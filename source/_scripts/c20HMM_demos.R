## ---------------
## Scripts for c20_Model_HMM
## Yun Yan <yanyun@whu.edu.cn>
## ---------------
library(ggplot2)
library(xkcd)
dna.bases <- c("A", "C", "G", "T")
seq.path <- "CTTCATGTGAAAGCAGACGTAAGTCA"

## ----------------------------
## Emission Prob Demo (Intron)
## ----------------------------
## Pre-settings
intron.emission.prob <- c(0.40, 0.10, 0.10, 0.40)
intron.emmision.matrix <- matrix(intron.emission.prob, nrow = 1, byrow = T)
colnames(intron.emmision.matrix) <- dna.bases
row.names(intron.emmision.matrix) <- c("I")
## Plot wheel of intron
intron.emmision.prob2 <- data.frame(Base = 
                                        paste(dna.bases, 
                                              as.character(intron.emission.prob),
                                              sep="="),
                                    Prob = intron.emission.prob)
p <- ggplot(intron.emmision.prob2, aes(x="", y=Prob, fill = Base)) + 
        geom_bar(width = 1) + 
        coord_polar("y") + 
        xlab('') + ylab('') + ggtitle('Intron Roulette Wheel') + 
        labs(fill="Base") + 
        theme_xkcd() + 
        theme(axis.ticks=element_blank(), text = element_text(size = 30)) 
p 
## Sample bases
set.seed(1026)
seq.len <- 26
sample.seq <- sample(x = dna.bases, prob = intron.emission.prob,
                     size = seq.len, replace = T)
cat(sample.seq, "\n")
table(sample.seq)

## ----------------------------
## Transition Prob Demo
## ----------------------------
state.bases <- c("E", "5", "I")
exon.transition.prob <- c(0.9, 0.1, 0)
exon.transition.matrix <- matrix(exon.transition.prob, nrow = 1, byrow = T)
colnames(exon.transition.matrix) <- state.bases
row.names(exon.transition.matrix) <- "E"
print(exon.transition.matrix)
set.seed(1026)
sample.state <- sample(x = state.bases, prob = exon.transition.prob,
                       size = 1, replace = T)
cat(sample.state)

## -----------------------------
## Generating sequence
## -----------------------------
## Construct transition matrix
S <- c(0,1,rep(0,3))
E <- c(0, 0.9, 0.1, 0, 0)
SS5 <- c(rep(0,3), 1, 0)
I <- c(rep(0,3), 0.9, 0.1)
N <- c(rep(0,4), 1)
myTransMtx <- rbind(S,E, SS5, I, N)
colnames(myTransMtx) <- rownames(myTransMtx)
print(myTransMtx)
## Construct emission matrix
S <- rep(0,4)
E <- rep(1/4, 4)
SS5 <- c(0.05, 0, 0.95, 0)
I <- c(0.4, 0.1, 0.1, 0.4)
N <- c(rep(0, 4))
myEmisMtx <- rbind(S,E, SS5, I, N)
colnames(myEmisMtx) <- c("A", "C", "G", "T")
print(myEmisMtx)
## Function to generate sequence
GenerateHMMSeq <- function(emmisionMx, transitionMx, len){
    state.bases <- row.names(emmisionMx)
    dna.bases   <- colnames(emmisionMx)
    out.seq.path <- rep(NA, len)
    out.state.path <- rep(NA, len)
    ## State of first base is "Start"
    state0 <- "S"
    cat("Start\n")
    
    for (i in 1:len) {
        ## Choose State
        state.i <- sample(x = state.bases, prob = transitionMx[state0, ], 
                          replace = T, size = 1)
        if(state.i == "N"){
            out.state.path[i] <- state.i
            cat("End\n")
            break
        }
        base.i  <- sample(x = dna.bases, prob = emmisionMx[state.i, ], 
                          replace = T, size = 1)
        cat("Base", i, "is", base.i, "under", state.i, "\n")
        out.state.path[i] <- state.i
        out.seq.path[i] <- base.i
        state0 <- state.i
    }
    cat("Seq Path:", na.omit(out.seq.path), "\n")
    cat("StatPath:", na.omit(out.state.path), "\n")
    return(1)
}

## --------------------------
## Compute Prob of state path
## --------------------------
ComputeStatePathProb <- function(seqPath, statePath, transMtx, emisMtx, log=TRUE){
    seqPath <- unlist(strsplit(toupper(seqPath), ""))
    statePath <- unlist(strsplit(statePath, ""))
    statePath[statePath=="5"] <- "SS5"
    ## First base
    ## previous is under Start State
    prob <-  transMtx["S", statePath[1]] * emisMtx[statePath[1], seqPath[1]]
    for (i in 2:length(seqPath)){
        state.i <- statePath[i]
        nucleotide.i <- seqPath[i]
        prob.i <- transMtx[statePath[i-1],state.i] * emisMtx[statePath[i], nucleotide.i]
        prob <- prob * prob.i
    }
    ## End State
    prob <- prob * transMtx[statePath[length(statePath)], "N"]
    if (log){
        return (log2(prob)/log2(exp(1))) 
    } else {
        return (prob)
    }
}

## --------------------
## Pick out the best State path
## --------------------
InspectorHMM <- function(obsSeq, transMtx, emisMtx ){
    sumProb <- 0
    allPath <- NULL
    allProb <- NULL
    print ("All State Path with non-zero Probability: ")
    for (i in 1:24) {
        path_i <- c(rep("E", i), "5", rep("I", 25-i))
        path_i <- paste(path_i, collapse="")
        allPath[i] <- path_i
        
        p <- computeProb(obsSeq, path_i,transMtx, emisMtx,log=FALSE)
        allProb[i] <- p 
        if (p != 0){
            print (path_i)
            print (p)
            sumProb <- sumProb + p 
        }
    }
    maxProb <- max(allProb, na.rm=TRUE)
    max <- which(allProb == maxProb)
    cat ("Best State Path and original Sequence", "\n")
    cat (allPath[max], "\n")
    cat (obsSeq,"\n")
    cat (paste("Maximum Probability: ", maxProb, sep=""),"\n")
    cat (paste("Maximum Prob (log):  ", log2(maxProb)/log2(exp(1)), sep=""),"\n")
    cat (paste("Posterior Decoding:  ", maxProb/sumProb, sep=""),"\n")
    return (allPath[max])
}

## ---------------------
## PTB CLIP-seq data
## ---------------------
## ---------------------
## Data retrieve and processing
## ---------------------

## Genome fasta
source("http://bioconductor.org/biocLite.R")
biocLite("BSgenome.Hsapiens.UCSC.hg18") 
# ~800Mb, or alternatively download and install the souce pkg
# my.cmd <- "axel http://www.bioconductor.org/packages/release/data/annotation/src/contrib/BSgenome.Hsapiens.UCSC.hg18_1.3.19.tar.gz" 
# system(my.cmd)
# install.packages("BSgenome.Hsapiens.UCSC.hg18_1.3.19.tar.gz", repos = NULL, type = "source")
library("BSgenome.Hsapiens.UCSC.hg18")
hg18.genome <- BSgenome.Hsapiens.UCSC.hg18

## Peak bed corrdinates
download.file(url = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE19nnn/GSE19323/suppl/GSE19323_ptb_peak_cluster.bed.gz", destfile = "ptb_peak_bed.gz")
my.cmd <- "gunzip ptb_peak_bed.gz" # file.name would be ptb_peak_bed
system(my.cmd)
ptb.peak <- read.table(file = "ptb_peak_bed", header = FALSE)
head(ptb.peak)
colnames(ptb.peak) <- c("Chr", "Start", "End", "Peak", "Height", "Strand")
neg.indices <- which(ptb.peak[, 6] == "-")

## Fetch CLIP-seq Cluster Sequence
library("GenomicRanges")
ptb.peak.GRobj <- GRanges(seqnames = Rle(as.character(ptb.peak[, 1])), 
                          ranges = IRanges(ptb.peak[, 2], 
                                           ptb.peak[, 3],
                                           names = ptb.peak[, 4]), 
                          strand = Rle(ptb.peak[, 6]))
ptb.peak.seq <- getSeq(hg18.genome, seqnames(ptb.peak.GRobj), 
                       start(ptb.peak.GRobj), 
                       end(ptb.peak.GRobj)) ## raw genome seq
ptb.peak.seq[neg.indices] <- reverseComplement(ptb.peak.rna[neg.indices]) ## transcript sequence 

## Cook to observations
library("plyr")
library("snowfall")
ptb.peak.seq <- as.list(as.character(ptb.peak.seq))
Seq2Obv <- function(dataList, size) {
    pass.indices <- sapply(dataList, nchar)
    dataList <- dataList[pass.indices >= size]
    dataList <- llply(dataList, function(ele) {
        ele <- as.character(dna2rna(DNAString(ele)))
        start.i <- 1
        ele <- substring(ele, 
                         seq(start.i, nchar(ele), size), 
                         seq((start.i + size - start.i), nchar(ele), size))
        return(ele[ele != ""])
    })
    dataList
}

Seq2ObvSuper <- function(dataList, size, startChar) {
    pass.indices <- sapply(dataList, nchar)
    dataList <- dataList[pass.indices >= (size + startChar - 1)]
    dataList <- sfLapply(dataList, function(ele) {
        ele <- as.character(dna2rna(DNAString(ele)))
        start.i <- startChar
        ele <- substring(ele, 
                         seq(start.i, nchar(ele), size), 
                         seq((size + start.i - 1), nchar(ele), size))
        return(ele[ele != ""])
    })
    dataList
}

ptb.peak.obv1 <- Seq2ObvSuper(ptb.peak.seq, size = 3, startChar = 1)
ptb.peak.obv2 <- Seq2ObvSuper(ptb.peak.seq, size = 3, startChar = 2)
ptb.peak.obv3 <- Seq2ObvSuper(ptb.peak.seq, size = 3, startChar = 3)
ptb.peak.obvFull <- c(ptb.peak.obv1, ptb.peak.obv2, ptb.peak.obv3)
## -----------------------------
## Estimate HMM of PTB
## -----------------------------
require("RHmm")
ptb.HMM.estimated <- HMMFit(ptb.peak.obvFull, dis = "DISCRETE", nStates = 2,
       control = list(init = "RANDOM"))
ptb.HMM.estimated.maxLLH <- ptb.HMM.estimated$LLH
for (i in 1:20) {
    ptb.HMM.estimated.each <- HMMFit(ptb.peak.obvFull, dis = "DISCRETE", 
                                     nStates = 2, 
                                     control = list(init = "RANDOM"))
    ptb.HMM.estimated.each.LLH <- ptb.HMM.estimated.each$LLH
    cat(ptb.HMM.estimated.maxLLH, "vs", ptb.HMM.estimated.each.LLH, "\n")
    if (ptb.HMM.estimated.each.LLH > ptb.HMM.estimated.maxLLH) {
        ptb.HMM.estimated <- ptb.HMM.estimated.each
        ptb.HMM.estimated.maxLLH <- ptb.HMM.estimated.each.LLH
    }
}
print(ptb.HMM.estimated)

ptb.HMM.transMx <- ptb.HMM.estimated$HMM$transMat
ptb.HMM.emisMx  <- ptb.HMM.estimated$HMM$distribution$proba
ptb.HMM.nStates <- ptb.HMM.estimated$HMM$distribution$nStates
ptb.HMM.nSymbols <- ptb.HMM.estimated$HMM$distribution$nLevels
ptb.HMM.symbols <- names(ptb.HMM.emisMx[[ptb.HMM.nStates]])
ptb.HMM.emisMx  <- matrix(unlist(ptb.HMM.emisMx), byrow = T, 
                          ncol = ptb.HMM.nSymbols, 
                          nrow = ptb.HMM.nStates)
colnames(ptb.HMM.emisMx) <- ptb.HMM.symbols
row.names(ptb.HMM.emisMx) <- c("State 1", "State 2")
# write.table(ptb.HMM.emisMx, file="ptb_emission.txt", col.names=T, quote=F, sep="\t")

## -----------------------------
## Visualize PTB HMM
## -----------------------------
ptb.HMM.emisMx.order <- ptb.HMM.emisMx[, order(ptb.HMM.emisMx[1, ], decreasing=F)]
library(reshape2)
ptb.HMM.emisMx.order <- as.data.frame(ptb.HMM.emisMx.order)
ptb.HMM.emisMx.order$State <- row.names(ptb.HMM.emisMx.order)
ptb.HMM.emisMx.order.melt <- melt(ptb.HMM.emisMx.order, id.vars = "State")
p <- ggplot(data = ptb.HMM.emisMx.order.melt, 
       aes(x = variable, y = value, fill = State)) + 
    geom_bar(stat="identity", position=position_dodge()) + coord_flip() + 
    theme_bw() + 
    theme(text = element_text(family = "xkcd", size = 20), 
          axis.text.x = element_text(angle = 90, hjust = 1), 
          legend.justification = c(0, 0), 
          legend.position=c(0.6, 0.7)) + 
    xlab("Triplet") + ylab("Emission Prob") 
    

    
## -----------------------------
## Construct PTB HMM
## -----------------------------
require("gtools")

rna.bases <- c("A", "C", "G", "U")
ptb.triplet <- permutations(n = 4, r = 3, v = rna.bases, 
                            repeats.allowed = TRUE)
ptb.triplet <- paste(ptb.triplet[, 1], ptb.triplet[, 2], ptb.triplet[, 3], 
                     sep = "")
ptb.triplet.len <- length(ptb.triplet)

ptb.state.bases <- c("Binding", "Nonbinding")
ptb.state.len <- length(ptb.state.bases)
ptb.init.state.prob <- c(0.5, 0.5)

ptb.trans.mx <- matrix(rep(rep(1/ptb.state.len, ptb.state.len), ptb.state.len), 
                       nrow = ptb.state.len, byrow = TRUE)
colnames(ptb.trans.mx) <- ptb.state.bases
row.names(ptb.trans.mx) <- ptb.state.bases

ptb.emis.mx <- matrix(rep(rep(1/ptb.triplet.len, ptb.triplet.len), ptb.state.len),
                      nrow = ptb.state.len, byrow = TRUE)
row.names(ptb.emis.mx) <- ptb.state.bases
colnames(ptb.emis.mx) <- ptb.triplet
ptb.emis.list <- split(ptb.emis.mx, 1:NROW(ptb.emis.mx))
names(ptb.emis.list) <- ptb.state.bases
## Init PTB HMM
ptb.HMM <- HMMSet(initProb = ptb.init.state.prob, transMat = ptb.trans.mx, 
                  dis = "DISCRETE", 
                  proba = ptb.emis.list,
                  labels = ptb.triplet)
