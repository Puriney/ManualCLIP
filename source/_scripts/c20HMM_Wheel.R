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
        opts(axis.ticks=theme_blank())
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

## --------------------------------------
## Force Compute Most Possible State Path
## --------------------------------------
InspectorHMM <- function(seqPath, transMtx, emisMtx ){
    sumProb <- 0
    allPath <- NULL
    allProb <- NULL
    print ("All State Path with non-zero Probability: ")
    for (i in 1:24) {
        path_i <- c(rep("E", i), "5", rep("I", 25-i))
        path_i <- paste(path_i, collapse="")
        allPath[i] <- path_i
        
        p <- computeProb(seqPath, path_i,transMtx, emisMtx,log=FALSE)
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
    cat (seqPath,"\n")
    cat (paste("Maximum Probability: ", maxProb, sep=""),"\n")
    cat (paste("Maximum Prob (log):  ", log2(maxProb)/log2(exp(1)), sep=""),"\n")
    cat (paste("Posterior Decoding:  ", maxProb/sumProb, sep=""),"\n")
    return (allPath[max])
}
