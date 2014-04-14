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
