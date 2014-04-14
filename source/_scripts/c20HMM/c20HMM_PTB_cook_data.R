
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