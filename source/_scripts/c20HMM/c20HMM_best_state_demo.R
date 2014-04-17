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