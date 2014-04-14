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