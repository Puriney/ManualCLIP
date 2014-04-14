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