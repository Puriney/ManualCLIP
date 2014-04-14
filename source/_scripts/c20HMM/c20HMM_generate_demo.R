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