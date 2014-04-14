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