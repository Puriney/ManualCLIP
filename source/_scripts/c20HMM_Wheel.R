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
intron.emmision.prob2 <- data.frame(Base = paste(dna.bases, as.character(intron.emission.prob), sep="="), Prob = intron.emission.prob)
p <- ggplot(intron.emmision.prob2, aes(x="", y=Prob, fill = Base)) + 
        geom_bar(width = 3) + 
        coord_polar("y") + 
        xlab('') + ylab('') + ggtitle('Intron Roulette Wheel') + labs(fill="Base") + 
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



