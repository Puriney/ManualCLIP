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