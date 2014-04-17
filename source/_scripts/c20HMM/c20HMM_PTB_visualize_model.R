## -----------------------------
## Visualize PTB HMM
## -----------------------------
ptb.HMM.emisMx.order <- ptb.HMM.emisMx[, order(ptb.HMM.emisMx[1, ], decreasing=F)]
library(reshape2)
ptb.HMM.emisMx.order <- as.data.frame(ptb.HMM.emisMx.order)
ptb.HMM.emisMx.order$State <- row.names(ptb.HMM.emisMx.order)
ptb.HMM.emisMx.order.melt <- melt(ptb.HMM.emisMx.order, id.vars = "State")
p <- ggplot(data = ptb.HMM.emisMx.order.melt, 
       aes(x = variable, y = value, fill = State)) + 
    geom_bar(stat="identity", position=position_dodge()) + coord_flip() + 
    theme_bw() + 
    theme(text = element_text(family = "xkcd", size = 20), 
          axis.text.x = element_text(angle = 90, hjust = 1), 
          legend.justification = c(0, 0), 
          legend.position=c(0.6, 0.7)) + 
    xlab("Triplet") + ylab("Emission Prob") 