#install.packages("devtools")
#library(devtools)
#devtools::install_github("kassambara/easyGgplot2", force = TRUE)
library(easyGgplot2)

setwd("Your project DEG directory")
directory <- getwd()

day0<- read.csv("WTDAY0-vs-KODAY0/SigDEGs.csv",header = T)
day0Up<- nrow(day0[day0$log2FoldChange > 0,])
day0Down<- nrow(day0[day0$log2FoldChange < 0,])

day3<- read.csv("WTDAY3-vs-KODAY3/SigDEGs.csv",header = T)
day3Up<- nrow(day3[day3$log2FoldChange > 0,])
day3Down<- nrow(day3[day3$log2FoldChange < 0,])

day6<- read.csv("WTDAY6-vs-KODAY6/SigDEGs.csv",header = T)
day6Up<- nrow(day6[day6$log2FoldChange > 0,])
day6Down<- nrow(day6[day6$log2FoldChange < 0,])

# Define the data frame for each comparison
df <- data.frame(DAY = c("Day 0","Day 3","Day 6"), 
                 Up = c(day0Up, day3Up, day6Up), 
                 Down= c(day0Down, day3Down, day6Down))

# png(file = "BarPlot.png", width = 1600, height = 1600, units = "px", res = 300)
pdf("BarPlot.pdf", width = 7, height = 14)
my_bar <- barplot(cbind(Down, Up) ~ DAY, data = df, border="black", width=c(1,1,1),
                  las=0 , 
                  col=c(rgb(0.0,0.0,1.0,1.0), rgb(1.0,0.0,0.0,1.0)), 
                  ylim=c(0,5000), 
                  space=0.3,
                  axes=T,
                  axisnames=T,
                  xlab = "",
                  ylab = "Number of DEGs",
                  main="" )

## Draw x-axis without labels.
axis(side = 1, at=c(0.8, 2.1, 3.4), labels = F)
axis(side = 1, at=c(0.8, 2.1, 3.4), labels = F)

# If need to Add vertical abline to separate the bars
# abline(v=c(1.3 , 2.5) , col="grey")

# Add the text 
text(my_bar, (df$Down/2), df$Down , cex=1.2, col="#FFFFFF") 
text(my_bar, (df$Down+df$Up/2), df$Up , cex=1.2, col="#FFFFFF") 
#text(my_bar, (df$Down+df$Up+100), paste0("Total: ",df$Up+df$Down) , cex=1.2, col="#000000") 

#Legend
legend("topleft", legend = c("Down","Up" ) , 
       col = c(rgb(0.0,0.0,1.0,1.0), rgb(1.0,0.0,0.0,1.0)) , 
       bty = "n", pch=15 , pt.cex = 2, cex = 1.2, horiz = FALSE, inset = c(0.05, 0.05))

axis.break(500,500,style="gap")

dev.off()

