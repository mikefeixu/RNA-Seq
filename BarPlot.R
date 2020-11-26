#install.packages("devtools")
#library(devtools)
#devtools::install_github("kassambara/easyGgplot2", force = TRUE)
library(easyGgplot2)

setwd("C:/Users/mikef/Dropbox (EinsteinMed)/RNA-Seq/")
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

df <- data.frame(DAY = c("Day 0","Day 3","Day 6"), 
                 Up = c(day0Up, day3Up, day6Up), 
                 Down= c(day0Down, day3Down, day6Down))

png(file = "BarPlot.png", width = 1600, height = 1600, units = "px", res = 300)
par(mar = c(6.1, 4.1, 4.1, 4.1), # change the margins
    lwd = 2, # increase the line thickness
    cex.axis = 1.2 # increase default axis label size
)

#A numerical vector of the form c(bottom, left, top, right) which gives
#the number of lines of margin to be specified on the four sides of the plot. The default is c(5, 4, 4, 2) + 0.1.
## Draw boxplot with no axes.
my_bar <- barplot(cbind(Down, Up) ~ DAY, data = df, border="black", width=c(1,1,1),
                  las=0 , 
                  col=c(rgb(0.0,0.0,1.0,1.0), rgb(1.0,0.0,0.0,1.0)), 
                  ylim=c(0,2500), 
                  space=0.3,
                  axes=T,
                  axisnames=T,
                  xlab = "",
                  ylab = "Number of DEGs",
                  main="" )

## Draw x-axis without labels.
axis(side = 1, at=c(0.8, 2.1, 3.4), labels = F)
axis(side = 1, at=c(0.8, 2.1, 3.4), labels = F)

# # Add abline
# abline(v=c(1.3 , 2.5) , col="grey")

# Add the text 
text(my_bar, (df$Down/2), df$Down , cex=1.2, col="#FFFFFF") 
text(my_bar, (df$Down+df$Up/2), df$Up , cex=1.2, col="#FFFFFF") 
text(my_bar, (df$Down+df$Up+100), paste0("Total: ",df$Up+df$Down) , cex=1.2, col="#000000") 

#Legend
legend("topleft", legend = c("Down","Up" ) , 
       col = c(rgb(0.0,0.0,1.0,1.0), rgb(1.0,0.0,0.0,1.0)) , 
       bty = "n", pch=15 , pt.cex = 2, cex = 1.2, horiz = FALSE, inset = c(0.05, 0.05))

dev.off()

# ## Draw y-axis.
# axis(side = 2,
#      ## Rotate labels perpendicular to y-axis.
#      las = 0,
#      ## Adjust y-axis label positions.
#      mgp = c(3, 0.75, 0), labels = F)

# ## Draw the x-axis labels.
# text(x = 1:length(dat),
#      ## Move labels to just below bottom of chart.
#      y = par("usr")[3] - 0.45,
#      ## Use names from the data list.
#      labels = names(dat),
#      ## Change the clipping region.
#      xpd = NA,
#      ## Rotate the labels by 35 degrees.
#      srt = 35,
#      ## Adjust the labels to almost 100% right-justified.
#      adj = 0.965,
#      ## Increase label size.
#      cex = 1.2)

#http://www.sthda.com/english/wiki/ggplot2-barplots-quick-start-guide-r-software-and-data-visualization