#!/usr/bin/Rscript

library(DESeq2)
library(gtools)
library(dplyr)

setwd("C:/Users/mikef/Dropbox (EinsteinMed)/RNA-Seq/WTDAY0-vs-KODAY0")
directory <- getwd()

comparison = "WTDAY0-vs-KODAY0"
refCond = sub("-vs-(.*)$", "", comparison)
FC = 0.0
negFC = FC*(-1)
Pvalue = 0.05
pvaluetype = "padj"

conditionfile <- "Conditions.txt"
Cond <- read.table(conditionfile, header=TRUE, row.names=1, fill=TRUE)
Cond$condition <- factor(Cond$condition)
Cond$batch <- factor(Cond$batch)
myCond <- Cond[mixedorder(row.names(Cond)),,drop=FALSE]
head(myCond)

countFiles <- grep("*.htseq.txt",list.files(directory),value=TRUE)
sampleNames <- sub(".htseq.txt","",countFiles)
countFiles <- countFiles[mixedorder(sampleNames),drop=FALSE]
sampleNames <- sub(".htseq.txt","",countFiles)

countData <- read.table(countFiles[1], row.names=1, check.names=FALSE, quote="\"", fill=TRUE, header=FALSE, sep="\t")
colnames(countData)[1] <- sampleNames[1]

countMatrix <- data.frame(matrix(ncol=0,nrow=nrow(countData)))
countMatrix <- cbind(countMatrix, countData)

for(i in 2:length(sampleNames)){
  countData <- read.table(countFiles[i], row.names=1, check.names=FALSE, quote="\"", fill=TRUE, header=FALSE, sep="\t")
  colnames(countData) <- sampleNames[i]
  countMatrix <- cbind(countMatrix, countData)
}

countMatrix <- countMatrix[-which(row.names(countMatrix) %in% 
                                    c("a", "__not_aligned", "__too_low_aQual", "__no_feature", 
                                      "__ambiguous", "__alignment_not_unique")),]

if (anyNA(myCond$batch)){
  dds = DESeqDataSetFromMatrix(countMatrix, myCond, design = ~ condition)
} else {
  dds = DESeqDataSetFromMatrix(countMatrix, myCond, design = ~ batch + condition)
}

dds$condition <- relevel(dds$condition, ref=refCond)
dds <- dds[rowSums(counts(dds)) > 10 * dim(myCond)[1]]

dds <- DESeq(dds)
res <- results(dds, alpha = 0.1)
res <- res[complete.cases(res),]
res

write.csv(countMatrix, file = "rawCounts.csv")

normalizedCounts <- counts(dds, normalized = T)
write.csv(normalizedCounts, file = "normalizedCounts.csv")

resOrdered <- res[order(res$padj),]
write.csv(as.data.frame(resOrdered), file="rawDEGs.csv")


SigDEGs <- subset(resOrdered, eval(as.name(pvaluetype)) < Pvalue & abs(log2FoldChange) > FC)
write.csv(as.data.frame(SigDEGs), file="rawSigDEGs.csv")
SigDEGsUp <- subset(resOrdered, eval(as.name(pvaluetype)) < Pvalue & (log2FoldChange) > FC)
SigDEGsDown <- subset(resOrdered, eval(as.name(pvaluetype)) < Pvalue & (log2FoldChange) < negFC)

SigDEGsList <- SigDEGs[0]
write.table(as.data.frame(SigDEGsList), "SigDEGsList.tsv", quote=FALSE, col.names=FALSE, sep="\t")
SigDEGsUpList <- SigDEGsUp[0]
write.table(as.data.frame(SigDEGsUpList), "SigDEGsUpList.tsv", quote=FALSE, col.names=FALSE, sep="\t")
SigDEGsDownList <- SigDEGsDown[0]
write.table(as.data.frame(SigDEGsDownList), "SigDEGsDownList.tsv", quote=FALSE, col.names=FALSE, sep="\t")


DEGs = read.table("rawDEGs.csv",header=T,sep=",")
SigDEGs = read.table("rawSigDEGs.csv", header=T,sep=",")
normalizedCounts = read.table("normalizedCounts.csv",header=T,sep=",",check.names=FALSE)
colnames(DEGs)[1] <- "ID"
colnames(SigDEGs)[1] <- "ID"
colnames(normalizedCounts)[1] <- "ID"
DEGNormalizedCounts = merge(DEGs,normalizedCounts,by="ID")
SigDEGNormalizedCounts = merge(SigDEGs,normalizedCounts,by="ID")
DEGNormalizedCounts = DEGNormalizedCounts[DEGNormalizedCounts$ID != "",]
SigDEGNormalizedCounts = SigDEGNormalizedCounts[SigDEGNormalizedCounts$ID != "",]

write.csv(DEGNormalizedCounts, file="DEGs.csv", row.names=F, quote=T)
write.csv(SigDEGNormalizedCounts, file="SigDEGs.csv", row.names=F, quote=T)

# Volcano Plots
png(file = "VolcanoPlot.png", width = 1600, height = 1600, units = "px", res = 300)
if ( pvaluetype == "padj" ){
  with(res, plot(log2FoldChange, -log10(padj), pch=20, xlim=c(min(res[,2]),max(res[,2])), ylim=c(0,25), main=comparison, cex.main=0.9, cex.lab=1, cex.axis=1))
  with(subset(res, padj < Pvalue & log2FoldChange > FC), points(log2FoldChange, -log10(padj), pch=20, col="red"))
  with(subset(res, padj < Pvalue & log2FoldChange < negFC), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
  dev.off()
} else{
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, xlim=c(min(res[,2]),max(res[,2])), main=comparison, cex.main=0.9, cex.lab=1, cex.axis=1))
  with(subset(res, pvalue < Pvalue & log2FoldChange > FC), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
  with(subset(res, pvalue < Pvalue & log2FoldChange < negFC), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
  dev.off()
}

# Box plots
png(file = "Boxplot_raw_data.png", width = 1600, height = 1600, units = "px", res = 300)
boxplot(log2(assays(dds)[["cooks"]]), range=0, las=2, main=comparison, pars = list(cex.axis=0.5, cex.main=0.9)) #raw data boxplot
dev.off()
png(file = "Boxplot_normalized_data.png", width = 1600, height = 1600, units = "px", res = 300)
boxplot(log2(counts(dds,normalized=TRUE)), range=20, las=2, main=comparison, pars = list(cex.axis=0.5, cex.main=0.9))
dev.off()








# Barplots

#install.packages("devtools")
#library(devtools)
#devtools::install_github("kassambara/easyGgplot2", force = TRUE)
####library(easyGgplot2)

setwd("C:/Users/Fei/Dropbox (EinsteinMed)/Fei/Qin-None-Stranded/")
#setwd("C:/Users/mikef/Dropbox (EinsteinMed)/RNA-Seq/")
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

# # Add abline
# abline(v=c(1.3 , 2.5) , col="grey")

# Add the text 
text(my_bar, (df$Down/2), df$Down , cex=1.2, col="#FFFFFF") 
text(my_bar, (df$Down+df$Up/2), df$Up , cex=1.2, col="#FFFFFF") 
#text(my_bar, (df$Down+df$Up+100), paste0("Total: ",df$Up+df$Down) , cex=1.2, col="#000000") 

#Legend
legend("topleft", legend = c("Down","Up" ) , 
       col = c(rgb(0.0,0.0,1.0,1.0), rgb(1.0,0.0,0.0,1.0)) , 
       bty = "n", pch=15 , pt.cex = 2, cex = 1.2, horiz = FALSE, inset = c(0.05, 0.05))

####axis.break(500,500,style="gap")
dev.off()
#http://www.sthda.com/english/wiki/ggplot2-barplots-quick-start-guide-r-software-and-data-visualization

