#!/usr/bin/Rscript
# Config comparison and projectdir before running the script.
# Create a folder for each comparion, add the counts files for samples from both groups to this folder
# Add Conditions.txt to each comparion folder, and update sample, group, and batch info for according comparions. Leave batch info blank if no batch effect expected.
# Add SelectedLabledGenes.csv to Selected_genes folder, and update it with the genes you want to lable in the MA plot and Volcano plot
library(DESeq2)
library(gtools)
library(dplyr)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
comparison <- "Your comparison" # Two groups connected by "-vs-": eg. "Control-vs-Treatment"
projectdir <- "Your project DEG directory" # Please use "/" or "\\" in the directory
workdir <- paste0(projectdir, comparison)
setwd(workdir) # Path to your comparison folder: eg. "./WTDAY0-vs-KODAY0"

refCond = sub("-vs-(.*)$", "", comparison)
FC = 1
negFC = FC*(-1)
Pvalue = 0.05
pvaluetype = "padj"

conditionfile <- "Conditions.txt"
Cond <- read.table(conditionfile, header=TRUE, row.names=1, fill=TRUE)
Cond$condition <- factor(Cond$condition)
Cond$batch <- factor(Cond$batch)
myCond <- Cond[mixedorder(row.names(Cond)),,drop=FALSE]
head(myCond)
sample_amount <- nrow(myCond)

countFiles <- grep("*.counts.txt",list.files(workdir),value=TRUE)
sampleNames <- sub(".counts.txt","",countFiles)
countFiles <- countFiles[mixedorder(sampleNames),drop=FALSE]
sampleNames <- sub(".counts.txt","",countFiles)

countData <- read.table(countFiles[1], row.names=1, check.names=FALSE, quote="\"", fill=TRUE, header=TRUE, sep="\t")
countData <- countData[, 5:ncol(countData), drop=FALSE]

countMatrix <- data.frame(matrix(ncol=0,nrow=nrow(countData)))
countMatrix <- cbind(countMatrix, countData)
colnames(countMatrix)[1] <- "length"
colnames(countMatrix)[2] <- sampleNames[1]

for(i in 2:length(sampleNames)){
  countData <- read.table(countFiles[i], row.names=1, check.names=FALSE, quote="\"", fill=TRUE, header=TRUE, sep="\t")
  countData <- countData[, 6:ncol(countData), drop=FALSE]
  colnames(countData) <- sampleNames[i]
  countMatrix <- cbind(countMatrix, countData)
}

countMatrix <- countMatrix[ ,2:ncol(countMatrix)]


if (anyNA(myCond$batch)){
  dds = DESeqDataSetFromMatrix(countMatrix, myCond, design = ~ condition)
} else {
  dds = DESeqDataSetFromMatrix(countMatrix, myCond, design = ~ batch + condition)
}

dds$condition <- relevel(dds$condition, ref=refCond)
dds <- dds[rowSums(counts(dds)) > 10 * dim(myCond)[1]]

dds <- DESeq(dds)
res <- results(dds)
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
  with(res, plot(log2FoldChange, -log10(padj), pch=20, xlim=c(min(res[,2]),max(res[,2])), ylim=c(0,250), main=comparison, cex.main=0.9, cex.lab=1, cex.axis=1))
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


# MA plots
# For help: https://github.com/kassambara/ggpubr/issues/70
library(ggpubr)
res<-res[order(res$padj, decreasing = T),]

# Get list of Genes for labling 
SelectedLabledGenes <- read.csv("../Selected_genes/SelectedLabledGenes.csv",header = F, sep = "\n")
SelectedLabledGenes <- SelectedLabledGenes$V1

# res <- subset(res, res@rownames %in% SelectedLabledGenes)  #Use subset for selecting genes on MA plot.
pma <- ggmaplot(res, main = comparison,
         fdr = 0.05, fc = 2, size = 0.4, alpha = 1,
         palette = c("red", "blue", "darkgray"),
         genenames = as.vector(res$name),
         legend = "top", top = 0,
         font.legend = "bold", 
         label.rectangle=T,
         label.select=SelectedLabledGenes,
         font.main = "bold",
         ylim=c(-5,5),
         ggtheme = ggplot2::theme_minimal())

pdf("MA.pdf")
pma
dev.off()
png(file = "MA.png", width = 1600, height = 1600, units = "px", res = 300)
pma
dev.off()

# Labeled Volcano Plots
# Get DEGs data
results <- read.table("rawDEGs.csv",header=T,sep=",")
results$genelabels <- factor(results$X, levels = SelectedLabledGenes)
results$Sig <- "NonSig"
results[results$log2FoldChange>=1, "Sig"]="Up"
results[results$log2FoldChange<=-1, "Sig"]="Down"
results$Sig <- as.factor(results$Sig)

VolcanoPlotLabled <- ggplot(results, aes(log2FoldChange, -log10(padj), label = genelabels, col = Sig)) + 
  geom_point() +
  geom_text_repel(col = "black", na.rm = TRUE, box.padding = unit(0.45, "lines"), hjust = 1, max.overlaps=10) + 
  scale_color_manual(values = c("blue", "grey","red")) +
  theme(legend.title = element_blank(), text = element_text(size = 12))+  
  ylim(c(0,300)) +
  xlim(c(-10,10)) +
  #geom_hline(yintercept = -log10(0.05), col = "green")+  #Uncomment if need the threshold line for -log10(padj).
  geom_vline(xintercept = c(1,-1), col = "purple")

png(file = "VolcanoPlotLabled.png", width = 1600, height = 1600, units = "px", res = 300)
VolcanoPlotLabled
dev.off()

pdf("VolcanoPlotLabled.pdf")
VolcanoPlotLabled
dev.off()

# Unlabeled Volcano Plots
png(file = "VolcanoPlot.png", width = 1600, height = 1600, units = "px", res = 300)
if ( pvaluetype == "padj" ){
  with(res, plot(log2FoldChange, -log10(padj), pch=20, xlim=c(min(res[,2]),max(res[,2])), ylim=c(0,300), main=comparison, cex.main=0.9, cex.lab=1, cex.axis=1))
  with(subset(res, padj < Pvalue & log2FoldChange > FC), points(log2FoldChange, -log10(padj), pch=20, col="red"))
  with(subset(res, padj < Pvalue & log2FoldChange < negFC), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
  abline(-log10(0.05), 0, col = "green")
  dev.off()
} else{
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, xlim=c(min(res[,2]),max(res[,2])), ylim=c(0,10), main=comparison, cex.main=0.9, cex.lab=1, cex.axis=1))
  with(subset(res, pvalue < Pvalue & log2FoldChange > FC), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
  with(subset(res, pvalue < Pvalue & log2FoldChange < negFC), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
  dev.off()
}
pdf("VolcanoPlot.pdf")
if ( pvaluetype == "padj" ){
  with(res, plot(log2FoldChange, -log10(padj), pch=20, xlim=c(-10,10), ylim=c(0,300),  main=comparison, cex.main=0.9, cex.lab=1, cex.axis=1))
  with(subset(res, padj < Pvalue & log2FoldChange > FC), points(log2FoldChange, -log10(padj), pch=20, col="red"))
  with(subset(res, padj < Pvalue & log2FoldChange < negFC), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
  abline(-log10(0.05), 0, col = "green")
  abline(v=1, b=90, col = "purple")
  abline(v=-1, b=90, col = "purple")
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


# Dispersion plot and fitting alternatives
png(file = "Dispersion_plot.png", width = 1600, height = 1600, units = "px", res = 300)
plotDispEsts(dds, main=comparison)
dev.off()


# PCA plot
rld <- rlog(dds)
png(file = "PCA.png", width = 1600, height = 1600, units = "px", res = 300)
pcaData <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
colnames(pcaData)[3:5] <- c("Group", "Condition", "Name")
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaplot <- ggplot(pcaData, aes(PC1, PC2, colour=Condition)) +
  geom_point(size=1) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() +
  ggtitle(comparison) +
  theme(plot.title = element_text(vjust = 5)) +
  # ylim(c(-20,20)) +
  # xlim(c(-3,3)) +
  theme(plot.title = element_text(hjust = 0.5))
plot(pcaplot)
dev.off()
pdf("PCA.pdf")
plot(pcaplot)
dev.off()


# Sample distance plots
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste( row.names(colData(rld)), rld$condition, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
png(file = "SampleDistances.png", width = 1600, height = 1600, units = "px", res = 300)
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists,clustering_distance_cols=sampleDists,col=colors, main=comparison, fontsize=5)
dev.off()
pdf("SampleDistances.pdf")
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists,clustering_distance_cols=sampleDists,col=colors, main=comparison, fontsize=5)
dev.off()
