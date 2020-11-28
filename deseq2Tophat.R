#!/usr/bin/Rscript

library(DESeq2)
library(gtools)
library(dplyr)

setwd("C:/Users/mikef/Dropbox (EinsteinMed)/RNA-Seq/WTDAY6-vs-KODAY6")
directory <- getwd()

comparison = "WTDAY6-vs-KODAY6"
refCond = sub("-vs-(.*)$", "", comparison)
FC = 0
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

#plotMA(dds)
