library(pheatmap)
#library(heatmaply)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("ComplexHeatmap")

#library(ComplexHeatmap)

setwd("C:/Users/mikef/Dropbox (EinsteinMed)/Fei/Reverse/")
directory <- getwd()

day0<- read.csv("WTDAY0-vs-KODAY0/SigDEGs.csv",header = T)
day0<- day0$ID

day3<- read.csv("WTDAY3-vs-KODAY3/SigDEGs.csv",header = T)
day3<- day3$ID

day6<- read.csv("WTDAY6-vs-KODAY6/SigDEGs.csv",header = T)
day6<- day6$ID

common <- union(union(day0, day3),day6)
#common <- day6

write.csv(common, "common.csv", row.names = F)


common <- read.csv("common.csv",header = T)
commonDEGs <-common$x

Counts = read.table("WT-vs-KO/normalizedCounts.csv",header=T,sep=",",check.names=FALSE, quote = "\"")
colnames(Counts)[1] = "ID"

zscore <- function(x){
  z <- (x - mean(x))/sd(x)
  z
}
for (i in 2:length(Counts[1,])){
  Counts[,i] = zscore(Counts[,i])
}

CommonTPMCounts <- Counts[Counts$ID %in% commonDEGs,]
rownames(CommonTPMCounts) <- CommonTPMCounts$ID
CommonTPMCounts <- CommonTPMCounts[,-1]
#CommonTPMCounts <- CommonTPMCounts[,c(2,6,7,8,9,10,11,12,13,3,4,5)]
annotation_col <-data.frame(c("Day0", "Day0","Day0","Day0",
                              "Day3", "Day3","Day3","Day3",
                              "Day6", "Day6","Day6","Day6"),
                            c("WT-1", "WT-2","KO-1","KO-2",
                              "WT-1", "WT-2","KO-1","KO-2",
                              "WT-1", "WT-2","KO-1","KO-2"))
colnames(annotation_col) <- c("DAY","Group")
annotation_col <- annotation_col[,c(2,1)]
row.names(annotation_col) <- colnames(CommonTPMCounts)
par(mar = c(5, 4, 4, 2), # change the margins
    lwd = 2, # increase the line thickness
    cex.axis = 1.0 # increase default axis label size
)

png(file = "SigDEGs.png", width = 1650, height = 1600, units = "px", res = 300)
pheatmap(CommonTPMCounts, border_color=NA, treeheight_row=0, scale="row", 
         cluster_rows=T, show_rownames=F, show_colnames=F, cluster_cols=F, 
         annotation_col=annotation_col,
         fontsize = 6, angle_col ="90", main="Heat Map of DEGs (Expression Z-Score)")
dev.off()


selected <- c("Pax6","Otx2","Fgf5","Nes",
              "Fgfr2","Brachyury","Sox17","Gata2","Gata4","Gata6",
            "Cdx2","Eomes","Krt7","Egfr")

write.csv(selected, "selected.csv", row.names = F)

#day6<- read.csv("WT-vs-KO/SigDEGs.csv",header = T)

#NormalizedCounts
NormalizedCounts = read.table("WT-vs-KO/normalizedCounts.csv",header=T,sep=",",check.names=FALSE, quote = "\"")
colnames(NormalizedCounts)[1] = "ID"
NormalizedCounts$WTDay0 = (NormalizedCounts[,2]+NormalizedCounts[,3])/2
NormalizedCounts$KODay0 = (NormalizedCounts[,4]+NormalizedCounts[,5])/2
NormalizedCounts$WTDay3 = (NormalizedCounts[,6]+NormalizedCounts[,7])/2
NormalizedCounts$KODay3 = (NormalizedCounts[,8]+NormalizedCounts[,9])/2
NormalizedCounts$WTDay6 = (NormalizedCounts[,10]+NormalizedCounts[,11])/2
NormalizedCounts$KODay6 = (NormalizedCounts[,12]+NormalizedCounts[,13])/2
NormalizedCounts = NormalizedCounts[,c(1,14:19)]

#TPM_values
NormalizedCounts = read.table("WT-vs-KO/TPM_values.csv",header=T,sep=",",check.names=FALSE, quote = "\"")
colnames(NormalizedCounts)[1] = "ID"
NormalizedCounts = NormalizedCounts[,c(1,2,6,7,8,9,10,11,12,13,3,4,5)]
NormalizedCounts$WTDay0 = (NormalizedCounts[,2]+NormalizedCounts[,3])/2
NormalizedCounts$KODay0 = (NormalizedCounts[,4]+NormalizedCounts[,5])/2
NormalizedCounts$WTDay3 = (NormalizedCounts[,6]+NormalizedCounts[,7])/2
NormalizedCounts$KODay3 = (NormalizedCounts[,8]+NormalizedCounts[,9])/2
NormalizedCounts$WTDay6 = (NormalizedCounts[,10]+NormalizedCounts[,11])/2
NormalizedCounts$KODay6 = (NormalizedCounts[,12]+NormalizedCounts[,13])/2
NormalizedCounts = NormalizedCounts[,c(1,14:19)]

#FPKM
NormalizedCounts = read.table("WT-vs-KO/genes.fpkm_tracking",header=T,sep="\t",check.names=FALSE, quote = "\"")
NormalizedCounts = NormalizedCounts[,c(1,10,14,18,22,26,30)]
colnames(NormalizedCounts) = c("ID","WTDay0","KODay0","WTDay3","KODay3","WTDay6","KODay6")
#NormalizedCounts = NormalizedCounts[, c(1,6,7)]


# zscore <- function(x){
#   z <- (x - mean(x))/sd(x)
#   z
# }
# for (i in 2:length(NormalizedCounts[1,])){
#   NormalizedCounts[,i] = zscore(NormalizedCounts[,i])
# }

selected <- read.csv("selected.csv",header = T)
selectedDEGs <-selected$x

SelectedNormalizedCounts <- NormalizedCounts[NormalizedCounts$ID %in% selectedDEGs,]
rownames(SelectedNormalizedCounts) <- SelectedNormalizedCounts$ID
SelectedNormalizedCounts <- SelectedNormalizedCounts[,-1]
SelectedNormalizedCounts <- SelectedNormalizedCounts[selectedDEGs,]
SelectedNormalizedCounts <- SelectedNormalizedCounts[!is.na(SelectedNormalizedCounts$WTDay6),]

annotation_col <-data.frame(c("Ectodermal", "Ectodermal", "Ectodermal", "Ectodermal", "Mensendodermal",
                              "Mensendodermal", "Mensendodermal", "Mensendodermal", "Mensendodermal",
                              "Trophectodermal", "Trophectodermal", "Trophectodermal", "Trophectodermal"))
colnames(annotation_col) <- c("Group")
row.names(annotation_col) <- row.names(SelectedNormalizedCounts)

zscore <- function(x){
  z <- (x - mean(x))/sd(x)
  z
}
for (i in 1:length(SelectedNormalizedCounts[1,])){
  SelectedNormalizedCounts[,i] = zscore(SelectedNormalizedCounts[,i])
}

library(RColorBrewer)
# Sets the minimum (0), the maximum (15), and the increasing steps (+1) for the color scale
# Note: if some of your genes are outside of this range, they will appear white on the heatmap
breaksList = seq(-1, 2, by = 0.1)

# Plots the first heatmap
pheatmap(expressionData[1:10, ], # Plots the first 10 genes of the dataset
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList) # Sets the breaks of the color scale as in breaksList

# Plots the second heatmap with the same color and breaks options
pheatmap(expressionData[20:30, ], # Plots the third 10 genes of the dataset
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
         breaks = breaksList)


png(file = "EBDay6.png", width = 850, height = 1600, units = "px", res = 300)
pheatmap(SelectedNormalizedCounts, border_color="black", treeheight_row=0, scale = "row",
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breaksList)),
         breaks = breaksList,
         cluster_rows=F, show_rownames=T, show_colnames=T, cluster_cols=F, annotation_row=annotation_col,
         fontsize = 6, angle_col ="0", main="EB Day 6")


dev.off()
