library(pheatmap)
library(heatmaply)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("ComplexHeatmap")

library(ComplexHeatmap)

setwd("C:/Users/mikef/Dropbox (EinsteinMed)/RNA-Seq/")
directory <- getwd()

day0<- read.csv("WTDAY0-vs-KODAY0/SigDEGs.csv",header = T)
day0<- day0$ID

day3<- read.csv("WTDAY3-vs-KODAY3/SigDEGs.csv",header = T)
day3<- day3$ID

day6<- read.csv("WTDAY6-vs-KODAY6/SigDEGs.csv",header = T)
day6<- day6$ID

common <- intersect(intersect(day0, day3),day6)
common <- day6


# common <- common[common != "Tmsb10" & common != "Myh9"
#                  & common != "Lgals1" & common != "Acsl4"
#                  & common != "Tpm1" & common != "Dsp"
#                  & common != "Gm21411" & common != "Erf"]
write.csv(common, "common.csv", row.names = F)


common <- read.csv("common.csv",header = T)
commonDEGs <-common$x

TPMCounts = read.table("WT-vs-KO/normalizedCounts.csv",header=T,sep=",",check.names=FALSE, quote = "\"")
colnames(TPMCounts)[1] = "ID"

zscore <- function(x){
  z <- (x - mean(x))/sd(x)
  z
}
for (i in 2:length(TPMCounts[1,])){
  TPMCounts[,i] = zscore(TPMCounts[,i])
}

CommonTPMCounts <- TPMCounts[TPMCounts$ID %in% commonDEGs,]
rownames(CommonTPMCounts) <- CommonTPMCounts$ID
CommonTPMCounts <- CommonTPMCounts[,-1]
#CommonTPMCounts <- CommonTPMCounts[,c(2,6,7,8,9,10,11,12,13,3,4,5)]
# colnames(CommonTPMCounts) <- c("Day0 WT-1", "Day0 WT-2","Day0 KO-1","Day0 KO-2",
#                                "Day3 WT-1", "Day3 WT-2","Day3 KO-1","Day3 KO-2",
#                                "Day6 WT-1", "Day6 WT-2","Day6 KO-1","Day6 KO-2")
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


selected <- c("Pax6","Oct2","Fgf5","Nestin","Cxxc5",
              "Fgfr2","Brachyury","Sox17","Gata2","Gata4","Gata6",
            "Cdx2","Eomes","Krt7","Egfr")

write.csv(selected, "selected.csv", row.names = F)




day6<- read.csv("WTDAY6-vs-KODAY6/SigDEGs.csv",header = T)
#day6<- day6[day6$ID %in% selected,]

NormalizedCounts = read.table("WTDAY6-vs-KODAY6/normalizedCounts.csv",header=T,sep=",",check.names=FALSE, quote = "\"")
colnames(NormalizedCounts)[1] = "ID"
NormalizedCounts$WT = (NormalizedCounts[,2]+NormalizedCounts[,3])/2
NormalizedCounts$KO = (NormalizedCounts[,4]+NormalizedCounts[,5])/2
NormalizedCounts = NormalizedCounts[,c(1,6:7)]

zscore <- function(x){
  z <- (x - mean(x))/sd(x)
  z
}
for (i in 2:length(NormalizedCounts[1,])){
  NormalizedCounts[,i] = zscore(NormalizedCounts[,i])
}

selected <- read.csv("selected.csv",header = T)
selectedDEGs <-selected$x

SelectedNormalizedCounts <- NormalizedCounts[NormalizedCounts$ID %in% selectedDEGs,]
rownames(SelectedNormalizedCounts) <- SelectedNormalizedCounts$ID
SelectedNormalizedCounts <- SelectedNormalizedCounts[,-1]
SelectedNormalizedCounts <- SelectedNormalizedCounts[selectedDEGs,]
SelectedNormalizedCounts <- SelectedNormalizedCounts[!is.na(SelectedNormalizedCounts$WT),]

annotation_col <-data.frame(c("Ectodermal", "Ectodermal","Ectodermal","Mensendodermal",
                              "Mensendodermal", "Mensendodermal","Mensendodermal","Mensendodermal",
                              "Trophectodermal", "Trophectodermal","Trophectodermal","Trophectodermal"))
colnames(annotation_col) <- c("Group")
row.names(annotation_col) <- row.names(SelectedNormalizedCounts)

png(file = "EBDay6.png", width = 850, height = 1600, units = "px", res = 300)
pheatmap(SelectedNormalizedCounts, border_color="black", treeheight_row=0, scale="row", 
         cluster_rows=F, show_rownames=T, show_colnames=T, cluster_cols=F, annotation_row=annotation_col,
         fontsize = 6, angle_col ="0", main="EB Day 6")
dev.off()
