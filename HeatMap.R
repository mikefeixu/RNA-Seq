library(pheatmap)
library(ComplexHeatmap)
library(RColorBrewer)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")

setwd("") # Path to your work directory.

directory <- getwd()

day0<- read.csv("WTDAY0-vs-KODAY0/SigDEGs.csv",header = T)
day0<- day0$ID

day3<- read.csv("WTDAY3-vs-KODAY3/SigDEGs.csv",header = T)
day3<- day3$ID

day6<- read.csv("WTDAY6-vs-KODAY6/SigDEGs.csv",header = T)
day6<- day6$ID

common <- union(union(day0, day3),day6)

write.csv(common, "common.csv", row.names = F)

common <- read.csv("common.csv",header = T)
commonDEGs <-common$x

Counts = read.table("WT-vs-KO/normalizedCounts.csv",header=T,sep=",",check.names=FALSE, quote = "\"")
colnames(Counts)[1] = "ID"

CommonTPMCounts <- Counts[Counts$ID %in% commonDEGs,]
rownames(CommonTPMCounts) <- CommonTPMCounts$ID
CommonTPMCounts <- CommonTPMCounts[,-1]

# Plot heatmap with package pheatmap
annotation_col <-data.frame(c("Day0", "Day0","Day0","Day0",
                              "Day3", "Day3","Day3","Day3",
                              "Day6", "Day6","Day6","Day6"),
                            c("WT-1", "WT-2","KO-1","KO-2",
                              "WT-1", "WT-2","KO-1","KO-2",
                              "WT-1", "WT-2","KO-1","KO-2"))
colnames(annotation_col) <- c("DAY","Group")
annotation_col <- annotation_col[,c(2,1)]
row.names(annotation_col) <- colnames(CommonTPMCounts)

png(file = "SigDEGs.png", width = 1650, height = 1600, units = "px", res = 300)
pheatmap(as.matrix(data.m), border_color=NA, treeheight_row=0, scale="row", 
         cluster_rows=T, show_rownames=F, show_colnames=T, cluster_cols=F, 
         annotation_col=annotation_col, use_raster=F,
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
         fontsize = 6, angle_col ="90", main="Heat Map of DEGs (Expression Z-Score)")
dev.off()

# Plot heatmap with package ComplexHeatmap
colnames(CommonTPMCounts) = c("WT-1", "WT-2","KO-1","KO-2",
                              "WT-1", "WT-2","KO-1","KO-2",
                              "WT-1", "WT-2","KO-1","KO-2")
da <- CommonTPMCounts
data <- apply(da,1,function(x){as.numeric(x)})
data <- t(data)
data.m<-apply(data,1,scale)
data.m<-t(data.m)
colnames(data.m)<- colnames(da)
rownames(data.m) <- rownames(da)

data.m<-data.m[rev(row.names(data.m)),]

ha0 = Heatmap(as.matrix(data.m),cluster_columns = FALSE,cluster_rows = TRUE,
              name = " Z score",
              row_dend_width = unit(12, "mm"), show_row_dend = FALSE,
              row_names_gp=gpar(fontsize = 10, col=c(rep("red",1), rep("black",51))), 
              row_title_gp = gpar(col = "black",fontsize = 8), show_row_names = FALSE,
              show_column_names = TRUE, column_names_gp = gpar(cex=1.2, fontsize=10,col= "black"),
              show_heatmap_legend = TRUE,column_names_side = "top",
              heatmap_legend_param=list(color_bar = "continuous", title_position = "topcenter",
                                        grid_height = unit(6, "mm"),grid_width = unit(4, "mm"),
                                        labels_gp = gpar(fontsize = 11,col="black")))

png(file = "SigDEGs_ComplexHeatmap.png", width = 1650, height = 1600, units = "px", res = 300)
ha0
dev.off()

# Plot heatmap for selected genes from normalized counts.
selected <- c("Pax6","Otx2","Fgf5","Nes",
              "Fgfr2","T","Sox17","Gata2","Gata4","Gata6",
            "Cdx2","Eomes","Krt7","Egfr")

write.csv(selected, "selected.csv", row.names = F)


#NormalizedCounts data processing
NormalizedCounts = read.table("normalizedCounts.csv",header=T,sep=",",check.names=FALSE, quote = "\"")
colnames(NormalizedCounts)[1] = "ID"
NormalizedCounts$WTDay0 = (NormalizedCounts[,2]+NormalizedCounts[,3])/2
NormalizedCounts$KODay0 = (NormalizedCounts[,4]+NormalizedCounts[,5])/2
NormalizedCounts$WTDay3 = (NormalizedCounts[,6]+NormalizedCounts[,7])/2
NormalizedCounts$KODay3 = (NormalizedCounts[,8]+NormalizedCounts[,9])/2
NormalizedCounts$WTDay6 = (NormalizedCounts[,10]+NormalizedCounts[,11])/2
NormalizedCounts$KODay6 = (NormalizedCounts[,12]+NormalizedCounts[,13])/2
NormalizedCounts = NormalizedCounts[,c(1,14:19)]

selected <- read.csv("selected.csv",header = T)
selectedDEGs <-selected$x

SelectedNormalizedCounts <- NormalizedCounts[NormalizedCounts$ID %in% selectedDEGs,]
rownames(SelectedNormalizedCounts) <- SelectedNormalizedCounts$ID
SelectedNormalizedCounts <- SelectedNormalizedCounts[,-1]
SelectedNormalizedCounts <- SelectedNormalizedCounts[selectedDEGs,]
SelectedNormalizedCounts <- SelectedNormalizedCounts[!is.na(SelectedNormalizedCounts$WTDay6),]
write.table(SelectedNormalizedCounts, "Pluripotency_genes.txt", row.names = T, sep="\t")

annotation_col <-data.frame(c("Ectodermal", "Ectodermal", "Ectodermal", "Ectodermal", "Mensendodermal", "Mensendodermal",
                              "Mensendodermal", "Mensendodermal", "Mensendodermal", "Mensendodermal",
                              "Trophectodermal", "Trophectodermal", "Trophectodermal", "Trophectodermal"))
colnames(annotation_col) <- c("Group")
row.names(annotation_col) <- row.names(SelectedNormalizedCounts)
write.table(SelectedNormalizedCounts, "Pluripotency_genes.txt", row.names = T, sep="\t")


# Plot heatmap for selected genes with package pheatmap
zscore <- function(x){
  z <- (x - mean(x))/sd(x)
  z
}
for (i in 1:length(SelectedNormalizedCounts[1,])){
  SelectedNormalizedCounts[,i] = zscore(SelectedNormalizedCounts[,i])
}

png(file = "EBDay6.png", width = 1600, height = 1600, units = "px", res = 300)
pheatmap(as.matrix(data.m), border_color=0, treeheight_row=0, scale = "row",
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breaksList)),
         breaks = breaksList,
         cluster_rows=F, show_rownames=T, show_colnames=T, cluster_cols=F, #annotation_row=annotation_col,
         fontsize = 6, angle_col ="0", main=" ")
dev.off()


# Plot heatmap for selected genes with package ComplexHeatmap
da0=read.table(file="Pluripotency_genes.txt",header=T,check.name=F,sep="\t",row.names=1)
da0=da0[,c(3,4,5,6)]
attach(da0)
da <-da0
rownames(da) <- rownames(da0)
data <- apply(da,1,function(x){as.numeric(x)})
data <- t(data)
data.m<-apply(data,1,scale)
data.m<-t(data.m)
colnames(data.m)<- colnames(da)
rownames(data.m) <- rownames(da)

ha1 = Heatmap(as.matrix(data.m),cluster_columns = FALSE,cluster_rows = TRUE,
              name = "Z score", clustering_method_rows = "average", clustering_distance_rows = "pearson",
              row_dend_width = unit(12, "mm"),show_row_dend = FALSE,
              row_names_gp=gpar(fontsize = 10,col=c(rep("red",1),rep("black",51))), 
              row_title_gp = gpar(col = "black",fontsize = 8), show_row_names = FALSE,
              show_column_names = TRUE, column_names_gp = gpar(cex=1.2, fontsize=10,col= "black"),
              show_heatmap_legend = TRUE,column_names_side = "top",
              heatmap_legend_param=list(color_bar = "continuous", title_position = "topcenter",
                                        grid_height = unit(6, "mm"),grid_width = unit(4, "mm"),
                                        labels_gp = gpar(fontsize = 11,col="black")))

png(file = "EB6.png", width = 1650, height = 1600, units = "px", res = 300)
ha1
dev.off()
