# Uncomment below lines to install pheatmap/ComplexHeatmap.
# Set project folder
# Set significant DEGs files
# Specify the normalizedCounts.csv for the comparison
# Define annotation column
# Copy the selectedGenes.tsv and update to your selected genes for small heatmap. Comment accordingly line if all Genes are included

library(pheatmap)
library(ComplexHeatmap) # ComplexHeatmap is used here for beautify the legends. It will replace pheatmap in the future.
library(RColorBrewer)

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("pheatmap")
# BiocManager::install("ComplexHeatmap")

setwd("<Your project directory>")

# Get all DEGs list from different comparisons to create a combined heatmap for all samples.
Counts = read.table("WT-vs-TKO/normalizedCounts.csv",header=T,sep=",",check.names=FALSE, quote = "\"")
colnames(Counts)[1] = "ID"

# Plot heatmap for selected genes from normalized counts.
selected <- read.csv("selectedGenes.tsv",header = F, sep = "\n")
selected <- selected$V1
write.csv(selected, "selected.csv", row.names = F)
selected <- read.csv("selected.csv",header = T)
selectedDEGs <-selected$x

# Filter Counts with selected genes. Comment this line if all genes are included.
SigDEGCounts <- Counts[Counts$ID %in% selectedDEGs,]

rownames(SigDEGCounts) <- SigDEGCounts$ID
SigDEGCounts <- SigDEGCounts[,-1]

# Plot heatmap with package pheatmap
# Define the annotation table for columns of SigDEGcounts: 
annotation_col <-data.frame(c("WT", "WT","WT",
                              "TKO", "TKO"),
                            c("WT1", "WT2","WT3",
                              "TKOc2", "TKOh7"))
colnames(annotation_col) <- c("Genotype","Group")
annotation_col <- annotation_col[,c(2,1)]
row.names(annotation_col) <- colnames(SigDEGCounts)
annotation_col <- annotation_col[,-c(1), drop = FALSE]

pdf("SigDEGs_HeatMap_Selected_mESC.pdf")
pheatmap(as.matrix(SigDEGCounts), border_color=NA, treeheight_row=0, scale="row",
         cluster_rows=T, show_rownames=T, show_colnames=T, cluster_cols=T,
         annotation_col=annotation_col, use_raster=F,
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
         fontsize = 6, angle_col ="90", main="Heat Map of DEGs (Expression Z-Score)")
dev.off()
