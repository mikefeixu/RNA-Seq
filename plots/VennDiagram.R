# Uncomment below lines to install Vennerable.
# Set project folder
# Set significant DEGs files

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("RBGL")
# BiocManager::install("graph")
# 
# install.packages("devtools");
# library(devtools);
# install_github("js229/Vennerable")

library(Vennerable)
setwd("<Your project directory>")
directory <- getwd()

WT<- read.csv("mESC_WT-vs-NPC_WT/SigDEGs.csv",header = T)
WT<- WT$ID

TKO<- read.csv("mESC_TKO-vs-NPC_TKO/SigDEGs.csv",header = T)
TKO<- TKO$ID

WT <- list(WT)
names(WT) = "WT"

TKO <- list(TKO)
names(TKO) = "TKO"


png(file = "Venn.png", width = 1600, height = 1600, units = "px", res = 300)
w <- Venn(Sets=c(WT,TKO))
plot(w)
dev.off()

pdf("Venn_WT_TKO.pdf")
plot(w)
dev.off()

mESC_Up<- read.table("mESC_WT-vs-mESC_TKO/SigDEGsUpList.tsv",header = F)
mESC_Up<- mESC_Up$V1

mESC_Down<- read.table("mESC_WT-vs-mESC_TKO/SigDEGsDownList.tsv",header = F)
mESC_Down<- mESC_Down$V1

NPC_Up<- read.table("NPC_WT-vs-NPC_TKO/SigDEGsUpList.tsv",header = F)
NPC_Up<- NPC_Up$V1

NPC_Down<- read.table("NPC_WT-vs-NPC_TKO/SigDEGsDownList.tsv",header = F)
NPC_Down<- NPC_Down$V1

mESC<- read.csv("mESC_WT-vs-mESC_TKO/SigDEGs.csv",header = T)
mESC<- mESC$ID

NPC<- read.csv("NPC_WT-vs-NPC_TKO/SigDEGs.csv",header = T)
NPC<- NPC$ID

NPC_Down_mESC <- setdiff(NPC_Down, mESC)
NPC_Up_mESC <- setdiff(NPC_Up, mESC)
NPC_and_mESC <- intersect(NPC, mESC)

write.csv(NPC_Down_mESC, "NPC_Down_mESC.csv",row.names = F)
write.csv(NPC_Up_mESC, "NPC_Up_mESC.csv",row.names = F)
write.csv(NPC_and_mESC, "NPC_and_mESC.csv",row.names = F)




