# Adjust the commented lines to tune the plots to an acceptable outlook. 
# Manual set the x-axis breaks/range by setting scale_x_continuous/xlim
# Set project folder
# Copy Selected_GO.tsv to GO folder, update with your selected GO Terms found from 
# DAVID https://david.ncifcrf.gov/home.jsp
library(ggplot2)
setwd("Path to your project folder")
selected_GO_terms <- "GO/Selected_GO.tsv"

GO <- read.csv("GO/Selected_GO.tsv", sep = "\t", header = T)
GO$Term <- gsub("^GO(.*)~","",GO$Term)
GO$EnrichmentFactor = GO$Count/GO$PopHits
GO$LogTransform = (-1)*log10(GO$PValue)

png(file = "Selected_GO.png", width = 2000, height = 1600, units = "px", res = 300)
ggplot(GO, aes(x = EnrichmentFactor, y = Term, label = "")) +
  geom_point(aes(size = Count, colour = LogTransform),alpha=.99) +
  #scale_x_continuous(breaks = seq(0.00, 0.20, 0.02)) +
  scale_y_discrete(limits=rev(GO$Term)) +
  scale_colour_gradient(low = "blue", high = "red") +
  labs(colour = "-log10(pvalue)", size = "Count") +
  labs(x = "Enrichment Factor", y = "GO Terms") +
  geom_text(hjust = 3, size = 2) +
  scale_size(range = c(1,10)) +
  #xlim(0.00, 0.15) +
  theme_bw()
dev.off()
