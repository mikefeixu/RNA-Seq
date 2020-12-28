setwd("C:/Users/Fei/Dropbox\ (EinsteinMed)/Fei/ChIP")
directory <- getwd()

GO<- read.csv("GO.tsv", sep = "\t", header = T)
GO$EnrichmentFactor = GO$Count/GO$PopHits
GO$LogTransform = (-1)*log10(GO$PValue)
#GO <- GO[-18,]

library(ggplot2)
png(file = "GO.png", width = 2000, height = 1600, units = "px", res = 300)
ggplot(GO, aes(x = EnrichmentFactor, y = Term, label = "")) +
  geom_point(aes(size = Count, colour = LogTransform),alpha=.99) +
  scale_x_continuous(breaks = seq(0.05, 0.3, 0.05)) +
  scale_y_discrete(limits=rev(GO$Term)) +
  #scale_colour_gradient(low = "deeppink", high = "yellow") +
  scale_colour_gradient(low = "blue", high = "red") +
  #scale_colour_gradient(low = "white", high = "black") +
  labs(colour = "-log10(pvalue)", size = "Count") +
  #ggtitle("Selected GO Terms") +
  labs(x = "Enrichment Factor", y = "GO Terms") +
  geom_text(hjust = 3, size = 2) +
  scale_size(range = c(1,10)) +
  theme_bw()
dev.off()
