setwd("C:/Users/mikef/Dropbox (EinsteinMed)/RNA-Seq/")
directory <- getwd()

GO<- read.csv("WTDAY6-vs-KODAY6/GO.tsv", sep = "\t", header = T)
GO$EnrichmentFactor = GO$Count/GO$PopHits
GO$LogTransform = (-1)*log10(GO$FDR)
#GO <- GO[-18,]

library(ggplot2)
ggplot(GO, aes(x = EnrichmentFactor, y = Term, label = "")) +
  geom_point(aes(size = Count, colour = LogTransform),alpha=.2) +
  scale_x_continuous(breaks = seq(0.1, 0.4, 0.01)) +
  scale_y_discrete(limits=rev(GO$Term)) +
  scale_colour_continuous(low = "red", high = "yellow") +
  labs(colour = "-log10(FDR)", size = "Count") +
  #ggtitle("Selected GO Terms") +
  labs(x = "Enrichment Factor", y = "GO Terms") +
  geom_text(hjust = 3, size = 2) +
  scale_size(range = c(1,10)) +
  theme_bw()
