library(vennerable)

setwd("C:/Users/Fei/Dropbox (EinsteinMed)/Fei/Qin-None-Stranded/")
directory <- getwd()

day0<- read.csv("WTDAY0-vs-KODAY0/SigDEGs.csv",header = T)
day0<- day0$ID

day3<- read.csv("WTDAY3-vs-KODAY3/SigDEGs.csv",header = T)
day3<- day3$ID

day6<- read.csv("WTDAY6-vs-KODAY6/SigDEGs.csv",header = T)
day6<- day6$ID

lday0 <- list(day0)
names(lday0) = " "
lday3 <- list(day3)
names(lday3) = "  "
lday6 <- list(day6)
names(lday6) = "   "

png(file = "Venn.png", width = 1600, height = 1600, units = "px", res = 300)
w <- Venn(Sets=c(lday0,lday3,lday6))
plot(w)
dev.off()

diffday6 <- setdiff(day6,day0)
diffday6 <- setdiff(diffday6,day3)
write.csv(diffday6, "diffday6.csv",row.names = F)
GOlist <- intersect(diffday6,day6)

write.csv(GOlist, "Golist.csv",row.names = F)

overrideTriple=T
venn.plot <- venn.diagram(
  x = list(I = day0,
           II = day3,
           III = day6),
  category.names = c(
    expression( bold('Day0') ),
    expression( bold('Day3') ),
    expression( bold('Day6') )
  ),
  filename = 'Fig3-1_triple_labels_sub_and_superscripts.tiff',
  output = TRUE,
  height = 3000,
  width = 3000,
  resolution = 300,
  compression = 'lzw',
  units = 'px',
  lwd = 1,
  lty = 'blank',
  fill = c('yellow', 'purple', 'green'),
  cex = 1,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.085, 0.085, 0.085),
  cat.fontfamily = "sans",
  scaled = F,
  euler.d = F,
  rotation = 1
)

#https://stackoverflow.com/questions/11727068/scaling-triple-venn-diagram-in-r-with-venndiagram-package
