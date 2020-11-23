setwd("C:/Users/mikef/Dropbox (EinsteinMed)/RNA-Seq/")
directory <- getwd()

day0<- read.csv("WTDAY0-vs-KODAY0/SigDEGs.csv",header = T)
day0<- day0$ID

day3<- read.csv("WTDAY3-vs-KODAY3/SigDEGs.csv",header = T)
day3<- day3$ID

day6<- read.csv("WTDAY6-vs-KODAY6/SigDEGs.csv",header = T)
day6<- day6$ID

library(VennDiagram)
#library(BioVenn)
#library(DescTools)
#library(gplots)
#library(nVennR)
#PlotVenn(x=list(day0, day3, day6), labels = list("day0", "day3", "day6"))
#biovenn <- draw.venn(day0, day3, day6, xtitle="Day0", ytitle="Day3", ztitle="Day6", t_c="#FFFFFF", subtitle="Venn Plot Day0, Day3, Day6", st_c="#FFFFFF", xt_c="#FFFFFF", yt_c="#FFFFFF", zt_c="#FFFFFF", nrtype="abs", nr_c="#FFFFFF", x_c="coral", y_c="coral", z_c="coral")

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
  lwd = 6,
  lty = 'blank',
  fill = c('yellow', 'purple', 'green'),
  cex = 3.5,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 3,
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
