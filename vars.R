require(ggplot2)
library(limma)
library(reshape)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(GenomicRanges)
library(knitr)
library(ComplexHeatmap)
library(DESeq2)
library(gplots)
library(data.table)
library(RColorBrewer)
library(corrplot)
library(DT)
library(cowplot)
library(circlize)
library(UpSetR) 
library(heatmap3)
library(gridExtra)
library(grid)
library(viridis)
library(dplyr)

DATADIR="A:/work/Kreidberg lab/Valerie/PAPER/DF/"
FIGDIR="A:/work/Kreidberg lab/Valerie/PAPER/FIG/"
RAWDIR="A:/work/Kreidberg lab/Valerie/PAPER/RAW/"
REPDIR="A:/work/Kreidberg lab/Valerie/PAPER/REPORTS/"

###############################################################################################
### sample name/color specification etc
###############################################################################################
if(T){
  
  l <- list()
  l[["control"]] <- c("VS868","VS679","VS582")
  l[["D7"]] <- c("VS1247","VS1249","VS36")
  l[["D14"]] <- c("VS838","VS839","VS840")
  l1 <- list()
  l1[["VS868"]] <- "control_R1"
  l1[["VS679"]] <- "control_R2"
  l1[["VS582"]] <- "control_R3"
  l1[["VS1247"]] <- "D7_R1"
  l1[["VS1249"]] <- "D7_R2"
  l1[["VS36"]] <- "D7_R3"
  l1[["VS838"]] <- "D14_R1"
  l1[["VS839"]] <- "D14_R2"
  l1[["VS840"]] <- "D14_R3"
  
  cols = c("control_R1"="darkolivegreen2", "control_R2"="darkolivegreen3", "control_R3"="darkolivegreen4", 
           "D7_R1"="firebrick1", "D7_R2"="firebrick3", "D7_R3"="firebrick4", 
           "D14_R1"="royalblue1", "D14_R2"="royalblue3", "D14_R3"="royalblue4",
           "VS868"="darkolivegreen2", "VS679"="darkolivegreen3", "VS582"="darkolivegreen4", 
           "VS1247"="firebrick1", "VS1249"="firebrick3", "VS36"="firebrick4", 
           "VS838"="royalblue1", "VS839"="royalblue3", "VS840"="royalblue4",
           "Control" = "darkolivegreen2","Day-7"="firebrick2","Day-14"="royalblue2",
           "control" = "darkolivegreen2","D7"="firebrick2","D14"="royalblue2",
           "Background" = "gray20")
  
  HMMCOLS = c("-"="black", "En-Pd"="gray", "En-Sp"="yellow3", "Hc-P"="salmon", "NS"="gray80", 
              "Pr-A"="green4", "Pr-F"="green", "Pr-W"="lightgreen", 
              "Tr-I"="lightblue", "Tr-P"="royalblue", "Tr-S"="royalblue4")
  
  libs <- data.frame(profile=c("VS1247", "VS1249", "VS36", "VS582", "VS679", "VS838",
                               "VS839", "VS840", "VS868"),
                     depths = c(153.1,160.3,45.2,185.6,156.5,31.8,39.4,48.5,157))
  libs$cov <- round((libs$depths*1e6*2*150)/2.5e9)
  libs$cov <- ifelse(libs$cov>=10,10,libs$cov)
  libs$hyperEdits <- c(9160,9049,3122,11078,6945,4247,2059,4574,1347)
  libs$unmap <- c(141958608,189244679,42447498,154495450,115871309,24739233,
    32654705,47923857,126216303)
  
  sampl <- c("VS1247"="D7_R1","VS1249"="D7_R2","VS36"="D7_R3","VS838"="D14_R1",
             "VS839"="D14_R2","VS840"="D14_R3","VS868"="Control_R1",
             "VS679"="Control_R2","VS582"="Control_R3")
  
  CUTOFFLEN=1
  FDR_CUTOFF=0.05
  
  gg_aes <- theme(legend.position = "none",
                  axis.text = element_text(size=16,color="black"),
                  axis.title = element_text(size=16,color="black"),
                  legend.text = element_text(size=14,colour = "black"),
                  strip.text = element_text(size=16,colour = "black"),
                  plot.title = element_text(size=20,colour = "black",hjust = 0.5),
                  legend.title = element_blank())
  
}