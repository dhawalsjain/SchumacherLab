source("functions.R")
source("vars.R")


INPUT="A:/work/mm10_Annotations/repeats.bed"

if(!file.exists( paste0(DATADIR,"repeats.RData")) ){
  cat("generating repeat file for mm10\n")
  mm10 <- read.delim(INPUT,header = F,stringsAsFactors = F)
  names(mm10) <- c("chr","start","end","name","width","strand")
  mm10 <- as(mm10,"GRanges")
  mm10$id <- paste0(1:length(mm10))
  alu <- mm10[mm10$name%in% c("B1_Mm","B1_Mur1","B1_Mur2","B1_Mur3","B1_Mur4","B1_Mus1","B1_Mus2","B1F","B1F1", "B1F2","B2_Mm1a","B2_Mm1t","B2_Mm2","B3","B3A","B4","B4A","BC1_Mm")]
  mm10 <- subset(mm10,!mm10$id%in%alu$id)
  mm10 <- reduce(mm10,ignore.strand=T)
  alu <- reduce(alu,ignore.strand=T)
  save(mm10,alu,file=paste0(DATADIR,"repeats.RData"))
  rm(mm10,alu)
}
