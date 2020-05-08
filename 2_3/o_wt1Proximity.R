rm(list = ls())
source("functions.R")
source("vars.R")

########################################################################################
## this script provides distance to nearest Wt1 peak

### liftover
if(F){
   
  a <- read.delim(paste0(RAWDIR,"wt1_mm9/D14_peaks.narrowPeak"),header=F,stringsAsFactors = F) 
  a <- a[,1:6]
  write.table(a,file=paste0(RAWDIR,"wt1_mm9/D14_peaks.bed"),sep="\t",quote = F,row.names = F,col.names = F)
   
  a <- read.delim(paste0(RAWDIR,"wt1_mm9/D9_peaks.narrowPeak"),header=F,stringsAsFactors = F) 
  a <- a[,1:6]
  write.table(a,file=paste0(RAWDIR,"wt1_mm9/D9_peaks.bed"),sep="\t",quote = F,row.names = F,col.names = F)
  a <- read.delim(paste0(RAWDIR,"wt1_mm9/WT_peaks.narrowPeak"),header=F,stringsAsFactors = F) 
  a <- a[,1:6]
  write.table(a,file=paste0(RAWDIR,"wt1_mm9/WT_peaks.bed"),sep="\t",quote = F,row.names = F,col.names = F)
  
   
   
}






