rm(list = ls())
source("functions.R")
source("vars.R")

###############################################################################################
### Editing and chromatin states

## 1. Use the chrome HMM states (15 state model) for mouse kidney from UCSC (mm10)
##    https://genome.ucsc.edu/cgi-bin/hgTables
## 2. E14, E15, E16 and P0 data
## 3. Focus on regions that have consistent chromotin state (i.e. same state across all 4). 
##    these are likely kidney-specific sequences with a given chromatin state
## 4. Chromatin states     
##      State 1 -       Dark Green - Promoter, Active (Pr-A)
##     State 2 -       Light Green - Promoter, Weak (Pr-W)
##     State 3 -       Light Grey - Promoter, Bivalent (Pr-B)
##     State 4 -       Green - Promoter, Flanking Region (Pr-F)
##     State 5 -       Bright Yellow - Enhancer, Strong TSS-distal (En-Sd)
##     State 6 -       Bright Yellow - Enhancer, Strong TSS-proximal (En-Sp)
##     State 7 -       Light Yellow - Enhancer, Weak (En-W)
##     State 8 -       Dark Grey - Enhancer, Poised TSS-distal (En-Pd)
##     State 9 -       Dark Grey - Enhancer, Poised TSS-proximal (En-Pp)
##     State 10 -       Dark Blue - Transcription, Strong (Tr-S)
##     State 11 -       Royal Blue - Transcription, Permissive (Tr-P)
##     State 12 -       Light Blue - Transcription, Initiation (Tr-I)
##     State 13 -       State 13 - Salmon - Heterochromatin, Polycomb-associated (Hc-P)
##     State 14 -       Pink - Heterochromatin, H3K9me3-associated (Hc-H)
##     State 15 -       White - No significant signal (Ns)
HMMCOLS = c("-"="black", "En-Pd"="gray", "En-Sp"="yellow3", "Hc-P"="salmon", "NS"="gray80", 
    "Pr-A"="green4", "Pr-F"="green", "Pr-W"="lightgreen", 
    "Tr-I"="lightblue", "Tr-P"="royalblue", "Tr-S"="royalblue4")


if(F){
  load(paste0(DATADIR,"EditingLevelReport_Anno.RData"))
  z <- imputed_edits
  hlpr$id <- gsub(":T>C|:A>G","",hlpr$idx)
  z <- z[match(hlpr$id,rownames(z)),]
  hlpr$cHMM_Final <- ifelse(is.na(hlpr$cHMM_Final),as.character('-'),as.character(hlpr$cHMM_Final))
  
  z$HMM <-hlpr$cHMM_Final
  z <- z[z$HMM!="-",]
  partition = as.factor(z$HMM)
  table(partition)
  rr <-hlpr[hlpr$id%in%rownames(z),]
  rr <- merge(rr,unique(hc_edits_anno[,c("id","ConsScore","hg19")]),by="id",all.x=T )
  rr <- rr[match(rownames(z),rr$id),]
  rr <- cbind(rr,z[,1:9])
  rr <- rr[,!names(rr)%in%c("id","cHMM_P0","cHMM_E16","cHMM_E15","cHMM_E14")]
  
  rownames(z) <- NULL
  col = colorRamp2(c(0, max(z[,1:9])), c("white","blue4"))
  Heatmap(z[,1:9],column_title = "Editing Dynamics\nby chromatin state",name = "scale",col = col,cluster_columns = F,
          split = partition, row_title_rot = 0,
          row_title_gp = gpar(col = HMMCOLS,fontsize=16,fontface="bold"))
  
  
  z <- cbind(z,rr[,c("idx","known","repeats","snpeff_uq","aaswap","ConsScore","hg19")])
  z <- z[,c("idx", "HMM", "known", "repeats","snpeff_uq", "aaswap", "ConsScore", "hg19",
            "D7_R1", "D7_R2", "D7_R3", "D14_R1", "D14_R2", "D14_R3", "control_R3", 
            "control_R2", "control_R1")]
  
}










