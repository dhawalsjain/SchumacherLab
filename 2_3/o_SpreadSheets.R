rm(list=ls())
source("functions.R")
source("vars.R")

###########################################3


if(F){
  load(paste0(DATADIR,"genes.RData"))
  load(paste0(DATADIR,"DEExpression_DESeq2.RData"))
  
  genes <- unique(gGenes[,c("chr", "strand", "gene", "geneStart", "geneEnd", "gene_id", "biotype")])
  genes$width <- abs(genes$geneStart-genes$geneEnd)
  cnt <- DESeq2::counts(outC[[2]],normalized=T) %>% as.data.frame
  names(cnt) <- unname(unlist(l1[names(cnt)]))
  cnt$id <- rownames(cnt)
  cnt <- merge(cnt,genes[,c("gene","gene_id","width")],by.x="id",by.y="gene_id")
  cnt <- unique(cnt)
  rownames(cnt) <- cnt$id
  cnt$id <- NULL
  cnt$width <-cnt$width/1000
  cnt[,1:9] <- round(cnt[,1:9]/cnt$width,2)
  cnt$width <- NULL
  cnt[,1:9] <- log2(cnt[,1:9]+1)
  cnt[,1:9] <- round(cnt[,1:9],2)
  cnt <- merge(cnt,genes,by="gene",all.x=T)
  
  ## FPKM
  cnt <- cnt[,c("chr", "strand", "geneStart", 
                "geneEnd", "gene_id", "biotype", "gene", "D7_R1", "D7_R2", "D7_R3", "D14_R1", "D14_R2", "D14_R3", 
                "control_R1", "control_R2", "control_R3")]
  
  ## DE genes D7/Control
  d7 <- out1[[1]]
  d7$baseMean <- d7$lfcSE <- d7$stat <- NULL
  names(d7)[c(1,4)] <- c("ensemble_id","FDR")
  d7 <- d7[!is.na(d7$FDR),]
  
  
  ## DE genes D14/Control
  d7 <- out0[[1]]
  d7$baseMean <- d7$lfcSE <- d7$stat <- NULL
  names(d7)[c(1,4)] <- c("ensemble_id","FDR")
  d7 <- d7[!is.na(d7$FDR),]
  d7
  
    ## DE genes D14/Control
  d7 <- out1[[1]]
  d7$baseMean <- d7$lfcSE <- d7$stat <- NULL
  names(d7)[c(1,4)] <- c("ensemble_id","FDR")
  d7 <- d7[!is.na(d7$FDR),]
  d7
  
  
  ## Editing levels
  load(paste0(DATADIR,"EditingLevelReport_Anno.RData"))
  sps <- imputed_edits
  sps$id <- rownames(sps)
  sps <- merge(sps,unique(hc_edits_anno[,c("id","V4","known","repeats","snpeff","aaswap")]))
  names(sps)[11:15] <- c("Edit","is_known","repeat","snpeff_annotation","AA_change")
  sps$in_Intron <- sps$in_3UTR <- sps$in_5UTR <- sps$in_spliceRegion <- sps$is_missense <- sps$is_intergenic <- sps$synony <- "no"
  sps[grep("intron",sps$snpeff_annotation),]$in_Intron <- as.character("yes")
  sps[grep("3_prime_UTR",sps$snpeff_annotation),]$in_3UTR <- as.character("yes")
  sps[grep("5_prime_UTR",sps$snpeff_annotation),]$in_5UTR <- as.character("yes")
  sps[grep("intergenic_region",sps$snpeff_annotation),]$is_intergenic <- as.character("yes")
  sps[grep("missense",sps$snpeff_annotation),]$is_missense <- as.character("yes")
  sps[grep("splice_",sps$snpeff_annotation),]$in_spliceRegion <- as.character("yes")
  sps[grep("synonymous",sps$snpeff_annotation),]$synony <- as.character("yes")
  
  
  
  ## DRE
  load(paste0(DATADIR,"DRE_DE_matrix.RData"))
  cn <- read.delim(paste0(DATADIR,"SnpEffector.txt"),header=T,stringsAsFactors = F)
  cn$Location <- gsub("-\\S*","",cn$Location)
  cn$REF <- ifelse(cn$Allele=="C",as.character("T"),as.character("G"))
  cn$idx <- paste0(cn$Location,":",cn$REF,">",cn$Allele)
  cn <- cn[,c("Consequence","SIFT","idx")] %>% unique
  cn$SIFT <- gsub("\\(\\S*\\)","",cn$SIFT)
  cn$SIFT <- gsub("deleterious_low_confidence|tolerated_low_confidence","tolerated",cn$SIFT)
  cn <- cn[,c("idx","SIFT")] %>% unique
  cn <- cn[cn$SIFT!="-",]
  res <- merge(res,cn,by="idx",all.x=T)
  ## D7/Control
  r1 <- res[res$adj.P.Val<0.05 & res$comparison=="D7_control",]
  r1 <- r1[,c("idx", "logFC","change", "comparison", "adj.P.Val",
              "known", "repeats","clustered","snpeff_uq",  
              "aaswap", "SIFT","control", "D7","D14", 
              "tracking_id", "gene", "expn_l2fc", "expn_fdr","biotype","hg19" )] %>% unique

  ## D14/Control
  r2 <- res[res$adj.P.Val<0.05 & res$comparison=="D14_control",]
  r2 <- r2[,c("idx", "logFC","change", "comparison", "adj.P.Val",
              "known", "repeats","clustered","snpeff_uq",  
              "aaswap", "SIFT","control", "D7","D14", 
              "tracking_id", "gene", "expn_l2fc", "expn_fdr","biotype","hg19" )] %>% unique
  
  ##
  load(paste0(DATADIR,"miRSeedOverlaps_filtermiRNAs.RData"))
  mirSeedOlp$names <- gsub("_\\D*","",mirSeedOlp$names)
  mirSeedOlp$check <- "yes"
  mirSeedOlp <- reshape2::dcast(mirSeedOlp,idx+names+mir~seq,value.var = "check")
  
  r1 <- merge(r1, unique(mirSeedOlp[,c("idx","mir","alt","ref")]), by="idx",all.x=T )
  r2 <- merge(r2, unique(mirSeedOlp[,c("idx","mir","alt","ref")]), by="idx",all.x=T )
  r1Seed <- r1
  r2Seed <- r2
  save(r1Seed,r2Seed,file=paste0(DATADIR,"miRNA-Seed_Analyses_TrxLevelInfo.RData"))
  
  load(paste0(DATADIR,"miRNA-Seed_Analyses_TrxLevelInfo.RData"))
  r1Seed <- r1Seed[!is.na(r1Seed$alt)  !is.na(r1Seed$ref),]
}










