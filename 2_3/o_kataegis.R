rm(list = ls())
source("functions.R")
source("vars.R")

########################################################################################
## this script plots intermutation distances for high confidence sites 


if(F){
  get.kataegis.dataframe <- function(df1){
    require(BSgenome.Mmusculus.UCSC.mm10)
    chrs <- data.frame(chr=seqnames(BSgenome.Mmusculus.UCSC.mm10),length=seqlengths(BSgenome.Mmusculus.UCSC.mm10))
    chrs <- chrs[-grep("_",chrs$chr),]
    chrs <- chrs[-grep("chrM",chrs$chr),]
    chrs$cumlen <- cumsum(as.numeric(chrs$length))
    chrs$cumlen <- c(0,cumsum(as.numeric(chrs$length))[1:20])
    chrs$chr <- gsub("chr","",chrs$chr)
    
    cf <- get.distance2next(df1)
    cf <- merge(cf,chrs[,c(1,3)],by="chr",all.x=T)
    cf$pos <- cf$start+cf$cumlen
    return(cf)
  }
  require(BSgenome.Mmusculus.UCSC.mm10)
  chrs <- data.frame(chr=seqnames(BSgenome.Mmusculus.UCSC.mm10),length=seqlengths(BSgenome.Mmusculus.UCSC.mm10))
  chrs <- chrs[-grep("_",chrs$chr),]
  chrs <- chrs[-grep("chrM",chrs$chr),]
  chrs$cumlen <- cumsum(as.numeric(chrs$length))
  chrs$cumlen <- c(0,cumsum(as.numeric(chrs$length))[1:20])
  chrs$chr <- gsub("chr","",chrs$chr)
  
  cf <- c()
  load(paste0(DATADIR,"EditingLevelReport_Anno.RData"))
  load(paste0(DATADIR,"DRE_DE_matrix.RData"))
  res <- res[,c("idx","biotype","gene")] %>% unique
  res <- res[!is.na(res$gene),]
  res <- data.table(res)
  funx <- function(x){
    paste0(unique(x),collapse = ",")
  }
  res <- res[,lapply(.SD, funx), by=list(idx)]
  res <- as.data.frame(res)
  
  df1 <- hc_edits_anno[hc_edits_anno$idx%in%hlpr[hlpr$control>0,]$idx,]
  df1 <- as(df1,"GRanges")
  df1 <- unique(df1[,c('idx','repeats','snpeff_uq',"known")]) %>% as.data.frame
  names(df1)[1] <- "chr"
  df1 <- merge(df1,res,by="idx",all.x=T)
  cf <- rbind(cf,cbind(get.kataegis.dataframe(df1),profile="control"))
  
  df1 <- hc_edits_anno[hc_edits_anno$idx%in%hlpr[hlpr$D7>0,]$idx,]
  df1 <- as(df1,"GRanges")
  df1 <- unique(df1[,c('idx','repeats','snpeff_uq',"known")]) %>% as.data.frame
  names(df1)[1] <- "chr"
  df1 <- merge(df1,res,by="idx",all.x=T)
  cf <- rbind(cf,cbind(get.kataegis.dataframe(df1),profile="D7"))
  
  df1 <- hc_edits_anno[hc_edits_anno$idx%in%hlpr[hlpr$D14>0,]$idx,]
  df1 <- as(df1,"GRanges")
  df1 <- unique(df1[,c('idx','repeats','snpeff_uq',"known")]) %>% as.data.frame
  names(df1)[1] <- "chr"
  df1 <- merge(df1,res,by="idx",all.x=T)
  cf <- rbind(cf,cbind(get.kataegis.dataframe(df1),profile="D14"))
  
  ## gene-wise FPKM
  load(paste0(DATADIR,"DEExpression_DESeq2.RData"))
  load(paste0(DATADIR,"genes.RData"))
  genes <- unique(gGenes[,1:6])
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
  cnt$D7_FPKM <- round(log2(apply(cnt[,1:3],1,gmean)+1),2)
  cnt$D14_FPKM <- round(log2(apply(cnt[,4:6],1,gmean)+1),2)
  cnt$Control_FPKM <- round(log2(apply(cnt[,7:9],1,gmean)+1),2)
  cnt <- cnt[,10:13]
  cnt <- merge(cnt,unique(genes[,c("chr", "strand", "geneStart", "geneEnd","gene")]),by="gene")
  cnt$start <- cnt$end <- cnt$geneStart
  cnt$geneEnd <- cnt$geneStart <- NULL
  cnt <- cnt[,c("chr", "start","end", "strand","gene", "D7_FPKM", "D14_FPKM", "Control_FPKM")]
  cnt <- get.kataegis.dataframe(cnt)
  
  z <- imputed_edits
  z$id <- rownames(z)
  z$D7 <- round(apply(z[,1:3],1,gmean),2)
  z$D14 <- round(apply(z[,4:6],1,gmean),2)
  z$Control <- round(apply(z[,7:9],1,gmean),2)
  z <- z[,10:13]
  uq <- as(hc_edits_anno,"GRanges")
  uq <- unique(uq[,"id"])
  uq <- as.data.frame(uq)
  z <- merge(z,uq,by="id")
  names(z)[5] <- "chr"
  z <- get.kataegis.dataframe(z)
  cntLevels <- z
  save(cf,chrs,cnt,cntLevels,file=paste0(DATADIR,"EditingSites_Kataegis.RData"))
}