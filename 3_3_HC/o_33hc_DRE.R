rm(list = ls())
source("functions.R")
source("vars.R")



###############################################################################################
### Differential RNA editing  data frame
###############################################################################################
## lot of manual entries
if(F){
  
  load(paste0(DATADIR,"33hc_EditingLevelReport_Anno.RData"))
  
  z <- imputed_edits
  z$D7 <- apply(z[,1:3],1,gmean)
  z$D14 <- apply(z[,4:6],1,gmean)
  z$control <- apply(z[,7:9],1,gmean)
  z$id <- rownames(z)
  design <- data.frame(WT=rep(1,6),MUvsWT= c(0,0,0,1,1,1))
  rownames(design) <- paste0("Array",1:nrow(design))
  
  mydeedit <- function(zz,design){
    zz <- zz[rowSums(zz)>0,]
    fit <- limma::lmFit(zz,design = design)
    fit <- limma::eBayes(fit,trend = F)
    cf <- limma::topTable(fit, coef="MUvsWT", adjust="BH",number = nrow(z))
    cf$id <- rownames(cf)
    cf$logFC <- cf$logFC*(-1)
    cf
  }
  
  dre <- c()
  names(z)[c(1:3,7:9)]
  cf <- mydeedit(z[,c(1:3,7:9)],design)
  cf <- merge(cf,z[,c("D14","D7","control","id")],by="id")
  cf$change <- cf$D7 - cf$control
  cf$comparison <- "D7_control"
  dre <- rbind(dre,cf)
  
  names(z)[c(4:9)]
  cf <- mydeedit(z[,c(4:9)],design)
  cf <- merge(cf,z[,c("D14","D7","control","id")],by="id")
  cf$change <- cf$D14 - cf$control
  cf$comparison <- "D14_control"
  dre <- rbind(dre,cf)
  rm(cf)
  
  hlpr$id <- gsub(":A>G|:T>C","",hlpr$idx)
  names(hlpr)[6:8] <- paste0("is_",names(hlpr)[6:8])
  dre <- merge(dre,hlpr,by="id")
  dre$chr <- gsub(":\\S*","",dre$id) %>% as.character
  dre$start <- dre$end <- gsub("\\S*:","",dre$id) %>% as.numeric
  dre <- as(dre,"GRanges")
  
  d2 <- reduce(resize(dre[dre$comparison=="D7_control"],101,"center"))
  d2 <- resize(d2,width(d2)-100,"center")
  d2 <- d2[width(d2)>50]
  d1 <- reduce(resize(dre[dre$comparison=="D14_control"],101,"center"))
  d1 <- resize(d1,width(d1)-100,"center")
  d1 <- d1[width(d1)>50]
  
  x <- dre[dre$comparison=="D14_control"]
  x$clustered <- countOverlaps(x,d1,ignore.strand=T)
  sum(x$clustered)
  y <- dre[dre$comparison=="D7_control"]
  y$clustered <- countOverlaps(y,d1,ignore.strand=T)
  sum(y$clustered)
  dre <- c(x,y)
  sum(dre$clustered)
  rm(d2,d1,x,y)
  rm(design,z)
  
  
  load(paste0(DATADIR,"DEExpression_DESeq2.RData"))
  names(out0[[1]])[c(3,7)] <- names(out1[[1]])[c(3,7)] <- c("expn_l2fc","expn_fdr")
  expn <- rbind(cbind(out0[[1]][,c(1,3,7:12)],comparison="D7_control"),
                cbind(out1[[1]][,c(1,3,7:12)],comparison="D14_control"))
  expn <- with(expn,GRanges(chr,IRanges(geneStart,geneEnd),strand,id,gene,expn_l2fc,expn_fdr,comparison))
  rm(out0,out1,outC)

  mycorrl <- function(dre,h){
    d1 <- dre[dre$comparison==h,]
    e1 <- expn[expn$comparison==h]
    e1$expn_fdr <- ifelse(is.na(e1$expn_fdr),1,e1$expn_fdr)
    o1 <- findOverlaps(d1,e1,ignore.strand=T) %>% as.data.frame
    o1$expn_fdr <- e1[o1$subjectHits]$expn_fdr
    o1 <- data.table(o1)
    o1 <- o1[,minv:=min(expn_fdr,na.rm=T),by=list(queryHits)]
    o1 <- data.frame(o1)
    o1 <- o1[o1$minv==o1$expn_fdr,]
    
    d1$tracking_id <- d1$gene <- d1$expn_l2fc <- d1$expn_fdr <- NA
    d1[o1$queryHits]$tracking_id <- e1[o1$subjectHits]$id
    d1[o1$queryHits]$gene <- e1[o1$subjectHits]$gene
    d1[o1$queryHits]$expn_l2fc <- e1[o1$subjectHits]$expn_l2fc
    d1[o1$queryHits]$expn_fdr <- e1[o1$subjectHits]$expn_fdr
    d1 <- as.data.frame(d1)
    return(d1)
  }
  res <- mycorrl(dre,h="D7_control")
  res <- rbind(res,mycorrl(dre,"D14_control"))
  
  load(paste0(DATADIR,"genes.RData"))
  res <- merge(res,unique(gGenes[,c("gene_id","biotype")]),by.x="tracking_id",by.y="gene_id",all.x=T)
  
  uq <- unique(hc_edits_anno[,c("idx","hg19")])
  res <- merge(res,uq, by="idx",all.x=T)
  
  save(res,file=paste0(DATADIR,"33hc_DRE_DE_matrix.RData"))
  rm(z,design,cf,dre,hlpr,d1,d2,x,y,expn,res)
  
  
}

#### GO enrichments
if(F){
  load(paste0(DATADIR,"DEExpression_DESeq2.RData"))
  ggg <- out0[[1]][,c("id","geneStart","geneEnd")]
  ggg$length <- abs(ggg$geneStart-ggg$geneEnd)
  ggg <- ggg[,c("id","length")] %>% unique
  cnt <- DESeq2::counts(outC[[2]],normalized=T) %>% as.data.frame
  cnt$id <- rownames(cnt)
  cnt <- merge(cnt,ggg,by="id",all.x=T)
  rownames(cnt) <- cnt$id
  cnt$id <- NULL
  cnt[,1:9] <- cnt[,1:9]*1000/cnt$length
  cnt$length <- NULL
  rm(out0,out1,outC,ggg)
  cnt <- apply(cnt,2,function(x){round(x*1e6/sum(x))})
  cnt <- as.data.frame(cnt)
  cnt <- cnt[rowSums(cnt)>0,]
  
  load(paste0(DATADIR,"genes.RData"))
  exp.trx <- gGenes[gGenes$gene_id%in%rownames(cnt),]
  exp.trx <- exp.trx$ucscID %>% unique
  rm(cnt)
  
  load(paste0(DATADIR,"33hc_DRE_DE_matrix.RData"))
  res <- res[res$adj.P.Val<0.05 & res$clustered>0,]
  out <- c()
  x <- res[res$comparison=="D7_control",]$tracking_id %>% unique
  x <- gGenes[gGenes$gene_id%in%x,]$ucscID
  y <- res[res$comparison=="D14_control",]$tracking_id %>% unique
  y <- gGenes[gGenes$gene_id%in%y,]$ucscID
  out <- rbind(
    cbind(mm9_GOStats.enrichment(x,p.val = 0.1,universe = exp.trx),comparison="D7_control",change="all"),
    cbind(mm9_GOStats.enrichment(y,p.val = 0.1,universe = exp.trx),comparison="D14_control",change="all")
  )
  for(f in c("up","dn")){
    if(f == "up"){
      x <- res[res$comparison=="D7_control" & res$logFC>0 ,]$tracking_id %>% unique
      y <- res[res$comparison=="D14_control" & res$logFC>0,]$tracking_id %>% unique
    }else{
      x <- res[res$comparison=="D7_control" & res$logFC<0 ,]$tracking_id %>% unique
      y <- res[res$comparison=="D14_control" & res$logFC<0,]$tracking_id %>% unique
    }
    x <- gGenes[gGenes$gene_id%in%x,]$ucscID
    y <- gGenes[gGenes$gene_id%in%y,]$ucscID
    out <- rbind(out,
      cbind(mm9_GOStats.enrichment(x,p.val = 0.1,universe = exp.trx),comparison="D7_control",change=f),
      cbind(mm9_GOStats.enrichment(y,p.val = 0.1,universe = exp.trx),comparison="D14_control",change=f)
    )
  }
  
  
  load(paste0(DATADIR,"33hc_DRE_DE_matrix.RData"))
  res <- res[res$adj.P.Val<0.05,]
  x <- res[res$comparison=="D7_control",]$tracking_id %>% unique
  y <- res[res$comparison=="D14_control",]$tracking_id %>% unique
  x <- gGenes[gGenes$gene_id%in%x,]$ucscID
  y <- gGenes[gGenes$gene_id%in%y,]$ucscID
  out1 <- rbind(
    cbind(mm9_GOStats.enrichment(x,p.val = 0.1,universe = exp.trx),comparison="D7_control",change="all"),
    cbind(mm9_GOStats.enrichment(y,p.val = 0.1,universe = exp.trx),comparison="D14_control",change="all")
  )
  for(f in c("up","dn")){
    if(f == "up"){
      x <- res[res$comparison=="D7_control" & res$logFC>0 ,]$tracking_id %>% unique
      y <- res[res$comparison=="D14_control" & res$logFC>0,]$tracking_id %>% unique
    }else{
      x <- res[res$comparison=="D7_control" & res$logFC<0 ,]$tracking_id %>% unique
      y <- res[res$comparison=="D14_control" & res$logFC<0,]$tracking_id %>% unique
    }
    x <- gGenes[gGenes$gene_id%in%x,]$ucscID
    y <- gGenes[gGenes$gene_id%in%y,]$ucscID
    out1 <- rbind(out1,
                 cbind(mm9_GOStats.enrichment(x,p.val = 0.1,universe = exp.trx),comparison="D7_control",change=f),
                 cbind(mm9_GOStats.enrichment(y,p.val = 0.1,universe = exp.trx),comparison="D14_control",change=f)
    )
  }
  
  
  save(out,out1,file=paste0(DATADIR,"33hc_DRE_DE_GOEnrichment.RData"))
  
}

#### GO enrichments by genomic features
if(F){
  load(paste0(DATADIR,"DEExpression_DESeq2.RData"))
  ggg <- out0[[1]][,c("id","geneStart","geneEnd")]
  ggg$length <- abs(ggg$geneStart-ggg$geneEnd)
  ggg <- ggg[,c("id","length")] %>% unique
  cnt <- DESeq2::counts(outC[[2]],normalized=T) %>% as.data.frame
  cnt$id <- rownames(cnt)
  cnt <- merge(cnt,ggg,by="id",all.x=T)
  rownames(cnt) <- cnt$id
  cnt$id <- NULL
  cnt[,1:9] <- cnt[,1:9]*1000/cnt$length
  cnt$length <- NULL
  rm(out0,out1,outC,ggg)
  cnt <- apply(cnt,2,function(x){round(x*1e6/sum(x))})
  cnt <- as.data.frame(cnt)
  cnt <- cnt[rowSums(cnt)>0,]
  
  load(paste0(DATADIR,"genes.RData"))
  exp.trx <- gGenes[gGenes$gene_id%in%rownames(cnt),]
  exp.trx <- exp.trx$ucscID %>% unique
  rm(cnt)
  
  load(paste0(DATADIR,"33hc_DRE_DE_matrix.RData"))
  res <- res[res$adj.P.Val<0.05,]
  
  out <- c()
  for(feat in c("5_prime_UTR","missense","splice","synonymous")  ){ levels(as.factor(res$snpeff_uq))
    for(f in c("up","dn")){
      cat(feat,"\t",f,"..\t")
      if(f == "up"){
        x <- res[res$comparison=="D7_control" & res$logFC>0 & res$snpeff_uq==feat,]$tracking_id %>% unique
        y <- res[res$comparison=="D14_control" & res$logFC>0 & res$snpeff_uq==feat,]$tracking_id %>% unique
      }else{
        x <- res[res$comparison=="D7_control" & res$logFC<0 & res$snpeff_uq==feat,]$tracking_id %>% unique
        y <- res[res$comparison=="D14_control" & res$logFC<0 & res$snpeff_uq==feat,]$tracking_id %>% unique
      }
      x <- gGenes[gGenes$gene_id%in%x,]$ucscID
      y <- gGenes[gGenes$gene_id%in%y,]$ucscID
      cat(length(x)," ", length(y), "\n")
      if(length(x)>10){
        out <- rbind(out,
                     cbind(mm9_GOStats.enrichment(x,p.val = 0.1,universe = exp.trx),comparison="D7_control",change=f,genomeAnno=feat))
      }
      if(length(y)>10){
        out <- rbind(out,
                     cbind(mm9_GOStats.enrichment(y,p.val = 0.1,universe = exp.trx),comparison="D14_control",change=f,genomeAnno=feat))
      }
      
    }
    rm(feat,f,x,y)
  }
  out4 <- out
  save(out4,file=paste0(DATADIR,"33hc_DRE_DE_GOEnrichment_byGenomicfeatures.RData"))
  
}




