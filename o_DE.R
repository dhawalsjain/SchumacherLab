rm(list=ls())
source("functions.R")
source("vars.R")


###############################################################################################
### Gene expression - various counts using HTSeq
###############################################################################################
WD=getwd()
setwd("A:/work/Kreidberg lab/Valerie/HTSeq")
if(!file.exists(paste0(DATADIR,"MiscCounts.RData"))){
  cntlist <- list()
  for(f in c("5utr","3utr","gene","bam")){
    cat(f,"\n")
    files <- Sys.glob(paste0("*.",f,".counts.txt"))
    cnt <- c()
    for(i in 1:length(files)){
      cat(files[i],"\n")
      d <- read.delim(files[i],header=F,stringsAsFactors = F)
      if(i==1){
        cnt <- d
      }else{
        cnt <- cbind(cnt,d$V2)
      }
      rm(d,i)
    }
    names(cnt)[2:ncol(cnt)] <- gsub(paste0(".",f,".counts.txt"),"",files)
    cntlist[[f]] <- cnt
    cat("\n")
  }
  save(cntlist,file=paste0(DATADIR,"MiscCounts.RData"))
}
setwd(WD)
## Note: bowtie2, -q 30
## Note: bowtie2, -q 0

###############################################################################################
### Gene expression - various counts using HTSeq (Sandrine's old data)
###############################################################################################
WD=getwd()
setwd("A:/work/Kreidberg lab/Sep2018_RNASeq/HTSeq_bt/")
if(!file.exists(paste0(DATADIR,"MiscCounts_WT1.RData"))){
  cntlist <- list()
  for(f in c("gene","tx")){
    cat(f,"\n")
    files <- Sys.glob(paste0("*.",f,"counts.txt"))
    cnt <- c()
    for(i in 1:length(files)){
      cat(files[i],"\n")
      d <- read.delim(files[i],header=F,stringsAsFactors = F)
      if(i==1){
        cnt <- d
      }else{
        cnt <- cbind(cnt,d$V3)
      }
      rm(d,i)
    }
    names(cnt)[3:ncol(cnt)] <- gsub(paste0(".",f,"counts.txt"),"",files)
    cntlist[[f]] <- cnt
    cat("\n")
  }
  save(cntlist,file=paste0(DATADIR,"MiscCounts_WT1.RData"))
}
setwd(WD)

## Note: bowtie2, -q 1
## Note: bowtie2, -q 30

###############################################################################################
### Differential gene expression
###############################################################################################
if(!file.exists(paste0(DATADIR,"DEExpression_DESeq2.RData"))){
  load(paste0(DATADIR,"MiscCounts.RData"))
  load(paste0(DATADIR,"genes.RData"))
  m <- cntlist[["bam"]]
  m <- m[-c((nrow(m)-4):nrow(m)),]
  rownames(m) <- m$V1
  names(m) <- gsub("[.]bam","",names(m))
  m$V1 <- NULL
  muq <- m[,grep("[.]q30",names(m))]
  names(muq) <- gsub("[.]q30","",names(muq))
  m <- m[,-grep("[.]q30",names(m))]
  
  genes <- unique(gGenes[,c(1:6,11)])
  muq <- subset(muq, rownames(muq)%in%genes$gene_id)
  
  type <- names(sampl)
  condition <- unname(sampl[type])
  condition <- gsub("_R\\S*","",condition)
  
  ## D7 v/s control
  out0 <- myDeSeq2(df2 = muq[,type[c(1:3,7:9)]],type = type[c(1:3,7:9)],condition = condition[c(1:3,7:9)])
  ## D14 vs control
  out1 <- myDeSeq2(df2 = muq[,type[c(4:9)]],type = type[c(4:9)],condition = condition[c(4:9)])
  
  ## D14 vs control
  outC <- myDeSeq2(df2 = muq[,type],type = type,condition = condition)
  
  out0[[1]] <- merge(out0[[1]],genes,by.x="id",by.y="gene_id")
  out1[[1]] <- merge(out1[[1]],genes,by.x="id",by.y="gene_id")
  outC[[1]] <- merge(outC[[1]],genes,by.x="id",by.y="gene_id")
  
  save(out0,out1,outC,file=paste0(DATADIR,"DEExpression_DESeq2.RData"))
  
}


###############################################################################################
### Differential gene expression _WT1 data
###############################################################################################
if(!file.exists(paste0(DATADIR,"DEExpression_DESeq2_WT1.RData"))){
  load(paste0(DATADIR,"MiscCounts_WT1.RData"))
  load(paste0(DATADIR,"genes.RData"))
  
  n <- cntlist[["gene"]]
  n <- n[-c((nrow(n)-4):nrow(n)),]
  rownames(n) <- n$V1
  n$V1 <- NULL
  nuq <- n[,grep("[.]q30",names(n))]
  names(nuq) <- gsub("[.]q30","",names(nuq))
  n <- n[,-grep("[.]q30",names(n))]
  n$V2 <- NULL
  
  genes <- unique(gGenes[,c(1:6,11)])
  nuq <- subset(nuq, rownames(nuq)%in%genes$gene_id)
  
  sampl <- c("WT1_M-48_S14"="Control_R1",  "WT1_M-85_S17"="Control_R2",  "WT1_M-93_S18"="Control_R3",
             "WT1_M-37_S10"="D7_R1",  "WT1_M-47_S13"="D7_R2",  "WT1_M-76_S15"="D7_R3",
             "WT1_M-39_S11"="D14_R1",  "WT1_M-45_S12"="D14_R2",  "WT1_M-82_S16"="D14_R3")
  names(nuq) <- names(n) <-   unname(sampl[names(nuq)])
  type <- names(nuq)
  condition <- gsub("_R\\S*","",type)
  
  ## D7 v/s control
  out2 <- myDeSeq2(df2 = nuq[,type[c(1,4,6,5,8,9)]],type = type[c(1,4,6,5,8,9)],condition = condition[c(1,4,6,5,8,9)])
  ## D14 vs control
  out3 <- myDeSeq2(df2 = nuq[,type[c(2,3,7,5,8,9)]],type = type[c(2,3,7,5,8,9)],condition = condition[c(2,3,7,5,8,9)])
  
  ## D14 vs control
  outCo <- myDeSeq2(df2 = nuq[,type],type = type,condition = condition)
  
  outCo[[1]] <- merge(outCo[[1]],genes,by.x="id",by.y="gene_id")
  out2[[1]] <- merge(out2[[1]],genes,by.x="id",by.y="gene_id")
  out3[[1]] <- merge(out3[[1]],genes,by.x="id",by.y="gene_id")
  
  save(out2,out3,outCo,file=paste0(DATADIR,"DEExpression_DESeq2_WT1.RData"))
  
  plotMA(out2[[2]])
}


###############################################################################################
### Differential gene expression, GO enrichment
###############################################################################################
INPUT="A:/work/Kreidberg lab/Valerie/Bam0_CuffDiff/genes.fpkm_tracking"
INPUT1="A:/work/mm9_Annotations/ucsc_genes.txt"
INPUT2="A:/work/Kreidberg lab/Valerie/Bam0_CuffDiff/genes.read_group_tracking"

if(!file.exists( paste0(DATADIR,"DEgenes_GOEnrichment_DESeq2.RData") )){
  
  load(paste0(DATADIR,"DEExpression_DESeq2.RData"))
  
  ## expressed genes
  d <- read.delim(INPUT,header=T,stringsAsFactors = F)
  d <- d[,c("gene_id","gene_short_name")] %>% unique
  z <- read.delim(INPUT1,header=T)
  z <- unique(z[,c(1,7)])
  names(z) <- c("ucsc","gene_short_name")
  d <- merge(d,z,by="gene_short_name")
  x <- read.delim(INPUT2,header=T,stringsAsFactors = F)
  x$condition <- paste0(x$condition,x$replicate)
  x <- x[,c(1,2,7)]
  x <- reshape2::dcast(x,tracking_id~condition,value.var = "FPKM")
  exp.genes <- x[rowSums(x[,2:10])>0,]$tracking_id
  exp.genes <- d[d$gene_id%in%exp.genes,]$ucsc %>% unique
  rm(x,z)
  
  ##
  d7 <- out0[[1]]
  d7 <- merge(d7,d,by.x="id",by.y="gene_id")
  d7 <- d7[!is.na(d7$padj) & d7$padj<0.05,]
  d14 <- out1[[1]]
  d14 <- merge(d14,d,by.x="id",by.y="gene_id")
  d14 <- d14[!is.na(d14$padj) & d14$padj<0.05,]
  d7$ucsc <- as.character(d7$ucsc)
  d14$ucsc <- as.character(d14$ucsc)
  exp.genes <- as.character(exp.genes)
  
  ##
  pl <- c()
  pl <- cbind(mm9_GOStats.enrichment(d7[d7[,3]>0,]$ucsc,p.val = 0.1,universe = exp.genes),comparison="D7_Down")
  pl <- rbind(pl, cbind(mm9_GOStats.enrichment(d7[d7[,3]<0,]$ucsc,p.val = 0.1,universe = exp.genes),comparison="D7_Up"))
  pl <- rbind(pl, cbind(mm9_GOStats.enrichment(d14[d14[,3]>0,]$ucsc,p.val = 0.1,universe = exp.genes),comparison="D14_Down"))
  pl <- rbind(pl, cbind(mm9_GOStats.enrichment(d14[d14[,3]<0,]$ucsc,p.val = 0.1,universe = exp.genes),comparison="D14_Up"))
  
  ##
  save(pl,file=paste0(DATADIR,"DEgenes_GOEnrichment_DESeq2.RData"))
  rm(pl,out0,out1,d,d7,d14,exp.genes)
}


###############################################################################################
### Old and New data comparison
###############################################################################################
if(F){
  library(sva)
  library(limma)
  
  load(paste0(DATADIR,"MiscCounts.RData"))
  m <- cntlist[["bam"]]
  m <- m[-c((nrow(m)-4):nrow(m)),]
  rownames(m) <- m$V1
  names(m) <- gsub("[.]bam","",names(m))
  m$V1 <- NULL
  muq <- m[,grep("[.]q30",names(m))]
  names(muq) <- gsub("[.]q30","",names(muq))
  m <- m[,-grep("[.]q30",names(m))]
  
  load(paste0(DATADIR,"MiscCounts_WT1.RData"))
  n <- cntlist[["gene"]]
  n <- n[-c((nrow(n)-4):nrow(n)),]
  rownames(n) <- n$V1
  n$V1 <- NULL
  nuq <- n[,grep("[.]q30",names(n))]
  names(nuq) <- gsub("[.]q30","",names(nuq))
  n <- n[,-grep("[.]q30",names(n))]
  n$V2 <- NULL
  
  sampl <- c("VS1247"="D7_R1","VS1249"="D7_R2","VS36"="D7_R3","VS838"="D14_R1",
             "VS839"="D14_R2","VS840"="D14_R3","VS868"="Control_R1",
             "VS679"="Control_R2","VS582"="Control_R3",
             "WT1_M-48_S14"="Control_R1",  "WT1_M-85_S17"="Control_R2",  "WT1_M-93_S18"="Control_R3",
             "WT1_M-37_S10"="D7_R1",  "WT1_M-47_S13"="D7_R2",  "WT1_M-76_S15"="D7_R3",
             "WT1_M-39_S11"="D14_R1",  "WT1_M-45_S12"="D14_R2",  "WT1_M-82_S16"="D14_R3")
  
  names(nuq) <-   unname(sampl[names(nuq)])
  names(muq) <-   unname(sampl[names(muq)])
  names(nuq) <- paste0(names(nuq),"B")
  nuq <- nuq[match(rownames(muq),rownames(nuq)),]
  rownames(nuq)==rownames(muq)
  z <- cbind(muq,nuq) 
  rm(m,n,cntlist)
  
  ### DESeq2
  type <- names(z)
  condition <- gsub("_\\S+","",type)
  batch=c(rep("new",9),rep("old",9))
  out <- myDeSeq2(df2 = z,type = type,condition = condition,batch)
  save(out,file = paste0(DATADIR,"WT1_ValerieCombibedDEExpn_DESeq2_batch.RData"))
  
  ### SVA
  z1 <- z[rowSums(z)>10,]
  info <- data.frame(sampleID = gsub("B$","",type), sample=condition, Batch=c(rep("New",9), rep("Old",9)))
  batch = info$Batch
  modcombat = model.matrix(~1, data=info)
  cf_sva = round(ComBat(dat=as.matrix(z1), batch=as.character(batch), mod=modcombat, par.prior=TRUE, prior.plots=TRUE))
  cf_sva <- as.data.frame(cf_sva)
  save(cf_sva,file = paste0(DATADIR,"WT1_ValerieCombibedDEExpn_DESeq2_SVA.RData"))
  
}


###############################################################################################
### miRNA gene
###############################################################################################
 ## new
if(F){
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
  cnt[,1:9] <- log2(cnt[,1:9]+1)
  cnt$txID <- rownames(cnt)
  g <- merge(cnt,unique(genes[,c("chr", "strand", "geneStart", "geneEnd","gene")]),by="gene")
  g <- g[,c( "chr", "strand", "geneStart", "geneEnd","gene","txID", 
             "D7_R1", "D7_R2", "D7_R3", "D14_R1", "D14_R2", "D14_R3", 
            "control_R1", "control_R2", "control_R3")]         
  write.table(g,file=paste0(DATADIR,"GenewiseFPKM_Log2_DESeq2.txt"),sep = "\t",quote = F,row.names = F)
  write.table(g[g$gene=="Mir193a",],file=paste0(DATADIR,"GenewiseFPKM_Log2_DESeq2_miR193a.txt"),sep = "\t",quote = F,row.names = F)
  
}
 ## WT1
if(F){
  load(paste0(DATADIR,"DEExpression_DESeq2_WT1.RData"))
  load(paste0(DATADIR,"genes.RData"))
  genes <- unique(gGenes[,1:6])
  genes$width <- abs(genes$geneStart-genes$geneEnd)
  cnt <- DESeq2::counts(outCo[[2]],normalized=T) %>% as.data.frame
  cnt$id <- rownames(cnt)
  cnt <- merge(cnt,unique(genes[,c("gene","gene_id","width")]),by.x="id",by.y="gene_id")
  cnt <- unique(cnt)
  rownames(cnt) <- cnt$id
  cnt$id <- NULL
  cnt$width <-cnt$width/1000
  cnt[,1:9] <- round(cnt[,1:9]/cnt$width,2)
  cnt$width <- NULL
  cnt[,1:9] <- log2(cnt[,1:9]+1)
  cnt$txID <- rownames(cnt)
  g <- merge(cnt,unique(genes[,c("chr", "strand", "geneStart", "geneEnd","gene")]),by="gene")
  g <- g[,c( "chr", "strand", "geneStart", "geneEnd","gene","txID", 
             "D7_R1", "D7_R2", "D7_R3", "D14_R1", "D14_R2", "D14_R3", 
             "Control_R1", "Control_R2", "Control_R3")]         
  names(g) <- gsub("D7_","D9_",names(g))
  write.table(g,file=paste0(DATADIR,"GenewiseFPKM_Log2_DESeq2_WT1.txt"),sep = "\t",quote = F,row.names = F)
  write.table(g[g$gene=="Mir193a",],file=paste0(DATADIR,"GenewiseFPKM_Log2_DESeq2_miR193a_WT1.txt"),sep = "\t",quote = F,row.names = F)
  
}             
             