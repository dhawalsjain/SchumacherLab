rm(list=ls())
source("functions.R")
source("vars.R")
MIRFAMILY="A:/work/mm10_Annotations/miR_Family_Info.txt"

###############################################################################################
### generate edited mRNA sequences for miR binding assessment
###############################################################################################
 
 ## transcript sequences (mRNA)
if(F){
  library(EnsDb.Mmusculus.v79)
  library(AnnotationHub)
  edb <- EnsDb.Mmusculus.v79
  dna <- ensembldb:::getGenomeTwoBitFile(edb)
  load(paste0(DATADIR,"33_EditingLevelReport_Anno.RData"))
  outr <- hc_edits_anno[grep("3_prime|5_prime",hc_edits_anno$snpeff),] 
  outr <- as(outr,"GRanges")
  outr <- outr[,1:3]
  outr <- unique(outr[,1:3])
  o1 <- genomeToTranscript(outr,edb) 
  gr <- c()
  for(i in 1:length(o1)){
    gr <- rbind(gr,cbind(values(o1[[i]]),start=start(o1[[i]]), end=start(o1[[i]]),
                         idx=outr[i]$idx, TxID=names(o1[[i]])))
  }
  Tx <- extractTranscriptSeqs(dna,exonsBy(edb, "tx", filter = TxIdFilter( gr$TxID ) ) )
  #rm(outr,o1,edb,dna,hc_edits_anno,hlpr,imputed_edits)
  
  gr <- as.data.frame(gr)
  gr <- subset(gr,gr$TxID!="")
  Tx <- Tx[match(gr$TxID,names(Tx))]
  g <- with(gr,GRanges(paste0(TxID,"_",start),IRanges(start,end),"*",idx))
  names(Tx) <- paste0(gr$TxID,"_",gr$start)
  g$REF <- as.data.frame(Tx[g])[,1]
  g$ALT <- NA
  g$ALT <- ifelse(g$REF == "A", as.character("G"),as.character(g$ALT))
  g$ALT <- ifelse(g$REF == "T", as.character("C"),as.character(g$ALT))
  seqlengths(g) <- width(Tx)
  
  altTX <- DNAStringSet()
  for(i in 1:length(Tx)){
    if(i%%100 ==0){cat(i,"\n")}
    gx <- subset(g,seqnames(g)==names(Tx)[i])
    altTX <- c(altTX,  replaceAt(Tx[i],at = IRanges(start(gx),end(gx)),gx$ALT))
    rm(gx)
  }
  width(Tx)==width(altTX)
  
  mirs <- read.delim(MIRFAMILY,header=T,stringsAsFactors = F)
  mirs <- mirs[grep("mmu",mirs$MiRBase.ID),]
  mirs$Mature.sequence <- gsub("U","T",mirs$Mature.sequence)
  #mirs <- mirs[mirs$MiRBase.Accession!="",]
  mseq <- DNAStringSet(as.character(mirs$Mature.sequence),start=1) 
  names(mseq) <- mirs$MiRBase.ID
  
  g$editloc <- start(g)
  g$altstart <- start(g)
  g$altcorstart <- g$editloc - g$altstart +1
  save(Tx,altTX,g,gr,file = paste0(DATADIR,"33_miRpreds_intermediateTmp.RData"))
  
}

### Seed-sequence matching for miRNAs 
if(F){
  load(paste0(DATADIR,"33_miRpreds_intermediateTmp.RData"))
  width(Tx)==width(altTX)
  nn <- names(Tx)
  Tx <- Tx[trim(resize(g,13,"center"))]
  altTX <- altTX[trim(resize(g,13,"center"))]
  names(Tx) <- names(altTX) <- nn
  
  mirs <- read.delim(MIRFAMILY,header=T,stringsAsFactors = F)
  mirs <- mirs[grep("mmu",mirs$MiRBase.ID),]
  mirs$Seed.m8 <- gsub("U","T",mirs$Seed.m8)
  mirs <- unique(mirs[,1:2])
  #mirs <- mirs[mirs$MiRBase.Accession!="",]
  mseq <- DNAStringSet(as.character(mirs$Seed.m8),start=1) 
  names(mseq) <- mirs$miR.family
  mseq <- reverseComplement(mseq)
  
  #grep(mseq[1],altTX)
  
  mirSeedOlp <- c()
  for(i in 1:length(mseq)){
    cat(i,": ", names(mseq)[i],"\n")
    k <- vmatchPattern(DNAString(as.character(mseq[i])),Tx,max.mismatch = 0) %>% unlist
    k <- g[seqnames(g)%in%names(k)] %>% as.data.frame
    if(nrow(k)>0){
      k$mir <- names(mseq)[i]
      k$seq <- "ref"
      mirSeedOlp<- rbind(mirSeedOlp,k)
    }
    rm(k)
    
    k <- (vmatchPattern(DNAString(as.character(mseq[i])),altTX,max.mismatch = 0)) %>% unlist
    k <- g[seqnames(g)%in%names(k)] %>% as.data.frame
    if(nrow(k)>0){
      k$mir <- names(mseq)[i]
      k$seq <- "alt"
      mirSeedOlp<- rbind(mirSeedOlp,k)
    }
    rm(k)
  }
  mirSeedOlp <- mirSeedOlp[,c("seqnames","mir","seq","idx")] %>% unique
  names(mirSeedOlp)[1] <- "names"
  save(mirSeedOlp,file=paste0(DATADIR,"33_miRSeedOverlaps_filtermiRNAs.RData"))

  
}

## genomic sequences for motif analysis (pre-mRNA)
## Used for Adar motif prediction analysis 
if(F){
  library(BSgenome.Mmusculus.UCSC.mm10)
  genome <- BSgenome.Mmusculus.UCSC.mm10
  
  load(paste0(DATADIR,"33_EditingLevelReport_Anno.RData"))
  load(paste0(DATADIR,"33_DRE_DE_matrix.RData"))
  res$seqnames <- paste0("chr",res$seqnames)
  res <- as(res,"GRanges")
  
  outr <- hc_edits_anno 
  outr <- outr[outr$seqnames%in%c(1:19,"X","Y"),]
  outr$seqnames <- paste0("chr",outr$seqnames)
  outr <- as(outr,"GRanges")
  outr <- outr[,c(1:3,14)]
  outr <- unique(outr)
  outr$d7_control <- countOverlaps(outr,res[res$adj.P.Val<0.05 & res$comparison=="D7_control"],ignore.strand=T)
  outr$d14_control <- countOverlaps(outr,res[res$adj.P.Val<0.05 & res$comparison=="D14_control"],ignore.strand=T)
  
  for(f in levels(as.factor(outr$snpeff_uq))){
    cat(f,"\n")
    o <- subset(outr, outr$snpeff_uq==f)
    seq <- getSeq(genome,resize(o,25,"center"))
    names(seq) <- o$id
    writeXStringSet(seq,file=paste0(DATADIR,"33_Edit_25bpContext_",f,".fa"))
    rm(o,seq,f)
  }
  
  
  hlpr$id <- gsub(":T>C|:A>G","",hlpr$idx)
  seq <- getSeq(genome,resize(outr,25,"center"))
  names(seq) <- outr$id
  writeXStringSet(seq,file=paste0(DATADIR,"33_Edit_25bpContext_all.fa"))
  
  h <-as.character(hlpr[hlpr$control>0,]$id) 
  h <- subset(h,h%in%names(seq))
  writeXStringSet(seq[h],file=paste0(DATADIR,"33_Edit_25bpContext_Control.fa"))
  
  h <-as.character(hlpr[hlpr$D7>0,]$id) 
  h <- subset(h,h%in%names(seq))
  writeXStringSet(seq[h],file=paste0(DATADIR,"33_Edit_25bpContext_D7.fa"))
  
  h <-as.character(hlpr[hlpr$D14>0,]$id) 
  h <- subset(h,h%in%names(seq))
  writeXStringSet(seq[h],file=paste0(DATADIR,"33_Edit_25bpContext_D14.fa"))
  
}

### Hyperedited sites -transcript sequences (mRNA) ## don't know what is this
if(F){
  
  library(EnsDb.Mmusculus.v79)
  library(AnnotationHub)
  edb <- EnsDb.Mmusculus.v79
  dna <- ensembldb:::getGenomeTwoBitFile(edb)
  load(paste0(DATADIR,"33_EditingLevelReport_Anno.RData"))
  load(paste0(DATADIR,"33_DRE_DE_matrix.RData"))
  
  outr <- hc_edits_anno[grep("3_prime|5_prime",hc_edits_anno$snpeff),] 
  outr <- outr[outr$idx%in%res[res$adj.P.Val<0.01,]$idx,]
  outr <- as(outr,"GRanges")
  outr <- outr[,1:3]
  outr <- unique(outr[,1:3])
  o1 <- genomeToTranscript(outr,edb) 
  gr <- c()
  for(i in 1:length(o1)){
    gr <- rbind(gr,cbind(values(o1[[i]]),start=start(o1[[i]]), end=start(o1[[i]]),
                         idx=outr[i]$idx, TxID=names(o1[[i]])))
  }
  gr <- gr[gr$TxID%in%expn$test_id,]
  Tx <- extractTranscriptSeqs(dna,exonsBy(edb, "tx", filter = TxIdFilter( gr$TxID ) ) )
  
  gr <- as.data.frame(gr)
  gr <- subset(gr,gr$TxID!="")
  g <- with(gr,GRanges(TxID,IRanges(start,end),"*",idx))
  g$REF <- as.data.frame(Tx[g])[,1]
  g$ALT <- NA
  g$ALT <- ifelse(g$REF == "A", as.character("G"),as.character(g$ALT))
  g$ALT <- ifelse(g$REF == "T", as.character("C"),as.character(g$ALT))
  g$id <- gsub(":A>G|:T>C","",g$idx)
  
  altTX <- DNAStringSet()
  for(i in 1:length(Tx)){
    if(i%%100 ==0){cat(i,"\n")}
    gx <- subset(g,seqnames(g)==names(Tx)[i])
    altTX <- c(altTX,  replaceAt(Tx[i],at = IRanges(start(gx),end(gx)),gx$ALT))
    rm(gx)
  }
  save(altTX,Tx,g,gr,paste0(DATADIR,"33_miRpreds_intermediateTmp_Clustered.RData"))
  
}


## Hyperedited/ reference sites -transcript sequences (pre-mRNA)
##  Used for 2ry structure prediction
if(F){
  library(BSgenome.Mmusculus.UCSC.mm10)
  genome <- BSgenome.Mmusculus.UCSC.mm10
  
  load(paste0(DATADIR,"33_EditingLevelReport_Anno.RData"))
  load(paste0(DATADIR,"33_DRE_DE_matrix.RData"))
  load(paste0(DATADIR,"genes.RData"))
  genes <- gGenes[,c("gene_id","strand")] %>% unique
  names(genes) <- c("tracking_id","gstrand")
  rm(gGenes)
  
  res$seqnames <- paste0("chr",res$seqnames)
  res$REF <- gsub(">\\S*","",gsub("\\S*:","",res$idx))
  res$ALT <- gsub("\\S*>","",res$idx)
  res <- merge(res,genes,by="tracking_id")
  res <- as(res,"GRanges")
  
  Cntr <- DNAStringSet()
  Exp <- DNAStringSet()
  nameS <- c()
  
  r1 <- res[res$comparison=="D7_control" & res$clustered>0 & !is.na(res$tracking_id),]
  r1 <- r1[r1$adj.P.Val<0.05,]
  r1 <- split(r1,r1$tracking_id)
  for(i in 1:length(r1)){
    if(min(r1[[i]]$adj.P.Val)<0.05 & length(r1[[i]])>1 ){
      cat(i,"\n")
      g <- GRanges(seqnames(r1[[i]])[1],IRanges(min(start(r1[[i]])), max(start(r1[[i]]))))
      g <- reduce(c(flank(g,100),g,flank(g,100,start = F)))
      gx <- r1[[i]]
      gx$REF <- as.character(getSeq(genome,gx))
      gx$ALT <- NA
      gx$ALT <- ifelse(gx$REF == "A", as.character("G"),as.character(gx$ALT))
      gx$ALT <- ifelse(gx$REF == "T", as.character("C"),as.character(gx$ALT))
      start(gx) <- (start(gx) - min(start(r1[[i]]))+1+100)
      end(gx) <- start(gx)
      sq <- getSeq(genome,g)
      sqC <- sqE <- sq
      gx1 <- gx[gx$control>gx$D7]
      gx2 <- gx[gx$control<gx$D7]
      
      if(length(gx1)>0){
        sqC <- replaceAt(sq,at = IRanges(start(gx1),end(gx1)),value = as.character(gx1$ALT) )
      }
      if(length(gx2)>0){
        sqE <- replaceAt(sq,at = IRanges(start(gx2),end(gx2)),value = as.character(gx2$ALT) )
      }
      
      if( unique(gx$gstrand) == "+" ){
        Cntr <- c(Cntr,complement(sqC) )
        Exp <- c(Exp,complement(sqE) )
      }else{
        Cntr <- c(Cntr,sqC)
        Exp <- c(Exp,sqE)
      }
      x <- paste0(unique(r1[[i]]$gene), "|",length(r1[[i]]),"|",min(r1[[i]]$adj.P.Val),"|",paste0(unique(r1[[i]]$snpeff_uq),collapse = "-"),"|",round(median(r1[[i]]$logFC),2),"|D7_control")
      nameS <- c(nameS, x)
      rm(g,gx,gx1,gx2,sq,sqC,sqE,x)    
    }
  }
  
  r1 <- res[res$comparison=="D14_control" & res$clustered>0 & !is.na(res$tracking_id),]
  r1 <- r1[r1$adj.P.Val<0.05,]
  r1 <- split(r1,r1$tracking_id)
  for(i in 1:length(r1)){
    if(min(r1[[i]]$adj.P.Val)<0.05 & length(r1[[i]])>1 ){
      cat(i,"\n")
      g <- GRanges(seqnames(r1[[i]])[1],IRanges(min(start(r1[[i]])), max(start(r1[[i]]))))
      g <- reduce(c(flank(g,100),g,flank(g,100,start = F)))
      gx <- r1[[i]]
      gx$REF <- as.character(getSeq(genome,gx))
      gx$ALT <- NA
      gx$ALT <- ifelse(gx$REF == "A", as.character("G"),as.character(gx$ALT))
      gx$ALT <- ifelse(gx$REF == "T", as.character("C"),as.character(gx$ALT))
      start(gx) <- (start(gx) - min(start(r1[[i]]))+1+100)
      end(gx) <- start(gx)
      sq <- getSeq(genome,g)
      sqC <- sqE <- sq
      gx1 <- gx[gx$control>gx$D14]
      gx2 <- gx[gx$control<gx$D14]
      
      if(length(gx1)>0){
        sqC <- replaceAt(sq,at = IRanges(start(gx1),end(gx1)),value = as.character(gx1$ALT) )
      }
      if(length(gx2)>0){
        sqE <- replaceAt(sq,at = IRanges(start(gx2),end(gx2)),value = as.character(gx2$ALT) )
      }
      
      if( unique(gx$gstrand) == "+" ){
        Cntr <- c(Cntr,complement(sqC) )
        Exp <- c(Exp,complement(sqE) )
      }else{
        Cntr <- c(Cntr,sqC)
        Exp <- c(Exp,sqE)
      }
      x <- paste0(unique(r1[[i]]$gene), "|",length(r1[[i]]),"|",min(r1[[i]]$adj.P.Val),"|",paste0(unique(r1[[i]]$snpeff_uq),collapse = "-"),"|",round(median(r1[[i]]$logFC),2),"|D14_control")
      nameS <- c(nameS, x)
      rm(g,gx,gx1,gx2,sq,sqC,sqE,x)    
    }
  }
  
  names(Cntr) <- names(Exp) <- nameS
  Cntr <- Cntr[width(Cntr)<1000]
  Exp <- Exp[width(Exp)<1000]
  
  save(Cntr,Exp,file=paste0(DATADIR,"33_EditingClusters_GenomicSequences.RData"))
  writeXStringSet(Cntr,file=paste0(DATADIR,"33_DRE_Clustered_RefSeqs_Genomic.fa"))
  writeXStringSet(Exp,file=paste0(DATADIR,"33_DRE_Clustered_AltSeqs_Genomic.fa"))
  
  
}


