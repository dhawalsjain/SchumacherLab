rm(list = ls())
source("functions.R")
source("vars.R")
UTR3="A:/work/mm10_Annotations/mm10_3UTR"
UTR5="A:/work/mm10_Annotations/mm10_5UTR"
INTRON="A:/work/mm10_Annotations/mm10_intorns"
EXON="A:/work/mm10_Annotations/mm10_exons"
UCE="A:/work/mm10_Annotations/uce.txt.txt"

###############################################################################################
### Editing level and data imputation and primary annotations
###############################################################################################
if(F){
  load(paste0(DATADIR,"CombinedEditsPool_q10.RData"))
  
  q <- unlist(l1) %>% as.data.frame
  q$profile <- rownames(q)
  names(q) <- c("sample_name","profile")
  libs <- merge(libs,q,by="profile")
  rm(q)
  
  out <- merge(out,libs[,c("profile","cov","sample_name")],all.x=T)
  
  ## High confidence edit calls
  o <- out[is.na(out$Insertion),]
  o <- o[is.na(o$Deletion),]
  o$H <- o$alt/o$depth
  #o <- o[o$H<1,]
  o <- o[o$H>0,]
  o <- o[o$V4%in%c("A>G","T>C"),]
  o <- o[o$depth>=10,]  ## need minimum 10 reads at the locus, no matter what the coverage is 
  
  ## data imputation
  z <- reshape2::dcast(o,id~profile,value.var = "H")
  rownames(z) <- z$id
  z$id <- NULL
  names(z) <- unname(unlist(l1[names(z)]))
  
  myfun <- function(x){
    if(sum(is.na(x))>0){
      return(rep(0,length(x)))
    }else{
      return(x)
    }
  }
  z1 <- apply(z[,grep("D7",names(z))],1,myfun) %>% t %>% as.data.frame
  z2 <- apply(z[,grep("D14",names(z))],1,myfun) %>% t %>% as.data.frame
  z3 <- apply(z[,grep("control",names(z))],1,myfun) %>% t %>% as.data.frame
  names(z1) <- names(z[,grep("D7",names(z))])
  names(z2) <- names(z[,grep("D14",names(z))])
  names(z3) <- names(z[,grep("control",names(z))])
  z <- cbind(z1,z2,z3)
  rm(z1,z2,z3)
  z <- round(z,3)
  z <- z[rowSums(z)>0,]
  z <- z[rowSums(z)<9,]
  
  ## known vs unknown
  o <- o[o$id%in%rownames(z),]
  db <- read.delim(paste0(DATADIR,"CombinedEditDB.bed"),header = F,stringsAsFactors = F)
  db$id <- paste0(db$V1,":",db$V2)
  o$known <- ifelse(o$id%in%db$id,1,0)
  
  ## repeat annotations
  load(paste0(DATADIR,"repeats.RData"))
  o <- with(o,GRanges(chr,IRanges(loc,loc),"*",id,V4,Insertion,Deletion,cov,sample_name,profile,H,known))
  o$repeats <- "Nonrepetitive"
  olap <- findOverlaps(o,alu,ignore.strand=T)
  o[queryHits(olap)]$repeats <- "Alu"
  olap <- findOverlaps(o,mm10,ignore.strand=T)
  o[queryHits(olap)]$repeats <- "Repetitive non-Alu"
  o <- as.data.frame(o)
  
  ## genomic annotations
  o$idx <- paste0(o$id,":",o$V4)
  regions <- read.delim(paste0(DATADIR,"CombinedPool.bed"),header = F,stringsAsFactors = F)
  regions$idx <- paste0(regions$V1,":",regions$V2,":",regions$V4)
  names(regions)[7:8] <- c("snpeff","aaswap")
  length(unique(regions$idx))==length(regions$idx)
  o <- merge(o,unique(regions[,c("idx","snpeff","aaswap")]),by="idx")
  
  hc_edits <- o
  imputed_edits <- z
  save(hc_edits,imputed_edits,file=paste0(DATADIR,"33_EditingLevelReport.RData"))
  
  rm(o,z,hc_edits,imputed_edits,db,alu,mm10,regions)
}

###############################################################################################
### Annotations for clusters, uniqueness etc 
###############################################################################################

if(F){
  myf <- function(path="A:/work/mm10_Annotations/mm10_3UTR"){
    z <- read.delim(path,header = F)
    z$V1 <- gsub("chr","",z$V1)
    z <- z[,c(1:3,6)]
    z <-with(z,GRanges(V1,IRanges(V2,V3),V6)) 
    z <- reduce(z,ignore.strand=T)
    z <- subset(z,seqnames(z)%in%c("X","Y",paste0("",1:22)))
    z
  }
  utr3 <- myf(UTR3)
  utr5 <- myf(UTR5)
  introns <- myf(INTRON)  
  exons <- myf(EXON)
  
  load(paste0(DATADIR,"33_EditingLevelReport.RData"))
  uqedits <- hc_edits[,c("seqnames","start","strand","V4")] %>% unique
  uqedits <- with(uqedits,GRanges(seqnames,IRanges(start,start),strand,edit=V4))
  uqedits$id <- paste0(seqnames(uqedits),":",start(uqedits))
  uqedits$idx <- paste0(uqedits$id,":",uqedits$edit)
  imp <- imputed_edits
  imp[imp==0]<- NA
  uqedits$intron <- countOverlaps(uqedits,introns,ignore.strand=T)
  uqedits$exon <- countOverlaps(uqedits,exons,ignore.strand=T)
  uqedits$utr3 <- countOverlaps(uqedits,utr3,ignore.strand=T)
  uqedits$utr5 <- countOverlaps(uqedits,utr5,ignore.strand=T)
  uqedits <- as.data.frame(uqedits)
  uqedits$intergenic <- apply(uqedits[,9:12],1,sum)
  uqedits$intergenic <- ifelse(uqedits$intergenic>0,100,200)
  uqedits$intergenic <- ifelse(uqedits$intergenic==200,1,0)
  uqedits <- as(uqedits,"GRanges")
  editC <- uqedits[uqedits$id%in%rownames(imp[!is.na(imp$control_R3),])] %>% as.data.frame
  editD7 <- uqedits[uqedits$id%in%rownames(imp[!is.na(imp$D7_R1),])] %>% as.data.frame
  editD14 <- uqedits[uqedits$id%in%rownames(imp[!is.na(imp$D14_R1),])] %>% as.data.frame
  rm(imp,uqedits)
  
  hlpr <- hc_edits[,c("idx","known","repeats","snpeff","aaswap")] %>% unique
  hlpr$control <- ifelse(hlpr$idx%in%editC$idx,1,0)
  hlpr$D7 <- ifelse(hlpr$idx%in%editD7$idx,1,0)
  hlpr$D14 <- ifelse(hlpr$idx%in%editD14$idx,1,0)
  hlpr$known <- ifelse(hlpr$known==0,as.character("novel"),as.character("known"))
  
  uce <- read.delim(UCE,header=T,stringsAsFactors = F)
  head(uce)
  uce$chr <- gsub("chr","",uce$chr)
  uce <- with(uce,GRanges(chr,IRanges(seq_start,seq_stop),"*",UCNE_name,UCNE_ID))
  qw <-  as(rbind(editC,editD14,editD7),"GRanges")
  qw <- unique(qw)
  olp <- findOverlaps(qw,resize(uce,width(uce)+101,"center"),ignore.strand=T) %>% as.data.frame
  olp$idx <- qw[olp$queryHits]$idx
  olp$UCNE_name <- uce[olp$subjectHits]$UCNE_name
  olp$UCNE_ID <- uce[olp$subjectHits]$UCNE_ID
  olp$queryHits <- olp$subjectHits <- NULL
  hlpr <- merge(hlpr,olp,by='idx',all.x=T)
  rm(uce,qw,olp)
  
  hc_edits_anno=hc_edits
  
  categories <- c("missense","synonymous","splice","3_prime_UTR", "5_prime_UTR","intron","intergenic_region")
  mm <- c()
  for(f in categories){
    mm <- rbind(mm,data.frame(idx=hc_edits_anno[grep(f,hc_edits_anno$snpeff),]$idx,anno=f,chk=1))
  }
  mm <- reshape2::dcast(mm,idx~anno,value.var = "chk")
  mm$anno <- NA
  for(f in categories){
    mm$anno <- ifelse(mm[,paste0(f)]>0 & is.na(mm$anno),as.character(f),mm$anno)
  }
  mm$snpeff_uq <- mm$anno
  mm <- mm[,c("idx","snpeff_uq")]
  idx <- unique(hc_edits_anno$idx)
  mm <- rbind(mm,data.frame(idx=idx[!idx%in%mm$idx],snpeff_uq=NA))
  
  hc_edits_anno <- merge(hc_edits_anno,mm[,c("idx","snpeff_uq")],by="idx")
  hlpr <- merge(hlpr,mm,by="idx")
  
  save(hc_edits_anno,imputed_edits,hlpr,file=paste0(DATADIR,"33_EditingLevelReport_Anno.RData"))
  rm(hc_edits,hc_edits_anno,hlpr,imputed_edits,editC,editD7,editD14)
  
}


## conservation according to 60Way MSA of vertebrates
if(F){
  load(paste0(DATADIR,"33_EditingLevelReport_Anno.RData"))
  m <- hc_edits_anno[,c(2:4,1)] %>% unique
  m$seqnames <- paste0("chr",m$seqnames)
  m$end <- m$start
  m$start <- m$start-1
  
  write.table(m[1:1000,],file=paste0(RAWDIR,"33_HC_EditingEvents1.bed"),quote = F,sep = "\t",row.names = F,col.names = F)
  write.table(m[1001:2000,],file=paste0(RAWDIR,"33_HC_EditingEvents2.bed"),quote = F,sep = "\t",row.names = F,col.names = F)
  write.table(m[2001:nrow(m),],file=paste0(RAWDIR,"33_HC_EditingEvents3.bed"),quote = F,sep = "\t",row.names = F,col.names = F)
  
  m$x <- paste0(m$seqnames,":",m$start,"-",m$end)
  write.table(m$x,file=paste0(RAWDIR,"33_HC_EditingEvents.liftOver"),quote = F,sep = "\t",row.names = F,col.names = F)
  
  
  m <- read.delim(paste0(RAWDIR,"33_HC_EditingEvents_PhastCons60Way.txt"),header=F,comment.char = "#",stringsAsFactors = F)
  m$V1 <- gsub("chr","",m$V1)
  m$id <- paste0(m$V1,":",m$V2)
  sum(m$id%in%hc_edits_anno$id)
  names(m)[10] <- "ConsScore"
  m <- m[,c("id","ConsScore")]
  m <- unique(m)
  length(unique(m$id))
  hc_edits_anno <- merge(hc_edits_anno,m,by="id",all.x=T)
  
  save(hc_edits_anno,imputed_edits,hlpr,file=paste0(DATADIR,"33_EditingLevelReport_Anno.RData"))
  
}


## conservation in human?
if(F){
  load(paste0(DATADIR,"33_EditingLevelReport_Anno.RData"))
  load(paste0(DATADIR,"hg19_RNAEdits_DB.RData"))
  
  a <- read.delim(paste0(RAWDIR,"33_HC_EditingEvents.liftover"),header=F,stringsAsFactors = F)
  f <- read.delim(paste0(RAWDIR,"33_liftoverFailure_mm102hg19.txt"),header=F,stringsAsFactors = F,comment.char = "#")
  b <- read.delim(paste0(RAWDIR,"33_liftover_mm102hg19.txt.bed"),header=F,stringsAsFactors = F)
  b$chr <- gsub("chr","",gsub(":\\S*","",b$V1))
  b$end <- gsub("\\S*-","",b$V1) %>% as.numeric
  b$start <- b$end-1
  b <- as(b,"GRanges")
  o <- findOverlaps(resize(b,5,"center"),hg19Edits,ignore.strand=T)
  b$tissue <- b$gene <- NA
  b[queryHits(o)]$gene <- as.character(hg19Edits[subjectHits(o)]$gene)
  b[queryHits(o)]$tissue <- as.character(hg19Edits[subjectHits(o)]$tissue)
  b$human <- paste0(b$V1,"|",b$gene,"|",b$tissue)
  
  a <- a[!a$V1%in%f$V1,] %>% as.data.frame
  a$hg19 <- b$human
  a$gene <- b$gene
  a <- a[!is.na(a$gene),]
  names(a)[1] <- "V1"
  a$chr <- gsub("chr","",gsub(":\\S*","",a$V1))
  a$end <- gsub("\\S*-","",a$V1) %>% as.numeric
  a$start <- a$end-1
  a <- as(a,"GRanges")
  a$gene <- a$V1 <- NULL
  
  hc_edits_anno <- as(hc_edits_anno,"GRanges")
  o <- findOverlaps(resize(a,5,"center"),hc_edits_anno,ignore.strand=T)
  hc_edits_anno$hg19 <- NA
  hc_edits_anno[subjectHits(o)]$hg19 <- as.character(a[queryHits(o)]$hg19)
  hc_edits_anno <- as.data.frame(hc_edits_anno)
  hc_edits_anno <- unique(hc_edits_anno)
  save(hc_edits_anno,imputed_edits,hlpr,file=paste0(DATADIR,"33_EditingLevelReport_Anno.RData"))

}


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
if(F){
  files <- Sys.glob(paste0(RAWDIR,"HMM/*.gz"))
  sts <- GRanges()
  for(f in files){
    cat(f,"\n")
    d <- read.delim(f,header=T,stringsAsFactors = F)
    f <- gsub(".txt.gz","",gsub(paste0(RAWDIR,"HMM/mm10_Renal_cHMM_"),"",f))
    d$stage <- f
    d$X.chrom <- gsub("chr","",d$X.chrom)
    d <- with(d,GRanges(X.chrom,IRanges(chromStart,chromEnd),"+",state=name,stage))
    sts <- c(sts,d)
    rm(d,f)
  }
  
  load(paste0(DATADIR,"33_EditingLevelReport_Anno.RData"))
  uq <- as(hc_edits_anno,"GRanges")
  uq <- unique(uq[,"idx"])
  uq$cHMM_E14 <- uq$cHMM_E15 <- uq$cHMM_E16 <- uq$cHMM_P0 <- NA 
  
  xx <- sts[sts$stage=="E14"]
  o <- findOverlaps(uq,xx,ignore.strand=T)
  uq[queryHits(o)]$cHMM_E14 <- xx[subjectHits(o)]$state
  
  xx <- sts[sts$stage=="E15"]
  o <- findOverlaps(uq,xx,ignore.strand=T)
  uq[queryHits(o)]$cHMM_E15 <- xx[subjectHits(o)]$state
  
  xx <- sts[sts$stage=="E16"]
  o <- findOverlaps(uq,xx,ignore.strand=T)
  uq[queryHits(o)]$cHMM_E16 <- xx[subjectHits(o)]$state
  
  xx <- sts[sts$stage=="P0"]
  o <- findOverlaps(uq,xx,ignore.strand=T)
  uq[queryHits(o)]$cHMM_P0 <- xx[subjectHits(o)]$state
  
  uq <- values(uq) %>% as.data.frame
  uq$cHMM_Final <- apply(uq[,2:5],1,function(x){
     if(length(unique(x))==1){
       return(unique(x))
     }else{
       return(NA)
     }
  })
  hlpr <- merge(hlpr,uq,by="idx",all.x=T)
  save(hc_edits_anno,imputed_edits,hlpr,file=paste0(DATADIR,"33_EditingLevelReport_Anno.RData"))
  
}

## enrichment
if(F){
  files <- Sys.glob(paste0(RAWDIR,"HMM/*.gz"))
  sts <- GRanges()
  for(f in files){
    cat(f,"\n")
    d <- read.delim(f,header=T,stringsAsFactors = F)
    f <- gsub(".txt.gz","",gsub(paste0(RAWDIR,"HMM/mm10_Renal_cHMM_"),"",f))
    d$stage <- f
    d$X.chrom <- gsub("chr","",d$X.chrom)
    d <- with(d,GRanges(X.chrom,IRanges(chromStart,chromEnd),"+",state=name,stage))
    sts <- c(sts,d)
    rm(d,f)
  }
  uq <- read.delim(paste0(RAWDIR,"mm10_5kRandom.bed"),header = F,stringsAsFactors = F)
  uq <- with(uq,GRanges(V1,IRanges(V2,V2),"+"))
  uq$cHMM_E14 <- uq$cHMM_E15 <- uq$cHMM_E16 <- uq$cHMM_P0 <- NA 
  
  xx <- sts[sts$stage=="E14"]
  o <- findOverlaps(uq,xx,ignore.strand=T)
  uq[queryHits(o)]$cHMM_E14 <- xx[subjectHits(o)]$state
  
  xx <- sts[sts$stage=="E15"]
  o <- findOverlaps(uq,xx,ignore.strand=T)
  uq[queryHits(o)]$cHMM_E15 <- xx[subjectHits(o)]$state
  
  xx <- sts[sts$stage=="E16"]
  o <- findOverlaps(uq,xx,ignore.strand=T)
  uq[queryHits(o)]$cHMM_E16 <- xx[subjectHits(o)]$state
  
  xx <- sts[sts$stage=="P0"]
  o <- findOverlaps(uq,xx,ignore.strand=T)
  uq[queryHits(o)]$cHMM_P0 <- xx[subjectHits(o)]$state
  
  uq <- values(uq) %>% as.data.frame
  uq$cHMM_Final <- apply(uq[,1:4],1,function(x){
    if(length(unique(x))==1){
      return(unique(x))
    }else{
      return(NA)
    }
  })
  uq <- uq[!is.na(uq$cHMM_Final),]
  
  load(paste0(DATADIR,"33_EditingLevelReport_Anno.RData"))
  hlpr <- hlpr[!is.na(hlpr$cHMM_Final),]
  pl <- rbind( cbind(as.data.frame(table(hlpr$cHMM_Final)),descr="Exp"),
               cbind(as.data.frame(table(uq$cHMM_Final)),descr="Random")
  )
  pl <- reshape2::dcast(pl,Var1~descr,value.var = "Freq")
  pl[is.na(pl)] <- 0
  
  pl$a <- pl$Exp
  pl$b <- nrow(hlpr)-pl$a
  pl$c <- pl$Random
  pl$d <- nrow(uq) - pl$c
  
  pl$p.val <- apply(pl[,4:7],1,function(x){
    -log10(fisher.test( matrix(x,nrow=2,byrow = T),alternative = "two.sided")$p.value)
  })
  pl$p.val <- ifelse(is.infinite(pl$p.val),100,pl$p.val)
  pl$status <- "enriched"
  pl[pl$a/nrow(hlpr) < pl$c/nrow(uq),]$status <- "depleted"
  cHMMEnrich <- pl
  save(cHMMEnrich,file=paste0(DATADIR,"33_cHMM_pValEnrichment.RData"))
  
}


## add distances to the  wt1 peaks from the edit sites
if(F){
  load(paste0(DATADIR,"33_EditingLevelReport_Anno.RData"))
  df1 <- unique(hc_edits_anno[,c('seqnames','start','end','idx','repeats','snpeff_uq',"known")]) %>% as.data.frame
  df1 <- as(df1,"GRanges")
  df1$D9_Wt1 <- df1$D14_Wt1 <- df1$Control_Wt1 <- NA
  
  gr <- read.delim(paste0(RAWDIR,"wt1_mm9/D9_mm10_Wt1peaks.bed"),header=F,stringsAsFactors = F)
  gr$V1 <- gsub("chr","",gr$V1)
  gr <- with(gr,GRanges(V1,IRanges(V2,V3),"*",score=V5))
  d <- distanceToNearest(df1,gr,ignore.strand=T) %>% as.data.frame
  df1[d$queryHits,]$D9_Wt1 <- d$distance
  
  gr <- read.delim(paste0(RAWDIR,"wt1_mm9/D14_mm10_Wt1peaks.bed"),header=F,stringsAsFactors = F)
  gr$V1 <- gsub("chr","",gr$V1)
  gr <- with(gr,GRanges(V1,IRanges(V2,V3),"*",score=V5))
  d <- distanceToNearest(df1,gr,ignore.strand=T) %>% as.data.frame
  df1[d$queryHits,]$D14_Wt1 <- d$distance
  
  gr <- read.delim(paste0(RAWDIR,"wt1_mm9/WT_mm10_Wt1peaks.bed"),header=F,stringsAsFactors = F)
  gr$V1 <- gsub("chr","",gr$V1)
  gr <- with(gr,GRanges(V1,IRanges(V2,V3),"*",score=V5))
  d <- distanceToNearest(df1,gr,ignore.strand=T) %>% as.data.frame
  df1[d$queryHits,]$Control_Wt1 <- d$distance
  
  df1 <- as.data.frame(df1)
  hlpr <- merge(hlpr,df1[,c("idx","Control_Wt1","D9_Wt1","D14_Wt1")],by="idx")
  save(hc_edits_anno,imputed_edits,hlpr,file=paste0(DATADIR,"33_EditingLevelReport_Anno.RData"))
  
}

## Distance to Wt1 peaks
if(F){
  myf <- function(path="A:/work/mm10_Annotations/mm10_3UTR"){
    z <- read.delim(path,header = F)
    z$V1 <- gsub("chr","",z$V1)
    z <- z[,c(1:3,6)]
    z <-with(z,GRanges(V1,IRanges(V2,V3),V6)) 
    z <- reduce(z,ignore.strand=T)
    z <- subset(z,seqnames(z)%in%c("X","Y",paste0("",1:22)))
    z
  }
  utr3 <- myf(UTR3)
  utr5 <- myf(UTR5)
  introns <- myf(INTRON)  
  exons <- myf(EXON)
  
  load(paste0(DATADIR,"33_EditingLevelReport_Anno.RData"))
  load(paste0(DATADIR,"33_DRE_DE_matrix.RData"))
  
  
  df1 <- read.delim(paste0(RAWDIR,"mm10_5kRandom.bed"),header = F,stringsAsFactors = F)
  df1 <- with(df1,GRanges(V1,IRanges(V2,V2),"+"))
  df1$D9_Wt1 <- df1$D14_Wt1 <- df1$Control_Wt1 <- NA
  
  gr <- read.delim(paste0(RAWDIR,"wt1_mm9/D9_mm10_Wt1peaks.bed"),header=F,stringsAsFactors = F)
  gr$V1 <- gsub("chr","",gr$V1)
  gr <- with(gr,GRanges(V1,IRanges(V2,V3),"*",score=V5))
  d <- distanceToNearest(df1,gr,ignore.strand=T) %>% as.data.frame
  df1[d$queryHits,]$D9_Wt1 <- d$distance
  
  gr <- read.delim(paste0(RAWDIR,"wt1_mm9/D14_mm10_Wt1peaks.bed"),header=F,stringsAsFactors = F)
  gr$V1 <- gsub("chr","",gr$V1)
  gr <- with(gr,GRanges(V1,IRanges(V2,V3),"*",score=V5))
  d <- distanceToNearest(df1,gr,ignore.strand=T) %>% as.data.frame
  df1[d$queryHits,]$D14_Wt1 <- d$distance
  
  gr <- read.delim(paste0(RAWDIR,"wt1_mm9/WT_mm10_Wt1peaks.bed"),header=F,stringsAsFactors = F)
  gr$V1 <- gsub("chr","",gr$V1)
  gr <- with(gr,GRanges(V1,IRanges(V2,V3),"*",score=V5))
  d <- distanceToNearest(df1,gr,ignore.strand=T) %>% as.data.frame
  df1[d$queryHits,]$Control_Wt1 <- d$distance
  df1 <- as.data.frame(df1)
  head(df1)
  
  df1 <- as(df1,"GRanges")
  df1$utr3 <- countOverlaps(df1,utr3,ignore.strand=T)
  df1$utr5 <- countOverlaps(df1,utr5,ignore.strand=T)
  df1$exon <- countOverlaps(df1,exons,ignore.strand=T)
  df1$intron <- countOverlaps(df1,introns,ignore.strand=T)
  df1 <- as.data.frame(df1)
  df1$anno <- NA
  df1[df1$utr3>0,]$anno <- "3UTR"
  df1[df1$utr5>0 & is.na(df1$anno),]$anno <- "5UTR"
  df1[df1$exon>0 & is.na(df1$anno),]$anno <- "Exon"
  df1[df1$intron>0 & is.na(df1$anno),]$anno <- "Intron"
  df1[is.na(df1$anno),]$anno <- "Intergenic"
  df1$utr3 <- df1$utr5 <- df1$exon <- df1$intron <- NULL
  
  uq <- as(hc_edits_anno,"GRanges") 
  uq <- uq[,"idx"] %>% unique
  uq$utr3 <- countOverlaps(uq,utr3,ignore.strand=T)
  uq$utr5 <- countOverlaps(uq,utr5,ignore.strand=T)
  uq$exon <- countOverlaps(uq,exons,ignore.strand=T)
  uq$intron <- countOverlaps(uq,introns,ignore.strand=T)
  uq <- as.data.frame(uq)
  uq$anno <- NA
  uq[uq$utr3>0,]$anno <- "3UTR"
  uq[uq$utr5>0 & is.na(uq$anno),]$anno <- "5UTR"
  uq[uq$exon>0 & is.na(uq$anno),]$anno <- "Exon"
  uq[uq$intron>0 & is.na(uq$anno),]$anno <- "Intron"
  uq[is.na(uq$anno),]$anno <- "Intergenic"
  table(uq$anno)
  
  df1 <- melt(df1,measure.vars = names(df1)[6:8])[,6:8]
  df1$set <- "random"
  df1$variable <- gsub("_Wt1","",df1$variable)
  df1$variable <- gsub("D9","D7",df1$variable)
  table(df1$variable)
  
  rd <- rbind(data.frame(variable="Control",value=hlpr[hlpr$control>0,]$Control_Wt1,set="edits",idx=hlpr[hlpr$control>0,]$idx),
              data.frame(variable="D7",value=hlpr[hlpr$D7>0,]$D9_Wt1,set="edits",idx=hlpr[hlpr$D7>0,]$idx),
              data.frame(variable="D14",value=hlpr[hlpr$D14>0,]$D14_Wt1,set="edits",idx=hlpr[hlpr$D14>0,]$idx))
  rd <- merge(rd,unique(uq[,c("idx","anno")]),by="idx")
  rd$idx <- NULL
  rd <- rd[,match( names(df1), names(rd))]
  rd <- rbind(rd,df1)
  
  save(rd,file=paste0(DATADIR,"33_Distance2Wt1Peaks_plotDF.RData"))
  
}
## Editing levels v/s read depth
if(F){
  load(paste0(DATADIR,"CombinedEditsPool_q10.RData"))
  k <- data.frame(name=names(l),profile=unlist(l))
  out <- merge(out,k,by="profile")
  out <- out[,c("id","name","depth","alt")] %>% unique
  out <- data.table(out)
  head(out)
  out = out[, lapply(.SD, sum, na.rm=TRUE), by=list(id,name)]
  out <- as.data.frame(out)
  rm(k)
  
  load(paste0(DATADIR,"33_EditingLevelReport_Anno.RData"))
  out <- out[out$id%in%rownames(imputed_edits),]
  oDepth <- out
  out$alt <- NULL
  
  z <- imputed_edits
  z$id <- rownames(z)
  z$D7 <- round(apply(z[,1:3],1,gmean),2)
  z$D14 <- round(apply(z[,4:6],1,gmean),2)
  z$control <- round(apply(z[,7:9],1,gmean),2)
  z <- z[,10:13]
  z <- melt(z,measure.vars = c("D7", "D14", "control"))
  z$idx <- paste0(z$id,"_",z$variable)
  out$idx <- paste0(out$id,"_",out$name)
  z$variable <- z$id <- NULL
  out$name <- out$id <-  NULL
  z <- merge(z,out,by="idx")
  z$profile <- gsub("\\S*_","",z$idx)
  
  eDepth <- z
  save(oDepth,eDepth,file=paste0(DATADIR,"33_readDepth_EditLevels.RData"))

}

### coplete editing annotation spreadsheet
if(T){
  load(paste0(DATADIR,"DEExpression_DESeq2.RData"))
  ggg <- out0[[1]][,c("id","geneStart","geneEnd","gene")]
  ggg$length <- abs(ggg$geneStart-ggg$geneEnd)
  ggg <- ggg[,c("id","length","gene")] %>% unique
  cnt <- DESeq2::counts(outC[[2]],normalized=F) %>% as.data.frame
  cnt$id <- rownames(cnt)
  cnt <- merge(cnt,ggg,by="id",all.x=T)
  rownames(cnt) <- cnt$id
  cnt$id <- NULL
  cnt[,1:9] <- cnt[,1:9]*1000/cnt$length
  cnt$length <- NULL
  cnt[,1:9] <- apply(cnt[,1:9],2,function(x){round(x*1e6/sum(x),2)})
  cnt <- as.data.frame(cnt)
  names(cnt)[1:9] <- unname(sampl[names(cnt)[1:9]])
  head(cnt)
  cnt$D7 <- round(apply(cnt[,1:3],1,gmean),2)
  cnt$D14 <- round(apply(cnt[,4:6],1,gmean),2)
  cnt$Control <- round(apply(cnt[,7:9],1,gmean),2)
  cnt <- cnt[,c("gene","Control","D7","D14")]
  cnt$gene_id <- rownames(cnt)
  cnt <- cnt[,c("gene","gene_id","Control","D7","D14")]
  cnt[,3:5] <- round(log2(cnt[,3:5]+1),2)
  names(cnt)[3:5] <- paste0(names(cnt[3:5]),"_log2Expn")
  rm(out0,out1,outC,ggg)
  
  load(paste0(DATADIR,"33_EditingLevelReport_Anno.RData"))
  z <- imputed_edits
  z$Control_edit <- round(apply(z[,7:9],1,gmean),2)
  z$D7_edit <- round(apply(z[,1:3],1,gmean),2)
  z$D14_edit <- round(apply(z[,4:6],1,gmean),2)
  z <- z[,10:12]
  z$id <- rownames(z)
  
  load(paste0(DATADIR,"genes.RData"))
  gGenes <- gGenes[,c(1:6,11)] %>% unique
  genes <- with(gGenes,GRanges(chr,IRanges(geneStart,geneEnd),strand,gene,gene_id,biotype))
  
  uq <- as(hc_edits_anno,"GRanges")
  uq <- unique(uq[,"id"])
  o <- findOverlaps(uq,genes,ignore.strand=T,select = "all")
  uq <- cbind(as.data.frame(uq[queryHits(o)]),as.data.frame(genes[subjectHits(o)])[,c("gene","gene_id","biotype")] )
  uq <- uq[!is.na(uq$gene_id),]
  uq <- unique(uq[,c("id","gene","gene_id","biotype")])
  
  h <- hlpr[,c("idx", "known", "repeats", "aaswap", "UCNE_name", "UCNE_ID", "snpeff_uq", 
               "cHMM_Final", "Control_Wt1", "D9_Wt1","D14_Wt1")] %>% unique
  h$id <- gsub(":A>G|:T>C","",h$idx)
  z <- merge(z,h,by="id")
  rm(h)
  
  z <- merge(z,uq,by="id",all.x=T)
  cnt$gene <- NULL
  z <- merge(z,cnt,by="gene_id",all.x=T)
  
  z <- z[,c("idx", "known", "repeats", "gene", "gene_id","biotype", "snpeff_uq",
            "Control_edit", "D7_edit", "D14_edit", "Control_log2Expn", "D7_log2Expn", "D14_log2Expn",
            "Control_Wt1", "D9_Wt1", "D14_Wt1",
            "aaswap", "UCNE_name", "UCNE_ID","cHMM_Final" )]
  
  z <- merge(z,unique(hc_edits_anno[,c("idx","ConsScore","hg19")]),by="idx")
  z$clusrd_in_Control <- z$clusrd_in_D7 <- z$clusrd_in_D14 <- 0
  
  h <- hlpr
  h$id <- gsub(":A>G|:T>C","",h$idx)
  h$chr <- gsub(":\\S*","",h$id)
  h$start <- h$end <- as.numeric(gsub("\\S*:","",h$id))
  h <- h[,c("chr","start","end","id","control","D7","D14","idx")]
  h <- as(h,"GRanges")
  h$idx <- as.character(h$idx)
  
  d2 <- reduce(resize(h[h$control>0],101,"center"))
  d2 <- resize(d2,width(d2)-100,"center")
  d2 <- d2[width(d2)>50]
  h1 <- h[h$control>0]
  h1$x <- countOverlaps(h1,d2,ignore.strand=T)
  h1 <- h1[h1$x>0]
  z[z$idx%in%h1$idx,]$clusrd_in_Control <- 1
  
  d2 <- reduce(resize(h[h$D7>0],101,"center"))
  d2 <- resize(d2,width(d2)-100,"center")
  d2 <- d2[width(d2)>50]
  h1 <- h[h$control>0]
  h1$x <- countOverlaps(h1,d2,ignore.strand=T)
  h1 <- h1[h1$x>0]
  z[z$idx%in%h1$idx,]$clusrd_in_D7 <- 1
  
  d2 <- reduce(resize(h[h$D14>0],101,"center"))
  d2 <- resize(d2,width(d2)-100,"center")
  d2 <- d2[width(d2)>50]
  h1 <- h[h$control>0]
  h1$x <- countOverlaps(h1,d2,ignore.strand=T)
  h1 <- h1[h1$x>0]
  z[z$idx%in%h1$idx,]$clusrd_in_D14 <- 1
  
  
  load(paste0(RAWDIR,"ensembl_g2go.RData"))
  z <- merge(z,g2go,by="gene_id",all.x=T)
  qq <- read.delim(paste0(RAWDIR,"Ensembl2MGI.txt"),header=T,stringsAsFactors = F)
  names(qq) <- c("gene_id","link")
  qq$link <- paste0("=HYPERLINK(\"http://www.informatics.jax.org/marker/",qq$link,"\",\"",qq$gene_id,"\")")
  z <- merge(z,qq,by="gene_id",all.x=T)
  
  z <- z[,c("idx", "known", "repeats", "gene", "gene_id", "biotype", "snpeff_uq", 
            "Control_edit", "D7_edit", "D14_edit", "Control_log2Expn", "D7_log2Expn", 
            "D14_log2Expn", "Control_Wt1", "D9_Wt1", "D14_Wt1", "aaswap", 
            "UCNE_name", "UCNE_ID", "cHMM_Final", "ConsScore", "hg19", "clusrd_in_Control", 
            "clusrd_in_D7", "clusrd_in_D14","GOAnno","link")]
  
  SPN <- z
  save(SPN,file=paste0(DATADIR,"33_Spreadsheet_HCEditsWithCompleteAnno.RData"))
  
}

### coplete editing annotation spreadsheet (all replicates)
if(T){
  load(paste0(DATADIR,"DEExpression_DESeq2.RData"))
  ggg <- out0[[1]][,c("id","geneStart","geneEnd","gene")]
  ggg$length <- abs(ggg$geneStart-ggg$geneEnd)
  ggg <- ggg[,c("id","length","gene")] %>% unique
  cnt <- DESeq2::counts(outC[[2]],normalized=F) %>% as.data.frame
  cnt$id <- rownames(cnt)
  cnt <- merge(cnt,ggg,by="id",all.x=T)
  rownames(cnt) <- cnt$id
  cnt$id <- NULL
  cnt[,1:9] <- cnt[,1:9]*1000/cnt$lengthO
  cnt$length <- NULL
  cnt[,1:9] <- apply(cnt[,1:9],2,function(x){round(x*1e6/sum(x),2)})
  cnt <- as.data.frame(cnt)
  names(cnt)[1:9] <- unname(sampl[names(cnt)[1:9]])
  head(cnt)
  cnt$D7 <- round(apply(cnt[,1:3],1,gmean),2)
  cnt$D14 <- round(apply(cnt[,4:6],1,gmean),2)
  cnt$Control <- round(apply(cnt[,7:9],1,gmean),2)
  #cnt <- cnt[,c("gene","Control","D7","D14")]
  cnt$gene_id <- rownames(cnt)
  cnt <- cnt[,c("gene_id","gene","Control_R1","Control_R2", "Control_R3", 
                "D7_R1", "D7_R2", "D7_R3", "D14_R1", "D14_R2", "D14_R3", "D7", "D14", "Control")]
  cnt[,3:ncol(cnt)] <- round(log2(cnt[,3:ncol(cnt)]+1),2)
  names(cnt)[3:ncol(cnt)] <- paste0(names(cnt[3:ncol(cnt)]),"_log2Expn")
  rm(out0,out1,outC,ggg)
  
  load(paste0(DATADIR,"33_EditingLevelReport_Anno.RData"))
  z <- imputed_edits
  z$Control_edit <- round(apply(z[,7:9],1,gmean),2)
  z$D7_edit <- round(apply(z[,1:3],1,gmean),2)
  z$D14_edit <- round(apply(z[,4:6],1,gmean),2)
  z$id <- rownames(z)
  names(z)[1:9] <- paste0(names(z)[1:9],"_edit")
  
  load(paste0(DATADIR,"genes.RData"))
  gGenes <- gGenes[,c(1:6,11)] %>% unique
  genes <- with(gGenes,GRanges(chr,IRanges(geneStart,geneEnd),strand,gene,gene_id,biotype))
  
  uq <- as(hc_edits_anno,"GRanges")
  uq <- unique(uq[,"id"])
  o <- findOverlaps(uq,genes,ignore.strand=T,select = "all")
  uq <- cbind(as.data.frame(uq[queryHits(o)]),as.data.frame(genes[subjectHits(o)])[,c("gene","gene_id","biotype")] )
  uq <- uq[!is.na(uq$gene_id),]
  uq <- unique(uq[,c("id","gene","gene_id","biotype")])
  
  h <- hlpr[,c("idx", "known", "repeats", "aaswap", "UCNE_name", "UCNE_ID", "snpeff_uq", 
               "cHMM_Final", "Control_Wt1", "D9_Wt1","D14_Wt1")] %>% unique
  h$id <- gsub(":A>G|:T>C","",h$idx)
  z <- merge(z,h,by="id")
  rm(h)
  
  z <- merge(z,uq,by="id",all.x=T)
  cnt$gene <- NULL
  z <- merge(z,cnt,by="gene_id",all.x=T)
  z <- merge(z,unique(hc_edits_anno[,c("idx","ConsScore","hg19")]),by="idx")
  
  z$clusrd_in_Control <- z$clusrd_in_D7 <- z$clusrd_in_D14 <- 0
  
  h <- hlpr
  h$id <- gsub(":A>G|:T>C","",h$idx)
  h$chr <- gsub(":\\S*","",h$id)
  h$start <- h$end <- as.numeric(gsub("\\S*:","",h$id))
  h <- h[,c("chr","start","end","id","control","D7","D14","idx")]
  h <- as(h,"GRanges")
  h$idx <- as.character(h$idx)
  d2 <- reduce(resize(h[h$control>0],101,"center"))
  d2 <- resize(d2,width(d2)-100,"center")
  d2 <- d2[width(d2)>50]
  h1 <- h[h$control>0]
  h1$x <- countOverlaps(h1,d2,ignore.strand=T)
  h1 <- h1[h1$x>0]
  z[z$idx%in%h1$idx,]$clusrd_in_Control <- 1
  d2 <- reduce(resize(h[h$D7>0],101,"center"))
  d2 <- resize(d2,width(d2)-100,"center")
  d2 <- d2[width(d2)>50]
  h1 <- h[h$control>0]
  h1$x <- countOverlaps(h1,d2,ignore.strand=T)
  h1 <- h1[h1$x>0]
  z[z$idx%in%h1$idx,]$clusrd_in_D7 <- 1
  d2 <- reduce(resize(h[h$D14>0],101,"center"))
  d2 <- resize(d2,width(d2)-100,"center")
  d2 <- d2[width(d2)>50]
  h1 <- h[h$control>0]
  h1$x <- countOverlaps(h1,d2,ignore.strand=T)
  h1 <- h1[h1$x>0]
  z[z$idx%in%h1$idx,]$clusrd_in_D14 <- 1
  
  dput(names(z))
  
  load(paste0(RAWDIR,"ensembl_g2go.RData"))
  z <- merge(z,g2go,by="gene_id",all.x=T)
  qq <- read.delim(paste0(RAWDIR,"Ensembl2MGI.txt"),header=T,stringsAsFactors = F)
  names(qq) <- c("gene_id","link")
  qq$link <- paste0("=HYPERLINK(\"http://www.informatics.jax.org/marker/",qq$link,"\",\"",qq$gene_id,"\")")
  z <- merge(z,qq,by="gene_id",all.x=T)
  z <- z[,c("idx", 
            "Control_edit", "D7_edit","D14_edit", 
            "Control_log2Expn","D7_log2Expn", "D14_log2Expn",  
            "gene_id", "biotype", "gene",  
            "known", "repeats", "aaswap", "UCNE_name", "UCNE_ID", "snpeff_uq", "cHMM_Final", 
            "Control_Wt1", "D9_Wt1", "D14_Wt1", 
            "ConsScore", "hg19", "clusrd_in_D14","clusrd_in_D7", "clusrd_in_Control",
            "control_R3_edit","control_R2_edit", "control_R1_edit", 
            "D7_R1_edit", "D7_R2_edit", "D7_R3_edit", 
            "D14_R1_edit", "D14_R2_edit", "D14_R3_edit", 
            "Control_R1_log2Expn", "Control_R2_log2Expn", "Control_R3_log2Expn", 
            "D7_R1_log2Expn", "D7_R2_log2Expn", "D7_R3_log2Expn", 
            "D14_R1_log2Expn", "D14_R2_log2Expn", "D14_R3_log2Expn","GOAnno","link")]
  
  SPNF <- z
  save(SPNF,file=paste0(DATADIR,"33_Spreadsheet_HCEditsWithCompleteAnnoFull.RData"))
  
  
  
}

### write spreadsheets
if(F){
  load(paste0(DATADIR,"Spreadsheet_HCEditsWithCompleteAnnoFull.RData"))
  write.table(SPNF,file=paste0(DATADIR,"Sheet_HCEditsWithCompleteAnnoFull_report.txt"),quote = F,sep = "\t",row.names = F)
  
  load(paste0(DATADIR,"Spreadsheet_HCEditsWithCompleteAnno.RData"))
  write.table(SPN,file=paste0(DATADIR,"Sheet_HCEditsWithCompleteAnno_report.txt"),quote = F,sep = "\t",row.names = F)
  
  
  load(paste0(DATADIR,"33_Spreadsheet_HCEditsWithCompleteAnnoFull.RData"))
  write.table(SPNF,file=paste0(DATADIR,"33_Sheet_HCEditsWithCompleteAnnoFull_report.txt"),quote = F,sep = "\t",row.names = F)
  
  load(paste0(DATADIR,"33_Spreadsheet_HCEditsWithCompleteAnno.RData"))
  write.table(SPN,file=paste0(DATADIR,"33_Sheet_HCEditsWithCompleteAnno_report.txt"),quote = F,sep = "\t",row.names = F)
  
  
}

##################################################################################################
### GO analysis of the eidting sites
##################################################################################################
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
  rm(cnt,out1,out2,outC,genes)
  
  load(paste0(DATADIR,"33_Spreadsheet_HCEditsWithCompleteAnno.RData"))
  
  x <- SPN[SPN$Control_edit>0 ,]$gene_id %>% unique
  x <- gGenes[gGenes$gene_id%in%x,]$ucscID
  y <- SPN[SPN$D7_edit>0 ,]$gene_id %>% unique
  y <- gGenes[gGenes$gene_id%in%y,]$ucscID
  z <- SPN[SPN$D14_edit>0 ,]$gene_id %>% unique
  z <- gGenes[gGenes$gene_id%in%z,]$ucscID
  out <- rbind(
    cbind(mm9_GOStats.enrichment(x,p.val = 0.1,universe = exp.trx),comparison="Control"),
    cbind(mm9_GOStats.enrichment(y,p.val = 0.1,universe = exp.trx),comparison="D7"),
    cbind(mm9_GOStats.enrichment(z,p.val = 0.1,universe = exp.trx),comparison="D14")
  )
  save(out,file=paste0(DATADIR,file="33_GOEnrichment_AllEditing.RData"))
  
  
  
}





## intron, exon annotations
if(F){
  myf <- function(path="A:/work/mm10_Annotations/mm10_3UTR"){
    z <- read.delim(path,header = F)
    z$V1 <- gsub("chr","",z$V1)
    z <- z[,c(1:3,6)]
    z <-with(z,GRanges(V1,IRanges(V2,V3),V6)) 
    z <- reduce(z,ignore.strand=T)
    z <- subset(z,seqnames(z)%in%c("X","Y",paste0("",1:22)))
    z
  }
  utr3 <- myf(UTR3)
  utr5 <- myf(UTR5)
  introns <- myf(INTRON)  
  exons <- myf(EXON)
  
  
  load(paste0(DATADIR,"33_EditingLevelReport_Anno.RData"))
  
  df1 <- as(hc_edits_anno,"GRanges")
  df1$utr3 <- countOverlaps(df1,utr3,ignore.strand=T)
  df1$utr5 <- countOverlaps(df1,utr5,ignore.strand=T)
  df1$exon <- countOverlaps(df1,exons,ignore.strand=T)
  df1$intron <- countOverlaps(df1,introns,ignore.strand=T)
  df1 <- as.data.frame(df1)
  df1$anno <- NA
  df1[df1$utr3>0,]$anno <- "3UTR"
  df1[df1$utr5>0 & is.na(df1$anno),]$anno <- "5UTR"
  df1[df1$exon>0 & is.na(df1$anno),]$anno <- "Exon"
  df1[df1$intron>0 & is.na(df1$anno),]$anno <- "Intron"
  df1[is.na(df1$anno),]$anno <- "Intergenic"
  df1$utr3 <- df1$utr5 <- df1$exon <- df1$intron <- NULL
  
  pl <- as.data.frame(table(df1$anno,df1$sample_name))
  pl$condition <- gsub("_R[123]","",pl$Var2)
  pl <- data.table(pl)
  pl <- pl[,tot:=sum(Freq),by=list(Var2,condition)]
  pl$Freq <- round(pl$Freq/pl$tot,2)
  pl <- pl[,mean:=mean(Freq),by=list(Var1,condition)]
  pl <- pl[,sd:=sd(Freq),by=list(Var1,condition)]
  pl <- pl[,ci:=1.96*sem(Freq),by=list(Var1,condition)]
  pl <- as.data.frame(pl)
  pl <- pl[,c("condition","Var1","mean","sd","ci")] %>% unique
  pl$condition <- factor(pl$condition,levels=c("control","D7","D14"))
  
  gl <- as.data.frame(table(df1$anno,df1$sample_name))
  gl$condition <- gsub("_R[123]","",gl$Var2)
  gl <- data.table(gl)
  gl <- gl[,mean:=mean(Freq),by=list(Var1,condition)]
  gl <- gl[,sd:=sd(Freq),by=list(Var1,condition)]
  gl <- gl[,ci:=1.96*sem(Freq),by=list(Var1,condition)]
  gl <- as.data.frame(gl)
  gl <- gl[,c("condition","Var1","mean","sd","ci")] %>% unique
  gl$condition <- factor(gl$condition,levels=c("control","D7","D14"))
  
  pllist <- list()
  
  pllist[[1]] <-  ggplot(pl,aes(x=Var1,y=mean,fill=condition,alpha=0.2))+
    geom_bar(stat="identity", color="black", position=position_dodge())+
    geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.2,position=position_dodge(.9))+
    xlab("")+ylab("fraction of editing sites")+
    scale_fill_manual(values = cols)+gg_aes+
    scale_y_continuous(expand = c(0,0))+
    theme(legend.position = "bottom")
  
  pllist[[2]] <-ggplot(gl,aes(x=Var1,y=mean,fill=condition,alpha=0.2))+
    geom_bar(stat="identity", color="black", position=position_dodge())+
    geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.2,position=position_dodge(.9))+
    xlab("")+ylab("number of editing sites")+
    scale_fill_manual(values = cols)+gg_aes+
    scale_y_continuous(expand = c(0,0))+
    theme(legend.position = "bottom")
  
  #pdf(paste0(FIGDIR,"Proportions of Editing cases_LOI.pdf"),width = 12,height = 6)
  grid.arrange(grobs=pllist,nrow=1)
  #dev.off()
  
  
}
## 


