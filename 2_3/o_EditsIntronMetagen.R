rm(list=ls())
source("functions.R")
source("vars.R")


###########################################################################3
INTRON="A:/work/mm10_Annotations/mm10_intorns"
EXON="A:/work/mm10_Annotations/mm10_exons"
UTR3="A:/work/mm10_Annotations/mm10_3UTR"
UTR5="A:/work/mm10_Annotations/mm10_5UTR"

if(F){
  myf <- function(path="A:/work/mm10_Annotations/mm10_3UTR"){
    z <- read.delim(path,header = F)
    z$V1 <- gsub("chr","",z$V1)
    z <- z[,c(1:4,6)]
    z <-with(z,GRanges(V1,IRanges(V2,V3),V6,id=V4)) 
    z <- subset(z,seqnames(z)%in%c("X","Y",paste0("",1:22)))
    z
  }
  myf_cnt <- function(gp,introns,WIDTH=100){
    vector.resizing <- function(x,final.len=100){
      y <- vector()
      len <- length(x)
      y <-spline(1:len,x,n=final.len)$y
      return(y)
    }
    
    OUTBIN=20
    OUTWIN=500
    OUTSTEPS=OUTWIN/OUTBIN
    
    #intp <- subsetByOverlaps(introns,gp,ignore.strand=T)
    intp <- introns
    orig <- intp
    strand(intp) <- "+"
    up <- flank(intp,OUTWIN)
    dn <- flank(intp,OUTWIN,start = F)
    
    intpl <- intp[width(intp)<WIDTH]
    intp <- intp[width(intp)>=WIDTH]
    id <- rep(intp$id,each=WIDTH)
    intp <- unlist(tile(intp,n = WIDTH))
    intp$id <- id
    intp$re <- countOverlaps(intp,gp,ignore.strand=T)
    
    id <- rep(as.character(intpl$id),width(intpl))
    intpl <- unlist(tile(intpl,width = 1))
    intpl$id <- id
    intpl$re <- countOverlaps(intpl,gp,ignore.strand=T)
    intpl <- data.frame(id=intpl$id, cnt=intpl$re)
    intpl <- data.table(intpl)
    intpl <- intpl[,lapply(.SD,vector.resizing),by=id] ## 
    intpl <- as.data.frame(intpl)   
    intpl$cnt <- ifelse(intpl$cnt<0,0,intpl$cnt)
    intpl$cnt <- round(intpl$cnt,2)
    
    id <- rep(as.character(up$id),each=OUTSTEPS)
    up <- unlist(tile(up,n = OUTSTEPS))
    up$id <- id
    up$re <- countOverlaps(up,gp,ignore.strand=T)
    
    id <- rep(as.character(dn$id),each=OUTSTEPS)
    dn <- unlist(tile(dn,n = OUTSTEPS))
    dn$id <- id
    dn$re <- countOverlaps(dn,gp,ignore.strand=T)
    
    m <- matrix(intpl$cnt,ncol = WIDTH,byrow = T)
    rownames(m) <- unique(intpl$id)
    n <- matrix(intp$re,ncol = WIDTH,byrow = T)
    rownames(n) <- unique(intp$id)
    m <- rbind(m,n)
    rm(n)
    
    u <- matrix(up$re,ncol = OUTSTEPS,byrow = T)
    rownames(u) <- unique(up$id)
    d <- matrix(dn$re,ncol = OUTSTEPS,byrow = T)
    rownames(d) <- unique(dn$id)
    
    d <- d[match(rownames(u),rownames(d)),]
    m <- m[match(rownames(u),rownames(m)),]
    m <- cbind(u,m,d)
    
    n <- subset(m,rownames(m)%in%orig[strand(orig)=="-"]$id)
    m <- subset(m,!rownames(m)%in%orig[strand(orig)=="-"]$id)
    n <- apply(n,1,rev) %>% t
    m <- rbind(m,n)
    return(as.data.frame(m))
  }
  myf_centered <- function(gp,introns,WIDTH=1001, bin=20){
    
    #intp <- subsetByOverlaps(introns,gp,ignore.strand=T)
    intp <- introns
    orig <- intp
    
    s <- resize(resize(intp,1,"start"),WIDTH,"center")
    e <- resize(resize(intp,1,"end"),WIDTH,"center")
    
    STEPS = round(WIDTH/bin)
    
    id <- rep(s$id,each=STEPS)
    s <- unlist(tile(s,n = STEPS))
    s$id <- id
    s$re <- countOverlaps(s,gp,ignore.strand=T)
    e <- unlist(tile(e,n = STEPS))
    e$id <- id
    e$re <- countOverlaps(e,gp,ignore.strand=T)
    
    m <- matrix(s$re,ncol = STEPS,byrow = T)
    rownames(m) <- unique(s$id)
    n <- matrix(e$re,ncol = STEPS,byrow = T)
    rownames(n) <- unique(e$id)
    
    mp <- subset(m,rownames(m)%in%orig[strand(orig)=="+"]$id)
    mn <- subset(m,!rownames(m)%in%orig[strand(orig)=="+"]$id)
    mn <- apply(mn,1,rev) %>% t
    m <- rbind(mp,mn)
    
    np <- subset(n,rownames(n)%in%orig[strand(orig)=="+"]$id)
    nn <- subset(n,!rownames(n)%in%orig[strand(orig)=="+"]$id)
    nn <- apply(nn,1,rev) %>% t
    n <- rbind(np,nn)
    
    m <- apply(m,1,function(x) if(sum(x)>0){round(x/sum(x),2)}else{x}) %>% t %>% as.data.frame
    n <- apply(n,1,function(x) if(sum(x)>0){round(x/sum(x),2)}else{x}) %>% t %>% as.data.frame
    
   l <- list(start=m,end=n)
   return(l)  
}
  sem <- function(x){
    x <- x[complete.cases(x)]
    1.96*sd(x,na.rm = T)/sqrt(length(x))
  }
  myf_GeneratePlotDf <- function(uq){
    pllist <- list()
    for(ff in c(INTRON,EXON,UTR3,UTR5)){
      cat(ff,"\n")
      ###
      introns <- myf(path = ff)  
      introns <- subsetByOverlaps(introns,resize(uq,2001,"center"),ignore.strand=T)
      
      l <- myf_centered(gp = unique(uq[grep("control",uq$sample_name)]),introns = introns,WIDTH = 1001,bin = 5)
      l1 <- myf_centered(gp = unique(uq[grep("D7",uq$sample_name)]),introns = introns,WIDTH = 1001,bin = 5)
      l2 <- myf_centered(gp = unique(uq[grep("D14",uq$sample_name)]),introns = introns,WIDTH = 1001,bin = 5)
      
      
      pl <- rbind(
        data.frame(pos=1:ncol(l[[1]]), mean=apply(l[[1]],2,mean), se=apply(l[[1]],2,sem),loc="start",profile="control"),
        data.frame(pos=1:ncol(l1[[1]]), mean=apply(l1[[1]],2,mean), se=apply(l1[[1]],2,sem),loc="start",profile="D7"),
        data.frame(pos=1:ncol(l2[[1]]), mean=apply(l2[[1]],2,mean), se=apply(l2[[1]],2,sem),loc="start",profile="D14"),
        data.frame(pos=1:ncol(l[[2]]), mean=apply(l[[2]],2,mean), se=apply(l[[2]],2,sem),loc="end",profile="control"),
        data.frame(pos=1:ncol(l1[[2]]), mean=apply(l1[[2]],2,mean), se=apply(l1[[2]],2,sem),loc="end",profile="D7"),
        data.frame(pos=1:ncol(l2[[2]]), mean=apply(l2[[2]],2,mean), se=apply(l2[[2]],2,sem),loc="end",profile="D14")
      )
      
      ff <- gsub("A:/work/mm10_Annotations/mm10_","",ff)
      pllist[[ paste0(ff," (",nrow(l[[1]]),")") ]] <- pl
      rm(ff,l1,l,l2,pl,introns)
    }
    return(pllist)
  }
  
  if(F){
    m <- myf_cnt(gp = unique(uq[grep("control",uq$sample_name)]),introns = introns)
    m <- apply(m,1,function(x) if(sum(x)>0){round(x*100/sum(x),2)}else{x}) %>% t %>% as.data.frame
    plot(1:ncol(m),apply(m,2,mean))
    abline(v = c(25,125))
    
    m <- myf_cnt(gp = unique(uq[grep("D7",uq$sample_name)]),introns = introns)
    m <- apply(m,1,function(x) if(sum(x)>0){round(x*100/sum(x),2)}else{x}) %>% t %>% as.data.frame
    points(1:ncol(m),apply(m,2,mean),col="red",pch=16)
    
    m <- myf_cnt(gp = unique(uq[grep("D14",uq$sample_name)]),introns = introns)
    m <- apply(m,1,function(x) if(sum(x)>0){round(x*100/sum(x),2)}else{x}) %>% t %>% as.data.frame
    points(1:ncol(m),apply(m,2,mean),col="green",pch=20)
  }
  
  load(paste0(DATADIR,"EditingLevelReport_Anno.RData"))
  uq <- as(hc_edits_anno,"GRanges")
  pllist <- myf_GeneratePlotDf(uq)
  save(pllist,file=paste0(DATADIR,"EditingMetagenePlots_IntronExonUTRs.RData"))
  
  
  uq <- read.delim(paste0(DATADIR,"CombinedEditDB.bed"),header = F,stringsAsFactors = F)
  uq <- with(uq,GRanges(V1,IRanges(V2,V3)))
  uq$sample_name <- "control"
  pllist <- myf_GeneratePlotDf(uq)
  save(pllist,file=paste0(DATADIR,"EditingMetagenePlots_IntronExonUTRs_DBCalls.RData"))
  
  
  load(paste0(DATADIR,"DRE_DE_matrix.RData"))
  res <- as(res,"GRanges")
  res <- res[res$adj.P.Val<0.05]
  
  uq <- res[res$comparison=="D7_control"]
  uq$sample_name <- "control"
  uq[uq$logFC>0]$sample_name <- "D7"
  pllist <- myf_GeneratePlotDf(uq)
  save(pllist,file=paste0(DATADIR,"EditingMetagenePlots_IntronExonUTRs_D7vsControl.RData"))
  
  uq <- res[res$comparison=="D14_control"]
  uq$sample_name <- "control"
  uq[uq$logFC>0]$sample_name <- "D14"
  pllist <- myf_GeneratePlotDf(uq)
  save(pllist,file=paste0(DATADIR,"EditingMetagenePlots_IntronExonUTRs_D14vsControl.RData"))
  

}



















