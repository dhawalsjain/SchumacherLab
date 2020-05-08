
####################################################################################################
#### base functions
####################################################################################################

### 1) Filter the editing calls for presence in 2/3 replicates
filt.function <- function(x){
  return(ifelse(sum(is.na(x[1:3])) <=CUTOFFLEN | sum(is.na(x[4:6])) <=CUTOFFLEN |sum(is.na(x[7:9]))<=CUTOFFLEN,1,0))
}

## 2) Impute the data 
impute.data <- function(z1){
  library(mice)
  z1$test <- apply(z1,1,filt.function)
  z1 <- z1[z1$test>0,]
  z0 <- mice(z1[,1:3],m = 5,maxit = 30)
  z2 <- complete(z0)
  z2$id <- rownames(z1)
  z2 <- reshape::melt(z2,measure.var=names(z2)[1:3])
  return(z2)
}

#### 3) Geometric mean
gmean <- function(x, na.rm=TRUE){
  exp(sum(log(x), na.rm=na.rm) / length(x))
}

### 4) GO enrichment
mm9_GOStats.enrichment <- function(genes,p.val,universe){
  require("GOstats")
  require("GSEABase")
  
  #d <- read.delim("A:/work/mm9_Annotations/mm9_MGI2ucsc.txt")
  #d <- unique(d[,c(2,3)])
  #b <- read.delim("A:/work/mm9_Annotations/mgigene_association.mgi",comment.char = "!",header=F)
  #b <- unique(b[,c(2,5,7)])
  #names(d)
  #names(b)[1] <- "MGI.ID"
  #d <- merge(d,b,by="MGI.ID",all.x=T)
  #d <- unique(d[,c(2:4)])
  #d$UCSC.ID <- ifelse(d$UCSC.ID=="",NA,as.character(d$UCSC.ID))
  #d <- d[complete.cases(d),]
  #write.table(d,file="A:/work/mm9_Annotations/mm9_GOTerms.txt",quote = F,sep = "\t",row.names = F)
  #frame <- read.delim("A:/work/mm9_Annotations/mm9_GOTerms.txt",header=T)
  #frame <- frame[complete.cases(frame),]
  #goframeData = frame[,c(2,3,1)]
  #names(goframeData) = c("go_id", "Evidence", "gene_id")
  #frame <- frame[!frame$V7%in%c("ND","IC","TAS","NAS"),]
  #goFrame=GOFrame(goframeData,organism="Mus musculus")
  #goAllFrame=GOAllFrame(goFrame)
  #gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
  #save(gsc,file="A:/work/mm9_Annotations/mm9 gene ontologies")
  
  load("A:/work/mm9_Annotations/mm9 gene ontologies")
  
  params <- GSEAGOHyperGParams(name="My Custom GSEA based annot Params",geneSetCollection=gsc,
                               geneIds = genes,universeGeneIds = universe,
                               ontology = "MF",pvalueCutoff = p.val,
                               conditional = TRUE,testDirection = "over")
  mf <- summary(hyperGTest(params))
  if(length(mf$Pvalue)>1){
    mf <- cbind(mf,Category="Molecular Function")
    names(mf)[1] <- c("ID")
  }
  
  params <- GSEAGOHyperGParams(name="My Custom GSEA based annot Params",geneSetCollection=gsc,
                               geneIds = genes,universeGeneIds = universe,
                               ontology = "BP",pvalueCutoff = p.val,
                               conditional = TRUE,testDirection = "over")
  bp <- summary(hyperGTest(params))
  if(length(bp$Pvalue)>1){
    bp <- cbind(bp,Category="Biological Process")
    names(bp)[1] <- c("ID")
  }
  
  params <- GSEAGOHyperGParams(name="My Custom GSEA based annot Params",geneSetCollection=gsc,
                               geneIds = genes,universeGeneIds = universe,
                               ontology = "CC",pvalueCutoff = p.val,
                               conditional = TRUE,testDirection = "over")
  cc <- summary(hyperGTest(params))
  if(length(cc$Pvalue)>1){
    cc <- cbind(cc,Category="Cellular Component")
    names(cc)[1] <- c("ID")
  }
  
  if(length(mf$Pvalue)==0){mf = data.frame(ID=NA,Pvalue=NA,OddsRatio=NA,ExpCount=NA,Count=NA,Size=NA,Term=NA,Category=NA)}
  if(length(bp$Pvalue)==0){bp = data.frame(ID=NA,Pvalue=NA,OddsRatio=NA,ExpCount=NA,Count=NA,Size=NA,Term=NA,Category=NA)}
  if(length(cc$Pvalue)==0){cc = data.frame(ID=NA,Pvalue=NA,OddsRatio=NA,ExpCount=NA,Count=NA,Size=NA,Term=NA,Category=NA)}
  g <- rbind(mf,bp,cc)
  g <- g[complete.cases(g),]
  g <- g[order(g$Pvalue),]
  
  return(g)
}

### 5) DESeq2 expression data frame
myDeSeq2 <- function(df2,type,condition,batch=NULL){
  require("DESeq2")
  
  if(!is.null(batch)){
    colData <- data.frame(condition,type,batch)
    row.names(colData) <- colData$type
    dds <- DESeqDataSetFromMatrix(countData = df2,colData = colData,design = ~ batch+condition)
  }else{
    colData <- data.frame(condition,type)
    row.names(colData) <- colData$type
    dds <- DESeqDataSetFromMatrix(countData = df2,colData = colData,design = ~ condition)
  }
  
  dds <- DESeq(dds)
  rld <- rlog(dds)
  res <- results(dds) %>% as.data.frame
  res$id <- rownames(res)
  j = names(dds@rowRanges@elementMetadata@listData)[12]
  j <- gsub("condition_","",j)
  j = unlist(strsplit(j,"_vs_"))
  res$log2FoldChange <- res$log2FoldChange*(-1)
  names(res)[2] <- paste0("log2 (",j[2],"/",j[1],")")
  
  return(list(res=res, object=dds,rlogs=rld))
}

## adopted function from some random online source
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

## 6) Venn
Venn2 <- function(set1, set2, names,plot=TRUE,return.df=T) {
  require(limma)
  stopifnot( length(names) == 2)
  set1 <- unique(set1)
  set2 <- unique(set2)
  
  # Form universe as union of all three sets
  universe <- sort( unique( c(set1, set2) ) )
  
  Counts <- matrix(0, nrow=length(universe), ncol=2)
  colnames(Counts) <- names
  
  for (i in 1:length(universe))
  {
    Counts[i,1] <- universe[i] %in% set1
    Counts[i,2] <- universe[i] %in% set2
  }
  
  #if(plot==TRUE){ 
  vennDiagram( vennCounts(Counts),circle.col = c(rgb(1,0,0,0.2),rgb(0,0,1,0.2)),cex=c(1.5,1.4,1.4)  ) 
  #}
  
  rownames(Counts) <- universe
  
  df1 <- list(AB=rownames(subset(Counts,Counts[,1]==1 & Counts[,2]==1)),
              A=rownames(subset(Counts,Counts[,1]==1 & Counts[,2]==0)),
              B=rownames(subset(Counts,Counts[,1]==0 & Counts[,2]==1))
  )
  
  names(df1) <- c(paste(names[1],names[2]), names[1],names[2])
  df2 <- c()
  for (name in names(df1)) {
    df2 <- rbind(df2, data.frame(condition=rep(name,length(df1[[name]])),tracking_id=df1[[name]]))
  }
  df2 <- unique(df2)
  if(return.df==T){
    return(df2)
  }
}

## 7) Vector resize
vector.resizing <- function(x,final.len=100){
  y <- vector()
  len <- length(x)
  y <-spline(1:len,x,n=final.len)$y
  return(y)
}

## 8) Dist2next
get.distance2next <- function(k){
  # k is a data.frame with 3 column headers: chr, start, end
  # additional columns are permisible
  
  k <- k[order(k$chr,k$start),]
  k$dist2next <- 0
  for(f in levels(as.factor(k$chr))){
    e <- k[k$chr==f,]$end
    s <- k[k$chr==f,]$start
    
    if(length(s)>1){
      s <- s[2:length(s)]
      s <- c(s,e[length(e)])
      k[k$chr==f,]$dist2next <- s-e
    }else if(length(s)==1){
      k[k$chr==f,]$dist2next <- 0
    }
    rm(s,e,f)
  }
  return(k)
}


sem <- function(x){
  sd(x, na.rm=TRUE) /  sqrt(length(x[!is.na(x)]))
}


## 