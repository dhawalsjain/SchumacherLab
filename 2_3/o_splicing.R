rm(list=ls())
source("functions.R")
source("vars.R")

###########################################3
### splicing analyses

########## Differential splicing with Suppa (generate expression file)
if(F){
  rm(list=ls())
  setwd(paste0(RAWDIR,"/HTSeq/star/"))
  files <- Sys.glob("*_star_1.txcounts_nua.txt")
  
  tx <-read.delim(paste0(RAWDIR,"/HTSeq/star/mm10_GRCm38ensemblTx.gtf"),header=F,stringsAsFactors = F)
  tx$tid <- gsub(";.*","",gsub(".*transcript_id ","",tx$V9)) %>% as.character
  tx$width <- abs(tx$V4-tx$V5)
  tx <- tx[,c("tid","width")] %>% unique
  
  cnt <- c()
  for(f in files){
    d <- read.delim(f,header = F,stringsAsFactors = F,comment.char = "_")
    cnt <- cbind(cnt,d$V3)
  }
  cnt <- as.data.frame(cnt)
  names(cnt) <- gsub("_star_1.txcounts_nua.txt","",files)
  cnt$tid <- d$V1
  
  cnt <- merge(cnt,tx, by ="tid",all.x=T)
  m <- cnt[,2:10]
  m[is.na(m)] <- 0
  cnt[,2:10] <- m
  cnt[,2:10] <- cnt[,2:10]*1000/cnt[,11]
  cnt <- cnt[complete.cases(cnt),]
  cnt[,2:10] <- apply(cnt[,2:10],2, function(x) { round(x*1e6/sum(x),2) })
  cnt$width <- NULL
  #write.table(cnt,file="Tx_TPM_expression_STAR.txt",sep = "\t",row.names = F,quote = F)
  
  write.table(cnt[,1:4],file="Tx_TPM_expression_STAR_D7.txt",sep = "\t",row.names = F,quote = F)
  write.table(cnt[,c(1,5:6,10)],file="Tx_TPM_expression_STAR_Control.txt",sep = "\t",row.names = F,quote = F)
  write.table(cnt[,c(1,7:9)],file="Tx_TPM_expression_STAR_D14.txt",sep = "\t",row.names = F,quote = F)
  
}

########## Differential splicing with Suppa (analyze output)
if(F){
  setwd(paste0(RAWDIR,"/HTSeq/star/"))
  files <- Sys.glob("*.dpsi")
  spl <- c()
  for(f in files){
    d <- read.delim(f,header = T,stringsAsFactors = F)
    names(d) <- c("change","pval")
    d$id <- rownames(d)
    d$test <- gsub("[.]dpsi","",f)
    spl <- rbind(spl,d)
    rm(f,d)
  }
  spl <- spl[is.finite(spl$change),]
  rownames(spl) <- NULL
  spl$ensembl <- gsub(";\\S*","",spl$id)
  spl$id <- gsub("\\S*;","",spl$id)
  spl$id <- gsub("-",":",spl$id)
  spl$num <- 1:nrow(spl)
  spl1 <- apply(spl,1,function(x){c(unlist(strsplit(x[3],":"))[2:6],x[c(1,2,4:6)]) } ) %>% t %>% as.data.frame
  spl1 <- melt(spl1,measure.vars = names(spl1)[2:5])
  spl1$value <- as.numeric(spl1$value)
  spl <- with(spl1,GRanges(id2,IRanges(value,value),"*",change,pval,test,ensembl,num))
  save(spl,file=paste0(DATADIR,"SplicingReport_SUPPA.RData"))
  
  

  
}

