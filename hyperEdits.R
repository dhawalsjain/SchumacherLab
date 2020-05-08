rm(list=ls())
source("functions.R")
source("vars.R")


############################################################################3

if(F){
  libs$sample <- unlist(l1[as.character(libs$profile)])
  libs$condition <- gsub("_\\S*","",libs$sample)
  libs$sample <- factor(libs$sample,levels=c("control_R1", "control_R2", "control_R3",
                                             "D7_R1", "D7_R2", "D7_R3",
                                             "D14_R1", "D14_R2", "D14_R3"))
  #libs$unmap <- libs$unmap/1e6
  
  ggplot(libs, aes(x=sample,y=unmap/1e6,fill=condition,col=condition,alpha=0.2))+
    geom_bar(stat = "identity",position="dodge")+
    xlab("")+ylab("number of unmapped mates \n in Millions")+
    scale_fill_manual(values = cols)+scale_color_manual(values=cols)+
    theme_bw()+gg_aes+
    theme(axis.text.x = element_text(angle = 90,hjust = 1))
  
  
  ggplot(libs, aes(x=sample,y=hyperEdits/(depths*1e6),fill=condition,col=condition,alpha=0.2))+
    geom_bar(stat = "identity",position="dodge")+
    xlab("")+ylab("number of hyperedited reads\n per Million")+
    scale_fill_manual(values = cols)+scale_color_manual(values=cols)+
    theme_bw()+gg_aes+
    theme(axis.text.x = element_text(angle = 90,hjust = 1))
  
  
}


if(F){
  files <- Sys.glob(paste0(RAWDIR,"hyperEdits/*.txt"))
  
  
  he <- GRanges()
  for(f in files){
    cat(f,"\n")
    d <- read.delim(f,header=F,stringsAsFactors = F)
    names(d) <- c("chr","start","Ref","Alt","id","mapq")
    d <- d[d$mapq>=10,]
    d$profile <- gsub(".a2g.sort.report.txt","",gsub("hyperEdits/","",gsub(RAWDIR,"",f)))
    d <- data.table(d)
    d <- unique(d)
    d <- as.data.frame(d)
    d <- with(d,GRanges(chr,IRanges(start,start),"*",id,Ref,Alt,mapq,profile))
    he <- c(he,d)
    rm(d,f)
  }
  
  he <- as.data.frame(he)
  he$seqnames <- paste0("chr",he$seqnames)
  he <- he[!he$seqnames%in%c("chrJH584304.1","chrGL456221.1","chrGL456383.1","chrMT"),]
  he <- as(he,"GRanges")
  library(BSgenome.Mmusculus.UCSC.mm10)
  genome <- BSgenome.Mmusculus.UCSC.mm10
  
  sq <- getSeq(genome,he)
  he$REF <- as.character(sq)
  rm(sq)
  he
  he <- he[he$REF%in%c("T","A")]
  
  he <- as.data.frame(he)
  libs$sample <- unlist(l1[as.character(libs$profile)])
  libs$condition <- gsub("_\\S*","",libs$sample)
  he <- merge(he,libs[,c("profile","condition","sample")],by="profile")
  
  x <- he[,c("id","condition","sample")]
  x$ch <- 1
  head(x)
  x <- data.table(x)
  x <- x[,tot:=sum(ch),by=list(condition,id)]
  y <- as.data.frame(x)
  y <- y[,c("id","condition","tot")] %>% unique
  
  ggplot(y, aes(x=tot,col=condition,fill=condition)) +geom_density(alpha=0.2)+
    xlim(0,100)+
    scale_fill_manual(values = cols)+
    scale_color_manual(values = cols)
  
  
  
  
  
  
  
  
  
  
}




############################################################################