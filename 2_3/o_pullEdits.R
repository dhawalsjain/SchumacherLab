source("functions.R")
source("vars.R")

###############################################################################################
### Combined pool of RNA editing sites
###############################################################################################
FILES=Sys.glob("A:/work/Kreidberg lab/Valerie/RNAEditing_Bowtie/*.editingSites.vcf")
INPUT="A:/work/Kreidberg lab/Valerie/editing_db/RediPortal.txt"
INPUT1="A:/work/Kreidberg lab/Valerie/editing_db/darned.txt"

if(!file.exists( paste0(DATADIR,"CombinedPool.bed")) ){
  
  ## our data
  a <- c()
  for(f in FILES){
    a <- rbind(a,read.delim(f,header = F,stringsAsFactors = F,comment.char = "#"))
  }
  
  ## Reduportal
  b <- read.delim(INPUT,header = T,comment.char = "#",stringsAsFactors = F)
  b <- b[b$Chr!="Chr",]
  b <- b[b$dbSNP=="-",]
  b$Chr <- gsub("chr","",b$Chr)
  
  
  ## DARNED
  d <- read.delim(INPUT1,header = T,stringsAsFactors = F,comment.char = "#")
  d$inrna <- ifelse(d$inrna==as.character("I"),as.character("G"),as.character(d$inrna))
  d$inrna <- ifelse(d$inrna==as.character("U"),as.character("T"),as.character(d$inrna))
  
  ## combined pool
  a$edit <- paste0(a$V4,">",a$V5)
  b$edit <- paste0(b$Ref,">",b$Ed)
  d$edit <- paste0(d$inchr,">",d$inrna)
  names(a)[1:2] <- names(b)[1:2] <- names(d)[1:2] <- c("chr","start")
  r <- rbind(a[,c(1:2,9)],b[,c(1:2,15)],d[,c(1:2,11)])
  r$end <- r$start
  r$score = 255
  r$strand = "+"
  r <- unique(r)
  r <- r[,c("chr","start","end","edit","score","strand")]
  write.table(r,file=paste0(DATADIR,"CombinedPool.bed"),quote = F,sep = "\t",row.names = F,col.names = F)
  
  ## known events
  q <- rbind(b[,c(1:2,15)],d[,c(1:2,11)])
  q$end <- q$start
  q$score = 255
  q$strand = "+"
  q <- unique(q)
  q <- q[,c("chr","start","end","edit","score","strand")]
  write.table(q,file=paste0(DATADIR,"CombinedEditDB.bed"),quote = F,sep = "\t",row.names = F,col.names = F)
  
  ## vcf for annotation
  ###CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	
  #1	4810987	.	T	C	17.61	.	ABHom=0;AC=2;AF=1;AN=2;BaseCounts=0,2,0,1;BaseQRankSum=0;DP=3;Dels=0;ExcessHet=3;FS=0;HaplotypeScore=0;MLEAC=2;MLEAF=1;MQ=34;MQ0=0;MQRankSum=0;OND=0;QD=5;ReadPosRankSum=0;SOR=0;GI=ENSMUSG00000025903:intron,ENSMUSG00000104217:intron;
  r1 <- r
  r1$REF <- gsub(">\\S*","",r1$edit)
  r1$ALT <- gsub("\\S*>","",r1$edit)
  r1$ID <- r1$QUAL <- r1$FILTER <- r1$INFO <- "."
  r1 <- r1[,c("chr","start","ID","REF","ALT","QUAL","FILTER","INFO")]
  write.table(r1,file=paste0(DATADIR,"CombinedPool.vcf"),quote = F,sep = "\t",row.names = F,col.names = F)
  rm(a,b,d,r,r1,q)
  
  ## snpeff annotation
  if(file.exists("CombinedPool.anno.vcf")){
    res <- read.delim("CombinedPool.anno.vcf",header = F,sep = "\t",comment.char = "#",stringsAsFactors = F)
    res$snpeff <- apply(res,1,function(x){
      q <- unique(unlist(strsplit(x[8],","))) ## annotation is the 8th column
      qq <- c()
      for(i in 1:length(q)){
        x <- unlist(strsplit(q[i],"\\|"))[7]
        y <- unlist(strsplit(q[i],"\\|"))[2]
        qq <- c(qq,paste0(y,"(",x,")"))
        rm(x,y)
      }
      qq <- gsub("_variant","",qq)
      qq <- paste0(unique(qq),collapse = ";")
      qq    
    })
    res$aaswap <- apply(res,1,function(x){
      q <- unique(unlist(strsplit(x[8],",")))
      qq <- c()
      for(i in 1:length(q)){
        x <- unlist(strsplit(q[i],"\\|")) [7]
        y <- gsub("p[.]","",unlist(strsplit(q[i],"\\|")) [11])
        if(y!="" & !is.na(y)){
          qq <- c(qq,paste0(y,"(",x,")"))
        } 
        rm(x,y,i)
      }
      if(is.null(qq)){
        return("NA")
      }else{
        qq <- paste0(unique(qq),collapse = ";")
        return(qq)
      }
    })
    r$snpeff <- res$snpeff
    r$aaswap <- res$aaswap
    write.table(r,file=paste0(DATADIR,"CombinedPoolAnno.bed"),quote = F,sep = "\t",row.names = F,col.names = F)
    rm(r,res)
  }else{
    cat("run SNPEFF program on 'CombinedPool.vcf' file\n")
    cat("java -Xmx4g -jar ~/snpEff.jar GRCm38.75 CombinedPool.vcf > CombinedPool.anno.vcf\n")
    stop("exiting!\n")
  }
  
}


###############################################################################################
### Calculate read depths at RNA editing sites
###############################################################################################
FILES <- Sys.glob("A:/work/Kreidberg lab/Valerie/RNAEditing_Bowtie/*.sort_1_btPileup.txt")

if(!file.exists(paste0(DATADIR,"CombinedEditsPool_q10.RData"))){
  regions <- read.delim(paste0(DATADIR,"CombinedPoolAnno.bed"),header = F,stringsAsFactors = F)
  regions$id <- paste0(regions$V1,":",regions$V2)
  out <- c()
  for(f in FILES){
    cat(f,"\n")
    z <- read.delim(f,header = T,stringsAsFactors = F)
    z[,5] <- z[,5]+z[,9]
    z[,6] <- z[,6]+z[,10]
    z[,7] <- z[,7]+z[,11]
    z[,8] <- z[,8]+z[,12]
    z$id <- paste0(z$chr,":",z$loc)
    z <- z[,c(1:8,13,14,16)]
    z <- subset(z,z$id%in%regions$id)
    
    z <- merge(z,unique(regions[,c("id","V4")]),by="id")
    z$REF <- as.character(gsub(">\\S*","",z$V4))
    z$ALT <- as.character(gsub("\\S*>","",z$V4))
    z <- z[z$ref==z$REF,]
    z$alt <- apply(z,1,function(x){
      if(x[14]=="A"){
        return(x[6])
      }else if(x[14]=="T"){
        return(x[7])
      }else if(x[14]=="C"){
        return(x[8])
      }else if(x[14]=="G"){
        return(x[9])
      }else{
        return(NA)
      }
    })
    z <- z[,c("chr","loc","id","V4","depth","alt","Insertion","Deletion")] %>% unique
    z$profile=gsub(".sort_1_btPileup.txt","",gsub(".*/","",f))
    z$alt <- as.numeric(z$alt)
    hist(z$alt/z$depth,breaks=100)
    out <- rbind(out,z)
    rm(z,f)
  }
  save(out, file=paste0(DATADIR,"CombinedEditsPool_q10.RData"))
  rm(regions,f,z,out)
}



#########################################################
## hg19
#########################################################
if(F){
  a <- read.delim(paste0(RAWDIR,"darned_hg19.txt"),header = T)
  b <- read.delim(paste0(RAWDIR,"radar_hg19.txt"),header=T) 
  
  a$inrna <- ifelse(a$inrna=="I",as.character("G"),as.character(a$inrna))
  a$edit <- paste0(a$inchr,">",a$inrna)
  b$chromosome <- gsub("chr","",b$chromosome)
  b$tissue <- "."
  
  a <- a[,c(1,2,6,9)]
  b <- b[,c(1:3,12)]
  names(a) <- names(b) <- c("chr","start","gene","tissue")
  gr <- rbind(a,b)
  rm(a,b)
  gr$end <- gr$start
  gr <- as(gr,"GRanges")
  gr <- unique(gr)
  hg19Edits <- gr
  save(hg19Edits,file=paste0(DATADIR,"hg19_RNAEdits_DB.RData"))
}



