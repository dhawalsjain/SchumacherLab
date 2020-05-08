source("functions.R")
source("vars.R")


###############################################################################################
### Editing dynamics over genomic regions
###############################################################################################
if(F){
  load(paste0(DATADIR,"EditingLevelReport_Anno.RData"))
  load(paste0(DATADIR,"DRE_DE_matrix.RData"))
  res <- merge(res,hlpr[,c("idx","snpeff_uq")],by="idx")
  res$col <- ifelse(res$adj.P.Val<0.05,0,1)
  
  cn <- read.delim(paste0(DATADIR,"SnpEffector.txt"),header=T,stringsAsFactors = F)
  cn$Location <- gsub("-\\S*","",cn$Location)
  cn$REF <- ifelse(cn$Allele=="C",as.character("T"),as.character("G"))
  cn$idx <- paste0(cn$Location,":",cn$REF,">",cn$Allele)
  cn <- cn[,c("Consequence","SIFT","idx")] %>% unique
  cn$SIFT <- gsub("\\(\\S*\\)","",cn$SIFT)
  table(cn$SIFT)
  cn$SIFT <- gsub("deleterious_low_confidence|tolerated_low_confidence","tolerated",cn$SIFT)
  cn <- cn[,c("idx","SIFT")] %>% unique
  cn <- cn[cn$SIFT!="-",]
  res <- merge(res,cn,by="idx",all.x=T)
  
  ## volcano plots
  ggplot(res,aes(x=change,y=-log10(adj.P.Val),col=factor(col)))+geom_point(size=2,alpha=0.5)+
    facet_wrap(~comparison,nrow=1,scales = "free")+
    scale_color_brewer(palette = "Dark2",direction = -1)+
    theme_bw()+theme(legend.position = "none")+gg_aes+
    xlab("Editing change")+ylab("FDR, -log10")
  
  ## spreadsheet
  r1 <- res[res$adj.P.Val<0.05,]
  r1 <- r1[,c("idx", "logFC","change", "comparison", "adj.P.Val",
              "known", "repeats","clustered","snpeff_uq",  
              "aaswap", "SIFT","control", "D7","D14", 
              "tracking_id", "gene", "expn_l2fc", "expn_fdr" )] %>% unique
  
  ## editing site associations
  r1 <- res[res$adj.P.Val<0.05,]
  r1$clustered <- ifelse(r1$clustered==0,as.character("loner"), as.character("clustered"))
  qq <- c()
  for(prof in levels(as.factor(r1$comparison))) {
    qwe <- r1[r1$comparison==prof & r1$logFC >0,]
    qq <- rbind(qq,
                data.frame(profile=prof, data.frame(table(qwe$clustered)),total=nrow(qwe),chnge="up") )
    qwe <- r1[r1$comparison==prof & r1$logFC <0,]
    qq <- rbind(qq,
                data.frame(profile=prof, data.frame(table(qwe$clustered)),total=nrow(qwe),chnge="down") )
    rm(qwe)
  }
  names(qq)[2:3] <- c("category","size")
  qq$category <- gsub("splice","splice_region",qq$category)
  qq$category <- gsub("_prime_","'",qq$category)
  qq$category <- gsub("_region","",qq$category)
  qq <- data.table(qq)
  qq <- qq[,totsize:=sum(size),by=list(profile,chnge)]
  qq <- as.data.frame(qq)
  qq$perc <- round(qq$size*100/qq$totsize,2)
  #qq$perc <- ifelse(qq$chnge=="up",qq$perc,qq$perc*(-1))
  
  ggplot(qq[qq$chnge=="up",], aes(y=perc, x=category,fill=category,alpha=0.2)) + geom_bar( stat="identity",position="dodge")+
    scale_fill_brewer(palette = "Dark2")+
    facet_wrap(~profile,nrow = 2)+coord_flip()+
    geom_hline(yintercept = 0,col="black")+
    theme_bw()+gg_aes+
    xlab("")+ylab("% of DRE")
  
  ggplot(qq[qq$chnge=="down",], aes(y=perc, x=category,fill=category,alpha=0.2)) + geom_bar( stat="identity",position="dodge")+
    scale_fill_brewer(palette = "Dark2")+
    facet_wrap(~profile,nrow = 2)+coord_flip()+
    geom_hline(yintercept = 0,col="black")+
    theme_bw()+gg_aes+
    xlab("")+ylab("% of DRE")
  
  ## consequence analyses
  r1 <- res[res$comparison=="D7_control" & !is.na(res$SIFT) & res$adj.P.Val<0.05,]
  r1 <- r1[,c("idx", "logFC","change", "comparison", "adj.P.Val",
              "known", "repeats","clustered","snpeff_uq",  
              "aaswap", "SIFT","control", "D7","D14", 
              "tracking_id", "gene", "expn_l2fc", "expn_fdr" )] %>% unique
  pie(table(r1$SIFT),col=c("red","green"))
  
  
  r1 <- res[res$comparison=="D14_control" & !is.na(res$SIFT) & res$adj.P.Val<0.05,]
  r1 <- r1[,c("idx", "logFC","change", "comparison", "adj.P.Val",
              "known", "repeats","clustered","snpeff_uq",  
              "aaswap", "SIFT","control", "D7","D14", 
              "tracking_id", "gene", "expn_l2fc", "expn_fdr" )] %>% unique
  pie(table(r1$SIFT),col=c("red","green"))
  
  
  ## DRE GO enrichment
  load(paste0(DATADIR,file="DRE_DE_GOEnrichment.RData"))
  z <- out[out$Count>=2,]
  z <- z[order(z$Pvalue,decreasing = F),]
  tms <- unique(z$Term)
  z <- out[out$Term%in%tms,]
  z$Pvalue <- -log10(z$Pvalue)
  z <- z[z$Category!='Cellular Component',]
  ggplot(z,aes(x=comparison,y=Term))+geom_point(aes(col=Pvalue,size=Count*100/Size))+
    xlab("")+ylab("")+theme_bw()+
    scale_color_gradientn( colours = c("blue", "green", "orange", "red"), limits = c( 1,max(z$Pvalue)) )+
    theme(axis.text = element_text(size=14,colour = "black"),
          axis.text.x = element_text(angle = 45,hjust = 1),
          strip.text = element_text(size=16,colour = "black"),
          legend.text = element_text(size=13,colour = "black"),
          legend.title = element_text(size=14,color = "black"))
  
  ###heatmap 
  require(circlize)
  require(ComplexHeatmap)
  ISOFORMCNT="A:/work/Kreidberg lab/Valerie/Bam0_CuffDiff/isoforms.read_group_tracking" ## isoform level diff expression by cuffdiff
  cnt <- read.delim(ISOFORMCNT,header=T,stringsAsFactors = F)
  cnt$condition <- paste0(cnt$condition,"_",cnt$replicate)
  cnt$condition <- gsub("_2","_3",cnt$condition)
  cnt$condition <- gsub("_1","_2",cnt$condition)
  cnt$condition <- gsub("_0","_1",cnt$condition)
  cnt$condition <- gsub("Control","control",cnt$condition)
  cnt <- cnt[,c("tracking_id","condition","FPKM")]
  cnt <- reshape2::dcast(cnt,tracking_id~condition,value.var = "FPKM")
  rownames(cnt) <- cnt$tracking_id
  cnt$tracking_id <- NULL
  cnt <- cnt[rowSums(cnt)>0,]
  cnt <- log2(cnt+1)
  cnt <- apply(cnt,1,function(x){round((x-mean(x))/sd(x),2)}) %>% t %>% as.data.frame
  col_rnorm = colorRamp2(c(min(cnt,na.rm=T), 0, max(cnt,na.rm=T)), c("white", "yellow2", "brown4"))
  res$snpeff_uq <- gsub("splice","splice_region",res$snpeff_uq)
  res$snpeff_uq <- gsub("_prime_","'",res$snpeff_uq)
  
  ht <- list()
  for(prof in c("D7_control","D14_control")){
    for(f in c("3'UTR","intron","splice_region","missense","synonymous")){
      r1 <- subset(res, res$clustered>0 & res$comparison==prof & res$snpeff_uq==f & res$logFC>0 & res$adj.P.Val<0.05)
      r1 <- unique(r1$tracking_id)
      r2 <- subset(res, res$clustered>0 & res$comparison==prof & res$snpeff_uq==f & res$logFC<0 & res$adj.P.Val<0.05)
      r2 <- unique(r2$tracking_id)
      m <- rbind(cnt[rownames(cnt)%in%r1,],cnt[rownames(cnt)%in%r2,])
      rownames(m) <- NULL
      ht[[paste0(prof,"_",f)]] <- 
        Heatmap(m,column_title = f,name = "scale",col = col_rnorm,cluster_columns = F,split = c(rep("Over\nEdited",length(r1)),rep("Under\nEdited",length(r2))))
      rm(r1,r2,m)
    }
  }
  ht[[1]]
  
  
  
}