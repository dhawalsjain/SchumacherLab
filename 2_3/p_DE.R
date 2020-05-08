rm(list=ls())
source("functions.R")
source("vars.R")


###############################################################################################
### QC
###############################################################################################
if(F){
  
  load(paste0(DATADIR,"DEExpression_DESeq2.RData"))
  d7 <- out0[[1]]
  d7 <- d7[!is.na(d7$padj) & d7$padj<0.05,]
  res <- out0[[1]]
  y <- c(nrow(res[!is.na(res$padj) & res$padj<FDR_CUTOFF & res[,3]<0,]),
         nrow(res[!is.na(res$padj) & res$padj<FDR_CUTOFF & res[,3]>0,]))
  #y <- c(nrow(res)-sum(y), y)
  names(y) <- c("up","down")
  h1 <- hist(res$pvalue[is.na(res$padj)], breaks=0:50/50, plot=FALSE)
  h2 <- hist(res$pvalue[!is.na(res$padj)], breaks=0:50/50, plot=FALSE)
  colori <- c(`do not pass`="khaki", `pass`="powderblue")
  
  
  par(mfrow=c(2,2))
  DESeq2::plotMA(out0[[2]],main=paste0("MA plot\n (QC plot)\n"))
  DESeq2::plotDispEsts(out0[[2]],main="Dispersion estimate\n (QC plot)")
  
  barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
          col = colori, space = 0, ylab="frequency",main="P-value distribution\n (QC plot)")
  text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
       adj = c(0.5,1.7), xpd=NA)
  legend("topright", fill=rev(colori), legend=rev(names(colori)))
  
  q <-barplot(y,main = paste0("Number of differential genes\n",gsub("log2","",names(res)[2])),ylab = "",cex.names  = 1.25,
              cex.axis = 1.25, ylim=c(0,max(y)+0.1*max(y)), col=c("lightblue","pink"))
  text(x = q, y = y, label = y, pos = 3,adj = -1, cex = 1.25, col = "red")
  
  
  
  d7$baseMean <- d7$lfcSE <- d7$stat <- NULL
  names(d7)[c(1,4)] <- c("ensemble_id","FDR")
  
}

### overlap analysis
if(F){
  load(paste0(DATADIR,"DEExpression_DESeq2.RData"))
  d7 <- out0[[1]]
  d7 <- d7[!is.na(d7$padj) & d7$padj<0.05,]
  d14 <- out1[[1]]
  d14 <- d14[!is.na(d14$padj) & d14$padj<0.05,]
  
    m <- cbind(comparison="D7",d7[,c(1,10,3)])
    names(m)[4] <- "log2FC"
    n <- cbind(comparison="D14",d14[,c(1,10,3)])
    names(n)[4] <- "log2FC"
    m <- rbind(m,n)
    rm(n)
    m$status <- ifelse(m$log2FC>0,as.character("down"),as.character("up"))
    m$log2FC=1
    m$comparison <- paste0(m$comparison,"_",m$status)
    m <- reshape2::dcast(m,id+gene~comparison,value.var = "log2FC")
    m[is.na(m)] <- 0
    
    upset(m, sets = names(m)[3:6], mb.ratio = c(0.5, 0.5), order.by = c("degree"),
          text.scale=rep(1.5,6),point.size=5,matrix.color = "steelblue4")
  
  
  
}

###############################################################################################
### ADAR gene FPKM
###############################################################################################
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
  
  go <- read.delim("A:/work/mm10_Annotations/GO_term_summary_20200123_210425.txt",row.names = NULL,stringsAsFactors = F)
  pods <- read.csv("A:/work/mm10_Annotations/scrna.mouse.podocyte.refined.csv",stringsAsFactors = F)
  pods <-pods$Podocytes.gene
  
  m <- cnt[cnt$gene=="Adar",]
  pl <- melt(m)
  pl$variable <- gsub("_2","_3",pl$variable)
  pl$variable <- gsub("_1","_2",pl$variable)
  pl$variable <- gsub("_0","_1",pl$variable)
  pl$variable <- gsub("Control","control",pl$variable)
  pl$gt <- gsub("_\\S*","",pl$variable)
  
  ggplot(pl, aes(x=gt, y=value, fill=gt,alpha=0.2)) + 
    geom_boxplot()+geom_jitter(width = 0.2)+
    scale_fill_manual(values=cols)+
    stat_compare_means(label = "p.signif", method = "t.test",ref.group = "control") +
    theme_bw()+xlab("")+ylab("editing ratio")+
    theme(axis.text = element_text(size=18,color="black"),
          axis.title = element_text(size=18,color="black"),
          legend.position = "none")
  rm(m,pl)
  
  
  m <- cnt[cnt$gene%in%pods,1:9]
  m <- m[rowSums(m)>0,]
  m <- t(apply(m,1,function(x) round((x-mean(x))/sd(x),3)))
  m <- round(m,3)
  x <- rownames(m)
  x <- cnt[match(x,rownames(cnt)),]$gene
  #x[-seq(1,length(x),2)] <- NA
  hmcols<-(colorRampPalette(c("yellow","orange","red4"))(256))
  heatmap3((m),tck=0.8,Colv = NA,cexRow = 0.7,labRow = x,
           symm = F,margins = c(5,20),cexCol = 0.7,col = hmcols,scale="none",useRaster = T)
  
  
  
  mygenes <- unique(go$MGI.Gene.Marker.ID)
  m <- cnt[cnt$gene%in%mygenes,1:9]
  m <- m[rowSums(m)>0,]
  m <- t(apply(m,1,function(x) round((x-mean(x))/sd(x),3)))
  m <- round(m,3)
  x <- rownames(m)
  x <- cnt[match(x,rownames(cnt)),]$gene
  x[-seq(1,length(x),5)] <- NA
  hmcols<-(colorRampPalette(c("yellow","orange","red4"))(256))
  heatmap3((m),tck=0.8,Colv = NA,cexRow = 0.7,labRow = x,
           symm = F,margins = c(5,20),cexCol = 0.7,col = hmcols,scale="none",useRaster = T)
  
  
}



###############################################################################################
### Gene Ontology
###############################################################################################
### DE genes
if(F){
  load(paste0(DATADIR,"DEgenes_GOEnrichment_DESeq2.RData"))
  
  z <- pl
  #z <- z[z$Pvalue<0.001,]
  z <- z[z$Count>=5,]
  z <- z[z$Category=='Biological Process',]
  z <- z[order(z$Pvalue,decreasing = F),]
  z <- z[1:25,]
  tms <- unique(z$Term)
  z <- pl[pl$Term%in%tms,]
  z$Pvalue <- -log10(z$Pvalue)
  
  ggplot(z,aes(x=comparison,y=Term))+geom_point(aes(col=Pvalue,size=Count*100/Size))+
    xlab("")+ylab("")+theme_bw()+
    scale_color_gradientn( colours = c("blue", "green", "orange", "red"), limits = c( 1,max(z$Pvalue)) )+
    theme(axis.text = element_text(size=14,colour = "black"),
          axis.text.x = element_text(angle = 45,hjust = 1),
          strip.text = element_text(size=16,colour = "black"),
          legend.text = element_text(size=13,colour = "black"),
          legend.title = element_text(size=14,color = "black"))
}

## DE (batch effect correction)
if(F){
  load(paste0(DATADIR,"WT1_ValerieCombibedDEExpn_DESeq2_batch.RData"))
  load(paste0(DATADIR,"WT1_ValerieCombibedDEExpn_DESeq2_SVA.RData"))
  load(paste0(DATADIR,"DEExpression_DESeq2.RData"))
  load(paste0(DATADIR,"DEExpression_DESeq2_WT1.RData"))
  expn <- read.delim("A:/work/Kreidberg lab/Valerie/Bam0_CuffDiff/isoform_exp.diff",header=T,stringsAsFactors = F)
  expn <- unique(expn[,2:3])
  
  names(cf_sva) <- c("D7_R1", "D7_R2", "D7_R3", "Control_R3", "Control_R2", "D14_R1", 
                     "D14_R2", "D14_R3", "Control_R1", "D9_S1", "D14_S1", "D14_S2", 
                     "D9_S2", "Control_S1", "D9_S3", "D14_S3", "Control_S2", 
                     "Control_S3")
  pca_data=prcomp(t(cf_sva))
  pca_data_perc=round(100*pca_data$sdev^2/sum(pca_data$sdev^2),1)
  df_pca_data=data.frame(PC1 = pca_data$x[,1], PC2 = pca_data$x[,2], 
                         sample = names(cf_sva), condition=gsub("_\\S+","",names(cf_sva)),
                         label=names(cf_sva),
                         batch=c(rep("New",9), rep("Old",9)))
  ggplot(df_pca_data, aes(PC1,PC2, color = batch,label=label))+geom_text_repel(box.padding = 0.1)+
    geom_point(size=8)+
    labs(x=paste0("PC1 (",pca_data_perc[1],")"), y=paste0("PC2 (",pca_data_perc[2],")"))+
    scale_color_brewer(palette = "Dark2")+
    theme(legend.position = "bottom")
  
  
  
  names(out0[[1]])[2] <- names(out1[[1]])[2] <- names(out2[[1]])[2] <- names(out3[[1]])[2] <- "log2FC"
  z <- c()
  z <- rbind(cbind(expt="new",test="D7_control",out0[[1]]),
             cbind(expt="new",test="D14_control",out1[[1]]),
             cbind(expt="old",test="D9_control",out2[[1]]),
             cbind(expt="old",test="D14_control",out3[[1]]))
  rm(out0,out1,out2,out3)
  
  z <- z[!is.na(z$padj) & z$padj<0.05,]
  z$test <- gsub("_control","",z$test)
  z$log2FCc <- ifelse(z$log2FC>0,as.character("up"),as.character("down"))
  z$var <- paste0(z$test,"_",z$expt,"_",z$log2FCc)
  z$test <- 1
  yy <- reshape2::dcast(z,id~var,value.var = "log2FC")
  yy <- merge(yy,expn,by.x="id",by.y="gene_id")
  yy[is.na(yy)] <- 0
  yy <- yy[,c(1,10,2:9)]
  z <- reshape2::dcast(z,id~var,value.var = "test")
  z[is.na(z)] <- 0
  
  
  upset(z, sets = names(z)[2:5], mb.ratio = c(0.5, 0.5), order.by = c("degree"),
        text.scale=rep(1.5,6),point.size=5,matrix.color = "steelblue4")
  
  upset(z, sets = names(z)[6:9], mb.ratio = c(0.5, 0.5), order.by = c("degree"),
        text.scale=rep(1.5,6),point.size=5,matrix.color = "steelblue4")
  
  upset(z, sets = names(z)[c(2,3,8,9)], mb.ratio = c(0.5, 0.5), order.by = c("degree"),
        text.scale=rep(1.5,6),point.size=5,matrix.color = "steelblue4")
  
  upset(z, sets = names(z)[c(4:7)], mb.ratio = c(0.5, 0.5), order.by = c("degree"),
        text.scale=rep(1.5,6),point.size=5,matrix.color = "steelblue4")
  
  
  z <- DESeq2::counts(out[[2]],normalized=T) %>% as.data.frame
  z <- z[rowSums(z)>0,]
  z <- log2(z+1)
  expn <- read.delim("A:/work/Kreidberg lab/Valerie/Bam0_CuffDiff/isoform_exp.diff",header=T,stringsAsFactors = F)
  expn <- unique(expn[,2:4])
  expn$chr <- gsub(":\\S*","",expn$locus) %>% as.character
  z$gene_id <- rownames(z)
  z <- merge(z,expn[,c("chr","gene_id")],by="gene_id")
  z$gene_id <- NULL
  z <- data.table(z)
  z <- z[, lapply(.SD, median, na.rm=TRUE), by=chr ] %>% as.data.frame
  names(z) <- c("chr", "D7_1", "D7_2", "D7_3", "Control_3", "Control_2", 
                "D14_1", "D14_2", "D14_3", "Control_1", "D9_1_wt1", "D14_1_wt1", 
                "D14_2_wt1", "D7_2_wt1", "Control_1_wt1", "D7_3_wt1", "D14_3_wt1", "Control_2_wt1", 
                "Control_3_wt1")
  z <- melt(z)
  z <- z[z$chr%in%c(1:19,"X"),]
  ggplot(z[z$chr%in%c(1:6),],aes(x=variable,y=value,fill=variable,alpha=0.2))+geom_bar(stat="identity",position="dodge")+
    facet_wrap(~chr,ncol=2,scales = "free_x")+
    coord_flip()+gg_aes
  
  
  
}


## DREs
if(F){
  load(paste0(DATADIR,file="DRE_DE_GOEnrichment.RData"))
  z <- out[out$Count>=2 & out$change=="all",]
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
  
  
  z <- out1[out1$Count>=2 & out1$change=="all",]
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
  
}


## DREs (3UTR or introns)
if(F){
  load(paste0(DATADIR,file="DRE_DE_GOEnrichment_3UTR.RData"))
  z <- out2[out2$Count>=3 & out2$change=="all",]
  z <- z[order(z$Pvalue,decreasing = F),]
  tms <- unique(z$Term)
  z <- out[out$Term%in%tms,]
  z$Pvalue <- -log10(z$Pvalue)
  z <- z[z$Category!='Cellular Component',]
  
  ggplot(z,aes(x=comparison,y=Term))+geom_point(aes(col=Pvalue,size=Count*100/Size))+
    xlab("")+ylab("")+theme_bw()+
    scale_color_gradientn( colours = c("blue", "green", "orange", "red"), limits = c( 1,max(z$Pvalue)) )+
    ggtitle("Editing changes within 3'UTR")+
    theme(axis.text = element_text(size=14,colour = "black"),
          axis.text.x = element_text(angle = 45,hjust = 1),
          strip.text = element_text(size=16,colour = "black"),
          legend.text = element_text(size=13,colour = "black"),
          legend.title = element_text(size=14,color = "black"),
          plot.title = element_text(size=18, hjust=0.5,colour = "black"))
  

  
  load(paste0(DATADIR,file="DRE_DE_GOEnrichment_intron.RData"))
  z <- out3[out3$Count>=7 & out3$change=="all",]
  z <- z[order(z$Pvalue,decreasing = F),]
  tms <- unique(z$Term)
  z <- out3[out3$Term%in%tms,]
  z$Pvalue <- -log10(z$Pvalue)
  z <- z[z$Category!='Cellular Component',]
  
  ggplot(z,aes(x=comparison,y=Term))+geom_point(aes(col=Pvalue,size=Count*100/Size))+
    xlab("")+ylab("")+theme_bw()+
    scale_color_gradientn( colours = c("blue", "green", "orange", "red"), limits = c( 1,max(z$Pvalue)) )+
    ggtitle("Editing changes within Introns")+
    theme(axis.text = element_text(size=14,colour = "black"),
          axis.text.x = element_text(angle = 45,hjust = 1),
          strip.text = element_text(size=16,colour = "black"),
          legend.text = element_text(size=13,colour = "black"),
          legend.title = element_text(size=14,color = "black"),
          plot.title = element_text(size=18, hjust=0.5,colour = "black"))
  
  
  
  
  load(paste0(DATADIR,file="DRE_DE_GOEnrichment_byGenomicfeatures.RData"))
  z <- out4[out4$Count>=7,]
  z <- z[order(z$Pvalue,decreasing = F),]
  tms <- unique(z$Term)
  z <- out4[out4$Term%in%tms,]
  z$Pvalue <- -log10(z$Pvalue)
  z <- z[z$Category!='Cellular Component',]
  
  ggplot(z,aes(x=comparison,y=Term))+geom_point(aes(col=Pvalue,size=Count*100/Size))+
    xlab("")+ylab("")+theme_bw()+
    scale_color_gradientn( colours = c("blue", "green", "orange", "red"), limits = c( 1,max(z$Pvalue)) )+
    ggtitle("Editing changes within Introns")+
    facet_wrap(~genomeAnno,nrow=1)+
    theme(axis.text = element_text(size=14,colour = "black"),
          axis.text.x = element_text(angle = 45,hjust = 1),
          strip.text = element_text(size=16,colour = "black"),
          legend.text = element_text(size=13,colour = "black"),
          legend.title = element_text(size=14,color = "black"),
          plot.title = element_text(size=18, hjust=0.5,colour = "black"))
  
  
}


