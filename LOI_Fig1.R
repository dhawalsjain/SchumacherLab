rm(list=ls())
source("functions.R")
source("vars.R")

pllist <- list()

### ADAR expression

gg_aes <- theme(legend.position = "none",
                axis.text = element_text(size=12,color="black"),
                axis.title = element_text(size=13,color="black"),
                legend.text = element_text(size=12,colour = "black"),
                strip.text = element_text(size=13,colour = "black"),
                plot.title = element_text(size=16,colour = "black",hjust = 0.5),
                legend.title = element_blank())



if(T){
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
  
  m <- cnt[cnt$gene=="Adar",]
  pl <- melt(m)
  pl$variable <- gsub("_2","_3",pl$variable)
  pl$variable <- gsub("_1","_2",pl$variable)
  pl$variable <- gsub("_0","_1",pl$variable)
  pl$variable <- gsub("Control","control",pl$variable)
  pl$gt <- gsub("_\\S*","",pl$variable)
  pl$gt <- factor(pl$gt,levels=c("control","D7","D14"))
  
  pllist[["AdarExpn"]] <- local({
    q <- ggplot(pl, aes(x=gt, y=value, fill=gt,alpha=0.2)) + 
    geom_boxplot()+geom_jitter(width = 0.2)+
    scale_fill_manual(values=cols)+
    stat_compare_means(label = "p.signif", method = "t.test",ref.group = "control") +
    theme_bw()+xlab("")+ylab("expression level\n log2 FPKM")+
    theme(axis.text = element_text(size=12,color="black"),
          axis.title = element_text(size=13,color="black"),
          legend.position = "none")
    print(q)
    })
  rm(m,pl,outC,out0,out1,cnt,genes)
  
}

  
## edit level, repeat, known, genomAnno, gneomAnnoClust
if(T){  
  categories <- c("3_prime_UTR", "5_prime_UTR", "intergenic_region", "intron", 
                  "missense", "synonymous", "splice")
  
  load(paste0(DATADIR,"33hc_EditingLevelReport_Anno.RData"))
  hc_edits_anno$condition <- gsub("_\\S*","",hc_edits_anno$sample_name)
  hc_edits_anno$condition <- factor(hc_edits_anno$condition,levels=c("control","D7","D14"))
  
  
  pllist[["EditLevels"]] <- local({
    q <- ggplot(hc_edits_anno, aes(x=condition,y=H,fill=condition,alpha=0.2)) + geom_violin() + theme_bw()+
    geom_boxplot(width=0.1, fill="white",outlier.size = NULL)+
    scale_color_manual(values=cols)+scale_fill_manual(values = cols)+
    stat_compare_means(label = "p.signif", method = "t.test",ref.group = "control")+
    xlab("")+ylab("editing levels\n")+gg_aes
    print(q)
    })
  
  
  pl <- rbind(cbind(hlpr[hlpr$control>0,],gt="control\n(1829)"),
              cbind(hlpr[hlpr$D7>0,],gt="D7\n(725)"),
              cbind(hlpr[hlpr$D14>0,],gt="D14\n(686)"))
  
  pllist[["known"]] <- local({
    q <- ggplot(pl, aes(fill=as.factor(known), x=gt,alpha=0.2)) + geom_bar(position="fill")+
      scale_fill_brewer(palette = "Dark2")+
      theme_bw()+gg_aes+theme(legend.position = "bottom")+
      xlab("")+ylab("fraction of edits")
    print(q)
  })  
  
  ## repeats
  pllist[["repeats"]] <- local({
    q <- ggplot(pl, aes(fill=as.factor(repeats), x=gt,alpha=0.2)) + geom_bar(position="fill")+
    scale_fill_brewer(palette = "Dark2")+
    theme_bw()+gg_aes+theme(legend.position = "bottom")+
    xlab("")+ylab("fraction of edits")
    print(q)
  })

  qq <- c()
  for(f in categories){
    for(prof in levels(pl$gt)) {
      qwe <- pl[pl$gt==prof,]
      qq <- rbind(qq,
                  data.frame(category=f,profile=prof, size=length(unique(grep(f,qwe$snpeff)))))
      rm(qwe)
    }
  }
  qq$category <- gsub("splice","splice_region",qq$category)
  qq$category <- gsub("_prime_","'",qq$category)
  qq$category <- gsub("_region","",qq$category)
  qq <- reshape2::dcast(qq,category~profile,value.var = "size")
  names(qq)[2:4] <- c("Control","D7","D14")
  qq[,2:4] <- apply(qq[,2:4],2,function(x) round(x*100/sum(x),1))
  
  pllist[["anno"]] <- local({
    q <- ggplot() + geom_col(aes(x = 2, y = Control, fill = category), data = qq) +
      geom_text(aes(x=2,label = paste0(Control,"%"), y=100-cumsum(Control)+Control*0.5 ),data = qq)+
      geom_col(aes(x = 3, y = D7, fill = category), data = qq, color = "black") +
      geom_text(aes(x=3,label = paste0(D7,"%"), y=100-cumsum(D7)+D7*0.6 ),data = qq)+
      geom_col(aes(x = 4, y = D14, fill = category), data = qq, color = "black") +
      geom_text(aes(x=4,label = paste0(D14,"%"), y=100-cumsum(D14)+D14*0.5 ),data = qq)+
      xlim(0, 4.5) + labs(x = NULL, y = NULL) + 
      scale_fill_brewer(palette = "Dark2")+
      theme(axis.ticks=element_blank(),
            axis.text=element_blank(),
            axis.title=element_blank(),
            axis.line = element_blank(),
            legend.title = element_blank(),
            plot.title = element_text(size=16,colour = "black",hjust = 0.5),
            legend.position = "bottom")+
      coord_polar(theta = "y") 
    print(q)
  })
  
  

  pl <- data.frame(profile=c("control","D7","D14"),
             sites=c("all","all","all"),
             n=c(sum(hlpr$control>0),sum(hlpr$D7>0), sum(hlpr$D14>0)))
  pl$profile <- factor(pl$profile,levels=c("control","D7","D14"))
  
  pllist[["Numbers"]] <- local({
    q <-ggplot(pl,aes(x=profile,y=n,fill=profile,alpha=0.2))+
      geom_bar(stat = "identity",position="dodge")+
      gg_aes+scale_fill_manual(values=cols)+
      xlab("")+ylab("number of sites")
    print(q)
  })
  
  
}

## GO
if(T){
  load(file=paste0(DATADIR,file="33hc_GOEnrichment_AllEditing.RData"))
  tms <- read.delim(file=paste0(DATADIR,file="SelectedGO_Valerie.txt"),header=F,stringsAsFactors = F)
  tms <- tms$V1
  z <- out[out$Pvalue<0.001,]
  z <- z[order(z$Pvalue,decreasing = F),]
  tms1 <- unique(z$Term)
  tms <- unique(c(tms,tms1))
  z <- out[out$Term%in%tms,]
  z$Pvalue <- -log10(z$Pvalue)
  z <- z[z$Category!='Cellular Component',]
  z$Term <- gsub("protein deubiquitination involved in ubiquitin-dependent protein catabolic process",
                 "protein deubiquitination",z$Term)
  
  pllist[["GO"]] <- local({
    q <- ggplot(z,aes(x=comparison,y=Term))+geom_point(aes(col=Pvalue,size=Count*100/Size))+
      xlab("")+ylab("")+theme_bw()+
      scale_color_gradientn( colours = c("blue", "green", "orange", "red"), limits = c( 1,max(z$Pvalue)) )+
      theme(axis.text = element_text(size=12,colour = "black"),
            axis.text.x = element_text(hjust = 1),
            strip.text = element_text(size=16,colour = "black"),
            legend.text = element_text(size=12,colour = "black"),
            legend.title = element_text(size=12,color = "black"),
            plot.title = element_text(size=16, hjust=0.5,colour = "black"))
    print(q)
  })
  
  pdf(paste0(FIGDIR,"AllEdit_SelectedGOTerms.pdf"),width = 8,height = 6)
  pllist[["GO"]]
  dev.off()
  
  
}


lay <- rbind(c(1,2,3),
             c(4,5,3),
             c(6,7,3))

pl <- list(pllist[[1]],pllist[[6]],pllist[[7]],pllist[[2]],pllist[[3]],pllist[[4]],pllist[[5]])
names(pl)

pdf(paste0(FIGDIR,"LOI_Fig1.pdf"),width = 12,height = 8)
grid.arrange(grobs = pl,widths=c(6,6,25),  layout_matrix = lay)
dev.off()



## correlation plot editing level vs expression
if(F){
  load(paste0(DATADIR,"33_Spreadsheet_HCEditsWithCompleteAnno.RData"))
  
  pdf(paste0(FIGDIR,"Test_ControlExpn vs EditingLevels.pdf"),width = 5,height = 5)
  ggplot(SPN,aes(x=Control_edit,y=Control_log2Expn))+geom_point()+
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0))+gg_aes
  dev.off()
  
  
}


## All edits
if(F){
  load(file=paste0(DATADIR,file="33hc_GOEnrichment_AllEditing.RData"))
  tms <- read.delim(file=paste0(DATADIR,file="SelectedGO_Valerie.txt"),header=F,stringsAsFactors = F)
  tms <- tms$V1
  z <- out[out$Pvalue<0.001,]
  z <- z[order(z$Pvalue,decreasing = F),]
  tms1 <- unique(z$Term)
  #tms <- unique(c(tms,tms1))
  z <- out[out$Term%in%tms,]
  #z$Pvalue <- -log10(z$Pvalue)
  z$Term <- gsub("protein deubiquitination involved in ubiquitin-dependent protein catabolic process",
                 "protein deubiquitination",z$Term)
  
    q <- ggplot(z,aes(x=comparison,y=Term))+geom_point(aes(col=-log10(Pvalue),size=Count*100/Size))+
      xlab("")+ylab("")+theme_bw()+
      scale_color_gradientn( colours = c("blue", "green", "orange", "red"), limits = c( 1,max(-log10(z$Pvalue)) ))+
      theme(axis.text = element_text(size=12,colour = "black"),
            axis.text.x = element_text(hjust = 1),
            strip.text = element_text(size=16,colour = "black"),
            legend.text = element_text(size=12,colour = "black"),
            legend.title = element_text(size=12,color = "black"),
            plot.title = element_text(size=16, hjust=0.5,colour = "black"))
  
  pdf(paste0(FIGDIR,"AllEdit_SelectedGOTerms.pdf"),width = 8,height = 6)
  q
  dev.off()
  
  
}


## 3'UTR
if(T){
  load(paste0(DATADIR,file="33hc_DRE_DE_GOEnrichment_byGenomicfeatures.RData"))
  out4$comparison <- paste0(out4$comparison,"_",out4$change)
  out4$comparison <- gsub("_control","",out4$comparison)
  out4$comparison <- gsub("_up","_overedited",out4$comparison)
  out4$comparison <- gsub("_dn","_underedited",out4$comparison)
  tms <- read.delim(file=paste0(DATADIR,file="SelectedGO_Valerie_3UTR.txt"),header=F,stringsAsFactors = F)
  tms <- tms$V1
  z <- out4[out4$Pvalue<0.0001,]
  z <- z[order(z$Pvalue,decreasing = F),]
  tms1 <- unique(z$Term)
  #tms <- unique(c(tms,tms1))
  z <- out4[out4$Term%in%tms,]
  #z$Pvalue <- -log10(z$Pvalue)
  z <- z[z$genomeAnno=="3_prime_UTR",]
  z$Term <- gsub("protein deubiquitination involved in ubiquitin-dependent protein catabolic process",
                 "protein deubiquitination",z$Term)
  z$comparison <-factor(z$comparison,levels=c("D7_overedited", "D7_underedited", "D14_overedited", "D14_underedited"))
  
 q <- ggplot(z,aes(x=comparison,y=Term))+geom_point(aes(col=-log10(Pvalue),size=Count*100/Size))+
      xlab("")+ylab("")+theme_bw()+
      scale_color_gradientn( colours = c("blue", "green", "orange", "red"), limits = c( 1,max( -log10(z$Pvalue)) ))+
      theme(axis.text = element_text(size=12,colour = "black"),
            axis.text.x = element_text(hjust = 1,angle = 45),
            strip.text = element_text(size=16,colour = "black"),
            legend.text = element_text(size=12,colour = "black"),
            legend.title = element_text(size=12,color = "black"),
            plot.title = element_text(size=16, hjust=0.5,colour = "black"))
  
  pdf(paste0(FIGDIR,"3UTREdit_SelectedGOTerms.pdf"),width = 8,height = 7)
  q
  dev.off()
  
  
}


## intron
if(T){
  load(paste0(DATADIR,file="33hc_DRE_DE_GOEnrichment_byGenomicfeatures.RData"))
  out4$comparison <- paste0(out4$comparison,"_",out4$change)
  out4$comparison <- gsub("_control","",out4$comparison)
  out4$comparison <- gsub("_up","_overedited",out4$comparison)
  out4$comparison <- gsub("_dn","_underedited",out4$comparison)
  tms <- read.delim(file=paste0(DATADIR,file="SelectedGO_Valerie_Intron.txt"),header=F,stringsAsFactors = F)
  tms <- tms$V1
  z <- out4[out4$Pvalue<0.1 & out4$genomeAnno=="intron" & out4$comparison=="D7_overedited",]
  z <- z[order(z$Pvalue,decreasing = F),]
  tms1 <- unique(z$Term)
  #tms <- unique(c(tms,tms1))
  z <- out4[out4$Term%in%tms,]
  #z$Pvalue <- -log10(z$Pvalue)
  z <- z[z$genomeAnno=="intron",]
  z$Term <- gsub("protein deubiquitination involved in ubiquitin-dependent protein catabolic process",
                 "protein deubiquitination",z$Term)
  z$comparison <-factor(z$comparison,levels=c("D7_overedited", "D7_underedited", "D14_overedited", "D14_underedited"))
  
  q <- ggplot(z,aes(x=comparison,y=Term))+geom_point(aes(col=-log10(Pvalue),size=Count*100/Size))+
    xlab("")+ylab("")+theme_bw()+
    scale_color_gradientn( colours = c("blue", "green", "orange", "red"), limits = c( 1,max( -log10(z$Pvalue)) ))+
    theme(axis.text = element_text(size=12,colour = "black"),
          axis.text.x = element_text(hjust = 1,angle = 45),
          strip.text = element_text(size=16,colour = "black"),
          legend.text = element_text(size=12,colour = "black"),
          legend.title = element_text(size=12,color = "black"),
          plot.title = element_text(size=16, hjust=0.5,colour = "black"))
  
  pdf(paste0(FIGDIR,"IntronEdit_SelectedGOTerms.pdf"),width = 8,height = 6)
  q
  dev.off()
  
  
}


