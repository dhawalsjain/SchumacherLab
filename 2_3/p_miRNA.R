rm(list=ls())
source("functions.R")
source("vars.R")

###############################################################################################
### miREnrich by seedOverlaps
###############################################################################################
if(F){
  load(paste0(DATADIR,"EditingLevelReport_Anno.RData"))
  load(paste0(DATADIR,"DRE_DE_matrix.RData"))
  load(paste0(DATADIR,"miRSeedOverlaps_filtermiRNAs.RData"))
  
  
  myplotdf <- function(r1,cutoff=3,titl=""){
    r2 <- table(r1$mir,r1$seq) %>% as.data.frame
    r2 <- reshape2::dcast(r2,Var1~Var2,value.var = "Freq")
    r2$a <- r2$ref
    r2$b <- nrow(r1) - r2$a
    r2$c <- r2$alt
    r2$d <- nrow(r1) - r2$c
    r2$pval <- apply(r2[,4:7],1,function(x){
      fisher.test(matrix(x,nrow = 2,byrow = T),alternative = "two.sided")$p.value
    })
    head(r2)
    r2$pval <- p.adjust(r2$pval,method = "fdr")
    r2$pval <- -log10(r2$pval)
    r2$pval <- ifelse(r2$a>r2$c,r2$pval, r2$pval*(-1))
    r2 <- r2[order(r2$pval,decreasing = T),]
    r2$x <- 1:nrow(r2)
    r2$col <- 0
    r2$col <- ifelse(r2$pval>cutoff,1,r2$col)
    r2$col <- ifelse(r2$pval< (-1)*cutoff,2,r2$col)
    r2$label <- NA
    r2$label <- ifelse(r2$col>0,as.character(r2$Var1),NA)
    gg <- ggplot(r2,aes(x=x,y=pval,col=col))+
      geom_rect(aes(xmin = -Inf,xmax = Inf,ymin = 0, ymax = Inf, fill = "lost"), alpha = .005,fill="yellow")+
      geom_rect(aes(xmin = -Inf,xmax = Inf,ymin = 0, ymax = -Inf, fill = "gain"), alpha = .01,fill="pink")+
      geom_point(size=2)+
      geom_label_repel(aes(label = label),box.padding   = 0.5,point.padding = 0.5,segment.color = 'grey50')+ 
      ylab("False discovery rate\n -log10")+
      ggtitle(titl)+
      theme_bw()+
      theme(axis.title = element_text(size=16,colour = "black"),
            axis.text = element_text(size=14,colour = "black"),
            plot.title = element_text(size=20,colour = "black",hjust = 0.5),
            legend.position = "none",
            axis.title.x = element_blank(),
            axis.text.x = element_blank())
    gg
  }
  
  r1 <- mirSeedOlp[mirSeedOlp$idx%in%hlpr[hlpr$control>0,]$idx,]
  p1 <- myplotdf(r1,titl = "Control")
  r1 <- mirSeedOlp[mirSeedOlp$idx%in%hlpr[hlpr$D7>0,]$idx,]
  p2 <- myplotdf(r1,titl = "Day-7")
  r1 <- mirSeedOlp[mirSeedOlp$idx%in%hlpr[hlpr$D14>0,]$idx,]
  p3 <- myplotdf(r1,titl = "Day-14")
  
  pdf(paste0(FIGDIR,"miRSeedBasedEnrichments.pdf"),width = 8,height = 10)
  grid.arrange(p1,p2,p3,nrow=3)
  dev.off()
  
  r1 <- res[res$adj.P.Val<0.05 & res$comparison=="D7_control",]
  r1 <- mirSeedOlp[mirSeedOlp$idx%in%r1$idx,]
  p1 <- myplotdf(r1,cutoff = 1.3,titl = "D7 vs control")
  r1 <- res[res$adj.P.Val<0.05 & res$comparison=="D14_control",]
  r1 <- mirSeedOlp[mirSeedOlp$idx%in%r1$idx,]
  p2 <- myplotdf(r1,cutoff = 1.3,titl = "D14 vs control")
  
  pdf(paste0(FIGDIR,"miRSeedBasedEnrichments_DRE.pdf"),width = 8,height = 10)
  grid.arrange(p1,p2,nrow=3)
  dev.off()
  
}


###############################################################################################
### miRNAvalidated  targets my miRTarscan
###############################################################################################
if(F){
  load(paste0(DATADIR,"DEExpression_DESeq2.RData"))
  ggg <- out0[[1]][,c("id","geneStart","geneEnd","gene")]
  ggg$length <- abs(ggg$geneStart-ggg$geneEnd)
  ggg <- ggg[,c("id","length","gene")] %>% unique
  cnt <- DESeq2::counts(outC[[2]],normalized=T) %>% as.data.frame
  cnt$id <- rownames(cnt)
  cnt <- merge(cnt,ggg,by="id",all.x=T)
  rownames(cnt) <- cnt$id
  cnt$id <- NULL
  cnt[,1:9] <- cnt[,1:9]*1000/cnt$length
  cnt$length <- NULL
  cnt[,1:9] <- apply(cnt[,1:9],2,function(x){round(x*1e6/sum(x),2)})
  cnt <- as.data.frame(cnt)
  cnt[,1:9] <- round(log2(cnt[,1:9]+1),2)
  names(cnt)[1:9] <- unname(sampl[names(cnt)[1:9]])
  cnt <- as.data.frame(cnt)
  rm(out0,out1,outC,ggg)
  
  load(paste0(DATADIR,"EditingLevelReport_Anno.RData"))
  load(paste0(DATADIR,"DRE_DE_matrix.RData"))
  load(paste0(DATADIR,"miRSeedOverlaps_filtermiRNAs.RData"))
  load(paste0(DATADIR,"genes.RData"))
  
  mirSeedOlp$txID <- gsub("_\\S*","",mirSeedOlp$names)
  mirSeedOlp <- merge(mirSeedOlp, unique(gGenes[,c("txID","gene")]),by="txID",all.x=T )
  mirSeedOlp$mir <- gsub("miR","mmu-miR",mirSeedOlp$mir)
  mirSeedOlp$test <- paste0(mirSeedOlp$gene,"|",mirSeedOlp$mir)
  
  
  tars <- read.delim("A:/work/mm10_Annotations/miRTarBase.txt",header=T,stringsAsFactors = F)
  tars$test <- paste0(tars$Target.Gene,"|",tars$miRNA)
  
  tars <- tars[tars$test%in%mirSeedOlp$test,]
  tars <- merge(tars,unique(mirSeedOlp[,c("idx","seq","test")]),by="test",all.x=T)
  ValTars <- reshape2::dcast(tars,test+Experiments+Support.Type+idx+Target.Gene~seq)
  
  ValTars <- merge(ValTars,cnt,by.x='Target.Gene',by.y="gene",all.x=T)
  
  z <- res[,c("idx","D7","D14","control")] %>% unique
  names(z)[2:4] <- paste0("editing_",names(z)[2:4])
  ValTars <- merge(ValTars,z,by="idx",all.x=T)
  ValTars[,8:16] <- apply(ValTars[,8:16],1,function(x) round((x-mean(x))/sd(x),2)) %>% t %>% as.data.frame
  
  ### plot
  m <- ValTars
  m$alt <- ifelse(m$alt>1,1,m$alt)
  m$ref <- ifelse(m$ref>1,1,m$ref)
  m <- m[m$ref!=m$alt,]
  rownames(m) <- 1:nrow(m)
  
  m <- m[order(m$alt,decreasing = T),]
  partition = ifelse(m$alt>0,as.character("gain"),as.character("loss"))
  col1 = colorRamp2(c(min(m[,8:16],na.rm=T), 0, max(m[,8:16],na.rm=T)), c("#fee8c8", "#fdbb84", "#e34a33"))
  col2 = colorRamp2(c(min(m[,17:19],na.rm=T), 0, max(m[,17:19],na.rm=T)), c("#efedf5", "#bcbddc", "#756bb1"))
  rownames(m) <- NULL
  
  
  ht <- Heatmap(partition, col = structure(2:7, names = paste0(levels(as.factor(partition)) )), name = "miRNA-binding",
          show_row_names = FALSE, width = unit(1, "mm"))+
  Heatmap(m[,8:16],column_title = "Expression",name = "Exp",col = col1,cluster_columns = F)+
  Heatmap(m[,17:19],column_title = "Editing",name = "Edit",col = col2,cluster_columns = F)
  
  
  draw(ht, split = partition, heatmap_legend_side = "bottom", gap = unit(2, "mm"))
  
  
  m$alt <- ifelse(m$alt>0, as.character("miRNA-binding GAIN"), as.character("miRNA-binding LOSS"))
  m$ref <- NULL
  names(m)[6] <- "miRNABinding"
  
  #datatable(m, extensions = 'Buttons',options = list(searching = TRUE,pageLength = 10,dom = 'Bfrtip', 
  #                                                   buttons = c('excel', "csv")))
  
}


###############################################################################################
###  Expression values for miRNA validated tagets, for which miRNA is edited 
###############################################################################################
if(F){
  load(paste0(DATADIR,"DEExpression_DESeq2.RData"))
  ggg <- out0[[1]][,c("id","geneStart","geneEnd","gene")]
  ggg$length <- abs(ggg$geneStart-ggg$geneEnd)
  ggg <- ggg[,c("id","length","gene")] %>% unique
  cnt <- DESeq2::counts(outC[[2]],normalized=T) %>% as.data.frame
  cnt$id <- rownames(cnt)
  cnt <- merge(cnt,ggg,by="id",all.x=T)
  rownames(cnt) <- cnt$id
  cnt$id <- NULL
  cnt[,1:9] <- cnt[,1:9]*1000/cnt$length
  cnt$length <- NULL
  cnt[,1:9] <- apply(cnt[,1:9],2,function(x){round(x*1e6/sum(x),2)})
  cnt <- as.data.frame(cnt)
  cnt[,1:9] <- round(log2(cnt[,1:9]+1),2)
  names(cnt)[1:9] <- unname(sampl[names(cnt)[1:9]])
  cnt <- as.data.frame(cnt)
  rm(out0,out1,outC,ggg)
  
  tars <- read.delim("A:/work/mm10_Annotations/miRTarBase.txt",header=T,stringsAsFactors = F)
  
  
  load(paste0(DATADIR,"EditedmiRNAs.RData"))
  load(paste0(DATADIR,"DRE_DE_matrix.RData"))
  tars <- tars[tars$miRNA%in%editedmiRs$name,]
  M <- merge(tars,cnt,by.x='Target.Gene',by.y="gene",all.x=T)
  M[,10:18] <- apply(M[,10:18],1,function(x) round((x-mean(x,na.rm=T))/sd(x,na.rm=T),2)) %>% t %>% as.data.frame
  M <- M[!is.na(M$Control_R1),]
  
  editedmiRs <- editedmiRs[,c("name","D7","D14","control")] %>% unique
  editedmiRs$D7 <- ifelse(is.na(editedmiRs$D7),0,1)
  editedmiRs$D14 <- ifelse(is.na(editedmiRs$D14),0,1)
  editedmiRs$control <- ifelse(is.na(editedmiRs$control),0,1)
  editedmiRs$name <- gsub("mir","miR",editedmiRs$name)
  
  M <- merge(M,editedmiRs,by.x="miRNA",by.y="name",all.x=T)
  rownames(M) <- paste0(M$miRNA,"|",M$Target.Gene,"|",1:nrow(M))
  
  M <- unique(M[,10:21])
  
  col1 = colorRamp2(c(min(M[,1:9],na.rm=T), 0, max(M[,1:9],na.rm=T)), c("#fee8c8", "#fdbb84", "#e34a33"))
  col2 = colorRamp2(c(0,1), c("#f0f0f0", "#636363"))
  
  Heatmap(M[,1:9],column_title = "Expression",name = "Exp",col = col1,cluster_columns = F,show_row_names = F)+
    Heatmap(M[,10:12],column_title = "Edited In Sample",name = "Edit",col = col2,cluster_columns = F)
  
  

}




###########################
#### MiRNA binding gain and loss summary across condiutions
###########################
if(F){
  load(paste0(DATADIR,"miRbinding_GainLoss_Timecourse.RData"))
  load(paste0(DATADIR,"miR_lossGain_dataframe.RData"))
  
  
  miRgain <- miRgain[order(miRgain$sd,decreasing = T),]
  miRloss <- miRloss[order(miRloss$sd,decreasing = T),]
  
  z <- miRgain[1:15,]
  z$sd <- NULL
  z <- melt(z,measure.vars = names(z)[2:4])
  z$status <- "gain in miRNA binding"
  z1 <- miRloss[1:15,]
  z1$sd <- NULL
  z1 <- melt(z1,measure.vars = names(z1)[2:4])
  z1$status <- "loss in miRNA binding"
  z <- rbind(z,z1)
  rm(z1)
  z$variable <- gsub("_\\S*","",z$variable)
  
  ggplot(z,aes(x=variable,y=miR))+#geom_point(aes(col=status,size=value,alpha=0.4))+
    xlab("")+ylab("")+theme_bw()+facet_wrap(~status,ncol = 1,scales = "free_y")+
    geom_text(aes(label=value,size=value,col=status),hjust=0, vjust=0)+
    theme(axis.text = element_text(size=14,colour = "black"),
          strip.text = element_text(size=16,colour = "black"),
          legend.position = "none")
  rm(z)
  
}



### miRNA-mRNA correlation plots based on Seedoverlaps using our data
if(F){
  
  load(paste0(DATADIR,"PlotDF_miRNA_mRNA_correlations_SeedBased.RData"))
  
  pllist <- list()
  for(f in levels(as.factor(td$idx))){
    pllist[[f]] <- local({
      td1 <- td[td$idx==f,]
      td1$variable <- factor(td1$variable,levels=c("control","D7","D14"))
      scaleFactor <- (max(td1$mean)) / max(td1$mean1)
      td1$mean1 <- td1$mean1*scaleFactor
      td1$scaleFactor <- scaleFactor
      q <- ggplot(td1,aes(x=variable))+
        geom_errorbar(aes(y=mean,ymin=mean-sem,ymax=mean+sem,col="red",group="red"), width=.2,position=position_dodge(0.01))+
        geom_errorbar(aes(y=mean1,ymin=mean1-sem1,ymax=mean1+sem1,col="blue",group="blue"), width=.2,position=position_dodge(0.01))+
        geom_line(aes(y=mean,group="red",col="red"),size=1.2)+  
        geom_line(aes(y=mean1,group="blue",col="blue"),size=1.2)+  
        scale_y_continuous(name="expression, log2 FPKM", sec.axis=sec_axis(~./unique(td1$scaleFactor), name="editing levels")) +
        gg_aes+xlab("")+
        theme(axis.title.y.left=element_text(color="cyan4"),
              axis.text.y.left=element_text(color="cyan4"),
              axis.title.y.right=element_text(color="red"),
              axis.text.y.right=element_text(color="red"))+
        ggtitle(paste0(unique(td1$gene)), subtitle = paste0( "Edit: ",unique(td1$idx), "\n",unique(td1$miRs)," (",unique(td1$ref),")\np= ", round(unique(td1$cor),2)  ))
      print(q)
    })
  }
  
  grid.arrange(grobs=pllist[12],nrow=2)
  
}


## miRNA-mRNA correlation plots based on seedoverlaps using mTarBase validated targets
if(F){
  load(paste0(DATADIR,"PlotDF_miRNA_mRNA_correlations_SeedBased_forTarBaseValTargets.RData"))
  td <- td[abs(td$cor)>0.7,]
  levels(as.factor(as.character(td$idx)))
  
  pllist <- list()
  for(f in levels(as.factor(td$idx))){
    td1 <- td[td$idx==f,]
    td1$variable <- factor(td1$variable,levels=c("control","D7","D14"))
    scaleFactor <- (max(td1$mean)) / max(td1$mean1)
    td1$mean1 <- td1$mean1*scaleFactor
    td1$scaleFactor <- scaleFactor
    q <- ggplot(td1,aes(x=variable))+
      geom_errorbar(aes(y=mean,ymin=mean-sem,ymax=mean+sem,col="red",group="red"), width=.2,position=position_dodge(0.01))+
      geom_errorbar(aes(y=mean1,ymin=mean1-sem1,ymax=mean1+sem1,col="blue",group="blue"), width=.2,position=position_dodge(0.01))+
      geom_line(aes(y=mean,group="red",col="red"),size=1.2)+  
      geom_line(aes(y=mean1,group="blue",col="blue"),size=1.2)+  
      scale_y_continuous(name="expression", sec.axis=sec_axis(~./unique(scaleFactor), name="editing levels")) +
      gg_aes+xlab("")+
      theme(axis.title.y.left=element_text(color="cyan4"),
            axis.text.y.left=element_text(color="cyan4"),
            axis.title.y.right=element_text(color="red"),
            axis.text.y.right=element_text(color="red"))+
      ggtitle(paste0(unique(td1$gene)), subtitle = paste0( "Edit: ",unique(td1$idx), "\n",unique(td1$miRs)," (",unique(td1$ref),")\np= ", round(unique(td1$cor),2)  ))
    pllist[[f]] <- print(q)
    rm(td1,q,f,scaleFactor)
    
  }
  
  grid.arrange(grobs=pllist[1:4],nrow=2)
  
}


### Correlation between expression of validated targets of edited miRNA and editing levels of the edited miRNA
if(F){
  load(paste0(DATADIR,"PlotDF_EditedMiRNA_CorrBn_valTarsNeditingLevels.RData"))
  td <- sq %>% unique
  td <- td[abs(td$cor)>0.5| td$Target.Gene=="Wdr76",]
  td$id <- paste0(td$idx,td$Target.Gene)
  
  f=levels(as.factor(as.character(td$idx)))
  
  pllist <- list()
  for(f in levels(as.factor(as.character(td$id))) ){
    td1 <- td[td$id==f,]
    td1$variable <- factor(td1$variable,levels=c("control","D7","D14"))
    scaleFactor <- (max(td1$mean)) / max(td1$mean1)
    td1$mean1 <- td1$mean1*scaleFactor
    td1$scaleFactor <- scaleFactor
    q <- ggplot(td1,aes(x=variable))+
      geom_errorbar(aes(y=mean,ymin=mean-sem,ymax=mean+sem,col="red",group="red"), width=.2,position=position_dodge(0.01))+
      geom_errorbar(aes(y=mean1,ymin=mean1-sem1,ymax=mean1+sem1,col="blue",group="blue"), width=.2,position=position_dodge(0.01))+
      geom_line(aes(y=mean,group="red",col="red"),size=1.2)+  
      geom_line(aes(y=mean1,group="blue",col="blue"),size=1.2)+  
      scale_y_continuous(name="expression", sec.axis=sec_axis(~./unique(td1$scaleFactor), name="editing levels")) +
      gg_aes+xlab("")+
      theme(axis.title.y.left=element_text(color="cyan4"),
            axis.text.y.left=element_text(color="cyan4"),
            axis.title.y.right=element_text(color="red"),
            axis.text.y.right=element_text(color="red"))+
      ggtitle(paste0(unique(td1$Target.Gene)), subtitle = paste0( "Edit: ",unique(td1$idx), "\n",unique(td1$miRs),"\np= ", round(unique(td1$cor),2)  ))
    pllist[[f]] <- q
    #rm(td1,q,f,scaleFactor)
  }
  
  grid.arrange(grobs=pllist,nrow=1)
  
}


