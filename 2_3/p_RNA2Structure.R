source("functions.R")
source("vars.R")


###############################################################################################
### RNA structure
###############################################################################################
if(F){
  load(paste0(DATADIR,"miRpreds_intermediateTmp_Clustered.RData"))
  trxn=names(Tx)[8]
  gx <- GRanges(trxn,IRanges( range(start(g[seqnames(g)==trxn]))[1],
                              range(start(g[seqnames(g)==trxn]))[2])  )
  gx <- resize(gx, width(gx)+201,"center")
  
  pl<- cnt[cnt$tracking_id==trxn,]
  pl$gt <- gsub("_\\S*","",pl$condition)
  ggplot(pl, aes(x=gt, y=FPKM, fill=gt)) + 
    geom_boxplot()+geom_jitter(width = 0.2)+
    stat_compare_means(label = "p.signif", method = "t.test",ref.group = "Control") +
    theme_bw()+xlab("")+
    theme(axis.text = element_text(size=18,color="black"),
          axis.title = element_text(size=18,color="black"),
          legend.position = "none")
  
  pl <- imputed_edits[rownames(imputed_edits)%in%g[seqnames(g)==trxn]$id,]
  if(nrow(pl)>1){
    heatmap(t(t(pl)))
  }else{
    pl <- melt(pl)
    pl$gt <- gsub("_\\S*","",pl$variable)
    ggplot(pl, aes(x=gt, y=value, fill=gt)) + 
      geom_boxplot()+geom_jitter(width = 0.2)+
      stat_compare_means(label = "p.signif", method = "t.test",ref.group = "Control") +
      theme_bw()+xlab("")+ylab("editing ratio")+
      theme(axis.text = element_text(size=18,color="black"),
            axis.title = element_text(size=18,color="black"),
            legend.position = "none")
  }
  
  
  
  rm(outr,o1,edb,dna,hc_edits_anno,hlpr,imputed_edits)
}



