rm(list=ls())
source("functions.R")
source("vars.R")

#############################################################################################
if(F){
  load(paste0(DATADIR,"EditingLevelReport_Anno.RData"))
  load(paste0(DATADIR,"DRE_DE_matrix.RData"))
  load(paste0(DATADIR,"SplicingReport_SUPPA.RData"))
  load(paste0(DATADIR,"genes.RData"))
  genes <- gGenes[,c("strand","gene","gene_id","biotype")] %>% unique
  rrr <- res[,c("idx","logFC","change","adj.P.Val","comparison")] %>% unique
  names(rrr) <- c("Editing_Event","DRE_Log2FC","DRE_change","DRE_FDR","comparison")
  
  r <- as(res,"GRanges")
  r1 <- r[r$comparison=="D7_control" & r$adj.P.Val<0.05]
  r2 <- r[r$comparison=="D14_control" & r$adj.P.Val<0.05]
  
  pl  <- c()
  q1 <- spl[grep("D7_Control",spl$test)]
  olp <- findOverlaps(r1,resize(q1,101,"center"),ignore.strand=T) %>% as.data.frame
  q1$gene <- q1$idx <- NA 
  q1[olp$subjectHits]$gene <- r1[olp$queryHits]$gene
  q1[olp$subjectHits]$idx <- r1[olp$queryHits]$idx
  q1 <- as.data.frame(q1)
  q1$comparison <- "D7_control"
  pl <- rbind(pl, q1)
  
  q1 <- spl[grep("D14_Control",spl$test)]
  olp <- findOverlaps(r2,resize(q1,101,"center"),ignore.strand=T) %>% as.data.frame
  q1$gene <- q1$idx <- NA 
  q1[olp$subjectHits]$gene <- r2[olp$queryHits]$gene
  q1[olp$subjectHits]$idx <- r2[olp$queryHits]$idx
  q1 <- as.data.frame(q1)
  q1$comparison <- "D14_control"
  pl <- rbind(pl, q1)
  rm(q1,r1,r2,r)
  
  #pl$gene <- ifelse(is.na(pl$gene),"",as.character(pl$gene))
  ## plots
  pl$test <- gsub("_\\S*","",pl$test)
  pl$label <- paste0(pl$gene,"(",pl$test,")")
  pl$label <- ifelse(is.na(pl$gene),as.character(""),as.character(pl$label))
  pl$change <- as.numeric(as.character(pl$change))
  pl$pval <- as.numeric(as.character(pl$pval))
  
  
  
  pl1 <- pl[pl$comparison=="D14_control",]
  pl2 <- pl1[pl1$pval<=0.1 & !is.na(pl1$gene),]
  gg <- ggplot(pl1,aes(x=change,y=-log10(pval)))+geom_point(size=0.25)+theme_bw()+
    geom_point(data = pl2,aes(x=change,y=-log10(pval),label=label),col="red",size=2.5,inherit.aes = F)+
    geom_text_repel(aes(x=change,y=-log10(pval),label=label),data = pl2,inherit.aes = F,point.padding = 0.2,
                    col="green4",size=8)+
    theme(axis.text = element_text(size=14,colour = "black"),
          axis.title = element_text(size=14,colour = "black"),
          legend.text = element_text(size=13,colour = "black"),
          legend.title = element_text(size=14,color = "black"),
          plot.title = element_text(size=18,colour = "black",hjust = 0.5))+
    xlab("change PSI")+ggtitle("D14 vs Control")
  gg
  
  
  pl1 <- pl[pl$comparison=="D7_control",]
  pl2 <- pl1[pl1$pval<=0.1 & !is.na(pl1$gene),]
  qq <- ggplot(pl1,aes(x=change,y=-log10(pval)))+geom_point(size=0.25)+theme_bw()+
    geom_point(data = pl2,aes(x=change,y=-log10(pval),label=label),col="red",size=2.5,inherit.aes = F)+
    geom_text_repel(aes(x=change,y=-log10(pval),label=label),data = pl2,inherit.aes = F,point.padding = 0.2,
                    col="green4",size=8)+
    theme(axis.text = element_text(size=14,colour = "black"),
          axis.title = element_text(size=14,colour = "black"),
          legend.text = element_text(size=13,colour = "black"),
          legend.title = element_text(size=14,color = "black"),
          plot.title = element_text(size=18,colour = "black",hjust = 0.5))+
    xlab("change PSI")+ggtitle("D7 vs Control")
  qq
  
  pl1 <- pl1[,c("seqnames", "start", "end", "change", "pval", "test", "ensembl", "idx")]
  names(pl1) <- c("chr","AS_start","AS_end","AS_change","Pvalue","Type_of_AS","gene_id", "Editing_Event")
  pl1 <- merge(pl1,genes,by="gene_id",all.x=T)
  pl1 <- merge(pl1,rrr[rrr$comparison=="D14_control",],by="Editing_Event",all.x=T)
  pl2 <- pl1[pl1$Pvalue<0.05 & !is.na(pl1$Editing_Event),]
  
  
}


#############################################################################################