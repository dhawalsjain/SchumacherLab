rm(list=ls())
source("functions.R")
source("vars.R")

###############################################################################################
### Editing in miRNA itself?
###############################################################################################
MIRFILE="A:/work/mm10_Annotations/mmu_mm10_miRNAGEnomicCoords.txt"

if(!file.exists(paste0(DATADIR,"EditedmiRNAs.RData"))){
  load(paste0(DATADIR,"DRE_DE_matrix.RData"))
  load(paste0(DATADIR,"EditingLevelReport_Anno.RData"))
  mirs <- read.delim(MIRFILE,header=F,comment.char = "#",stringsAsFactors = F)
  mirs$V1 <- gsub("chr","",mirs$V1)
  mirs$name <- gsub(";\\S*","",gsub("\\S*Name=","",mirs$V9))
  mirs <- with(mirs, GRanges(V1, IRanges(V4,V5),"+",anno=V3,name))
  
  hc <- as(hc_edits_anno,"GRanges")
  hc <- hc[,c("idx","sample_name")]
  
  mirs$control <- countOverlaps(mirs,unique(hc[hc$idx%in%hlpr[hlpr$control>0,]$idx]),ignore.strand=T)
  mirs$D7 <- countOverlaps(mirs,unique(hc[hc$idx%in%hlpr[hlpr$D7>0,]$idx]),ignore.strand=T)
  mirs$D14 <- countOverlaps(mirs,unique(hc[hc$idx%in%hlpr[hlpr$D14>0,]$idx]),ignore.strand=T)
  
  m <- as.data.frame(mirs)
  m <- m[rowSums(m[,8:10])>0,]
  m$control <- ifelse(m$control>0,as.character("edited"),NA)
  m$D7 <- ifelse(m$D7>0,as.character("edited"),NA)
  m$D14 <- ifelse(m$D14>0,as.character("edited"),NA)
  
  m <- as(m,"GRanges")
  o <- findOverlaps(m,hc,ignore.strand=T) %>% as.data.frame
  m <- as.data.frame(m)
  
  o <- cbind(m[o$queryHits,],idx=hc[o$subjectHits]$idx)
  o <- unique(o)
  
  editedmiRs <- o
  save(editedmiRs, file=paste0(DATADIR,"EditedmiRNAs.RData"))
}




###############################################################################################
###  Differential enrichment of miRNA targets
###############################################################################################
if(F){
  load(paste0(DATADIR,"EditingLevelReport_Anno.RData"))
  load(paste0(DATADIR,"DRE_DE_matrix.RData"))
  load(paste0(DATADIR,"miRSeedOverlaps_filtermiRNAs.RData"))
  
  r1 <- mirSeedOlp[mirSeedOlp$idx%in%hlpr[hlpr$control>0,]$idx,]
  r2 <- mirSeedOlp[mirSeedOlp$idx%in%hlpr[hlpr$D7>0,]$idx,]
  r3 <- mirSeedOlp[mirSeedOlp$idx%in%hlpr[hlpr$D14>0,]$idx,]
  
  r1 <- table(r1$mir,r1$seq) %>% as.data.frame
  r1 <- reshape2::dcast(r1, Var1~Var2,value.var = "Freq")
  r2 <- table(r2$mir,r2$seq) %>% as.data.frame
  r2 <- reshape2::dcast(r2, Var1~Var2,value.var = "Freq")
  r3 <- table(r3$mir,r3$seq) %>% as.data.frame
  r3 <- reshape2::dcast(r3, Var1~Var2,value.var = "Freq")
  
  names(r1) <- c("miR","control_gain","control_loss")
  names(r2) <- c("miR","D7_gain","D7_loss")
  names(r3) <- c("miR","D14_gain","D14_loss")
  
  myfun_enrich <- function(r1,r2,cond="gain"){
    if(cond == "gain"){
      x <- merge(r1[,c(1:2)],r2[,c(1:2)],by="miR",all.x=T,all.y=T)
    }else{
      x <- merge(r1[,c(1,3)],r2[,c(1,3)],by="miR",all.x=T,all.y=T)
    }
    
    x[,2] <- ifelse(is.na(x[,2]),0,x[,2])
    x[,3] <- ifelse(is.na(x[,3]),0,x[,3])
    x$baseT <- sum(x[,2])
    x$expT <- sum(x[,3])
    x$baseT <- x$baseT-x[,2]
    x$expT <- x$expT - x[,3]
    
    x$pval <- apply(x[,2:5],1,function(x){
      y <- fisher.test( matrix(x,nrow=2,byrow = T),alternative = "two.sided")$p.value
      return(y)
    })
    x$fdr <- p.adjust(x$pval,method = "fdr")
    x
  }
  
  d7ctr.gain <- myfun_enrich(r1,r2,cond = "gain")
  d7ctr.loss <- myfun_enrich(r1,r2,cond = "loss")
  d14ctr.gain <- myfun_enrich(r1,r3,cond = "gain")
  d14ctr.loss <- myfun_enrich(r1,r3,cond = "loss")
  
  p <- merge(d7ctr.gain[,1:3], d14ctr.gain[,c(1,3)],by="miR",all.x=T,all.y=T)
  p$sd <- apply(p[,2:4],1,sd,na.rm=T)
  q <- merge(d7ctr.loss[,1:3], d14ctr.loss[,c(1,3)],by="miR",all.x=T,all.y=T)
  q$sd <- apply(q[,2:4],1,sd,na.rm=T)
  
  
  miRgain <- p
  miRloss <- q
  save(miRgain,miRloss,file=paste0(DATADIR,"miRbinding_GainLoss_Timecourse.RData"))
  
  load(paste0(DATADIR,"genes.RData"))
  #gGenes <- gGenes[,c("chr", "strand", "gene", "geneStart", "geneEnd", "gene_id", "biotype")] %>% unique
  out <- reshape2::dcast(mirSeedOlp,names+mir+idx~seq)
  out$alt <- ifelse(is.na(out$alt),0,1)
  out$ref <- ifelse(is.na(out$ref),0,1)
  out <- out[rowSums(out[,4:5])<2,]
  out$alt <- ifelse(out$alt==0,"",as.character("gain"))
  out$ref <- ifelse(out$ref==0,"",as.character("loss"))
  out$txID <- gsub("_\\S*","",out$names)
  out <- merge(out, gGenes,by="txID",all.x=T)
  out <- out[,c( "gene", "mir", "idx", "alt", "ref", "biotype", "chr","txStart", "txEnd", "strand", 
                 "txID","gene_id")] %>% unique
  
  z <- imputed_edits
  z$id <- rownames(z)
  z$D7 <- round(apply(z[,1:3],1,gmean),2)
  z$D14 <- round(apply(z[,4:6],1,gmean),2)
  z$control <- round(apply(z[,7:9],1,gmean),2)
  z <- z[,10:13]
  names(z)[2:4] <- c("D7_editLevel","D14_editLevel","control_editLevel")
  uq <- hc_edits_anno[,c("id","idx","snpeff_uq","repeats","known","ConsScore","hg19")] %>% unique
  z <- merge(z,uq,by="id")
  z$id <- NULL
  out <- merge(out,z,by="idx")
  miRLGDF <- out
  save(miRLGDF,file=paste0(DATADIR,"miR_lossGain_dataframe.RData"))
  
}



### if there is a correlation between expression and editing for transcripts that have undergone RNAediting leading to miRNA binding loss or gain
if(F){
  load(paste0(DATADIR,"EditingLevelReport_Anno.RData"))
  load(paste0(DATADIR,"DRE_DE_matrix.RData"))
  load(paste0(DATADIR,"miRSeedOverlaps_filtermiRNAs.RData"))
  load(paste0(DATADIR,"miR_lossGain_dataframe.RData"))
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
  cq <- cnt
  cnt[,1:9] <- log2(cnt[,1:9]+1)
  rm(out0,out1,outC)
  names(cnt)[1:9] <- paste0(names(cnt)[1:9],"_Exp")
  
  z <- imputed_edits
  z$id <- rownames(z)
  
  miRLGDF$id <- gsub(":A>G|:T>C","",miRLGDF$idx)
  miRLGDF1 <- miRLGDF[,c("id","idx","gene","mir","alt","ref")]
  miRLGDF1 <- unique(miRLGDF1)
  
  td <- merge(miRLGDF1,z,by="id",all.x=T)
  td <- merge(td,cnt,by="gene")
  names(td)
  td <- td[,c("gene", "id", "idx", "mir", "alt", "ref", 
              "D7_R1", "D7_R2", "D7_R3", "D14_R1", "D14_R2", "D14_R3", 
              "control_R1","control_R2", "control_R3", 
              "D7_R1_Exp", "D7_R2_Exp", "D7_R3_Exp", 
              "D14_R1_Exp","D14_R2_Exp", "D14_R3_Exp", 
              "control_R1_Exp", "control_R2_Exp","control_R3_Exp")]
  td$cor <- apply(td[,7:24],1,function(x){cor(x[1:9],x[10:18])})
  td$ref <- ifelse(td$ref=="",as.character(td$alt),as.character(td$ref))
  td$alt <- NULL
  td <- td[!is.na(td$cor),]
  td$validity <- 0
  td[(td$ref=="gain" & td$cor<0|td$ref=="loss" & td$cor>0),]$validity <- 1
  save(td,file=paste0(DATADIR,"PlotDF_miRNA_mRNA_correlations_SeedBased_ALL.RData"))
  
  
  ### Add expression to mirLGDF
  cq$D7_Expn <- round(log2(apply(cq[,1:3],1,gmean)+1),2)
  cq$D14_Expn<- round(log2(apply(cq[,4:6],1,gmean)+1),2)
  cq$Control_Expn<- round(log2(apply(cq[,7:9],1,gmean)+1),2)
  cq <- cq[,11:13]
  cq$gene_id <- rownames(cq)
  miRLGDF <- merge(miRLGDF,cq,by="gene_id")
  miRLGDF <- merge(miRLGDF,unique(td[,c("idx","cor")]),by="idx",all.x=T)
  miRLGDF$ref <- ifelse(miRLGDF$ref=="",as.character(miRLGDF$alt),as.character(miRLGDF$ref))
  miRLGDF$alt <- NULL
  names(miRLGDF)[c(5,15)] <- c("miR_binding","annotation")  
  miRLGDF <- miRLGDF[,c("idx", "mir", "miR_binding",  "gene_id", "gene","txID","biotype", 
                        "chr", "txStart", "txEnd", "strand",  "D7_editLevel", 
                        "D14_editLevel", "control_editLevel", "D7_Expn", "D14_Expn", "Control_Expn",
                        "annotation", "repeats", "known", "ConsScore", "hg19","cor")] %>% unique
  miRLGDF <- miRLGDF[order(abs(miRLGDF$cor),decreasing = T),]
  
  load(file=paste0(DATADIR,"miRNA_geneExpression.RData"))
  mirExpn$D7_miRExpn <- round(log2(apply(mirExpn[,4:6],1,gmean)+1),2)
  mirExpn$D14_miRExpn<- round(log2(apply(mirExpn[,7:9],1,gmean)+1),2)
  mirExpn$Control_miRExpn<- round(log2(apply(mirExpn[,10:12],1,gmean)+1),2)
  mirExpn <- mirExpn[,c(1:3,13:15)] %>% unique
  mirExpn <- mirExpn[,c("miR_family","D7_miRExpn","D14_miRExpn","Control_miRExpn")] %>% unique
  miRLGDF <- merge(miRLGDF,mirExpn,by.x="mir",by.y="miR_family")
  miRLGDF <- miRLGDF[,c("mir", "idx", "miR_binding", "gene_id", "gene", "txID", "biotype", 
                        "chr", "txStart", "txEnd", "strand", "D7_editLevel", "D14_editLevel", 
                        "control_editLevel", "D7_Expn", "D14_Expn", "Control_Expn", "D7_miRExpn", 
                        "D14_miRExpn", "Control_miRExpn", "annotation", 
                        "repeats", "known", "ConsScore", "hg19", "cor")]
  save(miRLGDF,file=paste0(DATADIR,"miR_lossGain_dataframe_New.RData"))
  
  
  ## miRNA-mRNA influenced pairs, based on our data 
  td <- td[abs(td$cor)>0.8 & td$validity>0,]
  td$id <- NULL
  td <- data.table(td)
  td <- td[,miRs:=paste0(mir,collapse = ","),by=list(gene,idx,ref)]
  td <- as.data.frame(td)
  td$mir <- td$validity <- NULL
  td <- unique(td)
  td <- melt(td,measure.vars = names(td)[4:21])
  td$item <- "editing level"
  td[grep("_Exp",td$variable),]$item <- "expression"
  td$variable <- gsub("_\\S*","",td$variable)
  td <- data.table(td)
  td <- td[,mean:=mean(value),by=list(gene,idx,ref,miRs,variable,item)]
  td <- td[,sd:=sd(value),by=list(gene,idx,ref,miRs,variable,item)]
  td <- td[,sem:=sem(value),by=list(gene,idx,ref,miRs,variable,item)]
  td <- as.data.frame(td)
  td$value <- NULL
  td <- unique(td)
  td$id <- paste0(td$gene,td$idx,td$ref,td$miRs,td$variable)
  td1 <- td[td$item=="expression",]
  td2 <- td[td$item=="editing level",]
  names(td2)[8:10] <- c("mean1","sd1","sem1")
  td2 <- td2[,8:11]
  td <- merge(td1,td2,by="id")
  td$id <- NULL
  save(td,file=paste0(DATADIR,"PlotDF_miRNA_mRNA_correlations_SeedBased.RData"))
  
  
}


## miRNA-mRNA influenced pairs, based on validated miRTar base data
##  edit-levels and expression correlation pltos 
if(F){
  load(paste0(DATADIR,"PlotDF_miRNA_mRNA_correlations_SeedBased_ALL.RData"))
  load(paste0(DATADIR,"miRSeedOverlaps_filtermiRNAs.RData"))
  load(paste0(DATADIR,"genes.RData"))
  tars <- read.delim("A:/work/mm10_Annotations/miRTarBase.txt",header=T,stringsAsFactors = F)
  
  mirSeedOlp$txID <- gsub("_\\S*","",mirSeedOlp$names)
  mirSeedOlp <- merge(mirSeedOlp, unique(gGenes[,c("txID","gene")]),by="txID",all.x=T )
  mirSeedOlp$mir <- gsub("miR","mmu-miR",mirSeedOlp$mir)
  mirSeedOlp$test <- paste0(mirSeedOlp$gene,"|",mirSeedOlp$mir)
  tars$test <- paste0(tars$Target.Gene,"|",tars$miRNA)
  tars <- tars[tars$test%in%mirSeedOlp$test,]
  tars <- merge(tars,unique(mirSeedOlp[,c("idx","seq","test")]),by="test",all.x=T)
  ValTars <- reshape2::dcast(tars,test+Experiments+Support.Type+idx+Target.Gene~seq)
  ValTars$test <- gsub("mmu-","",ValTars$test)
  rm(tars)
  
  td$test <- paste0(td$gene,"|",td$mir)
  td <- td[td$test%in%ValTars$test,]
  
  td$test <- td$id <- NULL
  td <- unique(td)
  td <- data.table(td)
  td <- td[,miRs:=paste0(mir,collapse = ","),by=list(gene,idx,ref)]
  td <- as.data.frame(td)
  td$mir <- td$validity <- NULL
  td <- unique(td)
  td <- melt(td,measure.vars = names(td)[4:21])
  td$item <- "editing level"
  td[grep("_Exp",td$variable),]$item <- "expression"
  td$variable <- gsub("_\\S*","",td$variable)
  td <- data.table(td)
  td <- td[,mean:=mean(value),by=list(gene,idx,ref,miRs,variable,item)]
  td <- td[,sd:=sd(value),by=list(gene,idx,ref,miRs,variable,item)]
  td <- td[,sem:=sem(value),by=list(gene,idx,ref,miRs,variable,item)]
  td <- as.data.frame(td)
  td$value <- NULL
  td <- unique(td)
  td$id <- paste0(td$gene,td$idx,td$ref,td$miRs,td$variable)
  td1 <- td[td$item=="expression",]
  td2 <- td[td$item=="editing level",]
  names(td2)[8:10] <- c("mean1","sd1","sem1")
  td2 <- td2[,8:11]
  
  td <- merge(td1,td2,by="id")
  td$id <- NULL
  save(td,file=paste0(DATADIR,"PlotDF_miRNA_mRNA_correlations_SeedBased_forTarBaseValTargets.RData"))
  
  levels(as.factor(td$idx))
  
  
}

## miRNA-mRNA influenced pairs, based on validated miRTar base data for edited miRNAs
##  edit-levels and expression correlation pltos 

if(F){
  load(paste0(DATADIR,"EditingLevelReport_Anno.RData"))
  load(paste0(DATADIR,"EditedmiRNAs.RData"))
  load(paste0(DATADIR,"DRE_DE_matrix.RData"))
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
  rm(out0,out1,outC)
  names(cnt)[1:9] <- paste0(names(cnt)[1:9],"_Exp")
  
  z <- imputed_edits
  z$id <- rownames(z)
  
  e <- editedmiRs
  e$id <- gsub(":A>G|:T>C","",e$idx)
  e <- merge(e,z,by="id")
  e <- e[,c("id","idx","name","D7_R1", "D7_R2", "D7_R3", 
            "D14_R1", "D14_R2", "D14_R3", "control_R3", "control_R2", "control_R1")]
  tars <- read.delim("A:/work/mm10_Annotations/miRTarBase.txt",header=T,stringsAsFactors = F)
  tars <- tars[tars$miRNA%in%editedmiRs$name,]
  tars <- tars[,c("miRNA","Target.Gene","Experiments","Support.Type","References..PMID.")] %>% unique 
  
  sq <- merge(tars,cnt,by.x="Target.Gene",by.y="gene")
  sq <- merge(sq,e,by.x="miRNA",by.y="name")
  sq <- sq [,c("id", "idx","miRNA", "Target.Gene", "Experiments", "Support.Type", "References..PMID.", 
               "D7_R1_Exp", "D7_R2_Exp", "D7_R3_Exp", "D14_R1_Exp", "D14_R2_Exp", 
               "D14_R3_Exp", "control_R1_Exp", "control_R2_Exp", "control_R3_Exp", 
               "D7_R1", "D7_R2", "D7_R3", "D14_R1", "D14_R2", "D14_R3", 
               "control_R1", "control_R2", "control_R3")]
  
  sq$cor <- apply(sq[,8:25],1,function(x){cor(x[1:9],x[10:18])})
  sz <- sq 
  
  sq$References..PMID.<- "."
  sq$id <- NULL
  sq <- data.table(sq)
  sq <- sq[,miRs:=paste0(unique(miRNA),collapse = ","),by=list(idx,Target.Gene, Experiments)]
  sq <- as.data.frame(sq)
  sq$miRNA <-  NULL
  sq <- unique(sq)
  sq <- melt(sq,measure.vars = names(sq)[6:23])
  sq$item <- "editing level"
  sq[grep("_Exp",sq$variable),]$item <- "expression"
  sq$variable <- gsub("_\\S*","",sq$variable)
  
  sq <- data.table(sq)
  sq <- sq[,mean:=mean(value),by=list(Target.Gene,idx,miRs,variable,item)]
  sq <- sq[,sd:=sd(value),by=list(Target.Gene,idx,miRs,variable,item)]
  sq <- sq[,sem:=sem(value),by=list(Target.Gene,idx,miRs,variable,item)]
  sq <- as.data.frame(sq)
  sq$value <- NULL
  sq <- unique(sq)
  sq$id <- paste0(sq$Target.Gene,sq$idx,sq$miRs,sq$variable)
  
  sq1 <- sq[sq$item=="expression",]
  sq2 <- sq[sq$item=="editing level",]
  names(sq2)[10:12] <- c("mean1","sd1","sem1")
  sq2 <- sq2[,10:13]
  sq <- merge(unique(sq1),unique(sq2),by="id")
  sq$id <- NULL
  save(sq,file=paste0(DATADIR,"PlotDF_EditedMiRNA_CorrBn_valTarsNeditingLevels.RData"))
  
  
  
}


##########################################################################################################3
### miRNA predicted conserved targets
##########################################################################################################

if(F){
  mirs <- read.delim("A:/work/mm10_Annotations/miR_Family_Info.txt",header=T,stringsAsFactors = F)
  mirs <- mirs[grep("mmu",mirs$MiRBase.ID),]
  mirs$Mature.sequence <- gsub("U","T",mirs$Mature.sequence)
  #mirs <- mirs[mirs$MiRBase.Accession!="",]
  ## use the miRNA target binding sites predicted using TargetScan, mm10
  load(paste0(DATADIR,"EditingLevelReport_Anno.RData"))
  uqedits <- hc_edits_anno[,c("seqnames","start","strand","V4")] %>% unique
  uqedits <- with(uqedits,GRanges(seqnames,IRanges(start,start),strand,edit=V4))
  uqedits$id <- paste0(seqnames(uqedits),":",start(uqedits))
  imp <- imputed_edits
  imp[imp==0]<- NA
  editC <- uqedits[uqedits$id%in%rownames(imp[!is.na(imp$control_R3),])]
  editD7 <- uqedits[uqedits$id%in%rownames(imp[!is.na(imp$D7_R1),])]
  editD14 <- uqedits[uqedits$id%in%rownames(imp[!is.na(imp$D14_R1),])]
  rm(imp,uqedits)
  
  
  d <- read.delim("A:/work/mm10_Annotations/Predicted_Targets_miRNA.mm10.bed",header = F,stringsAsFactors = F)
  names(d)[1:6] <- c("chr","start","end","name","score","strand")
  d$chr <- gsub("chr","",d$chr)
  d$gene <- gsub(":\\S*","",d$name)
  d$miR <- gsub("\\S*:","",d$name)
  d <- as(d,"GRanges")
  
  o1 <- findOverlaps(editC,d,ignore.strand=T) %>% as.data.frame
  o1$mir <- d[o1$subjectHits]$miR
  o1$idx <- paste0(editC[o1$queryHits]$id,":",editC[o1$queryHits]$edit)
  o1$name <- d[o1$subjectHits]$name
  o1$cond <- "control"
  o2 <- findOverlaps(editD7,d,ignore.strand=T) %>% as.data.frame
  o2$mir <- d[o2$subjectHits]$miR
  o2$idx <- paste0(editD7[o2$queryHits]$id,":",editD7[o2$queryHits]$edit)
  o2$name <- d[o2$subjectHits]$name
  o2$cond <- "D7"
  o3 <- findOverlaps(editD14,d,ignore.strand=T) %>% as.data.frame
  o3$mir <- d[o3$subjectHits]$miR
  o3$idx <- paste0(editD14[o3$queryHits]$id,":",editD14[o3$queryHits]$edit)
  o3$name <- d[o3$subjectHits]$name
  o3$cond <- "D14"
  consv <- rbind(o1,o2,o3)
  save(consv,file=paste0(DATADIR,"ConservedMiRTargets_EnrichAlongEdits.RData"))
  
}

## plot consv
if(F){
  load(paste0(DATADIR,"EditingLevelReport_Anno.RData"))
  load(paste0(DATADIR,"DEExpression_DESeq2.RData"))
  load(paste0(DATADIR,"genes.RData"))
  load(paste0(DATADIR,"ConservedMiRTargets_EnrichAlongEdits.RData"))
  miRLGDF <- consv
  rm(consv)
  
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
  rm(out0,out1,outC)
  names(cnt)[1:9] <- paste0(names(cnt)[1:9],"_Exp")
  
  z <- imputed_edits
  z$id <- rownames(z)
  
  miRLGDF$id <- gsub(":A>G|:T>C","",miRLGDF$idx)
  miRLGDF$alt <- miRLGDF$ref <- "."
  miRLGDF$gene <- gsub(":\\S*","",miRLGDF$name)
  miRLGDF1 <- miRLGDF[,c("id","idx","gene","mir","alt","ref")]
  miRLGDF1 <- unique(miRLGDF1)
  
  td <- merge(miRLGDF1,z,by="id",all.x=T)
  td <- merge(td,cnt,by="gene")
  names(td)
  td <- td[,c("gene", "id", "idx", "mir", "alt", "ref", 
              "D7_R1", "D7_R2", "D7_R3", "D14_R1", "D14_R2", "D14_R3", 
              "control_R1","control_R2", "control_R3", 
              "D7_R1_Exp", "D7_R2_Exp", "D7_R3_Exp", 
              "D14_R1_Exp","D14_R2_Exp", "D14_R3_Exp", 
              "control_R1_Exp", "control_R2_Exp","control_R3_Exp")]
  td$cor <- apply(td[,7:24],1,function(x){cor(x[1:9],x[10:18])})
  td$alt <- NULL
  td <- td[!is.na(td$cor),]
  td$validity <- 0
  td[td$cor>0,]$validity <- 1
  ## miRNA-mRNA influenced pairs, based on our data 
  td <- td[abs(td$cor)>0.5 & td$validity>0,]
  td$id <- NULL
  td <- data.table(td)
  td <- td[,miRs:=paste0(mir,collapse = ","),by=list(gene,idx,ref)]
  td <- as.data.frame(td)
  td$mir <- td$validity <- NULL
  td <- unique(td)
  td <- melt(td,measure.vars = names(td)[4:21])
  
  td$item <- "editing level"
  td[grep("_Exp",td$variable),]$item <- "expression"
  td$variable <- gsub("_\\S*","",td$variable)
  td <- data.table(td)
  td <- td[,mean:=mean(value),by=list(gene,idx,ref,miRs,variable,item)]
  td <- td[,sd:=sd(value),by=list(gene,idx,ref,miRs,variable,item)]
  td <- td[,sem:=sem(value),by=list(gene,idx,ref,miRs,variable,item)]
  td <- as.data.frame(td)
  td$value <- NULL
  td <- unique(td)
  td$id <- paste0(td$gene,td$idx,td$ref,td$miRs,td$variable)
  td1 <- td[td$item=="expression",]
  td2 <- td[td$item=="editing level",]
  names(td2)[8:10] <- c("mean1","sd1","sem1")
  td2 <- td2[,8:11]
  td <- merge(td1,td2,by="id")
  td$id <- NULL
  save(td,file=paste0(DATADIR,"PlotDF_miRNA_mRNA_correlations_PredTargets.RData"))
  
}
  