rm(list=ls())
source("functions.R")
source("vars.R")


###############################################################################################
### initial plots cmparing editing annotations, per replicate
###############################################################################################

if(F){
  load(paste0(DATADIR,"EditingLevelReport_Anno.RData"))
  hc_edits_anno$condition <- gsub("_\\S*","",hc_edits_anno$sample_name)
  
  pl <- c()
  for(f in c("control","D7","D14")){
    r <- hc_edits_anno[grep(f,hc_edits_anno$sample_name),]
    r <- reshape2::dcast(r,idx~sample_name,value.var = "width")
    r$test <- apply(r[,2:4],1,function(x) sum(is.na(x)))
    r <- r[r$test<2,]
    r$test <- NULL
    pl <- rbind(pl,data.frame(perc=sum(is.na(r))*100/(nrow(r)*3),prof=f))
    rm(r)
  }
  ggplot(pl, aes(x=prof,y=perc,fill=prof,col=prof,alpha=0.2))+
    geom_bar(stat = "identity",position="dodge")+
    xlab("")+ylab("% of total sites\n that are imputed")+
    scale_fill_manual(values = cols)+scale_color_manual(values=cols)+
    theme_bw()+gg_aes
  
  
  
  libs$sample <- unlist(l1[as.character(libs$profile)])
  libs$condition <- gsub("_\\S*","",libs$sample)
  ggplot(libs, aes(x=sample,y=depths,fill=condition,col=condition,alpha=0.2))+
    geom_bar(stat = "identity",position="dodge")+
    xlab("")+ylab("number of sequences \n in Millions")+
    scale_fill_manual(values = cols)+scale_color_manual(values=cols)+
    theme_bw()+gg_aes+
    theme(axis.text.x = element_text(angle = 90,hjust = 1))
  
  
  ## original data
  
  ggplot(hc_edits_anno, aes(x=sample_name,y=H,fill=sample_name,alpha=0.2)) + geom_violin() + theme_bw()+
    geom_boxplot(width=0.1, fill="white",outlier.size = NULL)+
    scale_color_manual(values=cols)+scale_fill_manual(values = cols)+
    #stat_compare_means(label = "p.signif", method = "t.test",ref.group = "control_R1")+
    xlab("")+ylab("editing levels\n without imputation")+gg_aes+
    theme(axis.text.x = element_text(angle = 90))+
    ggtitle("individual samples")
  
  
  ## % of edits known
  pl <- table(hc_edits_anno$sample_name,hc_edits_anno$known) %>% as.data.frame
  pl$Var2 <- ifelse(pl$Var2==0, as.character("novel"),as.character("known"))
  pl$Var2 <- factor(pl$Var2,levels=c("novel","known"))
  ggplot(pl, aes(fill=Var2, y=Freq, x=Var1,alpha=0.2)) + geom_bar(position="stack", stat="identity")+
    scale_fill_brewer(palette = "Dark2")+
    theme_bw()+gg_aes+theme(axis.text.x = element_text(angle = 90),legend.position = "top")+
    xlab("")+ylab("number of edits")
  
  ## edits in Alu
  pl <- table(hc_edits_anno$repeats,hc_edits_anno$sample_name) %>% as.data.frame
  pl$Var1 <- factor(pl$Var1,levels = (c("Nonrepetitive", "Alu", "Repetitive non-Alu")))
  ggplot(pl, aes(fill=Var1, y=Freq, x=Var2,alpha=0.2)) + geom_bar(position="stack", stat="identity")+
    scale_fill_brewer(palette = "Dark2")+
    theme_bw()+gg_aes+theme(axis.text.x = element_text(angle = 90),
                            legend.position = "top")+
    xlab("")+ylab("number of edits")
  
  ## Genomic annotations, SnpEff
  categories <- c("3_prime_UTR", "5_prime_UTR", "intergenic_region", "intron", 
                  "missense", "splice_","synonymous")
  pl <- c()
  for(f in categories){
    for(prof in levels(hc_edits_anno$sample_name)) {
      qwe <- hc_edits_anno[hc_edits_anno$sample_name==prof,]
      pl <- rbind(pl,
                  data.frame(category=f,profile=prof, size=length(unique(grep(f,qwe$snpeff)))))
      rm(qwe)
    }
  }
  pl$category <- gsub("_prime_UTR","'UTR",pl$category)
  pl$category <- gsub("_region","",pl$category)
  pl$category  <- gsub("splice_","splice_region",pl$category)
  
  ggplot(pl, aes(fill=category, y=size, x=profile,alpha=0.5)) + geom_bar(position="fill", stat="identity")+
    scale_fill_brewer(palette = "Dark2")+
    theme_bw()+gg_aes+theme(axis.text.x = element_text(angle = 90),
                            legend.position = "top")+
    xlab("")+ylab("number of edits")
  
  
  
  
  pl <- c()
  for(prof in levels(hc_edits_anno$sample_name)) {
    qwe <- hc_edits_anno[hc_edits_anno$sample_name==prof,]
    pl <- rbind(pl,
                cbind(profile=prof,as.data.frame(table(qwe$snpeff_uq)) ))
    rm(qwe)
  }
  names(pl)[2:3] <- c("category","size")
  pl$category <- gsub("_prime_UTR","'UTR",pl$category)
  pl$category <- gsub("_region","",pl$category)
  pl$category  <- gsub("splice","splice_region",pl$category)
  
  ggplot(pl, aes(fill=category, y=size, x=profile,alpha=0.5)) + geom_bar(position="fill", stat="identity")+
    scale_fill_brewer(palette = "Dark2")+
    theme_bw()+gg_aes+theme(axis.text.x = element_text(angle = 90),
                            legend.position = "top")+
    xlab("")+ylab("number of edits")
  
  
  
  
  sps <- imputed_edits
  sps$id <- rownames(sps)
  sps <- merge(sps,unique(hc_edits_anno[,c("id","V4","known","repeats","snpeff","aaswap")]))
  names(sps)[11:15] <- c("Edit","is_known","repeat","snpeff_annotation","AA_change")
  sps$in_Intron <- sps$in_3UTR <- sps$in_5UTR <- sps$in_spliceRegion <- sps$is_missense <- sps$is_intergenic <- sps$synony <- "no"
  sps[grep("intron",sps$snpeff_annotation),]$in_Intron <- as.character("yes")
  sps[grep("3_prime_UTR",sps$snpeff_annotation),]$in_3UTR <- as.character("yes")
  sps[grep("5_prime_UTR",sps$snpeff_annotation),]$in_5UTR <- as.character("yes")
  sps[grep("intergenic_region",sps$snpeff_annotation),]$is_intergenic <- as.character("yes")
  sps[grep("missense",sps$snpeff_annotation),]$is_missense <- as.character("yes")
  sps[grep("splice_",sps$snpeff_annotation),]$in_spliceRegion <- as.character("yes")
  sps[grep("synonymous",sps$snpeff_annotation),]$synony <- as.character("yes")
  
}


##### Clustered sites
if(F){
  ### hyper-edited loci
  myf_getClusteredEdits <- function(Z1){
    unv <- reduce(resize(Z1,51,"center"),ignore.strand=T)
    unv <- unv[width(unv)>51]
    unv$id <- 1:length(unv)
    o <- findOverlaps(Z1,unv,ignore.strand=T) %>% as.data.frame
    Z1$clusterid <- NA
    Z1[o$queryHits]$clusterid <- unv[o$subjectHits]$id
    Z1 <- Z1[!is.na(Z1$clusterid)]
    Z1
  }
  hc <- c()
  
  Z <- as(hc_edits_anno,"GRanges")
  Z <- Z[grep("control",Z$sample_name)]
  a1 <- myf_getClusteredEdits(unique(Z))
  Z <- Z[Z$idx%in%a1$idx] %>% as.data.frame
  hc <- rbind(hc,Z)
  rm(Z,a1)
  Z <- as(hc_edits_anno,"GRanges")
  Z <- Z[grep("D7",Z$sample_name)]
  a1 <- myf_getClusteredEdits(unique(Z))
  Z <- Z[Z$idx%in%a1$idx] %>% as.data.frame
  hc <- rbind(hc,Z)
  rm(Z,a1)
  Z <- as(hc_edits_anno,"GRanges")
  Z <- Z[grep("D14",Z$sample_name)]
  a1 <- myf_getClusteredEdits(unique(Z))
  Z <- Z[Z$idx%in%a1$idx] %>% as.data.frame
  hc <- rbind(hc,Z)
  rm(Z,a1)
  
  
  ## original data
  ggplot(hc, aes(x=sample_name,y=H,fill=sample_name,alpha=0.2)) + geom_violin() + theme_bw()+
    geom_boxplot(width=0.1, fill="white",outlier.size = NULL)+
    scale_color_manual(values=cols)+scale_fill_manual(values = cols)+
    #stat_compare_means(label = "p.signif", method = "t.test",ref.group = "control_R1")+
    xlab("")+ylab("editing levels\n without imputation")+gg_aes+
    theme(axis.text.x = element_text(angle = 90))+
    ggtitle("individual samples")
  
  
  ## % of edits known
  pl <- table(hc$sample_name,hc$known) %>% as.data.frame
  pl$Var2 <- ifelse(pl$Var2==0, as.character("novel"),as.character("known"))
  pl$Var2 <- factor(pl$Var2,levels=c("novel","known"))
  ggplot(pl, aes(fill=Var2, y=Freq, x=Var1,alpha=0.2)) + geom_bar(position="stack", stat="identity")+
    scale_fill_brewer(palette = "Dark2")+
    theme_bw()+gg_aes+theme(axis.text.x = element_text(angle = 90),legend.position = "top")+
    xlab("")+ylab("number of edits")
  
  ## edits in Alu
  pl <- table(hc$repeats,hc$sample_name) %>% as.data.frame
  pl$Var1 <- factor(pl$Var1,levels = (c("Nonrepetitive", "Alu", "Repetitive non-Alu")))
  ggplot(pl, aes(fill=Var1, y=Freq, x=Var2,alpha=0.2)) + geom_bar(position="stack", stat="identity")+
    scale_fill_brewer(palette = "Dark2")+
    theme_bw()+gg_aes+theme(axis.text.x = element_text(angle = 90),
                            legend.position = "top")+
    xlab("")+ylab("number of edits")
  
  ## Genomic annotations, SnpEff
  categories <- c("3_prime_UTR", "5_prime_UTR", "intergenic_region", "intron", 
                  "missense", "splice_","synonymous")
  pl <- c()
  for(f in categories){
    for(prof in levels(hc$sample_name)) {
      qwe <- hc[hc$sample_name==prof,]
      pl <- rbind(pl,
                  data.frame(category=f,profile=prof, size=length(unique(grep(f,qwe$snpeff)))))
      rm(qwe)
    }
  }
  pl$category <- gsub("_prime_UTR","'UTR",pl$category)
  pl$category <- gsub("_region","",pl$category)
  pl$category  <- gsub("splice_","splice_region",pl$category)
  
  ggplot(pl, aes(fill=category, y=size, x=profile,alpha=0.5)) + geom_bar(position="fill", stat="identity")+
    scale_fill_brewer(palette = "Dark2")+
    theme_bw()+gg_aes+theme(axis.text.x = element_text(angle = 90),
                            legend.position = "top")+
    xlab("")+ylab("number of edits")
  
  
  pl <- c()
  for(prof in levels(hc$sample_name)) {
    qwe <- hc_edits_anno[hc$sample_name==prof,]
    pl <- rbind(pl,
                cbind(profile=prof,as.data.frame(table(qwe$snpeff_uq)) ))
    rm(qwe)
  }
  names(pl)[2:3] <- c("category","size")
  pl$category <- gsub("_prime_UTR","'UTR",pl$category)
  pl$category <- gsub("_region","",pl$category)
  pl$category  <- gsub("splice","splice_region",pl$category)
  
  ggplot(pl, aes(fill=category, y=size, x=profile,alpha=0.5)) + geom_bar(position="fill", stat="identity")+
    scale_fill_brewer(palette = "Dark2")+
    theme_bw()+gg_aes+theme(axis.text.x = element_text(angle = 90),
                            legend.position = "top")+
    xlab("")+ylab("number of edits")
  
}

###############################################################################################
### initial plots cmparing editing annotations, combined replicate
###############################################################################################
if(F){
  load(paste0(DATADIR,"EditingLevelReport_Anno.RData"))
  hc_edits_anno$condition <- gsub("_\\S*","",hc_edits_anno$sample_name)
  
  uqedits <- hc_edits_anno[,c("seqnames","start","strand","V4")] %>% unique
  uqedits <- with(uqedits,GRanges(seqnames,IRanges(start,start),strand,edit=V4))
  uqedits$id <- paste0(seqnames(uqedits),":",start(uqedits))
  uqedits$idx <- paste0(uqedits$id,":",uqedits$edit)
  imp <- imputed_edits
  imp[imp==0]<- NA
  editC <- uqedits[uqedits$id%in%rownames(imp[!is.na(imp$control_R3),])] %>% as.data.frame
  editD7 <- uqedits[uqedits$id%in%rownames(imp[!is.na(imp$D7_R1),])] %>% as.data.frame
  editD14 <- uqedits[uqedits$id%in%rownames(imp[!is.na(imp$D14_R1),])] %>% as.data.frame
  rm(imp,uqedits)
  
  
  ggplot(hc_edits_anno, aes(x=condition,y=H,fill=condition,alpha=0.2)) + geom_violin() + theme_bw()+
    geom_boxplot(width=0.1, fill="white",outlier.size = NULL)+
    scale_color_manual(values=cols)+scale_fill_manual(values = cols)+
    stat_compare_means(label = "p.signif", method = "t.test",ref.group = "control")+
    xlab("")+ylab("editing levels\n without imputation")+gg_aes+
    ggtitle("merged replicates")
  
  
  ## known/unknown
  pl <- rbind(cbind(hlpr[hlpr$control>0,],gt="control\n(3702)"),
              cbind(hlpr[hlpr$D7>0,],gt="D7\n(1677)"),
              cbind(hlpr[hlpr$D14>0,],gt="D14\n(1398)"))
  ggplot(pl, aes(fill=as.factor(known), x=gt,alpha=0.2)) + geom_bar(position="fill")+
    scale_fill_brewer(palette = "Dark2")+
    theme_bw()+gg_aes+theme(legend.position = "top")+
    xlab("")+ylab("fraction of edits")
  
  ## repeats
  ggplot(pl, aes(fill=as.factor(repeats), x=gt,alpha=0.2)) + geom_bar(position="fill")+
    scale_fill_brewer(palette = "Dark2")+
    theme_bw()+gg_aes+theme(legend.position = "top")+
    xlab("")+ylab("fraction of edits")
  
  ## combined annotations
  categories <- c("3_prime_UTR", "5_prime_UTR", "intergenic_region", "intron", 
                  "missense", "synonymous", "splice")
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
  ggplot(qq, aes(fill=category, y=size, x=profile,alpha=0.2)) + geom_bar(position="fill", stat="identity")+
    scale_fill_brewer(palette = "Dark2")+
    theme_bw()+gg_aes+theme(legend.position = "top")+
    xlab("")+ylab("number of edits")
  
  
  
  ### UpsetR plots
  m <- imputed_edits
  m <- ifelse(m>0,1,0) %>% as.data.frame
  m$D7 <- rowSums(m[,1:3])
  m$D14 <- rowSums(m[,4:6])
  m$Control <- rowSums(m[,7:9])
  m <- m[,10:12]
  m <- ifelse(m>0,1,0) %>% as.data.frame
  m$id <- rownames(m)
  m <- merge(m,unique(hc_edits_anno[,c("id","known","repeats")]),by="id")
  m$known <- ifelse(m$known==1,"known","novel")
  
  i<-0
  mylist<-list()
  vectorUniqueValue <- unique(m$known)
  while ( length(vectorUniqueValue)>0 ){
    i<-i+1
    mylist[[ i ]]<-list(query = elements, params = list("known",as.character(vectorUniqueValue)), 
                        active = T,query.name=as.character(vectorUniqueValue)[1])
    vectorUniqueValue<-vectorUniqueValue[-1]
  }
  upset(m, sets = names(m)[2:4], mb.ratio = c(0.5, 0.5), order.by = c("degree"),
        text.scale=rep(1.5,6),point.size=5,matrix.color = "steelblue4",
        queries = mylist,query.legend="bottom",color.pal = 1)
  rm(i,mylist,vectorUniqueValue)
  
  i<-0
  mylist<-list()
  vectorUniqueValue <- unique(m$repeats)
  while ( length(vectorUniqueValue)>0 ){
    i<-i+1
    mylist[[ i ]]<-list(query = elements, params = list("repeats",as.character(vectorUniqueValue)), 
                        active = T,query.name=as.character(vectorUniqueValue)[1])
    vectorUniqueValue<-vectorUniqueValue[-1]
  }
  upset(m, sets = names(m)[2:4], mb.ratio = c(0.5, 0.5), order.by = c("degree"),
        text.scale=rep(1.5,6),point.size=5,matrix.color = "steelblue4",
        queries = mylist,query.legend="bottom",color.pal = 1)
  rm(i,mylist,vectorUniqueValue)
  
}

## clustered
if(F){
  ### hyper-edited loci
  myf_getClusteredEdits <- function(Z1){
    unv <- reduce(resize(Z1,51,"center"),ignore.strand=T)
    unv <- unv[width(unv)>51]
    unv$id <- 1:length(unv)
    o <- findOverlaps(Z1,unv,ignore.strand=T) %>% as.data.frame
    Z1$clusterid <- NA
    Z1[o$queryHits]$clusterid <- unv[o$subjectHits]$id
    Z1 <- Z1[!is.na(Z1$clusterid)]
    Z1
  }
  hc <- c()
  hl <- c()
  
  Z <- as(hc_edits_anno,"GRanges")
  Z <- Z[grep("control",Z$sample_name)]
  a1 <- myf_getClusteredEdits(unique(Z))
  Z <- Z[Z$idx%in%a1$idx] %>% as.data.frame
  hc <- rbind(hc,Z)
  hl <- rbind(hl,unique(hlpr[hlpr$idx%in%Z$idx & hlpr$control>0,]))
  rm(Z,a1)
  Z <- as(hc_edits_anno,"GRanges")
  Z <- Z[grep("D7",Z$sample_name)]
  a1 <- myf_getClusteredEdits(unique(Z))
  Z <- Z[Z$idx%in%a1$idx] %>% as.data.frame
  hc <- rbind(hc,Z)
  hl <- rbind(hl,hlpr[hlpr$idx%in%Z$idx & hlpr$D7>0,])
  rm(Z,a1)
  Z <- as(hc_edits_anno,"GRanges")
  Z <- Z[grep("D14",Z$sample_name)]
  a1 <- myf_getClusteredEdits(unique(Z))
  Z <- Z[Z$idx%in%a1$idx] %>% as.data.frame
  hc <- rbind(hc,Z)
  hl <- rbind(hl,hlpr[hlpr$idx%in%Z$idx & hlpr$D14>0,])
  rm(Z,a1)
  hl <- unique(hl)
  
}

####################################################################################################
## Editing annotation with respect to transcript biotypes
####################################################################################################

if(F){
  load(paste0(DATADIR,"EditingLevelReport_Anno.RData"))
  load(paste0(DATADIR,"genes.RData"))
  hc_edits_anno$condition <- gsub("_\\S*","",hc_edits_anno$sample_name)
  Z <- as(hc_edits_anno,"GRanges")
  genes <- with(gGenes,GRanges(chr,IRanges(geneStart,geneEnd),strand,biotype,gene))
  genes <- unique(genes)
  o <- findOverlaps(Z,genes,ignore.strand=T)
  Z$biotype <- NA
  Z$gene <- NA
  Z[queryHits(o)]$biotype <- genes[subjectHits(o)]$biotype
  Z[queryHits(o)]$gene <- genes[subjectHits(o)]$gene
  hc_edits_anno <- as.data.frame(Z)
  Z <- unique(hc_edits_anno[,c("gene","biotype","idx")])
  hlpr <- merge(hlpr,Z,by="idx")
  
  pl <- rbind(
  cbind(as.data.frame(table(hlpr[hlpr$control>0,]$biotype)),profile="control"),
  cbind(as.data.frame(table(hlpr[hlpr$D7>0,]$biotype)),profile="D7"),
  cbind(as.data.frame(table(hlpr[hlpr$D14>0,]$biotype)),profile="D14"),
  cbind(as.data.frame(table(genes$biotype)),profile="Background"))
  
  pl <- data.table(pl)
  pl <- pl[,tot:=sum(Freq),by=list(profile)]
  pl <- as.data.frame(pl)
  pl$perc <- round(pl$Freq/pl$tot,2)
  pl <- pl[pl$perc>0,]
  
  ggplot(pl,aes(x=Var1,y=perc,fill=profile))+geom_bar(stat="identity",position="dodge")+
    gg_aes+ coord_flip()+scale_fill_manual(values = cols)+
    xlab("")+  theme(legend.position = "top")
  
  
}

 ## clustered events
if(F){
  load(paste0(DATADIR,"EditingLevelReport_Anno.RData"))
  load(paste0(DATADIR,"genes.RData"))
  hc_edits_anno$condition <- gsub("_\\S*","",hc_edits_anno$sample_name)
  Z <- as(hc_edits_anno,"GRanges")
  genes <- with(gGenes,GRanges(chr,IRanges(geneStart,geneEnd),strand,biotype,gene))
  genes <- unique(genes)
  o <- findOverlaps(Z,genes,ignore.strand=T)
  Z$biotype <- NA
  Z$gene <- NA
  Z[queryHits(o)]$biotype <- genes[subjectHits(o)]$biotype
  Z[queryHits(o)]$gene <- genes[subjectHits(o)]$gene
  hc_edits_anno <- as.data.frame(Z)
  Z <- unique(hc_edits_anno[,c("gene","biotype","idx")])
  hlpr <- merge(hlpr,Z,by="idx")
  
  
  Z <- as(hc_edits_anno,"GRanges")
  ## Control clusters
  uqc <- Z[grep("control",Z$sample_name),]
  uqc <- unique(uqc[,c("id","idx")])
  d1 <- reduce(resize(uqc,101,"center"))
  d1 <- resize(d1,width(d1)-100,"center")
  d1 <- d1[width(d1)>50]
  uqc$clustered <- countOverlaps(uqc,d1,ignore.strand=T)
  uqc <- uqc[uqc$clustered>0]
  hlpr$control_clustered <- ifelse(hlpr$idx%in%uqc$idx,1,0)
  sum(hlpr$control_clustered)
  ## D7 clusters
  uqc <- Z[grep("D7",Z$sample_name),]
  uqc <- unique(uqc[,c("id","idx")])
  d1 <- reduce(resize(uqc,101,"center"))
  d1 <- resize(d1,width(d1)-100,"center")
  d1 <- d1[width(d1)>50]
  uqc$clustered <- countOverlaps(uqc,d1,ignore.strand=T)
  uqc <- uqc[uqc$clustered>0]
  hlpr$D7_clustered <- ifelse(hlpr$idx%in%uqc$idx,1,0)
  ## D14 clusters
  uqc <- Z[grep("D14",Z$sample_name),]
  uqc <- unique(uqc[,c("id","idx")])
  d1 <- reduce(resize(uqc,101,"center"))
  d1 <- resize(d1,width(d1)-100,"center")
  d1 <- d1[width(d1)>50]
  uqc$clustered <- countOverlaps(uqc,d1,ignore.strand=T)
  uqc <- uqc[uqc$clustered>0]
  hlpr$D14_clustered <- ifelse(hlpr$idx%in%uqc$idx,1,0)
  
  pl <- rbind(
    cbind(as.data.frame(table(hlpr[hlpr$control>0 & hlpr$control_clustered>0,]$biotype)),profile="control"),
    cbind(as.data.frame(table(hlpr[hlpr$D7>0 & hlpr$D7_clustered>0,]$biotype)),profile="D7"),
    cbind(as.data.frame(table(hlpr[hlpr$D14>0 & hlpr$D14_clustered>0,]$biotype)),profile="D14"),
    cbind(as.data.frame(table(genes$biotype)),profile="Background"))
  
  pl <- data.table(pl)
  pl <- pl[,tot:=sum(Freq),by=list(profile)]
  pl <- as.data.frame(pl)
  pl$perc <- round(pl$Freq/pl$tot,2)
  pl <- pl[pl$perc>0,]
  
  ggplot(pl,aes(x=Var1,y=perc,fill=profile))+geom_bar(stat="identity",position="dodge")+
    gg_aes+ coord_flip()+scale_fill_manual(values = cols)+ylab("fraction of total edits")+
    xlab("")+  theme(legend.position = "top")
  
  
}


### correlation b/w editing levels and expression
if(F){
  load(paste0(DATADIR,"readDepth_EditLevels.RData"))
   
  oDepth$value <- oDepth$alt/oDepth$depth
    oDepth %>%
      ggplot(aes(log2(depth+1), value)) +geom_point() +
      stat_bin_hex(geom="point", aes(color=log2(..count..)), binwidth=c(.05,.05), size=1, shape=20, stroke=0)+
      scale_color_viridis(option="inferno",direction = 1)+
      geom_smooth(method = "lm") +facet_wrap(~name,nrow = 1)+gg_aes+
      xlab("expression level, (read depth on log2 scale)")+ylab("editing level")+
    stat_cor(method = "pearson", size=5,col="red",label.sep = "\n",
             label.x.npc = "right", label.y.npc = "top")
    
  
}



### Events with Wt1 peaks
if(F){
  load(paste0(DATADIR,"EditingLevelReport_Anno.RData"))
  load(paste0(DATADIR,"DRE_DE_matrix.RData"))
  load(paste0(DATADIR,"Distance2Wt1Peaks_plotDF.RData"))
  
  rd$set <- gsub("edits","RNA-edit site",rd$set)
  rd$variable <- factor(rd$variable,levels=c("Control","D7","D14"))
  
  ggplot(rd, aes(value,col=set)) + stat_ecdf(geom = "point")+ 
    scale_x_log10()+facet_wrap(~variable,nrow = 1,scales = "free_x")+
    geom_vline(xintercept = 1e4,col="gray")+
    xlab("distance to the Wt1 peak, log10")+ ylab("proportion of events")+
    theme_minimal()+gg_aes+
    theme(legend.position = "bottom")
  
  ggplot(rd, aes(value,col=set)) + stat_ecdf(geom = "point")+ 
    scale_x_log10()+facet_wrap(~anno+variable,nrow = 5)+
    geom_vline(xintercept = 1e4,col="gray")+
    xlab("distance to the Wt1 peak, log10")+ ylab("proportion of events")+
    theme_minimal()+gg_aes+
    theme(legend.position = "bottom")
  
}


### gneomic annotations by the chomatin states
if(F){
  load(paste0(DATADIR,"EditingLevelReport_Anno.RData"))
  
  pl <- c()
  pl <- rbind(pl, 
    cbind(profile="control", data.frame(table(hlpr[hlpr$control>0,]$snpeff_uq,hlpr[hlpr$control>0,]$cHMM_Final))),
    cbind(profile="D7", data.frame(table(hlpr[hlpr$D7>0,]$snpeff_uq,hlpr[hlpr$D7>0,]$cHMM_Final))),
    cbind(profile="D14", data.frame(table(hlpr[hlpr$D14>0,]$snpeff_uq,hlpr[hlpr$D14>0,]$cHMM_Final)))
  )
  pl$profile <- factor(pl$profile,levels=c("control","D7","D14"))
  pl <- pl[pl$Freq>0,]
  names(pl)[2] <- "annotation"
  
  p1 <- ggplot(pl,aes(x=Var2,y=Freq,fill=annotation))+geom_bar(stat = "identity")+
    facet_wrap(~profile)+
    theme(axis.text.x = element_text(angle=90,vjust=0),
          legend.position = "bottom",
          legend.title = element_blank())+
    scale_fill_brewer(palette = "Dark2")+
    scale_y_continuous(expand =c(0,0))+
    xlab("")+ylab("count")
  p2 <- ggplot(pl,aes(x=Var2,y=Freq,fill=annotation))+geom_bar(stat = "identity",position = "fill")+
    facet_wrap(~profile)+
    theme(axis.text.x = element_text(angle=90,vjust=0),
          legend.position = "bottom",
          legend.title = element_blank())+
    scale_fill_brewer(palette = "Dark2")+
    scale_y_continuous(expand =c(0,0))+
    xlab("")+ylab("fraction")
  grid.arrange(p1,p2,nrow=1)
  
  

}









