rm(list=ls())
source("functions.R")
source("vars.R")


###################################################################3333
####### 2/3
###################################################################3333

if(F){
  
  myf<- function(path="A:/work/Kreidberg lab/Valerie/PAPER/PSVisulaization/Ref/DRE_Clustered_RefSeqs_Genomic_OUT.fa"){
    ref <- read.delim(path,header=F,stringsAsFactors = F)
    rh <-  data.frame(id=ref[seq(1,1086,6),])
    rh$energy <- ref[seq(3,1086,6),]
    rh$energy <- gsub(".*?\\(","",rh$energy)
    rh$energy <- gsub("\\)","",rh$energy)
    rh$energy <- as.numeric(rh$energy)
    rh$id <- gsub(">","",rh$id)
    rm(ref)
    rh
  }
  
  ref <- myf(path="A:/work/Kreidberg lab/Valerie/PAPER/PSVisulaization/Ref/DRE_Clustered_RefSeqs_Genomic_OUT.fa")
  alt <- myf(path="A:/work/Kreidberg lab/Valerie/PAPER/PSVisulaization/Alt/DRE_Clustered_AltSeqs_Genomic_OUT.fa")
  names(ref)[2] <- "ref_energy"
  names(alt)[2] <- "alt_energy"
  
  rh <- merge(ref,alt,by="id")
  rh$pic <- paste0(gsub("\\|","_",rh$id),"_rel.jpg")
  
  rh$gene <- gsub("\\|\\S*","",rh$id)
  rh$id <- gsub("^\\S*?\\|","",rh$id)
  rh$ncluster <- gsub("\\|\\S*","",rh$id)
  rh$id <- gsub("^\\S*?\\|","",rh$id)
  rh$FDR <- gsub("\\|\\S*","",rh$id)
  rh$id <- gsub("^\\S*?\\|","",rh$id)
  rh$Anno <- gsub("\\|\\S*","",rh$id)
  rh$id <- gsub("^\\S*?\\|","",rh$id)
  rh$medainChange <- gsub("\\|\\S*","",rh$id)
  rh$id <- gsub("^\\S*?\\|","",rh$id)
  rh$comparison <- gsub("\\|\\S*","",rh$id)
  rh$id <- gsub("^\\S*?\\|","",rh$id)
  rh$id <- NULL
  
  rh <- rh[,c("gene", "ncluster", "FDR", "Anno", "medainChange", "comparison","ref_energy", "alt_energy", "pic")]
  rh$picAlt <- paste0("A:/work/Kreidberg lab/Valerie/PAPER/PSVisulaization/Alt/",rh$pic)
  rh$pic <- paste0("A:/work/Kreidberg lab/Valerie/PAPER/PSVisulaization/Ref/",rh$pic)
  
  
  
  
}

if(F){
  myf<- function(path="A:/work/Kreidberg lab/Valerie/PAPER/PSVisulaization/Ref/DRE_Clustered_RefSeqs_Genomic_OUT.fa"){
    ref <- read.delim(path,header=F,stringsAsFactors = F)
    rh <-  data.frame(id=ref[seq(1,nrow(ref),6),])
    rh$energy <- ref[seq(3,nrow(ref),6),]
    rh$energy <- gsub(".*?\\(","",rh$energy)
    rh$energy <- gsub("\\)","",rh$energy)
    rh$energy <- as.numeric(rh$energy)
    rh$id <- gsub(">","",rh$id)
    rm(ref)
    rh
  }
  
  ref <- myf(path="A:/work/Kreidberg lab/Valerie/PAPER/33hc_PSVisulaization/Cntr/33hc_DRE_Clustered_RefSeqs_Genomic_OUT.fa")
  alt <- myf(path="A:/work/Kreidberg lab/Valerie/PAPER/33hc_PSVisulaization/Exp/33hc_DRE_Clustered_AltSeqs_Genomic_OUT.fa")
  names(ref)[2] <- "ref_energy"
  names(alt)[2] <- "alt_energy"
  rh <- merge(ref,alt,by="id")
  rh$pic <- paste0(gsub("\\|","_",rh$id),"_rel.jpg")
  rh$gene <- gsub("\\|\\S*","",rh$id)
  rh$id <- gsub("^\\S*?\\|","",rh$id)
  rh$ncluster <- gsub("\\|\\S*","",rh$id)
  rh$id <- gsub("^\\S*?\\|","",rh$id)
  rh$FDR <- gsub("\\|\\S*","",rh$id)
  rh$id <- gsub("^\\S*?\\|","",rh$id)
  rh$Anno <- gsub("\\|\\S*","",rh$id)
  rh$id <- gsub("^\\S*?\\|","",rh$id)
  rh$medainChange <- gsub("\\|\\S*","",rh$id)
  rh$id <- gsub("^\\S*?\\|","",rh$id)
  rh$comparison <- gsub("\\|\\S*","",rh$id)
  rh$id <- gsub("^\\S*?\\|","",rh$id)
  rh$id <- NULL
  rh <- rh[,c("gene", "ncluster", "FDR", "Anno", "medainChange", "comparison","ref_energy", "alt_energy", "pic")]
  rhx <- rh
  
  load(paste0(DATADIR,"33hc_DRE_DE_matrix.RData"))
  r1 <- res[res$comparison=="D7_control" & res$clustered>0 & !is.na(res$tracking_id) & res$adj.P.Val<0.05,]
  r1 <- r1[r1$D7>r1$control,]
  rh1 <- rh[rh$comparison=="D7_control" & rh$gene%in%r1$gene,]
  r1 <- res[res$comparison=="D14_control" & res$clustered>0 & !is.na(res$tracking_id) & res$adj.P.Val<0.05,]
  r1 <- r1[r1$D14>r1$control,]
  rh2 <- rh[rh$comparison=="D14_control" & rh$gene%in%r1$gene,]
  rh <- rbind(rh1,rh2)
  rm(rh1,rh2,r1)  
  
  rh <- rh[rh$medainChange>0,c("gene","comparison","ref_energy","alt_energy")]
  #rh <- rh[rh$ref_energy>rh$alt_energy,]
  rh$comparison <-gsub("_control","",rh$comparison)
  rh <- reshape2::dcast(rh,gene+ref_energy ~ comparison, value.var = "alt_energy")
  names(rh)[2:4] <- c("Control","Day-14","Day-7")
  
  rh1 <- rh[,c(1,2,4)]
  rh1 <- rh1[complete.cases(rh1),]
  rh1[,3] <- ((rh1[,3]-rh1[,2])*100)/abs(rh1[,2])
  rh1[,2] <- 0
  rh1 <- melt(rh1,measure.vars = names(rh1)[2:3])
  p4 <- ggpaired(rh1, x = "variable", y = "value",color = "variable", 
                 add = "jitter",line.color = "gray", line.size = 0.4)+gg_aes+
    scale_color_manual(values=cols)+scale_fill_manual(values=cols)+
    geom_hline(yintercept = 0,col="red",linetype="dotted",size=1.2)+
    stat_compare_means(label = "p.signif", method = "t.test",paired = T,ref.group = "Control")+
    xlab("")+ylab("% change in energy\n[ (Control-D7)/Control ]")
  
  rh1 <- rh[,c(1,2,3)]
  rh1 <- rh1[complete.cases(rh1),]
  rh1[,3] <- ((rh1[,3]-rh1[,2])*100)/abs(rh1[,2])
  rh1[,2] <- 0
  rh1 <- melt(rh1,measure.vars = names(rh1)[2:3])
  p5 <- ggpaired(rh1, x = "variable", y = "value",color = "variable", 
            add = "jitter",line.color = "gray", line.size = 0.4)+gg_aes+
    scale_color_manual(values=cols)+scale_fill_manual(values=cols)+
    geom_hline(yintercept = 0,col="red",linetype="dotted",size=1.2)+
    stat_compare_means(label = "p.signif", method = "t.test",paired = T,ref.group = "Control")+
    xlab("")+ylab("% change in energy\n[ (Control-D14)/Control ]")
  
  p1 <- ggdraw() + draw_image(paste0("A:/work/Kreidberg lab/Valerie/PAPER/33hc_PSVisulaization/Cntr/",rhx$pic[50]),scale = 1)+
    ggtitle(paste0("Control\n",rhx$gene[50],"\n(",rhx$ref_energy[50],"  kcal/mol)" ) )+
    theme(plot.title = element_text(hjust = 0.5))
  p2 <- ggdraw() + draw_image(paste0("A:/work/Kreidberg lab/Valerie/PAPER/33hc_PSVisulaization/Exp/",rhx$pic[50]),scale = 1)+  
    ggtitle(paste0("Day-7\n",rhx$gene[50],"\n(",rhx$alt_energy[50],"  kcal/mol)" ) )+
    theme(plot.title = element_text(hjust = 0.5))
  p3 <- ggdraw() + draw_image(paste0("A:/work/Kreidberg lab/Valerie/PAPER/33hc_PSVisulaization/Cntr/",rhx$pic[50]),scale = 1)+  
    ggtitle(paste0("Day-14\n",rhx$gene[50],"\n(",rhx$ref_energy[50],"  kcal/mol)" ) )+
    theme(plot.title = element_text(hjust = 0.5))
  
  lay <- rbind(c(1,1,2,2,3,3),
               c(4,4,4,5,5,5))
  grid.arrange(p1,p2,p3,p4,p5,heights=c(5,4), layout_matrix = lay)
  
  
  
  
  
}



###################################################################3333
####### 3/3
###################################################################3333

if(F){
  
  myf<- function(path="A:/work/Kreidberg lab/Valerie/PAPER/PSVisulaization/Ref/DRE_Clustered_RefSeqs_Genomic_OUT.fa"){
    ref <- read.delim(path,header=F,stringsAsFactors = F)
    rh <-  data.frame(id=ref[seq(1,nrow(ref),6),])
    rh$energy <- ref[seq(3,nrow(ref),6),]
    rh$energy <- gsub(".*?\\(","",rh$energy)
    rh$energy <- gsub("\\)","",rh$energy)
    rh$energy <- as.numeric(rh$energy)
    rh$id <- gsub(">","",rh$id)
    rm(ref)
    rh
  }
  
  ref <- myf(path="A:/work/Kreidberg lab/Valerie/PAPER/33_PSVisulaization/Ref/33_DRE_Clustered_RefSeqs_Genomic_OUT.fa")
  alt <- myf(path="A:/work/Kreidberg lab/Valerie/PAPER/33_PSVisulaization/Alt/33_DRE_Clustered_AltSeqs_Genomic_OUT.fa")
  names(ref)[2] <- "ref_energy"
  names(alt)[2] <- "alt_energy"
  rh <- merge(ref,alt,by="id")
  rh$pic <- paste0(gsub("\\|","_",rh$id),"_rel.jpg")
  rh$gene <- gsub("\\|\\S*","",rh$id)
  rh$id <- gsub("^\\S*?\\|","",rh$id)
  rh$ncluster <- gsub("\\|\\S*","",rh$id)
  rh$id <- gsub("^\\S*?\\|","",rh$id)
  rh$FDR <- gsub("\\|\\S*","",rh$id)
  rh$id <- gsub("^\\S*?\\|","",rh$id)
  rh$Anno <- gsub("\\|\\S*","",rh$id)
  rh$id <- gsub("^\\S*?\\|","",rh$id)
  rh$medainChange <- gsub("\\|\\S*","",rh$id)
  rh$id <- gsub("^\\S*?\\|","",rh$id)
  rh$comparison <- gsub("\\|\\S*","",rh$id)
  rh$id <- gsub("^\\S*?\\|","",rh$id)
  rh$id <- NULL
  rh <- rh[,c("gene", "ncluster", "FDR", "Anno", "medainChange", "comparison","ref_energy", "alt_energy", "pic")]
  
  rh$picAlt <- paste0("A:/work/Kreidberg lab/Valerie/PAPER/33_PSVisulaization/Alt/",rh$pic)
  rh$pic <- paste0("A:/work/Kreidberg lab/Valerie/PAPER/33_PSVisulaization/Ref/",rh$pic)
  
  
  
  
}

if(F){
  myf<- function(path="A:/work/Kreidberg lab/Valerie/PAPER/PSVisulaization/Ref/DRE_Clustered_RefSeqs_Genomic_OUT.fa"){
    ref <- read.delim(path,header=F,stringsAsFactors = F)
    rh <-  data.frame(id=ref[seq(1,nrow(ref),6),])
    rh$energy <- ref[seq(3,nrow(ref),6),]
    rh$energy <- gsub(".*?\\(","",rh$energy)
    rh$energy <- gsub("\\)","",rh$energy)
    rh$energy <- as.numeric(rh$energy)
    rh$id <- gsub(">","",rh$id)
    rm(ref)
    rh
  }
  
  ref <- myf(path="A:/work/Kreidberg lab/Valerie/PAPER/33_PSVisulaization/Ref/33_DRE_Clustered_RefSeqs_Genomic_OUT.fa")
  alt <- myf(path="A:/work/Kreidberg lab/Valerie/PAPER/33_PSVisulaization/Alt/33_DRE_Clustered_AltSeqs_Genomic_OUT.fa")
  names(ref)[2] <- "ref_energy"
  names(alt)[2] <- "alt_energy"
  rh <- merge(ref,alt,by="id")
  rh$pic <- paste0(gsub("\\|","_",rh$id),"_rel.jpg")
  rh$gene <- gsub("\\|\\S*","",rh$id)
  rh$id <- gsub("^\\S*?\\|","",rh$id)
  rh$ncluster <- gsub("\\|\\S*","",rh$id)
  rh$id <- gsub("^\\S*?\\|","",rh$id)
  rh$FDR <- gsub("\\|\\S*","",rh$id)
  rh$id <- gsub("^\\S*?\\|","",rh$id)
  rh$Anno <- gsub("\\|\\S*","",rh$id)
  rh$id <- gsub("^\\S*?\\|","",rh$id)
  rh$medainChange <- gsub("\\|\\S*","",rh$id)
  rh$id <- gsub("^\\S*?\\|","",rh$id)
  rh$comparison <- gsub("\\|\\S*","",rh$id)
  rh$id <- gsub("^\\S*?\\|","",rh$id)
  rh$id <- NULL
  
  rh <- rh[,c("gene", "ncluster", "FDR", "Anno", "medainChange", "comparison","ref_energy", "alt_energy", "pic")]
  rhx <- rh
  
  load(paste0(DATADIR,"33_DRE_DE_matrix.RData"))
  r1 <- res[res$comparison=="D7_control" & res$clustered>0 & !is.na(res$tracking_id),]
  r1 <- r1[r1$D7<r1$control,]
  rh1 <- rh[rh$comparison=="D7_control" & rh$gene%in%r1$gene,]
  r1 <- res[res$comparison=="D14_control" & res$clustered>0 & !is.na(res$tracking_id),]
  r1 <- r1[r1$D14<r1$control,]
  rh2 <- rh[rh$comparison=="D14_control" & rh$gene%in%r1$gene,]
  rh <- rbind(rh1,rh2)
  rm(rh1,rh2,r1)  
  
  rh <- rh[,c("gene","comparison","ref_energy","alt_energy")]
  #rh <- rh[rh$ref_energy>rh$alt_energy,]
  rh$comparison <-gsub("_control","",rh$comparison)
  rh <- reshape2::dcast(rh,gene+ref_energy ~ comparison, value.var = "alt_energy")
  names(rh)[2] <- "control"
  
  rh1 <- rh[,c(1,2,4)]
  rh1 <- rh1[complete.cases(rh1),]
  rh1[,3] <- ((rh1[,3]-rh1[,2])*100)/abs(rh1[,2])
  rh1[,2] <- 0
  rh1 <- melt(rh1,measure.vars = names(rh1)[2:3])
  p4 <- ggpaired(rh1, x = "variable", y = "value",color = "variable", 
                 add = "jitter",line.color = "gray", line.size = 0.4)+gg_aes+
    scale_color_manual(values=cols)+scale_fill_manual(values=cols)+
    geom_hline(yintercept = 0,col="red",linetype="dotted",size=1.2)+
    stat_compare_means(label = "p.signif", method = "t.test",paired = T,ref.group = "control")+
    xlab("")+ylab("% change in energy\n[ (Control-D7)/Control ]")
  
  rh1 <- rh[,c(1,2,3)]
  rh1 <- rh1[complete.cases(rh1),]
  rh1[,3] <- ((rh1[,3]-rh1[,2])*100)/abs(rh1[,2])
  rh1[,2] <- 0
  rh1 <- melt(rh1,measure.vars = names(rh1)[2:3])
  p5 <- ggpaired(rh1, x = "variable", y = "value",color = "variable", 
                 add = "jitter",line.color = "gray", line.size = 0.4)+gg_aes+
    scale_color_manual(values=cols)+scale_fill_manual(values=cols)+
    geom_hline(yintercept = 0,col="red",linetype="dotted",size=1.2)+
    stat_compare_means(label = "p.signif", method = "t.test",paired = T,ref.group = "control")+
    xlab("")+ylab("% change in energy\n[ (Control-D14)/Control ]")
  
  p1 <- ggdraw() + draw_image(paste0("A:/work/Kreidberg lab/Valerie/PAPER/PSVisulaization/Ref/",rhx$pic[71]),scale = 1)+
    ggtitle("control\nPolrf3")+theme(plot.title = element_text(hjust = 0.5))
  p2 <- ggdraw() + draw_image(paste0("A:/work/Kreidberg lab/Valerie/PAPER/PSVisulaization/Alt/",rhx$pic[72]),scale = 1)+  
    ggtitle("Day-7\nPolrf3")+theme(plot.title = element_text(hjust = 0.5))
  p3 <- ggdraw() + draw_image(paste0("A:/work/Kreidberg lab/Valerie/PAPER/PSVisulaization/Alt/",rhx$pic[71]),scale = 1)+  
    ggtitle("Day-14\nPolrf3")+theme(plot.title = element_text(hjust = 0.5))
  
  lay <- rbind(c(1,1,2,2,3,3),
               c(4,4,4,5,5,5))
  grid.arrange(p1,p2,p3,p4,p5,heights=c(5,4), layout_matrix = lay)
  
  
  
  
  
}






