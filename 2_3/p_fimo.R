source("functions.R")
source("vars.R")


###############################################################################################
### FIMO motif analysis
###############################################################################################
if(F){
  fimo <- c()
  cntr <- read.delim(paste0(DATADIR,"meme/control_fimo/fimo.txt"),header=T,stringsAsFactors = F)
  cntr <- cntr[cntr$q.value<0.01 & cntr$strand=="+",]
  cntr <- unique(cntr[,1:2])
  fimo <- rbind(fimo,cbind(data.frame(table(cntr$X.pattern.name)),total=nrow(cntr),profile="control"))
  
  cntr <- read.delim(paste0(DATADIR,"meme/D7_fimo/fimo.txt"),header=T,stringsAsFactors = F)
  cntr <- cntr[cntr$q.value<0.01 & cntr$strand=="+",]
  cntr <- unique(cntr[,1:2])
  fimo <- rbind(fimo,cbind(data.frame(table(cntr$X.pattern.name)),total=nrow(cntr),profile="D7"))
  
  cntr <- read.delim(paste0(DATADIR,"meme/D14_fimo/fimo.txt"),header=T,stringsAsFactors = F)
  cntr <- cntr[cntr$q.value<0.01 & cntr$strand=="+",]
  cntr <- unique(cntr[,1:2])
  fimo <- rbind(fimo,cbind(data.frame(table(cntr$X.pattern.name)),total=nrow(cntr),profile="D14"))
  
  fimo$perc <- round(fimo$Freq*100/fimo$total  ,2)
  fimo  
  
  ggplot(fimo,aes(x=Var1,y=perc, fill=profile,alpha=0.2))+geom_bar(stat="identity",position = "dodge")+
    scale_fill_manual(values = cols)+
    xlab("")+ylab("% of total sites")+theme_bw()+gg_aes+
    theme(legend.position = "top")
  
  
}



###############################################################################################
### FIMO motif analysis _33
###############################################################################################
if(F){
  fimo <- c()
  cntr <- read.delim(paste0(DATADIR,"meme_33/control_fimo/fimo.txt"),header=T,stringsAsFactors = F)
  cntr <- cntr[cntr$q.value<0.01 & cntr$strand=="+",]
  cntr <- unique(cntr[,1:2])
  fimo <- rbind(fimo,cbind(data.frame(table(cntr$X.pattern.name)),total=nrow(cntr),profile="control"))
  
  cntr <- read.delim(paste0(DATADIR,"meme_33/D7_fimo/fimo.txt"),header=T,stringsAsFactors = F)
  cntr <- cntr[cntr$q.value<0.01 & cntr$strand=="+",]
  cntr <- unique(cntr[,1:2])
  fimo <- rbind(fimo,cbind(data.frame(table(cntr$X.pattern.name)),total=nrow(cntr),profile="D7"))
  
  cntr <- read.delim(paste0(DATADIR,"meme_33/D14_fimo/fimo.txt"),header=T,stringsAsFactors = F)
  cntr <- cntr[cntr$q.value<0.01 & cntr$strand=="+",]
  cntr <- unique(cntr[,1:2])
  fimo <- rbind(fimo,cbind(data.frame(table(cntr$X.pattern.name)),total=nrow(cntr),profile="D14"))
  
  fimo$perc <- round(fimo$Freq*100/fimo$total  ,2)
  fimo  
  
  ggplot(fimo,aes(x=Var1,y=perc, fill=profile,alpha=0.2))+geom_bar(stat="identity",position = "dodge")+
    scale_fill_manual(values = cols)+
    xlab("")+ylab("% of total sites")+theme_bw()+gg_aes+
    theme(legend.position = "top")
  
  
}