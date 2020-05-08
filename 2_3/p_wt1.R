rm(list = ls())
source("functions.R")
source("vars.R")


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
  
  myplfunc <- function(r,x,y,tit){
    q <- ggplot(r, aes(r[,x], log10(r[,y])) ) +geom_point() +
      stat_bin_hex(geom="point", aes(color=log2(..count..)), binwidth=c(.05,.05), size=1, shape=20, stroke=0)+
      scale_color_viridis(option="inferno",direction = 1)+
      geom_smooth(method = "lm") +gg_aes+
      ggtitle(tit)+
      xlab(x)+ylab("distance to Wt1 peak in Control\n (log10)")+
      stat_cor(method = "pearson", size=5,col="red",label.sep = "\n",
               label.x.npc = "right", label.y.npc = "top")
    q
  }
  
  
  pllist <- list()
  r <- res[res$comparison=="D7_control",]
  r <- merge(r,hlpr[,c("idx","Control_Wt1","D9_Wt1","D14_Wt1" )],by="idx")
  pllist[[1]] <- myplfunc(r,"change","Control_Wt1","Control-D7")
  r <- res[res$comparison=="D14_control",]
  r <- merge(r,hlpr[,c("idx","Control_Wt1","D9_Wt1","D14_Wt1" )],by="idx")
  pllist[[2]] <- myplfunc(r,"change","Control_Wt1","Control-D14")
  
  grid.arrange(grobs=pllist,nrow=1)
  

  
}
