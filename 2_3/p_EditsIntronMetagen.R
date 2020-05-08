
rm(list=ls())
source("functions.R")
source("vars.R")

######################33
if(F){
  myf_plot <- function(pllist){
    gb <- list()
    for(f in names(pllist)){
      rect1=rbind(data.frame(xmin=100,xmax=Inf,ymin=0,ymax=Inf,loc="start"),
                  data.frame(xmin=0,xmax=100,ymin=0,ymax=Inf,loc="end"))
      p <- ggplot(pllist[[f]],aes(x=pos,y=mean,ymin=mean-se,ymax=mean+se,col=profile))+
        geom_point(size=4,alpha=0.5)+geom_linerange()+
        ggtitle(f)+
        facet_wrap(~loc,nrow = 1,scales = "free_x")+
        geom_rect(data=rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="gray", alpha=0.3, inherit.aes = FALSE)+
        geom_vline(xintercept = 101,col="black")+
        ylab(paste0("row-scaled\n mean editing"))+xlab("")+
        scale_x_continuous(breaks = c(1,50,100,150,200),labels = c("-1kb","","feature","","+1kb"),expand = c(0,0))+
        scale_color_manual(values = cols)+theme_bw()+gg_aes+
        theme(strip.background = element_blank(),
              strip.placement =  "inside",
              legend.position = "bottom")
      gb[[f]] <- print(p)
    }
    gb
  }
  
  load(paste0(DATADIR,"EditingMetagenePlots_IntronExonUTRs.RData"))
  gb <- myf_plot(pllist)
  grid.arrange(grobs=gb,nrow=2)
  
  load(paste0(DATADIR,"EditingMetagenePlots_IntronExonUTRs_DBCalls.RData"))
  gb <- myf_plot(pllist)
  grid.arrange(grobs=gb,nrow=2)
  
  load(paste0(DATADIR,"EditingMetagenePlots_IntronExonUTRs_D7vsControl.RData"))
  gb <- myf_plot(pllist)
  grid.arrange(grobs=gb,nrow=2)
  
  
  load(paste0(DATADIR,"EditingMetagenePlots_IntronExonUTRs_D14vsControl.RData"))
  gb <- myf_plot(pllist)
  grid.arrange(grobs=gb,nrow=2)
  
  
}

