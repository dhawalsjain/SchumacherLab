rm(list = ls())
source("functions.R")
source("vars.R")

########################################################################################
## this script plots intermutation distances for high confidence sites 


if(F){
  load(paste0(DATADIR,"EditingSites_Kataegis.RData"))
  cf$known <- ifelse(cf$known==0, as.character("novel"),as.character("known"))
  
  rect1=rbind(data.frame(xmin=0,xmax=Inf,ymin=0,ymax=log10(50),loc="clustered"))
  
  p <- ggplot(cf,aes(y=log10(dist2next),x=pos,col=repeats)) + geom_point(alpha=0.4,size=1.5)+  
         geom_vline(xintercept = c(chrs$cumlen,chrs$cumlen[24]+chrs$length[24]),col="black") + xlab("") + ylab("Genomic Distance") +
         scale_y_continuous(expand = c(0, 0),breaks=c(0,2,4,6),labels=c("0", "0.1k", "10k","1Mio"))+
         scale_x_continuous(expand = c(0, 0),breaks=(chrs$cumlen+round(chrs$length/2)),labels=factor(chrs$chr) )+ #+ theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
         facet_wrap(~profile,ncol=1,scales = "free_y")+
         scale_fill_brewer(palette = "Dark2")+scale_color_brewer(palette = "Dark2")+
         geom_rect(data=rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="pink", alpha=0.3, inherit.aes = FALSE)+
         gg_aes+     
         theme(legend.position = "bottom",legend.key = element_rect(fill=NA),legend.title = element_blank(),  
               legend.text = element_text(colour="black", size = 15, face = "bold"),
               strip.background = element_blank(),strip.text = element_text(size=15,colour = "black"),
               plot.title = element_text(size = 16,hjust=0.5))+ 
         guides(colour = guide_legend(override.aes = list(alpha=1,size=5)))+
         ggtitle(label = paste("Editing sites"),subtitle = "(color by repeat annotations)")
  p
  
    
  p <- ggplot(cf,aes(y=log10(dist2next),x=pos,col=snpeff_uq)) + geom_point(alpha=0.4,size=1.5)+  
    geom_vline(xintercept = c(chrs$cumlen,chrs$cumlen[24]+chrs$length[24]),col="black") + xlab("") + ylab("Genomic Distance") +
    scale_y_continuous(expand = c(0, 0),breaks=c(0,2,4,6),labels=c("0", "0.1k", "10k","1Mio"))+
    scale_x_continuous(expand = c(0, 0),breaks=(chrs$cumlen+round(chrs$length/2)),labels=factor(chrs$chr) )+ #+ theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
    facet_wrap(~profile,ncol=1,scales = "free_y")+
    scale_fill_brewer(palette = "Dark2")+scale_color_brewer(palette = "Dark2")+
    geom_rect(data=rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="pink", alpha=0.3, inherit.aes = FALSE)+
    gg_aes+     
    theme(legend.position = "bottom",legend.key = element_rect(fill=NA),legend.title = element_blank(),  
          legend.text = element_text(colour="black", size = 15, face = "bold"),
          strip.background = element_blank(),strip.text = element_text(size=15,colour = "black"),
          plot.title = element_text(size = 16,hjust=0.5))+ 
    guides(colour = guide_legend(override.aes = list(alpha=1,size=5)))+
    ggtitle(label = paste("Editing sites"),subtitle = "(color by genomic annotations)") + theme(plot.title = element_text(size = 16,hjust=0.5))
  p
  
  
  p <- ggplot(cf,aes(y=log10(dist2next),x=pos,col=as.factor(known))) + geom_point(alpha=0.4,size=1.5)+  
    geom_vline(xintercept = c(chrs$cumlen,chrs$cumlen[24]+chrs$length[24]),col="black") + xlab("") + ylab("Genomic Distance") +
    scale_y_continuous(expand = c(0, 0),breaks=c(0,2,4,6),labels=c("0", "0.1k", "10k","1Mio"))+
    scale_x_continuous(expand = c(0, 0),breaks=(chrs$cumlen+round(chrs$length/2)),labels=factor(chrs$chr) )+ #+ theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
    facet_wrap(~profile,ncol=1,scales = "free_y")+
    scale_fill_brewer(palette = "Dark2")+scale_color_brewer(palette = "Dark2")+
    geom_rect(data=rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="pink", alpha=0.3, inherit.aes = FALSE)+
    gg_aes+     
    theme(legend.position = "bottom",legend.key = element_rect(fill=NA),legend.title = element_blank(),  
          legend.text = element_text(colour="black", size = 15, face = "bold"),
          strip.background = element_blank(),strip.text = element_text(size=15,colour = "black"),
          plot.title = element_text(size = 16,hjust=0.5))+ 
    guides(colour = guide_legend(override.aes = list(alpha=1,size=5)))+
    ggtitle(label = paste("Editing sites"),subtitle = "(color by novelty)") + theme(plot.title = element_text(size = 16,hjust=0.5))
  p

  cnt$col="1"
  cnt <- melt(cnt,measure.vars = c("D7_FPKM","D14_FPKM","Control_FPKM"))
  p <- ggplot(cnt,aes(y=value,x=pos,col=col)) + geom_point(alpha=0.4,size=0.8)+  
    geom_vline(xintercept = c(chrs$cumlen,chrs$cumlen[24]+chrs$length[24]),col="black") + xlab("") + ylab("Expression, log2 FPKM") +
    scale_y_continuous(expand = c(0, 0))+
    scale_x_continuous(expand = c(0, 0),breaks=(chrs$cumlen+round(chrs$length/2)),labels=factor(chrs$chr) )+ #+ theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
    facet_wrap(~variable,ncol=1,scales = "free_y")+
    scale_fill_brewer(palette = "Dark2")+scale_color_brewer(palette = "Dark2")+
    geom_rect(data=rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="pink", alpha=0.3, inherit.aes = FALSE)+
    gg_aes+     
    theme(legend.position = "bottom",legend.key = element_rect(fill=NA),legend.title = element_blank(),  
          legend.text = element_text(colour="black", size = 15, face = "bold"),
          strip.background = element_blank(),strip.text = element_text(size=15,colour = "black"),
          plot.title = element_text(size = 16,hjust=0.5))+ 
    guides(colour = guide_legend(override.aes = list(alpha=1,size=5)))+
    ggtitle(label = paste("Gene expression"),subtitle = "(log2 FPKM)") + theme(plot.title = element_text(size = 16,hjust=0.5))
  p
  

  cntLevels$col="1"
  cntLevels <- melt(cntLevels,measure.vars = c("D7","D14","Control"))
  p <- ggplot(cntLevels,aes(y=value,x=pos,col=col)) + geom_point(alpha=0.4,size=0.8)+  
    geom_vline(xintercept = c(chrs$cumlen,chrs$cumlen[24]+chrs$length[24]),col="black") + xlab("") + ylab("Editing level") +
    scale_y_continuous(expand = c(0, 0))+
    scale_x_continuous(expand = c(0, 0),breaks=(chrs$cumlen+round(chrs$length/2)),labels=factor(chrs$chr) )+ #+ theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
    facet_wrap(~variable,ncol=1,scales = "free_y")+
    scale_fill_brewer(palette = "Dark2")+scale_color_brewer(palette = "Dark2")+
    geom_rect(data=rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="pink", alpha=0.3, inherit.aes = FALSE)+
    gg_aes+ylim(0,1)+     
    theme(legend.position = "bottom",legend.key = element_rect(fill=NA),legend.title = element_blank(),  
          legend.text = element_text(colour="black", size = 15, face = "bold"),
          strip.background = element_blank(),strip.text = element_text(size=15,colour = "black"),
          plot.title = element_text(size = 16,hjust=0.5))+ 
    guides(colour = guide_legend(override.aes = list(alpha=1,size=5)))+
    ggtitle(label = paste("RNA editing level"),subtitle = "") + theme(plot.title = element_text(size = 16,hjust=0.5))
  p
  
    
}












