rm(list=ls())
source("functions.R")
source("vars.R")

############################################################################################################
#### produces genes data frame
############################################################################################################

TXFILE="A:/work/Kreidberg lab/Valerie/Bam0_CuffDiff/isoform_exp.diff" ## isoform level diff expression by cuffdiff
GENEFILE="A:/work/Kreidberg lab/Valerie/Bam0_CuffDiff/gene_exp.diff" ## gene level diff expression by cuffdiff
UCSC="A:/work/mm10_Annotations/Ensembl2ucsc.txt"
ANNO="A:/work/mm10_Annotations/mm10_Annotations.bed"
ENSEMBL = "A:/work/mm10_Annotations/EnsemblGeneLoci.txt.txt"

if(!file.exists(paste0(DATADIR,"genes.RData"))){
  
  
  ens2ucsc <- read.delim(UCSC,header=T,stringsAsFactors = F)
  names(ens2ucsc)[c(1,6)] <- c("ensembl","ucsc")
  ens2ucsc$ensembl <- gsub("[.]\\S*","",ens2ucsc$ensembl)
  ens2ucsc$chrom <- gsub("chr","",ens2ucsc$chrom)
  ens2ucsc <- unique(ens2ucsc)
  
  tx <- read.delim(TXFILE,header=T,stringsAsFactors = F)
  tx$chr <- gsub(":\\S*","",tx$locus) %>% as.character
  tx$start <-  gsub("-\\S*","",gsub("\\S*:","",tx$locus)) %>% as.numeric
  tx$end <-  gsub("\\S*-","",tx$locus) %>% as.numeric
  tx <- tx[,c("test_id","gene_id","gene")]
  tx <- unique(tx)
  tx <- merge(tx,ens2ucsc, by.x="test_id",by.y='ensembl')

  expn <- read.delim(GENEFILE,header=T,stringsAsFactors = F)
  expn$geneChr <- gsub(":\\S*","",expn$locus) %>% as.character
  expn$geneStart <-  gsub("-\\S*","",gsub("\\S*:","",expn$locus)) %>% as.numeric
  expn$geneEnd <-  gsub("\\S*-","",expn$locus) %>% as.numeric
  expn <- unique(expn[,c("geneChr","geneStart","geneEnd","gene_id")])
  genes <- merge(expn,tx,by="gene_id")
  genes <- genes[genes$geneChr==genes$chrom,]
  
  genes <- genes[,c("chrom", "strand", "gene","geneStart", "geneEnd","gene_id","test_id" ,"txStart", "txEnd", "ucsc")]
  names(genes)[c(1,7,10)] <- c("chr","txID","ucscID")
  gGenes <- genes
  
  
  anno <- read.delim(ANNO,header=F,stringsAsFactors = F)
  anno <- anno[,c(1,4)] %>% unique
  anno$gene_id <- gsub("\\|\\S*","",anno$V4)
  anno$V4 <- gsub("^\\S*?\\|","",anno$V4)
  anno$name <-  gsub("\\|\\S*","",anno$V4)
  anno$V4 <- gsub("^\\S*?\\|","",anno$V4)
  anno$biotype <-  gsub("\\|\\S*","",anno$V4)
  anno <- unique(anno[,c("gene_id","biotype")])
  gGenes <- merge(gGenes,anno,by="gene_id",all.x=T)
  gGenes <- gGenes[,c(2:6,1,7:11)] 
  save(gGenes,file=paste0(DATADIR,"genes.RData"))
  
}



if(F){
  g2go <- read.delim("A:/work/mm10_Annotations/ensembl2GO.txt",header=T)
  names(g2go) <- c("gene_id","GO")
  g2go <- g2go[g2go$GO!="",]
  g2go <- data.table(g2go)
  g2go <- g2go[,GOAnno:=paste0(GO,collapse = ","),by=list(gene_id)]
  g2go <- as.data.frame(g2go)
  g2go <- unique(g2go[,c(1,3)])
  save(g2go,file = paste0(RAWDIR,"ensembl_g2go.RData"))
  
  
}



