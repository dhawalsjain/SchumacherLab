rm(list=ls())
setwd("A:/work/scripts/SchumacherLab")
source("vars.R")



rmarkdown::render('R_DE_genes.Rmd',output_dir = REPDIR)


## 2/3 replicate based

rmarkdown::render('2_3/R_Editing_Annotations.Rmd',output_dir = REPDIR)

rmarkdown::render('2_3/R_Editing_Dynamics.Rmd',output_dir = REPDIR)

rmarkdown::render('2_3/R_miRNABindingChange.Rmd',output_dir = REPDIR)

rmarkdown::render('2_3/R_IntronicEdits.Rmd',output_dir = REPDIR)

#rmarkdown::render('2_3/R_RNAstructuresFinal.Rmd',output_dir = REPDIR)

rmarkdown::render('2_3/R_RNAstructures.Rmd',output_dir = REPDIR)

rmarkdown::render('2_3/R_RNAstructures1.Rmd',output_dir = REPDIR)

  rmarkdown::render('2_3/R_RNAstructures2.Rmd',output_dir = REPDIR)

rmarkdown::render('2_3/R_RNAstructures3.Rmd',output_dir = REPDIR)


#rmarkdown::render('2_3/R_Splicing.Rmd',output_dir = REPDIR)

#rmarkdown::render('2_3/R_SpreadSheets.Rmd',output_dir = REPDIR)



## 3/3 replicate based

rmarkdown::render('3_3/R_33_Editing_Annotations.Rmd',output_dir = REPDIR)

rmarkdown::render('3_3/R_33_Editing_Dynamics.Rmd',output_dir = REPDIR)

rmarkdown::render('3_3/R_33_miRNABindingChange.Rmd',output_dir = REPDIR)

rmarkdown::render('3_3/R_33_IntronicEdits.Rmd',output_dir = REPDIR)

rmarkdown::render('3_3/R_33_RNAstructuresFinal.Rmd',output_dir = REPDIR)

rmarkdown::render('3_3/R_33_RNAstructures.Rmd',output_dir = REPDIR)

rmarkdown::render('3_3/R_33_RNAstructures1.Rmd',output_dir = REPDIR)

rmarkdown::render('3_3/R_33_RNAstructures2.Rmd',output_dir = REPDIR)





## 3/3 hc replicate based

rmarkdown::render('3_3_HC/R_33hc_Editing_Annotations.Rmd',output_dir = REPDIR)

rmarkdown::render('3_3_HC/R_33hc_Editing_Dynamics.Rmd',output_dir = REPDIR)

rmarkdown::render('3_3_HC/R_33hc_miRNABindingChange.Rmd',output_dir = REPDIR)

rmarkdown::render('3_3_HC/R_33hc_IntronicEdits.Rmd',output_dir = REPDIR)

rmarkdown::render('3_3_HC/R_33hc_RNAstructuresFinal.Rmd',output_dir = REPDIR)

rmarkdown::render('3_3_HC/R_33hc_RNAstructures.Rmd',output_dir = REPDIR)

rmarkdown::render('3_3_HC/R_33hc_RNAstructures1.Rmd',output_dir = REPDIR)



