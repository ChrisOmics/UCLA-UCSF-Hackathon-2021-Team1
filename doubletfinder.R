
#Doublet Finder on total dataset
library(modes)
library(DoubletFinder)

for (sample in levels(as.factor(mgAVM02$sample_description))){
  LD1Ctrl = subset(mgAVM02,subset = sample_description == "Control_P5")
  sweep.res.listC <- paramSweep_v3(LD1Ctrl, PCs = 1:30, sct = FALSE)
  sweep.statsC <- summarizeSweep(sweep.res.listC, GT = FALSE)
  bcmvnC <- find.pK(sweep.statsC)
  homotypic.prop <- modelHomotypic(mgAVM02@meta.data$finalclusters) 
  nExp_poi <- round(0.075*nrow(mgAVM02@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  mgAVM02 = doubletFinder_v3(mgAVM02,PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj,reuse.pANN = F, sct = F)
  
}



setwd("~/Desktop")
