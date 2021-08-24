install.packages(c("BiocManager", 
                   "devtools",
                   "roxygen2",
                   'rversions'))
devtools::install_version("RcppAnnoy", version="0.0.18", repos="https://cloud.r-project.org")

BiocManager::install(c("Rhdf5lib",
                       "png",
                       "BiocGenerics",
                       "DelayedArray",
                       "DelayedMatrixStats",
                       "limma",
                       "S4Vectors",
                       "SingleCellExperiment",
                       "SummarizedExperiment",
                       "MAST",
                       "batchelor",
                       "chromVAR",
                       "motifmatchr",
                       "DirichletMultinomial",
                       "TFBSTools"),
                       update=TRUE, site_repository="https://cloud.r-project.org")
BiocManager::install("dittoSeq")
install.packages(c("usethis",
                   "reticulate",
                   "units",
                   "grr",
                   "Matrix.utils",
                   "pbmcapply",
                   "pheatmap",
                   "proxy",
                   "pscl",
                   "rsample",
                   "RhpcBLASctl",
                   "slam",
                   "spdep",
                   "speedglm",
                   "RcppProgress",
                   "RcppArmadillo",
                   "Rcpp",
                   "rlang",
                   "tibble",
                   "irlba",
                   "ggplot2",
                   "tidyr",
                   "cowplot",
                   "dplyr",
                   "tidyverse",
                   "ape",
                   "ggpubr"),
                 repos="https://cloud.r-project.org")

devtools::install_github("cole-trapnell-lab/leidenbase", upgrade=FALSE)
devtools::install_github("cole-trapnell-lab/monocle3", ref="1.0.0", upgrade=FALSE)
devtools::install_github("immunogenomics/harmony", ref="0.1", upgrade=FALSE)
install.packages("https://cran.r-project.org/src/contrib/Archive/spatstat/spatstat_1.64-1.tar.gz",
                 repos=NULL,type="source", INSTALL_opts = "--no-lock")
remotes::install_version("Seurat", version = "4.0.3", upgrade="never")
remotes::install_github("mojaveazure/seurat-disk")

# devtools::install_github("GreenleafLab/ArchR", ref="2be1294eb1fbff364fb538d71f4f545ee4384d09", repos = BiocManager::repositories())
# library(ArchR)
# ArchR::installExtraPackages()
