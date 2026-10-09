## R packages required for this repository (R 4.6.1 used for the manuscript).
## Install:  Rscript code/00_setup/R_packages.R
## The analysis machine used https://cloud.r-project.org; some country mirrors returned
## HTTP 403, so the upstream CRAN mirror is pinned here.

repos <- "https://cloud.r-project.org"

cran_reqs <- c(
  "survival", "survminer", "glmnet", "rms", "dcurves", "timeROC",
  "randomForestSRC", "CoxBoost", "survivalsvm", "xgboost", "c060", "stabsel",
  "data.table", "readr", "openxlsx", "ggplot2", "ggpubr", "pheatmap",
  "metafor", "mice", "VIM", "sva", "edgeR", "DESeq2", "circlize", "Matrix", "limma"
)

bioc_reqs <- c(
  "org.Hs.eg.db", "hgu133plus2.db", "clusterProfiler", "GSVA", "enrichplot",
  "SummarizedExperiment", "SingleCellExperiment", "BiocManager", "GEOquery",
  "biomaRt", "TCGAbiolinks", "MuSiC", "CellChat", "infercnv", "copykat"
)

missing_cran <- setdiff(cran_reqs, rownames(installed.packages()))
if (length(missing_cran)) install.packages(missing_cran, repos = repos)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = repos)
}
missing_bioc <- setdiff(bioc_reqs, rownames(installed.packages()))
if (length(missing_bioc)) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

## Notes
##  * `infercnv` needs JAGS and `copykat` a separate build step; neither was available on the
##    analysis machine, so the CNV check in stage 1 uses the infercnv support score already stored
##    in the processed single-cell object (see code/01_sc_anchor/07_cnv_malignant_qc.py).
##  * `CellChat` and `MuSiC` are used by stages 08_signalling and 07_state_composition.
##  * sessionInfo() of the analysis run: code/00_setup/sessionInfo_R.txt
