#!/usr/bin/env Rscript
# Stage 4 (external validation): apply the FROZEN CoxPH signature to CGGA-693 & CGGA-325.
# Risk score = sum(coef_g * expr_g) on log2(RSEM+1) expression, primary tumors with OS.
# C-index (manual Harrell), KM (within-cohort median high/low), log-rank.
## --- portable project roots: do not edit (see config/paths.R) ---------------
.start_dir <- tryCatch({
  .a <- commandArgs(trailingOnly = FALSE)
  dirname(normalizePath(sub("^--file=", "", .a[grep("^--file=", .a)][1])))
}, error = function(e) getwd(), warning = function(w) getwd())
.find_root <- function(d) {
  repeat {
    if (file.exists(file.path(d, "config", "paths.R"))) return(d)
    .p <- dirname(d)
    if (identical(.p, d)) stop("repository root not found -- set GLIOMA_ROOT")
    d <- .p
  }
}
REPO_ROOT <- if (nzchar(Sys.getenv("GLIOMA_ROOT"))) normalizePath(Sys.getenv("GLIOMA_ROOT")) else .find_root(.start_dir)
source(file.path(REPO_ROOT, "config", "paths.R"))
## ----------------------------------------------------------------------------

suppressMessages({library(survival); library(data.table)})
set.seed(20260909)
setwd(ROOT)
RES <- "results"; D <- file.path(DATA_ROOT, "CGGA")

# frozen signature
sig <- read.delim(file.path(RES, "signature_genes_coefs.tsv"), stringsAsFactors=FALSE)
sigmap <- setNames(sig$coefficient, sig$gene)
cat("frozen signature genes:", length(sigmap), "\n")

cind <- function(risk, t, e){ n<-length(risk); c<-0; d<-0
  for(i in 1:n){ if(!e[i]) next; ti<-t[i]
    for(j in 1:n){ if(i==j) next
      if(t[j]>ti || (t[j]==ti && !e[j])){ d<-d+1; if(risk[i]>risk[j]) c<-c+1 else if(risk[i]==risk[j]) c<-c+0.5 } } }
  if(d==0) NA else c/d }

for (name in c("693","325")) {
  cat("\n===== CGGA-", name, " =====\n", sep="")
  # extract the zip's .txt to a temp file (zip also holds __MACOSX resource entries)
  zf <- paste0(D, "/CGGA.mRNAseq_", name, ".RSEM-genes.20200506.txt.zip")
  tmp <- tempfile(fileext=".txt")
  system(paste0("unzip -p '", zf, "' 'CGGA.mRNAseq_", name, ".RSEM-genes.20200506.txt' > '", tmp, "'"))
  expr <- read.delim(tmp, sep="\t", check.names=FALSE)
  rownames(expr) <- expr[,1]; expr <- as.matrix(expr[,-1])
  clin <- read.delim(paste0(D, "/CGGA.mRNAseq_", name, "_clinical_clean.tsv"), stringsAsFactors=FALSE)
  # primary tumors with OS
  prim <- clin[clin$PRS_type == "Primary", ]
  prim$OS <- as.numeric(prim$OS); prim$Censor <- as.numeric(prim$Censor)
  keep <- !is.na(prim$OS) & prim$OS > 0
  clinp <- prim[keep, ]
  samples <- intersect(colnames(expr), clinp$CGGA_ID)
  expr <- expr[, samples]; clinp <- clinp[match(samples, clinp$CGGA_ID), ]
  cat("primary samples with OS:", nrow(clinp), "; events:", sum(clinp$Censor), "\n")

  # log2 transform
  X <- log2(as.matrix(expr) + 1)
  genes <- intersect(names(sigmap), rownames(X))
  cat("signature genes present:", length(genes), "/", length(sigmap), "\n")
  risk <- as.numeric(t(X[genes, , drop=FALSE]) %*% sigmap[genes])
  time <- clinp$OS; event <- clinp$Censor
  cat("C-index:", round(cind(risk, time, event), 4), "\n")

  # KM by within-cohort median
  med <- median(risk)
  grp <- ifelse(risk > med, "High", "Low")
  sfit <- survfit(Surv(time, event) ~ grp)
  lr <- survdiff(Surv(time, event) ~ grp)
  pval <- 1 - pchisq(lr$chisq, df=1)
  cat("median cutoff; High n=", sum(grp=="High"), " Low n=", sum(grp=="Low"),
      " log-rank p=", format(pval, digits=3), "\n", sep="")
  writeLines(paste("CGGA-", name, "\t", round(cind(risk,time,event),4), "\t", sum(grp=="High"),
                   sum(grp=="Low"), "\t", format(pval,digits=4), sep=""),
             paste0(RES, "/external_cgga_", name, "_result.txt"))
}
cat("DONE\n")
