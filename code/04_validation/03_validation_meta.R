#!/usr/bin/env Rscript
# Stage 4 (external validation summary + meta-analysis):
# Compute per-cohort C-index (bootstrap CI) for the RNA-seq cohorts and a random-effects
# meta-C-index (metafor). Consolidates the honest RNA-seq-centric validation framing.
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

suppressMessages({library(survival); library(metafor); library(data.table)})
set.seed(20260909)
setwd(ROOT)
RES <- "results"; D <- file.path(DATA_ROOT, "CGGA")

sig <- read.delim(file.path(RES, "signature_genes_coefs.tsv"), stringsAsFactors=FALSE)
sm <- setNames(sig$coefficient, sig$gene)

cind <- function(risk, t, e){ n<-length(risk); c<-0; d<-0
  for(i in 1:n){ if(!e[i]) next; ti<-t[i]
    for(j in 1:n){ if(i==j) next
      if(t[j]>ti || (t[j]==ti && !e[j])){ d<-d+1; if(risk[i]>risk[j]) c<-c+1 else if(risk[i]==risk[j]) c<-c+0.5 } } }
  if(d==0) NA else c/d }
boot_ci <- function(risk, t, e, B=300){
  ci <- numeric(B); n <- length(risk)
  for(b in 1:B){ idx<-sample(n, replace=TRUE); ci[b]<-cind(risk[idx], t[idx], e[idx]) }
  c(mean=mean(ci,na.rm=TRUE), lo=unname(quantile(ci,0.025,na.rm=TRUE)),
    hi=unname(quantile(ci,0.975,na.rm=TRUE)), sd=sd(ci,na.rm=TRUE))
}

rows <- list()
for (name in c("693","325")) {
  zf <- paste0(D, "/CGGA.mRNAseq_", name, ".RSEM-genes.20200506.txt.zip")
  tmp <- tempfile(fileext=".txt")
  system(paste0("unzip -p '", zf, "' 'CGGA.mRNAseq_", name, ".RSEM-genes.20200506.txt' > '", tmp, "'"))
  ex <- read.delim(tmp, sep="\t", check.names=FALSE); rownames(ex) <- ex[,1]; ex <- as.matrix(ex[,-1])
  clin <- read.delim(paste0(D, "/CGGA.mRNAseq_", name, "_clinical_clean.tsv"), stringsAsFactors=FALSE)
  prim <- clin[clin$PRS_type=="Primary", ]; prim$OS <- as.numeric(prim$OS); prim$Censor <- as.numeric(prim$Censor)
  keep <- !is.na(prim$OS) & prim$OS > 0; clinp <- prim[keep, ]
  samples <- intersect(colnames(ex), clinp$CGGA_ID); ex <- ex[, samples]; clinp <- clinp[match(samples, clinp$CGGA_ID), ]
  X <- log2(as.matrix(ex)+1); g <- intersect(names(sm), rownames(X))
  risk <- as.numeric(t(X[g,,drop=FALSE]) %*% sm[g])
  bc <- boot_ci(risk, clinp$OS, clinp$Censor)
  rows[[name]] <- data.frame(cohort=paste0("CGGA-",name), n=length(risk), events=sum(clinp$Censor),
                             genes=length(g), C_index=bc[["mean"]], lo=bc[["lo"]], hi=bc[["hi"]], se=bc[["sd"]])
  cat(sprintf("CGGA-%s: C=%.3f [%.3f-%.3f] n=%d\n", name, bc[["mean"]], bc[["lo"]], bc[["hi"]], length(risk)))
}
df <- do.call(rbind, rows)
# random-effects meta-analysis (pool C-index inverse-variance)
rma_fit <- tryCatch(rma(yi=df$C_index, sei=df$se, method="REML"), error=function(e) NULL)
if (!is.null(rma_fit)) {
  pooled <- data.frame(cohort="POOLED (REML)", n=sum(df$n), events=sum(df$events),
                       genes=NA, C_index=rma_fit$b[1], lo=rma_fit$ci.lb, hi=rma_fit$ci.ub,
                       se=rma_fit$se)
  df <- rbind(df, pooled)
  cat(sprintf("\nPOOLED RNA-seq C-index = %.3f [%.3f-%.3f], I2=%s\n",
              rma_fit$b[1], rma_fit$ci.lb, rma_fit$ci.ub, round(rma_fit$I2,1)))
}
write.table(df, file.path(RES, "external_validation_summary.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cat("Saved external_validation_summary.tsv\n")
cat("DONE\n")
