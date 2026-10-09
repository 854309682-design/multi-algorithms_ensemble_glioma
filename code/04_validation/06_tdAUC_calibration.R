#!/usr/bin/env Rscript
# Stage B2: time-dependent AUC (1/3/5 yr) for the frozen signature on the external cohorts,
# plus 1/3/5-yr calibration (training cohort, where the baseline hazard is defined).
# Input:  results/signature_pc_coefs.tsv, signature_pc.rds,
#         tcga_training_expr_log2cpm.tsv, tcga_training_clinical_mol.tsv ; <DATA_ROOT>/CGGA/*.zip
# Output: results/tdAUC_external.tsv, results/calibration_training_135.tsv
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

suppressMessages({library(survival); library(timeROC)})
set.seed(20260909)
setwd(ROOT)
RES <- "results"; D <- file.path(DATA_ROOT, "CGGA")
sig <- read.delim(file.path(RES, "signature_pc_coefs.tsv"), stringsAsFactors = FALSE)
sigmap <- setNames(sig$coefficient, sig$gene)
TIMES <- c(365, 1095, 1825)   # 1 / 3 / 5 years

# ---------- (1) external time-dependent AUC ----------
rows <- list()
for (name in c("693", "325")) {
  tr <- tempfile(fileext = ".txt")
  system(sprintf("unzip -p '%s/CGGA.mRNAseq_%s.RSEM-genes.zip' 'CGGA.mRNAseq_%s.RSEM-genes.20200506.txt' > '%s'", D, name, name, tr))
  ex <- read.delim(tr, sep = "\t", check.names = FALSE); rownames(ex) <- ex[, 1]; ex <- as.matrix(ex[, -1])
  tc <- tempfile(fileext = ".txt")
  system(sprintf("unzip -p '%s/CGGA.mRNAseq_%s_clinical.zip' 'CGGA.mRNAseq_%s_clinical.20200506.txt' > '%s'", D, name, name, tc))
  cl <- read.delim(tc, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE); names(cl)[1] <- "CGGA_ID"
  cl <- cl[trimws(cl$PRS_type) == "Primary", ]
  cl$OS <- suppressWarnings(as.numeric(cl$OS)); cl$Censor <- suppressWarnings(as.numeric(cl$Censor))
  cl <- cl[!is.na(cl$OS) & cl$OS > 0 & !is.na(cl$Censor), ]
  ss <- intersect(colnames(ex), trimws(cl$CGGA_ID)); X <- log2(ex[, ss, drop = FALSE] + 1); cl <- cl[match(ss, trimws(cl$CGGA_ID)), ]
  g <- intersect(names(sigmap), rownames(X)); risk <- as.numeric(t(X[g, , drop = FALSE]) %*% sigmap[g])
  roc <- timeROC(T = cl$OS, delta = cl$Censor, marker = risk, cause = 1, times = TIMES, iid = TRUE)
  auc <- as.numeric(roc$AUC); sdv <- as.numeric(roc$inference$vect_sd_1)
  cat(sprintf("\n== CGGA-%s == n=%d ev=%d | tdAUC\n", name, nrow(cl), sum(cl$Censor)))
  for (i in seq_along(TIMES)) cat(sprintf("  %d-yr AUC=%.3f [%.3f-%.3f]\n", TIMES[i]/365, auc[i], auc[i]-1.96*sdv[i], auc[i]+1.96*sdv[i]))
  rows[[name]] <- data.frame(cohort = paste0("CGGA-", name), horizon_yr = TIMES/365,
                             AUC = round(auc, 3), lo = round(auc - 1.96*sdv, 3), hi = round(auc + 1.96*sdv, 3))
}
write.table(do.call(rbind, rows), file.path(RES, "tdAUC_external.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---------- (2) calibration (training cohort, 1/3/5 yr, risk quartiles) ----------
expr <- read.delim(file.path(RES, "tcga_training_expr_log2cpm.tsv"), row.names = 1, check.names = FALSE)
clin <- read.delim(file.path(RES, "tcga_training_clinical_mol.tsv"), stringsAsFactors = FALSE)
expr <- expr[, intersect(colnames(expr), clin$sample_submitter_id)]
clin <- clin[match(colnames(expr), clin$sample_submitter_id), ]
lp <- as.numeric(t(as.matrix(expr[names(sigmap), ])) %*% sigmap)
f0 <- coxph(Surv(OS_time_days, OS_event) ~ lp, data = clin)
bh <- basehaz(f0, centered = FALSE)
H0 <- approx(bh$time, bh$hazard, xout = TIMES, rule = 2)$y
pred <- exp(-outer(exp(lp), H0))                      # predicted S(t) per sample x horizon
qbin <- cut(lp, quantile(lp, 0:4/4), include.lowest = TRUE, labels = FALSE)
sf <- survfit(Surv(OS_time_days, OS_event) ~ qbin, data = clin)
obs <- summary(sf, times = TIMES, extend = TRUE)$surv        # ordered: strata(qbin) x horizon
# assemble calibration table (predicted mean S(t) per risk quartile vs KM-observed)
predm <- sapply(seq_along(TIMES), function(i) tapply(pred[, i], qbin, mean))   # 4 x 3 (quartile x horizon)
obsm  <- matrix(obs, nrow = 4, ncol = length(TIMES), byrow = TRUE)             # 4 x 3
cal <- data.frame(quartile = rep(1:4, length(TIMES)), horizon_yr = rep(TIMES/365, each = 4),
                  predicted = round(as.numeric(predm), 3), observed = round(as.numeric(obsm), 3))
write.table(cal, file.path(RES, "calibration_training_135.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n== training calibration (predicted vs KM-observed S(t)) ==\n"); print(cal, row.names = FALSE)
cat("\nDONE -> results/tdAUC_external.tsv, results/calibration_training_135.tsv\n")
