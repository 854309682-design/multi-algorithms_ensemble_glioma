#!/usr/bin/env Rscript
# Stage B4: 3rd-cohort surrogate -- repeated patient-level hold-out inside TCGA with the
# FULL pipeline re-fit on each training split (univariate Cox -> LASSO -> freeze CoxPH),
# scored on the untouched test split. This is an honest internal-external validation when
# no additional independent RNA-seq cohort is available on disk (GSE16011 annotation blocked).
# Input:  results/tcga_training_expr_log2cpm.tsv, tcga_training_clinical_mol.tsv,
#         results/candidate_genes_single_cell.tsv ; <REPO_ROOT>/reference/protein_coding_genes.txt
# Output: results/tcga_holdout_validation.tsv
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

suppressMessages({library(survival); library(glmnet)})
set.seed(20260909)
setwd(ROOT)
RES <- "results"; REF <- REF_ROOT

cind <- function(risk, t, e) { n <- length(risk); c <- 0; d <- 0
  for (i in seq_len(n)) { if (!e[i]) next; ti <- t[i]
    for (j in seq_len(n)) { if (i == j) next
      if (t[j] > ti || (t[j] == ti && !e[j])) { d <- d + 1
        if (risk[i] > risk[j]) c <- c + 1 else if (risk[i] == risk[j]) c <- c + 0.5 } } }
  if (d == 0) NA else c / d }

expr <- read.delim(file.path(RES, "tcga_training_expr_log2cpm.tsv"), row.names = 1, check.names = FALSE)
clin <- read.delim(file.path(RES, "tcga_training_clinical_mol.tsv"), stringsAsFactors = FALSE)
cand <- read.delim(file.path(RES, "candidate_genes_single_cell.tsv"), stringsAsFactors = FALSE)
pc <- readLines(file.path(REF, "protein_coding_genes.txt"))
expr <- expr[, intersect(colnames(expr), clin$sample_submitter_id)]
clin <- clin[match(colnames(expr), clin$sample_submitter_id), ]
genes <- intersect(intersect(cand$gene, rownames(expr)), pc)
cat("candidate protein-coding genes:", length(genes), "\n")
Xall <- t(as.matrix(expr[genes, ])); time <- clin$OS_time_days; event <- clin$OS_event

R <- 100; res <- numeric(R)
for (r in 1:R) {
  # stratified 70/30 patient-level split
  tr <- unlist(lapply(split(seq_along(event), event), function(ix) { n <- round(0.7 * length(ix)); sample(ix, n) }))
  te <- setdiff(seq_along(event), tr)
  Xtr <- Xall[tr, , drop = FALSE]; ttr <- time[tr]; etr <- event[tr]
  p <- apply(Xtr, 2, function(x) tryCatch(summary(coxph(Surv(ttr, etr) ~ x))$coefficients[1, "Pr(>|z|)"], error = function(e) 1))
  keep <- which(p < 0.05); if (length(keep) < 5) keep <- order(p)[1:min(30, ncol(Xtr))]
  cv <- tryCatch(cv.glmnet(Xtr[, keep, drop = FALSE], Surv(ttr, etr), family = "cox", alpha = 1, nfolds = 5), error = function(e) NULL)
  if (is.null(cv)) { res[r] <- NA; next }
  cf <- as.matrix(coef(cv, s = "lambda.min")); sel <- rownames(cf)[cf[, 1] != 0]
  if (length(sel) < 2) { res[r] <- NA; next }
  df <- data.frame(time = ttr, event = etr, Xtr[, sel, drop = FALSE], check.names = FALSE)
  f <- tryCatch(coxph(Surv(time, event) ~ ., data = df), error = function(e) NULL)
  if (is.null(f)) { res[r] <- NA; next }
  lp <- predict(f, newdata = data.frame(Xall[te, sel, drop = FALSE], check.names = FALSE))
  res[r] <- cind(lp, time[te], event[te])
  if (r %% 10 == 0) cat("  iter", r, "done\n")
}
res <- res[!is.na(res)]
out <- data.frame(scheme = "TCGA repeated 70/30 holdout (pipeline re-fit)",
                  n_iter = length(res), cindex_mean = round(mean(res), 3), cindex_sd = round(sd(res), 3),
                  q025 = round(quantile(res, 0.025), 3), q500 = round(median(res), 3), q975 = round(quantile(res, 0.975), 3))
write.table(out, file.path(RES, "tcga_holdout_validation.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nDONE -> results/tcga_holdout_validation.tsv\n"); print(out, row.names = FALSE)
