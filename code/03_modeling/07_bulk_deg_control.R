#!/usr/bin/env Rscript
# Stage B1 (H1 control): a "naive bulk" signature derived WITHOUT single-cell anchoring,
# built with the identical screen -> LASSO -> freeze pipeline, and applied to the same
# external cohorts, to test whether single-cell anchoring beats bulk-derived features.
# Bulk features = top-300 genes from a limma contrast of survival extremes (short vs long OS).
# NOTE: DESeq2 needs raw counts; the count->sample manifest is no longer on disk, so we use
# limma on the log2-CPM training matrix (equivalent bulk-DEG control; documented in the paper).
# Input:  results/tcga_training_expr_log2cpm.tsv, tcga_training_clinical_mol.tsv,
#         signature_pc_coefs.tsv ; <DATA_ROOT>/CGGA/*.zip
# Output: results/bulk_deg_control_summary.tsv, results/bulk_deg_control_coefs.tsv
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

suppressMessages({library(survival); library(glmnet); library(limma)})
set.seed(20260909)
setwd(ROOT)
RES <- "results"; D <- file.path(DATA_ROOT, "CGGA")

cind <- function(risk, t, e) { ok <- !is.na(risk) & !is.na(t) & !is.na(e); risk <- risk[ok]; t <- t[ok]; e <- e[ok]
  n <- length(risk); if (n < 10) return(NA); c <- 0; d <- 0
  for (i in seq_len(n)) { if (!e[i]) next; ti <- t[i]
    for (j in seq_len(n)) { if (i == j) next
      if (t[j] > ti || (t[j] == ti && !e[j])) { d <- d + 1
        if (risk[i] > risk[j]) c <- c + 1 else if (risk[i] == risk[j]) c <- c + 0.5 } } }
  if (d == 0) NA else c / d }

expr <- read.delim(file.path(RES, "tcga_training_expr_log2cpm.tsv"), row.names = 1, check.names = FALSE)
clin <- read.delim(file.path(RES, "tcga_training_clinical_mol.tsv"), stringsAsFactors = FALSE)
expr <- expr[, intersect(colnames(expr), clin$sample_submitter_id)]
clin <- clin[match(colnames(expr), clin$sample_submitter_id), ]

# --- bulk-DEG control: limma, survival extremes (bottom vs top tertile of OS) ---
q <- quantile(clin$OS_time_days, c(1/3, 2/3))
grp <- ifelse(clin$OS_time_days <= q[1], "Short", ifelse(clin$OS_time_days >= q[2], "Long", NA))
keep <- !is.na(grp); cat("survival-extreme samples: Short", sum(grp == "Short", na.rm = TRUE), "Long", sum(grp == "Long", na.rm = TRUE), "\n")
fit <- eBayes(lmFit(expr[, keep], model.matrix(~ factor(grp[keep]))), trend = TRUE)
tt <- topTable(fit, coef = 2, number = 300, sort.by = "P")
genec <- rownames(tt); cat("bulk-DEG candidates:", length(genec), "\n")

# --- identical pipeline: univariate Cox -> LASSO -> freeze CoxPH ---
y <- Surv(clin$OS_time_days, clin$OS_event)
p <- apply(as.matrix(expr[genec, ]), 1, function(x) tryCatch(summary(coxph(y ~ x))$coefficients[1, "Pr(>|z|)"], error = function(e) 1))
keepg <- genec[p < 0.05]; cat("univariate p<0.05:", length(keepg), "\n")
X <- t(as.matrix(expr[keepg, ]))
cvf <- cv.glmnet(X, y, family = "cox", alpha = 1, nfolds = 10)
co <- as.matrix(coef(cvf, s = "lambda.min")); sel <- rownames(co)[co[, 1] != 0]
cat("lambda.min bulk signature genes:", length(sel), "\n")
df <- data.frame(y = y, t(as.matrix(expr[sel, ])), check.names = FALSE)
f <- coxph(y ~ ., data = df); bic <- coef(f)
cat("bulk in-sample C:", round(cind(predict(f, type = "lp"), clin$OS_time_days, clin$OS_event), 4), "\n")
write.table(data.frame(gene = names(bic), coefficient = bic, HR = exp(bic)),
            file.path(RES, "bulk_deg_control_coefs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# --- apply to external cohorts (frozen), compare with single-cell-anchored signature ---
sigsc <- setNames(read.delim(file.path(RES, "signature_pc_coefs.tsv"))$coefficient,
                  read.delim(file.path(RES, "signature_pc_coefs.tsv"))$gene)
res <- list()
for (name in c("693", "325")) {
  tr <- tempfile(fileext = ".txt")
  system(sprintf("unzip -p '%s/CGGA.mRNAseq_%s.RSEM-genes.zip' 'CGGA.mRNAseq_%s.RSEM-genes.20200506.txt' > '%s'", D, name, name, tr))
  ex <- read.delim(tr, sep = "\t", check.names = FALSE); rownames(ex) <- ex[, 1]; ex <- as.matrix(ex[, -1])
  clf <- sprintf("%s/CGGA.mRNAseq_%s_clinical_clean.tsv", D, name)
  if (file.exists(clf)) {
    cl <- read.delim(clf, stringsAsFactors = FALSE)
  } else {
    tc <- tempfile(fileext = ".txt")
    system(sprintf("unzip -p '%s/CGGA.mRNAseq_%s_clinical.zip' 'CGGA.mRNAseq_%s_clinical.20200506.txt' > '%s'", D, name, name, tc))
    cl <- read.delim(tc, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE); names(cl)[1] <- "CGGA_ID"
  }
  cl <- cl[trimws(cl$PRS_type) == "Primary", ]
  cl$OS <- suppressWarnings(as.numeric(cl$OS)); cl$Censor <- suppressWarnings(as.numeric(cl$Censor))
  cl <- cl[!is.na(cl$OS) & cl$OS > 0 & !is.na(cl$Censor), ]
  ss <- intersect(colnames(ex), trimws(cl$CGGA_ID)); Xe <- log2(ex[, ss, drop = FALSE] + 1); cl <- cl[match(ss, trimws(cl$CGGA_ID)), ]
  gb <- intersect(names(bic), rownames(Xe)); gs <- intersect(names(sigsc), rownames(Xe))
  rb <- as.numeric(t(Xe[gb, , drop = FALSE]) %*% bic[gb])
  rs <- as.numeric(t(Xe[gs, , drop = FALSE]) %*% sigsc[gs])
  cat(sprintf("\nCGGA-%s n=%d ev=%d | bulk C=%.3f (%d/%d genes) | single-cell C=%.3f (%d/%d genes)\n",
              name, nrow(cl), sum(cl$Censor), cind(rb, cl$OS, cl$Censor), length(gb), length(bic), cind(rs, cl$OS, cl$Censor), length(gs), length(sigsc)))
  res[[name]] <- data.frame(cohort = paste0("CGGA-", name), n = nrow(cl), events = sum(cl$Censor),
                            C_bulk_DEG = round(cind(rb, cl$OS, cl$Censor), 3),
                            C_single_cell = round(cind(rs, cl$OS, cl$Censor), 3))
}
out <- do.call(rbind, res)
write.table(out, file.path(RES, "bulk_deg_control_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nDONE -> results/bulk_deg_control_summary.tsv\n"); print(out, row.names = FALSE)
