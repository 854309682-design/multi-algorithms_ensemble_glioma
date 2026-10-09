#!/usr/bin/env Rscript
# Stage 3/4 sensitivity analyses (Fig 7):
#   (1) LASSO lambda sensitivity -- lambda.min vs lambda.1se (same cv.glmnet fit);
#   (2) training-cohort sensitivity -- refit the whole pipeline on TCGA GBM ONLY (drop LGG);
#   (3) validation-subset sensitivity -- the frozen primary signature on CGGA GBM-only samples.
# Each variant is frozen (CoxPH) and scored in-sample (TCGA) and on the two external RNA-seq
# cohorts (CGGA-693 / CGGA-325). Scoring = frozen coef x log2(RSEM+1), matching the PRIMARY external
# validation (no per-gene z-scoring), so the lambda.min variant reproduces the headline C-indexes.
# The lambda.min/all-TCGA variant is cross-checked against the primary results/signature_pc_coefs.tsv.
# Outputs: results/sensitivity_signature_variants.tsv, results/sensitivity_external.tsv,
#          results/sensitivity_external_gbmonly.tsv, results/figures/fig7_sensitivity.{pdf,jpg}
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
RES <- "results"; FIG <- "results/figures"
D <- file.path(DATA_ROOT, "CGGA")
REF <- REF_ROOT

cind <- function(risk, t, e) { ok <- !is.na(risk) & !is.na(t) & !is.na(e); risk <- risk[ok]; t <- t[ok]; e <- e[ok]
  n <- length(risk); if (n < 10) return(NA); c <- 0; d <- 0
  for (i in seq_len(n)) { if (!e[i]) next; ti <- t[i]
    for (j in seq_len(n)) { if (i == j) next
      if (t[j] > ti || (t[j] == ti && !e[j])) { d <- d + 1
        if (risk[i] > risk[j]) c <- c + 1 else if (risk[i] == risk[j]) c <- c + 0.5 } } }
  if (d == 0) NA else c / d }

# ---- training inputs -------------------------------------------------------------------------
expr <- read.delim(file.path(RES, "tcga_training_expr_log2cpm.tsv"), row.names = 1, check.names = FALSE)
clin <- read.delim(file.path(RES, "tcga_training_clinical.tsv"), stringsAsFactors = FALSE)
cand <- read.delim(file.path(RES, "candidate_genes_single_cell.tsv"), stringsAsFactors = FALSE)
pc <- readLines(file.path(REF, "protein_coding_genes.txt"))
expr <- expr[, intersect(colnames(expr), clin$sample_submitter_id)]
clin <- clin[match(colnames(expr), clin$sample_submitter_id), ]
cand_pc <- intersect(intersect(cand$gene, rownames(expr)), pc)
cat("protein-coding candidate genes:", length(cand_pc), "\n")

# ---- one pipeline run: univariate Cox -> cv.glmnet -> freeze CoxPH at a given lambda ---------
fit_variant <- function(keep_samples, lambda) {
  E <- expr[, keep_samples, drop = FALSE]; C <- clin[match(keep_samples, clin$sample_submitter_id), ]
  y <- Surv(C$OS_time_days, C$OS_event)
  mat <- t(as.matrix(E[cand_pc, ])); res <- list()
  for (g in cand_pc) { x <- mat[, g]; if (sd(x) == 0) next
    f <- tryCatch(coxph(y ~ x), error = function(e) NULL); if (is.null(f)) next
    res[[g]] <- summary(f)$coefficients[1, "Pr(>|z|)"] }
  p <- unlist(res); keepg <- names(p)[p < 0.05]
  X <- t(as.matrix(E[keepg, ]))
  set.seed(20260909); cvf <- cv.glmnet(X, y, family = "cox", alpha = 1, nfolds = 10)
  sel <- rownames(as.matrix(coef(cvf, s = lambda))); sel <- sel[as.matrix(coef(cvf, s = lambda))[, 1] != 0]
  sel <- intersect(sel, pc)
  df <- data.frame(y = y, t(as.matrix(E[sel, ])), check.names = FALSE)
  f <- coxph(y ~ ., data = df)
  list(genes = sel, coef = coef(f), insample = cind(predict(f, type = "lp"), C$OS_time_days, C$OS_event),
       n = nrow(df), events = sum(C$OS_event), n_uni = length(keepg))
}

all_s <- clin$sample_submitter_id
gbm_s <- clin$sample_submitter_id[clin$project == "TCGA-GBM"]
cat("samples: all =", length(all_s), "| GBM-only =", length(gbm_s), "\n")

v_all_min <- fit_variant(all_s, "lambda.min")     # primary (expect 29 genes)
v_all_1se <- fit_variant(all_s, "lambda.1se")     # lambda.1se
v_gbm_min <- fit_variant(gbm_s, "lambda.min")     # GBM-only training

# ---- sanity check vs the frozen primary signature -------------------------------------------
prim <- read.delim(file.path(RES, "signature_pc_coefs.tsv"), stringsAsFactors = FALSE)
cat("primary 29-gene reproduction:", length(intersect(prim$gene, v_all_min$genes)), "/",
    length(prim$gene), "genes identical\n")

# ---- external cohorts (same loader as 07) ----------------------------------------------------
load_cohort <- function(name) {
  tr <- tempfile(); tc <- tempfile()
  system(sprintf("unzip -p '%s/CGGA.mRNAseq_%s.RSEM-genes.zip' 'CGGA.mRNAseq_%s.RSEM-genes.20200506.txt' > '%s'", D, name, name, tr))
  system(sprintf("unzip -p '%s/CGGA.mRNAseq_%s_clinical.zip' 'CGGA.mRNAseq_%s_clinical.20200506.txt' > '%s'", D, name, name, tc))
  ex <- read.delim(tr, sep = "\t", check.names = FALSE); rownames(ex) <- ex[, 1]; ex <- as.matrix(ex[, -1])
  cl <- read.delim(tc, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE); names(cl)[1] <- "CGGA_ID"
  cl <- cl[trimws(cl$PRS_type) == "Primary", ]
  cl$OS <- suppressWarnings(as.numeric(cl$OS)); cl$Censor <- suppressWarnings(as.numeric(cl$Censor))
  cl <- cl[!is.na(cl$OS) & cl$OS > 0 & !is.na(cl$Censor), ]
  ss <- intersect(colnames(ex), trimws(cl$CGGA_ID)); ex <- log2(ex[, ss, drop = FALSE] + 1)
  cl <- cl[match(ss, trimws(cl$CGGA_ID)), ]
  g <- trimws(cl$Grade)
  cl$grade <- ifelse(grepl("IV", g), 4L, ifelse(grepl("III", g), 3L, ifelse(grepl("II\\b", g), 2L, NA_integer_)))
  list(X = ex, cl = cl)
}
coh <- list(`CGGA-693` = load_cohort("693"), `CGGA-325` = load_cohort("325"))

score <- function(sig, X) { g <- intersect(names(sig), rownames(X))
  list(risk = as.numeric(t(X[g, , drop = FALSE]) %*% sig[g]), used = length(g)) }

variants <- list(`lambda.min (all TCGA)` = v_all_min, `lambda.1se (all TCGA)` = v_all_1se,
                 `lambda.min (GBM only)` = v_gbm_min)

# table A -- signature variants
A <- do.call(rbind, lapply(names(variants), function(v) { s <- variants[[v]]
  data.frame(variant = v, n_genes = length(s$genes), n_uni_p05 = s$n_uni,
             training_n = s$n, training_events = s$events, in_sample_C = round(s$insample, 3)) }))
write.table(A, file.path(RES, "sensitivity_signature_variants.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# table B -- external C per variant
B <- do.call(rbind, lapply(names(variants), function(v) do.call(rbind, lapply(names(coh), function(cc) {
  sc <- score(variants[[v]]$coef, coh[[cc]]$X)
  data.frame(variant = v, cohort = cc, genes_used = sc$used,
             C_index = round(cind(sc$risk, coh[[cc]]$cl$OS, coh[[cc]]$cl$Censor), 3)) }))))
write.table(B, file.path(RES, "sensitivity_external.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# table C -- primary signature, external ALL vs GBM-only (grade IV)
Ctb <- do.call(rbind, lapply(names(coh), function(cc) { X <- coh[[cc]]$X; cl <- coh[[cc]]$cl
  sc <- score(v_all_min$coef, X); ok <- !is.na(cl$grade) & cl$grade == 4
  rbind(data.frame(cohort = cc, subset = "all", n = nrow(cl), events = sum(cl$Censor),
                   C_index = round(cind(sc$risk, cl$OS, cl$Censor), 3)),
        data.frame(cohort = cc, subset = "GBM only (grade IV)", n = sum(ok), events = sum(cl$Censor[ok]),
                   C_index = round(cind(sc$risk[ok], cl$OS[ok], cl$Censor[ok]), 3))) }))
write.table(Ctb, file.path(RES, "sensitivity_external_gbmonly.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nA)\n"); print(A, row.names = FALSE)
cat("\nB)\n"); print(B, row.names = FALSE)
cat("\nC)\n"); print(Ctb, row.names = FALSE)

# ---- Fig 7 -----------------------------------------------------------------------------------
draw <- function() {
  par(mfrow = c(1, 2), mar = c(6, 4, 3, 1))
  m <- matrix(c(A$in_sample_C, B$C_index), nrow = 3, byrow = FALSE,
              dimnames = list(c("Training", "CGGA-693", "CGGA-325"), A$variant))
  barplot(m, beside = TRUE, ylim = c(0, 1), col = c("grey25", "grey60", "firebrick"),
          ylab = "C-index", main = "(a)  Signature variants",
          legend.text = rownames(m), args.legend = list(x = "topright", bty = "n", cex = 0.8))
  abline(h = 0.5, lty = 3, col = "grey40", lwd = 0.8)
  m2 <- matrix(Ctb$C_index, nrow = 2)
  barplot(m2, beside = TRUE, ylim = c(0, 1), col = c("grey50", "firebrick"),
          names.arg = unique(Ctb$cohort), ylab = "C-index",
          main = "(b)  Primary 29-gene: all vs GBM-only",
          legend.text = unique(Ctb$subset), args.legend = list(x = "bottomright", bty = "n", cex = 0.8))
  abline(h = 0.5, lty = 3, col = "grey40", lwd = 0.8)
}
pdf(file.path(FIG, "fig7_sensitivity.pdf"), width = 10, height = 5); draw(); dev.off()
jpeg(file.path(FIG, "fig7_sensitivity.jpg"), width = 2000, height = 1000, res = 200); draw(); dev.off()
cat("\nDONE -> sensitivity_signature_variants.tsv / sensitivity_external.tsv / sensitivity_external_gbmonly.tsv + fig7_sensitivity.*\n")
