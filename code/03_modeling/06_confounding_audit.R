#!/usr/bin/env Rscript
# Stage C2 (rigor / confounding audit): is the single-cell-anchored signature an
# independent prognostic factor, or a proxy for WHO grade / IDH / MGMT / tumour purity?
# Input:  results/tcga_training_expr_log2cpm.tsv, results/tcga_training_clinical_mol.tsv,
#         results/signature_pc_coefs.tsv, results/candidate_genes_single_cell.tsv
# Output: results/confounding_audit.tsv, results/multivariable_cox.tsv,
#         results/risk_confounder_cor.tsv, results/clinical_vs_signature_cindex.tsv
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
RES <- "results"

cind <- function(risk, t, e) {
  if (length(unique(risk)) == 1) return(NA)
  n <- length(risk); cnt <- 0; den <- 0
  for (i in seq_len(n)) { if (!e[i]) next; ti <- t[i]
    for (j in seq_len(n)) { if (i == j) next
      if (t[j] > ti || (t[j] == ti && !e[j])) { den <- den + 1
        if (risk[i] > risk[j]) cnt <- cnt + 1 else if (risk[i] == risk[j]) cnt <- cnt + 0.5 } } }
  if (den == 0) NA else cnt / den
}

expr <- read.delim(file.path(RES, "tcga_training_expr_log2cpm.tsv"), row.names = 1, check.names = FALSE)
clin <- read.delim(file.path(RES, "tcga_training_clinical_mol.tsv"), stringsAsFactors = FALSE)
sig  <- read.delim(file.path(RES, "signature_pc_coefs.tsv"), stringsAsFactors = FALSE)
cand <- read.delim(file.path(RES, "candidate_genes_single_cell.tsv"), stringsAsFactors = FALSE)

expr <- expr[, intersect(colnames(expr), clin$sample_submitter_id)]
clin <- clin[match(colnames(expr), clin$sample_submitter_id), ]

# --- frozen risk score (29-gene protein-coding signature) ---
g <- intersect(sig$gene, rownames(expr))
risk <- as.numeric(t(as.matrix(expr[g, ])) %*% sig$coefficient[match(g, sig$gene)])
cat("signature genes used:", length(g), "/", nrow(sig), "\n")

d <- clin
d$risk <- as.numeric(scale(risk))
surv <- Surv(d$OS_time_days, d$OS_event)

# --- (a) risk vs confounders (Spearman) ---
cor_one <- function(x, nm) {
  ok <- !is.na(x)
  ct <- suppressWarnings(cor.test(d$risk[ok], x[ok], method = "spearman"))
  data.frame(confounder = nm, n = sum(ok), rho = unname(ct$estimate), p = ct$p.value)
}
cor_tab <- rbind(
  cor_one(d$age_at_index, "age"),
  cor_one(d$who_grade,    "WHO_grade"),
  cor_one(d$idh_mutant,   "IDH_mutant"),
  cor_one(d$mgmt_meth,    "MGMT_methylated"))
write.table(cor_tab, file.path(RES, "risk_confounder_cor.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# --- (b) purity proxy: PC1 of candidate matrix + leukocyte-marker score ---
cg <- intersect(cand$gene, rownames(expr))
X  <- t(as.matrix(expr[cg, ])); X <- X[, apply(X, 2, sd) > 0]
pc <- prcomp(X, center = TRUE, scale. = TRUE)
d$pc1 <- pc$x[, 1]
lmark <- intersect(c("PTPRC","CD3D","CD3E","CD2","CD14","CD68","CSF1R","AIF1","ITGAM","TYROBP","FCER1G","LYZ"), rownames(expr))
d$leuk <- rowMeans(scale(t(as.matrix(expr[lmark, ]))))
cat("PC1 var-explained:", round(100 * summary(pc)$importance[2, 1], 1), "%; leukocyte markers:", length(lmark), "\n")
purity_tab <- rbind(
  cor_one(d$pc1,  "PC1_candidate_genes"),
  cor_one(d$leuk, "leukocyte_marker_score"))
write.table(purity_tab, file.path(RES, "risk_confounder_cor.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE, append = TRUE, col.names = FALSE)

# --- (c) multivariable Cox (complete cases) ---
cc <- !is.na(d$who_grade) & !is.na(d$idh_mutant) & !is.na(d$mgmt_meth) & !is.na(d$age_at_index)
dd <- d[cc, ]
fit <- coxph(Surv(OS_time_days, OS_event) ~ risk + age_at_index + factor(who_grade) + idh_mutant + mgmt_meth, data = dd)
s <- summary(fit)
mv <- data.frame(term = rownames(s$coefficients), HR = s$conf.int[, "exp(coef)"],
                 lo = s$conf.int[, "lower .95"], hi = s$conf.int[, "upper .95"], p = s$coefficients[, "Pr(>|z|)"])
write.table(mv, file.path(RES, "multivariable_cox.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nMultivariable Cox (n=", nrow(dd), ", events=", sum(dd$OS_event), "):\n", sep = "")
print(mv, row.names = FALSE)

# --- (d) does candidate prognostic signal survive grade+IDH adjustment? ---
mat <- t(as.matrix(expr[cg, ]))
sub <- !is.na(d$who_grade) & !is.na(d$idh_mutant)
adj <- vapply(seq_along(cg), function(j) {
  x <- mat[sub, j]
  if (sd(x) == 0) return(NA_real_)
  f <- tryCatch(coxph(Surv(OS_time_days, OS_event) ~ x + factor(who_grade) + idh_mutant, data = d[sub, ]), error = function(e) NULL)
  if (is.null(f)) NA_real_ else summary(f)$coefficients["x", "Pr(>|z|)"]
}, numeric(1))
uni_raw <- vapply(seq_along(cg), function(j) {
  x <- mat[, j]
  if (sd(x) == 0) return(NA_real_)
  f <- tryCatch(coxph(Surv(OS_time_days, OS_event) ~ x, data = d), error = function(e) NULL)
  if (is.null(f)) NA_real_ else summary(f)$coefficients["x", "Pr(>|z|)"]
}, numeric(1))
hit <- data.frame(n_candidates = length(cg),
                  sig_unadj_p05 = sum(uni_raw < 0.05, na.rm = TRUE),
                  sig_adj_gradeIDH_p05 = sum(adj < 0.05, na.rm = TRUE))
write.table(hit, file.path(RES, "confounding_audit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nCandidate hit-rate unadjusted p<0.05:", hit$sig_unadj_p05, "/", length(cg),
    "  | grade+IDH-adjusted p<0.05:", hit$sig_adj_gradeIDH_p05, "/", length(cg), "\n")

# --- (e) discrimination: clinical vs signature vs combined (complete cases) ---
ci_risk <- cind(dd$risk, dd$OS_time_days, dd$OS_event)
fitc <- coxph(Surv(OS_time_days, OS_event) ~ age_at_index + factor(who_grade) + idh_mutant + mgmt_meth, data = dd)
ci_clin <- cind(predict(fitc), dd$OS_time_days, dd$OS_event)
fitb <- coxph(Surv(OS_time_days, OS_event) ~ risk + age_at_index + factor(who_grade) + idh_mutant + mgmt_meth, data = dd)
ci_both <- cind(predict(fitb), dd$OS_time_days, dd$OS_event)
cmp <- data.frame(model = c("signature_only","clinical_only(age+grade+IDH+MGMT)","signature+clinical"),
                  C_index = round(c(ci_risk, ci_clin, ci_both), 3), n = nrow(dd), events = sum(dd$OS_event))
write.table(cmp, file.path(RES, "clinical_vs_signature_cindex.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nDiscrimination (complete cases):\n"); print(cmp, row.names = FALSE)
cat("DONE\n")
