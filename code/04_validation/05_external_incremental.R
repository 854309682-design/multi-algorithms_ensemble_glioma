#!/usr/bin/env Rscript
# Stage C3 (rigor / decisive test): does the frozen single-cell-anchored signature add
# prognostic value BEYOND age + grade + IDH + 1p19q + MGMT in the external cohorts?
# Reports risk-only vs clinical-only vs risk+clinical C-index and the adjusted HR for risk.
# Input:  results/signature_pc_coefs.tsv ; <DATA_ROOT>/CGGA/CGGA.mRNAseq_{693,325}*.zip
# Output: results/external_incremental_value.tsv
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

suppressMessages({library(survival)})
setwd(ROOT)
RES <- "results"; D <- file.path(DATA_ROOT, "CGGA")

cind <- function(risk, t, e) { ok <- !is.na(risk) & !is.na(t) & !is.na(e); risk <- risk[ok]; t <- t[ok]; e <- e[ok]
  n <- length(risk); if (n < 10) return(NA); c <- 0; d <- 0
  for (i in seq_len(n)) { if (!e[i]) next; ti <- t[i]
    for (j in seq_len(n)) { if (i == j) next
      if (t[j] > ti || (t[j] == ti && !e[j])) { d <- d + 1
        if (risk[i] > risk[j]) c <- c + 1 else if (risk[i] == risk[j]) c <- c + 0.5 } } }
  if (d == 0) NA else c / d }

sig <- read.delim(file.path(RES, "signature_pc_coefs.tsv"), stringsAsFactors = FALSE)
sigmap <- setNames(sig$coefficient, sig$gene)

clean <- function(x) trimws(as.character(x))
mlab <- function(x, pos, neg) ifelse(clean(x) %in% pos, 1L, ifelse(clean(x) %in% neg, 0L, NA_integer_))

row_out <- list()
for (name in c("693", "325")) {
  zr <- file.path(D, paste0("CGGA.mRNAseq_", name, ".RSEM-genes.zip"))
  zc <- file.path(D, paste0("CGGA.mRNAseq_", name, "_clinical.zip"))
  tr <- tempfile(fileext = ".txt"); tc <- tempfile(fileext = ".txt")
  system(sprintf("unzip -p '%s' 'CGGA.mRNAseq_%s.RSEM-genes.20200506.txt' > '%s'", zr, name, tr))
  system(sprintf("unzip -p '%s' 'CGGA.mRNAseq_%s_clinical.20200506.txt' > '%s'", zc, name, tc))
  expr <- read.delim(tr, sep = "\t", check.names = FALSE); rownames(expr) <- expr[, 1]; expr <- as.matrix(expr[, -1])
  clin <- read.delim(tc, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  names(clin)[1] <- "CGGA_ID"

  clin <- clin[clean(clin$PRS_type) == "Primary", ]
  clin$OS <- suppressWarnings(as.numeric(clin$OS)); clin$Censor <- suppressWarnings(as.numeric(clin$Censor))
  clin <- clin[!is.na(clin$OS) & clin$OS > 0 & !is.na(clin$Censor), ]
  samples <- intersect(colnames(expr), clean(clin$CGGA_ID))
  X <- log2(expr[, samples, drop = FALSE] + 1); clin <- clin[match(samples, clean(clin$CGGA_ID)), ]

  genes <- intersect(names(sigmap), rownames(X))
  risk <- as.numeric(t(X[genes, , drop = FALSE]) %*% sigmap[genes])
  d <- data.frame(OS = clin$OS, ev = clin$Censor, risk = as.numeric(scale(risk)),
                  age = suppressWarnings(as.numeric(clean(clin$Age))),
                  grade = { g <- clean(clin$Grade)
                            ifelse(grepl("IV", g), 4L, ifelse(grepl("III", g), 3L, ifelse(grepl("II\\b", g), 2L, NA_integer_))) },
                  idh = mlab(clin$IDH_mutation_status, "Mutant", "Wildtype"),
                  codel = mlab(clin[["1p19q_codeletion_status"]], "Codel", "Non-codel"),
                  mgmt = mlab(clin[["MGMTp_methylation_status"]], c("methylated", "Methylated"), c("un-methylated", "Un-methylated")))
  d$grade[!d$grade %in% c(2, 3, 4)] <- NA

  ci_risk <- cind(d$risk, d$OS, d$ev)
  cc <- !is.na(d$age) & !is.na(d$grade) & !is.na(d$idh) & !is.na(d$codel) & !is.na(d$mgmt)
  dd <- d[cc, ]
  f_clin <- coxph(Surv(OS, ev) ~ age + factor(grade) + idh + codel + mgmt, data = dd)
  f_both <- coxph(Surv(OS, ev) ~ risk + age + factor(grade) + idh + codel + mgmt, data = dd)
  ci_clin <- cind(predict(f_clin), dd$OS, dd$ev)
  ci_both <- cind(predict(f_both), dd$OS, dd$ev)
  lrt <- anova(f_clin, f_both, test = "LRT")$`Pr(>|Chi|)`[2]
  s <- summary(f_both); hr <- s$conf.int["risk", ]; pv <- s$coefficients["risk", "Pr(>|z|)"]
  cat(sprintf("\n== CGGA-%s == n=%d ev=%d | sig genes %d/%d\n", name, nrow(d), sum(d$ev), length(genes), length(sigmap)))
  cat(sprintf("  C: risk-only %.3f | clinical-only %.3f | risk+clinical %.3f\n", ci_risk, ci_clin, ci_both))
  cat(sprintf("  adjusted risk HR=%.2f [%.2f-%.2f] p=%.2g ; added-value LRT p=%.2g\n", hr["exp(coef)"], hr["lower .95"], hr["upper .95"], pv, lrt))
  row_out[[name]] <- data.frame(cohort = paste0("CGGA-", name), n = nrow(d), events = sum(d$ev),
                                 C_risk = round(ci_risk, 3), C_clinical = round(ci_clin, 3), C_risk_clinical = round(ci_both, 3),
                                 adj_risk_HR = round(hr["exp(coef)"], 2), adj_risk_lo = round(hr["lower .95"], 2), adj_risk_hi = round(hr["upper .95"], 2),
                                 adj_risk_p = signif(pv, 3), added_value_LRT_p = signif(lrt, 3))
}
out <- do.call(rbind, row_out)
write.table(out, file.path(RES, "external_incremental_value.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nDONE -> results/external_incremental_value.tsv\n"); print(out, row.names = FALSE)
