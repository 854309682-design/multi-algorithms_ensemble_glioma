#!/usr/bin/env Rscript
# P1.4 - Does a state-composition model (bulk-deconvolved fractions) predict survival, and does
# it add value beyond clinical / the gene signature? Fit on TCGA, apply FROZEN to CGGA-693/325.
# Output: results/state_composition/state_model_performance.tsv, state_confounder_cor.tsv
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
set.seed(20260909)
setwd(ROOT)
RES <- "results"; SC <- file.path(RES, "state_composition"); D <- file.path(DATA_ROOT, "CGGA")

cind <- function(risk, t, e) { ok <- !is.na(risk) & !is.na(t) & !is.na(e); risk <- risk[ok]; t <- t[ok]; e <- e[ok]
  n <- length(risk); if (n < 10) return(NA); c <- 0; d <- 0
  for (i in seq_len(n)) { if (!e[i]) next; ti <- t[i]
    for (j in seq_len(n)) { if (i == j) next
      if (t[j] > ti || (t[j] == ti && !e[j])) { d <- d + 1
        if (risk[i] > risk[j]) c <- c + 1 else if (risk[i] == risk[j]) c <- c + 0.5 } } }
  if (d == 0) NA else c / d }

# states (drop near-empty Oligodendrocyte)
fr_t <- read.delim(file.path(SC, "state_fractions_tcga.tsv"), row.names = 1, check.names = FALSE)
st <- setdiff(colnames(fr_t), "Oligodendrocyte")
fr_t <- fr_t[, st]
mol <- read.delim(file.path(RES, "tcga_training_clinical_mol.tsv"), stringsAsFactors = FALSE)
mol <- mol[match(rownames(fr_t), mol$sample_submitter_id), ]
mol <- cbind(mol, fr_t)
num <- cbind(mol, fr_t)
dat <- data.frame(time = mol$OS_time_days, event = mol$OS_event, fr_t, check.names = FALSE)

# --- TCGA: fraction model (in-sample + patient-level 10-fold CV) ---
fmod <- coxph(Surv(time, event) ~ ., data = dat)
ci_in <- cind(predict(fmod), dat$time, dat$event)
fold <- integer(nrow(dat)); o <- sample(seq_len(nrow(dat))); for (k in 1:10) fold[o[seq(k, length(o), by = 10)]] <- k
cv <- sapply(1:10, function(k) { tr <- fold != k
  f <- coxph(Surv(time, event) ~ ., data = dat[tr, , drop = FALSE])
  cind(predict(f, newdata = dat[!tr, , drop = FALSE]), dat$time[!tr], dat$event[!tr]) })
cat("TCGA state-fraction model: in-sample C =", round(ci_in, 3), "| 10-fold CV C =", round(mean(cv, na.rm = TRUE), 3), "\n")

# --- confounder correlations (TCGA) ---
cc <- function(x, nm) { ok <- !is.na(x); r <- suppressWarnings(cor.test(predict(fmod)[ok], x[ok], method = "spearman"))
  data.frame(fraction_model_lp = nm, rho = round(unname(r$estimate), 3), p = signif(r$p.value, 3)) }
cor <- rbind(cc(mol$who_grade, "WHO_grade"), cc(mol$idh_mutant, "IDH_mutant"), cc(mol$mgmt_meth, "MGMT_methylated"),
             cc(mol$age_at_index, "age"))
write.table(cor, file.path(SC, "state_confounder_cor.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(cor, row.names = FALSE)

# --- frozen transfer to CGGA + clinical comparison ---
sig <- read.delim(file.path(RES, "signature_pc_coefs.tsv")); sigmap <- setNames(sig$coefficient, sig$gene)
tcga <- read.delim(file.path(RES, "tcga_training_expr_log2cpm.tsv"), row.names = 1, check.names = FALSE)
rows <- list()
for (name in c("693", "325")) {
  fr <- read.delim(file.path(SC, paste0("state_fractions_cgga", name, ".tsv")), row.names = 1, check.names = FALSE)[, st]
  tc <- tempfile(fileext = ".txt")
  system(sprintf("unzip -p '%s/CGGA.mRNAseq_%s_clinical.zip' 'CGGA.mRNAseq_%s_clinical.20200506.txt' > '%s'", D, name, name, tc))
  cl <- read.delim(tc, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE); names(cl)[1] <- "CGGA_ID"
  tr <- tempfile(fileext = ".txt")
  system(sprintf("unzip -p '%s/CGGA.mRNAseq_%s.RSEM-genes.zip' 'CGGA.mRNAseq_%s.RSEM-genes.20200506.txt' > '%s'", D, name, name, tr))
  ex <- read.delim(tr, sep = "\t", check.names = FALSE); rownames(ex) <- ex[, 1]; ex <- as.matrix(ex[, -1])
  cl <- cl[trimws(cl$PRS_type) == "Primary", ]; cl$OS <- suppressWarnings(as.numeric(cl$OS)); cl$Censor <- suppressWarnings(as.numeric(cl$Censor))
  cl <- cl[!is.na(cl$OS) & cl$OS > 0 & !is.na(cl$Censor), ]
  ss <- intersect(rownames(fr), trimws(cl$CGGA_ID)); fr <- fr[ss, , drop = FALSE]; cl <- cl[match(ss, trimws(cl$CGGA_ID)), ]
  clin <- data.frame(OS = cl$OS, ev = cl$Censor,
                     age = suppressWarnings(as.numeric(trimws(cl$Age))),
                     grade = { g <- trimws(cl$Grade); ifelse(grepl("IV", g), 4L, ifelse(grepl("III", g), 3L, ifelse(grepl("II\\b", g), 2L, NA))) },
                     idh = ifelse(trimws(cl$IDH_mutation_status) == "Mutant", 1L, ifelse(trimws(cl$IDH_mutation_status) == "Wildtype", 0L, NA)),
                     codel = ifelse(trimws(cl[["1p19q_codeletion_status"]]) == "Codel", 1L, ifelse(trimws(cl[["1p19q_codeletion_status"]]) == "Non-codel", 0L, NA)),
                     mgmt = ifelse(trimws(cl[["MGMTp_methylation_status"]]) %in% c("methylated","Methylated"), 1L,
                                   ifelse(trimws(cl[["MGMTp_methylation_status"]]) %in% c("un-methylated","Un-methylated"), 0L, NA)))
  clin <- cbind(clin, fr)
  cc2 <- !is.na(clin$age) & !is.na(clin$grade) & !is.na(clin$idh) & !is.na(clin$codel) & !is.na(clin$mgmt)
  dd <- clin[cc2, ]
  lp_state <- predict(fmod, newdata = dd[, st, drop = FALSE])                     # frozen state model
  g <- intersect(names(sigmap), rownames(ex)); rsig <- as.numeric(t(log2(ex[g, ss, drop = FALSE] + 1)) %*% sigmap[g])[cc2]
  fclin <- coxph(Surv(OS, ev) ~ age + factor(grade) + idh + codel + mgmt, data = dd)
  fcomb <- coxph(Surv(OS, ev) ~ lp_state + age + factor(grade) + idh + codel + mgmt, data = dd)
  lrt <- anova(fclin, fcomb, test = "LRT")$`Pr(>|Chi|)`[2]
  cat(sprintf("\n== CGGA-%s == n=%d ev=%d\n  state-fraction C=%.3f | gene-signature C=%.3f | clinical C=%.3f | clinical+state C=%.3f | added-value LRT p=%.2g\n",
              name, nrow(dd), sum(dd$ev), cind(lp_state, dd$OS, dd$ev), cind(rsig, dd$OS, dd$ev),
              cind(predict(fclin), dd$OS, dd$ev), cind(predict(fcomb), dd$OS, dd$ev), lrt))
  rows[[name]] <- data.frame(cohort = paste0("CGGA-", name), n = nrow(dd), events = sum(dd$ev),
                             C_state_fraction = round(cind(lp_state, dd$OS, dd$ev), 3),
                             C_gene_signature = round(cind(rsig, dd$OS, dd$ev), 3),
                             C_clinical = round(cind(predict(fclin), dd$OS, dd$ev), 3),
                             C_clinical_state = round(cind(predict(fcomb), dd$OS, dd$ev), 3),
                             added_value_LRT_p = signif(lrt, 3))
}
out <- do.call(rbind, rows)
write.table(out, file.path(SC, "state_model_performance.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nDONE -> results/state_composition/state_model_performance.tsv\n"); print(out, row.names = FALSE)
