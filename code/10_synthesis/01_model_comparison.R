#!/usr/bin/env Rscript
# P6.1 (synthesis) - Unified head-to-head comparison of all models in this study, within each cohort.
# Features: risk (29-gene signature), state (state-fraction model LP), signal (14-gene immune->
# malignant signalling score), clinical (age + WHO grade + IDH + MGMT, common to both cohorts).
# Each model evaluated by patient-level 10-fold CV C-index inside the cohort + incremental LRT vs clinical.
# Output: results/synthesis/model_comparison.tsv
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
RES <- "results"; SC <- file.path(RES, "state_composition"); SG <- file.path(RES, "signalling")
SYN <- file.path(RES, "synthesis"); dir.create(SYN, showWarnings = FALSE)
D <- file.path(DATA_ROOT, "CGGA")

cind <- function(risk, t, e) { ok <- !is.na(risk) & !is.na(t) & !is.na(e); risk <- risk[ok]; t <- t[ok]; e <- e[ok]
  n <- length(risk); if (n < 10) return(NA); c <- 0; d <- 0
  for (i in seq_len(n)) { if (!e[i]) next; ti <- t[i]
    for (j in seq_len(n)) { if (i == j) next
      if (t[j] > ti || (t[j] == ti && !e[j])) { d <- d + 1
        if (risk[i] > risk[j]) c <- c + 1 else if (risk[i] == risk[j]) c <- c + 0.5 } } }
  if (d == 0) NA else c / d }
cvC <- function(d, rhs, K = 10) { n <- nrow(d); fold <- integer(n); o <- sample(n)
  for (k in 1:K) fold[o[seq(k, length(o), by = K)]] <- k; lp <- rep(NA, n)
  for (k in 1:K) { tr <- fold != k
    f <- tryCatch(coxph(as.formula(paste("Surv(time,event) ~", rhs)), data = d[tr, , drop = FALSE]), error = function(e) NULL)
    if (!is.null(f)) lp[!tr] <- tryCatch(predict(f, newdata = d[!tr, , drop = FALSE]), error = function(e) NA) }
  cind(lp, d$time, d$event) }

sig <- read.delim(file.path(RES, "signature_pc_coefs.tsv")); sigmap <- setNames(sig$coefficient, sig$gene)
panel <- readRDS(file.path(SG, "sg_panels.rds"))$immune
st <- setdiff(colnames(read.delim(file.path(SC, "state_fractions_tcga.tsv"), row.names = 1, check.names = FALSE)), "Oligodendrocyte")
score_z <- function(ex, p) { g <- intersect(p, rownames(ex)); colMeans(t(scale(t(as.matrix(ex[g, , drop = FALSE])))), na.rm = TRUE) }

# state-fraction model coefficients (frozen, from TCGA)
tcga <- read.delim(file.path(RES, "tcga_training_expr_log2cpm.tsv"), row.names = 1, check.names = FALSE)
mol <- read.delim(file.path(RES, "tcga_training_clinical_mol.tsv"), stringsAsFactors = FALSE); mol <- mol[match(colnames(tcga), mol$sample_submitter_id), ]
fr_t <- read.delim(file.path(SC, "state_fractions_tcga.tsv"), row.names = 1, check.names = FALSE)[, st]
fstat <- coxph(Surv(mol$OS_time_days, mol$OS_event) ~ ., data = fr_t)

feat <- function(name) {
  fr <- read.delim(file.path(SC, paste0("state_fractions_cgga", name, ".tsv")), row.names = 1, check.names = FALSE)[, st]
  tc <- tempfile(); system(sprintf("unzip -p '%s/CGGA.mRNAseq_%s_clinical.zip' 'CGGA.mRNAseq_%s_clinical.20200506.txt' > '%s'", D, name, name, tc))
  cl <- read.delim(tc, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE); names(cl)[1] <- "CGGA_ID"
  tr <- tempfile(); system(sprintf("unzip -p '%s/CGGA.mRNAseq_%s.RSEM-genes.zip' 'CGGA.mRNAseq_%s.RSEM-genes.20200506.txt' > '%s'", D, name, name, tr))
  ex <- read.delim(tr, sep = "\t", check.names = FALSE); rownames(ex) <- ex[, 1]; ex <- log2(as.matrix(ex[, -1]) + 1)
  cl <- cl[trimws(cl$PRS_type) == "Primary", ]; cl$OS <- suppressWarnings(as.numeric(cl$OS)); cl$Censor <- suppressWarnings(as.numeric(cl$Censor))
  cl <- cl[!is.na(cl$OS) & cl$OS > 0 & !is.na(cl$Censor), ]
  ss <- intersect(colnames(ex), trimws(cl$CGGA_ID)); ex <- ex[, ss]; cl <- cl[match(ss, trimws(cl$CGGA_ID)), ]; fr <- fr[ss, ]
  gs <- intersect(names(sigmap), rownames(ex))
  d <- data.frame(time = cl$OS, event = cl$Censor,
                  age = suppressWarnings(as.numeric(trimws(cl$Age))),
                  grade = factor({ g <- trimws(cl$Grade); ifelse(grepl("IV", g), 4, ifelse(grepl("III", g), 3, ifelse(grepl("II\\b", g), 2, NA))) }),
                  idh = ifelse(trimws(cl$IDH_mutation_status) == "Mutant", 1, ifelse(trimws(cl$IDH_mutation_status) == "Wildtype", 0, NA)),
                  mgmt = ifelse(trimws(cl[["MGMTp_methylation_status"]]) %in% c("methylated","Methylated"), 1, ifelse(trimws(cl[["MGMTp_methylation_status"]]) %in% c("un-methylated","Un-methylated"), 0, NA)),
                  risk = as.numeric(t(ex[gs, , drop = FALSE]) %*% sigmap[gs]),
                  state = predict(fstat, newdata = fr), signal = score_z(ex, panel))
  d[complete.cases(d), ]
}
CLIN <- "age + grade + idh + mgmt"
models <- list(clinical = CLIN, gene_signature = "risk", state_composition = "state", signalling = "signal",
               clinical_signature = paste(CLIN, "+ risk"), clinical_state = paste(CLIN, "+ state"),
               clinical_signalling = paste(CLIN, "+ signal"), clinical_full = paste(CLIN, "+ risk + state + signal"))
rows <- list()
for (name in c("693","325")) {
  d <- feat(name)
  fclin <- coxph(as.formula(paste("Surv(time,event) ~", CLIN)), data = d)
  for (m in names(models)) {
    f <- coxph(as.formula(paste("Surv(time,event) ~", models[[m]])), data = d)
    lrt <- if (m == "clinical") NA else signif(anova(fclin, f, test = "LRT")$`Pr(>|Chi|)`[2], 3)
    rows[[paste(name, m)]] <- data.frame(cohort = paste0("CGGA-", name), model = m, n = nrow(d), events = sum(d$event),
                                         CV_C_index = round(cvC(d, models[[m]]), 3), vs_clinical_LRT_p = lrt)
  }
}
out <- do.call(rbind, rows)
write.table(out, file.path(SYN, "model_comparison.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("DONE -> results/synthesis/model_comparison.tsv\n"); print(out, row.names = FALSE)
