#!/usr/bin/env Rscript
# Stage 5 (clinical, extended) -- clinico-molecular nomogram + calibration + DCA for the
# frozen PRIMARY 29-gene protein-coding signature.
# Predictors: risk (frozen linear predictor) + age + WHO grade + IDH + MGMT (complete cases).
# Extends 01_nomogram_dca.R (risk + age only) after the molecular annotation join
# (code/02_bulk_prep/02_tcga_molecular_join.R). The risk + age model is refit on the same
# complete-case set so the two C-indices are directly comparable.
# Plots:  results/figures/05b_nomogram.pdf/.jpg, 05b_calibration_3yr.pdf/.jpg, 05b_dca_3yr.pdf/.jpg
# Tables: results/05b_nomogram_cox_models.tsv, results/05b_nomogram_cindex.tsv
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

library(survival)
library(rms)
library(dcurves)

set.seed(20260909)

setwd(ROOT)
RES <- "results"; FIG <- "results/figures"

dir.create(FIG, showWarnings = FALSE)

# Harrell's C via the project-standard manual concordance (cf. 03_modeling/06_confounding_audit.R)
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
sm <- setNames(sig$coefficient, sig$gene)
cat("primary signature genes:", length(sm), "(protein-coding)\n")

samples <- clin$sample_submitter_id
expr <- expr[, intersect(colnames(expr), samples)]
clin <- clin[match(colnames(expr), clin$sample_submitter_id), ]

g <- intersect(names(sm), rownames(expr))
cat("signature genes present in the training matrix:", length(g), "/", length(sm), "\n")

clin$risk  <- as.numeric(scale(t(as.matrix(expr[g, , drop = FALSE])) %*% sm[g]))  # per SD (matches the confounding audit)
clin$age   <- as.numeric(clin$age_at_index)
clin$grade <- factor(clin$who_grade)
clin$idh   <- factor(clin$idh_mutant, levels = c(0, 1), labels = c("WT", "mutant"))
clin$mgmt  <- factor(clin$mgmt_meth, levels = c(0, 1), labels = c("unmethylated", "methylated"))

keep <- !is.na(clin$risk) & !is.na(clin$age) & !is.na(clin$grade) & !is.na(clin$idh) &
        !is.na(clin$mgmt) & !is.na(clin$OS_time_days) & clin$OS_time_days > 0 & !is.na(clin$OS_event)
df <- clin[keep, ]
cat("extended nomogram dataset n =", nrow(df), " events =", sum(df$OS_event), "\n")
cat("grade distribution:\n"); print(table(df$grade))
cat("IDH:", paste(names(table(df$idh)), table(df$idh), collapse = " / "),
    " MGMT:", paste(names(table(df$mgmt)), table(df$mgmt), collapse = " / "), "\n")

dd <- datadist(df); options(datadist = "dd")

f_ref  <- cph(Surv(OS_time_days, OS_event) ~ risk + age, data = df, x = TRUE, y = TRUE, surv = TRUE)
f_full <- cph(Surv(OS_time_days, OS_event) ~ risk + age + grade + idh + mgmt,
              data = df, x = TRUE, y = TRUE, surv = TRUE)
cat("\n--- risk + age (reference) ---\n"); print(f_ref)
cat("\n--- risk + age + grade + IDH + MGMT (extended) ---\n"); print(summary(f_full))

hr_of <- function(f, tag) {
  cf <- coef(f); se <- sqrt(diag(vcov(f)))
  data.frame(model = tag, term = names(cf), HR = exp(cf),
             lo = exp(cf - 1.96 * se), hi = exp(cf + 1.96 * se),
             p = 2 * pnorm(-abs(cf / se)), row.names = NULL)
}
hr <- rbind(hr_of(f_ref, "risk+age"), hr_of(f_full, "risk+age+grade+IDH+MGMT"))
hr$HR <- round(hr$HR, 3); hr$lo <- round(hr$lo, 3); hr$hi <- round(hr$hi, 3)
hr$p <- signif(hr$p, 3)
write.table(hr, file.path(RES, "05b_nomogram_cox_models.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nHazard ratios:\n"); print(hr)

ci_tab <- data.frame(
  model = c("risk+age", "risk+age+grade+IDH+MGMT"),
  C_index = round(c(cind(predict(f_ref), df$OS_time_days, df$OS_event),
                    cind(predict(f_full), df$OS_time_days, df$OS_event)), 3),
  n = nrow(df), events = sum(df$OS_event))
write.table(ci_tab, file.path(RES, "05b_nomogram_cindex.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nDiscrimination:\n"); print(ci_tab)

# --- nomogram (1/3/5-yr OS from the extended model) ---
surv <- Survival(f_full)
nom <- nomogram(f_full, fun = list(function(x) 1 - surv(365, x),
                                   function(x) 1 - surv(1095, x),
                                   function(x) 1 - surv(1825, x)),
                funlabel = c("1-yr OS", "3-yr OS", "5-yr OS"))
pdf(file.path(FIG, "05b_nomogram.pdf"), width = 10, height = 7)
plot(nom)
dev.off()
jpeg(file.path(FIG, "05b_nomogram.jpg"), width = 2000, height = 1400, res = 220)
plot(nom)
dev.off()
cat("Nomogram saved.\n")

# --- 3-year calibration (transparent: KM-observed vs predicted 3-yr OS by risk decile) ---
S3  <- as.numeric(surv(1095, predict(f_full)))   # per-patient predicted 3-yr OS (nomogram machinery)
cat("distinct predicted 3-yr OS values:", length(unique(S3)), "of", length(S3), "\n")
grp <- cut(rank(S3, ties.method = "first"), breaks = 10, labels = FALSE)
cal <- do.call(rbind, lapply(sort(unique(grp)), function(k) {
  i <- which(grp == k)
  ss <- summary(survfit(Surv(OS_time_days, OS_event) ~ 1, data = df[i, ]), times = 1095, extend = TRUE)
  data.frame(bin = k, n = length(i), predicted = round(mean(S3[i]), 3),
             observed = round(as.numeric(ss$surv), 3), lo = round(as.numeric(ss$lower), 3),
             hi = round(as.numeric(ss$upper), 3))
}))
write.table(cal, file.path(RES, "05b_calibration_3yr.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n3-yr calibration (predicted vs KM-observed):\n"); print(cal)

cal_plot <- function() {
  plot(cal$predicted, cal$observed, xlim = c(0, 1), ylim = c(0, 1), pch = 19, cex = 1.1,
       xlab = "Predicted 3-yr OS", ylab = "Observed 3-yr OS (KM)",
       main = "Calibration (3-year OS)")
  segments(cal$predicted, cal$lo, cal$predicted, cal$hi, lwd = 1.2)
  abline(0, 1, lty = 3, col = "grey40")
}
pdf(file.path(FIG, "05b_calibration_3yr.pdf"), width = 6, height = 6); cal_plot(); dev.off()
jpeg(file.path(FIG, "05b_calibration_3yr.jpg"), width = 1200, height = 1200, res = 220); cal_plot(); dev.off()
cat("Calibration saved.\n")

# --- DCA (3-yr) ---
dfd <- data.frame(time = df$OS_time_days, event = df$OS_event, nomogram = 1 - S3)
dca_f <- dca(Surv(time, event) ~ nomogram, data = dfd, time = 1095)
p <- plot(dca_f) + ggplot2::theme_bw() +
     ggplot2::theme(panel.grid = ggplot2::element_blank()) +
     ggplot2::ggtitle("Decision curve analysis (3-year OS)")
dca_tab <- as.data.frame(dca_f$dca)
write.table(dca_tab, file.path(RES, "05b_dca_3yr.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("DCA columns:", paste(names(dca_tab), collapse = ", "), "\n")
cat("DCA rows:", nrow(dca_tab), "| distinct nomogram net-benefit values:",
    length(unique(dca_tab$net_benefit[dca_tab$variable == "nomogram"])), "\n")
pdf(file.path(FIG, "05b_dca_3yr.pdf"), width = 7, height = 6)
print(p)
dev.off()
jpeg(file.path(FIG, "05b_dca_3yr.jpg"), width = 1400, height = 1200, res = 220)
print(p)
dev.off()
cat("DCA saved.\n")

writeLines(paste0("NOMOG_DCA_MOL_DONE n=", nrow(df), " events=", sum(df$OS_event),
                  " C_full=", round(cind(predict(f_full), df$OS_time_days, df$OS_event), 3),
                  " C_ref=", round(cind(predict(f_ref), df$OS_time_days, df$OS_event), 3)),
           file.path(RES, "05b_nomogram_dca_note.txt"))
cat("DONE\n")
