#!/usr/bin/env Rscript
# Stage 5 (clinical): nomogram (1/3/5-y OS) + calibration + DCA for the frozen signature.
# Predictors: risk score (frozen LP) + age + grade. (IDH/MGMT to be added after molecular annotation.)
# Plots: 05_nomogram.pdf/png, 05_dca.pdf/png (PDF + JPEG, prefixed by script step id).
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

suppressMessages({library(survival); library(rms); library(dcurves)})

set.seed(20260909)

setwd(ROOT)
RES <- "results"; FIG <- "results/figures"

dir.create(FIG, showWarnings = FALSE)

expr <- read.delim(file.path(RES, "tcga_training_expr_log2cpm.tsv"), row.names = 1, check.names = FALSE)
clin <- read.delim(file.path(RES, "tcga_training_clinical.tsv"), stringsAsFactors = FALSE)
sig  <- read.delim(file.path(RES, "signature_genes_coefs.tsv"), stringsAsFactors = FALSE)
sm <- setNames(sig$coefficient, sig$gene)

samples <- clin$sample_submitter_id
expr <- expr[, intersect(colnames(expr), samples)]
clin <- clin[match(colnames(expr), clin$sample_submitter_id), ]

g <- intersect(names(sm), rownames(expr))
risk <- as.numeric(t(as.matrix(expr[g, , drop = FALSE])) %*% sm[g])

# build clinical covariates; grade is incomplete in TCGA clinical (G4/GBM mostly missing).
# Use risk + age (complete) for the nomogram; grade + IDH/MGMT to be added after molecular annotation.
clin$age <- as.numeric(clin$age_at_index)
clin$risk <- risk
df <- clin[!is.na(clin$age) & clin$OS_time_days > 0, ]
cat("nomogram dataset n =", nrow(df), " events =", sum(df$OS_event), "\n")
clin$grade <- NULL

# --- nomogram model ---
dd <- datadist(df); options(datadist = "dd")
f <- cph(Surv(OS_time_days, OS_event) ~ risk + age, data = df, x = TRUE, y = TRUE, surv = TRUE)
cat("Nomogram model\n"); print(summary(f))

surv <- Survival(f)
# nomogram at 1/3/5 years (in days)
nom <- nomogram(f, fun = list(function(x) 1 - surv(365, x),
                              function(x) 1 - surv(1095, x),
                              function(x) 1 - surv(1825, x)),
                funlabel = c("1-yr", "3-yr", "5-yr"))
pdf(file.path(FIG, "05_nomogram.pdf"), width = 9, height = 8)
plot(nom)
dev.off()
jpeg(file.path(FIG, "05_nomogram.jpg"), width = 1800, height = 1600, res = 220)
plot(nom)
dev.off()
cat("Nomogram saved.\n")

# --- calibration (risk + clinical vs observed) at 3yr ---
cal <- calibrate(f, u = 1095, cmethod = "KM", m = 60, B = 200)
# NOTE: rms::plot.calibrate draws its caption via title(sub = ...), which lands on the SAME
# margin line as the x-axis title and overprints "Predicted 3-yr survival". Suppress the
# built-in caption (subtitles = FALSE) and re-draw it two margin lines lower (same strings).
cN <- attr(cal, "n"); cD <- attr(cal, "d"); cP <- attr(cal, "p")
cM <- attr(cal, "m"); cB <- attr(cal, "B"); cWhat <- attr(cal, "what")
pdf(file.path(FIG, "05_calibration_3yr.pdf"), width = 6, height = 6)
par(mar = c(7.2, 4.6, 4, 2.2))
plot(cal, xlab = "Predicted 3-yr survival", ylab = "Observed", subtitles = FALSE)
mtext(paste0("n=", cN, " d=", cD, " p=", cP, ", ", cM, " subjects per group"),
      side = 1, line = 4.5, adj = 0, cex = 0.65, col = "grey30")
mtext(paste0("X - resampling optimism added, B=", cB),
      side = 1, line = 4.5, adj = 1, cex = 0.65, col = "grey30")
mtext("Gray: ideal", side = 1, line = 5.5, adj = 0, cex = 0.65, col = "grey30")
mtext(paste0("Based on ", cWhat), side = 1, line = 5.5, adj = 1, cex = 0.65, col = "grey30")
dev.off()
cat("Calibration saved.\n")

# --- DCA (predicted 3-yr probability vs treat-all vs treat-none) ---
# Convert Cox LP -> 3-yr event probability via rms survest, then run dcurves DCA.
p3 <- as.numeric(1 - survest(f, times = 1095)$surv)
dfd <- data.frame(time = df$OS_time_days, event = df$OS_event, p3 = p3)
dca_f <- dca(Surv(time, event) ~ p3, data = dfd, time = 1095)
pdf(file.path(FIG, "05_dca_3yr.pdf"), width = 7, height = 6)
plot(dca_f)
dev.off()
cat("DCA plots saved.\n")

# net benefit summary at a couple threshold probabilities
nb <- as.data.frame(dca_f$dca)
cat("DCA net benefit (risk model, 3-yr), head:\n")
print(head(nb, 6))
writeLines(paste0("NOMOG_DCA_DONE n=", nrow(df), " events=", sum(df$OS_event),
                  " model C=", round(concordance(f)$concordance, 3)),
           file.path(RES, "nomogram_dca_note.txt"))
cat("DONE\n")
