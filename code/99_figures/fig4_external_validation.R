#!/usr/bin/env Rscript
# Fig 4 — External validation (backbone): forest plot of C-index + KM curves (frozen signature).
# Plots: results/figures/fig4_forest.pdf/.jpg, fig4_km_cgga693.pdf/.jpg, fig4_km_cgga325.pdf/.jpg
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
RES <- "results"; FIG <- "results/figures"; DCGGA <- file.path(DATA_ROOT, "CGGA")

dir.create(FIG, showWarnings = FALSE)

sig <- read.delim(file.path(RES, "signature_pc_coefs.tsv"), stringsAsFactors = FALSE)
sm <- setNames(sig$coefficient, sig$gene)

# ---------- forest plot of C-index ----------
s <- read.delim(file.path(RES, "external_validation_summary_pc.tsv"), stringsAsFactors = FALSE)
s <- s[order(s$C_index), ]

pdf(file.path(FIG, "fig4_forest.pdf"), width = 7, height = 4)
par(mar = c(4, 9, 3, 2))
plot(s$C_index, seq_len(nrow(s)), xlim = c(0.4, 0.85), yaxt = "n", pch = 19, cex = 1.2,
     xlab = "C-index (95% CI)", ylab = "", main = "External validation (frozen signature)")
segments(s$lo, seq_len(nrow(s)), s$hi, seq_len(nrow(s)), lwd = 2)
abline(v = 0.5, lty = 3, col = "grey50")
abline(v = 0.760, lty = 2, col = "red")
axis(2, at = seq_len(nrow(s)), labels = s$cohort, las = 1)
dev.off()
jpeg(file.path(FIG, "fig4_forest.jpg"), width = 1400, height = 800, res = 200)
par(mar = c(4, 9, 3, 2))
plot(s$C_index, seq_len(nrow(s)), xlim = c(0.4, 0.85), yaxt = "n", pch = 19, cex = 1.2,
     xlab = "C-index (95% CI)", ylab = "", main = "External validation (frozen signature)")
segments(s$lo, seq_len(nrow(s)), s$hi, seq_len(nrow(s)), lwd = 2)
abline(v = 0.5, lty = 3, col = "grey50")
axis(2, at = seq_len(nrow(s)), labels = s$cohort, las = 1)
dev.off()

# ---------- KM curves per CGGA cohort (high/low by median) ----------
km_plot <- function(name, tag) {
  zf <- paste0(DCGGA, "/CGGA.mRNAseq_", name, ".RSEM-genes.zip")
  tmp <- tempfile(fileext = ".txt")
  system(paste0("unzip -p '", zf, "' 'CGGA.mRNAseq_", name, ".RSEM-genes.20200506.txt' > '", tmp, "'"))
  ex <- read.delim(tmp, sep = "\t", check.names = FALSE); rownames(ex) <- ex[, 1]; ex <- as.matrix(ex[, -1])
  clf <- paste0(DCGGA, "/CGGA.mRNAseq_", name, "_clinical_clean.tsv")
  if (file.exists(clf)) {
    clin <- read.delim(clf, stringsAsFactors = FALSE)
  } else {
    tc <- tempfile(fileext = ".txt")
    system(paste0("unzip -p '", DCGGA, "/CGGA.mRNAseq_", name, "_clinical.zip' 'CGGA.mRNAseq_", name, "_clinical.20200506.txt' > '", tc, "'"))
    clin <- read.delim(tc, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE); names(clin)[1] <- "CGGA_ID"
  }
  prim <- clin[clin$PRS_type == "Primary", ]; prim$OS <- as.numeric(prim$OS); prim$Censor <- as.numeric(prim$Censor)
  clinp <- prim[!is.na(prim$OS) & prim$OS > 0, ]
  samples <- intersect(colnames(ex), clinp$CGGA_ID); ex <- ex[, samples]; clinp <- clinp[match(samples, clinp$CGGA_ID), ]
  X <- log2(as.matrix(ex) + 1); g <- intersect(names(sm), rownames(X))
  risk <- as.numeric(t(X[g, , drop = FALSE]) %*% sm[g])
  grp <- factor(ifelse(risk > median(risk), "High", "Low"), levels = c("Low", "High"))
  fit <- survfit(Surv(clinp$OS, clinp$Censor) ~ grp)
  lr <- survdiff(Surv(clinp$OS, clinp$Censor) ~ grp); p <- 1 - pchisq(lr$chisq, 1)
  pdf(file.path(FIG, paste0("fig4_km_", tag, ".pdf")), width = 5, height = 5)
  plot(fit, col = c("blue", "red"), lwd = 2, xlab = "Overall survival (days)", ylab = "OS probability",
       main = paste0("CGGA-", name, " (log-rank p = ", format(p, digits = 2), ")"))
  legend("topright", c("Low risk", "High risk"), col = c("blue", "red"), lwd = 2, bty = "n")
  dev.off()
  jpeg(file.path(FIG, paste0("fig4_km_", tag, ".jpg")), width = 1000, height = 1000, res = 200)
  plot(fit, col = c("blue", "red"), lwd = 2, xlab = "Overall survival (days)", ylab = "OS probability",
       main = paste0("CGGA-", name, " (log-rank p = ", format(p, digits = 2), ")"))
  legend("topright", c("Low risk", "High risk"), col = c("blue", "red"), lwd = 2, bty = "n")
  dev.off()
  cat("KM saved for CGGA-", name, " p=", format(p, digits = 3), "\n", sep = "")
}
km_plot("693", "cgga693")
km_plot("325", "cgga325")
cat("DONE\n")
