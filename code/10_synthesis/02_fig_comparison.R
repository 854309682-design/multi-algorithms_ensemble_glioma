#!/usr/bin/env Rscript
# P6.2 - Figure: unified model comparison (10-fold CV C-index by model, per cohort).
# Outputs: results/figures/fig12_model_comparison.pdf + .jpg
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

setwd(ROOT)
SYN <- "results/synthesis"; FIG <- "results/figures"
d <- read.delim(file.path(SYN, "model_comparison.tsv"), stringsAsFactors = FALSE)
ord <- c("clinical","gene_signature","state_composition","signalling",
         "clinical_signature","clinical_state","clinical_signalling","clinical_full")
d$model <- factor(d$model, levels = ord)
m <- t(tapply(d$CV_C_index, list(d$model, d$cohort), identity)[ord, ])

draw <- function(dev) {
  if (dev == "pdf") pdf(file.path(FIG, "fig12_model_comparison.pdf"), width = 11.5, height = 4.6) else
                     jpeg(file.path(FIG, "fig12_model_comparison.jpg"), width = 2800, height = 1100, res = 300)
  layout(matrix(1:2, nrow = 1), widths = c(1, 0.28))
  lab <- c(clinical = "clinical", gene_signature = "gene sig", state_composition = "state comp",
           signalling = "signalling", clinical_signature = "clin+sig", clinical_state = "clin+state",
           clinical_signalling = "clin+signal", clinical_full = "clin+full")
  par(mar = c(8, 5, 3, 0.5), las = 1)
  bp <- barplot(m, beside = TRUE, ylim = c(0, 0.85), col = c("#4C72B0", "#C44E52"),
                ylab = "10-fold CV C-index", main = "(a)  Unified model comparison", axes = FALSE, axisnames = FALSE)
  abline(h = 0.5, lty = 3, col = "grey60"); axis(2, at = seq(0, 0.8, 0.1))
  axis(1, at = colMeans(bp), labels = lab[ord], tick = FALSE, line = -0.5, cex.axis = 0.75, las = 2)
  par(mar = c(8, 0.5, 3, 1))
  plot.new()
  legend("left", colnames(m), fill = c("#4C72B0", "#C44E52"), bty = "n", cex = 0.9)
  dev.off()
}
draw("pdf"); draw("jpg")
cat("saved fig12_model_comparison.pdf/.jpg\n")
