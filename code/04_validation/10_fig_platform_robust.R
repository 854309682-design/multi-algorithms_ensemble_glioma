#!/usr/bin/env Rscript
# P5.3 - Figure: platform-robust scoring (raw vs within-sample-rank srank vs per-gene-rank crank).
# Outputs: results/figures/fig11_platform_robust.pdf + .jpg
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
RES <- "results"; FIG <- file.path(RES, "figures")
d <- read.delim(file.path(RES, "platform_robust_scoring.tsv"), stringsAsFactors = FALSE)
rownames(d) <- paste0(d$cohort, ifelse(d$setting == "original", "", " (shift)"))
m <- as.matrix(d[, c("C_raw", "C_srank", "C_crank")])

draw <- function(dev) {
  if (dev == "pdf") pdf(file.path(FIG, "fig11_platform_robust.pdf"), width = 8.5, height = 4.5) else
                     jpeg(file.path(FIG, "fig11_platform_robust.jpg"), width = 2100, height = 1050, res = 300)
  par(mar = c(4, 5, 3, 1), las = 1)
  bp <- barplot(t(m), beside = TRUE, ylim = c(0, 0.9), col = c("#4C72B0", "#C44E52", "#55A868"),
                ylab = "C-index", main = "(d)  Platform-robust scoring (29-gene signature)", axes = FALSE,
                axisnames = FALSE, xpd = FALSE)
  abline(h = 0.5, lty = 3, col = "grey60"); axis(2, at = seq(0, 0.9, 0.1))
  axis(1, at = colMeans(bp), labels = rownames(m), tick = FALSE, line = 0.5, cex.axis = 0.85)
  legend("topright", c("raw (coef*expr)", "srank (within-sample rank)", "crank (per-gene rank)"),
         fill = c("#4C72B0", "#C44E52", "#55A868"), bty = "n", cex = 0.75)
  dev.off()
}
draw("pdf"); draw("jpg")
cat("saved fig11_platform_robust.pdf/.jpg\n")
