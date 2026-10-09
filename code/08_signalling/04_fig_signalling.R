#!/usr/bin/env Rscript
# P2.4 - Figure for the immune/TAM->malignant signalling axis (Fig 9).
# A: top L-R pairs targeting malignant states. B: C-index comparison. C: score HR by chemo status.
# Outputs: results/figures/fig9_signalling.pdf + .jpg
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
SG <- "results/signalling"; FIG <- "results/figures"; dir.create(FIG, showWarnings = FALSE)

lr <- read.delim(file.path(SG, "cc_lr_pairs.tsv"), stringsAsFactors = FALSE)
perf <- read.delim(file.path(SG, "signalling_score_performance.tsv"), stringsAsFactors = FALSE)
mal <- c("AC", "NPC", "MES", "OPC")
to_mal <- lr[lr$target %in% mal, ]
to_mal <- to_mal[order(-to_mal$prob), ]; to_mal <- to_mal[!duplicated(paste(to_mal$ligand, to_mal$receptor, to_mal$source)), ]
top <- head(to_mal, 15)
lab <- paste0(top$ligand, "-", top$receptor, "  (", top$source, "->", top$target, ")")
col <- ifelse(top$source %in% c("Immune", "Endothelial"), "#C44E52", "#4C72B0")

draw <- function(dev) {
  if (dev == "pdf") pdf(file.path(FIG, "fig9_signalling.pdf"), width = 13, height = 4.4) else
                     jpeg(file.path(FIG, "fig9_signalling.jpg"), width = 3300, height = 1100, res = 300)
  par(mfrow = c(1, 3), mar = c(5, 12, 3, 1), las = 1)
  # A
  barplot(rev(top$prob), horiz = TRUE, col = rev(col), names.arg = rev(lab), cex.names = 0.6,
          xlab = "communication probability", main = "(d)  Top L-R pairs -> malignant")
  legend("bottomright", c("malignant-derived", "immune/stroma-derived"), fill = c("#4C72B0", "#C44E52"), bty = "n", cex = 0.7)
  # B
  par(mar = c(5, 5, 3, 1))
  m <- as.matrix(perf[, c("C_score_immune", "C_clinical", "C_clinical_score")]); rownames(m) <- perf$cohort
  bp <- barplot(t(m), beside = TRUE, ylim = c(0, 0.9), col = c("#C44E52", "#55A868", "#8172B3"),
                ylab = "C-index", main = "(e)  External discrimination", axes = FALSE)
  abline(h = 0.5, lty = 3, col = "grey60"); axis(2, at = seq(0, 0.9, 0.1))
  legend("topright", c("signalling score", "clinical", "clinical+score"), fill = c("#C44E52", "#55A868", "#8172B3"), bty = "n", cex = 0.7)
  axis(1, at = colMeans(bp), labels = rownames(m), tick = FALSE, line = -0.5)
  # C - score HR by chemo status
  par(mar = c(5, 6, 3, 1))
  ylim <- range(c(perf$lo_chemo, perf$hi_chemo, perf$lo_nochemo, perf$hi_nochemo), na.rm = TRUE)
  plot(NA, xlim = c(0.5, 4.5), ylim = log(range(ylim)), log = "y", xaxt = "n", xlab = "", ylab = "HR per SD (log scale)",
       main = "(f)  Score HR: chemo vs untreated")
  abline(h = 1, lty = 3, col = "grey60")
  xs <- c(1, 2, 3.5, 4.5)
  for (i in seq_len(nrow(perf))) {
    xi <- xs[c(1,2) + 2*(i-1)]
    segments(xi[1], log(perf$hi_chemo[i]), xi[1], log(perf$lo_chemo[i]), lwd = 2, col = "#C44E52")
    points(xi[1], log(perf$HR_chemo[i]), pch = 19, col = "#C44E52", cex = 1.3)
    segments(xi[2], log(perf$hi_nochemo[i]), xi[2], log(perf$lo_nochemo[i]), lwd = 2, col = "#4C72B0")
    points(xi[2], log(perf$HR_nochemo[i]), pch = 19, col = "#4C72B0", cex = 1.3)
  }
  axis(1, at = xs, labels = rep(c("chemo", "no chemo"), nrow(perf)), cex.axis = 0.75)
  text(c(1.5, 4), par("usr")[4], labels = perf$cohort, xpd = NA, cex = 0.85, font = 2)
  legend("topright", c("chemo", "untreated"), col = c("#C44E52", "#4C72B0"), lwd = 2, bty = "n", cex = 0.7)
  dev.off()
}
draw("pdf"); draw("jpg")
cat("saved fig9_signalling.pdf/.jpg\n")
