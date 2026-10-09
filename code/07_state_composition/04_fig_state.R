#!/usr/bin/env Rscript
# P1.5 - Figure for the state-composition model (Fig 8).
# A: C-index comparison (state-fraction vs gene-signature vs clinical vs combined) on CGGA.
# B: confounding |rho| of the state model vs the gene signature (grade, IDH).
# C: mean malignant/myeloid state fractions by WHO grade (TCGA).
# Outputs: results/figures/fig8_state_composition.pdf + .jpg
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

suppressMessages(library(survival))
setwd(ROOT)
SC <- "results/state_composition"; FIG <- "results/figures"; dir.create(FIG, showWarnings = FALSE)

perf <- read.delim(file.path(SC, "state_model_performance.tsv"), stringsAsFactors = FALSE)
cor  <- read.delim(file.path(SC, "state_confounder_cor.tsv"), stringsAsFactors = FALSE)
fr   <- read.delim(file.path(SC, "state_fractions_tcga.tsv"), row.names = 1, check.names = FALSE)
mol  <- read.delim("results/tcga_training_clinical_mol.tsv", stringsAsFactors = FALSE)
mol  <- mol[match(rownames(fr), mol$sample_submitter_id), ]

draw <- function(dev) {
  if (dev == "pdf") pdf(file.path(FIG, "fig8_state_composition.pdf"), width = 13, height = 4.2) else
                     jpeg(file.path(FIG, "fig8_state_composition.jpg"), width = 3250, height = 1050, res = 300)
  layout(matrix(1:4, nrow = 1), widths = c(1, 1, 1, 0.5))
  par(mar = c(5, 5, 3, 1), las = 1)

  # A - C-index comparison
  m <- as.matrix(perf[, c("C_state_fraction", "C_gene_signature", "C_clinical", "C_clinical_state")])
  rownames(m) <- perf$cohort
  bp <- barplot(t(m), beside = TRUE, ylim = c(0, 0.9), col = c("#4C72B0", "#DD8452", "#55A868", "#C44E52"),
                xlab = "", ylab = "C-index", main = "(a)  External discrimination", axes = FALSE, axisnames = FALSE)
  abline(h = 0.5, lty = 3, col = "grey60"); axis(2, at = seq(0, 0.9, 0.1))
  legend("topright", c("state", "gene", "clinical", "clinical+state"), fill = c("#4C72B0", "#DD8452", "#55A868", "#C44E52"),
         bty = "n", cex = 0.7)
  axis(1, at = colMeans(bp), labels = rownames(m), tick = FALSE, line = -0.5)

  # B - confounding |rho|
  r_state <- abs(cor$rho); r_gene <- abs(c(0.757, 0.806))   # gene-signature rho (grade, IDH)
  bb <- barplot(rbind(r_state, r_gene), beside = TRUE, ylim = c(0, 1), col = c("#4C72B0", "#DD8452"),
                names.arg = c("grade", "IDH", "MGMT", "age"), ylab = "|Spearman rho|",
                main = "(b)  Confounding with risk score", axes = FALSE)
  axis(2, at = seq(0, 1, 0.2)); legend("topright", c("state model", "gene signature"), fill = c("#4C72B0", "#DD8452"), bty = "n", cex = 0.7)

  # C - mean fractions by WHO grade (TCGA)
  st <- c("AC", "OPC", "NPC", "MES", "TAM-IS", "myeloid_other", "microglia")
  fr2 <- fr[, st]; g <- mol$who_grade
  agg <- t(sapply(sort(unique(na.omit(g))), function(gi) colMeans(fr2[which(g == gi), , drop = FALSE])))
  rownames(agg) <- paste0("G", sort(unique(na.omit(g))))
  par(mar = c(5, 5, 3, 1))
  barplot(t(agg), beside = FALSE, col = rainbow(length(st)), ylim = c(0, 1),
          xlab = "WHO grade", ylab = "mean fraction", main = "(c)  State vs grade (TCGA)", axes = FALSE)
  axis(2, at = seq(0, 1, 0.2))
  par(mar = c(5, 0.5, 3, 1))
  plot.new()
  legend("left", st, fill = rainbow(length(st)), bty = "n", cex = 0.7)
  dev.off()
}
draw("pdf"); draw("jpg")
cat("saved fig8_state_composition.pdf/.jpg\n")
