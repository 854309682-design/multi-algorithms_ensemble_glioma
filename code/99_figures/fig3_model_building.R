#!/usr/bin/env Rscript
# Fig 3 -- model building: univariate Cox screen, LASSO-Cox shrinkage, >=6-learner CV.
# Inputs:  results/univariate_cox.tsv, lasso_genes.tsv, nested_cv_cindex.tsv, ensemble_cv_cindex.tsv
# Outputs: results/figures/fig3_model_building.pdf / .jpg, results/fig3_model_building_summary.tsv
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

set.seed(20260909)

setwd(ROOT)
RES <- "results"; FIG <- "results/figures"

uni  <- read.delim(file.path(RES, "univariate_cox.tsv"), stringsAsFactors = FALSE)
las  <- read.delim(file.path(RES, "lasso_genes.tsv"), stringsAsFactors = FALSE)
nest <- read.delim(file.path(RES, "nested_cv_cindex.tsv"), stringsAsFactors = FALSE)
ens  <- read.delim(file.path(RES, "ensemble_cv_cindex.tsv"), stringsAsFactors = FALSE)

uni$mlp <- -log10(uni$p)
sig <- uni$p < 0.05
cat("univariate tested:", nrow(uni), " p<0.05:", sum(sig), "\n")

nest$Cindex <- as.numeric(nest$Cindex); ens$Cindex <- as.numeric(ens$Cindex)
nest <- nest[!is.na(nest$Cindex), ]; ens <- ens[!is.na(ens$Cindex), ]
m <- merge(nest, ens, by = "learner", suffixes = c("_nested", "_cv"))
m$Cindex_nested <- as.numeric(m$Cindex_nested); m$Cindex_cv <- as.numeric(m$Cindex_cv)
m <- m[!is.na(m$Cindex_nested) & !is.na(m$Cindex_cv), ]   # drop SurvSVM (NA) and the stacking row
stk <- as.numeric(ens$Cindex[ens$learner == "stacking"])
o <- order(m$Cindex_nested)

draw <- function() {
  par(mfrow = c(1, 3), mar = c(4.2, 4.2, 3, 1))
  plot(uni$coef, uni$mlp, pch = 19, cex = 0.45, col = ifelse(sig, "firebrick", "grey72"),
       xlab = "univariate Cox coefficient", ylab = "-log10(p)",
       main = paste0("(a)  Univariate screen (", sum(sig), "/", nrow(uni), " p<0.05)"))
  abline(h = -log10(0.05), lty = 3, col = "grey40")
  bp <- barplot(las$n_genes, names.arg = las$lambda, col = "steelblue", ylim = c(0, 40),
                xlab = "LASSO-Cox lambda", ylab = "genes retained", main = "(b)  LASSO shrinkage")
  text(bp, las$n_genes, labels = las$n_genes, pos = 3, cex = 0.9)
  barplot(t(as.matrix(m[o, c("Cindex_nested", "Cindex_cv")])), beside = TRUE,
          names.arg = m$learner[o], col = c("grey35", "steelblue"), las = 2, cex.names = 0.8,
          ylim = c(0, 1), ylab = "C-index", main = "(c)  Learner CV (nested vs ensemble)")
  legend("topleft", c("nested (honest)", "ensemble"), fill = c("grey35", "steelblue"),
         bty = "n", cex = 0.75)
}
pdf(file.path(FIG, "fig3_model_building.pdf"), width = 13, height = 4.5); draw(); dev.off()
jpeg(file.path(FIG, "fig3_model_building.jpg"), width = 2900, height = 1000, res = 300); draw(); dev.off()

# NOTE: nested stacking (0.839) was only printed by 03_nested_cv.R and never persisted to a
# results file; it is taken verbatim from the frozen record (decision_log.md DEC-4, handoff Stage 3).
summ <- data.frame(
  item = c("genes_tested", "univariate_p05", "lasso_lambda_min", "lasso_lambda_1se",
           "best_single_learner", "best_single_C_nested", "stacking_C_nested_frozenDEC4",
           "stacking_C_ensemble_nonNested", "decision"),
  value = c(nrow(uni), sum(sig), las$n_genes[las$lambda == "min"], las$n_genes[las$lambda == "1se"],
            m$learner[order(-m$Cindex_nested)][1], round(max(m$Cindex_nested), 3),
            0.839, round(stk, 3),
            "best-single = CoxPH (pre-registered: nested stacking 0.839 < best 0.834 + 0.01)"))
write.table(summ, file.path(RES, "fig3_model_building_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(summ, row.names = FALSE)
cat("DONE\n")
