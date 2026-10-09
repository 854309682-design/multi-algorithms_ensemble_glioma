#!/usr/bin/env Rscript
# Fig 6 -- biology closure: signature-gene importance, single-cell state traceback, GO BP enrichment.
# Inputs:  results/signature_importance_states.tsv, results/signature_GO_BP.tsv
# Outputs: results/figures/fig6_biology_closure.pdf / .jpg, results/fig6_signature_state_counts.tsv
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

imp <- read.delim(file.path(RES, "signature_importance_states.tsv"), stringsAsFactors = FALSE)
go  <- read.delim(file.path(RES, "signature_GO_BP.tsv"), stringsAsFactors = FALSE)
cat("signature genes:", nrow(imp), " GO BP terms:", nrow(go), "\n")

imp$state[is.na(imp$state) | imp$state == ""] <- "unassigned"
sc <- as.data.frame(table(imp$state))
names(sc) <- c("state", "n_genes")
sc <- sc[order(-sc$n_genes), ]
write.table(sc, file.path(RES, "fig6_signature_state_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

top <- head(imp[order(-imp$importance), ], 15)
go  <- go[order(go$p.adjust), ]
gt  <- head(go, 8)
gt$Desc <- substr(gt$Description, 1, 42)

# colouring by state group (malignant state vs TAM vs hypoxia/prolif)
grp <- function(s) ifelse(grepl("TAM", s), "TAM", ifelse(grepl("Hypox", s), "Hypoxia",
             ifelse(grepl("Prolif", s), "Proliferation", "Malignant state")))
gcol <- c("Malignant state" = "#4363d8", "TAM" = "#3cb44b", "Hypoxia" = "#e6194B",
          "Proliferation" = "#f58231", "unassigned" = "grey60")

draw <- function() {
  par(mfrow = c(1, 3), mar = c(4.2, 8, 3, 1))
  bp <- barplot(rev(top$importance), names.arg = rev(top$gene), horiz = TRUE, las = 1,
                col = gcol[grp(rev(top$state))], cex.names = 0.8, xlab = "|coef| x SD",
                main = "(d)  Signature-gene importance")
  legend("bottomright", names(gcol)[1:4], fill = gcol[1:4], bty = "n", cex = 0.65)

  par(mar = c(4.2, 7, 3, 1))
  barplot(rev(sc$n_genes), names.arg = rev(sc$state), horiz = TRUE, las = 1, col = "steelblue",
          cex.names = 0.8, xlab = "signature genes", main = "(e)  Single-cell state traceback")

  par(mar = c(6, 4.2, 3, 4))
  pt <- seq(min(-log10(gt$p.adjust)), max(-log10(gt$p.adjust)), length.out = 50)
  ram <- colorRampPalette(c("grey85", "firebrick"))(50)
  plot(rev(gt$FoldEnrichment), seq_len(nrow(gt)), pch = 19, cex = rev(gt$Count) / 2.2,
       col = ram[findInterval(rev(-log10(gt$p.adjust)), pt)], yaxs = "i", xlim = c(0, max(gt$FoldEnrichment) * 1.1),
       yaxt = "n", ylab = "", xlab = "fold enrichment", main = "(f)  GO BP enrichment")
  axis(2, at = seq_len(nrow(gt)), labels = rev(gt$Desc), las = 1, cex.axis = 0.62)
  legend("bottomright", c("p.adj 0.2", "p.adj 1e-4"), pch = 19, col = ram[c(50, 1)], bty = "n", cex = 0.65, title = "-log10", title.cex = 0.65)
  mtext("dot size = gene count", side = 3, line = 0.2, cex = 0.6)
}
pdf(file.path(FIG, "fig6_biology_closure.pdf"), width = 15, height = 5.5); draw(); dev.off()
jpeg(file.path(FIG, "fig6_biology_closure.jpg"), width = 3300, height = 1250, res = 300); draw(); dev.off()
cat("DONE\n")
