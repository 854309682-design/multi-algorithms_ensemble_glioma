#!/usr/bin/env Rscript
# P1.3 - Deconvolve bulk transcriptomes into malignant + myeloid state fractions with MuSiC (v1.0.0).
# Reference: results/state_composition/ref_counts.mtx (+ ref_genes.txt, ref_cells.tsv).
# Bulk: TCGA training (log2CPM -> CPM) and CGGA-693/325 (RSEM -> linear).
# Output: results/state_composition/state_fractions_{tcga,cgga693,cgga325}.tsv (samples x states)
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

suppressMessages({library(MuSiC); library(SingleCellExperiment); library(Matrix)})
setwd(ROOT)
RES <- "results"; SC <- file.path(RES, "state_composition"); D <- file.path(DATA_ROOT, "CGGA")

sc <- as.matrix(readMM(file.path(SC, "ref_counts.mtx")))
genes <- readLines(file.path(SC, "ref_genes.txt"))
meta <- read.delim(file.path(SC, "ref_cells.tsv"), stringsAsFactors = FALSE)
rownames(sc) <- genes; colnames(sc) <- meta$cell
sce <- SingleCellExperiment(assays = list(counts = sc),
                            colData = DataFrame(cellType = meta$cellType, sampleID = meta$sampleID))
cat("reference:", nrow(sc), "genes x", ncol(sc), "cells;", length(unique(meta$cellType)), "states\n")

mus <- function(bulk, tag) {
  g <- intersect(rownames(bulk), rownames(sc))
  bulk <- bulk[g, , drop = FALSE]
  est <- music_prop(bulk.mtx = bulk, sc.sce = sce, clusters = "cellType", samples = "sampleID",
                    select.ct = unique(meta$cellType), verbose = FALSE)
  P <- as.data.frame(est$Est.prop.weighted)
  write.table(P, file.path(SC, paste0("state_fractions_", tag, ".tsv")), sep = "\t", quote = FALSE)
  cat(tag, ": deconvolved", nrow(P), "samples x", ncol(P), "states; genes used:", length(g), "\n")
}

# TCGA training (log2CPM -> CPM)
tcga <- read.delim(file.path(RES, "tcga_training_expr_log2cpm.tsv"), row.names = 1, check.names = FALSE)
mus(2^as.matrix(tcga), "tcga")

# CGGA-693 / 325 (RSEM -> linear)
for (name in c("693", "325")) {
  tr <- tempfile(fileext = ".txt")
  system(sprintf("unzip -p '%s/CGGA.mRNAseq_%s.RSEM-genes.zip' 'CGGA.mRNAseq_%s.RSEM-genes.20200506.txt' > '%s'", D, name, name, tr))
  ex <- read.delim(tr, sep = "\t", check.names = FALSE); rownames(ex) <- ex[, 1]; ex <- as.matrix(ex[, -1])
  mus(ex, paste0("cgga", name))
}
cat("DONE\n")
