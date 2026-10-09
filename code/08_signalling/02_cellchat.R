#!/usr/bin/env Rscript
# P2.2 - CellChat cell-cell communication on GSE103224 (malignant states + immune/stroma).
# Focus: signalling directed TO malignant states (esp. from Immune/TAM compartment).
# Output: results/signalling/cc_lr_pairs.tsv, cc_pathways.tsv, cc_net_weight.tsv
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

suppressMessages({library(CellChat); library(Matrix)})
set.seed(20260909)
setwd(ROOT)
SC <- "results/signalling"

expr <- as.matrix(readMM(file.path(SC, "cc_expr.mtx")))     # genes x cells (log-normalized)
genes <- readLines(file.path(SC, "cc_genes.txt"))
lab <- read.delim(file.path(SC, "cc_labels.tsv"), stringsAsFactors = FALSE)
rownames(expr) <- genes; colnames(expr) <- lab$cell
keep <- lab$label != "Oligodendrocyte"                       # only 23 cells -> drop
expr <- expr[, keep]; lab <- lab[keep, ]

dbg <- CellChatDB.human$geneInfo$Symbol
keepg <- rownames(expr) %in% dbg & Matrix::rowSums(expr > 0) >= 10
expr <- expr[keepg, ]
cat("CellChat input:", nrow(expr), "genes x", ncol(expr), "cells;", length(unique(lab$label)), "labels\n")

meta <- data.frame(labels = lab$label, row.names = lab$cell)
cellchat <- createCellChat(object = expr, meta = meta, group.by = "labels")
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

lr <- subsetCommunication(cellchat)                          # significant L-R pairs
write.table(lr, file.path(SC, "cc_lr_pairs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
pw <- subsetCommunication(cellchat, slot.name = "netP")
write.table(pw, file.path(SC, "cc_pathways.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
w <- as.data.frame(as.matrix(cellchat@net$weight))
write.table(data.frame(source = rownames(w), w, check.names = FALSE), file.path(SC, "cc_net_weight.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("significant L-R pairs:", nrow(lr), "| pathways:", nrow(pw), "\n")
mal <- c("AC", "NPC", "MES", "OPC")
to_mal <- lr[lr$target %in% mal, ]
cat("interactions targeting malignant states:", nrow(to_mal), "\n")
print(head(to_mal[order(-to_mal$prob), c("source", "target", "ligand", "receptor", "pathway_name", "prob")], 15), row.names = FALSE)
cat("DONE\n")
