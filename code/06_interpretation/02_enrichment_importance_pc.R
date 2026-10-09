#!/usr/bin/env Rscript
# Stage 6 (biology closure) -- regenerated on the PRIMARY 29-gene protein-coding signature.
# SHAP-style importance (coef x SD), single-cell state traceback, GO/KEGG enrichment.
# Overwrites signature_importance_states.tsv / signature_GO_BP.tsv / signature_KEGG.tsv;
# the earlier 34-gene versions are preserved in results/archive/*_34gene.tsv.
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

suppressMessages({library(clusterProfiler); library(org.Hs.eg.db)})
options(timeout = 30)   # KEGG REST may be unreachable offline; fail fast instead of hanging
set.seed(20260909)
setwd(ROOT)
RES <- "results"

expr <- read.delim(file.path(RES, "tcga_training_expr_log2cpm.tsv"), row.names = 1, check.names = FALSE)
sig  <- read.delim(file.path(RES, "signature_pc_coefs.tsv"), stringsAsFactors = FALSE)
cand <- read.delim(file.path(RES, "candidate_genes_single_cell.tsv"), stringsAsFactors = FALSE)
sm <- setNames(sig$coefficient, sig$gene)
cat("primary signature genes:", length(sm), "(protein-coding)\n")

# --- SHAP-style importance: |coef| * SD(gene) ---
genes <- intersect(names(sm), rownames(expr))
X <- as.matrix(expr[genes, ])
imp <- data.frame(gene = genes, coef = sm[genes], sd = apply(X, 1, sd),
                  importance = abs(sm[genes]) * apply(X, 1, sd),
                  state = cand$state[match(genes, cand$gene)], stringsAsFactors = FALSE)
imp <- imp[order(-imp$importance), ]
write.table(imp, file.path(RES, "signature_importance_states.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("signature gene importance (top 10):\n"); print(imp[1:10, c("gene", "importance", "state")])

# --- GO BP / KEGG enrichment ---
ids <- mapIds(org.Hs.eg.db, keys = imp$gene, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
ids <- ids[!is.na(ids)]
go <- enrichGO(unique(ids), OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
               pvalueCutoff = 0.1, qvalueCutoff = 0.2, readable = TRUE)
if (!is.null(go) && nrow(as.data.frame(go)) > 0) {
  write.table(as.data.frame(go), file.path(RES, "signature_GO_BP.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\nGO BP terms (top):\n"); print(head(as.data.frame(go)[, c("Description", "p.adjust", "Count")], 8))
}
kegg <- tryCatch(enrichKEGG(ids, organism = "hsa", pvalueCutoff = 0.1), error = function(e) NULL)
if (!is.null(kegg) && nrow(as.data.frame(kegg)) > 0) {
  write.table(as.data.frame(kegg), file.path(RES, "signature_KEGG.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\nKEGG terms (top):\n"); print(head(as.data.frame(kegg)[, c("Description", "p.adjust", "Count")], 8))
}
cat("DONE\n")
