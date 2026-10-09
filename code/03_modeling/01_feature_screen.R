#!/usr/bin/env Rscript
# Stage 3 (feature screen): univariate Cox (per candidate gene) -> LASSO-Cox
# Input:  results/tcga_training_expr_log2cpm.tsv, tcga_training_clinical.tsv,
#         results/candidate_genes_single_cell.tsv
# Output: results/univariate_cox.tsv, results/lasso_fit.rds, results/lasso_genes.tsv
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

suppressMessages({library(survival); library(glmnet)})
set.seed(20260909)
setwd(ROOT)
RES <- "results"

cat("Loading training data...\n")
expr <- read.delim(file.path(RES, "tcga_training_expr_log2cpm.tsv"), row.names=1, check.names=FALSE)
clin <- read.delim(file.path(RES, "tcga_training_clinical.tsv"), stringsAsFactors=FALSE)
cand <- read.delim(file.path(RES, "candidate_genes_single_cell.tsv"), stringsAsFactors=FALSE)

# align sample order
samples <- clin$sample_submitter_id
expr <- expr[, intersect(colnames(expr), samples)]
clin <- clin[match(colnames(expr), clin$sample_submitter_id), ]
cat("Training matrix:", nrow(expr), "genes x", ncol(expr), "samples; events:", sum(clin$OS_event), "\n")

# use genes present in both candidate set and expression
genes <- intersect(cand$gene, rownames(expr))
cat("candidate genes available in expr:", length(genes), "\n")

# build survival object
surv <- Surv(time=clin$OS_time_days, event=clin$OS_event)

# --- univariate Cox per gene ---
cat("Univariate Cox on", length(genes), "genes...\n")
mat <- t(as.matrix(expr[genes, ]))
res <- vector("list", length(genes))
names(res) <- genes
for (g in genes) {
  x <- mat[, g]
  if (sd(x) == 0) next
  fit <- tryCatch(coxph(surv ~ x), error=function(e) NULL)
  if (is.null(fit)) next
  s <- summary(fit)
  res[[g]] <- data.frame(gene=g, coef=coef(fit), HR=exp(coef(fit)),
                         se=s$coefficients[,"se(coef)"], p=s$coefficients[,"Pr(>|z|)"])
}
uni <- do.call(rbind, res)
uni <- uni[order(uni$p), ]
write.table(uni, file.path(RES, "univariate_cox.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cat("univariate results:", nrow(uni), "; p<0.05:", sum(uni$p < 0.05), "\n")

# --- LASSO-Cox on genes with p<0.05 (pre-registered threshold) ---
keepg <- uni$gene[uni$p < 0.05]
if (length(keepg) < 5) { cat("Too few genes after univariate screen; widening to FDR or top\n"); keepg <- head(uni$gene, 30) }
X <- t(as.matrix(expr[keepg, ]))
cat("LASSO on", length(keepg), "genes\n")
cvfit <- cv.glmnet(X, surv, family="cox", alpha=1, nfolds=10)
saveRDS(list(cvfit=cvfit, keepg=keepg), file.path(RES, "lasso_fit.rds"))

# lambda.min vs lambda.1se
coef_min <- as.matrix(coef(cvfit, s="lambda.min"))
coef_1se <- as.matrix(coef(cvfit, s="lambda.1se"))
sel_min <- rownames(coef_min)[coef_min[,1] != 0]
sel_1se <- rownames(coef_1se)[coef_1se[,1] != 0]
cat("lambda.min selected genes:", length(sel_min), "; lambda.1se:", length(sel_1se), "\n")
res_las <- data.frame(lambda=c("min","1se"), n_genes=c(length(sel_min), length(sel_1se)),
                      genes=c(paste(sel_min, collapse=","), paste(sel_1se, collapse=",")))
write.table(res_las, file.path(RES, "lasso_genes.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cat("DONE\n")
