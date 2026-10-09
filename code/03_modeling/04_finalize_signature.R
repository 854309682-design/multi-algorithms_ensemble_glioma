#!/usr/bin/env Rscript
# Stage 3 (finalize): freeze the prognostic signature.
# Model-selection rule (nested CV) -> BEST-SINGLE = CoxPH (stacking gain < 0.01).
# Features: lambda.min genes from a full-training LASSO. Refit CoxPH on ALL training data.
# Output (FROZEN, applied unchanged to external cohorts):
#   results/signature_model.rds, results/signature_genes_coefs.tsv, results/signature_genes.txt
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

expr <- read.delim(file.path(RES, "tcga_training_expr_log2cpm.tsv"), row.names=1, check.names=FALSE)
clin <- read.delim(file.path(RES, "tcga_training_clinical.tsv"), stringsAsFactors=FALSE)
samples <- clin$sample_submitter_id
expr <- expr[, intersect(colnames(expr), samples)]
clin  <- clin[match(colnames(expr), clin$sample_submitter_id), ]

# lambda.min genes
lasso <- readRDS(file.path(RES, "lasso_fit.rds"))
cm <- as.matrix(coef(lasso$cvfit, s="lambda.min")); genes <- rownames(cm)[cm[,1] != 0]
genes <- intersect(genes, rownames(expr))
cat("signature genes (lambda.min):", length(genes), "\n")

X <- t(as.matrix(expr[genes, ]))
y <- Surv(clin$OS_time_days, clin$OS_event)
df <- data.frame(y=y, X, check.names=FALSE)
f <- coxph(y ~ ., data=df)

# freeze
saveRDS(list(model=f, genes=genes, ref="CoxPH_lambda.min", date=Sys.time(),
             training_n=nrow(df), training_events=sum(clin$OS_event)),
        file.path(RES, "signature_model.rds"))
co <- coef(f)
write.table(data.frame(gene=names(co), coefficient=co, HR=exp(co)),
            file.path(RES, "signature_genes_coefs.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
writeLines(genes, file.path(RES, "signature_genes.txt"))

# in-sample C-index (apparent)
cind <- function(risk, t, e){ n<-length(risk); c<-0; d<-0
  for(i in 1:n){ if(!e[i]) next; ti<-t[i]
    for(j in 1:n){ if(i==j) next
      if(t[j]>ti || (t[j]==ti && !e[j])){ d<-d+1; if(risk[i]>risk[j]) c<-c+1 else if(risk[i]==risk[j]) c<-c+0.5 } } }
  if(d==0) NA else c/d }
lp <- predict(f, type="lp")
cat("Frozen CoxPH signature;\n  genes:", length(genes), "\n  in-sample C-index:", round(cind(lp, clin$OS_time_days, clin$OS_event), 4), "\n")
print(round(exp(coef(f)), 3))
cat("Saved signature_model.rds / signature_genes_coefs.tsv / signature_genes.txt\n")
cat("DONE\n")
