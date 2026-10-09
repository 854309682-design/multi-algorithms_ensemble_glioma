#!/usr/bin/env Rscript
# Stage 3 refinement: protein-coding-only signature.
# Restrict candidate genes to protein_coding, redo univariate Cox -> LASSO -> freeze CoxPH,
# then re-validate on CGGA-693 & CGGA-325. Saves signature_pc.rds / signature_pc_coefs.tsv.
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
RES <- "results"; D <- file.path(DATA_ROOT, "CGGA")
REF <- REF_ROOT   # dataset index/checksums/download scripts + reference gene list (shipped with the repository)

cind <- function(risk, t, e){ n<-length(risk); c<-0; d<-0
  for(i in 1:n){ if(!e[i]) next; ti<-t[i]
    for(j in 1:n){ if(i==j) next
      if(t[j]>ti || (t[j]==ti && !e[j])){ d<-d+1; if(risk[i]>risk[j]) c<-c+1 else if(risk[i]==risk[j]) c<-c+0.5 } } }
  if(d==0) NA else c/d }

expr <- read.delim(file.path(RES, "tcga_training_expr_log2cpm.tsv"), row.names=1, check.names=FALSE)
clin <- read.delim(file.path(RES, "tcga_training_clinical.tsv"), stringsAsFactors=FALSE)
cand <- read.delim(file.path(RES, "candidate_genes_single_cell.tsv"), stringsAsFactors=FALSE)
pc <- readLines(file.path(REF, "protein_coding_genes.txt"))

samples <- clin$sample_submitter_id
expr <- expr[, intersect(colnames(expr), samples)]
clin <- clin[match(colnames(expr), clin$sample_submitter_id), ]
y <- Surv(clin$OS_time_days, clin$OS_event)

genes <- intersect(intersect(cand$gene, rownames(expr)), pc)
cat("candidate genes (protein-coding) in expr:", length(genes), "\n")

# univariate Cox
mat <- t(as.matrix(expr[genes, ])); res <- list()
for (g in genes) { x <- mat[, g]; if (sd(x) == 0) next
  f <- tryCatch(coxph(y ~ x), error=function(e) NULL); if (is.null(f)) next
  res[[g]] <- summary(f)$coefficients[1, "Pr(>|z|)"] }
p <- unlist(res); keepg <- names(p)[p < 0.05]
cat("univariate p<0.05:", length(keepg), "\n")

# LASSO
X <- t(as.matrix(expr[keepg, ]))
cvf <- cv.glmnet(X, y, family="cox", alpha=1, nfolds=10)
sel <- rownames(as.matrix(coef(cvf, s="lambda.min")))[as.matrix(coef(cvf, s="lambda.min"))[,1] != 0]
sel <- intersect(sel, pc)
cat("lambda.min protein-coding signature genes:", length(sel), "\n")

# freeze
df <- data.frame(y=y, t(as.matrix(expr[sel, ])), check.names=FALSE)
f <- coxph(y ~ ., data=df)
co <- coef(f)
write.table(data.frame(gene=names(co), coefficient=co, HR=exp(co)), file.path(RES, "signature_pc_coefs.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
saveRDS(list(model=f, genes=sel), file.path(RES, "signature_pc.rds"))
cat("in-sample C:", round(cind(predict(f, type="lp"), clin$OS_time_days, clin$OS_event), 4), "\n")
sm <- co

# re-validate on CGGA
for (name in c("693","325")) {
  zf <- paste0(D, "/CGGA.mRNAseq_", name, ".RSEM-genes.20200506.txt.zip")
  tmp <- tempfile(fileext=".txt")
  system(paste0("unzip -p '", zf, "' 'CGGA.mRNAseq_", name, ".RSEM-genes.20200506.txt' > '", tmp, "'"))
  ex <- read.delim(tmp, sep="\t", check.names=FALSE); rownames(ex) <- ex[,1]; ex <- as.matrix(ex[,-1])
  cl <- read.delim(paste0(D, "/CGGA.mRNAseq_", name, "_clinical_clean.tsv"), stringsAsFactors=FALSE)
  prim <- cl[cl$PRS_type=="Primary", ]; prim$OS <- as.numeric(prim$OS); prim$Censor <- as.numeric(prim$Censor)
  clp <- prim[!is.na(prim$OS) & prim$OS>0, ]
  ss <- intersect(colnames(ex), clp$CGGA_ID); ex <- ex[, ss]; clp <- clp[match(ss, clp$CGGA_ID), ]
  Xc <- log2(as.matrix(ex)+1); g <- intersect(names(sm), rownames(Xc))
  risk <- as.numeric(t(Xc[g,,drop=FALSE]) %*% sm[g])
  cat(sprintf("  CGGA-%s: genes %d/%d  C=%.3f  (n=%d, ev=%d)\n", name, length(g), length(sm),
              cind(risk, clp$OS, clp$Censor), nrow(clp), sum(clp$Censor)))
}
cat("DONE\n")
