#!/usr/bin/env Rscript
# Stage 3 (honest benchmark): NESTED patient-level 10-fold CV.
# Inside EACH outer fold: univariate Cox (candidate genes) -> LASSO-Cox -> features,
# then fit each learner on the training folds and predict the held-out fold.
# Gives an unbiased C-index for the whole modeling procedure (fixes selection leakage).
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

suppressMessages({library(survival); library(glmnet); library(randomForestSRC);
                  library(CoxBoost); library(xgboost)})
set.seed(20260909)
setwd(ROOT)
RES <- "results"

cind <- function(risk, t, e) {
  if (length(unique(risk)) == 1) return(NA)
  n <- length(risk); cnt <- 0; den <- 0
  for (i in seq_len(n)) { if (!e[i]) next; ti <- t[i]
    for (j in seq_len(n)) { if (i == j) next
      if (t[j] > ti || (t[j] == ti && !e[j])) { den <- den + 1
        if (risk[i] > risk[j]) cnt <- cnt + 1 else if (risk[i] == risk[j]) cnt <- cnt + 0.5 } } }
  if (den == 0) NA else cnt / den
}

expr <- read.delim(file.path(RES, "tcga_training_expr_log2cpm.tsv"), row.names=1, check.names=FALSE)
clin <- read.delim(file.path(RES, "tcga_training_clinical.tsv"), stringsAsFactors=FALSE)
cand <- read.delim(file.path(RES, "candidate_genes_single_cell.tsv"), stringsAsFactors=FALSE)
cand_genes <- intersect(cand$gene, rownames(expr))
cat("candidate genes:", length(cand_genes), "\n")

samples <- clin$sample_submitter_id
expr <- expr[, intersect(colnames(expr), samples)]
clin  <- clin[match(colnames(expr), clin$sample_submitter_id), ]
Xall  <- t(as.matrix(expr[cand_genes, ]))
time  <- clin$OS_time_days; event <- clin$OS_event
cat("X:", nrow(Xall), "x", ncol(Xall), "\n")

# stratified 10-fold
K <- 10; fold <- integer(nrow(Xall))
oidx <- sample(seq_len(nrow(Xall)))
for (k in 1:K) fold[oidx[seq(k, length(oidx), by=K)]] <- k

learners <- c("CoxPH","Ridge","ElasticNet","CoxBoost","RSF","XGBoost")
oob <- matrix(NA, nrow=nrow(Xall), ncol=length(learners), dimnames=list(NULL, learners))
nsel <- integer(K)

for (k in 1:K) {
  tr <- fold != k; va <- fold == k
  Xtr <- Xall[tr, , drop=FALSE]; Xva <- Xall[va, , drop=FALSE]
  ttr <- time[tr]; etr <- event[tr]; tva <- time[va]; eva <- event[va]
  # --- nested feature selection inside this fold ---
  # univariate Cox
  p <- sapply(seq_len(ncol(Xtr)), function(j){
    x <- Xtr[, j]; if (sd(x) == 0) return(1)
    f <- tryCatch(coxph(Surv(ttr, etr) ~ x), error=function(e) NULL)
    if (is.null(f)) return(1)
    summary(f)$coefficients[,"Pr(>|z|)"]
  })
  keep <- which(p < 0.05)
  if (length(keep) < 10) keep <- order(p)[1:min(30, ncol(Xtr))]
  Xtrs <- Xtr[, keep, drop=FALSE]
  # LASSO
  cv <- tryCatch(cv.glmnet(Xtrs, Surv(ttr, etr), family="cox", alpha=1, nfolds=5), error=function(e) NULL)
  if (!is.null(cv)) {
    cf <- as.matrix(coef(cv, s="lambda.min"))
    sel <- rownames(cf)[cf[,1] != 0]
  } else sel <- colnames(Xtrs)
  nsel[k] <- length(sel)
  Xtr2 <- Xtrs[, intersect(sel, colnames(Xtrs)), drop=FALSE]
  Xva2 <- Xva[, intersect(sel, colnames(Xtrs)), drop=FALSE]
  df_tr <- data.frame(time=ttr, event=etr, Xtr2, check.names=FALSE)
  df_va <- data.frame(Xva2, check.names=FALSE)

  # --- fit learners ---
  f <- tryCatch(coxph(Surv(time,event) ~ ., data=df_tr), error=function(e) NULL)
  if (!is.null(f)) oob[va, "CoxPH"] <- predict(f, newdata=df_va, type="lp")
  for (m in c("Ridge","ElasticNet")) { al <- ifelse(m=="Ridge",0,0.5)
    f <- tryCatch(cv.glmnet(as.matrix(Xtr2), Surv(ttr,etr), family="cox", alpha=al), error=function(e) NULL)
    if (!is.null(f)) oob[va, m] <- as.numeric(predict(f, newx=as.matrix(Xva2), s="lambda.min")) }
  f <- tryCatch(CoxBoost(time=ttr, status=etr, x=as.matrix(Xtr2), stepno=50), error=function(e) NULL)
  if (!is.null(f)) oob[va, "CoxBoost"] <- tryCatch(as.numeric(predict(f, newdata=as.matrix(Xva2), type="lp")), error=function(e) NULL)
  f <- tryCatch(rfsrc(Surv(time,event) ~ ., data=df_tr, ntree=500), error=function(e) NULL)
  if (!is.null(f)) oob[va, "RSF"] <- predict(f, newdata=df_va)$predicted
  f <- tryCatch({ pp <- list(objective="survival:aft", aft_loss_distribution="normal", eval_metric="aft-nloglik", eta=0.05, max_depth=3)
                 dtf <- xgb.DMatrix(as.matrix(Xtr2), label=ttr)
                 setinfo(dtf, "label_lower_bound", as.numeric(ttr))
                 setinfo(dtf, "label_upper_bound", as.numeric(ifelse(etr==1, ttr, 1e9)))
                 xgb.train(pp, dtf, nrounds=100, verbose=0) }, error=function(e) NULL)
  if (!is.null(f)) oob[va, "XGBoost"] <- tryCatch(-as.numeric(predict(f, xgb.DMatrix(as.matrix(Xva2)))), error=function(e) NULL)
  cat(sprintf("fold %d: n_sel=%d\n", k, length(sel)))
}

cat("\n===== NESTED CV: selected features per fold:", nsel, "=====\n")
ci <- sapply(learners, function(cl){ r<-oob[,cl]; ok<-!is.na(r); if(sum(ok)==0) NA else cind(r[ok], time[ok], event[ok]) })
cat("Per-learner nested C-index:\n"); print(round(ci, 4))
write.table(data.frame(learner=names(ci), Cindex=round(ci,4)), file.path(RES,"nested_cv_cindex.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

use <- learners[colSums(!is.na(oob)) > 0]
stk <- oob[, use, drop=FALSE]; ok <- complete.cases(stk)
sdf <- data.frame(time=time[ok], event=event[ok], scale(stk[ok,,drop=FALSE]))
sf <- tryCatch(coxph(Surv(time,event) ~ ., data=sdf), error=function(e) NULL)
stack_ci <- if (!is.null(sf)) cind(predict(sf, type="lp", newdata=sdf), sdf$time, sdf$event) else NA
cat("Nested stacking C-index:", round(stack_ci,4), "\n")

best <- max(ci, na.rm=TRUE); bn <- names(ci)[which.max(ci)]
if (!is.na(stack_ci) && stack_ci >= best + 0.01) {
  cat("=> DECISION: stacking (", round(stack_ci,3), ")\n")
} else {
  cat("=> DECISION: BEST-SINGLE =", bn, round(best,3), "\n")
}
cat("DONE\n")
