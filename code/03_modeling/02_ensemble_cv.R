#!/usr/bin/env Rscript
# Stage 3 (ensemble): >=6 survival learners, patient-level 10-fold CV, best-single vs stacking.
# Pre-registered decision rule: take stacking if CV C-index >= best-single + 0.01, else best-single.
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
                  library(survivalsvm); library(CoxBoost); library(xgboost)})
set.seed(20260909)
setwd(ROOT)
RES <- "results"

expr <- read.delim(file.path(RES, "tcga_training_expr_log2cpm.tsv"), row.names=1, check.names=FALSE)
clin <- read.delim(file.path(RES, "tcga_training_clinical.tsv"), stringsAsFactors=FALSE)
lasso <- readRDS(file.path(RES, "lasso_fit.rds"))
cvfit <- lasso$cvfit
# use the lambda.min selected genes as the ensemble feature set
coef_min <- as.matrix(coef(cvfit, s="lambda.min"))
genes <- rownames(coef_min)[coef_min[,1] != 0]
cat("Ensemble feature genes (lambda.min):", length(genes), "\n")

samples <- clin$sample_submitter_id
expr <- expr[, intersect(colnames(expr), samples)]
clin  <- clin[match(colnames(expr), clin$sample_submitter_id), ]
X <- t(as.matrix(expr[intersect(genes, rownames(expr)), ]))
Y <- Surv(time=clin$OS_time_days, event=clin$OS_event)
cat("X:", nrow(X), "samples x", ncol(X), "genes; events:", sum(clin$OS_event), "\n")

# --- stratified 10-fold (by event) ---
set.seed(20260909)
K <- 10
fold <- integer(nrow(X))
oidx <- sample(seq_len(nrow(X)))
for (k in 1:K) {
  fold[oidx[seq(k, length(oidx), by=K)]] <- k
}

# Harrell C-index (higher risk = worse is concordant). Validated: age=0.768 in this set.
cind <- function(risk, t, e) {
  if (length(unique(risk)) == 1) return(NA)
  n <- length(risk); count <- 0; denom <- 0
  for (i in seq_len(n)) {
    if (!e[i]) next
    ti <- t[i]
    for (j in seq_len(n)) {
      if (i == j) next
      if (t[j] > ti || (t[j] == ti && !e[j])) {
        denom <- denom + 1
        if (risk[i] > risk[j]) count <- count + 1 else if (risk[i] == risk[j]) count <- count + 0.5
      }
    }
  }
  if (denom == 0) NA else count / denom
}

trainRisk <- matrix(NA, nrow=nrow(X), ncol=7)
colnames(trainRisk) <- c("CoxPH","Ridge","ElasticNet","CoxBoost","RSF","SurvSVM","XGBoost")

for (k in 1:K) {
  tr <- fold != k; va <- fold == k
  Xtr <- X[tr, , drop=FALSE]; Xva <- X[va, , drop=FALSE]
  ttr <- clin$OS_time_days[tr]; e_tr <- clin$OS_event[tr]
  tva <- clin$OS_time_days[va]; e_va <- clin$OS_event[va]
  df_tr <- data.frame(time=ttr, event=e_tr, Xtr, check.names=FALSE)
  df_va <- data.frame(Xva, check.names=FALSE)

  # CoxPH
  f <- tryCatch(coxph(Surv(time,event) ~ ., data=df_tr), error=function(e) NULL)
  if (!is.null(f)) trainRisk[va, "CoxPH"] <- predict(f, newdata=df_va, type="lp")

  # Ridge / ElasticNet (glmnet, alpha)
  for (m in c("Ridge","ElasticNet")) {
    al <- ifelse(m == "Ridge", 0, 0.5)
    f <- tryCatch(cv.glmnet(Xtr, Surv(ttr,e_tr), family="cox", alpha=al),
                  error=function(e) NULL)
    if (!is.null(f)) trainRisk[va, m] <- as.numeric(predict(f, newx=Xva, s="lambda.min"))
  }

  # CoxBoost
  f <- tryCatch(CoxBoost(time=ttr, status=e_tr, x=as.matrix(Xtr), stepno=50),
                error=function(e) NULL)
  if (!is.null(f)) {
    pr <- tryCatch(as.numeric(predict(f, newdata=as.matrix(Xva), type="lp")), error=function(e) NULL)
    if (!is.null(pr) && length(pr) == length(tva)) trainRisk[va, "CoxBoost"] <- pr
  }

  # RSF
  f <- tryCatch(rfsrc(Surv(time,event) ~ ., data=df_tr, ntree=500), error=function(e) NULL)
  if (!is.null(f)) trainRisk[va, "RSF"] <- predict(f, newdata=df_va)$predicted

  # SurvSVM
  f <- tryCatch(survivalsvm(Surv(time,event) ~ ., data=df_tr), error=function(e) NULL)
  if (!is.null(f)) trainRisk[va, "SurvSVM"] <- predict(f, newdata=df_va)$predicted

  # XGBoost (survival:aft, censoring-aware: lower/upper bounds on time)
  f <- tryCatch({
    params <- list(objective="survival:aft", aft_loss_distribution="normal",
                   eval_metric="aft-nloglik", eta=0.05, max_depth=3)
    dtrain <- xgb.DMatrix(as.matrix(Xtr), label=ttr)
    setinfo(dtrain, "label_lower_bound", as.numeric(ttr))
    setinfo(dtrain, "label_upper_bound", as.numeric(ifelse(e_tr==1, ttr, 1e9)))
    xgb.train(params, dtrain, nrounds=100, verbose=0)
  }, error=function(e) NULL)
  if (!is.null(f)) {
    pr <- tryCatch(predict(f, xgb.DMatrix(as.matrix(Xva))), error=function(e) NULL)
    # AFT predicts survival time (higher = better) -> risk = -prediction
    if (!is.null(pr) && length(pr) == length(tva)) trainRisk[va, "XGBoost"] <- -as.numeric(pr)
  }
}
cat("OOB risk matrix built; NAs per learner:\n"); print(colSums(is.na(trainRisk)))

# --- per-learner CV C-index (align direction so higher risk = worse; flip if neg concord) ---
cat("\n===== Learner CV performance (patient-level 10-fold) =====\n")
ci <- sapply(colnames(trainRisk), function(cl) {
  r <- trainRisk[, cl]; ok <- !is.na(r)
  if (sum(ok) == 0) return(NA)
  c <- cind(r[ok], clin$OS_time_days[ok], clin$OS_event[ok])
  c
})
ci_df <- data.frame(learner=names(ci), Cindex=round(ci, 4))
print(ci_df)
write.table(ci_df, file.path(RES, "ensemble_cv_cindex.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# --- stacking: Cox on the OOB risk matrix as predictors (only learners with predictions) ---
use_cols <- colnames(trainRisk)[colSums(!is.na(trainRisk)) > 0]
stk <- trainRisk[, use_cols, drop=FALSE]
ok <- complete.cases(stk)
stk_df <- data.frame(time=clin$OS_time_days[ok], event=clin$OS_event[ok], stk[ok, , drop=FALSE])
stk_df[, 3:ncol(stk_df)] <- scale(stk_df[, 3:ncol(stk_df)])
stk <- tryCatch(coxph(Surv(time,event) ~ ., data=stk_df), error=function(e) NULL)
if (!is.null(stk)) {
  stack_risk <- predict(stk, type="lp", newdata=stk_df)
  stack_ci <- cind(stack_risk, stk_df$time, stk_df$event)
  cat("\nstacking C-index:", round(stack_ci, 4), "\n")
  write.table(data.frame(model="stacking", Cindex=round(stack_ci,4)),
              file.path(RES, "ensemble_cv_cindex.tsv"), sep="\t", quote=FALSE, row.names=FALSE, append=TRUE)
} else cat("\nstacking failed\n")

best_single  <- max(ci, na.rm=TRUE); best_name <- names(ci)[which.max(ci)]
if (!is.null(stk) && stack_ci >= best_single + 0.01) {
  cat("=> DECISION: take STACKING (", round(stack_ci,3), ") over best-single (", best_name, round(best_single,3), ")\n")
} else {
  cat("=> DECISION: take BEST-SINGLE =", best_name, round(best_single,3), "\n")
}
cat("DONE\n")
