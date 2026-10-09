#!/usr/bin/env Rscript
# P5.2 (cross-cutting) - Platform-robust scoring. Three scorers on the SAME 29-gene signature:
#   raw   : sum(coef * expr)                                  (platform-sensitive)
#   srank : within-sample gene ranks, signed by coef           (invariant to per-sample scaling)
#   crank : within-cohort per-gene ranks -> z, coef-weighted   (invariant to per-gene monotone transform)
# Evaluated on RNA-seq (CGGA-693/325) and, when available, microarray (REMBRANDT).
# Also: a synthetic per-gene affine "platform shift" of CGGA shows raw breaks while rank survives.
# Output: results/platform_robust_scoring.tsv
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

suppressMessages({library(survival)})
set.seed(20260909)
setwd(ROOT)
RES <- "results"; D <- DATA_ROOT
sig <- read.delim(file.path(RES, "signature_pc_coefs.tsv"), stringsAsFactors = FALSE)
co <- setNames(sig$coefficient, sig$gene)

cind <- function(risk, t, e) { n <- length(risk); c <- 0; d <- 0
  for (i in seq_len(n)) { if (!e[i]) next; ti <- t[i]
    for (j in seq_len(n)) { if (i == j) next
      if (t[j] > ti || (t[j] == ti && !e[j])) { d <- d + 1
        if (risk[i] > risk[j]) c <- c + 1 else if (risk[i] == risk[j]) c <- c + 0.5 } } }
  if (d == 0) NA else c / d }

scorers <- function(X) {                      # X = genes x samples
  g <- intersect(names(co), rownames(X)); X <- as.matrix(X[g, , drop = FALSE]); w <- co[g]
  raw <- as.numeric(t(X) %*% w)
  R <- apply(X, 2, function(v) rank(v) / length(v))
  srank <- as.numeric(t(R) %*% sign(w))
  Z <- t(apply(X, 1, function(v) { r <- rank(v); (r - mean(r)) / sd(r) }))
  crank <- as.numeric(t(Z) %*% w)
  list(raw = raw, srank = srank, crank = crank, ng = length(g))
}

rows <- list()
# --- RNA-seq (CGGA) ---
cgga <- list()
for (name in c("693", "325")) {
  tc <- tempfile(); system(sprintf("unzip -p '%s/CGGA/CGGA.mRNAseq_%s_clinical.zip' 'CGGA.mRNAseq_%s_clinical.20200506.txt' > '%s'", D, name, name, tc))
  cl <- read.delim(tc, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE); names(cl)[1] <- "CGGA_ID"
  tr <- tempfile(); system(sprintf("unzip -p '%s/CGGA/CGGA.mRNAseq_%s.RSEM-genes.zip' 'CGGA.mRNAseq_%s.RSEM-genes.20200506.txt' > '%s'", D, name, name, tr))
  ex <- read.delim(tr, sep = "\t", check.names = FALSE); rownames(ex) <- ex[, 1]; ex <- log2(as.matrix(ex[, -1]) + 1)
  cl <- cl[trimws(cl$PRS_type) == "Primary", ]; cl$OS <- suppressWarnings(as.numeric(cl$OS)); cl$Censor <- suppressWarnings(as.numeric(cl$Censor))
  cl <- cl[!is.na(cl$OS) & cl$OS > 0 & !is.na(cl$Censor), ]
  ss <- intersect(colnames(ex), trimws(cl$CGGA_ID)); ex <- ex[, ss]; cl <- cl[match(ss, trimws(cl$CGGA_ID)), ]
  s <- scorers(ex); cgga[[name]] <- list(ex = ex, t = cl$OS, e = cl$Censor)
  rows[[paste0("CGGA-", name)]] <- data.frame(cohort = paste0("CGGA-", name), platform = "RNA-seq", setting = "original",
    n = nrow(cl), events = sum(cl$Censor), C_raw = round(cind(s$raw, cl$OS, cl$Censor), 3),
    C_srank = round(cind(s$srank, cl$OS, cl$Censor), 3), C_crank = round(cind(s$crank, cl$OS, cl$Censor), 3), genes = s$ng)
}
# --- synthetic per-SAMPLE affine "platform/normalisation shift" (CGGA-693) ---
# per-sample gain+offset breaks cross-sample comparability of the raw weighted sum,
# but leaves within-sample gene ranks (srank) invariant.
ex <- cgga[["693"]]$ex
k <- ncol(ex); a <- runif(k, 0.2, 5); b <- runif(k, -3, 3)
exp_shift <- sweep(ex, 2, a, "*"); exp_shift <- sweep(exp_shift, 2, b, "+")
s <- scorers(exp_shift)
rows[["CGGA-693 platform-shifted"]] <- data.frame(cohort = "CGGA-693", platform = "RNA-seq (synthetic shift)", setting = "per-sample affine",
  n = length(cgga[["693"]]$t), events = sum(cgga[["693"]]$e), C_raw = round(cind(s$raw, cgga[["693"]]$t, cgga[["693"]]$e), 3),
  C_srank = round(cind(s$srank, cgga[["693"]]$t, cgga[["693"]]$e), 3), C_crank = round(cind(s$crank, cgga[["693"]]$t, cgga[["693"]]$e), 3), genes = s$ng)

# --- microarray (REMBRANDT) if the raw files are present ---
rmb <- file.path(D, "GSE108474_REMBRANDT_GeneExpression.txt.gz")
if (file.exists(rmb)) {
  suppressMessages(library(hgu133plus2.db))
  exm <- as.matrix(read.delim(rmb, sep = "\t", check.names = FALSE, row.names = 1))
  ids <- mapIds(hgu133plus2.db, keys = rownames(exm), column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
  exm <- exm[!is.na(ids), ]; ids <- ids[!is.na(ids)]
  go <- unique(ids); ge <- matrix(NA, length(go), ncol(exm), dimnames = list(go, colnames(exm)))
  for (g in go) { p <- rownames(exm)[ids == g]; ge[g, ] <- if (length(p) == 1) exm[p, ] else apply(exm[p, , drop = FALSE], 2, max) }
  bios <- read.delim(file.path(D, "GSE108474_REMBRANDT_biospecimen_mapping_GEO.txt.gz"), stringsAsFactors = FALSE)
  s2s <- setNames(bios$SUBJECT_ID, sub("^X", "", bios$BIOSPECIMEN_ID))
  cl <- read.delim(file.path(D, "GSE108474_REMBRANDT_clinical.data.txt.gz"), stringsAsFactors = FALSE)
  cl <- cl[cl$TYPE == "PATIENT", ]; cl$OS_month <- suppressWarnings(as.numeric(cl$OVERALL_SURVIVAL_MONTHS)); cl$ev <- suppressWarnings(as.numeric(cl$EVENT_OS))
  c2s <- s2s[colnames(ge)]; m <- !is.na(c2s) & c2s %in% cl$SUBJECT_ID
  samp <- colnames(ge)[m]; sc <- cl[match(c2s[m], cl$SUBJECT_ID), ]
  td <- sc$OS_month * 30.4375; ev <- sc$ev; k <- !is.na(td) & td > 0 & !is.na(ev)
  s <- scorers(ge[, samp, drop = FALSE])
  rows[["REMBRANDT"]] <- data.frame(cohort = "REMBRANDT", platform = "microarray", setting = "original", n = sum(k), events = sum(ev[k]),
    C_raw = round(cind(s$raw[k], td[k], ev[k]), 3), C_srank = round(cind(s$srank[k], td[k], ev[k]), 3),
    C_crank = round(cind(s$crank[k], td[k], ev[k]), 3), genes = s$ng)
} else cat("NOTE: REMBRANDT raw files not found at", rmb, "- microarray arm skipped (see provisioning list).\n")

out <- do.call(rbind, rows)
write.table(out, file.path(RES, "platform_robust_scoring.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nDONE\n"); print(out, row.names = FALSE)
