#!/usr/bin/env Rscript
# Stage B3: benchmark >=3 PUBLISHED glioma signatures against the single-cell-anchored
# signature on the same external cohorts (CGGA-693 / CGGA-325).
#
# PUBLISHED SIGNATURES ARE EXTERNAL INPUT -- this script does NOT hard-code any gene list
# (that would risk transcribing a published list incorrectly). Drop curated tables into:
#     $GLIOMA_ROOT/reference/published_signatures/<study>.tsv
# with a header:  gene <tab> weight    (weight = published coefficient/HR-log; optional --
# if absent, a unit-weight z-score is used).  The `study` name is taken from the file name.
#
# Output: results/published_signature_comparison.tsv  (C-index per study per cohort)
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

suppressMessages(library(survival))
setwd(ROOT)
RES <- "results"; D <- file.path(DATA_ROOT, "CGGA")
PUB <- file.path(REF_ROOT, "published_signatures")

cind <- function(risk, t, e) { ok <- !is.na(risk) & !is.na(t) & !is.na(e); risk <- risk[ok]; t <- t[ok]; e <- e[ok]
  n <- length(risk); if (n < 10) return(NA); c <- 0; d <- 0
  for (i in seq_len(n)) { if (!e[i]) next; ti <- t[i]
    for (j in seq_len(n)) { if (i == j) next
      if (t[j] > ti || (t[j] == ti && !e[j])) { d <- d + 1
        if (risk[i] > risk[j]) c <- c + 1 else if (risk[i] == risk[j]) c <- c + 0.5 } } }
  if (d == 0) NA else c / d }

# our own signature for the comparison row
own <- read.delim(file.path(RES, "signature_pc_coefs.tsv"), stringsAsFactors = FALSE)
sigs <- list(`single-cell-anchored (this study)` = setNames(own$coefficient, own$gene))
if (dir.exists(PUB)) for (f in list.files(PUB, pattern = "\\.tsv$")) {
  t <- read.delim(file.path(PUB, f), stringsAsFactors = FALSE)
  w <- if ("weight" %in% names(t)) t$weight else rep(1, nrow(t))
  sigs[[sub("\\.tsv$", "", f)]] <- setNames(w, t$gene)
}
if (length(sigs) < 2) cat("!! No published signature tables found in", PUB, "-- only the in-house signature will be scored.\n")

# load cohorts once
cohorts <- list()
for (name in c("693", "325")) {
  tr <- tempfile(fileext = ".txt")
  system(sprintf("unzip -p '%s/CGGA.mRNAseq_%s.RSEM-genes.zip' 'CGGA.mRNAseq_%s.RSEM-genes.20200506.txt' > '%s'", D, name, name, tr))
  ex <- read.delim(tr, sep = "\t", check.names = FALSE); rownames(ex) <- ex[, 1]; ex <- as.matrix(ex[, -1])
  tc <- tempfile(fileext = ".txt")
  system(sprintf("unzip -p '%s/CGGA.mRNAseq_%s_clinical.zip' 'CGGA.mRNAseq_%s_clinical.20200506.txt' > '%s'", D, name, name, tc))
  cl <- read.delim(tc, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE); names(cl)[1] <- "CGGA_ID"
  cl <- cl[trimws(cl$PRS_type) == "Primary", ]
  cl$OS <- suppressWarnings(as.numeric(cl$OS)); cl$Censor <- suppressWarnings(as.numeric(cl$Censor))
  cl <- cl[!is.na(cl$OS) & cl$OS > 0 & !is.na(cl$Censor), ]
  ss <- intersect(colnames(ex), trimws(cl$CGGA_ID)); cohorts[[name]] <- list(X = log2(ex[, ss, drop = FALSE] + 1), cl = cl[match(ss, trimws(cl$CGGA_ID)), ])
}

rows <- list()
for (s in names(sigs)) {
  for (name in names(cohorts)) {
    X <- cohorts[[name]]$X; cl <- cohorts[[name]]$cl
    g <- intersect(names(sigs[[s]]), rownames(X)); if (length(g) < 3) next
    z <- t(scale(t(X[g, , drop = FALSE])))                 # per-gene z within cohort
    risk <- as.numeric(t(z) %*% sigs[[s]][g] / sum(abs(sigs[[s]][g])))
    rows[[paste(s, name)]] <- data.frame(study = s, cohort = paste0("CGGA-", name), genes_used = length(g),
                                          genes_total = length(sigs[[s]]), C_index = round(cind(risk, cl$OS, cl$Censor), 3))
  }
}
out <- do.call(rbind, rows)
write.table(out, file.path(RES, "published_signature_comparison.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nDONE -> results/published_signature_comparison.tsv\n"); print(out, row.names = FALSE)
