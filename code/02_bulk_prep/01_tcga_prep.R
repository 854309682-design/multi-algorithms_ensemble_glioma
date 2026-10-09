#!/usr/bin/env Rscript
# Stage 2 (bulk prep): TCGA-GBM/LGG training set
# - Read STAR-counts matrix (gene x sample)
# - Restrict to primary-tumor samples with valid OS (time>0, event 0/1)
# - TMM normalization (edgeR) + log2-CPM
# - Save normalized matrix + clinical (survival, age, gender, grade)
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

suppressMessages({library(edgeR); library(data.table)})
set.seed(20260909)
setwd(ROOT)

DATA <- file.path(DATA_ROOT, "TCGA")
OUT  <- file.path(ROOT, "results")

# --- counts matrix (gene x sample) ---
cat("Reading counts matrix...\n")
cnt <- fread(file.path(DATA, "TCGA_GBM_LGG_star_counts_matrix.tsv.gz"), sep="\t", data.table=FALSE)
rownames(cnt) <- cnt[,1]; cnt <- cnt[,-1]
cnt <- as.matrix(cnt)
cat("counts:", nrow(cnt), "genes x", ncol(cnt), "samples\n")

# --- sample -> case mapping + clinical ---
meta <- read.delim(file.path(DATA, "TCGA_GBM_LGG_sample_meta.tsv"), stringsAsFactors=FALSE)
clin <- read.delim(file.path(DATA, "tcga_clinical_survival.tsv"), stringsAsFactors=FALSE)
# survival per case
clin$OS_event <- as.integer(clin$OS_event)
clin <- clin[!is.na(clin$OS_event) & !is.na(clin$OS_time_days) & clin$OS_time_days > 0, ]

# keep primary-tumor samples; join to clinical by case_id
prim <- meta[meta$sample_type == "Primary Tumor", ]
m <- merge(prim, clin[, c("case_id","OS_time_days","OS_event","age_at_index","gender","tumor_grade")],
           by="case_id", all.x=TRUE)
m <- m[!is.na(m$OS_event), ]
cat("primary samples with OS:", nrow(m), "\n")

# intersect matrix columns with selected samples
keep_samples <- intersect(colnames(cnt), m$sample_submitter_id)
cnt <- cnt[, keep_samples]
m <- m[m$sample_submitter_id %in% keep_samples, ]
m <- m[match(colnames(cnt), m$sample_submitter_id), ]
cat("modeling samples:", ncol(cnt), "\n")

# --- TMM + log2-CPM ---
dge <- DGEList(counts=cnt)
dge <- calcNormFactors(dge)
x <- cpm(dge, log=TRUE)
cat("normalized matrix:", nrow(x), "x", ncol(x), "\n")

# save
write.table(x, file=file.path(OUT, "tcga_training_expr_log2cpm.tsv"), sep="\t", quote=FALSE)
write.table(m, file=file.path(OUT, "tcga_training_clinical.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cat("Saved tcga_training_expr_log2cpm.tsv + tcga_training_clinical.tsv\n")
cat("event count:", table(m$OS_event), "\n")
