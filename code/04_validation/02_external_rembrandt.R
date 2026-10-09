#!/usr/bin/env Rscript
# Stage 4 (external validation): REMBRANDT (GSE108474, GPL570).
# probe->gene via hgu133plus2.db (max probe expr per gene); samples->subjects via biospec mapping;
# clinical OS in months (converted to days, signature trained on days).
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

suppressMessages({library(survival); library(hgu133plus2.db)})
set.seed(20260909)
setwd(ROOT)
RES <- "results"; D <- DATA_ROOT
D2 <- gsub(" ", "", D)

sig <- read.delim(file.path(RES, "signature_pc_coefs.tsv"), stringsAsFactors=FALSE)
sigmap <- setNames(sig$coefficient, sig$gene)

cind <- function(risk, t, e){ n<-length(risk); c<-0; d<-0
  for(i in 1:n){ if(!e[i]) next; ti<-t[i]
    for(j in 1:n){ if(i==j) next
      if(t[j]>ti || (t[j]==ti && !e[j])){ d<-d+1; if(risk[i]>risk[j]) c<-c+1 else if(risk[i]==risk[j]) c<-c+0.5 } } }
  if(d==0) NA else c/d }

# --- expression (probe x sample) ---
cat("reading REMBRANDT expression...\n")
# header line has 541 sample names, data rows have probeID + 541 values -> use row.names=1
ex <- read.delim(paste0(D2, "/GSE108474_REMBRANDT_GeneExpression.txt.gz"), sep="\t",
                 check.names=FALSE, row.names=1)
ex <- as.matrix(ex)
cat("probes:", nrow(ex), " samples:", ncol(ex), "\n")

# --- probe -> gene (max expr probe per gene) ---
cat("mapping probes->genes (hgu133plus2)...\n")
ids <- mapIds(hgu133plus2.db, keys=rownames(ex), column="SYMBOL", keytype="PROBEID", multiVals="first")
has <- !is.na(ids)
ex <- ex[has, ]; ids <- ids[has]
# collapse: for each gene keep the probe with max mean expression
gene_order <- unique(ids)
gene_expr <- matrix(NA, nrow=length(gene_order), ncol=ncol(ex), dimnames=list(gene_order, colnames(ex)))
for (g in gene_order) {
  p <- rownames(ex)[ids == g]
  if (length(p) == 1) gene_expr[g, ] <- ex[p, ] else gene_expr[g, ] <- apply(ex[p, , drop=FALSE], 2, max)
}
cat("gene-level probes:", nrow(gene_expr), "\n")

# --- sample (biospecimen) -> subject ---
bios <- read.delim(paste0(D2, "/GSE108474_REMBRANDT_biospecimen_mapping_GEO.txt.gz"), stringsAsFactors=FALSE)
samp2subj <- setNames(bios$SUBJECT_ID, sub("^X", "", bios$BIOSPECIMEN_ID))  # sample -> subject

clin <- read.delim(paste0(D2, "/GSE108474_REMBRANDT_clinical.data.txt.gz"), stringsAsFactors=FALSE)
clin <- clin[clin$TYPE == "PATIENT", ]
clin$OS_month <- suppressWarnings(as.numeric(clin$OVERALL_SURVIVAL_MONTHS))
clin$ev <- suppressWarnings(as.numeric(clin$EVENT_OS))
cat("clinical patients:", nrow(clin), "; with OS:", sum(!is.na(clin$OS_month)), "\n")

# map each expression sample to a clinical subject (with OS)
col2subj <- samp2subj[colnames(gene_expr)]
matched <- !is.na(col2subj) & col2subj %in% clin$SUBJECT_ID
samples <- colnames(gene_expr)[matched]
subjects <- col2subj[matched]
subj_clin <- clin[match(subjects, clin$SUBJECT_ID), ]
cat("paired samples with clinical:", length(samples), "\n")

X <- gene_expr[, samples, drop=FALSE]
sig_genes <- intersect(names(sigmap), rownames(X))
cat("signature genes present:", length(sig_genes), "/", length(sigmap), "\n")
risk <- as.numeric(t(X[sig_genes, , drop=FALSE]) %*% sigmap[sig_genes])

time_days <- subj_clin$OS_month * 30.4375
event <- subj_clin$ev
keep <- !is.na(time_days) & time_days > 0 & !is.na(event)
cat("valid OS samples:", sum(keep), "; events:", sum(event[keep]), "\n")
C <- cind(risk[keep], time_days[keep], event[keep])
cat("REMBRANDT C-index:", round(C, 4), "\n")
# KM
med <- median(risk[keep]); grp <- ifelse(risk[keep] > med, "High", "Low")
lr <- survdiff(Surv(time_days[keep], event[keep]) ~ grp); p <- 1 - pchisq(lr$chisq, df=1)
cat("KM log-rank p:", format(p, digits=3), "\n")
writeLines(paste("REMBRANDT", round(C,4), sum(grp=="High"), sum(grp=="Low"), format(p,digits=4), sep="\t"),
           file.path(RES, "external_rembrandt_result.txt"))
cat("DONE\n")
