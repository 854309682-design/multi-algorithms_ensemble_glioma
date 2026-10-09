#!/usr/bin/env Rscript
# Stage 4 (external, 3rd cohort) -- GSE16011 (Gravendeel 2009; platform GPL8542).
#
# LIMITATION (established 2026-09-14): the GEO series matrix for GSE16011 contains NO
# overall-survival / event fields -- its four sample characteristics are tissue, gender,
# histology(+grade) and age at diagnosis -- and GSE16011 suppl holds only GSE16011_RAW.tar.
# Survival for this cohort lives in the Gravendeel 2009 J Clin Oncol supplementary tables,
# not in GEO. OS validation is therefore NOT possible from GEO alone.
#
# What this script CAN do (and does): apply the FROZEN 29-gene signature to GSE16011 probes
# (probe -> symbol via hgu133plus2.db; GPL8542 = re-annotated HG-U133 Plus 2.0) and test the
# expected grade ordering (IV > III > II) -- a cross-platform transfer / direction check,
# explicitly NOT an OS validation.
#
# Input:  <DATA_ROOT>/GSE16011_series_matrix.txt.gz, <DATA_ROOT>/GPL8542_family.soft.gz,
#         results/signature_pc_coefs.tsv
# Output: results/external_gse16011_grade_check.tsv, results/figures/figS_gse16011_grade.{pdf,jpg}
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

set.seed(20260909)
setwd(ROOT)
RES <- "results"; FIG <- "results/figures"; D <- DATA_ROOT

sig <- read.delim(file.path(RES, "signature_pc_coefs.tsv"), stringsAsFactors = FALSE)
sm <- setNames(sig$coefficient, sig$gene)

# ---- series matrix -------------------------------------------------------------------------
ln <- readLines(gzfile(file.path(D, "GSE16011_series_matrix.txt.gz")), warn = FALSE)
sq <- function(x) gsub('^"|"$', "", x)
gsm <- sq(strsplit(ln[grep("^!Sample_geo_accession", ln)[1]], "\t")[[1]])[-1]
gett <- function(key) { i <- grep(paste0("^", key, "\\t"), ln)[1]
  if (is.na(i)) return(NULL); sq(strsplit(ln[i], "\t")[[1]])[-1] }
chars <- ln[grep("^!Sample_characteristics", ln)]
ch <- t(vapply(chars, function(l) sq(strsplit(l, "\t")[[1]])[-1], character(length(gsm))))
histology <- apply(ch, 2, function(col) { h <- col[grepl("histology:", col)]
  if (length(h)) sub(".*histology:\\s*", "", h[1]) else NA_character_ })
grade <- ifelse(grepl("grade IV", histology), 4L,
                ifelse(grepl("grade III", histology), 3L,
                       ifelse(grepl("grade II", histology), 2L,
                              ifelse(grepl("grade I\\b", histology), 1L, NA_integer_))))
cat("samples:", length(gsm), "| with histology:", sum(!is.na(histology)),
    "| grade II/III/IV:", paste(table(grade), collapse = "/"), "\n")

tab_start <- grep('^"ID_REF"', ln)[1]
tab_end <- grep("^!series_matrix_table_end", ln)[1]
ex <- read.delim(text = ln[tab_start:(tab_end - 1)], row.names = 1, check.names = FALSE)
ex <- as.matrix(ex)
cat("series matrix:", nrow(ex), "probes x", ncol(ex), "samples\n")

# ---- probe -> gene symbol -------------------------------------------------------------------
# GPL8542 is the re-annotated Affymetrix HG-U133 Plus 2.0 platform whose 17,527 probe IDs are
# BrainArray-style "<entrez>_at" names (e.g. 10000_at = Entrez 10000 = AKT3) -- NOT standard Affy
# probe IDs, so hgu133plus2.db matches only ~1%. Map the numeric part as an Entrez Gene ID.
suppressMessages(library(org.Hs.eg.db))
entrez <- sub("^([0-9]+).*$", "\\1", rownames(ex))
ids <- mapIds(org.Hs.eg.db, keys = entrez, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
cat("probes mapped to a symbol:", sum(!is.na(ids)), "/", length(ids), "\n")
keep <- !is.na(ids) & ids != "" & ids != "---"
ex <- ex[keep, , drop = FALSE]; ids <- ids[keep]
go <- unique(ids); ge <- matrix(NA, length(go), ncol(ex), dimnames = list(go, colnames(ex)))
for (g in go) { p <- rownames(ex)[ids == g]
  ge[g, ] <- if (length(p) == 1) ex[p, ] else apply(ex[p, , drop = FALSE], 2, max) }
cat("gene-level matrix:", nrow(ge), "genes x", ncol(ge), "samples\n")

# ---- frozen signature score ----------------------------------------------------------------
g <- intersect(names(sm), rownames(ge))
cat("signature genes present:", length(g), "/", length(sm), "\n")
risk <- as.numeric(t(ge[g, , drop = FALSE]) %*% sm[g])

ok <- !is.na(grade) & grade %in% c(2, 3, 4)
kw <- kruskal.test(risk[ok], factor(grade[ok]))
gbm <- risk[ok][grade[ok] == 4]; lgg <- risk[ok][grade[ok] %in% c(2, 3)]
wt <- wilcox.test(gbm, lgg, alternative = "greater")
cat(sprintf("\nKruskal-Wallis across grade II/III/IV: chi2=%.1f p=%.3g\n", kw$statistic, kw$p.value))
cat(sprintf("GBM(IV) vs LGG(II-III) one-sided Wilcoxon: W=%.0f p=%.3g\n", wt$statistic, wt$p.value))
cat("median risk by grade:\n"); print(round(tapply(risk[ok], grade[ok], median), 3))

out <- data.frame(
  item = c("n_samples", "n_with_grade", "signature_genes_used", "kw_chi2", "kw_p",
           "gbm_vs_lgg_W", "gbm_vs_lgg_p", "median_risk_II", "median_risk_III", "median_risk_IV",
           "note"),
  value = c(length(gsm), sum(ok), length(g), round(kw$statistic, 1), signif(kw$p.value, 3),
            round(wt$statistic, 0), signif(wt$p.value, 3),
            round(median(risk[ok][grade[ok] == 2]), 3), round(median(risk[ok][grade[ok] == 3]), 3),
            round(median(risk[ok][grade[ok] == 4]), 3),
            "grade-direction check only; GSE16011 has no OS in GEO (survival is in the Gravendeel 2009 Cancer Research 69(23):9065 supplement)"))
write.table(out, file.path(RES, "external_gse16011_grade_check.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

draw <- function() {
  boxplot(risk[ok] ~ grade[ok], xlab = "WHO grade", ylab = "frozen signature risk score",
          main = sprintf("GSE16011 (GPL8542) grade-direction check (KW p=%.2g)", kw$p.value),
          cex.main = 0.8, col = c("grey80", "grey65", "firebrick"), outline = FALSE)
  abline(h = 0, lty = 3, col = "grey40")
}
pdf(file.path(FIG, "figS_gse16011_grade.pdf"), width = 5, height = 5); draw(); dev.off()
jpeg(file.path(FIG, "figS_gse16011_grade.jpg"), width = 1000, height = 1000, res = 220); draw(); dev.off()
cat("\nDONE -> results/external_gse16011_grade_check.tsv + figS_gse16011_grade.*\n")
