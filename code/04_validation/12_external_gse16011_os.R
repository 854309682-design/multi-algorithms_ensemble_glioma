#!/usr/bin/env Rscript
# Stage 4 (external, 3rd cohort) -- GSE16011 OS arm (Gravendeel 2009, *Cancer Research*; GPL8542).
#
# STATUS (2026-09-14): the GSE16011 OVERALL-SURVIVAL table is NOT obtainable in this environment.
#   - GEO GSE16011 carries no OS/event fields (verified: only 4 sample characteristics --
#     tissue / gender / histology(+grade) / age) and its suppl holds only GSE16011_RAW.tar.
#   - Survival lives in the Gravendeel 2009 paper's supplementary (Cancer Research
#     69(23):9065-9072; DOI 10.1158/0008-5472.CAN-09-2307) -- CLOSED access (Unpaywall
#     oa_status = closed, no PMC/repository copy), and the publisher page is Cloudflare-blocked.
#   See data_acquisition/README.md section 7 for the access routes.
#
# This script is the READY-TO-RUN OS arm: it auto-skips (exit 0, nothing written) while the table
# is absent and runs the full validation the moment the file lands.
#
# Expected input (TSV, curated by the user from the paper's supplementary):
#   <DATA_ROOT>/GSE16011_clinical_OS.tsv   columns:
#     sample      GSM id as in the series matrix (e.g. GSM405200)
#     os_months   overall-survival time in months (numeric)
#     event       1 = death, 0 = censored
#     (optional)  grade   WHO grade (II/III/IV) -- enables the grade-adjusted model
#
# Input:  <DATA_ROOT>/GSE16011_clinical_OS.tsv, <DATA_ROOT>/GSE16011_series_matrix.txt.gz,
#         results/signature_pc_coefs.tsv
# Output: results/external_gse16011_os.tsv, results/figures/figS_gse16011_km.{pdf,jpg}
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
OSF <- file.path(D, "GSE16011_clinical_OS.tsv")

if (!file.exists(OSF)) {
  cat("SKIP: GSE16011 OS table not present at\n  ", OSF, "\n")
  cat("      GSE16011 has no OS in GEO; survival is in the Gravendeel 2009 Cancer Res supplement.\n")
  cat("      Provision it (see data_acquisition/README.md section 7) and re-run. Nothing written.\n")
  quit(status = 0)
}

# Harrell C (manual -- survival::concordance has the opposite sign convention)
cind <- function(risk, t, e) { ok <- !is.na(risk) & !is.na(t) & !is.na(e); risk <- risk[ok]; t <- t[ok]; e <- e[ok]
  n <- length(risk); if (n < 10) return(NA); c <- 0; d <- 0
  for (i in seq_len(n)) { if (!e[i]) next; ti <- t[i]
    for (j in seq_len(n)) { if (i == j) next
      if (t[j] > ti || (t[j] == ti && !e[j])) { d <- d + 1
        if (risk[i] > risk[j]) c <- c + 1 else if (risk[i] == risk[j]) c <- c + 0.5 } } }
  if (d == 0) NA else c / d }

sig <- read.delim(file.path(RES, "signature_pc_coefs.tsv"), stringsAsFactors = FALSE)
sm <- setNames(sig$coefficient, sig$gene)

# ---- OS table -------------------------------------------------------------------------------
os <- read.delim(OSF, stringsAsFactors = FALSE, check.names = FALSE)
need <- c("sample", "os_months", "event")
if (!all(need %in% names(os))) stop("GSE16011 OS table must have columns: ", paste(need, collapse = ", "))
os$os_months <- suppressWarnings(as.numeric(os$os_months))
os$event <- suppressWarnings(as.numeric(os$event))
os <- os[!is.na(os$os_months) & os$os_months > 0 & !is.na(os$event), , drop = FALSE]
cat("OS table:", nrow(os), "samples with survival (", sum(os$event == 1), "events )\n")

# ---- series matrix -> gene-level matrix (same probe mapping as 11_external_gse16011_grade_check.R)
ln <- readLines(gzfile(file.path(D, "GSE16011_series_matrix.txt.gz")), warn = FALSE)
sq <- function(x) gsub('^"|"$', "", x)
gsm <- sq(strsplit(ln[grep("^!Sample_geo_accession", ln)[1]], "\t")[[1]])[-1]
tab_start <- grep('^"ID_REF"', ln)[1]; tab_end <- grep("^!series_matrix_table_end", ln)[1]
ex <- as.matrix(read.delim(text = ln[tab_start:(tab_end - 1)], row.names = 1, check.names = FALSE))
# GPL8542 probe IDs are BrainArray-style "<entrez>_at" -> map the numeric part as an Entrez ID
suppressMessages(library(org.Hs.eg.db))
entrez <- sub("^([0-9]+).*$", "\\1", rownames(ex))
ids <- mapIds(org.Hs.eg.db, keys = entrez, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
keep <- !is.na(ids) & ids != "" & ids != "---"
ex <- ex[keep, , drop = FALSE]; ids <- ids[keep]
go <- unique(ids); ge <- matrix(NA, length(go), ncol(ex), dimnames = list(go, colnames(ex)))
for (g in go) { p <- rownames(ex)[ids == g]
  ge[g, ] <- if (length(p) == 1) ex[p, ] else apply(ex[p, , drop = FALSE], 2, max) }
cat("gene-level matrix:", nrow(ge), "genes x", ncol(ge), "samples\n")

# ---- frozen score, matched to survival -------------------------------------------------------
g <- intersect(names(sm), rownames(ge))
cat("signature genes present:", length(g), "/", length(sm), "\n")
risk <- as.numeric(t(ge[g, , drop = FALSE]) %*% sm[g]); names(risk) <- colnames(ge)
ss <- intersect(gsm, os$sample)
if (length(ss) < 20) stop("only ", length(ss), " samples overlap between the series matrix and the OS table")
r <- risk[ss]; o <- os[match(ss, os$sample), ]
C <- cind(r, o$os_months, o$event)
sd_ <- survdiff(Surv(os_months, event) ~ I(r > median(r)), data = data.frame(os_months = o$os_months, event = o$event, r = r))
lr_p <- 1 - pchisq(sd_$chisq, length(sd_$n) - 1)
cat(sprintf("\nGSE16011 OS: n=%d events=%d  C-index=%.3f  log-rank p=%.3g\n", length(ss), sum(o$event == 1), C, lr_p))

adj <- "n/a"
if ("grade" %in% names(o)) {
  gr <- ifelse(grepl("IV", o$grade), 4L, ifelse(grepl("III", o$grade), 3L, ifelse(grepl("II\\b", o$grade), 2L, NA_integer_)))
  ok <- !is.na(gr)
  if (sum(ok) > 30) {
    fit <- coxph(Surv(o$os_months[ok], o$event[ok]) ~ scale(r[ok]) + factor(gr[ok]))
    adj <- sprintf("adj HR per SD %.2f (p=%.3g, n=%d)", exp(coef(fit)[1]), summary(fit)$coefficients[1, "Pr(>|z|)"], sum(ok))
    cat("grade-adjusted:", adj, "\n")
  }
}

out <- data.frame(item = c("n_samples", "n_events", "signature_genes_used", "C_index", "logrank_p", "grade_adjusted", "note"),
                  value = c(length(ss), sum(o$event == 1), length(g), round(C, 3), signif(lr_p, 3), adj,
                            "GSE16011 OS validation (frozen 29-gene signature); table curated from the Gravendeel 2009 Cancer Res supplement"))
write.table(out, file.path(RES, "external_gse16011_os.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

draw <- function() {
  plot(survfit(Surv(os_months, event) ~ I(r > median(r)), data = data.frame(os_months = o$os_months, event = o$event, r = r)),
       col = c("grey35", "firebrick"), lwd = 2, xlab = "months", ylab = "overall survival",
       main = sprintf("GSE16011 OS by signature (C=%.3f, log-rank p=%.2g)", C, lr_p))
  legend("topright", c("low risk", "high risk"), col = c("grey35", "firebrick"), lwd = 2, bty = "n")
}
pdf(file.path(FIG, "figS_gse16011_km.pdf"), width = 5, height = 5); draw(); dev.off()
jpeg(file.path(FIG, "figS_gse16011_km.jpg"), width = 1000, height = 1000, res = 220); draw(); dev.off()
cat("\nDONE -> results/external_gse16011_os.tsv + figS_gse16011_km.*\n")
