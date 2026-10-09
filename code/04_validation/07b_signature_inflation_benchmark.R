#!/usr/bin/env Rscript
# Stage B3 (methods capstone) -- CONFOUNDER-AWARE signature-inflation benchmark.
# Quantifies how much of each published gene set's apparent prognostic performance is explained
# by WHO grade / IDH, and compares it with the frozen single-cell-anchored signature.
# Per study x cohort: C-index, Spearman rho with grade and with IDH, and the grade+IDH-adjusted
# risk HR with the added-value LRT p (risk over clinical).
# Input:  results/signature_pc_coefs.tsv ;
#         <REPO_ROOT>/reference/published_signatures/*.tsv (gene<TAB>weight) ;
#         <DATA_ROOT>/CGGA/CGGA.mRNAseq_{693,325}*.zip
# Output: results/signature_inflation_benchmark.tsv
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
clean <- function(x) trimws(as.character(x))
mlab <- function(x, pos, neg) ifelse(clean(x) %in% pos, 1L, ifelse(clean(x) %in% neg, 0L, NA_integer_))

own <- read.delim(file.path(RES, "signature_pc_coefs.tsv"), stringsAsFactors = FALSE)
sigs <- list(`single-cell-anchored (this study)` = setNames(own$coefficient, own$gene))
for (f in list.files(PUB, pattern = "\\.tsv$")) {
  t <- read.delim(file.path(PUB, f), stringsAsFactors = FALSE)
  w <- if ("weight" %in% names(t)) t$weight else rep(1, nrow(t))
  sigs[[sub("\\.tsv$", "", f)]] <- setNames(w, t$gene)
}
cat("studies:", length(sigs), "\n")

coh <- list()
for (name in c("693", "325")) {
  tr <- tempfile(); tc <- tempfile()
  system(sprintf("unzip -p '%s/CGGA.mRNAseq_%s.RSEM-genes.zip' 'CGGA.mRNAseq_%s.RSEM-genes.20200506.txt' > '%s'", D, name, name, tr))
  system(sprintf("unzip -p '%s/CGGA.mRNAseq_%s_clinical.zip' 'CGGA.mRNAseq_%s_clinical.20200506.txt' > '%s'", D, name, name, tc))
  ex <- read.delim(tr, sep = "\t", check.names = FALSE); rownames(ex) <- ex[, 1]; ex <- as.matrix(ex[, -1])
  cl <- read.delim(tc, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE); names(cl)[1] <- "CGGA_ID"
  cl <- cl[clean(cl$PRS_type) == "Primary", ]
  cl$OS <- suppressWarnings(as.numeric(cl$OS)); cl$Censor <- suppressWarnings(as.numeric(cl$Censor))
  cl <- cl[!is.na(cl$OS) & cl$OS > 0 & !is.na(cl$Censor), ]
  ss <- intersect(colnames(ex), clean(cl$CGGA_ID)); ex <- log2(ex[, ss, drop = FALSE] + 1)
  cl <- cl[match(ss, clean(cl$CGGA_ID)), ]
  g <- clean(cl$Grade)
  cl$grade <- ifelse(grepl("IV", g), 4L, ifelse(grepl("III", g), 3L, ifelse(grepl("II\\b", g), 2L, NA_integer_)))
  cl$grade[!cl$grade %in% c(2, 3, 4)] <- NA
  cl$idh <- mlab(cl$IDH_mutation_status, "Mutant", "Wildtype")
  cl$age <- suppressWarnings(as.numeric(clean(cl$Age)))
  coh[[name]] <- list(X = ex, cl = cl)
}

rows <- list()
for (s in names(sigs)) {
  for (name in names(coh)) {
    X <- coh[[name]]$X; cl <- coh[[name]]$cl
    g <- intersect(names(sigs[[s]]), rownames(X)); if (length(g) < 3) next
    z <- t(scale(t(X[g, , drop = FALSE])))
    risk <- as.numeric(t(z) %*% sigs[[s]][g] / sum(abs(sigs[[s]][g])))
    C <- cind(risk, cl$OS, cl$Censor)
    rg <- suppressWarnings(cor.test(risk, cl$grade, method = "spearman")$estimate)
    ri <- suppressWarnings(cor.test(risk, cl$idh, method = "spearman")$estimate)
    ok <- !is.na(risk) & !is.na(cl$grade) & !is.na(cl$idh) & !is.na(cl$age)
    dd <- data.frame(OS = cl$OS[ok], ev = cl$Censor[ok], risk = as.numeric(scale(risk[ok])),
                     age = cl$age[ok], grade = cl$grade[ok], idh = cl$idh[ok])
    f_clin <- coxph(Surv(OS, ev) ~ age + factor(grade) + idh, data = dd)
    f_both <- coxph(Surv(OS, ev) ~ risk + age + factor(grade) + idh, data = dd)
    lrt <- tryCatch(anova(f_clin, f_both, test = "LRT")$`Pr(>|Chi|)`[2], error = function(e) NA)
    s2 <- summary(f_both)
    rows[[paste(s, name)]] <- data.frame(study = s, cohort = paste0("CGGA-", name), genes_used = length(g),
      genes_total = length(sigs[[s]]), C_index = round(C, 3),
      rho_grade = round(rg, 3), rho_IDH = round(ri, 3), n_adj = nrow(dd),
      adj_HR = round(s2$conf.int["risk", "exp(coef)"], 2),
      adj_p = signif(s2$coefficients["risk", "Pr(>|z|)"], 3), added_value_LRT_p = signif(lrt, 3))
  }
}
out <- do.call(rbind, rows)
write.table(out, file.path(RES, "signature_inflation_benchmark.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nDONE -> results/signature_inflation_benchmark.tsv\n")
cat("\nTop 15 by CGGA-693 C-index:\n")
print(head(out[order(-out$C_index), ], 15), row.names = FALSE)
