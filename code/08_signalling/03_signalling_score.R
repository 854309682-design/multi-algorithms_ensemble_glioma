#!/usr/bin/env Rscript
# P2.3 - Bulk signalling-module score from CellChat L-R pairs and its clinical value.
# Panel A: all genes in significant L-R pairs targeting malignant states ("incoming signal").
# Panel B: genes in pairs from Immune/Endothelial -> malignant ("microenvironment->malignant").
# Fit on TCGA, transfer FROZEN to CGGA-693/325; compare vs clinical; test treatment interaction.
# Output: results/signalling/signalling_score_performance.tsv, sg_panels.rds
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
RES <- "results"; SG <- file.path(RES, "signalling"); D <- file.path(DATA_ROOT, "CGGA")

cind <- function(risk, t, e) { ok <- !is.na(risk) & !is.na(t) & !is.na(e); risk <- risk[ok]; t <- t[ok]; e <- e[ok]
  n <- length(risk); if (n < 10) return(NA); c <- 0; d <- 0
  for (i in seq_len(n)) { if (!e[i]) next; ti <- t[i]
    for (j in seq_len(n)) { if (j == i) next
      if (t[j] > ti || (t[j] == ti && !e[j])) { d <- d + 1
        if (risk[i] > risk[j]) c <- c + 1 else if (risk[i] == risk[j]) c <- c + 0.5 } } }
  if (d == 0) NA else c / d }

lr <- read.delim(file.path(SG, "cc_lr_pairs.tsv"), stringsAsFactors = FALSE)
mal <- c("AC", "NPC", "MES", "OPC")
panel_all <- unique(c(lr$ligand[lr$target %in% mal], lr$receptor[lr$target %in% mal]))
panel_imm <- unique(c(lr$ligand[lr$target %in% mal & lr$source %in% c("Immune", "Endothelial")],
                      lr$receptor[lr$target %in% mal & lr$source %in% c("Immune", "Endothelial")]))
cat("panel_all:", length(panel_all), "genes | panel_imm:", length(panel_imm), "genes\n")
saveRDS(list(all = panel_all, immune = panel_imm), file.path(SG, "sg_panels.rds"))

score <- function(expr, panel) {                    # expr = genes x samples, log scale
  g <- intersect(panel, rownames(expr)); z <- t(scale(t(as.matrix(expr[g, , drop = FALSE]))))
  list(score = colMeans(z, na.rm = TRUE), n = length(g))
}

# TCGA (for fitting)
tcga <- read.delim(file.path(RES, "tcga_training_expr_log2cpm.tsv"), row.names = 1, check.names = FALSE)
mol <- read.delim(file.path(RES, "tcga_training_clinical_mol.tsv"), stringsAsFactors = FALSE)
mol <- mol[match(colnames(tcga), mol$sample_submitter_id), ]
sA <- score(tcga, panel_all); sB <- score(tcga, panel_imm)
dat <- data.frame(time = mol$OS_time_days, event = mol$OS_event, sA = sA$score, sB = sB$score)
cat("TCGA score A C (in-sample):", round(cind(dat$sA, dat$time, dat$event), 3),
    "| score B:", round(cind(dat$sB, dat$time, dat$event), 3), "\n")
stopifnot(!anyNA(dat))

rows <- list()
for (name in c("693", "325")) {
  tc <- tempfile(); system(sprintf("unzip -p '%s/CGGA.mRNAseq_%s_clinical.zip' 'CGGA.mRNAseq_%s_clinical.20200506.txt' > '%s'", D, name, name, tc))
  cl <- read.delim(tc, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE); names(cl)[1] <- "CGGA_ID"
  tr <- tempfile(); system(sprintf("unzip -p '%s/CGGA.mRNAseq_%s.RSEM-genes.zip' 'CGGA.mRNAseq_%s.RSEM-genes.20200506.txt' > '%s'", D, name, name, tr))
  ex <- read.delim(tr, sep = "\t", check.names = FALSE); rownames(ex) <- ex[, 1]; ex <- log2(as.matrix(ex[, -1]) + 1)
  cl <- cl[trimws(cl$PRS_type) == "Primary", ]; cl$OS <- suppressWarnings(as.numeric(cl$OS)); cl$Censor <- suppressWarnings(as.numeric(cl$Censor))
  cl <- cl[!is.na(cl$OS) & cl$OS > 0 & !is.na(cl$Censor), ]
  ss <- intersect(colnames(ex), trimws(cl$CGGA_ID)); ex <- ex[, ss]; cl <- cl[match(ss, trimws(cl$CGGA_ID)), ]
  sA2 <- score(ex, panel_all); sB2 <- score(ex, panel_imm)
  dd <- data.frame(OS = cl$OS, ev = cl$Censor, sA = sA2$score, sB = sB2$score,
                   age = suppressWarnings(as.numeric(trimws(cl$Age))),
                   grade = { g <- trimws(cl$Grade); ifelse(grepl("IV", g), 4L, ifelse(grepl("III", g), 3L, ifelse(grepl("II\\b", g), 2L, NA))) },
                   idh = ifelse(trimws(cl$IDH_mutation_status) == "Mutant", 1L, ifelse(trimws(cl$IDH_mutation_status) == "Wildtype", 0L, NA)),
                   codel = ifelse(trimws(cl[["1p19q_codeletion_status"]]) == "Codel", 1L, ifelse(trimws(cl[["1p19q_codeletion_status"]]) == "Non-codel", 0L, NA)),
                   mgmt = ifelse(trimws(cl[["MGMTp_methylation_status"]]) %in% c("methylated","Methylated"), 1L, ifelse(trimws(cl[["MGMTp_methylation_status"]]) %in% c("un-methylated","Un-methylated"), 0L, NA)),
                   chemo = suppressWarnings(as.numeric(trimws(cl$Chemo_status))), radio = suppressWarnings(as.numeric(trimws(cl$Radio_status))))
  cc <- complete.cases(dd[, c("OS","ev","sA","sB","age","grade","idh","codel","mgmt")]); dd <- dd[cc, ]
  fclin <- coxph(Surv(OS, ev) ~ age + factor(grade) + idh + codel + mgmt, data = dd)
  fb <- coxph(Surv(OS, ev) ~ sB + age + factor(grade) + idh + codel + mgmt, data = dd)
  lrtB <- anova(fclin, fb, test = "LRT")$`Pr(>|Chi|)`[2]
  int <- tryCatch({ dd2 <- dd[!is.na(dd$chemo), ]; summary(coxph(Surv(OS, ev) ~ sB * chemo, data = dd2))$coefficients["sB:chemo", "Pr(>|z|)"] }, error = function(e) NA)
  hrx <- function(sub) tryCatch({ h <- summary(coxph(Surv(OS, ev) ~ sB, data = sub))$conf.int["sB", ]; c(h["exp(coef)"], h["lower .95"], h["upper .95"]) }, error = function(e) c(NA, NA, NA))
  h1 <- hrx(dd[dd$chemo == 1 & !is.na(dd$chemo), ]); h0 <- hrx(dd[dd$chemo == 0 & !is.na(dd$chemo), ])
  cat(sprintf("\n== CGGA-%s == n=%d ev=%d\n  scoreA C=%.3f | scoreB C=%.3f | clinical C=%.3f | clinical+scoreB C=%.3f | LRT p=%.2g | scoreB x chemo p=%.3g\n",
              name, nrow(dd), sum(dd$ev), cind(dd$sA, dd$OS, dd$ev), cind(dd$sB, dd$OS, dd$ev),
              cind(predict(fclin), dd$OS, dd$ev), cind(predict(fb), dd$OS, dd$ev), lrtB, int))
  rows[[name]] <- data.frame(cohort = paste0("CGGA-", name), n = nrow(dd), events = sum(dd$ev),
                             C_score_all = round(cind(dd$sA, dd$OS, dd$ev), 3), C_score_immune = round(cind(dd$sB, dd$OS, dd$ev), 3),
                             C_clinical = round(cind(predict(fclin), dd$OS, dd$ev), 3),
                             C_clinical_score = round(cind(predict(fb), dd$OS, dd$ev), 3),
                             added_value_LRT_p = signif(lrtB, 3), score_x_chemo_p = signif(int, 3),
                             HR_chemo = round(h1[1], 2), lo_chemo = round(h1[2], 2), hi_chemo = round(h1[3], 2),
                             HR_nochemo = round(h0[1], 2), lo_nochemo = round(h0[2], 2), hi_nochemo = round(h0[3], 2))
}
out <- do.call(rbind, rows)
write.table(out, file.path(SG, "signalling_score_performance.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nDONE\n"); print(out, row.names = FALSE)
