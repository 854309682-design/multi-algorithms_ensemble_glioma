#!/usr/bin/env Rscript
# Stage 2 (bulk prep): attach TCGA molecular annotation (IDH / MGMT / WHO grade) to the
# training clinical table so H3 (multivariate Cox) and the full nomogram become possible.
# Input:  results/tcga_training_clinical.tsv
#         $GLIOMA_DATA/TCGA/TCGA_lgggbm_clinical.tsv  (cBioPortal-style)
# Output: results/tcga_training_clinical_mol.tsv  (+ coverage printed)
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
RES <- "results"; REF <- file.path(DATA_ROOT, "TCGA")

tr  <- read.delim(file.path(RES, "tcga_training_clinical.tsv"), stringsAsFactors = FALSE)
mol <- read.delim(file.path(REF, "TCGA_lgggbm_clinical.tsv"), stringsAsFactors = FALSE)

ne <- function(x) !is.na(x) & trimws(as.character(x)) != ""

m <- merge(tr, mol, by.x = "case_submitter_id", by.y = "case", all.x = TRUE, sort = FALSE)
m <- m[match(tr$sample_submitter_id, m$sample_submitter_id), ]

m$idh_mutant <- ifelse(!ne(m$idh), NA, ifelse(m$idh == "Mutant", 1L, 0L))
m$mgmt_meth  <- ifelse(!ne(m$mgmt), NA, ifelse(m$mgmt == "Methylated", 1L, 0L))
m$who_grade  <- ifelse(!ne(m$grade), NA, as.integer(sub("G", "", m$grade)))

out <- data.frame(case_id = m$case_id, sample_submitter_id = m$sample_submitter_id,
                  case_submitter_id = m$case_submitter_id, project = m$project,
                  OS_time_days = m$OS_time_days, OS_event = m$OS_event,
                  age_at_index = m$age_at_index, gender = m$gender,
                  idh_mutant = m$idh_mutant, mgmt_meth = m$mgmt_meth, who_grade = m$who_grade)
write.table(out, file.path(RES, "tcga_training_clinical_mol.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("training n:", nrow(out), " events:", sum(out$OS_event), "\n")
cat("idh covered   :", sum(!is.na(out$idh_mutant)), "/", nrow(out), " (mutant:", sum(out$idh_mutant == 1, na.rm = TRUE), ")\n")
cat("mgmt covered  :", sum(!is.na(out$mgmt_meth)), "/", nrow(out), " (methylated:", sum(out$mgmt_meth == 1, na.rm = TRUE), ")\n")
cat("grade covered :", sum(!is.na(out$who_grade)), "/", nrow(out), " (G2/3/4:", paste(table(out$who_grade), collapse = "/"), ")\n")
comp <- !is.na(out$idh_mutant) & !is.na(out$mgmt_meth) & !is.na(out$who_grade)
cat("complete idh+mgmt+grade:", sum(comp), " events:", sum(out$OS_event[comp]), "\n")
cat("DONE\n")
