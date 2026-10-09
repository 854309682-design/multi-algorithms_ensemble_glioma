## ---------------------------------------------------------------------------
## Resolved paths shared by every R script in this repository.
##
## Sourced automatically by the bootstrap block at the top of each script; it is
## not meant to be edited for a normal run. Override with environment variables:
##
##   GLIOMA_ROOT     repository root                        (default: auto-detected)
##   GLIOMA_DATA     raw-data tree                          (default: <sibling>/Datasets)
##   GLIOMA_SIBLING  directory holding companion projects   (default: dirname(ROOT))
##
## `REPO_ROOT` and `ROOT` are the same thing (the repository root); `ROOT` is kept
## because the analysis scripts use it as the working-directory anchor.
## ---------------------------------------------------------------------------

ROOT         <- REPO_ROOT
DATA_ROOT    <- normalizePath(Sys.getenv("GLIOMA_DATA", unset = file.path(dirname(REPO_ROOT), "Datasets")), mustWork = FALSE)
SIBLING_ROOT <- normalizePath(Sys.getenv("GLIOMA_SIBLING", unset = dirname(REPO_ROOT)), mustWork = FALSE)
REF_ROOT     <- file.path(ROOT, "reference")   # small reference tables shipped with the repo (+ published_signatures/)

if (!dir.exists(file.path(ROOT, "code"))) {
  stop("repository root does not contain code/ -- set GLIOMA_ROOT to the cloned repository")
}
