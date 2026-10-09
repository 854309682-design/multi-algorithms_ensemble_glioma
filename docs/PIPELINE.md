# Pipeline — run order and expected outputs

Every script is run from anywhere (`Rscript <script>` / `python <script>`); the bootstrap block
resolves the repository root itself. Outputs always land in `<ROOT>/results/`.
Stages must be run in the order below because each one consumes the previous stage's tables.

```bash
# R
Rscript code/03_modeling/05_signature_proteincoding.R
# Python (single-cell / spatial / figure assembly)
python code/09_spatial/01_ivygap_analysis.py
```

## 0. Environment

```bash
Rscript code/00_setup/R_packages.R                 # installs missing CRAN + Bioconductor packages
conda env create -f code/00_setup/environment.yml  # Python (scanpy / anndata / matplotlib / PyMuPDF)
```

Seeds: R `20260909`, Python `42`.

## 1. `01_sc_anchor` — single-cell discovery anchor (Python, scanpy)

| # | script | in → out |
|---|---|---|
| 01 | `01_preprocess_gse103224.py` | GSE103224 per-sample matrices → `<DATA_ROOT>/GSE103224/sc_processed.h5ad` |
| 02 | `02_annotation_gse103224.py` | `sc_processed.h5ad` → annotated object + per-cell labels |
| 03 | `03_metaprograms_gse103224.py` | malignant cells → NMF meta-programs |
| 04 | `04_tam_states_gse103224.py` | myeloid/TAM compartment → TAM states |
| 05 | `05_candidate_genes.py` | meta-programs + TAM states → `results/candidate_genes_single_cell.tsv` (537-gene candidate set) |
| 06 | `06_metaprograms_gse131928.py` | GSE131928 → cross-data meta-program consistency (Jaccard) |
| 07 | `07_cnv_malignant_qc.py` | `<DATA_ROOT>/processed/gse103224_*.h5ad` → `results/cnv_malignant_qc*.tsv` + FigS panels D/E |

## 2. `02_bulk_prep` — bulk training matrices

| script | in → out |
|---|---|
| `build_tcga_matrix.py` | TCGA STAR counts → TMM/log2-CPM training matrix + clinical table |
| `01_tcga_prep.R` | training matrix → `results/tcga_training_*` |
| `02_tcga_molecular_join.R` | cBioPortal-style clinical file → grade/IDH/MGMT annotation (`results/tcga_training_clinical_mol.tsv`) |
| `build_cgga_clinical.py` | CGGA clinical → analysis-ready clinical tables |

## 3. `03_modeling` — feature selection, model class, frozen signature

```bash
Rscript code/03_modeling/01_feature_screen.R          # univariate Cox screen (p < 0.05)
Rscript code/03_modeling/02_ensemble_cv.R             # 6-learner comparison (patient-level CV)
Rscript code/03_modeling/03_nested_cv.R               # nested CV (CoxPH vs stacking)
Rscript code/03_modeling/04_finalize_signature.R      # LASSO-Cox, 34-gene version
Rscript code/03_modeling/05_signature_proteincoding.R # PRIMARY: protein-coding refinement -> 29-gene frozen signature
Rscript code/03_modeling/06_confounding_audit.R       # grade/IDH confounding (C)
Rscript code/03_modeling/07_bulk_deg_control.R        # bulk-DEG control signature (H1)
Rscript code/03_modeling/08_sensitivity.R             # lambda.min vs lambda.1se, GBM-only training
```

Ground truth for the primary signature: `results/signature_pc_coefs.tsv` (+ `signature_pc.rds`).
Feature selection is **always inside** the cross-validation folds; the signature is never refit on a
validation cohort.

## 4. `04_validation` — external validation and controls

```bash
Rscript code/04_validation/01_external_cgga.R          # CGGA-693 / CGGA-325
Rscript code/04_validation/02_external_rembrandt.R     # GSE108474 microarray arm
Rscript code/04_validation/03_validation_meta.R        # random-effects meta-analysis (REML)
Rscript code/04_validation/04_validation_meta_pc.R     # meta-analysis of the protein-coding signature
Rscript code/04_validation/05_external_incremental.R   # incremental value over clinical / clinico-molecular
Rscript code/04_validation/06_tdAUC_calibration.R      # time-dependent AUC (1/3/5 yr) + calibration
Rscript code/04_validation/08_tcga_holdout.R           # repeated hold-out inside TCGA
Rscript code/04_validation/09_platform_robust.R        # rank-based (srank/crank) platform-robust scoring
Rscript code/04_validation/11_external_gse16011_grade_check.R  # GSE16011 direction check (BrainArray -> org.Hs.eg.db)
```

Published-signature benchmark (needs the MSigDB GMTs, see §5):

```bash
python code/04_validation/07a_gmt_to_signature_tsv.py   # MSigDB .gmt -> reference/published_signatures/*.tsv
python code/04_validation/07c_cgp_published_signatures.py  # 9 published signatures from C2:CGP
Rscript code/04_validation/07_published_signature_comparison.R
Rscript code/04_validation/07b_signature_inflation_benchmark.R
```

`12_external_gse16011_os.R` is the ready-to-run GSE16011 overall-survival arm. Overall-survival
data for GSE16011 are **not** in GEO (they live in the Gravendeel 2009 *Cancer Research* supplement,
DOI 10.1158/0008-5472.CAN-09-2307). The script **exits 0 without writing anything** until a table
`<DATA_ROOT>/GSE16011_clinical_OS.tsv` with columns `sample, OS_months, event` is supplied.

## 5. `05_clinical` … `10_synthesis`

```bash
Rscript code/05_clinical/02_nomogram_dca_mol.R         # PRIMARY clinico-molecular nomogram (risk+age+grade+IDH+MGMT) + calibration + DCA
Rscript code/06_interpretation/02_enrichment_importance_pc.R  # coefficient importance, state traceback, GO enrichment
python  code/07_state_composition/01_build_reference.py       # MuSiC reference (see README: companion object)
Rscript code/07_state_composition/02_deconvolve_music.R       # state fractions for TCGA/CGGA
Rscript code/07_state_composition/03_state_model.R            # state-composition model + confounding check
python  code/08_signalling/01_prep_cellchat.py                # CellChat input
Rscript code/08_signalling/02_cellchat.R                      # 186 L-R pairs / 122 pathways
Rscript code/08_signalling/03_signalling_score.R              # bulk signalling score + treatment interaction
python  code/09_spatial/01_ivygap_analysis.py                 # IvyGAP anatomical-structure scores
Rscript code/10_synthesis/01_model_comparison.R               # unified head-to-head model comparison
```

## 6. `99_figures` — manuscript figures

```bash
python  code/99_figures/fig1_sc_anchor.py
python  code/99_figures/fig2_states_candidates.py
Rscript code/99_figures/fig3_model_building.R
Rscript code/99_figures/fig4_external_validation.R
python  code/99_figures/fig5_clinical.py
Rscript code/99_figures/fig6_biology_closure.py
python  code/99_figures/assemble_manuscript_figures.py   # composes the 6 main + 4 supplementary composites
```

Plate scripts write `results/figures/*.pdf` + `*.jpg` (vector PDF **and** 300-dpi JPEG).
`assemble_manuscript_figures.py` combines the plates into `manuscript/Figures/Fig{1..6}` and
`FigS{1..4}` (vector-only composition via PyMuPDF) and also emits the 300-dpi JPEGs; panel letters
`(a)`, `(b)`, … are drawn inside the plate scripts, never at the composite level.

## 7. Numeric ground truth

Where this repository's scripts and any prose summary disagree, **the `results/` tables are
authoritative** (they are written directly by the analysis scripts). The tables underlying each
manuscript figure and table are listed in the Data-availability/reproducibility section of the
manuscript's supplementary information.

## 8. Known gaps / not shipped

- Raw and derived datasets are not redistributed (see `data_acquisition/README.md`).
- `reference/published_signatures/*.tsv` are regenerated by §5's `07a`/`07c` scripts from the
  MSigDB `h.all`, `c2.cp` and `c2.cgp` GMTs (`v2026.1.Hs`).
- The GSE16011 OS arm stays inert until the (closed-access) survival table is supplied.
- `code/10_synthesis/01_model_comparison.R` reads `results/ensemble_cv_cindex.tsv`, which contains a
  malformed row (a `SurvSVM` NA + a second header block); the script coerces and drops NAs — this
  reproduces the manuscript numbers.
