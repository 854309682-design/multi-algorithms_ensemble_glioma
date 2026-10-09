# multi-algorithms_ensemble_glioma

Analysis code for the manuscript **“A single-cell-anchored machine-learning prognostic signature
for glioma: multi-algorithm ensemble with multi-cohort external validation.”**

The repository contains the complete analysis pipeline: single-cell anchor discovery and
meta-program derivation, candidate-gene assembly, penalised-Cox modelling with a nested
cross-validated model-class comparison, multi-cohort external validation, confounding/incremental
value analyses, a clinico-molecular nomogram, state-composition and cell–cell-signalling
follow-ups, spatial check, and the figure-assembly scripts.

**Headline result.** A frozen **29-gene protein-coding Cox proportional-hazards signature** derived
from single-cell malignant-cell-state and TAM-state programs, externally validated in independent
RNA-seq cohorts at a **pooled C-index of 0.760 (95% CI 0.738–0.781, I² = 0)** (CGGA-693 0.761,
CGGA-325 0.759). Cross-platform (microarray) transfer is limited and is reported as a documented
limitation (REMBRANDT C = 0.483).

This is a **pure computational study**. All findings are exploratory/associational; no wet-lab
experiments are included.

---

## Repository layout

```
config/            path resolution shared by all scripts (paths.R, paths.py)
code/              00_setup 01_sc_anchor 02_bulk_prep 03_modeling 04_validation 05_clinical
                   06_interpretation 07_state_composition 08_signalling 09_spatial 10_synthesis 99_figures
reference/         small reference inputs shipped with the repo
                     protein_coding_genes.txt     (required by 03_modeling/05, 03_modeling/08, 04_validation/08)
                     published_signatures/        (generated: see docs/PIPELINE.md §5)
data_acquisition/  where to obtain each dataset, per-dataset download scripts + md5 checksums
docs/PIPELINE.md   run order, per-script description, seeds, expected outputs
```

## Paths — how this repository finds your data

No absolute paths are hard-coded. Every script starts with a short bootstrap block that

1. locates the repository root (walking up from the script until `config/paths.R` is found), and
2. sources `config/paths.{R,py}`, which exposes `ROOT` / `REPO_ROOT`, `DATA_ROOT`, `REF_ROOT`
   and `SIBLING_ROOT`.

Defaults and overrides:

| variable | meaning | default | override |
|---|---|---|---|
| repository root | the cloned repository (`ROOT`) | auto-detected | `GLIOMA_ROOT` |
| `DATA_ROOT` | tree holding the raw + intermediate datasets | `<parent of the clone>/Datasets` | `GLIOMA_DATA` |
| `REF_ROOT` | shipped reference tables (`reference/`) | `<ROOT>/reference` | — |
| `SIBLING_ROOT` | companion project objects (see below) | `<parent of the clone>` | `GLIOMA_SIBLING` |

So if you clone into `~/projects/multi-algorithms_ensemble_glioma` and keep the data in
`~/projects/Datasets`, **nothing has to be configured**. Otherwise:

```bash
export GLIOMA_DATA=/path/to/Datasets        # raw data tree
export GLIOMA_ROOT=/path/to/clone           # only needed when running from outside the clone
```

Scripts write their outputs to `<ROOT>/results/` (created on demand) — this directory is not
tracked by git. Each stage consumes the outputs of the previous ones, so run the pipeline in the
order documented in `docs/PIPELINE.md`.

### Companion project objects (one documented external dependency)

`code/07_state_composition/01_build_reference.py` builds the MuSiC reference from the GSE103224
states **and** the myeloid/TAM atlas of the companion study, expected at
`<SIBLING_ROOT>/Paper_1/objects/step05_myeloid_reclustered.h5ad`. Set `GLIOMA_SIBLING` to the
directory containing that project. The deconvolution results used in the manuscript are included
in the supplement; the script is provided so the reference construction is auditable.

## Environment

- **R 4.6.1** — `code/00_setup/R_packages.R` lists every CRAN/Bioconductor package used and installs
  the missing ones from `https://cloud.r-project.org`. The `sessionInfo()` of the analysis machine
  is in `code/00_setup/sessionInfo_R.txt`.
- **Python 3.11+** — `code/00_setup/environment.yml` (conda) lists the packages required by the
  single-cell, spatial and figure-assembly scripts (scanpy, anndata, harmonypy, numpy, pandas,
  scipy, scikit-learn, matplotlib, PyMuPDF).
- **Seeds** — R `20260909`, Python `42` (set in every script that involves clustering, resampling or
  model fitting).

## Data availability

All datasets are public and were re-used from the original publications; none are redistributed
here. Accessions, per-file download scripts and md5 checksums are in `data_acquisition/`
(see its README for the expected on-disk layout):

| role | dataset(s) |
|---|---|
| single-cell discovery anchor | GSE103224 |
| single-cell reproducibility | GSE131928 (Neftel et al.) |
| bulk training | TCGA-GBM/LGG (GDC STAR counts) |
| external validation (RNA-seq) | CGGA mRNAseq_693, CGGA mRNAseq_325 |
| external validation (microarray) | GSE108474 (REMBRANDT), GSE16011 / GPL8542 |
| spatial check | IvyGAP |

## Honest caveats (also in the manuscript)

- The signature is an **independent** prognostic factor but only **modestly incremental** over
  standard clinico-molecular variables (clinical-only C ≈ 0.79–0.80 vs signature-only 0.76).
- The candidate genes are largely **grade/IDH proxies**: 474/537 candidates are significant
  unadjusted vs 95/537 after grade+IDH adjustment.
- The single-cell-anchored candidate set beats an equivalent bulk-DEG control only modestly.
- **Cross-platform microarray transfer fails directionally** (REMBRANDT C = 0.483); rank-based
  scoring does not rescue it.
- The Stage-1 non-malignant single-cell classes in the discovery object are marker-score
  reconstructions (the original per-cell annotation object was not retained).

## Citation

Please cite the manuscript above (details to be finalised at publication). Code released under the
MIT licence — see `LICENSE`.
