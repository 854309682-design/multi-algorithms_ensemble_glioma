# Data acquisition — where the inputs come from and where they must sit

No dataset is redistributed in this repository. Everything below is public; please cite the
original studies. All accessions/roles are also machine-readable in `dataset_index.tsv`.

## Expected on-disk layout

`DATA_ROOT` defaults to `<DATA_ROOT>/` **next to the cloned repository** and can be overridden with
`export GLIOMA_DATA=/path/to/Datasets` (see the repository README). The pipeline expects:

```
<DATA_ROOT>/
├── GSE103224/                      # single-cell discovery anchor (8 per-sample matrices)
│   └── GSM*.filtered.matrix.txt.gz
├── GSE131928/                      # single-cell reproducibility (Neftel et al.)
├── GSE16011_series_matrix.txt.gz   # GSE16011 expression (GPL8542)
├── GPL8542_family.soft.gz          # GSE16011 platform annotation
├── GSE108474_REMBRANDT_GeneExpression.txt.gz
├── GSE108474_REMBRANDT_clinical.data.txt.gz
├── GSE108474_REMBRANDT_biospecimen_mapping_GEO.txt.gz
├── CGGA/                           # CGGA.mRNAseq_693 / _325 RSEM + clinical
├── IvyGAP/                         # IvyGAP expression + anatomical structure annotation
├── TCGA/                           # GDC STAR counts + clinical (built by tcga_*.sh/py)
├── processed/                      # pipeline intermediates written by stage 1 (h5ad)
│   └── GSE16011_clinical_OS.tsv    # optional: closed-access OS table for script 12 (see below)
└── Paper_1/…                       # companion-project objects for the MuSiC reference (see README)
```

## Sources

| dataset | role | source |
|---|---|---|
| GSE103224 | single-cell discovery anchor (8 high-grade glioma patients) | GEO accession GSE103224 |
| GSE131928 | single-cell reproducibility (Neftel et al.) | GEO accession GSE131928 |
| TCGA-GBM/LGG | bulk training (STAR counts) | GDC portal — `tcga_download.sh` + `tcga_manifest.py` |
| CGGA mRNAseq_693 / _325 | external validation (RNA-seq) | http://www.cgga.org.cn/download — `download_cgga.sh` |
| GSE108474 (REMBRANDT) | cross-platform validation (microarray) | GEO accession GSE108474 |
| GSE16011 (GPL8542) | cross-platform direction check | GEO accession GSE16011 |
| IvyGAP | spatial/anatomical-structure check | Ivy Glioblastoma Atlas Project |

## Scripts

```bash
bash download_geo.sh        # GEO series matrices used for the microarray arms
bash download_cgga.sh       # CGGA mRNAseq_693 / _325 (unzips the archives; results/ reads the .txt inside)
bash tcga_download.sh       # GDC STAR-count files for TCGA-GBM/LGG (needs tcga_manifest.py first)
python tcga_manifest.py     # builds the GDC file manifest
python tcga_clinical.py     # pulls the TCGA clinical/case JSON
md5sum -c checksums/geo_cgga.md5   # verify what you downloaded
```

All scripts honour `GLIOMA_DATA` and default to `<parent of this repository>/Datasets`.

## Data formats and traps (verified during the analysis)

- **GSE103224** — tab-delimited `<ENSG.ver>\t<symbol>\t<count>…`, **no header row, no barcode**;
  barcodes are constructed as `GSM<uid>_PJ###_CellN`. Only the `.txt.gz` files are matrices.
- **CGGA RSEM** — first column `Gene_Name`, then one column per `CGGA_XXXX` sample. Clinical holds
  `IDH_mutation_status, 1p19q_codeletion_status, MGMTp_methylation_status, Grade, Age, OS, Censor,
  PRS_type, Histology, Gender` and mixes primary **and recurrent** tumours (`PRS_type`).
  R cannot read the CGGA zips directly (`__MACOSX` entries) — `unzip -p` first.
- **GSE108474** — probe-level log2 matrix, rows are Affymetrix probes, clinical keyed by `SUBJECT_ID`
  (`EVENT_OS`, `OVERALL_SURVIVAL_MONTHS`, `WHO_GRADE`, `DISEASE_TYPE`).
- **GSE16011** — GEO series matrix on GPL8542. Map probes through the **BrainArray `<entrez>_at`**
  identifiers via `org.Hs.eg.db`, **not** `hgu133plus2.db`.
- **GSE16011 overall survival is not in GEO.** Script
  `code/04_validation/12_external_gse16011_os.R` needs `<DATA_ROOT>/GSE16011_clinical_OS.tsv`
  (columns `sample`, `OS_months`, `event`) derived from the Gravendeel 2009 *Cancer Research*
  supplement (DOI 10.1158/0008-5472.CAN-09-2307, closed access). Until it exists the script exits 0
  without writing anything.
