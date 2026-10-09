#!/usr/bin/env python3
"""
Build a gene-by-sample TCGA-GBM/LGG STAR-counts matrix (primary-tumor, one aliquot per case).

Inputs (<DATA_ROOT>/TCGA/):
  - tcga_gene_expr_manifest.json  (file_id/case_id/sample_submitter_id/sample_type/...)
  - GeneExpr/<uuid>.rna_seq.augmented_star_gene_counts.tsv  (per sample)
  - tcga_clinical_survival.tsv     (case-level OS; joined later)

Outputs:
  - <DATA_ROOT>/TCGA/TCGA_GBM_LGG_star_counts_matrix.tsv.gz   (rows=gene_id+gene_name, cols=sample)
  - <DATA_ROOT>/TCGA/TCGA_GBM_LGG_sample_meta.tsv              (sample -> case, primary status)
"""
# --- portable project roots: do not edit (see config/paths.py) ----------------
import os
import sys
import pathlib


def _find_root(_d):
    while True:
        if (_d / "config" / "paths.py").exists():
            return _d
        if _d == _d.parent:
            raise SystemExit("repository root not found -- set GLIOMA_ROOT")
        _d = _d.parent


_here = pathlib.Path(__file__).resolve().parent
REPO_ROOT = pathlib.Path(os.environ["GLIOMA_ROOT"]).resolve() if os.environ.get("GLIOMA_ROOT") else _find_root(_here)
sys.path.insert(0, str(REPO_ROOT / "config"))
from paths import ROOT, DATA_ROOT, REF_ROOT, SIBLING_ROOT  # noqa: E402
# -----------------------------------------------------------------------------

import json, os, gzip, csv, sys
import pandas as pd
from collections import defaultdict

TCGA = str(DATA_ROOT / "TCGA")
GEXPR = os.path.join(TCGA, "GeneExpr")

mani = json.load(open(os.path.join(TCGA, "tcga_gene_expr_manifest.json")))

# only primary-tumor samples; group by case and keep the aliquot with max total counts
prim = [e for e in mani if e.get("sample_type") == "Primary Tumor"]
print("primary-tumor file entries:", len(prim))

# read counts helper
def read_counts(path):
    df = pd.read_csv(path, sep="\t", comment="#")
    df = df[["gene_name", "unstranded"]].rename(columns={"unstranded": "count"})
    df = df[df["gene_name"].notna() & (df["gene_name"] != "")]
    df = df.groupby("gene_name", as_index=False)["count"].sum()  # collapse duplicate symbols
    return df.set_index("gene_name")["count"]

# choose one file per case (max total reads) among primary samples
by_case = defaultdict(list)
for e in prim:
    by_case[e["case_id"]].append(e)
chosen = []
for cid, es in by_case.items():
    # total counts = file_size proxy not read; compute real total lazily -> instead pick first, but better max file_size as pseudo-reads
    es = sorted(es, key=lambda e: e.get("file_size", 0), reverse=True)
    chosen.append(es[0])
print("one-aliquot-per-case chosen (primary):", len(chosen))

# merge (we read gene columns from each file and accumulate into a dict of dicts)
# Build via pandas concat of series (memory-efficient enough for ~800 x 23k)
series = {}
order = []
for i, e in enumerate(chosen):
    sid = e["sample_submitter_id"] or os.path.basename(e["file_name"])
    path = os.path.join(GEXPR, e["file_name"])
    if not os.path.exists(path):
        print("MISSING", e["file_name"]); continue
    s = read_counts(path)
    series[sid] = s
    order.append(sid)
    if (i+1) % 150 == 0:
        print("  read", i+1, "files")

mat = pd.concat(series.values(), axis=1, keys=order)
mat = mat.fillna(0.0).astype("int32")
mat.index.name = "gene_name"
mat = mat.sort_index()

# write output
outmat = os.path.join(TCGA, "TCGA_GBM_LGG_star_counts_matrix.tsv.gz")
mat.to_csv(outmat, sep="\t", compression="gzip")
print("matrix vars x samples:", mat.shape)

# sample metadata
meta = []
for e in chosen:
    meta.append({
        "sample_submitter_id": e["sample_submitter_id"],
        "case_id": e["case_id"],
        "case_submitter_id": e.get("case_submitter_id"),
        "sample_type": e["sample_type"],
        "project": e.get("project"),
        "file_id": e["file_id"],
    })
mdf = pd.DataFrame(meta)
mdf.to_csv(os.path.join(TCGA, "TCGA_GBM_LGG_sample_meta.tsv"), sep="\t", index=False)
print("saved meta:", mdf.shape)
print("DONE")
