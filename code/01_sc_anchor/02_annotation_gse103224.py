#!/usr/bin/env python3
"""
Stage 1 (part 2): marker-based cell-type annotation for GSE103224.
Scores canonical marker sets (scanpy tl.score_genes), assigns a coarse cell type,
and saves the annotated object. Also exports a raw-counts gene-symbol x cell matrix
for the CopyKAT CNV malignant call.
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

import os
import numpy as np
import scanpy as sc
import anndata as ad

BASE = str(DATA_ROOT / "GSE103224")
SCRIPT = str(REPO_ROOT / "code" / "01_sc_anchor")

sc.settings.verbosity = 2

adata = sc.read_h5ad(os.path.join(BASE, "sc_processed.h5ad"))
print("loaded:", adata.shape)

# gene symbol lookup (var_names are ENSG.ver_symbol); build reverse map symbol -> var_name
sym = adata.var["gene_symbol"].astype(str).to_dict()
var_by_symbol = {}
for v, s in zip(adata.var_names, adata.var["gene_symbol"].astype(str)):
    var_by_symbol.setdefault(s, v)

def score(markers, key):
    var_ids = []
    for g in markers:
        if g in var_by_symbol:
            var_ids.append(var_by_symbol[g])
    if not var_ids:
        print("  no var for", key); return 0
    try:
        sc.tl.score_genes(adata, var_ids, score_name=key, use_raw=False)
    except Exception as e:
        print("score fail", key, e); return 0
    return len(var_ids)

markers = {
    "T":        ["CD3D","CD3E","CD3G","IL7R","LCK","TRAC","CD2"],
    "NK":       ["NKG7","GNLY","KLRD1","GZMB","PRF1"],
    "B_plasma": ["CD79A","MS4A1","CD19","MZB1","JCHAIN","IGHG1"],
    "Myeloid":  ["CD68","CD163","C1QA","C1QB","C1QC","AIF1","ITGAM","FCGR3A","TMEM119","SPP1","MRC1","LYZ","S100A8"],
    "Oligo":    ["OLIG1","OLIG2","MBP","PLP1","SOX10","MOG","PDGFRA","MAG"],
    "Endo":     ["PECAM1","VWF","CLDN5","CDH5","FLT1","RAMP2"],
    "Tumor":    ["GFAP","SOX2","EGFR","NES","VIM","OLIG2","BCAN","CD44"],
}
for k, m in markers.items():
    n = score(m, f"score_{k}")
    print(f"  scored {k}: {n}/{len(m)} markers present")

# coarse assignment: assign the lineage with the MAX score across all 7 (incl Tumor).
# This prevents myeloids/immune/oligo/endo from leaking into the Tumor group.
lineages = ["T","NK","B_plasma","Myeloid","Oligo","Endo","Tumor"]
data = np.vstack([adata.obs[f"score_{k}"].values for k in lineages]).T  # cells x 7
idx = np.argmax(data, axis=1)
celltype = np.array(lineages)[idx]
adata.obs["celltype_marker"] = celltype
print("cell types:", adata.obs["celltype_marker"].value_counts().to_dict())

adata.write(os.path.join(BASE, "sc_annotated.h5ad"))
print("saved -> sc_annotated.h5ad")
print("cell type summary:")
print(adata.obs["celltype_marker"].value_counts())

