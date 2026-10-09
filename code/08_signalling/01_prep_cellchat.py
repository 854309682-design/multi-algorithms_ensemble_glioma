#!/usr/bin/env python3
"""P2.1 - Build the CellChat input from GSE103224: normalized expression + cell-state labels.

Labels: malignant AC/NPC/MES/OPC + Immune/Endothelial/Oligodendrocyte (score argmax).
Outputs (results/signalling/):
  cc_expr.mtx      genes x cells (log-normalized, sparse)
  cc_genes.txt     gene symbols
  cc_labels.tsv    cell, label
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
import pandas as pd
import anndata as ad
import scanpy as sc
import scipy.sparse as sp
from scipy.io import mmwrite

ROOT = str(REPO_ROOT)
OUT = os.path.join(REPO_ROOT, "results", "signalling")
os.makedirs(OUT, exist_ok=True)
MAX_PER = 500
rng = np.random.default_rng(42)

q = ad.read_h5ad(os.path.join(DATA_ROOT, "processed", "gse103224_qc.h5ad"))
st = ad.read_h5ad(os.path.join(DATA_ROOT, "processed", "gse103224_states.h5ad"))
st = st[q.obs_names]
obs = st.obs
mal = obs["malignant"].values.astype(int)
state = obs["state"].astype(str).values
scores = obs[["score_immune", "score_endothelial", "score_oligo"]].to_numpy()
nonmal = np.array(["Immune", "Endothelial", "Oligodendrocyte"])[np.argmax(scores, axis=1)]
lab = np.where(mal == 1, state, nonmal)
labels = ["AC", "NPC", "MES", "OPC", "Immune", "Endothelial", "Oligodendrocyte"]
keep = [i for i in range(q.n_obs) if lab[i] in labels]
q = q[keep]; lab = lab[keep]
sc.pp.normalize_total(q, target_sum=1e4); sc.pp.log1p(q)
print("cells:", q.n_obs, "labels:", dict(pd.Series(lab).value_counts()))

idx = []
for L in labels:
    w = np.where(lab == L)[0]
    idx.append(rng.choice(w, min(MAX_PER, len(w)), replace=False))
idx = np.sort(np.concatenate(idx))
X = sp.csr_matrix(q.X)[idx]                 # cells x genes
genes = np.array(q.var_names.astype(str))
mmwrite(os.path.join(OUT, "cc_expr.mtx"), sp.coo_matrix(X.T))   # genes x cells
pd.Series(genes).to_csv(os.path.join(OUT, "cc_genes.txt"), index=False, header=False)
pd.DataFrame({"cell": [f"c{i}" for i in range(len(idx))], "label": lab[idx]}).to_csv(
    os.path.join(OUT, "cc_labels.tsv"), sep="\t", index=False)
print("wrote cc_expr.mtx", X.T.shape, "; labels:", dict(pd.Series(lab[idx]).value_counts()))
