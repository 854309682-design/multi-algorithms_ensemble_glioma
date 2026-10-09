#!/usr/bin/env python3
"""P1.2 - Build the single-cell reference for bulk deconvolution (MuSiC).

States:
  malignant (GSE103224 states.h5ad): AC / NPC / MES / OPC
  non-malignant (GSE103224, score argmax): Immune / Endothelial / Oligodendrocyte
  myeloid (Paper_1 step05_myeloid_reclustered.h5ad): TAM-IS / myeloid_other / microglia

Outputs (results/state_composition/):
  ref_counts.mtx      genes x cells (raw counts, sparse, Matrix Market)
  ref_genes.txt       gene symbols (rows)
  ref_cells.tsv       cell, cellType, sampleID
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
import scipy.sparse as sp
from scipy.io import mmwrite

ROOT = str(REPO_ROOT)
OUT = os.path.join(REPO_ROOT, "results", "state_composition")
os.makedirs(OUT, exist_ok=True)
MAX_PER_TYPE = 500
rng = np.random.default_rng(42)

# bulk gene space (for intersection)
bulk_genes = pd.read_csv(os.path.join(REPO_ROOT, "results", "tcga_training_expr_log2cpm.tsv"),
                         sep="\t", index_col=0, usecols=[0]).index.astype(str).tolist()
bulk_genes = set(bulk_genes)

# ---- malignant + non-malignant from GSE103224 (full-gene qc.h5ad; labels from states.h5ad) ----
q = ad.read_h5ad(os.path.join(DATA_ROOT, "processed", "gse103224_qc.h5ad"))
st = ad.read_h5ad(os.path.join(DATA_ROOT, "processed", "gse103224_states.h5ad"))
st = st[q.obs_names]                                   # align to qc cell order (barcodes identical)
cnt = q.X
genes = np.array(q.var_names.astype(str))
obs = st.obs.copy()
mal = obs["malignant"].values.astype(int)
state = obs["state"].astype(str).values
sample = obs["sample"].astype(str).values
scores = obs[["score_immune", "score_endothelial", "score_oligo"]].to_numpy()
nonmal_lab = np.array(["Immune", "Endothelial", "Oligodendrocyte"])[np.argmax(scores, axis=1)]
ctype = np.where(mal == 1, state, nonmal_lab)
blocks = [(genes, cnt, ctype, sample)]

# ---- myeloid atlas from Paper_1 ----
m = ad.read_h5ad(os.path.join(SIBLING_ROOT, "Paper_1", "objects", "step05_myeloid_reclustered.h5ad"))
mcnt = m.layers["counts"]
mgenes = np.array(m.var_names.astype(str))
mctype = m.obs["state"].astype(str).values
msample = m.obs["dataset"].astype(str).values if "dataset" in m.obs else np.array(["myeloid"] * m.shape[0])
blocks.append((mgenes, mcnt, mctype, msample))

# ---- intersect to a common gene space present in bulk ----
common = set(blocks[0][0])
for g, _, _, _ in blocks[1:]:
    common &= set(g)
common &= bulk_genes
common = sorted(common)
print("common genes:", len(common))

mats, ctypes, samples = [], [], []
for g, c, ct, sp_ in blocks:
    idx = np.array([np.where(g == x)[0][0] for x in common])
    mats.append(sp.csr_matrix(c)[:, idx])          # cells x common genes
    ctypes.append(ct)
    samples.append(sp_)
X = sp.vstack(mats).tocsr()                        # cells x genes
ctypes = np.concatenate(ctypes)
samples = np.concatenate(samples)
print("reference cells before subsample:", X.shape[0], "states:", dict(pd.Series(ctypes).value_counts()))

# ---- subsample up to MAX_PER_TYPE per state ----
keep = []
for s in np.unique(ctypes):
    idx = np.where(ctypes == s)[0]
    if len(idx) > MAX_PER_TYPE:
        idx = rng.choice(idx, MAX_PER_TYPE, replace=False)
    keep.append(idx)
keep = np.sort(np.concatenate(keep))
X = X[keep]
ctypes = ctypes[keep]
samples = samples[keep]
print("reference cells after subsample:", X.shape[0])

# ---- write genes x cells ----
Xgc = X.T.tocoo()                                   # genes x cells
mmwrite(os.path.join(OUT, "ref_counts.mtx"), Xgc)
pd.Series(common).to_csv(os.path.join(OUT, "ref_genes.txt"), index=False, header=False)
pd.DataFrame({"cell": [f"c{i}" for i in range(X.shape[0])], "cellType": ctypes,
              "sampleID": samples}).to_csv(os.path.join(OUT, "ref_cells.tsv"), sep="\t", index=False)
print("wrote ref_counts.mtx (genes x cells =", Xgc.shape, "), ref_genes.txt, ref_cells.tsv")
