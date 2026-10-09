#!/usr/bin/env python3
"""
Stage 1 (part 1): GSE103224 scRNA load + QC + integration.
Loads the 8 per-sample filtered matrices, builds barcodes `PJ###_Cell####`,
filters cells/genes, normalizes to 1e4 + log1p, HVG(3000), PCA(50),
Harmony batch-correct by patient, Leiden clustering.
Output: <DATA_ROOT>/GSE103224/sc_processed.h5ad  (+ per-sample summary)
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

import os, glob, re
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import harmonypy as hm  # ensure available

BASE = str(DATA_ROOT / "GSE103224")
OUT = os.path.join(BASE, "sc_processed.h5ad")

sc.settings.verbosity = 3
sc.settings.seed = 20260909
np.random.seed(20260909)

files = sorted(glob.glob(os.path.join(BASE, "*.filtered.matrix.txt.gz")))
print("samples:", len(files))

adatas = []
for f in files:
    samp = os.path.basename(f).split(".")[0].replace("filtered", "").strip("_")  # e.g. GSM2758471_PJ016
    # parse patient id
    m = re.search(r"PJ(\d+)", os.path.basename(f))
    patient = f"PJ{m.group(1)}" if m else samp
    df = pd.read_csv(f, sep="\t", header=None)
    gene = df.iloc[:, 0].astype(str) + "_" + df.iloc[:, 1].astype(str)  # ENSG.ver_symbol
    symbol = df.iloc[:, 1].astype(str)
    counts = df.iloc[:, 2:].T  # cells x genes
    counts.index = [f"{patient}_Cell{n:04d}" for n in range(1, counts.shape[0] + 1)]
    counts.columns = gene.values
    # store symbol in var
    a = sc.AnnData(counts.astype("float32"))
    a.var["gene_symbol"] = symbol.values
    a.obs["patient"] = patient
    a.obs["sample"] = samp
    a.obs["n_genes_raw"] = (counts > 0).sum(axis=1).values
    adatas.append(a)
    print(f"  {patient}: cells={a.n_obs} genes={a.n_vars}")

adata = ad.concat(adatas, join="outer", label="sample_concat")
adata.obs_names_make_unique()
print("combined:", adata.shape)

# recover gene_symbol from var_names (encoded as ENSG.ver_symbol); fill NaN counts from outer join
def parse_sym(n):
    i = n.find("_")
    return n[i+1:] if i > 0 else n
adata.var["gene_symbol"] = [parse_sym(n) for n in adata.var_names]
adata.X = np.nan_to_num(adata.X, nan=0.0)

# QC
adata.var["mt"] = adata.var["gene_symbol"].astype(str).str.startswith("MT-").values
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
print("QC metrics computed")
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_cells(adata, max_genes=8000)
sc.pp.filter_genes(adata, min_cells=3)
print("after cell/gene filter:", adata.shape)

# normalize + log
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# HVG (use classic 'seurat' flavor — no scikit-misc needed; seurat_v3 would require scikit-misc)
sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat")
print("HVG:", adata.var["highly_variable"].sum())

# PCA on HVG
adata_pca = adata.copy()
sc.pp.pca(adata_pca, n_comps=50, use_highly_variable=True)

adatas = None  # free

# Harmony by patient
harmony_features = adata_pca.obsm["X_pca"]
obs = adata_pca.obs
meta = obs[["patient"]].reset_index().rename(columns={"index": "barcode"})
hm_out = hm.run_harmony(harmony_features, obs[["patient"]], ["patient"], verbose=False)
adata_pca.obsm["X_harmony"] = np.asarray(hm_out.Z_corr)

# neighborhood + Leiden
sc.pp.neighbors(adata_pca, use_rep="X_harmony", n_neighbors=15, n_pcs=50)
for res in [0.4, 1.0, 1.4, 2.0]:
    sc.tl.leiden(adata_pca, resolution=res, key_added=f"leiden_{res}", flavor="igraph", n_iterations=2, directed=False)
sc.tl.umap(adata_pca, min_dist=0.3)

# write
adata_pca.write(OUT)
print("saved ->", OUT)
print("final annotated object:", adata_pca.shape)
