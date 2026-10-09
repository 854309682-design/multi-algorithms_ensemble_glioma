#!/usr/bin/env python3
"""
Stage 1 (part 4): TAM / myeloid-state markers for GSE103224.
Subclusters the Myeloid compartment, rank_genes_groups to get per-state markers.
Output: <DATA_ROOT>/GSE103224/tam_state_markers.tsv
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
import scanpy as sc

BASE = str(DATA_ROOT / "GSE103224")
OUT = os.path.join(BASE, "tam_state_markers.tsv")

sc.settings.verbosity = 0
np.random.seed(20260909)

adata = sc.read_h5ad(os.path.join(BASE, "sc_annotated.h5ad"))
ad = adata[adata.obs["celltype_marker"] == "Myeloid"].copy()
print("myeloid cells:", ad.shape)

# re-derive embedding + neighbors on the myeloid subset
sc.pp.pca(ad, n_comps=30)
sc.pp.neighbors(ad, n_neighbors=15, n_pcs=30)
for res in [0.4, 0.6, 1.0]:
    sc.tl.leiden(ad, resolution=res, key_added=f"leiden_{res}", flavor="igraph", n_iterations=2, directed=False)
sc.tl.umap(ad, min_dist=0.3)

# use a moderated resolution for states (0.6)
res = "leiden_0.6"
print("state clusters @0.6:", ad.obs[res].value_counts().to_dict())
sc.tl.rank_genes_groups(ad, groupby=res, method="wilcoxon", use_raw=False)

# extract top-50 markers per state
rows = []
for st in np.sort(ad.obs[res].unique()):
    try:
        gs = sc.get.rank_genes_groups_df(ad, group=str(st)) if hasattr(sc.get, "rank_genes_groups_df") else None
    except Exception:
        gs = None
    if gs is None:
        # fallback: build manually
        gs = pd.DataFrame({"names": np.asarray(ad.uns["rank_genes_groups"]["names"])[int(st)],
                           "logfoldchanges": np.asarray(ad.uns["rank_genes_groups"]["logfoldchanges"])[int(st)],
                           "pvals_adj": np.asarray(ad.uns["rank_genes_groups"]["pvals_adj"])[int(st)]})
    top = gs.head(50)
    genes = [top["names"].values[i] for i in range(min(50, len(top)))]
    rows.append({"tam_state": f"TAM_{st}", "n_cells": int((ad.obs[res] == st).sum()),
                 "top_genes": ",".join(genes)})

mdf = pd.DataFrame(rows).sort_values("n_cells", ascending=False).reset_index(drop=True)
mdf.to_csv(OUT, sep="\t", index=False)
print("TAM states:", len(mdf), "->", OUT)
print(mdf[["tam_state","n_cells"]])
