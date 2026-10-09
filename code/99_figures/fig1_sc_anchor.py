#!/usr/bin/env python3
"""Fig 1 -- single-cell discovery anchor (GSE103224).

Panels: UMAP by cell class, UMAP by patient, UMAP by CNV-supported malignant call,
stacked cell-class composition per patient.

Cell class is reconstructed from the surviving on-disk objects (the original
<DATA_ROOT>/GSE103224/sc_annotated.h5ad and metaprograms.tsv are no longer on disk):
  malignant cells (integrated.h5ad `malignant` == 1, origin = celltype+infercnv_support)
    -> state label from states.h5ad (AC / NPC / MES / OPC)
  non-malignant cells -> argmax of score_immune / score_endothelial / score_oligo
    (identical rule to code/07_state_composition/01_build_reference.py)
Outputs: results/figures/fig1_sc_anchor.pdf / .jpg, results/fig1_cell_class_counts.tsv
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
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = str(REPO_ROOT)
PROC = os.path.join(DATA_ROOT, "processed")
FIG = os.path.join(REPO_ROOT, "results", "figures")
os.makedirs(FIG, exist_ok=True)
SEED = 20260909
np.random.seed(SEED)

integ = ad.read_h5ad(os.path.join(PROC, "gse103224_integrated.h5ad"))
st = ad.read_h5ad(os.path.join(PROC, "gse103224_states.h5ad"))
st = st[integ.obs_names]

mal = integ.obs["malignant"].values.astype(int)
state = st.obs["state"].astype(str).values
scores = st.obs[["score_immune", "score_endothelial", "score_oligo"]].to_numpy()
nonmal = np.array(["Immune", "Endothelial", "Oligodendrocyte"])[np.argmax(scores, axis=1)]
cls = np.where(mal == 1, state, nonmal)

umap = np.asarray(st.obsm["X_umap"])
sample = integ.obs["sample"].astype(str).values

classes = ["AC", "NPC", "MES", "OPC", "Immune", "Endothelial", "Oligodendrocyte"]
colors = ["#e6194B", "#4363d8", "#f58231", "#911eb4", "#3cb44b", "#42d4f4", "#808000"]
cmap = {c: colors[i] for i, c in enumerate(classes)}

tab = pd.crosstab(pd.Series(cls, name="cell_class"), pd.Series(sample, name="sample"))
tab.to_csv(os.path.join(REPO_ROOT, "results", "fig1_cell_class_counts.tsv"), sep="\t")
print(tab)
print("class totals:", dict(pd.Series(cls).value_counts()))

# --- panels A-C: UMAP ---
def umap_panel(ax, labels, title, lut=None, cmap_obj="tab20"):
    ax.set_facecolor("white")
    if lut is not None:
        for lab in lut:
            m = labels == lab
            ax.scatter(umap[m, 0], umap[m, 1], s=1.2, c=lut[lab], linewidths=0, label=lab, rasterized=False)
        ax.legend(markerscale=8, fontsize=6, frameon=False, loc="center left", bbox_to_anchor=(1, 0.5))
    else:
        cats = sorted(set(labels))
        cm = plt.get_cmap(cmap_obj, len(cats))
        for i, lab in enumerate(cats):
            m = labels == lab
            ax.scatter(umap[m, 0], umap[m, 1], s=1.2, c=[cm(i)], linewidths=0, label=lab, rasterized=False)
        ax.legend(markerscale=8, fontsize=6, frameon=False, loc="center left", bbox_to_anchor=(1, 0.5))
    ax.set_title(title, fontsize=9)
    ax.set_xlabel("UMAP1", fontsize=7)
    ax.set_ylabel("UMAP2", fontsize=7)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)

fig, axes = plt.subplots(2, 2, figsize=(11, 9))
umap_panel(axes[0, 0], cls, f"(a)  Cell class (n={len(cls):,})", lut=cmap)
umap_panel(axes[0, 1], sample, "(b)  Patient")
umap_panel(axes[1, 0], np.where(mal == 1, "malignant_candidate", "non-malignant"),
           "(c)  CNV-supported malignant call")

ax = axes[1, 1]
order = sorted(set(sample))
bottom = np.zeros(len(order))
for c in classes:
    vals = np.array([np.sum((cls == c) & (sample == s)) for s in order])
    ax.bar(order, vals, bottom=bottom, color=cmap[c], label=c, width=0.8)
    bottom += vals
ax.set_title("(d)  Cell class composition per patient", fontsize=9)
ax.set_ylabel("cells", fontsize=7)
ax.legend(markerscale=1, fontsize=6, frameon=False, ncol=2)
plt.tight_layout()
plt.savefig(os.path.join(FIG, "fig1_sc_anchor.pdf"), dpi=300, bbox_inches="tight")
plt.savefig(os.path.join(FIG, "fig1_sc_anchor.jpg"), dpi=300, bbox_inches="tight")
plt.close("all")
print("saved fig1_sc_anchor.pdf/.jpg")
