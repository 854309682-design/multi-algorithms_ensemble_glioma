#!/usr/bin/env python3
"""Fig 2 -- malignant states + single-cell candidate gene set (GSE103224).

The original NMF meta-program tables (<DATA_ROOT>/GSE103224/metaprograms.tsv,
programs_all.tsv) and sc_annotated.h5ad are no longer on disk. This panel is therefore
built from the surviving on-disk objects: malignant state labels (states.h5ad, CNV-supported
malignant call), the state-tagged candidate gene table, and the cross-data reproducibility
table produced by the GSE131928 NMF comparison.

Panel A  malignant-cell UMAP by state (AC / NPC / MES / OPC)
Panel B  candidate genes per source program (parsed from the candidate table `detail`)
Panel C  candidate genes by number of supporting programs
Panel D  cross-data reproducibility (top-gene Jaccard vs GSE131928) per meta-program
Outputs: results/figures/fig2_states_candidates.pdf/.jpg, results/fig2_candidate_program_counts.tsv
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
import re
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = str(REPO_ROOT)
RES = os.path.join(REPO_ROOT, "results")
FIG = os.path.join(RES, "figures")
SEED = 20260909
np.random.seed(SEED)

st = ad.read_h5ad(os.path.join(DATA_ROOT, "processed", "gse103224_states.h5ad"))
cand = pd.read_csv(os.path.join(RES, "candidate_genes_single_cell.tsv"), sep="\t")
cross = pd.read_csv(os.path.join(RES, "meta_program_crossdata_consistency.tsv"), sep="\t")

mal = st.obs["malignant"].values.astype(int) == 1
state = st.obs["state"].astype(str).values
umap = np.asarray(st.obsm["X_umap"])
print("malignant cells:", int(mal.sum()), "| states:", dict(pd.Series(state[mal]).value_counts()))

# --- B: genes per source program (tokenise `detail`) ---
prog_counts = {}
for detail in cand["detail"].astype(str):
    for tok in re.split(r"[;,]", detail):
        tok = tok.strip()
        if tok:
            prog_counts[tok] = prog_counts.get(tok, 0) + 1
pc = pd.Series(prog_counts).sort_values(ascending=False)
pc.to_frame("n_candidate_genes").to_csv(os.path.join(RES, "fig2_candidate_program_counts.tsv"), sep="\t")
print("source programs:", len(pc), "| top:", dict(pc.head(6)))

colors = {"AC": "#e6194B", "NPC": "#4363d8", "MES": "#f58231", "OPC": "#911eb4"}
fig, axes = plt.subplots(2, 2, figsize=(11, 9))

ax = axes[0, 0]
for s in ["AC", "NPC", "MES", "OPC"]:
    m = state == s
    ax.scatter(umap[m, 0], umap[m, 1], s=1.2, c=colors[s], linewidths=0, label=f"{s} (n={int(m.sum()):,})", rasterized=False)
ax.set_title("(e)  Malignant states (CNV-supported malignant call)", fontsize=9)
ax.legend(markerscale=8, fontsize=6, frameon=False)
ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("UMAP1", fontsize=7); ax.set_ylabel("UMAP2", fontsize=7)
for s in ax.spines.values():
    s.set_visible(False)

ax = axes[0, 1]
top = pc.head(18)[::-1]
ax.barh(range(len(top)), top.values, color="steelblue")
ax.set_yticks(range(len(top))); ax.set_yticklabels(top.index, fontsize=6)
ax.set_xlabel("candidate genes", fontsize=7)
ax.set_title("(f)  Candidate genes per source program", fontsize=9)

ax = axes[1, 0]
ns = cand["n_signals"].values
ax.hist(ns, bins=range(1, int(ns.max()) + 2), color="darkorange", edgecolor="white")
ax.set_xlabel("number of supporting programs (n_signals)", fontsize=7)
ax.set_ylabel("candidate genes", fontsize=7)
ax.set_title(f"(g)  Candidate set (n={len(cand)} genes)", fontsize=9)

ax = axes[1, 1]
cc = cross.sort_values(cross.columns[2], ascending=False)
ax.bar(range(len(cc)), cc.iloc[:, 2].values, color="seagreen")
ax.axhline(cc.iloc[:, 2].mean(), ls="--", c="#666666", lw=1)
ax.set_xticks([]); ax.set_xlim(-1, len(cc))
ax.set_xlabel(f"meta-programs (n={len(cc)}, sorted)", fontsize=7)
ax.set_ylabel("top-gene Jaccard vs GSE131928", fontsize=7)
ax.set_title(f"(h)  Cross-data reproducibility (mean {cc.iloc[:, 2].mean():.3f})", fontsize=9)

plt.tight_layout()
plt.savefig(os.path.join(FIG, "fig2_states_candidates.pdf"), dpi=300, bbox_inches="tight")
plt.savefig(os.path.join(FIG, "fig2_states_candidates.jpg"), dpi=300, bbox_inches="tight")
plt.close("all")
print("saved fig2_states_candidates.pdf/.jpg")
