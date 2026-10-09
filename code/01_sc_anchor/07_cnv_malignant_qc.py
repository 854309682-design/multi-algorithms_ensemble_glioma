#!/usr/bin/env python3
"""Stage 1 (CNV cross-check) -- QC / cross-check of the malignant call using the CNV support
already present in the on-disk object.

Context: the original plan called for a CopyKAT (inferCNV fallback) re-run. `copykat` is not
installable from the configured mirror and `infercnv` needs system JAGS, but the surviving
object <DATA_ROOT>/processed/gse103224_integrated.h5ad ALREADY carries a CNV-supported malignant
call (`malignant_origin = "celltype+infercnv_support"`, plus `cnv_score` and 61 `cnv_leiden`
clusters). This script therefore quantifies that CNV support instead of duplicating the caller:
  - CNV score by malignant status and by malignant state (AC/NPC/MES/OPC) vs non-malignant
  - per-patient malignant fraction
  - concordance of the CNV-supported call with the marker/state labels
Outputs: results/cnv_malignant_qc.tsv, results/figures/figS_cnv_support.pdf / .jpg
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
from scipy.stats import mannwhitneyu, kruskal

ROOT = str(REPO_ROOT)
PROC = os.path.join(DATA_ROOT, "processed")
RES = os.path.join(REPO_ROOT, "results")
FIG = os.path.join(RES, "figures")
np.random.seed(20260909)

integ = ad.read_h5ad(os.path.join(PROC, "gse103224_integrated.h5ad"))
st = ad.read_h5ad(os.path.join(PROC, "gse103224_states.h5ad"))
st = st[integ.obs_names]

mal = integ.obs["malignant"].values.astype(int)
cnv = np.asarray(integ.obs["cnv_score"].values, dtype=float)
state = st.obs["state"].astype(str).values
sample = integ.obs["sample"].astype(str).values
origin = str(integ.obs["malignant_origin"].astype(str).unique()[0])
print("malignant_origin:", origin, "| malignant:", int(mal.sum()), "/", len(mal))

u = mannwhitneyu(cnv[mal == 1], cnv[mal == 0], alternative="greater")
groups = ["AC", "NPC", "MES", "OPC", "not_malignant"]
by_state = [cnv[state == g] for g in groups]
kw = kruskal(*by_state)

rows = [{"comparison": "malignant vs non-malignant (CNV score)", "n_a": int((mal == 1).sum()),
         "n_b": int((mal == 0).sum()), "stat": float(u.statistic), "p": float(u.pvalue),
         "note": f"median malignant {np.median(cnv[mal == 1]):.4f} vs non-malignant {np.median(cnv[mal == 0]):.4f}"},
        {"comparison": "across states (AC/NPC/MES/OPC/not_malignant)", "n_a": len(groups),
         "n_b": 0, "stat": float(kw.statistic), "p": float(kw.pvalue),
         "note": "Kruskal-Wallis on CNV score"}]
for g in groups:
    rows.append({"comparison": f"cnv_score {g}", "n_a": int((state == g).sum()), "n_b": 0,
                 "stat": float(np.median(cnv[state == g])), "p": np.nan,
                 "note": f"median CNV score; IQR {np.percentile(cnv[state == g], 25):.3f}-{np.percentile(cnv[state == g], 75):.3f}"})
qc = pd.DataFrame(rows)
pd_tab = pd.crosstab(pd.Series(state, name="state"), pd.Series(sample, name="patient"))
pd_tab.to_csv(os.path.join(RES, "cnv_malignant_qc_by_patient.tsv"), sep="\t")
qc.to_csv(os.path.join(RES, "cnv_malignant_qc.tsv"), sep="\t", index=False)
print(qc.to_string(index=False))
print("\nmalignant fraction per patient:")
print((pd_tab.loc[[g for g in groups if g != "not_malignant"]].sum() / pd_tab.sum()).round(3).to_string())

fig, axes = plt.subplots(1, 3, figsize=(13, 4))
axes[0].hist(cnv[mal == 0], bins=60, alpha=0.75, label="non-malignant", color="#3cb44b")
axes[0].hist(cnv[mal == 1], bins=60, alpha=0.6, label="malignant", color="#e6194B")
axes[0].set_xlabel("CNV score", fontsize=8); axes[0].set_ylabel("cells", fontsize=8)
axes[0].set_title(f"(a)  CNV score by malignant call (MW p={u.pvalue:.1e})", fontsize=9)
axes[0].legend(fontsize=7, frameon=False)

bp = axes[1].boxplot(by_state, tick_labels=groups, patch_artist=True, showfliers=False)
for p, c in zip(bp["boxes"], ["#e6194B", "#4363d8", "#f58231", "#911eb4", "#999999"]):
    p.set_facecolor(c); p.set_alpha(0.7)
axes[1].set_ylabel("CNV score", fontsize=8)
axes[1].set_title(f"(b)  CNV score by state (KW p={kw.pvalue:.1e})", fontsize=9)
axes[1].tick_params(axis="x", labelsize=7)

frac = (pd_tab.loc[[g for g in groups if g != "not_malignant"]].sum() / pd_tab.sum())
axes[2].bar(range(len(frac)), frac.values, color="steelblue")
axes[2].set_xticks(range(len(frac))); axes[2].set_xticklabels(frac.index, fontsize=7)
axes[2].set_ylim(0, 1.05); axes[2].set_ylabel("fraction malignant_candidate", fontsize=8)
axes[2].set_title("(c)  Malignant fraction per patient", fontsize=9)

plt.tight_layout()
plt.savefig(os.path.join(FIG, "figS_cnv_support.pdf"), dpi=300, bbox_inches="tight")
plt.savefig(os.path.join(FIG, "figS_cnv_support.jpg"), dpi=300, bbox_inches="tight")
plt.close("all")
print("saved figS_cnv_support.pdf/.jpg")
