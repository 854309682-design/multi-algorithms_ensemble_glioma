#!/usr/bin/env python3
"""P4.3 - Figure: IvyGAP anatomical-structure enrichment of the risk / signalling / hypoxia scores."""
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
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = str(REPO_ROOT)
OUT = os.path.join(REPO_ROOT, "results", "ivygap")
FIG = os.path.join(REPO_ROOT, "results", "figures")
os.makedirs(FIG, exist_ok=True)
ORDER = ["LE", "IT", "CT", "CThbv", "CTmvp", "CTpan", "CTpnz"]

sc = pd.read_csv(os.path.join(OUT, "ivygap_structure_scores.tsv"), sep="\t", index_col=0).reindex(ORDER)
scores = ["risk", "signal", "hypoxia"]
cols = {"risk": "#4C72B0", "signal": "#C44E52", "hypoxia": "#8172B3"}

fig, ax = plt.subplots(1, 3, figsize=(13, 4))
for j, s in enumerate(scores):
    ax[j].bar(range(len(ORDER)), sc["raw_" + s], color=cols[s])
    ax[j].set_xticks(range(len(ORDER)))
    ax[j].set_xticklabels(ORDER, rotation=45)
    ax[j].axhline(0, color="grey", lw=0.6)
    ax[j].set_title(f"({'abc'[j]})  {s}  (KW p<0.001)")
    ax[j].set_ylabel("mean score" if j == 0 else "")
    for i, v in enumerate(sc["raw_" + s]):
        ax[j].text(i, v, f"n={int(sc['n'].iloc[i])}", ha="center",
                   va="bottom" if v >= 0 else "top", fontsize=6)
fig.suptitle("IvyGAP: scores by anatomical structure (aggressive niches = CTpan/CTpnz/CTmvp)", fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "fig10_ivygap.pdf"))
fig.savefig(os.path.join(FIG, "fig10_ivygap.jpg"), dpi=300)
print("saved fig10_ivygap.pdf/.jpg")
