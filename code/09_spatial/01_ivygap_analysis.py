#!/usr/bin/env python3
"""P4.2 - IvyGAP anatomical-structure validation of the risk / signalling / state scores.

Scores the IvyGAP region-annotated bulk RNA-seq (GSE107559) and tests whether the
single-cell-anchored risk and the immune/malignant signalling scores are higher in
aggressive/hypoxic niches (CTpan, CTpnz, CTmvp) than in leading edge (LE) / infiltrating
tumour (IT). Patient (tumour) effects removed by within-tumour centering.

Outputs (results/ivygap/):
  ivygap_structure_scores.tsv   mean score per structure (raw + tumour-centered)
  ivygap_stats.tsv              Kruskal-Wallis + within-tumour KW + pairwise vs LE
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
import gzip
import numpy as np
import pandas as pd
from scipy import stats

ROOT = str(REPO_ROOT)
IVY = os.path.join(DATA_ROOT, "IvyGAP")
RES = os.path.join(REPO_ROOT, "results")
OUT = os.path.join(RES, "ivygap")
os.makedirs(OUT, exist_ok=True)
REGION_ORDER = ["LE", "IT", "CT", "CThbv", "CTmvp", "CTpan", "CTpnz"]

# --- expression (FPKM) -> symbols ---
expr = pd.read_csv(os.path.join(IVY, "GSE107559_ivygap_fpkm_table.csv.gz"), index_col=0, compression="gzip")
rg = pd.read_csv(os.path.join(IVY, "GSE107559_ivygap_rows-genes.csv.gz"))
id2sym = dict(zip(rg["gene_id"].astype(str), rg["gene_symbol"].astype(str)))
expr.index = [id2sym.get(str(g), str(g)) for g in expr.index]
expr = expr[~expr.index.duplicated(keep="first")]
expr.columns = [str(c) for c in expr.columns]
expr = np.log2(expr + 1)

mp = pd.read_csv(os.path.join(OUT, "ivygap_sample_map.tsv"), sep="\t", dtype=str)
mp = mp[mp.well_id.isin(expr.columns)]
expr = expr[mp.well_id.tolist()]

def mean_z(genes):
    g = [x for x in genes if x in expr.index]
    sub = expr.loc[g]
    z = sub.sub(sub.mean(1), axis=0).div(sub.std(1).replace(0, np.nan), axis=0)
    return z.mean(0), len(g)

sig = pd.read_csv(os.path.join(RES, "signature_pc_coefs.tsv"), sep="\t")
risk, ng = mean_z(sig.gene.tolist())
lr = pd.read_csv(os.path.join(RES, "signalling", "cc_lr_pairs.tsv"), sep="\t")
mal = ["AC", "NPC", "MES", "OPC"]
imm = lr[lr.target.isin(mal) & lr.source.isin(["Immune", "Endothelial"])]
panel = sorted(set(imm.ligand) | set(imm.receptor))
sg, nsg = mean_z([p for p in panel if pd.notna(p)])
HYPOXIA = ["VEGFA","SLC2A1","LDHA","PGK1","PDK1","HK2","BNIP3","NDRG1","ENO1","PDK3","TPI1","P4HA1","ADM","MRPS17","ANKZF1"]
hyp, nh = mean_z(HYPOXIA)
print(f"genes matched: signature {ng}/29, signalling {nsg}/{len(panel)}, hypoxia {nh}/{len(HYPOXIA)}")

d = pd.DataFrame({"risk": risk.values, "signal": sg.values, "hypoxia": hyp.values}, index=expr.columns)
d = d.join(mp.set_index("well_id")[["structure", "tumor_id"]])
d["structure"] = d.structure.astype("category").cat.set_categories(REGION_ORDER, ordered=True)

# tumour-centered scores (remove patient effect)
for c in ["risk", "signal", "hypoxia"]:
    d[c + "_ctr"] = d[c] - d.groupby("tumor_id")[c].transform("mean")

mean_raw = d.groupby("structure", observed=True)[["risk", "signal", "hypoxia"]].mean().round(3)
mean_ctr = d.groupby("structure", observed=True)[["risk_ctr", "signal_ctr", "hypoxia_ctr"]].mean().round(3)
mean_ctr.columns = ["risk", "signal", "hypoxia"]
out = mean_raw.add_prefix("raw_").join(mean_ctr.add_prefix("ctr_"))
out["n"] = d.groupby("structure", observed=True).size()
out.to_csv(os.path.join(OUT, "ivygap_structure_scores.tsv"), sep="\t")
print("\nmean score by structure:\n", out.to_string())

rows = []
for c in ["risk", "signal", "hypoxia"]:
    grps = [d.loc[d.structure == s, c].dropna() for s in REGION_ORDER]
    kw = stats.kruskal(*grps)
    gctr = [d.loc[d.structure == s, c + "_ctr"].dropna() for s in REGION_ORDER]
    kwc = stats.kruskal(*gctr)
    # pairwise vs LE (raw and tumour-centered), CTpan + CTpnz + CTmvp
    for s in ["CTpan", "CTpnz", "CTmvp"]:
        a = d.loc[d.structure == "LE", c].dropna(); b = d.loc[d.structure == s, c].dropna()
        ac = d.loc[d.structure == "LE", c + "_ctr"].dropna(); bc = d.loc[d.structure == s, c + "_ctr"].dropna()
        rows.append({"score": c, "structure": s, "vs": "LE",
                     "MW_p_raw": round(stats.mannwhitneyu(b, a, alternative="greater").pvalue, 4),
                     "MW_p_tumour_centered": round(stats.mannwhitneyu(bc, ac, alternative="greater").pvalue, 4),
                     "mean_vs_LE_raw": round(b.mean() - a.mean(), 3)})
    rows.append({"score": c, "structure": "KW_all", "vs": "-", "MW_p_raw": round(kw.pvalue, 4),
                 "MW_p_tumour_centered": round(kwc.pvalue, 4), "mean_vs_LE_raw": np.nan})
st = pd.DataFrame(rows)
st.to_csv(os.path.join(OUT, "ivygap_stats.tsv"), sep="\t", index=False)
print("\nstructure enrichment stats:\n", st.to_string(index=False))
print("\nDONE")
