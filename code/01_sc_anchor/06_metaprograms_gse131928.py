#!/usr/bin/env python3
"""
Stage 1 (part 6): reproduce NMF meta-programs on GSE131928 (Neftel 2019, Smartseq2, IDHwt GBM)
and report cross-data consistency with the GSE103224 meta-programs (Jaccard).
Output: <DATA_ROOT>/GSE131928/metaprograms.tsv + consistency report in results/
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

import os, re
import numpy as np
import pandas as pd
from sklearn.decomposition import NMF
from collections import defaultdict

BASE = str(DATA_ROOT / "GSE131928")
OUT = os.path.join(BASE, "metaprograms.tsv")
F = os.path.join(BASE, "GSM3828672_Smartseq2_GBM_IDHwt_processed_TPM.tsv.gz")
GSE103 = str(DATA_ROOT / "GSE103224" / "metaprograms.tsv")
np.random.seed(20260909)

# load (first col = GENE symbol, rest = cells)
df = pd.read_csv(F, sep="\t")
df = df.set_index("GENE")
df = df.T  # cells x genes
print("GSE131928 shape (cells x genes):", df.shape)

# parse patient from cell name e.g. MGH101-P1-A04 -> MGH101-P1
def patient_of(c):
    m = re.match(r"^(M[A-Za-z0-9]+-\w+)", str(c))
    return m.group(1) if m else "unknown"
patients = [patient_of(c) for c in df.index]

# raw counts are TPM (positive); log1p for NMF (non-negative)
log = np.log1p(df.values.astype("float32"))

# minimal QC: drop near-empty genes (zeros in >95% cells) & housekeeping
sym = df.columns.values
expr_frac = (log > 0).mean(axis=0)
HK = ("RPL","RPS","MT-","MRPL","MRPS","MALAT1","GAPDH","ACTB","ACTG1","TMSB4X","TMSB10","UBB","UBC","PPIA")
keep = []
for i in range(log.shape[1]):
    g = sym[i]
    if not 0.05 <= expr_frac[i] <= 0.95: continue
    if not isinstance(g, str) or g == "": continue
    if g.startswith(HK): continue
    keep.append(i)
Xall = log[:, keep]
keep_sym = list(sym[keep])
print("genes kept:", Xall.shape[1])

# per-patient NMF (use patients with >=60 cells)
K_vals = [4, 5, 6]
programs = []
pats = pd.Series(patients)
for p in sorted(pats.unique()):
    idx = np.where(np.array(patients) == p)[0]
    if len(idx) < 60:
        print("  skip", p, len(idx)); continue
    Xp = (Xall[idx, :] - Xall[idx, :].min()) + 1e-4
    for K in K_vals:
        nmf = NMF(n_components=K, init="nndsvda", random_state=20260909, max_iter=1000, tol=1e-6)
        H = nmf.fit_transform(Xp)  # W (cells x K) -> use components via fit
        comp = nmf.components_  # K x genes
        for k in range(K):
            top_i = np.argsort(comp[k, :])[::-1][:50]
            genes = [keep_sym[i] for i in top_i]
            programs.append({"patient": p, "K": K, "program": f"{p}_K{K}_p{k}", "genes_str": ",".join(genes)})
    print("  ", p, "NMF done  (cells:", len(idx), ")")

pdf = pd.DataFrame(programs)
pdf.to_csv(OUT, sep="\t", index=False)
print("programs:", len(programs), "->", OUT)

# ---- cross-data consistency: GSE103224 meta-programs vs GSE131928 programs ----
mp103 = pd.read_csv(GSE103, sep="\t")
def jac(a, b):
    a=set(a); b=set(b); return len(a&b)/len(a|b) if (a|b) else 0
res = []
for _, r in mp103.iterrows():
    gs = r["top_genes"].split(",")
    best = 0; bestp = ""
    for _, pr in pdf.iterrows():
        j = jac(gs, pr["genes_str"].split(","))
        if j > best: best = j; bestp = pr["program"]
    res.append({"meta_program_103": r["meta_program"], "n_programs": r["n_programs"],
                "top_gene_jaccard_to_GSE131928": round(best, 3), "best_program_131928": bestp})
rep = pd.DataFrame(res).sort_values("top_gene_jaccard_to_GSE131928", ascending=False)
rep.to_csv(str(REPO_ROOT / "results" / "meta_program_crossdata_consistency.tsv"), sep="\t", index=False)
print("cross-data consistency (mean Jaccard):", round(rep["top_gene_jaccard_to_GSE131928"].mean(), 3))
print(rep.head(10))
