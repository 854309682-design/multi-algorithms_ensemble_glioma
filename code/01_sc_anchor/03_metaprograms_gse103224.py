#!/usr/bin/env python3
"""
Stage 1 (part 3): NMF meta-programs on the malignant (Tumor) subset of GSE103224.
For each patient, run NMF (sklearn) at K in {4,5,6}, extract top-50 genes per program,
then cluster all patient x program vectors by Jaccard overlap to get meta-programs.
Outputs:
  - <DATA_ROOT>/GSE103224/metaprograms.tsv   (meta-program id, size, top genes, contributing patients)
  - <DATA_ROOT>/GSE103224/programs_all.tsv   (per-patient program: patient, K, program, top genes list)
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

import os, itertools
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.decomposition import NMF
from collections import defaultdict

BASE = str(DATA_ROOT / "GSE103224")
OUT_META = os.path.join(BASE, "metaprograms.tsv")
OUT_PROG = os.path.join(BASE, "programs_all.tsv")

sc.settings.verbosity = 0
np.random.seed(20260909)

adata = sc.read_h5ad(os.path.join(BASE, "sc_annotated.h5ad"))
adata.obs["patient"] = adata.obs["patient"].astype(str)

# malignant cells only
mal = adata[adata.obs["celltype_marker"] == "Tumor"].copy()
print("malignant cells:", mal.shape)

# Gene set: top ~5000 genes by variance within malignant cells, expressed in 5-95% of malignant cells.
# Exclude ribosomal/mitochondrial/housekeeping genes (dominate NMF top-gene lists).
sub = mal.copy()
sym = sub.var["gene_symbol"].astype(str).to_dict()
Xx = sub.X.toarray() if hasattr(sub.X, "toarray") else np.asarray(sub.X)
expr_frac = (Xx > 0).mean(axis=0)
var_arr = Xx.var(axis=0)
HK = ("RPL","RPS","MT-","MRPL","MRPS","MALAT1","GAPDH","ACTB","ACTG1","TMSB4X","TMSB10",
      "UBB","UBC","PPIA","HSPA8","HSPA1A","HSPA1B","H3F3A","EEF1A1","EEF2","B2M")
cand = []
for i in range(sub.n_vars):
    g = sym[sub.var_names[i]]
    if not 0.05 <= expr_frac[i] <= 0.95: continue
    if g is None: continue
    if g.startswith(HK): continue
    if g == "" or g == "nan": continue
    cand.append(i)
cand = sorted(cand, key=lambda i: -var_arr[i])[:5000]
keep_idx = cand
X = Xx[:, keep_idx]
keep_sym = [sym[sub.var_names[i]] for i in keep_idx]
print("genes kept (top-5000 var, 5-95% expr, no housekeeping):", X.shape[1], "of", sub.n_vars)

# per-patient NMF
programs = []  # list of dict(patient, K, prog, genes_list)
K_vals = [4, 5, 6]
for patient in sorted(row for row in mal.obs["patient"].unique()):
    idx = np.where(mal.obs["patient"].values == patient)[0]
    Xp = X[idx, :]
    if Xp.shape[0] < 60:
        print(f"  skip {patient}: only {Xp.shape[0]} malignant cells"); continue
    Xp = (Xp - Xp.min()) + 1e-4  # ensure non-negative, shift
    for K in K_vals:
        nmf = NMF(n_components=K, init="nndsvda", random_state=20260909, max_iter=1000, tol=1e-6)
        W = nmf.fit_transform(Xp)
        H = nmf.components_  # K x genes
        # top-50 genes per program
        for k in range(K):
            top_idx = np.argsort(H[k, :])[::-1][:50]
            genes = [keep_sym[i] for i in top_idx]
            programs.append({"patient": patient, "K": K, "program": f"{patient}_K{K}_p{k}",
                             "genes": genes, "genes_str": ",".join(genes)})
    print(f"  {patient}: NMF K4/5/6 done")

pdf = pd.DataFrame(programs)
pdf.to_csv(OUT_PROG, sep="\t", index=False)
print("programs:", len(programs), "->", OUT_PROG)

# Jaccard meta-program clustering (graph: connect programs with Jaccard>=0.3)
gset = [set(p["genes_str"].split(",")) for p in programs]
G = np.zeros((len(programs), len(programs)))
for i in range(len(programs)):
    for j in range(i+1, len(programs)):
        inter = len(gset[i] & gset[j]); un = len(gset[i] | gset[j])
        jac = inter/un if un else 0
        G[i, j] = G[j, i] = jac

# connected components (simple BFS) at threshold 0.3
thr = 0.4
visited = np.zeros(len(programs), bool)
meta_id = 0
rows = []
for s in range(len(programs)):
    if visited[s]: continue
    comp = [s]; stack = [s]; visited[s] = True
    while stack:
        u = stack.pop()
        for v in range(len(programs)):
            if not visited[v] and G[u, v] >= thr:
                visited[v] = True; comp.append(v); stack.append(v)
    # top genes across component (frequency)
    freq = defaultdict(int)
    for c in comp:
        for g in gset[c]:
            freq[g] += 1
    top = [g for g, _ in sorted(freq.items(), key=lambda x: -x[1])[:40]]
    patients = sorted({programs[c]["patient"] for c in comp})
    rows.append({"meta_program": f"MP{meta_id}", "n_programs": len(comp),
                 "n_patients": len(patients), "patients": ";".join(patients),
                 "top_genes": ",".join(top)})
    meta_id += 1

mdf = pd.DataFrame(rows)
mdf = mdf.sort_values("n_programs", ascending=False).reset_index(drop=True)
mdf.to_csv(OUT_META, sep="\t", index=False)
print("meta-programs:", len(mdf), "->", OUT_META)
print("top meta-programs by size:")
print(mdf[["meta_program","n_programs","n_patients"]].head(10))
