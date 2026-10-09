#!/usr/bin/env python3
"""
Stage 1 (part 5): compile the single-cell-anchored candidate gene set.
- Load meta-programs (malignant states) + TAM state markers.
- Strip ENSG prefix -> gene symbol.
- Score each meta-program against curated GB state signatures (NPC/OPC/AC/MES/hypoxia/proliferation)
  and assign a state label by best overlap.
- Merge with TAM markers (labeled by TAM state).
- Intersect with the TCGA training gene space; also report overlap with CGGA-693 and the
  microarray platforms (GSE108474 gene-level, GSE16011 after probe mapping is later).
Output: results/candidate_genes_single_cell.tsv  (gene, source_state, n_signals)
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
import pandas as pd
from collections import defaultdict, Counter

BASE = str(DATA_ROOT / "GSE103224")
DATA = str(DATA_ROOT)
OUT = str(REPO_ROOT / "results" / "candidate_genes_single_cell.tsv")

def sym(g):
    # var_names "ENSG00000176410.7_DNAJC30" -> "DNAJC30"
    g = str(g)
    if "_" in g and g.startswith("ENSG"):
        return g.split("_", 1)[1]
    return g

# --- curated state signatures (published-gene sets, Neftel-like + hypoxia + proliferation) ---
SIG = {
    "NPC":  ["SOX2","PAX6","NES","FABP7","PTPRZ1","CD24","VIM","HES1","MKI67","EGFR"],
    "OPC":  ["OLIG1","OLIG2","PDGFRA","CSPG4","SOX8","NKX2-2","BCAN","GPR17","MAG","MBP"],
    "AC":   ["GFAP","AQP4","SLC1A3","ALDH1L1","AGT","SPARCL1","GJA1","APOE","CLU","S100B","FGF14"],
    "MES":  ["CHI3L1","ANXA1","ANXA2","CD44","SERPINE1","IGFBP7","VEGFA","TIMP1","FN1","TNC","LGALS1","ANXA5","TGFBI","PLAUR"],
    "Hypoxia": ["VEGFA","SLC2A1","PGK1","ENO1","ALDOA","TPI1","LDHA","CA9","HILPDA","ADM","BNIP3","NDRG1","P4HA1","HK2"],
    "Prolif": ["STMN1","TUBA1A","TOP2A","MKI67","PTTG1","UBE2C","NUSAP1","CCNB1","CDK1","CENPF","TUBA1B"],
}

# --- load meta-programs ---
mp = pd.read_csv(os.path.join(BASE, "metaprograms.tsv"), sep="\t")
mp["top_genes"] = mp["top_genes"].astype(str)
def sig_score(genes, siggenes):
    gs = set(genes)
    return len(gs & set(siggenes)) / max(1, min(len(siggenes), len(gs)))
def best_state(genes):
    scores = {k: sig_score(genes, v) for k, v in SIG.items()}
    return max(scores, key=scores.get), scores

# gene -> state label (majority vote across programs)
gene_state = defaultdict(list)
for _, r in mp.iterrows():
    genes = [sym(g) for g in r["top_genes"].split(",")]
    state, scores = best_state(genes)
    for g in genes:
        gene_state[g].append(state)

# --- load TAM markers ---
tam = pd.read_csv(os.path.join(BASE, "tam_state_markers.tsv"), sep="\t")
for _, r in tam.iterrows():
    for g in [sym(x) for x in r["top_genes"].split(",")]:
        gene_state[g].append(str(r["tam_state"]))

# --- build candidate table ---
rows = []
for g, states in gene_state.items():
    if not g or g in ("nan",""): continue
    c = Counter(states)
    rows.append({"gene": g, "state": "/".join([s for s,_ in c.most_common(2)] if len(c)>1 else [next(iter(c))]),
                 "n_signals": len(states), "detail": ";".join(states)})
cand = pd.DataFrame(rows)
cand = cand[cand["gene"].str.match(r"^[A-Z0-9]+$") & (~cand["gene"].str.match(r"^(ENSG|RP\d|AC\d|AL\d|CTD|AP\d|LINC|MIR|Y_RNA|U\d|Metazoa)"))].copy()
cand = cand.drop_duplicates("gene").reset_index(drop=True)
print("candidate genes (pre-bulk-intersection):", cand.shape)

# --- intersect with bulk training space (TCGA matrix gene symbols) ---
import gzip
tcga_genes = set()
with gzip.open(os.path.join(DATA, "TCGA/TCGA_GBM_LGG_star_counts_matrix.tsv.gz"), "rt") as fh:
    for line in fh:
        g = line.split("\t", 1)[0]
        tcga_genes.add(g.strip())
        if len(tcga_genes) > 70000: break
print("TCGA bulk gene space:", len(tcga_genes))
cand_b = cand[cand["gene"].isin(tcga_genes)].copy()
cand_b = cand_b.sort_values("n_signals", ascending=False).reset_index(drop=True)
cand_b.to_csv(OUT, sep="\t", index=False)
print("candidate genes intersected with TCGA:", cand_b.shape, "->", OUT)
print("state composition:")
print(cand_b["state"].value_counts().head(15))
