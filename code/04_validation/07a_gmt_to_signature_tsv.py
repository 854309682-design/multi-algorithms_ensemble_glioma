#!/usr/bin/env python3
"""Stage B3 (prep) -- convert selected MSigDB GMT gene sets into the per-study TSV format that
code/04_validation/07_published_signature_comparison.R consumes (gene<TAB>weight).

The P3 benchmark compares our frozen 29-gene signature against OTHER published gene SETS on the
same external cohorts. Per the agreed scope this uses MSigDB gene sets (score-based comparison)
rather than coefficient-bearing published glioma signatures.

Selection: HALLMARK sets (all) plus canonical sets whose names match a keyword list, capped so the
figure stays readable. Writes <SET_NAME>.tsv into <REPO_ROOT>/reference/published_signatures/.
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
import glob

PUB = str(REF_ROOT / "published_signatures")
GMTS = sorted(glob.glob(os.path.join(PUB, "*.gmt")))
KEY = re.compile(r"GLIOMA|GLIOBLASTOMA|HYPOXIA|GLYCOLYSIS|MESENCHYMAL|INFLAMMATORY|INTERFERON|"
                 r"^HALLMARK_E2F_TARGETS|^HALLMARK_G2M_CHECKPOINT|ANGIOGENESIS|EXTRACELLULAR_MATRIX|"
                 r"^REACTOME_CELL_CYCLE|ANTIGEN_PRESENTATION|IL6_JAK_STAT3|TGF_BETA|COMPLEMENT")
MAX_CANONICAL = 6

hits = {}
for p in GMTS:
    tag = "hallmark" if os.path.basename(p).startswith("h.all") else "canonical"
    n_can = 0
    with open(p) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            name, genes = parts[0], [g for g in parts[2:] if g]
            if tag == "hallmark":
                hits[name] = (genes, "hallmark")
            elif KEY.search(name) and n_can < MAX_CANONICAL:
                hits[name] = (genes, "canonical")
                n_can += 1

# drop any stale converted sets, then write
for f in glob.glob(os.path.join(PUB, "*_MSIGDB.tsv")):
    os.remove(f)
for name, (genes, tag) in sorted(hits.items()):
    safe = re.sub(r"[^A-Za-z0-9_]+", "_", name)[:80]
    with open(os.path.join(PUB, f"{safe}_MSIGDB.tsv"), "w") as out:
        out.write("gene\tweight\n")
        for g in genes:
            out.write(f"{g}\t1\n")
    print(f"{tag:9s} {len(genes):4d} genes  -> {safe}_MSIGDB.tsv")

print(f"\nconverted {len(hits)} gene sets from {len(GMTS)} GMT file(s)")
