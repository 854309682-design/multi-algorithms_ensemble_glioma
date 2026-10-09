#!/usr/bin/env python3
"""Stage B3 (prep, published-signature replay) -- extract a curated set of PUBLISHED glioma /
prognostic gene signatures from the MSigDB C2:CGP collection into the per-study TSV format that
code/04_validation/07_published_signature_comparison.R consumes (gene<TAB>weight).

Why C2:CGP: "chemical and genetic perturbations" holds the MSigDB-curated gene lists of the
canonical published signatures, so the benchmark can replay them WITHOUT hand-transcribing a
paper's supplementary table (which is exactly what 07_published_signature_comparison.R warns
against). MSigDB CGP sets are UNWEIGHTED -> weight = 1; the comparison script then scores each set
as a per-gene z-score mean, the standard published-signature replay.

Source: c2.cgp.v2026.1.Hs.symbols.gmt (MSigDB 2026.1.Hs), from
        https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/
Note:  Ragnum-2015 hypoxia is a PROSTATE signature with no MSigDB-curated set (its gene list would
        need the Oncotarget 2015 supplement); Buffa-2010 hypoxia is the second hypoxia metagene
        comparator instead.
Output: <study>.tsv in <REPO_ROOT>/reference/published_signatures/ (one per published signature).
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

PUB = str(REF_ROOT / "published_signatures")
CGP = os.path.join(PUB, "c2.cgp.v2026.1.Hs.symbols.gmt")

# MSigDB set name -> publication label used as the `study` name in the benchmark table
MAP = {
    "VERHAAK_GLIOBLASTOMA_PRONEURAL": "Verhaak-2010_Proneural",
    "VERHAAK_GLIOBLASTOMA_NEURAL": "Verhaak-2010_Neural",
    "VERHAAK_GLIOBLASTOMA_CLASSICAL": "Verhaak-2010_Classical",
    "VERHAAK_GLIOBLASTOMA_MESENCHYMAL": "Verhaak-2010_Mesenchymal",
    "WINTER_HYPOXIA_METAGENE": "Winter-2007_HypoxiaMetagene",
    "BUFFA_HYPOXIA_METAGENE": "Buffa-2010_HypoxiaMetagene",
    "YAMANAKA_GLIOBLASTOMA_SURVIVAL_UP": "Yamanaka-2009_GBM_Survival",
    "BEIER_GLIOMA_STEM_CELL_UP": "Beier-2007_GliomaStemCell",
    "ZHENG_GLIOBLASTOMA_PLASTICITY_UP": "Zheng-2021_GBM_Plasticity",
}

sets = {}
with open(CGP) as fh:
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 3 and parts[0] in MAP:
            sets[parts[0]] = [g for g in parts[2:] if g]

missing = [s for s in MAP if s not in sets]
if missing:
    raise SystemExit(f"!! set(s) not found in {os.path.basename(CGP)}: {missing}")

for label in MAP.values():
    stale = os.path.join(PUB, f"{label}.tsv")
    if os.path.exists(stale):
        os.remove(stale)

for name, label in MAP.items():
    genes = sets[name]
    with open(os.path.join(PUB, f"{label}.tsv"), "w") as out:
        out.write("gene\tweight\n")
        for g in genes:
            out.write(f"{g}\t1\n")
    print(f"{len(genes):4d} genes  {name:34s} -> {label}.tsv")

print(f"\nwrote {len(MAP)} published-signature tables into {PUB}")
