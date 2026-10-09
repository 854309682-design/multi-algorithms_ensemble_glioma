#!/usr/bin/env python3
"""
Build clean single-tab clinical/survival tables for CGGA-693 and CGGA-325
(outside-of-zip, normalized columns) and report Primary vs Recurrent breakdown.
Outputs -> <DATA_ROOT>/CGGA/<name>_clinical_clean.tsv
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

import zipfile, os, pandas as pd

CGGA = str(DATA_ROOT / "CGGA")

def read_clin(name):
    zp = f"{CGGA}/CGGA.mRNAseq_{name}_clinical.20200506.txt.zip"
    zf = zipfile.ZipFile(zp)
    txt = [n for n in zf.namelist() if n.endswith(".txt") and "__MACOSX" not in n][0]
    df = pd.read_csv(zf.open(txt), sep="\t")
    df.columns = [c.strip() for c in df.columns]
    return df

for name in ["693", "325"]:
    df = read_clin(name)
    idcol = df.columns[0]
    df = df.rename(columns={idcol: "CGGA_ID"})
    df["CGGA_ID"] = df["CGGA_ID"].astype(str).str.strip()
    # normalize PRS_type
    prs = df["PRS_type"].astype(str).str.strip().str.lower() if "PRS_type" in df.columns else None
    out = f"{CGGA}/CGGA.mRNAseq_{name}_clinical_clean.tsv"
    df.to_csv(out, sep="\t", index=False)
    print(f"CGGA-{name}: rows={len(df)} cols={len(df.columns)}")
    if prs is not None:
        print("  PRS_type:", prs.value_counts(dropna=False).to_dict())
    print("  ->", out)
