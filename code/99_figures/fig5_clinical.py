#!/usr/bin/env python3
"""Fig 5 -- clinical translation plate: composite of the clinico-molecular nomogram,
3-year calibration and 3-year DCA.

Composed from the SOURCE VECTOR PDFs (PyMuPDF show_pdf_page) -- never re-rasterised --
per the combined-figure rule. A 300-dpi JPEG preview is emitted alongside.
Inputs:  results/figures/05b_nomogram.pdf, 05b_calibration_3yr.pdf, 05b_dca_3yr.pdf
Outputs: results/figures/fig5_clinical.pdf / .jpg
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
import fitz

FIG = str(REPO_ROOT / "results" / "figures")
SRC = ["05b_nomogram.pdf", "05b_calibration_3yr.pdf", "05b_dca_3yr.pdf"]

docs = [fitz.open(os.path.join(FIG, f)) for f in SRC]
out = fitz.open()
W, H = 1020.0, 1180.0
page = out.new_page(width=W, height=H)
page.show_pdf_page(fitz.Rect(20, 20, 1000, 692), docs[0], 0)     # nomogram (top, full width)
page.insert_text((30, 46), "(a)", fontsize=16, fontname="hebo")
page.show_pdf_page(fitz.Rect(20, 712, 452, 1144), docs[1], 0)    # calibration (bottom left)
page.insert_text((30, 738), "(b)", fontsize=16, fontname="hebo")
page.show_pdf_page(fitz.Rect(488, 712, 992, 1144), docs[2], 0)   # DCA (bottom right)
page.insert_text((498, 738), "(c)", fontsize=16, fontname="hebo")
out.save(os.path.join(FIG, "fig5_clinical.pdf"))

pix = page.get_pixmap(dpi=300)
pix.save(os.path.join(FIG, "fig5_clinical.jpg"))
print("fig5_clinical.pdf/.jpg written from vector panels:", SRC)
for d in docs:
    d.close()
out.close()
