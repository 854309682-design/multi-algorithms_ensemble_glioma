#!/usr/bin/env python3
"""Assemble the submission figure set into manuscript/Figures/ as 6 combined MAIN figures and
4 combined SUPPLEMENTARY figures.

PANEL LETTERING. The multi-panel source plates now carry their own *flat* panel letters
((a), (b), (c) ...) drawn by their figure scripts, so this compositor adds NO further letters for them.
It only letters plates that it builds here from single-panel leaf PDFs: the Fig-2 external-validation
plate (forest + 2 Kaplan-Meier) and the Figure-S4 legacy-clinical plate.

Every combined plate is composed from its SOURCE VECTOR PDFs (PyMuPDF show_pdf_page) -- never
re-rasterised -- with 300-dpi JPEG previews; the source single-panel PDFs are kept in source_panels/.
Each composite carries its FIGURE TITLE (the manuscript legend lead-in) in a bold 13-pt band on top;
the band adds 30 pt to the page height.

Flat lettering per figure:
  Fig1 (a)-(h) | Fig2 (a)-(f) | Fig3 (a)-(f) | Fig4 (a)-(f) | Fig5 (a)-(d) | Fig6 (a)
  FigS1 (a)-(b) | FigS2 (a)-(c) | FigS3 (no letter) | FigS4 (a)-(c)
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
import glob
import shutil
import fitz

ROOT = str(REPO_ROOT)
SRC = os.path.join(ROOT, "results", "figures")
OUT = os.path.join(ROOT, "manuscript", "Figures")
PANELS = os.path.join(OUT, "source_panels")
INNER, GAP, PAD, W = 980.0, 22.0, 24.0, 1020.0
TITLE_H, TITLE_FS = 30.0, 13.0


def compose(rows, out_stem, outdir=OUT, title=None):
    """rows: list of rows; each row a list of (source_stem, width_fraction, letter|None).
    title: optional figure title drawn in a band on top (page height grows by TITLE_H)."""
    stems = [s for row in rows for s, _, _ in row]
    idx = {s: i for i, s in enumerate(stems)}
    docs = [fitz.open(os.path.join(SRC, s + ".pdf")) for s in stems]
    row_w, row_h = [], []
    for row in rows:
        widths = [f * INNER for _, f, _ in row]
        hs = [widths[k] * docs[idx[s]][0].rect.height / docs[idx[s]][0].rect.width
              for k, (s, _, _) in enumerate(row)]
        row_w.append(sum(widths) + GAP * (len(row) - 1))
        row_h.append(max(hs))
    band = TITLE_H if title else 0.0
    H = PAD + band + sum(row_h) + GAP * (len(rows) - 1) + PAD
    out = fitz.open()
    page = out.new_page(width=W, height=H)
    if title:
        page.insert_text((PAD, PAD + TITLE_FS), title, fontsize=TITLE_FS, fontname="hebo")
    y = PAD + band
    for r, row in enumerate(rows):
        x = PAD + (INNER - row_w[r]) / 2.0
        widths = [f * INNER for _, f, _ in row]
        for (s, _, letter), w in zip(row, widths):
            d = docs[idx[s]]
            ph = w * d[0].rect.height / d[0].rect.width
            yoff = (row_h[r] - ph) / 2.0
            page.show_pdf_page(fitz.Rect(x, y + yoff, x + w, y + yoff + ph), d, 0)
            if letter:
                page.insert_text((x + 8, y + yoff + 16), letter, fontsize=16, fontname="hebo")
            x += w + GAP
        y += row_h[r] + GAP
    out.save(os.path.join(outdir, out_stem + ".pdf"))
    page.get_pixmap(dpi=300).save(os.path.join(outdir, out_stem + ".jpg"))
    out.close()
    for d in docs:
        d.close()
    return H


# --- 1. (re)build the Fig-2 external-validation plate (forest + 2 KM) with flat letters (d) / (e) / (f) ----
compose([[("fig4_forest", 1.0, "(d)")],
         [("fig4_km_cgga693", 0.49, "(e)"), ("fig4_km_cgga325", 0.49, "(f)")]],
        "fig4_external_validation", outdir=SRC)

# --- 2. submission spec: (name, rows, title, legend);  letter=None => the plate letters itself ----
# Titles match the manuscript figure-legend lead-ins verbatim; each composite draws its title in a top band.
SPEC = [
    ("Fig1", [[("fig1_sc_anchor", 1.0, None)], [("fig2_states_candidates", 1.0, None)]],
     "Single-cell discovery anchor and candidate set (GSE103224)",
     "(a) UMAP by cell class, (b) by patient, (c) by CNV-supported malignant call, (d) cell-class composition per patient; (e) malignant-state UMAP, (f) candidate genes per source program, (g) the 537-gene candidate set, (h) cross-dataset program reproducibility."),
    ("Fig2", [[("fig3_model_building", 1.0, None)], [("fig4_external_validation", 0.62, None)]],
     "Model building and external validation",
     "(a) univariate screen, (b) LASSO path, (c) learner cross-validation; (d) cohort C-index forest plot, (e-f) per-cohort Kaplan-Meier."),
    ("Fig3", [[("fig5_clinical", 0.58, None)], [("fig6_biology_closure", 1.0, None)]],
     "Clinical translation and biology closure",
     "(a) nomogram, (b) 3-yr calibration, (c) 3-yr DCA; (d) signature-gene importance, (e) single-cell state traceback, (f) GO BP enrichment."),
    ("Fig4", [[("fig8_state_composition", 1.0, None)], [("fig9_signalling", 1.0, None)]],
     "Cell-state composition and immune-to-malignant signalling",
     "(a-c) state-composition model (discrimination, confounding, state vs grade); (d-f) signalling circuit, external discrimination, score HR by chemo."),
    ("Fig5", [[("fig10_ivygap", 1.0, None)], [("fig11_platform_robust", 0.70, None)]],
     "Spatial validation and platform-robust scoring",
     "(a-c) IvyGAP risk / signalling / hypoxia scores by anatomical structure; (d) platform-robust scoring."),
    ("Fig6", [[("fig12_model_comparison", 1.0, None)]],
     "Unified multi-omic model comparison",
     "(a) cross-validated discrimination of clinical, gene, state and signalling models and their combinations."),
    ("FigS1", [[("fig7_sensitivity", 1.0, None)]],
     "Sensitivity analyses",
     "(a) C-index by signature variant and cohort; (b) the primary signature on all vs GBM-only external samples."),
    ("FigS2", [[("figS_cnv_support", 1.0, None)]],
     "Single-cell CNV support",
     "(a) CNV score by malignant call, (b) by malignant state, (c) malignant fraction per patient."),
    ("FigS3", [[("figS_gse16011_grade", 0.55, None)]],
     "GSE16011 cross-platform direction check (GPL8542)",
     "The frozen score rises monotonically with WHO grade (GPL8542)."),
    ("FigS4", [[("05_nomogram", 0.55, "(a)")],
               [("05_calibration_3yr", 0.48, "(b)"), ("05_dca_3yr", 0.48, "(c)")]],
     "Legacy risk + age-only nomogram (sensitivity)",
     "(a) risk + age nomogram; (b) 3-yr calibration; (c) 3-yr DCA -- the legacy clinical model."),
]

os.makedirs(PANELS, exist_ok=True)
for f in glob.glob(os.path.join(OUT, "Fig*.pdf")) + glob.glob(os.path.join(OUT, "Fig*.jpg")):
    os.remove(f)

used = []
for name, rows, title, legend in SPEC:
    disp = name.replace("Fig", "Figure ")
    H = compose(rows, name, title=f"{disp}. {title}.")
    for s in [s for row in rows for s, _, _ in row]:
        if s not in used:
            used.append(s)
    print(f"  {name}.pdf/.jpg  <-  {[s for row in rows for s, _, _ in row]}  (page height {H:.0f} pt)")

for extra in ["fig4_forest", "fig4_km_cgga693", "fig4_km_cgga325"]:
    if extra not in used:
        used.append(extra)
for s in used:
    shutil.copy2(os.path.join(SRC, s + ".pdf"), os.path.join(PANELS, s + ".pdf"))

with open(os.path.join(OUT, "README.md"), "w") as fh:
    fh.write("# Paper 3 - submission figures (6 main + 4 supplementary)\n\n")
    fh.write("Combined figures assembled by `code/99_figures/assemble_manuscript_figures.py`. Each figure is\n")
    fh.write("a **vector PDF** and a 300-dpi **JPEG** preview, composed from source vector PDFs (PyMuPDF\n")
    fh.write("`show_pdf_page`), with a **single flat panel-letter sequence** ((a), (b), (c) ...) per figure and a\n")
    fh.write("**figure-title band on top** (the manuscript legend lead-in, bold 13 pt; +30 pt page height). The\n")
    fh.write("source single-panel PDFs are in `source_panels/`.\n\n")
    fh.write("| Figure | Letters | Panels | Title |\n|---|---|---|---|\n")
    LETTERS = {"Fig1": "(a)\u2013(h)", "Fig2": "(a)\u2013(f)", "Fig3": "(a)\u2013(f)", "Fig4": "(a)\u2013(f)",
               "Fig5": "(a)\u2013(d)", "Fig6": "(a)", "FigS1": "(a)\u2013(b)", "FigS2": "(a)\u2013(c)",
               "FigS3": "—", "FigS4": "(a)\u2013(c)"}
    for name, rows, title, legend in SPEC:
        fh.write(f"| {name} | {LETTERS[name]} | {' + '.join(s for row in rows for s, _, _ in row)} | {title} |\n")
    fh.write("\n## Legends\n\n")
    for name, rows, title, legend in SPEC:
        fh.write(f"**{name}. {title}.** {legend}\n\n")
    fh.write("**Note.** The GSE16011 overall-survival figure is absent (that analysis is blocked on a\n")
    fh.write("closed-access survival table; see `data_acquisition/README.md` section 7).\n")
print("wrote", os.path.join(OUT, "README.md"))
