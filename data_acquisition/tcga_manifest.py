#!/usr/bin/env python3
import json, urllib.request, urllib.parse, sys

API = "https://api.gdc.cancer.gov/files"
fields = "file_id,file_name,file_size,access,cases.case_id,cases.submitter_id,cases.samples.submitter_id,cases.samples.sample_type,cases.project.project_id"
filters = {
    "op": "and",
    "content": [
        {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-GBM", "TCGA-LGG"]}},
        {"op": "in", "content": {"field": "data_type", "value": ["Gene Expression Quantification"]}},
    ],
}
req = {
    "filters": filters,
    "format": "JSON",
    "size": 2000,
    "fields": fields,
}
data = json.dumps(req).encode()
req_obj = urllib.request.Request(API, data=data, headers={"Content-Type": "application/json"})
with urllib.request.urlopen(req_obj, timeout=60) as r:
    d = json.loads(r.read().decode())

hits = d["data"]["hits"]
total = d["data"]["pagination"]["total"]
print("total files:", total, "| returned:", len(hits))
sz = sum(h.get("file_size", 0) for h in hits)
print("sum size MB:", round(sz / 1024 / 1024, 1))

# Print distinct workflow_type values via file_name heuristics
import collections
wf = collections.Counter()
for h in hits:
    fn = h.get("file_name", "")
    if "star_gene_counts" in fn: wf["STAR_counts"] += 1
    elif "htseq" in fn.lower(): wf["HTSeq_counts"] += 1
    else: wf["other"] += 1
print("file_name classes:", dict(wf))

# Save manifest
mani = []
for h in hits:
    cs = h.get("cases", [{}])[0] if h.get("cases") else {}
    samp = (cs.get("samples") or [{}])[0] if cs else {}
    mani.append({
        "file_id": h["file_id"],
        "file_name": h.get("file_name", ""),
        "file_size": h.get("file_size", 0),
        "access": h.get("access"),
        "project": h.get("accepted_projects") or (cs.get("project") or {}).get("project_id"),
        "case_id": cs.get("case_id"),
        "case_submitter_id": cs.get("submitter_id"),
        "sample_submitter_id": samp.get("submitter_id"),
        "sample_type": samp.get("sample_type"),
    })
import os, pathlib, sys

# repository root = parent of data_acquisition/ ; DATA_ROOT honours GLIOMA_DATA (see config/paths.py)
REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "config"))
from paths import DATA_ROOT  # noqa: E402

out = str(DATA_ROOT / "TCGA" / "tcga_gene_expr_manifest.json")
os.makedirs(os.path.dirname(out), exist_ok=True)
with open(out, "w") as f:
    json.dump(mani, f, indent=1)
print("saved manifest ->", out, "entries:", len(mani))
