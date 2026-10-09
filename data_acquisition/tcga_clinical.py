#!/usr/bin/env python3
import json, urllib.request

API = "https://api.gdc.cancer.gov/cases"
fields = (
    "case_id,submitter_id,project.project_id,"
    "demographic.gender,demographic.age_at_index,demographic.race,demographic.vital_status,"
    "diagnoses.vital_status,diagnoses.days_to_death,diagnoses.days_to_last_followup,diagnoses.age_at_diagnosis,"
    "diagnoses.diagnosis_id,diagnoses.tissue_or_organ_of_origin,diagnoses.primary_diagnosis,"
    "diagnoses.last_known_disease_status,diagnoses.tumor_stage,diagnoses.morphology,"
    "samples.submitter_id,samples.sample_type,samples.is_ffpe"
)
filters = {"op": "in", "content": {"field": "project.project_id", "value": ["TCGA-GBM", "TCGA-LGG"]}}
req = {"filters": filters, "fields": fields, "size": 2000, "format": "JSON"}
req_obj = urllib.request.Request(API, data=json.dumps(req).encode(), headers={"Content-Type": "application/json"})
with urllib.request.urlopen(req_obj, timeout=90) as r:
    d = json.loads(r.read().decode())
hits = d["data"]["hits"]
print("cases fetched:", len(hits), "of", d["data"]["pagination"]["total"])
import os, pathlib, sys

# repository root = parent of data_acquisition/ ; DATA_ROOT honours GLIOMA_DATA (see config/paths.py)
REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "config"))
from paths import DATA_ROOT  # noqa: E402

out = str(DATA_ROOT / "TCGA" / "tcga_cases_clinical.json")
os.makedirs(os.path.dirname(out), exist_ok=True)
with open(out, "w") as f:
    json.dump(d, f, indent=1)
print("saved ->", out)

# quick summary
from collections import Counter
proj = Counter(h["project"]["project_id"] for h in hits)
print("by project:", dict(proj))
hists = Counter(h["diagnoses"][0].get("primary_diagnosis","") if h.get("diagnoses") else "" for h in hits)
print("top diagnoses:", hists.most_common(6))
