#!/usr/bin/env bash
set -euo pipefail
# Raw-data tree (override with GLIOMA_DATA; default: <DATA_ROOT>/ next to this repository)
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_ROOT="${GLIOMA_DATA:-$(dirname "$REPO_ROOT")/Datasets}"
MANIFEST="$DATA_ROOT/TCGA/tcga_gene_expr_manifest.json"
OUTDIR="$DATA_ROOT/TCGA/GeneExpr"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

# Build download list: file<tab>file_name
LIST=/tmp/tcga_filelist.txt
${PYTHON:-python3} - "$MANIFEST" "$LIST" <<'PY'
import json,sys
mani=json.load(open(sys.argv[1]))
with open(sys.argv[2],"w") as fo:
    for e in mani:
        fo.write(f"{e['file_id']}\t{e['file_name']}\n")
print("list lines:", mani.__len__())
PY

download_one() {
  local fid="$1" fn="$2"
  if [ -s "$fn" ]; then echo "SKIP $fn"; return 0; fi
  # tmp then atomic-ish move
  curl -sS --retry 3 --max-time 200 -o "$fn.tmp" -L "https://api.gdc.cancer.gov/data/$fid"
  if [ -s "$fn.tmp" ]; then mv "$fn.tmp" "$fn"; echo "OK $fn"; else echo "FAIL $fid"; rm -f "$fn.tmp"; fi
}
export -f download_one

cat "$LIST" | xargs -I{} -P 12 bash -c 'IFS=$'\''\t'\'' read -r fid fn <<< "{}"; download_one "$fid" "$fn"'

echo "ALL_DONE count=$(ls -1 | wc -l)"
