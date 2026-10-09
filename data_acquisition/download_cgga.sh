#!/usr/bin/env bash
set -euo pipefail
# Raw-data tree (override with GLIOMA_DATA; default: <DATA_ROOT>/ next to this repository)
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEST="${GLIOMA_DATA:-$(dirname "$REPO_ROOT")/Datasets}/CGGA"
mkdir -p "$DEST"
cd "$DEST"
BASE="https://www.cgga.org.cn/download?file="
declare -A files=(
  ["693_RSEM"]="download/20200506/CGGA.mRNAseq_693.RSEM-genes.20200506.txt.zip&type=mRNAseq_693&time=20200506,CGGA.mRNAseq_693.RSEM-genes.20200506.txt.zip"
  ["693_counts"]="download/20220620/CGGA.mRNAseq_693.Read_Counts-genes.20220620.txt.zip&type=mRNAseq_693_counts&time=20220620,CGGA.mRNAseq_693.Read_Counts-genes.20220620.txt.zip"
  ["693_clin"]="download/20200506/CGGA.mRNAseq_693_clinical.20200506.txt.zip&type=mRNAseq_693_clinical&time=20200506,CGGA.mRNAseq_693_clinical.20200506.txt.zip"
  ["325_RSEM"]="download/20200506/CGGA.mRNAseq_325.RSEM-genes.20200506.txt.zip&type=mRNAseq_325&time=20200506,CGGA.mRNAseq_325.RSEM-genes.20200506.txt.zip"
  ["325_counts"]="download/20220620/CGGA.mRNAseq_325.Read_Counts-genes.20220620.txt.zip&type=mRNAseq_325_counts&time=20220620,CGGA.mRNAseq_325.Read_Counts-genes.20220620.txt.zip"
  ["325_clin"]="download/20200506/CGGA.mRNAseq_325_clinical.20200506.txt.zip&type=mRNAseq_325_clinical&time=20200506,CGGA.mRNAseq_325_clinical.20200506.txt.zip"
)
for k in 693_RSEM 693_counts 693_clin 325_RSEM 325_counts 325_clin; do
  IFS=',' read -r q out <<< "${files[$k]}"
  echo "== $k -> $out =="
  curl -sSL --max-time 300 -o "$out" "${BASE}${q}" && echo "OK $out" || echo "FAIL $out"
done
echo "ALL_DONE"
