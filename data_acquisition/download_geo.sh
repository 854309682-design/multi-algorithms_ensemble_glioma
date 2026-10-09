#!/usr/bin/env bash
set -euo pipefail
# Raw-data tree (override with GLIOMA_DATA; default: <DATA_ROOT>/ next to this repository)
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEST="${GLIOMA_DATA:-$(dirname "$REPO_ROOT")/Datasets}"
mkdir -p "$DEST"
cd "$DEST"

echo "== GSE131928_RAW.tar =="
curl -sS --max-time 300 -o GSE131928_RAW.tar "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131928/suppl/GSE131928_RAW.tar" \
  && echo OK GSE131928_RAW.tar || echo "FAIL GSE131928_RAW.tar"

echo "== GSE16011_series_matrix =="
curl -sS --max-time 300 -o GSE16011_series_matrix.txt.gz "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE16nnn/GSE16011/matrix/GSE16011_series_matrix.txt.gz" \
  && echo OK GSE16011_series_matrix || echo "FAIL GSE16011_series_matrix"

echo "== GSE108474 files =="
curl -sS --max-time 600 -o GSE108474_REMBRANDT_GeneExpression.txt.gz "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE108nnn/GSE108474/suppl/GSE108474_REMBRANDT_GeneExpression.txt.gz" \
  && echo OK GSE108474_GeneExpression || echo "FAIL GSE108474_GeneExpression"
curl -sS --max-time 120 -o GSE108474_REMBRANDT_clinical.data.txt.gz "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE108nnn/GSE108474/suppl/GSE108474_REMBRANDT_clinical.data.txt.gz" \
  && echo OK GSE108474_clinical || echo "FAIL GSE108474_clinical"
curl -sS --max-time 120 -o GSE108474_REMBRANDT_biospecimen_mapping_GEO.txt.gz "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE108nnn/GSE108474/suppl/GSE108474_REMBRANDT_biospecimen_mapping_GEO.txt.gz" \
  && echo OK GSE108474_biospec || echo "FAIL GSE108474_biospec"

echo "ALL_DONE"
