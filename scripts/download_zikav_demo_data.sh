#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT_DIR="${PROJECT_ROOT}/demo_data/zikav_demo"

mkdir -p "${OUT_DIR}"

download() {
  local url="$1"
  local dest="$2"
  if [[ -s "${dest}" ]]; then
    echo "[skip] ${dest}"
    return 0
  fi
  echo "[download] ${dest}"
  curl -sS -L --fail --retry 3 --retry-delay 2 "${url}" -o "${dest}"
}

download "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/011/SRR25433011/SRR25433011_1.fastq.gz" "${OUT_DIR}/SRR25433011_1.fastq.gz"
download "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/011/SRR25433011/SRR25433011_2.fastq.gz" "${OUT_DIR}/SRR25433011_2.fastq.gz"
download "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/017/SRR25433117/SRR25433117_1.fastq.gz" "${OUT_DIR}/SRR25433117_1.fastq.gz"
download "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/017/SRR25433117/SRR25433117_2.fastq.gz" "${OUT_DIR}/SRR25433117_2.fastq.gz"

echo "[ok] zikav demo data ready in ${OUT_DIR}"
