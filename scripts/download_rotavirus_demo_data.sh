#!/usr/bin/env bash
set -euo pipefail

ROOT="/Users/wuhhh/Desktop/徐老师/代码/metagenomic"
OUT_DIR="${ROOT}/demo_data/rotavirus"

mkdir -p "${OUT_DIR}"

download() {
  local url="$1"
  local target="$2"
  if [[ -s "${target}" ]]; then
    echo "[skip] ${target}"
    return
  fi
  echo "[download] ${target}"
  curl -L --fail --retry 3 --retry-delay 2 -o "${target}" "${url}"
}

download "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR114/050/SRR11403250/SRR11403250_1.fastq.gz" "${OUT_DIR}/SRR11403250_1.fastq.gz"
download "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR114/050/SRR11403250/SRR11403250_2.fastq.gz" "${OUT_DIR}/SRR11403250_2.fastq.gz"
download "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR820/008/SRR8202298/SRR8202298_1.fastq.gz" "${OUT_DIR}/SRR8202298_1.fastq.gz"
download "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR820/008/SRR8202298/SRR8202298_2.fastq.gz" "${OUT_DIR}/SRR8202298_2.fastq.gz"
download "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR139/071/SRR13910671/SRR13910671_1.fastq.gz" "${OUT_DIR}/SRR13910671_1.fastq.gz"
download "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR139/071/SRR13910671/SRR13910671_2.fastq.gz" "${OUT_DIR}/SRR13910671_2.fastq.gz"

echo "[done] rotavirus demo data downloaded to ${OUT_DIR}"
