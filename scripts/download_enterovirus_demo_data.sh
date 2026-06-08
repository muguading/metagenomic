#!/usr/bin/env bash
set -euo pipefail

ROOT="/Users/wuhhh/Desktop/徐老师/代码/metagenomic"
OUT_DIR="${ROOT}/demo_data/enterovirus"

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

download "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR235/009/ERR2352279/ERR2352279_1.fastq.gz" "${OUT_DIR}/ERR2352279_1.fastq.gz"
download "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR235/009/ERR2352279/ERR2352279_2.fastq.gz" "${OUT_DIR}/ERR2352279_2.fastq.gz"
download "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR235/000/ERR2352270/ERR2352270_1.fastq.gz" "${OUT_DIR}/ERR2352270_1.fastq.gz"
download "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR235/000/ERR2352270/ERR2352270_2.fastq.gz" "${OUT_DIR}/ERR2352270_2.fastq.gz"
download "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR268/051/SRR26856151/SRR26856151_1.fastq.gz" "${OUT_DIR}/SRR26856151_1.fastq.gz"
download "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR268/051/SRR26856151/SRR26856151_2.fastq.gz" "${OUT_DIR}/SRR26856151_2.fastq.gz"
download "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR108/093/SRR10841593/SRR10841593_1.fastq.gz" "${OUT_DIR}/SRR10841593_1.fastq.gz"
download "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR108/093/SRR10841593/SRR10841593_2.fastq.gz" "${OUT_DIR}/SRR10841593_2.fastq.gz"

echo "[done] enterovirus demo data downloaded to ${OUT_DIR}"
