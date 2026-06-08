#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT_DIR="${PROJECT_ROOT}/demo_data/hepatovirus"

mkdir -p "${OUT_DIR}"

download() {
  local url="$1"
  local dest="$2"
  if [[ -s "${dest}" ]]; then
    echo "[skip] ${dest}"
    return 0
  fi
  echo "[download] ${dest}"
  curl -L --fail --retry 3 --retry-delay 2 -o "${dest}" "${url}"
}

download "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR147/038/ERR14794638/ERR14794638_1.fastq.gz" "${OUT_DIR}/ERR14794638_1.fastq.gz"
download "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR147/038/ERR14794638/ERR14794638_2.fastq.gz" "${OUT_DIR}/ERR14794638_2.fastq.gz"
download "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR158/009/ERR15860409/ERR15860409_1.fastq.gz" "${OUT_DIR}/ERR15860409_1.fastq.gz"
download "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR158/009/ERR15860409/ERR15860409_2.fastq.gz" "${OUT_DIR}/ERR15860409_2.fastq.gz"
download "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR159/066/ERR15993266/ERR15993266_1.fastq.gz" "${OUT_DIR}/ERR15993266_1.fastq.gz"
download "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR159/066/ERR15993266/ERR15993266_2.fastq.gz" "${OUT_DIR}/ERR15993266_2.fastq.gz"

echo "[done] hepatovirus demo data downloaded to ${OUT_DIR}"
