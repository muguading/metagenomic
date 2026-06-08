#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT_DIR="${PROJECT_ROOT}/demo_data/chikv_demo"

mkdir -p "${OUT_DIR}"

download() {
  local url="$1"
  local dest="$2"
  local tmp="${dest}.part"
  if [[ -s "${dest}" ]]; then
    echo "[skip] ${dest}"
    return 0
  fi
  echo "[download] ${dest}"
  rm -f "${tmp}"
  curl -sS -L --fail --retry 3 --retry-delay 2 "${url}" -o "${tmp}"
  mv "${tmp}" "${dest}"
}

download "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR313/003/ERR3131363/ERR3131363_1.fastq.gz" "${OUT_DIR}/ERR3131363_1.fastq.gz"
download "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR313/003/ERR3131363/ERR3131363_2.fastq.gz" "${OUT_DIR}/ERR3131363_2.fastq.gz"

echo "[ok] chikv demo data ready in ${OUT_DIR}"
