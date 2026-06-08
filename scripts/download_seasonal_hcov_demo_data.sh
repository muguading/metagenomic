#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT_DIR="${PROJECT_ROOT}/demo_data/seasonal_hcov"

mkdir -p "${OUT_DIR}"

download() {
  local url="$1"
  local dest="$2"
  if [[ -s "${dest}" ]]; then
    echo "[skip] ${dest}"
    return 0
  fi
  echo "[download] ${dest}"
  curl -sS -L "${url}" -o "${dest}"
}

download "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR191/010/SRR19134710/SRR19134710_1.fastq.gz" "${OUT_DIR}/SRR19134710_1.fastq.gz"
download "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR191/010/SRR19134710/SRR19134710_2.fastq.gz" "${OUT_DIR}/SRR19134710_2.fastq.gz"
download "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR325/002/SRR32515202/SRR32515202_1.fastq.gz" "${OUT_DIR}/SRR32515202_1.fastq.gz"
download "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR325/002/SRR32515202/SRR32515202_2.fastq.gz" "${OUT_DIR}/SRR32515202_2.fastq.gz"
download "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR191/011/SRR19134711/SRR19134711_1.fastq.gz" "${OUT_DIR}/SRR19134711_1.fastq.gz"
download "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR191/011/SRR19134711/SRR19134711_2.fastq.gz" "${OUT_DIR}/SRR19134711_2.fastq.gz"
download "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR162/057/SRR16227957/SRR16227957_1.fastq.gz" "${OUT_DIR}/SRR16227957_1.fastq.gz"
download "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR162/057/SRR16227957/SRR16227957_2.fastq.gz" "${OUT_DIR}/SRR16227957_2.fastq.gz"

echo "[ok] seasonal hcov demo data ready in ${OUT_DIR}"
