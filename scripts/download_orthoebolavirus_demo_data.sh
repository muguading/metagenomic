#!/usr/bin/env bash
set -euo pipefail

ROOT="/Users/wuhhh/Desktop/徐老师/代码/metagenomic"
OUT_DIR="${ROOT}/demo_data/orthoebolavirus_demo"
MANIFEST="${OUT_DIR}/manifest.tsv"

mkdir -p "${OUT_DIR}"

download() {
  local url="$1"
  local target="$2"
  if [[ -s "${target}" ]]; then
    echo "[skip] ${target}"
    return
  fi
  echo "[download] ${target}"
  curl -sS -L --fail --retry 3 --retry-delay 2 -o "${target}.part" "${url}"
  mv "${target}.part" "${target}"
}

RUN_ACCESSION="SRR2722780"
STUDY_ACCESSION="PRJNA298842"
SAMPLE_ACCESSION="SAMN04166791"
SCIENTIFIC_NAME="Zaire ebolavirus"
INSTRUMENT_MODEL="Illumina MiSeq"
LIBRARY_STRATEGY="RNA-Seq"
LIBRARY_LAYOUT="PAIRED"
READ1_URL="https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR272/000/SRR2722780/SRR2722780_1.fastq.gz"
READ2_URL="https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR272/000/SRR2722780/SRR2722780_2.fastq.gz"
READ1_BYTES="5380851"
READ2_BYTES="5762511"
READ_COUNT="61715"
BASE_COUNT="17588969"
SOURCE_DB="ENA"

download "${READ1_URL}" "${OUT_DIR}/${RUN_ACCESSION}_1.fastq.gz"
download "${READ2_URL}" "${OUT_DIR}/${RUN_ACCESSION}_2.fastq.gz"

{
  printf "run_accession\tstudy_accession\tsample_accession\tscientific_name\tinstrument_model\tlibrary_strategy\tlibrary_layout\tfastq_1_url\tfastq_2_url\tfastq_1_bytes\tfastq_2_bytes\tread_count\tbase_count\tsource_db\tnote\n"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${RUN_ACCESSION}" \
    "${STUDY_ACCESSION}" \
    "${SAMPLE_ACCESSION}" \
    "${SCIENTIFIC_NAME}" \
    "${INSTRUMENT_MODEL}" \
    "${LIBRARY_STRATEGY}" \
    "${LIBRARY_LAYOUT}" \
    "${READ1_URL}" \
    "${READ2_URL}" \
    "${READ1_BYTES}" \
    "${READ2_BYTES}" \
    "${READ_COUNT}" \
    "${BASE_COUNT}" \
    "${SOURCE_DB}" \
    "Public lightweight EBOV paired-end RNA-Seq run from ENA; intended for Orthoebolavirus local typing and Ebola Nextclade smoke testing."
} > "${MANIFEST}"

echo "[done] orthoebolavirus demo data downloaded to ${OUT_DIR}"
