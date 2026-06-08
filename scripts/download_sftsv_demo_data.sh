#!/usr/bin/env bash
set -euo pipefail

ROOT="/Users/wuhhh/Desktop/徐老师/代码/metagenomic"
OUT_DIR="${ROOT}/demo_data/sftsv_demo_single"
NGS_DIR="${OUT_DIR}/ngs"
MANIFEST="${OUT_DIR}/manifest.tsv"

mkdir -p "${NGS_DIR}"

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

RUN_ACCESSION="SRR14663545"
STUDY_ACCESSION="PRJNA732992"
SAMPLE_ACCESSION="SAMN19356952"
SCIENTIFIC_NAME="Severe fever with thrombocytopenia syndrome virus"
HOST="Homo sapiens"
COUNTRY="China: Wuhan"
INSTRUMENT_MODEL="NextSeq 550"
LIBRARY_LAYOUT="SINGLE"
READ1_URL="https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR146/045/SRR14663545/SRR14663545.fastq.gz"
READ2_URL=""
READ1_BYTES="57796487"
READ2_BYTES=""
READ_COUNT="1880523"
BASE_COUNT="93331568"
SOURCE_DB="ENA"

download "${READ1_URL}" "${NGS_DIR}/${RUN_ACCESSION}.fastq.gz"

{
  printf "run_accession\tstudy_accession\tsample_accession\tscientific_name\thost\tcountry\tinstrument_model\tlibrary_layout\tfastq_url\tfastq_bytes\tread_count\tbase_count\tsource_db\tnote\n"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${RUN_ACCESSION}" \
    "${STUDY_ACCESSION}" \
    "${SAMPLE_ACCESSION}" \
    "${SCIENTIFIC_NAME}" \
    "${HOST}" \
    "${COUNTRY}" \
    "${INSTRUMENT_MODEL}" \
    "${LIBRARY_LAYOUT}" \
    "${READ1_URL}" \
    "${READ1_BYTES}" \
    "${READ_COUNT}" \
    "${BASE_COUNT}" \
    "${SOURCE_DB}" \
    "Public Illumina single-end SFTSV demo run from ENA; lighter than the paired-end NovaSeq sample."
} > "${MANIFEST}"

echo "[done] SFTSV demo data downloaded to ${NGS_DIR}"
