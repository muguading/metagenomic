#!/usr/bin/env bash
set -euo pipefail

ROOT="/Users/wuhhh/Desktop/徐老师/代码/metagenomic"
OUT_DIR="${ROOT}/demo_data/orthohantavirus_demo"
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

RUN_ACCESSION="SRR22027719"
STUDY_ACCESSION="PRJNA893834"
SAMPLE_ACCESSION="SAMN31432885"
SCIENTIFIC_NAME="Orthohantavirus hantanense"
HOST="Apodemus agrarius"
COUNTRY="South Korea"
INSTRUMENT_MODEL="Illumina MiSeq"
LIBRARY_LAYOUT="PAIRED"
READ1_URL="https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR220/019/SRR22027719/SRR22027719_1.fastq.gz"
READ2_URL="https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR220/019/SRR22027719/SRR22027719_2.fastq.gz"
READ1_BYTES="3576878"
READ2_BYTES="3713731"
READ_COUNT="105499"
BASE_COUNT="31052335"
SOURCE_DB="ENA"

download "${READ1_URL}" "${NGS_DIR}/${RUN_ACCESSION}_1.fastq.gz"
download "${READ2_URL}" "${NGS_DIR}/${RUN_ACCESSION}_2.fastq.gz"

{
  printf "run_accession\tstudy_accession\tsample_accession\tscientific_name\thost\tcountry\tinstrument_model\tlibrary_layout\tfastq_1_url\tfastq_2_url\tfastq_1_bytes\tfastq_2_bytes\tread_count\tbase_count\tsource_db\tnote\n"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${RUN_ACCESSION}" \
    "${STUDY_ACCESSION}" \
    "${SAMPLE_ACCESSION}" \
    "${SCIENTIFIC_NAME}" \
    "${HOST}" \
    "${COUNTRY}" \
    "${INSTRUMENT_MODEL}" \
    "${LIBRARY_LAYOUT}" \
    "${READ1_URL}" \
    "${READ2_URL}" \
    "${READ1_BYTES}" \
    "${READ2_BYTES}" \
    "${READ_COUNT}" \
    "${BASE_COUNT}" \
    "${SOURCE_DB}" \
    "Public Illumina paired-end Hantaan orthohantavirus amplicon demo run from ENA; chosen to better match the current orthohantavirus typing workflow and stay lightweight for smoke testing."
} > "${MANIFEST}"

echo "[done] orthohantavirus demo data downloaded to ${NGS_DIR}"
