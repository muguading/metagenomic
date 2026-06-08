#!/usr/bin/env bash

set -euo pipefail

ROOT="/Users/wuhhh/Desktop/徐老师/代码/metagenomic"
OUT_DIR="${ROOT}/demo_data/astroviridae"
RAW_DIR="${OUT_DIR}/raw_interleaved"
MANIFEST="${OUT_DIR}/manifest.tsv"

mkdir -p "${RAW_DIR}"

RUNS=(
  "SRR8444451"
  "SRR8444452"
  "SRR8444453"
  "SRR8444454"
  "SRR8444455"
  "SRR8444456"
  "SRR8444457"
  "SRR8444458"
)

download() {
  local url="$1"
  local target="$2"
  if [[ -s "${target}" ]]; then
    echo "Exists: ${target}"
    return 0
  fi
  echo "Downloading ${url}"
  curl -L --fail --retry 3 --retry-delay 2 -o "${target}" "${url}"
}

fetch_run_meta() {
  local run="$1"
  curl -sS -L "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${run}&result=read_run&fields=run_accession,scientific_name,library_layout,instrument_model,fastq_ftp,fastq_bytes,read_count,base_count"
}

printf "run_accession\tscientific_name\tlibrary_layout\tinstrument_model\tfastq_url\tfastq_bytes\tread_count\tbase_count\traw_fastq\tfastq_mode\tnote\n" > "${MANIFEST}"

for run in "${RUNS[@]}"; do
  meta="$(fetch_run_meta "${run}")"
  line="$(printf '%s\n' "${meta}" | awk 'NR==2')"

  if [[ -z "${line}" ]]; then
    echo "Failed to fetch ENA metadata for ${run}" >&2
    exit 1
  fi

  scientific_name="$(printf '%s\n' "${line}" | awk -F '\t' '{print $2}')"
  library_layout="$(printf '%s\n' "${line}" | awk -F '\t' '{print $3}')"
  instrument_model="$(printf '%s\n' "${line}" | awk -F '\t' '{print $4}')"
  fastq_ftp="$(printf '%s\n' "${line}" | awk -F '\t' '{print $5}')"
  fastq_bytes="$(printf '%s\n' "${line}" | awk -F '\t' '{print $6}')"
  read_count="$(printf '%s\n' "${line}" | awk -F '\t' '{print $7}')"
  base_count="$(printf '%s\n' "${line}" | awk -F '\t' '{print $8}')"

  if [[ -z "${fastq_ftp}" ]]; then
    echo "ENA did not return FASTQ URL for ${run}" >&2
    exit 1
  fi

  fastq_url="https://${fastq_ftp}"
  raw_fastq="${RAW_DIR}/${run}.fastq.gz"
  download "${fastq_url}" "${raw_fastq}"

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${run}" \
    "${scientific_name}" \
    "${library_layout}" \
    "${instrument_model}" \
    "${fastq_url}" \
    "${fastq_bytes}" \
    "${read_count}" \
    "${base_count}" \
    "${raw_fastq}" \
    "single_fastq" \
    "ENA returns one FASTQ with headers ending in /1; use as single-end demo input" >> "${MANIFEST}"
done

echo "Done: ${OUT_DIR}"
