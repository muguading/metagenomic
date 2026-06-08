#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  cat <<'USAGE'
Usage:
  ./run_metagenome_analysis.sh <batch_input.tsv> [output_dir]

Batch TSV columns:
  样本名称	三代数据	二代数据左	二代数据右	物种信息
  sample_001	0	/path/to/sample_001_R1.fastq.gz	/path/to/sample_001_R2.fastq.gz	nolevel

Optional environment variables:
  CONDA_ROOT=/home/hpcdc/miniconda3
  CONDA_ENV=meta_main
  THREADS=10
  RUNFLOW=基因组组装,病毒组装,物种鉴定,元件预测
  RMHOST=norm
  ABUN=1
USAGE
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$SCRIPT_DIR}"
PIPELINE_SCRIPT="${PIPELINE_SCRIPT:-$PROJECT_ROOT/Bac_assemble_260112_newformat.py}"

INPUT_TSV="${1:-${INPUT_TSV:-}}"
OUTPUT_DIR="${2:-${OUTPUT_DIR:-$PROJECT_ROOT/outputs/metagenome_$(date +%Y%m%d_%H%M%S)}}"

CONDA_ROOT="${CONDA_ROOT:-/home/hpcdc/miniconda3}"
CONDA_ENV="${CONDA_ENV:-meta_main}"
THREADS="${THREADS:-10}"
RUNFLOW="${RUNFLOW:-基因组组装,病毒组装,物种鉴定,元件预测}"
RMHOST="${RMHOST:-norm}"
ABUN="${ABUN:-1}"

MIN_LONG_FILT="${MIN_LONG_FILT:-500}"
Q_FILT="${Q_FILT:-10}"
GENOME_LEN="${GENOME_LEN:-4m}"
ASM_TYPE="${ASM_TYPE:-shortasm}"
POLISH_TIMES="${POLISH_TIMES:-1}"
POLISH_SOFT="${POLISH_SOFT:-medaka}"

if [[ -z "$INPUT_TSV" ]]; then
  usage >&2
  exit 2
fi

if [[ ! -f "$INPUT_TSV" ]]; then
  echo "ERROR: batch TSV not found: $INPUT_TSV" >&2
  exit 2
fi

if [[ ! -f "$PIPELINE_SCRIPT" ]]; then
  echo "ERROR: pipeline script not found: $PIPELINE_SCRIPT" >&2
  exit 2
fi

if [[ -x "$CONDA_ROOT/bin/conda" ]]; then
  CONDA_BIN="$CONDA_ROOT/bin/conda"
elif [[ -x "$CONDA_ROOT/condabin/conda" ]]; then
  CONDA_BIN="$CONDA_ROOT/condabin/conda"
elif command -v conda >/dev/null 2>&1; then
  CONDA_BIN="$(command -v conda)"
else
  echo "ERROR: conda not found. Set CONDA_ROOT or make conda available in PATH." >&2
  exit 2
fi

mkdir -p "$OUTPUT_DIR"

echo "== Metagenome analysis =="
echo "Project root : $PROJECT_ROOT"
echo "Input TSV    : $INPUT_TSV"
echo "Output dir   : $OUTPUT_DIR"
echo "Conda        : $CONDA_BIN"
echo "Conda env    : $CONDA_ENV"
echo "Threads      : $THREADS"
echo "Runflow      : $RUNFLOW"
echo

cmd=(
  "$CONDA_BIN" run -n "$CONDA_ENV" --no-capture-output
  python -u "$PIPELINE_SCRIPT"
  --input "$INPUT_TSV"
  --analysis_target bacteria
  --inputtype fastq
  --minlongfilt "$MIN_LONG_FILT"
  --Qfilt "$Q_FILT"
  --barcodekit none
  --thread "$THREADS"
  --output "$OUTPUT_DIR"
  --fake_pip 0
  --method meta
  --long_type Nanopore
  --ref noref
  --gtf nogtf
  --genome_len "$GENOME_LEN"
  --asm_type "$ASM_TYPE"
  --polish_times "$POLISH_TIMES"
  --polish_soft "$POLISH_SOFT"
  --species nolevel
  --rmhost "$RMHOST"
  --runflow "$RUNFLOW"
  --abun "$ABUN"
  --rna 0
)

printf 'Running command:\n'
printf ' %q' "${cmd[@]}"
printf '\n\n'

"${cmd[@]}"
