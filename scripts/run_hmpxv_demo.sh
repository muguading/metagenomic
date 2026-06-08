#!/bin/zsh
set -euo pipefail

PROJECT_ROOT="/Users/wuhhh/Desktop/徐老师/代码/metagenomic"
PYTHON_BIN="/opt/homebrew/Caskroom/mambaforge/base/envs/ncov/bin/python3"
INPUT_DIR="/Users/wuhhh/Desktop/徐老师/代码/Baiyiapp_example/hmpxv/ngs"
OUTPUT_DIR="$PROJECT_ROOT/demo_data/hmpxv"

R1_FASTQ="$INPUT_DIR/hmpxv_ngs_R1.fastq.gz"
R2_FASTQ="$INPUT_DIR/hmpxv_ngs_R2.fastq.gz"

REFERENCE_CANDIDATES=(
  "$PROJECT_ROOT/database/virus/nextclade/hMPXV/sequences.fasta"
  "$PROJECT_ROOT/database/nextclade_db/hMPXV/sequences.fasta"
)

DATASET_CANDIDATES=(
  "$PROJECT_ROOT/database/virus/nextclade/hMPXV"
  "$PROJECT_ROOT/database/nextclade_db/hMPXV"
)

TAXONKIT_CANDIDATES=(
  "/opt/homebrew/bin/taxonkit"
  "/usr/local/bin/taxonkit"
  "/opt/homebrew/Caskroom/mambaforge/base/bin/taxonkit"
  "/Users/wuhhh/miniconda3/bin/taxonkit"
)

pick_first_existing() {
  local candidate
  for candidate in "$@"; do
    if [[ -e "$candidate" ]]; then
      printf '%s\n' "$candidate"
      return 0
    fi
  done
  return 1
}

if [[ ! -x "$PYTHON_BIN" ]]; then
  echo "未找到 ncov 环境 Python: $PYTHON_BIN" >&2
  exit 1
fi

if [[ ! -f "$R1_FASTQ" || ! -f "$R2_FASTQ" ]]; then
  echo "未找到猴痘二代数据，期望文件：" >&2
  echo "  $R1_FASTQ" >&2
  echo "  $R2_FASTQ" >&2
  exit 1
fi

REFERENCE_FASTA="$(pick_first_existing "${REFERENCE_CANDIDATES[@]}")" || {
  echo "未找到猴痘默认参考 sequences.fasta" >&2
  printf '尝试过的路径：\n' >&2
  printf '  %s\n' "${REFERENCE_CANDIDATES[@]}" >&2
  exit 1
}

NEXTCLADE_DATASET="$(pick_first_existing "${DATASET_CANDIDATES[@]}")" || {
  echo "未找到猴痘 Nextclade 数据集目录" >&2
  printf '尝试过的路径：\n' >&2
  printf '  %s\n' "${DATASET_CANDIDATES[@]}" >&2
  exit 1
}

if ! command -v taxonkit >/dev/null 2>&1; then
  TAXONKIT_BIN="$(pick_first_existing "${TAXONKIT_CANDIDATES[@]}")" || true
  if [[ -n "${TAXONKIT_BIN:-}" ]]; then
    export PATH="$(dirname "$TAXONKIT_BIN"):$PATH"
  fi
fi

if ! command -v taxonkit >/dev/null 2>&1; then
  echo "未找到 taxonkit，可先安装或把它加入 PATH 后再运行本脚本。" >&2
  echo "常见位置可检查：" >&2
  printf '  %s\n' "${TAXONKIT_CANDIDATES[@]}" >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

export META_NEXTCLADE_HMPXV_DATASET="$NEXTCLADE_DATASET"

cd "$PROJECT_ROOT"

exec "$PYTHON_BIN" Bac_assemble_260112_newformat.py \
  --input "$INPUT_DIR" \
  --analysis_target virus \
  --inputtype fastq \
  --thread 8 \
  --output "$OUTPUT_DIR" \
  --fake_pip 0 \
  --method freebayes \
  --long_type Nanopore \
  --ref "$REFERENCE_FASTA" \
  --gtf nogtf \
  --genome_len 200k \
  --asm_type shortref \
  --polish_times 1 \
  --polish_soft medaka \
  --species "Monkeypox virus" \
  --rmhost norm \
  --runflow "物种鉴定,基因组组装,分型鉴定" \
  --abun 1 \
  --rna 0
