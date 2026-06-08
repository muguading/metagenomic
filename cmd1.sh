#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="/Users/wuhhh/Desktop/徐老师/代码/metagenomic"
PYTHON_CANDIDATES=(
  "/opt/homebrew/Caskroom/mambaforge/base/envs/ncov/bin/python3"
  "$PROJECT_ROOT/.venv_web/bin/python"
  "$(command -v python3 2>/dev/null || true)"
)

INPUT_FASTA="${1:-$PROJECT_ROOT/database/virus/HIV/testdata/test.fasta}"
SAMPLE_PREFIX="${2:-hiv_cmd1}"
SPECIES_LABEL="${3:-HIV-1}"
WORKDIR="${4:-$PROJECT_ROOT/tmp_hiv_cmd1}"

pick_first_executable() {
  local candidate
  for candidate in "$@"; do
    if [[ -n "$candidate" && -x "$candidate" ]]; then
      printf '%s\n' "$candidate"
      return 0
    fi
  done
  return 1
}

PYTHON_BIN="$(pick_first_executable "${PYTHON_CANDIDATES[@]}")" || {
  echo "未找到可用 Python，请检查 ncov 环境或 python3。" >&2
  exit 1
}
export PYTHONPATH="$PROJECT_ROOT${PYTHONPATH:+:$PYTHONPATH}"

if [[ ! -f "$INPUT_FASTA" ]]; then
  echo "未找到输入 FASTA: $INPUT_FASTA" >&2
  exit 1
fi

rm -rf "$WORKDIR"
mkdir -p "$WORKDIR"
cp "$INPUT_FASTA" "$WORKDIR/$SAMPLE_PREFIX.consensus.fasta"

cd "$WORKDIR"

"$PYTHON_BIN" - <<PY
from metagenomic_refactor.virus_analysis import resolve_hiv_reference, virus_typing

pre = "${SAMPLE_PREFIX}"
species = "${SPECIES_LABEL}"
fasta = "${WORKDIR}/${SAMPLE_PREFIX}.consensus.fasta"

resolve_hiv_reference(pre, species=species, query_fasta=fasta)
virus_typing(pre, species)

print("HIV workflow finished")
print(f"summary: {pre}_serotype_result.tsv")
print(f"selection: {pre}_hiv_reference_selection/selection.tsv")
print(f"resistance: {pre}_hiv_resistance.tsv")
print(f"subtyping_json: {pre}_hiv_subtyping.json")
PY

echo
echo "结果目录: $WORKDIR"
echo "摘要结果: $WORKDIR/${SAMPLE_PREFIX}_serotype_result.tsv"
