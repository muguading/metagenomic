#!/bin/zsh
set -euo pipefail

IDS_JSON="${1:?ids json required}"
OUT_DIR="${2:?out dir required}"

mkdir -p "$OUT_DIR"

fetch_one() {
  local taxid="$1"
  local target="$OUT_DIR/${taxid}.json"
  local tmp="${target}.tmp"
  if [[ -s "$target" ]]; then
    return 0
  fi
  rm -f "$tmp"
  if curl --connect-timeout 10 --max-time 40 --retry 2 --retry-all-errors -sS -L "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon/${taxid}" > "$tmp" \
    && jq empty "$tmp" >/dev/null 2>&1; then
    mv "$tmp" "$target"
    return 0
  fi
  rm -f "$tmp"
  return 1
}

jq -r '.[]' "$IDS_JSON" | while IFS= read -r taxid; do
  [[ -n "$taxid" ]] || continue
  fetch_one "$taxid" || true
done
