#!/bin/zsh
set -euo pipefail

ROOT="${1:-/Users/wuhhh/Desktop/徐老师/代码/metagenomic}"
TMP_DIR="${KB_TAXONOMY_TMP_DIR:-$ROOT/tmp/kb_ncbi_taxonomy}"
QUERY_CACHE="$TMP_DIR/query_cache"
LINEAGE_CACHE="$TMP_DIR/lineage_cache"
MANIFEST="$TMP_DIR/manifest.json"
PENDING_IDS="$TMP_DIR/pending_lineage_ids.json"
MAPPING_OUT="$ROOT/database/knowledge_base/mappings/ncbi_taxonomy_links.json"

mkdir -p "$QUERY_CACHE" "$LINEAGE_CACHE"

fetch_json() {
  local url="$1"
  local target="$2"
  local tmp="${target}.tmp"
  rm -f "$tmp"
  if curl --connect-timeout 10 --max-time 40 --retry 2 --retry-all-errors -sS -L "$url" > "$tmp" && jq empty "$tmp" >/dev/null 2>&1; then
    mv "$tmp" "$target"
    return 0
  fi
  rm -f "$tmp"
  return 1
}

python3 "$ROOT/scripts/plan_kb_taxonomy_queries.py" \
  --project-root "$ROOT" \
  --manifest-out "$MANIFEST"

jq -r '.queries[] | [.key, .query] | @tsv' "$MANIFEST" | while IFS=$'\t' read -r key query; do
  target="$QUERY_CACHE/$key.json"
  if [[ -s "$target" ]]; then
    if jq empty "$target" >/dev/null 2>&1; then
      continue
    fi
    rm -f "$target"
  fi
  encoded="$(printf '%s' "$query" | jq -sRr @uri)"
  fetch_json "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon/${encoded}" "$target" || true
done

python3 "$ROOT/scripts/finalize_kb_ncbi_taxonomy.py" \
  --manifest "$MANIFEST" \
  --query-cache-dir "$QUERY_CACHE" \
  --lineage-cache-dir "$LINEAGE_CACHE" \
  --mapping-out "$MAPPING_OUT" \
  --pending-lineage-out "$PENDING_IDS"

if [[ -s "$PENDING_IDS" ]] && [[ "$(jq 'length' "$PENDING_IDS")" -gt 0 ]]; then
  jq -r '.[]' "$PENDING_IDS" | xargs -n 80 | while read -r batch; do
    key="$(printf '%s' "$batch" | shasum -a 1 | awk '{print $1}')"
    target="$LINEAGE_CACHE/$key.json"
    if [[ -s "$target" ]]; then
      if jq empty "$target" >/dev/null 2>&1; then
        continue
      fi
      rm -f "$target"
    fi
    ids="${batch// /,}"
    fetch_json "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon/${ids}" "$target" || true
  done
fi

python3 "$ROOT/scripts/finalize_kb_ncbi_taxonomy.py" \
  --manifest "$MANIFEST" \
  --query-cache-dir "$QUERY_CACHE" \
  --lineage-cache-dir "$LINEAGE_CACHE" \
  --mapping-out "$MAPPING_OUT" \
  --pending-lineage-out "$PENDING_IDS"
