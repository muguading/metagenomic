from __future__ import annotations

import argparse
import json
from pathlib import Path


RANK_KEY_MAP = {
    "DOMAIN": "domain",
    "SUPERKINGDOM": "domain",
    "KINGDOM": "kingdom",
    "PHYLUM": "phylum",
    "CLASS": "class",
    "ORDER": "order",
    "FAMILY": "family",
    "GENUS": "genus",
    "SPECIES": "species",
}


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def load_taxonomy_cache(directory: Path) -> dict[str, dict]:
    records: dict[str, dict] = {}
    if not directory.is_dir():
        return records
    for path in sorted(directory.glob("*.json")):
        try:
            payload = load_json(path)
        except (OSError, json.JSONDecodeError):
            continue
        if not isinstance(payload, dict):
            continue
        for node in payload.get("taxonomy_nodes") or []:
            taxonomy = node.get("taxonomy") or {}
            tax_id = str(taxonomy.get("tax_id") or "").strip()
            if tax_id:
                records[tax_id] = taxonomy
    return records


def select_taxonomy(record: dict, query_lookup: dict[str, dict]) -> tuple[dict | None, str | None]:
    for key in record.get("query_keys") or []:
        payload = query_lookup.get(key) or {}
        nodes = payload.get("taxonomy_nodes") or []
        if not nodes:
            continue
        taxonomy = nodes[0].get("taxonomy") or {}
        if str(taxonomy.get("tax_id") or "").strip():
            return taxonomy, key
    return None, None


def build_ncbi_taxonomy(taxonomy: dict, lineage_lookup: dict[str, dict], matched_query: str) -> dict:
    lineage_ids = [str(item) for item in (taxonomy.get("lineage") or []) if str(item).strip()]
    lineage_path: list[dict] = []
    classification: dict[str, dict] = {}
    for tax_id in lineage_ids:
        node = lineage_lookup.get(tax_id) or {}
        if not node:
            continue
        rank = str(node.get("rank") or "").upper()
        item = {
            "taxid": tax_id,
            "name": node.get("organism_name") or "",
            "rank": str(node.get("rank") or "").lower(),
        }
        lineage_path.append(item)
        rank_key = RANK_KEY_MAP.get(rank)
        if rank_key and item["name"]:
            classification[rank_key] = item

    terminal = {
        "taxid": str(taxonomy.get("tax_id") or "").strip(),
        "name": taxonomy.get("organism_name") or "",
        "rank": str(taxonomy.get("rank") or "").lower(),
    }
    rank_key = RANK_KEY_MAP.get(str(taxonomy.get("rank") or "").upper())
    if rank_key and terminal["name"]:
        classification[rank_key] = terminal

    return {
        "matched_query": matched_query,
        "scientific_name": terminal["name"],
        "rank": terminal["rank"],
        "blast_name": taxonomy.get("blast_name") or "",
        "common_name": taxonomy.get("common_name") or taxonomy.get("genbank_common_name") or "",
        "lineage_taxids": lineage_ids,
        "lineage_path": lineage_path,
        "classification": classification,
        "terminal": terminal,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--query-cache-dir", required=True)
    parser.add_argument("--lineage-cache-dir", required=True)
    parser.add_argument("--mapping-out", required=True)
    parser.add_argument("--pending-lineage-out", required=True)
    args = parser.parse_args()

    manifest = load_json(Path(args.manifest))
    query_cache_dir = Path(args.query_cache_dir)
    lineage_cache_dir = Path(args.lineage_cache_dir)
    query_lookup: dict[str, dict] = {}
    for path in sorted(query_cache_dir.glob("*.json")):
        try:
            query_lookup[path.stem] = load_json(path)
        except (OSError, json.JSONDecodeError):
            continue
    lineage_lookup = load_taxonomy_cache(lineage_cache_dir)

    entries: dict[str, dict] = {}
    pending_ids: set[str] = set()
    unresolved: list[str] = []
    linked_count = 0

    for record in manifest.get("records") or []:
        taxonomy, query_key = select_taxonomy(record, query_lookup)
        if not taxonomy or not query_key:
            unresolved.append(record.get("id") or record.get("species") or "unknown")
            continue
        linked_count += 1
        all_ids = [str(taxonomy.get("tax_id") or "").strip(), *[str(item) for item in (taxonomy.get("lineage") or []) if str(item).strip()]]
        for tax_id in all_ids:
            if tax_id and tax_id not in lineage_lookup:
                pending_ids.add(tax_id)
        matched_query = (query_lookup.get(query_key) or {}).get("taxonomy_nodes", [{}])[0].get("query", [""])[0]
        entries[record["id"]] = {
            "taxid": str(taxonomy.get("tax_id") or "").strip(),
            "ncbi_taxonomy": build_ncbi_taxonomy(taxonomy, lineage_lookup, matched_query),
        }

    Path(args.pending_lineage_out).write_text(
        json.dumps(sorted(pending_ids, key=lambda item: int(item) if item.isdigit() else item), ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    mapping = {
        "schema": "ncbi_taxonomy_links/v1",
        "source": {
            "label": "NCBI Datasets Taxonomy API",
            "url": "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/",
        },
        "entries": entries,
        "summary": {
            "linked_count": linked_count,
            "unresolved_count": len(unresolved),
        },
        "unresolved": unresolved,
    }
    Path(args.mapping_out).write_text(json.dumps(mapping, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"linked_count": linked_count, "unresolved_count": len(unresolved), "pending_lineage_count": len(pending_ids)}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
