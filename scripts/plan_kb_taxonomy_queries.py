from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path


def is_taxonomy_like_alias(value: str) -> bool:
    text = re.sub(r"\s+", " ", str(value or "").strip())
    if not text:
        return False
    if len(text) <= 4:
        return False
    if text.isupper():
        return False
    if " " not in text:
        return bool(re.fullmatch(r"[A-Z][a-z][A-Za-z-]{3,}", text))
    first = text.split(" ", 1)[0]
    if re.fullmatch(r"[A-Z]\.", first):
        return False
    return bool(re.fullmatch(r"[A-Z][a-zA-Z.-]+(?:\s+[A-Za-z0-9().-]+)+", text))


def build_candidates(species: str, aliases: list[str]) -> list[str]:
    raw_candidates: list[str] = [species, *[alias for alias in aliases if is_taxonomy_like_alias(alias)]]
    normalized: list[str] = []

    def push(value: str) -> None:
        text = re.sub(r"\s+", " ", str(value or "").strip())
        if not text or text in normalized:
            return
        normalized.append(text)

    for value in raw_candidates:
        push(value)
        push(re.sub(r"\bspp\.?$", "", value, flags=re.IGNORECASE).strip())
        push(re.sub(r"\bspp\.?$", "", value, flags=re.IGNORECASE).strip().split(" ")[0] if re.search(r"\bspp\.?$", value, flags=re.IGNORECASE) else "")
        push(re.sub(r"^Other\s+(.+?)\s+species$", r"\1", value, flags=re.IGNORECASE).strip())
        push(re.sub(r"^Other\s+(.+?)\s+Species$", r"\1", value, flags=re.IGNORECASE).strip())
        push(re.sub(r"^Other\s+pathogenic\s+(.+?)\s+spp\.?$", r"\1", value, flags=re.IGNORECASE).strip())
        push(re.sub(r"^Pathogenic\s+(.+)$", r"\1", value, flags=re.IGNORECASE).strip())
        push(re.sub(r"^Other\s+pathogenic\s+(.+)$", r"\1", value, flags=re.IGNORECASE).strip())
        push(re.sub(r"\s+complex$", "", value, flags=re.IGNORECASE).strip())
        push(re.sub(r"^(.+?)\s+serovar\s+.+$", r"\1", value, flags=re.IGNORECASE).strip())
        push(re.sub(r"^(.+?)\s+subsp\.\s+.+$", r"\1", value, flags=re.IGNORECASE).strip())
        parts = value.split()
        if len(parts) >= 2:
            push(" ".join(parts[:2]))

    return normalized


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", required=True)
    parser.add_argument("--manifest-out", required=True)
    args = parser.parse_args()

    project_root = Path(args.project_root)
    pathogen_dir = project_root / "database" / "knowledge_base" / "pathogens"
    records: list[dict] = []
    query_records: dict[str, dict] = {}

    for path in sorted(pathogen_dir.glob("*.json"), key=lambda item: item.name.lower()):
        payload = json.loads(path.read_text(encoding="utf-8"))
        species = str(payload.get("species") or "").strip()
        aliases = [str(item).strip() for item in payload.get("aliases") or [] if str(item).strip()]
        candidates = build_candidates(species, aliases)
        query_keys: list[str] = []
        for query in candidates:
            key = hashlib.sha1(query.encode("utf-8")).hexdigest()[:16]
            query_keys.append(key)
            query_records.setdefault(
                key,
                {
                    "key": key,
                    "query": query,
                },
            )
        records.append(
            {
                "id": payload.get("id") or path.stem,
                "file": str(path.relative_to(project_root)),
                "species": species,
                "aliases": aliases,
                "query_keys": query_keys,
            }
        )

    manifest = {
        "schema": "kb_ncbi_taxonomy_query_manifest/v1",
        "records": records,
        "queries": sorted(query_records.values(), key=lambda item: item["query"].lower()),
    }
    Path(args.manifest_out).write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"record_count": len(records), "query_count": len(query_records)}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
