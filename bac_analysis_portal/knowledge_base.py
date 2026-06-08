from __future__ import annotations

import json
from pathlib import Path


KNOWLEDGE_BASE_DIRNAME = "knowledge_base"


def _knowledge_base_root(project_root_text: str) -> Path:
    candidate = Path(project_root_text).expanduser()
    direct_root = candidate / KNOWLEDGE_BASE_DIRNAME
    nested_root = candidate / "database" / KNOWLEDGE_BASE_DIRNAME

    if candidate.name == "database":
        return direct_root
    if direct_root.is_dir() or (direct_root / "manifest.json").is_file():
        return direct_root
    return nested_root


def _read_json_file(path: Path) -> dict | list | None:
    if not path.is_file():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def _read_directory_objects(directory: Path) -> list[dict]:
    if not directory.is_dir():
        return []
    records: list[dict] = []
    for path in sorted(directory.glob("*.json"), key=lambda item: item.name.lower()):
        payload = _read_json_file(path)
        if isinstance(payload, dict):
            if not payload.get("pathogen_type"):
                payload["pathogen_type"] = "细菌"
            records.append(payload)
    return records


def _read_collection_entries(directory: Path) -> list[dict]:
    if not directory.is_dir():
        return []
    records: list[dict] = []
    for path in sorted(directory.glob("*.json"), key=lambda item: item.name.lower()):
        payload = _read_json_file(path)
        if not isinstance(payload, dict):
            continue
        entries = payload.get("entries")
        if isinstance(entries, list):
            records.extend([item for item in entries if isinstance(item, dict)])
    return records


def _read_vfdb_gene_links(root: Path) -> dict[str, list[dict]]:
    payload = _read_json_file(root / "mappings" / "vfdb_gene_links.json")
    if not isinstance(payload, dict):
        return {}
    entries = payload.get("entries")
    if not isinstance(entries, dict):
        return {}
    links: dict[str, list[dict]] = {}
    for rule_id, records in entries.items():
        if not isinstance(rule_id, str) or not isinstance(records, list):
            continue
        links[rule_id] = [record for record in records if isinstance(record, dict)]
    return links


def _read_card_gene_links(root: Path) -> dict[str, list[dict]]:
    payload = _read_json_file(root / "mappings" / "card_gene_links.json")
    if not isinstance(payload, dict):
        return {}
    entries = payload.get("entries")
    if not isinstance(entries, dict):
        return {}
    links: dict[str, list[dict]] = {}
    for rule_id, records in entries.items():
        if not isinstance(rule_id, str) or not isinstance(records, list):
            continue
        links[rule_id] = [record for record in records if isinstance(record, dict)]
    return links


def _read_ncbi_taxonomy_links(root: Path) -> dict[str, dict]:
    payload = (
        _read_json_file(root / "mappings" / "ncbi_taxonomy_links.generated.json")
        or _read_json_file(root / "mappings" / "ncbi_taxonomy_links.json")
    )
    if not isinstance(payload, dict):
        return {}
    entries = payload.get("entries")
    if not isinstance(entries, dict):
        return {}
    return {
        str(record_id): record
        for record_id, record in entries.items()
        if isinstance(record_id, str) and isinstance(record, dict)
    }


def _merge_gene_rule_links(gene_rules: list[dict], vfdb_links: dict[str, list[dict]], card_links: dict[str, list[dict]]) -> list[dict]:
    merged: list[dict] = []
    for rule in gene_rules:
        record = dict(rule)
        gene_type = str(record.get("gene_type") or "").upper()
        if gene_type == "VF":
            links = vfdb_links.get(str(record.get("id") or "").strip(), [])
            if links:
                record["vfdb_mappings"] = links
        elif gene_type == "ARG":
            links = card_links.get(str(record.get("id") or "").strip(), [])
            if links:
                record["card_mappings"] = links
        merged.append(record)
    return merged


def _merge_pathogen_taxonomy(pathogens: list[dict], taxonomy_links: dict[str, dict]) -> list[dict]:
    merged: list[dict] = []
    for item in pathogens:
        record = dict(item)
        taxonomy_record = taxonomy_links.get(str(record.get("id") or "").strip()) or {}
        if taxonomy_record.get("taxid"):
            record["taxid"] = str(taxonomy_record.get("taxid") or "").strip()
        if isinstance(taxonomy_record.get("ncbi_taxonomy"), dict):
            record["ncbi_taxonomy"] = taxonomy_record["ncbi_taxonomy"]
        merged.append(record)
    return merged


def _summarize_missing_keys(records: list[dict], required_keys: list[str]) -> list[str]:
    missing: list[str] = []
    for record in records:
        record_id = str(record.get("id") or record.get("species") or record.get("title") or "unknown").strip()
        absent = [key for key in required_keys if key not in record]
        if absent:
            missing.append(f"{record_id}: {', '.join(absent)}")
    return missing


def _summarize_missing_gene_links(records: list[dict]) -> list[str]:
    missing: list[str] = []
    for record in records:
        gene_type = str(record.get("gene_type") or "").upper()
        record_id = str(record.get("id") or record.get("gene_name") or "unknown").strip()
        if gene_type == "VF" and not record.get("vfdb_mappings"):
            missing.append(f"{record_id}: missing vfdb_mappings")
        elif gene_type == "ARG" and not record.get("card_mappings"):
            missing.append(f"{record_id}: missing card_mappings")
    return missing


def load_knowledge_base_bundle(project_root_text: str) -> dict:
    root = _knowledge_base_root(project_root_text)
    manifest = _read_json_file(root / "manifest.json")
    pathogens = _read_directory_objects(root / "pathogens")
    ncbi_taxonomy_links = _read_ncbi_taxonomy_links(root)
    pathogens = _merge_pathogen_taxonomy(pathogens, ncbi_taxonomy_links)
    gene_rules = _read_collection_entries(root / "genes")
    typing_rules = _read_collection_entries(root / "typing")
    vfdb_links = _read_vfdb_gene_links(root)
    card_links = _read_card_gene_links(root)
    gene_rules = _merge_gene_rule_links(gene_rules, vfdb_links, card_links)
    event_rules = _read_collection_entries(root / "rules" / "event_rules")
    downgrade_rules = _read_collection_entries(root / "rules" / "downgrade_rules")

    validation = {
        "pathogens": _summarize_missing_keys(
            pathogens,
            ["id", "species", "common_name", "clinical_significance", "public_health_significance"],
        ),
        "gene_rules": _summarize_missing_keys(
            gene_rules,
            ["id", "gene_name", "gene_type", "risk_level", "report_label"],
        ),
        "typing_rules": _summarize_missing_keys(
            typing_rules,
            ["id", "level", "broad_type", "serotype", "interpretation"],
        ),
        "event_rules": _summarize_missing_keys(
            event_rules,
            ["id", "title", "trigger", "interpretation"],
        ),
        "downgrade_rules": _summarize_missing_keys(
            downgrade_rules,
            ["id", "title", "condition", "penalty"],
        ),
        "gene_rule_links": _summarize_missing_gene_links(gene_rules),
    }
    validation["warnings"] = sum(len(items) for items in validation.values())

    return {
        "status": "ready" if root.is_dir() else "empty",
        "root": str(root),
        "manifest": manifest if isinstance(manifest, dict) else {},
        "collections": {
            "pathogens": pathogens,
            "gene_rules": gene_rules,
            "typing_rules": typing_rules,
            "event_rules": event_rules,
            "downgrade_rules": downgrade_rules,
        },
        "summary": {
            "pathogen_count": len(pathogens),
            "gene_rule_count": len(gene_rules),
            "typing_rule_count": len(typing_rules),
            "event_rule_count": len(event_rules),
            "downgrade_rule_count": len(downgrade_rules),
            "taxonomy_linked_pathogen_count": sum(1 for item in pathogens if str(item.get("taxid") or "").strip()),
            "vfdb_linked_gene_rule_count": sum(1 for item in gene_rules if isinstance(item.get("vfdb_mappings"), list) and item.get("vfdb_mappings")),
            "card_linked_gene_rule_count": sum(1 for item in gene_rules if isinstance(item.get("card_mappings"), list) and item.get("card_mappings")),
        },
        "validation": validation,
    }


def load_knowledge_base_summary(project_root_text: str) -> dict:
    bundle = load_knowledge_base_bundle(project_root_text)
    return {
        "status": bundle.get("status", "empty"),
        "root": bundle.get("root", ""),
        "manifest": bundle.get("manifest", {}),
        "summary": bundle.get("summary", {}),
        "validation": bundle.get("validation", {}),
    }
