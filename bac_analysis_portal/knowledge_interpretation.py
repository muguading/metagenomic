from __future__ import annotations

import re
from functools import lru_cache
from pathlib import Path

from .knowledge_base import load_knowledge_base_bundle
from .runtime_paths import _resolve_runtime_database_root

SALMONELLA_SEROVAR_ALIAS_MAP = {
    "伤寒沙门菌": ["Typhi", "S.Typhi", "Salmonella Typhi", "伤寒沙门氏菌"],
    "甲型副伤寒沙门菌": ["Paratyphi A", "S.Paratyphi A", "Salmonella Paratyphi A", "甲型副伤寒沙门氏菌"],
    "乙型副伤寒沙门菌": ["Paratyphi B", "S.Paratyphi B", "Salmonella Paratyphi B", "乙型副伤寒沙门氏菌"],
    "丙型副伤寒沙门菌": ["Paratyphi C", "S.Paratyphi C", "Salmonella Paratyphi C", "丙型副伤寒沙门氏菌"],
    "肠炎沙门菌": ["Enteritidis", "S.Enteritidis", "Salmonella Enteritidis", "肠炎沙门氏菌"],
    "鼠伤寒沙门菌": ["Typhimurium", "S.Typhimurium", "Salmonella Typhimurium", "鼠伤寒沙门氏菌"],
    "猪霍乱沙门菌": ["Choleraesuis", "Choleracsuis", "S.Choleraesuis", "S.Choleracsuis", "Salmonella Choleraesuis", "猪霍乱沙门氏菌"],
    "德尔卑沙门菌": ["Derby", "S.Derby", "Salmonella Derby", "德尔卑沙门氏菌"],
    "伦敦沙门菌": ["London", "S.London", "Salmonella London", "伦敦沙门氏菌"],
    "斯坦利沙门菌": ["Stanley", "S.Stanley", "Salmonella Stanley", "斯坦利沙门氏菌"],
    "山夫登堡沙门菌": ["Senftenberg", "S.Senftenberg", "Salmonella Senftenberg", "山夫登堡沙门氏菌"],
    "阿贡纳沙门菌": ["Agona", "S.Agona", "Salmonella Agona", "阿贡纳沙门氏菌"],
    "汤卜逊沙门菌": ["Thompson", "S.Thompson", "Salmonella Thompson", "汤卜逊沙门氏菌"],
    "罗森沙门菌": ["Rissen", "S.Rissen", "Salmonella Rissen", "罗森沙门氏菌"],
    "I,4,[5],12:i:-": ["1,4,[5],12:i:-", "S.1,4,[5],12:i:-", "单相鼠伤寒沙门菌", "单相鼠伤寒沙门氏菌"],
}


def _normalize_taxonomy_lookup_name(value: object) -> str:
    text = re.sub(r"\s+", " ", str(value or "").strip()).lower()
    if not text:
        return ""
    text = text.replace("，", " ").replace(",", " ")
    text = re.sub(r"\s+", " ", text).strip()
    return text

def _taxonomy_link_payload(pathogen: dict) -> dict | None:
    taxid = str(pathogen.get("taxid") or "").strip()
    ncbi = pathogen.get("ncbi_taxonomy") if isinstance(pathogen.get("ncbi_taxonomy"), dict) else {}
    if not taxid and not ncbi:
        return None
    classification = ncbi.get("classification") if isinstance(ncbi.get("classification"), dict) else {}
    terminal = ncbi.get("terminal") if isinstance(ncbi.get("terminal"), dict) else {}
    return {
        "pathogen_id": str(pathogen.get("id") or "").strip(),
        "species": str(pathogen.get("species") or "").strip(),
        "taxid": taxid,
        "scientific_name": str(ncbi.get("scientific_name") or terminal.get("name") or pathogen.get("species") or "").strip(),
        "rank": str(ncbi.get("rank") or terminal.get("rank") or "").strip(),
        "order": str(((classification.get("order") or {}) if isinstance(classification.get("order"), dict) else {}).get("name") or "").strip(),
        "family": str(((classification.get("family") or {}) if isinstance(classification.get("family"), dict) else {}).get("name") or "").strip(),
        "genus": str(((classification.get("genus") or {}) if isinstance(classification.get("genus"), dict) else {}).get("name") or "").strip(),
        "species_rank": str(((classification.get("species") or {}) if isinstance(classification.get("species"), dict) else {}).get("name") or "").strip(),
    }

def _normalize_serotype_lookup_text(value: object) -> str:
    text = str(value or "").strip().lower()
    if not text:
        return ""
    text = text.replace("（", "(").replace("）", ")")
    text = re.sub(r"\s+", "", text)
    return text

def _expand_salmonella_serovar_aliases(value: object) -> list[str]:
    raw = str(value or "").strip()
    if not raw:
        return []
    aliases = [raw]
    normalized_raw = _normalize_serotype_lookup_text(raw)
    for canonical, values in SALMONELLA_SEROVAR_ALIAS_MAP.items():
        candidates = [canonical, *values]
        normalized_candidates = [_normalize_serotype_lookup_text(item) for item in candidates]
        if normalized_raw and normalized_raw in normalized_candidates:
            aliases.extend(candidates)
            break
    seen: set[str] = set()
    output: list[str] = []
    for alias in aliases:
        text = str(alias or "").strip()
        key = _normalize_serotype_lookup_text(text)
        if text and key and key not in seen:
            seen.add(key)
            output.append(text)
    return output

def _tokenize_serotype_values(value: object) -> list[str]:
    text = str(value or "").strip()
    if not text:
        return []
    text = text.replace("，", "、").replace(",", "、").replace("；", "、").replace(";", "、")
    tokens = [item.strip() for item in re.split(r"[、/]", text) if item.strip()]
    normalized: list[str] = []
    for token in tokens:
        token_normalized = _normalize_serotype_lookup_text(token)
        if token_normalized and token_normalized not in normalized:
            normalized.append(token_normalized)
    return normalized

def _normalize_mlst_lookup_text(value: object) -> str:
    text = str(value or "").strip().upper()
    if not text or text in {"-", "--", "NONE", "FALSE"}:
        return ""
    text = re.sub(r"\s+", "", text)
    if re.fullmatch(r"\d+", text):
        return f"ST{text}"
    match = re.fullmatch(r"ST[-_]?(\d+)", text)
    if match:
        return f"ST{match.group(1)}"
    return text

def _build_kb_mlst_index(project_root_text: str) -> dict[str, list[dict]]:
    bundle = load_knowledge_base_bundle(project_root_text)
    collections = bundle.get("collections") if isinstance(bundle, dict) else {}
    pathogens = collections.get("pathogens") if isinstance(collections, dict) else []
    index: dict[str, list[dict]] = {}
    for pathogen in pathogens if isinstance(pathogens, list) else []:
        if not isinstance(pathogen, dict):
            continue
        mlst_rows = pathogen.get("mlst_associations")
        if not isinstance(mlst_rows, list) or not mlst_rows:
            continue
        record = {
            "pathogen_id": str(pathogen.get("id") or "").strip(),
            "species": str(pathogen.get("species") or "").strip(),
            "common_name": str(pathogen.get("common_name") or "").strip(),
            "mlst_associations": [item for item in mlst_rows if isinstance(item, dict)],
            "lineage_associations": [item for item in (pathogen.get("lineage_associations") or []) if isinstance(item, dict)],
        }
        candidate_names: list[str] = []
        for value in [
            pathogen.get("species"),
            *([item for item in (pathogen.get("aliases") or []) if isinstance(item, str)]),
        ]:
            normalized = _normalize_taxonomy_lookup_name(value)
            if normalized and normalized not in candidate_names:
                candidate_names.append(normalized)
        for name in candidate_names:
            index.setdefault(name, []).append(record)
    return index

def _lookup_kb_mlst_association(project_root_text: str, species_name: object, st_text: object) -> dict | None:
    normalized_species = _normalize_taxonomy_lookup_name(species_name)
    normalized_st = _normalize_mlst_lookup_text(st_text)
    if not normalized_species or not normalized_st or normalized_species in {"false", "-", "none"}:
        return None
    candidates = _build_kb_mlst_index(project_root_text).get(normalized_species, [])
    if not candidates:
        return None
    for candidate in candidates:
        for association in candidate.get("mlst_associations") or []:
            if _normalize_mlst_lookup_text(association.get("mlst")) != normalized_st:
                continue
            lineage_matches = []
            suffix = normalized_st[2:] if normalized_st.startswith("ST") else normalized_st
            for lineage in candidate.get("lineage_associations") or []:
                lineage_name = str(lineage.get("lineage") or "").strip()
                lineage_norm = re.sub(r"[^A-Z0-9]+", "", lineage_name.upper())
                if suffix and suffix in lineage_norm:
                    lineage_matches.append(lineage)
            return {
                "pathogen_id": candidate.get("pathogen_id"),
                "species": candidate.get("species"),
                "common_name": candidate.get("common_name"),
                "association": association,
                "lineage_matches": lineage_matches,
            }
    return None

def _build_mlst_knowledge_summary(project_root_text: str, rows: list[list[str]]) -> dict:
    summary_items: list[dict] = []
    seen: set[tuple[str, str]] = set()
    for row in rows:
        if not isinstance(row, list) or len(row) < 11:
            continue
        st_value = row[9] if len(row) > 9 else ""
        species_value = row[10] if len(row) > 10 else ""
        key = (_normalize_taxonomy_lookup_name(species_value), _normalize_mlst_lookup_text(st_value))
        if not key[0] or not key[1] or key in seen:
            continue
        seen.add(key)
        match = _lookup_kb_mlst_association(project_root_text, species_value, st_value)
        association = match.get("association") if isinstance(match, dict) else {}
        lineage_matches = match.get("lineage_matches") if isinstance(match, dict) else []
        if not association:
            continue
        lineage_text = "；".join(
            [str(item.get("lineage") or "").strip() for item in lineage_matches if str(item.get("lineage") or "").strip()]
        ) or "-"
        summary_items.append({
            "species": str((match or {}).get("species") or species_value or "-"),
            "st": str(association.get("mlst") or _normalize_mlst_lookup_text(st_value) or "-"),
            "lineage_type": str(association.get("lineage_type") or "-"),
            "lineage_text": lineage_text,
            "virulence": association.get("virulence_associations") or [],
            "resistance": association.get("resistance_associations") or [],
            "regional": association.get("regional_distribution") or [],
            "interpretation": str(association.get("interpretation") or "").strip(),
        })
    if not summary_items:
        return {"headline": "", "items": []}
    primary = summary_items[0]
    headline = f"当前 MLST 结果提示样本更接近 {primary['species']} 的 {primary['st']} 分型背景。"
    if primary.get("lineage_type") and primary["lineage_type"] != "-":
        headline += f" 该分型在知识库中归为“{primary['lineage_type']}”。"
    return {"headline": headline, "items": summary_items[:3]}

def _build_kb_serotype_index(project_root_text: str) -> dict[str, list[dict]]:
    bundle = load_knowledge_base_bundle(project_root_text)
    collections = bundle.get("collections") if isinstance(bundle, dict) else {}
    pathogens = collections.get("pathogens") if isinstance(collections, dict) else []
    index: dict[str, list[dict]] = {}
    for pathogen in pathogens if isinstance(pathogens, list) else []:
        if not isinstance(pathogen, dict):
            continue
        serotype_rows = pathogen.get("serotype_associations")
        if not isinstance(serotype_rows, list) or not serotype_rows:
            continue
        serotype_panels = pathogen.get("serotype_panels")
        record = {
            "pathogen_id": str(pathogen.get("id") or "").strip(),
            "species": str(pathogen.get("species") or "").strip(),
            "common_name": str(pathogen.get("common_name") or "").strip(),
            "serotype_panels": [item for item in serotype_panels if isinstance(item, dict)] if isinstance(serotype_panels, list) else [],
            "serotype_associations": [item for item in serotype_rows if isinstance(item, dict)],
        }
        candidate_names: list[str] = []
        for value in [
            pathogen.get("species"),
            *([item for item in (pathogen.get("aliases") or []) if isinstance(item, str)]),
        ]:
            normalized = _normalize_taxonomy_lookup_name(value)
            if normalized and normalized not in candidate_names:
                candidate_names.append(normalized)
        for name in candidate_names:
            index.setdefault(name, []).append(record)
    return index

def _lookup_kb_serotype_association(project_root_text: str, species_name: object, serotype_text: object) -> dict | None:
    normalized_species = _normalize_taxonomy_lookup_name(species_name)
    normalized_serotype = _normalize_serotype_lookup_text(serotype_text)
    if not normalized_species or not normalized_serotype or normalized_species in {"false", "-", "none"}:
        return None
    serotype_index = _build_kb_serotype_index(project_root_text)
    candidates = serotype_index.get(normalized_species, [])
    if not candidates:
        fuzzy_candidates: list[dict] = []
        for indexed_name, rows in serotype_index.items():
            if not indexed_name:
                continue
            if indexed_name == normalized_species or indexed_name in normalized_species or normalized_species in indexed_name:
                for row in rows:
                    if row not in fuzzy_candidates:
                        fuzzy_candidates.append(row)
        candidates = fuzzy_candidates
    if not candidates:
        return None
    best_match: dict | None = None
    best_score = -1
    for candidate in candidates:
        for association in candidate.get("serotype_associations") or []:
            candidate_values = [
                _normalize_serotype_lookup_text(association.get("serotype")),
                *[
                    _normalize_serotype_lookup_text(value)
                    for value in (association.get("serotype_aliases") or [])
                    if isinstance(value, str)
                ],
            ]
            if str(candidate.get("pathogen_id") or "") == "salmonella_enterica":
                candidate_values.extend(
                    _normalize_serotype_lookup_text(value)
                    for value in _expand_salmonella_serovar_aliases(association.get("serotype"))
                )
            candidate_values = [value for value in candidate_values if value]
            if not candidate_values:
                continue
            score = 0
            if normalized_serotype in candidate_values:
                score = 100
            elif len(normalized_serotype) >= 3 and any(
                len(value) >= 3 and (value in normalized_serotype or normalized_serotype in value)
                for value in candidate_values
            ):
                score = 60
            if score > best_score:
                best_score = score
                best_match = {
                    "pathogen_id": candidate.get("pathogen_id"),
                    "species": candidate.get("species"),
                    "common_name": candidate.get("common_name"),
                    "panels": candidate.get("serotype_panels") or [],
                    "association": association,
                }
    return best_match if best_score > 0 else None

def _build_kb_hepatitis_typing_index(project_root_text: str) -> dict[str, dict[str, list[dict]]]:
    bundle = load_knowledge_base_bundle(project_root_text)
    collections = bundle.get("collections") if isinstance(bundle, dict) else {}
    typing_rules = collections.get("typing_rules") if isinstance(collections, dict) else []
    by_broad: dict[str, list[dict]] = {}
    by_broad_serotype: dict[str, list[dict]] = {}
    for rule in typing_rules if isinstance(typing_rules, list) else []:
        if not isinstance(rule, dict):
            continue
        broad_type = str(rule.get("broad_type") or "").strip().upper()
        if broad_type not in {"HAV", "HBV", "HCV", "HDV", "HEV"}:
            continue
        by_broad.setdefault(broad_type, []).append(rule)
        candidate_values = [
            rule.get("serotype"),
            rule.get("subtype"),
            *[value for value in (rule.get("aliases") or []) if isinstance(value, str)],
        ]
        for value in candidate_values:
            normalized = _normalize_serotype_lookup_text(value)
            if normalized:
                by_broad_serotype.setdefault(f"{broad_type}:{normalized}", []).append(rule)
    return {"by_broad": by_broad, "by_broad_serotype": by_broad_serotype}

def _lookup_kb_hepatitis_typing_rule(project_root_text: str, broad_type: object, subtype_text: object = "") -> dict | None:
    normalized_broad = str(broad_type or "").strip().upper()
    if normalized_broad not in {"HAV", "HBV", "HCV", "HDV", "HEV"}:
        return None
    index = _build_kb_hepatitis_typing_index(project_root_text)
    normalized_subtype = _normalize_serotype_lookup_text(subtype_text)
    if normalized_subtype:
        candidates = index.get("by_broad_serotype", {}).get(f"{normalized_broad}:{normalized_subtype}", [])
        subtype_matches = [item for item in candidates if str(item.get("level") or "").strip() == "subtype"]
        if subtype_matches:
            return subtype_matches[0]
        if candidates:
            return candidates[0]
    broad_matches = [
        item
        for item in index.get("by_broad", {}).get(normalized_broad, [])
        if str(item.get("level") or "").strip() == "broad_type"
    ]
    return broad_matches[0] if broad_matches else None

def _hepatitis_typing_rule_to_summary_item(rule: dict, matched_on: object) -> dict:
    references = rule.get("references") if isinstance(rule.get("references"), list) else []
    reference_accessions = rule.get("reference_accessions") if isinstance(rule.get("reference_accessions"), list) else []
    return {
        "species": str(rule.get("species") or "").strip(),
        "common_name": str(rule.get("common_name") or "").strip(),
        "panel": str(rule.get("panel") or "-"),
        "serotype": str(rule.get("serotype") or rule.get("subtype") or rule.get("broad_type") or "").strip(),
        "matched_on": str(matched_on or rule.get("serotype") or "").strip(),
        "level": str(rule.get("level") or "").strip(),
        "broad_type": str(rule.get("broad_type") or "").strip(),
        "subtype": str(rule.get("subtype") or "").strip(),
        "reference_count": str(rule.get("reference_count") or len(reference_accessions) or len(references) or "0"),
        "reference_accessions": [str(item).strip() for item in reference_accessions if str(item).strip()],
        "references": [item for item in references if isinstance(item, dict)][:8],
        "virulence": [],
        "resistance": [],
        "regional": [],
        "interpretation": str(rule.get("interpretation") or "").strip(),
    }

def _build_kb_hiv_typing_index(project_root_text: str) -> dict[str, dict[str, list[dict]]]:
    bundle = load_knowledge_base_bundle(project_root_text)
    collections = bundle.get("collections") if isinstance(bundle, dict) else {}
    typing_rules = collections.get("typing_rules") if isinstance(collections, dict) else []
    by_broad: dict[str, list[dict]] = {}
    by_broad_serotype: dict[str, list[dict]] = {}
    for rule in typing_rules if isinstance(typing_rules, list) else []:
        if not isinstance(rule, dict):
            continue
        broad_type = str(rule.get("broad_type") or "").strip().upper()
        if broad_type not in {"HIV-1", "HIV-2"}:
            continue
        by_broad.setdefault(broad_type, []).append(rule)
        candidate_values = [
            rule.get("serotype"),
            rule.get("subtype"),
            *[value for value in (rule.get("aliases") or []) if isinstance(value, str)],
        ]
        for value in candidate_values:
            normalized = _normalize_serotype_lookup_text(value)
            if normalized:
                by_broad_serotype.setdefault(f"{broad_type}:{normalized}", []).append(rule)
    return {"by_broad": by_broad, "by_broad_serotype": by_broad_serotype}

def _lookup_kb_hiv_typing_rule(project_root_text: str, broad_type: object, subtype_text: object = "") -> dict | None:
    normalized_broad = str(broad_type or "").strip().upper()
    if normalized_broad not in {"HIV-1", "HIV-2"}:
        return None
    index = _build_kb_hiv_typing_index(project_root_text)
    normalized_subtype = _normalize_serotype_lookup_text(subtype_text)
    if normalized_subtype:
        candidates = index.get("by_broad_serotype", {}).get(f"{normalized_broad}:{normalized_subtype}", [])
        subtype_matches = [item for item in candidates if str(item.get("level") or "").strip() == "subtype"]
        if subtype_matches:
            return subtype_matches[0]
        if candidates:
            return candidates[0]
    broad_matches = [
        item
        for item in index.get("by_broad", {}).get(normalized_broad, [])
        if str(item.get("level") or "").strip() == "broad_type"
    ]
    return broad_matches[0] if broad_matches else None

def _hiv_typing_rule_to_summary_item(rule: dict, matched_on: object) -> dict:
    references = rule.get("references") if isinstance(rule.get("references"), list) else []
    reference_accessions = rule.get("reference_accessions") if isinstance(rule.get("reference_accessions"), list) else []
    regional = rule.get("regional_distribution") if isinstance(rule.get("regional_distribution"), list) else []
    return {
        "species": str(rule.get("species") or "").strip(),
        "common_name": str(rule.get("common_name") or "").strip(),
        "panel": str(rule.get("panel") or "-"),
        "serotype": str(rule.get("serotype") or rule.get("subtype") or rule.get("broad_type") or "").strip(),
        "matched_on": str(matched_on or rule.get("serotype") or "").strip(),
        "level": str(rule.get("level") or "").strip(),
        "broad_type": str(rule.get("broad_type") or "").strip(),
        "subtype": str(rule.get("subtype") or "").strip(),
        "reference_count": str(rule.get("reference_count") or len(reference_accessions) or len(references) or "0"),
        "reference_accessions": [str(item).strip() for item in reference_accessions if str(item).strip()],
        "references": [item for item in references if isinstance(item, dict)][:8],
        "virulence": [],
        "resistance": [],
        "regional": [str(item).strip() for item in regional if str(item).strip()],
        "interpretation": str(rule.get("interpretation") or "").strip(),
    }

def _build_kb_tb_typing_index(project_root_text: str) -> dict[str, list[dict]]:
    bundle = load_knowledge_base_bundle(project_root_text)
    collections = bundle.get("collections") if isinstance(bundle, dict) else {}
    typing_rules = collections.get("typing_rules") if isinstance(collections, dict) else []
    index: dict[str, list[dict]] = {}
    for rule in typing_rules if isinstance(typing_rules, list) else []:
        if not isinstance(rule, dict):
            continue
        if str(rule.get("pathogen_id") or "").strip() != "mycobacterium_tuberculosis":
            continue
        candidate_values = [
            rule.get("serotype"),
            rule.get("subtype"),
            *[value for value in (rule.get("aliases") or []) if isinstance(value, str)],
        ]
        if str(rule.get("level") or "").strip() == "broad_type":
            candidate_values.append(rule.get("broad_type"))
        for value in candidate_values:
            normalized = _normalize_serotype_lookup_text(value)
            if not normalized:
                continue
            bucket = index.setdefault(normalized, [])
            if rule not in bucket:
                bucket.append(rule)
    return index

def _tb_typing_rule_to_summary_item(rule: dict, matched_on: object) -> dict:
    regional = rule.get("regional_distribution") if isinstance(rule.get("regional_distribution"), list) else []
    markers = rule.get("key_markers") if isinstance(rule.get("key_markers"), list) else []
    return {
        "species": str(rule.get("species") or "").strip(),
        "common_name": str(rule.get("common_name") or "").strip(),
        "panel": str(rule.get("panel") or "-"),
        "serotype": str(rule.get("serotype") or rule.get("subtype") or rule.get("broad_type") or "").strip(),
        "matched_on": str(matched_on or rule.get("serotype") or rule.get("subtype") or "").strip(),
        "level": str(rule.get("level") or "").strip(),
        "broad_type": str(rule.get("broad_type") or "").strip(),
        "subtype": str(rule.get("subtype") or "").strip(),
        "reference_count": str(rule.get("reference_count") or "0"),
        "reference_accessions": [str(item).strip() for item in (rule.get("reference_accessions") or []) if str(item).strip()],
        "references": [item for item in (rule.get("references") or []) if isinstance(item, dict)][:8],
        "virulence": [],
        "resistance": [],
        "regional": [str(item).strip() for item in regional if str(item).strip()],
        "key_markers": [str(item).strip() for item in markers if str(item).strip()],
        "interpretation": str(rule.get("interpretation") or "").strip(),
    }

def _tb_typing_specificity(rule: dict) -> tuple[int, int]:
    level = str(rule.get("level") or "").strip().lower()
    level_rank = {
        "sublineage": 4,
        "lineage": 3,
        "subtype": 2,
        "broad_type": 1,
    }.get(level, 0)
    label = str(rule.get("subtype") or rule.get("serotype") or "").strip()
    dot_count = label.count(".")
    return (level_rank, dot_count)

def _build_tb_knowledge_summary(project_root_text: str, species_name: object, typing_values: list[object]) -> dict:
    normalized_species = _normalize_taxonomy_lookup_name(species_name)
    if normalized_species not in {
        "mycobacterium tuberculosis",
        "mtuberculosis",
        "mtb",
        "tb",
    }:
        return {"headline": "", "items": []}
    index = _build_kb_tb_typing_index(project_root_text)
    ordered_matches: list[tuple[str, dict]] = []
    seen_ids: set[str] = set()
    for value in typing_values:
        normalized = _normalize_serotype_lookup_text(value)
        if not normalized:
            continue
        for rule in index.get(normalized, []):
            rule_id = str(rule.get("id") or "").strip()
            if not rule_id or rule_id in seen_ids:
                continue
            seen_ids.add(rule_id)
            ordered_matches.append((str(value or "").strip(), rule))
    if not ordered_matches:
        broad_rule = next((
            rule for rule in index.get("mtbc", [])
            if isinstance(rule, dict) and str(rule.get("level") or "").strip() == "broad_type"
        ), None)
        if broad_rule:
            ordered_matches.append(("MTBC", broad_rule))
    if not ordered_matches:
        return {"headline": "", "items": []}

    best_match = max(ordered_matches, key=lambda item: _tb_typing_specificity(item[1]))
    best_rule = best_match[1]
    best_specificity = _tb_typing_specificity(best_rule)
    retained_matches = [
        item for item in ordered_matches
        if _tb_typing_specificity(item[1]) == best_specificity
    ]
    items = [_tb_typing_rule_to_summary_item(rule, matched_on) for matched_on, rule in retained_matches[:2]]
    if not items:
        return {"headline": "", "items": []}
    headline = "当前结核家系结果已命中知识库，可结合家系背景、地区分布和关键标记做流行病学解释。"
    top_labels = [item.get("serotype") for item in items if str(item.get("serotype") or "").strip()]
    if top_labels:
        headline = f"当前结核家系结果已命中知识库中的 {' / '.join(top_labels)} 条目，可直接补充分子流调背景解释。"
    return {"headline": headline, "items": items}

def _build_kb_pathogen_profile_index(project_root_text: str) -> dict[str, dict]:
    bundle = load_knowledge_base_bundle(project_root_text)
    collections = bundle.get("collections") if isinstance(bundle, dict) else {}
    pathogens = collections.get("pathogens") if isinstance(collections, dict) else []
    index: dict[str, dict] = {}
    for pathogen in pathogens if isinstance(pathogens, list) else []:
        if not isinstance(pathogen, dict):
            continue
        pathogen_id = str(pathogen.get("id") or "").strip()
        if pathogen_id:
            index[pathogen_id] = pathogen
    return index

def _lookup_kb_pathogen_profile(project_root_text: str, pathogen_id: str) -> dict | None:
    if not project_root_text or not pathogen_id:
        return None
    return _build_kb_pathogen_profile_index(project_root_text).get(str(pathogen_id).strip())

def _first_table_value(table: dict, candidates: list[str]) -> str:
    columns = table.get("columns") if isinstance(table, dict) else []
    rows = table.get("rows") if isinstance(table, dict) else []
    if not isinstance(columns, list) or not isinstance(rows, list) or not rows:
        return ""
    first_row = rows[0] if isinstance(rows[0], list) else []
    if not isinstance(first_row, list):
        return ""
    for column in candidates:
        if column not in columns:
            continue
        index = columns.index(column)
        if index >= len(first_row):
            continue
        value = str(first_row[index] or "").strip()
        if value and value not in {"-", "--", "False", "false", "none", "None"}:
            return value
    return ""

def _infer_priority_serotype_species(species_value: object, serotype_result: dict | None, checkm_info: dict | None = None) -> str:
    species_text = str(species_value or "").strip()
    normalized = _normalize_taxonomy_lookup_name(species_text)
    if normalized and normalized not in {"false", "-", "none", "--"}:
        return species_text
    checkm_candidates = [
        (checkm_info or {}).get("species_name") if isinstance(checkm_info, dict) else "",
        (checkm_info or {}).get("mlst_species_name") if isinstance(checkm_info, dict) else "",
    ]
    for value in checkm_candidates:
        text = str(value or "").strip()
        normalized_text = _normalize_taxonomy_lookup_name(text)
        if "salmonella" in normalized_text:
            return "Salmonella enterica"
        if normalized_text and normalized_text not in {"false", "-", "none", "--"}:
            return text
    serotype_table = {
        "columns": (serotype_result or {}).get("columns", []) if isinstance(serotype_result, dict) else [],
        "rows": (serotype_result or {}).get("rows", []) if isinstance(serotype_result, dict) else [],
    }
    serotype_clues = [
        _first_table_value(serotype_table, ["亚型全称", "菌种", "血清型注释信息(simple)", "血清型注释信息(details)", "血清型"]),
    ]
    for clue in serotype_clues:
        normalized_clue = _normalize_serotype_lookup_text(clue)
        if "salmonella" in normalized_clue or "沙门" in normalized_clue:
            return "Salmonella enterica"
        for alias_values in SALMONELLA_SEROVAR_ALIAS_MAP.values():
            if normalized_clue in {_normalize_serotype_lookup_text(item) for item in alias_values}:
                return "Salmonella enterica"
    return species_text

def _priority_serotype_candidates(serotype_value: object, serotype_result: dict | None) -> list[str]:
    candidates: list[str] = []

    def _push(value: object) -> None:
        text = str(value or "").strip()
        if not text or text in {"-", "--", "False", "false", "none", "None"}:
            return
        for alias in _expand_salmonella_serovar_aliases(text):
            key = _normalize_serotype_lookup_text(alias)
            if key and key not in {_normalize_serotype_lookup_text(item) for item in candidates}:
                candidates.append(alias)
        key = _normalize_serotype_lookup_text(text)
        if key and key not in {_normalize_serotype_lookup_text(item) for item in candidates}:
            candidates.append(text)
        for token in _tokenize_serotype_values(text):
            if token and token not in {_normalize_serotype_lookup_text(item) for item in candidates}:
                candidates.append(token)

    _push(serotype_value)
    columns = (serotype_result or {}).get("columns", []) if isinstance(serotype_result, dict) else []
    rows = (serotype_result or {}).get("rows", []) if isinstance(serotype_result, dict) else []
    first_row = rows[0] if isinstance(rows, list) and rows and isinstance(rows[0], list) else []
    for column in ["亚型全称", "菌种", "血清型注释信息(simple)", "血清型注释信息(details)", "血清型"]:
        if column in columns:
            index = columns.index(column)
            if index < len(first_row):
                _push(first_row[index])
    return candidates

def _enrich_priority_serotype_rows(project_root_text: str, table: dict, serotype_result: dict | None = None, checkm_info: dict | None = None) -> dict:
    columns = table.get("columns") if isinstance(table, dict) else []
    rows = table.get("rows") if isinstance(table, dict) else []
    if not isinstance(columns, list) or not isinstance(rows, list) or not columns:
        return {"columns": columns or [], "rows": rows or []}
    species_index = columns.index("物种") if "物种" in columns else -1
    serotype_index = columns.index("血清型") if "血清型" in columns else -1
    if species_index < 0 or serotype_index < 0:
        return {"columns": columns, "rows": rows}
    extra_columns = [
        "知识库血清型面板",
        "知识库命中血清型",
        "血清型-毒力关联",
        "血清型-耐药关联",
        "血清型-地域分布",
        "血清型知识库提示",
    ]
    enriched_rows: list[list[str]] = []
    for row in rows:
        current = list(row) if isinstance(row, list) else [str(row.get(column, "") or "") for column in columns]
        species_value = current[species_index] if species_index < len(current) else ""
        serotype_value = current[serotype_index] if serotype_index < len(current) else ""
        resolved_species = _infer_priority_serotype_species(species_value, serotype_result, checkm_info)
        if (not str(serotype_value or "").strip() or str(serotype_value or "").strip() in {"-", "--"}) and serotype_result:
            fallback_serotype = _first_table_value(serotype_result, ["亚型全称", "血清型注释信息(simple)", "血清型", "菌种"])
            if fallback_serotype:
                serotype_value = fallback_serotype
                if serotype_index < len(current):
                    current[serotype_index] = fallback_serotype
        if resolved_species and _normalize_taxonomy_lookup_name(species_value) in {"false", "-", "none", "--"} and species_index < len(current):
            current[species_index] = resolved_species
        match = None
        for candidate_serotype in _priority_serotype_candidates(serotype_value, serotype_result):
            match = _lookup_kb_serotype_association(project_root_text, resolved_species or species_value, candidate_serotype)
            if match:
                break
        association = match.get("association") if isinstance(match, dict) else {}
        panels = match.get("panels") if isinstance(match, dict) else []
        panel_text = "-"
        if isinstance(association, dict) and association.get("panel"):
            panel_text = str(association.get("panel") or "-")
        elif isinstance(panels, list) and panels:
            panel_values = []
            for item in panels:
                if not isinstance(item, dict):
                    continue
                title = str(item.get("panel") or "").strip()
                members = item.get("serotypes") if isinstance(item.get("serotypes"), list) else []
                if title:
                    panel_values.append(title)
                elif members:
                    panel_values.append(" / ".join(str(member).strip() for member in members if str(member).strip()))
            panel_text = "；".join([value for value in panel_values if value]) or "-"
        enriched_rows.append(current + [
            panel_text,
            str(association.get("serotype") or "-"),
            "；".join(association.get("virulence_associations") or []) if isinstance(association.get("virulence_associations"), list) and association.get("virulence_associations") else "-",
            "；".join(association.get("resistance_associations") or []) if isinstance(association.get("resistance_associations"), list) and association.get("resistance_associations") else "-",
            "；".join(association.get("regional_distribution") or []) if isinstance(association.get("regional_distribution"), list) and association.get("regional_distribution") else "-",
            str(association.get("interpretation") or "-"),
        ])
    return {"columns": columns + extra_columns, "rows": enriched_rows}

def _extract_viral_serotype_candidates(serotype_result: dict) -> list[str]:
    candidates: list[str] = []

    def _push(value: object) -> None:
        text = str(value or "").strip()
        if not text or text in {"-", "--", "none", "None"}:
            return
        if text not in candidates:
            candidates.append(text)
        for token in _tokenize_serotype_values(text):
            if token and token not in [_normalize_serotype_lookup_text(item) for item in candidates]:
                candidates.append(token)

    if not isinstance(serotype_result, dict):
        return candidates

    mode = str(serotype_result.get("mode") or "").strip()
    if mode == "hadv_typing":
        # For HAdV, the knowledge-summary primary target should be the final
        # total typing conclusion rather than PHF gene-level component types.
        _push(serotype_result.get("predicted_clade"))
        _push(serotype_result.get("predicted_serotype"))
        notes = str(serotype_result.get("notes") or "").strip()
        if notes:
            for pattern in [
                r"总分型：([^；]+)",
                r"自动选择参考型别：([^；]+)",
            ]:
                for matched in re.findall(pattern, notes):
                    _push(matched)
        return candidates

    ordered_keys = [
        "predicted_serotype",
        "predicted_subtype",
        "predicted_clade",
        "predicted_group",
        "predicted_lineage",
        "influenza_type",
        "ha_subtype",
        "na_subtype",
        "predicted_outbreak",
    ]
    if mode == "bandavirus_typing":
        ordered_keys = [
            "predicted_group",
            "predicted_clade",
            "predicted_lineage",
            "predicted_serotype",
            "predicted_subtype",
            "influenza_type",
            "ha_subtype",
            "na_subtype",
            "predicted_outbreak",
        ]

    for key in ordered_keys:
        _push(serotype_result.get(key))

    summary_cards = serotype_result.get("summary_cards")
    if isinstance(summary_cards, list):
        for item in summary_cards:
            if not isinstance(item, dict):
                continue
            label = str(item.get("label") or "").strip()
            value = item.get("value")
            if not label:
                continue
            if any(
                keyword in label
                for keyword in [
                    "亚型",
                    "分型",
                    "Lineage",
                    "Genotype",
                    "Clade",
                    "病毒属",
                    "病毒种",
                    "Penton",
                    "Hexon",
                    "Fiber",
                ]
            ) and not any(
                skip in label
                for skip in [
                    "覆盖",
                    "QC",
                    "参考序列",
                    "最近参考",
                    "注释文件",
                    "平均深度",
                    "支持 reads",
                    "覆盖碱基",
                ]
            ):
                _push(value)

    notes = str(serotype_result.get("notes") or "").strip()
    if notes:
        for pattern in [
            r"自动选择参考型别：([^；]+)",
            r"双位点分型：([^；]+)",
            r"RdRp/VP1：([^；]+)",
            r"大亚型：([^；]+)",
            r"大类分型：([^；]+)",
            r"S 子亚型：([^；]+)",
            r"A_F\(LMS\)：([^；]+)",
            r"CJ\(LMS\)：([^；]+)",
            r"G/P 分型：([^；]+)",
            r"组合分型：([^；]+)",
            r"总分型：([^；]+)",
            r"ORF2 分型：([^；]+)",
            r"病毒种：([^；]+)",
        ]:
            for matched in re.findall(pattern, notes):
                _push(matched)

    return candidates

def _build_viral_serotype_knowledge_summary(project_root_text: str, species_name: object, serotype_result: dict) -> dict:
    if not project_root_text or not isinstance(serotype_result, dict):
        return {"headline": "", "items": []}
    mode = str(serotype_result.get("mode") or "").strip()
    resolved_species_name = species_name
    if mode == "hiv_resistance":
        broad_type = str(serotype_result.get("predicted_group") or "").strip().upper()
        subtype = str(serotype_result.get("predicted_subtype") or serotype_result.get("predicted_clade") or "").strip()
        items: list[dict] = []
        broad_rule = _lookup_kb_hiv_typing_rule(project_root_text, broad_type, broad_type)
        subtype_rule = (
            _lookup_kb_hiv_typing_rule(project_root_text, broad_type, subtype)
            if broad_type == "HIV-1" and subtype and subtype not in {"-", "--", broad_type}
            else None
        )
        if isinstance(broad_rule, dict):
            items.append(_hiv_typing_rule_to_summary_item(broad_rule, broad_type))
        if isinstance(subtype_rule, dict) and subtype_rule.get("id") != (broad_rule or {}).get("id"):
            items.append(_hiv_typing_rule_to_summary_item(subtype_rule, subtype))
        if not items:
            return {"headline": "", "items": []}
        species_label = str((subtype_rule or broad_rule or {}).get("species") or species_name or "HIV").strip()
        if broad_type == "HIV-1" and subtype and subtype not in {"-", "--", broad_type}:
            headline = f"当前样本知识库先命中 {broad_type} 大亚型，并进一步关联到 HIV-1 {subtype} 子亚型 / CRF 背景。"
        elif broad_type:
            headline = f"当前样本知识库命中 {species_label} 的 {broad_type} 大亚型背景。"
        else:
            headline = f"当前样本命中 {species_label} 分型知识库。"
        return {"headline": headline, "items": items[:4]}
    if mode == "hepatovirus_typing":
        broad_type = str(serotype_result.get("predicted_group") or "").strip().upper()
        subtype = str(serotype_result.get("predicted_subtype") or serotype_result.get("predicted_clade") or "").strip()
        items: list[dict] = []
        broad_rule = _lookup_kb_hepatitis_typing_rule(project_root_text, broad_type, broad_type)
        subtype_rule = _lookup_kb_hepatitis_typing_rule(project_root_text, broad_type, subtype) if subtype and subtype not in {"-", "--"} else None
        if isinstance(broad_rule, dict):
            items.append(_hepatitis_typing_rule_to_summary_item(broad_rule, broad_type))
        if isinstance(subtype_rule, dict) and subtype_rule.get("id") != (broad_rule or {}).get("id"):
            items.append(_hepatitis_typing_rule_to_summary_item(subtype_rule, subtype))
        if not items:
            return {"headline": "", "items": []}
        species_label = str((subtype_rule or broad_rule or {}).get("species") or species_name or "肝炎病毒").strip()
        if broad_type and subtype and subtype not in {"-", "--", broad_type}:
            headline = f"当前样本知识库命中 {broad_type} 大亚型，并进一步关联到 {broad_type} {subtype} 子亚型/基因型背景。"
        elif broad_type:
            headline = f"当前样本知识库命中 {species_label} 的 {broad_type} 大亚型背景。"
        else:
            headline = f"当前样本命中 {species_label} 分型知识库。"
        return {"headline": headline, "items": items[:4]}
    normalized_species = _normalize_taxonomy_lookup_name(resolved_species_name)
    if not normalized_species or normalized_species in {"false", "-", "none"}:
        return {"headline": "", "items": []}
    candidates = _extract_viral_serotype_candidates(serotype_result)
    if not candidates:
        return {"headline": "", "items": []}

    items: list[dict] = []
    seen: set[tuple[str, str]] = set()
    orthohantavirus_hfrs_types = {"HTNV", "SEOV", "DOBV", "PUUV", "AMRV"}
    orthohantavirus_hps_types = {"ANDV", "BAYV", "BCCV", "CASV", "SNVV"}
    for candidate in candidates:
        match = _lookup_kb_serotype_association(project_root_text, resolved_species_name, candidate)
        association = match.get("association") if isinstance(match, dict) else {}
        if not isinstance(association, dict) or not association:
            continue
        matched_serotype = str(association.get("serotype") or "").strip()
        pathogen_id = str((match or {}).get("pathogen_id") or "").strip()
        interpretation_text = str(association.get("interpretation") or "").strip()
        virulence_values = list(association.get("virulence_associations") or []) if isinstance(association.get("virulence_associations"), list) else []
        regional_values = list(association.get("regional_distribution") or []) if isinstance(association.get("regional_distribution"), list) else []
        if pathogen_id == "orthohantavirus":
            profile_id = ""
            if matched_serotype.upper() in orthohantavirus_hfrs_types:
                profile_id = "hantaviruses_hfrs"
            elif matched_serotype.upper() in orthohantavirus_hps_types:
                profile_id = "hantaviruses_pulmonary_syndrome"
            profile = _lookup_kb_pathogen_profile(project_root_text, profile_id) if profile_id else None
            if isinstance(profile, dict):
                clinical_text = str(profile.get("clinical_significance") or "").strip()
                public_health_text = str(profile.get("public_health_significance") or "").strip()
                note_values = profile.get("interpretation_notes") if isinstance(profile.get("interpretation_notes"), list) else []
                extra_parts = [value for value in [clinical_text, public_health_text] if value]
                if extra_parts:
                    interpretation_text = "；".join([part for part in [interpretation_text, *extra_parts] if part])
                for value in note_values:
                    text = str(value or "").strip()
                    if text and text not in virulence_values:
                        virulence_values.append(text)
        key = (
            pathogen_id,
            matched_serotype or _normalize_serotype_lookup_text(candidate),
        )
        if key in seen:
            continue
        seen.add(key)
        items.append(
            {
                "species": str((match or {}).get("species") or resolved_species_name or "-"),
                "common_name": str((match or {}).get("common_name") or "").strip(),
                "panel": str(association.get("panel") or "-"),
                "serotype": matched_serotype or str(candidate).strip(),
                "matched_on": str(candidate).strip(),
                "virulence": virulence_values,
                "resistance": association.get("resistance_associations") or [],
                "regional": regional_values,
                "interpretation": interpretation_text,
            }
        )
        if len(items) >= 4:
            break

    if not items:
        return {"headline": "", "items": []}
    primary = items[0]
    primary_serotype = str(primary.get("serotype") or primary.get("matched_on") or "").strip()
    primary_species = str(primary.get("species") or species_name or "当前病毒").strip()
    headline = ""
    if mode == "hepatovirus_typing":
        broad_type = str(serotype_result.get("predicted_group") or "").strip()
        hav_subtype = str(serotype_result.get("predicted_clade") or "").strip()
        if broad_type and broad_type not in {"-", "--"} and hav_subtype and hav_subtype not in {"-", "--"}:
            headline = f"当前样本知识库首先命中 {broad_type} 大亚型，并进一步匹配到 HAV {hav_subtype} 子亚型背景。"
        elif broad_type and broad_type not in {"-", "--"}:
            headline = f"当前样本知识库首先命中 {broad_type} 大亚型背景。"
        elif primary_serotype:
            headline = f"当前样本知识库命中 {primary_species} 的 {primary_serotype} 分型背景。"
    elif primary_serotype:
        headline = f"当前病毒分型结果提示样本更接近 {primary_species} 的 {primary_serotype} 分型背景。"
    return {"headline": headline, "items": items}

def _build_kb_taxonomy_index(project_root_text: str) -> dict[str, list[dict]]:
    bundle = load_knowledge_base_bundle(project_root_text)
    collections = bundle.get("collections") if isinstance(bundle, dict) else {}
    pathogens = collections.get("pathogens") if isinstance(collections, dict) else []
    index: dict[str, list[dict]] = {}
    for pathogen in pathogens if isinstance(pathogens, list) else []:
        if not isinstance(pathogen, dict):
            continue
        payload = _taxonomy_link_payload(pathogen)
        if not payload:
            continue
        candidate_names: list[str] = []
        for value in [
            pathogen.get("species"),
            payload.get("scientific_name"),
            payload.get("genus"),
            payload.get("species_rank"),
            *([item for item in (pathogen.get("aliases") or []) if isinstance(item, str)]),
        ]:
            normalized = _normalize_taxonomy_lookup_name(value)
            if normalized and normalized not in candidate_names:
                candidate_names.append(normalized)
        for name in candidate_names:
            index.setdefault(name, []).append(payload)
    return index

def _choose_taxonomy_match(candidates: list[dict], query_name: str, genus_name: str = "") -> dict | None:
    normalized_query = _normalize_taxonomy_lookup_name(query_name)
    normalized_genus = _normalize_taxonomy_lookup_name(genus_name)
    ranked: list[tuple[int, dict]] = []
    for candidate in candidates:
        score = 0
        scientific_name = _normalize_taxonomy_lookup_name(candidate.get("scientific_name"))
        species_name = _normalize_taxonomy_lookup_name(candidate.get("species"))
        candidate_genus = _normalize_taxonomy_lookup_name(candidate.get("genus"))
        rank = str(candidate.get("rank") or "").lower()
        if normalized_query and scientific_name == normalized_query:
            score += 120
        if normalized_query and species_name == normalized_query:
            score += 110
        if normalized_query and candidate_genus == normalized_query:
            score += 60
        if normalized_genus and candidate_genus == normalized_genus:
            score += 40
        if rank == "species":
            score += 15
        elif rank == "genus":
            score += 5
        ranked.append((score, candidate))
    ranked.sort(key=lambda item: item[0], reverse=True)
    return ranked[0][1] if ranked and ranked[0][0] > 0 else None

def _lookup_kb_taxonomy(index: dict[str, list[dict]] | None, species_name: str = "", genus_name: str = "") -> dict | None:
    if not index:
        return None
    candidate_keys: list[str] = []
    for value in [species_name, genus_name]:
        normalized = _normalize_taxonomy_lookup_name(value)
        if normalized and normalized not in candidate_keys:
            candidate_keys.append(normalized)
    for key in candidate_keys:
        candidates = index.get(key) or []
        match = _choose_taxonomy_match(candidates, species_name or key, genus_name)
        if match:
            return match
    return None

def _normalize_kb_token(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value or "").strip().lower())

def _truncate_list(items: list[object], limit: int = 6) -> list[str]:
    values: list[str] = []
    for item in items:
        text = str(item or "").strip()
        if not text or text == "-" or text in values:
            continue
        values.append(text)
        if len(values) >= limit:
            break
    return values

def _joined_text(values: list[object], fallback: str = "--") -> str:
    cleaned = [str(value or "").strip() for value in values if str(value or "").strip() and str(value or "").strip() != "-"]
    return "、".join(cleaned) if cleaned else fallback

def _kb_gene_name_matches(detected_name: object, target_name: object) -> bool:
    detected = _normalize_kb_token(detected_name)
    target = _normalize_kb_token(target_name)
    if not detected or not target:
        return False
    if detected == target:
        return True
    if len(target) >= 4 and detected.startswith(target):
        return True
    if len(detected) >= 4 and target.startswith(detected):
        return True
    return False

def _build_kb_interpretation_context(project_root_text: str) -> dict:
    bundle = load_knowledge_base_bundle(project_root_text)
    collections = bundle.get("collections") if isinstance(bundle, dict) else {}
    pathogens = collections.get("pathogens") if isinstance(collections, dict) else []
    gene_rules = collections.get("gene_rules") if isinstance(collections, dict) else []
    event_rules = collections.get("event_rules") if isinstance(collections, dict) else []

    pathogen_by_taxid: dict[str, dict] = {}
    pathogen_by_name: dict[str, dict] = {}
    for pathogen in pathogens if isinstance(pathogens, list) else []:
        if not isinstance(pathogen, dict):
            continue
        taxid = str(pathogen.get("taxid") or "").strip()
        if taxid:
            pathogen_by_taxid[taxid] = pathogen
        names = [
            pathogen.get("species"),
            pathogen.get("common_name"),
            *((pathogen.get("aliases") or []) if isinstance(pathogen.get("aliases"), list) else []),
        ]
        ncbi = pathogen.get("ncbi_taxonomy") if isinstance(pathogen.get("ncbi_taxonomy"), dict) else {}
        terminal = ncbi.get("terminal") if isinstance(ncbi.get("terminal"), dict) else {}
        classification = ncbi.get("classification") if isinstance(ncbi.get("classification"), dict) else {}
        names.extend([
            ncbi.get("scientific_name"),
            terminal.get("name"),
            ((classification.get("genus") or {}) if isinstance(classification.get("genus"), dict) else {}).get("name"),
            ((classification.get("species") or {}) if isinstance(classification.get("species"), dict) else {}).get("name"),
        ])
        for name in names:
            normalized = _normalize_taxonomy_lookup_name(name)
            if normalized and normalized not in pathogen_by_name:
                pathogen_by_name[normalized] = pathogen

    return {
        "pathogen_by_taxid": pathogen_by_taxid,
        "pathogen_by_name": pathogen_by_name,
        "gene_rules": [item for item in gene_rules if isinstance(item, dict)],
        "event_rules": [item for item in event_rules if isinstance(item, dict)],
    }

def _rank_species_rows(species_taxonomy: dict, subspecies_taxonomy: dict | None = None) -> list[dict]:
    ranked: list[dict] = []
    for dataset in (species_taxonomy or {}, subspecies_taxonomy or {}):
        for row in dataset.get("rows") or []:
            if not isinstance(row, dict):
                continue
            species_name = str(row.get("种") or row.get("亚种") or row.get("NCBI种") or "").strip()
            if not species_name or species_name == "-":
                continue
            ranked.append(row)
    ranked.sort(
        key=lambda item: (
            int(item.get("序列数量数值") or 0),
            float(item.get("比例数值") or 0.0),
        ),
        reverse=True,
    )
    return ranked

def _match_kb_pathogen_profile(context: dict, row: dict | None) -> dict | None:
    if not isinstance(row, dict):
        return None
    pathogen_by_taxid = context.get("pathogen_by_taxid") if isinstance(context.get("pathogen_by_taxid"), dict) else {}
    pathogen_by_name = context.get("pathogen_by_name") if isinstance(context.get("pathogen_by_name"), dict) else {}
    taxid = str(row.get("NCBI TaxID") or "").strip()
    if taxid and taxid != "-" and taxid in pathogen_by_taxid:
        return pathogen_by_taxid[taxid]
    for value in [row.get("NCBI学名"), row.get("NCBI种"), row.get("种"), row.get("亚种"), row.get("属")]:
        normalized = _normalize_taxonomy_lookup_name(value)
        if normalized and normalized in pathogen_by_name:
            return pathogen_by_name[normalized]
    return None

def _collect_detected_gene_hits(elements: dict, gene_type: str) -> list[dict]:
    columns = elements.get("columns") or []
    rows = elements.get("rows") or []
    if not columns or not rows:
        return []
    gene_columns = ["基因名称", "毒力基因", "耐药基因", "Gene", "Best_Hit_ARO"]
    category_columns = ["VF分类", "耐药药物", "Drug Class", "Drug", "Resistance Mechanism"]
    hits: list[dict] = []
    for row in rows:
        if not isinstance(row, list):
            continue
        gene_name = ""
        for column in gene_columns:
            if column in columns:
                index = columns.index(column)
                gene_name = str(row[index] if index < len(row) else "").strip()
                if gene_name and gene_name != "-":
                    break
        if not gene_name or gene_name == "-":
            continue
        category = ""
        for column in category_columns:
            if column in columns:
                index = columns.index(column)
                category = str(row[index] if index < len(row) else "").strip()
                if category and category != "-":
                    break
        hits.append({
            "gene_type": gene_type,
            "gene_name": gene_name,
            "category": category or "未分类",
        })
    return hits

def _pathogen_name_pool(pathogen: dict | None) -> set[str]:
    if not isinstance(pathogen, dict):
        return set()
    names = {
        _normalize_taxonomy_lookup_name(pathogen.get("species")),
        _normalize_taxonomy_lookup_name(pathogen.get("common_name")),
    }
    for item in pathogen.get("aliases") or []:
        names.add(_normalize_taxonomy_lookup_name(item))
    return {item for item in names if item}

def _rule_applies_to_pathogen(rule: dict, pathogen: dict | None) -> bool:
    if not isinstance(rule, dict):
        return False
    species_scope = rule.get("species_scope")
    scope = rule.get("scope") if isinstance(rule.get("scope"), dict) else {}
    scope_species = scope.get("species")
    candidates = species_scope if isinstance(species_scope, list) else scope_species if isinstance(scope_species, list) else []
    if not candidates:
        return True
    pathogen_names = _pathogen_name_pool(pathogen)
    if not pathogen_names:
        return False
    for item in candidates:
        if _normalize_taxonomy_lookup_name(item) in pathogen_names:
            return True
    return False

def _evaluate_gene_rule_hits(context: dict, pathogen: dict | None, resistance_hits: list[dict], virulence_hits: list[dict]) -> list[dict]:
    hits_by_type = {"ARG": resistance_hits, "VF": virulence_hits}
    matched: list[dict] = []
    for rule in context.get("gene_rules") or []:
        gene_type = str(rule.get("gene_type") or "").upper()
        if gene_type not in hits_by_type or not _rule_applies_to_pathogen(rule, pathogen):
            continue
        detected_names = []
        for hit in hits_by_type[gene_type]:
            if _kb_gene_name_matches(hit.get("gene_name"), rule.get("gene_name")):
                detected_names.append(str(hit.get("gene_name") or "").strip())
        if not detected_names:
            continue
        matched.append({
            "id": str(rule.get("id") or "").strip(),
            "gene_type": gene_type,
            "gene_name": str(rule.get("gene_name") or "").strip(),
            "risk_level": str(rule.get("risk_level") or "").strip().lower() or "medium",
            "report_label": str(rule.get("report_label") or rule.get("gene_name") or "").strip(),
            "functional_meaning": str(rule.get("functional_meaning") or "").strip(),
            "evidence_strength": str(rule.get("evidence_strength") or "").strip().lower() or "medium",
            "matched_hits": _truncate_list(detected_names, 6),
        })
    matched.sort(key=lambda item: (item.get("risk_level") == "high", item.get("evidence_strength") == "high", len(item.get("matched_hits") or [])), reverse=True)
    return matched

def _evaluate_single_event_condition(condition: dict, resistance_hits: list[dict], virulence_hits: list[dict], has_mge_signal: bool) -> tuple[bool, list[str]]:
    if not isinstance(condition, dict):
        return False, []
    if condition.get("mobile_element_related") is True:
        return has_mge_signal, (["移动元件相关信号"] if has_mge_signal else [])
    gene_type = str(condition.get("gene_type") or "").upper()
    gene_names = condition.get("gene_name_in")
    if gene_type not in {"ARG", "VF"} or not isinstance(gene_names, list):
        return False, []
    source_hits = resistance_hits if gene_type == "ARG" else virulence_hits
    matched_hits: list[str] = []
    for target in gene_names:
        for hit in source_hits:
            if _kb_gene_name_matches(hit.get("gene_name"), target):
                matched_hits.append(str(hit.get("gene_name") or "").strip())
    return bool(matched_hits), _truncate_list(matched_hits, 8)

def _evaluate_event_rule_hits(context: dict, pathogen: dict | None, resistance_hits: list[dict], virulence_hits: list[dict], mge_monitoring: dict) -> list[dict]:
    has_mge_signal = int(((mge_monitoring.get("overview") or {}).get("total_hits")) or 0) > 0
    matched: list[dict] = []
    for rule in context.get("event_rules") or []:
        if not _rule_applies_to_pathogen(rule, pathogen):
            continue
        trigger = rule.get("trigger") if isinstance(rule.get("trigger"), dict) else {}
        all_conditions = trigger.get("all") if isinstance(trigger.get("all"), list) else []
        any_conditions = trigger.get("any") if isinstance(trigger.get("any"), list) else []
        matched_hits: list[str] = []
        success = True
        if all_conditions:
            for condition in all_conditions:
                condition_ok, condition_hits = _evaluate_single_event_condition(condition, resistance_hits, virulence_hits, has_mge_signal)
                if not condition_ok:
                    success = False
                    break
                matched_hits.extend(condition_hits)
        if success and any_conditions:
            any_success = False
            for condition in any_conditions:
                condition_ok, condition_hits = _evaluate_single_event_condition(condition, resistance_hits, virulence_hits, has_mge_signal)
                if condition_ok:
                    any_success = True
                    matched_hits.extend(condition_hits)
            success = any_success
        if not success:
            continue
        interpretation = rule.get("interpretation") if isinstance(rule.get("interpretation"), dict) else {}
        matched.append({
            "id": str(rule.get("id") or "").strip(),
            "title": str(rule.get("title") or "").strip(),
            "priority": str(rule.get("priority") or "").strip(),
            "summary": str(interpretation.get("summary") or "").strip(),
            "risk_level": str(interpretation.get("risk_level") or "").strip().lower() or "medium",
            "confidence": str(interpretation.get("confidence") or "").strip().lower() or "medium",
            "recommendation": [str(item).strip() for item in (rule.get("recommendation") or []) if str(item).strip()],
            "matched_hits": _truncate_list(matched_hits, 10),
        })
    matched.sort(key=lambda item: (item.get("risk_level") == "high", item.get("priority") == "P0", item.get("confidence") == "high"), reverse=True)
    return matched

def _summarize_typing_information(mlst_result: dict, serotype_result: dict) -> list[str]:
    fragments: list[str] = []
    rows = mlst_result.get("rows") or []
    columns = mlst_result.get("columns") or []
    if rows and columns:
        first_row = rows[0] if isinstance(rows[0], list) else []
        for key in ("ST型", "ST"):
            if key in columns:
                index = columns.index(key)
                value = str(first_row[index] if index < len(first_row) else "").strip()
                if value and value != "-":
                    fragments.append(f"MLST：{value}")
                    break
    if isinstance(serotype_result, dict):
        mode = str(serotype_result.get("mode") or "").strip()
        if mode == "sars_cov_2_nextclade":
            clade = str(serotype_result.get("predicted_clade") or "").strip()
            pango = str(serotype_result.get("pango_lineage") or "").strip()
            if clade and clade != "-":
                fragments.append(f"Nextclade：{clade}")
            if pango and pango != "-":
                fragments.append(f"Pango：{pango}")
            return fragments
        if mode == "monkeypox_nextclade":
            clade = str(serotype_result.get("predicted_clade") or "").strip()
            lineage = str(serotype_result.get("predicted_lineage") or "").strip()
            outbreak = str(serotype_result.get("predicted_outbreak") or "").strip()
            if clade and clade != "-":
                fragments.append(f"Nextclade：{clade}")
            if lineage and lineage != "-":
                fragments.append(f"Lineage：{lineage}")
            if outbreak and outbreak != "-":
                fragments.append(f"Outbreak：{outbreak}")
            return fragments
        predicted = str(serotype_result.get("predicted_serotype") or "").strip()
        if predicted and predicted != "-":
            fragments.append(f"血清型：{predicted}")
        else:
            rows = serotype_result.get("rows") or []
            columns = serotype_result.get("columns") or []
            first_row = rows[0] if rows and isinstance(rows[0], list) else []
            for key in ("血清型", "Serotype"):
                if key in columns:
                    index = columns.index(key)
                    value = str(first_row[index] if index < len(first_row) else "").strip()
                    if value and value != "-":
                        fragments.append(f"血清型：{value}")
                        break
    return fragments

def _interpretation_risk_label(pathogen: dict | None, gene_hits: list[dict], event_hits: list[dict], is_meta_method: bool) -> str:
    high_events = sum(1 for item in event_hits if str(item.get("risk_level") or "").lower() == "high")
    high_genes = sum(1 for item in gene_hits if str(item.get("risk_level") or "").lower() == "high")
    medium_genes = sum(1 for item in gene_hits if str(item.get("risk_level") or "").lower() == "medium")
    if high_events or high_genes >= 2:
        return "高风险"
    if high_genes or medium_genes >= 2:
        return "中风险"
    if isinstance(pathogen, dict) and pathogen.get("notifiable") and not is_meta_method:
        return "高风险"
    return "低风险"

def _derive_clinical_conclusion(dominant_species_name: str, pathogen: dict | None, event_hits: list[dict], is_meta_method: bool, supporting_species: list[dict]) -> str:
    profile_name = str((pathogen or {}).get("common_name") or (pathogen or {}).get("species") or dominant_species_name or "当前优势物种").strip()
    pathogen_type = str((pathogen or {}).get("pathogen_type") or "").strip().lower()
    is_viral_pathogen = pathogen_type == "病毒"
    if event_hits:
        summary = str(event_hits[0].get("summary") or "").strip()
        if is_meta_method:
            return f"当前样本以 {profile_name} 为主要优势病原线索，且出现“{summary}”相关证据，建议结合样本类型、临床表现与其他病原背景综合判断。"
        if is_viral_pathogen:
            return f"当前病毒结果以 {profile_name} 为主导检出对象，且出现“{summary}”相关证据，提示其在当前样本中具有明确分子判读意义。"
        return f"当前单菌结果以 {profile_name} 为主导检出病原，且出现“{summary}”相关证据，提示其临床意义较强。"
    if supporting_species and is_meta_method:
        return f"当前宏基因组结果以 {profile_name} 为主要优势病原线索，但同时存在其他高丰度物种背景，建议结合临床与流行病学资料谨慎解释。"
    if is_viral_pathogen:
        return (
            f"当前宏基因组结果以 {profile_name} 为主要病毒线索，宜结合覆盖度、分型结果和临床证据综合判读。"
            if is_meta_method
            else f"当前病毒结果以 {profile_name} 为主导判读对象，知识库提示其具有明确的病毒学和临床相关性。"
        )
    return (
        f"当前宏基因组结果以 {profile_name} 为主要病原线索，宜结合样本背景和临床证据综合判读。"
        if is_meta_method
        else f"当前单菌结果以 {profile_name} 为主导物种，知识库提示其具有明确临床相关性。"
    )

def _species_identity_key(row: dict | None) -> str:
    if not isinstance(row, dict):
        return ""
    for value in [row.get("NCBI学名"), row.get("NCBI种"), row.get("种")]:
        normalized = _normalize_taxonomy_lookup_name(value)
        if normalized and normalized != "-":
            return normalized
    return ""

def _build_knowledge_interpretation(
    project_root: Path,
    species_taxonomy: dict,
    subspecies_taxonomy: dict,
    resistance_elements: dict,
    virulence_elements: dict,
    mge_monitoring: dict,
    mlst_result: dict,
    serotype_result: dict,
    public_health_support: dict,
    is_meta_method: bool,
) -> dict:
    context = _build_kb_interpretation_context(str(_resolve_runtime_database_root()))
    ranked_rows = _rank_species_rows(species_taxonomy, subspecies_taxonomy)
    if not ranked_rows:
        return {"status": "empty"}

    dominant_row = ranked_rows[0]
    dominant_species_name = str(dominant_row.get("种") or dominant_row.get("亚种") or dominant_row.get("NCBI种") or "--").strip() or "--"
    pathogen = _match_kb_pathogen_profile(context, dominant_row)
    resistance_hits = _collect_detected_gene_hits(resistance_elements, "ARG")
    virulence_hits = _collect_detected_gene_hits(virulence_elements, "VF")
    matched_gene_rules = _evaluate_gene_rule_hits(context, pathogen, resistance_hits, virulence_hits)
    matched_event_rules = _evaluate_event_rule_hits(context, pathogen, resistance_hits, virulence_hits, mge_monitoring)
    dominant_identity = _species_identity_key(dominant_row)
    supporting_species = []
    for row in (species_taxonomy.get("rows") or [])[1:6]:
        ratio = float(row.get("比例数值") or 0.0)
        reads = int(row.get("序列数量数值") or 0)
        species_name = str(row.get("种") or row.get("亚种") or "").strip()
        if not species_name or species_name == "-" or ratio < 5:
            continue
        if dominant_identity and _species_identity_key(row) == dominant_identity:
            continue
        supporting_species.append({"species": species_name, "ratio": round(ratio, 2), "reads": reads})
        if len(supporting_species) >= 3:
            break
    intraspecies_signals = []
    for row in subspecies_taxonomy.get("rows") or []:
        ratio = float(row.get("比例数值") or 0.0)
        reads = int(row.get("序列数量数值") or 0)
        subspecies_name = str(row.get("亚种") or row.get("种") or "").strip()
        if not subspecies_name or subspecies_name == "-" or ratio < 1:
            continue
        intraspecies_signals.append({
            "label": subspecies_name,
            "ratio": round(ratio, 2),
            "reads": reads,
        })
        if len(intraspecies_signals) >= 4:
            break

    resistance_genes = _truncate_list([item.get("gene_name") for item in resistance_hits], 8)
    virulence_genes = _truncate_list([item.get("gene_name") for item in virulence_hits], 8)
    resistance_classes = _truncate_list([item.get("category") for item in resistance_hits], 8)
    infection_sites = [str(item).strip() for item in ((pathogen or {}).get("typical_infection_sites") or []) if str(item).strip()]
    syndrome_roles = [
        f"{str(item.get('syndrome') or '').strip()}（{str(item.get('role') or '相关病原').strip()}）"
        for item in ((pathogen or {}).get("syndrome_associations") or [])
        if isinstance(item, dict) and str(item.get("syndrome") or "").strip()
    ]
    risk_level = _interpretation_risk_label(pathogen, matched_gene_rules, matched_event_rules, is_meta_method)
    typing_fragments = _summarize_typing_information(mlst_result, serotype_result)
    top_event_summaries = _truncate_list([item.get("summary") for item in matched_event_rules], 3)
    top_rule_labels = _truncate_list([item.get("report_label") for item in matched_gene_rules], 6)
    profile_name = str((pathogen or {}).get("common_name") or (pathogen or {}).get("species") or dominant_species_name or "--").strip()
    pathogen_type = str((pathogen or {}).get("pathogen_type") or "").strip().lower()
    is_viral_pathogen = pathogen_type == "病毒"
    clinical_significance = str((pathogen or {}).get("clinical_significance") or "").strip()
    public_health_significance = str((pathogen or {}).get("public_health_significance") or "").strip()
    significance = (
        f"当前宏基因组结果中，{profile_name} 为主要优势病原线索。{clinical_significance}"
        if is_meta_method and clinical_significance
        else clinical_significance or f"当前结果以 {profile_name} 为主要判读对象。"
    )
    infection_hint_values = infection_sites or syndrome_roles
    infection_hint = _joined_text(infection_hint_values, "结合送检部位、临床综合征与宿主背景判断")
    treatment_hint = (
        f"当前耐药相关证据提示需重点关注 {_joined_text(resistance_classes, '相关耐药类别')}，建议结合药敏试验及感染部位综合制定治疗方案。"
        if resistance_classes
        else ("当前未见明确高风险耐药类别，但仍建议结合药敏试验、感染部位与宿主因素综合决策。" if not is_viral_pathogen else "当前结果以病毒分型、覆盖度和变异证据为主，药物相关解释仍需结合具体抗病毒方案、宿主状态及临床诊疗背景综合判断。")
    )
    virulence_hint = (
        f"检出 {len(virulence_hits)} 条毒力相关记录，重点提示包括 {_joined_text(top_event_summaries or top_rule_labels, '相关毒力因素')}。"
        if virulence_hits or matched_event_rules
        else ("当前未见明确高风险毒力组合，但仍应结合样本来源和临床综合征审慎判读。" if not is_viral_pathogen else "当前未见额外的病毒毒力规则命中，结果更适合结合分型、覆盖度和关键位点背景进行解释。")
    )
    clinical_recommendations = _truncate_list(
        [
            *((matched_event_rules[0].get("recommendation") or []) if matched_event_rules else []),
            *(
                [
                    "建议结合病毒分型、覆盖度、关键变异位点及样本来源综合确认结果解释边界。",
                    "建议结合患者症状、肝功能指标、流行病学接触史及影像学证据综合判断。",
                ]
                if is_viral_pathogen
                else [
                    "建议结合培养、药敏试验和样本来源综合确认病原学意义。",
                    "对于宏基因组结果，应结合背景菌和宿主状态避免把高丰度线索直接等同于致病菌结论。" if is_meta_method else "建议结合患者症状、炎症指标及影像学证据综合判断。",
                ]
            ),
        ],
        4,
    )
    cdc_recommendations = _truncate_list(
        [
            *((matched_event_rules[0].get("recommendation") or []) if matched_event_rules else []),
            "建议结合病例时空分布、接触史与同源性结果持续评估传播风险。",
            "如涉及重点病原或法定报告病原，建议按疾控流程及时上报和复核。",
        ],
        4,
    )
    support_evidence = []
    if clinical_significance:
        support_evidence.append({"manual": "病原体画像", "rule_type": "pathogen_profile", "basis": clinical_significance})
    if public_health_significance:
        support_evidence.append({"manual": "公共卫生画像", "rule_type": "public_health_profile", "basis": public_health_significance})
    for item in matched_event_rules[:2]:
        support_evidence.append({"manual": item.get("title") or "组合事件规则", "rule_type": "event_rule", "basis": item.get("summary") or ""})

    transmission_risk = "高" if (pathogen and pathogen.get("notifiable")) or risk_level == "高风险" else "中" if risk_level == "中风险" or matched_event_rules else "低"
    typing_summary = "；".join(typing_fragments) if typing_fragments else ("未获得稳定病毒分型结果" if is_viral_pathogen else "未获得稳定分型/血清型结果")
    reportable_label = (
        str(public_health_support.get("reporting_hint") or "").strip()
        or ("法定报告或重点报告病原体" if (pathogen and pathogen.get("notifiable")) else "重点监测病原体" if matched_event_rules or public_health_significance else "一般监测病原体")
    )

    return {
        "status": "ready",
        "dominant_species": {
            "species": dominant_species_name,
            "ratio": round(float(dominant_row.get("比例数值") or 0.0), 2),
            "reads": int(dominant_row.get("序列数量数值") or 0),
            "taxid": str(dominant_row.get("NCBI TaxID") or "-").strip() or "-",
            "scientific_name": str(dominant_row.get("NCBI学名") or "-").strip() or "-",
        },
        "pathogen_profile": {
            "id": str((pathogen or {}).get("id") or "").strip(),
            "species": str((pathogen or {}).get("species") or dominant_species_name).strip(),
            "common_name": str((pathogen or {}).get("common_name") or profile_name).strip(),
            "notifiable": bool((pathogen or {}).get("notifiable")),
            "pathogen_type": str((pathogen or {}).get("pathogen_type") or "").strip(),
        },
        "supporting_species": supporting_species,
        "intraspecies_signals": intraspecies_signals,
        "matched_gene_rules": matched_gene_rules,
        "matched_event_rules": matched_event_rules,
        "clinical": {
            "status": "ready",
            "speciesName": profile_name,
            "significance": significance,
            "commonLabel": (
                "属于常见临床病毒病原"
                if is_viral_pathogen and clinical_significance
                else ("属于常见临床致病菌" if clinical_significance else "已纳入病原知识库判读范围")
            ),
            "infectionHint": infection_hint,
            "resistanceCount": len(resistance_hits),
            "virulenceCount": len(virulence_hits),
            "resistanceGenes": resistance_genes,
            "resistanceClasses": resistance_classes,
            "virulenceGenes": virulence_genes,
            "riskLevel": risk_level,
            "treatmentHint": treatment_hint,
            "virulenceHint": virulence_hint,
            "conclusion": _derive_clinical_conclusion(dominant_species_name, pathogen, matched_event_rules, is_meta_method, supporting_species),
            "evidence": _truncate_list([
                f"主导物种：{dominant_species_name}（占比 {round(float(dominant_row.get('比例数值') or 0.0), 2):.2f}%）",
                f"知识库匹配：{profile_name}",
                *top_event_summaries,
                *[f"{item['report_label']}：{_joined_text(item.get('matched_hits') or [], item.get('gene_name') or '--')}" for item in matched_gene_rules[:3]],
            ], 6),
            "recommendations": clinical_recommendations,
        },
        "cdc": {
            "status": "ready",
            "speciesName": profile_name,
            "genus": str(dominant_row.get("NCBI属") or dominant_row.get("属") or "-").strip() or "-",
            "family": str(dominant_row.get("NCBI科") or dominant_row.get("科") or "-").strip() or "-",
            "reportableLabel": reportable_label,
            "pathogenNote": public_health_significance or clinical_significance or f"{profile_name} 已纳入知识库监测判读范围。",
            "resistanceCount": len(resistance_hits),
            "virulenceCount": len(virulence_hits),
            "resistanceGenes": resistance_genes,
            "virulenceGenes": virulence_genes,
            "resistanceClasses": resistance_classes,
            "hasEsbl": any("ctxm" in _normalize_kb_token(item) or "esbl" in _normalize_kb_token(item) for item in resistance_genes + resistance_classes),
            "hasCarbapenemase": any(token in _normalize_kb_token(item) for item in resistance_genes for token in ("kpc", "ndm", "oxa48", "oxa23", "vim", "imp")),
            "outbreakPotential": top_event_summaries[0] if top_event_summaries else ("提示存在较高传播处置价值，建议结合病例时空分布和接触史综合评估。" if pathogen and pathogen.get("notifiable") else "当前未形成明确暴发性分子证据，建议结合监测背景持续观察。"),
            "historyRelation": (
                f"当前样本已获得 {'，'.join(typing_fragments)}，建议结合本地历史株库与同源性结果进一步判断传播链。"
                if typing_fragments
                else ("当前尚缺乏稳定病毒分型结果，建议结合本地历史毒株数据库、流调线索和关键位点证据进一步比较。" if is_viral_pathogen else "当前尚缺乏稳定分型结果，建议结合本地历史菌株数据库和流调线索进一步比较。")
            ),
            "transmissionRisk": transmission_risk,
            "typingSummary": typing_summary,
            "monitoringAdvice": _joined_text(cdc_recommendations, "建议结合监测背景持续评估。"),
            "supportEvidence": support_evidence,
            "keyResistanceGenes": [
                {
                    "label": item.get("report_label") or item.get("gene_name") or "--",
                    "meaning": item.get("functional_meaning") or "--",
                    "status": "detected",
                    "matched_hits": item.get("matched_hits") or [],
                }
                for item in matched_gene_rules if item.get("gene_type") == "ARG"
            ][:6],
            "keyVirulenceGenes": [
                {
                    "label": item.get("report_label") or item.get("gene_name") or "--",
                    "meaning": item.get("functional_meaning") or "--",
                    "status": "detected",
                    "matched_hits": item.get("matched_hits") or [],
                }
                for item in matched_gene_rules if item.get("gene_type") == "VF"
            ][:6],
        },
    }
