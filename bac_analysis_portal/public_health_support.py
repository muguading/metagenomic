from __future__ import annotations

import json
from pathlib import Path

from .runtime_paths import _resolve_runtime_database_root

def _load_public_health_support_file(project_root_text: str, filename: str) -> dict:
    candidate = Path(project_root_text).expanduser()
    support_path = candidate / filename if candidate.name == "database" else candidate / "database" / filename
    if not support_path.is_file():
        fallback_path = candidate / filename
        if fallback_path.is_file():
            support_path = fallback_path
    if not support_path.is_file():
        return {"source": {}, "entries": []}
    try:
        payload = json.loads(support_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {"source": {}, "entries": []}
    if not isinstance(payload, dict):
        return {"source": {}, "entries": []}
    entries = payload.get("entries")
    if not isinstance(entries, list):
        payload["entries"] = []
    return payload

def _load_who_bppl_2024_support(project_root_text: str) -> dict:
    return _load_public_health_support_file(project_root_text, "who_bppl_2024_support.json")

def _load_china_cdc_bacteria_support(project_root_text: str) -> dict:
    return _load_public_health_support_file(project_root_text, "china_cdc_bacteria_support.json")

def _load_marker_annotations(payload: dict, entry_key: str) -> dict:
    annotations = payload.get("marker_annotations") if isinstance(payload, dict) else {}
    if not isinstance(annotations, dict):
        return {"key_serotypes": [], "key_resistance_genes": [], "key_virulence_genes": []}
    block = annotations.get(entry_key) if isinstance(annotations.get(entry_key), dict) else {}
    return {
        "key_serotypes": block.get("key_serotypes", []) if isinstance(block.get("key_serotypes"), list) else [],
        "key_resistance_genes": block.get("key_resistance_genes", []) if isinstance(block.get("key_resistance_genes"), list) else [],
        "key_virulence_genes": block.get("key_virulence_genes", []) if isinstance(block.get("key_virulence_genes"), list) else [],
    }

def _load_serotype_profiles(payload: dict, entry_key: str) -> list[dict]:
    profiles = payload.get("serotype_profiles") if isinstance(payload, dict) else {}
    if not isinstance(profiles, dict):
        return []
    matched = profiles.get(entry_key)
    return matched if isinstance(matched, list) else []

def _first_nonempty_from_indexes(row: object, columns: object, candidates: list[str]) -> object:
    if not isinstance(row, list) or not isinstance(columns, list):
        return ""
    for candidate in candidates:
        if candidate not in columns:
            continue
        index = columns.index(candidate)
        if index >= len(row):
            continue
        value = row[index]
        if value is None:
            continue
        text = str(value).strip()
        if text:
            return value
    return ""

def _collect_resistance_signals(rows: dict) -> tuple[str, str]:
    columns = rows.get("columns") if isinstance(rows, dict) else []
    table_rows = rows.get("rows") if isinstance(rows, dict) else []
    if not isinstance(columns, list) or not isinstance(table_rows, list):
        return "", ""
    gene_parts: list[str] = []
    class_parts: list[str] = []
    for row in table_rows:
        if not isinstance(row, list):
            continue
        gene_parts.extend(
            str(_first_nonempty_from_indexes(row, columns, ["基因名称", "耐药基因", "Gene", "Best_Hit_ARO"]) or "").split()
        )
        class_parts.extend(
            str(_first_nonempty_from_indexes(row, columns, ["耐药药物", "Drug Class", "Drug", "Resistance Mechanism"]) or "").split()
        )
    return " ".join(gene_parts).lower(), " ".join(class_parts).lower()

def _phenotype_supported(phenotype: str, gene_blob: str, class_blob: str) -> bool:
    phenotype_key = str(phenotype or "").strip().lower()
    if not phenotype_key:
        return False
    if phenotype_key == "carbapenem-resistant":
        return any(token in gene_blob for token in ("kpc", "ndm", "oxa", "vim", "imp", "ges")) or any(
            token in class_blob for token in ("carbapenem", "碳青霉烯")
        )
    if phenotype_key == "third-generation cephalosporin-resistant":
        return any(token in gene_blob for token in ("ctx-m", "shv", "tem", "ampc", "bla")) or any(
            token in class_blob for token in ("cephalosporin", "头孢", "beta-lactam", "β-内酰胺", "esbl")
        )
    if phenotype_key == "esbl":
        return any(token in gene_blob for token in ("ctx-m", "shv", "tem", "esbl")) or any(
            token in class_blob for token in ("esbl", "头孢", "cephalosporin", "beta-lactam", "β-内酰胺")
        )
    if phenotype_key == "fluoroquinolone-resistant":
        return any(token in gene_blob for token in ("qnr", "gyr", "parc", "pare")) or any(
            token in class_blob for token in ("fluoroquinolone", "quinolone", "喹诺酮", "ciprofloxacin", "levofloxacin")
        )
    if phenotype_key == "vancomycin-resistant":
        return any(token in gene_blob for token in ("vana", "vanb", "vanc")) or "vancomycin" in class_blob
    if phenotype_key == "methicillin-resistant":
        return any(token in gene_blob for token in ("meca", "mecc")) or any(
            token in class_blob for token in ("methicillin", "oxacillin")
        )
    if phenotype_key == "macrolide-resistant":
        return any(token in gene_blob for token in ("erm", "mef", "mph")) or "macrolide" in class_blob
    if phenotype_key == "ampicillin-resistant":
        return "ampicillin" in class_blob or any(token in gene_blob for token in ("tem", "rob", "bla"))
    if phenotype_key == "penicillin-resistant":
        return "penicillin" in class_blob or any(token in gene_blob for token in ("pbp",))
    if phenotype_key == "rifampicin-resistant":
        return any(token in gene_blob for token in ("rpob", "rif")) or "rifampicin" in class_blob
    if phenotype_key == "multidrug-resistant":
        keywords = (
            "carbapenem", "cephalosporin", "fluoroquinolone", "quinolone", "macrolide", "methicillin",
            "ampicillin", "penicillin", "vancomycin", "rifampicin", "头孢", "碳青霉烯", "喹诺酮", "大环内酯",
            "甲氧西林", "氨苄西林", "青霉素", "万古霉素", "利福平", "多重耐药", "mdr", "esbl",
        )
        return sum(1 for token in keywords if token in class_blob or token in gene_blob) >= 2
    return False

def _match_public_health_entries(entries: list, species_name: str, resistance_elements: dict, priority_key: str) -> list[dict]:
    species_text = str(species_name or "").strip().lower()
    if not species_text or not isinstance(entries, list):
        return []
    matched_entries: list[dict] = []
    for entry in entries:
        aliases = entry.get("aliases") if isinstance(entry, dict) else []
        if not isinstance(aliases, list):
            continue
        if any(str(alias or "").strip().lower() in species_text for alias in aliases):
            matched_entries.append(entry)
    if not matched_entries:
        return []
    gene_blob, class_blob = _collect_resistance_signals(resistance_elements)
    priority_rank = {
        "critical": 4,
        "high": 3,
        "medium": 2,
        "法定重点": 4,
        "食源性重点": 3,
        "医院感染重点": 3,
        "症候群重点": 2,
    }
    candidates: list[dict] = []
    for entry in matched_entries:
        phenotypes = entry.get("resistance_phenotypes") if isinstance(entry.get("resistance_phenotypes"), list) else []
        supported = [item for item in phenotypes if _phenotype_supported(str(item), gene_blob, class_blob)]
        candidate = dict(entry)
        candidate["supported_resistance_phenotypes"] = supported
        candidates.append(candidate)
    candidates.sort(
        key=lambda item: (
            len(item.get("supported_resistance_phenotypes") or []),
            priority_rank.get(str(item.get(priority_key) or "").strip(), 0),
        ),
        reverse=True,
    )
    return candidates

def _build_public_health_evidence(
    primary_source: str,
    primary_entry: dict,
    china_candidates: list[dict],
    who_candidates: list[dict],
) -> list[dict]:
    evidence: list[dict] = []
    display_name = str(primary_entry.get("display_name") or "").strip()
    surveillance_domain = str(primary_entry.get("surveillance_domain") or "").strip()
    china_priority = str(primary_entry.get("china_priority_group") or "").strip()
    monitoring_tags = primary_entry.get("monitoring_tags") if isinstance(primary_entry.get("monitoring_tags"), list) else []
    supported_phenotypes = primary_entry.get("supported_resistance_phenotypes") if isinstance(primary_entry.get("supported_resistance_phenotypes"), list) else []
    typing_hint = str(primary_entry.get("typing_hint") or "").strip()
    cluster_hint = str(primary_entry.get("cluster_hint") or "").strip()

    if primary_source == "china":
        if surveillance_domain == "食源性":
            evidence.append(
                {
                    "manual": "2023年食源性疾病监测工作手册-20221226.pdf",
                    "rule_type": "监测对象/暴发溯源",
                    "basis": f"{display_name or '该病原'}命中食源性监测域，需结合食源性病例监测、暴发监测和分子溯源规则判读。",
                }
            )
        if any("国家致病菌识别网" in str(tag) for tag in monitoring_tags) or surveillance_domain in {"呼吸道", "脑膜炎", "医院感染"}:
            evidence.append(
                {
                    "manual": "基于国家致病菌识别网的细菌性传染病实验室监测工作方案（2024版）.pdf",
                    "rule_type": "监测病原/上报时限/事件调查",
                    "basis": f"{display_name or '该病原'}属于{china_priority or '重点'}监测对象，适用症候群监测、重点病原上报和事件溯源调查规则。",
                }
            )
        if typing_hint or cluster_hint or supported_phenotypes:
            evidence.append(
                {
                    "manual": "国家致病菌识别网实验室监测技术手册（2022试用版）终稿(1).pdf",
                    "rule_type": "血清分型/耐药检测/PFGE-WGS",
                    "basis": typing_hint or cluster_hint or f"该病原存在{'; '.join(map(str, supported_phenotypes))}等重点耐药或分型判读依据。",
                }
            )
    if who_candidates:
        top_who = who_candidates[0]
        evidence.append(
            {
                "manual": "9789240093461-eng.pdf",
                "rule_type": "WHO优先病原体分层",
                "basis": f"{display_name or '该病原'}在 WHO 2024 中对应 {top_who.get('who_priority_group', '--')} 优先级，作为国际补充参考。",
            }
        )
    # 去重，避免同一本手册被重复塞两次
    deduped: list[dict] = []
    seen: set[tuple[str, str]] = set()
    for item in evidence:
        key = (str(item.get("manual") or ""), str(item.get("rule_type") or ""))
        if key in seen:
            continue
        seen.add(key)
        deduped.append(item)
    return deduped

def _build_public_health_support(
    project_root: Path,
    species_name: str,
    resistance_elements: dict,
) -> dict:
    database_root = _resolve_runtime_database_root()
    who_payload = _load_who_bppl_2024_support(str(database_root))
    china_payload = _load_china_cdc_bacteria_support(str(database_root))
    species_text = str(species_name or "").strip().lower()
    if not species_text:
        return {
            "status": "empty",
            "source": {"china": china_payload.get("source", {}), "who": who_payload.get("source", {})},
        }
    china_candidates = _match_public_health_entries(china_payload.get("entries"), species_name, resistance_elements, "china_priority_group")
    who_candidates = _match_public_health_entries(who_payload.get("entries"), species_name, resistance_elements, "who_priority_group")
    if not china_candidates and not who_candidates:
        return {
            "status": "unmatched",
            "source": {"china": china_payload.get("source", {}), "who": who_payload.get("source", {})},
            "species_name": species_name,
        }
    primary_source = "china" if china_candidates else "who"
    primary_entry = china_candidates[0] if china_candidates else who_candidates[0]
    marker_annotations = _load_marker_annotations(china_payload if primary_source == "china" else who_payload, str(primary_entry.get("key") or ""))
    serotype_profiles = _load_serotype_profiles(china_payload if primary_source == "china" else who_payload, str(primary_entry.get("key") or ""))
    phenotypes = primary_entry.get("resistance_phenotypes") if isinstance(primary_entry.get("resistance_phenotypes"), list) else []
    supported = primary_entry.get("supported_resistance_phenotypes") or []
    evidence = _build_public_health_evidence(primary_source, primary_entry, china_candidates, who_candidates)
    return {
        "status": "matched",
        "source": {
            "china": china_payload.get("source", {}),
            "who": who_payload.get("source", {}),
        },
        "primary_source": primary_source,
        "species_name": species_name,
        "display_name": primary_entry.get("display_name", species_name),
        "who_priority_group": (who_candidates[0] if who_candidates else {}).get("who_priority_group", ""),
        "china_priority_group": (china_candidates[0] if china_candidates else {}).get("china_priority_group", ""),
        "surveillance_domain": (china_candidates[0] if china_candidates else {}).get("surveillance_domain", ""),
        "monitoring_tags": (china_candidates[0] if china_candidates else {}).get("monitoring_tags", []),
        "resistance_phenotypes": phenotypes,
        "supported_resistance_phenotypes": supported,
        "public_health_significance": primary_entry.get("public_health_significance", ""),
        "transmission_context": primary_entry.get("transmission_context", ""),
        "outbreak_potential": primary_entry.get("outbreak_potential", ""),
        "reporting_hint": primary_entry.get("reporting_hint", ""),
        "cluster_hint": primary_entry.get("cluster_hint", ""),
        "typing_hint": primary_entry.get("typing_hint", ""),
        "resistance_focus": primary_entry.get("resistance_focus", ""),
        "key_serotypes": marker_annotations.get("key_serotypes", []),
        "serotype_profiles": serotype_profiles,
        "key_resistance_genes": marker_annotations.get("key_resistance_genes", []),
        "key_virulence_genes": marker_annotations.get("key_virulence_genes", []),
        "support_evidence": evidence,
        "domestic_matched_entries": [
            {
                "key": item.get("key", ""),
                "display_name": item.get("display_name", ""),
                "surveillance_domain": item.get("surveillance_domain", ""),
                "china_priority_group": item.get("china_priority_group", ""),
                "monitoring_tags": item.get("monitoring_tags", []),
                "resistance_phenotypes": item.get("resistance_phenotypes", []),
                "supported_resistance_phenotypes": item.get("supported_resistance_phenotypes", []),
            }
            for item in china_candidates
        ],
        "matched_entries": [
            {
                "key": item.get("key", ""),
                "display_name": item.get("display_name", ""),
                "who_priority_group": item.get("who_priority_group", ""),
                "resistance_phenotypes": item.get("resistance_phenotypes", []),
                "supported_resistance_phenotypes": item.get("supported_resistance_phenotypes", []),
            }
            for item in who_candidates
        ],
    }
