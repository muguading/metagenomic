from __future__ import annotations

import copy
import csv
import json
from pathlib import Path
from typing import Callable

from .admin_runtime import _load_conda_env_settings, _load_conda_root_setting
from .knowledge_interpretation import _normalize_taxonomy_lookup_name
from .parse_utils import _safe_float
from .report_cache import _report_cache_dir
from .report_sources import _resolve_report_source
from .runtime_config import resolve_runtime_env_name as _resolve_runtime_env_name
from .runtime_paths import _resolve_runtime_database_root
from .store import PortalStore, utc_now_iso
from .task_manager import AnalysisTaskManager, ValidationError

ADMIN_PATHOSOURCE_TRIGGER_RULE_DEFAULTS = {
    "enabled": False,
    "priority_only": False,
    "min_abundance_percent": 15.0,
    "min_support_reads": 800,
    "min_coverage_percent": 20.0,
    "allowed_sample_sources": [],
    "priority_species": [],
    "excluded_species": [],
    "max_reference_genomes": 30,
    "msa_method": "snippy",
    "tree_method": "ML",
    "auto_start": True,
}


def _normalize_text_list(raw_value: object) -> list[str]:
    if isinstance(raw_value, list):
        values = raw_value
    else:
        text = str(raw_value or "").replace("；", "\n").replace(";", "\n").replace("，", "\n").replace(",", "\n")
        values = text.splitlines()
    normalized: list[str] = []
    for item in values:
        text = str(item or "").strip()
        if text and text not in normalized:
            normalized.append(text)
    return normalized


def _normalize_admin_pathosource_trigger_rules(loaded: object) -> dict[str, object]:
    defaults = copy.deepcopy(ADMIN_PATHOSOURCE_TRIGGER_RULE_DEFAULTS)
    if not isinstance(loaded, dict):
        return defaults
    normalized = defaults
    normalized["enabled"] = bool(loaded.get("enabled", defaults["enabled"]))
    normalized["priority_only"] = bool(loaded.get("priority_only", defaults["priority_only"]))
    normalized["auto_start"] = bool(loaded.get("auto_start", defaults["auto_start"]))
    normalized["msa_method"] = str(loaded.get("msa_method", defaults["msa_method"]) or defaults["msa_method"]).strip() or defaults["msa_method"]
    normalized["tree_method"] = str(loaded.get("tree_method", defaults["tree_method"]) or defaults["tree_method"]).strip() or defaults["tree_method"]
    normalized["allowed_sample_sources"] = _normalize_text_list(loaded.get("allowed_sample_sources"))
    normalized["priority_species"] = _normalize_text_list(loaded.get("priority_species"))
    normalized["excluded_species"] = _normalize_text_list(loaded.get("excluded_species"))
    try:
        normalized["min_abundance_percent"] = max(0.0, min(100.0, float(loaded.get("min_abundance_percent", defaults["min_abundance_percent"]))))
    except (TypeError, ValueError):
        normalized["min_abundance_percent"] = defaults["min_abundance_percent"]
    try:
        normalized["min_support_reads"] = max(0, int(float(loaded.get("min_support_reads", defaults["min_support_reads"]))))
    except (TypeError, ValueError):
        normalized["min_support_reads"] = defaults["min_support_reads"]
    try:
        normalized["min_coverage_percent"] = max(0.0, min(100.0, float(loaded.get("min_coverage_percent", defaults["min_coverage_percent"]))))
    except (TypeError, ValueError):
        normalized["min_coverage_percent"] = defaults["min_coverage_percent"]
    try:
        normalized["max_reference_genomes"] = max(1, int(float(loaded.get("max_reference_genomes", defaults["max_reference_genomes"]))))
    except (TypeError, ValueError):
        normalized["max_reference_genomes"] = defaults["max_reference_genomes"]
    if normalized["msa_method"] not in {"snippy", "ska"}:
        normalized["msa_method"] = defaults["msa_method"]
    if normalized["tree_method"] not in {"ML", "NJ"}:
        normalized["tree_method"] = defaults["tree_method"]
    return normalized


def _load_admin_pathosource_trigger_rules(store: PortalStore) -> dict[str, object]:
    raw = str(store.get_setting("admin_pathosource_trigger_rules", "") or "").strip()
    if not raw:
        return copy.deepcopy(ADMIN_PATHOSOURCE_TRIGGER_RULE_DEFAULTS)
    try:
        loaded = json.loads(raw)
    except json.JSONDecodeError:
        return copy.deepcopy(ADMIN_PATHOSOURCE_TRIGGER_RULE_DEFAULTS)
    return _normalize_admin_pathosource_trigger_rules(loaded)


def _validate_admin_pathosource_trigger_rules(payload: dict[str, object]) -> dict[str, object]:
    normalized = _normalize_admin_pathosource_trigger_rules(payload)
    if normalized["priority_only"] and not normalized["priority_species"]:
        raise ValidationError("启用“仅重点病原触发”时，请至少填写一个重点病原物种。")
    return normalized


def _save_admin_pathosource_trigger_rules(store: PortalStore, payload: dict[str, object]) -> dict[str, object]:
    normalized = _validate_admin_pathosource_trigger_rules(payload)
    store.set_setting("admin_pathosource_trigger_rules", json.dumps(normalized, ensure_ascii=False))
    return normalized


def _sanitize_runtime_name(value: object, default: str = "item") -> str:
    text = str(value or "").strip() or default
    safe = "".join(ch if ch.isalnum() or ch in {"-", "_", "."} else "_" for ch in text)
    return safe.strip("._-") or default


def _normalize_species_match_key(value: object) -> str:
    return _normalize_taxonomy_lookup_name(value)


def _pathosource_builtin_species_alias(species_name: str) -> str:
    normalized = _normalize_species_match_key(species_name)
    alias_map = {
        "salmonella enterica": "salmonella",
        "salmonella bongori": "salmonella",
        "escherichia coli": "E_coli",
        "shigella sonnei": "Shigella",
        "shigella flexneri": "Shigella",
        "vibrio parahaemolyticus": "Parahemolyticus",
        "vibrio cholerae": "cholerae",
        "yersinia enterocolitica": "Y_enterocolitica",
        "campylobacter jejuni": "Campylobacter",
        "campylobacter coli": "Campylobacter",
        "neisseria meningitidis": "Nmen",
        "listeria monocytogenes": "Lmono",
        "klebsiella pneumoniae": "Kpne",
        "streptococcus suis": "Suare",
        "bacillus cereus": "Bcere",
        "brucella melitensis": "Brucella",
        "brucella abortus": "Brucella",
        "haemophilus influenzae": "HPinf",
    }
    if normalized in alias_map:
        return alias_map[normalized]
    genus = normalized.split(" ", 1)[0] if normalized else ""
    genus_alias_map = {
        "salmonella": "salmonella",
        "escherichia": "E_coli",
        "shigella": "Shigella",
        "vibrio": "cholerae",
        "yersinia": "Y_enterocolitica",
        "campylobacter": "Campylobacter",
        "neisseria": "Nmen",
        "listeria": "Lmono",
        "klebsiella": "Kpne",
        "streptococcus": "Suare",
        "bacillus": "Bcere",
        "brucella": "Brucella",
        "haemophilus": "HPinf",
    }
    return alias_map.get(normalized) or genus_alias_map.get(genus) or species_name


def _estimate_pathosource_trigger_coverage_percent(assembly_coverage: dict) -> float:
    points = assembly_coverage.get("points") if isinstance(assembly_coverage, dict) else []
    numeric_points = []
    for point in points if isinstance(points, list) else []:
        value = _safe_float(point)
        if value is not None:
            numeric_points.append(value)
    if not numeric_points:
        return 0.0
    covered = sum(1 for value in numeric_points if value > 0)
    return round((covered / len(numeric_points)) * 100, 2) if numeric_points else 0.0


def _resolve_auto_pathosource_current_fasta(report_dir: Path, sample_name: str) -> Path | None:
    candidates = [
        report_dir / "tmp_combine.fa",
        report_dir / "tmp_combine.fasta",
        report_dir / "megahit_output" / "final.contigs.fa",
        report_dir / "viral_assembly" / "megahit_output" / "final.contigs.fa",
        report_dir / f"{sample_name}.fa",
        report_dir / f"{sample_name}.fasta",
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    return None


def _select_auto_pathosource_history_records(
    store: PortalStore,
    taxid: str,
    *,
    exclude_task_id: str = "",
    max_items: int = 30,
) -> list[dict[str, object]]:
    target_taxid = str(taxid or "").strip()
    if not target_taxid or target_taxid == "-":
        return []
    matches: list[dict[str, object]] = []
    for item in store.list_sample_library():
        final_fasta_path = Path(str(item.get("final_fasta_path") or "").strip()).expanduser()
        if not final_fasta_path.is_file():
            continue
        if exclude_task_id and str(item.get("task_id") or "").strip() == exclude_task_id:
            continue
        item_taxid = str(item.get("taxid") or "").strip()
        if item_taxid != target_taxid:
            continue
        matches.append(
            {
                "sample_key": str(item.get("sample_key") or "").strip(),
                "sample_name": str(item.get("sample_name") or item.get("sample_key") or "").strip(),
                "task_name": str(item.get("task_name") or item.get("owner") or "library").strip() or "library",
                "final_fasta_path": str(final_fasta_path.resolve()),
                "species_name": str(item.get("species_name") or "").strip(),
                "taxid": item_taxid,
                "updated_at": str(item.get("updated_at") or item.get("imported_at") or "").strip(),
            }
        )
    matches.sort(key=lambda item: str(item.get("updated_at") or ""), reverse=True)
    return matches[: max(1, int(max_items or 1))]


def _write_auto_pathosource_input_sheet(
    cache_dir: Path,
    sample_name: str,
    current_fasta: Path,
    history_records: list[dict[str, object]],
) -> Path:
    cache_dir.mkdir(parents=True, exist_ok=True)
    target = cache_dir / f"{_sanitize_runtime_name(sample_name, 'sample')}_pathosource_input.tsv"
    rows = [["样本名称", "单端数据", "右端数据"]]
    rows.append([sample_name, str(current_fasta), "none"])
    seen_names = {sample_name}
    for index, item in enumerate(history_records, start=1):
        base_name = str(item.get("sample_name") or item.get("sample_key") or f"history_{index}").strip() or f"history_{index}"
        history_name = base_name
        suffix = 1
        while history_name in seen_names:
            suffix += 1
            history_name = f"{base_name}_{suffix}"
        seen_names.add(history_name)
        rows.append([history_name, str(item.get("final_fasta_path") or "").strip(), "none"])
    with target.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerows(rows)
    return target


def _pathosource_trigger_rule_summary(rules: dict[str, object]) -> dict[str, object]:
    return {
        "enabled": bool(rules.get("enabled")),
        "priority_only": bool(rules.get("priority_only")),
        "auto_start": bool(rules.get("auto_start")),
        "min_abundance_percent": float(rules.get("min_abundance_percent") or 0.0),
        "min_support_reads": int(rules.get("min_support_reads") or 0),
        "min_coverage_percent": float(rules.get("min_coverage_percent") or 0.0),
        "max_reference_genomes": int(rules.get("max_reference_genomes") or 0),
        "allowed_sample_sources": list(rules.get("allowed_sample_sources") or []),
        "priority_species": list(rules.get("priority_species") or []),
        "excluded_species": list(rules.get("excluded_species") or []),
        "msa_method": str(rules.get("msa_method") or ""),
        "tree_method": str(rules.get("tree_method") or ""),
    }


def _pathosource_trigger_check(
    key: str,
    label: str,
    status: str,
    detail: str,
    *,
    value: object = "",
    threshold: object = "",
) -> dict[str, object]:
    return {
        "key": key,
        "label": label,
        "status": status,
        "detail": detail,
        "value": value,
        "threshold": threshold,
    }


def _build_auto_pathosource_trigger_payload(
    status: str,
    reason: str,
    *,
    rules: dict[str, object] | None = None,
    sample_name: str = "",
    sample_source: str = "",
    candidate: dict[str, object] | None = None,
    checks: list[dict[str, object]] | None = None,
    evaluated_candidates: list[dict[str, object]] | None = None,
    child_task: dict[str, object] | None = None,
    history_count: int | None = None,
    current_fasta: str = "",
    input_sheet: str = "",
    output_dir: str = "",
    coverage_percent: float | None = None,
    extra: dict[str, object] | None = None,
) -> dict[str, object]:
    payload: dict[str, object] = {
        "status": status,
        "reason": reason,
        "sample_name": sample_name,
        "sample_source": sample_source,
        "checks": checks or [],
        "evaluated_candidates": evaluated_candidates or [],
    }
    if rules is not None:
        payload["rules"] = _pathosource_trigger_rule_summary(rules)
    if candidate:
        payload["trigger_species"] = candidate.get("species_name", "")
        payload["trigger_taxid"] = candidate.get("taxid", "")
        payload["candidate"] = candidate
    if child_task:
        payload["child_task"] = child_task
        payload["child_task_id"] = child_task.get("id", "")
        payload["child_task_name"] = child_task.get("name", "")
    if history_count is not None:
        payload["history_count"] = history_count
    if current_fasta:
        payload["current_fasta"] = current_fasta
    if input_sheet:
        payload["input_sheet"] = input_sheet
    if output_dir:
        payload["output_dir"] = output_dir
    if coverage_percent is not None:
        payload["coverage_percent"] = round(float(coverage_percent), 2)
    if extra:
        payload.update(extra)
    return payload


def _pick_auto_pathosource_candidate(
    species_taxonomy: dict,
    rules: dict[str, object],
    *,
    coverage_percent: float,
    sample_source: str = "",
) -> tuple[dict[str, object] | None, str, dict[str, object]]:
    rows = species_taxonomy.get("rows") if isinstance(species_taxonomy, dict) else []
    if not isinstance(rows, list) or not rows:
        return None, "未读取到可用于自动溯源判定的物种结果。", {
            "checks": [
                _pathosource_trigger_check("species_results", "物种结果", "failed", "未读取到种水平分类结果，无法执行自动溯源判定。"),
            ],
            "evaluated_candidates": [],
        }
    allowed_sources = {_normalize_species_match_key(item) for item in (rules.get("allowed_sample_sources") or []) if str(item).strip()}
    if allowed_sources and _normalize_species_match_key(sample_source) not in allowed_sources:
        return None, f"当前样本来源“{sample_source or '-'}”不在允许自动溯源的样本来源范围内。", {
            "checks": [
                _pathosource_trigger_check(
                    "sample_source",
                    "样本来源",
                    "failed",
                    f"当前样本来源“{sample_source or '-'}”不在允许范围内。",
                    value=sample_source or "-",
                    threshold="、".join(str(item) for item in rules.get("allowed_sample_sources") or []),
                ),
            ],
            "evaluated_candidates": [],
        }
    excluded = {_normalize_species_match_key(item) for item in (rules.get("excluded_species") or []) if str(item).strip()}
    priority = {_normalize_species_match_key(item) for item in (rules.get("priority_species") or []) if str(item).strip()}
    min_ratio = float(rules.get("min_abundance_percent") or 0.0)
    min_reads = int(rules.get("min_support_reads") or 0)
    min_coverage = float(rules.get("min_coverage_percent") or 0.0)
    ranked_rows = sorted(
        [row for row in rows if isinstance(row, dict)],
        key=lambda item: (
            float(item.get("比例数值") or 0.0),
            int(item.get("序列数量数值") or 0),
        ),
        reverse=True,
    )
    source_check = _pathosource_trigger_check(
        "sample_source",
        "样本来源",
        "passed",
        "样本来源未限制，或当前样本来源在允许范围内。",
        value=sample_source or "-",
        threshold="、".join(str(item) for item in rules.get("allowed_sample_sources") or []) or "不限",
    )
    evaluated_candidates: list[dict[str, object]] = []
    for row in ranked_rows:
        species_name = str(row.get("种") or row.get("NCBI种") or "").strip()
        if not species_name or species_name == "-":
            continue
        normalized_species = _normalize_species_match_key(species_name)
        ratio = float(row.get("比例数值") or 0.0)
        reads = int(row.get("序列数量数值") or 0)
        candidate_checks = [source_check]
        rejection_reasons: list[str] = []
        if normalized_species in excluded:
            rejection_reasons.append("在排除物种名单中")
            candidate_checks.append(_pathosource_trigger_check("excluded_species", "排除名单", "failed", "候选物种在排除名单中。", value=species_name))
        else:
            candidate_checks.append(_pathosource_trigger_check("excluded_species", "排除名单", "passed", "候选物种未命中排除名单。", value=species_name))
        if bool(rules.get("priority_only")) and normalized_species not in priority:
            rejection_reasons.append("未在重点病原名单中")
            candidate_checks.append(_pathosource_trigger_check("priority_species", "重点病原", "failed", "已启用仅重点病原触发，但候选物种不在名单中。", value=species_name))
        elif bool(rules.get("priority_only")):
            candidate_checks.append(_pathosource_trigger_check("priority_species", "重点病原", "passed", "候选物种在重点病原名单中。", value=species_name))
        else:
            candidate_checks.append(_pathosource_trigger_check("priority_species", "重点病原", "passed", "未启用仅重点病原触发。", value="不限"))
        if ratio < min_ratio:
            rejection_reasons.append(f"丰度 {ratio:.2f}% < {min_ratio:.2f}%")
            candidate_checks.append(_pathosource_trigger_check("min_abundance_percent", "相对丰度", "failed", "候选物种相对丰度低于阈值。", value=f"{ratio:.2f}%", threshold=f"{min_ratio:.2f}%"))
        else:
            candidate_checks.append(_pathosource_trigger_check("min_abundance_percent", "相对丰度", "passed", "候选物种相对丰度达到阈值。", value=f"{ratio:.2f}%", threshold=f"{min_ratio:.2f}%"))
        if reads < min_reads:
            rejection_reasons.append(f"支持 reads {reads} < {min_reads}")
            candidate_checks.append(_pathosource_trigger_check("min_support_reads", "支持 reads", "failed", "候选物种支持 reads 低于阈值。", value=reads, threshold=min_reads))
        else:
            candidate_checks.append(_pathosource_trigger_check("min_support_reads", "支持 reads", "passed", "候选物种支持 reads 达到阈值。", value=reads, threshold=min_reads))
        if coverage_percent < min_coverage:
            rejection_reasons.append(f"覆盖度 {coverage_percent:.2f}% < {min_coverage:.2f}%")
            candidate_checks.append(_pathosource_trigger_check("min_coverage_percent", "参考覆盖度", "failed", "当前报告覆盖度低于阈值。", value=f"{coverage_percent:.2f}%", threshold=f"{min_coverage:.2f}%"))
        else:
            candidate_checks.append(_pathosource_trigger_check("min_coverage_percent", "参考覆盖度", "passed", "当前报告覆盖度达到阈值。", value=f"{coverage_percent:.2f}%", threshold=f"{min_coverage:.2f}%"))
        candidate_summary = {
            "species_name": species_name,
            "ratio": round(ratio, 2),
            "reads": reads,
            "coverage_percent": round(coverage_percent, 2),
            "taxid": str(row.get("NCBI TaxID") or "").strip(),
            "is_priority": normalized_species in priority,
            "decision": "rejected" if rejection_reasons else "selected",
            "reason": "；".join(rejection_reasons) if rejection_reasons else "达到自动溯源触发条件。",
        }
        if len(evaluated_candidates) < 6:
            evaluated_candidates.append(candidate_summary)
        if rejection_reasons:
            continue
        return {
            "species_name": species_name,
            "ratio": round(ratio, 2),
            "reads": reads,
            "coverage_percent": round(coverage_percent, 2),
            "genus": str(row.get("属") or "").strip(),
            "taxid": str(row.get("NCBI TaxID") or "").strip(),
            "is_priority": normalized_species in priority,
        }, "", {"checks": candidate_checks, "evaluated_candidates": evaluated_candidates}
    return None, "当前物种结果尚未达到自动溯源阈值。", {
        "checks": [
            source_check,
            _pathosource_trigger_check("candidate", "候选物种", "failed", "没有物种同时满足名单、丰度、reads 和覆盖度条件。"),
        ],
        "evaluated_candidates": evaluated_candidates,
    }

# admin runtime constants moved to bac_analysis_portal.admin_runtime.

# admin runtime helpers moved to bac_analysis_portal.admin_runtime.

def maybe_auto_trigger_pathosource_for_meta_task(
    task: dict[str, object],
    payload: dict[str, object],
    *,
    store: PortalStore,
    task_manager: AnalysisTaskManager,
    project_root: Path,
    owner_fallback: Callable[[], str],
    selected_sample: str = "",
) -> dict[str, object]:
    params = task.get("params") if isinstance(task.get("params"), dict) else {}
    workstation_key = str(params.get("workstation_key") or "").strip().lower()
    method = str(params.get("method") or "").strip().lower()
    if workstation_key != "metagenome" and method != "meta":
        payload["auto_pathosource_trigger"] = _build_auto_pathosource_trigger_payload(
            "skipped",
            "当前任务不是宏基因组任务，未执行自动溯源判定。",
        )
        return payload
    if str(task.get("status") or "").strip().upper() != "SUCCEEDED":
        payload["auto_pathosource_trigger"] = _build_auto_pathosource_trigger_payload(
            "skipped",
            "任务尚未完成，当前不会自动派生溯源子任务。",
        )
        return payload
    rules = _load_admin_pathosource_trigger_rules(store)
    if not bool(rules.get("enabled")):
        payload["auto_pathosource_trigger"] = _build_auto_pathosource_trigger_payload(
            "disabled",
            "后台未启用自动溯源触发规则。",
            rules=rules,
        )
        return payload

    report_task = payload.get("task") if isinstance(payload.get("task"), dict) else {}
    report_source = _resolve_report_source(task, selected_sample)
    if not report_source.get("available"):
        payload["auto_pathosource_trigger"] = _build_auto_pathosource_trigger_payload(
            "skipped",
            str(report_source.get("reason") or "当前报告尚未定位到服务器结果目录。").strip(),
            rules=rules,
        )
        return payload
    report_dir = Path(report_source.get("report_dir") or "").expanduser().resolve() if report_source.get("report_dir") else None
    sample_name = str(report_source.get("selected_sample") or report_task.get("sample_name") or selected_sample or "").strip()
    report_output_dir = str(report_task.get("output_dir") or params.get("output_dir") or (str(report_dir) if report_dir else "")).strip()
    if not sample_name or report_dir is None or not report_dir.is_dir():
        payload["auto_pathosource_trigger"] = _build_auto_pathosource_trigger_payload(
            "skipped",
            "当前报告尚未定位到有效的样本名称或报告目录。",
            rules=rules,
            sample_name=sample_name,
        )
        return payload

    trigger_key = f"{sample_name}::{selected_sample or sample_name}"
    auto_pathosource_state = task.get("auto_pathosource") if isinstance(task.get("auto_pathosource"), dict) else {}
    existing_state = auto_pathosource_state.get(trigger_key) if isinstance(auto_pathosource_state, dict) else None
    if isinstance(existing_state, dict):
        if str(existing_state.get("status") or "").strip() == "pending":
            payload["auto_pathosource_trigger"] = _build_auto_pathosource_trigger_payload(
                "pending",
                str(existing_state.get("reason") or "正在创建自动溯源子任务，请稍后刷新。").strip(),
                rules=rules,
                sample_name=sample_name,
                candidate={"species_name": existing_state.get("trigger_species", ""), "taxid": existing_state.get("trigger_taxid", "")},
                checks=existing_state.get("checks") if isinstance(existing_state.get("checks"), list) else [],
                history_count=int(existing_state.get("history_count") or 0),
            )
            return payload
        child_task_id = str(existing_state.get("child_task_id") or "").strip()
        if child_task_id:
            try:
                child_task = task_manager.get_task(child_task_id, log_lines=0, owner=None)
            except KeyError:
                child_task = None
            if child_task is not None:
                payload["auto_pathosource_trigger"] = _build_auto_pathosource_trigger_payload(
                    "linked",
                    str(existing_state.get("reason") or "已存在自动创建的溯源子任务。").strip(),
                    rules=rules,
                    sample_name=sample_name,
                    candidate={"species_name": existing_state.get("trigger_species", ""), "taxid": existing_state.get("trigger_taxid", "")},
                    checks=existing_state.get("checks") if isinstance(existing_state.get("checks"), list) else [],
                    history_count=int(existing_state.get("history_count") or 0),
                    input_sheet=str(existing_state.get("input_sheet") or ""),
                    output_dir=str(existing_state.get("output_dir") or ""),
                    child_task={
                        "id": child_task_id,
                        "name": child_task.get("name"),
                        "status": child_task.get("status"),
                        "output_dir": str((child_task.get("params") or {}).get("output_dir") or ""),
                    },
                )
                return payload
        previous_status = str(existing_state.get("status") or "").strip()
        if previous_status in {"failed", "blocked", "would_trigger"}:
            payload["auto_pathosource_trigger"] = _build_auto_pathosource_trigger_payload(
                previous_status,
                str(existing_state.get("reason") or "自动溯源触发已留下历史状态。").strip(),
                rules=rules,
                sample_name=sample_name,
                candidate={"species_name": existing_state.get("trigger_species", ""), "taxid": existing_state.get("trigger_taxid", "")},
                checks=existing_state.get("checks") if isinstance(existing_state.get("checks"), list) else [],
                history_count=int(existing_state.get("history_count") or 0),
                current_fasta=str(existing_state.get("current_fasta") or ""),
                input_sheet=str(existing_state.get("input_sheet") or ""),
                output_dir=str(existing_state.get("output_dir") or ""),
            )
            return payload

    sections = payload.get("sections") if isinstance(payload.get("sections"), dict) else {}
    taxonomy_section = sections.get("taxonomy") if isinstance(sections.get("taxonomy"), dict) else {}
    species_taxonomy = taxonomy_section.get("species_taxonomy") if isinstance(taxonomy_section.get("species_taxonomy"), dict) else {"rows": []}
    coverage_section = sections.get("coverage") if isinstance(sections.get("coverage"), dict) else {"points": []}
    coverage_percent = _estimate_pathosource_trigger_coverage_percent(coverage_section)
    sample_source = str(params.get("sample_source") or report_task.get("sample_source") or "").strip()
    candidate, candidate_reason, candidate_evaluation = _pick_auto_pathosource_candidate(
        species_taxonomy,
        rules,
        coverage_percent=coverage_percent,
        sample_source=sample_source,
    )
    if candidate is None:
        payload["auto_pathosource_trigger"] = _build_auto_pathosource_trigger_payload(
            "not_triggered",
            candidate_reason,
            rules=rules,
            sample_name=sample_name,
            sample_source=sample_source,
            checks=candidate_evaluation.get("checks") if isinstance(candidate_evaluation.get("checks"), list) else [],
            evaluated_candidates=candidate_evaluation.get("evaluated_candidates") if isinstance(candidate_evaluation.get("evaluated_candidates"), list) else [],
            coverage_percent=coverage_percent,
        )
        return payload

    checks = candidate_evaluation.get("checks") if isinstance(candidate_evaluation.get("checks"), list) else []
    evaluated_candidates = candidate_evaluation.get("evaluated_candidates") if isinstance(candidate_evaluation.get("evaluated_candidates"), list) else []
    if not bool(rules.get("auto_start")):
        payload["auto_pathosource_trigger"] = _build_auto_pathosource_trigger_payload(
            "would_trigger",
            "已命中自动溯源阈值，但后台关闭了自动启动子任务，仅生成触发建议。",
            rules=rules,
            sample_name=sample_name,
            sample_source=sample_source,
            candidate=candidate,
            checks=checks,
            evaluated_candidates=evaluated_candidates,
            coverage_percent=coverage_percent,
        )
        return payload

    current_fasta = _resolve_auto_pathosource_current_fasta(report_dir, sample_name)
    if current_fasta is None:
        payload["auto_pathosource_trigger"] = _build_auto_pathosource_trigger_payload(
            "blocked",
            "已命中自动溯源阈值，但当前报告目录中没有找到可用于 PathoSource 的 fasta 文件。",
            rules=rules,
            sample_name=sample_name,
            sample_source=sample_source,
            candidate=candidate,
            checks=checks + [
                _pathosource_trigger_check("current_fasta", "当前样本 fasta", "failed", "报告目录中未找到 tmp_combine.fa、final.contigs.fa 或样本同名 fasta。"),
            ],
            evaluated_candidates=evaluated_candidates,
            coverage_percent=coverage_percent,
        )
        return payload

    candidate_taxid = str(candidate.get("taxid") or "").strip()
    if not candidate_taxid or candidate_taxid == "-":
        payload["auto_pathosource_trigger"] = _build_auto_pathosource_trigger_payload(
            "blocked",
            "已命中自动溯源阈值，但当前候选物种缺少有效 TaxID，未创建子任务。",
            rules=rules,
            sample_name=sample_name,
            sample_source=sample_source,
            candidate=candidate,
            checks=checks + [
                _pathosource_trigger_check("trigger_taxid", "候选 TaxID", "failed", "候选物种没有可用于检索历史样本的 TaxID。", value=candidate_taxid or "-"),
            ],
            evaluated_candidates=evaluated_candidates,
            current_fasta=str(current_fasta),
            coverage_percent=coverage_percent,
        )
        return payload

    history_records = _select_auto_pathosource_history_records(
        store,
        candidate_taxid,
        exclude_task_id=str(task.get("id") or ""),
        max_items=int(rules.get("max_reference_genomes") or 30),
    )
    if not history_records:
        payload["auto_pathosource_trigger"] = _build_auto_pathosource_trigger_payload(
            "blocked",
            f"已命中自动溯源阈值，但样本库中暂未找到 TaxID {candidate_taxid} 的历史 fasta，未创建子任务。",
            rules=rules,
            sample_name=sample_name,
            sample_source=sample_source,
            candidate=candidate,
            checks=checks + [
                _pathosource_trigger_check("current_fasta", "当前样本 fasta", "passed", "已找到当前样本可用于 PathoSource 的 fasta。", value=str(current_fasta)),
                _pathosource_trigger_check("history_records", "历史株", "failed", f"样本库中没有 TaxID {candidate_taxid} 的可用历史 fasta。", value=0, threshold=f"1-{int(rules.get('max_reference_genomes') or 30)}"),
            ],
            evaluated_candidates=evaluated_candidates,
            current_fasta=str(current_fasta),
            history_count=0,
            coverage_percent=coverage_percent,
        )
        return payload

    checks = checks + [
        _pathosource_trigger_check("current_fasta", "当前样本 fasta", "passed", "已找到当前样本可用于 PathoSource 的 fasta。", value=str(current_fasta)),
        _pathosource_trigger_check("history_records", "历史株", "passed", f"已找到 {len(history_records)} 条同 TaxID 历史 fasta。", value=len(history_records), threshold=f"1-{int(rules.get('max_reference_genomes') or 30)}"),
    ]
    next_state = dict(auto_pathosource_state or {})
    next_state[trigger_key] = {
        "status": "pending",
        "trigger_species": str(candidate.get("species_name") or ""),
        "trigger_taxid": candidate_taxid,
        "history_count": len(history_records),
        "created_at": utc_now_iso(),
        "reason": "已命中自动溯源阈值，正在创建 PathoSource 子任务。",
        "checks": checks,
    }
    task_manager.update_task_fields(str(task.get("id") or ""), {"auto_pathosource": next_state}, owner=None)

    cache_dir = _report_cache_dir(report_dir) / "auto_pathosource"
    input_sheet = _write_auto_pathosource_input_sheet(
        cache_dir,
        sample_name,
        current_fasta,
        history_records,
    )
    preferred_msa_method = str(rules.get("msa_method") or "snippy").strip() or "snippy"
    reference_path = str(history_records[0].get("final_fasta_path") or "").strip()
    child_species_value = _pathosource_builtin_species_alias(str(candidate.get("species_name") or ""))
    if preferred_msa_method == "snippy" and not reference_path:
        preferred_msa_method = "ska"
    child_output_dir = (Path(report_output_dir) / "auto_pathosource" / _sanitize_runtime_name(str(candidate.get("species_name") or sample_name), "pathosource")).resolve()
    child_output_dir.mkdir(parents=True, exist_ok=True)
    child_payload = {
        "task_name": f"{task.get('name') or sample_name}_auto_trace_{_sanitize_runtime_name(str(candidate.get('species_name') or sample_name), 'species')}",
        "input_path": str(input_sheet),
        "output_dir": str(child_output_dir),
        "thread": max(4, min(16, int(params.get("thread") or 8))),
        "species": child_species_value,
        "ref": reference_path if preferred_msa_method == "snippy" and reference_path else "False",
        "pathosource_ref": reference_path if preferred_msa_method == "snippy" and reference_path else "False",
        "pathosource_cgmlstana": "no",
        "pathosource_gubbins": "yes",
        "pathosource_msamethod": preferred_msa_method,
        "pathosource_treemethod": str(rules.get("tree_method") or "ML"),
        "pathosource_bootstrap": 1000,
        "pathosource_mltype": "MFP",
        "pathosource_mode": "P",
        "pathosource_cgmlstversion": "none",
        "workstation_key": "pathosource",
    }
    try:
        child_task = task_manager.create_task(
            child_payload,
            owner=str(task.get("owner") or owner_fallback() or ""),
            owner_group=str(task.get("owner_group") or ""),
            pipeline_script=str((project_root / "PathoSource.py").resolve()),
            pipeline_python=_resolve_runtime_env_name(project_root, store.get_setting("pipeline_python", "base")),
            database_root=str(_resolve_runtime_database_root()),
            conda_root=_load_conda_root_setting(store),
            max_concurrent_tasks=int(store.get_setting("max_concurrent_tasks", "2") or "2"),
            conda_envs=_load_conda_env_settings(store),
            extra_task_fields={
                "parent_task_id": str(task.get("id") or ""),
                "trigger_context": {
                    "source": "metagenome_auto_pathosource",
                    "sample_name": sample_name,
                    "trigger_species": str(candidate.get("species_name") or ""),
                    "trigger_taxid": candidate_taxid,
                    "input_sheet": str(input_sheet),
                    "history_count": len(history_records),
                    "current_fasta": str(current_fasta),
                },
            },
        )
    except Exception as exc:
        next_state[trigger_key] = {
            "status": "failed",
            "trigger_species": str(candidate.get("species_name") or ""),
            "trigger_taxid": candidate_taxid,
            "history_count": len(history_records),
            "created_at": utc_now_iso(),
            "reason": f"自动创建 PathoSource 子任务失败：{exc}",
            "checks": checks,
            "current_fasta": str(current_fasta),
            "input_sheet": str(input_sheet),
            "output_dir": str(child_output_dir),
        }
        task_manager.update_task_fields(str(task.get("id") or ""), {"auto_pathosource": next_state}, owner=None)
        payload["auto_pathosource_trigger"] = _build_auto_pathosource_trigger_payload(
            "failed",
            f"已命中自动溯源阈值，但创建 PathoSource 子任务失败：{exc}",
            rules=rules,
            sample_name=sample_name,
            sample_source=sample_source,
            candidate=candidate,
            checks=checks + [
                _pathosource_trigger_check("child_task", "子任务创建", "failed", str(exc)),
            ],
            evaluated_candidates=evaluated_candidates,
            current_fasta=str(current_fasta),
            input_sheet=str(input_sheet),
            output_dir=str(child_output_dir),
            history_count=len(history_records),
            coverage_percent=coverage_percent,
        )
        return payload
    next_state[trigger_key] = {
        "status": "created",
        "child_task_id": str(child_task.get("id") or ""),
        "child_task_name": str(child_task.get("name") or ""),
        "trigger_species": str(candidate.get("species_name") or ""),
        "trigger_taxid": candidate_taxid,
        "history_count": len(history_records),
        "created_at": utc_now_iso(),
        "reason": "命中宏基因组自动溯源阈值，已创建 PathoSource 子任务。",
        "checks": checks,
        "current_fasta": str(current_fasta),
        "input_sheet": str(input_sheet),
        "output_dir": str(child_output_dir),
    }
    task_manager.update_task_fields(str(task.get("id") or ""), {"auto_pathosource": next_state}, owner=None)
    payload["auto_pathosource_trigger"] = _build_auto_pathosource_trigger_payload(
        "created",
        "已命中自动溯源阈值，并创建 PathoSource 子任务。",
        rules=rules,
        sample_name=sample_name,
        sample_source=sample_source,
        candidate=candidate,
        checks=checks + [
            _pathosource_trigger_check("child_task", "子任务创建", "passed", "PathoSource 子任务已创建。", value=str(child_task.get("id") or "")),
        ],
        evaluated_candidates=evaluated_candidates,
        child_task={
            "id": child_task.get("id"),
            "name": child_task.get("name"),
            "status": child_task.get("status"),
            "output_dir": str(child_payload.get("output_dir") or ""),
        },
        current_fasta=str(current_fasta),
        input_sheet=str(input_sheet),
        output_dir=str(child_output_dir),
        history_count=len(history_records),
        coverage_percent=coverage_percent,
    )
    return payload
