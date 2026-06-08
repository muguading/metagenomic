from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def analytics_snapshot_path(task_dir: Path) -> Path:
    return task_dir / "task_analytics.json"


def read_task_analytics_snapshot(task_dir: Path) -> dict[str, Any] | None:
    path = analytics_snapshot_path(task_dir)
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None
    return payload if isinstance(payload, dict) else None


def write_task_analytics_snapshot(task_dir: Path, snapshot: dict[str, Any]) -> None:
    path = analytics_snapshot_path(task_dir)
    tmp_path = path.with_suffix(f"{path.suffix}.tmp")
    tmp_path.write_text(json.dumps(snapshot, ensure_ascii=False, indent=2), encoding="utf-8")
    tmp_path.replace(path)


def resolve_queue_task_workstation_key(task: dict[str, Any]) -> str:
    params = task.get("params") if isinstance(task.get("params"), dict) else {}
    explicit = str(params.get("workstation_key") or "").strip().lower()
    if explicit and explicit != "default":
        return explicit
    if str(params.get("method") or "").strip().lower() == "meta":
        return "metagenome"
    if str(params.get("analysis_target") or "").strip().lower() == "virus":
        return "virus"
    return "bacteria"


def get_queue_quality_band_meta(band: str = "") -> dict[str, str]:
    normalized = str(band or "").strip().lower()
    mapping = {
        "success": {"label": "整体稳定", "note": "当前质控或结果基础整体平稳。"},
        "warning": {"label": "建议复核", "note": "建议回看关键质量指标或结果基础。"},
        "danger": {"label": "需重点关注", "note": "当前结果提示需要优先处理的质量风险。"},
        "neutral": {"label": "待生成", "note": "当前还没有足够的摘要用于判断质量层级。"},
    }
    return mapping.get(normalized, mapping["neutral"])


def parse_percent_number(display_value: Any) -> float | None:
    raw = str(display_value or "").replace("%", "").strip()
    if not raw or raw == "--":
        return None
    try:
        value = float(raw)
    except ValueError:
        return None
    return value if value == value else None


def build_queue_analytics_snapshot(report: dict[str, Any]) -> dict[str, Any]:
    metrics = report.get("overview_metrics") if isinstance(report.get("overview_metrics"), list) else []
    metric_map = {
        str(metric.get("key") or ""): metric
        for metric in metrics
        if isinstance(metric, dict) and str(metric.get("key") or "")
    }
    task_meta = report.get("task") if isinstance(report.get("task"), dict) else {}
    task_like = {"params": task_meta}
    workstation_key = resolve_queue_task_workstation_key(task_like)
    quality_band = "neutral"
    dominant_species = ""
    bacteria_quality: dict[str, float] | None = None

    if workstation_key == "metagenome":
        sections = report.get("sections") if isinstance(report.get("sections"), dict) else {}
        binning_results = sections.get("binning_results") if isinstance(sections.get("binning_results"), dict) else {}
        quality_summary = binning_results.get("quality", {}).get("summary", {}) if isinstance(binning_results.get("quality"), dict) else {}
        avg_completeness = _to_number(quality_summary.get("avg_completeness"))
        avg_contamination = _to_number(quality_summary.get("avg_contamination"))
        if avg_completeness is not None and avg_contamination is not None:
            quality_band = "success" if avg_completeness >= 80 and avg_contamination < 10 else "warning"
        dominant_species = _first_metric_display(metric_map.get("meta_species_mge"), 0)
    elif workstation_key == "virus":
        q20 = parse_percent_number(_first_metric_display(metric_map.get("q_metrics"), 0))
        q30 = parse_percent_number(_first_metric_display(metric_map.get("q_metrics"), 1))
        quality_band = "danger" if _lt(q20, 90) or _lt(q30, 80) else "success"
        dominant_species = (
            _first_metric_display(metric_map.get("virus_taxonomy"), 1)
            or _first_metric_display(metric_map.get("virus_taxonomy"), 3)
            or _first_metric_display(metric_map.get("species_estimation"), 0)
        )
    elif workstation_key == "community":
        sections = report.get("sections") if isinstance(report.get("sections"), dict) else {}
        community = sections.get("community") if isinstance(sections.get("community"), dict) else {}
        summary = community.get("summary") if isinstance(community.get("summary"), dict) else {}
        sample_count = _to_int(metric_map.get("community_samples", {}).get("display")) or _to_int(summary.get("sample_count")) or 0
        if sample_count >= 3 and str(summary.get("group_column") or "").strip():
            quality_band = "success"
        elif sample_count > 0:
            quality_band = "warning"
    elif workstation_key == "pathosource":
        sections = report.get("sections") if isinstance(report.get("sections"), dict) else {}
        pathosource = sections.get("pathosource") if isinstance(sections.get("pathosource"), dict) else {}
        snp_matrix = pathosource.get("snp_matrix") if isinstance(pathosource.get("snp_matrix"), dict) else {}
        distance_summary = snp_matrix.get("summary") if isinstance(snp_matrix.get("summary"), dict) else {}
        min_distance = _to_number(distance_summary.get("min_distance"))
        quality_band = "danger" if min_distance is not None and min_distance <= 20 else "neutral"
        dominant_species = str(task_meta.get("species") or task_meta.get("pathosource_species") or "").strip()
    else:
        q20 = parse_percent_number(_first_metric_display(metric_map.get("q_metrics"), 0))
        q30 = parse_percent_number(_first_metric_display(metric_map.get("q_metrics"), 1))
        quality_band = "danger" if _lt(q20, 90) or _lt(q30, 80) else "success"
        dominant_species = _first_metric_display(metric_map.get("species_estimation"), 0)
        checkm_metric = metric_map.get("checkm_metrics") if isinstance(metric_map.get("checkm_metrics"), dict) else {}
        completeness = parse_percent_number(_first_metric_display(checkm_metric, 0))
        contamination = parse_percent_number(_first_metric_display(checkm_metric, 1))
        if completeness is not None and contamination is not None:
            bacteria_quality = {"completeness": completeness, "contamination": contamination}

    if dominant_species in {"--", "False"}:
        dominant_species = ""

    return {
        "ready": True,
        "taskId": str(task_meta.get("id") or ""),
        "taskName": str(task_meta.get("task_name") or task_meta.get("name") or ""),
        "workstationKey": workstation_key,
        "qualityBand": quality_band or "neutral",
        "qualityLabel": get_queue_quality_band_meta(quality_band).get("label", ""),
        "dominantSpecies": dominant_species,
        "bacteriaQuality": bacteria_quality,
    }


def build_pending_queue_analytics_snapshot(task: dict[str, Any]) -> dict[str, Any]:
    workstation_key = resolve_queue_task_workstation_key(task if isinstance(task, dict) else {})
    quality_band = "neutral"
    return {
        "ready": False,
        "taskId": str(task.get("id") or ""),
        "taskName": str(task.get("name") or ""),
        "workstationKey": workstation_key,
        "qualityBand": quality_band,
        "qualityLabel": get_queue_quality_band_meta(quality_band).get("label", ""),
        "dominantSpecies": "",
        "bacteriaQuality": None,
    }


def _first_metric_display(metric: Any, index: int) -> str:
    if not isinstance(metric, dict):
        return ""
    items = metric.get("items")
    if not isinstance(items, list) or index >= len(items) or index < 0:
        return ""
    item = items[index]
    if not isinstance(item, dict):
        return ""
    return str(item.get("display") or "").strip()


def _to_number(value: Any) -> float | None:
    try:
        if value in (None, ""):
            return None
        return float(str(value).strip())
    except (TypeError, ValueError):
        return None


def _to_int(value: Any) -> int | None:
    numeric = _to_number(value)
    return int(numeric) if numeric is not None else None


def _lt(value: float | None, target: float) -> bool:
    return value is not None and value < target
