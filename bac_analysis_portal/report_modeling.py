from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .sample_modeling_service import RISK_LEVEL_LABELS
from .store import PortalStore


def _safe_json(value: Any, fallback: Any) -> Any:
    if isinstance(value, (dict, list)):
        return value
    if not isinstance(value, str) or not value.strip():
        return fallback
    try:
        return json.loads(value)
    except json.JSONDecodeError:
        return fallback


def build_modeling_risk_section(project_root: Path, *, sample_name: str = "", report_dir: Path | str = "") -> dict[str, Any]:
    store = PortalStore(db_path=project_root / "bac_analysis_portal.sqlite3", project_root=project_root)
    score = store.find_latest_sample_modeling_score_for_report(sample_name=sample_name, report_dir=str(report_dir or ""))
    if not score:
        return {
            "status": "empty",
            "headline": "暂无样本库建模评分",
            "summary": "当前报告尚未匹配到样本库建模结果。可在样本库建模面板生成训练集并运行基线模型后回写评分。",
            "items": [],
        }
    level = str(score.get("risk_level") or "routine").strip() or "routine"
    label = str(score.get("risk_label") or RISK_LEVEL_LABELS.get(level, "常规")).strip()
    explanations = _safe_json(score.get("explanation_json"), [])
    run_id = str(score.get("run_id") or "").strip()
    run: dict[str, Any] = {}
    if run_id:
        try:
            run = store.get_sample_modeling_run(run_id)
        except KeyError:
            run = {}
    metrics = _safe_json(run.get("metrics_json"), {})
    created_at = str(run.get("created_at") or score.get("created_at") or "").strip()
    date_part = created_at[:10].replace("-", "") if created_at else "unknown"
    model_version = f"MV-{date_part}-{run_id[-6:]}" if run_id else ""
    algorithm = str(run.get("algorithm") or "").strip()
    disclaimer = "本评分来自样本库基线模型，仅用于复核优先级排序和结果解释辅助；不得单独作为临床诊断、治疗决策或公共卫生处置依据。"
    return {
        "status": "ready",
        "run_id": run_id,
        "model_version": model_version,
        "algorithm": algorithm,
        "model_mode": str(metrics.get("mode") or "").strip(),
        "model_supervised": bool(metrics.get("supervised")),
        "model_sample_count": int(run.get("sample_count") or 0),
        "model_labeled_count": int(run.get("labeled_count") or 0),
        "sample_key": score.get("sample_key") or "",
        "risk_score": round(float(score.get("risk_score") or 0), 1),
        "risk_level": level,
        "risk_label": label,
        "headline": f"样本库建模评分：{label}",
        "summary": f"基于当前样本库基线模型，该样本风险分值为 {round(float(score.get('risk_score') or 0), 1)}。",
        "disclaimer": disclaimer,
        "items": [str(item) for item in explanations if str(item).strip()],
        "feature_snapshot": _safe_json(score.get("feature_snapshot_json"), {}),
        "created_at": score.get("created_at") or "",
    }
