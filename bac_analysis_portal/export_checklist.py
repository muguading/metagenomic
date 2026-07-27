from __future__ import annotations

from typing import Any


def _text(value: Any, default: str = "") -> str:
    text = str(value or "").strip()
    return text or default


def _item(item_id: str, label: str, detail: str, state: str, required: bool = True) -> dict[str, Any]:
    return {
        "id": item_id,
        "label": label,
        "detail": detail,
        "state": state,
        "required": required,
    }


def _count_rows(section: object) -> int:
    if isinstance(section, dict):
        rows = section.get("rows")
        return len(rows) if isinstance(rows, list) else 0
    return 0


def _task_closure(payload: dict[str, Any]) -> dict[str, Any]:
    task = payload.get("task") if isinstance(payload.get("task"), dict) else {}
    closure = task.get("closure_status") if isinstance(task.get("closure_status"), dict) else {}
    return closure


def _closure_step_state(closure: dict[str, Any], step_id: str) -> str:
    steps = closure.get("steps") if isinstance(closure.get("steps"), list) else []
    for step in steps:
        if not isinstance(step, dict):
            continue
        if str(step.get("id") or "").strip() == step_id:
            return str(step.get("state") or "").strip()
    return ""


def _closure_step_done(closure: dict[str, Any], step_id: str) -> bool:
    return _closure_step_state(closure, step_id) == "done"


def build_export_checklist(payload: dict[str, Any]) -> dict[str, Any]:
    task = payload.get("task") if isinstance(payload.get("task"), dict) else {}
    sections = payload.get("sections") if isinstance(payload.get("sections"), dict) else {}
    closure = sections.get("workflow_closure") if isinstance(sections.get("workflow_closure"), dict) else {}
    task_closure = _task_closure(payload)
    public_health = sections.get("public_health_support") if isinstance(sections.get("public_health_support"), dict) else {}
    rv = sections.get("resistance_virulence") if isinstance(sections.get("resistance_virulence"), dict) else {}
    resistance_count = _count_rows((rv.get("resistance_elements") or {}) if isinstance(rv, dict) else {})
    virulence_count = _count_rows((rv.get("virulence_elements") or {}) if isinstance(rv, dict) else {})
    task_id = _text(task.get("id"))
    sample_name = _text(task.get("sample_name") or task.get("name") or task_id, "当前样本")
    risk_level = _text(closure.get("risk_level"), "待评估")
    public_health_matched = _text(public_health.get("status")) == "matched"
    has_task_closure = bool(task_closure)
    task_closure_label = _text(task_closure.get("label"), "待确认")
    task_closure_summary = _text(task_closure.get("summary"))
    review_done = _closure_step_done(task_closure, "review_report")
    database_done = _closure_step_done(task_closure, "import_sample")
    history_done = _closure_step_done(task_closure, "compare_history")
    archive_ready = _text(task_closure.get("state")) in {"needs_archive_export", "closed"} and review_done and database_done and history_done

    report_review_state = "ready" if review_done else ("attention" if has_task_closure else "manual")
    public_health_state = (
        "ready"
        if has_task_closure and review_done
        else ("attention" if public_health_matched or risk_level in {"高", "中"} else "ready")
    )
    database_trace_state = "ready" if database_done else ("attention" if has_task_closure else "manual")
    history_state = "ready" if history_done else ("attention" if has_task_closure else "manual")
    archive_state = "ready" if archive_ready else ("attention" if has_task_closure else "ready")

    items = [
        _item(
            "report_review",
            "报告结论已人工复核",
            (
                f"当前任务闭环状态：{task_closure_label}；{task_closure_summary or '导出前须完成双人复核签核。'}"
                if has_task_closure
                else f"导出前应确认 {sample_name} 的物种/分型/质控结论与正文一致。"
            ),
            report_review_state,
            True,
        ),
        _item(
            "evidence_review",
            "关键证据已检查",
            f"当前报告包含耐药条目 {resistance_count} 项、毒力条目 {virulence_count} 项；导出前应确认表格和摘要口径一致。",
            "ready" if resistance_count or virulence_count else "manual",
            True,
        ),
        _item(
            "public_health_review",
            "公共卫生意义已确认",
            "命中重点/法定/优先病原时，需确认是否触发流调、会商或上报流程。",
            public_health_state,
            True,
        ),
        _item(
            "database_trace",
            "样本库沉淀已完成",
            (
                f"当前任务已入库 {int(task_closure.get('imported_sample_count') or 0)} 份样本；导出版可回溯到样本库记录。"
                if database_done
                else "导出版应能回溯到任务和样本库记录；如尚未入库，应先完成样本库沉淀。"
            ),
            database_trace_state,
            True,
        ),
        _item(
            "history_compare",
            "历史对照已确认",
            (
                "已完成同地区、同型别或同物种历史样本筛查，并形成处置留痕。"
                if history_done
                else "正式导出前应完成历史样本筛查，并记录是否发现聚集、传播或关联线索。"
            ),
            history_state,
            True,
        ),
        _item(
            "export_archive",
            "导出用途与留痕信息已确认",
            (
                "前置闭环动作已完成，本次正式导出会自动记录格式、时间、任务编号和操作者。"
                if archive_ready
                else "只有报告复核、样本入库和历史对照均完成后，才能形成正式归档导出记录。"
            ),
            archive_state,
            True,
        ),
    ]
    blocking_count = sum(1 for item in items if item["required"] and item["state"] in {"manual", "attention"})
    return {
        "status": "ready",
        "title": "导出前疾控复核清单",
        "summary": "导出不是截图保存，而是形成可会商、可归档、可追溯的正式材料；以下项目应在导出前逐项确认。",
        "readiness": "needs_review" if blocking_count else "ready",
        "blocking_count": blocking_count,
        "items": items,
    }
