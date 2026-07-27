from __future__ import annotations

from typing import Any


def _text(value: Any, default: str = "") -> str:
    text = str(value or "").strip()
    return text or default


def _evidence(label: str, value: str, state: str = "available") -> dict[str, str]:
    return {"label": label, "value": value, "state": state}


def _action(action_id: str, label: str, href: str = "", state: str = "recommended") -> dict[str, str]:
    return {"id": action_id, "label": label, "href": href, "state": state}


def build_sample_closure_trace(
    record: dict[str, Any],
    *,
    version_events: list[dict[str, Any]] | None = None,
    task_closure_status: dict[str, Any] | None = None,
) -> dict[str, Any]:
    task_id = _text(record.get("task_id"))
    task_name = _text(record.get("task_name"), "未关联任务")
    imported_at = _text(record.get("imported_at"))
    updated_at = _text(record.get("updated_at"))
    report_dir = _text(record.get("report_dir"))
    final_fasta = _text(record.get("final_fasta_path"))
    source = _text(record.get("source_submission_id"))
    report_href = f"/tasks/{task_id}/result-page" if task_id else ""
    task_closure = task_closure_status if isinstance(task_closure_status, dict) else {}
    task_closure_state = _text(task_closure.get("state"))
    task_closure_label = _text(task_closure.get("label"), "未同步")
    task_closed = task_closure_state in {"closed", "failure_archived", "stopped_archived"}

    missing: list[str] = []
    if not task_id:
        missing.append("来源任务")
    if not report_dir:
        missing.append("报告目录")
    if not final_fasta:
        missing.append("Final FASTA")
    if not imported_at:
        missing.append("入库时间")

    if missing:
        state = "incomplete"
        label = "证据链待补齐"
        summary = f"该样本已入库，但缺少 {('、'.join(missing))}；后续追溯、复核或上报前应先补齐。"
    else:
        state = "traceable"
        label = "证据链完整"
        summary = "该样本已关联来源任务、报告目录、Final FASTA 和入库时间，可回到报告复核并开展历史对照。"
        if task_closure_state and not task_closed:
            summary += f" 来源任务当前仍为“{task_closure_label}”，尚不能作为最终交付证据。"

    evidence = [
        _evidence("来源任务", task_name if task_id else "未关联", "available" if task_id else "missing"),
        _evidence("入库时间", imported_at or "未记录", "available" if imported_at else "missing"),
        _evidence("报告目录", report_dir or "未记录", "available" if report_dir else "missing"),
        _evidence("Final FASTA", final_fasta or "未记录", "available" if final_fasta else "missing"),
    ]
    if source:
        evidence.append(_evidence("提交记录", source))
    if updated_at:
        evidence.append(_evidence("最后更新", updated_at))
    if task_closure_state:
        evidence.append(_evidence("来源任务闭环", task_closure_label, "available" if task_closed else "pending"))

    actions = [
        _action("open_report", "回到报告复核", report_href, "available" if report_href else "disabled"),
        _action(
            "compare_history",
            "历史对照",
            f"/workstation?tab=queue&task={task_id}&closure_action=compare_history" if task_id else "/workstation?tab=database",
            "available",
        ),
        _action("complete_metadata", "补齐主档", "", "recommended"),
        _action("export_archive", "导出归档材料", "", "recommended"),
    ]
    if task_id:
        actions.insert(1, _action("open_task", "查看来源任务", f"/workstation?tab=queue&task={task_id}", "available"))

    return {
        "state": state,
        "label": label,
        "summary": summary,
        "task_id": task_id,
        "task_name": task_name,
        "task_closure_state": task_closure_state,
        "task_closure_label": task_closure_label if task_closure_state else "",
        "task_closed": task_closed,
        "report_href": report_href,
        "missing": missing,
        "evidence": evidence,
        "actions": actions,
        "recent_events": _build_recent_events(version_events or []),
    }


def _build_recent_events(events: list[dict[str, Any]]) -> list[dict[str, str]]:
    recent: list[dict[str, str]] = []
    for event in events[:5]:
        action = _text(event.get("action"))
        if not action:
            continue
        recent.append(
            {
                "action": action,
                "summary": _text(event.get("summary")),
                "operator": _text(event.get("operator")),
                "version_label": _text(event.get("version_label")),
                "created_at": _text(event.get("created_at")),
            }
        )
    return recent
