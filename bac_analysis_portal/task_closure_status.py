from __future__ import annotations

import json
from datetime import datetime
from typing import Any


def _task_status(task: dict[str, Any]) -> str:
    return str(task.get("status") or "").strip().upper()


def _step(step_id: str, label: str, state: str, href: str = "") -> dict[str, str]:
    return {"id": step_id, "label": label, "state": state, "href": href}


def _confirm_action(action_id: str, label: str, detail: str) -> dict[str, str]:
    return {"id": action_id, "label": label, "detail": detail}


def _event_is_after_current_run(event: dict[str, Any], finished_at: object) -> bool:
    cutoff = str(finished_at or "").strip()
    if not cutoff:
        return True
    event_time = str(event.get("created_at") or "").strip()
    if not event_time:
        return False
    try:
        return datetime.fromisoformat(event_time.replace("Z", "+00:00")).timestamp() >= datetime.fromisoformat(cutoff.replace("Z", "+00:00")).timestamp()
    except ValueError:
        return False


def _completed_closure_actions(events: list[dict[str, Any]], *, task_owner: str, finished_at: object = "") -> set[str]:
    action_map = {
        "确认报告复核": "review_report",
        "确认历史对照": "compare_history",
        "确认导出归档": "archive_export",
    }
    completed: set[str] = set()
    for event in events:
        if str(event.get("outcome") or "").strip() != "success":
            continue
        if not _event_is_after_current_run(event, finished_at):
            continue
        action = str(event.get("action") or "").strip()
        action_id = action_map.get(action)
        if not action_id:
            continue
        if action_id == "compare_history" and _event_request_value(event, "comparison_outcome") not in {"no_signal", "signal_found"}:
            continue
        if action_id == "review_report":
            reviewer = str(event.get("username") or "").strip()
            if not task_owner or not reviewer or reviewer == task_owner:
                continue
            if _event_request_value(event, "review_outcome") not in {"approved", "approved_with_notes"}:
                continue
        completed.add(action_id)
    return completed


def _event_request_value(event: dict[str, Any], key: str) -> str:
    return str(_summary_payload(event.get("request_summary")).get(key) or "").strip()


def _latest_history_comparison_outcome(events: list[dict[str, Any]], *, finished_at: object = "") -> str:
    for event in events:
        if str(event.get("action") or "").strip() != "确认历史对照":
            continue
        if str(event.get("outcome") or "").strip() != "success":
            continue
        if not _event_is_after_current_run(event, finished_at):
            continue
        return _event_request_value(event, "comparison_outcome")
    return ""


def _responsibility_trace(events: list[dict[str, Any]], *, task_owner: str, finished_at: object = "") -> dict[str, Any]:
    reviewer = ""
    history_reviewer = ""
    exporter = ""
    for event in events:
        if str(event.get("outcome") or "").strip() != "success":
            continue
        if not _event_is_after_current_run(event, finished_at):
            continue
        action = str(event.get("action") or "").strip()
        username = str(event.get("username") or "").strip()
        if action == "确认报告复核" and username and username != task_owner and not reviewer:
            reviewer = username
        elif action == "确认历史对照" and username and not history_reviewer:
            history_reviewer = username
        elif action == "确认导出归档" and username and not exporter:
            exporter = username
    return {
        "analyst": task_owner,
        "reviewer": reviewer,
        "history_reviewer": history_reviewer,
        "exporter": exporter,
        "separation_verified": bool(task_owner and reviewer and task_owner != reviewer),
    }


def _latest_event_for_action(events: list[dict[str, Any]], action_name: str, *, finished_at: object = "") -> dict[str, Any]:
    for event in events:
        if str(event.get("action") or "").strip() != action_name:
            continue
        if str(event.get("outcome") or "").strip() != "success":
            continue
        if not _event_is_after_current_run(event, finished_at):
            continue
        return event
    return {}


def _delivery_evidence(
    *,
    imported_count: int,
    responsibility: dict[str, Any],
    events: list[dict[str, Any]],
    finished_at: object = "",
) -> list[dict[str, str]]:
    review_event = _latest_event_for_action(events, "确认报告复核", finished_at=finished_at)
    history_event = _latest_event_for_action(events, "确认历史对照", finished_at=finished_at)
    export_event = _latest_event_for_action(events, "确认导出归档", finished_at=finished_at)
    history_outcome = _event_request_value(history_event, "comparison_outcome")
    history_label = {
        "no_signal": "未发现聚集或传播线索",
        "signal_found": "发现聚集或传播线索",
    }.get(history_outcome, "已确认")
    export_format = _event_request_value(export_event, "format").upper()
    return [
        {
            "id": "report_review",
            "label": "报告复核",
            "value": str(responsibility.get("reviewer") or review_event.get("username") or "已确认"),
            "detail": _event_note(review_event) or "物种、分型、质控及公共卫生解释已人工复核。",
        },
        {
            "id": "database_import",
            "label": "样本入库",
            "value": f"{imported_count} 份样本",
            "detail": "结果已沉淀到样本库，可用于后续追溯、统计和历史对照。",
        },
        {
            "id": "history_compare",
            "label": "历史对照",
            "value": str(responsibility.get("history_reviewer") or history_event.get("username") or "已确认"),
            "detail": _event_note(history_event) or history_label,
        },
        {
            "id": "archive_export",
            "label": "正式导出",
            "value": str(responsibility.get("exporter") or export_event.get("username") or "已留痕"),
            "detail": f"正式导出格式：{export_format}" if export_format else (_event_note(export_event) or "报告导出归档已留痕。"),
        },
    ]


def build_task_closure_status(
    task: dict[str, Any],
    *,
    result_exists: bool = False,
    imported_sample_count: int = 0,
    audit_events: list[dict[str, Any]] | None = None,
    eligible_reviewers: list[str] | None = None,
    report_issue: str = "",
) -> dict[str, Any]:
    task_id = str(task.get("id") or "").strip()
    task_owner = str(task.get("owner") or "").strip()
    status = _task_status(task)
    result_href = f"/tasks/{task_id}/result-page" if task_id and result_exists else ""
    database_href = f"/workstation?tab=database&task={task_id}" if task_id else "/workstation?tab=database"
    imported_count = max(0, int(imported_sample_count or 0))
    raw_events = audit_events or []
    recent_events = _build_recent_events(raw_events)
    finished_at = task.get("finished_at")
    completed_actions = _completed_closure_actions(raw_events, task_owner=task_owner, finished_at=finished_at)
    responsibility = _responsibility_trace(raw_events, task_owner=task_owner, finished_at=finished_at)
    responsibility["eligible_reviewers"] = list(eligible_reviewers or [])
    responsibility["eligible_reviewer_count"] = len(eligible_reviewers or [])
    review_done = "review_report" in completed_actions
    history_done = "compare_history" in completed_actions
    archive_done = "archive_export" in completed_actions
    latest_history_outcome = _latest_history_comparison_outcome(raw_events, finished_at=finished_at)
    failure_archived = any(
        str(event.get("action") or "").strip() == "确认失败归档"
        and str(event.get("outcome") or "").strip() == "success"
        and _event_is_after_current_run(event, finished_at)
        for event in raw_events
    )
    stopped_archived = any(
        str(event.get("action") or "").strip() == "确认停止归档"
        and str(event.get("outcome") or "").strip() == "success"
        and _event_is_after_current_run(event, finished_at)
        for event in raw_events
    )
    failure_diagnosis = task.get("failure_diagnosis") if isinstance(task.get("failure_diagnosis"), dict) else {}
    failure_label = str(failure_diagnosis.get("label") or "").strip()
    failure_summary = str(failure_diagnosis.get("summary") or "").strip()

    if status in {"QUEUED", "RUNNING", "PAUSED"}:
        return {
            "state": "analysis_active",
            "label": "分析中",
            "summary": "任务尚未完成，当前闭环重点是观察进度、日志和人工确认节点。",
            "next_action": {"id": "watch_progress", "label": "查看进度", "href": ""},
            "imported_sample_count": imported_count,
            "recent_events": recent_events,
            "responsibility": responsibility,
            "confirmable_actions": [],
            "steps": [
                _step("analysis", "完成分析", "current"),
                _step("review_report", "报告复核", "pending"),
                _step("import_sample", "样本入库", "pending"),
                _step("compare_history", "历史对照", "pending"),
                _step("archive_export", "导出归档", "pending"),
            ],
        }

    if status == "FAILED":
        if failure_archived:
            return {
                "state": "failure_archived",
                "label": "失败已归档",
                "summary": "任务失败已完成原因复核，并确认不再重跑；失败证据与处置依据已留痕。",
                "next_action": {"id": "review_trace", "label": "查看留痕", "href": ""},
                "imported_sample_count": imported_count,
                "recent_events": recent_events,
                "responsibility": responsibility,
                "confirmable_actions": [],
                "steps": [
                    _step("analysis", "分析失败", "blocked"),
                    _step("failure_review", "失败原因复核", "done"),
                    _step("failure_disposition", "处置决定", "done"),
                    _step("failure_archive", "失败归档", "done"),
                ],
            }
        return {
            "state": "analysis_failed",
            "label": "失败待复核",
            "summary": (
                f"任务失败，系统初判为“{failure_label}”：{failure_summary}"
                if failure_label and failure_summary
                else "任务失败，需回看日志后实际重新运行，或填写原因并确认接受失败归档。"
            ),
            "next_action": {"id": "review_log", "label": "查看日志", "href": ""},
            "imported_sample_count": imported_count,
            "recent_events": recent_events,
            "responsibility": responsibility,
            "confirmable_actions": [
                _confirm_action("accept_failure", "接受失败并归档", "确认无需继续重跑，并记录失败原因、影响判断及处置依据。")
            ],
            "steps": [
                _step("analysis", "完成分析", "blocked"),
                _step("review_report", "报告复核", "pending"),
                _step("import_sample", "样本入库", "pending"),
                _step("compare_history", "历史对照", "pending"),
                _step("archive_export", "导出归档", "pending"),
            ],
        }

    if status == "STOPPED":
        if stopped_archived:
            return {
                "state": "stopped_archived",
                "label": "停止已归档",
                "summary": "任务已停止，并完成停止原因复核与终止处置留痕。",
                "next_action": {"id": "review_trace", "label": "查看留痕", "href": ""},
                "imported_sample_count": imported_count,
                "recent_events": recent_events,
                "responsibility": responsibility,
                "confirmable_actions": [],
                "steps": [
                    _step("analysis", "任务停止", "blocked"),
                    _step("stop_review", "停止原因复核", "done"),
                    _step("stop_disposition", "终止决定", "done"),
                    _step("stop_archive", "停止归档", "done"),
                ],
            }
        return {
            "state": "analysis_stopped",
            "label": "已停止",
            "summary": "任务已停止；应实际重建任务，或填写原因并确认停止归档。",
            "next_action": {"id": "decide_rebuild", "label": "确认是否重建", "href": ""},
            "imported_sample_count": imported_count,
            "recent_events": recent_events,
            "responsibility": responsibility,
            "confirmable_actions": [
                _confirm_action("accept_stop", "确认停止并归档", "确认无需继续运行，并记录停止原因、影响判断及终止依据。")
            ],
            "steps": [
                _step("analysis", "完成分析", "blocked"),
                _step("review_report", "报告复核", "pending"),
                _step("import_sample", "样本入库", "pending"),
                _step("compare_history", "历史对照", "pending"),
                _step("archive_export", "导出归档", "pending"),
            ],
        }

    if not result_exists:
        return {
            "state": "report_missing",
            "label": "待生成报告",
            "summary": str(report_issue or "任务已完成，但尚未识别可浏览报告；需先确认输出目录和结果文件。"),
            "next_action": {"id": "open_output", "label": "检查输出", "href": ""},
            "imported_sample_count": imported_count,
            "recent_events": recent_events,
            "responsibility": responsibility,
            "confirmable_actions": [],
            "steps": [
                _step("analysis", "完成分析", "done"),
                _step("review_report", "报告复核", "blocked"),
                _step("import_sample", "样本入库", "pending"),
                _step("compare_history", "历史对照", "pending"),
                _step("archive_export", "导出归档", "pending"),
            ],
        }

    if not review_done and eligible_reviewers is not None and not eligible_reviewers:
        return {
            "state": "reviewer_unavailable",
            "label": "缺复核人",
            "summary": "报告已就绪，但当前没有符合责任分离要求的复核账号；需先配置另一名管理员或对应组管理员。",
            "next_action": {"id": "configure_reviewer", "label": "配置复核人", "href": "/workstation?tab=admin&admin_section=admin-users-section"},
            "imported_sample_count": imported_count,
            "recent_events": recent_events,
            "responsibility": responsibility,
            "confirmable_actions": [],
            "steps": [
                _step("analysis", "完成分析", "done"),
                _step("review_report", "报告复核", "blocked", result_href),
                _step("import_sample", "样本入库", "done" if imported_count > 0 else "pending", database_href),
                _step("compare_history", "历史对照", "pending"),
                _step("archive_export", "导出归档", "pending"),
            ],
        }

    if not review_done:
        return {
            "state": "needs_report_review",
            "label": "待复核",
            "summary": "报告已就绪，但尚未由非任务创建人完成复核签核；双人复核完成后才能继续推进历史对照与归档。",
            "next_action": {"id": "review_report", "label": "确认报告复核", "href": result_href},
            "imported_sample_count": imported_count,
            "recent_events": recent_events,
            "responsibility": responsibility,
            "confirmable_actions": [
                _confirm_action("review_report", "确认报告复核完成", "须由非任务创建人确认物种、分型、质控及公共卫生解释已人工复核。")
            ],
            "steps": [
                _step("analysis", "完成分析", "done"),
                _step("review_report", "报告复核", "current", result_href),
                _step("import_sample", "样本入库", "done" if imported_count > 0 else "pending", database_href),
                _step("compare_history", "历史对照", "pending"),
                _step("archive_export", "导出归档", "pending"),
            ],
        }

    if imported_count <= 0:
        return {
            "state": "needs_database_import",
            "label": "待入库",
            "summary": "报告已完成复核，但样本尚未沉淀到样本库；历史对照和证据追溯还缺关键锚点。",
            "next_action": {"id": "import_sample", "label": "导入数据库", "href": database_href},
            "imported_sample_count": imported_count,
            "recent_events": recent_events,
            "responsibility": responsibility,
            "confirmable_actions": [],
            "steps": [
                _step("analysis", "完成分析", "done"),
                _step("review_report", "报告复核", "done", result_href),
                _step("import_sample", "样本入库", "current", database_href),
                _step("compare_history", "历史对照", "pending"),
                _step("archive_export", "导出归档", "pending"),
            ],
        }

    if not history_done:
        escalation_pending = latest_history_outcome == "needs_escalation"
        return {
            "state": "needs_history_compare",
            "label": "待升级处置" if escalation_pending else "待历史对照",
            "summary": (
                "历史对照已判定需升级复核或处置；完成升级处置并记录结论后，需重新确认历史对照。"
                if escalation_pending
                else f"报告已复核且已有 {imported_count} 份样本入库；需确认已完成同地区、同型别或同物种历史对照。"
            ),
            "next_action": {"id": "compare_history", "label": "历史对照", "href": "/workstation?tab=database"},
            "imported_sample_count": imported_count,
            "recent_events": recent_events,
            "responsibility": responsibility,
            "confirmable_actions": [
                _confirm_action("compare_history", "确认历史对照完成", "确认已完成历史样本筛查，并记录是否发现聚集或传播线索。")
            ],
            "steps": [
                _step("analysis", "完成分析", "done"),
                _step("review_report", "报告复核", "done", result_href),
                _step("import_sample", "样本入库", "done", database_href),
                _step("compare_history", "历史对照", "current", "/workstation?tab=database"),
                _step("archive_export", "导出归档", "pending"),
            ],
        }

    if not archive_done:
        return {
            "state": "needs_archive_export",
            "label": "待归档",
            "summary": "报告复核、样本入库和历史对照均已确认；请从报告页完成正式导出，系统将自动留痕并闭环。",
            "next_action": {"id": "archive_export", "label": "导出归档", "href": result_href},
            "imported_sample_count": imported_count,
            "recent_events": recent_events,
            "responsibility": responsibility,
            "confirmable_actions": [],
            "steps": [
                _step("analysis", "完成分析", "done"),
                _step("review_report", "报告复核", "done", result_href),
                _step("import_sample", "样本入库", "done", database_href),
                _step("compare_history", "历史对照", "done", "/workstation?tab=database"),
                _step("archive_export", "导出归档", "current", result_href),
            ],
        }

    return {
        "state": "closed",
        "label": "已闭环",
        "summary": "报告复核、样本入库、历史对照和导出归档均已确认并留痕。",
        "next_action": {"id": "review_trace", "label": "查看留痕", "href": ""},
        "imported_sample_count": imported_count,
        "delivery_evidence": _delivery_evidence(
            imported_count=imported_count,
            responsibility=responsibility,
            events=raw_events,
            finished_at=finished_at,
        ),
        "recent_events": recent_events,
        "responsibility": responsibility,
        "confirmable_actions": [],
        "steps": [
            _step("analysis", "完成分析", "done"),
            _step("review_report", "报告复核", "done", result_href),
            _step("import_sample", "样本入库", "done", database_href),
            _step("compare_history", "历史对照", "done", "/workstation?tab=database"),
            _step("archive_export", "导出归档", "done", result_href),
        ],
    }


def _build_recent_events(events: list[dict[str, Any]]) -> list[dict[str, str]]:
    recent: list[dict[str, str]] = []
    for event in events[:5]:
        action = str(event.get("action") or "").strip()
        if not action:
            continue
        recent.append(
            {
                "action": action,
                "operator": str(event.get("username") or "").strip(),
                "outcome": str(event.get("outcome") or "").strip(),
                "created_at": str(event.get("created_at") or "").strip(),
                "detail": _event_note(event),
            }
        )
    return recent


def _event_note(event: dict[str, Any]) -> str:
    request_payload = _summary_payload(event.get("request_summary"))
    response_payload = _summary_payload(event.get("response_summary"), nested=False)
    note = str(request_payload.get("note") or "").strip()
    diagnosis_label = str(request_payload.get("diagnosis_label") or response_payload.get("diagnosis_label") or "").strip()
    comparison_labels = {
        "no_signal": "未发现聚集或传播线索",
        "signal_found": "发现聚集或传播线索",
        "needs_escalation": "需升级复核或处置",
    }
    comparison_outcome = str(request_payload.get("comparison_outcome") or response_payload.get("comparison_outcome") or "").strip()
    comparison_label = str(response_payload.get("comparison_outcome_label") or comparison_labels.get(comparison_outcome) or "").strip()
    review_labels = {
        "approved": "复核通过",
        "approved_with_notes": "复核通过但需关注",
        "needs_reanalysis": "需重分析或补充复核",
    }
    review_outcome = str(request_payload.get("review_outcome") or response_payload.get("review_outcome") or "").strip()
    review_label = str(response_payload.get("review_outcome_label") or review_labels.get(review_outcome) or "").strip()
    if note:
        prefix = diagnosis_label or comparison_label or review_label
        return f"{prefix}；{note}" if prefix else note
    if diagnosis_label:
        return diagnosis_label
    if comparison_label:
        return comparison_label
    if review_label:
        return review_label
    export_format = str(request_payload.get("format") or "").strip().upper()
    return f"正式导出格式：{export_format}" if export_format else ""


def _summary_payload(raw_value: object, *, nested: bool = True) -> dict[str, Any]:
    raw = str(raw_value or "").strip()
    if not raw:
        return {}
    try:
        summary = json.loads(raw)
    except (TypeError, ValueError):
        return {}
    if not isinstance(summary, dict):
        return {}
    if not nested:
        return summary
    payload = summary.get("json")
    return payload if isinstance(payload, dict) else {}
