from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Callable

from .identity import UserIdentity
from .sample_library_manager import SampleLibraryManager
from .store import PortalStore
from .task_manager import AnalysisTaskManager, ValidationError


CLOSURE_ACTION_LABELS = {
    "review_report": "报告复核完成",
    "compare_history": "历史对照完成",
}

REPORT_EXPORT_FORMATS = {"pdf", "word", "html"}
REPORT_REVIEW_OUTCOMES = {
    "approved": "复核通过",
    "approved_with_notes": "复核通过但需关注",
    "needs_reanalysis": "需重分析或补充复核",
}
COMPLETED_REPORT_REVIEW_OUTCOMES = {"approved", "approved_with_notes"}
HISTORY_COMPARISON_OUTCOMES = {
    "no_signal": "未发现聚集或传播线索",
    "signal_found": "发现聚集或传播线索",
    "needs_escalation": "需升级复核或处置",
}
COMPLETED_HISTORY_COMPARISON_OUTCOMES = {"no_signal", "signal_found"}


def _next_step(action_id: str, label: str, detail: str = "") -> dict[str, str]:
    return {"id": action_id, "label": label, "detail": detail}


def _review_confirmation_guidance(review_outcome: str, *, imported_count: int) -> tuple[str, dict[str, str]]:
    if review_outcome == "needs_reanalysis":
        return (
            "已记录“需重分析或补充复核”，后续入库、历史对照和正式归档仍保持阻塞。",
            _next_step("decide_rebuild", "重建或补充复核", "请先修正报告依据，再由非任务创建人重新完成复核签核。"),
        )
    if imported_count <= 0:
        return (
            "报告复核已留痕，下一步请将结果沉淀到样本库。",
            _next_step("import_sample", "导入样本库", "样本入库后才能形成历史对照和证据追溯锚点。"),
        )
    return (
        "报告复核已留痕，样本库已有记录，下一步请确认历史对照结论。",
        _next_step("compare_history", "确认历史对照", "请记录是否发现聚集、传播或关联线索。"),
    )


def _history_confirmation_guidance(comparison_outcome: str) -> tuple[str, dict[str, str]]:
    if comparison_outcome == "needs_escalation":
        return (
            "已记录“需升级复核或处置”，正式归档导出仍保持阻塞。",
            _next_step("compare_history", "升级处置后重新确认历史对照", "完成会商、流调或复核处置后，再记录最终历史对照结论。"),
        )
    if comparison_outcome == "signal_found":
        return (
            "历史对照已留痕，并记录发现聚集或传播线索；下一步请从报告页完成正式导出归档。",
            _next_step("archive_export", "正式导出归档", "导出记录会作为任务交付闭环的最后一步。"),
        )
    return (
        "历史对照已留痕，未发现聚集或传播线索；下一步请从报告页完成正式导出归档。",
        _next_step("archive_export", "正式导出归档", "导出记录会作为任务交付闭环的最后一步。"),
    )


def _database_import_guidance(result: dict) -> dict:
    imported_count = int(result.get("imported_count") or 0)
    skipped_count = int(result.get("skipped_count") or 0)
    if imported_count > 0:
        message = f"已导入 {imported_count} 个样本"
        if skipped_count:
            message += f"，跳过 {skipped_count} 个"
        message += "；下一步请确认历史对照结论。"
        return {
            **result,
            "message": message,
            "next_action": _next_step("compare_history", "确认历史对照", "请筛查同地区、同型别或同物种历史样本，并记录是否发现聚集或传播线索。"),
        }
    message = "没有成功导入样本"
    if skipped_count:
        message += f"，跳过 {skipped_count} 个"
    message += "；请先检查输出文件和可导入样本清单。"
    return {
        **result,
        "message": message,
        "next_action": _next_step("open_output", "检查输出文件", "确认 final.fasta、报告目录或 meta bin fasta 是否已生成。"),
    }


def _event_is_after_current_run(event: dict, finished_at: object) -> bool:
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


def _successful_closure_action_labels(events: list[dict], *, task_owner: str, finished_at: object = "") -> set[str]:
    completed: set[str] = set()
    for item in events:
        if str(item.get("outcome") or "").strip() != "success":
            continue
        if not _event_is_after_current_run(item, finished_at):
            continue
        action = str(item.get("action") or "").strip()
        if action == "确认历史对照":
            if _audit_request_value(item, "comparison_outcome") not in COMPLETED_HISTORY_COMPARISON_OUTCOMES:
                continue
        if action == "确认报告复核":
            reviewer = str(item.get("username") or "").strip()
            if not task_owner or not reviewer or reviewer == task_owner:
                continue
            if _audit_request_value(item, "review_outcome") not in COMPLETED_REPORT_REVIEW_OUTCOMES:
                continue
        completed.add(action)
    return completed


def _audit_request_value(event: dict, key: str) -> str:
    raw = str(event.get("request_summary") or "").strip()
    if not raw:
        return ""
    try:
        summary = json.loads(raw)
    except (TypeError, ValueError):
        return ""
    payload = summary.get("json") if isinstance(summary, dict) else None
    return str(payload.get(key) or "").strip() if isinstance(payload, dict) else ""


@dataclass(frozen=True)
class TaskLifecycleService:
    store: PortalStore
    task_manager: AnalysisTaskManager
    sample_manager: SampleLibraryManager
    ensure_can_view_task: Callable[[dict], None]
    ensure_can_modify_task: Callable[[dict], None]
    ensure_can_control_task: Callable[[dict], None]
    open_in_file_manager: Callable[[Path], None]

    def get_visible(self, task_id: str, *, log_lines: int = 0) -> dict:
        task = self.task_manager.get_task(task_id, log_lines=log_lines, owner=None)
        self.ensure_can_view_task(task)
        return task

    def task_directory(self, task_id: str) -> Path:
        return self.task_manager.task_root / str(task_id)

    def delete(self, task_id: str) -> None:
        task = self.task_manager.get_task(task_id, log_lines=0, owner=None)
        self.ensure_can_modify_task(task)
        self.task_manager.delete_task(task_id, owner=None)

    def pause(self, task_id: str) -> dict:
        return self._control(task_id, self.task_manager.pause_task)

    def resume(self, task_id: str) -> dict:
        return self._control(task_id, self.task_manager.resume_task)

    def stop(self, task_id: str) -> dict:
        return self._control(task_id, self.task_manager.stop_task)

    def open_output(self, task_id: str) -> Path:
        task = self.get_visible(task_id)
        output_dir = str((task.get("params") or {}).get("output_dir", "") or "").strip()
        if not output_dir:
            raise ValidationError("当前任务没有输出目录")
        target = Path(output_dir).expanduser().resolve()
        if not target.exists():
            raise ValidationError(f"输出目录不存在: {target}")
        self.open_in_file_manager(target)
        return target

    def import_to_database(self, task_id: str, payload: dict, *, identity: UserIdentity) -> dict:
        task = self.get_visible(task_id)
        requested_scope = str(payload.get("library_scope") or "").strip()
        library_scope = (
            requested_scope
            if identity.role == "admin" and requested_scope in {"main", "personal"}
            else ("main" if identity.role == "admin" else "personal")
        )
        if str((task.get("params") or {}).get("method") or "").strip().lower() == "meta":
            selected_bins = payload.get("selected_bins") or []
            result = self.sample_manager.import_meta_task_bins(
                task,
                selected_bins=selected_bins if isinstance(selected_bins, list) else [],
                library_scope=library_scope,
            )
            return _database_import_guidance(result)
        return _database_import_guidance(self.sample_manager.import_task_samples(task, library_scope=library_scope))

    def preview_database_import(self, task_id: str) -> dict:
        task = self.get_visible(task_id)
        if str((task.get("params") or {}).get("method") or "").strip().lower() != "meta":
            raise ValidationError("当前只有 meta 任务需要预览 bin 导入清单")
        return self.sample_manager.preview_meta_task_bins(task)

    def confirm_closure_action(self, task_id: str, action_id: str, payload: dict, *, identity: UserIdentity) -> dict:
        action_id = str(action_id or "").strip()
        if action_id not in CLOSURE_ACTION_LABELS:
            raise ValidationError("不支持的闭环动作")
        task = self.task_manager.get_task(task_id, log_lines=0, owner=None)
        self.ensure_can_modify_task(task)
        if str(task.get("status") or "").strip().upper() != "SUCCEEDED":
            raise ValidationError("仅已完成任务可以确认闭环动作")
        if action_id == "review_report":
            task_owner = str(task.get("owner") or "").strip()
            if not task_owner:
                raise ValidationError("任务缺少分析执行人，无法完成双人复核签核")
            if task_owner == identity.username:
                raise ValidationError("报告复核必须由非任务创建人完成")
            review_outcome = str(payload.get("review_outcome") or "").strip()
            if review_outcome not in REPORT_REVIEW_OUTCOMES:
                raise ValidationError("请选择报告复核结论")
        else:
            review_outcome = ""
        events = self.store.list_audit_logs_for_target("task", task_id, limit=50)
        completed_actions = _successful_closure_action_labels(
            events,
            task_owner=str(task.get("owner") or "").strip(),
            finished_at=task.get("finished_at"),
        )
        if action_id == "compare_history":
            if self.store.count_sample_library_records_by_task_id(task_id) <= 0:
                raise ValidationError("请先将任务结果导入样本库，再确认历史对照")
            if "确认报告复核" not in completed_actions:
                raise ValidationError("请先确认报告复核完成")
            comparison_outcome = str(payload.get("comparison_outcome") or "").strip()
            if comparison_outcome not in HISTORY_COMPARISON_OUTCOMES:
                raise ValidationError("请选择历史对照结论")
        else:
            comparison_outcome = ""
        note = str(payload.get("note") or "").strip()
        if not note:
            raise ValidationError("请填写本次闭环确认结论")
        if action_id == "review_report":
            message, next_action = _review_confirmation_guidance(
                review_outcome,
                imported_count=self.store.count_sample_library_records_by_task_id(task_id),
            )
        else:
            message, next_action = _history_confirmation_guidance(comparison_outcome)
        return {
            "status": "confirmed",
            "action_id": action_id,
            "label": CLOSURE_ACTION_LABELS[action_id],
            "message": message,
            "next_action": next_action,
            "note": note,
            "review_outcome": review_outcome,
            "review_outcome_label": REPORT_REVIEW_OUTCOMES.get(review_outcome, ""),
            "comparison_outcome": comparison_outcome,
            "comparison_outcome_label": HISTORY_COMPARISON_OUTCOMES.get(comparison_outcome, ""),
            "task_id": task_id,
        }

    def record_report_export(self, task_id: str, payload: dict) -> dict:
        task = self.task_manager.get_task(task_id, log_lines=0, owner=None)
        self.ensure_can_modify_task(task)
        if str(task.get("status") or "").strip().upper() != "SUCCEEDED":
            raise ValidationError("仅已完成任务可以记录正式报告导出")
        export_format = str(payload.get("format") or "").strip().lower()
        if export_format not in REPORT_EXPORT_FORMATS:
            raise ValidationError("不支持的报告导出格式")
        events = self.store.list_audit_logs_for_target("task", task_id, limit=50)
        completed_actions = _successful_closure_action_labels(
            events,
            task_owner=str(task.get("owner") or "").strip(),
            finished_at=task.get("finished_at"),
        )
        if "确认报告复核" not in completed_actions or "确认历史对照" not in completed_actions:
            raise ValidationError("请先完成报告复核、样本入库和历史对照，再执行正式导出归档")
        if self.store.count_sample_library_records_by_task_id(task_id) <= 0:
            raise ValidationError("请先将任务结果导入样本库，再执行正式导出归档")
        return {
            "status": "recorded",
            "action_id": "archive_export",
            "label": "导出归档完成",
            "message": "正式报告导出已留痕，任务已具备交付闭环证据。",
            "next_action": _next_step("review_trace", "查看闭环留痕", "可回到任务队列核查报告复核、样本入库、历史对照和导出记录。"),
            "format": export_format,
            "task_id": task_id,
        }

    def record_failure_disposition(self, task_id: str, payload: dict, *, identity: UserIdentity) -> dict:
        task = self.task_manager.get_task(task_id, log_lines=0, owner=None)
        self.ensure_can_modify_task(task)
        task_status = str(task.get("status") or "").strip().upper()
        decision = str(payload.get("decision") or "").strip()
        expected_decision = {"FAILED": "accept_failure", "STOPPED": "accept_stop"}.get(task_status)
        if not expected_decision:
            raise ValidationError("仅失败或已停止任务可以确认异常归档")
        if decision != expected_decision:
            raise ValidationError("不支持的异常处置决定")
        note = str(payload.get("note") or "").strip()
        if not note:
            raise ValidationError("请填写异常原因及不再继续运行的处置依据")
        diagnosis = task.get("failure_diagnosis") if isinstance(task.get("failure_diagnosis"), dict) else {}
        label = "失败归档完成" if task_status == "FAILED" else "停止归档完成"
        return {
            "status": "recorded",
            "decision": decision,
            "label": label,
            "message": f"{label}，异常原因和处置依据已留痕。",
            "next_action": _next_step("review_trace", "查看异常留痕", "可回到任务队列核查失败/停止原因、影响判断和处置依据。"),
            "note": note,
            "diagnosis_category": str(payload.get("diagnosis_category") or diagnosis.get("category") or "").strip(),
            "diagnosis_label": str(payload.get("diagnosis_label") or diagnosis.get("label") or "").strip(),
            "operator": identity.username,
            "task_id": task_id,
        }

    def _control(self, task_id: str, operation: Callable[..., dict]) -> dict:
        task = self.task_manager.get_task(task_id, log_lines=0, owner=None)
        self.ensure_can_control_task(task)
        return operation(task_id, owner=None, max_concurrent_tasks=self._max_concurrent_tasks())

    def _max_concurrent_tasks(self) -> int:
        return int(self.store.get_setting("max_concurrent_tasks", "2") or "2")
