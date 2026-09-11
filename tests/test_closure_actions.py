from __future__ import annotations

import unittest
from unittest.mock import Mock

from flask import Flask

from bac_analysis_portal.audit_hooks import _meta
from bac_analysis_portal.identity import UserIdentity
from bac_analysis_portal.task_lifecycle_service import TaskLifecycleService
from bac_analysis_portal.task_manager import ValidationError


class ClosureActionTests(unittest.TestCase):
    def _service(self, *, imported_count: int = 1, audit_events: list[dict] | None = None, status: str = "SUCCEEDED") -> TaskLifecycleService:
        store = Mock()
        store.count_sample_library_records_by_task_id.return_value = imported_count
        store.list_audit_logs_for_target.return_value = audit_events or []
        task_manager = Mock()
        task_manager.get_task.return_value = {"id": "task-1", "status": status, "owner": "analyst", "params": {}}
        return TaskLifecycleService(
            store=store,
            task_manager=task_manager,
            sample_manager=Mock(),
            ensure_can_view_task=Mock(),
            ensure_can_modify_task=Mock(),
            ensure_can_control_task=Mock(),
            open_in_file_manager=Mock(),
        )

    def _identity(self, username: str = "reviewer") -> UserIdentity:
        return UserIdentity(username=username, role="group_admin", group_name="cdc")

    def test_report_review_can_be_confirmed_for_completed_task(self) -> None:
        result = self._service().confirm_closure_action(
            "task-1",
            "review_report",
            {"review_outcome": "approved", "note": "已复核"},
            identity=self._identity(),
        )
        self.assertEqual(result["status"], "confirmed")
        self.assertEqual(result["action_id"], "review_report")
        self.assertEqual(result["note"], "已复核")
        self.assertEqual(result["review_outcome_label"], "复核通过")
        self.assertEqual(result["next_action"]["id"], "compare_history")
        self.assertIn("下一步", result["message"])

    def test_database_import_returns_history_comparison_guidance(self) -> None:
        service = self._service()
        service.sample_manager.import_task_samples.return_value = {
            "status": "ok",
            "imported_count": 1,
            "skipped_count": 0,
            "items": [{"sample_name": "sample-a"}],
            "skipped": [],
        }
        result = service.import_to_database("task-1", {}, identity=self._identity())
        self.assertEqual(result["next_action"]["id"], "compare_history")
        self.assertIn("下一步请确认历史对照结论", result["message"])

    def test_database_import_without_records_returns_output_check_guidance(self) -> None:
        service = self._service()
        service.sample_manager.import_task_samples.return_value = {
            "status": "ok",
            "imported_count": 0,
            "skipped_count": 1,
            "items": [],
            "skipped": [{"sample_name": "sample-a", "reason": "缺少 final.fasta"}],
        }
        result = service.import_to_database("task-1", {}, identity=self._identity())
        self.assertEqual(result["next_action"]["id"], "open_output")
        self.assertIn("检查输出文件", result["next_action"]["label"])

    def test_report_review_requires_structured_outcome(self) -> None:
        with self.assertRaises(ValidationError):
            self._service().confirm_closure_action("task-1", "review_report", {"note": "已复核"}, identity=self._identity())

    def test_reanalysis_review_does_not_unlock_history_confirmation(self) -> None:
        with self.assertRaises(ValidationError):
            self._service(
                imported_count=1,
                audit_events=[
                    {
                        "action": "确认报告复核",
                        "username": "reviewer",
                        "outcome": "success",
                        "request_summary": '{"json": {"review_outcome": "needs_reanalysis", "note": "质控不达标"}}',
                    }
                ],
            ).confirm_closure_action(
                "task-1",
                "compare_history",
                {"comparison_outcome": "no_signal", "note": "未发现聚集"},
                identity=self._identity(),
            )

    def test_task_creator_cannot_review_own_report(self) -> None:
        with self.assertRaises(ValidationError):
            self._service().confirm_closure_action(
                "task-1",
                "review_report",
                {"note": "本人复核"},
                identity=self._identity("analyst"),
            )

    def test_history_confirmation_requires_import_and_report_review(self) -> None:
        with self.assertRaises(ValidationError):
            self._service(imported_count=0).confirm_closure_action("task-1", "compare_history", {"note": "未发现聚集"}, identity=self._identity())
        with self.assertRaises(ValidationError):
            self._service(imported_count=1).confirm_closure_action("task-1", "compare_history", {"note": "未发现聚集"}, identity=self._identity())
        result = self._service(
            imported_count=1,
            audit_events=[
                {
                    "action": "确认报告复核",
                    "username": "reviewer",
                    "outcome": "success",
                    "request_summary": '{"json": {"review_outcome": "approved", "note": "已复核"}}',
                }
            ],
        ).confirm_closure_action(
            "task-1",
            "compare_history",
            {"comparison_outcome": "no_signal", "note": "未发现聚集"},
            identity=self._identity(),
        )
        self.assertEqual(result["action_id"], "compare_history")
        self.assertEqual(result["comparison_outcome_label"], "未发现聚集或传播线索")
        self.assertEqual(result["next_action"]["id"], "archive_export")

    def test_history_confirmation_requires_structured_outcome(self) -> None:
        service = self._service(
            imported_count=1,
            audit_events=[
                {
                    "action": "确认报告复核",
                    "username": "reviewer",
                    "outcome": "success",
                    "request_summary": '{"json": {"review_outcome": "approved", "note": "已复核"}}',
                }
            ],
        )
        with self.assertRaises(ValidationError):
            service.confirm_closure_action("task-1", "compare_history", {"note": "已查看历史样本"}, identity=self._identity())

    def test_legacy_unstructured_review_audit_does_not_unlock_history(self) -> None:
        with self.assertRaises(ValidationError):
            self._service(
                imported_count=1,
                audit_events=[{"action": "确认报告复核", "username": "reviewer", "outcome": "success"}],
            ).confirm_closure_action(
                "task-1",
                "compare_history",
                {"comparison_outcome": "no_signal", "note": "未发现聚集"},
                identity=self._identity(),
            )

    def test_failed_review_audit_does_not_unlock_history_confirmation(self) -> None:
        with self.assertRaises(ValidationError):
            self._service(
                imported_count=1,
                audit_events=[{"action": "确认报告复核", "outcome": "failed"}],
            ).confirm_closure_action("task-1", "compare_history", {"note": "未发现聚集"}, identity=self._identity())

    def test_creator_review_audit_does_not_unlock_downstream_actions(self) -> None:
        creator_review = [{"action": "确认报告复核", "username": "analyst", "outcome": "success"}]
        with self.assertRaises(ValidationError):
            self._service(imported_count=1, audit_events=creator_review).confirm_closure_action(
                "task-1",
                "compare_history",
                {"note": "未发现聚集"},
                identity=self._identity(),
            )
        with self.assertRaises(ValidationError):
            self._service(
                imported_count=1,
                audit_events=creator_review + [{"action": "确认历史对照", "username": "reviewer", "outcome": "success"}],
            ).record_report_export("task-1", {"format": "pdf"})

    def test_closure_confirmation_requires_auditable_conclusion(self) -> None:
        with self.assertRaises(ValidationError):
            self._service().confirm_closure_action("task-1", "review_report", {"note": " "}, identity=self._identity())

    def test_archive_export_cannot_be_manually_confirmed(self) -> None:
        with self.assertRaises(ValidationError):
            self._service(
                imported_count=1,
                audit_events=[
                    {"action": "确认报告复核", "outcome": "success"},
                    {"action": "确认历史对照", "outcome": "success"},
                ],
            ).confirm_closure_action("task-1", "archive_export", {}, identity=self._identity())

    def test_report_export_requires_completed_prior_steps(self) -> None:
        with self.assertRaises(ValidationError):
            self._service(imported_count=1).record_report_export("task-1", {"format": "pdf"})
        with self.assertRaises(ValidationError):
            self._service(
                imported_count=0,
                audit_events=[
                    {"action": "确认报告复核", "outcome": "success"},
                    {"action": "确认历史对照", "outcome": "success"},
                ],
            ).record_report_export("task-1", {"format": "pdf"})
        result = self._service(
            imported_count=1,
            audit_events=[
                {
                    "action": "确认报告复核",
                    "username": "reviewer",
                    "outcome": "success",
                    "request_summary": '{"json": {"review_outcome": "approved", "note": "已复核"}}',
                },
                {
                    "action": "确认历史对照",
                    "username": "reviewer",
                    "outcome": "success",
                    "request_summary": '{"json": {"comparison_outcome": "no_signal", "note": "未发现聚集"}}',
                },
            ],
        ).record_report_export("task-1", {"format": "word"})
        self.assertEqual(result["status"], "recorded")
        self.assertEqual(result["format"], "word")
        self.assertEqual(result["next_action"]["id"], "review_trace")
        self.assertIn("交付闭环", result["message"])

    def test_review_that_needs_reanalysis_returns_blocking_guidance(self) -> None:
        result = self._service().confirm_closure_action(
            "task-1",
            "review_report",
            {"review_outcome": "needs_reanalysis", "note": "覆盖度不足"},
            identity=self._identity(),
        )
        self.assertEqual(result["next_action"]["id"], "decide_rebuild")
        self.assertIn("仍保持阻塞", result["message"])

    def test_history_escalation_returns_blocking_guidance(self) -> None:
        result = self._service(
            imported_count=1,
            audit_events=[
                {
                    "action": "确认报告复核",
                    "username": "reviewer",
                    "outcome": "success",
                    "request_summary": '{"json": {"review_outcome": "approved", "note": "已复核"}}',
                }
            ],
        ).confirm_closure_action(
            "task-1",
            "compare_history",
            {"comparison_outcome": "needs_escalation", "note": "发现疑似聚集"},
            identity=self._identity(),
        )
        self.assertEqual(result["next_action"]["id"], "compare_history")
        self.assertIn("正式归档导出仍保持阻塞", result["message"])

    def test_report_export_rejects_unknown_format(self) -> None:
        with self.assertRaises(ValidationError):
            self._service().record_report_export("task-1", {"format": "zip"})

    def test_report_export_is_blocked_while_history_comparison_needs_escalation(self) -> None:
        events = [
            {
                "action": "确认报告复核",
                "username": "reviewer",
                "outcome": "success",
                "request_summary": '{"json": {"review_outcome": "approved", "note": "已复核"}}',
            },
            {
                "action": "确认历史对照",
                "username": "reviewer",
                "outcome": "success",
                "request_summary": '{"json": {"comparison_outcome": "needs_escalation", "note": "发现疑似聚集"}}',
            },
        ]
        with self.assertRaises(ValidationError):
            self._service(imported_count=1, audit_events=events).record_report_export("task-1", {"format": "pdf"})

    def test_failed_task_can_be_archived_only_with_explicit_reason(self) -> None:
        with self.assertRaises(ValidationError):
            self._service(status="FAILED").record_failure_disposition(
                "task-1",
                {"decision": "accept_failure", "note": ""},
                identity=self._identity(),
            )
        result = self._service(status="FAILED").record_failure_disposition(
            "task-1",
            {"decision": "accept_failure", "note": "输入样本损坏，无法补采"},
            identity=self._identity(),
        )
        self.assertEqual(result["status"], "recorded")
        self.assertEqual(result["decision"], "accept_failure")
        self.assertEqual(result["next_action"]["id"], "review_trace")

    def test_failure_disposition_returns_server_side_diagnosis_for_audit_trace(self) -> None:
        service = self._service(status="FAILED")
        service.task_manager.get_task.return_value["failure_diagnosis"] = {
            "category": "dependency_missing",
            "label": "运行环境依赖缺失",
        }
        result = service.record_failure_disposition(
            "task-1",
            {"decision": "accept_failure", "note": "已确认无需重跑"},
            identity=self._identity(),
        )
        self.assertEqual(result["diagnosis_category"], "dependency_missing")
        self.assertEqual(result["diagnosis_label"], "运行环境依赖缺失")

    def test_stopped_task_can_be_archived_with_explicit_reason(self) -> None:
        result = self._service(status="STOPPED").record_failure_disposition(
            "task-1",
            {"decision": "accept_stop", "note": "用户确认取消，样本需重新送检"},
            identity=self._identity(),
        )
        self.assertEqual(result["status"], "recorded")
        self.assertEqual(result["decision"], "accept_stop")
        self.assertEqual(result["label"], "停止归档完成")

    def test_successful_task_cannot_be_archived_as_failure(self) -> None:
        with self.assertRaises(ValidationError):
            self._service().record_failure_disposition(
                "task-1",
                {"decision": "accept_failure", "note": "不适用"},
                identity=self._identity(),
            )

    def test_audit_hook_names_actual_report_export_as_archive_action(self) -> None:
        app = Flask(__name__)
        with app.test_request_context("/api/tasks/task-1/report-exports", method="POST"):
            module, action, target_type, target_id = _meta(
                "/api/tasks/task-1/report-exports",
                "POST",
            )
        self.assertEqual((module, action, target_type, target_id), ("任务", "确认导出归档", "task", "task-1"))

    def test_audit_hook_does_not_name_manual_archive_attempt_as_export(self) -> None:
        app = Flask(__name__)
        path = "/api/tasks/task-1/closure-actions/archive_export"
        with app.test_request_context(path, method="POST"):
            module, action, target_type, target_id = _meta(path, "POST")
        self.assertEqual((module, action, target_type, target_id), ("任务", "确认闭环动作", "task", "task-1"))

    def test_audit_hook_names_failure_disposition(self) -> None:
        app = Flask(__name__)
        path = "/api/tasks/task-1/failure-disposition"
        with app.test_request_context(path, method="POST"):
            module, action, target_type, target_id = _meta(path, "POST")
        self.assertEqual((module, action, target_type, target_id), ("任务", "确认失败归档", "task", "task-1"))

    def test_audit_hook_names_stopped_disposition(self) -> None:
        app = Flask(__name__)
        path = "/api/tasks/task-1/failure-disposition"
        with app.test_request_context(path, method="POST", json={"decision": "accept_stop"}):
            module, action, target_type, target_id = _meta(path, "POST")
        self.assertEqual((module, action, target_type, target_id), ("任务", "确认停止归档", "task", "task-1"))


if __name__ == "__main__":
    unittest.main()
