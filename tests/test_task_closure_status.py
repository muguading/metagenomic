from __future__ import annotations

import unittest

from bac_analysis_portal.task_closure_status import build_task_closure_status


class TaskClosureStatusTests(unittest.TestCase):
    def test_completed_task_without_database_import_points_to_import_step(self) -> None:
        closure = build_task_closure_status(
            {"id": "task-1", "status": "SUCCEEDED", "owner": "analyst"},
            result_exists=True,
            imported_sample_count=0,
        )
        self.assertEqual(closure["state"], "needs_report_review")
        self.assertEqual(closure["next_action"]["id"], "review_report")
        self.assertEqual(closure["confirmable_actions"][0]["id"], "review_report")

    def test_completed_task_with_imported_samples_points_to_history_comparison(self) -> None:
        closure = build_task_closure_status(
            {"id": "task-1", "status": "SUCCEEDED", "owner": "analyst"},
            result_exists=True,
            imported_sample_count=2,
            audit_events=[
                {
                    "action": "确认报告复核",
                    "username": "admin",
                    "outcome": "success",
                    "created_at": "2026-06-08T08:00:00",
                    "request_summary": '{"json": {"review_outcome": "approved", "note": "分型与质控已复核"}}',
                },
            ],
        )
        self.assertEqual(closure["state"], "needs_history_compare")
        self.assertEqual(closure["imported_sample_count"], 2)
        self.assertEqual(closure["next_action"]["id"], "compare_history")
        self.assertEqual(closure["recent_events"][0]["operator"], "admin")
        self.assertEqual(closure["recent_events"][0]["detail"], "复核通过；分型与质控已复核")
        self.assertTrue(closure["responsibility"]["separation_verified"])

    def test_reanalysis_review_keeps_task_in_report_review_state(self) -> None:
        closure = build_task_closure_status(
            {"id": "task-1", "status": "SUCCEEDED", "owner": "analyst"},
            result_exists=True,
            imported_sample_count=1,
            audit_events=[
                {
                    "action": "确认报告复核",
                    "username": "admin",
                    "outcome": "success",
                    "request_summary": '{"json": {"review_outcome": "needs_reanalysis", "note": "质控不达标"}}',
                },
            ],
        )
        self.assertEqual(closure["state"], "needs_report_review")
        self.assertEqual(closure["recent_events"][0]["detail"], "需重分析或补充复核；质控不达标")

    def test_task_is_closed_only_after_all_confirmations(self) -> None:
        closure = build_task_closure_status(
            {"id": "task-1", "status": "SUCCEEDED", "owner": "analyst"},
            result_exists=True,
            imported_sample_count=1,
            audit_events=[
                {"action": "确认导出归档", "username": "admin", "outcome": "success"},
                {
                    "action": "确认历史对照",
                    "username": "admin",
                    "outcome": "success",
                    "request_summary": '{"json": {"comparison_outcome": "no_signal", "note": "未发现聚集"}}',
                },
                {
                    "action": "确认报告复核",
                    "username": "admin",
                    "outcome": "success",
                    "request_summary": '{"json": {"review_outcome": "approved", "note": "已复核"}}',
                },
            ],
        )
        self.assertEqual(closure["state"], "closed")
        self.assertTrue(all(step["state"] == "done" for step in closure["steps"]))
        self.assertEqual([item["id"] for item in closure["delivery_evidence"]], ["report_review", "database_import", "history_compare", "archive_export"])
        self.assertEqual(closure["delivery_evidence"][1]["value"], "1 份样本")

    def test_archive_trace_exposes_actual_export_format(self) -> None:
        closure = build_task_closure_status(
            {"id": "task-1", "status": "SUCCEEDED", "owner": "analyst"},
            result_exists=True,
            imported_sample_count=1,
            audit_events=[
                {
                    "action": "确认导出归档",
                    "username": "admin",
                    "outcome": "success",
                    "request_summary": '{"json": {"format": "pdf"}}',
                },
                {
                    "action": "确认历史对照",
                    "username": "admin",
                    "outcome": "success",
                    "request_summary": '{"json": {"comparison_outcome": "no_signal", "note": "未发现聚集"}}',
                },
                {
                    "action": "确认报告复核",
                    "username": "admin",
                    "outcome": "success",
                    "request_summary": '{"json": {"review_outcome": "approved", "note": "已复核"}}',
                },
            ],
        )
        self.assertEqual(closure["recent_events"][0]["detail"], "正式导出格式：PDF")
        archive_evidence = [item for item in closure["delivery_evidence"] if item["id"] == "archive_export"][0]
        self.assertEqual(archive_evidence["detail"], "正式导出格式：PDF")

    def test_archive_step_requires_actual_report_export_instead_of_manual_confirmation(self) -> None:
        closure = build_task_closure_status(
            {"id": "task-1", "status": "SUCCEEDED", "owner": "analyst"},
            result_exists=True,
            imported_sample_count=1,
            audit_events=[
                {
                    "action": "确认历史对照",
                    "username": "admin",
                    "outcome": "success",
                    "request_summary": '{"json": {"comparison_outcome": "no_signal", "note": "未发现聚集"}}',
                },
                {
                    "action": "确认报告复核",
                    "username": "admin",
                    "outcome": "success",
                    "request_summary": '{"json": {"review_outcome": "approved", "note": "已复核"}}',
                },
            ],
        )
        self.assertEqual(closure["state"], "needs_archive_export")
        self.assertEqual(closure["confirmable_actions"], [])
        self.assertIn("报告页完成正式导出", closure["summary"])

    def test_creator_review_audit_does_not_satisfy_separation_of_duties(self) -> None:
        closure = build_task_closure_status(
            {"id": "task-1", "status": "SUCCEEDED", "owner": "analyst"},
            result_exists=True,
            imported_sample_count=1,
            audit_events=[
                {"action": "确认报告复核", "username": "analyst", "outcome": "success"},
            ],
        )
        self.assertEqual(closure["state"], "needs_report_review")
        self.assertFalse(closure["responsibility"]["separation_verified"])
        self.assertEqual(closure["responsibility"]["reviewer"], "")

    def test_anonymous_review_audit_does_not_satisfy_separation_of_duties(self) -> None:
        closure = build_task_closure_status(
            {"id": "task-1", "status": "SUCCEEDED", "owner": "analyst"},
            result_exists=True,
            imported_sample_count=1,
            audit_events=[{"action": "确认报告复核", "username": "", "outcome": "success"}],
        )
        self.assertEqual(closure["state"], "needs_report_review")
        self.assertFalse(closure["responsibility"]["separation_verified"])

    def test_missing_eligible_reviewer_is_an_explicit_blocker(self) -> None:
        closure = build_task_closure_status(
            {"id": "task-1", "status": "SUCCEEDED", "owner": "admin"},
            result_exists=True,
            imported_sample_count=0,
            eligible_reviewers=[],
        )
        self.assertEqual(closure["state"], "reviewer_unavailable")
        self.assertEqual(closure["next_action"]["id"], "configure_reviewer")
        self.assertIn("admin-users-section", closure["next_action"]["href"])
        self.assertEqual(closure["steps"][1]["state"], "blocked")

    def test_missing_report_exposes_specific_diagnostic_reason(self) -> None:
        closure = build_task_closure_status(
            {"id": "task-1", "status": "SUCCEEDED", "owner": "analyst"},
            result_exists=False,
            report_issue="输出目录存在，但其中没有识别到结果文件",
        )
        self.assertEqual(closure["state"], "report_missing")
        self.assertEqual(closure["summary"], "输出目录存在，但其中没有识别到结果文件")

    def test_failed_task_blocks_closure_before_report_review(self) -> None:
        closure = build_task_closure_status(
            {"id": "task-1", "status": "FAILED", "owner": "analyst"},
            result_exists=False,
            imported_sample_count=0,
        )
        self.assertEqual(closure["state"], "analysis_failed")
        self.assertEqual(closure["steps"][0]["state"], "blocked")
        self.assertEqual(closure["confirmable_actions"][0]["id"], "accept_failure")

    def test_failed_task_surfaces_structured_diagnosis_in_closure_summary(self) -> None:
        closure = build_task_closure_status(
            {
                "id": "task-1",
                "status": "FAILED",
                "owner": "analyst",
                "failure_diagnosis": {
                    "label": "运行环境依赖缺失",
                    "summary": "分析程序缺少运行依赖，当前失败不能用于判断样本质量。",
                },
            },
        )
        self.assertIn("运行环境依赖缺失", closure["summary"])
        self.assertIn("不能用于判断样本质量", closure["summary"])

    def test_failed_task_becomes_resolved_only_after_failure_archive_audit(self) -> None:
        closure = build_task_closure_status(
            {"id": "task-1", "status": "FAILED", "owner": "analyst"},
            audit_events=[
                {
                    "action": "确认失败归档",
                    "username": "admin",
                    "outcome": "success",
                    "request_summary": '{"json": {"decision": "accept_failure", "note": "样本无法补采"}}',
                }
            ],
        )
        self.assertEqual(closure["state"], "failure_archived")
        self.assertEqual(closure["recent_events"][0]["detail"], "样本无法补采")
        self.assertEqual(closure["confirmable_actions"], [])

    def test_failure_archive_trace_preserves_diagnosis_and_human_conclusion(self) -> None:
        closure = build_task_closure_status(
            {"id": "task-1", "status": "FAILED", "owner": "analyst"},
            audit_events=[
                {
                    "action": "确认失败归档",
                    "username": "admin",
                    "outcome": "success",
                    "request_summary": '{"json": {"diagnosis_label": "运行环境依赖缺失", "note": "已转系统管理员处理"}}',
                }
            ],
        )
        self.assertEqual(closure["recent_events"][0]["detail"], "运行环境依赖缺失；已转系统管理员处理")

    def test_history_comparison_trace_preserves_structured_outcome(self) -> None:
        closure = build_task_closure_status(
            {"id": "task-1", "status": "SUCCEEDED", "owner": "analyst"},
            result_exists=True,
            imported_sample_count=1,
            audit_events=[
                {
                    "action": "确认历史对照",
                    "username": "admin",
                    "outcome": "success",
                    "request_summary": '{"json": {"comparison_outcome": "signal_found", "note": "同区同期发现近似样本"}}',
                }
            ],
        )
        self.assertEqual(closure["recent_events"][0]["detail"], "发现聚集或传播线索；同区同期发现近似样本")

    def test_history_escalation_blocks_archive_until_comparison_is_reconfirmed(self) -> None:
        closure = build_task_closure_status(
            {"id": "task-1", "status": "SUCCEEDED", "owner": "analyst"},
            result_exists=True,
            imported_sample_count=1,
            audit_events=[
                {
                    "action": "确认历史对照",
                    "username": "admin",
                    "outcome": "success",
                    "request_summary": '{"json": {"comparison_outcome": "needs_escalation", "note": "疑似聚集需复核"}}',
                },
                {
                    "action": "确认报告复核",
                    "username": "reviewer",
                    "outcome": "success",
                    "request_summary": '{"json": {"review_outcome": "approved", "note": "已复核"}}',
                },
            ],
        )
        self.assertEqual(closure["state"], "needs_history_compare")
        self.assertEqual(closure["label"], "待升级处置")
        self.assertIn("重新确认历史对照", closure["summary"])

    def test_old_failure_archive_does_not_apply_to_a_new_failed_run(self) -> None:
        closure = build_task_closure_status(
            {"id": "task-1", "status": "FAILED", "owner": "analyst", "finished_at": "2026-06-08T10:00:00+00:00"},
            audit_events=[
                {
                    "action": "确认失败归档",
                    "username": "admin",
                    "outcome": "success",
                    "created_at": "2026-06-08T09:00:00+00:00",
                }
            ],
        )
        self.assertEqual(closure["state"], "analysis_failed")

    def test_old_review_does_not_unlock_a_new_successful_run(self) -> None:
        closure = build_task_closure_status(
            {"id": "task-1", "status": "SUCCEEDED", "owner": "analyst", "finished_at": "2026-06-08T10:00:00+00:00"},
            result_exists=True,
            imported_sample_count=1,
            audit_events=[
                {
                    "action": "确认报告复核",
                    "username": "reviewer",
                    "outcome": "success",
                    "created_at": "2026-06-08T09:00:00+00:00",
                }
            ],
        )
        self.assertEqual(closure["state"], "needs_report_review")

    def test_stopped_task_can_be_archived_with_audit_trace(self) -> None:
        closure = build_task_closure_status(
            {"id": "task-1", "status": "STOPPED", "owner": "analyst"},
            audit_events=[
                {
                    "action": "确认停止归档",
                    "username": "admin",
                    "outcome": "success",
                    "request_summary": '{"json": {"decision": "accept_stop", "note": "用户取消并重新送检"}}',
                }
            ],
        )
        self.assertEqual(closure["state"], "stopped_archived")
        self.assertEqual(closure["recent_events"][0]["detail"], "用户取消并重新送检")

    def test_stopped_task_exposes_archive_action_before_disposition(self) -> None:
        closure = build_task_closure_status(
            {"id": "task-1", "status": "STOPPED", "owner": "analyst"},
        )
        self.assertEqual(closure["state"], "analysis_stopped")
        self.assertEqual(closure["confirmable_actions"][0]["id"], "accept_stop")


if __name__ == "__main__":
    unittest.main()
