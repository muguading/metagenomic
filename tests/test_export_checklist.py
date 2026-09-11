from __future__ import annotations

import unittest

from bac_analysis_portal.export_checklist import build_export_checklist


class ExportChecklistTests(unittest.TestCase):
    def test_priority_report_requires_manual_public_health_review_before_export(self) -> None:
        checklist = build_export_checklist(
            {
                "task": {"id": "task-1", "name": "Priority sample"},
                "sections": {
                    "workflow_closure": {"risk_level": "高"},
                    "public_health_support": {"status": "matched"},
                    "resistance_virulence": {
                        "resistance_elements": {"rows": [["blaKPC"]]},
                        "virulence_elements": {"rows": [["rmpA"]]},
                    },
                },
            }
        )
        self.assertEqual(checklist["status"], "ready")
        self.assertEqual(checklist["readiness"], "needs_review")
        self.assertGreaterEqual(checklist["blocking_count"], 1)
        labels = [item["label"] for item in checklist["items"]]
        self.assertIn("公共卫生意义已确认", labels)

    def test_checklist_keeps_archive_traceability_item(self) -> None:
        checklist = build_export_checklist({"task": {"id": "task-2"}, "sections": {}})
        action_ids = [item["id"] for item in checklist["items"]]
        self.assertIn("database_trace", action_ids)
        self.assertIn("history_compare", action_ids)
        self.assertIn("export_archive", action_ids)

    def test_task_closure_state_blocks_export_until_review_import_and_history_are_done(self) -> None:
        checklist = build_export_checklist(
            {
                "task": {
                    "id": "task-3",
                    "name": "Needs review",
                    "closure_status": {
                        "state": "needs_report_review",
                        "label": "待复核",
                        "summary": "报告已就绪，但尚未由非任务创建人完成复核签核。",
                        "imported_sample_count": 0,
                        "steps": [
                            {"id": "analysis", "state": "done"},
                            {"id": "review_report", "state": "current"},
                            {"id": "import_sample", "state": "pending"},
                            {"id": "compare_history", "state": "pending"},
                            {"id": "archive_export", "state": "pending"},
                        ],
                    },
                },
                "sections": {"workflow_closure": {"risk_level": "低"}},
            }
        )
        states = {item["id"]: item["state"] for item in checklist["items"]}
        self.assertEqual(states["report_review"], "attention")
        self.assertEqual(states["database_trace"], "attention")
        self.assertEqual(states["history_compare"], "attention")
        self.assertEqual(states["export_archive"], "attention")
        self.assertGreaterEqual(checklist["blocking_count"], 4)

    def test_task_closure_state_marks_archive_ready_after_required_steps(self) -> None:
        checklist = build_export_checklist(
            {
                "task": {
                    "id": "task-4",
                    "name": "Ready archive",
                    "closure_status": {
                        "state": "needs_archive_export",
                        "label": "待归档",
                        "summary": "报告复核、样本入库和历史对照均已确认。",
                        "imported_sample_count": 1,
                        "steps": [
                            {"id": "analysis", "state": "done"},
                            {"id": "review_report", "state": "done"},
                            {"id": "import_sample", "state": "done"},
                            {"id": "compare_history", "state": "done"},
                            {"id": "archive_export", "state": "current"},
                        ],
                    },
                },
                "sections": {
                    "workflow_closure": {"risk_level": "低"},
                    "public_health_support": {"status": "empty"},
                    "resistance_virulence": {
                        "resistance_elements": {"rows": [["blaTEM"]]},
                        "virulence_elements": {"rows": [["fimH"]]},
                    },
                },
            }
        )
        states = {item["id"]: item["state"] for item in checklist["items"]}
        self.assertEqual(states["report_review"], "ready")
        self.assertEqual(states["database_trace"], "ready")
        self.assertEqual(states["history_compare"], "ready")
        self.assertEqual(states["export_archive"], "ready")
        self.assertEqual(checklist["readiness"], "ready")
        self.assertEqual(checklist["blocking_count"], 0)


if __name__ == "__main__":
    unittest.main()
