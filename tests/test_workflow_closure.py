from __future__ import annotations

import unittest
from pathlib import Path

from bac_analysis_portal.result_page import build_result_back_target
from bac_analysis_portal.workflow_closure import build_workflow_closure


class WorkflowClosureTests(unittest.TestCase):
    def test_high_priority_public_health_match_requires_database_and_history_actions(self) -> None:
        closure = build_workflow_closure(
            {
                "task": {"id": "task-1", "params": {"method": "spades"}},
                "overview_metrics": [
                    {"key": "species_estimation", "items": [{"display": "Klebsiella pneumoniae"}]},
                    {"key": "total_bases", "display": "1.2 Gb"},
                ],
                "sections": {
                    "public_health_support": {
                        "status": "matched",
                        "china_priority_group": "医院感染重点",
                    },
                    "resistance_virulence": {
                        "resistance_elements": {"rows": [["blaKPC"], ["blaNDM"], ["oqxA"]]},
                        "virulence_elements": {"rows": [["rmpA"]]},
                    },
                },
            }
        )
        self.assertEqual(closure["risk_level"], "高")
        action_ids = [item["id"] for item in closure["actions"]]
        self.assertIn("import_sample", action_ids)
        self.assertIn("compare_history", action_ids)
        self.assertIn("trace_source", action_ids)
        self.assertIn("export_cdc_report", action_ids)
        hrefs = {item["id"]: item["href"] for item in closure["actions"]}
        self.assertIn("tab=queue", hrefs["review_report"])
        self.assertIn("closure_action=import_sample", hrefs["import_sample"])
        self.assertIn("closure_action=compare_history", hrefs["compare_history"])

    def test_low_priority_result_still_guides_report_review_and_archiving(self) -> None:
        closure = build_workflow_closure(
            {
                "task": {"id": "task-2", "params": {"method": "spades"}},
                "overview_metrics": [{"key": "species_estimation", "items": [{"display": "Background organism"}]}],
                "sections": {
                    "public_health_support": {"status": "unmatched"},
                    "resistance_virulence": {
                        "resistance_elements": {"rows": []},
                        "virulence_elements": {"rows": []},
                    },
                },
            }
        )
        self.assertEqual(closure["risk_level"], "低")
        action_ids = [item["id"] for item in closure["actions"]]
        self.assertEqual(action_ids[:2], ["review_report", "import_sample"])
        self.assertNotIn("trace_source", action_ids)

    def test_frontend_runtime_contains_closure_panel_renderer(self) -> None:
        runtime = (Path(__file__).parents[1] / "bac_analysis_portal" / "static" / "report_runtime.js").read_text(encoding="utf-8")
        self.assertIn("function renderWorkflowClosurePanel", runtime)
        self.assertGreaterEqual(runtime.count("renderWorkflowClosurePanel(data)"), 2)
        self.assertIn("function renderExportChecklistPanel", runtime)
        self.assertGreaterEqual(runtime.count("renderExportChecklistPanel(data)"), 2)
        self.assertIn("function confirmReportExportReadiness", runtime)
        self.assertIn("function buildReportExportChecklistSummaryMarkup", runtime)
        self.assertIn("function recordReportExport", runtime)
        self.assertIn("/report-exports", runtime)
        self.assertIn('const CLOSURE_SYNC_STORAGE_KEY = "bac-closure-sync-event"', runtime)
        self.assertIn("function notifyPortalClosureSync", runtime)
        self.assertIn("function showReportToast", runtime)
        self.assertIn("window.localStorage.setItem", runtime)
        self.assertIn('notifyPortalClosureSync(taskId, "archive_export", { format })', runtime)
        self.assertIn("showReportToast(exportRecord.message", runtime)
        self.assertGreaterEqual(runtime.count("await recordReportExport(format)"), 3)
        self.assertIn('const exportWindow = window.open("about:blank", "_blank");', runtime)
        self.assertIn("await recordReportExport(format);\n    downloadBlob", runtime)
        self.assertIn("report-document-review-summary", runtime)

    def test_report_result_back_links_preserve_task_focus(self) -> None:
        package_root = Path(__file__).parents[1] / "bac_analysis_portal"
        template = (package_root / "templates" / "result_report.html").read_text(encoding="utf-8")
        result_page = (package_root / "result_page.py").read_text(encoding="utf-8")
        queue_target = build_result_back_target("task-1")
        database_target = build_result_back_target("task-1", return_to="database", sample_key="main::task-1::sample-a")
        self.assertEqual(queue_target["label"], "返回任务列表")
        self.assertIn("tab=queue&task=task-1", queue_target["href"])
        self.assertEqual(database_target["label"], "返回样本列表")
        self.assertIn("tab=database", database_target["href"])
        self.assertIn("sample_key=main%3A%3Atask-1%3A%3Asample-a", database_target["href"])
        self.assertIn("back_href", template)
        self.assertIn("back_label", template)
        self.assertIn("data-back-label", template)
        self.assertIn("onclick='window.location.href={{", template)
        self.assertIn("build_result_back_target", result_page)
        runtime = (package_root / "static" / "report_runtime.js").read_text(encoding="utf-8")
        self.assertIn("backNode && !backNode.dataset.backLabel", runtime)


if __name__ == "__main__":
    unittest.main()
