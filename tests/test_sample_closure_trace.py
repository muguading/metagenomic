from __future__ import annotations

import unittest

from bac_analysis_portal.sample_closure_trace import build_sample_closure_trace


class SampleClosureTraceTests(unittest.TestCase):
    def test_traceable_sample_links_back_to_report_and_history_comparison(self) -> None:
        trace = build_sample_closure_trace(
            {
                "task_id": "task-1",
                "task_name": "Outbreak isolate",
                "report_dir": "/tmp/report",
                "final_fasta_path": "/tmp/report/sample.final.fasta",
                "imported_at": "2026-06-08T08:00:00Z",
                "updated_at": "2026-06-08T09:00:00Z",
            },
            version_events=[
                {
                    "action": "publish",
                    "summary": "个人库样本经审核发布入主数据库",
                    "operator": "admin",
                    "version_label": "main-publish-20260608",
                    "created_at": "2026-06-08T08:00:00Z",
                }
            ],
        )
        self.assertEqual(trace["state"], "traceable")
        self.assertEqual(trace["report_href"], "/tasks/task-1/result-page")
        action_ids = [item["id"] for item in trace["actions"]]
        self.assertIn("open_report", action_ids)
        self.assertIn("compare_history", action_ids)
        hrefs = {item["id"]: item["href"] for item in trace["actions"]}
        self.assertEqual(hrefs["open_task"], "/workstation?tab=queue&task=task-1")
        self.assertIn("closure_action=compare_history", hrefs["compare_history"])
        self.assertEqual(trace["recent_events"][0]["action"], "publish")

    def test_traceable_sample_surfaces_source_task_closure_state(self) -> None:
        trace = build_sample_closure_trace(
            {
                "task_id": "task-1",
                "task_name": "任务一",
                "report_dir": "/tmp/report",
                "final_fasta_path": "/tmp/report/sample.final.fasta",
                "imported_at": "2026-06-08T08:00:00Z",
            },
            task_closure_status={"state": "needs_report_review", "label": "待复核"},
        )
        self.assertEqual(trace["task_closure_state"], "needs_report_review")
        self.assertFalse(trace["task_closed"])
        self.assertIn("尚不能作为最终交付证据", trace["summary"])
        self.assertIn("来源任务闭环", [item["label"] for item in trace["evidence"]])

    def test_incomplete_sample_reports_missing_evidence(self) -> None:
        trace = build_sample_closure_trace(
            {
                "sample_name": "sample-a",
                "task_name": "",
                "report_dir": "",
                "final_fasta_path": "",
                "imported_at": "",
            }
        )
        self.assertEqual(trace["state"], "incomplete")
        self.assertIn("来源任务", trace["missing"])
        self.assertIn("报告目录", trace["missing"])
        self.assertIn("Final FASTA", trace["missing"])
        self.assertFalse(trace["report_href"])


if __name__ == "__main__":
    unittest.main()
