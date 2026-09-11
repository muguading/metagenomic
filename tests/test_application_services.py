from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import Mock

from bac_analysis_portal.admin_user_service import AdminUserService
from bac_analysis_portal.application_errors import AuthorizationError
from bac_analysis_portal.dataset_service import DatasetService
from bac_analysis_portal.identity import UserIdentity
from bac_analysis_portal.portal_helpers import task_report_availability
from bac_analysis_portal.portal_service import PortalService, _eligible_reviewer_usernames
from bac_analysis_portal.queue_maintenance_service import QueueMaintenanceService
from bac_analysis_portal.sample_library_manager import SampleLibraryManager
from bac_analysis_portal import serotype_reports
from bac_analysis_portal.task_manager import ValidationError
from bac_analysis_portal.task_manager import AnalysisTaskManager
from bac_analysis_portal.task_submission_service import TaskSubmissionService


class ApplicationServiceTests(unittest.TestCase):
    def test_task_serialization_preserves_portal_only_task_type(self) -> None:
        manager = AnalysisTaskManager(project_root=Path("/tmp/project"), task_root=Path("/tmp/tasks"), python_executable="python")
        task = {
            "id": "task-1",
            "task_type": "nextstrain_auspice_build",
            "nextstrain_cli": "/opt/nextstrain",
            "log_path": "",
        }
        serialized = manager._serialize_task(task, include_log=False)
        self.assertEqual(serialized["task_type"], "nextstrain_auspice_build")
        self.assertEqual(serialized["nextstrain_cli"], "/opt/nextstrain")

    def test_influenza_report_dependencies_are_wired(self) -> None:
        self.assertTrue(callable(serotype_reports._read_influenza_resistance_annotation_table))
        self.assertTrue(callable(serotype_reports._build_influenza_mutation_display_table))

    def test_dynamic_report_source_is_recognized_as_browsable_report(self) -> None:
        with TemporaryDirectory() as directory:
            report_dir = Path(directory)
            (report_dir / "summary.tsv").write_text("sum_len\n100\n", encoding="utf-8")
            availability = task_report_availability(
                {
                    "id": "task-1",
                    "status": "SUCCEEDED",
                    "params": {"input_path": directory, "output_dir": directory, "analysis_target": "virus"},
                }
            )
        self.assertTrue(availability["available"])
        self.assertEqual(availability["mode"], "dynamic")

    def test_eligible_reviewers_respect_role_group_and_task_ownership(self) -> None:
        users = [
            {"username": "analyst", "role": "group_admin", "group_name": "Bio"},
            {"username": "bio-reviewer", "role": "group_admin", "group_name": "Bio"},
            {"username": "other-group", "role": "group_admin", "group_name": "Other"},
            {"username": "regular-user", "role": "user", "group_name": "Bio"},
            {"username": "admin", "role": "admin", "group_name": ""},
        ]
        reviewers = _eligible_reviewer_usernames({"owner": "analyst", "owner_group": "Bio"}, users)
        self.assertEqual(reviewers, ["bio-reviewer", "admin"])

    def test_task_submission_rejects_unlicensed_module_before_creating_task(self) -> None:
        store = Mock()
        store.get_user.return_value = {"allowed_modules": ["bacteria"]}
        task_manager = Mock()
        service = TaskSubmissionService(
            project_root=Path("/tmp/project"),
            store=store,
            task_manager=task_manager,
            ensure_can_view_task=Mock(),
            ensure_can_modify_task=Mock(),
            ensure_user_has_virus_permission=Mock(),
            resolve_requested_virus_permission=Mock(return_value=""),
            resolve_pipeline_script=Mock(),
            resolve_runtime_env_name=Mock(),
            demo_type_to_virus_permission={},
        )
        with self.assertRaises(AuthorizationError):
            service.create({"workstation_key": "virus"}, identity=UserIdentity("user", "user", "group"))
        task_manager.create_task.assert_not_called()

    def test_task_submission_defaults_virus_tasks_to_ncov_runtime(self) -> None:
        store = Mock()
        store.get_user.return_value = {"allowed_modules": ["virus"]}
        store.get_setting.side_effect = lambda _key, default="": default
        task_manager = Mock()
        task_manager.create_task.return_value = {"id": "task-1"}
        service = TaskSubmissionService(
            project_root=Path("/tmp/project"),
            store=store,
            task_manager=task_manager,
            ensure_can_view_task=Mock(),
            ensure_can_modify_task=Mock(),
            ensure_user_has_virus_permission=Mock(),
            resolve_requested_virus_permission=Mock(return_value="ncov"),
            resolve_pipeline_script=Mock(return_value="/tmp/project/Bac_assemble_260112_newformat.py"),
            resolve_runtime_env_name=Mock(side_effect=lambda _root, value: value),
            demo_type_to_virus_permission={},
        )
        service.create(
            {
                "workstation_key": "virus",
                "analysis_target": "virus",
                "pipeline_script": "Bac_assemble_260112_newformat.py",
                "species": "SARS-CoV-2",
            },
            identity=UserIdentity("user", "user", "group"),
        )
        self.assertEqual(task_manager.create_task.call_args.kwargs["pipeline_python"], "ncov")

    def test_dataset_service_rejects_invalid_auspice_json_before_persisting(self) -> None:
        operations = Mock()
        service = DatasetService(
            project_root=Path("/tmp/project"),
            store=Mock(),
            task_manager=Mock(),
            sample_manager=Mock(),
            operations=operations,
            nextstrain_build_task_type="nextstrain",
            can_view_task=Mock(),
            can_view_uploaded_auspice=Mock(),
            ensure_can_view_uploaded_auspice=Mock(),
            load_uploaded_auspice_metadata=Mock(),
            uploaded_auspice_dataset_path=Mock(),
            open_in_file_manager=Mock(),
        )
        with self.assertRaises(ValidationError):
            service.create_auspice_upload(
                raw_bytes=b"not-json",
                original_filename="data.json",
                name="",
                identity=UserIdentity("user", "user", "group"),
            )
        operations.store_uploaded_auspice.assert_not_called()

    def test_admin_user_service_protects_default_admin(self) -> None:
        service = AdminUserService(store=Mock())
        with self.assertRaises(ValidationError):
            service.update("admin", {"role": "user"})
        with self.assertRaises(ValidationError):
            service.delete("admin", current_username="someone")

    def test_queue_maintenance_coordinates_refresh_and_reconcile(self) -> None:
        store = Mock()
        store.get_setting.return_value = "3"
        task_manager = Mock()
        QueueMaintenanceService(store=store, task_manager=task_manager).maintain_once()
        task_manager.refresh_monitored_tasks.assert_called_once_with()
        task_manager.reconcile_queue.assert_called_once_with(3)

    def test_portal_service_adds_closure_status_to_visible_tasks(self) -> None:
        store = Mock()
        store.get_setting.return_value = "2"
        store.list_users.return_value = [
            {"username": "admin", "role": "admin", "group_name": ""},
            {"username": "reviewer", "role": "admin", "group_name": ""},
        ]
        store.count_sample_library_records_by_task_id.return_value = 1
        store.list_audit_logs_for_target.return_value = [
            {"action": "确认导出归档", "username": "admin", "outcome": "success", "created_at": "2026-06-08T10:00:00"},
            {
                "action": "确认历史对照",
                "username": "admin",
                "outcome": "success",
                "created_at": "2026-06-08T09:00:00",
                "request_summary": '{"json": {"comparison_outcome": "no_signal", "note": "未发现聚集"}}',
            },
            {
                "action": "确认报告复核",
                "username": "admin",
                "outcome": "success",
                "created_at": "2026-06-08T08:00:00",
                "request_summary": '{"json": {"review_outcome": "approved", "note": "已复核"}}',
            },
        ]
        task_manager = Mock()
        task_manager.list_tasks.return_value = [
            {"id": "task-1", "name": "Task 1", "status": "SUCCEEDED", "owner": "analyst", "params": {}}
        ]
        service = PortalService(
            project_root=Path("/tmp/project"),
            store=store,
            task_manager=task_manager,
            cpu_cache={},
            can_view_task=Mock(return_value=True),
            is_portal_only_task=Mock(return_value=False),
            task_report_availability=Mock(return_value={"available": True, "mode": "dynamic", "reason": ""}),
        )
        tasks = service.list_visible_tasks()
        self.assertEqual(tasks[0]["closure_status"]["state"], "closed")
        self.assertEqual(tasks[0]["closure_status"]["imported_sample_count"], 1)
        self.assertEqual(tasks[0]["closure_status"]["recent_events"][0]["action"], "确认导出归档")

    def test_portal_service_can_enrich_single_task_with_closure_status(self) -> None:
        store = Mock()
        store.list_users.return_value = [
            {"username": "analyst", "role": "user", "group_name": "Bio"},
            {"username": "reviewer", "role": "group_admin", "group_name": "Bio"},
        ]
        store.count_sample_library_records_by_task_id.return_value = 0
        store.list_audit_logs_for_target.return_value = []
        task_manager = Mock()
        service = PortalService(
            project_root=Path("/tmp/project"),
            store=store,
            task_manager=task_manager,
            cpu_cache={},
            can_view_task=Mock(return_value=True),
            is_portal_only_task=Mock(return_value=False),
            task_report_availability=Mock(return_value={"available": True, "mode": "dynamic", "reason": ""}),
        )
        task = {"id": "task-2", "name": "Task 2", "status": "SUCCEEDED", "owner": "analyst", "owner_group": "Bio", "params": {}}
        enriched = service.enrich_task_with_closure_status(task)
        self.assertIs(enriched, task)
        self.assertEqual(task["closure_status"]["state"], "needs_report_review")
        self.assertEqual(task["closure_status"]["responsibility"]["eligible_reviewers"], ["reviewer"])

    def test_sample_library_manager_adds_closure_trace_to_visible_records(self) -> None:
        store = Mock()
        store.list_sample_library_by_scope.return_value = [
            {
                "sample_key": "main::task-1::sample-a",
                "sample_name": "sample-a",
                "task_id": "task-1",
                "task_name": "Task 1",
                "owner": "user",
                "owner_group": "group",
                "library_scope": "main",
                "visibility_scope": "group",
                "report_dir": "/tmp/report",
                "final_fasta_path": "/tmp/report/sample-a.final.fasta",
                "imported_at": "2026-06-08T08:00:00Z",
                "updated_at": "2026-06-08T09:00:00Z",
                "custom_metadata_json": "[]",
            }
        ]
        store.find_pending_submission_for_sample.return_value = None
        store.list_sample_library_version_logs_for_sample.return_value = [
            {
                "action": "publish",
                "summary": "个人库样本经审核发布入主数据库",
                "operator": "admin",
                "version_label": "main-publish-20260608",
                "created_at": "2026-06-08T08:00:00Z",
            }
        ]
        manager = SampleLibraryManager(
            store=store,
            resolve_report_source=Mock(),
            build_report_payload=Mock(),
            resolve_report_sample_name=Mock(),
            human_bp=Mock(return_value=""),
        )
        manager.list_metadata_templates = Mock(return_value=[])
        records = manager.list_visible(scope="main", role="user", username="user", group_name="group")
        self.assertEqual(records[0]["closure_trace"]["state"], "traceable")
        self.assertEqual(records[0]["closure_trace"]["report_href"], "/tasks/task-1/result-page")
        self.assertEqual(records[0]["closure_trace"]["recent_events"][0]["operator"], "admin")

    def test_sample_closure_trace_uses_live_task_report_availability(self) -> None:
        store = Mock()
        store.list_sample_library_by_scope.return_value = [
            {
                "sample_key": "main::task-1::sample-a",
                "sample_name": "sample-a",
                "task_id": "task-1",
                "task_name": "Task 1",
                "owner": "user",
                "owner_group": "group",
                "library_scope": "main",
                "visibility_scope": "group",
                "report_dir": "",
                "final_fasta_path": "/tmp/report/sample-a.final.fasta",
                "imported_at": "2026-06-08T08:00:00Z",
                "updated_at": "2026-06-08T09:00:00Z",
                "custom_metadata_json": "[]",
            }
        ]
        store.find_pending_submission_for_sample.return_value = None
        store.list_sample_library_version_logs_for_sample.return_value = []
        store.list_audit_logs_for_target.return_value = []
        store.count_sample_library_records_by_task_id.return_value = 1
        task_manager = Mock()
        task_manager.get_task.return_value = {"id": "task-1", "name": "Task 1", "status": "SUCCEEDED", "owner": "user", "params": {}}
        manager = SampleLibraryManager(
            store=store,
            resolve_report_source=Mock(return_value={"available": True, "report_dir": "/tmp/report"}),
            build_report_payload=Mock(),
            resolve_report_sample_name=Mock(),
            human_bp=Mock(return_value=""),
            task_manager=task_manager,
        )
        manager.list_metadata_templates = Mock(return_value=[])
        records = manager.list_visible(scope="main", role="user", username="user", group_name="group")
        trace = records[0]["closure_trace"]
        self.assertEqual(trace["report_href"], "/tasks/task-1/result-page")
        self.assertEqual(trace["task_closure_state"], "needs_report_review")
        self.assertNotEqual(trace["task_closure_state"], "report_missing")

    def test_dataset_service_delegates_valid_upload_to_operations(self) -> None:
        operations = Mock()
        operations.store_uploaded_auspice.return_value = {"id": "upload-1"}
        operations.serialize_uploaded_auspice_dataset.return_value = {"upload_id": "upload-1"}
        service = DatasetService(
            project_root=Path("/tmp/project"),
            store=Mock(),
            task_manager=Mock(),
            sample_manager=Mock(),
            operations=operations,
            nextstrain_build_task_type="nextstrain",
            can_view_task=Mock(),
            can_view_uploaded_auspice=Mock(),
            ensure_can_view_uploaded_auspice=Mock(),
            load_uploaded_auspice_metadata=Mock(),
            uploaded_auspice_dataset_path=Mock(),
            open_in_file_manager=Mock(),
        )
        result = service.create_auspice_upload(
            raw_bytes=b'{"meta":{"title":"demo"}}',
            original_filename="data.json",
            name="",
            identity=UserIdentity("user", "user", "group"),
        )
        self.assertEqual(result, {"upload_id": "upload-1"})
        operations.store_uploaded_auspice.assert_called_once()


if __name__ == "__main__":
    unittest.main()
