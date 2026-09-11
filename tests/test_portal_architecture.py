from __future__ import annotations

import inspect
from dataclasses import fields
from pathlib import Path
import unittest

from bac_analysis_portal.application import create_app
from bac_analysis_portal.app_services import APP_SERVICES_EXTENSION_KEY, get_app_services
from bac_analysis_portal.route_registry import ROUTE_REGISTRARS


class PortalArchitectureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.app = create_app({"TESTING": True, "SECRET_KEY": "test-secret"})

    def test_route_modules_only_depend_on_the_flask_app_at_registration(self) -> None:
        for registrar in ROUTE_REGISTRARS:
            with self.subTest(registrar=registrar.__name__):
                self.assertEqual(list(inspect.signature(registrar).parameters), ["app"])

    def test_application_services_are_installed_and_available_in_context(self) -> None:
        self.assertIn(APP_SERVICES_EXTENSION_KEY, self.app.extensions)
        with self.app.app_context():
            services = get_app_services()
        self.assertIsNotNone(services.portal_service)
        self.assertIsNotNone(services.queue_maintenance_service)

    def test_route_registry_preserves_the_full_route_surface(self) -> None:
        self.assertGreaterEqual(len(self.app.url_map._rules), 123)
        endpoints = {rule.endpoint for rule in self.app.url_map.iter_rules()}
        for endpoint in {"login", "list_tasks", "create_task", "task_report_data", "confirm_task_closure_action", "record_task_report_export", "record_task_failure_disposition", "filesystem"}:
            self.assertIn(endpoint, endpoints)

    def test_report_payload_exposes_workflow_closure_section(self) -> None:
        source = (Path(__file__).parents[1] / "bac_analysis_portal" / "report_payload.py").read_text(encoding="utf-8")
        self.assertIn('"workflow_closure"', source)
        self.assertIn("build_workflow_closure(payload)", source)
        self.assertIn('"export_checklist"', source)
        self.assertIn("build_export_checklist(payload)", source)

    def test_report_data_recomputes_export_checklist_with_live_closure_status(self) -> None:
        source = (Path(__file__).parents[1] / "bac_analysis_portal" / "task_report_routes.py").read_text(encoding="utf-8")
        self.assertIn("services.portal_service.enrich_task_with_closure_status(task)", source)
        self.assertIn('payload_task["closure_status"]', source)
        self.assertIn('sections["export_checklist"] = build_export_checklist(payload)', source)

    def test_dynamic_reports_use_the_platform_report_page(self) -> None:
        source = (Path(__file__).parents[1] / "bac_analysis_portal" / "task_lifecycle_routes.py").read_text(encoding="utf-8")
        self.assertIn('url_for("task_result_page"', source)
        self.assertIn('report_availability.get("mode") == "static"', source)

    def test_factory_configuration_overrides_defaults(self) -> None:
        self.assertTrue(self.app.testing)
        self.assertEqual(self.app.secret_key, "test-secret")

    def test_application_factory_contains_no_business_route_decorators(self) -> None:
        package_root = Path(__file__).parents[1] / "bac_analysis_portal"
        for filename in ("app.py", "application.py"):
            with self.subTest(filename=filename):
                source = (package_root / filename).read_text(encoding="utf-8")
                self.assertNotIn("@app.get(", source)
                self.assertNotIn("@app.post(", source)
                self.assertNotIn("@app.put(", source)
                self.assertNotIn("@app.delete(", source)

    def test_package_modules_do_not_depend_on_compatibility_entrypoint(self) -> None:
        package_root = Path(__file__).parents[1] / "bac_analysis_portal"
        offenders = []
        for module_path in package_root.glob("*.py"):
            if module_path.name == "app.py":
                continue
            if "from .app import" in module_path.read_text(encoding="utf-8"):
                offenders.append(module_path.name)
        self.assertEqual(offenders, [])

    def test_route_modules_do_not_write_directly_to_disk(self) -> None:
        package_root = Path(__file__).parents[1] / "bac_analysis_portal"
        forbidden = (".mkdir(", ".write_text(", ".write_bytes(", ".unlink(", ".rename(")
        offenders = []
        for module_path in package_root.glob("*_routes.py"):
            source = module_path.read_text(encoding="utf-8")
            if any(token in source for token in forbidden):
                offenders.append(module_path.name)
        self.assertEqual(offenders, [])

    def test_domain_operations_do_not_depend_on_flask_session(self) -> None:
        package_root = Path(__file__).parents[1] / "bac_analysis_portal"
        offenders = []
        for pattern in ("*operations.py", "*service.py"):
            for module_path in package_root.glob(pattern):
                if "from flask import session" in module_path.read_text(encoding="utf-8"):
                    offenders.append(module_path.name)
        self.assertEqual(sorted(set(offenders)), [])

    def test_application_services_do_not_depend_on_routes_or_composition_root(self) -> None:
        package_root = Path(__file__).parents[1] / "bac_analysis_portal"
        offenders = []
        for module_path in package_root.glob("*service.py"):
            source = module_path.read_text(encoding="utf-8")
            if "_routes import" in source or "from .application import" in source or "from .app_services import" in source:
                offenders.append(module_path.name)
        self.assertEqual(offenders, [])

    def test_cross_cutting_services_are_grouped_by_capability(self) -> None:
        with self.app.app_context():
            services = get_app_services()
        for grouped_name in ("access", "filesystem", "runtime", "batch_inputs"):
            self.assertIsNotNone(getattr(services, grouped_name))
        for legacy_name in (
            "login_required",
            "admin_required",
            "can_view_task",
            "open_in_file_manager",
            "create_directory",
            "resolve_pipeline_script",
            "scan_fastq_directory_for_batch_rows",
        ):
            self.assertFalse(hasattr(services, legacy_name), legacy_name)

    def test_app_services_exposes_no_raw_callbacks_or_infrastructure(self) -> None:
        with self.app.app_context():
            services = get_app_services()
        self.assertEqual([field.name for field in fields(services) if callable(getattr(services, field.name))], [])
        for forbidden in ("store", "task_manager", "cpu_cache", "reference_operations", "dataset_operations", "admin_operations"):
            self.assertFalse(hasattr(services, forbidden), forbidden)

    def test_batch_import_orchestration_is_not_implemented_in_routes(self) -> None:
        package_root = Path(__file__).parents[1] / "bac_analysis_portal"
        forbidden = ("create_batch_import_run(", "update_batch_import_run(", "_parse_database_batch_upload(")
        offenders = []
        for module_path in package_root.glob("*_routes.py"):
            source = module_path.read_text(encoding="utf-8")
            if any(token in source for token in forbidden):
                offenders.append(module_path.name)
        self.assertEqual(offenders, [])

    def test_task_command_routes_delegate_to_application_services(self) -> None:
        package_root = Path(__file__).parents[1] / "bac_analysis_portal"
        for filename in ("task_submission_routes.py", "task_lifecycle_routes.py"):
            with self.subTest(filename=filename):
                source = (package_root / filename).read_text(encoding="utf-8")
                self.assertNotIn("services.store", source)
                self.assertNotIn("services.task_manager", source)
                self.assertNotIn("task_manager.", source)

    def test_dataset_routes_delegate_to_dataset_service(self) -> None:
        source = (
            Path(__file__).parents[1] / "bac_analysis_portal" / "dataset_routes.py"
        ).read_text(encoding="utf-8")
        for forbidden in ("services.store", "services.task_manager", "task_manager.", "write_nextstrain_manifest"):
            self.assertNotIn(forbidden, source)

    def test_routes_do_not_create_analysis_tasks_directly(self) -> None:
        package_root = Path(__file__).parents[1] / "bac_analysis_portal"
        offenders = []
        for module_path in package_root.glob("*_routes.py"):
            if "task_manager.create_task(" in module_path.read_text(encoding="utf-8"):
                offenders.append(module_path.name)
        self.assertEqual(offenders, [])

    def test_routes_do_not_access_store_or_task_manager_directly(self) -> None:
        package_root = Path(__file__).parents[1] / "bac_analysis_portal"
        offenders = []
        for module_path in package_root.glob("*_routes.py"):
            source = module_path.read_text(encoding="utf-8")
            if "services.store" in source or "services.task_manager" in source:
                offenders.append(module_path.name)
        self.assertEqual(offenders, [])
        with self.app.app_context():
            services = get_app_services()
        self.assertFalse(hasattr(services, "store"))
        self.assertFalse(hasattr(services, "task_manager"))


if __name__ == "__main__":
    unittest.main()
