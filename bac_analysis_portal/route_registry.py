from __future__ import annotations

from flask import Flask

from .admin_audit_routes import register_admin_audit_routes
from .admin_settings_routes import register_admin_settings_routes
from .admin_monitor_routes import register_admin_monitor_routes
from .admin_user_routes import register_admin_user_routes
from .auspice_routes import register_auspice_routes
from .batch_input_routes import register_batch_input_routes
from .export_routes import register_export_routes
from .filesystem_routes import register_filesystem_routes
from .knowledge_routes import register_knowledge_routes
from .modeling_routes import register_modeling_routes
from .portal_routes import register_portal_routes
from .sample_library_routes import register_sample_library_routes
from .sample_modeling_routes import register_sample_modeling_routes
from .reference_database_routes import register_reference_database_routes
from .dataset_routes import register_dataset_routes
from .task_lifecycle_routes import register_task_lifecycle_routes
from .task_report_routes import register_task_report_routes
from .task_submission_routes import register_task_submission_routes


ROUTE_REGISTRARS = (
    register_knowledge_routes,
    register_portal_routes,
    register_filesystem_routes,
    register_task_lifecycle_routes,
    register_auspice_routes,
    register_task_report_routes,
    register_export_routes,
    register_task_submission_routes,
    register_batch_input_routes,
    register_admin_user_routes,
    register_admin_audit_routes,
    register_admin_settings_routes,
    register_admin_monitor_routes,
    register_sample_library_routes,
    register_sample_modeling_routes,
    register_modeling_routes,
    register_reference_database_routes,
    register_dataset_routes,
)


def register_app_routes(app: Flask) -> None:
    for registrar in ROUTE_REGISTRARS:
        registrar(app)
