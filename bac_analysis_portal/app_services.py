from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from flask import Flask, current_app

from .batch_import_service import BatchImportService
from .service_groups import AccessControlServices, AuspiceServices, BatchInputServices, FilesystemServices, ReportServices, RuntimeServices
from .task_submission_service import TaskSubmissionService
from .task_lifecycle_service import TaskLifecycleService
from .dataset_service import DatasetService
from .admin_monitor_service import AdminMonitorService
from .portal_service import PortalService
from .admin_user_service import AdminUserService
from .admin_audit_service import AdminAuditService
from .knowledge_service import KnowledgeService
from .reference_database_service import ReferenceDatabaseService
from .admin_settings_service import AdminSettingsService
from .queue_maintenance_service import QueueMaintenanceService
from .sample_library_manager import SampleLibraryManager


APP_SERVICES_EXTENSION_KEY = "bac_analysis_portal.services"


@dataclass(frozen=True)
class AppServices:
    """Application composition root shared by route modules."""

    project_root: Path
    access: AccessControlServices
    filesystem: FilesystemServices
    runtime: RuntimeServices
    batch_inputs: BatchInputServices
    reports: ReportServices
    auspice: AuspiceServices
    sample_library: SampleLibraryManager
    batch_import_service: BatchImportService
    task_submission_service: TaskSubmissionService
    task_lifecycle_service: TaskLifecycleService
    dataset_service: DatasetService
    admin_monitor_service: AdminMonitorService
    portal_service: PortalService
    admin_user_service: AdminUserService
    admin_audit_service: AdminAuditService
    knowledge_service: KnowledgeService
    reference_database_service: ReferenceDatabaseService
    admin_settings_service: AdminSettingsService
    queue_maintenance_service: QueueMaintenanceService


def install_app_services(app: Flask, services: AppServices) -> None:
    app.extensions[APP_SERVICES_EXTENSION_KEY] = services


def get_app_services(app: Flask | None = None) -> AppServices:
    active_app = app if app is not None else current_app
    services = active_app.extensions.get(APP_SERVICES_EXTENSION_KEY)
    if not isinstance(services, AppServices):
        raise RuntimeError("AppServices have not been installed on this Flask application")
    return services
