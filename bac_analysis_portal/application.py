from __future__ import annotations

import os
from functools import partial
from pathlib import Path
from typing import Any, Mapping

from flask import Flask

from .access_control import (
    DEMO_TYPE_TO_VIRUS_PERMISSION,
    add_no_cache_headers,
    admin_required,
    can_view_task,
    can_view_uploaded_auspice,
    current_username,
    ensure_can_control_task,
    ensure_can_modify_task,
    ensure_can_view_task,
    ensure_can_view_uploaded_auspice,
    ensure_user_has_virus_permission,
    is_logged_in,
    login_required,
    resolve_requested_virus_permission,
)
from .admin_operations import AdminOperations, detect_offline_update_sources, resolve_admin_path, resolve_admin_root, run_offline_update, run_online_update
from .app_lifecycle import register_error_handlers, start_queue_watch_daemon
from .app_services import AppServices, install_app_services
from .audit_hooks import register_audit_hooks
from .auto_pathosource_service import _load_admin_pathosource_trigger_rules, _save_admin_pathosource_trigger_rules, maybe_auto_trigger_pathosource_for_meta_task
from .batch_scanner import scan_fastq_directory_for_batch_rows
from .batch_input_operations import write_batch_input
from .batch_import_service import BatchImportService
from .dataset_operations import (
    DatasetOperations,
    NEXTSTRAIN_BUILD_TASK_TYPE,
    NEXTSTRAIN_CLI_PATH,
    is_portal_only_task,
    load_uploaded_auspice_metadata,
    normalize_nextstrain_dataset_prefix,
    parse_sample_location_json,
    resolve_nextstrain_dataset_main_json,
    resolve_nextstrain_dataset_sidecar,
    serialize_nextstrain_build_task,
    serialize_uploaded_auspice_dataset,
    uploaded_auspice_dataset_path,
    uploaded_auspice_metadata_path,
    uploaded_auspice_root,
    list_uploaded_auspice_metadata,
    store_uploaded_auspice,
    delete_uploaded_auspice,
    write_nextstrain_manifest,
)
from .filesystem_operations import create_directory, rename_path
from .import_templates import _load_batch_import_precheck, _precheck_reference_batch_upload, _store_batch_import_precheck
from .portal_helpers import open_in_file_manager, resolve_task_result_html, task_report_availability
from .reference_database_ops import _batch_import_reference_records, _import_reference_records_from_source, _register_reference_database_record, _sanitize_host_slug
from .reference_download_manager import ReferenceDownloadManager
from .reference_operations import ReferenceOperations
from .report_formatters import _human_bp
from .report_payload import _build_report_payload
from .report_sources import _resolve_report_sample_name, _resolve_report_source
from .route_registry import register_app_routes
from .runtime_config import env_flag, resolve_pipeline_script, resolve_runtime_env_name
from .runtime_hooks import register_runtime_hooks
from .sample_library_manager import SampleLibraryManager
from .service_groups import AccessControlServices, AuspiceServices, BatchInputServices, FilesystemServices, ReportServices, RuntimeServices
from .store import PortalStore
from .task_manager import AnalysisTaskManager
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
from .security import configure_session_security, register_security_hooks, resolve_portal_mode, validate_security_config


def create_app(config: Mapping[str, Any] | None = None) -> Flask:
    package_root = Path(__file__).resolve().parent
    project_root = package_root.parent
    precheck_cache: dict[str, dict[str, object]] = {}
    app = Flask(__name__, template_folder=str(package_root / "templates"), static_folder=str(package_root / "static"))
    portal_mode = resolve_portal_mode(config)
    app.config.from_mapping(
        PORTAL_MODE=portal_mode,
        SECRET_KEY=os.environ.get("PORTAL_SECRET_KEY", "") if portal_mode == "production" else os.environ.get("PORTAL_SECRET_KEY", "bac-analysis-portal-dev-key"),
        INITIAL_ADMIN_PASSWORD=os.environ.get("PORTAL_INITIAL_ADMIN_PASSWORD", "") if portal_mode == "production" else os.environ.get("PORTAL_INITIAL_ADMIN_PASSWORD", "admin123"),
        PORTAL_COOKIE_SECURE=env_flag("PORTAL_COOKIE_SECURE", default=portal_mode == "production"),
        SEND_FILE_MAX_AGE_DEFAULT=0,
        ENABLE_KNOWLEDGE_BASE_TEST_PANEL=env_flag("PORTAL_ENABLE_KNOWLEDGE_BASE_PANEL", default=True),
        PROJECT_ROOT=str(project_root),
    )
    if config:
        app.config.from_mapping(config)
    validate_security_config(app)
    configure_session_security(app)

    task_manager = AnalysisTaskManager.from_project_root(project_root)
    store = PortalStore.from_project_root(
        project_root,
        initial_admin_password=str(app.config["INITIAL_ADMIN_PASSWORD"]),
        rotate_weak_admin=portal_mode == "production",
    )
    app.config.pop("INITIAL_ADMIN_PASSWORD", None)
    reference_downloads = ReferenceDownloadManager(store=store, project_root=project_root)
    sample_library = SampleLibraryManager(
        store=store,
        resolve_report_source=_resolve_report_source,
        build_report_payload=_build_report_payload,
        resolve_report_sample_name=_resolve_report_sample_name,
        human_bp=_human_bp,
        task_manager=task_manager,
    )
    app.config["PORTAL_STORE"] = store

    register_audit_hooks(app, store)
    register_security_hooks(app)
    register_runtime_hooks(app, project_root=project_root, store=store, task_manager=task_manager, is_logged_in=is_logged_in)
    reference_operations = ReferenceOperations(
        list_download_jobs=reference_downloads.list_jobs,
        delete_download_job=reference_downloads.delete,
        start_download_job=reference_downloads.start,
        import_from_source=_import_reference_records_from_source,
        register_record=_register_reference_database_record,
        batch_import=_batch_import_reference_records,
        precheck_batch=_precheck_reference_batch_upload,
        store_precheck=lambda **kwargs: _store_batch_import_precheck(precheck_cache, **kwargs),
        load_precheck=lambda **kwargs: _load_batch_import_precheck(precheck_cache, **kwargs),
    )
    dataset_operations = DatasetOperations(
        nextstrain_cli_path=NEXTSTRAIN_CLI_PATH,
        serialize_nextstrain_build_task=serialize_nextstrain_build_task,
        uploaded_auspice_root=uploaded_auspice_root,
        uploaded_auspice_metadata_path=uploaded_auspice_metadata_path,
        serialize_uploaded_auspice_dataset=serialize_uploaded_auspice_dataset,
        parse_sample_location_json=parse_sample_location_json,
        sanitize_slug=_sanitize_host_slug,
        list_uploaded_auspice_metadata=list_uploaded_auspice_metadata,
        store_uploaded_auspice=store_uploaded_auspice,
        delete_uploaded_auspice=delete_uploaded_auspice,
        write_nextstrain_manifest=write_nextstrain_manifest,
    )
    admin_operations = AdminOperations(
        load_pathosource_trigger_rules=lambda: _load_admin_pathosource_trigger_rules(store),
        save_pathosource_trigger_rules=lambda payload: _save_admin_pathosource_trigger_rules(store, payload),
        resolve_admin_root=resolve_admin_root,
        resolve_admin_path=resolve_admin_path,
        detect_offline_update_sources=detect_offline_update_sources,
        run_online_update=run_online_update,
        run_offline_update=run_offline_update,
    )
    cpu_cache: dict[str, float] = {}
    install_app_services(
        app,
        AppServices(
            project_root=project_root,
            access=AccessControlServices(
                login_required=login_required,
                admin_required=admin_required,
                is_logged_in=is_logged_in,
                can_view_task=can_view_task,
                ensure_can_view_task=ensure_can_view_task,
                ensure_can_modify_task=ensure_can_modify_task,
                ensure_can_control_task=ensure_can_control_task,
                can_view_uploaded_auspice=can_view_uploaded_auspice,
                ensure_can_view_uploaded_auspice=ensure_can_view_uploaded_auspice,
                demo_type_to_virus_permission=DEMO_TYPE_TO_VIRUS_PERMISSION,
                ensure_user_has_virus_permission=ensure_user_has_virus_permission,
                resolve_requested_virus_permission=resolve_requested_virus_permission,
            ),
            filesystem=FilesystemServices(
                open_in_file_manager=open_in_file_manager,
                create_directory=create_directory,
                rename_path=rename_path,
            ),
            runtime=RuntimeServices(
                resolve_pipeline_script=resolve_pipeline_script,
                resolve_runtime_env_name=resolve_runtime_env_name,
            ),
            batch_inputs=BatchInputServices(
                scan_fastq_directory_for_batch_rows=scan_fastq_directory_for_batch_rows,
                write_batch_input=write_batch_input,
            ),
            reports=ReportServices(
                resolve_task_result_html=resolve_task_result_html,
                task_report_availability=task_report_availability,
                resolve_report_source=_resolve_report_source,
                maybe_auto_trigger_pathosource_for_meta_task=partial(
                    maybe_auto_trigger_pathosource_for_meta_task,
                    store=store,
                    task_manager=task_manager,
                    project_root=project_root,
                    owner_fallback=current_username,
                ),
            ),
            auspice=AuspiceServices(
                nextstrain_build_task_type=NEXTSTRAIN_BUILD_TASK_TYPE,
                add_no_cache_headers=add_no_cache_headers,
                normalize_nextstrain_dataset_prefix=normalize_nextstrain_dataset_prefix,
                load_uploaded_auspice_metadata=load_uploaded_auspice_metadata,
                uploaded_auspice_dataset_path=uploaded_auspice_dataset_path,
                resolve_nextstrain_dataset_main_json=resolve_nextstrain_dataset_main_json,
                resolve_nextstrain_dataset_sidecar=resolve_nextstrain_dataset_sidecar,
            ),
            sample_library=sample_library,
            batch_import_service=BatchImportService(
                store=store,
                project_root=project_root,
                sample_manager=sample_library,
                reference_operations=reference_operations,
            ),
            task_submission_service=TaskSubmissionService(
                project_root=project_root,
                store=store,
                task_manager=task_manager,
                ensure_can_view_task=ensure_can_view_task,
                ensure_can_modify_task=ensure_can_modify_task,
                ensure_user_has_virus_permission=ensure_user_has_virus_permission,
                resolve_requested_virus_permission=resolve_requested_virus_permission,
                resolve_pipeline_script=resolve_pipeline_script,
                resolve_runtime_env_name=resolve_runtime_env_name,
                demo_type_to_virus_permission=DEMO_TYPE_TO_VIRUS_PERMISSION,
            ),
            task_lifecycle_service=TaskLifecycleService(
                store=store,
                task_manager=task_manager,
                sample_manager=sample_library,
                ensure_can_view_task=ensure_can_view_task,
                ensure_can_modify_task=ensure_can_modify_task,
                ensure_can_control_task=ensure_can_control_task,
                open_in_file_manager=open_in_file_manager,
            ),
            dataset_service=DatasetService(
                project_root=project_root,
                store=store,
                task_manager=task_manager,
                sample_manager=sample_library,
                operations=dataset_operations,
                nextstrain_build_task_type=NEXTSTRAIN_BUILD_TASK_TYPE,
                can_view_task=can_view_task,
                can_view_uploaded_auspice=can_view_uploaded_auspice,
                ensure_can_view_uploaded_auspice=ensure_can_view_uploaded_auspice,
                load_uploaded_auspice_metadata=load_uploaded_auspice_metadata,
                uploaded_auspice_dataset_path=uploaded_auspice_dataset_path,
                open_in_file_manager=open_in_file_manager,
            ),
            admin_monitor_service=AdminMonitorService(
                project_root=project_root,
                store=store,
                task_manager=task_manager,
                resolve_pipeline_script=resolve_pipeline_script,
                resolve_runtime_env_name=resolve_runtime_env_name,
            ),
            portal_service=PortalService(
                project_root=project_root,
                store=store,
                task_manager=task_manager,
                cpu_cache=cpu_cache,
                can_view_task=can_view_task,
                is_portal_only_task=is_portal_only_task,
                task_report_availability=task_report_availability,
            ),
            admin_user_service=AdminUserService(store=store),
            admin_audit_service=AdminAuditService(store=store),
            knowledge_service=KnowledgeService(store=store),
            reference_database_service=ReferenceDatabaseService(
                store=store, project_root=project_root, operations=reference_operations
            ),
            admin_settings_service=AdminSettingsService(
                project_root=project_root,
                store=store,
                operations=admin_operations,
                resolve_pipeline_script=resolve_pipeline_script,
                resolve_runtime_env_name=resolve_runtime_env_name,
                create_directory=create_directory,
            ),
            queue_maintenance_service=QueueMaintenanceService(store=store, task_manager=task_manager),
        ),
    )
    register_app_routes(app)
    register_error_handlers(app)
    start_queue_watch_daemon(app)
    return app
