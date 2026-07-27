from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable


@dataclass(frozen=True)
class AccessControlServices:
    login_required: Callable
    admin_required: Callable
    is_logged_in: Callable[[], bool]
    can_view_task: Callable[[dict], bool]
    ensure_can_view_task: Callable[[dict], None]
    ensure_can_modify_task: Callable[[dict], None]
    ensure_can_control_task: Callable[[dict], None]
    can_view_uploaded_auspice: Callable[[dict], bool]
    ensure_can_view_uploaded_auspice: Callable[[dict], None]
    demo_type_to_virus_permission: dict[str, str]
    ensure_user_has_virus_permission: Callable[[dict, str], None]
    resolve_requested_virus_permission: Callable[[Any], str]


@dataclass(frozen=True)
class FilesystemServices:
    open_in_file_manager: Callable[[Path], None]
    create_directory: Callable[..., Path]
    rename_path: Callable[..., Path]


@dataclass(frozen=True)
class RuntimeServices:
    resolve_pipeline_script: Callable[[Path, str], str]
    resolve_runtime_env_name: Callable[[Path, str], str]


@dataclass(frozen=True)
class BatchInputServices:
    scan_fastq_directory_for_batch_rows: Callable[[Path, str], list[list[str]]]
    write_batch_input: Callable[[Path, list[list[str]]], Path]


@dataclass(frozen=True)
class ReportServices:
    resolve_task_result_html: Callable[[dict], Path | None]
    task_report_availability: Callable[[dict], dict[str, object]]
    resolve_report_source: Callable[[dict, str], dict]
    maybe_auto_trigger_pathosource_for_meta_task: Callable[..., dict]


@dataclass(frozen=True)
class AuspiceServices:
    nextstrain_build_task_type: str
    add_no_cache_headers: Callable
    normalize_nextstrain_dataset_prefix: Callable
    load_uploaded_auspice_metadata: Callable
    uploaded_auspice_dataset_path: Callable
    resolve_nextstrain_dataset_main_json: Callable
    resolve_nextstrain_dataset_sidecar: Callable
