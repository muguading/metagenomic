from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable

from .admin_runtime import _load_conda_env_settings, _load_conda_root_setting, _resolve_workstation_module
from .application_errors import AuthorizationError
from .identity import UserIdentity
from .runtime_paths import _resolve_runtime_database_root
from .store import PortalStore
from .task_manager import AnalysisTaskManager, ValidationError
from .taxonomy_reports import _build_community_merge_bundle, _is_metagenome_task_record, _parse_community_merge_items


@dataclass(frozen=True)
class TaskSubmissionService:
    project_root: Path
    store: PortalStore
    task_manager: AnalysisTaskManager
    ensure_can_view_task: Callable[[dict], None]
    ensure_can_modify_task: Callable[[dict], None]
    ensure_user_has_virus_permission: Callable[[dict, str], None]
    resolve_requested_virus_permission: Callable[[object], str]
    resolve_pipeline_script: Callable[[Path, str], str]
    resolve_runtime_env_name: Callable[[Path, str], str]
    demo_type_to_virus_permission: dict[str, str]

    def create(self, payload: dict, *, identity: UserIdentity) -> dict:
        requested_workstation = _resolve_workstation_module(payload.get("workstation_key"))
        user = self.store.get_user(identity.username)
        if requested_workstation["key"] not in (user.get("allowed_modules") or []):
            raise AuthorizationError("当前账号未开通该分析模块权限")
        payload = self._prepare_community_merge(payload)
        analysis_target = str(payload.get("analysis_target") or "").strip().lower()
        if requested_workstation["key"] == "virus" or analysis_target == "virus":
            self.ensure_user_has_virus_permission(user, self.resolve_requested_virus_permission(payload.get("species")))
        workspace_root = self._workspace_root()
        pipeline_script_value = str(payload.get("pipeline_script") or "").strip() or self.store.get_setting(
            "pipeline_script", requested_workstation["pipeline_script"]
        )
        return self.task_manager.create_task(
            payload,
            owner=identity.username,
            owner_group=identity.group_name,
            **self._runtime_arguments(workspace_root, pipeline_script_value, payload=payload),
        )

    def rerun(self, task_id: str, payload: dict, *, identity: UserIdentity) -> dict:
        task = self.task_manager.get_task(task_id, log_lines=0, owner=None)
        self.ensure_can_modify_task(task)
        payload = self._prepare_community_merge(
            payload,
            fallback_name=str(task.get("name") or "community_merge"),
            workstation_fallback=str((task.get("params") or {}).get("workstation_key") or ""),
        )
        workspace_root = self._workspace_root()
        pipeline_script_value = str(payload.get("pipeline_script") or task.get("pipeline_script") or "").strip() or self.store.get_setting(
            "pipeline_script", "Bac_assemble_260112_newformat.py"
        )
        return self.task_manager.rerun_task(
            task_id,
            payload,
            owner=None,
            owner_group=identity.group_name,
            **self._runtime_arguments(workspace_root, pipeline_script_value, payload=payload),
        )

    def create_demo(self, payload: dict, *, identity: UserIdentity) -> dict:
        user = self.store.get_user(identity.username)
        demo_type = str(payload.get("demo_type") or "fastq").strip().lower()
        demo_virus_key = self.demo_type_to_virus_permission.get(demo_type, "")
        if demo_virus_key:
            if "virus" not in (user.get("allowed_modules") or []):
                raise AuthorizationError("当前账号未开通该分析模块权限")
            self.ensure_user_has_virus_permission(user, demo_virus_key)
        return self.task_manager.create_demo_task(
            owner=identity.username,
            owner_group=identity.group_name,
            demo_type=demo_type,
            workspace_root=self._workspace_root(),
            database_root=_resolve_runtime_database_root(),
        )

    def _prepare_community_merge(
        self, payload: dict, fallback_name: str = "community_merge", workstation_fallback: str = ""
    ) -> dict:
        if (
            str(payload.get("workstation_key") or workstation_fallback).strip() != "community"
            or str(payload.get("community_input_source") or "").strip() != "merge_tasks"
        ):
            return payload
        merge_entries = []
        for item in _parse_community_merge_items(payload.get("community_merge_tasks")):
            task = self.task_manager.get_task(item["task_id"], log_lines=0, owner=None)
            self.ensure_can_view_task(task)
            if not _is_metagenome_task_record(task):
                raise ValidationError(f"任务 {task.get('name') or item['task_id']} 不是宏基因组任务，不能用于群落汇总")
            merge_entries.append({"task": task, "group": item["group"]})
        return {
            **payload,
            **_build_community_merge_bundle(
                project_root=self.project_root,
                task_name=str(payload.get("task_name") or fallback_name).strip(),
                group_column=str(payload.get("community_group_column") or "Group").strip() or "Group",
                merge_entries=merge_entries,
            ),
        }

    def _workspace_root(self) -> Path:
        return Path(self.store.get_setting("workspace_root", str(self.project_root))).expanduser().resolve()

    def _pipeline_python_setting(self, payload: dict) -> str:
        workstation_key = str(payload.get("workstation_key") or "").strip().lower()
        analysis_target = str(payload.get("analysis_target") or "").strip().lower()
        if workstation_key == "virus" or analysis_target == "virus":
            return "ncov"
        return self.store.get_setting("pipeline_python", "base")

    def _runtime_arguments(self, workspace_root: Path, pipeline_script_value: str, *, payload: dict | None = None) -> dict:
        payload = payload or {}
        return {
            "pipeline_script": self.resolve_pipeline_script(workspace_root, pipeline_script_value),
            "pipeline_python": self.resolve_runtime_env_name(
                workspace_root, self._pipeline_python_setting(payload)
            ),
            "database_root": str(_resolve_runtime_database_root()),
            "conda_root": _load_conda_root_setting(self.store),
            "max_concurrent_tasks": int(self.store.get_setting("max_concurrent_tasks", "2") or "2"),
            "conda_envs": _load_conda_env_settings(self.store),
        }
