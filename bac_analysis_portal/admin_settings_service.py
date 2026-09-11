from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable

from .admin_operations import AdminOperations
from .admin_runtime import (
    CONDA_ENV_SETTINGS,
    _detect_conda_root,
    _list_conda_env_names,
    _load_conda_env_settings,
    _resolve_home_and_desktop_paths,
)
from .nextclade_database import NextcladeDatabaseManager
from .runtime_paths import _normalize_database_root_path
from .store import PortalStore
from .task_manager import ValidationError


@dataclass(frozen=True)
class AdminSettingsService:
    project_root: Path
    store: PortalStore
    operations: AdminOperations
    resolve_pipeline_script: Callable[[Path, str], str]
    resolve_runtime_env_name: Callable[[Path, str], str]
    create_directory: Callable[..., Path]

    def get_settings(self) -> dict:
        detected_conda_root = _detect_conda_root()
        conda_root_value = str(
            self.store.get_setting("conda_root", str(detected_conda_root) if detected_conda_root else "") or ""
        ).strip()
        conda_root_path = Path(conda_root_value).expanduser().resolve() if conda_root_value else detected_conda_root
        return {
            "workspace_root": self.store.get_setting("workspace_root", str(self.project_root)),
            "pipeline_script": self.store.get_setting("pipeline_script", "Bac_assemble_260112_newformat.py"),
            "pipeline_python": self.resolve_runtime_env_name(
                self.project_root, self.store.get_setting("pipeline_python", "base")
            ),
            "conda_root": str(conda_root_path) if conda_root_path else "",
            "detected_conda_envs": _list_conda_env_names(conda_root_path),
            "database_root": self.store.get_setting("database_root", str(self.project_root)),
            "max_concurrent_tasks": int(self.store.get_setting("max_concurrent_tasks", "2") or "2"),
            "conda_envs": _load_conda_env_settings(self.store),
            "conda_env_fields": CONDA_ENV_SETTINGS,
            "pathosource_trigger_rules": self.operations.load_pathosource_trigger_rules(),
        }

    def list_conda_envs(self, root_value: str) -> dict:
        conda_root = Path(root_value).expanduser().resolve() if root_value else _detect_conda_root()
        if conda_root is None or not conda_root.is_dir():
            raise ValidationError(f"Conda 安装路径不存在: {root_value or conda_root or ''}")
        return {"conda_root": str(conda_root), "detected_conda_envs": _list_conda_env_names(conda_root)}

    def browse_filesystem(self, *, selector: str, root: str, path: str) -> dict:
        browse_root = self.operations.resolve_admin_root(root)
        current_path = self.operations.resolve_admin_path(browse_root, path)
        if not current_path.is_dir():
            raise ValidationError(f"只能浏览目录: {current_path}")
        try:
            children = sorted(current_path.iterdir(), key=lambda item: (not item.is_dir(), item.name.lower()))
        except PermissionError as exc:
            raise ValidationError(f"没有访问目录权限: {current_path}") from exc
        items = []
        for child in children:
            if child.name.startswith(".") and child.name not in {".", ".."}:
                continue
            item_type = "directory" if child.is_dir() else "file"
            if selector == "workspace_root" and item_type != "directory":
                continue
            try:
                modified_at = child.stat().st_mtime
            except OSError:
                modified_at = None
            items.append({"name": child.name, "type": item_type, "path": str(child.resolve()), "modified_at": modified_at})
        home_path, desktop_path = _resolve_home_and_desktop_paths()
        return {
            "root": str(browse_root),
            "current_path": str(current_path),
            "parent_path": str(current_path.parent.resolve()) if current_path != browse_root else "",
            "home_path": str(home_path),
            "desktop_path": str(desktop_path),
            "selector": selector,
            "items": items,
        }

    def create_folder(self, *, current_path: str, name: str) -> Path:
        parent = Path(str(current_path or "").strip() or "/").expanduser().resolve()
        return self.create_directory(parent, name)

    def update_settings(self, payload: dict) -> dict:
        workspace_root = str(payload.get("workspace_root", "")).strip()
        script_path = str(payload.get("pipeline_script", "")).strip()
        pipeline_python = str(payload.get("pipeline_python", "")).strip()
        conda_root = str(payload.get("conda_root", "")).strip()
        database_root = str(payload.get("database_root", "")).strip()
        if not all([workspace_root, script_path, pipeline_python, database_root]):
            raise ValidationError("部署基准目录、脚本路径、运行环境和数据库部署目录均不能为空")
        try:
            max_task_count = max(1, int(str(payload.get("max_concurrent_tasks", "")).strip() or "2"))
        except ValueError as exc:
            raise ValidationError("最大同时运行任务数量必须是整数") from exc
        workspace_root_path = Path(workspace_root).expanduser().resolve()
        conda_root_path = Path(conda_root).expanduser().resolve() if conda_root else _detect_conda_root()
        database_root_path = _normalize_database_root_path(Path(database_root))
        if not workspace_root_path.is_dir():
            raise ValidationError(f"部署基准目录不存在: {workspace_root_path}")
        if conda_root and (conda_root_path is None or not conda_root_path.is_dir()):
            raise ValidationError(f"Conda 安装路径不存在: {conda_root_path}")
        if not database_root_path.is_dir():
            raise ValidationError(f"数据库部署目录不存在: {database_root_path}")
        candidate = Path(self.resolve_pipeline_script(workspace_root_path, script_path))
        runtime_env = self.resolve_runtime_env_name(workspace_root_path, pipeline_python)
        raw_conda_envs = payload.get("conda_envs") or {}
        conda_envs = {}
        for item in CONDA_ENV_SETTINGS:
            raw_value = str(
                raw_conda_envs.get(item["key"], self.store.get_setting(f"conda_env_{item['key']}", item["default"])) or ""
            ).strip()
            if not raw_value:
                raise ValidationError(f"{item['label']} 环境名不能为空")
            conda_envs[item["key"]] = self.resolve_runtime_env_name(workspace_root_path, raw_value)
        if not candidate.is_file():
            raise ValidationError(f"脚本路径不存在: {candidate}")
        settings = {
            "workspace_root": str(workspace_root_path),
            "pipeline_script": str(candidate.relative_to(workspace_root_path)),
            "pipeline_python": runtime_env,
            "conda_root": str(conda_root_path) if conda_root_path else "",
            "database_root": database_root,
            "max_concurrent_tasks": str(max_task_count),
        }
        for key, value in settings.items():
            self.store.set_setting(key, value)
        for key, value in conda_envs.items():
            self.store.set_setting(f"conda_env_{key}", value)
        return {
            **settings,
            "max_concurrent_tasks": max_task_count,
            "detected_conda_envs": _list_conda_env_names(conda_root_path),
            "conda_envs": conda_envs,
            "conda_env_fields": CONDA_ENV_SETTINGS,
        }

    def workspace_root(self) -> Path:
        return Path(self.store.get_setting("workspace_root", str(self.project_root))).expanduser().resolve()

    def load_pathosource_trigger_rules(self) -> dict:
        return self.operations.load_pathosource_trigger_rules()

    def save_pathosource_trigger_rules(self, payload: dict) -> dict:
        if not isinstance(payload, dict):
            raise ValidationError("触发规则配置格式不正确")
        return self.operations.save_pathosource_trigger_rules(payload)

    def update_sources(self) -> list[dict]:
        return self.operations.detect_offline_update_sources(self.workspace_root())

    def run_online_update(self) -> dict:
        return self.operations.run_online_update(self.workspace_root())

    def run_offline_update(self, source_path: str) -> dict:
        return self.operations.run_offline_update(
            self.workspace_root(), Path(str(source_path or "").strip()).expanduser().resolve()
        )

    def check_nextclade_datasets(self) -> dict[str, object]:
        return NextcladeDatabaseManager(project_root=self.project_root, store=self.store).check()

    def update_nextclade_datasets(self, directories: object = None) -> dict[str, object]:
        return NextcladeDatabaseManager(project_root=self.project_root, store=self.store).update(directories)
