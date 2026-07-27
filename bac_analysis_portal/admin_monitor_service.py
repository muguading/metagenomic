from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable

from .admin_runtime import (
    _build_admin_realtime_monitor_payload,
    _load_conda_env_settings,
    _load_conda_root_setting,
    _resolve_workstation_module,
)
from .identity import UserIdentity
from .runtime_paths import _resolve_runtime_database_root
from .store import PortalStore
from .task_manager import AnalysisTaskManager, ValidationError


@dataclass(frozen=True)
class AdminMonitorService:
    project_root: Path
    store: PortalStore
    task_manager: AnalysisTaskManager
    resolve_pipeline_script: Callable[[Path, str], str]
    resolve_runtime_env_name: Callable[[Path, str], str]

    def create_watch_tasks(self, payload: dict, *, identity: UserIdentity) -> list[dict]:
        monitor_name = str(payload.get("monitor_name", "") or "").strip()
        input_path = self._optional_directory(payload.get("input_path"), "默认测序下机目录")
        output_root = self._optional_directory(payload.get("output_root"), "默认分析输出根目录")
        modules = payload.get("modules") if isinstance(payload.get("modules"), list) else []
        if not modules:
            raise ValidationError("请至少选择一个需要运行的模块")
        try:
            stable_minutes = max(1, int(payload.get("watch_stable_minutes", 30)))
            poll_minutes = max(1, int(payload.get("watch_poll_minutes", 5)))
            max_samples = max(0, int(payload.get("watch_max_samples", 0)))
        except (TypeError, ValueError) as exc:
            raise ValidationError("监听时长、轮询间隔和样本上限必须是整数") from exc

        workspace_root = Path(self.store.get_setting("workspace_root", str(self.project_root))).expanduser().resolve()
        pipeline_python = self.resolve_runtime_env_name(
            workspace_root, self.store.get_setting("pipeline_python", "base")
        )
        conda_envs = _load_conda_env_settings(self.store)
        created_items = []
        for module in modules:
            if not isinstance(module, dict):
                continue
            module_key = str(module.get("key", "") or "").strip().lower()
            workstation = _resolve_workstation_module(module_key)
            resolved_input = self._optional_directory(module.get("input_path"), workstation["label"]) or input_path
            resolved_output = self._optional_directory(module.get("output_root"), workstation["label"]) or output_root
            if resolved_input is None:
                raise ValidationError(f"{workstation['label']} 未配置有效的监控目录")
            if resolved_output is None:
                raise ValidationError(f"{workstation['label']} 未配置有效的输出目录")
            monitor_payload = _build_admin_realtime_monitor_payload(
                module_key=module_key,
                preset_key=str(module.get("preset", "") or "").strip(),
                preset_overrides=module.get("preset_overrides") if isinstance(module.get("preset_overrides"), dict) else {},
                monitor_name=monitor_name,
                input_path=resolved_input,
                output_root=resolved_output,
                watch_stable_minutes=stable_minutes,
                watch_poll_minutes=poll_minutes,
                watch_max_samples=max_samples,
            )
            pipeline_script_value = self.store.get_setting("pipeline_script", workstation["pipeline_script"])
            created = self.task_manager.create_task(
                monitor_payload,
                owner=identity.username or "admin",
                owner_group=identity.group_name,
                pipeline_script=self.resolve_pipeline_script(
                    workspace_root, str(pipeline_script_value or workstation["pipeline_script"])
                ),
                pipeline_python=pipeline_python,
                database_root=str(_resolve_runtime_database_root()),
                conda_root=_load_conda_root_setting(self.store),
                max_concurrent_tasks=int(self.store.get_setting("max_concurrent_tasks", "2") or "2"),
                conda_envs=conda_envs,
            )
            created_items.append(
                {
                    "id": created.get("id"),
                    "name": created.get("name"),
                    "workstation_key": module_key,
                    "watch_label": ((created.get("watch") or {}) if isinstance(created.get("watch"), dict) else {}).get("label")
                    or "监听中 · 等待测序数据写入",
                }
            )
        if not created_items:
            raise ValidationError("没有生成任何实时监控任务")
        return created_items

    @staticmethod
    def _optional_directory(value: object, label: str) -> Path | None:
        text = str(value or "").strip()
        if not text:
            return None
        path = Path(text).expanduser().resolve()
        if not path.is_dir():
            raise ValidationError(f"{label}不存在: {path}")
        return path
