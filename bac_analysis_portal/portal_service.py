from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable

from .admin_runtime import _first_allowed_module, _resolve_workstation_module
from .store import PortalStore
from .system_status import _collect_server_status
from .task_closure_status import build_task_closure_status
from .task_manager import AnalysisTaskManager, ValidationError


@dataclass(frozen=True)
class PortalService:
    project_root: Path
    store: PortalStore
    task_manager: AnalysisTaskManager
    cpu_cache: dict[str, float]
    can_view_task: Callable[[dict], bool]
    is_portal_only_task: Callable[[dict], bool]
    task_report_availability: Callable[[dict], dict[str, object]]

    def authenticate(self, username: str, password: str) -> dict:
        try:
            existing_user = self.store.get_user(username)
        except KeyError:
            existing_user = None
        if existing_user and existing_user.get("is_expired"):
            raise ValidationError("账号已过有效期，请联系管理员续期")
        user = self.store.authenticate(username, password)
        if user is None:
            raise ValidationError("用户名或密码错误")
        return user

    def get_user(self, username: str) -> dict:
        return self.store.get_user(username)

    def default_module(self, username: str) -> str:
        return _first_allowed_module(self.get_user(username))

    def resolve_workstation(self, requested_module: object, username: str) -> tuple[dict, bool]:
        workstation = _resolve_workstation_module(requested_module)
        user = self.get_user(username)
        return workstation, workstation["key"] in (user.get("allowed_modules") or [])

    def health(self) -> dict:
        return {
            "status": "ok",
            "workspace_root": self.store.get_setting("workspace_root", str(self.project_root)),
            "script_path": self.store.get_setting("pipeline_script", "Bac_assemble_260112_newformat.py"),
        }

    def list_visible_tasks(self) -> list[dict]:
        self.task_manager.reconcile_queue(int(self.store.get_setting("max_concurrent_tasks", "2") or "2"))
        users = self.store.list_users()
        visible_tasks = [
            task
            for task in self.task_manager.list_tasks()
            if self.can_view_task(task) and not self.is_portal_only_task(task)
        ]
        for task in visible_tasks:
            self.enrich_task_with_closure_status(task, users=users)
        return visible_tasks

    def enrich_task_with_closure_status(self, task: dict, *, users: list[dict] | None = None) -> dict:
        report_availability = self.task_report_availability(task)
        result_exists = bool(report_availability.get("available"))
        task_id = str(task.get("id") or "")
        imported_count = self.store.count_sample_library_records_by_task_id(task_id)
        audit_events = self.store.list_audit_logs_for_target("task", task_id, limit=50)
        task["closure_status"] = build_task_closure_status(
            task,
            result_exists=result_exists,
            imported_sample_count=imported_count,
            audit_events=audit_events,
            eligible_reviewers=_eligible_reviewer_usernames(task, users if users is not None else self.store.list_users()),
            report_issue=str(report_availability.get("reason") or ""),
        )
        return task

    def server_status(self) -> dict:
        return _collect_server_status(self.project_root, self.cpu_cache)


def _eligible_reviewer_usernames(task: dict, users: list[dict]) -> list[str]:
    owner = str(task.get("owner") or "").strip()
    owner_group = str(task.get("owner_group") or "").strip()
    if not owner:
        return []
    eligible: list[str] = []
    for user in users:
        username = str(user.get("username") or "").strip()
        role = str(user.get("role") or "").strip()
        group_name = str(user.get("group_name") or "").strip()
        if not username or username == owner or bool(user.get("is_expired")):
            continue
        if role == "admin" or (role == "group_admin" and owner_group and group_name == owner_group):
            eligible.append(username)
    return eligible
