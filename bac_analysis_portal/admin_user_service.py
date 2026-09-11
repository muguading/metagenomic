from __future__ import annotations

from dataclasses import dataclass

from .store import PortalStore
from .task_manager import ValidationError


@dataclass(frozen=True)
class AdminUserService:
    store: PortalStore

    def list_users(self) -> list[dict]:
        return self.store.list_users()

    def create(self, payload: dict) -> dict:
        username = str(payload.get("username", "")).strip()
        password = str(payload.get("password", ""))
        if not username or not password:
            raise ValidationError("用户名和密码不能为空")
        return self.store.create_user(
            username=username,
            password=password,
            role=str(payload.get("role", "user")).strip(),
            display_name=str(payload.get("display_name", "")).strip(),
            group_name=str(payload.get("group_name", "")).strip(),
            allowed_viruses=payload.get("allowed_viruses") if isinstance(payload.get("allowed_viruses"), list) else None,
            module_expirations=payload.get("module_expirations") if isinstance(payload.get("module_expirations"), dict) else None,
            account_expires_at=str(payload.get("account_expires_at", "")).strip() or None,
            allowed_modules=payload.get("allowed_modules") if isinstance(payload.get("allowed_modules"), list) else None,
        )

    def update(self, username: str, payload: dict) -> dict:
        if username == "admin" and str(payload.get("role", "")).strip() in {"user", "group_admin"}:
            raise ValidationError("默认管理员账号不能降级")
        return self.store.update_user(
            username,
            new_username=str(payload.get("username", "")).strip() or None,
            role=str(payload.get("role", "")).strip() or None,
            group_name=str(payload.get("group_name", "")).strip() if "group_name" in payload else None,
            display_name=str(payload.get("display_name", "")).strip() if "display_name" in payload else None,
            new_password=str(payload.get("password", "")).strip() or None,
            allowed_viruses=payload.get("allowed_viruses") if isinstance(payload.get("allowed_viruses"), list) else None,
            module_expirations=payload.get("module_expirations") if isinstance(payload.get("module_expirations"), dict) else None,
            account_expires_at=str(payload.get("account_expires_at", "")).strip() if "account_expires_at" in payload else None,
            allowed_modules=payload.get("allowed_modules") if isinstance(payload.get("allowed_modules"), list) else None,
        )

    def delete(self, username: str, *, current_username: str) -> None:
        if current_username == username:
            raise ValidationError("不能删除当前登录用户")
        if username == "admin":
            raise ValidationError("默认管理员账号不能删除")
        self.store.delete_user(username)
