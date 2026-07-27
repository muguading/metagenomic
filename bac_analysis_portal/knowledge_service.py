from __future__ import annotations

from dataclasses import dataclass

from .application_errors import AuthorizationError
from .identity import UserIdentity
from .knowledge_base import load_knowledge_base_bundle, load_knowledge_base_summary
from .runtime_paths import _resolve_runtime_database_root
from .store import PortalStore
from .task_manager import ValidationError
from .virus_report_templates import (
    _list_virus_report_templates,
    _reset_virus_report_template_override,
    _rollback_virus_report_template_override,
    _save_virus_report_template_override,
)


@dataclass(frozen=True)
class KnowledgeService:
    store: PortalStore

    def summary(self) -> dict:
        return load_knowledge_base_summary(str(_resolve_runtime_database_root()))

    def bundle(self) -> dict:
        return load_knowledge_base_bundle(str(_resolve_runtime_database_root()))

    def list_virus_templates(self) -> list[dict]:
        return _list_virus_report_templates(self.store)

    def update_virus_template(self, template_id: str, payload: dict, *, identity: UserIdentity) -> dict:
        self._ensure_admin(identity)
        if not isinstance(payload, dict):
            raise ValidationError("模板内容格式不正确")
        return self._result(_save_virus_report_template_override(self.store, template_id.strip(), payload, identity.username))

    def reset_virus_template(self, template_id: str, *, identity: UserIdentity) -> dict:
        self._ensure_admin(identity)
        return self._result(_reset_virus_report_template_override(self.store, template_id.strip(), identity.username))

    def rollback_virus_template(self, template_id: str, payload: dict, *, identity: UserIdentity) -> dict:
        self._ensure_admin(identity)
        if not isinstance(payload, dict):
            raise ValidationError("回滚版本格式不正确")
        event_id = str(payload.get("eventId") or "").strip()
        if not event_id:
            raise ValidationError("缺少回滚版本 ID")
        return self._result(
            _rollback_virus_report_template_override(self.store, template_id.strip(), event_id, identity.username)
        )

    def _result(self, item: dict) -> dict:
        return {"item": item, "items": self.list_virus_templates()}

    @staticmethod
    def _ensure_admin(identity: UserIdentity) -> None:
        if identity.role != "admin":
            raise AuthorizationError("Administrator permission required")
