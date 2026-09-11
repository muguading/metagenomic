from __future__ import annotations

import json
import os
import uuid
from dataclasses import dataclass
from pathlib import Path

from .application_errors import AuthorizationError
from .identity import UserIdentity
from .knowledge_base import load_knowledge_base_bundle
from .runtime_paths import _resolve_runtime_database_root
from .store import PortalStore
from .task_manager import ValidationError
from .virus_report_templates import (
    _list_virus_report_templates,
    _reset_virus_report_template_override,
    _rollback_virus_report_template_override,
    _save_virus_report_template_override,
)


PATHONET_KNOWLEDGE_BASE_ENV = "META_PATHONET_KNOWLEDGE_BASE"
PATHONET_DEFAULT_DOCUMENT = {
    "schema_version": "v1",
    "id": "pathonet_typing",
    "title": "PathoNet 重点血清型与毒力基因规则",
    "description": "供 PathoNet 在运行时标记重点关注血清型和毒力基因。",
}


def _pathonet_knowledge_base_path(database_root: Path) -> Path:
    configured_path = str(os.environ.get(PATHONET_KNOWLEDGE_BASE_ENV) or "").strip()
    if configured_path:
        return Path(configured_path).expanduser().resolve()
    return database_root / "knowledge_base" / "pathonet" / "pathonet_typing.json"


def _normalize_pathonet_values(value: object, *, field: str, entry_index: int) -> list[str]:
    if not isinstance(value, list):
        raise ValidationError(f"第 {entry_index + 1} 条 PathoNet 规则的 {field} 必须为列表")
    if len(value) > 250:
        raise ValidationError(f"第 {entry_index + 1} 条 PathoNet 规则的 {field} 最多保留 250 项")
    normalized: list[str] = []
    for item in value:
        if not isinstance(item, str):
            raise ValidationError(f"第 {entry_index + 1} 条 PathoNet 规则的 {field} 只能包含文本")
        text = item.strip()
        if not text:
            continue
        if len(text) > 200:
            raise ValidationError(f"第 {entry_index + 1} 条 PathoNet 规则的 {field} 单项不能超过 200 个字符")
        if text not in normalized:
            normalized.append(text)
    return normalized


def _normalize_pathonet_entries(value: object) -> list[dict[str, object]]:
    if not isinstance(value, list):
        raise ValidationError("PathoNet 规则必须包含 entries 列表")
    if len(value) > 200:
        raise ValidationError("PathoNet 规则最多保留 200 个物种")
    normalized: list[dict[str, object]] = []
    seen_species: set[str] = set()
    for index, item in enumerate(value):
        if not isinstance(item, dict):
            raise ValidationError(f"第 {index + 1} 条 PathoNet 规则必须为对象")
        species = str(item.get("species") or "").strip()
        if not species:
            raise ValidationError(f"第 {index + 1} 条 PathoNet 规则缺少运行物种键")
        if len(species) > 120:
            raise ValidationError(f"第 {index + 1} 条 PathoNet 规则的运行物种键不能超过 120 个字符")
        if species in seen_species:
            raise ValidationError(f"PathoNet 规则存在重复的运行物种键：{species}")
        seen_species.add(species)
        normalized.append(
            {
                "species": species,
                "serotype": _normalize_pathonet_values(item.get("serotype"), field="重点血清型", entry_index=index),
                "vfgene": _normalize_pathonet_values(item.get("vfgene"), field="重点毒力基因", entry_index=index),
            }
        )
    return normalized


def _read_pathonet_document(path: Path) -> dict[str, object]:
    if not path.is_file():
        return {**PATHONET_DEFAULT_DOCUMENT, "entries": []}
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValidationError(f"无法读取 PathoNet 知识库：{path}") from exc
    if not isinstance(payload, dict):
        raise ValidationError("PathoNet 知识库根节点必须为对象")
    return payload


def _write_pathonet_document(path: Path, entries: list[dict[str, object]]) -> None:
    document = {**PATHONET_DEFAULT_DOCUMENT, "entries": entries}
    temporary_path = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary_path.write_text(json.dumps(document, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        temporary_path.replace(path)
    except OSError as exc:
        raise ValidationError(f"PathoNet 知识库保存失败：{path}") from exc
    finally:
        if temporary_path.exists():
            temporary_path.unlink(missing_ok=True)


@dataclass(frozen=True)
class KnowledgeService:
    store: PortalStore

    def summary(self) -> dict:
        bundle = self.bundle()
        return {
            "status": bundle.get("status", "empty"),
            "root": bundle.get("root", ""),
            "manifest": bundle.get("manifest", {}),
            "summary": bundle.get("summary", {}),
            "validation": bundle.get("validation", {}),
        }

    def bundle(self) -> dict:
        database_root = _resolve_runtime_database_root()
        bundle = load_knowledge_base_bundle(str(database_root))
        pathonet_path = _pathonet_knowledge_base_path(database_root)
        validation = bundle.setdefault("validation", {})
        try:
            document = _read_pathonet_document(pathonet_path)
            entries = _normalize_pathonet_entries(document.get("entries"))
            pathonet_validation: list[str] = []
        except ValidationError as exc:
            entries = []
            pathonet_validation = [str(exc)]
        collections = bundle.setdefault("collections", {})
        collections["pathonet_rules"] = entries
        validation["pathonet_rules"] = pathonet_validation
        validation["warnings"] = sum(
            len(items) for key, items in validation.items() if key != "warnings" and isinstance(items, list)
        )
        summary = bundle.setdefault("summary", {})
        summary["pathonet_rule_count"] = len(entries)
        bundle["pathonet_configuration"] = {
            "source": "environment" if str(os.environ.get(PATHONET_KNOWLEDGE_BASE_ENV) or "").strip() else "database_root",
            "path": str(pathonet_path),
        }
        return bundle

    def update_pathonet_rules(self, payload: dict, *, identity: UserIdentity) -> dict:
        self._ensure_admin(identity)
        if not isinstance(payload, dict):
            raise ValidationError("PathoNet 规则内容格式不正确")
        entries = _normalize_pathonet_entries(payload.get("entries"))
        _write_pathonet_document(_pathonet_knowledge_base_path(_resolve_runtime_database_root()), entries)
        return self.bundle()

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
