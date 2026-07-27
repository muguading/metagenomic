from __future__ import annotations

from pathlib import Path

from flask import current_app

def _default_project_root() -> Path:
    return Path(__file__).resolve().parent.parent

def _normalize_database_root_path(candidate: Path) -> Path:
    resolved = candidate.expanduser().resolve()
    nested_database = resolved / "database"
    if resolved.name != "database" and nested_database.is_dir():
        return nested_database
    return resolved

def _resolve_runtime_database_root() -> Path:
    project_root = _default_project_root()
    store = None
    try:
        store = current_app.config.get("PORTAL_STORE")
        configured_project_root = current_app.config.get("PROJECT_ROOT")
        if configured_project_root:
            project_root = Path(str(configured_project_root)).expanduser().resolve()
    except RuntimeError:
        store = None

    if store is None:
        return _normalize_database_root_path(project_root)

    raw_value = str(store.get_setting("database_root", "") or "").strip()
    if not raw_value:
        return _normalize_database_root_path(project_root)
    try:
        candidate = _normalize_database_root_path(Path(raw_value))
    except OSError:
        return _normalize_database_root_path(project_root)
    return candidate if candidate.is_dir() else _normalize_database_root_path(project_root)
