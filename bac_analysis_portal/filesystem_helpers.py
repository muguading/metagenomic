from __future__ import annotations

from pathlib import Path

from .task_manager import ValidationError

def _resolve_browser_path(base_root: Path, path_arg: str) -> Path:
    raw = str(path_arg or "").strip()
    if not raw:
        return base_root
    candidate = Path(raw).expanduser()
    return candidate.resolve() if candidate.is_absolute() else (base_root / candidate).resolve()


def _to_browser_path(base_root: Path, path: Path) -> str:
    if path == base_root:
        return ""
    try:
        return str(path.relative_to(base_root))
    except ValueError:
        return str(path)


def _is_within_root(base_root: Path, candidate: Path) -> bool:
    try:
        candidate.relative_to(base_root)
        return True
    except ValueError:
        return False


def _assert_within_root(base_root: Path, candidate: Path) -> None:
    try:
        candidate.relative_to(base_root)
    except ValueError as exc:
        raise ValidationError(f"路径超出允许范围: {candidate}") from exc


def _resolve_optional_existing_path(project_root: Path, raw_value: object) -> str:
    text = str(raw_value or "").strip()
    if not text:
        return ""
    candidate = Path(text).expanduser()
    candidate = candidate.resolve() if candidate.is_absolute() else (project_root / candidate).resolve()
    if not candidate.exists():
        raise ValidationError(f"文件不存在: {candidate}")
    return str(candidate)
