from __future__ import annotations

from pathlib import Path

from .filesystem_helpers import _assert_within_root
from .task_manager import ValidationError


def create_directory(parent: Path, name: str, *, allowed_root: Path | None = None) -> Path:
    normalized_name = str(name or "").strip()
    if not normalized_name:
        raise ValidationError("文件夹名称不能为空")
    if any(token in normalized_name for token in ("/", "\\")):
        raise ValidationError("文件夹名称不能包含路径分隔符")
    if not parent.is_dir():
        raise ValidationError(f"只能在目录下新建文件夹: {parent}")
    target = (parent / normalized_name).resolve()
    if allowed_root is not None:
        _assert_within_root(allowed_root.resolve(), target)
    try:
        target.mkdir(exist_ok=False)
    except PermissionError as exc:
        raise ValidationError(f"没有写入目录权限: {parent}") from exc
    return target


def rename_path(source: Path, new_name: str, *, allowed_root: Path) -> Path:
    normalized_name = str(new_name or "").strip()
    if not source.exists():
        raise ValidationError(f"待重命名目标不存在: {source}")
    if not normalized_name:
        raise ValidationError("新名称不能为空")
    if any(token in normalized_name for token in ("/", "\\")):
        raise ValidationError("新名称不能包含路径分隔符")
    target = (source.parent / normalized_name).resolve()
    _assert_within_root(allowed_root.resolve(), target)
    source.rename(target)
    return target
