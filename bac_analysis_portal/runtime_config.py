from __future__ import annotations

import os
from pathlib import Path

from .filesystem_helpers import _assert_within_root
from .task_manager import ValidationError


def env_flag(name: str, default: bool = False) -> bool:
    raw_value = str(os.environ.get(name, "") or "").strip().lower()
    return default if not raw_value else raw_value not in {"0", "false", "no", "off"}


def resolve_runtime_env_name(project_root: Path, raw_value: str) -> str:
    raw = str(raw_value or "").strip()
    if not raw:
        raise ValidationError("运行环境 Conda 环境名不能为空")
    if "/" not in raw and "\\" not in raw:
        return raw
    candidate = Path(raw).expanduser()
    candidate = candidate.resolve() if candidate.is_absolute() else (project_root / candidate).resolve()
    if candidate.name.lower().endswith(".py") or "bac_assemble_260112_newformat.py" in candidate.name.lower():
        raise ValidationError("运行环境应填写 Conda 环境名，或选择该环境中的 python 可执行文件，不能填写分析脚本路径")
    parts = list(candidate.parts)
    if "envs" in parts and parts.index("envs") + 1 < len(parts):
        return parts[parts.index("envs") + 1]
    if candidate.name.lower().startswith("python"):
        return "base"
    if parts:
        return parts[-1]
    raise ValidationError("无法识别运行环境，请填写 Conda 环境名，或选择 envs/<环境名>/bin/python")


def resolve_pipeline_script(project_root: Path, script_path: str) -> str:
    raw = str(script_path or "").strip()
    if not raw:
        raise ValidationError("脚本路径不能为空")
    candidate = Path(raw).expanduser()
    candidate = candidate.resolve() if candidate.is_absolute() else (project_root / candidate).resolve()
    _assert_within_root(project_root.resolve(), candidate)
    return str(candidate)
