from __future__ import annotations

import platform
import subprocess
from pathlib import Path

from .report_sources import _resolve_report_source
from .task_manager import ValidationError


def resolve_task_result_html(task: dict) -> Path | None:
    input_path = str((task.get("params") or {}).get("input_path", "")).strip()
    if not input_path:
        return None
    candidate = Path(input_path).expanduser().resolve()
    search_dir = candidate if candidate.is_dir() else candidate.parent
    candidates = sorted(search_dir.glob("*_bacgenome.html")) if search_dir.is_dir() else []
    return candidates[0] if candidates else None


def task_has_embedded_result(task: dict) -> bool:
    return str(task_report_availability(task).get("mode") or "") == "dynamic"


def task_report_availability(task: dict) -> dict[str, object]:
    static_result = resolve_task_result_html(task)
    if static_result is not None:
        return {"available": True, "mode": "static", "path": static_result, "reason": ""}
    report_source = _resolve_report_source(task)
    if bool(report_source.get("available")):
        return {"available": True, "mode": "dynamic", "path": None, "reason": ""}
    return {
        "available": False,
        "mode": "missing",
        "path": None,
        "reason": str(report_source.get("reason") or "尚未识别可浏览报告。"),
    }


def open_in_file_manager(target: Path) -> None:
    system_name = platform.system().lower()
    command = ["open", str(target)] if system_name == "darwin" else (["explorer", str(target)] if system_name == "windows" else ["xdg-open", str(target)])
    try:
        subprocess.Popen(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, start_new_session=True)
    except FileNotFoundError as exc:
        raise ValidationError("当前系统缺少可用的文件管理器打开命令") from exc
    except OSError as exc:
        raise ValidationError(f"无法打开分析文件夹: {target}") from exc
