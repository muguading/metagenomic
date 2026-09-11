from __future__ import annotations

import os
import platform
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

def _collect_server_status(project_root: Path, cpu_cache: dict[str, float]) -> dict:
    memory = _read_memory_status()
    disk = _read_disk_status(project_root)
    cpu = _read_cpu_status(cpu_cache)
    return {
        "hostname": platform.node() or "-",
        "platform": platform.platform(),
        "machine": platform.machine() or "-",
        "python": sys.version.split()[0],
        "cpu_count": os.cpu_count() or 0,
        "sampled_at": datetime.now().isoformat(timespec="seconds"),
        "cpu": cpu,
        "memory": memory,
        "disk": disk,
    }


def _read_cpu_status(cpu_cache: dict[str, float]) -> dict:
    load_avg = _read_load_average()
    percent, source = _read_cpu_percent(cpu_cache)
    detail = f"{os.cpu_count() or 0} 核 / 负载 {load_avg}" if load_avg != "-" else f"{os.cpu_count() or 0} 核"
    return {
        "percent": percent,
        "source": source,
        "detail": detail,
        "load_average": load_avg,
    }


def _read_cpu_percent(cpu_cache: dict[str, float]) -> tuple[float | None, str]:
    proc_stat = Path("/proc/stat")
    if proc_stat.is_file():
        try:
            with proc_stat.open("r", encoding="utf-8") as handle:
                first_line = handle.readline().strip()
            parts = first_line.split()
            if len(parts) >= 5 and parts[0] == "cpu":
                values = [float(item) for item in parts[1:]]
                idle = values[3] + (values[4] if len(values) > 4 else 0.0)
                total = sum(values)
                previous_total = cpu_cache.get("total")
                previous_idle = cpu_cache.get("idle")
                cpu_cache["total"] = total
                cpu_cache["idle"] = idle
                if previous_total is not None and previous_idle is not None and total > previous_total:
                    idle_delta = idle - previous_idle
                    total_delta = total - previous_total
                    percent = max(0.0, min(100.0, (1.0 - idle_delta / total_delta) * 100.0))
                    return round(percent, 1), "proc-stat"
        except OSError:
            pass

    try:
        load = os.getloadavg()[0]
        cpu_total = max(os.cpu_count() or 1, 1)
        percent = max(0.0, min(100.0, load / cpu_total * 100.0))
        return round(percent, 1), "load-average"
    except (AttributeError, OSError):
        return None, "unavailable"


def _read_load_average() -> str:
    try:
        load1, load5, load15 = os.getloadavg()
        return f"{load1:.2f} / {load5:.2f} / {load15:.2f}"
    except (AttributeError, OSError):
        return "-"


def _read_memory_status() -> dict:
    meminfo = Path("/proc/meminfo")
    if meminfo.is_file():
        values: dict[str, int] = {}
        try:
            for line in meminfo.read_text(encoding="utf-8", errors="ignore").splitlines():
                if ":" not in line:
                    continue
                key, raw_value = line.split(":", 1)
                amount = raw_value.strip().split()[0]
                if amount.isdigit():
                    values[key] = int(amount) * 1024
            total = values.get("MemTotal", 0)
            available = values.get("MemAvailable", values.get("MemFree", 0))
            used = max(total - available, 0)
            percent = (used / total * 100.0) if total else None
            return {
                "total": total,
                "used": used,
                "free": available,
                "percent": round(percent, 1) if percent is not None else None,
                "total_human": _human_bytes(total),
                "used_human": _human_bytes(used),
                "free_human": _human_bytes(available),
            }
        except OSError:
            pass

    if platform.system() == "Darwin":
        darwin_memory = _read_memory_status_darwin()
        if darwin_memory is not None:
            return darwin_memory

    total = _read_total_memory_fallback()
    return {
        "total": total,
        "used": 0,
        "free": total,
        "percent": 0.0 if total else None,
        "total_human": _human_bytes(total),
        "used_human": _human_bytes(0),
        "free_human": _human_bytes(total),
    }


def _read_memory_status_darwin() -> dict | None:
    total = _read_total_memory_fallback()
    if total <= 0:
        return None

    try:
        vm_stat_output = subprocess.check_output(["vm_stat"], text=True, stderr=subprocess.DEVNULL)
    except (OSError, subprocess.SubprocessError):
        return None

    page_size = _read_vm_stat_page_size(vm_stat_output) or os.sysconf("SC_PAGE_SIZE")
    page_values: dict[str, int] = {}
    for raw_line in vm_stat_output.splitlines():
        line = raw_line.strip()
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        digits = "".join(ch for ch in value if ch.isdigit())
        if digits:
            page_values[key.strip()] = int(digits)

    free = (
        page_values.get("Pages free", 0)
        + page_values.get("Pages inactive", 0)
        + page_values.get("Pages speculative", 0)
    ) * int(page_size)
    free = max(0, min(free, total))
    used = max(total - free, 0)
    percent = round(used / total * 100.0, 1) if total else None
    return {
        "total": total,
        "used": used,
        "free": free,
        "percent": percent,
        "total_human": _human_bytes(total),
        "used_human": _human_bytes(used),
        "free_human": _human_bytes(free),
    }


def _read_vm_stat_page_size(vm_stat_output: str) -> int | None:
    for raw_line in vm_stat_output.splitlines():
        line = raw_line.strip()
        if "page size of" not in line:
            continue
        start = line.find("page size of")
        if start < 0:
            continue
        suffix = line[start + len("page size of") :]
        digits = "".join(ch for ch in suffix if ch.isdigit())
        if digits:
            return int(digits)
    return None


def _read_total_memory_fallback() -> int:
    try:
        page_size = os.sysconf("SC_PAGE_SIZE")
        page_count = os.sysconf("SC_PHYS_PAGES")
        return int(page_size) * int(page_count)
    except (ValueError, OSError, AttributeError):
        return 0


def _read_disk_status(project_root: Path) -> dict:
    usage = shutil.disk_usage(project_root)
    used_percent = (usage.used / usage.total * 100.0) if usage.total else None
    return {
        "path": str(project_root),
        "total": usage.total,
        "used": usage.used,
        "free": usage.free,
        "used_percent": round(used_percent, 1) if used_percent is not None else None,
        "total_human": _human_bytes(usage.total),
        "used_human": _human_bytes(usage.used),
        "free_human": _human_bytes(usage.free),
    }


def _human_bytes(value: int | float) -> str:
    size = float(value or 0)
    units = ["B", "KB", "MB", "GB", "TB", "PB"]
    for unit in units:
        if size < 1024.0 or unit == units[-1]:
            if unit == "B":
                return f"{int(size)} {unit}"
            return f"{size:.1f} {unit}"
        size /= 1024.0
    return "0 B"
