from __future__ import annotations

from .parse_utils import _safe_float, _safe_int

def _human_count(value: object) -> str:
    number = _safe_int(value)
    if number is None:
        return "--"
    return f"{number:,}"

def _human_bp(value: object) -> str:
    number = _safe_int(value)
    if number is None:
        return "--"
    units = [("bp", 1), ("Kb", 10**3), ("Mb", 10**6), ("Gb", 10**9)]
    unit = "bp"
    divisor = 1
    for unit_name, base in units:
        if number >= base:
            unit = unit_name
            divisor = base
    scaled = number / divisor
    if divisor == 1:
        return f"{number} {unit}"
    return f"{scaled:.2f} {unit}"

def _coerce_percent_value(value: object) -> float | None:
    numeric = _safe_float(value)
    if numeric is None:
        return None
    if 0 <= numeric <= 1:
        return round(numeric * 100, 2)
    return round(numeric, 2)

def _display_percent(value: object) -> str:
    numeric = _coerce_percent_value(value)
    if numeric is None:
        return "--"
    return f"{numeric:.2f}%"

def _display_percent_points(value: object) -> str:
    numeric = _safe_float(value)
    if numeric is None:
        return "--"
    return f"{numeric:.2f}%"
