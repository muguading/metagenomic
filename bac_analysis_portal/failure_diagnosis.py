from __future__ import annotations

import re
from pathlib import Path
from typing import Any


_DIAGNOSIS_RULES = (
    {
        "category": "dependency_missing",
        "label": "运行环境依赖缺失",
        "patterns": (
            re.compile(r"ModuleNotFoundError:\s*No module named ['\"]([^'\"]+)['\"]", re.I),
            re.compile(r"command not found:\s*(\S+)", re.I),
        ),
        "summary": "分析程序缺少运行依赖，当前失败不能用于判断样本质量。",
        "recommendation": "先由系统管理员修复运行环境或依赖路径，再重新运行任务。",
        "impact": "系统运行环境异常；当前任务未形成可用于疾控判读的分析结果。",
        "rerun_recommended": True,
    },
    {
        "category": "input_invalid",
        "label": "输入文件异常",
        "patterns": (
            re.compile(r"(?:input|fastq|fasta|file).*(?:not found|does not exist|empty|invalid)", re.I),
            re.compile(r"(?:not found|does not exist|empty|invalid).*(?:input|fastq|fasta|file)", re.I),
            re.compile(r"No such file or directory:\s*['\"]?([^'\"\n]+)", re.I),
            re.compile(r"unexpected end of file", re.I),
        ),
        "summary": "输入文件缺失、为空或格式异常，分析无法继续。",
        "recommendation": "核对样本文件完整性、命名和路径；修正或重新获取数据后重建任务。",
        "impact": "当前样本没有形成有效分析结果，不能据此作出阴性或排除判断。",
        "rerun_recommended": True,
    },
    {
        "category": "resource_exhausted",
        "label": "计算资源不足",
        "patterns": (
            re.compile(r"out of memory|cannot allocate memory|killed process|oom", re.I),
            re.compile(r"no space left on device|disk quota exceeded", re.I),
        ),
        "summary": "任务因内存或磁盘资源不足中断。",
        "recommendation": "释放或扩充计算资源后重新运行；重跑前核对输出目录剩余空间。",
        "impact": "属于计算资源异常，当前任务结果不完整，不能进入复核和归档主流程。",
        "rerun_recommended": True,
    },
    {
        "category": "permission_denied",
        "label": "文件权限异常",
        "patterns": (re.compile(r"permission denied|operation not permitted", re.I),),
        "summary": "分析程序没有读取输入或写入输出所需的权限。",
        "recommendation": "由系统管理员修复目录权限后重新运行任务。",
        "impact": "属于系统访问权限异常，当前任务未形成完整结果。",
        "rerun_recommended": True,
    },
    {
        "category": "pipeline_error",
        "label": "分析流程异常",
        "patterns": (
            re.compile(r"traceback \(most recent call last\)", re.I),
            re.compile(r"\berror\b|\bfailed\b|exception", re.I),
        ),
        "summary": "分析流程发生未归类异常，需要结合日志证据进一步复核。",
        "recommendation": "由分析人员或系统管理员复核日志；确认原因修复后优先重新运行。",
        "impact": "当前任务结果可能不完整，不能直接用于公共卫生判读。",
        "rerun_recommended": True,
    },
)


def build_failure_diagnosis(task: dict[str, Any], log_path: Path) -> dict[str, Any]:
    status = str(task.get("status") or "").strip().upper()
    if status == "STOPPED":
        return {
            "category": "manually_stopped",
            "label": "任务人工停止",
            "confidence": "confirmed",
            "summary": "任务在完成分析前被停止，需要确认停止原因及是否重建。",
            "recommendation": "确认停止原因；仍需结果时重建任务，否则填写终止依据并归档。",
            "impact": "当前任务没有形成完整可判读结果。",
            "evidence": [],
            "rerun_recommended": False,
        }
    if status != "FAILED":
        return {}

    log_text = ""
    if log_path.is_file():
        log_text = log_path.read_text(encoding="utf-8", errors="replace")
    for rule in _DIAGNOSIS_RULES:
        evidence = _match_evidence(log_text, rule["patterns"])
        if not evidence:
            continue
        return {
            "category": rule["category"],
            "label": rule["label"],
            "confidence": "high" if rule["category"] != "pipeline_error" else "medium",
            "summary": rule["summary"],
            "recommendation": rule["recommendation"],
            "impact": rule["impact"],
            "evidence": evidence,
            "rerun_recommended": rule["rerun_recommended"],
        }
    return {
        "category": "unknown",
        "label": "原因待人工确认",
        "confidence": "low",
        "summary": "日志中未识别到明确失败模式，需要人工复核。",
        "recommendation": "查看完整日志并记录最终判断；原因未确认前不建议直接接受失败归档。",
        "impact": "当前任务没有形成可确认的完整结果。",
        "evidence": [],
        "rerun_recommended": False,
    }


def _match_evidence(log_text: str, patterns: tuple[re.Pattern[str], ...]) -> list[str]:
    lines = [line.strip() for line in log_text.splitlines() if line.strip()]
    matches: list[str] = []
    for line in reversed(lines):
        if any(pattern.search(line) for pattern in patterns):
            compact = line if len(line) <= 240 else f"{line[:237]}..."
            if compact not in matches:
                matches.append(compact)
        if len(matches) >= 3:
            break
    return list(reversed(matches))
