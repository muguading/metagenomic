#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


PLACEHOLDER_PARTS = (
    "replace-with",
    "<",
    ">",
    "请替换",
)

REQUIRED_PRODUCTION_ENV = (
    "PORTAL_MODE",
    "PORTAL_SECRET_KEY",
    "PORTAL_INITIAL_ADMIN_PASSWORD",
    "PORTAL_DB_PATH",
    "BAC_ANALYSIS_TASK_ROOT",
)

OPTIONAL_ASSET_PATHS = (
    "META_KRAKEN_DB",
    "META_VIRUS_KRAKEN_DB",
    "META_VIRSORTER2_DB",
    "META_CHECKV_DB",
    "META_GENOMAD_DB",
    "META_MOBILEOG_DB",
    "META_MOBILEOG_META",
)

DEFAULT_CONDA_ENV_NAMES = (
    "base",
    "genomad_aux",
    "amr_aux",
    "host_filter",
    "mag_aux",
    "cm210",
    "sistr_hicap",
    "qiime2",
    "microeco",
    "genovi",
    "plasflow",
    "longread_aux",
    "chewie",
    "ncov",
    "choleraefinder",
    "report_env",
    "Vlib",
    "tNGS",
)

FULL_ANALYSIS_CONDA_ENV_NAMES = (
    "web_runtime",
    "meta_main",
    "genomad_aux",
    "amr_aux",
    "host_filter",
    "mag_aux",
    "cm210",
    "sistr_hicap",
    "qiime2",
    "microeco",
    "genovi",
    "plasflow",
    "longread_aux",
    "chewie",
    "ncov",
    "choleraefinder",
    "report_env",
    "Vlib",
    "tNGS",
)

FULL_ANALYSIS_ENV_BINARIES = (
    ("web_runtime", "python"),
    ("web_runtime", "gunicorn"),
    ("meta_main", "python"),
    ("meta_main", "fastp"),
    ("meta_main", "kraken2"),
    ("meta_main", "spades.py"),
    ("genomad_aux", "virsorter"),
    ("genomad_aux", "checkv"),
    ("qiime2", "python"),
    ("microeco", "R"),
    ("report_env", "R"),
)


@dataclass(frozen=True)
class Check:
    name: str
    status: str
    message: str
    detail: str = ""

    def as_dict(self) -> dict[str, str]:
        return {
            "name": self.name,
            "status": self.status,
            "message": self.message,
            "detail": self.detail,
        }


def parse_env_file(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    if not path.is_file():
        return values
    for line_number, raw_line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("export "):
            line = line[len("export "):].strip()
        if "=" not in line:
            values[f"__INVALID_LINE_{line_number}"] = raw_line
            continue
        key, value = line.split("=", 1)
        key = key.strip()
        value = value.strip().strip("'").strip('"')
        if key:
            values[key] = value
    return values


def merged_env(env_file: Path | None) -> dict[str, str]:
    values = dict(os.environ)
    if env_file is not None:
        values.update(parse_env_file(env_file))
    return values


def has_placeholder(value: str) -> bool:
    text = str(value or "").strip()
    return not text or any(part in text for part in PLACEHOLDER_PARTS)


def resolve_path(raw_value: str, project_root: Path) -> Path:
    path = Path(str(raw_value or "").strip()).expanduser()
    return path if path.is_absolute() else project_root / path


def check_writable_dir(path: Path) -> tuple[bool, str]:
    try:
        path.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(prefix=".deploy-check-", dir=path, delete=True) as handle:
            handle.write(b"ok")
            handle.flush()
        return True, "writable"
    except OSError as exc:
        return False, str(exc)


def list_conda_envs(conda_root: Path | None) -> set[str]:
    envs: set[str] = set()
    if conda_root is None:
        return envs
    if (conda_root / "bin" / "python").is_file() or (conda_root / "python.exe").is_file():
        envs.add("base")
    envs_dir = conda_root / "envs"
    if envs_dir.is_dir():
        envs.update(child.name for child in envs_dir.iterdir() if child.is_dir())
    return envs


def detect_conda_root(values: dict[str, str]) -> Path | None:
    raw_root = values.get("CONDA_ROOT") or values.get("META_CONDA_ROOT") or values.get("PORTAL_CONDA_ROOT") or ""
    if raw_root.strip():
        return Path(raw_root).expanduser()
    conda_exe = values.get("CONDA_EXE") or shutil.which("conda") or ""
    if conda_exe:
        candidate = Path(conda_exe).expanduser()
        try:
            resolved = candidate.resolve()
        except OSError:
            resolved = candidate
        if resolved.name in {"conda", "conda.exe"}:
            parent = resolved.parent
            if parent.name in {"bin", "condabin", "Scripts"}:
                return parent.parent
    return None


def requested_conda_envs(values: dict[str, str], extra_envs: Iterable[str]) -> list[str]:
    envs = {item.strip() for item in extra_envs if item.strip()}
    for key in ("CONDA_ENV", "PIPELINE_CONDA_ENV"):
        if values.get(key, "").strip():
            envs.add(values[key].strip())
    raw = values.get("PORTAL_REQUIRED_CONDA_ENVS", "")
    if raw.strip():
        envs.update(item.strip() for item in raw.split(",") if item.strip())
    return sorted(envs)


def disk_free_gb(path: Path) -> float:
    target = path
    while not target.exists() and target.parent != target:
        target = target.parent
    usage = shutil.disk_usage(target)
    return usage.free / (1024 ** 3)


def build_checks(args: argparse.Namespace) -> list[Check]:
    project_root = Path(args.project_root).expanduser().resolve()
    env_file = Path(args.env_file).expanduser().resolve() if args.env_file else None
    values = merged_env(env_file)
    full_analysis = args.profile == "full-analysis"
    strict_assets = args.strict_assets or full_analysis
    checks: list[Check] = []

    if env_file is None:
        checks.append(Check("env_file", "warn", "未指定 env 文件，将只检查当前进程环境。"))
    elif env_file.is_file():
        checks.append(Check("env_file", "ok", f"已读取环境文件: {env_file}"))
    else:
        checks.append(Check("env_file", "fail", f"环境文件不存在: {env_file}"))

    invalid_lines = sorted(key for key in values if key.startswith("__INVALID_LINE_"))
    if invalid_lines:
        checks.append(Check("env_syntax", "fail", "环境文件存在无法解析的行。", ", ".join(invalid_lines)))
    else:
        checks.append(Check("env_syntax", "ok", "环境文件语法可解析。"))

    version = sys.version_info
    if version >= (3, 10):
        checks.append(Check("python_version", "ok", f"Python {version.major}.{version.minor}.{version.micro}"))
    else:
        checks.append(Check("python_version", "fail", f"Python 版本过低: {version.major}.{version.minor}.{version.micro}", "需要 Python 3.10+"))

    for relative_path in ("requirements-web.txt", "bac_analysis_portal/app.py", "bac_analysis_portal/application.py"):
        target = project_root / relative_path
        checks.append(
            Check(
                f"source:{relative_path}",
                "ok" if target.is_file() else "fail",
                f"源码文件存在: {relative_path}" if target.is_file() else f"源码文件缺失: {relative_path}",
            )
        )

    portal_mode = values.get("PORTAL_MODE", "").strip()
    if portal_mode == "production":
        for key in REQUIRED_PRODUCTION_ENV:
            value = values.get(key, "")
            if has_placeholder(value):
                checks.append(Check(f"env:{key}", "fail", f"{key} 未配置有效生产值。"))
            else:
                checks.append(Check(f"env:{key}", "ok", f"{key} 已配置。"))
        if values.get("PORTAL_INITIAL_ADMIN_PASSWORD", "") == "admin123":
            checks.append(Check("env:admin_password", "fail", "生产环境不能使用演示管理员密码 admin123。"))
    elif portal_mode:
        checks.append(Check("env:PORTAL_MODE", "warn", f"当前不是生产模式: {portal_mode}"))
    else:
        checks.append(Check("env:PORTAL_MODE", "fail", "PORTAL_MODE 未配置。"))

    cookie_secure = values.get("PORTAL_COOKIE_SECURE", "").strip()
    if cookie_secure in {"0", "1", "true", "false", "True", "False"}:
        checks.append(Check("env:PORTAL_COOKIE_SECURE", "ok", f"PORTAL_COOKIE_SECURE={cookie_secure}"))
    else:
        checks.append(Check("env:PORTAL_COOKIE_SECURE", "warn", "PORTAL_COOKIE_SECURE 未明确配置为 0 或 1。"))

    db_path = resolve_path(values.get("PORTAL_DB_PATH", "bac_analysis_portal.sqlite3"), project_root)
    ok, detail = check_writable_dir(db_path.parent)
    checks.append(Check("path:portal_db_parent", "ok" if ok else "fail", f"数据库父目录: {db_path.parent}", detail))

    task_root = resolve_path(values.get("BAC_ANALYSIS_TASK_ROOT", "analysis_tasks"), project_root)
    ok, detail = check_writable_dir(task_root)
    checks.append(Check("path:task_root", "ok" if ok else "fail", f"任务目录: {task_root}", detail))

    for name, path in (("portal_db_parent", db_path.parent), ("task_root", task_root)):
        try:
            free_gb = disk_free_gb(path)
        except OSError as exc:
            checks.append(Check(f"disk:{name}", "fail", f"无法读取磁盘空间: {path}", str(exc)))
            continue
        status = "ok" if free_gb >= args.min_free_gb else "fail"
        checks.append(Check(f"disk:{name}", status, f"{name} 可用空间 {free_gb:.1f} GB", f"最低要求 {args.min_free_gb:.1f} GB"))

    database_root_raw = values.get("META_DATABASE_ROOT", "")
    if database_root_raw.strip():
        database_root = resolve_path(database_root_raw, project_root)
        if database_root.is_dir():
            readable = os.access(database_root, os.R_OK)
            checks.append(Check("path:META_DATABASE_ROOT", "ok" if readable else "fail", f"数据库目录: {database_root}", "readable" if readable else "not readable"))
        else:
            checks.append(Check("path:META_DATABASE_ROOT", "fail" if strict_assets else "warn", f"数据库目录不存在: {database_root}"))
    else:
        checks.append(Check("path:META_DATABASE_ROOT", "fail" if full_analysis else "warn", "META_DATABASE_ROOT 未配置；Web-only 部署可以暂时接受。"))

    for key in OPTIONAL_ASSET_PATHS:
        raw_value = values.get(key, "")
        if not raw_value.strip():
            checks.append(Check(f"asset:{key}", "fail" if full_analysis else "warn", f"{key} 未配置。"))
            continue
        target = resolve_path(raw_value, project_root)
        exists = target.exists()
        checks.append(
            Check(
                f"asset:{key}",
                "ok" if exists else ("fail" if strict_assets else "warn"),
                f"{key}: {target}" if exists else f"{key} 指向的资产不存在: {target}",
            )
        )

    conda_root = detect_conda_root(values)
    if conda_root is None:
        checks.append(Check("conda:root", "fail" if full_analysis else "warn", "未检测到 Conda；Web-only 部署可以运行，真实分析不可运行。"))
    elif conda_root.is_dir():
        checks.append(Check("conda:root", "ok", f"Conda 根目录: {conda_root}"))
        conda_exe = conda_root / "bin" / "conda"
        checks.append(
            Check(
                "conda:exe",
                "ok" if conda_exe.is_file() and os.access(conda_exe, os.X_OK) else ("fail" if full_analysis else "warn"),
                f"Conda 可执行文件: {conda_exe}",
            )
        )
        conda_envs = list_conda_envs(conda_root)
        if conda_envs:
            checks.append(Check("conda:env_scan", "ok", f"检测到 {len(conda_envs)} 个 Conda 环境。"))
        else:
            checks.append(Check("conda:env_scan", "fail" if full_analysis else "warn", "Conda 根目录存在，但未扫描到 envs。"))
        required_envs = requested_conda_envs(values, args.conda_env)
        if args.check_default_conda_envs:
            required_envs = sorted(set(required_envs).union(DEFAULT_CONDA_ENV_NAMES))
        if full_analysis:
            required_envs = sorted(set(required_envs).union(FULL_ANALYSIS_CONDA_ENV_NAMES))
        for env_name in required_envs:
            checks.append(
                Check(
                    f"conda:env:{env_name}",
                    "ok" if env_name in conda_envs else "fail",
                    f"Conda 环境存在: {env_name}" if env_name in conda_envs else f"Conda 环境缺失: {env_name}",
                )
            )
        if full_analysis:
            for env_name, binary_name in FULL_ANALYSIS_ENV_BINARIES:
                target = conda_root / "envs" / env_name / "bin" / binary_name
                checks.append(
                    Check(
                        f"conda:bin:{env_name}:{binary_name}",
                        "ok" if target.is_file() and os.access(target, os.X_OK) else "fail",
                        f"关键工具存在: {target}" if target.exists() else f"关键工具缺失: {target}",
                    )
                )
    else:
        checks.append(Check("conda:root", "fail", f"Conda 根目录不存在: {conda_root}"))

    compose_file = project_root / "deployment" / "docker-compose.portal.yml"
    dockerfile = project_root / "deployment" / "Dockerfile.portal"
    checks.append(Check("docker:compose_file", "ok" if compose_file.is_file() else "fail", f"Compose 文件: {compose_file}"))
    checks.append(Check("docker:dockerfile", "ok" if dockerfile.is_file() else "fail", f"Dockerfile: {dockerfile}"))

    return checks


def summarize(checks: list[Check]) -> dict[str, object]:
    counts = {"ok": 0, "warn": 0, "fail": 0}
    for check in checks:
        counts[check.status] = counts.get(check.status, 0) + 1
    return {
        "status": "fail" if counts.get("fail", 0) else "warn" if counts.get("warn", 0) else "ok",
        "counts": counts,
        "checks": [check.as_dict() for check in checks],
    }


def print_text_report(result: dict[str, object]) -> None:
    status = str(result["status"])
    counts = result["counts"]
    print(f"Deployment check: {status.upper()}  ok={counts['ok']} warn={counts['warn']} fail={counts['fail']}")
    print()
    for item in result["checks"]:
        marker = {"ok": "OK", "warn": "WARN", "fail": "FAIL"}[item["status"]]
        line = f"[{marker}] {item['name']}: {item['message']}"
        if item.get("detail"):
            line += f" ({item['detail']})"
        print(line)


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Check Pathogen Workbench deployment readiness.")
    parser.add_argument("--project-root", default=str(Path(__file__).resolve().parents[1]), help="Project root to inspect.")
    parser.add_argument("--env-file", default="deployment/portal.env", help="Environment file to inspect.")
    parser.add_argument("--profile", choices=("web", "full-analysis"), default="web", help="Deployment profile to validate.")
    parser.add_argument("--json", action="store_true", help="Print machine-readable JSON.")
    parser.add_argument("--strict-assets", action="store_true", help="Treat missing optional database assets as failures.")
    parser.add_argument("--min-free-gb", type=float, default=5.0, help="Minimum free disk space for state and task paths.")
    parser.add_argument("--conda-env", action="append", default=[], help="Required Conda environment name. May be repeated.")
    parser.add_argument("--check-default-conda-envs", action="store_true", help="Check the default full analysis Conda environment set.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv or sys.argv[1:])
    checks = build_checks(args)
    result = summarize(checks)
    if args.json:
        print(json.dumps(result, ensure_ascii=False, indent=2))
    else:
        print_text_report(result)
    return 1 if result["status"] == "fail" else 0


if __name__ == "__main__":
    raise SystemExit(main())
