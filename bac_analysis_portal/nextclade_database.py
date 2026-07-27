from __future__ import annotations

from contextlib import contextmanager
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import threading
from typing import Iterator
from uuid import uuid4

from .admin_runtime import _detect_conda_root
from .runtime_paths import _normalize_database_root_path
from .store import PortalStore
from .task_manager import ValidationError


NEXTCLADE_ENV_NAME = "ncov"
NEXTCLADE_DATASET_ROOTS = (
    Path("virus") / "nextclade_db",
    Path("nextclade_db"),
    Path("virus") / "nextclade",
)
NEXTCLADE_LIST_TIMEOUT_SECONDS = 120
NEXTCLADE_DOWNLOAD_TIMEOUT_SECONDS = 900
_NEXTCLADE_UPDATE_LOCK = threading.Lock()


def _normalized_text(value: object) -> str:
    return " ".join(str(value or "").strip().casefold().split())


def _configured_database_root(project_root: Path, store: PortalStore) -> Path:
    raw_value = str(store.get_setting("database_root", str(project_root)) or "").strip()
    candidate = Path(raw_value).expanduser() if raw_value else project_root
    try:
        resolved = _normalize_database_root_path(candidate)
    except OSError:
        resolved = _normalize_database_root_path(project_root)
    if resolved.is_dir():
        return resolved
    fallback = _normalize_database_root_path(project_root)
    if fallback.is_dir():
        return fallback
    raise ValidationError(f"数据库部署目录不存在: {resolved}")


def _resolve_dataset_root(database_root: Path) -> Path:
    for relative in NEXTCLADE_DATASET_ROOTS:
        candidate = database_root / relative
        try:
            if candidate.is_dir() and any(candidate.glob("*/pathogen.json")):
                return candidate.resolve()
        except OSError:
            continue
    searched = "、".join(str(database_root / item) for item in NEXTCLADE_DATASET_ROOTS)
    raise ValidationError(f"未找到 Nextclade 数据库目录，已检查: {searched}")


def _standard_conda_roots() -> tuple[Path, ...]:
    home = Path.home()
    return (
        home / "miniconda3",
        home / "anaconda3",
        home / "miniforge3",
        home / "mambaforge",
        Path("/opt/miniconda3"),
        Path("/opt/homebrew/Caskroom/miniconda/base"),
        Path("/opt/homebrew/Caskroom/mambaforge/base"),
    )


def _resolve_ncov_nextclade(store: PortalStore) -> Path:
    candidates: list[Path] = []
    configured = str(store.get_setting("conda_root", "") or "").strip()
    if configured:
        candidates.append(Path(configured).expanduser())
    detected = _detect_conda_root()
    if detected is not None:
        candidates.append(detected)
    candidates.extend(_standard_conda_roots())

    seen: set[str] = set()
    for root in candidates:
        key = str(root)
        if key in seen:
            continue
        seen.add(key)
        executables = (
            root / "envs" / NEXTCLADE_ENV_NAME / "bin" / "nextclade",
            root / "envs" / NEXTCLADE_ENV_NAME / "Scripts" / "nextclade.exe",
        )
        for executable in executables:
            if executable.is_file() and os.access(executable, os.X_OK):
                return executable.resolve()
    raise ValidationError("未找到 ncov Conda 环境中的 Nextclade，请先在后台管理配置正确的 Conda 安装路径")


def _run_nextclade(executable: Path, arguments: list[str], *, timeout: int) -> subprocess.CompletedProcess[str]:
    try:
        completed = subprocess.run(
            [str(executable), *arguments],
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        raise ValidationError("Nextclade 操作超时，请检查网络连接后重试") from exc
    except OSError as exc:
        raise ValidationError(f"无法运行 ncov 环境中的 Nextclade: {exc}") from exc
    if completed.returncode != 0:
        detail = (completed.stderr or completed.stdout or "Nextclade 返回未知错误").strip()
        raise ValidationError(f"Nextclade 操作失败: {detail}")
    return completed


def _read_json_object(path: Path) -> dict[str, object]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValidationError(f"无法读取数据集元数据 {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ValidationError(f"数据集元数据格式无效: {path}")
    return payload


def _load_remote_datasets(executable: Path) -> list[dict[str, object]]:
    completed = _run_nextclade(
        executable,
        ["dataset", "list", "--json"],
        timeout=NEXTCLADE_LIST_TIMEOUT_SECONDS,
    )
    try:
        payload = json.loads(completed.stdout)
    except json.JSONDecodeError as exc:
        raise ValidationError("Nextclade 返回的数据集列表不是有效 JSON") from exc
    if not isinstance(payload, list):
        raise ValidationError("Nextclade 返回的数据集列表格式无效")
    return [item for item in payload if isinstance(item, dict)]


def _metadata_signature(payload: dict[str, object]) -> tuple[str, str, str]:
    attributes = payload.get("attributes")
    if not isinstance(attributes, dict):
        attributes = {}
    return (
        _normalized_text(attributes.get("name")),
        _normalized_text(attributes.get("reference accession")),
        _normalized_text(attributes.get("reference name")),
    )


def _dataset_aliases(payload: dict[str, object]) -> set[str]:
    aliases = {_normalized_text(payload.get("path"))}
    shortcuts = payload.get("shortcuts")
    if isinstance(shortcuts, list):
        aliases.update(_normalized_text(item) for item in shortcuts)
    aliases.discard("")
    return aliases


def _match_remote_dataset(
    local_payload: dict[str, object], remote_datasets: list[dict[str, object]]
) -> dict[str, object] | None:
    local_aliases = _dataset_aliases(local_payload)
    alias_matches = [item for item in remote_datasets if local_aliases & _dataset_aliases(item)]
    if len(alias_matches) == 1:
        return alias_matches[0]

    local_signature = _metadata_signature(local_payload)
    signature_matches = [
        item
        for item in remote_datasets
        if _metadata_signature(item) == local_signature and all(local_signature)
    ]
    if len(signature_matches) == 1:
        return signature_matches[0]
    return None


def _version_tag(payload: dict[str, object]) -> str:
    version = payload.get("version")
    return str(version.get("tag") or "").strip() if isinstance(version, dict) else ""


def _display_name(payload: dict[str, object], fallback: str) -> str:
    attributes = payload.get("attributes")
    if isinstance(attributes, dict):
        name = str(attributes.get("name") or "").strip()
        if name:
            return name
    return fallback


def _summarize(items: list[dict[str, object]]) -> dict[str, int]:
    return {
        "total_count": len(items),
        "latest_count": sum(item.get("status") == "latest" for item in items),
        "outdated_count": sum(item.get("status") == "outdated" for item in items),
        "unmatched_count": sum(item.get("status") == "unmatched" for item in items),
        "failed_count": sum(item.get("status") in {"invalid", "update_failed"} for item in items),
    }


def _scan_local_datasets(dataset_root: Path, remote_datasets: list[dict[str, object]]) -> list[dict[str, object]]:
    items: list[dict[str, object]] = []
    for metadata_path in sorted(dataset_root.glob("*/pathogen.json"), key=lambda item: item.parent.name.casefold()):
        directory = metadata_path.parent.name
        try:
            local_payload = _read_json_object(metadata_path)
        except ValidationError as exc:
            items.append(
                {
                    "directory": directory,
                    "display_name": directory,
                    "local_tag": "",
                    "latest_tag": "",
                    "dataset_path": "",
                    "status": "invalid",
                    "updatable": False,
                    "message": str(exc),
                }
            )
            continue
        remote = _match_remote_dataset(local_payload, remote_datasets)
        local_tag = _version_tag(local_payload)
        if remote is None:
            items.append(
                {
                    "directory": directory,
                    "display_name": _display_name(local_payload, directory),
                    "local_tag": local_tag,
                    "latest_tag": "",
                    "dataset_path": "",
                    "status": "unmatched",
                    "updatable": False,
                    "message": "无法在当前 Nextclade 在线列表中唯一匹配该数据集",
                }
            )
            continue
        latest_tag = _version_tag(remote)
        is_latest = bool(local_tag and latest_tag and local_tag == latest_tag)
        items.append(
            {
                "directory": directory,
                "display_name": _display_name(local_payload, directory),
                "local_tag": local_tag,
                "latest_tag": latest_tag,
                "dataset_path": str(remote.get("path") or "").strip(),
                "status": "latest" if is_latest else "outdated",
                "updatable": bool(not is_latest and latest_tag and remote.get("path")),
                "message": "已是最新兼容版本" if is_latest else "存在可更新的兼容版本",
            }
        )
    if not items:
        raise ValidationError(f"Nextclade 数据库目录中没有数据集: {dataset_root}")
    return items


def _nextclade_version(executable: Path) -> str:
    completed = _run_nextclade(executable, ["--version"], timeout=30)
    return completed.stdout.strip() or "未知"


@contextmanager
def _update_guard() -> Iterator[None]:
    if not _NEXTCLADE_UPDATE_LOCK.acquire(blocking=False):
        raise ValidationError("已有 Nextclade 数据库更新正在执行，请稍后重试")
    try:
        yield
    finally:
        _NEXTCLADE_UPDATE_LOCK.release()


def _download_and_replace(
    executable: Path,
    dataset_root: Path,
    item: dict[str, object],
) -> None:
    directory = str(item.get("directory") or "").strip()
    dataset_path = str(item.get("dataset_path") or "").strip()
    expected_tag = str(item.get("latest_tag") or "").strip()
    target = dataset_root / directory
    if not directory or target.parent.resolve() != dataset_root.resolve() or not target.is_dir() or target.is_symlink():
        raise ValidationError(f"数据集目录无效: {directory or '-'}")
    if not dataset_path or not expected_tag:
        raise ValidationError(f"数据集 {directory} 缺少在线路径或版本信息")

    staging_root = Path(tempfile.mkdtemp(prefix=".nextclade_update_", dir=dataset_root))
    download_root = staging_root / "dataset"
    backup = dataset_root / f".{directory}.backup-{uuid4().hex}"
    try:
        _run_nextclade(
            executable,
            [
                "dataset",
                "get",
                "--name",
                dataset_path,
                "--tag",
                expected_tag,
                "--output-dir",
                str(download_root),
            ],
            timeout=NEXTCLADE_DOWNLOAD_TIMEOUT_SECONDS,
        )
        downloaded_metadata = _read_json_object(download_root / "pathogen.json")
        if _version_tag(downloaded_metadata) != expected_tag:
            raise ValidationError(f"数据集 {directory} 下载版本校验失败")

        target.rename(backup)
        try:
            download_root.rename(target)
        except OSError:
            backup.rename(target)
            raise
        shutil.rmtree(backup, ignore_errors=True)
    except OSError as exc:
        if backup.is_dir() and not target.exists():
            backup.rename(target)
        raise ValidationError(f"替换数据集 {directory} 失败: {exc}") from exc
    finally:
        shutil.rmtree(staging_root, ignore_errors=True)


class NextcladeDatabaseManager:
    def __init__(self, *, project_root: Path, store: PortalStore) -> None:
        self.project_root = project_root
        self.store = store

    def check(self) -> dict[str, object]:
        database_root = _configured_database_root(self.project_root, self.store)
        dataset_root = _resolve_dataset_root(database_root)
        executable = _resolve_ncov_nextclade(self.store)
        remote_datasets = _load_remote_datasets(executable)
        items = _scan_local_datasets(dataset_root, remote_datasets)
        return {
            "checked_at": datetime.now(timezone.utc).isoformat(),
            "database_root": str(dataset_root),
            "environment": NEXTCLADE_ENV_NAME,
            "nextclade_version": _nextclade_version(executable),
            "items": items,
            "summary": _summarize(items),
        }

    def update(self, directories: object = None) -> dict[str, object]:
        with _update_guard():
            status = self.check()
            items = status["items"]
            assert isinstance(items, list)
            requested: set[str] | None = None
            if directories is not None:
                if not isinstance(directories, list) or not all(isinstance(item, str) for item in directories):
                    raise ValidationError("待更新数据集列表格式不正确")
                requested = {item.strip() for item in directories if item.strip()}
                known = {str(item.get("directory") or "") for item in items}
                unknown = sorted(requested - known)
                if unknown:
                    raise ValidationError(f"包含未知数据集目录: {', '.join(unknown)}")

            executable = _resolve_ncov_nextclade(self.store)
            dataset_root = Path(str(status["database_root"]))
            update_results: list[dict[str, object]] = []
            for item in items:
                directory = str(item.get("directory") or "")
                if item.get("status") != "outdated" or (requested is not None and directory not in requested):
                    continue
                try:
                    _download_and_replace(executable, dataset_root, item)
                except ValidationError as exc:
                    item["status"] = "update_failed"
                    item["updatable"] = True
                    item["message"] = str(exc)
                    update_results.append({"directory": directory, "status": "failed", "message": str(exc)})
                    continue
                item["local_tag"] = item.get("latest_tag")
                item["status"] = "latest"
                item["updatable"] = False
                item["message"] = "更新完成"
                update_results.append({"directory": directory, "status": "updated", "message": "更新完成"})

            if not update_results:
                raise ValidationError("没有可更新的 Nextclade 数据集")
            status["updated_at"] = datetime.now(timezone.utc).isoformat()
            status["update_results"] = update_results
            status["summary"] = _summarize(items)
            status["updated_count"] = sum(item["status"] == "updated" for item in update_results)
            status["update_failed_count"] = sum(item["status"] == "failed" for item in update_results)
            return status
