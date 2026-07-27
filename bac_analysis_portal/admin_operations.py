from __future__ import annotations

from dataclasses import dataclass
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from typing import Callable

from .task_manager import ValidationError


UPDATE_REPO_URL = "https://github.com/muguading/metagenomic.git"
UPDATE_EXCLUDE_DIRS = {".git", "__pycache__", ".pytest_cache", ".mypy_cache", ".ruff_cache", ".venv", ".venv_web", ".venv_genome_db_test", "analysis_tasks", "generated_batch_inputs", "node_modules"}
UPDATE_EXCLUDE_FILES = {"bac_analysis_portal.sqlite3", ".DS_Store"}


def resolve_admin_root(root_arg: str) -> Path:
    candidate = Path(str(root_arg or "").strip()).expanduser().resolve() if str(root_arg or "").strip() else Path.home().resolve()
    if not candidate.is_dir():
        raise ValidationError(f"目录不存在: {candidate}")
    return candidate


def resolve_admin_path(root: Path, path_arg: str) -> Path:
    raw = str(path_arg or "").strip()
    candidate = Path(raw).expanduser() if raw else root
    return candidate.resolve() if candidate.is_absolute() else (root / candidate).resolve()


def _is_update_candidate_root(path: Path) -> bool:
    try:
        return path.is_dir() and (path / "Bac_assemble_260112_newformat.py").is_file() and (path / "bac_analysis_portal").is_dir()
    except OSError:
        return False


def _collect_update_files(root: Path) -> list[str]:
    files = []
    for current_root, dirnames, filenames in os.walk(root):
        dirnames[:] = [name for name in dirnames if name not in UPDATE_EXCLUDE_DIRS]
        files.extend(str((Path(current_root) / name).relative_to(root)) for name in filenames if name not in UPDATE_EXCLUDE_FILES and not name.endswith((".pyc", ".pyo")))
    return files


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _calculate_update_diff(source_root: Path, target_root: Path) -> dict:
    added, changed, skipped, copy_files = 0, 0, 0, []
    for relative in sorted(_collect_update_files(source_root)):
        source, target = source_root / relative, target_root / relative
        if not target.exists():
            added += 1
            copy_files.append(relative)
        elif _file_sha256(source) != _file_sha256(target):
            changed += 1
            copy_files.append(relative)
        else:
            skipped += 1
    return {"added_files": added, "changed_files": changed, "skipped_files": skipped, "copy_files": copy_files}


def _apply_update_tree(source_root: Path, target_root: Path) -> dict:
    if not target_root.is_dir():
        raise ValidationError(f"部署目录不存在: {target_root}")
    diff = _calculate_update_diff(source_root, target_root)
    for relative in diff["copy_files"]:
        target = target_root / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source_root / relative, target)
    return {**{key: diff[key] for key in ("added_files", "changed_files", "skipped_files")}, "target_path": str(target_root), "restart_required": True}


def detect_offline_update_sources(workspace_root: Path) -> list[dict]:
    roots = []
    for base in (Path("/Volumes"), Path("/media"), Path("/run/media"), Path("/mnt")):
        if base.is_dir():
            try:
                roots.extend(child.resolve() for child in base.iterdir() if child.is_dir())
            except OSError:
                pass
    candidates = []
    for mount_root in roots:
        queue = [(mount_root, 0)]
        while queue:
            current, depth = queue.pop(0)
            if _is_update_candidate_root(current) and current != workspace_root:
                diff = _calculate_update_diff(current, workspace_root)
                candidates.append({"path": str(current), "label": f"{mount_root.name} · {current.name} · 新增{diff['added_files']} / 更新{diff['changed_files']}", "mount_root": str(mount_root), **{key: diff[key] for key in ("added_files", "changed_files", "skipped_files")}})
                continue
            if depth < 3:
                try:
                    queue.extend((child, depth + 1) for child in current.iterdir() if child.is_dir() and not child.name.startswith("."))
                except OSError:
                    pass
    return candidates


def run_online_update(workspace_root: Path) -> dict:
    with tempfile.TemporaryDirectory(prefix="metagenomic_update_") as tmpdir:
        clone_root = Path(tmpdir) / "repo"
        completed = subprocess.run(["git", "clone", "--depth", "1", UPDATE_REPO_URL, str(clone_root)], capture_output=True, text=True)
        if completed.returncode != 0 or not _is_update_candidate_root(clone_root):
            raise ValidationError(f"在线更新失败: {(completed.stderr or completed.stdout or '更新源无效').strip()}")
        return {**_apply_update_tree(clone_root, workspace_root), "mode": "online", "source_path": str(clone_root), "repo_url": UPDATE_REPO_URL}


def run_offline_update(workspace_root: Path, source_path: Path) -> dict:
    if not _is_update_candidate_root(source_path):
        raise ValidationError("离线更新源不包含有效的部署目录结构")
    return {**_apply_update_tree(source_path, workspace_root), "mode": "offline", "source_path": str(source_path)}

@dataclass(frozen=True)
class AdminOperations:
    load_pathosource_trigger_rules: Callable[[], dict[str, object]]
    save_pathosource_trigger_rules: Callable[[dict[str, object]], dict[str, object]]
    resolve_admin_root: Callable[[str], Path]
    resolve_admin_path: Callable[[Path, str], Path]
    detect_offline_update_sources: Callable[[Path], list[dict]]
    run_online_update: Callable[[Path], dict]
    run_offline_update: Callable[[Path, Path], dict]
