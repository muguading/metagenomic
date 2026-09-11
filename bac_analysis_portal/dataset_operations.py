from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Callable
from datetime import datetime
import uuid

NEXTSTRAIN_BUILD_TASK_TYPE = "nextstrain_auspice_build"
NEXTSTRAIN_CLI_PATH = Path("/Users/wuhhh/.nextstrain/cli-standalone/nextstrain")
AUSPICE_UPLOAD_KIND = "uploaded_auspice_json"


def is_portal_only_task(task: dict) -> bool:
    return str(task.get("task_type") or "").strip() == NEXTSTRAIN_BUILD_TASK_TYPE


def parse_sample_location_json(raw_value: object) -> dict[str, str]:
    if isinstance(raw_value, dict):
        payload = raw_value
    else:
        try:
            payload = json.loads(str(raw_value or "").strip())
        except (json.JSONDecodeError, TypeError):
            return {}
    return {str(key): str(value or "").strip() for key, value in payload.items()} if isinstance(payload, dict) else {}


def normalize_nextstrain_dataset_prefix(prefix: str) -> tuple[str, str] | None:
    parts = [part for part in str(prefix or "").strip().strip("/").split("/") if part]
    if len(parts) >= 2 and parts[0] in {"nextstrain", "uploaded"}:
        return ("task" if parts[0] == "nextstrain" else "uploaded", parts[1])
    if len(parts) >= 3 and parts[0] == "auspice" and parts[1] in {"nextstrain", "uploaded"}:
        return ("task" if parts[1] == "nextstrain" else "uploaded", parts[2])
    return None


def uploaded_auspice_root(project_root: Path) -> Path:
    return (project_root / "public" / "auspice-uploads").resolve()


def uploaded_auspice_metadata_path(project_root: Path, upload_id: str) -> Path:
    return uploaded_auspice_root(project_root) / "metadata" / f"{upload_id}.json"


def uploaded_auspice_dataset_path(project_root: Path, upload_id: str) -> Path:
    return uploaded_auspice_root(project_root) / "datasets" / f"{upload_id}.json"


def load_uploaded_auspice_metadata(project_root: Path, upload_id: str) -> dict:
    path = uploaded_auspice_metadata_path(project_root, upload_id)
    if not path.is_file():
        raise KeyError(f"Uploaded auspice dataset not found: {upload_id}")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise KeyError(f"Uploaded auspice dataset is invalid: {upload_id}")
    return payload


def resolve_nextstrain_dataset_main_json(project_root: Path, task: dict) -> Path:
    return Path(str((task.get("params") or {}).get("output_dir") or "")).expanduser().resolve() / "auspice" / "zika.json"


def resolve_nextstrain_dataset_sidecar(project_root: Path, task: dict, sidecar_type: str) -> Path:
    output_dir = Path(str((task.get("params") or {}).get("output_dir") or "")).expanduser().resolve() / "auspice"
    suffix = {"root-sequence": "_root-sequence", "tip-frequencies": "_tip-frequencies", "measurements": "_measurements"}.get(sidecar_type, "")
    return output_dir / f"zika{suffix}.json"


def serialize_nextstrain_build_task(project_root: Path, task: dict) -> dict:
    params = task.get("params") if isinstance(task.get("params"), dict) else {}
    dataset_json = resolve_nextstrain_dataset_main_json(project_root, task)
    manifest = {}
    manifest_path = Path(str(params.get("output_dir") or "")).expanduser().resolve() / "build_manifest.json"
    if manifest_path.is_file():
        try:
            loaded = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest = loaded if isinstance(loaded, dict) else {}
        except json.JSONDecodeError:
            pass
    return {"entry_kind": "build", "task_id": str(task.get("id") or "").strip(), "task_name": str(task.get("name") or "").strip(), "status": str(task.get("status") or "").strip().upper(), "created_at": task.get("created_at"), "started_at": task.get("started_at"), "finished_at": task.get("finished_at"), "workspace_dir": str(params.get("output_dir") or "").strip(), "reference_genbank": str(params.get("ref") or "").strip(), "sample_count": int(manifest.get("sample_count") or len(manifest.get("samples") or [])), "species_name": str(manifest.get("species_name") or "").strip(), "dataset_ready": dataset_json.is_file(), "dataset_path": f"/nextstrain/{str(task.get('id') or '').strip()}", "dataset_json_path": str(dataset_json), "log_path": str(task.get("log_path") or "").strip()}


def serialize_uploaded_auspice_dataset(project_root: Path, metadata: dict) -> dict:
    upload_id = str(metadata.get("id") or "").strip()
    dataset_json = uploaded_auspice_dataset_path(project_root, upload_id)
    return {"entry_kind": "upload", "upload_id": upload_id, "task_id": "", "task_name": str(metadata.get("name") or metadata.get("original_filename") or upload_id).strip(), "status": "UPLOADED", "created_at": metadata.get("created_at"), "started_at": metadata.get("created_at"), "finished_at": metadata.get("created_at"), "workspace_dir": "", "reference_genbank": "", "sample_count": int(metadata.get("sample_count") or 0), "species_name": str(metadata.get("species_name") or "").strip(), "dataset_ready": dataset_json.is_file(), "dataset_path": f"/uploaded/{upload_id}", "dataset_json_path": str(dataset_json), "log_path": "", "owner": str(metadata.get("owner") or "").strip(), "owner_group": str(metadata.get("owner_group") or "").strip(), "original_filename": str(metadata.get("original_filename") or "").strip(), "source_kind": AUSPICE_UPLOAD_KIND}


def list_uploaded_auspice_metadata(project_root: Path) -> list[dict]:
    metadata_dir = uploaded_auspice_root(project_root) / "metadata"
    items = []
    if not metadata_dir.is_dir():
        return items
    for metadata_file in sorted(metadata_dir.glob("*.json"), key=lambda item: item.name, reverse=True):
        try:
            metadata = json.loads(metadata_file.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        if isinstance(metadata, dict):
            items.append(metadata)
    return items


def store_uploaded_auspice(project_root: Path, raw_bytes: bytes, original_filename: str, name: str, owner: str, owner_group: str) -> dict:
    payload = json.loads(raw_bytes.decode("utf-8"))
    meta_payload = payload.get("meta") if isinstance(payload.get("meta"), dict) else {}
    inferred_name = str(meta_payload.get("title") or meta_payload.get("description") or "").strip()
    upload_id = f"auspice_upload_{datetime.now().strftime('%Y%m%d%H%M%S')}_{uuid.uuid4().hex[:8]}"
    upload_root = uploaded_auspice_root(project_root)
    metadata_dir, dataset_dir = upload_root / "metadata", upload_root / "datasets"
    metadata_dir.mkdir(parents=True, exist_ok=True)
    dataset_dir.mkdir(parents=True, exist_ok=True)
    (dataset_dir / f"{upload_id}.json").write_bytes(raw_bytes)
    metadata = {
        "id": upload_id,
        "name": name or inferred_name or Path(original_filename).stem,
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "owner": owner,
        "owner_group": owner_group,
        "original_filename": original_filename,
        "species_name": "",
        "sample_count": 0,
        "sha256": hashlib.sha256(raw_bytes).hexdigest(),
    }
    uploaded_auspice_metadata_path(project_root, upload_id).write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    return metadata


def delete_uploaded_auspice(project_root: Path, upload_id: str) -> None:
    for path in (uploaded_auspice_metadata_path(project_root, upload_id), uploaded_auspice_dataset_path(project_root, upload_id)):
        if path.exists():
            path.unlink()


def write_nextstrain_manifest(project_root: Path, workspace_name: str, manifest: dict) -> tuple[Path, Path]:
    target_root = (project_root / "public" / "nextstrain-builds").resolve()
    workspace_dir = target_root / workspace_name
    manifest_dir = target_root / "manifests"
    manifest_dir.mkdir(parents=True, exist_ok=True)
    workspace_dir.parent.mkdir(parents=True, exist_ok=True)
    manifest_path = manifest_dir / f"{workspace_dir.name}.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    return manifest_path, workspace_dir


@dataclass(frozen=True)
class DatasetOperations:
    nextstrain_cli_path: Path
    serialize_nextstrain_build_task: Callable[[Path, dict], dict]
    uploaded_auspice_root: Callable[[Path], Path]
    uploaded_auspice_metadata_path: Callable[[Path, str], Path]
    serialize_uploaded_auspice_dataset: Callable[[Path, dict], dict]
    parse_sample_location_json: Callable[[object], dict[str, str]]
    sanitize_slug: Callable[[str], str]
    list_uploaded_auspice_metadata: Callable[[Path], list[dict]]
    store_uploaded_auspice: Callable[[Path, bytes, str, str, str, str], dict]
    delete_uploaded_auspice: Callable[[Path, str], None]
    write_nextstrain_manifest: Callable[[Path, str, dict], tuple[Path, Path]]
