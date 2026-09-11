from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Callable

from .admin_runtime import _load_conda_env_settings, _load_conda_root_setting
from .dataset_operations import DatasetOperations
from .identity import UserIdentity
from .sample_library_manager import SampleLibraryManager
from .store import PortalStore
from .task_manager import AnalysisTaskManager, ValidationError, read_json


@dataclass(frozen=True)
class DatasetService:
    project_root: Path
    store: PortalStore
    task_manager: AnalysisTaskManager
    sample_manager: SampleLibraryManager
    operations: DatasetOperations
    nextstrain_build_task_type: str
    can_view_task: Callable[[dict], bool]
    can_view_uploaded_auspice: Callable[[dict], bool]
    ensure_can_view_uploaded_auspice: Callable[[dict], None]
    load_uploaded_auspice_metadata: Callable[[Path, str], dict]
    uploaded_auspice_dataset_path: Callable[[Path, str], Path]
    open_in_file_manager: Callable[[Path], None]

    def list_nextstrain_builds(self) -> list[dict]:
        self.task_manager.reconcile_queue(self._max_concurrent_tasks())
        items = []
        for task_file in sorted(self.task_manager.task_root.glob("*/task.json"), key=lambda item: item.parent.name, reverse=True):
            try:
                task = read_json(task_file)
            except Exception:
                continue
            if str(task.get("task_type") or "").strip() == self.nextstrain_build_task_type and self.can_view_task(task):
                items.append(self.operations.serialize_nextstrain_build_task(self.project_root, task))
        return items

    def list_auspice_uploads(self) -> list[dict]:
        return [
            self.operations.serialize_uploaded_auspice_dataset(self.project_root, metadata)
            for metadata in self.operations.list_uploaded_auspice_metadata(self.project_root)
            if self.can_view_uploaded_auspice(metadata)
        ]

    def create_auspice_upload(
        self, *, raw_bytes: bytes, original_filename: str, name: str, identity: UserIdentity
    ) -> dict:
        if not original_filename.lower().endswith(".json"):
            raise ValidationError("目前只支持上传 .json 格式的 Auspice 数据集")
        try:
            payload = json.loads(raw_bytes.decode("utf-8"))
        except Exception as exc:
            raise ValidationError("上传文件不是有效的 JSON") from exc
        if not raw_bytes or not isinstance(payload, dict):
            raise ValidationError("Auspice 数据集必须是非空 JSON 对象")
        metadata = self.operations.store_uploaded_auspice(
            self.project_root, raw_bytes, original_filename, name, identity.username, identity.group_name
        )
        return self.operations.serialize_uploaded_auspice_dataset(self.project_root, metadata)

    def open_auspice_upload(self, upload_id: str) -> Path:
        metadata = self.load_uploaded_auspice_metadata(self.project_root, upload_id)
        self.ensure_can_view_uploaded_auspice(metadata)
        target = self.uploaded_auspice_dataset_path(self.project_root, upload_id).parent
        if not target.exists():
            raise ValidationError(f"上传目录不存在: {target}")
        self.open_in_file_manager(target)
        return target

    def delete_auspice_upload(self, upload_id: str) -> None:
        metadata = self.load_uploaded_auspice_metadata(self.project_root, upload_id)
        self.ensure_can_view_uploaded_auspice(metadata)
        self.operations.delete_uploaded_auspice(self.project_root, upload_id)

    def create_nextstrain_build(self, payload: dict, *, identity: UserIdentity) -> tuple[dict, dict]:
        sample_keys = payload.get("sample_keys") if isinstance(payload.get("sample_keys"), list) else []
        if len(sample_keys) < 2:
            raise ValidationError("请至少勾选 2 个样本后再构建 Nextstrain 数据集")
        reference_genbank = Path(str(payload.get("reference_genbank_path") or "").strip()).expanduser()
        if not reference_genbank.is_file():
            raise ValidationError("参考 GenBank 文件不存在，请重新选择")
        if not self.operations.nextstrain_cli_path.is_file():
            raise ValidationError(f"Nextstrain CLI 不存在: {self.operations.nextstrain_cli_path}")
        records, species_candidates = [], []
        for raw_key in sample_keys:
            try:
                record = self.sample_manager.get_visible(
                    str(raw_key or "").strip(),
                    role=identity.role,
                    username=identity.username,
                    group_name=identity.group_name,
                )
            except KeyError as exc:
                raise ValidationError(f"样本不可见或不存在: {raw_key}") from exc
            if not Path(str(record.get("final_fasta_path") or "")).expanduser().is_file():
                raise ValidationError(f"样本缺少可读取的 final.fasta: {record.get('sample_name') or raw_key}")
            records.append(record)
            species = str(record.get("species_name") or record.get("mlst_species_name") or "").strip()
            if species and species not in species_candidates:
                species_candidates.append(species)
        if len(records) < 2 or len(species_candidates) > 1:
            raise ValidationError("当前样本不足 2 个或混入多个物种，无法构建时间树")
        timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
        build_name = str(payload.get("build_name") or "").strip() or f"{species_candidates[0] if species_candidates else 'nextstrain'}_{timestamp}"
        workspace_name = self.operations.sanitize_slug(build_name) or f"nextstrain_{timestamp}"
        manifest = {
            "build_name": build_name,
            "species_name": species_candidates[0] if species_candidates else "",
            "created_by": identity.username,
            "samples": [
                {
                    **{
                        key: str(record.get(key) or "").strip()
                        for key in (
                            "sample_key", "sample_name", "final_fasta_path", "genome_id", "sample_alias", "country",
                            "collection_date", "owner", "host_info", "note",
                        )
                    },
                    "species_name": str(record.get("species_name") or record.get("mlst_species_name") or "").strip(),
                    "location": self.operations.parse_sample_location_json(record.get("location_json")),
                }
                for record in records
            ],
        }
        manifest_path, workspace_dir = self.operations.write_nextstrain_manifest(
            self.project_root, workspace_name, manifest
        )
        created = self.task_manager.create_task(
            {
                "task_name": build_name,
                "input_path": str(manifest_path),
                "output_dir": str(workspace_dir),
                "analysis_target": "virus",
                "inputtype": "fasta",
                "thread": 1,
                "method": "bwa",
                "asm_type": "shortref",
                "species": species_candidates[0] if species_candidates else "Virus",
                "ref": str(reference_genbank.resolve()),
                "runflow": "分型鉴定",
                "workstation_key": "virus",
            },
            owner=identity.username,
            owner_group=identity.group_name,
            pipeline_script=str((self.project_root / "scripts" / "run_nextstrain_auspice_build.py").resolve()),
            pipeline_python="base",
            max_concurrent_tasks=self._max_concurrent_tasks(),
            database_root=str(self.project_root / "database"),
            conda_root=_load_conda_root_setting(self.store),
            conda_envs=_load_conda_env_settings(self.store),
            extra_task_fields={"task_type": self.nextstrain_build_task_type, "nextstrain_cli": str(self.operations.nextstrain_cli_path)},
        )
        build = {
            "task_id": str(created.get("id") or "").strip(),
            "dataset_path": f"/nextstrain/{str(created.get('id') or '').strip()}",
            "workspace_dir": str(workspace_dir),
        }
        return created, build

    def _max_concurrent_tasks(self) -> int:
        return int(self.store.get_setting("max_concurrent_tasks", "2") or "2")
