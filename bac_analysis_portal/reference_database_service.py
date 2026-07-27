from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .reference_database_ops import (
    _build_reference_database_index,
    _build_single_reference_index,
    _delete_reference_database_record,
    _queue_pathogen_cgmlst_panels,
    _save_uploaded_reference_genome,
    _save_uploaded_task_reference,
    _update_reference_database_record,
)
from .reference_operations import ReferenceOperations
from .store import PortalStore


@dataclass(frozen=True)
class ReferenceDatabaseService:
    store: PortalStore
    project_root: Path
    operations: ReferenceOperations

    def list_download_jobs(self, category: str) -> list[dict]:
        return self.operations.list_download_jobs(category)

    def delete_download_job(self, job_id: str, category: str) -> None:
        self.operations.delete_download_job(job_id, category)

    def start_download_job(self, *, category: str, payload: dict, owner: str) -> dict:
        return self.operations.start_download_job(category=category, payload=payload, owner=owner)

    def list_records(self, category: str) -> list[dict]:
        return self.store.list_host_database(category)

    def get_record(self, host_key: str) -> dict:
        return self.store.get_host_database_record(host_key)

    def list_panels(self, category: str = "pathogen", panel_type: str = "cgmlst") -> list[dict]:
        return self.store.list_reference_panels(category, panel_type)

    def import_from_source(self, *, category: str, owner: str, **kwargs) -> dict:
        return self.operations.import_from_source(
            store=self.store, project_root=self.project_root, category=category, owner=owner, **kwargs
        )

    def register_record(self, *, category: str, owner: str, **kwargs) -> dict:
        return self.operations.register_record(
            store=self.store, project_root=self.project_root, category=category, owner=owner, **kwargs
        )

    def save_uploaded_genome(self, *, category: str, upload: Any) -> Path:
        return _save_uploaded_reference_genome(project_root=self.project_root, category=category, upload=upload)

    def save_uploaded_task_reference(self, *, upload: Any, owner: str) -> Path:
        return _save_uploaded_task_reference(project_root=self.project_root, upload=upload, owner=owner)

    def update_record(self, *, category: str, host_key: str, payload: dict) -> dict:
        return _update_reference_database_record(
            store=self.store, project_root=self.project_root, host_key=host_key, category=category, payload=payload
        )

    def build_index(self, category: str) -> dict:
        return _build_reference_database_index(store=self.store, project_root=self.project_root, category=category)

    def build_single_index(self, *, category: str, host_key: str) -> dict:
        return _build_single_reference_index(
            store=self.store, project_root=self.project_root, category=category, host_key=host_key
        )

    def build_selected_indexes(self, host_keys: list) -> dict:
        built_items = [
            self.build_single_index(category="pathogen", host_key=str(host_key).strip())
            for host_key in host_keys
            if str(host_key or "").strip()
        ]
        return {
            "status": "built",
            "message": f"已完成 {len(built_items)} 条选中病原参考的索引构建",
            "built_items": built_items,
            "items": self.list_records("pathogen"),
        }

    def delete_record(self, *, category: str, host_key: str) -> dict:
        return _delete_reference_database_record(
            store=self.store, project_root=self.project_root, host_key=host_key, category=category
        )

    def queue_cgmlst_panels(self, *, host_keys: list, threshold: str, threads: int, owner: str) -> dict:
        return _queue_pathogen_cgmlst_panels(
            store=self.store,
            project_root=self.project_root,
            host_keys=host_keys,
            threshold=threshold,
            threads=threads,
            owner=owner,
        )
