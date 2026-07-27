from __future__ import annotations

import hashlib
import json
import uuid
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from .filesystem_helpers import _resolve_optional_existing_path
from .identity import UserIdentity
from .import_templates import (
    _build_database_import_metadata_items,
    _parse_database_batch_upload,
    _precheck_database_batch_upload,
)
from .reference_operations import ReferenceOperations
from .sample_library_manager import SampleLibraryManager
from .store import PortalStore
from .task_manager import ValidationError


SAMPLE_IMPORT_TEXT_FIELDS = (
    "pathogen_type", "species_name", "mlst_species_name", "mlst_st", "serotype_result", "q20_rate", "q30_rate",
    "completeness", "contamination", "contig_count", "plasmid_count", "resistance_count", "virulence_count",
    "resistance_gene_hits", "virulence_gene_hits", "resistance_mge_hits", "virulence_mge_hits", "genome_id", "taxid",
    "description", "sample_source", "collection_date", "gender", "host_info", "sample_type", "sequencing_method", "note",
)


@dataclass(frozen=True)
class BatchImportService:
    store: PortalStore
    project_root: Path
    sample_manager: SampleLibraryManager
    reference_operations: ReferenceOperations

    def list_runs(self, *, identity: UserIdentity, import_type: str = "", category: str = "", limit: int = 100) -> list[dict]:
        items = self.store.list_batch_import_runs(import_type=import_type, category=category, limit=limit)
        if identity.role != "admin":
            items = [item for item in items if str(item.get("operator") or "") == identity.username]
        for item in items:
            for key in ("precheck_json", "result_json"):
                raw_value = str(item.get(key) or "").strip()
                try:
                    item[key[:-5]] = json.loads(raw_value) if raw_value else {}
                except json.JSONDecodeError:
                    item[key[:-5]] = {}
        return items

    def precheck_reference(self, *, category: str, filename: str, content: bytes, identity: UserIdentity) -> dict:
        result = self.reference_operations.precheck_batch(
            project_root=self.project_root, category=category, filename=filename, content=content
        )
        return self._record_precheck(
            kind="reference", category=category, filename=filename, content=content, result=result, identity=identity
        )

    def import_reference(self, *, category: str, precheck_id: str, identity: UserIdentity) -> dict:
        filename, content, batch_id = self._load_precheck(
            precheck_id=precheck_id, kind="reference", category=category, identity=identity
        )
        result = self.reference_operations.batch_import(
            store=self.store,
            project_root=self.project_root,
            category=category,
            owner=identity.username,
            filename=filename,
            content=content,
        )
        result["batch_id"] = batch_id
        self._complete_run(batch_id, result)
        return result

    def precheck_samples(self, *, filename: str, content: bytes, identity: UserIdentity) -> dict:
        result = _precheck_database_batch_upload(
            project_root=self.project_root,
            sample_library_manager=self.sample_manager,
            filename=filename,
            content=content,
        )
        return self._record_precheck(
            kind="database", category="sample", filename=filename, content=content, result=result, identity=identity
        )

    def import_samples(self, *, precheck_id: str, requested_scope: str, identity: UserIdentity) -> dict:
        filename, content, batch_id = self._load_precheck(
            precheck_id=precheck_id, kind="database", category="sample", identity=identity
        )
        library_scope = (
            requested_scope
            if identity.role == "admin" and requested_scope in {"main", "personal"}
            else ("main" if identity.role == "admin" else "personal")
        )
        imported_items, skipped = [], []
        metadata_templates = self.sample_manager.list_metadata_templates()
        for index, row in enumerate(_parse_database_batch_upload(filename, content), start=2):
            sample_name = str(row.get("sample_name") or "").strip()
            final_fasta_path = str(row.get("final_fasta_path") or "").strip()
            if not sample_name and not final_fasta_path:
                continue
            try:
                province, city, district, detail = (
                    str(row.get(key) or "").strip() for key in ("province", "city", "district", "location_detail")
                )
                country = str(row.get("country") or "").strip() or ("中国" if any([province, city, district]) else "")
                metadata_items = _build_database_import_metadata_items(row, metadata_templates)
                self.sample_manager.validate_metadata_items(metadata_items)
                result = self.sample_manager.import_local_sample(
                    owner=identity.username,
                    owner_group=identity.group_name,
                    library_scope=library_scope,
                    sample_name=sample_name,
                    final_fasta_path=_resolve_optional_existing_path(self.project_root, final_fasta_path),
                    country=country,
                    location_json=json.dumps(
                        {"province": province, "city": city, "district": district, "detail": detail}, ensure_ascii=False
                    ),
                    custom_metadata_json=json.dumps(metadata_items, ensure_ascii=False),
                    **{key: str(row.get(key) or "").strip() for key in SAMPLE_IMPORT_TEXT_FIELDS},
                )
                imported_items.extend(result.get("items") or [])
            except Exception as exc:
                skipped.append({"row": str(index), "sample_name": sample_name or "-", "reason": str(exc)})
        result = {
            "status": "ok",
            "batch_id": batch_id,
            "imported_count": len(imported_items),
            "skipped_count": len(skipped),
            "items": imported_items,
            "skipped": skipped,
        }
        self._complete_run(batch_id, result)
        return result

    def _record_precheck(
        self,
        *,
        kind: str,
        category: str,
        filename: str,
        content: bytes,
        result: dict,
        identity: UserIdentity,
    ) -> dict:
        batch_id = f"batch-{uuid.uuid4().hex}"
        self.store.create_batch_import_run(
            {
                "batch_id": batch_id,
                "import_type": kind,
                "category": category,
                "filename": filename,
                "operator": identity.username,
                "role": identity.role,
                "group_name": identity.group_name,
                "status": "precheck_passed" if result.get("can_import") else "precheck_failed",
                "precheck": result,
                "issue_count": int(result.get("issue_count") or 0),
                "content_hash": hashlib.sha256(content).hexdigest(),
            }
        )
        result["batch_id"] = batch_id
        return self.reference_operations.store_precheck(
            kind=kind,
            category=category,
            owner=identity.username,
            filename=filename,
            content=content,
            result=result,
            batch_id=batch_id,
        )

    def _load_precheck(self, *, precheck_id: str, kind: str, category: str, identity: UserIdentity) -> tuple[str, bytes, str]:
        if not precheck_id:
            raise ValidationError("请先完成批量导入预检，通过后再开始正式导入")
        precheck = self.reference_operations.load_precheck(
            precheck_id=precheck_id, kind=kind, category=category, owner=identity.username
        )
        if precheck is None:
            raise ValidationError("请先完成批量导入预检，通过后再开始正式导入")
        return precheck

    def _complete_run(self, batch_id: str, result: dict) -> None:
        if not batch_id:
            return
        skipped_count = int(result.get("skipped_count") or 0)
        self.store.update_batch_import_run(
            batch_id,
            {
                "status": "completed_with_skips" if skipped_count else "completed",
                "result": result,
                "imported_count": int(result.get("imported_count") or 0),
                "skipped_count": skipped_count,
                "issue_count": skipped_count,
                "completed_at": datetime.utcnow().isoformat(timespec="seconds") + "Z",
            },
        )
