from __future__ import annotations

import threading
import time
import uuid
from pathlib import Path

from .reference_database_ops import (
    _download_remote_host_genome,
    _fetch_ncbi_taxid_reference_candidates,
    _register_reference_database_record,
)
from .store import PortalStore
from .task_manager import ValidationError


class ReferenceDownloadManager:
    def __init__(self, *, store: PortalStore, project_root: Path) -> None:
        self.store = store
        self.project_root = project_root
        self.jobs: dict[str, dict[str, object]] = {}
        self.lock = threading.Lock()

    def list_jobs(self, category: str) -> list[dict[str, object]]:
        cutoff = time.time() - 3600
        with self.lock:
            for job_id in [key for key, item in self.jobs.items() if item.get("status") in {"success", "failed"} and float(item.get("updated_ts") or 0) < cutoff]:
                self.jobs.pop(job_id, None)
            items = [dict(item) for item in self.jobs.values() if str(item.get("category") or "") == category]
        return sorted(items, key=lambda item: float(item.get("updated_ts") or 0), reverse=True)

    def update(self, job_id: str, **patch: object) -> None:
        with self.lock:
            if job_id in self.jobs:
                self.jobs[job_id].update(patch)
                self.jobs[job_id]["updated_ts"] = time.time()

    def delete(self, job_id: str, category: str) -> None:
        with self.lock:
            current = self.jobs.get(job_id)
            if current is None:
                raise ValidationError("参考下载任务不存在或已清除")
            if str(current.get("category") or "") != category:
                raise ValidationError("参考下载任务分类不匹配")
            if str(current.get("status") or current.get("__jobStatus") or "").lower() in {"downloading", "importing"}:
                raise ValidationError("任务仍在处理中，暂不支持直接删除")
            self.jobs.pop(job_id, None)

    def start(self, *, category: str, payload: dict[str, object], owner: str) -> dict[str, object]:
        provider, query = str(payload.get("provider") or "").strip().lower(), str(payload.get("query") or "").strip()
        ncbi_mode = str(payload.get("ncbi_mode") or "datasets").strip().lower()
        host_name, genome_name, taxid = (str(payload.get(key) or "").strip() for key in ("host_name", "genome_name", "taxid"))
        is_taxid_batch = category == "pathogen" and provider == "ncbi" and ncbi_mode == "taxid_refs"
        effective_taxid = query if is_taxid_batch else taxid
        source_label = str(payload.get("source_label") or ("NCBI TaxID 精选" if is_taxid_batch else provider.upper())).strip()
        job_id = f"refjob-{uuid.uuid4().hex[:12]}"
        now = time.time()
        job = {
            "job_id": job_id, "host_key": f"job:{job_id}", "db_category": category, "category": category,
            "host_name": host_name, "genome_name": genome_name or host_name or (f"TaxID {effective_taxid} 精选参考" if is_taxid_batch else query),
            "taxid": effective_taxid, "source_type": provider, "source_label": source_label,
            "source_accession": query if query.startswith(("GCF_", "GCA_")) else "", "status": "downloading",
            "progress_percent": 1, "message": "等待开始下载", "owner": owner, "created_ts": now, "updated_ts": now,
            "__pending": True, "__jobStatus": "downloading", "__jobPercent": 1,
        }
        with self.lock:
            self.jobs[job_id] = job
        threading.Thread(target=self._run, args=(job_id, category, payload, owner, is_taxid_batch, source_label), name=f"reference-download-{job_id}", daemon=True).start()
        return job

    def _run(self, job_id: str, category: str, payload: dict[str, object], owner: str, is_taxid_batch: bool, source_label: str) -> None:
        provider, query = str(payload.get("provider") or "").strip().lower(), str(payload.get("query") or "").strip()
        description = str(payload.get("description") or "").strip()

        def progress(stage: str, percent: int, message: str) -> None:
            status = "downloading" if stage == "downloading" else "importing"
            self.update(job_id, status=status, progress_percent=percent, message=message, __jobStatus=status, __jobPercent=percent)

        try:
            candidates = _fetch_ncbi_taxid_reference_candidates(query, max_items=int(payload.get("max_genomes") or 20)) if is_taxid_batch else [payload]
            if not candidates:
                raise ValidationError(f"TaxID {query} 没有筛选出可下载的高质量参考基因组")
            records = []
            for index, candidate in enumerate(candidates, start=1):
                accession = str(candidate.get("source_accession") or query).strip()
                display_name = str(candidate.get("genome_name") or payload.get("genome_name") or payload.get("host_name") or accession).strip()
                downloaded = _download_remote_host_genome(project_root=self.project_root, provider=provider, query=accession, host_name=display_name, category=category, ncbi_mode="datasets" if is_taxid_batch else str(payload.get("ncbi_mode") or "datasets"), progress_callback=progress)
                progress("importing", min(99, int(index * 100 / len(candidates))), "下载完成，正在写入参考基因组数据库")
                records.append(_register_reference_database_record(store=self.store, project_root=self.project_root, category=category, host_name=str(payload.get("host_name") or downloaded.get("host_name") or candidate.get("host_name") or "").strip(), genome_name=str(payload.get("genome_name") or downloaded.get("genome_name") or display_name).strip(), taxid=str(payload.get("taxid") or downloaded.get("taxid") or candidate.get("taxid") or query).strip(), source_type=provider, source_label=source_label, source_accession=str(downloaded.get("source_accession") or accession).strip(), source_url=str(downloaded.get("source_url") or "").strip(), source_fasta_path=Path(str(downloaded["fasta_path"])), description=description, owner=owner))
            self.update(job_id, status="success", message=f"imported:{len(records)}", progress_percent=100, __jobStatus="success", __jobPercent=100)
        except Exception as exc:
            self.update(job_id, status="failed", progress_percent=0, message=str(exc), __jobStatus="failed", __jobPercent=100)
