from __future__ import annotations

import gzip
import http.client
import json
import os
import platform
import shlex
import shutil
import subprocess
import threading
import time
import uuid
import zipfile
from datetime import datetime
from pathlib import Path
from urllib.parse import quote, urlencode
from urllib.request import Request, urlopen

from metagenomic_refactor.common import resolve_conda_env_name

from .admin_runtime import _load_conda_root_setting, _resolve_conda_exe_from_root
from .filesystem_helpers import _resolve_optional_existing_path
from .knowledge_interpretation import _normalize_taxonomy_lookup_name
from .import_templates import _parse_database_batch_upload
from .store import PortalStore
from .task_manager import ValidationError


def _register_reference_database_record(
    *,
    store: PortalStore,
    project_root: Path,
    category: str,
    host_name: str,
    genome_name: str,
    taxid: str,
    source_type: str,
    source_label: str,
    source_accession: str,
    source_url: str,
    source_fasta_path: Path,
    description: str,
    owner: str,
) -> dict[str, object]:
    host_name, genome_name = str(host_name or "").strip(), str(genome_name or "").strip()
    if not host_name:
        raise ValidationError("宿主名称不能为空" if category == "host" else "病原名称不能为空")
    if not genome_name:
        raise ValidationError("基因组名称不能为空")
    if not source_fasta_path.is_file():
        raise ValidationError(f"宿主基因组文件不存在: {source_fasta_path}")
    slug = _sanitize_host_slug(genome_name)
    host_root = _reference_database_root(project_root, category)
    fasta_target = host_root / "genomes" / f"{slug}.fa"
    _write_normalized_fasta(source_fasta_path, fasta_target)
    return store.upsert_host_database_record(
        {
            "host_key": slug, "db_category": category, "host_name": host_name, "genome_name": genome_name,
            "taxid": str(taxid or "").strip(), "source_type": source_type, "source_label": source_label,
            "source_accession": source_accession, "source_url": source_url, "description": description,
            "fasta_path": str(fasta_target), "index_prefix": str(host_root / "index" / slug / "reference"),
            "index_status": "pending", "index_message": "已导入，等待构建索引", "owner": owner,
        }
    )


def _batch_import_reference_records(*, store: PortalStore, project_root: Path, category: str, owner: str, filename: str, content: bytes) -> dict[str, object]:
    imported_items, skipped = [], []
    display_name = "宿主" if category == "host" else "病原"
    for index, row in enumerate(_parse_database_batch_upload(filename, content), start=2):
        host_name, genome_name, fasta_path = (str(row.get(key) or "").strip() for key in ("host_name", "genome_name", "fasta_path"))
        if host_name in {"宿主名称", "病原名称", "物种名称"} and genome_name == "基因组名称":
            continue
        if not host_name and not genome_name and not fasta_path:
            continue
        try:
            imported_items.append(
                _register_reference_database_record(
                    store=store, project_root=project_root, category=category, host_name=host_name, genome_name=genome_name,
                    taxid=str(row.get("taxid") or "").strip(), source_type="local",
                    source_label=str(row.get("source_label") or "批量导入").strip(),
                    source_accession=str(row.get("source_accession") or "").strip(), source_url="",
                    source_fasta_path=Path(_resolve_optional_existing_path(project_root, fasta_path)),
                    description=str(row.get("description") or "").strip(), owner=owner,
                )
            )
        except Exception as exc:
            skipped.append({"row": str(index), "host_name": host_name or f"{display_name}名称缺失", "genome_name": genome_name or "-", "reason": str(exc)})
    return {"status": "ok", "imported_count": len(imported_items), "skipped_count": len(skipped), "items": imported_items, "skipped": skipped}


def _looks_like_fasta_path(path: Path) -> bool:
    return path.name.lower().endswith((".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz"))


def _strip_fasta_suffix(name: str) -> str:
    lowered = name.lower()
    for suffix in (".fa.gz", ".fasta.gz", ".fna.gz", ".fasta", ".fna", ".fa", ".gz"):
        if lowered.endswith(suffix):
            return name[: -len(suffix)] or name
    return Path(name).stem


def _collect_reference_import_sources(project_root: Path, raw_value: object) -> list[Path]:
    resolved = _resolve_optional_existing_path(project_root, raw_value)
    if not resolved:
        raise ValidationError("请提供参考基因组路径")
    candidate = Path(resolved)
    if candidate.is_file() and _looks_like_fasta_path(candidate):
        return [candidate]
    if candidate.is_dir():
        paths = sorted((item for item in candidate.rglob("*") if item.is_file() and _looks_like_fasta_path(item)), key=lambda item: str(item).lower())
        if paths:
            return paths
        raise ValidationError(f"目录中未找到 FASTA 文件: {candidate}")
    if candidate.is_file() and candidate.name.lower().endswith(".txt"):
        paths = []
        for raw_line in candidate.read_text(encoding="utf-8", errors="ignore").splitlines():
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            item = Path(line).expanduser()
            item = item.resolve() if item.is_absolute() else (project_root / item).resolve()
            if not item.is_file() or not _looks_like_fasta_path(item):
                raise ValidationError(f"list.txt 中存在无效 FASTA 路径: {item}")
            paths.append(item)
        if paths:
            return paths
    raise ValidationError("仅支持单个 FASTA、包含 FASTA 的目录，或每行一个 FASTA 路径的 list.txt")


def _import_reference_records_from_source(
    *,
    store: PortalStore,
    project_root: Path,
    category: str,
    host_name: str,
    genome_name: str,
    taxid: str,
    source_type: str,
    source_label: str,
    source_accession: str,
    source_url: str,
    source_path_value: object,
    description: str,
    owner: str,
) -> dict[str, object]:
    paths = _collect_reference_import_sources(project_root, source_path_value)
    items = [
        _register_reference_database_record(
            store=store, project_root=project_root, category=category, host_name=host_name,
            genome_name=genome_name if len(paths) == 1 and genome_name else _strip_fasta_suffix(path.name),
            taxid=taxid, source_type=source_type, source_label=source_label, source_accession=source_accession,
            source_url=source_url, source_fasta_path=path, description=description, owner=owner,
        )
        for path in paths
    ]
    return {"status": "ok", "imported_count": len(items), "items": items, "source_mode": "single" if len(items) == 1 else "multi"}

def _reference_database_root(project_root: Path, category: str) -> Path:
    category_name = "host_database" if category == "host" else "pathogen_database"
    return project_root / category_name

def _sanitize_host_slug(value: str) -> str:
    text = str(value or "").strip()
    if not text:
        text = "host_genome"
    normalized = "".join(ch if ch.isalnum() or ch in {"-", "_"} else "_" for ch in text)
    while "__" in normalized:
        normalized = normalized.replace("__", "_")
    return normalized.strip("._") or "host_genome"

def _normalize_ncbi_accession(value: str) -> str:
    text = str(value or "").strip()
    if text.startswith(("GCF_", "GCA_")) and "." not in text:
        parts = text.rsplit("_", 1)
        if len(parts) == 2 and parts[1].isdigit():
            return f"{parts[0]}.{parts[1]}"
    return text

def _resolve_datasets_binary(project_root: Path) -> Path | None:
    scripts_dir = project_root / "scripts"
    system_name = platform.system().lower()
    candidates: list[Path] = []
    if system_name == "linux":
        candidates.append(scripts_dir / "datasets_linux")
    elif system_name == "darwin":
        candidates.append(scripts_dir / "datasets_macos")
    candidates.append(scripts_dir / "datasets")
    which_path = shutil.which("datasets")
    if which_path:
        candidates.append(Path(which_path))
    for candidate in candidates:
        if not candidate.is_file():
            continue
        try:
            mode = candidate.stat().st_mode
            if not os.access(candidate, os.X_OK):
                candidate.chmod(mode | 0o111)
        except OSError:
            pass
        if os.access(candidate, os.X_OK):
            return candidate
    return None

def _save_uploaded_reference_genome(*, project_root: Path, category: str, upload) -> Path:
    filename = Path(str(upload.filename or "").strip()).name
    if not filename:
        raise ValidationError("上传文件名不能为空")
    lowered = filename.lower()
    if not lowered.endswith((".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz")):
        raise ValidationError("仅支持 .fa / .fasta / .fna 及其 gz 压缩文件")
    target_root = _reference_database_root(project_root, category) / "uploads"
    target_root.mkdir(parents=True, exist_ok=True)
    stem = _sanitize_host_slug(Path(filename).stem)
    timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
    suffixes = "".join(Path(filename).suffixes) or Path(filename).suffix or ".fa"
    target = target_root / f"{stem}_{timestamp}{suffixes}"
    upload.save(target)
    try:
        _validate_fasta_content(target)
    except ValidationError:
        target.unlink(missing_ok=True)
        raise
    return target

def _save_uploaded_task_reference(*, project_root: Path, upload, owner: str) -> Path:
    filename = Path(str(upload.filename or "").strip()).name
    if not filename:
        raise ValidationError("上传文件名不能为空")
    lowered = filename.lower()
    if not lowered.endswith((".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz")):
        raise ValidationError("仅支持 .fa / .fasta / .fna 及其 gz 压缩文件")
    owner_slug = _sanitize_host_slug(str(owner or "anonymous")) or "anonymous"
    target_root = project_root / "analysis_tasks" / "_task_reference_uploads" / owner_slug
    target_root.mkdir(parents=True, exist_ok=True)
    stem = _sanitize_host_slug(Path(filename).stem) or "reference"
    timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
    suffixes = "".join(Path(filename).suffixes) or Path(filename).suffix or ".fa"
    target = target_root / f"{stem}_{timestamp}{suffixes}"
    upload.save(target)
    try:
        _validate_fasta_content(target)
    except ValidationError:
        target.unlink(missing_ok=True)
        raise
    return target

def _resolve_reference_fasta_for_update(project_root: Path, category: str, current_record: dict[str, object], payload: dict[str, object]) -> str | None:
    uploaded_fasta_path = str(payload.get("uploaded_fasta_path") or "").strip()
    if uploaded_fasta_path:
        source_path = Path(uploaded_fasta_path)
        host_root = _reference_database_root(project_root, category)
        target = host_root / "genomes" / f"{_sanitize_host_slug(str(payload.get('genome_name') or payload.get('host_name') or current_record.get('genome_name') or current_record.get('host_name') or current_record.get('host_key') or 'reference'))}.fa"
        _write_normalized_fasta(source_path, target)
        return str(target)
    if "fasta_path" not in payload:
        return None
    raw_path = str(payload.get("fasta_path") or "").strip()
    if not raw_path:
        return str(current_record.get("fasta_path") or "")
    source_path = Path(_resolve_optional_existing_path(project_root, raw_path))
    host_root = _reference_database_root(project_root, category)
    target = host_root / "genomes" / f"{_sanitize_host_slug(str(payload.get('genome_name') or payload.get('host_name') or current_record.get('genome_name') or current_record.get('host_name') or current_record.get('host_key') or 'reference'))}.fa"
    _write_normalized_fasta(source_path, target)
    return str(target)

def _write_normalized_fasta(source_path: Path, destination_path: Path) -> None:
    _validate_fasta_content(source_path)
    destination_path.parent.mkdir(parents=True, exist_ok=True)
    if source_path.resolve() == destination_path.resolve():
        return
    if source_path.suffix.lower() == ".gz":
        with gzip.open(source_path, "rt", encoding="utf-8", errors="ignore") as src, destination_path.open("w", encoding="utf-8") as dest:
            shutil.copyfileobj(src, dest)
        return
    try:
        if destination_path.exists():
            destination_path.unlink()
        os.link(source_path, destination_path)
    except OSError:
        shutil.copy2(source_path, destination_path)

def _validate_fasta_content(path: Path) -> None:
    if not path.is_file() or path.stat().st_size <= 0:
        raise ValidationError("FASTA 文件为空")
    opener = gzip.open if path.name.lower().endswith(".gz") else open
    try:
        with opener(path, "rt", encoding="utf-8", errors="ignore") as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if not line:
                    continue
                if not line.startswith(">"):
                    raise ValidationError("FASTA 文件格式错误：首个有效行必须以 > 开头")
                return
    except OSError as exc:
        raise ValidationError(f"FASTA 文件无法读取: {path}") from exc
    raise ValidationError("FASTA 文件为空")

def _extract_ncbi_taxid(raw: object) -> str:
    if raw is None:
        return ""
    if isinstance(raw, (str, int)):
        text = str(raw).strip()
        return text if text.isdigit() else ""
    if isinstance(raw, dict):
        for key in ("tax_id", "taxid", "taxonomy_id"):
            value = _extract_ncbi_taxid(raw.get(key))
            if value:
                return value
        for value in raw.values():
            found = _extract_ncbi_taxid(value)
            if found:
                return found
    if isinstance(raw, list):
        for item in raw:
            found = _extract_ncbi_taxid(item)
            if found:
                return found
    return ""

def _find_nested_text(raw: object, keys: tuple[str, ...]) -> str:
    if isinstance(raw, dict):
        for key in keys:
            value = raw.get(key)
            text = str(value or "").strip()
            if text:
                return text
        for value in raw.values():
            found = _find_nested_text(value, keys)
            if found:
                return found
    if isinstance(raw, list):
        for item in raw:
            found = _find_nested_text(item, keys)
            if found:
                return found
    return ""

def _extract_ncbi_organism_name(report: dict[str, object]) -> str:
    candidates = [
        report.get("organism", {}).get("organism_name") if isinstance(report.get("organism"), dict) else None,
        report.get("organism", {}).get("scientific_name") if isinstance(report.get("organism"), dict) else None,
        report.get("organism_name"),
        report.get("scientific_name"),
        report.get("current_organism_name"),
    ]
    for value in candidates:
        text = str(value or "").strip()
        if text:
            return text
    return _find_nested_text(report, ("organism_name", "scientific_name", "current_organism_name"))

def _extract_ncbi_assembly_name(report: dict[str, object]) -> str:
    strain_name = _find_nested_text(
        report,
        (
            "strain",
            "isolate",
            "cultivar",
            "breed",
            "ecotype",
            "serovar",
            "bioproject_lineage",
        ),
    )
    if strain_name:
        return strain_name
    assembly_info = report.get("assembly_info")
    if isinstance(assembly_info, dict):
        for key in ("assembly_name", "assembly_level", "assembly_type"):
            text = str(assembly_info.get(key) or "").strip()
            if text:
                return text
    for key in ("assembly_name", "display_name"):
        text = str(report.get(key) or "").strip()
        if text:
            return text
    return _find_nested_text(report, ("assembly_name", "display_name", "submitter_name"))

def _normalize_ncbi_species_name(organism_name: str, assembly_name: str) -> str:
    organism_text = str(organism_name or "").strip()
    assembly_text = str(assembly_name or "").strip()
    if not organism_text:
        return ""
    if assembly_text and organism_text.lower().endswith(assembly_text.lower()):
        trimmed = organism_text[: len(organism_text) - len(assembly_text)].strip(" ,;:-_/")
        if trimmed:
            return trimmed
    for marker in (" strain ", " substr. ", " subsp. ", " isolate "):
        index = organism_text.lower().find(marker)
        if index > 0:
            trimmed = organism_text[:index].strip()
            if trimmed:
                return trimmed
    return organism_text

def _extract_ncbi_biosource_value(raw: object, allowed_types: tuple[str, ...] = ("strain", "isolate")) -> str:
    biosource = raw.get("biosource") if isinstance(raw, dict) else {}
    infraspecies = biosource.get("infraspecieslist") if isinstance(biosource, dict) else None
    if isinstance(infraspecies, list):
        for item in infraspecies:
            if not isinstance(item, dict):
                continue
            sub_type = str(item.get("sub_type") or item.get("subtype") or "").strip().lower()
            sub_value = str(item.get("sub_value") or item.get("subname") or item.get("value") or "").strip()
            if sub_type in allowed_types and sub_value:
                return sub_value
    return ""

def _score_ncbi_assembly_quality(doc: dict[str, object]) -> tuple[int, int, int, int, int, str, str]:
    refseq_category = str(doc.get("refseq_category") or "").strip().lower()
    assembly_level = str(doc.get("assemblylevel") or doc.get("assembly_level") or "").strip().lower()
    assembly_status = str(doc.get("assemblystatus") or doc.get("versionstatus") or "").strip().lower()
    accession = _normalize_ncbi_accession(str(doc.get("assemblyaccession") or doc.get("assembly_accession") or "").strip())
    release_date = str(doc.get("seqreleasedate") or doc.get("submissiondate") or doc.get("lastupdatedate") or "").strip()
    refseq_score = 3 if refseq_category == "reference genome" else 2 if refseq_category == "representative genome" else 0
    level_score = {
        "complete genome": 4,
        "chromosome": 3,
        "scaffold": 2,
        "contig": 1,
    }.get(assembly_level, 0)
    accession_score = 1 if accession.startswith("GCF_") else 0
    latest_score = 1 if "latest" in assembly_status or "current" in assembly_status else 0
    preferred_score = 1 if refseq_score > 0 or level_score >= 3 else 0
    return (preferred_score, refseq_score, level_score, accession_score, latest_score, release_date, accession)

def _build_ncbi_reference_dedup_key(doc: dict[str, object]) -> str:
    species_name = _normalize_taxonomy_lookup_name(doc.get("host_name") or "")
    strain_name = _normalize_taxonomy_lookup_name(doc.get("strain_name") or "")
    genome_name = _normalize_taxonomy_lookup_name(doc.get("genome_name") or "")
    if species_name and strain_name:
        return f"{species_name}::{strain_name}"
    if species_name and genome_name:
        return f"{species_name}::{genome_name}"
    return genome_name or species_name

def _fetch_ncbi_taxid_reference_candidates(taxid: str, *, max_items: int = 20) -> list[dict[str, str]]:
    taxid = str(taxid or "").strip()
    if not taxid.isdigit():
        raise ValidationError("TaxID 必须是纯数字")
    limit = max(1, min(200, int(max_items or 20)))
    search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urlencode(
        {
            "db": "assembly",
            "term": f"txid{taxid}[Organism:exp]",
            "retmode": "json",
            "retmax": str(min(max(limit * 12, 100), 500)),
        }
    )
    try:
        with urlopen(Request(search_url, headers={"User-Agent": "bac-analysis-portal/1.0"}), timeout=60) as response:
            search_payload = json.loads(response.read().decode("utf-8", errors="ignore"))
    except Exception as exc:
        raise ValidationError(f"查询 TaxID {taxid} 的 NCBI 装配列表失败: {exc}") from exc
    id_list = (((search_payload or {}).get("esearchresult") or {}).get("idlist") or [])
    if not id_list:
        raise ValidationError(f"NCBI 中没有找到 TaxID {taxid} 对应的装配记录")

    docs: list[dict[str, object]] = []
    for start in range(0, len(id_list), 100):
        batch = id_list[start : start + 100]
        summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?" + urlencode(
            {
                "db": "assembly",
                "id": ",".join(str(item) for item in batch),
                "retmode": "json",
            }
        )
        try:
            with urlopen(Request(summary_url, headers={"User-Agent": "bac-analysis-portal/1.0"}), timeout=60) as response:
                summary_payload = json.loads(response.read().decode("utf-8", errors="ignore"))
        except Exception:
            continue
        result = (summary_payload or {}).get("result") or {}
        for item_id in batch:
            doc = result.get(str(item_id)) or {}
            if isinstance(doc, dict) and doc:
                docs.append(doc)

    if not docs:
        raise ValidationError(f"TaxID {taxid} 的装配摘要读取失败")

    candidates: list[dict[str, object]] = []
    for doc in docs:
        accession = _normalize_ncbi_accession(str(doc.get("assemblyaccession") or doc.get("assembly_accession") or "").strip())
        if not accession.startswith(("GCF_", "GCA_")):
            continue
        organism_name = _extract_ncbi_organism_name(doc)
        assembly_name = _extract_ncbi_assembly_name(doc)
        species_name = _normalize_ncbi_species_name(organism_name, assembly_name) or organism_name
        strain_name = _extract_ncbi_biosource_value(doc)
        refseq_category = str(doc.get("refseq_category") or "").strip().lower()
        assembly_level = str(doc.get("assemblylevel") or doc.get("assembly_level") or "").strip().lower()
        score = _score_ncbi_assembly_quality(doc)
        candidates.append(
            {
                "source_accession": accession,
                "host_name": species_name.strip(),
                "genome_name": (strain_name or assembly_name or accession).strip(),
                "taxid": str(doc.get("taxid") or taxid).strip(),
                "strain_name": strain_name.strip(),
                "assembly_level": assembly_level,
                "refseq_category": refseq_category,
                "quality_score": score,
            }
        )

    if not candidates:
        raise ValidationError(f"TaxID {taxid} 没有可下载的参考基因组 accession")

    candidates.sort(key=lambda item: item.get("quality_score") or (), reverse=True)

    selected: list[dict[str, str]] = []
    seen_accessions: set[str] = set()
    seen_groups: set[str] = set()
    quality_tiers = (
        lambda item: (item.get("quality_score") or (0,))[0] >= 1 and str(item.get("assembly_level") or "") in {"complete genome", "chromosome"},
        lambda item: (item.get("quality_score") or (0,))[0] >= 1,
        lambda item: str(item.get("assembly_level") or "") in {"complete genome", "chromosome", "scaffold"},
        lambda item: True,
    )
    for predicate in quality_tiers:
        for item in candidates:
            accession = str(item.get("source_accession") or "").strip()
            if not accession or accession in seen_accessions:
                continue
            dedup_key = _build_ncbi_reference_dedup_key(item)
            if dedup_key and dedup_key in seen_groups:
                continue
            if not predicate(item):
                continue
            selected.append(
                {
                    "source_accession": accession,
                    "host_name": str(item.get("host_name") or "").strip(),
                    "genome_name": str(item.get("genome_name") or accession).strip(),
                    "taxid": str(item.get("taxid") or taxid).strip(),
                }
            )
            seen_accessions.add(accession)
            if dedup_key:
                seen_groups.add(dedup_key)
            if len(selected) >= limit:
                return selected
    return selected[:limit]

def _fetch_ncbi_accession_metadata(accession: str) -> dict[str, str]:
    accession = str(accession or "").strip()
    if not accession:
        return {}
    search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urlencode({
        "db": "assembly",
        "term": f"{accession}[Assembly Accession]",
        "retmode": "json",
    })
    try:
        with urlopen(Request(search_url, headers={"User-Agent": "bac-analysis-portal/1.0"}), timeout=60) as response:
            search_payload = json.loads(response.read().decode("utf-8", errors="ignore"))
    except Exception:
        return {}
    id_list = (((search_payload or {}).get("esearchresult") or {}).get("idlist") or [])
    if not id_list:
        return {}
    summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?" + urlencode({
        "db": "assembly",
        "id": str(id_list[0]),
        "retmode": "json",
    })
    try:
        with urlopen(Request(summary_url, headers={"User-Agent": "bac-analysis-portal/1.0"}), timeout=60) as response:
            summary_payload = json.loads(response.read().decode("utf-8", errors="ignore"))
    except Exception:
        return {}
    result = (summary_payload or {}).get("result") or {}
    doc = result.get(str(id_list[0])) or {}
    species_name = str(doc.get("speciesname") or doc.get("organism") or "").strip()
    assembly_name = str(doc.get("assemblyname") or "").strip()
    taxid = str(doc.get("taxid") or "").strip()
    return {
        "host_name": species_name,
        "genome_name": _extract_ncbi_biosource_value(doc) or assembly_name,
        "taxid": taxid,
        "source_accession": accession,
    }

def _read_ncbi_download_metadata(download_target: Path) -> dict[str, str]:
    if not zipfile.is_zipfile(download_target):
        return {}
    try:
        with zipfile.ZipFile(download_target) as archive:
            report_names = [name for name in archive.namelist() if name.endswith("assembly_data_report.jsonl")]
            for report_name in report_names:
                with archive.open(report_name) as handle:
                    for raw_line in handle:
                        line = raw_line.decode("utf-8", errors="ignore").strip()
                        if not line:
                            continue
                        try:
                            report = json.loads(line)
                        except json.JSONDecodeError:
                            continue
                        organism_name = _extract_ncbi_organism_name(report)
                        taxid = _extract_ncbi_taxid(report)
                        assembly_name = _extract_ncbi_assembly_name(report)
                        species_name = _normalize_ncbi_species_name(organism_name, assembly_name)
                        if organism_name or taxid or assembly_name:
                            return {
                                "host_name": species_name or organism_name,
                                "taxid": taxid,
                                "genome_name": assembly_name,
                            }
    except OSError:
        return {}
    return {}

def _validate_downloaded_reference_file(download_target: Path) -> None:
    if zipfile.is_zipfile(download_target):
        return
    head = download_target.read_bytes()[:4096]
    text = head.decode("utf-8", errors="ignore").strip()
    if not text:
        raise ValidationError("下载结果为空")
    if text.startswith(">"):
        return
    if text.startswith("{") or text.startswith("["):
        try:
            payload = json.loads(text)
            error_text = json.dumps(payload, ensure_ascii=False)[:400]
        except json.JSONDecodeError:
            error_text = text[:200]
        raise ValidationError(f"NCBI 未返回基因组文件，返回内容为: {error_text}")
    raise ValidationError(f"下载结果不是 FASTA 文件，返回内容为: {text[:200]}")

def _download_with_resume(*, source_url: str, destination: Path, progress_callback=None, max_attempts: int = 5) -> None:
    downloaded_bytes = destination.stat().st_size if destination.exists() else 0
    total_bytes = 0
    chunk_size = 1024 * 1024
    last_error: Exception | None = None

    for attempt in range(1, max_attempts + 1):
        request_obj = Request(source_url, headers={"User-Agent": "bac-analysis-portal/1.0"})
        if downloaded_bytes > 0:
            request_obj.add_header("Range", f"bytes={downloaded_bytes}-")
        try:
            with urlopen(request_obj, timeout=180) as response:
                content_range = str(response.headers.get("Content-Range") or "").strip()
                if downloaded_bytes > 0 and not content_range:
                    downloaded_bytes = 0
                    if destination.exists():
                        destination.unlink(missing_ok=True)
                with destination.open("ab" if downloaded_bytes > 0 else "wb") as handle:
                    if content_range and "/" in content_range:
                        try:
                            total_bytes = int(content_range.rsplit("/", 1)[1])
                        except (TypeError, ValueError):
                            total_bytes = total_bytes
                    elif response.headers.get("Content-Length"):
                        try:
                            response_length = int(response.headers.get("Content-Length") or 0)
                        except (TypeError, ValueError):
                            response_length = 0
                        total_bytes = downloaded_bytes + response_length if downloaded_bytes > 0 else response_length

                    while True:
                        try:
                            chunk = response.read(chunk_size)
                        except http.client.IncompleteRead as exc:
                            chunk = exc.partial or b""
                            if chunk:
                                handle.write(chunk)
                                downloaded_bytes += len(chunk)
                                if progress_callback is not None:
                                    percent = 10
                                    if total_bytes > 0:
                                        percent = max(2, min(94, int(downloaded_bytes * 94 / total_bytes)))
                                    progress_callback("downloading", percent, f"下载连接中断，正在重试（第 {attempt} 次）")
                            raise
                        if not chunk:
                            break
                        handle.write(chunk)
                        downloaded_bytes += len(chunk)
                        if progress_callback is not None:
                            percent = 10
                            if total_bytes > 0:
                                percent = max(2, min(94, int(downloaded_bytes * 94 / total_bytes)))
                            progress_callback("downloading", percent, "正在下载参考基因组")
            return
        except http.client.IncompleteRead as exc:
            last_error = exc
            time.sleep(min(2 * attempt, 8))
            continue
        except Exception as exc:
            last_error = exc
            if attempt >= max_attempts:
                break
            time.sleep(min(2 * attempt, 8))
            continue

    raise ValidationError(f"下载宿主基因组失败: {last_error}")

def _download_with_datasets(
    *,
    datasets_binary: Path,
    accession: str,
    destination: Path,
    progress_callback=None,
    max_attempts: int = 3,
) -> None:
    last_error = ""
    for attempt in range(1, max_attempts + 1):
        if destination.exists():
            destination.unlink(missing_ok=True)
        command = [
            str(datasets_binary),
            "download",
            "genome",
            "accession",
            accession,
            "--filename",
            str(destination),
            "--include",
            "genome",
        ]
        process = subprocess.Popen(
            command,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            text=True,
        )
        last_size = -1
        stable_rounds = 0
        while True:
            return_code = process.poll()
            current_size = destination.stat().st_size if destination.exists() else 0
            if current_size != last_size:
                stable_rounds = 0
                last_size = current_size
                if progress_callback is not None:
                    progress_callback("downloading", 10, f"NCBI datasets 下载中（第 {attempt}/{max_attempts} 次），已接收 {current_size} bytes")
            else:
                stable_rounds += 1
                if progress_callback is not None and stable_rounds % 5 == 0:
                    progress_callback("downloading", 10, f"NCBI datasets 下载中（第 {attempt}/{max_attempts} 次），请稍候")
            if return_code is not None:
                if return_code != 0:
                    last_error = f"datasets 退出码 {return_code}"
                    break
                if not destination.exists():
                    last_error = f"datasets 未生成输出文件 {destination}"
                    break
                return
            time.sleep(1)
        if attempt < max_attempts:
            time.sleep(min(2 * attempt, 6))
    raise ValidationError(f"下载宿主基因组失败: {last_error or 'datasets 多次重试后仍失败'}")

def _download_remote_host_genome(*, project_root: Path, provider: str, query: str, host_name: str, category: str = "host", ncbi_mode: str = "datasets", progress_callback=None) -> dict[str, str]:
    provider = str(provider or "").strip().lower()
    ncbi_mode = str(ncbi_mode or "datasets").strip().lower()
    query = str(query or "").strip()
    if not query:
        raise ValidationError("请填写 accession 或下载链接")
    host_root = _reference_database_root(project_root, category)
    downloads_dir = host_root / "downloads"
    downloads_dir.mkdir(parents=True, exist_ok=True)
    source_url = query if query.startswith(("http://", "https://")) else ""
    source_accession = ""
    if provider == "ncbi" and not source_url:
        source_accession = _normalize_ncbi_accession(query)
    if provider in {"gtdb", "ensembl"} and not source_url:
        label = "Ensembl" if provider == "ensembl" else "GTDB"
        raise ValidationError(f"{label} 下载请提供可直接访问的 FASTA 或压缩包链接")

    slug_seed = host_name or source_accession or (Path(source_url.split("?")[0]).stem if source_url else "") or provider
    slug = _sanitize_host_slug(slug_seed)
    if provider == "ncbi" and source_accession and ncbi_mode == "datasets":
        datasets_binary = _resolve_datasets_binary(project_root)
        if datasets_binary is None:
            raise ValidationError("未找到可用的 datasets 软件，请确认 scripts/datasets_linux 或 scripts/datasets_macos 存在且可执行")
        download_target = downloads_dir / f"{slug}.ncbi_dataset.zip"
        if progress_callback is not None:
            progress_callback("downloading", 5, "正在通过 NCBI datasets 下载参考基因组")
        _download_with_datasets(
            datasets_binary=datasets_binary,
            accession=source_accession,
            destination=download_target,
            progress_callback=progress_callback,
        )
    else:
        if provider == "ncbi" and source_accession:
            source_url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{quote(source_accession)}/download?include_annotation_type=GENOME_FASTA"
        download_target = downloads_dir / f"{slug}.download"
        _download_with_resume(source_url=source_url, destination=download_target, progress_callback=progress_callback)
        _validate_downloaded_reference_file(download_target)

    ncbi_metadata = _read_ncbi_download_metadata(download_target) if provider == "ncbi" else {}

    if zipfile.is_zipfile(download_target):
        if progress_callback is not None:
            progress_callback("importing", 96, "下载完成，正在提取 FASTA")
        with zipfile.ZipFile(download_target) as archive:
            fasta_members = [name for name in archive.namelist() if name.lower().endswith((".fa", ".fasta", ".fna"))]
            if not fasta_members:
                raise ValidationError("下载文件中未找到 FASTA")
            chosen = fasta_members[0]
            extracted = downloads_dir / f"{slug}{Path(chosen).suffix or '.fa'}"
            with archive.open(chosen) as src, extracted.open("wb") as dest:
                shutil.copyfileobj(src, dest)
            download_target.unlink(missing_ok=True)
            download_target = extracted
    else:
        if progress_callback is not None:
            progress_callback("importing", 96, "下载完成，正在整理参考文件")
        suffix = Path(source_url.split("?")[0]).suffix or ".fa"
        normalized = downloads_dir / f"{slug}{suffix}"
        download_target.replace(normalized)
        download_target = normalized

    inferred_name = host_name or ncbi_metadata.get("host_name") or source_accession or Path(download_target).stem
    if provider == "ncbi" and (not ncbi_metadata.get("host_name") or not ncbi_metadata.get("genome_name")) and source_accession:
        accession_metadata = _fetch_ncbi_accession_metadata(source_accession)
        for key, value in accession_metadata.items():
            if value and not ncbi_metadata.get(key):
                ncbi_metadata[key] = value
        inferred_name = host_name or ncbi_metadata.get("host_name") or source_accession or Path(download_target).stem
    return {
        "host_name": inferred_name,
        "genome_name": str(ncbi_metadata.get("genome_name") or "").strip(),
        "taxid": str(ncbi_metadata.get("taxid") or "").strip(),
        "source_accession": source_accession,
        "source_url": source_url,
        "fasta_path": str(download_target),
    }

def _append_fasta_to_handle(handle, fasta_path: Path) -> None:
    with fasta_path.open("r", encoding="utf-8", errors="ignore") as source:
        content = source.read().strip()
        if not content:
            return
        handle.write(content)
        handle.write("\n")

def _run_host_index_command(command: list[str], workdir: Path) -> tuple[bool, str]:
    tool_name = command[0]
    if shutil.which(tool_name) is None:
        return False, f"{tool_name} 未安装，已跳过"
    completed = subprocess.run(command, cwd=workdir, capture_output=True, text=True)
    if completed.returncode != 0:
        error_text = (completed.stderr or completed.stdout or "").strip() or f"{tool_name} 执行失败"
        return False, error_text
    return True, f"{tool_name} 已完成"

def _reference_panel_root(project_root: Path, category: str, panel_type: str = "cgmlst") -> Path:
    return _reference_database_root(project_root, category) / "panels" / panel_type

def _resolve_chewbbaca_command(store: PortalStore | None = None) -> list[str]:
    env_name = "genomad_aux"
    if store is not None:
        env_name = resolve_conda_env_name(str(store.get_setting("conda_env_vfind", env_name) or env_name).strip() or env_name)
    raw = str(os.environ.get("CHEWBBACA_BIN") or "").strip()
    if raw:
        return shlex.split(raw)
    conda_root = _load_conda_root_setting(store) if store is not None else ""
    return [_resolve_conda_exe_from_root(conda_root), "run", "-n", env_name, "--no-capture-output", "chewBBACA.py"]

def _resolve_chewbbaca_training_file(species_slug: str) -> Path | None:
    specific_env = str(os.environ.get(f"CHEWBBACA_PTF_{species_slug.upper()}") or "").strip()
    generic_env = str(os.environ.get("CHEWBBACA_PTF") or "").strip()
    for candidate in (specific_env, generic_env):
        if not candidate:
            continue
        path = Path(candidate).expanduser().resolve()
        if path.is_file():
            return path
    return None

def _run_reference_panel_command(command: list[str], workdir: Path) -> None:
    executable = command[0]
    if shutil.which(executable) is None and not Path(executable).is_file():
        raise ValidationError(f"{executable} 未安装或不可执行")
    completed = subprocess.run(command, cwd=workdir, capture_output=True, text=True)
    if completed.returncode != 0:
        error_text = (completed.stderr or completed.stdout or "").strip() or "命令执行失败"
        raise ValidationError(error_text)

def _link_or_copy_reference_fasta(source: Path, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() or target.is_symlink():
        target.unlink()
    try:
        target.symlink_to(source)
    except OSError:
        shutil.copy2(source, target)

def _materialize_cgmlst_loci(schema_dir: Path, loci_list_path: Path, target_dir: Path) -> int:
    if not loci_list_path.is_file():
        return 0
    target_dir.mkdir(parents=True, exist_ok=True)
    copied = 0
    names = [line.strip() for line in loci_list_path.read_text(encoding="utf-8", errors="ignore").splitlines() if line.strip()]
    for locus in names:
        candidates = [
            schema_dir / locus,
            schema_dir / f"{locus}.fasta",
            schema_dir / f"{locus}.fa",
            schema_dir / f"{locus}.fna",
        ]
        for candidate in candidates:
            if not candidate.is_file():
                continue
            shutil.copy2(candidate, target_dir / candidate.name)
            copied += 1
            break
    return copied

def _run_single_cgmlst_panel_build(
    *,
    store: PortalStore,
    project_root: Path,
    panel_key: str,
    species_name: str,
    selected_records: list[dict[str, object]],
    threshold: str,
    threads: int,
    owner: str,
) -> None:
    species_slug = _sanitize_host_slug(species_name or panel_key)
    panel_root = _reference_panel_root(project_root, "pathogen", "cgmlst") / species_slug / panel_key
    input_dir = panel_root / "input_genomes"
    schema_dir = panel_root / "schema_seed"
    allele_call_dir = panel_root / "allele_call"
    cgmlst_dir = panel_root / "cgmlst_panel"
    loci_dir = cgmlst_dir / "loci"
    panel_root.mkdir(parents=True, exist_ok=True)
    input_dir.mkdir(parents=True, exist_ok=True)
    selected_host_keys = [str(row.get("host_key") or "").strip() for row in selected_records if str(row.get("host_key") or "").strip()]
    store.upsert_reference_panel(
        {
            "panel_key": panel_key,
            "db_category": "pathogen",
            "panel_type": "cgmlst",
            "species_name": species_name,
            "species_slug": species_slug,
            "selected_host_keys_json": selected_host_keys,
            "genome_count": len(selected_records),
            "schema_dir": str(schema_dir),
            "cgmlst_dir": str(cgmlst_dir),
            "threshold": threshold,
            "status": "running",
            "message": "正在准备参考基因组并启动 chewBBACA",
            "owner": owner,
        }
    )
    try:
        if len(selected_records) < 2:
            raise ValidationError("构建 cgMLST panel 至少需要 2 条同物种参考基因组")
        for row in selected_records:
            fasta_path = Path(str(row.get("fasta_path") or "")).expanduser().resolve()
            if not fasta_path.is_file():
                raise ValidationError(f"参考 FASTA 不存在: {fasta_path}")
            target = input_dir / f"{_sanitize_host_slug(str(row.get('genome_name') or row.get('host_key') or fasta_path.stem))}.fa"
            _link_or_copy_reference_fasta(fasta_path, target)

        chewbbaca = _resolve_chewbbaca_command(store)
        ptf_path = _resolve_chewbbaca_training_file(species_slug)
        create_schema_cmd = chewbbaca + [
            "CreateSchema",
            "-i",
            str(input_dir),
            "-o",
            str(schema_dir),
            "--cpu",
            str(threads),
        ]
        if ptf_path is not None:
            create_schema_cmd.extend(["--ptf", str(ptf_path)])
        _run_reference_panel_command(create_schema_cmd, panel_root)

        allele_call_cmd = chewbbaca + [
            "AlleleCall",
            "-i",
            str(input_dir),
            "-g",
            str(schema_dir),
            "-o",
            str(allele_call_dir),
            "--cpu",
            str(threads),
        ]
        _run_reference_panel_command(allele_call_cmd, panel_root)

        extract_cmd = chewbbaca + [
            "ExtractCgMLST",
            "-i",
            str(allele_call_dir / "results_alleles.tsv"),
            "-o",
            str(cgmlst_dir),
            "--t",
            str(threshold),
        ]
        _run_reference_panel_command(extract_cmd, panel_root)

        loci_list_candidates = sorted(cgmlst_dir.glob("cgMLSTschema*.txt"))
        loci_list_path = loci_list_candidates[0] if loci_list_candidates else Path("")
        copied_count = _materialize_cgmlst_loci(schema_dir, loci_list_path, loci_dir) if loci_list_candidates else 0
        results_alleles_path = allele_call_dir / "results_alleles.tsv"
        message = f"cgMLST panel 构建完成，纳入 {len(selected_records)} 条参考；提取 {copied_count} 个 loci"
        store.upsert_reference_panel(
            {
                "panel_key": panel_key,
                "db_category": "pathogen",
                "panel_type": "cgmlst",
                "species_name": species_name,
                "species_slug": species_slug,
                "selected_host_keys_json": selected_host_keys,
                "genome_count": len(selected_records),
                "schema_dir": str(schema_dir),
                "cgmlst_dir": str(cgmlst_dir),
                "loci_list_path": str(loci_list_path) if loci_list_candidates else "",
                "results_alleles_path": str(results_alleles_path) if results_alleles_path.is_file() else "",
                "threshold": threshold,
                "status": "built",
                "message": message,
                "owner": owner,
            }
        )
    except Exception as exc:
        store.upsert_reference_panel(
            {
                "panel_key": panel_key,
                "db_category": "pathogen",
                "panel_type": "cgmlst",
                "species_name": species_name,
                "species_slug": species_slug,
                "selected_host_keys_json": selected_host_keys,
                "genome_count": len(selected_records),
                "schema_dir": str(schema_dir),
                "cgmlst_dir": str(cgmlst_dir),
                "threshold": threshold,
                "status": "failed",
                "message": str(exc),
                "owner": owner,
            }
        )

def _queue_pathogen_cgmlst_panels(
    *,
    store: PortalStore,
    project_root: Path,
    host_keys: list[object],
    threshold: str,
    threads: int,
    owner: str,
) -> dict[str, object]:
    selected_keys = [str(item or "").strip() for item in host_keys if str(item or "").strip()]
    if not selected_keys:
        raise ValidationError("请先勾选病原参考基因组")
    all_records = {
        str(item.get("host_key") or ""): item
        for item in store.list_host_database("pathogen")
    }
    selected_records = [all_records[key] for key in selected_keys if key in all_records]
    if not selected_records:
        raise ValidationError("未找到对应的病原参考记录")

    grouped: dict[str, list[dict[str, object]]] = {}
    for row in selected_records:
        species_name = str(row.get("host_name") or "").strip()
        if not species_name:
            raise ValidationError(f"参考记录缺少病原名称，无法构建 cgMLST panel: {row.get('host_key')}")
        grouped.setdefault(species_name, []).append(row)

    queued_items: list[dict[str, object]] = []
    timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
    for species_name, rows in grouped.items():
        species_slug = _sanitize_host_slug(species_name)
        panel_key = f"{species_slug}_{timestamp}_{uuid.uuid4().hex[:6]}"
        record = store.upsert_reference_panel(
            {
                "panel_key": panel_key,
                "db_category": "pathogen",
                "panel_type": "cgmlst",
                "species_name": species_name,
                "species_slug": species_slug,
                "selected_host_keys_json": [str(row.get("host_key") or "") for row in rows],
                "genome_count": len(rows),
                "threshold": threshold,
                "status": "pending",
                "message": "已加入构建队列",
                "owner": owner,
            }
        )
        queued_items.append(record)
        threading.Thread(
            target=_run_single_cgmlst_panel_build,
            kwargs={
                "store": store,
                "project_root": project_root,
                "panel_key": panel_key,
                "species_name": species_name,
                "selected_records": rows,
                "threshold": threshold,
                "threads": threads,
                "owner": owner,
            },
            daemon=True,
            name=f"cgmlst-panel-{panel_key}",
        ).start()

    return {
        "status": "queued",
        "message": f"已按物种提交 {len(queued_items)} 个 cgMLST panel 构建任务",
        "items": queued_items,
    }

def _build_single_reference_index(*, store: PortalStore, project_root: Path, category: str, host_key: str) -> dict[str, object]:
    record = store.get_host_database_record(host_key)
    fasta_path = Path(str(record.get("fasta_path") or ""))
    if str(record.get("db_category") or "") != category:
        raise ValidationError("参考记录分类不匹配")
    if not fasta_path.is_file():
        raise ValidationError("该参考记录缺少可用的 FASTA 文件")

    index_dir = _reference_database_root(project_root, category) / "index" / _sanitize_host_slug(str(record.get("host_key") or host_key))
    index_dir.mkdir(parents=True, exist_ok=True)
    stem = _sanitize_host_slug(str(record.get("host_key") or host_key))
    prefix = index_dir / "reference"
    results: list[str] = []
    for command in (
        ["samtools", "faidx", str(fasta_path)],
        ["minimap2", "-d", str(prefix.with_suffix(".mmi")), str(fasta_path)],
        ["bowtie2-build", str(fasta_path), str(prefix)],
    ):
        ok, message = _run_host_index_command(command, index_dir)
        results.append(message)

    built_files = [str(path) for path in [fasta_path.with_suffix(fasta_path.suffix + ".fai"), prefix.with_suffix(".mmi")] if path.exists()]
    built_files.extend(str(path) for path in index_dir.glob("reference*.bt2*"))
    status = "built" if any("已完成" in message for message in results) else "prepared"
    message = "；".join(results)
    updated = store.update_host_database_record(
        host_key,
        index_prefix=str(prefix),
        index_status=status,
        index_message=message,
    )
    return {
        "status": status,
        "host_key": host_key,
        "index_prefix": str(prefix),
        "built_files": built_files,
        "message": message,
        "record": updated,
    }

def _build_reference_database_index(*, store: PortalStore, project_root: Path, category: str) -> dict[str, object]:
    records = store.list_host_database(category)
    available_records = [row for row in records if Path(str(row.get("fasta_path") or "")).is_file()]
    if not available_records:
        raise ValidationError("当前参考数据库中暂无可用于构建索引的 FASTA")
    built_items = [
        _build_single_reference_index(store=store, project_root=project_root, category=category, host_key=str(row.get("host_key") or ""))
        for row in available_records
    ]
    return {
        "status": "built",
        "message": f"已完成 {len(built_items)} 条参考记录的索引构建",
        "built_items": built_items,
        "items": store.list_host_database(category),
    }

def _update_reference_database_record(*, store: PortalStore, project_root: Path, host_key: str, category: str, payload: dict[str, object]) -> dict[str, object]:
    current = store.get_host_database_record(host_key)
    if str(current.get("db_category") or "") != category:
        raise ValidationError("参考记录分类不匹配")
    next_host_name = str(payload.get("host_name") or current.get("host_name") or "").strip()
    next_genome_name = str(payload.get("genome_name") or current.get("genome_name") or current.get("host_name") or "").strip()
    next_fasta_path = _resolve_reference_fasta_for_update(project_root, category, current, payload)
    fasta_changed = next_fasta_path is not None and next_fasta_path != str(current.get("fasta_path") or "")
    return store.update_host_database_record(
        host_key,
        host_name=next_host_name,
        genome_name=next_genome_name,
        taxid=str(payload.get("taxid") or current.get("taxid") or "").strip(),
        source_label=str(payload.get("source_label") or current.get("source_label") or "").strip(),
        source_accession=str(payload.get("source_accession") or current.get("source_accession") or "").strip(),
        source_url=str(payload.get("source_url") or current.get("source_url") or "").strip(),
        description=str(payload.get("description") or current.get("description") or "").strip(),
        fasta_path=next_fasta_path if next_fasta_path is not None else None,
        index_prefix="" if fasta_changed else None,
        index_status="pending" if fasta_changed else None,
        index_message="FASTA 已更新，请重新构建索引" if fasta_changed else None,
    )

def _delete_reference_database_record(*, store: PortalStore, project_root: Path, host_key: str, category: str) -> dict[str, object]:
    current = store.get_host_database_record(host_key)
    if str(current.get("db_category") or "") != category:
        raise ValidationError("参考记录分类不匹配")
    host_root = _reference_database_root(project_root, category)
    downloads_dir = host_root / "downloads"
    fasta_path = Path(str(current.get("fasta_path") or "")).expanduser()
    index_prefix = Path(str(current.get("index_prefix") or "")).expanduser()
    index_dir = index_prefix.parent if str(index_prefix) else (host_root / "index" / _sanitize_host_slug(host_key))
    if fasta_path.is_file():
        fasta_path.unlink(missing_ok=True)
        fai_path = fasta_path.with_suffix(f"{fasta_path.suffix}.fai")
        fai_path.unlink(missing_ok=True)
    download_stems = {
        _sanitize_host_slug(str(host_key or "")),
        _sanitize_host_slug(str(current.get("genome_name") or "")),
        _sanitize_host_slug(str(current.get("host_name") or "")),
        _sanitize_host_slug(str(current.get("source_accession") or "")),
    }
    for stem in {item for item in download_stems if item}:
        for path in downloads_dir.glob(f"{stem}*"):
            if path.is_file():
                path.unlink(missing_ok=True)
    if index_dir.exists():
        shutil.rmtree(index_dir, ignore_errors=True)
    store.delete_host_database_record(host_key)
    return {"status": "deleted", "host_key": host_key}
