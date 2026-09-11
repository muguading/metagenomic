from __future__ import annotations

import csv
import io
import re
import secrets
import time
import zipfile
from datetime import datetime
from pathlib import Path
from urllib.parse import quote

from flask import Response, jsonify, request, send_file, session

from .app_services import get_app_services
from .export_utils import (
    _build_delimited_bytes,
    _build_xlsx_bytes,
    _build_xlsx_workbook_bytes,
    _normalize_export_columns,
    _normalize_export_rows,
    _normalize_export_sheets,
    _sanitize_export_filename,
)
from .import_templates import _extract_batch_upload_headers, _parse_database_batch_upload
from .task_manager import ValidationError


ANALYSIS_EXPORT_TYPES = {
    "fasta": {
        "label": "FASTA",
        "patterns": ("{sample}.final.fasta",),
    },
    "qc": {
        "label": "质控结果",
        "patterns": ("{sample}.fastp2.json", "{sample}.fastp2.log"),
    },
    "assembly": {
        "label": "组装结果",
        "patterns": ("Assem_info.tsv", "Assem_info1.tsv", "{sample}.checkm.tsv"),
    },
    "resistance_virulence": {
        "label": "耐药毒力结果",
        "patterns": ("Assem_abricate_CARD.tsv", "Assem_abricate_VFDB.tsv"),
    },
    "mge": {
        "label": "MGE结果",
        "patterns": ("*MGE*", "*mge*", "*genomad*", "*geNomad*"),
        "any_match": True,
    },
    "serotype": {
        "label": "血清型鉴定结果",
        "patterns": ("{sample}_serotype_result.tsv", "{sample}.keblo.tsv", "{sample}.pathonet_result.tsv"),
    },
}


def _tsv_bytes(columns: list[str], rows: list[list[object]]) -> bytes:
    buffer = io.StringIO()
    writer = csv.writer(buffer, delimiter="\t", lineterminator="\n")
    writer.writerow(columns)
    writer.writerows(rows)
    return buffer.getvalue().encode("utf-8-sig")


def _safe_zip_part(value: object, fallback: str = "unknown") -> str:
    text = str(value or "").strip() or fallback
    return "".join(ch if ch.isalnum() or ch in ("-", "_", ".") else "_" for ch in text).strip("._") or fallback


def _replace_archive_sample_name(filename: str, source_sample_name: str, export_sample_name: str) -> str:
    """Replace a source sample-name prefix in an artifact filename when present."""
    source = _safe_zip_part(source_sample_name, "")
    target = _safe_zip_part(export_sample_name)
    if source and target != source:
        if filename.startswith(f"{source}."):
            return f"{target}{filename[len(source):]}"
        if filename.startswith(f"{source}_"):
            return f"{target}{filename[len(source):]}"
    return filename


def _analysis_archive_name(
    artifact_type: str,
    sample_name: str,
    path: Path,
    task_id: str = "",
    *,
    source_sample_name: str = "",
    include_task_prefix: bool = False,
) -> str:
    safe_sample = _safe_zip_part(sample_name)
    safe_task = _safe_zip_part(task_id, "")
    filename = _replace_archive_sample_name(path.name, source_sample_name, sample_name)
    if safe_sample and not filename.startswith(f"{safe_sample}.") and not filename.startswith(f"{safe_sample}_"):
        filename = f"{safe_sample}_{filename}"
    if include_task_prefix and safe_task and not filename.startswith(f"{safe_task}_"):
        filename = f"{safe_task}_{filename}"
    return str(Path(artifact_type) / filename)


def _duplicate_sample_names(raw_samples: list[object]) -> set[str]:
    sample_tasks: dict[str, set[str]] = {}
    for item in raw_samples:
        if not isinstance(item, dict):
            continue
        sample_name = str(item.get("sample_name") or "").strip()
        task_id = str(item.get("task_id") or "").strip()
        if not sample_name or not task_id:
            continue
        sample_tasks.setdefault(sample_name, set()).add(task_id)
    return {sample_name for sample_name, task_ids in sample_tasks.items() if len(task_ids) > 1}


def _export_sample_name(sample_item: dict) -> str:
    return str(sample_item.get("export_sample_name") or sample_item.get("sample_name") or "").strip()


def _duplicate_export_sample_names(raw_samples: list[object]) -> set[str]:
    return _duplicate_sample_names(
        [
            {
                "task_id": item.get("task_id"),
                "sample_name": _export_sample_name(item),
            }
            for item in raw_samples
            if isinstance(item, dict)
        ]
    )


def _is_ncov_or_monkeypox_task(task: dict) -> bool:
    """Identify virus tasks whose exported consensus FASTA header follows the sample label."""
    params = task.get("params") if isinstance(task.get("params"), dict) else {}
    values = [
        task.get("name"),
        task.get("demo_type"),
        task.get("pipeline_script"),
        params.get("species"),
        params.get("ref"),
    ]
    text = " ".join(str(value or "").lower() for value in values)
    markers = (
        "sars-cov-2",
        "sars cov 2",
        "2019-ncov",
        "covid-19",
        "ncov",
        "新冠",
        "新型冠状",
        "monkeypox",
        "mpox",
        "hmpxv",
        "猴痘",
    )
    return any(marker in text for marker in markers)


def _rewrite_first_fasta_contig_name(content: bytes, sample_name: str) -> bytes:
    """Return a FASTA copy whose first record identifier is the safe exported sample name."""
    target = _safe_zip_part(sample_name).encode("utf-8")
    header_start = content.find(b">")
    while header_start >= 0 and header_start not in {0} and content[header_start - 1:header_start] not in {b"\n", b"\r"}:
        header_start = content.find(b">", header_start + 1)
    if header_start < 0:
        return content
    line_end = content.find(b"\n", header_start)
    if line_end < 0:
        line_end = len(content)
    old_header = content[header_start + 1:line_end].rstrip(b"\r")
    description_start = len(old_header)
    for index, character in enumerate(old_header):
        if character in b" \t":
            description_start = index
            break
    description = old_header[description_start:]
    return content[:header_start] + b">" + target + description + content[line_end:]


_ILLUMINA_LIBRARY_SUFFIX_RE = re.compile(r"_S\d+_L\d{3}_\d{3}$", re.IGNORECASE)


def _canonical_meta_sequencing_name(value: object) -> str:
    """Remove the standard Illumina sample/lane/read suffix used in FASTQ folder names."""
    return _ILLUMINA_LIBRARY_SUFFIX_RE.sub("", str(value or "").strip()).strip()


def _prepare_meta_name_mapping(rows: list[dict[str, str]], sequencing_column: str, sample_column: str) -> dict[str, object]:
    mapping: dict[str, str] = {}
    duplicates: set[str] = set()
    skipped_rows = 0
    for row in rows:
        sequencing_name = str(row.get(sequencing_column) or "").strip()
        sample_name = str(row.get(sample_column) or "").strip()
        if not sequencing_name or not sample_name:
            skipped_rows += 1
            continue
        existing = mapping.get(sequencing_name)
        if existing is not None and existing != sample_name:
            duplicates.add(sequencing_name)
            continue
        mapping[sequencing_name] = sample_name
    for sequencing_name in duplicates:
        mapping.pop(sequencing_name, None)

    normalized_mapping: dict[str, str] = {}
    ambiguous_normalized_names: set[str] = set()
    for sequencing_name, sample_name in mapping.items():
        normalized_name = _canonical_meta_sequencing_name(sequencing_name)
        if not normalized_name:
            continue
        existing = normalized_mapping.get(normalized_name)
        if existing is not None and existing != sample_name:
            ambiguous_normalized_names.add(normalized_name)
            continue
        normalized_mapping[normalized_name] = sample_name
    for normalized_name in ambiguous_normalized_names:
        normalized_mapping.pop(normalized_name, None)

    return {
        "mapping": mapping,
        "normalized_mapping": normalized_mapping,
        "mapped_count": len(mapping),
        "library_suffix_normalized_count": len(normalized_mapping),
        "skipped_rows": skipped_rows,
        "duplicate_sequencing_names": sorted(duplicates),
        "ambiguous_normalized_sequencing_names": sorted(ambiguous_normalized_names),
    }


def _iter_report_sample_dirs(report_source: dict, requested_sample: str) -> list[tuple[str, Path]]:
    mode = str(report_source.get("mode") or "single")
    report_dir = Path(report_source["report_dir"]).resolve()
    requested = str(requested_sample or "").strip()
    if mode == "multi":
        root_dir = Path(report_source.get("root_dir") or report_dir).resolve()
        if requested:
            sample_dir = (root_dir / requested).resolve()
            if sample_dir.is_dir() and root_dir in sample_dir.parents:
                return [(requested, sample_dir)]
        selected = str(report_source.get("selected_sample") or "").strip() or report_dir.name
        return [(selected, report_dir)]
    sample_name = requested or str(report_source.get("selected_sample") or "").strip() or report_dir.name
    return [(sample_name, report_dir)]


def _find_analysis_artifacts(sample_dir: Path, sample_name: str, artifact_type: str) -> tuple[list[Path], list[str]]:
    spec = ANALYSIS_EXPORT_TYPES[artifact_type]
    found: list[Path] = []
    missing: list[str] = []
    seen: set[Path] = set()
    pattern_labels: list[str] = []
    for pattern_template in spec["patterns"]:
        pattern = pattern_template.format(sample=sample_name)
        pattern_labels.append(pattern)
        matches: list[Path] = []
        for candidate in sample_dir.glob(pattern):
            resolved = candidate.resolve()
            if resolved.is_file():
                matches.append(resolved)
            elif resolved.is_dir():
                matches.extend(path.resolve() for path in resolved.rglob("*") if path.is_file())
        matches = sorted(matches, key=lambda item: str(item).lower())
        safe_matches = [path for path in matches if sample_dir == path.parent or sample_dir in path.parents]
        if not safe_matches:
            if not spec.get("any_match"):
                missing.append(pattern)
            continue
        for path in safe_matches:
            if path in seen:
                continue
            seen.add(path)
            found.append(path)
    if spec.get("any_match") and not found:
        missing.append(" / ".join(pattern_labels))
    return found, missing


def register_export_routes(app) -> None:
    services = get_app_services(app)
    login_required = services.access.login_required
    metadata_import_cache: dict[str, dict[str, object]] = {}

    def _purge_expired_metadata_imports() -> None:
        expires_before = time.monotonic() - 30 * 60
        for import_id, item in list(metadata_import_cache.items()):
            if float(item.get("created_at", 0) or 0) < expires_before:
                metadata_import_cache.pop(import_id, None)

    def _get_metadata_import(import_id: str) -> dict[str, object]:
        _purge_expired_metadata_imports()
        item = metadata_import_cache.get(import_id)
        if not item or str(item.get("username") or "") != str(session.get("username") or ""):
            raise ValidationError("Meta 导入记录不存在或已过期，请重新导入文件")
        return item

    @app.post("/api/export/table")
    @login_required
    def export_table():
        payload = request.get_json(force=True)
        title = str(payload.get("title", "")).strip() or "结果表"
        export_format = str(payload.get("format", "csv")).strip().lower()
        columns = _normalize_export_columns(payload.get("columns", []))
        rows = _normalize_export_rows(payload.get("rows", []))
        sheets = _normalize_export_sheets(payload.get("sheets", []))
        filename_root = _sanitize_export_filename(str(payload.get("filename", "")).strip() or title)

        if export_format == "csv":
            content = _build_delimited_bytes(columns, rows, ",")
            mimetype = "text/csv; charset=utf-8"
            extension = "csv"
        elif export_format == "tsv":
            content = _build_delimited_bytes(columns, rows, "\t")
            mimetype = "text/tab-separated-values; charset=utf-8"
            extension = "tsv"
        elif export_format == "xlsx":
            content = _build_xlsx_workbook_bytes(sheets) if sheets else _build_xlsx_bytes(title, columns, rows)
            mimetype = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            extension = "xlsx"
        else:
            raise ValidationError("仅支持导出 csv、tsv 或 xlsx")

        filename = f"{filename_root}.{extension}"
        return Response(
            content,
            mimetype=mimetype,
            headers={
                "Content-Disposition": f"attachment; filename*=UTF-8''{quote(filename)}",
                "Cache-Control": "no-store",
            },
        )

    @app.post("/api/export/sample-meta/preview")
    @login_required
    def preview_sample_meta_import():
        upload = request.files.get("file")
        if upload is None or not str(upload.filename or "").strip():
            raise ValidationError("请先选择 Meta 文件")
        filename = Path(str(upload.filename or "")).name
        if Path(filename).suffix.lower() not in {".csv", ".tsv", ".xlsx"}:
            raise ValidationError("Meta 文件只支持 csv、tsv 或 xlsx")
        content = upload.read()
        if len(content) > 10 * 1024 * 1024:
            raise ValidationError("Meta 文件不能超过 10 MB")
        headers = _extract_batch_upload_headers(filename, content)
        rows = _parse_database_batch_upload(filename, content)
        if not headers:
            raise ValidationError("Meta 文件缺少表头")
        if not rows:
            raise ValidationError("Meta 文件没有可用于匹配的数据行")
        import_id = secrets.token_urlsafe(18)
        _purge_expired_metadata_imports()
        metadata_import_cache[import_id] = {
            "created_at": time.monotonic(),
            "username": str(session.get("username") or ""),
            "filename": filename,
            "headers": headers,
            "rows": rows,
        }
        return jsonify({
            "import_id": import_id,
            "filename": filename,
            "headers": headers,
            "row_count": len(rows),
        })

    @app.post("/api/export/sample-meta/mapping")
    @login_required
    def prepare_sample_meta_mapping():
        payload = request.get_json(force=True) or {}
        import_id = str(payload.get("import_id") or "").strip()
        sequencing_column = str(payload.get("sequencing_column") or "").strip()
        sample_column = str(payload.get("sample_column") or "").strip()
        item = _get_metadata_import(import_id)
        headers = item.get("headers") if isinstance(item.get("headers"), list) else []
        if sequencing_column not in headers or sample_column not in headers:
            raise ValidationError("请选择 Meta 文件中的测序名称列和样本名称列")
        if sequencing_column == sample_column:
            raise ValidationError("测序名称列和样本名称列不能相同")
        rows = item.get("rows") if isinstance(item.get("rows"), list) else []
        result = _prepare_meta_name_mapping(rows, sequencing_column, sample_column)
        return jsonify({
            "filename": item.get("filename"),
            "sequencing_column": sequencing_column,
            "sample_column": sample_column,
            **result,
        })

    @app.post("/api/tasks/batch-analysis-export")
    @login_required
    def export_batch_analysis_results():
        payload = request.get_json(force=True) or {}
        raw_samples = payload.get("samples")
        if not isinstance(raw_samples, list) or not raw_samples:
            raise ValidationError("请选择至少一个样本")
        raw_types = payload.get("artifact_types")
        if not isinstance(raw_types, list) or not raw_types:
            raise ValidationError("请选择至少一种分析结果")
        artifact_types = []
        for item in raw_types:
            key = str(item or "").strip()
            if key not in ANALYSIS_EXPORT_TYPES:
                raise ValidationError(f"不支持的分析结果类型: {key or '-'}")
            if key not in artifact_types:
                artifact_types.append(key)

        exported_at = datetime.now().strftime("%Y%m%d_%H%M%S")
        zip_buffer = io.BytesIO()
        manifest_rows: list[list[object]] = []
        missing_rows: list[list[object]] = []
        seen_requests: set[tuple[str, str]] = set()
        zip_names: set[str] = set()
        duplicate_sample_names = _duplicate_export_sample_names(raw_samples)

        with zipfile.ZipFile(zip_buffer, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            for sample_item in raw_samples:
                if not isinstance(sample_item, dict):
                    continue
                task_id = str(sample_item.get("task_id") or "").strip()
                sample_name = str(sample_item.get("sample_name") or "").strip()
                export_sample_name = _export_sample_name(sample_item)
                if not task_id:
                    continue
                request_key = (task_id, sample_name)
                if request_key in seen_requests:
                    continue
                seen_requests.add(request_key)

                task = services.task_lifecycle_service.get_visible(task_id, log_lines=0)
                if str(task.get("status") or "").upper() != "SUCCEEDED":
                    raise ValidationError(f"{task.get('name') or task_id} 尚未完成，不能导出分析结果")
                report_source = services.reports.resolve_report_source(task, sample_name)
                if not report_source.get("available"):
                    missing_rows.append([
                        task_id,
                        task.get("name") or task_id,
                        sample_name or "-",
                        "结果目录",
                        str(report_source.get("reason") or "无法定位结果目录"),
                    ])
                    continue

                for resolved_sample, sample_dir in _iter_report_sample_dirs(report_source, sample_name):
                    source_sample_label = sample_name or resolved_sample or sample_dir.name
                    sample_label = export_sample_name or source_sample_label
                    for artifact_type in artifact_types:
                        found_paths, missing_patterns = _find_analysis_artifacts(sample_dir, source_sample_label, artifact_type)
                        type_label = ANALYSIS_EXPORT_TYPES[artifact_type]["label"]
                        for pattern in missing_patterns:
                            missing_rows.append([task_id, task.get("name") or task_id, sample_label, type_label, pattern])
                        for path in found_paths:
                            archive_name = _analysis_archive_name(
                                artifact_type,
                                sample_label,
                                path,
                                task_id,
                                source_sample_name=source_sample_label,
                                include_task_prefix=sample_label in duplicate_sample_names,
                            )
                            if archive_name in zip_names:
                                archive_path = Path(archive_name)
                                archive_name = str(archive_path.with_name(f"{len(zip_names)}_{archive_path.name}"))
                            zip_names.add(archive_name)
                            if artifact_type == "fasta" and sample_label != source_sample_label and _is_ncov_or_monkeypox_task(task):
                                exported_content = _rewrite_first_fasta_contig_name(path.read_bytes(), sample_label)
                                archive.writestr(archive_name, exported_content)
                                exported_size = len(exported_content)
                            else:
                                archive.write(path, archive_name)
                                exported_size = path.stat().st_size
                            manifest_rows.append([
                                task_id,
                                task.get("name") or task_id,
                                sample_label,
                                type_label,
                                str(path),
                                archive_name,
                                exported_size,
                            ])

            archive.writestr(
                "manifest.tsv",
                _tsv_bytes(
                    ["任务ID", "任务名称", "样本名称", "结果类型", "原始路径", "压缩包路径", "文件大小"],
                    manifest_rows,
                ),
            )
            archive.writestr(
                "missing_files.tsv",
                _tsv_bytes(
                    ["任务ID", "任务名称", "样本名称", "结果类型", "缺失文件或原因"],
                    missing_rows,
                ),
            )

        zip_buffer.seek(0)
        filename = f"analysis_results_{exported_at}.zip"
        response = send_file(
            zip_buffer,
            mimetype="application/zip",
            as_attachment=True,
            download_name=filename,
            max_age=0,
        )
        response.headers["Content-Disposition"] = f"attachment; filename*=UTF-8''{quote(filename)}"
        response.headers["Cache-Control"] = "no-store"
        return response
