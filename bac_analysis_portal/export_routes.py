from __future__ import annotations

import csv
import io
import zipfile
from datetime import datetime
from pathlib import Path
from urllib.parse import quote

from flask import Response, request, send_file

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


def _analysis_archive_name(artifact_type: str, sample_name: str, path: Path, task_id: str = "", *, include_task_prefix: bool = False) -> str:
    safe_sample = _safe_zip_part(sample_name)
    safe_task = _safe_zip_part(task_id, "")
    filename = path.name
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
        duplicate_sample_names = _duplicate_sample_names(raw_samples)

        with zipfile.ZipFile(zip_buffer, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            for sample_item in raw_samples:
                if not isinstance(sample_item, dict):
                    continue
                task_id = str(sample_item.get("task_id") or "").strip()
                sample_name = str(sample_item.get("sample_name") or "").strip()
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
                    sample_label = sample_name or resolved_sample or sample_dir.name
                    for artifact_type in artifact_types:
                        found_paths, missing_patterns = _find_analysis_artifacts(sample_dir, sample_label, artifact_type)
                        type_label = ANALYSIS_EXPORT_TYPES[artifact_type]["label"]
                        for pattern in missing_patterns:
                            missing_rows.append([task_id, task.get("name") or task_id, sample_label, type_label, pattern])
                        for path in found_paths:
                            archive_name = _analysis_archive_name(
                                artifact_type,
                                sample_label,
                                path,
                                task_id,
                                include_task_prefix=sample_label in duplicate_sample_names,
                            )
                            if archive_name in zip_names:
                                archive_path = Path(archive_name)
                                archive_name = str(archive_path.with_name(f"{len(zip_names)}_{archive_path.name}"))
                            zip_names.add(archive_name)
                            archive.write(path, archive_name)
                            manifest_rows.append([
                                task_id,
                                task.get("name") or task_id,
                                sample_label,
                                type_label,
                                str(path),
                                archive_name,
                                path.stat().st_size,
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
