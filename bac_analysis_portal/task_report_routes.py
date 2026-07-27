from __future__ import annotations

from pathlib import Path

from flask import abort, jsonify, request, send_file

from .app_services import get_app_services
from .assembly_reports import _build_cgview_json_from_gbk
from .export_checklist import build_export_checklist
from .report_payload import _build_report_payload
from .runtime_paths import _resolve_runtime_database_root
from .task_analytics import (
    build_pending_queue_analytics_snapshot,
    build_queue_analytics_snapshot,
    read_task_analytics_snapshot,
    write_task_analytics_snapshot,
)
from .task_manager import ValidationError


def register_task_report_routes(app) -> None:
    services = get_app_services(app)
    lifecycle = services.task_lifecycle_service
    login_required = services.access.login_required
    resolve_report_source = services.reports.resolve_report_source
    maybe_auto_trigger_pathosource_for_meta_task = services.reports.maybe_auto_trigger_pathosource_for_meta_task

    @app.get("/api/tasks/<task_id>/report-data")
    @login_required
    def task_report_data(task_id: str):
        task = lifecycle.get_visible(task_id)
        selected_sample = str(request.args.get("sample", "") or "").strip()
        payload = _build_report_payload(task, selected_sample=selected_sample)
        services.portal_service.enrich_task_with_closure_status(task)
        payload_task = payload.get("task") if isinstance(payload.get("task"), dict) else {}
        payload_task["closure_status"] = task.get("closure_status") or {}
        payload["task"] = payload_task
        sections = payload.get("sections") if isinstance(payload.get("sections"), dict) else {}
        sections["export_checklist"] = build_export_checklist(payload)
        payload["sections"] = sections
        payload = maybe_auto_trigger_pathosource_for_meta_task(task, payload, selected_sample=selected_sample)
        return jsonify(payload)

    @app.get("/api/tasks/<task_id>/analytics-snapshot")
    @login_required
    def task_analytics_snapshot(task_id: str):
        task = lifecycle.get_visible(task_id)
        task_dir = lifecycle.task_directory(task_id)
        cached = read_task_analytics_snapshot(task_dir)
        if cached:
            return jsonify(cached)
        if str(task.get("status") or "").upper() != "SUCCEEDED":
            return jsonify(build_pending_queue_analytics_snapshot(task))
        payload = _build_report_payload(task)
        snapshot = build_queue_analytics_snapshot(payload)
        write_task_analytics_snapshot(task_dir, snapshot)
        return jsonify(snapshot)

    @app.get("/api/tasks/<task_id>/report-asset/<path:asset_name>")
    @login_required
    def task_report_asset(task_id: str, asset_name: str):
        task = lifecycle.get_visible(task_id)
        selected_sample = str(request.args.get("sample", "") or "").strip()
        render_as = str(request.args.get("render_as", "") or "").strip().lower()
        report_source = resolve_report_source(task, selected_sample)
        if not report_source.get("available"):
            raise ValidationError(str(report_source.get("reason") or "服务器结果目录尚未就绪。"))
        report_dir = Path(report_source["report_dir"]).resolve()
        database_root = _resolve_runtime_database_root()
        if asset_name == "__ncov_ref.fna":
            ncov_ref = database_root / "virus" / "ncov" / "ref.fna"
            if not ncov_ref.is_file():
                abort(404)
            response = send_file(str(ncov_ref), mimetype="text/plain; charset=utf-8", conditional=True)
            response.headers["Cache-Control"] = "no-store, no-cache, must-revalidate, max-age=0"
            return response
        if asset_name == "__ncov_ref.fna.fai":
            ncov_ref_fai = database_root / "virus" / "ncov" / "ref.fna.fai"
            if not ncov_ref_fai.is_file():
                abort(404)
            response = send_file(str(ncov_ref_fai), mimetype="text/plain; charset=utf-8", conditional=True)
            response.headers["Cache-Control"] = "no-store, no-cache, must-revalidate, max-age=0"
            return response
        if asset_name == "__ncov_genomic.gff":
            ncov_gff = database_root / "virus" / "ncov" / "genomic.gff"
            if not ncov_gff.is_file():
                abort(404)
            response = send_file(str(ncov_gff), mimetype="text/plain; charset=utf-8", conditional=True)
            response.headers["Cache-Control"] = "no-store, no-cache, must-revalidate, max-age=0"
            return response
        if asset_name == "__hmpxv_genome_annotation.gff3":
            hmpxv_db_root = database_root / "virus" / "nextclade" / "hMPXV"
            if not hmpxv_db_root.is_dir():
                hmpxv_db_root = database_root / "nextclade_db" / "hMPXV"
            hmpxv_gff = hmpxv_db_root / "genome_annotation.gff3"
            if not hmpxv_gff.is_file():
                abort(404)
            response = send_file(str(hmpxv_gff), mimetype="text/plain; charset=utf-8", conditional=True)
            response.headers["Cache-Control"] = "no-store, no-cache, must-revalidate, max-age=0"
            return response
        target = (report_dir / asset_name).resolve()
        if report_dir not in target.parents or not target.is_file():
            abort(404)
        target_suffix = target.suffix.lower()
        if target_suffix == ".gbk":
            if render_as == "cgview-json":
                try:
                    target = _build_cgview_json_from_gbk(report_dir, target, target.stem)
                except RuntimeError as exc:
                    raise ValidationError(f"CGView 图谱构建失败: {exc}") from exc
                response = send_file(str(target), mimetype="application/json")
            else:
                response = send_file(str(target), mimetype="text/plain; charset=utf-8")
        elif target_suffix == ".json":
            response = send_file(str(target), mimetype="application/json")
        elif target_suffix in {".bam", ".bai", ".fa", ".fasta", ".fai", ".gff", ".gff3", ".vcf", ".html", ".js", ".svg", ".png", ".csv", ".tsv"}:
            mimetype = {
                ".bam": "application/octet-stream",
                ".bai": "application/octet-stream",
                ".fa": "text/plain; charset=utf-8",
                ".fasta": "text/plain; charset=utf-8",
                ".fai": "text/plain; charset=utf-8",
                ".gff": "text/plain; charset=utf-8",
                ".gff3": "text/plain; charset=utf-8",
                ".vcf": "text/plain; charset=utf-8",
                ".html": "text/html; charset=utf-8",
                ".js": "application/javascript; charset=utf-8",
                ".svg": "image/svg+xml",
                ".png": "image/png",
                ".csv": "text/csv; charset=utf-8",
                ".tsv": "text/tab-separated-values; charset=utf-8",
            }.get(target_suffix, "application/octet-stream")
            response = send_file(str(target), mimetype=mimetype, conditional=True)
        else:
            abort(403)
        response.headers["Cache-Control"] = "no-store, no-cache, must-revalidate, max-age=0"
        response.headers["Pragma"] = "no-cache"
        response.headers["Expires"] = "0"
        return response
