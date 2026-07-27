from __future__ import annotations

import json

from flask import Response, jsonify, request, session

from .app_services import get_app_services
from .identity import UserIdentity
from .filesystem_helpers import _resolve_optional_existing_path
from .import_templates import (
    _build_database_import_template_text,
    _build_database_import_template_xlsx,
)
from .task_manager import ValidationError


def register_sample_library_routes(app) -> None:
    services = get_app_services(app)
    login_required = services.access.login_required
    admin_required = services.access.admin_required
    batch_imports = services.batch_import_service

    def manager():
        return services.sample_library

    def identity() -> dict[str, str]:
        return {
            "role": str(session.get("role") or ""),
            "username": str(session.get("username") or ""),
            "group_name": str(session.get("group_name") or ""),
        }

    def batch_identity() -> UserIdentity:
        return UserIdentity.from_mapping(identity())

    @app.get("/api/database/samples")
    @login_required
    def list_database_samples():
        return jsonify({"items": manager().list_visible(scope=str(request.args.get("scope", "main") or "main").strip(), **identity())})

    @app.get("/api/database/submissions")
    @login_required
    def list_database_submissions():
        return jsonify({"items": manager().list_submissions(**identity())})

    @app.get("/api/database/version-logs")
    @login_required
    def list_database_version_logs():
        return jsonify({"items": manager().list_version_logs(**identity())})

    @app.get("/api/database/releases")
    @login_required
    def list_database_releases():
        return jsonify({"items": manager().list_release_versions(**identity())})

    @app.get("/api/database/metadata-templates")
    @login_required
    def list_database_metadata_templates():
        return jsonify({"items": manager().list_metadata_templates()})

    @app.put("/api/database/metadata-templates")
    @login_required
    def save_database_metadata_templates():
        try:
            items = manager().save_metadata_templates((request.get_json(force=True).get("items") or []))
        except ValueError as exc:
            raise ValidationError(str(exc)) from exc
        return jsonify({"items": items})

    @app.get("/api/database/samples/<path:sample_key>")
    @login_required
    def get_database_sample(sample_key: str):
        try:
            return jsonify(manager().get_visible(sample_key, **identity()))
        except KeyError as exc:
            raise KeyError(str(exc)) from exc

    @app.put("/api/database/samples/<path:sample_key>")
    @login_required
    def update_database_sample(sample_key: str):
        payload = request.get_json(force=True)
        editable_fields = {
            "genome_id", "pathogen_type", "sample_alias", "taxid", "mlst_st", "serotype_result",
            "resistance_gene_hits", "virulence_gene_hits", "resistance_mge_hits", "virulence_mge_hits",
            "description", "sample_source", "collection_date", "gender", "country", "host_info",
            "location_json", "sample_type", "sequencing_method", "genome_length", "note", "visibility_scope",
            "custom_metadata_json", "metadata_templates",
        }
        try:
            updated = manager().update_visible(
                sample_key,
                **identity(),
                **{key: payload.get(key) for key in editable_fields},
            )
        except (PermissionError, ValueError) as exc:
            raise ValidationError(str(exc)) from exc
        return jsonify(updated)

    @app.post("/api/database/samples/batch-update")
    @login_required
    def batch_update_database_samples():
        payload = request.get_json(force=True)
        try:
            updated = manager().batch_update_visible(
                payload.get("sample_keys") or [],
                **identity(),
                custom_metadata_json=payload.get("custom_metadata_json"),
                metadata_templates=payload.get("metadata_templates"),
            )
        except (PermissionError, ValueError) as exc:
            raise ValidationError(str(exc)) from exc
        return jsonify(updated)

    @app.delete("/api/database/samples/<path:sample_key>")
    @login_required
    def delete_database_sample(sample_key: str):
        try:
            manager().delete_visible(sample_key, **identity())
        except PermissionError as exc:
            raise ValidationError(str(exc)) from exc
        return jsonify({"status": "deleted", "sample_key": sample_key})

    @app.post("/api/database/samples/<path:sample_key>/submit")
    @login_required
    def submit_database_sample(sample_key: str):
        try:
            return jsonify(manager().submit_personal_to_main(sample_key, **identity()))
        except PermissionError as exc:
            raise ValidationError(str(exc)) from exc

    @app.post("/api/database/submissions/<request_id>/review")
    @admin_required
    def review_database_submission(request_id: str):
        payload = request.get_json(force=True)
        return jsonify(
            manager().review_submission(
                request_id,
                action=str(payload.get("action") or "").strip(),
                admin_username=str(session.get("username") or ""),
                note=str(payload.get("note") or "").strip(),
            )
        )

    @app.post("/api/database/releases")
    @admin_required
    def publish_database_release():
        payload = request.get_json(force=True)
        try:
            return jsonify(
                manager().publish_release_version(
                    version_label=str(payload.get("version_label") or "").strip(),
                    note=str(payload.get("note") or "").strip(),
                    **identity(),
                )
            )
        except (PermissionError, ValueError) as exc:
            raise ValidationError(str(exc)) from exc

    @app.post("/api/database/local-import")
    @login_required
    def import_local_sample_to_database():
        payload = request.get_json(force=True)
        requested_scope = str(payload.get("library_scope") or "").strip()
        is_admin = str(session.get("role") or "") == "admin"
        library_scope = requested_scope if is_admin and requested_scope in {"main", "personal"} else ("main" if is_admin else "personal")
        text_fields = (
            "sample_name", "pathogen_type", "species_name", "mlst_species_name", "mlst_st", "serotype_result",
            "q20_rate", "q30_rate", "completeness", "contamination", "contig_count", "plasmid_count",
            "resistance_count", "virulence_count", "resistance_gene_hits", "virulence_gene_hits",
            "resistance_mge_hits", "virulence_mge_hits", "genome_id", "taxid", "description", "sample_source",
            "collection_date", "gender", "country", "host_info", "location_json", "sample_type",
            "sequencing_method", "note",
        )
        return jsonify(
            manager().import_local_sample(
                owner=str(session.get("username") or ""),
                owner_group=str(session.get("group_name") or ""),
                library_scope=library_scope,
                final_fasta_path=_resolve_optional_existing_path(services.project_root, payload.get("final_fasta_path")),
                custom_metadata_json=json.dumps(payload.get("custom_metadata_json") or [], ensure_ascii=False),
                **{key: str(payload.get(key) or "").strip() for key in text_fields},
            )
        )

    @app.get("/api/batch-import-runs")
    @login_required
    def list_batch_import_runs():
        items = batch_imports.list_runs(
            identity=batch_identity(),
            import_type=str(request.args.get("import_type") or "").strip(),
            category=str(request.args.get("category") or "").strip(),
            limit=request.args.get("limit", default=100, type=int),
        )
        return jsonify({"items": items})

    @app.post("/api/database/batch-import")
    @login_required
    def import_batch_samples_to_database():
        return jsonify(
            batch_imports.import_samples(
                precheck_id=str(request.form.get("precheck_id") or "").strip(),
                requested_scope=str(request.form.get("library_scope") or request.args.get("library_scope") or "").strip(),
                identity=batch_identity(),
            )
        )

    @app.post("/api/database/batch-import/precheck")
    @login_required
    def precheck_batch_samples_to_database():
        upload = request.files.get("file")
        if upload is None or not upload.filename:
            raise ValidationError("请先选择批量导入文件")
        filename, content = str(upload.filename or ""), upload.read()
        return jsonify(batch_imports.precheck_samples(filename=filename, content=content, identity=batch_identity()))

    @app.get("/api/database/import-template")
    @login_required
    def download_database_import_template():
        file_format = str(request.args.get("format", "xlsx") or "xlsx").strip().lower()
        if file_format not in {"xlsx", "csv", "tsv"}:
            raise ValidationError("模板格式只支持 xlsx、csv 或 tsv")
        if file_format == "xlsx":
            content, mimetype, suffix = _build_database_import_template_xlsx(manager().list_metadata_templates()), "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", "xlsx"
        else:
            content = _build_database_import_template_text("," if file_format == "csv" else "\t").encode("utf-8-sig")
            mimetype, suffix = ("text/csv" if file_format == "csv" else "text/tab-separated-values"), file_format
        return Response(content, mimetype=mimetype, headers={"Content-Disposition": f'attachment; filename="sample_database_import_template.{suffix}"'})
