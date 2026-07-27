from __future__ import annotations

from pathlib import Path

from flask import Response, jsonify, request, session

from .app_services import get_app_services
from .identity import UserIdentity
from .import_templates import _build_reference_import_template_text, _build_reference_import_template_xlsx
from .task_manager import ValidationError


def register_reference_database_routes(app) -> None:
    services = get_app_services(app)
    login_required = services.access.login_required
    admin_required = services.access.admin_required
    database = services.reference_database_service
    batch_imports = services.batch_import_service

    def owner() -> str:
        return str(session.get("username") or "")

    def identity() -> UserIdentity:
        return UserIdentity(
            username=owner(),
            role=str(session.get("role") or ""),
            group_name=str(session.get("group_name") or ""),
        )

    def import_local(category: str):
        payload = request.get_json(force=True)
        return jsonify(
            database.import_from_source(
                category=category,
                host_name=str(payload.get("host_name") or "").strip(),
                genome_name=str(payload.get("genome_name") or "").strip(),
                taxid=str(payload.get("taxid") or "").strip(), source_type="local",
                source_label=str(payload.get("source_label") or "本地导入").strip(),
                source_accession=str(payload.get("source_accession") or "").strip(), source_url="",
                source_path_value=payload.get("fasta_path"), description=str(payload.get("description") or "").strip(),
                owner=owner(),
            )
        )

    def upload_local(category: str):
        upload = request.files.get("file")
        if upload is None or not upload.filename:
            raise ValidationError("请先上传 FASTA 文件")
        uploaded_path = database.save_uploaded_genome(category=category, upload=upload)
        return jsonify(
            database.register_record(
                category=category,
                host_name=str(request.form.get("host_name") or "").strip(),
                genome_name=str(request.form.get("genome_name") or "").strip() or Path(upload.filename).stem,
                taxid=str(request.form.get("taxid") or "").strip(), source_type="local",
                source_label=str(request.form.get("source_label") or "本地上传").strip(),
                source_accession=str(request.form.get("source_accession") or "").strip(), source_url="",
                source_fasta_path=uploaded_path, description=str(request.form.get("description") or "").strip(), owner=owner(),
            )
        )

    def import_batch(category: str):
        return jsonify(
            batch_imports.import_reference(
                category=category,
                precheck_id=str(request.form.get("precheck_id") or "").strip(),
                identity=identity(),
            )
        )

    def remote_import(category: str):
        payload = request.get_json(force=True)
        provider = str(payload.get("provider") or "").strip().lower()
        allowed_providers = {"ncbi", "ensembl"} if category == "host" else {"ncbi", "gtdb"}
        allowed_modes = {"datasets", "api"} if category == "host" else {"datasets", "api", "taxid_refs"}
        if provider not in allowed_providers:
            raise ValidationError(f"provider 只支持 {' 或 '.join(sorted(allowed_providers))}")
        ncbi_mode = str(payload.get("ncbi_mode") or "datasets").strip().lower()
        if ncbi_mode not in allowed_modes:
            raise ValidationError(f"ncbi_mode 只支持 {'、'.join(sorted(allowed_modes))}")
        query = str(payload.get("query") or "").strip()
        if not query:
            raise ValidationError("请填写 accession 或下载链接")
        if category == "pathogen" and provider == "ncbi" and ncbi_mode == "taxid_refs":
            if not query.isdigit():
                raise ValidationError("TaxID 下载模式下请输入纯数字 TaxID")
            payload["taxid"] = query
            try:
                payload["max_genomes"] = max(1, min(200, int(payload.get("max_genomes") or 20)))
            except (TypeError, ValueError) as exc:
                raise ValidationError("最多导入参考数必须是 1-200 的整数") from exc
        return jsonify({"status": "queued", "job": database.start_download_job(category=category, payload=payload, owner=owner())})

    @app.get("/api/reference-database/import-template")
    @login_required
    def download_reference_database_import_template():
        category = str(request.args.get("category", "host") or "host").strip().lower()
        if category not in {"host", "pathogen"}:
            raise ValidationError("category 只支持 host 或 pathogen")
        file_format = str(request.args.get("format", "xlsx") or "xlsx").strip().lower()
        if file_format not in {"xlsx", "csv", "tsv"}:
            raise ValidationError("模板格式只支持 xlsx、csv 或 tsv")
        if file_format == "xlsx":
            content = _build_reference_import_template_xlsx(category)
            mimetype, suffix = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", "xlsx"
        else:
            delimiter = "," if file_format == "csv" else "\t"
            content = _build_reference_import_template_text(category, delimiter).encode("utf-8-sig")
            mimetype, suffix = ("text/csv" if file_format == "csv" else "text/tab-separated-values"), file_format
        filename = "host_reference_import_template" if category == "host" else "pathogen_reference_import_template"
        return Response(content, mimetype=mimetype, headers={"Content-Disposition": f'attachment; filename="{filename}.{suffix}"'})

    @app.get("/api/host-database/records")
    @login_required
    def list_host_database_records():
        return jsonify({"items": database.list_records("host")})

    @app.get("/api/host-database/remote-import-jobs")
    @login_required
    def list_host_database_remote_import_jobs():
        return jsonify({"items": database.list_download_jobs("host")})

    @app.delete("/api/host-database/remote-import-jobs/<job_id>")
    @admin_required
    def delete_host_database_remote_import_job(job_id: str):
        database.delete_download_job(job_id, "host")
        return jsonify({"status": "deleted", "job_id": job_id})

    @app.get("/api/host-database/records/<path:host_key>")
    @login_required
    def get_host_database_record(host_key: str):
        return jsonify(database.get_record(host_key))

    @app.get("/api/pathogen-database/records")
    @login_required
    def list_pathogen_database_records():
        return jsonify({"items": database.list_records("pathogen")})

    @app.get("/api/pathogen-database/remote-import-jobs")
    @login_required
    def list_pathogen_database_remote_import_jobs():
        return jsonify({"items": database.list_download_jobs("pathogen")})

    @app.delete("/api/pathogen-database/remote-import-jobs/<job_id>")
    @admin_required
    def delete_pathogen_database_remote_import_job(job_id: str):
        database.delete_download_job(job_id, "pathogen")
        return jsonify({"status": "deleted", "job_id": job_id})

    @app.get("/api/pathogen-database/records/<path:host_key>")
    @login_required
    def get_pathogen_database_record(host_key: str):
        return jsonify(database.get_record(host_key))

    @app.get("/api/pathogen-database/cgmlst-panels")
    @login_required
    def list_pathogen_cgmlst_panels():
        return jsonify({"items": database.list_panels("pathogen", "cgmlst")})

    @app.post("/api/task-references/upload")
    @login_required
    def upload_task_reference():
        upload = request.files.get("file")
        if upload is None or not upload.filename:
            raise ValidationError("请先上传参考 FASTA 文件")
        saved_path = database.save_uploaded_task_reference(upload=upload, owner=str(session.get("username") or ""))
        return jsonify({"status": "ok", "path": str(saved_path), "name": Path(str(upload.filename or "")).name})

    @app.post("/api/host-database/local-import")
    @admin_required
    def import_local_host_genome():
        return import_local("host")

    @app.post("/api/pathogen-database/local-import")
    @admin_required
    def import_local_pathogen_genome():
        return import_local("pathogen")

    @app.post("/api/host-database/upload-import")
    @admin_required
    def upload_local_host_genome():
        return upload_local("host")

    @app.post("/api/pathogen-database/upload-import")
    @admin_required
    def upload_local_pathogen_genome():
        return upload_local("pathogen")

    @app.post("/api/host-database/batch-import")
    @admin_required
    def import_batch_host_genomes():
        return import_batch("host")

    @app.post("/api/pathogen-database/batch-import")
    @admin_required
    def import_batch_pathogen_genomes():
        return import_batch("pathogen")

    @app.post("/api/reference-database/batch-import/precheck")
    @admin_required
    def precheck_batch_reference_genomes():
        category = str(request.form.get("category") or request.args.get("category") or "host").strip().lower()
        if category not in {"host", "pathogen"}:
            raise ValidationError("category 只支持 host 或 pathogen")
        upload = request.files.get("file")
        if upload is None or not upload.filename:
            raise ValidationError("请先选择批量导入文件")
        filename, content = str(upload.filename or ""), upload.read()
        return jsonify(batch_imports.precheck_reference(category=category, filename=filename, content=content, identity=identity()))

    @app.post("/api/host-database/remote-import")
    @admin_required
    def import_remote_host_genome():
        return remote_import("host")

    @app.post("/api/pathogen-database/remote-import")
    @admin_required
    def import_remote_pathogen_genome():
        return remote_import("pathogen")

    def update_record(category: str, host_key: str):
        if request.files:
            payload = dict(request.form)
            upload = request.files.get("file")
            if upload is not None and upload.filename:
                payload["uploaded_fasta_path"] = str(database.save_uploaded_genome(category=category, upload=upload))
        else:
            payload = request.get_json(force=True)
        return jsonify(database.update_record(host_key=host_key, category=category, payload=payload))

    @app.put("/api/host-database/records/<path:host_key>")
    @admin_required
    def update_host_database_record(host_key: str):
        return update_record("host", host_key)

    @app.put("/api/pathogen-database/records/<path:host_key>")
    @admin_required
    def update_pathogen_database_record(host_key: str):
        return update_record("pathogen", host_key)

    @app.post("/api/host-database/build-index")
    @admin_required
    def build_host_database_index():
        return jsonify(database.build_index("host"))

    @app.post("/api/host-database/records/<path:host_key>/build-index")
    @admin_required
    def build_single_host_database_index(host_key: str):
        return jsonify(database.build_single_index(category="host", host_key=host_key))

    @app.delete("/api/host-database/records/<path:host_key>")
    @admin_required
    def delete_host_database_record(host_key: str):
        return jsonify(database.delete_record(category="host", host_key=host_key))

    @app.post("/api/pathogen-database/build-index")
    @admin_required
    def build_pathogen_database_index():
        payload = request.get_json(silent=True) or {}
        host_keys = payload.get("host_keys") or []
        if not host_keys:
            return jsonify(database.build_index("pathogen"))
        return jsonify(database.build_selected_indexes(host_keys))

    @app.post("/api/pathogen-database/cgmlst-panels/build")
    @admin_required
    def build_pathogen_cgmlst_panels():
        payload = request.get_json(force=True) or {}
        return jsonify(
            database.queue_cgmlst_panels(
                host_keys=payload.get("host_keys") or [],
                threshold=str(payload.get("threshold") or "0.95").strip() or "0.95",
                threads=max(1, int(payload.get("threads") or 8)), owner=str(session.get("username") or ""),
            )
        )

    @app.post("/api/pathogen-database/records/<path:host_key>/build-index")
    @admin_required
    def build_single_pathogen_database_index(host_key: str):
        return jsonify(database.build_single_index(category="pathogen", host_key=host_key))

    @app.delete("/api/pathogen-database/records/<path:host_key>")
    @admin_required
    def delete_pathogen_database_record(host_key: str):
        return jsonify(database.delete_record(category="pathogen", host_key=host_key))
