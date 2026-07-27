from __future__ import annotations

from pathlib import Path

from flask import jsonify, request, session

from .app_services import get_app_services
from .identity import UserIdentity
from .task_manager import ValidationError


def register_dataset_routes(app) -> None:
    services = get_app_services(app)
    login_required = services.access.login_required
    datasets = services.dataset_service

    def identity() -> UserIdentity:
        return UserIdentity(
            username=str(session.get("username") or ""),
            role=str(session.get("role") or ""),
            group_name=str(session.get("group_name") or ""),
        )

    @app.get("/api/database/nextstrain-builds")
    @login_required
    def list_database_nextstrain_builds():
        return jsonify({"items": datasets.list_nextstrain_builds()})

    @app.get("/api/database/auspice-uploads")
    @login_required
    def list_database_auspice_uploads():
        return jsonify({"items": datasets.list_auspice_uploads()})

    @app.post("/api/database/auspice-uploads")
    @login_required
    def create_database_auspice_upload():
        upload = request.files.get("file")
        if upload is None or not str(upload.filename or "").strip():
            raise ValidationError("请先选择一个 Auspice JSON 文件")
        original_filename = Path(str(upload.filename or "").strip()).name
        item = datasets.create_auspice_upload(
            raw_bytes=upload.read(),
            original_filename=original_filename,
            name=str(request.form.get("name") or "").strip(),
            identity=identity(),
        )
        return jsonify({"item": item})

    @app.post("/api/database/auspice-uploads/<upload_id>/open-output")
    @login_required
    def open_database_auspice_upload_output(upload_id: str):
        target = datasets.open_auspice_upload(upload_id)
        return jsonify({"status": "opened", "path": str(target)})

    @app.delete("/api/database/auspice-uploads/<upload_id>")
    @login_required
    def delete_database_auspice_upload(upload_id: str):
        datasets.delete_auspice_upload(upload_id)
        return jsonify({"status": "deleted", "upload_id": upload_id})

    @app.post("/api/database/nextstrain-builds")
    @login_required
    def create_database_nextstrain_build():
        task, build = datasets.create_nextstrain_build(request.get_json(force=True), identity=identity())
        return jsonify({"task": task, "build": build})
