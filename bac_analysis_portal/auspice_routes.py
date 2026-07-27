from __future__ import annotations

from pathlib import Path

from flask import Response, abort, jsonify, request, send_file, send_from_directory, url_for

from .app_services import get_app_services
def register_auspice_routes(app) -> None:
    services = get_app_services(app)
    project_root = services.project_root
    lifecycle = services.task_lifecycle_service
    login_required = services.access.login_required
    nextstrain_build_task_type = services.auspice.nextstrain_build_task_type
    add_no_cache_headers = services.auspice.add_no_cache_headers
    normalize_nextstrain_dataset_prefix = services.auspice.normalize_nextstrain_dataset_prefix
    load_uploaded_auspice_metadata = services.auspice.load_uploaded_auspice_metadata
    ensure_can_view_uploaded_auspice = services.access.ensure_can_view_uploaded_auspice
    uploaded_auspice_dataset_path = services.auspice.uploaded_auspice_dataset_path
    resolve_nextstrain_dataset_main_json = services.auspice.resolve_nextstrain_dataset_main_json
    resolve_nextstrain_dataset_sidecar = services.auspice.resolve_nextstrain_dataset_sidecar

    def _auspice_index_response() -> Response:
        build_root = (project_root / "public" / "auspice-us" / "dist").resolve()
        index_path = build_root / "index.html"
        if not index_path.is_file():
            abort(503, description="本地 Auspice 构建尚未就绪。")
        html = index_path.read_text(encoding="utf-8", errors="ignore")
        html = html.replace("/favicon.png", url_for("static", filename="app_icon.png"))
        return add_no_cache_headers(Response(html, mimetype="text/html"))

    @app.get("/public/<path:filename>")
    @login_required
    def public_asset(filename: str):
        public_root = (project_root / "public").resolve()
        response = send_from_directory(str(public_root), filename)
        response.headers["Cache-Control"] = "no-store, no-cache, must-revalidate, max-age=0"
        response.headers["Pragma"] = "no-cache"
        response.headers["Expires"] = "0"
        return response

    @app.get("/auspice")
    @app.get("/auspice/")
    @login_required
    def local_auspice_page():
        return _auspice_index_response()

    @app.get("/uploaded/<upload_id>")
    @app.get("/uploaded/<upload_id>/")
    @app.get("/auspice/uploaded/<upload_id>")
    @app.get("/auspice/uploaded/<upload_id>/")
    @login_required
    def uploaded_auspice_page(upload_id: str):
        metadata = load_uploaded_auspice_metadata(project_root, upload_id)
        ensure_can_view_uploaded_auspice(metadata)
        dataset_path = uploaded_auspice_dataset_path(project_root, upload_id)
        if not dataset_path.is_file():
            abort(404)
        return _auspice_index_response()

    @app.get("/nextstrain/<task_id>")
    @app.get("/nextstrain/<task_id>/")
    @login_required
    def nextstrain_dataset_page(task_id: str):
        try:
            task = lifecycle.get_visible(task_id)
        except KeyError:
            abort(404)
        if str(task.get("task_type") or "").strip() != nextstrain_build_task_type:
            abort(404)
        return _auspice_index_response()

    @app.get("/dist/<path:filename>")
    @login_required
    def local_auspice_asset(filename: str):
        build_root = (project_root / "public" / "auspice-us" / "dist").resolve()
        return add_no_cache_headers(send_from_directory(str(build_root), filename))

    def _local_auspice_available_response():
        prefix = str(request.args.get("prefix", "") or "").strip()
        resolved = normalize_nextstrain_dataset_prefix(prefix)
        if not resolved:
            return jsonify({"datasets": []})
        resolved_type, dataset_id = resolved
        if resolved_type == "task":
            try:
                task = lifecycle.get_visible(dataset_id)
            except KeyError:
                return jsonify({"datasets": []})
            if str(task.get("task_type") or "").strip() != nextstrain_build_task_type:
                return jsonify({"datasets": []})
            dataset_json = resolve_nextstrain_dataset_main_json(project_root, task)
            normalized_prefix = f"nextstrain/{dataset_id}"
        elif resolved_type == "uploaded":
            metadata = load_uploaded_auspice_metadata(project_root, dataset_id)
            ensure_can_view_uploaded_auspice(metadata)
            dataset_json = uploaded_auspice_dataset_path(project_root, dataset_id)
            normalized_prefix = f"uploaded/{dataset_id}"
        else:
            return jsonify({"datasets": []})
        if not dataset_json.is_file():
            return jsonify({"datasets": []})
        return jsonify({"datasets": [{"request": normalized_prefix, "build": normalized_prefix}]})

    @app.get("/getAvailable")
    @app.get("/charon/getAvailable")
    @login_required
    def local_auspice_available():
        return _local_auspice_available_response()

    def _local_auspice_dataset_response():
        prefix = str(request.args.get("prefix", "") or "").strip()
        dataset_type = str(request.args.get("type", "") or "").strip().lower()
        resolved = normalize_nextstrain_dataset_prefix(prefix)
        if not resolved:
            abort(404)
        resolved_type, dataset_id = resolved
        if resolved_type == "task":
            try:
                task = lifecycle.get_visible(dataset_id)
            except KeyError:
                abort(404)
            if str(task.get("task_type") or "").strip() != nextstrain_build_task_type:
                abort(404)
            dataset_path = resolve_nextstrain_dataset_sidecar(project_root, task, dataset_type)
        elif resolved_type == "uploaded":
            metadata = load_uploaded_auspice_metadata(project_root, dataset_id)
            ensure_can_view_uploaded_auspice(metadata)
            dataset_path = uploaded_auspice_dataset_path(project_root, dataset_id)
        else:
            abort(404)
        if not dataset_path.is_file():
            abort(404)
        response = send_file(str(dataset_path), mimetype="application/json", conditional=True)
        response.headers["Cache-Control"] = "no-store, no-cache, must-revalidate, max-age=0"
        response.headers["Pragma"] = "no-cache"
        response.headers["Expires"] = "0"
        return response

    @app.get("/getDataset")
    @app.get("/charon/getDataset")
    @login_required
    def local_auspice_dataset():
        return _local_auspice_dataset_response()
