from __future__ import annotations

from flask import jsonify, request

from .app_services import get_app_services


def register_admin_settings_routes(app) -> None:
    services = get_app_services(app)
    admin_required = services.access.admin_required
    settings = services.admin_settings_service

    @app.get("/api/admin/settings")
    @admin_required
    def admin_settings():
        return jsonify(settings.get_settings())

    @app.get("/api/admin/conda-envs")
    @admin_required
    def admin_conda_envs():
        return jsonify(settings.list_conda_envs(str(request.args.get("root", "") or "").strip()))

    @app.get("/api/admin/pathosource-trigger-rules")
    @admin_required
    def admin_pathosource_trigger_rules():
        return jsonify(settings.load_pathosource_trigger_rules())

    @app.put("/api/admin/pathosource-trigger-rules")
    @admin_required
    def update_admin_pathosource_trigger_rules():
        return jsonify(settings.save_pathosource_trigger_rules(request.get_json(force=True)))

    @app.get("/api/admin/filesystem")
    @admin_required
    def admin_filesystem():
        return jsonify(
            settings.browse_filesystem(
                selector=request.args.get("selector", default="workspace_root", type=str),
                root=request.args.get("root", default="", type=str),
                path=request.args.get("path", default="", type=str),
            )
        )

    @app.post("/api/admin/filesystem/mkdir")
    @admin_required
    def admin_filesystem_mkdir():
        payload = request.get_json(force=True)
        target = settings.create_folder(current_path=payload.get("path", ""), name=payload.get("name", ""))
        return jsonify({"status": "created", "path": str(target.resolve())}), 201

    @app.put("/api/admin/settings")
    @admin_required
    def update_admin_settings():
        return jsonify(settings.update_settings(request.get_json(force=True)))

    @app.get("/api/admin/update/sources")
    @admin_required
    def admin_update_sources():
        return jsonify({"items": settings.update_sources()})

    @app.post("/api/admin/update/online")
    @admin_required
    def admin_online_update():
        return jsonify(settings.run_online_update())

    @app.post("/api/admin/update/offline")
    @admin_required
    def admin_offline_update():
        return jsonify(settings.run_offline_update(request.get_json(force=True).get("source_path", "")))

    @app.get("/api/admin/nextclade-datasets")
    @admin_required
    def admin_nextclade_datasets():
        return jsonify(settings.check_nextclade_datasets())

    @app.post("/api/admin/nextclade-datasets/update")
    @admin_required
    def admin_update_nextclade_datasets():
        payload = request.get_json(silent=True) or {}
        return jsonify(settings.update_nextclade_datasets(payload.get("directories")))
