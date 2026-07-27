from __future__ import annotations

from flask import abort, current_app, jsonify, request, session

from .app_services import get_app_services
from .identity import UserIdentity


def register_knowledge_routes(app) -> None:
    services = get_app_services(app)
    login_required = services.access.login_required
    knowledge = services.knowledge_service

    def identity() -> UserIdentity:
        return UserIdentity(
            username=str(session.get("username") or ""),
            role=str(session.get("role") or ""),
            group_name=str(session.get("group_name") or ""),
        )

    def no_cache(payload):
        response = jsonify(payload)
        response.headers["Cache-Control"] = "no-store, max-age=0"
        response.headers["Pragma"] = "no-cache"
        response.headers["Expires"] = "0"
        return response

    def ensure_panel_enabled() -> None:
        if not current_app.config.get("ENABLE_KNOWLEDGE_BASE_TEST_PANEL", False):
            abort(404)

    @app.get("/api/knowledge-base/summary")
    @login_required
    def knowledge_base_summary():
        ensure_panel_enabled()
        return no_cache(knowledge.summary())

    @app.get("/api/knowledge-base/bundle")
    @login_required
    def knowledge_base_bundle():
        ensure_panel_enabled()
        return no_cache(knowledge.bundle())

    @app.get("/api/report-templates/virus")
    @login_required
    def list_virus_report_templates():
        return no_cache({"items": knowledge.list_virus_templates()})

    @app.put("/api/report-templates/virus/<template_id>")
    @login_required
    def update_virus_report_template(template_id: str):
        return jsonify(knowledge.update_virus_template(template_id, request.get_json(force=True), identity=identity()))

    @app.delete("/api/report-templates/virus/<template_id>")
    @login_required
    def reset_virus_report_template(template_id: str):
        return jsonify(knowledge.reset_virus_template(template_id, identity=identity()))

    @app.post("/api/report-templates/virus/<template_id>/rollback")
    @login_required
    def rollback_virus_report_template(template_id: str):
        return jsonify(knowledge.rollback_virus_template(template_id, request.get_json(force=True), identity=identity()))
