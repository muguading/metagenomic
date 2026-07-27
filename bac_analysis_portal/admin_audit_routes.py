from __future__ import annotations

from flask import jsonify, request

from .app_services import get_app_services


def register_admin_audit_routes(app) -> None:
    services = get_app_services(app)
    admin_required = services.access.admin_required
    audit = services.admin_audit_service

    @app.get("/api/admin/audit-logs")
    @admin_required
    def admin_audit_logs():
        return jsonify(
            audit.query(
                username=str(request.args.get("username", "") or "").strip(),
                module=str(request.args.get("module", "") or "").strip(),
                action=str(request.args.get("action", "") or "").strip(),
                outcome=str(request.args.get("outcome", "") or "").strip(),
                search=str(request.args.get("search", "") or "").strip(),
                limit=request.args.get("limit", default=500, type=int),
            )
        )
