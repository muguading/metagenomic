from __future__ import annotations

from flask import jsonify, request, session

from .app_services import get_app_services
from .identity import UserIdentity


def register_admin_monitor_routes(app) -> None:
    services = get_app_services(app)
    admin_required = services.access.admin_required
    monitor = services.admin_monitor_service

    @app.post("/api/admin/watch-tasks")
    @admin_required
    def create_admin_watch_tasks():
        identity = UserIdentity(
            username=str(session.get("username") or "admin"),
            role=str(session.get("role") or ""),
            group_name=str(session.get("group_name") or ""),
        )
        return jsonify({"items": monitor.create_watch_tasks(request.get_json(force=True), identity=identity)}), 201
