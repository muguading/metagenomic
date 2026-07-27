from __future__ import annotations

from flask import jsonify, request, session

from .app_services import get_app_services
def register_admin_user_routes(app) -> None:
    services = get_app_services(app)
    admin_required = services.access.admin_required
    users = services.admin_user_service

    @app.get("/api/admin/users")
    @admin_required
    def list_users():
        return jsonify({"items": users.list_users()})

    @app.post("/api/admin/users")
    @admin_required
    def create_user():
        return jsonify(users.create(request.get_json(force=True))), 201

    @app.put("/api/admin/users/<username>")
    @admin_required
    def update_user(username: str):
        updated = users.update(username, request.get_json(force=True))
        if session.get("username") == username:
            session["username"] = updated["username"]
            session["role"] = updated["role"]
            session["group_name"] = updated.get("group_name", "")
        return jsonify(updated)

    @app.delete("/api/admin/users/<username>")
    @admin_required
    def delete_user(username: str):
        users.delete(username, current_username=str(session.get("username") or ""))
        return jsonify({"status": "deleted", "username": username})
