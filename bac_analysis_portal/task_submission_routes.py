from __future__ import annotations

from flask import jsonify, request, session

from .app_services import get_app_services
from .identity import UserIdentity


def register_task_submission_routes(app) -> None:
    services = get_app_services(app)
    login_required = services.access.login_required
    submissions = services.task_submission_service

    def identity() -> UserIdentity:
        return UserIdentity(
            username=str(session.get("username") or ""),
            role=str(session.get("role") or ""),
            group_name=str(session.get("group_name") or ""),
        )

    @app.post("/api/tasks")
    @login_required
    def create_task():
        return jsonify(submissions.create(request.get_json(force=True), identity=identity())), 201

    @app.post("/api/tasks/<task_id>/rerun")
    @login_required
    def rerun_task(task_id: str):
        return jsonify(submissions.rerun(task_id, request.get_json(force=True), identity=identity())), 200

    @app.post("/api/tasks/demo")
    @login_required
    def create_demo_task():
        return jsonify(submissions.create_demo(request.get_json(silent=True) or {}, identity=identity())), 201
