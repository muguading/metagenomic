from __future__ import annotations

from flask import Response, jsonify, render_template, request, session, url_for

from .app_services import get_app_services
from .community_reports import _is_community_report_task
from .identity import UserIdentity
from .report_cache import REPORT_CACHE_VERSION
from .report_payload import _build_report_payload
from .result_page import _inject_result_preview_style, build_result_back_target


def register_task_lifecycle_routes(app) -> None:
    services = get_app_services(app)
    login_required = services.access.login_required
    lifecycle = services.task_lifecycle_service
    resolve_task_result_html = services.reports.resolve_task_result_html
    task_report_availability = services.reports.task_report_availability

    def identity() -> UserIdentity:
        return UserIdentity(
            username=str(session.get("username") or ""),
            role=str(session.get("role") or ""),
            group_name=str(session.get("group_name") or ""),
        )

    @app.get("/api/tasks/<task_id>")
    @login_required
    def get_task(task_id: str):
        log_lines = request.args.get("log_lines", default=120, type=int)
        task = lifecycle.get_visible(task_id, log_lines=log_lines)
        report_availability = task_report_availability(task)
        result_html = report_availability.get("path")
        task["result_exists"] = bool(report_availability.get("available"))
        task["result_url"] = (
            url_for("task_result_view", task_id=task_id)
            if report_availability.get("mode") == "static"
            else (url_for("task_result_page", task_id=task_id) if task["result_exists"] else "")
        )
        task["result_name"] = result_html.name if result_html else ""
        return jsonify(task)

    @app.delete("/api/tasks/<task_id>")
    @login_required
    def delete_task(task_id: str):
        lifecycle.delete(task_id)
        return jsonify({"status": "deleted", "id": task_id})

    @app.post("/api/tasks/<task_id>/pause")
    @login_required
    def pause_task(task_id: str):
        return jsonify(lifecycle.pause(task_id))

    @app.post("/api/tasks/<task_id>/resume")
    @login_required
    def resume_task(task_id: str):
        return jsonify(lifecycle.resume(task_id))

    @app.post("/api/tasks/<task_id>/stop")
    @login_required
    def stop_task(task_id: str):
        return jsonify(lifecycle.stop(task_id))

    @app.post("/api/tasks/<task_id>/open-output")
    @login_required
    def open_task_output(task_id: str):
        target = lifecycle.open_output(task_id)
        return jsonify({"status": "opened", "path": str(target)})

    @app.post("/api/tasks/<task_id>/database-import")
    @login_required
    def import_task_to_database(task_id: str):
        return jsonify(lifecycle.import_to_database(task_id, request.get_json(silent=True) or {}, identity=identity()))

    @app.get("/api/tasks/<task_id>/database-import-preview")
    @login_required
    def preview_task_database_import(task_id: str):
        return jsonify(lifecycle.preview_database_import(task_id))

    @app.post("/api/tasks/<task_id>/closure-actions/<action_id>")
    @login_required
    def confirm_task_closure_action(task_id: str, action_id: str):
        return jsonify(lifecycle.confirm_closure_action(task_id, action_id, request.get_json(silent=True) or {}, identity=identity()))

    @app.post("/api/tasks/<task_id>/report-exports")
    @login_required
    def record_task_report_export(task_id: str):
        return jsonify(lifecycle.record_report_export(task_id, request.get_json(silent=True) or {}))

    @app.post("/api/tasks/<task_id>/failure-disposition")
    @login_required
    def record_task_failure_disposition(task_id: str):
        return jsonify(lifecycle.record_failure_disposition(task_id, request.get_json(silent=True) or {}, identity=identity()))

    @app.get("/api/tasks/<task_id>/result-view")
    @login_required
    def task_result_view(task_id: str):
        task = lifecycle.get_visible(task_id)
        result_html = resolve_task_result_html(task)
        if result_html is None or not result_html.is_file():
            if _is_community_report_task(task):
                payload = _build_report_payload(task)
                return render_template("community_result_view.html", report=payload)
            raise KeyError(f"Result not found for task: {task_id}")
        html = result_html.read_text(encoding="utf-8", errors="ignore")
        html = _inject_result_preview_style(html, task_id=task_id)
        return Response(html, mimetype="text/html")

    @app.get("/tasks/<task_id>/result-page")
    @login_required
    def task_result_page(task_id: str):
        task = lifecycle.get_visible(task_id)
        selected_sample = str(request.args.get("sample", "") or "").strip()
        back_target = build_result_back_target(
            task_id,
            return_to=str(request.args.get("return_to", "") or ""),
            sample_key=str(request.args.get("sample_key", "") or ""),
            database_section=str(request.args.get("database_section", "") or ""),
        )
        return render_template(
            "result_report.html",
            task=task,
            selected_sample=selected_sample,
            back_href=back_target["href"],
            back_label=back_target["label"],
            report_static_version=REPORT_CACHE_VERSION,
        )
