from __future__ import annotations

from flask import jsonify, redirect, render_template, request, session, url_for

from .app_services import get_app_services
from .task_manager import ASM_METHOD_OPTIONS
from .security import login_rate_limit_key


def register_portal_routes(app) -> None:
    services = get_app_services(app)
    project_root = services.project_root
    login_required = services.access.login_required
    is_logged_in = services.access.is_logged_in
    portal = services.portal_service

    @app.get("/login")
    def login():
        if is_logged_in():
            return redirect(url_for("index", module=portal.default_module(str(session.get("username") or ""))))
        return render_template("login.html")

    @app.post("/login")
    def login_post():
        payload = request.get_json(force=True) if request.is_json else request.form
        username = str(payload.get("username", "")).strip()
        password = str(payload.get("password", ""))
        limiter = app.extensions["portal_login_rate_limiter"]
        rate_limit_key = login_rate_limit_key(username)
        if limiter.is_blocked(rate_limit_key):
            return jsonify({"error": "登录失败次数过多，请稍后再试"}), 429
        try:
            user = portal.authenticate(username, password)
        except Exception:
            limiter.record_failure(rate_limit_key)
            raise
        limiter.clear(rate_limit_key)
        session.clear()
        session["username"] = user["username"]
        session["role"] = user["role"]
        session["group_name"] = user.get("group_name", "")
        session.permanent = True
        return jsonify({"status": "ok", "user": user})

    @app.post("/logout")
    @login_required
    def logout():
        session.clear()
        return jsonify({"status": "logged_out"})

    @app.get("/")
    @login_required
    def launcher():
        return redirect(url_for("index", module=portal.default_module(str(session.get("username") or ""))))

    @app.get("/workstation")
    @login_required
    def index():
        username = str(session.get("username") or "")
        workstation, is_allowed = portal.resolve_workstation(request.args.get("module"), username)
        if not is_allowed:
            return redirect(url_for("index", module=portal.default_module(username)))
        return render_template(
            "index.html",
            project_root=str(project_root),
            script_path=workstation["pipeline_script"],
            workstation_key=workstation["key"],
            workstation_label=workstation["label"],
            workstation_subtitle=workstation["subtitle"],
            knowledge_base_enabled=bool(app.config.get("ENABLE_KNOWLEDGE_BASE_TEST_PANEL", False)),
        )

    @app.get("/api/health")
    @login_required
    def health():
        return jsonify(portal.health())

    @app.get("/api/session")
    @login_required
    def current_session():
        return jsonify(portal.get_user(session["username"]))

    @app.get("/api/tasks")
    @login_required
    def list_tasks():
        return jsonify({"items": portal.list_visible_tasks()})

    @app.get("/api/asm-options")
    @login_required
    def asm_options():
        return jsonify({"items": ASM_METHOD_OPTIONS})

    @app.get("/api/server-status")
    @login_required
    def server_status():
        return jsonify(portal.server_status())
