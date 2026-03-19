from __future__ import annotations

from pathlib import Path

from flask import Flask, jsonify, redirect, render_template, request, session, url_for
from functools import wraps

from .task_manager import ASM_METHOD_OPTIONS, AnalysisTaskManager, ValidationError
from .store import PortalStore


def create_app() -> Flask:
    package_root = Path(__file__).resolve().parent
    project_root = package_root.parent
    app = Flask(
        __name__,
        template_folder=str(package_root / "templates"),
        static_folder=str(package_root / "static"),
    )
    app.secret_key = "bac-analysis-portal-dev-key"
    task_manager = AnalysisTaskManager.from_project_root(project_root)
    store = PortalStore.from_project_root(project_root)

    @app.before_request
    def protect_routes():
        if request.endpoint in {"login", "login_post", "static"}:
            return None
        if not _is_logged_in():
            if request.path.startswith("/api/"):
                return jsonify({"error": "Authentication required"}), 401
            return redirect(url_for("login"))
        return None

    @app.get("/login")
    def login():
        if _is_logged_in():
            return redirect(url_for("index"))
        return render_template("login.html")

    @app.post("/login")
    def login_post():
        payload = request.get_json(force=True) if request.is_json else request.form
        user = store.authenticate(str(payload.get("username", "")).strip(), str(payload.get("password", "")))
        if user is None:
            raise ValidationError("用户名或密码错误")
        session["username"] = user["username"]
        session["role"] = user["role"]
        return jsonify({"status": "ok", "user": user})

    @app.post("/logout")
    @login_required
    def logout():
        session.clear()
        return jsonify({"status": "logged_out"})

    @app.get("/")
    @login_required
    def index():
        return render_template(
            "index.html",
            project_root=str(project_root),
            script_path=store.get_setting("pipeline_script", "Bac_assemble_260112_newformat.py"),
        )

    @app.get("/api/health")
    @login_required
    def health():
        return jsonify(
            {
                "status": "ok",
                "workspace_root": store.get_setting("workspace_root", str(project_root)),
                "script_path": store.get_setting("pipeline_script", "Bac_assemble_260112_newformat.py"),
            }
        )

    @app.get("/api/session")
    @login_required
    def current_session():
        return jsonify(store.get_user(session["username"]))

    @app.get("/api/tasks")
    @login_required
    def list_tasks():
        owner = None if session.get("role") == "admin" else session["username"]
        return jsonify({"items": task_manager.list_tasks(owner=owner)})

    @app.get("/api/asm-options")
    @login_required
    def asm_options():
        return jsonify({"items": ASM_METHOD_OPTIONS})

    @app.get("/api/filesystem")
    @login_required
    def filesystem():
        path_arg = request.args.get("path", default="", type=str)
        selector = request.args.get("selector", default="input", type=str)
        base_root = project_root.resolve()
        current_path = _resolve_browser_path(base_root, path_arg)
        if not current_path.is_dir():
            raise ValidationError(f"只能浏览目录: {current_path}")

        items = []
        for child in sorted(current_path.iterdir(), key=lambda item: (not item.is_dir(), item.name.lower())):
            if child.name.startswith(".") and child.name not in {".", ".."}:
                continue
            item_type = "directory" if child.is_dir() else "file"
            if selector == "output" and item_type != "directory":
                continue
            items.append(
                {
                    "name": child.name,
                    "type": item_type,
                    "path": _to_browser_path(base_root, child),
                }
            )

        parent_relative = ""
        if current_path != base_root:
            parent_relative = _to_browser_path(base_root, current_path.parent)

        return jsonify(
            {
                "root": str(base_root),
                "current_path": str(current_path),
                "relative_path": _to_browser_path(base_root, current_path),
                "parent_relative_path": parent_relative,
                "selector": selector,
                "items": items,
            }
        )

    @app.post("/api/filesystem/mkdir")
    @login_required
    def filesystem_mkdir():
        payload = request.get_json(force=True)
        base_root = project_root.resolve()
        parent = _resolve_browser_path(base_root, payload.get("path", ""))
        if not parent.is_dir():
            raise ValidationError(f"只能在目录下新建文件夹: {parent}")
        name = str(payload.get("name", "")).strip()
        if not name:
            raise ValidationError("文件夹名称不能为空")
        if any(token in name for token in ["/", "\\"]):
            raise ValidationError("文件夹名称不能包含路径分隔符")
        target = (parent / name).resolve()
        _assert_within_root(base_root, target)
        target.mkdir(exist_ok=False)
        return jsonify({"status": "ok", "path": _to_browser_path(base_root, target)})

    @app.post("/api/filesystem/rename")
    @login_required
    def filesystem_rename():
        payload = request.get_json(force=True)
        base_root = project_root.resolve()
        source = _resolve_browser_path(base_root, payload.get("path", ""))
        if not source.exists():
            raise ValidationError(f"待重命名目标不存在: {source}")
        new_name = str(payload.get("name", "")).strip()
        if not new_name:
            raise ValidationError("新名称不能为空")
        if any(token in new_name for token in ["/", "\\"]):
            raise ValidationError("新名称不能包含路径分隔符")
        target = (source.parent / new_name).resolve()
        _assert_within_root(base_root, target)
        source.rename(target)
        return jsonify({"status": "ok", "path": _to_browser_path(base_root, target)})

    @app.get("/api/tasks/<task_id>")
    @login_required
    def get_task(task_id: str):
        log_lines = request.args.get("log_lines", default=120, type=int)
        owner = None if session.get("role") == "admin" else session["username"]
        return jsonify(task_manager.get_task(task_id, log_lines=log_lines, owner=owner))

    @app.post("/api/tasks")
    @login_required
    def create_task():
        payload = request.get_json(force=True)
        workspace_root = Path(store.get_setting("workspace_root", str(project_root))).expanduser().resolve()
        created = task_manager.create_task(
            payload,
            owner=session["username"],
            pipeline_script=_resolve_pipeline_script(
                workspace_root,
                store.get_setting("pipeline_script", "Bac_assemble_260112_newformat.py"),
            ),
        )
        return jsonify(created), 201

    @app.get("/api/admin/settings")
    @admin_required
    def admin_settings():
        return jsonify(
            {
                "workspace_root": store.get_setting("workspace_root", str(project_root)),
                "pipeline_script": store.get_setting("pipeline_script", "Bac_assemble_260112_newformat.py"),
            }
        )

    @app.get("/api/admin/filesystem")
    @admin_required
    def admin_filesystem():
        selector = request.args.get("selector", default="workspace_root", type=str)
        root_arg = request.args.get("root", default="", type=str)
        path_arg = request.args.get("path", default="", type=str)
        browse_root = _resolve_admin_root(root_arg)
        current_path = _resolve_admin_path(browse_root, path_arg)
        if not current_path.is_dir():
            raise ValidationError(f"只能浏览目录: {current_path}")

        items = []
        try:
            children = sorted(current_path.iterdir(), key=lambda item: (not item.is_dir(), item.name.lower()))
        except PermissionError as exc:
            raise ValidationError(f"没有访问目录权限: {current_path}") from exc

        for child in children:
            if child.name.startswith(".") and child.name not in {".", ".."}:
                continue
            item_type = "directory" if child.is_dir() else "file"
            if selector == "workspace_root" and item_type != "directory":
                continue
            items.append({"name": child.name, "type": item_type, "path": str(child.resolve())})

        parent_path = ""
        if current_path != browse_root:
            parent_path = str(current_path.parent.resolve())

        return jsonify(
            {
                "root": str(browse_root),
                "current_path": str(current_path),
                "parent_path": parent_path,
                "selector": selector,
                "items": items,
            }
        )

    @app.put("/api/admin/settings")
    @admin_required
    def update_admin_settings():
        payload = request.get_json(force=True)
        workspace_root = str(payload.get("workspace_root", "")).strip()
        script_path = str(payload.get("pipeline_script", "")).strip()
        if not workspace_root:
            raise ValidationError("部署基准目录不能为空")
        if not script_path:
            raise ValidationError("脚本路径不能为空")
        workspace_root_path = Path(workspace_root).expanduser().resolve()
        if not workspace_root_path.is_dir():
            raise ValidationError(f"部署基准目录不存在: {workspace_root_path}")
        candidate = Path(_resolve_pipeline_script(workspace_root_path, script_path))
        if not candidate.is_file():
            raise ValidationError(f"脚本路径不存在: {candidate}")
        store.set_setting("workspace_root", str(workspace_root_path))
        store.set_setting("pipeline_script", str(candidate.relative_to(workspace_root_path)))
        return jsonify(
            {
                "workspace_root": str(workspace_root_path),
                "pipeline_script": str(candidate.relative_to(workspace_root_path)),
            }
        )

    @app.get("/api/admin/users")
    @admin_required
    def list_users():
        return jsonify({"items": store.list_users()})

    @app.post("/api/admin/users")
    @admin_required
    def create_user():
        payload = request.get_json(force=True)
        username = str(payload.get("username", "")).strip()
        password = str(payload.get("password", ""))
        role = str(payload.get("role", "user")).strip()
        if not username or not password:
            raise ValidationError("用户名和密码不能为空")
        created = store.create_user(
            username=username,
            password=password,
            role=role,
            display_name=str(payload.get("display_name", "")).strip(),
        )
        return jsonify(created), 201

    @app.put("/api/admin/users/<username>")
    @admin_required
    def update_user(username: str):
        payload = request.get_json(force=True)
        updated = store.update_user(
            username,
            role=str(payload.get("role", "")).strip() or None,
            display_name=str(payload.get("display_name", "")).strip() if "display_name" in payload else None,
            new_password=str(payload.get("password", "")).strip() or None,
        )
        return jsonify(updated)

    @app.errorhandler(ValidationError)
    def handle_validation(error: ValidationError):
        return jsonify({"error": str(error)}), 400

    @app.errorhandler(KeyError)
    def handle_missing(error: KeyError):
        return jsonify({"error": str(error)}), 404

    @app.errorhandler(Exception)
    def handle_unknown(error: Exception):
        app.logger.exception("Bac analysis portal error")
        return jsonify({"error": str(error)}), 500

    return app


def login_required(view_func):
    @wraps(view_func)
    def wrapped(*args, **kwargs):
        if not _is_logged_in():
            return jsonify({"error": "Authentication required"}), 401
        return view_func(*args, **kwargs)

    return wrapped


def admin_required(view_func):
    @wraps(view_func)
    def wrapped(*args, **kwargs):
        if not _is_logged_in():
            return jsonify({"error": "Authentication required"}), 401
        if session.get("role") != "admin":
            return jsonify({"error": "Administrator permission required"}), 403
        return view_func(*args, **kwargs)

    return wrapped


def _is_logged_in() -> bool:
    return bool(session.get("username") and session.get("role"))


def _resolve_browser_path(base_root: Path, path_arg: str) -> Path:
    raw = str(path_arg or "").strip()
    if not raw:
        return base_root
    candidate = Path(raw).expanduser()
    candidate = candidate.resolve() if candidate.is_absolute() else (base_root / candidate).resolve()
    _assert_within_root(base_root, candidate)
    return candidate


def _to_browser_path(base_root: Path, path: Path) -> str:
    if path == base_root:
        return ""
    return str(path.relative_to(base_root))


def _assert_within_root(base_root: Path, candidate: Path) -> None:
    try:
        candidate.relative_to(base_root)
    except ValueError as exc:
        raise ValidationError(f"路径超出允许范围: {candidate}") from exc


def _resolve_pipeline_script(project_root: Path, script_path: str) -> str:
    raw = str(script_path or "").strip()
    if not raw:
        raise ValidationError("脚本路径不能为空")
    candidate = Path(raw).expanduser()
    candidate = candidate.resolve() if candidate.is_absolute() else (project_root / candidate).resolve()
    _assert_within_root(project_root.resolve(), candidate)
    return str(candidate)

def _resolve_admin_root(root_arg: str) -> Path:
    raw = str(root_arg or "").strip()
    if not raw:
        return Path.home().resolve()
    candidate = Path(raw).expanduser().resolve()
    if not candidate.is_dir():
        raise ValidationError(f"目录不存在: {candidate}")
    return candidate


def _resolve_admin_path(root: Path, path_arg: str) -> Path:
    raw = str(path_arg or "").strip()
    if not raw:
        return root
    candidate = Path(raw).expanduser()
    return candidate.resolve() if candidate.is_absolute() else (root / candidate).resolve()


if __name__ == "__main__":
    app = create_app()
    app.run(host="0.0.0.0", port=5055, debug=True)
