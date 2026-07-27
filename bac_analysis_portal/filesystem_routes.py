from __future__ import annotations

from flask import jsonify, request

from .admin_runtime import _resolve_home_and_desktop_paths
from .app_services import get_app_services
from .filesystem_helpers import (
    _assert_within_root,
    _is_within_root,
    _resolve_browser_path,
    _to_browser_path,
)
from .task_manager import ValidationError


def register_filesystem_routes(app) -> None:
    services = get_app_services(app)
    project_root = services.project_root
    login_required = services.access.login_required

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
            try:
                modified_at = child.stat().st_mtime
            except OSError:
                modified_at = None
            items.append(
                {
                    "name": child.name,
                    "type": item_type,
                    "path": _to_browser_path(base_root, child),
                    "modified_at": modified_at,
                }
            )

        parent_relative = ""
        if current_path != current_path.parent:
            parent_relative = _to_browser_path(base_root, current_path.parent)

        home_path, desktop_path = _resolve_home_and_desktop_paths()
        return jsonify(
            {
                "root": str(base_root),
                "current_path": str(current_path),
                "relative_path": _to_browser_path(base_root, current_path),
                "parent_relative_path": parent_relative,
                "home_path": str(home_path),
                "desktop_path": str(desktop_path),
                "selector": selector,
                "within_root": _is_within_root(base_root, current_path),
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
        target = services.filesystem.create_directory(parent, payload.get("name", ""), allowed_root=base_root)
        return jsonify({"status": "ok", "path": _to_browser_path(base_root, target)})

    @app.post("/api/filesystem/rename")
    @login_required
    def filesystem_rename():
        payload = request.get_json(force=True)
        base_root = project_root.resolve()
        source = _resolve_browser_path(base_root, payload.get("path", ""))
        target = services.filesystem.rename_path(source, payload.get("name", ""), allowed_root=base_root)
        return jsonify({"status": "ok", "path": _to_browser_path(base_root, target)})
