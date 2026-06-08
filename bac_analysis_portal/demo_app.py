from __future__ import annotations

from pathlib import Path

from flask import Response, abort, current_app, jsonify, redirect, request, send_from_directory, url_for

from .app import create_app


HTML_INJECTION = """
<link rel="stylesheet" href="/static/demo.css">
<script>window.__BAC_PORTAL_DEMO__ = {"mode":"lite","name":"客户演示版"};</script>
<script defer src="/static/demo.js"></script>
"""

READ_ONLY_ALLOWED_POSTS = {
    "/login",
    "/logout",
    "/api/tasks",
    "/api/tasks/demo",
}

READ_ONLY_BLOCKED_PREFIXES = (
    "/api/host-database",
    "/api/pathogen-database",
    "/api/reference-database",
    "/api/knowledge-base",
)

READ_ONLY_BLOCKED_EXACT = {
    "/api/database/local-import",
    "/api/database/batch-import",
    "/api/database/batch-update",
    "/api/database/metadata-templates",
    "/api/database/import-template",
    "/api/export-table",
}


def _inject_demo_assets(response: Response) -> Response:
    content_type = str(response.headers.get("Content-Type") or "").lower()
    if response.status_code != 200 or "text/html" not in content_type:
        return response
    try:
        body = response.get_data(as_text=True)
    except Exception:
        return response
    if "/demo-static/demo.css" in body:
        return response
    if "</head>" in body:
        body = body.replace("</head>", f"{HTML_INJECTION}\n</head>", 1)
    elif "</body>" in body:
        body = body.replace("</body>", f"{HTML_INJECTION}\n</body>", 1)
    response.set_data(body)
    response.headers["Cache-Control"] = "no-store, max-age=0"
    return response


def _is_blocked_request() -> bool:
    path = request.path or ""
    if path.startswith("/demo-static/"):
        return False
    if request.method in {"POST", "PUT", "PATCH", "DELETE"}:
        if path in READ_ONLY_ALLOWED_POSTS:
            return False
        if path in READ_ONLY_BLOCKED_EXACT:
            return True
        if path.startswith("/api/tasks/") and path not in {"/api/tasks/demo"}:
            return True
        if path.startswith("/api/database/"):
            return True
    if any(path.startswith(prefix) for prefix in READ_ONLY_BLOCKED_PREFIXES):
        return True
    return False


def create_demo_app():
    app = create_app()
    app.config["ENABLE_KNOWLEDGE_BASE_TEST_PANEL"] = False

    @app.before_request
    def demo_read_only_guard():
        if not _is_blocked_request():
            return None
        if request.path.startswith("/api/"):
            return jsonify({"error": "当前为客户演示版，已关闭该操作。"}), 403
        return redirect(url_for("index", module="bacteria"))

    @app.after_request
    def demo_after_request(response: Response):
        return _inject_demo_assets(response)

    return app


if __name__ == "__main__":
    demo_app = create_demo_app()
    demo_app.run(host="127.0.0.1", port=5060, debug=False)
