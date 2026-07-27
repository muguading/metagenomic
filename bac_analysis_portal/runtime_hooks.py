from __future__ import annotations

import os
from pathlib import Path
from typing import Callable

from flask import jsonify, redirect, request, url_for

from .store import PortalStore
from .task_manager import AnalysisTaskManager


def register_runtime_hooks(
    app,
    *,
    project_root: Path,
    store: PortalStore,
    task_manager: AnalysisTaskManager,
    is_logged_in: Callable[[], bool],
) -> None:
    def sync_runtime_task_manager() -> None:
        raw_value = str(store.get_setting("workspace_root", "") or "").strip()
        try:
            runtime_root = Path(raw_value).expanduser().resolve() if raw_value else project_root
        except OSError:
            runtime_root = project_root
        if not runtime_root.is_dir():
            runtime_root = project_root
        runtime_task_root = Path(os.environ.get("BAC_ANALYSIS_TASK_ROOT", runtime_root / "analysis_tasks")).expanduser()
        runtime_task_root.mkdir(parents=True, exist_ok=True)
        task_manager.project_root = runtime_root
        task_manager.task_root = runtime_task_root

    @app.before_request
    def protect_routes():
        sync_runtime_task_manager()
        if request.endpoint in {"login", "login_post", "static"}:
            return None
        if not is_logged_in():
            if request.path.startswith("/api/"):
                return jsonify({"error": "Authentication required"}), 401
            return redirect(url_for("login"))
        return None
