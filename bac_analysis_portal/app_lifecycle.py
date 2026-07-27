from __future__ import annotations

import os
import threading
import time
import uuid

from flask import Flask, jsonify
from werkzeug.exceptions import HTTPException

from .app_services import get_app_services
from .application_errors import AuthorizationError
from .task_manager import ValidationError


def register_error_handlers(app: Flask) -> None:
    @app.errorhandler(ValidationError)
    def handle_validation(error: ValidationError):
        return jsonify({"error": str(error)}), 400

    @app.errorhandler(ValueError)
    def handle_value_error(error: ValueError):
        return jsonify({"error": str(error)}), 400

    @app.errorhandler(KeyError)
    def handle_missing(error: KeyError):
        return jsonify({"error": str(error)}), 404

    @app.errorhandler(AuthorizationError)
    def handle_authorization(error: AuthorizationError):
        return jsonify({"error": str(error)}), 403

    @app.errorhandler(Exception)
    def handle_unknown(error: Exception):
        if isinstance(error, HTTPException):
            return error
        app.logger.exception("Bac analysis portal error")
        return jsonify({"error": "服务器内部错误"}), 500


def start_queue_watch_daemon(app: Flask) -> None:
    if app.testing:
        return
    if app.debug and os.environ.get("WERKZEUG_RUN_MAIN") != "true":
        return

    services = get_app_services(app)
    owner_id = f"{os.getpid()}-{uuid.uuid4().hex}"

    def _worker() -> None:
        while True:
            try:
                if services.queue_maintenance_service.acquire_lease(owner_id, ttl_seconds=180):
                    services.queue_maintenance_service.maintain_once()
            except Exception:
                app.logger.exception("Queue watch daemon error")
            time.sleep(60)

    threading.Thread(target=_worker, name="task-queue-watch-daemon", daemon=True).start()
