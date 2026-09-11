from __future__ import annotations

import os
import secrets
import threading
import time
from collections import defaultdict, deque
from dataclasses import dataclass, field
from datetime import timedelta
from typing import Any, Mapping

from flask import Flask, jsonify, request, session


PORTAL_MODES = {"production", "development", "demo", "test"}
UNSAFE_METHODS = {"POST", "PUT", "PATCH", "DELETE"}


def resolve_portal_mode(config: Mapping[str, Any] | None = None) -> str:
    configured = str((config or {}).get("PORTAL_MODE") or "").strip().lower()
    if not configured and bool((config or {}).get("TESTING")):
        configured = "test"
    mode = configured or str(os.environ.get("PORTAL_MODE") or "production").strip().lower()
    if mode not in PORTAL_MODES:
        raise RuntimeError(f"PORTAL_MODE must be one of: {', '.join(sorted(PORTAL_MODES))}")
    return mode


def validate_security_config(app: Flask) -> None:
    mode = str(app.config["PORTAL_MODE"])
    if mode != "production":
        return
    if not str(app.config.get("SECRET_KEY") or "").strip():
        raise RuntimeError("PORTAL_SECRET_KEY is required when PORTAL_MODE=production")
    password = str(app.config.get("INITIAL_ADMIN_PASSWORD") or "")
    if not password:
        raise RuntimeError("PORTAL_INITIAL_ADMIN_PASSWORD is required when PORTAL_MODE=production")
    if password == "admin123":
        raise RuntimeError("PORTAL_INITIAL_ADMIN_PASSWORD must not use the demo password")


def csrf_token() -> str:
    token = str(session.get("_csrf_token") or "")
    if not token:
        token = secrets.token_urlsafe(32)
        session["_csrf_token"] = token
    return token


@dataclass
class LoginRateLimiter:
    max_attempts: int = 5
    window_seconds: int = 900
    _failures: dict[str, deque[float]] = field(default_factory=lambda: defaultdict(deque))
    _lock: threading.Lock = field(default_factory=threading.Lock)

    def _active(self, key: str, now: float) -> deque[float]:
        failures = self._failures[key]
        threshold = now - self.window_seconds
        while failures and failures[0] <= threshold:
            failures.popleft()
        return failures

    def is_blocked(self, key: str, *, now: float | None = None) -> bool:
        with self._lock:
            return len(self._active(key, time.time() if now is None else now)) >= self.max_attempts

    def record_failure(self, key: str, *, now: float | None = None) -> None:
        timestamp = time.time() if now is None else now
        with self._lock:
            self._active(key, timestamp).append(timestamp)

    def clear(self, key: str) -> None:
        with self._lock:
            self._failures.pop(key, None)


def login_rate_limit_key(username: str) -> str:
    ip_address = str(request.headers.get("X-Forwarded-For", request.remote_addr or "")).split(",")[0].strip()
    return f"{ip_address}:{str(username or '').strip().lower()}"


def register_security_hooks(app: Flask) -> None:
    app.extensions["portal_login_rate_limiter"] = LoginRateLimiter(
        max_attempts=int(app.config.get("LOGIN_RATE_LIMIT_ATTEMPTS", 5)),
        window_seconds=int(app.config.get("LOGIN_RATE_LIMIT_WINDOW_SECONDS", 900)),
    )
    app.jinja_env.globals["csrf_token"] = csrf_token

    @app.before_request
    def enforce_csrf():
        if request.method not in UNSAFE_METHODS:
            csrf_token()
            return None
        expected = str(session.get("_csrf_token") or "")
        supplied = str(request.headers.get("X-CSRF-Token") or request.form.get("_csrf_token") or "")
        if expected and secrets.compare_digest(expected, supplied):
            return None
        return jsonify({"error": "CSRF token missing or invalid"}), 403


def configure_session_security(app: Flask) -> None:
    app.config.update(
        SESSION_COOKIE_HTTPONLY=True,
        SESSION_COOKIE_SAMESITE="Lax",
        SESSION_COOKIE_SECURE=bool(app.config.get("PORTAL_COOKIE_SECURE")),
        PERMANENT_SESSION_LIFETIME=timedelta(hours=8),
    )
