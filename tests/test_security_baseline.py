from __future__ import annotations

from flask import Flask
from werkzeug.security import check_password_hash

from bac_analysis_portal.application import create_app
from bac_analysis_portal.security import LoginRateLimiter, configure_session_security, resolve_portal_mode, validate_security_config
from bac_analysis_portal.store import PortalStore


def test_portal_mode_defaults_to_production_and_testing_selects_test(monkeypatch) -> None:
    monkeypatch.delenv("PORTAL_MODE", raising=False)
    assert resolve_portal_mode() == "production"
    assert resolve_portal_mode({"TESTING": True}) == "test"
    assert resolve_portal_mode({"PORTAL_MODE": "demo"}) == "demo"
    monkeypatch.setenv("PORTAL_MODE", "development")
    assert resolve_portal_mode() == "development"


def test_production_requires_secret_and_strong_initial_password() -> None:
    app = Flask(__name__)
    app.config.update(PORTAL_MODE="production", SECRET_KEY="", INITIAL_ADMIN_PASSWORD="")
    try:
        validate_security_config(app)
    except RuntimeError as error:
        assert "PORTAL_SECRET_KEY" in str(error)
    else:
        raise AssertionError("production mode accepted a missing secret")

    app.config.update(SECRET_KEY="strong-secret", INITIAL_ADMIN_PASSWORD="admin123")
    try:
        validate_security_config(app)
    except RuntimeError as error:
        assert "demo password" in str(error)
    else:
        raise AssertionError("production mode accepted the demo password")


def test_session_security_defaults() -> None:
    app = Flask(__name__)
    app.config["PORTAL_COOKIE_SECURE"] = True
    configure_session_security(app)
    assert app.config["SESSION_COOKIE_HTTPONLY"] is True
    assert app.config["SESSION_COOKIE_SAMESITE"] == "Lax"
    assert app.config["SESSION_COOKIE_SECURE"] is True


def test_csrf_rejects_missing_token_and_accepts_valid_token() -> None:
    app = create_app({"TESTING": True, "SECRET_KEY": "test-secret"})
    assert "INITIAL_ADMIN_PASSWORD" not in app.config
    client = app.test_client()

    rejected = client.post("/login", json={"username": "admin", "password": "admin123"})
    assert rejected.status_code == 403

    client.get("/login")
    with client.session_transaction() as session:
        token = session["_csrf_token"]
    accepted = client.post(
        "/login",
        json={"username": "admin", "password": "admin123"},
        headers={"X-CSRF-Token": token},
    )
    assert accepted.status_code == 200


def test_login_rate_limiter_blocks_and_resets() -> None:
    limiter = LoginRateLimiter(max_attempts=2, window_seconds=10)
    limiter.record_failure("key", now=100)
    assert limiter.is_blocked("key", now=101) is False
    limiter.record_failure("key", now=102)
    assert limiter.is_blocked("key", now=103) is True
    assert limiter.is_blocked("key", now=113) is False
    limiter.record_failure("key", now=114)
    limiter.clear("key")
    assert limiter.is_blocked("key", now=115) is False


def test_production_store_bootstrap_and_weak_password_rotation(tmp_path) -> None:
    store = PortalStore(db_path=tmp_path / "portal.sqlite3", project_root=tmp_path)
    store.initialize(initial_admin_password="admin123")
    assert store.authenticate("admin", "admin123") is not None

    store.initialize(initial_admin_password="replacement-password", rotate_weak_admin=True)
    assert store.authenticate("admin", "admin123") is None
    assert store.authenticate("admin", "replacement-password") is not None
    assert check_password_hash(store.get_user_with_hash("admin")["password_hash"], "replacement-password")


def test_login_endpoint_rate_limits_repeated_failures() -> None:
    app = create_app({"TESTING": True, "SECRET_KEY": "test-secret"})
    client = app.test_client()
    client.get("/login")
    with client.session_transaction() as session:
        token = session["_csrf_token"]
    headers = {"X-CSRF-Token": token}

    for _ in range(5):
        response = client.post("/login", json={"username": "admin", "password": "wrong"}, headers=headers)
        assert response.status_code == 400
    blocked = client.post("/login", json={"username": "admin", "password": "wrong"}, headers=headers)
    assert blocked.status_code == 429


def test_unknown_errors_return_generic_message() -> None:
    app = create_app({"TESTING": True, "PROPAGATE_EXCEPTIONS": False, "SECRET_KEY": "test-secret"})

    @app.get("/raise-test-error")
    def raise_test_error():
        raise RuntimeError("sensitive internal detail")

    client = app.test_client()
    client.get("/login")
    with client.session_transaction() as session:
        token = session["_csrf_token"]
    login = client.post(
        "/login",
        json={"username": "admin", "password": "admin123"},
        headers={"X-CSRF-Token": token},
    )
    assert login.status_code == 200
    response = client.get("/raise-test-error")
    assert response.status_code == 500
    assert response.get_json() == {"error": "服务器内部错误"}
