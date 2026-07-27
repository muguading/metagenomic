from __future__ import annotations

import pytest

from conftest import create_portal_user, login_as, write_task_record


pytestmark = pytest.mark.release


def test_release_unauthenticated_api_requests_are_rejected(release_app_client) -> None:
    _app, client = release_app_client

    response = client.get("/api/tasks")

    assert response.status_code == 401
    assert response.get_json()["error"] == "Authentication required"


def test_release_regular_user_cannot_view_other_group_task_but_admin_can(release_app_client) -> None:
    app, client = release_app_client
    write_task_record(app, "release-private-task", status="SUCCEEDED", owner="admin", owner_group="admin-group")
    create_portal_user(app, username="analyst", password="analyst-pass", role="user", group_name="qa")

    user_headers = login_as(client, "analyst", "analyst-pass")
    forbidden = client.get("/api/tasks/release-private-task", headers=user_headers)
    client.post("/logout", headers=user_headers)

    admin_headers = login_as(client)
    allowed = client.get("/api/tasks/release-private-task", headers=admin_headers)

    assert forbidden.status_code == 404
    assert allowed.status_code == 200
    assert allowed.get_json()["id"] == "release-private-task"


def test_release_admin_endpoints_reject_regular_users(release_app_client) -> None:
    app, client = release_app_client
    create_portal_user(app, username="analyst", password="analyst-pass", role="user", group_name="qa")
    headers = login_as(client, "analyst", "analyst-pass")

    response = client.get("/api/host-database/records", headers=headers)
    admin_only = client.delete("/api/host-database/remote-import-jobs/job-1", headers=headers)

    assert response.status_code == 200
    assert admin_only.status_code == 403
    assert admin_only.get_json()["error"] == "Administrator permission required"


def test_release_csrf_rejects_missing_token_on_login(release_app_client) -> None:
    _app, client = release_app_client

    response = client.post("/login", json={"username": "admin", "password": "admin123"})

    assert response.status_code == 403


def test_release_report_asset_rejects_disallowed_extensions_and_path_traversal(release_app_client, tmp_path) -> None:
    app, client = release_app_client
    headers = login_as(client)
    report_dir = tmp_path / "report"
    report_dir.mkdir()
    (report_dir / "summary.tsv").write_text("sample\treads\nS1\t10\n", encoding="utf-8")
    (report_dir / "secret.py").write_text("print('no')\n", encoding="utf-8")
    write_task_record(app, "release-assets", status="SUCCEEDED", output_dir=report_dir)

    disallowed = client.get("/api/tasks/release-assets/report-asset/secret.py", headers=headers)
    traversal = client.get("/api/tasks/release-assets/report-asset/..%2Fsecret.py", headers=headers)

    assert disallowed.status_code == 403
    assert traversal.status_code in {403, 404}
