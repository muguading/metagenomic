from __future__ import annotations

import threading

import pytest

pytest.importorskip("playwright.sync_api")


pytestmark = pytest.mark.release


@pytest.fixture
def release_live_server(release_app_client):
    app, _client = release_app_client
    werkzeug = pytest.importorskip("werkzeug.serving")
    server = werkzeug.make_server("127.0.0.1", 0, app)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield f"http://127.0.0.1:{server.server_port}"
    finally:
        server.shutdown()
        thread.join(timeout=5)


def test_release_browser_can_login_open_queue_and_refresh_result_page(page, release_live_server) -> None:
    page.goto(f"{release_live_server}/login")
    page.fill("#login-username", "admin")
    page.fill("#login-password", "admin123")
    page.click("#login-form button[type='submit']")
    page.wait_for_url("**/workstation**")
    page.locator("[data-tab-target='queue-tab']").click()
    page.wait_for_selector("#queue-tab.active")
    page.reload()
    page.wait_for_selector("#queue-tab.active")
    assert "黄浦区公共卫生病原基因数据库和生物信息学分析系统" in page.title()
