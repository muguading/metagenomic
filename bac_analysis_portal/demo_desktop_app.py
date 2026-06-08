from __future__ import annotations

import time
from pathlib import Path

from werkzeug.serving import make_server

from .demo_app import create_demo_app
from .desktop_app import (
    CONNECTION_WINDOW_HTML,
    DEFAULT_LOCAL_TARGET,
    ConnectionController,
    _build_brand_image_src,
    _is_local_target,
    _normalize_remote_target,
    _read_desktop_config,
    _resolve_local_server_binding,
    _write_desktop_config,
)


APP_NAME = "黄浦区公共卫生病原基因数据库和生物信息学分析系统 Demo"


class DemoPortalServerThread:
    def __init__(self, host: str, port: int):
        self.host = host
        self.port = port
        self.app = create_demo_app()
        self.server = make_server(host, port, self.app, threaded=True)
        self.context = self.app.app_context()
        self.context.push()

    def start(self) -> None:
        import threading

        self._thread = threading.Thread(target=self.server.serve_forever, name="bac-analysis-demo-server", daemon=True)
        self._thread.start()

    def shutdown(self) -> None:
        try:
            self.server.shutdown()
        finally:
            self.context.pop()


class DemoConnectionController(ConnectionController):
    def submit(self, target: str, username: str, password: str) -> dict:
        raw_target = str(target or "").strip()
        username_text = str(username or "").strip()
        password_text = str(password or "")
        if not username_text or not password_text:
            return {"ok": False, "error": "请输入用户名和密码。"}
        self.pending_credentials = {"username": username_text, "password": password_text}
        try:
            normalized_target = _normalize_remote_target(raw_target)
        except ValueError as error:
            return {"ok": False, "error": str(error)}
        _write_desktop_config({"server_url": raw_target.strip() or DEFAULT_LOCAL_TARGET})

        if not _is_local_target(normalized_target):
            base_url = normalized_target.rstrip("/")
            return {"ok": True, "url": f"{base_url}/login"}

        try:
            host, port = _resolve_local_server_binding(raw_target or DEFAULT_LOCAL_TARGET)
            server = DemoPortalServerThread(host, port)
            server.start()
            time.sleep(0.6)
            self.server = server
            return {"ok": True, "url": f"http://{host}:{port}/login"}
        except Exception as error:
            if self.server is not None:
                try:
                    self.server.shutdown()
                except Exception:
                    pass
                self.server = None
            return {"ok": False, "error": f"演示版本机服务启动失败：{error}"}


def launch_demo_desktop_app() -> None:
    try:
        import webview
    except Exception as error:  # pragma: no cover
        raise RuntimeError("未安装 pywebview，请先执行 `pip install -r requirements-web.txt`。") from error

    config = _read_desktop_config()
    initial_target = str(config.get("server_url") or DEFAULT_LOCAL_TARGET)
    project_root = Path(__file__).resolve().parent.parent
    brand_image_src = _build_brand_image_src(project_root)
    selector_html = (
        CONNECTION_WINDOW_HTML
        .replace("__INITIAL_TARGET__", initial_target.replace("\\", "\\\\").replace('"', "&quot;"))
        .replace("__BRAND_IMAGE_SRC__", brand_image_src)
    )
    icon_candidates = [
        project_root / "bac_analysis_portal" / "static" / "app_icon.png",
        project_root / "bac_analysis_portal" / "static" / "favicon.png",
    ]
    icon_path = next((str(path) for path in icon_candidates if path.is_file()), None)
    controller = DemoConnectionController()

    window = webview.create_window(
        APP_NAME,
        html=selector_html,
        js_api=controller,
        width=780,
        height=620,
        min_size=(720, 560),
        text_select=True,
        confirm_close=True,
        background_color="#edf2f8",
    )
    controller.bind_window(window)
    window.events.loaded += lambda: controller.autofill_login_if_needed()
    if icon_path and hasattr(window, "set_icon"):
        try:
            window.set_icon(icon_path)
        except Exception:
            pass

    try:
        webview.start(private_mode=False)
    finally:
        if controller.server is not None:
            controller.server.shutdown()


if __name__ == "__main__":
    launch_demo_desktop_app()
