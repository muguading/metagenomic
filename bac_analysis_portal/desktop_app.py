from __future__ import annotations

import json
import socket
import threading
import time
from contextlib import closing
from base64 import b64encode
from pathlib import Path
from urllib.parse import urlparse

from flask import Flask, jsonify, request
from werkzeug.serving import make_server

from .app import create_app

APP_NAME = "黄浦区公共卫生病原基因数据库和生物信息学分析系统"
DESKTOP_CONFIG_PATH = Path.home() / ".bac_analysis_portal_desktop.json"
DEFAULT_LOCAL_TARGET = "127.0.0.1:5055"
CONNECTION_WINDOW_HTML = """
<!doctype html>
<html lang="zh-CN">
  <head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>连接工作站</title>
    <style>
      :root {
        color-scheme: light;
        --ink: #1c2a44;
        --muted: #6b7891;
        --line: rgba(35, 54, 88, 0.12);
        --panel: rgba(249, 251, 254, 0.96);
        --brand: #355c94;
        --brand-strong: #27446d;
        --white: #ffffff;
      }
      * { box-sizing: border-box; }
      body {
        margin: 0;
        min-height: 100vh;
        font-family: "SF Pro Text", "PingFang SC", "Helvetica Neue", sans-serif;
        color: var(--ink);
        background: linear-gradient(180deg, rgba(247, 250, 253, 0.96), rgba(238, 243, 249, 0.98));
      }
      .login-shell {
        min-height: 100vh;
        display: grid;
        grid-template-columns: minmax(0, 1.08fr) minmax(420px, 0.92fr);
        align-items: stretch;
      }
      .login-brand-panel {
        display: grid;
        place-items: center;
        padding: clamp(28px, 4vw, 52px);
      }
      .login-brand-surface {
        width: min(100%, 760px);
        min-height: min(72vh, 820px);
        display: grid;
        place-items: center;
      }
      .login-brand-image {
        width: min(100%, 700px);
        max-height: min(74vh, 820px);
        object-fit: contain;
        display: block;
        filter: drop-shadow(0 18px 46px rgba(22, 40, 66, 0.12));
      }
      .login-card {
        display: flex;
        flex-direction: column;
        justify-content: center;
        padding: clamp(28px, 4vw, 54px);
        background: linear-gradient(180deg, rgba(252, 253, 255, 0.88), rgba(247, 250, 253, 0.92));
        border-left: 1px solid rgba(214, 224, 237, 0.6);
        backdrop-filter: blur(10px);
      }
      .eyebrow {
        margin: 0 0 10px;
        font-size: 0.82rem;
        letter-spacing: 0.18em;
        text-transform: uppercase;
        color: #41679d;
        font-weight: 700;
      }
      h1 {
        margin: 0;
        font-size: clamp(2rem, 4vw, 3rem);
        line-height: 1.08;
        font-family: "Newsreader", "Source Serif 4", serif;
      }
      .hero-copy {
        margin: 12px 0 0;
        font-size: 1rem;
        line-height: 1.72;
        color: var(--muted);
        max-width: 560px;
      }
      .login-form {
        display: grid;
        gap: 14px;
        margin-top: 24px;
        max-width: 560px;
      }
      .login-form label {
        display: grid;
        gap: 8px;
      }
      .login-form label span {
        font-size: 0.94rem;
        font-weight: 700;
        color: #44546f;
      }
      .field input,
      .login-form input {
        width: 100%;
        min-height: 52px;
        padding: 0 16px;
        border-radius: 16px;
        border: 1px solid rgba(35, 54, 88, 0.14);
        background: rgba(255, 255, 255, 0.96);
        font-size: 1rem;
        color: var(--ink);
        outline: none;
      }
      .field input:focus {
        border-color: rgba(53, 92, 148, 0.34);
        box-shadow: 0 0 0 4px rgba(53, 92, 148, 0.08);
      }
      .field {
        display: grid;
        gap: 8px;
      }
      .field label {
        display: block;
        font-size: 0.94rem;
        font-weight: 700;
        color: #485674;
      }
      .hint {
        color: #7a879d;
        font-size: 0.88rem;
      }
      .actions {
        display: flex;
        justify-content: flex-end;
        gap: 10px;
        margin-top: 6px;
      }
      button {
        min-height: 46px;
        padding: 0 18px;
        border-radius: 999px;
        font-size: 0.95rem;
        font-weight: 700;
        cursor: pointer;
        border: 1px solid transparent;
      }
      .ghost {
        background: rgba(255, 255, 255, 0.84);
        border-color: rgba(35, 54, 88, 0.1);
        color: #5a6780;
      }
      .primary {
        background: var(--brand);
        color: #fff;
      }
      .primary:hover { background: var(--brand-strong); }
      .error {
        display: none;
        margin-top: 14px;
        padding: 12px 14px;
        border-radius: 14px;
        background: rgba(154, 54, 54, 0.1);
        color: #873b3b;
        font-size: 0.92rem;
        line-height: 1.6;
      }
      .error.visible { display: block; }
      .login-note {
        margin-top: 20px;
        display: inline-flex;
        gap: 10px;
        align-items: center;
        flex-wrap: wrap;
        color: var(--muted);
        font-size: 0.92rem;
      }
      .login-note code {
        padding: 6px 10px;
        border-radius: 999px;
        background: rgba(255, 255, 255, 0.84);
        border: 1px solid rgba(35, 54, 88, 0.08);
        color: #355c94;
        font-size: 0.88rem;
        font-weight: 700;
      }
      @media (max-width: 640px) {
        .login-shell {
          grid-template-columns: 1fr;
        }
        .login-brand-panel {
          padding: 20px 20px 0;
        }
        .login-brand-surface {
          min-height: auto;
        }
        .login-brand-image {
          width: min(100%, 420px);
          max-height: 220px;
        }
        .login-card {
          padding: 26px 20px 28px;
          border-left: none;
          border-top: 1px solid rgba(214, 224, 237, 0.6);
        }
        .mode-grid { grid-template-columns: 1fr; }
        .actions { flex-direction: column-reverse; }
        button { width: 100%; }
      }
    </style>
  </head>
  <body>
    <main class="login-shell">
      <section class="login-brand-panel" aria-label="品牌展示">
        <div class="login-brand-surface">
          <img class="login-brand-image" src="__BRAND_IMAGE_SRC__" alt="黄浦区公共卫生病原基因数据库和生物信息学分析系统">
        </div>
      </section>
      <section class="login-card">
        <p class="eyebrow">Clinical Assembly Console</p>
        <h1>连接分析工作站</h1>
        <p class="hero-copy">统一在一个入口完成工作站地址连接与账号登录。默认填写本机地址，若分析系统部署在服务器上，直接改成服务器 IP 或域名即可。</p>
        <form class="login-form" onsubmit="event.preventDefault(); submitConnect();">
          <div class="field" id="server-field">
            <label for="server-input">工作站地址</label>
            <input id="server-input" type="text" autocomplete="off" placeholder="例如：127.0.0.1:5055 或 192.168.1.10:5055">
            <div class="hint">默认本机地址；若分析系统部署在服务器，请改为服务器 IP:端口 或完整 URL。</div>
          </div>
          <label>
            <span>用户名</span>
            <input id="username-input" type="text" autocomplete="username" placeholder="例如：admin">
          </label>
          <label>
            <span>密码</span>
            <input id="password-input" type="password" autocomplete="current-password" placeholder="请输入登录密码">
          </label>
          <div id="connection-error" class="error" role="alert"></div>
          <div class="actions">
            <button class="ghost" type="button" onclick="cancelConnect()">取消</button>
            <button class="primary" id="connect-button" type="submit">进入软件</button>
          </div>
        </form>
        <div class="login-note">
          <span>默认管理员</span>
          <code>admin / admin123</code>
        </div>
      </section>
    </main>
    <script>
      const initialTarget = "__INITIAL_TARGET__";
      const desktopBridgeUrl = "__DESKTOP_BRIDGE_URL__";

      const errorBox = document.getElementById("connection-error");
      const serverInput = document.getElementById("server-input");
      const usernameInput = document.getElementById("username-input");
      const passwordInput = document.getElementById("password-input");
      const connectButton = document.getElementById("connect-button");
      serverInput.value = initialTarget;

      function showError(message) {
        errorBox.textContent = message;
        errorBox.classList.add("visible");
      }

      function clearError() {
        errorBox.textContent = "";
        errorBox.classList.remove("visible");
      }

      function updateConnectButtonState() {
        connectButton.disabled = false;
        connectButton.textContent = "进入软件";
      }

      async function postDesktopBridge(path, payload) {
        const response = await fetch(`${desktopBridgeUrl}${path}`, {
          method: "POST",
          headers: {
            "Content-Type": "application/json",
          },
          body: JSON.stringify(payload || {}),
        });
        const result = await response.json().catch(() => null);
        if (!response.ok) {
          throw new Error((result && result.error) || "桌面接口调用失败，请稍后重试。");
        }
        return result || {};
      }

      async function submitConnect() {
        clearError();
        connectButton.disabled = true;
        connectButton.textContent = "连接中...";
        try {
          const result = await postDesktopBridge("/submit", {
            target: serverInput.value,
            username: usernameInput.value,
            password: passwordInput.value,
          });
          if (!result || !result.ok) {
            showError((result && result.error) || "连接失败，请检查地址或本机运行环境。");
            return;
          }
          window.location.replace(result.url);
        } catch (error) {
          showError((error && error.message) || "连接失败，请稍后重试。");
        } finally {
          updateConnectButtonState();
        }
      }

      async function cancelConnect() {
        try {
          await postDesktopBridge("/cancel", {});
        } catch (error) {
          window.close();
        }
      }

      [serverInput, usernameInput, passwordInput].forEach((input) => input.addEventListener("keydown", (event) => {
        if (event.key === "Enter") {
          event.preventDefault();
          submitConnect();
        }
      }));

      updateConnectButtonState();

    </script>
  </body>
</html>
"""


def _find_free_port(host: str = "127.0.0.1") -> int:
    with closing(socket.socket(socket.AF_INET, socket.SOCK_STREAM)) as sock:
        sock.bind((host, 0))
        sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        return int(sock.getsockname()[1])


class PortalServerThread(threading.Thread):
    def __init__(self, host: str, port: int):
        super().__init__(name="bac-analysis-portal-server", daemon=True)
        self.host = host
        self.port = port
        self.app = create_app()
        self.server = make_server(host, port, self.app, threaded=True)
        self.context = self.app.app_context()
        self.context.push()

    def run(self) -> None:
        self.server.serve_forever()

    def shutdown(self) -> None:
        try:
            self.server.shutdown()
        finally:
            self.context.pop()


class DesktopBridgeThread(threading.Thread):
    def __init__(self, host: str, port: int, controller: "ConnectionController"):
        super().__init__(name="bac-analysis-desktop-bridge", daemon=True)
        self.host = host
        self.port = port
        self.controller = controller
        self.app = Flask("bac_analysis_portal_desktop_bridge")
        self._register_routes()
        self.server = make_server(host, port, self.app, threaded=True)
        self.context = self.app.app_context()
        self.context.push()

    def _register_routes(self) -> None:
        @self.app.after_request
        def add_cors_headers(response):
            response.headers["Access-Control-Allow-Origin"] = "*"
            response.headers["Access-Control-Allow-Headers"] = "Content-Type"
            response.headers["Access-Control-Allow-Methods"] = "POST, OPTIONS"
            return response

        @self.app.route("/submit", methods=["POST", "OPTIONS"])
        def submit():
            if request.method == "OPTIONS":
                return ("", 204)
            payload = request.get_json(silent=True) or {}
            result = self.controller.submit(
                str(payload.get("target") or ""),
                str(payload.get("username") or ""),
                str(payload.get("password") or ""),
            )
            status_code = 200 if result.get("ok") else 400
            return jsonify(result), status_code

        @self.app.route("/cancel", methods=["POST", "OPTIONS"])
        def cancel():
            if request.method == "OPTIONS":
                return ("", 204)
            threading.Thread(target=self.controller.cancel, daemon=True).start()
            return jsonify({"ok": True})

    def run(self) -> None:
        self.server.serve_forever()

    def shutdown(self) -> None:
        try:
            self.server.shutdown()
        finally:
            self.context.pop()


def _read_desktop_config() -> dict:
    try:
        if DESKTOP_CONFIG_PATH.is_file():
            return json.loads(DESKTOP_CONFIG_PATH.read_text(encoding="utf-8"))
    except Exception:
        pass
    return {}


def _write_desktop_config(payload: dict) -> None:
    try:
        DESKTOP_CONFIG_PATH.write_text(
            json.dumps(payload, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
    except Exception:
        pass


def _normalize_remote_target(raw_value: str) -> str:
    value = (raw_value or "").strip()
    if not value:
        raise ValueError("请输入服务器 IP 或完整访问地址。")
    if not value.startswith(("http://", "https://")):
        value = f"http://{value}"
    parsed = urlparse(value)
    if not parsed.netloc:
        raise ValueError("服务器地址格式不正确，请填写如 192.168.1.10:5055。")
    return value


def _is_local_target(url: str) -> bool:
    parsed = urlparse(url)
    hostname = str(parsed.hostname or "").strip().lower()
    return hostname in {"127.0.0.1", "localhost", "::1"}


def _resolve_local_server_binding(raw_target: str) -> tuple[str, int]:
    normalized = _normalize_remote_target(raw_target)
    parsed = urlparse(normalized)
    host = "127.0.0.1"
    port = int(parsed.port or 5055)
    return host, port


def _build_brand_image_src(project_root: Path) -> str:
    logo_path = project_root / "bac_analysis_portal" / "static" / "image_v1.png"
    if not logo_path.is_file():
        return ""
    try:
        encoded = b64encode(logo_path.read_bytes()).decode("ascii")
        return f"data:image/png;base64,{encoded}"
    except Exception:
        return ""


class ConnectionController:
    def __init__(self) -> None:
        self.window = None
        self.server: PortalServerThread | None = None
        self.bridge: DesktopBridgeThread | None = None
        self.pending_credentials: dict[str, str] | None = None

    def bind_window(self, window: object) -> None:
        self.window = window

    def bind_bridge(self, bridge: DesktopBridgeThread) -> None:
        self.bridge = bridge

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
            server = PortalServerThread(host, port)
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
            return {"ok": False, "error": f"本机服务启动失败：{error}"}

    def cancel(self) -> None:
        if self.window is not None:
            try:
                self.window.destroy()
            except Exception:
                pass

    def autofill_login_if_needed(self) -> None:
        if self.window is None or not self.pending_credentials:
            return
        try:
            current_url = str(self.window.get_current_url() or "")
        except Exception:
            current_url = ""
        if "/login" not in current_url:
            return
        payload = json.dumps(self.pending_credentials, ensure_ascii=False)
        script = f"""
            (function() {{
              const creds = {payload};
              const username = document.getElementById('login-username');
              const password = document.getElementById('login-password');
              const form = document.getElementById('login-form');
              if (!username || !password || !form) {{
                return false;
              }}
              username.value = creds.username || '';
              password.value = creds.password || '';
              form.dispatchEvent(new Event('submit', {{ bubbles: true, cancelable: true }}));
              return true;
            }})();
        """
        try:
            time.sleep(0.15)
            self.window.evaluate_js(script)
            self.pending_credentials = None
        except Exception:
            pass


def launch_desktop_app() -> None:
    try:
        import webview
    except Exception as error:  # pragma: no cover - runtime dependency check
        raise RuntimeError(
            "未安装 pywebview，请先执行 `pip install -r requirements-web.txt`。"
        ) from error

    config = _read_desktop_config()
    initial_target = str(config.get("server_url") or DEFAULT_LOCAL_TARGET)
    project_root = Path(__file__).resolve().parent.parent
    brand_image_src = _build_brand_image_src(project_root)
    bridge_host = "127.0.0.1"
    bridge_port = _find_free_port(bridge_host)
    controller = ConnectionController()
    bridge = DesktopBridgeThread(bridge_host, bridge_port, controller)
    bridge.start()
    controller.bind_bridge(bridge)
    selector_html = (
        CONNECTION_WINDOW_HTML
        .replace("__INITIAL_TARGET__", initial_target.replace("\\", "\\\\").replace('"', "&quot;"))
        .replace("__DESKTOP_BRIDGE_URL__", f"http://{bridge_host}:{bridge_port}")
        .replace("__BRAND_IMAGE_SRC__", brand_image_src)
    )
    icon_candidates = [
        project_root / "bac_analysis_portal" / "static" / "app_icon.png",
        project_root / "bac_analysis_portal" / "static" / "favicon.png",
    ]
    icon_path = next((str(path) for path in icon_candidates if path.is_file()), None)
    window = webview.create_window(
        APP_NAME,
        html=selector_html,
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
        if controller.bridge is not None:
            controller.bridge.shutdown()
        if controller.server is not None:
            controller.server.shutdown()


__all__ = ["launch_desktop_app"]
