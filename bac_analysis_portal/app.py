from __future__ import annotations

import csv
import gzip
import io
import json
import os
import platform
import shutil
import sys
import time
import zipfile
from datetime import datetime
from pathlib import Path
from urllib.parse import quote
from xml.sax.saxutils import escape as xml_escape

from flask import Flask, Response, jsonify, redirect, render_template, request, session, url_for
from functools import wraps

from .task_manager import ASM_METHOD_OPTIONS, AnalysisTaskManager, ValidationError
from .store import PortalStore


def create_app() -> Flask:
    package_root = Path(__file__).resolve().parent
    project_root = package_root.parent
    cpu_cache: dict[str, float] = {}
    app = Flask(
        __name__,
        template_folder=str(package_root / "templates"),
        static_folder=str(package_root / "static"),
    )
    app.secret_key = "bac-analysis-portal-dev-key"
    task_manager = AnalysisTaskManager.from_project_root(project_root)
    store = PortalStore.from_project_root(project_root)

    @app.before_request
    def protect_routes():
        if request.endpoint in {"login", "login_post", "static"}:
            return None
        if not _is_logged_in():
            if request.path.startswith("/api/"):
                return jsonify({"error": "Authentication required"}), 401
            return redirect(url_for("login"))
        return None

    @app.get("/login")
    def login():
        if _is_logged_in():
            return redirect(url_for("index"))
        return render_template("login.html")

    @app.post("/login")
    def login_post():
        payload = request.get_json(force=True) if request.is_json else request.form
        user = store.authenticate(str(payload.get("username", "")).strip(), str(payload.get("password", "")))
        if user is None:
            raise ValidationError("用户名或密码错误")
        session["username"] = user["username"]
        session["role"] = user["role"]
        session["group_name"] = user.get("group_name", "")
        return jsonify({"status": "ok", "user": user})

    @app.post("/logout")
    @login_required
    def logout():
        session.clear()
        return jsonify({"status": "logged_out"})

    @app.get("/")
    @login_required
    def index():
        return render_template(
            "index.html",
            project_root=str(project_root),
            script_path=store.get_setting("pipeline_script", "Bac_assemble_260112_newformat.py"),
        )

    @app.get("/api/health")
    @login_required
    def health():
        return jsonify(
            {
                "status": "ok",
                "workspace_root": store.get_setting("workspace_root", str(project_root)),
                "script_path": store.get_setting("pipeline_script", "Bac_assemble_260112_newformat.py"),
            }
        )

    @app.get("/api/session")
    @login_required
    def current_session():
        return jsonify(store.get_user(session["username"]))

    @app.get("/api/tasks")
    @login_required
    def list_tasks():
        tasks = task_manager.list_tasks()
        visible = [task for task in tasks if _can_view_task(task)]
        return jsonify({"items": visible})

    @app.get("/api/asm-options")
    @login_required
    def asm_options():
        return jsonify({"items": ASM_METHOD_OPTIONS})

    @app.get("/api/server-status")
    @login_required
    def server_status():
        return jsonify(_collect_server_status(project_root, cpu_cache))

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
            items.append(
                {
                    "name": child.name,
                    "type": item_type,
                    "path": _to_browser_path(base_root, child),
                }
            )

        parent_relative = ""
        if current_path != current_path.parent:
            parent_relative = _to_browser_path(base_root, current_path.parent)

        return jsonify(
            {
                "root": str(base_root),
                "current_path": str(current_path),
                "relative_path": _to_browser_path(base_root, current_path),
                "parent_relative_path": parent_relative,
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
        name = str(payload.get("name", "")).strip()
        if not name:
            raise ValidationError("文件夹名称不能为空")
        if any(token in name for token in ["/", "\\"]):
            raise ValidationError("文件夹名称不能包含路径分隔符")
        target = (parent / name).resolve()
        _assert_within_root(base_root, target)
        target.mkdir(exist_ok=False)
        return jsonify({"status": "ok", "path": _to_browser_path(base_root, target)})

    @app.post("/api/filesystem/rename")
    @login_required
    def filesystem_rename():
        payload = request.get_json(force=True)
        base_root = project_root.resolve()
        source = _resolve_browser_path(base_root, payload.get("path", ""))
        if not source.exists():
            raise ValidationError(f"待重命名目标不存在: {source}")
        new_name = str(payload.get("name", "")).strip()
        if not new_name:
            raise ValidationError("新名称不能为空")
        if any(token in new_name for token in ["/", "\\"]):
            raise ValidationError("新名称不能包含路径分隔符")
        target = (source.parent / new_name).resolve()
        _assert_within_root(base_root, target)
        source.rename(target)
        return jsonify({"status": "ok", "path": _to_browser_path(base_root, target)})

    @app.get("/api/tasks/<task_id>")
    @login_required
    def get_task(task_id: str):
        log_lines = request.args.get("log_lines", default=120, type=int)
        task = task_manager.get_task(task_id, log_lines=log_lines, owner=None)
        _ensure_can_view_task(task)
        result_html = _resolve_task_result_html(task)
        task["result_exists"] = result_html is not None
        task["result_url"] = url_for("task_result_view", task_id=task_id) if result_html else ""
        task["result_name"] = result_html.name if result_html else ""
        return jsonify(task)

    @app.delete("/api/tasks/<task_id>")
    @login_required
    def delete_task(task_id: str):
        task = task_manager.get_task(task_id, log_lines=0, owner=None)
        _ensure_can_modify_task(task)
        task_manager.delete_task(task_id, owner=None)
        return jsonify({"status": "deleted", "id": task_id})

    @app.get("/api/tasks/<task_id>/result-view")
    @login_required
    def task_result_view(task_id: str):
        task = task_manager.get_task(task_id, log_lines=0, owner=None)
        _ensure_can_view_task(task)
        result_html = _resolve_task_result_html(task)
        if result_html is None or not result_html.is_file():
            raise KeyError(f"Result not found for task: {task_id}")
        html = result_html.read_text(encoding="utf-8", errors="ignore")
        html = _inject_result_preview_style(html)
        return Response(html, mimetype="text/html")

    @app.get("/tasks/<task_id>/result-page")
    @login_required
    def task_result_page(task_id: str):
        task = task_manager.get_task(task_id, log_lines=0, owner=None)
        _ensure_can_view_task(task)
        return render_template("result_report.html", task=task)

    @app.get("/api/tasks/<task_id>/report-data")
    @login_required
    def task_report_data(task_id: str):
        task = task_manager.get_task(task_id, log_lines=0, owner=None)
        _ensure_can_view_task(task)
        return jsonify(_build_report_payload(task))

    @app.post("/api/export/table")
    @login_required
    def export_table():
        payload = request.get_json(force=True)
        title = str(payload.get("title", "")).strip() or "结果表"
        export_format = str(payload.get("format", "csv")).strip().lower()
        columns = _normalize_export_columns(payload.get("columns", []))
        rows = _normalize_export_rows(payload.get("rows", []))
        filename_root = _sanitize_export_filename(str(payload.get("filename", "")).strip() or title)

        if export_format == "csv":
            content = _build_delimited_bytes(columns, rows, ",")
            mimetype = "text/csv; charset=utf-8"
            extension = "csv"
        elif export_format == "tsv":
            content = _build_delimited_bytes(columns, rows, "\t")
            mimetype = "text/tab-separated-values; charset=utf-8"
            extension = "tsv"
        elif export_format == "xlsx":
            content = _build_xlsx_bytes(title, columns, rows)
            mimetype = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            extension = "xlsx"
        else:
            raise ValidationError("仅支持导出 csv、tsv 或 xlsx")

        filename = f"{filename_root}.{extension}"
        return Response(
            content,
            mimetype=mimetype,
            headers={
                "Content-Disposition": f"attachment; filename*=UTF-8''{quote(filename)}",
                "Cache-Control": "no-store",
            },
        )

    @app.post("/api/tasks")
    @login_required
    def create_task():
        payload = request.get_json(force=True)
        workspace_root = Path(store.get_setting("workspace_root", str(project_root))).expanduser().resolve()
        created = task_manager.create_task(
            payload,
            owner=session["username"],
            owner_group=str(session.get("group_name", "") or ""),
            pipeline_script=_resolve_pipeline_script(
                workspace_root,
                store.get_setting("pipeline_script", "Bac_assemble_260112_newformat.py"),
            ),
            pipeline_python=_resolve_runtime_python_path(
                workspace_root,
                store.get_setting("pipeline_python", sys.executable),
            ),
        )
        return jsonify(created), 201

    @app.post("/api/tasks/demo")
    @login_required
    def create_demo_task():
        created = task_manager.create_demo_task(owner=session["username"], owner_group=str(session.get("group_name", "") or ""))
        return jsonify(created), 201

    @app.post("/api/batch-inputs")
    @login_required
    def create_batch_input():
        payload = request.get_json(force=True)
        rows = payload.get("rows", [])
        if not isinstance(rows, list) or not rows:
            raise ValidationError("批量输入不能为空")

        normalized_rows: list[list[str]] = []
        for index, row in enumerate(rows, start=1):
            if not isinstance(row, dict):
                raise ValidationError(f"第 {index} 行格式不正确")
            sample_name = str(row.get("sample_name", "")).strip()
            species = str(row.get("species", "")).strip()
            third_gen = _resolve_optional_existing_path(project_root, row.get("third_gen"))
            short_left = _resolve_optional_existing_path(project_root, row.get("short_left"))
            short_right = _resolve_optional_existing_path(project_root, row.get("short_right"))
            if not sample_name:
                raise ValidationError(f"第 {index} 行样本名称不能为空")
            if not any([third_gen, short_left, short_right]):
                raise ValidationError(f"第 {index} 行至少需要填写一项测序数据")
            normalized_rows.append([sample_name, third_gen, short_left, short_right, species])

        batch_dir = project_root / "generated_batch_inputs"
        batch_dir.mkdir(parents=True, exist_ok=True)
        filename = f"batch_input_{datetime.now().strftime('%Y%m%d%H%M%S')}.tsv"
        target = batch_dir / filename
        with target.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="	")
            writer.writerow(["sample_name", "third_gen", "short_left", "short_right", "species"])
            writer.writerows(normalized_rows)

        return jsonify(
            {
                "path": _to_browser_path(project_root.resolve(), target.resolve()),
                "absolute_path": str(target.resolve()),
                "rows": len(normalized_rows),
            }
        ), 201

    @app.get("/api/admin/settings")
    @admin_required
    def admin_settings():
        return jsonify(
            {
                "workspace_root": store.get_setting("workspace_root", str(project_root)),
                "pipeline_script": store.get_setting("pipeline_script", "Bac_assemble_260112_newformat.py"),
                "pipeline_python": store.get_setting("pipeline_python", sys.executable),
                "max_concurrent_tasks": int(store.get_setting("max_concurrent_tasks", "2") or "2"),
            }
        )

    @app.get("/api/admin/filesystem")
    @admin_required
    def admin_filesystem():
        selector = request.args.get("selector", default="workspace_root", type=str)
        root_arg = request.args.get("root", default="", type=str)
        path_arg = request.args.get("path", default="", type=str)
        browse_root = _resolve_admin_root(root_arg)
        current_path = _resolve_admin_path(browse_root, path_arg)
        if not current_path.is_dir():
            raise ValidationError(f"只能浏览目录: {current_path}")

        items = []
        try:
            children = sorted(current_path.iterdir(), key=lambda item: (not item.is_dir(), item.name.lower()))
        except PermissionError as exc:
            raise ValidationError(f"没有访问目录权限: {current_path}") from exc

        for child in children:
            if child.name.startswith(".") and child.name not in {".", ".."}:
                continue
            item_type = "directory" if child.is_dir() else "file"
            if selector == "workspace_root" and item_type != "directory":
                continue
            items.append({"name": child.name, "type": item_type, "path": str(child.resolve())})

        parent_path = ""
        if current_path != browse_root:
            parent_path = str(current_path.parent.resolve())

        return jsonify(
            {
                "root": str(browse_root),
                "current_path": str(current_path),
                "parent_path": parent_path,
                "selector": selector,
                "items": items,
            }
        )

    @app.put("/api/admin/settings")
    @admin_required
    def update_admin_settings():
        payload = request.get_json(force=True)
        workspace_root = str(payload.get("workspace_root", "")).strip()
        script_path = str(payload.get("pipeline_script", "")).strip()
        pipeline_python = str(payload.get("pipeline_python", "")).strip()
        max_concurrent_tasks = str(payload.get("max_concurrent_tasks", "")).strip() or "2"
        if not workspace_root:
            raise ValidationError("部署基准目录不能为空")
        if not script_path:
            raise ValidationError("脚本路径不能为空")
        if not pipeline_python:
            raise ValidationError("运行环境 Python 路径不能为空")
        try:
            max_task_count = max(1, int(max_concurrent_tasks))
        except ValueError as exc:
            raise ValidationError("最大同时运行任务数量必须是整数") from exc
        workspace_root_path = Path(workspace_root).expanduser().resolve()
        if not workspace_root_path.is_dir():
            raise ValidationError(f"部署基准目录不存在: {workspace_root_path}")
        candidate = Path(_resolve_pipeline_script(workspace_root_path, script_path))
        runtime_python = Path(_resolve_runtime_python_path(workspace_root_path, pipeline_python))
        if not candidate.is_file():
            raise ValidationError(f"脚本路径不存在: {candidate}")
        if not runtime_python.is_file():
            raise ValidationError(f"运行环境 Python 不存在: {runtime_python}")
        store.set_setting("workspace_root", str(workspace_root_path))
        store.set_setting("pipeline_script", str(candidate.relative_to(workspace_root_path)))
        store.set_setting("pipeline_python", str(runtime_python))
        store.set_setting("max_concurrent_tasks", str(max_task_count))
        return jsonify(
            {
                "workspace_root": str(workspace_root_path),
                "pipeline_script": str(candidate.relative_to(workspace_root_path)),
                "pipeline_python": str(runtime_python),
                "max_concurrent_tasks": max_task_count,
            }
        )

    @app.get("/api/admin/users")
    @admin_required
    def list_users():
        return jsonify({"items": store.list_users()})

    @app.post("/api/admin/users")
    @admin_required
    def create_user():
        payload = request.get_json(force=True)
        username = str(payload.get("username", "")).strip()
        password = str(payload.get("password", ""))
        role = str(payload.get("role", "user")).strip()
        if not username or not password:
            raise ValidationError("用户名和密码不能为空")
        created = store.create_user(
            username=username,
            password=password,
            role=role,
            display_name=str(payload.get("display_name", "")).strip(),
            group_name=str(payload.get("group_name", "")).strip(),
        )
        return jsonify(created), 201

    @app.put("/api/admin/users/<username>")
    @admin_required
    def update_user(username: str):
        payload = request.get_json(force=True)
        if username == "admin" and str(payload.get("role", "")).strip() in {"user", "group_admin"}:
            raise ValidationError("默认管理员账号不能降级")
        updated = store.update_user(
            username,
            new_username=str(payload.get("username", "")).strip() or None,
            role=str(payload.get("role", "")).strip() or None,
            group_name=str(payload.get("group_name", "")).strip() if "group_name" in payload else None,
            display_name=str(payload.get("display_name", "")).strip() if "display_name" in payload else None,
            new_password=str(payload.get("password", "")).strip() or None,
        )
        if session.get("username") == username:
            session["username"] = updated["username"]
            session["role"] = updated["role"]
            session["group_name"] = updated.get("group_name", "")
        return jsonify(updated)

    @app.delete("/api/admin/users/<username>")
    @admin_required
    def delete_user(username: str):
        if session.get("username") == username:
            raise ValidationError("不能删除当前登录用户")
        if username == "admin":
            raise ValidationError("默认管理员账号不能删除")
        store.delete_user(username)
        return jsonify({"status": "deleted", "username": username})

    @app.errorhandler(ValidationError)
    def handle_validation(error: ValidationError):
        return jsonify({"error": str(error)}), 400

    @app.errorhandler(ValueError)
    def handle_value_error(error: ValueError):
        return jsonify({"error": str(error)}), 400

    @app.errorhandler(KeyError)
    def handle_missing(error: KeyError):
        return jsonify({"error": str(error)}), 404

    @app.errorhandler(Exception)
    def handle_unknown(error: Exception):
        app.logger.exception("Bac analysis portal error")
        return jsonify({"error": str(error)}), 500

    return app


def login_required(view_func):
    @wraps(view_func)
    def wrapped(*args, **kwargs):
        if not _is_logged_in():
            return jsonify({"error": "Authentication required"}), 401
        return view_func(*args, **kwargs)

    return wrapped


def admin_required(view_func):
    @wraps(view_func)
    def wrapped(*args, **kwargs):
        if not _is_logged_in():
            return jsonify({"error": "Authentication required"}), 401
        if session.get("role") != "admin":
            return jsonify({"error": "Administrator permission required"}), 403
        return view_func(*args, **kwargs)

    return wrapped


def _is_logged_in() -> bool:
    return bool(session.get("username") and session.get("role"))


def _can_view_task(task: dict) -> bool:
    role = str(session.get("role") or "")
    group_name = str(session.get("group_name") or "")
    if role == "admin":
        return True
    if group_name and str(task.get("owner_group") or "") == group_name:
        return str(task.get("owner_group") or "") == group_name and bool(group_name)
    return False


def _ensure_can_view_task(task: dict) -> None:
    if not _can_view_task(task):
        raise KeyError(f"Task not found: {task.get('id', '-')}")


def _ensure_can_modify_task(task: dict) -> None:
    role = str(session.get("role") or "")
    group_name = str(session.get("group_name") or "")
    if role == "admin":
        return
    if role == "group_admin" and group_name and str(task.get("owner_group") or "") == group_name:
        return
    raise ValidationError("只有管理员或 group 管理可以删除组内任务")


def _resolve_browser_path(base_root: Path, path_arg: str) -> Path:
    raw = str(path_arg or "").strip()
    if not raw:
        return base_root
    candidate = Path(raw).expanduser()
    return candidate.resolve() if candidate.is_absolute() else (base_root / candidate).resolve()


def _to_browser_path(base_root: Path, path: Path) -> str:
    if path == base_root:
        return ""
    try:
        return str(path.relative_to(base_root))
    except ValueError:
        return str(path)


def _is_within_root(base_root: Path, candidate: Path) -> bool:
    try:
        candidate.relative_to(base_root)
        return True
    except ValueError:
        return False


def _assert_within_root(base_root: Path, candidate: Path) -> None:
    try:
        candidate.relative_to(base_root)
    except ValueError as exc:
        raise ValidationError(f"路径超出允许范围: {candidate}") from exc


def _resolve_optional_existing_path(project_root: Path, raw_value: object) -> str:
    text = str(raw_value or "").strip()
    if not text:
        return ""
    candidate = Path(text).expanduser()
    candidate = candidate.resolve() if candidate.is_absolute() else (project_root / candidate).resolve()
    if not candidate.exists():
        raise ValidationError(f"文件不存在: {candidate}")
    return str(candidate)


def _resolve_task_result_html(task: dict) -> Path | None:
    input_path = str((task.get("params") or {}).get("input_path", "")).strip()
    if not input_path:
        return None
    candidate = Path(input_path).expanduser().resolve()
    search_dir = candidate if candidate.is_dir() else candidate.parent
    if not search_dir.is_dir():
        return None
    matches = sorted(search_dir.glob("*_bacgenome.html"))
    return matches[0] if matches else None


def _resolve_runtime_python_path(project_root: Path, python_path: str) -> str:
    raw = str(python_path or "").strip()
    if not raw:
        raise ValidationError("运行环境 Python 路径不能为空")
    candidate = Path(raw).expanduser()
    candidate = candidate.resolve() if candidate.is_absolute() else (project_root / candidate).resolve()
    return str(candidate)


def _resolve_pipeline_script(project_root: Path, script_path: str) -> str:
    raw = str(script_path or "").strip()
    if not raw:
        raise ValidationError("脚本路径不能为空")
    candidate = Path(raw).expanduser()
    candidate = candidate.resolve() if candidate.is_absolute() else (project_root / candidate).resolve()
    _assert_within_root(project_root.resolve(), candidate)
    return str(candidate)

def _resolve_admin_root(root_arg: str) -> Path:
    raw = str(root_arg or "").strip()
    if not raw:
        return Path.home().resolve()
    candidate = Path(raw).expanduser().resolve()
    if not candidate.is_dir():
        raise ValidationError(f"目录不存在: {candidate}")
    return candidate


def _resolve_admin_path(root: Path, path_arg: str) -> Path:
    raw = str(path_arg or "").strip()
    if not raw:
        return root
    candidate = Path(raw).expanduser()
    return candidate.resolve() if candidate.is_absolute() else (root / candidate).resolve()


def _collect_server_status(project_root: Path, cpu_cache: dict[str, float]) -> dict:
    memory = _read_memory_status()
    disk = _read_disk_status(project_root)
    cpu = _read_cpu_status(cpu_cache)
    return {
        "hostname": platform.node() or "-",
        "platform": platform.platform(),
        "machine": platform.machine() or "-",
        "python": sys.version.split()[0],
        "cpu_count": os.cpu_count() or 0,
        "sampled_at": datetime.now().isoformat(timespec="seconds"),
        "cpu": cpu,
        "memory": memory,
        "disk": disk,
    }


def _read_cpu_status(cpu_cache: dict[str, float]) -> dict:
    load_avg = _read_load_average()
    percent, source = _read_cpu_percent(cpu_cache)
    detail = f"{os.cpu_count() or 0} 核 / 负载 {load_avg}" if load_avg != "-" else f"{os.cpu_count() or 0} 核"
    return {
        "percent": percent,
        "source": source,
        "detail": detail,
        "load_average": load_avg,
    }


def _read_cpu_percent(cpu_cache: dict[str, float]) -> tuple[float | None, str]:
    proc_stat = Path("/proc/stat")
    if proc_stat.is_file():
        try:
            with proc_stat.open("r", encoding="utf-8") as handle:
                first_line = handle.readline().strip()
            parts = first_line.split()
            if len(parts) >= 5 and parts[0] == "cpu":
                values = [float(item) for item in parts[1:]]
                idle = values[3] + (values[4] if len(values) > 4 else 0.0)
                total = sum(values)
                previous_total = cpu_cache.get("total")
                previous_idle = cpu_cache.get("idle")
                cpu_cache["total"] = total
                cpu_cache["idle"] = idle
                if previous_total is not None and previous_idle is not None and total > previous_total:
                    idle_delta = idle - previous_idle
                    total_delta = total - previous_total
                    percent = max(0.0, min(100.0, (1.0 - idle_delta / total_delta) * 100.0))
                    return round(percent, 1), "proc-stat"
        except OSError:
            pass

    try:
        load = os.getloadavg()[0]
        cpu_total = max(os.cpu_count() or 1, 1)
        percent = max(0.0, min(100.0, load / cpu_total * 100.0))
        return round(percent, 1), "load-average"
    except (AttributeError, OSError):
        return None, "unavailable"


def _read_load_average() -> str:
    try:
        load1, load5, load15 = os.getloadavg()
        return f"{load1:.2f} / {load5:.2f} / {load15:.2f}"
    except (AttributeError, OSError):
        return "-"


def _read_memory_status() -> dict:
    meminfo = Path("/proc/meminfo")
    if meminfo.is_file():
        values: dict[str, int] = {}
        try:
            for line in meminfo.read_text(encoding="utf-8", errors="ignore").splitlines():
                if ":" not in line:
                    continue
                key, raw_value = line.split(":", 1)
                amount = raw_value.strip().split()[0]
                if amount.isdigit():
                    values[key] = int(amount) * 1024
            total = values.get("MemTotal", 0)
            available = values.get("MemAvailable", values.get("MemFree", 0))
            used = max(total - available, 0)
            percent = (used / total * 100.0) if total else None
            return {
                "total": total,
                "used": used,
                "free": available,
                "percent": round(percent, 1) if percent is not None else None,
                "total_human": _human_bytes(total),
                "used_human": _human_bytes(used),
                "free_human": _human_bytes(available),
            }
        except OSError:
            pass

    total = _read_total_memory_fallback()
    return {
        "total": total,
        "used": 0,
        "free": total,
        "percent": 0.0 if total else None,
        "total_human": _human_bytes(total),
        "used_human": _human_bytes(0),
        "free_human": _human_bytes(total),
    }


def _read_total_memory_fallback() -> int:
    try:
        page_size = os.sysconf("SC_PAGE_SIZE")
        page_count = os.sysconf("SC_PHYS_PAGES")
        return int(page_size) * int(page_count)
    except (ValueError, OSError, AttributeError):
        return 0


def _read_disk_status(project_root: Path) -> dict:
    usage = shutil.disk_usage(project_root)
    used_percent = (usage.used / usage.total * 100.0) if usage.total else None
    return {
        "path": str(project_root),
        "total": usage.total,
        "used": usage.used,
        "free": usage.free,
        "used_percent": round(used_percent, 1) if used_percent is not None else None,
        "total_human": _human_bytes(usage.total),
        "used_human": _human_bytes(usage.used),
        "free_human": _human_bytes(usage.free),
    }


def _human_bytes(value: int | float) -> str:
    size = float(value or 0)
    units = ["B", "KB", "MB", "GB", "TB", "PB"]
    for unit in units:
        if size < 1024.0 or unit == units[-1]:
            if unit == "B":
                return f"{int(size)} {unit}"
            return f"{size:.1f} {unit}"
        size /= 1024.0
    return "0 B"


def _normalize_export_columns(columns: object) -> list[str]:
    if not isinstance(columns, list):
        raise ValidationError("导出表头格式不正确")
    return [str(value or "") for value in columns]


def _normalize_export_rows(rows: object) -> list[list[str]]:
    if not isinstance(rows, list):
        raise ValidationError("导出表格内容格式不正确")
    normalized: list[list[str]] = []
    for row in rows:
        if not isinstance(row, list):
            raise ValidationError("导出表格内容格式不正确")
        normalized.append([str(value or "") for value in row])
    return normalized


def _sanitize_export_filename(name: str) -> str:
    safe = "".join(char if char.isalnum() or char in {"-", "_"} else "_" for char in name.strip())
    return safe.strip("_") or "export"


def _build_delimited_bytes(columns: list[str], rows: list[list[str]], delimiter: str) -> bytes:
    buffer = io.StringIO()
    writer = csv.writer(buffer, delimiter=delimiter, lineterminator="\n")
    if columns:
        writer.writerow(columns)
    writer.writerows(rows)
    return buffer.getvalue().encode("utf-8-sig")


def _xlsx_column_name(index: int) -> str:
    result = ""
    current = index + 1
    while current > 0:
        current, remainder = divmod(current - 1, 26)
        result = chr(65 + remainder) + result
    return result


def _build_xlsx_bytes(title: str, columns: list[str], rows: list[list[str]]) -> bytes:
    sheet_name = _sanitize_export_filename(title)[:31] or "Sheet1"

    def build_row_xml(row_number: int, values: list[str]) -> str:
        cells = []
        for index, value in enumerate(values):
            reference = f"{_xlsx_column_name(index)}{row_number}"
            text = xml_escape(str(value or ""))
            cells.append(
                f'<c r="{reference}" t="inlineStr"><is><t xml:space="preserve">{text}</t></is></c>'
            )
        return f'<row r="{row_number}">{"".join(cells)}</row>'

    sheet_rows = []
    current_row = 1
    if columns:
        sheet_rows.append(build_row_xml(current_row, columns))
        current_row += 1
    for row in rows:
        sheet_rows.append(build_row_xml(current_row, row))
        current_row += 1

    workbook_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        f'<sheets><sheet name="{xml_escape(sheet_name)}" sheetId="1" r:id="rId1"/></sheets>'
        '</workbook>'
    )
    sheet_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        f'<sheetData>{"".join(sheet_rows)}</sheetData>'
        '</worksheet>'
    )
    content_types_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
        '<Default Extension="xml" ContentType="application/xml"/>'
        '<Override PartName="/xl/workbook.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
        '<Override PartName="/xl/worksheets/sheet1.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
        '<Override PartName="/xl/styles.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>'
        '</Types>'
    )
    rels_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" '
        'Target="xl/workbook.xml"/>'
        '</Relationships>'
    )
    workbook_rels_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" '
        'Target="worksheets/sheet1.xml"/>'
        '<Relationship Id="rId2" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" '
        'Target="styles.xml"/>'
        '</Relationships>'
    )
    styles_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        '<fonts count="1"><font><sz val="11"/><name val="Calibri"/></font></fonts>'
        '<fills count="1"><fill><patternFill patternType="none"/></fill></fills>'
        '<borders count="1"><border/></borders>'
        '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>'
        '<cellXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/></cellXfs>'
        '<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>'
        '</styleSheet>'
    )

    output = io.BytesIO()
    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("[Content_Types].xml", content_types_xml)
        archive.writestr("_rels/.rels", rels_xml)
        archive.writestr("xl/workbook.xml", workbook_xml)
        archive.writestr("xl/_rels/workbook.xml.rels", workbook_rels_xml)
        archive.writestr("xl/worksheets/sheet1.xml", sheet_xml)
        archive.writestr("xl/styles.xml", styles_xml)
    return output.getvalue()


def _build_report_payload(task: dict) -> dict:
    params = task.get("params", {})
    report_dir = _resolve_report_directory(task)
    sample_name = _resolve_report_sample_name(task, report_dir)
    summary_info = _read_summary_metrics(report_dir / "summary.tsv")
    checkm_info = _read_checkm_metrics(report_dir / f"{sample_name}.checkm.tsv") if sample_name else {}
    assembly_profile = _read_assembly_profile(report_dir / "Assem_info.tsv")
    contig_depth_relationship = _read_contig_depth_relationship(report_dir / "Assem_info.tsv")
    fastp_info = _read_fastp_metrics(report_dir / f"{sample_name}.fastp2.json") if sample_name else {}
    assembly_summary = _read_tsv_rows(report_dir / f"{sample_name}.assemble.result.tsv") if sample_name else {"columns": [], "rows": []}
    assembly_coverage = _read_coverage_profile(report_dir / f"{sample_name}_ngs.per-base.bed.gz") if sample_name else {}
    contig_annotation = _read_tsv_rows(report_dir / "flye_output" / "assembly_info.txt")
    checkm_quality = _read_checkm2_quality(report_dir / "checkm2_out" / "quality_report.tsv")
    gene_annotation_summary = _read_tsv_rows(report_dir / f"{sample_name}.genefun_summary.tsv") if sample_name else {"columns": [], "rows": []}
    gene_length_distribution = _read_gene_length_distribution(report_dir / f"{sample_name}_gene_raw_sum.tsv") if sample_name else {"status": "empty", "points": []}
    rv_summary = _read_tsv_rows(report_dir / "Assem_info1.tsv")
    virulence_elements = _merge_annotation_with_summary(
        report_dir / "Assem_abricate_VFDB.tsv",
        report_dir / "VFDB_summary.tsv",
        mode="virulence",
    )
    virulence_relationship = _build_category_gene_relationship(
        virulence_elements,
        left_key="VF分类",
        right_key="基因名称",
        label="VF 分类与毒力基因关系图",
    )
    resistance_elements = _merge_annotation_with_summary(
        report_dir / "Assem_abricate_CARD.tsv",
        report_dir / "CARD_summary.tsv",
        mode="resistance",
    )
    resistance_relationship = _build_category_gene_relationship(
        resistance_elements,
        left_key="耐药药物",
        right_key="基因名称",
        label="耐药药物与耐药基因关系图",
        split_delimiters=[";"],
    )
    rv_overview = _build_resistance_virulence_summary(rv_summary, virulence_elements, resistance_elements)
    species_taxonomy = _read_taxonomy_list(report_dir / f"{sample_name}_2.list.txt", terminal_column="种") if sample_name else {"rows": [], "rank_options": []}
    subspecies_taxonomy = _read_taxonomy_list(report_dir / f"{sample_name}_2.list2.txt", terminal_column="亚种") if sample_name else {"rows": [], "rank_options": []}
    taxonomy_abundance = _build_taxonomy_abundance(species_taxonomy, subspecies_taxonomy)
    taxonomy_risk_summary = _build_taxonomy_risk_summary(species_taxonomy, subspecies_taxonomy)
    assembly_taxonomy = _read_assembly_taxonomy(
        report_dir / "Assem_info1.tsv",
        report_dir / f"{sample_name}_assem.kraken2.txt",
    ) if sample_name else {"columns": [], "rows": []}
    mlst_result = _read_mlst_result(report_dir / f"{sample_name}.mlst_Stat.txt") if sample_name else {"columns": [], "rows": [], "gene_show_map": {}, "default_gene": ""}
    serotype_result = _read_tsv_rows(report_dir / f"{sample_name}_serotype_result.tsv") if sample_name else {"columns": [], "rows": []}
    priority_serotype = _read_tsv_rows(report_dir / f"{sample_name}.pathonet_result.tsv") if sample_name else {"columns": [], "rows": []}
    return {
        "task": {
            "id": task.get("id"),
            "name": task.get("name"),
            "status": task.get("status"),
            "owner": task.get("owner"),
            "group": task.get("owner_group", ""),
            "created_at": task.get("created_at"),
            "started_at": task.get("started_at"),
            "finished_at": task.get("finished_at"),
            "input_path": params.get("input_path", ""),
            "output_dir": params.get("output_dir", ""),
            "asm_type": params.get("asm_type", ""),
            "method": params.get("method", ""),
            "species": params.get("species", ""),
            "sample_name": sample_name,
        },
        "overview_metrics": [
            {"key": "total_bases", "label": "总测序数据量", "type": "single", "value": summary_info.get("sum_len"), "display": _human_bp(summary_info.get("sum_len")), "unit": "bp"},
            {
                "key": "assembly_profile",
                "label": "组装情况",
                "type": "assembly_profile",
                "contig_count": assembly_profile.get("contig_count"),
                "plasmid_count": assembly_profile.get("plasmid_count"),
                "total_count": assembly_profile.get("total_count"),
                "total_length": assembly_profile.get("total_length"),
            },
            {
                "key": "q_metrics",
                "label": "Q20 / Q30",
                "type": "paired",
                "items": [
                    {"label": "Q20", "display": _display_percent(summary_info.get("q20_rate"))},
                    {"label": "Q30", "display": _display_percent(summary_info.get("q30_rate"))},
                ],
            },
            {
                "key": "checkm_metrics",
                "label": "完整性 / 污染率",
                "type": "paired",
                "items": [
                    {"label": "完整性", "display": _display_percent(checkm_info.get("completeness"))},
                    {"label": "污染率", "display": _display_percent(checkm_info.get("contamination"))},
                ],
            },
            {
                "key": "species_estimation",
                "label": "物种预估",
                "type": "paired",
                "items": [
                    {"label": "物种名称", "display": checkm_info.get("species_name") or "--"},
                    {"label": "MLST 物种名称", "display": checkm_info.get("mlst_species_name") or "--"},
                ],
            },
        ],
        "sections": {
            "overview": {"status": "empty"},
            "raw_qc": {
                "status": "ready" if fastp_info else "empty",
                "paired_end": {
                    "left": fastp_info.get("read1", {"status": "empty"}),
                    "right": fastp_info.get("read2", {"status": "empty"}),
                },
                "fastp": fastp_info.get("summary", {"status": "empty", "plots": []}),
            },
            "species_identification": {
                "status": "ready" if species_taxonomy.get("rows") or subspecies_taxonomy.get("rows") or assembly_taxonomy.get("rows") else "empty",
                "species": species_taxonomy,
                "subspecies": subspecies_taxonomy,
                "abundance": taxonomy_abundance,
                "risk_summary": taxonomy_risk_summary,
                "assembly_taxonomy": assembly_taxonomy,
            },
            "assembly": {
                "status": "empty",
                "summary": assembly_summary,
                "coverage": assembly_coverage,
                "contig_annotation": contig_annotation,
                "contig_depth_relationship": contig_depth_relationship,
                "checkm": checkm_quality,
                "gene_annotation_summary": gene_annotation_summary,
                "gene_length_distribution": gene_length_distribution,
            },
            "resistance_virulence": {
                "status": "empty",
                "overview": rv_overview,
                "summary": rv_summary,
                "virulence_elements": virulence_elements,
                "virulence_relationship": virulence_relationship,
                "resistance_elements": resistance_elements,
                "resistance_relationship": resistance_relationship,
            },
            "mlst": {
                "status": "ready" if mlst_result.get("rows") else "empty",
                "columns": mlst_result.get("columns", []),
                "rows": mlst_result.get("rows", []),
                "gene_show_map": mlst_result.get("gene_show_map", {}),
                "default_gene": mlst_result.get("default_gene", ""),
            },
            "serotype": {
                "status": "ready" if serotype_result.get("rows") else "empty",
                "columns": serotype_result.get("columns", []),
                "rows": serotype_result.get("rows", []),
            },
            "priority_serotype": {
                "status": "ready" if priority_serotype.get("rows") else "empty",
                "columns": priority_serotype.get("columns", []),
                "rows": priority_serotype.get("rows", []),
            },
        },
    }


def _resolve_report_directory(task: dict) -> Path:
    input_path = str((task.get("params") or {}).get("input_path", "")).strip()
    if not input_path:
        return Path(task.get("project_root") or ".").resolve()
    candidate = Path(input_path).expanduser().resolve()
    return candidate if candidate.is_dir() else candidate.parent


def _resolve_report_sample_name(task: dict, report_dir: Path) -> str:
    params = task.get("params") or {}
    task_name = str(task.get("name") or params.get("task_name") or "").strip()
    input_path = str(params.get("input_path") or "").strip()
    input_name = Path(input_path).name if input_path else ""
    if input_name == "fastq" and (report_dir / "Men-IGT.fastp2.json").is_file():
        return "Men-IGT"
    for pattern in ("*.checkm.tsv", "*.fastp2.json", "*_bacgenome.html"):
        matches = sorted(report_dir.glob(pattern))
        if matches:
            name = matches[0].name
            if name.endswith(".checkm.tsv"):
                return name[:-10]
            if name.endswith(".fastp2.json"):
                return name[:-11]
            if name.endswith("_bacgenome.html"):
                return name[:-14]
    if task_name and task_name != "demo_fastq":
        return task_name
    return Path(input_name).stem if input_name else ""


def _read_summary_metrics(path: Path) -> dict:
    if not path.is_file():
        return {}
    total_sum_len = 0
    weighted_q20 = 0.0
    weighted_q30 = 0.0
    try:
        with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                value = _safe_int(row.get("sum_len"))
                if value is not None:
                    total_sum_len += value
                    q20 = _safe_float(row.get("Q20(%)"))
                    q30 = _safe_float(row.get("Q30(%)"))
                    if q20 is not None:
                        weighted_q20 += q20 * value
                    if q30 is not None:
                        weighted_q30 += q30 * value
    except OSError:
        return {}
    return {
        "sum_len": total_sum_len or None,
        "q20_rate": round(weighted_q20 / total_sum_len, 2) if total_sum_len else None,
        "q30_rate": round(weighted_q30 / total_sum_len, 2) if total_sum_len else None,
    }


def _read_checkm_metrics(path: Path) -> dict:
    if not path.is_file():
        return {}
    try:
        with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="	")
            first = next(reader, None)
            if not first:
                return {}
            return {
                "contamination": first.get("污染率") or first.get("contamination") or None,
                "completeness": first.get("完整性") or first.get("completeness") or None,
                "species_name": first.get("物种名称") or first.get("species_name") or None,
                "mlst_species_name": first.get("mlst 物种名称") or first.get("mlst_species_name") or None,
            }
    except OSError:
        return {}


def _read_assembly_profile(path: Path) -> dict:
    if not path.is_file():
        return {"contig_count": None, "plasmid_count": None, "total_count": None, "total_length": None}
    contig_count = 0
    plasmid_count = 0
    total_count = 0
    total_length = 0
    try:
        with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
            reader = csv.reader(handle, delimiter="	")
            next(reader, None)
            for row in reader:
                if len(row) < 5:
                    continue
                total_count += 1
                length_value = _safe_int(row[1] if len(row) > 1 else None)
                if length_value is not None:
                    total_length += length_value
                genome_type = str(row[4]).lower()
                if "plasmid" in genome_type:
                    plasmid_count += 1
                else:
                    contig_count += 1
    except OSError:
        return {"contig_count": None, "plasmid_count": None, "total_count": None, "total_length": None}
    return {
        "contig_count": contig_count,
        "plasmid_count": plasmid_count,
        "total_count": total_count,
        "total_length": total_length or None,
    }


def _read_contig_depth_relationship(path: Path) -> dict:
    if not path.is_file():
        return {"status": "empty", "points": []}
    points: list[dict[str, object]] = []
    scatter_points: list[dict[str, object]] = []
    try:
        with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="	")
            for row in reader:
                contig_name = str(row.get("序列名称", "")).strip()
                depth = _safe_float(row.get("平均深度"))
                length = _safe_int(row.get("序列长度"))
                raw_type = str(row.get("基因组/质粒", "")).strip()
                if not contig_name or depth is None:
                    continue
                seq_type = "质粒" if "plasmid" in raw_type.lower() else "基因组"
                points.append({
                    "name": contig_name,
                    "depth": round(depth, 2),
                    "type": seq_type,
                })
                if length is not None and length > 0:
                    scatter_points.append({
                        "name": contig_name,
                        "depth": round(depth, 2),
                        "length": length,
                        "type": seq_type,
                    })
    except OSError:
        return {"status": "empty", "points": []}
    return {
        "status": "ready" if points else "empty",
        "label": "基因组/质粒与平均深度关系图",
        "x_label": "序列类型",
        "y_label": "平均深度",
        "points": points,
        "length_depth_scatter": {
            "status": "ready" if scatter_points else "empty",
            "label": "Contig长度与平均测序深度散点图",
            "x_label": "Contig长度(bp)",
            "y_label": "平均测序深度",
            "points": scatter_points,
        },
    }


def _read_coverage_profile(path: Path, target_points: int = 800) -> dict:
    if not path.is_file():
        return {"status": "empty", "points": []}
    depths: list[int] = []
    contigs: set[str] = set()
    try:
        with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 4:
                    continue
                contigs.add(parts[0])
                depth = _safe_float(parts[3])
                if depth is None:
                    continue
                depths.append(int(round(depth)))
    except OSError:
        return {"status": "empty", "points": []}
    if not depths:
        return {"status": "empty", "points": []}

    total_bases = len(depths)
    if total_bases <= target_points:
        points = [round(value, 2) for value in depths]
        x_values = list(range(1, total_bases + 1))
    else:
        chunk_size = max(1, (total_bases + target_points - 1) // target_points)
        points = []
        x_values = []
        for start in range(0, total_bases, chunk_size):
            chunk = depths[start:start + chunk_size]
            if not chunk:
                continue
            points.append(round(sum(chunk) / len(chunk), 2))
            x_values.append(start + 1)
    max_depth = max(points) if points else None
    mean_depth = round(sum(depths) / total_bases, 2) if depths else None
    x_ticks = [1, max(1, total_bases // 2), total_bases]
    return {
        "status": "ready",
        "label": "基因组覆盖度",
        "x_label": "基因组位置",
        "y_label": "测序深度",
        "points": points,
        "x_values": x_values,
        "total_bases": total_bases,
        "contig_count": len(contigs),
        "max_depth": max_depth,
        "mean_depth": mean_depth,
        "x_ticks": x_ticks,
    }


def _parse_mlst_line(line: str) -> list[str]:
    parts: list[str] = []
    current: list[str] = []
    depth = 0
    for char in line.strip():
        if char == ',' and depth == 0:
            parts.append(''.join(current).strip())
            current = []
            continue
        if char == '(':
            depth += 1
        elif char == ')' and depth > 0:
            depth -= 1
        current.append(char)
    if current:
        parts.append(''.join(current).strip())
    return parts


def _read_mlst_result(path: Path) -> dict:
    if not path.is_file():
        return {"columns": [], "rows": [], "gene_show_map": {}, "default_gene": ""}
    raw = _read_tsv_rows(path)
    if not raw["columns"] or not raw["rows"]:
        return {"columns": [], "rows": [], "gene_show_map": {}, "default_gene": ""}

    rows: list[list[str]] = []
    gene_show_map: dict[str, str] = {}
    default_gene = ""

    for record in (
        {raw["columns"][index]: row[index] if index < len(row) else "" for index in range(len(raw["columns"]))}
        for row in raw["rows"]
    ):
        host_gene_display = str(record.get("管家基因", "")).strip()
        allele_no = str(record.get("管家基因序号", "")).strip()
        if not host_gene_display:
            continue
        host_gene_key = host_gene_display if "_" in host_gene_display else (
            f"{host_gene_display}_{allele_no}" if allele_no else host_gene_display
        )
        display_gene_name = host_gene_display.rsplit("_", 1)[0] if "_" in host_gene_display else host_gene_display
        derived_allele = host_gene_display.rsplit("_", 1)[1] if "_" in host_gene_display else ""
        show_path = path.parent / f"{host_gene_key}_gene_show.txt"
        show_text = show_path.read_text(encoding="utf-8", errors="ignore") if show_path.is_file() else ""
        if show_text and not default_gene:
            default_gene = host_gene_key
        gene_show_map[host_gene_key] = show_text
        rows.append([
            display_gene_name,
            allele_no or derived_allele,
            str(record.get("序列名称", "")).strip(),
            str(record.get("起始位置", "")).strip(),
            str(record.get("终止位置", "")).strip(),
            str(record.get("比对起始位置", "")).strip(),
            str(record.get("比对终止位置", "")).strip(),
            str(record.get("比对长度", "")).strip(),
            str(record.get("一致性%", "")).strip(),
            str(record.get("序列分型(ST)", "")).strip(),
            str(record.get("物种信息", "")).strip(),
            host_gene_display,
            host_gene_key,
        ])

    return {
        "columns": ["Host Gene", "等位基因", "序列名称", "起始位置", "终止位置", "比对起始位置", "比对终止位置", "比对长度", "一致性%", "序列分型(ST)", "物种信息", "Host Gene 展示", "Host Gene ID"],
        "rows": rows,
        "gene_show_map": gene_show_map,
        "default_gene": default_gene or (rows[0][12] if rows else ""),
    }


def _read_taxonomy_list(path: Path, terminal_column: str) -> dict:
    raw = _read_tsv_rows(path)
    if not raw["columns"] or not raw["rows"]:
        return {"rows": [], "rank_options": []}
    rank_options = [column for column in ["界", "门", "纲", "目", "科", "属"] if column in raw["columns"]]
    rows = []
    for row in raw["rows"]:
        record = {raw["columns"][index]: row[index] if index < len(row) else "" for index in range(len(raw["columns"]))}
        record["比例数值"] = _safe_float(record.get("比例")) or 0.0
        record["序列数量数值"] = _safe_int(record.get("序列数量")) or 0
        rows.append(record)
    return {"rows": rows, "rank_options": rank_options, "terminal_column": terminal_column}


def _build_taxonomy_abundance(species_taxonomy: dict, subspecies_taxonomy: dict) -> dict:
    rank_order = ["界", "门", "纲", "目", "科", "属", "种", "亚种"]
    rank_sources = {
        "界": species_taxonomy,
        "门": species_taxonomy,
        "纲": species_taxonomy,
        "目": species_taxonomy,
        "科": species_taxonomy,
        "属": species_taxonomy,
        "种": species_taxonomy,
        "亚种": subspecies_taxonomy,
    }
    palette = [
        "#526a86", "#76834f", "#8a6654", "#6d6481", "#4e7b75",
        "#9b7a3f", "#8a4d47", "#5d7c83", "#7b6d5a", "#697789", "#8d8d8d",
    ]
    ranks: list[dict] = []
    for rank in rank_order:
        dataset = rank_sources.get(rank) or {}
        rows = dataset.get("rows") or []
        if not rows:
            continue
        groups: dict[str, dict[str, float | int]] = {}
        for row in rows:
            name = str(row.get(rank, "")).strip() or "未注释"
            record = groups.setdefault(name, {"ratio": 0.0, "reads": 0})
            record["ratio"] += float(row.get("比例数值") or 0.0)
            record["reads"] += int(row.get("序列数量数值") or 0)
        ranked = sorted(groups.items(), key=lambda item: (item[1]["ratio"], item[1]["reads"]), reverse=True)
        segments = []
        for index, (name, values) in enumerate(ranked):
            segments.append({
                "name": name,
                "ratio": round(float(values["ratio"]), 2),
                "reads": int(values["reads"]),
                "color": palette[index % len(palette)],
            })
        ranks.append({
            "rank": rank,
            "segments": segments,
            "total_ratio": round(sum(segment["ratio"] for segment in segments), 2),
        })
    return {"status": "ready" if ranks else "empty", "ranks": ranks}


def _build_taxonomy_risk_summary(species_taxonomy: dict, subspecies_taxonomy: dict) -> dict:
    datasets = [species_taxonomy or {}, subspecies_taxonomy or {}]
    pathogenicity: dict[str, dict[str, float | int]] = {}
    hazard: dict[str, dict[str, float | int]] = {}
    kingdom_groups: dict[str, dict[str, float | int]] = {
        "细菌": {"reads": 0, "records": 0},
        "病毒": {"reads": 0, "records": 0},
        "真菌": {"reads": 0, "records": 0},
    }
    total_reads = 0
    total_records = 0

    def _normalize_kingdom(value: str) -> str | None:
        text = value.strip().lower()
        if not text or text == "-":
            return None
        if any(token in text for token in ("细菌", "bacteria", "eubacteria")):
            return "细菌"
        if any(token in text for token in ("病毒", "virus", "viruses")):
            return "病毒"
        if any(token in text for token in ("真菌", "fungi", "fungus", "mycota")):
            return "真菌"
        return None

    for dataset in datasets:
        for row in dataset.get("rows") or []:
            reads = int(row.get("序列数量数值") or 0)
            ratio = float(row.get("比例数值") or 0.0)
            total_reads += reads
            total_records += 1

            kingdom_label = _normalize_kingdom(str(row.get("界", "")).strip())
            if kingdom_label:
                kingdom_groups[kingdom_label]["reads"] += reads
                kingdom_groups[kingdom_label]["records"] += 1

            pathogenic_label = str(row.get("致病性", "")).strip()
            if pathogenic_label and pathogenic_label != "-":
                record = pathogenicity.setdefault(pathogenic_label, {"reads": 0, "ratio": 0.0, "records": 0})
                record["reads"] += reads
                record["ratio"] += ratio
                record["records"] += 1

            hazard_label = str(row.get("危害程度等级", "")).strip()
            if hazard_label and hazard_label != "-":
                record = hazard.setdefault(hazard_label, {"reads": 0, "ratio": 0.0, "records": 0})
                record["reads"] += reads
                record["ratio"] += ratio
                record["records"] += 1

    def _format_rank(mapping: dict[str, dict[str, float | int]]) -> list[dict]:
        return [
            {
                "label": label,
                "reads": int(values["reads"]),
                "ratio": round(float(values["ratio"]), 2),
                "records": int(values["records"]),
            }
            for label, values in sorted(mapping.items(), key=lambda item: (item[1]["reads"], item[1]["ratio"]), reverse=True)
        ]

    pathogenicity_rows = _format_rank(pathogenicity)
    hazard_rows = _format_rank(hazard)
    dominant_pathogenicity = pathogenicity_rows[0]["label"] if pathogenicity_rows else ""
    dominant_hazard = hazard_rows[0]["label"] if hazard_rows else ""

    if pathogenicity_rows or hazard_rows:
        summary_parts = []
        if dominant_pathogenicity:
            summary_parts.append(f"当前序列物种鉴定结果以“{dominant_pathogenicity}”类型为主")
        if dominant_hazard:
            summary_parts.append(f"危害程度以“{dominant_hazard}”为主要等级")
        if total_reads:
            summary_parts.append(f"纳入统计的分类结果覆盖 {total_reads} 条序列")
        narrative = "，".join(summary_parts) + "。"
    else:
        narrative = "当前序列物种鉴定结果未提供可汇总的致病性或危害程度等级信息。"

    return {
        "status": "ready" if pathogenicity_rows or hazard_rows else "empty",
        "kingdom_summary": [
            {
                "label": label,
                "reads": int(values["reads"]),
                "records": int(values["records"]),
                "ratio": round((int(values["reads"]) / total_reads) * 100, 2) if total_reads else 0.0,
            }
            for label, values in kingdom_groups.items()
        ],
        "pathogenicity": pathogenicity_rows,
        "hazard": hazard_rows,
        "narrative": narrative,
        "total_reads": total_reads,
        "total_records": total_records,
    }


def _read_assembly_taxonomy(assem_info_path: Path, kraken_report_path: Path) -> dict:
    assem_info = _read_tsv_rows(assem_info_path)
    if not assem_info["columns"] or not assem_info["rows"]:
        return {"columns": [], "rows": []}

    rank_mapping = {
        "D": "界",
        "P": "门",
        "C": "纲",
        "O": "目",
        "F": "科",
        "G": "属",
        "S": "种",
    }
    tracked_ranks = ["界", "门", "纲", "目", "科", "属", "种"]
    lineage_by_taxid: dict[str, dict[str, str]] = {}
    lineage_stack: dict[str, str] = {}

    if kraken_report_path.is_file():
        try:
            with kraken_report_path.open("r", encoding="utf-8", errors="ignore") as handle:
                for line in handle:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 6:
                        continue
                    rank_code = str(parts[3]).strip()
                    taxid = str(parts[4]).strip()
                    raw_name = str(parts[5]).rstrip()
                    clean_name = raw_name.strip()
                    if not taxid or not clean_name:
                        continue

                    normalized_rank = rank_code[0] if rank_code else ""
                    if normalized_rank not in rank_mapping:
                        continue

                    current_rank = rank_mapping[normalized_rank]
                    current_index = tracked_ranks.index(current_rank)
                    for stale_rank in tracked_ranks[current_index:]:
                        lineage_stack.pop(stale_rank, None)
                    lineage_stack[current_rank] = clean_name

                    lineage_by_taxid[taxid] = {rank: lineage_stack.get(rank, "") for rank in tracked_ranks}
        except OSError:
            lineage_by_taxid = {}

    assem_records = [
        {assem_info["columns"][index]: row[index] if index < len(row) else "" for index in range(len(assem_info["columns"]))}
        for row in assem_info["rows"]
    ]

    columns = [
        "序列名称",
        "序列长度",
        "平均深度",
        "是否成环",
        "基因组/质粒",
        "质粒分型",
        "taxid",
        "物种名称",
        "界",
        "门",
        "纲",
        "目",
        "科",
        "属",
    ]
    rows: list[list[str]] = []
    for record in assem_records:
        taxid = str(record.get("taxid", "")).strip()
        lineage = lineage_by_taxid.get(taxid, {})
        rows.append([
            str(record.get("序列名称", "")).strip(),
            str(record.get("序列长度", "")).strip(),
            str(record.get("平均深度", "")).strip(),
            str(record.get("是否成环", "")).strip(),
            str(record.get("基因组/质粒", "")).strip(),
            str(record.get("质粒分型", "")).strip(),
            taxid,
            str(record.get("物种名称", "")).strip(),
            str(lineage.get("界", "")).strip(),
            str(lineage.get("门", "")).strip(),
            str(lineage.get("纲", "")).strip(),
            str(lineage.get("目", "")).strip(),
            str(lineage.get("科", "")).strip(),
            str(lineage.get("属", "")).strip(),
        ])
    return {"columns": columns, "rows": rows}


def _read_gene_length_distribution(path: Path) -> dict:
    raw = _read_tsv_rows(path)
    if not raw["columns"] or not raw["rows"]:
        return {"status": "empty", "points": []}
    x_values: list[str] = []
    points: list[int] = []
    for row in raw["rows"]:
        label = row[0] if row else ""
        value = _safe_int(row[1] if len(row) > 1 else None)
        if not label or value is None:
            continue
        x_values.append(str(label))
        points.append(value)
    if not points:
        return {"status": "empty", "points": []}
    return {
        "status": "ready",
        "label": "基因长度与数量分布",
        "x_label": raw["columns"][0] if raw["columns"] else "基因长度范围",
        "y_label": raw["columns"][1] if len(raw["columns"]) > 1 else "Gene数量",
        "x_values": x_values,
        "points": points,
        "max_count": max(points),
    }


def _build_resistance_virulence_summary(rv_summary: dict, virulence_elements: dict, resistance_elements: dict) -> dict:
    rv_rows = rv_summary.get("rows") or []
    rv_columns = rv_summary.get("columns") or []
    virulence_rows = virulence_elements.get("rows") or []
    resistance_rows = resistance_elements.get("rows") or []

    def _count_by_column(rows: list[list[str]], columns: list[str], target: str, split_delimiters: list[str] | None = None) -> list[dict]:
        if target not in columns:
            return []
        index = columns.index(target)
        counts: dict[str, int] = {}
        for row in rows:
            raw = str(row[index] if index < len(row) else "").strip()
            if not raw or raw == '-':
                continue
            parts = [raw]
            if split_delimiters:
                parts = [raw]
                for delimiter in split_delimiters:
                    expanded = []
                    for item in parts:
                        expanded.extend(item.split(delimiter))
                    parts = expanded
            for item in [part.strip() for part in parts if part.strip() and part.strip() != '-']:
                counts[item] = counts.get(item, 0) + 1
        return [
            {"label": label, "count": count}
            for label, count in sorted(counts.items(), key=lambda item: item[1], reverse=True)[:8]
        ]

    def _rv_metric(column_name: str) -> int:
        if column_name not in rv_columns:
            return 0
        idx = rv_columns.index(column_name)
        total = 0
        for row in rv_rows:
            try:
                total += int(float(str(row[idx]).strip() or '0'))
            except ValueError:
                continue
        return total

    virulence_top = _count_by_column(virulence_rows, virulence_elements.get("columns") or [], "VF分类")
    resistance_top = _count_by_column(resistance_rows, resistance_elements.get("columns") or [], "耐药药物", split_delimiters=[';'])
    virulence_genes = len({str(row[(virulence_elements.get("columns") or []).index("基因名称")]).strip() for row in virulence_rows if "基因名称" in (virulence_elements.get("columns") or []) and str(row[(virulence_elements.get("columns") or []).index("基因名称")]).strip() not in {'', '-'}}) if virulence_rows else 0
    resistance_genes = len({str(row[(resistance_elements.get("columns") or []).index("基因名称")]).strip() for row in resistance_rows if "基因名称" in (resistance_elements.get("columns") or []) and str(row[(resistance_elements.get("columns") or []).index("基因名称")]).strip() not in {'', '-'}}) if resistance_rows else 0

    resistance_note = f"当前共检出 {len(resistance_rows)} 条耐药元件记录、{resistance_genes} 个耐药基因，主要集中于“{resistance_top[0]['label']}”相关类别。" if resistance_top else f"当前共检出 {len(resistance_rows)} 条耐药元件记录、{resistance_genes} 个耐药基因。"
    virulence_note = f"当前共检出 {len(virulence_rows)} 条毒力元件记录、{virulence_genes} 个毒力基因，主要集中于“{virulence_top[0]['label']}”类别。" if virulence_top else f"当前共检出 {len(virulence_rows)} 条毒力元件记录、{virulence_genes} 个毒力基因。"

    return {
        "status": "ready" if rv_rows or virulence_rows or resistance_rows else "empty",
        "resistance": {
            "note": resistance_note,
            "hit_count": len(resistance_rows),
            "gene_count": resistance_genes,
            "top_categories": resistance_top,
            "summary_count": _rv_metric("耐药基因数量"),
        },
        "virulence": {
            "note": virulence_note,
            "hit_count": len(virulence_rows),
            "gene_count": virulence_genes,
            "top_categories": virulence_top,
            "summary_count": _rv_metric("毒力基因数量"),
        },
    }


def _merge_annotation_with_summary(detail_path: Path, summary_path: Path, mode: str) -> dict:
    detail = _read_tsv_rows(detail_path)
    summary = _read_tsv_rows(summary_path)
    if not detail["columns"] and not summary["columns"]:
        return {"columns": [], "rows": []}

    detail_records = [
        {detail["columns"][index]: row[index] if index < len(row) else "" for index in range(len(detail["columns"]))}
        for row in detail["rows"]
    ]
    summary_records = [
        {summary["columns"][index]: row[index] if index < len(row) else "" for index in range(len(summary["columns"]))}
        for row in summary["rows"]
    ]
    summary_index: dict[tuple[str, str], dict] = {}
    for record in summary_records:
        key = (
            str(record.get("基因名称", "")).strip(),
            str(record.get("片段名称", "")).strip(),
        )
        summary_index[key] = record

    if mode == "virulence":
        columns = [
            "Contig名称", "物种名称", "taxid", "基因名称", "覆盖度%", "一致性%", "产物",
            "VF分类", "VF名称", "起始碱基", "终止碱基", "正负链",
            "覆盖度(>0)%", "覆盖度(>10)%", "覆盖度(>100)%", "平均深度", "最低深度", "最高深度",
        ]
    else:
        columns = [
            "Contig名称", "物种名称", "taxid", "基因名称", "覆盖度%", "一致性%", "产物",
            "耐药药物", "起始碱基", "终止碱基", "正负链",
            "覆盖度(>0)%", "覆盖度(>10)%", "覆盖度(>100)%", "平均深度", "最低深度", "最高深度",
        ]

    rows = []
    for record in detail_records:
        summary_record = summary_index.get(
            (
                str(record.get("基因名称", "")).strip(),
                str(record.get("Contig名称", "")).strip(),
            ),
            {},
        )
        merged = []
        for column in columns:
            if column in record:
                merged.append(record.get(column, ""))
            elif column in summary_record:
                merged.append(summary_record.get(column, ""))
            else:
                merged.append("")
        rows.append(merged)

    if rows:
        return {"columns": columns, "rows": rows}

    fallback_columns = summary["columns"] or detail["columns"]
    fallback_rows = summary["rows"] or detail["rows"]
    return {"columns": fallback_columns, "rows": fallback_rows}


def _build_category_gene_relationship(data: dict, left_key: str, right_key: str, label: str, split_delimiters: list[str] | None = None) -> dict:
    columns = data.get("columns") or []
    rows = data.get("rows") or []
    if not columns or not rows:
        return {"status": "empty", "nodes_left": [], "nodes_right": [], "links": []}
    left_index = columns.index(left_key) if left_key in columns else -1
    right_index = columns.index(right_key) if right_key in columns else -1
    if left_index < 0 or right_index < 0:
        return {"status": "empty", "nodes_left": [], "nodes_right": [], "links": []}
    split_delimiters = split_delimiters or []
    link_counts: dict[tuple[str, str], int] = {}
    left_totals: dict[str, int] = {}
    right_totals: dict[str, int] = {}
    for row in rows:
        left_value = str(row[left_index] if left_index < len(row) else "").strip()
        right_value = str(row[right_index] if right_index < len(row) else "").strip()
        if not left_value or not right_value or left_value == "-" or right_value == "-":
            continue
        left_values = [left_value]
        for delimiter in split_delimiters:
            if delimiter in left_value:
                left_values = [item.strip() for item in left_value.split(delimiter) if item.strip()]
                break
        for left_item in left_values:
            key = (left_item, right_value)
            link_counts[key] = link_counts.get(key, 0) + 1
            left_totals[left_item] = left_totals.get(left_item, 0) + 1
            right_totals[right_value] = right_totals.get(right_value, 0) + 1
    if not link_counts:
        return {"status": "empty", "nodes_left": [], "nodes_right": [], "links": []}
    top_left = sorted(left_totals.items(), key=lambda item: (-item[1], item[0]))[:8]
    top_right = sorted(right_totals.items(), key=lambda item: (-item[1], item[0]))[:12]
    keep_left = {name for name, _ in top_left}
    keep_right = {name for name, _ in top_right}
    links = [
        {"source": left, "target": right, "value": count}
        for (left, right), count in sorted(link_counts.items(), key=lambda item: (-item[1], item[0][0], item[0][1]))
        if left in keep_left and right in keep_right
    ]
    return {
        "status": "ready" if links else "empty",
        "label": label,
        "left_label": left_key,
        "right_label": right_key,
        "nodes_left": [{"name": name, "value": value} for name, value in top_left],
        "nodes_right": [{"name": name, "value": value} for name, value in top_right],
        "links": links,
    }


def _read_tsv_rows(path: Path) -> dict:
    if not path.is_file():
        return {"columns": [], "rows": []}
    try:
        with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
            reader = csv.reader(handle, delimiter="	")
            rows = list(reader)
    except OSError:
        return {"columns": [], "rows": []}
    if not rows:
        return {"columns": [], "rows": []}
    header = rows[0]
    data_rows = rows[1:]
    return {"columns": header, "rows": data_rows}


def _read_checkm2_quality(path: Path) -> dict:
    raw = _read_tsv_rows(path)
    if not raw["columns"]:
        return {"columns": [], "rows": []}
    rename_map = {
        "Name": "样本名称",
        "Completeness": "完整性",
        "Contamination": "污染率",
        "Completeness_Model_Used": "完整性模型",
        "Translation_Table_Used": "翻译表",
        "Coding_Density": "编码密度",
        "Contig_N50": "Contig N50",
        "Average_Gene_Length": "平均基因长度",
        "Genome_Size": "基因组大小",
        "GC_Content": "GC 含量",
        "Total_Coding_Sequences": "编码序列数",
        "Total_Contigs": "Contig 总数",
        "Max_Contig_Length": "最大 Contig 长度",
        "Additional_Notes": "附加说明",
    }
    keep_columns = [
        "Name", "Completeness", "Contamination", "Genome_Size",
        "GC_Content", "Total_Contigs", "Contig_N50", "Max_Contig_Length",
        "Coding_Density", "Total_Coding_Sequences", "Completeness_Model_Used", "Additional_Notes",
    ]
    indices = [raw["columns"].index(col) for col in keep_columns if col in raw["columns"]]
    columns = [rename_map.get(raw["columns"][index], raw["columns"][index]) for index in indices]
    rows = []
    for row in raw["rows"]:
        trimmed = []
        for index in indices:
            value = row[index] if index < len(row) else ""
            if raw["columns"][index] in {"Completeness", "Contamination"} and value not in {"", "None", None}:
                try:
                    value = f"{float(value):.2f}%"
                except ValueError:
                    pass
            elif raw["columns"][index] == "GC_Content" and value not in {"", "None", None}:
                try:
                    numeric = float(value)
                    value = f"{numeric * 100:.2f}%" if numeric <= 1 else f"{numeric:.2f}%"
                except ValueError:
                    pass
            elif raw["columns"][index] == "Genome_Size" and value not in {"", "None", None}:
                value = _human_bp(value)
            elif raw["columns"][index] == "Coding_Density" and value not in {"", "None", None}:
                try:
                    numeric = float(value)
                    value = f"{numeric * 100:.2f}%" if numeric <= 1 else f"{numeric:.2f}%"
                except ValueError:
                    pass
            trimmed.append(value)
        rows.append(trimmed)
    return {"columns": columns, "rows": rows}


def _read_fastp_metrics(path: Path) -> dict:
    if not path.is_file():
        return {}
    try:
        payload = json.loads(path.read_text(encoding="utf-8", errors="ignore"))
    except (OSError, json.JSONDecodeError):
        return {}
    summary = payload.get("summary", {})
    before = summary.get("before_filtering", {})
    after = summary.get("after_filtering", {})
    filtering = payload.get("filtering_result", {})
    duplication = payload.get("duplication", {})
    read1_before = payload.get("read1_before_filtering", {})
    read2_before = payload.get("read2_before_filtering", {})
    read1_after = payload.get("read1_after_filtering", {})
    read2_after = payload.get("read2_after_filtering", {})
    return {
        "summary": {
            "status": "ready",
            "sequencing": summary.get("sequencing", ""),
            "before_filtering": before,
            "after_filtering": after,
            "filtering_result": filtering,
            "duplication_rate": duplication.get("rate"),
            "insert_size": payload.get("insert_size", {}),
            "adapter_cutting": payload.get("adapter_cutting", {}),
            "base_distribution": {
                "read1": read1_after.get("content_curves", {}),
                "read2": read2_after.get("content_curves", {}),
            },
        },
        "read1": {
            "status": "ready",
            "label": "R1",
            "before": read1_before,
            "after": read1_after,
            "quality_curves": read1_before.get("quality_curves", {}),
            "content_curves": read1_before.get("content_curves", {}),
            "before_summary": _summarize_fastp_read_block(read1_before),
            "after_summary": _summarize_fastp_read_block(read1_before),
        },
        "read2": {
            "status": "ready",
            "label": "R2",
            "before": read2_before,
            "after": read2_after,
            "quality_curves": read2_before.get("quality_curves", {}),
            "content_curves": read2_before.get("content_curves", {}),
            "before_summary": _summarize_fastp_read_block(read2_before),
            "after_summary": _summarize_fastp_read_block(read2_before),
        },
    }


def _summarize_fastp_read_block(section: dict) -> dict:
    total_reads = _safe_int(section.get("total_reads"))
    total_bases = _safe_int(section.get("total_bases"))
    q20_bases = _safe_int(section.get("q20_bases"))
    q30_bases = _safe_int(section.get("q30_bases"))
    content_curves = section.get("content_curves", {}) or {}
    gc_curve = content_curves.get("GC", []) or []
    mean_length = round(total_bases / total_reads, 2) if total_reads and total_bases else None
    q20_rate = round(q20_bases / total_bases, 6) if q20_bases is not None and total_bases else None
    q30_rate = round(q30_bases / total_bases, 6) if q30_bases is not None and total_bases else None
    gc_content = round(sum((_safe_float(value) or 0) for value in gc_curve) / len(gc_curve), 6) if gc_curve else None
    return {
        "total_reads": total_reads,
        "total_bases": total_bases,
        "mean_length": mean_length,
        "q20_rate": q20_rate,
        "q30_rate": q30_rate,
        "gc_content": gc_content,
    }


def _safe_int(value: object) -> int | None:
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return None


def _safe_float(value: object) -> float | None:
    try:
        return float(str(value))
    except (TypeError, ValueError):
        return None


def _human_bp(value: object) -> str:
    number = _safe_int(value)
    if number is None:
        return "--"
    units = [("bp", 1), ("Kb", 10**3), ("Mb", 10**6), ("Gb", 10**9)]
    unit = "bp"
    divisor = 1
    for unit_name, base in units:
        if number >= base:
            unit = unit_name
            divisor = base
    scaled = number / divisor
    if divisor == 1:
        return f"{number} {unit}"
    return f"{scaled:.2f} {unit}"


def _display_percent(value: object) -> str:
    numeric = _safe_float(value)
    if numeric is None:
        return "--"
    return f"{numeric:.2f}%"


def _inject_result_preview_style(html: str) -> str:
    preview_style = """
<style id="portal-result-preview-style">
html, body {
  margin: 0 !important;
  padding: 0 !important;
  width: 100% !important;
  min-height: 100% !important;
  overflow-x: auto !important;
  background: #ffffff !important;
}
body {
  font-size: 16px !important;
}
.container-fluid.main-container,
div.main-container {
  width: min(1680px, calc(100vw - 96px)) !important;
  max-width: min(1680px, calc(100vw - 96px)) !important;
  margin: 0 auto !important;
  padding: 18px 24px 32px !important;
  box-sizing: border-box !important;
}
.container-fluid.main-container > .row {
  display: flex !important;
  align-items: flex-start !important;
  gap: 24px !important;
  margin: 0 !important;
}
.container-fluid.main-container > .row::before,
.container-fluid.main-container > .row::after {
  display: none !important;
}
.container-fluid.main-container > .row > [class*="col-"] {
  float: none !important;
  padding-left: 0 !important;
  padding-right: 0 !important;
}
.container-fluid.main-container > .row > .col-xs-12.col-sm-4.col-md-3 {
  flex: 0 0 280px !important;
  width: 280px !important;
  max-width: 280px !important;
}
.container-fluid.main-container > .row > .toc-content {
  flex: 1 1 auto !important;
  width: auto !important;
  max-width: none !important;
  min-width: 0 !important;
  padding: 0 !important;
}
div.tocify, #section-TOC {
  position: sticky !important;
  top: 82px !important;
  width: 280px !important;
  max-width: 280px !important;
  max-height: calc(100vh - 100px) !important;
  margin: 0 !important;
  overflow: auto !important;
}
.html-widget, .plotly, .datatables, table, img, svg, canvas {
  max-width: 100% !important;
}
pre {
  white-space: pre-wrap !important;
  word-break: break-word !important;
}
@media (max-width: 1024px) {
  .container-fluid.main-container,
  div.main-container {
    width: calc(100vw - 28px) !important;
    max-width: calc(100vw - 28px) !important;
    padding: 12px !important;
  }
  .container-fluid.main-container > .row {
    display: block !important;
  }
  .container-fluid.main-container > .row > .col-xs-12.col-sm-4.col-md-3,
  .container-fluid.main-container > .row > .toc-content {
    width: 100% !important;
    max-width: none !important;
  }
  div.tocify, #section-TOC {
    position: relative !important;
    top: auto !important;
    width: 100% !important;
    max-width: none !important;
    max-height: none !important;
    margin-bottom: 16px !important;
  }
}
</style>
"""
    if 'id="portal-result-preview-style"' in html:
        return html
    if '</head>' in html:
        return html.replace('</head>', f'{preview_style}</head>', 1)
    return preview_style + html


def _inject_result_page_chrome(html: str, task: dict) -> str:
    title = str(task.get("name") or task.get("id") or "任务结果")
    status = str(task.get("status") or "未知")
    chrome_style = f"""
<style id="portal-result-page-chrome">
body {{
  padding-top: 64px !important;
  scroll-padding-top: 76px !important;
}}
#portal-result-page-bar {{
  position: fixed;
  top: 0;
  left: 0;
  right: 0;
  z-index: 9999;
  height: 64px;
  display: flex;
  align-items: center;
  justify-content: space-between;
  gap: 12px;
  padding: 0 18px;
  border-bottom: 1px solid rgba(35, 54, 88, 0.12);
  background: rgba(246, 248, 252, 0.98);
  backdrop-filter: blur(10px);
  box-sizing: border-box;
}}
#portal-result-page-left {{
  display: flex;
  align-items: center;
  gap: 12px;
  min-width: 0;
}}
#portal-result-float-back {{
  position: fixed;
  top: 12px;
  left: 18px;
  z-index: 10002;
  border: 1px solid rgba(35, 54, 88, 0.14);
  background: #ffffff;
  color: #17233a;
  border-radius: 14px;
  padding: 10px 14px;
  font: inherit;
  cursor: pointer;
  box-shadow: 0 10px 26px rgba(23, 35, 58, 0.12);
}}
#portal-result-title {{
  font-size: 20px;
  font-weight: 700;
  color: #17233a;
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
}}
#portal-result-status {{
  display: inline-flex;
  align-items: center;
  justify-content: center;
  min-width: 96px;
  height: 40px;
  padding: 0 14px;
  border-radius: 999px;
  background: #355c94;
  color: #fff;
  font-weight: 700;
  font-size: 14px;
}}
@media (max-width: 720px) {{
  body {{
    padding-top: 84px !important;
  }}
  #portal-result-page-bar {{
    height: 84px;
    padding: 10px 12px;
    align-items: flex-start;
  }}
  #portal-result-page-left {{
    flex-direction: column;
    align-items: flex-start;
    gap: 8px;
    padding-top: 38px;
  }}
  #portal-result-float-back {{
    top: 10px;
    left: 12px;
    padding: 9px 12px;
  }}
}}
</style>
<script>
window.addEventListener('DOMContentLoaded', function () {{
  var backButton = document.getElementById('portal-result-float-back');
  if (backButton) {{
    backButton.addEventListener('click', function () {{
      window.location.href = '/?tab=queue'
    }});
  }}
}});
</script>
"""
    chrome_body = f"""
<button id="portal-result-float-back" type="button">返回任务列表</button>
<div id="portal-result-page-bar">
  <div id="portal-result-page-left">
    <div id="portal-result-title">{title}</div>
  </div>
  <div id="portal-result-status">{status}</div>
</div>
"""
    if '</head>' in html and 'portal-result-page-chrome' not in html:
        html = html.replace('</head>', chrome_style + '\n</head>', 1)
    if '<body>' in html and 'portal-result-page-bar' not in html:
        html = html.replace('<body>', '<body>\n' + chrome_body, 1)
    return html


def _build_result_placeholder_page(task: dict) -> str:
    title = str(task.get("name") or task.get("id") or "任务结果")
    status = str(task.get("status") or "未知")
    return f"""<!doctype html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{title}</title>
  <style>
    body {{ margin: 0; font-family: sans-serif; background: #f3f6fb; color: #17233a; }}
    .bar {{ display:flex; justify-content:space-between; align-items:center; padding:12px 18px; border-bottom:1px solid rgba(35,54,88,.12); background:rgba(246,248,252,.98); }}
    .left {{ display:flex; align-items:center; gap:12px; }}
    button {{ border:1px solid rgba(35,54,88,.14); background:#fff; border-radius:12px; padding:10px 14px; cursor:pointer; }}
    .status {{ display:inline-flex; align-items:center; justify-content:center; min-width:96px; height:40px; padding:0 14px; border-radius:999px; background:#6a7488; color:#fff; font-weight:700; }}
    .body {{ padding:40px 24px; }}
  </style>
</head>
<body>
  <div class="bar">
    <div class="left"><button onclick="window.location.href='/?tab=queue'">返回任务列表</button><strong>{title}</strong></div>
    <div class="status">{status}</div>
  </div>
  <div class="body">当前没有可展示的结果页面。</div>
</body>
</html>"""

if __name__ == "__main__":
    app = create_app()
    app.run(host="0.0.0.0", port=5055, debug=True)
