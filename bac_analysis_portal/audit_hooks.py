from __future__ import annotations

import json
import uuid
from datetime import datetime

from flask import Response, g, request, session

from .store import PortalStore


def _sanitize(value: object, *, depth: int = 0) -> object:
    if depth > 3:
        return "..."
    if isinstance(value, dict):
        return {
            str(key): "***" if any(token in str(key).lower() for token in {"password", "token", "secret", "hash"}) else _sanitize(item, depth=depth + 1)
            for key, item in value.items()
        }
    if isinstance(value, (list, tuple)):
        return [_sanitize(item, depth=depth + 1) for item in value[:30]]
    text = str(value)
    return f"{text[:240]}..." if len(text) > 240 else text


def _request_summary() -> str:
    payload: dict[str, object] = {}
    json_data = request.get_json(silent=True)
    if isinstance(json_data, (dict, list)):
        payload["json"] = _sanitize(json_data)
    if request.form:
        payload["form"] = _sanitize(request.form.to_dict(flat=True))
    if request.files:
        payload["files"] = [str(name or "") for name in request.files.keys()]
    if request.args:
        payload["args"] = _sanitize(request.args.to_dict(flat=True))
    return json.dumps(payload, ensure_ascii=False, sort_keys=True) if payload else ""


def _response_summary(response: Response) -> str:
    try:
        data = response.get_json(silent=True) if response.is_json else None
        if isinstance(data, dict):
            focus = {
                key: data.get(key)
                for key in ("status", "error", "message", "id", "username", "path", "decision", "diagnosis_category", "diagnosis_label", "comparison_outcome", "comparison_outcome_label", "review_outcome", "review_outcome_label")
                if key in data
            }
            return json.dumps(_sanitize(focus), ensure_ascii=False, sort_keys=True)
    except Exception:
        pass
    return ""


def _identity() -> tuple[str, str, str]:
    actor = getattr(g, "audit_actor", None) or {}
    username = str(actor.get("username") or session.get("username") or "").strip()
    role = str(actor.get("role") or session.get("role") or "").strip()
    group_name = str(actor.get("group_name") or session.get("group_name") or "").strip()
    if username:
        return username, role, group_name
    payload = request.get_json(silent=True)
    fallback = str(payload.get("username") or "").strip() if isinstance(payload, dict) else str(request.form.get("username") or "").strip()
    return fallback, role, group_name


def _meta(path: str, method: str) -> tuple[str, str, str, str]:
    segments = [segment for segment in path.split("/") if segment]
    if "tasks" in segments:
        index = segments.index("tasks")
        target_id = segments[index + 1] if index + 1 < len(segments) else ""
        closure_action_map = {
            "review_report": "确认报告复核",
            "compare_history": "确认历史对照",
        }
        action = next((label for token, label in {"database-import": "导入数据库", "pause": "暂停任务", "resume": "恢复任务", "stop": "停止任务", "rerun": "重新运行", "rebuild": "重建任务"}.items() if token in segments), "任务操作")
        if "closure-actions" in segments:
            closure_index = segments.index("closure-actions")
            closure_action = segments[closure_index + 1] if closure_index + 1 < len(segments) else ""
            action = closure_action_map.get(closure_action, "确认闭环动作")
        if "report-exports" in segments:
            action = "确认导出归档"
        if "failure-disposition" in segments:
            payload = request.get_json(silent=True)
            action = "确认停止归档" if isinstance(payload, dict) and payload.get("decision") == "accept_stop" else "确认失败归档"
        if method == "POST" and segments[-1] == "tasks":
            action = "创建任务"
        elif method == "DELETE":
            action = "删除任务"
        return "任务", action, "task", target_id
    if "database" in segments and "samples" in segments:
        index = segments.index("samples")
        return "样本数据库", "样本库操作", "sample", segments[index + 1] if index + 1 < len(segments) else ""
    if "modeling" in segments:
        target_id = ""
        if "datasets" in segments:
            action = "创建建模数据集" if method == "POST" else "建模数据集操作"
            target_id = segments[segments.index("datasets") + 1] if segments.index("datasets") + 1 < len(segments) else ""
        elif "feature-sets" in segments:
            action = "创建特征方案" if method == "POST" else "特征方案操作"
        elif "train" in segments:
            action = "启动模型训练"
        elif "models" in segments:
            index = segments.index("models")
            target_id = segments[index + 1] if index + 1 < len(segments) else ""
            if "activate" in segments:
                action = "启用模型"
            elif "deactivate" in segments:
                action = "停用模型"
            elif "predict" in segments:
                action = "执行模型预测"
            elif method == "DELETE":
                action = "删除模型"
            else:
                action = "模型库操作"
        elif "predictions" in segments:
            action = "预测结果操作"
            target_id = segments[segments.index("predictions") + 1] if segments.index("predictions") + 1 < len(segments) else ""
        else:
            action = "建模平台操作"
        return "模型构建与预测", action, "modeling", target_id
    if "host-database" in segments:
        return "参考库", "宿主数据库操作", "reference", request.view_args.get("host_key", "") if request.view_args else ""
    if "pathogen-database" in segments:
        return "参考库", "病原数据库操作", "reference", request.view_args.get("host_key", "") if request.view_args else ""
    if "admin" in segments:
        if "users" in segments:
            return "后台管理", "用户管理操作", "user", request.view_args.get("username", "") if request.view_args else ""
        if "settings" in segments:
            return "后台管理", "脚本设置修改", "setting", ""
        if "update" in segments:
            return "后台管理", "系统更新", "update", ""
    if path in {"/login", "/logout"}:
        return "会话", "登录" if path == "/login" else "退出登录", "session", ""
    return "系统", "非查看操作", "", ""


def register_audit_hooks(app, store: PortalStore) -> None:
    @app.before_request
    def capture_audit_context() -> None:
        g.audit_actor = {
            "username": str(session.get("username") or "").strip(),
            "role": str(session.get("role") or "").strip(),
            "group_name": str(session.get("group_name") or "").strip(),
        }

    @app.after_request
    def write_audit_log(response: Response) -> Response:
        if request.method in {"GET", "HEAD"}:
            response.headers.update({"Cache-Control": "no-store, no-cache, must-revalidate, max-age=0", "Pragma": "no-cache", "Expires": "0"})
        if request.method in {"GET", "HEAD", "OPTIONS"} or request.path.startswith("/static/") or request.path == "/maketree":
            return response
        try:
            username, role, group_name = _identity()
            module, action, target_type, target_id = _meta(request.path, request.method)
            store.record_audit_log(
                {
                    "event_id": f"audit-{uuid.uuid4().hex}",
                    "username": username, "role": role, "group_name": group_name,
                    "method": request.method, "path": request.path, "endpoint": str(request.endpoint or ""),
                    "module": module, "action": action, "target_type": target_type, "target_id": target_id,
                    "status_code": int(response.status_code), "outcome": "success" if response.status_code < 400 else "failed",
                    "ip_address": str(request.headers.get("X-Forwarded-For", request.remote_addr or "")).split(",")[0].strip(),
                    "user_agent": str(request.user_agent.string or "").strip()[:255],
                    "request_summary": _request_summary(), "response_summary": _response_summary(response),
                    "created_at": datetime.utcnow().isoformat(),
                }
            )
        except Exception:
            app.logger.exception("audit logging failed")
        return response
