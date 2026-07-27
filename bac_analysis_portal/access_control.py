from __future__ import annotations

from functools import wraps
import re

from flask import Response, current_app, jsonify, session

from .store import DEFAULT_ALLOWED_VIRUSES, PortalStore
from .task_manager import ValidationError

VIRUS_PERMISSION_ALIASES = {
    "ncov": ["sars-cov-2", "新型冠状病毒", "新冠病毒", "新冠", "2019-ncov"],
    "flu": ["influenza virus", "influenza a virus", "influenza b virus", "流感病毒", "甲型流感病毒", "乙型流感病毒", "甲流", "乙流"],
    "rsv": ["respiratory syncytial virus", "human respiratory syncytial virus", "orthopneumovirus hominis", "呼吸道合胞病毒", "rsv"],
    "hmpv": ["human metapneumovirus", "metapneumovirus", "人偏肺病毒", "hmpv"],
    "hpiv": ["human parainfluenza virus", "parainfluenza virus", "副流感病毒", "hpiv"],
    "hadv": ["adenovirus", "human adenovirus", "腺病毒"],
    "rhinovirus": ["human rhinovirus", "rhinovirus", "鼻病毒", "hrv"],
    "seasonal_hcov": ["human coronavirus", "coronavirus", "人冠状病毒", "冠状病毒"],
    "mpox": ["monkeypox virus", "mpox virus", "human monkeypox virus", "猴痘病毒", "猴痘", "mpox", "hmpxv"],
    "denv": ["dengue virus", "dengue", "denv", "登革热病毒", "登革热"],
    "zikav": ["zika virus", "zika", "zikv", "zikav", "寨卡病毒", "寨卡"],
    "chikv": ["chikungunya virus", "chikungunya", "chikv", "基孔肯雅病毒", "基孔肯雅"],
    "bandavirus": ["bandavirus dabieense", "bandavirus", "sftsv", "大别班达病毒", "班达病毒", "发热伴血小板减少综合征病毒"],
    "orthohantavirus": ["orthohantavirus", "hantavirus", "汉坦病毒", "汉他病毒"],
    "orthoebolavirus": ["orthoebolavirus", "ebola virus", "ebolavirus", "ebov", "埃博拉病毒", "正埃博拉病毒"],
    "norovirus": ["norovirus", "诺如病毒"], "rotavirus": ["rotavirus a", "rotavirus", "a组轮状病毒", "轮状病毒"],
    "astroviridae": ["astrovirus", "星状病毒"], "enterovirus": ["human enterovirus", "enterovirus", "coxsackievirus", "echovirus", "肠道病毒", "柯萨奇病毒", "埃可病毒"],
    "hepatovirus": ["hepatovirus", "hepatitis a virus", "hav", "甲肝病毒", "甲型肝炎病毒", "甲型肝炎"],
    "hiv": ["hiv", "hiv-1", "hiv1", "human immunodeficiency virus", "human immunodeficiency virus 1", "艾滋病病毒", "艾滋病毒"],
    "sapovirus": ["sapovirus", "札如病毒"],
}
DEMO_TYPE_TO_VIRUS_PERMISSION = {
    key: key for key in (
        "ncov", "flu", "rsv", "hmpv", "hpiv", "hadv", "norovirus", "rotavirus", "astroviridae",
        "bandavirus", "denv", "zikav", "chikv", "rhinovirus", "enterovirus", "hepatovirus", "hiv",
        "orthohantavirus", "orthoebolavirus", "seasonal_hcov",
    )
}
DEMO_TYPE_TO_VIRUS_PERMISSION["hmpxv"] = "mpox"


def _normalize_virus_permission_lookup(value: object) -> str:
    return re.sub(r"[\s,_/]+", " ", str(value or "").strip().lower().replace("（", "(").replace("）", ")")).strip()


def resolve_requested_virus_permission(species_name: object) -> str:
    normalized = _normalize_virus_permission_lookup(species_name)
    for key, aliases in VIRUS_PERMISSION_ALIASES.items():
        if any(normalized == _normalize_virus_permission_lookup(alias) or _normalize_virus_permission_lookup(alias) in normalized for alias in aliases):
            return key
    return ""


def ensure_user_has_virus_permission(user: dict[str, object], virus_key: str) -> None:
    if not virus_key:
        return
    allowed = user.get("allowed_viruses") if isinstance(user, dict) else None
    normalized = {str(item or "").strip().lower() for item in (allowed if isinstance(allowed, list) else DEFAULT_ALLOWED_VIRUSES) if str(item or "").strip()}
    if virus_key not in normalized:
        raise ValidationError("当前账号未开通该病毒类型权限")


def login_required(view_func):
    @wraps(view_func)
    def wrapped(*args, **kwargs):
        if not is_logged_in():
            return jsonify({"error": "Authentication required"}), 401
        return view_func(*args, **kwargs)

    return wrapped


def admin_required(view_func):
    @wraps(view_func)
    def wrapped(*args, **kwargs):
        if not is_logged_in():
            return jsonify({"error": "Authentication required"}), 401
        if session.get("role") != "admin":
            return jsonify({"error": "Administrator permission required"}), 403
        return view_func(*args, **kwargs)

    return wrapped


def is_logged_in() -> bool:
    username = str(session.get("username") or "").strip()
    if not (username and session.get("role")):
        return False
    try:
        user = get_store().get_user(username)
    except KeyError:
        session.clear()
        return False
    if user.get("is_expired"):
        session.clear()
        return False
    return True


def current_username() -> str:
    return str(session.get("username") or "").strip()


def get_store() -> PortalStore:
    store = current_app.config.get("PORTAL_STORE")
    if not isinstance(store, PortalStore):
        raise RuntimeError("Portal store is not configured")
    return store


def can_view_task(task: dict) -> bool:
    role = str(session.get("role") or "")
    group_name = str(session.get("group_name") or "")
    username = str(session.get("username") or "")
    if role == "admin":
        return True
    if username and str(task.get("owner") or "") == username:
        return True
    if group_name and str(task.get("owner_group") or "") == group_name:
        return True
    return False


def ensure_can_view_task(task: dict) -> None:
    if not can_view_task(task):
        raise KeyError(f"Task not found: {task.get('id', '-')}")


def ensure_can_modify_task(task: dict) -> None:
    role = str(session.get("role") or "")
    group_name = str(session.get("group_name") or "")
    if role == "admin":
        return
    if role == "group_admin" and group_name and str(task.get("owner_group") or "") == group_name:
        return
    raise ValidationError("只有管理员或 group 管理可以删除组内任务")


def ensure_can_control_task(task: dict) -> None:
    role = str(session.get("role") or "")
    group_name = str(session.get("group_name") or "")
    username = str(session.get("username") or "")
    if role == "admin":
        return
    if role == "group_admin" and group_name and str(task.get("owner_group") or "") == group_name:
        return
    if username and username == str(task.get("owner") or ""):
        return
    raise ValidationError("只有任务本人、管理员或 group 管理可以控制任务")


def can_view_uploaded_auspice(metadata: dict) -> bool:
    return (
        str(session.get("role") or "") == "admin"
        or bool(str(session.get("username") or "") and str(session.get("username") or "") == str(metadata.get("owner") or "").strip())
        or bool(str(session.get("group_name") or "") and str(session.get("group_name") or "") == str(metadata.get("owner_group") or "").strip())
    )


def ensure_can_view_uploaded_auspice(metadata: dict) -> None:
    if not can_view_uploaded_auspice(metadata):
        raise KeyError(f"Uploaded auspice dataset not found: {metadata.get('id', '-')}")


def add_no_cache_headers(response: Response) -> Response:
    response.headers["Cache-Control"] = "no-store, no-cache, must-revalidate, max-age=0"
    response.headers["Pragma"] = "no-cache"
    response.headers["Expires"] = "0"
    return response
