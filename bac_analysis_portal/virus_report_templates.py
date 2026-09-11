from __future__ import annotations

import copy
import json
import uuid
from datetime import datetime

from flask import current_app

from .store import PortalStore
from .task_manager import ValidationError

VIRUS_REPORT_TEMPLATE_SETTING_KEY = "virus_report_templates"
VIRUS_REPORT_TEMPLATE_HISTORY_SETTING_KEY = "virus_report_template_history"
DEFAULT_VIRUS_REPORT_TEMPLATES: list[dict[str, object]] = [
    {
        "id": "virus_influenza",
        "name": "流感病毒",
        "category": "respiratory",
        "modes": ["influenza_typing"],
        "clinicalRisk": "中风险",
        "cdcRisk": "中",
        "evidenceBasis": "流感类型、HA/NA 分型、segment 组成、覆盖度与突变注释",
        "clinicalMeaning": "流感病毒检出具有明确呼吸道感染相关性，但病情轻重、传染期和治疗决策仍需结合症状、采样时间、抗原/核酸复核及基础疾病综合判断。",
        "cdcMeaning": "流感结果应纳入呼吸道传染病季节性监测，重点关注亚型变化、聚集性病例和重症病例比例。",
        "clinicalRecommendations": ["建议结合发病时间窗、重症风险因素和当地诊疗规范评估抗病毒治疗时机。", "建议结合呼吸道症状、抗原/核酸复核结果和采样质量判断检出结果的临床相关性。", "如为重症、聚集性或特殊人群样本，建议优先复核分型、segment 覆盖度和关键突变位点。"],
        "cdcRecommendations": ["建议纳入流感季节性监测，关注亚型变化、重症比例和聚集性病例。", "建议与本地哨点监测、疫苗株背景和历史序列进行比对，判断是否存在异常谱系变化。", "如出现学校、养老机构或医院聚集性病例，建议补充分型/同源性证据并按规范处置。"],
    },
    {
        "id": "virus_monkeypox",
        "name": "猴痘病毒",
        "category": "contact-transmitted",
        "modes": ["monkeypox_nextclade"],
        "clinicalRisk": "中风险",
        "cdcRisk": "中",
        "evidenceBasis": "Nextclade clade/lineage、覆盖度、QC 与突变位点",
        "clinicalMeaning": "猴痘病毒检出需结合皮疹、发热、暴露史与采样部位判断临床相关性；分型结果可辅助追踪传播背景，但不替代临床诊断流程。",
        "cdcMeaning": "猴痘结果具有公共卫生追踪价值，建议结合个案调查、接触者管理和本地历史序列开展传播链评估。",
        "clinicalRecommendations": ["建议结合皮疹部位、病程阶段、暴露史和采样类型评估检出结果的临床意义。", "建议必要时复核关键位点和覆盖度，避免因低覆盖或混样造成谱系解释偏差。", "临床处置仍应依据现行诊疗规范和感染防控要求执行。"],
        "cdcRecommendations": ["建议结合个案调查、接触者管理和活动轨迹评估传播链。", "建议将 clade/lineage 与本地及上级平台历史序列比对，判断是否为既有传播链延续。", "如存在聚集性或跨区域关联，应补充同源性分析和暴露网络信息。"],
    },
    {
        "id": "virus_hepatitis",
        "name": "肝炎病毒",
        "category": "hepatitis",
        "modes": ["hepatovirus_typing"],
        "clinicalRisk": "中风险",
        "cdcRisk": "中",
        "evidenceBasis": "肝炎病毒 broad 大亚型、子亚型/基因型参考竞争、覆盖度与突变注释",
        "clinicalMeaning": "肝炎病毒分型结果可辅助判断病毒类别与分子流行病学背景；临床诊断仍需结合肝功能、血清学标志物、病毒载量、病程阶段和既往感染/免疫史综合判断。",
        "cdcMeaning": "肝炎病毒结果应重点关注传播途径、感染来源和同型别聚集情况；不同大亚型对应的流调问题不同，不应只按通用病毒阳性处理。",
        "clinicalRecommendations": ["建议结合肝功能、血清学标志物、病毒载量和临床病程判断活动性感染与疾病阶段。", "建议核对 broad 大亚型与子亚型/基因型是否一致，必要时复核覆盖度和关键参考株选择。", "如涉及治疗或随访，应由临床结合指南、既往感染史和免疫状态综合决策。"],
        "cdcRecommendations": ["建议按 HAV/HBV/HCV/HDV/HEV 不同传播特点分别补充流调信息，不要混用同一处置路径。", "建议关注同一单位、家庭或共同暴露场景中的同型别聚集信号。", "如用于暴发研判，应补充采样时间、暴露史和更高分辨率的序列比较证据。"],
    },
    {
        "id": "virus_natural_focus",
        "name": "自然疫源性病毒",
        "category": "natural-focus",
        "modes": ["bandavirus_typing", "orthohantavirus_typing", "orthoebolavirus_typing", "ebola_nextclade"],
        "clinicalRisk": "中风险",
        "cdcRisk": "中",
        "evidenceBasis": "自然疫源性病毒分型、分段参考选择、覆盖度与潜在重配提示",
        "clinicalMeaning": "该类病毒结果应结合发热、出血倾向、肾损伤或血小板变化等临床表现判断；分段分型结果可辅助判断自然疫源背景和潜在暴露来源。",
        "cdcMeaning": "该类病毒具备自然疫源性监测意义，建议结合病例暴露史、地域来源、媒介/宿主线索和本地监测资料综合研判。",
        "clinicalRecommendations": ["建议结合发热、血小板、肾功能、出血倾向及流行病学暴露史综合判断临床相关性。", "建议复核分段分型是否一致，若多片段证据不一致，应谨慎解释潜在重配或混合信号。", "如病情进展或暴露史明确，建议结合规范检测和临床专科意见动态评估。"],
        "cdcRecommendations": ["建议补充病例居住地、活动地、野外或农田暴露、动物或媒介接触信息。", "建议结合宿主和媒介监测资料判断是否存在自然疫源地活跃信号。", "如同一区域出现多例同型别结果，建议开展时空聚集和传播风险复核。"],
    },
    {
        "id": "virus_vector_borne",
        "name": "虫媒病毒",
        "category": "vector-borne",
        "modes": ["denv_nextclade", "zikav_nextclade", "chikv_nextclade"],
        "clinicalRisk": "低-中风险",
        "cdcRisk": "中",
        "evidenceBasis": "虫媒病毒 clade/lineage、覆盖度、QC 与突变位点",
        "clinicalMeaning": "虫媒病毒检出需结合发热、皮疹、关节痛、出血表现、旅行史和采样时间窗判断临床相关性；分型结果更适合用于输入来源和传播背景分析。",
        "cdcMeaning": "该结果具备虫媒病毒监测意义，建议结合旅行史、媒介密度、病例时空分布和本地输入/本地传播背景研判。",
        "clinicalRecommendations": ["建议结合旅行史、蚊媒暴露史、发病时间窗和血清学/核酸复核结果判断临床意义。", "建议关注采样时间对核酸检出率的影响，必要时补充血清学或复采证据。", "如出现重症表现，应结合当地诊疗规范和实验室确认结果及时评估。"],
        "cdcRecommendations": ["建议核查旅行史、活动轨迹和发病地，区分输入病例与本地传播风险。", "建议结合媒介密度、季节和周边病例分布判断是否需要强化媒介控制。", "如出现同区域同时间多例，应补充分型/系统发育比较和现场流调证据。"],
    },
    {
        "id": "virus_bloodborne",
        "name": "血源性病毒",
        "category": "bloodborne",
        "modes": ["hiv_resistance"],
        "clinicalRisk": "中风险",
        "cdcRisk": "中",
        "evidenceBasis": "血源性病毒分型、耐药或关键位点、覆盖度与参考选择结果",
        "clinicalMeaning": "血源性或慢性感染相关病毒结果应结合确认试验、病毒载量、免疫状态和既往诊疗史综合解释；测序分型可辅助耐药、传播背景和随访管理。",
        "cdcMeaning": "该结果可用于传播网络和重点人群监测线索，但不应替代确认试验、个案管理和规范报告流程。",
        "clinicalRecommendations": ["建议结合确认试验、病毒载量、免疫状态和既往治疗史判断临床意义。", "如涉及耐药或亚型解释，应复核覆盖度、关键位点和参考选择结果。", "诊疗决策应依据现行临床指南和专科评估，不以单一测序报告替代。"],
        "cdcRecommendations": ["建议结合个案管理、传播风险评估和重点人群监测资料进行解释。", "如用于传播网络分析，应补充匿名化流调信息和更高分辨率序列比较。", "涉及报告管理时，应按现行规范和确认试验结果执行。"],
    },
    {
        "id": "virus_respiratory",
        "name": "呼吸道病毒",
        "category": "respiratory",
        "modes": ["rsv_nextclade", "hmpv_nextclade", "hpiv_typing", "hadv_typing", "rhinovirus_typing", "seasonal_hcov_typing"],
        "clinicalRisk": "低-中风险",
        "cdcRisk": "低",
        "evidenceBasis": "呼吸道病毒分型、覆盖度、QC 与突变位点",
        "clinicalMeaning": "呼吸道病毒检出需结合症状、采样部位、病程阶段和共感染背景判断临床意义；分型结果可辅助判断流行株背景和院内/社区传播线索。",
        "cdcMeaning": "呼吸道病毒结果适合纳入季节性和聚集性监测，重点关注同型别病例聚集、特殊机构暴发和重症病例。",
        "clinicalRecommendations": ["建议结合呼吸道症状、病程阶段、采样质量和共感染证据判断临床相关性。", "如为重症或免疫低下患者，建议复核覆盖度和关键变异位点。", "抗病毒或感染控制决策应结合当地诊疗规范、病原确认和患者风险因素。"],
        "cdcRecommendations": ["建议关注同型别病例在学校、养老机构、医院等场景中的聚集信号。", "建议与同期本地呼吸道病原监测数据对照，判断是否存在流行株变化。", "如出现异常重症或聚集性事件，应补充分型和同源性证据。"],
    },
    {
        "id": "virus_enteric",
        "name": "肠道/胃肠炎相关病毒",
        "category": "enteric",
        "modes": ["norovirus_typing", "enterovirus_typing", "rotavirus_typing", "astroviridae_typing"],
        "clinicalRisk": "低-中风险",
        "cdcRisk": "低",
        "evidenceBasis": "肠道病毒或胃肠炎相关病毒分型、覆盖度与参考竞争结果",
        "clinicalMeaning": "肠道病毒或胃肠炎相关病毒检出需结合腹泻、呕吐、发热、采样时间和暴露史判断临床相关性；分型结果可辅助食品、水源或机构聚集事件研判。",
        "cdcMeaning": "该结果适合纳入肠道传染病或胃肠炎聚集性监测，重点关注共同暴露、机构传播和同型别聚集。",
        "clinicalRecommendations": ["建议结合胃肠道症状、采样时间、脱水程度和共感染结果判断临床意义。", "如结果用于个案诊疗，需结合病程和其他病原检测排除偶然携带或残留核酸。", "建议复核分型和覆盖度，尤其是用于聚集性事件解释时。"],
        "cdcRecommendations": ["建议补充共同就餐、水源、托幼/学校/养老机构暴露史。", "建议关注同型别病例聚集，并结合环境或食品样本结果进行综合判断。", "如涉及暴发调查，应补充采样时间轴和序列同源性比较。"],
    },
]


def _virus_report_template_category(serotype_result: dict) -> str:
    mode = str(serotype_result.get("mode") or "").strip()
    if mode == "influenza_typing":
        return "respiratory"
    if mode == "monkeypox_nextclade":
        return "contact-transmitted"
    if mode == "hepatovirus_typing":
        return "hepatitis"
    if mode in {"bandavirus_typing", "orthohantavirus_typing"}:
        return "natural-focus"
    if mode in {"denv_nextclade", "zikav_nextclade", "chikv_nextclade"}:
        return "vector-borne"
    if mode == "hiv_resistance":
        return "bloodborne"
    if mode in {"rsv_nextclade", "hmpv_nextclade", "hpiv_typing", "hadv_typing", "rhinovirus_typing", "seasonal_hcov_typing"}:
        return "respiratory"
    if mode in {"norovirus_typing", "enterovirus_typing", "rotavirus_typing", "astroviridae_typing"}:
        return "enteric"
    return "general"

def _entry_matches_virus_report_mode(entry: dict, mode: str) -> bool:
    match = entry.get("match") if isinstance(entry.get("match"), dict) else {}
    modes = match.get("modes")
    if not isinstance(modes, list):
        return False
    return mode in {str(item or "").strip() for item in modes}

def _normalize_report_template_list(value: object) -> list[str]:
    if not isinstance(value, list):
        return []
    return [str(item or "").strip() for item in value if str(item or "").strip()]

def _sanitize_virus_report_template_payload(payload: dict, base: dict | None = None) -> dict:
    base = base or {}
    sanitized = dict(base)
    for key in ["name", "category", "clinicalRisk", "cdcRisk", "evidenceBasis", "clinicalMeaning", "cdcMeaning"]:
        value = str(payload.get(key, sanitized.get(key, "")) or "").strip()
        if key in {"clinicalMeaning", "cdcMeaning"} and len(value) > 1200:
            raise ValidationError("报告意义文本过长，请控制在 1200 字以内")
        sanitized[key] = value
    for key in ["clinicalRecommendations", "cdcRecommendations"]:
        source = payload.get(key, sanitized.get(key, []))
        items = _normalize_report_template_list(source)
        if len(items) > 8:
            raise ValidationError("建议条目最多保留 8 条")
        if any(len(item) > 300 for item in items):
            raise ValidationError("单条建议请控制在 300 字以内")
        sanitized[key] = items
    sanitized["id"] = str(base.get("id") or payload.get("id") or "").strip()
    sanitized["modes"] = _normalize_report_template_list(base.get("modes") or payload.get("modes"))
    sanitized["source"] = str(payload.get("source") or base.get("source") or "default").strip()
    return sanitized

def _load_virus_report_template_overrides(store: PortalStore | None = None) -> dict[str, dict]:
    active_store = store
    if active_store is None:
        try:
            active_store = current_app.config.get("PORTAL_STORE")
        except RuntimeError:
            active_store = None
    if active_store is None:
        return {}
    raw_value = str(active_store.get_setting(VIRUS_REPORT_TEMPLATE_SETTING_KEY, "[]") or "[]").strip()
    try:
        payload = json.loads(raw_value)
    except json.JSONDecodeError:
        return {}
    if not isinstance(payload, list):
        return {}
    overrides: dict[str, dict] = {}
    for item in payload:
        if not isinstance(item, dict):
            continue
        template_id = str(item.get("id") or "").strip()
        if template_id:
            overrides[template_id] = item
    return overrides

def _list_virus_report_templates(store: PortalStore | None = None) -> list[dict]:
    overrides = _load_virus_report_template_overrides(store)
    history_map = _load_virus_report_template_history(store)
    templates: list[dict] = []
    for base in DEFAULT_VIRUS_REPORT_TEMPLATES:
        template_id = str(base.get("id") or "").strip()
        merged = dict(base)
        if isinstance(overrides.get(template_id), dict):
            merged.update(overrides[template_id])
            merged["source"] = "custom"
        else:
            merged["source"] = "default"
        template = _sanitize_virus_report_template_payload(merged, base)
        template["history"] = history_map.get(template_id, [])[:8]
        templates.append(template)
    return templates

def _load_virus_report_template_history(store: PortalStore | None = None) -> dict[str, list[dict]]:
    active_store = store
    if active_store is None:
        try:
            active_store = current_app.config.get("PORTAL_STORE")
        except RuntimeError:
            active_store = None
    if active_store is None:
        return {}
    raw_value = str(active_store.get_setting(VIRUS_REPORT_TEMPLATE_HISTORY_SETTING_KEY, "{}") or "{}").strip()
    try:
        payload = json.loads(raw_value)
    except json.JSONDecodeError:
        return {}
    if not isinstance(payload, dict):
        return {}
    history: dict[str, list[dict]] = {}
    for template_id, records in payload.items():
        if not isinstance(template_id, str) or not isinstance(records, list):
            continue
        cleaned: list[dict] = []
        for record in records:
            if isinstance(record, dict):
                cleaned.append(record)
        history[template_id] = cleaned[:20]
    return history

def _record_virus_report_template_history(store: PortalStore, template_id: str, action: str, snapshot: dict, actor: str = "") -> None:
    history = _load_virus_report_template_history(store)
    records = history.get(template_id, [])
    event = {
        "eventId": f"vrt_{datetime.now().strftime('%Y%m%d%H%M%S')}_{uuid.uuid4().hex[:8]}",
        "action": action,
        "actor": str(actor or "").strip() or "system",
        "createdAt": datetime.now().isoformat(timespec="seconds"),
        "snapshot": _sanitize_virus_report_template_payload(snapshot, snapshot),
    }
    history[template_id] = [event] + records
    history[template_id] = history[template_id][:20]
    store.set_setting(VIRUS_REPORT_TEMPLATE_HISTORY_SETTING_KEY, json.dumps(history, ensure_ascii=False, indent=2))

def _save_virus_report_template_override(store: PortalStore, template_id: str, payload: dict, actor: str = "", history_action: str = "save") -> dict:
    base = next((item for item in DEFAULT_VIRUS_REPORT_TEMPLATES if str(item.get("id") or "").strip() == template_id), None)
    if not isinstance(base, dict):
        raise ValidationError("报告模板不存在")
    sanitized = _sanitize_virus_report_template_payload(payload, base)
    sanitized["source"] = "custom"
    overrides = _load_virus_report_template_overrides(store)
    overrides[template_id] = sanitized
    ordered = [overrides[str(item.get("id") or "").strip()] for item in DEFAULT_VIRUS_REPORT_TEMPLATES if str(item.get("id") or "").strip() in overrides]
    store.set_setting(VIRUS_REPORT_TEMPLATE_SETTING_KEY, json.dumps(ordered, ensure_ascii=False, indent=2))
    if history_action:
        _record_virus_report_template_history(store, template_id, history_action, sanitized, actor)
    return sanitized

def _reset_virus_report_template_override(store: PortalStore, template_id: str, actor: str = "") -> dict:
    base = next((item for item in DEFAULT_VIRUS_REPORT_TEMPLATES if str(item.get("id") or "").strip() == template_id), None)
    if not isinstance(base, dict):
        raise ValidationError("报告模板不存在")
    overrides = _load_virus_report_template_overrides(store)
    overrides.pop(template_id, None)
    ordered = [overrides[str(item.get("id") or "").strip()] for item in DEFAULT_VIRUS_REPORT_TEMPLATES if str(item.get("id") or "").strip() in overrides]
    store.set_setting(VIRUS_REPORT_TEMPLATE_SETTING_KEY, json.dumps(ordered, ensure_ascii=False, indent=2))
    sanitized = _sanitize_virus_report_template_payload(dict(base, source="default"), base)
    _record_virus_report_template_history(store, template_id, "reset", sanitized, actor)
    return sanitized

def _rollback_virus_report_template_override(store: PortalStore, template_id: str, event_id: str, actor: str = "") -> dict:
    history = _load_virus_report_template_history(store)
    record = next((item for item in history.get(template_id, []) if str(item.get("eventId") or "") == event_id), None)
    if not isinstance(record, dict) or not isinstance(record.get("snapshot"), dict):
        raise ValidationError("未找到可回滚的模板版本")
    snapshot = dict(record["snapshot"])
    snapshot["source"] = "custom"
    item = _save_virus_report_template_override(store, template_id, snapshot, actor, "rollback")
    return item

def _build_virus_report_template_summary(project_root_text: str, serotype_result: dict) -> dict:
    if not isinstance(serotype_result, dict):
        return {"status": "empty"}
    entries = _list_virus_report_templates()
    if not entries:
        return {"status": "empty"}
    mode = str(serotype_result.get("mode") or "").strip()
    category = _virus_report_template_category(serotype_result)
    selected = next((entry for entry in entries if isinstance(entry, dict) and mode in set(entry.get("modes") or [])), None)
    if selected is None:
        selected = next((entry for entry in entries if isinstance(entry, dict) and str(entry.get("category") or "").strip() == category), None)
    if not isinstance(selected, dict):
        return {"status": "empty"}
    return {
        "status": "ready",
        "source": str(selected.get("source") or "default").strip(),
        "id": str(selected.get("id") or "").strip(),
        "category": str(selected.get("category") or category or "general").strip(),
        "clinicalRisk": str(selected.get("clinicalRisk") or "").strip(),
        "cdcRisk": str(selected.get("cdcRisk") or "").strip(),
        "evidenceBasis": str(selected.get("evidenceBasis") or "").strip(),
        "clinicalMeaning": str(selected.get("clinicalMeaning") or "").strip(),
        "cdcMeaning": str(selected.get("cdcMeaning") or "").strip(),
        "clinicalRecommendations": _normalize_report_template_list(selected.get("clinicalRecommendations")),
        "cdcRecommendations": _normalize_report_template_list(selected.get("cdcRecommendations")),
    }
