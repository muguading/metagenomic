from __future__ import annotations

from typing import Any


def _text(value: Any, default: str = "") -> str:
    text = str(value or "").strip()
    return text or default


def _first_metric_display(metrics: list[dict], key: str, default: str = "--") -> str:
    for item in metrics:
        if item.get("key") != key:
            continue
        if _text(item.get("display")):
            return _text(item.get("display"), default)
        nested = item.get("items") if isinstance(item.get("items"), list) else []
        for nested_item in nested:
            if _text(nested_item.get("display")):
                return _text(nested_item.get("display"), default)
    return default


def _dominant_species(payload: dict) -> str:
    metrics = payload.get("overview_metrics") if isinstance(payload.get("overview_metrics"), list) else []
    for item in metrics:
        if item.get("key") != "species_estimation":
            continue
        for nested in item.get("items") or []:
            display = _text(nested.get("display"))
            if display and display != "--":
                return display
    knowledge = ((payload.get("sections") or {}).get("knowledge_interpretation") or {})
    dominant = knowledge.get("dominant_species") if isinstance(knowledge.get("dominant_species"), dict) else {}
    return _text(dominant.get("species") or dominant.get("scientific_name"), "未明确")


def _resistance_count(sections: dict) -> int:
    rows = (((sections.get("resistance_virulence") or {}).get("resistance_elements") or {}).get("rows") or [])
    return len(rows) if isinstance(rows, list) else 0


def _virulence_count(sections: dict) -> int:
    rows = (((sections.get("resistance_virulence") or {}).get("virulence_elements") or {}).get("rows") or [])
    return len(rows) if isinstance(rows, list) else 0


def _risk_level(public_health: dict, resistance_count: int, virulence_count: int, is_meta: bool) -> tuple[str, str]:
    support_status = _text(public_health.get("status"))
    priority = _text(public_health.get("china_priority_group") or public_health.get("who_priority_group")).lower()
    if support_status == "matched" and any(token in priority for token in ("critical", "法定重点", "high", "医院感染重点")):
        return "高", "命中重点监测/优先病原，并存在公共卫生知识库支持。"
    if resistance_count >= 10 or virulence_count >= 10:
        return "高", "耐药或毒力条目数量较高，建议优先复核并纳入重点监测。"
    if support_status == "matched" or resistance_count >= 3 or virulence_count >= 3 or is_meta:
        return "中", "存在需要结合流调、历史样本和质控结果继续判读的监测信号。"
    return "低", "当前未见明确高优先级公共卫生信号，建议按常规监测归档。"


def _action(action_id: str, label: str, reason: str, href: str = "", priority: str = "recommended") -> dict:
    return {"id": action_id, "label": label, "reason": reason, "href": href, "priority": priority}


def _build_actions(task: dict, risk: str, public_health: dict, resistance_count: int, virulence_count: int) -> list[dict]:
    task_id = _text(task.get("id"))
    actions = [
        _action(
            "review_report",
            "复核报告结论",
            "先确认物种、分型、耐药毒力与质控是否支持当前研判。",
            f"/workstation?tab=queue&task={task_id}&closure_action=review_report" if task_id else "",
            "required",
        ),
        _action(
            "import_sample",
            "入库沉淀样本",
            "将本次结果纳入样本库，后续才能做历史对照、趋势跟踪和证据追溯。",
            f"/workstation?tab=queue&task={task_id}&closure_action=import_sample" if task_id else "/workstation?tab=queue",
            "required" if risk in {"高", "中"} else "recommended",
        ),
    ]
    if risk in {"高", "中"} or public_health.get("status") == "matched":
        actions.append(
            _action(
                "compare_history",
                "开展历史对照",
                "比对同地区、同型别或同物种历史样本，判断是否存在聚集或传播链线索。",
                f"/workstation?tab=queue&task={task_id}&closure_action=compare_history" if task_id else "/workstation?tab=queue",
                "required" if risk == "高" else "recommended",
            )
        )
    if risk == "高" or resistance_count >= 3 or virulence_count >= 3:
        actions.append(
            _action(
                "trace_source",
                "评估是否溯源",
                "耐药、毒力或重点病原信号较强时，应考虑 PathoSource/同源性分析。",
                "/workstation?module=pathosource",
                "recommended",
            )
        )
    actions.append(
        _action(
            "export_cdc_report",
            "导出疾控报告",
            "形成可用于会商、归档或上报前复核的标准化材料。",
            "",
            "recommended",
        )
    )
    return actions


def build_workflow_closure(payload: dict) -> dict:
    task = payload.get("task") if isinstance(payload.get("task"), dict) else {}
    sections = payload.get("sections") if isinstance(payload.get("sections"), dict) else {}
    metrics = payload.get("overview_metrics") if isinstance(payload.get("overview_metrics"), list) else []
    public_health = sections.get("public_health_support") if isinstance(sections.get("public_health_support"), dict) else {}
    resistance_count = _resistance_count(sections)
    virulence_count = _virulence_count(sections)
    is_meta = _text(task.get("method") or (task.get("params") or {}).get("method")).lower() == "meta"
    risk, risk_reason = _risk_level(public_health, resistance_count, virulence_count, is_meta)
    species_name = _dominant_species(payload)
    evidence = [
        {"label": "主要病原", "value": species_name},
        {"label": "总数据量/覆盖", "value": _first_metric_display(metrics, "total_bases")},
        {"label": "耐药条目", "value": str(resistance_count)},
        {"label": "毒力条目", "value": str(virulence_count)},
    ]
    if public_health.get("status") == "matched":
        evidence.append(
            {
                "label": "公共卫生知识库",
                "value": _text(
                    public_health.get("china_priority_group")
                    or public_health.get("who_priority_group")
                    or public_health.get("surveillance_domain"),
                    "已命中",
                ),
            }
        )
    return {
        "status": "ready",
        "title": "疾控处置闭环",
        "risk_level": risk,
        "risk_reason": risk_reason,
        "summary": f"当前样本主要指向 {species_name}；系统建议按“复核结论 → 入库沉淀 → 历史对照 → 报告导出”的顺序完成闭环。",
        "evidence": evidence,
        "actions": _build_actions(task, risk, public_health, resistance_count, virulence_count),
        "audit_hint": "上述动作均应保留操作者、时间、任务版本、参考库版本和导出记录，形成单样本证据链。",
    }
