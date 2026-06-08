from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any, Callable
from uuid import uuid4

from .store import PortalStore


ReportSourceResolver = Callable[[dict, str], dict[str, Any]]
ReportPayloadBuilder = Callable[[dict, str], dict[str, Any]]
SampleNameResolver = Callable[[dict, Path], str]
HumanBpFormatter = Callable[[object], str]


DEFAULT_SAMPLE_METADATA_TEMPLATES: list[dict[str, Any]] = [
    {
        "key": "case_id",
        "label": "病例/事件编号",
        "type": "text",
        "group": "case",
        "requirement": "required",
        "placeholder": "例如 HP-CDC-2026-001",
        "help_text": "建议使用统一的病例或事件编号，便于跨样本追踪。",
        "dictionary_name": "病例编号规则",
    },
    {
        "key": "patient_id",
        "label": "患者/个案编号",
        "type": "text",
        "group": "case",
        "requirement": "recommended",
        "placeholder": "例如 PID-240301",
        "help_text": "个案编号可与 LIS、住院号或内部个案台账保持一致。",
        "dictionary_name": "个案编号规则",
    },
    {
        "key": "surveillance_source",
        "label": "监测来源",
        "type": "select",
        "options": ["门急诊", "住院", "ICU", "发热门诊", "社区监测", "环境监测", "食品监测", "污水监测", "动物源", "其他"],
        "group": "epi",
        "requirement": "required",
        "help_text": "优先使用标准来源词表，避免出现“住院部”“病房”等自由文本混写。",
        "dictionary_name": "监测来源标准词表",
    },
    {
        "key": "suspected_syndrome",
        "label": "疑似症候群",
        "type": "select",
        "options": ["脑膜炎/脑膜脑炎", "肺炎/呼吸道感染", "败血症/血流感染", "腹泻/肠道感染", "皮肤软组织感染", "泌尿生殖道感染", "医院感染监测", "环境异常事件", "不明原因感染", "其他"],
        "group": "epi",
        "requirement": "required",
        "help_text": "按监测专题选择最接近的症候群，便于后续统计和专题报告。",
        "dictionary_name": "症候群标准词表",
    },
    {
        "key": "submitting_unit",
        "label": "送检单位",
        "type": "text",
        "options": ["黄浦区疾控中心", "区中心医院", "瑞金医院黄浦分院", "社区卫生服务中心", "学校卫生室", "第三方检测机构"],
        "group": "case",
        "requirement": "recommended",
        "placeholder": "例如 黄浦区疾控中心",
        "help_text": "填写报告或标本送检来源单位，便于回溯。",
        "dictionary_name": "送检单位标准词表",
    },
    {
        "key": "collection_unit",
        "label": "采样单位",
        "type": "text",
        "options": ["黄浦区疾控采样组", "瑞金医院黄浦分院", "社区采样点", "学校采样点", "环境监测点位", "病区采样小组"],
        "group": "case",
        "requirement": "recommended",
        "placeholder": "例如 瑞金医院黄浦分院",
        "help_text": "若与送检单位不同，建议单独记录。",
        "dictionary_name": "采样单位标准词表",
    },
    {
        "key": "ward_department",
        "label": "科室/病区",
        "type": "text",
        "options": ["ICU", "急诊科", "发热门诊", "呼吸内科", "感染科", "神经内科", "儿科", "环境采样区"],
        "group": "sampling",
        "requirement": "recommended",
        "placeholder": "例如 ICU / 神经内科",
        "help_text": "院感监测样本建议补齐病区或重点科室。",
        "dictionary_name": "科室病区标准词表",
    },
    {
        "key": "specimen_category",
        "label": "标本类别",
        "type": "select",
        "options": ["血液", "脑脊液", "呼吸道", "粪便/肛拭子", "尿液", "脓液/分泌物", "组织/活检", "环境拭子", "污水/水体", "食品", "其他"],
        "group": "sampling",
        "requirement": "required",
        "help_text": "建议统一使用标准标本类别，避免同义写法造成统计分散。",
        "dictionary_name": "标本类别标准词表",
    },
    {
        "key": "collection_site",
        "label": "采样地点",
        "type": "location",
        "group": "sampling",
        "requirement": "required",
        "options": ["门诊采样点", "急诊采样点", "ICU监测点", "病区治疗室", "污水处理点", "环境高频接触面", "学校食堂", "社区采样点"],
        "help_text": "中国样本建议补齐省/市/区县，detail 可写医院、街道、点位名称。",
        "dictionary_name": "采样地点详细点位词表",
    },
    {
        "key": "cluster_status",
        "label": "聚集性状态",
        "type": "select",
        "options": ["散发", "聚集性", "暴发相关", "待判定"],
        "group": "epi",
        "requirement": "recommended",
        "help_text": "用于区分常规散发样本与聚集/暴发相关样本。",
        "dictionary_name": "聚集性状态词表",
    },
    {
        "key": "epidemiology_link",
        "label": "流行病学关联",
        "type": "text",
        "group": "epi",
        "requirement": "recommended",
        "placeholder": "例如 与 2026 春季 ICU 事件关联",
        "help_text": "可记录接触史、共同暴露、事件编号等关键信息。",
    },
    {
        "key": "traditional_result",
        "label": "传统检测结果",
        "type": "text",
        "group": "lab",
        "requirement": "recommended",
        "placeholder": "例如 培养阳性，PCR 阳性，药敏待回报",
        "help_text": "可填写培养、PCR、抗原、血清学或药敏摘要。",
    },
]


def _parse_metadata_items(raw_value: object) -> list[dict[str, Any]]:
    text = str(raw_value or "").strip()
    if not text:
        return []
    try:
        items = json.loads(text)
    except Exception:
        return []
    return [item for item in items if isinstance(item, dict)] if isinstance(items, list) else []


def _extract_mlst_st_from_payload(payload: object) -> str:
    sections = payload.get("sections", {}) if isinstance(payload, dict) else {}
    mlst_section = sections.get("mlst", {}) if isinstance(sections.get("mlst"), dict) else {}
    columns = mlst_section.get("columns") if isinstance(mlst_section.get("columns"), list) else []
    rows = mlst_section.get("rows") if isinstance(mlst_section.get("rows"), list) else []
    if not columns or not rows:
        return ""
    target_index = next(
        (idx for idx, label in enumerate(columns) if str(label or "").strip() in {"序列分型(ST)", "序列分型", "ST"}),
        -1,
    )
    if target_index < 0:
        return ""
    seen: list[str] = []
    for row in rows:
        if not isinstance(row, list) or target_index >= len(row):
            continue
        value = str(row[target_index] or "").strip()
        if not value or value == "-":
            continue
        normalized = value if value.upper().startswith("ST") else f"ST{value}"
        if normalized not in seen:
            seen.append(normalized)
    return " / ".join(seen[:3])


def _clean_typing_token(value: object) -> str:
    text = str(value or "").strip()
    if not text or text in {"-", "--", "/", "-/-", "- / -"}:
        return ""
    return text


def _append_unique_typing(parts: list[str], value: object, *, prefix: str = "") -> None:
    cleaned = _clean_typing_token(value)
    if not cleaned:
        return
    display = f"{prefix}{cleaned}" if prefix else cleaned
    normalized = display.lower().replace(" ", "")
    existing = [item.lower().replace(" ", "") for item in parts]
    if normalized in existing:
        return
    if not prefix and len(cleaned) <= 3:
        for item in existing:
            if item.endswith(f"-{normalized}") or item.endswith(f"/{normalized}"):
                return
    parts.append(display)


def _summary_card_value(section: dict[str, Any], label: str) -> str:
    summary_cards = section.get("summary_cards") if isinstance(section.get("summary_cards"), list) else []
    for item in summary_cards:
        if not isinstance(item, dict):
            continue
        if str(item.get("label") or "").strip() != label:
            continue
        return _clean_typing_token(item.get("value"))
    return ""


def _extract_structured_serotype_from_section(section: dict[str, Any]) -> str:
    mode = str(section.get("mode") or "").strip()
    parts: list[str] = []
    if mode == "influenza_typing":
        _append_unique_typing(parts, section.get("predicted_serotype"))
        if not parts:
            influenza_type = _clean_typing_token(section.get("influenza_type"))
            ha_subtype = _clean_typing_token(section.get("ha_subtype"))
            na_subtype = _clean_typing_token(section.get("na_subtype"))
            if influenza_type and ha_subtype and na_subtype:
                parts.append(f"{influenza_type}({ha_subtype}{na_subtype})")
    elif mode in {"rsv_nextclade", "hmpv_nextclade", "denv_nextclade", "zikav_nextclade", "chikv_nextclade"}:
        _append_unique_typing(parts, section.get("predicted_lineage"))
        _append_unique_typing(parts, section.get("predicted_serotype"))
        _append_unique_typing(parts, section.get("predicted_clade"))
    elif mode == "monkeypox_nextclade":
        _append_unique_typing(parts, section.get("predicted_lineage"))
        _append_unique_typing(parts, section.get("predicted_clade"))
        _append_unique_typing(parts, section.get("predicted_outbreak"))
    elif mode == "hadv_typing":
        _append_unique_typing(parts, section.get("predicted_clade"))
        _append_unique_typing(parts, _summary_card_value(section, "Penton"), prefix="Penton:")
        _append_unique_typing(parts, _summary_card_value(section, "Hexon"), prefix="Hexon:")
        _append_unique_typing(parts, _summary_card_value(section, "Fiber"), prefix="Fiber:")
    elif mode == "norovirus_typing":
        _append_unique_typing(parts, section.get("predicted_clade"))
    elif mode == "enterovirus_typing":
        _append_unique_typing(parts, section.get("predicted_clade"))
        _append_unique_typing(parts, section.get("predicted_group"))
    elif mode == "hepatovirus_typing":
        _append_unique_typing(parts, section.get("predicted_group"))
        _append_unique_typing(parts, section.get("predicted_subtype"))
        _append_unique_typing(parts, section.get("predicted_clade"))
    elif mode == "bandavirus_typing":
        _append_unique_typing(parts, section.get("predicted_lineage"))
        _append_unique_typing(parts, section.get("predicted_clade"))
        _append_unique_typing(parts, section.get("predicted_group"))
    elif mode == "orthohantavirus_typing":
        _append_unique_typing(parts, section.get("predicted_clade"))
        _append_unique_typing(parts, section.get("predicted_group"))
    elif mode == "ebola_nextclade":
        _append_unique_typing(parts, section.get("predicted_lineage"))
        _append_unique_typing(parts, section.get("predicted_clade"))
        _append_unique_typing(parts, section.get("predicted_serotype"))
    elif mode == "orthoebolavirus_typing":
        _append_unique_typing(parts, section.get("predicted_clade"))
        _append_unique_typing(parts, section.get("predicted_group"))
    elif mode == "astroviridae_typing":
        _append_unique_typing(parts, section.get("predicted_clade"))
        _append_unique_typing(parts, section.get("predicted_lineage"))
        _append_unique_typing(parts, section.get("predicted_group"))
    elif mode == "rhinovirus_typing":
        _append_unique_typing(parts, section.get("predicted_clade"))
        _append_unique_typing(parts, section.get("predicted_group"))
    elif mode == "seasonal_hcov_typing":
        _append_unique_typing(parts, section.get("predicted_subtype"))
        _append_unique_typing(parts, section.get("predicted_clade"))
    elif mode == "rotavirus_typing":
        _append_unique_typing(parts, section.get("predicted_subtype"))
        _append_unique_typing(parts, section.get("predicted_group"))
    elif mode == "hpiv_typing":
        _append_unique_typing(parts, section.get("predicted_clade"))
    else:
        for key in (
            "predicted_subtype",
            "predicted_clade",
            "predicted_serotype",
            "predicted_group",
            "predicted_lineage",
            "predicted_outbreak",
        ):
            _append_unique_typing(parts, section.get(key))
    return " / ".join(parts[:4])


VIRUS_TASK_SPECIES_DISPLAY_NAMES: dict[str, str] = {
    "sars-cov-2": "新冠病毒",
    "severe acute respiratory syndrome coronavirus 2": "新冠病毒",
    "influenza virus": "流感病毒",
    "influenza a virus": "甲型流感病毒",
    "influenza b virus": "乙型流感病毒",
    "respiratory syncytial virus": "呼吸道合胞病毒",
    "human respiratory syncytial virus": "呼吸道合胞病毒",
    "human metapneumovirus": "人偏肺病毒",
    "metapneumovirus": "人偏肺病毒",
    "human parainfluenza virus": "人副流感病毒",
    "parainfluenza virus": "人副流感病毒",
    "adenovirus": "腺病毒",
    "human adenovirus": "腺病毒",
    "human rhinovirus": "鼻病毒",
    "rhinovirus": "鼻病毒",
    "human coronavirus": "季节性冠状病毒",
    "coronavirus": "季节性冠状病毒",
    "monkeypox virus": "猴痘病毒",
    "mpox virus": "猴痘病毒",
    "dengue virus": "登革热病毒",
    "zika virus": "寨卡病毒",
    "chikungunya virus": "基孔肯雅病毒",
    "bandavirus dabieense": "大别班达病毒",
    "orthohantavirus": "汉坦病毒",
    "hantavirus": "汉坦病毒",
    "severe fever with thrombocytopenia syndrome virus": "大别班达病毒",
    "norovirus": "诺如病毒",
    "rotavirus a": "轮状病毒",
    "rotavirus": "轮状病毒",
    "astrovirus": "星状病毒",
    "human enterovirus": "肠道病毒",
    "enterovirus": "肠道病毒",
    "sapovirus": "札如病毒",
}


def _resolve_import_species_name(task_species: object, fallback_species_name: object, *, analysis_target: object = "") -> str:
    fallback = str(fallback_species_name or "").strip()
    task_species_text = str(task_species or "").strip()
    if _normalize_pathogen_type(analysis_target) != "virus":
        return fallback or task_species_text
    if not task_species_text:
        return fallback
    if any("\u4e00" <= char <= "\u9fff" for char in task_species_text):
        return task_species_text
    normalized = task_species_text.lower()
    return VIRUS_TASK_SPECIES_DISPLAY_NAMES.get(normalized, task_species_text)


def _infer_viral_suspected_syndrome(species_name: object, serotype_section: object) -> str:
    species = str(species_name or "").strip().lower()
    section = serotype_section if isinstance(serotype_section, dict) else {}
    mode = str(section.get("mode") or "").strip()
    predicted_clade = str(section.get("predicted_clade") or "").strip().upper()
    if mode in {
        "influenza_typing",
        "rsv_nextclade",
        "hmpv_nextclade",
        "hpiv_typing",
        "hadv_typing",
        "rhinovirus_typing",
        "seasonal_hcov_typing",
        "sars_cov_2_nextclade",
    }:
        return "肺炎/呼吸道感染"
    if mode == "enterovirus_typing":
        return "脑膜炎/脑膜脑炎" if "EV-A71" in predicted_clade else "腹泻/肠道感染"
    if mode in {"norovirus_typing", "rotavirus_typing", "astroviridae_typing"}:
        return "腹泻/肠道感染"
    if mode == "bandavirus_typing":
        return "败血症/血流感染"
    if mode == "orthohantavirus_typing":
        return "不明原因感染"
    if mode in {"denv_nextclade", "zikav_nextclade", "chikv_nextclade", "ebola_nextclade"}:
        return "不明原因感染"
    if mode == "monkeypox_nextclade":
        return "皮肤软组织感染"
    if any(token in species for token in ["influenza", "syncytial", "metapneumovirus", "parainfluenza", "adenovirus", "rhinovirus", "coronavirus", "sars-cov-2"]):
        return "肺炎/呼吸道感染"
    if any(token in species for token in ["norovirus", "rotavirus", "astrovirus", "enterovirus"]):
        return "腹泻/肠道感染"
    if any(token in species for token in ["bandavirus", "sftsv", "dabie"]):
        return "败血症/血流感染"
    if any(token in species for token in ["orthohantavirus", "hantavirus", "汉坦", "汉他"]):
        return "不明原因感染"
    if any(token in species for token in ["orthoebolavirus", "ebolavirus", "ebola", "埃博拉"]):
        return "不明原因感染"
    if any(token in species for token in ["dengue", "zika", "chikungunya"]):
        return "不明原因感染"
    if any(token in species for token in ["mpox", "monkeypox"]):
        return "皮肤软组织感染"
    return ""


def _normalize_pathogen_type(value: object) -> str:
    text = str(value or "").strip().lower()
    if text in {"virus", "viral", "病毒"}:
        return "virus"
    if text in {"bacteria", "bacterial", "细菌"}:
        return "bacteria"
    return ""


def _infer_pathogen_type(species_name: object, *, analysis_target: object = "", mlst_st: object = "") -> str:
    explicit = _normalize_pathogen_type(analysis_target)
    if explicit:
        return explicit
    if str(mlst_st or "").strip():
        return "bacteria"
    species = str(species_name or "").strip().lower()
    if any(token in species for token in [
        "influenza virus",
        "respiratory syncytial virus",
        "metapneumovirus",
        "parainfluenza virus",
        "adenovirus",
        "rhinovirus",
        "coronavirus",
        "sars-cov-2",
        "norovirus",
        "rotavirus",
        "astrovirus",
        "enterovirus",
        "bandavirus",
        "orthohantavirus",
        "hantavirus",
        "sftsv",
        "dengue virus",
        "zika virus",
        "chikungunya virus",
        "mpox",
        "monkeypox virus",
    ]):
        return "virus"
    return "bacteria"


def _build_task_import_metadata_items(species_name: object, serotype_section: object) -> list[dict[str, Any]]:
    syndrome = _infer_viral_suspected_syndrome(species_name, serotype_section)
    if not syndrome:
        return []
    syndrome_options = next(
        (
            item.get("options") or []
            for item in DEFAULT_SAMPLE_METADATA_TEMPLATES
            if str(item.get("key") or "").strip() == "suspected_syndrome"
        ),
        [],
    )
    return [
        {
            "key": "suspected_syndrome",
            "label": "疑似症候群",
            "type": "select",
            "options": syndrome_options,
            "value": syndrome,
        }
    ]


def _extract_serotype_from_payload(payload: object) -> str:
    sections = payload.get("sections", {}) if isinstance(payload, dict) else {}
    serotype_section = sections.get("serotype") if isinstance(sections.get("serotype"), dict) else {}
    structured_value = _extract_structured_serotype_from_section(serotype_section)
    if structured_value:
        return structured_value
    tables: list[dict[str, Any]] = []
    if isinstance(sections.get("priority_serotype"), dict):
        tables.append(sections["priority_serotype"])
    if isinstance(sections.get("serotype"), dict):
        tables.append(sections["serotype"])

    candidates = {"知识库命中血清型", "血清型", "血清群", "亚型预测", "predicted_serotype", "Serotype", "New_serotype"}
    for table in tables:
        columns = table.get("columns") if isinstance(table.get("columns"), list) else []
        rows = table.get("rows") if isinstance(table.get("rows"), list) else []
        if not columns or not rows:
            continue
        for idx, label in enumerate(columns):
            if str(label or "").strip() not in candidates:
                continue
            for row in rows:
                if not isinstance(row, list) or idx >= len(row):
                    continue
                value = str(row[idx] or "").strip()
                if value and value != "-":
                    return value
    return ""


def _find_table_column_index(columns: list[object], candidates: tuple[str, ...]) -> int:
    labels = [str(label or "").strip() for label in (columns or [])]
    for candidate in candidates:
        if candidate in labels:
            return labels.index(candidate)
    return -1


def _summarize_detected_genes(table: object, *, limit: int = 8) -> str:
    if not isinstance(table, dict):
        return ""
    columns = table.get("columns") if isinstance(table.get("columns"), list) else []
    rows = table.get("rows") if isinstance(table.get("rows"), list) else []
    gene_index = _find_table_column_index(columns, ("基因名称", "gene_name", "Gene"))
    if gene_index < 0 or not rows:
        return ""
    seen: list[str] = []
    for row in rows:
        if not isinstance(row, list) or gene_index >= len(row):
            continue
        gene_name = str(row[gene_index] or "").strip()
        if not gene_name or gene_name == "-" or gene_name in seen:
            continue
        seen.append(gene_name)
    if not seen:
        return ""
    if len(seen) <= limit:
        return "、".join(seen)
    return f"{'、'.join(seen[:limit])} 等 {len(seen)} 项"


def _summarize_mge_gene_hits(table: object, *, limit: int = 6) -> str:
    if not isinstance(table, dict):
        return ""
    columns = table.get("columns") if isinstance(table.get("columns"), list) else []
    rows = table.get("rows") if isinstance(table.get("rows"), list) else []
    gene_index = _find_table_column_index(columns, ("基因名称", "gene_name", "Gene"))
    mge_index = _find_table_column_index(columns, ("关联元件类型", "元件类型", "MGE类型"))
    if gene_index < 0 or mge_index < 0 or not rows:
        return ""
    seen: list[str] = []
    for row in rows:
        if not isinstance(row, list) or gene_index >= len(row) or mge_index >= len(row):
            continue
        gene_name = str(row[gene_index] or "").strip()
        mge_type = str(row[mge_index] or "").strip()
        if not gene_name or gene_name == "-":
            continue
        normalized_types = [item.strip() for item in mge_type.split("、") if item.strip() and item.strip() != "未识别"]
        if not normalized_types:
            continue
        summary = f"{gene_name}（{'/'.join(normalized_types)}）"
        if summary in seen:
            continue
        seen.append(summary)
    if not seen:
        return ""
    if len(seen) <= limit:
        return "、".join(seen)
    return f"{'、'.join(seen[:limit])} 等 {len(seen)} 项"


@dataclass
class SampleLibraryManager:
    store: PortalStore
    resolve_report_source: ReportSourceResolver
    build_report_payload: ReportPayloadBuilder
    resolve_report_sample_name: SampleNameResolver
    human_bp: HumanBpFormatter

    def list_visible(self, *, scope: str, role: str, username: str, group_name: str) -> list[dict[str, Any]]:
        rows = self.store.list_sample_library_by_scope(scope)
        return [self._enrich_record_permissions(row, role=role, username=username, group_name=group_name) for row in rows if self._can_view_record(row, role=role, username=username, group_name=group_name)]

    def get_visible(self, sample_key: str, *, role: str, username: str, group_name: str) -> dict[str, Any]:
        record = self.store.get_sample_library_record(sample_key)
        if not self._can_view_record(record, role=role, username=username, group_name=group_name):
            raise KeyError(f"样本不存在: {sample_key}")
        return self._enrich_record_permissions(record, role=role, username=username, group_name=group_name)

    def list_submissions(self, *, role: str, username: str, group_name: str) -> list[dict[str, Any]]:
        rows = self.store.list_sample_library_submissions()
        visible: list[dict[str, Any]] = []
        for row in rows:
            if role == "admin":
                visible.append({**row, "can_review": row.get("status") == "pending"})
                continue
            if row.get("owner") == username:
                visible.append({**row, "can_review": False})
        return visible

    def list_version_logs(self, *, role: str, username: str, group_name: str) -> list[dict[str, Any]]:
        rows = self.store.list_sample_library_version_logs(limit=60)
        if role == "admin":
            return rows
        return [row for row in rows if str(row.get("owner") or "") == username]

    def list_release_versions(self, *, role: str, username: str, group_name: str) -> list[dict[str, Any]]:
        rows = self.store.list_sample_library_release_versions(scope="main", limit=40)
        if role == "admin":
            return rows
        return []

    def publish_release_version(
        self,
        *,
        version_label: str,
        note: str,
        role: str,
        username: str,
        group_name: str,
    ) -> dict[str, Any]:
        if role != "admin":
            raise PermissionError("只有管理员可以发布数据库版本")
        normalized_label = str(version_label or "").strip() or self._default_release_label()
        current_snapshot = self._build_release_snapshot(scope="main")
        previous_release = next(iter(self.store.list_sample_library_release_versions(scope="main", limit=1)), None)
        previous_snapshot = self._parse_release_snapshot(previous_release.get("snapshot_json") if previous_release else "[]")
        changes = self._compare_release_snapshots(previous_snapshot, current_snapshot)
        summary = (
            "当前主数据库与上一版本完全一致，更新日志为空。"
            if changes["change_count"] <= 0
            else f"新增 {changes['added_count']} 条、更新 {changes['updated_count']} 条、删除 {changes['deleted_count']} 条。"
        )
        release = self.store.create_sample_library_release_version(
            {
                "release_id": f"release::{uuid4().hex}",
                "version_label": normalized_label,
                "scope": "main",
                "summary": summary,
                "note": str(note or "").strip(),
                "operator": username,
                "change_count": changes["change_count"],
                "added_count": changes["added_count"],
                "updated_count": changes["updated_count"],
                "deleted_count": changes["deleted_count"],
                "changes_json": json.dumps(changes["changes"], ensure_ascii=False),
                "snapshot_json": json.dumps(current_snapshot, ensure_ascii=False),
                "created_at": self._now(),
            }
        )
        self.store.create_sample_library_version_log(
            {
                "event_id": f"sample-version::{uuid4().hex}",
                "version_label": normalized_label,
                "request_id": "",
                "sample_key": "",
                "sample_name": normalized_label,
                "owner": "",
                "owner_group": "",
                "action": "release",
                "summary": summary,
                "operator": username,
                "payload_json": json.dumps(
                    {
                        "release_id": release["release_id"],
                        "scope": "main",
                        "note": str(note or "").strip(),
                        "change_count": changes["change_count"],
                        "added_count": changes["added_count"],
                        "updated_count": changes["updated_count"],
                        "deleted_count": changes["deleted_count"],
                    },
                    ensure_ascii=False,
                ),
                "created_at": release["created_at"],
            }
        )
        return release

    def list_metadata_templates(self) -> list[dict[str, Any]]:
        return self._merge_default_metadata_templates(self.store.list_sample_library_metadata_templates())

    def save_metadata_templates(self, metadata_templates: list[dict[str, Any]]) -> list[dict[str, Any]]:
        normalized_templates = self._build_shared_metadata_templates(metadata_templates)
        self.store.upsert_sample_library_metadata_templates(normalized_templates)
        return self.list_metadata_templates()

    def validate_metadata_items(
        self,
        metadata_items: list[dict[str, Any]] | None,
        *,
        existing_metadata_json: str | None = None,
        partial: bool = False,
    ) -> None:
        template_map = {
            str(item.get("key") or "").strip(): item
            for item in self.list_metadata_templates()
            if str(item.get("key") or "").strip()
        }
        merged_items: dict[str, dict[str, Any]] = {}
        for item in _parse_metadata_items(existing_metadata_json):
            key = str(item.get("key") or "").strip()
            if not key:
                continue
            merged_items[key] = item
        for item in metadata_items or []:
            if not isinstance(item, dict):
                continue
            key = str(item.get("key") or "").strip()
            if not key:
                continue
            merged_items[key] = item

        errors: list[str] = []
        for key, template in template_map.items():
            requirement = str(template.get("requirement") or "").strip()
            if requirement != "required":
                continue
            if partial and key not in merged_items:
                continue
            current = merged_items.get(key) or {}
            field_type = str(current.get("type") or template.get("type") or "text").strip() or "text"
            value = current.get("value")
            if field_type == "location":
                parts = value if isinstance(value, dict) else {}
                is_valid = any(str(parts.get(part) or "").strip() for part in ("province", "city", "district", "detail"))
            else:
                is_valid = bool(str(value or "").strip())
            if not is_valid:
                errors.append(f"{template.get('label') or key} 为必填字段")

        for key, item in merged_items.items():
            template = template_map.get(key) or {}
            field_type = str(item.get("type") or template.get("type") or "text").strip() or "text"
            if field_type != "select":
                continue
            options = item.get("options") if isinstance(item.get("options"), list) else template.get("options")
            options = [str(option).strip() for option in (options or []) if str(option).strip()]
            value = str(item.get("value") or "").strip()
            if value and options and value not in options:
                errors.append(f"{template.get('label') or item.get('label') or key} 不在标准词表内")

        if errors:
            raise ValueError("；".join(errors))

    def batch_update_visible(
        self,
        sample_keys: list[str],
        *,
        role: str,
        username: str,
        group_name: str,
        custom_metadata_json: list[dict[str, Any]] | str | None = None,
        metadata_templates: list[dict[str, Any]] | None = None,
    ) -> dict[str, Any]:
        normalized_keys = [str(item or "").strip() for item in sample_keys if str(item or "").strip()]
        if not normalized_keys:
            raise ValueError("未选择样本")
        incoming_items = _parse_metadata_items(custom_metadata_json) if isinstance(custom_metadata_json, str) else (
            [
                item for item in (custom_metadata_json or [])
                if isinstance(item, dict)
            ]
            if custom_metadata_json is not None else []
        )
        update_map: dict[str, dict[str, Any]] = {}
        for item in incoming_items:
            key = str(item.get("key") or "").strip()
            if not key:
                continue
            normalized = {
                "key": key,
                "label": str(item.get("label") or key).strip(),
                "type": str(item.get("type") or "text").strip() or "text",
                "options": [str(option).strip() for option in (item.get("options") or []) if str(option).strip()]
                if isinstance(item.get("options"), list)
                else [],
                "value": (
                    {
                        "province": str((item.get("value") or {}).get("province") or "").strip(),
                        "city": str((item.get("value") or {}).get("city") or "").strip(),
                        "district": str((item.get("value") or {}).get("district") or "").strip(),
                        "detail": str((item.get("value") or {}).get("detail") or "").strip(),
                    }
                    if str(item.get("type") or "text").strip() == "location"
                    else str(item.get("value") or "").strip()
                ),
            }
            if normalized["type"] == "location":
                if not any(str(part).strip() for part in normalized["value"].values()):
                    continue
            elif not str(normalized["value"] or "").strip():
                continue
            update_map[key] = normalized
        if not update_map:
            raise ValueError("没有可批量更新的主档字段")
        self.validate_metadata_items(list(update_map.values()), partial=True)

        updated_count = 0
        for sample_key in normalized_keys:
            record = self.get_visible(sample_key, role=role, username=username, group_name=group_name)
            if not self._can_edit_record(record, role=role, username=username):
                continue
            existing_items = {
                str(item.get("key") or "").strip(): item
                for item in _parse_metadata_items(record.get("custom_metadata_json"))
                if isinstance(item, dict) and str(item.get("key") or "").strip()
            }
            for key, item in update_map.items():
                existing_items[key] = item
            self.store.update_sample_library_record(
                record["sample_key"],
                custom_metadata_json=json.dumps(list(existing_items.values()), ensure_ascii=False),
            )
            updated_count += 1
        if metadata_templates is not None:
            normalized_templates = self._build_shared_metadata_templates(metadata_templates)
            self.store.upsert_sample_library_metadata_templates(normalized_templates)
        return {"updated_count": updated_count, "selected_count": len(normalized_keys)}

    def update_visible(
        self,
        sample_key: str,
        *,
        role: str,
        username: str,
        group_name: str,
        genome_id: str | None = None,
        pathogen_type: str | None = None,
        sample_alias: str | None = None,
        taxid: str | None = None,
        mlst_st: str | None = None,
        serotype_result: str | None = None,
        resistance_gene_hits: str | None = None,
        virulence_gene_hits: str | None = None,
        resistance_mge_hits: str | None = None,
        virulence_mge_hits: str | None = None,
        description: str | None = None,
        sample_source: str | None = None,
        collection_date: str | None = None,
        gender: str | None = None,
        country: str | None = None,
        host_info: str | None = None,
        location_json: str | None = None,
        sample_type: str | None = None,
        sequencing_method: str | None = None,
        genome_length: str | None = None,
        note: str | None = None,
        visibility_scope: str | None = None,
        custom_metadata_json: list[dict[str, str]] | str | None = None,
        metadata_templates: list[dict[str, Any]] | None = None,
    ) -> dict[str, Any]:
        record = self.get_visible(sample_key, role=role, username=username, group_name=group_name)
        if not self._can_edit_record(record, role=role, username=username):
            raise PermissionError("当前用户无权修改此样本")
        metadata_json = None
        if custom_metadata_json is not None:
            if isinstance(custom_metadata_json, str):
                self.validate_metadata_items(_parse_metadata_items(custom_metadata_json))
                metadata_json = custom_metadata_json
            else:
                normalized_metadata = [
                    {
                        "key": str(item.get("key") or "").strip(),
                        "label": str(item.get("label") or item.get("key") or "").strip(),
                        "type": str(item.get("type") or "text").strip() or "text",
                        "options": [
                            str(option).strip()
                            for option in (item.get("options") or [])
                            if str(option).strip()
                        ] if isinstance(item.get("options"), list) else [],
                        "value": (
                            {
                                "province": str((item.get("value") or {}).get("province") or "").strip(),
                                "city": str((item.get("value") or {}).get("city") or "").strip(),
                                "district": str((item.get("value") or {}).get("district") or "").strip(),
                                "detail": str((item.get("value") or {}).get("detail") or "").strip(),
                            }
                            if str(item.get("type") or "text").strip() == "location"
                            else str(item.get("value") or "").strip()
                        ),
                    }
                    for item in custom_metadata_json
                    if isinstance(item, dict) and (
                        str(item.get("key") or "").strip()
                        or str(item.get("label") or "").strip()
                        or str(item.get("value") or "").strip()
                        or (isinstance(item.get("value"), dict) and any(str(v).strip() for v in item.get("value").values()))
                    )
                ]
                self.validate_metadata_items(normalized_metadata)
                metadata_json = json.dumps(normalized_metadata, ensure_ascii=False)
        updated = self.store.update_sample_library_record(
            record["sample_key"],
            genome_id=genome_id,
            pathogen_type=_normalize_pathogen_type(pathogen_type),
            sample_alias=sample_alias,
            taxid=taxid,
            mlst_st=mlst_st,
            serotype_result=serotype_result,
            resistance_gene_hits=resistance_gene_hits,
            virulence_gene_hits=virulence_gene_hits,
            resistance_mge_hits=resistance_mge_hits,
            virulence_mge_hits=virulence_mge_hits,
            description=description,
            sample_source=sample_source,
            collection_date=collection_date,
            gender=gender,
            country=country,
            host_info=host_info,
            location_json=location_json,
            sample_type=sample_type,
            sequencing_method=sequencing_method,
            genome_length=genome_length,
            note=note,
            visibility_scope=visibility_scope,
            custom_metadata_json=metadata_json,
        )
        if metadata_templates is not None:
            normalized_templates = self._build_shared_metadata_templates(metadata_templates)
            self.store.upsert_sample_library_metadata_templates(normalized_templates)
        return self.get_visible(updated["sample_key"], role=role, username=username, group_name=group_name)

    def delete_visible(self, sample_key: str, *, role: str, username: str, group_name: str) -> None:
        record = self.get_visible(sample_key, role=role, username=username, group_name=group_name)
        if not self._can_delete_record(record, role=role, username=username):
            raise PermissionError("当前用户无权删除此样本")
        self.store.delete_sample_library_record(sample_key)

    def submit_personal_to_main(self, sample_key: str, *, role: str, username: str, group_name: str) -> dict[str, Any]:
        record = self.get_visible(sample_key, role=role, username=username, group_name=group_name)
        if str(record.get("library_scope") or "") != "personal":
            raise ValueError("只有个人数据库样本可以提交至主数据库")
        owner = str(record.get("owner") or "").strip()
        if role != "admin" and owner != username:
            raise PermissionError("只能提交自己的个人数据库样本")
        existing = self.store.find_pending_submission_for_sample(sample_key)
        owner_group = str(record.get("owner_group") or "").strip() or group_name
        payload_json = json.dumps(record, ensure_ascii=False)
        if existing:
            return self.store.upsert_sample_library_submission(
                {
                    **existing,
                    "sample_name": record.get("sample_name") or "",
                    "owner": owner or username,
                    "owner_group": owner_group,
                    "payload_json": payload_json,
                    "status": "pending",
                    "review_note": "",
                    "reviewed_by": "",
                    "reviewed_at": "",
                }
            )
        return self.store.upsert_sample_library_submission(
            {
                "request_id": f"submission::{uuid4().hex}",
                "personal_sample_key": sample_key,
                "sample_name": record.get("sample_name") or "",
                "owner": owner or username,
                "owner_group": owner_group,
                "payload_json": payload_json,
                "status": "pending",
                "review_note": "",
                "reviewed_by": "",
                "reviewed_at": "",
            }
        )

    def review_submission(self, request_id: str, *, action: str, admin_username: str, note: str = "") -> dict[str, Any]:
        submission = self.store.get_sample_library_submission(request_id)
        if action not in {"approve", "reject"}:
            raise ValueError("action must be approve or reject")
        if submission.get("status") != "pending":
            raise ValueError("当前审核申请已处理")
        payload = json.loads(str(submission.get("payload_json") or "{}") or "{}")
        if action == "approve":
            main_record = dict(payload)
            main_record["sample_key"] = f"main::{request_id}"
            main_record["library_scope"] = "main"
            main_record["visibility_scope"] = "group"
            main_record["source_submission_id"] = request_id
            self.store.upsert_sample_library_record(main_record)
            version_stamp = self._now().replace(":", "").replace("-", "").replace("T", "-")[:15]
            self.store.create_sample_library_version_log(
                {
                    "event_id": f"sample-version::{uuid4().hex}",
                    "version_label": f"main-publish-{version_stamp}",
                    "request_id": request_id,
                    "sample_key": main_record["sample_key"],
                    "sample_name": submission.get("sample_name") or main_record.get("sample_name") or "",
                    "owner": submission.get("owner") or "",
                    "owner_group": submission.get("owner_group") or "",
                    "action": "publish",
                    "summary": "个人库样本经审核发布入主数据库",
                    "operator": admin_username,
                    "payload_json": json.dumps(
                        {
                            "source_submission_id": request_id,
                            "personal_sample_key": submission.get("personal_sample_key") or "",
                            "published_sample_key": main_record["sample_key"],
                            "review_note": note,
                        },
                        ensure_ascii=False,
                    ),
                    "created_at": self._now(),
                }
            )
        return self.store.upsert_sample_library_submission(
            {
                **submission,
                "status": "approved" if action == "approve" else "rejected",
                "review_note": note,
                "reviewed_by": admin_username,
                "reviewed_at": self._now(),
            }
        )

    def import_task_samples(self, task: dict, *, library_scope: str = "personal") -> dict[str, Any]:
        report_source = self.resolve_report_source(task, "")
        sample_names = list(report_source.get("samples") or [])
        if report_source.get("mode") != "multi":
            fallback_name = self.resolve_report_sample_name(task, report_source["report_dir"])
            sample_names = [fallback_name] if fallback_name else [""]

        imported_rows: list[dict[str, Any]] = []
        skipped: list[dict[str, str]] = []
        for sample_name in sample_names:
            payload = self.build_report_payload(task, sample_name)
            sample_task = payload.get("task") or {}
            task_params = task.get("params") or {}
            normalized_sample = str(sample_task.get("sample_name") or sample_name or "").strip()
            current_report_source = self.resolve_report_source(task, normalized_sample)
            report_dir = Path(str(current_report_source.get("report_dir") or report_source.get("report_dir") or "")).resolve()
            if not normalized_sample:
                skipped.append({"sample_name": sample_name or "-", "reason": "未识别样本名称"})
                continue

            final_fasta = report_dir / f"{normalized_sample}.final.fasta"
            if not final_fasta.is_file():
                skipped.append({"sample_name": normalized_sample, "reason": f"缺少 {normalized_sample}.final.fasta"})
                continue

            metrics = {item.get("key"): item for item in payload.get("overview_metrics", [])}
            species_items = metrics.get("species_estimation", {}).get("items", [])
            q_items = metrics.get("q_metrics", {}).get("items", [])
            checkm_items = metrics.get("checkm_metrics", {}).get("items", [])
            assembly_metric = metrics.get("assembly_profile", {})
            rv_sections = payload.get("sections", {}).get("resistance_virulence", {})
            mge_sections = payload.get("sections", {}).get("mge_monitoring", {})
            serotype_section = payload.get("sections", {}).get("serotype", {})
            genome_length = self._read_fasta_total_length(final_fasta) or assembly_metric.get("total_length")
            fallback_species_name = str(species_items[0].get("display") if len(species_items) > 0 else "") or ""
            species_name = _resolve_import_species_name(
                task_params.get("species"),
                fallback_species_name,
                analysis_target=task_params.get("analysis_target"),
            )
            metadata_items = _build_task_import_metadata_items(species_name, serotype_section)

            record = self.store.upsert_sample_library_record(
                {
                    "sample_key": f"{library_scope}::{task.get('id')}::{normalized_sample}",
                    "genome_id": normalized_sample,
                    "sample_name": normalized_sample,
                    "task_id": str(task.get("id") or ""),
                    "task_name": str(task.get("name") or ""),
                    "owner": str(task.get("owner") or ""),
                    "owner_group": str(task.get("owner_group") or ""),
                    "report_dir": str(report_dir),
                    "output_dir": str(task_params.get("output_dir", "") or ""),
                    "final_fasta_path": str(final_fasta),
                    "species_name": species_name,
                    "pathogen_type": _infer_pathogen_type(
                        species_name,
                        analysis_target=task_params.get("analysis_target"),
                        mlst_st=_extract_mlst_st_from_payload(payload),
                    ),
                    "taxid": "",
                    "mlst_species_name": str(species_items[1].get("display") if len(species_items) > 1 else "") or "",
                    "mlst_st": _extract_mlst_st_from_payload(payload),
                    "serotype_result": _extract_serotype_from_payload(payload),
                    "genome_length": str(genome_length or ""),
                    "q20_rate": str(q_items[0].get("display") if len(q_items) > 0 else "") or "",
                    "q30_rate": str(q_items[1].get("display") if len(q_items) > 1 else "") or "",
                    "completeness": str(checkm_items[0].get("display") if len(checkm_items) > 0 else "") or "",
                    "contamination": str(checkm_items[1].get("display") if len(checkm_items) > 1 else "") or "",
                    "contig_count": str(assembly_metric.get("contig_count") if assembly_metric.get("contig_count") is not None else ""),
                    "plasmid_count": str(assembly_metric.get("plasmid_count") if assembly_metric.get("plasmid_count") is not None else ""),
                    "total_length": self.human_bp(assembly_metric.get("total_length")),
                    "resistance_count": str(len(rv_sections.get("resistance_elements", {}).get("rows", []) or [])),
                    "virulence_count": str(len(rv_sections.get("virulence_elements", {}).get("rows", []) or [])),
                    "resistance_gene_hits": _summarize_detected_genes(rv_sections.get("resistance_elements")),
                    "virulence_gene_hits": _summarize_detected_genes(rv_sections.get("virulence_elements")),
                    "resistance_mge_hits": _summarize_mge_gene_hits(mge_sections.get("resistance")),
                    "virulence_mge_hits": _summarize_mge_gene_hits(mge_sections.get("virulence")),
                    "description": "",
                    "gender": "",
                    "country": "",
                    "location_json": "",
                    "sample_type": str(task_params.get("inputtype", "") or ""),
                    "sequencing_method": str(task_params.get("method", "") or ""),
                    "custom_metadata_json": json.dumps(metadata_items, ensure_ascii=False),
                    "library_scope": library_scope,
                    "visibility_scope": "public" if library_scope == "main" and not str(task.get("owner_group") or "").strip() else "group",
                    "source_submission_id": "",
                }
            )
            imported_rows.append(record)

        return {
            "status": "ok",
            "imported_count": len(imported_rows),
            "skipped_count": len(skipped),
            "items": imported_rows,
            "skipped": skipped,
        }

    def preview_meta_task_bins(self, task: dict) -> dict[str, Any]:
        report_source = self.resolve_report_source(task, "")
        report_dir = Path(str(report_source.get("report_dir") or "")).resolve()
        payload = self.build_report_payload(task, "")
        binning_section = ((payload.get("sections") or {}).get("binning_results") or {})
        quality_table = ((binning_section.get("quality") or {}).get("table") or {})
        taxonomy_table = ((binning_section.get("taxonomy") or {}).get("table") or {})
        quality_columns = quality_table.get("columns") or []
        quality_rows = quality_table.get("rows") or []
        taxonomy_columns = taxonomy_table.get("columns") or []
        taxonomy_rows = taxonomy_table.get("rows") or []

        taxonomy_map: dict[str, str] = {}
        if "Bin名称" in taxonomy_columns and "种" in taxonomy_columns:
            name_index = taxonomy_columns.index("Bin名称")
            species_index = taxonomy_columns.index("种")
            for row in taxonomy_rows:
                bin_name = str(row[name_index] if name_index < len(row) else "").strip()
                species_name = str(row[species_index] if species_index < len(row) else "").strip()
                if bin_name:
                    taxonomy_map[bin_name] = species_name or "-"

        items: list[dict[str, Any]] = []
        if {"Bin名称", "完整性", "污染率"}.issubset(set(quality_columns)):
            name_index = quality_columns.index("Bin名称")
            completeness_index = quality_columns.index("完整性")
            contamination_index = quality_columns.index("污染率")
            for row in quality_rows:
                bin_name = str(row[name_index] if name_index < len(row) else "").strip()
                if not bin_name:
                    continue
                fasta_path = self._resolve_meta_bin_fasta(report_dir, bin_name)
                items.append(
                    {
                        "bin_name": bin_name,
                        "completeness": str(row[completeness_index] if completeness_index < len(row) else "").strip(),
                        "contamination": str(row[contamination_index] if contamination_index < len(row) else "").strip(),
                        "species_name": taxonomy_map.get(bin_name, "-"),
                        "fasta_path": str(fasta_path) if fasta_path else "",
                        "available": bool(fasta_path and fasta_path.is_file()),
                    }
                )
        return {"status": "ok", "items": items, "report_dir": str(report_dir)}

    def import_meta_task_bins(self, task: dict, *, selected_bins: list[str], library_scope: str = "personal") -> dict[str, Any]:
        preview = self.preview_meta_task_bins(task)
        preview_map = {str(item.get("bin_name") or "").strip(): item for item in preview.get("items") or []}
        selected = [str(item or "").strip() for item in selected_bins if str(item or "").strip()]
        params = task.get("params") or {}
        imported_rows: list[dict[str, Any]] = []
        skipped: list[dict[str, str]] = []
        report_dir = Path(str(preview.get("report_dir") or "")).resolve() if preview.get("report_dir") else Path(".").resolve()

        for bin_name in selected:
            item = preview_map.get(bin_name)
            if not item:
                skipped.append({"sample_name": bin_name, "reason": "未找到该 bin 的预览信息"})
                continue
            fasta_path = Path(str(item.get("fasta_path") or "")).expanduser()
            if not fasta_path.is_file():
                skipped.append({"sample_name": bin_name, "reason": "对应 bin fasta 不存在"})
                continue
            genome_length = self._read_fasta_total_length(fasta_path)
            record = self.store.upsert_sample_library_record(
                {
                    "sample_key": f"{library_scope}::{task.get('id')}::{bin_name}",
                    "genome_id": bin_name,
                    "sample_name": bin_name,
                    "task_id": str(task.get("id") or ""),
                    "task_name": str(task.get("name") or ""),
                    "owner": str(task.get("owner") or ""),
                    "owner_group": str(task.get("owner_group") or ""),
                    "report_dir": str(report_dir),
                    "output_dir": str(params.get("output_dir") or ""),
                    "final_fasta_path": str(fasta_path.resolve()),
                    "species_name": str(item.get("species_name") or "").strip(),
                    "taxid": "",
                    "mlst_species_name": "",
                    "genome_length": str(genome_length or ""),
                    "q20_rate": "",
                    "q30_rate": "",
                    "completeness": str(item.get("completeness") or "").strip(),
                    "contamination": str(item.get("contamination") or "").strip(),
                    "contig_count": "",
                    "plasmid_count": "",
                    "total_length": self.human_bp(genome_length),
                    "resistance_count": "",
                    "virulence_count": "",
                    "description": "",
                    "gender": "",
                    "country": "",
                    "location_json": "",
                    "sample_type": "fasta",
                    "sequencing_method": str(params.get("method") or "").strip(),
                    "custom_metadata_json": "[]",
                    "library_scope": library_scope,
                    "visibility_scope": "public" if library_scope == "main" and not str(task.get("owner_group") or "").strip() else "group",
                    "source_submission_id": "",
                }
            )
            imported_rows.append(record)

        return {
            "status": "ok",
            "imported_count": len(imported_rows),
            "skipped_count": len(skipped),
            "items": imported_rows,
            "skipped": skipped,
        }

    def import_local_sample(
        self,
        *,
        owner: str,
        owner_group: str,
        library_scope: str = "personal",
        sample_name: str,
        final_fasta_path: str,
        pathogen_type: str = "",
        species_name: str = "",
        mlst_species_name: str = "",
        mlst_st: str = "",
        serotype_result: str = "",
        resistance_gene_hits: str = "",
        virulence_gene_hits: str = "",
        resistance_mge_hits: str = "",
        virulence_mge_hits: str = "",
        q20_rate: str = "",
        q30_rate: str = "",
        completeness: str = "",
        contamination: str = "",
        contig_count: str = "",
        plasmid_count: str = "",
        resistance_count: str = "",
        virulence_count: str = "",
        genome_id: str = "",
        taxid: str = "",
        description: str = "",
        sample_source: str = "",
        collection_date: str = "",
        gender: str = "",
        country: str = "",
        host_info: str = "",
        location_json: str = "",
        sample_type: str = "",
        sequencing_method: str = "",
        note: str = "",
        custom_metadata_json: str = "[]",
    ) -> dict[str, Any]:
        normalized_sample = str(sample_name or "").strip()
        if not normalized_sample:
            raise ValueError("样本名称不能为空")
        final_fasta = Path(str(final_fasta_path or "").strip()).expanduser().resolve()
        if not final_fasta.is_file():
            raise ValueError(f"final.fasta 不存在: {final_fasta}")
        genome_length = self._read_fasta_total_length(final_fasta)
        record = self.store.upsert_sample_library_record(
            {
                "sample_key": f"{library_scope}::local::{final_fasta}",
                "genome_id": str(genome_id or normalized_sample).strip(),
                "sample_name": normalized_sample,
                "task_id": "",
                "task_name": "本地导入",
                "owner": str(owner or "").strip(),
                "owner_group": str(owner_group or "").strip(),
                "report_dir": str(final_fasta.parent),
                "output_dir": str(final_fasta.parent),
                "final_fasta_path": str(final_fasta),
                "species_name": str(species_name or "").strip(),
                "pathogen_type": _normalize_pathogen_type(pathogen_type) or _infer_pathogen_type(species_name, mlst_st=mlst_st),
                "taxid": str(taxid or "").strip(),
                "mlst_species_name": str(mlst_species_name or "").strip(),
                "mlst_st": str(mlst_st or "").strip(),
                "serotype_result": str(serotype_result or "").strip(),
                "resistance_gene_hits": str(resistance_gene_hits or "").strip(),
                "virulence_gene_hits": str(virulence_gene_hits or "").strip(),
                "resistance_mge_hits": str(resistance_mge_hits or "").strip(),
                "virulence_mge_hits": str(virulence_mge_hits or "").strip(),
                "genome_length": str(genome_length or ""),
                "q20_rate": str(q20_rate or "").strip(),
                "q30_rate": str(q30_rate or "").strip(),
                "completeness": str(completeness or "").strip(),
                "contamination": str(contamination or "").strip(),
                "contig_count": str(contig_count or "").strip(),
                "plasmid_count": str(plasmid_count or "").strip(),
                "total_length": self.human_bp(genome_length),
                "resistance_count": str(resistance_count or "").strip(),
                "virulence_count": str(virulence_count or "").strip(),
                "description": str(description or "").strip(),
                "gender": str(gender or "").strip(),
                "country": str(country or "").strip(),
                "location_json": str(location_json or "").strip(),
                "sample_type": str(sample_type or "").strip(),
                "sequencing_method": str(sequencing_method or "").strip(),
                "custom_metadata_json": str(custom_metadata_json or "[]").strip() or "[]",
                "sample_alias": "",
                "sample_source": str(sample_source or "").strip(),
                "collection_date": str(collection_date or "").strip(),
                "host_info": str(host_info or "").strip(),
                "note": str(note or "").strip(),
                "library_scope": library_scope,
                "visibility_scope": "public" if library_scope == "main" and not str(owner_group or "").strip() else "group",
                "source_submission_id": "",
            }
        )
        return {
            "status": "ok",
            "imported_count": 1,
            "skipped_count": 0,
            "items": [record],
            "skipped": [],
        }

    def _can_view_record(self, record: dict[str, Any], *, role: str, username: str, group_name: str) -> bool:
        if role == "admin":
            return True
        scope = str(record.get("library_scope") or "main")
        if scope == "main":
            owner_group = str(record.get("owner_group") or "").strip()
            if str(record.get("visibility_scope") or "group") == "public":
                return True
            return bool(group_name) and owner_group == group_name
        return str(record.get("owner") or "") == username

    def _can_edit_record(self, record: dict[str, Any], *, role: str, username: str) -> bool:
        scope = str(record.get("library_scope") or "main")
        owner = str(record.get("owner") or "")
        if scope == "main":
            return role == "admin"
        if role == "admin":
            return owner == username
        return owner == username

    def _can_delete_record(self, record: dict[str, Any], *, role: str, username: str) -> bool:
        return self._can_edit_record(record, role=role, username=username)

    def _enrich_record_permissions(self, record: dict[str, Any], *, role: str, username: str, group_name: str) -> dict[str, Any]:
        record = self._apply_demo_sample_metadata(record)
        scope = str(record.get("library_scope") or "main")
        owner = str(record.get("owner") or "")
        can_submit_to_main = scope == "personal" and (owner == username or role == "admin")
        pending_submission = self.store.find_pending_submission_for_sample(record["sample_key"]) if can_submit_to_main else None
        metadata_completion = self._build_metadata_completion(record.get("custom_metadata_json"))
        return {
            **record,
            "can_edit": self._can_edit_record(record, role=role, username=username),
            "can_delete": self._can_delete_record(record, role=role, username=username),
            "can_submit_to_main": can_submit_to_main,
            "pending_submission_status": pending_submission.get("status") if pending_submission else "",
            "pending_submission_id": pending_submission.get("request_id") if pending_submission else "",
            **metadata_completion,
        }

    def _build_metadata_completion(self, raw_metadata_json: object) -> dict[str, Any]:
        template_map = {
            str(item.get("key") or "").strip(): item
            for item in self.list_metadata_templates()
            if str(item.get("key") or "").strip()
        }
        items = {
            str(item.get("key") or "").strip(): item
            for item in _parse_metadata_items(raw_metadata_json)
            if str(item.get("key") or "").strip()
        }
        missing_required: list[str] = []
        missing_recommended: list[str] = []

        def _has_value(item: dict[str, Any], template: dict[str, Any]) -> bool:
            field_type = str(item.get("type") or template.get("type") or "text").strip() or "text"
            value = item.get("value")
            if field_type == "location":
                parts = value if isinstance(value, dict) else {}
                return any(str(parts.get(part) or "").strip() for part in ("province", "city", "district", "detail"))
            return bool(str(value or "").strip())

        for key, template in template_map.items():
            requirement = str(template.get("requirement") or "").strip()
            if requirement not in {"required", "recommended"}:
                continue
            current = items.get(key) or {}
            if _has_value(current, template):
                continue
            label = str(template.get("label") or key).strip()
            if requirement == "required":
                missing_required.append(label)
            else:
                missing_recommended.append(label)

        if missing_required:
            status = "missing_required"
            label = "缺关键字段"
            tone = "danger"
        elif missing_recommended:
            status = "missing_recommended"
            label = f"缺 {len(missing_recommended)} 项"
            tone = "warning"
        else:
            status = "complete"
            label = "主档完整"
            tone = "success"

        return {
            "metadata_completion_status": status,
            "metadata_completion_label": label,
            "metadata_completion_tone": tone,
            "metadata_missing_required_count": len(missing_required),
            "metadata_missing_recommended_count": len(missing_recommended),
            "metadata_missing_required_fields": missing_required,
            "metadata_missing_recommended_fields": missing_recommended,
            "metadata_missing_summary": "；".join(missing_required or missing_recommended),
        }

    def _apply_demo_sample_metadata(self, record: dict[str, Any]) -> dict[str, Any]:
        if not self._looks_like_demo_record(record):
            return record
        existing_items = _parse_metadata_items(record.get("custom_metadata_json"))
        if existing_items:
            return record

        sample_name = str(record.get("sample_name") or "").strip()
        species_name = str(record.get("species_name") or "").strip()
        location = self._parse_record_location(record.get("location_json"))
        sample_source = str(record.get("sample_source") or "").strip()
        collection_date = str(record.get("collection_date") or "").strip()
        seed = sum(ord(ch) for ch in sample_name) % 97
        city = location.get("city") or location.get("province") or "上海"
        district = location.get("district") or "示范区"
        detail = location.get("detail") or f"{city}示例监测点"
        syndrome = self._suggest_demo_syndrome(species_name, sample_source)
        surveillance_source = self._suggest_demo_source(sample_source, species_name, seed)
        cluster_status = self._suggest_demo_cluster_status(sample_source, seed)
        specimen_category = self._suggest_demo_specimen_category(sample_source)
        collection_unit = f"{city}疾控采样组"
        submitting_unit = f"{city}{district}示例送检单位"
        department = self._suggest_demo_department(syndrome, surveillance_source)
        sample_token = sample_name.split("_")[0].replace("GCF", "").replace("GCA", "").replace("|", "")[:10] or "DEMO"
        case_code = f"DEMO-{collection_date.replace('-', '') or '202603'}-{seed:02d}"
        patient_code = f"PAT-{sample_token or f'{seed:02d}'}"
        traditional_result = f"{species_name or '病原'} 传统检测提示阳性，建议与测序结果联合解读。"
        epidemiology_link = f"演示样本：{city}{district}区域监测链路，用于展示样本主档与结果统计联动。"
        metadata_items = [
            {"key": "case_id", "label": "病例/事件编号", "type": "text", "value": case_code},
            {"key": "patient_id", "label": "患者/个案编号", "type": "text", "value": patient_code},
            {
                "key": "surveillance_source",
                "label": "监测来源",
                "type": "select",
                "options": [item["options"] for item in DEFAULT_SAMPLE_METADATA_TEMPLATES if item["key"] == "surveillance_source"][0],
                "value": surveillance_source,
            },
            {
                "key": "suspected_syndrome",
                "label": "疑似症候群",
                "type": "select",
                "options": [item["options"] for item in DEFAULT_SAMPLE_METADATA_TEMPLATES if item["key"] == "suspected_syndrome"][0],
                "value": syndrome,
            },
            {"key": "submitting_unit", "label": "送检单位", "type": "text", "value": submitting_unit},
            {"key": "collection_unit", "label": "采样单位", "type": "text", "value": collection_unit},
            {"key": "ward_department", "label": "科室/病区", "type": "text", "value": department},
            {
                "key": "specimen_category",
                "label": "标本类别",
                "type": "select",
                "options": [item["options"] for item in DEFAULT_SAMPLE_METADATA_TEMPLATES if item["key"] == "specimen_category"][0],
                "value": specimen_category,
            },
            {
                "key": "collection_site",
                "label": "采样地点",
                "type": "location",
                "value": {
                    "province": location.get("province") or "",
                    "city": location.get("city") or "",
                    "district": location.get("district") or "",
                    "detail": detail,
                },
            },
            {
                "key": "cluster_status",
                "label": "聚集性状态",
                "type": "select",
                "options": [item["options"] for item in DEFAULT_SAMPLE_METADATA_TEMPLATES if item["key"] == "cluster_status"][0],
                "value": cluster_status,
            },
            {"key": "epidemiology_link", "label": "流行病学关联", "type": "text", "value": epidemiology_link},
            {"key": "traditional_result", "label": "传统检测结果", "type": "text", "value": traditional_result},
        ]
        return {
            **record,
            "custom_metadata_json": json.dumps(metadata_items, ensure_ascii=False),
        }

    def _looks_like_demo_record(self, record: dict[str, Any]) -> bool:
        text = " ".join(
            [
                str(record.get("note") or ""),
                str(record.get("description") or ""),
                str(record.get("report_dir") or ""),
                str(record.get("sample_name") or ""),
            ]
        ).lower()
        if "自动生成示例" in text or "demo_data" in text or "示例样本" in text:
            return True
        sample_name = str(record.get("sample_name") or "").strip().lower()
        return sample_name.startswith(("gcf_", "gca_", "rna"))

    def _parse_record_location(self, raw_value: object) -> dict[str, str]:
        text = str(raw_value or "").strip()
        if not text:
            return {"province": "", "city": "", "district": "", "detail": ""}
        try:
            parsed = json.loads(text)
        except Exception:
            return {"province": "", "city": "", "district": "", "detail": text}
        if not isinstance(parsed, dict):
            return {"province": "", "city": "", "district": "", "detail": ""}
        return {
            "province": str(parsed.get("province") or "").strip(),
            "city": str(parsed.get("city") or "").strip(),
            "district": str(parsed.get("district") or "").strip(),
            "detail": str(parsed.get("detail") or "").strip(),
        }

    def _suggest_demo_syndrome(self, species_name: str, sample_source: str) -> str:
        species = species_name.lower()
        source = sample_source.lower()
        if "环境" in sample_source or "污水" in sample_source or "食品" in sample_source:
            return "环境异常事件"
        if "血" in sample_source:
            return "败血症/血流感染"
        if "脑膜炎奈瑟" in species_name or "meningitidis" in species:
            return "脑膜炎/脑膜脑炎"
        if any(name in species for name in ["klebsiella", "acinetobacter", "pseudomonas", "staphylococcus aureus"]):
            return "医院感染监测"
        if "escherichia coli" in species or "coli" in species:
            return "腹泻/肠道感染" if "粪" in sample_source else "泌尿生殖道感染"
        return "不明原因感染"

    def _suggest_demo_source(self, sample_source: str, species_name: str, seed: int) -> str:
        if "环境" in sample_source:
            return "环境监测"
        if "污水" in sample_source:
            return "污水监测"
        if "食品" in sample_source:
            return "食品监测"
        if "脑膜炎奈瑟" in species_name or "Neisseria meningitidis" in species_name:
            return "发热门诊" if seed % 2 == 0 else "住院"
        options = ["门急诊", "住院", "ICU", "社区监测"]
        return options[seed % len(options)]

    def _suggest_demo_cluster_status(self, sample_source: str, seed: int) -> str:
        if "环境" in sample_source:
            return "聚集性"
        options = ["散发", "聚集性", "待判定", "暴发相关"]
        return options[seed % len(options)]

    def _suggest_demo_specimen_category(self, sample_source: str) -> str:
        if "血" in sample_source:
            return "血液"
        if "环境" in sample_source:
            return "环境拭子"
        if "污水" in sample_source or "水" in sample_source:
            return "污水/水体"
        if "食品" in sample_source:
            return "食品"
        return "其他"

    def _suggest_demo_department(self, syndrome: str, surveillance_source: str) -> str:
        if syndrome == "脑膜炎/脑膜脑炎":
            return "神经感染组"
        if syndrome == "败血症/血流感染":
            return "感染科 / 重症医学科"
        if syndrome == "医院感染监测":
            return "医院感染管理科"
        if "环境" in surveillance_source or "污水" in surveillance_source or "食品" in surveillance_source:
            return "公共卫生监测组"
        return "临床微生物室"

    def _now(self) -> str:
        from .store import utc_now_iso
        return utc_now_iso()

    def _default_release_label(self) -> str:
        stamp = self._now().replace("T", " ").replace("Z", "")
        return f"SampleDB {stamp}"

    def _build_release_snapshot(self, *, scope: str) -> list[dict[str, Any]]:
        rows = self.store.list_sample_library_by_scope(scope)
        snapshot = [self._normalize_release_record(row) for row in rows]
        snapshot.sort(key=lambda item: (str(item.get("sample_name") or ""), str(item.get("sample_key") or "")))
        return snapshot

    def _normalize_release_record(self, record: dict[str, Any]) -> dict[str, Any]:
        payload = {
            key: value
            for key, value in record.items()
            if key not in {"imported_at", "updated_at", "source_submission_id"}
        }
        payload["location_json"] = self._normalize_json_string(payload.get("location_json"), fallback={})
        payload["custom_metadata_json"] = self._normalize_json_string(payload.get("custom_metadata_json"), fallback=[])
        for key, value in list(payload.items()):
            if isinstance(value, str):
                payload[key] = value.strip()
        return payload

    def _normalize_json_string(self, raw_value: Any, *, fallback: object) -> str:
        text = str(raw_value or "").strip()
        if not text:
            return json.dumps(fallback, ensure_ascii=False, sort_keys=True)
        try:
            parsed = json.loads(text)
        except Exception:
            return json.dumps(fallback, ensure_ascii=False, sort_keys=True)
        return json.dumps(parsed, ensure_ascii=False, sort_keys=True)

    def _parse_release_snapshot(self, raw_snapshot: str | None) -> list[dict[str, Any]]:
        try:
            parsed = json.loads(str(raw_snapshot or "[]") or "[]")
        except Exception:
            parsed = []
        if not isinstance(parsed, list):
            return []
        normalized: list[dict[str, Any]] = []
        for item in parsed:
            if isinstance(item, dict):
                normalized.append(self._normalize_release_record(item))
        return normalized

    def _compare_release_snapshots(
        self,
        previous_snapshot: list[dict[str, Any]],
        current_snapshot: list[dict[str, Any]],
    ) -> dict[str, Any]:
        previous_map = {str(item.get("sample_key") or ""): item for item in previous_snapshot if str(item.get("sample_key") or "").strip()}
        current_map = {str(item.get("sample_key") or ""): item for item in current_snapshot if str(item.get("sample_key") or "").strip()}
        added_keys = sorted(set(current_map.keys()) - set(previous_map.keys()))
        deleted_keys = sorted(set(previous_map.keys()) - set(current_map.keys()))
        updated_items: list[dict[str, Any]] = []
        for sample_key in sorted(set(current_map.keys()) & set(previous_map.keys())):
            previous_item = previous_map[sample_key]
            current_item = current_map[sample_key]
            changed_fields = sorted(
                [
                    key
                    for key in set(previous_item.keys()) | set(current_item.keys())
                    if previous_item.get(key) != current_item.get(key)
                ]
            )
            if changed_fields:
                updated_items.append(
                    {
                        "sample_key": sample_key,
                        "sample_name": current_item.get("sample_name") or previous_item.get("sample_name") or sample_key,
                        "changed_fields": changed_fields,
                    }
                )
        return {
            "change_count": len(added_keys) + len(deleted_keys) + len(updated_items),
            "added_count": len(added_keys),
            "updated_count": len(updated_items),
            "deleted_count": len(deleted_keys),
            "changes": {
                "added": [
                    {
                        "sample_key": key,
                        "sample_name": current_map[key].get("sample_name") or key,
                    }
                    for key in added_keys
                ],
                "updated": updated_items,
                "deleted": [
                    {
                        "sample_key": key,
                        "sample_name": previous_map[key].get("sample_name") or key,
                    }
                    for key in deleted_keys
                ],
            },
        }

    def _read_fasta_total_length(self, path: Path) -> int | None:
        if not path.is_file():
            return None
        total = 0
        try:
            with path.open("r", encoding="utf-8", errors="ignore") as handle:
                for line in handle:
                    if line.startswith(">"):
                        continue
                    total += len(line.strip())
        except OSError:
            return None
        return total or None

    def _resolve_meta_bin_fasta(self, report_dir: Path, bin_name: str) -> Path | None:
        direct_candidates = [
            report_dir / "BASALT_out" / "meta_drep_out" / "binning_genomes" / f"{bin_name}.fa",
            report_dir / "BASALT_out" / "meta_drep_out" / "binning_genomes" / f"{bin_name}.fasta",
        ]
        for candidate in direct_candidates:
            if candidate.is_file():
                return candidate.resolve()
        for suffix in ("fa", "fasta", "fna"):
            matches = sorted(report_dir.rglob(f"{bin_name}.{suffix}"))
            for candidate in matches:
                if candidate.is_file():
                    return candidate.resolve()
        return None

    def _build_shared_metadata_templates(self, metadata_templates: list[dict[str, Any]]) -> list[dict[str, Any]]:
        template_map: dict[str, dict[str, Any]] = {}

        for item in self._merge_default_metadata_templates(metadata_templates):
            if not isinstance(item, dict):
                continue
            key = str(item.get("key") or "").strip()
            if not key:
                continue
            template_map[key] = {
                "key": key,
                "label": str(item.get("label") or key).strip(),
                "type": str(item.get("type") or "text").strip() or "text",
                "options": [str(option).strip() for option in (item.get("options") or []) if str(option).strip()]
                if isinstance(item.get("options"), list)
                else [],
                "group": str(item.get("group") or "").strip(),
                "requirement": str(item.get("requirement") or "").strip(),
                "placeholder": str(item.get("placeholder") or "").strip(),
                "help_text": str(item.get("help_text") or "").strip(),
                "dictionary_name": str(item.get("dictionary_name") or "").strip(),
            }

        for record in self.store.list_sample_library():
            raw_metadata = str(record.get("custom_metadata_json") or "").strip()
            if not raw_metadata:
                continue
            try:
                items = json.loads(raw_metadata)
            except Exception:
                continue
            if not isinstance(items, list):
                continue
            for item in items:
                if not isinstance(item, dict):
                    continue
                key = str(item.get("key") or "").strip()
                if not key:
                    continue
                template_map[key] = {
                    "key": key,
                    "label": str(item.get("label") or template_map.get(key, {}).get("label") or key).strip(),
                    "type": str(item.get("type") or template_map.get(key, {}).get("type") or "text").strip() or "text",
                    "options": [str(option).strip() for option in (item.get("options") or []) if str(option).strip()]
                    if isinstance(item.get("options"), list)
                    else template_map.get(key, {}).get("options", []),
                    "group": str(item.get("group") or template_map.get(key, {}).get("group") or "").strip(),
                    "requirement": str(item.get("requirement") or template_map.get(key, {}).get("requirement") or "").strip(),
                    "placeholder": str(item.get("placeholder") or template_map.get(key, {}).get("placeholder") or "").strip(),
                    "help_text": str(item.get("help_text") or template_map.get(key, {}).get("help_text") or "").strip(),
                    "dictionary_name": str(item.get("dictionary_name") or template_map.get(key, {}).get("dictionary_name") or "").strip(),
                }

        return list(template_map.values())

    def _merge_default_metadata_templates(self, items: list[dict[str, Any]] | None) -> list[dict[str, Any]]:
        merged: dict[str, dict[str, Any]] = {}
        for position, item in enumerate(DEFAULT_SAMPLE_METADATA_TEMPLATES):
            merged[str(item["key"])] = {
                "key": str(item["key"]),
                "label": str(item.get("label") or item["key"]),
                "type": str(item.get("type") or "text"),
                "options": list(item.get("options") or []),
                "group": str(item.get("group") or "").strip(),
                "requirement": str(item.get("requirement") or "").strip(),
                "placeholder": str(item.get("placeholder") or "").strip(),
                "help_text": str(item.get("help_text") or "").strip(),
                "dictionary_name": str(item.get("dictionary_name") or "").strip(),
                "position": position,
            }
        for item in items or []:
            if not isinstance(item, dict):
                continue
            key = str(item.get("key") or "").strip()
            if not key:
                continue
            current = merged.get(key, {"key": key, "label": key, "type": "text", "options": [], "position": len(merged)})
            merged[key] = {
                "key": key,
                "label": str(item.get("label") or current["label"]).strip(),
                "type": str(item.get("type") or current["type"]).strip() or "text",
                "options": [str(option).strip() for option in (item.get("options") or current["options"]) if str(option).strip()],
                "group": str(item.get("group") or current.get("group") or "").strip(),
                "requirement": str(item.get("requirement") or current.get("requirement") or "").strip(),
                "placeholder": str(item.get("placeholder") or current.get("placeholder") or "").strip(),
                "help_text": str(item.get("help_text") or current.get("help_text") or "").strip(),
                "dictionary_name": str(item.get("dictionary_name") or current.get("dictionary_name") or "").strip(),
                "position": int(item.get("position")) if str(item.get("position") or "").strip().isdigit() else current["position"],
                "updated_at": item.get("updated_at", current.get("updated_at", "")),
            }
        return sorted(merged.values(), key=lambda item: (int(item.get("position") or 0), str(item.get("label") or "")))
