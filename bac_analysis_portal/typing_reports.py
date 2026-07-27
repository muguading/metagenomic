from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from .knowledge_interpretation import _build_mlst_knowledge_summary, _build_tb_knowledge_summary, _normalize_mlst_lookup_text
from .table_io import _read_tsv_rows

TB_DRUG_NAME_ZH = {
    "Isoniazid": "异烟肼",
    "Rifampicin": "利福平",
    "Ethambutol": "乙胺丁醇",
    "Pyrazinamide": "吡嗪酰胺",
    "Levofloxacin": "左氧氟沙星",
    "Moxifloxacin": "莫西沙星",
    "Amikacin": "阿米卡星",
    "Kanamycin": "卡那霉素",
    "Capreomycin": "卷曲霉素",
    "Linezolid": "利奈唑胺",
    "Bedaquiline": "贝达喹啉",
    "Clofazimine": "氯法齐明",
    "Delamanid": "德拉马尼",
    "Ethionamide": "乙硫异烟胺",
    "Streptomycin": "链霉素",
    "Pretomanid": "普托马尼",
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

def _read_mlst_result(path: Path, project_root_text: str = "") -> dict:
    if not path.is_file():
        return {"columns": [], "rows": [], "gene_show_map": {}, "default_gene": "", "knowledge_summary": {"headline": "", "items": []}}
    raw = _read_tsv_rows(path)
    if not raw["columns"] or not raw["rows"]:
        return {"columns": [], "rows": [], "gene_show_map": {}, "default_gene": "", "knowledge_summary": {"headline": "", "items": []}}

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

    knowledge_summary = _build_mlst_knowledge_summary(project_root_text, rows) if project_root_text else {"headline": "", "items": []}
    return {
        "columns": ["Host Gene", "等位基因", "序列名称", "起始位置", "终止位置", "比对起始位置", "比对终止位置", "比对长度", "一致性%", "序列分型(ST)", "物种信息", "Host Gene 展示", "Host Gene ID"],
        "rows": rows,
        "gene_show_map": gene_show_map,
        "default_gene": default_gene or (rows[0][12] if rows else ""),
        "knowledge_summary": knowledge_summary,
    }

def _looks_like_neisseria_meningitidis_text(value: object) -> bool:
    text = str(value or "").strip().lower()
    if not text:
        return False
    return (
        "neisseria meningitidis" in text
        or "meningitidis" in text
        or "脑膜炎奈瑟" in text
        or "流脑" in text
    )

def _looks_like_tuberculosis_text(value: object) -> bool:
    text = str(value or "").strip().lower()
    if not text:
        return False
    return (
        "mycobacterium tuberculosis" in text
        or "mycobacterium_tuberculosis" in text
        or "mycobacterium tuberculosis complex" in text
        or "m. tuberculosis" in text
        or "结核分枝杆菌" in text
        or "结核杆菌" in text
    )

def _is_mlst_neisseria_meningitidis(mlst_result: dict, checkm_info: dict) -> bool:
    for value in (
        checkm_info.get("species_name"),
        checkm_info.get("mlst_species_name"),
    ):
        if _looks_like_neisseria_meningitidis_text(value):
            return True

    columns = mlst_result.get("columns", []) if isinstance(mlst_result, dict) else []
    rows = mlst_result.get("rows", []) if isinstance(mlst_result, dict) else []
    species_index = columns.index("物种信息") if "物种信息" in columns else -1
    st_index = columns.index("序列分型(ST)") if "序列分型(ST)" in columns else -1
    for row in rows if isinstance(rows, list) else []:
        if species_index >= 0 and species_index < len(row) and _looks_like_neisseria_meningitidis_text(row[species_index]):
            return True
        if st_index >= 0 and st_index < len(row):
            st_text = _normalize_mlst_lookup_text(row[st_index])
            if st_text == "ST4821":
                return True
    return False

def _build_neisseria_amr_summary(call_rows: list[dict[str, str]]) -> dict:
    positive = [row for row in call_rows if str(row.get("present") or "").strip().lower() == "yes"]
    review = [
        row for row in call_rows
        if str(row.get("status") or "").strip().lower() in {"manual_review_recommended", "position_not_covered", "gene_not_found_in_faa"}
    ]
    highlights: list[str] = []
    if positive:
        grouped: dict[str, list[str]] = {}
        for row in positive:
            gene = str(row.get("gene") or "-").strip()
            grouped.setdefault(gene, []).append(str(row.get("protein_change") or row.get("nucleotide_change") or "-").strip())
        for gene, changes in grouped.items():
            uniq = []
            for change in changes:
                if change and change not in uniq:
                    uniq.append(change)
            highlights.append(f"{gene}: {', '.join(uniq)}")
    headline = "未检出明确命中的脑膜炎奈瑟菌已知耐药位点。"
    if positive:
        headline = f"检出 {len(positive)} 个脑膜炎奈瑟菌已知耐药相关位点。"
    elif review:
        headline = "未检出明确阳性位点，但部分目标基因需要人工复核或序列补充确认。"
    positive_keys = {
        (
            str(row.get("gene") or "").strip().lower(),
            str(row.get("protein_change") or row.get("nucleotide_change") or "").strip(),
        )
        for row in positive
    }
    review_keys = {
        (
            str(row.get("gene") or "").strip().lower(),
            str(row.get("protein_change") or row.get("nucleotide_change") or "").strip(),
        )
        for row in review
    }
    interpretation_items: list[str] = []
    classic_peni = {
        ("pena", "F504L"),
        ("pena", "A510V"),
        ("pena", "I515V"),
        ("pena", "H541N"),
        ("pena", "I566V"),
    }
    if classic_peni.issubset(positive_keys):
        interpretation_items.append("命中经典 penA 五位点组合，提示青霉素敏感性下降。")
    if ("pena", "A549T") in positive_keys:
        interpretation_items.append("额外命中 penA A549T，进一步支持青霉素非敏感背景。")
    if ("pena", "N512Y") in positive_keys:
        interpretation_items.append("命中 penA N512Y，提示第三代头孢敏感性下降风险，建议结合完整 penA 单倍型或药敏结果进一步判读。")
    if ("gyra", "T91I") in positive_keys or ("gyra", "T91F") in positive_keys:
        interpretation_items.append("命中 gyrA T91 位点耐药突变，支持环丙沙星耐药。")
    elif any(key in positive_keys for key in {("gyra", "D95A"), ("gyra", "D95G"), ("gyra", "D95N"), ("gyra", "D95Y")}):
        interpretation_items.append("命中 gyrA D95 位点变异，提示喹诺酮敏感性下降或耐药风险。")
    if any(key in positive_keys for key in {("rpob", "H552Y"), ("rpob", "H552N"), ("rpob", "H552R"), ("rpob", "S548F"), ("rpob", "S557F"), ("rpob", "D542N")}):
        interpretation_items.append("命中 rpoB 利福平耐药热点位点，提示利福平耐药风险。")
    if any(key in positive_keys for key in {("folp", "F31L"), ("folp", "G194C"), ("folp", "R228S"), ("folp", "195_196insGS")}):
        interpretation_items.append("命中 folP 磺胺类相关位点，提示磺胺类耐药风险。")
    if ("folp", "195_196insGS") in review_keys:
        interpretation_items.append("folP 195_196insGS 当前为待确认状态，需要结合 CDS 核酸或蛋白比对复核；若确认阳性，提示磺胺类耐药。")
    if not interpretation_items and review:
        interpretation_items.append("当前主要为待确认位点，建议结合更完整序列或药敏试验进行复核。")
    return {
        "headline": headline,
        "highlights": highlights,
        "positive_count": len(positive),
        "review_count": len(review),
        "interpretation_items": interpretation_items,
    }

def _translate_neisseria_amr_drug_class(value: str) -> str:
    mapping = {
        "beta_lactam": "β-内酰胺类",
        "fluoroquinolone": "氟喹诺酮类",
        "rifamycin": "利福平类",
        "sulfonamide": "磺胺类",
        "folate_pathway": "叶酸通路相关",
        "macrolide": "大环内酯类",
    }
    return mapping.get(str(value or "").strip().lower(), str(value or "").strip() or "-")

def _translate_neisseria_amr_drug(value: str) -> str:
    mapping = {
        "penicillin": "青霉素",
        "ampicillin": "氨苄西林",
        "ceftriaxone_or_cefotaxime": "头孢曲松/头孢噻肟",
        "cefotaxime": "头孢噻肟",
        "ceftriaxone": "头孢曲松",
        "ciprofloxacin": "环丙沙星",
        "levofloxacin": "左氧氟沙星",
        "rifampin": "利福平",
        "sulfamethoxazole": "磺胺甲噁唑",
        "sulfonamide": "磺胺类",
        "azithromycin": "阿奇霉素",
    }
    key = str(value or "").strip().lower()
    return mapping.get(key, str(value or "").strip() or "-")

def _translate_neisseria_amr_association(value: str) -> str:
    mapping = {
        "reduced_susceptibility": "敏感性下降",
        "reduced_susceptibility_or_resistance": "敏感性下降或耐药",
        "resistance": "耐药",
    }
    return mapping.get(str(value or "").strip().lower(), str(value or "").strip() or "-")

def _translate_neisseria_amr_note(value: str) -> str:
    text = str(value or "").strip()
    if not text:
        return "-"
    replacements = [
        ("Part of the classic PenI-associated PBP2 pattern; usually interpreted together with additional penA substitutions.", "属于经典 PenI 相关 PBP2 组合位点，通常需要结合其他 penA 位点共同判读。"),
        ("Highlighted in Chen Mingliang team papers as an additional penA substitution associated with penicillin nonsusceptibility, including nonclassic penA184-like patterns.", "陈明亮团队文献强调该位点可作为青霉素非敏感的补充 penA 位点，也见于非经典 penA184 样模式。"),
        ("Reported among penA changes linked to reduced third-generation cephalosporin susceptibility; interpret in penA haplotype context rather than as a standalone marker.", "文献将其列为与第三代头孢敏感性下降相关的 penA 位点，更适合放在 penA 单倍型背景下整体判读。"),
        ("Observed in penA795/FC428-like cefotaxime-resistant backgrounds described in recent Chinese meningococcal reports.", "该位点见于近期中国文献报道的 penA795/FC428-like 头孢噻肟耐药背景。"),
        ("Added from Chen Mingliang team's cefotaxime-resistant reports; part of FC428-like or NEIS1753_5058-associated cephalosporin-resistance combinations.", "该位点来自陈明亮团队头孢噻肟耐药文献，属于 FC428-like 或 NEIS1753_5058 相关头孢耐药组合的一部分。"),
        ("Core ciprofloxacin resistance-associated QRDR substitution in Neisseria meningitidis.", "脑膜炎奈瑟菌喹诺酮耐药最核心的 QRDR 位点之一。"),
        ("Often reported with elevated ciprofloxacin MIC; may occur with other QRDR substitutions.", "常与环丙沙星 MIC 升高相关，也可能与其他 QRDR 位点共同出现。"),
        ("Typically considered an accessory QRDR substitution; interpret together with gyrA.", "通常属于辅助性 QRDR 位点，建议与 gyrA 位点联合判读。"),
        ("Rifampin resistance hotspot in meningococci.", "脑膜炎奈瑟菌利福平耐药热点位点。"),
        ("Classic sulfonamide-associated folP substitution.", "经典磺胺类相关 folP 位点。"),
        ("Candidate macrolide-associated regulator change; evidence in meningococci is limited.", "候选性大环内酯相关调控位点，但在脑膜炎奈瑟菌中的证据仍有限。"),
        ("A Gly-Ser insertion between codons 195 and 196 is linked to elevated sulfonamide MIC.", "195-196 位密码子之间的 Gly-Ser 插入与磺胺类 MIC 升高相关。"),
        ("This target is recorded as a nucleotide-level change (195_196insGS); direct confirmation needs CDS nucleotide extraction/alignment.", "该位点按核酸变化记录，需结合 CDS 核酸提取或序列比对进一步确认。"),
        ("Insertion-style changes should be confirmed on the CDS nucleotide alignment.", "插入类变化建议结合 CDS 核酸比对进一步确认。"),
        ("Nucleotide-level targets require CDS-level nucleotide alignment.", "核酸层面的位点需要结合 CDS 核酸比对判定。"),
        ("Selected protein is too short to cover amino-acid position ", "当前选中的蛋白序列长度不足，无法覆盖氨基酸位点 "),
        ("Multiple candidate proteins found for ", "检测到多个候选蛋白："),
    ]
    for source, target in replacements:
        text = text.replace(source, target)
    return text

def _build_neisseria_amr_section(report_dir: Path, sample_name: str, mlst_result: dict, checkm_info: dict, project_root: Path) -> dict:
    if not _is_mlst_neisseria_meningitidis(mlst_result, checkm_info):
        return {"status": "empty", "headline": "", "columns": [], "rows": [], "source_label": ""}
    result_path = report_dir / f"{sample_name}.neisseria_amr_calls.csv"
    if not result_path.is_file():
        return {
            "status": "empty",
            "headline": "MLST/物种结果提示为脑膜炎奈瑟菌，但主分析流程尚未产出耐药位点识别结果文件。",
            "columns": [],
            "rows": [],
            "source_label": "脑膜炎奈瑟菌耐药位点库（文献整理版）",
        }
    raw_calls: list[dict[str, str]] = []
    try:
        with result_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            raw_calls = [{str(key): str(value or "") for key, value in row.items()} for row in reader]
    except OSError:
        raw_calls = []
    if not raw_calls:
        return {
            "status": "empty",
            "headline": "脑膜炎奈瑟菌耐药位点结果文件存在，但暂未读取到有效记录。",
            "columns": [],
            "rows": [],
            "source_label": "脑膜炎奈瑟菌耐药位点库（文献整理版）",
        }

    columns = ["基因", "位点", "药物/药物类别", "关联解释", "判定", "观测氨基酸", "状态", "置信度", "说明"]
    rows: list[list[str]] = []
    status_label_map = {
        "called": "已判定",
        "manual_review_recommended": "建议复核",
        "position_not_covered": "位点未覆盖",
        "gene_not_found_in_faa": "目标基因未找到",
        "unsupported_target_format": "规则格式待补充",
    }
    for call in raw_calls:
        gene = str(call.get("gene") or "").strip()
        site = str(call.get("protein_change") or call.get("nucleotide_change") or "-").strip()
        drug = " / ".join([
            item for item in [
                _translate_neisseria_amr_drug(str(call.get("drug") or "").strip()),
                _translate_neisseria_amr_drug_class(str(call.get("drug_class") or "").strip()),
            ] if item and item != "-"
        ])
        present = str(call.get("present") or "").strip().lower()
        status = str(call.get("status") or "").strip()
        if present == "no":
            continue
        if present == "yes":
            verdict = "命中"
        elif present == "no":
            verdict = "未命中"
        else:
            verdict = "待确认"
        rows.append([
            gene or "-",
            site or "-",
            drug or "-",
            _translate_neisseria_amr_association(str(call.get("association") or "-").strip()),
            verdict,
            str(call.get("observed_aa") or "-").strip(),
            status_label_map.get(status, status or "-"),
            str(call.get("confidence") or "-").strip(),
            _translate_neisseria_amr_note(str(call.get("notes") or "-").strip()),
        ])

    summary = _build_neisseria_amr_summary(raw_calls)
    return {
        "status": "ready",
        "headline": summary.get("headline", ""),
        "highlights": summary.get("highlights", []),
        "interpretation_items": summary.get("interpretation_items", []),
        "positive_count": summary.get("positive_count", 0),
        "review_count": summary.get("review_count", 0),
        "columns": columns,
        "rows": rows,
        "source_label": "脑膜炎奈瑟菌耐药位点库（文献整理版）",
    }

def _read_tb_summary(report_dir: Path) -> dict:
    path = report_dir / "tb_analysis" / "tb_summary.json"
    if not path.is_file():
        return {}
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    return payload if isinstance(payload, dict) else {}

def _translate_tb_drug_name(value: object) -> str:
    text = str(value or "").strip()
    if not text:
        return "-"
    return TB_DRUG_NAME_ZH.get(text, text)

def _classify_tb_resistance_grade(focus_calls: list[dict[str, Any]]) -> dict[str, Any]:
    resistant_drugs = {
        str(item.get("drug") or "").strip()
        for item in focus_calls
        if isinstance(item, dict) and str(item.get("drug") or "").strip()
    }
    first_line_hits = {
        drug for drug in resistant_drugs
        if drug in {"Isoniazid", "Rifampicin", "Ethambutol", "Pyrazinamide", "Streptomycin"}
    }
    fluoroquinolone_hits = {
        drug for drug in resistant_drugs
        if drug in {"Levofloxacin", "Moxifloxacin"}
    }
    group_a_hits = {
        drug for drug in resistant_drugs
        if drug in {"Bedaquiline", "Linezolid"}
    }
    has_rif = "Rifampicin" in resistant_drugs
    has_inh = "Isoniazid" in resistant_drugs

    if has_rif and fluoroquinolone_hits and group_a_hits:
        label = "XDR-TB"
        tone = "high"
        reason = "已命中利福平耐药，并同时命中氟喹诺酮类和贝达喹啉/利奈唑胺相关重点耐药证据。"
    elif has_rif and fluoroquinolone_hits:
        label = "Pre-XDR-TB"
        tone = "high"
        reason = "已命中利福平耐药，并同时命中左氧氟沙星/莫西沙星等氟喹诺酮类重点耐药证据。"
    elif has_rif and has_inh:
        label = "MDR-TB"
        tone = "high"
        reason = "已同时命中利福平和异烟肼重点耐药证据。"
    elif has_rif:
        label = "RR-TB"
        tone = "watch"
        reason = "已命中利福平重点耐药证据，可按利福平耐药相关结核重点关注。"
    elif resistant_drugs:
        label = "DR-TB"
        tone = "watch"
        if first_line_hits:
            reason = f"已命中非利福平重点耐药证据，主要涉及{ '、'.join(_translate_tb_drug_name(item) for item in sorted(first_line_hits)) }。"
        else:
            reason = "已命中重点耐药证据，但未形成 RR/MDR/Pre-XDR/XDR 的更高分级。"
    else:
        label = "敏感"
        tone = "safe"
        reason = "当前未检出 WHO 1/2 级重点耐药证据，暂不支持耐药分级升级。"

    return {
        "label": label,
        "tone": tone,
        "reason": reason,
        "resistant_drugs": sorted(resistant_drugs),
        "fluoroquinolone_hits": sorted(fluoroquinolone_hits),
        "group_a_hits": sorted(group_a_hits),
    }

def _build_tb_amr_section(report_dir: Path, sample_name: str, checkm_info: dict) -> dict:
    if not any(
        _looks_like_tuberculosis_text(checkm_info.get(key))
        for key in ("species_name", "mlst_species_name")
    ):
        return {"status": "empty", "headline": "", "columns": [], "rows": [], "source_label": ""}
    payload = _read_tb_summary(report_dir)
    if not payload:
        return {
            "status": "empty",
            "headline": "物种结果提示为结核分枝杆菌，但当前尚未产出 H37Rv 有参 SNP 与 mutation catalogue 判读文件。",
            "columns": [],
            "rows": [],
            "source_label": "WHO mutation catalogue 2023 + H37Rv 参考",
        }
    catalogue = payload.get("catalogue") if isinstance(payload.get("catalogue"), dict) else {}
    summary = catalogue.get("summary") if isinstance(catalogue.get("summary"), dict) else {}
    raw_rows = catalogue.get("matched_rows") if isinstance(catalogue.get("matched_rows"), list) else []
    columns = ["药物", "药物分层", "基因", "Catalogue 突变", "样本匹配层级", "样本核酸变化", "样本氨基酸变化", "判读结论", "最终分级", "说明"]
    rows = [
        [
            _translate_tb_drug_name(item.get("药物")),
            str(item.get("药物分层") or "-"),
            str(item.get("基因") or "-"),
            str(item.get("catalogue突变") or "-"),
            str(item.get("样本匹配层级") or "-"),
            str(item.get("样本碱基变化") or "-"),
            str(item.get("样本氨基酸变化") or "-"),
            str(item.get("判读结论") or "-"),
            str(item.get("最终分级") or "-"),
            str(item.get("注释") or "-"),
        ]
        for item in raw_rows if isinstance(item, dict)
    ]
    drug_calls = summary.get("drug_calls") if isinstance(summary.get("drug_calls"), list) else []
    highlights = [
        f"{_translate_tb_drug_name(item.get('drug'))}：{str(item.get('verdict') or '-')}（{', '.join([str(x) for x in (item.get('mutations') or []) if str(x).strip()][:3])}）"
        for item in drug_calls if isinstance(item, dict)
    ]
    focus_calls = [
        {
            **item,
            "drug": _translate_tb_drug_name(item.get("drug")),
            "drug_en": str(item.get("drug") or "").strip(),
            "evidence_grade": str(item.get("grade") or "").strip() or "-",
        }
        for item in (summary.get("focus_drug_calls") if isinstance(summary.get("focus_drug_calls"), list) else [])
        if isinstance(item, dict)
    ]
    other_calls = [
        {
            **item,
            "drug": _translate_tb_drug_name(item.get("drug")),
            "drug_en": str(item.get("drug") or "").strip(),
            "evidence_grade": str(item.get("grade") or "").strip() or "-",
        }
        for item in (summary.get("other_drug_calls") if isinstance(summary.get("other_drug_calls"), list) else [])
        if isinstance(item, dict)
    ]
    tb_resistance_grade = _classify_tb_resistance_grade(
        [item for item in (summary.get("focus_drug_calls") if isinstance(summary.get("focus_drug_calls"), list) else []) if isinstance(item, dict)]
    )
    interpretation_items = [
        str(item).strip()
        for item in (summary.get("interpretation_items") or [])
        if str(item).strip()
    ]
    for en_name, zh_name in TB_DRUG_NAME_ZH.items():
        interpretation_items = [item.replace(en_name, zh_name) for item in interpretation_items]
    interpretation_items.insert(0, f"基于 WHO 1/2 级重点耐药证据，当前样本可归为 {tb_resistance_grade['label']}：{tb_resistance_grade['reason']}")
    focus_count = int(summary.get("focus_drug_count") or 0)
    total_drug_count = int(summary.get("total_drug_count") or 0)
    matched_variant_count = int(summary.get("matched_variant_count") or 0)
    grade_label = str(tb_resistance_grade.get("label") or "").strip()
    grade_reason = str(tb_resistance_grade.get("reason") or "").strip()
    headline = f"基于 H37Rv 有参 SNP 与 WHO mutation catalogue，当前样本结核耐药分级为 {grade_label}"
    if grade_reason:
        headline += f"；{grade_reason}"
    if focus_count:
        headline += f" 共识别到 {focus_count} 个药物存在 WHO 1/2 级重点耐药证据。"
    else:
        headline += "，当前未识别到 WHO 1/2 级重点耐药证据。"
    return {
        "status": "ready",
        "headline": headline,
        "highlights": highlights,
        "interpretation_items": interpretation_items,
        "positive_count": focus_count,
        "review_count": total_drug_count,
        "focus_calls": focus_calls,
        "other_calls": other_calls,
        "matched_variant_count": matched_variant_count,
        "resistance_grade": tb_resistance_grade,
        "columns": columns,
        "rows": rows,
        "source_label": "WHO mutation catalogue 2023（version 2）+ H37Rv 有参 SNP 判读",
        "organism_label": "结核分枝杆菌",
        "title": "结核分枝杆菌耐药位点判读",
        "tag_label": "WHO catalogue",
    }

def _build_tb_typing_mlst_section(report_dir: Path, sample_name: str, checkm_info: dict) -> dict | None:
    if not any(
        _looks_like_tuberculosis_text(checkm_info.get(key))
        for key in ("species_name", "mlst_species_name")
    ):
        return None
    payload = _read_tb_summary(report_dir)
    if not payload:
        return {
            "columns": ["项目", "结果", "说明"],
            "rows": [],
            "gene_show_map": {},
            "default_gene": "",
            "knowledge_summary": {
                "headline": "已识别为结核分枝杆菌，但当前尚未读取到 tb-profiler 分型摘要文件。",
                "items": [],
            },
            "title": "结核分枝杆菌分型摘要",
            "tag_label": "tb-profiler",
            "empty_message": "尚未生成可展示的结核分型摘要。",
            "detail_empty_message": "当前没有附加的位点级详情可供展开。",
            "generic_detail_note": "当前展示的是基于 tb-profiler 汇总的结核分型摘要。",
        }
    tbprofiler = payload.get("tbprofiler") if isinstance(payload.get("tbprofiler"), dict) else {}
    lineage = str(tbprofiler.get("predicted_lineage") or "-").strip() or "-"
    family = str(tbprofiler.get("family") or "-").strip() or "-"
    spoligotype = str(tbprofiler.get("spoligotype") or "-").strip() or "-"
    drug_type = str(tbprofiler.get("drug_type") or "-").strip() or "-"
    rows = [
        ["Lineage", lineage, "tb-profiler 主家系判定"],
        ["Family / Spoligotype", family, "对应 spoligotype 家族归属"],
        ["Spoligotype", spoligotype, "若数据库未返回则显示 -"],
        ["DR type", drug_type, "tb-profiler 内置耐药类别标签，仅作补充展示"],
    ]
    summary_items = [
        {
            "lineage_text": lineage if lineage != "-" else "",
            "virulence": [],
            "resistance": [f"tb-profiler DR type：{drug_type}"] if drug_type != "-" else [],
            "regional": [family] if family != "-" else [],
            "interpretation": "结核样本不使用常规 MLST 管家基因分型，这里改为展示 tb-profiler 的 lineage / family 摘要。",
        }
    ]
    return {
        "columns": ["项目", "结果", "说明"],
        "rows": rows,
        "gene_show_map": {},
        "default_gene": "",
        "knowledge_summary": {
            "headline": "已切换为结核专用分型摘要，直接展示 tb-profiler 的家系与家族结果。",
            "items": summary_items,
        },
        "title": "结核分枝杆菌分型摘要",
        "tag_label": "tb-profiler",
        "empty_message": "未读取到可展示的结核分型结果。",
        "detail_empty_message": "当前没有附加的位点级详情可供展开。",
        "generic_detail_note": "当前展示的是基于 tb-profiler 汇总的结核分型摘要。",
    }

def _build_tb_serotype_section(report_dir: Path, sample_name: str, checkm_info: dict) -> dict | None:
    if not any(
        _looks_like_tuberculosis_text(checkm_info.get(key))
        for key in ("species_name", "mlst_species_name")
    ):
        return None
    payload = _read_tb_summary(report_dir)
    if not payload:
        return {
            "status": "empty",
            "mode": "tb_profiler",
            "columns": [],
            "rows": [],
            "summary_cards": [],
            "notes": "已识别为结核分枝杆菌，但尚未读取到 tb-profiler 家系分析结果。",
        }
    tbprofiler = payload.get("tbprofiler") if isinstance(payload.get("tbprofiler"), dict) else {}
    lineage = str(tbprofiler.get("predicted_lineage") or "-").strip() or "-"
    family = str(tbprofiler.get("family") or "-").strip() or "-"
    spoligotype = str(tbprofiler.get("spoligotype") or "-").strip() or "-"
    raw_lineage_items = tbprofiler.get("raw_lineage_items") if isinstance(tbprofiler.get("raw_lineage_items"), list) else []
    main_lineage = "-"
    sub_lineage = lineage
    main_family = "-"
    rd_value = "-"
    if raw_lineage_items:
        lineage_values = [
            str(item.get("lineage") or "").strip()
            for item in raw_lineage_items
            if isinstance(item, dict) and str(item.get("lineage") or "").strip()
        ]
        family_values = [
            str(item.get("family") or "").strip()
            for item in raw_lineage_items
            if isinstance(item, dict) and str(item.get("family") or "").strip()
        ]
        rd_values = [
            str(item.get("rd") or "").strip()
            for item in raw_lineage_items
            if isinstance(item, dict) and str(item.get("rd") or "").strip()
        ]
        if lineage_values:
            main_lineage = lineage_values[0]
            sub_lineage = lineage_values[-1]
        if family_values:
            main_family = family_values[0]
        if rd_values:
            rd_value = rd_values[-1]
    if main_lineage == "-" and lineage != "-":
        main_lineage = lineage.split(".", 1)[0] if "." in lineage else lineage
    if main_family == "-" and family != "-":
        main_family = family
    if spoligotype == "-" and main_family != "-":
        spoligotype = main_family
    knowledge_summary = _build_tb_knowledge_summary(
        str(Path(__file__).resolve().parent.parent),
        checkm_info.get("species_name") or checkm_info.get("mlst_species_name") or "Mycobacterium tuberculosis",
        [
            "MTBC",
            main_lineage,
            sub_lineage,
            lineage,
            family,
            spoligotype,
            rd_value,
        ],
    )
    summary_cards = [
        {"label": "主家系", "value": main_lineage},
        {"label": "子家系", "value": sub_lineage},
        {"label": "Spoligotype", "value": spoligotype},
        {"label": "RD", "value": rd_value},
    ]
    return {
        "status": "ready",
        "mode": "tb_profiler",
        "predicted_lineage": lineage,
        "predicted_serotype": family,
        "predicted_group": "MTB",
        "summary_cards": summary_cards,
        "notes": "基于 TB 环境中的 tb-profiler 完成结核分枝杆菌家系分析；耐药结论单独基于 H37Rv 有参 SNP 与 WHO mutation catalogue 组织展示。",
        "knowledge_summary": knowledge_summary,
        "columns": ["项目", "结果"],
        "rows": [
            ["主家系", main_lineage],
            ["子家系", sub_lineage],
            ["Spoligotype", spoligotype],
            ["RD", rd_value],
        ],
    }
