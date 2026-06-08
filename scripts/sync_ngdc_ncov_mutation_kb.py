#!/usr/bin/env python3
from __future__ import annotations

import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable
from urllib.request import Request, urlopen


PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = PROJECT_ROOT / "database" / "virus" / "ncov" / "ngdc_mutation_kb"
RAW_DIR = OUTPUT_DIR / "raw"
USER_AGENT = "metagenomic-portal/1.0 (+https://ngdc.cncb.ac.cn/ncov/knowledge/mutation)"
ORF1A_AA_LENGTH = 4401

ENDPOINTS = {
    "mutation_effect_transmission": {
        "url": "https://ngdc.cncb.ac.cn/ncov/rest/mutationEffect/list/1",
        "section_key": "transmission",
        "section_label": "传播/感染力",
    },
    "mutation_effect_antibody": {
        "url": "https://ngdc.cncb.ac.cn/ncov/rest/mutationEffect/list/2",
        "section_key": "antibody",
        "section_label": "抗体逃逸",
    },
    "mutation_effect_drug": {
        "url": "https://ngdc.cncb.ac.cn/ncov/rest/mutationEffect/list/3",
        "section_key": "drug",
        "section_label": "药物耐药",
    },
    "mutation_effect_tcell": {
        "url": "https://ngdc.cncb.ac.cn/ncov/rest/mutationEffect/list/4",
        "section_key": "tcell",
        "section_label": "T 细胞表位",
    },
    "pathogenicity_prediction": {
        "url": "https://ngdc.cncb.ac.cn/ncov/api/kg/path/list",
        "section_key": "pathogenicity",
        "section_label": "致病性/分子机制",
    },
    "structure_stability": {
        "url": "https://ngdc.cncb.ac.cn/ncov/api/kg/structure/list",
        "section_key": "structure",
        "section_label": "结构稳定性",
    },
    "mutation_overview": {
        "url": "https://ngdc.cncb.ac.cn/ncov/api/kg/overview",
        "section_key": "overview",
        "section_label": "知识概览",
    },
}

EFFECT_TRANSLATIONS = {
    "Increase": "增加",
    "Decrease": "降低",
    "Yes": "是",
    "No": "否",
    "NA": "无明确结论",
    "Drug resistance": "耐药",
    "Destabilizing": "降低稳定性",
    "Stabilizing": "提高稳定性",
    "SNP": "单核苷酸变异",
    "MutPred2": "MutPred2 预测",
    "CD4": "CD4",
    "CD8": "CD8",
}

METHOD_TRANSLATIONS = {
    "Pseudovirus entry assay": "假病毒入侵实验",
    "Pseudovirus neutralization assays": "假病毒中和实验",
    "Structure model and clinical trial": "结构模型与临床试验",
    "In vitro expriment": "体外实验",
    "in vitro expriment": "体外实验",
    "Molecular dynamics simulation": "分子动力学模拟",
    "Predicted by SAMMBE-3D": "SAMMBE-3D 预测",
    "NMR spectroscopy": "核磁共振实验",
    "Cell culture": "细胞培养实验",
    "Cell culture selection": "细胞培养筛选",
    "Experiment_verfied": "实验验证",
    "IEDB database": "IEDB 数据库",
    "Structure model": "结构模型",
    "Clinical trial": "临床试验",
}

DETAIL_REPLACEMENTS = [
    ("Increased viral entry", "病毒入侵能力增加"),
    ("Significantly resistant to neutralization", "对中和作用显著耐受"),
    ("Increase neutralization escape", "中和逃逸能力增加"),
    ("Affect neutralization by mAb", "影响单克隆抗体中和"),
    ("Resistant to Remdesivir", "对瑞德西韦耐药"),
    ("Mutation associated with potential small molecule inhibitor resistance", "与潜在小分子抑制剂耐药相关"),
    ("The mutation increased the production of viruslike particles", "该突变增加了类病毒颗粒产生"),
    ("abolished HMA binding", "消除了 HMA 结合能力"),
]


def fetch_json(url: str) -> list[dict]:
    request = Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
    with urlopen(request, timeout=60) as response:
        payload = json.loads(response.read().decode("utf-8", errors="ignore"))
    return payload if isinstance(payload, list) else []


def normalize_gene_name(gene: object) -> str:
    text = str(gene or "").strip()
    if not text:
        return ""
    normalized = text.lower().replace(" ", "")
    mapping = {
        "orf1ab": "ORF1ab",
        "orf1a": "ORF1a",
        "orf1b": "ORF1b",
        "orf3a": "ORF3a",
        "orf6": "ORF6",
        "orf7a": "ORF7a",
        "orf7b": "ORF7b",
        "orf8": "ORF8",
        "orf9b": "ORF9b",
        "orf10": "ORF10",
        "spike": "S",
    }
    return mapping.get(normalized, text)


def replace_first_number(text: str, number: int) -> str:
    return re.sub(r"\d+", str(number), text, count=1)


def build_aa_aliases(gene: object, mutation: object) -> list[str]:
    gene_name = normalize_gene_name(gene)
    change = str(mutation or "").strip()
    if not gene_name or not change:
        return []
    aliases = {f"{gene_name}:{change}"}
    match = re.search(r"(\d+)", change)
    if not match:
        return sorted(aliases)
    position = int(match.group(1))
    if gene_name == "ORF1ab":
        if position <= ORF1A_AA_LENGTH:
            aliases.add(f"ORF1a:{replace_first_number(change, position)}")
        else:
            aliases.add(f"ORF1b:{replace_first_number(change, position - ORF1A_AA_LENGTH)}")
    elif gene_name == "ORF1a":
        aliases.add(f"ORF1ab:{replace_first_number(change, position)}")
    elif gene_name == "ORF1b":
        aliases.add(f"ORF1ab:{replace_first_number(change, position + ORF1A_AA_LENGTH)}")
    return sorted(aliases)


def compact_text(value: object) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    return re.sub(r"\s+", " ", text)


def translate_phrase(text: object, mapping: dict[str, str]) -> str:
    raw = compact_text(text)
    if not raw:
        return ""
    return mapping.get(raw, raw)


def translate_detail_text(text: object) -> str:
    raw = compact_text(text)
    if not raw:
        return ""
    translated = raw
    for source, target in DETAIL_REPLACEMENTS:
        translated = translated.replace(source, target)
    translated = translated.replace("May decrease the eff", "可能降低效应")
    translated = translated.replace("Based on the changes of Gibbs free energy", "基于 Gibbs 自由能变化分析")
    translated = translated.replace("It may casue the loss of efficacy for some vaccines.", "可能导致部分疫苗效力下降。")
    translated = translated.replace("It reduces spike sensitivity to neutralization", "会降低刺突蛋白对中和作用的敏感性")
    return translated


def normalize_contents_codes(value: object) -> str:
    mapping = {
        "A": "抗体",
        "I": "感染力",
        "S": "结构稳定性",
        "P": "致病性",
        "T": "T 细胞表位",
        "D": "药物耐药",
    }
    codes = [mapping.get(item.strip(), item.strip()) for item in str(value or "").split(";") if item.strip()]
    deduped: list[str] = []
    for code in codes:
        if code and code not in deduped:
            deduped.append(code)
    return "；".join(deduped)


def normalize_record(dataset_name: str, meta: dict, row: dict) -> dict | None:
    gene = normalize_gene_name(row.get("gene"))
    mutation = compact_text(
        row.get("aminoAcidMutation")
        or row.get("subsitution")
        or row.get("mutation")
    )
    aa_keys = build_aa_aliases(gene, mutation)
    if not gene or not mutation or not aa_keys:
        return None

    section_key = str(meta.get("section_key") or "").strip()
    section_label = str(meta.get("section_label") or "").strip()
    effect = ""
    detail = ""
    method = ""
    pmid = ""
    extra = {}

    if section_key == "transmission":
        effect = compact_text(row.get("transmissionEffect"))
        detail = compact_text(row.get("transmissionDetail"))
        method = compact_text(row.get("transmissionMethod"))
        pmid = compact_text(row.get("transmissionPmid"))
    elif section_key == "antibody":
        effect = compact_text(row.get("antibodyEffect"))
        detail = compact_text(row.get("antibodyDetail"))
        method = compact_text(row.get("antibodyMethod"))
        pmid = compact_text(row.get("antibodyPmid"))
    elif section_key == "drug":
        effect = compact_text(row.get("drugEffect"))
        detail = compact_text(row.get("drugDetail"))
        method = compact_text(row.get("drugMethod"))
        pmid = compact_text(row.get("drugPmid"))
    elif section_key == "tcell":
        effect = compact_text(row.get("tcellEffect"))
        detail = compact_text(row.get("tcellDescription"))
        method = compact_text(row.get("tcellType"))
        pmid = compact_text(row.get("tcellPmid"))
        extra = {
            "epitope": compact_text(row.get("tcellEpitope")),
            "hla": compact_text(row.get("tcellHla")),
            "epitope_position": compact_text(row.get("tcellPosition")),
        }
    elif section_key == "pathogenicity":
        effect = compact_text(row.get("mutPred2Score"))
        detail = compact_text(row.get("molecularMechanisms"))
        method = "MutPred2"
        extra = {
            "motif_information": compact_text(row.get("motifInformation")),
        }
    elif section_key == "structure":
        effect = compact_text(row.get("effect"))
        detail = compact_text(row.get("ddg"))
        method = compact_text(row.get("mutationType"))
    elif section_key == "overview":
        effect = normalize_contents_codes(row.get("contents"))
        detail = ""
        method = compact_text(row.get("resources"))

    return {
        "record_id": row.get("id"),
        "dataset": dataset_name,
        "source": "NGDC",
        "section_key": section_key,
        "section_label": section_label,
        "gene": gene,
        "mutation": mutation,
        "aa_keys": aa_keys,
        "primary_aa_key": aa_keys[0],
        "genomic_position": compact_text(row.get("variant") or row.get("nucPos")),
        "nuc_change": compact_text(row.get("nucChange")),
        "effect": effect,
        "effect_zh": translate_phrase(effect, EFFECT_TRANSLATIONS),
        "detail": detail,
        "detail_zh": translate_detail_text(detail),
        "method": method,
        "method_zh": translate_phrase(method, METHOD_TRANSLATIONS),
        "pmid": pmid,
        "extra": extra,
        "raw": row,
    }


def dedupe_records(records: Iterable[dict]) -> list[dict]:
    seen: set[tuple[str, str, str, str, str]] = set()
    deduped: list[dict] = []
    for record in records:
        key = (
            str(record.get("section_key") or ""),
            str(record.get("gene") or ""),
            str(record.get("mutation") or ""),
            str(record.get("effect") or ""),
            str(record.get("detail") or ""),
        )
        if key in seen:
            continue
        seen.add(key)
        deduped.append(record)
    return deduped


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    RAW_DIR.mkdir(parents=True, exist_ok=True)

    raw_payloads: dict[str, list[dict]] = {}
    normalized_records: list[dict] = []
    counts: dict[str, int] = {}

    for dataset_name, meta in ENDPOINTS.items():
        rows = fetch_json(meta["url"])
        raw_payloads[dataset_name] = rows
        counts[dataset_name] = len(rows)
        (RAW_DIR / f"{dataset_name}.json").write_text(
            json.dumps(rows, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
        for row in rows:
            if not isinstance(row, dict):
                continue
            normalized = normalize_record(dataset_name, meta, row)
            if normalized:
                normalized_records.append(normalized)

    normalized_records = dedupe_records(normalized_records)
    by_aa_key: dict[str, list[dict]] = {}
    for record in normalized_records:
        for aa_key in record.get("aa_keys") or []:
            by_aa_key.setdefault(str(aa_key), []).append(record)

    for aa_key, items in list(by_aa_key.items()):
        by_aa_key[aa_key] = dedupe_records(items)

    payload = {
        "source": "NGDC SARS-CoV-2 mutation knowledge",
        "source_page": "https://ngdc.cncb.ac.cn/ncov/knowledge/mutation",
        "downloaded_at": datetime.now(timezone.utc).isoformat(),
        "endpoint_counts": counts,
        "record_count": len(normalized_records),
        "by_aa_key": by_aa_key,
    }
    (OUTPUT_DIR / "knowledge_index.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(json.dumps({
        "output_dir": str(OUTPUT_DIR),
        "record_count": len(normalized_records),
        "endpoint_counts": counts,
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
