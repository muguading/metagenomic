#!/usr/bin/env python3
"""Classify TB resistance grade from WHO catalogue drug calls.

This is a standalone copy of the report-page logic used for TB resistance
grading. It intentionally depends only on the Python standard library so it can
be moved into another project.

Expected input can be one of:
1. A metagenomic tb_summary.json containing catalogue.summary.focus_drug_calls.
2. A JSON list of drug call objects, each with a "drug" field.
3. A JSON object with "focus_drug_calls" or "drugs".

Only WHO 1/2 resistance evidence should be passed in. When using tb_summary.json,
that filtering has already been done by the postprocess script.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any


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

ZH_TO_TB_DRUG_NAME = {value: key for key, value in TB_DRUG_NAME_ZH.items()}

FIRST_LINE_DRUGS = {
    "Isoniazid",
    "Rifampicin",
    "Ethambutol",
    "Pyrazinamide",
    "Streptomycin",
}
FLUOROQUINOLONE_DRUGS = {"Levofloxacin", "Moxifloxacin"}
GROUP_A_DRUGS = {"Bedaquiline", "Linezolid"}
FOCUS_GRADES = {"1) Assoc w R", "2) Assoc w R - Interim"}


def normalize_drug_name(value: object) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    return ZH_TO_TB_DRUG_NAME.get(text, text)


def translate_drug_name(value: object) -> str:
    text = str(value or "").strip()
    if not text:
        return "-"
    return TB_DRUG_NAME_ZH.get(text, text)


def extract_focus_calls(payload: Any) -> list[dict[str, Any]]:
    """Extract WHO 1/2 focus drug calls from supported JSON shapes."""
    if isinstance(payload, list):
        return [item for item in payload if isinstance(item, dict)]

    if not isinstance(payload, dict):
        return []

    catalogue = payload.get("catalogue")
    if isinstance(catalogue, dict):
        summary = catalogue.get("summary")
        if isinstance(summary, dict) and isinstance(summary.get("focus_drug_calls"), list):
            return [item for item in summary["focus_drug_calls"] if isinstance(item, dict)]

    for key in ("focus_drug_calls", "drugs", "drug_calls"):
        calls = payload.get(key)
        if isinstance(calls, list):
            return [item for item in calls if isinstance(item, dict)]

    return []


def filter_focus_calls(calls: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Keep only WHO 1/2 resistance evidence when grade fields are present."""
    filtered = []
    for item in calls:
        grade = str(item.get("grade") or item.get("最终分级") or "").strip()
        if not grade or grade in FOCUS_GRADES:
            filtered.append(item)
    return filtered


def classify_tb_resistance_grade(focus_calls: list[dict[str, Any]]) -> dict[str, Any]:
    resistant_drugs = {
        normalize_drug_name(item.get("drug") or item.get("药物"))
        for item in focus_calls
        if isinstance(item, dict)
    }
    resistant_drugs.discard("")

    first_line_hits = sorted(drug for drug in resistant_drugs if drug in FIRST_LINE_DRUGS)
    fluoroquinolone_hits = sorted(
        drug for drug in resistant_drugs if drug in FLUOROQUINOLONE_DRUGS
    )
    group_a_hits = sorted(drug for drug in resistant_drugs if drug in GROUP_A_DRUGS)

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
            reason = (
                "已命中非利福平重点耐药证据，主要涉及"
                + "、".join(translate_drug_name(drug) for drug in first_line_hits)
                + "。"
            )
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
        "resistant_drugs_zh": [translate_drug_name(drug) for drug in sorted(resistant_drugs)],
        "fluoroquinolone_hits": fluoroquinolone_hits,
        "fluoroquinolone_hits_zh": [
            translate_drug_name(drug) for drug in fluoroquinolone_hits
        ],
        "group_a_hits": group_a_hits,
        "group_a_hits_zh": [translate_drug_name(drug) for drug in group_a_hits],
    }


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Classify TB resistance grade from WHO 1/2 catalogue drug calls."
    )
    parser.add_argument("input_json", type=Path, help="tb_summary.json or drug-call JSON")
    parser.add_argument(
        "--no-grade-filter",
        action="store_true",
        help="Do not filter calls by grade/最终分级 before classification.",
    )
    parser.add_argument(
        "--pretty",
        action="store_true",
        help="Pretty-print JSON output.",
    )
    args = parser.parse_args()

    try:
        payload = load_json(args.input_json)
    except (OSError, json.JSONDecodeError) as exc:
        print(f"failed to read JSON: {exc}", file=sys.stderr)
        return 2

    focus_calls = extract_focus_calls(payload)
    if not args.no_grade_filter:
        focus_calls = filter_focus_calls(focus_calls)

    result = classify_tb_resistance_grade(focus_calls)
    result["focus_call_count"] = len(focus_calls)

    json.dump(
        result,
        sys.stdout,
        ensure_ascii=False,
        indent=2 if args.pretty else None,
        sort_keys=True,
    )
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
