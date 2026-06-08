from __future__ import annotations

import csv
import json
import re
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
HIV_DB = ROOT / "database" / "virus" / "HIV"
REGA_REF_DIR = HIV_DB / "rega_reference_genomes"
OUT_DIR = ROOT / "database" / "knowledge_base" / "typing"
OUT_PATH = OUT_DIR / "hiv_typing.json"

REFERENCE_MANIFEST = REGA_REF_DIR / "reference_manifest.tsv"
SUPPLEMENT_MANIFEST = REGA_REF_DIR / "supplement_manifest.tsv"

HIV1_SPECIES = "Human immunodeficiency virus 1"
HIV2_SPECIES = "Human immunodeficiency virus 2"
HIV1_COMMON = "艾滋病病毒1型"
HIV2_COMMON = "艾滋病病毒2型"
PATHOGEN_ID = "human_immunodeficiency_virus"

PANEL_BROAD = "HIV 大亚型知识库"
PANEL_SUBTYPE = "HIV-1 亚型 / CRF 知识库"

MANUAL_GROUPS = {
    "A6": {
        "aliases": ["A6", "Subtype A6", "subtype A6", "HIV-1 A6", "HIV-1 subtype A6"],
        "subtype_type": "HIV-1 subtype-related lineage",
        "regional_distribution": ["东欧和中亚传播网络中常见，可视作 A1 相关地方性谱系背景。"],
        "interpretation": "A6 常被视作与 A1 相关的地方性谱系，适合用于补充东欧/中亚相关传播背景；结果仍应结合系统发育和重组证据保守解释。",
    }
}

SPECIAL_INTERPRETATIONS = {
    "B": "HIV-1 B 是全球研究和耐药解释中最常见的参考背景之一，在欧美及男男性传播网络中长期占主导；适合作为耐药、传播链和参考株选择的核心分型标签。",
    "C": "HIV-1 C 是全球流行最广的 HIV-1 纯亚型之一，在南部非洲和部分亚洲地区负担突出；命中后适合补充区域流行背景与传播网络解释。",
    "A1": "HIV-1 A1 属于常见纯亚型之一，在东非、东欧及部分跨区域传播网络中较常见；适合作为谱系来源和流行病学关联的补充分层。",
    "A2": "HIV-1 A2 属于较少见但可稳定识别的纯亚型，适合作为 HIV-1 内部谱系差异和输入性背景的补充标签。",
    "D": "HIV-1 D 与较高病毒复制活性和东非流行背景常被共同讨论；结果应结合宿主、治疗史和地区流调资料解释。",
    "F1": "HIV-1 F1 属于较少见但临床和分子流调中可识别的纯亚型，常见于南美和部分欧洲传播网络。",
    "F2": "HIV-1 F2 属于低频纯亚型，更适合作为分子流调和输入性背景记录，而不应单独放大临床严重性。",
    "G": "HIV-1 G 常与西非及相关传播网络背景共同出现，可作为地区传播和重组来源分析的重要标签。",
    "CRF01_AE": "CRF01_AE 是亚洲地区最重要的 HIV-1 循环重组型之一，在中国和东南亚传播网络中极为常见；命中后应优先结合本地流行谱和重组背景解释。",
    "CRF02_AG": "CRF02_AG 是西非和全球输入性样本中最常见的 HIV-1 循环重组型之一，常用于描述跨区域传播和重组来源背景。",
    "CRF07_BC": "CRF07_BC 是中国传播网络中极具代表性的 HIV-1 循环重组型之一，适合用于补充本土传播谱系和分子流调背景。",
    "CRF08_BC": "CRF08_BC 是中国常见的 HIV-1 循环重组型之一，常与西南地区及特定传播网络相关；结果应结合地区流调背景解释。",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def sanitize_id(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", text.strip().lower()).strip("_")


def normalize_group(raw_group: str) -> str:
    text = str(raw_group or "").strip()
    if text.startswith("Subtype "):
        return text.split(" ", 1)[1].strip()
    return text.strip()


def subtype_type_for(group: str) -> str:
    if group.startswith("CRF"):
        return "circulating recombinant form (CRF)"
    if group in MANUAL_GROUPS:
        return str(MANUAL_GROUPS[group].get("subtype_type") or "HIV-1 subtype-related lineage")
    return "HIV-1 subtype"


def aliases_for(group: str) -> list[str]:
    values = [
        group,
        f"HIV-1 {group}",
        f"HIV-1 subtype {group}",
        f"Subtype {group}",
        f"subtype {group}",
    ]
    if group.startswith("CRF"):
        trimmed = group[3:]
        values.extend([f"CRF {trimmed}", group.replace("_", "-"), trimmed, trimmed.replace("_", "-")])
    if group in MANUAL_GROUPS:
        values.extend(MANUAL_GROUPS[group].get("aliases") or [])
    seen: list[str] = []
    for value in values:
        item = str(value or "").strip()
        if item and item not in seen:
            seen.append(item)
    return seen


def interpretation_for(group: str) -> str:
    if group in MANUAL_GROUPS:
        return str(MANUAL_GROUPS[group].get("interpretation") or "")
    if group in SPECIAL_INTERPRETATIONS:
        return SPECIAL_INTERPRETATIONS[group]
    if group.startswith("CRF"):
        return f"{group} 属于 HIV-1 循环重组型（CRF），适合作为重组背景、参考株选择和分子流行病学关联解释依据。"
    return f"HIV-1 {group} 属于可稳定识别的纯亚型，适合作为参考株选择、传播网络和区域流行背景的补充分型标签。"


def regional_distribution_for(group: str) -> list[str]:
    if group in MANUAL_GROUPS:
        return [str(item).strip() for item in (MANUAL_GROUPS[group].get("regional_distribution") or []) if str(item).strip()]
    if group == "CRF01_AE":
        return ["东南亚和中国传播网络中常见。"]
    if group == "CRF02_AG":
        return ["西非及相关输入性传播网络中常见。"]
    if group == "CRF07_BC":
        return ["中国本土传播网络中常见。"]
    if group == "CRF08_BC":
        return ["中国西南相关传播网络中常见。"]
    if group == "B":
        return ["欧美及部分全球传播网络中常见。"]
    if group == "C":
        return ["南部非洲和部分亚洲地区常见。"]
    return []


def build_fasta_index() -> dict[str, str]:
    index: dict[str, str] = {}
    for path in sorted(REGA_REF_DIR.rglob("*.fasta")):
        accession_match = re.search(r"__([A-Z]{1,4}\d{5,8}(?:\.\d+)?)\.fasta$", path.name)
        if not accession_match:
            continue
        accession = accession_match.group(1).split(".", 1)[0]
        index.setdefault(accession, str(path.resolve()))
    return index


def clean_reference(accession: str, title: str, fasta_path: str) -> dict[str, str]:
    return {
        "accession": accession,
        "title": title,
        "fasta_path": fasta_path,
    }


def sort_group_key(group: str) -> tuple[int, str]:
    if group == "A6":
        return (2, group)
    if group.startswith("CRF"):
        return (3, group)
    return (1, group)


def build_entries() -> list[dict[str, object]]:
    fasta_index = build_fasta_index()
    by_group: dict[str, list[dict[str, str]]] = defaultdict(list)

    for row in read_tsv(REFERENCE_MANIFEST):
        if str(row.get("download_status") or "").strip().lower() != "downloaded":
            continue
        accession = str(row.get("note") or "").strip()
        group = normalize_group(str(row.get("group") or ""))
        if not accession or not group:
            continue
        by_group[group].append(
            {
                "accession": accession,
                "title": str(row.get("label") or "").strip(),
                "fasta_path": fasta_index.get(accession, ""),
            }
        )

    for row in read_tsv(SUPPLEMENT_MANIFEST):
        accession = str(row.get("accession") or "").strip()
        group = normalize_group(str(row.get("group") or row.get("subtype") or ""))
        if not accession or not group:
            continue
        by_group[group].append(
            {
                "accession": accession,
                "title": str(row.get("title") or "").strip(),
                "fasta_path": str(row.get("path") or "").strip() or fasta_index.get(accession, ""),
            }
        )

    groups = sorted(set(by_group.keys()) | set(MANUAL_GROUPS.keys()), key=sort_group_key)
    reference_rows = [record for records in by_group.values() for record in records]

    entries: list[dict[str, object]] = [
        {
            "id": "hiv_hiv1",
            "level": "broad_type",
            "broad_type": "HIV-1",
            "typing_source": "REGA-like HIV typing",
            "species": HIV1_SPECIES,
            "common_name": HIV1_COMMON,
            "pathogen_id": PATHOGEN_ID,
            "panel": PANEL_BROAD,
            "serotype": "HIV-1",
            "aliases": ["HIV-1", HIV1_SPECIES, HIV1_COMMON],
            "subtype_panel": PANEL_SUBTYPE,
            "subtypes": groups,
            "reference_count": len(reference_rows),
            "reference_accessions": sorted({row["accession"] for row in reference_rows if row["accession"]}),
            "interpretation": "HIV-1 命中后应继续结合纯亚型、CRF 和重组信号进行分子流调与参考株选择解释；耐药分析建议再回到 HXB2 坐标体系进行判读。",
        },
        {
            "id": "hiv_hiv2",
            "level": "broad_type",
            "broad_type": "HIV-2",
            "typing_source": "Broad coverage typing",
            "species": HIV2_SPECIES,
            "common_name": HIV2_COMMON,
            "pathogen_id": PATHOGEN_ID,
            "panel": PANEL_BROAD,
            "serotype": "HIV-2",
            "aliases": ["HIV-2", HIV2_SPECIES, HIV2_COMMON],
            "reference_count": 0,
            "reference_accessions": [],
            "interpretation": "HIV-2 与 HIV-1 在流行范围、传播效率和耐药解释流程上存在明显差异；在当前本地流程中，应优先停留在 broad typing 层并避免直接套用 HIV-1 子亚型/耐药规则。",
        },
    ]

    for group in groups:
        refs = by_group.get(group, [])
        accessions = []
        seen_accessions: set[str] = set()
        for row in refs:
            accession = str(row.get("accession") or "").strip()
            if accession and accession not in seen_accessions:
                seen_accessions.add(accession)
                accessions.append(accession)
        entries.append(
            {
                "id": f"hiv1_{sanitize_id(group)}",
                "level": "subtype",
                "broad_type": "HIV-1",
                "typing_source": "REGA-like HIV typing",
                "species": HIV1_SPECIES,
                "common_name": HIV1_COMMON,
                "pathogen_id": PATHOGEN_ID,
                "panel": PANEL_SUBTYPE,
                "serotype": group,
                "subtype": group,
                "subtype_type": subtype_type_for(group),
                "aliases": aliases_for(group),
                "reference_count": len(accessions),
                "reference_accessions": accessions,
                "references": [
                    clean_reference(
                        str(row.get("accession") or "").strip(),
                        str(row.get("title") or "").strip(),
                        str(row.get("fasta_path") or "").strip(),
                    )
                    for row in refs
                ],
                "regional_distribution": regional_distribution_for(group),
                "interpretation": interpretation_for(group),
            }
        )

    return entries


def main() -> None:
    payload = {
        "schema_version": "v1",
        "id": "hiv_typing",
        "title": "HIV 大亚型与子亚型知识库",
        "description": "关联 HIV-1 / HIV-2 broad typing、HIV-1 纯亚型与常见 CRF、参考 accession、参考 FASTA 和分型解释提示。",
        "entries": build_entries(),
    }
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {OUT_PATH}")
    print(f"entries {len(payload['entries'])}")


if __name__ == "__main__":
    main()
