from __future__ import annotations

import csv
import json
import re
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
HEPATITIS_DB = ROOT / "database" / "virus" / "Hepatovirus"
OUT_DIR = ROOT / "database" / "knowledge_base" / "typing"
OUT_PATH = OUT_DIR / "hepatitis_virus_typing.json"

BROAD_PROFILES = {
    "HAV": {
        "typing_source": "A",
        "species": "Hepatitis A virus",
        "common_name": "甲型肝炎病毒",
        "pathogen_id": "hepatitis_a_virus",
        "panel": "肝炎病毒 HAV 大亚型及子亚型知识库",
        "subtype_label": "HAV 子亚型",
        "interpretation": "HAV 命中后应进一步结合 IA/IB/IIA/IIB/IIIA/IIIB 等子亚型进行分子流调解释。",
    },
    "HBV": {
        "typing_source": "B",
        "species": "Hepatitis B virus",
        "common_name": "乙型肝炎病毒",
        "pathogen_id": "hepatitis_b_virus",
        "panel": "肝炎病毒 HBV 基因型知识库",
        "subtype_label": "HBV genotype",
        "interpretation": "HBV genotype 可作为慢性乙肝分子流行病学和参考株选择的补充分型背景。",
    },
    "HCV": {
        "typing_source": "C",
        "species": "Hepatitis C virus",
        "common_name": "丙型肝炎病毒",
        "pathogen_id": "hepatitis_c_virus",
        "panel": "肝炎病毒 HCV genotype/subtype 知识库",
        "subtype_label": "HCV genotype/subtype",
        "interpretation": "HCV genotype/subtype 是丙肝分子流调和治疗背景解释的重要分层标签。",
    },
    "HDV": {
        "typing_source": "D",
        "species": "Hepatitis D virus",
        "common_name": "丁型肝炎病毒",
        "pathogen_id": "hepatitis_d_virus",
        "panel": "肝炎病毒 HDV genotype 知识库",
        "subtype_label": "HDV genotype",
        "interpretation": "HDV genotype 应结合 HBV 背景解释，用于丁肝谱系和区域传播背景补充。",
    },
    "HEV": {
        "typing_source": "E",
        "species": "Hepatitis E virus",
        "common_name": "戊型肝炎病毒",
        "pathogen_id": "hepatitis_e_virus",
        "panel": "肝炎病毒 HEV genotype/subtype 知识库",
        "subtype_label": "HEV genotype/subtype",
        "interpretation": "HEV genotype/subtype 可提示人兽共患、食源性暴露和区域传播背景；结果需结合样本类型与覆盖度解释。",
    },
}

MANIFESTS = {
    "HAV": HEPATITIS_DB / "HAV_subtypes" / "hav_subtype_complete_genomes_manifest.tsv",
    "HBV": HEPATITIS_DB / "typingB_reference_genomes" / "typingB_reference_genomes_manifest.tsv",
    "HCV": HEPATITIS_DB / "typingC_reference_genomes" / "typingC_reference_genomes_manifest.tsv",
    "HDV": HEPATITIS_DB / "typingD_reference_genomes" / "typingD_reference_genomes_manifest.tsv",
    "HEV": HEPATITIS_DB / "typingE_reference_genomes" / "typingE_reference_genomes_manifest.tsv",
}


def read_manifest(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def infer_subtype(row: dict[str, str], broad_type: str) -> str:
    genotype = str(row.get("genotype") or "").strip()
    if genotype:
        return genotype
    text = " ".join(str(row.get(key) or "") for key in ("virus_name", "header", "abbrev", "species"))
    match = re.search(r"\bgenotype\s*([0-9]+[A-Za-z]?|[A-Za-z])\b", text, flags=re.IGNORECASE)
    if match:
        return match.group(1)
    match = re.search(r"\bsubtype\s*([0-9]+[A-Za-z]?|[A-Za-z])\b", text, flags=re.IGNORECASE)
    if match:
        return match.group(1)
    abbrev = str(row.get("abbrev") or "").strip()
    upper = abbrev.upper()
    if broad_type == "HBV" and upper.startswith("HBV-"):
        return abbrev.split("-", 1)[1].strip()
    if broad_type in {"HCV", "HDV"} and upper.startswith(broad_type):
        return abbrev[len(broad_type):].strip("-_ ")
    return abbrev or broad_type


def clean_reference(row: dict[str, str]) -> dict[str, str]:
    return {
        "accession": str(row.get("accession") or "").strip(),
        "isolate": str(row.get("isolate") or "").strip(),
        "header": str(row.get("header") or "").strip(),
        "virus_name": str(row.get("virus_name") or "").strip(),
        "species": str(row.get("species") or "").strip(),
        "sequence_length": str(row.get("sequence_length") or "").strip(),
        "fasta_path": str(row.get("fasta_path") or "").strip(),
        "gff3_path": str(row.get("gff3_path") or "").strip(),
    }


def aliases_for(broad_type: str, subtype: str) -> list[str]:
    values = [
        subtype,
        f"{broad_type} {subtype}",
        f"{broad_type}-{subtype}",
        f"{broad_type}{subtype}",
        f"genotype {subtype}",
        f"subtype {subtype}",
    ]
    seen: list[str] = []
    for value in values:
        if value and value not in seen:
            seen.append(value)
    return seen


def build_entries() -> list[dict[str, object]]:
    entries: list[dict[str, object]] = []
    for broad_type, profile in BROAD_PROFILES.items():
        rows = [row for row in read_manifest(MANIFESTS[broad_type]) if str(row.get("status") or "").strip().lower() in {"", "ok"}]
        by_subtype: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row in rows:
            subtype = infer_subtype(row, broad_type).strip()
            if not subtype:
                continue
            by_subtype[subtype].append(row)
        entries.append(
            {
                "id": f"hepatitis_{broad_type.lower()}",
                "level": "broad_type",
                "broad_type": broad_type,
                "typing_source": profile["typing_source"],
                "species": profile["species"],
                "common_name": profile["common_name"],
                "pathogen_id": profile["pathogen_id"],
                "panel": f"肝炎病毒 {broad_type} 大亚型知识库",
                "serotype": broad_type,
                "aliases": [broad_type, profile["species"], profile["common_name"]],
                "subtype_panel": profile["panel"],
                "subtypes": sorted(by_subtype.keys(), key=lambda item: (len(item), item)),
                "reference_count": len(rows),
                "reference_accessions": sorted({str(row.get("accession") or "").strip() for row in rows if str(row.get("accession") or "").strip()}),
                "interpretation": profile["interpretation"],
            }
        )
        for subtype in sorted(by_subtype.keys(), key=lambda item: (len(item), item)):
            subtype_rows = by_subtype[subtype]
            accessions = sorted({str(row.get("accession") or "").strip() for row in subtype_rows if str(row.get("accession") or "").strip()})
            entries.append(
                {
                    "id": f"hepatitis_{broad_type.lower()}_{re.sub(r'[^A-Za-z0-9]+', '_', subtype).strip('_').lower()}",
                    "level": "subtype",
                    "broad_type": broad_type,
                    "typing_source": profile["typing_source"],
                    "species": profile["species"],
                    "common_name": profile["common_name"],
                    "pathogen_id": profile["pathogen_id"],
                    "panel": profile["panel"],
                    "serotype": subtype,
                    "subtype": subtype,
                    "subtype_type": profile["subtype_label"],
                    "aliases": aliases_for(broad_type, subtype),
                    "reference_count": len(subtype_rows),
                    "reference_accessions": accessions,
                    "references": [clean_reference(row) for row in subtype_rows],
                    "interpretation": f"{broad_type} {subtype} 属于{profile['subtype_label']}；可作为 {profile['common_name']} 的分子分型、参考株选择和流行病学关联解释依据。",
                }
            )
    return entries


def main() -> None:
    payload = {
        "schema_version": "v1",
        "id": "hepatitis_virus_typing",
        "title": "肝炎病毒大亚型与子亚型知识库",
        "description": "关联 HAV/HBV/HCV/HDV/HEV 大亚型、子亚型/基因型、参考 accession、参考 FASTA/GFF3 与报告判读提示。",
        "entries": build_entries(),
    }
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {OUT_PATH}")
    print(f"entries {len(payload['entries'])}")


if __name__ == "__main__":
    main()
