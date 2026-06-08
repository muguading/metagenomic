from __future__ import annotations

import csv
import json
from pathlib import Path

from bac_analysis_portal.store import PortalStore


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "public" / "zika-tutorial" / "data"
METADATA_PATH = DATA_DIR / "metadata.tsv"
SEQUENCES_PATH = DATA_DIR / "sequences.fasta"
OUTPUT_FASTA_DIR = DATA_DIR / "sample_library_fastas"


def parse_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    current_name = ""
    chunks: list[str] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_name:
                records[current_name] = "".join(chunks)
            current_name = line[1:].strip()
            chunks = []
            continue
        chunks.append(line)
    if current_name:
        records[current_name] = "".join(chunks)
    return records


def metadata_items(row: dict[str, str]) -> list[dict[str, object]]:
    location_value = {
        "province": "",
        "city": str(row.get("division") or "").strip(),
        "district": "",
        "detail": str(row.get("city") or "").strip(),
    }
    items = [
        {"key": "case_id", "label": "病例/事件编号", "type": "text", "options": [], "value": str(row.get("accession") or "").strip()},
        {"key": "surveillance_source", "label": "监测来源", "type": "select", "options": [], "value": "公开数据库"},
        {"key": "suspected_syndrome", "label": "疑似症候群", "type": "select", "options": [], "value": "其他"},
        {"key": "specimen_category", "label": "标本类别", "type": "select", "options": [], "value": "其他"},
        {"key": "collection_site", "label": "采样地点", "type": "location", "options": [], "value": location_value},
        {"key": "cluster_status", "label": "聚集性状态", "type": "select", "options": [], "value": "待判定"},
        {"key": "epidemiology_link", "label": "流行病学关联", "type": "text", "options": [], "value": str(row.get("region") or "").strip()},
        {"key": "traditional_result", "label": "传统检测结果", "type": "text", "options": [], "value": str(row.get("db") or "").strip()},
        {"key": "virus", "label": "病毒名称", "type": "text", "options": [], "value": str(row.get("virus") or "").strip()},
        {"key": "accession", "label": "公开库编号", "type": "text", "options": [], "value": str(row.get("accession") or "").strip()},
        {"key": "segment", "label": "片段", "type": "text", "options": [], "value": str(row.get("segment") or "").strip()},
        {"key": "authors", "label": "作者", "type": "text", "options": [], "value": str(row.get("authors") or "").strip()},
        {"key": "source_url", "label": "来源链接", "type": "text", "options": [], "value": str(row.get("url") or "").strip()},
        {"key": "paper_title", "label": "文献标题", "type": "text", "options": [], "value": str(row.get("title") or "").strip()},
        {"key": "journal", "label": "期刊", "type": "text", "options": [], "value": str(row.get("journal") or "").strip()},
        {"key": "paper_url", "label": "文献链接", "type": "text", "options": [], "value": str(row.get("paper_url") or "").strip()},
    ]
    return items


def main() -> None:
    if not METADATA_PATH.is_file():
        raise FileNotFoundError(f"Missing metadata file: {METADATA_PATH}")
    if not SEQUENCES_PATH.is_file():
        raise FileNotFoundError(f"Missing fasta file: {SEQUENCES_PATH}")

    OUTPUT_FASTA_DIR.mkdir(parents=True, exist_ok=True)
    store = PortalStore.from_project_root(PROJECT_ROOT)
    sequences = parse_fasta(SEQUENCES_PATH)

    with METADATA_PATH.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        inserted = 0
        missing_sequences: list[str] = []
        for row in reader:
            strain = str(row.get("strain") or "").strip()
            if not strain:
                continue
            sequence = sequences.get(strain)
            if not sequence:
                missing_sequences.append(strain)
                continue

            fasta_path = OUTPUT_FASTA_DIR / f"{strain.replace('/', '_')}.fasta"
            fasta_path.write_text(f">{strain}\n{sequence}\n", encoding="utf-8")

            collection_date = str(row.get("date") or "").strip().replace("XX", "01")
            location_json = json.dumps(
                {
                    "country": str(row.get("country") or "").strip(),
                    "province": "",
                    "city": str(row.get("division") or "").strip(),
                    "district": "",
                    "detail": str(row.get("city") or "").strip(),
                    "region": str(row.get("region") or "").strip(),
                },
                ensure_ascii=False,
            )
            metadata_json = json.dumps(metadata_items(row), ensure_ascii=False)

            store.upsert_sample_library_record(
                {
                    "sample_key": f"main::zika-tutorial::{strain}",
                    "genome_id": str(row.get("accession") or "").strip(),
                    "sample_name": strain,
                    "task_name": "zika_tutorial_dataset",
                    "owner": "admin",
                    "owner_group": "virus_tutorial",
                    "final_fasta_path": str(fasta_path.resolve()),
                    "species_name": "Zika virus",
                    "pathogen_type": "virus",
                    "genome_length": str(len(sequence)),
                    "description": f"Imported from public/zika-tutorial/data ({row.get('db') or 'public source'})",
                    "country": str(row.get("country") or "").strip(),
                    "location_json": location_json,
                    "sample_type": "public_reference",
                    "sequencing_method": "public_dataset",
                    "custom_metadata_json": metadata_json,
                    "sample_alias": str(row.get("accession") or "").strip(),
                    "sample_source": "public/zika-tutorial/data",
                    "collection_date": collection_date,
                    "host_info": str(row.get("authors") or "").strip(),
                    "note": f"{row.get('title') or ''} | {row.get('journal') or ''}".strip(" |"),
                    "library_scope": "main",
                    "visibility_scope": "group",
                }
            )
            inserted += 1

    print(f"Imported {inserted} Zika tutorial samples into sample_library.")
    if missing_sequences:
        print("Missing sequences for:")
        for item in missing_sequences:
            print(f"  - {item}")


if __name__ == "__main__":
    main()
