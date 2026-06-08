from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


ROOT = Path(__file__).resolve().parent.parent
REF_DIR = ROOT / "database" / "virus" / "norovirus" / "cdc_typing_refs"
BACKUP_DIR = REF_DIR / "backup_refs"


def sanitize_tree_label(value: str) -> str:
    sanitized = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value or "").strip())
    return sanitized.strip("_") or "sample"


def build_backup_refs(gene_key: str) -> dict[str, int]:
    source_tsv = REF_DIR / f"cdc_norovirus_{gene_key}_refs.tsv"
    source_fasta = REF_DIR / f"cdc_norovirus_{gene_key}_refs.fasta"
    output_tsv = BACKUP_DIR / f"cdc_norovirus_{gene_key}_backup_refs.tsv"
    output_fasta = BACKUP_DIR / f"cdc_norovirus_{gene_key}_backup_refs.fasta"

    with source_tsv.open("r", encoding="utf-8-sig", errors="ignore", newline="") as handle:
        rows = [dict(row) for row in csv.DictReader(handle, delimiter="\t") if row]
    records_by_accession: dict[str, SeqRecord] = {}
    for record in SeqIO.parse(str(source_fasta), "fasta"):
        accession = str(record.id).split("_", 1)[0].strip()
        if accession and accession not in records_by_accession:
            records_by_accession[accession] = record

    subtype_counts: dict[str, int] = defaultdict(int)
    selected_rows: list[dict[str, str]] = []
    selected_records: list[SeqRecord] = []
    for row in rows:
        subtype = str(row.get("subtype") or "").strip()
        accession = str(row.get("accession") or "").strip()
        if not subtype or not accession:
            continue
        if subtype_counts[subtype] >= 3:
            continue
        template = records_by_accession.get(accession)
        if template is None:
            continue
        subtype_counts[subtype] += 1
        rank = subtype_counts[subtype]
        selected_rows.append(
            {
                **row,
                "backup_rank": str(rank),
                "selection_rule": "per_subtype_top3_in_source_order",
            }
        )
        canonical_id = sanitize_tree_label(f"{accession}_{subtype}_{str(row.get('gene') or gene_key).strip()}")
        selected_records.append(
            SeqRecord(
                template.seq,
                id=canonical_id,
                name=canonical_id,
                description=f"{canonical_id} {template.description}".strip(),
            )
        )

    BACKUP_DIR.mkdir(parents=True, exist_ok=True)
    with output_tsv.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = ["source", "gene", "subtype", "accession", "label", "source_url", "backup_rank", "selection_rule"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(selected_rows)
    SeqIO.write(selected_records, str(output_fasta), "fasta")
    return {
        "subtype_count": len(subtype_counts),
        "record_count": len(selected_rows),
    }


def main() -> None:
    summary_rows: list[list[str]] = []
    for gene_key in ("rdrp", "vp1"):
        stats = build_backup_refs(gene_key)
        summary_rows.append([gene_key.upper(), str(stats["subtype_count"]), str(stats["record_count"])])
    summary_path = BACKUP_DIR / "cdc_norovirus_backup_ref_counts.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["gene", "subtype_count", "record_count"])
        writer.writerows(summary_rows)


if __name__ == "__main__":
    main()
