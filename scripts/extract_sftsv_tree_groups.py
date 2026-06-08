from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

import openpyxl


DEFAULT_XLSX = Path("database/virus/bandavirus/A_Fgroup/SFTSV_tree_sample_groups.xlsx")
DEFAULT_DB_CSV = Path("database/virus/bandavirus/bandavirus_db.csv")
DEFAULT_DB_FASTA = Path("database/virus/bandavirus/bandavirus_db.fasta")
DEFAULT_OUT_DIR = Path("database/virus/bandavirus/A_Fgroup/selected_refs")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract SFTSV tree-grouped references from local bandavirus_db resources."
    )
    parser.add_argument("--xlsx", type=Path, default=DEFAULT_XLSX, help="Input grouping workbook")
    parser.add_argument("--db-csv", type=Path, default=DEFAULT_DB_CSV, help="bandavirus_db metadata CSV")
    parser.add_argument("--db-fasta", type=Path, default=DEFAULT_DB_FASTA, help="bandavirus_db FASTA")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR, help="Output directory")
    return parser.parse_args()


def normalize_accession(value: object) -> str:
    return str(value or "").strip().split("/", 1)[0].split(".", 1)[0]


def safe_name(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in {"_", "-", "."} else "_" for ch in value).strip("_") or "unknown"


def wrap_sequence(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[index:index + width] for index in range(0, len(sequence), width))


def read_fasta_map(path: Path) -> dict[str, tuple[str, str]]:
    fasta_map: dict[str, tuple[str, str]] = {}
    header = ""
    accession = ""
    seq_lines: list[str] = []
    for raw_line in path.open(encoding="utf-8", errors="ignore"):
        line = raw_line.rstrip("\n")
        if line.startswith(">"):
            if header:
                fasta_map[accession] = (header, "".join(seq_lines).upper())
            header = line[1:].strip()
            accession = header.split()[0].split("|")[0].split(".", 1)[0].strip()
            seq_lines = []
        else:
            seq_lines.append(line.strip())
    if header:
        fasta_map[accession] = (header, "".join(seq_lines).upper())
    return fasta_map


def load_db_rows(path: Path) -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            accession = (row.get("Accession") or "").strip().split(".", 1)[0]
            rows[accession] = row
    return rows


def load_group_rows(path: Path) -> list[dict[str, str]]:
    wb = openpyxl.load_workbook(path, data_only=True)
    ws = wb["All_segments"]
    rows: list[dict[str, str]] = []
    for raw_row in ws.iter_rows(min_row=2, values_only=True):
        if not any(raw_row):
            continue
        segment, sample_id, group = raw_row[:3]
        rows.append(
            {
                "segment": str(segment or "").strip(),
                "sample_id": str(sample_id or "").strip(),
                "group": str(group or "").strip(),
                "accession": normalize_accession(sample_id),
            }
        )
    return rows


def normalize_key(text: str) -> str:
    return "".join(ch for ch in text.lower() if ch.isalnum())


def resolve_accession(
    row: dict[str, str],
    db_rows: dict[str, dict[str, str]],
) -> tuple[str | None, str]:
    accession = row["accession"]
    if accession in db_rows:
        return accession, "direct"

    manual_map = {
        ("L", "NC018136"): "NC_018136",
        ("M", "PQ761121"): "PP761121",
        ("M", "PQ761122"): "PP761122",
        ("M", "KC473534"): "KC473541",
        ("M", "KF356640"): "KF356540",
        ("M", "OQ670930"): "JQ670930",
        ("S", "AB898541"): "AB985541",
        ("L", "SFTSZJK1"): "PX226566",
        ("M", "SFTSZJK1"): "PX226564",
        ("S", "SFTSZJK1"): "PX226562",
        ("L", "SFTSZJK2"): "PX226567",
        ("M", "SFTSZJK2"): "PX226565",
        ("S", "SFTSZJK2"): "PX226563",
    }
    mapped = manual_map.get((row["segment"], accession))
    if mapped and mapped in db_rows:
        return mapped, "manual_map"

    sample_id = row["sample_id"]
    parts = sample_id.split("/")
    label = parts[1] if len(parts) > 1 else parts[0]
    label_norm = normalize_key(label)
    candidates: list[str] = []
    for candidate_acc, db_row in db_rows.items():
        if (db_row.get("Segment") or "").strip() != row["segment"]:
            continue
        haystacks = [
            db_row.get("Accession", ""),
            db_row.get("Isolate", ""),
            db_row.get("GenBank_Title", ""),
        ]
        hay_norm = normalize_key(" ".join(haystacks))
        if label_norm and label_norm in hay_norm:
            candidates.append(candidate_acc)
    if len(candidates) == 1:
        return candidates[0], "label_match"

    return None, "unresolved"


def build_header(group_row: dict[str, str], db_row: dict[str, str], sequence: str) -> str:
    isolate = (db_row.get("Isolate") or "").strip()
    return " | ".join(
        [
            f"segment={group_row['segment']}",
            f"group={group_row['group']}",
            f"accession={group_row.get('resolved_accession', group_row['accession'])}",
            f"sample_id={group_row['sample_id']}",
            f"length={len(sequence)}",
            f"organism={db_row.get('Organism_Name', '').strip()}",
            f"isolate={isolate}",
        ]
    )


def main() -> int:
    args = parse_args()
    xlsx_path = args.xlsx.expanduser().resolve()
    db_csv_path = args.db_csv.expanduser().resolve()
    db_fasta_path = args.db_fasta.expanduser().resolve()
    out_dir = args.out_dir.expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    segment_group_dir = out_dir / "segment_groups"
    segment_group_dir.mkdir(parents=True, exist_ok=True)

    group_rows = load_group_rows(xlsx_path)
    db_rows = load_db_rows(db_csv_path)
    fasta_map = read_fasta_map(db_fasta_path)

    manifest_rows: list[dict[str, str]] = []
    combined_entries: list[str] = []
    grouped_entries: dict[tuple[str, str], list[str]] = defaultdict(list)
    grouped_sample_ids: dict[tuple[str, str], list[str]] = defaultdict(list)
    missing_rows: list[dict[str, str]] = []

    for row in group_rows:
        resolved_accession, resolve_mode = resolve_accession(row, db_rows)
        if not resolved_accession or resolved_accession not in fasta_map:
            missing_rows.append(
                {
                    "segment": row["segment"],
                    "group": row["group"],
                    "sample_id": row["sample_id"],
                    "requested_accession": row["accession"],
                    "resolve_mode": resolve_mode,
                }
            )
            continue
        db_row = db_rows[resolved_accession]
        _, sequence = fasta_map[resolved_accession]
        row = dict(row)
        row["resolved_accession"] = resolved_accession
        header = build_header(row, db_row, sequence)
        fasta_entry = f">{header}\n{wrap_sequence(sequence)}"
        combined_entries.append(fasta_entry)
        grouped_entries[(row["segment"], row["group"])].append(fasta_entry)
        grouped_sample_ids[(row["segment"], row["group"])].append(row["sample_id"])

        manifest_row = dict(db_row)
        manifest_row.update(
            {
                "Tree_Segment": row["segment"],
                "Tree_Group": row["group"],
                "Tree_Sample_ID": row["sample_id"],
                "Tree_Requested_Accession": row["accession"],
                "Tree_Resolved_Accession": resolved_accession,
                "Resolve_Mode": resolve_mode,
                "Selected_From": "bandavirus_db.fasta",
                "Selected_Header": header,
                "Sequence_Length_From_FASTA": str(len(sequence)),
            }
        )
        manifest_rows.append(manifest_row)

    if not combined_entries:
        raise SystemExit("No records could be resolved from local bandavirus_db resources.")

    combined_fasta_path = out_dir / "SFTSV_tree_group_segments.fasta"
    combined_fasta_path.write_text("\n".join(combined_entries).strip() + "\n", encoding="utf-8")

    manifest_columns = list(next(iter(db_rows.values())).keys()) + [
        "Tree_Segment",
        "Tree_Group",
        "Tree_Sample_ID",
        "Tree_Requested_Accession",
        "Tree_Resolved_Accession",
        "Resolve_Mode",
        "Selected_From",
        "Selected_Header",
        "Sequence_Length_From_FASTA",
    ]
    manifest_path = out_dir / "SFTSV_tree_group_segments_manifest.tsv"
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=manifest_columns, delimiter="\t")
        writer.writeheader()
        for row in manifest_rows:
            writer.writerow({column: row.get(column, "") for column in manifest_columns})

    summary_rows: list[dict[str, str]] = []
    for key in sorted(grouped_entries):
        segment, group = key
        file_name = f"{safe_name(segment)}_{safe_name(group)}.fasta"
        fasta_path = segment_group_dir / file_name
        fasta_path.write_text("\n".join(grouped_entries[key]).strip() + "\n", encoding="utf-8")
        summary_rows.append(
            {
                "segment": segment,
                "group": group,
                "record_count": str(len(grouped_entries[key])),
                "sample_ids": " ; ".join(grouped_sample_ids[key]),
                "fasta_path": str(fasta_path.resolve()),
            }
        )

    summary_path = out_dir / "SFTSV_tree_group_summary.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["segment", "group", "record_count", "sample_ids", "fasta_path"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    missing_path = out_dir / "SFTSV_tree_group_missing.tsv"
    with missing_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["segment", "group", "sample_id", "requested_accession", "resolve_mode"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(missing_rows)

    counts = Counter((row["segment"], row["group"]) for row in group_rows)
    print(f"records\t{len(group_rows)}")
    print(f"resolved_records\t{len(manifest_rows)}")
    print(f"missing_records\t{len(missing_rows)}")
    print(f"groups\t{len(counts)}")
    print(f"combined_fasta\t{combined_fasta_path}")
    print(f"manifest\t{manifest_path}")
    print(f"summary\t{summary_path}")
    print(f"missing\t{missing_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
