from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import openpyxl


DEFAULT_XLSX = Path("database/virus/bandavirus/CJ_group/CJ_group.xlsx")
DEFAULT_DB_CSV = Path("database/virus/bandavirus/bandavirus_db.csv")
DEFAULT_DB_FASTA = Path("database/virus/bandavirus/bandavirus_db.fasta")
DEFAULT_OUT_DIR = Path("database/virus/bandavirus/CJ_group/selected_refs")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract CJ_group bandavirus segment references from local bandavirus_db resources."
    )
    parser.add_argument("--xlsx", type=Path, default=DEFAULT_XLSX, help="Input CJ_group workbook")
    parser.add_argument("--db-csv", type=Path, default=DEFAULT_DB_CSV, help="bandavirus_db metadata CSV")
    parser.add_argument("--db-fasta", type=Path, default=DEFAULT_DB_FASTA, help="bandavirus_db FASTA")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR, help="Output directory")
    return parser.parse_args()


def normalize_accession(value: object) -> str:
    return str(value or "").strip().split(".", 1)[0]


def normalize_date(value: object) -> str:
    if isinstance(value, datetime):
        return value.strftime("%Y-%m-%d")
    if value is None:
        return ""
    text = str(value).strip()
    return text


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
            accession = normalize_accession(header.split()[0].split("|")[0])
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
            rows[normalize_accession(row.get("Accession", ""))] = row
    return rows


def load_cj_rows(path: Path) -> list[dict[str, str]]:
    wb = openpyxl.load_workbook(path, data_only=True)
    ws = wb.active
    rows: list[dict[str, str]] = []
    for raw_row in ws.iter_rows(min_row=4, values_only=True):
        if not any(raw_row):
            continue
        viral_strain, isolation_source, isolation_date, geographic_origin, l_acc, l_geno, m_acc, m_geno, s_acc, s_geno = raw_row[:10]
        strain = str(viral_strain or "").strip()
        common = {
            "viral_strain": strain,
            "isolation_source": str(isolation_source or "").strip(),
            "isolation_date": normalize_date(isolation_date),
            "geographic_origin": str(geographic_origin or "").strip(),
        }
        for segment, accession, genotype in (
            ("L", l_acc, l_geno),
            ("M", m_acc, m_geno),
            ("S", s_acc, s_geno),
        ):
            row = dict(common)
            row["segment"] = segment
            row["accession"] = normalize_accession(accession)
            row["genotype"] = str(genotype or "").strip()
            rows.append(row)
    return rows


def build_segment_header(row: dict[str, str], db_row: dict[str, str], sequence: str) -> str:
    return " | ".join(
        [
            f"strain={row['viral_strain']}",
            f"segment={row['segment']}",
            f"genotype={row['genotype']}",
            f"accession={row['accession']}",
            f"length={len(sequence)}",
            f"organism={db_row.get('Organism_Name', '').strip()}",
            f"isolate={db_row.get('Isolate', '').strip()}",
        ]
    )


def main() -> int:
    args = parse_args()
    xlsx_path = args.xlsx.expanduser().resolve()
    db_csv_path = args.db_csv.expanduser().resolve()
    db_fasta_path = args.db_fasta.expanduser().resolve()
    out_dir = args.out_dir.expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    strain_dir = out_dir / "strains"
    subtype_dir = out_dir / "subtypes"
    strain_dir.mkdir(parents=True, exist_ok=True)
    subtype_dir.mkdir(parents=True, exist_ok=True)

    cj_rows = load_cj_rows(xlsx_path)
    db_rows = load_db_rows(db_csv_path)
    fasta_map = read_fasta_map(db_fasta_path)

    missing = [row["accession"] for row in cj_rows if row["accession"] not in db_rows or row["accession"] not in fasta_map]
    if missing:
        raise SystemExit(f"Missing accessions in local db: {sorted(set(missing))}")

    manifest_rows: list[dict[str, str]] = []
    combined_entries: list[str] = []
    strain_entries: dict[str, list[str]] = defaultdict(list)
    subtype_entries: dict[tuple[str, str], list[str]] = defaultdict(list)
    subtype_manifest_counter: dict[tuple[str, str], list[str]] = defaultdict(list)

    for row in cj_rows:
        db_row = db_rows[row["accession"]]
        _, sequence = fasta_map[row["accession"]]
        header = build_segment_header(row, db_row, sequence)
        fasta_entry = f">{header}\n{wrap_sequence(sequence)}"
        combined_entries.append(fasta_entry)
        strain_entries[row["viral_strain"]].append(fasta_entry)
        subtype_entries[(row["segment"], row["genotype"])].append(fasta_entry)
        subtype_manifest_counter[(row["segment"], row["genotype"])].append(row["accession"])

        manifest_row = dict(db_row)
        manifest_row.update(
            {
                "CJ_Viral_Strain": row["viral_strain"],
                "CJ_Isolation_Source": row["isolation_source"],
                "CJ_Isolation_Date": row["isolation_date"],
                "CJ_Geographic_Origin": row["geographic_origin"],
                "CJ_Segment": row["segment"],
                "CJ_Genotype": row["genotype"],
                "Selected_From": "bandavirus_db.fasta",
                "Selected_Header": header,
                "Sequence_Length_From_FASTA": str(len(sequence)),
            }
        )
        manifest_rows.append(manifest_row)

    combined_fasta_path = out_dir / "CJ_group_segment_genomes.fasta"
    combined_fasta_path.write_text("\n".join(combined_entries).strip() + "\n", encoding="utf-8")

    manifest_columns = list(next(iter(db_rows.values())).keys()) + [
        "CJ_Viral_Strain",
        "CJ_Isolation_Source",
        "CJ_Isolation_Date",
        "CJ_Geographic_Origin",
        "CJ_Segment",
        "CJ_Genotype",
        "Selected_From",
        "Selected_Header",
        "Sequence_Length_From_FASTA",
    ]
    manifest_path = out_dir / "CJ_group_segment_genomes_manifest.tsv"
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=manifest_columns, delimiter="\t")
        writer.writeheader()
        for row in manifest_rows:
            writer.writerow({column: row.get(column, "") for column in manifest_columns})

    for strain, entries in strain_entries.items():
        safe_name = strain.replace("/", "_")
        (strain_dir / f"{safe_name}.fasta").write_text("\n".join(entries).strip() + "\n", encoding="utf-8")

    subtype_summary_rows: list[dict[str, str]] = []
    for (segment, genotype), entries in sorted(subtype_entries.items()):
        file_name = f"{segment}_{genotype}.fasta"
        (subtype_dir / file_name).write_text("\n".join(entries).strip() + "\n", encoding="utf-8")
        subtype_summary_rows.append(
            {
                "segment": segment,
                "genotype": genotype,
                "record_count": str(len(entries)),
                "accessions": ",".join(subtype_manifest_counter[(segment, genotype)]),
                "fasta_path": str((subtype_dir / file_name).resolve()),
            }
        )

    subtype_summary_path = out_dir / "CJ_group_subtype_manifest.tsv"
    with subtype_summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["segment", "genotype", "record_count", "accessions", "fasta_path"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(subtype_summary_rows)

    print(f"records\t{len(cj_rows)}")
    print(f"strains\t{len(strain_entries)}")
    print(f"subtype_groups\t{len(subtype_entries)}")
    print(f"combined_fasta\t{combined_fasta_path}")
    print(f"manifest\t{manifest_path}")
    print(f"subtype_manifest\t{subtype_summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
