from __future__ import annotations

import argparse
import csv
from pathlib import Path


def _read_fasta_records(path: Path) -> dict[str, tuple[str, str]]:
    records: dict[str, tuple[str, str]] = {}
    header = ""
    chunks: list[str] = []

    def _flush() -> None:
        nonlocal header, chunks
        if not header:
            return
        accession = header.split()[0].lstrip(">").strip()
        sequence = "".join(chunks).strip().upper()
        if accession and sequence:
            records[accession] = (header.lstrip(">"), sequence)
        header = ""
        chunks = []

    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            text = line.strip()
            if not text:
                continue
            if text.startswith(">"):
                _flush()
                header = text
                chunks = []
            elif header:
                chunks.append(text)
    _flush()
    return records


def main() -> int:
    parser = argparse.ArgumentParser(description="Deduplicate human rhinovirus complete genomes by species_group + type_label")
    parser.add_argument(
        "--in-dir",
        type=Path,
        default=Path("database/virus/rhinovirus/full_genomes"),
        help="Input directory containing combined fasta and manifest",
    )
    args = parser.parse_args()

    in_dir = args.in_dir.expanduser().resolve()
    manifest_path = in_dir / "human_rhinovirus_A_B_C_manifest.tsv"
    fasta_path = in_dir / "human_rhinovirus_A_B_C_complete_genomes.fasta"
    if not manifest_path.is_file() or not fasta_path.is_file():
        raise SystemExit("Missing input manifest or fasta.")

    fasta_records = _read_fasta_records(fasta_path)
    typed_best: dict[tuple[str, str], dict[str, str]] = {}
    untyped_rows: list[dict[str, str]] = []

    with manifest_path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            accession = str(row.get("accession") or "").strip()
            species_group = str(row.get("species_group") or "").strip()
            type_label = str(row.get("type_label") or "").strip()
            if not accession or accession not in fasta_records:
                continue
            if not type_label:
                untyped_rows.append(dict(row))
                continue
            key = (species_group, type_label)
            existing = typed_best.get(key)
            current_length = int(row.get("sequence_length") or 0)
            current_header = str(row.get("header") or "")
            if existing is None:
                typed_best[key] = dict(row)
                continue
            existing_length = int(existing.get("sequence_length") or 0)
            existing_header = str(existing.get("header") or "")
            if (current_length, current_header, accession) > (existing_length, existing_header, str(existing.get("accession") or "")):
                typed_best[key] = dict(row)

    typed_rows = sorted(
        typed_best.values(),
        key=lambda row: (
            str(row.get("species_group") or ""),
            int(str(row.get("type_label") or "999999")) if str(row.get("type_label") or "").isdigit() else 999999,
            str(row.get("accession") or ""),
        ),
    )
    untyped_rows = sorted(
        untyped_rows,
        key=lambda row: (
            str(row.get("species_group") or ""),
            str(row.get("accession") or ""),
        ),
    )

    typed_manifest_path = in_dir / "human_rhinovirus_A_B_C_dedup_by_type_manifest.tsv"
    typed_fasta_path = in_dir / "human_rhinovirus_A_B_C_dedup_by_type.fasta"
    untyped_manifest_path = in_dir / "human_rhinovirus_A_B_C_untyped_manifest.tsv"
    untyped_fasta_path = in_dir / "human_rhinovirus_A_B_C_untyped.fasta"

    columns = ["species_group", "species_name", "accession", "type_label", "sequence_length", "header", "title"]
    with typed_manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for row in typed_rows:
            writer.writerow({column: row.get(column, "") for column in columns})

    with typed_fasta_path.open("w", encoding="utf-8") as handle:
        for row in typed_rows:
            accession = str(row.get("accession") or "")
            header, sequence = fasta_records[accession]
            handle.write(f">{header}\n")
            for index in range(0, len(sequence), 80):
                handle.write(sequence[index:index + 80] + "\n")

    with untyped_manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for row in untyped_rows:
            writer.writerow({column: row.get(column, "") for column in columns})

    with untyped_fasta_path.open("w", encoding="utf-8") as handle:
        for row in untyped_rows:
            accession = str(row.get("accession") or "")
            if accession not in fasta_records:
                continue
            header, sequence = fasta_records[accession]
            handle.write(f">{header}\n")
            for index in range(0, len(sequence), 80):
                handle.write(sequence[index:index + 80] + "\n")

    print(f"[ok] typed representatives: {len(typed_rows)}")
    print(f"[ok] untyped sequences: {len(untyped_rows)}")
    print(f"[ok] wrote {typed_fasta_path}")
    print(f"[ok] wrote {typed_manifest_path}")
    print(f"[ok] wrote {untyped_fasta_path}")
    print(f"[ok] wrote {untyped_manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
