#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DB_DIR = ROOT / "database" / "virus" / "Hepatovirus"
OUT_DIR = DB_DIR / "broad_reference_genomes"


def load_manifest(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t") if row]


def broad_type_for_source(source: str, row: dict[str, str]) -> str:
    source = str(source or "").strip().upper()
    if source == "A":
        return str(row.get("abbrev") or "").strip()
    if source == "B":
        return "HBV"
    if source == "C":
        return "HCV"
    if source == "D":
        return "HDV"
    if source == "E":
        return "HEV"
    return str(row.get("abbrev") or "").strip()


def species_label_for_source(source: str, row: dict[str, str]) -> str:
    source = str(source or "").strip().upper()
    if source == "A":
        broad = str(row.get("abbrev") or "").strip().upper()
        if broad == "HAV":
            return "Hepatitis A virus"
        return str(row.get("species") or row.get("virus_name") or "Hepatovirus").strip()
    if source == "B":
        return "Hepatitis B virus"
    if source == "C":
        return "Hepatitis C virus"
    if source == "D":
        return "Hepatitis D virus"
    if source == "E":
        return "Hepatitis E virus"
    return str(row.get("species") or row.get("virus_name") or "").strip()


def source_manifest(source: str) -> Path:
    source = str(source or "").strip().upper()
    if source == "A":
        return (DB_DIR / "reference_genomes" / "hepatovirus_typing_reference_genomes_manifest.tsv").resolve()
    return (DB_DIR / f"typing{source}_reference_genomes" / f"typing{source}_reference_genomes_manifest.tsv").resolve()


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    combined_fasta = OUT_DIR / "hepatitis_broad_reference_genomes.fasta"
    combined_manifest = OUT_DIR / "hepatitis_broad_reference_genomes_manifest.tsv"

    merged_rows: list[dict[str, str]] = []
    fasta_paths: list[Path] = []
    for source in ["A", "B", "C", "D", "E"]:
        manifest_path = source_manifest(source)
        if not manifest_path.is_file():
            raise SystemExit(f"missing manifest: {manifest_path}")
        for row in load_manifest(manifest_path):
            fasta_path = Path(str(row.get("fasta_path") or "").strip())
            if not fasta_path.is_file():
                continue
            merged = dict(row)
            merged["typing_source"] = source
            merged["broad_type"] = broad_type_for_source(source, row)
            merged["species_label"] = species_label_for_source(source, row)
            merged_rows.append(merged)
            fasta_paths.append(fasta_path)

    merged_rows.sort(key=lambda item: (str(item.get("typing_source") or ""), str(item.get("accession") or "")))
    with combined_fasta.open("w", encoding="utf-8") as out_handle:
        for row in merged_rows:
            fasta_path = Path(str(row.get("fasta_path") or "").strip())
            out_handle.write(fasta_path.read_text(encoding="utf-8", errors="ignore").rstrip() + "\n")

    fieldnames = [
        "typing_source",
        "broad_type",
        "species_label",
        "genus",
        "species",
        "virus_name",
        "isolate",
        "accession",
        "available_sequence",
        "abbrev",
        "header",
        "sequence_length",
        "fasta_path",
        "gff3_path",
        "fasta_status",
        "gff3_status",
        "status",
    ]
    with combined_manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in merged_rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})

    print(f"[done] merged {len(merged_rows)} references -> {OUT_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
