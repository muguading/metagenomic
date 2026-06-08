#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path("/Users/wuhhh/Desktop/徐老师/代码/metagenomic")
FULL_GENOMES_DIR = ROOT / "database/virus/rhinovirus/full_genomes"
REFERENCE_DIR = ROOT / "database/virus/rhinovirus/reference_genomes"
INPUT_FASTA = FULL_GENOMES_DIR / "human_rhinovirus_A_B_C_dedup_by_type.fasta"
INPUT_MANIFEST = FULL_GENOMES_DIR / "human_rhinovirus_A_B_C_dedup_by_type_manifest.tsv"
OUTPUT_FASTA = REFERENCE_DIR / "human_rhinovirus_representative_genomes.fasta"
OUTPUT_MANIFEST = REFERENCE_DIR / "human_rhinovirus_representative_genomes.tsv"
GFF_DIR = REFERENCE_DIR / "gff3"


def parse_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    header: str | None = None
    chunks: list[str] = []
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records[header] = "".join(chunks)
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        records[header] = "".join(chunks)
    return records


def normalized_type_name(species_group: str, type_label: str) -> str:
    return f"RV-{species_group}{type_label}"


def normalized_record_id(accession_root: str, species_group: str, type_label: str) -> str:
    return f"{accession_root}_{normalized_type_name(species_group, type_label)}"


def main() -> None:
    REFERENCE_DIR.mkdir(parents=True, exist_ok=True)
    records = parse_fasta(INPUT_FASTA)

    with INPUT_MANIFEST.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)

    manifest_rows: list[dict[str, str]] = []

    with OUTPUT_FASTA.open("w") as fasta_out:
        for row in rows:
            accession = row["accession"]
            accession_root = accession.split(".", 1)[0]
            header = row["header"]
            sequence = records[header]
            species_group = row["species_group"]
            type_label = row["type_label"]
            species_name = row["species_name"]
            normalized_type = normalized_type_name(species_group, type_label)
            record_id = normalized_record_id(accession_root, species_group, type_label)
            gff_path = GFF_DIR / f"{accession_root}.gff3"
            gff_status = "downloaded" if gff_path.exists() else "missing"

            fasta_out.write(f">{record_id} {species_name} type {normalized_type} accession {accession}\n")
            for start in range(0, len(sequence), 80):
                fasta_out.write(sequence[start : start + 80] + "\n")

            manifest_rows.append(
                {
                    "record_id": record_id,
                    "species_group": species_group,
                    "species_name": species_name,
                    "type_label": type_label,
                    "normalized_type": normalized_type,
                    "accession": accession,
                    "accession_root": accession_root,
                    "sequence_length": str(len(sequence)),
                    "source_header": header,
                    "source_title": row["title"],
                    "fasta_path": str(OUTPUT_FASTA),
                    "gff_path": str(gff_path),
                    "gff_status": gff_status,
                }
            )

    fieldnames = [
        "record_id",
        "species_group",
        "species_name",
        "type_label",
        "normalized_type",
        "accession",
        "accession_root",
        "sequence_length",
        "source_header",
        "source_title",
        "fasta_path",
        "gff_path",
        "gff_status",
    ]
    with OUTPUT_MANIFEST.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)


if __name__ == "__main__":
    main()
