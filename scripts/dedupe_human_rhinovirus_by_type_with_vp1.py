#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path("/Users/wuhhh/Desktop/徐老师/代码/metagenomic")
FULL_GENOMES_DIR = ROOT / "database/virus/rhinovirus/full_genomes"
INPUT_FASTA = FULL_GENOMES_DIR / "human_rhinovirus_A_B_C_complete_genomes.fasta"
INPUT_MANIFEST = FULL_GENOMES_DIR / "human_rhinovirus_A_B_C_manifest.tsv"
INPUT_GFF_DIR = FULL_GENOMES_DIR / "gff3_all"
WITH_VP1_FASTA = FULL_GENOMES_DIR / "human_rhinovirus_A_B_C_with_vp1.fasta"
WITH_VP1_MANIFEST = FULL_GENOMES_DIR / "human_rhinovirus_A_B_C_with_vp1_manifest.tsv"
DEDUP_FASTA = FULL_GENOMES_DIR / "human_rhinovirus_A_B_C_dedup_by_type_with_vp1.fasta"
DEDUP_MANIFEST = FULL_GENOMES_DIR / "human_rhinovirus_A_B_C_dedup_by_type_with_vp1_manifest.tsv"
NO_VP1_TYPED_MANIFEST = FULL_GENOMES_DIR / "human_rhinovirus_A_B_C_typed_without_vp1_manifest.tsv"


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


def parse_attributes(raw: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in raw.split(";"):
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        attrs[key] = value
    return attrs


def gff_has_vp1(gff_path: Path) -> bool:
    if not gff_path.exists():
        return False
    with gff_path.open() as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                continue
            feature_type = fields[2]
            attrs = parse_attributes(fields[8])
            combined = " ".join(
                [
                    attrs.get("product", ""),
                    attrs.get("gene", ""),
                    attrs.get("Name", ""),
                    attrs.get("Note", ""),
                ]
            ).upper()
            if feature_type in {"gene", "CDS", "mature_protein_region_of_CDS"} and ("VP1" in combined or "1D" in combined):
                return True
    return False


def write_fasta(path: Path, rows: list[dict[str, str]], records: dict[str, str]) -> None:
    with path.open("w") as handle:
        for row in rows:
            header = row["header"]
            sequence = records[header]
            handle.write(f">{header}\n")
            for idx in range(0, len(sequence), 80):
                handle.write(sequence[idx : idx + 80] + "\n")


def write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    records = parse_fasta(INPUT_FASTA)
    with INPUT_MANIFEST.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    with_vp1_rows: list[dict[str, str]] = []
    typed_without_vp1_rows: list[dict[str, str]] = []

    for row in rows:
        accession = str(row.get("accession") or "").strip()
        accession_root = accession.split(".", 1)[0]
        gff_path = INPUT_GFF_DIR / f"{accession_root}.gff3"
        type_label = str(row.get("type_label") or "").strip()
        enriched = dict(row)
        enriched["accession_root"] = accession_root
        enriched["gff_path"] = str(gff_path)
        enriched["vp1_in_gff"] = "yes" if gff_has_vp1(gff_path) else "no"
        if enriched["vp1_in_gff"] == "yes":
            with_vp1_rows.append(enriched)
        elif type_label:
            typed_without_vp1_rows.append(enriched)

    write_fasta(WITH_VP1_FASTA, with_vp1_rows, records)
    write_manifest(WITH_VP1_MANIFEST, with_vp1_rows)
    write_manifest(NO_VP1_TYPED_MANIFEST, typed_without_vp1_rows)

    dedup_map: dict[tuple[str, str], dict[str, str]] = {}
    for row in with_vp1_rows:
        key = (row["species_group"], row["type_label"])
        current = dedup_map.get(key)
        if current is None or int(row["sequence_length"]) > int(current["sequence_length"]):
            dedup_map[key] = row

    dedup_rows = sorted(
        dedup_map.values(),
        key=lambda row: (
            row["species_group"],
            0 if str(row["type_label"]).isdigit() else 1,
            int(row["type_label"]) if str(row["type_label"]).isdigit() else str(row["type_label"]),
        ),
    )
    write_fasta(DEDUP_FASTA, dedup_rows, records)
    write_manifest(DEDUP_MANIFEST, dedup_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
