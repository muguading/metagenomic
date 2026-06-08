#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Dict, Iterable, Iterator, Tuple


HA_RE = re.compile(r"(H\d+)", re.IGNORECASE)
NA_RE = re.compile(r"(N\d+)", re.IGNORECASE)
CSV_NAME_RE = re.compile(r"\(([^()]+)\(([A-Za-z0-9,\-/]+)\)\)\s*$")
FASTA_HEADER_RE = re.compile(r"^(\S+)\s+(.+?)\s+(\d{4})(?:/\S*)?\s+\d+\s+\((HA|NA)\)$", re.IGNORECASE)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build influenza A HA/NA subtype reference FASTAs from flu.csv and Influnza_HANA.fa."
    )
    parser.add_argument(
        "--source-dir",
        default="virusdatabase/influenza",
        help="Directory containing flu.csv and Influnza_HANA.fa (default: %(default)s).",
    )
    parser.add_argument(
        "--database-root",
        default="database",
        help="Reference database root. Output will be written to <root>/virus/influenza_a/ (default: %(default)s).",
    )
    return parser.parse_args()


def iter_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    header = None
    seq_parts = []
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        yield header, "".join(seq_parts)


def normalize_accession(raw: str) -> str:
    return (raw or "").strip().split(".")[0].upper()


def extract_subtype(serotype: str, protein: str) -> str:
    text = (serotype or "").strip().upper()
    if not text:
        return ""
    if protein.upper() == "HA":
        match = HA_RE.search(text)
    else:
        match = NA_RE.search(text)
    return match.group(1).upper() if match else ""


def normalize_strain(raw: str) -> str:
    return re.sub(r"\s+", "", (raw or "").strip()).upper()


def extract_year(raw: str) -> str:
    match = re.search(r"(\d{4})", raw or "")
    return match.group(1) if match else ""


def parse_csv_name(raw_name: str) -> str:
    text = (raw_name or "").strip()
    match = CSV_NAME_RE.search(text)
    if match:
        return normalize_strain(match.group(1))
    return normalize_strain(text)


def load_annotations(csv_path: Path) -> Dict[Tuple[str, str, str], Dict[str, str]]:
    annotations: Dict[Tuple[str, str, str], Dict[str, str]] = {}
    with csv_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            protein = str(row.get("protein") or "").strip().upper()
            serotype = str(row.get("serotype") or "").strip()
            if protein not in {"HA", "NA"}:
                continue
            subtype = extract_subtype(serotype, protein)
            if not subtype:
                continue
            strain = parse_csv_name(str(row.get("name") or ""))
            year = extract_year(str(row.get("date") or "")) or extract_year(str(row.get("name") or ""))
            if not strain or not year:
                continue
            key = (strain, year, protein)
            existing = annotations.get(key)
            if existing and existing.get("subtype") != subtype:
                # Keep the first stable mapping when noisy metadata disagree.
                continue
            annotations[key] = {
                "protein": protein,
                "serotype": serotype,
                "subtype": subtype,
                "strain": strain,
                "year": year,
                "host": str(row.get("host") or "").strip(),
                "name": str(row.get("name") or "").strip(),
                "country": str(row.get("country") or "").strip(),
                "date": str(row.get("date") or "").strip(),
            }
    return annotations


def classify_record(
    header: str, sequence: str, annotations: Dict[Tuple[str, str, str], Dict[str, str]]
) -> Tuple[str, str] | None:
    if not sequence:
        return None
    match = FASTA_HEADER_RE.match(header)
    if not match:
        return None
    accession_raw, name, year, protein_hint = match.groups()
    accession = normalize_accession(accession_raw)
    strain = normalize_strain(name)
    protein = protein_hint.upper()
    ann = annotations.get((strain, year, protein))
    if not ann:
        return None
    subtype = ann["subtype"]
    clean_name = " ".join(name.split()) if name else accession
    title = f"{accession}|{subtype}|{clean_name}"
    return protein, f">{title}\n{sequence}\n"


def write_records(path: Path, records: Iterable[str]) -> int:
    count = 0
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for record in records:
            handle.write(record)
            count += 1
    return count


def main() -> int:
    args = parse_args()
    source_dir = Path(args.source_dir).expanduser().resolve()
    csv_path = source_dir / "flu.csv"
    fasta_path = source_dir / "Influnza_HANA.fa"
    if not csv_path.is_file():
        raise SystemExit(f"Missing source CSV: {csv_path}")
    if not fasta_path.is_file():
        raise SystemExit(f"Missing source FASTA: {fasta_path}")

    annotations = load_annotations(csv_path)
    ha_records = []
    na_records = []
    for header, sequence in iter_fasta(fasta_path):
        classified = classify_record(header, sequence, annotations)
        if not classified:
            continue
        protein, formatted = classified
        if protein == "HA":
            ha_records.append(formatted)
        elif protein == "NA":
            na_records.append(formatted)

    database_root = Path(args.database_root).expanduser().resolve()
    out_dir = database_root / "virus" / "influenza_a"
    ha_path = out_dir / "ha_subtypes.fa"
    na_path = out_dir / "na_subtypes.fa"
    ha_count = write_records(ha_path, ha_records)
    na_count = write_records(na_path, na_records)

    print(f"HA refs: {ha_count} -> {ha_path}")
    print(f"NA refs: {na_count} -> {na_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
