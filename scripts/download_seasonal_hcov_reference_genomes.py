#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import sys
import urllib.parse
import urllib.request
from pathlib import Path


EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

REFERENCE_GENOMES = [
    {
        "virus": "HCoV-229E",
        "species": "Human coronavirus 229E",
        "accession": "NC_002645.1",
    },
    {
        "virus": "HCoV-NL63",
        "species": "Human coronavirus NL63",
        "accession": "NC_005831.2",
    },
    {
        "virus": "HCoV-OC43",
        "species": "Human coronavirus OC43",
        "accession": "NC_006213.1",
    },
    {
        "virus": "HCoV-HKU1",
        "species": "Human coronavirus HKU1",
        "accession": "NC_006577.2",
    },
]


def fetch_fasta(accession: str) -> str:
    params = {
        "db": "nuccore",
        "id": accession,
        "rettype": "fasta",
        "retmode": "text",
    }
    url = f"{EUTILS_BASE}?{urllib.parse.urlencode(params)}"
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "metagenomic-seasonal-hcov-downloader/1.0"},
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        return response.read().decode("utf-8", errors="ignore")


def parse_fasta(text: str) -> tuple[str, str]:
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines or not lines[0].startswith(">"):
        raise ValueError("Invalid FASTA response from NCBI")
    header = lines[0][1:].strip()
    sequence = "".join(lines[1:]).upper()
    if not sequence:
        raise ValueError("Empty FASTA sequence from NCBI")
    return header, sequence


def wrap_sequence(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[index:index + width] for index in range(0, len(sequence), width))


def main() -> int:
    parser = argparse.ArgumentParser(description="Download four seasonal human coronavirus reference genomes")
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("database/virus/seasonal_coronavirus/reference_genomes"),
        help="Output directory",
    )
    args = parser.parse_args()

    out_dir = args.out_dir.expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    combined_fasta_path = out_dir / "seasonal_hcov_reference_genomes.fasta"
    manifest_path = out_dir / "seasonal_hcov_reference_genomes_manifest.tsv"

    manifest_rows: list[dict[str, str]] = []
    combined_entries: list[str] = []

    for record in REFERENCE_GENOMES:
        accession = record["accession"]
        virus = record["virus"]
        fasta_text = fetch_fasta(accession)
        header, sequence = parse_fasta(fasta_text)
        fasta_entry = f">{header}\n{wrap_sequence(sequence)}\n"
        per_accession_path = out_dir / f"{accession}.fasta"
        per_accession_path.write_text(fasta_entry, encoding="utf-8")
        combined_entries.append(fasta_entry.rstrip())
        manifest_rows.append(
            {
                "virus": virus,
                "species": record["species"],
                "accession": accession,
                "sequence_length": str(len(sequence)),
                "header": header,
                "fasta_path": per_accession_path.name,
                "source": f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}",
            }
        )
        print(f"[downloaded] {virus} {accession}", file=sys.stderr)

    combined_fasta_path.write_text("\n".join(combined_entries) + "\n", encoding="utf-8")
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["virus", "species", "accession", "sequence_length", "header", "fasta_path", "source"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(manifest_rows)

    print(f"[ok] wrote {combined_fasta_path}", file=sys.stderr)
    print(f"[ok] wrote {manifest_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
