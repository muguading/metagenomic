#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

from docx import Document


ROOT = Path("/Users/wuhhh/Desktop/徐老师/代码/metagenomic")
DOCX_PATH = ROOT / "database/virus/seasonal_coronavirus/HCoV_OC43_genomes/temi_a_2019560_sm7144.docx"
OUTDIR = ROOT / "database/virus/seasonal_coronavirus/HCoV_OC43_genomes"
METADATA_TSV = OUTDIR / "HCoV_OC43_docx_accessions_expanded.tsv"
MISSING_TSV = OUTDIR / "HCoV_OC43_docx_download_missing.tsv"


def expand_accession(value: str) -> list[str]:
    raw = value.strip()
    if "-" not in raw:
        return [raw]

    left, right = raw.split("-", 1)
    prefix = "".join(ch for ch in left if not ch.isdigit())
    left_num = "".join(ch for ch in left if ch.isdigit())
    right_num = "".join(ch for ch in right if ch.isdigit())
    if not prefix or not left_num or not right_num:
        return [raw]

    width = len(left_num)
    if len(right_num) < width:
        right_num = left_num[: width - len(right_num)] + right_num

    start = int(left_num)
    end = int(right_num)
    if end < start:
        return [raw]

    return [f"{prefix}{number:0{width}d}" for number in range(start, end + 1)]


def parse_docx_rows(path: Path) -> list[dict[str, str]]:
    doc = Document(path)
    table = doc.tables[0]
    rows: list[dict[str, str]] = []

    for row_index, row in enumerate(table.rows[2:], start=3):
        cells = [cell.text.strip().replace("\n", " ") for cell in row.cells]
        if len(cells) < 9:
            continue
        strain_name = cells[6]
        accession_raw = cells[7]
        collection_year = cells[8]
        if not strain_name or not accession_raw:
            continue

        expanded = expand_accession(accession_raw)
        for accession in expanded:
            rows.append(
                {
                    "docx_row": str(row_index),
                    "strain_name": strain_name,
                    "accession_raw": accession_raw,
                    "accession": accession,
                    "collection_year": collection_year,
                    "source_docx": str(path),
                }
            )
    return rows


def ncbi_efetch_text(accession: str, rettype: str) -> str:
    params = urllib.parse.urlencode(
        {
            "db": "nuccore",
            "id": accession,
            "rettype": rettype,
            "retmode": "text",
        }
    )
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?{params}"
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": "metagenomic-codex/1.0",
            "Accept": "text/plain",
        },
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        return response.read().decode("utf-8", errors="replace")


def write_text(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    rows = parse_docx_rows(DOCX_PATH)

    with METADATA_TSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "docx_row",
                "strain_name",
                "accession_raw",
                "accession",
                "collection_year",
                "source_docx",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    missing: list[dict[str, str]] = []
    downloaded_fasta = 0
    downloaded_gff3 = 0

    for row in rows:
        accession = row["accession"]
        fasta_path = OUTDIR / f"{accession}.fasta"
        gff_path = OUTDIR / f"{accession}.gff3"

        if not fasta_path.exists():
            try:
                fasta_text = ncbi_efetch_text(accession, "fasta")
                if not fasta_text.lstrip().startswith(">"):
                    raise ValueError("not_fasta")
                write_text(fasta_path, fasta_text)
                downloaded_fasta += 1
            except (urllib.error.URLError, TimeoutError, ValueError) as exc:
                missing.append(
                    {
                        **row,
                        "artifact": "fasta",
                        "reason": str(exc),
                    }
                )
            time.sleep(0.34)

        if not gff_path.exists():
            try:
                gff_text = ncbi_efetch_text(accession, "gff3")
                if "##gff-version" not in gff_text:
                    raise ValueError("not_gff3")
                write_text(gff_path, gff_text)
                downloaded_gff3 += 1
            except (urllib.error.URLError, TimeoutError, ValueError) as exc:
                missing.append(
                    {
                        **row,
                        "artifact": "gff3",
                        "reason": str(exc),
                    }
                )
            time.sleep(0.34)

    if missing:
        with MISSING_TSV.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "docx_row",
                    "strain_name",
                    "accession_raw",
                    "accession",
                    "collection_year",
                    "source_docx",
                    "artifact",
                    "reason",
                ],
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerows(missing)
    elif MISSING_TSV.exists():
        MISSING_TSV.unlink()

    summary = {
        "rows_in_docx": len(rows),
        "unique_accessions": len({row["accession"] for row in rows}),
        "downloaded_fasta": downloaded_fasta,
        "downloaded_gff3": downloaded_gff3,
        "missing_records": len(missing),
        "metadata_tsv": str(METADATA_TSV),
        "missing_tsv": str(MISSING_TSV) if missing else "",
    }
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
