#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import http.client
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import openpyxl


EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
SVIEWER_URL = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi"
DEFAULT_XLSX = Path("database/virus/Orthoebolavirus/typing.xlsx")
DEFAULT_OUT_DIR = Path("database/virus/Orthoebolavirus/reference_genomes")
REQUEST_SLEEP = 0.34


def request_text(url: str, retries: int = 5) -> str:
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        request = urllib.request.Request(url, headers={"User-Agent": "metagenomic-orthoebolavirus-downloader/1.0"})
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                return response.read().decode("utf-8", errors="ignore")
        except (urllib.error.URLError, TimeoutError, http.client.IncompleteRead, ConnectionError) as exc:
            last_error = exc
            if attempt == retries:
                break
            time.sleep(1.5 * attempt)
    assert last_error is not None
    raise last_error


def clean(value: object) -> str:
    return str(value or "").strip()


def safe_id(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())


def load_rows(xlsx_path: Path) -> list[dict[str, str]]:
    workbook = openpyxl.load_workbook(xlsx_path, data_only=True)
    worksheet = workbook.active
    rows: list[dict[str, str]] = []
    for raw_row in worksheet.iter_rows(min_row=2, values_only=True):
        genus, species, virus_name, isolate, accession, available, abbrev = raw_row[:7]
        if not clean(accession):
            continue
        row = {
            "genus": clean(genus),
            "species": clean(species),
            "virus_name": clean(virus_name),
            "isolate": clean(isolate),
            "accession": clean(accession).split(".")[0],
            "available_sequence": clean(available),
            "abbrev": clean(abbrev),
        }
        if row["genus"] != "Orthoebolavirus":
            continue
        row["record_id"] = f"{safe_id(row['abbrev'])}_{row['accession']}"
        rows.append(row)
    return rows


def fetch_fasta(accession: str, email: str = "", api_key: str = "") -> str:
    params = {"db": "nuccore", "id": accession, "rettype": "fasta", "retmode": "text"}
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key
    return request_text(f"{EUTILS_BASE}/efetch.fcgi?{urllib.parse.urlencode(params)}")


def fetch_gff3(accession: str, email: str = "", api_key: str = "") -> str:
    params = {"id": accession, "db": "nuccore", "report": "gff3", "retmode": "text"}
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key
    return request_text(f"{SVIEWER_URL}?{urllib.parse.urlencode(params)}")


def parse_fasta(fasta_text: str) -> tuple[str, str]:
    header = ""
    sequence_parts: list[str] = []
    for line in fasta_text.splitlines():
        if line.startswith(">"):
            if header:
                raise ValueError("Expected one FASTA record, got multiple records")
            header = line[1:].strip()
        elif header:
            sequence_parts.append(line.strip())
    sequence = "".join(sequence_parts).upper()
    if not header or not sequence:
        raise ValueError("Downloaded FASTA is empty")
    return header, sequence


def wrap_sequence(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[index:index + width] for index in range(0, len(sequence), width))


def build_header(row: dict[str, str], source_header: str) -> str:
    parts = [
        row["record_id"],
        f"abbrev={row['abbrev']}",
        f"accession={row['accession']}",
        f"species={row['species']}",
        f"virus_name={row['virus_name']}",
    ]
    if row["isolate"]:
        parts.append(f"isolate={row['isolate']}")
    parts.append(f"source_header={source_header}")
    return " | ".join(parts)


def normalize_gff(gff_text: str, seqid: str) -> tuple[str, int]:
    lines: list[str] = []
    feature_count = 0
    for raw_line in gff_text.splitlines():
        line = raw_line.rstrip("\n")
        if not line:
            continue
        if line.startswith("#"):
            if line.startswith("##sequence-region"):
                parts = line.split()
                if len(parts) >= 4:
                    line = f"##sequence-region {seqid} {parts[-2]} {parts[-1]}"
            lines.append(line)
            continue
        columns = line.split("\t")
        if len(columns) != 9:
            continue
        columns[0] = seqid
        lines.append("\t".join(columns))
        feature_count += 1
    text = "\n".join(lines).strip()
    return (text + "\n" if text else ""), feature_count


def main() -> int:
    parser = argparse.ArgumentParser(description="Download Orthoebolavirus genomes and GFF3 annotations from typing.xlsx")
    parser.add_argument("--xlsx", type=Path, default=DEFAULT_XLSX)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--email", default="")
    parser.add_argument("--api-key", default="")
    args = parser.parse_args()

    out_dir = args.out_dir.resolve()
    fasta_dir = out_dir / "fasta"
    gff_dir = out_dir / "gff3"
    gff_norm_dir = out_dir / "gff3_normalized"
    for directory in (fasta_dir, gff_dir, gff_norm_dir):
        directory.mkdir(parents=True, exist_ok=True)

    rows = load_rows(args.xlsx)
    manifest_rows: list[dict[str, str]] = []
    combined_entries: list[str] = []

    for row in rows:
        fasta_path = fasta_dir / f"{row['record_id']}.fasta"
        gff_path = gff_dir / f"{row['record_id']}.gff3"
        gff_norm_path = gff_norm_dir / f"{row['record_id']}.gff3"

        fasta_status = "cached"
        if fasta_path.exists() and fasta_path.stat().st_size > 0:
            source_header, sequence = parse_fasta(fasta_path.read_text(encoding="utf-8", errors="ignore"))
            header = source_header
        else:
            source_header, sequence = parse_fasta(fetch_fasta(row["accession"], email=args.email, api_key=args.api_key))
            time.sleep(REQUEST_SLEEP)
            header = build_header(row, source_header)
            fasta_path.write_text(f">{header}\n{wrap_sequence(sequence)}\n", encoding="utf-8")
            fasta_status = "downloaded"

        gff_status = "cached"
        if gff_path.exists() and gff_path.stat().st_size > 0:
            gff_text = gff_path.read_text(encoding="utf-8", errors="ignore")
        else:
            gff_text = fetch_gff3(row["accession"], email=args.email, api_key=args.api_key).strip()
            time.sleep(REQUEST_SLEEP)
            if "##gff-version" not in gff_text:
                raise RuntimeError(f"NCBI did not return valid GFF3 for {row['accession']}")
            gff_path.write_text(gff_text + "\n", encoding="utf-8")
            gff_status = "downloaded"

        normalized_gff, feature_count = normalize_gff(gff_text, row["record_id"])
        gff_norm_path.write_text(normalized_gff, encoding="utf-8")
        combined_entries.append(f">{header}\n{wrap_sequence(sequence)}")

        manifest_row = dict(row)
        manifest_row.update(
            {
                "sequence_length": str(len(sequence)),
                "gff_feature_count": str(feature_count),
                "fasta_status": fasta_status,
                "gff_status": gff_status,
                "fasta_path": str(fasta_path),
                "gff_path": str(gff_norm_path),
                "header": header,
            }
        )
        manifest_rows.append(manifest_row)
        print(f"[done] {row['abbrev']} {row['accession']} fasta={fasta_status} gff={gff_status}", file=sys.stderr)

    (out_dir / "orthoebolavirus_typing_reference_genomes.fasta").write_text("\n".join(combined_entries) + "\n", encoding="utf-8")
    manifest_path = out_dir / "orthoebolavirus_typing_reference_genomes_manifest.tsv"
    columns = [
        "genus",
        "species",
        "virus_name",
        "isolate",
        "abbrev",
        "available_sequence",
        "accession",
        "record_id",
        "sequence_length",
        "gff_feature_count",
        "fasta_status",
        "gff_status",
        "fasta_path",
        "gff_path",
        "header",
    ]
    with manifest_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(manifest_rows)
    print(f"[manifest] {manifest_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
