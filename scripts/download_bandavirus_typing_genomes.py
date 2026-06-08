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
DEFAULT_XLSX = Path("database/virus/bandavirus/typing.xlsx")
DEFAULT_OUT_DIR = Path("database/virus/bandavirus/reference_genomes")


def _request(url: str, retries: int = 5, sleep_seconds: float = 1.5) -> str:
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        request = urllib.request.Request(
            url,
            headers={"User-Agent": "metagenomic-bandavirus-downloader/1.0"},
        )
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                return response.read().decode("utf-8", errors="ignore")
        except (urllib.error.URLError, TimeoutError, http.client.IncompleteRead, ConnectionError) as exc:
            last_error = exc
            if attempt == retries:
                break
            wait_seconds = sleep_seconds * attempt
            print(f"[retry] request failed on attempt {attempt}/{retries}: {exc}", file=sys.stderr)
            time.sleep(wait_seconds)
    assert last_error is not None
    raise last_error


def _parse_accession_cell(value: object) -> tuple[str, str]:
    text = str(value or "").strip()
    match = re.match(r"^([LMS]):\s*([A-Za-z0-9_.]+)\s*;?$", text)
    if not match:
        raise ValueError(f"Unexpected accession cell: {text!r}")
    return match.group(1), match.group(2)


def load_manifest_rows(xlsx_path: Path) -> list[dict[str, str]]:
    wb = openpyxl.load_workbook(xlsx_path, data_only=True)
    ws = wb.active

    current_meta: dict[str, str] | None = None
    rows: list[dict[str, str]] = []
    for raw_row in ws.iter_rows(min_row=2, values_only=True):
        genus, species, virus_name, isolate, accession_cell, available_sequence, abbrev = raw_row[:7]
        if species or virus_name or isolate or available_sequence or abbrev:
            current_meta = {
                "genus": str(genus or "").strip(),
                "species": str(species or "").strip(),
                "virus_name": str(virus_name or "").strip(),
                "isolate": str(isolate or "").strip(),
                "available_sequence": str(available_sequence or "").strip(),
                "abbrev": str(abbrev or "").strip(),
            }
        if not current_meta or not accession_cell:
            continue
        segment, accession = _parse_accession_cell(accession_cell)
        row = dict(current_meta)
        row["segment"] = segment
        row["accession"] = accession
        rows.append(row)
    return rows


def efetch_fasta(accession: str, email: str = "", api_key: str = "") -> str:
    params = {
        "db": "nuccore",
        "id": accession,
        "rettype": "fasta",
        "retmode": "text",
    }
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key
    url = f"{EUTILS_BASE}/efetch.fcgi?{urllib.parse.urlencode(params)}"
    return _request(url)


def parse_single_fasta(fasta_text: str) -> tuple[str, str]:
    header = ""
    seq_lines: list[str] = []
    for line in fasta_text.splitlines():
        if line.startswith(">"):
            if header:
                raise ValueError("Expected one FASTA record, got multiple.")
            header = line[1:].strip()
        elif header:
            seq_lines.append(line.strip())
    sequence = "".join(seq_lines).upper()
    if not header or not sequence:
        raise ValueError("Downloaded FASTA is empty.")
    return header, sequence


def wrap_fasta_sequence(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[index:index + width] for index in range(0, len(sequence), width))


def build_header(row: dict[str, str], original_header: str) -> str:
    parts = [
        row["abbrev"],
        row["segment"],
        row["accession"],
        row["species"],
        row["virus_name"],
    ]
    if row["isolate"]:
        parts.append(f"isolate={row['isolate']}")
    parts.append(f"source_header={original_header}")
    return " | ".join(parts)


def main() -> int:
    parser = argparse.ArgumentParser(description="Download representative bandavirus segment genomes from typing.xlsx")
    parser.add_argument("--xlsx", type=Path, default=DEFAULT_XLSX, help="Path to bandavirus typing workbook")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR, help="Output directory")
    parser.add_argument("--email", default="", help="Optional email for NCBI E-utilities")
    parser.add_argument("--api-key", default="", help="Optional NCBI API key")
    args = parser.parse_args()

    xlsx_path = args.xlsx.expanduser().resolve()
    out_dir = args.out_dir.expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest_rows = load_manifest_rows(xlsx_path)
    if not manifest_rows:
        raise SystemExit(f"No accession rows found in {xlsx_path}")

    combined_fasta_path = out_dir / "bandavirus_typing_reference_segments.fasta"
    manifest_path = out_dir / "bandavirus_typing_reference_segments_manifest.tsv"

    combined_entries: list[str] = []
    final_rows: list[dict[str, str]] = []

    for row in manifest_rows:
        type_dir = out_dir / row["abbrev"]
        type_dir.mkdir(parents=True, exist_ok=True)
        file_name = f"{row['abbrev']}_{row['segment']}_{row['accession']}.fasta"
        fasta_path = type_dir / file_name

        if fasta_path.exists() and fasta_path.stat().st_size > 0:
            fasta_text = fasta_path.read_text(encoding="utf-8", errors="ignore")
            original_header, sequence = parse_single_fasta(fasta_text)
            rewritten_header = original_header
            if not original_header.startswith(f"{row['abbrev']} | {row['segment']} | {row['accession']}"):
                rewritten_header = build_header(row, original_header)
                fasta_text = f">{rewritten_header}\n{wrap_fasta_sequence(sequence)}\n"
                fasta_path.write_text(fasta_text, encoding="utf-8")
            combined_entries.append(fasta_text.strip())
            row_copy = dict(row)
            row_copy["sequence_length"] = str(len(sequence))
            row_copy["fasta_path"] = str(fasta_path)
            row_copy["header"] = rewritten_header
            final_rows.append(row_copy)
            print(f"[skip] reuse {fasta_path.name}", file=sys.stderr)
            continue

        raw_fasta = efetch_fasta(row["accession"], email=args.email, api_key=args.api_key)
        time.sleep(0.34)
        original_header, sequence = parse_single_fasta(raw_fasta)
        rewritten_header = build_header(row, original_header)
        fasta_text = f">{rewritten_header}\n{wrap_fasta_sequence(sequence)}\n"
        fasta_path.write_text(fasta_text, encoding="utf-8")
        combined_entries.append(fasta_text.strip())

        row_copy = dict(row)
        row_copy["sequence_length"] = str(len(sequence))
        row_copy["fasta_path"] = str(fasta_path)
        row_copy["header"] = rewritten_header
        final_rows.append(row_copy)
        print(f"[done] {row['abbrev']} {row['segment']} {row['accession']}", file=sys.stderr)

    combined_fasta_path.write_text("\n".join(combined_entries).strip() + "\n", encoding="utf-8")

    columns = [
        "genus",
        "species",
        "virus_name",
        "isolate",
        "abbrev",
        "available_sequence",
        "segment",
        "accession",
        "sequence_length",
        "fasta_path",
        "header",
    ]
    final_rows.sort(key=lambda item: (item["abbrev"], item["segment"], item["accession"]))
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for row in final_rows:
            writer.writerow({column: row.get(column, "") for column in columns})

    summary_rows: dict[str, list[str]] = {}
    for row in final_rows:
        summary_rows.setdefault(row["abbrev"], []).append(row["segment"])
    for abbrev, segments in sorted(summary_rows.items()):
        ordered = ",".join(sorted(segments))
        print(f"[summary] {abbrev}: {ordered}", file=sys.stderr)

    print(f"[ok] wrote {combined_fasta_path}", file=sys.stderr)
    print(f"[ok] wrote {manifest_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
