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
DEFAULT_XLSX = Path("database/virus/orthohantavirus/typing.xlsx")
DEFAULT_OUT_DIR = Path("database/virus/orthohantavirus/reference_genomes")
REQUEST_SLEEP = 0.34


def _request(url: str, retries: int = 6, sleep_seconds: float = 1.5) -> str:
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        request = urllib.request.Request(
            url,
            headers={"User-Agent": "metagenomic-orthohantavirus-downloader/1.0"},
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


def _clean_field(value: object) -> str:
    return str(value or "").strip()


def _primary_abbrev(abbrev: str, accession: str) -> str:
    parts = [part.strip() for part in abbrev.split(";") if part.strip()]
    primary = parts[0] if parts else accession
    return re.sub(r"[^A-Za-z0-9._-]+", "_", primary)


def load_manifest_rows(xlsx_path: Path) -> list[dict[str, str]]:
    workbook = openpyxl.load_workbook(xlsx_path, data_only=True)
    worksheet = workbook.active

    current_meta: dict[str, str] | None = None
    rows: list[dict[str, str]] = []
    for raw_row in worksheet.iter_rows(min_row=2, values_only=True):
        genus, species, virus_name, isolate, accession_cell, available_sequence, abbrev = raw_row[:7]
        genus_text = _clean_field(genus)
        accession_text = _clean_field(accession_cell)
        if genus_text and genus_text != "Orthohantavirus":
            current_meta = None
            continue
        if any(value not in (None, "") for value in (species, virus_name, isolate, available_sequence, abbrev)):
            current_meta = {
                "genus": genus_text,
                "species": _clean_field(species),
                "virus_name": _clean_field(virus_name),
                "isolate": _clean_field(isolate),
                "available_sequence": _clean_field(available_sequence),
                "abbrev": _clean_field(abbrev),
            }
        if not current_meta or not accession_text:
            continue
        if not re.match(r"^[LMS]:\s*[A-Za-z0-9_.]+\s*;?$", accession_text):
            continue
        segment, accession = _parse_accession_cell(accession_text)
        row = dict(current_meta)
        row["segment"] = segment
        row["accession"] = accession
        row["accession_root"] = accession.split(".")[0]
        row["abbrev_primary"] = _primary_abbrev(row["abbrev"], row["accession_root"])
        row["record_id"] = f"{row['abbrev_primary']}_{segment}_{row['accession_root']}"
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


def fetch_gff3(accession: str, email: str = "", api_key: str = "") -> str:
    params = {
        "id": accession,
        "db": "nuccore",
        "report": "gff3",
        "retmode": "text",
    }
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key
    url = f"{SVIEWER_URL}?{urllib.parse.urlencode(params)}"
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
        row["record_id"],
        f"abbrev={row['abbrev']}",
        f"segment={row['segment']}",
        f"accession={row['accession']}",
    ]
    if row["species"]:
        parts.append(f"species={row['species']}")
    if row["virus_name"]:
        parts.append(f"virus_name={row['virus_name']}")
    if row["isolate"]:
        parts.append(f"isolate={row['isolate']}")
    parts.append(f"source_header={original_header}")
    return " | ".join(parts)


def normalize_gff_text(gff_text: str, seqid: str) -> tuple[str, int]:
    output_lines: list[str] = []
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
            output_lines.append(line)
            continue
        columns = line.split("\t")
        if len(columns) != 9:
            continue
        columns[0] = seqid
        output_lines.append("\t".join(columns))
        feature_count += 1
    normalized_text = "\n".join(output_lines).strip()
    if normalized_text:
        normalized_text += "\n"
    return normalized_text, feature_count


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Download Orthohantavirus typing reference segments and matching GFF3 annotations from typing.xlsx"
    )
    parser.add_argument("--xlsx", type=Path, default=DEFAULT_XLSX, help="Path to orthohantavirus typing workbook")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR, help="Output directory")
    parser.add_argument("--email", default="", help="Optional email for NCBI E-utilities")
    parser.add_argument("--api-key", default="", help="Optional NCBI API key")
    args = parser.parse_args()

    xlsx_path = args.xlsx.expanduser().resolve()
    out_dir = args.out_dir.expanduser().resolve()
    gff_raw_dir = out_dir / "gff3"
    gff_norm_dir = out_dir / "gff3_normalized"
    out_dir.mkdir(parents=True, exist_ok=True)
    gff_raw_dir.mkdir(parents=True, exist_ok=True)
    gff_norm_dir.mkdir(parents=True, exist_ok=True)

    manifest_rows = load_manifest_rows(xlsx_path)
    if not manifest_rows:
        raise SystemExit(f"No accession rows found in {xlsx_path}")

    combined_fasta_path = out_dir / "orthohantavirus_typing_reference_segments.fasta"
    manifest_path = out_dir / "orthohantavirus_typing_reference_segments_manifest.tsv"

    combined_entries: list[str] = []
    final_rows: list[dict[str, str]] = []
    failures: list[tuple[str, str]] = []

    for row in manifest_rows:
        type_dir = out_dir / row["abbrev_primary"]
        type_dir.mkdir(parents=True, exist_ok=True)
        fasta_path = type_dir / f"{row['record_id']}.fasta"
        gff_raw_path = gff_raw_dir / f"{row['record_id']}.gff3"
        gff_norm_path = gff_norm_dir / f"{row['record_id']}.gff3"

        row_copy = dict(row)
        row_copy.update(
            {
                "sequence_length": "",
                "fasta_path": str(fasta_path),
                "gff_raw_path": str(gff_raw_path),
                "gff_path": str(gff_norm_path),
                "header": "",
                "gff_feature_count": "",
                "fasta_status": "",
                "gff_status": "",
                "status": "ok",
                "note": "",
            }
        )

        try:
            fasta_status = "cached"
            if fasta_path.exists() and fasta_path.stat().st_size > 0:
                fasta_text = fasta_path.read_text(encoding="utf-8", errors="ignore")
                original_header, sequence = parse_single_fasta(fasta_text)
                rewritten_header = original_header
                if not original_header.startswith(f"{row['record_id']} |"):
                    rewritten_header = build_header(row, original_header)
                    fasta_text = f">{rewritten_header}\n{wrap_fasta_sequence(sequence)}\n"
                    fasta_path.write_text(fasta_text, encoding="utf-8")
                    fasta_status = "rewritten"
            else:
                raw_fasta = efetch_fasta(row["accession"], email=args.email, api_key=args.api_key)
                time.sleep(REQUEST_SLEEP)
                original_header, sequence = parse_single_fasta(raw_fasta)
                rewritten_header = build_header(row, original_header)
                fasta_text = f">{rewritten_header}\n{wrap_fasta_sequence(sequence)}\n"
                fasta_path.write_text(fasta_text, encoding="utf-8")
                fasta_status = "downloaded"

            gff_status = "cached"
            raw_gff_text = ""
            if gff_raw_path.exists() and gff_raw_path.stat().st_size > 0:
                raw_gff_text = gff_raw_path.read_text(encoding="utf-8", errors="ignore")
            else:
                raw_gff_text = fetch_gff3(row["accession"], email=args.email, api_key=args.api_key).strip()
                time.sleep(REQUEST_SLEEP)
                if "##gff-version" in raw_gff_text:
                    gff_raw_path.write_text(raw_gff_text + "\n", encoding="utf-8")
                    gff_status = "downloaded"
                else:
                    raw_gff_text = ""
                    gff_status = "missing"

            feature_count = 0
            if raw_gff_text and "##gff-version" in raw_gff_text:
                normalized_gff_text, feature_count = normalize_gff_text(raw_gff_text, row["record_id"])
                if normalized_gff_text:
                    gff_norm_path.write_text(normalized_gff_text, encoding="utf-8")
                    if gff_status == "cached" and not gff_norm_path.exists():
                        gff_status = "normalized"
                else:
                    gff_status = "empty_after_normalization"
            else:
                if gff_norm_path.exists():
                    gff_norm_path.unlink()

            combined_entries.append(fasta_text.strip())
            row_copy.update(
                {
                    "sequence_length": str(len(sequence)),
                    "header": rewritten_header,
                    "gff_feature_count": str(feature_count),
                    "fasta_status": fasta_status,
                    "gff_status": gff_status,
                    "status": "ok" if gff_status != "missing" else "gff_missing",
                    "note": "" if gff_status != "missing" else "NCBI did not return a valid GFF3 record",
                }
            )
            final_rows.append(row_copy)
            print(
                f"[done] {row['abbrev_primary']} {row['segment']} {row['accession']} "
                f"(fasta={fasta_status}, gff={gff_status})",
                file=sys.stderr,
            )
        except Exception as exc:  # noqa: BLE001
            row_copy["status"] = "error"
            row_copy["note"] = str(exc)
            final_rows.append(row_copy)
            failures.append((row["record_id"], str(exc)))
            print(f"[fail] {row['record_id']}: {exc}", file=sys.stderr)

    combined_fasta_path.write_text("\n".join(combined_entries).strip() + "\n", encoding="utf-8")

    columns = [
        "genus",
        "species",
        "virus_name",
        "isolate",
        "abbrev",
        "abbrev_primary",
        "available_sequence",
        "segment",
        "accession",
        "accession_root",
        "record_id",
        "sequence_length",
        "gff_feature_count",
        "fasta_status",
        "gff_status",
        "status",
        "fasta_path",
        "gff_raw_path",
        "gff_path",
        "header",
        "note",
    ]
    final_rows.sort(key=lambda item: (item["abbrev_primary"], item["segment"], item["accession"]))
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for row in final_rows:
            writer.writerow({column: row.get(column, "") for column in columns})

    summary_rows: dict[str, list[str]] = {}
    for row in final_rows:
        if row.get("status") != "error":
            summary_rows.setdefault(row["abbrev_primary"], []).append(row["segment"])
    for abbrev, segments in sorted(summary_rows.items()):
        ordered = ",".join(sorted(segments))
        print(f"[summary] {abbrev}: {ordered}", file=sys.stderr)

    print(f"[ok] wrote {combined_fasta_path}", file=sys.stderr)
    print(f"[ok] wrote {manifest_path}", file=sys.stderr)

    if failures:
        print(f"[warn] {len(failures)} records failed", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
