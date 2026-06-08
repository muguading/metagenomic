from __future__ import annotations

import argparse
import csv
import http.client
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path


SPECIES = [
    ("A", "Rhinovirus A"),
    ("B", "Rhinovirus B"),
    ("C", "Rhinovirus C"),
]

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def _request(url: str, retries: int = 5, sleep_seconds: float = 1.5) -> str:
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        request = urllib.request.Request(
            url,
            headers={"User-Agent": "metagenomic-rhinovirus-downloader/1.0"},
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


def _batched(items: list[str], size: int) -> list[list[str]]:
    return [items[index:index + size] for index in range(0, len(items), size)]


def esearch_ids(term: str, email: str = "", api_key: str = "") -> list[str]:
    params = {
        "db": "nuccore",
        "term": term,
        "retmax": "10000",
        "retmode": "json",
    }
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key
    url = f"{EUTILS_BASE}/esearch.fcgi?{urllib.parse.urlencode(params)}"
    payload = json.loads(_request(url))
    id_list = payload.get("esearchresult", {}).get("idlist", [])
    return [str(item).strip() for item in id_list if str(item).strip()]


def efetch_fasta(id_batch: list[str], email: str = "", api_key: str = "") -> str:
    params = {
        "db": "nuccore",
        "id": ",".join(id_batch),
        "rettype": "fasta",
        "retmode": "text",
    }
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key
    url = f"{EUTILS_BASE}/efetch.fcgi?{urllib.parse.urlencode(params)}"
    return _request(url)


def _extract_type_label(header: str, species_group: str) -> str:
    text = str(header or "")
    patterns = [
        rf"human rhinovirus {species_group}(?: type)?[-\s]*([0-9]+)",
        rf"rhinovirus {species_group}(?: type)?[-\s]*([0-9]+)",
        rf"\bhrv[-\s]?{species_group}[-\s]*([0-9]+)\b",
        rf"\brv[-\s]?{species_group}[-\s]*([0-9]+)\b",
        rf"\b{species_group}[-\s]*([0-9]+)\b",
    ]
    for pattern in patterns:
        matched = re.search(pattern, text, flags=re.IGNORECASE)
        if matched:
            return matched.group(1)
    return ""


def _parse_fasta_records(fasta_text: str, species_group: str, species_name: str) -> tuple[list[dict[str, str]], list[str]]:
    rows: list[dict[str, str]] = []
    fasta_entries: list[str] = []
    current_header = ""
    current_seq_lines: list[str] = []

    def _flush() -> None:
        nonlocal current_header, current_seq_lines
        if not current_header or not current_seq_lines:
            return
        sequence = "".join(line.strip() for line in current_seq_lines).upper()
        accession = current_header.split()[0].lstrip(">").strip()
        if len(sequence) < 6000 or len(sequence) > 9000:
            current_header = ""
            current_seq_lines = []
            return
        type_label = _extract_type_label(current_header, species_group)
        rows.append(
            {
                "species_group": species_group,
                "species_name": species_name,
                "accession": accession,
                "type_label": type_label,
                "sequence_length": str(len(sequence)),
                "header": current_header.lstrip(">"),
                "title": " ".join(current_header.split()[1:]).strip(),
            }
        )
        wrapped = "\n".join(sequence[index:index + 80] for index in range(0, len(sequence), 80))
        fasta_entries.append(f"{current_header}\n{wrapped}")
        current_header = ""
        current_seq_lines = []

    for line in fasta_text.splitlines():
        if line.startswith(">"):
            _flush()
            current_header = line.strip()
            current_seq_lines = []
        elif current_header:
            current_seq_lines.append(line.strip())
    _flush()
    return rows, fasta_entries


def main() -> int:
    parser = argparse.ArgumentParser(description="Download human rhinovirus A/B/C complete genomes from NCBI")
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("database/virus/rhinovirus/full_genomes"),
        help="Output directory",
    )
    parser.add_argument("--email", default="", help="Optional email for NCBI E-utilities")
    parser.add_argument("--api-key", default="", help="Optional NCBI API key")
    args = parser.parse_args()

    out_dir = args.out_dir.expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = out_dir / "human_rhinovirus_A_B_C_manifest.tsv"
    combined_fasta_path = out_dir / "human_rhinovirus_A_B_C_complete_genomes.fasta"

    manifest_rows: list[dict[str, str]] = []
    combined_entries: list[str] = []
    seen_accessions: set[str] = set()

    for species_group, species_name in SPECIES:
        species_fasta_path = out_dir / f"human_rhinovirus_{species_group}_complete_genomes.fasta"
        if species_fasta_path.exists() and species_fasta_path.stat().st_size > 0:
            fasta_text = species_fasta_path.read_text(encoding="utf-8", errors="ignore")
            species_rows, species_entries = _parse_fasta_records(fasta_text, species_group, species_name)
            if species_rows:
                for row, entry in zip(species_rows, species_entries):
                    accession = row["accession"]
                    if accession in seen_accessions:
                        continue
                    seen_accessions.add(accession)
                    manifest_rows.append(row)
                    combined_entries.append(entry)
                print(f"[skip] {species_name}: reuse existing {species_fasta_path.name} ({len(species_rows)} records)", file=sys.stderr)
                continue

        term = (
            f"\"{species_name}\"[Organism] "
            f"AND \"complete genome\"[Title] "
            f"NOT patent[All Fields]"
        )
        print(f"[search] {species_name}", file=sys.stderr)
        ids = esearch_ids(term, email=args.email, api_key=args.api_key)
        time.sleep(0.34)
        if not ids:
            print(f"[warn] no records found for {species_name}", file=sys.stderr)
            continue

        species_rows: list[dict[str, str]] = []
        species_entries: list[str] = []
        for batch in _batched(ids, 200):
            fasta_text = efetch_fasta(batch, email=args.email, api_key=args.api_key)
            time.sleep(0.34)
            batch_rows, batch_entries = _parse_fasta_records(fasta_text, species_group, species_name)
            species_rows.extend(batch_rows)
            species_entries.extend(batch_entries)

        if not species_rows:
            continue

        unique_rows: list[dict[str, str]] = []
        unique_entries: list[str] = []
        local_seen: set[str] = set()
        for row, entry in zip(species_rows, species_entries):
            accession = row["accession"]
            if accession in local_seen:
                continue
            local_seen.add(accession)
            unique_rows.append(row)
            unique_entries.append(entry)

        unique_rows.sort(
            key=lambda row: (
                int(row.get("type_label") or 999999) if str(row.get("type_label") or "").isdigit() else 999999,
                row.get("accession") or "",
            )
        )
        species_fasta_path.write_text("\n".join(unique_entries).strip() + "\n", encoding="utf-8")
        for row, entry in zip(unique_rows, unique_entries):
            accession = row["accession"]
            if accession in seen_accessions:
                continue
            seen_accessions.add(accession)
            manifest_rows.append(row)
            combined_entries.append(entry)
        print(f"[done] {species_name}: {len(unique_rows)} records", file=sys.stderr)

    if not manifest_rows:
        raise SystemExit("No rhinovirus complete genomes were downloaded.")

    manifest_rows.sort(
        key=lambda row: (
            row.get("species_group") or "",
            int(row.get("type_label") or 999999) if str(row.get("type_label") or "").isdigit() else 999999,
            row.get("accession") or "",
        )
    )
    combined_fasta_path.write_text("\n".join(combined_entries).strip() + "\n", encoding="utf-8")
    manifest_columns = ["species_group", "species_name", "accession", "type_label", "sequence_length", "header", "title"]
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=manifest_columns, delimiter="\t")
        writer.writeheader()
        for row in manifest_rows:
            writer.writerow({column: row.get(column, "") for column in manifest_columns})

    print(f"[ok] wrote {combined_fasta_path}", file=sys.stderr)
    print(f"[ok] wrote {manifest_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
