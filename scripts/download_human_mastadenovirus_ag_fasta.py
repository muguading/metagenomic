from __future__ import annotations

import argparse
import csv
import http.client
import json
import re
import sys
import time
import urllib.parse
import urllib.error
import urllib.request
from pathlib import Path


SPECIES = [
    "Human mastadenovirus A",
    "Human mastadenovirus B",
    "Human mastadenovirus C",
    "Human mastadenovirus D",
    "Human mastadenovirus E",
    "Human mastadenovirus F",
    "Human mastadenovirus G",
]

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def _request(url: str, retries: int = 5, sleep_seconds: float = 1.5) -> str:
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        request = urllib.request.Request(
            url,
            headers={
                "User-Agent": "metagenomic-hadv-downloader/1.0",
            },
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


def _parse_species_suffix(species_name: str) -> str:
    return species_name.rsplit(" ", 1)[-1].strip().upper()


def _extract_type_label(header: str) -> str:
    text = str(header or "")
    patterns = [
        r"human adenovirus(?: type)?\s*([0-9]+)",
        r"hadv[-\s]?([0-9]+)",
        r"mastadenovirus [a-z]\s*([0-9]+)",
        r"\btype\s*([0-9]+)\b",
    ]
    for pattern in patterns:
        matched = re.search(pattern, text, flags=re.IGNORECASE)
        if matched:
            return matched.group(1)
    return ""


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


def _parse_fasta_records(fasta_text: str, species_suffix: str, species_name: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    current_header = ""
    current_seq_lines: list[str] = []
    for line in fasta_text.splitlines():
        if line.startswith(">"):
            if current_header and current_seq_lines:
                accession = current_header.split()[0].lstrip(">").strip()
                type_label = _extract_type_label(current_header)
                rows.append(
                    {
                        "species_group": species_suffix,
                        "species_name": species_name,
                        "accession": accession,
                        "type_label": type_label,
                        "header": current_header.lstrip(">"),
                        "sequence_length": str(sum(len(item.strip()) for item in current_seq_lines)),
                        "title": " ".join(current_header.split()[1:]).strip(),
                    }
                )
            current_header = line.strip()
            current_seq_lines = []
        elif current_header:
            current_seq_lines.append(line.strip())
    if current_header and current_seq_lines:
        accession = current_header.split()[0].lstrip(">").strip()
        type_label = _extract_type_label(current_header)
        rows.append(
            {
                "species_group": species_suffix,
                "species_name": species_name,
                "accession": accession,
                "type_label": type_label,
                "header": current_header.lstrip(">"),
                "sequence_length": str(sum(len(item.strip()) for item in current_seq_lines)),
                "title": " ".join(current_header.split()[1:]).strip(),
            }
        )
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description="Download Human mastadenovirus A-G complete genome FASTA from NCBI")
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("database/virus/hadv/full_genomes"),
        help="Output directory",
    )
    parser.add_argument(
        "--email",
        default="",
        help="Optional email for NCBI E-utilities",
    )
    parser.add_argument(
        "--api-key",
        default="",
        help="Optional NCBI API key",
    )
    args = parser.parse_args()

    out_dir = args.out_dir.expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = out_dir / "human_mastadenovirus_A_G_manifest.tsv"
    combined_fasta_path = out_dir / "human_mastadenovirus_A_G_complete_genomes.fasta"

    manifest_rows: list[dict[str, str]] = []
    combined_chunks: list[str] = []

    for species_name in SPECIES:
        species_suffix = _parse_species_suffix(species_name)
        species_fasta_path = out_dir / f"human_mastadenovirus_{species_suffix}_complete_genomes.fasta"
        if species_fasta_path.exists() and species_fasta_path.stat().st_size > 0:
            fasta_text = species_fasta_path.read_text(encoding="utf-8", errors="ignore").strip()
            species_manifest_rows = _parse_fasta_records(fasta_text, species_suffix, species_name)
            if species_manifest_rows:
                manifest_rows.extend(species_manifest_rows)
                combined_chunks.append(fasta_text)
                print(f"[skip] {species_name}: reuse existing {species_fasta_path.name} ({len(species_manifest_rows)} records)", file=sys.stderr)
                continue

        term = f"\"{species_name}\"[Organism] AND \"complete genome\"[Title]"
        print(f"[search] {species_name}", file=sys.stderr)
        ids = esearch_ids(term, email=args.email, api_key=args.api_key)
        time.sleep(0.34)
        if not ids:
            print(f"[warn] no records found for {species_name}", file=sys.stderr)
            continue

        species_manifest_rows: list[dict[str, str]] = []
        species_fasta_chunks: list[str] = []

        for batch in _batched(ids, 200):
            fasta_text = efetch_fasta(batch, email=args.email, api_key=args.api_key)
            time.sleep(0.34)

            species_manifest_rows.extend(_parse_fasta_records(fasta_text, species_suffix, species_name))
            species_fasta_chunks.append(fasta_text.strip())

        if not species_manifest_rows:
            continue
        species_manifest_rows.sort(
            key=lambda row: (
                int(row.get("type_label") or 999999) if str(row.get("type_label") or "").isdigit() else 999999,
                row.get("accession") or "",
            )
        )
        species_fasta_path.write_text("\n".join(chunk for chunk in species_fasta_chunks if chunk).strip() + "\n", encoding="utf-8")
        manifest_rows.extend(species_manifest_rows)
        combined_chunks.extend(chunk for chunk in species_fasta_chunks if chunk)
        print(f"[done] {species_name}: {len(species_manifest_rows)} records", file=sys.stderr)

    if not manifest_rows:
        raise SystemExit("No FASTA records were downloaded.")

    combined_fasta_path.write_text("\n".join(combined_chunks).strip() + "\n", encoding="utf-8")
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
