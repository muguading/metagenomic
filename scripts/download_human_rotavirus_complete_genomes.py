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
    (
        "A",
        "Human rotavirus A",
        '("Human rotavirus A"[Organism] OR "Rotavirus A"[Organism]) AND (human[All Fields] OR Homo sapiens[Organism]) AND ("complete genome"[Title] OR "complete coding sequence"[Title] OR "complete cds"[Title]) NOT patent[All Fields]',
    ),
    (
        "B",
        "Human rotavirus B",
        '("Human rotavirus B"[Organism] OR "Rotavirus B"[Organism]) AND (human[All Fields] OR Homo sapiens[Organism]) AND ("complete genome"[Title] OR "complete coding sequence"[Title] OR "complete cds"[Title]) NOT patent[All Fields]',
    ),
    (
        "C",
        "Human rotavirus C",
        '("Human rotavirus C"[Organism] OR "Rotavirus C"[Organism]) AND (human[All Fields] OR Homo sapiens[Organism]) AND ("complete genome"[Title] OR "complete coding sequence"[Title] OR "complete cds"[Title]) NOT patent[All Fields]',
    ),
]

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def _request(url: str, retries: int = 5, sleep_seconds: float = 1.5) -> str:
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        request = urllib.request.Request(
            url,
            headers={"User-Agent": "metagenomic-rotavirus-downloader/1.0"},
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
        "retmax": "20000",
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


def _looks_like_human_rotavirus(header: str) -> bool:
    text = str(header or "").lower()
    excluded = [
        "porcine",
        "swine",
        "bovine",
        "ovine",
        "equine",
        "murine",
        "rat ",
        "rabbit",
        "feline",
        "canine",
        "simian",
        "bat ",
        "avian",
        "chicken",
        "turkey",
        "deer",
        "pangolin",
        "sea lion",
        "shrew",
        "marmot",
        "rodent",
    ]
    if any(token in text for token in excluded):
        return False
    return "human rotavirus" in text or "homo sapiens" in text or "/hu/" in text or "/human/" in text


def _extract_g_genotype(header: str) -> str:
    text = str(header or "")
    patterns = [
        r"\b(G\d{1,2})P\[\d{1,2}\]\b",
        r"\b(G\d{1,2})\b",
    ]
    for pattern in patterns:
        matched = re.search(pattern, text, flags=re.IGNORECASE)
        if matched:
            return matched.group(1).upper()
    return ""


def _extract_p_genotype(header: str) -> str:
    text = str(header or "")
    patterns = [
        r"\bG\d{1,2}(P\[\d{1,2}\])\b",
        r"\b(P\[\d{1,2}\])\b",
    ]
    for pattern in patterns:
        matched = re.search(pattern, text, flags=re.IGNORECASE)
        if matched:
            return matched.group(1).upper()
    return ""


def _extract_segment_label(header: str) -> str:
    text = str(header or "")
    patterns = [
        r"\bsegment\s+([1-9]|10|11)\b",
        r"\b(VP[1-7]|NSP[1-6])\b",
    ]
    for pattern in patterns:
        matched = re.search(pattern, text, flags=re.IGNORECASE)
        if matched:
            value = matched.group(1).upper()
            if value.isdigit():
                return f"segment_{value}"
            return value
    return ""


def _extract_isolate_label(header: str) -> str:
    text = str(header or "")
    patterns = [
        r"\bRVA/[A-Za-z0-9_\-./]+",
        r"\bRVB/[A-Za-z0-9_\-./]+",
        r"\bRVC/[A-Za-z0-9_\-./]+",
        r"\b(Hu/[A-Za-z0-9_\-./]+)",
        r"\bisolate\s+([A-Za-z0-9_\-./]+)",
        r"\bstrain\s+([A-Za-z0-9_\-./]+)",
    ]
    for pattern in patterns:
        matched = re.search(pattern, text, flags=re.IGNORECASE)
        if matched:
            return matched.group(0 if pattern.startswith(r"\bR") else 1).strip()
    return ""


def _parse_fasta_records(
    fasta_text: str,
    species_group: str,
    species_name: str,
) -> tuple[list[dict[str, str]], list[str]]:
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
        if len(sequence) < 500:
            current_header = ""
            current_seq_lines = []
            return
        if not _looks_like_human_rotavirus(current_header):
            current_header = ""
            current_seq_lines = []
            return
        rows.append(
            {
                "species_group": species_group,
                "species_name": species_name,
                "accession": accession,
                "g_genotype": _extract_g_genotype(current_header) if species_group == "A" else "",
                "p_genotype": _extract_p_genotype(current_header) if species_group == "A" else "",
                "segment_label": _extract_segment_label(current_header),
                "isolate_label": _extract_isolate_label(current_header),
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
    parser = argparse.ArgumentParser(description="Download human rotavirus A/B/C complete genome records from NCBI")
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("database/virus/Rotavirus/full_genomes"),
        help="Output directory",
    )
    parser.add_argument("--email", default="", help="Optional email for NCBI E-utilities")
    parser.add_argument("--api-key", default="", help="Optional NCBI API key")
    args = parser.parse_args()

    out_dir = args.out_dir.expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = out_dir / "human_rotavirus_A_B_C_manifest.tsv"
    combined_fasta_path = out_dir / "human_rotavirus_A_B_C_complete_genomes.fasta"

    manifest_rows: list[dict[str, str]] = []
    combined_entries: list[str] = []
    seen_accessions: set[str] = set()

    for species_group, species_name, term in SPECIES:
        species_fasta_path = out_dir / f"human_rotavirus_{species_group}_complete_genomes.fasta"
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
            print(f"[warn] {species_name}: fetched ids but no human complete-genome-like records survived filtering", file=sys.stderr)
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
                row.get("g_genotype") or "ZZZ",
                row.get("p_genotype") or "ZZZ",
                row.get("segment_label") or "ZZZ",
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
        raise SystemExit("No human rotavirus records were downloaded.")

    manifest_rows.sort(
        key=lambda row: (
            row.get("species_group") or "",
            row.get("g_genotype") or "ZZZ",
            row.get("p_genotype") or "ZZZ",
            row.get("segment_label") or "ZZZ",
            row.get("accession") or "",
        )
    )
    combined_fasta_path.write_text("\n".join(combined_entries).strip() + "\n", encoding="utf-8")
    manifest_columns = [
        "species_group",
        "species_name",
        "accession",
        "g_genotype",
        "p_genotype",
        "segment_label",
        "isolate_label",
        "sequence_length",
        "header",
        "title",
    ]
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
