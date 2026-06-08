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


EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ALL_TYPE_NUMBERS = list(range(1, 114))
MANUAL_REFERENCES = {
    53: ("HAdV-D53", "FJ169625"),
    63: ("HAdV-D63", "JN935766"),
    68: ("HAdV-B68", "JN860678"),
    73: ("HAdV-D73", "KY618676"),
    76: ("HAdV-B76", "KF633445"),
    77: ("HAdV-B77", "KF268328"),
    78: ("HAdV-B78", "KT970441"),
    79: ("HAdV-B79", "LC177352"),
    80: ("HAdV-D80", "KY618679"),
    81: ("HAdV-D81", "AB765926"),
    82: ("HAdV-D82", "LC066535"),
    85: ("HAdV-D85", "LC314153"),
    86: ("HAdV-D86", "KX868297"),
    87: ("HAdV-D87", "MF476841"),
    91: ("HAdV-D91", "KF268208"),
    93: ("HAdV-D93", "KF268334"),
    94: ("HAdV-D94", "KF268201"),
    96: ("HAdV-D96", "KF268327"),
    97: ("HAdV-D97", "KF268320"),
    99: ("HAdV-D99", "KF268211"),
    104: ("HAdV-C104", "MH558113"),
    106: ("HAdV-B106", "ON393912"),
    108: ("HAdV-C108", "ON054624"),
    109: ("HAdV-D109", "OM830314"),
    111: ("HAdV-D111", "LC652931"),
    112: ("HAdV-D112", "OQ679041"),
    113: ("HAdV-D113", "MW694832"),
}


def _request(url: str, retries: int = 5, sleep_seconds: float = 1.5) -> str:
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        request = urllib.request.Request(
            url,
            headers={"User-Agent": "metagenomic-hadv-supplement/1.0"},
        )
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                return response.read().decode("utf-8", errors="ignore")
        except (urllib.error.URLError, TimeoutError, http.client.IncompleteRead, ConnectionError) as exc:
            last_error = exc
            if attempt == retries:
                break
            print(f"[retry] request failed on attempt {attempt}/{retries}: {exc}", file=sys.stderr)
            time.sleep(sleep_seconds * attempt)
    assert last_error is not None
    raise last_error


def _read_fasta(path: Path) -> dict[str, tuple[str, str]]:
    records: dict[str, tuple[str, str]] = {}
    header = ""
    seq_lines: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header:
                    seq_id = header.split()[0]
                    records[seq_id] = (header, "".join(seq_lines))
                header = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line.strip())
    if header:
        seq_id = header.split()[0]
        records[seq_id] = (header, "".join(seq_lines))
    return records


def _read_manifest(path: Path) -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            accession = row.get("accession", "").strip()
            if accession:
                rows[accession] = row
    return rows


def _extract_type_number(type_label: str) -> int | None:
    matched = re.search(r"(\d+)$", type_label or "")
    return int(matched.group(1)) if matched else None


def _header_matches_type(header: str, type_number: int) -> bool:
    patterns = [
        rf"Human adenovirus(?: type)?\s*{type_number}(?!\d)",
        rf"HAdV-(?:[A-G])?{type_number}(?!\d)",
        rf"hAdV-(?:[A-G])?{type_number}(?!\d)",
    ]
    return any(re.search(pattern, header, flags=re.IGNORECASE) for pattern in patterns)


def _extract_type_label(header: str, accession: str, manifest_row: dict[str, str] | None = None) -> str:
    matched = re.search(r"HAdV-([A-G])(\d+)", header, flags=re.IGNORECASE)
    if matched:
        return f"HAdV-{matched.group(1).upper()}{matched.group(2)}"
    matched = re.search(r"Human mastadenovirus\s+([A-G])", header, flags=re.IGNORECASE)
    if matched:
        number = _extract_type_number(header) or _extract_type_number(accession)
        if number is not None:
            return f"HAdV-{matched.group(1).upper()}{number}"
    if manifest_row:
        group = (manifest_row.get("species_group") or "").strip().upper()
        number = _extract_type_number(header) or _extract_type_number(accession)
        if group and number is not None:
            return f"HAdV-{group}{number}"
    return ""


def _candidate_score(accession: str, header: str, type_number: int) -> tuple[int, int, str]:
    score = 0
    if re.search(rf"HAdV-[A-G]{type_number}(?!\d)", header, flags=re.IGNORECASE):
        score += 100
    if re.search(rf"Human adenovirus(?: type)?\s*{type_number}(?!\d)", header, flags=re.IGNORECASE):
        score += 60
    if re.search(r"\[P\d+H\d+F\d+\]", header):
        score += 20
    if accession.startswith(("NC_", "AC_", "AP_")):
        score += 200
    elif accession.startswith(("JN", "JQ", "HQ", "KF", "OM", "ON", "OQ", "OR", "PP", "PX")):
        score += 40
    return (score, -len(header), accession)


def _wrap_sequence(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[index:index + width] for index in range(0, len(sequence), width))


def _search_ncbi_ids(type_number: int, email: str = "", api_key: str = "") -> list[str]:
    terms = [
        f"\"Human adenovirus {type_number}\"[Title]",
        f"\"Human adenovirus type {type_number}\"[Title]",
        f"\"HAdV-{type_number}\"[Title]",
    ] + [f"\"HAdV-{letter}{type_number}\"[Title]" for letter in "ABCDEFG"]
    params = {
        "db": "nuccore",
        "term": f"({' OR '.join(terms)}) AND \"complete genome\"[Title]",
        "retmax": "20",
        "retmode": "json",
    }
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key
    url = f"{EUTILS_BASE}/esearch.fcgi?{urllib.parse.urlencode(params)}"
    payload = json.loads(_request(url))
    return [item for item in payload.get("esearchresult", {}).get("idlist", []) if item]


def _efetch_fasta_entries(ids: list[str], email: str = "", api_key: str = "") -> list[tuple[str, str]]:
    if not ids:
        return []
    params = {
        "db": "nuccore",
        "id": ",".join(ids),
        "rettype": "fasta",
        "retmode": "text",
    }
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key
    url = f"{EUTILS_BASE}/efetch.fcgi?{urllib.parse.urlencode(params)}"
    text = _request(url).strip()
    entries: list[tuple[str, str]] = []
    header = ""
    seq_lines: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if header:
                entries.append((header, "".join(seq_lines)))
            header = line[1:].strip()
            seq_lines = []
        else:
            seq_lines.append(line.strip())
    if header:
        entries.append((header, "".join(seq_lines)))
    return entries


def _efetch_single_accession(accession: str, email: str = "", api_key: str = "") -> tuple[str, str] | None:
    entries = _efetch_fasta_entries([accession], email=email, api_key=api_key)
    if not entries:
        return None
    return entries[0]


def main() -> int:
    parser = argparse.ArgumentParser(description="Supplement HAdV representative genomes up to all available complete-genome types")
    parser.add_argument(
        "--existing-manifest",
        type=Path,
        default=Path("database/virus/hadv/reference_genomes/hadv_phf_typing_db_representative_genomes.tsv"),
    )
    parser.add_argument(
        "--existing-fasta",
        type=Path,
        default=Path("database/virus/hadv/reference_genomes/hadv_phf_typing_db_representative_genomes.fasta"),
    )
    parser.add_argument(
        "--downloaded-fasta",
        type=Path,
        default=Path("database/virus/hadv/full_genomes/human_mastadenovirus_A_G_complete_genomes.fasta"),
    )
    parser.add_argument(
        "--downloaded-manifest",
        type=Path,
        default=Path("database/virus/hadv/full_genomes/human_mastadenovirus_A_G_manifest.tsv"),
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("database/virus/hadv/reference_genomes"),
    )
    parser.add_argument("--email", default="", help="Optional email for NCBI E-utilities")
    parser.add_argument("--api-key", default="", help="Optional NCBI API key")
    args = parser.parse_args()

    existing_rows = list(_read_manifest(args.existing_manifest.resolve()).values())
    existing_type_numbers = {number for number in (_extract_type_number(row.get("type_label", "")) for row in existing_rows) if number is not None}
    missing_type_numbers = [number for number in ALL_TYPE_NUMBERS if number not in existing_type_numbers]

    downloaded_records = _read_fasta(args.downloaded_fasta.resolve())
    downloaded_manifest_raw = _read_manifest(args.downloaded_manifest.resolve())
    downloaded_manifest = {accession.split(".")[0]: row for accession, row in downloaded_manifest_raw.items()}
    downloaded_by_base_accession = {seq_id.split(".")[0]: (header, sequence) for seq_id, (header, sequence) in downloaded_records.items()}

    supplement_rows: list[dict[str, str]] = []
    supplement_entries: list[str] = []
    unresolved_rows: list[dict[str, str]] = []

    for type_number in missing_type_numbers:
        manual = MANUAL_REFERENCES.get(type_number)
        if manual:
            manual_type_label, manual_accession = manual
            selected_header = ""
            selected_sequence = ""
            source = "downloaded_full_genomes"
            if manual_accession in downloaded_by_base_accession:
                selected_header, selected_sequence = downloaded_by_base_accession[manual_accession]
            else:
                print(f"[fetch] {manual_type_label} <- {manual_accession}", file=sys.stderr)
                fetched = _efetch_single_accession(manual_accession, email=args.email, api_key=args.api_key)
                time.sleep(0.34)
                if fetched is None:
                    unresolved_rows.append(
                        {
                            "type_number": str(type_number),
                            "status": "manual_accession_missing",
                            "note": f"Representative accession {manual_accession} could not be fetched",
                        }
                    )
                    continue
                selected_header, selected_sequence = fetched
                source = "ncbi_efetch"
            supplement_entries.append(f">{selected_header}\n{_wrap_sequence(selected_sequence)}\n")
            supplement_rows.append(
                {
                    "accession": manual_accession,
                    "type_label": manual_type_label,
                    "sequence_length": str(len(selected_sequence)),
                    "source": source,
                    "full_genome_header": selected_header,
                }
            )
            print(f"[add] {manual_type_label} <- {manual_accession}", file=sys.stderr)
            continue

        candidates: list[tuple[tuple[int, int, str], str, str, dict[str, str] | None, str]] = []

        for accession, (header, sequence) in downloaded_by_base_accession.items():
            if _header_matches_type(header, type_number):
                manifest_row = downloaded_manifest.get(accession)
                candidates.append((_candidate_score(accession, header, type_number), accession, header, manifest_row, sequence))

        source = "downloaded_full_genomes"
        if not candidates:
            print(f"[search] type {type_number}", file=sys.stderr)
            ids = _search_ncbi_ids(type_number, email=args.email, api_key=args.api_key)
            time.sleep(0.34)
            entries = _efetch_fasta_entries(ids, email=args.email, api_key=args.api_key)
            time.sleep(0.34)
            for header, sequence in entries:
                accession = header.split()[0].split(".")[0]
                if _header_matches_type(header, type_number):
                    candidates.append((_candidate_score(accession, header, type_number), accession, header, None, sequence))
            source = "ncbi_efetch"

        if not candidates:
            unresolved_rows.append(
                {
                    "type_number": str(type_number),
                    "status": "unresolved",
                    "note": "No complete genome representative found with current exact-title NCBI search",
                }
            )
            continue

        _, accession, header, manifest_row, sequence = max(candidates, key=lambda item: item[0])
        type_label = _extract_type_label(header, str(type_number), manifest_row=manifest_row)
        if not type_label:
            unresolved_rows.append(
                {
                    "type_number": str(type_number),
                    "status": "unresolved_label",
                    "note": f"Matched complete genome but could not infer species letter from header: {accession}",
                }
            )
            continue

        supplement_entries.append(f">{header}\n{_wrap_sequence(sequence)}\n")
        supplement_rows.append(
            {
                "accession": accession,
                "type_label": type_label,
                "sequence_length": str(len(sequence)),
                "source": source,
                "full_genome_header": header,
            }
        )
        print(f"[add] {type_label} <- {accession}", file=sys.stderr)

    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    supplement_fasta = out_dir / "hadv_additional_representative_genomes.fasta"
    supplement_tsv = out_dir / "hadv_additional_representative_genomes.tsv"
    merged_fasta = out_dir / "hadv_representative_genomes_expanded.fasta"
    merged_tsv = out_dir / "hadv_representative_genomes_expanded.tsv"
    unresolved_tsv = out_dir / "hadv_representative_genomes_unresolved.tsv"

    supplement_fasta.write_text("".join(supplement_entries), encoding="utf-8")
    with supplement_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["accession", "type_label", "sequence_length", "source", "full_genome_header"], delimiter="\t")
        writer.writeheader()
        for row in supplement_rows:
            writer.writerow(row)

    existing_fasta_text = args.existing_fasta.resolve().read_text(encoding="utf-8")
    merged_fasta.write_text(existing_fasta_text + "".join(supplement_entries), encoding="utf-8")
    merged_rows = existing_rows + supplement_rows
    with merged_tsv.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = sorted({key for row in merged_rows for key in row.keys()}, key=lambda item: ["accession", "type_label", "sequence_length", "source", "fiber_header", "hexon_header", "penton_header", "full_genome_header"].index(item) if item in ["accession", "type_label", "sequence_length", "source", "fiber_header", "hexon_header", "penton_header", "full_genome_header"] else 100)
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in merged_rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})

    with unresolved_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["type_number", "status", "note"], delimiter="\t")
        writer.writeheader()
        for row in unresolved_rows:
            writer.writerow(row)

    print(f"[ok] wrote {supplement_fasta}", file=sys.stderr)
    print(f"[ok] wrote {supplement_tsv}", file=sys.stderr)
    print(f"[ok] wrote {merged_fasta}", file=sys.stderr)
    print(f"[ok] wrote {merged_tsv}", file=sys.stderr)
    print(f"[ok] wrote {unresolved_tsv}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
