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
DISCOVERY_QUERY = '"Human rotavirus A"[Organism] AND ("VP7"[Title] OR "VP4"[Title]) AND "complete cds"[Title]'


def _request(url: str, retries: int = 5, sleep_seconds: float = 1.5) -> str:
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        request = urllib.request.Request(
            url,
            headers={"User-Agent": "metagenomic-rotavirusA-subtype-downloader/1.0"},
        )
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                return response.read().decode("utf-8", errors="ignore")
        except (urllib.error.URLError, TimeoutError, http.client.IncompleteRead, ConnectionError) as exc:
            last_error = exc
            if attempt == retries:
                break
            time.sleep(sleep_seconds * attempt)
    assert last_error is not None
    raise last_error


def _batched(items: list[str], size: int) -> list[list[str]]:
    return [items[index:index + size] for index in range(0, len(items), size)]


def esearch_ids(term: str, retmax: int = 1000, email: str = "", api_key: str = "") -> list[str]:
    params = {
        "db": "nuccore",
        "term": term,
        "retmax": str(retmax),
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


def _extract_combo_and_isolate(header: str) -> tuple[str, str]:
    text = str(header or "")
    combo_match = re.search(r"(G\d{1,2}P\[\d{1,2}\])", text, flags=re.IGNORECASE)
    isolate_match = re.search(
        r"(RVA(?:/[A-Za-z0-9_.-]+){3,}/G\d{1,2}P\[\d{1,2}\]|RVA-[A-Za-z0-9_.-]+_[0-9]{4}_G\d{1,2}P\[\d{1,2}\])",
        text,
        flags=re.IGNORECASE,
    )
    combo = combo_match.group(1).upper() if combo_match else ""
    isolate = isolate_match.group(1) if isolate_match else ""
    return combo, isolate


def _extract_segment_label(header: str) -> str:
    text = str(header or "")
    patterns = [
        r"\b(VP[1-7])\b",
        r"\b(NSP[1-6])\b",
        r"\bsegment\s+([1-9]|10|11)\b",
    ]
    for pattern in patterns:
        matched = re.search(pattern, text, flags=re.IGNORECASE)
        if matched:
            value = matched.group(1).upper()
            if value.isdigit():
                return f"segment_{value}"
            return value
    return ""


def _parse_fasta_records(fasta_text: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    current_header = ""
    current_seq_lines: list[str] = []

    def _flush() -> None:
        nonlocal current_header, current_seq_lines
        if not current_header or not current_seq_lines:
            return
        sequence = "".join(line.strip() for line in current_seq_lines).upper()
        accession = current_header.split()[0].lstrip(">").strip()
        combo, isolate = _extract_combo_and_isolate(current_header)
        rows.append(
            {
                "accession": accession,
                "combo": combo,
                "isolate_label": isolate,
                "segment_label": _extract_segment_label(current_header),
                "sequence_length": str(len(sequence)),
                "header": current_header.lstrip(">"),
                "sequence": sequence,
            }
        )
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
    return rows


def _wrap_sequence(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[index:index + width] for index in range(0, len(sequence), width))


def discover_candidate_isolates(limit: int, email: str = "", api_key: str = "") -> dict[str, list[str]]:
    ids = esearch_ids(DISCOVERY_QUERY, retmax=limit, email=email, api_key=api_key)
    candidates: dict[str, list[str]] = {}
    for batch in _batched(ids, 200):
        records = _parse_fasta_records(efetch_fasta(batch, email=email, api_key=api_key))
        time.sleep(0.34)
        for record in records:
            combo = record["combo"]
            isolate = record["isolate_label"]
            if not combo or not isolate:
                continue
            bucket = candidates.setdefault(combo, [])
            if isolate not in bucket:
                bucket.append(isolate)
    return candidates


def fetch_complete_isolate_records(isolate_label: str, email: str = "", api_key: str = "") -> list[dict[str, str]]:
    term = f'"{isolate_label}"[All Fields] AND "Human rotavirus A"[Organism] AND "complete cds"[Title]'
    ids = esearch_ids(term, retmax=50, email=email, api_key=api_key)
    if not ids:
        return []
    rows: list[dict[str, str]] = []
    for batch in _batched(ids, 50):
        rows.extend(_parse_fasta_records(efetch_fasta(batch, email=email, api_key=api_key)))
        time.sleep(0.34)
    deduped: dict[str, dict[str, str]] = {}
    for row in rows:
        deduped[row["accession"]] = row
    return list(deduped.values())


def choose_best_segment_set(records: list[dict[str, str]]) -> list[dict[str, str]]:
    by_segment: dict[str, list[dict[str, str]]] = {}
    for row in records:
        segment = row.get("segment_label") or ""
        if not segment:
            continue
        by_segment.setdefault(segment, []).append(row)
    chosen: list[dict[str, str]] = []
    for segment, rows in sorted(by_segment.items()):
        rows.sort(
            key=lambda row: (
                -int(row.get("sequence_length") or 0),
                row.get("accession") or "",
            )
        )
        chosen.append(rows[0])
    return chosen


def main() -> int:
    parser = argparse.ArgumentParser(description="Download representative complete segment sets for human rotavirus A G/P subtypes")
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("database/virus/Rotavirus/subtype_complete_genomes"),
        help="Output directory",
    )
    parser.add_argument("--discovery-limit", type=int, default=400, help="How many VP7/VP4 records to scan for subtype discovery")
    parser.add_argument("--email", default="", help="Optional email for NCBI E-utilities")
    parser.add_argument("--api-key", default="", help="Optional NCBI API key")
    args = parser.parse_args()

    out_dir = args.out_dir.expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = out_dir / "rotavirus_a_subtype_complete_genomes_manifest.tsv"
    subtype_manifest_path = out_dir / "rotavirus_a_subtype_complete_genomes_summary.tsv"
    combined_fasta_path = out_dir / "rotavirus_a_subtype_complete_genomes.fasta"

    for stale_file in out_dir.glob("*.fasta"):
        stale_file.unlink(missing_ok=True)
    for stale_file in (manifest_path, subtype_manifest_path):
        stale_file.unlink(missing_ok=True)

    candidates = discover_candidate_isolates(args.discovery_limit, email=args.email, api_key=args.api_key)
    if not candidates:
        raise SystemExit("No rotavirus A subtype representatives discovered.")

    segment_rows: list[dict[str, str]] = []
    subtype_rows: list[dict[str, str]] = []
    combined_entries: list[str] = []

    for combo, isolate_candidates in sorted(candidates.items()):
        isolate_label = ""
        records: list[dict[str, str]] = []
        unique_segments: list[str] = []
        for candidate in isolate_candidates:
            candidate_records = fetch_complete_isolate_records(candidate, email=args.email, api_key=args.api_key)
            candidate_best = choose_best_segment_set(candidate_records)
            candidate_segments = sorted({row.get("segment_label") or "" for row in candidate_best if row.get("segment_label")})
            if len(candidate_segments) > len(unique_segments):
                isolate_label = candidate
                records = candidate_best
                unique_segments = candidate_segments
            if len(candidate_segments) >= 11:
                isolate_label = candidate
                records = candidate_best
                unique_segments = candidate_segments
                break
        if not records:
            print(f"[warn] {combo}: no complete records found among {len(isolate_candidates)} isolate candidates", file=sys.stderr)
            continue
        if len(unique_segments) < 11:
            print(
                f"[warn] {combo}: best isolate {isolate_label} only has {len(unique_segments)} unique segments, skip as incomplete",
                file=sys.stderr,
            )
            continue
        records.sort(key=lambda row: (row.get("segment_label") or "ZZZ", row.get("accession") or ""))
        safe_combo = combo.replace("[", "_").replace("]", "").replace("/", "_")
        safe_isolate = "".join(ch if ch.isalnum() or ch in {"-", "_", "."} else "_" for ch in isolate_label)[:120]
        per_subtype_fasta = out_dir / f"{safe_combo}__{safe_isolate}.fasta"
        subtype_entries: list[str] = []
        for row in records:
            subtype_entries.append(f">{row['header']}\n{_wrap_sequence(row['sequence'])}")
            segment_rows.append(
                {
                    "combo": combo,
                    "isolate_label": isolate_label,
                    "segment_label": row.get("segment_label", ""),
                    "accession": row["accession"],
                    "sequence_length": row["sequence_length"],
                    "header": row["header"],
                    "fasta_file": per_subtype_fasta.name,
                }
            )
        per_subtype_fasta.write_text("\n".join(subtype_entries).strip() + "\n", encoding="utf-8")
        combined_entries.extend(subtype_entries)
        subtype_rows.append(
            {
                "combo": combo,
                "isolate_label": isolate_label,
                "segment_count": str(len(records)),
                "segments": ",".join(unique_segments),
                "fasta_file": per_subtype_fasta.name,
            }
        )
        print(f"[done] {combo}: {len(records)} unique segments from {isolate_label}", file=sys.stderr)

    if not segment_rows:
        raise SystemExit("No complete subtype segment sets were downloaded.")

    combined_fasta_path.write_text("\n".join(combined_entries).strip() + "\n", encoding="utf-8")

    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["combo", "isolate_label", "segment_label", "accession", "sequence_length", "header", "fasta_file"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(segment_rows)

    with subtype_manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["combo", "isolate_label", "segment_count", "segments", "fasta_file"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(subtype_rows)

    print(f"[ok] wrote {combined_fasta_path}", file=sys.stderr)
    print(f"[ok] wrote {manifest_path}", file=sys.stderr)
    print(f"[ok] wrote {subtype_manifest_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
