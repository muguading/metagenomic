#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

from download_rega_hiv_reference_genomes import COMBINED_FASTA_PATH, FASTA_DIR, VALID_HEADER_MARKERS, fetch_fasta
from download_rega_hiv_reference_manifest import MANIFEST_PATH, OUTPUT_DIR


SUPPLEMENT_DIR = OUTPUT_DIR / "supplement_fasta"
SUPPLEMENT_MANIFEST_PATH = OUTPUT_DIR / "supplement_manifest.tsv"
SUPPLEMENT_SUMMARY_PATH = OUTPUT_DIR / "supplement_summary.json"
ESEARCH_URL = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    "?db=nucleotide&retmode=json&retmax=50&term={term}"
)
ESUMMARY_URL = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    "?db=nucleotide&retmode=json&id={ids}"
)
CURATED_SUPPLEMENTS = {
    "Subtype A2": ["ON902752"],
    "CRF04_CPX": ["AF049337"],
    "CRF05_DF": ["AF193253", "AF076998"],
    "CRF10_CD": ["AF289548", "AF289549", "AF289550"],
}


def ncbi_json(url: str) -> dict:
    with urllib.request.urlopen(url, timeout=60) as response:
        return json.loads(response.read().decode("utf-8"))


def normalize_token(value: str) -> str:
    return "".join(ch for ch in value.upper() if ch.isalnum())


def group_query(group: str) -> tuple[str, str]:
    if group.startswith("Subtype "):
        code = group.split()[-1]
        term = f'(Human immunodeficiency virus 1[Organism]) AND ("subtype {code}"[All Fields]) AND (complete genome[All Fields])'
        return code, term
    code = group
    spaced = code.replace("_", " ")
    term = (
        f'(Human immunodeficiency virus 1[Organism]) AND '
        f'("{code}"[All Fields] OR "{spaced}"[All Fields]) AND '
        f'(complete genome[All Fields])'
    )
    return code, term


def subtype_tag(summary: dict) -> str:
    subname = summary.get("subname", "")
    if not subname:
        return ""
    return subname.split("|", 1)[0].strip()


def header_title(fasta: str) -> str:
    return fasta.splitlines()[0].lstrip(">").strip()


def is_valid_candidate(summary: dict, desired_code: str) -> bool:
    if summary.get("organism") != "Human immunodeficiency virus 1":
        return False
    if summary.get("completeness", "").lower() != "complete":
        return False
    title = summary.get("title", "").lower()
    if "complete genome" not in title:
        return False
    if not any(marker in title or marker in summary.get("organism", "").lower() for marker in VALID_HEADER_MARKERS):
        return False

    desired = normalize_token(desired_code)
    subtype = subtype_tag(summary)
    subtype_norm = normalize_token(subtype)
    if not subtype_norm:
        return False
    if subtype_norm != desired:
        return False
    return True


def existing_accessions() -> set[str]:
    accessions: set[str] = set()
    with MANIFEST_PATH.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["download_status"] == "downloaded" and row["note"]:
                accessions.add(row["note"].split()[0].split(";")[0])
    return accessions


def load_missing_groups() -> dict[str, int]:
    missing: dict[str, int] = {}
    with MANIFEST_PATH.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["download_status"] != "downloaded":
                missing[row["group"]] = missing.get(row["group"], 0) + 1
    return missing


def candidate_summaries(group: str) -> list[dict]:
    desired_code, term = group_query(group)
    url = ESEARCH_URL.format(term=urllib.parse.quote(term, safe=""))
    ids = ncbi_json(url)["esearchresult"]["idlist"]
    if not ids:
        return []
    time.sleep(0.34)
    summary_url = ESUMMARY_URL.format(ids=",".join(ids))
    payload = ncbi_json(summary_url)["result"]
    summaries = [payload[uid] for uid in payload["uids"]]
    valid = [summary for summary in summaries if is_valid_candidate(summary, desired_code)]
    return valid


def rebuild_combined_fasta() -> None:
    chunks: list[str] = []
    for path in sorted(FASTA_DIR.glob("*.fasta")):
        chunks.append(path.read_text(encoding="utf-8"))
    for path in sorted(SUPPLEMENT_DIR.glob("*.fasta")):
        chunks.append(path.read_text(encoding="utf-8"))
    COMBINED_FASTA_PATH.write_text("".join(chunks), encoding="utf-8")


def main() -> int:
    SUPPLEMENT_DIR.mkdir(parents=True, exist_ok=True)
    for stale in SUPPLEMENT_DIR.glob("*.fasta"):
        stale.unlink()

    missing = load_missing_groups()
    seen = existing_accessions()
    rows: list[dict[str, str]] = []
    summary = {"downloaded": [], "insufficient": []}

    for group in sorted(missing):
        need = missing[group]
        chosen = 0
        for candidate in candidate_summaries(group):
            accession = candidate["accessionversion"].split(".", 1)[0]
            if accession in seen:
                continue
            fasta = fetch_fasta(accession)
            path = SUPPLEMENT_DIR / f"{group.replace('/', '_')}__supplement__{accession}.fasta"
            path.write_text(fasta, encoding="utf-8")
            seen.add(accession)
            rows.append(
                {
                    "group": group,
                    "accession": accession,
                    "accession_version": candidate["accessionversion"],
                    "subtype": subtype_tag(candidate),
                    "title": candidate["title"],
                    "path": str(path),
                }
            )
            summary["downloaded"].append(
                {
                    "group": group,
                    "accession": accession,
                    "title": candidate["title"],
                }
            )
            chosen += 1
            time.sleep(0.34)
            if chosen >= need:
                break

        for accession in CURATED_SUPPLEMENTS.get(group, []):
            if chosen >= need:
                break
            if accession in seen:
                continue
            fasta = fetch_fasta(accession)
            path = SUPPLEMENT_DIR / f"{group.replace('/', '_')}__supplement__{accession}.fasta"
            path.write_text(fasta, encoding="utf-8")
            seen.add(accession)
            rows.append(
                {
                    "group": group,
                    "accession": accession,
                    "accession_version": accession,
                    "subtype": group.replace("Subtype ", ""),
                    "title": header_title(fasta),
                    "path": str(path),
                }
            )
            summary["downloaded"].append(
                {
                    "group": group,
                    "accession": accession,
                    "title": header_title(fasta),
                }
            )
            chosen += 1
            time.sleep(0.34)
        if chosen < need:
            summary["insufficient"].append({"group": group, "needed": need, "downloaded": chosen})

    with SUPPLEMENT_MANIFEST_PATH.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["group", "accession", "accession_version", "subtype", "title", "path"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    SUPPLEMENT_SUMMARY_PATH.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    rebuild_combined_fasta()
    print(
        json.dumps(
            {
                "supplement_manifest": str(SUPPLEMENT_MANIFEST_PATH),
                "supplement_summary": str(SUPPLEMENT_SUMMARY_PATH),
                "downloaded": len(summary["downloaded"]),
                "insufficient": len(summary["insufficient"]),
            },
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
