#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import time
import urllib.error
import urllib.request
from pathlib import Path

from download_rega_hiv_reference_manifest import MANIFEST_PATH, OUTPUT_DIR


FASTA_DIR = OUTPUT_DIR / "fasta"
SUMMARY_PATH = OUTPUT_DIR / "download_summary.json"
COMBINED_FASTA_PATH = OUTPUT_DIR / "rega_hiv_reference_genomes.fasta"
EFETCH_URL = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    "?db=nucleotide&id={accession}&rettype=fasta&retmode=text"
)
VALID_HEADER_MARKERS = (
    "human immunodeficiency virus",
    "hiv-1",
    "hiv1",
    "hiv type 1",
)


def fetch_fasta(accession: str) -> str:
    url = EFETCH_URL.format(accession=accession)
    with urllib.request.urlopen(url, timeout=30) as response:
        payload = response.read().decode("utf-8", errors="replace").strip()
    if not payload.startswith(">"):
        raise ValueError(f"{accession} did not return FASTA")
    header = payload.splitlines()[0].lower()
    if not any(marker in header for marker in VALID_HEADER_MARKERS):
        raise ValueError(f"{accession} returned a non-HIV sequence")
    return payload + "\n"


def load_manifest() -> list[dict[str, str]]:
    with MANIFEST_PATH.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_manifest(rows: list[dict[str, str]]) -> None:
    fieldnames = ["group", "label", "accession_candidates", "download_status", "note"]
    with MANIFEST_PATH.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    FASTA_DIR.mkdir(parents=True, exist_ok=True)
    for stale_fasta in FASTA_DIR.glob("*.fasta"):
        stale_fasta.unlink()
    rows = load_manifest()
    summary = {
        "downloaded": [],
        "unresolved": [],
        "failed": [],
    }
    combined_chunks: list[str] = []

    for row in rows:
        candidates = [item for item in row["accession_candidates"].split(",") if item]
        if not candidates:
            row["download_status"] = "unresolved"
            summary["unresolved"].append(
                {
                    "group": row["group"],
                    "label": row["label"],
                    "reason": row["note"] or "No accession candidate",
                }
            )
            continue

        fasta_text = ""
        chosen_accession = ""
        failures: list[str] = []
        for accession in candidates:
            try:
                fasta_text = fetch_fasta(accession)
                chosen_accession = accession
                break
            except (urllib.error.URLError, ValueError) as exc:
                failures.append(f"{accession}: {exc}")
                time.sleep(0.34)

        if not fasta_text:
            row["download_status"] = "failed"
            row["note"] = "; ".join(failures)
            summary["failed"].append(
                {
                    "group": row["group"],
                    "label": row["label"],
                    "candidates": candidates,
                    "errors": failures,
                }
            )
            continue

        fasta_path = FASTA_DIR / f"{row['group'].replace('/', '_')}__{row['label']}__{chosen_accession}.fasta"
        fasta_path.write_text(fasta_text, encoding="utf-8")
        combined_chunks.append(fasta_text)
        row["download_status"] = "downloaded"
        row["note"] = chosen_accession
        summary["downloaded"].append(
            {
                "group": row["group"],
                "label": row["label"],
                "accession": chosen_accession,
                "path": str(fasta_path),
            }
        )
        time.sleep(0.34)

    COMBINED_FASTA_PATH.write_text("".join(combined_chunks), encoding="utf-8")
    SUMMARY_PATH.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    write_manifest(rows)
    print(json.dumps(
        {
            "manifest": str(MANIFEST_PATH),
            "combined_fasta": str(COMBINED_FASTA_PATH),
            "summary": str(SUMMARY_PATH),
            "downloaded": len(summary["downloaded"]),
            "failed": len(summary["failed"]),
            "unresolved": len(summary["unresolved"]),
        },
        ensure_ascii=False,
    ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
