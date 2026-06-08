#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations

import csv
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path


ROOT = Path("/Users/wuhhh/Desktop/徐老师/代码/metagenomic")
MISSING_TSV = ROOT / "database" / "virus" / "enterovirus" / "reference_genomes" / "abcd_vp1" / "enterovirus_abcd_subtype_completeness.tsv"
REF_DIR = ROOT / "database" / "virus" / "enterovirus" / "reference_genomes"
GFF_DIR = REF_DIR / "gff3"
REPORT_TSV = ROOT / "database" / "virus" / "enterovirus" / "reference_genomes" / "abcd_vp1" / "download_missing_complete_genomes_report.tsv"

ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
GFF_URL = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi"
REQUEST_SLEEP = 0.34


def fetch_text(url: str, params: dict[str, str], timeout: int = 60) -> str:
    query = urllib.parse.urlencode(params)
    full_url = f"{url}?{query}"
    request = urllib.request.Request(
        full_url,
        headers={"User-Agent": "metagenomic-enterovirus-abcd-downloader/1.0"},
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read().decode("utf-8", errors="replace")


def fetch_json(url: str, params: dict[str, str], timeout: int = 60) -> dict:
    return json.loads(fetch_text(url, params, timeout=timeout))


def load_missing_targets() -> list[dict[str, str]]:
    targets = []
    with MISSING_TSV.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["status"] == "no_complete_genome":
                targets.append(row)
    return targets


def build_queries(row: dict[str, str]) -> list[str]:
    virus_name = row["virus_name"].strip()
    abbrev = row["abbrev"].strip()
    queries = []
    if virus_name:
        queries.append(f'"{virus_name}"[Title] AND ("complete genome"[Title] OR "complete sequence"[Title])')
        queries.append(f'"{virus_name}"[All Fields] AND ("complete genome"[Title] OR "complete sequence"[Title])')
    if abbrev.startswith("RV-"):
        suffix = abbrev.replace("RV-", "")
        queries.append(f'"rhinovirus {suffix}"[Title] AND ("complete genome"[Title] OR "complete sequence"[Title])')
        queries.append(f'"human rhinovirus {suffix}"[Title] AND ("complete genome"[Title] OR "complete sequence"[Title])')
    elif abbrev.startswith("EV-"):
        suffix = abbrev.replace("EV-", "")
        queries.append(f'"enterovirus {suffix}"[Title] AND ("complete genome"[Title] OR "complete sequence"[Title])')
    elif abbrev.startswith("CVA"):
        suffix = abbrev.replace("CVA", "")
        queries.append(f'"coxsackievirus A{suffix}"[Title] AND ("complete genome"[Title] OR "complete sequence"[Title])')
    elif abbrev.startswith("CVB"):
        suffix = abbrev.replace("CVB", "")
        queries.append(f'"coxsackievirus B{suffix}"[Title] AND ("complete genome"[Title] OR "complete sequence"[Title])')
    elif abbrev.startswith("E"):
        queries.append(f'"echovirus {abbrev[1:]}"[Title] AND ("complete genome"[Title] OR "complete sequence"[Title])')
    return list(dict.fromkeys(query for query in queries if query))


def esearch_ids(term: str, retmax: int = 10) -> list[str]:
    payload = fetch_json(
        ESEARCH_URL,
        {
            "db": "nuccore",
            "retmode": "json",
            "term": term,
            "retmax": str(retmax),
        },
    )
    time.sleep(REQUEST_SLEEP)
    return payload.get("esearchresult", {}).get("idlist", [])


def esummary_one(uid: str) -> dict[str, str]:
    payload = fetch_json(
        ESUMMARY_URL,
        {"db": "nuccore", "id": uid, "retmode": "json"},
    )
    time.sleep(REQUEST_SLEEP)
    result = payload.get("result", {})
    return result.get(uid, {})


def looks_complete(title: str, length: int) -> bool:
    text = (title or "").lower()
    if "nearly complete" in text or "near complete" in text:
        return False
    if "complete genome" in text:
        return True
    if "complete sequence" in text and length >= 6000:
        return True
    return False


def download_fasta(accession: str, out_path: Path) -> None:
    fasta_text = fetch_text(
        EFETCH_URL,
        {"db": "nuccore", "id": accession, "rettype": "fasta", "retmode": "text"},
    )
    if not fasta_text.lstrip().startswith(">"):
        raise RuntimeError(f"{accession}: invalid FASTA response")
    out_path.write_text(fasta_text, encoding="utf-8")
    time.sleep(REQUEST_SLEEP)


def download_gff(accession: str, out_path: Path) -> None:
    gff_text = fetch_text(
        GFF_URL,
        {"id": accession, "db": "nuccore", "report": "gff3", "retmode": "text"},
    )
    if "##gff-version" not in gff_text:
        raise RuntimeError(f"{accession}: invalid GFF3 response")
    out_path.write_text(gff_text, encoding="utf-8")
    time.sleep(REQUEST_SLEEP)


def choose_candidate(row: dict[str, str]) -> dict[str, str] | None:
    seen_uids = set()
    for term in build_queries(row):
        for uid in esearch_ids(term):
            if uid in seen_uids:
                continue
            seen_uids.add(uid)
            summary = esummary_one(uid)
            accession = str(summary.get("caption") or "").strip()
            title = str(summary.get("title") or "").strip()
            length = int(summary.get("slen") or summary.get("length") or 0)
            if not accession:
                continue
            if not looks_complete(title, length):
                continue
            return {
                "uid": uid,
                "accession": accession,
                "title": title,
                "length": str(length),
                "term": term,
            }
    return None


def main() -> int:
    targets = load_missing_targets()
    rows = []
    for row in targets:
        result = {
            "big_group": row["big_group"],
            "abbrev": row["abbrev"],
            "virus_name": row["virus_name"],
            "species": row["species"],
            "status": "",
            "download_accession": "",
            "download_title": "",
            "download_length": "",
            "query_term": "",
            "fasta_path": "",
            "gff_path": "",
            "note": "",
        }
        try:
            candidate = choose_candidate(row)
            if candidate is None:
                result["status"] = "not_found"
                result["note"] = "NCBI 未检到符合 complete genome/complete sequence 条件的候选"
                rows.append(result)
                continue
            accession = candidate["accession"]
            fasta_path = REF_DIR / f"{accession}.fasta"
            gff_path = GFF_DIR / f"{accession}.gff3"
            if not fasta_path.exists() or fasta_path.stat().st_size == 0:
                download_fasta(accession, fasta_path)
            if not gff_path.exists() or gff_path.stat().st_size == 0:
                download_gff(accession, gff_path)
            result.update(
                {
                    "status": "downloaded",
                    "download_accession": accession,
                    "download_title": candidate["title"],
                    "download_length": candidate["length"],
                    "query_term": candidate["term"],
                    "fasta_path": str(fasta_path),
                    "gff_path": str(gff_path),
                    "note": "",
                }
            )
        except Exception as exc:  # noqa: BLE001
            result["status"] = "error"
            result["note"] = str(exc)
        rows.append(result)
        print(f"{result['status']}\t{result['abbrev']}\t{result['download_accession'] or '-'}")

    with REPORT_TSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "big_group",
                "abbrev",
                "virus_name",
                "species",
                "status",
                "download_accession",
                "download_title",
                "download_length",
                "query_term",
                "fasta_path",
                "gff_path",
                "note",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    downloaded = sum(1 for row in rows if row["status"] == "downloaded")
    print(f"Downloaded {downloaded} complete genomes. Report: {REPORT_TSV}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
