from __future__ import annotations

import argparse
import csv
import http.client
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path


EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def _request(url: str, retries: int = 5, sleep_seconds: float = 1.5) -> str:
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        request = urllib.request.Request(
            url,
            headers={"User-Agent": "metagenomic-rhinovirus-gff-downloader/1.0"},
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


def _fetch_gff3(accession: str, email: str = "", api_key: str = "") -> str:
    params = {
        "db": "nuccore",
        "id": accession,
        "rettype": "gff3",
        "retmode": "text",
    }
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key
    url = f"{EUTILS_BASE}/efetch.fcgi?{urllib.parse.urlencode(params)}"
    return _request(url)


def main() -> int:
    parser = argparse.ArgumentParser(description="Download GFF3 annotations for rhinovirus representative genomes")
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("database/virus/rhinovirus/full_genomes/human_rhinovirus_A_B_C_dedup_by_type_manifest.tsv"),
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("database/virus/rhinovirus/reference_genomes/gff3"),
    )
    parser.add_argument("--email", default="", help="Optional email for NCBI E-utilities")
    parser.add_argument("--api-key", default="", help="Optional NCBI API key")
    args = parser.parse_args()

    manifest_path = args.manifest.resolve()
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    with manifest_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    status_rows: list[dict[str, str]] = []
    for row in rows:
        accession = str(row.get("accession") or "").strip()
        species_group = str(row.get("species_group") or "").strip()
        type_label = str(row.get("type_label") or "").strip()
        header = str(row.get("header") or "").strip()
        out_path = out_dir / f"{accession.split('.')[0]}.gff3"

        if out_path.is_file() and out_path.stat().st_size > 0:
            status_rows.append(
                {
                    "accession": accession,
                    "species_group": species_group,
                    "type_label": type_label,
                    "status": "existing",
                    "gff_path": str(out_path),
                    "header": header,
                }
            )
            continue

        print(f"[fetch] HRV-{species_group}{type_label} <- {accession}", file=sys.stderr)
        text = _fetch_gff3(accession, email=args.email, api_key=args.api_key).strip()
        time.sleep(0.34)
        if text.startswith("##gff-version 3"):
            out_path.write_text(text + "\n", encoding="utf-8")
            status = "downloaded"
        else:
            status = "missing"
        status_rows.append(
            {
                "accession": accession,
                "species_group": species_group,
                "type_label": type_label,
                "status": status,
                "gff_path": str(out_path) if status != "missing" else "",
                "header": header,
            }
        )

    summary_path = out_dir / "download_status.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["accession", "species_group", "type_label", "status", "gff_path", "header"],
            delimiter="\t",
        )
        writer.writeheader()
        for row in status_rows:
            writer.writerow(row)

    print(f"[ok] wrote {summary_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
