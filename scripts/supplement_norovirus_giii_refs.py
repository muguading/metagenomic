from __future__ import annotations

import csv
import textwrap
import urllib.parse
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REF_DIR = ROOT / "database/virus/norovirus/cdc_typing_refs"


SUPPLEMENTS = {
    "rdrp": [
        {
            "subtype": "GIII.P1",
            "accession": "AJ011099",
            "fetch_id": "AJ011099.1",
            "label": "Jena_AJ011099",
        },
        {
            "subtype": "GIII.P2",
            "accession": "AF097917",
            "fetch_id": "AF097917.5",
            "label": "Newbury2_AF097917",
        },
        {
            "subtype": "GIII.P2",
            "accession": "MK159169",
            "fetch_id": "MK159169.1",
            "label": "Bo-BET-17-18-CH_MK159169",
        },
    ],
    "vp1": [
        {
            "subtype": "GIII.1",
            "accession": "AJ011099",
            "fetch_id": "AJ011099.1",
            "label": "Jena_AJ011099",
        },
        {
            "subtype": "GIII.2",
            "accession": "AF097917",
            "fetch_id": "AF097917.5",
            "label": "Newbury2_AF097917",
        },
        {
            "subtype": "GIII.3",
            "accession": "EU193658",
            "fetch_id": "EU193658.4",
            "label": "Norsewood30_EU193658",
        },
        {
            "subtype": "GIII.4",
            "accession": "MK159169",
            "fetch_id": "MK159169.1",
            "label": "Bo-BET-17-18-CH_MK159169",
        },
    ],
}


def fetch_fasta(fetch_id: str) -> tuple[str, str]:
    params = urllib.parse.urlencode(
        {
            "db": "nuccore",
            "id": fetch_id,
            "rettype": "fasta",
            "retmode": "text",
        }
    )
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?{params}"
    with urllib.request.urlopen(url, timeout=60) as response:
        text = response.read().decode("utf-8")
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines or not lines[0].startswith(">"):
        raise RuntimeError(f"Failed to fetch FASTA for {fetch_id}")
    header = lines[0][1:].strip()
    seq = "".join(lines[1:]).upper()
    if not seq:
        raise RuntimeError(f"Empty sequence for {fetch_id}")
    return header, seq


def wrap_fasta(seq: str, width: int = 80) -> str:
    return "\n".join(textwrap.wrap(seq, width))


def append_manifest_rows(tsv_path: Path, gene_name: str, rows: list[dict[str, str]]) -> int:
    existing_keys: set[tuple[str, str]] = set()
    with tsv_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            existing_keys.add((str(row.get("accession") or ""), str(row.get("subtype") or "")))

    added = 0
    with tsv_path.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["source", "gene", "subtype", "accession", "label", "source_url"],
            delimiter="\t",
        )
        for row in rows:
            key = (row["accession"], row["subtype"])
            if key in existing_keys:
                continue
            writer.writerow(
                {
                    "source": "NCBI supplement for GIII",
                    "gene": gene_name,
                    "subtype": row["subtype"],
                    "accession": row["accession"],
                    "label": row["label"],
                    "source_url": f"https://www.ncbi.nlm.nih.gov/nuccore/{row['fetch_id']}",
                }
            )
            existing_keys.add(key)
            added += 1
    return added


def append_fasta_entries(fasta_path: Path, gene_key: str, rows: list[dict[str, str]]) -> int:
    existing_headers: set[str] = set()
    with fasta_path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                existing_headers.add(line[1:].split()[0])

    entries: list[str] = []
    added = 0
    gene_suffix = gene_key.upper()
    for row in rows:
        header, seq = fetch_fasta(row["fetch_id"])
        record_id = f"{row['accession']}_{row['subtype']}_{gene_suffix}"
        if record_id in existing_headers:
            continue
        entries.append(f">{record_id} {header}\n{wrap_fasta(seq)}\n")
        existing_headers.add(record_id)
        added += 1

    if entries:
        with fasta_path.open("a", encoding="utf-8") as handle:
            if fasta_path.stat().st_size > 0:
                handle.write("\n")
            handle.write("\n".join(entry.rstrip("\n") for entry in entries))
            handle.write("\n")
    return added


def main() -> int:
    summary: list[str] = []
    for gene_key, rows in SUPPLEMENTS.items():
        gene_name = "RdRp" if gene_key == "rdrp" else "VP1"
        tsv_path = REF_DIR / f"cdc_norovirus_{gene_key}_refs.tsv"
        fasta_path = REF_DIR / f"cdc_norovirus_{gene_key}_refs.fasta"
        added_manifest = append_manifest_rows(tsv_path, gene_name, rows)
        added_fasta = append_fasta_entries(fasta_path, gene_key, rows)
        summary.append(
            f"{gene_key}\tmanifest_added={added_manifest}\tfasta_added={added_fasta}"
        )
    print("\n".join(summary))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
