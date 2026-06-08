from __future__ import annotations

import argparse
import csv
import http.client
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path


EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
PHF_DB_PATHS = {
    "fiber": Path("database/virus/hadv/blastn_db_fiber/hadv_types_ref_fiber.fa"),
    "hexon": Path("database/virus/hadv/blastn_db_hexon/hadv_types_ref_hexon.fa"),
    "penton": Path("database/virus/hadv/blastn_db_penton/hadv_types_ref_penton.fa"),
}


def _request(url: str, retries: int = 5, sleep_seconds: float = 1.5) -> str:
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        request = urllib.request.Request(
            url,
            headers={"User-Agent": "metagenomic-hadv-reference-downloader/1.0"},
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


def _parse_reference_accessions() -> dict[str, dict[str, str]]:
    refs: dict[str, dict[str, str]] = {}
    for gene_name, path in PHF_DB_PATHS.items():
        with path.resolve().open(encoding="utf-8") as handle:
            for line in handle:
                if not line.startswith(">"):
                    continue
                header = line[1:].strip()
                seq_id = header.split()[0]
                accession = seq_id.split("_HAdV-")[0].split(".")[0]
                matched = re.search(r"(HAdV-[A-Z]\d+)", header, flags=re.IGNORECASE)
                type_label = matched.group(1).upper() if matched else ""
                row = refs.setdefault(
                    accession,
                    {
                        "accession": accession,
                        "type_label": type_label,
                        "fiber_header": "",
                        "hexon_header": "",
                        "penton_header": "",
                    },
                )
                row[f"{gene_name}_header"] = header
                if type_label and not row.get("type_label"):
                    row["type_label"] = type_label
    return refs


def _efetch_fasta(accession: str, email: str = "", api_key: str = "") -> str:
    params = {
        "db": "nuccore",
        "id": accession,
        "rettype": "fasta",
        "retmode": "text",
    }
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key
    url = f"{EUTILS_BASE}/efetch.fcgi?{urllib.parse.urlencode(params)}"
    return _request(url)


def _wrap_sequence(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[index:index + width] for index in range(0, len(sequence), width))


def main() -> int:
    parser = argparse.ArgumentParser(description="Download representative full genomes for the HAdV PHF typing database")
    parser.add_argument(
        "--existing-fasta",
        type=Path,
        default=Path("database/virus/hadv/full_genomes/human_mastadenovirus_A_G_complete_genomes.fasta"),
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("database/virus/hadv/reference_genomes"),
    )
    parser.add_argument("--email", default="", help="Optional email for NCBI E-utilities")
    parser.add_argument("--api-key", default="", help="Optional NCBI API key")
    args = parser.parse_args()

    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    ref_rows = _parse_reference_accessions()
    existing_records = _read_fasta(args.existing_fasta.resolve())
    existing_by_accession = {seq_id.split(".")[0]: (header, sequence) for seq_id, (header, sequence) in existing_records.items()}

    fasta_entries: list[str] = []
    manifest_rows: list[dict[str, str]] = []
    missing_accessions: list[str] = []

    for accession in sorted(ref_rows):
        header = ""
        sequence = ""
        source = "downloaded_full_genomes"
        if accession in existing_by_accession:
            header, sequence = existing_by_accession[accession]
        else:
            print(f"[fetch] {accession}", file=sys.stderr)
            fasta_text = _efetch_fasta(accession, email=args.email, api_key=args.api_key).strip()
            time.sleep(0.34)
            if not fasta_text.startswith(">"):
                missing_accessions.append(accession)
                continue
            lines = fasta_text.splitlines()
            header = lines[0][1:].strip()
            sequence = "".join(line.strip() for line in lines[1:])
            source = "ncbi_efetch"
        fasta_entries.append(f">{header}\n{_wrap_sequence(sequence)}\n")
        manifest_rows.append(
            {
                "accession": accession,
                "type_label": ref_rows[accession].get("type_label", ""),
                "sequence_length": str(len(sequence)),
                "source": source,
                "fiber_header": ref_rows[accession].get("fiber_header", ""),
                "hexon_header": ref_rows[accession].get("hexon_header", ""),
                "penton_header": ref_rows[accession].get("penton_header", ""),
                "full_genome_header": header,
            }
        )

    out_fasta = out_dir / "hadv_phf_typing_db_representative_genomes.fasta"
    out_manifest = out_dir / "hadv_phf_typing_db_representative_genomes.tsv"
    out_fasta.write_text("".join(fasta_entries), encoding="utf-8")
    with out_manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "accession",
                "type_label",
                "sequence_length",
                "source",
                "fiber_header",
                "hexon_header",
                "penton_header",
                "full_genome_header",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in manifest_rows:
            writer.writerow(row)

    if missing_accessions:
        print("[warn] missing accessions: " + ",".join(missing_accessions), file=sys.stderr)

    print(f"[ok] wrote {out_fasta}", file=sys.stderr)
    print(f"[ok] wrote {out_manifest}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
