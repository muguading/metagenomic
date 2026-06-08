#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path


ROOT = Path("/Users/wuhhh/Desktop/徐老师/代码/metagenomic")
DEFAULT_QUERY = ROOT / "database" / "virus" / "astroviridae" / "astro_db.fasta"
DEFAULT_REF_FASTA = ROOT / "database" / "virus" / "astroviridae" / "reference_genomes" / "astroviridae_orf2_capsid_references.fasta"
DEFAULT_REF_MANIFEST = ROOT / "database" / "virus" / "astroviridae" / "reference_genomes" / "astroviridae_orf2_capsid_references.tsv"
DEFAULT_OUT_DIR = ROOT / "tmp" / "astroviridae_orf2_typing_demo"
DEFAULT_BLASTN = Path("/opt/homebrew/Caskroom/mambaforge/base/envs/ncov/bin/blastn")

OUTFMT_FIELDS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "slen",
]


def resolve_blastn(explicit: str) -> str:
    candidate = str(explicit or "").strip()
    if candidate:
        return candidate
    if DEFAULT_BLASTN.is_file():
        return str(DEFAULT_BLASTN)
    return "blastn"


def iter_fasta(path: Path):
    header = ""
    seq_lines: list[str] = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header:
                    yield header, "".join(seq_lines)
                header = line[1:].strip()
                seq_lines = []
                continue
            if header:
                seq_lines.append(line.strip())
    if header:
        yield header, "".join(seq_lines)


def read_query_records(path: Path, limit: int) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    for index, (header, sequence) in enumerate(iter_fasta(path), start=1):
        query_id = header.split()[0]
        records.append({"id": query_id, "header": header, "sequence": sequence})
        if limit > 0 and index >= limit:
            break
    return records


def write_subset_fasta(records: list[dict[str, str]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as handle:
        for record in records:
            sequence = record["sequence"]
            wrapped = "\n".join(sequence[index:index + 80] for index in range(0, len(sequence), 80))
            handle.write(f">{record['header']}\n{wrapped}\n")


def load_manifest(path: Path) -> dict[str, dict[str, str]]:
    mapping: dict[str, dict[str, str]] = {}
    with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            accession = str(row.get("accession") or "").strip()
            header_token = str(row.get("header") or "").strip().split()[0].split(".", 1)[0]
            for key in {accession, header_token}:
                if key:
                    mapping[key] = row
    return mapping


def run_blastn(blastn_bin: str, query_fasta: Path, ref_fasta: Path, out_path: Path, threads: int) -> None:
    cmd = [
        blastn_bin,
        "-task",
        "blastn",
        "-query",
        str(query_fasta),
        "-subject",
        str(ref_fasta),
        "-outfmt",
        "6 " + " ".join(OUTFMT_FIELDS),
        "-evalue",
        "1e-20",
        "-max_target_seqs",
        "50",
        "-dust",
        "no",
        "-num_threads",
        str(max(1, int(threads or 1))),
        "-out",
        str(out_path),
    ]
    subprocess.run(cmd, check=True)


def normalize_subject_id(subject_id: str) -> str:
    return str(subject_id or "").strip().split(".", 1)[0]


def parse_blast_hits(path: Path, manifest: dict[str, dict[str, str]]) -> tuple[list[dict[str, str]], dict[str, dict[str, str]]]:
    all_rows: list[dict[str, str]] = []
    best_by_query: dict[str, dict[str, str]] = {}
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != len(OUTFMT_FIELDS):
                continue
            row = dict(zip(OUTFMT_FIELDS, parts))
            subject_key = normalize_subject_id(row["sseqid"])
            meta = manifest.get(subject_key, {})
            try:
                align_len = int(float(row["length"]))
                subject_len = int(float(row["slen"]))
                pident = float(row["pident"])
                bitscore = float(row["bitscore"])
            except (TypeError, ValueError):
                continue
            subject_cov = min(100.0, (align_len / subject_len) * 100.0) if subject_len else 0.0
            result = {
                "query_id": row["qseqid"],
                "subject_id": row["sseqid"],
                "subject_accession": subject_key,
                "genus": str(meta.get("genus") or "").strip(),
                "species": str(meta.get("species") or "").strip(),
                "subtype": str(meta.get("abbrev") or "").strip(),
                "virus_name": str(meta.get("virus_name") or "").strip(),
                "reference_isolate": str(meta.get("isolate") or "").strip(),
                "pident": f"{pident:.2f}",
                "align_length": str(align_len),
                "subject_coverage_pct": f"{subject_cov:.2f}",
                "bitscore": f"{bitscore:.1f}",
                "evalue": str(row["evalue"]),
                "reference_header": str(meta.get("header") or ""),
            }
            all_rows.append(result)
            current_score = (bitscore, subject_cov, pident, align_len)
            previous = best_by_query.get(result["query_id"])
            if previous is None or current_score > previous.get("_score", (0.0, 0.0, 0.0, 0)):
                result["_score"] = current_score
                best_by_query[result["query_id"]] = result
    return all_rows, best_by_query


def build_summary_rows(records: list[dict[str, str]], best_by_query: dict[str, dict[str, str]]) -> list[dict[str, str]]:
    summary_rows: list[dict[str, str]] = []
    for record in records:
        query_id = record["id"]
        best = best_by_query.get(query_id)
        sequence_length = len(record["sequence"])
        if best is None:
            summary_rows.append(
                {
                    "query_id": query_id,
                    "sequence_length": str(sequence_length),
                    "genus": "",
                    "species": "",
                    "subtype": "",
                    "virus_name": "",
                    "subject_accession": "",
                    "subject_coverage_pct": "0.00",
                    "pident": "0.00",
                    "bitscore": "0.0",
                    "status": "unassigned",
                    "note": "未命中 ORF2 参考库",
                }
            )
            continue
        subject_cov = float(best["subject_coverage_pct"])
        pident = float(best["pident"])
        if subject_cov >= 80.0 and pident >= 80.0:
            status = "typed"
            note = "ORF2 高置信命中"
        elif subject_cov >= 60.0 and pident >= 70.0:
            status = "review"
            note = "ORF2 中等命中，建议复核"
        else:
            status = "weak"
            note = "ORF2 弱命中，建议复核"
        summary_rows.append(
            {
                "query_id": query_id,
                "sequence_length": str(sequence_length),
                "genus": best["genus"],
                "species": best["species"],
                "subtype": best["subtype"],
                "virus_name": best["virus_name"],
                "subject_accession": best["subject_accession"],
                "subject_coverage_pct": best["subject_coverage_pct"],
                "pident": best["pident"],
                "bitscore": best["bitscore"],
                "status": status,
                "note": note,
            }
        )
    return summary_rows


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: ("" if key.startswith("_") else row.get(key, "")) for key in fieldnames})


def main() -> int:
    parser = argparse.ArgumentParser(description="Type astroviridae sequences by local ORF2 capsid reference blastn")
    parser.add_argument("--input-fasta", "--input", dest="input_fasta", type=Path, default=DEFAULT_QUERY)
    parser.add_argument("--reference-fasta", type=Path, default=DEFAULT_REF_FASTA)
    parser.add_argument("--reference-manifest", type=Path, default=DEFAULT_REF_MANIFEST)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--blastn-bin", type=str, default="")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--limit", type=int, default=0, help="Only type the first N sequences; use 0 for all")
    args = parser.parse_args()

    input_fasta = args.input_fasta.resolve()
    ref_fasta = args.reference_fasta.resolve()
    ref_manifest = args.reference_manifest.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if not input_fasta.is_file():
        raise FileNotFoundError(f"Input FASTA not found: {input_fasta}")
    if not ref_fasta.is_file():
        raise FileNotFoundError(f"Reference FASTA not found: {ref_fasta}")
    if not ref_manifest.is_file():
        raise FileNotFoundError(f"Reference manifest not found: {ref_manifest}")

    records = read_query_records(input_fasta, args.limit)
    if not records:
        raise RuntimeError(f"No FASTA records found in {input_fasta}")

    manifest = load_manifest(ref_manifest)
    subset_fasta = output_dir / "query_subset.fasta"
    blast_out = output_dir / "orf2_blastn.tsv"
    summary_tsv = output_dir / "typing_summary.tsv"
    all_hits_tsv = output_dir / "typing_all_hits.tsv"

    write_subset_fasta(records, subset_fasta)
    run_blastn(resolve_blastn(args.blastn_bin), subset_fasta, ref_fasta, blast_out, args.threads)
    all_rows, best_by_query = parse_blast_hits(blast_out, manifest)
    summary_rows = build_summary_rows(records, best_by_query)

    write_tsv(
        summary_tsv,
        summary_rows,
        [
            "query_id",
            "sequence_length",
            "genus",
            "species",
            "subtype",
            "virus_name",
            "subject_accession",
            "subject_coverage_pct",
            "pident",
            "bitscore",
            "status",
            "note",
        ],
    )
    write_tsv(
        all_hits_tsv,
        all_rows,
        [
            "query_id",
            "subject_id",
            "subject_accession",
            "genus",
            "species",
            "subtype",
            "virus_name",
            "reference_isolate",
            "pident",
            "align_length",
            "subject_coverage_pct",
            "bitscore",
            "evalue",
            "reference_header",
        ],
    )

    typed_count = sum(1 for row in summary_rows if row["status"] == "typed")
    print(f"Input FASTA: {input_fasta}")
    print(f"Processed sequences: {len(records)}")
    print(f"Typed sequences: {typed_count}")
    print(f"Summary TSV: {summary_tsv}")
    print(f"All hits TSV: {all_hits_tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
