#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path

from Bio import SeqIO


ROOT = Path("/Users/wuhhh/Desktop/徐老师/代码/metagenomic")
DEFAULT_QUERY = ROOT / "ev_db.fasta"
DEFAULT_REF_FASTA = ROOT / "database" / "virus" / "enterovirus" / "reference_genomes" / "abcd_vp1" / "enterovirus_abcd_vp1.fasta"
DEFAULT_REF_MANIFEST = ROOT / "database" / "virus" / "enterovirus" / "reference_genomes" / "abcd_vp1" / "enterovirus_abcd_vp1.tsv"
DEFAULT_OUT_DIR = ROOT / "tmp" / "enterovirus_vp1_typing_demo"
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


def load_manifest(path: Path) -> dict[str, dict[str, str]]:
    mapping: dict[str, dict[str, str]] = {}
    with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            accession = str(row.get("accession") or "").strip()
            accession_full = str(row.get("accession_full") or "").strip()
            keys = {
                accession,
                accession_full.split(".", 1)[0] if accession_full else "",
                str(row.get("header") or "").split()[0].split(".", 1)[0],
            }
            for key in keys:
                key = str(key or "").strip()
                if key:
                    mapping[key] = row
    return mapping


def write_subset_fasta(records: list, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, str(out_path), "fasta")


def read_query_records(path: Path, limit: int) -> list:
    records = list(SeqIO.parse(str(path), "fasta"))
    if limit > 0:
        return records[:limit]
    return records


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
                "big_group": str(meta.get("big_group") or "").strip(),
                "subtype": str(meta.get("abbrev") or "").strip(),
                "virus_name": str(meta.get("virus_name") or "").strip(),
                "pident": f"{pident:.2f}",
                "align_length": str(align_len),
                "subject_coverage_pct": f"{subject_cov:.2f}",
                "bitscore": f"{bitscore:.1f}",
                "evalue": str(row["evalue"]),
                "reference_header": str(meta.get("header") or ""),
            }
            all_rows.append(result)
            previous = best_by_query.get(result["query_id"])
            current_score = (bitscore, subject_cov, pident, align_len)
            if previous is None:
                result["_score"] = current_score
                best_by_query[result["query_id"]] = result
                continue
            previous_score = previous.get("_score") or (0.0, 0.0, 0.0, 0)
            if current_score > previous_score:
                result["_score"] = current_score
                best_by_query[result["query_id"]] = result
    return all_rows, best_by_query


def build_summary_rows(records: list, best_by_query: dict[str, dict[str, str]]) -> list[dict[str, str]]:
    summary_rows: list[dict[str, str]] = []
    for record in records:
        query_id = str(record.id)
        best = best_by_query.get(query_id)
        sequence_length = len(str(record.seq))
        if best is None:
            summary_rows.append(
                {
                    "query_id": query_id,
                    "sequence_length": str(sequence_length),
                    "big_group": "",
                    "subtype": "",
                    "virus_name": "",
                    "subject_accession": "",
                    "subject_coverage_pct": "0.00",
                    "pident": "0.00",
                    "bitscore": "0.0",
                    "status": "unassigned",
                    "note": "未命中 VP1 参考库",
                }
            )
            continue
        subject_cov = float(best["subject_coverage_pct"])
        pident = float(best["pident"])
        if subject_cov >= 80.0 and pident >= 80.0:
            status = "typed"
            note = "VP1 高置信命中"
        elif subject_cov >= 60.0 and pident >= 70.0:
            status = "review"
            note = "VP1 中等命中，建议复核"
        else:
            status = "weak"
            note = "VP1 弱命中，建议复核"
        summary_rows.append(
            {
                "query_id": query_id,
                "sequence_length": str(sequence_length),
                "big_group": best["big_group"],
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
            clean = {key: ("" if key.startswith("_") else row.get(key, "")) for key in fieldnames}
            writer.writerow(clean)


def main() -> int:
    parser = argparse.ArgumentParser(description="Type enterovirus complete genomes by local VP1 reference blastn")
    parser.add_argument("--input-fasta", "--input", dest="input_fasta", type=Path, default=DEFAULT_QUERY)
    parser.add_argument("--reference-fasta", type=Path, default=DEFAULT_REF_FASTA)
    parser.add_argument("--reference-manifest", type=Path, default=DEFAULT_REF_MANIFEST)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--blastn-bin", type=str, default="")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--limit", type=int, default=100, help="Only type the first N sequences; use 0 for all")
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
    blast_out = output_dir / "vp1_blastn.tsv"
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
            "big_group",
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
            "big_group",
            "subtype",
            "virus_name",
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
