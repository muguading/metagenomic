#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


DEFAULT_FASTA = Path("database/virus/rsv/rsva_db.fasta")
DEFAULT_META = Path("database/virus/rsv/rsva_db.tsv")
DEFAULT_OUT_FASTA = Path("database/virus/rsv/rsva_db_c99.fasta")
DEFAULT_OUT_META = Path("database/virus/rsv/rsva_db_c99.tsv")
DEFAULT_SUMMARY = Path("database/virus/rsv/dedup_filtc99.summary.tsv")
DEFAULT_MMSEQS = Path("/opt/homebrew/Caskroom/mambaforge/base/envs/ncov/bin/mmseqs")


@dataclass(frozen=True)
class FastaRecord:
    fasta_id: str
    accession: str
    description: str
    sequence: str


@dataclass(frozen=True)
class MetaSequence:
    accession: str
    meta_key: str
    genotype: str
    meta_index: int
    fasta: FastaRecord


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Group RSV-A sequences by genotype, deduplicate each group at 0.99 identity with mmseqs, and export merged FASTA/TSV/summary.",
    )
    parser.add_argument("--fasta", type=Path, default=DEFAULT_FASTA, help="Input FASTA.")
    parser.add_argument("--meta", type=Path, default=DEFAULT_META, help="Input metadata TSV.")
    parser.add_argument("--out-fasta", type=Path, default=DEFAULT_OUT_FASTA, help="Output deduplicated FASTA.")
    parser.add_argument("--out-meta", type=Path, default=DEFAULT_OUT_META, help="Output deduplicated metadata TSV.")
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY, help="Output summary TSV.")
    parser.add_argument("--mmseqs", type=Path, default=DEFAULT_MMSEQS, help="Path to mmseqs binary.")
    parser.add_argument("--min-seq-id", type=float, default=0.99, help="mmseqs minimum sequence identity. Default: 0.99")
    parser.add_argument("--coverage", type=float, default=0.99, help="mmseqs coverage threshold (-c). Default: 0.99")
    parser.add_argument("--threads", type=int, default=8, help="mmseqs threads. Default: 8")
    return parser.parse_args()


def normalize_genotype(value: str | None) -> str:
    text = str(value or "").strip()
    return text if text else "Unknown"


def build_meta_key(row: dict[str, str]) -> str:
    seq_name = str(row.get("seqName") or "").strip()
    if seq_name:
        fasta_id, accession, _description = parse_fasta_header(seq_name)
        return accession or fasta_id
    accession = str(row.get("Accession") or "").strip()
    if accession:
        return accession.split(".", 1)[0]
    return ""


def parse_fasta_header(header: str) -> tuple[str, str, str]:
    first_token = header.split()[0]
    fasta_id = first_token.strip()
    accession = fasta_id.split("|", 1)[0].split(".", 1)[0]
    description = header[len(first_token):].strip()
    return fasta_id, accession, description


def read_fasta(path: Path) -> dict[str, FastaRecord]:
    records: dict[str, FastaRecord] = {}
    header = ""
    seq_chunks: list[str] = []
    with path.open(encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    fasta_id, accession, description = parse_fasta_header(header)
                    records[accession] = FastaRecord(
                        fasta_id=fasta_id,
                        accession=accession,
                        description=description,
                        sequence="".join(seq_chunks).upper(),
                    )
                header = line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line)
    if header:
        fasta_id, accession, description = parse_fasta_header(header)
        records[accession] = FastaRecord(
            fasta_id=fasta_id,
            accession=accession,
            description=description,
            sequence="".join(seq_chunks).upper(),
        )
    return records


def read_metadata(path: Path, fasta_records: dict[str, FastaRecord]) -> tuple[list[str], list[dict[str, str]], list[MetaSequence]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        has_accession = "Accession" in fieldnames
        has_seqname = "seqName" in fieldnames
        has_genotype = "Genotype" in fieldnames or "clade" in fieldnames
        if (not has_accession and not has_seqname) or not has_genotype:
            raise ValueError("Metadata TSV must contain Accession or seqName, plus Genotype or clade columns.")
        rows = list(reader)
    matched: list[MetaSequence] = []
    for index, row in enumerate(rows):
        accession = build_meta_key(row)
        if not accession:
            continue
        fasta = fasta_records.get(accession)
        if fasta is None:
            continue
        matched.append(
            MetaSequence(
                accession=accession,
                meta_key=accession,
                genotype=normalize_genotype(row.get("Genotype") or row.get("clade")),
                meta_index=index,
                fasta=fasta,
            )
        )
    return fieldnames, rows, matched


def write_group_fasta(path: Path, records: list[MetaSequence]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for record in records:
            header = record.fasta.fasta_id
            if record.fasta.description:
                header = f"{header} {record.fasta.description}"
            handle.write(f">{header}\n")
            sequence = record.fasta.sequence
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def read_rep_accessions(rep_fasta: Path) -> set[str]:
    accessions: set[str] = set()
    with rep_fasta.open(encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith(">"):
                _fasta_id, accession, _description = parse_fasta_header(line[1:].strip())
                accessions.add(accession)
    return accessions


def run_mmseqs_group(
    mmseqs_bin: Path,
    records: list[MetaSequence],
    min_seq_id: float,
    coverage: float,
    threads: int,
    tmp_root: Path,
) -> set[str]:
    if len(records) <= 1:
        return {record.accession for record in records}
    group_fasta = tmp_root / "input.fasta"
    cluster_prefix = tmp_root / "cluster"
    mmseqs_tmp = tmp_root / "mmseqs_tmp"
    write_group_fasta(group_fasta, records)
    cmd = [
        str(mmseqs_bin),
        "easy-cluster",
        str(group_fasta),
        str(cluster_prefix),
        str(mmseqs_tmp),
        "--dbtype",
        "2",
        "--min-seq-id",
        str(min_seq_id),
        "--cov-mode",
        "5",
        "-c",
        str(coverage),
        "--similarity-type",
        "2",
        "--cluster-mode",
        "2",
        "--threads",
        str(max(1, threads)),
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    rep_fasta = cluster_prefix.with_name(cluster_prefix.name + "_rep_seq.fasta")
    if not rep_fasta.is_file():
        raise FileNotFoundError(f"mmseqs representative FASTA not found: {rep_fasta}")
    return read_rep_accessions(rep_fasta)


def write_output_fasta(path: Path, records: list[MetaSequence]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for record in records:
            header = record.fasta.fasta_id
            if record.fasta.description:
                header = f"{header} {record.fasta.description}"
            handle.write(f">{header}\n")
            sequence = record.fasta.sequence
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def write_output_meta(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_summary(path: Path, rows: list[dict[str, str | int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["Genotype", "Input_Sequences", "Retained_Sequences", "Removed_Sequences", "Retention_Rate"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    mmseqs_bin = shutil.which(str(args.mmseqs)) if not args.mmseqs.is_absolute() else str(args.mmseqs)
    if not mmseqs_bin or not Path(mmseqs_bin).is_file():
        raise FileNotFoundError(f"mmseqs not found: {args.mmseqs}")

    fasta_records = read_fasta(args.fasta)
    fieldnames, meta_rows, matched_records = read_metadata(args.meta, fasta_records)
    if not matched_records:
        raise ValueError("No metadata records could be matched to FASTA accessions.")

    genotype_groups: dict[str, list[MetaSequence]] = defaultdict(list)
    for record in matched_records:
        genotype_groups[record.genotype].append(record)

    kept_accessions: set[str] = set()
    summary_rows: list[dict[str, str | int]] = []
    tmp_parent = Path.cwd() / "tmp"
    tmp_parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="rsva_mmseqs_c99_", dir=str(tmp_parent)) as tmp_dir:
        tmp_root = Path(tmp_dir)
        for genotype in sorted(genotype_groups, key=lambda value: (value == "Unknown", value)):
            group_records = genotype_groups[genotype]
            group_dir = tmp_root / sanitize_filename(genotype)
            group_dir.mkdir(parents=True, exist_ok=True)
            retained = run_mmseqs_group(
                Path(mmseqs_bin),
                group_records,
                min_seq_id=args.min_seq_id,
                coverage=args.coverage,
                threads=args.threads,
                tmp_root=group_dir,
            )
            kept_accessions.update(retained)
            input_count = len(group_records)
            retained_count = len(retained)
            summary_rows.append(
                {
                    "Genotype": genotype,
                    "Input_Sequences": input_count,
                    "Retained_Sequences": retained_count,
                    "Removed_Sequences": input_count - retained_count,
                    "Retention_Rate": f"{(retained_count / input_count):.4f}" if input_count else "0.0000",
                }
            )

    retained_rows = [
        row for row in meta_rows
        if build_meta_key(row) in kept_accessions
    ]
    matched_by_accession = {record.meta_key: record for record in matched_records}
    retained_records = [
        matched_by_accession[build_meta_key(row)]
        for row in retained_rows
    ]

    write_output_fasta(args.out_fasta, retained_records)
    write_output_meta(args.out_meta, fieldnames, retained_rows)
    total_input = len(matched_records)
    total_retained = len(retained_rows)
    summary_rows.append(
        {
            "Genotype": "__TOTAL__",
            "Input_Sequences": total_input,
            "Retained_Sequences": total_retained,
            "Removed_Sequences": total_input - total_retained,
            "Retention_Rate": f"{(total_retained / total_input):.4f}" if total_input else "0.0000",
        }
    )
    write_summary(args.summary, summary_rows)


def sanitize_filename(value: str) -> str:
    text = value.strip() or "Unknown"
    return "".join(ch if ch.isalnum() or ch in {"_", "-", "."} else "_" for ch in text)


if __name__ == "__main__":
    main()
