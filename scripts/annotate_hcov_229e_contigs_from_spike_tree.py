#!/usr/bin/env python3
from __future__ import annotations

import csv
import re
from collections import Counter
from pathlib import Path


PROJECT_ROOT = Path("/Users/wuhhh/Desktop/徐老师/代码/metagenomic")
GENOMES_DIR = PROJECT_ROOT / "database/virus/seasonal_coronavirus/HCoV_229E_genomes"
TREE_TABLE = PROJECT_ROOT / "database/virus/seasonal_coronavirus/HCoV_229E_resolved_accessions.tsv"
OUT_CONTIG = GENOMES_DIR / "HCoV_229E_contig_spike_tree_annotations.tsv"
OUT_SAMPLE = GENOMES_DIR / "HCoV_229E_sample_spike_tree_annotations.tsv"


GENE_PATTERNS: list[tuple[str, str]] = [
    ("surface glycoprotein (S) gene", "spike"),
    ("nucleocapsid protein (N) gene", "nucleocapsid"),
    ("RNA-dependent RNA polymerase (RdRp) gene", "rdrp"),
]


def normalize_space(text: str) -> str:
    return re.sub(r"\s+", " ", text.strip())


def infer_sample_name_from_tree_label(label: str) -> str:
    parts = normalize_space(label).split(" ")
    if len(parts) < 4:
        return ""
    return " ".join(parts[1:-2])


def infer_gene_and_sample_from_header(header: str) -> tuple[str, str]:
    text = header[1:] if header.startswith(">") else header
    text = normalize_space(text)
    isolate_marker = " isolate "
    if isolate_marker not in text:
        return "unknown", ""
    after_isolate = text.split(isolate_marker, 1)[1]
    if ", complete genome" in after_isolate:
        sample = after_isolate.split(", complete genome", 1)[0].strip()
        return "complete_genome", sample
    for marker, gene in GENE_PATTERNS:
        token = f" {marker},"
        if token in after_isolate:
            sample = after_isolate.split(token, 1)[0].strip()
            return gene, sample
    sample = after_isolate.split(",", 1)[0].strip()
    return "unknown", sample


def load_tree_annotations() -> dict[str, dict[str, str]]:
    mapping: dict[str, dict[str, str]] = {}
    with TREE_TABLE.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            tree_label = normalize_space(row["#NAME?"])
            sample_name = infer_sample_name_from_tree_label(tree_label)
            if not sample_name:
                continue
            mapping[sample_name] = {
                "tree_label": tree_label,
                "genogroup": normalize_space(row["genogroup"]),
            }
    return mapping


def build_rows() -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    tree_map = load_tree_annotations()
    contig_rows: list[dict[str, str]] = []
    sample_map: dict[str, dict[str, str]] = {}

    for fasta_path in sorted(GENOMES_DIR.glob("*.fasta")):
        header = ""
        seq_len = 0
        with fasta_path.open() as handle:
            for line in handle:
                line = line.rstrip("\n")
                if line.startswith(">"):
                    header = line
                elif line:
                    seq_len += len(line.strip())
        accession = fasta_path.stem
        accession_version = header[1:].split(" ", 1)[0] if header else accession
        gene_name, sample_name = infer_gene_and_sample_from_header(header)
        tree_info = tree_map.get(sample_name)
        if tree_info:
            tree_label = tree_info["tree_label"]
            genogroup = tree_info["genogroup"]
            status = "matched"
        else:
            tree_label = ""
            genogroup = "unresolved"
            status = "sample_not_in_spike_tree"

        row = {
            "accession": accession,
            "accession_version": accession_version,
            "sample_name": sample_name,
            "gene_name": gene_name,
            "sequence_length_nt": str(seq_len),
            "genogroup": genogroup,
            "tree_label": tree_label,
            "match_status": status,
            "fasta_header": header,
            "source_file": fasta_path.name,
        }
        contig_rows.append(row)

        sample_entry = sample_map.setdefault(
            sample_name,
            {
                "sample_name": sample_name,
                "genogroup": genogroup,
                "tree_label": tree_label,
                "match_status": status,
                "genes_present": "",
                "accessions": "",
            },
        )
        sample_entry.setdefault("_genes", []).append(gene_name)
        sample_entry.setdefault("_accessions", []).append(accession)

    sample_rows: list[dict[str, str]] = []
    for sample_name in sorted(sample_map):
        entry = sample_map[sample_name]
        genes = [gene for gene in entry.pop("_genes") if gene]
        accessions = entry.pop("_accessions")
        entry["genes_present"] = ",".join(sorted(genes))
        entry["accessions"] = ",".join(accessions)
        sample_rows.append(entry)

    return contig_rows, sample_rows


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    contig_rows, sample_rows = build_rows()
    write_tsv(
        OUT_CONTIG,
        contig_rows,
        [
            "accession",
            "accession_version",
            "sample_name",
            "gene_name",
            "sequence_length_nt",
            "genogroup",
            "tree_label",
            "match_status",
            "fasta_header",
            "source_file",
        ],
    )
    write_tsv(
        OUT_SAMPLE,
        sample_rows,
        [
            "sample_name",
            "genogroup",
            "tree_label",
            "match_status",
            "genes_present",
            "accessions",
        ],
    )

    gene_counts = Counter(row["gene_name"] for row in contig_rows)
    matched_counts = Counter(row["match_status"] for row in contig_rows)
    print(f"Wrote {len(contig_rows)} contig annotations to {OUT_CONTIG}")
    print(f"Wrote {len(sample_rows)} sample annotations to {OUT_SAMPLE}")
    print("Gene counts:")
    for gene_name in sorted(gene_counts):
        print(f"  {gene_name}\t{gene_counts[gene_name]}")
    print("Match status counts:")
    for status in sorted(matched_counts):
        print(f"  {status}\t{matched_counts[status]}")


if __name__ == "__main__":
    main()
