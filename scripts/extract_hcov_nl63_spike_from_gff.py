#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path("/Users/wuhhh/Desktop/徐老师/代码/metagenomic")
GENOME_DIR = ROOT / "database/virus/seasonal_coronavirus/HCoV_NL63_genomes"
OUTPUT_FASTA = GENOME_DIR / "HCoV_NL63_spike_genes.fasta"
OUTPUT_TSV = GENOME_DIR / "HCoV_NL63_spike_genes.tsv"


def parse_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    header: str | None = None
    chunks: list[str] = []
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records[header] = "".join(chunks)
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
    if header is not None:
        records[header] = "".join(chunks)
    return records


def parse_attributes(raw: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in raw.split(";"):
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        attrs[key] = value
    return attrs


def reverse_complement(seq: str) -> str:
    table = str.maketrans(
        {
            "A": "T",
            "C": "G",
            "G": "C",
            "T": "A",
            "R": "Y",
            "Y": "R",
            "M": "K",
            "K": "M",
            "B": "V",
            "V": "B",
            "D": "H",
            "H": "D",
            "N": "N",
            "a": "t",
            "c": "g",
            "g": "c",
            "t": "a",
            "r": "y",
            "y": "r",
            "m": "k",
            "k": "m",
            "b": "v",
            "v": "b",
            "d": "h",
            "h": "d",
            "n": "n",
        }
    )
    return seq.translate(table)[::-1]


def extract_subseq(seq: str, start: int, end: int, strand: str) -> str:
    fragment = seq[start - 1 : end]
    if strand == "-":
        return reverse_complement(fragment)
    return fragment


def is_spike_feature(feature_type: str, attrs: dict[str, str]) -> bool:
    if feature_type not in {"CDS", "gene"}:
        return False
    gene = attrs.get("gene", "").upper()
    name = attrs.get("Name", "").upper()
    product = attrs.get("product", "").upper()
    note = attrs.get("Note", "").upper()
    combined = " ".join([gene, name, product, note])
    return gene == "S" or name == "S" or "SPIKE" in combined


def find_spike_feature(gff_path: Path) -> tuple[int, int, str, str, str] | None:
    gene_hit: tuple[int, int, str, str, str] | None = None
    cds_hit: tuple[int, int, str, str, str] | None = None

    with gff_path.open() as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                continue
            seqid, source, feature_type, start, end, score, strand, phase, attributes = fields
            attrs = parse_attributes(attributes)
            if not is_spike_feature(feature_type, attrs):
                continue
            gene_name = attrs.get("gene", "").strip()
            if not gene_name and is_spike_feature(feature_type, attrs):
                gene_name = "S"
            hit = (
                int(start),
                int(end),
                strand,
                gene_name,
                attrs.get("product", attrs.get("Note", "")),
            )
            if feature_type == "CDS":
                cds_hit = hit
            elif gene_hit is None:
                gene_hit = hit

    return cds_hit or gene_hit


def wrap_sequence(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def main() -> None:
    fasta_paths = sorted(
        path for path in GENOME_DIR.glob("*.fasta") if path.name != OUTPUT_FASTA.name
    )
    rows: list[dict[str, str]] = []
    missing: list[str] = []

    with OUTPUT_FASTA.open("w") as fasta_out:
        for fasta_path in fasta_paths:
            accession = fasta_path.stem
            gff_path = GENOME_DIR / f"{accession}.gff3"
            if not gff_path.exists():
                missing.append(f"{accession}\tmissing_gff3")
                continue

            records = parse_fasta(fasta_path)
            if not records:
                missing.append(f"{accession}\tempty_fasta")
                continue

            header, genome_seq = next(iter(records.items()))
            spike_hit = find_spike_feature(gff_path)
            if spike_hit is None:
                missing.append(f"{accession}\tmissing_spike_annotation")
                continue

            start, end, strand, gene_name, product = spike_hit
            spike_seq = extract_subseq(genome_seq, start, end, strand)
            fasta_header = (
                f"{accession}|gene={gene_name or 'S'}|product={product or 'spike'}"
                f"|start={start}|end={end}|strand={strand}"
            )
            fasta_out.write(f">{fasta_header}\n{wrap_sequence(spike_seq)}\n")

            rows.append(
                {
                    "accession": accession,
                    "fasta_header": header,
                    "gff_path": str(gff_path),
                    "start": str(start),
                    "end": str(end),
                    "strand": strand,
                    "length_nt": str(len(spike_seq)),
                    "gene_name": gene_name,
                    "product": product,
                }
            )

    with OUTPUT_TSV.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "accession",
                "fasta_header",
                "gff_path",
                "start",
                "end",
                "strand",
                "length_nt",
                "gene_name",
                "product",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    if missing:
        missing_path = GENOME_DIR / "HCoV_NL63_spike_genes.missing.tsv"
        missing_path.write_text("accession\treason\n" + "\n".join(missing) + "\n")
    else:
        missing_path = GENOME_DIR / "HCoV_NL63_spike_genes.missing.tsv"
        if missing_path.exists():
            missing_path.unlink()

    print(f"Extracted {len(rows)} spike genes to {OUTPUT_FASTA}")
    print(f"Wrote manifest to {OUTPUT_TSV}")
    if missing:
        print(f"Missing annotations for {len(missing)} accessions: {missing_path}")


if __name__ == "__main__":
    main()
