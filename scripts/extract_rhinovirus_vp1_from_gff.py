#!/usr/bin/env python3
from __future__ import annotations

import csv
import urllib.request
from pathlib import Path


ROOT = Path("/Users/wuhhh/Desktop/徐老师/代码/metagenomic")
REFERENCE_DIR = ROOT / "database/virus/rhinovirus/reference_genomes"
GENOME_FASTA = REFERENCE_DIR / "human_rhinovirus_representative_genomes.fasta"
MANIFEST = REFERENCE_DIR / "human_rhinovirus_representative_genomes.tsv"
OUTPUT_FASTA = REFERENCE_DIR / "human_rhinovirus_vp1_representative_genes.fasta"
OUTPUT_MANIFEST = REFERENCE_DIR / "human_rhinovirus_vp1_representative_genes.tsv"
GENBANK_CACHE_DIR = REFERENCE_DIR / "genbank_cache"


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
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
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


def is_vp1_feature(feature_type: str, attrs: dict[str, str]) -> bool:
    if feature_type != "mature_protein_region_of_CDS":
        return False
    product = attrs.get("product", "").upper()
    note = attrs.get("Note", "").upper()
    combined = f"{product} {note}"
    return "VP1" in combined or product == "1D" or "1D (VP1)" in combined


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
            "N": "N",
            "n": "n",
        }
    )
    return seq.translate(table)[::-1]


def extract_subseq(seq: str, start: int, end: int, strand: str) -> str:
    frag = seq[start - 1 : end]
    if strand == "-":
        return reverse_complement(frag)
    return frag


def read_gff_vp1_feature(gff_path: Path) -> tuple[int, int, str, str, str] | None:
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
            if not is_vp1_feature(feature_type, attrs):
                continue
            return (
                int(start),
                int(end),
                strand,
                attrs.get("product", ""),
                attrs.get("Note", ""),
            )
    return None


def fetch_genbank_text(accession: str) -> str:
    GENBANK_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_path = GENBANK_CACHE_DIR / f"{accession}.gb"
    if cache_path.exists():
        return cache_path.read_text()
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=nucleotide&id={accession}&rettype=gbwithparts&retmode=text"
    )
    with urllib.request.urlopen(url, timeout=60) as response:
        text = response.read().decode("utf-8")
    cache_path.write_text(text)
    return text


def parse_genbank_qualifier_value(raw: str) -> str:
    raw = raw.strip()
    if raw.startswith('"') and raw.endswith('"'):
        return raw[1:-1]
    return raw


def parse_genbank_location(raw: str) -> tuple[int, int, str] | None:
    text = raw.strip()
    strand = "+"
    if text.startswith("complement(") and text.endswith(")"):
        strand = "-"
        text = text[len("complement(") : -1]
    if ".." not in text:
        return None
    start_raw, end_raw = text.split("..", 1)
    start = "".join(ch for ch in start_raw if ch.isdigit())
    end = "".join(ch for ch in end_raw if ch.isdigit())
    if not start or not end:
        return None
    return int(start), int(end), strand


def read_genbank_vp1_feature(accession: str) -> tuple[int, int, str, str, str] | None:
    text = fetch_genbank_text(accession)
    lines = text.splitlines()
    in_features = False
    current_key: str | None = None
    current_location = ""
    qualifiers: dict[str, str] = {}
    features: list[tuple[str, str, dict[str, str]]] = []

    for line in lines:
        if line.startswith("FEATURES             Location/Qualifiers"):
            in_features = True
            continue
        if not in_features:
            continue
        if line.startswith("ORIGIN"):
            break
        if len(line) >= 21 and line[:5] == "     " and line[5:21].strip():
            if current_key is not None:
                features.append((current_key, current_location.strip(), qualifiers))
            current_key = line[5:21].strip()
            current_location = line[21:].rstrip()
            qualifiers = {}
            continue
        if current_key is None:
            continue
        if line.startswith(" " * 21 + "/"):
            body = line[21:].strip()
            if "=" in body:
                key, value = body[1:].split("=", 1)
                qualifiers[key] = parse_genbank_qualifier_value(value)
            else:
                qualifiers[body[1:]] = ""
            continue
        if line.startswith(" " * 21):
            continuation = line[21:].rstrip()
            if current_location and not current_location.endswith(","):
                current_location += continuation.strip()
            elif qualifiers:
                last_key = next(reversed(qualifiers))
                qualifiers[last_key] += continuation.strip().strip('"')
            else:
                current_location += continuation.strip()

    if current_key is not None:
        features.append((current_key, current_location.strip(), qualifiers))

    for key, location, qualifiers in features:
        if key != "mat_peptide":
            continue
        product = qualifiers.get("product", "")
        note = qualifiers.get("note", "")
        combined = f"{product} {note}".upper()
        if "VP1" not in combined and product.upper() != "1D":
            continue
        coords = parse_genbank_location(location)
        if coords is None:
            continue
        start, end, strand = coords
        return start, end, strand, product, note
    return None


def main() -> None:
    genomes = parse_fasta(GENOME_FASTA)

    with MANIFEST.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    output_rows: list[dict[str, str]] = []
    missing: list[str] = []

    with OUTPUT_FASTA.open("w") as fasta_out:
        for row in rows:
            record_id = row["record_id"]
            accession = row["accession"]
            accession_root = row["accession_root"]
            gff_path = Path(row["gff_path"])
            genome_seq = genomes[record_id + f" {row['species_name']} type {row['normalized_type']} accession {accession}"]

            vp1_hit = read_gff_vp1_feature(gff_path)
            annotation_source = "gff3"
            if vp1_hit is None:
                vp1_hit = read_genbank_vp1_feature(accession)
                annotation_source = "genbank"

            if vp1_hit is None:
                missing.append(record_id)
                continue

            start, end, strand, product, note = vp1_hit
            vp1_seq = extract_subseq(genome_seq, start, end, strand)
            vp1_record_id = f"{record_id}_VP1"
            fasta_out.write(
                f">{vp1_record_id} {row['species_name']} {row['normalized_type']} VP1 accession {accession} "
                f"coords {start}-{end} strand {strand}\n"
            )
            for idx in range(0, len(vp1_seq), 80):
                fasta_out.write(vp1_seq[idx : idx + 80] + "\n")

            output_rows.append(
                {
                    "vp1_record_id": vp1_record_id,
                    "record_id": record_id,
                    "species_group": row["species_group"],
                    "species_name": row["species_name"],
                    "type_label": row["type_label"],
                    "normalized_type": row["normalized_type"],
                    "accession": accession,
                    "accession_root": accession_root,
                    "vp1_start": str(start),
                    "vp1_end": str(end),
                    "strand": strand,
                    "vp1_length_nt": str(len(vp1_seq)),
                    "product": product,
                    "note": note,
                    "annotation_source": annotation_source,
                    "genome_fasta_path": str(GENOME_FASTA),
                    "gff_path": str(gff_path),
                }
            )

    fieldnames = [
        "vp1_record_id",
        "record_id",
        "species_group",
        "species_name",
        "type_label",
        "normalized_type",
        "accession",
        "accession_root",
        "vp1_start",
        "vp1_end",
        "strand",
        "vp1_length_nt",
        "product",
        "note",
        "annotation_source",
        "genome_fasta_path",
        "gff_path",
    ]
    with OUTPUT_MANIFEST.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)

    if missing:
        raise SystemExit(f"Missing VP1 features for {len(missing)} records: {', '.join(missing[:10])}")


if __name__ == "__main__":
    main()
