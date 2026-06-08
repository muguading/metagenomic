#!/usr/bin/env python3
from __future__ import annotations

import csv
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


ROOT = Path("/Users/wuhhh/Desktop/徐老师/代码/metagenomic")
BASE_DIR = ROOT / "database/virus/seasonal_coronavirus/HCoV_OC43_genomes"
REF_DIR = BASE_DIR / "spike_tree_reference_genomes"

OC_SPIKE_TSV = BASE_DIR / "HCoV_OC43_spike_genes.tsv"
OC_SPIKE_FASTA = BASE_DIR / "HCoV_OC43_spike_genes.fasta"

REF_SPIKE_TSV = BASE_DIR / "HCoV_OC43_spike_tree_reference_spike_genes.tsv"
REF_SPIKE_FASTA = BASE_DIR / "HCoV_OC43_spike_tree_reference_spike_genes.fasta"
REF_META_TSV = BASE_DIR / "HCoV_OC43_spike_tree_reference_accessions.tsv"

MERGED_SPIKE_TSV = BASE_DIR / "HCoV_OC43_spike_genes_merged_with_references.tsv"
MERGED_SPIKE_FASTA = BASE_DIR / "HCoV_OC43_spike_genes_merged_with_references.fasta"


REFERENCE_ACCESSIONS: list[dict[str, str]] = [
    {"accession": "NC_006213.1", "tree_label": "ATCC VR-759", "tree_genotype": "A", "note": "reference from Spike tree; figure label omits underscore"},
    {"accession": "AY391777.1", "tree_label": "Human coronavirus OC43", "tree_genotype": "A", "note": "reference from Spike tree"},
    {"accession": "AY585229.1", "tree_label": "serotype OC43-Paris", "tree_genotype": "A", "note": "reference from Spike tree"},
    {"accession": "KU131570.1", "tree_label": "HCoV-OC43/UK/London/2011", "tree_genotype": "E", "note": "reference from Spike tree"},
    {"accession": "MT118671.1", "tree_label": "2018 8596", "tree_genotype": "E", "note": "reference from Spike tree"},
    {"accession": "KY967360.1", "tree_label": "HCoV OC43/Seattle/USA/SC2476/2015", "tree_genotype": "E", "note": "reference from Spike tree"},
    {"accession": "KF923888.1", "tree_label": "2145A/2010", "tree_genotype": "B", "note": "reference from Spike tree"},
    {"accession": "KF923886.1", "tree_label": "1908A/2010", "tree_genotype": "B", "note": "reference from Spike tree"},
    {"accession": "KF923889.1", "tree_label": "1926/2006", "tree_genotype": "B", "note": "reference from Spike tree"},
    {"accession": "KF923896.1", "tree_label": "3074A/2012", "tree_genotype": "B", "note": "reference from Spike tree"},
    {"accession": "KF923905.1", "tree_label": "229/2005", "tree_genotype": "C", "note": "reference from Spike tree"},
    {"accession": "JN129834.1", "tree_label": "HK04/01", "tree_genotype": "C", "note": "reference from Spike tree"},
    {"accession": "AY903460.1", "tree_label": "19572 Belgium 2004", "tree_genotype": "C", "note": "reference from Spike tree"},
    {"accession": "KF923900.1", "tree_label": "3647/2006", "tree_genotype": "C", "note": "reference from Spike tree"},
    {"accession": "KY554973.1", "tree_label": "N07-1689B 116X", "tree_genotype": "C", "note": "reference from Spike tree"},
    {"accession": "KY674920.1", "tree_label": "N09-595B", "tree_genotype": "D", "note": "reference from Spike tree"},
    {"accession": "KF923914.1", "tree_label": "5508/2007", "tree_genotype": "D", "note": "reference from Spike tree"},
    {"accession": "KF923923.1", "tree_label": "1892A/2008", "tree_genotype": "D", "note": "reference from Spike tree"},
    {"accession": "KF923897.1", "tree_label": "3269A/2012", "tree_genotype": "F", "note": "reference from Spike tree"},
    {"accession": "KF923903.1", "tree_label": "12691/2012", "tree_genotype": "F", "note": "reference from Spike tree"},
    {"accession": "KX538965.1", "tree_label": "MY-U208/12", "tree_genotype": "F", "note": "reference from Spike tree"},
    {"accession": "KX538974.1", "tree_label": "MY-U945/12", "tree_genotype": "F", "note": "reference from Spike tree"},
    {"accession": "MK303620.1", "tree_label": "MDS2", "tree_genotype": "G", "note": "reference from Spike tree"},
    {"accession": "KF923904.1", "tree_label": "12694/2012", "tree_genotype": "G", "note": "reference from Spike tree"},
    {"accession": "MG197710.1", "tree_label": "BJ-124", "tree_genotype": "G", "note": "reference from Spike tree"},
    {"accession": "KX538978.1", "tree_label": "MY-U1758/13", "tree_genotype": "G", "note": "reference from Spike tree"},
    {"accession": "KX538979.1", "tree_label": "MY-U1975/13", "tree_genotype": "G", "note": "reference from Spike tree"},
    {"accession": "MN306043.1", "tree_label": "HCoV OC43/Seattle/USA/SC0841/2019", "tree_genotype": "H", "note": "reference from Spike tree"},
    {"accession": "KP198610.1", "tree_label": "2058A/10", "tree_genotype": "unresolved", "note": "visible on Spike tree but bracket assignment is unclear"},
    {"accession": "MF374985.2", "tree_label": "HCoV-OC43/USA/TCNP 00212-2017", "tree_genotype": "I", "note": "reference from Spike tree"},
    {"accession": "KY967358.1", "tree_label": "HCoV OC43/Seattle/USA/SC2770/2015", "tree_genotype": "I", "note": "reference from Spike tree"},
    {"accession": "MN306041.1", "tree_label": "HCoV OC43/Seattle/USA/SC0810/2019", "tree_genotype": "Novel genotype K", "note": "reference from Spike tree"},
    {"accession": "MT118683.1", "tree_label": "2019 0522", "tree_genotype": "Novel genotype K", "note": "reference from Spike tree"},
    {"accession": "MT118692.1", "tree_label": "2019 1540", "tree_genotype": "Novel genotype K", "note": "reference from Spike tree"},
    {"accession": "MN306042.1", "tree_label": "HCoV OC43/Seattle/USA/SC0839/2019", "tree_genotype": "Novel genotype K", "note": "reference from Spike tree"},
    {"accession": "MZ450972.1", "tree_label": "OC43/Seattle/USA/20959/2021", "tree_genotype": "Novel genotype K", "note": "reference from Spike tree"},
]

DOWNLOAD_WORKERS = 6


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
    frag = seq[start - 1 : end]
    return reverse_complement(frag) if strand == "-" else frag


def wrap_sequence(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def fetch_text(accession: str, rettype: str) -> str:
    ids_to_try = [accession]
    if "." in accession:
        ids_to_try.append(accession.split(".", 1)[0])

    last_error: Exception | None = None
    for accession_id in ids_to_try:
        params = urllib.parse.urlencode(
            {"db": "nuccore", "id": accession_id, "rettype": rettype, "retmode": "text"}
        )
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?{params}"
        for attempt in range(5):
            req = urllib.request.Request(url, headers={"User-Agent": "metagenomic-codex/1.0"})
            try:
                with urllib.request.urlopen(req, timeout=120) as response:
                    text = response.read().decode("utf-8", errors="replace")
                if text.startswith("Error: Failed to understand id"):
                    raise ValueError(text.strip())
                return text
            except Exception as exc:  # network retries for flaky NCBI responses
                last_error = exc
                time.sleep(1.5 * (attempt + 1))
    if last_error is not None:
        raise last_error
    raise RuntimeError(f"failed to fetch {accession} {rettype}")


def download_one(entry: dict[str, str]) -> None:
    accession = entry["accession"]
    fasta_path = REF_DIR / f"{accession}.fasta"
    gff_path = REF_DIR / f"{accession}.gff3"

    if not fasta_path.exists():
        fasta_text = fetch_text(accession, "fasta")
        fasta_path.write_text(fasta_text, encoding="utf-8")
        time.sleep(0.2)
    if not gff_path.exists():
        gff_text = fetch_text(accession, "gff3")
        gff_path.write_text(gff_text, encoding="utf-8")
        time.sleep(0.2)


def ensure_reference_downloads() -> None:
    REF_DIR.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=DOWNLOAD_WORKERS) as executor:
        futures = [executor.submit(download_one, entry) for entry in REFERENCE_ACCESSIONS]
        for future in as_completed(futures):
            future.result()


def parse_fasta(path: Path) -> tuple[str, str]:
    header = ""
    chunks: list[str] = []
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith(">"):
                if header:
                    break
                header = line[1:].strip()
            elif header:
                chunks.append(line.strip())
    return header, "".join(chunks)


def inspect_gff(path: Path) -> dict[str, str]:
    region_attrs: dict[str, str] = {}
    spike_hit: dict[str, str] | None = None
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                continue
            seqid, source, feature_type, start, end, score, strand, phase, attributes = fields
            attrs = parse_attributes(attributes)
            if feature_type == "region" and not region_attrs:
                region_attrs = attrs

            product = attrs.get("product", "").upper()
            gene = attrs.get("gene", "").upper()
            name = attrs.get("Name", "").upper()
            note = attrs.get("Note", "").upper()
            combined = " ".join([gene, name, product, note])
            if feature_type in {"gene", "CDS", "sequence_feature"} and (
                gene == "S" or name == "S" or "SPIKE" in combined
                or product == "S"
            ):
                gene_name = attrs.get("gene", attrs.get("Name", "")).strip()
                if "SPIKE" in combined:
                    gene_name = "S"
                spike_hit = {
                    "feature_type": feature_type,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "gene": gene_name,
                    "product": attrs.get("product", ""),
                    "partial": attrs.get("partial", "false"),
                }
                if feature_type == "CDS":
                    break

    result = {
        "region_strain": region_attrs.get("strain", ""),
        "region_collection_date": region_attrs.get("collection-date", ""),
        "region_country": region_attrs.get("country", ""),
    }
    if spike_hit:
        result.update(spike_hit)
    return result


def build_reference_spike_outputs() -> None:
    with REF_META_TSV.open("w", newline="", encoding="utf-8") as meta_handle, \
        REF_SPIKE_TSV.open("w", newline="", encoding="utf-8") as tsv_handle, \
        REF_SPIKE_FASTA.open("w", encoding="utf-8") as fasta_handle:
        meta_writer = csv.DictWriter(
            meta_handle,
            fieldnames=["accession", "tree_label", "tree_genotype", "note"],
            delimiter="\t",
        )
        meta_writer.writeheader()
        meta_writer.writerows(REFERENCE_ACCESSIONS)

        spike_writer = csv.DictWriter(
            tsv_handle,
            fieldnames=[
                "accession",
                "spike_fasta_header",
                "tree_label",
                "tree_genotype",
                "tree_genotype_note",
                "region_strain",
                "region_collection_date",
                "region_country",
                "spike_feature_type",
                "spike_start",
                "spike_end",
                "spike_strand",
                "spike_gene_name",
                "spike_product",
                "spike_partial",
                "spike_length_nt",
                "gff_path",
                "fasta_path",
                "source_set",
            ],
            delimiter="\t",
        )
        spike_writer.writeheader()

        for entry in REFERENCE_ACCESSIONS:
            accession = entry["accession"]
            fasta_path = REF_DIR / f"{accession}.fasta"
            gff_path = REF_DIR / f"{accession}.gff3"
            fasta_header, genome_seq = parse_fasta(fasta_path)
            gff_info = inspect_gff(gff_path)
            if not gff_info.get("start"):
                continue

            spike_seq = extract_subseq(
                genome_seq,
                int(gff_info["start"]),
                int(gff_info["end"]),
                gff_info["strand"],
            )
            spike_header = (
                f"{accession}|tree_genotype={entry['tree_genotype']}|tree_label={entry['tree_label']}"
                f"|start={gff_info['start']}|end={gff_info['end']}|partial={gff_info.get('partial', 'false')}"
            )
            fasta_handle.write(f">{spike_header}\n{wrap_sequence(spike_seq)}\n")

            spike_writer.writerow(
                {
                    "accession": accession,
                    "spike_fasta_header": spike_header,
                    "tree_label": entry["tree_label"],
                    "tree_genotype": entry["tree_genotype"],
                    "tree_genotype_note": entry["note"],
                    "region_strain": gff_info.get("region_strain", ""),
                    "region_collection_date": gff_info.get("region_collection_date", ""),
                    "region_country": gff_info.get("region_country", ""),
                    "spike_feature_type": gff_info.get("feature_type", ""),
                    "spike_start": gff_info.get("start", ""),
                    "spike_end": gff_info.get("end", ""),
                    "spike_strand": gff_info.get("strand", ""),
                    "spike_gene_name": gff_info.get("gene", ""),
                    "spike_product": gff_info.get("product", ""),
                    "spike_partial": gff_info.get("partial", ""),
                    "spike_length_nt": str(len(spike_seq)),
                    "gff_path": str(gff_path),
                    "fasta_path": str(fasta_path),
                    "source_set": "reference_non_oc_label",
                }
            )


def merge_outputs() -> None:
    merged_rows: list[dict[str, str]] = []
    for path in [OC_SPIKE_TSV, REF_SPIKE_TSV]:
        with path.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                row = dict(row)
                if "source_set" not in row:
                    row["source_set"] = "study_oc43"
                merged_rows.append(row)

    merged_fieldnames = [
        "accession",
        "spike_fasta_header",
        "docx_strain_name",
        "docx_accession_raw",
        "docx_collection_year",
        "tree_label",
        "tree_genotype",
        "tree_genotype_note",
        "region_strain",
        "region_collection_date",
        "region_country",
        "spike_feature_type",
        "spike_start",
        "spike_end",
        "spike_strand",
        "spike_gene_name",
        "spike_product",
        "spike_partial",
        "spike_length_nt",
        "gff_path",
        "fasta_path",
        "source_set",
    ]

    with MERGED_SPIKE_TSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=merged_fieldnames, delimiter="\t")
        writer.writeheader()
        for row in merged_rows:
            writer.writerow({field: row.get(field, "") for field in merged_fieldnames})

    with MERGED_SPIKE_FASTA.open("w", encoding="utf-8") as out_handle:
        for path in [OC_SPIKE_FASTA, REF_SPIKE_FASTA]:
            out_handle.write(path.read_text(encoding="utf-8"))


def main() -> None:
    ensure_reference_downloads()
    build_reference_spike_outputs()
    merge_outputs()
    print(f"Wrote reference metadata: {REF_META_TSV}")
    print(f"Wrote reference spike TSV: {REF_SPIKE_TSV}")
    print(f"Wrote reference spike FASTA: {REF_SPIKE_FASTA}")
    print(f"Wrote merged spike TSV: {MERGED_SPIKE_TSV}")
    print(f"Wrote merged spike FASTA: {MERGED_SPIKE_FASTA}")


if __name__ == "__main__":
    main()
