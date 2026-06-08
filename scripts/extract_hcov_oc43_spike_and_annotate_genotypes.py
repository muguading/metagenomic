#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path("/Users/wuhhh/Desktop/徐老师/代码/metagenomic")
GENOME_DIR = ROOT / "database/virus/seasonal_coronavirus/HCoV_OC43_genomes"
DOCX_META_TSV = GENOME_DIR / "HCoV_OC43_docx_accessions_expanded.tsv"
SPIKE_FASTA = GENOME_DIR / "HCoV_OC43_spike_genes.fasta"
SPIKE_TSV = GENOME_DIR / "HCoV_OC43_spike_genes.tsv"
CONTIG_TSV = GENOME_DIR / "HCoV_OC43_contig_genotypes.tsv"


# Manual typing map curated from Supplement Figure 1B (Spike tree)
# in temi_a_2019560_sm7144.docx.
GENOTYPE_BY_ACCESSION: dict[str, tuple[str, str]] = {
    "MW532110": ("D", "assigned from Spike tree"),
    "MW532111": ("D", "assigned from Spike tree"),
    "MW532112": ("C", "assigned from Spike tree"),
    "MW532113": ("I", "assigned from Spike tree"),
    "MW532114": ("I", "assigned from Spike tree"),
    "MW532115": ("I", "assigned from Spike tree"),
    "MW532116": ("I", "assigned from Spike tree"),
    "MW532117": ("Novel genotype K", "assigned from Spike tree"),
    "MW532118": ("Novel genotype K", "assigned from Spike tree"),
    "MW532108": ("H", "assigned from Spike tree"),
    "MW532109": ("H", "assigned from Spike tree; figure label shows China/12/2018"),
    "MW532119": ("unresolved", "not clearly shown in Supplement Figure 1B Spike tree"),
    "OK318939": ("Novel genotype J", "assigned from Spike tree"),
    "OK318940": ("Novel genotype J", "assigned from Spike tree"),
    "OK318941": ("Novel genotype J", "assigned from Spike tree"),
    "OK391221": ("Novel genotype J", "assigned from Spike tree"),
    "OK318942": ("Novel genotype J", "assigned from Spike tree"),
    "OK391222": ("Novel genotype J", "assigned from Spike tree"),
    "OK391223": ("Novel genotype J", "assigned from Spike tree"),
    "OK391224": ("Novel genotype J", "assigned from Spike tree"),
    "OK391225": ("Novel genotype J", "assigned from Spike tree"),
    "OK318944": ("Novel genotype J", "assigned from Spike tree"),
    "OK391226": ("Novel genotype J", "assigned from Spike tree"),
    "OK318918": ("Novel genotype K", "assigned from Spike tree"),
    "OK318945": ("Novel genotype K", "assigned from Spike tree"),
    "OK391238": ("Novel genotype K", "assigned from Spike tree"),
    "OK391231": ("Novel genotype K", "assigned from Spike tree; expanded from docx range OK391231-37"),
    "OK391232": ("Novel genotype K", "assigned from Spike tree; expanded from docx range OK391231-37"),
    "OK391233": ("Novel genotype K", "assigned from Spike tree; expanded from docx range OK391231-37"),
    "OK391234": ("Novel genotype K", "assigned from Spike tree; expanded from docx range OK391231-37"),
    "OK391235": ("Novel genotype K", "assigned from Spike tree; expanded from docx range OK391231-37"),
    "OK391236": ("Novel genotype K", "assigned from Spike tree; expanded from docx range OK391231-37"),
    "OK391237": ("Novel genotype K", "assigned from Spike tree; expanded from docx range OK391231-37"),
    "OK391227": ("Novel genotype K", "assigned from Spike tree"),
    "OK391229": ("Novel genotype K", "assigned from Spike tree"),
    "OK318946": ("Novel genotype K", "assigned from Spike tree"),
    "OK318917": ("Novel genotype K", "assigned from Spike tree"),
    "OK318947": ("Novel genotype K", "assigned from Spike tree"),
    "OK500297": ("Novel genotype K", "assigned from Spike tree"),
    "OK500298": ("Novel genotype K", "assigned from Spike tree"),
    "OK500299": ("unlabeled", "present between Novel genotype K and genotype I in Spike tree; no explicit bracket label"),
    "OK500300": ("Novel genotype K", "assigned from Spike tree"),
    "OK500301": ("unlabeled", "present between Novel genotype K and genotype I in Spike tree; no explicit bracket label"),
    "OK500302": ("Novel genotype K", "assigned from Spike tree"),
    "OK500303": ("Novel genotype K", "assigned from Spike tree; figure label shows China/67/2020"),
    "OK500304": ("Novel genotype K", "assigned from Spike tree; figure label shows China/70/2020"),
    "OK500305": ("Novel genotype K", "assigned from Spike tree; figure label shows China/73/2020"),
    "OK391228": ("Novel genotype K", "assigned from Spike tree"),
    "OK391230": ("Novel genotype K", "assigned from Spike tree"),
}


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


def wrap_sequence(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def read_docx_meta(path: Path) -> dict[str, dict[str, str]]:
    meta: dict[str, dict[str, str]] = {}
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            meta[row["accession"]] = row
    return meta


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
                gene == "S"
                or name == "S"
                or "SPIKE" in combined
            ):
                gene_name = attrs.get("gene", attrs.get("Name", "")).strip()
                if "SPIKE" in combined:
                    gene_name = "S"
                spike_hit = {
                    "seqid": seqid,
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
        "seqid": region_attrs.get("ID", "").split(":")[0] if region_attrs else "",
        "region_strain": region_attrs.get("strain", ""),
        "region_country": region_attrs.get("country", ""),
        "region_collection_date": region_attrs.get("collection-date", ""),
    }
    if spike_hit:
        result.update(spike_hit)
    return result


def main() -> None:
    docx_meta = read_docx_meta(DOCX_META_TSV)
    fasta_paths = sorted(
        path
        for path in GENOME_DIR.glob("*.fasta")
        if path.name not in {SPIKE_FASTA.name}
    )

    contig_rows: list[dict[str, str]] = []
    spike_rows: list[dict[str, str]] = []

    with SPIKE_FASTA.open("w") as fasta_out:
        for fasta_path in fasta_paths:
            accession = fasta_path.stem
            gff_path = GENOME_DIR / f"{accession}.gff3"
            if not gff_path.exists():
                continue

            fasta_header, genome_seq = parse_fasta(fasta_path)
            gff_info = inspect_gff(gff_path)
            meta = docx_meta.get(accession, {})
            genotype, genotype_note = GENOTYPE_BY_ACCESSION.get(
                accession, ("unknown", "no manual assignment found from Spike tree")
            )

            contig_row = {
                "accession": accession,
                "fasta_header": fasta_header,
                "docx_strain_name": meta.get("strain_name", ""),
                "docx_accession_raw": meta.get("accession_raw", ""),
                "docx_collection_year": meta.get("collection_year", ""),
                "region_strain": gff_info.get("region_strain", ""),
                "region_collection_date": gff_info.get("region_collection_date", ""),
                "region_country": gff_info.get("region_country", ""),
                "tree_genotype": genotype,
                "tree_genotype_note": genotype_note,
                "has_spike_feature": "yes" if gff_info.get("start") else "no",
                "spike_feature_type": gff_info.get("feature_type", ""),
                "spike_start": gff_info.get("start", ""),
                "spike_end": gff_info.get("end", ""),
                "spike_strand": gff_info.get("strand", ""),
                "spike_gene_name": gff_info.get("gene", ""),
                "spike_product": gff_info.get("product", ""),
                "spike_partial": gff_info.get("partial", ""),
                "gff_path": str(gff_path),
                "fasta_path": str(fasta_path),
            }
            contig_rows.append(contig_row)

            if not gff_info.get("start"):
                continue

            spike_seq = extract_subseq(
                genome_seq,
                int(gff_info["start"]),
                int(gff_info["end"]),
                gff_info["strand"],
            )
            spike_header = (
                f"{accession}|tree_genotype={genotype}|strain={contig_row['docx_strain_name'] or contig_row['region_strain']}"
                f"|start={gff_info['start']}|end={gff_info['end']}|partial={gff_info.get('partial', 'false')}"
            )
            fasta_out.write(f">{spike_header}\n{wrap_sequence(spike_seq)}\n")

            spike_rows.append(
                {
                    "accession": accession,
                    "spike_fasta_header": spike_header,
                    "docx_strain_name": contig_row["docx_strain_name"],
                    "docx_accession_raw": contig_row["docx_accession_raw"],
                    "docx_collection_year": contig_row["docx_collection_year"],
                    "region_strain": contig_row["region_strain"],
                    "region_collection_date": contig_row["region_collection_date"],
                    "region_country": contig_row["region_country"],
                    "tree_genotype": contig_row["tree_genotype"],
                    "tree_genotype_note": contig_row["tree_genotype_note"],
                    "spike_feature_type": contig_row["spike_feature_type"],
                    "spike_start": contig_row["spike_start"],
                    "spike_end": contig_row["spike_end"],
                    "spike_strand": contig_row["spike_strand"],
                    "spike_gene_name": contig_row["spike_gene_name"],
                    "spike_product": contig_row["spike_product"],
                    "spike_partial": contig_row["spike_partial"],
                    "spike_length_nt": str(len(spike_seq)),
                    "gff_path": contig_row["gff_path"],
                    "fasta_path": contig_row["fasta_path"],
                }
            )

    with CONTIG_TSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "accession",
                "fasta_header",
                "docx_strain_name",
                "docx_accession_raw",
                "docx_collection_year",
                "region_strain",
                "region_collection_date",
                "region_country",
                "tree_genotype",
                "tree_genotype_note",
                "has_spike_feature",
                "spike_feature_type",
                "spike_start",
                "spike_end",
                "spike_strand",
                "spike_gene_name",
                "spike_product",
                "spike_partial",
                "gff_path",
                "fasta_path",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(contig_rows)

    with SPIKE_TSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "accession",
                "spike_fasta_header",
                "docx_strain_name",
                "docx_accession_raw",
                "docx_collection_year",
                "region_strain",
                "region_collection_date",
                "region_country",
                "tree_genotype",
                "tree_genotype_note",
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
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(spike_rows)

    print(f"Wrote contig genotype table: {CONTIG_TSV}")
    print(f"Wrote spike table: {SPIKE_TSV}")
    print(f"Wrote spike fasta: {SPIKE_FASTA}")
    print(f"Contigs processed: {len(contig_rows)}")
    print(f"Spike sequences extracted: {len(spike_rows)}")


if __name__ == "__main__":
    main()
