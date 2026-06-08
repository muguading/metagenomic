#!/usr/bin/env python3
from __future__ import annotations

import csv
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


ROOT = Path("/Users/wuhhh/Desktop/徐老师/代码/metagenomic")
BASE_DIR = ROOT / "database/virus/seasonal_coronavirus/HCoV_HKU1_genomes"
DOWNLOAD_DIR = BASE_DIR / "tree_reference_genomes"

ACCESSION_TSV = BASE_DIR / "HCoV_HKU1_tree_accessions.tsv"
SPIKE_TSV = BASE_DIR / "HCoV_HKU1_spike_genes.tsv"
SPIKE_FASTA = BASE_DIR / "HCoV_HKU1_spike_genes.fasta"
MISSING_TSV = BASE_DIR / "HCoV_HKU1_spike_missing.tsv"

DOWNLOAD_WORKERS = 6


TREE_ACCESSIONS: list[dict[str, str]] = [
    {"accession": "KY983584", "tree_label": "SC2628/USA/2015", "subtype": "A"},
    {"accession": "LC315650", "tree_label": "Tokyo/SGH-15/2014/JPN/2014", "subtype": "A"},
    {"accession": "KF686344", "tree_label": "HKU1/human/USA/HKU1-15/2009", "subtype": "A"},
    {"accession": "KT779555", "tree_label": "BJ01-p3/CHN/2009", "subtype": "A"},
    {"accession": "KT779556", "tree_label": "BJ01-p9/CHN/2009", "subtype": "A"},
    {"accession": "DQ415900", "tree_label": "N23/CHN/2005", "subtype": "A"},
    {"accession": "DQ415901", "tree_label": "N24/CHN/2005", "subtype": "A"},
    {"accession": "DQ415909", "tree_label": "N13/CHN/2004", "subtype": "A"},
    {"accession": "DQ415904", "tree_label": "N6/CHN/2004", "subtype": "A"},
    {"accession": "DQ415905", "tree_label": "N7/CHN/2004", "subtype": "A"},
    {"accession": "DQ415906", "tree_label": "N9/CHN/2004", "subtype": "A"},
    {"accession": "DQ415907", "tree_label": "N10/CHN/2004", "subtype": "A"},
    {"accession": "DQ415908", "tree_label": "N11/CHN/2004", "subtype": "A"},
    {"accession": "DQ415910", "tree_label": "N14/CHN/2004", "subtype": "A"},
    {"accession": "DQ415914", "tree_label": "N18/CHN/2004", "subtype": "A"},
    {"accession": "DQ415896", "tree_label": "N19/CHN/2004", "subtype": "A"},
    {"accession": "AY597011", "tree_label": "HKU1 genotype A/CHN/2004", "subtype": "A"},
    {"accession": "DQ415903", "tree_label": "N3/CHN/2003", "subtype": "A"},
    {"accession": "KF686345", "tree_label": "HKU1/human/USA/HKU1-20/2010", "subtype": "A"},
    {"accession": "KF430201", "tree_label": "HKU1/human/USA/HKU1-18/2010", "subtype": "A"},
    {"accession": "KF686341", "tree_label": "HKU1/human/USA/HKU1-10/2010", "subtype": "A"},
    {"accession": "KF686339", "tree_label": "HKU1/human/USA/HKU1-3/2009", "subtype": "A"},
    {"accession": "KF686340", "tree_label": "HKU1/human/USA/HKU1-5/2009", "subtype": "A"},
    {"accession": "KY674942", "tree_label": "N09-1627B/USA/2016", "subtype": "A"},
    {"accession": "KY674943", "tree_label": "N09-1605B/USA/2016", "subtype": "A"},
    {"accession": "KY674941", "tree_label": "N09-1663B/USA/2016", "subtype": "A"},
    {"accession": "KF430202", "tree_label": "HKU1/human/USA/HKU1-7/2010", "subtype": "A"},
    {"accession": "KF686346", "tree_label": "HKU1/human/USA/HKU1-12/2010", "subtype": "A"},
    {"accession": "KF686342", "tree_label": "HKU1/human/USA/HKU1-11/2009", "subtype": "A"},
    {"accession": "KF686343", "tree_label": "HKU1/human/USA/HKU1-13/2010", "subtype": "A"},
    {"accession": "HM034837", "tree_label": "Caen1/FRA/2005", "subtype": "A"},
    {"accession": "DQ415912", "tree_label": "N16/CHN/2004", "subtype": "C"},
    {"accession": "DQ415897", "tree_label": "N20/CHN/2004", "subtype": "C"},
    {"accession": "DQ415913", "tree_label": "N17/CHN/2004", "subtype": "C"},
    {"accession": "DQ415899", "tree_label": "N22/CHN/2005", "subtype": "C"},
    {"accession": "DQ339101", "tree_label": "N5P8/CHN/2004", "subtype": "C"},
    {"accession": "DQ415898", "tree_label": "N21/CHN/2004", "subtype": "C"},
    {"accession": "MK167038", "tree_label": "SC2521/USA/2017", "subtype": "B"},
    {"accession": "AY884001", "tree_label": "HKU1 genotype B/CHN/2003", "subtype": "B"},
    {"accession": "MH940245", "tree_label": "SI17244/THA/2017", "subtype": "B"},
    {"accession": "DQ415911", "tree_label": "N15/CHN/2004", "subtype": "B"},
    {"accession": "DQ415902", "tree_label": "N25/CHN/2005", "subtype": "B"},
    {"accession": "KF686338", "tree_label": "HKU1/human/USA/HKU1-1/2005", "subtype": "B"},
    {"accession": "KY674921", "tree_label": "N08-87/USA/2016", "subtype": "B"},
    {"accession": "LC315651", "tree_label": "Tokyo/SGH-18/2016/JPN/2016", "subtype": "B"},
    {"accession": "LC654447", "tree_label": "Fukushima H815/JPN/2020", "subtype": "B"},
    {"accession": "LC654448", "tree_label": "Fukushima H821/JPN/2020", "subtype": "B"},
    {"accession": "LC654449", "tree_label": "Fukushima O943/JPN/2020", "subtype": "B"},
]


def fetch_text(accession: str, rettype: str) -> str:
    params = urllib.parse.urlencode(
        {"db": "nuccore", "id": accession, "rettype": rettype, "retmode": "text"}
    )
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?{params}"
    last_error: Exception | None = None
    for attempt in range(5):
        request = urllib.request.Request(url, headers={"User-Agent": "metagenomic-codex/1.0"})
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                text = response.read().decode("utf-8", errors="replace")
            if text.startswith("Error: Failed to understand id"):
                raise ValueError(text.strip())
            return text
        except Exception as exc:
            last_error = exc
            time.sleep(1.5 * (attempt + 1))
    if last_error is not None:
        raise last_error
    raise RuntimeError(f"failed to fetch {accession} {rettype}")


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


def read_genbank_spike(accession: str) -> dict[str, str] | None:
    text = fetch_text(accession, "gbwithparts")
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
                qualifiers[key] = value.strip().strip('"')
            else:
                qualifiers[body[1:]] = ""
            continue
        if line.startswith(" " * 21):
            continuation = line[21:].rstrip()
            if qualifiers:
                last_key = next(reversed(qualifiers))
                qualifiers[last_key] += continuation.strip().strip('"')
            else:
                current_location += continuation.strip()

    if current_key is not None:
        features.append((current_key, current_location.strip(), qualifiers))

    for key, location, qualifiers in features:
        gene = qualifiers.get("gene", "").upper()
        product = qualifiers.get("product", "").upper()
        note = qualifiers.get("note", "").upper()
        combined = " ".join([gene, product, note])
        if key not in {"gene", "CDS", "misc_feature"}:
            continue
        if gene == "S" or "SPIKE" in combined or product == "S":
            coords = parse_genbank_location(location)
            if coords is None:
                continue
            start, end, strand = coords
            return {
                "feature_type": key,
                "start": str(start),
                "end": str(end),
                "strand": strand,
                "gene": "S",
                "product": qualifiers.get("product", ""),
                "partial": "false",
            }
    return None


def download_one(entry: dict[str, str]) -> None:
    accession = entry["accession"]
    fasta_path = DOWNLOAD_DIR / f"{accession}.fasta"
    gff_path = DOWNLOAD_DIR / f"{accession}.gff3"
    if not fasta_path.exists():
        fasta_path.write_text(fetch_text(accession, "fasta"), encoding="utf-8")
        time.sleep(0.2)
    if not gff_path.exists():
        gff_path.write_text(fetch_text(accession, "gff3"), encoding="utf-8")
        time.sleep(0.2)


def ensure_downloads() -> None:
    DOWNLOAD_DIR.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=DOWNLOAD_WORKERS) as executor:
        futures = [executor.submit(download_one, entry) for entry in TREE_ACCESSIONS]
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
            "A": "T", "C": "G", "G": "C", "T": "A",
            "R": "Y", "Y": "R", "M": "K", "K": "M",
            "B": "V", "V": "B", "D": "H", "H": "D", "N": "N",
            "a": "t", "c": "g", "g": "c", "t": "a",
            "r": "y", "y": "r", "m": "k", "k": "m",
            "b": "v", "v": "b", "d": "h", "h": "d", "n": "n",
        }
    )
    return seq.translate(table)[::-1]


def extract_subseq(seq: str, start: int, end: int, strand: str) -> str:
    fragment = seq[start - 1 : end]
    return reverse_complement(fragment) if strand == "-" else fragment


def wrap_sequence(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


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
                gene == "S" or name == "S" or "SPIKE" in combined or product == "S"
            ):
                gene_name = attrs.get("gene", attrs.get("Name", "")).strip()
                if "SPIKE" in combined or product == "S":
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
    elif path.stem:
        fallback = read_genbank_spike(path.stem)
        if fallback:
            result.update(fallback)
    return result


def main() -> None:
    ensure_downloads()

    with ACCESSION_TSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["accession", "tree_label", "subtype"], delimiter="\t")
        writer.writeheader()
        writer.writerows(TREE_ACCESSIONS)

    missing: list[dict[str, str]] = []
    with SPIKE_TSV.open("w", newline="", encoding="utf-8") as tsv_handle, SPIKE_FASTA.open("w", encoding="utf-8") as fasta_handle:
        writer = csv.DictWriter(
            tsv_handle,
            fieldnames=[
                "accession",
                "spike_fasta_header",
                "tree_label",
                "subtype",
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
            ],
            delimiter="\t",
        )
        writer.writeheader()

        for entry in TREE_ACCESSIONS:
            accession = entry["accession"]
            fasta_path = DOWNLOAD_DIR / f"{accession}.fasta"
            gff_path = DOWNLOAD_DIR / f"{accession}.gff3"
            fasta_header, genome_seq = parse_fasta(fasta_path)
            gff_info = inspect_gff(gff_path)
            if not gff_info.get("start"):
                missing.append(entry)
                continue
            spike_seq = extract_subseq(genome_seq, int(gff_info["start"]), int(gff_info["end"]), gff_info["strand"])
            spike_header = (
                f"{accession}|subtype={entry['subtype']}|tree_label={entry['tree_label']}"
                f"|start={gff_info['start']}|end={gff_info['end']}|partial={gff_info.get('partial', 'false')}"
            )
            fasta_handle.write(f">{spike_header}\n{wrap_sequence(spike_seq)}\n")
            writer.writerow(
                {
                    "accession": accession,
                    "spike_fasta_header": spike_header,
                    "tree_label": entry["tree_label"],
                    "subtype": entry["subtype"],
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
                }
            )

    if missing:
        with MISSING_TSV.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=["accession", "tree_label", "subtype"], delimiter="\t")
            writer.writeheader()
            writer.writerows(missing)
    elif MISSING_TSV.exists():
        MISSING_TSV.unlink()

    print(f"Wrote accession table: {ACCESSION_TSV}")
    print(f"Wrote spike TSV: {SPIKE_TSV}")
    print(f"Wrote spike FASTA: {SPIKE_FASTA}")
    print(f"Accessions in tree: {len(TREE_ACCESSIONS)}")
    print(f"Spike sequences extracted: {len(TREE_ACCESSIONS) - len(missing)}")
    print(f"Missing spike records: {len(missing)}")


if __name__ == "__main__":
    main()
