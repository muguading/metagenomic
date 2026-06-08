#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

from Bio import SeqIO


ROOT = Path("/Users/wuhhh/Desktop/徐老师/代码/metagenomic")
DB_DIR = ROOT / "database" / "bacteria" / "yersinia_enterocolitica" / "o_antigen"
TMPDB_DIR = DB_DIR / "tmpdb"


EMBL_MARKERS = [
    {
        "source": "FR729477.1.embl",
        "marker_id": "YE_O3_oag_export",
        "serotype": "O:3",
        "gene": "oag_export",
        "locus_tag": "Y11_09840",
        "start": 1022852,
        "end": 1024108,
        "strand": "+",
        "product": "membrane protein involved in the export of O-antigen, teichoic acid lipoteichoic acids",
    },
    {
        "source": "FR729477.1.embl",
        "marker_id": "YE_O3_gtr_like",
        "serotype": "O:3",
        "gene": "gtr_like",
        "locus_tag": "Y11_09850",
        "start": 1024101,
        "end": 1025057,
        "strand": "+",
        "product": "beta-1,3-glucosyltransferase",
    },
    {
        "source": "AM286415.1.embl",
        "marker_id": "YE_O8_wzz",
        "serotype": "O:8",
        "gene": "wzz",
        "locus_tag": "YE3070",
        "start": 3339012,
        "end": 3340103,
        "strand": "-",
        "product": "O-antigen chain length determinant",
    },
    {
        "source": "AM286415.1.embl",
        "marker_id": "YE_O8_wbcJ",
        "serotype": "O:8",
        "gene": "wbcJ",
        "locus_tag": "YE3074",
        "start": 3344007,
        "end": 3344972,
        "strand": "-",
        "product": "GDP-fucose synthetase",
    },
    {
        "source": "AM286415.1.embl",
        "marker_id": "YE_O8_wzx",
        "serotype": "O:8",
        "gene": "wzx",
        "locus_tag": "YE3082",
        "start": 3352027,
        "end": 3353316,
        "strand": "-",
        "product": "Wzx flippase",
    },
    {
        "source": "CP002246.1.embl",
        "marker_id": "YE_O9_gmd",
        "serotype": "O:9",
        "gene": "gmd",
        "locus_tag": "YE105_C1504",
        "start": 1774513,
        "end": 1775634,
        "strand": "+",
        "product": "GDP-mannose 4,6-dehydratase",
    },
    {
        "source": "CP002246.1.embl",
        "marker_id": "YE_O9_per",
        "serotype": "O:9",
        "gene": "per",
        "locus_tag": "YE105_C1505",
        "start": 1775650,
        "end": 1776735,
        "strand": "+",
        "product": "Perosamine synthetase, Per protein",
    },
    {
        "source": "CP002246.1.embl",
        "marker_id": "YE_O9_wzm_like",
        "serotype": "O:9",
        "gene": "wzm_like",
        "locus_tag": "YE105_C1506",
        "start": 1776791,
        "end": 1777573,
        "strand": "+",
        "product": "lipopolysaccharide transport system permease protein",
    },
]


O527_MARKERS = [
    {
        "marker_id": "YE_O527_gt2",
        "serotype": "O:5,27",
        "gene": "gt2",
        "seqid": "NZ_CIFH01000063.1",
        "locus_tag": "AP310_RS19300",
        "start": 11392,
        "end": 12153,
        "strand": "-",
        "product": "glycosyltransferase family 2 protein",
    },
    {
        "marker_id": "YE_O527_rfbC",
        "serotype": "O:5,27",
        "gene": "rfbC",
        "seqid": "NZ_CIFH01000063.1",
        "locus_tag": "AP310_RS19305",
        "start": 12187,
        "end": 12735,
        "strand": "-",
        "product": "dTDP-4-dehydrorhamnose 3,5-epimerase",
    },
    {
        "marker_id": "YE_O527_rfbD",
        "serotype": "O:5,27",
        "gene": "rfbD",
        "seqid": "NZ_CIFH01000063.1",
        "locus_tag": "AP310_RS19310",
        "start": 12740,
        "end": 13618,
        "strand": "-",
        "product": "dTDP-4-dehydrorhamnose reductase",
    },
]


def _extract_record_map(path: Path, fmt: str) -> dict[str, object]:
    return {record.id: record for record in SeqIO.parse(path, fmt)}


def _extract_fasta_map(path: Path) -> dict[str, object]:
    return {record.id: record for record in SeqIO.parse(path, "fasta")}


def _slice_sequence(record, start: int, end: int, strand: str) -> str:
    seq = record.seq[start - 1 : end]
    if strand == "-":
        seq = seq.reverse_complement()
    return str(seq)


def main() -> int:
    embl_cache: dict[str, dict[str, object]] = {}
    fasta_records = _extract_fasta_map(TMPDB_DIR / "o527_genomic.fna")
    out_fasta = DB_DIR / "reference_o_antigen_markers.fasta"
    lines: list[str] = []

    for marker in EMBL_MARKERS:
        source = marker["source"]
        if source not in embl_cache:
            embl_cache[source] = _extract_record_map(TMPDB_DIR / source, "embl")
        record_id = source.replace(".embl", "")
        record = embl_cache[source][record_id]
        seq = _slice_sequence(record, marker["start"], marker["end"], marker["strand"])
        header = (
            f">{marker['marker_id']} serotype={marker['serotype']} gene={marker['gene']} "
            f"locus_tag={marker['locus_tag']} source={record_id} "
            f"coords={marker['start']}-{marker['end']} strand={marker['strand']} "
            f"product=\"{marker['product']}\""
        )
        lines.extend([header, seq])

    for marker in O527_MARKERS:
        record = fasta_records[marker["seqid"]]
        seq = _slice_sequence(record, marker["start"], marker["end"], marker["strand"])
        header = (
            f">{marker['marker_id']} serotype={marker['serotype']} gene={marker['gene']} "
            f"locus_tag={marker['locus_tag']} source={marker['seqid']} "
            f"coords={marker['start']}-{marker['end']} strand={marker['strand']} "
            f"product=\"{marker['product']}\""
        )
        lines.extend([header, seq])

    out_fasta.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(out_fasta)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
