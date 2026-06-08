#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import subprocess
from collections import defaultdict
from pathlib import Path


ROOT = Path("/Users/wuhhh/Desktop/徐老师/代码/metagenomic")
DEFAULT_REF_DIR = ROOT / "database" / "bacteria" / "yersinia_enterocolitica" / "o_antigen"
DEFAULT_REF_FASTA = DEFAULT_REF_DIR / "reference_o_antigen_markers.fasta"
DEFAULT_REF_MANIFEST = DEFAULT_REF_DIR / "reference_o_antigen_markers.tsv"
DEFAULT_OUT_DIR = ROOT / "tmp" / "yersinia_enterocolitica_o_antigen_typing"
DEFAULT_BLASTN = ROOT / "soft" / "ncbi-blast" / "bin" / "blastn"

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


def read_fasta_lengths(path: Path) -> tuple[dict[str, int], int]:
    record_lengths: dict[str, int] = {}
    current_id = ""
    current_length = 0
    total_length = 0
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id:
                    record_lengths[current_id] = current_length
                    total_length += current_length
                current_id = line[1:].split()[0]
                current_length = 0
                continue
            current_length += len(line)
        if current_id:
            record_lengths[current_id] = current_length
            total_length += current_length
    return record_lengths, total_length


def parse_bool(value: object) -> bool:
    normalized = str(value or "").strip().lower()
    return normalized in {"1", "true", "yes", "y", "required", "core"}


def load_manifest(path: Path) -> tuple[dict[str, dict[str, str]], dict[str, list[dict[str, str]]]]:
    marker_meta: dict[str, dict[str, str]] = {}
    serotype_markers: dict[str, list[dict[str, str]]] = defaultdict(list)
    with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required_columns = {"marker_id", "serotype", "gene"}
        if not reader.fieldnames or not required_columns.issubset(set(reader.fieldnames)):
            raise ValueError(
                f"Reference manifest must contain columns {sorted(required_columns)}: {path}"
            )
        for row in reader:
            marker_id = str(row.get("marker_id") or "").strip()
            serotype = str(row.get("serotype") or "").strip()
            gene = str(row.get("gene") or "").strip()
            if not marker_id or not serotype or not gene:
                continue
            normalized = {
                "marker_id": marker_id,
                "serotype": serotype,
                "gene": gene,
                "role": str(row.get("role") or "").strip() or ("required" if parse_bool(row.get("required")) else "accessory"),
                "required": "1" if parse_bool(row.get("required")) or str(row.get("role") or "").strip().lower() == "required" else "0",
                "group_id": str(row.get("group_id") or "").strip(),
                "note": str(row.get("note") or "").strip(),
            }
            marker_meta[marker_id] = normalized
            serotype_markers[serotype].append(normalized)
    if not marker_meta:
        raise ValueError(f"No valid marker rows found in manifest: {path}")
    return marker_meta, serotype_markers


def run_blastn(
    blastn_bin: str,
    query_fasta: Path,
    ref_fasta: Path,
    out_path: Path,
    threads: int,
) -> None:
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
        "200",
        "-dust",
        "no",
        "-num_threads",
        str(max(1, int(threads or 1))),
        "-out",
        str(out_path),
    ]
    subprocess.run(cmd, check=True)


def parse_blast_hits(
    path: Path,
    marker_meta: dict[str, dict[str, str]],
    min_hit_pident: float,
    min_hit_coverage: float,
) -> list[dict[str, str]]:
    hits: list[dict[str, str]] = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != len(OUTFMT_FIELDS):
                continue
            row = dict(zip(OUTFMT_FIELDS, parts))
            marker_id = str(row["sseqid"] or "").strip().split()[0]
            meta = marker_meta.get(marker_id)
            if meta is None:
                continue
            try:
                pident = float(row["pident"])
                align_len = int(float(row["length"]))
                ref_len = int(float(row["slen"]))
                bitscore = float(row["bitscore"])
            except (TypeError, ValueError):
                continue
            if ref_len <= 0:
                continue
            coverage = align_len / ref_len
            if pident < min_hit_pident or coverage < min_hit_coverage:
                continue
            hits.append(
                {
                    "query_id": row["qseqid"],
                    "marker_id": marker_id,
                    "serotype": meta["serotype"],
                    "gene": meta["gene"],
                    "role": meta["role"],
                    "required": meta["required"],
                    "group_id": meta["group_id"],
                    "pident": f"{pident:.2f}",
                    "align_length": str(align_len),
                    "reference_length": str(ref_len),
                    "coverage": f"{coverage:.4f}",
                    "bitscore": f"{bitscore:.1f}",
                    "evalue": str(row["evalue"]),
                    "qstart": row["qstart"],
                    "qend": row["qend"],
                    "sstart": row["sstart"],
                    "send": row["send"],
                    "note": meta["note"],
                }
            )
    return hits


def best_hit_per_marker(hits: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    best: dict[str, dict[str, str]] = {}
    for hit in hits:
        marker_id = hit["marker_id"]
        score = (
            float(hit["bitscore"]),
            float(hit["coverage"]),
            float(hit["pident"]),
            int(hit["align_length"]),
        )
        previous = best.get(marker_id)
        if previous is None:
            hit["_score"] = score
            best[marker_id] = hit
            continue
        previous_score = previous.get("_score") or (0.0, 0.0, 0.0, 0)
        if score > previous_score:
            hit["_score"] = score
            best[marker_id] = hit
    return best


def build_candidate_rows(
    sample_id: str,
    serotype_markers: dict[str, list[dict[str, str]]],
    marker_best_hits: dict[str, dict[str, str]],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for serotype, markers in sorted(serotype_markers.items()):
        expected_total = len(markers)
        required_markers = [marker for marker in markers if marker["required"] == "1"]
        required_total = len(required_markers)
        detected_hits: list[dict[str, str]] = []
        detected_required_hits: list[dict[str, str]] = []
        missing_required: list[str] = []
        for marker in markers:
            marker_id = marker["marker_id"]
            best_hit = marker_best_hits.get(marker_id)
            if best_hit is not None:
                detected_hits.append(best_hit)
                if marker["required"] == "1":
                    detected_required_hits.append(best_hit)
            elif marker["required"] == "1":
                missing_required.append(marker["gene"])
        detected_total = len(detected_hits)
        required_detected = len(detected_required_hits)
        expected_fraction = (detected_total / expected_total) if expected_total else 0.0
        required_fraction = (required_detected / required_total) if required_total else 0.0
        mean_pident = sum(float(hit["pident"]) for hit in detected_hits) / detected_total if detected_hits else 0.0
        mean_coverage = sum(float(hit["coverage"]) for hit in detected_hits) / detected_total if detected_hits else 0.0
        mean_bitscore = sum(float(hit["bitscore"]) for hit in detected_hits) / detected_total if detected_hits else 0.0
        # Favor the amount of independent locus support over raw alignment
        # length-derived bitscore averages. This avoids shorter serotype panels
        # outranking richer clusters solely because one marker is longer.
        score = (
            required_detected * 1000.0
            + detected_total * 100.0
            + required_fraction * 100.0
            + expected_fraction * 50.0
            + mean_pident
            + (mean_coverage * 10.0)
        )
        rows.append(
            {
                "sample_id": sample_id,
                "serotype": serotype,
                "expected_markers": str(expected_total),
                "required_markers": str(required_total),
                "detected_markers": str(detected_total),
                "detected_required_markers": str(required_detected),
                "expected_fraction": f"{expected_fraction:.4f}",
                "required_fraction": f"{required_fraction:.4f}",
                "mean_pident": f"{mean_pident:.2f}",
                "mean_coverage": f"{mean_coverage:.4f}",
                "mean_bitscore": f"{mean_bitscore:.1f}",
                "score": f"{score:.4f}",
                "detected_genes": ";".join(sorted({hit["gene"] for hit in detected_hits})),
                "missing_required_genes": ";".join(missing_required),
            }
        )
    rows.sort(
        key=lambda row: (
            float(row["score"]),
            float(row["required_fraction"]),
            float(row["expected_fraction"]),
            float(row["mean_pident"]),
            float(row["mean_coverage"]),
        ),
        reverse=True,
    )
    return rows


def classify_top_candidate(
    candidates: list[dict[str, str]],
    min_required_fraction: float,
    min_expected_fraction: float,
    min_pident: float,
    min_coverage: float,
    suspected_min_required_fraction: float,
    suspected_min_expected_fraction: float,
    suspected_min_pident: float,
    suspected_min_coverage: float,
    min_score_delta: float,
) -> tuple[str, str]:
    if not candidates:
        return "unassigned", "未找到任何 O 抗候选命中"
    top = candidates[0]
    second = candidates[1] if len(candidates) > 1 else None
    top_score = float(top["score"])
    second_score = float(second["score"]) if second is not None else 0.0
    score_delta = top_score - second_score
    top_required_fraction = float(top["required_fraction"])
    top_expected_fraction = float(top["expected_fraction"])
    top_pident = float(top["mean_pident"])
    top_coverage = float(top["mean_coverage"])

    if (
        top_required_fraction >= min_required_fraction
        and top_expected_fraction >= min_expected_fraction
        and top_pident >= min_pident
        and top_coverage >= min_coverage
        and (second is None or score_delta >= min_score_delta or top_required_fraction > float(second["required_fraction"]))
    ):
        return "typed", "满足高置信 O 抗分型阈值"

    if (
        top_required_fraction >= suspected_min_required_fraction
        and top_expected_fraction >= suspected_min_expected_fraction
        and top_pident >= suspected_min_pident
        and top_coverage >= suspected_min_coverage
    ):
        if second is not None and score_delta < min_score_delta:
            return "ambiguous", "首位与次位候选得分接近，建议人工复核"
        return "review", "达到低一级证据阈值，建议结合参考库质量和位点完整性复核"

    return "weak", "命中不足以支持可靠 O 抗分型"


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Type Yersinia enterocolitica O-antigen from a local marker reference panel"
    )
    parser.add_argument("--input-fasta", "--input", dest="input_fasta", type=Path, required=True)
    parser.add_argument("--sample-id", type=str, default="", help="Optional sample name; defaults to input FASTA stem")
    parser.add_argument("--reference-fasta", type=Path, default=DEFAULT_REF_FASTA)
    parser.add_argument("--reference-manifest", type=Path, default=DEFAULT_REF_MANIFEST)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--blastn-bin", type=str, default="")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--min-hit-pident", type=float, default=70.0)
    parser.add_argument("--min-hit-coverage", type=float, default=0.50)
    parser.add_argument("--min-pident", type=float, default=90.0)
    parser.add_argument("--min-coverage", type=float, default=0.85)
    parser.add_argument("--min-required-fraction", type=float, default=1.0)
    parser.add_argument("--min-expected-fraction", type=float, default=0.60)
    parser.add_argument("--suspected-min-pident", type=float, default=80.0)
    parser.add_argument("--suspected-min-coverage", type=float, default=0.60)
    parser.add_argument("--suspected-min-required-fraction", type=float, default=0.50)
    parser.add_argument("--suspected-min-expected-fraction", type=float, default=0.40)
    parser.add_argument("--min-score-delta", type=float, default=25.0)
    args = parser.parse_args()

    input_fasta = args.input_fasta.resolve()
    ref_fasta = args.reference_fasta.resolve()
    ref_manifest = args.reference_manifest.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if not input_fasta.is_file():
        raise FileNotFoundError(f"Input FASTA not found: {input_fasta}")
    if not ref_fasta.is_file():
        raise FileNotFoundError(
            f"Reference FASTA not found: {ref_fasta}. Please populate the Yersinia O-antigen reference panel first."
        )
    if not ref_manifest.is_file():
        raise FileNotFoundError(
            f"Reference manifest not found: {ref_manifest}. Please populate the Yersinia O-antigen reference panel first."
        )

    sample_id = str(args.sample_id or input_fasta.stem).strip() or input_fasta.stem
    contig_lengths, assembly_length = read_fasta_lengths(input_fasta)
    if not contig_lengths:
        raise RuntimeError(f"No FASTA records found in {input_fasta}")

    marker_meta, serotype_markers = load_manifest(ref_manifest)
    blast_out = output_dir / f"{sample_id}.o_antigen.blastn.tsv"
    marker_hits_tsv = output_dir / f"{sample_id}.o_antigen.marker_hits.tsv"
    candidates_tsv = output_dir / f"{sample_id}.o_antigen.candidates.tsv"
    summary_tsv = output_dir / f"{sample_id}.o_antigen.summary.tsv"
    summary_json = output_dir / f"{sample_id}.o_antigen.summary.json"

    run_blastn(
        blastn_bin=resolve_blastn(args.blastn_bin),
        query_fasta=input_fasta,
        ref_fasta=ref_fasta,
        out_path=blast_out,
        threads=args.threads,
    )
    raw_hits = parse_blast_hits(
        blast_out,
        marker_meta,
        min_hit_pident=args.min_hit_pident,
        min_hit_coverage=args.min_hit_coverage,
    )
    marker_best_hits = best_hit_per_marker(raw_hits)
    best_hit_rows = sorted(
        marker_best_hits.values(),
        key=lambda row: (
            row["serotype"],
            row["gene"],
            -float(row["bitscore"]),
        ),
    )
    candidate_rows = build_candidate_rows(sample_id, serotype_markers, marker_best_hits)
    top_candidate = candidate_rows[0] if candidate_rows else {}
    status, note = classify_top_candidate(
        candidate_rows,
        min_required_fraction=args.min_required_fraction,
        min_expected_fraction=args.min_expected_fraction,
        min_pident=args.min_pident,
        min_coverage=args.min_coverage,
        suspected_min_required_fraction=args.suspected_min_required_fraction,
        suspected_min_expected_fraction=args.suspected_min_expected_fraction,
        suspected_min_pident=args.suspected_min_pident,
        suspected_min_coverage=args.suspected_min_coverage,
        min_score_delta=args.min_score_delta,
    )
    second_score = float(candidate_rows[1]["score"]) if len(candidate_rows) > 1 else 0.0
    top_score = float(top_candidate.get("score", "0") or 0.0)
    summary_row = {
        "sample_id": sample_id,
        "input_fasta": str(input_fasta),
        "assembly_length": str(assembly_length),
        "contig_count": str(len(contig_lengths)),
        "predicted_o_serotype": str(top_candidate.get("serotype") or ""),
        "status": status,
        "note": note,
        "score": f"{top_score:.4f}",
        "score_delta_vs_next": f"{(top_score - second_score):.4f}",
        "expected_markers": str(top_candidate.get("expected_markers") or "0"),
        "required_markers": str(top_candidate.get("required_markers") or "0"),
        "detected_markers": str(top_candidate.get("detected_markers") or "0"),
        "detected_required_markers": str(top_candidate.get("detected_required_markers") or "0"),
        "expected_fraction": str(top_candidate.get("expected_fraction") or "0"),
        "required_fraction": str(top_candidate.get("required_fraction") or "0"),
        "mean_pident": str(top_candidate.get("mean_pident") or "0"),
        "mean_coverage": str(top_candidate.get("mean_coverage") or "0"),
        "detected_genes": str(top_candidate.get("detected_genes") or ""),
        "missing_required_genes": str(top_candidate.get("missing_required_genes") or ""),
        "reference_fasta": str(ref_fasta),
        "reference_manifest": str(ref_manifest),
    }

    write_tsv(
        marker_hits_tsv,
        best_hit_rows,
        [
            "query_id",
            "marker_id",
            "serotype",
            "gene",
            "role",
            "required",
            "group_id",
            "pident",
            "align_length",
            "reference_length",
            "coverage",
            "bitscore",
            "evalue",
            "qstart",
            "qend",
            "sstart",
            "send",
            "note",
        ],
    )
    write_tsv(
        candidates_tsv,
        candidate_rows,
        [
            "sample_id",
            "serotype",
            "expected_markers",
            "required_markers",
            "detected_markers",
            "detected_required_markers",
            "expected_fraction",
            "required_fraction",
            "mean_pident",
            "mean_coverage",
            "mean_bitscore",
            "score",
            "detected_genes",
            "missing_required_genes",
        ],
    )
    write_tsv(
        summary_tsv,
        [summary_row],
        [
            "sample_id",
            "input_fasta",
            "assembly_length",
            "contig_count",
            "predicted_o_serotype",
            "status",
            "note",
            "score",
            "score_delta_vs_next",
            "expected_markers",
            "required_markers",
            "detected_markers",
            "detected_required_markers",
            "expected_fraction",
            "required_fraction",
            "mean_pident",
            "mean_coverage",
            "detected_genes",
            "missing_required_genes",
            "reference_fasta",
            "reference_manifest",
        ],
    )
    summary_json.write_text(
        json.dumps(
            {
                "sample_id": sample_id,
                "summary": summary_row,
                "top_candidates": candidate_rows[:10],
                "marker_hits": [
                    {key: value for key, value in row.items() if not key.startswith("_")}
                    for row in best_hit_rows
                ],
            },
            ensure_ascii=False,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    print(f"Sample: {sample_id}")
    print(f"Predicted O serotype: {summary_row['predicted_o_serotype'] or '-'}")
    print(f"Status: {status}")
    print(f"Summary TSV: {summary_tsv}")
    print(f"Candidate TSV: {candidates_tsv}")
    print(f"Marker hits TSV: {marker_hits_tsv}")
    print(f"Summary JSON: {summary_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
