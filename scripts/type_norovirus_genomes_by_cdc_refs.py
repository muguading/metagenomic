from __future__ import annotations

import argparse
import csv
import re
import subprocess
from collections import defaultdict
from pathlib import Path


GENE_CONFIGS = {
    "rdrp": {
        "ref_fasta": Path("database/virus/norovirus/cdc_typing_refs/cdc_norovirus_rdrp_refs.fasta"),
        "manifest": Path("database/virus/norovirus/cdc_typing_refs/cdc_norovirus_rdrp_refs.tsv"),
    },
    "vp1": {
        "ref_fasta": Path("database/virus/norovirus/cdc_typing_refs/cdc_norovirus_vp1_refs.fasta"),
        "manifest": Path("database/virus/norovirus/cdc_typing_refs/cdc_norovirus_vp1_refs.tsv"),
    },
}

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
]


def _read_fasta(path: Path) -> dict[str, tuple[str, str]]:
    records: dict[str, tuple[str, str]] = {}
    header = ""
    seq_lines: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header:
                    seq_id = header.split()[0]
                    records[seq_id] = (header, "".join(seq_lines))
                header = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line.strip())
    if header:
        seq_id = header.split()[0]
        records[seq_id] = (header, "".join(seq_lines))
    return records


def _read_manifest(path: Path, gene_name: str) -> dict[str, dict[str, str]]:
    mapping: dict[str, dict[str, str]] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            accession = str(row.get("accession") or "").strip()
            subtype = str(row.get("subtype") or "").strip()
            if not accession or not subtype:
                continue
            subject_id = f"{accession}_{subtype}_{gene_name.upper()}"
            mapping[subject_id] = row
            mapping[accession.split(".", 1)[0]] = row
    return mapping


def _count_subjects_per_subtype(path: Path) -> dict[str, int]:
    counts: dict[str, int] = defaultdict(int)
    seen: set[tuple[str, str]] = set()
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            subtype = _normalize_type(str(row.get("subtype") or ""))
            accession = str(row.get("accession") or "").strip()
            if not subtype or not accession:
                continue
            key = (subtype, accession.split(".", 1)[0])
            if key in seen:
                continue
            seen.add(key)
            counts[subtype] += 1
    return counts


def _normalize_type(value: str) -> str:
    return str(value or "").strip().upper().replace("_", ".")


def _extract_genogroup(type_label: str) -> str:
    matched = re.match(r"^(G[IVX]+)", _normalize_type(type_label))
    return matched.group(1) if matched else ""


def _build_blast_db(makeblastdb_bin: Path, fasta_path: Path, db_prefix: Path) -> None:
    expected = db_prefix.with_suffix(".nhr")
    if expected.exists():
        return
    db_prefix.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(makeblastdb_bin),
        "-in",
        str(fasta_path),
        "-dbtype",
        "nucl",
        "-out",
        str(db_prefix),
    ]
    subprocess.run(cmd, check=True)


def _run_blastn(blastn_bin: Path, query_fasta: Path, db_prefix: Path, out_path: Path, threads: int) -> None:
    cmd = [
        str(blastn_bin),
        "-task",
        "blastn",
        "-query",
        str(query_fasta),
        "-db",
        str(db_prefix),
        "-outfmt",
        "6 " + " ".join(OUTFMT_FIELDS),
        "-evalue",
        "1e-20",
        "-max_target_seqs",
        "50",
        "-dust",
        "no",
        "-num_threads",
        str(max(1, threads)),
        "-out",
        str(out_path),
    ]
    subprocess.run(cmd, check=True)


def _summarize_hits_by_subtype(
    blast_tsv: Path,
    subject_meta: dict[str, dict[str, str]],
    subject_records: dict[str, tuple[str, str]],
    subtype_subject_totals: dict[str, int],
) -> dict[str, dict[str, dict[str, str]]]:
    summarized_hits: dict[str, dict[str, dict[str, object]]] = defaultdict(dict)
    with blast_tsv.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            values = line.split("\t")
            row = dict(zip(OUTFMT_FIELDS, values))
            qseqid = row["qseqid"]
            subject_id = row["sseqid"]
            meta = subject_meta.get(subject_id) or subject_meta.get(subject_id.split("_", 1)[0]) or {}
            subtype = _normalize_type(meta.get("subtype") or "")
            if not subtype:
                continue
            ref_sequence = subject_records.get(subject_id, ("", ""))[1]
            ref_length = len(ref_sequence) if ref_sequence else 0
            align_length = int(row["length"])
            qcov_ref = (align_length / ref_length) if ref_length else 0.0
            candidate = {
                "subtype": subtype,
                "genogroup": _extract_genogroup(subtype),
                "subject_id": subject_id,
                "accession": str(meta.get("accession") or ""),
                "label": str(meta.get("label") or ""),
                "pident": row["pident"],
                "align_length": row["length"],
                "mismatch": row["mismatch"],
                "gapopen": row["gapopen"],
                "qstart": row["qstart"],
                "qend": row["qend"],
                "sstart": row["sstart"],
                "send": row["send"],
                "evalue": row["evalue"],
                "bitscore": row["bitscore"],
                "ref_length": str(ref_length),
                "qcov_ref": f"{qcov_ref:.4f}",
            }
            bucket = summarized_hits[qseqid].get(subtype)
            if bucket is None:
                summarized_hits[qseqid][subtype] = {
                    "best_hit": candidate,
                    "hit_count": 1,
                    "subject_ids": {subject_id},
                    "sum_pident": float(candidate["pident"]),
                    "sum_qcov_ref": float(candidate["qcov_ref"]),
                    "max_pident": float(candidate["pident"]),
                    "max_qcov_ref": float(candidate["qcov_ref"]),
                }
                continue
            previous = bucket["best_hit"]
            previous_score = (
                float(previous["bitscore"]),
                float(previous["pident"]),
                float(previous["qcov_ref"]),
                int(previous["align_length"]),
            )
            candidate_score = (
                float(candidate["bitscore"]),
                float(candidate["pident"]),
                float(candidate["qcov_ref"]),
                int(candidate["align_length"]),
            )
            if candidate_score > previous_score:
                bucket["best_hit"] = candidate
            bucket["hit_count"] = int(bucket["hit_count"]) + 1
            cast_subjects = bucket["subject_ids"]
            assert isinstance(cast_subjects, set)
            cast_subjects.add(subject_id)
            bucket["sum_pident"] = float(bucket["sum_pident"]) + float(candidate["pident"])
            bucket["sum_qcov_ref"] = float(bucket["sum_qcov_ref"]) + float(candidate["qcov_ref"])
            bucket["max_pident"] = max(float(bucket["max_pident"]), float(candidate["pident"]))
            bucket["max_qcov_ref"] = max(float(bucket["max_qcov_ref"]), float(candidate["qcov_ref"]))

    normalized: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)
    for qseqid, subtype_buckets in summarized_hits.items():
        for subtype, bucket in subtype_buckets.items():
            best_hit = dict(bucket["best_hit"])
            hit_count = int(bucket["hit_count"])
            subject_count = len(bucket["subject_ids"])
            subject_count_total = int(subtype_subject_totals.get(subtype, 0))
            subtype_hit_fraction = (subject_count / subject_count_total) if subject_count_total else 0.0
            best_hit.update(
                {
                    "hit_count": str(hit_count),
                    "subject_count_hit": str(subject_count),
                    "subject_count_total": str(subject_count_total),
                    "mean_pident": f"{float(bucket['sum_pident']) / hit_count:.3f}",
                    "mean_qcov_ref": f"{float(bucket['sum_qcov_ref']) / hit_count:.4f}",
                    "max_pident": f"{float(bucket['max_pident']):.3f}",
                    "max_qcov_ref": f"{float(bucket['max_qcov_ref']):.4f}",
                    "subtype_hit_fraction": f"{subtype_hit_fraction:.4f}",
                }
            )
            normalized[qseqid][subtype] = best_hit
    return normalized


def _pick_top_hit(subtype_hits: dict[str, dict[str, str]]) -> dict[str, str]:
    if not subtype_hits:
        return {}
    ranked = sorted(
        subtype_hits.values(),
        key=lambda item: (
            float(item["bitscore"]),
            float(item["pident"]),
            float(item["qcov_ref"]),
            float(item.get("subtype_hit_fraction", "0")),
            int(item.get("subject_count_hit", "0")),
            int(item["align_length"]),
        ),
        reverse=True,
    )
    return ranked[0]


def _top2_distance(subtype_hits: dict[str, dict[str, str]]) -> tuple[str, str]:
    ranked = sorted(
        subtype_hits.values(),
        key=lambda item: (
            float(item["bitscore"]),
            float(item["pident"]),
            float(item["qcov_ref"]),
            float(item.get("subtype_hit_fraction", "0")),
            int(item["align_length"]),
        ),
        reverse=True,
    )
    if len(ranked) < 2:
        return "", ""
    bitscore_delta = float(ranked[0]["bitscore"]) - float(ranked[1]["bitscore"])
    pident_delta = float(ranked[0]["pident"]) - float(ranked[1]["pident"])
    qcov_delta = float(ranked[0]["qcov_ref"]) - float(ranked[1]["qcov_ref"])
    support_delta = float(ranked[0].get("subtype_hit_fraction", "0")) - float(ranked[1].get("subtype_hit_fraction", "0"))
    return f"{bitscore_delta:.2f}", f"{pident_delta:.2f}", f"{qcov_delta:.4f}", f"{support_delta:.4f}"


def _evaluate_confident_hit(
    subtype_hits: dict[str, dict[str, str]],
    min_pident: float,
    suspected_min_pident: float,
    min_qcov_ref: float,
    suspected_min_qcov_ref: float,
    min_subtype_hit_fraction: float,
    min_subject_count_hit: int,
    min_bitscore_delta: float,
    min_pident_delta: float,
) -> tuple[dict[str, str], str]:
    top_hit = _pick_top_hit(subtype_hits)
    if not top_hit:
        return {}, "no_hit"
    top_pident = float(top_hit["pident"])
    top_qcov_ref = float(top_hit["qcov_ref"])
    if top_pident < suspected_min_pident:
        return {}, "low_identity"
    if top_qcov_ref < suspected_min_qcov_ref:
        return {}, "low_coverage"
    suspected_reasons: list[str] = []
    if top_pident < min_pident:
        suspected_reasons.append("identity")
    if top_qcov_ref < min_qcov_ref:
        suspected_reasons.append("coverage")
    if suspected_reasons:
        return top_hit, f"suspected_{'_and_'.join(suspected_reasons)}"
    subject_count_total = int(top_hit.get("subject_count_total", "0"))
    subject_count_hit = int(top_hit.get("subject_count_hit", "0"))
    subtype_hit_fraction = float(top_hit.get("subtype_hit_fraction", "0"))
    if subject_count_total <= 2:
        if subject_count_hit < 1:
            return {}, "low_consistency"
    else:
        if subtype_hit_fraction < min_subtype_hit_fraction and subject_count_hit < min_subject_count_hit:
            return {}, "low_consistency"
    ranked = sorted(
        subtype_hits.values(),
        key=lambda item: (
            float(item["bitscore"]),
            float(item["pident"]),
            float(item["qcov_ref"]),
            float(item.get("subtype_hit_fraction", "0")),
            int(item.get("subject_count_hit", "0")),
            int(item["align_length"]),
        ),
        reverse=True,
    )
    if len(ranked) >= 2:
        bitscore_delta = float(ranked[0]["bitscore"]) - float(ranked[1]["bitscore"])
        pident_delta = float(ranked[0]["pident"]) - float(ranked[1]["pident"])
        same_genogroup = ranked[0].get("genogroup") and ranked[0].get("genogroup") == ranked[1].get("genogroup")
        if same_genogroup and bitscore_delta < min_bitscore_delta and pident_delta < min_pident_delta:
            return {}, "ambiguous_top2"
    return top_hit, "accepted"


def _write_gene_detail(
    out_path: Path,
    query_records: dict[str, tuple[str, str]],
    best_hits_by_query: dict[str, dict[str, dict[str, str]]],
    gene_name: str,
) -> None:
    columns = [
        "query_id",
        "query_length",
        f"{gene_name}_subtype",
        f"{gene_name}_genogroup",
        "subject_id",
        "accession",
        "label",
        "pident",
        "align_length",
        "ref_length",
        "qcov_ref",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "hit_count",
        "subject_count_hit",
        "subject_count_total",
        "mean_pident",
        "mean_qcov_ref",
        "subtype_hit_fraction",
        "top2_bitscore_delta",
        "top2_pident_delta",
        "top2_qcov_delta",
        "top2_support_delta",
    ]
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for query_id in sorted(query_records):
            top_hit = _pick_top_hit(best_hits_by_query.get(query_id, {}))
            bit_delta, pid_delta, qcov_delta, support_delta = _top2_distance(best_hits_by_query.get(query_id, {}))
            writer.writerow(
                {
                    "query_id": query_id,
                    "query_length": len(query_records[query_id][1]),
                    f"{gene_name}_subtype": top_hit.get("subtype", ""),
                    f"{gene_name}_genogroup": top_hit.get("genogroup", ""),
                    "subject_id": top_hit.get("subject_id", ""),
                    "accession": top_hit.get("accession", ""),
                    "label": top_hit.get("label", ""),
                    "pident": top_hit.get("pident", ""),
                    "align_length": top_hit.get("align_length", ""),
                    "ref_length": top_hit.get("ref_length", ""),
                    "qcov_ref": top_hit.get("qcov_ref", ""),
                    "qstart": top_hit.get("qstart", ""),
                    "qend": top_hit.get("qend", ""),
                    "sstart": top_hit.get("sstart", ""),
                    "send": top_hit.get("send", ""),
                    "evalue": top_hit.get("evalue", ""),
                    "bitscore": top_hit.get("bitscore", ""),
                    "hit_count": top_hit.get("hit_count", ""),
                    "subject_count_hit": top_hit.get("subject_count_hit", ""),
                    "subject_count_total": top_hit.get("subject_count_total", ""),
                    "mean_pident": top_hit.get("mean_pident", ""),
                    "mean_qcov_ref": top_hit.get("mean_qcov_ref", ""),
                    "subtype_hit_fraction": top_hit.get("subtype_hit_fraction", ""),
                    "top2_bitscore_delta": bit_delta,
                    "top2_pident_delta": pid_delta,
                    "top2_qcov_delta": qcov_delta,
                    "top2_support_delta": support_delta,
                }
            )


def _compose_dual_type(rdrp_hit: dict[str, str], vp1_hit: dict[str, str]) -> str:
    rdrp_type = rdrp_hit.get("subtype", "")
    vp1_type = vp1_hit.get("subtype", "")
    if rdrp_type and vp1_type:
        return f"{rdrp_type}_{vp1_type}"
    return rdrp_type or vp1_type or ""


def _compose_status(
    rdrp_hit: dict[str, str],
    vp1_hit: dict[str, str],
    rdrp_reason: str,
    vp1_reason: str,
) -> tuple[str, str]:
    suspected_reasons = {
        "suspected_identity",
        "suspected_coverage",
        "suspected_identity_and_coverage",
    }
    if rdrp_hit and vp1_hit:
        if rdrp_hit.get("genogroup") and vp1_hit.get("genogroup") and rdrp_hit["genogroup"] != vp1_hit["genogroup"]:
            return "discordant", "RdRp and VP1 genogroups are discordant"
        if rdrp_reason in suspected_reasons or vp1_reason in suspected_reasons:
            return "suspected", "At least one typing call is a suspected assignment (identity 60-85% and/or coverage 60-90%)"
        return "dual_typed", "Both RdRp and VP1 typing calls are available"
    if rdrp_hit:
        if rdrp_reason in suspected_reasons:
            return "suspected", "Only RdRp provides a suspected assignment (identity 60-85% and/or coverage 60-90%)"
        return "rdrp_only", "Only RdRp typing call is available"
    if vp1_hit:
        if vp1_reason in suspected_reasons:
            return "suspected", "Only VP1 provides a suspected assignment (identity 60-85% and/or coverage 60-90%)"
        return "vp1_only", "Only VP1 typing call is available"
    return "unassigned", "No confident CDC reference hit found"


def main() -> int:
    parser = argparse.ArgumentParser(description="Type norovirus genomes using CDC RdRp and VP1 reference panels")
    parser.add_argument(
        "--query-fasta",
        type=Path,
        default=Path("database/virus/norovirus/full_genomes/human_norovirus_complete_genomes.fasta"),
    )
    parser.add_argument(
        "--blastn-bin",
        type=Path,
        default=Path("/opt/homebrew/Caskroom/mambaforge/base/envs/ncov/bin/blastn"),
    )
    parser.add_argument(
        "--makeblastdb-bin",
        type=Path,
        default=Path("/opt/homebrew/Caskroom/mambaforge/base/envs/ncov/bin/makeblastdb"),
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=8,
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("database/virus/norovirus/typing_by_cdc_refs"),
    )
    parser.add_argument("--min-pident", type=float, default=85.0)
    parser.add_argument("--suspected-min-pident", type=float, default=60.0)
    parser.add_argument("--min-qcov-ref", type=float, default=0.90)
    parser.add_argument("--suspected-min-qcov-ref", type=float, default=0.60)
    parser.add_argument("--min-subtype-hit-fraction", type=float, default=0.30)
    parser.add_argument("--min-subject-count-hit", type=int, default=2)
    parser.add_argument("--min-bitscore-delta", type=float, default=20.0)
    parser.add_argument("--min-pident-delta", type=float, default=0.5)
    args = parser.parse_args()

    query_fasta = args.query_fasta.resolve()
    blastn_bin = args.blastn_bin.resolve()
    makeblastdb_bin = args.makeblastdb_bin.resolve()
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    db_dir = out_dir / "blastdb"
    db_dir.mkdir(parents=True, exist_ok=True)

    query_records = _read_fasta(query_fasta)
    gene_hits: dict[str, dict[str, dict[str, dict[str, str]]]] = {}

    for gene_name, config in GENE_CONFIGS.items():
        ref_fasta = config["ref_fasta"].resolve()
        ref_manifest = config["manifest"].resolve()
        ref_records = _read_fasta(ref_fasta)
        ref_meta = _read_manifest(ref_manifest, gene_name)
        subtype_subject_totals = _count_subjects_per_subtype(ref_manifest)
        db_prefix = db_dir / f"cdc_norovirus_{gene_name}"
        _build_blast_db(makeblastdb_bin, ref_fasta, db_prefix)
        blast_out = out_dir / f"{gene_name}.blastn.tsv"
        _run_blastn(blastn_bin, query_fasta, db_prefix, blast_out, args.threads)
        best_hits_by_query = _summarize_hits_by_subtype(blast_out, ref_meta, ref_records, subtype_subject_totals)
        gene_hits[gene_name] = best_hits_by_query
        _write_gene_detail(out_dir / f"{gene_name}.typing.tsv", query_records, best_hits_by_query, gene_name)

    summary_columns = [
        "query_id",
        "query_length",
        "rdrp_type",
        "rdrp_genogroup",
        "rdrp_subject",
        "rdrp_pident",
        "rdrp_qcov_ref",
        "rdrp_subtype_hit_fraction",
        "rdrp_subject_count_hit",
        "rdrp_subject_count_total",
        "rdrp_accept_reason",
        "vp1_type",
        "vp1_genogroup",
        "vp1_subject",
        "vp1_pident",
        "vp1_qcov_ref",
        "vp1_subtype_hit_fraction",
        "vp1_subject_count_hit",
        "vp1_subject_count_total",
        "vp1_accept_reason",
        "dual_type",
        "typing_status",
        "typing_note",
    ]
    summary_path = out_dir / "norovirus_dual_typing.summary.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_columns, delimiter="\t")
        writer.writeheader()
        for query_id in sorted(query_records):
            rdrp_hit, rdrp_reason = _evaluate_confident_hit(
                gene_hits["rdrp"].get(query_id, {}),
                min_pident=args.min_pident,
                suspected_min_pident=args.suspected_min_pident,
                min_qcov_ref=args.min_qcov_ref,
                suspected_min_qcov_ref=args.suspected_min_qcov_ref,
                min_subtype_hit_fraction=args.min_subtype_hit_fraction,
                min_subject_count_hit=args.min_subject_count_hit,
                min_bitscore_delta=args.min_bitscore_delta,
                min_pident_delta=args.min_pident_delta,
            )
            vp1_hit, vp1_reason = _evaluate_confident_hit(
                gene_hits["vp1"].get(query_id, {}),
                min_pident=args.min_pident,
                suspected_min_pident=args.suspected_min_pident,
                min_qcov_ref=args.min_qcov_ref,
                suspected_min_qcov_ref=args.suspected_min_qcov_ref,
                min_subtype_hit_fraction=args.min_subtype_hit_fraction,
                min_subject_count_hit=args.min_subject_count_hit,
                min_bitscore_delta=args.min_bitscore_delta,
                min_pident_delta=args.min_pident_delta,
            )
            status, note = _compose_status(rdrp_hit, vp1_hit, rdrp_reason, vp1_reason)
            writer.writerow(
                {
                    "query_id": query_id,
                    "query_length": len(query_records[query_id][1]),
                    "rdrp_type": rdrp_hit.get("subtype", ""),
                    "rdrp_genogroup": rdrp_hit.get("genogroup", ""),
                    "rdrp_subject": rdrp_hit.get("subject_id", ""),
                    "rdrp_pident": rdrp_hit.get("pident", ""),
                    "rdrp_qcov_ref": rdrp_hit.get("qcov_ref", ""),
                    "rdrp_subtype_hit_fraction": rdrp_hit.get("subtype_hit_fraction", ""),
                    "rdrp_subject_count_hit": rdrp_hit.get("subject_count_hit", ""),
                    "rdrp_subject_count_total": rdrp_hit.get("subject_count_total", ""),
                    "rdrp_accept_reason": rdrp_reason,
                    "vp1_type": vp1_hit.get("subtype", ""),
                    "vp1_genogroup": vp1_hit.get("genogroup", ""),
                    "vp1_subject": vp1_hit.get("subject_id", ""),
                    "vp1_pident": vp1_hit.get("pident", ""),
                    "vp1_qcov_ref": vp1_hit.get("qcov_ref", ""),
                    "vp1_subtype_hit_fraction": vp1_hit.get("subtype_hit_fraction", ""),
                    "vp1_subject_count_hit": vp1_hit.get("subject_count_hit", ""),
                    "vp1_subject_count_total": vp1_hit.get("subject_count_total", ""),
                    "vp1_accept_reason": vp1_reason,
                    "dual_type": _compose_dual_type(rdrp_hit, vp1_hit),
                    "typing_status": status,
                    "typing_note": note,
                }
            )

    print(summary_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
