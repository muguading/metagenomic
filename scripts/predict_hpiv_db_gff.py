#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO


SUBTYPES = ("HPIV1", "HPIV2", "HPIV3", "HPIV4a", "HPIV4b")


@dataclass
class GffFeature:
    seqid: str
    source: str
    feature_type: str
    start: int
    end: int
    score: str
    strand: str
    phase: str
    attributes: str


@dataclass
class PafAlignment:
    query_name: str
    query_length: int
    query_start: int
    query_end: int
    strand: str
    target_name: str
    target_length: int
    target_start: int
    target_end: int
    matches: int
    block_length: int
    mapq: int
    cg: str


def parse_args() -> argparse.Namespace:
    project_root = Path(__file__).resolve().parents[1]
    default_db_root = project_root / "database" / "virus" / "hpiv"
    parser = argparse.ArgumentParser(
        description="根据 HPIV reference fasta + gff3，为各 subtype 数据库中的每条样本序列批量预测 GFF3。"
    )
    parser.add_argument("--db-root", type=Path, default=default_db_root, help="HPIV 数据库目录")
    parser.add_argument("--out-root", type=Path, default=default_db_root / "predicted_gff", help="输出目录")
    parser.add_argument("--summary", type=Path, default=default_db_root / "predicted_gff" / "hpiv_gff_transfer_summary.tsv", help="汇总 TSV 路径")
    parser.add_argument("--subtypes", nargs="*", default=list(SUBTYPES), help="要处理的 subtype 列表，如 HPIV1 HPIV3")
    parser.add_argument("--minimap2", type=str, default="", help="minimap2 路径，默认自动探测")
    parser.add_argument("--min-mapq", type=int, default=20, help="接受的最小 MAPQ")
    parser.add_argument("--min-query-coverage", type=float, default=0.7, help="最小 query 覆盖比例")
    parser.add_argument("--limit-per-subtype", type=int, default=0, help="每个 subtype 仅测试前 N 条序列，0 为全部")
    return parser.parse_args()


def resolve_minimap2(configured: str) -> str:
    candidates = [
        configured.strip() if configured else "",
        shutil.which("minimap2") or "",
        "/opt/homebrew/bin/minimap2",
        "/usr/local/bin/minimap2",
        str((Path(__file__).resolve().parents[1] / "soft" / "minimap2" / "minimap2").resolve()),
    ]
    for candidate in candidates:
        if candidate and Path(candidate).is_file():
            return str(Path(candidate).resolve())
    raise FileNotFoundError("未找到 minimap2，请通过 --minimap2 指定路径。")


def parse_gff(gff_path: Path) -> tuple[list[str], list[GffFeature]]:
    headers: list[str] = []
    features: list[GffFeature] = []
    with gff_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.strip():
                continue
            if line.startswith("#"):
                headers.append(line.rstrip("\n"))
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            try:
                start = int(parts[3])
                end = int(parts[4])
            except ValueError:
                continue
            features.append(
                GffFeature(
                    seqid=parts[0],
                    source=parts[1],
                    feature_type=parts[2],
                    start=start,
                    end=end,
                    score=parts[5],
                    strand=parts[6],
                    phase=parts[7],
                    attributes=parts[8],
                )
            )
    return headers, features


def parse_paf_line(line: str) -> PafAlignment | None:
    parts = line.rstrip("\n").split("\t")
    if len(parts) < 12:
        return None
    tags = {}
    for item in parts[12:]:
        fields = item.split(":", 2)
        if len(fields) == 3:
            tags[fields[0]] = fields[2]
    cg = tags.get("cg", "")
    if not cg:
        return None
    return PafAlignment(
        query_name=parts[0],
        query_length=int(parts[1]),
        query_start=int(parts[2]),
        query_end=int(parts[3]),
        strand=parts[4],
        target_name=parts[5],
        target_length=int(parts[6]),
        target_start=int(parts[7]),
        target_end=int(parts[8]),
        matches=int(parts[9]),
        block_length=int(parts[10]),
        mapq=int(parts[11]),
        cg=cg,
    )


def choose_best_alignment(alignments: list[PafAlignment], min_mapq: int, min_query_coverage: float) -> PafAlignment | None:
    filtered = [
        item for item in alignments
        if item.mapq >= min_mapq and item.query_length > 0 and ((item.query_end - item.query_start) / item.query_length) >= min_query_coverage
    ]
    if not filtered:
        return None
    return max(
        filtered,
        key=lambda item: (
            item.matches,
            item.block_length,
            item.mapq,
            (item.query_end - item.query_start),
        ),
    )


def parse_cigar(cigar: str) -> list[tuple[int, str]]:
    return [(int(length), op) for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)]


def map_ref_position(aln: PafAlignment, ref_pos_1based: int) -> int | None:
    ref_cursor = aln.target_start + 1
    query_consumed = 0
    for length, op in parse_cigar(aln.cg):
        if op in {"M", "=", "X"}:
            if ref_cursor <= ref_pos_1based < ref_cursor + length:
                delta = ref_pos_1based - ref_cursor
                query_offset = query_consumed + delta
                if aln.strand == "+":
                    return aln.query_start + 1 + query_offset
                return aln.query_end - query_offset
            ref_cursor += length
            query_consumed += length
            continue
        if op in {"D", "N"}:
            if ref_cursor <= ref_pos_1based < ref_cursor + length:
                return None
            ref_cursor += length
            continue
        if op in {"I", "S", "H", "P"}:
            query_consumed += length
            continue
    return None


def transfer_feature(feature: GffFeature, aln: PafAlignment, target_seqid: str) -> GffFeature | None:
    mapped_start = map_ref_position(aln, feature.start)
    mapped_end = map_ref_position(aln, feature.end)
    if mapped_start is None or mapped_end is None:
        return None
    new_start = min(mapped_start, mapped_end)
    new_end = max(mapped_start, mapped_end)
    return GffFeature(
        seqid=target_seqid,
        source=feature.source,
        feature_type=feature.feature_type,
        start=new_start,
        end=new_end,
        score=feature.score,
        strand=feature.strand,
        phase=feature.phase,
        attributes=feature.attributes,
    )


def sanitize_filename(name: str) -> str:
    text = str(name or "").strip()
    text = re.sub(r"[^\w.\-]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("_")
    return text or "sample"


def write_gff(
    output_path: Path,
    target_seqid: str,
    target_length: int,
    header_lines: list[str],
    features: list[GffFeature],
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write("##gff-version 3\n")
        for header in header_lines:
            if header.startswith("##gff-version"):
                continue
            if header.startswith("##sequence-region"):
                continue
            handle.write(f"{header}\n")
        handle.write(f"##sequence-region {target_seqid} 1 {target_length}\n")
        for feature in sorted(features, key=lambda item: (item.start, item.end, item.feature_type)):
            handle.write(
                "\t".join(
                    [
                        feature.seqid,
                        feature.source,
                        feature.feature_type,
                        str(feature.start),
                        str(feature.end),
                        feature.score,
                        feature.strand,
                        feature.phase,
                        feature.attributes,
                    ]
                )
                + "\n"
            )


def load_fasta_records(fasta_path: Path, limit: int = 0) -> list:
    records = list(SeqIO.parse(str(fasta_path), "fasta"))
    if limit > 0:
        return records[:limit]
    return records


def run_minimap2(minimap2_path: str, ref_path: Path, query_path: Path, paf_path: Path) -> None:
    with paf_path.open("w", encoding="utf-8") as handle:
        subprocess.run(
            [
                minimap2_path,
                "-x",
                "asm5",
                "-c",
                "--eqx",
                str(ref_path),
                str(query_path),
            ],
            stdout=handle,
            stderr=subprocess.DEVNULL,
            check=True,
        )


def main() -> None:
    args = parse_args()
    db_root = args.db_root.expanduser().resolve()
    out_root = args.out_root.expanduser().resolve()
    summary_path = args.summary.expanduser().resolve()
    minimap2_path = resolve_minimap2(args.minimap2)
    out_root.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)

    summary_rows: list[dict[str, object]] = []

    for subtype in args.subtypes:
        subtype_name = str(subtype).strip()
        ref_fasta = db_root / f"{subtype_name}.fna"
        ref_gff = db_root / f"{subtype_name}.gff3"
        db_fasta = db_root / f"{subtype_name}_db.fasta"
        subtype_out_dir = out_root / subtype_name

        if not ref_fasta.is_file() or not ref_gff.is_file() or not db_fasta.is_file():
            summary_rows.append(
                {
                    "subtype": subtype_name,
                    "sample": "-",
                    "status": "missing_input",
                    "mapped_features": 0,
                    "total_features": 0,
                    "query_coverage": 0.0,
                    "mapq": 0,
                    "output_gff": "",
                }
            )
            continue

        header_lines, ref_features = parse_gff(ref_gff)
        query_records = load_fasta_records(db_fasta, args.limit_per_subtype)
        if not query_records:
            summary_rows.append(
                {
                    "subtype": subtype_name,
                    "sample": "-",
                    "status": "empty_fasta",
                    "mapped_features": 0,
                    "total_features": len(ref_features),
                    "query_coverage": 0.0,
                    "mapq": 0,
                    "output_gff": "",
                }
            )
            continue

        with tempfile.TemporaryDirectory(prefix=f"{subtype_name}_gff_transfer_") as tmp_dir:
            query_path = Path(tmp_dir) / f"{subtype_name}_queries.fasta"
            SeqIO.write(query_records, str(query_path), "fasta")
            paf_path = Path(tmp_dir) / f"{subtype_name}.paf"
            run_minimap2(minimap2_path, ref_fasta, query_path, paf_path)

            alignments_by_query: dict[str, list[PafAlignment]] = {}
            with paf_path.open("r", encoding="utf-8", errors="ignore") as handle:
                for line in handle:
                    aln = parse_paf_line(line)
                    if aln is None:
                        continue
                    alignments_by_query.setdefault(aln.query_name, []).append(aln)

        for record in query_records:
            record_id = str(record.id).strip()
            best = choose_best_alignment(
                alignments_by_query.get(record_id, []),
                min_mapq=args.min_mapq,
                min_query_coverage=args.min_query_coverage,
            )
            if best is None:
                summary_rows.append(
                    {
                        "subtype": subtype_name,
                        "sample": record_id,
                        "status": "no_valid_alignment",
                        "mapped_features": 0,
                        "total_features": len(ref_features),
                        "query_coverage": 0.0,
                        "mapq": 0,
                        "output_gff": "",
                    }
                )
                continue

            transferred_features = []
            for feature in ref_features:
                mapped = transfer_feature(feature, best, record_id)
                if mapped is not None:
                    transferred_features.append(mapped)

            output_path = subtype_out_dir / f"{sanitize_filename(record_id)}.gff3"
            write_gff(output_path, record_id, len(record.seq), header_lines, transferred_features)
            summary_rows.append(
                {
                    "subtype": subtype_name,
                    "sample": record_id,
                    "status": "ok" if transferred_features else "no_mapped_features",
                    "mapped_features": len(transferred_features),
                    "total_features": len(ref_features),
                    "query_coverage": round((best.query_end - best.query_start) / best.query_length, 6) if best.query_length else 0.0,
                    "mapq": best.mapq,
                    "output_gff": str(output_path),
                }
            )

    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "subtype",
                "sample",
                "status",
                "mapped_features",
                "total_features",
                "query_coverage",
                "mapq",
                "output_gff",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"已完成 HPIV GFF 迁移预测，结果目录: {out_root}")
    print(f"汇总文件: {summary_path}")


if __name__ == "__main__":
    main()
