#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path


DB_NAME = "SARSCoV2_WuhanHu1"


def resolve_bcftools_path() -> str:
    candidates = [
        shutil.which("bcftools"),
        "/opt/homebrew/Caskroom/mambaforge/base/envs/ncov/bin/bcftools",
        "/opt/homebrew/Caskroom/mambaforge/base/envs/Pathogen/bin/bcftools",
        "/opt/homebrew/bin/bcftools",
        "/usr/local/bin/bcftools",
    ]
    for candidate in candidates:
        if candidate and Path(candidate).is_file():
            return str(Path(candidate).resolve())
    raise FileNotFoundError("未找到 bcftools，可执行文件不存在于 PATH 或预设的 conda 环境路径中。")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Annotate SARS-CoV-2 VCF with snpEff and export a mutation table."
    )
    parser.add_argument("vcf", type=Path, help="Input VCF path, e.g. snps.raw.vcf")
    parser.add_argument(
        "--snpeff",
        type=Path,
        default=Path("snpEff/exec/snpeff"),
        help="Path to snpEff launcher",
    )
    parser.add_argument(
        "--ncov-db-dir",
        type=Path,
        default=Path("database/virus/ncov"),
        help="Directory containing ref.fna and genome.gff",
    )
    parser.add_argument(
        "--out-prefix",
        type=Path,
        default=None,
        help="Output prefix; defaults to input basename without .vcf",
    )
    return parser.parse_args()


def parse_info_field(info: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in info.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
        else:
            result[item] = "true"
    return result


def max_homopolymer_run(sequence: str | None) -> int:
    text = str(sequence or "").strip().upper()
    if not text:
        return 0
    best = 1
    current = 1
    for index in range(1, len(text)):
        if text[index] == text[index - 1]:
            current += 1
            if current > best:
                best = current
        else:
            current = 1
    return best


def is_poly_variant(ref: str | None, alt: str | None) -> bool:
    return max(max_homopolymer_run(ref), max_homopolymer_run(alt)) >= 5


def first_float(value: str | None, default: float = 0.0) -> float:
    if not value:
        return default
    token = value.split(",")[0]
    try:
        return float(token)
    except ValueError:
        return default


def first_int(value: str | None, default: int = 0) -> int:
    if not value:
        return default
    token = value.split(",")[0]
    try:
        return int(float(token))
    except ValueError:
        return default


def choose_best_ann(ann_value: str | None) -> dict[str, str]:
    empty = {
        "allele": "",
        "annotation": "",
        "impact": "",
        "gene_name": "",
        "gene_id": "",
        "feature_type": "",
        "feature_id": "",
        "biotype": "",
        "rank": "",
        "hgvs_c": "",
        "hgvs_p": "",
        "cdna": "",
        "cds": "",
        "aa": "",
        "distance": "",
        "messages": "",
    }
    if not ann_value:
        return empty

    impact_rank = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3, "": 4}
    best_parts: list[str] | None = None
    best_key: tuple[int, int, int] | None = None

    for entry in ann_value.split(","):
        parts = entry.split("|")
        parts += [""] * (16 - len(parts))
        impact = parts[2]
        feature_type = parts[5]
        annotation = parts[1]
        key = (
            impact_rank.get(impact, 9),
            0 if feature_type == "transcript" else 1,
            0 if "intergenic" not in annotation else 1,
        )
        if best_key is None or key < best_key:
            best_key = key
            best_parts = parts[:16]

    if best_parts is None:
        return empty

    return {
        "allele": best_parts[0],
        "annotation": best_parts[1],
        "impact": best_parts[2],
        "gene_name": best_parts[3],
        "gene_id": best_parts[4],
        "feature_type": best_parts[5],
        "feature_id": best_parts[6],
        "biotype": best_parts[7],
        "rank": best_parts[8],
        "hgvs_c": best_parts[9],
        "hgvs_p": best_parts[10],
        "cdna": best_parts[11],
        "cds": best_parts[12],
        "aa": best_parts[13],
        "distance": best_parts[14],
        "messages": best_parts[15],
    }


def ensure_snpeff_db(snpeff_path: Path, ncov_db_dir: Path, project_root: Path) -> None:
    ref_fasta = ncov_db_dir / "ref.fna"
    genome_gff = ncov_db_dir / "genome.gff"
    if not genome_gff.exists():
        alt_gff = ncov_db_dir / "genomic.gff"
        if alt_gff.exists():
            genome_gff = alt_gff
        else:
            raise FileNotFoundError(f"GFF not found: {genome_gff}")
    if not ref_fasta.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {ref_fasta}")

    snpeff_root = snpeff_path.resolve().parents[1]
    config_path = snpeff_root / "snpEff.config"
    db_dir = snpeff_root / "data" / DB_NAME
    db_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(ref_fasta, db_dir / "sequences.fa")
    shutil.copy2(genome_gff, db_dir / "genes.gff")

    config_text = config_path.read_text(encoding="utf-8")
    genome_entry = f"{DB_NAME}.genome : SARS-CoV-2 Wuhan-Hu-1"
    ref_entry = f"{DB_NAME}.reference : {ref_fasta.resolve()}"
    updates: list[str] = []
    if genome_entry not in config_text:
        updates.append(genome_entry)
    if ref_entry not in config_text:
        updates.append(ref_entry)
    if updates:
        with config_path.open("a", encoding="utf-8") as handle:
            handle.write("\n" + "\n".join(updates) + "\n")

    predictor_bin = db_dir / "snpEffectPredictor.bin"
    if predictor_bin.exists():
        return

    subprocess.run(
        [
            str(snpeff_path.resolve()),
            "build",
            "-gff3",
            "-noCheckCds",
            "-noCheckProtein",
            DB_NAME,
        ],
        cwd=project_root,
        check=True,
    )


def annotate_vcf(snpeff_path: Path, input_vcf: Path, output_vcf: Path, project_root: Path) -> None:
    with output_vcf.open("w", encoding="utf-8") as handle:
        subprocess.run(
            [
                str(snpeff_path.resolve()),
                "-noStats",
                DB_NAME,
                str(input_vcf.resolve()),
            ],
            cwd=project_root,
            check=True,
            stdout=handle,
        )


def resolve_variant_maf(info_map: dict[str, str], dp: int) -> float:
    ao = first_float(info_map.get("AO"), 0.0)
    ro = first_float(info_map.get("RO"), 0.0)
    allele_depth = ao + ro
    if allele_depth > 0:
        return round((ao / allele_depth), 6)
    maf = first_float(info_map.get("AF"), -1.0)
    if maf >= 0:
        return maf
    return round((ao / dp), 6) if dp > 0 else 0.0


def classify_variant_quality(qual: float, dp: int, maf: float, ref: str | None, alt: str | None) -> str:
    maf_threshold = 0.75 if is_poly_variant(ref, alt) else 0.1
    return "高质量突变" if qual > 10 and dp > 10 and maf > maf_threshold else "低质量突变"


def build_high_quality_consensus(
    *,
    annotated_vcf: Path,
    ref_fasta: Path,
    high_quality_vcf: Path,
    high_quality_vcf_gz: Path,
    consensus_fasta: Path,
) -> None:
    bcftools = resolve_bcftools_path()
    with annotated_vcf.open("r", encoding="utf-8") as in_handle, high_quality_vcf.open("w", encoding="utf-8") as out_handle:
        for line in in_handle:
            if not line:
                continue
            if line.startswith("#"):
                out_handle.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            ref = fields[3]
            alt = fields[4]
            qual = first_float(fields[5], 0.0)
            info_map = parse_info_field(fields[7])
            dp = first_int(info_map.get("DP"), 0)
            maf = resolve_variant_maf(info_map, dp)
            if classify_variant_quality(qual, dp, maf, ref, alt) == "高质量突变":
                out_handle.write(line)
    subprocess.run(
        [
            bcftools,
            "view",
            "-Oz",
            "-o",
            str(high_quality_vcf_gz),
            str(high_quality_vcf),
        ],
        check=True,
    )
    subprocess.run(
        [
            bcftools,
            "index",
            "-f",
            str(high_quality_vcf_gz),
        ],
        check=True,
    )
    with consensus_fasta.open("w", encoding="utf-8") as handle:
        subprocess.run(
            [
                bcftools,
                "consensus",
                "-f",
                str(ref_fasta),
                str(high_quality_vcf_gz),
            ],
            check=True,
            stdout=handle,
        )


def export_tables(annotated_vcf: Path, out_tsv: Path, out_json: Path) -> dict[str, int]:
    rows = []
    total = 0
    high_quality = 0

    with annotated_vcf.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            total += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            chrom, pos, record_id, ref, alt, qual_raw, filt, info = fields[:8]
            info_map = parse_info_field(info)
            ann = choose_best_ann(info_map.get("ANN"))

            qual = first_float(qual_raw, 0.0)
            dp = first_int(info_map.get("DP"), 0)
            maf = resolve_variant_maf(info_map, dp)
            quality_label = classify_variant_quality(qual, dp, maf, ref, alt)
            if quality_label == "高质量突变":
                high_quality += 1

            variant_type = info_map.get("TYPE", "")
            row = {
                "染色体": chrom,
                "位置": int(pos),
                "ID": record_id,
                "参考碱基": ref,
                "突变碱基": alt,
                "核苷酸突变": f"{ref}{pos}{alt}",
                "变异类型": variant_type,
                "质量值QUAL": qual,
                "测序深度DP": dp,
                "突变频率MAF": round(maf, 6),
                "质量分层": quality_label,
                "注释效应": ann["annotation"],
                "影响等级": ann["impact"],
                "基因": ann["gene_name"],
                "基因ID": ann["gene_id"],
                "转录本ID": ann["feature_id"],
                "HGVS.c": ann["hgvs_c"],
                "HGVS.p": ann["hgvs_p"],
                "氨基酸位点": ann["aa"],
                "警告信息": ann["messages"],
                "FILTER": filt,
            }
            rows.append(row)

    fieldnames = [
        "染色体",
        "位置",
        "ID",
        "参考碱基",
        "突变碱基",
        "核苷酸突变",
        "变异类型",
        "质量值QUAL",
        "测序深度DP",
        "突变频率MAF",
        "质量分层",
        "注释效应",
        "影响等级",
        "基因",
        "基因ID",
        "转录本ID",
        "HGVS.c",
        "HGVS.p",
        "氨基酸位点",
        "警告信息",
        "FILTER",
    ]

    with out_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    with out_json.open("w", encoding="utf-8") as handle:
        json.dump(
            {
                "status": "ready",
                "source_vcf": str(annotated_vcf),
                "total_variants": total,
                "high_quality_variants": high_quality,
                "low_quality_variants": total - high_quality,
                "rows": rows,
            },
            handle,
            ensure_ascii=False,
            indent=2,
        )

    return {
        "total": total,
        "high_quality": high_quality,
        "low_quality": total - high_quality,
    }


def main() -> int:
    args = parse_args()
    project_root = Path(__file__).resolve().parents[1]
    input_vcf = args.vcf.resolve()
    snpeff_path = (project_root / args.snpeff).resolve() if not args.snpeff.is_absolute() else args.snpeff.resolve()
    ncov_db_dir = (project_root / args.ncov_db_dir).resolve() if not args.ncov_db_dir.is_absolute() else args.ncov_db_dir.resolve()

    if args.out_prefix:
        out_prefix = args.out_prefix.resolve() if args.out_prefix.is_absolute() else (project_root / args.out_prefix).resolve()
    else:
        suffix = input_vcf.name[:-4] if input_vcf.name.endswith(".vcf") else input_vcf.name
        out_prefix = input_vcf.with_name(suffix)

    annotated_vcf = out_prefix.with_name(out_prefix.name + ".ann.vcf")
    out_tsv = out_prefix.with_name(out_prefix.name + ".mutation_table.tsv")
    out_json = out_prefix.with_name(out_prefix.name + ".mutation_table.json")
    high_quality_vcf = out_prefix.with_name(out_prefix.name + ".high_quality.vcf")
    high_quality_vcf_gz = out_prefix.with_name(out_prefix.name + ".high_quality.vcf.gz")
    consensus_fasta = out_prefix.with_name(out_prefix.name + ".high_quality.consensus.fasta")
    ref_fasta = ncov_db_dir / "ref.fna"

    ensure_snpeff_db(snpeff_path, ncov_db_dir, project_root)
    annotate_vcf(snpeff_path, input_vcf, annotated_vcf, project_root)
    stats = export_tables(annotated_vcf, out_tsv, out_json)
    build_high_quality_consensus(
        annotated_vcf=annotated_vcf,
        ref_fasta=ref_fasta,
        high_quality_vcf=high_quality_vcf,
        high_quality_vcf_gz=high_quality_vcf_gz,
        consensus_fasta=consensus_fasta,
    )

    print(
        json.dumps(
            {
                "status": "ready",
                "annotated_vcf": str(annotated_vcf),
                "high_quality_vcf": str(high_quality_vcf),
                "high_quality_vcf_gz": str(high_quality_vcf_gz),
                "consensus_fasta": str(consensus_fasta),
                "mutation_table_tsv": str(out_tsv),
                "mutation_table_json": str(out_json),
                **stats,
            },
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
