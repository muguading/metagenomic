#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


AA3_TO_1 = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Gln": "Q",
    "Glu": "E",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
    "Ter": "*",
    "Stop": "*",
}

GRADE_TO_CALL = {
    "1) Assoc w R": "耐药相关",
    "2) Assoc w R - Interim": "耐药相关（临时）",
    "3) Uncertain significance": "意义未明",
    "4) Not assoc w R - Interim": "未证实耐药（临时）",
    "5) Not assoc w R": "未证实耐药",
}
GRADE_PRIORITY = {
    "1) Assoc w R": 0,
    "2) Assoc w R - Interim": 1,
    "3) Uncertain significance": 2,
    "4) Not assoc w R - Interim": 3,
    "5) Not assoc w R": 4,
}
PREDICTION_PRIORITY = {
    "耐药相关": 0,
    "耐药相关（临时）": 1,
    "意义未明": 2,
    "未证实耐药（临时）": 3,
    "未证实耐药": 4,
}
TB_DRUG_LINE_MAP = {
    "Isoniazid": "一线药物",
    "Rifampicin": "一线药物",
    "Ethambutol": "一线药物",
    "Pyrazinamide": "一线药物",
    "Levofloxacin": "二线药物",
    "Moxifloxacin": "二线药物",
    "Amikacin": "二线药物",
    "Kanamycin": "二线药物",
    "Capreomycin": "二线药物",
    "Linezolid": "二线药物",
    "Bedaquiline": "二线药物",
    "Clofazimine": "二线药物",
    "Delamanid": "二线药物",
    "Ethionamide": "二线药物",
    "Streptomycin": "二线药物",
}
TB_NAME_TOKENS = (
    "mycobacterium tuberculosis",
    "mycobacterium_tuberculosis",
    "mycobacterium tuberculosis complex",
    "m. tuberculosis",
    "mtb",
    "结核分枝杆菌",
    "结核杆菌",
)


def compact_text(value: object) -> str:
    return str(value or "").strip()


def classify_tb_drug_line(drug_name: object) -> str:
    return TB_DRUG_LINE_MAP.get(compact_text(drug_name), "未分层")


def _parse_gff_attributes(field: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for chunk in str(field or "").strip().split(";"):
        item = chunk.strip()
        if not item or "=" not in item:
            continue
        key, value = item.split("=", 1)
        key = key.strip()
        value = value.strip()
        if key:
            attributes[key] = value
    return attributes


def prepare_snpeff_reference_gff(source_gff: Path, output_gff: Path) -> bool:
    if not source_gff.is_file() or source_gff.stat().st_size == 0:
        return False
    cds_rows: list[tuple[str, str, int, int, str, str, str, dict[str, str]]] = []
    gene_rows: list[tuple[str, str, int, int, str, str, dict[str, str]]] = []
    sequence_regions: list[str] = []
    seen_regions: set[str] = set()
    with source_gff.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith("##sequence-region"):
                if line not in seen_regions:
                    sequence_regions.append(line)
                    seen_regions.add(line)
                continue
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 9:
                continue
            seqid, source, feature_type, start, end, score, strand, phase, attrs = parts
            if feature_type not in {"CDS", "gene"}:
                continue
            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                continue
            attributes = _parse_gff_attributes(attrs)
            normalized_source = source or "RefSeq"
            normalized_score = score or "."
            normalized_strand = strand or "+"
            if feature_type == "CDS":
                cds_rows.append((seqid, normalized_source, start_i, end_i, normalized_score, normalized_strand, phase or "0", attributes))
            else:
                gene_rows.append((seqid, normalized_source, start_i, end_i, normalized_score, normalized_strand, attributes))
    if not cds_rows and not gene_rows:
        shutil.copy2(source_gff, output_gff)
        return True
    output_gff.parent.mkdir(parents=True, exist_ok=True)
    with output_gff.open("w", encoding="utf-8") as handle:
        handle.write("##gff-version 3\n")
        for line in sequence_regions:
            handle.write(f"{line}\n")
        for seqid, source, start_i, end_i, score, strand, attrs in gene_rows:
            gene_id = attrs.get("ID") or attrs.get("locus_tag") or attrs.get("gene") or f"gene_{seqid}_{start_i}_{end_i}"
            gene_name = attrs.get("gene") or attrs.get("Name") or gene_id
            gene_attrs = f"ID={gene_id};Name={gene_name};gene={gene_name}"
            handle.write(f"{seqid}\t{source}\tgene\t{start_i}\t{end_i}\t{score}\t{strand}\t.\t{gene_attrs}\n")
        for seqid, source, start_i, end_i, score, strand, phase, attrs in cds_rows:
            parent_id = attrs.get("Parent") or attrs.get("gene") or attrs.get("locus_tag") or attrs.get("ID") or f"gene_{seqid}_{start_i}_{end_i}"
            cds_id = attrs.get("ID") or f"{parent_id}.cds"
            gene_name = attrs.get("gene") or attrs.get("Name") or parent_id
            product = attrs.get("product") or gene_name
            cds_attrs = f"ID={cds_id};Parent={parent_id};Name={gene_name};gene={gene_name};product={product}"
            handle.write(f"{seqid}\t{source}\tCDS\t{start_i}\t{end_i}\t{score}\t{strand}\t{phase}\t{cds_attrs}\n")
    return True


def normalize_tb_text(value: object) -> str:
    return re.sub(r"\s+", " ", compact_text(value).lower())


def looks_like_tb_species(value: object) -> bool:
    text = normalize_tb_text(value)
    return any(token in text for token in TB_NAME_TOKENS)


def run_shell(command: str, *, cwd: Path, log_path: Path) -> None:
    with log_path.open("a", encoding="utf-8") as log_handle:
        log_handle.write(f"$ {command}\n")
        log_handle.flush()
        result = subprocess.run(
            command,
            shell=True,
            cwd=str(cwd),
            stdout=log_handle,
            stderr=log_handle,
            text=True,
            check=False,
        )
    if result.returncode != 0:
        raise RuntimeError(f"命令执行失败: {command}")


def run_command(command: list[str], *, cwd: Path, log_path: Path, env: dict[str, str] | None = None) -> None:
    with log_path.open("a", encoding="utf-8") as log_handle:
        log_handle.write(f"$ {' '.join(shlex.quote(part) for part in command)}\n")
        log_handle.flush()
        result = subprocess.run(
            command,
            cwd=str(cwd),
            stdout=log_handle,
            stderr=log_handle,
            text=True,
            env=env,
            check=False,
        )
    if result.returncode != 0:
        raise RuntimeError(f"命令执行失败: {' '.join(command)}")


def read_first_tsv_row(path: Path) -> dict[str, str]:
    if not path.is_file():
        return {}
    with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        row = next(reader, None)
    return {str(key): compact_text(value) for key, value in (row or {}).items()}


def detect_inputs(report_dir: Path, sample_name: str) -> dict[str, str]:
    pair_r1 = report_dir / f"{sample_name}.R1.fastq.gz"
    pair_r2 = report_dir / f"{sample_name}.R2.fastq.gz"
    single_r1 = report_dir / f"{sample_name}.R1.fastq.gz"
    long_fastq = report_dir / f"{sample_name}.final.fastq"
    final_fasta = report_dir / f"{sample_name}.final.fasta"
    if pair_r1.is_file() and pair_r2.is_file():
        return {"mode": "paired_short", "read1": str(pair_r1), "read2": str(pair_r2)}
    if single_r1.is_file():
        return {"mode": "single_short", "read1": str(single_r1)}
    if long_fastq.is_file():
        return {"mode": "long", "read1": str(long_fastq)}
    if final_fasta.is_file():
        return {"mode": "fasta", "fasta": str(final_fasta)}
    raise FileNotFoundError(f"{sample_name} 未找到可用于结核后处理的 reads 或 fasta")


def build_snpeff_db(work_dir: Path, ref_fasta: Path, ref_gff: Path, log_path: Path, project_root: Path) -> None:
    ref_dir = work_dir / "ref"
    genomes_dir = work_dir / "genomes"
    ref_dir.mkdir(parents=True, exist_ok=True)
    genomes_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(ref_fasta, genomes_dir / "ref.fa")
    prepared = prepare_snpeff_reference_gff(ref_gff, ref_dir / "genes.gff")
    if not prepared:
        shutil.copy2(ref_gff, ref_dir / "genes.gff")
    (work_dir / "snpEff.config").write_text("ref.genome : ref\n", encoding="utf-8")
    snpeff_jar = (project_root / "snpEff" / "snpEff.jar").resolve()
    if not snpeff_jar.is_file():
        raise FileNotFoundError(f"未找到 snpEff.jar: {snpeff_jar}")
    run_command(
        [
            "java",
            "-jar",
            str(snpeff_jar),
            "build",
            "-noCheckCds",
            "-noCheckProtein",
            "-gff3",
            "ref",
            "-c",
            "snpEff.config",
            "-dataDir",
            ".",
        ],
        cwd=work_dir,
        log_path=log_path,
    )


def call_variants(
    *,
    work_dir: Path,
    sample_name: str,
    inputs: dict[str, str],
    threads: int,
    platform: str,
    log_path: Path,
    project_root: Path,
) -> tuple[Path, Path, Path, Path]:
    ref_fasta = (project_root / "database" / "bacteria" / "tb" / "GCF_000195955.2_ASM19595v2_genomic.fna").resolve()
    ref_gff = (project_root / "database" / "bacteria" / "tb" / "genomic.gff").resolve()
    if not ref_fasta.is_file() or not ref_gff.is_file():
        raise FileNotFoundError("未找到 H37Rv 参考基因组或注释文件")
    build_snpeff_db(work_dir, ref_fasta, ref_gff, log_path, project_root)
    bam_path = work_dir / "ref.mapping.bam"
    genomes_ref = work_dir / "genomes" / "ref.fa"
    mode = inputs.get("mode", "")
    if mode == "paired_short":
        run_shell(
            f"bwa index {shlex.quote(str(genomes_ref))}",
            cwd=work_dir,
            log_path=log_path,
        )
        run_shell(
            "bwa mem -t {threads} {ref} {r1} {r2} | samtools sort -o {bam}".format(
                threads=threads,
                ref=shlex.quote(str(genomes_ref)),
                r1=shlex.quote(inputs["read1"]),
                r2=shlex.quote(inputs["read2"]),
                bam=shlex.quote(str(bam_path)),
            ),
            cwd=work_dir,
            log_path=log_path,
        )
    elif mode == "single_short":
        run_shell(
            f"bwa index {shlex.quote(str(genomes_ref))}",
            cwd=work_dir,
            log_path=log_path,
        )
        run_shell(
            "bwa mem -t {threads} {ref} {r1} | samtools sort -o {bam}".format(
                threads=threads,
                ref=shlex.quote(str(genomes_ref)),
                r1=shlex.quote(inputs["read1"]),
                bam=shlex.quote(str(bam_path)),
            ),
            cwd=work_dir,
            log_path=log_path,
        )
    elif mode == "long":
        preset = "map-ont" if platform == "nanopore" else "map-pb"
        run_shell(
            "minimap2 -ax {preset} {ref} {reads} -t {threads} | samtools sort -o {bam}".format(
                preset=preset,
                ref=shlex.quote(str(genomes_ref)),
                reads=shlex.quote(inputs["read1"]),
                threads=threads,
                bam=shlex.quote(str(bam_path)),
            ),
            cwd=work_dir,
            log_path=log_path,
        )
    else:
        raise RuntimeError("缺少可用于 H37Rv 有参 SNP 识别的 reads")
    run_command(["samtools", "index", str(bam_path)], cwd=work_dir, log_path=log_path)
    run_command(["samtools", "faidx", str(genomes_ref)], cwd=work_dir, log_path=log_path)
    run_shell(
        f"fasta_generate_regions.py {shlex.quote(str(genomes_ref))}.fai 200000 > ref.txt",
        cwd=work_dir,
        log_path=log_path,
    )
    run_shell(
        (
            "freebayes-parallel ref.txt {threads} -p 2 -P 0 -C 2 -F 0.05 "
            "--min-coverage 10 --min-repeat-entropy 1.0 -q 30 -m 30 --strict-vcf "
            "-f {ref} {bam} > snps.raw.vcf"
        ).format(
            threads=threads,
            ref=shlex.quote(str(genomes_ref)),
            bam=shlex.quote(str(bam_path)),
        ),
        cwd=work_dir,
        log_path=log_path,
    )
    run_shell(
        (
            "bcftools view --include "
            "\"QUAL>=20 && FMT/DP>=10 && (FMT/AO+FMT/RO)>0 && (FMT/AO)/(FMT/AO+FMT/RO)>=0.9\" "
            "snps.raw.vcf | "
            "bcftools annotate --remove "
            "'^INFO/TYPE,^INFO/DP,^INFO/RO,^INFO/AO,^INFO/AB,^FORMAT/GT,^FORMAT/DP,^FORMAT/RO,^FORMAT/AO,^FORMAT/QR,^FORMAT/QA,^FORMAT/GL' "
            "> snps.filt1.vcf"
        ),
        cwd=work_dir,
        log_path=log_path,
    )
    filtered_vcf = work_dir / "snps.filt1.vcf"
    filtered_vcfgz = work_dir / "snps.filt1.vcf.gz"
    run_command(
        ["bcftools", "view", "-Oz", "-o", str(filtered_vcfgz), str(filtered_vcf)],
        cwd=work_dir,
        log_path=log_path,
    )
    run_command(
        ["bcftools", "index", "-f", str(filtered_vcfgz)],
        cwd=work_dir,
        log_path=log_path,
    )
    snpeff_jar = (project_root / "snpEff" / "snpEff.jar").resolve()
    annotated_vcf = work_dir / "snps.anno.vcf"
    with annotated_vcf.open("w", encoding="utf-8") as handle:
        result = subprocess.run(
            [
                "java",
                "-jar",
                str(snpeff_jar),
                "ann",
                "-c",
                "snpEff.config",
                "-dataDir",
                ".",
                "ref",
                str(filtered_vcfgz),
            ],
            cwd=str(work_dir),
            stdout=handle,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
    if result.returncode != 0:
        raise RuntimeError(f"snpEff 注释失败: {result.stderr}")
    anno_tsv = work_dir / f"{sample_name}.anno.tsv"
    export_annotated_vcf_table(annotated_vcf, anno_tsv)
    consensus_fasta = work_dir / f"{sample_name}.tbprofiler.consensus.fa"
    run_command(
        [
            "bcftools",
            "consensus",
            "-f",
            str(genomes_ref),
            "-o",
            str(consensus_fasta),
            str(filtered_vcfgz),
        ],
        cwd=work_dir,
        log_path=log_path,
    )
    return bam_path, filtered_vcfgz, anno_tsv, consensus_fasta


def parse_info_map(info_field: str) -> dict[str, str]:
    info_map: dict[str, str] = {}
    for chunk in compact_text(info_field).split(";"):
        if not chunk:
            continue
        if "=" in chunk:
            key, value = chunk.split("=", 1)
            info_map[key] = value
        else:
            info_map[chunk] = ""
    return info_map


def pick_primary_ann(info_map: dict[str, str]) -> list[str]:
    ann_value = compact_text(info_map.get("ANN"))
    if not ann_value:
        return []
    first = ann_value.split(",")[0]
    return first.split("|")


def export_annotated_vcf_table(annotated_vcf: Path, output_tsv: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with annotated_vcf.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, _vid, ref, alt, _qual, _flt, info = parts[:8]
            info_map = parse_info_map(info)
            ann_fields = pick_primary_ann(info_map)
            effect = ann_fields[1] if len(ann_fields) > 1 else ""
            impact = ann_fields[2] if len(ann_fields) > 2 else ""
            gene = ann_fields[3] if len(ann_fields) > 3 else ""
            hgvs_c = ann_fields[9] if len(ann_fields) > 9 else ""
            hgvs_p = ann_fields[10] if len(ann_fields) > 10 else ""
            row = {
                "片段名称": compact_text(chrom),
                "突变位置": compact_text(pos),
                "参考碱基": compact_text(ref),
                "突变碱基": compact_text(alt),
                "影响基因": compact_text(gene),
                "突变类型": compact_text(effect),
                "突变影响": compact_text(impact),
                "碱基变化": compact_text(hgvs_c),
                "氨基酸变化": compact_text(hgvs_p),
            }
            rows.append(row)
    with output_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["片段名称", "突变位置", "参考碱基", "突变碱基", "影响基因", "突变类型", "突变影响", "碱基变化", "氨基酸变化"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)
    return rows


def normalize_gene_name(value: object) -> str:
    return compact_text(value).replace(" ", "").lower()


def _aa_to_one_letter(token: str) -> str:
    token = compact_text(token)
    if not token:
        return token
    if len(token) == 1:
        return token
    return AA3_TO_1.get(token[:1].upper() + token[1:].lower(), token)


def normalize_protein_change(value: object) -> str:
    text = compact_text(value)
    if not text:
        return ""
    text = re.sub(r"^p\.", "", text)
    text = text.replace("(", "").replace(")", "")
    text = text.replace("Ter", "*")
    match = re.match(r"([A-Za-z\*]+)(\d+)([A-Za-z\*=]+)", text)
    if match:
        ref_aa = _aa_to_one_letter(match.group(1))
        pos = match.group(2)
        alt_aa = _aa_to_one_letter(match.group(3))
        return f"p.{ref_aa}{pos}{alt_aa}"
    return f"p.{text}" if not text.startswith("p.") else text


def normalize_catalogue_mutation(value: object) -> str:
    text = compact_text(value)
    if not text:
        return ""
    if text.startswith("p.") or re.match(r"^[A-Za-z\*]{1,3}\d+[A-Za-z\*=]{1,3}$", text):
        return normalize_protein_change(text)
    return text.replace(" ", "")


def load_catalogue_index(catalogue_path: Path) -> dict[tuple[str, str], list[dict[str, str]]]:
    index: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    with catalogue_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            gene = normalize_gene_name(row.get("gene"))
            mutation = normalize_catalogue_mutation(row.get("mutation"))
            if not gene or not mutation:
                continue
            payload = {str(key): compact_text(value) for key, value in row.items()}
            index[(gene, mutation)].append(payload)
    return index


def read_annotation_rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [{str(key): compact_text(value) for key, value in row.items()} for row in reader]


def match_catalogue(
    annotation_rows: list[dict[str, str]],
    catalogue_index: dict[tuple[str, str], list[dict[str, str]]],
) -> list[dict[str, str]]:
    matched_rows: list[dict[str, str]] = []
    seen: set[tuple[str, str, str, str]] = set()
    for row in annotation_rows:
        gene = normalize_gene_name(row.get("影响基因"))
        if not gene:
            continue
        sample_changes = [
            ("核酸变化", compact_text(row.get("碱基变化")).replace(" ", "")),
            ("蛋白变化", normalize_protein_change(row.get("氨基酸变化"))),
        ]
        for match_level, mutation in sample_changes:
            if not mutation:
                continue
            for catalogue_row in catalogue_index.get((gene, mutation), []):
                dedupe_key = (
                    compact_text(catalogue_row.get("drug")),
                    compact_text(catalogue_row.get("gene")),
                    compact_text(catalogue_row.get("mutation")),
                    compact_text(row.get("突变位置")),
                )
                if dedupe_key in seen:
                    continue
                seen.add(dedupe_key)
                final_grade = compact_text(catalogue_row.get("FINAL CONFIDENCE GRADING"))
                matched_rows.append(
                    {
                        "药物": compact_text(catalogue_row.get("drug")),
                        "药物分层": classify_tb_drug_line(catalogue_row.get("drug")),
                        "基因": compact_text(catalogue_row.get("gene")),
                        "catalogue突变": compact_text(catalogue_row.get("mutation")),
                        "样本匹配层级": match_level,
                        "样本碱基变化": compact_text(row.get("碱基变化")),
                        "样本氨基酸变化": compact_text(row.get("氨基酸变化")),
                        "突变位置": compact_text(row.get("突变位置")),
                        "参考碱基": compact_text(row.get("参考碱基")),
                        "突变碱基": compact_text(row.get("突变碱基")),
                        "突变类型": compact_text(row.get("突变类型")),
                        "突变影响": compact_text(row.get("突变影响")),
                        "最终分级": final_grade,
                        "判读结论": GRADE_TO_CALL.get(final_grade, final_grade or "-"),
                        "注释": compact_text(catalogue_row.get("Comment")),
                        "variant": compact_text(catalogue_row.get("variant")),
                        "effect": compact_text(catalogue_row.get("effect")),
                    }
                )
    matched_rows.sort(
        key=lambda item: (
            PREDICTION_PRIORITY.get(compact_text(item.get("判读结论")), 99),
            compact_text(item.get("药物")).lower(),
            compact_text(item.get("基因")).lower(),
            compact_text(item.get("catalogue突变")),
        )
    )
    return matched_rows


def build_resistance_summary(matched_rows: list[dict[str, str]]) -> dict[str, object]:
    by_drug: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in matched_rows:
        by_drug[compact_text(row.get("药物"))].append(row)
    drug_calls: list[dict[str, object]] = []
    focus_drugs = 0
    other_catalogue_drugs = 0
    for drug, rows in sorted(by_drug.items(), key=lambda item: item[0].lower()):
        best = sorted(
            rows,
            key=lambda item: (
                GRADE_PRIORITY.get(compact_text(item.get("最终分级")), 99),
                compact_text(item.get("catalogue突变")),
            ),
        )[0]
        verdict = compact_text(best.get("判读结论")) or "-"
        final_grade = compact_text(best.get("最终分级")) or "-"
        is_focus = final_grade in {"1) Assoc w R", "2) Assoc w R - Interim"}
        if is_focus:
            focus_drugs += 1
        else:
            other_catalogue_drugs += 1
        drug_calls.append(
            {
                "drug": drug,
                "drug_line": classify_tb_drug_line(drug),
                "verdict": verdict,
                "grade": final_grade,
                "is_focus": is_focus,
                "mutations": [
                    compact_text(item.get("catalogue突变"))
                    for item in rows
                    if compact_text(item.get("catalogue突变"))
                ],
            }
        )
    focus_calls = [item for item in drug_calls if item.get("is_focus")]
    other_calls = [item for item in drug_calls if not item.get("is_focus")]
    if drug_calls:
        if focus_calls:
            headline = (
                f"基于 H37Rv 有参 SNP 与 WHO mutation catalogue，"
                f"共识别到 {focus_drugs} 个药物存在 WHO 1/2 级重点耐药证据。"
            )
        else:
            headline = (
                f"基于 H37Rv 有参 SNP 与 WHO mutation catalogue，"
                f"当前未识别到 WHO 1/2 级重点耐药证据。"
            )
    else:
        headline = "基于 H37Rv 有参 SNP 与 WHO mutation catalogue，当前未匹配到已收录的耐药判读条目。"
    interpretation_items: list[str] = []
    if focus_calls:
        interpretation_items.append(
            "重点关注的耐药药物（WHO 1/2 级）："
            + ", ".join([str(item.get("drug") or "-") for item in focus_calls])
            + "。"
        )
        first_line_focus = [str(item.get("drug") or "-") for item in focus_calls if str(item.get("drug_line") or "") == "一线药物"]
        second_line_focus = [str(item.get("drug") or "-") for item in focus_calls if str(item.get("drug_line") or "") == "二线药物"]
        if first_line_focus:
            interpretation_items.append("其中一线药物包括：" + ", ".join(first_line_focus) + "。")
        if second_line_focus:
            interpretation_items.append("其中二线药物包括：" + ", ".join(second_line_focus) + "。")
    highlights = [
        f"{item['drug']}: {item['verdict']}（{', '.join(item['mutations'][:3])}）"
        for item in focus_calls[:8]
    ]
    return {
        "headline": headline,
        "drug_calls": drug_calls,
        "focus_drug_calls": focus_calls,
        "other_drug_calls": other_calls,
        "total_drug_count": len(drug_calls),
        "focus_drug_count": focus_drugs,
        "other_catalogue_drug_count": other_catalogue_drugs,
        "matched_variant_count": len(matched_rows),
        "interpretation_items": interpretation_items,
        "highlights": highlights,
    }


def load_tbprofiler_json(tbprofiler_dir: Path, sample_name: str) -> dict:
    candidates = [
        tbprofiler_dir / f"{sample_name}.results.json",
        tbprofiler_dir / "results" / f"{sample_name}.results.json",
    ]
    for path in candidates:
        if path.is_file():
            try:
                return json.loads(path.read_text(encoding="utf-8"))
            except json.JSONDecodeError:
                return {}
    return {}


def first_present(payload: dict, keys: list[str]) -> object:
    for key in keys:
        value = payload.get(key)
        if value not in (None, "", [], {}):
            return value
    return ""


def parse_tbprofiler_summary(payload: dict) -> dict[str, object]:
    lineage = ""
    family = ""
    lineage_items = payload.get("lineage") if isinstance(payload.get("lineage"), list) else []
    if lineage_items and all(isinstance(item, dict) for item in lineage_items):
        lineage_labels = [compact_text(item.get("lineage")) for item in lineage_items if compact_text(item.get("lineage"))]
        family_labels = [compact_text(item.get("family")) for item in lineage_items if compact_text(item.get("family"))]
        lineage = lineage_labels[-1] if lineage_labels else ""
        family = family_labels[-1] if family_labels else ""
    else:
        lineage_value = first_present(payload, ["sublin", "sublineage", "lineage", "main_lin", "main_lineage"])
        if isinstance(lineage_value, list):
            lineage = "; ".join([compact_text(item) for item in lineage_value if compact_text(item)])
        elif isinstance(lineage_value, dict):
            lineage = "; ".join([compact_text(value) for value in lineage_value.values() if compact_text(value)])
        else:
            lineage = compact_text(lineage_value)
    if not lineage:
        lineage = compact_text(first_present(payload, ["main_lineage", "main_lin"]))
    if not family:
        family = compact_text(first_present(payload, ["spoligotype", "family", "lineage_family"]))
    spoligotype = compact_text(first_present(payload, ["spoligotype", "spoligo"]))
    drug_type = compact_text(first_present(payload, ["drtype", "drug_resistance", "resistance_profile"]))
    summary_rows = []
    if lineage:
        summary_rows.append(["Lineage", lineage])
    if family:
        summary_rows.append(["Family/Spoligotype", family])
    if spoligotype and spoligotype != family:
        summary_rows.append(["Spoligotype", spoligotype])
    if drug_type:
        summary_rows.append(["tb-profiler DR type", drug_type])
    return {
        "predicted_lineage": lineage or "-",
        "family": family or "-",
        "spoligotype": spoligotype or "-",
        "drug_type": drug_type or "-",
        "lineage_rows": summary_rows,
        "raw_lineage_items": lineage_items,
    }


def run_tbprofiler(
    *,
    bam_path: Path,
    consensus_fasta: Path,
    tbprofiler_dir: Path,
    sample_name: str,
    threads: int,
    platform: str,
    env: dict[str, str],
    log_path: Path,
) -> dict[str, object]:
    tbprofiler_dir.mkdir(parents=True, exist_ok=True)
    configured_bin = compact_text(env.get("TB_PROFILER_BIN"))
    candidates = []
    if configured_bin:
        candidates.append(shlex.split(configured_bin))
    for fallback in (
        ["conda", "run", "-n", compact_text(env.get("TB_PROFILER_ENV")) or "TB", "--no-capture-output", "tb-profiler"],
        ["conda", "run", "-n", "TB_ONT", "--no-capture-output", "tb-profiler"],
        ["tb-profiler"],
    ):
        if fallback not in candidates:
            candidates.append(fallback)
    last_error = ""
    for base_cmd in candidates:
        cmd = base_cmd + [
            "profile",
            "--fasta",
            str(consensus_fasta),
            "--platform",
            platform,
            "--prefix",
            sample_name,
            "--dir",
            str(tbprofiler_dir),
            "--txt",
            "--csv",
            "--threads",
            str(max(1, threads)),
        ]
        try:
            run_command(cmd, cwd=tbprofiler_dir, log_path=log_path, env=env)
            return parse_tbprofiler_summary(load_tbprofiler_json(tbprofiler_dir, sample_name))
        except Exception as exc:  # noqa: BLE001
            last_error = str(exc)
            continue
    return {
        "predicted_lineage": "-",
        "family": "-",
        "spoligotype": "-",
        "drug_type": "-",
        "lineage_rows": [["tb-profiler", f"运行失败: {last_error or '未找到可用命令'}"]],
        "raw_lineage_items": [],
    }


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = [
        "药物",
        "药物分层",
        "基因",
        "catalogue突变",
        "样本匹配层级",
        "样本碱基变化",
        "样本氨基酸变化",
        "突变位置",
        "参考碱基",
        "突变碱基",
        "突变类型",
        "突变影响",
        "最终分级",
        "判读结论",
        "注释",
        "variant",
        "effect",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description="结核分枝杆菌后处理：H37Rv 有参 SNP + WHO catalogue 判读 + tb-profiler 家系分析")
    parser.add_argument("--report-dir", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--platform", choices=["illumina", "nanopore", "pacbio"], default="illumina")
    args = parser.parse_args()

    report_dir = args.report_dir.expanduser().resolve()
    project_root = args.project_root.expanduser().resolve()
    sample_name = compact_text(args.sample)
    tb_root = report_dir / "tb_analysis"
    ref_call_dir = tb_root / "reference_call"
    tbprofiler_dir = tb_root / "tbprofiler"
    tb_root.mkdir(parents=True, exist_ok=True)
    ref_call_dir.mkdir(parents=True, exist_ok=True)
    log_path = tb_root / "tb_postprocess.log"

    checkm_row = read_first_tsv_row(report_dir / f"{sample_name}.checkm.tsv")
    if not any(looks_like_tb_species(checkm_row.get(key)) for key in ("物种名称", "species_name", "mlst 物种名称", "mlst_species_name")):
        raise RuntimeError(f"{sample_name} 当前物种结果不是结核分枝杆菌，跳过 TB 专用后处理")

    env = os.environ.copy()
    inputs = detect_inputs(report_dir, sample_name)
    bam_path, filtered_vcf, anno_tsv, consensus_fasta = call_variants(
        work_dir=ref_call_dir,
        sample_name=sample_name,
        inputs=inputs,
        threads=max(1, args.threads),
        platform=args.platform,
        log_path=log_path,
        project_root=project_root,
    )
    tbprofiler_summary = run_tbprofiler(
        bam_path=bam_path,
        consensus_fasta=consensus_fasta,
        tbprofiler_dir=tbprofiler_dir,
        sample_name=sample_name,
        threads=max(1, args.threads),
        platform=args.platform,
        env=env,
        log_path=log_path,
    )
    catalogue_path = (
        project_root
        / "database"
        / "bacteria"
        / "tb"
        / "mutation-catalogue-2023-main"
        / "Final Result Files"
        / "WHO-UCN-TB-2023.6-eng_catalogue_master_file.txt"
    ).resolve()
    if not catalogue_path.is_file():
        raise FileNotFoundError(f"未找到 WHO catalogue 文件: {catalogue_path}")
    annotation_rows = read_annotation_rows(anno_tsv)
    catalogue_index = load_catalogue_index(catalogue_path)
    matched_rows = match_catalogue(annotation_rows, catalogue_index)
    resistance_summary = build_resistance_summary(matched_rows)
    catalogue_match_tsv = tb_root / "tb_catalogue_matches.tsv"
    write_tsv(catalogue_match_tsv, matched_rows)
    payload = {
        "status": "ready",
        "sample_name": sample_name,
        "reference": {
            "fasta": str((project_root / "database" / "bacteria" / "tb" / "GCF_000195955.2_ASM19595v2_genomic.fna").resolve()),
            "gff": str((project_root / "database" / "bacteria" / "tb" / "genomic.gff").resolve()),
            "filtered_vcf": str(filtered_vcf),
            "annotation_tsv": str(anno_tsv),
            "bam": str(bam_path),
            "consensus_fasta": str(consensus_fasta),
            "input_mode": inputs.get("mode", ""),
        },
        "tbprofiler": tbprofiler_summary,
        "catalogue": {
            "source": str(catalogue_path),
            "matched_rows": matched_rows,
            "summary": resistance_summary,
            "match_table": str(catalogue_match_tsv),
        },
    }
    (tb_root / "tb_summary.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"status": "ok", "tb_root": str(tb_root), "matched_rows": len(matched_rows)}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
