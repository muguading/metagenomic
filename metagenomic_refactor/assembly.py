from __future__ import annotations

import os
import csv
import re
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from metagenomic_refactor.common import conda_run_command, conda_run_prefix, run_command
from metagenomic_refactor.context import get_runtime_context, update_runtime_context
from metagenomic_refactor.mag_binning import (
    MagBinningConfig,
    MagSample,
    export_legacy_binning_layout,
    run_mag_binning,
)
from metagenomic_refactor.taxonomy import exreadsID1, proc_kra1
from metagenomic_refactor.virus_analysis import (
    build_astroviridae_orf2_phylogeny_assets,
    build_enterovirus_vp1_phylogeny_assets,
    build_hpiv_coverage_assets,
    build_norovirus_gene_phylogeny_assets,
    build_rhinovirus_vp1_phylogeny_assets,
    build_seasonal_hcov_spike_phylogeny_assets,
    prepare_astroviridae_reference_annotation,
    prepare_astroviridae_sample_annotation,
    prepare_enterovirus_reference_annotation,
    prepare_enterovirus_sample_annotation,
    resolve_hepatovirus_reference,
    resolve_hiv_reference,
    prepare_norovirus_reference_annotation,
    prepare_monkeypox_reference_annotation,
    prepare_snpeff_reference_gff,
    resolve_bandavirus_reference,
    resolve_orthohantavirus_reference,
    prepare_rhinovirus_reference_annotation,
    prepare_rhinovirus_sample_annotation,
    prepare_seasonal_hcov_reference_annotation,
    prepare_seasonal_hcov_sample_annotation,
    resolve_denv_reference,
    resolve_astroviridae_reference,
    resolve_enterovirus_reference,
    resolve_hadv_reference,
    resolve_hpiv_reference,
    resolve_influenza_reference,
    resolve_norovirus_reference,
    resolve_orthoebolavirus_reference,
    resolve_rotavirus_reference,
    resolve_rhinovirus_reference,
    resolve_seasonal_hcov_reference,
    resolve_rsv_reference,
    run_astroviridae_consensus_typing,
    run_bandavirus_consensus_typing,
    run_enterovirus_consensus_typing,
    run_hepatovirus_consensus_typing,
    run_orthoebolavirus_consensus_typing,
    run_orthohantavirus_consensus_typing,
    run_rotavirus_consensus_typing,
    run_rhinovirus_consensus_typing,
    run_norovirus_consensus_typing,
)


def _database_root() -> Path:
    raw = str(os.environ.get("META_DATABASE_ROOT") or "").strip()
    if raw:
        root = Path(raw).expanduser().resolve()
        return root if root.name == "database" else root / "database"
    return Path(__file__).resolve().parent.parent / "database"


def _database_path(*parts: str) -> str:
    return str(_database_root().joinpath(*parts))


def _checkm2_database_path() -> str:
    return str(Path(os.environ.get("META_CHECKM2_DB", _database_path("checkm2", "uniref100.KO.1.dmnd"))).expanduser().resolve())


def _staramr_database_path() -> str:
    return str(Path(os.environ.get("META_STARAMR_DB", _database_path("staramr"))).expanduser().resolve())


def _safe_path_token(value) -> str:
    token = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value or "").strip())
    return token.strip("._") or "unknown"


LEGACY_CALLBACKS = {}


def register_assembly_callbacks(**kwargs):
    LEGACY_CALLBACKS.update(kwargs)


def _resolve_external_command(configured_command, fallback):
    command = str(configured_command or fallback).strip()
    return shlex.split(command) if command else [fallback]


def _is_influenza_species(species: str | None) -> bool:
    normalized = str(species or "").strip().lower()
    if not normalized:
        return False
    if "parainfluenza" in normalized or "副流感" in normalized:
        return False
    return "influenza" in normalized or "流感" in normalized


def _filter_influenza_reads(source_fastq: str, target_fastq: Path, min_length: int = 200, logf=None) -> str:
    source_path = Path(source_fastq)
    if not source_path.is_file() or source_path.stat().st_size == 0:
        return str(source_fastq)
    target_fastq.parent.mkdir(parents=True, exist_ok=True)
    run_command(
        f"seqkit seq -m {min_length} {shlex.quote(str(source_path))} > {shlex.quote(str(target_fastq))}",
        logf=logf,
    )
    if target_fastq.is_file() and target_fastq.stat().st_size > 0:
        return str(target_fastq)
    return str(source_fastq)


def _choose_reference_variant_caller(method: str, asmt: str) -> str:
    normalized = str(method or "").strip().lower()
    if normalized in {"freebayes", "clair3"}:
        return normalized
    return "clair3" if str(asmt or "").strip() == "longref" else "freebayes"


def _existing_input_path(path_value) -> str:
    candidate = str(path_value or "").strip()
    if not candidate or candidate == "0":
        return ""
    return candidate if Path(candidate).is_file() else ""


def _normalize_assembly_inputs(inf, fq1, fq2) -> tuple[str, str, str]:
    return (
        _existing_input_path(inf),
        _existing_input_path(fq1),
        _existing_input_path(fq2),
    )


def _is_sars_cov_2_species(species: str | None) -> bool:
    normalized = str(species or "").strip().lower()
    return (
        normalized in {"sars-cov-2", "新型冠状病毒", "新冠病毒", "新冠"}
        or "sars-cov-2" in normalized
    )


def _is_monkeypox_species(species: str | None) -> bool:
    normalized = str(species or "").strip().lower()
    return (
        normalized in {"monkeypox virus", "mpox", "猴痘", "猴痘病毒"}
        or "monkeypox" in normalized
        or "mpox" in normalized
        or "猴痘" in normalized
    )


def _is_rsv_species(species: str | None) -> bool:
    normalized = str(species or "").strip().lower()
    return (
        normalized in {
            "respiratory syncytial virus",
            "respiratory syncytial virus a",
            "respiratory syncytial virus b",
            "human respiratory syncytial virus",
            "human respiratory syncytial virus a",
            "human respiratory syncytial virus b",
            "orthopneumovirus hominis",
            "呼吸道合胞病毒",
            "rsv",
            "rsv-a",
            "rsv-b",
        }
        or "respiratory syncytial" in normalized
        or "合胞病毒" in normalized
    )


def _is_denv_species(species: str | None) -> bool:
    normalized = str(species or "").strip().lower()
    return (
        normalized in {
            "dengue virus",
            "dengue",
            "登革热病毒",
            "登革热",
            "denv",
            "dengue virus 1",
            "dengue virus 2",
            "dengue virus 3",
            "dengue virus 4",
            "denv1",
            "denv2",
            "denv3",
            "denv4",
            "denv-1",
            "denv-2",
            "denv-3",
            "denv-4",
        }
        or "dengue" in normalized
        or "登革热" in normalized
    )


def _is_zika_species(species: str | None) -> bool:
    normalized = str(species or "").strip().lower()
    return (
        normalized in {
            "zika virus",
            "zika",
            "zikv",
            "zikav",
            "寨卡病毒",
            "寨卡",
        }
        or "zika" in normalized
        or "寨卡" in normalized
    )


def _is_chikv_species(species: str | None) -> bool:
    normalized = str(species or "").strip().lower()
    return (
        normalized in {"chikungunya virus", "chikungunya", "chikv", "基孔肯雅病毒", "基孔肯雅"}
        or "chikungunya" in normalized
        or "基孔肯雅" in normalized
    )


def _is_hpiv_species(species: str | None) -> bool:
    normalized = str(species or "").strip().lower()
    return (
        normalized in {
            "human parainfluenza virus",
            "parainfluenza virus",
            "副流感病毒",
            "hpiv",
            "human parainfluenza virus 1",
            "human parainfluenza virus 2",
            "human parainfluenza virus 3",
            "human parainfluenza virus 4a",
            "human parainfluenza virus 4b",
            "hpiv1",
            "hpiv2",
            "hpiv3",
            "hpiv4a",
            "hpiv4b",
            "hpiv-1",
            "hpiv-2",
            "hpiv-3",
            "hpiv-4a",
            "hpiv-4b",
        }
        or "parainfluenza" in normalized
        or "副流感" in normalized
    )


def _is_hadv_species(species: str | None) -> bool:
    normalized = str(species or "").strip().lower()
    return (
        normalized in {
            "human adenovirus",
            "human mastadenovirus",
            "mastadenovirus hominis",
            "adenovirus",
            "hadv",
            "人腺病毒",
            "腺病毒",
        }
        or normalized.startswith("hadv")
        or "adenovirus" in normalized
        or "mastadenovirus" in normalized
        or "腺病毒" in normalized
    )


def _is_norovirus_species(species: str | None) -> bool:
    normalized = str(species or "").strip().lower()
    return (
        normalized in {
            "norovirus",
            "human norovirus",
            "norwalk virus",
            "诺如病毒",
            "诺瓦克病毒",
            "noro",
        }
        or "norovirus" in normalized
        or "norwalk" in normalized
        or "诺如" in normalized
    )


def _is_hepatovirus_species(species: str | None) -> bool:
    normalized = str(species or "").strip().lower()
    return (
        normalized in {
            "hepatovirus",
            "hepatitis a virus",
            "hepatitis b virus",
            "hepatitis c virus",
            "hepatitis d virus",
            "hepatitis e virus",
            "hepatovirus ahepa",
            "hav",
            "hbv",
            "hcv",
            "hdv",
            "hev",
            "甲肝病毒",
            "乙肝病毒",
            "丙肝病毒",
            "丁肝病毒",
            "戊肝病毒",
            "甲型肝炎病毒",
            "乙型肝炎病毒",
            "丙型肝炎病毒",
            "丁型肝炎病毒",
            "戊型肝炎病毒",
            "甲型肝炎",
            "乙型肝炎",
            "丙型肝炎",
            "丁型肝炎",
            "戊型肝炎",
        }
        or "hepatovirus" in normalized
        or re.search(r"\bhepatitis\s+[abcde]\s+virus\b", normalized) is not None
        or normalized in {"hav", "hbv", "hcv", "hdv", "hev"}
        or any(token in normalized for token in ("甲肝", "乙肝", "丙肝", "丁肝", "戊肝", "肝炎病毒"))
    )


def _is_hiv_species(species: str | None) -> bool:
    normalized = str(species or "").strip().lower()
    return (
        normalized in {
            "hiv",
            "hiv-1",
            "hiv1",
            "hiv-2",
            "hiv2",
            "human immunodeficiency virus",
            "human immunodeficiency virus 1",
            "human immunodeficiency virus 2",
            "艾滋病病毒",
            "艾滋病毒",
        }
        or "human immunodeficiency virus" in normalized
        or re.search(r"\bhiv(?:[-\s]?[12])?\b", normalized) is not None
        or "艾滋" in normalized
    )


def _is_rhinovirus_species(species: str | None) -> bool:
    normalized = str(species or "").strip().lower()
    return (
        normalized in {
            "human rhinovirus",
            "rhinovirus",
            "鼻病毒",
            "hrv",
            "rv",
            "human rhinovirus a",
            "human rhinovirus b",
            "human rhinovirus c",
        }
        or "rhinovirus" in normalized
        or "鼻病毒" in normalized
    )


def _is_enterovirus_species(species: str | None) -> bool:
    normalized = str(species or "").strip().lower()
    return (
        normalized not in {"", "n/a"}
        and "rhinovirus" not in normalized
        and "鼻病毒" not in normalized
        and (
            normalized in {
                "human enterovirus",
                "enterovirus",
                "human enterovirus a",
                "human enterovirus b",
                "human enterovirus c",
                "human enterovirus d",
                "enterovirus a",
                "enterovirus b",
                "enterovirus c",
                "enterovirus d",
                "肠道病毒",
            }
            or "enterovirus" in normalized
            or "肠道病毒" in normalized
        )
    )


def _is_bandavirus_species(species: str | None) -> bool:
    normalized = str(species or "").strip().lower()
    return (
        normalized not in {"", "n/a"}
        and (
            normalized in {
                "bandavirus",
                "bandavirus dabieense",
                "severe fever with thrombocytopenia syndrome virus",
                "sftsv",
                "heartland virus",
                "hunter island group virus",
                "guertu virus",
                "bhanja virus",
                "lone star virus",
                "kismaayo virus",
                "razdan virus",
                "zwiesel bat bandavirus",
                "蜱传播班达病毒",
                "班达病毒",
                "发热伴血小板减少综合征病毒",
            }
            or "bandavirus" in normalized
            or "sftsv" in normalized
            or "thrombocytopenia syndrome virus" in normalized
            or "班达病毒" in normalized
            or "血小板减少综合征病毒" in normalized
        )
    )


def _is_orthohantavirus_species(species: str | None) -> bool:
    normalized = str(species or "").strip().lower()
    if normalized in {"", "n/a"} or _is_bandavirus_species(species):
        return False
    return (
        normalized in {
            "orthohantavirus",
            "hantavirus",
            "汉坦病毒",
            "汉他病毒",
        }
        or "orthohantavirus" in normalized
        or "hantavirus" in normalized
        or "汉坦" in normalized
        or "汉他" in normalized
    )


def _is_orthoebolavirus_species(species: str | None) -> bool:
    normalized = str(species or "").strip().lower()
    if normalized in {"", "n/a"}:
        return False
    return (
        normalized in {
            "orthoebolavirus",
            "ebola virus",
            "ebolavirus",
            "ebov",
            "sudan virus",
            "sudv",
            "bundibugyo virus",
            "bdbv",
            "reston virus",
            "restv",
            "tai forest virus",
            "taï forest virus",
            "tafv",
            "bombali virus",
            "bomv",
            "埃博拉病毒",
            "正埃博拉病毒",
        }
        or "orthoebolavirus" in normalized
        or "ebolavirus" in normalized
        or "ebola" in normalized
        or "埃博拉" in normalized
    )


def _is_astroviridae_species(species: str | None) -> bool:
    normalized = str(species or "").strip().lower()
    return (
        normalized in {
            "astroviridae",
            "astrovirus",
            "human astrovirus",
            "avian astrovirus",
            "mamastrovirus",
            "avastrovirus",
            "星状病毒",
        }
        or "astrovirus" in normalized
        or "mamastrovirus" in normalized
        or "avastrovirus" in normalized
        or "星状病毒" in normalized
    )


def _is_rotavirus_species(species: str | None) -> bool:
    normalized = str(species or "").strip().lower()
    return (
        normalized in {
            "human rotavirus",
            "rotavirus",
            "轮状病毒",
            "human rotavirus a",
            "human rotavirus b",
            "human rotavirus c",
            "rotavirus a",
            "rotavirus b",
            "rotavirus c",
            "轮状病毒a",
            "轮状病毒b",
            "轮状病毒c",
            "a组轮状病毒",
            "b组轮状病毒",
            "c组轮状病毒",
            "rva",
            "rvb",
            "rvc",
        }
        or "rotavirus" in normalized
        or "轮状病毒" in normalized
    )


def _is_seasonal_hcov_species(species: str | None) -> bool:
    normalized = str(species or "").strip().lower()
    if not normalized:
        return False
    if _is_sars_cov_2_species(species):
        return False
    return (
        normalized in {
            "human coronavirus",
            "seasonal human coronavirus",
            "seasonal coronavirus",
            "人冠状病毒",
            "季节性冠状病毒",
            "human coronavirus 229e",
            "human coronavirus nl63",
            "human coronavirus oc43",
            "human coronavirus hku1",
            "hcov-229e",
            "hcov-nl63",
            "hcov-oc43",
            "hcov-hku1",
            "229e",
            "nl63",
            "oc43",
            "hku1",
        }
        or ("coronavirus" in normalized and any(token in normalized for token in ["229e", "nl63", "oc43", "hku1"]))
        or ("冠状病毒" in normalized and any(token in normalized for token in ["229e", "nl63", "oc43", "hku1"]))
    )


def _allows_auto_reference_selection(species: str | None) -> bool:
    return _is_influenza_species(species) or _is_rsv_species(species) or _is_denv_species(species) or _is_zika_species(species) or _is_chikv_species(species) or _is_hpiv_species(species) or _is_hadv_species(species) or _is_norovirus_species(species) or _is_rhinovirus_species(species) or _is_enterovirus_species(species) or _is_astroviridae_species(species) or _is_rotavirus_species(species) or _is_seasonal_hcov_species(species) or _is_bandavirus_species(species) or _is_orthohantavirus_species(species) or _is_orthoebolavirus_species(species) or _is_hepatovirus_species(species) or _is_hiv_species(species)


def _resolve_project_snpeff_paths() -> tuple[Path, Path]:
    project_root = Path(__file__).resolve().parents[1]
    snpeff_launcher = project_root / "snpEff" / "exec" / "snpeff"
    snpeff_jar = project_root / "snpEff" / "snpEff.jar"
    return snpeff_launcher, snpeff_jar


def _run_ncov_raw_vcf_annotation() -> None:
    project_root = Path(__file__).resolve().parents[1]
    script_path = project_root / "scripts" / "annotate_ncov_vcf.py"
    snpeff_path, _ = _resolve_project_snpeff_paths()
    ncov_db_dir = project_root / "database" / "virus" / "ncov"
    raw_vcf = Path("snps.raw.vcf")
    if not raw_vcf.is_file() or raw_vcf.stat().st_size == 0:
        raise FileNotFoundError("未生成 snps.raw.vcf，无法执行新冠突变注释。")
    cmd = [
        sys.executable,
        str(script_path),
        str(raw_vcf.resolve()),
        "--snpeff",
        str(snpeff_path),
        "--ncov-db-dir",
        str(ncov_db_dir),
    ]
    subprocess.run(cmd, check=True)


def _max_homopolymer_run(sequence: str | None) -> int:
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


def _is_monkeypox_low_quality_indel(ref: str | None, alt: str | None) -> bool:
    ref_text = str(ref or "").strip().upper()
    alt_text = str(alt or "").strip().upper()
    if not ref_text or not alt_text or len(ref_text) == len(alt_text):
        return False
    return max(_max_homopolymer_run(ref_text), _max_homopolymer_run(alt_text)) > 6


def _first_info_token(value: str | None, default: str = "") -> str:
    text = str(value or "").strip()
    if not text:
        return default
    return text.split(",")[0].strip() or default


def _first_info_int(value: str | None, default: int = 0) -> int:
    token = _first_info_token(value, "")
    try:
        return int(float(token))
    except ValueError:
        return default


def _first_info_float(value: str | None, default: float = 0.0) -> float:
    token = _first_info_token(value, "")
    try:
        return float(token)
    except ValueError:
        return default


def _parse_info_field(info: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in str(info or "").split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
        else:
            result[item] = ""
    return result


def _filter_monkeypox_raw_vcf_to_filt1(raw_vcf: Path, output_vcf: Path) -> None:
    if not raw_vcf.is_file():
        raise FileNotFoundError(f"未找到原始 VCF: {raw_vcf}")
    with raw_vcf.open("r", encoding="utf-8", errors="ignore") as in_handle, output_vcf.open("w", encoding="utf-8") as out_handle:
        for line in in_handle:
            if not line:
                continue
            if line.startswith("#"):
                out_handle.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            ref = fields[3]
            alt = _first_info_token(fields[4], "")
            qual = _first_info_float(fields[5], 0.0)
            info_map = _parse_info_field(fields[7])
            fmt_keys = fields[8].split(":")
            sample_values = fields[9].split(":")
            sample_map = {fmt_keys[i]: sample_values[i] if i < len(sample_values) else "" for i in range(len(fmt_keys))}
            dp = _first_info_int(info_map.get("DP"), 0)
            if dp <= 0:
                dp = _first_info_int(sample_map.get("DP"), 0)
            ao = _first_info_float(info_map.get("AO"), 0.0)
            if ao <= 0:
                ao = _first_info_float(sample_map.get("AO"), 0.0)
            ro = _first_info_float(info_map.get("RO"), 0.0)
            if ro < 0:
                ro = 0.0
            if ro <= 0:
                ro = _first_info_float(sample_map.get("RO"), 0.0)
            total_allele_depth = ao + ro
            if total_allele_depth <= 0 and dp > 0:
                total_allele_depth = float(dp)
            maf = (ao / total_allele_depth) if total_allele_depth > 0 else 0.0
            if _is_monkeypox_low_quality_indel(ref, alt):
                continue
            if qual > 10 and dp > 10 and maf > 0.1:
                out_handle.write(line)


def run_meta_viral_assembly(Pre, threads):
    runtime = get_runtime_context()
    if str(runtime.method or "").strip() != "meta":
        return
    with open("viral_assembly.log", "w") as logf:
        _run_meta_viral_assembly(Pre, threads, logf)


def _run_meta_viral_assembly(Pre, threads, logf):
    runtime = get_runtime_context()
    resources = runtime.resources
    viral_dir = Path("viral_assembly")
    viral_dir.mkdir(exist_ok=True)
    prefix = f"{Pre}_2" if Path(f"{Pre}_2.report.txt").is_file() else Pre
    report_path = Path(f"{prefix}.report.txt")
    out_path = Path(f"{prefix}.out.txt")
    read1_path = Path(f"{Pre}.R1.fastq.gz")
    read2_path = Path(f"{Pre}.R2.fastq.gz")
    if not report_path.is_file() or not out_path.is_file() or not read1_path.is_file():
        _write_viral_assembly_placeholder(
            viral_dir,
            reason="缺少病毒 reads 提取所需的 kraken2 结果或原始样本 reads，未执行病毒组装。",
        )
        return

    virus_taxids = [int(item) for item in proc_kra1(str(report_path), 10239, "D")]
    if not virus_taxids:
        _write_viral_assembly_placeholder(
            viral_dir,
            reason="未在 kraken2 分类结果中检出可提取的病毒域 reads，跳过病毒组装。",
        )
        return

    exreadsID1(virus_taxids, str(out_path), str(read1_path), str(read2_path) if read2_path.is_file() else 0)
    extracted_r1 = Path("10239.1.fastq")
    extracted_r2 = Path("10239.2.fastq")
    viral_r1 = viral_dir / "virus_reads_R1.fastq"
    viral_r2 = viral_dir / "virus_reads_R2.fastq"
    if extracted_r1.is_file():
        shutil.move(str(extracted_r1), viral_r1)
    if extracted_r2.is_file():
        shutil.move(str(extracted_r2), viral_r2)
    if not viral_r1.is_file() or viral_r1.stat().st_size == 0:
        _write_viral_assembly_placeholder(
            viral_dir,
            reason="病毒 reads 提取后为空，未生成病毒组装结果。",
        )
        return

    megahit_dir = viral_dir / "megahit_output"
    if not (megahit_dir / "final.contigs.fa").is_file():
        cmd = [
            *conda_run_prefix("mag_aux"),
            "megahit",
            "-1",
            str(viral_r1),
        ]
        if viral_r2.is_file() and viral_r2.stat().st_size > 0:
            cmd.extend(["-2", str(viral_r2)])
        cmd.extend(["-o", str(megahit_dir)])
        run_command(" ".join(shlex.quote(part) for part in cmd), logf=logf)

    contig_path = megahit_dir / "final.contigs.fa"
    if not contig_path.is_file() or contig_path.stat().st_size == 0:
        raise RuntimeError("病毒组装失败：viral_assembly/megahit_output/final.contigs.fa 未生成。")

    virsorter_out = viral_dir / "virsorter2"
    if not (virsorter_out / "final-viral-score.tsv").is_file():
        virsorter_cmd = _resolve_external_command(
            resources.virsorter2_bin if resources is not None else "",
            "virsorter",
        ) + [
            "run",
            "-w",
            str(virsorter_out),
            "-i",
            str(contig_path),
            "-j",
            str(threads),
            "all",
            "--include-groups",
            "dsDNAphage,ssDNA,RNA",
        ]
        if resources is not None and str(resources.virsorter2_db or "").strip():
            virsorter_cmd.extend(["--db-dir", str(resources.virsorter2_db)])
        run_command(" ".join(shlex.quote(part) for part in virsorter_cmd), logf=logf)

    checkv_out = viral_dir / "checkv"
    if not (checkv_out / "contamination.tsv").is_file():
        checkv_cmd = _resolve_external_command(
            resources.checkv_bin if resources is not None else "",
            "checkv",
        ) + [
            "end_to_end",
            str(contig_path),
            str(checkv_out),
            "-t",
            str(threads),
        ]
        if resources is not None and str(resources.checkv_db or "").strip():
            checkv_cmd.extend(["-d", str(resources.checkv_db)])
        run_command(" ".join(shlex.quote(part) for part in checkv_cmd), logf=logf)

    summary_df = _build_viral_contig_summary(contig_path, virsorter_out, checkv_out)
    summary_path = viral_dir / "viral_contig_summary.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)
    retained = summary_df.loc[summary_df["retained"] == "yes", "contig_id"].astype(str).tolist()
    retained_fasta = viral_dir / "viral_retained_contigs.fa"
    _write_filtered_fasta(contig_path, retained_fasta, retained)
    with (viral_dir / "viral_summary.tsv").open("w", encoding="utf-8") as handle:
        handle.write("指标\t值\n")
        handle.write(f"病毒候选contig数\t{summary_df.shape[0]}\n")
        handle.write(f"最终保留contig数\t{len(retained)}\n")
        handle.write(f"总保留长度\t{int(summary_df.loc[summary_df['retained'] == 'yes', 'contig_length'].sum())}\n")


def _write_viral_assembly_placeholder(viral_dir: Path, reason: str) -> None:
    summary_path = viral_dir / "viral_contig_summary.tsv"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        [{
            "contig_id": "-",
            "contig_length": 0,
            "virsorter2_score": "",
            "virsorter2_group": "",
            "hallmark": 0,
            "viral_genes": 0,
            "host_genes": 0,
            "checkv_quality": "",
            "completeness": "",
            "contamination": "",
            "retained": "no",
            "retention_reason": reason,
        }],
    ).to_csv(summary_path, sep="\t", index=False)
    with (viral_dir / "viral_summary.tsv").open("w", encoding="utf-8") as handle:
        handle.write("指标\t值\n")
        handle.write("病毒候选contig数\t0\n")
        handle.write("最终保留contig数\t0\n")
        handle.write(f"说明\t{reason}\n")
    (viral_dir / "viral_retained_contigs.fa").write_text("", encoding="utf-8")


def _read_contig_lengths(contig_path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    with contig_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            lengths[str(record.id)] = len(record.seq)
    return lengths


def _load_virsorter2_scores(virsorter_out: Path) -> pd.DataFrame:
    score_path = virsorter_out / "final-viral-score.tsv"
    if not score_path.is_file() or score_path.stat().st_size == 0:
        return pd.DataFrame(columns=["contig_id", "virsorter2_score", "virsorter2_group", "hallmark"])
    score_df = pd.read_table(score_path)
    if score_df.empty:
        return pd.DataFrame(columns=["contig_id", "virsorter2_score", "virsorter2_group", "hallmark"])
    rename_map = {}
    if "seqname" in score_df.columns:
        rename_map["seqname"] = "contig_id"
    if "max_score" in score_df.columns:
        rename_map["max_score"] = "virsorter2_score"
    if "max_score_group" in score_df.columns:
        rename_map["max_score_group"] = "virsorter2_group"
    score_df = score_df.rename(columns=rename_map)
    score_df["contig_id"] = score_df["contig_id"].astype(str).str.split("|").str[0]
    if "hallmark" not in score_df.columns:
        score_df["hallmark"] = 0
    keep_columns = ["contig_id", "virsorter2_score", "virsorter2_group", "hallmark"]
    return score_df[[column for column in keep_columns if column in score_df.columns]].drop_duplicates("contig_id")


def _load_checkv_contamination(checkv_out: Path) -> pd.DataFrame:
    contamination_path = checkv_out / "contamination.tsv"
    if not contamination_path.is_file() or contamination_path.stat().st_size == 0:
        return pd.DataFrame(columns=["contig_id", "total_genes", "viral_genes", "host_genes", "provirus"])
    checkv_df = pd.read_table(contamination_path)
    if "contig_id" not in checkv_df.columns:
        return pd.DataFrame(columns=["contig_id", "total_genes", "viral_genes", "host_genes", "provirus"])
    keep = [column for column in ["contig_id", "total_genes", "viral_genes", "host_genes", "provirus"] if column in checkv_df.columns]
    return checkv_df[keep].drop_duplicates("contig_id")


def _load_checkv_quality(checkv_out: Path) -> pd.DataFrame:
    quality_path = checkv_out / "quality_summary.tsv"
    if not quality_path.is_file() or quality_path.stat().st_size == 0:
        return pd.DataFrame(columns=["contig_id", "checkv_quality", "completeness", "contamination"])
    quality_df = pd.read_table(quality_path)
    rename_map = {}
    if "contig_id" not in quality_df.columns and "genome" in quality_df.columns:
        rename_map["genome"] = "contig_id"
    if "checkv_quality" not in quality_df.columns and "miuvig_quality" in quality_df.columns:
        rename_map["miuvig_quality"] = "checkv_quality"
    quality_df = quality_df.rename(columns=rename_map)
    for missing in ["checkv_quality", "completeness", "contamination"]:
        if missing not in quality_df.columns:
            quality_df[missing] = ""
    keep = ["contig_id", "checkv_quality", "completeness", "contamination"]
    return quality_df[[column for column in keep if column in quality_df.columns]].drop_duplicates("contig_id")


def _viral_retention_decision(row: pd.Series) -> tuple[str, str]:
    hallmark = int(float(row.get("hallmark") or 0))
    viral_genes = int(float(row.get("viral_genes") or 0))
    host_genes = int(float(row.get("host_genes") or 0))
    length = int(float(row.get("contig_length") or 0))
    quality = str(row.get("checkv_quality") or "").strip()
    if quality in {"Complete", "High-quality", "Medium-quality"}:
        return "yes", f"CheckV={quality}"
    if hallmark > 2:
        return "yes", "VirSorter2 hallmark>2"
    if viral_genes > 0:
        return "yes", "CheckV viral_genes>0"
    if viral_genes == 0 and host_genes == 0:
        return "yes", "viral_genes=0 且 host_genes=0"
    if viral_genes == 0 and host_genes == 1 and length >= 10000:
        return "yes", "host_genes=1 且 contig>=10kb"
    return "no", "未通过 VirSorter2/CheckV 过滤规则"


def _build_viral_contig_summary(contig_path: Path, virsorter_out: Path, checkv_out: Path) -> pd.DataFrame:
    contig_lengths = _read_contig_lengths(contig_path)
    base_df = pd.DataFrame(
        [{"contig_id": contig_id, "contig_length": length} for contig_id, length in contig_lengths.items()],
    )
    if base_df.empty:
        return pd.DataFrame(columns=[
            "contig_id", "contig_length", "virsorter2_score", "virsorter2_group", "hallmark",
            "viral_genes", "host_genes", "checkv_quality", "completeness", "contamination",
            "retained", "retention_reason",
        ])
    merged = base_df.merge(_load_virsorter2_scores(virsorter_out), on="contig_id", how="left")
    merged = merged.merge(_load_checkv_contamination(checkv_out), on="contig_id", how="left")
    merged = merged.merge(_load_checkv_quality(checkv_out), on="contig_id", how="left")
    for column in ["virsorter2_score", "hallmark", "viral_genes", "host_genes", "completeness", "contamination"]:
        if column not in merged.columns:
            merged[column] = 0
    for column in ["virsorter2_group", "checkv_quality"]:
        if column not in merged.columns:
            merged[column] = ""
    merged["hallmark"] = pd.to_numeric(merged["hallmark"], errors="coerce").fillna(0).astype(int)
    merged["viral_genes"] = pd.to_numeric(merged["viral_genes"], errors="coerce").fillna(0).astype(int)
    merged["host_genes"] = pd.to_numeric(merged["host_genes"], errors="coerce").fillna(0).astype(int)
    decisions = merged.apply(_viral_retention_decision, axis=1)
    merged["retained"] = decisions.apply(lambda item: item[0])
    merged["retention_reason"] = decisions.apply(lambda item: item[1])
    return merged[[
        "contig_id", "contig_length", "virsorter2_score", "virsorter2_group", "hallmark",
        "viral_genes", "host_genes", "checkv_quality", "completeness", "contamination",
        "retained", "retention_reason",
    ]].sort_values(["retained", "contig_length"], ascending=[False, False])


def _write_filtered_fasta(source_fasta: Path, output_fasta: Path, keep_ids: list[str]) -> None:
    keep_set = {str(item).strip() for item in keep_ids if str(item).strip()}
    with source_fasta.open("r", encoding="utf-8", errors="ignore") as in_handle, output_fasta.open("w", encoding="utf-8") as out_handle:
        for record in SeqIO.parse(in_handle, "fasta"):
            if str(record.id) in keep_set:
                SeqIO.write(record, out_handle, "fasta")


def _resolve_hiv_hxb2_reference() -> Path:
    return Path(__file__).resolve().parents[1] / "database" / "virus" / "HIV" / "HXB2_K03455.fasta"


def _build_hiv_hxb2_consensus(inf, fq1, fq2, threads, Pre, pts, pst, method, asmt, logf=None) -> str:
    runtime = get_runtime_context()
    if runtime.analysis_target != "virus" or not _is_hiv_species(runtime.species):
        return ""
    hxb2_ref = _resolve_hiv_hxb2_reference()
    if not hxb2_ref.is_file():
        if logf is not None:
            logf.write(f"[HIV_HXB2_CONSENSUS_SKIP] missing reference: {hxb2_ref}\n")
            logf.flush()
        return ""

    output_path = Path(f"{Pre}.hxb2.consensus.fasta").resolve()
    workdir = Path(f"{Pre}_hiv_hxb2_consensus_work")
    inf_abs = str(Path(inf).resolve()) if str(inf or "").strip() else ""
    fq1_abs = str(Path(fq1).resolve()) if str(fq1 or "").strip() else ""
    fq2_abs = str(Path(fq2).resolve()) if str(fq2 or "").strip() else ""
    if workdir.exists():
        shutil.rmtree(workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    previous_state = {
        "ref": runtime.ref,
        "gtf": runtime.gtf,
        "species": runtime.species,
    }
    original_cwd = Path.cwd()
    try:
        update_runtime_context(ref=str(hxb2_ref.resolve()), gtf="nogtf")
        os.chdir(workdir)
        reassm_fun(inf_abs, fq1_abs, fq2_abs, threads, Pre, pts, pst, method, asmt, logf, str(output_path))
    except Exception as exc:
        if logf is not None:
            logf.write(f"[HIV_HXB2_CONSENSUS_FAIL] {exc}\n")
            logf.flush()
        return ""
    finally:
        os.chdir(original_cwd)
        update_runtime_context(**previous_state)

    if output_path.is_file() and output_path.stat().st_size > 0:
        if logf is not None:
            logf.write(f"[HIV_HXB2_CONSENSUS_READY] {output_path}\n")
            logf.flush()
        return str(output_path)
    if logf is not None:
        logf.write("[HIV_HXB2_CONSENSUS_FAIL] output fasta missing or empty\n")
        logf.flush()
    return ""


def asb_func(inf, fq1, fq2, threads, Pre, lelID, pts, pst, method, asmt="longasm", ref="noref", gtf="nogtf", tryref=False):
    runtime = get_runtime_context()
    with open("asb.log", "w") as f:
        outputfa = f"{Pre}.consensus.fasta"
        polifa = f"{Pre}.polish.fasta"
        hiv_reference_method = ""
        fq1, fq2, _ = prepare_assembly_inputs(inf, fq1, fq2, threads, Pre, lelID, method)
        inf, fq1, fq2 = _normalize_assembly_inputs(inf, fq1, fq2)

        if not os.path.isfile(f"assm_{method}_{pts}_ok"):
            methodlist = method.split(",")
            if len(methodlist) == 1:
                run_method = methodlist[0]
                if run_method in ["flye", "canu", "unicycler", "masurca", "spades", "raven", "wtdbg2", "miniasm", "meta"]:
                    denovo_asb(inf, fq1, fq2, threads, Pre, pts, pst, run_method, asmt, f, outputfa)
                else:
                    if ref != "noref" or _allows_auto_reference_selection(runtime.species):
                        reassm_fun(inf, fq1, fq2, threads, Pre, pts, pst, run_method, asmt, f, outputfa)
                        hiv_reference_method = run_method
                    else:
                        print("有参组装未提供参考基因组")
                        sys.exit()
            elif len(methodlist) == 2:
                denovomethod = methodlist[0]
                reasmmethod = methodlist[1]
                denovo_asb(inf, fq1, fq2, threads, Pre, pts, pst, denovomethod, asmt, f, outputfa)
                if ref != "noref" or _allows_auto_reference_selection(runtime.species):
                    reassm_fun(inf, fq1, fq2, threads, Pre, pts, pst, reasmmethod, asmt, f, outputfa="noforce")
                    hiv_reference_method = reasmmethod
                else:
                    print("有参组装未提供参考基因组")
                    sys.exit()
    if method != "meta":
        if os.path.isfile(outputfa):
            finalfa = f"{Pre}.final.fasta"
            with open("Anno1.log", "w") as f:
                if runtime.analysis_target != 'virus':
                    if os.path.isfile(polifa):
                        renamefa(polifa, finalfa)
                    else:
                        renamefa(outputfa, finalfa)
                else:
                    subprocess.run(f'cp {outputfa} {finalfa}',shell=True)
                
                map_assembly_reads(finalfa, inf, fq1, fq2, threads, Pre, asmt, runtime.long_type, f)
                if hiv_reference_method and runtime.analysis_target == 'virus' and _is_hiv_species(runtime.species):
                    _build_hiv_hxb2_consensus(inf, fq1, fq2, threads, Pre, pts, pst, hiv_reference_method, asmt, f)
                if runtime.analysis_target != 'virus':
                    build_assembly_info(finalfa, Pre, method)
                    flyedb, plasmidlist = enhance_plasmid_results(finalfa, Pre)
                    if os.path.isfile(finalfa):
                        write_finalfasta_stats(finalfa, Pre)
                    annotate_with_prokka(finalfa, Pre, threads)
                    export_plasmid_gbk_and_cgview(Pre, plasmidlist)
                    run_checkm_and_write_summary(Pre, finalfa, threads, runtime.species, flyedb)
                    run_genovi_summary(Pre)
        else:
            raise Exception(f"{method}组装报错可能原因:\n1.数据量过少，提高数据量后继续分析\n2.样本中含有一定量的杂菌污染")
    time.sleep(0.5)


def Annotate_func(Pre, threads):
    runtime = get_runtime_context()
    with open("Anno.log", "w") as f1:
        finalfa = f"{Pre}.final.fasta"
        renamefa(finalfa, finalfa)
        if runtime.analysis_target != 'virus':
            if not os.path.isfile(f"{Pre}_prokka/{Pre}.tsv"):
                if os.path.isfile(finalfa):
                    write_finalfasta_stats(finalfa, Pre)
            
                annotate_with_prokka(finalfa, Pre, threads, "Anno.log")
                flyedb = build_fasta_assembly_info(finalfa, Pre, runtime.method)
                flyedb, plasmidlist = enhance_plasmid_results(finalfa, Pre, "Anno.log")
                export_plasmid_gbk_and_cgview(Pre, plasmidlist)
                run_checkm_and_write_summary(Pre, finalfa, threads, runtime.species, flyedb)
                run_genovi_summary(Pre)
            write_gene_summaries(Pre)
        time.sleep(5)


def polish_func(Pre, ptimes, threads, psoft="medaka", input_fastq: str | None = None):
    print(f"开始抛光 抛光软件: {psoft} 抛光次数: {ptimes}")
    ptimes = str(ptimes)
    medaka_cmd = conda_run_command("longread_aux", "medaka_consensus")
    try:
        with open("polish.log", "w") as f:
            for i in range(1, int(ptimes) + 1):
                input_file = input_fastq if i == 1 and input_fastq else (f"{Pre}.final.fastq" if i == 1 else f"medaka_output{i-1}/consensus.fasta")
                output_dir = f"medaka_output{i}"
                output_file = f"{output_dir}/consensus.fasta"
                cmd = f"{medaka_cmd} -i {input_file} -d {Pre}.consensus.fasta -o {output_dir} -t {threads} > medaka.log"
                subprocess.run(cmd, shell=True, stdout=f, stderr=f)
                if i == int(ptimes):
                    subprocess.run(f"seqkit seq -w0 {output_file} > {Pre}.polish.fasta", shell=True, stdout=f, stderr=f)
        if os.path.getsize(f"{Pre}.polish.fasta") == 0:
            subprocess.run(f"seqkit seq -w0 {Pre}.consensus.fasta > {Pre}.polish.fasta", shell=True)
    except Exception as e:
        print(f"抛光过程出现错误: {e}")
    print("抛光结束")


def _export_influenza_reference_outputs(
    pre: str,
    threads: int,
    outputfa: str,
    variant_caller: str,
    long_reads: str | None = None,
    polish_times: str | int = 1,
    polish_soft: str = "medaka",
    logf=None,
):
    flu_root = Path("wf_flu")
    consensus_dir = flu_root / "consensus"
    typing_dir = flu_root / "typing"
    consensus_dir.mkdir(parents=True, exist_ok=True)
    typing_dir.mkdir(parents=True, exist_ok=True)

    try:
        if Path("ref.mapping.bam").is_file():
            run_command(
                f"samtools depth -aa -d 0 ref.mapping.bam > {shlex.quote(str(flu_root / 'coverage.depth.tsv'))}",
                logf=logf,
                check=False,
            )
    except Exception:
        pass

    final_consensus = Path(outputfa)
    if variant_caller == "clair3" and long_reads:
        try:
            polish_func(pre, polish_times or 1, threads, polish_soft, input_fastq=long_reads)
            polished = Path(f"{pre}.polish.fasta")
            if polished.is_file() and polished.stat().st_size > 0:
                final_consensus = polished
        except Exception as exc:
            if logf is not None:
                logf.write(f"[WF_FLU_POLISH_FAIL] {exc}\n")
                logf.flush()

    if final_consensus.is_file() and final_consensus.stat().st_size > 0:
        subprocess.run(
            f"cp {shlex.quote(str(final_consensus))} {shlex.quote(str(consensus_dir / 'consensus.fasta'))}",
            shell=True,
            stdout=logf,
            stderr=logf,
        )
        subprocess.run(
            f"cp {shlex.quote(str(final_consensus))} {shlex.quote(str(outputfa))}",
            shell=True,
            stdout=logf,
            stderr=logf,
        )

    abricate_bin = shutil.which("abricate")
    if abricate_bin and Path(outputfa).is_file():
        subprocess.run(
            f"{abricate_bin} --db insaflu --threads {threads} {shlex.quote(str(outputfa))} > {shlex.quote(str(typing_dir / 'insaflu.typing.tsv'))}",
            shell=True,
            stdout=logf,
            stderr=logf,
        )


def wait_for_file(filepath, cinterval=2):
    while not os.path.exists(filepath):
        time.sleep(1)
    lastsize = 1
    while True:
        current_size = os.path.getsize(filepath)
        if current_size == lastsize:
            break
        lastsize = current_size
        time.sleep(cinterval)
    return True


def rebinning():
    infpath = "BASALT_out/meta_drep_out/dereplicated_genomes/"
    outpath = "BASALT_out/meta_drep_out/binning_genomes/"
    if not os.path.isdir(outpath):
        os.makedirs(outpath)
    else:
        for stale in os.listdir(outpath):
            stale_path = os.path.join(outpath, stale)
            if os.path.isdir(stale_path):
                shutil.rmtree(stale_path)
            else:
                os.remove(stale_path)
    falist = [i for i in os.listdir(infpath) if i.endswith(".fa")]
    n = 1
    open("binning_name.tsv", "w").write("oldname\tnewname\n")
    if falist:
        for MAG in falist:
            MAG = f"{infpath}/{MAG}"
            oldname = MAG.split("/")[-1]
            oldname = re.sub(r"\.(fa|fasta|fna)(\.gz)?$", "", oldname)
            outMAG = f"{outpath}/MAG_{n}.fa"
            open("binning_name.tsv", "a").write(f"{oldname}\tMAG_{n}\n")
            open(outMAG, "w").write("")
            with open(MAG) as f:
                m = 1
                for line in f:
                    line = line.strip()
                    if line.startswith(">"):
                        open(outMAG, "a").write(f">{m}\n")
                        m += 1
                    else:
                        open(outMAG, "a").write(f"{line}\n")
            n += 1


def combinebin(refinedir, ofa):
    open(ofa, "w").write("")
    list1 = [i for i in os.listdir(refinedir) if i.endswith("fa")]
    for i in list1:
        filen = f"{refinedir}/{i}"
        newname = i.replace(".fa", "")
        with open(filen) as f:
            for line in f:
                if line.startswith(">"):
                    tmp_contig = line.strip().replace(">", "")
                    open(ofa, "a").write(f">{newname}_{tmp_contig}\n")
                else:
                    open(ofa, "a").write(line)


def bingtdbtk_fun():
    inf = "gtdbtk_out/gtdbtk.bac120.summary.tsv"
    with open("gtdbtk.log", "w") as f:
        if not os.path.isfile(inf):
            subprocess.run(conda_run_command("mag_aux", "gtdbtk classify_wf --genome_dir BASALT_out/meta_drep_out/dereplicated_genomes --out_dir gtdbtk_out -x .fa --cpus 10 --force"), shell=True, stdout=f, stderr=f)


def bincheckm2_fun():
    inf = "bin_checkm2out/quality_report.tsv"
    checkm2_db = _checkm2_database_path()
    with open("bincheckm2.log", "w") as f:
        if not os.path.isfile(inf):
            subprocess.run(conda_run_command("cm210", f"checkm2 predict --thread 10 --input BASALT_out/meta_drep_out/binning_genomes/ --output-directory bin_checkm2out -x .fa --database_path {shlex.quote(checkm2_db)}"), shell=True, stdout=f, stderr=f)


def binvfdrdb():
    inf1 = "bin_vfdb.tsv"
    inf2 = "bin_card.tsv"
    inf3 = "binning_rgi_new.txt"
    inf4 = "staramr_result/plasmidfinder.tsv"
    staramr_db = shlex.quote(_staramr_database_path())
    with open("binning.log", "w") as f:
        if not os.path.isfile(inf1):
            subprocess.run("abricate BASALT_out/meta_drep_out/binning_genomes/MAG_*.fa --db vfdb > bin_vfdb.tsv ", shell=True, stdout=f, stderr=f)
        if not os.path.isfile(inf2):
            subprocess.run("abricate BASALT_out/meta_drep_out/binning_genomes/MAG_*.fa --db card > bin_card.tsv ", shell=True, stdout=f, stderr=f)
        if not os.path.isfile(inf3):
            subprocess.run(conda_run_command("amr_aux", "rgi main -i tmp_combine.fa --clean --include_loose -o binning_rgi_new -n 20 -g PYRODIGAL --low_quality"), shell=True)
            subprocess.run("rm -r binning_rgi_new.json", shell=True)
        if not os.path.isfile(inf4):
            subprocess.run(f"staramr search -d {staramr_db} BASALT_out/meta_drep_out/binning_genomes/*.fa -o staramr_result -n 10", shell=True, stdout=f, stderr=f)


def meta_plasmid(Pre):
    staramr_db = shlex.quote(_staramr_database_path())
    with open("plasmid.log", "w") as f:
        if not os.path.isfile(f"{Pre}_plaspredict.tsv") or not os.path.isfile("staramr_result/plasmidfinder.tsv"):
            subprocess.run(conda_run_command("plasflow", f"PlasFlow.py --input tmp_combine.fa --output {Pre}_plaspredict.tsv"), shell=True, stdout=f, stderr=f)
            subprocess.run(f"staramr search -d {staramr_db} BASALT_out/meta_drep_out/binning_genomes/*.fa -o staramr_result -n 10", shell=True, stdout=f, stderr=f)
        plasmiddb = pd.read_table("staramr_result/plasmidfinder.tsv")
        rawindexname = plasmiddb.index.name
        if plasmiddb.shape[0] > 1:
            plasmiddb = plasmiddb.groupby("Contig").apply(_join_plasmid)
            plasmiddb.index.name = rawindexname
        plasmiddb["contig_name"] = plasmiddb.apply(lambda x: f"{x['Isolate ID']}_{x['Contig']}", axis=1)
        plasmiddb = plasmiddb[["contig_name", "Plasmid"]]
        plasflowdb = pd.read_table(f"{Pre}_plaspredict.tsv")
        plasflowdb = plasflowdb[["contig_name", "label"]]
        plasdb = plasflowdb.merge(plasmiddb, on="contig_name", how="left").fillna("-")
        plasdb.to_csv(f"{Pre}_meta_plaspredict.tsv", sep="\t", index=False)


def meta_tpm():
    with open("coverm.log", "w") as f:
        if not os.path.isfile("meta_tpm.tsv"):
            subprocess.run(conda_run_command("mag_aux", "coverm genome -1 2.1.fastq -2 2.2.fastq -d  BASALT_out/meta_drep_out/binning_genomes/ -x .fa --min-read-percent-identity 95 --min-read-aligned-percent 75 -m tpm -o meta_tpm.tsv -t 10"), shell=True, stderr=f, stdout=f)
            subprocess.run(conda_run_command("mag_aux", "coverm contig -1 2.1.fastq -2 2.2.fastq -r tmp_combine.fa --min-read-percent-identity 95 --min-read-aligned-percent 75 -m metabat -o meta_contig_metabat.tsv -t 10"), shell=True, stderr=f, stdout=f)


def binning_result(Pre):
    rebinning()
    combinebin("BASALT_out/meta_drep_out/binning_genomes", "tmp_combine.fa")
    bingtdbtk_fun()
    bincheckm2_fun()
    binvfdrdb()
    meta_plasmid(Pre)
    meta_tpm()
    vfdb = pd.read_table("bin_vfdb.tsv")
    argdict = {"ARG": [], "contig_name": [], "Name": [], "AR Gene(abricate)": [], "AR Gene(rgi)": [], "AR Gene(resfinder)": []}
    vfdb["Name"] = vfdb["#FILE"].str.split("/").str[-1].str.split(".").str[0]
    vfdb["contig_name"] = vfdb.apply(lambda x: f"{x['Name']}_{x['SEQUENCE']}", axis=1)
    vfdb = vfdb[["contig_name", "GENE"]]
    vfdb.rename(columns={"GENE": "VF Gene"}, inplace=True)
    carddb = pd.read_table("bin_card.tsv")
    carddb["tmpgene"] = carddb["GENE"].str.lower()
    rgidb = pd.read_table("binning_rgi_new.txt")
    rgidb = rgidb[rgidb["Cut_Off"].isin(["Strict", "Perfect"])]
    rgidb["tmpgene"] = rgidb["Best_Hit_ARO"].str.lower()
    rgidb["contig_name"] = rgidb["Contig"]
    resdb = pd.read_table("staramr_result/resfinder.tsv")
    resdb["tmpgene"] = resdb["Gene"].str.lower()
    resdb["contig_name"] = resdb.apply(lambda x: f"{x['Isolate ID']}_{x['Contig']}", axis=1)
    carddb["Name"] = carddb["#FILE"].str.split("/").str[-1].str.split(".").str[0]
    carddb["contig_name"] = carddb.apply(lambda x: f"{x['Name']}_{x['SEQUENCE']}", axis=1)
    carddb = carddb[["contig_name", "GENE"]]
    carddb["tmpgene"] = carddb["GENE"].str.lower()
    carddb.rename(columns={"GENE": "AR Gene"}, inplace=True)
    arglist = list(set(resdb["tmpgene"].tolist() + rgidb["tmpgene"].tolist() + carddb["tmpgene"].tolist()))
    for argene in arglist:
        argdict["ARG"].append(argene)
        contiglist = []
        if argene in carddb["tmpgene"].tolist():
            argdict["AR Gene(abricate)"].append("+")
            contiglist.append(carddb.loc[carddb["tmpgene"] == argene]["contig_name"].tolist()[0])
        else:
            argdict["AR Gene(abricate)"].append("-")
        if argene in resdb["tmpgene"].tolist():
            argdict["AR Gene(resfinder)"].append("+")
            contiglist.append(resdb.loc[resdb["tmpgene"] == argene]["contig_name"].tolist()[0])
        else:
            argdict["AR Gene(resfinder)"].append("-")
        if argene in rgidb["tmpgene"].tolist():
            argdict["AR Gene(rgi)"].append("+")
            contiglist.append(rgidb.loc[rgidb["tmpgene"] == argene]["contig_name"].tolist()[0])
        else:
            argdict["AR Gene(rgi)"].append("-")
        if contiglist:
            tmpName = contiglist[0]
            argdict["contig_name"].append(tmpName)
            argdict["Name"].append("_".join(tmpName.split("_")[:2]))
        else:
            argdict["contig_name"].append("-")
    argdb = pd.DataFrame(argdict)
    Alldb = pd.read_table(f"{Pre}_meta_plaspredict.tsv")
    Alldb = Alldb.merge(vfdb, on="contig_name", how="left").merge(argdb, on="contig_name", how="left").fillna("-")
    Alldb["Name"] = Alldb["contig_name"].str.split("_").str[:2].str.join("_")
    gtdbdb = pd.read_table("gtdbtk_out/gtdbtk.bac120.summary.tsv")
    gtdbdb["Name"] = gtdbdb["user_genome"]
    lvlist = ["D", "P", "C", "O", "F", "G", "S"]
    for lv in lvlist:
        gtdbdb[lv] = gtdbdb["classification"].str.split(";").str[lvlist.index(lv)].str.split("__").str[1]
    gtdbdb = gtdbdb[["Name", "D", "P", "C", "O", "F", "G", "S"]]
    binnamedb = pd.read_table("binning_name.tsv")
    gtdbdb = gtdbdb.merge(binnamedb, left_on="Name", right_on="oldname")[["newname", "D", "P", "C", "O", "F", "G", "S"]].rename(columns={"newname": "Name"})
    tpmdb = pd.read_table("meta_tpm.tsv")
    tpmdb.columns = ["Name", Pre]
    Alldb = Alldb.merge(gtdbdb, on="Name", how="left").merge(tpmdb, on="Name", how="left")
    Alldb.to_csv("meta_plas_vf_card.tsv", sep="\t", index=False)


def run_meta_mag_binning(fq1, fq2, threads, Pre, log_handle):
    sample = MagSample(
        sample=Pre,
        contigs=os.path.abspath("megahit_output/final.contigs.fa"),
        fastq1=os.path.abspath(fq1),
        fastq2=os.path.abspath(fq2) if fq2 else None,
    )
    config = MagBinningConfig(
        outdir=os.path.abspath("codex_mag_binning"),
        threads=threads,
        min_contig_len=1500,
        score_threshold=0.1,
        force=False,
    )
    summary_path = run_mag_binning([sample], config)
    print(f"MAG分箱完成，汇总文件: {summary_path}", file=log_handle)
    log_handle.flush()
    export_legacy_binning_layout(
        os.path.join(config.outdir, sample.sample),
        os.path.abspath("BASALT_out"),
    )


def denovo_asb(inf, fq1, fq2, threads, Pre, pts, pst, method, asmt, f, outputfa):
    runtime = get_runtime_context()
    assembly_long_type = runtime.long_type
    assembly_genome_len = runtime.genome_len
    if asmt == "longasm":
        if method == "flye":
            if assembly_long_type == "Nanopore":
                readsq = float(os.popen(f"nanoq -i {inf} -s -t 5 -vvv|grep 'Mean read quality'|cut -d ':' -f2").read().strip())
                per = 100 - 10 ** (float(readsq) / -10) * 100
                if per > 95:
                    subprocess.run(f"flye --nano-hq {inf} -g {assembly_genome_len} -o flye_output -t {threads} -i 3", shell=True, stdout=f, stderr=f)
                else:
                    subprocess.run(f"flye --nano-raw {inf} -g {assembly_genome_len} -o flye_output -t {threads} -i 3", shell=True, stdout=f, stderr=f)
            elif assembly_long_type == "PacBio_CLR":
                subprocess.run(f"flye --pacbio-raw {inf} -g {assembly_genome_len} -o flye_output -t {threads} -i 3", shell=True, stdout=f, stderr=f)
            elif assembly_long_type == "PacBio_CCS":
                subprocess.run(f"flye --pacbio-hifi {inf} -g {assembly_genome_len} -o flye_output -t {threads} -i 3", shell=True, stdout=f, stderr=f)
            subprocess.run(f"cp flye_output/assembly.fasta {outputfa}", shell=True)
            flyedb = pd.read_table("flye_output/assembly_info.txt", names=["seq_name", "length", "cov", "circ", "repeat", "mult", "alt_group", "graph_path"], skiprows=1)
            flyedb.rename(columns={"seq_name": "序列名称", "length": "序列长度", "cov": "平均深度", "circ": "是否成环"}, inplace=True)
            flyedb[["序列名称", "序列长度", "平均深度", "是否成环"]].to_csv("flye_output/assembly_info.txt", sep="\t", index=False)
        elif method == "miniasm":
            if assembly_long_type == "Nanopore":
                subprocess.run(f"minimap2 -x ava-ont -t {threads} {inf} {inf} | gzip -1 > {Pre}_reads.paf.gz", shell=True, stdout=f, stderr=f)
            elif assembly_long_type == "PacBio_CLR":
                subprocess.run(f"minimap2 -x map-pb -t {threads} {inf} {inf} | gzip -1 > {Pre}_reads.paf.gz", shell=True, stdout=f, stderr=f)
            elif assembly_long_type == "PacBio_CCS":
                subprocess.run(f"minimap2 -x map-hifi -t {threads} {inf} {inf} | gzip -1 > {Pre}_reads.paf.gz", shell=True, stdout=f, stderr=f)
            subprocess.run(f"miniasm -f {inf} {Pre}_reads.paf.gz > {Pre}_reads.gfa", shell=True, stdout=f, stderr=f)
            subprocess.run(f"gfatools gfa2fa {Pre}_reads.gfa > {Pre}_miniasm.fa", shell=True, stdout=f, stderr=f)
            subprocess.run(f"cp {Pre}_miniasm.fa {outputfa}", shell=True)
        elif method == "wtdbg2":
            subprocess.run(f"perl ~/biosoft/wtdbg2/wtdbg2.pl -t {threads} -x ont -g {assembly_genome_len} -o wtdbg2 {inf}", shell=True, stdout=f, stderr=f)
            subprocess.run(f"cp wtdbg2.cns.fa {outputfa}", shell=True)
        elif method == "canu":
            if assembly_long_type == "Nanopore":
                subprocess.run(f"time canu -d canu -p canu genomeSize={assembly_genome_len} maxThreads={threads} -nanopore-raw {inf} >canu.log", shell=True, stdout=f, stderr=f)
            else:
                subprocess.run(f"time canu -d canu -p canu genomeSize={assembly_genome_len} maxThreads={threads} -pacbio-raw {inf} >canu.log", shell=True, stdout=f, stderr=f)
            subprocess.run(f"cp canu/canu.contigs.fasta {outputfa}", shell=True)
        elif method == "unicycler":
            subprocess.run(f"unicycler -t {threads} -l {inf} -o unicycler", shell=True, stdout=f, stderr=f)
            subprocess.run(f"cp unicycler/assembly.fasta {outputfa}", shell=True)
        elif method == "raven":
            subprocess.run(f"raven {inf} -t {threads} > raven.fasta", shell=True, stdout=f, stderr=f)
            subprocess.run(f"cp raven.fasta {outputfa}", shell=True)
        else:
            print(f"请确认传入参数method是否正确，可选[flye,canu,wtdbg2,unicycler,miniasm],您输入的为：{method}")
        if os.path.isfile(outputfa):
            polish_func(Pre, pts, threads, pst)
    elif asmt == "shortasm":
        if method == "spades":
            if fq2:
                subprocess.run(f"spades.py --pe1-1 {fq1} --pe1-2 {fq2} -t {threads} -o spades_output --cov-cutoff 8 --isolate", shell=True, stdout=f, stderr=f)
            else:
                subprocess.run(f"spades.py -s {fq1} -t {threads} -o spades_output --isolate --cov-cutoff 8", shell=True, stdout=f, stderr=f)
            if os.path.isfile("spades_output/contigs.fasta"):
                subprocess.run(f"cp spades_output/contigs.fasta {outputfa}", shell=True, stdout=f, stderr=f)
        elif method == "masurca":
            if fq2:
                subprocess.run(f"masurca -i {fq1},{fq2} -t {threads} -o masurca_output", shell=True, stdout=f, stderr=f)
            else:
                subprocess.run(f"masurca -i {fq1} -t {threads} -o masurca_output", shell=True, stdout=f, stderr=f)
            CAdir = os.popen("ls -d CA").read().split("\n")[0]
            if os.path.isfile(f"{CAdir}/primary.genome.scf.fasta"):
                subprocess.run(f"cp {CAdir}/primary.genome.scf.fasta {outputfa}", shell=True)
            else:
                if os.path.isfile(f"{CAdir}/scaffolds.ref.fa"):
                    subprocess.run(f"cp {CAdir}/scaffolds.ref.fa {outputfa}", shell=True)
                else:
                    print("组装失败")
                    sys.exit()
        elif method == "meta":
            if not os.path.isfile("megahit_output/final.contigs.fa"):
                if fq2:
                    subprocess.run(conda_run_command("mag_aux", f"megahit -1 {fq1} -2 {fq2} -t 10 -o megahit_output"), shell=True, stdout=f, stderr=f)
                else:
                    subprocess.run(conda_run_command("mag_aux", f"megahit -1 {fq1} -t 10 -o megahit_output"), shell=True, stdout=f, stderr=f)
            if os.path.isfile("megahit_output/final.contigs.fa") and os.path.getsize("megahit_output/final.contigs.fa") != 0:
                if not os.path.isdir("BASALT_out"):
                    os.makedirs("BASALT_out")
                run_meta_mag_binning(fq1, fq2, threads, Pre, f)
                binning_result(Pre)
            else:
                print("宏基因组组装失败")
                sys.exit()
    elif asmt == "shortlongasm":
        if method == "unicycler":
            if fq2:
                subprocess.run(f"unicycler -1 {fq1} -2 {fq2} -l {inf} -t {threads} -o unicycler", shell=True, stdout=f, stderr=f)
            elif fq1:
                subprocess.run(f"unicycler -1 {fq1} -l {inf} -t {threads} -o unicycler", shell=True, stdout=f, stderr=f)
            if os.path.isfile("unicycler/assembly.fasta"):
                subprocess.run(f"cp unicycler/assembly.fasta {outputfa}", shell=True)
        elif method == "masurca":
            if fq2:
                subprocess.run(f"masurca -i {fq1},{fq2} -t {threads} -o masurca_output -r {inf}", shell=True, stdout=f, stderr=f)
                CAdir = os.popen("ls -d CA.*").read().split("\n")[0]
                if os.path.isfile(f"{CAdir}/primary.genome.scf.fasta"):
                    subprocess.run(f"cp {CAdir}/primary.genome.scf.fasta {outputfa}", shell=True)
                elif os.path.isfile(f"{CAdir}/scaffolds.ref.fa"):
                    subprocess.run(f"cp {CAdir}/scaffolds.ref.fa {outputfa}", shell=True)
            else:
                subprocess.run(f"masurca -i {fq1} -t {threads} -o masurca_output -r {inf}", shell=True, stdout=f, stderr=f)
                CAdir = os.popen("ls -d CA.*").read().split("\n")[0]
                if os.path.isfile(f"{CAdir}/primary.genome.scf.fasta"):
                    subprocess.run(f"cp {CAdir}/primary.genome.scf.fasta {outputfa}", shell=True)
                elif os.path.isfile(f"{CAdir}/scaffolds.ref.fa"):
                    subprocess.run(f"cp {CAdir}/scaffolds.ref.fa {outputfa}", shell=True)


def reassm_fun(inf, fq1, fq2, threads, Pre, pts, pst, method, asmt, f, outputfa):
    runtime = get_runtime_context()
    inf, fq1, fq2 = _normalize_assembly_inputs(inf, fq1, fq2)
    runtime_ref = str(runtime.ref or "").strip()
    if _is_influenza_species(runtime.species) and (not runtime_ref or runtime_ref == "noref" or not Path(runtime_ref).is_file()):
        runtime_ref = resolve_influenza_reference(
            Pre,
            runtime.species,
            runtime_ref,
            single_fastq=inf,
            fq1=fq1,
            fq2=fq2,
            long_type=runtime.long_type,
            threads=threads,
            logf=f,
        )
        if runtime_ref:
            update_runtime_context(ref=runtime_ref)
    if _is_rsv_species(runtime.species) and (not runtime_ref or runtime_ref == "noref" or not Path(runtime_ref).is_file()):
        rsv_selection = resolve_rsv_reference(
            Pre,
            runtime.species,
            runtime_ref,
            single_fastq=inf,
            fq1=fq1,
            fq2=fq2,
            long_type=runtime.long_type,
            threads=threads,
            logf=f,
        )
        runtime_ref = str(rsv_selection.get("reference_path") or "").strip()
        runtime_gtf_selected = str(rsv_selection.get("gff_path") or "").strip() or "nogtf"
        runtime_species_selected = str(rsv_selection.get("species_label") or runtime.species).strip() or runtime.species
        if runtime_ref:
            update_runtime_context(ref=runtime_ref, gtf=runtime_gtf_selected, species=runtime_species_selected)
    if _is_denv_species(runtime.species) and (not runtime_ref or runtime_ref == "noref" or not Path(runtime_ref).is_file()):
        denv_selection = resolve_denv_reference(
            Pre,
            runtime.species,
            runtime_ref,
            single_fastq=inf,
            fq1=fq1,
            fq2=fq2,
            long_type=runtime.long_type,
            threads=threads,
            logf=f,
        )
        runtime_ref = str(denv_selection.get("reference_path") or "").strip()
        runtime_gtf_selected = str(denv_selection.get("gff_path") or "").strip() or "nogtf"
        runtime_species_selected = str(denv_selection.get("species_label") or runtime.species).strip() or runtime.species
        if runtime_ref:
            update_runtime_context(ref=runtime_ref, gtf=runtime_gtf_selected, species=runtime_species_selected)
    if _is_zika_species(runtime.species) and (not runtime_ref or runtime_ref == "noref" or not Path(runtime_ref).is_file()):
        zika_ref = PROJECT_ROOT / "database" / "nextclade_db" / "zikav" / "reference.fasta"
        zika_gtf = PROJECT_ROOT / "database" / "nextclade_db" / "zikav" / "genome_annotation.gff3"
        runtime_ref = str(zika_ref.resolve()) if zika_ref.is_file() else runtime_ref
        runtime_gtf_selected = str(zika_gtf.resolve()) if zika_gtf.is_file() else str(runtime.gtf or "nogtf").strip() or "nogtf"
        if runtime_ref and runtime_ref != "noref" and Path(runtime_ref).is_file():
            update_runtime_context(ref=runtime_ref, gtf=runtime_gtf_selected, species=runtime.species)
    if _is_chikv_species(runtime.species) and (not runtime_ref or runtime_ref == "noref" or not Path(runtime_ref).is_file()):
        chikv_ref = PROJECT_ROOT / "database" / "nextclade_db" / "chikv" / "reference.fasta"
        chikv_gtf = PROJECT_ROOT / "database" / "nextclade_db" / "chikv" / "genome_annotation.gff3"
        runtime_ref = str(chikv_ref.resolve()) if chikv_ref.is_file() else runtime_ref
        runtime_gtf_selected = str(chikv_gtf.resolve()) if chikv_gtf.is_file() else str(runtime.gtf or "nogtf").strip() or "nogtf"
        if runtime_ref and runtime_ref != "noref" and Path(runtime_ref).is_file():
            update_runtime_context(ref=runtime_ref, gtf=runtime_gtf_selected, species=runtime.species)
    if _is_hpiv_species(runtime.species) and (not runtime_ref or runtime_ref == "noref" or not Path(runtime_ref).is_file()):
        hpiv_selection = resolve_hpiv_reference(
            Pre,
            runtime.species,
            runtime_ref,
            single_fastq=inf,
            fq1=fq1,
            fq2=fq2,
            long_type=runtime.long_type,
            threads=threads,
            logf=f,
        )
        runtime_ref = str(hpiv_selection.get("reference_path") or "").strip()
        runtime_gtf_selected = str(hpiv_selection.get("gff_path") or "").strip() or "nogtf"
        runtime_species_selected = str(hpiv_selection.get("species_label") or runtime.species).strip() or runtime.species
        if runtime_ref:
            update_runtime_context(ref=runtime_ref, gtf=runtime_gtf_selected, species=runtime_species_selected)
    if _is_hadv_species(runtime.species) and (not runtime_ref or runtime_ref == "noref" or not Path(runtime_ref).is_file()):
        hadv_selection = resolve_hadv_reference(
            Pre,
            runtime.species,
            runtime_ref,
            single_fastq=inf,
            fq1=fq1,
            fq2=fq2,
            long_type=runtime.long_type,
            threads=threads,
            logf=f,
        )
        runtime_ref = str(hadv_selection.get("reference_path") or "").strip()
        runtime_gtf_selected = str(hadv_selection.get("gff_path") or "").strip() or "nogtf"
        runtime_species_selected = str(hadv_selection.get("species_label") or runtime.species).strip() or runtime.species
        if runtime_ref:
            update_runtime_context(ref=runtime_ref, gtf=runtime_gtf_selected, species=runtime_species_selected)
    if _is_norovirus_species(runtime.species) and (not runtime_ref or runtime_ref == "noref" or not Path(runtime_ref).is_file()):
        norovirus_selection = resolve_norovirus_reference(
            Pre,
            runtime.species,
            runtime_ref,
            single_fastq=inf,
            fq1=fq1,
            fq2=fq2,
            long_type=runtime.long_type,
            threads=threads,
            logf=f,
        )
        runtime_ref = str(norovirus_selection.get("reference_path") or "").strip()
        runtime_gtf_selected = str(norovirus_selection.get("gff_path") or "").strip() or "nogtf"
        runtime_species_selected = str(norovirus_selection.get("species_label") or runtime.species).strip() or runtime.species
        if runtime_ref:
            update_runtime_context(ref=runtime_ref, gtf=runtime_gtf_selected, species=runtime_species_selected)
    if _is_hepatovirus_species(runtime.species) and (not runtime_ref or runtime_ref == "noref" or not Path(runtime_ref).is_file()):
        hepatovirus_selection = resolve_hepatovirus_reference(
            Pre,
            runtime.species,
            runtime_ref,
            single_fastq=inf,
            fq1=fq1,
            fq2=fq2,
            long_type=runtime.long_type,
            threads=threads,
            logf=f,
        )
        runtime_ref = str(hepatovirus_selection.get("reference_path") or "").strip()
        runtime_gtf_selected = str(hepatovirus_selection.get("gff_path") or "").strip() or "nogtf"
        runtime_species_selected = str(hepatovirus_selection.get("species_label") or runtime.species).strip() or runtime.species
        if runtime_ref:
            update_runtime_context(ref=runtime_ref, gtf=runtime_gtf_selected, species=runtime_species_selected)
    if _is_hiv_species(runtime.species) and (not runtime_ref or runtime_ref == "noref" or not Path(runtime_ref).is_file()):
        hiv_selection = resolve_hiv_reference(
            Pre,
            runtime.species,
            runtime_ref,
            single_fastq=inf,
            fq1=fq1,
            fq2=fq2,
            long_type=runtime.long_type,
            threads=threads,
            logf=f,
        )
        runtime_ref = str(hiv_selection.get("reference_path") or "").strip()
        runtime_gtf_selected = str(hiv_selection.get("gff_path") or "").strip() or "nogtf"
        runtime_species_selected = str(hiv_selection.get("species_label") or runtime.species).strip() or runtime.species
        if runtime_ref:
            update_runtime_context(ref=runtime_ref, gtf=runtime_gtf_selected, species=runtime_species_selected)
    if _is_rhinovirus_species(runtime.species) and (not runtime_ref or runtime_ref == "noref" or not Path(runtime_ref).is_file()):
        rhinovirus_selection = resolve_rhinovirus_reference(
            Pre,
            runtime.species,
            runtime_ref,
            single_fastq=inf,
            fq1=fq1,
            fq2=fq2,
            long_type=runtime.long_type,
            threads=threads,
            logf=f,
        )
        runtime_ref = str(rhinovirus_selection.get("reference_path") or "").strip()
        runtime_gtf_selected = str(rhinovirus_selection.get("gff_path") or "").strip() or "nogtf"
        runtime_species_selected = str(rhinovirus_selection.get("species_label") or runtime.species).strip() or runtime.species
        if runtime_ref:
            update_runtime_context(ref=runtime_ref, gtf=runtime_gtf_selected, species=runtime_species_selected)
    if _is_enterovirus_species(runtime.species) and (not runtime_ref or runtime_ref == "noref" or not Path(runtime_ref).is_file()):
        enterovirus_selection = resolve_enterovirus_reference(
            Pre,
            runtime.species,
            runtime_ref,
            single_fastq=inf,
            fq1=fq1,
            fq2=fq2,
            long_type=runtime.long_type,
            threads=threads,
            logf=f,
        )
        runtime_ref = str(enterovirus_selection.get("reference_path") or "").strip()
        runtime_gtf_selected = str(enterovirus_selection.get("gff_path") or "").strip() or "nogtf"
        runtime_species_selected = str(enterovirus_selection.get("species_label") or runtime.species).strip() or runtime.species
        if runtime_ref:
            update_runtime_context(ref=runtime_ref, gtf=runtime_gtf_selected, species=runtime_species_selected)
    if _is_bandavirus_species(runtime.species) and (not runtime_ref or runtime_ref == "noref" or not Path(runtime_ref).is_file()):
        bandavirus_selection = resolve_bandavirus_reference(
            Pre,
            runtime.species,
            runtime_ref,
            single_fastq=inf,
            fq1=fq1,
            fq2=fq2,
            long_type=runtime.long_type,
            threads=threads,
            logf=f,
        )
        runtime_ref = str(bandavirus_selection.get("reference_path") or "").strip()
        runtime_gtf_selected = str(bandavirus_selection.get("gff_path") or "").strip() or "nogtf"
        runtime_species_selected = str(bandavirus_selection.get("species_label") or runtime.species).strip() or runtime.species
        if runtime_ref:
            update_runtime_context(ref=runtime_ref, gtf=runtime_gtf_selected, species=runtime_species_selected)
    if _is_orthohantavirus_species(runtime.species) and (not runtime_ref or runtime_ref == "noref" or not Path(runtime_ref).is_file()):
        orthohantavirus_selection = resolve_orthohantavirus_reference(
            Pre,
            runtime.species,
            runtime_ref,
            single_fastq=inf,
            fq1=fq1,
            fq2=fq2,
            long_type=runtime.long_type,
            threads=threads,
            logf=f,
        )
        runtime_ref = str(orthohantavirus_selection.get("reference_path") or "").strip()
        runtime_gtf_selected = str(orthohantavirus_selection.get("gff_path") or "").strip() or "nogtf"
        runtime_species_selected = str(orthohantavirus_selection.get("species_label") or runtime.species).strip() or runtime.species
        if runtime_ref:
            update_runtime_context(ref=runtime_ref, gtf=runtime_gtf_selected, species=runtime_species_selected)
    if _is_orthoebolavirus_species(runtime.species) and (not runtime_ref or runtime_ref == "noref" or not Path(runtime_ref).is_file()):
        orthoebolavirus_selection = resolve_orthoebolavirus_reference(
            Pre,
            runtime.species,
            runtime_ref,
            single_fastq=inf,
            fq1=fq1,
            fq2=fq2,
            long_type=runtime.long_type,
            threads=threads,
            logf=f,
        )
        runtime_ref = str(orthoebolavirus_selection.get("reference_path") or "").strip()
        runtime_gtf_selected = str(orthoebolavirus_selection.get("gff_path") or "").strip() or "nogtf"
        runtime_species_selected = str(orthoebolavirus_selection.get("species_label") or runtime.species).strip() or runtime.species
        if runtime_ref:
            update_runtime_context(ref=runtime_ref, gtf=runtime_gtf_selected, species=runtime_species_selected)
    if _is_astroviridae_species(runtime.species) and (not runtime_ref or runtime_ref == "noref" or not Path(runtime_ref).is_file()):
        astroviridae_selection = resolve_astroviridae_reference(
            Pre,
            runtime.species,
            runtime_ref,
            single_fastq=inf,
            fq1=fq1,
            fq2=fq2,
            long_type=runtime.long_type,
            threads=threads,
            logf=f,
        )
        runtime_ref = str(astroviridae_selection.get("reference_path") or "").strip()
        runtime_gtf_selected = str(astroviridae_selection.get("gff_path") or "").strip() or "nogtf"
        runtime_species_selected = str(astroviridae_selection.get("species_label") or runtime.species).strip() or runtime.species
        if runtime_ref:
            update_runtime_context(ref=runtime_ref, gtf=runtime_gtf_selected, species=runtime_species_selected)
    if _is_rotavirus_species(runtime.species) and (not runtime_ref or runtime_ref == "noref" or not Path(runtime_ref).is_file()):
        rotavirus_selection = resolve_rotavirus_reference(
            Pre,
            runtime.species,
            runtime_ref,
            single_fastq=inf,
            fq1=fq1,
            fq2=fq2,
            long_type=runtime.long_type,
            threads=threads,
            logf=f,
        )
        runtime_ref = str(rotavirus_selection.get("reference_path") or "").strip()
        runtime_gtf_selected = str(rotavirus_selection.get("gff_path") or "").strip() or "nogtf"
        runtime_species_selected = str(rotavirus_selection.get("species_label") or runtime.species).strip() or runtime.species
        if runtime_ref:
            update_runtime_context(ref=runtime_ref, gtf=runtime_gtf_selected, species=runtime_species_selected)
    if _is_seasonal_hcov_species(runtime.species) and (not runtime_ref or runtime_ref == "noref" or not Path(runtime_ref).is_file()):
        seasonal_hcov_selection = resolve_seasonal_hcov_reference(
            Pre,
            runtime.species,
            runtime_ref,
            single_fastq=inf,
            fq1=fq1,
            fq2=fq2,
            long_type=runtime.long_type,
            threads=threads,
            logf=f,
        )
        runtime_ref = str(seasonal_hcov_selection.get("reference_path") or "").strip()
        runtime_gtf_selected = str(seasonal_hcov_selection.get("gff_path") or "").strip() or "nogtf"
        runtime_species_selected = str(seasonal_hcov_selection.get("species_label") or runtime.species).strip() or runtime.species
        if runtime_ref:
            update_runtime_context(ref=runtime_ref, gtf=runtime_gtf_selected, species=runtime_species_selected)
    if not runtime_ref or runtime_ref == "noref" or not Path(runtime_ref).is_file():
        raise FileNotFoundError(f"未找到可用于有参组装的参考基因组：{runtime_ref or 'noref'}")
    runtime_gtf = runtime.gtf
    if _is_rsv_species(runtime.species) and (not str(runtime_gtf or "").strip() or str(runtime_gtf).strip() == "nogtf"):
        runtime_ref_lower = str(runtime_ref or "").strip().lower()
        selected_type = "A" if "rsv_a" in runtime_ref_lower else ("B" if "rsv_b" in runtime_ref_lower else "")
        if selected_type:
            rsv_selection = resolve_rsv_reference(
                Pre,
                f"Respiratory syncytial virus {selected_type}",
                runtime_ref,
                single_fastq=inf,
                fq1=fq1,
                fq2=fq2,
                long_type=runtime.long_type,
                threads=threads,
                logf=f,
            )
            runtime_gtf = str(rsv_selection.get("gff_path") or "").strip() or runtime_gtf
            update_runtime_context(gtf=runtime_gtf)
    if _is_denv_species(runtime.species) and (not str(runtime_gtf or "").strip() or str(runtime_gtf).strip() == "nogtf"):
        runtime_ref_lower = str(runtime_ref or "").strip().lower()
        selected_type = ""
        for denv_type in ["1", "2", "3", "4"]:
            if f"denv{denv_type}" in runtime_ref_lower:
                selected_type = denv_type
                break
        if selected_type:
            denv_selection = resolve_denv_reference(
                Pre,
                f"Dengue virus {selected_type}",
                runtime_ref,
                single_fastq=inf,
                fq1=fq1,
                fq2=fq2,
                long_type=runtime.long_type,
                threads=threads,
                logf=f,
            )
            runtime_gtf = str(denv_selection.get("gff_path") or "").strip() or runtime_gtf
            update_runtime_context(gtf=runtime_gtf)
    if _is_zika_species(runtime.species) and (not str(runtime_gtf or "").strip() or str(runtime_gtf).strip() == "nogtf"):
        zika_gtf = PROJECT_ROOT / "database" / "nextclade_db" / "zikav" / "genome_annotation.gff3"
        if zika_gtf.is_file():
            runtime_gtf = str(zika_gtf.resolve())
            update_runtime_context(gtf=runtime_gtf)
    if _is_chikv_species(runtime.species) and (not str(runtime_gtf or "").strip() or str(runtime_gtf).strip() == "nogtf"):
        chikv_gtf = PROJECT_ROOT / "database" / "nextclade_db" / "chikv" / "genome_annotation.gff3"
        if chikv_gtf.is_file():
            runtime_gtf = str(chikv_gtf.resolve())
            update_runtime_context(gtf=runtime_gtf)
    if _is_hpiv_species(runtime.species) and (not str(runtime_gtf or "").strip() or str(runtime_gtf).strip() == "nogtf"):
        runtime_ref_lower = str(runtime_ref or "").strip().lower()
        selected_type = ""
        for hpiv_type, token in [("1", "hpiv1"), ("2", "hpiv2"), ("3", "hpiv3"), ("4A", "hpiv4a"), ("4B", "hpiv4b")]:
            if token in runtime_ref_lower:
                selected_type = hpiv_type
                break
        if selected_type:
            hpiv_selection = resolve_hpiv_reference(
                Pre,
                f"Human parainfluenza virus {selected_type}",
                runtime_ref,
                single_fastq=inf,
                fq1=fq1,
                fq2=fq2,
                long_type=runtime.long_type,
                threads=threads,
                logf=f,
            )
            runtime_gtf = str(hpiv_selection.get("gff_path") or "").strip() or runtime_gtf
            update_runtime_context(gtf=runtime_gtf)
    if _is_hadv_species(runtime.species) and (not str(runtime_gtf or "").strip() or str(runtime_gtf).strip() == "nogtf"):
        selection_path = Path(f"{Pre}_hadv_reference_selection") / "selection.tsv"
        if selection_path.is_file():
            with selection_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            if rows:
                selected_gff = str(rows[0].get("gff_path") or "").strip()
                if selected_gff and selected_gff != "nogtf" and Path(selected_gff).is_file():
                    runtime_gtf = selected_gff
                    update_runtime_context(gtf=runtime_gtf)
    variant_caller = _choose_reference_variant_caller(method, asmt)
    if not os.path.isdir("genomes"):
        os.makedirs("genomes")
    subprocess.run(f"seqkit seq {runtime_ref} > genomes/ref.fa", shell=True)
    subprocess.run("samtools faidx genomes/ref.fa ", shell=True)
    if runtime_gtf == "nogtf" and _is_monkeypox_species(runtime.species):
        generated_gtf = prepare_monkeypox_reference_annotation(
            Pre,
            Path("genomes/ref.fa"),
            Path.cwd(),
            logf=f,
        )
        if generated_gtf is not None:
            runtime_gtf = str(generated_gtf.resolve())
            update_runtime_context(gtf=runtime_gtf)
    if runtime_gtf == "nogtf" and _is_norovirus_species(runtime.species):
        generated_gtf = prepare_norovirus_reference_annotation(
            Pre,
            Path("genomes/ref.fa"),
            Path.cwd(),
            logf=f,
        )
        if generated_gtf is not None:
            runtime_gtf = str(generated_gtf.resolve())
            update_runtime_context(gtf=runtime_gtf)
            selection_path = Path(f"{Pre}_norovirus_reference_selection") / "selection.tsv"
            if selection_path.is_file():
                with selection_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
                    rows = list(csv.DictReader(handle, delimiter="\t"))
                    fieldnames = list(rows[0].keys()) if rows else []
                if rows and "gff_path" in fieldnames:
                    rows[0]["gff_path"] = runtime_gtf
                    with selection_path.open("w", encoding="utf-8", newline="") as handle:
                        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
                        writer.writeheader()
                        writer.writerows(rows)
    if _is_rhinovirus_species(runtime.species):
        selection_path = Path(f"{Pre}_rhinovirus_reference_selection") / "selection.tsv"
        species_group = ""
        rows = []
        fieldnames: list[str] = []
        if selection_path.is_file():
            with selection_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
                fieldnames = list(rows[0].keys()) if rows else []
            if rows:
                species_group = str(rows[0].get("species_group") or "").strip().upper()
        if species_group in {"A", "B", "C"}:
            generated_gtf = prepare_rhinovirus_reference_annotation(
                Pre,
                Path("genomes/ref.fa"),
                Path.cwd(),
                species_group,
                logf=f,
            )
            if generated_gtf is not None:
                runtime_gtf = str(generated_gtf.resolve())
                update_runtime_context(gtf=runtime_gtf)
                if selection_path.is_file() and rows and "gff_path" in fieldnames:
                    rows[0]["gff_path"] = runtime_gtf
                    with selection_path.open("w", encoding="utf-8", newline="") as handle:
                        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
                        writer.writeheader()
                        writer.writerows(rows)
    if _is_enterovirus_species(runtime.species):
        selection_path = Path(f"{Pre}_enterovirus_reference_selection") / "selection.tsv"
        big_group = ""
        rows = []
        fieldnames: list[str] = []
        if selection_path.is_file():
            with selection_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
                fieldnames = list(rows[0].keys()) if rows else []
            if rows:
                big_group = str(rows[0].get("big_group") or "").strip().upper()
        if big_group in {"A", "B", "C", "D"}:
            generated_gtf = prepare_enterovirus_reference_annotation(
                Pre,
                Path("genomes/ref.fa"),
                Path.cwd(),
                big_group,
                logf=f,
            )
            if generated_gtf is not None:
                runtime_gtf = str(generated_gtf.resolve())
                update_runtime_context(gtf=runtime_gtf)
                if selection_path.is_file() and rows and "gff_path" in fieldnames:
                    rows[0]["gff_path"] = runtime_gtf
                    with selection_path.open("w", encoding="utf-8", newline="") as handle:
                        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
                        writer.writeheader()
                        writer.writerows(rows)
    if _is_astroviridae_species(runtime.species):
        selection_path = Path(f"{Pre}_astroviridae_reference_selection") / "selection.tsv"
        astro_genus = ""
        rows = []
        fieldnames: list[str] = []
        if selection_path.is_file():
            with selection_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
                fieldnames = list(rows[0].keys()) if rows else []
            if rows:
                astro_genus = str(rows[0].get("genus") or "").strip()
        if astro_genus in {"Mamastrovirus", "Avastrovirus"}:
            generated_gtf = prepare_astroviridae_reference_annotation(
                Pre,
                Path("genomes/ref.fa"),
                Path.cwd(),
                astro_genus,
                logf=f,
            )
            if generated_gtf is not None:
                runtime_gtf = str(generated_gtf.resolve())
                update_runtime_context(gtf=runtime_gtf)
                if selection_path.is_file() and rows and "gff_path" in fieldnames:
                    rows[0]["gff_path"] = runtime_gtf
                    with selection_path.open("w", encoding="utf-8", newline="") as handle:
                        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
                        writer.writeheader()
                        writer.writerows(rows)
    if _is_seasonal_hcov_species(runtime.species) and (not str(runtime_gtf or "").strip() or str(runtime_gtf).strip() == "nogtf"):
        selection_path = Path(f"{Pre}_seasonal_hcov_reference_selection") / "selection.tsv"
        hcov_type = ""
        rows = []
        fieldnames: list[str] = []
        if selection_path.is_file():
            with selection_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
                fieldnames = list(rows[0].keys()) if rows else []
            if rows:
                hcov_type = str(rows[0].get("hcov_type") or "").strip()
        if hcov_type:
            generated_gtf = prepare_seasonal_hcov_reference_annotation(
                Pre,
                Path("genomes/ref.fa"),
                Path.cwd(),
                hcov_type,
                logf=f,
            )
            if generated_gtf is not None:
                runtime_gtf = str(generated_gtf.resolve())
                update_runtime_context(gtf=runtime_gtf)
                if selection_path.is_file() and rows and "gff_path" in fieldnames:
                    rows[0]["gff_path"] = runtime_gtf
                    with selection_path.open("w", encoding="utf-8", newline="") as handle:
                        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
                        writer.writeheader()
                        writer.writerows(rows)
    if runtime_gtf == "nogtf":
        ifAvcf = 0
        print("snp位点不进行额外注释")
        sys.stdout.flush()
    else:
        _, snpeff_jar = _resolve_project_snpeff_paths()
        if not os.path.isdir("ref"):
            os.makedirs("ref")
        prepared_gff = prepare_snpeff_reference_gff(Path(runtime_gtf), Path("ref/genes.gff"))
        if (not prepared_gff) and os.path.isfile(runtime_gtf):
            shutil.copy2(runtime_gtf, "ref/genes.gff")
        open("snpEff.config", "w").write("ref.genome : ref")
        subprocess.run(
            f"java -jar {shlex.quote(str(snpeff_jar))} build -noCheckCds -noCheckProtein -gff3 ref -c snpEff.config -dataDir ./",
            shell=True,
            stdout=f,
            stderr=f,
        )
        if not os.path.isfile("ref/snpEffectPredictor.bin"):
            print("gff与fa不匹配，snp位点不进行额外注释")
            ifAvcf = 0
            sys.stdout.flush()
        else:
            print("注释文件正常，snp位点根据注释文件进行注释")
            ifAvcf = 1
            sys.stdout.flush()

    influenza_mode = _is_influenza_species(runtime.species) and asmt in {"shortref", "longref"}
    influenza_reads = inf
    if influenza_mode:
        flu_root = Path("wf_flu")
        flu_reads_dir = flu_root / "reads"
        flu_reads_dir.mkdir(parents=True, exist_ok=True)
        if asmt == "longref":
            influenza_reads = _filter_influenza_reads(inf, flu_reads_dir / "influenza.filtered.fastq", 200, logf=f)
            inf = influenza_reads
        else:
            if fq1:
                fq1 = _filter_influenza_reads(fq1, flu_reads_dir / "influenza.R1.filtered.fastq", 50, logf=f)
            if fq2:
                fq2 = _filter_influenza_reads(fq2, flu_reads_dir / "influenza.R2.filtered.fastq", 50, logf=f)

    if "long" in asmt:
        if not inf:
            raise FileNotFoundError(f"{Pre} 缺少可用长读输入文件，无法执行 {asmt} 有参组装")
        subprocess.run(f"minimap2 -ax map-ont genomes/ref.fa {inf} -t {threads} |samtools sort -o ref.mapping.bam", shell=True)
    elif "short" in asmt:
        with open("mapping.log", "w") as mapplog:
            subprocess.run("bwa index genomes/ref.fa", shell=True, stdout=mapplog, stderr=mapplog)
            print(fq1,fq2)
            if fq2:
                subprocess.run(f"bwa mem -t {threads} genomes/ref.fa {fq1} {fq2}  |samtools sort -o ref.mapping.bam", shell=True, stdout=mapplog, stderr=mapplog)
            else:
                subprocess.run(f"bwa mem  -t {threads} genomes/ref.fa {fq1}  |samtools sort -o ref.mapping.bam", shell=True, stdout=mapplog, stderr=mapplog)
    if runtime.analysis_target != 'virus':
        subprocess.run("samtools index ref.mapping.bam", shell=True)
        subprocess.run(f"mosdepth -b 1000 ref_map ref.mapping.bam -t {threads}", shell=True)
        subprocess.run("gunzip ref_map.regions.bed.gz", shell=True)
        subprocess.run(f"mosdepth -b 1 ref_map1 ref.mapping.bam -t {threads}", shell=True)
        subprocess.run("gunzip ref_map1.regions.bed.gz", shell=True)
        subprocess.run("cp ref_map1.regions.bed ref.regions.bed", shell=True)  
    else:
        subprocess.run("samtools index ref.mapping.bam", shell=True)
        subprocess.run("samtools depth -aa ref.mapping.bam > ref_map.pre.regions.bed", shell=True)
        tmppd = pd.read_table('ref_map.pre.regions.bed',header=None)
        tmppd['oldpos'] = tmppd[1]-1
        tmppd[[0,'oldpos',1,2]].to_csv('ref_map.regions.bed',index=False,sep='\t',header=False)
        subprocess.run("cp ref_map.regions.bed ref.regions.bed", shell=True)    
    def _outfun(x):
        tmpdict = {}
        ofname = x["GeneName"].tolist()[0]
        x["start"] = x.reset_index().index + 1
        x["end"] = x.reset_index().index + 2
        x[["Chrom", "start", "end", "Depth"]].to_csv(f"geneDepth/{ofname}.tsv", sep="\t", header=False, index=False)
        tmpdict["片段名称"] = x["Chrom"].tolist()[0]
        tmpdict["起始位置"] = x["start"].min()
        tmpdict["终止位置"] = x["start"].max()
        tmpdict["覆盖度(>0)%"] = round(x[x["Depth"] > 0].shape[0] / x.shape[0], 4) * 100
        tmpdict["覆盖度(>10)%"] = round(x[x["Depth"] > 10].shape[0] / x.shape[0], 4) * 100
        tmpdict["覆盖度(>100)%"] = round(x[x["Depth"] > 100].shape[0] / x.shape[0], 4) * 100
        tmpdict["平均深度"] = round(x["Depth"].mean(), 2)
        tmpdict["最低深度"] = x["Depth"].min()
        tmpdict["最高深度"] = x["Depth"].max()
        return pd.DataFrame(tmpdict, index=[0]).round(2)

    if os.path.isfile("ref/genes.gff"):
        open("geneNamelist.txt", "w").write("")
        if not os.path.isdir("geneDepth"):
            os.makedirs("geneDepth")
        if runtime.analysis_target != 'virus':
            with open("ref/genes.gff") as f1:
                for line in f1:
                    if not line.startswith("#"):
                        line = line.strip().split("\t")
                        if line[2] == "gene":
                            if "gene=" in line[8]:
                                gName = line[8].split("gene=")[1].split(";")[0].split("/")[0]
                                open(f"geneDepth/{gName}.bed", "w").write(f"{line[0]}\t{line[3]}\t{line[4]}\t{gName}\n")
                                open("geneNamelist.txt", "a").write(f"{gName}\n")
                            else:
                                if "ID=" in line[8]:
                                    gName = line[8].split("ID=")[1].split(";")[0]
                                    open(f"geneDepth/{gName}.bed", "w").write(f"{line[0]}\t{line[3]}\t{line[4]}\t{gName}\n")
                                    open("geneNamelist.txt", "a").write(f"{gName}\n")
            subprocess.run("cat geneDepth/*.bed > All_gene.bed", shell=True)
            subprocess.run("bedtools intersect -a All_gene.bed -b ref_map1.regions.bed -wb > ref_map1.Anno.regions.bed", shell=True)
            tmpd = pd.read_table("ref_map1.Anno.regions.bed", header=None, names=["Chrom", "start", "end", "GeneName", "c1", "s1", "e1", "Depth"])
            if tmpd.shape[0] != 0:
                newtd = tmpd.groupby("GeneName").apply(_outfun).reset_index(level=0)
                newtd.rename(columns={"GeneName": "基因名称"}, inplace=True)
                newtd.to_csv("gene_summary.tsv", sep="\t", index=False)
            else:
                print("gff与基因组不匹配，无法展示各个基因区段覆盖度")
                sys.stdout.flush()
        
        #subprocess.run(f"python /data1/shanghai_pip/meta_genome/soft/IGV_js/IGV_new.py -r genomes/ref.fa -m ref.mapping.bam -o ./ -s {Pre}", shell=True)
        #subprocess.run("cp /data1/shanghai_pip/meta_genome/soft/IGV_js/igv.min.js ./", shell=True)
    #subprocess.run("mosdepth -b 1 -t 10 ref ref.mapping.bam", shell=True)
    #subprocess.run("gunzip ref.regions.bed.gz", shell=True)
    refbeddb = pd.read_table("ref.regions.bed", header=None)
    refbeddb[refbeddb[3] == 0].to_csv("mask.bed", index=False, header=False, sep="\t")
    if variant_caller == "freebayes":
        if runtime.analysis_target != 'virus':
            subprocess.run("fasta_generate_regions.py genomes/ref.fa.fai 200000 > ref.txt", shell=True)
        else:
            subprocess.run("fasta_generate_regions.py genomes/ref.fa.fai 1000 > ref.txt", shell=True)
        with open("snps.raw.vcf", "w", encoding="utf-8") as raw_vcf_handle:
            subprocess.run(
                [
                    "freebayes-parallel",
                    "ref.txt",
                    str(threads),
                    "-p",
                    "2",
                    "-P",
                    "0",
                    "-C",
                    "2",
                    "-F",
                    "0.05",
                    "--min-coverage",
                    "10",
                    "--min-repeat-entropy",
                    "1.0",
                    "-q",
                    "30",
                    "-m",
                    "30",
                    "--strict-vcf",
                    "-f",
                    "genomes/ref.fa",
                    "ref.mapping.bam",
                ],
                stdout=raw_vcf_handle,
                stderr=f,
                check=False,
            )
        if _is_sars_cov_2_species(runtime.species):
            print("检测到新冠病毒样本，开始对 snps.raw.vcf 执行 snpEff 注释并生成突变位点表")
            sys.stdout.flush()
            _run_ncov_raw_vcf_annotation()
        if _is_monkeypox_species(runtime.species):
            _filter_monkeypox_raw_vcf_to_filt1(Path("snps.raw.vcf"), Path("snps.filt1.vcf"))
        else:
            subprocess.run("bcftools view --include 'QUAL>=20 && FMT/DP>=10 && (FMT/AO+FMT/RO)>0 && (FMT/AO)/(FMT/AO+FMT/RO)>=0.9' snps.raw.vcf  | bcftools annotate --remove '^INFO/TYPE,^INFO/DP,^INFO/RO,^INFO/AO,^INFO/AB,^FORMAT/GT,^FORMAT/DP,^FORMAT/RO,^FORMAT/AO,^FORMAT/QR,^FORMAT/QA,^FORMAT/GL' > snps.filt1.vcf", shell=True)
        if ifAvcf:
            _, snpeff_jar = _resolve_project_snpeff_paths()
            subprocess.run(f"java -jar {shlex.quote(str(snpeff_jar))} ann -c snpEff.config -Datadir . ref snps.filt1.vcf > snps.anno.vcf", shell=True)
            skinums = int(os.popen("grep '##' snps.anno.vcf |wc -l").read())
            annodb = pd.read_table("snps.anno.vcf", skiprows=skinums)
            annodb["参考碱基"] = annodb["REF"] + "(" + annodb["unknown"].str.split("|").str[-1].str.strip().str.split(":").str[-4] + ")"
            annodb["突变碱基"] = annodb["ALT"] + "(" + annodb["unknown"].str.split("|").str[-1].str.strip().str.split(":").str[-3] + ")"
            annodb["突变类型"] = annodb["INFO"].str.split("|").str[1]
            annodb["突变影响"] = annodb["INFO"].str.split("|").str[2]
            annodb["影响基因"] = annodb["INFO"].str.split("|").str[3]
            annodb["碱基变化"] = annodb["INFO"].str.split("|").str[9]
            annodb["氨基酸变化"] = annodb["INFO"].str.split("|").str[10]
            annodb["突变位置"] = annodb["POS"]
            annodb["片段名称"] = annodb["#CHROM"]
            annodb = annodb[["片段名称", "突变位置", "参考碱基", "突变碱基", "影响基因", "突变类型", "突变影响", "碱基变化", "氨基酸变化"]]
            annodb["氨基酸变化"] = annodb.apply(lambda x: "-" if x["氨基酸变化"] == "" else x["氨基酸变化"], axis=1)
            annodb.fillna("-", inplace=True)
            annodb.to_csv(f"{Pre}.anno.tsv", sep="\t", index=False)
        else:
            subprocess.run("ln -s snps.filt1.vcf snps.anno.vcf", shell=True)
            skinums = int(os.popen("grep '##' snps.anno.vcf |wc -l").read())
            annodb = pd.read_table("snps.anno.vcf", skiprows=skinums)
            annodb["参考碱基"] = annodb["REF"] + "(" + annodb["unknown"].str.split("|").str[-1].str.strip().str.split(":").str[-5] + ")"
            annodb["突变碱基"] = annodb["ALT"] + "(" + annodb["unknown"].str.split("|").str[-1].str.strip().str.split(":").str[-3] + ")"
            annodb["突变位置"] = annodb["POS"]
            annodb["片段名称"] = annodb["#CHROM"]
            annodb = annodb[["片段名称", "突变位置", "参考碱基", "突变碱基"]]
            annodb.to_csv(f"{Pre}.anno.tsv", sep="\t", index=False)
    elif variant_caller == "clair3":
        with open("clair3.log", "w") as clg:
            clair3_cmd = f"run_clair3.sh --bam_fn=ref.mapping.bam --ref_fn=genomes/ref.fa --threads={threads} --platform=ont --model_path=${{CONDA_PREFIX}}/bin/models/r941_prom_sup_g5014 --output=./ --include_all_ctgs --enable_long_indel --snp_min_af=0.05"
            subprocess.run(conda_run_command("longread_aux", f"bash -lc {shlex.quote(clair3_cmd)}"), shell=True, stdout=clg, stderr=clg)
        subprocess.run("samtools view -h -F 2308 ref.mapping.bam |samtools sort -o ref.filter.bam", shell=True)
        subprocess.run("samtools index ref.filter.bam", shell=True)
        subprocess.run(f"perbase base-depth ref.filter.bam -t {threads} > perbase.bed", shell=True)
        pd.read_table("perbase.bed")
        if ifAvcf:
            _, snpeff_jar = _resolve_project_snpeff_paths()
            subprocess.run(f"java -jar {shlex.quote(str(snpeff_jar))} ann -c snpEff.config -Datadir . ref merge_output.vcf.gz > snps.anno.vcf", shell=True)
            skinums = int(os.popen("grep '##' snps.anno.vcf |wc -l").read())
            annodb = pd.read_table("snps.anno.vcf", skiprows=skinums)
            annodb["参考碱基"] = annodb["REF"] + "(" + annodb["SAMPLE"].str.split("|").str[-1].str.strip().str.split(":").str[-2].str.split(",").str[0] + ")"
            annodb["突变碱基"] = annodb["ALT"] + "(" + annodb["SAMPLE"].str.split("|").str[-1].str.strip().str.split(":").str[-2].str.split(",").str[1] + ")"
            annodb["突变类型"] = annodb["INFO"].str.split("|").str[1]
            annodb["突变影响"] = annodb["INFO"].str.split("|").str[2]
            annodb["影响基因"] = annodb["INFO"].str.split("|").str[3]
            annodb["碱基变化"] = annodb["INFO"].str.split("|").str[9]
            annodb["氨基酸变化"] = annodb["INFO"].str.split("|").str[10]
            annodb["突变位置"] = annodb["POS"]
            annodb["片段名称"] = annodb["#CHROM"]
            annodb = annodb[["片段名称", "突变位置", "参考碱基", "突变碱基", "影响基因", "突变类型", "突变影响", "碱基变化", "氨基酸变化"]]
            annodb["氨基酸变化"] = annodb.apply(lambda x: "-" if x["氨基酸变化"] == "" else x["氨基酸变化"], axis=1)
            annodb.to_csv(f"{Pre}.anno.tsv", sep="\t", index=False)
        else:
            subprocess.run("gunzip merge_output.vcf.gz -c > snps.anno.vcf", shell=True)
            skinums = int(os.popen("grep '##' snps.anno.vcf |wc -l").read())
            annodb = pd.read_table("snps.anno.vcf", skiprows=skinums)
            annodb["参考碱基"] = annodb["REF"] + "(" + annodb["SAMPLE"].str.split("|").str[-1].str.strip().str.split(":").str[-2].str.split(",").str[0] + ")"
            annodb["突变碱基"] = annodb["ALT"] + "(" + annodb["SAMPLE"].str.split("|").str[-1].str.strip().str.split(":").str[-2].str.split(",").str[1] + ")"
            annodb["突变位置"] = annodb["POS"]
            annodb["片段名称"] = annodb["#CHROM"]
            annodb = annodb[["片段名称", "突变位置", "参考碱基", "突变碱基"]]
            annodb.to_csv(f"{Pre}.anno.tsv", sep="\t", index=False)
    subprocess.run("bcftools convert -Oz -o snps.anno.vcf.gz snps.anno.vcf", shell=True)
    subprocess.run("bcftools index -f snps.anno.vcf.gz", shell=True)
    if outputfa != "noforce":
        subprocess.run(f"bcftools consensus -f genomes/ref.fa -o {outputfa} snps.anno.vcf.gz -m mask.bed", shell=True)
        if _is_hpiv_species(runtime.species):
            build_hpiv_coverage_assets(
                Pre,
                runtime.species,
                single_fastq=inf,
                fq1=fq1,
                fq2=fq2,
                long_type=runtime.long_type,
                threads=threads,
                logf=f,
            )
        if _is_norovirus_species(runtime.species):
            run_norovirus_consensus_typing(
                Pre,
                Path(outputfa),
                logf=f,
            )
            if os.path.isfile("ref/genes.gff"):
                build_norovirus_gene_phylogeny_assets(
                    Pre,
                    Path(outputfa),
                    Path("ref/genes.gff"),
                    logf=f,
                )
        if _is_hepatovirus_species(runtime.species):
            run_hepatovirus_consensus_typing(
                Pre,
                Path(outputfa),
                logf=f,
            )
        if _is_rhinovirus_species(runtime.species):
            rhinovirus_consensus = run_rhinovirus_consensus_typing(
                Pre,
                Path(outputfa),
                logf=f,
            )
            species_group = str(rhinovirus_consensus.get("species_group") or "").strip().upper()
            if species_group in {"A", "B", "C"}:
                sample_gff = prepare_rhinovirus_sample_annotation(
                    Pre,
                    Path(outputfa),
                    Path.cwd(),
                    species_group,
                    logf=f,
                )
                if sample_gff is not None and Path(sample_gff).is_file():
                    build_rhinovirus_vp1_phylogeny_assets(
                        Pre,
                        Path(outputfa),
                        Path(sample_gff),
                        species_group=species_group,
                        logf=f,
                    )
        if _is_enterovirus_species(runtime.species):
            enterovirus_consensus = run_enterovirus_consensus_typing(
                Pre,
                Path(outputfa),
                logf=f,
            )
            big_group = str(enterovirus_consensus.get("big_group") or "").strip().upper()
            if big_group in {"A", "B", "C", "D"}:
                sample_gff = prepare_enterovirus_sample_annotation(
                    Pre,
                    Path(outputfa),
                    Path.cwd(),
                    big_group,
                    logf=f,
                )
                if sample_gff is not None and Path(sample_gff).is_file():
                    build_enterovirus_vp1_phylogeny_assets(
                        Pre,
                        Path(outputfa),
                        Path(sample_gff),
                        big_group=big_group,
                        logf=f,
                    )
        if _is_bandavirus_species(runtime.species):
            run_bandavirus_consensus_typing(
                Pre,
                Path(outputfa),
                logf=f,
            )
        if _is_orthohantavirus_species(runtime.species):
            run_orthohantavirus_consensus_typing(
                Pre,
                Path(outputfa),
                logf=f,
            )
        if _is_orthoebolavirus_species(runtime.species):
            run_orthoebolavirus_consensus_typing(
                Pre,
                Path(outputfa),
                logf=f,
            )
        if _is_astroviridae_species(runtime.species):
            astroviridae_consensus = run_astroviridae_consensus_typing(
                Pre,
                Path(outputfa),
                logf=f,
            )
            astro_genus = str(astroviridae_consensus.get("genus") or "").strip()
            if astro_genus in {"Mamastrovirus", "Avastrovirus"}:
                sample_gff = prepare_astroviridae_sample_annotation(
                    Pre,
                    Path(outputfa),
                    Path.cwd(),
                    astro_genus,
                    logf=f,
                )
                if sample_gff is not None and Path(sample_gff).is_file():
                    build_astroviridae_orf2_phylogeny_assets(
                        Pre,
                        Path(outputfa),
                        Path(sample_gff),
                        genus=astro_genus,
                        logf=f,
                    )
        if _is_rotavirus_species(runtime.species):
            run_rotavirus_consensus_typing(
                Pre,
                Path(outputfa),
                logf=f,
            )
        if _is_seasonal_hcov_species(runtime.species):
            selection_path = Path(f"{Pre}_seasonal_hcov_reference_selection") / "selection.tsv"
            hcov_type = ""
            if selection_path.is_file():
                with selection_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
                    row = next(csv.DictReader(handle, delimiter="\t"), None) or {}
                hcov_type = str(row.get("hcov_type") or "").strip()
            if hcov_type:
                sample_gff = prepare_seasonal_hcov_sample_annotation(
                    Pre,
                    Path(outputfa),
                    Path.cwd(),
                    hcov_type,
                    logf=f,
                )
                if sample_gff is not None and Path(sample_gff).is_file():
                    build_seasonal_hcov_spike_phylogeny_assets(
                        Pre,
                        Path(outputfa),
                        Path(sample_gff),
                        hcov_type,
                        logf=f,
                    )
        if influenza_mode:
            _export_influenza_reference_outputs(
                Pre,
                threads,
                outputfa,
                variant_caller,
                long_reads=influenza_reads if asmt == "longref" else None,
                polish_times=pts,
                polish_soft=pst,
                logf=f,
            )


def renamefa(inf, ofn):
    subprocess.run(f"seqkit sort -l -r {inf}|seqkit fx2tab > tmpfa.tab", shell=True)
    afile = pd.read_table("tmpfa.tab", header=None)
    afile["contignum"] = afile.index + 1
    runtime = get_runtime_context()
    if runtime.analysis_target != "virus":
        if afile.loc[afile[1].str.len() > 1000, :].shape[0] > 0:
            afile = afile.loc[afile[1].str.len() > 1000, :]
    afile[3] = afile.apply(lambda x: f"contig_{x.contignum}", axis=1)
    afile[[0, 3]].to_csv("transname.tsv", sep="\t", index=False, header=False)
    with open(ofn, "w", encoding="utf-8") as handle:
        for _, row in afile.iterrows():
            handle.write(f">contig_{row.contignum}\n")
            handle.write(f"{str(row[1]).strip()}\n")


def _nonempty_file(path: str | Path) -> bool:
    try:
        return Path(path).is_file() and Path(path).stat().st_size > 0
    except OSError:
        return False


def _extracted_reads_ready(read_prefix: str, expect_read2: bool) -> bool:
    if not _nonempty_file(f"{read_prefix}.1.fastq"):
        return False
    return not expect_read2 or _nonempty_file(f"{read_prefix}.2.fastq")


def _reuse_extracted_reads(read_prefix: str, expect_read2: bool, logf=None):
    if logf is not None:
        read2_msg = f" and {read_prefix}.2.fastq" if expect_read2 else ""
        logf.write(f"[skip] reuse existing non-empty {read_prefix}.1.fastq{read2_msg}; skip kraken2 read extraction\n")
        logf.flush()
    Path("kk2_ok").touch()
    return f"{read_prefix}.1.fastq", f"{read_prefix}.2.fastq" if expect_read2 else 0


def prepare_assembly_inputs(inf, fq1, fq2, threads, Pre, lelID, method):
    runtime = get_runtime_context()
    krdb = runtime.krdb
    is_virus = runtime.analysis_target == "virus"
    use_bracken = not is_virus
    default_taxid = "10239" if is_virus else "2"
    inf_path = str(inf).strip() if inf not in {None, 0} else ""
    fq1_path = str(fq1).strip() if fq1 not in {None, 0} else ""
    fq2_path = str(fq2).strip() if fq2 not in {None, 0} else ""
    has_inf = bool(inf_path) and Path(inf_path).is_file()
    has_fq1 = bool(fq1_path) and Path(fq1_path).is_file()
    has_fq2 = bool(fq2_path) and Path(fq2_path).is_file()
    expected_read_prefix = str(lelID).split(",", 1)[0] if lelID != "nolevel" else default_taxid
    extracted_reads_ready = _extracted_reads_ready(expected_read_prefix, has_fq2)
    with open("tmpkk2.log", "w") as kkf:
        if extracted_reads_ready:
            _reuse_extracted_reads(expected_read_prefix, has_fq2, kkf)
        else:
            if has_inf:
                print(inf_path)
                subprocess.run(f"kraken2 --db {krdb} --threads {threads} --output {Pre}.out.txt --report {Pre}.report.txt {inf_path}", shell=True, stdout=kkf, stderr=kkf)
                if use_bracken:
                    subprocess.run(f"bracken -d {krdb} -o {Pre}.bracken1.txt -w {Pre}.bracken2.txt -l S -t 10  -i {Pre}.report.txt", shell=True, stdout=kkf, stderr=kkf)

            if has_fq1 and has_fq2:
                subprocess.run(f"kraken2 --db {krdb} --threads {threads} --output {Pre}_2.out.txt --report {Pre}_2.report.txt {fq1_path} {fq2_path}", shell=True, stdout=kkf, stderr=kkf)
                if use_bracken:
                    subprocess.run(f"bracken -d {krdb} -o {Pre}_2.bracken1.txt -w {Pre}_2.bracken2.txt -l S -t 10  -i {Pre}_2.report.txt", shell=True, stdout=kkf, stderr=kkf)
                    _run_bracken_sub(krdb, f"{Pre}_2.report.txt", f"{Pre}_2", kkf)
            elif has_fq1:
                subprocess.run(f"kraken2 --db {krdb} --threads {threads} --output {Pre}_2.out.txt --report {Pre}_2.report.txt {fq1_path}", shell=True, stdout=kkf, stderr=kkf)
                if use_bracken:
                    subprocess.run(f"bracken -d {krdb} -o {Pre}_2.bracken1.txt -w {Pre}_2.bracken2.txt -l S -t 10  -i {Pre}_2.report.txt", shell=True, stdout=kkf, stderr=kkf)
                    _run_bracken_sub(krdb, f"{Pre}_2.report.txt", f"{Pre}_2", kkf)

    pre2 = f"{Pre}_2" if (has_fq1 or has_fq2) else Pre
    if lelID != "nolevel":
        level = lelID.split(",")[1]
        krakenfile = f"{pre2}.report.txt"
        tkid = lelID.split(",")[0]
        if _extracted_reads_ready(tkid, bool(fq2)):
            fq1, fq2 = _reuse_extracted_reads(tkid, bool(fq2))
        else:
            taxlist1 = [int(i) for i in proc_kra1(krakenfile, tkid, level)]
            if fq2:
                exreadsID1(taxlist1, f"{pre2}.out.txt", f"{Pre}.R1.fastq.gz", f"{Pre}.R2.fastq.gz")
                fq1 = f"{tkid}.1.fastq"
                fq2 = f"{tkid}.2.fastq"
            else:
                exreadsID1(taxlist1, f"{pre2}.out.txt", f"{Pre}.R1.fastq.gz", 0)
                fq1 = f"{tkid}.1.fastq"
                fq2 = 0
    elif method == "meta":
        if _extracted_reads_ready(default_taxid, bool(fq2)):
            fq1, fq2 = _reuse_extracted_reads(default_taxid, bool(fq2))
        else:
            taxlist1 = [int(i) for i in proc_kra1(f"{pre2}.report.txt", default_taxid, "D")]
            if fq2:
                exreadsID1(taxlist1, f"{pre2}.out.txt", f"{Pre}.R1.fastq.gz", f"{Pre}.R2.fastq.gz")
                fq1 = f"{default_taxid}.1.fastq"
                fq2 = f"{default_taxid}.2.fastq"
            else:
                exreadsID1(taxlist1, f"{pre2}.out.txt", f"{Pre}.R1.fastq.gz", 0)
                fq1 = f"{default_taxid}.1.fastq"
                fq2 = 0

    return fq1, fq2, pre2


def map_assembly_reads(finalfa, inf, fq1, fq2, threads, Pre, asmt, long_type, logf):
    inf, fq1, fq2 = _normalize_assembly_inputs(inf, fq1, fq2)
    if asmt in ["longasm", "longref"]:
        if not inf:
            raise FileNotFoundError(f"{Pre} 缺少可用长读输入文件，无法回贴组装 reads")
        if long_type == "Nanopore":
            run_command(f"minimap2 -ax map-ont {finalfa} {inf} -t {threads} |samtools sort -o {Pre}.sorted.bam", logf=logf)
        elif long_type == "PacBio_CLR":
            run_command(f"minimap2 -ax map-pb {finalfa} {inf} -t {threads} |samtools sort -o {Pre}.sorted.bam", logf=logf)
        elif long_type == "PacBio_CCS":
            run_command(f"minimap2 -ax map-hifi {finalfa} {inf} -t {threads} |samtools sort -o {Pre}.sorted.bam", logf=logf)
    elif asmt in ["shortasm", "shortref"]:
        run_command(f"bwa index {finalfa}", logf=logf)
        if fq2:
            run_command(f"minimap2 -ax sr {finalfa} {fq1} {fq2} -t {threads}|samtools sort -o {Pre}_ngs.sorted.bam", logf=logf)
        else:
            run_command(f"minimap2 -ax sr {finalfa} {fq1} -t {threads}|samtools sort -o {Pre}_ngs.sorted.bam", logf=logf)
    else:
        if not inf:
            raise FileNotFoundError(f"{Pre} 缺少可用长读输入文件，无法回贴组装 reads")
        if long_type == "Nanopore":
            run_command(f"minimap2 -ax map-ont {finalfa} {inf} -t {threads} |samtools sort -o {Pre}.sorted.bam", logf=logf)
        elif long_type == "PacBio_CLR":
            run_command(f"minimap2 -ax map-pb {finalfa} {inf} -t {threads} |samtools sort -o {Pre}.sorted.bam", logf=logf)
        elif long_type == "PacBio_CCS":
            run_command(f"minimap2 -ax map-hifi {finalfa} {inf} -t {threads} |samtools sort -o {Pre}.sorted.bam", logf=logf)
        if fq2:
            run_command(f"bwa mem {finalfa} {fq1} {fq2} -t {threads}|samtools sort -o {Pre}_ngs.sorted.bam", logf=logf)
        elif fq1:
            run_command(f"bwa mem {finalfa} {fq1} -t {threads}|samtools sort -o {Pre}_ngs.sorted.bam", logf=logf)

    _generate_depth_files(Pre, logf)


def build_assembly_info(finalfa, Pre, method):
    canudb = None
    if method == "canu":
        subprocess.run("seqkit fx2tab canu/canu.contigs.fasta -n > canu.txt", shell=True)
        canudb = pd.read_table("canu.txt", sep=" ", header=None)
        canudb["len"] = canudb[1].str.replace("len=", "").str.strip().astype("int")
        canudb = canudb.sort_values("len", ascending=False)
        canudb["index"] = canudb.index + 1
        canudb["contig"] = "contig_" + canudb["index"].astype("str")
        canudb["cir"] = canudb[6].str.replace("suggestCircular=", "")

    if not os.path.exists("flye_output"):
        os.makedirs("flye_output")
        subprocess.run(f"seqkit fx2tab {finalfa} -n -l > flye_output/tmp.stat.txt", shell=True)
        flyedb = pd.read_table("flye_output/tmp.stat.txt", header=None)
        flyedb["是否成环"] = "-"
        flyedb.rename(columns={0: "序列名称", 1: "序列长度"}, inplace=True)
        if method == "canu" and canudb is not None:
            for ctg in flyedb["序列名称"].tolist():
                flyedb.loc[flyedb["序列名称"] == ctg, "是否成环"] = canudb.loc[canudb["contig"] == ctg, "cir"].tolist()[0]

        mosdb = pd.read_table(f"{Pre}.regions.bed", header=None) if os.path.isfile(f"{Pre}.regions.bed") else pd.read_table(f"{Pre}_ngs.regions.bed", header=None)
        flyedb = mosdb.groupby(0).agg("mean").merge(flyedb, left_on=0, right_on="序列名称")
        flyedb.rename(columns={3: "平均深度"}, inplace=True)
        flyedb = flyedb[["序列名称", "序列长度", "平均深度", "是否成环"]]
        flyedb["平均深度"] = flyedb["平均深度"].round()
        flyedb.sort_values("序列长度", axis=0, inplace=True, ascending=False)
        flyedb.to_csv("flye_output/assembly_info.txt", sep="\t", index=False)
        return

    transdb = pd.read_table("transname.tsv", header=None, names=["oldname", "newname"])
    subprocess.run(f"seqkit fx2tab {finalfa} -n -l > flye_output/tmp.stat.txt", shell=True)
    newcontigdb = pd.read_table("flye_output/tmp.stat.txt", header=None)
    flyedb = pd.read_table("flye_output/assembly_info.txt")
    flyedb = flyedb.merge(transdb, left_on="序列名称", right_on="oldname")
    flyedb = flyedb.merge(newcontigdb, left_on="newname", right_on=0)
    flyedb["序列名称"] = flyedb["newname"]
    flyedb["序列长度"] = flyedb[1]
    flyedb = flyedb.sort_values("序列长度", ascending=False)
    flyedb = flyedb[["序列名称", "序列长度", "平均深度", "是否成环"]]
    flyedb.to_csv("flye_output/assembly_info.txt", sep="\t", index=False)


def write_finalfasta_stats(finalfa, Pre):
    subprocess.run(f"seqkit stat -a -T -G N {finalfa} > ragtag_sum.tsv", shell=True)
    ragtagdb = pd.read_table("ragtag_sum.tsv")
    ragtagdb = ragtagdb[["num_seqs", "sum_len", "max_len", "min_len", "sum_gap", "N50", "GC(%)"]]
    ragtagdb["样本名称"] = Pre
    ragtagdb.rename(columns={"num_seqs": "contig数量", "sum_len": "总长度", "min_len": "最小contig长度", "max_len": "最大contig长度"}, inplace=True)
    ragtagdb["N比例(%)"] = ((ragtagdb["sum_gap"] / ragtagdb["总长度"]) * 100).round(2)
    ragtagdb = ragtagdb[["样本名称", "contig数量", "总长度", "最大contig长度", "最小contig长度", "N50", "GC(%)", "N比例(%)"]]
    ragtagdb.to_csv("finalfasta.tsv", index=False, sep="\t")


def build_fasta_assembly_info(finalfa, Pre, method):
    canudb = None
    if method == "canu":
        subprocess.run("seqkit fx2tab canu/canu.contigs.fasta -n > canu.txt", shell=True)
        canudb = pd.read_table("canu.txt", sep=" ", header=None)
        canudb["len"] = canudb[1].str.replace("len=", "").str.strip().astype("int")
        canudb = canudb.sort_values("len", ascending=False)
        canudb["index"] = canudb.index + 1
        canudb["contig"] = "contig_" + canudb["index"].astype("str")
        canudb["cir"] = canudb[6].str.replace("suggestCircular=", "")

    if not os.path.exists("flye_output"):
        os.makedirs("flye_output")
    subprocess.run(f"seqkit fx2tab {finalfa} -n -l > flye_output/tmp.stat.txt", shell=True)
    flyedb = pd.read_table("flye_output/tmp.stat.txt", header=None)
    flyedb["是否成环"] = "-"
    flyedb.rename(columns={0: "序列名称", 1: "序列长度"}, inplace=True)
    if method == "canu" and canudb is not None:
        for ctg in flyedb["序列名称"].tolist():
            flyedb.loc[flyedb["序列名称"] == ctg, "是否成环"] = canudb.loc[canudb["contig"] == ctg, "cir"].tolist()[0]
    flyedb["平均深度"] = "-"
    flyedb = flyedb[["序列名称", "序列长度", "平均深度", "是否成环"]]
    flyedb.sort_values("序列长度", axis=0, inplace=True, ascending=False)
    flyedb.to_csv("flye_output/assembly_info.txt", sep="\t", index=False)
    return flyedb


def annotate_with_prokka(finalfa, Pre, threads, log_path="asb.log"):
    with open(log_path, "a") as logf:
        subprocess.run(f"prokka --force --outdir {Pre}_prokka --prefix {Pre} --addgenes --cpus {threads} {finalfa}", shell=True, stdout=logf, stderr=logf)

    prkskpn = int(os.popen(f"grep '##sequence-region' {Pre}_prokka/{Pre}.gff|wc -l").read()) + 1
    prkfan = int(os.popen(f"grep -n '#FASTA' {Pre}_prokka/{Pre}.gff").read().split(":")[0])
    prkskfn = int(os.popen(f"cat {Pre}_prokka/{Pre}.gff|wc -l").read()) - prkfan + 1
    prokkadb = pd.read_table(
        f"{Pre}_prokka/{Pre}.gff",
        skiprows=prkskpn,
        skipfooter=prkskfn,
        engine="python",
        header=None,
        names=["染色体", "数据库", "类型", "起始位置", "终止位置", "t1", "链方向", "t2", "注释1"],
    )
    prokkadb = prokkadb[["染色体", "类型", "起始位置", "终止位置", "链方向", "注释1"]]
    prokkadb["注释1"] = prokkadb["注释1"].fillna("").astype(str)
    prokkadb["基因名称"] = prokkadb["注释1"].str.extract(r"Name=([^;]+)")
    prokkadb["locus标签"] = prokkadb["注释1"].str.extract(r"ID=([^;]+)")
    prokkadb.fillna("-", inplace=True)
    prokkadb = prokkadb[["染色体", "类型", "起始位置", "终止位置", "链方向", "基因名称", "locus标签"]]
    prokkadb.to_csv(f"{Pre}.prokka.tsv", sep="\t", index=False)

    monthchin = int(os.popen(f"grep '月' {Pre}_prokka/{Pre}.gbk|wc -l").read().strip())
    if monthchin:
        if os.popen(f"head -n 1 {Pre}_prokka/{Pre}.gbk").read().strip()[72] == "月":
            subprocess.run(f"sed -i 's/月/  /g' {Pre}_prokka/{Pre}.gbk", shell=True)
        else:
            subprocess.run(f"sed -i 's/月/ /g' {Pre}_prokka/{Pre}.gbk", shell=True)
    subprocess.run(f"cp {Pre}_prokka/{Pre}.gbk tt.gbk", shell=True)


def enhance_plasmid_results(finalfa, Pre, log_path="asb.log"):
    staramr_db = shlex.quote(_staramr_database_path())
    with open(log_path, "a") as logf:
        subprocess.run(conda_run_command("plasflow", f"PlasFlow.py --input {finalfa} --output {Pre}_plaspredict.tsv"), shell=True, stdout=logf, stderr=logf)
        subprocess.run(f"staramr search -d {staramr_db} {shlex.quote(str(finalfa))} -o staramr_result -n 30", shell=True, stdout=logf, stderr=logf)

    plasmiddb = pd.read_table("staramr_result/plasmidfinder.tsv")
    rawindexname = plasmiddb.index.name
    if plasmiddb.shape[0] > 1:
        plasmiddb = plasmiddb.groupby("Contig").apply(_join_plasmid)
        plasmiddb.index.name = rawindexname

    plasflowdb = pd.read_table(f"{Pre}_plaspredict.tsv")[["contig_name", "label"]]
    flyedb = pd.read_table("flye_output/assembly_info.txt")
    flyedb = flyedb.merge(plasflowdb, left_on="序列名称", right_on="contig_name").drop("contig_name", axis=1)
    flyedb = flyedb.merge(plasmiddb, left_on="序列名称", right_on="Contig", how="left").drop("Contig", axis=1)
    flyedb = flyedb.rename(columns={"label": "基因组/质粒", "Plasmid": "质粒分型"})
    flyedb.fillna("-", inplace=True)
    flyedb.to_csv("tmp1.tsv", sep="\t", index=False)
    flyedb["序列长度"] = flyedb["序列长度"].astype("int")
    plasmidlist = flyedb.loc[
        (flyedb["基因组/质粒"].str.contains("plasmid")) | ((flyedb["质粒分型"] != "-") & (flyedb["序列长度"] < 1000000)),
        "序列名称",
    ].tolist()
    flyedb["占比"] = (flyedb["序列长度"] / flyedb["序列长度"].sum()).round(2)
    flyedb = flyedb[["序列名称", "序列长度", "平均深度", "是否成环", "基因组/质粒", "质粒分型", "占比"]]
    flyedb.to_csv("flye_output/assembly_info.txt", sep="\t", index=False)
    return flyedb, plasmidlist


def export_plasmid_gbk_and_cgview(Pre, plasmidlist):
    input_gbk = f"{Pre}_prokka/{Pre}.gbk"
    if len(plasmidlist) != 0:
        for contig_id_to_extract in plasmidlist:
            for record in SeqIO.parse(input_gbk, "genbank"):
                if record.id == contig_id_to_extract:
                    with open(f"{Pre}_prokka/{contig_id_to_extract}.gbk", "w") as output_handle:
                        SeqIO.write(record, output_handle, "genbank")
                    with open("cgview.log", "a") as cgvf:
                        subprocess.run(f"ruby /data/deploy/bio-elite/bio/script/cgview_builder_cli.rb -s {Pre}_prokka/{contig_id_to_extract}.gbk -o {contig_id_to_extract}.json -c /data/deploy/bio-elite/bio/script/CGview.yaml -n {contig_id_to_extract}", shell=True, stdout=cgvf, stderr=cgvf)

        records = list(SeqIO.parse(input_gbk, "genbank"))
        filtered_records = [record for record in records if record.id not in plasmidlist]
        with open(f"{Pre}_prokka/main.gbk", "w") as output_handle:
            SeqIO.write(filtered_records, output_handle, "genbank")
        with open("cgview.log", "a") as cgvf:
            subprocess.run(f"ruby /data/deploy/bio-elite/bio/script/cgview_builder_cli.rb -s {Pre}_prokka/main.gbk -o main.json -c /data/deploy/bio-elite/bio/script/CGview.yaml -n main", shell=True, stdout=cgvf, stderr=cgvf)
    else:
        with open("cgview.log", "a") as cgvf:
            subprocess.run(f"ruby /data/deploy/bio-elite/bio/script/cgview_builder_cli.rb -s {Pre}_prokka/{Pre}.gbk -o main.json -c /data/deploy/bio-elite/bio/script/CGview.yaml -n {Pre}", shell=True, stdout=cgvf, stderr=cgvf)


def run_checkm_and_write_summary(Pre, finalfa, threads, species, flyedb):
    if not os.path.isdir(f"{Pre}_bin_genome_out"):
        os.makedirs(f"{Pre}_bin_genome_out")
    subprocess.run(f"cp {finalfa} {Pre}_bin_genome_out/{Pre}.fna", shell=True)

    checkm2_db = _checkm2_database_path()
    with open("checkm2", "w") as cmf:
        subprocess.run(conda_run_command("cm210", f"checkm2 predict -i  {Pre}_bin_genome_out -x fna -o checkm2_out -t {threads} --force --database_path {shlex.quote(checkm2_db)}"), shell=True, stdout=cmf, stderr=cmf)
        if os.path.isfile("checkm2_out/quality_report.tsv"):
            checkmdb = pd.read_table("checkm2_out/quality_report.tsv")
            checkmdb["样本名称"] = Pre
            checkmdb.rename(columns={"Completeness": "完整性", "Contamination": "污染率"}, inplace=True)
            checkmdb["物种名称"] = species
            checkmdb[["样本名称", "物种名称", "污染率", "完整性"]].to_csv(f"{Pre}.checkm.tsv", sep="\t", index=False)

    subprocess.run(f"seqkit stat {finalfa} -b -T -a > {Pre}.fasum.tsv", shell=True)
    fadb = pd.read_table(f"{Pre}.fasum.tsv")
    assdict = {
        "Contig数量": fadb["num_seqs"].tolist()[0],
        "N50长度": fadb["N50"].tolist()[0],
        "组装基因组长度": fadb["sum_len"].tolist()[0],
        "最长片段长度": fadb["max_len"].tolist()[0],
        "污染率": "-",
        "完整性": "-",
        "主基因组是否成环": "-",
    }
    if os.path.isfile(f"{Pre}.checkm.tsv"):
        cdb = pd.read_table(f"{Pre}.checkm.tsv")
        assdict["污染率"] = round(float(cdb["污染率"].tolist()[0]), 2)
        assdict["完整性"] = round(float(cdb["完整性"].tolist()[0]), 2)
    assdb = pd.DataFrame(assdict, index=[0])
    assdb["样本名称"] = Pre
    assdb = assdb[["样本名称", "Contig数量", "N50长度", "最长片段长度", "主基因组是否成环", "污染率", "完整性"]]
    assdb["主基因组是否成环"] = flyedb["是否成环"].tolist()[0]
    assdb.to_csv(f"{Pre}.assemble.result.tsv", sep="\t", index=False)


def run_genovi_summary(Pre):
    with open("genovi.log", "w") as fv:
        try:
            subprocess.run(conda_run_command("genovi", f"genovi -i {Pre}_prokka/{Pre}.gbk -o {Pre}_genovi -s draft"), shell=True, stdout=fv, stderr=fv)
            cogdb = pd.read_table(f"{Pre}_genovi/{Pre}_genovi_COG_Classification.csv", sep=",", header=1)
            cogdb.iloc[-2, :].to_csv("Cog_summary.tsv", sep="\t", header=False)
        except Exception:
            pass


def write_gene_summaries(Pre):
    tmpgenedb = pd.read_table(f"{Pre}_prokka/{Pre}.tsv")
    tmpgenedb = tmpgenedb[tmpgenedb.ftype == "gene"]
    tmpgenedb_dict = {}
    for minlg in range(0, 2000, 100):
        tmpgenedb_dict[f"{minlg}-{minlg+100}"] = sum(tmpgenedb.length_bp.between(minlg, minlg + 100))
    tmpgenedb_dict[">2000"] = sum(tmpgenedb.length_bp >= 2000)
    genedb = pd.DataFrame(tmpgenedb_dict, index=["Gene数量"]).T
    genedb["范围"] = genedb.index
    genedb = genedb[["范围", "Gene数量"]]
    genedb.to_csv(f"{Pre}_gene_raw_sum.tsv", sep="\t", index=False)

    gene_fundb = pd.read_table(f"{Pre}_prokka/{Pre}.txt", sep=":")
    gene_fundb.index = gene_fundb.organism
    gene_fundb = gene_fundb.drop("organism", axis=1).T
    gene_fundb.to_csv(f"{Pre}.genefun_summary.tsv", sep="\t", index=False)

    if os.path.isfile(f"{Pre}.uniqgene.fasta"):
        subprocess.run(f"seqkit fx2tab -n -l {Pre}.uniqgene.fasta > {Pre}.uniqgene.tsv", shell=True)
        tmpgeneudb = pd.read_table(f"{Pre}.uniqgene.tsv", names=["Genename", "length"])
        tmpgeneudb_dict = {}
        for minlg in range(0, 2000, 100):
            tmpgeneudb_dict[f"{minlg}-{minlg+100}"] = sum(tmpgeneudb.length.between(minlg, minlg + 100))
        tmpgeneudb_dict[">2000"] = sum(tmpgeneudb.length >= 2000)
        geneudb = pd.DataFrame(tmpgeneudb_dict, index=["Gene数量"]).T
        geneudb["范围"] = geneudb.index
        geneudb = geneudb[["范围", "Gene数量"]]
        geneudb.to_csv(f"{Pre}_gene_uniq_sum.tsv", sep="\t", index=False)

    subprocess.run(f"grep CDS {Pre}_prokka/{Pre}.gff > {Pre}.CDS.gff", shell=True)
    with open(f"{Pre}.CDS.gff") as f:
        open(f"{Pre}.Contig_gene.tsv", "w").write("Contig\tCDS\n")
        for line in f:
            line = line.strip().split("\t")
            if line[2] == "CDS":
                contig = line[0]
                cdsid = line[8].split("locus_tag=")[-1].split(";")[0]
                open(f"{Pre}.Contig_gene.tsv", "a").write(f"{contig}\t{cdsid}\n")

    subprocess.run(f"grep 'gene' {Pre}_prokka/{Pre}.gff > {Pre}.gene.gff", shell=True)
    with open(f"{Pre}.gene.gff") as f:
        open(f"{Pre}.gene.bed", "w").write("")
        for line in f:
            line = line.strip().split("\t")
            if line[2] == "gene":
                contig = line[0]
                cdsid = line[8].split("locus_tag=")[-1].split(";")[0]
                start = line[3]
                end = line[4]
                open(f"{Pre}.gene.bed", "a").write(f"{contig}\t{start}\t{end}\t{cdsid}\n")
    subprocess.run(f"samtools faidx {Pre}_prokka/{Pre}.fna", shell=True)
    subprocess.run(f"bedtools getfasta -fi {Pre}_prokka/{Pre}.fna -bed {Pre}.gene.bed -name > tmp_{Pre}.gene.fasta", shell=True)
    subprocess.run(f"cut -d ':' -f1  tmp_{Pre}.gene.fasta > {Pre}.gene.fasta", shell=True)


def getvfID(tmpdb):
    matches = re.findall(r"VF\d+", tmpdb["产物"])
    if len(matches) > 0:
        if len(matches) > 1:
            return "|".join([i for i in matches])
        return matches[0]
    return "-"


def rgi_fun(Pre):
    subprocess.run(conda_run_command("amr_aux", f"rgi main -i {Pre}.final.fasta -o {Pre}.rgi --clean --include_loose --low_quality -g PYRODIGAL"), shell=True)
    rgidb = pd.read_table(f"{Pre}.rgi.txt")
    rgidb = rgidb[["Contig", "Start", "Stop", "Orientation", "Cut_Off", "Pass_Bitscore", "Best_Hit_ARO", "Best_Identities", "Model_type", "SNPs_in_Best_Hit_ARO", "AMR Gene Family", "Drug Class"]]
    rgidb["序列名称"] = rgidb["Contig"].str.split("_").str[:2].str.join("_")
    rgidb.rename(columns={"Start": "序列起始", "Stop": "序列终止", "Orientation": "正负链", "Cut_Off": "过滤标准", "Pass_Bitscore": "比对得分", "Best_Hit_ARO": "耐药基因", "Best_Identities": "一致性%", "Model_type": "基因数据库", "SNPs_in_Best_Hit_ARO": "基因突变", "AMR Gene Family": "耐药基因家族", "Drug Class": "药物类别"}, inplace=True)
    rgidb.drop("Contig", axis=1, inplace=True)
    rgidb = rgidb[["序列名称", "序列起始", "序列终止", "正负链", "过滤标准", "比对得分", "耐药基因", "一致性%", "基因数据库", "基因突变", "耐药基因家族", "药物类别"]]
    rgidb.sort_values("比对得分", ascending=False, inplace=True)
    rgidb.to_csv(f"{Pre}.rgi.tsv", sep="\t", index=False)
    rgidb.to_csv(f"{Pre}.rgi.bed", sep="\t", index=False, header=False)


def outfun(Pre, x, typeF):
    ts = x["start"].min()
    te = x["end"].max()
    ofname = x["GeneName"].tolist()[0]
    file_token = _safe_path_token(ofname)
    Spename = x["Species"].tolist()[0]
    x["start"] = x.reset_index().index + 1
    x["end"] = x.reset_index().index + 2
    tsv_path = Path("geneDepth") / f"{file_token}_{typeF}.tsv"
    bed_path = Path("geneDepth") / f"{file_token}_{typeF}.bed"
    fasta_path = Path("geneDepth") / f"{file_token}_{typeF}.fasta"
    x[["Chrom", "start", "end", "Depth"]].to_csv(tsv_path, sep="\t", header=False, index=False)
    bed_path.write_text(f"{x['Chrom'].tolist()[0]}\t{ts}\t{te}\t{Pre}_{ofname}\t{Pre}_{ofname}\t{x['strand'].tolist()[0]}")
    subprocess.run(f"bedtools getfasta -fi {shlex.quote(f'{Pre}.final.fasta')} -bed {shlex.quote(str(bed_path))} -name -s > {shlex.quote(str(fasta_path))}", shell=True)
    tmpdict = {"片段名称": x["Chrom"].tolist()[0], "物种名称": Spename, "起始位置": x["start"].min(), "终止位置": x["end"].max(), "覆盖度(>0)%": round(x[x["Depth"] > 0].shape[0] / x.shape[0], 4) * 100, "覆盖度(>10)%": round(x[x["Depth"] > 10].shape[0] / x.shape[0], 4) * 100, "覆盖度(>100)%": round(x[x["Depth"] > 100].shape[0] / x.shape[0], 4) * 100, "平均深度": round(x["Depth"].mean(), 2), "最低深度": x["Depth"].min(), "最高深度": x["Depth"].max()}
    return pd.DataFrame(tmpdict, index=[0]).round(2)


def _run_bracken_sub(krdb, report_path, prefix, kkf):
    testbrkdb = pd.read_table(report_path, header=None)
    if "S4" in testbrkdb[3]:
        level = "S3"
    elif "S3" in testbrkdb[3]:
        level = "S2"
    else:
        level = "S1"
    subprocess.run(f"bracken -d {krdb} -o {prefix}_Sub.bracken1.txt -w {prefix}_Sub.bracken2.txt -l {level} -t 10  -i {report_path}", shell=True, stdout=kkf, stderr=kkf)


def _write_contig_beds(prefix):
    if not os.path.isdir("Contigbedfile"):
        os.makedirs("Contigbedfile")
    contigbed = pd.read_table(f"{prefix}.regions.bed", header=None)
    contig1bed = pd.read_table(f"{prefix}_1.regions.bed", header=None)
    contigbed.groupby(0).apply(lambda x: x.to_csv(f"Contigbedfile/{x[0].tolist()[0]}.bed", index=False, header=False, sep="\t"))
    contig1bed.groupby(0).apply(lambda x: x.to_csv(f"Contigbedfile/{x[0].tolist()[0]}_dis1.bed", index=False, header=False, sep="\t"))
    for bed in os.listdir("Contigbedfile"):
        if not bed.endswith("_dis1.bed"):
            ttPre = bed.replace(".bed", "")
            if int(os.popen(f"cat Contigbedfile/{bed}|wc -l").read()) < 10:
                subprocess.run(f"mv Contigbedfile/{ttPre}_dis1.bed Contigbedfile/{ttPre}.bed ", shell=True)
    contig1bed[contig1bed[3] == 0].to_csv("mask.bed", sep="\t", index=False, header=False)


def _generate_depth_files(Pre, logf):
    runtime = get_runtime_context()
    if runtime.analysis_target != 'virus':
        if os.path.isfile(f"{Pre}.sorted.bam"):
            subprocess.run(f"samtools index {Pre}.sorted.bam", shell=True, stdout=logf, stderr=logf)
            subprocess.run(f"mosdepth -b 1000 {Pre} {Pre}.sorted.bam", shell=True)
            subprocess.run(f"mosdepth -b 1 {Pre}_1 {Pre}.sorted.bam", shell=True)
            subprocess.run(f"gunzip -f {Pre}.regions.bed.gz", shell=True)
            subprocess.run(f"gunzip -f {Pre}_1.regions.bed.gz", shell=True)
            _write_contig_beds(Pre)

        if os.path.isfile(f"{Pre}_ngs.sorted.bam"):
            subprocess.run(f"samtools index {Pre}_ngs.sorted.bam", shell=True, stdout=logf, stderr=logf)
            subprocess.run(f"mosdepth -b 1000 {Pre}_ngs {Pre}_ngs.sorted.bam", shell=True)
            subprocess.run(f"mosdepth -b 1 {Pre}_ngs_1 {Pre}_ngs.sorted.bam", shell=True)
            subprocess.run(f"gunzip -f {Pre}_ngs.regions.bed.gz", shell=True)
            subprocess.run(f"gunzip -f {Pre}_ngs_1.regions.bed.gz", shell=True)
            _write_contig_beds(f"{Pre}_ngs")
    else:
        if os.path.isfile(f"{Pre}.sorted.bam"):
            subprocess.run(f'samtools depth -aa {Pre}.sorted.bam > {Pre}.pre.regions.bed',shell=True)
            tmppd = pd.read_table(f'{Pre}.pre.regions.bed',header=None)
            tmppd['oldpos'] = tmppd[1]-1
            tmppd[[0,'oldpos',1,2]].to_csv(f'{Pre}.regions.bed',index=False,sep='\t',header=False)
        else:
            subprocess.run(f'samtools depth -aa {Pre}_ngs.sorted.bam > {Pre}.pre.regions.bed',shell=True)
            tmppd = pd.read_table(f'{Pre}.pre.regions.bed',header=None)
            tmppd['oldpos'] = tmppd[1]-1
            tmppd[[0,'oldpos',1,2]].to_csv(f'{Pre}.regions.bed',index=False,sep='\t',header=False)



def _join_plasmid(group):
    plasmids = "|".join(group["Plasmid"])
    newdb = group
    newdb["Plasmid"] = plasmids
    return newdb.iloc[0, :]
