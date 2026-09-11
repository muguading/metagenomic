#!/usr/bin/env python3
"""Backfill missing meta_plas_vf_card.tsv files under fastq_analysis.

Example:
    python scripts/backfill_meta_plas_vf_card.py /path/to/output/fastq_analysis
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class BackfillTask:
    sample: str
    sample_dir: Path
    status: str
    contigs: Path | None = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "扫描 output/fastq_analysis 下缺失 meta_plas_vf_card.tsv 的宏基因组样本目录，"
            "复用已有中间结果，并补算缺失的 PlasFlow/staramr/GTDB-Tk/CheckM2/CoverM 后生成该表。"
        )
    )
    parser.add_argument(
        "fastq_analysis",
        type=Path,
        help="输出目录中的 fastq_analysis 路径，例如 /data/run1/fastq_analysis。",
    )
    parser.add_argument(
        "--sample",
        help="只处理指定样本目录名；默认扫描 fastq_analysis 下所有样本目录。",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="即使 meta_plas_vf_card.tsv 已存在也重新生成。",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="只列出将处理/跳过的样本，不运行分析或写入文件。",
    )
    parser.add_argument(
        "--reuse-only",
        action="store_true",
        help="只根据现有结果拼表，不补跑 PlasFlow/staramr/GTDB-Tk/CheckM2/CoverM。",
    )
    return parser.parse_args()


def _has_meta_inputs(sample_dir: Path) -> bool:
    return (
        (sample_dir / "megahit_output" / "final.contigs.fa").is_file()
        or (sample_dir / "codex_mag_binning").is_dir()
        or (sample_dir / "BASALT_out").is_dir()
    )


def _sample_dirs(fastq_analysis: Path, sample: str | None) -> list[Path]:
    if sample:
        candidate = fastq_analysis / sample
        if not candidate.is_dir():
            raise FileNotFoundError(f"样本目录不存在：{candidate}")
        return [candidate]
    if _has_meta_inputs(fastq_analysis):
        return [fastq_analysis]
    ignored = {"barout", "basecaller_outputs", "__pycache__"}
    return sorted(
        path
        for path in fastq_analysis.iterdir()
        if path.is_dir() and path.name not in ignored and _has_meta_inputs(path)
    )


def _find_filtered_contigs(sample_dir: Path, sample: str) -> Path | None:
    roots = [
        sample_dir / "codex_mag_binning" / sample / "filtered_contigs",
        sample_dir / "codex_mag_binning",
    ]
    candidates: list[Path] = []
    for root in roots:
        if not root.is_dir():
            continue
        candidates.extend(root.rglob("*.min1500.fa"))
        candidates.extend(root.rglob("*.fa"))
        candidates.extend(root.rglob("*.fasta"))
        candidates.extend(root.rglob("*.fna"))
    candidates = [
        path
        for path in candidates
        if path.is_file() and path.stat().st_size > 0 and "binning_genomes" not in path.parts
    ]
    if not candidates:
        return None
    return sorted(candidates, key=lambda path: (".min1500." not in path.name, len(path.parts), path.name))[0]


def _ensure_final_contigs(sample_dir: Path, sample: str) -> Path | None:
    final_contigs = sample_dir / "megahit_output" / "final.contigs.fa"
    if final_contigs.is_file() and final_contigs.stat().st_size > 0:
        return final_contigs
    fallback = _find_filtered_contigs(sample_dir, sample)
    if fallback is None:
        return None
    final_contigs.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(fallback, final_contigs)
    return final_contigs


def _ensure_legacy_derep_input(sample_dir: Path, sample: str, contigs: Path) -> Path:
    derep_dir = sample_dir / "BASALT_out" / "meta_drep_out" / "dereplicated_genomes"
    derep_dir.mkdir(parents=True, exist_ok=True)
    has_fasta = any(
        path.is_file() and path.suffix.lower() in {".fa", ".fasta", ".fna"}
        for path in derep_dir.iterdir()
    )
    if has_fasta:
        return derep_dir
    shutil.copy2(contigs, derep_dir / f"{sample}.fa")
    return derep_dir


def discover_tasks(fastq_analysis: Path, sample: str | None, force: bool) -> list[BackfillTask]:
    fastq_analysis = fastq_analysis.expanduser().resolve()
    if not fastq_analysis.is_dir():
        raise FileNotFoundError(f"fastq_analysis 路径不存在：{fastq_analysis}")
    tasks: list[BackfillTask] = []
    for sample_dir in _sample_dirs(fastq_analysis, sample):
        sample_name = sample or sample_dir.name
        output = sample_dir / "meta_plas_vf_card.tsv"
        if output.is_file() and output.stat().st_size > 0 and not force:
            tasks.append(BackfillTask(sample_name, sample_dir, "skip_exists", output))
            continue
        contigs = _ensure_final_contigs(sample_dir, sample_name)
        if contigs is None:
            tasks.append(BackfillTask(sample_name, sample_dir, "skip_no_contigs"))
            continue
        tasks.append(BackfillTask(sample_name, sample_dir, "run", contigs))
    return tasks


def run_task(task: BackfillTask, reuse_only: bool = False) -> None:
    sys.path.insert(0, str(PROJECT_ROOT))
    from metagenomic_refactor.assembly import binning_result

    if task.contigs is not None:
        _ensure_legacy_derep_input(task.sample_dir, task.sample, task.contigs)
    old_cwd = Path.cwd()
    try:
        os.chdir(task.sample_dir)
        binning_result(
            task.sample,
            run_missing_tools=not reuse_only,
            progress=lambda message: print(f"  {task.sample} {message}", flush=True),
        )
    finally:
        os.chdir(old_cwd)
    output = task.sample_dir / "meta_plas_vf_card.tsv"
    if not output.is_file() or output.stat().st_size == 0:
        raise RuntimeError(f"未生成有效 meta_plas_vf_card.tsv：{output}")


def main() -> int:
    args = parse_args()
    tasks = discover_tasks(args.fastq_analysis, args.sample, args.force)
    if not tasks:
        print("未发现可处理的样本目录。")
        return 0

    runnable = [task for task in tasks if task.status == "run"]
    skipped_exists = [task for task in tasks if task.status == "skip_exists"]
    skipped_no_contigs = [task for task in tasks if task.status == "skip_no_contigs"]
    print(
        f"发现 {len(tasks)} 个样本目录：待补跑 {len(runnable)}，"
        f"已存在跳过 {len(skipped_exists)}，缺少 contigs 跳过 {len(skipped_no_contigs)}。"
    )
    for task in tasks:
        label = {
            "run": "补跑",
            "skip_exists": "已存在",
            "skip_no_contigs": "缺少contigs",
        }.get(task.status, task.status)
        print(f"[{label}] {task.sample}: {task.sample_dir}")
        if task.contigs:
            print(f"  contigs: {task.contigs}")
    if args.dry_run:
        return 0

    failures: list[tuple[BackfillTask, Exception]] = []
    for task in runnable:
        try:
            run_task(task, reuse_only=args.reuse_only)
            print(f"  完成：{task.sample_dir / 'meta_plas_vf_card.tsv'}")
        except Exception as error:
            failures.append((task, error))
            print(f"  失败：{task.sample}: {error}", file=sys.stderr)

    print(f"补跑完成：成功 {len(runnable) - len(failures)}，失败 {len(failures)}。")
    if failures:
        for task, error in failures:
            print(f"- {task.sample_dir}: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
