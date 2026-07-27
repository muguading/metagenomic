#!/usr/bin/env python3
"""Rerun MGE results for completed MAG tasks.

This script deliberately does *not* rerun binning, GTDB-Tk, AMR, plasmid, or
abundance analyses.  It finds completed task directories and reruns the
``genomad_mge.py`` MGE workflow in the same task directories.

CheckM2 backfilling is intentionally opt-in via ``--with-checkm2``.

Example:
    python scripts/run_binning_result.py --tasks-dir /path/to/completed_tasks
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
QUALITY_COLUMNS = ("Completeness", "Contamination")
OUTPUT_COLUMNS = {"Completeness": "完整性", "Contamination": "污染率"}
DEPLOY_ROOT_DEFAULT = Path(os.environ.get("META_DEPLOY_ROOT", PROJECT_ROOT))
GENOMAD_DB_DEFAULT = Path(os.environ.get("META_GENOMAD_DB", "database/genomad_db"))
MOBILEOG_DB_DEFAULT = Path(os.environ.get("META_MOBILEOG_DB", "database/beatrix/mobileOG-db.dmnd"))
MOBILEOG_META_DEFAULT = Path(
    os.environ.get("META_MOBILEOG_META", "database/beatrix/mobileOG-db-beatrix-1.6-All.csv")
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "在已完成 MAG 任务目录中重跑 genomad_mge.py；默认不处理 CheckM2，"
            "不改写 meta_plas_vf_card.tsv。"
        )
    )
    parser.add_argument(
        "--tasks-dir", type=Path, required=True,
        help="包含一个或多个已完成任务目录的根目录；会递归查找 meta_plas_vf_card.tsv。",
    )
    parser.add_argument(
        "--deploy-root", type=Path, default=DEPLOY_ROOT_DEFAULT,
        help="部署根目录；相对数据库路径都按此目录解析（默认：META_DEPLOY_ROOT 或脚本所在项目根目录）。",
    )
    parser.add_argument("--threads", type=int, default=10, help="MGE/可选 CheckM2 线程数（默认：10）。")
    parser.add_argument("--conda-exe", default=os.environ.get("CONDA_EXE", "conda"), help="conda 可执行文件。")
    parser.add_argument("--checkm2-env", default="cm210", help="CheckM2 conda 环境（默认：cm210）。")
    parser.add_argument("--mge-env", default=os.environ.get("META_MGE_ENV", "genomad_aux"), help="MGE/mobileOG conda 环境（默认：genomad_aux）。")
    parser.add_argument(
        "--checkm2-db", type=Path,
        default=Path(os.environ.get("META_CHECKM2_DB", PROJECT_ROOT / "database/checkm2/uniref100.KO.1.dmnd")),
        help="CheckM2 .dmnd 数据库；默认读取 META_CHECKM2_DB 或项目数据库。",
    )
    parser.add_argument(
        "--genomad-db", type=Path,
        default=GENOMAD_DB_DEFAULT,
        help="geNomad 数据库目录；默认 database/genomad_db，相对路径按 --deploy-root 解析。",
    )
    parser.add_argument(
        "--mobileog-db", type=Path,
        default=MOBILEOG_DB_DEFAULT,
        help="mobileOG DIAMOND 数据库；默认 database/beatrix/mobileOG-db.dmnd，相对路径按 --deploy-root 解析。",
    )
    parser.add_argument(
        "--mobileog-meta", type=Path,
        default=MOBILEOG_META_DEFAULT,
        help="mobileOG 注释表；默认 database/beatrix/mobileOG-db-beatrix-1.6-All.csv，相对路径按 --deploy-root 解析。",
    )
    parser.add_argument("--sample", help="单任务时手动指定样本前缀；批量任务通常自动推断。")
    parser.add_argument(
        "--with-checkm2", action="store_true",
        help="额外补跑/复用 CheckM2 并回填 meta_plas_vf_card.tsv；默认不启用。",
    )
    parser.add_argument("--skip-mge", action="store_true", help="不重跑 genomad_mge.py。")
    parser.add_argument("--force", action="store_true", help="配合 --with-checkm2 使用：即使已有质量报告也重新运行 CheckM2。")
    parser.add_argument(
        "--no-force-mge", action="store_true",
        help="MGE 默认会重跑并覆盖样本 geNomad 子目录；加此参数则复用已有 geNomad 汇总。",
    )
    parser.add_argument("--dry-run", action="store_true", help="只列出将处理的任务，不运行或写入文件。")
    parser.add_argument("--no-backup", action="store_true", help="回填前不备份原 meta_plas_vf_card.tsv。")
    return parser.parse_args()


def resolve_deploy_path(path: Path, deploy_root: Path) -> Path:
    path = path.expanduser()
    if path.is_absolute():
        return path.resolve()
    return (deploy_root / path).resolve()


def discover_tasks(root: Path) -> list[Path]:
    root = root.expanduser().resolve()
    if not root.is_dir():
        raise FileNotFoundError(f"任务根目录不存在：{root}")
    return sorted({table.parent for table in root.rglob("meta_plas_vf_card.tsv")})


def run_checkm2(task_dir: Path, args: argparse.Namespace, genome_dir: Path, quality_report: Path) -> None:
    command = [
        args.conda_exe, "run", "--no-capture-output", "-n", args.checkm2_env,
        "checkm2", "predict", "--thread", str(args.threads), "--input", str(genome_dir),
        "--output-directory", str(quality_report.parent), "-x", ".fa", "--force",
        "--database_path", str(args.checkm2_db),
    ]
    print(f"  + {' '.join(command)}", flush=True)
    with (task_dir / "bincheckm2.log").open("a", encoding="utf-8") as log:
        subprocess.run(command, cwd=task_dir, check=True, stdout=log, stderr=subprocess.STDOUT)
    if not quality_report.is_file() or quality_report.stat().st_size == 0:
        raise RuntimeError(f"CheckM2 未生成有效质量报告：{quality_report}")


def infer_sample_name(task_dir: Path, explicit_sample: str | None) -> str:
    if explicit_sample:
        return explicit_sample
    candidates: set[str] = set()
    patterns = {
        "*_meta_plaspredict.tsv": "_meta_plaspredict.tsv",
        "*.mge_risk_summary.tsv": ".mge_risk_summary.tsv",
        "*.integrated_mge_summary.tsv": ".integrated_mge_summary.tsv",
        "*.QC2.summary.tsv": ".QC2.summary.tsv",
        "*.final.json": ".final.json",
    }
    for pattern, suffix in patterns.items():
        for path in task_dir.glob(pattern):
            name = path.name
            if name.endswith(suffix):
                candidates.add(name[: -len(suffix)])
    if len(candidates) == 1:
        return next(iter(candidates))
    if not candidates:
        raise RuntimeError(
            f"无法自动推断样本前缀：{task_dir}。请确认存在 *_meta_plaspredict.tsv/"
            "*.mge_risk_summary.tsv，或用 --sample 指定。"
        )
    raise RuntimeError(f"样本前缀不唯一：{sorted(candidates)}。请用 --sample 指定。")


def run_mge(task_dir: Path, args: argparse.Namespace, sample: str) -> None:
    fasta_path = task_dir / "tmp_combine.fa"
    if not fasta_path.is_file():
        raise RuntimeError(f"缺少 MGE 输入 fasta：{fasta_path}")

    command = [
        sys.executable,
        "-c",
        (
            "from pathlib import Path\n"
            "from metagenomic_refactor.context import RuntimeContext, set_runtime_context\n"
            "from metagenomic_refactor.genomad_mge import "
            "GeNomadConfig, GeNomadSample, run_genomad_sample, _run_mobileog, "
            "integrate_mge_tables, summarize_mge_risk\n"
            "sample = __import__('os').environ['BACKFILL_MGE_SAMPLE']\n"
            "threads = int(__import__('os').environ['BACKFILL_MGE_THREADS'])\n"
            "cwd = Path.cwd()\n"
            "fasta = cwd / 'tmp_combine.fa'\n"
            "db = Path(__import__('os').environ['META_GENOMAD_DB'])\n"
            "env = __import__('os').environ.get('META_MGE_ENV', 'genomad_aux')\n"
            "force = __import__('os').environ.get('BACKFILL_MGE_FORCE') == '1'\n"
            "set_runtime_context(RuntimeContext(ofn=sample, runflow='', method='meta', rmhost='', tspeabun=''))\n"
            "cfg = GeNomadConfig(outdir=cwd, database=db, threads=threads, conda_env=env, force=force)\n"
            "run_genomad_sample(GeNomadSample(sample=sample, fasta=fasta), cfg)\n"
            "_run_mobileog(sample, fasta, threads, cwd, cfg)\n"
            "integrate_mge_tables(cwd, sample)\n"
            "summarize_mge_risk(sample, cwd, cfg)\n"
        ),
    ]
    env = os.environ.copy()
    env["META_GENOMAD_DB"] = str(args.genomad_db)
    env["META_MOBILEOG_DB"] = str(args.mobileog_db)
    env["META_MOBILEOG_META"] = str(args.mobileog_meta)
    env["META_MGE_ENV"] = args.mge_env
    env["BACKFILL_MGE_SAMPLE"] = sample
    env["BACKFILL_MGE_THREADS"] = str(args.threads)
    env["BACKFILL_MGE_FORCE"] = "0" if args.no_force_mge else "1"
    env["PYTHONPATH"] = os.pathsep.join(
        [str(PROJECT_ROOT), env["PYTHONPATH"]] if env.get("PYTHONPATH") else [str(PROJECT_ROOT)]
    )
    print(f"  + 在 {task_dir} 重跑 genomad_mge.py：sample={sample}", flush=True)
    with (task_dir / "genomad_mge.backfill.log").open("a", encoding="utf-8") as log:
        subprocess.run(command, cwd=task_dir, env=env, check=True, stdout=log, stderr=subprocess.STDOUT)
    for output_name in (f"{sample}.integrated_mge_summary.tsv", f"{sample}.mge_risk_summary.tsv"):
        output_path = task_dir / output_name
        if not output_path.is_file() or output_path.stat().st_size == 0:
            raise RuntimeError(f"MGE 未生成有效结果：{output_path}")


def merge_quality_report(task_dir: Path, quality_report: Path, backup: bool) -> tuple[int, int]:
    import pandas as pd

    result_path = task_dir / "meta_plas_vf_card.tsv"
    result = pd.read_table(result_path)
    quality = pd.read_table(quality_report)
    missing = {"Name", *QUALITY_COLUMNS}.difference(quality.columns)
    if missing:
        raise RuntimeError(f"CheckM2 报告缺少必要列 {sorted(missing)}：{quality_report}")
    if "Name" not in result.columns:
        raise RuntimeError(f"汇总表缺少 Name 列，无法匹配 MAG：{result_path}")

    quality = quality[["Name", *QUALITY_COLUMNS]].drop_duplicates("Name")
    quality = quality.rename(columns=OUTPUT_COLUMNS)
    result = result.drop(columns=[column for column in OUTPUT_COLUMNS.values() if column in result.columns])
    merged = result.merge(quality, on="Name", how="left")
    matched = int(merged["完整性"].notna().sum())
    total = len(merged)
    if matched == 0:
        raise RuntimeError("CheckM2 的 MAG 名称与 meta_plas_vf_card.tsv 的 Name 列完全无法匹配，未写入结果。")
    merged[["完整性", "污染率"]] = merged[["完整性", "污染率"]].fillna("-")

    if backup:
        backup_path = result_path.with_name("meta_plas_vf_card.before_checkm2_merge.tsv")
        if not backup_path.exists():
            shutil.copy2(result_path, backup_path)
    merged.to_csv(result_path, sep="\t", index=False)
    return matched, total


def process_task(task_dir: Path, args: argparse.Namespace) -> None:
    genome_dir = task_dir / "BASALT_out/meta_drep_out/binning_genomes"
    quality_report = task_dir / "bin_checkm2out/quality_report.tsv"
    sample = infer_sample_name(task_dir, args.sample) if not args.skip_mge else ""

    checkm2_state = "跳过"
    if args.with_checkm2:
        if not genome_dir.is_dir():
            raise RuntimeError(f"缺少 binning_genomes 目录：{genome_dir}")
        if not list(genome_dir.glob("*.fa")):
            raise RuntimeError(f"binning_genomes 中没有 .fa MAG 文件：{genome_dir}")
        checkm2_state = "重跑" if args.force or not quality_report.is_file() else "复用"
    mge_state = "跳过" if args.skip_mge else "重跑"
    print(f"[CheckM2:{checkm2_state} MGE:{mge_state}] {task_dir}")
    if args.dry_run:
        if sample:
            print(f"  样本前缀：{sample}")
        return
    if args.with_checkm2:
        if args.force or not quality_report.is_file() or quality_report.stat().st_size == 0:
            run_checkm2(task_dir, args, genome_dir, quality_report)
        matched, total = merge_quality_report(task_dir, quality_report, backup=not args.no_backup)
        print(f"  已回填完整性/污染率：{matched}/{total} 条 contig 记录匹配 CheckM2 MAG。")
    if not args.skip_mge:
        run_mge(task_dir, args, sample)
        print(f"  已重跑 MGE：{sample}.integrated_mge_summary.tsv / {sample}.mge_risk_summary.tsv")


def main() -> int:
    args = parse_args()
    if args.threads < 1:
        raise ValueError("--threads 必须大于 0")
    args.deploy_root = args.deploy_root.expanduser().resolve()
    args.checkm2_db = args.checkm2_db.expanduser().resolve()
    args.genomad_db = resolve_deploy_path(args.genomad_db, args.deploy_root)
    args.mobileog_db = resolve_deploy_path(args.mobileog_db, args.deploy_root)
    args.mobileog_meta = resolve_deploy_path(args.mobileog_meta, args.deploy_root)
    if args.with_checkm2:
        try:
            import pandas  # noqa: F401
        except ImportError as error:
            raise RuntimeError("该脚本需要 pandas；请在包含 pandas 的 Python/conda 环境中运行。") from error

    tasks = discover_tasks(args.tasks_dir)
    if not tasks:
        raise RuntimeError("未找到 meta_plas_vf_card.tsv；请确认 --tasks-dir 指向已完成任务的父目录。")
    needs_checkm2 = args.with_checkm2 and (
        args.force or any(
        not (task / "bin_checkm2out/quality_report.tsv").is_file() or
        (task / "bin_checkm2out/quality_report.tsv").stat().st_size == 0
        for task in tasks
        )
    )
    if needs_checkm2 and not args.dry_run and not args.checkm2_db.is_file():
        raise FileNotFoundError(f"CheckM2 数据库不存在：{args.checkm2_db}")
    if not args.skip_mge and not args.dry_run and not args.genomad_db.is_dir():
        raise FileNotFoundError(f"geNomad 数据库目录不存在：{args.genomad_db}")
    if not args.skip_mge and not args.dry_run and not args.mobileog_db.exists():
        raise FileNotFoundError(f"mobileOG DIAMOND 数据库不存在：{args.mobileog_db}")
    if not args.skip_mge and not args.dry_run and not args.mobileog_meta.is_file():
        raise FileNotFoundError(f"mobileOG 注释表不存在：{args.mobileog_meta}")
    print(f"发现 {len(tasks)} 个待补救任务。")
    failures: list[tuple[Path, Exception]] = []
    for task in tasks:
        try:
            process_task(task, args)
        except Exception as error:  # Continue so a batch can salvage the remaining tasks.
            failures.append((task, error))
            print(f"  失败：{error}", file=sys.stderr)
    print(f"完成：成功 {len(tasks) - len(failures)}，失败 {len(failures)}。")
    if failures:
        print("失败任务：", file=sys.stderr)
        for task, error in failures:
            print(f"- {task}: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
