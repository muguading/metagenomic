#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import multiprocessing
import shutil
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


DEFAULT_RUNFLOW = "基因组组装,物种鉴定,结构变异检测,功能注释,元件预测,耐药与毒力,mlst检验,血清型检验"
DEFAULT_SPECIES = "Neisseria meningitidis"
DEFAULT_METHOD = "spades,freebayes"
DEFAULT_ASM_TYPE = "shortasm"
DEFAULT_REF = "Nmen"
DEFAULT_GTF = "nogtf"
DEFAULT_GENOME_LEN = "2.2m"
DEFAULT_POLISH_TIMES = "1"
DEFAULT_POLISH_SOFT = "medaka"
DEFAULT_BARCODEKIT = "none"
DEFAULT_MIN_LEN = "500"
DEFAULT_MIN_Q = "10"
DEFAULT_RMHOST = "norm"
DEFAULT_ABUN = "1"
DEFAULT_RNA = "0"
DEFAULT_IFANNO = "Anno"
DEFAULT_LONG_TYPE = "Nanopore"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="脑膜炎奈瑟菌专用全流程脚本：从原始数据质控到最终报告输出，并补充耐药位点判读。"
    )
    parser.add_argument("-i", "--input", required=True, help="输入文件路径")
    parser.add_argument("-o", "--output", required=True, help="输出目录")
    parser.add_argument(
        "--inputtype",
        "-tp",
        default="fastq",
        help="输入文件类型，常见为 fastq / barcode_fastq / fasta / @",
    )
    parser.add_argument("--thread", "-t", type=int, default=10, help="线程数")
    parser.add_argument(
        "--barcodekit",
        "-bk",
        default=DEFAULT_BARCODEKIT,
        help="只有 barcode 数据时才需要填写，其他情况保持默认。",
    )
    parser.add_argument(
        "--skip-serotype",
        action="store_true",
        help="跳过血清型检验；默认会跑脑膜炎奈瑟菌血清型。",
    )
    parser.add_argument(
        "--skip-amr-postprocess",
        action="store_true",
        help="跳过脑膜炎奈瑟菌专属耐药位点补充判读。",
    )
    return parser.parse_args()


def _normalize_delivery_args(args: argparse.Namespace) -> argparse.Namespace:
    args.analysis_target = "bacteria"
    args.minlongfilt = DEFAULT_MIN_LEN
    args.Qfilt = DEFAULT_MIN_Q
    args.fake_pip = 0
    args.method = DEFAULT_METHOD
    args.long_type = DEFAULT_LONG_TYPE
    args.ref = DEFAULT_REF
    args.gtf = DEFAULT_GTF
    args.genome_len = DEFAULT_GENOME_LEN
    args.asm_type = DEFAULT_ASM_TYPE
    args.polish_times = DEFAULT_POLISH_TIMES
    args.polish_soft = DEFAULT_POLISH_SOFT
    args.species = DEFAULT_SPECIES
    args.ifanno = DEFAULT_IFANNO
    args.rmhost = DEFAULT_RMHOST
    args.abun = DEFAULT_ABUN
    args.rna = DEFAULT_RNA
    if args.skip_serotype:
        args.runflow = "基因组组装,物种鉴定,结构变异检测,功能注释,元件预测,耐药与毒力,mlst检验"
    else:
        args.runflow = DEFAULT_RUNFLOW
    return args


def _looks_like_neisseria_meningitidis(value: object) -> bool:
    text = str(value or "").strip().lower()
    if not text:
        return False
    return (
        "neisseria meningitidis" in text
        or "meningitidis" in text
        or "脑膜炎奈瑟" in text
        or "流脑" in text
    )


def _read_first_tsv_row(path: Path) -> dict[str, str]:
    if not path.is_file():
        return {}
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            row = next(reader, None)
            return {str(k): str(v or "") for k, v in (row or {}).items()}
    except OSError:
        return {}


def _iter_report_dirs(output_dir: Path):
    seen: set[Path] = set()
    for mlst_path in sorted(output_dir.rglob("*.mlst_Stat.txt")):
        report_dir = mlst_path.parent.resolve()
        if report_dir in seen:
            continue
        seen.add(report_dir)
        sample_name = mlst_path.name.removesuffix(".mlst_Stat.txt")
        yield report_dir, sample_name


def _find_gbk(report_dir: Path, sample_name: str) -> Path | None:
    for candidate in (
        report_dir / f"{sample_name}_prokka" / f"{sample_name}.gbk",
        report_dir / f"{sample_name}_prokka" / "main.gbk",
        report_dir / f"{sample_name}.gbk",
    ):
        if candidate.is_file():
            return candidate
    return None


def _copy_amr_call_into_result_dir(report_dir: Path, sample_name: str, output_path: Path) -> None:
    result_dir = report_dir / f"{sample_name}_genome_complete_result" / "5.Fun_Element"
    if not result_dir.exists():
        return
    result_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(output_path, result_dir / output_path.name)


def _postprocess_neisseria_amr(project_root: Path, output_dir: Path) -> list[Path]:
    script_path = project_root / "scripts" / "check_neisseria_meningitidis_amr_sites.py"
    site_table = project_root / "database" / "NM_mutate" / "neisseria_meningitidis_snp_amr_associations_literature_updated.csv"
    if not script_path.is_file():
        raise FileNotFoundError(f"脑膜炎奈瑟菌耐药位点脚本不存在: {script_path}")
    if not site_table.is_file():
        raise FileNotFoundError(f"脑膜炎奈瑟菌耐药位点数据库不存在: {site_table}")

    generated: list[Path] = []
    for report_dir, sample_name in _iter_report_dirs(output_dir):
        checkm_row = _read_first_tsv_row(report_dir / f"{sample_name}.checkm.tsv")
        if not any(
            [
                _looks_like_neisseria_meningitidis(checkm_row.get("物种名称")),
                _looks_like_neisseria_meningitidis(checkm_row.get("species_name")),
                _looks_like_neisseria_meningitidis(checkm_row.get("mlst 物种名称")),
                _looks_like_neisseria_meningitidis(checkm_row.get("mlst_species_name")),
            ]
        ):
            continue
        gbk_path = _find_gbk(report_dir, sample_name)
        if gbk_path is None:
            print(f"[WARN] {sample_name} 未找到 GBK，跳过脑膜炎奈瑟耐药位点补充判读。", file=sys.stderr)
            continue
        output_path = report_dir / f"{sample_name}.neisseria_amr_calls.csv"
        command = [
            sys.executable,
            str(script_path),
            str(gbk_path),
            "--site-table",
            str(site_table),
            "--output",
            str(output_path),
        ]
        result = subprocess.run(
            command,
            cwd=str(project_root),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"{sample_name} 脑膜炎奈瑟耐药位点补充判读失败。\n"
                f"命令: {' '.join(command)}\n"
                f"输出:\n{result.stdout}"
            )
        generated.append(output_path)
        _copy_amr_call_into_result_dir(report_dir, sample_name, output_path)
    return generated


def _build_runtime(argv: argparse.Namespace):
    from metagenomic_refactor.annotation import AnnoFun, DrugFinder, VFDR
    from metagenomic_refactor.assembly import Annotate_func, asb_func, run_meta_viral_assembly
    from metagenomic_refactor.config import build_pipeline_config
    from metagenomic_refactor.context import RuntimeContext, set_runtime_context
    from metagenomic_refactor.genomad_mge import AnnoEle
    from metagenomic_refactor.strain_typing import mlst_only, mlst_serotype, serotype_only
    from metagenomic_refactor.taxonomy import kk2
    from metagenomic_refactor.virus_analysis import nextclade_identify, virus_typing
    from metagenomic_refactor.workflow import register_workflow_callbacks

    pipeline_config = build_pipeline_config(argv, multiprocessing.cpu_count())
    set_runtime_context(
        RuntimeContext(
            ofn=pipeline_config.ofn,
            runflow=pipeline_config.runflow,
            method=pipeline_config.method,
            rmhost=pipeline_config.rmhost,
            tspeabun=pipeline_config.tspeabun,
            nt=pipeline_config.nt,
            krdb=pipeline_config.krdb,
            wkdir=pipeline_config.ofn,
            long_type=pipeline_config.long_type,
            analysis_target=pipeline_config.analysis_target,
            species=pipeline_config.species,
            genome_len=pipeline_config.genome_len,
            ref=pipeline_config.ref,
            gtf=pipeline_config.gtf,
            vfmeta=pipeline_config.vfmeta,
            resources=pipeline_config.resources,
        )
    )
    register_workflow_callbacks(
        asb_func=asb_func,
        Annotate_func=Annotate_func,
        kk2=kk2,
        DrugFinder=DrugFinder,
        AnnoFun=AnnoFun,
        AnnoEle=AnnoEle,
        VFDR=VFDR,
        mlst_only=mlst_only,
        serotype_only=serotype_only,
        mlst_serotype=mlst_serotype,
        virus_nextclade_identify=nextclade_identify,
        virus_typing=virus_typing,
        meta_viral_assembly=run_meta_viral_assembly,
    )
    return pipeline_config


def _run_standalone_pipeline(pipeline_config) -> None:
    from metagenomic_refactor.runner import RunnerConfig, run_pipeline_entry
    from metagenomic_refactor.workflow import main_process

    output_dir = Path(pipeline_config.ofn)
    output_dir.mkdir(parents=True, exist_ok=True)
    run_pipeline_entry(
        RunnerConfig(
            raw_input=pipeline_config.inf,
            inf=pipeline_config.inf,
            intype=pipeline_config.intype,
            ofn=pipeline_config.ofn,
            barkit=pipeline_config.barkit,
            tmpfake=pipeline_config.tmpfake,
            fastq1=pipeline_config.fastq1,
            fastq2=pipeline_config.fastq2,
            nt=pipeline_config.nt,
            llid=pipeline_config.species if pipeline_config.species != "False" else "nolevel",
            mmethod=pipeline_config.mmethod,
            minl=pipeline_config.minl,
            minQ=pipeline_config.minQ,
            asm_type=pipeline_config.asm_type,
            ptimes=pipeline_config.ptimes,
            psoft=pipeline_config.psoft,
            rnalib=pipeline_config.rnalib,
            ref=pipeline_config.ref,
            gtf=pipeline_config.gtf,
        ),
        main_process=main_process,
    )


def _collect_reports(output_dir: Path) -> list[Path]:
    return sorted(output_dir.rglob("*_bacgenome.html"))


def main() -> int:
    args = _normalize_delivery_args(parse_args())

    print("==> 启动脑膜炎奈瑟菌独立分析脚本")
    print(f"项目根目录: {PROJECT_ROOT}")
    print(f"物种: {args.species}")
    print(f"输入类型: {args.inputtype}")
    print(f"运行模块: {args.runflow}")
    sys.stdout.flush()

    pipeline_config = _build_runtime(args)
    _run_standalone_pipeline(pipeline_config)

    amr_outputs: list[Path] = []
    output_dir = Path(pipeline_config.ofn).resolve()
    if not args.skip_amr_postprocess:
        print("==> 开始补充脑膜炎奈瑟菌耐药位点判读")
        sys.stdout.flush()
        amr_outputs = _postprocess_neisseria_amr(PROJECT_ROOT, output_dir)

    report_paths = _collect_reports(output_dir)
    print("==> 脑膜炎奈瑟菌分析完成")
    if report_paths:
        print("报告文件:")
        for path in report_paths:
            print(f"  - {path}")
    if amr_outputs:
        print("脑膜炎奈瑟菌耐药位点结果:")
        for path in amr_outputs:
            print(f"  - {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
