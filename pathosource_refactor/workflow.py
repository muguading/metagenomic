from __future__ import annotations

import os
import subprocess
import sys
import time

from pathosource_refactor.common import format_seconds
from pathosource_refactor.context import get_runtime_context, update_runtime_context
from pathosource_refactor.phylogeny import Pair_dis, build_tree, consensus_fun, snippy
from pathosource_refactor.typing import cgmlst, mlst


BUILTIN_REFS = [
    "salmonella",
    "E_coli",
    "Shigella",
    "Parahemolyticus",
    "cholerae",
    "Y_enterocolitica",
    "Campylobacter",
    "HPinf",
    "Brucella",
    "Lmono",
    "Kpne",
    "Suare",
    "Bcere",
    "Nmen",
]


def _resolve_reference(ref: str) -> str:
    runtime = get_runtime_context()
    if os.path.isfile(ref):
        return os.path.abspath(ref)
    if ref in BUILTIN_REFS:
        return runtime.speciesdb.loc[runtime.speciesdb["Species"] == ref, "reference"].tolist()[0]
    tref = f"{runtime.resources.ref_fadb_dir}/{ref}_genomic.fna.gz"
    subprocess.run(f"seqkit seq {tref} > tmp_ref.fa", shell=True)
    return f"{os.getcwd()}/tmp_ref.fa"


def _prepare_workspace() -> None:
    runtime = get_runtime_context()
    sample_root = f"{runtime.run_path}/fastq_analysis/sample"
    subprocess.run(f"rm -rf {sample_root}", shell=True)
    os.makedirs(sample_root, exist_ok=True)
    with open(f"{runtime.run_path}/sample_result.txt", "w") as handle:
        handle.write("sample")
    os.chdir(sample_root)
    update_runtime_context(run_path=runtime.run_path)


def _print_progress(step: int, total_step: int, msg: str, runtime_text: str | None = None) -> None:
    text = f"task_step：{step}/{total_step}\t样本进度：1/1\t样本：sample\t {msg}"
    if runtime_text is not None:
        text += f"\t运行时间:{runtime_text}"
    print(text)
    sys.stdout.flush()


def main_process(inputin, ref, inputout=0) -> None:
    runtime = get_runtime_context()
    _prepare_workspace()
    reference = _resolve_reference(ref)
    subprocess.run(f"seqkit seq {reference} > ref.fa", shell=True)
    subprocess.run(f"sed -i 's/\t/\t/g' {runtime.input_path}", shell=True)
    _print_progress(0, 4, "任务开始")
    start_time = time.time()

    snippy(inputout, inputin, "ref.fa")
    use_gubbins = 1 if str(runtime.gubbins or "yes").strip().lower() in {"yes", "y", "1", "true"} else 0
    consensus_fun("Samplelist.txt", runtime.msamethod, "ref.fa", use_gubbins, "test1")
    step_time = time.time()
    _print_progress(1, 4, "生成一致性序列已结束", format_seconds(step_time - start_time))

    build_tree("clean.core.aln", runtime.treemethod, subs=runtime.mltype, threads=runtime.threads, bs=runtime.bootstrap)
    step2_time = time.time()
    _print_progress(2, 4, "构建进化树已结束", format_seconds(step2_time - step_time))

    Pair_dis()
    step3_time = time.time()
    _print_progress(3, 4, "聚类计算已结束", format_seconds(step3_time - step2_time))

    mlst()
    if runtime.cgmlstana != "no":
        cgmlst("sample", runtime.species, runtime.cgmlstversion)
    step4_time = time.time()
    _print_progress(4, 4, "mlst&cgmlst分析已结束", format_seconds(step4_time - step3_time))
