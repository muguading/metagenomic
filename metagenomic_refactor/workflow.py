from __future__ import annotations

import os
import subprocess
import sys
import time
import traceback

from metagenomic_refactor.common import format_seconds, is_fastq_input
from metagenomic_refactor.context import get_runtime_context, update_runtime_context
from metagenomic_refactor.qc import QC_func
from metagenomic_refactor.report import combine_func
from metagenomic_refactor.strain_typing import level2Spe
from metagenomic_refactor.virus_analysis import (
    _is_hepatovirus,
    detect_influenza_type,
    prepare_influenza_reference_set,
    resolve_hepatovirus_reference,
)


CALLBACKS = {}


def register_workflow_callbacks(**kwargs):
    CALLBACKS.update(kwargs)


def ensure_sample_result_file():
    runtime = get_runtime_context()
    result_fp = f"{runtime.ofn}/sample_result.txt"
    if not os.path.isfile(result_fp):
        with open(result_fp, "w") as f:
            f.write("样本名称\tContig数量\tN50长度\t最长片段长度\t主基因组是否成环\t污染率\t完整性\n")


def print_progress(step, total_step, sam_num, all_num, sample, msg, runtime=None):
    text = f"task_step：{step}/{total_step}\t样本进度：{sam_num}/{all_num}\t样本：{sample}\t{msg}"
    if runtime is not None:
        text += f"\t运行时间:{runtime}"
    print(text)
    sys.stdout.flush()


def append_sample_result(pre):
    runtime = get_runtime_context()
    result_tsv = f"{pre}.assemble.result.tsv"
    if os.path.isfile(result_tsv):
        subprocess.run(f"sed -n '2p' {result_tsv} >> {runtime.ofn}/sample_result.txt", shell=True)


def is_assembly_output_ready(pre: str, method: str) -> bool:
    if str(method or "").strip() == "meta":
        meta_output = "tmp_combine.fa"
        return os.path.isfile(meta_output) and os.path.getsize(meta_output) > 0
    final_fasta = f"{pre}.final.fasta"
    return os.path.isfile(final_fasta) and os.path.getsize(final_fasta) > 0


def get_flow_list():
    runtime = get_runtime_context()
    raw_flows = [x.strip() for x in runtime.runflow.split(",") if x.strip()]
    expanded: list[str] = []
    for flow in raw_flows:
        if flow == "mlst与血清型":
            for item in ("mlst检验", "血清型检验"):
                if item not in expanded:
                    expanded.append(item)
            continue
        if flow not in expanded:
            expanded.append(flow)
    return expanded


def get_read_mode(pre):
    ont = os.path.isfile(f"{pre}.final.fastq")
    r1 = os.path.isfile(f"{pre}.R1.fastq.gz")
    r2 = os.path.isfile(f"{pre}.R2.fastq.gz")
    if ont:
        if r1 and r2:
            return "ont_pe"
        if r1:
            return "ont_se"
        return "ont"
    if r1 and r2:
        return "pe"
    if r1:
        return "se"
    return "none"


def infer_species(pre, llid):
    try:
        if llid != "nolevel":
            return level2Spe(pre, llid.split(",")[0])
    except Exception:
        pass
    return 0


def _sample_skip_flag_path(pre: str) -> str:
    return f"{pre}.skip_remaining.txt"


def _clear_sample_skip_flag(pre: str) -> None:
    path = _sample_skip_flag_path(pre)
    if os.path.isfile(path):
        os.remove(path)


def _read_sample_skip_reason(pre: str) -> str:
    path = _sample_skip_flag_path(pre)
    if not os.path.isfile(path):
        return ""
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as handle:
            return handle.read().strip()
    except OSError:
        return ""


def _is_influenza_scope() -> bool:
    runtime = get_runtime_context()
    species = str(runtime.species or "").strip().lower()
    if runtime.analysis_target != "virus":
        return False
    if "parainfluenza" in species or "副流感" in species:
        return False
    return "influenza" in species or "流感" in species


def _is_sars_cov_2_scope() -> bool:
    runtime = get_runtime_context()
    species = str(runtime.species or "").strip().lower()
    if runtime.analysis_target != "virus":
        return False
    return (
        species == "sars-cov-2"
        or species == "新型冠状病毒"
        or species == "新冠病毒"
        or species == "新冠"
        or "sars-cov-2" in species
    )


def run_fake_pipeline(pre, sam_num, all_num):
    print_progress(0, 5, sam_num, all_num, pre, "开始假流程分析")
    subprocess.run("cp -r /data1/shanghai_pip/meta_genome/fakefile/Sample1/* ./", shell=True)
    subprocess.run(f"rename 's/Sample1/{pre}/g' *", shell=True)
    subprocess.run(f"rename 's/Sample1/{pre}/g' */*", shell=True)
    time.sleep(5)
    print_progress(1, 5, sam_num, all_num, pre, "数据质控结束")
    time.sleep(5)
    print_progress(2, 5, sam_num, all_num, pre, "序列鉴定结束")
    time.sleep(5)
    print_progress(3, 5, sam_num, all_num, pre, "基因组完成图组装结束")
    time.sleep(5)
    print_progress(4, 5, sam_num, all_num, pre, "基因注释结束")
    time.sleep(5)
    print_progress(5, 5, sam_num, all_num, pre, "数据合并结束")


def run_fastq_qc_and_assembly(infile, fastq1, fastq2, threads, llid, mm, pre, ml, mq, asmtype, pts, pst, rnalib, ref, gtf, flowlist, sam_num, all_num, start_time):
    runtime = get_runtime_context()
    _clear_sample_skip_flag(pre)
    assembly_selected = "基因组组装" in runtime.runflow
    preidentified_species = False
    if _is_influenza_scope():
        if not os.path.isfile("QC_ok"):
            QC_func(infile, fastq1, fastq2, mq, ml, pre, rnalib, threads, runtime.method)
        if not os.path.isfile("kk2_ok"):
            run_species_identification(pre, threads)
        if "物种鉴定" in flowlist:
            flowlist.remove("物种鉴定")
        preidentified_species = True
        if asmtype in {"shortref", "longref"}:
            flu_result = prepare_influenza_reference_set(
                pre,
                runtime.species,
                runtime.ref or ref,
                single_fastq=infile,
                fq1=fastq1,
                fq2=fastq2,
                long_type=runtime.long_type,
                threads=threads,
            )
            detected_type = str(flu_result.get("influenza_type") or "Other").strip()
            if detected_type in {"Influenza A virus", "Influenza B virus"}:
                update_runtime_context(species=detected_type, ref=str(flu_result.get("reference_path") or "").strip())
                ref = str(flu_result.get("reference_path") or "").strip()
            else:
                detected_type = "Other"
        else:
            detected_type = detect_influenza_type(pre, runtime.species)
            if detected_type != "-":
                update_runtime_context(species=detected_type)
        if detected_type == "Other":
            t1 = time.time()
            runtime_text = format_seconds(t1 - start_time)
            print_progress(1, 1, sam_num, all_num, pre, "物种鉴定结束，未判定为甲/乙流感，终止后续组装与分型", runtime_text)
            return {"flowstep": 1, "flownum": 1, "timepoint": t1, "species_only_complete": True}
        if detected_type in {"Influenza C virus", "Influenza D virus"}:
            if "virus_typing" in CALLBACKS:
                CALLBACKS["virus_typing"](pre, detected_type)
            t1 = time.time()
            runtime_text = format_seconds(t1 - start_time)
            print_progress(1, 1, sam_num, all_num, pre, f"物种鉴定结束，判定为{detected_type}，跳过后续组装与分型", runtime_text)
            return {"flowstep": 1, "flownum": 1, "timepoint": t1, "species_only_complete": True}

    if not is_assembly_output_ready(pre, runtime.method):
        flownum = len(flowlist) + 1
        flowstep = 1
        if "基因组组装" in flowlist:
            flowlist.remove("基因组组装")

        if not os.path.isfile("QC_ok"):
            QC_func(infile, fastq1, fastq2, mq, ml, pre, rnalib, threads, runtime.method)

        if "基因组组装" in runtime.runflow:
            if os.path.isfile(f"{pre}.R2.fastq.gz"):
                CALLBACKS["asb_func"](f"{pre}.final.fastq", f"{pre}.R1.fastq.gz", f"{pre}.R2.fastq.gz", threads, pre, llid, pts, pst, runtime.method, asmtype, ref, gtf)
            elif os.path.isfile(f"{pre}.R1.fastq.gz"):
                CALLBACKS["asb_func"](f"{pre}.final.fastq", f"{pre}.R1.fastq.gz", 0, threads, pre, llid, pts, pst, runtime.method, asmtype, ref, gtf)
            else:
                CALLBACKS["asb_func"](f"{pre}.final.fastq", 0, 0, threads, pre, llid, pts, pst, runtime.method, asmtype, ref, gtf)

            skip_reason = _read_sample_skip_reason(pre)
            if skip_reason:
                t1 = time.time()
                runtime_text = format_seconds(t1 - start_time)
                print_progress(1, 1, sam_num, all_num, pre, f"基因分型支持 reads 过少，跳过当前样本：{skip_reason}", runtime_text)
                return {"flowstep": 1, "flownum": 1, "timepoint": t1, "species_only_complete": True}

            if runtime.method != "meta":
                CALLBACKS["Annotate_func"](pre, threads)

        t1 = time.time()
        runtime_text = format_seconds(t1 - start_time)
        finish_message = "数据质控&物种鉴定&组装结束" if preidentified_species else "数据质控&组装结束"
        print_progress(flowstep, flownum, sam_num, all_num, pre, finish_message, runtime_text)
        return {"flowstep": flowstep, "flownum": flownum, "timepoint": t1, "species_only_complete": False}

    if "基因组组装" in flowlist:
        flowlist.remove("基因组组装")
    flownum = len(flowlist) + (1 if assembly_selected else 0) + 1
    flowstep = 1 if assembly_selected else 0
    t1 = time.time()
    if assembly_selected:
        runtime_text = format_seconds(t1 - start_time)
        print_progress(flowstep, flownum, sam_num, all_num, pre, "组装结果已完成，进行后面的分析", runtime_text)
    return {"flowstep": flowstep, "flownum": flownum, "timepoint": t1, "species_only_complete": False}


def run_species_identification(pre, threads):
    mode = get_read_mode(pre)
    if mode == "ont":
        print("三代鉴定")
        sys.stdout.flush()
        CALLBACKS["kk2"](f"{pre}.final.fastq", 0, 0, threads, pre)
    elif mode == "ont_pe":
        print("三代鉴定+二代双端鉴定")
        sys.stdout.flush()
        CALLBACKS["kk2"](f"{pre}.final.fastq", f"{pre}.R1.fastq.gz", f"{pre}.R2.fastq.gz", threads, pre)
    elif mode == "ont_se":
        print("三代鉴定+二代单端鉴定")
        sys.stdout.flush()
        CALLBACKS["kk2"](f"{pre}.final.fastq", f"{pre}.R1.fastq.gz", 0, threads, pre)
    elif mode == "pe":
        print("二代双端鉴定")
        sys.stdout.flush()
        CALLBACKS["kk2"](0, f"{pre}.R1.fastq.gz", f"{pre}.R2.fastq.gz", threads, pre)
    elif mode == "se":
        print("二代单端鉴定")
        sys.stdout.flush()
        CALLBACKS["kk2"](0, f"{pre}.R1.fastq.gz", 0, threads, pre)
    else:
        print(f"{pre} 未检测到可用于物种鉴定的 reads")
        sys.stdout.flush()


def run_fastq_flow(flow, pre, threads, llid):
    runtime = get_runtime_context()
    if flow == "物种鉴定":
        if not os.path.isfile("kk2_ok"):
            run_species_identification(pre, threads)
    elif flow == "结构变异检测":
        if not os.path.isfile("SV_ok"):
            CALLBACKS["DrugFinder"](pre, threads)
    elif flow == "功能注释":
        CALLBACKS["AnnoFun"](pre, threads)
    elif flow == "元件预测":
        CALLBACKS["AnnoEle"](pre, threads)
    elif flow == "病毒组装" and runtime.method == "meta":
        CALLBACKS["meta_viral_assembly"](pre, threads)
    elif flow == "耐药与毒力":
        CALLBACKS["VFDR"](pre, threads)
    elif flow == "mlst检验":
        species = infer_species(pre, llid)
        CALLBACKS["mlst_only"](pre, species)
    elif flow == "血清型检验":
        species = infer_species(pre, llid)
        CALLBACKS["serotype_only"](pre, species)
    elif flow == "分型鉴定" and runtime.analysis_target == "virus":
        CALLBACKS["virus_typing"](pre, runtime.species)


def run_fastq_pipeline(infile, fastq1, fastq2, threads, llid, mm, pre, ml, mq, asmtype, all_num, sam_num, pts, pst, rnalib, ref, gtf, flowlist, start_time):
    init_result = run_fastq_qc_and_assembly(infile, fastq1, fastq2, threads, llid, mm, pre, ml, mq, asmtype, pts, pst, rnalib, ref, gtf, flowlist, sam_num, all_num, start_time)
    flowstep = int(init_result["flowstep"])
    flownum = int(init_result["flownum"])
    if init_result.get("species_only_complete"):
        return
    for flow in flowlist:
        flowstep += 1
        step_start = time.time()
        try:
            run_fastq_flow(flow, pre, threads, llid)
        except Exception as e:
            print(f"{pre} {flow} 失败: {e}")
            sys.stdout.flush()
        runtime_text = format_seconds(time.time() - step_start)
        print_progress(flowstep, flownum, sam_num, all_num, pre, f"{flow}已结束", runtime_text)

    flowstep += 1
    try:
        combine_func(pre)
    except Exception as e:
        print(f"{pre} 数据合并失败: {e}")
        traceback.print_exc()
        sys.stdout.flush()
    total_runtime = format_seconds(time.time() - start_time)
    if "基因组组装" in get_runtime_context().runflow:
        append_sample_result(pre)
    print_progress(flowstep, flownum, sam_num, all_num, pre, "数据合并结束", total_runtime)


def run_fasta_init(pre, threads, flowlist, sam_num, all_num, start_time):
    if _is_sars_cov_2_scope():
        flowlist[:] = [item for item in flowlist if item not in {"基因组组装", "物种鉴定"}]
        rt_flag = 0
        flowstep = 0
        flownum = len(flowlist) + 1
        return rt_flag, flowstep, flownum
    runtime = get_runtime_context()
    prokka_dir = f"{runtime.ofn}/fastq_analysis/{pre}/{pre}_prokka"
    if not os.path.isdir(prokka_dir):
        rt_flag = 1
        flownum = len(flowlist) + 2
        flowstep = 1
        CALLBACKS["Annotate_func"](pre, threads)
        t1 = time.time()
        runtime_text = format_seconds(t1 - start_time)
        print_progress(flowstep, flownum, sam_num, all_num, pre, "基因预测已完成", runtime_text)
        return rt_flag, flowstep, flownum
    rt_flag = 0
    flowstep = 0
    flownum = len(flowlist) + 1
    return rt_flag, flowstep, flownum


def run_fasta_flow(flow, pre, threads, llid):
    runtime = get_runtime_context()
    if flow == "物种鉴定" and runtime.analysis_target == "virus":
        if _is_sars_cov_2_scope():
            return
        if _is_hepatovirus(runtime.species):
            final_fasta = f"{pre}.final.fasta"
            selection = resolve_hepatovirus_reference(
                pre,
                species=runtime.species,
                requested_ref=str(runtime.ref or "").strip(),
                query_fasta=final_fasta,
                threads=threads,
            )
            selected_ref = str(selection.get("reference_path") or "").strip()
            selected_gtf = str(selection.get("gff_path") or "").strip() or "nogtf"
            selected_species = str(selection.get("species_label") or runtime.species).strip() or runtime.species
            if selected_ref:
                update_runtime_context(ref=selected_ref, gtf=selected_gtf, species=selected_species)
            return
        CALLBACKS["virus_nextclade_identify"](pre, runtime.species)
    elif flow == "功能注释":
        CALLBACKS["AnnoFun"](pre, threads)
    elif flow == "元件预测":
        CALLBACKS["AnnoEle"](pre, threads)
    elif flow == "病毒组装" and runtime.method == "meta":
        CALLBACKS["meta_viral_assembly"](pre, threads)
    elif flow == "耐药与毒力":
        CALLBACKS["VFDR"](pre, threads, "fasta")
    elif flow == "mlst检验":
        species = infer_species(pre, llid)
        CALLBACKS["mlst_only"](pre, species)
    elif flow == "血清型检验":
        species = infer_species(pre, llid)
        CALLBACKS["serotype_only"](pre, species)
    elif flow == "分型鉴定" and runtime.analysis_target == "virus":
        CALLBACKS["virus_typing"](pre, runtime.species)


def run_fasta_pipeline(infile, threads, llid, mm, pre, all_num, sam_num, flowlist, start_time):
    rt_flag, flowstep, flownum = run_fasta_init(pre, threads, flowlist, sam_num, all_num, start_time)
    for flow in flowlist:
        flowstep += 1
        step_start = time.time()
        try:
            run_fasta_flow(flow, pre, threads, llid)
        except Exception as e:
            print(f"{pre} {flow} 失败: {e}")
            sys.stdout.flush()
        runtime_text = format_seconds(time.time() - step_start)
        print_progress(flowstep, flownum, sam_num, all_num, pre, f"{flow}已结束", runtime_text)

    flowstep += 1
    try:
        combine_func(pre)
    except Exception as e:
        print(f"{pre} 数据合并失败: {e}")
        traceback.print_exc()
        sys.stdout.flush()
    total_runtime = format_seconds(time.time() - start_time)
    if rt_flag:
        append_sample_result(pre)
    print_progress(flowstep, flownum, sam_num, all_num, pre, "数据合并结束", total_runtime)


def main_process(infile, fastq1, fastq2, threads, llid, mm, Pre, ml, mq, asmtype, Allnum, Samnum, pts, pst, rnalib, iffake=0, ref="noref", gtf="nogtf"):
    start_time = time.time()
    runtime = get_runtime_context()
    update_runtime_context(
        species=runtime.base_species or runtime.species,
        ref=runtime.base_ref or ref,
        gtf=runtime.base_gtf or gtf,
    )
    if iffake:
        run_fake_pipeline(Pre, Samnum, Allnum)
        return
    ensure_sample_result_file()
    flowlist = get_flow_list()
    total_for_start = len(flowlist) + 1
    print_progress(0, total_for_start, Samnum, Allnum, Pre, "开始分析")
    if is_fastq_input(infile, fastq1, fastq2):
        run_fastq_pipeline(infile, fastq1, fastq2, threads, llid, mm, Pre, ml, mq, asmtype, Allnum, Samnum, pts, pst, rnalib, ref, gtf, flowlist, start_time)
    else:
        run_fasta_pipeline(infile, threads, llid, mm, Pre, Allnum, Samnum, flowlist, start_time)
