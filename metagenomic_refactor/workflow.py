from __future__ import annotations

import os
import subprocess
import sys
import time

from metagenomic_refactor.common import format_seconds, is_fastq_input
from metagenomic_refactor.context import get_runtime_context
from metagenomic_refactor.qc import QC_func
from metagenomic_refactor.report import combine_func


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


def get_flow_list():
    runtime = get_runtime_context()
    return [x.strip() for x in runtime.runflow.split(",") if x.strip()]


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
            return CALLBACKS["level2Spe"](pre, llid.split(",")[0])
    except Exception:
        pass
    return 0


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
    if not os.path.isfile(f"{pre}.final.fasta"):
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

            if runtime.method != "meta":
                CALLBACKS["Annotate_func"](pre, threads)

        t1 = time.time()
        runtime_text = format_seconds(t1 - start_time)
        print_progress(flowstep, flownum, sam_num, all_num, pre, "数据质控&组装结束", runtime_text)
        return flowstep, flownum, t1

    flownum = len(flowlist)
    flowstep = 0
    t1 = time.time()
    return flowstep, flownum, t1


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
    elif flow == "耐药与毒力":
        CALLBACKS["VFDR"](pre, threads)
    elif flow == "mlst与血清型":
        species = infer_species(pre, llid)
        CALLBACKS["mlst_serotype"](pre, species)


def run_fastq_pipeline(infile, fastq1, fastq2, threads, llid, mm, pre, ml, mq, asmtype, all_num, sam_num, pts, pst, rnalib, ref, gtf, flowlist, start_time):
    flowstep, flownum, _ = run_fastq_qc_and_assembly(infile, fastq1, fastq2, threads, llid, mm, pre, ml, mq, asmtype, pts, pst, rnalib, ref, gtf, flowlist, sam_num, all_num, start_time)
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
    combine_func(pre)
    total_runtime = format_seconds(time.time() - start_time)
    if "基因组组装" in get_runtime_context().runflow:
        append_sample_result(pre)
    print_progress(flowstep, flownum, sam_num, all_num, pre, "数据合并结束", total_runtime)


def run_fasta_init(pre, threads, flowlist, sam_num, all_num, start_time):
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
    if flow == "功能注释":
        CALLBACKS["AnnoFun"](pre, threads)
    elif flow == "元件预测":
        CALLBACKS["AnnoEle"](pre, threads)
    elif flow == "耐药与毒力":
        CALLBACKS["VFDR"](pre, threads, "fasta")
    elif flow == "mlst与血清型":
        species = infer_species(pre, llid)
        CALLBACKS["mlst_serotype"](pre, species)


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
    combine_func(pre)
    total_runtime = format_seconds(time.time() - start_time)
    if rt_flag:
        append_sample_result(pre)
    print_progress(flowstep, flownum, sam_num, all_num, pre, "数据合并结束", total_runtime)


def main_process(infile, fastq1, fastq2, threads, llid, mm, Pre, ml, mq, asmtype, Allnum, Samnum, pts, pst, rnalib, iffake=0, ref="noref", gtf="nogtf"):
    start_time = time.time()
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
