from __future__ import annotations

import os
import subprocess
import sys
from dataclasses import dataclass

import pandas as pd

from metagenomic_refactor.common import basecaller, check_input, is_fasta, is_fastq
from metagenomic_refactor.context import update_runtime_context


@dataclass
class RunnerConfig:
    raw_input: str
    inf: str
    intype: str
    ofn: str
    barkit: str
    tmpfake: int
    fastq1: int
    fastq2: int
    nt: int
    llid: str
    mmethod: str
    minl: str
    minQ: str
    asm_type: str
    ptimes: str
    psoft: str
    rnalib: str
    ref: str
    gtf: str


def _set_wkdir(path: str) -> None:
    update_runtime_context(wkdir=path)


def _init_analysis_root(ofn: str) -> str:
    os.makedirs("fastq_analysis", exist_ok=True)
    with open("Sample_result.txt", "w") as f:
        f.write("样本名称\t序列数量\t碱基数量\t运行时间\n")
    root = os.path.join(ofn, "fastq_analysis")
    sample_list = os.path.join(root, "Samplelist.txt")
    if not os.path.isfile(sample_list):
        open(sample_list, "w").close()
    return root


def _record_sample(ofn: str, sample: str) -> None:
    with open(f"{ofn}/fastq_analysis/Samplelist.txt", "a") as f:
        f.write(f"{sample}\n")


def _ensure_sample_dir(root: str, sample: str) -> str:
    sample_dir = os.path.join(root, sample)
    os.makedirs(sample_dir, exist_ok=True)
    return sample_dir


def _run_sample(main_process, cfg: RunnerConfig, infile, fastq1, fastq2, sample, anum, snum, llid=None):
    target_llid = cfg.llid if llid is None else llid
    main_process(
        infile,
        fastq1,
        fastq2,
        cfg.nt,
        target_llid,
        cfg.mmethod,
        sample,
        cfg.minl,
        cfg.minQ,
        cfg.asm_type,
        anum,
        snum,
        cfg.ptimes,
        cfg.psoft,
        cfg.rnalib,
        cfg.tmpfake,
        cfg.ref,
        cfg.gtf,
    )


def _print_failure(sample: str, snum: int, anum: int) -> None:
    print(f"{sample}数据量不足,或者基因组组装失败")
    print(f"样本进度：{snum}/{anum}\t样本：{sample}\t数据分析中断")
    sys.stdout.flush()


def _is_fasta_suffix(path: str) -> bool:
    return path.endswith(("fa", "fas", "fna", "fasta", "fa.gz", "fas.gz", "fasta.gz", "fna.gz"))


def _run_non_list_mode(cfg: RunnerConfig, main_process) -> None:
    protype = check_input(cfg.inf, cfg.intype)
    print(protype)
    sys.stdout.flush()
    default_llid = cfg.llid or "nolevel"

    if protype == "fqdir":
        sample = "sample1"
        anum = 1
        snum = 1
        sample_dir = _ensure_sample_dir("fastq_analysis", sample)
        _record_sample(cfg.ofn, sample)
        os.chdir(sample_dir)
        _set_wkdir(os.getcwd())
        subprocess.run(f"cat {cfg.inf}/*.f*q* |seqkit rmdup -i |seqkit seq > {sample}.raw.fastq", shell=True)
        if os.path.getsize(f"{sample}.raw.fastq") == 0:
            print("文件夹内没有三代fastq格式文件")
            raise SystemExit()
        try:
            print("三代fastq文件夹分析模式")
            print(f"{sample}.raw.fastq", cfg.fastq1, cfg.fastq2, cfg.nt, cfg.mmethod, sample, cfg.minl, cfg.minQ, cfg.asm_type, anum, snum, cfg.ptimes, cfg.psoft)
            sys.stdout.flush()
            _run_sample(main_process, cfg, f"{sample}.raw.fastq", cfg.fastq1, cfg.fastq2, sample, anum, snum, default_llid)
        except Exception:
            print("fqdir分析失败")
        os.chdir(cfg.ofn)
        _set_wkdir(cfg.ofn)
    elif protype == "fqfile":
        sample = "sample1"
        anum = 1
        snum = 1
        sample_dir = _ensure_sample_dir("fastq_analysis", sample)
        _record_sample(cfg.ofn, sample)
        os.chdir("fastq_analysis")
        if cfg.tmpfake:
            subprocess.run("cp -r /data1/shanghai_pip/meta_genome/fake_result/* ./", shell=True)
        os.chdir(cfg.ofn)
        os.chdir(sample_dir)
        _set_wkdir(os.getcwd())
        subprocess.run(f"ln -s {cfg.inf} ./{sample}.raw.fastq", shell=True)
        if os.path.getsize(f"{sample}.raw.fastq") == 0:
            print("输入三代fastq格式文件有误，请核对")
            raise SystemExit()
        try:
            print("三代fastq文件分析模式")
            print(f"{sample}.raw.fastq", cfg.fastq1, cfg.fastq2, cfg.nt, cfg.mmethod, sample, cfg.minl, cfg.minQ, cfg.asm_type, anum, snum, cfg.ptimes, cfg.psoft)
            sys.stdout.flush()
            _run_sample(main_process, cfg, f"{sample}.raw.fastq", cfg.fastq1, cfg.fastq2, sample, anum, snum, default_llid)
        except Exception:
            print("fqfile分析失败")
        os.chdir(cfg.ofn)
        _set_wkdir(cfg.ofn)
    elif protype == "f5dir":
        os.chdir("fastq_analysis")
        _set_wkdir(os.getcwd())
        basecaller("fast5", cfg.intype, cfg.inf)
        if cfg.barkit != "none":
            subprocess.run(f"/home/dell/biosoft/ont-guppy/bin/guppy_barcoder -i basecaller_outputs/pass -s barout -x auto --barcode_kits {cfg.barkit}", shell=True)
            anum = len(os.listdir("barout"))
            snum = 1
            root = os.getcwd()
            for item in os.listdir("barout"):
                if item.startswith("barcode"):
                    sample_dir = _ensure_sample_dir(root, item)
                    os.chdir(sample_dir)
                    _set_wkdir(os.getcwd())
                    _record_sample(cfg.ofn, item)
                    subprocess.run(f"cat {cfg.ofn}/fastq_analysis/barout/{item}/*.f*q*|seqkit seq > {item}.raw.fastq", shell=True)
                    try:
                        _run_sample(main_process, cfg, f"{item}.raw.fastq", cfg.fastq1, cfg.fastq2, item, anum, snum, default_llid)
                    except Exception:
                        print("f5dir分析失败")
                        _print_failure(item, snum, anum)
                    snum += 1
                    os.chdir(root)
                    _set_wkdir(root)
            os.chdir(cfg.ofn)
            _set_wkdir(cfg.ofn)
        else:
            sample = "sample1"
            anum = 1
            snum = 1
            os.makedirs(sample, exist_ok=True)
            os.chdir(sample)
            _set_wkdir(os.getcwd())
            _record_sample(cfg.ofn, sample)
            subprocess.run(f"cat {cfg.ofn}/fastq_analysis/basecaller_outputs/pass/*.f*q*|seqkit seq > {sample}.raw.fastq", shell=True)
            try:
                _run_sample(main_process, cfg, f"{sample}.raw.fastq", cfg.fastq1, cfg.fastq2, sample, anum, snum, default_llid)
            except Exception:
                print("数据量不足")
            os.chdir(cfg.ofn)
            _set_wkdir(cfg.ofn)
            subprocess.run("cat */*/Sample_result.txt >> Sample_result.txt", shell=True)
    elif protype == "bardir":
        os.chdir("fastq_analysis")
        root = os.getcwd()
        _set_wkdir(root)
        anum = len([i for i in os.listdir(cfg.inf) if i.startswith("barcode")])
        snum = 1
        print("输入文件为barcode文件夹")
        for item in os.listdir(cfg.inf):
            if item.startswith("barcode"):
                try:
                    sample_dir = _ensure_sample_dir(root, item)
                    os.chdir(sample_dir)
                    _set_wkdir(os.getcwd())
                    _record_sample(cfg.ofn, item)
                    if len([fq for fq in os.listdir(f"{cfg.inf}/{item}") if is_fastq(f"{cfg.inf}/{item}/{fq}")]):
                        subprocess.run(f"cat {cfg.inf}/{item}/*.f*q*|seqkit rmdup -i |seqkit seq > {item}.raw.fastq", shell=True)
                    else:
                        subprocess.run(f"touch {item}.raw.fastq", shell=True)
                    _run_sample(main_process, cfg, f"{item}.raw.fastq", cfg.fastq1, cfg.fastq2, item, anum, snum, default_llid)
                except Exception:
                    _print_failure(item, snum, anum)
                snum += 1
                os.chdir(root)
                _set_wkdir(root)
        os.chdir(cfg.ofn)
        _set_wkdir(cfg.ofn)
    elif protype == "barfqdir":
        os.chdir("fastq_analysis")
        root = os.getcwd()
        _set_wkdir(root)
        os.makedirs("fastq_out", exist_ok=True)
        subprocess.run(f"/home/dell/biosoft/ont-guppy/bin/guppy_barcoder -i {cfg.inf} -s barout -x auto --barcode_kits {cfg.barkit}", shell=True)
        anum = len([i for i in os.listdir("barout") if i.startswith("barcode")])
        snum = 1
        for item in os.listdir("barout"):
            if item.startswith("barcode"):
                try:
                    sample_dir = _ensure_sample_dir(root, item)
                    os.chdir(sample_dir)
                    _set_wkdir(os.getcwd())
                    _record_sample(cfg.ofn, item)
                    if len([fq for fq in os.listdir(f"{cfg.ofn}/fastq_analysis/barout/{item}/") if is_fastq(f"{cfg.ofn}/fastq_analysis/barout/{item}/{fq}")]):
                        subprocess.run(f"cat {cfg.ofn}/fastq_analysis/barout/{item}/*.f*q*|seqkit seq > {item}.raw.fastq", shell=True)
                    else:
                        subprocess.run(f"touch {item}.raw.fastq", shell=True)
                    _run_sample(main_process, cfg, f"{item}.raw.fastq", cfg.fastq1, cfg.fastq2, item, anum, snum, default_llid)
                except Exception:
                    _print_failure(item, snum, anum)
                snum += 1
                os.chdir(root)
                _set_wkdir(root)
        os.chdir(cfg.ofn)
        _set_wkdir(cfg.ofn)
    elif protype == "barfqfile":
        os.chdir("fastq_analysis")
        root = os.getcwd()
        _set_wkdir(root)
        os.makedirs("fastq_out", exist_ok=True)
        subprocess.run(f"ln -s {cfg.inf} fastq_out", shell=True)
        if cfg.barkit != "none":
            subprocess.run(f"/home/dell/biosoft/ont-guppy/bin/guppy_barcoder -i fastq_out/ -s barout -x auto --barcode_kits {cfg.barkit}", shell=True)
            anum = len([i for i in os.listdir("barout") if i.startswith("barcode")])
            snum = 1
            for item in os.listdir("barout"):
                if item.startswith("barcode"):
                    try:
                        sample_dir = _ensure_sample_dir(root, item)
                        os.chdir(sample_dir)
                        _set_wkdir(os.getcwd())
                        _record_sample(cfg.ofn, item)
                        subprocess.run(f"cat {cfg.ofn}/fastq_analysis/barout/{item}/*.f*q*|seqkit rmdup -i |seqkit seq > {item}.raw.fastq", shell=True)
                        _run_sample(main_process, cfg, f"{item}.raw.fastq", cfg.fastq1, cfg.fastq2, item, anum, snum, default_llid)
                    except Exception:
                        _print_failure(item, snum, anum)
                    snum += 1
                    os.chdir(root)
                    _set_wkdir(root)
        os.chdir(cfg.ofn)
        _set_wkdir(cfg.ofn)
    elif protype == "fadir":
        os.chdir("fastq_analysis")
        root = os.getcwd()
        _set_wkdir(root)
        faendlist = ["fa", "fasta", "fna", "fas", "fa.gz", "fasta.gz", "fna.gz", "fas.gz"]
        anum = 0
        for tfile in os.listdir(cfg.inf):
            endssuf = [i for i in faendlist if tfile.endswith(i)]
            if endssuf:
                anum += 1
        snum = 1
        for tfile in os.listdir(cfg.inf):
            endssuf = [i for i in faendlist if tfile.endswith(i)]
            if endssuf:
                sample = tfile.replace(f".{endssuf[0]}", "")
                _record_sample(cfg.ofn, sample)
                sample_dir = _ensure_sample_dir(root, sample)
                os.chdir(sample_dir)
                _set_wkdir(os.getcwd())
                subprocess.run(f"seqkit seq -w 0  {cfg.inf}/{tfile} > {sample}.final.fasta", shell=True)
                try:
                    _run_sample(main_process, cfg, f"{cfg.inf}/{tfile}", cfg.fastq1, cfg.fastq2, sample, anum, snum, default_llid)
                except Exception:
                    _print_failure(sample, snum, anum)
                snum += 1
                os.chdir(cfg.ofn)
                _set_wkdir(cfg.ofn)
    elif protype == "fafile":
        os.chdir("fastq_analysis")
        root = os.getcwd()
        _set_wkdir(root)
        sample = "sample1"
        anum = 1
        snum = 1
        _record_sample(cfg.ofn, sample)
        sample_dir = _ensure_sample_dir(root, sample)
        os.chdir(sample_dir)
        _set_wkdir(os.getcwd())
        subprocess.run(f"seqkit seq -w 0  {cfg.inf} > {sample}.final.fasta", shell=True)
        _run_sample(main_process, cfg, cfg.inf, cfg.fastq1, cfg.fastq2, sample, anum, snum, default_llid)
        os.chdir(cfg.ofn)
        _set_wkdir(cfg.ofn)
    elif protype == "pod5":
        os.chdir("fastq_analysis")
        root = os.getcwd()
        _set_wkdir(root)
        if not os.path.isdir(f"{cfg.ofn}/fastq_analysis/basecaller_outputs/"):
            basecaller("pod5", cfg.intype, cfg.inf)
        if cfg.barkit != "none":
            subprocess.run(f"/home/dell/biosoft/ont-guppy/bin/guppy_barcoder -i basecaller_outputs/*.f*q* -s barout -x auto --barcode_kits {cfg.barkit}", shell=True)
            anum = len(os.listdir("barout"))
            snum = 1
            for item in os.listdir("barout"):
                if item.startswith("barcode"):
                    sample_dir = _ensure_sample_dir(root, item)
                    os.chdir(sample_dir)
                    _set_wkdir(os.getcwd())
                    _record_sample(cfg.ofn, item)
                    if os.path.isdir(f"{cfg.ofn}/fastq_analysis/basecaller_outputs/"):
                        subprocess.run(f"cat {cfg.ofn}/fastq_analysis/barout/{item}/*.f*q*|seqkit seq > {item}.raw.fastq", shell=True)
                    try:
                        _run_sample(main_process, cfg, f"{item}.raw.fastq", cfg.fastq1, cfg.fastq2, item, anum, snum, default_llid)
                    except Exception:
                        _print_failure(item, snum, anum)
                    snum += 1
                    os.chdir(root)
                    _set_wkdir(root)
            os.chdir(cfg.ofn)
            _set_wkdir(cfg.ofn)
        else:
            sample = "sample1"
            anum = 1
            snum = 1
            os.makedirs(sample, exist_ok=True)
            os.chdir(sample)
            _set_wkdir(os.getcwd())
            _record_sample(cfg.ofn, sample)
            if os.path.isdir(f"{cfg.ofn}/fastq_analysis/basecaller_outputs/"):
                subprocess.run(f"cat {cfg.ofn}/fastq_analysis/basecaller_outputs/*.f*q*|seqkit seq > {sample}.raw.fastq", shell=True)
            try:
                _run_sample(main_process, cfg, f"{sample}.raw.fastq", cfg.fastq1, cfg.fastq2, sample, anum, snum, default_llid)
            except Exception:
                print("数据量不足")
            os.chdir(cfg.ofn)
            _set_wkdir(cfg.ofn)
            subprocess.run("cat */*/Sample_result.txt >> Sample_result.txt", shell=True)


def _run_list_mode(cfg: RunnerConfig, main_process) -> None:
    os.chdir("fastq_analysis")
    root = os.getcwd()
    _set_wkdir(root)
    snum = 1
    listdb = pd.read_table(cfg.inf)
    listdb["样本名称"] = listdb["样本名称"].astype("str")
    listdb.fillna(0, inplace=True)
    anum = listdb.shape[0]
    print("list分析模式")
    print(listdb)
    sys.stdout.flush()
    for line in listdb.index.tolist():
        print(line)
        sys.stdout.flush()
        pre = listdb.iloc[line,].样本名称
        tinf = listdb.iloc[line,].三代数据
        fastq1 = listdb.iloc[line,].二代数据左
        fastq2 = listdb.iloc[line,].二代数据右
        llid = listdb.iloc[line,].物种信息
        print(pre, tinf, fastq1, fastq2, llid)
        sys.stdout.flush()
        if tinf and cfg.intype == "fastq":
            if not _is_fasta_suffix(str(tinf)):
                protype = check_input(tinf, "fastq")
                print(protype)
                try:
                    sample_dir = _ensure_sample_dir(root, pre)
                    _record_sample(cfg.ofn, pre)
                    os.chdir(sample_dir)
                    _set_wkdir(os.getcwd())
                    if protype == "fqdir":
                        subprocess.run(f"cat {tinf}/*.f*q* |seqkit rmdup -i |seqkit seq > {pre}.raw.fastq", shell=True)
                        if os.path.getsize(f"{pre}.raw.fastq") == 0:
                            print(f"{pre}文件夹内没有三代fastq格式文件")
                            raise SystemExit()
                    elif protype == "fqfile":
                        subprocess.run(f"cat {tinf}|seqkit seq > {pre}.raw.fastq", shell=True)
                        if os.path.getsize(f"{pre}.raw.fastq") == 0:
                            print(f"{pre}三代fastq格式文件有误，请核对")
                            raise SystemExit()
                    _run_sample(main_process, cfg, f"{pre}.raw.fastq", fastq1, fastq2, pre, anum, snum, llid)
                except Exception:
                    _print_failure(pre, snum, anum)
                snum += 1
                os.chdir(cfg.ofn)
                _set_wkdir(cfg.ofn)
            else:
                protype = check_input(tinf, "fasta")
                if protype == "fadir":
                    faendlist = ["fa", "fasta", "fna", "fas", "fa.gz", "fasta.gz", "fna.gz", "fas.gz"]
                    fadir_anum = 0
                    for tfile in os.listdir(cfg.inf):
                        if [i for i in faendlist if tfile.endswith(i)]:
                            fadir_anum += 1
                    snum = 1
                    for tfile in os.listdir(cfg.inf):
                        endssuf = [i for i in faendlist if tfile.endswith(i)]
                        if endssuf:
                            sample = tfile.replace(f".{endssuf[0]}", "")
                            _record_sample(cfg.ofn, sample)
                            sample_dir = _ensure_sample_dir(root, sample)
                            os.chdir(sample_dir)
                            _set_wkdir(os.getcwd())
                            subprocess.run(f"seqkit seq -i -w 0  {tinf} > {sample}.final.fasta", shell=True)
                            try:
                                _run_sample(main_process, cfg, f"{sample}.final.fasta", fastq1, fastq2, sample, fadir_anum, snum, llid)
                            except Exception:
                                _print_failure(sample, snum, fadir_anum)
                            snum += 1
                            os.chdir(cfg.ofn)
                            _set_wkdir(cfg.ofn)
                elif protype == "fafile":
                    _record_sample(cfg.ofn, pre)
                    sample_dir = _ensure_sample_dir(root, pre)
                    os.chdir(sample_dir)
                    _set_wkdir(os.getcwd())
                    subprocess.run(f"seqkit seq -i -w 0 {tinf} > {pre}.final.fasta", shell=True)
                    _run_sample(main_process, cfg, f"{pre}.final.fasta", fastq1, fastq2, pre, anum, snum, llid)
                    os.chdir(cfg.ofn)
                    _set_wkdir(cfg.ofn)
                    snum += 1
        else:
            _record_sample(cfg.ofn, pre)
            sample_dir = _ensure_sample_dir(root, pre)
            os.chdir(sample_dir)
            _set_wkdir(os.getcwd())
            subprocess.run(f"touch {pre}.raw.fastq", shell=True)
            _run_sample(main_process, cfg, 0, fastq1, fastq2, pre, anum, snum, llid)
            snum += 1
            os.chdir(cfg.ofn)
            _set_wkdir(cfg.ofn)


def run_pipeline_entry(cfg: RunnerConfig, main_process) -> None:
    root = _init_analysis_root(cfg.ofn)
    os.chdir(cfg.ofn)
    _set_wkdir(cfg.ofn)
    raw_input = cfg.raw_input
    if (os.path.isfile(raw_input) and os.popen(f"less {raw_input}|head -n 1").read().strip().startswith("@")) or os.path.isdir(raw_input) or is_fasta(raw_input):
        _run_non_list_mode(cfg, main_process)
    else:
        _run_list_mode(cfg, main_process)
    os.chdir(cfg.ofn)
    _set_wkdir(cfg.ofn)
