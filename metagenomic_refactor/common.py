from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import glob
from Bio import SeqIO


def is_fasta(filename):
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)


def is_fastq(file_path):
    if file_path != 0:
        suflist = ["fastq", "fq", "fastq.gz", "fq.gz"]
        tlist = [i for i in suflist if str(file_path).endswith(i)]
        return len(tlist) > 0
    return False


def run_cmd(cmd, logf=None):
    subprocess.run(cmd, shell=True, stdout=logf, stderr=logf)


def copy_pattern(patterns, dest):
    dest = Path(dest)
    dest.mkdir(parents=True, exist_ok=True)

    for p in patterns:
        for f in glob.glob(p):
            try:
                shutil.copy(f, dest)
            except Exception:
                pass


def format_seconds(seconds):
    total_seconds = int(seconds)
    hours = total_seconds // 3600
    minutes = (total_seconds % 3600) // 60
    seconds = total_seconds % 60
    return f"{hours}小时{minutes}分钟{seconds}秒"


def file_exists(path):
    return bool(path) and str(path) != "0" and os.path.isfile(path)


def is_fastq_input(infile, fastq1, fastq2):
    return is_fastq(infile) or is_fastq(fastq1) or is_fastq(fastq2)


def check_input(infile, intp):
    if intp == "barcode_fastq":
        if os.path.isdir(infile):
            for i in os.listdir(infile):
                if i.startswith("barcode"):
                    return "bardir"
            print("文件夹下未检测到barcode文件夹")
            raise SystemExit()
        print("请输入barcode_fastq文件夹")
        raise SystemExit()

    if intp == "fastq":
        if os.path.isdir(infile):
            return "fqdir"
        if os.path.isfile(infile):
            return "fqfile"

    elif ".cfg" in intp:
        if os.path.isdir(infile):
            return "f5dir"
        print("请输入f5dir文件夹")
        raise SystemExit()

    elif "@" in intp:
        if os.path.isdir(infile):
            return "pod5"
        print("请输入pod5文件夹")
        raise SystemExit()

    elif intp == "fasta":
        if os.path.isdir(infile):
            return "fadir"
        return "fafile"

    return None


def get_free_gpu_memory():
    try:
        result = subprocess.run(
            ["nvidia-smi", "--query-gpu=memory.free", "--format=csv,noheader,nounits"],
            stdout=subprocess.PIPE,
        )
        memory_free = result.stdout.decode("utf-8").strip()
        return memory_free or "noGPU"
    except Exception:
        return "noGPU"


def basecaller(inputtype, params, inf, ofn="basecaller_outputs"):
    freem = get_free_gpu_memory()
    if freem == "noGPU":
        return

    freem = int(freem)
    if not os.path.isdir(ofn):
        os.makedirs(ofn)

    if os.path.isfile(inf):
        os.makedirs(f"{ofn}/tmp")
        subprocess.run(f"cp {inf} {ofn}/tmp", shell=True)
        inf = f"{ofn}/tmp"

    if inputtype == "fast5":
        print("下机数据为fast5模式")
        print(f"guppy_basecaller -r -i {inf} -s {ofn} -x auto -c {params}")
        if freem / 1.2 > 5500:
            subprocess.run(f"guppy_basecaller -r -i {inf} -s {ofn} -x auto -c {params}", shell=True)
        else:
            subprocess.run(
                f"guppy_basecaller -r -i {inf} -s {ofn} -x auto -c {params} --chunk_size 1000",
                shell=True,
            )
        subprocess.run(f"seqkit seq {ofn}/pass/*.fastq > {ofn}/basecaller.fastq", shell=True)
        return

    print("下机数据为pod5模式")
    qdict = {"fast@": 8, "hac@": 9, "sup@": 10}
    configdir = "/home/dell/biosoft/dorado"
    paramlist = params.split("_")
    if len(paramlist) == 4:
        _, _, _, modeln = paramlist
        qthod = qdict.get(modeln)
        subprocess.run(
            f"{configdir}/bin/dorado basecaller -r {configdir}/{params}v3.3 {inf} --min-qscore {qthod} --emit-fastq > {ofn}/basecaller.fastq",
            shell=True,
        )
    else:
        _, _, _, cspeed, modeln = paramlist
        qthod = qdict.get(modeln)
        if cspeed == "400bps":
            subprocess.run(
                f"{configdir}/bin/dorado basecaller -r {configdir}/{params}v4.1.0 {inf} --min-qscore {qthod} --emit-fastq  > {ofn}/basecaller.fastq",
                shell=True,
            )
            if os.path.getsize(f"{ofn}/basecaller.fastq") == 0:
                subprocess.run(
                    f"{configdir}/bin/dorado basecaller -r {configdir}/{params}v4.2.0 {inf} --min-qscore {qthod} --emit-fastq  > {ofn}/basecaller.fastq",
                    shell=True,
                )
        else:
            subprocess.run(
                f"{configdir}/bin/dorado basecaller -r {configdir}/{params}v4.1.0 {inf} --min-qscore {qthod} --emit-fastq  > {ofn}/basecaller.fastq",
                shell=True,
            )
