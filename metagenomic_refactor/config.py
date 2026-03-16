from __future__ import annotations

import os
import subprocess
from dataclasses import dataclass

import pandas as pd


@dataclass
class PipelineConfig:
    inf: str
    ofn: str
    sumpath: str
    nt: int
    rnalib: str
    krdb: str
    intype: str
    minl: str
    minQ: str
    barkit: str
    tmpfake: int
    long_type: str
    method: str
    asm_type: str
    ref: str
    gtf: str
    runflow: str
    rmhost: str
    tspeabun: str
    speciesdb: pd.DataFrame
    IfAnno: str
    sc1: str
    refdict: str
    genome_len: str
    ptimes: str
    psoft: str
    fastq1: int
    fastq2: int
    species: str
    speciesrefdb: pd.DataFrame
    vfmeta: pd.DataFrame
    config_out: str
    mmethod: str


def build_pipeline_config(argv, cpu_count: int) -> PipelineConfig:
    if argv.input:
        inf = os.path.abspath(argv.input)
    else:
        raise SystemExit("未检测到输入数据")

    ofn = os.path.abspath(argv.output)
    sumpath = f"{ofn}/fastq_analysis"

    nt = argv.thread
    if nt > cpu_count:
        nt = cpu_count

    krdb = "/home/dell/kraken2_custom_202101_24G"
    inputtypedic = {"1": ".cfg", "2": "fastq", "3": "barcode_fastq", "4": "fasta", "5": "@"}
    intype = inputtypedic.get(argv.inputtype, argv.inputtype)
    speciesdb = pd.read_table("/data1/shanghai_pip/meta_genome/pathotable.tsv")

    ref = argv.ref
    refdict = ""
    if "ref" in argv.asm_type:
        if os.path.isfile(ref):
            ref = os.path.abspath(ref)
        elif ref in ["salmonella", "E_coli", "Shigella", "Parahemolyticus", "cholerae", "Y_enterocolitica", "Campylobacter", "Brucella", "Lmono", "Kpne", "Suare", "Bcere", "Nmen", "HPinf"]:
            ref = speciesdb.loc[speciesdb["Species"] == ref, "reference"].tolist()[0]
        else:
            tref = f"/data1/shanghai_pip/meta_genome/database/fadb/{ref}_genomic.fna.gz"
            subprocess.run(f"seqkit seq {tref} > tmp_ref.fa", shell=True)
            ref = f"{os.getcwd()}/tmp_ref.fa"

    speciesrefdb = pd.read_table("/data1/shanghai_pip/meta_genome/gc_gtdbmeta.tsv")
    vfmeta = pd.read_table("/data1/shanghai_pip/meta_genome/database/vfdb/VFs_meta.tsv", encoding="Windows-1252")

    config_out = f"""输入文件\t{inf}
输出文件\t{ofn}
输入文件类型\t{intype}
线程数:\t{nt}
最小序列长度:\t{argv.minlongfilt}
最小序列质量:\t{argv.Qfilt}
试剂盒:\t{argv.barcodekit}
测序数据类型:\t{argv.asm_type}
长读长类型:\t{argv.long_type}
组装方法:\t{argv.method}
预估基因组大小:\t{argv.genome_len}
polish软件:\t{argv.polish_soft}
polish次数:\t{argv.polish_times}
物种:\t{argv.species}
参考基因组:\t{ref}
运行模块:\t{argv.runflow}
"""

    return PipelineConfig(
        inf=inf,
        ofn=ofn,
        sumpath=sumpath,
        nt=nt,
        rnalib=argv.rna,
        krdb=krdb,
        intype=intype,
        minl=argv.minlongfilt,
        minQ=argv.Qfilt,
        barkit=argv.barcodekit,
        tmpfake=argv.fake_pip,
        long_type=argv.long_type,
        method=argv.method,
        asm_type=argv.asm_type,
        ref=ref,
        gtf=argv.gtf,
        runflow=argv.runflow,
        rmhost=argv.rmhost,
        tspeabun=argv.abun,
        speciesdb=speciesdb,
        IfAnno=argv.ifAnno,
        sc1="/data/deploy/TB_soft/other_soft/3_kreport2krona.py",
        refdict=refdict,
        genome_len=argv.genome_len,
        ptimes=argv.polish_times,
        psoft=argv.polish_soft,
        fastq1=0,
        fastq2=0,
        species=argv.species,
        speciesrefdb=speciesrefdb,
        vfmeta=vfmeta,
        config_out=config_out,
        mmethod=argv.method,
    )
