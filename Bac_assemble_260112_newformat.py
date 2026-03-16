#!/home/dell/miniconda3/envs/TB_ONT/bin/python
# 规定了应使用的Python解释器的路径
import pandas as pd
import os
import subprocess
import argparse
import sys 
import time
import re
import multiprocessing
from metagenomic_refactor.config import build_pipeline_config
__author__='wsh'
__version__='1.0.0'
#----version new!------
#组装模块1.0.0
__date__='20260306'
parser = argparse.ArgumentParser(description='纯菌+宏测组装分析脚本')      # 创建一个名为'纯菌+宏测组装分析脚本'的解析器（对象）。通过parser.add_argument定义命令行参数， parser.parse_args() 获取并解析传入的参数
parser.add_argument('--input','-i',type=str,default=False,help='输入文件路径')             # -i,输入文件路径。type=str，输入的值会被解析为'字符串'
parser.add_argument('--inputtype','-tp',type=str,default='fastq',help='输入文件类型')      # -tp,默认为fastq
parser.add_argument('--minlongfilt','-ml',type=str,default=500,help='最短序列长度')        # -ml,默认最短长度为500
parser.add_argument('--Qfilt','-mq',type=str,default=10,help='最小序列质量')               # -mq,默认最小序列质量10
parser.add_argument('--barcodekit','-bk',type=str,default='none',help='拆分试剂盒')       # -bk,默认不拆分试剂盒
parser.add_argument('--thread','-t',type=int,default=10,help='线程数量')                  # -t,默认线程数量10。type=int,输入的值会被解析为'整数'
parser.add_argument('--output','-o',type=str,default='tmpdir1',help='输出文件')           # -o,默认输出文件夹为tmpdir1
parser.add_argument('--fake_pip','-f',type=int,default=0,help='是否是假流程')             # -f,默认为0，非假流程
parser.add_argument('--method','-m',type=str,default='spades,freebayes',help='组装软件')  # -m,输入组装软件，默认组装软件为spades。可以输入两个spades,freebayes（二代）或clair3（三代）.有参组装输入：freebayes
parser.add_argument('--long_type','-lt',type=str,default='Nanopore',help='数据类型')      # -lt,默认数据类型为Nanopore
parser.add_argument('--ref','-r',type=str,default='noref',help='输入参考基因组')           # -r,输入参考基因组路径/名称
parser.add_argument('--gtf','-gtf',type=str,default='nogtf',help='输入注释文件')           # -gtf，输入注释的.gff文件。gff文件比gtf更加简洁，GCF：RefSeq标准化和改进的基因组
parser.add_argument('--genome_len','-gl',type=str,default='4m',help='预估基因组大小')      # -gl,默认基因组大小4m
parser.add_argument('--asm_type','-at',type=str,default='shortasm',help='组装方式')        # -at,组装类型。默认为shortasm.包括：'shortref'，'longref'，'longasm'，'shortasm'，'shortlongasm'，'shortlongasm'
parser.add_argument('--polish_times','-pt',type=str,default='1',help='抛光次数')          # -pt,默认抛光次数1次
parser.add_argument('--polish_soft','-ps',type=str,default='medaka',help='抛光软件')      # -ps,默认抛光软件medaka
parser.add_argument('--species','-species',type=str,default='False',help='物种信息')    # -species,默认物种信息为False
parser.add_argument('--ifAnno','-ifanno',type=str,default='Anno',help='是否注释')         # -ifanno，默认进行注释
parser.add_argument('--rmhost','-rmh',type=str,default='norm',help='是否去宿主')          # -rmh，默认不去宿主.norm不去宿主，其他情况去
parser.add_argument('--runflow','-rf',type=str,default='All',help='运行节点')             # -rf，运行节点:基因组组装，物种鉴定，结构变异检测，功能注释，元件预测，耐药与毒力，mlst与血清型
parser.add_argument('--abun','-abun',type=str,default='1',help='过滤阈值')             # -abum,默认过滤阈值0.85
parser.add_argument('--rna','-rna',type=str,default='0',help='是否rna建库')              # rna为0不是rna建库。rna != 0且-rmh ！= norm，进入rna分析
argv = parser.parse_args()   # 解析输入的参数并存储到argv中
command_line = " ".join(sys.argv)  
print(command_line)  # 打印命令行参数
sys.stdout.flush()   # 刷新标准输出

pipeline_config = build_pipeline_config(argv, multiprocessing.cpu_count())

# Legacy functions still read this subset of globals directly.
inf = pipeline_config.inf
ofn = pipeline_config.ofn
nt = pipeline_config.nt
rnalib = pipeline_config.rnalib
Krdb = pipeline_config.krdb
intype = pipeline_config.intype
minl = pipeline_config.minl
minQ = pipeline_config.minQ
barkit = pipeline_config.barkit
tmpfake = pipeline_config.tmpfake
long_type = pipeline_config.long_type
method = pipeline_config.method
asm_type = pipeline_config.asm_type
ref = pipeline_config.ref
gtf = pipeline_config.gtf
runflow = pipeline_config.runflow
rmhost = pipeline_config.rmhost
tspeabun = pipeline_config.tspeabun
genome_len = pipeline_config.genome_len
ptimes = pipeline_config.ptimes
psoft = pipeline_config.psoft
species = pipeline_config.species
vfmeta = pipeline_config.vfmeta

print(pipeline_config.config_out)  # 打印传入的全部参数
sys.stdout.flush()

# 创建并进入输出文件夹
if not os.path.isdir(ofn): # isdir,是否为路径/文件夹
    os.makedirs(ofn)       # 创建路径/文件夹
os.chdir(ofn)              # 切换到（ofn）路径下

from metagenomic_refactor.assembly import (
    Annotate_func,
    asb_func,
)
from metagenomic_refactor.annotation import AnnoEle, AnnoFun, DrugFinder, VFDR
from metagenomic_refactor.context import RuntimeContext, set_runtime_context
from metagenomic_refactor.runner import RunnerConfig, run_pipeline_entry
from metagenomic_refactor.taxonomy import kk2
from metagenomic_refactor.typing import mlst_serotype
from metagenomic_refactor.workflow import main_process, register_workflow_callbacks

set_runtime_context(
    RuntimeContext(
        ofn=ofn,
        runflow=runflow,
        method=method,
        rmhost=rmhost,
        tspeabun=tspeabun,
        nt=nt,
        krdb=Krdb,
        wkdir=ofn,
        long_type=long_type,
        species=species,
        genome_len=genome_len,
        ref=ref,
        gtf=gtf,
        vfmeta=vfmeta,
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
    mlst_serotype=mlst_serotype,
)

def main():
    run_pipeline_entry(
        RunnerConfig(
            raw_input=argv.input,
            inf=inf,
            intype=intype,
            ofn=ofn,
            barkit=barkit,
            tmpfake=tmpfake,
            fastq1=pipeline_config.fastq1,
            fastq2=pipeline_config.fastq2,
            nt=nt,
            llid=species if species != 'False' else 'nolevel',
            mmethod=pipeline_config.mmethod,
            minl=minl,
            minQ=minQ,
            asm_type=asm_type,
            ptimes=ptimes,
            psoft=psoft,
            rnalib=rnalib,
            ref=ref,
            gtf=gtf,
        ),
        main_process=main_process,
    )
    print('所有分析进程已经完成')


if __name__ == "__main__":
    main()
