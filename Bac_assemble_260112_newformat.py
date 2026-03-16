#!/home/dell/miniconda3/envs/TB_ONT/bin/python
# 规定了应使用的Python解释器的路径
import pandas as pd
import os
import subprocess
import argparse
import sys 
import numpy as np
import time
import re
from Bio import SeqIO
import multiprocessing
from Bio import Phylo
import itertools
import glob
import math
import json
import pathogenprofiler as pp
import csv
import ast
from typing import Iterable, Set, Tuple, Dict, List
import pytaxonkit
import shutil
import logging
from dataclasses import dataclass
from typing import Dict, Any, Optional
import glob
from pathlib import Path
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

# 将输入路径转换为绝对路径
if argv.input:
    inf = os.path.abspath(argv.input)
else:
    print('未检测到输入数据')
    sys.exit()    # 退出执行
ofn = os.path.abspath(argv.output)  # 将输出路径转换为绝对路径
sumpath = f'{ofn}/fastq_analysis'   # 总结分析数据的保存路径

# 合理控制调用CPU数目
nt = argv.thread                           # nt存储输入的线程数
num_threads = multiprocessing.cpu_count()  # 获取当前可用cpu数
if nt > num_threads:
    nt = num_threads
rnalib = argv.rna
Krdb='/home/dell/kraken2_custom_202101_24G'
# 获取输入文件类型   
inputtypedic = {'1':'.cfg','2':'fastq','3':'barcode_fastq','4':'fasta','5':'@'}  # 定义一个字典inputtypedic
intype = inputtypedic.get(argv.inputtype,argv.inputtype)  # 根据填写的数字（键）找到对应的文件类型(值）或直接选取输入的文字
minl = argv.minlongfilt
minQ = argv.Qfilt
barkit = argv.barcodekit
tmpfake = argv.fake_pip
long_type = argv.long_type
method = argv.method
asm_type = argv.asm_type
ref = argv.ref
gtf = argv.gtf
runflow = argv.runflow
rmhost = argv.rmhost
tspeabun = argv.abun        
speciesdb = pd.read_table('/data1/shanghai_pip/meta_genome/pathotable.tsv')  # speciesdb = ref等的预设文件夹
IfAnno = argv.ifAnno
sc1='/data/deploy/TB_soft/other_soft/3_kreport2krona.py'  # 定义脚本路径
# ref值传入
refdict = ''
if 'ref' in asm_type:           # 如果进行有参组装
    if os.path.isfile(ref):     # isfile,是否存在文件  
        ref = os.path.abspath(ref) # 转化为绝对路径
    elif ref in ['salmonella','E_coli','Shigella','Parahemolyticus','cholerae','Y_enterocolitica','Campylobacter','Brucella','Lmono','Kpne','Suare','Bcere','Nmen','HPinf']:
        ref = speciesdb.loc[speciesdb['Species']==ref,'reference'].tolist()[0]        # 从speciesdb中提取对应物种'ref'的'reference'列的值
    else:
        tref = f'/data1/shanghai_pip/meta_genome/database/fadb/{ref}_genomic.fna.gz'  # 自己传入参考序列时，放置的路径
        # 将输入的文件，解压并复制到tmp_ref.fa
        subprocess.run(f'seqkit seq {tref} > tmp_ref.fa',shell=True) # subprocess.run,在python内部执行外部命令。f’是一种字符串格式化方式，subprocess.run('seqkit seq '+tref+' > tmp_ref.fa',shell=True)
        ref = f'{os.getcwd()}/tmp_ref.fa'  # 生成绝对路径
genome_len = argv.genome_len
ptimes = argv.polish_times
psoft = argv.polish_soft
fastq1 = 0
fastq2 = 0
species = argv.species
speciesrefdb = pd.read_table('/data1/shanghai_pip/meta_genome/gc_gtdbmeta.tsv')  # speciesrefdb = 传入GC、GTDB数据库
vfmeta = pd.read_table('/data1/shanghai_pip/meta_genome/database/vfdb/VFs_meta.tsv', encoding='Windows-1252')   # vfmeta = VFDB数据库

# 输出配置参数
config_out = f'''输入文件\t{inf}
输出文件\t{ofn}
输入文件类型\t{intype}
线程数:\t{nt}
最小序列长度:\t{minl}
最小序列质量:\t{minQ}
试剂盒:\t{barkit}
测序数据类型:\t{asm_type}
长读长类型:\t{long_type}
组装方法:\t{method}
预估基因组大小:\t{genome_len}
polish软件:\t{psoft}
polish次数:\t{ptimes}
物种:\t{species}
参考基因组:\t{ref}
运行模块:\t{runflow}
'''
mmethod = method   
print(config_out)  # 打印传入的全部参数
sys.stdout.flush()

# 创建并进入输出文件夹
if not os.path.isdir(ofn): # isdir,是否为路径/文件夹
    os.makedirs(ofn)       # 创建路径/文件夹
os.chdir(ofn)              # 切换到（ofn）路径下


#--   module -----
# 定义is_fasta函数和is_fastq函数：判断文件格式
def is_fasta(filename):
    with open(filename, "r") as handle:  # r代表只读，w可修改
        fasta = SeqIO.parse(handle, "fasta")  # 使用 SeqIO.parse 按照 FASTA 格式解析文件
        return any(fasta)    # 至少有一条符合条件的序列返回true
def is_fastq(file_path): 
    if file_path != 0:
        suflist = ['fastq','fq','fastq.gz','fq.gz']
        tlist = [i for i in suflist if file_path.endswith(i)]   # 通过检查文件的扩展名来判断文件类型
        if len(tlist) > 0: 
            return True
        else:
            return False
    else:
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
            except:
                pass

# ---------- logging ----------
def get_logger(name="qc", level=logging.INFO):
    logger = logging.getLogger(name)
    if not logger.handlers:
        logger.setLevel(level)
        h = logging.StreamHandler()
        fmt = logging.Formatter("[%(asctime)s] %(levelname)s %(name)s: %(message)s")
        h.setFormatter(fmt)
        logger.addHandler(h)
    return logger


# ---------- helpers ----------
def safe_read_json(path: str, logger: logging.Logger) -> Optional[Dict[str, Any]]:
    try:
        with open(path) as f:
            return json.load(f)
    except FileNotFoundError:
        logger.warning(f"JSON not found: {path}")
    except json.JSONDecodeError as e:
        logger.error(f"JSON parse error: {path} ({e})")
    return None

# ---------- logging ----------
def get_logger(name="qc", level=logging.INFO):
    logger = logging.getLogger(name)
    if not logger.handlers:
        logger.setLevel(level)
        h = logging.StreamHandler()
        fmt = logging.Formatter("[%(asctime)s] %(levelname)s %(name)s: %(message)s")
        h.setFormatter(fmt)
        logger.addHandler(h)
    return logger


# ---------- helpers ----------
def safe_read_json(path: str, logger: logging.Logger) -> Optional[Dict[str, Any]]:
    try:
        with open(path) as f:
            return json.load(f)
    except FileNotFoundError:
        logger.warning(f"JSON not found: {path}")
    except json.JSONDecodeError as e:
        logger.error(f"JSON parse error: {path} ({e})")
    return None
def read_fastqc_data(path: str, logger: logging.Logger) -> Optional[Dict[str, str]]:
    """
    读取 fastqc_data.txt 形成 key->value 字典。
    """
    try:
        res = {}
        with open(path) as f:
            for line in f:
                if "\t" not in line:
                    continue
                k, v = line.rstrip("\n").split("\t")[:2]
                res[k] = v
        return res
    except FileNotFoundError:
        logger.warning(f"fastqc_data.txt not found: {path}")
        return None


def normalize_fastqc_images(fastqc_dir: str, logger: logging.Logger):
    """
    将 fastqc Images 内部文件名标准化（避免不同版本 fastqc 图片命名差异导致报告找不到图）。
    """
    imgdir = os.path.join(fastqc_dir, "Images")
    if not os.path.isdir(imgdir):
        logger.warning(f"Images dir missing: {imgdir}")
        return

    mapping = {
        "per_base_quality.png": "per_base_sequence_quality.png",
        "per_sequence_quality.png": "per_sequence_quality_scores.png",
        "duplication_levels.png": "sequence_duplication_levels.png",
    }

    for old, new in mapping.items():
        src = os.path.join(imgdir, old)
        dst = os.path.join(imgdir, new)
        if os.path.exists(src):
            try:
                shutil.move(src, dst)
            except Exception as e:
                logger.error(f"Move failed {src} -> {dst}: {e}")


def normalize_summary_txt(path: str, logger: logging.Logger):
    """
    fastqc summary.txt 通常第一行和最后一行是特殊条目，很多人会 drop。
    这里做成鲁棒版：>=9行 drop[0,8]；否则 drop[0,last]
    """
    if not os.path.exists(path):
        logger.warning(f"summary.txt not found: {path}")
        return

    try:
        df = pd.read_table(path, names=["status", "names", "tt"])
        if df.shape[0] == 0:
            logger.warning(f"summary.txt empty: {path}")
            return

        if df.shape[0] >= 9:
            df = df.drop([0, 8])
        else:
            df = df.drop([0, df.shape[0] - 1])

        df.to_csv(path, sep="\t", index=False, header=False)
    except Exception as e:
        logger.error(f"Normalize summary.txt failed: {path} ({e})")


def to_gb_str(x: Any) -> str:
    try:
        return f"{round(float(x) / 1e9, 2)} G"
    except Exception:
        return "NA"


def safe_float(x: Any, ndigits=2) -> Any:
    try:
        return round(float(x), ndigits)
    except Exception:
        return "NA"


def safe_int(x: Any) -> Any:
    try:
        return int(x)
    except Exception:
        return "NA"


def pick(d: Dict[str, Any], *keys, default=None):
    """
    连续取嵌套key：pick(d,'a','b','c')
    """
    cur = d
    for k in keys:
        if not isinstance(cur, dict) or k not in cur:
            return default
        cur = cur[k]
    return cur


# ---------- config ----------
@dataclass(frozen=True)
class QCTarget:
    name: str                 # Nano / R1 / R2
    trigger_dir: str          # Nano_qc / R1_qc / R2_qc
    raw_fastqc_dir: str       # raw_fastqc / raw.R1_fastqc / raw.R2_fastqc
    final_fastqc_dir: str     # {Pre}.final_fastqc / {Pre}.R1_fastqc / {Pre}.R2_fastqc
    json_path: str            # json文件路径模板（含Pre）
    mode: str                 # "nano" / "read1" / "read2"
    out_tsv: str              # 输出tsv文件名


# ---------- core ----------
def build_qc_table(
    raw_fastqc: Dict[str, str],
    final_fastqc: Dict[str, str],
    before: Dict[str, Any],
    after: Dict[str, Any],
) -> pd.DataFrame:
    qc = {"过滤前": {}, "过滤后": {}}

    # fastqc 基础字段
    qc["过滤前"]["总序列数"] = safe_int(raw_fastqc.get("Total Sequences"))
    qc["过滤前"]["总碱基数"] = raw_fastqc.get("Total Bases", "NA")
    qc["过滤前"]["低质量序列"] = safe_int(raw_fastqc.get("Sequences flagged as poor quality"))
    qc["过滤前"]["序列长度"] = raw_fastqc.get("Sequence length", "NA")
    qc["过滤前"]["GC%"] = raw_fastqc.get("%GC", "NA")

    qc["过滤后"]["总序列数"] = safe_int(final_fastqc.get("Total Sequences"))
    qc["过滤后"]["总碱基数"] = final_fastqc.get("Total Bases", "NA")
    qc["过滤后"]["低质量序列"] = safe_int(final_fastqc.get("Sequences flagged as poor quality"))
    qc["过滤后"]["序列长度"] = final_fastqc.get("Sequence length", "NA")
    qc["过滤后"]["GC%"] = final_fastqc.get("%GC", "NA")

    # q20/q30
    qc["过滤前"]["q20数据量"] = to_gb_str(before.get("q20_bases"))
    qc["过滤前"]["q20比例"] = safe_float(before.get("q20_rate"))

    qc["过滤前"]["q30数据量"] = to_gb_str(before.get("q30_bases"))
    qc["过滤前"]["q30比例"] = safe_float(before.get("q30_rate"))

    qc["过滤后"]["q20数据量"] = to_gb_str(after.get("q20_bases"))
    qc["过滤后"]["q20比例"] = safe_float(after.get("q20_rate"))

    qc["过滤后"]["q30数据量"] = to_gb_str(after.get("q30_bases"))
    qc["过滤后"]["q30比例"] = safe_float(after.get("q30_rate"))

    return pd.DataFrame(qc).T


def extract_before_after(fpjson: Dict[str, Any], mode: str, logger: logging.Logger):
    """
    mode:
      - nano: fpjson['summary']['before_filtering'], fpjson['summary']['after_filtering']
      - read1/read2: fpjson['read1_before_filtering'] ...
    """
    if mode == "nano":
        before = pick(fpjson, "summary", "before_filtering", default={}) or {}
        after = pick(fpjson, "summary", "after_filtering", default={}) or {}
        return before, after

    if mode in ("read1", "read2"):
        before = fpjson.get(f"{mode}_before_filtering", {}) or {}
        after = fpjson.get(f"{mode}_after_filtering", {}) or {}
        # 二代 fastp2.json 没有 q20_rate/q30_rate 的话，兜底用 bases/total_bases 计算
        for tag in (before, after):
            tb = tag.get("total_bases") or 0
            if "q20_rate" not in tag and tb:
                tag["q20_rate"] = (tag.get("q20_bases", 0) / tb)
            if "q30_rate" not in tag and tb:
                tag["q30_rate"] = (tag.get("q30_bases", 0) / tb)
        return before, after

    logger.error(f"Unknown mode: {mode}")
    return {}, {}


def run_one_target(pre: str, t: QCTarget, logger: logging.Logger) -> bool:
    """
    单个target的QC汇总。成功返回True。
    """
    if not os.path.isdir(t.trigger_dir):
        logger.info(f"Skip {t.name}: trigger dir not exists ({t.trigger_dir})")
        return False

    raw_fastqc_data_path = os.path.join(t.raw_fastqc_dir, "fastqc_data.txt")
    final_fastqc_data_path = os.path.join(t.final_fastqc_dir.format(pre=pre), "fastqc_data.txt")

    raw_fastqc = read_fastqc_data(raw_fastqc_data_path, logger)
    final_fastqc = read_fastqc_data(final_fastqc_data_path, logger)

    if raw_fastqc is None or final_fastqc is None:
        logger.warning(f"{t.name}: fastqc_data missing, skip table.")
        return False

    json_path = t.json_path.format(pre=pre)
    fpjson = safe_read_json(json_path, logger)
    if fpjson is None:
        logger.warning(f"{t.name}: fastp json missing, skip table.")
        return False

    before, after = extract_before_after(fpjson, t.mode, logger)
    df = build_qc_table(raw_fastqc, final_fastqc, before, after)

    out_path = t.out_tsv.format(pre=pre)
    df.to_csv(out_path, sep="\t")
    logger.info(f"{t.name}: write {out_path}")

    # 标准化图片与summary
    normalize_fastqc_images(t.raw_fastqc_dir, logger)
    normalize_fastqc_images(t.final_fastqc_dir.format(pre=pre), logger)

    normalize_summary_txt(os.path.join(t.raw_fastqc_dir, "summary.txt"), logger)
    normalize_summary_txt(os.path.join(t.final_fastqc_dir.format(pre=pre), "summary.txt"), logger)

    return True


def summaryfastqc_prod(pre: str, logger: Optional[logging.Logger] = None) -> Dict[str, bool]:
    """
    生产级入口：自动检测 Nano / R1 / R2 三种情况，分别输出对应TSV。
    返回各模块是否成功：{'Nano':True/False, 'R1':..., 'R2':...}
    """
    logger = logger or get_logger("summaryfastqc", logging.INFO)

    targets = [
        QCTarget(
            name="Nano",
            trigger_dir="Nano_qc",
            raw_fastqc_dir="raw_fastqc",
            final_fastqc_dir="{pre}.final_fastqc",
            json_path="{pre}.filter.fastp.json",   # 你原来用了 raw + filter 两个json；生产上推荐只用filter这个（含before/after）
            mode="nano",
            out_tsv="NanoFastqc.tsv",
        ),
        QCTarget(
            name="R1",
            trigger_dir="R1_qc",
            raw_fastqc_dir="raw.R1_fastqc",
            final_fastqc_dir="{pre}.R1_fastqc",
            json_path="{pre}.fastp2.json",
            mode="read1",
            out_tsv="R1_Fastqc.tsv",
        ),
        QCTarget(
            name="R2",
            trigger_dir="R2_qc",
            raw_fastqc_dir="raw.R2_fastqc",
            final_fastqc_dir="{pre}.R2_fastqc",
            json_path="{pre}.fastp2.json",
            mode="read2",
            out_tsv="R2_Fastqc.tsv",
        ),
    ]

    results = {}
    for t in targets:
        results[t.name] = run_one_target(pre, t, logger)

    return results

# 定义proc_kra函数：指定一个分类等级（lel:G、S等）和物种ID（tax），返回原先指定的物种ID + 其下分类所有的物种ID 
def proc_kra(kraken,tax,lel):   # raw_{Pre}.report.txt，tmpfile.taxonomy_id.tolist()[0]，'G'
    tmplist = [tax] 
    if lel in ['R','D','K','P','C','O','F','G','S']:    
        rawlist = ['R','D','K','P','C','O','F','G','S']
    else:  # 亚种
        rawlist = ['S1','S2','S3','S4','S5','S6']
    if tax != 0:  # =0表示未被鉴定的序列
        kradb = pd.read_table(kraken,header=None)  
        #kradb = kradb[(kradb[3]==lel)&(kradb[4]==tax)]
        tmpindex = kradb[(kradb[3]==lel)&(kradb[4]==tax)].index.tolist()[0]+1  # 选取同时满足lel和tax的行，+1为真正所在的行
        if tmpindex <= kradb.shape[0]-1:  # 确认tmpindex不为最后一行
            def getlindex(tmpindex):
                tmpl = kradb.iloc[tmpindex,3]
                for tl in rawlist:
                    if tl in tmpl:
                        if tl == tmpl:
                            tmpl = tl
                            tmlindex = rawlist.index(tmpl)
                        else:
                            tmpl = tl
                            tmlindex = rawlist.index(tmpl)+1  # 适用于lel=S1的情况。整个内函数最后实现，['R','D','K','P','C','O','F','G','S'] lel对应的索引+1
                        return tmlindex
            while getlindex(tmpindex) > rawlist.index(lel):   # while，只要条件成立会一直运行
                tmplist.append(kradb.iloc[tmpindex,4])        # .append()，用于将一个元素追加到列表的末尾
                tmpindex+=1
    return tmplist
# 定义exreadsID函数：根据taxlist中的物种ID，生成对应的reads ID，再使用 seqkit 提取
def exreadsID(Pre,taxlist,kraresult,fq1,fq2=0):  # Pre,taxlist1,raw_{Pre}.out.txt（包含物种ID及 reads ID）
    Maintax = taxlist[0]  # Maintax = taxlist列表中的第一个元素
    kraredb = pd.read_table(kraresult,header=None)
    pd.DataFrame(kraredb[kraredb[2].isin(taxlist)][1].unique()).to_csv(f'{Maintax}_tfqID.txt',sep='\n',index=False,header=False) # 选出所有taxlist中对应raw_{Pre}.out.txt行的 reads ID
    if os.popen(f'''head -n 1 {Maintax}_tfqID.txt''').read().strip().endswith('/1') or os.popen(f'''head -n 1 {Maintax}_tfqID.txt''').read().strip().endswith('/2'):  # 提取华大数据read ID，分成_fq1ID.txt和_fq2ID.txt
        subprocess.run(f'''cut -f1 -d '/' {Maintax}_tfqID.txt|sort -u|sed 's/$/\/1/' > {Maintax}_fq1ID.txt''',shell=True)
        subprocess.run(f'''cut -f1 -d '/' {Maintax}_tfqID.txt|sort -u|sed 's/$/\/2/' > {Maintax}_fq2ID.txt''',shell=True)
    else:  # 提取因美纳数据read ID，分成_fq1ID.txt和_fq2ID.txt
        subprocess.run(f'ln -s {Maintax}_tfqID.txt {Maintax}_fq1ID.txt',shell=True)
        subprocess.run(f'ln -s {Maintax}_tfqID.txt {Maintax}_fq2ID.txt',shell=True)
    if fq1.endswith('final.fastq'):  # 三代数据直接根据.txt进行序列提取
        subprocess.run(f'seqkit grep -n -f {Maintax}_fq1ID.txt {fq1} > {Pre}_t.final.fastq',shell=True)
    else: # 二代数据进行提取
        subprocess.run(f'seqkit grep -n -f {Maintax}_fq1ID.txt {fq1} > {Pre}.R1.fastq',shell=True)
        subprocess.run(f'gzip -f  {Pre}.R1.fastq',shell=True)
        if fq2:
            subprocess.run(f'seqkit grep -n -f {Maintax}_fq2ID.txt {fq2} > {Pre}.R2.fastq',shell=True)
            subprocess.run(f'gzip -f {Pre}.R2.fastq',shell=True)

# 定义ngs_qc函数：对所有二代以.fastq.gz结尾的文件用seqkit进行统计，生成{Pre}.QC2.summary.tsv
def ngs_qc(Pre):
    if not os.path.isfile('summary.tsv') or os.path.getsize('summary.tsv') == 0:
        subprocess.run(f'seqkit stat *_t.R*.fastq.gz -T -a -b > summary.tsv',shell=True)
    afile = pd.read_table('summary.tsv')
    afile = afile[['file','format','type','num_seqs','sum_len','min_len','avg_len','max_len','N50']]
    afile['file'] = afile['file'].str.split('.').str[1]
    afile.rename(columns={'num_seqs':'序列总数','sum_len':'碱基总数','min_len':'最短序列长度','avg_len':'平均序列长度','max_len':'最长序列长度'},inplace=True)
    afile.to_csv(f'{Pre}.QC2.summary.tsv',sep='\t',index=False)
# 定义QC_func函数：对二代或/和三代数据进行质控,并统计质控信息
def QC_func(inf,fq1,fq2,minq,minl,Pre,rnalib,threads,method):
    #1.QC Rawfastq 2.trim adapter barcode 3.filter low quality 4.rm Host contam 5.rm common contam 6.summary QC result
    #---检查三代测序数量
    if inf:
        if method != 'meta':
            subprocess.run(f'rasusa reads --bases 2gb -o {Pre}.sub.fastq {inf}',shell=True)
        else:
            subprocess.run(f'rasusa reads --bases 20gb -o {Pre}.sub.fastq {inf}',shell=True)
        subprocess.run(f'''fastp --in1 {Pre}.sub.fastq \\
                    --out1 {Pre}.clean.fastq \\
                    --thread {threads} \\
                    --length_required=50 \\
                    --n_base_limit=6 \\
                    --compression=6 \\
                    -Q \\
                    -A \\
                    --qualified_quality_phred=10 \\
                    --json {Pre}.raw.fastp.json \\
                    2> {Pre}.raw.fastp.log''',shell=True)  # 主要为了获得.json文件的统计信息.-A,-Q,不清洗。
        with open(f'QC.log','w') as f:     # 文件对象被赋值给变量f                
            subprocess.run(f'porechop -i {Pre}.clean.fastq -o {Pre}.trim.fastq -t {threads} --no_split --check_reads 2000',shell=True,stdout=f,stderr=f)  # 利用porechop去除三代序列中的适配子。输出：{Pre}.trim.fastq
        subprocess.run(f'seqkit seq -m {minl} -Q {minq} {Pre}.trim.fastq > {Pre}.filter.fastq',shell=True)  # 利用seqkit控制最短序列长度和质量，输出：{Pre}.filter.fastq
        if rmhost == 'norm':
            subprocess.run(f'ln -s  {Pre}.filter.fastq {Pre}.rmhost.fastq',shell=True)   # 默认rmhost = norm不进行去宿主过程。ln -s，创建一个符号链接Sample1.rmhost.fastq指向Sample1.filter.fastq
        else:
            # -n,指向hostile环境；--aligner，指定使用的比对工具minimap2
            subprocess.run(f'/home/dell/miniconda3/bin/conda run -n hostile hostile  clean --fastq1 {Pre}.filter.fastq --threads {threads} --aligner minimap2 --index /data/deploy/meta_new/Database/Host_Ref/hostile/human-t2t-hla.argos-bacteria-985_rs-viral-202401_ml-phage-202401.mmi > rmcont.json',shell=True)
            subprocess.run(f'seqkit seq {Pre}.filter.clean.fastq.gz > {Pre}.rmhost.fastq',shell=True)
            subprocess.run(f'rm  {Pre}.filter.clean.fastq.gz',shell=True)   # 删除过程文件
        subprocess.run(f'''fastp --in1  {Pre}.rmhost.fastq \\
                    --out1 {Pre}.final.fastq \\
                    --thread {threads} \\
                    --length_required=50 \\
                    --qualified_quality_phred=10 \\
                    --n_base_limit=6 \\
                    --compression=6 \\
                    -Q \\
                    -A \\
                    --json {Pre}.filter.fastp.json \\
                    2> {Pre}.filter.fastp.log''',shell=True)  # 再次使用fastp获取质控后的.json 文件统计信息                 
        subprocess.run(f'seqkit stat *.fastq -T -a -b > summary.tsv',shell=True)       # sekit统计所有fastq文件的基本信息                   
        # 提取统计信息并计算平均质量值
        afile = pd.read_table('summary.tsv')
        afile = afile[['file','format','type','num_seqs','sum_len','min_len','avg_len','max_len','N50']]  # 外层方括号：用于 DataFrame 的列选择；内层方括号：创建列名列表
        afile['file'] = afile['file'].str.split('.').str[1]  # afile['file'] = 将afile中的file列按照'.'进行分割，取第二部分
        afile.index = afile['file']
        afile.drop(['format','type','file'],axis=1,inplace=True)     # axis=1 表示删除列
        
        # 获取平均质量值，并添加到新列中
        calqdict = {}     # 创建一个空字典
        def get_meanQ(tlist):  
            newlist = [(1-10 ** (-i / 10)) for i in tlist]
            newq = round(-10 * math.log10(1-sum(newlist)/len(newlist)),2)
            return newq
        for status in afile.index:
            subprocess.run(f'seqkit fx2tab -n -i -q {Pre}.{status}.fastq > tmp.tab',shell=True)  # -fx2tab,将 FASTA/FASTQ 文件转换为制表符分隔的文本文件
            calqdb = pd.read_table('tmp.tab',header=None)
            calqdict[status] = get_meanQ(calqdb[1].tolist())           
        calqdb = pd.DataFrame(calqdict,index=['平均质量值']).T
        calqdb['file'] = calqdb.index
        calqdb.reset_index(drop=True,inplace=True)
        afile.reset_index(inplace=True)
        afile = afile.merge(calqdb,on='file')
        
        # {Pre}.QC.summary.tsv = 汇总不同过滤状态下的序列信息
        afile.rename(columns={'num_seqs':'序列总数','sum_len':'碱基总数','min_len':'最短序列长度','avg_len':'平均序列长度','max_len':'最长序列长度'},inplace=True)
        afile['file'] = afile['file'].replace({'raw':'原始序列','clean':'去除过短/质量过差','filter':'质量/长度过滤','trim':'去除接头/barcode后序列','rmhost':'去宿主后序列','final':'最终序列'}) #如果 afile['file'] 是 'raw'替换成 '原始序列'
        order = ['原始序列','去除过短/质量过差','去除接头/barcode后序列','质量/长度过滤','去宿主后序列','最终序列']
        afile = afile.set_index('file').reindex(order).reset_index()
        afile = afile.rename(columns={'file':'过滤状态'})
        afile.to_csv(f'{Pre}.QC.summary.tsv',sep='\t',index=False)
        # 清理中间文件
        #subprocess.run(f'rm {Pre}.clean.fastq {Pre}.trim.fastq {Pre}.filter.fastq {Pre}.rmhost.fastq',shell=True)
        
       # 用fastqc评估.raw和.final文件,最后输出到Nano_qc 文件夹中 
        subprocess.run(f'ln -s {inf} ./raw.fastq',shell=True)                
        if not os.path.isdir('Nano_qc'):
            os.makedirs('Nano_qc')
        with open('qc.log','a') as f:
            subprocess.run(f'fastqc raw.fastq {Pre}.final.fastq -o Nano_qc -t {threads}',shell=True,stdout=f,stderr=f)
            subprocess.run(f'''unzip -o 'Nano_qc/*.zip' ''',shell=True,stdout=f,stderr=f)
        subprocess.run(f'rm raw.fastq',shell=True)  
        
        #---当第一个属水平超过85% 就用全部数据进行组装，当低于这个数值是提取丰度第一的属的数据进行组装。结果记录在rawkk2.log
        with open('rawkk2.log','w') as rawkkf:
            subprocess.run(f'kraken2 --db /data1/shanghai_pip/meta_genome/database/kraken2_custom_202101_24G --threads {threads} --output raw_{Pre}.out.txt --report raw_{Pre}.report.txt {Pre}.final.fastq',shell=True,stdout=rawkkf,stderr=rawkkf)
            subprocess.run(f'bracken -d /data1/shanghai_pip/meta_genome/database/kraken2_custom_202101_24G -o raw_{Pre}.bracken1.txt -w raw_{Pre}.bracken2.txt -l S -t 10  -i raw_{Pre}.report.txt',shell=True,stdout=rawkkf,stderr=rawkkf)
            tmpfile = pd.read_table(f'raw_{Pre}.bracken1.txt',sep='\t')
            tmpfile = tmpfile[tmpfile['name']!='Homo']     # 去除宿主基因序列                           
            ONTSpe = tmpfile.name.tolist()[0]    # ONTSpe 和 tkid 分别获取 第一行 的'物种名称'和'分类ID'。可根据需要，更改提取的行数
            tkid = tmpfile.taxonomy_id.tolist()[0]  
            krakenfile = f'raw_{Pre}.report.txt'  # 存储 Kraken2 的分类报告文件路径
            tkidl = os.popen(f'grep {tkid} {krakenfile}|cut -f4').read().strip()   # tkid1 = krakenfile中tkid对应行的第四列信息，即种属水平
            if float(tmpfile.loc[tmpfile['name']==ONTSpe,'fraction_total_reads'].tolist()[0]) < float(tspeabun):       # 从bracken1.txt获取ONTSpe对应的物种丰度与阈值（0.85）进行比较
                 #taxlist1 = proc_kra(krakenfile,tkid,tkidl)                                                           
                 taxlist1 = proc_kra(krakenfile,517,'G')                  # taxlist1 = 包含所有对应物种层级（属），及其以下全部层级的物种ID.自己输入时注意G和前面数字层级对应
                 exreadsID(Pre,taxlist1,f'raw_{Pre}.out.txt',f'{Pre}.final.fastq')                                     # 根据taxlist中的物种ID，生成对应的reads ID，再使用 seqkit 提取
                 subprocess.run(f'mv {Pre}_t.final.fastq {Pre}.final.fastq',shell=True)                                     # {tkid}.1.fastq 的文件重命名(替换） {Pre}.final.fastq
    else:  # 无三代数据
        pass
    
    ## 二代数据    
    if fq1 and fq2:   # 开始处理二代数据
        if not os.path.isfile(f'{Pre}_sub.R1.fastq') or not os.path.isfile(f'{Pre}_sub.R2.fastq'):
        #-下采样，rasusa工具，--bases提取数据量
            if method != 'meta':
                subprocess.run(f'rasusa reads --bases 10gb -o {Pre}_sub.R1.fastq -o {Pre}_sub.R2.fastq {fq1} {fq2}',shell=True)      # R1+R2=10GB
            else:
                subprocess.run(f'rasusa reads --bases 100gb -o {Pre}_sub.R1.fastq -o {Pre}_sub.R2.fastq {fq1} {fq2}',shell=True)
        subprocess.run(f'seqkit stat -T {Pre}_sub.R1.fastq {Pre}_sub.R2.fastq > raw_summary.tsv',shell=True)
        if not os.path.isfile(f'{Pre}.fastp2.json'):
            subprocess.run(f'''fastp --in1 {Pre}_sub.R1.fastq \\
        --out1 {Pre}_t.R1.fastq.gz \\
        --in2 {Pre}_sub.R2.fastq \\
        --out2 {Pre}_t.R2.fastq.gz \\
        --thread {threads} \\
        --length_required=50 \\
        --n_base_limit=6 \\
        --compression=6 \\
        --detect_adapter_for_pe \\
        --json {Pre}.fastp2.json \\
        2> {Pre}.fastp2.log     
        ''',shell=True )   
        ngs_qc(Pre)    # 调用ngs_qc函数，对质控文件进行统计
        if not os.path.isdir('R1_qc'):
            os.makedirs('R1_qc')  
        if not os.path.isdir('R2_qc'):
            os.makedirs('R2_qc')
        with open('qc.log','a') as f:
            subprocess.run(f'ln -s  {Pre}_sub.R1.fastq  raw.R1.fastq',shell=True)        # 可以通过 raw.R1.fastq 和 raw.R2.fastq 来访问原始下采样后的 R1 和 R2 FASTQ 文件
            subprocess.run(f'ln -s  {Pre}_sub.R2.fastq raw.R2.fastq',shell=True)
            if not os.path.isfile(f'{Pre}.R1.fastq.gz'):
                if rmhost == 'norm':
                    subprocess.run(f'ln -s {Pre}_t.R1.fastq.gz  {Pre}.R1.fastq.gz;ln -s {Pre}_t.R2.fastq.gz  {Pre}.R2.fastq.gz',shell=True)
                else:
                    if rnalib == '0':
                        subprocess.run(f'/home/dell/miniconda3/bin/conda run -n hostile hostile clean --fastq1 {Pre}_t.R1.fastq.gz --threads {threads} --fastq2 {Pre}_t.R2.fastq.gz --aligner bowtie2 --index /data/deploy/meta_new/Database/Host_Ref/hostile/human-t2t-hla.argos-bacteria-985_rs-viral-202401_ml-phage-202401 > rmcont.json',shell=True,stdout=f,stderr=f)
                        subprocess.run(f'mv  {Pre}_t.R1.clean_1.fastq.gz  {Pre}.R1.fastq.gz',shell=True)
                        subprocess.run(f'mv  {Pre}_t.R2.clean_2.fastq.gz  {Pre}.R2.fastq.gz',shell=True)      # 将去宿主后的文件重命名(移动）为{Pre}.R1.fastq.gz
                    else:
                        print('测试rna建库')
                        if not os.path.isfile(f'kneaddata_out/{Pre}_t.R1_kneaddata_paired_1.fastq'):  # # 去除人源+rRNA序列
                            subprocess.run(f'/home/dell/miniconda3/bin/conda run -n kneaddata kneaddata --bypass-trf --bypass-trim -i1 {Pre}_t.R1.fastq.gz  -i2 {Pre}_t.R2.fastq.gz -o kneaddata_out -t {threads} -db /data/Ref/human_hg38_refMrna -db /data/Ref/SILVA_128_LSUParc_SSUParc_ribosomal_RNA' ,shell=True,stdout=f,stderr=f)  ## SILVA 第 128 版rRNA数据库
                        subprocess.run(f'pigz -p 10 kneaddata_out/{Pre}_t.R1_kneaddata_paired_1.fastq -c > {Pre}.R1.fastq.gz',shell=True,stdout=f,stderr=f)
                        subprocess.run(f'pigz -p 10 kneaddata_out/{Pre}_t.R1_kneaddata_paired_2.fastq -c > {Pre}.R2.fastq.gz',shell=True,stdout=f,stderr=f)
            if not os.path.isfile(f'R1_qc/raw.R1_fastqc.html'):   
                subprocess.run(f'fastqc raw.R1.fastq {Pre}.R1.fastq.gz -o R1_qc -t {threads}',shell=True,stdout=f,stderr=f)   # raw.R1.fastq为下采样后的数据，{Pre}.R1.fastq.gz为质控后的数据
                subprocess.run(f'fastqc raw.R2.fastq {Pre}.R2.fastq.gz -o R2_qc -t {threads}',shell=True,stdout=f,stderr=f)
                subprocess.run(f'''unzip -o 'R1_qc/*.zip' ''',shell=True,stdout=f,stderr=f)
                subprocess.run(f'''unzip -o 'R2_qc/*.zip' ''',shell=True,stdout=f,stderr=f)   
            subprocess.run(f'fastp -i {Pre}.R1.fastq.gz -I {Pre}.R2.fastq.gz -o tt.1.fq -O tt.2.fq -w {threads} --json {Pre}.final.json',shell=True,stdout=f,stderr=f)  
            subprocess.run(f'rm tt.*.fq',shell=True)  
        #    with open('rawkk2.log','w') as rawkkf:
        #    subprocess.run(f'kraken2 --db /data1/shanghai_pip/meta_genome/database/kraken2_custom_202101_24G --threads {threads} --output raw_{Pre}.out.txt --report raw_{Pre}.report.txt {Pre}.R1.fastq.gz {Pre}.R2.fastq.gz',shell=True,stdout=rawkkf,stderr=rawkkf)
        #    subprocess.run(f'bracken -d /data1/shanghai_pip/meta_genome/database/kraken2_custom_202101_24G -o raw_{Pre}.bracken1.txt -w raw_{Pre}.bracken2.txt -l G -t 10  -i raw_{Pre}.report.txt',shell=True,stdout=rawkkf,stderr=rawkkf)
        #    tmpfile = pd.read_table(f'raw_{Pre}.bracken1.txt')
        #    tmpfile = tmpfile[tmpfile['name']!='Homo']                                                                   # 去除宿主基因序列
        #    NGSspe = tmpfile.name.tolist()[0]
        #    tkid = tmpfile.taxonomy_id.tolist()[0]
        #    krakenfile = f'raw_{Pre}.report.txt'
        #    tkidl = os.popen(f'grep {tkid} {krakenfile}|cut -f4').read().strip()
        #    #if not float(tmpfile.loc[tmpfile['name']==NGSspe,'fraction_total_reads'].tolist()[0]) > 0.85:            
        #    if float(tmpfile.loc[tmpfile['name']==NGSspe,'fraction_total_reads'].tolist()[0]) < float(tspeabun):
        #        taxlist1 = proc_kra(krakenfile,2093,'G')
        #        exreadsID(Pre,taxlist1,f'raw_{Pre}.out.txt',f'{Pre}.R1.fastq.gz',f'{Pre}.R2.fastq.gz')              


    # 针对二代单端的处理
    elif not fq2 and fq1: 
        if not os.path.isfile(f'{Pre}_sub.R1.fastq'):
            if method != 'meta':
                subprocess.run(f'rasusa reads --bases 20gb -o {Pre}_sub.R1.fastq  {fq1}',shell=True)      # R1+R2=10GB
            else:
                subprocess.run(f'rasusa reads --bases 20gb -o {Pre}_sub.R1.fastq  {fq1}',shell=True)
        subprocess.run(f'seqkit stat -T {Pre}_sub.R1.fastq > raw_summary.tsv',shell=True)
        if not os.path.isfile(f'{Pre}_t.R1.fastq.gz'):
            subprocess.run(f'''fastp --in1  {Pre}_sub.R1.fastq \\
        --out1 {Pre}_t.R1.fastq.gz \\
        --thread {threads} \\
        --length_required=50 \\
        --qualified_quality_phred=10 \\
        --n_base_limit=6 \\
        --compression=6 \\
        --json {Pre}.fastp2.json \\                
        2> {Pre}.fastp2.log''',shell=True)         # 单端数据的质控相较于双端，少了检测接头--detect_adapter_for_pe，多了一步--qualified_quality_phred=10检测测序质量
        ngs_qc(Pre)
        if not os.path.isdir('R1_qc'):
            os.makedirs('R1_qc')
        with open('qc.log','a') as f:
            subprocess.run(f'ln -s  {Pre}_sub.R1.fastq raw.R1.fastq',shell=True)   # 下采样的{Pre}_sub.R1.fastq定向为raw.R1.fastq
            if rmhost == 'norm': 
                subprocess.run(f'ln -s {Pre}_t.R1.fastq.gz {Pre}.R1.fastq.gz',shell=True)
            else:
                if rnalib == '0':
                        subprocess.run(f'/home/dell/miniconda3/bin/conda run -n hostile hostile clean --fastq1 {Pre}_t.R1.fastq.gz --threads {threads}  --aligner bowtie2 --index /data/deploy/meta_new/Database/Host_Ref/hostile/human-t2t-hla.argos-bacteria-985_rs-viral-202401_ml-phage-202401 > rmcont.json',shell=True,stdout=f,stderr=f)
                        subprocess.run(f'mv  {Pre}_t.R1.clean.fastq.gz  {Pre}.R1.fastq.gz',shell=True)
                else:
                    print('测试rna建库')
                    if not os.path.isfile(f'kneaddata_out/{Pre}_t.R1_kneaddata_paired_1.fastq'):
                        subprocess.run(f'/home/dell/miniconda3/bin/conda run -n kneaddata kneaddata --bypass-trf --bypass-trim -i1 {Pre}_t.R1.fastq.gz  -o kneaddata_out -t {threads} -db /data/Ref/human_hg38_refMrna -db /data/Ref/SILVA_128_LSUParc_SSUParc_ribosomal_RNA' ,shell=True,stdout=f,stderr=f)
                    subprocess.run(f'pigz -p 10 kneaddata_out/{Pre}_t.R1_kneaddata_paired_1.fastq -c > {Pre}.R1.fastq.gz',shell=True,stdout=f,stderr=f)
            subprocess.run(f'fastqc raw.R1.fastq {Pre}.R1.fastq.gz -o R1_qc -t {threads}',shell=True,stdout=f,stderr=f)
            subprocess.run(f'''unzip -o 'R1_qc/*.zip' ''',shell=True,stdout=f,stderr=f)
            subprocess.run(f'fastp -i {Pre}.R1.fastq.gz -o tt.1.fq -w {threads} --json {Pre}.final.json',shell=True,stdout=f,stderr=f)  
            subprocess.run(f'rm tt.*.fq',shell=True)  
    summaryfastqc_prod(Pre)   # 调用函数

    open('QC_ok','w').write('已跑过')
    time.sleep(0.5)

# 定义polish_func函数：medaka对三代数据进行抛光。第一次以原始文件（{Pre}.final.fastq）作为输入，后续以上一次生成的共识序列 medaka_output{i-1}/consensus.fasta，都以-d 拼接序列（强制命名） {Pre}.consensus.fasta为参考
def polish_func(Pre, ptimes, threads, psoft='medaka'):     # Pre,pts,threads,pst(-ps,传入）
    print(f'开始抛光 抛光软件: {psoft} 抛光次数: {ptimes}')
    ptimes = str(ptimes)
    medaka_cmd = "/home/dell/miniconda3/bin/conda run -n medaka medaka_consensus"   
    try:
        with open('polish.log', 'w') as f:
            for i in range(1, int(ptimes) + 1):   # 从1到ptimes进行循环
                input_file = f"{Pre}.final.fastq" if i == 1 else f"medaka_output{i-1}/consensus.fasta"
                output_dir = f"medaka_output{i}"
                output_file = f"{output_dir}/consensus.fasta"                
                cmd = f"{medaka_cmd} -i {input_file} -d {Pre}.consensus.fasta -o {output_dir} -t {threads} > medaka.log"
                subprocess.run(cmd, shell=True, stdout=f, stderr=f)
                if i == int(ptimes): # 达到抛光次数
                    subprocess.run(f'seqkit seq -w0 {output_file} > {Pre}.polish.fasta', shell=True, stdout=f, stderr=f)       
        if os.path.getsize(f'{Pre}.polish.fasta') == 0:
            subprocess.run(f'seqkit seq -w0 {Pre}.consensus.fasta > {Pre}.polish.fasta', shell=True)  # 如果抛光后无数据，保留原始的拼接数据
    
    except Exception as e:
        print(f"抛光过程出现错误: {e}")
    
    print(f'抛光结束')
def wait_for_file(filepath,cinterval=2):
    while not os.path.exists(filepath):
        times.sleep(1)
    lastsize = 1
    while True:
        curruent_size = os.path.getsize(filepath)
        if filepath == lastsize:
            break
        lastsize =  curruent_size
        time.sleep(cinterval)
    return True


#定义rebinning函数：去冗余 MAG “标准化命名 + 重写 contig 头”
def rebinning():
    infpath = f'BASALT_out/meta_drep_out/dereplicated_genomes/'
    outpath = f'BASALT_out/meta_drep_out/binning_genomes/'
    if not os.path.isdir(outpath):
        os.makedirs(outpath)
    falist = [i for i in os.listdir(infpath) if i.endswith('.fa')]
    n=1
    open('binning_name.tsv','w').write(f'oldname\tnewname\n')
    if falist:
        for MAG in falist:
            MAG = f'{infpath}/{MAG}'
            oldname = MAG.split('/')[-1]
            oldname = re.sub(r'\.(fa|fasta|fna)(\.gz)?$', '', oldname)
            outMAG = f'{outpath}/MAG_{n}.fa'
            open(f'binning_name.tsv','a').write(f'{oldname}\tMAG_{n}\n')
            open(outMAG,'w').write('')
            with open(MAG) as f:
                m=1
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        open(outMAG,'a').write(f'>{m}\n')
                        m+=1
                    else:
                        open(outMAG,'a').write(f'{line}\n')
            n+=1

# 定义combinebin函数：把所有 MAG 的 contig 放进同一个 fasta 文件，同时给 contig 加上“来源 MAG 前缀”
def combinebin(refinedir,ofa):
    open(ofa,'w').write('')
    list1 = [i for i in os.listdir(refinedir) if i.endswith('fa')]
    for i in list1:
        filen = f'{refinedir}/{i}'
        newname = i.replace('.fa','')
        #print(newname)
        n=1
        with open(filen) as f:
            for line in f:
                if line.startswith('>'):
                    tmp_contig = line.strip().replace('>','')
                    open(ofa,'a').write(f'>{newname}_{tmp_contig}\n')
                    n+=1
                else:
                    open(ofa,'a').write(line)
# 定义bingtdbtk函数：用GTDB-Tk对MAG进行物种鉴定
def bingtdbtk_fun():
    inf = f'gtdbtk_out/gtdbtk.bac120.summary.tsv'
    with open('gtdbtk.log','w') as f:
        if not os.path.isfile(inf):
            subprocess.run(f'/home/dell/miniconda3/bin/conda run -n gtdbtk --no-capture-output gtdbtk classify_wf --genome_dir BASALT_out/meta_drep_out/dereplicated_genomes --out_dir gtdbtk_out -x .fa --cpus 10 --force',shell=True,stdout=f,stderr=f)
# 定义bincheckm2_fun函数：如果 bin_checkm2out/quality_report.tsv 不存在，就跑 CheckM2
def bincheckm2_fun():
    inf = 'bin_checkm2out/quality_report.tsv'
    with open('bincheckm2.log','w') as f:
        if not os.path.isfile(inf):
            subprocess.run(f'/home/dell/miniconda3/bin/conda run -n cm210 --no-capture-output checkm2 predict --thread 10 --input BASALT_out/meta_drep_out/binning_genomes/ --output-directory bin_checkm2out -x .fa',shell=True,stdout=f,stderr=f)
#定义binvfdrdb函数：耐药、毒力、质粒（部分）扫描
def binvfdrdb():
    inf1 = f'bin_vfdb.tsv'
    inf2 = f'bin_card.tsv'
    inf3 = 'binning_rgi_new.txt'
    inf4 = 'staramr_result/plasmidfinder.tsv'
    with open('binning.log','w') as f:
        if not os.path.isfile(inf1):
            subprocess.run(f'abricate BASALT_out/meta_drep_out/binning_genomes/MAG_*.fa --db vfdb > bin_vfdb.tsv ',shell=True,stdout=f,stderr=f)
        if not os.path.isfile(inf2):
            subprocess.run(f'abricate BASALT_out/meta_drep_out/binning_genomes/MAG_*.fa --db card > bin_card.tsv ',shell=True,stdout=f,stderr=f)
        if not os.path.isfile(inf3):
            subprocess.run(f'/home/dell/miniconda3/bin/conda run --no-capture-output  -n RGI_new rgi main -i tmp_combine.fa --clean --include_loose -o binning_rgi_new -n 20 -g PYRODIGAL --low_quality',shell=True)
            subprocess.run(f'rm -r binning_rgi_new.json',shell=True)
        if not os.path.isfile(inf4):  #plasmidfinder.tsv、resfinder.tsv
            subprocess.run(f'staramr search BASALT_out/meta_drep_out/binning_genomes/*.fa -o staramr_result -n 10',shell=True,stdout=f,stderr=f)

def join_plasmid(group):     # 定义join_plasmid函数：将每个 contig 对应的质粒预测结果（PlasFlow+plasmidfinder）进行合并
    plasmids = '|'.join(group['Plasmid'])  # 将每个contig可能包含的多个质粒预测结果（plasmid 列的值），用'|'进行分割合并
    newdb = group
    newdb['Plasmid'] = plasmids
    newdb = newdb.iloc[0,:]   # 经过上面一步后,同一个contig的多行同时包含所有预测结果，因此保留第一个即可
    return newdb   
def meta_plasmid(Pre):  # 跑PlasFlow
    with open('plasmid.log','w') as f:
        if not os.path.isfile(f'{Pre}_plaspredict.tsv') or not os.path.isfile(f'staramr_result/plasmidfinder.tsv'):
            subprocess.run(f'/home/dell/miniconda3/bin/conda run --no-capture-output -n plasflow PlasFlow.py --input tmp_combine.fa --output {Pre}_plaspredict.tsv',shell=True,stdout=f,stderr=f) 
            subprocess.run(f'staramr search BASALT_out/meta_drep_out/binning_genomes/*.fa -o staramr_result -n 10',shell=True,stdout=f,stderr=f)
        plasmiddb = pd.read_table('staramr_result/plasmidfinder.tsv')  # 读取staramr的质粒（分型）预测结果
        rawindexname = plasmiddb.index.name   # 保存 plasmiddb 索引名称
        if plasmiddb.shape[0] >1:   # plasmiddb行数＞1
            plasmiddb = plasmiddb.groupby('Contig').apply(join_plasmid)  # 按照 contig名称 进行分组，调用join_plasmid函数
            plasmiddb.index.name = rawindexname  # 将索引名称改回原来的
        plasmiddb['contig_name'] = plasmiddb.apply(lambda x:f'''{x['Isolate ID']}_{x['Contig']}''',axis=1)
        plasmiddb = plasmiddb[['contig_name','Plasmid']]
        plasflowdb = pd.read_table(f'{Pre}_plaspredict.tsv') # 读取 PlasFlow 预测结果文件
        plasflowdb = plasflowdb[['contig_name','label']]  
        plasdb = plasflowdb.merge(plasmiddb,on='contig_name',how='left').fillna('-')
        plasdb.to_csv(f'{Pre}_meta_plaspredict.tsv',sep='\t',index=False)

# CoverM：reads：2.1.fastq + 2.2.fastq；参考：binning_genomes/*.fa；指标：-m tpm；输出：meta_tpm.tsv
def meta_tpm():
    with open('coverm.log','w') as f:
        if not os.path.isfile('meta_tpm.tsv'):
            subprocess.run(f'/home/dell/miniconda3/bin/conda run -n coverm coverm genome -1 2.1.fastq -2 2.2.fastq -d  BASALT_out/meta_drep_out/binning_genomes/ -x .fa --min-read-percent-identity 95 --min-read-aligned-percent 75 -m tpm -o meta_tpm.tsv -t 10',shell=True,stderr=f,stdout=f)

def binning_result(Pre):
    #1.MAG整理 2.物种鉴定 3.checkm2 4.耐药毒力 5.质粒鉴定
    rebinning()
    combinebin('BASALT_out/meta_drep_out/binning_genomes','tmp_combine.fa')
    bingtdbtk_fun()
    bincheckm2_fun()
    binvfdrdb()
    meta_plasmid(Pre)
    meta_tpm()
    vfdb = pd.read_table('bin_vfdb.tsv')
    argdict = {'ARG':[],
    'contig_name':[],
    'Name':[],
    'AR Gene(abricate)':[],
    'AR Gene(rgi)':[],
    'AR Gene(resfinder)':[]}
    vfdb['Name'] = vfdb['#FILE'].str.split('/').str[-1].str.split('.').str[0]
    vfdb['contig_name'] = vfdb.apply(lambda x:f'''{x['Name']}_{x['SEQUENCE']}''',axis=1)
    vfdb = vfdb[['contig_name','GENE']]
    vfdb.rename(columns={'GENE':'VF Gene'},inplace=True)
    carddb = pd.read_table('bin_card.tsv')
    carddb['tmpgene'] = carddb['GENE'].str.lower()
    rgidb = pd.read_table('binning_rgi_new.txt')
    rgidb = rgidb[rgidb['Cut_Off'].isin(['Strict','Perfect'])]
    rgidb['tmpgene'] = rgidb['Best_Hit_ARO'].str.lower()
    rgidb['contig_name'] = rgidb['Contig']
    resdb = pd.read_table('staramr_result/resfinder.tsv')
    resdb['tmpgene'] = resdb['Gene'].str.lower()
    resdb['contig_name'] = resdb.apply(lambda x:f'''{x['Isolate ID']}_{x['Contig']}''',axis=1)
    carddb['Name'] = carddb['#FILE'].str.split('/').str[-1].str.split('.').str[0]
    carddb['contig_name'] = carddb.apply(lambda x:f'''{x['Name']}_{x['SEQUENCE']}''',axis=1)
    carddb = carddb[['contig_name','GENE']]
    carddb['tmpgene'] = carddb['GENE'].str.lower()
    carddb.rename(columns={'GENE':'AR Gene'},inplace=True)
    arglist = list(set(resdb['tmpgene'].tolist()+rgidb['tmpgene'].tolist()+carddb['tmpgene'].tolist()))
    for argene in arglist:
        argdict['ARG'].append(argene)
        contiglist = []
        if argene in carddb['tmpgene'].tolist():
            argdict['AR Gene(abricate)'].append('+')
            contiglist.append(carddb.loc[carddb['tmpgene']==argene]['contig_name'].tolist()[0])
        else:
            argdict['AR Gene(abricate)'].append('-')
        if argene in resdb['tmpgene'].tolist():
            argdict['AR Gene(resfinder)'].append('+')
            contiglist.append(resdb.loc[resdb['tmpgene']==argene]['contig_name'].tolist()[0])
        else:
            argdict['AR Gene(resfinder)'].append('-')
        if argene in rgidb['tmpgene'].tolist():
            argdict['AR Gene(rgi)'].append('+')
            contiglist.append(rgidb.loc[rgidb['tmpgene']==argene]['contig_name'].tolist()[0])
        else:
            argdict['AR Gene(rgi)'].append('-')
        if contiglist:
            tmpName = contiglist[0]
            argdict['contig_name'].append(tmpName)
            argdict['Name'].append('_'.join(tmpName.split('_')[:2]))
        else:
            argdict['contig_name'].append('-')
    lengths = {k: len(v) for k, v in argdict.items()}
    argdb = pd.DataFrame(argdict)
    Alldb = pd.read_table(f'{Pre}_meta_plaspredict.tsv')
    Alldb = Alldb.merge(vfdb,on='contig_name',how='left').merge(argdb,on='contig_name',how='left').fillna('-')
    Alldb['Name'] = Alldb['contig_name'].str.split('_').str[:2].str.join('_')
    gtdbdb = pd.read_table('gtdbtk_out/gtdbtk.bac120.summary.tsv')
    gtdbdb['Name'] = gtdbdb['user_genome']
    lvlist = ['D','P','C','O','F','G','S']
    for lv in lvlist:
        gtdbdb[lv] = gtdbdb['classification'].str.split(';').str[lvlist.index(lv)].str.split('__').str[1]
    gtdbdb = gtdbdb[['Name','D','P','C','O','F','G','S']]
    binnamedb = pd.read_table('binning_name.tsv')
    gtdbdb = gtdbdb.merge(binnamedb,left_on='Name',right_on='oldname')[['newname','D','P','C','O','F','G','S']].rename(columns={'newname':'Name'})
    tpmdb = pd.read_table('meta_tpm.tsv')
    tpmdb.columns = ['Name',Pre]
    Alldb = Alldb.merge(gtdbdb,on='Name',how='left').merge(tpmdb,on='Name',how='left')
    Alldb.to_csv('meta_plas_vf_card.tsv',sep='\t',index=False)

# 定义denovo_asb函数：无参组装
def denovo_asb(inf,fq1,fq2,threads,Pre,pts,pst,method,asmt,f,outputfa):
    if asmt == 'longasm':   # 组装方式：'longasm'，由-at传入              
        if method == "flye":                        
            if long_type == 'Nanopore':         # 数据类型:'Nanopore'，-lt传入
                readsq = float(os.popen(f'''nanoq -i {inf} -s -t 5 -vvv|grep 'Mean read quality'|cut -d ':' -f2''').read().strip()) # 通过 nanoq 工具对输入 FASTQ 文件进行统计，提取了平均质量值
                per=100-10**(float(readsq)/-10)*100         # 正确率（per）
                if per > 95:
                    subprocess.run(f"flye --nano-hq {inf} -g {genome_len} -o flye_output -t {threads} -i 3",shell=True,stdout=f,stderr=f)  # 正确率＞95： --nano-hq
                else:
                    subprocess.run(f"flye --nano-raw {inf} -g {genome_len} -o flye_output -t {threads} -i 3",shell=True,stdout=f,stderr=f) # 正确率≤95：--nano-raw
            elif long_type == 'PacBio_CLR':     # PacBio_CLR数据类型，--pacbio-raw
                subprocess.run(f"flye --pacbio-raw {inf} -g {genome_len} -o flye_output -t {threads} -i 3",shell=True,stdout=f,stderr=f)   
            elif long_type == 'PacBio_CCS':     # PacBio_CCS数据类型，--pacbio-hifi
                subprocess.run(f"flye --pacbio-hifi {inf} -g {genome_len} -o flye_output -t {threads} -i 3",shell=True,stdout=f,stderr=f)  
            subprocess.run(f"cp flye_output/assembly.fasta {outputfa}",shell=True)
            flyedb = pd.read_table(f'flye_output/assembly_info.txt',names=['seq_name','length','cov','circ','repeat','mult','alt_group','graph_path'],skiprows=1) # skiprows=1，跳过第一行
            flyedb.rename(columns={'seq_name':'序列名称','length':'序列长度','cov':'平均深度','circ':'是否成环'},inplace=True)
            flyedb[['序列名称','序列长度','平均深度','是否成环']].to_csv(f'flye_output/assembly_info.txt',sep='\t',index=False)  # 从 flyedb 数据框中选择特定的列 保存到 assembly_info.txt                      
        # 其他三代拼接软件canu,wtdbg2,unicycler,miniasm
        elif method == 'miniasm':      # miniasm拼接
            if long_type == 'Nanopore':
                subprocess.run(f"minimap2 -x ava-ont -t {threads} {inf} {inf} | gzip -1 > {Pre}_reads.paf.gz",shell=True,stdout=f,stderr=f)
            elif long_type == 'PacBio_CLR':
                subprocess.run(f"minimap2 -x map-pb -t {threads} {inf} {inf} | gzip -1 > {Pre}_reads.paf.gz",shell=True,stdout=f,stderr=f)
            elif long_type == 'PacBio_CCS':
                subprocess.run(f"minimap2 -x map-hifi -t {threads} {inf} {inf} | gzip -1 > {Pre}_reads.paf.gz",shell=True,stdout=f,stderr=f)
            subprocess.run(f"miniasm -f {inf} {Pre}_reads.paf.gz > {Pre}_reads.gfa",shell=True,stdout=f,stderr=f)
            subprocess.run(f"gfatools gfa2fa {Pre}_reads.gfa > {Pre}_miniasm.fa",shell=True,stdout=f,stderr=f)
            subprocess.run(f"cp {Pre}_miniasm.fa {outputfa}",shell=True)
        elif method == "wtdbg2":       # wtdbg2拼接
            subprocess.run(f"perl ~/biosoft/wtdbg2/wtdbg2.pl -t {threads} -x ont -g {genome_len} -o wtdbg2 {inf}",shell=True,stdout=f,stderr=f)
            subprocess.run(f"cp wtdbg2.cns.fa {outputfa}",shell=True)
        elif method == "canu":         # canu拼接
            if long_type == 'Nanopore':
                subprocess.run(f"time canu -d canu -p canu genomeSize={genome_len} maxThreads={threads} -nanopore-raw {inf} >canu.log",shell=True,stdout=f,stderr=f)
            else:
                subprocess.run(f"time canu -d canu -p canu genomeSize={genome_len} maxThreads={threads} -pacbio-raw {inf} >canu.log",shell=True,stdout=f,stderr=f)
            subprocess.run(f"cp canu/canu.contigs.fasta {outputfa}",shell=True)                        
        elif method == "unicycler":     # unicyler拼接
            subprocess.run(f"unicycler -t {threads} -l {inf} -o unicycler",shell=True,stdout=f,stderr=f)
            subprocess.run(f"cp unicycler/assembly.fasta {outputfa}",shell=True)                        
        elif method == 'raven':   # raven拼接
            subprocess.run(f"raven {inf} -t {threads} > raven.fasta",shell=True,stdout=f,stderr=f)
            subprocess.run(f"cp raven.fasta {outputfa}",shell=True)                          
        else:
            print(f'请确认传入参数method是否正确，可选[flye,canu,wtdbg2,unicycler,miniasm],您输入的为：{method}')                   
        if os.path.isfile(outputfa):                # 成功获得三代拼接序列
            polish_func(Pre,pts,threads,pst)        # 调用polish_func抛光函数
            
    # 二代数据拼接，spades、masurca                      
    elif asmt == 'shortasm':       # 由-at传入  
        if method == 'spades':
            if fq2:
                subprocess.run(f"spades.py --pe1-1 {fq1} --pe1-2 {fq2} -t {threads} -o spades_output --cov-cutoff 8 --isolate",shell=True,stdout=f,stderr=f)  # --isolate：单一菌株的基因组
            else:
                subprocess.run(f"spades.py -s {fq1} -t {threads} -o spades_output --isolate --cov-cutoff 8",shell=True,stdout=f,stderr=f)
            if os.path.isfile('spades_output/contigs.fasta'):
                subprocess.run(f'cp spades_output/contigs.fasta {outputfa}',shell=True,stdout=f,stderr=f)                           
        
        elif method == 'masurca':
            if fq2:
                subprocess.run(f'masurca -i {fq1},{fq2} -t {threads} -o masurca_output',shell=True,stdout=f,stderr=f)
            else:
                subprocess.run(f'masurca -i {fq1} -t {threads} -o masurca_output',shell=True,stdout=f,stderr=f)
            CAdir = os.popen('ls -d CA').read().split('\n')[0]
            if os.path.isfile(f'{CAdir}/primary.genome.scf.fasta'):
                subprocess.run(f'cp {CAdir}/primary.genome.scf.fasta {outputfa}',shell=True)
            else:
                if os.path.isfile(f'{CAdir}/scaffolds.ref.fa'):
                    subprocess.run(f'cp {CAdir}/scaffolds.ref.fa {outputfa}',shell=True)
                else:
                    print('组装失败')
                    sys.exit()
        elif method == 'meta':  # megahit宏组装
            if not os.path.isfile(f'megahit_output/final.contigs.fa'):
                if fq2:
                    subprocess.run(f'/home/dell/miniconda3/bin/conda run -n BASALT megahit -1 {fq1} -2 {fq2} -o megahit_output',shell=True,stdout=f,stderr=f)
                else:
                    subprocess.run(f'/home/dell/miniconda3/bin/conda run -n BASALT megahit -1 {fq1} -o megahit_output',shell=True,stdout=f,stderr=f)
            if os.path.isfile(f'megahit_output/final.contigs.fa') and os.path.getsize('megahit_output/final.contigs.fa') != 0:
                #---运行时间1-2h
                tmpwkdir = os.getcwd()
                #---判断drep文件
                if not os.path.isdir(f'BASALT_out'):
                    os.makedirs(f'BASALT_out')
                #os.chdir('BASALT')
                subprocess.run(f'ln -s {tmpwkdir}/{fq1} ./BASALT_out',shell=True)
                if fq2:
                    subprocess.run(f'ln -s {tmpwkdir}/{fq2} ./BASALT_out',shell=True)
                subprocess.run(f'ln -s {tmpwkdir}/megahit_output/final.contigs.fa ./BASALT_out',shell=True)
                os.chdir('BASALT_out')
                if os.path.isdir('BestBinset_outlier_refined'):
                    if not [i for i in os.listdir('BestBinset_outlier_refined')]:
                        if fq2:     #--module autobinning自动分箱；--module refinement优化
                            subprocess.run(f'/home/dell/miniconda3/bin/conda run -n BASALT  BASALT -a final.contigs.fa -s {fq1},{fq2} -t {threads} --module autobinning -q checkm2',shell=True,stdout=f,stderr=f)
                            subprocess.run(f'/home/dell/miniconda3/bin/conda run -n BASALT  BASALT -a final.contigs.fa -s {fq1},{fq2} -t {threads} --module refinement -q checkm2',shell=True,stdout=f,stderr=f)
                        else:
                            subprocess.run(f'/home/dell/miniconda3/bin/conda run -n BASALT  BASALT -a final.contigs.fa -s {fq1} -t {threads} --module autobinning -q checkm2',shell=True,stdout=f,stderr=f)
                            subprocess.run(f'/home/dell/miniconda3/bin/conda run -n BASALT  BASALT -a final.contigs.fa -s {fq1} -t {threads} --module refinement -q checkm2',shell=True,stdout=f,stderr=f)
                else:
                        if fq2:
                            subprocess.run(f'/home/dell/miniconda3/bin/conda run -n BASALT  BASALT -a final.contigs.fa -s {fq1},{fq2} -t {threads} --module autobinning -q checkm2',shell=True,stdout=f,stderr=f)
                            subprocess.run(f'/home/dell/miniconda3/bin/conda run -n BASALT  BASALT -a final.contigs.fa -s {fq1},{fq2} -t {threads} --module refinement -q checkm2',shell=True,stdout=f,stderr=f)
                        else:
                            subprocess.run(f'/home/dell/miniconda3/bin/conda run -n BASALT  BASALT -a final.contigs.fa -s {fq1} -t {threads} --module autobinning -q checkm2',shell=True,stdout=f,stderr=f)
                            subprocess.run(f'/home/dell/miniconda3/bin/conda run -n BASALT  BASALT -a final.contigs.fa -s {fq1} -t {threads} --module refinement -q checkm2',shell=True,stdout=f,stderr=f)

                #---检测文件，合并fasta
                if os.path.isdir('BestBinset_outlier_refined'):
                    binlist = [i for i in os.listdir('BestBinset_outlier_refined') if i.endswith('fa') or i.endswith('fasta')] 
                    if len(binlist)>1:
                        bincheckmdb = pd.read_table('BestBinset_outlier_refined/quality_report.tsv',header=None) 
                        #dRep--去重，--ignoreGenomeQuality不考虑质量
                        subprocess.run(f'/home/dell/miniconda3/bin/conda run -n BASALT dRep dereplicate meta_drep_out -g BestBinset_outlier_refined/*.fa* --ignoreGenomeQuality -p {threads}',shell=True,stdout=f,stderr=f)
                        if os.path.isdir('meta_drep_out/dereplicated_genomes'):
                            dreplist = [i for i in os.listdir('meta_drep_out/dereplicated_genomes') if i.endswith('fasta') or i.endswith('fa')]
                            if not dreplist:
                                subprocess.run(f'cp BestBinset_outlier_refined/*.fa* meta_drep_out/dereplicated_genomes/',shell=True)
                        else:
                            os.makedirs('meta_drep_out/dereplicated_genomes')
                            subprocess.run(f'cp BestBinset_outlier_refined/*.fa* meta_drep_out/dereplicated_genomes/',shell=True)

                    else:
                        if not os.path.isdir(f'meta_drep_out/dereplicated_genomes'):
                            os.makedirs('meta_drep_out/dereplicated_genomes')
                        subprocess.run(f'cp BestBinset_outlier_refined/*.fa* meta_drep_out/dereplicated_genomes/tmp.fa',shell=True)

                    #combinebin('meta_drep_out/dereplicated_genomes',outputfa)
                    #subprocess.run(f'cp {outputfa} {tmpwkdir}',shell=True)
                else:
                    print('refine binning derep失败')
                    subprocess.run(f'cp {tmpwkdir}/megahit_output/final.contigs.fa {tmpwkdir}/{outputfa}',shell=True)
                os.chdir(tmpwkdir)
                binning_result(Pre)

            else:
                print('宏基因组组装失败')
                sys.exit()  #宏组装的合并拼接文件为consensus.fasta文件，过滤短片段后的最终文件为final.fasta
                    
# 二代数据和三代数据混合拼接：unicycler、masurca
    elif asmt == 'shortlongasm':    # 由-at传入 
        if method == 'unicycler':   
            if fq2:
                subprocess.run(f'unicycler -1 {fq1} -2 {fq2} -l {inf} -t {threads} -o unicycler',shell=True,stdout=f,stderr=f)  # -l，指定长读长数据
            elif fq1:
                subprocess.run(f'unicycler -1 {fq1} -l {inf} -t {threads} -o unicycler',shell=True,stdout=f,stderr=f)
            if os.path.isfile(f'unicycler/assembly.fasta'):
                subprocess.run(f"cp unicycler/assembly.fasta {outputfa}",shell=True)
        
        elif method == 'masurca':
            if fq2:
                subprocess.run(f'masurca -i {fq1},{fq2} -t {threads} -o masurca_output -r {inf}',shell=True,stdout=f,stderr=f)  # -r，指定长读长数据
                CAdir = os.popen('ls -d CA.*').read().split('\n')[0]                       
                if os.path.isfile(f'{CAdir}/primary.genome.scf.fasta'):
                    subprocess.run(f'cp {CAdir}/primary.genome.scf.fasta {outputfa}',shell=True)
                elif os.path.isfile(f'{CAdir}/scaffolds.ref.fa'):
                    subprocess.run(f'cp {CAdir}/scaffolds.ref.fa {outputfa}',shell=True)
            else:
                subprocess.run(f'masurca -i {fq1} -t {threads} -o masurca_output -r {inf}',shell=True,stdout=f,stderr=f)
                CAdir = os.popen('ls -d CA.*').read().split('\n')[0]
                if os.path.isfile(f'{CAdir}/primary.genome.scf.fasta'):
                    subprocess.run(f'cp {CAdir}/primary.genome.scf.fasta {outputfa}',shell=True)
                elif os.path.isfile(f'{CAdir}/scaffolds.ref.fa'):
                    subprocess.run(f'cp {CAdir}/scaffolds.ref.fa {outputfa}',shell=True)

# 定义reassm_fun函数：有参组装
def reassm_fun(inf,fq1,fq2,threads,Pre,pts,pst,method,asmt,f,outputfa):
    if not os.path.isdir('genomes'):
        os.makedirs('genomes')
    subprocess.run(f'seqkit seq {ref} > genomes/ref.fa',shell=True)  # 使用 seqkit 将 {ref} 文件格式化输出到 genomes/ref.fa
    subprocess.run(f'samtools faidx genomes/ref.fa ',shell=True)    # 使用samtools为genomes/ref.fa 创建.fai索引文件               
    # 是否进行额外注释
    if gtf == 'nogtf':
        ifAvcf = 0
        print('snp位点不进行额外注释')
        sys.stdout.flush()
    else:  # 进行注释    
        if not os.path.isdir('ref'):
            os.makedirs('ref')
        subprocess.run(f'cp {gtf} ref/genes.gff',shell=True)  # 将传入的 GFF 文件复制到 ref/genes.gff
        open('snpEff.config','w').write('ref.genome : ref')   # 创建并写入 snpEff 配置文件
        subprocess.run(f'java -jar /home/dell/biosoft/snpEff/snpEff.jar build -gff3 ref -c snpEff.config -dataDir ./',shell=True)  # 使用snpEff工具构建GFF3格式的注释文件
        if not os.path.isfile('ref/snpEffectPredictor.bin'):  # 检查snpEff构建注释文件是否成功
            print('gff与fa不匹配，snp位点不进行额外注释') # 即，{ref}与{gtf}不符
            ifAvcf = 0
            sys.stdout.flush()
        else:
            print('注释文件正常，snp位点根据注释文件进行注释')            
            ifAvcf = 1
            sys.stdout.flush()
            
    # 长序列使用minimap2进行比对，短序列使用bwa进行比对。最后生成ref.mapping.bam文件
    if 'long' in asmt:
        subprocess.run(f'minimap2 -ax map-ont genomes/ref.fa {inf} -t {threads} |samtools sort -o ref.mapping.bam',shell=True) # -ax map-ont 适用于 Oxford Nanopore比对    
    elif 'short' in asmt:
        with open('mapping.log','w') as mapplog:
            subprocess.run(f'bwa index genomes/ref.fa',shell=True)  #构建索引
            if fq2:
                subprocess.run(f'bwa mem genomes/ref.fa {fq1} {fq2} -t {threads} |samtools sort -o ref.mapping.bam',shell=True,stdout=mapplog,stderr=mapplog)
            else:
                subprocess.run(f'bwa mem genomes/ref.fa {fq1}  -t {threads} |samtools sort -o ref.mapping.bam',shell=True,stdout=mapplog,stderr=mapplog)                
    # 使用mosdepth工具，分别使用步长1000bp和1bp计算覆盖度              
    subprocess.run(f'samtools index ref.mapping.bam',shell=True)  # 使用 samtools index 为 BAM 文件生成索引
    subprocess.run(f'mosdepth -b 1000 ref_map ref.mapping.bam -t {threads}',shell=True) # 使用 mosdepth 计算在每个 1000 bp 的区间内的覆盖度
    subprocess.run(f'gunzip ref_map.regions.bed.gz',shell=True)                         # gunzip解压         
    subprocess.run(f'mosdepth -b 1 ref_map1 ref.mapping.bam -t {threads}',shell=True)   # 使用 mosdepth 计算以每个碱基为单位的覆盖度
    subprocess.run(f'gunzip ref_map1.regions.bed.gz',shell=True)                        # gunzip解压              
    # 定义outfun内函数：基因覆盖深度统计
    def outfun(x):
        tmpdict = {}  # 建立一个空字典
        ofname = x['GeneName'].tolist()[0]
        x['start'] = x.reset_index().index+1
        x['end'] = x.reset_index().index+2
        x[['Chrom','start','end','Depth']].to_csv(f'geneDepth/{ofname}.tsv',sep='\t',header=False,index=False)
        tmpdict['片段名称'] = x['Chrom'].tolist()[0]
        tmpdict['起始位置'] = x['start'].min()
        tmpdict['终止位置'] = x['start'].max()
        tmpdict['覆盖度(>0)%'] =  round(x[x['Depth']>0].shape[0]/x.shape[0],4)*100 # 覆盖度大于0%的百分比
        tmpdict['覆盖度(>10)%'] =  round(x[x['Depth']>10].shape[0]/x.shape[0],4)*100  # 覆盖度大于10%的覆盖度
        tmpdict['覆盖度(>100)%'] =  round(x[x['Depth']>100].shape[0]/x.shape[0],4)*100 # 覆盖度大于100%的覆盖度
        tmpdict['平均深度'] = round(x['Depth'].mean(),2)
        tmpdict['最低深度'] = x['Depth'].min()
        tmpdict['最高深度'] = x['Depth'].max()
        return pd.DataFrame(tmpdict,index=[0]).round(2)                      
    # 读取 GFF 文件，提取每个基因的 序列ID、起始位置、终止位置、基因名，每个基因生成单独的bed文件               
    if os.path.isfile('ref/genes.gff'):        # 成功接收传入.gff文件
        open('geneNamelist.txt','w').write('')
        if not os.path.isdir('geneDepth'):
            os.makedirs('geneDepth')
        with open('ref/genes.gff') as f:
            for line in f:
                if not line.startswith('#'):   # 不选入开头几行注释
                    line = line.strip().split('\t')  # 将每一行数据按tab分割成列表
                    if line[2] == 'gene':      # 选出第3列是'gene'对应的行
                        if 'gene=' in line[8]: # 该行第9个字段包含'gene='
                            gName = line[8].split('gene=')[1].split(';')[0].split('/')[0] #先按 gene= 切分，取第二部分；再将结果按 ';'切分，取第一部分；最后再按'/'切分，取第一部分，获得基因名称
                            open(f'geneDepth/{gName}.bed','w').write(f'{line[0]}\t{line[3]}\t{line[4]}\t{gName}\n')  # 该基因的 序列ID（line[0]）、起始位置（line[3]）、终止位置（line[4]）及基因名写入对应的 .bed 文件
                            open('geneNamelist.txt','a').write(f'{gName}\n')   # 以追加模式将 'gName' 写入geneNamelist.txt
                        else:
                            if 'ID=' in line[8]: # 该行第9个字段包含'ID='
                                gName = line[8].split('ID=')[1].split(';')[0]  #先按 gene= 切分，取第二部分；再将结果按 “;” 切分，取第一部分
                                open(f'geneDepth/{gName}.bed','w').write(f'{line[0]}\t{line[3]}\t{line[4]}\t{gName}\n') 
                                open('geneNamelist.txt','a').write(f'{gName}\n')                 
        subprocess.run(f'cat geneDepth/*.bed > All_gene.bed',shell=True)  # 合并多个BED 文件，生成All_gene.bed文件
        subprocess.run(f'bedtools intersect -a All_gene.bed -b ref_map1.regions.bed -wb > ref_map1.Anno.regions.bed',shell=True) # 使用bedtools取All_gene.bed和ref_map1.regions.bed文件的交集
        # 调用outfun基因深度统计函数汇总基因覆盖度等信息
        tmpd = pd.read_table('ref_map1.Anno.regions.bed',header=None,names=['Chrom','start','end','GeneName','c1','s1','e1','Depth']) # 为合并后的文件指定列的名称
        if tmpd.shape[0] != 0:
            newtd = tmpd.groupby('GeneName').apply(outfun).reset_index(level=0)   # groupby（）：相同基因名称 (GeneName) 的行分为一组。应用 outfun 函数，返回每个基因的统计信息
            newtd.rename(columns={'GeneName':'基因名称'},inplace=True)             # 将列'GeneName'重命名为'基因名称'
            newtd.to_csv('gene_summary.tsv',sep='\t',index=False)
        else:
            print('gff与基因组不匹配，无法展示各个基因区段覆盖度')
            sys.stdout.flush()
    subprocess.run(f'python /data1/shanghai_pip/meta_genome/soft/IGV_js/IGV_new.py -r genomes/ref.fa -m ref.mapping.bam -o ./ -s {Pre}',shell=True) #调用IGV_new.py脚本，1000bp步长的比对文件：ref.mapping.bam，参考基因组：genomes/ref.fa
    subprocess.run(f'cp /data1/shanghai_pip/meta_genome/soft/IGV_js/igv.min.js ./',shell=True)  # 将 IGV 在浏览器中可视化显示基因组所需的 JavaScript 文件 igv.min.js 复制到当前目录
    subprocess.run(f'mosdepth -b 1 -t 10 ref ref.mapping.bam',shell=True)                         # 使用 mosdepth 计算每个碱基的覆盖度，输出ref.regions.bed.gz
    subprocess.run(f'gunzip ref.regions.bed.gz',shell=True)
    refbeddb = pd.read_table('ref.regions.bed',header=None)    
    refbeddb[refbeddb[3]==0].to_csv('mask.bed',index=False,header=False,sep='\t')  #  将覆盖度为 0 的区域保存为 mask.bed 文件              
    
    # 针对二代序列，freebayes变异检测工具
    if method == 'freebayes':
        subprocess.run(f'fasta_generate_regions.py genomes/ref.fa.fai 200000 > ref.txt',shell=True)   # 调用外脚本fasta_generate_regions.py,将基因组索引文件划分为多个200kb大小的的文件，方便将任务分成多个并行子任务
        # 使用 freebayes 对 ref.mapping.bam 进行变异检测，生成snps.raw.vcf | 使用 bcftools 进行过滤和删除不必要的字段，生成 snps.filt1.vcf
        subprocess.run(f'freebayes-parallel ref.txt {threads} -p 2 -P 0 -C 2 -F 0.05 --min-coverage 10 --min-repeat-entropy 1.0 -q 30 -m 30 --strict-vcf -f genomes/ref.fa ref.mapping.bam > snps.raw.vcf',shell=True)
        subprocess.run(f'''bcftools view --include 'QUAL>=20 && FMT/DP>=10 && (FMT/AO)/(FMT/DP)>=0.9' snps.raw.vcf  | bcftools annotate --remove '^INFO/TYPE,^INFO/DP,^INFO/RO,^INFO/AO,^INFO/AB,^FORMAT/GT,^FORMAT/DP,^FORMAT/RO,^FORMAT/AO,^FORMAT/QR,^FORMAT/QA,^FORMAT/GL' > snps.filt1.vcf''',shell=True)
        if ifAvcf:  # 进行注释
            subprocess.run(f'''java -jar /home/dell/biosoft/snpEff/snpEff.jar ann  -c snpEff.config -Datadir . ref snps.filt1.vcf > snps.anno.vcf''',shell=True) # 使用 snpEff 工具对 ref snps.filt1.vcf 进行注释，生成 snps.anno.vcf
            skinums = int(os.popen('''grep '##' snps.anno.vcf |wc -l''').read())
            annodb = pd.read_table('snps.anno.vcf',skiprows=skinums)  # annodb = 跳过了前面以'##'开头行的 snps.anno.vcf 文件
            annodb['参考碱基'] = annodb['REF']+"("+annodb['unknown'].str.split('|').str[-1].str.strip().str.split(':').str[-4] + ")"
            annodb['突变碱基'] = annodb['ALT']+"("+annodb['unknown'].str.split('|').str[-1].str.strip().str.split(':').str[-3] + ")"
            annodb['突变类型'] = annodb['INFO'].str.split('|').str[1]
            annodb['突变影响'] = annodb['INFO'].str.split('|').str[2]
            annodb['影响基因'] = annodb['INFO'].str.split('|').str[3]
            annodb['碱基变化'] = annodb['INFO'].str.split('|').str[9] 
            annodb['氨基酸变化'] = annodb['INFO'].str.split('|').str[10]
            annodb['突变位置'] = annodb['POS']
            annodb['片段名称'] = annodb['#CHROM']
            annodb = annodb[['片段名称','突变位置','参考碱基','突变碱基','影响基因','突变类型','突变影响','碱基变化','氨基酸变化']]
            annodb['氨基酸变化'] = annodb.apply(lambda x:'-' if x['氨基酸变化'] == '' else x['氨基酸变化'],axis=1)
            annodb.fillna('-',inplace=True) # 将所有的 NaN（缺失值）替换为 '-'
            annodb.to_csv(f'{Pre}.anno.tsv',sep='\t',index=False) # 最后生成{Pre}.anno.tsv表格
        else:  # 不进行注释
            subprocess.run(f'ln -s snps.filt1.vcf snps.anno.vcf',shell=True)
            skinums = int(os.popen('''grep '##' snps.anno.vcf |wc -l''').read())
            annodb = pd.read_table('snps.anno.vcf',skiprows=skinums)
            annodb['参考碱基'] = annodb['REF']+"("+annodb['unknown'].str.split('|').str[-1].str.strip().str.split(':').str[-5] + ")"
            annodb['突变碱基'] = annodb['ALT']+"("+annodb['unknown'].str.split('|').str[-1].str.strip().str.split(':').str[-3] + ")"
            annodb['突变位置'] = annodb['POS']
            annodb['片段名称'] = annodb['#CHROM']
            annodb = annodb[['片段名称','突变位置','参考碱基','突变碱基']]
            annodb.to_csv(f'{Pre}.anno.tsv',sep='\t',index=False)
    
    # 针对三代序列，clair3变异检测工具
    elif method == 'clair3':
        with open('clair3.log','w') as clg:   # 目前的--platform 和 --model_path参数仅适用于Nanopore
            subprocess.run(f'conda run -n clair3 --no-capture-output run_clair3.sh --bam_fn=ref.mapping.bam --ref_fn=genomes/ref.fa --threads={threads} --platform=ont --model_path=/home/dell/miniconda3/envs/clair3/bin/models/r941_prom_sup_g5014 --output=./ --include_all_ctgs --enable_long_indel --snp_min_af=0.05',shell=True,stdout=clg,stderr=clg)
        subprocess.run(f'samtools view -h -F 2308 ref.mapping.bam |samtools sort -o ref.filter.bam',shell=True)   # 使用 samtools 过滤并排序 ref.mapping.bam 文件,输出 ref.filter.bam 
        subprocess.run(f'samtools index ref.filter.bam',shell=True)       #使用 samtools index 创建 ref.filter.bam 文件的索引
        subprocess.run(f'perbase base-depth ref.filter.bam -t {threads} > perbase.bed',shell=True)  # 使用 perbase 工具计算 ref.filter.bam 中每个碱基的覆盖深度，结果保存到 perbase.bed 文件中
        perbasedb = pd.read_table('perbase.bed')
        if ifAvcf:
            subprocess.run(f'''java -jar /home/dell/biosoft/snpEff/snpEff.jar ann  -c snpEff.config -Datadir . ref merge_output.vcf.gz > snps.anno.vcf''',shell=True)
            skinums = int(os.popen('''grep '##' snps.anno.vcf |wc -l''').read())
            annodb = pd.read_table('snps.anno.vcf',skiprows=skinums)
            annodb['参考碱基'] = annodb['REF']+"("+annodb['SAMPLE'].str.split('|').str[-1].str.strip().str.split(':').str[-2].str.split(',').str[0]+")"
            annodb['突变碱基'] = annodb['ALT']+"("+annodb['SAMPLE'].str.split('|').str[-1].str.strip().str.split(':').str[-2].str.split(',').str[1]+")"
            annodb['突变类型'] = annodb['INFO'].str.split('|').str[1]
            annodb['突变影响'] = annodb['INFO'].str.split('|').str[2]
            annodb['影响基因'] = annodb['INFO'].str.split('|').str[3]
            annodb['碱基变化'] = annodb['INFO'].str.split('|').str[9] 
            annodb['氨基酸变化'] = annodb['INFO'].str.split('|').str[10]
            annodb['突变位置'] = annodb['POS']
            annodb['片段名称'] = annodb['#CHROM']
            annodb = annodb[['片段名称','突变位置','参考碱基','突变碱基','影响基因','突变类型','突变影响','碱基变化','氨基酸变化']]
            annodb['氨基酸变化'] = annodb.apply(lambda x:'-' if x['氨基酸变化'] == '' else x['氨基酸变化'],axis=1)
            annodb.to_csv(f'{Pre}.anno.tsv',sep='\t',index=False)
            #subprocess.run(f'per',shell=True)
        else:
            subprocess.run(f'gunzip merge_output.vcf.gz -c > snps.anno.vcf',shell=True)
            skinums = int(os.popen('''grep '##' snps.anno.vcf |wc -l''').read())
            annodb = pd.read_table('snps.anno.vcf',skiprows=skinums)
            annodb['参考碱基'] = annodb['REF']+"("+annodb['SAMPLE'].str.split('|').str[-1].str.strip().str.split(':').str[-2].str.split(',').str[0]+")"
            annodb['突变碱基'] = annodb['ALT']+"("+annodb['SAMPLE'].str.split('|').str[-1].str.strip().str.split(':').str[-2].str.split(',').str[1]+")"
            annodb['突变位置'] = annodb['POS']
            annodb['片段名称'] = annodb['#CHROM']
            annodb = annodb[['片段名称','突变位置','参考碱基','突变碱基']]
            annodb.to_csv(f'{Pre}.anno.tsv',sep='\t',index=False)
    # --将测序深度低于10x的区域标记为N
    subprocess.run(f'bcftools convert -Oz -o snps.anno.vcf.gz snps.anno.vcf',shell=True)  # bcftools convert 压缩文件
    subprocess.run(f'bcftools index -f snps.anno.vcf.gz',shell=True)                      # 创建索引
    if outputfa !='noforce':
        subprocess.run(f'bcftools consensus -f genomes/ref.fa -o {outputfa} snps.anno.vcf.gz -m mask.bed',shell=True)     # bcftools consensus，生成共识序列

# 定义renamefa函数：将输入的fasta文件序列从长到短排列，保留＞1000bp的片段并统一命名
def renamefa(inf,ofn):  # (polifa/outputfa/finalfa,finalfa) 
    #1.如果长度没有>1000则不过滤
    subprocess.run(f'seqkit sort -l -r {inf}|seqkit fx2tab > tmpfa.tab',shell=True) # seqkit对输入的fasta文件按序列长度递减排序 | 转换为制表符分隔的文本格式
    afile = pd.read_table('tmpfa.tab',header=None)  
    afile['contignum'] = afile.index+1           # 添加一列'contignum'=行索引+1，即加一列行数
    if afile.loc[afile[1].str.len()>1000,].shape[0] > 0:
        afile = afile.loc[afile[1].str.len()>1000,]  # 选择第二列字符串长度大于1000的行
    afile[3] = afile.apply(lambda x: f'contig_{x.contignum}',axis=1)  # 添加第四列，每行命名格式：contig_contignum
    afile[[0,3]].to_csv('transname.tsv',sep='\t',index=False,header=False)   # 原始序列ID（第一列）+ 新生成的序列ID（第四列） = transname.tsv
    afile[0] = afile.apply(lambda x: f'>contig_{x.contignum}',axis=1)        # 将第一列重命名为>contig_contignum,即fasta的标准序列ID格式  
    afile = afile[[0,1]]      # 保留第一列和第二列
    afile.to_csv(ofn,sep='\n',index=False,header=False)   # finalfa = 改换contig名称的fasta的标准格式

# 定义compareid1函数：
def compareid1(level1,level2,rawlist):
    rlevel1 = [i for i in ['R','D','K','P','C','O','F','G','S'] if i in level1][0]
    rlevel2 = [i for i in ['R','D','K','P','C','O','F','G','S'] if i in level2][0]

    if rawlist.index(rlevel1) == rawlist.index(rlevel2):
        if level1 > level2:
            return 1 
        elif level1 == level2:
            return 0
        else:
            return -1
    elif rawlist.index(rlevel1) > rawlist.index(rlevel2):
        return 1
    else:
        return -1
# 定义proc_kra1函数：  
def proc_kra1(kraken,tax,lel='S'):
    tmplist = [tax]
    if [i for i in ['R','D','K','P','C','O','F','G','S'] if i in lel]:    
        rawlist = ['R','D','K','P','C','O','F','G','S']
    else:
        rawlist = ['S1','S2','S3','S4','S5','S6']
    if tax != 0:
        kradb = pd.read_table(kraken,header=None)
        #kradb = kradb[(kradb[3]==lel)&(kradb[4]==tax)]
        kradb[4] = kradb[4].astype('str')
        tmpindex = kradb[(kradb[3]==lel)&(kradb[4]==str(tax))].index.tolist()[0]+1
        if tmpindex <=  kradb.shape[0]-1:
            def getlindex(tmpindex):
                tmpl = kradb.iloc[tmpindex,3]
                for tl in rawlist:
                    if tl in tmpl:
                        if tl == tmpl:
                            tmpl = tl
                            tmlindex = rawlist.index(tmpl)
                        else:
                            tmpl = tl
                            tmlindex = rawlist.index(tmpl)+1
                        return tmlindex
            while compareid1(kradb.iloc[tmpindex,3],lel,rawlist) == 1 and tmpindex <= kradb.shape[0]-2:
                tmplist.append(kradb.iloc[tmpindex,4])
                tmpindex+=1
    return tmplist
# 定义exreadsID1函数： 
def exreadsID1(taxlist,kraresult,fq1,fq2=0):
    Maintax = taxlist[0]
    #kraredb = pd.read_table(kraresult,header=None)
    kraredb = pd.read_csv(kraresult, header=None, usecols=[1, 2], dtype={1:'str',2:'int32'},sep='\t')
    tmp2db = kraredb[kraredb[2].isin(taxlist)]
    tmp1db = pd.DataFrame(tmp2db[1].unique())
    tmp2db.to_csv(f'{Maintax}.id.tsv',sep='\t',index=False)
    #tmp1db.to_csv(f'{Maintax}_fqID.txt',index=False,header=False)
    pd.DataFrame(tmp1db).to_csv(f'{Maintax}_fqID.txt', index=False, header=False)
    subprocess.run(f'head -n 1 {Maintax}_fqID.txt > tt.txt',shell=True)
    subprocess.run(f'''cut -d '/' -f1 {Maintax}_fqID.txt|sort -u > {Maintax}.listID.txt''',shell=True)
    if os.popen(f'''head -n 1 tt.txt''').read().strip().endswith('/1') or os.popen(f'''head -n 1 tt.txt''').read().strip().endswith('/2'):
        subprocess.run(f'''sed 's/$/\/1/' {Maintax}.listID.txt > {Maintax}.listID1.txt''',shell=True)
        subprocess.run(f'''sed 's/$/\/2/' {Maintax}.listID.txt > {Maintax}.listID2.txt''',shell=True)
        subprocess.run(f'seqkit grep -f {Maintax}.listID1.txt {fq1} > {Maintax}.1.fastq',shell=True)
        if fq2:
            subprocess.run(f'seqkit grep -f {Maintax}.listID2.txt {fq2} > {Maintax}.2.fastq',shell=True)
    else:
        subprocess.run(f'seqkit grep -f {Maintax}.listID.txt {fq1} > {Maintax}.1.fastq',shell=True)
        if fq2:
            subprocess.run(f'seqkit grep -f {Maintax}.listID.txt {fq2} > {Maintax}.2.fastq',shell=True)

# 定义asb_func函数：组装
def asb_func(inf,fq1,fq2,threads,Pre,lelID,pts,pst,method,asmt='longasm',ref='noref',gtf='nogtf',tryref=False):
    with open('asb.log','w') as f:
        outputfa = f'{Pre}.consensus.fasta'   # 拼接文件
        polifa = f'{Pre}.polish.fasta'        # 最终抛光文件 
        with open('tmpkk2.log','w') as kkf:
            if not os.path.isfile('2.1.fastq'):
                if inf:         
                    subprocess.run(f'kraken2 --db {Krdb} --threads {threads} --output {Pre}.out.txt --report {Pre}.report.txt {inf}',shell=True,stdout=kkf,stderr=kkf)
                    subprocess.run(f'bracken -d {Krdb} -o {Pre}.bracken1.txt -w {Pre}.bracken2.txt -l S -t 10  -i {Pre}.report.txt',shell=True,stdout=kkf,stderr=kkf) # {Pre}.bracken1.txt：丰度重估结果 {Pre}.bracken2.txt：详细的丰度信息，包括不同分类层级的读数数量            
                if fq1 and fq2:   
                    subprocess.run(f'kraken2 --db {Krdb} --threads {threads} --output {Pre}_2.out.txt --report {Pre}_2.report.txt {fq1} {fq2}',shell=True,stdout=kkf,stderr=kkf)
                    subprocess.run(f'bracken -d {Krdb} -o {Pre}_2.bracken1.txt -w {Pre}_2.bracken2.txt -l S -t 10  -i {Pre}_2.report.txt',shell=True,stdout=kkf,stderr=kkf)
                    testbrkdb = pd.read_table(f'{Pre}_2.report.txt',header=None)
                    if 'S4' in testbrkdb[3]:
                        subprocess.run(f'bracken -d {Krdb} -o {Pre}_2_Sub.bracken1.txt -w {Pre}_2_Sub.bracken2.txt -l S3 -t 10  -i {Pre}_2.report.txt',shell=True,stdout=kkf,stderr=kkf)
                    else:
                        if 'S3' in testbrkdb[3]:
                            subprocess.run(f'bracken -d {Krdb} -o {Pre}_2_Sub.bracken1.txt -w {Pre}_2_Sub.bracken2.txt -l S2 -t 10  -i {Pre}_2.report.txt',shell=True,stdout=kkf,stderr=kkf)
                        else:
                            subprocess.run(f'bracken -d {Krdb} -o {Pre}_2_Sub.bracken1.txt -w {Pre}_2_Sub.bracken2.txt -l S1 -t 10  -i {Pre}_2.report.txt',shell=True,stdout=kkf,stderr=kkf)                
                elif fq1:       
                    subprocess.run(f'kraken2 --db {Krdb} --threads {threads} --output {Pre}_2.out.txt --report {Pre}_2.report.txt {fq1}',shell=True,stdout=kkf,stderr=kkf)
                    subprocess.run(f'bracken -d {Krdb} -o {Pre}_2.bracken1.txt -w {Pre}_2.bracken2.txt -l S -t 10  -i {Pre}_2.report.txt',shell=True,stdout=kkf,stderr=kkf)
                    testbrkdb = pd.read_table(f'{Pre}_2.report.txt',header=None)
                    if 'S4' in testbrkdb[3]:
                        subprocess.run(f'bracken -d {Krdb} -o {Pre}_2_Sub.bracken1.txt -w {Pre}_2_Sub.bracken2.txt -l S3 -t 10  -i {Pre}_2.report.txt',shell=True,stdout=kkf,stderr=kkf)
                    else:
                        if 'S3' in testbrkdb[3]:
                            subprocess.run(f'bracken -d {Krdb} -o {Pre}_2_Sub.bracken1.txt -w {Pre}_2_Sub.bracken2.txt -l S2 -t 10  -i {Pre}_2.report.txt',shell=True,stdout=kkf,stderr=kkf)
                        else:
                            subprocess.run(f'bracken -d {Krdb} -o {Pre}_2_Sub.bracken1.txt -w {Pre}_2_Sub.bracken2.txt -l S1 -t 10  -i {Pre}_2.report.txt',shell=True,stdout=kkf,stderr=kkf)
    
        # 分类水平数据提取
        if fq1 or fq2:
            Pre2=f'{Pre}_2'
        else:
            Pre2 = Pre 
        if lelID != 'nolevel':
            level = lelID.split(',')[1]
            krakenfile = f'{Pre2}.report.txt'
            tkid= lelID.split(',')[0]
            tkidl = os.popen(f'grep {tkid} {krakenfile}|cut -f4').read().strip()
            taxlist1 = proc_kra1(krakenfile,tkid,level)
            taxlist1 = [int(i) for i in taxlist1]
            if fq2:
                exreadsID1(taxlist1,f'{Pre2}.out.txt',f'{Pre}.R1.fastq.gz',f'{Pre}.R2.fastq.gz')
                fq1 = f'{tkid}.1.fastq'
                fq2 = f'{tkid}.2.fastq'
            else:
                exreadsID1(taxlist1,f'{Pre2}.out.txt',f'{Pre}.R1.fastq.gz',0)
                fq1 = f'{tkid}.1.fastq'
                fq2 = 0
        else:  # meta + nolevel提取细菌属
            if method == 'meta':
                if not os.path.isfile(f'2.1.fastq'):
                    level = 'D'
                    krakenfile = f'{Pre2}.report.txt'
                    tkid= 2
                    tkidl = 'D'
                    taxlist1 = proc_kra1(krakenfile,tkid,level)
                    taxlist1 = [int(i) for i in taxlist1]
                    if fq2:
                        exreadsID1(taxlist1,f'{Pre2}.out.txt',f'{Pre}.R1.fastq.gz',f'{Pre}.R2.fastq.gz')
                        fq1 = f'{tkid}.1.fastq'
                        fq2 = f'{tkid}.2.fastq'
                    else:
                        exreadsID1(taxlist1,f'{Pre2}.out.txt',f'{Pre}.R1.fastq.gz',0)
                        fq1 = f'{tkid}.1.fastq'
                        fq2 = 0

        if not os.path.isfile(f'assm_{method}_{pts}_ok'):   # 如果没有组装完成文件            
            methodlist = method.split(',')
            if len(methodlist) == 1:
                method = methodlist[0]
                if method in ['flye','canu','unicycler','masurca','spades','raven','wtdbg2','miniasm','meta']:
                    denovo_asb(inf,fq1,fq2,threads,Pre,pts,pst,method,asmt,f,outputfa)
                else:
                    if ref != 'noref':
                        reassm_fun(inf,fq1,fq2,threads,Pre,pts,pst,method,asmt,f,outputfa)
                    else:
                        print('有参组装未提供参考基因组')
                        sys.exit()               
            elif len(methodlist) == 2:
                denovomethod = methodlist[0]
                reasmmethod = methodlist[1]
                denovo_asb(inf,fq1,fq2,threads,Pre,pts,pst,denovomethod,asmt,f,outputfa)
                if ref != 'noref':
                    reassm_fun(inf,fq1,fq2,threads,Pre,pts,pst,reasmmethod,asmt,f,outputfa='noforce')
                else:
                    print('有参组装未提供参考基因组')
                    sys.exit()
    if method != 'meta':
        if os.path.isfile(outputfa):
            finalfa = f'{Pre}.final.fasta'
            with open('Anno1.log','w') as f:
                if os.path.isfile(polifa): # polifa = {Pre}.polish.fasta三代拼接后的最终抛光文件
                    renamefa(polifa,finalfa)
                else:
                    renamefa(outputfa,finalfa)               
                
                 # 根据三代不同序列类型，运用 minimap2 的不同映射命令            
                if asmt == 'longasm' or asmt == 'longref':       
                    if long_type == 'Nanopore':
                        subprocess.run(f'minimap2 -ax map-ont {finalfa} {inf} -t {threads} |samtools sort -o {Pre}.sorted.bam',shell=True,stdout=f,stderr=f)
                    elif long_type == 'PacBio_CLR':
                        subprocess.run(f'minimap2 -ax map-pb {finalfa} {inf} -t {threads} |samtools sort -o {Pre}.sorted.bam',shell=True,stdout=f,stderr=f)
                    elif long_type == 'PacBio_CCS':
                        subprocess.run(f'minimap2 -ax map-hifi {finalfa} {inf} -t {threads} |samtools sort -o {Pre}.sorted.bam',shell=True,stdout=f,stderr=f)
                 # 二代数据，运用 bwa 的映射命令                     
                elif asmt == 'shortasm' or asmt == 'shortref':
                    subprocess.run(f'bwa index {finalfa}',shell=True)
                    if fq2:
                        subprocess.run(f'bwa mem {finalfa} {fq1} {fq2} -t {threads}|samtools sort -o {Pre}_ngs.sorted.bam',shell=True)
                    else:
                        subprocess.run(f'bwa mem {finalfa} {fq1} -t {threads}|samtools sort -o {Pre}_ngs.sorted.bam',shell=True)
                # asmt == 'shortlongasm'
                else:
                    if long_type == 'Nanopore':
                        subprocess.run(f'minimap2 -ax map-ont {finalfa} {inf} -t {threads} |samtools sort -o {Pre}.sorted.bam',shell=True,stdout=f,stderr=f)
                    elif long_type == 'PacBio_CLR':
                        subprocess.run(f'minimap2 -ax map-pb {finalfa} {inf} -t {threads} |samtools sort -o {Pre}.sorted.bam',shell=True,stdout=f,stderr=f)
                    elif long_type == 'PacBio_CCS':
                        subprocess.run(f'minimap2 -ax map-hifi {finalfa} {inf} -t {threads} |samtools sort -o {Pre}.sorted.bam',shell=True,stdout=f,stderr=f)
                    if fq2:
                        subprocess.run(f'bwa mem {finalfa} {fq1} {fq2} -t {threads}|samtools sort -o {Pre}_ngs.sorted.bam',shell=True)
                    elif not fq2:
                        subprocess.run(f'bwa mem {finalfa} {fq1} -t {threads}|samtools sort -o {Pre}_ngs.sorted.bam',shell=True)
                
                # 三代数据映射的 {Pre}.sorted.bam 文件，生成1000bp窗口和1bp窗口的覆盖度信息                   
                if os.path.isfile(f'{Pre}.sorted.bam'):
                    subprocess.run(f'samtools index {Pre}.sorted.bam',shell=True,stdout=f,stderr=f)
                    subprocess.run(f'mosdepth -b 1000 {Pre} {Pre}.sorted.bam',shell=True)
                    subprocess.run(f'mosdepth -b 1 {Pre}_1 {Pre}.sorted.bam',shell=True)
                    subprocess.run(f'gunzip -f {Pre}.regions.bed.gz',shell=True)
                    subprocess.run(f'gunzip -f {Pre}_1.regions.bed.gz',shell=True) 
                    # 根据 contig 名称分组分别保存到单独的 .bed 文件   
                    if not os.path.isdir('Contigbedfile'):
                        os.makedirs('Contigbedfile')
                    contigbed = pd.read_table(f'{Pre}.regions.bed',header=None)
                    contig1bed = pd.read_table(f'{Pre}_1.regions.bed',header=None)
                    contigbed.groupby(0).apply(lambda x:x.to_csv(f'Contigbedfile/{x[0].tolist()[0]}.bed',index=False,header=False,sep='\t'))
                    contig1bed.groupby(0).apply(lambda x:x.to_csv(f'Contigbedfile/{x[0].tolist()[0]}_dis1.bed',index=False,header=False,sep='\t'))
                    for bed in os.listdir('Contigbedfile'):
                        if not bed.endswith('_dis1.bed'):
                            ttPre = bed.replace('.bed','')  # ttPre = 去除 .bed文件 的后缀.bed，只保留contig名称
                            if int(os.popen(f'cat Contigbedfile/{bed}|wc -l').read()) < 10:  # 如果 {ttPre}.bed文件 行数小于10。用 {ttPre}_dis1.bed 文件替换 {ttPre}.bed 文件
                                subprocess.run(f'mv Contigbedfile/{ttPre}_dis1.bed Contigbedfile/{ttPre}.bed ',shell=True)   
                    contig1bed[contig1bed[3]==0].to_csv('mask.bed',sep='\t',index=False,header=False)   # mask.bed = contig1bed文件 中覆覆盖为0的区域
                
                # 二代数据映射的BAM文件，生成1000bp窗口以及1bp窗口的覆盖度信息等            
                if os.path.isfile(f'{Pre}_ngs.sorted.bam'):
                    subprocess.run(f'samtools index {Pre}_ngs.sorted.bam',shell=True,stdout=f,stderr=f)
                    subprocess.run(f'mosdepth -b 1000 {Pre}_ngs {Pre}_ngs.sorted.bam',shell=True)
                    subprocess.run(f'mosdepth -b 1 {Pre}_ngs_1 {Pre}_ngs.sorted.bam',shell=True)
                    subprocess.run(f'gunzip -f {Pre}_ngs.regions.bed.gz',shell=True)
                    subprocess.run(f'gunzip -f {Pre}_ngs_1.regions.bed.gz',shell=True)
                    if not os.path.isdir('Contigbedfile'):
                        os.makedirs('Contigbedfile')
                    contigbed = pd.read_table(f'{Pre}_ngs.regions.bed',header=None)
                    contig1bed = pd.read_table(f'{Pre}_ngs_1.regions.bed',header=None)
                    contigbed.groupby(0).apply(lambda x:x.to_csv(f'Contigbedfile/{x[0].tolist()[0]}.bed',index=False,header=False,sep='\t'))
                    contig1bed.groupby(0).apply(lambda x:x.to_csv(f'Contigbedfile/{x[0].tolist()[0]}_dis1.bed',index=False,header=False,sep='\t'))
                    for bed in os.listdir('Contigbedfile'):
                        if not bed.endswith('_dis1.bed'):
                            ttPre = bed.replace('.bed','')
                            if int(os.popen(f'cat Contigbedfile/{bed}|wc -l').read()) < 10:
                                subprocess.run(f'mv Contigbedfile/{ttPre}_dis1.bed Contigbedfile/{ttPre}.bed ',shell=True)
                    contig1bed[contig1bed[3]==0].to_csv('mask.bed',sep='\t',index=False,header=False)       # mask.bed 文件包含 contig1bed文件 中覆覆盖为0的区域
                
                if method == 'canu':   # 三代拼接软件canu的特殊处理，它具有额外信息如：长度和是否成环。flye也有，但输出格式以flye为模板，因此无需更改格式
                    subprocess.run(f'seqkit fx2tab canu/canu.contigs.fasta -n > canu.txt',shell=True)
                    canudb = pd.read_table('canu.txt',sep=' ',header=None)
                    canudb['len'] = canudb[1].str.replace('len=','').str.strip().astype('int')
                    canudb = canudb.sort_values('len',ascending=False)
                    canudb['index'] = canudb.index+1
                    canudb['contig'] = 'contig' + '_' + canudb['index'].astype('str')
                    canudb['cir'] = canudb[6].str.replace('suggestCircular=','')  #从第 6 列提取是否建议环状（circular）的信息
                    #canudb[['contig','cir']].to_csv(f'canu_sum.tsv',sep='\t',index=False)   # 将序列对应的是否建议成环、长度提取到canu_sum.tsv
                if not os.path.exists('flye_output'):
                    os.makedirs('flye_output')
                    subprocess.run(f'seqkit fx2tab {finalfa} -n -l > flye_output/tmp.stat.txt',shell=True)  # 输出每条序列的名称和长度到tmp.stat.txt
                    flyedb = pd.read_table('flye_output/tmp.stat.txt',header=None)                          
                    flyedb['是否成环'] = '-'                                                                 # 新增一列名为'是否成环',初始为'-' 
                    flyedb.rename(columns = {0:'序列名称',1:'序列长度'},inplace=True)                         # 列索引为 0 的列重命名为 序列名称，列索引为 1 的列重命名为 序列长度
                    if method == 'canu': 
                        for ctg in flyedb['序列名称'].tolist():
                            flyedb.loc[flyedb['序列名称']==ctg,'是否成环'] = canudb.loc[canudb['contig']==ctg,'cir'].tolist()[0]
                    
                    if os.path.isfile(f'{Pre}.regions.bed'):    # {Pre}.regions.bed，三代数据的1000bp比对文件
                        mosdb = pd.read_table(f'{Pre}.regions.bed',header=None)
                    else: 
                        mosdb = pd.read_table(f'{Pre}_ngs.regions.bed',header=None)  # 二代
                    # 将覆盖度文件（mosdb）与组装结果（flyedb）合并，以便生成包含每个 contig 的名称、长度、平均深度和是否成环信息的汇总文件assembly_info.txt
                    flyedb = mosdb.groupby(0).agg('mean').merge(flyedb,left_on=0,right_on='序列名称')
                    flyedb.rename(columns={3:'平均深度'},inplace=True)
                    flyedb = flyedb[['序列名称','序列长度','平均深度','是否成环']]
                    flyedb.平均深度 = flyedb.平均深度.round()
                    flyedb.sort_values('序列长度',axis=0,inplace=True,ascending=False)
                    flyedb.to_csv(f'flye_output/assembly_info.txt',sep='\t',index=False)
                else: # 已经存在 flye_output/assembly_info.txt 文件的情况下，更新汇总文件的内容 
                    transdb = pd.read_table('transname.tsv',header=None,names=['oldname','newname'])
                    subprocess.run(f'seqkit fx2tab {finalfa} -n -l > flye_output/tmp.stat.txt',shell=True)
                    newcontigdb = pd.read_table('flye_output/tmp.stat.txt',header=None)
                    flyedb = pd.read_table('flye_output/assembly_info.txt')
                    flyedb = flyedb.merge(transdb,left_on='序列名称',right_on='oldname')
                    flyedb = flyedb.merge(newcontigdb,left_on='newname',right_on=0)
                    flyedb['序列名称'] = flyedb['newname']
                    flyedb['序列长度'] = flyedb[1]
                    flyedb = flyedb.sort_values('序列长度',ascending=False)
                    flyedb = flyedb[['序列名称','序列长度','平均深度','是否成环']]
                    flyedb.to_csv(f'flye_output/assembly_info.txt',sep='\t',index=False)
                
                #---1.0.1 质粒鉴定及分型
                with open('asb.log','a') as f:  # 使用plasflow对序列进行基因组/质粒预测，生成{Pre}_plaspredict.tsv； 使用staramr工具的plasmidfinder进行质粒鉴定生成staramr_result/plasmidfinder.tsv 
                    subprocess.run(f'/home/dell/miniconda3/bin/conda run --no-capture-output -n plasflow PlasFlow.py --input {finalfa} --output {Pre}_plaspredict.tsv',shell=True,stdout=f,stderr=f) 
                    subprocess.run(f'staramr search  {finalfa} -o staramr_result -n 30',shell=True,stdout=f,stderr=f)                               
                def join_plasmid(group):     # 定义join_plasmid函数：将每个 contig 对应的质粒预测结果进行合并
                    plasmids = '|'.join(group['Plasmid'])  # 将每个contig可能包含的多个质粒预测结果（plasmid 列的值），用'|'进行分割合并
                    newdb = group
                    newdb['Plasmid'] = plasmids
                    newdb = newdb.iloc[0,:]   # 经过上面一步后,同一个contig的多行同时包含所有预测结果，因此保留第一个即可
                    return newdb              
                plasmiddb = pd.read_table('staramr_result/plasmidfinder.tsv')  # 读取staramr的质粒（分型）预测结果
                rawindexname = plasmiddb.index.name   # 保存 plasmiddb 索引名称
                if plasmiddb.shape[0] >1:   # plasmiddb行数＞1
                    plasmiddb = plasmiddb.groupby('Contig').apply(join_plasmid)  # 按照 contig名称 进行分组，调用join_plasmid函数
                    plasmiddb.index.name = rawindexname  # 将索引名称改回原来的
                plasflowdb = pd.read_table(f'{Pre}_plaspredict.tsv') # 读取 PlasFlow 预测结果文件
                plasflowdb = plasflowdb[['contig_name','label']]     # 仅保留'contig_name','label'两列
                flyedb = pd.read_table('flye_output/assembly_info.txt')
                flyedb = flyedb.merge(plasflowdb,left_on='序列名称',right_on='contig_name').drop('contig_name',axis=1)  # 按照'序列名称'和'contig_name'配对进行合并，最后去除'contig_name'列
                flyedb = flyedb.merge(plasmiddb,left_on='序列名称',right_on='Contig',how='left').drop('Contig',axis=1)  # 按照contig和序列名称进行配对，contig无但序列名称有的也将会被保留
                flyedb = flyedb.rename(columns={'label':'基因组/质粒','Plasmid':'质粒分型'})
                #---提取将质粒gbk和主基因组gbk区分
                flyedb.fillna('-',inplace=True)  #用 '-' 来填充 flyedb 数据框中所有的缺失值（NaN）
                sys.stdout.flush()
                flyedb.to_csv('tmp1.tsv',sep='\t',index=False)
                flyedb['序列长度'] = flyedb['序列长度'].astype('int')
                # 选择plasflow预测为plasmid 或者 staramr预测不为'—'，序列长度小于 1,000,000 的记录。将筛选后的序列名称列转换为列表存储在 plasmidlist 中
                plasmidlist = flyedb.loc[(flyedb['基因组/质粒'].str.contains('plasmid')) | (flyedb['质粒分型']!='-') & (flyedb['序列长度'] < 1000000),'序列名称'].tolist()
                flyedb['占比'] = flyedb['序列长度']/flyedb['序列长度'].sum() # 新加一列'占比':计算每条 contig 占总序列长度的比例
                flyedb['占比'] =flyedb['占比'].round(2)                    # 保留小数点后两位
                flyedb = flyedb[['序列名称','序列长度','平均深度','是否成环','基因组/质粒','质粒分型','占比']]
                flyedb.to_csv(f'flye_output/assembly_info.txt',sep='\t',index=False)
                
                # 利用 seqkit 工具对最终组装的 finalfa 文件进行统计，保存为 finalfasta.tsv 文件
                if os.path.isfile(finalfa):
                    subprocess.run(f'''seqkit stat -a -T -G N {finalfa} > ragtag_sum.tsv''',shell=True)
                    ragtagdb = pd.read_table('ragtag_sum.tsv')
                    ragtagdb = ragtagdb[['num_seqs','sum_len','max_len','min_len','sum_gap','N50','GC(%)']]
                    ragtagdb['样本名称'] = Pre
                    ragtagdb.rename(columns={'num_seqs':'contig数量','sum_len':'总长度','min_len':'最小contig长度','max_len':'最大contig长度'},inplace=True)
                    ragtagdb['N比例(%)'] = ((ragtagdb['sum_gap']/ragtagdb['总长度'])*100).round(2)
                    ragtagdb = ragtagdb[['样本名称','contig数量','总长度','最大contig长度','最小contig长度','N50','GC(%)','N比例(%)']]
                    ragtagdb.to_csv('finalfasta.tsv',index=False,sep='\t')
                with open('asb.log','a') as f:  # 调用 Prokka 对基因组进行注释，整理并提取关键的注释信息
                    subprocess.run(f'prokka --force --outdir {Pre}_prokka --prefix {Pre} --addgenes --cpus {threads} {finalfa}',shell=True,stdout=f,stderr=f)  # 调用 Prokka 对（finalfa）进行注释，结果输出到 {Pre}_prokka 目录中
                    #prokkadb = pd.read_table(f'{Pre}_prokka/{Pre}.tsv',names=['注释标签','类型','长度_bp','基因名','EC号','COG号','功能描述']).drop(0,axis=0)
                    prkskpn = int(os.popen(f'''grep '##sequence-region' {Pre}_prokka/{Pre}.gff|wc -l''').read())+1      # 统计需要跳过的行数
                    prkfan = int(os.popen(f'''grep -n '#FASTA' {Pre}_prokka/{Pre}.gff''').read().split(':')[0])         # 确定 .gff 文件中 '#FASTA' 标记的行号,表示序列部分的开始
                    prkskfn = int(os.popen(f'''cat {Pre}_prokka/{Pre}.gff|wc -l''').read()) - prkfan+1                  # 计算出序列部分的总行数
                    # 从 Prokka 注释生成的 .gff 文件中读取非序列部分的数据，提取关键信息到{Pre}.prokka.tsv文件
                    prokkadb = pd.read_table(f'{Pre}_prokka/{Pre}.gff',skiprows=prkskpn,skipfooter=prkskfn,engine='python',header=None,names=['染色体','数据库','类型','起始位置','终止位置','t1','链方向','t2','注释1'])
                    prokkadb = prokkadb[['染色体','类型','起始位置','终止位置','链方向','注释1']]      
                    prokkadb['基因名称'] = prokkadb['注释1'].str.split('Name=').str[1].str.split(';').str[0]
                    prokkadb['locus标签'] = prokkadb['注释1'].str.split('ID=').str[1].str.split(';').str[0]  # 从'注释1'中分别提取出'基因名称'和'locus标签'
                    prokkadb.fillna('-',inplace=True)
                    prokkadb = prokkadb[['染色体','类型','起始位置','终止位置','链方向','基因名称','locus标签']]
                    prokkadb.to_csv(f'{Pre}.prokka.tsv',sep='\t',index=False)
                    
                    # 对prokka注释后的.gbk文件中可能包含的中文字符 '月’ 用空格进行替代
                    monthchin = int(os.popen(f'''grep '月' {Pre}_prokka/{Pre}.gbk|wc -l''').read().strip())
                    if monthchin:
                        if os.popen(f'head -n 1 {Pre}_prokka/{Pre}.gbk').read().strip()[72] == '月':
                            subprocess.run(f'''sed -i 's/月/  /g' {Pre}_prokka/{Pre}.gbk''',shell=True)
                        else:
                            subprocess.run(f'''sed -i 's/月/ /g' {Pre}_prokka/{Pre}.gbk''',shell=True)
                    subprocess.run(f'cp {Pre}_prokka/{Pre}.gbk tt.gbk',shell=True)
                    
                    # 提取{Pre}.gff文件中类型为 repeat_region 的重复区域，输出{Pre}.repeat.tsv。该部分的重复区域由prokka注释结果产生，参考意义不大。该段在元件预测中有
                    #ignum = int(os.popen(f'''cat {Pre}_prokka/{Pre}.gff|grep '##' |wc -l ''').read().strip())-1
                    #afile = pd.read_table(f'{Pre}_prokka/{Pre}.gff',skiprows=ignum,low_memory=False,header=None)
                    #afile = afile.loc[~afile[3].isna()]
                    #repeatfile = afile[afile[2]=='repeat_region']  # prokka注释的repeat_region并不一定准确
                    #repeatfile = afile[[0,1,2,3,4,8]]
                    #repeatfile.rename(columns={0:'Contig名称',1:'数据库',2:'类型',3:'序列开始',4:'序列结尾',8:'结果注释'},inplace=True)
                    #repeatfile.to_csv(f'{Pre}.repeat.tsv',sep='\t',index=False)
                    
                # plasmidlist中的序列与prokka注释的基因结合,记录在{contig_id_to_extract}.gbk文件中
                if len(plasmidlist) != 0:  # 质粒序列
                    input_gbk = f'{Pre}_prokka/{Pre}.gbk'
                    records = SeqIO.parse(input_gbk, "genbank")  # 使用 Biopython 的 SeqIO 模块来解析 GenBank 格式的文件
                    for contig_id_to_extract in plasmidlist:  # 遍历 plasmidlist 中的每个 contig ID
                        for record in records:   # 遍历 records 中的每个记录对象
                            if record.id == contig_id_to_extract:   # 两者相等时
                                with open(f'{Pre}_prokka/{contig_id_to_extract}.gbk', "w") as output_handle:  
                                    SeqIO.write(record, output_handle, "genbank")   # 保存提取到的 record 信息至 {Pre}_prokka/{contig_id_to_extract}.gbk             
                                with open('cgview.log','a') as cgvf: #调用 cgview_builder_cli.rb 脚本，.gbk 转换成用于可视化的 .json 格式   
                                    subprocess.run(f'ruby /data/deploy/bio-elite/bio/script/cgview_builder_cli.rb -s {Pre}_prokka/{contig_id_to_extract}.gbk -o {contig_id_to_extract}.json -c /data/deploy/bio-elite/bio/script/CGview.yaml -n {contig_id_to_extract}',shell=True,stdout=cgvf,stderr=cgvf)
                    
                    
                    records = list(SeqIO.parse(input_gbk, "genbank"))
                    filtered_records = [record for record in records if record.id not in plasmidlist]
                    with open(f'{Pre}_prokka/main.gbk', "w") as output_handle:  # 去除包含在plasmidlist中的序列ID 生成主基因组 main.gbk 文件
                        SeqIO.write(filtered_records, output_handle, "genbank")
                    with open('cgview.log','a') as cgvf:        
                        subprocess.run(f'ruby /data/deploy/bio-elite/bio/script/cgview_builder_cli.rb -s {Pre}_prokka/main.gbk -o main.json -c /data/deploy/bio-elite/bio/script/CGview.yaml -n main',shell=True,stdout=cgvf,stderr=cgvf)
                else: # 没有质粒序列时的处理方式
                    with open('cgview.log','a') as cgvf:
                        subprocess.run(f'ruby /data/deploy/bio-elite/bio/script/cgview_builder_cli.rb -s {Pre}_prokka/{Pre}.gbk -o main.json -c /data/deploy/bio-elite/bio/script/CGview.yaml -n {Pre}',shell=True,stdout=cgvf,stderr=cgvf)
                if not os.path.isdir(f'{Pre}_bin_genome_out'):
                    os.makedirs(f'{Pre}_bin_genome_out')
                
                subprocess.run(f'cp {finalfa} {Pre}_bin_genome_out/{Pre}.fna',shell=True)
                # seqkit + checkm: 确定对应物种的完整性、污染率 = {Pre}.assemble.result.tsv
                with open('checkm2','w') as cmf:  
                    subprocess.run(f'/home/dell/miniconda3/bin/conda run --no-capture-output  -n cm2  checkm2 predict -i  {Pre}_bin_genome_out -x fna -o checkm2_out -t {threads} --force --database_path /data1/shanghai_pip/meta_genome/uniref100.KO.1.dmnd',shell=True,stdout=cmf,stderr=cmf)  # checkm2进行质量评估
                    if os.path.isfile('checkm2_out/quality_report.tsv'):    # 将quality_report.tsv转化成可读性更强的{Pre}.checkm.tsv
                        checkmdb = pd.read_table('checkm2_out/quality_report.tsv')
                        checkmdb['样本名称'] = Pre
                        checkmdb.rename(columns={'Completeness':'完整性','Contamination':'污染率'},inplace=True)
                        checkmdb['物种名称'] = species 
                        checkmdb[['样本名称','物种名称','污染率','完整性']].to_csv(f'{Pre}.checkm.tsv',sep='\t',index=False)
            subprocess.run(f'seqkit stat {Pre}.final.fasta -b -T -a > {Pre}.fasum.tsv',shell=True)
            fadb = pd.read_table(f'{Pre}.fasum.tsv')
            fadb['Sample'] = fadb['file'].str.replace('.final.fasta','')
            Assdict = {}
            Assdict['Contig数量'] = fadb['num_seqs'].tolist()[0]
            Assdict['N50长度'] = fadb['N50'].tolist()[0]
            Assdict['组装基因组长度'] = fadb['sum_len'].tolist()[0]
            Assdict['最长片段长度'] = fadb['max_len'].tolist()[0]
            Assdict['污染率'] = '-'
            Assdict['完整性'] = '-'
            Assdict['主基因组是否成环'] = '-'
            if os.path.isfile(f'{Pre}.checkm.tsv'):
                cdb = pd.read_table(f'{Pre}.checkm.tsv')
                Assdict['污染率'] = round(float(cdb['污染率'].tolist()[0]),2)
                Assdict['完整性'] = round(float(cdb['完整性'].tolist()[0]),2)
            Assdb = pd.DataFrame(Assdict,index=[0])
            Assdb['样本名称'] = Pre
            Assdb = Assdb[['样本名称','Contig数量','N50长度','最长片段长度','主基因组是否成环','污染率','完整性']]
            Assdb['主基因组是否成环'] = flyedb['是否成环'].tolist()[0]
            Assdb.to_csv(f'{Pre}.assemble.result.tsv',sep='\t',index=False)
            # genovi中的genovi工具对输入文件进行功能注释
            with open('genovi.log','w') as fv:
                try:  # 使用genovi 对 {Pre}_prokka/{Pre}.gbk进行注释，主要获取 COG 分类
                    subprocess.run(f'/home/dell/miniconda3/bin/conda run -n genovi genovi -i {Pre}_prokka/{Pre}.gbk -o {Pre}_genovi -s draft',shell=True,stdout=fv,stderr=fv)
                    cogdb = pd.read_table(f'{Pre}_genovi/{Pre}_genovi_COG_Classification.csv',sep=',',header=1)
                    cogdb.iloc[-2,:].to_csv('Cog_summary.tsv',sep='\t',header=False)
                except:
                    pass

        else:
            raise Exception(f"{method}组装报错可能原因:\n1.数据量过少，提高数据量后继续分析\n2.样本中含有一定量的杂菌污染")

    time.sleep(0.5)
   
# 定义Annotate_func函数：
def Annotate_func(Pre,threads):
    with open('Anno.log','w') as f1:
        finalfa = f'{Pre}.final.fasta'
        renamefa(finalfa,finalfa)   #给传入fasta用
        if not os.path.isfile(f'{Pre}_prokka/{Pre}.tsv'):  # 未调用asb_func函数
            if os.path.isfile(finalfa):  
                # 利用 seqkit 工具对最终组装的 finalfa 文件进行统计，保存为 finalfasta.tsv 文件
                subprocess.run(f'''seqkit stat -a -T -G N {Pre}.final.fasta > ragtag_sum.tsv''',shell=True)  # 通过 seqkit 工具进行统计分析
                ragtagdb = pd.read_table('ragtag_sum.tsv')
                ragtagdb = ragtagdb[['num_seqs','sum_len','max_len','min_len','sum_gap','N50','GC(%)']]
                ragtagdb['样本名称'] = Pre
                ragtagdb.rename(columns={'num_seqs':'contig数量','sum_len':'总长度','min_len':'最小contig长度','max_len':'最大contig长度'},inplace=True)
                ragtagdb['N比例(%)'] = ((ragtagdb['sum_gap']/ragtagdb['总长度'])*100).round(2)
                ragtagdb = ragtagdb[['样本名称','contig数量','总长度','最大contig长度','最小contig长度','N50','GC(%)','N比例(%)']]
                ragtagdb.to_csv('finalfasta.tsv',index=False,sep='\t')          
            subprocess.run(f'prokka --force --outdir {Pre}_prokka --prefix {Pre} --addgenes --cpus {threads} {Pre}.final.fasta',shell=True,stdout=f1,stderr=f1)  
            prkskpn = int(os.popen(f'''grep '##sequence-region' {Pre}_prokka/{Pre}.gff|wc -l''').read())+1  # 统计需要跳过的行数
            prkfan = int(os.popen(f'''grep -n '#FASTA' {Pre}_prokka/{Pre}.gff''').read().split(':')[0])     # 确定 .gff 文件中 '#FASTA' 标记的行号,表示序列部分的开始
            prkskfn = int(os.popen(f'''cat {Pre}_prokka/{Pre}.gff|wc -l''').read()) - prkfan+1              # 计算出序列部分的总行数
            prokkadb = pd.read_table(f'{Pre}_prokka/{Pre}.gff',skiprows=prkskpn,skipfooter=prkskfn,engine='python',header=None,names=['染色体','数据库','类型','起始位置','终止位置','t1','链方向','t2','注释1'])
            prokkadb = prokkadb[['染色体','类型','起始位置','终止位置','链方向','注释1']]
            prokkadb['基因名称'] = prokkadb['注释1'].str.split('Name=').str[1].str.split(';').str[0]
            prokkadb['locus标签'] = prokkadb['注释1'].str.split('ID=').str[1].str.split(';').str[0]
            prokkadb.fillna('-',inplace=True)
            prokkadb = prokkadb[['染色体','类型','起始位置','终止位置','链方向','基因名称','locus标签']]
            prokkadb.to_csv(f'{Pre}.prokka.tsv',sep='\t',index=False)
            
            # 对prokka注释后的.gbk文件中可能包含的中文字符 '月’ 用空格进行替代
            monthchin = int(os.popen(f'''grep '月' {Pre}_prokka/{Pre}.gbk|wc -l''').read().strip())
            if monthchin:
                if '月' in os.popen(f'head -n 1 {Pre}_prokka/{Pre}.gbk').read().strip(): 
                    subprocess.run(f'''sed -i 's/月/  /g' {Pre}_prokka/{Pre}.gbk''',shell=True)
                else:
                    subprocess.run(f'''sed -i 's/月/ /g' {Pre}_prokka/{Pre}.gbk''',shell=True)
            subprocess.run(f'cp {Pre}_prokka/{Pre}.gbk tt.gbk',shell=True) 
            
            # 提取{Pre}.gff文件中类型为 repeat_region 的重复区域，输出{Pre}.repeat.tsv。被注释掉的原因为‘prokka注释的repeat_region并不一定准确“，后续有其他软件替代其功能
            #ignum = int(os.popen(f'''cat {Pre}_prokka/{Pre}.gff|grep '##' |wc -l ''').read().strip())-1
            #afile = pd.read_table(f'{Pre}_prokka/{Pre}.gff',skiprows=ignum,low_memory=False,header=None)
            #afile = afile.loc[~afile[3].isna()]
            #repeatfile = afile[afile[2]=='repeat_region']
            #repeatfile = afile[[0,1,2,3,4,8]]
            #repeatfile.rename(columns={0:'Contig名称',1:'数据库',2:'类型',3:'序列开始',4:'序列结尾',8:'结果注释'},inplace=True)
            #repeatfile.to_csv(f'{Pre}.repeat.tsv',sep='\t',index=False)

            if method == 'canu':   # 三代拼接软件canu的特殊处理，它具有额外信息如：长度和是否成环。flye也有，但输出格式以flye为模板，因此无需更改格式
                subprocess.run(f'seqkit fx2tab canu/canu.contigs.fasta -n > canu.txt',shell=True)
                canudb = pd.read_table('canu.txt',sep=' ',header=None)
                canudb['len'] = canudb[1].str.replace('len=','').str.strip().astype('int')
                canudb = canudb.sort_values('len',ascending=False)
                canudb['index'] = canudb.index+1
                canudb['contig'] = 'contig' + '_' + canudb['index'].astype('str')
                canudb['cir'] = canudb[6].str.replace('suggestCircular=','')  #从第 6 列提取是否建议环状（circular）的信息
                #canudb[['contig','cir']].to_csv(f'canu_sum.tsv',sep='\t',index=False)   # 将序列对应的是否建议成环、长度提取到canu_sum.tsv
            if not os.path.exists('flye_output'):
                os.makedirs('flye_output')
                subprocess.run(f'seqkit fx2tab {finalfa} -n -l > flye_output/tmp.stat.txt',shell=True)  # 输出每条序列的名称和长度到tmp.stat.txt
                flyedb = pd.read_table('flye_output/tmp.stat.txt',header=None)                          
                flyedb['是否成环'] = '-'                                                                 # 新增一列名为'是否成环',初始为'-' 
                flyedb.rename(columns = {0:'序列名称',1:'序列长度'},inplace=True)                         # 列索引为 0 的列重命名为 序列名称，列索引为 1 的列重命名为 序列长度
                if method == 'canu': 
                    for ctg in flyedb['序列名称'].tolist():
                        flyedb.loc[flyedb['序列名称']==ctg,'是否成环'] = canudb.loc[canudb['contig']==ctg,'cir'].tolist()[0]
            flyedb['平均深度'] = '-'                                      # 相比于fastq文件，无法获得'平均深度'信息
            flyedb = flyedb[['序列名称','序列长度','平均深度','是否成环']]
            flyedb.sort_values('序列长度',axis=0,inplace=True,ascending=False)
            flyedb.to_csv(f'flye_output/assembly_info.txt',sep='\t',index=False)
            #---1.0.1 质粒鉴定及分型
            with open('Anno.log','a') as f:  # 使用plasflow对序列进行基因组/质粒预测，生成{Pre}_plaspredict.tsv； 使用staramr工具的plasmidfinder 进行质粒鉴定生成staramr_result/plasmidfinder.tsv 
                subprocess.run(f'/home/dell/miniconda3/bin/conda run --no-capture-output -n plasflow PlasFlow.py --input {finalfa} --output {Pre}_plaspredict.tsv',shell=True,stdout=f,stderr=f) 
                subprocess.run(f'staramr search  {finalfa} -o staramr_result -n 30',shell=True,stdout=f,stderr=f)                               
            def join_plasmid(group):     # 定义join_plasmid函数：将每个 contig 对应的质粒预测结果进行合并
                plasmids = '|'.join(group['Plasmid'])  # 将每个contig可能包含的多个质粒预测结果（plasmid 列的值），用'|'进行分割合并
                newdb = group
                newdb['Plasmid'] = plasmids
                newdb = newdb.iloc[0,:]   # 经过上面一步后,同一个contig的多行同时包含所有预测结果，因此保留第一个即可
                return newdb              
            plasmiddb = pd.read_table('staramr_result/plasmidfinder.tsv')  # 读取staramr的质粒（分型）预测结果
            rawindexname = plasmiddb.index.name   # 保存 plasmiddb 索引名称
            if plasmiddb.shape[0] >1:   # plasmiddb行数＞1
                plasmiddb = plasmiddb.groupby('Contig').apply(join_plasmid)  # 按照 contig名称 进行分组，调用join_plasmid函数
                plasmiddb.index.name = rawindexname  # 将索引名称改回原来的
            plasflowdb = pd.read_table(f'{Pre}_plaspredict.tsv') # 读取 PlasFlow 预测结果文件
            plasflowdb = plasflowdb[['contig_name','label']]     # 仅保留'contig_name','label'两列
            flyedb = pd.read_table('flye_output/assembly_info.txt')
            flyedb = flyedb.merge(plasflowdb,left_on='序列名称',right_on='contig_name').drop('contig_name',axis=1)  # 按照'序列名称'和'contig_name'配对进行合并，最后去除'contig_name'列
            flyedb = flyedb.merge(plasmiddb,left_on='序列名称',right_on='Contig',how='left').drop('Contig',axis=1)  # 按照contig和序列名称进行配对，contig无但序列名称有的也将会被保留
            flyedb = flyedb.rename(columns={'label':'基因组/质粒','Plasmid':'质粒分型'}) 
            #---提取将质粒gbk和主基因组gbk区分
            flyedb.fillna('-',inplace=True)     # 用 '-' 来填充 flyedb 数据框中所有的缺失值（NaN）
            sys.stdout.flush()
            flyedb.to_csv('tmp1.tsv',sep='\t',index=False)
            flyedb['序列长度'] = flyedb['序列长度'].astype('int')
            # 选择plasflow预测为plasmid 或者 staramr预测不为'—'，序列长度小于 1,000,000 的记录。将筛选后的序列名称列转换为列表存储在 plasmidlist 中
            plasmidlist = flyedb.loc[(flyedb['基因组/质粒'].str.contains('plasmid')) | (flyedb['质粒分型']!='-') & (flyedb['序列长度'] < 1000000),'序列名称'].tolist()
            flyedb['占比'] = flyedb['序列长度']/flyedb['序列长度'].sum()  # 新加一列'占比':计算每条 contig 占总序列长度的比例
            flyedb['占比'] =flyedb['占比'].round(2)                     # 保留小数点后两位
            flyedb = flyedb[['序列名称','序列长度','平均深度','是否成环','基因组/质粒','质粒分型','占比']]
            flyedb.to_csv(f'flye_output/assembly_info.txt',sep='\t',index=False)
             
            # plasmidlist中的序列与prokka注释的基因结合,记录在{contig_id_to_extract}.gbk文件中
            if len(plasmidlist) != 0:
                input_gbk = f'{Pre}_prokka/{Pre}.gbk'
                records = SeqIO.parse(input_gbk, "genbank")
                for contig_id_to_extract in plasmidlist:      # 遍历 plasmidlist 中的每个 contig ID
                    for record in records:
                        if record.id == contig_id_to_extract:
                            with open(f'{Pre}_prokka/{contig_id_to_extract}.gbk', "w") as output_handle:
                                SeqIO.write(record, output_handle, "genbank")
                            with open('cgview.log','a') as cgvf: #调用 cgview_builder_cli.rb 脚本，.gbk 转换成用于可视化的 .json 格式   
                                subprocess.run(f'ruby /data/deploy/bio-elite/bio/script/cgview_builder_cli.rb -s {Pre}_prokka/{contig_id_to_extract}.gbk -o {contig_id_to_extract}.json -c /data/deploy/bio-elite/bio/script/CGview.yaml -n {contig_id_to_extract}',shell=True,stdout=cgvf,stderr=cgvf)
                
                records = list(SeqIO.parse(input_gbk, "genbank"))
                filtered_records = [record for record in records if record.id not in plasmidlist]
                with open(f'{Pre}_prokka/main.gbk', "w") as output_handle:
                    SeqIO.write(filtered_records, output_handle, "genbank")
                with open('cgview.log','a') as cgvf:        
                    subprocess.run(f'ruby /data/deploy/bio-elite/bio/script/cgview_builder_cli.rb -s {Pre}_prokka/main.gbk -o main.json -c /data/deploy/bio-elite/bio/script/CGview.yaml -n main',shell=True,stdout=cgvf,stderr=cgvf)
            else: # 没有质粒序列时的处理方式
                subprocess.run(f'ruby /data/deploy/bio-elite/bio/script/cgview_builder_cli.rb -s {Pre}_prokka/{Pre}.gbk -o main.json -c /data/deploy/bio-elite/bio/script/CGview.yaml -n {Pre}',shell=True)
            if not os.path.isdir(f'{Pre}_bin_genome_out'):
                os.makedirs(f'{Pre}_bin_genome_out')
            
            subprocess.run(f'cp {finalfa} {Pre}_bin_genome_out/{Pre}.fna',shell=True)          
            # 1.判断物种 2.确定对应物种的完整性和污染率
            with open('checkm2','w') as cmf:
                subprocess.run(f'/home/dell/miniconda3/bin/conda run --no-capture-output  -n cm2 checkm2 predict -i  {Pre}_bin_genome_out -x fna -o checkm2_out -t {threads} --force --database_path /data1/shanghai_pip/meta_genome/uniref100.KO.1.dmnd',shell=True,stdout=cmf,stderr=cmf)
                if os.path.isfile('checkm2_out/quality_report.tsv'):
                    checkmdb = pd.read_table('checkm2_out/quality_report.tsv')
                    checkmdb['样本名称'] = Pre
                    checkmdb.rename(columns={'Completeness':'完整性','Contamination':'污染率'},inplace=True)
                    checkmdb['物种名称'] = species
                    checkmdb[['样本名称','物种名称','污染率','完整性']].to_csv(f'{Pre}.checkm.tsv',sep='\t',index=False)
            subprocess.run(f'seqkit stat {Pre}.final.fasta -b -T -a > {Pre}.fasum.tsv',shell=True)  # seqkit统计结果输出到{Pre}.fasum.tsv
            fadb = pd.read_table(f'{Pre}.fasum.tsv')
            fadb['Sample'] = fadb['file'].str.replace('.final.fasta','')
            Assdict = {}
            Assdict['Contig数量'] = fadb['num_seqs'].tolist()[0]
            Assdict['N50长度'] = fadb['N50'].tolist()[0]
            Assdict['组装基因组长度'] = fadb['sum_len'].tolist()[0]
            Assdict['最长片段长度'] = fadb['max_len'].tolist()[0]
            Assdict['污染率'] = '-'
            Assdict['完整性'] = '-'
            Assdict['主基因组是否成环'] = '-'
            if os.path.isfile(f'{Pre}.checkm.tsv'):
                cdb = pd.read_table(f'{Pre}.checkm.tsv')
                Assdict['污染率'] = round(float(cdb['污染率'].tolist()[0]),2)
                Assdict['完整性'] = round(float(cdb['完整性'].tolist()[0]),2)
            Assdb = pd.DataFrame(Assdict,index=[0])
            Assdb['样本名称'] = Pre
            Assdb = Assdb[['样本名称','Contig数量','N50长度','最长片段长度','主基因组是否成环','污染率','完整性']]
            Assdb['主基因组是否成环'] = flyedb['是否成环'].tolist()[0]
            Assdb.to_csv(f'{Pre}.assemble.result.tsv',sep='\t',index=False)   # 以上代码均在asb_func的后半段出现 
            # genovi中的genovi工具对输入文件进行功能注释
            with open('genovi.log','w') as fv:
                try:  # 使用genovi 对 {Pre}_prokka/{Pre}.gbk进行注释，主要获取 COG 分类
                    subprocess.run(f'/home/dell/miniconda3/bin/conda run -n genovi genovi -i {Pre}_prokka/{Pre}.gbk -o {Pre}_genovi -s draft',shell=True,stdout=fv,stderr=fv)
                    cogdb = pd.read_table(f'{Pre}_genovi/{Pre}_genovi_COG_Classification.csv',sep=',',header=1)
                    cogdb.iloc[-2,:].to_csv('Cog_summary.tsv',sep='\t',header=False)
                except:
                    pass
            
        # 选出{Pre}.tsv中'基因'相关的行并进行长度统计，保存统计结果到{Pre}_gene_raw_sum.tsv
        tmpgenedb = pd.read_table(f'{Pre}_prokka/{Pre}.tsv')
        tmpgenedb = tmpgenedb[tmpgenedb.ftype == 'gene']    # 保留'gene'相关的数据
        tmpgenedb_dict = {}             
        for minlg in range(0,2000,100):
            tmpgenedb_dict[f'{minlg}-{minlg+100}'] =  sum(tmpgenedb.length_bp.between(minlg,minlg+100))
        tmpgenedb_dict['>2000'] = sum(tmpgenedb.length_bp >= 2000)  # 以100bp为区间，统计0~2000不同长度范围内的基因数量，> 2000的单独统计
        genedb = pd.DataFrame(tmpgenedb_dict,index=['Gene数量']).T
        genedb['范围'] = genedb.index
        genedb = genedb[['范围','Gene数量']]
        genedb.to_csv(f'{Pre}_gene_raw_sum.tsv',sep='\t',index=False)
        # 读取{Pre}.txt，行列转置后，生成{Pre}.genefun_summary.tsv
        gene_fundb = pd.read_table(f'{Pre}_prokka/{Pre}.txt',sep=':')  # 读取{Pre}.txt以':'分隔符
        gene_fundb.index = gene_fundb.organism                         # 将 organism(物种）列设为行索引,行列转置，生成{Pre}.genefun_summary.tsv
        gene_fundb = gene_fundb.drop('organism',axis=1).T
        gene_fundb.to_csv(f'{Pre}.genefun_summary.tsv',sep='\t',index=False)
        # 统计{Pre}.uniqgene.fasta （经过筛选的核心基因或其他具有特定功能的基因）的长度分布，生成{Pre}.uniqgene.tsv；计算基因长度分布，生成{Pre}_gene_uniq_sum.tsv。uniqgene主要用于宏基因组，本流程中不会启用
        if os.path.isfile(f'{Pre}.uniqgene.fasta'):
            subprocess.run(f'seqkit fx2tab -n -l {Pre}.uniqgene.fasta > {Pre}.uniqgene.tsv',shell=True)  # -n 表示提取序列名称，-l 表示提取序列长度
            tmpgeneudb = pd.read_table(f'{Pre}.uniqgene.tsv',names=['Genename','length'])
            tmpgeneudb_dict = {}
            for minlg in range(0,2000,100):
                tmpgeneudb_dict[f'{minlg}-{minlg+100}'] =  sum(tmpgeneudb.length.between(minlg,minlg+100))
            tmpgeneudb_dict['>2000'] = sum(tmpgeneudb.length >= 2000)    
            geneudb = pd.DataFrame(tmpgeneudb_dict,index=['Gene数量']).T
            geneudb['范围'] = geneudb.index
            geneudb = geneudb[['范围','Gene数量']]
            geneudb.to_csv(f'{Pre}_gene_uniq_sum.tsv',sep='\t',index=False)
        # Contig : CDS / GENE
        subprocess.run(f'grep CDS {Pre}_prokka/{Pre}.gff > {Pre}.CDS.gff',shell=True)  # 筛选出所有包含 'CDS' 的行
        with open(f'{Pre}.CDS.gff') as f:
            open(f'{Pre}.Contig_gene.tsv','w').write(f'Contig\tCDS\n')
            for line in f:
                line = line.strip().split('\t')
                if line[2] == 'CDS':
                    contig = line[0]
                    cdsid = line[8].split('locus_tag=')[-1].split(';')[0] # 从第9列（属性信息）中提取出 'locus_tag='（即基因的唯一标识符）
                    open(f'{Pre}.Contig_gene.tsv','a').write(f'{contig}\t{cdsid}\n')  # 将提取的 Contig 名称 和 CDS ID 保存到 Pre.Contig_gene.tsv 文件中           
        subprocess.run(f'''grep 'gene' {Pre}_prokka/{Pre}.gff > {Pre}.gene.gff''',shell=True)  # 筛选出所有包含 'gene' 的行
        with open(f'{Pre}.gene.gff') as f:
            open(f'{Pre}.gene.bed','w').write(f'')
            for line in f:
                line = line.strip().split('\t')
                if line[2] == 'gene':
                    contig = line[0]  
                    cdsid = line[8].split('locus_tag=')[-1].split(';')[0]  
                    start = line[3]
                    end = line[4]
                    open(f'{Pre}.gene.bed','a').write(f'{contig}\t{start}\t{end}\t{cdsid}\n')   #  contig名称、基因起始位置、终止位置、基因ID 追加写入输出文件 {Pre}.gene.bed 中
        subprocess.run(f'samtools faidx {Pre}_prokka/{Pre}.fna',shell=True)  # 构建索引
        subprocess.run(f'bedtools getfasta -fi {Pre}_prokka/{Pre}.fna -bed {Pre}.gene.bed -name > tmp_{Pre}.gene.fasta',shell=True) # 根据 Pre.gene.bed 文件中列出的基因区域提取相应的基因序列
        subprocess.run(f'''cut -d ':' -f1  tmp_{Pre}.gene.fasta > {Pre}.gene.fasta ''',shell=True) # 整理后输出为{Pre}.gene.fasta
        time.sleep(5)

#定义DrugFinder函数：三代数据使用 minimap2 和 cuteSV ，二代数据使用 bwa 和 delly ，进行比对和结构变异检测。
def DrugFinder(Pre,threads):
    # 针对三代测序数据，使用 minimap2 比对到参考基因组中，cuteSV 检测结构变异（SV），snpEff 对SV进行注释，IGV_new.py 脚本生成 IGV 兼容的可视化文件
    if os.path.isfile(f'{Pre}.final.fastq'):  
        subprocess.run(f'minimap2 -ax map-ont /data/deploy/TB_soft/ref/TB/ref.fa {Pre}.final.fastq -t {threads}|samtools sort -o 2.CuteSV/{Pre}.ref.sort.bam',shell=True,stdout=f,stderr=f)
        subprocess.run(f'samtools index 2.CuteSV/{Pre}.ref.sort.bam',shell=True,stdout=f,stderr=f)
        subprocess.run(f'cuteSV 2.CuteSV/{Pre}.ref.sort.bam /data/deploy/TB_soft/ref/TB/ref.fa 2.CuteSV/{Pre}.cuteSV.vcf  2.CuteSV/ --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --threads {threads} -s 10',shell=True,stdout=f,stderr=f)
        subprocess.run(f'/home/dell/biosoft/snpEff/scripts/snpEff ann -noLog -noStats -no-downstream -no-upstream -no-utr -c reference/snpeff.config -dataDir . ref 2.CuteSV/{Pre}.cuteSV.vcf > 2.CuteSV/{Pre}.anno.vcf',shell=True,stdout=f,stderr=f)
        subprocess.run(f'python /data/deploy/TB_soft/other_soft/IGV_js/IGV_new.py -r ref.fa -m 2.CuteSV/{Pre}.ref.sort.bam -o 2.CuteSV/ -s {Pre}',shell=True,stderr=f,stdout=f)
        subprocess.run(f'cp /data/deploy/TB_soft/other_soft/IGV_js/igv.min.js 2.CuteSV/',shell=True)
    # 针对二代测序数据，使用 bwa mem 比对到参考基因组中，delly 检测结构变异（SV），snpEff 对SV进行注释，IGV_new.py 脚本生成 IGV 兼容的可视化文件。生成2.CuteSV/{Pre}.delly.anno.vcf
    if os.path.isfile(f'{Pre}.R1.fastq.gz'):
        if os.path.isfile(f'{Pre}.R2.fastq.gz'):
            subprocess.run(f'bwa mem /data/deploy/TB_soft/ref/TB/ref.fa {Pre}.R1.fastq.gz {Pre}.R2.fastq.gz -t {threads}|samtools sort -o 2.CuteSV/{Pre}.ngs_ref.sort.bam',shell=True,stdout=f,stderr=f)
        else:
            subprocess.run(f'bwa mem /data/deploy/TB_soft/ref/TB/ref.fa {Pre}.R1.fastq.gz -t {threads}|samtools sort -o 2.CuteSV/{Pre}.ngs_ref.sort.bam',shell=True,stdout=f,stderr=f)
        subprocess.run(f'samtools index 2.CuteSV/{Pre}.ngs_ref.sort.bam',shell=True,stdout=f,stderr=f)
        subprocess.run(f'delly call -o 2.CuteSV/{Pre}.delly.bcf -s 3 -q 20 -g ref.fa 2.CuteSV/{Pre}.ngs_ref.sort.bam',shell=True)
        subprocess.run(f'bcftools view 2.CuteSV/{Pre}.delly.bcf > 2.CuteSV/{Pre}.delly.vcf',shell=True)
        subprocess.run(f'''bcftools view -i "FILTER='PASS' & GT='1/1'" 2.CuteSV/{Pre}.delly.vcf > 2.CuteSV/{Pre}.delly.filt.vcf''',shell=True)
        subprocess.run(f'/home/dell/biosoft/snpEff/scripts/snpEff ann -noLog -noStats -no-downstream -no-upstream -no-utr -c reference/snpeff.config -dataDir . ref  2.CuteSV/{Pre}.delly.filt.vcf > 2.CuteSV/{Pre}.delly.anno.vcf',shell=True)
    open('SV_ok','w').write('')

# 定义getinfo函数：细菌序列比对毒力和耐药基因
def getinfo(Pre,threads=10):
    #1.card 2.vfdb
    refdict = {'card':'/home/dell/miniconda3/envs/TB_ONT/db/card/sequences','vfdb':'/data/deploy/meta_genome/database/vfdb.fasta'}
    for db in ['card','vfdb']:
        ref = refdict.get(db)
        if not os.path.isfile(f'2.{db}.sorted.bam') or os.popen(f'samtools view 2.{db}.sorted.bam|wc -l').read().strip() == '0':
            if os.path.isfile(f'2.2.fastq'):
                subprocess.run(f'minimap2 -ax sr {ref} 2.1.fastq 2.2.fastq -t 10 |samtools sort -o 2.{db}.sorted.bam',shell=True)
            else:
                subprocess.run(f'minimap2 -ax sr {ref} 2.1.fastq -t 10 |samtools sort -o 2.{db}.sorted.bam',shell=True)
            subprocess.run(f'samtools index 2.{db}.sorted.bam',shell=True)
            subprocess.run(f'mosdepth -b1 {db} 2.{db}.sorted.bam -n -t {threads}',shell=True)
            subprocess.run(f'gunzip {db}.regions.bed.gz -f',shell=True)
            subprocess.run(f'samtools idxstat 2.{db}.sorted.bam > {db}.stat.tsv',shell=True)
        dbfile = pd.read_table(f'{db}.regions.bed',header=None)
        coninfo = pd.read_table(f'{db}.stat.tsv',header=None,usecols=[0,2])
        coninfo.columns = [0,'card_subreads']
        #method1 
        depdb = pd.DataFrame(dbfile.groupby(0).apply(lambda x:round(sum(x[3])/x.shape[0],2)).reset_index(name='card_dep'))
        covdb = pd.DataFrame(dbfile.groupby(0).apply(lambda x:round(sum(x[3]>0)/x.shape[0],2)).reset_index(name='card_cov'))
        #subdb = pd.DataFrame(afile.groupby(0).apply(lambda x:sum(x[3])/x.shape[0]).reset_index(name='card_subreads'))
        rawdb = depdb.merge(covdb,on=0).merge(coninfo,on=0)
        rawdb = rawdb.sort_values('card_subreads',ascending=False)
        if db == 'card':
            rawdb = rawdb.loc[(rawdb['card_cov']>=0.1) & (rawdb['card_subreads']>10),:]
        else:
            rawdb = rawdb.loc[(rawdb['card_cov']>=0.01) & (rawdb['card_subreads']>10),:].head(50)
        if db == 'card':
            metadb = pd.read_table('/data/deploy/meta_genome/database/aro_index.tsv',usecols=['Model Name','AMR Gene Family','Drug Class','Resistance Mechanism'])
            rawdb['Model Name'] = rawdb[0].str.split('~~~').str[1]
            rawdb = rawdb.merge(metadb,on = 'Model Name')
            rawdb.columns = ['片段名称','平均深度','覆盖率','支持序列数','Model','基因家族','耐药分类','耐药机制']
            rawdb['耐药基因'] = rawdb['片段名称'].str.split('~~~').str[1]
            rawdb[['耐药基因','平均深度','覆盖率','支持序列数','基因家族','耐药分类','耐药机制']].to_csv(f'2.card.tsv',sep='\t',index=False)
        else:
            metadb = pd.read_table('/data/deploy/meta_genome/database/VFs_meta.tsv',encoding='Windows-1252',usecols=['VFID','Bacteria','Function','Mechanism'])
            contigdb = pd.read_table('/data/deploy/meta_genome/database/vfdb.contig.tsv')
            rawdb = rawdb.merge(contigdb,left_on=0,right_on='Contig Name')
            rawdb = rawdb.merge(metadb,on='VFID')
            rawdb['毒力基因'] = rawdb[0].str.split('~~~').str[1]
            rawdb.columns = ['片段名称','平均深度','覆盖率','支持序列数','Contig','VFID','菌株','毒力功能','毒力机制','毒力基因']
            rawdb[['毒力基因','平均深度','覆盖率','支持序列数','VFID','菌株','毒力功能','毒力机制']].to_csv(f'2.vfdb.tsv',sep='\t',index=False)
# 定义getCovDep函数：没有组装时的抽序列,2D序列比对毒力和耐药基因
def getCovDep(Pre,Pre2):
    #1.提取细菌序列 2.组装？ 3.blast? 4.将reads比对到blast/直接比对card vfdb数据库 5.获取测序深度和覆盖度
    if not os.path.isfile('2.1.fastq'):
        level = 'D'
        krakenfile = f'{Pre2}.report.txt'
        tkid='2'
        tkidl = os.popen(f'grep {tkid} {krakenfile}|cut -f4').read().strip()
        taxlist1 = proc_kra1(krakenfile,tkid,level)
        taxlist1 = [int(i) for i in taxlist1]
        if os.path.isfile(f'{Pre}.R2.fastq.gz') and os.path.getsize(f'{Pre}.R2.fastq.gz') != 0:
            exreadsID1(taxlist1,f'{Pre2}.out.txt',f'{Pre}.R1.fastq.gz',f'{Pre}.R2.fastq.gz')
        else:
            exreadsID1(taxlist1,f'{Pre2}.out.txt',f'{Pre}.R1.fastq.gz',0)
    getinfo(Pre)
# 定义kk2函数：对二代和/或三代数据进行物种鉴定、整理表格及物种丰度可视化
def kk2(inf,fq1,fq2,threads,Pre):   
    with open('kk2.log','w') as kkf:
        if inf:         #对三代测序数据进行物种鉴定，并提取丰度最高的物种到ONTSpe。待完善
            if not os.path.isfile(f'{Pre}.list.txt'):
                if not os.path.isfile(f'{Pre}.report.txt'):
                    subprocess.run(f'kraken2 --db {Krdb} --threads {threads} --output {Pre}.out.txt --report {Pre}.report.txt {inf}',shell=True,stdout=kkf,stderr=kkf)
                subprocess.run(f'bracken -d {Krdb} -o {Pre}.bracken1.txt -w {Pre}.bracken2.txt -l S -t 10  -i {Pre}.report.txt',shell=True,stdout=kkf,stderr=kkf) 
                subprocess.run(f'bracken -d {Krdb} -o {Pre}.bracken1.txt -w {Pre}.bracken2.txt -l S3 -t 10  -i {Pre}.report.txt',shell=True,stdout=kkf,stderr=kkf) # {Pre}.bracken1.txt：丰度重估结果 {Pre}.bracken2.txt：详细的丰度信息，包括不同分类层级的读数数量
            tmpfile = pd.read_table(f'{Pre}.bracken1.txt')
            ONTSpe = tmpfile.name.tolist()[0]
            try:
                getCovDep(Pre,Pre)
            except:
                pass
            
        if fq1 and fq2:   #对二代测序数据进行物种鉴定，并提取丰度最高的物种到ngsSpe
            if not os.path.isfile(f'{Pre}_2.list.txt') and not os.path.isfile(f'{Pre}_2.report.txt'):
                subprocess.run(f'kraken2 --db {Krdb} --threads {threads} --output {Pre}_2.out.txt --report {Pre}_2.report.txt {fq1} {fq2}',shell=True,stdout=kkf,stderr=kkf)
                subprocess.run(f'bracken -d {Krdb} -o {Pre}_2.bracken1.txt -w {Pre}_2.bracken2.txt -l S -t 10  -i {Pre}_2.report.txt',shell=True,stdout=kkf,stderr=kkf)
                testbrkdb = pd.read_table(f'{Pre}_2.report.txt',header=None)
                if 'S4' in testbrkdb[3]:
                    subprocess.run(f'bracken -d {Krdb} -o {Pre}_2_Sub.bracken1.txt -w {Pre}_2_Sub.bracken2.txt -l S3 -t 10  -i {Pre}_2.report.txt',shell=True,stdout=kkf,stderr=kkf)
                else:
                    if 'S3' in testbrkdb[3]:
                        subprocess.run(f'bracken -d {Krdb} -o {Pre}_2_Sub.bracken1.txt -w {Pre}_2_Sub.bracken2.txt -l S2 -t 10  -i {Pre}_2.report.txt',shell=True,stdout=kkf,stderr=kkf)
                    else:
                        subprocess.run(f'bracken -d {Krdb} -o {Pre}_2_Sub.bracken1.txt -w {Pre}_2_Sub.bracken2.txt -l S1 -t 10  -i {Pre}_2.report.txt',shell=True,stdout=kkf,stderr=kkf)
            tmpfile2 = pd.read_table(f'{Pre}_2.bracken1.txt')
            ngsSpe = tmpfile2.name.tolist()[0]
            getCovDep(Pre,f'{Pre}_2')
            
        elif fq1:        #对二代单端测序数据进行物种鉴定，并提取丰度最高的物种到ngsSpe
            if not os.path.isfile(f'{Pre}_2.list.txt') and not os.path.isfile(f'{Pre}_2.report.txt'):
                subprocess.run(f'kraken2 --db {Krdb} --threads {threads} --output {Pre}_2.out.txt --report {Pre}_2.report.txt {fq1}',shell=True,stdout=kkf,stderr=kkf)
                subprocess.run(f'bracken -d {Krdb} -o {Pre}_2.bracken1.txt -w {Pre}_2.bracken2.txt -l S -t 10  -i {Pre}_2.report.txt',shell=True,stdout=kkf,stderr=kkf)
                testbrkdb = pd.read_table(f'{Pre}_2.report.txt',header=None)
                if 'S4' in testbrkdb[3]:
                    subprocess.run(f'bracken -d {Krdb} -o {Pre}_2_Sub.bracken1.txt -w {Pre}_2_Sub.bracken2.txt -l S3 -t 10  -i {Pre}_2.report.txt',shell=True,stdout=kkf,stderr=kkf)
                else:
                    if 'S3' in testbrkdb[3]:
                        subprocess.run(f'bracken -d {Krdb} -o {Pre}_2_Sub.bracken1.txt -w {Pre}_2_Sub.bracken2.txt -l S2 -t 10  -i {Pre}_2.report.txt',shell=True,stdout=kkf,stderr=kkf)
                    else:
                        subprocess.run(f'bracken -d {Krdb} -o {Pre}_2_Sub.bracken1.txt -w {Pre}_2_Sub.bracken2.txt -l S1 -t 10  -i {Pre}_2.report.txt',shell=True,stdout=kkf,stderr=kkf)
            tmpfile2 = pd.read_table(f'{Pre}_2.bracken1.txt')
            ngsSpe = tmpfile2.name.tolist()[0]
            try:
                getCovDep(Pre,f'{Pre}_2')
            except:
                pass
     

        if 'ngsSpe' in dir() and 'ONTSpe' in dir():  # 二代和三代都有时，比较最高丰度物种是否一致
            if ngsSpe != ONTSpe:
                 raise Exception(f'二三代不是同一菌种测序数据')
        if 'ONTSpe' in dir():
            tmpfile = tmpfile[['name','taxonomy_id','taxonomy_lvl','new_est_reads','fraction_total_reads']]
        else:
            tmpfile = tmpfile2[['name','taxonomy_id','taxonomy_lvl','new_est_reads','fraction_total_reads']]
        tmpfile.rename(columns={'name':'物种','taxonomy_id':'taxid','taxonomy_lvl':'水平','new_est_reads':'序列数量','fraction_total_reads':'相对丰度'},inplace=True)
        sanofile = pd.read_table('/data/Ref/Meta_anno/taxa_info_20210508.txt')
        sanofile['taxid'] = sanofile['taxid'].astype('int')  # taxid 列的数据类型转换为整数类型
        kran_dic = {}                                                               
        kran_dic1 = {}                                                               
        
        # 定义一个内部函数kran_summ：物种序列＞2条时才会被选入，进行一系列的可视化文件准备
        def kran_summ(ofn1,kfile):  # {Pre}.list.txt/{Pre}.list.txt , {Pre}.bracken2.txt
            with open(kfile) as f:
                tad,tap,tac,tap,taf,tag = '-'*6
                for line in f:
                    line = line.strip().split('\t')
                    line[5] = line[5].strip()  # 种类名称
                    prop = line[0].strip()     # 丰度
                    snum = line[1].strip()     # 归属于该分类单元的总读数数量
                    taxid = line[4].strip()    # 分类 ID
                    if line[5] == 'unclassified':
                        #-----[物种英文，匹配到的reads数量,百分比]，将未分类的物种信息写入文件
                        #f1.write(f'{line[5]}\t{line[1]}\t{line[0]}\t-\t-\t-\t-\t-\t-\t-\t-\n')
                        #f2.write(f'{line[5]}\t{line[1]}\t{line[0]}\t-\t-\t-\t-\t-\t-\t-\t-\n')
                        pass
                    else:
                        if line[3] == 'D': tad =line[5]  # tad 保存界（Domain）的分类名称
                        elif line[3] =='P': tap = line[5] # tap：保存门（Phylum）的分类名称
                        elif line[3] =='C': tac = line[5] # tac：保存纲（Class）的分类名称
                        elif line[3] =='O': tao = line[5] # tao：保存目（Order）的分类名称
                        elif line[3] =='F': taf =line[5]  # taf：保存科（Family）的分类名称
                        elif line[3] =='G': tag = line[5] # tag：保存属（Genus）的分类名称
                        elif line[3] =='S' and int(line[1]) > 2:   # 种（Species）的序列数量＞2才会被选入
                            kran_dic[line[5]]= {'D': tad,'P':tap,'C':tac,'O':tao,'F':taf,'G':tag,'S':line[5],'比例':prop,'序列数量':snum,'taxid':taxid}
                        elif line[3] in ['S1','S2','S3'] and int(line[1]) > 2:   # 亚种（SubSpecies）的序列数量＞2才会被选入
                            kran_dic1[line[5]]= {'D': tad,'P':tap,'C':tac,'O':tao,'F':taf,'G':tag,'亚种':line[5],'比例':prop,'序列数量':snum,'taxid':taxid}
            tmpdb = pd.DataFrame(kran_dic).T
            tmpdb.sort_values('比例',inplace=True,ascending=False)  # 按照'比例'列的值进行从大到小排序
            tmpdb['taxid'] = tmpdb['taxid'].astype('int')
            tmpdb = tmpdb.merge(sanofile,on='taxid',how='left')   # 依据 taxid 列进行合并，保留 tmpdb 中的所有行。增加致病菌信息
            tmpdb.rename(columns={'D':'界','P':'门','C':'纲','O':'目','F':'科','G':'属','S':'种'},inplace=True)
            tmpdb.to_csv(ofn1,sep='\t',index=False)     # 输出为{Pre}.list.txt  
            tmpdb1 = pd.DataFrame(kran_dic1).T
            if tmpdb1.shape[0] >0 :
                tmpdb1.sort_values('比例',inplace=True,ascending=False)  # 按照'比例'列的值进行从大到小排序
                tmpdb1['taxid'] = tmpdb1['taxid'].astype('int')
                tmpdb1 = tmpdb1.merge(sanofile,on='taxid',how='left')   # 依据 taxid 列进行合并，保留 tmpdb 中的所有行。增加致病菌信息
                tmpdb1.rename(columns={'D':'界','P':'门','C':'纲','O':'目','F':'科','G':'属','S':'种'},inplace=True)
                tmpdb1.to_csv(ofn1,sep='\t',index=False)     # 输出为{Pre}.list.txt

        
        # 传入{Pre}*.bracken2.txt，调用kran_summ、3_kreport2krona.py等进行可视化，split_kraken统计
        if os.path.isfile(f'{Pre}.bracken2.txt'):
            kran_summ(f'{Pre}.list.txt',f'{Pre}.bracken2.txt') # 结果文件{Pre}.list.txt
            subprocess.run(f'/data1/shanghai_pip/meta_genome/soft/3_kreport2krona.py -r {Pre}.bracken2.txt -o {Pre}.krona.txt',shell=True,stdout=kkf,stderr=kkf) # {Pre}.bracken2.txt 结果转换为 {Pre}.krona.txt
            subprocess.run(f'ktImportText {Pre}.krona.txt -o {Pre}.krona.html',shell=True)  # 使用ktImportText生成{Pre}.krona.html
            #split_kraken(Pre,'./')  # 该功能在R中'Spedb <- Spedb[,c('taxname','中文名','序列数量','比例','致病源性','可能引起的疾病','G')]',有更简便的解决方案 
        if os.path.isfile(f'{Pre}_2.bracken2.txt'):
            kran_summ(f'{Pre}_2.list.txt',f'{Pre}_2.bracken2.txt')
            kran_summ(f'{Pre}_2.list2.txt',f'{Pre}_2_Sub.bracken2.txt')
            subprocess.run(f'/data1/shanghai_pip/meta_genome/soft/3_kreport2krona.py -r {Pre}_2.bracken2.txt -o {Pre}_2.krona.txt',shell=True,stdout=kkf,stderr=kkf)
            subprocess.run(f'ktImportText {Pre}_2.krona.txt -o {Pre}_2.krona.html',shell=True,stdout=kkf,stderr=kkf)
            #split_kraken(Pre,'./')
    tmpdb = pd.read_table(f'{Pre}_2.list.txt')
    #---物种数---
    Summarydb_K = pd.DataFrame(tmpdb['界'].value_counts()).reset_index()
    Summarydb_P = pd.DataFrame(tmpdb['危害程度等级'].value_counts()).reset_index()
    print(Summarydb_K)
    print(Summarydb_P)
    #---序列数---
    Summarydb_Kreads  = tmpdb.groupby('界').apply(lambda x:x['序列数量'].sum()).reset_index(name='序列数')
    print(Summarydb_Kreads)
    plantdb = pd.read_csv('/data/Ref/Meta_anno/All.sort.csv',usecols=['taxonId','type'])
    tmpdb = tmpdb.merge(plantdb,left_on = 'taxid',right_on='taxonId',how='left').fillna('-')
    tmpdb.drop('taxonId',inplace=True,axis=1)
    tmpdb.to_csv(f'{Pre}.anno.tsv',sep='\t',index=False)
    #---汇总----
    if os.path.isfile(f'{Pre}_2.report.txt'):
        rawrpdb = pd.read_table(f'{Pre}_2.report.txt',header=None)
    else:
        rawrpdb = pd.read_table(f'{Pre}.report.txt',header=None)
    rawrpdb[5] = rawrpdb[5].str.strip()
    rawrpdb1 = rawrpdb.copy()
    rawrpdb = rawrpdb.loc[rawrpdb[5].isin(['unclassified','root'])]
    summarydict = {}
    if 'unclassified' in rawrpdb[5].tolist():
        ureads =  rawrpdb.loc[rawrpdb[5]=='unclassified',1].tolist()[0]
    else:
        ureads = 0
    if 'root' in rawrpdb[5].tolist():
        creads =  rawrpdb.loc[rawrpdb[5]=='root',1].tolist()[0]
    else:
        creads = 0
    if os.path.isfile(f'{Pre}_2.list2.txt'):
        Sdb = pd.read_table(f'{Pre}_2.list.txt')
        Subdb = pd.read_table(f'{Pre}_2.list2.txt')
    else:
        Sdb = pd.read_table(f'{Pre}.list.txt')
        Subdb = pd.read_table(f'{Pre}.list2.txt')
    Prodb = pd.read_table('/data/Ref/Meta_anno/AllSpeProid_rank.txt',header=None,names=['Taxid','Type'])
    Prodb = Prodb.loc[Prodb['Type']=='Species',:]
    Prodb1 = rawrpdb1.loc[rawrpdb1[5].isin(Prodb['Taxid'].tolist()),:]
    Wormdb = pd.read_table('/data/Ref/Meta_anno/wormbase.tsv')
    Wormdb1 = rawrpdb1.loc[rawrpdb1[5].isin(Wormdb['taxid'].tolist()),:]
    summarydict1 = {}
    summarydict['有效序列'] = rawrpdb[1].sum()
    summarydict['未识别序列'] = ureads
    summarydict['未识别序列比例'] = f'{round(ureads/(ureads+creads)*100,2)}%'
    summarydict['可识别序列'] = creads
    summarydict['可识别序列比例'] = f'{round(creads/(ureads+creads)*100,2)}%'
    summarydict['校正识别序列数(种)'] = Sdb['序列数量'].sum()
    summarydict['校正识别序列数(亚种)'] = Subdb['序列数量'].sum()
    summarydict['细菌'] = Sdb.loc[Sdb['界']=='Bacteria','序列数量'].sum()
    summarydict['病毒'] = Sdb.loc[Sdb['界']=='Viruses','序列数量'].sum()
    summarydict1['细菌'] = Sdb.loc[Sdb['界']=='Bacteria',].shape[0]
    summarydict1['病毒'] = Sdb.loc[Sdb['界']=='Viruses',].shape[0]
    if '真菌' in rawrpdb[5].tolist()[0]:
        summarydict['真菌'] = rawrpdb1.loc[rawrpdb1[5]=='Fungi',1].tolist()[0]
        summarydict1['真菌'] = rawrpdb1.loc[rawrpdb1[5]=='Fungi',].shape[0]
    else:
        summarydict['真菌'] = 0
        summarydict1['真菌'] = 0
    if '古菌' in rawrpdb[5].tolist()[0]:
        summarydict['古菌'] = rawrpdb1.loc[rawrpdb1[5]=='Archaea',1].tolist()[0]
    else:
        summarydict['古菌'] = 0
    summarydict['原生动物'] = Prodb1[2].sum()
    summarydict['寄生虫'] = Wormdb1[2].sum()
    summarydict1['寄生虫'] = Wormdb1.shape[0]
    #--宿主
    hostdb = pd.read_table('summary.tsv')
    hostdb1 = pd.read_table('R1_Fastqc.tsv')
    hostrate = round((hostdb['num_seqs'][0]-hostdb1['总序列数'][1])/hostdb['num_seqs'][0],4)*100
    summarydict1['宿主'] = hostrate
    pd.DataFrame(summarydict,index=[0]).to_csv(f'Summary_kraken.csv')
    pd.DataFrame(summarydict1,index=[0]).to_csv(f'Summary_kraken1.csv')
    open('kk2_ok','w').write('')   # 创建一个名为 kk2_ok 的空文件，作为处理完成的标志
 
# 定义getvfID函数：从tmpdb 数据框的'产物'列中提取与'VF'相关的 ID
def getvfID(tmpdb):
    matches = re.findall(r'VF\d+', tmpdb['产物'])  # matches = tmpdb['产物'] 列中正则匹配的字符串
    if len(matches) >0:
        if len(matches) > 1:
            newID = '|'.join([i for i in matches])
        else:
            newID = matches[0]
    else:
        newID = '-'
    return newID     # 例，newID = ['vfID1'|'vfID2']
# 定义rgi_fun函数：通过RGI进行抗药性基因的鉴定
def rgi_fun(Pre):
    subprocess.run(f'/home/dell/miniconda3/bin/conda run -n RGI_new rgi main -i {Pre}.final.fasta -o {Pre}.rgi --clean --include_loose --low_quality -g PYRODIGAL',shell=True) # 调用 RGI 工具的主程序来进行抗药性基因的预测
    rgidb = pd.read_table(f'{Pre}.rgi.txt')
    rgidb = rgidb[['Contig','Start','Stop','Orientation','Cut_Off','Pass_Bitscore','Best_Hit_ARO','Best_Identities','Model_type','SNPs_in_Best_Hit_ARO','AMR Gene Family','Drug Class']]
    rgidb['序列名称'] = rgidb['Contig'].str.split('_').str[:2].str.join('_')
    rgidb.rename(columns={'Start':'序列起始','Stop':'序列终止','Orientation':'正负链','Cut_Off':'过滤标准','Pass_Bitscore':'比对得分','Best_Hit_ARO':'耐药基因','Best_Identities':'一致性%','Model_type':'基因数据库','SNPs_in_Best_Hit_ARO':'基因突变','AMR Gene Family':'耐药基因家族','Drug Class':'药物类别'},inplace=True)
    rgidb.drop('Contig',axis=1,inplace=True)
    rgidb = rgidb[['序列名称','序列起始','序列终止','正负链','过滤标准','比对得分','耐药基因','一致性%','基因数据库','基因突变','耐药基因家族','药物类别']]
    rgidb.sort_values('比对得分',ascending=False,inplace=True)  # 比对得分 对数据进行降序排序
    rgidb.to_csv(f'{Pre}.rgi.tsv',sep='\t',index=False)       
    rgidb.to_csv(f'{Pre}.rgi.bed',sep='\t',index=False,header=False)         # 结果分别保存为{Pre}.rgi.tsv和{Pre}.rgi.bed  
# 定义outfun函数：assem_vfdr中调用
def outfun(Pre,x,typeF):
    tmpdict = {}
    ts = x['start'].min() # 获取基因的最小起始位置
    te = x['end'].max()   # 获取基因的最大终止位置
    ofname = x['GeneName'].tolist()[0] # ofname = 'GeneName'的第一个元素
    Spename = x['Species'].tolist()[0]
    x['start'] = x.reset_index().index+1 
    x['end'] = x.reset_index().index+2  # 将'start'和'end'列的值更新为,从1开始的索引值 和 从2开始的索引值
    x[['Chrom','start','end','Depth']].to_csv(f'geneDepth/{ofname}_{typeF}.tsv',sep='\t',header=False,index=False)  # 提取'Chrom','start','end','Depth'列输出到geneDepth/{ofname}_{typeF}.tsv
    open(f'geneDepth/{ofname}_{typeF}.bed','w').write(f'''{x['Chrom'].tolist()[0]}\t{ts}\t{te}\t{Pre}_{ofname}\t{Pre}_{ofname}\t{x['strand'].tolist()[0]}''') # 写入相关行，生成geneDepth/{ofname}_{typeF}.bed
    #  使用bedtools getfasta，以bed文件为参考，提取fasta基因序列。-name：保留基因的名称作为序列的 ID，-s：处理正负链，输出为
    subprocess.run(f'''bedtools getfasta -fi {wkdir}/{Pre}/{Pre}.final.fasta -bed "geneDepth/{ofname}_{typeF}.bed" -name -s > "geneDepth/{ofname}_{typeF}.fasta" ''',shell=True)
    tmpdict['片段名称'] = x['Chrom'].tolist()[0]
    tmpdict['物种名称'] = Spename
    tmpdict['起始位置'] = x['start'].min()
    tmpdict['终止位置'] = x['end'].max()
    tmpdict['覆盖度(>0)%'] =  round(x[x['Depth']>0].shape[0]/x.shape[0],4)*100
    tmpdict['覆盖度(>10)%'] =  round(x[x['Depth']>10].shape[0]/x.shape[0],4)*100
    tmpdict['覆盖度(>100)%'] =  round(x[x['Depth']>100].shape[0]/x.shape[0],4)*100
    tmpdict['平均深度'] = round(x['Depth'].mean(),2)
    tmpdict['最低深度'] = x['Depth'].min()
    tmpdict['最高深度'] = x['Depth'].max()
    return pd.DataFrame(tmpdict,index=[0]).round(2)
# 定义hAMRCom(函数：
def hAMRCom(Pre):
    subprocess.run(f'abricate --db card {Pre}.final.fasta > {Pre}.abricate.tsv',shell=True)
    subprocess.run(f'/home/dell/miniconda3/bin/conda run -n hamronization hamronize abricate {Pre}.abricate.tsv --format tsv --analysis_software_version 1.0.1 --reference_database_version 20250207 > {Pre}.hamr.abricate.tsv',shell=True)
    subprocess.run(f'/home/dell/miniconda3/bin/conda run -n RGI rgi main -i {Pre}.final.fasta -o {Pre}.rgi --clean --include_loose',shell=True)
    subprocess.run(f'/home/dell/miniconda3/bin/conda run -n hamronization hamronize rgi {Pre}.rgi.txt  --format tsv --analysis_software_version 6.0.3 --reference_database_version 20250207 --input_file_name {Pre}.final > {Pre}.hamr.rgi.tsv',shell=True)
    subprocess.run(f'/home/dell/miniconda3/bin/conda run -n hamronization run_resfinder.py -ifa {Pre}.final.fasta -o {Pre}_resfinder -acq',shell=True)
    subprocess.run(f'/home/dell/miniconda3/bin/conda run -n hamronization hamronize resfinder {Pre}_resfinder/ResFinder_results_tab.txt  --format tsv --analysis_software_version 4.6.0 --reference_database_version 20250207 --input_file_name {Pre}.final > {Pre}.hamr.resfinder.tsv',shell=True)
    subprocess.run(f'/home/dell/miniconda3/bin/conda run -n hamronization amrfinder -n {Pre}.final.fasta -o {Pre}.amrfinder.tsv -t 10',shell=True)
    subprocess.run(f'/home/dell/miniconda3/bin/conda run -n hamronization hamronize amrfinderplus {Pre}.amrfinder.tsv  --format tsv --analysis_software_version 4.0.15 --reference_database_version 20250207 --input_file_name {Pre}.final > {Pre}.hamr.amrfinderplus.tsv',shell=True)
    subprocess.run(f'hamronize summarize *.hamr*.tsv -t interactive -o {Pre}.hamr.html',shell=True)
# 定义assem_vfdr函数：
def assem_vfdr(Pre,inty='fastq'):
    print('组装与耐药毒力分析开始')
    if inty == 'fastq':
        if os.path.isfile(f'{Pre}_1.regions.bed'):   # 三代数据以1bp窗口进行比对的文件
            subprocess.run(f'ln -s {Pre}_1.regions.bed {Pre}_assem.regions.bed',shell=True)
        else:                                        # 二代数据以1bp窗口进行比对的文件  
            subprocess.run(f'ln -s {Pre}_ngs_1.regions.bed {Pre}_assem.regions.bed',shell=True)
    #---20260104 修改物种鉴定从kraken2到gtdbtk
    if not os.path.isfile(f"{Pre}_assem.kraken2.txt"):
        subprocess.run(f'kraken2 --db {Krdb} --threads 10 --output {Pre}_assem.txt --report {Pre}_assem.kraken2.txt {Pre}.final.fasta',shell=True) # 运行kraken2对{Pre}.final.fasta进行物种鉴定
    subprocess.run(f'abricate  --minid 50 --mincov 50 --db card --threads 10 --quiet  {Pre}.final.fasta > Assem_abricate_CARD.txt',shell=True)
    subprocess.run(f'abricate  --minid 50 --mincov 50 --db vfdb --threads 10 --quiet  {Pre}.final.fasta > Assem_abricate_VFDB.txt',shell=True) # 运行abricate再次进行card和vfdb注释，无门槛
    asvfdb = pd.read_table(f'Assem_abricate_VFDB.txt')
    asdrdb = pd.read_table(f'Assem_abricate_CARD.txt')
    ask2db = pd.read_table(f'{Pre}_assem.txt',names=['type1','contig','taxid','readsL','taxidlist'])  # type1：C 分类成功，U 未分类；contig ；taxid；readsL：contig长度；taxidlist：物种分类的详细路径
    ask2db1 = pd.read_table(f'{Pre}_assem.kraken2.txt',names=['Abundance','ReadsC','ReadsO','taxid','SciName']) # Abundance：相对丰度百分比；ReadsC：能分类到的全部序列数；ReadsO：只能到这个层级的序列数
    ask2db1['Name'] = ask2db1['SciName'].str.strip()
    ask2db = ask2db.merge(ask2db1,on='taxid',how='left') # ask2db = {Pre}_assem.txt和{Pre}_assem.kraken2.txt按照'taxid'进行合并，左链接即保留{Pre}_assem.txt
    asdrdb = asdrdb.merge(ask2db,left_on='SEQUENCE',right_on='contig',how='left') # asdrdb = 按照asdrdb的'SEQUENCE'列和ask2db的'contig'列进行合并
    asvfdb = asvfdb.merge(ask2db,left_on='SEQUENCE',right_on='contig',how='left') # asvfdb = 按照asvfdb的'SEQUENCE'列和ask2db的'contig'列进行合并
    asvfdb = asvfdb.fillna('-')  # 所有的 NaN 值（缺失值）替换为 '-'
    asvfdb = asvfdb[['SEQUENCE','Name','taxid','GENE','%COVERAGE','%IDENTITY','PRODUCT','START', 'END', 'STRAND']]
    asvfdb.rename(columns={'SEQUENCE':'Contig名称','START':'起始碱基','END':'终止碱基','STRAND':'正负链','GENE':'基因名称','%COVERAGE':'覆盖度%','%IDENTITY':'一致性%','PRODUCT':'产物'},inplace=True)
    if asvfdb.shape[0] > 0:
        asvfdb['VFID'] = asvfdb.apply(lambda x:getvfID(x),axis=1) # 新增VFID列 = 对每一行调用函数 getvfID()，从产物列提取ID
    else:
        asvfdb['VFID'] = '-'
    
    asvfdb = asvfdb.merge(vfmeta,on='VFID',how='left')
    asvfdb.rename(columns={'Name':'物种名称','VF_Name':'VF名称','VF_FullName':'VF全称','Bacteria':'物种来源','VFcategory':'VF分类','Characteristics':'特征','Structure':'结构','Function':'功能','Mechanism':'机制','Reference':'文献来源'},inplace=True)
    asvfdb = asvfdb[['Contig名称','物种名称','taxid','基因名称','覆盖度%','一致性%','产物','VF分类','VF名称','起始碱基','终止碱基','正负链']]  # 添加列名并选取特定的列
    asdrdb = asdrdb[['SEQUENCE','Name','taxid','GENE','%COVERAGE','%IDENTITY','PRODUCT','RESISTANCE','START', 'END', 'STRAND']]
    asdrdb.rename(columns={'SEQUENCE':'Contig名称','START':'起始碱基','END':'终止碱基','STRAND':'正负链','GENE':'基因名称','%COVERAGE':'覆盖度%','%IDENTITY':'一致性%','PRODUCT':'产物','RESISTANCE':'耐药药物','Name':'物种名称'},inplace=True)
    asvfdb['基因名称'] = asvfdb['基因名称'].str.replace("'",'').str.replace('/','_').str.replace(' ','').replace('(','_').replace(')','_')
    asdrdb['基因名称'] = asdrdb['基因名称'].str.replace("'",'').str.replace('/','_').str.replace(' ','').replace('(','_').replace(')','_')
    asvfdb.to_csv('Assem_abricate_VFDB.tsv',sep='\t',index=False)
    asdrdb.to_csv('Assem_abricate_CARD.tsv',sep='\t',index=False)  # 整理后，输出Assem_abricate_VFDB.tsv和Assem_abricate_CARD.tsv
    
    if not os.path.isdir('geneDepth'):
        os.makedirs('geneDepth')
    if asvfdb.shape[0] > 0:
        asvfdb[['Contig名称','起始碱基','终止碱基','基因名称','正负链','物种名称']].to_csv('Assem_abricate_VFDB.bed',header=False,index=False,sep='\t')
        if inty == 'fastq':
            subprocess.run(f'bedtools intersect -a Assem_abricate_VFDB.bed -b {Pre}_assem.regions.bed -wb > Assem_abricate_VFDB.depth.bed',shell=True) # Assem_abricate_VFDB.depth.bed = Assem_abricate_VFDB.bed和{Pre}_assem.regions.bed文件的交集  
            asvfddb = pd.read_table('Assem_abricate_VFDB.depth.bed',header=None,names=['Chrom','start','end','GeneName','strand','Species','c1','s1','e1','Depth'])  # 添加列名
            asvfddb = asvfddb.groupby('GeneName').apply(lambda x:outfun(Pre,x,'vfdb')).reset_index(level=0) # 按 GeneName 进行分组，再调用outfun函数。reset_index(level=0)：重设索引，以便将 GeneName 作为列而不是索引
            asvfddb.rename(columns={'GeneName':'基因名称'},inplace=True)
            asvfddb.to_csv(f'VFDB_summary.tsv',sep='\t',index=False)    # outfun函数的统计结果，合并输出为VFDB_summary.tsv
    else:
        open('VFDB_summary.tsv','w').write(f'基因名称\t片段名称\t起始位置\t终止位置\t覆盖度(>0)%\t覆盖度(>10)%\t覆盖度(>100)%\t平均深度\t最低深度\t最高深度\n')
    
    if asdrdb.shape[0] > 0:
        if inty == 'fastq':
            asdrdb[['Contig名称','起始碱基','终止碱基','基因名称','正负链','物种名称']].to_csv('Assem_abricate_CARD.bed',header=False,index=False,sep='\t')
            subprocess.run(f'bedtools intersect -a Assem_abricate_CARD.bed -b {Pre}_assem.regions.bed -wb > Assem_abricate_CARD.depth.bed',shell=True)
            asdrddb = pd.read_table('Assem_abricate_CARD.depth.bed',header=None,names=['Chrom','start','end','GeneName','strand','Species','c1','s1','e1','Depth'])
            asdrddb = asdrddb.groupby('GeneName').apply(lambda x:outfun(Pre,x,'card')).reset_index(level=0)
            asdrddb.rename(columns={'GeneName':'基因名称'},inplace=True)
            asdrddb.to_csv(f'CARD_summary.tsv',sep='\t',index=False) # outfun函数的统计结果，输出为CARD_summary.tsv
    else:
        open('CARD_summary.tsv','w').write(f'基因名称\t片段名称\t起始位置\t终止位置\t覆盖度(>0)%\t覆盖度(>10)%\t覆盖度(>100)%\t平均深度\t最低深度\t最高深度\n')
    
    #with open('split.log','w') as spf:
    #    subprocess.run(f'''seqkit split -i --by-id-prefix '' -O  ass_split {Pre}.final''',shell=True,stdout=spf,stderr=spf)
    #    subprocess.run(f'seqkit stat -T ass_split/*.fasta > assembly_info.txt',shell=True,stdout=spf,stderr=spf)  
    assdb = pd.read_table(f'flye_output/assembly_info.txt')  # 序列的统计信息表
    assdb = assdb.merge(ask2db,how='left',left_on='序列名称',right_on='contig')
    assdb = assdb[['序列名称','序列长度','平均深度','是否成环','基因组/质粒','质粒分型','taxid','Name']]
    assdb['毒力基因数量'] = assdb.apply(lambda x:asvfdb[asvfdb['Contig名称']==x['序列名称']].shape[0],axis=1) # 查找每条contig对应的'毒力基因'数量
    assdb['毒力基因'] = assdb.apply(lambda x:','.join(asvfdb.loc[asvfdb['Contig名称']==x['序列名称'],'基因名称'].tolist()) if asvfdb.loc[asvfdb['Contig名称']==x['序列名称']].shape[0] > 0 else '-' ,axis=1) # 查找与该序列片段对应的所有毒力基因名称，以逗号分隔
    assdb['耐药基因数量'] = assdb.apply(lambda x:asdrdb[asdrdb['Contig名称']==x['序列名称']].shape[0],axis=1) # 查找每条contig对应的'耐药基因'数量
    assdb['耐药基因'] = assdb.apply(lambda x:','.join(asdrdb.loc[asdrdb['Contig名称']==x['序列名称'],'基因名称'].tolist()) if asdrdb.loc[asdrdb['Contig名称']==x['序列名称']].shape[0] > 0 else '-' ,axis=1)
    assdb = assdb.rename(columns={'contig':'Contig名称','length':'片段长度','Name':'物种名称'}) # 
    assdb.to_csv('Assem_info1.tsv',sep='\t',index=False)
    assdb[['序列名称','序列长度','平均深度','是否成环','基因组/质粒','质粒分型','物种名称','毒力基因数量','毒力基因','耐药基因数量','耐药基因']].to_csv('Assem_info.tsv',sep='\t',index=False) # 输出为Assem_info.tsv
    #rgi_fun(Pre)
    #20250207
    #hAMRCom(Pre)
   
# 定义VFDR函数：耐药与毒力基因
def VFDR(Pre,threads,intype='fastq'):
    with open('vfdr.log','w') as f1:
        # VFDB毒力基因
        subprocess.run(f'abricate --db vfdb {Pre}.final.fasta --threads {threads} --minid 50 --mincov 50 >{Pre}.vfdb.tsv',shell=True,stdout=f1,stderr=f1) # abricate进行VFDB注释,最小相似性和最小覆盖度阈值均为 50%,结果保存到 {Pre}.vfdb.tsv
        vfdb = pd.read_table(f'{Pre}.vfdb.tsv')
        vfdb = vfdb[['SEQUENCE','START','END','STRAND','GENE','%COVERAGE','%IDENTITY','PRODUCT']] 
        vfdb.rename(columns={'SEQUENCE':'Contig名称','START':'起始碱基','END':'终止碱基','STRAND':'正负链','GENE':'基因名称','%COVERAGE':'覆盖度%','%IDENTITY':'一致性%','PRODUCT':'产物'},inplace=True) # 添加并重命名列名
        vfdb['VFID'] = vfdb.apply(lambda x:getvfID(x),axis=1)   # 新增VFID列 = 对每一行调用函数 getvfID()，从产物列提取ID
        vfdb = vfdb.merge(vfmeta,on = 'VFID',how='left') 
        vfdb.rename(columns={'VF_Name':'VF名称','VF_FullName':'VF全称','Bacteria':'物种来源','VFcategory':'VF分类','Characteristics':'特征','Structure':'结构','Function':'功能','Mechanism':'机制','Reference':'文献来源'},inplace=True)
        vfdb.to_csv(f'{Pre}.vfdb.tsv',sep='\t',index=False)   # 合并VFs_meta.tsv数据库，重命名列名，输出为{Pre}.vfdb.tsv
        subprocess.run(f'''awk -v OFS='\t' '{{print $1,$2,$3,$5}}' {Pre}.vfdb.tsv|sed '1d' > vfdb.bed''',shell=True) # vfdb.bed = 选取{Pre}.vfdb.tsv数据库的1，2，3，5列，去除表头
        if int(os.popen('cat vfdb.bed|wc -l').read()) >0:
            vfdbdb = pd.read_table('vfdb.bed',header=None)
            vfdbdb[3] = vfdbdb.apply(lambda x:x[3] if len(x[3])<10 else f'{x[3][0:5]}..{x[3][-5:]}',axis=1)
            vfdbdb.to_csv('vfdb.bed',sep='\t',index=False,header=False)  # 基因名称（x[3]）的长度小于 10，保留原名称。≥ 10，处理为 ABCDE..MNOP
        else:
            subprocess.run(f'rm vfdb.bed',shell=True)
        
        # card耐药基因    
        subprocess.run(f'abricate --db card {Pre}.final.fasta --threads {threads} --minid 50 --mincov 50 >{Pre}.card.tsv',shell=True,stdout=f1,stderr=f1) # abricate进行card注释,最小相似性和最小覆盖度阈值均为 50%,结果保存到 {Pre}.card.tsv
        card = pd.read_table(f'{Pre}.card.tsv')
        card = card[['SEQUENCE','START','END','STRAND','GENE','%COVERAGE','%IDENTITY','PRODUCT','RESISTANCE']]  # 添加并重命名列名
        card.rename(columns={'SEQUENCE':'Contig名称','START':'起始碱基','END':'终止碱基','STRAND':'正负链','GENE':'基因名称','%COVERAGE':'覆盖度%','%IDENTITY':'一致性%','PRODUCT':'产物','RESISTANCE':'耐药药物'},inplace=True) 
        card.to_csv(f'{Pre}.card.tsv',sep='\t',index=False)  # 输出为{Pre}.card.tsv
        subprocess.run(f'''awk -v OFS='\t' '{{print $1,$2,$3,$5}}' {Pre}.card.tsv|sed '1d' > card.bed''',shell=True)    # card.bed = 选取{Pre}.card.tsv数据库的1，2，3，5列，去除表头
        if int(os.popen('cat card.bed|wc -l').read()) >0:
            carddb = pd.read_table('card.bed',header=None)
            carddb[3] = carddb.apply(lambda x:x[3] if len(x[3])<10 else f'{x[3][0:5]}..{x[3][-5:]}',axis=1)  # 基因名称（x[3]）的长度小于 10，保留原名称。≥ 10，处理为 ABCDE..MNOP
            carddb.to_csv('card.bed',sep='\t',index=False,header=False)
        else:
            subprocess.run('rm card.bed',shell=True)
        assem_vfdr(Pre,intype)   # 调用assem_vfdr函数
        tmpbind = '/'.join(os.getcwd().split('/')[:2])
        # 20260104
        if method != 'meta':
            subprocess.run(f'''/home/dell/miniconda3/bin/singularity exec --bind {tmpbind}:{tmpbind} /home/dell/biosoft/mummer2circos.simg mummer2circos -q {os.getcwd()}/{Pre}.final.fasta -r {os.getcwd()}/{Pre}.final.fasta -gb {os.getcwd()}/{Pre}_prokka/{Pre}.gbk -l  -o '{Pre}_raw' ''',shell=True,stdout=f1,stderr=f1)
            for feature in ['card','vfdb','rgi']:     # card/vfdb.bed，{Pre}.card/vfdb.tsv在鉴定时有门槛
                if os.path.isfile(f'{feature}.bed'):  # mummer2circos圈图：总共生成了四张圈图，标注了card基因的圈图在Rmarkdown中进行展示
                    tmpbind = '/'.join(os.getcwd().split('/')[:2])
                    subprocess.run(f'''/home/dell/miniconda3/bin/singularity exec --bind {tmpbind}:{tmpbind} /home/dell/biosoft/mummer2circos.simg mummer2circos -q {os.getcwd()}/{Pre}.final.fasta -r {os.getcwd()}/{Pre}.final.fasta -gb {os.getcwd()}/{Pre}_prokka/{Pre}.gbk -l -lf {feature}.bed -o '{Pre}_{feature}' ''',shell=True,stdout=f1,stderr=f1)
                    

def bp_vaccine(Pre):  # 百日咳疫苗基因型——BLASTn
    tmpdir = '/data1/shanghai_pip/meta_genome/database/BIGsdb/bordetella'
    tmpdict = {}
    #20250529 根据基因长度修改比对长度阈值 将prn鉴定错误,以及fim没有比对结果的报错去除
    genethrdict = {'23S_rRNA':2800,'fhaB-2400_5550':3100,'fim2':600,'fim3':600,'prn':2600,'ptxA':800,'ptxB':350,'ptxC':650,'ptxD':350,'ptxE':350,'ptxP':150,'tcfA':1900}
    for i in os.listdir(tmpdir):
        if i.endswith('.fas'):
            tgene = i.replace('.fas','') #  tgene = .fas前面部分，即保留纯基因名称
            bitlen = genethrdict.get(tgene,100)
            subprocess.run(f'''blastn -db {tmpdir}/{tgene} -query {Pre}.final.fasta -out {tgene}.blast.out -num_threads 10 -evalue 1e-5 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' -max_target_seqs 10  -perc_identity 90 -max_hsps 1''',shell=True)
            if os.path.getsize(f'{tgene}.blast.out') != 0:
                tmpdb = pd.read_table(f'{tgene}.blast.out',header=None)  # tmpdb = 运行BLASTn比对生成的结果文件{tgene}.blast.out
            #20250324 修正blast比对过短导致的鉴定错
                if tmpdb.shape[0]>0:
                    tmpdb = tmpdb.loc[tmpdb[3]>bitlen]
                    tmpdict[tgene] = {}
                    if tmpdb.shape[0] > 0:
                        tmpdict[tgene] = {'基因名称':tgene,'Contig名称':tmpdb.iloc[0,0],'起始位置':tmpdb.iloc[0,6],'终止位置':tmpdb.iloc[0,7],'分型':tmpdb.iloc[0,1],'一致性':tmpdb.iloc[0,2],'差异碱基数量':tmpdb.iloc[0,4]}
                    else:
                        tmpdict[tgene] = {'基因名称':tgene,'Contig名称':'-','起始位置':'-','终止位置':'-','分型':'-','一致性':'-','差异碱基数量':'-'}
            else:
                tmpdict[tgene] = {'基因名称':tgene,'Contig名称':'-','起始位置':'-','终止位置':'-','分型':'-','一致性':'-','差异碱基数量':'-'}
    newdb = pd.DataFrame(tmpdict).T.reset_index(drop=True)
    newdb.to_csv(f'{Pre}_scheme.tsv',sep='\t',index=False) # 保存每个基因的比对信息到{Pre}_scheme.tsv。对应Rmarkdown中的'特殊分型'

def serotype_HI(pre):  # 流感嗜血亚型预测——HICAP
    subprocess.run(f'''cut -d '' -f 1 {pre}.final.fasta > {pre}.fasta''',shell=True)  # 每一行的第一部分提取出来
    subprocess.run(f'/home/dell/miniconda3/bin/conda run -n GTDBtk python /home/dell/biosoft/hicap-master/hicap-runner.py -q  {pre}.final.fasta -o ./',shell=True) # 运行HICAP
    if os.path.isfile(f'{pre}.tsv'):
        hidb = pd.read_table(f'{pre}.tsv')
        hidb['样本名称'] = pre
        hidb.rename(columns={'predicted_serotype':'亚型预测','genes_identified':'检测基因','IS1016_hits':'IS1016数量'},inplace=True)
        hidb = hidb[['样本名称','亚型预测','检测基因','IS1016数量']]
        hidb.to_csv(f'{pre}_serotype_result.tsv',sep = '\t',index = False)
        sers = hidb['亚型预测'].tolist()[0]
    else:
        open(f'{pre}_serotype_result.tsv','w').write(f'样本名称\t亚型预测\t检测基因\tIS1016数量\n')
        open(f'{pre}_serotype_result.tsv','a').write(f'{pre}\t-\t-\t-')
        sers = '-'
    return sers   # 对应Rmarkdown中的'特殊分型'

#20260130。假结核和鼠疫耶尔森菌的区分
def serotype_ys(pre):   #耶尔森菌
    #1.PLA|PO2088|opgG|inv 2.97plus
    YS_dict = {'样本名':pre,'PLA':'-','YPO2088':'-','inv':'-','opgG':'-','97_predict':'Unknown','物种预测':'未知耶尔森'}
    subprocess.run(f'cat {pre}.final.fasta|seqkit amplicon -p /data/deploy/meta_genome/YS_primer.tsv --bed > YS_primer.bed',shell=True)
    if os.path.isfile(f'YS_primer.bed') and os.path.getsize('YS_primer.bed') != 0:
        primerdb = pd.read_table('YS_primer.bed',header=None)
        genelist = primerdb[3].tolist()
        for targene in genelist:
            YS_dict[targene] = '+'
    subprocess.run(f'python /data/deploy/meta_genome/Identify_Y.pestis-main/Identify_Y.pestis/Identify_Y.pestis_from_dir.py {pre}.final.fasta YS_97.txt',shell=True)
    if os.path.isfile(f'YS_97.txt') and os.path.getsize('YS_97.txt') != 0:
        YS97db = pd.read_table('YS_97.txt')
        YS_dict['97_predict'] = YS97db['Is_Ypestis'].tolist()[0]
    if YS_dict['PLA']=='+' and YS_dict['YPO2088'] == '+' and YS_dict['97_predict'] =='Yes':
        YS_dict['物种预测'] = '鼠疫耶尔森'
    elif YS_dict['opgG']=='+' and YS_dict['inv'] == '+' and YS_dict['97_predict'] =='No':
        YS_dict['物种预测'] = '假结核耶尔森'
    s1 = pd.DataFrame(YS_dict,index=[0])
    s1.to_csv(f'{pre}_serotype_result.tsv',sep = '\t',index =0)
    sers = s1['物种预测'].tolist()[0]

    return sers

def serotype_B(pre):   # 大肠+致贺——ectyper
    subprocess.run(f'/home/dell/miniconda3/bin/conda run -n cm210 ectyper -i ./{pre}.final.fasta -c {nt} -o {pre}_serotype --verify --pathotype',shell=True)
    s1 = pd.read_table(f'{pre}_serotype/output.tsv')
    serotype_s = '/data/deploy/bio-elite/bio/load_file/pathogenic/ecoli_seotype_anno.txt'
    s1 =s1[['Name', 'O-type', 'H-type', 'Serotype', 'Species','Pathotype','StxSubtypes']]
    s1.columns=['样本名','O抗原','H抗原','血清型','物种','分型','Stx亚型']
    s1['样本名']=pre
    s1['志贺分型'] = '-'
    if os.path.isfile('flex_out.tsv'):
        subprocess.run(f'rm flex_out.tsv',shell=True)
    if s1['物种'].tolist()[0] == 'Shigella flexneri':
        subprocess.run(f'cat {pre}.final.fasta |seqkit amplicon -p /data/test/cptyper_test/flex_primer.tsv --bed -o flex_out.tsv',shell=True)
        flexdb = pd.read_table('flex_out.tsv',header=None)
        flexplist = set(flexdb[3].tolist())
        flextype = '-'
        if flexplist:
            if 'wzx1' in flexplist:
                if flexplist == set(['wzx1','gtrI']):
                    flextype = '1a'
                elif flexplist == set(['wzx1','gtrI','oac']):
                    flextype = '1b'
                elif flexplist == set(['wzx1','gtrI','oac','gtrIC']):
                    flextype = '1c'
                elif flexplist == set(['wzx1','gtrII']):
                    flextype = '2a'
                elif flexplist == set(['wzx1','gtrII','gtrX']):
                    flextype = '2b' 
                elif flexplist == set(['wzx1','oac','gtrX']):
                    flextype = '3a'
                elif flexplist == set(['wzx1','oac']):
                    flextype = '3b'
                elif flexplist == set(['wzx1','gtrIV']):
                    flextype = '4a'
                elif flexplist == set(['wzx1','gtrIV','oac']):
                    flextype = '4b'
                elif flexplist == set(['wzx1','gtrV']):
                    flextype = '5a'
                elif flexplist == set(['wzx1','gtrX']):
                    flextype = 'X或Xv'
                elif flexplist == set(['wzx1']):
                    flextype = 'Y'
                else:
                    flextype = '-'
            else:
                if flexplist == set(['wzx6']):
                    flextype = 'F6'           
        s1['志贺分型'] = flextype
    s1.to_csv(f'{pre}_serotype_result.tsv',sep = '\t',index =0)
    sers = s1['血清型'].tolist()[0]

    return sers

def serotype_A(pre): # 沙门氏菌血清型——sistr
    #os.system(f"mkdir {run_path}/05.sero_type")
    #subprocess.run(f'/data/deploy/meta_genome/soft/SeqSero2-master/bin/SeqSero2_package.py -t 4 -m k -i {pre}.final.fasta -d sero_type',shell=True) # 运行seqsero2血清型鉴定
    subprocess.run(f'/home/dell/miniconda3/bin/conda run --no-capture-output -n GTDBtk sistr -i {pre}.final.fasta {pre} -f tab -o {pre}_result',shell=True) # 运行sistr进行血清型鉴定
    serotype_s = '/data/deploy/bio-elite/bio/load_file/pathogenic/salmonella_52seotype_v2.txt'
    sero_info = pd.read_table(serotype_s)
    sero = pd.read_table(f'{pre}_result.tab') # 读取 sistr 的结果文件 {pre}_result.tab
    sero=sero.rename(columns={'o_antigen':'O抗原','h1':'H1相抗原(fliC)','h2':'H2相抗原(fljB)'}) # 重命名
    sero['O抗原']=sero['O抗原'].astype('str')
    serox = sero.merge(sero_info,on=['O抗原','H1相抗原(fliC)'],how='left') # serox = salmonella_52seotype_v2.txt 与 sero合并
    serox.fillna('-',inplace=True)
    sero_result = pd.DataFrame()
    sero_result['样本']=pre
    sero_result[['O抗原','H1相抗原(fliC)','H2相抗原(fljB)']]=serox[['O抗原','H1相抗原(fliC)','H2相抗原(fljB)_y']]
    sero_result[['菌种','血清型','血清型注释信息(simple)','血清型注释信息(details)']] = serox[['serovar','serogroup','simple_description','details_description']]
    sero_result['抗原组成'] = sero_result['O抗原'] + ':'+ sero_result['H1相抗原(fliC)'] + ":" +sero_result['H2相抗原(fljB)']
    sero_result['样本']=pre
    sero_result['亚型全称'] = serox['name']
    sero_result.to_csv(f'{pre}_serotype_result.tsv',sep = '\t',index =0) # 结果保存为{pre}_serotype_result.tsv
    sers = sero_result['血清型'].tolist()[0]
    return sers

def serotype_kb(pre):  #克雷伯菌血清型分型——kleborate
    subprocess.run(f'kleborate --all -o results.txt -a {pre}.final.fasta > {pre}.keblo.tsv',shell=True) # 运行kleborate进行血清型分型，输出结果到{pre}.keblo.tsv
    kledb = pd.read_table(f'{pre}.keblo.tsv')
    kledb['样本名称'] = pre
    kledb.rename(columns={'virulence_score':'毒力得分','resistance_score':'耐药得分','Yersiniabactin':'耶尔森菌素','Colibactin':'大肠菌素','Bla_chr':'氨苄类耐药SHV等位基因','SHV_mutations':'SHV耐药突变','wzi':'wzi荚膜预测'},inplace=True)
    kledb['KO血清型'] = kledb['K_locus'].tolist()[0] + '|' + kledb['O_locus'].tolist()[0] # 添加新列'KO血清型' =  K_locus 和 O_locus 两列的第一个值用竖线'|'连接而成
    kledb = kledb[['样本名称','ST','毒力得分','耐药得分','耶尔森菌素','大肠菌素','氨苄类耐药SHV等位基因','SHV耐药突变','wzi荚膜预测','KO血清型']]
    kledb.to_csv(f'{pre}_serotype_result.tsv',sep = '\t',index = False) # 结果保存到 {pre}_serotype_result.tsv
    sers = kledb['KO血清型'].tolist()[0]
    return sers 
    
def serotype_nm(pre):  # 奈瑟氏菌亚型预测——PMGA 
    subprocess.run(f'/home/dell/miniconda3/bin/conda run --no-capture-output -n RGI pmga {pre}.final.fasta  -t 10 --force --blastdir /data1/shanghai_pip/meta_genome/database/pmga/',shell=True) # 运行PMGA 
    nmdb = pd.read_table(f'pmga/{pre}.finalsta.txt')
    nmdb['样本名称'] = pre
    nmdb.rename(columns={'prediction':'亚型预测','genes_present':'验证基因集','notes':'注释'},inplace=True)
    nmdb = nmdb[['样本名称','亚型预测','验证基因集','注释']]
    nmdb.to_csv(f'{pre}_serotype_result.tsv',sep = '\t',index = False)
    sers = nmdb['亚型预测'].tolist()[0]
    return sers
    
def serotype_D(pre):  # 副溶血弧菌血清型——VPsero
    serotype_s='/data/deploy/bio-elite/bio/load_file/pathogenic/vp_serotype_anno.txt'
    subprocess.run(f'mkdir VPsero;cp ./{pre}.final.fasta VPsero',shell=True)  # 将{pre}.final.fasta复制到VPsero
    subprocess.run(f'python /home/dell/biosoft/VPsero-master/program.py -i VPsero -o my_out_put_2  -n {nt}',shell=True)  # 调用 VPsero 工具
    s1 =pd.read_excel('my_out_put_2/serotype_predict/04.predict_result/all_strain_predict_result.xlsx')
    s1 =s1[['O_Spec_Gene', 'K_Spec_Gene','Predict_O_sero', 'Predict_K_sero', 'New_serotype']]
    s1.columns=['O血清型基因','K血清型基因','O血清型','K血清型','血清型类型']
    s1['血清型']= s1['O血清型']+':'+ s1['K血清型']
    s1.to_csv(f'{pre}_serotype_result.tsv',sep = '\t',index =0) # 输出结果为{pre}_serotype_result.tsv
    subprocess.run('cp my_out_put_i/serotype_predict/04.predict_result/all_strain_predict_result.xlsx {pre}_strain_predict_result.xlsx')
    sers = s1['血清型'].tolist()[0]
    return sers

def serotype_groupA(pre):     #GAS，A链emm分型——emm_typing.py
    subprocess.run(f'python /home/dell/biosoft/emm_typing/emm_typing/emm_typing.py -f {pre}.final.fasta',shell=True) # 调用emm_typing.py 脚本
    gasdb = pd.read_table(f'emm_results.tab')
    gasdb.rename(columns={'emm-type':'亚型预测','pident':'一致性','Isolate':'样本名称','length':'比对长度'},inplace=True)
    gasdb = gasdb[['样本名称','亚型预测','一致性','比对长度']]
    gasdb.to_csv(f'{pre}_serotype_result.tsv',sep = '\t',index = False)  # 结果输出为{pre}_serotype_result.tsv
    sers = gasdb['亚型预测'].tolist()[0]
    return sers

def serotype_E(pre):     # 霍乱亚型预测——自行编写的PCR脚本
    serodict = {}
    #subprocess.run(f'/home/dell/miniconda3/bin/conda run -n GTDBtk python /home/dell/biosoft/vista/vista.py search --data-path /home/dell/biosoft/vista/data/ --cpus 10 {pre}.final.fasta > {pre}.serotype.json',shell=True)
    subprocess.run(f'''/home/dell/miniconda3/bin/conda run -n  choleraefinder python /data/test/test_VP/choleraefinder/choleraefinder.py --input {pre}.final.fasta -o ./ -p /data/test/test_VP/choleraefinder_db/  -t 0.95 -l 0.95 -q''',shell=True)
    serojson = json.load(open(f'data_CholeraeFinder.json'))
    serodict['样本名称'] = pre
    serodict['血清型'] = serojson['choleraefinder']['typing_cholerae']['serogroup']
    serodict['生物型'] = serojson['choleraefinder']['typing_cholerae']['biotype']
    serodb = pd.DataFrame(serodict,index=[0])
    serodb.to_csv(f'{pre}_serotype_result.tsv',sep='\t',index=False)   # 结果保存到{pre}_serotype_result.tsv
    sers = serojson['choleraefinder']['typing_cholerae']['serogroup']
    return sers

def serotype_MLVA(pre):    #布鲁氏菌MLVA分析——MLVA_finder.py
    if not os.path.isdir(f'{pre}_mlvafafile'):
        os.makedirs(f'{pre}_mlvafafile')
    subprocess.run(f'cp {pre}.final.fasta {pre}_mlvafafile',shell=True)     # 调用MLVA_finder.py，引物文件 Brucella_primers.txt
    subprocess.run(f'python /data1/shanghai_pip/meta_genome/database/MLVA_finder/MLVA_finder.py -i {pre}_mlvafafile -o ./ -p /data1/shanghai_pip/meta_genome/database/MLVA_finder/data_test/primers/Brucella_primers.txt',shell=True)
    mlvadb = pd.read_table(f'{pre}_mlvafafile_output.csv',sep=',')
    mlvadb['样本名称'] = pre
    mlvadb.rename(columns={'primer':'引物','position1':'起始位置','position2':'终止位置','size':'扩增片段大小','allele':'重复基因数量'},inplace=True)
    mlvadb = mlvadb[['样本名称','引物','起始位置','终止位置','扩增片段大小','重复基因数量']]
    mlvadb.to_csv(f'{pre}_serotype_result.tsv',sep = '\t',index = False)
    sers = ';'.join(mlvadb['allele'].tolist())
    return sers

def serotype_st(pre):  # 金葡家系分析——SALTY
    if not os.path.isdir(f'{pre}_st_fafile'):
        os.makedirs(f'{pre}_st_fafile')
    subprocess.run(f'cp {pre}.final.fasta {pre}_st_fafile',shell=True)
    subprocess.run(f'/home/dell/miniconda3/bin/conda run --no-capture-output -n RGI salty -i {pre}_st_fafile -o {pre}_st_fafile -t 10',shell=True)
    stdb = pd.read_table(f'{pre}_st_fafile/summaryReport.txt')
    stdb['样本名称'] = pre
    stdb.rename(columns={'Lineage':'家系','SACOL1908':'SACOL1908基因座等位基因','SACOL0451':'SACOL0451基因座等位基因','SACOL2725':'SACOL2725基因座等位基因'},inplace=True)  
    stdb = stdb[['样本名称','家系','SACOL1908基因座等位基因','SACOL0451基因座等位基因','SACOL2725基因座等位基因']]   
    stdb.to_csv(f'{pre}_serotype_result.tsv',sep = '\t',index = False)   # 输出到{pre}_serotype_result.tsv
    sers = stdb['家系'].tolist()[0]
    return sers

def serotype_lm(pre):    # 单增李斯特血清型——lissero
    subprocess.run(f'/home/dell/miniconda3/bin/conda run --no-capture-output -n RGI lissero {pre}.final.fasta > Lm_sero.tsv ',shell=True)
    stdb = pd.read_table(f'Lm_sero.tsv')
    stdb['样本名称'] = pre
    stdb.rename(columns={'SEROTYPE':'血清型'},inplace=True)
    stdb = stdb[['样本名称','血清型','PRS','LMO0737','LMO1118','ORF2110','ORF2819']]  
    stdb.to_csv(f'{pre}_serotype_result.tsv',sep = '\t',index = False)  # 结果保存为{pre}_serotype_result.tsv
    sers = stdb['血清型'].tolist()[0]
    return sers

def serotype_bt(pre):   # 蜡样芽孢杆菌毒力因子——btyper3
    subprocess.run(f'btyper3 -i {pre}.final.fasta -o Bcere_sero',shell=True)  # 调用btyper3
    btdb = pd.read_table(f'Bcere_sero/btyper3_final_results/{pre}.final_final_results.txt')
    btdb['样本名称'] = pre
    btdb.rename(columns={'final_taxon_names':'物种名称','anthrax_toxin(genes)':'nthrax基因集','emetic_toxin_cereulide(genes)':'cereulide基因集','diarrheal_toxin_Nhe(genes)':'Nhe基因集','diarrheal_toxin_Hbl(genes)':'Hbl基因集','diarrheal_toxin_CytK(top_hit)':'CytK基因集','sphingomyelinase_Sph(gene)':'Sph基因集','capsule_Cap(genes)':'Cap基因集'},inplace=True)
    btdb = btdb[['样本名称','物种名称','nthrax基因集','cereulide基因集','Nhe基因集','Hbl基因集','CytK基因集','Sph基因集','Cap基因集']]
    btdb.to_csv(f'{pre}_serotype_result.tsv',sep = '\t',index = False)
    sers = btdb['物种名称'].tolist()[0]
    return sers

def serotype_SS(Pre):   # 猪链血清型
    #1.构建cps库 2.获取突变位点
    database='/data/deploy/meta_genome/database/SS_sero/'
    blastref=f'{database}/Ssuis_Serotyping'
    cpsref=f'{database}/Ssuis_cps2K.fasta'
    with open('sero.log','w') as serof:
        subprocess.run(f'blastn -query {Pre}.final.fasta -db {blastref} -out {Pre}.cps.out.tsv -outfmt 6 -perc_identity 90 ',shell=True)
        Serodb = pd.read_table(f'{Pre}.cps.out.tsv',header=None)
        lengdb = pd.read_table(f'{database}/Ssuis_Serotyping.tsv')
        Serodb = Serodb.merge(lengdb,left_on=1,right_on='Serotype')
        Serodb['Perc'] = Serodb[3]/Serodb['Length']
        Serodb = Serodb.loc[Serodb['Perc']>0.9]
        if Serodb.shape[0] > 0:
            Serotype = Serodb['Serotype'].tolist()[0]
        else:
            Serotype = 'notype'
        #print(Serotype)
        Serotype = Serotype.replace('cps-','')
        if Serotype == '1' or Serotype == '2':
            subprocess.run(f'nucmer --maxmatch -b 200 -c 65 -d 0.12 -g 90 -l 20 {cpsref} {Pre}.final.fasta -p {Pre}',shell=True,stdout=serof,stderr=serof)
            subprocess.run(f'show-snps {Pre}.delta -T > {Pre}.snps.out',shell=True)
            snpdb = pd.read_table(f'{Pre}.snps.out',skiprows=2)
            snplist = snpdb['[P2]'].tolist()
            if Serotype == '1':
                if '483' in snplist:
                    pass
                else:
                    Serotype = '14'
            else:
                if '483' in snplist:
                    Serotype = '1/2'
                else:
                    pass

        #otherResult
        if os.path.isfile(f'{Pre}.SsuisChara.tsv'):
            subprocess.run(f'rm -r {Pre}.SsuisChara.tsv',shell=True)
        subprocess.run(f'python /data/deploy/meta_genome/database/SsuisChara/SsuisChara.py -i {Pre}.final.fasta -o {Pre}.SsuisChara.tsv',shell=True)
        SsuisCdb = pd.read_table(f'{Pre}.SsuisChara.tsv')
        SsuisCdb = SsuisCdb[['human infection potential','AMRG_level','aminoglycoside','macrolide','tetracycline']]
        SsuisCdb['样本名称'] = Pre
        SsuisCdb['血清型'] = Serotype
        SsuisCdb = SsuisCdb.rename(columns={'human infection potential':'感染等级','aminoglycoside':'氨基糖苷类','macrolide':'大环内酯类','tetracycline':'四环素类','AMRG_level':'耐药数量'})
        SsuisCdb = SsuisCdb[['样本名称','血清型','感染等级','耐药数量','氨基糖苷类','大环内酯类','四环素类']]
        SsuisCdb.to_csv(f'{Pre}_serotype_result.tsv',sep='\t',index=False)
        return Serotype

#20250317 add mlva results for BP.目前脚本VNTR3b始终为0
def len2mlvacopy(PrimerN,PrimerL):
    primerdict = {'VNTR1':[383,15,9],'VNTR3a':[135,5,7],'VNTR3b':[0,0,0],'VNTR4':[232,12,9],'VNTR5':[143,6,7],'VNTR6':[234,9,11]}
    modlen,modcopylen,modcopynum = primerdict[PrimerN]
    if PrimerN != 'VNTR3b':
        PrimerCopy = round(modcopynum-(modlen-PrimerL)/modcopylen)
    else:
        PrimerCopy = 0
    return PrimerCopy
def bp_mlva(Pre,primer='/data/test/mlva/mlva_primer.tsv',typetable='/data1/shanghai_pip/meta_genome/BPMLVA.table.tsv'):
    subprocess.run(f'cat {Pre}.final.fasta |seqkit amplicon -p {primer} --bed > {Pre}.raw.mlva.tsv',shell=True)
    rawdb = pd.read_table(f'{Pre}.raw.mlva.tsv',header=None,names=['Chrom','Startpos','Endpos','PrimerName','Mismatch','Strand','Sequence'])
    finaldict = {}
    for PrimerN in rawdb['PrimerName'].tolist():
        tmpdb  = rawdb.loc[rawdb['PrimerName']==PrimerN,]
        if PrimerN == 'VNTR3':
            if tmpdb.shape[0] == 1:
                PrimerN = 'VNTR3a'
                finaldict[PrimerN] = len2mlvacopy(PrimerN,tmpdb['Endpos'].tolist()[0]-tmpdb['Startpos'].tolist()[0])
                finaldict['VNTR3b'] = len2mlvacopy('VNTR3b',tmpdb['Endpos'].tolist()[0]-tmpdb['Startpos'].tolist()[0])
            else:
                #20250520判断N个亚型是否一致,多个VNTR3a时都计算拷贝数，相同写入一个值，否则‘；’拼接多个值。VNTR3b始终为0
                VNTR3list = []
                tmpdb = tmpdb.reset_index()
                #finaldict[PrimerN] = len2mlvacopy(PrimerN,tmpdb['Endpos'].tolist()[0]-tmpdb['Startpos'].tolist()[0])
                for tmpi in tmpdb.index:
                    Endpos = tmpdb.loc[tmpdb.index==tmpi,'Endpos'].tolist()[0]
                    Startpos = tmpdb.loc[tmpdb.index==tmpi,'Startpos'].tolist()[0]
                    tmptype = len2mlvacopy('VNTR3a',tmpdb['Endpos'].tolist()[0]-tmpdb['Startpos'].tolist()[0])
                    VNTR3list.append(tmptype)
                if len(set(VNTR3list)) == 1:
                    finaldict['VNTR3a'] = VNTR3list[0]

                else:
                    finaldict['VNTR3a'] = ';'.join(list(set(VNTR3list)))
                finaldict['VNTR3b'] = 0
        elif PrimerN == 'VNTR4': # 出现多次VNTR4时，标记为"未知"
            if tmpdb.shape[0] ==1:
                finaldict['VNTR4'] = len2mlvacopy(PrimerN,tmpdb['Endpos'].tolist()[0]-tmpdb['Startpos'].tolist()[0])
            else:
                finaldict['VNTR4'] = "未知"
        else:
            finaldict[PrimerN] = len2mlvacopy(PrimerN,tmpdb['Endpos'].tolist()[0]-tmpdb['Startpos'].tolist()[0])
    typedb = pd.read_table(typetable)
    typedb = typedb.astype('str')
    typedb['combine'] = typedb.apply(lambda x: x['VNTR1']+'_'+x['VNTR3a']+'_'+x['VNTR3b']+'_'+x['VNTR4']+'_'+x['VNTR5']+'_'+x['VNTR6'],axis=1) #9_7_0_9_7_11
    typedb.to_csv('tt1.tsv',sep='\t',index=False)
    finaldb = pd.DataFrame(finaldict,index=[0]).reset_index(drop=True)
    finaldb = finaldb.astype('str')
    if finaldb.shape[1] == 6:
        finaldb['combine'] =finaldb.apply(lambda x:x['VNTR1']+'_'+x['VNTR3a']+'_'+x['VNTR3b']+'_'+x['VNTR4']+'_'+x['VNTR5']+'_'+x['VNTR6'],axis=1)
        finaldb = finaldb.merge(typedb,on='combine',how='left')
        finaldb['样本名称'] = Pre
        finaldb.rename(columns = {'VNTR1_x':'VNTR1','VNTR3a_x':'VNTR3a','VNTR3b_x':'VNTR3b','VNTR4_x':'VNTR4','VNTR5_x':'VNTR5','VNTR6_x':'VNTR6'},inplace=True)
        finaldb = finaldb[['样本名称','MT','VNTR1','VNTR3a','VNTR3b','VNTR4','VNTR5','VNTR6']]
        finaldb.to_csv(f'{Pre}.mlva.tsv',sep='\t',index=False)
    else:
        print(f'{Pre} {finaldb}')

#20250714
def len2mlvacopy_mp(PrimerN,PrimerL):
    primerdict = {'Mpn13':[428,16,4],'Mpn14':[378,21,4],'Mpn15':[192,21,4],'Mpn16':[447,47,4]}
    modlen,modcopylen,modcopynum = primerdict[PrimerN]
    PrimerCopy = math.ceil(modcopynum-(modlen-PrimerL)/modcopylen)
    return PrimerCopy
def mp_mlva(Pre,primer='/data/test/mlva/mlva_primer_mp.tsv'):
    subprocess.run(f'cat {Pre}.final.fasta |seqkit amplicon  -p {primer} --bed > {Pre}.raw.mlva.tsv',shell=True)
    rawdb = pd.read_table(f'{Pre}.raw.mlva.tsv',header=None,names=['Chrom','Startpos','Endpos','PrimerName','Mismatch','Strand','Sequence'])
    finaldict = {}
    for PrimerN in rawdb['PrimerName'].tolist():
        tmpdb  = rawdb.loc[rawdb['PrimerName']==PrimerN,]
        finaldict[PrimerN] = len2mlvacopy_mp(PrimerN,tmpdb['Endpos'].tolist()[0]-tmpdb['Startpos'].tolist()[0])
    finaldb = pd.DataFrame(finaldict,index=[0]).reset_index(drop=True)
    finaldb = finaldb.astype('str')
    print(finaldb)
    if finaldb.shape[1] == 4:
        finaldb['combine'] =finaldb.apply(lambda x:x['Mpn13']+'_'+x['Mpn14']+'_'+x['Mpn15']+'_'+x['Mpn16'],axis=1)
    #    finaldb = finaldb.merge(typedb,on='combine',how='left')
        finaldb['样本名称'] = Pre
    #    finaldb.rename(columns = {'VNTR1_x':'VNTR1','VNTR3a_x':'VNTR3a','VNTR3b_x':'VNTR3b','VNTR4_x':'VNTR4','VNTR5_x':'VNTR5','VNTR6_x':'VNTR6'},inplace=True)
        finaldb = finaldb[['样本名称','combine','Mpn13','Mpn14','Mpn15','Mpn16']]
        finaldb.to_csv(f'{Pre}.mlva.tsv',sep='\t',index=False)
    else:
        print(f'{Pre} {finaldb}')
#20250317 A2037G AMR for BP
#20250429 解决连续碱基识别错误
def bp_2037(Pre):
    #1.mapping rrn 2.call snp 3.collect data
    with open('bp2037.log','w') as f: 
        if os.path.isfile(f'{Pre}.R2.fastq.gz'):
            subprocess.run(f'bwa mem /data1/shanghai_pip/meta_genome/rrn.fasta {Pre}.R1.fastq.gz {Pre}.R2.fastq.gz -t 10 |samtools sort -o {Pre}.rrn.sorted.bam',shell=True,stdout=f,stderr=f)
        else:
            subprocess.run(f'bwa mem /data1/shanghai_pip/meta_genome/rrn.fasta {Pre}.R1.fastq.gz -t 10 |samtools sort -o {Pre}.rrn.sorted.bam',shell=True,stdout=f,stderr=f)
        subprocess.run(f'samtools index {Pre}.rrn.sorted.bam',stdout=f,stderr=f,shell=True)
        subprocess.run(f'freebayes -f /data1/shanghai_pip/meta_genome/rrn.fasta {Pre}.rrn.sorted.bam > {Pre}.rrn.vcf',stdout=f,stderr=f,shell=True)
        subprocess.run(f'/home/dell/miniconda3/envs/PathoSource/bin/vt normalize {Pre}.rrn.vcf -r /data1/shanghai_pip/meta_genome/rrn.fasta -o {Pre}.rrn.filt1.vcf',shell=True,stdout=f,stderr=f) # 最简+左对齐原则
        skip_rows = int(os.popen(f'''grep '##' {Pre}.rrn.filt1.vcf|wc -l ''').read())
        rrndb = pd.read_table(f'{Pre}.rrn.filt1.vcf',skiprows=skip_rows)
        if rrndb.shape[0] > 0:
            rrndb['GT'] = rrndb['unknown'].str.split(':').str[0]
            rrndb = rrndb[rrndb['GT']!='0/0']         # 0/0代表软件预测突变不可信
            if rrndb.shape[0]>0:
                rrndb['POS'] = rrndb['POS'].astype('str')
                rrndb['symbol'] = rrndb['REF']+rrndb['POS']+rrndb['ALT']
                if any(rrndb['symbol']=='A2037G'):
                    rrndb.to_csv(f'{Pre}.2037.tsv',sep='\t',index=False)

def fa_2037(Pre):
    with open('2037.log','w') as f:
        subprocess.run(f'nucmer /data1/shanghai_pip/meta_genome/rrn.fasta {Pre}.final.fasta',shell=True,stdout=f,stderr=f)
        if int(os.popen(f'show-snps out.delta|grep -w 2037|wc -l').read()) >= 1:
            subprocess.run(f'show-snps out.delta > {Pre}.2037.tsv',shell=True)

#物种名称：传入的taxid根据tmp.txt转化为相应的mlst名称，如果一个taxid命中多行用kraken2再鉴定一次。两种弧菌的kraken2名称，需要更正
def level2Spe(Pre,taxid):
    #----get taxid vs narmal pathogen
    id2dict = {'Vibrio parahaemolyticus':'vparahaemolyticus','Vibrio cholerae':'vcholerae'}
    taxid = str(taxid)
    Spedb = pd.read_table('/data/test/level2Spe/tmp.txt')
    Spedb['idlist'] = Spedb['idlist'].str.split(';')
    nSpedb = Spedb[Spedb['idlist'].apply(lambda x:taxid in x)]
    if nSpedb.shape[0] == 1:
        nSpe = nSpedb['Spe'].tolist()[0]
    else:
        if nSpedb.shape[0] > 1:
            #nSpe = ';'.join(nSpedb['Spe'].tolist())
            if not os.path.isfile(f'{Pre}.kraken2.report.txt'):
                subprocess.run(f'kraken2 --db /home/dell/kraken2_custom_202101_24G {Pre}.final.fasta --report {Pre}.kraken2.report.txt --output {Pre}.kraken2.txt -t 10',shell=True)
            fak2db = pd.read_table(f'{Pre}.kraken2.report.txt',header=None)
            nSpe = id2dict.get(fak2db.loc[fak2db[3]=='S'][5].str.strip().tolist()[0],0)
        else:
            nSpe = 0
    return nSpe


#20251111，弯曲菌血清型及鉴定
Prim = Tuple[str, str]   # (mix, primer_name), e.g. ("beta","Mu_HS31")
# canonical mixes
A = "alpha"; B = "beta"; G = "gamma"; D = "delta"
# ---- Primer → Mix mapping (from Table 3) ----
PRIMER_TO_MIX: Dict[str, str] = {
    # Alpha
    "Mu_HS2": A, "Mu_HS3": A, "Mu_HS4A": A, "Mu_HS6": A, "Mu_HS10": A, "Mu_HS15": A,
    "Mu_HS41": A, "Mu_HS53": A, "Mu_HS19": A, "Mu_HS63": A, "Mu_HS33": A,
    # Beta
    "Mu_HS1": B, "Mu_HS4B": B, "Mu_HS8": B, "Mu_HS23/36": B, "Mu_HS42": B, "Mu_HS57": B,
    "Mu_HS12": B, "Mu_HS27": B, "Mu_HS21": B, "Mu_HS31": B,
    # Gamma
    "Mu_HS44": G, "Mu_HS45": G, "Mu_HS29": G, "Mu_HS22": G, "Mu_HS9": G, "Mu_HS37": G,
    "Mu_HS18": G, "lpxA": G,
    # Delta
    "Mu_HS58": D, "Mu_HS52": D, "Mu_HS60": D, "Mu_HS55": D, "Mu_HS32": D, "Mu_HS11": D,
    "Mu_HS40": D, "Mu_HS38": D,
}

# ---- Aliases (resolve Mu_HS5 vs Mu_HS31) ----
ALIASES: Dict[str, str] = {
    "Mu_HS5": "Mu_HS31",
}
def norm_primer(p: str) -> str:
    return ALIASES.get(p, p)

class Rule:
    def __init__(self,
                 name: str,
                 must_have: Set[Prim],
                 must_not: Set[Prim] = None,
                 allowed_extra: Set[Prim] = None):
        self.name = name
        self.must_have = {(m, norm_primer(p)) for m, p in (must_have or set())}
        self.must_not = {(m, norm_primer(p)) for m, p in (must_not or set())}
        self.allowed_extra = {(m, norm_primer(p)) for m, p in (allowed_extra or set())}

# Rules based on Table 4
A = "alpha"; B = "beta"; G = "gamma"; D = "delta"
RULES: List[Rule] = [
    Rule("HS1",  {(B,"Mu_HS1")}),
    Rule("HS2",  {(A,"Mu_HS2")}),
    Rule("HS3",  {(A,"Mu_HS3")}),
    Rule("HS4 complex", {(A,"Mu_HS4A")}),
    Rule("CG8486 (HS4 complex member)", {(B,"Mu_HS4B")}),
    Rule("HS5",  {(B,"Mu_HS31"), (G,"Mu_HS45")}),
    Rule("HS6",  {(A,"Mu_HS6")}),
    Rule("HS7",  {(A,"Mu_HS6")}),
    Rule("HS8",  {(B,"Mu_HS8")}),
    Rule("HS9",  {(G,"Mu_HS9")}),
    Rule("HS10", {(A,"Mu_HS10")}),
    Rule("HS11", {(D,"Mu_HS11")}),
    Rule("HS12", {(B,"Mu_HS12")}),
    Rule("HS13", {(A,"Mu_HS4A")}),
    Rule("HS15", {(A,"Mu_HS15")}, allowed_extra={(D,"Mu_HS58")}),
    Rule("HS16", {(A,"Mu_HS4A"), (B,"Mu_HS4B")}, allowed_extra={(D,"Mu_HS52")}),
    Rule("HS17", {(B,"Mu_HS8")}),
    Rule("HS18", {(G,"Mu_HS18")}),
    Rule("HS19", {(A,"Mu_HS19")}),
    Rule("HS21", {(B,"Mu_HS21")}),
    Rule("HS22", {(G,"Mu_HS22")}),
    Rule("HS23", {(B,"Mu_HS23/36")}),
    Rule("HS27", {(B,"Mu_HS27")}),
    Rule("HS29", {(G,"Mu_HS29")}),
    Rule("HS31", {(B,"Mu_HS31")}, allowed_extra={(A,"Mu_HS15")}),
    Rule("HS32", {(D,"Mu_HS32"), (G,"Mu_HS45")}, allowed_extra={(B,"Mu_HS8")}),
    Rule("HS33", {(A,"Mu_HS33")}),
    Rule("HS35", {(A,"Mu_HS33")}),
    Rule("HS36", {(B,"Mu_HS23/36")}),
    Rule("HS37", {(G,"Mu_HS37")}),
    Rule("HS38", {(D,"Mu_HS38")}),
    Rule("HS40", {(D,"Mu_HS40")}),
    Rule("HS41", {(A,"Mu_HS41")}),
    Rule("HS42", {(B,"Mu_HS42")}),
    Rule("HS43", {(A,"Mu_HS4A")}),
    Rule("HS44", {(G,"Mu_HS44")}),
    Rule("HS45", {(G,"Mu_HS45")}, must_not={(B,"Mu_HS31"), (D,"Mu_HS32"), (D,"Mu_HS60")}),
    Rule("HS50", {(A,"Mu_HS4A")}),
    Rule("HS52", {(D,"Mu_HS52")}),
    Rule("HS53", {(A,"Mu_HS53")}),
    Rule("HS55", {(D,"Mu_HS55")}),
    Rule("HS57", {(B,"Mu_HS57")}),
    Rule("HS58", {(D,"Mu_HS58")}, allowed_extra={(A,"Mu_HS15")}),
    Rule("HS60", {(D,"Mu_HS60"), (G,"Mu_HS45")}),
    Rule("HS62", {(A,"Mu_HS4A")}),
    Rule("HS63", {(A,"Mu_HS63")}),
    Rule("HS64", {(A,"Mu_HS4A"), (B,"Mu_HS4B")}),
    Rule("HS65", {(A,"Mu_HS4A")}),
]

LPXA = (G, "lpxA")
def primers_to_pairs(primer_names: Iterable[str]) -> Set[Prim]:
    pairs: Set[Prim] = set()
    for raw in primer_names:
        p = norm_primer(str(raw).strip())
        mix = PRIMER_TO_MIX.get(p)
        if mix:
            pairs.add((mix, p))
    return pairs
def call_capsule_types(detected: Iterable[Prim], require_lpxA: bool = False) -> List[str]:
    det: Set[Prim] = {(m, norm_primer(p)) for m, p in detected}

    if require_lpxA and LPXA not in det:
        return ["Uninterpretable (lpxA negative)"]

    calls: List[str] = []
    for r in RULES:
        if not r.must_have.issubset(det):
            continue
        if any((m,p) in det for (m,p) in r.must_not):
            continue
        calls.append(r.name)

    if "HS4 complex" in calls:
        #if any(x in calls for x in ["HS13","HS43","HS50","HS62","HS65","HS16","HS64"]):
        #    calls = [c for c in calls if c != "HS4 complex"]
        calls = ['HS4']

    return sorted(set(calls)) or ["Untypeable"]
def call_capsule_types_from_primernames(primer_names: Iterable[str], require_lpxA: bool = False) -> List[str]:
    pairs = primers_to_pairs(primer_names)
    return call_capsule_types(pairs, require_lpxA=require_lpxA)
def serotype_Cb(Pre):   #空肠弯曲菌
    sers = '-'
    subprocess.run(f' cat {Pre}.final.fasta |seqkit amplicon -p /data/test/CB_type/CB_primer.tsv --bed > {Pre}.CB_primer.tsv',shell=True)
    serodict = {'样本名称':Pre,'血清型':'-','物种可靠':'否'}
    if os.path.isfile(f'{Pre}.CB_primer.tsv') and os.path.getsize(f'{Pre}.CB_primer.tsv') != 0:
        CBdb = pd.read_table(f'{Pre}.CB_primer.tsv',header=None)
        if CBdb.shape[0] > 0:
            #1.判断物种 2.判断血清型
            if 'lpxA' in CBdb[3].tolist():
                serodict['物种可靠'] = '是'
            sers = call_capsule_types_from_primernames(CBdb[3].tolist())[0]
            sers = sers.replace('HS','HS:')
            serodict['血清型'] = sers
    serodb = pd.DataFrame(serodict,index=[0])
    serodb.to_csv(f'{Pre}_serotype_result.tsv',sep='\t',index=False)
    return sers


#20251112，特殊毒力基因抽取表格
def PathoNet(Pre,species):
    PathoSamdict = {'样本名称':Pre,'物种':species,'血清型':'-','毒力基因':'-'}
    PathoNetdict = {'vcholerae':{'serotype':['O1','O139'],
    'vfgene':['ctxA','ctxB']},
    'senterica':{'serotype':['S.Typhi','S.Paratyphi A','S.Paratyphi B','S.Paratyphi C','S.Enteritidis','S.Typhimurium','S.Choleracsuis','S.Derby','S.London','S.Stanley','S.Calabar','S.Agona','S.Thompson','S.Rissen','S.enterica subsp. enterica serovar Typhimurium monophasic variant'],'vfgene':[]},
    'campylobacter':{'serotype':['HS:1','HS:2','HS:4','HS:19','HS:23','HS:41','HS:44'],'vfgene':['hcp','virB','ciaB','ggt','cdtA','cdtB','ctdC','cgtA','cgtB','wlaN','cstII']},
    'klebsiella':{'serotype':['K1','K2','K5','K20','K54','K57'],'vfgene':[]},
    'ecoli':{'serotype':['O2','O45','O103','O111','O121','O145','O157'],'vfgene':['stx1A','stx1B','stx2A','stx2B','stxA']},
    'Shigella':{'serotype':['1a','1b','1c','2a','2b','3a','3b','4a','4b','5a','5b','X','Xv','F6','Y'],'vfgene':['stx1A','stx1B','stx2A','stx2B','stxA']}}
    if species == 'campylobacter':
        pathodict = PathoNetdict[species]
        if os.path.isfile(f'{Pre}.vfdb.tsv') and os.path.getsize(f'{Pre}.vfdb.tsv') != 0:
            cpvfdb = pd.read_table(f'{Pre}.vfdb.tsv')
            tarvflist = [i for i in pathodict['vfgene'] if i in  cpvfdb['基因名称'].tolist()]
            if tarvflist:
                PathoSamdict['毒力基因'] = ';'.join(tarvflist)
        if os.path.isfile(f'{Pre}_serotype_result.tsv') and os.path.getsize(f'{Pre}_serotype_result.tsv') != 0:
            cpserodb = pd.read_table(f'{Pre}_serotype_result.tsv')
            if cpserodb['血清型'].tolist()[0] in pathodict['serotype']:
                PathoSamdict['血清型'] = f'''{cpserodb['血清型'].tolist()[0]}(重点关注)'''
            else:
                PathoSamdict['血清型'] = cpserodb['血清型'].tolist()[0]
    if species == 'klebsiella':
        pathodict = PathoNetdict[species]
        if os.path.isfile(f'{Pre}.vfdb.tsv') and os.path.getsize(f'{Pre}.vfdb.tsv') != 0:
            klvfdb = pd.read_table(f'{Pre}.vfdb.tsv')
            tarvflist = [i for i in pathodict['vfgene'] if i in  klvfdb['基因名称'].tolist()]
            if tarvflist:
                PathoSamdict['毒力基因'] = ';'.join(tarvflist)
        if os.path.isfile(f'{Pre}_serotype_result.tsv') and os.path.getsize(f'{Pre}_serotype_result.tsv') != 0:
            cpserodb = pd.read_table(f'{Pre}_serotype_result.tsv')
            if cpserodb['KO血清型'].tolist()[0].split('|')[0].replace('KL','K') in pathodict['serotype']:
                PathoSamdict['血清型'] = f'''{cpserodb['KO血清型'].tolist()[0]}.split('|')[0].replace('KL','K')(重点关注)'''
            else:
                PathoSamdict['血清型'] = cpserodb['KO血清型'].tolist()[0].split('|')[0].replace('KL','K')

    if species == 'senterica':
        pathodict = PathoNetdict[species]
        if os.path.isfile(f'{Pre}.vfdb.tsv') and os.path.getsize(f'{Pre}.vfdb.tsv') != 0:
            salvfdb = pd.read_table(f'{Pre}.vfdb.tsv')
            tarvflist = [i for i in pathodict['vfgene'] if i in  salvfdb['基因名称'].tolist()]
            if tarvflist:
                PathoSamdict['毒力基因'] = ';'.join(tarvflist)
        if os.path.isfile(f'{Pre}_serotype_result.tsv') and os.path.getsize(f'{Pre}_serotype_result.tsv') != 0:
            cpserodb = pd.read_table(f'{Pre}_serotype_result.tsv')
            if cpserodb['亚型全称'].tolist()[0] in pathodict['serotype']:
                PathoSamdict['血清型'] = f'''{cpserodb['亚型全称'].tolist()[0]}(重点关注)'''
            else:
                PathoSamdict['血清型'] = cpserodb['亚型全称'].tolist()[0]

    if species == 'vcholerae':
        pathodict = PathoNetdict[species]
        if os.path.isfile(f'{Pre}.vfdb.tsv') and os.path.getsize(f'{Pre}.vfdb.tsv') != 0:
            vchovfdb = pd.read_table(f'{Pre}.vfdb.tsv')
            tarvflist = [i for i in pathodict['vfgene'] if i in vchovfdb['基因名称'].tolist()]
            if tarvflist:
                PathoSamdict['毒力基因'] = ';'.join(tarvflist)
        if os.path.isfile(f'{Pre}_serotype_result.tsv') and os.path.getsize(f'{Pre}_serotype_result.tsv') != 0:
            cpserodb = pd.read_table(f'{Pre}_serotype_result.tsv')
            cpserodb = cpserodb.fillna('-')
            if cpserodb['血清型'].tolist()[0] in pathodict['serotype']:
                PathoSamdict['血清型'] = f'''{cpserodb['血清型'].tolist()[0]}(重点关注)'''
            else:
                PathoSamdict['血清型'] = cpserodb['血清型'].tolist()[0]
    if species == 'ecoli_achtman_4':
        ecodb = pd.read_table(f'{Pre}_serotype_result.tsv')
        if ecodb['物种'].tolist()[0] == 'Escherichia coli':
            pathodict = PathoNetdict['ecoli']
            stype = ecodb['O抗原'].tolist()[0]
            if stype in pathodict['serotype']:
                PathoSamdict['血清型'] = f'{stype}(重点关注)' 
            else:
                PathoSamdict['血清型'] = f'{stype}' 
        elif 'Shigella' in ecodb['物种'].tolist()[0]:
            pathodict = PathoNetdict['Shigella']
            stype = ecodb['志贺分型'].tolist()[0]
            if stype in pathodict['serotype']:
                PathoSamdict['血清型'] = f'{stype}(重点关注)'
            else:
                PathoSamdict['血清型'] = f'{stype}'
        PathoSamdict['物种'] = ecodb['物种'].tolist()[0]
        ecovfdb = pd.read_table(f'{Pre}.vfdb.tsv')
        tarvflist = [i for i in pathodict['vfgene'] if i in ecovfdb['基因名称'].tolist()] 
        if tarvflist:
            PathoSamdict['毒力基因'] = ';'.join(tarvflist)


    PathoSamdb = pd.DataFrame(PathoSamdict,index=[0])
    PathoSamdb.to_csv(f'{Pre}.pathonet_result.tsv',sep='\t',index=False)
        


def extract_SpeID(row):
     if pd.isna(row['FullLineageRanks']) or pd.isna(row['FullLineage']):
         return "noSpe"

     ranks = row['FullLineageRanks'].split(';')

     if 'species' not in ranks:
         return "noSpe"

     idx = ranks.index('species')
     lineage = row['FullLineageTaxIDs'].split(';')

     if idx < len(lineage):
         return lineage[idx]
     else:
         return "noSpe"

def extract_Spe(row):
    if pd.isna(row['FullLineageRanks']) or pd.isna(row['FullLineage']):
        return "noSpe"

    ranks = row['FullLineageRanks'].split(';')
    if 'species' not in ranks:
        return "noSpe"

    idx = ranks.index('species')
    lineage = row['FullLineage'].split(';')
    if idx < len(lineage):
        return lineage[idx]
    else:
        return "noSpe"

def is_non_numeric_in_bracket(x):
    if pd.isna(x):
        return False
    m = re.search(r'\((.*?)\)', str(x))
    if m:
        return not m.group(1).isdigit()
    return False


# mlst分型
def mlst_serotype(Pre,tSpe):
    #20260210
    kradb1 = pd.read_table(f'{Pre}_assem.kraken2.txt',header=None)
    #krspe =  kradb1.loc[kradb1[3]=='S',5].tolist()[0].strip()
    #krspeid =  kradb1.loc[kradb1[3]=='S',4].tolist()[0]
    cmdb = pd.read_table(f'{Pre}.checkm.tsv')
    Asdb = pd.read_table('Assem_info1.tsv')
    cfile = pytaxonkit.lineage(Asdb['taxid'].tolist())
    cfile['Species'] = cfile.apply(extract_SpeID, axis=1)
    cfile['SpeciesName'] = cfile.apply(extract_Spe, axis=1)
    tcfile = cfile.merge(Asdb,left_on='TaxID',right_on='taxid').drop_duplicates().groupby('Species').sum('序列长度').reset_index().sort_values('序列长度',ascending=False)
    krspeid = tcfile.loc[tcfile['Species']!='noSpe']['Species'].tolist()[0]
    #krspe = kradb1.loc[kradb1[4].astype('str')==krspeid,5].tolist()[0].strip()
    krspe = cfile.loc[cfile['Species']==krspeid,'SpeciesName'].tolist()[0]
    krspeidlist = cfile.loc[cfile['Species']==str(krspeid),'TaxID'].tolist()
    mainper = int(Asdb.loc[Asdb['taxid'].isin(krspeidlist),'序列长度'].sum()/Asdb.序列长度.sum()*100)
    if 'noSpe' in tcfile['Species'].tolist():
        noSpeper = round(tcfile.loc[tcfile['Species']=='noSpe','序列长度'].tolist()[0]/tcfile.序列长度.sum()*100,2)
        cmdb['物种名称'] = f'''{krspe}({mainper}%) noSpe({noSpeper}%)'''
    else:
        cmdb['物种名称'] = f'''{krspe}({mainper}%)'''
    cmdb.to_csv(f'{Pre}.checkm.tsv',sep='\t',index=False)

    subprocess.run(f"mlst --quiet --csv {Pre}.final.fasta > {Pre}_mlst.csv",shell=True) # 运行 mlst 工具进行分型
    mlst_gene = pd.read_table(f"{Pre}_mlst.csv",sep = ",",header=None) 
    if mlst_gene.shape[1] <= 4:  # mlst鉴定失败
        if tSpe:
            subprocess.run(f"mlst --scheme {tSpe} --quiet --csv {Pre}.final.fasta > {Pre}_mlst.csv",shell=True)
        else:
            subprocess.run(f"mlst --scheme ecoli --quiet --csv {Pre}.final.fasta > {Pre}_mlst.csv",shell=True)
    mlst_gene = pd.read_table(f"{Pre}_mlst.csv",sep = ",",header=None)
    sub = mlst_gene.iloc[:,3:]
    if sub.applymap(is_non_numeric_in_bracket).sum(axis=1).tolist()[0] < 3:
        mlst_B = mlst_gene.iloc[0,1]  # mlst_B = 第一行第二个（物种信息）
    else:
        mlst_B = 'UnKnown'
    mlst_st = mlst_gene.iloc[0,2] # mlst_st = 第一行第三个（序列分型ST）
    genes = mlst_gene.iloc[0,].tolist() # gene = 第一行
    for x in genes:
        x = str(x)
        if x.find("(") >=0:   # 若字符串 x 中是否包含左括号'('，则返回索引 >=0
            gene = x.split("(")[0]  # gene = 括号前的部分
            gene_num= x.split("(")[1].split(")")[0] # gene_num = 括号中的部分
            try:
                #gene_num = re.findall(r"\d+",gene_num)[0] # gene_num = 提取gene_num的数字部分
                os.system(f"seqkit grep -p {gene}_{gene_num} /home/dell/miniconda3/envs/TB_ONT/db/pubmlst/{mlst_B}/{gene}.tfa >> mlst_all.fa") # 根据{gene}_{gene_num}提取对应{gene}.tfa序列追加到 mlst_all.fa 文件中
            except:
                pass      
    with open('mlst.log','w') as mlstg:
        subprocess.run(f"dnadiff mlst_all.fa {Pre}.final.fasta",shell=True,stdout=mlstg,stderr=mlstg)  # 运行dnadiff比较样本fasta和mlst官方基因序列的差异
        subprocess.run(f"show-coords -lTH out.delta|sort -k7nr > mlst.coords",shell=True,stdout=mlstg,stderr=mlstg)   # mlst.coords = 运行show-coords提取 dnadiff生成的out.delta文件中的比对信息
    if os.path.isfile('mlst.coords') and os.path.getsize('mlst.coords') != 0:
        mlst_gene = pd.read_table("mlst.coords",header=None) # 起始位置，参考序列比对的起始位置； 比对起始位置，查询序列比对的起始位置； 基因长度：参考序列的长度； 比对长度：比对的片段在参考序列中所覆盖的长度
        mlst_gene.columns= ["起始位置","终止位置","比对起始位置","比对终止位置","基因长度","比对长度","一致性%","基因长度","序列长度","管家基因","序列名称"] 
        sss = mlst_gene.pop("管家基因")
        mlst_gene.insert(0,"管家基因",sss) # 调整'管家基因'列的位置到第一列
        sss = mlst_gene.pop("序列名称")
        mlst_gene.insert(3,"序列名称",sss) # 调整'序列名称'列的位置到第四列
        mlst_gene["序列分型(ST)"]= mlst_st
        mlst_gene["物种信息"]=mlst_B   # 添加"序列分型(ST)"和"物种信息"列
        sss = mlst_gene.pop("序列长度") # 将 "序列长度" 列提取到'sss',并删除原表格列
        mlst_gene['管家基因序号'] = mlst_gene['管家基因'].str.split('_').str[1]
        mlst_gene = mlst_gene[['管家基因','管家基因序号','起始位置','终止位置','序列名称','比对起始位置','比对终止位置','基因长度','比对长度','一致性%','序列分型(ST)','物种信息']]
        mlst_gene.to_csv(f"{Pre}.mlst_Stat.txt",index=0,sep='\t')  # 结果输出到{Pre}.mlst_Stat.txt        
        ##添加首行给前端确定字符间距
        #c_contig= mlst_gene["序列名称"].tolist()
        c_gene = mlst_gene["管家基因"].tolist()  # c_gene = '管家基因'列的值
        for x in c_gene:
            y = mlst_gene[mlst_gene["管家基因"]==x ]["序列名称"].tolist()[0] # y = mlst_gene 中与 x 相同的 "管家基因" 对应的 "序列名称"
            os.system(f"show-aligns out.delta {x} {y}|sed -n '/-- Alignments/,/--   END/p' > {x}_gene_show.txt") # 提取out.delta中{x}和{y}比对结果的 -- Alignments 到 -- 
    
    cdb = pd.read_table(f'{Pre}.checkm.tsv')
    cdb['mlst 物种名称'] = mlst_B
    print(mlst_B,Pre)
    cdb[['样本名称','物种名称','mlst 物种名称','污染率','完整性']].to_csv(f'{Pre}.checkm.tsv',sep = '\t',index=False)
    if mlst_B == 'bordetella_3':  # 百日咳疫苗基因型——BLASTn
        bp_vaccine(Pre)       
        if os.path.isfile(f'{Pre}.R1.fastq.gz'):
            bp_2037(Pre)
        else:
            fa_2037(Pre)
        bp_mlva(Pre) 
    elif mlst_B == 'klebsiella': # 肺炎克雷伯菌血清型——kleborate
        serotype_kb(Pre)
    elif 'ecoli' in mlst_B:      # 大肠杆菌及志贺血清型——ectyper.在mlst鉴定中两者名称一致     
        serotype_B(Pre)
    elif 'senterica' in mlst_B:  # 沙门氏菌血清型——sistr
        serotype_A(Pre)
    elif 'hinfluenzae' in mlst_B:  # 流感嗜血亚型预测——HICAP
        serotype_HI(Pre)
    #elif 'neisseria' in mlst_B:   # 奈瑟氏菌亚型预测——PMGA          
        #serotype_nm(Pre)
    elif 'vparahaemolyticus' in mlst_B:   # 副溶血弧菌亚型预测——VPsero       
        serotype_D(Pre)
    elif 'spyogenes' in mlst_B:          # 化脓链球菌亚型预测——emm_typing.py
        serotype_groupA(Pre)
    elif 'vcholerae' in mlst_B:            # 霍乱亚型预测——vista    
        serotype_E(Pre)
    #elif 'brucella' in mlst_B:         # 布鲁氏菌亚型预测——mlva  
        #serotype_MLVA(Pre)
    elif 'saureus' in mlst_B:        # 金葡亚型预测——salty    
        serotype_st(Pre)
    elif 'listeria' in mlst_B:     # 单增李斯特亚型预测——lissero       
        serotype_lm(Pre)
    elif 'bcereus' in mlst_B:   # 芽孢杆菌亚型预测——PMGA          
        serotype_bt(Pre)
    elif 'mpneumoniae' in mlst_B:
        mp_mlva(Pre)
    elif 'ssuis' in mlst_B:    #猪链球菌
        serotype_SS(Pre)
    elif 'campylobacter' in mlst_B:   #空肠弯曲菌
        serotype_Cb(Pre)
    if os.path.isfile(f'{Pre}_2.report.txt'):
        kradb = pd.read_table(f'{Pre}_2.report.txt',header=None)
    else:
        kradb = pd.read_table(f'{Pre}_assem.kraken2.txt',header=None)
    if kradb.loc[kradb[3]=='G',5].tolist()[0].strip() == 'Yersinia':
        serotype_ys(Pre)
    if tSpe:    # 前期有强制指定
        PathoNet(Pre,tSpe)
    else:
        PathoNet(Pre,mlst_B)

# 定义AnnoEle函数：元件预测等注释
def AnnoEle(Pre,threads):  # 1.minced-CRISPR 2.PhySpy-噬菌体 3.Dimob.pl-基因岛 #4.signalp6-信号肽 #5.bakta—ncRNA 6.mefinder - fasta移动元件 7.diamond blastp - uniq.faa移动元件  8.R - GO富集分析
    # 1.minced-CRISPR
    with open('AnnoElog','w') as f1:
        subprocess.run(f'minced {Pre}_prokka/{Pre}.fna {Pre}_CRISPRs.txt {Pre}_CRISPRs.gff',shell=True,stderr=f1,stdout=f1)  # 使用 minced 工具对 CRISPR 进行预测
        subprocess.run(f'''grep 'minced:0.4.2' {Pre}_CRISPRs.gff > {Pre}_CRISPRs.tsv''',shell=True)   # 从 {Pre}_CRISPRs.gff 中筛选出包含 'minced:0.4.2' 的行，保存到{Pre}_CRISPRs.tsv
        CRIdb = pd.read_table(f'{Pre}_CRISPRs.tsv',header=None,names=['序列ID','软件版本','片段类型','开始位置','终止位置','得分','a','b','结果注释']) # 添加列名
        CRIdb = CRIdb[['序列ID','软件版本','片段类型','开始位置','终止位置','得分','结果注释']]  # 保留需要的列
        CRIdb.to_csv(f'{Pre}.CRISPR.tsv',sep='\t',index=False)  # 保存为 {Pre}.CRISPR.tsv
        ignum = int(os.popen(f'''cat {Pre}_prokka/{Pre}.gff|grep '##' |wc -l ''').read().strip())-1 # 计算需要跳过的行
        afile = pd.read_table(f'{Pre}_prokka/{Pre}.gff',skiprows=ignum,low_memory=False,header=None) # 跳过这些行后读取
        afile = afile.loc[~afile[3].isna()]   # 去除第4列（基因组特征的起始位置）为空的行
        repeatfile = afile[afile[2]=='repeat_region']  # 筛选重复区域repeat_region
        repeatfile = afile[[0,1,2,3,4,8]]
        repeatfile.rename(columns={0:'Contig名称',1:'数据库',2:'类型',3:'序列开始',4:'序列结尾',8:'结果注释'},inplace=True) # 添加列名
        repeatfile.to_csv(f'{Pre}.repeat.tsv',sep='\t',index=False)  # 保存到{Pre}.repeat.tsv
    # 2.PhySpy-噬菌体
    with open('physpy.log','a') as phyf:
        subprocess.run(f'PhiSpy.py {Pre}_prokka/{Pre}.gbk -o {Pre}_PhySpy --threads {threads}',shell=True,stdout=phyf,stderr=phyf) # 运行PhiSpy.py - 噬菌体
    if int(os.popen(f'cat {Pre}_PhySpy/prophage_coordinates.tsv|wc -l').read().strip()) >=1:  # 统计行数至少＞1
        bfile = pd.read_table(f'{Pre}_PhySpy/prophage_coordinates.tsv',header=None)
        bfile.rename(columns={0:'噬菌体名称',1:'Contig名称',2:'Contig起始位置',3:'Contig终止位置',4:'attL起始位置',5:'attL终止位置',6:'attR起始位置',7:'attR终止位置',8:'attL序列',9:'attR序列',10:'原因'},inplace=True)
        bfile.to_csv(f'{Pre}.phi.tsv',index=False,sep='\t') # 添加列名后，保存到{Pre}.phi.tsv
    # 3.Dimob.pl-基因岛
    with open('AnnoElog','a') as f1:
        subprocess.run(f'/data/deploy/meta_genome/soft/islandpath/Dimob.pl {Pre}_prokka/{Pre}.gbk {Pre}',shell=True,stdout=f1,stderr=f1)  # 运行Dimob.pl - 基因岛
        isdb = pd.read_table(f'{Pre}_annot.tsv')
        isdb.rename(columns={'##GI_id':'GI号','sequence':'序列ID','start':'起始位置','end':' 终止位置','strand':'正负链','orf_name':'开放阅读框','annotation':'结果注释'},inplace=True)
        isdb.to_csv(f'{Pre}.annot.tsv',sep='\t',index=False) # 添加列名后，保存到{Pre}.annot.tsv
        # 4.signalp6-信号肽
        #subprocess.run(f'signalp6 -fasta {Pre}_prokka/{Pre}.faa --output_dir {Pre}_sigIP --bsize 20',shell=True,stdout=f1,stderr=f1)   # 运行signalp6 - 信号肽
        #sigdb = pd.read_table(f'{Pre}_sigIP/region_output.gff3',skiprows=1,sep='\t',header=None) # 跳过第一行
        #sigdb = sigdb[[0,2,3,4]]
        #sigdb['faaID'] = sigdb[0].str.split(' ').str[0]
        #sigdb['序列信息'] = sigdb[0].str.split(' ').str[1:].str.join(' ')
        #sigdb.rename(columns={2:'信号肽区域',3:'预测起始位置',4:'预测终止位置'},inplace=True)
        #sigdb1 = pd.read_table(f'{Pre}_sigIP/prediction_results.txt',skiprows=1,sep='\t')
        #sigdb1['faaID'] = sigdb1['# ID'].str.split(' ').str[0]
        #sigdb = sigdb.merge(sigdb1,on='faaID')
        #sigdb.rename(columns={'Prediction':'预测信号肽'},inplace=True)
        #sigdb = sigdb[['faaID','预测起始位置','预测终止位置','序列信息','信号肽区域','预测信号肽']]
        #sigdb.sort_values(['faaID','预测起始位置'],ascending=True,inplace=True)
        #sigdb.to_csv(f'{Pre}.sigIP.tsv',sep='\t',index=False)  # 保存到{Pre}.sigIP.tsv
    # 5.bakta—ncRNA(非编码RNA)
    #with open('bakta.log','w') as bakf:  # 运行bakta进行基因注释，ncRNA（非编码 RNA）
        #subprocess.run(f'/home/dell/miniconda3/bin/conda run -n GTDBtk bakta --db /data/deploy/meta_genome/soft/baktk/db-light --threads {threads} {Pre}.final.fasta --output {Pre}_batka --prefix {Pre}',shell=True,stdout=bakf,stderr=bakf)
        #if int(os.popen(f'cat {Pre}_batka/{Pre}.gff3|grep ncRNA|wc -l').read()) > 0:
            #subprocess.run(f'grep ncRNA {Pre}_batka/{Pre}.gff3 > {Pre}.ncRNAt.tsv',shell=True)
            #ncdb = pd.read_table(f'{Pre}.ncRNAt.tsv',names=['序列名称','数据库','鉴定类别','起始位置','终止位置','evalue','正负链','tmp1','Info'])
            #ncdb['基因名称'] = ncdb['Info'].str.split(';').str[3].str.replace('gene=','')
            #ncdb['基因详情'] = ncdb['Info'].str.split(';').str[1].str.replace('Name=','')
            #ncdb = ncdb[['序列名称','数据库','鉴定类别','起始位置','终止位置','evalue','正负链','基因名称','基因详情']]
            #ncdb.to_csv(f'{Pre}.ncRNA.tsv',sep='\t',index=False)  # 保存到{Pre}.ncRNA.tsv 
    # 6.mefinder - fasta移动元件 
    try:
        with open('AnnoElog','a') as f1:
            subprocess.run(f'mefinder find --contig {Pre}.final.fasta  {Pre} -t {threads}',shell=True,stdout=f1,stderr=f1) # 运行mefinder
            medb = pd.read_table(f'{Pre}.csv',sep=',',skiprows=5)
            medb = medb[['contig','start','end','name','type','allele_len','e_value','identity','coverage']]
            medb.rename(columns={'name':'元件名称','type':'元件类型','allele_len':'元件长度','identity':'一致性','coverage':'覆盖度','contig':'序列名称','start':'起始位置','end':'终止位置'},inplace=True)
            medb.to_csv(f'{Pre}.medb.tsv',sep='\t',index=False)  # 添加列名后，保存到{Pre}.medb.tsv
    except:
        print('不是完整基因组，移动元件鉴定失败')       
    # 7.diamond blastp - uniq.faa移动元件
    with open('mobileOG.log','w') as mbf:
        mgemt = pd.read_table('/data/deploy/meta_genome/database/beatrix/mobileOG-db-beatrix-1.6-All.csv',sep=',',low_memory=False)  # mgemt = 读取移动元件数据库文件,以','分隔
        subprocess.run(f'diamond blastp -q {Pre}.faa --db /data/deploy/meta_genome/database/beatrix/mobileOG-db  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --out {Pre}.mgeblast.tsv',shell=True,stdout=mbf,stderr=mbf)
        mgedb = pd.read_table(f'{Pre}.mgeblast.tsv',names=['序列名称','参考基因组名称','相似性(%)','长度','差异数量','空缺数量','序列起始','序列终止','参考起始','参考终止','evalue','比对得分'])
        mgedb = mgedb[(mgedb['相似性(%)']>25) & (mgedb['evalue']<1e-5)]  # 保留相似性大于 25% 且 evalue 小于 1e-5 的比对结果
        mgedb['ID'] = mgedb['参考基因组名称'].str.split('|').str[0]
        mgedb['类型'] = mgedb['参考基因组名称'].str.split('|').str[3:].str.join('|')
        mgedb = mgedb.merge(mgemt,left_on='ID',right_on='mobileOG Entry Name') # mobileOG-db-beatrix-1.6-All.csv数据库信息 与 {Pre}.mgeblast.tsv按照ID合并
        mgedb = mgedb[['序列名称','参考基因组名称','类型','Manual Annotation','Name','相似性(%)','长度','差异数量','空缺数量','序列起始','序列终止','参考起始','参考终止','evalue','比对得分']]
        mgedb.rename(columns={'Manual Annotation':'注释','Name':'元件名称'},inplace=True) 
        mgedb.to_csv(f'{Pre}.medb.tsv',sep='\t',index=False) # 添加列名后，保存到{Pre}.medb.tsv
    # 8. R - GO富集分析  
    with open('enrichplot.log','w') as enf:
        subprocess.run(f'/home/dell/miniconda3/bin/conda run -n report_env Rscript /data/deploy/meta_genome/GO.R {os.getcwd()}',shell=True,stdout=enf,stderr=enf) # 使用 Rscript 运行 R 脚本，GO 富集分析

def combine_func(Pre):

    result_dir = Path(f"{Pre}_genome_complete_result")
    anno_dir = Path(f"{Pre}_anno_sum")

    with open("combine.log", "w") as logf:

        # ------------------------------------------------
        # 清理旧结果
        # ------------------------------------------------
        if result_dir.exists():
            shutil.rmtree(result_dir)

        result_dir.mkdir()

        # ------------------------------------------------
        # 1 Data QC summary
        # ------------------------------------------------
        copy_pattern(
            ["*QC*.summary.tsv", "*.fastp*.json"],
            result_dir / "1.DataSum"
        )

        copy_pattern(
            ["*qc"],
            result_dir / "1.DataSum"
        )

        # ------------------------------------------------
        # 2 Species identification
        # ------------------------------------------------
        copy_pattern(
            ["*.list.txt", "*.krona.html"],
            result_dir / "2.Spereads"
        )

        # ------------------------------------------------
        # 3 Assembly
        # ------------------------------------------------
        copy_pattern(
            [
                "flye_output/assembly_info.txt",
                "*.checkm.tsv",
                "*_raw.png",
                "*.final.fasta"
            ],
            result_dir / "3.Assemble"
        )

        # ------------------------------------------------
        # 4 Repeat
        # ------------------------------------------------
        copy_pattern(
            ["*.repeat.tsv", "*.mummer.*"],
            result_dir / "4.Repeat"
        )

        # ------------------------------------------------
        # 5 Functional element
        # ------------------------------------------------
        copy_pattern(
            [
                "*.rgi.tsv",
                "*.card.tsv",
                "*.vfdb.tsv",
                "*.phi.tsv",
                "*.CRISPR.tsv",
                "*.annot.tsv",
                "*.sigIP.tsv",
                "*.medb.tsv",
                "*.ncRNA.tsv"
            ],
            result_dir / "5.Fun_Element"
        )

        # ------------------------------------------------
        # 7 MLST
        # ------------------------------------------------
        copy_pattern(
            ["*mlst_Stat.txt", "*gene_show.txt"],
            result_dir / "7.Mlst"
        )

        # ------------------------------------------------
        # Annotation summary
        # ------------------------------------------------
        copy_pattern(
            [
                "*.swiss.tsv",
                "*.GO.tsv",
                "*.KEGG.tsv",
                "*.PFAM.tsv",
                "*.Cog.tsv",
                "*.CAZy.tsv",
                "*.rgi.tsv",
                "*.card.tsv",
                "*.vfdb.tsv",
                "*.phi.tsv",
                "*.CRISPR.tsv",
                "*.annot.tsv",
                "*.sigIP.tsv",
                "*.medb.tsv",
                "*.ncRNA.tsv"
            ],
            anno_dir
        )

        run_cmd(
            f"tar -zcf {Pre}.anno.tar.gz {anno_dir}",
            logf
        )

        # ------------------------------------------------
        # R report
        # ------------------------------------------------
        if method != "meta":

            run_cmd(
                f"""
                /home/dell/miniconda3/bin/conda run -n report_env \
                Rscript /data1/shanghai_pip/meta_genome/report_asb.R \
                {os.getcwd()}
                """,
                logf
            )

        else:

            run_cmd(
                f"""
                /home/dell/miniconda3/bin/conda run -n report_env \
                Rscript /data/deploy/meta_genome/report_db/report_meta.R \
                {os.getcwd()}
                """,
                logf
            )

        # ------------------------------------------------
        # meta_out
        # ------------------------------------------------
        meta_dir = Path("meta_out")
        meta_dir.mkdir(exist_ok=True)

        copy_pattern(
            [
                "Summary_kraken.csv",
                f"{Pre}*list*.txt",
                "2.vfdb.tsv",
                "2.card.tsv",
                "summary.tsv",
                "test.fastp*.json"
            ],
            meta_dir
        )

        # ------------------------------------------------
        # 删除中间文件
        # ------------------------------------------------
        run_cmd(
            "rm -rf 2.listID*.txt 2_fqID.txt 2.*.sorted.bam "
            "2.id.tsv *_t.R*.fastq.gz *_sub.R*.fastq *regions.bed",
            logf
        )

        # ------------------------------------------------
        # reads rename
        # ------------------------------------------------
        if Path("2.1.fastq").exists():
            Path("2.1.fastq").rename(f"{Pre}D2.1.fastq")
            run_cmd(f"pigz -p 10 {Pre}D2.1.fastq", logf)

        if Path("2.2.fastq").exists():
            Path("2.2.fastq").rename(f"{Pre}D2.2.fastq")
            run_cmd(f"pigz -p 10 {Pre}D2.2.fastq", logf)

# 定义format_seconds函数：将秒转化成/小时/分钟/秒
def format_seconds(seconds):
    total_seconds = int(seconds)     # 将秒数保留整数部分
    hours = total_seconds // 3600
    minutes = (total_seconds % 3600) // 60   # %代表取余数
    seconds = total_seconds % 60
    return f"{hours}小时{minutes}分钟{seconds}秒"

# 主函数
# =========================================================
# 基础工具函数
# =========================================================

def ensure_sample_result_file():
    result_fp = f"{ofn}/sample_result.txt"
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
    result_tsv = f"{pre}.assemble.result.tsv"
    if os.path.isfile(result_tsv):
        subprocess.run(
            f"sed -n '2p' {result_tsv} >> {ofn}/sample_result.txt",
            shell=True
        )


def file_exists(path):
    return bool(path) and str(path) != "0" and os.path.isfile(path)


def is_fastq_input(infile, fastq1, fastq2):
    return is_fastq(infile) or is_fastq(fastq1) or is_fastq(fastq2)


def get_flow_list():
    return [x.strip() for x in runflow.split(",") if x.strip()]


def get_read_mode(pre):
    """
    返回：
        ont       : 只有三代
        ont_se    : 三代 + 二代单端
        ont_pe    : 三代 + 二代双端
        se        : 只有二代单端
        pe        : 只有二代双端
        none      : 未识别到reads
    """
    ont = os.path.isfile(f"{pre}.final.fastq")
    r1 = os.path.isfile(f"{pre}.R1.fastq.gz")
    r2 = os.path.isfile(f"{pre}.R2.fastq.gz")

    if ont:
        if r1 and r2:
            return "ont_pe"
        elif r1:
            return "ont_se"
        else:
            return "ont"

    if r1 and r2:
        return "pe"
    elif r1:
        return "se"
    else:
        return "none"


def infer_species(pre, llid):
    try:
        if llid != "nolevel":
            return level2Spe(pre, llid.split(",")[0])
    except Exception:
        pass
    return 0


# =========================================================
# fake 流程
# =========================================================

def run_fake_pipeline(pre, sam_num, all_num):
    print_progress(0, 5, sam_num, all_num, pre, "开始假流程分析")

    subprocess.run(
        "cp -r /data1/shanghai_pip/meta_genome/fakefile/Sample1/* ./",
        shell=True
    )
    subprocess.run(
        f"""rename 's/Sample1/{pre}/g' *""",
        shell=True
    )
    subprocess.run(
        f"""rename 's/Sample1/{pre}/g' */*""",
        shell=True
    )

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


# =========================================================
# fastq 主线
# =========================================================

def run_fastq_qc_and_assembly(
    infile, fastq1, fastq2, threads, llid, mm, pre,
    ml, mq, asmtype, pts, pst, rnalib, ref, gtf,
    flowlist, sam_num, all_num, start_time
):
    """
    返回:
        flowstep, flownum, anchor_time
    """
    if not os.path.isfile(f"{pre}.final.fasta"):
        flownum = len(flowlist) + 1
        flowstep = 1

        if "基因组组装" in flowlist:
            flowlist.remove("基因组组装")

        if not os.path.isfile("QC_ok"):
            QC_func(infile, fastq1, fastq2, mq, ml, pre, rnalib, threads, method)

        if "基因组组装" in runflow:
            if os.path.isfile(f"{pre}.R2.fastq.gz"):
                asb_func(
                    f"{pre}.final.fastq",
                    f"{pre}.R1.fastq.gz",
                    f"{pre}.R2.fastq.gz",
                    threads, pre, llid, pts, pst, method, asmtype, ref, gtf
                )
            elif os.path.isfile(f"{pre}.R1.fastq.gz"):
                asb_func(
                    f"{pre}.final.fastq",
                    f"{pre}.R1.fastq.gz",
                    0,
                    threads, pre, llid, pts, pst, method, asmtype, ref, gtf
                )
            else:
                asb_func(
                    f"{pre}.final.fastq",
                    0,
                    0,
                    threads, pre, llid, pts, pst, method, asmtype, ref, gtf
                )

            if method != "meta":
                Annotate_func(pre, threads)

        t1 = time.time()
        runtime = format_seconds(t1 - start_time)
        print_progress(flowstep, flownum, sam_num, all_num, pre, "数据质控&组装结束", runtime)
        return flowstep, flownum, t1

    else:
        flownum = len(flowlist)
        flowstep = 0
        t1 = time.time()
        return flowstep, flownum, t1


def run_species_identification(pre, threads):
    mode = get_read_mode(pre)

    if mode == "ont":
        print("三代鉴定")
        sys.stdout.flush()
        kk2(f"{pre}.final.fastq", 0, 0, threads, pre)

    elif mode == "ont_pe":
        print("三代鉴定+二代双端鉴定")
        sys.stdout.flush()
        kk2(f"{pre}.final.fastq", f"{pre}.R1.fastq.gz", f"{pre}.R2.fastq.gz", threads, pre)

    elif mode == "ont_se":
        print("三代鉴定+二代单端鉴定")
        sys.stdout.flush()
        kk2(f"{pre}.final.fastq", f"{pre}.R1.fastq.gz", 0, threads, pre)

    elif mode == "pe":
        print("二代双端鉴定")
        sys.stdout.flush()
        kk2(0, f"{pre}.R1.fastq.gz", f"{pre}.R2.fastq.gz", threads, pre)

    elif mode == "se":
        print("二代单端鉴定")
        sys.stdout.flush()
        kk2(0, f"{pre}.R1.fastq.gz", 0, threads, pre)

    else:
        print(f"{pre} 未检测到可用于物种鉴定的 reads")
        sys.stdout.flush()


def run_fastq_flow(flow, pre, threads, llid):
    if flow == "物种鉴定":
        if not os.path.isfile("kk2_ok"):
            run_species_identification(pre, threads)

    elif flow == "结构变异检测":
        if not os.path.isfile("SV_ok"):
            DrugFinder(pre, threads)

    elif flow == "功能注释":
        AnnoFun(pre, threads)

    elif flow == "元件预测":
        AnnoEle(pre, threads)

    elif flow == "耐药与毒力":
        VFDR(pre, threads)

    elif flow == "mlst与血清型":
        species = infer_species(pre, llid)
        mlst_serotype(pre, species)


def run_fastq_pipeline(
    infile, fastq1, fastq2, threads, llid, mm, pre,
    ml, mq, asmtype, all_num, sam_num,
    pts, pst, rnalib, ref, gtf,
    flowlist, start_time
):
    flowstep, flownum, last_time = run_fastq_qc_and_assembly(
        infile, fastq1, fastq2, threads, llid, mm, pre,
        ml, mq, asmtype, pts, pst, rnalib, ref, gtf,
        flowlist, sam_num, all_num, start_time
    )

    for flow in flowlist:
        flowstep += 1
        step_start = time.time()

        try:
            run_fastq_flow(flow, pre, threads, llid)
        except Exception as e:
            print(f"{pre} {flow} 失败: {e}")
            sys.stdout.flush()

        runtime = format_seconds(time.time() - step_start)
        print_progress(flowstep, flownum, sam_num, all_num, pre, f"{flow}已结束", runtime)
        last_time = time.time()

    flowstep += 1
    combine_func(pre)

    total_runtime = format_seconds(time.time() - start_time)
    if "基因组组装" in runflow:
        append_sample_result(pre)

    print_progress(flowstep, flownum, sam_num, all_num, pre, "数据合并结束", total_runtime)


# =========================================================
# fasta 主线
# =========================================================

def run_fasta_init(pre, threads, flowlist, sam_num, all_num, start_time):
    """
    返回:
        rt_flag, flowstep, flownum
    """
    prokka_dir = f"{ofn}/fastq_analysis/{pre}/{pre}_prokka"

    if not os.path.isdir(prokka_dir):
        rt_flag = 1
        flownum = len(flowlist) + 2
        flowstep = 1

        Annotate_func(pre, threads)

        t1 = time.time()
        runtime = format_seconds(t1 - start_time)
        print_progress(flowstep, flownum, sam_num, all_num, pre, "基因预测已完成", runtime)
        return rt_flag, flowstep, flownum

    else:
        rt_flag = 0
        flowstep = 0
        flownum = len(flowlist) + 1
        return rt_flag, flowstep, flownum


def run_fasta_flow(flow, pre, threads, llid):
    if flow == "功能注释":
        AnnoFun(pre, threads)

    elif flow == "元件预测":
        AnnoEle(pre, threads)

    elif flow == "耐药与毒力":
        VFDR(pre, threads, "fasta")

    elif flow == "mlst与血清型":
        species = infer_species(pre, llid)
        mlst_serotype(pre, species)


def run_fasta_pipeline(
    infile, threads, llid, mm,
    pre, all_num, sam_num,
    flowlist, start_time
):
    rt_flag, flowstep, flownum = run_fasta_init(
        pre, threads, flowlist, sam_num, all_num, start_time
    )

    for flow in flowlist:
        flowstep += 1
        step_start = time.time()

        try:
            run_fasta_flow(flow, pre, threads, llid)
        except Exception as e:
            print(f"{pre} {flow} 失败: {e}")
            sys.stdout.flush()

        runtime = format_seconds(time.time() - step_start)
        print_progress(flowstep, flownum, sam_num, all_num, pre, f"{flow}已结束", runtime)

    flowstep += 1
    combine_func(pre)

    total_runtime = format_seconds(time.time() - start_time)
    if rt_flag:
        append_sample_result(pre)

    print_progress(flowstep, flownum, sam_num, all_num, pre, "数据合并结束", total_runtime)


# =========================================================
# 主函数
# =========================================================

def main_process(
    infile, fastq1, fastq2, threads, llid, mm, Pre,
    ml, mq, asmtype, Allnum, Samnum, pts, pst,
    rnalib, iffake=tmpfake, ref=ref, gtf=gtf
):
    start_time = time.time()

    if iffake:
        run_fake_pipeline(Pre, Samnum, Allnum)
        return

    ensure_sample_result_file()

    flowlist = get_flow_list()
    total_for_start = len(flowlist) + 1
    print_progress(0, total_for_start, Samnum, Allnum, Pre, "开始分析")

    if is_fastq_input(infile, fastq1, fastq2):
        run_fastq_pipeline(
            infile, fastq1, fastq2, threads, llid, mm, Pre,
            ml, mq, asmtype, Allnum, Samnum,
            pts, pst, rnalib, ref, gtf,
            flowlist, start_time
        )
    else:
        run_fasta_pipeline(
            infile, threads, llid, mm,
            Pre, Allnum, Samnum,
            flowlist, start_time
        )


# 定义 check_input 函数：确定输入的文件类型
def check_input(infile,intp):   # 输入文件的绝对路径（inf）， 输入文件类型（intype）       
    if intp == 'barcode_fastq':   # 已拆分,需要cat在一起的数据放在同一个文件夹内
        if os.path.isdir(infile): # 检查输入的是否为文件夹
            for i in os.listdir(infile):
                if i.startswith('barcode'): # 检查文件名称是否以'barcode'开头
                    return 'bardir'
                    break      # 找到一个以'barcode'开头的文件就结束循环
            else:
                print('文件夹下未检测到barcode文件夹')
                sys.exit()      
        else:
            print('请输入barcode_fastq文件夹')
            sys.exit()    
    elif intp == 'fastq':
        #if barkit == 'none':
        if os.path.isdir(infile):  # fqdir:fq文件夹
            return 'fqdir'
        elif os.path.isfile(infile):  # fqfile:fq文件
            return 'fqfile'
        #else:
        #    if os.path.isdir(infile):
        #        return 'barfqdir'
        #    elif os.path.isfile(infile):
        #        return 'barfqfile'
    
    elif '.cfg' in intp:      
        if os.path.isdir(infile):
            return 'f5dir'
        else:
            print('请输入f5dir文件夹')
            sys.exit()           
    elif '@' in intp:
        if os.path.isdir(infile):
            return 'pod5'
        else:
            print('请输入pod5文件夹')
            sys.exit()           
    elif intp == 'fasta':
        if os.path.isdir(infile):
            return 'fadir'         # fadir:fa文件夹
        else:
            return 'fafile'        # fafile:fa文件

# 定义get_free_gpu_memory函数：使用 NVIDIA 的 nvidia-smi 工具获取可用的 GPU 内存量
def get_free_gpu_memory():
    try:
        result = subprocess.run(["nvidia-smi", "--query-gpu=memory.free", "--format=csv,noheader,nounits"], stdout=subprocess.PIPE) # stdout=subprocess.PIPE用于捕获命令的标准输出
        memory_free = result.stdout.decode('utf-8').strip()  # 去除空白字符
        if memory_free:
            return memory_free
        else:
            return "noGPU"
    except Exception as e:
        return "noGPU"
        
# 定义basecaller函数：guppy_basecaller将fast5转变为fastq输出；dorado将pod5转变为fastq输出
def basecaller(inputtype,params,inf,ofn='basecaller_outputs'):  # 传入'fast5'或'pod5',intype,inf
    freem = get_free_gpu_memory()  # 获取可用GPU量
    if freem != 'noGPU':           # 如果我有GPU
        freem = int(freem)
        if not os.path.isdir(ofn):
            os.makedirs(ofn)
        # 如果我传入的inf是一个文件，将他变为文件夹（处理的必须是个文件夹），在本流程中不会被调用    
        if os.path.isfile(inf):  
            os.makedirs(f'{ofn}/tmp')
            subprocess.run(f'cp {inf} {ofn}/tmp',shell=True)  # 将 inf 所指定的文件复制到 ofn 目录下的 tmp 子目录中
            inf = f'{ofn}/tmp'                               # 为inf重新赋值
            
        if inputtype == 'fast5':                           
            print('下机数据为fast5模式')
            print(f'guppy_basecaller -r -i {inf} -s {ofn} -x auto -c {params}')
            sys.stdout.flush()
            if freem/1.2 > 5500:   # 根据可用GPU量，利用--chunk_size 1000控制调用量
                subprocess.run(f'guppy_basecaller -r -i {inf} -s {ofn} -x auto -c {params}',shell=True)
            else:
                subprocess.run(f'guppy_basecaller -r -i {inf} -s {ofn} -x auto -c {params} --chunk_size 1000',shell=True)
            subprocess.run(f'seqkit seq {ofn}/pass/*.fastq > {ofn}/basecaller.fastq',shell=True)   # guppy_basecaller会生成一个pass文件夹，一个未通过的文件夹

        else:
            print('下机数据为pod5模式')  
            qdict = {'fast@':8,'hac@':9,'sup@':10}   # 键fast@，值8
            configdir = '/home/dell/biosoft/dorado'  
            paramlist = params.split('_')            # 以'_ '为分隔符，将 params 字符串分割成一个列表
            #1.判断文件格式 2.输出。 dorado不会生成pass文件夹，只有guppy才会生成
            if len(paramlist) == 4:  # 如果他有四个元素
                seqtype,stype,v1,modeln = paramlist   # 将 paramlist 中的四个元素分别赋值给 seqtype、stype、v1 和 modeln 这四个变量
                qthod = qdict.get(modeln)
                subprocess.run(f'{configdir}/bin/dorado basecaller -r {configdir}/{params}v3.3 {inf} --min-qscore {qthod} --emit-fastq > {ofn}/basecaller.fastq',shell=True)  # --min-qscore {qthod}：设置了最小的质量分数（qthod）参数
            else:
                seqtype,stype,v1,cspeed,modeln = paramlist
                qthod = qdict.get(modeln)
                if cspeed == '400bps':
                    subprocess.run(f'{configdir}/bin/dorado basecaller -r {configdir}/{params}v4.1.0 {inf} --min-qscore {qthod} --emit-fastq  > {ofn}/basecaller.fastq',shell=True)
                    if os.path.getsize(f'{ofn}/basecaller.fastq') == 0:
                        subprocess.run(f'{configdir}/bin/dorado basecaller -r {configdir}/{params}v4.2.0 {inf} --min-qscore {qthod} --emit-fastq  > {ofn}/basecaller.fastq',shell=True)
                else:
                    subprocess.run(f'{configdir}/bin/dorado basecaller -r {configdir}/{params}v4.1.0 {inf} --min-qscore {qthod} --emit-fastq  > {ofn}/basecaller.fastq',shell=True)

legacy_asb_func = asb_func

from metagenomic_refactor.assembly import (
    asb_func,
    outfun,
    register_assembly_callbacks,
    rgi_fun,
)
from metagenomic_refactor.common import (
    basecaller,
    check_input,
    copy_pattern,
    file_exists,
    format_seconds,
    get_free_gpu_memory,
    is_fasta,
    is_fastq,
    is_fastq_input,
    run_cmd,
)
from metagenomic_refactor.context import RuntimeContext, set_runtime_context
from metagenomic_refactor.qc import (
    QC_func,
    get_logger,
    ngs_qc,
    normalize_fastqc_images,
    read_fastqc_data,
    safe_read_json,
    summaryfastqc_prod,
)
from metagenomic_refactor.report import combine_func
from metagenomic_refactor.taxonomy import (
    exreadsID1,
    getCovDep,
    getinfo,
    kk2,
    proc_kra1,
)
from metagenomic_refactor.workflow import (
    append_sample_result,
    ensure_sample_result_file,
    get_flow_list,
    get_read_mode,
    infer_species,
    main_process,
    print_progress,
    register_workflow_callbacks,
    run_fake_pipeline,
    run_fasta_flow,
    run_fasta_init,
    run_fasta_pipeline,
    run_fastq_flow,
    run_fastq_pipeline,
    run_fastq_qc_and_assembly,
    run_species_identification,
)

set_runtime_context(
    RuntimeContext(
        ofn=ofn,
        runflow=runflow,
        method=method,
        rmhost=rmhost,
        tspeabun=tspeabun,
        krdb=Krdb,
    )
)

register_assembly_callbacks(legacy_asb_func=legacy_asb_func)

register_workflow_callbacks(
    asb_func=asb_func,
    Annotate_func=Annotate_func,
    kk2=kk2,
    DrugFinder=DrugFinder,
    AnnoFun=AnnoFun,
    AnnoEle=AnnoEle,
    VFDR=VFDR,
    mlst_serotype=mlst_serotype,
    level2Spe=level2Spe,
)

#----                    main process                    ----  目前在ofn路径下，Snum = 开始分析第几个样本，Anum = 总共有多少个样本
# 三代数据+fasta：1.argv.input是一个文件且第一行以'@'开始，即fastq文件 2.是一个文件夹 3.调用is_fasta函数            
if os.path.isfile(argv.input) and os.popen(f'less {argv.input}|head -n 1').read().strip().startswith('@') or os.path.isdir(argv.input) or is_fasta(argv.input):
    protype = check_input(inf,intype)   # 调用check_input函数，输入：输入文件的绝对路径（inf）和 输入文件类型（intype）
    print(protype)  
    sys.stdout.flush()    
    
    if protype == 'fqdir': # 单个样本的三代数据文件夹，三代数据可以边测序边分析因此可以有多个FQ文件                                                                         
        Sam = 'sample1'    # 用于创建工作目录名称
        Anum = 1
        Snum = 1
        os.makedirs(f'fastq_analysis/{Sam}')  # 在当前工作目录下创建一个新的目录
        open(f'{ofn}/fastq_analysis/Samplelist.txt','a').write(f'{Sam}\n')  # 以追加模式，将sam值一个一行写入Samplelist.txt
        open('Sample_result.txt','w').write(f'样本名称\t序列数量\t碱基数量\t运行时间\n')  # 以写入模式，创建（覆盖）Sample_result.txt，并写入
        os.chdir(f'fastq_analysis/{Sam}')  # 切换工作目录。f'使得字符串中可以直接嵌入变量和表达式
        wkdir = os.getcwd()                # 获取当前工作目录的绝对路径，储存到wkdir
        subprocess.run(f'cat {inf}/*.f*q* |seqkit rmdup -i |seqkit seq > {Sam}.raw.fastq',shell=True)  # seqkit rmdup -i，去重保留第一个序列ID（考虑到三代测序中途获取FQ文件，后续又再次选择全部FQ文件造成数据重复的情况）
        if os.path.getsize(f'{Sam}.raw.fastq') == 0:  # 如果文件大小为0
            print('文件夹内没有三代fastq格式文件')
            sys.exit()
        try:
            print('三代fastq文件夹分析模式')
            print(f'{Sam}.raw.fastq',fastq1,fastq2,nt,mmethod,Sam,minl,minQ,asm_type,Anum,Snum,ptimes,psoft)
            sys.stdout.flush()
            main_process(f'{Sam}.raw.fastq',fastq1,fastq2,nt,llid,mmethod,Sam,minl,minQ,asm_type,Anum,Snum,ptimes,psoft,rnalib,tmpfake,ref,gtf)  # 调用并传递参数给main_process函数
        except:
            print('fqdir分析失败')
        os.chdir(ofn)  
    elif protype == 'fqfile':   # 单个三代数据文件
        Sam = 'sample1'
        Anum = 1
        Snum = 1
        if not os.path.isdir(f'fastq_analysis/{Sam}'):
            os.makedirs(f'fastq_analysis/{Sam}')
        open(f'{ofn}/fastq_analysis/Samplelist.txt','a').write(f'{Sam}\n')
        open('Sample_result.txt','w').write(f'样本名称\t序列数量\t碱基数量\t运行时间\n')
        #--测试假流程脚本---
        os.chdir('fastq_analysis') 
        if tmpfake:
            subprocess.run(f'cp -r /data1/shanghai_pip/meta_genome/fake_result/* ./',shell=True)  # -r，复制目录及其内容
        os.chdir(ofn)
        #---
        os.chdir(f'fastq_analysis/{Sam}')
        wkdir = os.getcwd()
        #subprocess.run(f'cat {inf}|seqkit seq > {Sam}.raw.fastq',shell=True)
        subprocess.run(f'ln -s {inf} ./{Sam}.raw.fastq',shell=True)        # ln -s ，创建一个名为{Sam}.raw.fastq的符号链接指向{inf}。 -s，软链接，可以跨盘；-h，硬链接，删除原始文件不会失效
        if os.path.getsize(f'{Sam}.raw.fastq') == 0:
            print('输入三代fastq格式文件有误，请核对')
            sys.exit()
        try:
            print('三代fastq文件分析模式')
            print(f'{Sam}.raw.fastq',fastq1,fastq2,nt,mmethod,Sam,minl,minQ,asm_type,Anum,Snum,ptimes,psoft)
            sys.stdout.flush()
            main_process(f'{Sam}.raw.fastq',fastq1,fastq2,nt,llid,mmethod,Sam,minl,minQ,asm_type,Anum,Snum,ptimes,psoft,rnalib,tmpfake,ref,gtf)
        except:
            print('fqfile分析失败')
        os.chdir(ofn)

    elif protype == 'f5dir':
        os.makedirs(f'fastq_analysis')
        open('Sample_result.txt','w').write(f'样本名称\t序列数量\t碱基数量\t运行时间\n')
        os.chdir('fastq_analysis')
        wkdir = os.getcwd()
        basecaller('fast5',intype,inf)    # 调用basecaller函数：guppy_basecaller将fast5转变为fastq输出
        #subprocess.run(f'/home/dell/biosoft/ont-guppy/bin/guppy_basecaller -r -i {inf} -s fastq_out -x auto -c dna_r9.4.1_450bps_hac.cfg',shell=True)  # basecaller函数已包含其功能，多了GPU的合理调用
        if barkit != 'none':  # 如果需要拆分试剂盒
            subprocess.run(f'/home/dell/biosoft/ont-guppy/bin/guppy_barcoder -i basecaller_outputs/pass -s barout -x auto --barcode_kits {barkit}',shell=True) 
            Anum = len([i for i in os.listdir('barout')])    # Anum = 存储在barout中被拆分出来的样本数
            Snum = 1
            for i in os.listdir('barout'):
                if i.startswith('barcode'):
                    os.makedirs(i)         
                    os.chdir(f'{wkdir}/{i}')
                    Sam = i
                    open(f'{ofn}/fastq_analysis/Samplelist.txt','a').write(f'{Sam}\n')
                    subprocess.run(f'cat {ofn}/fastq_analysis/barout/{i}/*.f*q*|seqkit seq > {Sam}.raw.fastq',shell=True)  # 将 {i}目录下的所有 fastq 文件合并
                    try:
                        main_process(f'{Sam}.raw.fastq',fastq1,fastq2,nt,llid,mmethod,Sam,minl,minQ,asm_type,Anum,Snum,ptimes,psoft,tmpfake,ref,gtf)
                        Snum+=1    # Snum的值＋1
                    except:
                        print(f'f5dir分析失败')
                        print(f'task_step：{flowstep}/{flownum}\t样本进度：{Snum}/{Anum}\t样本：{Sam}\t数据分析中断')
                        Snum+=1
                    os.chdir(wkdir)   
            os.chdir(ofn)
        else:
            Sam = 'sample1'
            Anum = 1
            Snum = 1
            os.makedirs(Sam)
            os.chdir(Sam)
            open(f'{ofn}/fastq_analysis/Samplelist.txt','a').write(f'{Sam}\n')
            open('Sample_result.txt','w').write(f'样本名称\t序列数量\t碱基数量\t运行时间\n')
            subprocess.run(f'cat {ofn}/fastq_analysis/basecaller_outputs/pass/*.f*q*|seqkit seq > {Sam}.raw.fastq',shell=True)
            try:
                main_process(f'{Sam}.raw.fastq',fastq1,fastq2,nt,llid,mmethod,Sam,minl,minQ,asm_type,Anum,Snum,ptimes,psoft,rnalib,tmpfake,ref,gtf)
            except:
                print('数据量不足')
            os.chdir(ofn)
            subprocess.run(f'cat */*/Sample_result.txt >> Sample_result.txt',shell=True)     # 合并所有的Sample_result.txt文件
            
    elif protype == 'bardir':   # barcode_fastq
        os.makedirs(f'fastq_analysis')
        open('Sample_result.txt','w').write(f'样本名称\t序列数量\t碱基数量\t运行时间\n')
        os.chdir('fastq_analysis')
        wkdir = os.getcwd()
        Anum = len([i for i in os.listdir(inf) if i.startswith('barcode')])  # Anum= 存储在inf中以'barcode'开头的样本数
        Snum = 1
        print('输入文件为barcode文件夹')
        for i in os.listdir(inf):
            if i.startswith('barcode'):
                Sam = i        
                try:
                    os.makedirs(i)
                    os.chdir(f'{wkdir}/{i}')       
                    open(f'{ofn}/fastq_analysis/Samplelist.txt','a').write(f'{Sam}\n')
                    if len([fq for fq in os.listdir(f'{inf}/{i}') if is_fastq(f'{inf}/{i}/{fq}')]):  # 遍历 {inf}/{i} 目录中的所有文件和目录，将符合is_fastq函数的值进行返回，返回的值长度＞0
                        subprocess.run(f'cat {inf}/{i}/*.f*q*|seqkit rmdup -i |seqkit seq > {Sam}.raw.fastq',shell=True)   # seqkit rmdup -i , 移除重复的序列
                    else:
                        print('hello')
                        subprocess.run(f'touch {Sam}.raw.fastq',shell=True)  # touch：新建一个空文件夹或刷新原文件夹
                    main_process(f'{Sam}.raw.fastq',fastq1,fastq2,nt,llid,mmethod,Sam,minl,minQ,asm_type,Anum,Snum,ptimes,psoft,rnalib,tmpfake,ref,gtf) 
                    Snum+=1
                    os.chdir(wkdir)
                except:
                    print(f'bardir{i}数据量不足,或者基因组组装失败')
                    print(f'task_step：{flowstep}/{flownum}\t样本进度：{Snum}/{Anum}\t样本：{Sam}\t数据分析中断')
                    Snum+=1  
                    os.chdir(wkdir)                    
        os.chdir(ofn)         
    
    elif protype == 'barfqdir':   # 未拆分的bardir，暂未启用
        os.makedirs(f'fastq_analysis')
        open('Sample_result.txt','w').write(f'样本名称\t序列数量\t碱基数量\t运行时间\n')
        os.chdir('fastq_analysis')
        wkdir = os.getcwd()
        os.makedirs('fastq_out')
        subprocess.run(f'/home/dell/biosoft/ont-guppy/bin/guppy_barcoder -i {inf} -s barout -x auto --barcode_kits {barkit}',shell=True)
        Anum = len([i for i in os.listdir('barout') if i.startswith('barcode')])
        Snum = 1
        for i in os.listdir('barout'):
            if i.startswith('barcode'):
                Sam = i
                try:
                    os.makedirs(i)
                    os.chdir(f'{wkdir}/{i}')
                    open(f'{ofn}/fastq_analysis/Samplelist.txt','a').write(f'{Sam}\n')
                    if len([fq for fq in os.listdir(f'{ofn}/fastq_analysis/barout/{i}/') if is_fastq(f'{ofn}/fastq_analysis/barout/{i}/{fq}')]):
                        subprocess.run(f'cat {ofn}/fastq_analysis/barout/{i}/*.f*q*|seqkit seq > {Sam}.raw.fastq',shell=True)
                    else:
                        print('hello')
                        subprocess.run(f'touch {Sam}.raw.fastq',shell=True)  # touch：新建一个空文件夹或刷新原文件夹
                    main_process(f'{Sam}.raw.fastq',fastq1,fastq2,nt,llid,mmethod,Sam,minl,minQ,asm_type,Anum,Snum,ptimes,psoft,rnalib,tmpfake,ref,gtf)
                    Snum+=1
                    os.chdir(wkdir)
                except:
                    os.chdir(wkdir)
                    print(f'barfqdir{i}数据量不足,或者基因组组装失败')
                    print(f'task_step：{flowstep}/{flownum}\t样本进度：{Snum}/{Anum}\t样本：{Sam}\t数据分析中断')
                    Snum+=1
        os.chdir(ofn)
    elif protype == 'barfqfile':  # 暂未启用
        os.makedirs(f'fastq_analysis')
        open('Sample_result.txt','w').write(f'样本名称\t序列数量\t碱基数量\t运行时间\n')
        os.chdir('fastq_analysis')
        wkdir = os.getcwd()
        os.makedirs('fastq_out')
        subprocess.run(f'ln -s {inf} fastq_out',shell=True)   #  subprocess.run(f'cp {inf} fastq_out',shell=True)
        if barkit != 'none':
            subprocess.run(f'/home/dell/biosoft/ont-guppy/bin/guppy_barcoder -i fastq_out/ -s barout -x auto --barcode_kits {barkit}',shell=True)
            Anum = len([i for i in os.listdir('barout') if i.startswith('barcode')])
            Snum = 1
            for i in os.listdir('barout'):
                if i.startswith('barcode'):
                    Sam = i
                    try:
                        os.makedirs(i)
                        os.chdir(f'{wkdir}/{i}')
                        open(f'{ofn}/fastq_analysis/Samplelist.txt','a').write(f'{Sam}\n')
                        subprocess.run(f'cat {ofn}/fastq_analysis/barout/{i}/*.f*q*|seqkit rmdup -i |seqkit seq > {Sam}.raw.fastq',shell=True)            
                        main_process(f'{Sam}.raw.fastq',fastq1,fastq2,nt,llid,mmethod,Sam,minl,minQ,asm_type,Anum,Snum,ptimes,psoft,rnalib,tmpfake,ref,gtf)
                        Snum+=1
                        os.chdir(wkdir)
                    except:
                        os.chdir(wkdir)
                        print(f'{i}数据量不足,或者基因组组装失败')
                        print(f'task_step：{flowstep}/{flownum}\t样本进度：{Snum}/{Anum}\t样本：{Sam}\t数据分析中断')
                        Snum+=1
        os.chdir(ofn)

    elif protype == 'fadir':  
        os.makedirs(f'fastq_analysis') 
        open('Sample_result.txt','w').write(f'样本名称\t序列数量\t碱基数量\t运行时间\n')
        os.chdir('fastq_analysis')
        wkdir = os.getcwd()
        faendlist = ['fa','fasta','fna','fas','fa.gz','fasta.gz','fna.gz','fas.gz']
        Anum = 0   
        for tfile in os.listdir(inf):
            endssuf = [i for i in faendlist if tfile.endswith(i)]  # endssuf = 为i赋值faendlist中的值，如果tfile中的值以i结尾，这个i值会被提取出来
            if sum(endssuf) > 0:
                Anum+=1                # 计数总共有几个符合后缀名称的文件
        Snum = 1 
        for tfile in os.listdir(inf):
            endssuf = [i for i in faendlist if tfile.endswith(i)]
            if sum(endssuf) > 0:
                Sam = tfile.replace(f'.{endssuf[0]}','')  # sam = tfile的 '.'+'endssuf'列表中第一个值 被替换为'空格'，即只保留文件名称前面部分
                open(f'{ofn}/fastq_analysis/Samplelist.txt','a').write(f'{Sam}\n')
                os.makedirs(f'{wkdir}/{Sam}')
                os.chdir(f'{wkdir}/{Sam}')
                subprocess.run(f'seqkit seq -w 0  {inf}/{tfile} > {Sam}.final.fasta',shell=True)  # 指定输出序列的每行宽度默认为60,-w 0 表示不对序列换行
                try:
                    main_process(f'{inf}/{tfile}',fastq1,fastq2,nt,llid,mmethod,Sam,minl,minQ,asm_type,Anum,Snum,ptimes,psoft,rnalib,tmpfake,ref,gtf)
                except:
                    print(f'fadir{i}数据量不足,或者基因组组装失败')
                    print(f'task_step：{flowstep}/{flownum}\t样本进度：{Snum}/{Anum}\t样本：{Sam}\t数据分析中断')
                Snum+=1
                os.chdir(ofn)
    elif protype == 'fafile':
        os.makedirs(f'fastq_analysis')
        open('Sample_result.txt','w').write(f'样本名称\t序列数量\t碱基数量\t运行时间\n')
        os.chdir('fastq_analysis')
        wkdir = os.getcwd()
        Anum = 1
        Snum = 1
        Sam = 'sample1' 
        open(f'{ofn}/fastq_analysis/Samplelist.txt','a').write(f'{Sam}\n')
        os.makedirs(f'{wkdir}/{Sam}')
        os.chdir(f'{wkdir}/{Sam}')
        subprocess.run(f'seqkit seq -w 0  {inf} > {Sam}.final.fasta',shell=True)
        main_process(inf,fastq1,fastq2,nt,llid,mmethod,Sam,minl,minQ,asm_type,Anum,Snum,ptimes,psoft,rnalib,tmpfake,ref,gtf)
        os.chdir(ofn)
    
    elif protype == 'pod5':
        if not os.path.isdir('fastq_analysis'):
            os.makedirs(f'fastq_analysis')
        open('Sample_result.txt','w').write(f'样本名称\t序列数量\t碱基数量\t运行时间\n') 
        os.chdir('fastq_analysis')
        wkdir = os.getcwd()
        if not os.path.isdir(f'{ofn}/fastq_analysis/basecaller_outputs/'): 
            basecaller('pod5',intype,inf)
        if barkit != 'none':
            subprocess.run(f'/home/dell/biosoft/ont-guppy/bin/guppy_barcoder -i basecaller_outputs/*.f*q*-s barout -x auto --barcode_kits {barkit}',shell=True)  
            Anum = len([i for i in os.listdir('barout')])
            Snum = 1
            for i in os.listdir('barout'):
                if i.startswith('barcode'):
                    os.makedirs(i)
                    os.chdir(f'{wkdir}/{i}')
                    Sam = i
                    open(f'{ofn}/fastq_analysis/Samplelist.txt','a').write(f'{Sam}\n')
                    if os.path.isdir(f'{ofn}/fastq_analysis/basecaller_outputs/'):
                        subprocess.run(f'cat {ofn}/fastq_analysis/barout/{i}/*.f*q*|seqkit seq > {Sam}.raw.fastq',shell=True)
                    try:
                        main_process(f'{Sam}.raw.fastq',fastq1,fastq2,nt,llid,mmethod,Sam,minl,minQ,asm_type,Anum,Snum,ptimes,psoft,rnalib,tmpfake,ref,gtf)
                        Snum+=1
                    except:
                        print(f'{i}数据量不足,或者基因组组装失败')
                        print(f'task_step：{flowstep}/{flownum}\t样本进度：{Snum}/{Anum}\t样本：{Sam}\t数据分析中断')
                        Snum+=1
                    os.chdir(wkdir)
            os.chdir(ofn)
        else:    # pod5单个样本
            Sam = 'sample1'
            Anum = 1
            Snum = 1
            os.makedirs(Sam)
            os.chdir(Sam)
            open(f'{ofn}/fastq_analysis/Samplelist.txt','a').write(f'{Sam}\n')
            open('Sample_result.txt','w').write(f'样本名称\t序列数量\t碱基数量\t运行时间\n')
            if os.path.isdir(f'{ofn}/fastq_analysis/basecaller_outputs/'):
                subprocess.run(f'cat {ofn}/fastq_analysis/basecaller_outputs/*.f*q*|seqkit seq > {Sam}.raw.fastq',shell=True)
            try:
                main_process(f'{Sam}.raw.fastq',fastq1,fastq2,nt,llid,mmethod,Sam,minl,minQ,asm_type,Anum,Snum,ptimes,psoft,rnalib,tmpfake,ref,gtf)
            except:
                print('数据量不足')
            os.chdir(ofn)
            subprocess.run(f'cat */*/Sample_result.txt >> Sample_result.txt',shell=True)

else:  # 传入list文本，fastq数据。list文本从左往右依次为：样本名称、三代数据、二代数据左、二代数据右、物种信息。
       # 注，三代数据列为：三代+fasta;每一行：三代数据传入的只能是单个样本文件/文件夹，fasta数据可以是多个样本组成的文件夹
    if not os.path.isdir('fastq_analysis'):
        os.makedirs(f'fastq_analysis')
    print('list分析模式')
    os.chdir('fastq_analysis')
    wkdir = os.getcwd()
    Snum = 1 
    listdb = pd.read_table(inf) 
    listdb['样本名称'] = listdb['样本名称'].astype('str')  # 将listdb的样本名称列的内容转化为字符串格式。使纯数字的样本名称可行
    listdb.fillna(0,inplace=True)  # fillna:将缺失值填充为0；inplace=True,在原表格上进行
    Anum = listdb.shape[0]         # shape[0]获取行数；shape[1]获取列数
    print(listdb)
    for line in listdb.index.tolist():  # listdb.index:行索引（不包括列名）; tolist：将 Index 对象转换为一个普通的 Python 列表，即[0,1,2]
        print(line)
        sys.stdout.flush()
        Pre =  listdb.iloc[line,].样本名称   
        tinf = listdb.iloc[line,].三代数据
        fastq1 = listdb.iloc[line,].二代数据左
        fastq2 = listdb.iloc[line,].二代数据右
        llid = listdb.iloc[line,].物种信息
       
        ## 三代fastq
        print(Pre,tinf,fastq1,fastq2,llid)       
        if tinf and intype == 'fastq':   # 如果thif有数据且文件类型为fastq
            if not (tinf.endswith('fa') or tinf.endswith('fas') or tinf.endswith('fna') or tinf.endswith('fasta') or tinf.endswith('fa.gz') or tinf.endswith('fas.gz') or  tinf.endswith('fasta.gz') or tinf.endswith('fna.gz')):  # 如果tinf不以这些结尾
                protype = check_input(tinf,'fastq')   # 调用check_input函数
                print(protype)
                try:
                    if protype == 'fqdir':   # 三代数据文件夹
                        Sam = Pre   
                        open(f'{ofn}/fastq_analysis/Samplelist.txt','a').write(f'{Sam}\n')
                        if not os.path.isdir(f'{wkdir}/{Sam}'):
                            os.makedirs(f'{wkdir}/{Sam}')
                        os.chdir(f'{wkdir}/{Sam}')
                        subprocess.run(f'cat {tinf}/*.f*q* |seqkit rmdup -i |seqkit seq > {Sam}.raw.fastq',shell=True)  # 合并序列
                        if os.path.getsize(f'{Sam}.raw.fastq') == 0:
                            print(f'{Sam}文件夹内没有三代fastq格式文件')
                            sys.exit()
                        main_process(f'{Sam}.raw.fastq',fastq1,fastq2,nt,llid,mmethod,Sam,minl,minQ,asm_type,Anum,Snum,ptimes,psoft,rnalib,tmpfake,ref,gtf)
                        Snum+=1
                        os.chdir(ofn)
                    elif protype == 'fqfile':
                        Sam = Pre
                        open(f'{ofn}/fastq_analysis/Samplelist.txt','a').write(f'{Sam}\n')
                        if not os.path.isdir(f'{wkdir}/{Sam}'):
                            os.makedirs(f'{wkdir}/{Sam}')
                        open('Sample_result.txt','w').write(f'样本名称\t序列数量\t碱基数量\t运行时间\n')
                        os.chdir(f'{wkdir}/{Sam}')
                        subprocess.run(f'cat {tinf}|seqkit seq > {Sam}.raw.fastq',shell=True)
                        if os.path.getsize(f'{Sam}.raw.fastq') == 0:
                            print(f'{Sam}三代fastq格式文件有误，请核对')
                            sys.exit()
                        main_process(f'{Sam}.raw.fastq',fastq1,fastq2,nt,llid,mmethod,Sam,minl,minQ,asm_type,Anum,Snum,ptimes,psoft,rnalib,tmpfake,ref,gtf)
                        Snum+=1
                        os.chdir(ofn)
                except:
                    os.chdir(ofn)
                    print(f'{Pre}数据量不足,或者基因组组装失败')
                    print(f'task_step：{flowstep}/{flownum}\t样本进度：{Snum}/{Anum}\t样本：{Pre}\t数据分析中断')
                    Snum+=1
            
            ##fasta
            else:      # 和三代数据共用 intype == 'fastq'，可以更好地做到fasta+fastq数据的一起上传
                protype = check_input(tinf,'fasta')
                if protype == 'fadir':
                    Anum = 0
                    faendlist = ['fa','fasta','fna','fas','fa.gz','fasta.gz','fna.gz','fas.gz']
                    for tfile in os.listdir(inf):
                        endssuf = [i for i in faendlist if tfile.endswith(i)]   # endssuf = 为i赋值faendlist中的值，如果tfile中的值以i结尾，这个i值会被提取出来。
                        if sum(endssuf) > 0:
                            Anum+=1
                    Snum = 1
                    for tfile in os.listdir(inf):
                        endssuf = [i for i in faendlist if tfile.endswith(i)]
                        if sum(endssuf) > 0:
                            Sam = tfile.replace(f'.{endssuf[0]}','')
                            open(f'{ofn}/fastq_analysis/Samplelist.txt','a').write(f'{Sam}\n')
                            if not os.path.isdir(f'{wkdir}/{Sam}'):
                                os.makedirs(f'{wkdir}/{Sam}')
                            open('Sample_result.txt','w').write(f'样本名称\t序列数量\t碱基数量\t运行时间\n')
                            os.chdir(f'{wkdir}/{Sam}')
                            subprocess.run(f'seqkit seq -i -w 0  {tinf} > {Sam}.final.fasta',shell=True)
                            try:
                                main_process(f'{Sam}.final.fasta',fastq1,fastq2,nt,llid,mmethod,Sam,minl,minQ,asm_type,Anum,Snum,ptimes,psoft,rnalib,tmpfake,ref,gtf)
                            except:
                                print(f'{Sam}数据量不足,或者基因组组装失败')
                                print(f'task_step：{flowstep}/{flownum}\t样本进度：{Snum}/{Anum}\t样本：{Sam}\t数据分析中断')
                            Snum+=1
                            os.chdir(ofn)
                elif protype == 'fafile':                
                    Sam = Pre
                    open(f'{ofn}/fastq_analysis/Samplelist.txt','a').write(f'{Sam}\n')
                    if not os.path.isdir(f'{wkdir}/{Sam}'):
                        os.makedirs(f'{wkdir}/{Sam}')
                    open('Sample_result.txt','w').write(f'样本名称\t序列数量\t碱基数量\t运行时间\n')
                    os.chdir(f'{wkdir}/{Sam}')
                    subprocess.run(f'seqkit seq -i -w 0 {tinf} > {Sam}.final.fasta',shell=True)
                    main_process(f'{Sam}.final.fasta',fastq1,fastq2,nt,llid,mmethod,Sam,minl,minQ,asm_type,Anum,Snum,ptimes,psoft,rnalib,tmpfake,ref,gtf)
                    os.chdir(ofn)
                    Snum +=1
        
        ##二代fastq
        else:   # tinf = 0
            Sam = Pre
            open(f'{ofn}/fastq_analysis/Samplelist.txt','a').write(f'{Sam}\n')
            if not os.path.isdir(f'{wkdir}/{Sam}'):
                os.makedirs(f'{wkdir}/{Sam}')
            os.chdir(f'{wkdir}/{Sam}')
            subprocess.run(f'touch {Sam}.raw.fastq',shell=True)
            main_process(0,fastq1,fastq2,nt,llid,mmethod,Sam,minl,minQ,asm_type,Anum,Snum,ptimes,psoft,rnalib,tmpfake,ref,gtf)
            Snum+=1
            os.chdir(ofn)            
print('所有分析进程已经完成')
