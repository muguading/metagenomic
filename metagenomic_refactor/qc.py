from __future__ import annotations

import json
import logging
import math
import os
import shutil
import subprocess
import time
from dataclasses import dataclass
from typing import Any, Dict, Optional

import pandas as pd

from metagenomic_refactor.context import get_runtime_context


def get_logger(name="qc", level=logging.INFO):
    logger = logging.getLogger(name)
    if not logger.handlers:
        logger.setLevel(level)
        h = logging.StreamHandler()
        fmt = logging.Formatter("[%(asctime)s] %(levelname)s %(name)s: %(message)s")
        h.setFormatter(fmt)
        logger.addHandler(h)
    return logger


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
    cur = d
    for k in keys:
        if not isinstance(cur, dict) or k not in cur:
            return default
        cur = cur[k]
    return cur


@dataclass(frozen=True)
class QCTarget:
    name: str
    trigger_dir: str
    raw_fastqc_dir: str
    final_fastqc_dir: str
    json_path: str
    mode: str
    out_tsv: str


def build_qc_table(raw_fastqc, final_fastqc, before, after) -> pd.DataFrame:
    qc = {"过滤前": {}, "过滤后": {}}
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
    if mode == "nano":
        before = pick(fpjson, "summary", "before_filtering", default={}) or {}
        after = pick(fpjson, "summary", "after_filtering", default={}) or {}
        return before, after

    if mode in ("read1", "read2"):
        before = fpjson.get(f"{mode}_before_filtering", {}) or {}
        after = fpjson.get(f"{mode}_after_filtering", {}) or {}
        for tag in (before, after):
            tb = tag.get("total_bases") or 0
            if "q20_rate" not in tag and tb:
                tag["q20_rate"] = tag.get("q20_bases", 0) / tb
            if "q30_rate" not in tag and tb:
                tag["q30_rate"] = tag.get("q30_bases", 0) / tb
        return before, after

    logger.error(f"Unknown mode: {mode}")
    return {}, {}


def run_one_target(pre: str, t: QCTarget, logger: logging.Logger) -> bool:
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

    normalize_fastqc_images(t.raw_fastqc_dir, logger)
    normalize_fastqc_images(t.final_fastqc_dir.format(pre=pre), logger)
    normalize_summary_txt(os.path.join(t.raw_fastqc_dir, "summary.txt"), logger)
    normalize_summary_txt(os.path.join(t.final_fastqc_dir.format(pre=pre), "summary.txt"), logger)
    return True


def summaryfastqc_prod(pre: str, logger: Optional[logging.Logger] = None) -> Dict[str, bool]:
    logger = logger or get_logger("summaryfastqc", logging.INFO)
    targets = [
        QCTarget("Nano", "Nano_qc", "raw_fastqc", "{pre}.final_fastqc", "{pre}.filter.fastp.json", "nano", "NanoFastqc.tsv"),
        QCTarget("R1", "R1_qc", "raw.R1_fastqc", "{pre}.R1_fastqc", "{pre}.fastp2.json", "read1", "R1_Fastqc.tsv"),
        QCTarget("R2", "R2_qc", "raw.R2_fastqc", "{pre}.R2_fastqc", "{pre}.fastp2.json", "read2", "R2_Fastqc.tsv"),
    ]
    results = {}
    for t in targets:
        results[t.name] = run_one_target(pre, t, logger)
    return results


def proc_kra(kraken, tax, lel):
    tmplist = [tax]
    if lel in ["R", "D", "K", "P", "C", "O", "F", "G", "S"]:
        rawlist = ["R", "D", "K", "P", "C", "O", "F", "G", "S"]
    else:
        rawlist = ["S1", "S2", "S3", "S4", "S5", "S6"]
    if tax != 0:
        kradb = pd.read_table(kraken, header=None)
        tmpindex = kradb[(kradb[3] == lel) & (kradb[4] == tax)].index.tolist()[0] + 1
        if tmpindex <= kradb.shape[0] - 1:
            def getlindex(tmpindex):
                tmpl = kradb.iloc[tmpindex, 3]
                for tl in rawlist:
                    if tl in tmpl:
                        if tl == tmpl:
                            return rawlist.index(tl)
                        return rawlist.index(tl) + 1
                return 0

            while getlindex(tmpindex) > rawlist.index(lel):
                tmplist.append(kradb.iloc[tmpindex, 4])
                tmpindex += 1
    return tmplist


def exreadsID(Pre, taxlist, kraresult, fq1, fq2=0):
    Maintax = taxlist[0]
    kraredb = pd.read_table(kraresult, header=None)
    pd.DataFrame(kraredb[kraredb[2].isin(taxlist)][1].unique()).to_csv(
        f"{Maintax}_tfqID.txt", sep="\n", index=False, header=False
    )
    head_id = os.popen(f"head -n 1 {Maintax}_tfqID.txt").read().strip()
    if head_id.endswith("/1") or head_id.endswith("/2"):
        subprocess.run(f"cut -f1 -d '/' {Maintax}_tfqID.txt|sort -u|sed 's/$/\\/1/' > {Maintax}_fq1ID.txt", shell=True)
        subprocess.run(f"cut -f1 -d '/' {Maintax}_tfqID.txt|sort -u|sed 's/$/\\/2/' > {Maintax}_fq2ID.txt", shell=True)
    else:
        subprocess.run(f"ln -s {Maintax}_tfqID.txt {Maintax}_fq1ID.txt", shell=True)
        subprocess.run(f"ln -s {Maintax}_tfqID.txt {Maintax}_fq2ID.txt", shell=True)
    if str(fq1).endswith("final.fastq"):
        subprocess.run(f"seqkit grep -n -f {Maintax}_fq1ID.txt {fq1} > {Pre}_t.final.fastq", shell=True)
    else:
        subprocess.run(f"seqkit grep -n -f {Maintax}_fq1ID.txt {fq1} > {Pre}.R1.fastq", shell=True)
        subprocess.run(f"gzip -f  {Pre}.R1.fastq", shell=True)
        if fq2:
            subprocess.run(f"seqkit grep -n -f {Maintax}_fq2ID.txt {fq2} > {Pre}.R2.fastq", shell=True)
            subprocess.run(f"gzip -f {Pre}.R2.fastq", shell=True)


def ngs_qc(Pre):
    if not os.path.isfile("summary.tsv") or os.path.getsize("summary.tsv") == 0:
        subprocess.run("seqkit stat *_t.R*.fastq.gz -T -a -b > summary.tsv", shell=True)
    afile = pd.read_table("summary.tsv")
    afile = afile[["file", "format", "type", "num_seqs", "sum_len", "min_len", "avg_len", "max_len", "N50"]]
    afile["file"] = afile["file"].str.split(".").str[1]
    afile.rename(columns={"num_seqs": "序列总数", "sum_len": "碱基总数", "min_len": "最短序列长度", "avg_len": "平均序列长度", "max_len": "最长序列长度"}, inplace=True)
    afile.to_csv(f"{Pre}.QC2.summary.tsv", sep="\t", index=False)


def QC_func(inf, fq1, fq2, minq, minl, Pre, rnalib, threads, method):
    runtime = get_runtime_context()
    if inf:
        if method != "meta":
            subprocess.run(f"rasusa reads --bases 2gb -o {Pre}.sub.fastq {inf}", shell=True)
        else:
            subprocess.run(f"rasusa reads --bases 20gb -o {Pre}.sub.fastq {inf}", shell=True)
        subprocess.run(
            f"""fastp --in1 {Pre}.sub.fastq \
                    --out1 {Pre}.clean.fastq \
                    --thread {threads} \
                    --length_required=50 \
                    --n_base_limit=6 \
                    --compression=6 \
                    -Q \
                    -A \
                    --qualified_quality_phred=10 \
                    --json {Pre}.raw.fastp.json \
                    2> {Pre}.raw.fastp.log""",
            shell=True,
        )
        with open("QC.log", "w") as f:
            subprocess.run(f"porechop -i {Pre}.clean.fastq -o {Pre}.trim.fastq -t {threads} --no_split --check_reads 2000", shell=True, stdout=f, stderr=f)
        subprocess.run(f"seqkit seq -m {minl} -Q {minq} {Pre}.trim.fastq > {Pre}.filter.fastq", shell=True)
        if runtime.rmhost == "norm":
            subprocess.run(f"ln -s  {Pre}.filter.fastq {Pre}.rmhost.fastq", shell=True)
        else:
            subprocess.run(
                f"/home/dell/miniconda3/bin/conda run -n hostile hostile  clean --fastq1 {Pre}.filter.fastq --threads {threads} --aligner minimap2 --index /data/deploy/meta_new/Database/Host_Ref/hostile/human-t2t-hla.argos-bacteria-985_rs-viral-202401_ml-phage-202401.mmi > rmcont.json",
                shell=True,
            )
            subprocess.run(f"seqkit seq {Pre}.filter.clean.fastq.gz > {Pre}.rmhost.fastq", shell=True)
            subprocess.run(f"rm  {Pre}.filter.clean.fastq.gz", shell=True)
        subprocess.run(
            f"""fastp --in1  {Pre}.rmhost.fastq \
                    --out1 {Pre}.final.fastq \
                    --thread {threads} \
                    --length_required=50 \
                    --qualified_quality_phred=10 \
                    --n_base_limit=6 \
                    --compression=6 \
                    -Q \
                    -A \
                    --json {Pre}.filter.fastp.json \
                    2> {Pre}.filter.fastp.log""",
            shell=True,
        )
        subprocess.run("seqkit stat *.fastq -T -a -b > summary.tsv", shell=True)
        afile = pd.read_table("summary.tsv")
        afile = afile[["file", "format", "type", "num_seqs", "sum_len", "min_len", "avg_len", "max_len", "N50"]]
        afile["file"] = afile["file"].str.split(".").str[1]
        afile.index = afile["file"]
        afile.drop(["format", "type", "file"], axis=1, inplace=True)

        def get_meanQ(tlist):
            newlist = [(1 - 10 ** (-i / 10)) for i in tlist]
            return round(-10 * math.log10(1 - sum(newlist) / len(newlist)), 2)

        calqdict = {}
        for status in afile.index:
            subprocess.run(f"seqkit fx2tab -n -i -q {Pre}.{status}.fastq > tmp.tab", shell=True)
            calqdb = pd.read_table("tmp.tab", header=None)
            calqdict[status] = get_meanQ(calqdb[1].tolist())
        calqdb = pd.DataFrame(calqdict, index=["平均质量值"]).T
        calqdb["file"] = calqdb.index
        calqdb.reset_index(drop=True, inplace=True)
        afile.reset_index(inplace=True)
        afile = afile.merge(calqdb, on="file")
        afile.rename(columns={"num_seqs": "序列总数", "sum_len": "碱基总数", "min_len": "最短序列长度", "avg_len": "平均序列长度", "max_len": "最长序列长度"}, inplace=True)
        afile["file"] = afile["file"].replace({"raw": "原始序列", "clean": "去除过短/质量过差", "filter": "质量/长度过滤", "trim": "去除接头/barcode后序列", "rmhost": "去宿主后序列", "final": "最终序列"})
        order = ["原始序列", "去除过短/质量过差", "去除接头/barcode后序列", "质量/长度过滤", "去宿主后序列", "最终序列"]
        afile = afile.set_index("file").reindex(order).reset_index()
        afile = afile.rename(columns={"file": "过滤状态"})
        afile.to_csv(f"{Pre}.QC.summary.tsv", sep="\t", index=False)

        subprocess.run(f"ln -s {inf} ./raw.fastq", shell=True)
        if not os.path.isdir("Nano_qc"):
            os.makedirs("Nano_qc")
        with open("qc.log", "a") as f:
            subprocess.run(f"fastqc raw.fastq {Pre}.final.fastq -o Nano_qc -t {threads}", shell=True, stdout=f, stderr=f)
            subprocess.run("unzip -o 'Nano_qc/*.zip' ", shell=True, stdout=f, stderr=f)
        subprocess.run("rm raw.fastq", shell=True)

        with open("rawkk2.log", "w") as rawkkf:
            subprocess.run(f"kraken2 --db /data1/shanghai_pip/meta_genome/database/kraken2_custom_202101_24G --threads {threads} --output raw_{Pre}.out.txt --report raw_{Pre}.report.txt {Pre}.final.fastq", shell=True, stdout=rawkkf, stderr=rawkkf)
            subprocess.run(f"bracken -d /data1/shanghai_pip/meta_genome/database/kraken2_custom_202101_24G -o raw_{Pre}.bracken1.txt -w raw_{Pre}.bracken2.txt -l S -t 10  -i raw_{Pre}.report.txt", shell=True, stdout=rawkkf, stderr=rawkkf)
            tmpfile = pd.read_table(f"raw_{Pre}.bracken1.txt", sep="\t")
            tmpfile = tmpfile[tmpfile["name"] != "Homo"]
            ONTSpe = tmpfile.name.tolist()[0]
            tkid = tmpfile.taxonomy_id.tolist()[0]
            krakenfile = f"raw_{Pre}.report.txt"
            tkidl = os.popen(f"grep {tkid} {krakenfile}|cut -f4").read().strip()
            if float(tmpfile.loc[tmpfile["name"] == ONTSpe, "fraction_total_reads"].tolist()[0]) < float(runtime.tspeabun):
                taxlist1 = proc_kra(krakenfile, 517, "G")
                exreadsID(Pre, taxlist1, f"raw_{Pre}.out.txt", f"{Pre}.final.fastq")
                subprocess.run(f"mv {Pre}_t.final.fastq {Pre}.final.fastq", shell=True)

    if fq1 and fq2:
        if not os.path.isfile(f"{Pre}_sub.R1.fastq") or not os.path.isfile(f"{Pre}_sub.R2.fastq"):
            if method != "meta":
                subprocess.run(f"rasusa reads --bases 10gb -o {Pre}_sub.R1.fastq -o {Pre}_sub.R2.fastq {fq1} {fq2}", shell=True)
            else:
                subprocess.run(f"rasusa reads --bases 100gb -o {Pre}_sub.R1.fastq -o {Pre}_sub.R2.fastq {fq1} {fq2}", shell=True)
        subprocess.run(f"seqkit stat -T {Pre}_sub.R1.fastq {Pre}_sub.R2.fastq > raw_summary.tsv", shell=True)
        if not os.path.isfile(f"{Pre}.fastp2.json"):
            subprocess.run(
                f"""fastp --in1 {Pre}_sub.R1.fastq \
        --out1 {Pre}_t.R1.fastq.gz \
        --in2 {Pre}_sub.R2.fastq \
        --out2 {Pre}_t.R2.fastq.gz \
        --thread {threads} \
        --length_required=50 \
        --n_base_limit=6 \
        --compression=6 \
        --detect_adapter_for_pe \
        --json {Pre}.fastp2.json \
        2> {Pre}.fastp2.log""",
                shell=True,
            )
        ngs_qc(Pre)
        if not os.path.isdir("R1_qc"):
            os.makedirs("R1_qc")
        if not os.path.isdir("R2_qc"):
            os.makedirs("R2_qc")
        with open("qc.log", "a") as f:
            subprocess.run(f"ln -s  {Pre}_sub.R1.fastq  raw.R1.fastq", shell=True)
            subprocess.run(f"ln -s  {Pre}_sub.R2.fastq raw.R2.fastq", shell=True)
            if not os.path.isfile(f"{Pre}.R1.fastq.gz"):
                if runtime.rmhost == "norm":
                    subprocess.run(f"ln -s {Pre}_t.R1.fastq.gz  {Pre}.R1.fastq.gz;ln -s {Pre}_t.R2.fastq.gz  {Pre}.R2.fastq.gz", shell=True)
                else:
                    if rnalib == "0":
                        subprocess.run(f"/home/dell/miniconda3/bin/conda run -n hostile hostile clean --fastq1 {Pre}_t.R1.fastq.gz --threads {threads} --fastq2 {Pre}_t.R2.fastq.gz --aligner bowtie2 --index /data/deploy/meta_new/Database/Host_Ref/hostile/human-t2t-hla.argos-bacteria-985_rs-viral-202401_ml-phage-202401 > rmcont.json", shell=True, stdout=f, stderr=f)
                        subprocess.run(f"mv  {Pre}_t.R1.clean_1.fastq.gz  {Pre}.R1.fastq.gz", shell=True)
                        subprocess.run(f"mv  {Pre}_t.R2.clean_2.fastq.gz  {Pre}.R2.fastq.gz", shell=True)
                    else:
                        print("测试rna建库")
                        if not os.path.isfile(f"kneaddata_out/{Pre}_t.R1_kneaddata_paired_1.fastq"):
                            subprocess.run(f"/home/dell/miniconda3/bin/conda run -n kneaddata kneaddata --bypass-trf --bypass-trim -i1 {Pre}_t.R1.fastq.gz  -i2 {Pre}_t.R2.fastq.gz -o kneaddata_out -t {threads} -db /data/Ref/human_hg38_refMrna -db /data/Ref/SILVA_128_LSUParc_SSUParc_ribosomal_RNA", shell=True, stdout=f, stderr=f)
                        subprocess.run(f"pigz -p 10 kneaddata_out/{Pre}_t.R1_kneaddata_paired_1.fastq -c > {Pre}.R1.fastq.gz", shell=True, stdout=f, stderr=f)
                        subprocess.run(f"pigz -p 10 kneaddata_out/{Pre}_t.R1_kneaddata_paired_2.fastq -c > {Pre}.R2.fastq.gz", shell=True, stdout=f, stderr=f)
            if not os.path.isfile("R1_qc/raw.R1_fastqc.html"):
                subprocess.run(f"fastqc raw.R1.fastq {Pre}.R1.fastq.gz -o R1_qc -t {threads}", shell=True, stdout=f, stderr=f)
                subprocess.run(f"fastqc raw.R2.fastq {Pre}.R2.fastq.gz -o R2_qc -t {threads}", shell=True, stdout=f, stderr=f)
                subprocess.run("unzip -o 'R1_qc/*.zip' ", shell=True, stdout=f, stderr=f)
                subprocess.run("unzip -o 'R2_qc/*.zip' ", shell=True, stdout=f, stderr=f)
            subprocess.run(f"fastp -i {Pre}.R1.fastq.gz -I {Pre}.R2.fastq.gz -o tt.1.fq -O tt.2.fq -w {threads} --json {Pre}.final.json", shell=True, stdout=f, stderr=f)
            subprocess.run("rm tt.*.fq", shell=True)

    elif not fq2 and fq1:
        if not os.path.isfile(f"{Pre}_sub.R1.fastq"):
            subprocess.run(f"rasusa reads --bases 20gb -o {Pre}_sub.R1.fastq  {fq1}", shell=True)
        subprocess.run(f"seqkit stat -T {Pre}_sub.R1.fastq > raw_summary.tsv", shell=True)
        if not os.path.isfile(f"{Pre}_t.R1.fastq.gz"):
            subprocess.run(
                f"""fastp --in1  {Pre}_sub.R1.fastq \
        --out1 {Pre}_t.R1.fastq.gz \
        --thread {threads} \
        --length_required=50 \
        --qualified_quality_phred=10 \
        --n_base_limit=6 \
        --compression=6 \
        --json {Pre}.fastp2.json \
        2> {Pre}.fastp2.log""",
                shell=True,
            )
        ngs_qc(Pre)
        if not os.path.isdir("R1_qc"):
            os.makedirs("R1_qc")
        with open("qc.log", "a") as f:
            subprocess.run(f"ln -s  {Pre}_sub.R1.fastq raw.R1.fastq", shell=True)
            if runtime.rmhost == "norm":
                subprocess.run(f"ln -s {Pre}_t.R1.fastq.gz {Pre}.R1.fastq.gz", shell=True)
            else:
                if rnalib == "0":
                    subprocess.run(f"/home/dell/miniconda3/bin/conda run -n hostile hostile clean --fastq1 {Pre}_t.R1.fastq.gz --threads {threads}  --aligner bowtie2 --index /data/deploy/meta_new/Database/Host_Ref/hostile/human-t2t-hla.argos-bacteria-985_rs-viral-202401_ml-phage-202401 > rmcont.json", shell=True, stdout=f, stderr=f)
                    subprocess.run(f"mv  {Pre}_t.R1.clean.fastq.gz  {Pre}.R1.fastq.gz", shell=True)
                else:
                    print("测试rna建库")
                    if not os.path.isfile(f"kneaddata_out/{Pre}_t.R1_kneaddata_paired_1.fastq"):
                        subprocess.run(f"/home/dell/miniconda3/bin/conda run -n kneaddata kneaddata --bypass-trf --bypass-trim -i1 {Pre}_t.R1.fastq.gz  -o kneaddata_out -t {threads} -db /data/Ref/human_hg38_refMrna -db /data/Ref/SILVA_128_LSUParc_SSUParc_ribosomal_RNA", shell=True, stdout=f, stderr=f)
                    subprocess.run(f"pigz -p 10 kneaddata_out/{Pre}_t.R1_kneaddata_paired_1.fastq -c > {Pre}.R1.fastq.gz", shell=True, stdout=f, stderr=f)
            subprocess.run(f"fastqc raw.R1.fastq {Pre}.R1.fastq.gz -o R1_qc -t {threads}", shell=True, stdout=f, stderr=f)
            subprocess.run("unzip -o 'R1_qc/*.zip' ", shell=True, stdout=f, stderr=f)
            subprocess.run(f"fastp -i {Pre}.R1.fastq.gz -o tt.1.fq -w {threads} --json {Pre}.final.json", shell=True, stdout=f, stderr=f)
            subprocess.run("rm tt.*.fq", shell=True)

    summaryfastqc_prod(Pre)
    open("QC_ok", "w").write("已跑过")
    time.sleep(0.5)
