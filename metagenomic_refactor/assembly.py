from __future__ import annotations

import re
import subprocess

import pandas as pd

from metagenomic_refactor.context import get_runtime_context


LEGACY_CALLBACKS = {}


def register_assembly_callbacks(**kwargs):
    LEGACY_CALLBACKS.update(kwargs)


def asb_func(inf, fq1, fq2, threads, Pre, lelID, pts, pst, method, asmt="longasm", ref="noref", gtf="nogtf", tryref=False):
    return LEGACY_CALLBACKS["legacy_asb_func"](inf, fq1, fq2, threads, Pre, lelID, pts, pst, method, asmt, ref, gtf, tryref)


def getvfID(tmpdb):
    matches = re.findall(r"VF\d+", tmpdb["产物"])
    if len(matches) > 0:
        if len(matches) > 1:
            return "|".join([i for i in matches])
        return matches[0]
    return "-"


def rgi_fun(Pre):
    subprocess.run(f"/home/dell/miniconda3/bin/conda run -n RGI_new rgi main -i {Pre}.final.fasta -o {Pre}.rgi --clean --include_loose --low_quality -g PYRODIGAL", shell=True)
    rgidb = pd.read_table(f"{Pre}.rgi.txt")
    rgidb = rgidb[["Contig", "Start", "Stop", "Orientation", "Cut_Off", "Pass_Bitscore", "Best_Hit_ARO", "Best_Identities", "Model_type", "SNPs_in_Best_Hit_ARO", "AMR Gene Family", "Drug Class"]]
    rgidb["序列名称"] = rgidb["Contig"].str.split("_").str[:2].str.join("_")
    rgidb.rename(columns={"Start": "序列起始", "Stop": "序列终止", "Orientation": "正负链", "Cut_Off": "过滤标准", "Pass_Bitscore": "比对得分", "Best_Hit_ARO": "耐药基因", "Best_Identities": "一致性%", "Model_type": "基因数据库", "SNPs_in_Best_Hit_ARO": "基因突变", "AMR Gene Family": "耐药基因家族", "Drug Class": "药物类别"}, inplace=True)
    rgidb.drop("Contig", axis=1, inplace=True)
    rgidb = rgidb[["序列名称", "序列起始", "序列终止", "正负链", "过滤标准", "比对得分", "耐药基因", "一致性%", "基因数据库", "基因突变", "耐药基因家族", "药物类别"]]
    rgidb.sort_values("比对得分", ascending=False, inplace=True)
    rgidb.to_csv(f"{Pre}.rgi.tsv", sep="\t", index=False)
    rgidb.to_csv(f"{Pre}.rgi.bed", sep="\t", index=False, header=False)


def outfun(Pre, x, typeF):
    runtime = get_runtime_context()
    ts = x["start"].min()
    te = x["end"].max()
    ofname = x["GeneName"].tolist()[0]
    Spename = x["Species"].tolist()[0]
    x["start"] = x.reset_index().index + 1
    x["end"] = x.reset_index().index + 2
    x[["Chrom", "start", "end", "Depth"]].to_csv(f"geneDepth/{ofname}_{typeF}.tsv", sep="\t", header=False, index=False)
    open(f"geneDepth/{ofname}_{typeF}.bed", "w").write(f"{x['Chrom'].tolist()[0]}\t{ts}\t{te}\t{Pre}_{ofname}\t{Pre}_{ofname}\t{x['strand'].tolist()[0]}")
    subprocess.run(f"bedtools getfasta -fi {runtime.wkdir}/{Pre}/{Pre}.final.fasta -bed geneDepth/{ofname}_{typeF}.bed -name -s > geneDepth/{ofname}_{typeF}.fasta", shell=True)
    tmpdict = {"片段名称": x["Chrom"].tolist()[0], "物种名称": Spename, "起始位置": x["start"].min(), "终止位置": x["end"].max(), "覆盖度(>0)%": round(x[x["Depth"] > 0].shape[0] / x.shape[0], 4) * 100, "覆盖度(>10)%": round(x[x["Depth"] > 10].shape[0] / x.shape[0], 4) * 100, "覆盖度(>100)%": round(x[x["Depth"] > 100].shape[0] / x.shape[0], 4) * 100, "平均深度": round(x["Depth"].mean(), 2), "最低深度": x["Depth"].min(), "最高深度": x["Depth"].max()}
    return pd.DataFrame(tmpdict, index=[0]).round(2)
