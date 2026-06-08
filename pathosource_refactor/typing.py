from __future__ import annotations

import os
import subprocess

import pandas as pd
import pymysql

from pathosource_refactor.common import checkfile
from pathosource_refactor.context import get_runtime_context


def mlst() -> None:
    runtime = get_runtime_context()
    with open(f"{runtime.run_path}/task.log", "a") as log_handle:
        subprocess.run("mlst *.fasta >mlst.txt", shell=True, stdout=log_handle, stderr=log_handle)
    headerlist = ["FILE", "物种名称", "ST分型"]
    vi_mlst = pd.read_table("mlst.txt", header=None)
    allenum = vi_mlst.shape[1] - 3
    for i in range(1, allenum + 1):
        headerlist.append(f"Allele {i}")
    vi_mlst.columns = headerlist
    vi_mlst["物种名称"] = runtime.species
    vi_mlst["FILE"] = vi_mlst["FILE"].str.replace(".fasta", "", regex=False)
    vi_mlst = vi_mlst.rename(columns={"FILE": "样本名"})
    vi_mlst.to_csv("sample_MLST_summary.txt", index=0, sep="\t")
    vi_mlst.drop("ST分型", axis=1, inplace=True)
    vi_mlst.drop("物种名称", axis=1, inplace=True)
    for i in range(1, vi_mlst.shape[1]):
        vi_mlst.iloc[:, i] = vi_mlst.iloc[:, i].str.extract(r"\((.*?)\)")
    vi_mlst.dropna(inplace=True)
    vi_mlst.to_csv("mlst.vcf", sep="\t", index=0)
    if vi_mlst.shape[0] >= 3:
        try:
            subprocess.run(f"grapetree -p mlst.vcf -n {runtime.threads} > mlst.nwk", shell=True)
        except Exception:
            print("样本无差异，无法绘制mlst最小生成树")


def cgmlst(pre: str, tspe: str, version: str) -> None:
    runtime = get_runtime_context()
    resources = runtime.resources
    speciesdict = {
        "salmonella": "Salmonella enterica",
        "E_coli": "Escherichia coli",
        "Shigella": "Shigella flexneri",
        "Parahemolyticus": "Vibrio parahaemolyticus",
        "cholerae": "Vibrio cholerae",
        "Y_enterocolitica": "Yersinia enterocolitica",
        "Campylobacter": "Campylobacter coli",
        "Brucella": "Brucella spp",
        "Lmono": "Listeria monocytogenes",
        "Kpne": "Klebsiella pneumoniae",
        "HPinf": "Haemophilus influenzae",
        "Suare": "Staphylococcus aureus",
        "Bcere": "Bacillus cereus",
        "Nmen": "Neisseria meningitidis",
    }
    query = "SELECT * FROM gtdb"
    connection = pymysql.connect(
        host=resources.mysql_host,
        user=resources.mysql_user,
        password=resources.mysql_password,
        database=resources.mysql_database,
    )
    metafile = pd.read_sql(query, connection)
    connection.close()
    metafile.fillna("-", inplace=True)
    if "DB" in metafile.columns:
        metafile = metafile[metafile["DB"] != "-"]
    if "GCA" in tspe or "GCF" in tspe:
        spegcf = tspe
        species = metafile.loc[metafile["ncbi_acc"] == spegcf, "s"].tolist()[0]
    else:
        species = speciesdict.get(tspe, tspe)
    tspe = species.replace(" ", "_")
    print(f"cgmlst Spe:{tspe}")
    if tspe != "" and version != "none":
        folder = f"{tspe}_cgmlst"
        if not os.path.isdir(folder):
            os.makedirs(folder)
        with open("Samplelist.txt") as handle:
            for line in handle:
                sam = line.split("\t")[0]
                tfile = line.split("\t")[1].strip()
                ff = checkfile(tfile)
                if ff == "fasta":
                    infa = line.split("\t")[1].strip()
                else:
                    infa = f"{sam}/snps.consensus.fa"
                subprocess.run(f"cp {infa} {folder}/{sam}.fasta", shell=True)
        with open("cgmlstcall.log", "w") as log_handle:
            if version == "cgmlstorg":
                subprocess.run(f"conda run -n chewie chewBBACA.py AlleleCall -i {folder} -g {resources.cgmlst_org_db}/{tspe}_cgmlst --cpu 10 -o cgmlst_result", shell=True, stdout=log_handle, stderr=log_handle)
                scheme_num = len([i for i in os.listdir(f"{resources.cgmlst_org_db}/{tspe}_cgmlst") if i.endswith("fasta")])
            elif version == "selforg":
                subprocess.run(
                    f"conda run -n chewie chewBBACA.py AlleleCall -i {folder} -g {resources.cgmlst_self_org_db}/{tspe}_Scheme/schema_seed/ --gl {resources.cgmlst_self_org_db}/{tspe}_Scheme_cgMLST/cgMLSTschema95.txt --cpu 10 -o cgmlst_result",
                    shell=True,
                    stdout=log_handle,
                    stderr=log_handle,
                )
                scheme_num = int(os.popen(f"cat {resources.cgmlst_self_org_db}/{tspe}_Scheme_cgMLST/cgMLSTschema95.txt|wc -l").read().strip())
            else:
                tmpgl = [i for i in os.listdir(f"{resources.cgmlst_other_org_db}/{tspe}/{version}/cgmlst_Scheme_cgMLST/") if i.endswith(".txt")][0]
                subprocess.run(
                    f"conda run -n chewie chewBBACA.py AlleleCall -i {folder} -g {resources.cgmlst_other_org_db}/{tspe}/{version}/cgmlst_Scheme/schema_seed/ --gl {resources.cgmlst_other_org_db}/{tspe}/{version}/cgmlst_Scheme_cgMLST/{tmpgl} --cpu 10 -o cgmlst_result",
                    shell=True,
                    stdout=log_handle,
                    stderr=log_handle,
                )
                scheme_num = int(os.popen(f"cat {resources.cgmlst_other_org_db}/{tspe}/{version}/cgmlst_Scheme_cgMLST/{tmpgl}|wc -l").read().strip())

        cgsumdb = pd.read_table("cgmlst_result/results_statistics.tsv")
        cgsumdb1 = cgsumdb.copy()
        cgsumdb1["识别到基因数量"] = cgsumdb.apply(lambda x: x[1:-4].sum(), axis=1)
        cgsumdb1["总等位基因数量"] = scheme_num
        cgsumdb1["样本名称"] = cgsumdb["FILE"]
        cgsumdb1["cgST"] = "-"
        cgsumdb1["匹配等位基因数量"] = cgsumdb["EXC"].tolist()[0]
        cgsumdb1["检测率%"] = cgsumdb1.apply(lambda x: round(x["识别到基因数量"] / x["总等位基因数量"], 2), axis=1)
        cgsumdb1["匹配率%"] = cgsumdb1.apply(lambda x: round(x["匹配等位基因数量"] / x["总等位基因数量"], 2), axis=1)
        cgsumdb1 = cgsumdb1[["样本名称", "总等位基因数量", "识别到基因数量", "检测率%", "cgST", "匹配等位基因数量", "匹配率%"]]
        cgsumdb1.to_csv(f"{pre}_cgmlst_summary_result.tsv", index=False, sep="\t")
        subprocess.run("grapetree -p cgmlst_result/results_alleles.tsv -n 10 > cgmlst.nwk", shell=True)
