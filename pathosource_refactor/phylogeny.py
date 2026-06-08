from __future__ import annotations

import os
import subprocess
import sys
from itertools import combinations

import pandas as pd
from Bio import AlignIO

from pathosource_refactor.common import checkfile, select_cols
from pathosource_refactor.context import get_runtime_context


def snippy(input_file_path: str, meta_path: str | None, ref: str) -> None:
    runtime = get_runtime_context()
    resources = runtime.resources
    name_list = []
    path_list = []
    tmpdict = {}
    print(meta_path)
    if meta_path is not None:
        meta = pd.read_table(meta_path, dtype="str")
    else:
        meta = pd.DataFrame({})
    input_file = pd.read_table(input_file_path, dtype="str")
    if meta.shape[0] > 0:
        for sam in meta.iloc[:, 0].tolist():
            database_name = meta.loc[meta["样本名"] == sam, "任务名"].tolist()[0]
            if meta.loc[meta["样本名"] == sam, "是否内置"].tolist()[0] == "yes":
                print(runtime.mode)
                if runtime.mode == "G":
                    leftpath = os.popen(f'find {resources.pathogen_db_root}/{runtime.species}/ -name "{sam}*.fna" ').read().strip()
                else:
                    leftpath = os.popen(f'find {resources.pathogen_snippy_fasta_root}/ -name "{sam}*.f*a" ').read().strip()
                tmpdict[sam] = {"left": leftpath, "right": "none"}
            else:
                if runtime.mode == "G":
                    genome_path = f"/data/outputs2/{database_name}/{sam}.fa"
                    if os.path.isfile(genome_path):
                        tmpdict[sam] = {"left": genome_path, "right": "none"}
                else:
                    tmpdict[sam] = {"left": f"{resources.self_db_root}/{database_name}/{sam}.fa", "right": "none"}

    if input_file.shape[0] > 0:
        for sam in input_file.iloc[:, 0].tolist():
            leftpath = input_file.loc[input_file["样本名称"] == sam, "单端数据"].tolist()[0]
            rightpath = input_file.loc[input_file["样本名称"] == sam, "右端数据"].tolist()[0]
            tmpdict[sam] = {"left": leftpath, "right": rightpath}
    samdb = pd.DataFrame(tmpdict).T
    df_out = samdb.apply(select_cols, axis=1, result_type="expand")
    df_out.insert(0, "SampleName", df_out.index)
    print(df_out)
    df_out.to_csv("Samplelist.txt", index=False, header=False, sep="\t")
    with open("Samplelist.txt") as handle:
        for line in handle:
            sam, fq1, fq2 = line.rstrip("\n").split("\t")
            if checkfile(fq1) == "fasta":
                subprocess.run(f"ln -s {fq1} ./{sam}.raw.fasta", shell=True)


def get_info(pre: str, sdb: str):
    dbdir = f"{get_runtime_context().resources.self_db_root}/{sdb}"
    gff = 0
    if os.path.isfile(f"{dbdir}/{pre}.gff3"):
        gff = f"{dbdir}/{pre}.gff3"
    return gff


def convcf(filen: str, pre: str) -> None:
    rows = int(os.popen("cat " + filen + " |grep '^\\#'|wc -l ").read()) - 1
    table = pd.read_table(filen, skiprows=rows)
    table = table.iloc[:, 9:]
    df2 = pd.DataFrame(table.values.T, index=table.columns, columns=table.index)
    df2.replace(1, 2, inplace=True)
    df2.replace(0, 1, inplace=True)
    df2.to_csv(f"{pre}.vcf", sep="\t")


def consensus_fun(inf: str, method: str, ref: str, ifgubbin: int, db: str) -> None:
    runtime = get_runtime_context()
    with open("consensus.log", "w") as log_handle:
        open("snippylist.txt", "w").write("")
        with open(inf, "r") as sample_handle:
            for line in sample_handle:
                open("snippylist.txt", "a").write(f"{line.strip()}\n")
        inf = "snippylist.txt"
        if method == "snippy":
            print(f"使用{method}生成一致性序列")
            if os.path.isfile(ref):
                subprocess.run(f"snippy-multi {inf} --reference {ref} > runme.sh", shell=True, stdout=log_handle, stderr=log_handle)
                subprocess.run("sh runme.sh", shell=True, stdout=log_handle, stderr=log_handle)
                subprocess.run("snippy-clean_full_aln core.full.aln > clean.full.aln", shell=True, stdout=log_handle, stderr=log_handle)
                if ifgubbin:
                    subprocess.run("run_gubbins.py -p gubbins clean.full.aln", shell=True, stdout=log_handle, stderr=log_handle)
                    subprocess.run("snp-sites -c gubbins.filtered_polymorphic_sites.fasta > clean.core.aln", shell=True, stdout=log_handle, stderr=log_handle)
                else:
                    subprocess.run("snp-sites -c clean.full.aln > clean.core.aln", shell=True, stdout=log_handle, stderr=log_handle)
            else:
                raise ValueError("使用 snippy 必须输入参考基因组")
        elif method == "ska2":
            print(f"使用{method}生成一致性序列")
            subprocess.run(f"ska build -f {inf} -o merged.ska", shell=True, stdout=log_handle, stderr=log_handle)
            subprocess.run("ska align -o core.aln merged.ska.skf --threads 10 --ambig-mask", shell=True, stdout=log_handle, stderr=log_handle)
            subprocess.run("ska align -m 0 -o full.aln merged.ska.skf --threads 10", shell=True, stdout=log_handle, stderr=log_handle)
            subprocess.run("snp-sites -c core.aln > clean.core.aln", shell=True)
            subprocess.run("mv full.aln clean.full.aln", shell=True)
        else:
            print(f"使用{method}生成一致性序列")
            if not os.path.isdir("gff"):
                os.makedirs("gff")
            with open(inf) as sample_handle:
                for line in sample_handle:
                    sam, fafile = line.strip().split("\t")
                    gff, skf = [0, 0]
                    if not gff:
                        subprocess.run(
                            f"conda run --no-capture-output -n GTDBtk bakta --threads 10 {fafile} --db {runtime.resources.bakta_db} --skip-trna --skip-tmrna --skip-ncrna --skip-pseudo --skip-rrna --skip-ncrna --skip-ncrna-region --skip-crispr -o gff -p {sam} --force",
                            shell=True,
                            stderr=log_handle,
                            stdout=log_handle,
                        )
                    else:
                        subprocess.run(f"ln -s {runtime.resources.self_db_root}/{db}/{sam}.gff3 gff/{sam}.gff3", shell=True)
            subprocess.run("roary -e -p 10 gff/*.gff3", shell=True, stderr=log_handle, stdout=log_handle)
            subprocess.run("snp-sites -c core_gene_alignment.aln > clean.core.aln", shell=True)
            subprocess.run("mv core_gene_alignment.aln clean.full.aln", shell=True)
        subprocess.run("snp-sites -v clean.core.aln > core.vcf", shell=True)
        convcf("core.vcf", "convert")
        subprocess.run(f"grapetree -p convert.vcf  -n {runtime.threads} > grapetree.nwk", shell=True)

        if not os.path.isfile("core.tab"):
            tmpclist = ["#CHROM", "POS", "REF"]
            coredb = pd.read_table("core.vcf", skiprows=3)
            slist = coredb.columns.tolist()[9 : coredb.shape[1]]
            for tsam in slist:
                coredb[tsam] = coredb.apply(lambda x: x["ALT"] if x[tsam] != 0 else x["REF"], axis=1)
            tmpclist.extend(coredb.columns.tolist()[9 : coredb.shape[1]])
            coredb = coredb.loc[:, tmpclist]
            coredb.rename(columns={"#CHROM": "CHR"}, inplace=True)
            coredb.to_csv("core.tab", sep="\t", index=False)
        tmpfile = pd.read_table("core.tab")
        rndict = {"CHR": "染色体", "POS": "变异位点位置", "REF": "参考位点碱基"}
        tmpfile.rename(columns=rndict, inplace=True)
        tmpfile.to_csv("Mutate.tsv", sep="\t", index=False)


def check_msa_fasta(fasta_file: str) -> str:
    try:
        alignment = AlignIO.read(fasta_file, "fasta")
    except Exception as exc:
        return f"无法读取文件: {exc}"
    length = alignment.get_alignment_length()
    for record in alignment:
        if len(record.seq) != length:
            return "序列长度不一致"
    valid_chars = set("ATCGatcg-N")
    for record in alignment:
        if not set(record.seq).issubset(valid_chars):
            return "发现非法字符"
    return "MSA FASTA格式正常"


def build_tree(inf: str, method: str, subs: str = "MIX+MFP", threads: int = 10, bs: int = 1000) -> None:
    if check_msa_fasta(inf) != "MSA FASTA格式正常":
        raise ValueError("输入的msa文件格式有问题")
    subprocess.run(f"seqkit grep -v -p 'Reference' {inf} > rmref.tmp.aln", shell=True)
    subprocess.run("snp-sites -c rmref.tmp.aln > rmref.core.aln", shell=True)
    inf = "rmref.core.aln"
    with open("buildtree.log", "w") as log_handle:
        if method == "NJ":
            subprocess.run(f"/data1/shanghai_pip/meta_genome/soft/soft/rapidNJ-master/bin/rapidnj {inf} -i fa -b {bs} -c {threads}|sed 1d > clean.core.aln.contree", shell=True)
        elif method == "MP":
            subprocess.run(f"/data1/shanghai_pip/meta_genome/soft/soft/mpboot-avx -s {inf} -B {bs}", shell=True, stderr=log_handle, stdout=log_handle)
            subprocess.run(f"mv {inf}.tree  clean.core.aln.contree", shell=True)
        elif method == "ML":
            subprocess.run(f"iqtree2 -m {subs} -s {inf} -T {threads} -B {bs} -alrt 1000", shell=True, stderr=log_handle, stdout=log_handle)
            subprocess.run(f"mv {inf}.treefile clean.core.aln.contree", shell=True)
        elif method == "Bayes":
            pass
        subprocess.run("sed -i \"s/'//g\" clean.core.aln.contree", shell=True)


def generate_combinations(lst):
    return list(combinations(lst, 2))


def merge_lists_with_common_elements(lst):
    merged_list = []
    for sublist in lst:
        if not any(set(sublist) & set(item) for item in merged_list):
            merged_list.append(sublist)
        else:
            merged = False
            for i, existing in enumerate(merged_list):
                if set(sublist) & set(existing):
                    merged_list[i] = list(set(sublist) | set(existing))
                    merged = True
            if not merged:
                merged_list.append(sublist)
    return merged_list


def Pair_dis(nsnp: int = 10, binmod=0, dism: str = "ANI") -> None:
    nsnp = int(nsnp)
    tmpdb = pd.read_table("Samplelist.txt", header=None, dtype="str")
    tmpdb[[0]].to_csv("Samlist.txt", sep="\t", index=False, header=False)
    samlist = tmpdb[0].tolist()
    subprocess.run("seqkit split -i clean.full.aln -O ./ --by-id-prefix 'split_full' ", shell=True)
    subprocess.run("seqkit split -i clean.core.aln -O ./ --by-id-prefix 'split_core' ", shell=True)
    if os.path.isfile("genomelist.txt"):
        subprocess.run("rm genomelist.txt", shell=True)
    if os.path.isfile("genome_corelist.txt"):
        subprocess.run("rm genome_corelist.txt", shell=True)
    for item in os.listdir():
        if item.startswith("split_full") and item.endswith("aln") and not item.endswith("Reference.aln"):
            open("genomelist.txt", "a").write(f"{os.getcwd()}/{item}\n")
        elif item.startswith("split_core") and item.endswith("aln") and not item.endswith("Reference.aln"):
            open("genome_corelist.txt", "a").write(f"{os.getcwd()}/{item}\n")
    subprocess.run("seqkit grep -v -p Reference clean.core.aln > core_dref.aln", shell=True)
    subprocess.run("seqkit grep -v -p Reference clean.full.aln > core_dref.full.aln", shell=True)
    subprocess.run("snp-dists core_dref.aln > dis.mat.txt", shell=True)
    subprocess.run("snp-dists -m core_dref.aln > t_dis.tsv", shell=True)
    df = pd.read_table("t_dis.tsv", names=["A", "B", "D"])
    df = df[df["A"] != df["B"]]
    df["AB_tuple"] = df[["A", "B"]].apply(lambda x: sorted(list(x)), axis=1)
    df = df[~df["AB_tuple"].duplicated()]
    if not binmod:
        disdict = {"0": 0, "1-10": 0, "10-100": 0, "100-1000": 0, "1000+": 0}
        disdict["0"] = len([i for i in df["D"].tolist() if i == 0])
        disdict["1-10"] = len([i for i in df["D"].tolist() if i > 0 and i <= 10])
        disdict["10-100"] = len([i for i in df["D"].tolist() if i > 10 and i <= 100])
        disdict["100-1000"] = len([i for i in df["D"].tolist() if i > 100 and i <= 1000])
        disdict["1000+"] = len([i for i in df["D"].tolist() if i > 1000])
    else:
        disdict = {}
        binlist = sorted([int(i) for i in str(binmod).split(",")])
        if len(binlist) >= 2:
            for modnum, mod in enumerate(binlist):
                if modnum == 0:
                    disdict[f"<={mod}"] = len([i for i in df["D"].tolist() if i <= mod])
                elif modnum == len(binlist) - 1:
                    oldmod = int(binlist[modnum - 1])
                    disdict[f"{oldmod}-{mod}"] = len([i for i in df["D"].tolist() if i > oldmod and i <= mod])
                    disdict[f"{mod}+"] = len([i for i in df["D"].tolist() if i > mod])
                else:
                    oldmod = int(binlist[modnum - 1])
                    disdict[f"{oldmod}-{mod}"] = len([i for i in df["D"].tolist() if i > oldmod and i <= mod])
        else:
            mod = binlist[0]
            disdict[f"<={mod}"] = len([i for i in df["D"].tolist() if i <= mod])
            disdict[f">{mod}"] = len([i for i in df["D"].tolist() if i > mod])
    pd.DataFrame(disdict, index=["0"]).to_csv("dis_bin.tsv", sep="\t", index=False)
    clusterlist = []
    for tmpi in df.index:
        tmpdf = df.loc[df.index == tmpi, :]
        if tmpdf["D"].tolist()[0] <= nsnp:
            clusterlist.append([tmpdf["A"].tolist()[0], tmpdf["B"].tolist()[0]])
    merlist = merge_lists_with_common_elements(clusterlist)
    clusternum = 1
    if os.path.isfile("Cluster.tsv"):
        subprocess.run("rm Cluster.tsv", shell=True)
    open("Cluster.tsv", "w").write("聚类名称\t样本数量\t聚类样本\t最大snp差异\t最小snp差异\t平均snp差异\n")
    for cluster in merlist:
        snplist = []
        all_clu = generate_combinations(cluster)
        for sam1, sam2 in all_clu:
            tmpdisdf = df.loc[((df["A"] == sam1) & (df["B"] == sam2)) | ((df["A"] == sam2) & (df["B"] == sam1)), :]
            snplist.append(tmpdisdf["D"].tolist()[0])
        open("Cluster.tsv", "a").write(f"Cluster{clusternum}\t{len(cluster)}\t{','.join(cluster)}\t{max(snplist)}\t{min(snplist)}\t{sum(snplist)/len(snplist)}\n")
        clusternum += 1
    if dism in ["TN93", "SNP"]:
        subprocess.run("java -jar ~/miniconda3/envs/PathoSource/bin/SeqRuler.jar -i core_dref.aln -d {0} -o t_Gdis.txt -a average -c 10".format(dism), shell=True)
        misdf = pd.read_table("t_Gdis.txt", sep=",")
        misdf.to_csv("Gdis.txt", sep="\t", index=False)
        rmisdf = misdf[["Target", "Source", "Distance"]]
        rmisdf.columns = ["Source", "Target", "Distance"]
        misdf = pd.concat([misdf, rmisdf]).reset_index(drop=True)
        for tsam in samlist:
            misdf.loc[len(misdf)] = [tsam, tsam, 0]
        misdf["xlab"] = "-"
        misdf["ylab"] = "-"
        for tindex in misdf.index:
            xsam = misdf.loc[misdf.index == tindex, "Source"].tolist()[0]
            ysam = misdf.loc[misdf.index == tindex, "Target"].tolist()[0]
            misdf.loc[misdf.index == tindex, "xlab"] = samlist.index(xsam)
            misdf.loc[misdf.index == tindex, "ylab"] = samlist.index(ysam)
        misdf[["Source", "Target", "Distance", "xlab", "ylab"]].to_csv("Gdis_core.txt", sep="\t", index=False)
    elif dism == "ANI":
        with open("ani.log", "w") as anl:
            if int(os.popen("seqkit stat clean.full.aln -T |cut -f5|tail -n1").read().strip()) > 1000000:
                subprocess.run("skani dist --rl genomelist.txt  --ql genomelist.txt > Full_ANI.txt ", shell=True, stdout=anl, stderr=anl)
            else:
                subprocess.run("skani dist --rl genomelist.txt  --ql genomelist.txt > Full_ANI.txt -m 200 ", shell=True, stdout=anl, stderr=anl)
        misdf = pd.read_table("Full_ANI.txt", dtype="str")
        misdf["Source"] = misdf["Ref_name"]
        misdf["Target"] = misdf["Query_name"]
        misdf["Distance"] = misdf["ANI"]
        misdf = misdf[misdf["Source"] != misdf["Target"]]
        misdf["AB_tuple"] = misdf[["Source", "Target"]].apply(lambda x: sorted(list(x)), axis=1)
        misdf = misdf[~misdf["AB_tuple"].duplicated()]
        misdf = misdf[["Source", "Target", "Distance"]]
        misdf.to_csv("Gdis.txt", sep="\t", index=False)
        rmisdf = misdf[["Target", "Source", "Distance"]]
        rmisdf.columns = ["Source", "Target", "Distance"]
        misdf = pd.concat([misdf, rmisdf]).reset_index(drop=True)
        for tsam in samlist:
            misdf.loc[len(misdf)] = [tsam, tsam, 100]
        misdf["xlab"] = "-"
        misdf["ylab"] = "-"
        for tindex in misdf.index:
            xsam = misdf.loc[misdf.index == tindex, "Source"].tolist()[0]
            ysam = misdf.loc[misdf.index == tindex, "Target"].tolist()[0]
            misdf.loc[misdf.index == tindex, "xlab"] = samlist.index(xsam)
            misdf.loc[misdf.index == tindex, "ylab"] = samlist.index(ysam)
        misdf[["Source", "Target", "Distance", "xlab", "ylab"]].to_csv("Gdis_full.txt", sep="\t", index=False)
