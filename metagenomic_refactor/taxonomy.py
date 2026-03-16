from __future__ import annotations

import os
import subprocess

import pandas as pd

from metagenomic_refactor.context import get_runtime_context


def compareid1(level1, level2, rawlist=None):
    rawlist = rawlist or ["R", "D", "K", "P", "C", "O", "F", "G", "S"]
    rlevel1 = [i for i in rawlist if i in level1][0]
    rlevel2 = [i for i in rawlist if i in level2][0]
    if rawlist.index(rlevel1) == rawlist.index(rlevel2):
        if level1 > level2:
            return 1
        if level1 == level2:
            return 0
        return -1
    if rawlist.index(rlevel1) > rawlist.index(rlevel2):
        return 1
    return -1


def proc_kra1(kraken, tax, lel="S"):
    tmplist = [tax]
    if [i for i in ["R", "D", "K", "P", "C", "O", "F", "G", "S"] if i in lel]:
        rawlist = ["R", "D", "K", "P", "C", "O", "F", "G", "S"]
    else:
        rawlist = ["S1", "S2", "S3", "S4", "S5", "S6"]
    if tax != 0:
        kradb = pd.read_table(kraken, header=None)
        kradb[4] = kradb[4].astype("str")
        tmpindex = kradb[(kradb[3] == lel) & (kradb[4] == str(tax))].index.tolist()[0] + 1
        if tmpindex <= kradb.shape[0] - 1:
            while compareid1(kradb.iloc[tmpindex, 3], lel, rawlist) == 1 and tmpindex <= kradb.shape[0] - 2:
                tmplist.append(kradb.iloc[tmpindex, 4])
                tmpindex += 1
    return tmplist


def exreadsID1(taxlist, kraresult, fq1, fq2=0):
    Maintax = taxlist[0]
    kraredb = pd.read_csv(kraresult, header=None, usecols=[1, 2], dtype={1: "str", 2: "int32"}, sep="\t")
    tmp2db = kraredb[kraredb[2].isin(taxlist)]
    tmp1db = pd.DataFrame(tmp2db[1].unique())
    tmp2db.to_csv(f"{Maintax}.id.tsv", sep="\t", index=False)
    pd.DataFrame(tmp1db).to_csv(f"{Maintax}_fqID.txt", index=False, header=False)
    subprocess.run(f"head -n 1 {Maintax}_fqID.txt > tt.txt", shell=True)
    subprocess.run(f"cut -d '/' -f1 {Maintax}_fqID.txt|sort -u > {Maintax}.listID.txt", shell=True)
    if os.popen("head -n 1 tt.txt").read().strip().endswith("/1") or os.popen("head -n 1 tt.txt").read().strip().endswith("/2"):
        subprocess.run(f"sed 's/$/\\/1/' {Maintax}.listID.txt > {Maintax}.listID1.txt", shell=True)
        subprocess.run(f"sed 's/$/\\/2/' {Maintax}.listID.txt > {Maintax}.listID2.txt", shell=True)
        subprocess.run(f"seqkit grep -f {Maintax}.listID1.txt {fq1} > {Maintax}.1.fastq", shell=True)
        if fq2:
            subprocess.run(f"seqkit grep -f {Maintax}.listID2.txt {fq2} > {Maintax}.2.fastq", shell=True)
    else:
        subprocess.run(f"seqkit grep -f {Maintax}.listID.txt {fq1} > {Maintax}.1.fastq", shell=True)
        if fq2:
            subprocess.run(f"seqkit grep -f {Maintax}.listID.txt {fq2} > {Maintax}.2.fastq", shell=True)


def getinfo(Pre, threads=10):
    refdict = {"card": "/home/dell/miniconda3/envs/TB_ONT/db/card/sequences", "vfdb": "/data/deploy/meta_genome/database/vfdb.fasta"}
    for db in ["card", "vfdb"]:
        ref = refdict.get(db)
        if not os.path.isfile(f"2.{db}.sorted.bam") or os.popen(f"samtools view 2.{db}.sorted.bam|wc -l").read().strip() == "0":
            if os.path.isfile("2.2.fastq"):
                subprocess.run(f"minimap2 -ax sr {ref} 2.1.fastq 2.2.fastq -t 10 |samtools sort -o 2.{db}.sorted.bam", shell=True)
            else:
                subprocess.run(f"minimap2 -ax sr {ref} 2.1.fastq -t 10 |samtools sort -o 2.{db}.sorted.bam", shell=True)
            subprocess.run(f"samtools index 2.{db}.sorted.bam", shell=True)
            subprocess.run(f"mosdepth -b1 {db} 2.{db}.sorted.bam -n -t {threads}", shell=True)
            subprocess.run(f"gunzip {db}.regions.bed.gz -f", shell=True)
            subprocess.run(f"samtools idxstat 2.{db}.sorted.bam > {db}.stat.tsv", shell=True)
        dbfile = pd.read_table(f"{db}.regions.bed", header=None)
        coninfo = pd.read_table(f"{db}.stat.tsv", header=None, usecols=[0, 2])
        coninfo.columns = [0, "card_subreads"]
        depdb = pd.DataFrame(dbfile.groupby(0).apply(lambda x: round(sum(x[3]) / x.shape[0], 2)).reset_index(name="card_dep"))
        covdb = pd.DataFrame(dbfile.groupby(0).apply(lambda x: round(sum(x[3] > 0) / x.shape[0], 2)).reset_index(name="card_cov"))
        rawdb = depdb.merge(covdb, on=0).merge(coninfo, on=0).sort_values("card_subreads", ascending=False)
        if db == "card":
            rawdb = rawdb.loc[(rawdb["card_cov"] >= 0.1) & (rawdb["card_subreads"] > 10), :]
            metadb = pd.read_table("/data/deploy/meta_genome/database/aro_index.tsv", usecols=["Model Name", "AMR Gene Family", "Drug Class", "Resistance Mechanism"])
            rawdb["Model Name"] = rawdb[0].str.split("~~~").str[1]
            rawdb = rawdb.merge(metadb, on="Model Name")
            rawdb.columns = ["片段名称", "平均深度", "覆盖率", "支持序列数", "Model", "基因家族", "耐药分类", "耐药机制"]
            rawdb["耐药基因"] = rawdb["片段名称"].str.split("~~~").str[1]
            rawdb[["耐药基因", "平均深度", "覆盖率", "支持序列数", "基因家族", "耐药分类", "耐药机制"]].to_csv("2.card.tsv", sep="\t", index=False)
        else:
            rawdb = rawdb.loc[(rawdb["card_cov"] >= 0.01) & (rawdb["card_subreads"] > 10), :].head(50)
            metadb = pd.read_table("/data/deploy/meta_genome/database/VFs_meta.tsv", encoding="Windows-1252", usecols=["VFID", "Bacteria", "Function", "Mechanism"])
            contigdb = pd.read_table("/data/deploy/meta_genome/database/vfdb.contig.tsv")
            rawdb = rawdb.merge(contigdb, left_on=0, right_on="Contig Name")
            rawdb = rawdb.merge(metadb, on="VFID")
            rawdb["毒力基因"] = rawdb[0].str.split("~~~").str[1]
            rawdb.columns = ["片段名称", "平均深度", "覆盖率", "支持序列数", "Contig", "VFID", "菌株", "毒力功能", "毒力机制", "毒力基因"]
            rawdb[["毒力基因", "平均深度", "覆盖率", "支持序列数", "VFID", "菌株", "毒力功能", "毒力机制"]].to_csv("2.vfdb.tsv", sep="\t", index=False)


def getCovDep(Pre, Pre2):
    if not os.path.isfile("2.1.fastq"):
        taxlist1 = proc_kra1(f"{Pre2}.report.txt", "2", "D")
        taxlist1 = [int(i) for i in taxlist1]
        if os.path.isfile(f"{Pre}.R2.fastq.gz") and os.path.getsize(f"{Pre}.R2.fastq.gz") != 0:
            exreadsID1(taxlist1, f"{Pre2}.out.txt", f"{Pre}.R1.fastq.gz", f"{Pre}.R2.fastq.gz")
        else:
            exreadsID1(taxlist1, f"{Pre2}.out.txt", f"{Pre}.R1.fastq.gz", 0)
    getinfo(Pre)


def run_bracken_sub(report_path, prefix, krdb, kkf):
    testbrkdb = pd.read_table(report_path, header=None)
    if "S4" in testbrkdb[3]:
        level = "S3"
    elif "S3" in testbrkdb[3]:
        level = "S2"
    else:
        level = "S1"
    subprocess.run(f"bracken -d {krdb} -o {prefix}_Sub.bracken1.txt -w {prefix}_Sub.bracken2.txt -l {level} -t 10  -i {report_path}", shell=True, stdout=kkf, stderr=kkf)


def kk2(inf, fq1, fq2, threads, Pre):
    krdb = get_runtime_context().krdb
    with open("kk2.log", "w") as kkf:
        if inf:
            if not os.path.isfile(f"{Pre}.list.txt"):
                if not os.path.isfile(f"{Pre}.report.txt"):
                    subprocess.run(f"kraken2 --db {krdb} --threads {threads} --output {Pre}.out.txt --report {Pre}.report.txt {inf}", shell=True, stdout=kkf, stderr=kkf)
                subprocess.run(f"bracken -d {krdb} -o {Pre}.bracken1.txt -w {Pre}.bracken2.txt -l S -t 10  -i {Pre}.report.txt", shell=True, stdout=kkf, stderr=kkf)
                subprocess.run(f"bracken -d {krdb} -o {Pre}.bracken1.txt -w {Pre}.bracken2.txt -l S3 -t 10  -i {Pre}.report.txt", shell=True, stdout=kkf, stderr=kkf)
            tmpfile = pd.read_table(f"{Pre}.bracken1.txt")
            ONTSpe = tmpfile.name.tolist()[0]
            try:
                getCovDep(Pre, Pre)
            except Exception:
                pass

        if fq1 and fq2:
            if not os.path.isfile(f"{Pre}_2.list.txt") and not os.path.isfile(f"{Pre}_2.report.txt"):
                subprocess.run(f"kraken2 --db {krdb} --threads {threads} --output {Pre}_2.out.txt --report {Pre}_2.report.txt {fq1} {fq2}", shell=True, stdout=kkf, stderr=kkf)
                subprocess.run(f"bracken -d {krdb} -o {Pre}_2.bracken1.txt -w {Pre}_2.bracken2.txt -l S -t 10  -i {Pre}_2.report.txt", shell=True, stdout=kkf, stderr=kkf)
                run_bracken_sub(f"{Pre}_2.report.txt", f"{Pre}_2", krdb, kkf)
            tmpfile2 = pd.read_table(f"{Pre}_2.bracken1.txt")
            ngsSpe = tmpfile2.name.tolist()[0]
            getCovDep(Pre, f"{Pre}_2")
        elif fq1:
            if not os.path.isfile(f"{Pre}_2.list.txt") and not os.path.isfile(f"{Pre}_2.report.txt"):
                subprocess.run(f"kraken2 --db {krdb} --threads {threads} --output {Pre}_2.out.txt --report {Pre}_2.report.txt {fq1}", shell=True, stdout=kkf, stderr=kkf)
                subprocess.run(f"bracken -d {krdb} -o {Pre}_2.bracken1.txt -w {Pre}_2.bracken2.txt -l S -t 10  -i {Pre}_2.report.txt", shell=True, stdout=kkf, stderr=kkf)
                run_bracken_sub(f"{Pre}_2.report.txt", f"{Pre}_2", krdb, kkf)
            tmpfile2 = pd.read_table(f"{Pre}_2.bracken1.txt")
            ngsSpe = tmpfile2.name.tolist()[0]
            try:
                getCovDep(Pre, f"{Pre}_2")
            except Exception:
                pass

        if "ngsSpe" in dir() and "ONTSpe" in dir() and ngsSpe != ONTSpe:
            raise Exception("二三代不是同一菌种测序数据")

        tmpfile = tmpfile[["name", "taxonomy_id", "taxonomy_lvl", "new_est_reads", "fraction_total_reads"]] if "ONTSpe" in dir() else tmpfile2[["name", "taxonomy_id", "taxonomy_lvl", "new_est_reads", "fraction_total_reads"]]
        tmpfile.rename(columns={"name": "物种", "taxonomy_id": "taxid", "taxonomy_lvl": "水平", "new_est_reads": "序列数量", "fraction_total_reads": "相对丰度"}, inplace=True)
        tmpfile.to_csv(f"{Pre}.taxonomy_summary.tsv", sep="\t", index=False)
        open("kk2_ok", "w").write("")
