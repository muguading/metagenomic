from __future__ import annotations

import os
import re
import subprocess

import pandas as pd

from metagenomic_refactor.assembly import outfun, run_genovi_summary, write_gene_summaries
from metagenomic_refactor.common import conda_base_bin, conda_run_command
from metagenomic_refactor.context import get_runtime_context


def AnnoFun(Pre, threads):
    # Compatibility shim: the original standalone AnnoFun implementation is no
    # longer present in the workspace, so we preserve the workflow entrypoint
    # and reuse the functional-annotation outputs still generated in the
    # refactored pipeline.
    with open("AnnoFun.log", "a") as f1:
        if os.path.isfile(f"{Pre}_prokka/{Pre}.gbk"):
            run_genovi_summary(Pre)
        if os.path.isfile(f"{Pre}_prokka/{Pre}.tsv") and os.path.isfile(f"{Pre}_prokka/{Pre}.txt"):
            write_gene_summaries(Pre)
        else:
            f1.write(f"{Pre}\tmissing prokka outputs for AnnoFun compatibility mode\n")


def DrugFinder(Pre, threads):
    with open("SV.log", "w") as f:
        if os.path.isfile(f"{Pre}.final.fastq"):
            subprocess.run(f"minimap2 -ax map-ont /data/deploy/TB_soft/ref/TB/ref.fa {Pre}.final.fastq -t {threads}|samtools sort -o 2.CuteSV/{Pre}.ref.sort.bam", shell=True, stdout=f, stderr=f)
            subprocess.run(f"samtools index 2.CuteSV/{Pre}.ref.sort.bam", shell=True, stdout=f, stderr=f)
            subprocess.run(f"cuteSV 2.CuteSV/{Pre}.ref.sort.bam /data/deploy/TB_soft/ref/TB/ref.fa 2.CuteSV/{Pre}.cuteSV.vcf  2.CuteSV/ --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --threads {threads} -s 10", shell=True, stdout=f, stderr=f)
            subprocess.run(f"/home/dell/biosoft/snpEff/scripts/snpEff ann -noLog -noStats -no-downstream -no-upstream -no-utr -c reference/snpeff.config -dataDir . ref 2.CuteSV/{Pre}.cuteSV.vcf > 2.CuteSV/{Pre}.anno.vcf", shell=True, stdout=f, stderr=f)
            subprocess.run(f"python /data/deploy/TB_soft/other_soft/IGV_js/IGV_new.py -r ref.fa -m 2.CuteSV/{Pre}.ref.sort.bam -o 2.CuteSV/ -s {Pre}", shell=True, stdout=f, stderr=f)
            subprocess.run("cp /data/deploy/TB_soft/other_soft/IGV_js/igv.min.js 2.CuteSV/", shell=True)
        if os.path.isfile(f"{Pre}.R1.fastq.gz"):
            if os.path.isfile(f"{Pre}.R2.fastq.gz"):
                subprocess.run(f"bwa mem /data/deploy/TB_soft/ref/TB/ref.fa {Pre}.R1.fastq.gz {Pre}.R2.fastq.gz -t {threads}|samtools sort -o 2.CuteSV/{Pre}.ngs_ref.sort.bam", shell=True, stdout=f, stderr=f)
            else:
                subprocess.run(f"bwa mem /data/deploy/TB_soft/ref/TB/ref.fa {Pre}.R1.fastq.gz -t {threads}|samtools sort -o 2.CuteSV/{Pre}.ngs_ref.sort.bam", shell=True, stdout=f, stderr=f)
            subprocess.run(f"samtools index 2.CuteSV/{Pre}.ngs_ref.sort.bam", shell=True, stdout=f, stderr=f)
            subprocess.run(f"delly call -o 2.CuteSV/{Pre}.delly.bcf -s 3 -q 20 -g ref.fa 2.CuteSV/{Pre}.ngs_ref.sort.bam", shell=True)
            subprocess.run(f"bcftools view 2.CuteSV/{Pre}.delly.bcf > 2.CuteSV/{Pre}.delly.vcf", shell=True)
            subprocess.run(f"""bcftools view -i "FILTER='PASS' & GT='1/1'" 2.CuteSV/{Pre}.delly.vcf > 2.CuteSV/{Pre}.delly.filt.vcf""", shell=True)
            subprocess.run(f"/home/dell/biosoft/snpEff/scripts/snpEff ann -noLog -noStats -no-downstream -no-upstream -no-utr -c reference/snpeff.config -dataDir . ref  2.CuteSV/{Pre}.delly.filt.vcf > 2.CuteSV/{Pre}.delly.anno.vcf", shell=True)
    open("SV_ok", "w").write("")


def getvfID(value):
    matches = re.findall(r"VF\d+", str(value))
    if len(matches) > 0:
        if len(matches) > 1:
            newID = "|".join([i for i in matches])
        else:
            newID = matches[0]
    else:
        newID = "-"
    return newID


def hAMRCom(Pre):
    subprocess.run(f"abricate --db card {Pre}.final.fasta > {Pre}.abricate.tsv", shell=True)
    subprocess.run(f"{conda_run_command('amr_aux', f'hamronize abricate {Pre}.abricate.tsv --format tsv --analysis_software_version 1.0.1 --reference_database_version 20250207')} > {Pre}.hamr.abricate.tsv", shell=True)
    subprocess.run(conda_run_command("amr_aux", f"rgi main -i {Pre}.final.fasta -o {Pre}.rgi --clean --include_loose"), shell=True)
    subprocess.run(f"{conda_run_command('amr_aux', f'hamronize rgi {Pre}.rgi.txt  --format tsv --analysis_software_version 6.0.3 --reference_database_version 20250207 --input_file_name {Pre}.final')} > {Pre}.hamr.rgi.tsv", shell=True)
    subprocess.run(conda_run_command("amr_aux", f"run_resfinder.py -ifa {Pre}.final.fasta -o {Pre}_resfinder -acq"), shell=True)
    subprocess.run(f"{conda_run_command('amr_aux', f'hamronize resfinder {Pre}_resfinder/ResFinder_results_tab.txt  --format tsv --analysis_software_version 4.6.0 --reference_database_version 20250207 --input_file_name {Pre}.final')} > {Pre}.hamr.resfinder.tsv", shell=True)
    subprocess.run(conda_run_command("amr_aux", f"amrfinder -n {Pre}.final.fasta -o {Pre}.amrfinder.tsv -t 10"), shell=True)
    subprocess.run(f"hamronize summarize *.hamr*.tsv -t interactive -o {Pre}.hamr.html", shell=True)


def assem_vfdr(Pre, inty="fastq"):
    runtime = get_runtime_context()
    runtime_krdb = runtime.krdb
    runtime_vfmeta = runtime.vfmeta
    print("组装与耐药毒力分析开始")
    if inty == "fastq":
        if os.path.isfile(f"{Pre}_1.regions.bed"):
            subprocess.run(f"ln -s {Pre}_1.regions.bed {Pre}_assem.regions.bed", shell=True)
        else:
            subprocess.run(f"ln -s {Pre}_ngs_1.regions.bed {Pre}_assem.regions.bed", shell=True)
    if not os.path.isfile(f"{Pre}_assem.kraken2.txt"):
        subprocess.run(f"kraken2 --db {runtime_krdb} --threads 10 --output {Pre}_assem.txt --report {Pre}_assem.kraken2.txt {Pre}.final.fasta", shell=True)
    subprocess.run(f"abricate  --minid 50 --mincov 50 --db card --threads 10 --quiet  {Pre}.final.fasta > Assem_abricate_CARD.txt", shell=True)
    subprocess.run(f"abricate  --minid 50 --mincov 50 --db vfdb --threads 10 --quiet  {Pre}.final.fasta > Assem_abricate_VFDB.txt", shell=True)
    asvfdb = pd.read_table("Assem_abricate_VFDB.txt")
    asdrdb = pd.read_table("Assem_abricate_CARD.txt")
    ask2db = pd.read_table(f"{Pre}_assem.txt", names=["type1", "contig", "taxid", "readsL", "taxidlist"])
    ask2db1 = pd.read_table(f"{Pre}_assem.kraken2.txt", names=["Abundance", "ReadsC", "ReadsO", "taxid", "SciName"])
    ask2db1["Name"] = ask2db1["SciName"].str.strip()
    ask2db = ask2db.merge(ask2db1, on="taxid", how="left")
    asdrdb = asdrdb.merge(ask2db, left_on="SEQUENCE", right_on="contig", how="left")
    asvfdb = asvfdb.merge(ask2db, left_on="SEQUENCE", right_on="contig", how="left").fillna("-")
    asvfdb = asvfdb[["SEQUENCE", "Name", "taxid", "GENE", "%COVERAGE", "%IDENTITY", "PRODUCT", "START", "END", "STRAND"]]
    asvfdb.rename(columns={"SEQUENCE": "Contig名称", "START": "起始碱基", "END": "终止碱基", "STRAND": "正负链", "GENE": "基因名称", "%COVERAGE": "覆盖度%", "%IDENTITY": "一致性%", "PRODUCT": "产物"}, inplace=True)
    if asvfdb.shape[0] > 0:
        asvfdb["VFID"] = asvfdb["产物"].fillna("").astype(str).apply(getvfID)
    else:
        asvfdb["VFID"] = "-"
    asvfdb = asvfdb.merge(runtime_vfmeta, on="VFID", how="left")
    asvfdb.rename(columns={"Name": "物种名称", "VF_Name": "VF名称", "VF_FullName": "VF全称", "Bacteria": "物种来源", "VFcategory": "VF分类", "Characteristics": "特征", "Structure": "结构", "Function": "功能", "Mechanism": "机制", "Reference": "文献来源"}, inplace=True)
    asvfdb = asvfdb[["Contig名称", "物种名称", "taxid", "基因名称", "覆盖度%", "一致性%", "产物", "VF分类", "VF名称", "起始碱基", "终止碱基", "正负链"]]
    asdrdb = asdrdb[["SEQUENCE", "Name", "taxid", "GENE", "%COVERAGE", "%IDENTITY", "PRODUCT", "RESISTANCE", "START", "END", "STRAND"]]
    asdrdb.rename(columns={"SEQUENCE": "Contig名称", "START": "起始碱基", "END": "终止碱基", "STRAND": "正负链", "GENE": "基因名称", "%COVERAGE": "覆盖度%", "%IDENTITY": "一致性%", "PRODUCT": "产物", "RESISTANCE": "耐药药物", "Name": "物种名称"}, inplace=True)
    asvfdb["基因名称"] = (
        asvfdb["基因名称"].astype(str)
        .str.replace("'", "", regex=False)
        .str.replace("/", "_", regex=False)
        .str.replace(" ", "", regex=False)
        .str.replace("(", "_", regex=False)
        .str.replace(")", "_", regex=False)
    )
    asdrdb["基因名称"] = (
        asdrdb["基因名称"].astype(str)
        .str.replace("'", "", regex=False)
        .str.replace("/", "_", regex=False)
        .str.replace(" ", "", regex=False)
        .str.replace("(", "_", regex=False)
        .str.replace(")", "_", regex=False)
    )
    asvfdb.to_csv("Assem_abricate_VFDB.tsv", sep="\t", index=False)
    asdrdb.to_csv("Assem_abricate_CARD.tsv", sep="\t", index=False)
    if not os.path.isdir("geneDepth"):
        os.makedirs("geneDepth")
    if asvfdb.shape[0] > 0:
        asvfdb[["Contig名称", "起始碱基", "终止碱基", "基因名称", "正负链", "物种名称"]].to_csv("Assem_abricate_VFDB.bed", header=False, index=False, sep="\t")
        if inty == "fastq":
            subprocess.run(f"bedtools intersect -a Assem_abricate_VFDB.bed -b {Pre}_assem.regions.bed -wb > Assem_abricate_VFDB.depth.bed", shell=True)
            asvfddb = pd.read_table("Assem_abricate_VFDB.depth.bed", header=None, names=["Chrom", "start", "end", "GeneName", "strand", "Species", "c1", "s1", "e1", "Depth"])
            asvfddb = asvfddb.groupby("GeneName").apply(lambda x: outfun(Pre, x, "vfdb")).reset_index(level=0)
            asvfddb.rename(columns={"GeneName": "基因名称"}, inplace=True)
            asvfddb.to_csv("VFDB_summary.tsv", sep="\t", index=False)
    else:
        open("VFDB_summary.tsv", "w").write("基因名称\t片段名称\t起始位置\t终止位置\t覆盖度(>0)%\t覆盖度(>10)%\t覆盖度(>100)%\t平均深度\t最低深度\t最高深度\n")
    if asdrdb.shape[0] > 0:
        if inty == "fastq":
            asdrdb[["Contig名称", "起始碱基", "终止碱基", "基因名称", "正负链", "物种名称"]].to_csv("Assem_abricate_CARD.bed", header=False, index=False, sep="\t")
            subprocess.run("bedtools intersect -a Assem_abricate_CARD.bed -b {Pre}_assem.regions.bed -wb > Assem_abricate_CARD.depth.bed".format(Pre=Pre), shell=True)
            asdrddb = pd.read_table("Assem_abricate_CARD.depth.bed", header=None, names=["Chrom", "start", "end", "GeneName", "strand", "Species", "c1", "s1", "e1", "Depth"])
            asdrddb = asdrddb.groupby("GeneName").apply(lambda x: outfun(Pre, x, "card")).reset_index(level=0)
            asdrddb.rename(columns={"GeneName": "基因名称"}, inplace=True)
            asdrddb.to_csv("CARD_summary.tsv", sep="\t", index=False)
    else:
        open("CARD_summary.tsv", "w").write("基因名称\t片段名称\t起始位置\t终止位置\t覆盖度(>0)%\t覆盖度(>10)%\t覆盖度(>100)%\t平均深度\t最低深度\t最高深度\n")
    assdb = pd.read_table("flye_output/assembly_info.txt")
    assdb = assdb.merge(ask2db, how="left", left_on="序列名称", right_on="contig")
    assdb = assdb[["序列名称", "序列长度", "平均深度", "是否成环", "基因组/质粒", "质粒分型", "taxid", "Name"]]
    assdb["毒力基因数量"] = assdb.apply(lambda x: asvfdb[asvfdb["Contig名称"] == x["序列名称"]].shape[0], axis=1)
    assdb["毒力基因"] = assdb.apply(lambda x: ",".join(asvfdb.loc[asvfdb["Contig名称"] == x["序列名称"], "基因名称"].tolist()) if asvfdb.loc[asvfdb["Contig名称"] == x["序列名称"]].shape[0] > 0 else "-", axis=1)
    assdb["耐药基因数量"] = assdb.apply(lambda x: asdrdb[asdrdb["Contig名称"] == x["序列名称"]].shape[0], axis=1)
    assdb["耐药基因"] = assdb.apply(lambda x: ",".join(asdrdb.loc[asdrdb["Contig名称"] == x["序列名称"], "基因名称"].tolist()) if asdrdb.loc[asdrdb["Contig名称"] == x["序列名称"]].shape[0] > 0 else "-", axis=1)
    assdb = assdb.rename(columns={"contig": "Contig名称", "length": "片段长度", "Name": "物种名称"})
    assdb.to_csv("Assem_info1.tsv", sep="\t", index=False)
    assdb[["序列名称", "序列长度", "平均深度", "是否成环", "基因组/质粒", "质粒分型", "物种名称", "毒力基因数量", "毒力基因", "耐药基因数量", "耐药基因"]].to_csv("Assem_info.tsv", sep="\t", index=False)


def VFDR(Pre, threads, intype="fastq"):
    runtime = get_runtime_context()
    runtime_method = runtime.method
    runtime_vfmeta = runtime.vfmeta
    with open("vfdr.log", "w") as f1:
        subprocess.run(f"abricate --db vfdb {Pre}.final.fasta --threads {threads} --minid 50 --mincov 50 >{Pre}.vfdb.tsv", shell=True, stdout=f1, stderr=f1)
        vfdb = pd.read_table(f"{Pre}.vfdb.tsv")
        vfdb = vfdb[["SEQUENCE", "START", "END", "STRAND", "GENE", "%COVERAGE", "%IDENTITY", "PRODUCT"]]
        vfdb.rename(columns={"SEQUENCE": "Contig名称", "START": "起始碱基", "END": "终止碱基", "STRAND": "正负链", "GENE": "基因名称", "%COVERAGE": "覆盖度%", "%IDENTITY": "一致性%", "PRODUCT": "产物"}, inplace=True)
        vfdb["VFID"] = vfdb["产物"].fillna("").astype(str).apply(getvfID)
        vfdb = vfdb.merge(runtime_vfmeta, on="VFID", how="left")
        vfdb.rename(columns={"VF_Name": "VF名称", "VF_FullName": "VF全称", "Bacteria": "物种来源", "VFcategory": "VF分类", "Characteristics": "特征", "Structure": "结构", "Function": "功能", "Mechanism": "机制", "Reference": "文献来源"}, inplace=True)
        vfdb.to_csv(f"{Pre}.vfdb.tsv", sep="\t", index=False)
        subprocess.run(f"""awk -v OFS='\t' '{{print $1,$2,$3,$5}}' {Pre}.vfdb.tsv|sed '1d' > vfdb.bed""", shell=True)
        if int(os.popen("cat vfdb.bed|wc -l").read()) > 0:
            vfdbdb = pd.read_table("vfdb.bed", header=None)
            vfdbdb[3] = vfdbdb.apply(lambda x: x[3] if len(x[3]) < 10 else f"{x[3][0:5]}..{x[3][-5:]}", axis=1)
            vfdbdb.to_csv("vfdb.bed", sep="\t", index=False, header=False)
        else:
            subprocess.run("rm vfdb.bed", shell=True)
        subprocess.run(f"abricate --db card {Pre}.final.fasta --threads {threads} --minid 50 --mincov 50 >{Pre}.card.tsv", shell=True, stdout=f1, stderr=f1)
        card = pd.read_table(f"{Pre}.card.tsv")
        card = card[["SEQUENCE", "START", "END", "STRAND", "GENE", "%COVERAGE", "%IDENTITY", "PRODUCT", "RESISTANCE"]]
        card.rename(columns={"SEQUENCE": "Contig名称", "START": "起始碱基", "END": "终止碱基", "STRAND": "正负链", "GENE": "基因名称", "%COVERAGE": "覆盖度%", "%IDENTITY": "一致性%", "PRODUCT": "产物", "RESISTANCE": "耐药药物"}, inplace=True)
        card.to_csv(f"{Pre}.card.tsv", sep="\t", index=False)
        subprocess.run(f"""awk -v OFS='\t' '{{print $1,$2,$3,$5}}' {Pre}.card.tsv|sed '1d' > card.bed""", shell=True)
        if int(os.popen("cat card.bed|wc -l").read()) > 0:
            carddb = pd.read_table("card.bed", header=None)
            carddb[3] = carddb.apply(lambda x: x[3] if len(x[3]) < 10 else f"{x[3][0:5]}..{x[3][-5:]}", axis=1)
            carddb.to_csv("card.bed", sep="\t", index=False, header=False)
        else:
            subprocess.run("rm card.bed", shell=True)
        assem_vfdr(Pre, intype)
        tmpbind = "/".join(os.getcwd().split("/")[:2])
        if runtime_method != "meta":
            subprocess.run(f"""{conda_base_bin("singularity")} exec --bind {tmpbind}:{tmpbind} /home/dell/biosoft/mummer2circos.simg mummer2circos -q {os.getcwd()}/{Pre}.final.fasta -r {os.getcwd()}/{Pre}.final.fasta -gb {os.getcwd()}/{Pre}_prokka/{Pre}.gbk -l  -o '{Pre}_raw' """, shell=True, stdout=f1, stderr=f1)
            for feature in ["card", "vfdb", "rgi"]:
                if os.path.isfile(f"{feature}.bed"):
                    subprocess.run(f"""{conda_base_bin("singularity")} exec --bind {tmpbind}:{tmpbind} /home/dell/biosoft/mummer2circos.simg mummer2circos -q {os.getcwd()}/{Pre}.final.fasta -r {os.getcwd()}/{Pre}.final.fasta -gb {os.getcwd()}/{Pre}_prokka/{Pre}.gbk -l -lf {feature}.bed -o '{Pre}_{feature}' """, shell=True, stdout=f1, stderr=f1)


def AnnoEle(Pre, threads):
    with open("AnnoElog", "w") as f1:
        subprocess.run(f"minced {Pre}_prokka/{Pre}.fna {Pre}_CRISPRs.txt {Pre}_CRISPRs.gff", shell=True, stderr=f1, stdout=f1)
        subprocess.run(f"""grep 'minced:0.4.2' {Pre}_CRISPRs.gff > {Pre}_CRISPRs.tsv""", shell=True)
        CRIdb = pd.read_table(f"{Pre}_CRISPRs.tsv", header=None, names=["序列ID", "软件版本", "片段类型", "开始位置", "终止位置", "得分", "a", "b", "结果注释"])
        CRIdb = CRIdb[["序列ID", "软件版本", "片段类型", "开始位置", "终止位置", "得分", "结果注释"]]
        CRIdb.to_csv(f"{Pre}.CRISPR.tsv", sep="\t", index=False)
    with open("mobileOG.log", "w") as mbf:
        mgemt = pd.read_table("/data/deploy/meta_genome/database/beatrix/mobileOG-db-beatrix-1.6-All.csv", sep=",", low_memory=False)
        subprocess.run(f"diamond blastp -q {Pre}.faa --db /data/deploy/meta_genome/database/beatrix/mobileOG-db  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --out {Pre}.mgeblast.tsv", shell=True, stdout=mbf, stderr=mbf)
        mgedb = pd.read_table(f"{Pre}.mgeblast.tsv", names=["序列名称", "参考基因组名称", "相似性(%)", "长度", "差异数量", "空缺数量", "序列起始", "序列终止", "参考起始", "参考终止", "evalue", "比对得分"])
        mgedb = mgedb[(mgedb["相似性(%)"] > 25) & (mgedb["evalue"] < 1e-5)]
        mgedb["ID"] = mgedb["参考基因组名称"].str.split("|").str[0]
        mgedb["类型"] = mgedb["参考基因组名称"].str.split("|").str[3:].str.join("|")
        mgedb = mgedb.merge(mgemt, left_on="ID", right_on="mobileOG Entry Name")
        mgedb = mgedb[["序列名称", "参考基因组名称", "类型", "Manual Annotation", "Name", "相似性(%)", "长度", "差异数量", "空缺数量", "序列起始", "序列终止", "参考起始", "参考终止", "evalue", "比对得分"]]
        mgedb.rename(columns={"Manual Annotation": "注释", "Name": "元件名称"}, inplace=True)
        mgedb.to_csv(f"{Pre}.medb.tsv", sep="\t", index=False)
