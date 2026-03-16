from __future__ import annotations

import os
import re
import subprocess
import sys
import time

import pandas as pd
from Bio import SeqIO

from metagenomic_refactor.context import get_runtime_context
from metagenomic_refactor.taxonomy import exreadsID1, proc_kra1


LEGACY_CALLBACKS = {}


def register_assembly_callbacks(**kwargs):
    LEGACY_CALLBACKS.update(kwargs)


def asb_func(inf, fq1, fq2, threads, Pre, lelID, pts, pst, method, asmt="longasm", ref="noref", gtf="nogtf", tryref=False):
    runtime = get_runtime_context()
    with open("asb.log", "w") as f:
        outputfa = f"{Pre}.consensus.fasta"
        polifa = f"{Pre}.polish.fasta"
        fq1, fq2, _ = prepare_assembly_inputs(inf, fq1, fq2, threads, Pre, lelID, method)

        if not os.path.isfile(f"assm_{method}_{pts}_ok"):
            methodlist = method.split(",")
            if len(methodlist) == 1:
                run_method = methodlist[0]
                if run_method in ["flye", "canu", "unicycler", "masurca", "spades", "raven", "wtdbg2", "miniasm", "meta"]:
                    denovo_asb(inf, fq1, fq2, threads, Pre, pts, pst, run_method, asmt, f, outputfa)
                else:
                    if ref != "noref":
                        reassm_fun(inf, fq1, fq2, threads, Pre, pts, pst, run_method, asmt, f, outputfa)
                    else:
                        print("有参组装未提供参考基因组")
                        sys.exit()
            elif len(methodlist) == 2:
                denovomethod = methodlist[0]
                reasmmethod = methodlist[1]
                denovo_asb(inf, fq1, fq2, threads, Pre, pts, pst, denovomethod, asmt, f, outputfa)
                if ref != "noref":
                    reassm_fun(inf, fq1, fq2, threads, Pre, pts, pst, reasmmethod, asmt, f, outputfa="noforce")
                else:
                    print("有参组装未提供参考基因组")
                    sys.exit()
    if method != "meta":
        if os.path.isfile(outputfa):
            finalfa = f"{Pre}.final.fasta"
            with open("Anno1.log", "w") as f:
                if os.path.isfile(polifa):
                    renamefa(polifa, finalfa)
                else:
                    renamefa(outputfa, finalfa)
                map_assembly_reads(finalfa, inf, fq1, fq2, threads, Pre, asmt, runtime.long_type, f)
                build_assembly_info(finalfa, Pre, method)
                flyedb, plasmidlist = enhance_plasmid_results(finalfa, Pre)
                if os.path.isfile(finalfa):
                    write_finalfasta_stats(finalfa, Pre)
                annotate_with_prokka(finalfa, Pre, threads)
                export_plasmid_gbk_and_cgview(Pre, plasmidlist)
                run_checkm_and_write_summary(Pre, finalfa, threads, runtime.species, flyedb)
                run_genovi_summary(Pre)
        else:
            raise Exception(f"{method}组装报错可能原因:\n1.数据量过少，提高数据量后继续分析\n2.样本中含有一定量的杂菌污染")
    time.sleep(0.5)


def Annotate_func(Pre, threads):
    runtime = get_runtime_context()
    with open("Anno.log", "w") as f1:
        finalfa = f"{Pre}.final.fasta"
        renamefa(finalfa, finalfa)
        if not os.path.isfile(f"{Pre}_prokka/{Pre}.tsv"):
            if os.path.isfile(finalfa):
                write_finalfasta_stats(finalfa, Pre)
            annotate_with_prokka(finalfa, Pre, threads, "Anno.log")
            flyedb = build_fasta_assembly_info(finalfa, Pre, runtime.method)
            flyedb, plasmidlist = enhance_plasmid_results(finalfa, Pre, "Anno.log")
            export_plasmid_gbk_and_cgview(Pre, plasmidlist)
            run_checkm_and_write_summary(Pre, finalfa, threads, runtime.species, flyedb)
            run_genovi_summary(Pre)
        write_gene_summaries(Pre)
        time.sleep(5)


def polish_func(Pre, ptimes, threads, psoft="medaka"):
    print(f"开始抛光 抛光软件: {psoft} 抛光次数: {ptimes}")
    ptimes = str(ptimes)
    medaka_cmd = "/home/dell/miniconda3/bin/conda run -n medaka medaka_consensus"
    try:
        with open("polish.log", "w") as f:
            for i in range(1, int(ptimes) + 1):
                input_file = f"{Pre}.final.fastq" if i == 1 else f"medaka_output{i-1}/consensus.fasta"
                output_dir = f"medaka_output{i}"
                output_file = f"{output_dir}/consensus.fasta"
                cmd = f"{medaka_cmd} -i {input_file} -d {Pre}.consensus.fasta -o {output_dir} -t {threads} > medaka.log"
                subprocess.run(cmd, shell=True, stdout=f, stderr=f)
                if i == int(ptimes):
                    subprocess.run(f"seqkit seq -w0 {output_file} > {Pre}.polish.fasta", shell=True, stdout=f, stderr=f)
        if os.path.getsize(f"{Pre}.polish.fasta") == 0:
            subprocess.run(f"seqkit seq -w0 {Pre}.consensus.fasta > {Pre}.polish.fasta", shell=True)
    except Exception as e:
        print(f"抛光过程出现错误: {e}")
    print("抛光结束")


def wait_for_file(filepath, cinterval=2):
    while not os.path.exists(filepath):
        time.sleep(1)
    lastsize = 1
    while True:
        current_size = os.path.getsize(filepath)
        if current_size == lastsize:
            break
        lastsize = current_size
        time.sleep(cinterval)
    return True


def rebinning():
    infpath = "BASALT_out/meta_drep_out/dereplicated_genomes/"
    outpath = "BASALT_out/meta_drep_out/binning_genomes/"
    if not os.path.isdir(outpath):
        os.makedirs(outpath)
    falist = [i for i in os.listdir(infpath) if i.endswith(".fa")]
    n = 1
    open("binning_name.tsv", "w").write("oldname\tnewname\n")
    if falist:
        for MAG in falist:
            MAG = f"{infpath}/{MAG}"
            oldname = MAG.split("/")[-1]
            oldname = re.sub(r"\.(fa|fasta|fna)(\.gz)?$", "", oldname)
            outMAG = f"{outpath}/MAG_{n}.fa"
            open("binning_name.tsv", "a").write(f"{oldname}\tMAG_{n}\n")
            open(outMAG, "w").write("")
            with open(MAG) as f:
                m = 1
                for line in f:
                    line = line.strip()
                    if line.startswith(">"):
                        open(outMAG, "a").write(f">{m}\n")
                        m += 1
                    else:
                        open(outMAG, "a").write(f"{line}\n")
            n += 1


def combinebin(refinedir, ofa):
    open(ofa, "w").write("")
    list1 = [i for i in os.listdir(refinedir) if i.endswith("fa")]
    for i in list1:
        filen = f"{refinedir}/{i}"
        newname = i.replace(".fa", "")
        with open(filen) as f:
            for line in f:
                if line.startswith(">"):
                    tmp_contig = line.strip().replace(">", "")
                    open(ofa, "a").write(f">{newname}_{tmp_contig}\n")
                else:
                    open(ofa, "a").write(line)


def bingtdbtk_fun():
    inf = "gtdbtk_out/gtdbtk.bac120.summary.tsv"
    with open("gtdbtk.log", "w") as f:
        if not os.path.isfile(inf):
            subprocess.run("/home/dell/miniconda3/bin/conda run -n gtdbtk --no-capture-output gtdbtk classify_wf --genome_dir BASALT_out/meta_drep_out/dereplicated_genomes --out_dir gtdbtk_out -x .fa --cpus 10 --force", shell=True, stdout=f, stderr=f)


def bincheckm2_fun():
    inf = "bin_checkm2out/quality_report.tsv"
    with open("bincheckm2.log", "w") as f:
        if not os.path.isfile(inf):
            subprocess.run("/home/dell/miniconda3/bin/conda run -n cm210 --no-capture-output checkm2 predict --thread 10 --input BASALT_out/meta_drep_out/binning_genomes/ --output-directory bin_checkm2out -x .fa", shell=True, stdout=f, stderr=f)


def binvfdrdb():
    inf1 = "bin_vfdb.tsv"
    inf2 = "bin_card.tsv"
    inf3 = "binning_rgi_new.txt"
    inf4 = "staramr_result/plasmidfinder.tsv"
    with open("binning.log", "w") as f:
        if not os.path.isfile(inf1):
            subprocess.run("abricate BASALT_out/meta_drep_out/binning_genomes/MAG_*.fa --db vfdb > bin_vfdb.tsv ", shell=True, stdout=f, stderr=f)
        if not os.path.isfile(inf2):
            subprocess.run("abricate BASALT_out/meta_drep_out/binning_genomes/MAG_*.fa --db card > bin_card.tsv ", shell=True, stdout=f, stderr=f)
        if not os.path.isfile(inf3):
            subprocess.run("/home/dell/miniconda3/bin/conda run --no-capture-output  -n RGI_new rgi main -i tmp_combine.fa --clean --include_loose -o binning_rgi_new -n 20 -g PYRODIGAL --low_quality", shell=True)
            subprocess.run("rm -r binning_rgi_new.json", shell=True)
        if not os.path.isfile(inf4):
            subprocess.run("staramr search BASALT_out/meta_drep_out/binning_genomes/*.fa -o staramr_result -n 10", shell=True, stdout=f, stderr=f)


def meta_plasmid(Pre):
    with open("plasmid.log", "w") as f:
        if not os.path.isfile(f"{Pre}_plaspredict.tsv") or not os.path.isfile("staramr_result/plasmidfinder.tsv"):
            subprocess.run(f"/home/dell/miniconda3/bin/conda run --no-capture-output -n plasflow PlasFlow.py --input tmp_combine.fa --output {Pre}_plaspredict.tsv", shell=True, stdout=f, stderr=f)
            subprocess.run("staramr search BASALT_out/meta_drep_out/binning_genomes/*.fa -o staramr_result -n 10", shell=True, stdout=f, stderr=f)
        plasmiddb = pd.read_table("staramr_result/plasmidfinder.tsv")
        rawindexname = plasmiddb.index.name
        if plasmiddb.shape[0] > 1:
            plasmiddb = plasmiddb.groupby("Contig").apply(_join_plasmid)
            plasmiddb.index.name = rawindexname
        plasmiddb["contig_name"] = plasmiddb.apply(lambda x: f"{x['Isolate ID']}_{x['Contig']}", axis=1)
        plasmiddb = plasmiddb[["contig_name", "Plasmid"]]
        plasflowdb = pd.read_table(f"{Pre}_plaspredict.tsv")
        plasflowdb = plasflowdb[["contig_name", "label"]]
        plasdb = plasflowdb.merge(plasmiddb, on="contig_name", how="left").fillna("-")
        plasdb.to_csv(f"{Pre}_meta_plaspredict.tsv", sep="\t", index=False)


def meta_tpm():
    with open("coverm.log", "w") as f:
        if not os.path.isfile("meta_tpm.tsv"):
            subprocess.run("/home/dell/miniconda3/bin/conda run -n coverm coverm genome -1 2.1.fastq -2 2.2.fastq -d  BASALT_out/meta_drep_out/binning_genomes/ -x .fa --min-read-percent-identity 95 --min-read-aligned-percent 75 -m tpm -o meta_tpm.tsv -t 10", shell=True, stderr=f, stdout=f)


def binning_result(Pre):
    rebinning()
    combinebin("BASALT_out/meta_drep_out/binning_genomes", "tmp_combine.fa")
    bingtdbtk_fun()
    bincheckm2_fun()
    binvfdrdb()
    meta_plasmid(Pre)
    meta_tpm()
    vfdb = pd.read_table("bin_vfdb.tsv")
    argdict = {"ARG": [], "contig_name": [], "Name": [], "AR Gene(abricate)": [], "AR Gene(rgi)": [], "AR Gene(resfinder)": []}
    vfdb["Name"] = vfdb["#FILE"].str.split("/").str[-1].str.split(".").str[0]
    vfdb["contig_name"] = vfdb.apply(lambda x: f"{x['Name']}_{x['SEQUENCE']}", axis=1)
    vfdb = vfdb[["contig_name", "GENE"]]
    vfdb.rename(columns={"GENE": "VF Gene"}, inplace=True)
    carddb = pd.read_table("bin_card.tsv")
    carddb["tmpgene"] = carddb["GENE"].str.lower()
    rgidb = pd.read_table("binning_rgi_new.txt")
    rgidb = rgidb[rgidb["Cut_Off"].isin(["Strict", "Perfect"])]
    rgidb["tmpgene"] = rgidb["Best_Hit_ARO"].str.lower()
    rgidb["contig_name"] = rgidb["Contig"]
    resdb = pd.read_table("staramr_result/resfinder.tsv")
    resdb["tmpgene"] = resdb["Gene"].str.lower()
    resdb["contig_name"] = resdb.apply(lambda x: f"{x['Isolate ID']}_{x['Contig']}", axis=1)
    carddb["Name"] = carddb["#FILE"].str.split("/").str[-1].str.split(".").str[0]
    carddb["contig_name"] = carddb.apply(lambda x: f"{x['Name']}_{x['SEQUENCE']}", axis=1)
    carddb = carddb[["contig_name", "GENE"]]
    carddb["tmpgene"] = carddb["GENE"].str.lower()
    carddb.rename(columns={"GENE": "AR Gene"}, inplace=True)
    arglist = list(set(resdb["tmpgene"].tolist() + rgidb["tmpgene"].tolist() + carddb["tmpgene"].tolist()))
    for argene in arglist:
        argdict["ARG"].append(argene)
        contiglist = []
        if argene in carddb["tmpgene"].tolist():
            argdict["AR Gene(abricate)"].append("+")
            contiglist.append(carddb.loc[carddb["tmpgene"] == argene]["contig_name"].tolist()[0])
        else:
            argdict["AR Gene(abricate)"].append("-")
        if argene in resdb["tmpgene"].tolist():
            argdict["AR Gene(resfinder)"].append("+")
            contiglist.append(resdb.loc[resdb["tmpgene"] == argene]["contig_name"].tolist()[0])
        else:
            argdict["AR Gene(resfinder)"].append("-")
        if argene in rgidb["tmpgene"].tolist():
            argdict["AR Gene(rgi)"].append("+")
            contiglist.append(rgidb.loc[rgidb["tmpgene"] == argene]["contig_name"].tolist()[0])
        else:
            argdict["AR Gene(rgi)"].append("-")
        if contiglist:
            tmpName = contiglist[0]
            argdict["contig_name"].append(tmpName)
            argdict["Name"].append("_".join(tmpName.split("_")[:2]))
        else:
            argdict["contig_name"].append("-")
    argdb = pd.DataFrame(argdict)
    Alldb = pd.read_table(f"{Pre}_meta_plaspredict.tsv")
    Alldb = Alldb.merge(vfdb, on="contig_name", how="left").merge(argdb, on="contig_name", how="left").fillna("-")
    Alldb["Name"] = Alldb["contig_name"].str.split("_").str[:2].str.join("_")
    gtdbdb = pd.read_table("gtdbtk_out/gtdbtk.bac120.summary.tsv")
    gtdbdb["Name"] = gtdbdb["user_genome"]
    lvlist = ["D", "P", "C", "O", "F", "G", "S"]
    for lv in lvlist:
        gtdbdb[lv] = gtdbdb["classification"].str.split(";").str[lvlist.index(lv)].str.split("__").str[1]
    gtdbdb = gtdbdb[["Name", "D", "P", "C", "O", "F", "G", "S"]]
    binnamedb = pd.read_table("binning_name.tsv")
    gtdbdb = gtdbdb.merge(binnamedb, left_on="Name", right_on="oldname")[["newname", "D", "P", "C", "O", "F", "G", "S"]].rename(columns={"newname": "Name"})
    tpmdb = pd.read_table("meta_tpm.tsv")
    tpmdb.columns = ["Name", Pre]
    Alldb = Alldb.merge(gtdbdb, on="Name", how="left").merge(tpmdb, on="Name", how="left")
    Alldb.to_csv("meta_plas_vf_card.tsv", sep="\t", index=False)


def denovo_asb(inf, fq1, fq2, threads, Pre, pts, pst, method, asmt, f, outputfa):
    runtime = get_runtime_context()
    assembly_long_type = runtime.long_type
    assembly_genome_len = runtime.genome_len
    if asmt == "longasm":
        if method == "flye":
            if assembly_long_type == "Nanopore":
                readsq = float(os.popen(f"nanoq -i {inf} -s -t 5 -vvv|grep 'Mean read quality'|cut -d ':' -f2").read().strip())
                per = 100 - 10 ** (float(readsq) / -10) * 100
                if per > 95:
                    subprocess.run(f"flye --nano-hq {inf} -g {assembly_genome_len} -o flye_output -t {threads} -i 3", shell=True, stdout=f, stderr=f)
                else:
                    subprocess.run(f"flye --nano-raw {inf} -g {assembly_genome_len} -o flye_output -t {threads} -i 3", shell=True, stdout=f, stderr=f)
            elif assembly_long_type == "PacBio_CLR":
                subprocess.run(f"flye --pacbio-raw {inf} -g {assembly_genome_len} -o flye_output -t {threads} -i 3", shell=True, stdout=f, stderr=f)
            elif assembly_long_type == "PacBio_CCS":
                subprocess.run(f"flye --pacbio-hifi {inf} -g {assembly_genome_len} -o flye_output -t {threads} -i 3", shell=True, stdout=f, stderr=f)
            subprocess.run(f"cp flye_output/assembly.fasta {outputfa}", shell=True)
            flyedb = pd.read_table("flye_output/assembly_info.txt", names=["seq_name", "length", "cov", "circ", "repeat", "mult", "alt_group", "graph_path"], skiprows=1)
            flyedb.rename(columns={"seq_name": "序列名称", "length": "序列长度", "cov": "平均深度", "circ": "是否成环"}, inplace=True)
            flyedb[["序列名称", "序列长度", "平均深度", "是否成环"]].to_csv("flye_output/assembly_info.txt", sep="\t", index=False)
        elif method == "miniasm":
            if assembly_long_type == "Nanopore":
                subprocess.run(f"minimap2 -x ava-ont -t {threads} {inf} {inf} | gzip -1 > {Pre}_reads.paf.gz", shell=True, stdout=f, stderr=f)
            elif assembly_long_type == "PacBio_CLR":
                subprocess.run(f"minimap2 -x map-pb -t {threads} {inf} {inf} | gzip -1 > {Pre}_reads.paf.gz", shell=True, stdout=f, stderr=f)
            elif assembly_long_type == "PacBio_CCS":
                subprocess.run(f"minimap2 -x map-hifi -t {threads} {inf} {inf} | gzip -1 > {Pre}_reads.paf.gz", shell=True, stdout=f, stderr=f)
            subprocess.run(f"miniasm -f {inf} {Pre}_reads.paf.gz > {Pre}_reads.gfa", shell=True, stdout=f, stderr=f)
            subprocess.run(f"gfatools gfa2fa {Pre}_reads.gfa > {Pre}_miniasm.fa", shell=True, stdout=f, stderr=f)
            subprocess.run(f"cp {Pre}_miniasm.fa {outputfa}", shell=True)
        elif method == "wtdbg2":
            subprocess.run(f"perl ~/biosoft/wtdbg2/wtdbg2.pl -t {threads} -x ont -g {assembly_genome_len} -o wtdbg2 {inf}", shell=True, stdout=f, stderr=f)
            subprocess.run(f"cp wtdbg2.cns.fa {outputfa}", shell=True)
        elif method == "canu":
            if assembly_long_type == "Nanopore":
                subprocess.run(f"time canu -d canu -p canu genomeSize={assembly_genome_len} maxThreads={threads} -nanopore-raw {inf} >canu.log", shell=True, stdout=f, stderr=f)
            else:
                subprocess.run(f"time canu -d canu -p canu genomeSize={assembly_genome_len} maxThreads={threads} -pacbio-raw {inf} >canu.log", shell=True, stdout=f, stderr=f)
            subprocess.run(f"cp canu/canu.contigs.fasta {outputfa}", shell=True)
        elif method == "unicycler":
            subprocess.run(f"unicycler -t {threads} -l {inf} -o unicycler", shell=True, stdout=f, stderr=f)
            subprocess.run(f"cp unicycler/assembly.fasta {outputfa}", shell=True)
        elif method == "raven":
            subprocess.run(f"raven {inf} -t {threads} > raven.fasta", shell=True, stdout=f, stderr=f)
            subprocess.run(f"cp raven.fasta {outputfa}", shell=True)
        else:
            print(f"请确认传入参数method是否正确，可选[flye,canu,wtdbg2,unicycler,miniasm],您输入的为：{method}")
        if os.path.isfile(outputfa):
            polish_func(Pre, pts, threads, pst)
    elif asmt == "shortasm":
        if method == "spades":
            if fq2:
                subprocess.run(f"spades.py --pe1-1 {fq1} --pe1-2 {fq2} -t {threads} -o spades_output --cov-cutoff 8 --isolate", shell=True, stdout=f, stderr=f)
            else:
                subprocess.run(f"spades.py -s {fq1} -t {threads} -o spades_output --isolate --cov-cutoff 8", shell=True, stdout=f, stderr=f)
            if os.path.isfile("spades_output/contigs.fasta"):
                subprocess.run(f"cp spades_output/contigs.fasta {outputfa}", shell=True, stdout=f, stderr=f)
        elif method == "masurca":
            if fq2:
                subprocess.run(f"masurca -i {fq1},{fq2} -t {threads} -o masurca_output", shell=True, stdout=f, stderr=f)
            else:
                subprocess.run(f"masurca -i {fq1} -t {threads} -o masurca_output", shell=True, stdout=f, stderr=f)
            CAdir = os.popen("ls -d CA").read().split("\n")[0]
            if os.path.isfile(f"{CAdir}/primary.genome.scf.fasta"):
                subprocess.run(f"cp {CAdir}/primary.genome.scf.fasta {outputfa}", shell=True)
            else:
                if os.path.isfile(f"{CAdir}/scaffolds.ref.fa"):
                    subprocess.run(f"cp {CAdir}/scaffolds.ref.fa {outputfa}", shell=True)
                else:
                    print("组装失败")
                    sys.exit()
        elif method == "meta":
            if not os.path.isfile("megahit_output/final.contigs.fa"):
                if fq2:
                    subprocess.run(f"/home/dell/miniconda3/bin/conda run -n BASALT megahit -1 {fq1} -2 {fq2} -o megahit_output", shell=True, stdout=f, stderr=f)
                else:
                    subprocess.run(f"/home/dell/miniconda3/bin/conda run -n BASALT megahit -1 {fq1} -o megahit_output", shell=True, stdout=f, stderr=f)
            if os.path.isfile("megahit_output/final.contigs.fa") and os.path.getsize("megahit_output/final.contigs.fa") != 0:
                tmpwkdir = os.getcwd()
                if not os.path.isdir("BASALT_out"):
                    os.makedirs("BASALT_out")
                subprocess.run(f"ln -s {tmpwkdir}/{fq1} ./BASALT_out", shell=True)
                if fq2:
                    subprocess.run(f"ln -s {tmpwkdir}/{fq2} ./BASALT_out", shell=True)
                subprocess.run(f"ln -s {tmpwkdir}/megahit_output/final.contigs.fa ./BASALT_out", shell=True)
                os.chdir("BASALT_out")
                if os.path.isdir("BestBinset_outlier_refined"):
                    if not [i for i in os.listdir("BestBinset_outlier_refined")]:
                        if fq2:
                            subprocess.run(f"/home/dell/miniconda3/bin/conda run -n BASALT  BASALT -a final.contigs.fa -s {fq1},{fq2} -t {threads} --module autobinning -q checkm2", shell=True, stdout=f, stderr=f)
                            subprocess.run(f"/home/dell/miniconda3/bin/conda run -n BASALT  BASALT -a final.contigs.fa -s {fq1},{fq2} -t {threads} --module refinement -q checkm2", shell=True, stdout=f, stderr=f)
                        else:
                            subprocess.run(f"/home/dell/miniconda3/bin/conda run -n BASALT  BASALT -a final.contigs.fa -s {fq1} -t {threads} --module autobinning -q checkm2", shell=True, stdout=f, stderr=f)
                            subprocess.run(f"/home/dell/miniconda3/bin/conda run -n BASALT  BASALT -a final.contigs.fa -s {fq1} -t {threads} --module refinement -q checkm2", shell=True, stdout=f, stderr=f)
                else:
                    if fq2:
                        subprocess.run(f"/home/dell/miniconda3/bin/conda run -n BASALT  BASALT -a final.contigs.fa -s {fq1},{fq2} -t {threads} --module autobinning -q checkm2", shell=True, stdout=f, stderr=f)
                        subprocess.run(f"/home/dell/miniconda3/bin/conda run -n BASALT  BASALT -a final.contigs.fa -s {fq1},{fq2} -t {threads} --module refinement -q checkm2", shell=True, stdout=f, stderr=f)
                    else:
                        subprocess.run(f"/home/dell/miniconda3/bin/conda run -n BASALT  BASALT -a final.contigs.fa -s {fq1} -t {threads} --module autobinning -q checkm2", shell=True, stdout=f, stderr=f)
                        subprocess.run(f"/home/dell/miniconda3/bin/conda run -n BASALT  BASALT -a final.contigs.fa -s {fq1} -t {threads} --module refinement -q checkm2", shell=True, stdout=f, stderr=f)
                if os.path.isdir("BestBinset_outlier_refined"):
                    binlist = [i for i in os.listdir("BestBinset_outlier_refined") if i.endswith("fa") or i.endswith("fasta")]
                    if len(binlist) > 1:
                        pd.read_table("BestBinset_outlier_refined/quality_report.tsv", header=None)
                        subprocess.run(f"/home/dell/miniconda3/bin/conda run -n BASALT dRep dereplicate meta_drep_out -g BestBinset_outlier_refined/*.fa* --ignoreGenomeQuality -p {threads}", shell=True, stdout=f, stderr=f)
                        if os.path.isdir("meta_drep_out/dereplicated_genomes"):
                            dreplist = [i for i in os.listdir("meta_drep_out/dereplicated_genomes") if i.endswith("fasta") or i.endswith("fa")]
                            if not dreplist:
                                subprocess.run("cp BestBinset_outlier_refined/*.fa* meta_drep_out/dereplicated_genomes/", shell=True)
                        else:
                            os.makedirs("meta_drep_out/dereplicated_genomes")
                            subprocess.run("cp BestBinset_outlier_refined/*.fa* meta_drep_out/dereplicated_genomes/", shell=True)
                    else:
                        if not os.path.isdir("meta_drep_out/dereplicated_genomes"):
                            os.makedirs("meta_drep_out/dereplicated_genomes")
                        subprocess.run("cp BestBinset_outlier_refined/*.fa* meta_drep_out/dereplicated_genomes/tmp.fa", shell=True)
                else:
                    print("refine binning derep失败")
                    subprocess.run(f"cp {tmpwkdir}/megahit_output/final.contigs.fa {tmpwkdir}/{outputfa}", shell=True)
                os.chdir(tmpwkdir)
                binning_result(Pre)
            else:
                print("宏基因组组装失败")
                sys.exit()
    elif asmt == "shortlongasm":
        if method == "unicycler":
            if fq2:
                subprocess.run(f"unicycler -1 {fq1} -2 {fq2} -l {inf} -t {threads} -o unicycler", shell=True, stdout=f, stderr=f)
            elif fq1:
                subprocess.run(f"unicycler -1 {fq1} -l {inf} -t {threads} -o unicycler", shell=True, stdout=f, stderr=f)
            if os.path.isfile("unicycler/assembly.fasta"):
                subprocess.run(f"cp unicycler/assembly.fasta {outputfa}", shell=True)
        elif method == "masurca":
            if fq2:
                subprocess.run(f"masurca -i {fq1},{fq2} -t {threads} -o masurca_output -r {inf}", shell=True, stdout=f, stderr=f)
                CAdir = os.popen("ls -d CA.*").read().split("\n")[0]
                if os.path.isfile(f"{CAdir}/primary.genome.scf.fasta"):
                    subprocess.run(f"cp {CAdir}/primary.genome.scf.fasta {outputfa}", shell=True)
                elif os.path.isfile(f"{CAdir}/scaffolds.ref.fa"):
                    subprocess.run(f"cp {CAdir}/scaffolds.ref.fa {outputfa}", shell=True)
            else:
                subprocess.run(f"masurca -i {fq1} -t {threads} -o masurca_output -r {inf}", shell=True, stdout=f, stderr=f)
                CAdir = os.popen("ls -d CA.*").read().split("\n")[0]
                if os.path.isfile(f"{CAdir}/primary.genome.scf.fasta"):
                    subprocess.run(f"cp {CAdir}/primary.genome.scf.fasta {outputfa}", shell=True)
                elif os.path.isfile(f"{CAdir}/scaffolds.ref.fa"):
                    subprocess.run(f"cp {CAdir}/scaffolds.ref.fa {outputfa}", shell=True)


def reassm_fun(inf, fq1, fq2, threads, Pre, pts, pst, method, asmt, f, outputfa):
    runtime = get_runtime_context()
    runtime_ref = runtime.ref
    runtime_gtf = runtime.gtf
    if not os.path.isdir("genomes"):
        os.makedirs("genomes")
    subprocess.run(f"seqkit seq {runtime_ref} > genomes/ref.fa", shell=True)
    subprocess.run("samtools faidx genomes/ref.fa ", shell=True)
    if runtime_gtf == "nogtf":
        ifAvcf = 0
        print("snp位点不进行额外注释")
        sys.stdout.flush()
    else:
        if not os.path.isdir("ref"):
            os.makedirs("ref")
        subprocess.run(f"cp {runtime_gtf} ref/genes.gff", shell=True)
        open("snpEff.config", "w").write("ref.genome : ref")
        subprocess.run("java -jar /home/dell/biosoft/snpEff/snpEff.jar build -gff3 ref -c snpEff.config -dataDir ./", shell=True)
        if not os.path.isfile("ref/snpEffectPredictor.bin"):
            print("gff与fa不匹配，snp位点不进行额外注释")
            ifAvcf = 0
            sys.stdout.flush()
        else:
            print("注释文件正常，snp位点根据注释文件进行注释")
            ifAvcf = 1
            sys.stdout.flush()
    if "long" in asmt:
        subprocess.run(f"minimap2 -ax map-ont genomes/ref.fa {inf} -t {threads} |samtools sort -o ref.mapping.bam", shell=True)
    elif "short" in asmt:
        with open("mapping.log", "w") as mapplog:
            subprocess.run("bwa index genomes/ref.fa", shell=True)
            if fq2:
                subprocess.run(f"bwa mem genomes/ref.fa {fq1} {fq2} -t {threads} |samtools sort -o ref.mapping.bam", shell=True, stdout=mapplog, stderr=mapplog)
            else:
                subprocess.run(f"bwa mem genomes/ref.fa {fq1}  -t {threads} |samtools sort -o ref.mapping.bam", shell=True, stdout=mapplog, stderr=mapplog)
    subprocess.run("samtools index ref.mapping.bam", shell=True)
    subprocess.run(f"mosdepth -b 1000 ref_map ref.mapping.bam -t {threads}", shell=True)
    subprocess.run("gunzip ref_map.regions.bed.gz", shell=True)
    subprocess.run(f"mosdepth -b 1 ref_map1 ref.mapping.bam -t {threads}", shell=True)
    subprocess.run("gunzip ref_map1.regions.bed.gz", shell=True)

    def _outfun(x):
        tmpdict = {}
        ofname = x["GeneName"].tolist()[0]
        x["start"] = x.reset_index().index + 1
        x["end"] = x.reset_index().index + 2
        x[["Chrom", "start", "end", "Depth"]].to_csv(f"geneDepth/{ofname}.tsv", sep="\t", header=False, index=False)
        tmpdict["片段名称"] = x["Chrom"].tolist()[0]
        tmpdict["起始位置"] = x["start"].min()
        tmpdict["终止位置"] = x["start"].max()
        tmpdict["覆盖度(>0)%"] = round(x[x["Depth"] > 0].shape[0] / x.shape[0], 4) * 100
        tmpdict["覆盖度(>10)%"] = round(x[x["Depth"] > 10].shape[0] / x.shape[0], 4) * 100
        tmpdict["覆盖度(>100)%"] = round(x[x["Depth"] > 100].shape[0] / x.shape[0], 4) * 100
        tmpdict["平均深度"] = round(x["Depth"].mean(), 2)
        tmpdict["最低深度"] = x["Depth"].min()
        tmpdict["最高深度"] = x["Depth"].max()
        return pd.DataFrame(tmpdict, index=[0]).round(2)

    if os.path.isfile("ref/genes.gff"):
        open("geneNamelist.txt", "w").write("")
        if not os.path.isdir("geneDepth"):
            os.makedirs("geneDepth")
        with open("ref/genes.gff") as f1:
            for line in f1:
                if not line.startswith("#"):
                    line = line.strip().split("\t")
                    if line[2] == "gene":
                        if "gene=" in line[8]:
                            gName = line[8].split("gene=")[1].split(";")[0].split("/")[0]
                            open(f"geneDepth/{gName}.bed", "w").write(f"{line[0]}\t{line[3]}\t{line[4]}\t{gName}\n")
                            open("geneNamelist.txt", "a").write(f"{gName}\n")
                        else:
                            if "ID=" in line[8]:
                                gName = line[8].split("ID=")[1].split(";")[0]
                                open(f"geneDepth/{gName}.bed", "w").write(f"{line[0]}\t{line[3]}\t{line[4]}\t{gName}\n")
                                open("geneNamelist.txt", "a").write(f"{gName}\n")
        subprocess.run("cat geneDepth/*.bed > All_gene.bed", shell=True)
        subprocess.run("bedtools intersect -a All_gene.bed -b ref_map1.regions.bed -wb > ref_map1.Anno.regions.bed", shell=True)
        tmpd = pd.read_table("ref_map1.Anno.regions.bed", header=None, names=["Chrom", "start", "end", "GeneName", "c1", "s1", "e1", "Depth"])
        if tmpd.shape[0] != 0:
            newtd = tmpd.groupby("GeneName").apply(_outfun).reset_index(level=0)
            newtd.rename(columns={"GeneName": "基因名称"}, inplace=True)
            newtd.to_csv("gene_summary.tsv", sep="\t", index=False)
        else:
            print("gff与基因组不匹配，无法展示各个基因区段覆盖度")
            sys.stdout.flush()
    subprocess.run(f"python /data1/shanghai_pip/meta_genome/soft/IGV_js/IGV_new.py -r genomes/ref.fa -m ref.mapping.bam -o ./ -s {Pre}", shell=True)
    subprocess.run("cp /data1/shanghai_pip/meta_genome/soft/IGV_js/igv.min.js ./", shell=True)
    subprocess.run("mosdepth -b 1 -t 10 ref ref.mapping.bam", shell=True)
    subprocess.run("gunzip ref.regions.bed.gz", shell=True)
    refbeddb = pd.read_table("ref.regions.bed", header=None)
    refbeddb[refbeddb[3] == 0].to_csv("mask.bed", index=False, header=False, sep="\t")
    if method == "freebayes":
        subprocess.run("fasta_generate_regions.py genomes/ref.fa.fai 200000 > ref.txt", shell=True)
        subprocess.run(f"freebayes-parallel ref.txt {threads} -p 2 -P 0 -C 2 -F 0.05 --min-coverage 10 --min-repeat-entropy 1.0 -q 30 -m 30 --strict-vcf -f genomes/ref.fa ref.mapping.bam > snps.raw.vcf", shell=True)
        subprocess.run("bcftools view --include 'QUAL>=20 && FMT/DP>=10 && (FMT/AO)/(FMT/DP)>=0.9' snps.raw.vcf  | bcftools annotate --remove '^INFO/TYPE,^INFO/DP,^INFO/RO,^INFO/AO,^INFO/AB,^FORMAT/GT,^FORMAT/DP,^FORMAT/RO,^FORMAT/AO,^FORMAT/QR,^FORMAT/QA,^FORMAT/GL' > snps.filt1.vcf", shell=True)
        if ifAvcf:
            subprocess.run("java -jar /home/dell/biosoft/snpEff/snpEff.jar ann  -c snpEff.config -Datadir . ref snps.filt1.vcf > snps.anno.vcf", shell=True)
            skinums = int(os.popen("grep '##' snps.anno.vcf |wc -l").read())
            annodb = pd.read_table("snps.anno.vcf", skiprows=skinums)
            annodb["参考碱基"] = annodb["REF"] + "(" + annodb["unknown"].str.split("|").str[-1].str.strip().str.split(":").str[-4] + ")"
            annodb["突变碱基"] = annodb["ALT"] + "(" + annodb["unknown"].str.split("|").str[-1].str.strip().str.split(":").str[-3] + ")"
            annodb["突变类型"] = annodb["INFO"].str.split("|").str[1]
            annodb["突变影响"] = annodb["INFO"].str.split("|").str[2]
            annodb["影响基因"] = annodb["INFO"].str.split("|").str[3]
            annodb["碱基变化"] = annodb["INFO"].str.split("|").str[9]
            annodb["氨基酸变化"] = annodb["INFO"].str.split("|").str[10]
            annodb["突变位置"] = annodb["POS"]
            annodb["片段名称"] = annodb["#CHROM"]
            annodb = annodb[["片段名称", "突变位置", "参考碱基", "突变碱基", "影响基因", "突变类型", "突变影响", "碱基变化", "氨基酸变化"]]
            annodb["氨基酸变化"] = annodb.apply(lambda x: "-" if x["氨基酸变化"] == "" else x["氨基酸变化"], axis=1)
            annodb.fillna("-", inplace=True)
            annodb.to_csv(f"{Pre}.anno.tsv", sep="\t", index=False)
        else:
            subprocess.run("ln -s snps.filt1.vcf snps.anno.vcf", shell=True)
            skinums = int(os.popen("grep '##' snps.anno.vcf |wc -l").read())
            annodb = pd.read_table("snps.anno.vcf", skiprows=skinums)
            annodb["参考碱基"] = annodb["REF"] + "(" + annodb["unknown"].str.split("|").str[-1].str.strip().str.split(":").str[-5] + ")"
            annodb["突变碱基"] = annodb["ALT"] + "(" + annodb["unknown"].str.split("|").str[-1].str.strip().str.split(":").str[-3] + ")"
            annodb["突变位置"] = annodb["POS"]
            annodb["片段名称"] = annodb["#CHROM"]
            annodb = annodb[["片段名称", "突变位置", "参考碱基", "突变碱基"]]
            annodb.to_csv(f"{Pre}.anno.tsv", sep="\t", index=False)
    elif method == "clair3":
        with open("clair3.log", "w") as clg:
            subprocess.run(f"conda run -n clair3 --no-capture-output run_clair3.sh --bam_fn=ref.mapping.bam --ref_fn=genomes/ref.fa --threads={threads} --platform=ont --model_path=/home/dell/miniconda3/envs/clair3/bin/models/r941_prom_sup_g5014 --output=./ --include_all_ctgs --enable_long_indel --snp_min_af=0.05", shell=True, stdout=clg, stderr=clg)
        subprocess.run("samtools view -h -F 2308 ref.mapping.bam |samtools sort -o ref.filter.bam", shell=True)
        subprocess.run("samtools index ref.filter.bam", shell=True)
        subprocess.run(f"perbase base-depth ref.filter.bam -t {threads} > perbase.bed", shell=True)
        pd.read_table("perbase.bed")
        if ifAvcf:
            subprocess.run("java -jar /home/dell/biosoft/snpEff/snpEff.jar ann  -c snpEff.config -Datadir . ref merge_output.vcf.gz > snps.anno.vcf", shell=True)
            skinums = int(os.popen("grep '##' snps.anno.vcf |wc -l").read())
            annodb = pd.read_table("snps.anno.vcf", skiprows=skinums)
            annodb["参考碱基"] = annodb["REF"] + "(" + annodb["SAMPLE"].str.split("|").str[-1].str.strip().str.split(":").str[-2].str.split(",").str[0] + ")"
            annodb["突变碱基"] = annodb["ALT"] + "(" + annodb["SAMPLE"].str.split("|").str[-1].str.strip().str.split(":").str[-2].str.split(",").str[1] + ")"
            annodb["突变类型"] = annodb["INFO"].str.split("|").str[1]
            annodb["突变影响"] = annodb["INFO"].str.split("|").str[2]
            annodb["影响基因"] = annodb["INFO"].str.split("|").str[3]
            annodb["碱基变化"] = annodb["INFO"].str.split("|").str[9]
            annodb["氨基酸变化"] = annodb["INFO"].str.split("|").str[10]
            annodb["突变位置"] = annodb["POS"]
            annodb["片段名称"] = annodb["#CHROM"]
            annodb = annodb[["片段名称", "突变位置", "参考碱基", "突变碱基", "影响基因", "突变类型", "突变影响", "碱基变化", "氨基酸变化"]]
            annodb["氨基酸变化"] = annodb.apply(lambda x: "-" if x["氨基酸变化"] == "" else x["氨基酸变化"], axis=1)
            annodb.to_csv(f"{Pre}.anno.tsv", sep="\t", index=False)
        else:
            subprocess.run("gunzip merge_output.vcf.gz -c > snps.anno.vcf", shell=True)
            skinums = int(os.popen("grep '##' snps.anno.vcf |wc -l").read())
            annodb = pd.read_table("snps.anno.vcf", skiprows=skinums)
            annodb["参考碱基"] = annodb["REF"] + "(" + annodb["SAMPLE"].str.split("|").str[-1].str.strip().str.split(":").str[-2].str.split(",").str[0] + ")"
            annodb["突变碱基"] = annodb["ALT"] + "(" + annodb["SAMPLE"].str.split("|").str[-1].str.strip().str.split(":").str[-2].str.split(",").str[1] + ")"
            annodb["突变位置"] = annodb["POS"]
            annodb["片段名称"] = annodb["#CHROM"]
            annodb = annodb[["片段名称", "突变位置", "参考碱基", "突变碱基"]]
            annodb.to_csv(f"{Pre}.anno.tsv", sep="\t", index=False)
    subprocess.run("bcftools convert -Oz -o snps.anno.vcf.gz snps.anno.vcf", shell=True)
    subprocess.run("bcftools index -f snps.anno.vcf.gz", shell=True)
    if outputfa != "noforce":
        subprocess.run(f"bcftools consensus -f genomes/ref.fa -o {outputfa} snps.anno.vcf.gz -m mask.bed", shell=True)


def renamefa(inf, ofn):
    subprocess.run(f"seqkit sort -l -r {inf}|seqkit fx2tab > tmpfa.tab", shell=True)
    afile = pd.read_table("tmpfa.tab", header=None)
    afile["contignum"] = afile.index + 1
    if afile.loc[afile[1].str.len() > 1000, :].shape[0] > 0:
        afile = afile.loc[afile[1].str.len() > 1000, :]
    afile[3] = afile.apply(lambda x: f"contig_{x.contignum}", axis=1)
    afile[[0, 3]].to_csv("transname.tsv", sep="\t", index=False, header=False)
    afile[0] = afile.apply(lambda x: f">contig_{x.contignum}", axis=1)
    afile = afile[[0, 1]]
    afile.to_csv(ofn, sep="\n", index=False, header=False)


def prepare_assembly_inputs(inf, fq1, fq2, threads, Pre, lelID, method):
    runtime = get_runtime_context()
    krdb = runtime.krdb
    with open("tmpkk2.log", "w") as kkf:
        if not os.path.isfile("2.1.fastq"):
            if inf:
                subprocess.run(f"kraken2 --db {krdb} --threads {threads} --output {Pre}.out.txt --report {Pre}.report.txt {inf}", shell=True, stdout=kkf, stderr=kkf)
                subprocess.run(f"bracken -d {krdb} -o {Pre}.bracken1.txt -w {Pre}.bracken2.txt -l S -t 10  -i {Pre}.report.txt", shell=True, stdout=kkf, stderr=kkf)

            if fq1 and fq2:
                subprocess.run(f"kraken2 --db {krdb} --threads {threads} --output {Pre}_2.out.txt --report {Pre}_2.report.txt {fq1} {fq2}", shell=True, stdout=kkf, stderr=kkf)
                subprocess.run(f"bracken -d {krdb} -o {Pre}_2.bracken1.txt -w {Pre}_2.bracken2.txt -l S -t 10  -i {Pre}_2.report.txt", shell=True, stdout=kkf, stderr=kkf)
                _run_bracken_sub(krdb, f"{Pre}_2.report.txt", f"{Pre}_2", kkf)
            elif fq1:
                subprocess.run(f"kraken2 --db {krdb} --threads {threads} --output {Pre}_2.out.txt --report {Pre}_2.report.txt {fq1}", shell=True, stdout=kkf, stderr=kkf)
                subprocess.run(f"bracken -d {krdb} -o {Pre}_2.bracken1.txt -w {Pre}_2.bracken2.txt -l S -t 10  -i {Pre}_2.report.txt", shell=True, stdout=kkf, stderr=kkf)
                _run_bracken_sub(krdb, f"{Pre}_2.report.txt", f"{Pre}_2", kkf)

    pre2 = f"{Pre}_2" if (fq1 or fq2) else Pre
    if lelID != "nolevel":
        level = lelID.split(",")[1]
        krakenfile = f"{pre2}.report.txt"
        tkid = lelID.split(",")[0]
        taxlist1 = [int(i) for i in proc_kra1(krakenfile, tkid, level)]
        if fq2:
            exreadsID1(taxlist1, f"{pre2}.out.txt", f"{Pre}.R1.fastq.gz", f"{Pre}.R2.fastq.gz")
            fq1 = f"{tkid}.1.fastq"
            fq2 = f"{tkid}.2.fastq"
        else:
            exreadsID1(taxlist1, f"{pre2}.out.txt", f"{Pre}.R1.fastq.gz", 0)
            fq1 = f"{tkid}.1.fastq"
            fq2 = 0
    elif method == "meta" and not os.path.isfile("2.1.fastq"):
        taxlist1 = [int(i) for i in proc_kra1(f"{pre2}.report.txt", 2, "D")]
        if fq2:
            exreadsID1(taxlist1, f"{pre2}.out.txt", f"{Pre}.R1.fastq.gz", f"{Pre}.R2.fastq.gz")
            fq1 = "2.1.fastq"
            fq2 = "2.2.fastq"
        else:
            exreadsID1(taxlist1, f"{pre2}.out.txt", f"{Pre}.R1.fastq.gz", 0)
            fq1 = "2.1.fastq"
            fq2 = 0

    return fq1, fq2, pre2


def map_assembly_reads(finalfa, inf, fq1, fq2, threads, Pre, asmt, long_type, logf):
    if asmt in ["longasm", "longref"]:
        if long_type == "Nanopore":
            subprocess.run(f"minimap2 -ax map-ont {finalfa} {inf} -t {threads} |samtools sort -o {Pre}.sorted.bam", shell=True, stdout=logf, stderr=logf)
        elif long_type == "PacBio_CLR":
            subprocess.run(f"minimap2 -ax map-pb {finalfa} {inf} -t {threads} |samtools sort -o {Pre}.sorted.bam", shell=True, stdout=logf, stderr=logf)
        elif long_type == "PacBio_CCS":
            subprocess.run(f"minimap2 -ax map-hifi {finalfa} {inf} -t {threads} |samtools sort -o {Pre}.sorted.bam", shell=True, stdout=logf, stderr=logf)
    elif asmt in ["shortasm", "shortref"]:
        subprocess.run(f"bwa index {finalfa}", shell=True)
        if fq2:
            subprocess.run(f"bwa mem {finalfa} {fq1} {fq2} -t {threads}|samtools sort -o {Pre}_ngs.sorted.bam", shell=True)
        else:
            subprocess.run(f"bwa mem {finalfa} {fq1} -t {threads}|samtools sort -o {Pre}_ngs.sorted.bam", shell=True)
    else:
        if long_type == "Nanopore":
            subprocess.run(f"minimap2 -ax map-ont {finalfa} {inf} -t {threads} |samtools sort -o {Pre}.sorted.bam", shell=True, stdout=logf, stderr=logf)
        elif long_type == "PacBio_CLR":
            subprocess.run(f"minimap2 -ax map-pb {finalfa} {inf} -t {threads} |samtools sort -o {Pre}.sorted.bam", shell=True, stdout=logf, stderr=logf)
        elif long_type == "PacBio_CCS":
            subprocess.run(f"minimap2 -ax map-hifi {finalfa} {inf} -t {threads} |samtools sort -o {Pre}.sorted.bam", shell=True, stdout=logf, stderr=logf)
        if fq2:
            subprocess.run(f"bwa mem {finalfa} {fq1} {fq2} -t {threads}|samtools sort -o {Pre}_ngs.sorted.bam", shell=True)
        elif fq1:
            subprocess.run(f"bwa mem {finalfa} {fq1} -t {threads}|samtools sort -o {Pre}_ngs.sorted.bam", shell=True)

    _generate_depth_files(Pre, logf)


def build_assembly_info(finalfa, Pre, method):
    canudb = None
    if method == "canu":
        subprocess.run("seqkit fx2tab canu/canu.contigs.fasta -n > canu.txt", shell=True)
        canudb = pd.read_table("canu.txt", sep=" ", header=None)
        canudb["len"] = canudb[1].str.replace("len=", "").str.strip().astype("int")
        canudb = canudb.sort_values("len", ascending=False)
        canudb["index"] = canudb.index + 1
        canudb["contig"] = "contig_" + canudb["index"].astype("str")
        canudb["cir"] = canudb[6].str.replace("suggestCircular=", "")

    if not os.path.exists("flye_output"):
        os.makedirs("flye_output")
        subprocess.run(f"seqkit fx2tab {finalfa} -n -l > flye_output/tmp.stat.txt", shell=True)
        flyedb = pd.read_table("flye_output/tmp.stat.txt", header=None)
        flyedb["是否成环"] = "-"
        flyedb.rename(columns={0: "序列名称", 1: "序列长度"}, inplace=True)
        if method == "canu" and canudb is not None:
            for ctg in flyedb["序列名称"].tolist():
                flyedb.loc[flyedb["序列名称"] == ctg, "是否成环"] = canudb.loc[canudb["contig"] == ctg, "cir"].tolist()[0]

        mosdb = pd.read_table(f"{Pre}.regions.bed", header=None) if os.path.isfile(f"{Pre}.regions.bed") else pd.read_table(f"{Pre}_ngs.regions.bed", header=None)
        flyedb = mosdb.groupby(0).agg("mean").merge(flyedb, left_on=0, right_on="序列名称")
        flyedb.rename(columns={3: "平均深度"}, inplace=True)
        flyedb = flyedb[["序列名称", "序列长度", "平均深度", "是否成环"]]
        flyedb["平均深度"] = flyedb["平均深度"].round()
        flyedb.sort_values("序列长度", axis=0, inplace=True, ascending=False)
        flyedb.to_csv("flye_output/assembly_info.txt", sep="\t", index=False)
        return

    transdb = pd.read_table("transname.tsv", header=None, names=["oldname", "newname"])
    subprocess.run(f"seqkit fx2tab {finalfa} -n -l > flye_output/tmp.stat.txt", shell=True)
    newcontigdb = pd.read_table("flye_output/tmp.stat.txt", header=None)
    flyedb = pd.read_table("flye_output/assembly_info.txt")
    flyedb = flyedb.merge(transdb, left_on="序列名称", right_on="oldname")
    flyedb = flyedb.merge(newcontigdb, left_on="newname", right_on=0)
    flyedb["序列名称"] = flyedb["newname"]
    flyedb["序列长度"] = flyedb[1]
    flyedb = flyedb.sort_values("序列长度", ascending=False)
    flyedb = flyedb[["序列名称", "序列长度", "平均深度", "是否成环"]]
    flyedb.to_csv("flye_output/assembly_info.txt", sep="\t", index=False)


def write_finalfasta_stats(finalfa, Pre):
    subprocess.run(f"seqkit stat -a -T -G N {finalfa} > ragtag_sum.tsv", shell=True)
    ragtagdb = pd.read_table("ragtag_sum.tsv")
    ragtagdb = ragtagdb[["num_seqs", "sum_len", "max_len", "min_len", "sum_gap", "N50", "GC(%)"]]
    ragtagdb["样本名称"] = Pre
    ragtagdb.rename(columns={"num_seqs": "contig数量", "sum_len": "总长度", "min_len": "最小contig长度", "max_len": "最大contig长度"}, inplace=True)
    ragtagdb["N比例(%)"] = ((ragtagdb["sum_gap"] / ragtagdb["总长度"]) * 100).round(2)
    ragtagdb = ragtagdb[["样本名称", "contig数量", "总长度", "最大contig长度", "最小contig长度", "N50", "GC(%)", "N比例(%)"]]
    ragtagdb.to_csv("finalfasta.tsv", index=False, sep="\t")


def build_fasta_assembly_info(finalfa, Pre, method):
    canudb = None
    if method == "canu":
        subprocess.run("seqkit fx2tab canu/canu.contigs.fasta -n > canu.txt", shell=True)
        canudb = pd.read_table("canu.txt", sep=" ", header=None)
        canudb["len"] = canudb[1].str.replace("len=", "").str.strip().astype("int")
        canudb = canudb.sort_values("len", ascending=False)
        canudb["index"] = canudb.index + 1
        canudb["contig"] = "contig_" + canudb["index"].astype("str")
        canudb["cir"] = canudb[6].str.replace("suggestCircular=", "")

    if not os.path.exists("flye_output"):
        os.makedirs("flye_output")
    subprocess.run(f"seqkit fx2tab {finalfa} -n -l > flye_output/tmp.stat.txt", shell=True)
    flyedb = pd.read_table("flye_output/tmp.stat.txt", header=None)
    flyedb["是否成环"] = "-"
    flyedb.rename(columns={0: "序列名称", 1: "序列长度"}, inplace=True)
    if method == "canu" and canudb is not None:
        for ctg in flyedb["序列名称"].tolist():
            flyedb.loc[flyedb["序列名称"] == ctg, "是否成环"] = canudb.loc[canudb["contig"] == ctg, "cir"].tolist()[0]
    flyedb["平均深度"] = "-"
    flyedb = flyedb[["序列名称", "序列长度", "平均深度", "是否成环"]]
    flyedb.sort_values("序列长度", axis=0, inplace=True, ascending=False)
    flyedb.to_csv("flye_output/assembly_info.txt", sep="\t", index=False)
    return flyedb


def annotate_with_prokka(finalfa, Pre, threads, log_path="asb.log"):
    with open(log_path, "a") as logf:
        subprocess.run(f"prokka --force --outdir {Pre}_prokka --prefix {Pre} --addgenes --cpus {threads} {finalfa}", shell=True, stdout=logf, stderr=logf)

    prkskpn = int(os.popen(f"grep '##sequence-region' {Pre}_prokka/{Pre}.gff|wc -l").read()) + 1
    prkfan = int(os.popen(f"grep -n '#FASTA' {Pre}_prokka/{Pre}.gff").read().split(":")[0])
    prkskfn = int(os.popen(f"cat {Pre}_prokka/{Pre}.gff|wc -l").read()) - prkfan + 1
    prokkadb = pd.read_table(
        f"{Pre}_prokka/{Pre}.gff",
        skiprows=prkskpn,
        skipfooter=prkskfn,
        engine="python",
        header=None,
        names=["染色体", "数据库", "类型", "起始位置", "终止位置", "t1", "链方向", "t2", "注释1"],
    )
    prokkadb = prokkadb[["染色体", "类型", "起始位置", "终止位置", "链方向", "注释1"]]
    prokkadb["基因名称"] = prokkadb["注释1"].str.split("Name=").str[1].str.split(";").str[0]
    prokkadb["locus标签"] = prokkadb["注释1"].str.split("ID=").str[1].str.split(";").str[0]
    prokkadb.fillna("-", inplace=True)
    prokkadb = prokkadb[["染色体", "类型", "起始位置", "终止位置", "链方向", "基因名称", "locus标签"]]
    prokkadb.to_csv(f"{Pre}.prokka.tsv", sep="\t", index=False)

    monthchin = int(os.popen(f"grep '月' {Pre}_prokka/{Pre}.gbk|wc -l").read().strip())
    if monthchin:
        if os.popen(f"head -n 1 {Pre}_prokka/{Pre}.gbk").read().strip()[72] == "月":
            subprocess.run(f"sed -i 's/月/  /g' {Pre}_prokka/{Pre}.gbk", shell=True)
        else:
            subprocess.run(f"sed -i 's/月/ /g' {Pre}_prokka/{Pre}.gbk", shell=True)
    subprocess.run(f"cp {Pre}_prokka/{Pre}.gbk tt.gbk", shell=True)


def enhance_plasmid_results(finalfa, Pre, log_path="asb.log"):
    with open(log_path, "a") as logf:
        subprocess.run(f"/home/dell/miniconda3/bin/conda run --no-capture-output -n plasflow PlasFlow.py --input {finalfa} --output {Pre}_plaspredict.tsv", shell=True, stdout=logf, stderr=logf)
        subprocess.run(f"staramr search  {finalfa} -o staramr_result -n 30", shell=True, stdout=logf, stderr=logf)

    plasmiddb = pd.read_table("staramr_result/plasmidfinder.tsv")
    rawindexname = plasmiddb.index.name
    if plasmiddb.shape[0] > 1:
        plasmiddb = plasmiddb.groupby("Contig").apply(_join_plasmid)
        plasmiddb.index.name = rawindexname

    plasflowdb = pd.read_table(f"{Pre}_plaspredict.tsv")[["contig_name", "label"]]
    flyedb = pd.read_table("flye_output/assembly_info.txt")
    flyedb = flyedb.merge(plasflowdb, left_on="序列名称", right_on="contig_name").drop("contig_name", axis=1)
    flyedb = flyedb.merge(plasmiddb, left_on="序列名称", right_on="Contig", how="left").drop("Contig", axis=1)
    flyedb = flyedb.rename(columns={"label": "基因组/质粒", "Plasmid": "质粒分型"})
    flyedb.fillna("-", inplace=True)
    flyedb.to_csv("tmp1.tsv", sep="\t", index=False)
    flyedb["序列长度"] = flyedb["序列长度"].astype("int")
    plasmidlist = flyedb.loc[
        (flyedb["基因组/质粒"].str.contains("plasmid")) | ((flyedb["质粒分型"] != "-") & (flyedb["序列长度"] < 1000000)),
        "序列名称",
    ].tolist()
    flyedb["占比"] = (flyedb["序列长度"] / flyedb["序列长度"].sum()).round(2)
    flyedb = flyedb[["序列名称", "序列长度", "平均深度", "是否成环", "基因组/质粒", "质粒分型", "占比"]]
    flyedb.to_csv("flye_output/assembly_info.txt", sep="\t", index=False)
    return flyedb, plasmidlist


def export_plasmid_gbk_and_cgview(Pre, plasmidlist):
    input_gbk = f"{Pre}_prokka/{Pre}.gbk"
    if len(plasmidlist) != 0:
        for contig_id_to_extract in plasmidlist:
            for record in SeqIO.parse(input_gbk, "genbank"):
                if record.id == contig_id_to_extract:
                    with open(f"{Pre}_prokka/{contig_id_to_extract}.gbk", "w") as output_handle:
                        SeqIO.write(record, output_handle, "genbank")
                    with open("cgview.log", "a") as cgvf:
                        subprocess.run(f"ruby /data/deploy/bio-elite/bio/script/cgview_builder_cli.rb -s {Pre}_prokka/{contig_id_to_extract}.gbk -o {contig_id_to_extract}.json -c /data/deploy/bio-elite/bio/script/CGview.yaml -n {contig_id_to_extract}", shell=True, stdout=cgvf, stderr=cgvf)

        records = list(SeqIO.parse(input_gbk, "genbank"))
        filtered_records = [record for record in records if record.id not in plasmidlist]
        with open(f"{Pre}_prokka/main.gbk", "w") as output_handle:
            SeqIO.write(filtered_records, output_handle, "genbank")
        with open("cgview.log", "a") as cgvf:
            subprocess.run(f"ruby /data/deploy/bio-elite/bio/script/cgview_builder_cli.rb -s {Pre}_prokka/main.gbk -o main.json -c /data/deploy/bio-elite/bio/script/CGview.yaml -n main", shell=True, stdout=cgvf, stderr=cgvf)
    else:
        with open("cgview.log", "a") as cgvf:
            subprocess.run(f"ruby /data/deploy/bio-elite/bio/script/cgview_builder_cli.rb -s {Pre}_prokka/{Pre}.gbk -o main.json -c /data/deploy/bio-elite/bio/script/CGview.yaml -n {Pre}", shell=True, stdout=cgvf, stderr=cgvf)


def run_checkm_and_write_summary(Pre, finalfa, threads, species, flyedb):
    if not os.path.isdir(f"{Pre}_bin_genome_out"):
        os.makedirs(f"{Pre}_bin_genome_out")
    subprocess.run(f"cp {finalfa} {Pre}_bin_genome_out/{Pre}.fna", shell=True)

    with open("checkm2", "w") as cmf:
        subprocess.run(f"/home/dell/miniconda3/bin/conda run --no-capture-output  -n cm2 checkm2 predict -i  {Pre}_bin_genome_out -x fna -o checkm2_out -t {threads} --force --database_path /data1/shanghai_pip/meta_genome/uniref100.KO.1.dmnd", shell=True, stdout=cmf, stderr=cmf)
        if os.path.isfile("checkm2_out/quality_report.tsv"):
            checkmdb = pd.read_table("checkm2_out/quality_report.tsv")
            checkmdb["样本名称"] = Pre
            checkmdb.rename(columns={"Completeness": "完整性", "Contamination": "污染率"}, inplace=True)
            checkmdb["物种名称"] = species
            checkmdb[["样本名称", "物种名称", "污染率", "完整性"]].to_csv(f"{Pre}.checkm.tsv", sep="\t", index=False)

    subprocess.run(f"seqkit stat {finalfa} -b -T -a > {Pre}.fasum.tsv", shell=True)
    fadb = pd.read_table(f"{Pre}.fasum.tsv")
    assdict = {
        "Contig数量": fadb["num_seqs"].tolist()[0],
        "N50长度": fadb["N50"].tolist()[0],
        "组装基因组长度": fadb["sum_len"].tolist()[0],
        "最长片段长度": fadb["max_len"].tolist()[0],
        "污染率": "-",
        "完整性": "-",
        "主基因组是否成环": "-",
    }
    if os.path.isfile(f"{Pre}.checkm.tsv"):
        cdb = pd.read_table(f"{Pre}.checkm.tsv")
        assdict["污染率"] = round(float(cdb["污染率"].tolist()[0]), 2)
        assdict["完整性"] = round(float(cdb["完整性"].tolist()[0]), 2)
    assdb = pd.DataFrame(assdict, index=[0])
    assdb["样本名称"] = Pre
    assdb = assdb[["样本名称", "Contig数量", "N50长度", "最长片段长度", "主基因组是否成环", "污染率", "完整性"]]
    assdb["主基因组是否成环"] = flyedb["是否成环"].tolist()[0]
    assdb.to_csv(f"{Pre}.assemble.result.tsv", sep="\t", index=False)


def run_genovi_summary(Pre):
    with open("genovi.log", "w") as fv:
        try:
            subprocess.run(f"/home/dell/miniconda3/bin/conda run -n genovi genovi -i {Pre}_prokka/{Pre}.gbk -o {Pre}_genovi -s draft", shell=True, stdout=fv, stderr=fv)
            cogdb = pd.read_table(f"{Pre}_genovi/{Pre}_genovi_COG_Classification.csv", sep=",", header=1)
            cogdb.iloc[-2, :].to_csv("Cog_summary.tsv", sep="\t", header=False)
        except Exception:
            pass


def write_gene_summaries(Pre):
    tmpgenedb = pd.read_table(f"{Pre}_prokka/{Pre}.tsv")
    tmpgenedb = tmpgenedb[tmpgenedb.ftype == "gene"]
    tmpgenedb_dict = {}
    for minlg in range(0, 2000, 100):
        tmpgenedb_dict[f"{minlg}-{minlg+100}"] = sum(tmpgenedb.length_bp.between(minlg, minlg + 100))
    tmpgenedb_dict[">2000"] = sum(tmpgenedb.length_bp >= 2000)
    genedb = pd.DataFrame(tmpgenedb_dict, index=["Gene数量"]).T
    genedb["范围"] = genedb.index
    genedb = genedb[["范围", "Gene数量"]]
    genedb.to_csv(f"{Pre}_gene_raw_sum.tsv", sep="\t", index=False)

    gene_fundb = pd.read_table(f"{Pre}_prokka/{Pre}.txt", sep=":")
    gene_fundb.index = gene_fundb.organism
    gene_fundb = gene_fundb.drop("organism", axis=1).T
    gene_fundb.to_csv(f"{Pre}.genefun_summary.tsv", sep="\t", index=False)

    if os.path.isfile(f"{Pre}.uniqgene.fasta"):
        subprocess.run(f"seqkit fx2tab -n -l {Pre}.uniqgene.fasta > {Pre}.uniqgene.tsv", shell=True)
        tmpgeneudb = pd.read_table(f"{Pre}.uniqgene.tsv", names=["Genename", "length"])
        tmpgeneudb_dict = {}
        for minlg in range(0, 2000, 100):
            tmpgeneudb_dict[f"{minlg}-{minlg+100}"] = sum(tmpgeneudb.length.between(minlg, minlg + 100))
        tmpgeneudb_dict[">2000"] = sum(tmpgeneudb.length >= 2000)
        geneudb = pd.DataFrame(tmpgeneudb_dict, index=["Gene数量"]).T
        geneudb["范围"] = geneudb.index
        geneudb = geneudb[["范围", "Gene数量"]]
        geneudb.to_csv(f"{Pre}_gene_uniq_sum.tsv", sep="\t", index=False)

    subprocess.run(f"grep CDS {Pre}_prokka/{Pre}.gff > {Pre}.CDS.gff", shell=True)
    with open(f"{Pre}.CDS.gff") as f:
        open(f"{Pre}.Contig_gene.tsv", "w").write("Contig\tCDS\n")
        for line in f:
            line = line.strip().split("\t")
            if line[2] == "CDS":
                contig = line[0]
                cdsid = line[8].split("locus_tag=")[-1].split(";")[0]
                open(f"{Pre}.Contig_gene.tsv", "a").write(f"{contig}\t{cdsid}\n")

    subprocess.run(f"grep 'gene' {Pre}_prokka/{Pre}.gff > {Pre}.gene.gff", shell=True)
    with open(f"{Pre}.gene.gff") as f:
        open(f"{Pre}.gene.bed", "w").write("")
        for line in f:
            line = line.strip().split("\t")
            if line[2] == "gene":
                contig = line[0]
                cdsid = line[8].split("locus_tag=")[-1].split(";")[0]
                start = line[3]
                end = line[4]
                open(f"{Pre}.gene.bed", "a").write(f"{contig}\t{start}\t{end}\t{cdsid}\n")
    subprocess.run(f"samtools faidx {Pre}_prokka/{Pre}.fna", shell=True)
    subprocess.run(f"bedtools getfasta -fi {Pre}_prokka/{Pre}.fna -bed {Pre}.gene.bed -name > tmp_{Pre}.gene.fasta", shell=True)
    subprocess.run(f"cut -d ':' -f1  tmp_{Pre}.gene.fasta > {Pre}.gene.fasta", shell=True)


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
    ts = x["start"].min()
    te = x["end"].max()
    ofname = x["GeneName"].tolist()[0]
    Spename = x["Species"].tolist()[0]
    x["start"] = x.reset_index().index + 1
    x["end"] = x.reset_index().index + 2
    x[["Chrom", "start", "end", "Depth"]].to_csv(f"geneDepth/{ofname}_{typeF}.tsv", sep="\t", header=False, index=False)
    open(f"geneDepth/{ofname}_{typeF}.bed", "w").write(f"{x['Chrom'].tolist()[0]}\t{ts}\t{te}\t{Pre}_{ofname}\t{Pre}_{ofname}\t{x['strand'].tolist()[0]}")
    subprocess.run(f"bedtools getfasta -fi {Pre}.final.fasta -bed geneDepth/{ofname}_{typeF}.bed -name -s > geneDepth/{ofname}_{typeF}.fasta", shell=True)
    tmpdict = {"片段名称": x["Chrom"].tolist()[0], "物种名称": Spename, "起始位置": x["start"].min(), "终止位置": x["end"].max(), "覆盖度(>0)%": round(x[x["Depth"] > 0].shape[0] / x.shape[0], 4) * 100, "覆盖度(>10)%": round(x[x["Depth"] > 10].shape[0] / x.shape[0], 4) * 100, "覆盖度(>100)%": round(x[x["Depth"] > 100].shape[0] / x.shape[0], 4) * 100, "平均深度": round(x["Depth"].mean(), 2), "最低深度": x["Depth"].min(), "最高深度": x["Depth"].max()}
    return pd.DataFrame(tmpdict, index=[0]).round(2)


def _run_bracken_sub(krdb, report_path, prefix, kkf):
    testbrkdb = pd.read_table(report_path, header=None)
    if "S4" in testbrkdb[3]:
        level = "S3"
    elif "S3" in testbrkdb[3]:
        level = "S2"
    else:
        level = "S1"
    subprocess.run(f"bracken -d {krdb} -o {prefix}_Sub.bracken1.txt -w {prefix}_Sub.bracken2.txt -l {level} -t 10  -i {report_path}", shell=True, stdout=kkf, stderr=kkf)


def _write_contig_beds(prefix):
    if not os.path.isdir("Contigbedfile"):
        os.makedirs("Contigbedfile")
    contigbed = pd.read_table(f"{prefix}.regions.bed", header=None)
    contig1bed = pd.read_table(f"{prefix}_1.regions.bed", header=None)
    contigbed.groupby(0).apply(lambda x: x.to_csv(f"Contigbedfile/{x[0].tolist()[0]}.bed", index=False, header=False, sep="\t"))
    contig1bed.groupby(0).apply(lambda x: x.to_csv(f"Contigbedfile/{x[0].tolist()[0]}_dis1.bed", index=False, header=False, sep="\t"))
    for bed in os.listdir("Contigbedfile"):
        if not bed.endswith("_dis1.bed"):
            ttPre = bed.replace(".bed", "")
            if int(os.popen(f"cat Contigbedfile/{bed}|wc -l").read()) < 10:
                subprocess.run(f"mv Contigbedfile/{ttPre}_dis1.bed Contigbedfile/{ttPre}.bed ", shell=True)
    contig1bed[contig1bed[3] == 0].to_csv("mask.bed", sep="\t", index=False, header=False)


def _generate_depth_files(Pre, logf):
    if os.path.isfile(f"{Pre}.sorted.bam"):
        subprocess.run(f"samtools index {Pre}.sorted.bam", shell=True, stdout=logf, stderr=logf)
        subprocess.run(f"mosdepth -b 1000 {Pre} {Pre}.sorted.bam", shell=True)
        subprocess.run(f"mosdepth -b 1 {Pre}_1 {Pre}.sorted.bam", shell=True)
        subprocess.run(f"gunzip -f {Pre}.regions.bed.gz", shell=True)
        subprocess.run(f"gunzip -f {Pre}_1.regions.bed.gz", shell=True)
        _write_contig_beds(Pre)

    if os.path.isfile(f"{Pre}_ngs.sorted.bam"):
        subprocess.run(f"samtools index {Pre}_ngs.sorted.bam", shell=True, stdout=logf, stderr=logf)
        subprocess.run(f"mosdepth -b 1000 {Pre}_ngs {Pre}_ngs.sorted.bam", shell=True)
        subprocess.run(f"mosdepth -b 1 {Pre}_ngs_1 {Pre}_ngs.sorted.bam", shell=True)
        subprocess.run(f"gunzip -f {Pre}_ngs.regions.bed.gz", shell=True)
        subprocess.run(f"gunzip -f {Pre}_ngs_1.regions.bed.gz", shell=True)
        _write_contig_beds(f"{Pre}_ngs")


def _join_plasmid(group):
    plasmids = "|".join(group["Plasmid"])
    newdb = group
    newdb["Plasmid"] = plasmids
    return newdb.iloc[0, :]
