from __future__ import annotations

import os
import shutil
from pathlib import Path

from metagenomic_refactor.common import copy_pattern, run_cmd
from metagenomic_refactor.context import get_runtime_context


def combine_func(Pre):
    runtime = get_runtime_context()
    result_dir = Path(f"{Pre}_genome_complete_result")
    anno_dir = Path(f"{Pre}_anno_sum")

    with open("combine.log", "w") as logf:
        if result_dir.exists():
            shutil.rmtree(result_dir)
        result_dir.mkdir()

        copy_pattern(["*QC*.summary.tsv", "*.fastp*.json"], result_dir / "1.DataSum")
        copy_pattern(["*qc"], result_dir / "1.DataSum")
        copy_pattern(["*.list.txt", "*.krona.html"], result_dir / "2.Spereads")
        copy_pattern(["flye_output/assembly_info.txt", "*.checkm.tsv", "*_raw.png", "*.final.fasta"], result_dir / "3.Assemble")
        copy_pattern(["*.repeat.tsv", "*.mummer.*"], result_dir / "4.Repeat")
        copy_pattern(["*.rgi.tsv", "*.card.tsv", "*.vfdb.tsv", "*.phi.tsv", "*.CRISPR.tsv", "*.annot.tsv", "*.sigIP.tsv", "*.medb.tsv", "*.ncRNA.tsv"], result_dir / "5.Fun_Element")
        copy_pattern(["*mlst_Stat.txt", "*gene_show.txt"], result_dir / "7.Mlst")
        copy_pattern(["*.swiss.tsv", "*.GO.tsv", "*.KEGG.tsv", "*.PFAM.tsv", "*.Cog.tsv", "*.CAZy.tsv", "*.rgi.tsv", "*.card.tsv", "*.vfdb.tsv", "*.phi.tsv", "*.CRISPR.tsv", "*.annot.tsv", "*.sigIP.tsv", "*.medb.tsv", "*.ncRNA.tsv"], anno_dir)
        run_cmd(f"tar -zcf {Pre}.anno.tar.gz {anno_dir}", logf)

        if runtime.method != "meta":
            run_cmd(f"/home/dell/miniconda3/bin/conda run -n report_env Rscript /data1/shanghai_pip/meta_genome/report_asb.R {os.getcwd()}", logf)
        else:
            run_cmd(f"/home/dell/miniconda3/bin/conda run -n report_env Rscript /data/deploy/meta_genome/report_db/report_meta.R {os.getcwd()}", logf)

        meta_dir = Path("meta_out")
        meta_dir.mkdir(exist_ok=True)
        copy_pattern(["Summary_kraken.csv", f"{Pre}*list*.txt", "2.vfdb.tsv", "2.card.tsv", "summary.tsv", "test.fastp*.json"], meta_dir)
        run_cmd("rm -rf 2.listID*.txt 2_fqID.txt 2.*.sorted.bam 2.id.tsv *_t.R*.fastq.gz *_sub.R*.fastq *regions.bed", logf)

        if Path("2.1.fastq").exists():
            Path("2.1.fastq").rename(f"{Pre}D2.1.fastq")
            run_cmd(f"pigz -p 10 {Pre}D2.1.fastq", logf)
        if Path("2.2.fastq").exists():
            Path("2.2.fastq").rename(f"{Pre}D2.2.fastq")
            run_cmd(f"pigz -p 10 {Pre}D2.2.fastq", logf)
