from __future__ import annotations

import os
import shutil
from pathlib import Path

from metagenomic_refactor.common import CommandExecutionError, conda_run_command, copy_pattern, run_cmd
from metagenomic_refactor.context import get_runtime_context


def combine_func(Pre):
    runtime = get_runtime_context()
    result_dir = Path(f"{Pre}_genome_complete_result")
    anno_dir = Path(f"{Pre}_anno_sum")
    resources = runtime.resources

    with open("combine.log", "w") as logf:
        if result_dir.exists():
            shutil.rmtree(result_dir)
        result_dir.mkdir()

        copy_pattern(["*QC*.summary.tsv", "*.fastp*.json"], result_dir / "1.DataSum")
        copy_pattern(["*qc"], result_dir / "1.DataSum")
        copy_pattern(["*.list.txt", "*.krona.html", "*_nextclade/*.tsv", "*_nextclade/*.log"], result_dir / "2.Spereads")
        copy_pattern(["flye_output/assembly_info.txt", "*.checkm.tsv", "*_raw.png", "*.final.fasta"], result_dir / "3.Assemble")
        copy_pattern(["viral_assembly/*.tsv", "viral_assembly/*.fa", "viral_assembly/megahit_output/*.fa", "viral_assembly/virsorter2/*.tsv", "viral_assembly/checkv/*.tsv"], result_dir / "3.Assemble")
        copy_pattern(["*.repeat.tsv", "*.mummer.*"], result_dir / "4.Repeat")
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
                "*.ncRNA.tsv",
                "snps.raw.vcf",
                "snps.raw.ann.vcf",
                "snps.raw.mutation_table.tsv",
                "snps.raw.mutation_table.json",
                "*_hiv_resistance.tsv",
                "*_hiv_resistance.json",
            ],
            result_dir / "5.Fun_Element",
        )
        copy_pattern(["*mlst_Stat.txt", "*gene_show.txt"], result_dir / "7.Mlst")
        copy_pattern(["*.swiss.tsv", "*.GO.tsv", "*.KEGG.tsv", "*.PFAM.tsv", "*.Cog.tsv", "*.CAZy.tsv", "*.rgi.tsv", "*.card.tsv", "*.vfdb.tsv", "*.phi.tsv", "*.CRISPR.tsv", "*.annot.tsv", "*.sigIP.tsv", "*.medb.tsv", "*.ncRNA.tsv"], anno_dir)
        run_cmd(f"tar -zcf {Pre}.anno.tar.gz {anno_dir}", logf)
        #if runtime.analysis_target != 'virus':
        #    if runtime.method != "meta":
        #        report_script = resources.report_asb_script if resources is not None else "/data1/shanghai_pip/meta_genome/report_asb.R"
        #    else:
        #        report_script = resources.report_meta_script if resources is not None else "/data/deploy/meta_genome/report_db/report_meta.R"
        #    if not Path(report_script).is_file():
        #        raise FileNotFoundError(f"报告脚本不存在: {report_script}，请检查部署环境或配置环境变量")
        #    try:
        #        run_cmd(conda_run_command("report_env", f"Rscript {report_script} {os.getcwd()}"), logf)
        #    except CommandExecutionError as exc:
        #        raise RuntimeError(f"报告生成失败，请查看 {Path('combine.log').resolve()}；失败命令: {exc.command}") from exc

        meta_dir = Path("meta_out")
        meta_dir.mkdir(exist_ok=True)
        copy_pattern(["Summary_kraken.csv", f"{Pre}*list*.txt", "2.vfdb.tsv", "2.card.tsv", "summary.tsv", "test.fastp*.json"], meta_dir)
        run_cmd("rm -rf 2.listID*.txt 2_fqID.txt 2.*.sorted.bam 2.id.tsv *_t.R*.fastq.gz *_sub.R*.fastq", logf)

        if Path("2.1.fastq").exists():
            Path("2.1.fastq").rename(f"{Pre}D2.1.fastq")
            run_cmd(f"pigz -p 10 {Pre}D2.1.fastq", logf)
        if Path("2.2.fastq").exists():
            Path("2.2.fastq").rename(f"{Pre}D2.2.fastq")
            run_cmd(f"pigz -p 10 {Pre}D2.2.fastq", logf)
