from __future__ import annotations

import argparse
import os
from dataclasses import dataclass, field

import pandas as pd


def _env_path(name: str, default: str) -> str:
    return os.path.abspath(os.environ.get(name, default))


@dataclass
class ResourcePaths:
    species_table: str
    fameta_table: str
    self_db_root: str
    pathogen_db_root: str
    pathogen_snippy_fasta_root: str
    ref_fadb_dir: str
    bakta_db: str
    cgmlst_org_db: str
    cgmlst_self_org_db: str
    cgmlst_other_org_db: str
    mysql_host: str
    mysql_user: str
    mysql_password: str
    mysql_database: str


@dataclass
class PipelineConfig:
    input_path: str
    species: str
    threads: int
    output_dir: str
    ref: str
    meta: str | None
    cgmlstana: str
    gubbins: str
    msamethod: str
    treemethod: str
    bootstrap: int
    mltype: str
    mode: str
    cgmlstversion: str
    speciesdb: pd.DataFrame
    fametadb: pd.DataFrame
    config_out: str
    resources: ResourcePaths = field(repr=False)


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Phylogenetic_Tree")
    parser.add_argument("--species", "-s", type=str, default="salmonella", help="可选数据类型:[salmonella,E_coli,Shigella,Parahemolyticus,cholerae,Y_enterocolitica,Campylobacter]")
    parser.add_argument("--input", "-i", type=str, default=False, help="样本名\\tR1(fa)\\tR2")
    parser.add_argument("--threads", "-t", type=int, default=10, help="线程数量,默认线程：10")
    parser.add_argument("--output", "-o", type=str, help="输出文件夹，默认在当前运行路径下生成：output文件夹")
    parser.add_argument("--ref", "-r", type=str, default=False, help="参考基因组")
    parser.add_argument("--meta", "-m", type=str, help="自建库的列表信息")
    parser.add_argument("--cgmlstana", "-c", type=str, default="no", help="可选[yes,no]")
    parser.add_argument("--gubbins", "-gb", type=str, default="yes", help="是否使用 gubbins，可选[yes,no]")
    parser.add_argument("--msamethod", "-mr", type=str, default="snippy", help="可选[snippy,roray,ska]")
    parser.add_argument("--treemethod", "-tr", type=str, default="ML", help="可选[ML,NJ,MP]")
    parser.add_argument("--Bootstrap", "-bs", type=int, default=1000, help="bootstrap 次数")
    parser.add_argument("--mltype", "-mt", type=str, default="MFP", help="碱基替代先验模型")
    parser.add_argument("--mode", "-md", type=str, default="P", help="溯源数据库")
    parser.add_argument("--cgmlst", "-cg", type=str, default="none", help="cgmlst版本")
    return parser


def build_resource_paths() -> ResourcePaths:
    return ResourcePaths(
        species_table=_env_path("PATHOSOURCE_SPECIES_TABLE", "/data1/shanghai_pip/meta_genome/pathotable.tsv"),
        fameta_table=_env_path("PATHOSOURCE_FAMETA_TABLE", "/data1/shanghai_pip/meta_genome/gc_gtdbmeta.tsv"),
        self_db_root=_env_path("PATHOSOURCE_SELF_DB_ROOT", "/data/deploy/bio-elite/bio/load_file/pathogenic_Self_DB"),
        pathogen_db_root=_env_path("PATHOSOURCE_PATHOGEN_DB_ROOT", "/data/deploy/bio-elite/bio/load_file/pathodb/Pathogen_DB"),
        pathogen_snippy_fasta_root=_env_path("PATHOSOURCE_SNIPPY_FASTA_ROOT", "/data/deploy/bio-elite/bio/load_file/pathogenic/snippy/fasta"),
        ref_fadb_dir=_env_path("PATHOSOURCE_REF_FADB_DIR", "/data/deploy/meta_genome/database/fadb"),
        bakta_db=_env_path("PATHOSOURCE_BAKTA_DB", "/data1/shanghai_pip/meta_genome/database/baktk/db-light"),
        cgmlst_org_db=_env_path("PATHOSOURCE_CGMLST_ORG_DB", "/data/deploy/bio-elite/bio/load_file/cgmlstdb/cgmlstorgdb"),
        cgmlst_self_org_db=_env_path("PATHOSOURCE_CGMLST_SELF_DB", "/data/deploy/bio-elite/bio/load_file/cgmlstdb/cgmlstorgdb"),
        cgmlst_other_org_db=_env_path("PATHOSOURCE_CGMLST_OTHER_DB", "/data/deploy/bio-elite/bio/load_file/pathocgmlst"),
        mysql_host=os.environ.get("PATHOSOURCE_MYSQL_HOST", "127.0.0.1"),
        mysql_user=os.environ.get("PATHOSOURCE_MYSQL_USER", "baiyi"),
        mysql_password=os.environ.get("PATHOSOURCE_MYSQL_PASSWORD", "baiyi123@+1s"),
        mysql_database=os.environ.get("PATHOSOURCE_MYSQL_DATABASE", "baiyi"),
    )


def build_pipeline_config(argv, cpu_count: int) -> PipelineConfig:
    resources = build_resource_paths()
    if not argv.input:
        raise SystemExit("未检测到输入数据")
    input_path = os.path.abspath(argv.input)
    output_dir = os.path.abspath(argv.output or "output")
    threads = min(int(argv.threads), cpu_count)
    speciesdb = pd.read_table(resources.species_table)
    fametadb = pd.read_table(resources.fameta_table)
    msa_alias = {"ska": "ska2", "ska2": "ska2", "roary": "roray", "roray": "roray"}
    msamethod = msa_alias.get(argv.msamethod, argv.msamethod)
    config_out = f""" --input {input_path}
--threads {threads}
--output {output_dir}
--species {argv.species}
--meta {argv.meta}
--ref {argv.ref}
--cgmlstana {argv.cgmlstana}
--gubbins {argv.gubbins}
--treemethod {argv.treemethod}
--msamethod {msamethod}
--Bootstrap {argv.Bootstrap}
--mltype {argv.mltype}
--cgmlst {argv.cgmlst}
"""
    return PipelineConfig(
        input_path=input_path,
        species=argv.species,
        threads=threads,
        output_dir=output_dir,
        ref=argv.ref,
        meta=argv.meta,
        cgmlstana=argv.cgmlstana,
        gubbins=str(argv.gubbins or "yes").strip().lower() or "yes",
        msamethod=msamethod,
        treemethod=argv.treemethod,
        bootstrap=argv.Bootstrap,
        mltype=argv.mltype,
        mode=argv.mode,
        cgmlstversion=argv.cgmlst,
        speciesdb=speciesdb,
        fametadb=fametadb,
        config_out=config_out,
        resources=resources,
    )
