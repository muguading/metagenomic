from __future__ import annotations

import os
import subprocess
from dataclasses import dataclass, field

import pandas as pd

from metagenomic_refactor.common import conda_env_path, conda_run_command


@dataclass
class ResourcePaths:
    kraken_db: str
    virus_kraken_db: str
    species_table: str
    species_ref_table: str
    vf_meta_table: str
    krona_script: str
    taxa_info: str
    all_sort: str
    all_spe_proid_rank: str
    wormbase: str
    card_sequences: str
    vfdb_fasta: str
    aro_index: str
    vfdb_meta_annotation: str
    vfdb_contig: str
    report_asb_script: str
    report_meta_script: str
    virsorter2_bin: str
    virsorter2_db: str
    checkv_bin: str
    checkv_db: str


@dataclass
class PipelineConfig:
    inf: str
    ofn: str
    sumpath: str
    nt: int
    rnalib: str
    krdb: str
    intype: str
    analysis_target: str
    minl: str
    minQ: str
    barkit: str
    tmpfake: int
    long_type: str
    method: str
    asm_type: str
    ref: str
    gtf: str
    runflow: str
    rmhost: str
    tspeabun: str
    speciesdb: pd.DataFrame
    IfAnno: str
    sc1: str
    refdict: str
    genome_len: str
    ptimes: str
    psoft: str
    fastq1: int
    fastq2: int
    species: str
    speciesrefdb: pd.DataFrame
    vfmeta: pd.DataFrame
    config_out: str
    mmethod: str
    resources: ResourcePaths = field(repr=False)


def _env_path(name: str, default: str) -> str:
    return os.path.abspath(os.environ.get(name, default))


def _database_root() -> str:
    raw = str(os.environ.get("META_DATABASE_ROOT") or "").strip()
    if raw:
        candidate = os.path.abspath(os.path.expanduser(raw))
        if os.path.basename(candidate.rstrip(os.sep)) == "database":
            return candidate
        return os.path.join(candidate, "database")
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "database"))


def _database_path(*parts: str) -> str:
    return os.path.join(_database_root(), *parts)


def _read_table_if_exists(path: str, *, encoding: str | None = None) -> pd.DataFrame:
    if path and os.path.isfile(path):
        kwargs = {"encoding": encoding} if encoding else {}
        return pd.read_table(path, **kwargs)
    return pd.DataFrame()


def build_resource_paths() -> ResourcePaths:
    return ResourcePaths(
        kraken_db=_env_path("META_KRAKEN_DB", _database_path("kraken2_custom_202101_24G")),
        virus_kraken_db=_env_path("META_VIRUS_KRAKEN_DB", _database_path("virus", "k2_viral_20260226")),
        species_table=_env_path("META_SPECIES_TABLE", _database_path("pathotable.tsv")),
        species_ref_table=_env_path("META_SPECIES_REF_TABLE", _database_path("gc_gtdbmeta.tsv")),
        vf_meta_table=_env_path("META_VF_META_TABLE", _database_path("vfdb", "VFs_meta.tsv")),
        krona_script=_env_path("META_KRONA_SCRIPT", _database_path("TB_soft", "other_soft", "3_kreport2krona.py")),
        taxa_info=_env_path("META_TAXA_INFO", _database_path("Meta_anno", "taxa_info_20210508.txt")),
        all_sort=_env_path("META_ALL_SORT", _database_path("Meta_anno", "All.sort.csv")),
        all_spe_proid_rank=_env_path("META_ALL_SPE_PROID_RANK", _database_path("Meta_anno", "AllSpeProid_rank.txt")),
        wormbase=_env_path("META_WORMBASE", _database_path("Meta_anno", "wormbase.tsv")),
        card_sequences=_env_path("META_CARD_SEQUENCES", conda_env_path("meta_main", "db", "card", "sequences")),
        vfdb_fasta=_env_path("META_VFDB_FASTA", _database_path("vfdb.fasta")),
        aro_index=_env_path("META_ARO_INDEX", _database_path("aro_index.tsv")),
        vfdb_meta_annotation=_env_path("META_VFDB_META_ANNOTATION", _database_path("VFs_meta.tsv")),
        vfdb_contig=_env_path("META_VFDB_CONTIG", _database_path("vfdb.contig.tsv")),
        report_asb_script=_env_path("META_REPORT_ASB_SCRIPT", _database_path("report_asb.R")),
        report_meta_script=_env_path("META_REPORT_META_SCRIPT", _database_path("report_db", "report_meta.R")),
        virsorter2_bin=os.environ.get("META_VIRSORTER2_BIN", conda_run_command("genomad_aux", "virsorter")),
        virsorter2_db=_env_path("META_VIRSORTER2_DB", _database_path("virsorter2")),
        checkv_bin=os.environ.get("META_CHECKV_BIN", conda_run_command("genomad_aux", "checkv")),
        checkv_db=_env_path("META_CHECKV_DB", _database_path("checkvDB", "checkv-db-v1.5")),
    )


def build_pipeline_config(argv, cpu_count: int) -> PipelineConfig:
    resources = build_resource_paths()
    if argv.input:
        inf = os.path.abspath(argv.input)
    else:
        raise SystemExit("未检测到输入数据")

    ofn = os.path.abspath(argv.output)
    sumpath = f"{ofn}/fastq_analysis"

    nt = argv.thread
    if nt > cpu_count:
        nt = cpu_count

    krdb = resources.virus_kraken_db if argv.analysis_target == "virus" else resources.kraken_db
    inputtypedic = {"1": ".cfg", "2": "fastq", "3": "barcode_fastq", "4": "fasta", "5": "@"}
    intype = inputtypedic.get(argv.inputtype, argv.inputtype)
    speciesdb = _read_table_if_exists(resources.species_table)

    ref = argv.ref
    refdict = ""
    if "ref" in argv.asm_type:
        if str(ref).strip() in {"", "noref", "None", "none"}:
            ref = "noref"
        elif os.path.isfile(ref):
            ref = os.path.abspath(ref)
        elif ref in ["salmonella", "E_coli", "Shigella", "Parahemolyticus", "cholerae", "Y_enterocolitica", "Campylobacter", "Brucella", "Lmono", "Kpne", "Suare", "Bcere", "Nmen", "HPinf"]:
            if speciesdb.empty:
                raise SystemExit(f"未找到物种参考映射表: {resources.species_table}")
            ref = speciesdb.loc[speciesdb["Species"] == ref, "reference"].tolist()[0]
        else:
            tref = _database_path("fadb", f"{ref}_genomic.fna.gz")
            subprocess.run(f"seqkit seq {tref} > tmp_ref.fa", shell=True)
            ref = f"{os.getcwd()}/tmp_ref.fa"

    if argv.analysis_target == "virus":
        speciesrefdb = _read_table_if_exists(resources.species_ref_table)
        vfmeta = _read_table_if_exists(resources.vf_meta_table, encoding="Windows-1252")
    else:
        speciesrefdb = pd.read_table(resources.species_ref_table)
        vfmeta = pd.read_table(resources.vf_meta_table, encoding="Windows-1252")

    config_out = f"""输入文件\t{inf}
输出文件\t{ofn}
分析对象\t{argv.analysis_target}
输入文件类型\t{intype}
线程数:\t{nt}
最小序列长度:\t{argv.minlongfilt}
最小序列质量:\t{argv.Qfilt}
试剂盒:\t{argv.barcodekit}
测序数据类型:\t{argv.asm_type}
长读长类型:\t{argv.long_type}
组装方法:\t{argv.method}
预估基因组大小:\t{argv.genome_len}
polish软件:\t{argv.polish_soft}
polish次数:\t{argv.polish_times}
物种:\t{argv.species}
参考基因组:\t{ref}
运行模块:\t{argv.runflow}
"""

    return PipelineConfig(
        inf=inf,
        ofn=ofn,
        sumpath=sumpath,
        nt=nt,
        rnalib=argv.rna,
        krdb=krdb,
        intype=intype,
        analysis_target=argv.analysis_target,
        minl=argv.minlongfilt,
        minQ=argv.Qfilt,
        barkit=argv.barcodekit,
        tmpfake=argv.fake_pip,
        long_type=argv.long_type,
        method=argv.method,
        asm_type=argv.asm_type,
        ref=ref,
        gtf=argv.gtf,
        runflow=argv.runflow,
        rmhost=argv.rmhost,
        tspeabun=argv.abun,
        speciesdb=speciesdb,
        IfAnno=argv.ifAnno,
        sc1=resources.krona_script,
        refdict=refdict,
        genome_len=argv.genome_len,
        ptimes=argv.polish_times,
        psoft=argv.polish_soft,
        fastq1=0,
        fastq2=0,
        species=argv.species,
        speciesrefdb=speciesrefdb,
        vfmeta=vfmeta,
        config_out=config_out,
        mmethod=argv.method,
        resources=resources,
    )
