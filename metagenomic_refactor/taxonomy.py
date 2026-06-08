from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pandas as pd

from metagenomic_refactor.common import conda_env_path, run_command
from metagenomic_refactor.context import get_runtime_context

TAXA_INFO_PATH = "/data/Ref/Meta_anno/taxa_info_20210508.txt"
ALL_SORT_PATH = "/data/Ref/Meta_anno/All.sort.csv"
ALL_SPE_PROID_RANK_PATH = "/data/Ref/Meta_anno/AllSpeProid_rank.txt"
WORMBASE_PATH = "/data/Ref/Meta_anno/wormbase.tsv"
KRONA_SCRIPT_PATH = "/data/deploy/TB_soft/other_soft/3_kreport2krona.py"


def _resource_path(name: str, fallback: str) -> str:
    runtime = get_runtime_context()
    if runtime.resources and hasattr(runtime.resources, name):
        value = getattr(runtime.resources, name)
        if value:
            return value
    return fallback


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


def exreadsID1(taxlist, kraresult, fq1, fq2=0, kkf=None):
    Maintax = taxlist[0]
    kraredb = pd.read_csv(kraresult, header=None, usecols=[1, 2], dtype={1: "str", 2: "int32"}, sep="\t")
    tmp2db = kraredb[kraredb[2].isin(taxlist)]
    tmp1db = pd.DataFrame(tmp2db[1].unique())
    tmp2db.to_csv(f"{Maintax}.id.tsv", sep="\t", index=False)
    pd.DataFrame(tmp1db).to_csv(f"{Maintax}_fqID.txt", index=False, header=False)
    run_command(f"head -n 1 {Maintax}_fqID.txt > tt.txt", check=True, logf=kkf)
    run_command(f"cut -d '/' -f1 {Maintax}_fqID.txt|sort -u > {Maintax}.listID.txt", check=True, logf=kkf)
    if os.popen("head -n 1 tt.txt").read().strip().endswith("/1") or os.popen("head -n 1 tt.txt").read().strip().endswith("/2"):
        run_command(f"sed 's/$/\\/1/' {Maintax}.listID.txt > {Maintax}.listID1.txt", check=True, logf=kkf)
        run_command(f"sed 's/$/\\/2/' {Maintax}.listID.txt > {Maintax}.listID2.txt", check=True, logf=kkf)
        run_command(f"seqkit grep -f {Maintax}.listID1.txt {fq1} > {Maintax}.1.fastq", check=True, logf=kkf)
        if fq2:
            run_command(f"seqkit grep -f {Maintax}.listID2.txt {fq2} > {Maintax}.2.fastq", check=True, logf=kkf)
    else:
        run_command(f"seqkit grep -f {Maintax}.listID.txt {fq1} > {Maintax}.1.fastq", check=True, logf=kkf)
        if fq2:
            run_command(f"seqkit grep -f {Maintax}.listID.txt {fq2} > {Maintax}.2.fastq", check=True, logf=kkf)


def getinfo(Pre, read_prefix="2", threads=10, kkf=None):
    refdict = {
        "card": _resource_path("card_sequences", conda_env_path("meta_main", "db", "card", "sequences")),
        "vfdb": _resource_path("vfdb_fasta", "/data/deploy/meta_genome/database/vfdb.fasta"),
    }
    read1 = f"{read_prefix}.1.fastq"
    read2 = f"{read_prefix}.2.fastq"
    for db in ["card", "vfdb"]:
        ref = refdict.get(db)
        if not os.path.isfile(f"2.{db}.sorted.bam") or os.popen(f"samtools view 2.{db}.sorted.bam|wc -l").read().strip() == "0":
            if os.path.isfile(read2):
                run_command(f"minimap2 -ax sr {ref} {read1} {read2} -t 10 |samtools sort -o 2.{db}.sorted.bam", check=True, logf=kkf)
            else:
                run_command(f"minimap2 -ax sr {ref} {read1} -t 10 |samtools sort -o 2.{db}.sorted.bam", check=True, logf=kkf)
            run_command(f"samtools index 2.{db}.sorted.bam", check=True, logf=kkf)
            run_command(f"mosdepth -b1 {db} 2.{db}.sorted.bam -n -t {threads}", check=True, logf=kkf)
            run_command(f"gunzip {db}.regions.bed.gz -f", check=True, logf=kkf)
            run_command(f"samtools idxstat 2.{db}.sorted.bam > {db}.stat.tsv", check=True, logf=kkf)
        dbfile = pd.read_table(f"{db}.regions.bed", header=None)
        coninfo = pd.read_table(f"{db}.stat.tsv", header=None, usecols=[0, 2])
        coninfo.columns = [0, "card_subreads"]
        depdb = pd.DataFrame(dbfile.groupby(0).apply(lambda x: round(sum(x[3]) / x.shape[0], 2)).reset_index(name="card_dep"))
        covdb = pd.DataFrame(dbfile.groupby(0).apply(lambda x: round(sum(x[3] > 0) / x.shape[0], 2)).reset_index(name="card_cov"))
        rawdb = depdb.merge(covdb, on=0).merge(coninfo, on=0).sort_values("card_subreads", ascending=False)
        if db == "card":
            rawdb = rawdb.loc[(rawdb["card_cov"] >= 0.1) & (rawdb["card_subreads"] > 10), :]
            metadb = pd.read_table(_resource_path("aro_index", "/data/deploy/meta_genome/database/aro_index.tsv"), usecols=["Model Name", "AMR Gene Family", "Drug Class", "Resistance Mechanism"])
            rawdb["Model Name"] = rawdb[0].str.split("~~~").str[1]
            rawdb = rawdb.merge(metadb, on="Model Name")
            rawdb.columns = ["片段名称", "平均深度", "覆盖率", "支持序列数", "Model", "基因家族", "耐药分类", "耐药机制"]
            rawdb["耐药基因"] = rawdb["片段名称"].str.split("~~~").str[1]
            rawdb[["耐药基因", "平均深度", "覆盖率", "支持序列数", "基因家族", "耐药分类", "耐药机制"]].to_csv("2.card.tsv", sep="\t", index=False)
        else:
            rawdb = rawdb.loc[(rawdb["card_cov"] >= 0.01) & (rawdb["card_subreads"] > 10), :].head(50)
            metadb = pd.read_table(_resource_path("vfdb_meta_annotation", "/data/deploy/meta_genome/database/VFs_meta.tsv"), encoding="Windows-1252", usecols=["VFID", "Bacteria", "Function", "Mechanism"])
            contigdb = pd.read_table(_resource_path("vfdb_contig", "/data/deploy/meta_genome/database/vfdb.contig.tsv"))
            rawdb = rawdb.merge(contigdb, left_on=0, right_on="Contig Name")
            rawdb = rawdb.merge(metadb, on="VFID")
            rawdb["毒力基因"] = rawdb[0].str.split("~~~").str[1]
            rawdb.columns = ["片段名称", "平均深度", "覆盖率", "支持序列数", "Contig", "VFID", "菌株", "毒力功能", "毒力机制", "毒力基因"]
            rawdb[["毒力基因", "平均深度", "覆盖率", "支持序列数", "VFID", "菌株", "毒力功能", "毒力机制"]].to_csv("2.vfdb.tsv", sep="\t", index=False)


def getCovDep(Pre, Pre2, kkf=None):
    runtime = get_runtime_context()
    maintax = "10239" if runtime.analysis_target == "virus" else "2"
    leveltm = 'R1' if runtime.analysis_target == "virus" else "D"
    if not os.path.isfile(f"{maintax}.1.fastq"):
        taxlist1 = proc_kra1(f"{Pre2}.report.txt", maintax, leveltm)
        taxlist1 = [int(i) for i in taxlist1]
        if os.path.isfile(f"{Pre}.R2.fastq.gz") and os.path.getsize(f"{Pre}.R2.fastq.gz") != 0:
            exreadsID1(taxlist1, f"{Pre2}.out.txt", f"{Pre}.R1.fastq.gz", f"{Pre}.R2.fastq.gz", kkf=kkf)
        else:
            exreadsID1(taxlist1, f"{Pre2}.out.txt", f"{Pre}.R1.fastq.gz", 0, kkf=kkf)
    if runtime.analysis_target != "virus":
        getinfo(Pre, read_prefix=maintax, threads=10, kkf=kkf)


def run_bracken_sub(report_path, prefix, krdb, kkf):
    testbrkdb = pd.read_table(report_path, header=None)
    if "S4" in testbrkdb[3]:
        level = "S3"
    elif "S3" in testbrkdb[3]:
        level = "S2"
    else:
        level = "S1"
    run_command(f"bracken -d {krdb} -o {prefix}_Sub.bracken1.txt -w {prefix}_Sub.bracken2.txt -l {level} -t 10  -i {report_path}", logf=kkf)


def _read_table_if_exists(path: str | Path, **kwargs) -> pd.DataFrame:
    path = Path(path)
    if not path.is_file():
        return pd.DataFrame()
    read_kwargs = dict(kwargs)
    if "sep" not in read_kwargs and "delimiter" not in read_kwargs:
        read_kwargs["sep"] = "," if path.suffix.lower() == ".csv" else "\t"
    return pd.read_csv(path, **read_kwargs)


def _write_empty_table(path: str | Path, columns: list[str]) -> None:
    pd.DataFrame(columns=columns).to_csv(path, sep="\t", index=False)


def _ontology_outputs_complete(inf, fq1, fq2, Pre: str) -> bool:
    runtime = get_runtime_context()
    use_bracken = runtime.analysis_target != "virus"
    required_files = [Path(f"{Pre}.taxonomy_summary.tsv"), Path("Summary_kraken.csv"), Path("Summary_kraken1.csv"), Path(f"{Pre}.anno.tsv")]
    if inf:
        required_files.extend([Path(f"{Pre}.report.txt"), Path(f"{Pre}.out.txt"), Path(f"{Pre}.list.txt"), Path(f"{Pre}.list2.txt")])
        if use_bracken:
            required_files.extend([Path(f"{Pre}.bracken1.txt"), Path(f"{Pre}.bracken2.txt")])
    if fq1:
        required_files.extend([Path(f"{Pre}_2.report.txt"), Path(f"{Pre}_2.out.txt"), Path(f"{Pre}_2.list.txt"), Path(f"{Pre}_2.list2.txt")])
        if use_bracken:
            required_files.extend([Path(f"{Pre}_2.bracken1.txt"), Path(f"{Pre}_2.bracken2.txt")])
            if fq2:
                required_files.extend([Path(f"{Pre}_2_Sub.bracken1.txt"), Path(f"{Pre}_2_Sub.bracken2.txt")])
    return all(path.is_file() for path in required_files)


def _kran_summ(output_species: str, output_subspecies: str, bracken_report: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    sanofile = _read_table_if_exists(_resource_path("taxa_info", TAXA_INFO_PATH))
    if not sanofile.empty and "taxid" in sanofile.columns:
        sanofile["taxid"] = pd.to_numeric(sanofile["taxid"], errors="coerce").astype("Int64")

    kran_dic: dict[str, dict] = {}
    kran_dic1: dict[str, dict] = {}
    tad = tap = tac = tao = taf = tag = "-"

    with open(bracken_report, "r", encoding="utf-8", errors="ignore") as f:
        for raw_line in f:
            line = raw_line.strip().split("\t")
            if len(line) < 6:
                continue
            line[5] = line[5].strip()
            prop = line[0].strip()
            snum = line[1].strip()
            taxid = line[4].strip()
            level = line[3]
            if line[5] == "unclassified":
                continue
            if level == "D":
                tad = line[5]
            elif level == "P":
                tap = line[5]
            elif level == "C":
                tac = line[5]
            elif level == "O":
                tao = line[5]
            elif level == "F":
                taf = line[5]
            elif level == "G":
                tag = line[5]
            elif level == "S" and int(float(snum)) > 2:
                kran_dic[line[5]] = {"D": tad, "P": tap, "C": tac, "O": tao, "F": taf, "G": tag, "S": line[5], "比例": prop, "序列数量": snum, "taxid": taxid}
            elif level in ["S1", "S2", "S3"] and int(float(snum)) > 2:
                kran_dic1[line[5]] = {"D": tad, "P": tap, "C": tac, "O": tao, "F": taf, "G": tag, "亚种": line[5], "比例": prop, "序列数量": snum, "taxid": taxid}

    species_df = pd.DataFrame(kran_dic).T if kran_dic else pd.DataFrame(columns=["D", "P", "C", "O", "F", "G", "S", "比例", "序列数量", "taxid"])
    subspecies_df = pd.DataFrame(kran_dic1).T if kran_dic1 else pd.DataFrame(columns=["D", "P", "C", "O", "F", "G", "亚种", "比例", "序列数量", "taxid"])

    if not species_df.empty:
        species_df["比例"] = pd.to_numeric(species_df["比例"], errors="coerce")
        species_df["序列数量"] = pd.to_numeric(species_df["序列数量"], errors="coerce")
        species_df["taxid"] = pd.to_numeric(species_df["taxid"], errors="coerce").astype("Int64")
        species_df.sort_values("比例", inplace=True, ascending=False)
        if not sanofile.empty and "taxid" in sanofile.columns:
            species_df = species_df.merge(sanofile, on="taxid", how="left")

    if not subspecies_df.empty:
        subspecies_df["比例"] = pd.to_numeric(subspecies_df["比例"], errors="coerce")
        subspecies_df["序列数量"] = pd.to_numeric(subspecies_df["序列数量"], errors="coerce")
        subspecies_df["taxid"] = pd.to_numeric(subspecies_df["taxid"], errors="coerce").astype("Int64")
        subspecies_df.sort_values("比例", inplace=True, ascending=False)
        if not sanofile.empty and "taxid" in sanofile.columns:
            subspecies_df = subspecies_df.merge(sanofile, on="taxid", how="left")

    rename_map = {"D": "界", "P": "门", "C": "纲", "O": "目", "F": "科", "G": "属", "S": "种"}
    species_df.rename(columns=rename_map, inplace=True)
    subspecies_df.rename(columns=rename_map, inplace=True)

    species_df.to_csv(output_species, sep="\t", index=False)
    subspecies_df.to_csv(output_subspecies, sep="\t", index=False)
    return species_df, subspecies_df


def _safe_ratio(numerator: float, denominator: float) -> str:
    if denominator == 0:
        return "0.0%"
    return f"{round(numerator / denominator * 100, 2)}%"


def _build_summary_outputs(Pre: str) -> None:
    list_path = Path(f"{Pre}_2.list.txt") if Path(f"{Pre}_2.list.txt").is_file() else Path(f"{Pre}.list.txt")
    list2_path = Path(f"{Pre}_2.list2.txt") if Path(f"{Pre}_2.list2.txt").is_file() else Path(f"{Pre}.list2.txt")
    report_path = Path(f"{Pre}_2.report.txt") if Path(f"{Pre}_2.report.txt").is_file() else Path(f"{Pre}.report.txt")
    if not list_path.is_file() or not report_path.is_file():
        return

    tmpdb = pd.read_table(list_path)
    if tmpdb.empty:
        return

    plantdb = _read_table_if_exists(_resource_path("all_sort", ALL_SORT_PATH), usecols=["taxonId", "type"])
    if not plantdb.empty and "taxid" in tmpdb.columns:
        tmpdb = tmpdb.merge(plantdb, left_on="taxid", right_on="taxonId", how="left").fillna("-")
        if "taxonId" in tmpdb.columns:
            tmpdb.drop("taxonId", inplace=True, axis=1)
    tmpdb.to_csv(f"{Pre}.anno.tsv", sep="\t", index=False)

    rawrpdb1 = pd.read_table(report_path, header=None)
    rawrpdb1[5] = rawrpdb1[5].astype(str).str.strip()
    rawrpdb = rawrpdb1.loc[rawrpdb1[5].isin(["unclassified", "root"])]

    ureads = rawrpdb.loc[rawrpdb[5] == "unclassified", 1].tolist()[0] if "unclassified" in rawrpdb[5].tolist() else 0
    creads = rawrpdb.loc[rawrpdb[5] == "root", 1].tolist()[0] if "root" in rawrpdb[5].tolist() else 0

    Sdb = pd.read_table(list_path)
    Subdb = pd.read_table(list2_path) if list2_path.is_file() else pd.DataFrame(columns=["序列数量"])
    Prodb = _read_table_if_exists(_resource_path("all_spe_proid_rank", ALL_SPE_PROID_RANK_PATH), header=None, names=["Taxid", "Type"])
    if not Prodb.empty:
        Prodb = Prodb.loc[Prodb["Type"] == "Species", :]
        Prodb1 = rawrpdb1.loc[rawrpdb1[5].isin(Prodb["Taxid"].astype(str).tolist()), :]
    else:
        Prodb1 = pd.DataFrame()
    Wormdb = _read_table_if_exists(_resource_path("wormbase", WORMBASE_PATH))
    Wormdb1 = rawrpdb1.loc[rawrpdb1[5].isin(Wormdb["taxid"].astype(str).tolist()), :] if not Wormdb.empty and "taxid" in Wormdb.columns else pd.DataFrame()

    summarydict = {}
    summarydict1 = {}
    total_reads = ureads + creads
    summarydict["有效序列"] = rawrpdb[1].sum() if not rawrpdb.empty else total_reads
    summarydict["未识别序列"] = ureads
    summarydict["未识别序列比例"] = _safe_ratio(ureads, total_reads)
    summarydict["可识别序列"] = creads
    summarydict["可识别序列比例"] = _safe_ratio(creads, total_reads)
    summarydict["校正识别序列数(种)"] = pd.to_numeric(Sdb.get("序列数量"), errors="coerce").sum() if "序列数量" in Sdb.columns else 0
    summarydict["校正识别序列数(亚种)"] = pd.to_numeric(Subdb.get("序列数量"), errors="coerce").sum() if "序列数量" in Subdb.columns else 0
    summarydict["细菌"] = pd.to_numeric(Sdb.loc[Sdb["界"] == "Bacteria", "序列数量"], errors="coerce").sum() if "界" in Sdb.columns else 0
    summarydict["病毒"] = pd.to_numeric(Sdb.loc[Sdb["界"] == "Viruses", "序列数量"], errors="coerce").sum() if "界" in Sdb.columns else 0
    summarydict1["细菌"] = Sdb.loc[Sdb["界"] == "Bacteria", :].shape[0] if "界" in Sdb.columns else 0
    summarydict1["病毒"] = Sdb.loc[Sdb["界"] == "Viruses", :].shape[0] if "界" in Sdb.columns else 0

    fungi_rows = rawrpdb1.loc[rawrpdb1[5] == "Fungi", :]
    summarydict["真菌"] = fungi_rows[1].sum() if not fungi_rows.empty else 0
    summarydict1["真菌"] = fungi_rows.shape[0]
    archaea_rows = rawrpdb1.loc[rawrpdb1[5] == "Archaea", :]
    summarydict["古菌"] = archaea_rows[1].sum() if not archaea_rows.empty else 0
    summarydict["原生动物"] = Prodb1[2].sum() if not Prodb1.empty and 2 in Prodb1.columns else 0
    summarydict["寄生虫"] = Wormdb1[2].sum() if not Wormdb1.empty and 2 in Wormdb1.columns else 0
    summarydict1["寄生虫"] = Wormdb1.shape[0]

    if Path("summary.tsv").is_file() and Path("R1_Fastqc.tsv").is_file():
        hostdb = pd.read_table("summary.tsv")
        hostdb1 = pd.read_table("R1_Fastqc.tsv")
        try:
            hostrate = round((hostdb["num_seqs"][0] - hostdb1["总序列数"][1]) / hostdb["num_seqs"][0], 4) * 100
        except Exception:
            hostrate = 0
        summarydict1["宿主"] = hostrate

    pd.DataFrame(summarydict, index=[0]).to_csv("Summary_kraken.csv")
    pd.DataFrame(summarydict1, index=[0]).to_csv("Summary_kraken1.csv")


def _generate_taxonomy_outputs(Pre: str, kkf) -> None:
    runtime = get_runtime_context()
    use_bracken = runtime.analysis_target != "virus"
    krona_script_path = _resource_path("krona_script", KRONA_SCRIPT_PATH)
    if use_bracken and Path(f"{Pre}.bracken2.txt").is_file():
        _kran_summ(f"{Pre}.list.txt", f"{Pre}.list2.txt", f"{Pre}.bracken2.txt")
        if Path(krona_script_path).is_file():
            run_command(f"{krona_script_path} -r {Pre}.bracken2.txt -o {Pre}.krona.txt", logf=kkf)
            run_command(f"ktImportText {Pre}.krona.txt -o {Pre}.krona.html", logf=kkf)
    elif Path(f"{Pre}.report.txt").is_file():
        _kran_summ(f"{Pre}.list.txt", f"{Pre}.list2.txt", f"{Pre}.report.txt")
        if Path(krona_script_path).is_file():
            run_command(f"{krona_script_path} -r {Pre}.report.txt -o {Pre}.krona.txt", logf=kkf)
            run_command(f"ktImportText {Pre}.krona.txt -o {Pre}.krona.html", logf=kkf)
    if use_bracken and Path(f"{Pre}_2.bracken2.txt").is_file():
        _kran_summ(f"{Pre}_2.list.txt", f"{Pre}_2.list2.txt", f"{Pre}_2.bracken2.txt")
        if Path(f"{Pre}_2_Sub.bracken2.txt").is_file():
            _kran_summ(f"{Pre}_2_sub.list.txt", f"{Pre}_2.list2.txt", f"{Pre}_2_Sub.bracken2.txt")
        else:
            _write_empty_table(f"{Pre}_2.list2.txt", ["界", "门", "纲", "目", "科", "属", "亚种", "比例", "序列数量", "taxid"])
        if Path(krona_script_path).is_file():
            run_command(f"{krona_script_path} -r {Pre}_2.bracken2.txt -o {Pre}_2.krona.txt", logf=kkf)
            run_command(f"ktImportText {Pre}_2.krona.txt -o {Pre}_2.krona.html", logf=kkf)
    elif Path(f"{Pre}_2.report.txt").is_file():
        _kran_summ(f"{Pre}_2.list.txt", f"{Pre}_2.list2.txt", f"{Pre}_2.report.txt")
        _write_empty_table(f"{Pre}_2_sub.list.txt", ["界", "门", "纲", "目", "科", "属", "亚种", "比例", "序列数量", "taxid"])
        if Path(krona_script_path).is_file():
            run_command(f"{krona_script_path} -r {Pre}_2.report.txt -o {Pre}_2.krona.txt", logf=kkf)
            run_command(f"ktImportText {Pre}_2.krona.txt -o {Pre}_2.krona.html", logf=kkf)
    _build_summary_outputs(Pre)


def kk2(inf, fq1, fq2, threads, Pre):
    runtime = get_runtime_context()
    krdb = runtime.krdb
    use_bracken = runtime.analysis_target != "virus"
    with open("kk2.log", "w") as kkf:
        if _ontology_outputs_complete(inf, fq1, fq2, Pre):
            kkf.write(f"[skip] {Pre} 物种鉴定结果已存在，跳过 kraken2/bracken 重复执行\n")
            open("kk2_ok", "w").write("")
            return
        if inf:
            print('?????')
            if not os.path.isfile(f"{Pre}.list.txt"):
                if not os.path.isfile(f"{Pre}.report.txt"):
                    run_command(f"kraken2 --db {krdb} --threads {threads} --output {Pre}.out.txt --report {Pre}.report.txt {inf}", logf=kkf)
                if use_bracken:
                    run_command(f"bracken -d {krdb} -o {Pre}.bracken1.txt -w {Pre}.bracken2.txt -l S -t 10  -i {Pre}.report.txt", logf=kkf)
                    run_command(f"bracken -d {krdb} -o {Pre}.bracken1.txt -w {Pre}.bracken2.txt -l S3 -t 10  -i {Pre}.report.txt", logf=kkf)
            tmpfile = pd.read_table(f"{Pre}.bracken1.txt") if use_bracken else pd.read_table(f"{Pre}.report.txt", header=None, names=["percent", "reads_clade", "reads_taxon", "taxonomy_lvl", "taxonomy_id", "name"])
            ONTSpe = tmpfile.name.tolist()[0]
            try:
                getCovDep(Pre, Pre, kkf=kkf)
            except Exception:
                pass

        if fq1 and fq2:
            print('?????')
            if not os.path.isfile(f"{Pre}_2.list.txt") and not os.path.isfile(f"{Pre}_2.report.txt"):
                run_command(f"kraken2 --db {krdb} --threads {threads} --output {Pre}_2.out.txt --report {Pre}_2.report.txt {fq1} {fq2}", logf=kkf)
                if use_bracken:
                    run_command(f"bracken -d {krdb} -o {Pre}_2.bracken1.txt -w {Pre}_2.bracken2.txt -l S -t 10  -i {Pre}_2.report.txt", logf=kkf)
                    run_bracken_sub(f"{Pre}_2.report.txt", f"{Pre}_2", krdb, kkf)
            tmpfile2 = pd.read_table(f"{Pre}_2.bracken1.txt") if use_bracken else pd.read_table(f"{Pre}_2.report.txt", header=None, names=["percent", "reads_clade", "reads_taxon", "taxonomy_lvl", "taxonomy_id", "name"])
            ngsSpe = tmpfile2.name.tolist()[0]
            getCovDep(Pre, f"{Pre}_2", kkf=kkf)
        elif fq1:
            print('?????')
            if not os.path.isfile(f"{Pre}_2.list.txt") and not os.path.isfile(f"{Pre}_2.report.txt"):
                run_command(f"kraken2 --db {krdb} --threads {threads} --output {Pre}_2.out.txt --report {Pre}_2.report.txt {fq1}", logf=kkf)
                if use_bracken:
                    run_command(f"bracken -d {krdb} -o {Pre}_2.bracken1.txt -w {Pre}_2.bracken2.txt -l S -t 10  -i {Pre}_2.report.txt", logf=kkf)
                    run_bracken_sub(f"{Pre}_2.report.txt", f"{Pre}_2", krdb, kkf)
            tmpfile2 = pd.read_table(f"{Pre}_2.bracken1.txt") if use_bracken else pd.read_table(f"{Pre}_2.report.txt", header=None, names=["percent", "reads_clade", "reads_taxon", "taxonomy_lvl", "taxonomy_id", "name"])
            ngsSpe = tmpfile2.name.tolist()[0]
            try:
                getCovDep(Pre, f"{Pre}_2", kkf=kkf)
            except Exception:
                pass

        if "ngsSpe" in dir() and "ONTSpe" in dir() and ngsSpe != ONTSpe:
            raise Exception("二三代不是同一菌种测序数据")

        if use_bracken:
            taxonomy_df = tmpfile[["name", "taxonomy_id", "taxonomy_lvl", "new_est_reads", "fraction_total_reads"]] if "ONTSpe" in dir() else tmpfile2[["name", "taxonomy_id", "taxonomy_lvl", "new_est_reads", "fraction_total_reads"]]
        else:
            source_df = tmpfile if "ONTSpe" in dir() else tmpfile2
            taxonomy_df = source_df[["name", "taxonomy_id", "taxonomy_lvl", "reads_clade", "percent"]].copy()
            taxonomy_df.rename(columns={"reads_clade": "new_est_reads", "percent": "fraction_total_reads"}, inplace=True)
        taxonomy_df = taxonomy_df.copy()
        taxonomy_df.rename(columns={"name": "物种", "taxonomy_id": "taxid", "taxonomy_lvl": "水平", "new_est_reads": "序列数量", "fraction_total_reads": "相对丰度"}, inplace=True)
        taxonomy_df.to_csv(f"{Pre}.taxonomy_summary.tsv", sep="\t", index=False)
        _generate_taxonomy_outputs(Pre, kkf)
        open("kk2_ok", "w").write("")
