from __future__ import annotations

import csv
from pathlib import Path

import pytest

from metagenomic_refactor.genomad_mge import (
    GeNomadMGEError,
    _find_summary_table,
    _genomad_prefix,
    _resolve_annoele_fasta_path,
    _resolve_mobileog_query_faa,
    _tool_prefix,
    aggregate_genomad_results,
    build_parser,
    integrate_mge_tables,
    GeNomadConfig,
    load_manifest,
    load_samples,
    summarize_mge_risk,
)
from metagenomic_refactor.context import RuntimeContext, set_runtime_context


def _write_fasta(path: Path, name: str = "contig1") -> Path:
    path.write_text(f">{name}\nACGTACGT\n", encoding="utf-8")
    return path


def test_load_manifest_reads_required_columns(tmp_path: Path) -> None:
    fasta = _write_fasta(tmp_path / "sample1.fa")
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "sample\tfasta\n"
        f"sample A\t{fasta}\n",
        encoding="utf-8",
    )

    samples = load_manifest(manifest)

    assert len(samples) == 1
    assert samples[0].sample == "sample_A"
    assert samples[0].fasta == fasta.resolve()


def test_load_samples_requires_either_manifest_or_single_sample(tmp_path: Path) -> None:
    with pytest.raises(GeNomadMGEError):
        load_samples(manifest=None, sample=None, fasta=None)


def test_find_summary_table_prefers_matching_file(tmp_path: Path) -> None:
    sample_dir = tmp_path / "genomad"
    summary_dir = sample_dir / "sample_summary"
    summary_dir.mkdir(parents=True)
    plasmid = summary_dir / "sample_plasmid_summary.tsv"
    plasmid.write_text("seq_name\tlength\ncontig1\t1000\n", encoding="utf-8")

    found = _find_summary_table(sample_dir, "plasmid")

    assert found == plasmid


def test_aggregate_genomad_results_merges_plasmid_and_provirus(tmp_path: Path) -> None:
    sample_dir = tmp_path / "sample1" / "genomad"
    summary_dir = sample_dir / "sample1_summary"
    summary_dir.mkdir(parents=True)
    (summary_dir / "sample1_plasmid_summary.tsv").write_text(
        "seq_name\tlength\ttopology\tplasmid_score\ttaxonomy\n"
        "plasmid_contig\t12000\tcircular\t0.98\tPlasmid taxon\n",
        encoding="utf-8",
    )
    (summary_dir / "sample1_plasmid_genes.tsv").write_text(
        "gene\tstart\tend\n"
        "plasmid_contig_1\t20\t400\n"
        "plasmid_contig_2\t800\t2400\n",
        encoding="utf-8",
    )
    (summary_dir / "sample1_provirus_summary.tsv").write_text(
        "seq_name\tlength\tcoordinates\tprovirus_score\ttaxonomy\n"
        "host_contig\t45000\t100-5000\t0.91\tCaudoviricetes\n",
        encoding="utf-8",
    )
    outdir = tmp_path / "out"

    merged = aggregate_genomad_results(sample_dir, "sample1", outdir)

    with merged.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    assert len(rows) == 2
    assert rows[0]["mge_type"] == "plasmid"
    assert rows[0]["sequence_id"] == "plasmid_contig"
    assert rows[0]["start"] == "20"
    assert rows[0]["end"] == "2400"
    assert rows[1]["mge_type"] == "provirus"
    assert rows[1]["start"] == "100"
    assert rows[1]["end"] == "5000"


def test_build_parser_accepts_single_sample_mode() -> None:
    parser = build_parser()

    args = parser.parse_args(
        [
            "--fasta",
            "assembly.fa",
            "--sample",
            "S1",
            "--outdir",
            "outdir",
            "--database",
            "dbdir",
        ]
    )

    assert args.fasta == "assembly.fa"
    assert args.sample == "S1"


def test_genomad_prefix_uses_fixed_genomad_env() -> None:
    assert _genomad_prefix() == [
        "/home/dell/miniconda3/condabin/conda",
        "run",
        "--no-capture-output",
        "-n",
        "genomad",
    ]


def test_tool_prefix_uses_cfg_env_for_non_genomad_tools(tmp_path: Path) -> None:
    cfg = GeNomadConfig(outdir=tmp_path, database=tmp_path, conda_env="mge_tools")

    assert _tool_prefix(cfg) == [
        "/home/dell/miniconda3/condabin/conda",
        "run",
        "--no-capture-output",
        "-n",
        "mge_tools",
    ]


def test_resolve_annoele_fasta_path_uses_tmp_combine_for_meta(tmp_path: Path) -> None:
    set_runtime_context(
        RuntimeContext(
            ofn=str(tmp_path),
            runflow="元件预测",
            method="meta",
            rmhost="norm",
            tspeabun="1",
        )
    )

    assert _resolve_annoele_fasta_path("sample1", tmp_path) == tmp_path / "tmp_combine.fa"


def test_resolve_annoele_fasta_path_uses_final_fasta_for_non_meta(tmp_path: Path) -> None:
    set_runtime_context(
        RuntimeContext(
            ofn=str(tmp_path),
            runflow="元件预测",
            method="spades,freebayes",
            rmhost="norm",
            tspeabun="1",
        )
    )

    assert _resolve_annoele_fasta_path("sample1", tmp_path) == tmp_path / "sample1.final.fasta"


def test_resolve_mobileog_query_faa_uses_same_input_fasta_as_genomad(tmp_path: Path) -> None:
    input_fasta = tmp_path / "tmp_combine.fa"
    input_fasta.write_text(">contig1\nACGT\n", encoding="utf-8")

    resolved = _resolve_mobileog_query_faa("sample1", input_fasta, tmp_path)

    assert resolved == input_fasta


def test_integrate_mge_tables_merges_genomad_and_mobileog(tmp_path: Path) -> None:
    outdir = tmp_path
    (outdir / "sample1.genomad_mge_summary.tsv").write_text(
        "sample\tmge_type\tsequence_id\tlength\tstart\tend\tscore\ttopology\ttaxonomy\tsource_tsv\n"
        "sample1\tplasmid\tcontigA\t10000\t1\t10000\t0.99\tcircular\tPlasmid\tgenomad.tsv\n",
        encoding="utf-8",
    )
    (outdir / "sample1.mobileog.tsv").write_text(
        "序列名称\t参考基因组名称\t类型\t注释\t元件名称\t相似性(%)\t长度\t差异数量\t空缺数量\t序列起始\t序列终止\t参考起始\t参考终止\tevalue\t比对得分\n"
        "contigC\tMOG001|x|x|IS\tIS transposase\tIS1\t88\t300\t1\t0\t50\t349\t1\t300\t1e-20\t250\n",
        encoding="utf-8",
    )

    merged = integrate_mge_tables(outdir, "sample1")

    with merged.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    assert len(rows) == 2
    assert {row["method"] for row in rows} == {"geNomad", "mobileOG"}


def test_integrate_mge_tables_supports_nested_genomad_and_top_level_mobileog(tmp_path: Path) -> None:
    nested = tmp_path / "A01"
    nested.mkdir()
    (nested / "A01.genomad_mge_summary.tsv").write_text(
        "sample\tmge_type\tsequence_id\tlength\tstart\tend\tscore\ttopology\ttaxonomy\tsource_tsv\n"
        "A01\tplasmid\tMAG_35_46\t6670\t1\t6670\t1.0\tNo terminal repeats\tPlasmid\tgenomad.tsv\n",
        encoding="utf-8",
    )
    (tmp_path / "A01.mobileog.tsv").write_text(
        "序列名称\t参考基因组名称\t类型\t注释\t元件名称\t相似性(%)\t长度\t差异数量\t空缺数量\t序列起始\t序列终止\t参考起始\t参考终止\tevalue\t比对得分\n"
        "MAG_35_1\tMOG001|x|x|IS\tIS\tanno\tname\t88\t300\t1\t0\t50\t349\t1\t300\t1e-20\t250\n",
        encoding="utf-8",
    )

    merged = integrate_mge_tables(tmp_path, "A01")

    with merged.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    assert len(rows) == 2
    assert rows[0]["sequence_id"]
    assert {row["method"] for row in rows} == {"geNomad", "mobileOG"}


def test_summarize_mge_risk_reports_level_a_for_non_meta_boundary_overlap(tmp_path: Path) -> None:
    (tmp_path / "sample1.integrated_mge_summary.tsv").write_text(
        "sample\tmethod\tmge_type\tsequence_id\tstart\tend\tlength\tscore\tannotation\tsource_tsv\n"
        "sample1\tgeNomad\tplasmid\tcontigA\t100\t1000\t901\t1.0\tplasmid boundary\tgenomad.tsv\n"
        "sample1\tmobileOG\tcore\tcontigA\t1200\t1500\t301\t250\tintegrase\tmobileog.tsv\n",
        encoding="utf-8",
    )
    (tmp_path / "sample1.card.tsv").write_text(
        "Contig名称\t起始碱基\t终止碱基\t正负链\t基因名称\t覆盖度%\t一致性%\t产物\t耐药药物\n"
        "contigA\t200\t500\t+\tblaA\t99\t99\tbeta-lactamase\tbeta-lactam\n",
        encoding="utf-8",
    )
    (tmp_path / "sample1.vfdb.tsv").write_text(
        "Contig名称\t起始碱基\t终止碱基\t正负链\t基因名称\t覆盖度%\t一致性%\t产物\tVFID\n",
        encoding="utf-8",
    )
    prokka_dir = tmp_path / "sample1_prokka"
    prokka_dir.mkdir()
    (prokka_dir / "sample1.gff").write_text(
        "##gff-version 3\n"
        "contigA\tprokka\tCDS\t200\t500\t.\t+\t0\tID=cds1;gene=blaA;product=beta-lactamase\n"
        "contigA\tprokka\tCDS\t1200\t1500\t.\t+\t0\tID=cds2;gene=int1;product=integrase\n",
        encoding="utf-8",
    )
    set_runtime_context(
        RuntimeContext(
            ofn=str(tmp_path),
            runflow="元件预测",
            method="spades,freebayes",
            rmhost="norm",
            tspeabun="1",
        )
    )

    out = summarize_mge_risk("sample1", tmp_path, GeNomadConfig(outdir=tmp_path, database=tmp_path))

    with out.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    assert len(rows) == 1
    assert rows[0]["risk_level"] == "A"
    assert rows[0]["within_mge_boundary"] == "yes"


def test_summarize_mge_risk_uses_meta_bin_hit_sequence_names(tmp_path: Path) -> None:
    (tmp_path / "sample1.integrated_mge_summary.tsv").write_text(
        "sample\tmethod\tmge_type\tsequence_id\tstart\tend\tlength\tscore\tannotation\tsource_tsv\n"
        "sample1\tmobileOG\tcore\tMAG_1_contig_3\t400\t900\t501\t200\tintegrase\tmobileog.tsv\n",
        encoding="utf-8",
    )
    (tmp_path / "bin_card.tsv").write_text(
        "#FILE\tSEQUENCE\tSTART\tEND\tGENE\tPRODUCT\n"
        "/tmp/MAG_1.fa\tcontig_3\t100\t300\tblaA\tbeta-lactamase\n",
        encoding="utf-8",
    )
    (tmp_path / "bin_vfdb.tsv").write_text(
        "#FILE\tSEQUENCE\tSTART\tEND\tGENE\tPRODUCT\n",
        encoding="utf-8",
    )
    (tmp_path / "tmp_combine.genes.gff").write_text(
        "##gff-version 3\n"
        "MAG_1_contig_3\tProdigal\tCDS\t100\t300\t.\t+\t0\tID=orf1;gene=blaA;product=beta-lactamase\n"
        "MAG_1_contig_3\tProdigal\tCDS\t400\t900\t.\t+\t0\tID=orf2;gene=int1;product=integrase\n",
        encoding="utf-8",
    )
    set_runtime_context(
        RuntimeContext(
            ofn=str(tmp_path),
            runflow="元件预测",
            method="meta",
            rmhost="norm",
            tspeabun="1",
        )
    )

    out = summarize_mge_risk("sample1", tmp_path, GeNomadConfig(outdir=tmp_path, database=tmp_path))

    with out.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    assert len(rows) == 1
    assert rows[0]["sequence_id"] == "MAG_1_contig_3"
    assert rows[0]["nearest_core_gene"] == "core"
