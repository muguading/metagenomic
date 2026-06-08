from __future__ import annotations

import csv
import subprocess
from pathlib import Path

import pytest

from metagenomic_refactor.meta_quant import (
    MetaQuantError,
    MetaQuantConfig,
    MetaQuantSample,
    _coverm_prefix,
    _predict_genes_for_bins,
    _ribodetector_prefix,
    _server_has_gpu,
    build_parser,
    compute_fpm_tpm_scatter_stats,
    load_coverm_table,
    load_manifest,
    merge_quant_tables,
    render_fpm_tpm_scatter,
)


def _write_fastq(path: Path) -> Path:
    path.write_text("@r1\nACGT\n+\n!!!!\n", encoding="utf-8")
    return path


def _write_fasta(path: Path, name: str = "seq1") -> Path:
    path.write_text(f">{name}\nACGTACGT\n", encoding="utf-8")
    return path


def test_load_manifest_reads_required_columns(tmp_path: Path) -> None:
    bins_dir = tmp_path / "bins"
    bins_dir.mkdir()
    _write_fasta(bins_dir / "bin.1.fa")
    dna1 = _write_fastq(tmp_path / "dna.R1.fastq")
    dna2 = _write_fastq(tmp_path / "dna.R2.fastq")
    rna1 = _write_fastq(tmp_path / "rna.R1.fastq")
    rna2 = _write_fastq(tmp_path / "rna.R2.fastq")
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "\t".join(["sample", "bins_dir", "dna_fastq1", "dna_fastq2", "rna_fastq1", "rna_fastq2"])
        + "\n"
        + "\t".join(["sample A", str(bins_dir), str(dna1), str(dna2), str(rna1), str(rna2)])
        + "\n",
        encoding="utf-8",
    )

    samples = load_manifest(manifest)

    assert len(samples) == 1
    assert samples[0].sample == "sample_A"
    assert samples[0].bins_dir == bins_dir.resolve()


def test_load_manifest_rejects_missing_fastq_column(tmp_path: Path) -> None:
    bins_dir = tmp_path / "bins"
    bins_dir.mkdir()
    _write_fasta(bins_dir / "bin.1.fa")
    dna1 = _write_fastq(tmp_path / "dna.R1.fastq")
    dna2 = _write_fastq(tmp_path / "dna.R2.fastq")
    rna1 = _write_fastq(tmp_path / "rna.R1.fastq")
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "\t".join(["sample", "bins_dir", "dna_fastq1", "dna_fastq2", "rna_fastq1", "rna_fastq2"])
        + "\n"
        + "\t".join(["sample1", str(bins_dir), str(dna1), str(dna2), str(rna1), ""])
        + "\n",
        encoding="utf-8",
    )

    with pytest.raises(MetaQuantError):
        load_manifest(manifest)


def test_load_coverm_table_reads_count_and_tpm(tmp_path: Path) -> None:
    table = tmp_path / "coverm.tsv"
    table.write_text(
        "Contig\tCount\tTPM\n"
        "bin1|gene1\t10\t100.5\n"
        "bin1|gene2\t20\t200.25\n",
        encoding="utf-8",
    )

    result = load_coverm_table(table)

    assert result["bin1|gene1"]["count"] == 10
    assert result["bin1|gene2"]["tpm"] == 200.25


def test_merge_quant_tables_calculates_dna_fpm_and_keeps_rna_tpm(tmp_path: Path) -> None:
    gene_meta = tmp_path / "gene_metadata.tsv"
    gene_meta.write_text(
        "gene_id\tbin_id\tsource_gene_id\tnt_length\n"
        "bin1|gene1\tbin1\tgene1\t120\n"
        "bin1|gene2\tbin1\tgene2\t240\n",
        encoding="utf-8",
    )
    dna_coverm = tmp_path / "dna.tsv"
    dna_coverm.write_text(
        "Contig\tCount\tTPM\n"
        "bin1|gene1\t30\t300\n"
        "bin1|gene2\t70\t700\n",
        encoding="utf-8",
    )
    rna_coverm = tmp_path / "rna.tsv"
    rna_coverm.write_text(
        "Contig\tCount\tTPM\n"
        "bin1|gene1\t5\t400\n"
        "bin1|gene2\t15\t600\n",
        encoding="utf-8",
    )
    output = tmp_path / "merged.tsv"

    merged_path = merge_quant_tables(gene_meta, dna_coverm, rna_coverm, output)

    with merged_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    assert len(rows) == 2
    assert rows[0]["gene_id"] == "bin1|gene1"
    assert float(rows[0]["dna_fpm"]) == 300000.0
    assert float(rows[1]["dna_fpm"]) == 700000.0
    assert float(rows[0]["rna_tpm"]) == 400.0
    assert float(rows[1]["rna_count"]) == 15.0


def test_compute_fpm_tpm_scatter_stats_returns_perfect_r2_for_linear_data(tmp_path: Path) -> None:
    gene_quant = tmp_path / "gene_quant.tsv"
    gene_quant.write_text(
        "gene_id\tbin_id\tsource_gene_id\tnt_length\tdna_count\tdna_fpm\trna_count\trna_tpm\n"
        "g1\tb1\tg1\t100\t10\t9\t10\t9\n"
        "g2\tb1\tg2\t100\t10\t99\t10\t99\n"
        "g3\tb1\tg3\t100\t10\t999\t10\t999\n",
        encoding="utf-8",
    )

    stats = compute_fpm_tpm_scatter_stats(gene_quant)

    assert stats.point_count == 3
    assert stats.r2 == 1.0


def test_render_fpm_tpm_scatter_writes_stats_file(tmp_path: Path) -> None:
    gene_quant = tmp_path / "gene_quant.tsv"
    gene_quant.write_text(
        "gene_id\tbin_id\tsource_gene_id\tnt_length\tdna_count\tdna_fpm\trna_count\trna_tpm\n"
        "g1\tb1\tg1\t100\t10\t10\t10\t20\n"
        "g2\tb1\tg2\t100\t20\t20\t20\t40\n",
        encoding="utf-8",
    )
    png_path = tmp_path / "scatter.png"
    stats_path = tmp_path / "scatter.tsv"

    try:
        render_fpm_tpm_scatter(gene_quant, png_path, stats_path)
    except MetaQuantError as exc:
        if "matplotlib" in str(exc):
            pytest.skip("matplotlib not available")
        raise

    assert png_path.is_file()
    assert stats_path.is_file()
    content = stats_path.read_text(encoding="utf-8")
    assert "r2" in content
    assert "log10(DNA FPM + 1)" in content


def test_build_parser_supports_skip_ribodetector_flag() -> None:
    parser = build_parser()

    args = parser.parse_args(
        [
            "--manifest",
            "samples.tsv",
            "--outdir",
            "outdir",
            "--skip-ribodetector",
        ]
    )

    assert args.skip_ribodetector is True


def test_server_has_gpu_returns_true_when_nvidia_smi_succeeds(monkeypatch: pytest.MonkeyPatch) -> None:
    def fake_run(*args, **kwargs):
        return subprocess.CompletedProcess(args=args[0], returncode=0, stdout="GPU 0: Test\n", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)

    assert _server_has_gpu() is True


def test_server_has_gpu_returns_false_when_nvidia_smi_missing(monkeypatch: pytest.MonkeyPatch) -> None:
    def fake_run(*args, **kwargs):
        raise OSError("nvidia-smi not found")

    monkeypatch.setattr(subprocess, "run", fake_run)

    assert _server_has_gpu() is False


def test_ribodetector_prefix_uses_requested_conda_path() -> None:
    assert _ribodetector_prefix() == [
        "/home/dell/miniconda3/condabin/conda",
        "run",
        "--no-capture-output",
        "-n",
        "Ribodetector",
    ]


def test_coverm_prefix_uses_fixed_coverm_env() -> None:
    assert _coverm_prefix() == [
        "/home/dell/miniconda3/condabin/conda",
        "run",
        "--no-capture-output",
        "-n",
        "coverm",
    ]


def test_predict_genes_for_bins_reuses_existing_predictions(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    bins_dir = tmp_path / "bins"
    bins_dir.mkdir()
    _write_fasta(bins_dir / "bin1.fa", "contig1")
    sample_out = tmp_path / "sample_out"
    pred_dir = sample_out / "gene_catalog" / "predictions" / "bin1"
    pred_dir.mkdir(parents=True)
    (pred_dir / "bin1.genes.fna").write_text(">gene1\nATGC\n", encoding="utf-8")
    (pred_dir / "bin1.proteins.faa").write_text(">gene1\nM\n", encoding="utf-8")
    (pred_dir / "bin1.genes.gff").write_text("##gff-version 3\n", encoding="utf-8")

    sample = MetaQuantSample(
        sample="sample1",
        bins_dir=bins_dir,
        dna_fastq1=tmp_path / "dna1.fq",
        dna_fastq2=tmp_path / "dna2.fq",
        rna_fastq1=tmp_path / "rna1.fq",
        rna_fastq2=tmp_path / "rna2.fq",
    )
    cfg = MetaQuantConfig(outdir=tmp_path / "out")

    def fail_run_command(*args, **kwargs):
        raise AssertionError("prodigal should not be called when prediction outputs already exist")

    monkeypatch.setattr("metagenomic_refactor.meta_quant._run_command", fail_run_command)

    combined_cds, gene_meta = _predict_genes_for_bins(sample, cfg, sample_out)

    assert combined_cds.is_file()
    assert gene_meta.is_file()
    assert "bin1|gene1" in combined_cds.read_text(encoding="utf-8")
