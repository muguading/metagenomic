from __future__ import annotations

import subprocess
from pathlib import Path

from scripts.backfill_meta_plas_vf_card import discover_tasks, run_task
from metagenomic_refactor.assembly import binning_result


def test_binning_result_writes_meta_table_when_mag_bins_are_missing(
    monkeypatch, tmp_path: Path
) -> None:
    project_root = Path(__file__).resolve().parents[1]
    sample_binning = project_root / "test_data" / "codex_mag_binning"
    contigs = tmp_path / "megahit_output" / "final.contigs.fa"
    contigs.parent.mkdir()
    contigs.write_text(
        (
            sample_binning / "MH073106" / "filtered_contigs" / "MH073106.min1500.fa"
        ).read_text(encoding="utf-8"),
        encoding="utf-8",
    )

    def fake_run(*args, **kwargs):
        return subprocess.CompletedProcess(args=args[0] if args else "", returncode=0)

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr("metagenomic_refactor.assembly.subprocess.run", fake_run)

    binning_result("MH073106")

    meta_table = tmp_path / "meta_plas_vf_card.tsv"
    assert meta_table.is_file()
    content = meta_table.read_text(encoding="utf-8")
    assert content.startswith("contig_name\tlabel\tPlasmid")
    assert "MAG_1_" in content
    assert (tmp_path / "binning_status.tsv").is_file()


def test_backfill_script_discovers_fastq_analysis_sample_and_generates_table(
    monkeypatch, tmp_path: Path, capsys
) -> None:
    project_root = Path(__file__).resolve().parents[1]
    sample_binning = project_root / "test_data" / "codex_mag_binning"
    sample_dir = tmp_path / "fastq_analysis" / "MH073106"
    filtered = sample_dir / "codex_mag_binning" / "MH073106" / "filtered_contigs"
    filtered.mkdir(parents=True)
    (filtered / "MH073106.min1500.fa").write_text(
        (
            sample_binning / "MH073106" / "filtered_contigs" / "MH073106.min1500.fa"
        ).read_text(encoding="utf-8"),
        encoding="utf-8",
    )

    def fake_run(*args, **kwargs):
        return subprocess.CompletedProcess(args=args[0] if args else "", returncode=0)

    monkeypatch.setattr("metagenomic_refactor.assembly.subprocess.run", fake_run)

    tasks = discover_tasks(tmp_path / "fastq_analysis", sample=None, force=False)

    assert len(tasks) == 1
    assert tasks[0].status == "run"
    assert (sample_dir / "megahit_output" / "final.contigs.fa").is_file()

    run_task(tasks[0])

    assert (sample_dir / "meta_plas_vf_card.tsv").is_file()
    captured = capsys.readouterr()
    assert "MH073106 [1/9] 准备 MAG/contig 输入" in captured.out
    assert "MH073106 [9/9] 合并生成 meta_plas_vf_card.tsv" in captured.out
