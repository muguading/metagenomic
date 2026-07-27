from __future__ import annotations

import subprocess

from metagenomic_refactor.mag_binning import MagBinningError
from metagenomic_refactor.runner import RunnerConfig, _run_list_mode


def test_list_mode_continues_after_short_read_sample_failure(
    monkeypatch, tmp_path
) -> None:
    analysis_root = tmp_path / "fastq_analysis"
    analysis_root.mkdir()
    batch_input = tmp_path / "samples.tsv"
    batch_input.write_text(
        "样本名称\t三代数据\t二代数据左\t二代数据右\t物种信息\n"
        "failed_sample\t\tfailed_R1.fastq.gz\tfailed_R2.fastq.gz\t\n"
        "next_sample\t\tnext_R1.fastq.gz\tnext_R2.fastq.gz\t\n",
        encoding="utf-8",
    )
    cfg = RunnerConfig(
        raw_input=str(batch_input),
        inf=str(batch_input),
        intype="fastq",
        ofn=str(tmp_path),
        barkit="none",
        tmpfake=0,
        fastq1=0,
        fastq2=0,
        nt=4,
        llid="",
        mmethod="meta",
        minl="0",
        minQ="0",
        asm_type="shortasm",
        ptimes="1",
        psoft="none",
        rnalib="none",
        ref="noref",
        gtf="nogtf",
    )
    processed_samples: list[str] = []

    def fake_main_process(*args) -> None:
        sample = args[6]
        processed_samples.append(sample)
        if sample == "failed_sample":
            raise MagBinningError("mag_binning failed")

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr("metagenomic_refactor.runner._set_wkdir", lambda path: None)
    monkeypatch.setattr(
        "metagenomic_refactor.runner.subprocess.run",
        lambda *args, **kwargs: subprocess.CompletedProcess(args=args, returncode=0),
    )

    _run_list_mode(cfg, fake_main_process)

    assert processed_samples == ["failed_sample", "next_sample"]
