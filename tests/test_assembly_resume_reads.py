from __future__ import annotations

import sys
import types

sys.modules.setdefault("pytaxonkit", types.SimpleNamespace())

from metagenomic_refactor.assembly import _extracted_reads_ready


def test_extracted_reads_ready_accepts_nonempty_single_end(tmp_path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    (tmp_path / "2.1.fastq").write_text("@r1\nACGT\n+\n!!!!\n")

    assert _extracted_reads_ready("2", expect_read2=False)


def test_extracted_reads_ready_rejects_empty_primary_read(tmp_path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    (tmp_path / "2.1.fastq").write_text("")

    assert not _extracted_reads_ready("2", expect_read2=False)


def test_extracted_reads_ready_requires_nonempty_read2_for_paired_input(tmp_path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    (tmp_path / "2.1.fastq").write_text("@r1\nACGT\n+\n!!!!\n")
    (tmp_path / "2.2.fastq").write_text("")

    assert not _extracted_reads_ready("2", expect_read2=True)

    (tmp_path / "2.2.fastq").write_text("@r2\nTGCA\n+\n!!!!\n")

    assert _extracted_reads_ready("2", expect_read2=True)
