"""Tests for virus-specific fastp trimming options."""

import metagenomic_refactor.context as runtime_context

from metagenomic_refactor.context import RuntimeContext
from metagenomic_refactor.qc import _virus_fastp_trim_options


def _runtime(analysis_target: str) -> RuntimeContext:
    return RuntimeContext(
        ofn=".",
        runflow="",
        method="",
        rmhost="norm",
        tspeabun="1",
        analysis_target=analysis_target,
    )


def test_virus_fastp_trims_twenty_bases_from_each_end(monkeypatch) -> None:
    monkeypatch.setattr(runtime_context, "runtime", _runtime("virus"))

    assert _virus_fastp_trim_options(paired_end=True) == (
        "--trim_front1 20 --trim_tail1 20 --trim_front2 20 --trim_tail2 20"
    )
    assert _virus_fastp_trim_options(paired_end=False) == "--trim_front1 20 --trim_tail1 20"


def test_non_virus_fastp_keeps_existing_trimming_behavior(monkeypatch) -> None:
    monkeypatch.setattr(runtime_context, "runtime", _runtime("bacteria"))

    assert _virus_fastp_trim_options(paired_end=True) == ""
