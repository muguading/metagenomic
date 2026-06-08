from __future__ import annotations

import sys
import types

sys.modules.setdefault("pytaxonkit", types.SimpleNamespace())

from metagenomic_refactor.context import RuntimeContext, set_runtime_context
from metagenomic_refactor.workflow import get_flow_list, is_assembly_output_ready


def test_get_flow_list_expands_legacy_mlst_and_serotype_flow() -> None:
    set_runtime_context(
        RuntimeContext(
            ofn=".",
            runflow="功能注释,mlst与血清型,耐药与毒力",
            method="spades,freebayes",
            rmhost="norm",
            tspeabun="1",
        )
    )

    assert get_flow_list() == ["功能注释", "mlst检验", "血清型检验", "耐药与毒力"]


def test_get_flow_list_keeps_new_flow_names() -> None:
    set_runtime_context(
        RuntimeContext(
            ofn=".",
            runflow="mlst检验,血清型检验",
            method="spades,freebayes",
            rmhost="norm",
            tspeabun="1",
        )
    )

    assert get_flow_list() == ["mlst检验", "血清型检验"]


def test_assembly_output_ready_requires_nonempty_final_fasta(tmp_path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    (tmp_path / "sample.final.fasta").write_text("")

    assert not is_assembly_output_ready("sample", "spades")

    (tmp_path / "sample.final.fasta").write_text(">contig1\nACGT\n")

    assert is_assembly_output_ready("sample", "spades")
