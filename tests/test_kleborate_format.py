import sys
import types
from pathlib import Path


def test_new_kleborate_output_is_normalised_to_legacy_serotype_format():
    sys.modules.setdefault("pytaxonkit", types.SimpleNamespace())
    from metagenomic_refactor.strain_typing import _normalise_kleborate_table, _read_kleborate_result_table

    raw = _read_kleborate_result_table(Path("test_data/result"))
    normalised = _normalise_kleborate_table(raw, "HPKP26004_S6_L001_001")

    assert normalised.columns.tolist() == [
        "样本名称",
        "ST",
        "毒力得分",
        "耐药得分",
        "耶尔森菌素",
        "大肠菌素",
        "氨苄类耐药SHV等位基因",
        "SHV耐药突变",
        "wzi荚膜预测",
        "KO血清型",
    ]
    assert normalised.iloc[0].to_dict() == {
        "样本名称": "HPKP26004_S6_L001_001",
        "ST": "ST23",
        "毒力得分": "5",
        "耐药得分": "0",
        "耶尔森菌素": "ybt 1; ICEKp10",
        "大肠菌素": "clb 2",
        "氨苄类耐药SHV等位基因": "SHV-11^",
        "SHV耐药突变": "35Q",
        "wzi荚膜预测": "wzi1",
        "KO血清型": "KL1|OL2α.2",
    }
