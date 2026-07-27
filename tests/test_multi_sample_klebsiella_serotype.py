from bac_analysis_portal.report_sources import _format_multi_sample_bacterial_serotype, _read_multi_sample_typing_call


def test_multi_sample_bacterial_serotype_uses_klebsiella_ko_serotype():
    row = {
        "样本名称": "HPKP26004_S6_L001_001",
        "ST": "ST23",
        "毒力得分": "5",
        "耐药得分": "0",
        "KO血清型": "KL1|OL2α.2",
    }

    assert _format_multi_sample_bacterial_serotype(row) == "KL1|OL2α.2"


def test_multi_sample_bacterial_serotype_can_build_klebsiella_ko_from_loci():
    row = {
        "klebsiella_pneumo_complex__kaptive__K_locus": "KL1",
        "klebsiella_pneumo_complex__kaptive__O_locus": "OL2α.2",
    }

    assert _format_multi_sample_bacterial_serotype(row) == "KL1|OL2α.2"


def test_multi_sample_typing_call_falls_back_to_available_serotype_file(tmp_path):
    sample_dir = tmp_path / "HPKP26006_S8_L001_001"
    sample_dir.mkdir()
    (sample_dir / "actual_prefix_serotype_result.tsv").write_text(
        "样本名称\tST\tKO血清型\n"
        "HPKP26006_S8_L001_001\tST11\tKL64|OL2α.1\n",
        encoding="utf-8",
    )

    result = _read_multi_sample_typing_call(
        sample_dir,
        "HPKP26006_S8_L001_001",
        {"species_name": "Klebsiella pneumoniae(97%)"},
        is_virus=False,
    )

    assert result["serotype"] == "KL64|OL2α.1"


def test_multi_sample_typing_call_ignores_placeholder_exact_serotype_file(tmp_path):
    sample = "HPKP26006_S8_L001_001"
    sample_dir = tmp_path / sample
    sample_dir.mkdir()
    (sample_dir / f"{sample}_serotype_result.tsv").write_text(
        "样本名称\t血清型\n"
        f"{sample}\t-\n",
        encoding="utf-8",
    )
    (sample_dir / "rerun_kleborate_serotype_result.tsv").write_text(
        "样本名称\tST\tKO血清型\n"
        f"{sample}\tST11\tKL64|OL2α.1\n",
        encoding="utf-8",
    )

    result = _read_multi_sample_typing_call(
        sample_dir,
        sample,
        {"species_name": "Klebsiella pneumoniae(97%)"},
        is_virus=False,
    )

    assert result["serotype"] == "KL64|OL2α.1"
