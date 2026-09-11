from pathlib import Path

from conftest import login_as


APP_JS = Path(__file__).resolve().parents[1] / "bac_analysis_portal" / "static" / "app.js"
REPORT_RUNTIME = Path(__file__).resolve().parents[1] / "bac_analysis_portal" / "static" / "report_runtime.js"


def test_influenza_submission_allows_workflow_selected_reference(release_app_client, release_data) -> None:
    _app, client = release_app_client
    headers = login_as(client)
    response = client.post(
        "/api/tasks",
        json={
            "workstation_key": "virus",
            "task_name": "influenza_auto_reference",
            "input_path": str(release_data["fasta"]),
            "inputtype": "fasta",
            "output_dir": str(release_data["output_root"]),
            "asm_type": "shortref",
            "method": "bwa",
            "species": "Influenza virus",
            "thread": 2,
        },
        headers=headers,
    )
    assert response.status_code == 201, response.get_data(as_text=True)
    assert response.get_json()["params"]["ref"] == "noref"


def test_submission_ui_makes_reference_optional_only_for_influenza() -> None:
    source = APP_JS.read_text(encoding="utf-8")
    assert "function doesVirusRequireReference()" in source
    assert "return isVirusWorkstation() && !isInfluenzaSpecies();" in source
    assert "if (doesVirusRequireReference()) {\n      required.push(\"ref\");" in source
    assert "流感参考由工作流自动选择" in source


def test_influenza_subtype_maps_to_matching_ha_na_nextclade_datasets(tmp_path, monkeypatch) -> None:
    from metagenomic_refactor import virus_analysis

    root = tmp_path / "database"
    datasets = root / "nextclade_db"
    for name in ("nextstrain_flu_h1n1_ha", "nextstrain_flu_h1n1_na"):
        (datasets / name).mkdir(parents=True)
        (datasets / name / "pathogen.json").write_text("{}", encoding="utf-8")
    monkeypatch.setattr(virus_analysis, "_database_root", lambda: root)

    assert [item.name for item in virus_analysis._resolve_influenza_segment_nextclade_datasets("H1N1", "HA")] == ["nextstrain_flu_h1n1_ha"]
    assert [item.name for item in virus_analysis._resolve_influenza_segment_nextclade_datasets("H1N1", "NA")] == ["nextstrain_flu_h1n1_na"]


def test_influenza_na_selection_aggregates_reads_across_references() -> None:
    from metagenomic_refactor.virus_analysis import _select_influenza_subtype_reference

    common = {"influenza_type": "Influenza A virus", "segment_group": "NA", "source": "na_subtype"}
    rows = [
        {**common, "reference_id": "A_NA_N1__single", "subtype": "N1", "coverage_pct": 20.58, "mean_depth": 1.59, "covered_bases": 297, "mapped_reads": 64},
        {**common, "reference_id": "A_NA_N2__first", "subtype": "N2", "coverage_pct": 20.14, "mean_depth": 0.30, "covered_bases": 284, "mapped_reads": 40},
        {**common, "reference_id": "A_NA_N2__second", "subtype": "N2", "coverage_pct": 18.65, "mean_depth": 0.20, "covered_bases": 263, "mapped_reads": 35},
    ]

    selected = _select_influenza_subtype_reference(rows, "Influenza A virus", "NA")

    assert selected is not None
    assert selected["subtype"] == "N2"
    assert selected["reference_id"] == "A_NA_N2__first"


def test_influenza_subtype_selection_prioritizes_meaningfully_higher_coverage() -> None:
    from metagenomic_refactor.virus_analysis import _select_influenza_subtype_reference

    common = {"influenza_type": "Influenza A virus", "segment_group": "HA", "source": "ha_subtype"}
    rows = [
        {**common, "reference_id": "A_HA_H1__reads", "subtype": "H1", "coverage_pct": 4.98, "mean_depth": 16.35, "covered_bases": 92, "mapped_reads": 847},
        {**common, "reference_id": "A_HA_H3__coverage", "subtype": "H3", "coverage_pct": 21.71, "mean_depth": 1.05, "covered_bases": 375, "mapped_reads": 816},
    ]

    selected = _select_influenza_subtype_reference(rows, "Influenza A virus", "HA")

    assert selected is not None
    assert selected["subtype"] == "H3"
    assert selected["reference_id"] == "A_HA_H3__coverage"


def test_influenza_pairwise_candidates_use_one_reference_per_leading_subtype() -> None:
    from metagenomic_refactor.virus_analysis import _select_influenza_pairwise_candidates

    common = {"influenza_type": "Influenza A virus", "segment_group": "NA", "source": "na_subtype"}
    rows = [
        {**common, "reference_id": "A_NA_N1__best", "subtype": "N1", "coverage_pct": 25.16, "mean_depth": 3.13, "mapped_reads": 127},
        {**common, "reference_id": "A_NA_N2__best", "subtype": "N2", "coverage_pct": 15.25, "mean_depth": 0.93, "mapped_reads": 37},
        {**common, "reference_id": "A_NA_N2__lower", "subtype": "N2", "coverage_pct": 14.54, "mean_depth": 2.34, "mapped_reads": 73},
    ]

    candidates = _select_influenza_pairwise_candidates(rows, "Influenza A virus", "NA")

    assert [row["reference_id"] for row in candidates] == ["A_NA_N1__best", "A_NA_N2__best"]


def test_influenza_report_includes_ready_nextclade_segment_results(tmp_path) -> None:
    from bac_analysis_portal.serotype_reports import _read_influenza_typing_section

    workflow_dir = tmp_path / "wf_flu"
    workflow_dir.mkdir()
    (workflow_dir / "typing_summary.tsv").write_text(
        "status\tinfluenza_type\tha_subtype\tna_subtype\tsubtype_call\treference_path\n"
        "ready\tInfluenza A virus\tH3\tN2\tH3N2\t/tmp/test1.final_segments.fa\n",
        encoding="utf-8",
    )
    nextclade_dir = workflow_dir / "nextclade"
    nextclade_dir.mkdir()
    (nextclade_dir / "segment_analysis.tsv").write_text(
        "segment\tsubtype_call\tdataset\tclade\tstatus\tdetail\tqc_json\n"
        "HA\tH3N2\tnextstrain_flu_h3n2_ha_CY163680\t3C.2a1b.2a\tready\tgood\ttest1.nextstrain_flu_h3n2_ha_CY163680.json\n"
        "HA\tH3N2\tnextstrain_flu_h3n2_ha_EPI1857216\t3C.2a1b.2a\tready\tgood\ttest1.nextstrain_flu_h3n2_ha_EPI1857216.json\n"
        "NA\tH3N2\tnextstrain_flu_h3n2_na_EPI1857215\t-\tready\tgood\ttest1.nextstrain_flu_h3n2_na_EPI1857215.json\n",
        encoding="utf-8",
    )
    for dataset in ("nextstrain_flu_h3n2_ha_CY163680", "nextstrain_flu_h3n2_na_EPI1857215"):
        (nextclade_dir / f"test1.{dataset}.json").write_text(
            '{"results": [{"qc": {"overallStatus": "good", "overallScore": 0, "missingData": {"status": "good", "score": 0}}}]}',
            encoding="utf-8",
        )

    section = _read_influenza_typing_section(tmp_path, "test1")

    assert section is not None
    assert section["nextclade_segments"]["columns"] == ["segment", "subtype_call", "dataset", "clade", "status", "detail", "qc_json"]
    assert section["nextclade_segments"]["rows"][0][3] == "3C.2a1b.2a"
    assert section["summary_cards"][1]["value"] == "H3（Nextclade：3C.2a1b.2a）"
    assert section["nextclade_qc_results"][0]["qc_result"]["overallStatus"] == "good"
    assert "已完成 2 个 HA/NA 节段的 Nextclade 分型" in section["notes"]


def test_influenza_report_runtime_renders_nextclade_only_when_results_are_ready() -> None:
    source = REPORT_RUNTIME.read_text(encoding="utf-8")

    assert 'id="influenza-typing-nextclade"' in source
    assert "流感 Nextclade 节段分型" in source
    assert "节段 Nextclade QC结果" in source
    assert "const hasNextcladeResults" in source
    assert "const hasInfluenzaNextclade" in source


def test_influenza_subtype_display_joins_distinct_nextclade_clades() -> None:
    from bac_analysis_portal.serotype_reports import _format_influenza_subtype_with_nextclade

    table = {
        "columns": ["segment", "clade", "status"],
        "rows": [["HA", "A.1", "ready"], ["HA", "B.2", "ready"], ["HA", "A.1", "ready"]],
    }

    assert _format_influenza_subtype_with_nextclade("H3", "HA", table) == "H3（Nextclade：A.1|B.2）"
