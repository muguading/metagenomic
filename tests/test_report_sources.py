from __future__ import annotations

import gzip
from pathlib import Path

from bac_analysis_portal.report_payload import _build_report_payload
from bac_analysis_portal.report_sources import _read_multi_sample_depth_coverage, _read_multi_sample_mean_depth, _resolve_report_source


def test_multi_sample_report_source_prefers_ready_sample_directory(tmp_path) -> None:
    fastq_root = tmp_path / "outputs" / "fastq_analysis"
    empty_sample = fastq_root / "S1"
    ready_sample = fastq_root / "S2"
    empty_sample.mkdir(parents=True)
    ready_sample.mkdir()
    (ready_sample / "summary.tsv").write_text("sum_len\n100\n", encoding="utf-8")

    source = _resolve_report_source({
        "id": "task-1",
        "name": "batch",
        "status": "RUNNING",
        "params": {"output_dir": str(tmp_path / "outputs")},
    })

    assert source["available"] is True
    assert source["mode"] == "multi"
    assert source["selected_sample"] == "S2"
    assert source["report_dir"] == ready_sample.resolve()


def test_multi_sample_report_source_reports_pending_when_no_sample_has_results(tmp_path) -> None:
    fastq_root = tmp_path / "outputs" / "fastq_analysis"
    (fastq_root / "S1").mkdir(parents=True)
    (fastq_root / "S2").mkdir()

    source = _resolve_report_source({
        "id": "task-1",
        "name": "batch",
        "status": "RUNNING",
        "params": {"output_dir": str(tmp_path / "outputs")},
    })

    assert source["available"] is False
    assert source["mode"] == "multi"
    assert source["samples"] == ["S1", "S2"]
    assert "还没有样本产出可展示结果文件" in source["reason"]


def test_multi_sample_report_source_skips_unready_preferred_sample(tmp_path) -> None:
    fastq_root = tmp_path / "outputs" / "fastq_analysis"
    unready_preferred = fastq_root / "S1"
    ready_sample = fastq_root / "SAMPLE01"
    unready_preferred.mkdir(parents=True)
    ready_sample.mkdir()
    (ready_sample / "nextclade_output").mkdir()
    (ready_sample / "nextclade_output" / "nextclade.tsv").write_text("seqName\tclade\nSAMPLE01\t20A\n", encoding="utf-8")

    source = _resolve_report_source({
        "id": "task-1",
        "name": "batch",
        "status": "RUNNING",
        "params": {"output_dir": str(tmp_path / "outputs"), "sample_name": "S1"},
    })

    assert source["available"] is True
    assert source["selected_sample"] == "SAMPLE01"
    assert source["report_dir"] == ready_sample.resolve()


def test_multi_sample_report_source_honors_explicit_unready_sample(tmp_path) -> None:
    fastq_root = tmp_path / "outputs" / "fastq_analysis"
    explicit_sample = fastq_root / "S1"
    ready_sample = fastq_root / "SAMPLE01"
    explicit_sample.mkdir(parents=True)
    ready_sample.mkdir()
    (ready_sample / "summary.tsv").write_text("sum_len\n100\n", encoding="utf-8")

    source = _resolve_report_source(
        {
            "id": "task-1",
            "name": "batch",
            "status": "RUNNING",
            "params": {"output_dir": str(tmp_path / "outputs"), "sample_name": "SAMPLE01"},
        },
        selected_sample="S1",
    )

    assert source["available"] is True
    assert source["selected_sample"] == "S1"
    assert source["report_dir"] == explicit_sample.resolve()


def test_ncov_multi_sample_demo_skips_pending_sample_and_builds_payload(tmp_path) -> None:
    project_root = Path(__file__).resolve().parent.parent
    ncov_demo = project_root / "demo_data" / "ncov"
    fastq_root = tmp_path / "ncov_multi" / "fastq_analysis"
    pending_sample = fastq_root / "000_pending"
    ready_sample = fastq_root / "SAMPLE01"
    pending_sample.mkdir(parents=True)
    ready_sample.symlink_to(ncov_demo, target_is_directory=True)

    task = {
        "id": "local-ncov-multi",
        "name": "local-ncov-multi",
        "status": "RUNNING",
        "owner": "tester",
        "owner_group": "",
        "created_at": "",
        "started_at": "",
        "finished_at": "",
        "params": {
            "output_dir": str(tmp_path / "ncov_multi"),
            "sample_name": "000_pending",
            "analysis_target": "virus",
            "species": "SARS-CoV-2",
            "workstation_key": "virus",
            "method": "",
            "asm_type": "",
        },
    }

    source = _resolve_report_source(task)
    assert source["available"] is True
    assert source["selected_sample"] == "SAMPLE01"

    payload = _build_report_payload(task)
    assert payload["task"]["report_mode"] == "multi"
    assert payload["task"]["sample_name"] == "SAMPLE01"
    assert payload["task"]["multi_sample_summary"]["ready_count"] == 1
    assert payload["task"]["multi_sample_summary"]["table"]["rows"][1]["ready"] is True
    assert payload["task"]["multi_sample_summary"]["table"]["rows"][1]["species_name"] == "SARS-CoV-2"
    assert payload["task"]["multi_sample_summary"]["table"]["rows"][1]["coverage_1x"] == "100.00%"
    assert payload["task"]["multi_sample_summary"]["table"]["rows"][1]["coverage_10x"] == "100.00%"
    assert payload["task"]["multi_sample_summary"]["table"]["rows"][1]["coverage_100x"] == "100.00%"
    assert payload["task"]["multi_sample_summary"]["table"]["rows"][1]["mean_depth"] == "2524.40"
    assert payload["sections"]["assembly"]["coverage"]["status"] == "ready"
    assert payload["sections"]["serotype"]["status"] == "ready"


def test_other_virus_multi_sample_demos_build_batch_overview(tmp_path) -> None:
    project_root = Path(__file__).resolve().parent.parent
    cases = [
        ("rsv_demo", "Respiratory syncytial virus B", "B.D.1", "100.00%", "99.99%", "99.33%"),
        ("hiv_demo", "HIV-1", "HIV-1 / B", "96.71%", "95.90%", "94.38%"),
        ("hmpv_demo", "Human metapneumovirus", "A2.2.2", "93.43%", "93.08%", "87.30%"),
        ("denv_demo", "Dengue virus 4", "4II_B.1.3", "98.75%", "97.99%", "86.83%"),
    ]

    for demo_name, species, typing, coverage_1x, coverage_10x, coverage_100x in cases:
        fastq_root = tmp_path / demo_name / "outputs" / "fastq_analysis"
        pending_sample = fastq_root / "000_pending"
        ready_sample = fastq_root / "SAMPLE01"
        pending_sample.mkdir(parents=True)
        ready_sample.symlink_to(project_root / "demo_data" / demo_name, target_is_directory=True)

        task = {
            "id": f"local-{demo_name}-multi",
            "name": f"local-{demo_name}-multi",
            "status": "RUNNING",
            "owner": "tester",
            "owner_group": "",
            "created_at": "",
            "started_at": "",
            "finished_at": "",
            "params": {
                "output_dir": str(tmp_path / demo_name / "outputs"),
                "sample_name": "000_pending",
                "analysis_target": "virus",
                "species": species,
                "workstation_key": "virus",
                "method": "",
                "asm_type": "",
            },
        }

        source = _resolve_report_source(task)
        assert source["available"] is True
        assert source["selected_sample"] == "SAMPLE01"

        payload = _build_report_payload(task)
        summary = payload["task"]["multi_sample_summary"]
        ready_row = summary["table"]["rows"][1]
        assert payload["task"]["report_mode"] == "multi"
        assert summary["sample_count"] == 2
        assert summary["ready_count"] == 1
        assert ready_row["ready"] is True
        assert ready_row["species_name"] == species
        assert ready_row["typing"] == typing
        assert ready_row["coverage_1x"] == coverage_1x
        assert ready_row["coverage_10x"] == coverage_10x
        assert ready_row["coverage_100x"] == coverage_100x


def test_multi_sample_mean_depth_falls_back_to_per_base_bed_gz(tmp_path) -> None:
    report_dir = tmp_path / "SAMPLE01"
    report_dir.mkdir()
    with gzip.open(report_dir / "ref_map.per-base.bed.gz", "wt", encoding="utf-8") as handle:
        handle.write("NC_045512.2\t0\t2\t10\n")
        handle.write("NC_045512.2\t2\t5\t20\n")

    assert _read_multi_sample_mean_depth(report_dir, "SAMPLE01") == "16.00"
    assert _read_multi_sample_depth_coverage(report_dir, "SAMPLE01") == {
        "coverage_1x": "100.00%",
        "coverage_10x": "100.00%",
        "coverage_100x": "0.00%",
    }


def test_multi_sample_depth_coverage_falls_back_to_mosdepth_dist(tmp_path) -> None:
    report_dir = tmp_path / "SAMPLE01"
    report_dir.mkdir()
    (report_dir / "ref_map.mosdepth.global.dist.txt").write_text(
        "total\t100\t0.25\n"
        "total\t10\t0.80\n"
        "total\t1\t0.95\n",
        encoding="utf-8",
    )

    assert _read_multi_sample_depth_coverage(report_dir, "SAMPLE01") == {
        "coverage_1x": "95.00%",
        "coverage_10x": "80.00%",
        "coverage_100x": "25.00%",
    }
