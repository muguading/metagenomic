from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import Mock

from bac_analysis_portal.auto_pathosource_service import (
    _pick_auto_pathosource_candidate,
    maybe_auto_trigger_pathosource_for_meta_task,
)


def test_pick_auto_pathosource_candidate_explains_rejected_thresholds() -> None:
    candidate, reason, evaluation = _pick_auto_pathosource_candidate(
        {
            "rows": [
                {"种": "Klebsiella pneumoniae", "NCBI TaxID": "573", "比例数值": 8.0, "序列数量数值": 900},
                {"种": "Escherichia coli", "NCBI TaxID": "562", "比例数值": 3.0, "序列数量数值": 1200},
            ]
        },
        {
            "enabled": True,
            "priority_only": False,
            "min_abundance_percent": 15.0,
            "min_support_reads": 800,
            "min_coverage_percent": 20.0,
            "allowed_sample_sources": [],
            "priority_species": [],
            "excluded_species": [],
        },
        coverage_percent=35.0,
        sample_source="血液",
    )
    assert candidate is None
    assert reason == "当前物种结果尚未达到自动溯源阈值。"
    assert evaluation["checks"][1]["status"] == "failed"
    assert evaluation["evaluated_candidates"][0]["species_name"] == "Klebsiella pneumoniae"
    assert "丰度 8.00% < 15.00%" in evaluation["evaluated_candidates"][0]["reason"]


def test_auto_pathosource_respects_auto_start_disabled(tmp_path: Path) -> None:
    report_dir = tmp_path / "report"
    report_dir.mkdir()
    (report_dir / "summary.tsv").write_text("sum_len\n100\n", encoding="utf-8")
    (report_dir / "tmp_combine.fa").write_text(">sample\nATGC\n", encoding="utf-8")

    rules = {
        "enabled": True,
        "priority_only": False,
        "auto_start": False,
        "min_abundance_percent": 15.0,
        "min_support_reads": 800,
        "min_coverage_percent": 20.0,
        "allowed_sample_sources": [],
        "priority_species": [],
        "excluded_species": [],
        "max_reference_genomes": 30,
        "msa_method": "snippy",
        "tree_method": "ML",
    }
    store = Mock()
    store.get_setting.side_effect = lambda key, default="": json.dumps(rules, ensure_ascii=False) if key == "admin_pathosource_trigger_rules" else default
    store.list_sample_library.return_value = []
    task_manager = Mock()
    task = {
        "id": "task-1",
        "name": "meta-task",
        "status": "SUCCEEDED",
        "owner": "analyst",
        "params": {
            "workstation_key": "metagenome",
            "method": "meta",
            "output_dir": str(report_dir),
            "sample_source": "血液",
        },
    }
    payload = {
        "task": {"sample_name": "sample", "output_dir": str(report_dir)},
        "sections": {
            "taxonomy": {
                "species_taxonomy": {
                    "rows": [
                        {"种": "Klebsiella pneumoniae", "NCBI TaxID": "573", "比例数值": 24.0, "序列数量数值": 1800}
                    ]
                }
            },
            "coverage": {"points": [1, 1, 1, 0]},
        },
    }

    result = maybe_auto_trigger_pathosource_for_meta_task(
        task,
        payload,
        store=store,
        task_manager=task_manager,
        project_root=Path(__file__).parents[1],
        owner_fallback=lambda: "admin",
    )

    trigger = result["auto_pathosource_trigger"]
    assert trigger["status"] == "would_trigger"
    assert trigger["trigger_species"] == "Klebsiella pneumoniae"
    assert trigger["rules"]["auto_start"] is False
    assert any(item["key"] == "min_abundance_percent" and item["status"] == "passed" for item in trigger["checks"])
    task_manager.create_task.assert_not_called()
    task_manager.update_task_fields.assert_not_called()
