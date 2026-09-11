from __future__ import annotations

from pathlib import Path
from unittest.mock import Mock

from bac_analysis_portal.modeling_service import ModelingPlatformService, _prediction_interpretation
from bac_analysis_portal.store import PortalStore


def make_store(tmp_path: Path) -> PortalStore:
    store = PortalStore(db_path=tmp_path / "portal.sqlite3", project_root=tmp_path)
    store.initialize(initial_admin_password="test-password")
    return store


def test_modeling_dataset_uses_visible_sample_snapshot(tmp_path: Path) -> None:
    store = make_store(tmp_path)
    sample_library = Mock()
    sample_library.list_visible.return_value = [
        {
            "sample_key": "s1",
            "sample_name": "S1",
            "species_name": "Klebsiella pneumoniae",
            "collection_date": "2026-01-02",
            "country": "上海",
            "sample_type": "blood",
            "sample_source": "clinical",
            "sequencing_method": "illumina",
            "mlst_st": "ST11",
            "serotype_result": "K64",
            "resistance_gene_hits": "blaKPC-2",
            "virulence_gene_hits": "iucA",
            "custom_metadata_json": "[]",
        }
    ]
    service = ModelingPlatformService(store=store, sample_library=sample_library, project_root=tmp_path)

    dataset = service.create_dataset(
        scope="main",
        role="admin",
        username="admin",
        group_name="",
        payload={"dataset_name": "训练集", "filters": {"pathogen": "Klebsiella"}},
    )

    assert dataset["sample_ids"] == ["s1"]
    assert dataset["summary"]["sample_count"] == 1
    assert dataset["filters"]["pathogen"] == "Klebsiella"


def test_modeling_feature_set_records_schema_version(tmp_path: Path) -> None:
    store = make_store(tmp_path)
    service = ModelingPlatformService(store=store, sample_library=Mock(), project_root=tmp_path)

    feature_set = service.create_feature_set(
        username="analyst",
        payload={"feature_set_name": "基础特征", "version": "2", "feature_schema": {"features": ["region", "ST", "AMR genes"]}},
    )

    assert feature_set["version"] == "2"
    assert feature_set["feature_schema"]["features"] == ["region", "ST", "AMR genes"]
    assert service.list_feature_sets()["items"][0]["feature_set_id"] == feature_set["feature_set_id"]


def test_trend_dataset_uses_historical_window_to_label_future_change(tmp_path: Path) -> None:
    store = make_store(tmp_path)
    service = ModelingPlatformService(store=store, sample_library=Mock(), project_root=tmp_path)
    rows = []
    for month in range(1, 25):
        hit_count = 2 if month <= 8 or month >= 17 else 8
        for index in range(10):
            rows.append(
                {
                    "collection_date": f"{2024 + ((month - 1) // 12)}-{((month - 1) % 12) + 1:02d}-15",
                    "serotype_result": "Enteritidis" if index < hit_count else "Typhimurium",
                    "resistance_gene_hits": "blaTEM-1B" if index < hit_count else "",
                    "virulence_gene_hits": "",
                }
            )

    x_rows, labels, timeline = service._build_trend_dataset(
        rows,
        signal_type="serotype",
        target_name="Enteritidis",
        lookback_months=3,
        forecast_months=1,
        regular_threshold=0.2,
    )

    assert len(x_rows) == len(labels) == len(timeline)
    assert "常规变化" in labels
    assert any(label in {"明显上升", "明显下降"} for label in labels)
    assert all("historical_rate" in vector and "rate_change" in vector for vector in x_rows)


def test_prediction_interpretation_separates_confidence_from_business_risk() -> None:
    excluded = _prediction_interpretation("modeling_review_label", "排除", 0.96)
    unknown_business_rule = _prediction_interpretation("predict_st", "ST11", 0.99)
    trend_alert = _prediction_interpretation("trend_regular_change", "明显上升", 0.82)

    assert excluded["confidence_level"] == "high"
    assert excluded["risk_level"] == "low"
    assert unknown_business_rule["confidence_level"] == "high"
    assert unknown_business_rule["risk_level"] == "undefined"
    assert trend_alert["risk_level"] == "high"
    assert "明显高于" in trend_alert["risk_reason"]


def test_trend_evidence_identifies_changed_indicators_and_deviation_points(tmp_path: Path) -> None:
    service = ModelingPlatformService(store=make_store(tmp_path), sample_library=Mock(), project_root=tmp_path)
    rows = []
    for month in range(1, 13):
        hit_count = 9 if month == 12 else 2
        for index in range(10):
            rows.append(
                {
                    "collection_date": f"2025-{month:02d}-15",
                    "serotype_result": "Enteritidis" if index < hit_count else "Typhimurium",
                    "resistance_gene_hits": "",
                    "virulence_gene_hits": "",
                }
            )

    evidence = service._trend_evidence(rows, signal_type="serotype", target_name="Enteritidis", lookback_months=3, min_monthly_samples=20)

    assert {item["name"] for item in evidence["indicator_changes"]} == {"目标检出率", "目标阳性数", "样本总量", "月度波动性"}
    assert evidence["summary"]["severe_count"] == 1
    assert evidence["severe_points"][0]["month"] == "2025-12"
    assert evidence["severe_points"][0]["rate"] == 0.9
    assert evidence["severe_points"][0]["sample_size_warning"] is True
    assert evidence["series"][-1]["normal_upper"] < evidence["series"][-1]["rate"]


def test_trained_model_saves_full_classification_evaluation(tmp_path: Path) -> None:
    store = make_store(tmp_path)
    rows = [
        {
            "sample_key": f"s{index}",
            "sample_name": f"S{index}",
            "species_name": "Salmonella enterica",
            "collection_date": f"2025-{(index % 12) + 1:02d}-15",
            "country": "上海" if index % 2 else "江苏",
            "sample_type": "food",
            "sample_source": "食品",
            "sequencing_method": "illumina",
            "mlst_st": "ST11" if index % 2 else "ST34",
            "serotype_result": "Enteritidis" if index % 2 else "Typhimurium",
            "resistance_gene_hits": "blaTEM-1B" if index % 2 else "tetA",
            "virulence_gene_hits": "invA",
            "custom_metadata_json": "[]",
        }
        for index in range(80)
    ]
    sample_library = Mock()
    sample_library.list_visible.return_value = rows
    service = ModelingPlatformService(store=store, sample_library=sample_library, project_root=tmp_path)
    dataset = service.create_dataset(scope="main", role="admin", username="admin", group_name="", payload={"dataset_name": "训练集"})
    feature_set = service.create_feature_set(username="admin", payload={"feature_set_name": "基础特征", "feature_schema": {"features": ["region", "serotype", "AMR genes"], "missing_policy": "warn"}})

    result = service.train(
        scope="main",
        role="admin",
        username="admin",
        group_name="",
        payload={
            "model_name": "ST 分类模型",
            "dataset_id": dataset["dataset_id"],
            "feature_set_id": feature_set["feature_set_id"],
            "target_type": "predict_st",
            "algorithm": "logistic_regression",
            "params": {"test_size": 0.25, "random_state": 42, "split_mode": "random"},
        },
    )

    model = result["model"]
    assert model["metrics"]["macro_f1"] == 1.0
    assert len(model["metrics"]["roc_curve"]) >= 2
    assert len(model["metrics"]["pr_curve"]) >= 2
    assert set(model["class_metrics"]) == {"ST11", "ST34"}
    assert model["confusion_matrix"]["labels"] == ["ST11", "ST34"]
    assert model["training_summary"]["sample_count"] == 80
    assert model["feature_schema"]["feature_count"] == 3
    assert model["reliability_warnings"]
    assert "模型评估报告" in service.export_model_report(model["model_id"])


def test_trend_prediction_returns_auditable_change_evidence(tmp_path: Path) -> None:
    store = make_store(tmp_path)
    rows = []
    for month in range(1, 31):
        rate_count = 2 if month <= 10 or month > 20 else 8
        year = 2023 + ((month - 1) // 12)
        month_of_year = ((month - 1) % 12) + 1
        for index in range(10):
            rows.append(
                {
                    "sample_key": f"trend-{month}-{index}",
                    "sample_name": f"Trend {month}-{index}",
                    "species_name": "Salmonella enterica",
                    "collection_date": f"{year}-{month_of_year:02d}-15",
                    "country": "上海",
                    "sample_type": "food",
                    "sample_source": "食品",
                    "sequencing_method": "illumina",
                    "mlst_st": "ST11",
                    "serotype_result": "Enteritidis" if index < rate_count else "Typhimurium",
                    "resistance_gene_hits": "",
                    "virulence_gene_hits": "",
                    "custom_metadata_json": "[]",
                }
            )
    sample_library = Mock()
    sample_library.list_visible.return_value = rows
    service = ModelingPlatformService(store=store, sample_library=sample_library, project_root=tmp_path)
    dataset = service.create_dataset(scope="main", role="admin", username="admin", group_name="", payload={"dataset_name": "趋势训练集"})
    feature_set = service.create_feature_set(username="admin", payload={"feature_set_name": "趋势特征", "feature_schema": {"features": ["collection_date"]}})
    trained = service.train(
        scope="main",
        role="admin",
        username="admin",
        group_name="",
        payload={
            "model_name": "血清型趋势模型",
            "dataset_id": dataset["dataset_id"],
            "feature_set_id": feature_set["feature_set_id"],
            "target_type": "trend_regular_change",
            "target_name": "Enteritidis",
            "algorithm": "logistic_regression",
            "trend_signal_type": "serotype",
            "lookback_months": 3,
            "forecast_months": 1,
            "regular_threshold": 0.2,
            "params": {"test_size": 0.2, "split_mode": "time"},
        },
    )

    prediction = service.predict(
        model_id=trained["model"]["model_id"],
        scope="main",
        role="admin",
        username="admin",
        group_name="",
        payload={"sample_ids": [row["sample_key"] for row in rows]},
    )
    result = prediction["results"][0]

    assert result["trend_evidence"]["indicator_changes"]
    assert result["trend_evidence"]["series"]
    assert "deviation_points" in result["trend_evidence"]
    assert result["trend_evidence"]["method"].startswith("按月检出率")


def test_trend_signal_scan_trains_without_target_and_returns_candidate_signals(tmp_path: Path) -> None:
    store = make_store(tmp_path)
    rows = []
    for month in range(1, 19):
        year = 2025 + ((month - 1) // 12)
        month_of_year = ((month - 1) % 12) + 1
        amr_hits = 2 if month <= 12 else 12
        serotype_hits = 12 if month <= 12 else 3
        for index in range(20):
            rows.append(
                {
                    "sample_key": f"scan-{month}-{index}",
                    "sample_name": f"Scan {month}-{index}",
                    "species_name": "Salmonella enterica",
                    "collection_date": f"{year}-{month_of_year:02d}-15",
                    "country": "上海",
                    "sample_type": "food",
                    "sample_source": "食品",
                    "sequencing_method": "illumina",
                    "mlst_st": "ST11" if index < amr_hits else "ST34",
                    "serotype_result": "Enteritidis" if index < serotype_hits else "Typhimurium",
                    "resistance_gene_hits": "blaTEM-1B" if index < amr_hits else "rareGene" if month == 1 and index == 19 else "",
                    "virulence_gene_hits": "invA" if index < 3 else "",
                    "custom_metadata_json": "[]",
                }
            )
    sample_library = Mock()
    sample_library.list_visible.return_value = rows
    service = ModelingPlatformService(store=store, sample_library=sample_library, project_root=tmp_path)
    dataset = service.create_dataset(scope="main", role="admin", username="admin", group_name="", payload={"dataset_name": "扫描训练集"})
    feature_set = service.create_feature_set(username="admin", payload={"feature_set_name": "扫描配置", "feature_schema": {"features": ["collection_date"]}})

    trained = service.train(
        scope="main",
        role="admin",
        username="admin",
        group_name="",
        payload={
            "model_name": "自动偏移扫描",
            "dataset_id": dataset["dataset_id"],
            "feature_set_id": feature_set["feature_set_id"],
            "target_type": "trend_signal_scan",
            "target_name": "",
            "scan_signal_types": ["amr_gene", "st", "serotype", "vf_gene"],
            "lookback_months": 3,
            "forecast_months": 1,
            "regular_threshold": 0.1,
            "min_monthly_samples": 10,
            "min_hit_count": 5,
            "severe_z_threshold": 3,
        },
    )
    prediction = service.predict(
        model_id=trained["model"]["model_id"],
        scope="main",
        role="admin",
        username="admin",
        group_name="",
        payload={"sample_ids": [row["sample_key"] for row in rows]},
    )
    signal_ids = {item["signal_id"] for item in prediction["results"]}

    assert trained["model"]["algorithm"] == "rule_based_trend_scan"
    assert trained["model"]["config"]["scan_signal_types"] == ["amr_gene", "st", "serotype", "vf_gene"]
    assert "amr_gene::blaTEM-1B" in signal_ids
    assert "serotype::Enteritidis" in signal_ids
    assert "amr_gene::rareGene" not in signal_ids
    assert all("series" in item and "indicator_changes" in item["trend_evidence"] for item in prediction["results"])
