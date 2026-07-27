from __future__ import annotations

import json
import math
import pickle
import platform
import re
import statistics
import time
from html import escape as html_escape
from collections import Counter
from dataclasses import dataclass
from importlib.metadata import PackageNotFoundError, version as package_version
from pathlib import Path
from typing import Any
from uuid import uuid4

from .sample_library_manager import SampleLibraryManager
from .store import PortalStore, utc_now_iso


META_FEATURES = [
    "collection_date",
    "region",
    "sample_type",
    "source_category",
    "hospital_or_lab",
    "project_id",
    "batch_id",
    "sequencing_platform",
    "patient_group",
    "age_group",
    "gender",
]

TYPING_FEATURES = ["ST", "MLST", "cgMLST cluster", "serotype", "SNP cluster", "phylogenetic cluster"]
GENE_FEATURES = ["AMR genes", "VF genes", "plasmid replicons", "mobile elements"]
TREND_SCAN_TYPES = {
    "amr_gene": "耐药基因",
    "vf_gene": "毒力基因",
    "st": "ST/MLST",
    "serotype": "血清型",
}
STAT_FEATURES = [
    "type_historical_frequency",
    "amr_historical_frequency",
    "vf_historical_frequency",
    "region_detection_rate",
    "type_growth_rate",
    "amr_growth_rate",
    "vf_growth_rate",
]

FEATURE_OPTIONS = {
    "meta": META_FEATURES,
    "typing": TYPING_FEATURES,
    "genes": GENE_FEATURES,
    "statistics": STAT_FEATURES,
}

ALGORITHM_OPTIONS = {
    "logistic_regression": {"label": "Logistic Regression", "package": "sklearn"},
    "random_forest": {"label": "Random Forest", "package": "sklearn"},
    "gradient_boosting": {"label": "Gradient Boosting", "package": "sklearn"},
    "svm": {"label": "SVM", "package": "sklearn"},
    "knn": {"label": "KNN", "package": "sklearn"},
    "xgboost": {"label": "XGBoost", "package": "xgboost"},
    "lightgbm": {"label": "LightGBM", "package": "lightgbm"},
}

DISCLAIMER = "模型输出仅用于数据趋势分析与人工复核参考，不作为最终业务处置依据。"
METRICS_DISCLAIMER = "高指标不代表可直接采用，需要结合数据质量和业务背景。"
MODEL_MANAGER_ROLES = {"admin", "group_admin"}

BUSINESS_RISK_LABELS = {
    "high": "重点复核",
    "focus": "持续关注",
    "routine": "常规",
    "low": "低关注",
    "review": "待复核",
    "undefined": "未定义",
}


def _confidence_level(probability: float | None) -> str:
    if probability is None:
        return "unknown"
    if probability >= 0.8:
        return "high"
    if probability >= 0.55:
        return "medium"
    return "low"


def _prediction_interpretation(target_type: str, predicted_label: str, probability: float | None) -> dict[str, str]:
    """Keep model confidence separate from target-specific business risk."""
    target_type = str(target_type or "").strip()
    predicted_label = str(predicted_label or "").strip()
    confidence = _confidence_level(probability)
    risk_level = "undefined"
    risk_reason = "该预测目标尚未配置关注级别规则，不能仅凭模型置信度判定。"

    if target_type in {"trend_regular_change", "trend_signal_scan"}:
        if predicted_label in {"明显上升", "异常上升"}:
            risk_level = "high"
            risk_reason = "模型判断后续检出率可能明显高于历史常规波动范围。"
        elif predicted_label in {"明显下降", "异常下降"}:
            risk_level = "review"
            risk_reason = "模型判断后续检出率可能明显下降，建议复核采样量、检测策略和数据完整性。"
        elif predicted_label in {"明显变化", "异常变化"}:
            risk_level = "focus"
            risk_reason = "模型判断后续变化可能超出历史常规波动范围，但旧模型未区分变化方向。"
        else:
            risk_level = "routine"
            risk_reason = "模型判断后续变化仍处于历史常规波动范围。"
    elif target_type in {"modeling_review_label", ""}:
        mapping = {
            "异常信号": ("high", "模型识别到需重点复核的信号，建议优先人工复核。"),
            "重点关注": ("focus", "模型识别到持续关注信号，建议持续跟踪。"),
            "常规": ("routine", "模型判断该样本属于常规监测范围。"),
            "排除": ("low", "模型判断该样本不属于当前持续关注范围。"),
        }
        risk_level, risk_reason = mapping.get(predicted_label, (risk_level, risk_reason))
    elif target_type == "risk_level":
        normalized = predicted_label.lower()
        mapping = {"高": "high", "重点复核": "high", "high": "high", "中": "focus", "持续关注": "focus", "medium": "focus", "低": "low", "低关注": "low", "low": "low"}
        risk_level = mapping.get(normalized, "review")
        risk_reason = "关注级别来自模型预测标签，仍需结合业务证据复核。"

    confidence_text = {
        "high": "模型对该结论置信度较高。",
        "medium": "模型对该结论置信度中等，建议结合业务信息复核。",
        "low": "模型对该结论置信度较低，应优先人工复核。",
        "unknown": "当前算法未提供可解释的预测概率。",
    }[confidence]
    return {
        "confidence_level": confidence,
        "risk_level": risk_level,
        "risk_label": BUSINESS_RISK_LABELS[risk_level],
        "risk_reason": risk_reason,
        "prediction_explanation": f"{risk_reason}{confidence_text}",
    }


class EncodedClassifier:
    """Makes optional classifiers that require numeric targets work with platform labels."""

    def __init__(self, estimator: Any):
        self.estimator = estimator
        self.label_encoder = None
        self.classes_ = None

    def fit(self, x: Any, y: Any) -> "EncodedClassifier":
        from sklearn.preprocessing import LabelEncoder

        self.label_encoder = LabelEncoder()
        encoded = self.label_encoder.fit_transform(y)
        self.estimator.fit(x, encoded)
        self.classes_ = self.label_encoder.classes_
        return self

    def predict(self, x: Any) -> Any:
        encoded = self.estimator.predict(x)
        return self.label_encoder.inverse_transform(encoded.astype(int))

    def predict_proba(self, x: Any) -> Any:
        return self.estimator.predict_proba(x)


def _json_loads(value: Any, fallback: Any) -> Any:
    if isinstance(value, (dict, list)):
        return value
    if not isinstance(value, str) or not value.strip():
        return fallback
    try:
        return json.loads(value)
    except json.JSONDecodeError:
        return fallback


def _json_dumps(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True)


def _metadata_items(record: dict[str, Any]) -> list[dict[str, Any]]:
    parsed = _json_loads(record.get("custom_metadata_json"), [])
    return [item for item in parsed if isinstance(item, dict)] if isinstance(parsed, list) else []


def _metadata_value(record: dict[str, Any], key: str) -> str:
    normalized = str(key or "").strip()
    for item in _metadata_items(record):
        if str(item.get("key") or "").strip() == normalized:
            value = item.get("value")
            if isinstance(value, dict):
                return " / ".join(str(value.get(part) or "").strip() for part in ("province", "city", "district", "detail") if str(value.get(part) or "").strip())
            return str(value or "").strip()
    return ""


def _location_region(record: dict[str, Any]) -> str:
    parsed = _json_loads(record.get("location_json"), {})
    if isinstance(parsed, dict):
        for key in ("province", "city", "district", "detail"):
            text = str(parsed.get(key) or "").strip()
            if text:
                return text
    return str(record.get("country") or "").strip()


def _split_hits(value: Any) -> list[str]:
    text = str(value or "").strip()
    if not text:
        return []
    parsed = _json_loads(text, None)
    values: list[str] = []
    if isinstance(parsed, list):
        for item in parsed:
            if isinstance(item, dict):
                values.append(str(item.get("gene") or item.get("name") or item.get("subject") or "").strip())
            else:
                values.append(str(item or "").strip())
    else:
        values = [part.strip() for part in re.split(r"[,;|\n]+", text) if part.strip()]
    return sorted({item for item in values if item})


def _month_index(value: Any) -> int | None:
    text = str(value or "").strip()
    match = re.match(r"^(\d{4})-(\d{2})", text)
    if not match:
        return None
    year, month = int(match.group(1)), int(match.group(2))
    return year * 12 + month - 1 if 1 <= month <= 12 else None


def _month_label(index: int) -> str:
    year, month_zero = divmod(index, 12)
    return f"{year:04d}-{month_zero + 1:02d}"


def _sample_field(record: dict[str, Any], key: str) -> str:
    mapping = {
        "collection_date": "collection_date",
        "region": "",
        "sample_type": "sample_type",
        "source_category": "sample_source",
        "hospital_or_lab": "",
        "project_id": "task_id",
        "batch_id": "task_name",
        "sequencing_platform": "sequencing_method",
        "patient_group": "",
        "age_group": "",
        "gender": "gender",
        "ST": "mlst_st",
        "MLST": "mlst_species_name",
        "serotype": "serotype_result",
    }
    if key == "region":
        return _location_region(record)
    column = mapping.get(key, "")
    return str(record.get(column) or "").strip() if column else _metadata_value(record, key)


def _summarize_samples(rows: list[dict[str, Any]]) -> dict[str, Any]:
    dates = sorted(str(row.get("collection_date") or "").strip() for row in rows if str(row.get("collection_date") or "").strip())
    missing = {
        "collection_date": sum(1 for row in rows if not str(row.get("collection_date") or "").strip()),
        "region": sum(1 for row in rows if not _location_region(row)),
        "pathogen": sum(1 for row in rows if not str(row.get("species_name") or "").strip()),
        "amr_genes": sum(1 for row in rows if not _split_hits(row.get("resistance_gene_hits"))),
        "vf_genes": sum(1 for row in rows if not _split_hits(row.get("virulence_gene_hits"))),
    }
    return {
        "sample_count": len(rows),
        "time_span": {"start": dates[0] if dates else "", "end": dates[-1] if dates else ""},
        "region_distribution": dict(Counter(_location_region(row) or "未填写" for row in rows).most_common(20)),
        "pathogen_distribution": dict(Counter(str(row.get("species_name") or "未填写").strip() for row in rows).most_common(20)),
        "missing_values": missing,
    }


def _available_packages() -> dict[str, bool]:
    status: dict[str, bool] = {}
    for name in {"sklearn", "xgboost", "lightgbm"}:
        try:
            __import__(name)
            status[name] = True
        except Exception:
            status[name] = False
    return status


def _serialize_dataset(row: dict[str, Any]) -> dict[str, Any]:
    return {**row, "filters": _json_loads(row.get("filters_json"), {}), "sample_ids": _json_loads(row.get("sample_ids_json"), []), "summary": _json_loads(row.get("summary_json"), {})}


def _serialize_feature_set(row: dict[str, Any]) -> dict[str, Any]:
    return {**row, "feature_schema": _json_loads(row.get("feature_schema_json"), {})}


def _serialize_model(row: dict[str, Any]) -> dict[str, Any]:
    return {
        **row,
        "dataset_snapshot": _json_loads(row.get("dataset_snapshot_json"), {}),
        "config": _json_loads(row.get("config_json"), {}),
        "metrics": _json_loads(row.get("metrics_json"), {}),
        "class_metrics": _json_loads(row.get("class_metrics_json"), {}),
        "confusion_matrix": _json_loads(row.get("confusion_matrix_json"), []),
        "feature_importance": _json_loads(row.get("feature_importance_json"), []),
        "training_summary": _json_loads(row.get("training_summary_json"), {}),
        "validation_summary": _json_loads(row.get("validation_summary_json"), {}),
        "reliability_warnings": _json_loads(row.get("reliability_warnings_json"), []),
        "feature_schema": _json_loads(row.get("feature_schema_json"), {}),
        "train_params": _json_loads(row.get("train_params_json"), {}),
        "split_strategy": _json_loads(row.get("split_strategy_json"), {}),
    }


def _serialize_job(row: dict[str, Any]) -> dict[str, Any]:
    return {**row, "logs": _json_loads(row.get("logs_json"), []), "config": _json_loads(row.get("config_json"), {})}


def _serialize_prediction(row: dict[str, Any]) -> dict[str, Any]:
    return {**row, "sample_ids": _json_loads(row.get("sample_ids_json"), []), "results": _json_loads(row.get("results_json"), [])}


@dataclass
class ModelingPlatformService:
    store: PortalStore
    sample_library: SampleLibraryManager
    project_root: Path

    def options(self) -> dict[str, Any]:
        packages = _available_packages()
        algorithms = []
        for key, meta in ALGORITHM_OPTIONS.items():
            package = meta["package"]
            algorithms.append({**meta, "key": key, "available": bool(packages.get(package))})
        return {
            "feature_options": FEATURE_OPTIONS,
            "algorithms": algorithms,
            "disclaimer": DISCLAIMER,
            "metrics_disclaimer": METRICS_DISCLAIMER,
        }

    def _prediction_record_with_interpretation(self, row: dict[str, Any]) -> dict[str, Any]:
        serialized = _serialize_prediction(row)
        try:
            model = self.get_model(str(serialized.get("model_id") or ""))
        except KeyError:
            model = {}
        target_type = str(model.get("target_type") or "")
        normalized_results = []
        for result in serialized.get("results") or []:
            item = dict(result)
            item_target_type = target_type if item.get("signal_id") else "trend_regular_change" if item.get("signal_type") else target_type
            probability = item.get("predicted_probability")
            numeric_probability = float(probability) if isinstance(probability, (int, float)) else None
            interpretation = _prediction_interpretation(item_target_type, str(item.get("predicted_label") or ""), numeric_probability)
            for key, value in interpretation.items():
                item.setdefault(key, value)
            # Legacy records used confidence thresholds as risk. Replace that misleading value.
            if "prediction_explanation" not in result:
                item.update(interpretation)
            normalized_results.append(item)
        serialized["results"] = normalized_results
        serialized["target_type"] = target_type
        serialized["target_name"] = str(model.get("target_name") or "")
        serialized["model_name"] = str(model.get("model_name") or serialized.get("model_id") or "")
        return serialized

    def search_samples(self, *, scope: str, role: str, username: str, group_name: str, filters: dict[str, Any]) -> dict[str, Any]:
        rows = self.sample_library.list_visible(scope=scope, role=role, username=username, group_name=group_name)

        def text_filter(row: dict[str, Any], filter_key: str, getter) -> bool:
            wanted = str(filters.get(filter_key) or "").strip().lower()
            return not wanted or wanted in str(getter(row) or "").strip().lower()

        def bool_presence(row: dict[str, Any], filter_key: str, getter) -> bool:
            raw = filters.get(filter_key)
            if raw in (None, "", "any"):
                return True
            expected = str(raw).lower() in {"1", "true", "yes", "present"}
            return bool(getter(row)) == expected

        start = str(filters.get("collection_start") or "").strip()
        end = str(filters.get("collection_end") or "").strip()
        min_q30 = str(filters.get("min_q30") or "").strip()
        try:
            min_q30_value = float(min_q30) if min_q30 else None
        except ValueError:
            min_q30_value = None

        filtered = []
        for row in rows:
            date = str(row.get("collection_date") or "").strip()
            if start and date and date < start:
                continue
            if end and date and date > end:
                continue
            if min_q30_value is not None:
                try:
                    q30 = float(str(row.get("q30_rate") or "0").replace("%", ""))
                except ValueError:
                    q30 = 0.0
                if q30 < min_q30_value:
                    continue
            checks = [
                text_filter(row, "pathogen", lambda item: item.get("species_name")),
                text_filter(row, "project_id", lambda item: item.get("task_id") or _metadata_value(item, "project_id")),
                text_filter(row, "batch_id", lambda item: item.get("task_name") or _metadata_value(item, "batch_id")),
                text_filter(row, "region", _location_region),
                text_filter(row, "sample_type", lambda item: item.get("sample_type")),
                text_filter(row, "source_category", lambda item: item.get("sample_source")),
                text_filter(row, "sequencing_platform", lambda item: item.get("sequencing_method")),
                bool_presence(row, "has_typing", lambda item: item.get("mlst_st") or item.get("serotype_result")),
                bool_presence(row, "has_amr", lambda item: _split_hits(item.get("resistance_gene_hits"))),
                bool_presence(row, "has_vf", lambda item: _split_hits(item.get("virulence_gene_hits"))),
            ]
            if all(checks):
                filtered.append(row)
        return {"items": filtered[:500], "sample_ids": [row["sample_key"] for row in filtered], "summary": _summarize_samples(filtered), "filters": filters}

    def create_dataset(self, *, scope: str, role: str, username: str, group_name: str, payload: dict[str, Any]) -> dict[str, Any]:
        filters = payload.get("filters") if isinstance(payload.get("filters"), dict) else payload
        search = self.search_samples(scope=scope, role=role, username=username, group_name=group_name, filters=filters)
        dataset = self.store.create_modeling_dataset(
            {
                "dataset_id": f"dataset::{uuid4().hex}",
                "dataset_name": str(payload.get("dataset_name") or payload.get("name") or "未命名训练数据集").strip(),
                "scope": scope,
                "filters_json": _json_dumps(filters),
                "sample_ids_json": _json_dumps(search["sample_ids"]),
                "summary_json": _json_dumps(search["summary"]),
                "created_by": username,
                "created_at": utc_now_iso(),
            }
        )
        return _serialize_dataset(dataset)

    def list_datasets(self) -> dict[str, Any]:
        return {"items": [_serialize_dataset(row) for row in self.store.list_modeling_datasets()]}

    def get_dataset(self, dataset_id: str) -> dict[str, Any]:
        return _serialize_dataset(self.store.get_modeling_dataset(dataset_id))

    def create_feature_set(self, *, username: str, payload: dict[str, Any]) -> dict[str, Any]:
        schema = payload.get("feature_schema") if isinstance(payload.get("feature_schema"), dict) else {}
        if not schema:
            schema = {
                "features": payload.get("features") if isinstance(payload.get("features"), list) else META_FEATURES[:4] + ["ST", "serotype", "AMR genes", "VF genes"],
                "missing_policy": str(payload.get("missing_policy") or "warn"),
            }
        row = self.store.create_modeling_feature_set(
            {
                "feature_set_id": f"features::{uuid4().hex}",
                "feature_set_name": str(payload.get("feature_set_name") or payload.get("name") or "未命名特征方案").strip(),
                "version": str(payload.get("version") or "1").strip() or "1",
                "feature_schema_json": _json_dumps(schema),
                "created_by": username,
                "created_at": utc_now_iso(),
            }
        )
        return _serialize_feature_set(row)

    def list_feature_sets(self) -> dict[str, Any]:
        return {"items": [_serialize_feature_set(row) for row in self.store.list_modeling_feature_sets()], **self.options()}

    def _samples_by_ids(self, sample_ids: list[str], *, scope: str, role: str, username: str, group_name: str) -> list[dict[str, Any]]:
        visible = self.sample_library.list_visible(scope=scope, role=role, username=username, group_name=group_name)
        by_id = {str(row.get("sample_key") or ""): row for row in visible}
        return [by_id[item] for item in sample_ids if item in by_id]

    def _feature_keys(self, feature_set: dict[str, Any]) -> list[str]:
        schema = _serialize_feature_set(feature_set)["feature_schema"]
        features = schema.get("features") if isinstance(schema, dict) else []
        return [str(item or "").strip() for item in features if str(item or "").strip()] or META_FEATURES[:4]

    def _feature_vector(self, record: dict[str, Any], feature_keys: list[str]) -> tuple[dict[str, Any], list[str]]:
        vector: dict[str, Any] = {}
        missing: list[str] = []
        for key in feature_keys:
            if key == "AMR genes":
                hits = _split_hits(record.get("resistance_gene_hits"))
                vector["amr_gene_count"] = len(hits)
                vector["amr_genes"] = "|".join(hits[:30])
                if not hits:
                    missing.append(key)
            elif key == "VF genes":
                hits = _split_hits(record.get("virulence_gene_hits"))
                vector["vf_gene_count"] = len(hits)
                vector["vf_genes"] = "|".join(hits[:30])
                if not hits:
                    missing.append(key)
            elif key == "plasmid replicons":
                vector["plasmid_count"] = str(record.get("plasmid_count") or "")
                if not vector["plasmid_count"]:
                    missing.append(key)
            elif key == "mobile elements":
                hits = _split_hits(record.get("resistance_mge_hits")) + _split_hits(record.get("virulence_mge_hits"))
                vector["mobile_element_count"] = len(hits)
                if not hits:
                    missing.append(key)
            elif key in STAT_FEATURES:
                vector[key] = 0.0
                missing.append(key)
            else:
                value = _sample_field(record, key)
                vector[key] = value
                if not value:
                    missing.append(key)
        return vector, missing

    def _target_label(self, record: dict[str, Any], target_type: str, target_name: str) -> str:
        target_type = str(target_type or "").strip()
        target_name = str(target_name or "").strip()
        if target_type == "binary_amr_gene":
            return "yes" if target_name and target_name in _split_hits(record.get("resistance_gene_hits")) else "no"
        if target_type == "binary_vf_gene":
            return "yes" if target_name and target_name in _split_hits(record.get("virulence_gene_hits")) else "no"
        if target_type in {"predict_st", "target_st"}:
            return str(record.get("mlst_st") or "").strip()
        if target_type in {"predict_serotype", "target_serotype"}:
            return str(record.get("serotype_result") or "").strip()
        if target_type == "risk_level":
            return _metadata_value(record, target_name or "risk_level")
        if target_name:
            return _sample_field(record, target_name) or _metadata_value(record, target_name)
        return _metadata_value(record, "modeling_review_label")

    def _trend_target_present(self, record: dict[str, Any], signal_type: str, target_name: str) -> bool:
        if signal_type == "serotype":
            return str(record.get("serotype_result") or "").strip().lower() == target_name.lower()
        if signal_type == "st":
            return str(record.get("mlst_st") or "").strip().lower() == target_name.lower()
        if signal_type == "vf_gene":
            return target_name in _split_hits(record.get("virulence_gene_hits"))
        return target_name in _split_hits(record.get("resistance_gene_hits"))

    def _trend_candidate_names(self, rows: list[dict[str, Any]], signal_type: str) -> list[str]:
        names: set[str] = set()
        for row in rows:
            if signal_type == "amr_gene":
                names.update(_split_hits(row.get("resistance_gene_hits")))
            elif signal_type == "vf_gene":
                names.update(_split_hits(row.get("virulence_gene_hits")))
            elif signal_type == "st":
                name = str(row.get("mlst_st") or "").strip()
                if name:
                    names.add(name)
            elif signal_type == "serotype":
                name = str(row.get("serotype_result") or "").strip()
                if name:
                    names.add(name)
        return sorted(names)

    def _trend_vector(
        self,
        month_rows: dict[int, list[dict[str, Any]]],
        *,
        cutoff: int,
        lookback_months: int,
        signal_type: str,
        target_name: str,
    ) -> dict[str, Any]:
        current = [
            row
            for month in range(cutoff - lookback_months, cutoff)
            for row in month_rows.get(month, [])
        ]
        previous = [
            row
            for month in range(cutoff - (lookback_months * 2), cutoff - lookback_months)
            for row in month_rows.get(month, [])
        ]
        current_hits = sum(self._trend_target_present(row, signal_type, target_name) for row in current)
        previous_hits = sum(self._trend_target_present(row, signal_type, target_name) for row in previous)
        current_rate = current_hits / len(current) if current else 0.0
        previous_rate = previous_hits / len(previous) if previous else current_rate
        return {
            "historical_rate": round(current_rate, 6),
            "previous_rate": round(previous_rate, 6),
            "rate_change": round(current_rate - previous_rate, 6),
            "historical_sample_count": len(current),
            "previous_sample_count": len(previous),
            "historical_target_hits": current_hits,
            "lookback_months": lookback_months,
            "cutoff_month_sin": round(math.sin((cutoff % 12) * math.pi / 6), 6),
            "cutoff_month_cos": round(math.cos((cutoff % 12) * math.pi / 6), 6),
            "signal_type": signal_type,
            "target_name": target_name,
        }

    def _build_trend_dataset(
        self,
        rows: list[dict[str, Any]],
        *,
        signal_type: str,
        target_name: str,
        lookback_months: int,
        forecast_months: int,
        regular_threshold: float,
    ) -> tuple[list[dict[str, Any]], list[str], list[dict[str, Any]]]:
        if not target_name:
            raise ValueError("趋势预测必须填写目标耐药基因或血清型名称。")
        month_rows: dict[int, list[dict[str, Any]]] = {}
        for row in rows:
            month = _month_index(row.get("collection_date"))
            if month is not None:
                month_rows.setdefault(month, []).append(row)
        if not month_rows:
            raise ValueError("训练数据缺少有效采样日期，无法构建时间窗口。")
        first_month, last_month = min(month_rows), max(month_rows)
        x_rows: list[dict[str, Any]] = []
        y_rows: list[str] = []
        timeline: list[dict[str, Any]] = []
        start = first_month + lookback_months
        stop = last_month - forecast_months + 1
        for cutoff in range(start, stop + 1):
            vector = self._trend_vector(
                month_rows,
                cutoff=cutoff,
                lookback_months=lookback_months,
                signal_type=signal_type,
                target_name=target_name,
            )
            future = [
                row
                for month in range(cutoff, cutoff + forecast_months)
                for row in month_rows.get(month, [])
            ]
            if vector["historical_sample_count"] <= 0 or not future:
                continue
            future_hits = sum(self._trend_target_present(row, signal_type, target_name) for row in future)
            future_rate = future_hits / len(future)
            change = future_rate - float(vector["historical_rate"])
            if abs(change) <= regular_threshold:
                label = "常规变化"
            elif change > 0:
                label = "明显上升"
            else:
                label = "明显下降"
            x_rows.append(vector)
            y_rows.append(label)
            timeline.append(
                {
                    "collection_date": _month_label(cutoff),
                    "historical_rate": vector["historical_rate"],
                    "future_rate": round(future_rate, 6),
                    "rate_change": round(change, 6),
                    "label": label,
                    "future_sample_count": len(future),
                }
            )
        return x_rows, y_rows, timeline

    def _latest_trend_vector(
        self,
        rows: list[dict[str, Any]],
        *,
        signal_type: str,
        target_name: str,
        lookback_months: int,
    ) -> tuple[dict[str, Any], str]:
        month_rows: dict[int, list[dict[str, Any]]] = {}
        for row in rows:
            month = _month_index(row.get("collection_date"))
            if month is not None:
                month_rows.setdefault(month, []).append(row)
        if not month_rows:
            raise ValueError("待预测数据缺少有效采样日期。")
        cutoff = max(month_rows) + 1
        vector = self._trend_vector(
            month_rows,
            cutoff=cutoff,
            lookback_months=lookback_months,
            signal_type=signal_type,
            target_name=target_name,
        )
        if vector["historical_sample_count"] <= 0:
            raise ValueError("历史窗口内没有可用于趋势预测的样本。")
        return vector, _month_label(cutoff)

    def _trend_evidence(
        self,
        rows: list[dict[str, Any]],
        *,
        signal_type: str,
        target_name: str,
        lookback_months: int,
        min_monthly_samples: int = 10,
        severe_z: float = 3.0,
    ) -> dict[str, Any]:
        month_rows: dict[int, list[dict[str, Any]]] = {}
        for row in rows:
            month = _month_index(row.get("collection_date"))
            if month is not None:
                month_rows.setdefault(month, []).append(row)
        series: list[dict[str, Any]] = []
        baseline_window = max(6, lookback_months * 2)
        for month in sorted(month_rows):
            samples = month_rows[month]
            hits = sum(self._trend_target_present(row, signal_type, target_name) for row in samples)
            rate = hits / len(samples) if samples else 0.0
            history = [float(item["rate"]) for item in series[-baseline_window:]]
            baseline_mean = statistics.mean(history) if history else rate
            baseline_std = statistics.pstdev(history) if len(history) >= 2 else 0.0
            tolerance = max(0.05, baseline_std * 2)
            lower = max(0.0, baseline_mean - tolerance)
            upper = min(1.0, baseline_mean + tolerance)
            deviation = rate - baseline_mean
            deviation_score = abs(deviation) / max(0.025, baseline_std) if history else 0.0
            severity = "severe" if len(history) >= 3 and deviation_score >= severe_z else "warning" if len(history) >= 3 and (rate < lower or rate > upper) else "normal"
            series.append(
                {
                    "month": _month_label(month),
                    "sample_count": len(samples),
                    "hit_count": hits,
                    "rate": round(rate, 6),
                    "baseline_mean": round(baseline_mean, 6),
                    "normal_lower": round(lower, 6),
                    "normal_upper": round(upper, 6),
                    "deviation": round(deviation, 6),
                    "deviation_score": round(deviation_score, 3),
                    "severity": severity,
                    "direction": "上升" if deviation > 0 else "下降" if deviation < 0 else "平稳",
                    "sample_size_warning": len(samples) < min_monthly_samples,
                }
            )

        current = series[-lookback_months:]
        previous = series[-(lookback_months * 2):-lookback_months]

        def total(items: list[dict[str, Any]], key: str) -> float:
            return float(sum(float(item.get(key) or 0) for item in items))

        current_samples = total(current, "sample_count")
        previous_samples = total(previous, "sample_count")
        current_hits = total(current, "hit_count")
        previous_hits = total(previous, "hit_count")
        current_rate = current_hits / current_samples if current_samples else 0.0
        previous_rate = previous_hits / previous_samples if previous_samples else current_rate
        current_volatility = statistics.pstdev([float(item["rate"]) for item in current]) if len(current) >= 2 else 0.0
        previous_volatility = statistics.pstdev([float(item["rate"]) for item in previous]) if len(previous) >= 2 else 0.0

        def indicator(name: str, current_value: float, previous_value: float, unit: str, *, scale: float = 1.0) -> dict[str, Any]:
            change = current_value - previous_value
            relative = change / abs(previous_value) if previous_value else None
            return {
                "name": name,
                "current": round(current_value * scale, 3),
                "previous": round(previous_value * scale, 3),
                "change": round(change * scale, 3),
                "relative_change": round(relative, 3) if relative is not None else None,
                "unit": unit,
                "direction": "上升" if change > 0 else "下降" if change < 0 else "平稳",
                "impact": round(abs(change) * scale, 3),
            }

        indicators = [
            indicator("目标检出率", current_rate, previous_rate, "百分点", scale=100),
            indicator("目标阳性数", current_hits, previous_hits, "例"),
            indicator("样本总量", current_samples, previous_samples, "份"),
            indicator("月度波动性", current_volatility, previous_volatility, "百分点", scale=100),
        ]
        indicators.sort(key=lambda item: float(item["impact"]), reverse=True)
        deviation_points = sorted(
            (item for item in series if item["severity"] != "normal"),
            key=lambda item: (item["severity"] == "severe", float(item["deviation_score"])),
            reverse=True,
        )
        severe_points = [item for item in deviation_points if item["severity"] == "severe"]
        latest = series[-1] if series else {}
        return {
            "method": f"按月检出率；滚动 {baseline_window} 个月基线；常规区间为基线均值 ± 2σ，至少保留 ±5 个百分点容差；显著偏移阈值为 {severe_z:g}σ。",
            "series": series,
            "indicator_changes": indicators,
            "deviation_points": deviation_points,
            "severe_points": severe_points,
            "summary": {
                "months": len(series),
                "deviation_count": len(deviation_points),
                "severe_count": len(severe_points),
                "latest_rate": latest.get("rate"),
                "latest_normal_lower": latest.get("normal_lower"),
                "latest_normal_upper": latest.get("normal_upper"),
                "largest_change_indicator": indicators[0]["name"] if indicators else "",
                "low_sample_months": sum(bool(item["sample_size_warning"]) for item in series),
                "min_monthly_samples": min_monthly_samples,
            },
        }

    def _trend_scan_types(self, payload: dict[str, Any]) -> list[str]:
        raw = payload.get("scan_signal_types")
        if raw is None:
            raw = payload.get("trend_scan_types")
        if isinstance(raw, str):
            candidates = [part.strip() for part in re.split(r"[,;|\s]+", raw) if part.strip()]
        elif isinstance(raw, list):
            candidates = [str(part or "").strip() for part in raw if str(part or "").strip()]
        else:
            candidates = []
        selected = [item for item in candidates if item in TREND_SCAN_TYPES]
        return selected or ["amr_gene", "st", "serotype"]

    def _trend_signal_result(
        self,
        *,
        model_row: dict[str, Any],
        evidence: dict[str, Any],
        signal_type: str,
        signal_name: str,
        lookback_months: int,
        forecast_months: int,
        regular_threshold: float,
        severe_z: float,
        applicability_warnings: list[str],
    ) -> dict[str, Any]:
        indicators = evidence.get("indicator_changes") if isinstance(evidence.get("indicator_changes"), list) else []
        rate_change_item = next((item for item in indicators if item.get("name") == "目标检出率"), {})
        rate_change = float(rate_change_item.get("change") or 0) / 100
        current_rate = float(rate_change_item.get("current") or 0) / 100
        previous_rate = float(rate_change_item.get("previous") or 0) / 100
        severe_points = evidence.get("severe_points") if isinstance(evidence.get("severe_points"), list) else []
        deviation_points = evidence.get("deviation_points") if isinstance(evidence.get("deviation_points"), list) else []
        if abs(rate_change) <= regular_threshold and not severe_points:
            predicted = "常规变化"
        elif rate_change > 0:
            predicted = "明显上升"
        else:
            predicted = "明显下降"
        interpretation = _prediction_interpretation("trend_signal_scan", predicted, None)
        if predicted == "明显上升":
            interpretation.update(
                {
                    "risk_level": "high" if severe_points else "focus",
                    "risk_label": BUSINESS_RISK_LABELS["high" if severe_points else "focus"],
                    "risk_reason": f"{TREND_SCAN_TYPES.get(signal_type, signal_type)} {signal_name} 检出率较前一窗口上升 {rate_change * 100:.1f} 个百分点。",
                }
            )
        elif predicted == "明显下降":
            interpretation.update(
                {
                    "risk_level": "review",
                    "risk_label": BUSINESS_RISK_LABELS["review"],
                    "risk_reason": f"{TREND_SCAN_TYPES.get(signal_type, signal_type)} {signal_name} 检出率较前一窗口下降 {abs(rate_change) * 100:.1f} 个百分点，需复核采样量和检测策略。",
                }
            )
        else:
            interpretation.update(
                {
                    "risk_level": "routine",
                    "risk_label": BUSINESS_RISK_LABELS["routine"],
                    "risk_reason": f"{TREND_SCAN_TYPES.get(signal_type, signal_type)} {signal_name} 未超出当前偏移阈值。",
                }
            )
        interpretation["prediction_explanation"] = f"{interpretation['risk_reason']} 已扫描历史常规区间、偏移点和小样本提示。"
        return {
            "sample_id": f"signal::{signal_type}::{signal_name}",
            "sample_name": f"{signal_name} 偏移信号",
            "signal_id": f"{signal_type}::{signal_name}",
            "signal_type": signal_type,
            "signal_type_label": TREND_SCAN_TYPES.get(signal_type, signal_type),
            "signal_name": signal_name,
            "predicted_label": predicted,
            "predicted_probability": None,
            **interpretation,
            "current_rate": round(current_rate, 6),
            "previous_rate": round(previous_rate, 6),
            "rate_change": round(rate_change, 6),
            "deviation_points": deviation_points,
            "severe_points": severe_points,
            "series": evidence.get("series", []),
            "trend_evidence": evidence,
            "top_features": indicators[:4],
            "missing_features": [],
            "model_id": model_row["model_id"],
            "model_version": model_row["version"],
            "prediction_time": utc_now_iso(),
            "history_window_months": lookback_months,
            "forecast_window_months": forecast_months,
            "regular_threshold": regular_threshold,
            "severe_z_threshold": severe_z,
            "applicability_warnings": applicability_warnings,
        }

    def _scan_trend_signals(
        self,
        rows: list[dict[str, Any]],
        *,
        model_row: dict[str, Any],
        config: dict[str, Any],
        applicability_warnings: list[str],
    ) -> list[dict[str, Any]]:
        lookback_months = int(config.get("lookback_months") or 3)
        forecast_months = int(config.get("forecast_months") or 1)
        regular_threshold = float(config.get("regular_threshold") or 0.1)
        min_monthly_samples = int(config.get("min_monthly_samples") or 10)
        min_hit_count = int(config.get("min_hit_count") or 5)
        severe_z = float(config.get("severe_z_threshold") or 3)
        target_filter = str(config.get("target_name") or "").strip().lower()
        results: list[dict[str, Any]] = []
        for signal_type in self._trend_scan_types(config):
            for signal_name in self._trend_candidate_names(rows, signal_type):
                if target_filter and target_filter not in signal_name.lower():
                    continue
                evidence = self._trend_evidence(
                    rows,
                    signal_type=signal_type,
                    target_name=signal_name,
                    lookback_months=lookback_months,
                    min_monthly_samples=min_monthly_samples,
                    severe_z=severe_z,
                )
                series = evidence.get("series") if isinstance(evidence.get("series"), list) else []
                months_with_samples = sum(1 for point in series if int(point.get("sample_count") or 0) >= min_monthly_samples)
                months_with_hits = sum(1 for point in series if int(point.get("hit_count") or 0) > 0)
                total_hits = sum(int(point.get("hit_count") or 0) for point in series)
                if len(series) < max(lookback_months * 2, 6) or months_with_samples < max(lookback_months + 1, 4):
                    continue
                if total_hits < min_hit_count or months_with_hits < 2:
                    continue
                result = self._trend_signal_result(
                    model_row=model_row,
                    evidence=evidence,
                    signal_type=signal_type,
                    signal_name=signal_name,
                    lookback_months=lookback_months,
                    forecast_months=forecast_months,
                    regular_threshold=regular_threshold,
                    severe_z=severe_z,
                    applicability_warnings=applicability_warnings,
                )
                if result["predicted_label"] != "常规变化" or result["severe_points"] or result["deviation_points"]:
                    results.append(result)
        return sorted(
            results,
            key=lambda item: (len(item.get("severe_points") or []), abs(float(item.get("rate_change") or 0)), len(item.get("deviation_points") or [])),
            reverse=True,
        )

    def _build_estimator(self, algorithm: str, params: dict[str, Any]) -> Any:
        try:
            from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
            from sklearn.feature_extraction import DictVectorizer
            from sklearn.neighbors import KNeighborsClassifier
            from sklearn.pipeline import make_pipeline
            from sklearn.svm import SVC
            from sklearn.linear_model import LogisticRegression
        except Exception as exc:
            raise ValueError("当前环境未安装 scikit-learn，无法训练第一阶段模型。") from exc
        seed = int(params.get("random_state") or 42)
        if algorithm == "xgboost":
            try:
                from xgboost import XGBClassifier
            except Exception as exc:
                raise ValueError("当前环境未安装 XGBoost。") from exc
            estimator = EncodedClassifier(
                XGBClassifier(
                    n_estimators=int(params.get("n_estimators") or 160),
                    max_depth=int(params.get("max_depth") or 6),
                    learning_rate=float(params.get("learning_rate") or 0.08),
                    random_state=seed,
                    eval_metric="mlogloss",
                )
            )
        elif algorithm == "lightgbm":
            try:
                from lightgbm import LGBMClassifier
            except Exception as exc:
                raise ValueError("当前环境未安装 LightGBM。") from exc
            estimator = EncodedClassifier(
                LGBMClassifier(
                    n_estimators=int(params.get("n_estimators") or 160),
                    max_depth=int(params.get("max_depth") or -1),
                    learning_rate=float(params.get("learning_rate") or 0.08),
                    random_state=seed,
                    verbosity=-1,
                )
            )
        elif algorithm == "random_forest":
            estimator = RandomForestClassifier(n_estimators=int(params.get("n_estimators") or 120), class_weight=params.get("class_weight") or None, random_state=seed)
        elif algorithm == "gradient_boosting":
            estimator = GradientBoostingClassifier(random_state=seed)
        elif algorithm == "svm":
            estimator = SVC(probability=True, class_weight=params.get("class_weight") or None, random_state=seed)
        elif algorithm == "knn":
            estimator = KNeighborsClassifier(n_neighbors=int(params.get("n_neighbors") or 5))
        else:
            estimator = LogisticRegression(max_iter=int(params.get("max_iter") or 1000), class_weight=params.get("class_weight") or None, random_state=seed)
        return make_pipeline(DictVectorizer(sparse=False), estimator)

    def _feature_importance(self, model: Any) -> list[dict[str, Any]]:
        try:
            vectorizer = model.steps[0][1]
            estimator = model.steps[-1][1]
            feature_names = list(vectorizer.get_feature_names_out())
            if isinstance(estimator, EncodedClassifier):
                estimator = estimator.estimator
            values = getattr(estimator, "feature_importances_", None)
            if values is None:
                coefficients = getattr(estimator, "coef_", None)
                if coefficients is not None:
                    import numpy as np

                    values = np.mean(np.abs(coefficients), axis=0)
            if values is None:
                return []
            ranked = sorted(
                ({"feature": str(name), "importance": round(float(value), 6)} for name, value in zip(feature_names, values)),
                key=lambda item: item["importance"],
                reverse=True,
            )
            return ranked[:30]
        except Exception:
            return []

    def _package_version(self, algorithm: str) -> str:
        package = "xgboost" if algorithm == "xgboost" else "lightgbm" if algorithm == "lightgbm" else "scikit-learn"
        try:
            return package_version(package)
        except PackageNotFoundError:
            return ""

    def _feature_schema_summary(self, feature_set: dict[str, Any], feature_keys: list[str]) -> dict[str, Any]:
        schema = _serialize_feature_set(feature_set)["feature_schema"]
        return {
            **schema,
            "feature_set_name": str(feature_set.get("feature_set_name") or ""),
            "feature_set_version": str(feature_set.get("version") or ""),
            "feature_count": len(feature_keys),
            "meta_count": sum(key in META_FEATURES for key in feature_keys),
            "typing_count": sum(key in TYPING_FEATURES for key in feature_keys),
            "amr_feature_count": sum(key == "AMR genes" for key in feature_keys),
            "vf_feature_count": sum(key == "VF genes" for key in feature_keys),
            "other_genomic_count": sum(key in GENE_FEATURES and key not in {"AMR genes", "VF genes"} for key in feature_keys),
            "categorical_encoding": "DictVectorizer one-hot encoding",
            "normalization": "none",
            "feature_selection": "none",
        }

    def _reliability_warnings(
        self,
        *,
        total_count: int,
        test_count: int,
        label_distribution: Counter,
        split_mode: str,
        target_type: str,
        missing_sample_count: int,
        train_accuracy: float,
        test_accuracy: float,
    ) -> list[dict[str, str]]:
        warnings: list[dict[str, str]] = []
        if total_count < 100 or test_count < 30:
            warnings.append({"level": "warning", "message": "该模型测试集样本量较少，结果仅供参考。"})
        rare = [str(label) for label, count in label_distribution.items() if count < 20]
        if rare:
            warnings.append({"level": "warning", "message": f"类别 {', '.join(rare)} 样本量不足，评估结果可能不稳定。"})
        counts = list(label_distribution.values())
        if counts and min(counts) / max(counts) < 0.25:
            warnings.append({"level": "warning", "message": "目标类别分布不均衡，Accuracy 可能偏高，建议关注 Recall、Macro F1 和 PR-AUC。"})
        if target_type == "trend_regular_change" and split_mode != "time":
            warnings.append({"level": "danger", "message": "该模型未使用时间切分，不建议用于未来趋势分析。"})
        if train_accuracy - test_accuracy > 0.15:
            warnings.append({"level": "warning", "message": "训练集和测试集指标差距较大，可能存在过拟合。"})
        if total_count and missing_sample_count / total_count > 0.3:
            warnings.append({"level": "warning", "message": "超过 30% 的训练样本存在特征缺失，输出稳定性可能受影响。"})
        if not warnings:
            warnings.append({"level": "success", "message": "未发现明显可靠性问题，可用于后续复核验证。"})
        warnings.append({"level": "info", "message": DISCLAIMER})
        return warnings

    def train(self, *, scope: str, role: str, username: str, group_name: str, payload: dict[str, Any]) -> dict[str, Any]:
        if role not in MODEL_MANAGER_ROLES:
            raise PermissionError("当前用户无权启动模型训练")
        now = utc_now_iso()
        job_id = f"train::{uuid4().hex}"
        logs = ["训练任务已创建。"]
        config = dict(payload)
        self.store.upsert_modeling_train_job({"job_id": job_id, "model_id": "", "status": "running", "logs_json": _json_dumps(logs), "config_json": _json_dumps(config), "created_by": username, "created_at": now, "updated_at": now})
        try:
            dataset = self.get_dataset(str(payload.get("dataset_id") or ""))
            feature_set = self.store.get_modeling_feature_set(str(payload.get("feature_set_id") or ""))
            feature_keys = self._feature_keys(feature_set)
            sample_ids = dataset.get("sample_ids") or []
            rows = self._samples_by_ids(sample_ids, scope=scope, role=role, username=username, group_name=group_name)
            target_type = str(payload.get("target_type") or "modeling_review_label").strip()
            target_name = str(payload.get("target_name") or "").strip()
            params = payload.get("params") if isinstance(payload.get("params"), dict) else {}
            if target_type == "trend_signal_scan":
                lookback_months = int(payload.get("lookback_months") or 3)
                forecast_months = int(payload.get("forecast_months") or 1)
                regular_threshold = float(payload.get("regular_threshold") or 0.1)
                min_monthly_samples = int(payload.get("min_monthly_samples") or 10)
                min_hit_count = int(payload.get("min_hit_count") or 5)
                severe_z = float(payload.get("severe_z_threshold") or 3)
                if lookback_months not in {3, 6}:
                    raise ValueError("历史窗口仅支持 3 个月或 6 个月。")
                if forecast_months not in {1, 3}:
                    raise ValueError("后续预测窗口仅支持 1 个月或 3 个月。")
                if not 0 < regular_threshold <= 1:
                    raise ValueError("偏移阈值必须在 0 到 1 之间。")
                scan_types = self._trend_scan_types(payload)
                config.update(
                    {
                        "scan_signal_types": scan_types,
                        "target_name": target_name,
                        "lookback_months": lookback_months,
                        "forecast_months": forecast_months,
                        "regular_threshold": regular_threshold,
                        "min_monthly_samples": min_monthly_samples,
                        "min_hit_count": min_hit_count,
                        "severe_z_threshold": severe_z,
                        "candidate_rule": "从 AMR/VF 基因命中、mlst_st 和 serotype_result 自动枚举候选指标。",
                        "deviation_rule": "按月检出率滚动基线，超过偏移阈值或历史常规区间时输出偏移信号。",
                    }
                )
                row_summary = _summarize_samples(rows)
                detected_candidates = {
                    signal_type: len(self._trend_candidate_names(rows, signal_type))
                    for signal_type in scan_types
                }
                now_done = utc_now_iso()
                model_id = f"model::{uuid4().hex}"
                version = f"v{now[:10].replace('-', '')}-{model_id[-6:]}"
                artifact_dir = self.project_root / "bac_analysis_portal" / "model_artifacts"
                artifact_dir.mkdir(parents=True, exist_ok=True)
                artifact_path = artifact_dir / f"{model_id.replace(':', '_')}.pkl"
                with artifact_path.open("wb") as handle:
                    pickle.dump({"model": None, "feature_keys": [], "target_type": target_type, "target_name": target_name, "config": config}, handle)
                metrics = {
                    "algorithm": "rule_based_trend_scan",
                    "candidate_counts": detected_candidates,
                    "scan_signal_types": scan_types,
                    "sample_count": len(rows),
                    "lookback_months": lookback_months,
                    "forecast_months": forecast_months,
                    "regular_threshold": regular_threshold,
                    "min_monthly_samples": min_monthly_samples,
                    "min_hit_count": min_hit_count,
                    "severe_z_threshold": severe_z,
                    "split_mode": "time",
                    "disclaimer": METRICS_DISCLAIMER,
                }
                training_summary = {
                    "dataset_name": dataset.get("dataset_name"),
                    "dataset_id": dataset.get("dataset_id"),
                    "sample_count": len(rows),
                    "train_sample_count": len(rows),
                    "validation_sample_count": 0,
                    "test_sample_count": 0,
                    "time_span": row_summary.get("time_span", {}),
                    "region_distribution": row_summary.get("region_distribution", {}),
                    "sample_type_distribution": dict(Counter(str(row.get("sample_type") or "未填写") for row in rows)),
                    "pathogen_distribution": row_summary.get("pathogen_distribution", {}),
                    "label_distribution": {},
                    "missing_values": row_summary.get("missing_values", {}),
                    "excluded_sample_count": 0,
                    "excluded_reasons": {},
                    "candidate_counts": detected_candidates,
                }
                feature_schema = {
                    "feature_set_name": str(feature_set.get("feature_set_name") or ""),
                    "feature_set_version": str(feature_set.get("version") or ""),
                    "feature_count": len(scan_types),
                    "scan_signal_types": scan_types,
                    "candidate_rule": config["candidate_rule"],
                    "missing_policy": "低样本和低命中候选自动过滤",
                    "categorical_encoding": "none",
                    "normalization": "monthly detection rate",
                    "feature_selection": "minimum sample and hit filters",
                }
                model_row = self.store.create_modeling_model(
                    {
                        "model_id": model_id,
                        "model_name": str(payload.get("model_name") or "自动偏移扫描模型").strip(),
                        "version": version,
                        "pathogen": str(payload.get("pathogen") or next(iter(row_summary.get("pathogen_distribution", {})), "")).strip(),
                        "target_type": target_type,
                        "target_name": target_name,
                        "algorithm": "rule_based_trend_scan",
                        "feature_set_id": str(feature_set.get("feature_set_id") or ""),
                        "feature_set_version": str(feature_set.get("version") or ""),
                        "dataset_id": str(dataset.get("dataset_id") or ""),
                        "dataset_snapshot_json": _json_dumps({"sample_ids": sample_ids, "summary": dataset.get("summary"), "created_at": dataset.get("created_at")}),
                        "config_json": _json_dumps(config),
                        "metrics_json": _json_dumps(metrics),
                        "class_metrics_json": _json_dumps({}),
                        "confusion_matrix_json": _json_dumps({"labels": [], "matrix": []}),
                        "feature_importance_json": _json_dumps([]),
                        "training_summary_json": _json_dumps(training_summary),
                        "validation_summary_json": _json_dumps({"training_seconds": 0, "algorithm_version": "rule-based", "feature_importance_supported": False}),
                        "reliability_warnings_json": _json_dumps([{"level": "info", "message": "该模型为规则型偏移扫描配置，用于发现候选指标的趋势偏移信号。"}, {"level": "info", "message": DISCLAIMER}]),
                        "feature_schema_json": _json_dumps(feature_schema),
                        "train_params_json": _json_dumps(config),
                        "split_strategy_json": _json_dumps({"mode": "time", "time_split": True, "cross_validation": "none", "independent_test_set": False}),
                        "artifact_path": str(artifact_path),
                        "status": "active",
                        "created_by": username,
                        "created_at": now,
                        "updated_at": now_done,
                    }
                )
                logs.append("自动偏移扫描配置已保存，可在预测阶段自动枚举候选信号。")
                job = self.store.upsert_modeling_train_job({"job_id": job_id, "model_id": model_id, "status": "completed", "logs_json": _json_dumps(logs), "config_json": _json_dumps(config), "created_by": username, "created_at": now, "updated_at": now_done})
                return {"job": _serialize_job(job), "model": _serialize_model(model_row)}
            x_rows: list[dict[str, Any]] = []
            y_rows: list[str] = []
            model_rows: list[dict[str, Any]] = []
            missing_by_sample: dict[str, list[str]] = {}
            if target_type == "trend_regular_change":
                signal_type = str(payload.get("trend_signal_type") or "amr_gene").strip()
                lookback_months = int(payload.get("lookback_months") or 3)
                forecast_months = int(payload.get("forecast_months") or 1)
                regular_threshold = float(payload.get("regular_threshold") or 0.1)
                if lookback_months not in {3, 6}:
                    raise ValueError("历史窗口仅支持 3 个月或 6 个月。")
                if forecast_months not in {1, 3}:
                    raise ValueError("后续预测窗口仅支持 1 个月或 3 个月。")
                if not 0 < regular_threshold <= 1:
                    raise ValueError("常规变化阈值必须在 0 到 1 之间。")
                x_rows, y_rows, model_rows = self._build_trend_dataset(
                    rows,
                    signal_type=signal_type,
                    target_name=target_name,
                    lookback_months=lookback_months,
                    forecast_months=forecast_months,
                    regular_threshold=regular_threshold,
                )
                feature_keys = list(x_rows[0].keys()) if x_rows else []
                config.update(
                    {
                        "trend_signal_type": signal_type,
                        "lookback_months": lookback_months,
                        "forecast_months": forecast_months,
                        "regular_threshold": regular_threshold,
                    }
                )
                params = {**params, "split_mode": "time"}
            else:
                for row in rows:
                    label = self._target_label(row, target_type, target_name)
                    if not label:
                        continue
                    vector, missing = self._feature_vector(row, feature_keys)
                    x_rows.append(vector)
                    y_rows.append(label)
                    model_rows.append(row)
                    missing_by_sample[str(row.get("sample_key") or "")] = missing
            if len(x_rows) < 2 or len(set(y_rows)) < 2:
                raise ValueError("可训练窗口/样本不足，或历史数据中只有一种变化类别。请扩大时间范围或调整常规变化阈值。")
            algorithm = str(payload.get("algorithm") or "logistic_regression").strip()
            model = self._build_estimator(algorithm, params)
            try:
                from sklearn.metrics import (
                    accuracy_score,
                    average_precision_score,
                    classification_report,
                    confusion_matrix,
                    f1_score,
                    precision_recall_curve,
                    precision_score,
                    recall_score,
                    roc_curve,
                    roc_auc_score,
                )
                from sklearn.model_selection import train_test_split
                from sklearn.preprocessing import label_binarize
            except Exception as exc:
                raise ValueError("当前环境未安装 scikit-learn 评估组件。") from exc
            test_size = float(params.get("test_size") or 0.2)
            random_state = int(params.get("random_state") or 42)
            split_mode = str(params.get("split_mode") or "random").strip()
            indices = list(range(len(x_rows)))
            if split_mode == "time":
                ordered = sorted(indices, key=lambda idx: str(model_rows[idx].get("collection_date") or ""))
                cut = max(1, min(len(ordered) - 1, int(len(ordered) * (1 - test_size))))
                train_idx, test_idx = ordered[:cut], ordered[cut:]
            else:
                train_idx, test_idx = train_test_split(indices, test_size=test_size, random_state=random_state, stratify=y_rows if min(Counter(y_rows).values()) > 1 else None)
            training_started = time.perf_counter()
            model.fit([x_rows[idx] for idx in train_idx], [y_rows[idx] for idx in train_idx])
            predictions = list(model.predict([x_rows[idx] for idx in test_idx]))
            train_predictions = list(model.predict([x_rows[idx] for idx in train_idx]))
            training_seconds = round(time.perf_counter() - training_started, 4)
            actual = [y_rows[idx] for idx in test_idx]
            train_actual = [y_rows[idx] for idx in train_idx]
            labels = sorted(set(y_rows))
            class_report = classification_report(actual, predictions, labels=labels, output_dict=True, zero_division=0)
            class_metrics = {
                label: {
                    "precision": round(float(class_report.get(label, {}).get("precision", 0)), 4),
                    "recall": round(float(class_report.get(label, {}).get("recall", 0)), 4),
                    "f1": round(float(class_report.get(label, {}).get("f1-score", 0)), 4),
                    "support": int(class_report.get(label, {}).get("support", 0)),
                }
                for label in labels
            }
            matrix = confusion_matrix(actual, predictions, labels=labels).tolist()
            train_accuracy = float(accuracy_score(train_actual, train_predictions))
            test_accuracy = float(accuracy_score(actual, predictions))
            metrics = {
                "accuracy": round(test_accuracy, 4),
                "precision": round(float(precision_score(actual, predictions, average="weighted", zero_division=0)), 4),
                "recall": round(float(recall_score(actual, predictions, average="weighted", zero_division=0)), 4),
                "f1": round(float(f1_score(actual, predictions, average="weighted", zero_division=0)), 4),
                "macro_f1": round(float(f1_score(actual, predictions, average="macro", zero_division=0)), 4),
                "weighted_f1": round(float(f1_score(actual, predictions, average="weighted", zero_division=0)), 4),
                "macro_precision": round(float(precision_score(actual, predictions, average="macro", zero_division=0)), 4),
                "macro_recall": round(float(recall_score(actual, predictions, average="macro", zero_division=0)), 4),
                "train_accuracy": round(train_accuracy, 4),
                "train_sample_count": len(train_idx),
                "test_sample_count": len(test_idx),
                "label_distribution": dict(Counter(y_rows)),
                "split_mode": split_mode,
                "missing_feature_samples": sum(1 for items in missing_by_sample.values() if items),
                "disclaimer": METRICS_DISCLAIMER,
            }
            try:
                probabilities = model.predict_proba([x_rows[idx] for idx in test_idx])
                binary_actual = label_binarize(actual, classes=labels)
                if len(labels) == 2:
                    positive = probabilities[:, 1]
                    metrics["roc_auc"] = round(float(roc_auc_score(binary_actual.ravel(), positive)), 4)
                    metrics["pr_auc"] = round(float(average_precision_score(binary_actual.ravel(), positive)), 4)
                    false_positive_rate, true_positive_rate, _ = roc_curve(binary_actual.ravel(), positive)
                    curve_precision, curve_recall, _ = precision_recall_curve(binary_actual.ravel(), positive)
                    metrics["roc_curve"] = [
                        {"x": round(float(x), 6), "y": round(float(y), 6)}
                        for x, y in zip(false_positive_rate, true_positive_rate)
                    ]
                    metrics["pr_curve"] = [
                        {"x": round(float(x), 6), "y": round(float(y), 6)}
                        for x, y in zip(curve_recall, curve_precision)
                    ]
                    metrics["positive_label"] = labels[1]
                    metrics["positive_prevalence"] = round(float(sum(binary_actual.ravel()) / len(binary_actual)), 4)
                else:
                    metrics["roc_auc"] = round(float(roc_auc_score(binary_actual, probabilities, average="macro", multi_class="ovr")), 4)
                    metrics["pr_auc"] = round(float(average_precision_score(binary_actual, probabilities, average="macro")), 4)
                    top_k = min(3, len(labels))
                    top_predictions = [[labels[index] for index in row.argsort()[-top_k:]] for row in probabilities]
                    metrics["top_k_accuracy"] = round(sum(actual_label in predicted for actual_label, predicted in zip(actual, top_predictions)) / len(actual), 4)
            except Exception:
                pass
            if target_type == "trend_regular_change":
                metrics["time_backtest"] = model_rows
                metrics["lookback_months"] = config["lookback_months"]
                metrics["forecast_months"] = config["forecast_months"]
                metrics["regular_threshold"] = config["regular_threshold"]
            row_summary = _summarize_samples(rows)
            label_distribution = Counter(y_rows)
            training_summary = {
                "dataset_name": dataset.get("dataset_name"),
                "dataset_id": dataset.get("dataset_id"),
                "sample_count": len(x_rows),
                "train_sample_count": len(train_idx),
                "validation_sample_count": 0,
                "test_sample_count": len(test_idx),
                "time_span": row_summary.get("time_span", {}),
                "region_distribution": row_summary.get("region_distribution", {}),
                "sample_type_distribution": dict(Counter(str(row.get("sample_type") or "未填写") for row in rows)),
                "pathogen_distribution": row_summary.get("pathogen_distribution", {}),
                "label_distribution": dict(label_distribution),
                "class_balance_ratio": round(min(label_distribution.values()) / max(label_distribution.values()), 4),
                "missing_values": row_summary.get("missing_values", {}),
                "excluded_sample_count": max(0, len(rows) - len(x_rows)),
                "excluded_reasons": {"missing_target": max(0, len(rows) - len(x_rows))},
            }
            dated_train = sorted(str(model_rows[idx].get("collection_date") or "") for idx in train_idx if str(model_rows[idx].get("collection_date") or ""))
            dated_test = sorted(str(model_rows[idx].get("collection_date") or "") for idx in test_idx if str(model_rows[idx].get("collection_date") or ""))
            split_strategy = {
                "mode": split_mode,
                "time_split": split_mode == "time",
                "test_size": test_size,
                "random_state": random_state,
                "train_time_range": {"start": dated_train[0] if dated_train else "", "end": dated_train[-1] if dated_train else ""},
                "test_time_range": {"start": dated_test[0] if dated_test else "", "end": dated_test[-1] if dated_test else ""},
                "cross_validation": "none",
                "independent_test_set": False,
            }
            reliability_warnings = self._reliability_warnings(
                total_count=len(x_rows),
                test_count=len(test_idx),
                label_distribution=label_distribution,
                split_mode=split_mode,
                target_type=target_type,
                missing_sample_count=metrics["missing_feature_samples"],
                train_accuracy=train_accuracy,
                test_accuracy=test_accuracy,
            )
            validation_summary = {
                "training_seconds": training_seconds,
                "algorithm_version": self._package_version(algorithm),
                "machine": platform.platform(),
                "feature_importance_supported": bool(self._feature_importance(model)),
            }
            feature_importance = self._feature_importance(model)
            feature_schema = self._feature_schema_summary(feature_set, feature_keys)
            model_id = f"model::{uuid4().hex}"
            version = f"v{now[:10].replace('-', '')}-{model_id[-6:]}"
            artifact_dir = self.project_root / "bac_analysis_portal" / "model_artifacts"
            artifact_dir.mkdir(parents=True, exist_ok=True)
            artifact_path = artifact_dir / f"{model_id.replace(':', '_')}.pkl"
            with artifact_path.open("wb") as handle:
                pickle.dump({"model": model, "feature_keys": feature_keys, "target_type": target_type, "target_name": target_name, "config": config}, handle)
            model_row = self.store.create_modeling_model(
                {
                    "model_id": model_id,
                    "model_name": str(payload.get("model_name") or "未命名模型").strip(),
                    "version": version,
                    "pathogen": str(payload.get("pathogen") or next(iter(row_summary.get("pathogen_distribution", {})), "")).strip(),
                    "target_type": target_type,
                    "target_name": target_name,
                    "algorithm": algorithm,
                    "feature_set_id": str(feature_set.get("feature_set_id") or ""),
                    "feature_set_version": str(feature_set.get("version") or ""),
                    "dataset_id": str(dataset.get("dataset_id") or ""),
                    "dataset_snapshot_json": _json_dumps({"sample_ids": sample_ids, "summary": dataset.get("summary"), "created_at": dataset.get("created_at")}),
                    "config_json": _json_dumps(config),
                    "metrics_json": _json_dumps(metrics),
                    "class_metrics_json": _json_dumps(class_metrics),
                    "confusion_matrix_json": _json_dumps({"labels": labels, "matrix": matrix}),
                    "feature_importance_json": _json_dumps(feature_importance),
                    "training_summary_json": _json_dumps(training_summary),
                    "validation_summary_json": _json_dumps(validation_summary),
                    "reliability_warnings_json": _json_dumps(reliability_warnings),
                    "feature_schema_json": _json_dumps(feature_schema),
                    "train_params_json": _json_dumps(params),
                    "split_strategy_json": _json_dumps(split_strategy),
                    "artifact_path": str(artifact_path),
                    "status": "active",
                    "created_by": username,
                    "created_at": now,
                    "updated_at": utc_now_iso(),
                }
            )
            logs.append("训练完成，模型与评估指标已保存。")
            job = self.store.upsert_modeling_train_job({"job_id": job_id, "model_id": model_id, "status": "completed", "logs_json": _json_dumps(logs), "config_json": _json_dumps(config), "created_by": username, "created_at": now, "updated_at": utc_now_iso()})
            return {"job": _serialize_job(job), "model": _serialize_model(model_row)}
        except Exception as exc:
            logs.append(str(exc))
            job = self.store.upsert_modeling_train_job({"job_id": job_id, "model_id": "", "status": "failed", "logs_json": _json_dumps(logs), "config_json": _json_dumps(config), "created_by": username, "created_at": now, "updated_at": utc_now_iso()})
            return {"job": _serialize_job(job), "error": str(exc)}

    def get_train_job(self, job_id: str) -> dict[str, Any]:
        return _serialize_job(self.store.get_modeling_train_job(job_id))

    def list_models(self) -> dict[str, Any]:
        return {"items": [_serialize_model(row) for row in self.store.list_modeling_models()], **self.options()}

    def get_model(self, model_id: str) -> dict[str, Any]:
        return _serialize_model(self.store.get_modeling_model(model_id))

    def set_model_status(self, model_id: str, status: str, *, role: str) -> dict[str, Any]:
        if role not in MODEL_MANAGER_ROLES:
            raise PermissionError("当前用户无权启用或停用模型")
        return _serialize_model(self.store.update_modeling_model_status(model_id, status))

    def delete_model(self, model_id: str, *, role: str) -> dict[str, Any]:
        if role != "admin":
            raise PermissionError("只有管理员可以删除模型")
        self.store.delete_modeling_model(model_id)
        return {"deleted": True, "model_id": model_id}

    def export_model_report(self, model_id: str) -> str:
        model = self.get_model(model_id)
        metrics = model.get("metrics") or {}
        training = model.get("training_summary") or {}
        schema = model.get("feature_schema") or {}
        warnings = model.get("reliability_warnings") or []
        class_rows = "".join(
            f"<tr><td>{html_escape(str(label))}</td><td>{values.get('precision', '-')}</td><td>{values.get('recall', '-')}</td><td>{values.get('f1', '-')}</td><td>{values.get('support', '-')}</td></tr>"
            for label, values in (model.get("class_metrics") or {}).items()
        )
        importance_rows = "".join(
            f"<tr><td>{html_escape(str(item.get('feature') or ''))}</td><td>{item.get('importance', '-')}</td></tr>"
            for item in (model.get("feature_importance") or [])[:20]
        )
        warning_rows = "".join(f"<li>{html_escape(str(item.get('message') or item))}</li>" for item in warnings)
        matrix = model.get("confusion_matrix") or {}
        matrix_labels = matrix.get("labels") if isinstance(matrix.get("labels"), list) else []
        matrix_rows = matrix.get("matrix") if isinstance(matrix.get("matrix"), list) else []
        matrix_html = "".join(
            f"<tr><th>{html_escape(str(matrix_labels[index] if index < len(matrix_labels) else index))}</th>{''.join(f'<td>{html_escape(str(value))}</td>' for value in row)}</tr>"
            for index, row in enumerate(matrix_rows)
        )
        parameter_labels = {
            "class_weight": "类别不平衡处理",
            "random_state": "随机种子",
            "split_mode": "切分方式",
            "test_size": "测试集比例",
            "n_estimators": "迭代次数",
            "max_depth": "最大树深度",
            "learning_rate": "学习率",
            "n_neighbors": "邻居数量",
            "max_iter": "最大迭代次数",
        }

        def parameter_value(key: str, value: Any) -> str:
            text = str(value or "").strip()
            if key == "split_mode":
                return "按时间切分" if text == "time" else "随机切分" if text == "random" else text or "-"
            if key == "class_weight":
                return "自动平衡类别权重" if text == "balanced" else "不处理" if not text else text
            return "未使用" if text == "none" else text or "-"

        parameter_rows = "".join(
            f"<tr><th>{html_escape(parameter_labels.get(key, key.replace('_', ' ')))}</th><td>{html_escape(parameter_value(key, value))}</td></tr>"
            for key, value in (model.get("train_params") or {}).items()
        )
        return f"""<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><title>{html_escape(str(model.get('model_name') or '模型评估报告'))}</title>
<style>body{{font:14px/1.6 Arial,sans-serif;color:#172033;max-width:1100px;margin:32px auto;padding:0 24px}}h1,h2{{color:#0f2747}}table{{width:100%;border-collapse:collapse;margin:12px 0 24px}}th,td{{border:1px solid #dbe3ee;padding:8px;text-align:left}}th{{background:#f3f6fa}}.notice{{padding:12px;background:#fff7ed;border:1px solid #fed7aa}}</style></head><body>
<h1>模型评估报告</h1><h2>模型基本信息</h2><table><tr><th>模型名称</th><td>{html_escape(str(model.get('model_name') or ''))}</td><th>模型ID</th><td>{html_escape(str(model.get('model_id') or ''))}</td></tr><tr><th>版本</th><td>{html_escape(str(model.get('version') or ''))}</td><th>算法</th><td>{html_escape(str(model.get('algorithm') or ''))}</td></tr><tr><th>预测目标</th><td>{html_escape(str(model.get('target_name') or model.get('target_type') or ''))}</td><th>适用病原</th><td>{html_escape(str(model.get('pathogen') or ''))}</td></tr></table>
<h2>训练数据</h2><table><tr><th>数据集</th><td>{html_escape(str(training.get('dataset_name') or model.get('dataset_id') or ''))}</td><th>样本总数</th><td>{training.get('sample_count', '-')}</td></tr><tr><th>训练样本</th><td>{training.get('train_sample_count', '-')}</td><th>测试样本</th><td>{training.get('test_sample_count', '-')}</td></tr></table>
<h2>特征信息</h2><table><tr><th>特征集</th><td>{html_escape(str(schema.get('feature_set_name') or model.get('feature_set_id') or ''))}</td><th>特征总数</th><td>{schema.get('feature_count', '-')}</td></tr></table>
<h2>核心指标</h2><table><tr><th>Accuracy</th><th>Precision</th><th>Recall</th><th>F1</th><th>Macro F1</th><th>Weighted F1</th><th>ROC-AUC</th><th>PR-AUC</th></tr><tr><td>{metrics.get('accuracy', '-')}</td><td>{metrics.get('precision', '-')}</td><td>{metrics.get('recall', '-')}</td><td>{metrics.get('f1', '-')}</td><td>{metrics.get('macro_f1', '-')}</td><td>{metrics.get('weighted_f1', '-')}</td><td>{metrics.get('roc_auc', '-')}</td><td>{metrics.get('pr_auc', '-')}</td></tr></table>
<h2>分类报告</h2><table><tr><th>类别</th><th>Precision</th><th>Recall</th><th>F1</th><th>Support</th></tr>{class_rows}</table>
<h2>混淆矩阵</h2><table><tr><th>真实 \\ 预测</th>{''.join(f'<th>{html_escape(str(label))}</th>' for label in matrix_labels)}</tr>{matrix_html or '<tr><td>旧模型未保存混淆矩阵</td></tr>'}</table>
<h2>训练参数</h2><table>{parameter_rows or '<tr><td>没有额外训练参数</td></tr>'}</table>
<h2>特征重要性</h2><table><tr><th>特征</th><th>重要性</th></tr>{importance_rows or '<tr><td colspan="2">当前算法不支持特征重要性解释</td></tr>'}</table>
<h2>可靠性提示</h2><ul>{warning_rows}</ul><p class="notice">本模型输出仅用于数据趋势分析与人工复核参考，不作为最终业务处置依据。</p></body></html>"""

    def predict(self, *, model_id: str, scope: str, role: str, username: str, group_name: str, payload: dict[str, Any]) -> dict[str, Any]:
        model_row = self.get_model(model_id)
        sample_ids = payload.get("sample_ids") if isinstance(payload.get("sample_ids"), list) else []
        if not sample_ids:
            filters = payload.get("filters") if isinstance(payload.get("filters"), dict) else {}
            sample_ids = self.search_samples(scope=scope, role=role, username=username, group_name=group_name, filters=filters)["sample_ids"]
        rows = self._samples_by_ids([str(item) for item in sample_ids], scope=scope, role=role, username=username, group_name=group_name)
        artifact = Path(str(model_row.get("artifact_path") or ""))
        if not artifact.exists():
            raise ValueError("模型文件不存在，无法执行预测。")
        with artifact.open("rb") as handle:
            payload_model = pickle.load(handle)
        estimator = payload_model["model"]
        feature_keys = list(payload_model.get("feature_keys") or [])
        artifact_config = payload_model.get("config") if isinstance(payload_model.get("config"), dict) else {}
        results = []
        applicability_warnings: list[str] = []
        expected_pathogen = str(model_row.get("pathogen") or "").strip()
        row_pathogens = sorted({str(row.get("species_name") or "").strip() for row in rows if str(row.get("species_name") or "").strip()})
        if expected_pathogen and any(pathogen != expected_pathogen for pathogen in row_pathogens):
            raise ValueError(f"预测前校验失败：模型适用病原为 {expected_pathogen}，待预测样本包含 {', '.join(row_pathogens)}。")
        training_summary = model_row.get("training_summary") or {}
        covered_types = set((training_summary.get("sample_type_distribution") or {}).keys())
        outside_types = sorted({str(row.get("sample_type") or "未填写") for row in rows} - covered_types) if covered_types else []
        if outside_types:
            applicability_warnings.append(f"样本类型超出训练数据覆盖范围：{', '.join(outside_types)}")
        covered_regions = set((training_summary.get("region_distribution") or {}).keys())
        outside_regions = sorted({_location_region(row) or "未填写" for row in rows} - covered_regions) if covered_regions else []
        if outside_regions:
            applicability_warnings.append(f"地区超出训练数据覆盖范围：{', '.join(outside_regions)}")
        training_time_span = training_summary.get("time_span") if isinstance(training_summary.get("time_span"), dict) else {}
        training_end = str(training_time_span.get("end") or "")
        prediction_dates = sorted(str(row.get("collection_date") or "") for row in rows if str(row.get("collection_date") or ""))
        if training_end and prediction_dates and prediction_dates[-1] > training_end:
            applicability_warnings.append(f"待预测样本时间晚于训练数据范围（训练截止 {training_end}），请谨慎解释结果。")
        if str(payload_model.get("target_type") or "") == "trend_signal_scan":
            results = self._scan_trend_signals(
                rows,
                model_row=model_row,
                config=artifact_config,
                applicability_warnings=applicability_warnings,
            )
            if not results:
                results.append(
                    {
                        "sample_id": "signal::none",
                        "sample_name": "未发现明显偏移信号",
                        "signal_id": "none",
                        "signal_type": "summary",
                        "signal_type_label": "扫描汇总",
                        "signal_name": "未发现明显偏移信号",
                        "predicted_label": "常规变化",
                        "predicted_probability": None,
                        **_prediction_interpretation("trend_signal_scan", "常规变化", None),
                        "current_rate": 0,
                        "previous_rate": 0,
                        "rate_change": 0,
                        "deviation_points": [],
                        "severe_points": [],
                        "series": [],
                        "trend_evidence": {"series": [], "indicator_changes": [], "deviation_points": [], "severe_points": [], "summary": {"deviation_count": 0, "severe_count": 0}},
                        "top_features": [],
                        "missing_features": [],
                        "model_id": model_row["model_id"],
                        "model_version": model_row["version"],
                        "prediction_time": utc_now_iso(),
                        "applicability_warnings": applicability_warnings,
                    }
                )
        elif str(payload_model.get("target_type") or "") == "trend_regular_change":
            signal_type = str(artifact_config.get("trend_signal_type") or "amr_gene")
            target_name = str(payload_model.get("target_name") or model_row.get("target_name") or "")
            lookback_months = int(artifact_config.get("lookback_months") or 3)
            forecast_months = int(artifact_config.get("forecast_months") or 1)
            vector, cutoff_month = self._latest_trend_vector(
                rows,
                signal_type=signal_type,
                target_name=target_name,
                lookback_months=lookback_months,
            )
            trend_evidence = self._trend_evidence(
                rows,
                signal_type=signal_type,
                target_name=target_name,
                lookback_months=lookback_months,
            )
            predicted = str(estimator.predict([vector])[0])
            probability = None
            if hasattr(estimator, "predict_proba"):
                try:
                    probability = round(float(max(estimator.predict_proba([vector])[0])), 4)
                except Exception:
                    probability = None
            interpretation = _prediction_interpretation("trend_regular_change", predicted, probability)
            results.append(
                {
                    "sample_id": f"trend::{signal_type}::{target_name}",
                    "sample_name": f"{target_name} 后续变化趋势",
                    "predicted_label": predicted,
                    "predicted_probability": probability,
                    **interpretation,
                    "top_features": {
                        "historical_rate": vector["historical_rate"],
                        "previous_rate": vector["previous_rate"],
                        "rate_change": vector["rate_change"],
                        "historical_sample_count": vector["historical_sample_count"],
                    },
                    "missing_features": [],
                    "model_id": model_row["model_id"],
                    "model_version": model_row["version"],
                    "prediction_time": utc_now_iso(),
                    "history_window_months": lookback_months,
                    "forecast_window_months": forecast_months,
                    "forecast_start_month": cutoff_month,
                    "signal_type": signal_type,
                    "target_name": target_name,
                    "trend_evidence": trend_evidence,
                }
            )
        else:
            critical_missing: list[str] = []
            for row in rows:
                vector, missing = self._feature_vector(row, feature_keys)
                if feature_keys and len(missing) / len(feature_keys) > 0.5:
                    critical_missing.append(str(row.get("sample_key") or ""))
                predicted = str(estimator.predict([vector])[0])
                probability = None
                if hasattr(estimator, "predict_proba"):
                    try:
                        probability = round(float(max(estimator.predict_proba([vector])[0])), 4)
                    except Exception:
                        probability = None
                interpretation = _prediction_interpretation(str(payload_model.get("target_type") or ""), predicted, probability)
                results.append(
                    {
                        "sample_id": row.get("sample_key"),
                        "sample_name": row.get("sample_name"),
                        "predicted_label": predicted,
                        "predicted_probability": probability,
                        **interpretation,
                        "top_features": feature_keys[:8],
                        "missing_features": missing,
                        "model_id": model_row["model_id"],
                        "model_version": model_row["version"],
                        "prediction_time": utc_now_iso(),
                        "applicability_warnings": applicability_warnings,
                    }
                )
            if critical_missing:
                raise ValueError(f"预测前校验失败：{len(critical_missing)} 个样本缺失特征比例超过 50%，示例：{', '.join(critical_missing[:5])}。")
        record = self.store.create_modeling_prediction(
            {
                "prediction_id": f"prediction::{uuid4().hex}",
                "model_id": model_row["model_id"],
                "model_version": model_row["version"],
                "sample_ids_json": _json_dumps([row.get("sample_key") for row in rows]),
                "results_json": _json_dumps(results),
                "created_by": username,
                "created_at": utc_now_iso(),
            }
        )
        return {**self._prediction_record_with_interpretation(record), "disclaimer": DISCLAIMER}

    def list_predictions(self) -> dict[str, Any]:
        return {"items": [self._prediction_record_with_interpretation(row) for row in self.store.list_modeling_predictions()], "disclaimer": DISCLAIMER}

    def get_prediction(self, prediction_id: str) -> dict[str, Any]:
        return {**self._prediction_record_with_interpretation(self.store.get_modeling_prediction(prediction_id)), "disclaimer": DISCLAIMER}
