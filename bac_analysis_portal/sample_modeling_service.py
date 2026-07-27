from __future__ import annotations

import json
import math
import statistics
from dataclasses import dataclass
from typing import Any
from uuid import uuid4

from .sample_library_manager import SampleLibraryManager
from .store import PortalStore, utc_now_iso


MODELING_LABEL_FIELD = "modeling_review_label"
MODELING_LABEL_OPTIONS = ["常规", "重点关注", "异常信号", "排除"]
SUPERVISED_MIN_LABELED_SAMPLES = 30
SUPERVISED_RECOMMENDED_LABELED_SAMPLES = 50

FEATURE_COLUMNS = [
    "q30_rate",
    "completeness",
    "contamination",
    "contig_count",
    "resistance_count",
    "virulence_count",
    "has_mlst",
    "has_serotype",
    "has_collection_date",
    "has_sample_source",
    "pathogen_is_virus",
    "pathogen_is_bacteria",
]

RISK_LEVEL_LABELS = {
    "routine": "常规",
    "focus": "重点关注",
    "abnormal": "异常信号",
    "review": "需人工复核",
}

ALGORITHM_OPTIONS = {
    "auto_baseline": "自动选择基线",
    "rule_baseline": "规则基线",
    "centroid_baseline": "质心监督基线",
}


def _parse_float(value: Any) -> float | None:
    text = str(value or "").strip().replace("%", "")
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _parse_metadata_items(raw_value: Any) -> list[dict[str, Any]]:
    if isinstance(raw_value, list):
        return [item for item in raw_value if isinstance(item, dict)]
    if not isinstance(raw_value, str) or not raw_value.strip():
        return []
    try:
        parsed = json.loads(raw_value)
    except json.JSONDecodeError:
        return []
    return [item for item in parsed if isinstance(item, dict)] if isinstance(parsed, list) else []


def _metadata_value(record: dict[str, Any], field_key: str) -> str:
    for item in _parse_metadata_items(record.get("custom_metadata_json")):
        if str(item.get("key") or "").strip() == field_key:
            return str(item.get("value") or "").strip()
    return ""


def _normalize_algorithm(value: Any) -> str:
    text = str(value or "auto_baseline").strip()
    return text if text in ALGORITHM_OPTIONS else "auto_baseline"


def _normalize_sample_filter(value: Any) -> str:
    text = str(value or "all").strip()
    return text if text in {"all", "selected", "modeling_test", "labeled_only", "virus", "bacteria"} else "all"


def _feature_snapshot(record: dict[str, Any]) -> dict[str, float]:
    pathogen_type = str(record.get("pathogen_type") or "").strip().lower()
    snapshot = {
        "q30_rate": _parse_float(record.get("q30_rate")) or 0.0,
        "completeness": _parse_float(record.get("completeness")) or 0.0,
        "contamination": _parse_float(record.get("contamination")) or 0.0,
        "contig_count": _parse_float(record.get("contig_count")) or 0.0,
        "resistance_count": _parse_float(record.get("resistance_count")) or 0.0,
        "virulence_count": _parse_float(record.get("virulence_count")) or 0.0,
        "has_mlst": 1.0 if str(record.get("mlst_st") or "").strip() else 0.0,
        "has_serotype": 1.0 if str(record.get("serotype_result") or "").strip() else 0.0,
        "has_collection_date": 1.0 if str(record.get("collection_date") or "").strip() else 0.0,
        "has_sample_source": 1.0 if str(record.get("sample_source") or "").strip() else 0.0,
        "pathogen_is_virus": 1.0 if pathogen_type == "virus" else 0.0,
        "pathogen_is_bacteria": 1.0 if pathogen_type == "bacteria" else 0.0,
    }
    return snapshot


def _rule_score(record: dict[str, Any], species_counts: dict[str, int]) -> tuple[float, list[str]]:
    features = _feature_snapshot(record)
    score = 18.0
    reasons: list[str] = []
    q30 = features["q30_rate"]
    if q30 and q30 < 90:
        score += 18
        reasons.append(f"Q30 为 {q30:g}%，低于 90% 的复核阈值")
    elif q30 and q30 < 95:
        score += 8
        reasons.append(f"Q30 为 {q30:g}%，处于质量关注区间")
    completeness = features["completeness"]
    if completeness and completeness < 80:
        score += 14
        reasons.append(f"完整度为 {completeness:g}%，低于 80%")
    contamination = features["contamination"]
    if contamination and contamination > 8:
        score += 16
        reasons.append(f"污染度为 {contamination:g}%，高于 8%")
    resistance = features["resistance_count"]
    virulence = features["virulence_count"]
    if resistance >= 8:
        score += 16
        reasons.append(f"耐药命中 {int(resistance)} 项，建议结合基因明细复核")
    elif resistance >= 3:
        score += 8
        reasons.append(f"耐药命中 {int(resistance)} 项")
    if virulence >= 50:
        score += 14
        reasons.append(f"毒力命中 {int(virulence)} 项，提示毒力谱偏高")
    elif virulence >= 15:
        score += 6
        reasons.append(f"毒力命中 {int(virulence)} 项")
    species_name = str(record.get("species_name") or "").strip()
    if species_name and species_counts.get(species_name, 0) == 1:
        score += 6
        reasons.append("样本库中该物种仅出现 1 次，属于稀有条目")
    if not str(record.get("sample_source") or "").strip():
        score += 4
        reasons.append("缺少样本来源，影响专题模型解释")
    if not str(record.get("collection_date") or "").strip():
        score += 4
        reasons.append("缺少采样日期，无法进行时间偏离建模")
    if not reasons:
        reasons.append("现有质控和注释字段未触发明显风险阈值")
    return min(100.0, round(score, 1)), reasons[:5]


def _risk_level(score: float) -> str:
    if score >= 72:
        return "review"
    if score >= 56:
        return "abnormal"
    if score >= 38:
        return "focus"
    return "routine"


def _standardize(rows: list[dict[str, Any]]) -> dict[str, tuple[float, float]]:
    stats: dict[str, tuple[float, float]] = {}
    for column in FEATURE_COLUMNS:
        values = [float(row["features"].get(column) or 0) for row in rows]
        mean = statistics.fmean(values) if values else 0.0
        stdev = statistics.pstdev(values) if len(values) > 1 else 1.0
        stats[column] = (mean, stdev or 1.0)
    return stats


def _distance(features: dict[str, float], centroid: dict[str, float], stats: dict[str, tuple[float, float]]) -> float:
    total = 0.0
    for column in FEATURE_COLUMNS:
        mean, stdev = stats[column]
        left = (float(features.get(column) or 0) - mean) / stdev
        right = (float(centroid.get(column) or 0) - mean) / stdev
        total += (left - right) ** 2
    return math.sqrt(total)


def _distribution(values: list[Any]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for value in values:
        key = str(value or "").strip() or "未填写"
        counts[key] = counts.get(key, 0) + 1
    return dict(sorted(counts.items(), key=lambda item: (-item[1], item[0])))


def _rate(part: int, total: int) -> float:
    if total <= 0:
        return 0.0
    return round((part / total) * 100, 1)


def _model_version_label(run: dict[str, Any] | None) -> str:
    if not run:
        return ""
    created_at = str(run.get("created_at") or "").strip()
    date_part = created_at[:10].replace("-", "") if created_at else "unknown"
    suffix = str(run.get("run_id") or "")[-6:] or "000000"
    return f"MV-{date_part}-{suffix}"


def _safe_json_loads(value: Any, fallback: Any) -> Any:
    if isinstance(value, (dict, list)):
        return value
    if not isinstance(value, str) or not value.strip():
        return fallback
    try:
        return json.loads(value)
    except json.JSONDecodeError:
        return fallback


@dataclass
class SampleModelingService:
    store: PortalStore
    sample_library: SampleLibraryManager

    def ensure_label_templates(self) -> dict[str, Any]:
        current = self.sample_library.list_metadata_templates()
        existing_keys = {str(item.get("key") or "").strip() for item in current}
        next_items = list(current)
        added = False
        if MODELING_LABEL_FIELD not in existing_keys:
            next_items.append(
                {
                    "key": MODELING_LABEL_FIELD,
                    "label": "建模复核标签",
                    "type": "select",
                    "options": MODELING_LABEL_OPTIONS,
                    "requirement": "recommended",
                    "description": "用于样本库建模闭环的人工复核标签。",
                }
            )
            added = True
        if added:
            current = self.sample_library.save_metadata_templates(next_items)
        return {"items": current, "label_field": MODELING_LABEL_FIELD, "options": MODELING_LABEL_OPTIONS, "added": added}

    def _label_fields(self) -> list[dict[str, Any]]:
        fields = [{"key": MODELING_LABEL_FIELD, "label": "建模复核标签", "type": "select"}]
        existing = {MODELING_LABEL_FIELD}
        for item in self.sample_library.list_metadata_templates():
            key = str(item.get("key") or "").strip()
            if not key or key in existing:
                continue
            field_type = str(item.get("type") or item.get("field_type") or "").strip()
            if field_type not in {"select", "text"}:
                continue
            label = str(item.get("label") or key).strip()
            fields.append({"key": key, "label": label, "type": field_type})
            existing.add(key)
        return fields

    def _filter_samples(
        self,
        samples: list[dict[str, Any]],
        *,
        sample_filter: str,
        sample_keys: list[str] | None,
        label_field: str,
    ) -> list[dict[str, Any]]:
        normalized_filter = _normalize_sample_filter(sample_filter)
        selected_keys = {str(item or "").strip() for item in sample_keys or [] if str(item or "").strip()}
        rows = list(samples)
        if normalized_filter == "selected":
            rows = [row for row in rows if str(row.get("sample_key") or "") in selected_keys]
        elif normalized_filter == "modeling_test":
            rows = [row for row in rows if str(row.get("sample_key") or "").startswith("modeling-test-")]
        elif normalized_filter == "labeled_only":
            rows = [row for row in rows if _metadata_value(row, label_field)]
        elif normalized_filter in {"virus", "bacteria"}:
            rows = [row for row in rows if str(row.get("pathogen_type") or "").strip().lower() == normalized_filter]
        return rows

    def build_dataset(
        self,
        *,
        scope: str,
        role: str,
        username: str,
        group_name: str,
        label_field: str = MODELING_LABEL_FIELD,
        sample_filter: str = "all",
        sample_keys: list[str] | None = None,
    ) -> dict[str, Any]:
        label_field = str(label_field or MODELING_LABEL_FIELD).strip() or MODELING_LABEL_FIELD
        samples = self.sample_library.list_visible(scope=scope, role=role, username=username, group_name=group_name)
        samples = self._filter_samples(samples, sample_filter=sample_filter, sample_keys=sample_keys, label_field=label_field)
        species_counts: dict[str, int] = {}
        for record in samples:
            species = str(record.get("species_name") or "").strip()
            if species:
                species_counts[species] = species_counts.get(species, 0) + 1
        rows: list[dict[str, Any]] = []
        missing = {"collection_date": 0, "sample_source": 0, "label": 0}
        sample_names: list[str] = []
        for record in samples:
            sample_names.append(str(record.get("sample_name") or "").strip())
            label = _metadata_value(record, label_field)
            if not label:
                missing["label"] += 1
            if not str(record.get("collection_date") or "").strip():
                missing["collection_date"] += 1
            if not str(record.get("sample_source") or "").strip():
                missing["sample_source"] += 1
            score, reasons = _rule_score(record, species_counts)
            rows.append(
                {
                    "sample_key": record.get("sample_key"),
                    "sample_name": record.get("sample_name"),
                    "species_name": record.get("species_name"),
                    "pathogen_type": record.get("pathogen_type"),
                    "label": label,
                    "features": _feature_snapshot(record),
                    "rule_score": score,
                    "rule_level": _risk_level(score),
                    "rule_explanations": reasons,
                }
            )
        labeled_count = len([row for row in rows if row.get("label")])
        label_distribution = _distribution([row.get("label") for row in rows if row.get("label")])
        duplicate_sample_names = sorted([
            name for name, count in _distribution([name for name in sample_names if name]).items() if count > 1
        ])
        quality = {
            "label_distribution": label_distribution,
            "pathogen_distribution": _distribution([row.get("pathogen_type") for row in rows]),
            "sample_source_distribution": _distribution([record.get("sample_source") for record in samples]),
            "species_distribution": _distribution([row.get("species_name") for row in rows]),
            "missing_rates": {
                "collection_date": _rate(missing["collection_date"], len(rows)),
                "sample_source": _rate(missing["sample_source"], len(rows)),
                "label": _rate(missing["label"], len(rows)),
            },
            "sample_source_completeness": _rate(len(rows) - missing["sample_source"], len(rows)),
            "collection_date_completeness": _rate(len(rows) - missing["collection_date"], len(rows)),
            "label_completeness": _rate(labeled_count, len(rows)),
            "duplicate_sample_names": duplicate_sample_names[:20],
            "label_class_count": len(label_distribution),
            "supervised_min_labeled_samples": SUPERVISED_MIN_LABELED_SAMPLES,
            "supervised_recommended_labeled_samples": SUPERVISED_RECOMMENDED_LABELED_SAMPLES,
        }
        can_train_supervised = labeled_count >= SUPERVISED_MIN_LABELED_SAMPLES and len(label_distribution) >= 2
        return {
            "label_field": MODELING_LABEL_FIELD,
            "active_label_field": label_field,
            "sample_filter": _normalize_sample_filter(sample_filter),
            "selected_sample_count": len([item for item in sample_keys or [] if str(item or "").strip()]),
            "feature_columns": FEATURE_COLUMNS,
            "sample_count": len(rows),
            "labeled_count": labeled_count,
            "missing": missing,
            "quality": quality,
            "can_train_supervised": can_train_supervised,
            "supervised_gate": {
                "passed": can_train_supervised,
                "min_labeled_samples": SUPERVISED_MIN_LABELED_SAMPLES,
                "recommended_labeled_samples": SUPERVISED_RECOMMENDED_LABELED_SAMPLES,
                "labeled_count": labeled_count,
                "label_class_count": len(label_distribution),
                "message": (
                    "已达到监督基线最低门槛。"
                    if can_train_supervised
                    else f"监督基线暂不启用：至少需要 {SUPERVISED_MIN_LABELED_SAMPLES} 条有效人工标签且不少于 2 个标签类别；建议达到 {SUPERVISED_RECOMMENDED_LABELED_SAMPLES} 条后再作为稳定模型使用。"
                ),
            },
            "rows": rows,
            "preview": rows[:12],
        }

    def _serialize_run(self, run: dict[str, Any] | None) -> dict[str, Any] | None:
        if not run:
            return None
        metrics = _safe_json_loads(run.get("metrics_json"), {})
        feature_columns = _safe_json_loads(run.get("feature_columns_json"), [])
        dataset_preview = _safe_json_loads(run.get("dataset_preview_json"), [])
        return {
            **run,
            "model_version": _model_version_label(run),
            "metrics": metrics if isinstance(metrics, dict) else {},
            "feature_columns": feature_columns if isinstance(feature_columns, list) else [],
            "dataset_preview": dataset_preview if isinstance(dataset_preview, list) else [],
        }

    def readiness(
        self,
        *,
        scope: str,
        role: str,
        username: str,
        group_name: str,
        label_field: str = MODELING_LABEL_FIELD,
        sample_filter: str = "all",
        sample_keys: list[str] | None = None,
        algorithm: str = "auto_baseline",
    ) -> dict[str, Any]:
        algorithm = _normalize_algorithm(algorithm)
        dataset = self.build_dataset(scope=scope, role=role, username=username, group_name=group_name, label_field=label_field, sample_filter=sample_filter, sample_keys=sample_keys)
        runs = [self._serialize_run(run) for run in self.store.list_sample_modeling_runs(limit=6)]
        scores = self.store.list_sample_modeling_scores(limit=12)
        sample_count = int(dataset["sample_count"])
        labeled_count = int(dataset["labeled_count"])
        missing = dataset["missing"]
        readiness_score = 0
        if sample_count >= 20:
            readiness_score += 30
        elif sample_count >= 8:
            readiness_score += 18
        readiness_score += min(30, labeled_count * 5)
        readiness_score += 15 if missing["sample_source"] == 0 and sample_count else 0
        readiness_score += 15 if missing["collection_date"] == 0 and sample_count else 0
        readiness_score += 10 if scores else 0
        gate = dataset["supervised_gate"]
        return {
            "config": {
                "sample_filter": dataset["sample_filter"],
                "selected_sample_count": dataset["selected_sample_count"],
                "label_field": dataset["active_label_field"],
                "algorithm": algorithm,
            },
            "label_fields": self._label_fields(),
            "algorithm_options": [{"key": key, "label": label} for key, label in ALGORITHM_OPTIONS.items()],
            "sample_filter_options": [
                {"key": "all", "label": "全部可见样本"},
                {"key": "selected", "label": "样本列表已勾选"},
                {"key": "modeling_test", "label": "模型训练测试样本"},
                {"key": "labeled_only", "label": "已有标签样本"},
                {"key": "virus", "label": "病毒样本"},
                {"key": "bacteria", "label": "细菌样本"},
            ],
            "dataset": {key: dataset[key] for key in ("sample_count", "labeled_count", "missing", "quality", "can_train_supervised", "supervised_gate", "active_label_field", "sample_filter", "selected_sample_count", "feature_columns", "preview")},
            "readiness_score": min(100, readiness_score),
            "recommended_next_step": "可以训练可解释监督基线" if dataset["can_train_supervised"] else gate["message"],
            "latest_run": runs[0] if runs else None,
            "runs": [run for run in runs if run],
            "scores": scores,
        }

    def train_baseline(
        self,
        *,
        scope: str,
        role: str,
        username: str,
        group_name: str,
        label_field: str = MODELING_LABEL_FIELD,
        sample_filter: str = "all",
        sample_keys: list[str] | None = None,
        algorithm: str = "auto_baseline",
    ) -> dict[str, Any]:
        requested_algorithm = _normalize_algorithm(algorithm)
        dataset = self.build_dataset(scope=scope, role=role, username=username, group_name=group_name, label_field=label_field, sample_filter=sample_filter, sample_keys=sample_keys)
        rows = dataset["rows"]
        if not rows:
            raise ValueError("当前建模配置没有匹配到可训练样本")
        species_counts = {str(row.get("species_name") or ""): 0 for row in rows}
        for row in rows:
            species = str(row.get("species_name") or "")
            if species:
                species_counts[species] = species_counts.get(species, 0) + 1
        algorithm = "rule_baseline"
        metrics = {
            "mode": "规则基线",
            "supervised": False,
            "warning": dataset["supervised_gate"]["message"],
            "supervised_gate": dataset["supervised_gate"],
        }
        centroids: dict[str, dict[str, float]] = {}
        stats = _standardize(rows)
        use_supervised = requested_algorithm == "centroid_baseline" or (requested_algorithm == "auto_baseline" and dataset["can_train_supervised"])
        if requested_algorithm == "centroid_baseline" and not dataset["can_train_supervised"]:
            raise ValueError(dataset["supervised_gate"]["message"])
        if use_supervised:
            algorithm = "centroid_baseline"
            labels = sorted({str(row.get("label") or "") for row in rows if row.get("label")})
            for label in labels:
                labeled_rows = [row for row in rows if row.get("label") == label]
                centroids[label] = {
                    column: statistics.fmean(float(item["features"].get(column) or 0) for item in labeled_rows)
                    for column in FEATURE_COLUMNS
                }
            metrics = {
                "mode": "质心监督基线",
                "supervised": True,
                "label_count": len(labels),
                "labels": labels,
                "supervised_gate": dataset["supervised_gate"],
                "training_config": {
                    "sample_filter": dataset["sample_filter"],
                    "selected_sample_count": dataset["selected_sample_count"],
                    "label_field": dataset["active_label_field"],
                    "algorithm": algorithm,
                },
            }
        elif requested_algorithm == "auto_baseline":
            algorithm = "rule_baseline"
        metrics["training_config"] = {
            "sample_filter": dataset["sample_filter"],
            "selected_sample_count": dataset["selected_sample_count"],
            "label_field": dataset["active_label_field"],
            "algorithm": algorithm,
        }
        run_id = f"modeling::{uuid4().hex}"
        score_rows = []
        for row in rows:
            score = float(row.get("rule_score") or 0)
            explanations = list(row.get("rule_explanations") or [])
            predicted_label = ""
            if centroids:
                distances = {
                    label: _distance(row["features"], centroid, stats)
                    for label, centroid in centroids.items()
                }
                predicted_label = min(distances, key=distances.get)
                score = max(score, 42.0 if predicted_label == "重点关注" else score)
                score = max(score, 62.0 if predicted_label == "异常信号" else score)
                explanations.insert(0, f"监督基线预测为“{predicted_label}”")
            level = _risk_level(score)
            score_rows.append(
                {
                    "sample_key": row["sample_key"],
                    "sample_name": row["sample_name"],
                    "risk_score": score,
                    "risk_level": level,
                    "risk_label": predicted_label or RISK_LEVEL_LABELS[level],
                    "explanation_json": json.dumps(explanations[:5], ensure_ascii=False),
                    "feature_snapshot_json": json.dumps(row["features"], ensure_ascii=False, sort_keys=True),
                }
            )
        run = self.store.create_sample_modeling_run(
            {
                "run_id": run_id,
                "run_type": "baseline",
                "algorithm": algorithm,
                "label_field": dataset["active_label_field"],
                "sample_count": dataset["sample_count"],
                "labeled_count": dataset["labeled_count"],
                "feature_columns_json": json.dumps(FEATURE_COLUMNS, ensure_ascii=False),
                "metrics_json": json.dumps(metrics, ensure_ascii=False),
                "dataset_preview_json": json.dumps(dataset["preview"], ensure_ascii=False),
                "status": "completed",
                "message": "已生成训练集并完成样本库风险评分回写。",
                "created_by": username,
                "created_at": utc_now_iso(),
            }
        )
        scores = self.store.replace_sample_modeling_scores(run_id, score_rows)
        serialized_run = self._serialize_run(run) or run
        return {"run": serialized_run, "scores": scores, "dataset": {key: dataset[key] for key in ("sample_count", "labeled_count", "missing", "quality", "can_train_supervised", "supervised_gate", "active_label_field", "sample_filter", "selected_sample_count")}}
