from __future__ import annotations

import csv
import json
import math
import random
import re
import zipfile
from pathlib import Path

from .admin_runtime import _resolve_workstation_module
from .parse_utils import _safe_float, _safe_int
from .table_io import _read_tsv_rows

def _find_column_index(columns: list, candidates: list[str]) -> int:
    normalized = {str(column or "").strip().lower(): index for index, column in enumerate(columns)}
    for candidate in candidates:
        key = str(candidate or "").strip().lower()
        if key in normalized:
            return normalized[key]
    return -1

def _cell_at(row: list, index: int) -> str:
    if index < 0 or index >= len(row):
        return ""
    return str(row[index] or "").strip()

def _first_row_mapping(table: dict) -> dict[str, str]:
    columns = table.get("columns") or []
    rows = table.get("rows") or []
    first = rows[0] if rows and isinstance(rows[0], list) else []
    return {str(column or "").strip(): _cell_at(first, index) for index, column in enumerate(columns)}

def _first_table_cell(row: list, columns: list, candidates: list[str]) -> str:
    index = _find_column_index(columns, candidates)
    return _cell_at(row, index)

def _first_mapping_value(mapping: dict, candidates: list[str]) -> str:
    normalized = {str(key or "").strip().lower(): str(value or "").strip() for key, value in (mapping or {}).items()}
    for candidate in candidates:
        value = normalized.get(str(candidate or "").strip().lower(), "")
        if value and value != "-":
            return value
    return ""

def _is_community_report_task(task: dict) -> bool:
    params = task.get("params") or {}
    workstation_key = _resolve_workstation_module(params.get("workstation_key")).get("key", "")
    pipeline_script = str(task.get("pipeline_script") or "").strip()
    return workstation_key == "community" or pipeline_script.endswith("CommunityAnalysis.py")

def _read_community_summary(path: Path) -> dict:
    if not path.is_file():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}

def _read_jsonl_records(path: Path) -> tuple[dict, list[dict]]:
    if not path.is_file():
        return {}, []
    header: dict = {}
    rows: list[dict] = []
    try:
        with path.open("r", encoding="utf-8", errors="ignore") as handle:
            for index, line in enumerate(handle):
                raw = str(line or "").strip()
                if not raw:
                    continue
                try:
                    record = json.loads(raw)
                except json.JSONDecodeError:
                    continue
                if index == 0 and isinstance(record, dict) and "fields" in record:
                    header = record
                    continue
                if isinstance(record, dict):
                    rows.append(record)
    except OSError:
        return {}, []
    return header, rows

def _read_community_taxonomy_preview(report_dir: Path, preview_size: int = 10) -> dict:
    taxonomy_path = report_dir / "taxonomy_export" / "taxonomy.tsv"
    raw = _read_tsv_rows(taxonomy_path)
    columns = raw.get("columns", [])
    rows = raw.get("rows", [])
    if not columns or not rows:
        return {"summary": {}, "columns": [], "rows": []}
    column_labels = {
        "Feature ID": "特征 ID",
        "Taxon": "分类注释",
        "Confidence": "置信度",
    }
    preview_columns = [column_labels.get(column, column) for column in columns[:3]]
    preview_rows = [row[:3] for row in rows[:preview_size]]
    confidences = [_safe_float(row[2] if len(row) > 2 else None) for row in rows]
    valid_confidences = [value for value in confidences if value is not None]
    return {
        "summary": {
            "feature_count": len(rows),
            "avg_confidence": round(sum(valid_confidences) / len(valid_confidences), 4) if valid_confidences else None,
        },
        "columns": preview_columns,
        "rows": preview_rows,
    }

def _community_taxonomy_rank_name(level: int) -> str:
    return {
        2: "门",
        3: "纲",
        4: "目",
        5: "科",
        6: "属",
        7: "种",
    }.get(level, f"Level {level}")

def _community_taxonomy_label(raw: str, level: int) -> str:
    parts = [part.strip() for part in str(raw or "").split(";") if part.strip()]
    target = parts[-1] if parts else str(raw or "").strip()
    target = re.sub(r"^[a-z]__+", "", target).strip("_")
    if not target:
        target = "未注释"
    return f"{_community_taxonomy_rank_name(level)}:{target}"

def _is_community_taxonomy_column(value: str) -> bool:
    text = str(value or "").strip()
    if not text or text == "index":
        return False
    if text.startswith("d__"):
        return True
    if text.startswith("Unassigned"):
        return True
    return False

def _read_metadata_group_lookup(metadata_path: Path, sample_id_column: str, group_column: str) -> dict[str, str]:
    if not metadata_path.is_file():
        return {}
    try:
        with metadata_path.open("r", encoding="utf-8", errors="ignore") as handle:
            lines = [line.rstrip("\n") for line in handle if line.strip()]
    except OSError:
        return {}
    if not lines:
        return {}
    reader = csv.DictReader(lines, delimiter="\t")
    lookup: dict[str, str] = {}
    for row in reader:
        if any(str(value or "").strip().startswith("#q2:") for value in row.values()):
            continue
        sample_id = str(row.get(sample_id_column) or "").strip()
        if not sample_id:
            continue
        lookup[sample_id] = str(row.get(group_column) or "未分组").strip() or "未分组"
    return lookup

def _find_qzv_member(names: list[str], suffix: str) -> str:
    for name in names:
        if name.endswith(suffix):
            return name
    raise FileNotFoundError(f"未在 qzv 中找到 {suffix}")

def _read_qzv_text(path: Path, suffix: str) -> str:
    with zipfile.ZipFile(path) as archive:
        member = _find_qzv_member(archive.namelist(), suffix)
        return archive.read(member).decode("utf-8", errors="ignore")

def _read_qzv_tsv_rows(path: Path, suffix: str) -> tuple[list[str], list[list[str]]]:
    raw_text = _read_qzv_text(path, suffix)
    lines = [line for line in raw_text.splitlines() if line.strip()]
    if not lines:
        return [], []
    reader = csv.reader(lines, delimiter="\t")
    rows = list(reader)
    if not rows:
        return [], []
    columns = [str(item or "").strip() for item in rows[0]]
    data_rows: list[list[str]] = []
    for row in rows[1:]:
        if row and str(row[0] or "").strip().startswith("#q2:"):
            continue
        data_rows.append([str(item or "").strip() for item in row[: len(columns)]])
    return columns, data_rows

def _read_community_demux_preview(report_dir: Path, preview_size: int = 12) -> dict:
    qzv_path = report_dir / "demux.qzv"
    columns: list[str] = []
    rows: list[list[str]] = []
    if qzv_path.is_file():
        try:
            columns, rows = _read_qzv_tsv_rows(qzv_path, "data/per-sample-fastq-counts.tsv")
        except Exception:
            columns, rows = [], []
    if not columns or not rows:
        demux_path = report_dir / "community_demux_summary.tsv"
        raw = _read_tsv_rows(demux_path)
        columns = raw.get("columns", [])
        rows = raw.get("rows", [])
    if not columns or not rows:
        return {"summary": {}, "columns": [], "rows": []}
    column_labels = {
        "sample ID": "样本 ID",
        "sample_id": "样本 ID",
        "forward sequence count": "正向 reads 数",
        "forward_count": "正向 reads 数",
        "reverse sequence count": "反向 reads 数",
        "reverse_count": "反向 reads 数",
    }
    columns = [column_labels.get(column, column) for column in columns]
    forward_values = [_safe_int(row[1] if len(row) > 1 else None) for row in rows]
    reverse_values = [_safe_int(row[2] if len(row) > 2 else None) for row in rows]
    valid_forward = [value for value in forward_values if value is not None]
    valid_reverse = [value for value in reverse_values if value is not None]
    preview_rows = rows[:preview_size]
    return {
        "summary": {
            "sample_count": len(rows),
            "forward_total": sum(valid_forward) if valid_forward else None,
            "reverse_total": sum(valid_reverse) if valid_reverse else None,
        },
        "columns": columns,
        "rows": preview_rows,
    }

def _read_community_denoise_preview(report_dir: Path, preview_size: int = 12) -> dict:
    qzv_path = report_dir / "denoising-stats-dada2.qzv"
    if not qzv_path.is_file():
        return {"summary": {}, "columns": [], "rows": []}
    try:
        columns, rows = _read_qzv_tsv_rows(qzv_path, "data/metadata.tsv")
    except Exception:
        return {"summary": {}, "columns": [], "rows": []}
    if not columns or not rows:
        return {"summary": {}, "columns": [], "rows": []}
    filtered_pct_idx = columns.index("percentage of input passed filter") if "percentage of input passed filter" in columns else -1
    non_chimeric_pct_idx = columns.index("percentage of input non-chimeric") if "percentage of input non-chimeric" in columns else -1
    merged_idx = columns.index("merged") if "merged" in columns else -1
    preview_columns = [
        "sample-id",
        "input",
        "filtered",
        "merged",
        "non-chimeric",
        "percentage of input non-chimeric",
    ]
    preview_indexes = [columns.index(name) for name in preview_columns if name in columns]
    preview_rows = [[row[index] if index < len(row) else "" for index in preview_indexes] for row in rows[:preview_size]]
    valid_filtered_pct = [_safe_float(row[filtered_pct_idx]) for row in rows] if filtered_pct_idx >= 0 else []
    valid_non_chimeric_pct = [_safe_float(row[non_chimeric_pct_idx]) for row in rows] if non_chimeric_pct_idx >= 0 else []
    valid_merged = [_safe_int(row[merged_idx]) for row in rows] if merged_idx >= 0 else []
    clean_filtered_pct = [value for value in valid_filtered_pct if value is not None]
    clean_non_chimeric_pct = [value for value in valid_non_chimeric_pct if value is not None]
    clean_merged = [value for value in valid_merged if value is not None]
    column_labels = {
        "sample-id": "样本 ID",
        "input": "输入 reads",
        "filtered": "过滤后 reads",
        "denoised": "去噪后 reads",
        "merged": "合并后 reads",
        "non-chimeric": "非嵌合 reads",
        "percentage of input passed filter": "过滤保留率(%)",
        "percentage of input merged": "合并保留率(%)",
        "percentage of input non-chimeric": "非嵌合保留率(%)",
    }
    return {
        "summary": {
            "sample_count": len(rows),
            "avg_pass_filter_pct": round(sum(clean_filtered_pct) / len(clean_filtered_pct), 2) if clean_filtered_pct else None,
            "avg_non_chimeric_pct": round(sum(clean_non_chimeric_pct) / len(clean_non_chimeric_pct), 2) if clean_non_chimeric_pct else None,
            "merged_total": sum(clean_merged) if clean_merged else None,
        },
        "columns": [column_labels.get(columns[index], columns[index]) for index in preview_indexes],
        "rows": preview_rows,
    }

def _read_community_taxa_abundance_summary(report_dir: Path, metadata_groups: dict[str, str] | None = None, levels: tuple[int, ...] = (2, 4, 5, 6), top_n: int = 15, sample_limit: int = 120) -> dict:
    qzv_path = report_dir / "taxa-barplot.qzv"
    if not qzv_path.is_file():
        return {"summary": {}, "columns": [], "rows": [], "levels": {}}
    summary_rows: list[list[str]] = []
    level_counts: dict[str, int] = {}
    level_map: dict[str, dict[str, object]] = {}
    try:
        with zipfile.ZipFile(qzv_path) as archive:
            names = archive.namelist()
            for level in levels:
                try:
                    member = _find_qzv_member(names, f"data/level-{level}.csv")
                except FileNotFoundError:
                    continue
                text = archive.read(member).decode("utf-8", errors="ignore")
                reader = csv.reader(line for line in text.splitlines() if line.strip())
                rows = list(reader)
                if len(rows) < 2:
                    continue
                header = rows[0]
                value_columns = []
                for index, column in enumerate(header[1:], start=1):
                    text_column = str(column or "").strip()
                    if not text_column:
                        continue
                    if not _is_community_taxonomy_column(text_column):
                        continue
                    value_columns.append((index, text_column))
                if not value_columns:
                    continue
                taxa_metrics: list[tuple[str, float, int]] = []
                for index, taxon in value_columns:
                    total = 0.0
                    present = 0
                    sample_totals = 0
                    for row in rows[1:]:
                        numeric_values = []
                        for value_index, _taxon in value_columns:
                            numeric_values.append(_safe_float(row[value_index] if value_index < len(row) else None) or 0.0)
                        sample_sum = sum(numeric_values)
                        if sample_sum <= 0:
                            continue
                        sample_totals += 1
                        value = _safe_float(row[index] if index < len(row) else None) or 0.0
                        relative = (value / sample_sum) * 100 if sample_sum else 0.0
                        total += relative
                        if value > 0:
                            present += 1
                    if sample_totals <= 0:
                        continue
                    taxa_metrics.append((taxon, total / sample_totals, present))
                taxa_metrics.sort(key=lambda item: item[1], reverse=True)
                top_rows = taxa_metrics[:top_n]
                level_label = _community_taxonomy_rank_name(level)
                level_counts[level_label] = len(top_rows)
                level_rows: list[list[str]] = []
                sample_series: list[dict[str, object]] = []
                for row in rows[1:1 + sample_limit]:
                    sample_name = str(row[0] or "").strip()
                    if not sample_name:
                        continue
                    sample_numeric = {}
                    sample_sum = 0.0
                    for index, taxon in value_columns:
                        value = _safe_float(row[index] if index < len(row) else None) or 0.0
                        sample_numeric[taxon] = value
                        sample_sum += value
                    if sample_sum <= 0:
                        continue
                    segments = []
                    for taxon, _relative_mean, _present in top_rows:
                        value = sample_numeric.get(taxon, 0.0)
                        ratio = (value / sample_sum) * 100 if sample_sum else 0.0
                        if ratio <= 0:
                            continue
                        segments.append({
                            "label": _community_taxonomy_label(taxon, level),
                            "ratio": round(ratio, 2),
                        })
                    other_ratio = max(0.0, 100.0 - sum(item["ratio"] for item in segments))
                    if other_ratio > 0.01:
                        segments.append({"label": "其他", "ratio": round(other_ratio, 2)})
                    sample_series.append({
                        "sample": sample_name,
                        "group": (metadata_groups or {}).get(sample_name, "未分组"),
                        "segments": segments,
                    })
                for taxon, relative_mean, present in top_rows:
                    row = [
                        level_label,
                        _community_taxonomy_label(taxon, level),
                        f"{relative_mean:.2f}%",
                        str(present),
                    ]
                    summary_rows.append(row)
                    level_rows.append(row[1:])
                level_map[level_label] = {
                    "level": level_label,
                    "columns": ["分类单元", "平均相对丰度", "检出样本数"],
                    "rows": level_rows,
                    "sample_series": sample_series,
                }
    except OSError:
        return {"summary": {}, "columns": [], "rows": [], "levels": {}}
    return {
        "summary": {
            "level_count": len(level_counts),
            "row_count": len(summary_rows),
        },
        "columns": ["分类水平", "分类单元", "平均相对丰度", "检出样本数"],
        "rows": summary_rows,
        "levels": level_map,
    }

def _quantile(sorted_values: list[float], fraction: float) -> float:
    if not sorted_values:
        return 0.0
    if len(sorted_values) == 1:
        return float(sorted_values[0])
    position = (len(sorted_values) - 1) * fraction
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return float(sorted_values[lower])
    weight = position - lower
    return float(sorted_values[lower] * (1 - weight) + sorted_values[upper] * weight)

def _boxplot_stats(values: list[float]) -> dict[str, float | int]:
    ordered = sorted(float(value) for value in values)
    if not ordered:
        return {"n": 0, "min": 0.0, "q1": 0.0, "median": 0.0, "q3": 0.0, "max": 0.0, "mean": 0.0}
    return {
        "n": len(ordered),
        "min": round(ordered[0], 4),
        "q1": round(_quantile(ordered, 0.25), 4),
        "median": round(_quantile(ordered, 0.5), 4),
        "q3": round(_quantile(ordered, 0.75), 4),
        "max": round(ordered[-1], 4),
        "mean": round(sum(ordered) / len(ordered), 4),
    }

def _average_ranks(values: list[float]) -> list[float]:
    ordered = sorted(enumerate(values), key=lambda item: item[1])
    ranks = [0.0] * len(values)
    index = 0
    while index < len(ordered):
        end = index + 1
        while end < len(ordered) and ordered[end][1] == ordered[index][1]:
            end += 1
        average_rank = (index + 1 + end) / 2
        for cursor in range(index, end):
            original_index = ordered[cursor][0]
            ranks[original_index] = average_rank
        index = end
    return ranks

def _kruskal_wallis_h(groups: list[list[float]]) -> float:
    valid_groups = [group for group in groups if group]
    if len(valid_groups) < 2:
        return 0.0
    all_values = [float(value) for group in valid_groups for value in group]
    total_n = len(all_values)
    if total_n <= 1:
        return 0.0
    ranks = _average_ranks(all_values)
    cursor = 0
    rank_sums = []
    group_sizes = []
    for group in valid_groups:
        group_size = len(group)
        group_ranks = ranks[cursor:cursor + group_size]
        rank_sums.append(sum(group_ranks))
        group_sizes.append(group_size)
        cursor += group_size
    statistic = (12 / (total_n * (total_n + 1))) * sum((rank_sum ** 2) / size for rank_sum, size in zip(rank_sums, group_sizes)) - 3 * (total_n + 1)
    value_counts: dict[float, int] = {}
    for value in all_values:
        value_counts[value] = value_counts.get(value, 0) + 1
    tie_correction = 1 - sum(count ** 3 - count for count in value_counts.values()) / (total_n ** 3 - total_n) if total_n > 1 else 1
    if tie_correction > 0:
        statistic /= tie_correction
    return float(max(statistic, 0.0))

def _permutation_kruskal_pvalue(groups: list[list[float]], permutations: int = 1200, seed: int = 42) -> float | None:
    valid_groups = [group for group in groups if group]
    if len(valid_groups) < 2:
        return None
    flattened = [float(value) for group in valid_groups for value in group]
    if len(flattened) <= 1:
        return None
    group_sizes = [len(group) for group in valid_groups]
    observed = _kruskal_wallis_h(valid_groups)
    rng = random.Random(seed)
    exceed_count = 0
    shuffled = list(flattened)
    for _ in range(permutations):
        rng.shuffle(shuffled)
        cursor = 0
        permuted_groups = []
        for size in group_sizes:
            permuted_groups.append(shuffled[cursor:cursor + size])
            cursor += size
        if _kruskal_wallis_h(permuted_groups) >= observed - 1e-12:
            exceed_count += 1
    return round((exceed_count + 1) / (permutations + 1), 4)

def _normal_survival_probability(z_score: float) -> float:
    return 0.5 * math.erfc(z_score / math.sqrt(2.0))

def _mann_whitney_pairwise_pvalue(group_a: list[float], group_b: list[float]) -> float | None:
    values_a = [float(value) for value in group_a]
    values_b = [float(value) for value in group_b]
    if not values_a or not values_b:
        return None
    merged = values_a + values_b
    ranks = _average_ranks(merged)
    n1 = len(values_a)
    n2 = len(values_b)
    rank_sum_a = sum(ranks[:n1])
    u1 = rank_sum_a - (n1 * (n1 + 1)) / 2
    u2 = n1 * n2 - u1
    u_stat = min(u1, u2)
    mean_u = (n1 * n2) / 2
    value_counts: dict[float, int] = {}
    for value in merged:
        value_counts[value] = value_counts.get(value, 0) + 1
    total_n = n1 + n2
    tie_term = sum(count ** 3 - count for count in value_counts.values())
    variance = (n1 * n2 / 12) * ((total_n + 1) - tie_term / max(total_n * (total_n - 1), 1))
    if variance <= 0:
        return None
    z_score = (abs(u_stat - mean_u) - 0.5) / math.sqrt(variance)
    p_value = min(1.0, max(0.0, 2 * _normal_survival_probability(abs(z_score))))
    return round(p_value, 4)

def _community_alpha_metric_label(metric_key: str) -> str:
    labels = {
        "shannon": "Shannon 指数",
        "observed_features": "观察到的特征数",
        "pielou_evenness": "Pielou 均匀度",
        "simpson": "Simpson 指数",
        "faith_pd": "Faith 系统发育多样性",
        "chao1": "Chao1 丰富度",
        "goods_coverage": "Good's coverage",
    }
    return labels.get(metric_key, metric_key.replace("_", " ").strip() or metric_key)

def _read_alpha_metric_export_rows(export_dir: Path) -> dict[str, float]:
    target = export_dir / "alpha-diversity.tsv"
    if not target.is_file():
        return {}
    rows = _read_tsv_rows(target)
    sample_idx = rows["columns"].index("sample-id") if "sample-id" in rows["columns"] else 0
    value_idx = rows["columns"].index("alpha_diversity") if "alpha_diversity" in rows["columns"] else 1
    result: dict[str, float] = {}
    for row in rows["rows"]:
        if sample_idx >= len(row) or value_idx >= len(row):
            continue
        sample_id = str(row[sample_idx] or "").strip()
        value = _safe_float(row[value_idx])
        if sample_id and value is not None:
            result[sample_id] = round(float(value), 4)
    return result

def _read_community_alpha_rarefaction_summary(report_dir: Path, metadata_groups: dict[str, str] | None = None, preferred_depth: int | None = None, preview_size: int = 18) -> dict:
    qzv_path = report_dir / "alpha-rarefaction.qzv"
    if not qzv_path.is_file():
        return {"summary": {}, "sample_columns": [], "sample_rows": [], "group_columns": [], "group_rows": [], "boxplots": {}, "pairwise": {}, "rarefaction": {}}

    def _select_depth_columns(header: list[str], target_depth: int | None) -> tuple[int | None, list[int]]:
        depth_map: dict[int, list[int]] = {}
        for index, column in enumerate(header[1:], start=1):
            text = str(column or "").strip()
            if not text.startswith("depth-"):
                continue
            depth_text = text.split("_iter-", 1)[0].removeprefix("depth-")
            depth_value = int(_safe_float(depth_text) or 0)
            if depth_value <= 0:
                continue
            depth_map.setdefault(depth_value, []).append(index)
        if not depth_map:
            return None, []
        available_depths = sorted(depth_map)
        if target_depth and any(depth == target_depth for depth in available_depths):
            chosen_depth = target_depth
        elif target_depth:
            chosen_depth = min(available_depths, key=lambda depth: (abs(depth - target_depth), -depth))
        else:
            chosen_depth = max(available_depths)
        return chosen_depth, depth_map.get(chosen_depth, [])

    def _build_rarefaction_curves(rows: list[list[str]], header: list[str], limit_groups: int = 8) -> dict:
        depth_map: dict[int, list[int]] = {}
        for index, column in enumerate(header[1:], start=1):
            text = str(column or "").strip()
            if not text.startswith("depth-"):
                continue
            depth_text = text.split("_iter-", 1)[0].removeprefix("depth-")
            depth_value = int(_safe_float(depth_text) or 0)
            if depth_value <= 0:
                continue
            depth_map.setdefault(depth_value, []).append(index)
        if not depth_map:
            return {}
        grouped_points: dict[str, dict[int, list[float]]] = {}
        grouped_counts: dict[str, int] = {}
        for row in rows:
            sample_id = str(row[0] or "").strip()
            if not sample_id:
                continue
            group_name = (metadata_groups or {}).get(sample_id, "未分组")
            grouped_counts[group_name] = grouped_counts.get(group_name, 0) + 1
            group_bucket = grouped_points.setdefault(group_name, {})
            for depth in sorted(depth_map):
                values = [_safe_float(row[index] if index < len(row) else None) for index in depth_map[depth]]
                valid_values = [value for value in values if value is not None]
                if not valid_values:
                    continue
                group_bucket.setdefault(depth, []).append(sum(valid_values) / len(valid_values))
        curves = []
        max_y = 0.0
        ranked_groups = sorted(grouped_points, key=lambda group: (-grouped_counts.get(group, 0), group))[:limit_groups]
        for group_name in ranked_groups:
            points = []
            for depth in sorted(grouped_points[group_name]):
                values = grouped_points[group_name][depth]
                mean_value = round(sum(values) / len(values), 4)
                max_y = max(max_y, mean_value)
                points.append({"x": depth, "y": mean_value})
            if points:
                curves.append({
                    "group": group_name,
                    "n": grouped_counts.get(group_name, 0),
                    "points": points,
                })
        return {
            "label": "Observed Features 稀释曲线",
            "x_label": "测序深度",
            "y_label": "观察到的特征数",
            "curves": curves,
            "max_y": round(max_y, 4) if max_y else None,
        }

    metrics = ["shannon", "observed_features", "pielou_evenness", "simpson", "faith_pd", "chao1", "goods_coverage"]
    sample_map: dict[str, dict[str, object]] = {}
    chosen_depth: int | None = None
    rarefaction_chart: dict[str, object] = {}
    try:
        with zipfile.ZipFile(qzv_path) as archive:
            names = archive.namelist()
            for metric_key in metrics:
                try:
                    member = _find_qzv_member(names, f"data/{metric_key}.csv")
                except FileNotFoundError:
                    continue
                text = archive.read(member).decode("utf-8", errors="ignore")
                reader = csv.reader(line for line in text.splitlines() if line.strip())
                rows = list(reader)
                if len(rows) < 2:
                    continue
                if metric_key == "observed_features":
                    rarefaction_chart = _build_rarefaction_curves(rows[1:], rows[0]) or {}
                metric_depth, depth_indexes = _select_depth_columns(rows[0], preferred_depth)
                if chosen_depth is None:
                    chosen_depth = metric_depth
                if not depth_indexes:
                    continue
                for row in rows[1:]:
                    sample_id = str(row[0] or "").strip()
                    if not sample_id:
                        continue
                    values = [_safe_float(row[index] if index < len(row) else None) for index in depth_indexes]
                    valid_values = [value for value in values if value is not None]
                    if not valid_values:
                        continue
                    entry = sample_map.setdefault(sample_id, {
                        "sample": sample_id,
                        "group": (metadata_groups or {}).get(sample_id, "未分组"),
                    })
                    entry[metric_key] = round(sum(valid_values) / len(valid_values), 4)
    except (OSError, FileNotFoundError, zipfile.BadZipFile):
        return {"summary": {}, "sample_columns": [], "sample_rows": [], "group_columns": [], "group_rows": [], "boxplots": {}, "pairwise": {}, "rarefaction": {}}

    for metric in metrics:
        exported_values = _read_alpha_metric_export_rows(report_dir / f"alpha-{metric}_export")
        for sample_id, value in exported_values.items():
            entry = sample_map.setdefault(sample_id, {
                "sample": sample_id,
                "group": (metadata_groups or {}).get(sample_id, "未分组"),
            })
            entry[metric] = value

    available_metrics = [metric for metric in metrics if any(entry.get(metric) is not None for entry in sample_map.values())]
    sample_rows_source = [entry for entry in sample_map.values() if any(entry.get(metric) is not None for metric in available_metrics)]
    if not sample_rows_source:
        return {"summary": {}, "sample_columns": [], "sample_rows": [], "group_columns": [], "group_rows": [], "boxplots": {}, "pairwise": {}, "rarefaction": rarefaction_chart}

    sample_rows_source.sort(key=lambda item: (str(item.get("group") or "未分组"), str(item.get("sample") or "")))
    group_map: dict[str, list[dict[str, object]]] = {}
    for entry in sample_rows_source:
        group_name = str(entry.get("group") or "未分组")
        group_map.setdefault(group_name, []).append(entry)

    group_rows: list[list[str]] = []
    metric_group_boxes: dict[str, list[dict[str, object]]] = {metric: [] for metric in available_metrics}
    metric_group_values: dict[str, dict[str, list[float]]] = {metric: {} for metric in available_metrics}
    for group_name, items in sorted(group_map.items(), key=lambda item: item[0]):
        shannon_values = [float(item["shannon"]) for item in items if item.get("shannon") is not None]
        observed_values = [float(item["observed_features"]) for item in items if item.get("observed_features") is not None]
        shannon_mean = round(sum(shannon_values) / len(shannon_values), 4) if shannon_values else None
        observed_mean = round(sum(observed_values) / len(observed_values), 2) if observed_values else None
        group_rows.append([
            group_name,
            str(len(items)),
            f"{shannon_mean:.3f}" if shannon_mean is not None else "--",
            f"{observed_mean:.2f}" if observed_mean is not None else "--",
        ])
        for metric in available_metrics:
            metric_values = [float(item[metric]) for item in items if item.get(metric) is not None]
            if metric_values:
                metric_group_values[metric][group_name] = metric_values
                metric_group_boxes[metric].append({"label": group_name, **_boxplot_stats(metric_values)})

    shannon_all = [float(item["shannon"]) for item in sample_rows_source if item.get("shannon") is not None]
    observed_all = [float(item["observed_features"]) for item in sample_rows_source if item.get("observed_features") is not None]
    metric_pvalues = {
        metric: _permutation_kruskal_pvalue([
            values for _group, values in sorted(metric_group_values[metric].items(), key=lambda item: item[0])
        ])
        for metric in available_metrics
    }
    def _pairwise_rows(group_values: dict[str, list[float]]) -> list[list[str]]:
        rows: list[list[str]] = []
        ordered_groups = sorted(group_values)
        for left_index, left_name in enumerate(ordered_groups):
            for right_name in ordered_groups[left_index + 1:]:
                p_value = _mann_whitney_pairwise_pvalue(group_values[left_name], group_values[right_name])
                rows.append([
                    left_name,
                    right_name,
                    str(len(group_values[left_name])),
                    str(len(group_values[right_name])),
                    f"{p_value:.4f}" if p_value is not None else "--",
                    "显著" if p_value is not None and p_value < 0.05 else "不显著",
                ])
        rows.sort(key=lambda item: (_safe_float(item[4]) if item[4] != "--" else 99, item[0], item[1]))
        return rows

    pairwise_map = {
        metric: {
            "columns": ["分组 A", "分组 B", "A 组样本数", "B 组样本数", "P 值", "显著性"],
            "rows": _pairwise_rows(metric_group_values[metric]),
            "test": "Mann-Whitney U",
        }
        for metric in available_metrics
    }
    sample_rows = [
        [
            str(entry.get("sample") or "--"),
            str(entry.get("group") or "未分组"),
            f"{float(entry['shannon']):.3f}" if entry.get("shannon") is not None else "--",
            f"{float(entry['observed_features']):.2f}" if entry.get("observed_features") is not None else "--",
        ]
        for entry in sample_rows_source[:preview_size]
    ]
    return {
        "summary": {
            "selected_depth": chosen_depth,
            "sample_count": len(sample_rows_source),
            "group_count": len(group_rows),
            "shannon_mean": round(sum(shannon_all) / len(shannon_all), 4) if shannon_all else None,
            "observed_features_mean": round(sum(observed_all) / len(observed_all), 2) if observed_all else None,
            "shannon_pvalue": metric_pvalues.get("shannon"),
            "observed_features_pvalue": metric_pvalues.get("observed_features"),
            "available_metrics": available_metrics,
        },
        "sample_columns": ["样本 ID", "分组", "Shannon 指数", "观察到的特征数"],
        "sample_rows": sample_rows,
        "group_columns": ["分组", "样本数", "平均 Shannon 指数", "平均观察到的特征数"],
        "group_rows": group_rows,
        "boxplots": {
            metric: {
                "label": f"{_community_alpha_metric_label(metric)}箱线图",
                "tab_label": _community_alpha_metric_label(metric),
                "x_label": "分组",
                "y_label": _community_alpha_metric_label(metric),
                "groups": metric_group_boxes[metric],
                "p_value": metric_pvalues.get(metric),
                "significant": bool(metric_pvalues.get(metric) is not None and metric_pvalues.get(metric) < 0.05),
                "test": "置换 Kruskal-Wallis",
            }
            for metric in available_metrics
        },
        "pairwise": pairwise_map,
        "rarefaction": {
            **(rarefaction_chart or {}),
            "suggested_depth": chosen_depth,
        },
    }

def _read_community_biomarker_summary(report_dir: Path, preview_size: int = 12) -> dict:
    base_dir = report_dir / "microeco_biomarker"
    lefse_raw = _read_tsv_rows(base_dir / "lefse_diff.tsv")
    rf_raw = _read_tsv_rows(base_dir / "rf_importance.tsv")
    assets: list[dict[str, str]] = []
    for label, file_name in [
        ("LEfSe 差异表", "lefse_diff.tsv"),
        ("LEfSe 条形图 PNG", "lefse_barplot.png"),
        ("LEfSe 条形图 PDF", "lefse_barplot.pdf"),
        ("RF 特征重要性", "rf_importance.tsv"),
        ("RF 特征重要性 PNG", "rf_importance.png"),
        ("RF 特征重要性 PDF", "rf_importance.pdf"),
        ("Biomarker 运行摘要", "run_summary.txt"),
    ]:
        asset_path = base_dir / file_name
        if asset_path.exists():
            assets.append({"label": label, "status": "ready", "path": str(asset_path)})
    if not lefse_raw.get("rows") and not rf_raw.get("rows"):
        return {"summary": {}, "columns": [], "rows": [], "assets": assets}

    def _normalize_row_map(raw: dict) -> list[dict[str, str]]:
        columns = [str(col or "").strip() for col in raw.get("columns", [])]
        rows = raw.get("rows", [])
        normalized_rows: list[dict[str, str]] = []
        for row in rows:
            normalized_rows.append({
                columns[index] if index < len(columns) else f"col_{index + 1}": str(row[index] if index < len(row) else "").strip()
                for index in range(len(columns))
            })
        return normalized_rows

    lefse_rows = _normalize_row_map(lefse_raw)
    rf_rows = _normalize_row_map(rf_raw)
    lefse_preview: list[list[str]] = []
    lefse_significant = 0
    lefse_lda_max = None
    for row in lefse_rows:
        taxon = str(row.get("Taxa") or row.get("taxa") or row.get("Taxon") or "").strip()
        group = str(row.get("Group") or row.get("Comparison") or row.get("group") or "").strip() or "--"
        pvalue = _safe_float(row.get("P.unadj") or row.get("P.adj") or row.get("Pvalue") or row.get("P"))
        lda_score = _safe_float(row.get("LDA") or row.get("LDA_score") or row.get("Score") or row.get("Importance"))
        if pvalue is not None and pvalue <= 0.05:
            lefse_significant += 1
        if lda_score is not None:
            lefse_lda_max = max(lefse_lda_max, lda_score) if lefse_lda_max is not None else lda_score
        if taxon and len(lefse_preview) < preview_size:
            lefse_preview.append([
                "LEfSe",
                taxon,
                group,
                f"{pvalue:.4f}" if pvalue is not None else "--",
                f"{lda_score:.4f}" if lda_score is not None else "--",
            ])
    rf_preview: list[list[str]] = []
    rf_top_importance = None
    for row in rf_rows[:preview_size]:
        taxon = str(row.get("Taxa") or row.get("taxa") or row.get("Feature") or "").strip()
        importance = _safe_float(row.get("Importance") or row.get("importance"))
        method = str(row.get("Method") or "RF").strip() or "RF"
        if importance is not None:
            rf_top_importance = max(rf_top_importance, importance) if rf_top_importance is not None else importance
        if taxon:
            rf_preview.append([
                "RF",
                taxon,
                method,
                "--",
                f"{importance:.4f}" if importance is not None else "--",
            ])
    preview_rows = (lefse_preview[: max(1, preview_size // 2)] + rf_preview[: max(1, preview_size - max(1, preview_size // 2))])[:preview_size]
    return {
        "summary": {
            "lefse_feature_count": len(lefse_rows),
            "lefse_significant_count": lefse_significant,
            "lefse_lda_max": round(lefse_lda_max, 4) if lefse_lda_max is not None else None,
            "rf_feature_count": len(rf_rows),
            "rf_top_importance": round(rf_top_importance, 4) if rf_top_importance is not None else None,
            "preview_mode": "lefse_rf",
        },
        "columns": ["方法", "分类单元", "分组/模型", "P值", "效应值/重要性"],
        "rows": preview_rows,
        "lefse": {
            "columns": ["方法", "分类单元", "分组", "P值", "LDA"],
            "rows": lefse_preview[:preview_size],
        },
        "rf": {
            "columns": ["方法", "分类单元", "模型", "P值", "重要性"],
            "rows": rf_preview[:preview_size],
        },
        "assets": assets,
    }

def _read_community_network_summary(report_dir: Path, preview_size: int = 16, graph_node_limit: int = 120) -> dict:
    base_dir = report_dir / "microeco_network"
    empty_result = {
        "summary": {},
        "nodes": [],
        "edges": [],
        "node_preview": {"columns": [], "rows": []},
        "edge_preview": {"columns": [], "rows": []},
        "module_preview": {"columns": [], "rows": []},
        "role_preview": {"columns": [], "rows": []},
        "eigen_preview": {"columns": [], "rows": []},
        "assets": [],
    }
    if not base_dir.is_dir():
        return empty_result

    summary_rows = _read_tsv_rows(base_dir / "network_summary.tsv")
    node_raw = _read_tsv_rows(base_dir / "node_table.tsv")
    edge_raw = _read_tsv_rows(base_dir / "edge_table.tsv")
    module_raw = _read_tsv_rows(base_dir / "module_summary.tsv")
    role_raw = _read_tsv_rows(base_dir / "role_summary.tsv")
    eigen_raw = _read_tsv_rows(base_dir / "eigen_summary.tsv")

    def _rows_to_maps(raw: dict) -> list[dict[str, str]]:
        columns = [str(col or "").strip() for col in raw.get("columns", [])]
        rows = raw.get("rows", [])
        normalized_rows: list[dict[str, str]] = []
        for row in rows:
            normalized_rows.append({
                columns[index] if index < len(columns) else f"col_{index + 1}": str(row[index] if index < len(row) else "").strip()
                for index in range(len(columns))
            })
        return normalized_rows

    summary_map = {
        str(row[0] if len(row) > 0 else "").strip(): str(row[1] if len(row) > 1 else "").strip()
        for row in summary_rows.get("rows", [])
        if row
    }
    node_rows = _rows_to_maps(node_raw)
    edge_rows = _rows_to_maps(edge_raw)
    module_rows = _rows_to_maps(module_raw)
    role_rows = _rows_to_maps(role_raw)
    eigen_rows = _rows_to_maps(eigen_raw)

    node_rows_sorted = sorted(
        node_rows,
        key=lambda item: (_safe_float(item.get("degree")) or 0, _safe_float(item.get("Abundance")) or 0),
        reverse=True,
    )
    graph_node_rows = node_rows_sorted[:graph_node_limit]
    graph_node_names = {str(item.get("name") or "").strip() for item in graph_node_rows if str(item.get("name") or "").strip()}

    nodes = [
        {
            "id": str(item.get("name") or "").strip(),
            "label": str(item.get("Genus") or item.get("name") or "").strip() or str(item.get("name") or "--").strip(),
            "module": str(item.get("module") or "未分模块").strip() or "未分模块",
            "degree": _safe_float(item.get("degree")) or 0.0,
            "abundance": _safe_float(item.get("Abundance")) or 0.0,
            "phylum": str(item.get("Phylum") or "未注释").strip() or "未注释",
            "genus": str(item.get("Genus") or "未注释").strip() or "未注释",
            "role": str(item.get("taxa_roles") or "未分类").strip() or "未分类",
            "z": _safe_float(item.get("z")),
            "p": _safe_float(item.get("p")),
        }
        for item in graph_node_rows
        if str(item.get("name") or "").strip()
    ]
    edges = [
        {
            "source": str(item.get("node1") or "").strip(),
            "target": str(item.get("node2") or "").strip(),
            "label": str(item.get("label") or "").strip() or "+",
            "weight": _safe_float(item.get("weight")) or 0.0,
        }
        for item in edge_rows
        if str(item.get("node1") or "").strip() in graph_node_names and str(item.get("node2") or "").strip() in graph_node_names
    ]

    assets: list[dict[str, str]] = []
    for label, file_name in [
        ("网络摘要", "network_summary.tsv"),
        ("节点属性表", "node_table.tsv"),
        ("边属性表", "edge_table.tsv"),
        ("模块统计", "module_summary.tsv"),
        ("节点角色统计", "role_summary.tsv"),
        ("模块特征向量", "eigen_summary.tsv"),
        ("运行摘要", "run_summary.txt"),
        ("门水平正相关汇总", "phylum_links_positive.tsv"),
        ("门水平负相关汇总", "phylum_links_negative.tsv"),
    ]:
        asset_path = base_dir / file_name
        if asset_path.exists():
            assets.append({"label": label, "status": "ready", "path": str(asset_path)})

    positive_edge_count = _safe_int(summary_map.get("positive_edge_count"))
    negative_edge_count = _safe_int(summary_map.get("negative_edge_count"))
    module_count = _safe_int(summary_map.get("module_count"))
    node_count = _safe_int(summary_map.get("node_count"))
    edge_count = _safe_int(summary_map.get("edge_count"))
    avg_degree = _safe_float(summary_map.get("avg_degree"))
    density = _safe_float(summary_map.get("density"))
    modularity = _safe_float(summary_map.get("modularity"))
    mean_abs_weight = _safe_float(summary_map.get("mean_abs_weight"))

    return {
        "summary": {
            "status": summary_map.get("status") or "",
            "taxa_level": summary_map.get("taxa_level") or "",
            "cor_method": summary_map.get("cor_method") or "",
            "cor_cut": _safe_float(summary_map.get("cor_cut")),
            "p_thres": _safe_float(summary_map.get("p_thres")),
            "filter_thres": _safe_float(summary_map.get("filter_thres")),
            "node_count": node_count,
            "edge_count": edge_count,
            "positive_edge_count": positive_edge_count,
            "negative_edge_count": negative_edge_count,
            "module_count": module_count,
            "avg_degree": round(avg_degree, 4) if avg_degree is not None else None,
            "density": round(density, 6) if density is not None else None,
            "modularity": round(modularity, 6) if modularity is not None else None,
            "mean_abs_weight": round(mean_abs_weight, 4) if mean_abs_weight is not None else None,
            "preview_nodes": len(nodes),
            "preview_edges": len(edges),
        },
        "nodes": nodes,
        "edges": edges,
        "node_preview": {
            "columns": ["节点", "模块", "Degree", "丰度%", "门", "属", "角色"],
            "rows": [
                [
                    str(item.get("name") or "--"),
                    str(item.get("module") or "--"),
                    str(item.get("degree") or "--"),
                    str(item.get("Abundance") or "--"),
                    str(item.get("Phylum") or "--"),
                    str(item.get("Genus") or "--"),
                    str(item.get("taxa_roles") or "--"),
                ]
                for item in node_rows_sorted[:preview_size]
            ],
        },
        "edge_preview": {
            "columns": ["节点 1", "节点 2", "方向", "权重"],
            "rows": [
                [
                    str(item.get("node1") or "--"),
                    str(item.get("node2") or "--"),
                    str(item.get("label") or "--"),
                    str(item.get("weight") or "--"),
                ]
                for item in sorted(edge_rows, key=lambda row: abs(_safe_float(row.get("weight")) or 0), reverse=True)[:preview_size]
            ],
        },
        "module_preview": {
            "columns": ["模块", "节点数", "平均 Degree", "最高 Degree", "总丰度%"],
            "rows": [
                [
                    str(item.get("module") or "--"),
                    str(item.get("node_count") or "--"),
                    str(item.get("mean_degree") or "--"),
                    str(item.get("max_degree") or "--"),
                    str(item.get("total_abundance") or "--"),
                ]
                for item in module_rows[:preview_size]
            ],
        },
        "role_preview": {
            "columns": ["节点角色", "节点数", "平均 Degree"],
            "rows": [
                [
                    str(item.get("taxa_role") or "--"),
                    str(item.get("node_count") or "--"),
                    str(item.get("mean_degree") or "--"),
                ]
                for item in role_rows[:preview_size]
            ],
        },
        "eigen_preview": {
            "columns": ["模块", "解释度/特征值"],
            "rows": [
                [
                    str(item.get("module") or item.get("row.names") or "--"),
                    str(item.get("variance_explained") or item.get("PC1") or item.get("value") or "--"),
                ]
                for item in eigen_rows[:preview_size]
            ],
        },
        "assets": assets,
    }

def _read_community_beta_summary(report_dir: Path, group_column: str = "SampleType", preview_size: int = 12) -> dict:
    base_dir = report_dir / "microeco_beta"
    empty_result = {
        "summary": {},
        "assets": [],
        "pcoa_points": [],
        "pcoa_preview": {"columns": [], "rows": []},
        "nmds_points": [],
        "nmds_preview": {"columns": [], "rows": []},
        "distance_matrix": {"samples": [], "groups": {}, "rows": [], "min": None, "max": None},
        "permanova": {"columns": [], "rows": []},
        "anosim": {"columns": [], "rows": []},
        "group_distances": {"columns": [], "rows": []},
        "group_counts": {"columns": [], "rows": []},
        "dispersion": {"summary": {}, "lines": []},
    }
    if not base_dir.is_dir():
        return empty_result

    def _read_named_rows(path: Path, rename_map: dict[str, str] | None = None, limit: int = preview_size) -> dict:
        raw = _read_tsv_rows(path)
        columns = [str(col or "").strip() for col in raw.get("columns", [])]
        rows = raw.get("rows", [])
        if not columns:
            return {"columns": [], "rows": []}
        if columns and not columns[0]:
            columns[0] = "index"
        mapped_columns = [rename_map.get(column, column) if rename_map else column for column in columns]
        trimmed_rows = []
        for row in rows[:limit]:
            trimmed_rows.append([str(row[index] if index < len(row) else "").strip() for index in range(len(columns))])
        return {"columns": mapped_columns, "rows": trimmed_rows}

    pcoa_raw = _read_tsv_rows(base_dir / "pcoa_scores.tsv")
    pcoa_columns = [str(col or "").strip() for col in pcoa_raw.get("columns", [])]
    pcoa_rows = pcoa_raw.get("rows", [])
    pcoa_preview = {"columns": [], "rows": []}
    pcoa_points: list[dict[str, object]] = []
    sample_count = 0
    group_values: set[str] = set()
    if pcoa_columns and pcoa_rows:
        if not pcoa_columns[0]:
            pcoa_columns[0] = "SampleID"
        column_lookup = {column: index for index, column in enumerate(pcoa_columns)}
        pco1_index = column_lookup.get("PCo1")
        pco2_index = column_lookup.get("PCo2")
        sample_index = column_lookup.get("SampleID", 0)
        group_index = column_lookup.get(group_column)
        if pco1_index is not None and pco2_index is not None:
            preview_rows: list[list[str]] = []
            for row in pcoa_rows:
                sample_name = str(row[sample_index] if sample_index < len(row) else "").strip()
                group_name = str(row[group_index] if group_index is not None and group_index < len(row) else "未分组").strip() or "未分组"
                if sample_name:
                    sample_count += 1
                if group_name:
                    group_values.add(group_name)
                pco1_value = _safe_float(row[pco1_index] if pco1_index < len(row) else None)
                pco2_value = _safe_float(row[pco2_index] if pco2_index < len(row) else None)
                if pco1_value is not None and pco2_value is not None:
                    pcoa_points.append({
                        "sample": sample_name or "--",
                        "group": group_name,
                        "x": round(pco1_value, 6),
                        "y": round(pco2_value, 6),
                    })
                preview_rows.append([
                    sample_name or "--",
                    group_name,
                    str(row[pco1_index] if pco1_index < len(row) else "").strip(),
                    str(row[pco2_index] if pco2_index < len(row) else "").strip(),
                ])
            pcoa_preview = {"columns": ["样本 ID", "分组", "PCo1", "PCo2"], "rows": preview_rows[:preview_size]}

    nmds_raw = _read_tsv_rows(base_dir / "nmds_scores.tsv")
    nmds_columns = [str(col or "").strip() for col in nmds_raw.get("columns", [])]
    nmds_rows = nmds_raw.get("rows", [])
    nmds_preview = {"columns": [], "rows": []}
    nmds_points: list[dict[str, object]] = []
    nmds_stress = None
    if nmds_columns and nmds_rows:
        if not nmds_columns[0]:
            nmds_columns[0] = "SampleID"
        column_lookup = {column: index for index, column in enumerate(nmds_columns)}
        axis1_index = column_lookup.get("MDS1", column_lookup.get("NMDS1"))
        axis2_index = column_lookup.get("MDS2", column_lookup.get("NMDS2"))
        sample_index = column_lookup.get("SampleID", 0)
        group_index = column_lookup.get(group_column)
        if "Stress" in nmds_columns:
            stress_index = column_lookup.get("Stress")
            if stress_index is not None and nmds_rows:
                nmds_stress = _safe_float(nmds_rows[0][stress_index] if stress_index < len(nmds_rows[0]) else None)
        if axis1_index is not None and axis2_index is not None:
            preview_rows: list[list[str]] = []
            for row in nmds_rows:
                sample_name = str(row[sample_index] if sample_index < len(row) else "").strip()
                group_name = str(row[group_index] if group_index is not None and group_index < len(row) else "未分组").strip() or "未分组"
                axis1_value = _safe_float(row[axis1_index] if axis1_index < len(row) else None)
                axis2_value = _safe_float(row[axis2_index] if axis2_index < len(row) else None)
                if axis1_value is not None and axis2_value is not None:
                    nmds_points.append({
                        "sample": sample_name or "--",
                        "group": group_name,
                        "x": round(axis1_value, 6),
                        "y": round(axis2_value, 6),
                    })
                preview_rows.append([
                    sample_name or "--",
                    group_name,
                    str(row[axis1_index] if axis1_index < len(row) else "").strip(),
                    str(row[axis2_index] if axis2_index < len(row) else "").strip(),
                ])
            nmds_preview = {"columns": ["样本 ID", "分组", "NMDS1", "NMDS2"], "rows": preview_rows[:preview_size]}
    if nmds_stress is None:
        run_summary_path = base_dir / "run_summary.txt"
        if run_summary_path.is_file():
            try:
                for line in run_summary_path.read_text(encoding="utf-8", errors="ignore").splitlines():
                    if not str(line).startswith("nmds_stress:"):
                        continue
                    nmds_stress = _safe_float(str(line).split(":", 1)[1].strip())
                    break
            except OSError:
                pass

    permanova = _read_named_rows(
        base_dir / "permanova.tsv",
        {
            "index": "来源",
            "Df": "自由度",
            "SumOfSqs": "平方和",
            "R2": "R2",
            "F": "F值",
            "Pr(>F)": "P值",
            "Significance": "显著性",
        },
    )
    anosim = _read_named_rows(
        base_dir / "anosim.tsv",
        {
            "index": "序号",
            "Test": "检验",
            "permutations": "置换次数",
            "statistic.R": "R统计量",
            "p.value": "P值",
            "Significance": "显著性",
        },
    )
    group_distances = _read_named_rows(
        base_dir / "within_group_distance_stats.tsv",
        {
            "index": "序号",
            "Comparison": "比较",
            "Measure": "指标",
            "Method": "方法",
            "Group": "分组",
            "P.unadj": "未校正P值",
            "P.adj": "校正P值",
            "Significance": "显著性",
        },
    )
    group_counts = _read_named_rows(
        base_dir / "group_sample_counts.tsv",
        {
            group_column: "分组",
            "sample_count": "样本数",
        },
        limit=preview_size * 2,
    )
    distance_matrix = {"samples": [], "groups": {}, "rows": [], "min": None, "max": None}
    within_group_raw = _read_tsv_rows(base_dir / "within_group_distances.tsv")
    within_group_columns = [str(col or "").strip() for col in within_group_raw.get("columns", [])]
    within_group_rows = within_group_raw.get("rows", [])
    if within_group_columns and within_group_rows:
        if within_group_columns and not within_group_columns[0]:
            within_group_columns[0] = "index"
        column_lookup = {column: index for index, column in enumerate(within_group_columns)}
        value_index = column_lookup.get("Value")
        group_index = column_lookup.get(group_column)
        grouped_values: dict[str, list[float]] = {}
        if value_index is not None and group_index is not None:
            for row in within_group_rows:
                group_name = str(row[group_index] if group_index < len(row) else "").strip() or "未分组"
                value = _safe_float(row[value_index] if value_index < len(row) else None)
                if value is None:
                    continue
                grouped_values.setdefault(group_name, []).append(float(value))
    matrix_raw = _read_tsv_rows(base_dir / "bray_distance_matrix.tsv")
    matrix_columns = [str(col or "").strip() for col in matrix_raw.get("columns", [])]
    matrix_rows = matrix_raw.get("rows", [])
    if matrix_columns and matrix_rows:
        if matrix_columns and not matrix_columns[0]:
            matrix_columns = matrix_columns[1:]
        samples = [str(value or "").strip() for value in matrix_columns if str(value or "").strip()]
        pcoa_group_lookup = {str(point.get("sample") or ""): str(point.get("group") or "未分组") for point in pcoa_points}
        matrix_values: list[list[float]] = []
        valid_values: list[float] = []
        for row in matrix_rows:
            numeric_row: list[float] = []
            for index in range(len(samples)):
                source_index = index + 1
                value = _safe_float(row[source_index] if source_index < len(row) else None)
                numeric = float(value) if value is not None else 0.0
                numeric_row.append(round(numeric, 6))
                valid_values.append(numeric)
            matrix_values.append(numeric_row)
        distance_matrix = {
            "samples": samples,
            "groups": {sample: pcoa_group_lookup.get(sample, "未分组") for sample in samples},
            "rows": matrix_values,
            "min": round(min(valid_values), 6) if valid_values else None,
            "max": round(max(valid_values), 6) if valid_values else None,
        }

    permanova_r2 = None
    permanova_p = None
    if permanova["rows"]:
        first_row = permanova["rows"][0]
        if "R2" in permanova["columns"]:
            permanova_r2 = _safe_float(first_row[permanova["columns"].index("R2")])
        if "P值" in permanova["columns"]:
            permanova_p = _safe_float(first_row[permanova["columns"].index("P值")])

    anosim_r = None
    anosim_p = None
    if anosim["rows"]:
        first_row = anosim["rows"][0]
        if "R统计量" in anosim["columns"]:
            anosim_r = _safe_float(first_row[anosim["columns"].index("R统计量")])
        if "P值" in anosim["columns"]:
            anosim_p = _safe_float(first_row[anosim["columns"].index("P值")])

    betadisper_summary: dict[str, object] = {}
    betadisper_lines: list[str] = []
    betadisper_path = base_dir / "betadisper.txt"
    if betadisper_path.is_file():
        try:
            raw_text = betadisper_path.read_text(encoding="utf-8", errors="ignore")
            all_lines = [str(line or "").rstrip() for line in raw_text.splitlines()]
            betadisper_lines = [line for line in all_lines if line.strip()][:preview_size]
            for line in all_lines:
                stripped = line.strip()
                if not stripped.startswith("Groups"):
                    continue
                parts = re.split(r"\s+", stripped)
                if len(parts) >= 7:
                    betadisper_summary = {
                        "f_value": _safe_float(parts[4]),
                        "p_value": _safe_float(parts[6]),
                        "significance": parts[7] if len(parts) >= 8 else "",
                    }
                    break
        except OSError:
            betadisper_lines = []

    assets: list[dict[str, str]] = []
    for label, file_name in [
        ("PCoA 排序图", "pcoa_plot.png"),
        ("NMDS 排序图", "nmds_plot.png"),
        ("组内距离分布图", "within_group_distance_plot.png"),
        ("PCoA 排序图 PDF", "pcoa_plot.pdf"),
        ("NMDS 排序图 PDF", "nmds_plot.pdf"),
        ("组内距离分布图 PDF", "within_group_distance_plot.pdf"),
        ("Bray-Curtis 距离矩阵", "bray_distance_matrix.tsv"),
        ("PCoA 坐标表", "pcoa_scores.tsv"),
        ("NMDS 坐标表", "nmds_scores.tsv"),
        ("PERMANOVA 结果", "permanova.tsv"),
        ("ANOSIM 结果", "anosim.tsv"),
        ("Betadisper 结果", "betadisper.txt"),
        ("组内距离统计", "within_group_distance_stats.tsv"),
    ]:
        asset_path = base_dir / file_name
        if asset_path.exists():
            assets.append({"label": label, "status": "ready", "path": str(asset_path)})

    return {
        "summary": {
            "sample_count": sample_count or None,
            "group_count": len(group_values) or None,
            "measure": "bray",
            "permanova_r2": round(permanova_r2, 4) if permanova_r2 is not None else None,
            "permanova_p": permanova_p,
            "anosim_r": round(anosim_r, 4) if anosim_r is not None else None,
            "anosim_p": anosim_p,
            "nmds_stress": round(nmds_stress, 4) if nmds_stress is not None else None,
            "betadisper_f": betadisper_summary.get("f_value"),
            "betadisper_p": betadisper_summary.get("p_value"),
            "betadisper_significance": betadisper_summary.get("significance") or "",
        },
        "assets": assets,
        "pcoa_points": pcoa_points,
        "pcoa_preview": pcoa_preview,
        "nmds_points": nmds_points,
        "nmds_preview": nmds_preview,
        "distance_matrix": distance_matrix,
        "permanova": permanova,
        "anosim": anosim,
        "group_distances": group_distances,
        "group_counts": group_counts,
        "dispersion": {"summary": betadisper_summary, "lines": betadisper_lines},
    }

def _build_community_report_payload(*, task: dict, report_dir: Path, report_source: dict) -> dict:
    summary = _read_community_summary(report_dir / "community_summary.json")
    workflow_mode = str(summary.get("workflow_mode") or "abundance").strip()
    metadata_summary = summary.get("metadata_summary") or {}
    taxonomy_summary = summary.get("taxonomy_summary") or {}
    demux_summary = summary.get("demux_summary") or {}
    demux_preview = _read_community_demux_preview(report_dir)
    denoise_preview = _read_community_denoise_preview(report_dir)
    taxonomy_preview = _read_community_taxonomy_preview(report_dir)
    differential_summary = _read_community_biomarker_summary(report_dir)
    input_summary = summary.get("input_summary") or {}
    parameters = summary.get("parameters") or {}
    modules = summary.get("modules") if isinstance(summary.get("modules"), list) else []
    commands = summary.get("commands") if isinstance(summary.get("commands"), list) else []
    outputs = summary.get("outputs") if isinstance(summary.get("outputs"), list) else []
    task_params = task.get("params") or {}
    metadata_path_text = str(task_params.get("metadata") or "").strip()
    metadata_lookup = _read_metadata_group_lookup(
        Path(metadata_path_text),
        str(metadata_summary.get("sample_id_column") or "sample-id"),
        str(metadata_summary.get("group_column") or "SampleType"),
    ) if metadata_path_text else {}
    taxa_abundance = _read_community_taxa_abundance_summary(report_dir, metadata_groups=metadata_lookup)
    alpha_summary = _read_community_alpha_rarefaction_summary(
        report_dir,
        metadata_groups=metadata_lookup,
        preferred_depth=int(_safe_float(demux_summary.get("suggested_sampling_depth")) or 0) or None,
    )
    beta_summary = _read_community_beta_summary(
        report_dir,
        group_column=str(metadata_summary.get("group_column") or task_params.get("group_column") or "SampleType"),
    )
    network_summary = _read_community_network_summary(report_dir)
    known_outputs = [
        ("demux.qzv", report_dir / "demux.qzv"),
        ("denoising-stats-dada2.qzv", report_dir / "denoising-stats-dada2.qzv"),
        ("taxonomy.qzv", report_dir / "taxonomy.qzv"),
        ("taxa-barplot.qzv", report_dir / "taxa-barplot.qzv"),
        ("alpha-rarefaction.qzv", report_dir / "alpha-rarefaction.qzv"),
        ("ancombc2-sampletype.qzv", report_dir / "ancombc2-sampletype.qzv"),
        ("taxonomy.tsv", report_dir / "taxonomy_export" / "taxonomy.tsv"),
        ("ancombc2_export", report_dir / "ancombc2_export"),
        ("microeco_beta", report_dir / "microeco_beta"),
        ("microeco_biomarker", report_dir / "microeco_biomarker"),
        ("microeco_network", report_dir / "microeco_network"),
        ("lefse_diff.tsv", report_dir / "microeco_biomarker" / "lefse_diff.tsv"),
        ("lefse_barplot.png", report_dir / "microeco_biomarker" / "lefse_barplot.png"),
        ("lefse_barplot.pdf", report_dir / "microeco_biomarker" / "lefse_barplot.pdf"),
        ("rf_importance.tsv", report_dir / "microeco_biomarker" / "rf_importance.tsv"),
        ("rf_importance.png", report_dir / "microeco_biomarker" / "rf_importance.png"),
        ("rf_importance.pdf", report_dir / "microeco_biomarker" / "rf_importance.pdf"),
        ("network_summary.tsv", report_dir / "microeco_network" / "network_summary.tsv"),
        ("node_table.tsv", report_dir / "microeco_network" / "node_table.tsv"),
        ("edge_table.tsv", report_dir / "microeco_network" / "edge_table.tsv"),
        ("module_summary.tsv", report_dir / "microeco_network" / "module_summary.tsv"),
    ]
    for metric_qzv in sorted(report_dir.glob("alpha-*.qzv")):
        known_outputs.append((metric_qzv.name, metric_qzv))
    for metric_export in sorted(report_dir.glob("alpha-*_export")):
        known_outputs.append((metric_export.name, metric_export))
    output_map: dict[str, dict] = {}
    for item in outputs:
        if not isinstance(item, dict):
            continue
        path_value = str(item.get("path") or "").strip()
        label_value = str(item.get("label") or "").strip() or (Path(path_value).name if path_value else "")
        status_value = str(item.get("status") or "").strip().lower()
        candidate_path = Path(path_value) if path_value else None
        if candidate_path and candidate_path.exists():
            status_value = "ready"
        elif status_value in {"planned", "needs-demux-artifact", ""}:
            status_value = "missing"
        normalized = {
            "label": label_value or "--",
            "status": status_value or "missing",
            "path": path_value,
        }
        output_map[normalized["label"]] = normalized
    for label, output_path in known_outputs:
        output_map[label] = {
            "label": label,
            "status": "ready" if output_path.exists() else "missing",
            "path": str(output_path),
        }
    outputs = list(output_map.values())
    ready_outputs = sum(1 for item in outputs if str(item.get("status") or "").strip().lower() == "ready")
    return {
        "task": {
            "id": task.get("id"),
            "name": task.get("name"),
            "status": task.get("status"),
            "owner": task.get("owner"),
            "group": task.get("owner_group", ""),
            "created_at": task.get("created_at"),
            "started_at": task.get("started_at"),
            "finished_at": task.get("finished_at"),
            "input_path": task_params.get("input_path", ""),
            "output_dir": task_params.get("output_dir", ""),
            "asm_type": "",
            "method": "community",
            "analysis_target": "bacteria",
            "species": "",
            "sample_name": "",
            "sample_display_name": task.get("name", ""),
            "samples": [],
            "report_mode": "cohort",
            "report_kind": "community_meta_ecology",
            "workstation_key": "community",
        },
        "overview_metrics": [
            {"key": "community_samples", "label": "纳入样本", "type": "single", "display": str(metadata_summary.get("sample_count") or "--")},
            {"key": "community_metadata_columns", "label": "元数据字段", "type": "single", "display": str(len(metadata_summary.get("metadata_columns") or []) or "--")},
            {"key": "community_taxa", "label": "Taxonomy 条目", "type": "single", "display": str(taxonomy_summary.get("feature_count") or "--")},
            {"key": "community_modules", "label": "分析模块", "type": "single", "display": str(len(modules) or "--")},
            {"key": "community_outputs", "label": "就绪输出", "type": "single", "display": str(ready_outputs or "--")},
        ],
        "sections": {
            "community": {
                "summary": {
                    "sample_count": metadata_summary.get("sample_count"),
                    "group_column": metadata_summary.get("group_column", ""),
                    "sample_id_column": metadata_summary.get("sample_id_column", ""),
                    "input_type": input_summary.get("path_type", ""),
                    "input_entry_count": input_summary.get("entry_count"),
                    "workflow_mode": workflow_mode,
                    "taxonomy_feature_count": taxonomy_summary.get("feature_count"),
                    "taxonomy_feature_column": taxonomy_summary.get("feature_column", ""),
                    "demux_sample_count": demux_summary.get("sample_count"),
                    "demux_trunc_len_f": demux_summary.get("suggested_trunc_len_f"),
                    "demux_trunc_len_r": demux_summary.get("suggested_trunc_len_r"),
                    "demux_sampling_depth": demux_summary.get("suggested_sampling_depth"),
                    "taxonomy_level": parameters.get("taxonomy_level", ""),
                    "normalization": parameters.get("normalization", ""),
                },
                "metadata": {
                    "columns": ["字段", "值"],
                    "rows": [
                        ["流程模式", "QIIME2 amplicon" if workflow_mode == "amplicon" else "丰度表模式"],
                        ["样本 ID 列", str(metadata_summary.get("sample_id_column") or "--")],
                        ["分组列", str(metadata_summary.get("group_column") or "--")],
                        ["样本数", str(metadata_summary.get("sample_count") or "--")],
                        ["元数据字段", ", ".join(metadata_summary.get("metadata_columns") or []) or "--"],
                        ["输入路径类型", str(input_summary.get("path_type") or "--")],
                        ["输入条目数", str(input_summary.get("entry_count") or "--")],
                        ["taxonomy 条目数", str(taxonomy_summary.get("feature_count") or "--")],
                        ["taxonomy 主键列", str(taxonomy_summary.get("feature_column") or "--")],
                        ["demux 样本数", str(demux_summary.get("sample_count") or "--")],
                        ["建议 trunc-len-f / r", f"{demux_summary.get('suggested_trunc_len_f') or '--'} / {demux_summary.get('suggested_trunc_len_r') or '--'}"],
                        ["建议 rarefaction depth", str(demux_summary.get("suggested_sampling_depth") or "--")],
                        ["标准化方式", str(parameters.get("normalization") or "--")],
                        ["统计层级", str(parameters.get("taxonomy_level") or "--")],
                    ],
                },
                "demux_preview": {
                    "summary": demux_preview.get("summary") or {},
                    "columns": demux_preview.get("columns") or [],
                    "rows": demux_preview.get("rows") or [],
                },
                "denoise_preview": {
                    "summary": denoise_preview.get("summary") or {},
                    "columns": denoise_preview.get("columns") or [],
                    "rows": denoise_preview.get("rows") or [],
                },
                "modules": {
                    "columns": ["模块", "状态", "运行环境", "说明"],
                    "rows": [
                        [
                            str(item.get("label") or "--"),
                            str(item.get("status") or "--"),
                            str(item.get("runtime") or "--"),
                            str(item.get("description") or "--"),
                        ]
                        for item in modules
                    ],
                },
                "commands": {
                    "columns": ["模块", "运行环境", "命令"],
                    "rows": [
                        [
                            str(item.get("module") or "--"),
                            str(item.get("runtime") or "--"),
                            str(item.get("command") or "--"),
                        ]
                        for item in commands
                    ],
                },
                "outputs": {
                    "columns": ["输出项", "状态", "路径"],
                    "rows": [
                        [
                            str(item.get("label") or "--"),
                            str(item.get("status") or "--"),
                            str(item.get("path") or "--"),
                        ]
                        for item in outputs
                    ],
                },
                "taxonomy_preview": {
                    "summary": taxonomy_preview.get("summary") or {},
                    "columns": taxonomy_preview.get("columns") or [],
                    "rows": taxonomy_preview.get("rows") or [],
                },
                "taxa_abundance": {
                    "summary": taxa_abundance.get("summary") or {},
                    "columns": taxa_abundance.get("columns") or [],
                    "rows": taxa_abundance.get("rows") or [],
                    "levels": taxa_abundance.get("levels") or {},
                },
                "differential": {
                    "summary": differential_summary.get("summary") or {},
                    "columns": differential_summary.get("columns") or [],
                    "rows": differential_summary.get("rows") or [],
                    "lefse": differential_summary.get("lefse") or {"columns": [], "rows": []},
                    "rf": differential_summary.get("rf") or {"columns": [], "rows": []},
                    "assets": differential_summary.get("assets") or [],
                },
                "alpha": {
                    "summary": alpha_summary.get("summary") or {},
                    "sample_columns": alpha_summary.get("sample_columns") or [],
                    "sample_rows": alpha_summary.get("sample_rows") or [],
                    "group_columns": alpha_summary.get("group_columns") or [],
                    "group_rows": alpha_summary.get("group_rows") or [],
                    "boxplots": alpha_summary.get("boxplots") or {},
                    "pairwise": alpha_summary.get("pairwise") or {},
                    "rarefaction": alpha_summary.get("rarefaction") or {},
                },
                "beta": {
                    "summary": beta_summary.get("summary") or {},
                    "assets": beta_summary.get("assets") or [],
                    "pcoa_points": beta_summary.get("pcoa_points") or [],
                    "pcoa_preview": beta_summary.get("pcoa_preview") or {"columns": [], "rows": []},
                    "nmds_points": beta_summary.get("nmds_points") or [],
                    "nmds_preview": beta_summary.get("nmds_preview") or {"columns": [], "rows": []},
                    "distance_matrix": beta_summary.get("distance_matrix") or {"samples": [], "groups": {}, "rows": [], "min": None, "max": None},
                    "permanova": beta_summary.get("permanova") or {"columns": [], "rows": []},
                    "anosim": beta_summary.get("anosim") or {"columns": [], "rows": []},
                    "group_distances": beta_summary.get("group_distances") or {"columns": [], "rows": []},
                    "group_counts": beta_summary.get("group_counts") or {"columns": [], "rows": []},
                    "dispersion": beta_summary.get("dispersion") or {"summary": {}, "lines": []},
                },
                "network": {
                    "summary": network_summary.get("summary") or {},
                    "nodes": network_summary.get("nodes") or [],
                    "edges": network_summary.get("edges") or [],
                    "node_preview": network_summary.get("node_preview") or {"columns": [], "rows": []},
                    "edge_preview": network_summary.get("edge_preview") or {"columns": [], "rows": []},
                    "module_preview": network_summary.get("module_preview") or {"columns": [], "rows": []},
                    "role_preview": network_summary.get("role_preview") or {"columns": [], "rows": []},
                    "eigen_preview": network_summary.get("eigen_preview") or {"columns": [], "rows": []},
                    "assets": network_summary.get("assets") or [],
                },
                "notes": summary.get("notes") if isinstance(summary.get("notes"), list) else [],
                "report_source": {
                    key: str(value) if isinstance(value, Path) else value
                    for key, value in report_source.items()
                },
            },
        },
    }
