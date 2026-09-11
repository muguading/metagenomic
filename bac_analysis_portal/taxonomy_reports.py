from __future__ import annotations

import math
import re
import uuid
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any

from .knowledge_interpretation import _lookup_kb_taxonomy
from .parse_utils import _safe_float, _safe_int
from .report_sources import _resolve_report_sample_name, _resolve_report_source
from .table_io import _read_tsv_rows, write_tsv
from .task_manager import ValidationError, write_json

def _read_taxonomy_list(path: Path, terminal_column: str, taxonomy_index: dict[str, list[dict]] | None = None) -> dict:
    raw = _read_tsv_rows(path)
    if not raw["columns"] or not raw["rows"]:
        return {"rows": [], "rank_options": []}
    rank_options = [column for column in ["界", "门", "纲", "目", "科", "属"] if column in raw["columns"]]
    rows = []
    for row in raw["rows"]:
        record = {raw["columns"][index]: row[index] if index < len(row) else "" for index in range(len(raw["columns"]))}
        record["比例数值"] = _safe_float(record.get("比例")) or 0.0
        record["序列数量数值"] = _safe_int(record.get("序列数量")) or 0
        matched_taxonomy = _lookup_kb_taxonomy(
            taxonomy_index,
            species_name=str(record.get("种") or record.get(terminal_column) or "").strip(),
            genus_name=str(record.get("属") or "").strip(),
        )
        if matched_taxonomy:
            record["NCBI TaxID"] = matched_taxonomy.get("taxid") or "-"
            record["NCBI学名"] = matched_taxonomy.get("scientific_name") or "-"
            record["NCBI分类等级"] = matched_taxonomy.get("rank") or "-"
            record["NCBI目"] = matched_taxonomy.get("order") or "-"
            record["NCBI科"] = matched_taxonomy.get("family") or "-"
            record["NCBI属"] = matched_taxonomy.get("genus") or "-"
            record["NCBI种"] = matched_taxonomy.get("species_rank") or "-"
        else:
            record["NCBI TaxID"] = "-"
            record["NCBI学名"] = "-"
            record["NCBI分类等级"] = "-"
            record["NCBI目"] = "-"
            record["NCBI科"] = "-"
            record["NCBI属"] = "-"
            record["NCBI种"] = "-"
        rows.append(record)
    return {"rows": rows, "rank_options": rank_options, "terminal_column": terminal_column}

def _parse_community_merge_items(raw_value: object) -> list[dict[str, str]]:
    items: list[dict[str, str]] = []
    if not isinstance(raw_value, list):
        return items
    for item in raw_value:
        if not isinstance(item, dict):
            continue
        task_id = str(item.get("task_id") or "").strip()
        group = str(item.get("group") or "").strip() or "未分组"
        if task_id:
            items.append({"task_id": task_id, "group": group})
    return items

def _is_metagenome_task_record(task: dict[str, Any]) -> bool:
    params = task.get("params") if isinstance(task, dict) else {}
    workstation_key = str((params or {}).get("workstation_key") or "").strip().lower()
    method = str((params or {}).get("method") or "").strip().lower()
    return workstation_key == "metagenome" or method == "meta"

def _community_merge_sample_id(base_text: str, seen: set[str], index: int) -> str:
    normalized = re.sub(r"[^\w.\-]+", "_", str(base_text or "").strip(), flags=re.UNICODE).strip("._")
    if not normalized:
        normalized = f"sample_{index}"
    candidate = normalized
    suffix = 2
    while candidate in seen:
        candidate = f"{normalized}_{suffix}"
        suffix += 1
    seen.add(candidate)
    return candidate

def _community_taxonomy_value(row: dict[str, Any], *keys: str) -> str:
    for key in keys:
        value = str(row.get(key) or "").strip()
        if value and value != "-":
            return value
    return ""

def _build_community_merge_bundle(
    *,
    project_root: Path,
    task_name: str,
    group_column: str,
    merge_entries: list[dict[str, Any]],
) -> dict[str, Any]:
    if len(merge_entries) < 2:
        raise ValidationError("群落分析的宏基因组任务汇总模式至少需要选择 2 个任务")
    bundle_root = project_root / "tmp" / "community_merge_inputs"
    bundle_root.mkdir(parents=True, exist_ok=True)
    bundle_name = _community_merge_sample_id(task_name or "community_merge", set(), 1)
    bundle_dir = bundle_root / f"{bundle_name}_{datetime.now().strftime('%Y%m%d%H%M%S')}_{uuid.uuid4().hex[:6]}"
    bundle_dir.mkdir(parents=True, exist_ok=True)

    metadata_path = bundle_dir / "community_metadata.tsv"
    taxonomy_path = bundle_dir / "community_taxonomy.tsv"
    abundance_path = bundle_dir / "community_abundance.tsv"
    manifest_path = bundle_dir / "community_merge_manifest.json"

    seen_sample_ids: set[str] = set()
    abundance_by_sample: dict[str, dict[str, float]] = {}
    taxonomy_lookup: dict[str, dict[str, str]] = {}
    species_totals: defaultdict[str, float] = defaultdict(float)
    metadata_rows: list[list[str]] = []
    manifest_items: list[dict[str, Any]] = []

    for index, entry in enumerate(merge_entries, start=1):
        task = entry["task"]
        task_id = str(task.get("id") or "").strip()
        task_display_name = str(task.get("name") or task_id or f"meta_{index}").strip()
        report_source = _resolve_report_source(task)
        if not report_source.get("available"):
            raise ValidationError(str(report_source.get("reason") or f"宏基因组任务 {task_display_name} 尚未定位到服务器结果目录"))
        report_dir = Path(report_source["report_dir"]).resolve()
        sample_name = str(report_source.get("selected_sample") or _resolve_report_sample_name(task, report_dir) or task_display_name).strip()
        species_payload = _read_taxonomy_list(report_dir / f"{sample_name}_2.list.txt", terminal_column="种")
        taxonomy_rows = species_payload.get("rows") if isinstance(species_payload, dict) else []
        if not taxonomy_rows:
            raise ValidationError(f"宏基因组任务 {task_display_name} 尚未生成可用于群落汇总的物种列表")

        sample_id = _community_merge_sample_id(sample_name or task_display_name, seen_sample_ids, index)
        sample_weights: defaultdict[str, float] = defaultdict(float)
        for row in taxonomy_rows:
            if not isinstance(row, dict):
                continue
            species_name = _community_taxonomy_value(row, "种", "Species", "species")
            if not species_name:
                continue
            reads = int(row.get("序列数量数值") or 0)
            ratio = float(row.get("比例数值") or 0.0)
            weight = float(reads) if reads > 0 else max(ratio, 0.0)
            if weight <= 0:
                continue
            sample_weights[species_name] += weight
            if species_name not in taxonomy_lookup:
                taxonomy_lookup[species_name] = {
                    "Taxon": species_name,
                    "Kingdom": _community_taxonomy_value(row, "界", "Kingdom", "kingdom"),
                    "Phylum": _community_taxonomy_value(row, "门", "Phylum", "phylum"),
                    "Class": _community_taxonomy_value(row, "纲", "Class", "class"),
                    "Order": _community_taxonomy_value(row, "目", "Order", "order"),
                    "Family": _community_taxonomy_value(row, "科", "Family", "family"),
                    "Genus": _community_taxonomy_value(row, "属", "Genus", "genus"),
                    "Species": species_name,
                    "TaxID": _community_taxonomy_value(row, "NCBI TaxID", "TaxID", "taxid"),
                }
        if not sample_weights:
            raise ValidationError(f"宏基因组任务 {task_display_name} 的物种列表为空，无法用于群落汇总")

        total_weight = sum(sample_weights.values())
        normalized_weights = {
            species_name: round(weight / total_weight, 6)
            for species_name, weight in sample_weights.items()
        }
        abundance_by_sample[sample_id] = normalized_weights
        for species_name, weight in normalized_weights.items():
            species_totals[species_name] += weight

        group_value = str(entry.get("group") or "").strip() or "未分组"
        metadata_rows.append([sample_id, group_value, task_display_name, task_id])
        manifest_items.append(
            {
                "task_id": task_id,
                "task_name": task_display_name,
                "sample_id": sample_id,
                "group": group_value,
                "report_dir": str(report_dir),
                "source_species_count": len(sample_weights),
            }
        )

    ordered_species = sorted(species_totals.keys(), key=lambda key: (-species_totals[key], key.lower()))
    ordered_sample_ids = [row[0] for row in metadata_rows]
    write_tsv(
        abundance_path,
        ["Taxon", *ordered_sample_ids],
        [
            [species_name, *[f"{abundance_by_sample.get(sample_id, {}).get(species_name, 0.0):.6f}" for sample_id in ordered_sample_ids]]
            for species_name in ordered_species
        ],
    )
    write_tsv(
        taxonomy_path,
        ["Taxon", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "TaxID"],
        [
            [
                taxonomy_lookup[species_name].get("Taxon", ""),
                taxonomy_lookup[species_name].get("Kingdom", ""),
                taxonomy_lookup[species_name].get("Phylum", ""),
                taxonomy_lookup[species_name].get("Class", ""),
                taxonomy_lookup[species_name].get("Order", ""),
                taxonomy_lookup[species_name].get("Family", ""),
                taxonomy_lookup[species_name].get("Genus", ""),
                taxonomy_lookup[species_name].get("Species", ""),
                taxonomy_lookup[species_name].get("TaxID", ""),
            ]
            for species_name in ordered_species
        ],
    )
    write_tsv(
        metadata_path,
        ["sample-id", group_column, "SourceTask", "SourceTaskId"],
        metadata_rows,
    )
    write_json(
        manifest_path,
        {
            "created_at": datetime.now().isoformat(timespec="seconds"),
            "task_name": task_name,
            "group_column": group_column,
            "generated_input_path": str(bundle_dir),
            "generated_abundance_path": str(abundance_path),
            "generated_taxonomy_path": str(taxonomy_path),
            "generated_metadata_path": str(metadata_path),
            "items": manifest_items,
        },
    )
    return {
        "input_path": str(bundle_dir),
        "community_metadata": str(metadata_path),
        "community_taxonomy": str(taxonomy_path),
        "community_merge_tasks": [{"task_id": item["task_id"], "group": item["group"]} for item in manifest_items],
        "community_generated_abundance": str(abundance_path),
        "community_generated_manifest": str(manifest_path),
    }

def _build_taxonomy_abundance(species_taxonomy: dict, subspecies_taxonomy: dict) -> dict:
    rank_order = ["界", "门", "纲", "目", "科", "属", "种", "亚种"]
    rank_sources = {
        "界": species_taxonomy,
        "门": species_taxonomy,
        "纲": species_taxonomy,
        "目": species_taxonomy,
        "科": species_taxonomy,
        "属": species_taxonomy,
        "种": species_taxonomy,
        "亚种": subspecies_taxonomy,
    }
    palette = [
        "#526a86", "#76834f", "#8a6654", "#6d6481", "#4e7b75",
        "#9b7a3f", "#8a4d47", "#5d7c83", "#7b6d5a", "#697789", "#8d8d8d",
    ]
    ranks: list[dict] = []
    for rank in rank_order:
        dataset = rank_sources.get(rank) or {}
        rows = dataset.get("rows") or []
        if not rows:
            continue
        groups: dict[str, dict[str, float | int]] = {}
        for row in rows:
            name = str(row.get(rank, "")).strip() or "未注释"
            record = groups.setdefault(name, {"ratio": 0.0, "reads": 0})
            record["ratio"] += float(row.get("比例数值") or 0.0)
            record["reads"] += int(row.get("序列数量数值") or 0)
        ranked = sorted(groups.items(), key=lambda item: (item[1]["ratio"], item[1]["reads"]), reverse=True)
        segments = []
        for index, (name, values) in enumerate(ranked):
            segments.append({
                "name": name,
                "ratio": round(float(values["ratio"]), 2),
                "reads": int(values["reads"]),
                "color": palette[index % len(palette)],
            })
        ranks.append({
            "rank": rank,
            "segments": segments,
            "total_ratio": round(sum(segment["ratio"] for segment in segments), 2),
        })
    return {"status": "ready" if ranks else "empty", "ranks": ranks}

def _build_taxonomy_rarefaction(species_taxonomy: dict, subspecies_taxonomy: dict) -> dict:
    def _build_points(dataset: dict, terminal_column: str) -> dict[str, object]:
        rows = dataset.get("rows") or []
        counts: list[int] = []
        for row in rows:
            reads = int(row.get("序列数量数值") or 0)
            label = str(row.get(terminal_column, "")).strip()
            if reads <= 0 or not label or label == "-":
                continue
            counts.append(reads)
        if not counts:
            return {"points": [], "final_expected": None, "plateau_state": "empty", "tail_gain": None}
        total_reads = sum(counts)
        max_points = 32
        if total_reads <= max_points:
            sample_sizes = list(range(1, total_reads + 1))
        else:
            sample_sizes = sorted({max(1, round(total_reads * (index / (max_points - 1)))) for index in range(1, max_points + 1)})
            if sample_sizes[0] != 1:
                sample_sizes.insert(0, 1)
            if sample_sizes[-1] != total_reads:
                sample_sizes.append(total_reads)

        lgamma_total = math.lgamma(total_reads + 1)

        def _prob_not_seen(sample_size: int, taxon_reads: int) -> float:
            if sample_size <= 0:
                return 1.0
            if sample_size > total_reads - taxon_reads:
                return 0.0
            log_prob = (
                math.lgamma(total_reads - taxon_reads + 1)
                - math.lgamma(sample_size + 1)
                - math.lgamma(total_reads - taxon_reads - sample_size + 1)
                - (lgamma_total - math.lgamma(sample_size + 1) - math.lgamma(total_reads - sample_size + 1))
            )
            return math.exp(log_prob)

        points: list[dict[str, object]] = []
        for sample_size in sample_sizes:
            expected_richness = 0.0
            for taxon_reads in counts:
                expected_richness += 1.0 - _prob_not_seen(sample_size, taxon_reads)
            points.append({
                "x": int(sample_size),
                "y": round(expected_richness, 2),
            })
        final_expected = float(points[-1]["y"]) if points else None
        if len(points) >= 2 and final_expected:
            start_index = max(0, int(len(points) * 0.8) - 1)
            tail_gain = float(points[-1]["y"]) - float(points[start_index]["y"])
            tail_ratio = tail_gain / max(final_expected, 1.0)
            if tail_ratio <= 0.03:
                plateau_state = "stable"
            elif tail_ratio <= 0.1:
                plateau_state = "approaching"
            else:
                plateau_state = "rising"
        else:
            tail_gain = None
            plateau_state = "unknown"
        return {
            "points": points,
            "final_expected": round(final_expected, 2) if final_expected is not None else None,
            "plateau_state": plateau_state,
            "tail_gain": round(tail_gain, 2) if tail_gain is not None else None,
        }

    species_curve = _build_points(species_taxonomy or {}, "种")
    subspecies_curve = _build_points(subspecies_taxonomy or {}, "亚种")
    species_points = species_curve["points"]
    subspecies_points = subspecies_curve["points"]
    notes: list[str] = []
    if species_curve["final_expected"] is not None:
        species_note = f"种水平终点期望分类数约 {species_curve['final_expected']}"
        if species_curve["plateau_state"] == "stable":
            species_note += "，曲线已接近平稳"
        elif species_curve["plateau_state"] == "approaching":
            species_note += "，曲线趋于平缓"
        elif species_curve["plateau_state"] == "rising":
            species_note += "，曲线仍有上升空间"
        notes.append(species_note)
    if subspecies_curve["final_expected"] is not None:
        subspecies_note = f"亚种水平终点期望分类数约 {subspecies_curve['final_expected']}"
        if subspecies_curve["plateau_state"] == "stable":
            subspecies_note += "，曲线已接近平稳"
        elif subspecies_curve["plateau_state"] == "approaching":
            subspecies_note += "，曲线趋于平缓"
        elif subspecies_curve["plateau_state"] == "rising":
            subspecies_note += "，曲线仍有上升空间"
        notes.append(subspecies_note)
    return {
        "status": "ready" if species_points or subspecies_points else "empty",
        "label": "分类稀释曲线",
        "x_label": "累计序列数量",
        "y_label": "累计检出分类数",
        "species_points": species_points,
        "subspecies_points": subspecies_points,
        "species_final_expected": species_curve["final_expected"],
        "subspecies_final_expected": subspecies_curve["final_expected"],
        "species_plateau_state": species_curve["plateau_state"],
        "subspecies_plateau_state": subspecies_curve["plateau_state"],
        "note": "；".join(notes) + "。" if notes else "当前分类稀释曲线尚不足以给出稳定性判读。",
    }

def _is_virus_taxonomy_row(row: dict | None) -> bool:
    if not isinstance(row, dict):
        return False
    for value in (
        row.get("界"),
        row.get("种"),
        row.get("亚种"),
        row.get("NCBI学名"),
        row.get("NCBI科"),
        row.get("NCBI属"),
    ):
        text = str(value or "").strip().lower()
        if text and any(token in text for token in ("病毒", "virus", "viridae", "virales")):
            return True
    return False

def _is_fungus_taxonomy_row(row: dict | None) -> bool:
    if not isinstance(row, dict):
        return False
    for value in (
        row.get("界"),
        row.get("种"),
        row.get("亚种"),
        row.get("NCBI学名"),
        row.get("NCBI科"),
        row.get("NCBI属"),
    ):
        text = str(value or "").strip().lower()
        if text and any(token in text for token in ("真菌", "fungi", "fungus", "mycota", "mycetaceae", "mycetales", "mycetes")):
            return True
    return False

def _build_taxonomy_identity_payload(row: dict | None) -> dict | None:
    if not isinstance(row, dict):
        return None
    species_name = str(row.get("种") or row.get("亚种") or row.get("NCBI种") or "").strip()
    if not species_name or species_name == "-":
        return None
    return {
        "species": species_name,
        "ratio": round(float(row.get("比例数值") or 0.0), 2),
        "reads": int(row.get("序列数量数值") or 0),
        "taxid": str(row.get("NCBI TaxID") or "").strip() or "-",
        "scientific_name": str(row.get("NCBI学名") or "").strip() or "-",
        "rank": str(row.get("NCBI分类等级") or "").strip() or "-",
        "order": str(row.get("NCBI目") or "").strip() or "-",
        "family": str(row.get("NCBI科") or "").strip() or "-",
        "genus": str(row.get("NCBI属") or "").strip() or "-",
        "species_rank": str(row.get("NCBI种") or "").strip() or "-",
    }

def _extract_dominant_virus_taxonomy(species_taxonomy: dict, subspecies_taxonomy: dict | None = None) -> dict | None:
    ranked_rows: list[dict] = []
    for dataset in (species_taxonomy or {}, subspecies_taxonomy or {}):
        for row in dataset.get("rows") or []:
            if _is_virus_taxonomy_row(row):
                ranked_rows.append(row)
    ranked_rows.sort(
        key=lambda item: (
            int(item.get("序列数量数值") or 0),
            float(item.get("比例数值") or 0.0),
        ),
        reverse=True,
    )
    return _build_taxonomy_identity_payload(ranked_rows[0]) if ranked_rows else None

def _extract_dominant_fungus_taxonomy(species_taxonomy: dict, subspecies_taxonomy: dict | None = None) -> dict | None:
    ranked_rows: list[dict] = []
    for dataset in (species_taxonomy or {}, subspecies_taxonomy or {}):
        for row in dataset.get("rows") or []:
            if _is_fungus_taxonomy_row(row):
                ranked_rows.append(row)
    ranked_rows.sort(
        key=lambda item: (
            int(item.get("序列数量数值") or 0),
            float(item.get("比例数值") or 0.0),
        ),
        reverse=True,
    )
    return _build_taxonomy_identity_payload(ranked_rows[0]) if ranked_rows else None

def _build_taxonomy_risk_summary(species_taxonomy: dict, subspecies_taxonomy: dict) -> dict:
    datasets = [species_taxonomy or {}, subspecies_taxonomy or {}]
    pathogenicity: dict[str, dict[str, float | int]] = {}
    hazard: dict[str, dict[str, float | int]] = {}
    kingdom_groups: dict[str, dict[str, float | int]] = {
        "细菌": {"reads": 0, "records": 0},
        "病毒": {"reads": 0, "records": 0},
        "真菌": {"reads": 0, "records": 0},
    }
    total_reads = 0
    total_records = 0

    def _normalize_kingdom(value: str) -> str | None:
        text = value.strip().lower()
        if not text or text == "-":
            return None
        if any(token in text for token in ("细菌", "bacteria", "eubacteria")):
            return "细菌"
        if any(token in text for token in ("病毒", "virus", "viruses")):
            return "病毒"
        if any(token in text for token in ("真菌", "fungi", "fungus", "mycota")):
            return "真菌"
        return None

    for dataset in datasets:
        for row in dataset.get("rows") or []:
            reads = int(row.get("序列数量数值") or 0)
            ratio = float(row.get("比例数值") or 0.0)
            total_reads += reads
            total_records += 1

            kingdom_label = _normalize_kingdom(str(row.get("界", "")).strip())
            if kingdom_label:
                kingdom_groups[kingdom_label]["reads"] += reads
                kingdom_groups[kingdom_label]["records"] += 1

            pathogenic_label = str(row.get("致病性", "")).strip()
            if pathogenic_label and pathogenic_label != "-":
                record = pathogenicity.setdefault(pathogenic_label, {"reads": 0, "ratio": 0.0, "records": 0})
                record["reads"] += reads
                record["ratio"] += ratio
                record["records"] += 1

            hazard_label = str(row.get("危害程度等级", "")).strip()
            if hazard_label and hazard_label != "-":
                record = hazard.setdefault(hazard_label, {"reads": 0, "ratio": 0.0, "records": 0})
                record["reads"] += reads
                record["ratio"] += ratio
                record["records"] += 1

    def _format_rank(mapping: dict[str, dict[str, float | int]]) -> list[dict]:
        return [
            {
                "label": label,
                "reads": int(values["reads"]),
                "ratio": round(float(values["ratio"]), 2),
                "records": int(values["records"]),
            }
            for label, values in sorted(mapping.items(), key=lambda item: (item[1]["reads"], item[1]["ratio"]), reverse=True)
        ]

    pathogenicity_rows = _format_rank(pathogenicity)
    hazard_rows = _format_rank(hazard)
    dominant_pathogenicity = pathogenicity_rows[0]["label"] if pathogenicity_rows else ""
    dominant_hazard = hazard_rows[0]["label"] if hazard_rows else ""

    if pathogenicity_rows or hazard_rows:
        summary_parts = []
        if dominant_pathogenicity:
            summary_parts.append(f"当前序列物种鉴定结果以“{dominant_pathogenicity}”类型为主")
        if dominant_hazard:
            summary_parts.append(f"危害程度以“{dominant_hazard}”为主要等级")
        if total_reads:
            summary_parts.append(f"纳入统计的分类结果覆盖 {total_reads} 条序列")
        narrative = "，".join(summary_parts) + "。"
    else:
        narrative = "当前序列物种鉴定结果未提供可汇总的致病性或危害程度等级信息。"

    return {
        "status": "ready" if pathogenicity_rows or hazard_rows else "empty",
        "dominant_virus": _extract_dominant_virus_taxonomy(species_taxonomy, subspecies_taxonomy),
        "dominant_fungus": _extract_dominant_fungus_taxonomy(species_taxonomy, subspecies_taxonomy),
        "kingdom_summary": [
            {
                "label": label,
                "reads": int(values["reads"]),
                "records": int(values["records"]),
                "ratio": round((int(values["reads"]) / total_reads) * 100, 2) if total_reads else 0.0,
            }
            for label, values in kingdom_groups.items()
        ],
        "pathogenicity": pathogenicity_rows,
        "hazard": hazard_rows,
        "narrative": narrative,
        "total_reads": total_reads,
        "total_records": total_records,
    }

def _build_taxonomy_interpretation(species_taxonomy: dict, checkm_info: dict | None = None) -> dict:
    rows = list(species_taxonomy.get("rows") or [])
    valid_rows = []
    for row in rows:
        species_name = str(row.get("种", "")).strip()
        if not species_name or species_name == "-":
            continue
        ratio = float(row.get("比例数值") or 0.0)
        reads = int(row.get("序列数量数值") or 0)
        genus_name = str(row.get("属", "")).strip() or (species_name.split()[0] if species_name else "")
        valid_rows.append({
            "species": species_name,
            "genus": genus_name,
            "ratio": ratio,
            "reads": reads,
        })
    valid_rows.sort(key=lambda item: (item["reads"], item["ratio"]), reverse=True)

    if not valid_rows:
        empty_card = {
            "status": "empty",
            "tone": "neutral",
            "badge": "暂无判读",
            "headline": "未检出足够的物种分类结果",
            "summary": "当前分类结果不足，暂无法对相近物种混淆或单菌/混菌状态给出可靠提示。",
            "metrics": [],
            "evidence": [],
        }
        return {"status": "empty", "confusion_hint": dict(empty_card), "mixture_hint": dict(empty_card)}

    top = valid_rows[0]
    second = valid_rows[1] if len(valid_rows) > 1 else None
    same_genus_rows = [row for row in valid_rows[1:] if row["genus"] and row["genus"] == top["genus"]]
    same_genus_ratio = round(sum(row["ratio"] for row in same_genus_rows), 2)
    same_genus_count = len(same_genus_rows)
    species_over_five = [row for row in valid_rows if row["ratio"] >= 5]
    top_two_ratio = round(top["ratio"] + (second["ratio"] if second else 0.0), 2)
    competing_ratio = round(second["ratio"], 2) if second else 0.0
    competing_reads = int(second["reads"]) if second else 0
    competing_species = second["species"] if second else "--"
    checkm_species = str((checkm_info or {}).get("species_name") or "").strip()
    mlst_species = str((checkm_info or {}).get("mlst_species_name") or "").strip()

    confusion_tone = "success"
    confusion_badge = "较稳定"
    confusion_headline = "未见明显相近物种混淆"
    confusion_summary = f"当前分类结果以 {top['species']} 为主导，相近物种干扰信号不强。"
    confusion_evidence = [
        f"主导物种为 {top['species']}，序列占比约 {top['ratio']:.2f}%。",
    ]
    if second:
        confusion_evidence.append(f"第二位物种为 {second['species']}，占比约 {second['ratio']:.2f}%。")
    if same_genus_count and second and second["genus"] == top["genus"]:
        if second["ratio"] >= 10 or (top["reads"] and second["reads"] / max(top["reads"], 1) >= 0.35) or same_genus_count >= 3:
            confusion_tone = "warning"
            confusion_badge = "需要关注"
            confusion_headline = "存在相近物种混淆风险"
            confusion_summary = (
                f"{top['genus']} 属内同时出现多个高占比物种，"
                f"其中 {competing_species} 与主导物种接近，需结合其他证据进一步确认。"
            )
        elif second["ratio"] >= 5 or same_genus_ratio >= 12:
            confusion_tone = "attention"
            confusion_badge = "提示"
            confusion_headline = "需关注近缘物种干扰"
            confusion_summary = (
                f"{top['genus']} 属内还存在一定比例的近缘物种信号，"
                f"当前更适合保守解读为主导物种附近的近缘分类结果。"
            )
        confusion_evidence.append(
            f"同属候选物种共 {same_genus_count} 个，累计占比约 {same_genus_ratio:.2f}%。"
        )
    elif checkm_species and mlst_species and checkm_species != mlst_species:
        confusion_tone = "attention"
        confusion_badge = "提示"
        confusion_headline = "主流程物种结果存在差异"
        confusion_summary = (
            f"分类主导物种为 {top['species']}，但 CheckM 与 MLST 返回的物种名称不完全一致，建议结合装配和 MLST 结果复核。"
        )
        confusion_evidence.append(f"CheckM 物种：{checkm_species or '--'}；MLST 物种：{mlst_species or '--'}。")

    confusion_metrics = [
        {"label": "主导物种", "value": top["species"]},
        {"label": "次高物种", "value": competing_species},
        {"label": "同属累计占比", "value": f"{same_genus_ratio:.2f}%"},
    ]

    mixture_tone = "success"
    mixture_badge = "倾向单菌"
    mixture_headline = "单菌信号较明确"
    mixture_summary = f"{top['species']} 占比明显，当前结果整体更接近单菌样本。"
    mixture_evidence = [
        f"主导物种 {top['species']} 占比约 {top['ratio']:.2f}%。",
        f"占比 ≥5% 的物种数为 {len(species_over_five)} 个。",
    ]

    if top["ratio"] < 45 or len(species_over_five) >= 3 or (second and second["ratio"] >= 20 and top_two_ratio <= 85):
        mixture_tone = "warning"
        mixture_badge = "疑似混菌"
        mixture_headline = "存在混菌信号"
        mixture_summary = "多个物种占比同时较高，当前结果更像混菌样本或显著背景混入。"
    elif top["ratio"] < 65 or (second and second["ratio"] >= 12) or len(species_over_five) == 2:
        mixture_tone = "attention"
        mixture_badge = "需关注"
        mixture_headline = "单菌/混菌边界不够清晰"
        mixture_summary = "样本存在次高物种背景，当前更适合解释为主导物种明显，但仍需关注可能的混杂信号。"

    if second:
        mixture_evidence.append(
            f"第二位物种 {second['species']} 占比约 {second['ratio']:.2f}%（{competing_reads} 条序列）。"
        )
    mixture_evidence.append(f"前两位物种合计占比约 {top_two_ratio:.2f}%。")

    mixture_metrics = [
        {"label": "主导物种占比", "value": f"{top['ratio']:.2f}%"},
        {"label": "次高物种占比", "value": f"{competing_ratio:.2f}%"},
        {"label": "≥5% 物种数", "value": str(len(species_over_five))},
    ]

    return {
        "status": "ready",
        "confusion_hint": {
            "status": "ready",
            "tone": confusion_tone,
            "badge": confusion_badge,
            "headline": confusion_headline,
            "summary": confusion_summary,
            "metrics": confusion_metrics,
            "evidence": confusion_evidence,
        },
        "mixture_hint": {
            "status": "ready",
            "tone": mixture_tone,
            "badge": mixture_badge,
            "headline": mixture_headline,
            "summary": mixture_summary,
            "metrics": mixture_metrics,
            "evidence": mixture_evidence,
        },
    }

def _extract_dominant_species(species_taxonomy: dict) -> dict:
    ranked = []
    for row in species_taxonomy.get("rows") or []:
        species_name = str(row.get("种", "")).strip()
        if not species_name or species_name == "-":
            continue
        ranked.append({
            "species": species_name,
            "ratio": float(row.get("比例数值") or 0.0),
            "reads": int(row.get("序列数量数值") or 0),
        })
    ranked.sort(key=lambda item: (item["reads"], item["ratio"]), reverse=True)
    return ranked[0] if ranked else {"species": "--", "ratio": 0.0, "reads": 0}

def _read_assembly_taxonomy(assem_info_path: Path, kraken_report_path: Path, taxonomy_index: dict[str, list[dict]] | None = None) -> dict:
    assem_info = _read_tsv_rows(assem_info_path)
    if not assem_info["columns"] or not assem_info["rows"]:
        return {"columns": [], "rows": []}

    rank_mapping = {
        "D": "界",
        "P": "门",
        "C": "纲",
        "O": "目",
        "F": "科",
        "G": "属",
        "S": "种",
    }
    tracked_ranks = ["界", "门", "纲", "目", "科", "属", "种"]
    lineage_by_taxid: dict[str, dict[str, str]] = {}
    lineage_stack: dict[str, str] = {}

    if kraken_report_path.is_file():
        try:
            with kraken_report_path.open("r", encoding="utf-8", errors="ignore") as handle:
                for line in handle:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 6:
                        continue
                    rank_code = str(parts[3]).strip()
                    taxid = str(parts[4]).strip()
                    raw_name = str(parts[5]).rstrip()
                    clean_name = raw_name.strip()
                    if not taxid or not clean_name:
                        continue

                    normalized_rank = rank_code[0] if rank_code else ""
                    if normalized_rank not in rank_mapping:
                        continue

                    current_rank = rank_mapping[normalized_rank]
                    current_index = tracked_ranks.index(current_rank)
                    for stale_rank in tracked_ranks[current_index:]:
                        lineage_stack.pop(stale_rank, None)
                    lineage_stack[current_rank] = clean_name

                    lineage_by_taxid[taxid] = {rank: lineage_stack.get(rank, "") for rank in tracked_ranks}
        except OSError:
            lineage_by_taxid = {}

    assem_records = [
        {assem_info["columns"][index]: row[index] if index < len(row) else "" for index in range(len(assem_info["columns"]))}
        for row in assem_info["rows"]
    ]

    columns = [
        "序列名称",
        "序列长度",
        "平均深度",
        "是否成环",
        "基因组/质粒",
        "质粒分型",
        "taxid",
        "物种名称",
        "NCBI TaxID",
        "NCBI学名",
        "NCBI分类等级",
        "界",
        "门",
        "纲",
        "目",
        "科",
        "属",
        "NCBI目",
        "NCBI科",
        "NCBI属",
        "NCBI种",
    ]
    rows: list[list[str]] = []
    for record in assem_records:
        taxid = str(record.get("taxid", "")).strip()
        lineage = lineage_by_taxid.get(taxid, {})
        species_name = str(record.get("物种名称", "")).strip()
        matched_taxonomy = _lookup_kb_taxonomy(
            taxonomy_index,
            species_name=species_name,
            genus_name=str(lineage.get("属", "")).strip(),
        ) or {}
        rows.append([
            str(record.get("序列名称", "")).strip(),
            str(record.get("序列长度", "")).strip(),
            str(record.get("平均深度", "")).strip(),
            str(record.get("是否成环", "")).strip(),
            str(record.get("基因组/质粒", "")).strip(),
            str(record.get("质粒分型", "")).strip(),
            taxid,
            species_name,
            str(matched_taxonomy.get("taxid") or "-"),
            str(matched_taxonomy.get("scientific_name") or "-"),
            str(matched_taxonomy.get("rank") or "-"),
            str(lineage.get("界", "")).strip(),
            str(lineage.get("门", "")).strip(),
            str(lineage.get("纲", "")).strip(),
            str(lineage.get("目", "")).strip(),
            str(lineage.get("科", "")).strip(),
            str(lineage.get("属", "")).strip(),
            str(matched_taxonomy.get("order") or "-"),
            str(matched_taxonomy.get("family") or "-"),
            str(matched_taxonomy.get("genus") or "-"),
            str(matched_taxonomy.get("species_rank") or "-"),
        ])
    return {"columns": columns, "rows": rows}
