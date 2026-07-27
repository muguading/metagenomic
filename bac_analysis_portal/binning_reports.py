from __future__ import annotations

from pathlib import Path

from .parse_utils import _safe_float, _safe_int
from .report_formatters import _human_bp
from .report_sources import _read_fasta_assembly_summary
from .table_io import _read_tsv_rows

def _read_checkm2_quality(path: Path) -> dict:
    raw = _read_tsv_rows(path)
    if not raw["columns"]:
        return {"columns": [], "rows": []}
    rename_map = {
        "Name": "样本名称",
        "Completeness": "完整性",
        "Contamination": "污染率",
        "Completeness_Model_Used": "完整性模型",
        "Translation_Table_Used": "翻译表",
        "Coding_Density": "编码密度",
        "Contig_N50": "Contig N50",
        "Average_Gene_Length": "平均基因长度",
        "Genome_Size": "基因组大小",
        "GC_Content": "GC 含量",
        "Total_Coding_Sequences": "编码序列数",
        "Total_Contigs": "Contig 总数",
        "Max_Contig_Length": "最大 Contig 长度",
        "Additional_Notes": "附加说明",
    }
    keep_columns = [
        "Name", "Completeness", "Contamination", "Genome_Size",
        "GC_Content", "Total_Contigs", "Contig_N50", "Max_Contig_Length",
        "Coding_Density", "Total_Coding_Sequences", "Completeness_Model_Used", "Additional_Notes",
    ]
    indices = [raw["columns"].index(col) for col in keep_columns if col in raw["columns"]]
    columns = [rename_map.get(raw["columns"][index], raw["columns"][index]) for index in indices]
    rows = []
    for row in raw["rows"]:
        trimmed = []
        for index in indices:
            value = row[index] if index < len(row) else ""
            if raw["columns"][index] in {"Completeness", "Contamination"} and value not in {"", "None", None}:
                try:
                    value = f"{float(value):.2f}%"
                except ValueError:
                    pass
            elif raw["columns"][index] == "GC_Content" and value not in {"", "None", None}:
                try:
                    numeric = float(value)
                    value = f"{numeric * 100:.2f}%" if numeric <= 1 else f"{numeric:.2f}%"
                except ValueError:
                    pass
            elif raw["columns"][index] == "Genome_Size" and value not in {"", "None", None}:
                value = _human_bp(value)
            elif raw["columns"][index] == "Coding_Density" and value not in {"", "None", None}:
                try:
                    numeric = float(value)
                    value = f"{numeric * 100:.2f}%" if numeric <= 1 else f"{numeric:.2f}%"
                except ValueError:
                    pass
            trimmed.append(value)
        rows.append(trimmed)
    return {"columns": columns, "rows": rows}

def _bin_quality_tier(completeness: float | None, contamination: float | None) -> str:
    comp = float(completeness or 0.0)
    cont = float(contamination or 0.0)
    if comp >= 90 and cont <= 5:
        return "高质量"
    if comp >= 50 and cont <= 10:
        return "中质量"
    return "低质量"

def _build_binning_quality_section(path: Path) -> dict:
    raw = _read_tsv_rows(path)
    if not raw["columns"] or not raw["rows"]:
        return {"status": "empty", "summary": {}, "charts": {}, "table": {"columns": [], "rows": []}}
    name_map = _read_binning_name_map(path.parent.parent)
    records: list[dict[str, object]] = []
    for row in raw["rows"]:
        record = {raw["columns"][index]: row[index] if index < len(row) else "" for index in range(len(raw["columns"]))}
        completeness = _safe_float(record.get("Completeness"))
        contamination = _safe_float(record.get("Contamination"))
        genome_size = _safe_int(record.get("Genome_Size"))
        total_contigs = _safe_int(record.get("Total_Contigs"))
        contig_n50 = _safe_int(record.get("Contig_N50"))
        gc_content = _safe_float(record.get("GC_Content"))
        tier = _bin_quality_tier(completeness, contamination)
        records.append(
            {
                "bin_id": name_map.get(str(record.get("Name") or "").strip(), str(record.get("Name") or "").strip()),
                "completeness": completeness,
                "contamination": contamination,
                "tier": tier,
                "genome_size": genome_size,
                "total_contigs": total_contigs,
                "contig_n50": contig_n50,
                "gc_content": gc_content,
                "coding_density": _safe_float(record.get("Coding_Density")),
                "coding_sequences": _safe_int(record.get("Total_Coding_Sequences")),
                "model": str(record.get("Completeness_Model_Used") or "").strip(),
                "notes": str(record.get("Additional_Notes") or "").strip(),
            }
        )
    records = [item for item in records if item.get("bin_id")]
    if not records:
        return {"status": "empty", "summary": {}, "charts": {}, "table": {"columns": [], "rows": []}}

    completeness_buckets = [
        ("<50", lambda v: v is not None and v < 50),
        ("50-70", lambda v: v is not None and 50 <= v < 70),
        ("70-90", lambda v: v is not None and 70 <= v < 90),
        ("≥90", lambda v: v is not None and v >= 90),
    ]
    contamination_buckets = [
        ("0-5", lambda v: v is not None and 0 <= v <= 5),
        ("5-10", lambda v: v is not None and 5 < v <= 10),
        ("10-20", lambda v: v is not None and 10 < v <= 20),
        (">20", lambda v: v is not None and v > 20),
    ]
    quality_counts = {"高质量": 0, "中质量": 0, "低质量": 0}
    for item in records:
        quality_counts[str(item["tier"])] += 1

    summary = {
        "total_bins": len(records),
        "hq_bins": quality_counts["高质量"],
        "mq_bins": quality_counts["中质量"],
        "lq_bins": quality_counts["低质量"],
        "avg_completeness": round(sum(float(item["completeness"] or 0.0) for item in records) / len(records), 2),
        "avg_contamination": round(sum(float(item["contamination"] or 0.0) for item in records) / len(records), 2),
    }
    charts = {
        "completeness": {
            "label": "bin 完整性分布",
            "x_label": "完整性区间",
            "y_label": "bin 数量",
            "x_values": [label for label, _ in completeness_buckets],
            "points": [sum(1 for item in records if predicate(item["completeness"])) for label, predicate in completeness_buckets],
        },
        "contamination": {
            "label": "bin 污染率分布",
            "x_label": "污染率区间",
            "y_label": "bin 数量",
            "x_values": [label for label, _ in contamination_buckets],
            "points": [sum(1 for item in records if predicate(item["contamination"])) for label, predicate in contamination_buckets],
        },
        "quality_tier": {
            "label": "bin 质量分层",
            "x_label": "质量等级",
            "y_label": "bin 数量",
            "x_values": ["高质量", "中质量", "低质量"],
            "points": [quality_counts["高质量"], quality_counts["中质量"], quality_counts["低质量"]],
        },
    }
    table_columns = ["Bin名称", "质量等级", "完整性", "污染率", "基因组大小", "Contig总数", "Contig N50", "GC含量", "编码密度", "编码序列数", "完整性模型", "附加说明"]
    table_rows: list[list[str]] = []
    for item in records:
        gc_value = item["gc_content"]
        gc_display = "-" if gc_value is None else (f"{float(gc_value) * 100:.2f}%" if float(gc_value) <= 1 else f"{float(gc_value):.2f}%")
        density_value = item["coding_density"]
        density_display = "-" if density_value is None else (f"{float(density_value) * 100:.2f}%" if float(density_value) <= 1 else f"{float(density_value):.2f}%")
        table_rows.append(
            [
                str(item["bin_id"]),
                str(item["tier"]),
                "-" if item["completeness"] is None else f"{float(item['completeness']):.2f}%",
                "-" if item["contamination"] is None else f"{float(item['contamination']):.2f}%",
                _human_bp(item["genome_size"]) if item["genome_size"] else "-",
                str(item["total_contigs"] or "-"),
                str(item["contig_n50"] or "-"),
                gc_display,
                density_display,
                str(item["coding_sequences"] or "-"),
                str(item["model"] or "-"),
                str(item["notes"] or "-"),
            ]
        )
    return {
        "status": "ready",
        "summary": summary,
        "charts": charts,
        "table": {"columns": table_columns, "rows": table_rows},
    }

def _parse_gtdb_classification(value: str) -> dict[str, str]:
    result = {"domain": "-", "phylum": "-", "class": "-", "order": "-", "family": "-", "genus": "-", "species": "-"}
    for part in str(value or "").split(";"):
        text = part.strip()
        if not text or "__" not in text:
            continue
        prefix, label = text.split("__", 1)
        label = label.strip() or "-"
        mapping = {
            "d": "domain",
            "p": "phylum",
            "c": "class",
            "o": "order",
            "f": "family",
            "g": "genus",
            "s": "species",
        }
        target = mapping.get(prefix.strip().lower())
        if target:
            result[target] = label
    return result

def _read_binning_name_map(base_dir: Path) -> dict[str, str]:
    raw = _read_tsv_rows(base_dir / "binning_name.tsv")
    columns = raw.get("columns") or []
    rows = raw.get("rows") or []
    if not columns or not rows or "oldname" not in columns or "newname" not in columns:
        return {}
    old_index = columns.index("oldname")
    new_index = columns.index("newname")
    mapping: dict[str, str] = {}
    for row in rows:
        old_name = str(row[old_index] if old_index < len(row) else "").strip()
        new_name = str(row[new_index] if new_index < len(row) else "").strip()
        if old_name and new_name:
            mapping[old_name] = new_name
    return mapping

def _build_binning_taxonomy_section(path: Path) -> dict:
    raw = _read_tsv_rows(path)
    if not raw["columns"] or not raw["rows"]:
        return {"status": "empty", "summary": {}, "charts": {}, "table": {"columns": [], "rows": []}}
    name_map = _read_binning_name_map(path.parent.parent)
    records: list[dict[str, str]] = []
    for row in raw["rows"]:
        record = {raw["columns"][index]: row[index] if index < len(row) else "" for index in range(len(raw["columns"]))}
        lineage = _parse_gtdb_classification(str(record.get("classification") or ""))
        original_bin_id = str(record.get("user_genome") or "").strip()
        records.append(
            {
                "bin_id": name_map.get(original_bin_id, original_bin_id),
                "domain": lineage["domain"],
                "phylum": lineage["phylum"],
                "genus": lineage["genus"],
                "species": lineage["species"],
                "classification": str(record.get("classification") or "").strip(),
                "reference": str(record.get("closest_genome_reference") or "").strip() or "-",
                "ani": str(record.get("closest_genome_ani") or "").strip() or "-",
                "af": str(record.get("closest_genome_af") or "").strip() or "-",
                "method": str(record.get("classification_method") or "").strip() or "-",
                "warning": str(record.get("warnings") or "").strip() or "-",
            }
        )
    records = [item for item in records if item.get("bin_id")]
    if not records:
        return {"status": "empty", "summary": {}, "charts": {}, "table": {"columns": [], "rows": []}}

    def top_counts(key: str, limit: int = 8) -> tuple[list[str], list[int]]:
        counts: dict[str, int] = {}
        for item in records:
            label = str(item.get(key) or "-").strip() or "-"
            counts[label] = counts.get(label, 0) + 1
        ranked = sorted(counts.items(), key=lambda x: (-x[1], x[0]))[:limit]
        return [name for name, _ in ranked], [count for _, count in ranked]

    phylum_labels, phylum_counts = top_counts("phylum", 8)
    genus_labels, genus_counts = top_counts("genus", 8)
    method_labels, method_counts = top_counts("method", 6)
    classified_count = sum(1 for item in records if str(item.get("species") or "-") not in {"-", "Unclassified"})
    summary = {
        "total_bins": len(records),
        "classified_bins": classified_count,
        "unclassified_bins": len(records) - classified_count,
        "top_phylum": phylum_labels[0] if phylum_labels else "-",
        "top_genus": genus_labels[0] if genus_labels else "-",
    }
    charts = {
        "phylum": {
            "label": "bin 门水平分布",
            "x_label": "门",
            "y_label": "bin 数量",
            "x_values": phylum_labels,
            "points": phylum_counts,
        },
        "genus": {
            "label": "bin 属水平分布",
            "x_label": "属",
            "y_label": "bin 数量",
            "x_values": genus_labels,
            "points": genus_counts,
        },
        "method": {
            "label": "GTDB-Tk 分类方法统计",
            "x_label": "分类方法",
            "y_label": "bin 数量",
            "x_values": method_labels,
            "points": method_counts,
        },
    }
    table_columns = ["Bin名称", "门", "属", "种", "参考基因组", "ANI", "AF", "分类方法", "完整分类结果", "警告信息"]
    table_rows = [
        [
            item["bin_id"],
            item["phylum"],
            item["genus"],
            item["species"],
            item["reference"],
            item["ani"],
            item["af"],
            item["method"],
            item["classification"] or "-",
            item["warning"],
        ]
        for item in records
    ]
    return {
        "status": "ready",
        "summary": summary,
        "charts": charts,
        "table": {"columns": table_columns, "rows": table_rows},
    }

def _build_meta_viral_assembly_section(base_dir: Path) -> dict:
    summary_raw = _read_tsv_rows(base_dir / "viral_summary.tsv")
    contig_table = _read_tsv_rows(base_dir / "viral_contig_summary.tsv")
    retained_summary = _read_fasta_assembly_summary(base_dir / "viral_retained_contigs.fa")
    raw_contigs_summary = _read_fasta_assembly_summary(base_dir / "megahit_output" / "final.contigs.fa")
    if not contig_table["columns"] and not summary_raw["rows"] and raw_contigs_summary.get("status") == "empty":
        return {"status": "empty", "summary": {}, "table": {"columns": [], "rows": []}}

    summary_map: dict[str, str] = {}
    for row in summary_raw.get("rows", []):
        if len(row) >= 2:
            summary_map[str(row[0]).strip()] = str(row[1]).strip()

    total_contigs = _safe_int(summary_map.get("病毒候选contig数"))
    retained_contigs = _safe_int(summary_map.get("最终保留contig数"))
    retained_length = _safe_int(summary_map.get("总保留长度"))
    note = str(summary_map.get("说明") or "").strip()
    rows = contig_table.get("rows") or []
    kept_rows = [row for row in rows if len(row) >= 11 and str(row[10]).strip().lower() == "yes"]
    quality_index = contig_table["columns"].index("checkv_quality") if "checkv_quality" in contig_table["columns"] else -1
    best_quality = "-"
    if quality_index >= 0 and kept_rows:
        qualities = [str(row[quality_index]).strip() for row in kept_rows if quality_index < len(row) and str(row[quality_index]).strip()]
        for label in ["Complete", "High-quality", "Medium-quality", "Low-quality", "Not-determined"]:
            if label in qualities:
                best_quality = label
                break
    summary = {
        "candidate_contigs": total_contigs if total_contigs is not None else len(rows),
        "retained_contigs": retained_contigs if retained_contigs is not None else len(kept_rows),
        "retained_length": retained_length if retained_length is not None else retained_summary.get("total_length"),
        "raw_contigs": raw_contigs_summary.get("contig_count"),
        "best_quality": best_quality,
        "note": note or "-",
    }
    display_columns = {
        "contig_id": "Contig",
        "contig_length": "长度",
        "virsorter2_score": "VirSorter2得分",
        "virsorter2_group": "病毒组别",
        "hallmark": "Hallmark基因",
        "viral_genes": "病毒基因数",
        "host_genes": "宿主基因数",
        "checkv_quality": "CheckV质量",
        "completeness": "完整性",
        "contamination": "污染率",
        "retained": "是否保留",
        "retention_reason": "保留依据",
    }
    if contig_table["columns"]:
        remapped_columns = [display_columns.get(column, column) for column in contig_table["columns"]]
    else:
        remapped_columns = []
    return {
        "status": "ready",
        "summary": summary,
        "table": {"columns": remapped_columns, "rows": rows},
    }
