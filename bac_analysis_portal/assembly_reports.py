from __future__ import annotations

import csv
import gzip
import json
import os
import shlex
import subprocess
from pathlib import Path
from typing import Any

from .parse_utils import _safe_float, _safe_int
from .report_artifacts import _resolve_report_artifact_path
from .report_sources import _read_fasta_assembly_summary
from .table_io import _read_tsv_rows

def _discover_cgview_assets(report_dir: Path, sample_name: str) -> dict:
    if not report_dir.is_dir() or not sample_name:
        return {"status": "empty", "summary": {"map_count": 0}, "maps": []}
    prokka_dir = report_dir / f"{sample_name}_prokka"
    if not prokka_dir.is_dir():
        return {"status": "empty", "summary": {"map_count": 0}, "maps": []}
    maps: list[dict[str, str]] = []
    main_gbk = prokka_dir / "main.gbk"
    fallback_main_gbk = prokka_dir / f"{sample_name}.gbk"
    if main_gbk.is_file():
        maps.append({
            "key": "main",
            "label": "主基因组环形图",
            "asset_name": f"{prokka_dir.name}/main.gbk",
            "role": "main",
        })
    elif fallback_main_gbk.is_file():
        maps.append({
            "key": sample_name,
            "label": "主基因组环形图",
            "asset_name": f"{prokka_dir.name}/{sample_name}.gbk",
            "role": "main",
        })
    for path in sorted(prokka_dir.glob("*.gbk"), key=lambda item: item.name.lower()):
        if path.name.startswith("."):
            continue
        if path.name in {"main.gbk", f"{sample_name}.gbk"}:
            continue
        maps.append({
            "key": path.stem,
            "label": f"{path.stem} 环形图",
            "asset_name": f"{prokka_dir.name}/{path.name}",
            "role": "plasmid" if "plasmid" in path.stem.lower() else "contig",
        })

    return {
        "status": "ready" if maps else "empty",
        "summary": {"map_count": len(maps)},
        "maps": maps,
    }

def _build_cgview_json_from_gbk(report_dir: Path, gbk_path: Path, map_name: str) -> Path:
    cache_dir = report_dir / ".portal_report_cache" / "cgview"
    cache_dir.mkdir(parents=True, exist_ok=True)
    json_path = cache_dir / f"{map_name}.cgview.json"
    if json_path.is_file():
        try:
            if json_path.stat().st_mtime >= gbk_path.stat().st_mtime:
                return json_path
        except OSError:
            pass

    builder_command = str(
        os.environ.get("CGVIEW_BUILDER_BIN")
        or "ruby /data/deploy/bio-elite/bio/script/cgview_builder_cli.rb"
    ).strip()
    config_path = str(
        os.environ.get("CGVIEW_CONFIG_PATH")
        or "/data/deploy/bio-elite/bio/script/CGview.yaml"
    ).strip()
    command = shlex.split(builder_command) + [
        "-s", str(gbk_path),
        "-o", str(json_path),
        "-c", config_path,
        "-n", map_name,
    ]
    completed = subprocess.run(
        command,
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode == 0 and json_path.is_file():
        return json_path

    builder_error = completed.stderr.strip() or completed.stdout.strip() or "CGView builder 执行失败"
    try:
        return _build_cgview_json_with_biopython(report_dir, gbk_path, map_name)
    except Exception as fallback_error:
        raise RuntimeError(f"{builder_error}; Python fallback 也失败: {fallback_error}") from fallback_error
    return json_path

def _build_cgview_json_with_biopython(report_dir: Path, gbk_path: Path, map_name: str) -> Path:
    try:
        from Bio import BiopythonParserWarning, SeqIO
    except ImportError as error:
        raise RuntimeError("Biopython 不可用，无法执行 CGView JSON 回退构建") from error

    import warnings

    cache_dir = report_dir / ".portal_report_cache" / "cgview"
    cache_dir.mkdir(parents=True, exist_ok=True)
    json_path = cache_dir / f"{map_name}.cgview.json"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", BiopythonParserWarning)
        record = SeqIO.read(str(gbk_path), "genbank")

    contig_name = str(getattr(record, "id", "") or map_name or gbk_path.stem).strip() or gbk_path.stem
    sequence = str(getattr(record, "seq", "") or "").lower()
    contig_length = len(sequence)
    if contig_length <= 0:
        raise RuntimeError(f"GenBank 文件 {gbk_path.name} 不包含可绘制的序列长度")

    features: list[dict[str, Any]] = []
    for feature in getattr(record, "features", []) or []:
        if not feature or getattr(feature, "type", "") == "source":
            continue
        try:
            start = int(feature.location.start) + 1
            stop = int(feature.location.end)
        except Exception:
            continue
        if start <= 0 or stop <= 0:
            continue
        start, stop = sorted((start, stop))
        qualifiers = getattr(feature, "qualifiers", {}) or {}
        name = (
            next(iter(qualifiers.get("gene", [])), "")
            or next(iter(qualifiers.get("locus_tag", [])), "")
            or next(iter(qualifiers.get("product", [])), "")
            or str(getattr(feature, "type", "") or "feature")
        )
        sanitized_qualifiers = {
            str(key): [str(item) for item in value]
            for key, value in qualifiers.items()
            if isinstance(value, list)
        }
        features.append({
            "start": start,
            "stop": stop,
            "strand": -1 if getattr(feature.location, "strand", 1) == -1 else 1,
            "name": str(name).strip() or str(getattr(feature, "type", "") or "feature"),
            "type": str(getattr(feature, "type", "") or "misc_feature"),
            "contig": contig_name,
            "source": "genbank-features",
            "legend": str(getattr(feature, "type", "") or "misc_feature"),
            "qualifiers": sanitized_qualifiers,
        })

    if not features:
        features.append({
            "start": 1,
            "stop": contig_length,
            "strand": 1,
            "name": contig_name,
            "type": "misc_feature",
            "contig": contig_name,
            "source": "genbank-features",
            "legend": "misc_feature",
            "qualifiers": {"note": ["fallback genome span"]},
        })

    payload = {
        "cgview": {
            "name": map_name,
            "settings": {
                "format": "circular",
                "backgroundColor": "rgb(255,255,255)",
                "geneticCode": 11,
            },
            "backbone": {},
            "ruler": {},
            "dividers": {},
            "annotation": {},
            "sequence": {
                "contigs": [{
                    "name": contig_name,
                    "length": contig_length,
                    "seq": sequence,
                }],
            },
            "legend": {"visible": True},
            "tracks": [{
                "name": "Features",
                "separateFeaturesBy": "strand",
                "position": "both",
                "dataType": "feature",
                "dataMethod": "source",
                "dataKeys": "genbank-features",
            }],
            "captions": [],
            "version": "1.7.0",
            "features": features,
        },
    }
    json_path.write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")
    return json_path

def _fallback_assembly_profile(assembly_profile: dict, contig_annotation: dict, assembly_summary: dict) -> dict:
    if assembly_profile.get("contig_count") is not None:
        return assembly_profile
    columns = contig_annotation.get("columns") or []
    rows = contig_annotation.get("rows") or []
    if columns and rows:
        length_index = columns.index("序列长度") if "序列长度" in columns else -1
        genome_type_index = columns.index("基因组/质粒") if "基因组/质粒" in columns else -1
        contig_count = len(rows)
        plasmid_count = 0
        total_length = 0
        has_total_length = False
        for row in rows:
            if genome_type_index >= 0 and genome_type_index < len(row):
                raw_type = str(row[genome_type_index] or "")
                if "plasmid" in raw_type.lower():
                    plasmid_count += 1
            if length_index >= 0 and length_index < len(row):
                length_value = _safe_int(row[length_index])
                if length_value is not None:
                    total_length += length_value
                    has_total_length = True
        return {
            "contig_count": contig_count,
            "plasmid_count": plasmid_count,
            "total_count": contig_count,
            "total_length": total_length if has_total_length else None,
        }
    columns = assembly_summary.get("columns") or []
    rows = assembly_summary.get("rows") or []
    if not columns or not rows:
        return assembly_profile
    try:
        contig_index = columns.index("Contig数量")
    except ValueError:
        return assembly_profile
    first_row = rows[0] if rows else []
    contig_count = _safe_int(first_row[contig_index] if contig_index < len(first_row) else None)
    if contig_count is None:
        return assembly_profile
    return {
        "contig_count": contig_count,
        "plasmid_count": 0,
        "total_count": contig_count,
        "total_length": assembly_profile.get("total_length"),
    }

def _read_contig_depth_relationship(path: Path) -> dict:
    if not path.is_file():
        return {"status": "empty", "points": []}
    points: list[dict[str, object]] = []
    scatter_points: list[dict[str, object]] = []
    try:
        with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="	")
            for row in reader:
                contig_name = str(row.get("序列名称", "")).strip()
                depth = _safe_float(row.get("平均深度"))
                length = _safe_int(row.get("序列长度"))
                raw_type = str(row.get("基因组/质粒", "")).strip()
                if not contig_name or depth is None:
                    continue
                seq_type = "质粒" if "plasmid" in raw_type.lower() else "基因组"
                points.append({
                    "name": contig_name,
                    "depth": round(depth, 2),
                    "type": seq_type,
                })
                if length is not None and length > 0:
                    scatter_points.append({
                        "name": contig_name,
                        "depth": round(depth, 2),
                        "length": length,
                        "type": seq_type,
                    })
    except OSError:
        return {"status": "empty", "points": []}
    return {
        "status": "ready" if points else "empty",
        "label": "基因组/质粒与平均深度关系图",
        "x_label": "序列类型",
        "y_label": "平均深度",
        "points": points,
        "length_depth_scatter": {
            "status": "ready" if scatter_points else "empty",
            "label": "Contig长度与平均测序深度散点图",
            "x_label": "Contig长度(bp)",
            "y_label": "平均测序深度",
            "points": scatter_points,
        },
    }

def _read_coverage_profile(path: Path, target_points: int = 800) -> dict:
    if not path.is_file():
        return {"status": "empty", "points": []}
    depths: list[int] = []
    contigs: set[str] = set()
    contig_depths: dict[str, list[int]] = {}
    opener = gzip.open if path.suffix == ".gz" else open
    try:
        with opener(path, "rt", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 4:
                    continue
                contig_name = str(parts[0]).strip() or "unknown"
                contigs.add(contig_name)
                depth = _safe_float(parts[3])
                if depth is None:
                    continue
                depth_value = int(round(depth))
                depths.append(depth_value)
                contig_depths.setdefault(contig_name, []).append(depth_value)
    except OSError:
        return {"status": "empty", "points": []}
    if not depths:
        return {"status": "empty", "points": []}

    def _compress_depths(raw_depths: list[int], max_points: int) -> tuple[list[float], list[int]]:
        total = len(raw_depths)
        if total <= max_points:
            return [round(value, 2) for value in raw_depths], list(range(1, total + 1))
        chunk_size = max(1, (total + max_points - 1) // max_points)
        compressed_points: list[float] = []
        compressed_x: list[int] = []
        for start in range(0, total, chunk_size):
            chunk = raw_depths[start:start + chunk_size]
            if not chunk:
                continue
            compressed_points.append(round(sum(chunk) / len(chunk), 2))
            compressed_x.append(start + 1)
        return compressed_points, compressed_x

    total_bases = len(depths)
    points, x_values = _compress_depths(depths, target_points)
    max_depth = max(points) if points else None
    mean_depth = round(sum(depths) / total_bases, 2) if depths else None
    covered_bases = sum(1 for depth in depths if depth > 0)
    covered_10x_bases = sum(1 for depth in depths if depth >= 10)
    covered_100x_bases = sum(1 for depth in depths if depth >= 100)
    x_ticks = [1, max(1, total_bases // 2), total_bases]
    segment_summaries: list[dict[str, object]] = []
    for contig_name, contig_values in sorted(contig_depths.items(), key=lambda item: item[0]):
        contig_total = len(contig_values)
        contig_points, contig_x_values = _compress_depths(contig_values, min(target_points, 240))
        contig_covered = sum(1 for depth in contig_values if depth > 0)
        contig_covered_10x = sum(1 for depth in contig_values if depth >= 10)
        contig_covered_100x = sum(1 for depth in contig_values if depth >= 100)
        segment_summaries.append({
            "name": contig_name,
            "label": contig_name,
            "points": contig_points,
            "x_values": contig_x_values,
            "x_ticks": [1, max(1, contig_total // 2), contig_total],
            "total_bases": contig_total,
            "max_depth": max(contig_values) if contig_values else None,
            "mean_depth": round(sum(contig_values) / contig_total, 2) if contig_total else None,
            "coverage_fraction": round(contig_covered / contig_total, 6) if contig_total else None,
            "coverage_10x_fraction": round(contig_covered_10x / contig_total, 6) if contig_total else None,
            "coverage_100x_fraction": round(contig_covered_100x / contig_total, 6) if contig_total else None,
        })
    return {
        "status": "ready",
        "label": "基因组覆盖度",
        "x_label": "基因组位置",
        "y_label": "测序深度",
        "points": points,
        "x_values": x_values,
        "total_bases": total_bases,
        "contig_count": len(contigs),
        "max_depth": max_depth,
        "mean_depth": mean_depth,
        "coverage_fraction": round(covered_bases / total_bases, 6) if total_bases else None,
        "coverage_10x_fraction": round(covered_10x_bases / total_bases, 6) if total_bases else None,
        "coverage_100x_fraction": round(covered_100x_bases / total_bases, 6) if total_bases else None,
        "x_ticks": x_ticks,
        "segments": segment_summaries,
    }

def _normalize_bandavirus_assembly_coverage(assembly_coverage: dict, serotype_result: dict) -> dict:
    if not isinstance(assembly_coverage, dict) or str(assembly_coverage.get("status") or "").strip() != "ready":
        return assembly_coverage
    coverage_segments = assembly_coverage.get("segments")
    if not isinstance(coverage_segments, list) or not coverage_segments:
        return assembly_coverage

    segment_order = {"L": 0, "M": 1, "S": 2}
    af_table = serotype_result.get("af_segment_typing") if isinstance(serotype_result.get("af_segment_typing"), dict) else {}
    cj_table = serotype_result.get("cj_segment_typing") if isinstance(serotype_result.get("cj_segment_typing"), dict) else {}
    af_columns = af_table.get("columns") if isinstance(af_table.get("columns"), list) else []
    af_rows = af_table.get("rows") if isinstance(af_table.get("rows"), list) else []
    cj_columns = cj_table.get("columns") if isinstance(cj_table.get("columns"), list) else []
    cj_rows = cj_table.get("rows") if isinstance(cj_table.get("rows"), list) else []

    def _table_row_map(columns: list, rows: list, key_name: str) -> dict[str, dict[str, str]]:
        if key_name not in columns:
            return {}
        key_index = columns.index(key_name)
        mapping: dict[str, dict[str, str]] = {}
        for row in rows:
            if not isinstance(row, list):
                continue
            key = str(row[key_index] if key_index < len(row) else "").strip()
            if not key:
                continue
            row_map = {
                str(columns[index]): str(row[index] if index < len(row) else "").strip()
                for index in range(len(columns))
            }
            mapping[key] = row_map
        return mapping

    af_by_segment = _table_row_map(af_columns, af_rows, "segment")
    cj_by_segment = _table_row_map(cj_columns, cj_rows, "segment")

    normalized_segments: list[dict] = []
    for raw_item in coverage_segments:
        if not isinstance(raw_item, dict):
            continue
        item = dict(raw_item)
        raw_name = str(item.get("name") or item.get("label") or "").strip()
        segment_symbol = ""
        accession = ""
        if raw_name:
            parts = raw_name.split("_")
            tail = str(parts[-1]).strip().upper() if parts else ""
            if tail in segment_order:
                segment_symbol = tail
                accession = "_".join(parts[:-1]).strip()
        if not segment_symbol:
            continue
        af_info = af_by_segment.get(segment_symbol, {})
        cj_info = cj_by_segment.get(segment_symbol, {})
        item["segment"] = segment_symbol
        item["name"] = segment_symbol
        item["label"] = f"{segment_symbol} 片段"
        item["accession"] = str(af_info.get("accession") or accession or "").strip()
        item["af_group"] = str(af_info.get("af_group") or "").strip()
        item["cj_group"] = str(cj_info.get("cj_group") or "").strip()
        item["sort_order"] = segment_order.get(segment_symbol, 99)
        normalized_segments.append(item)

    if len(normalized_segments) < 2:
        return assembly_coverage

    normalized_segments.sort(key=lambda item: (int(item.get("sort_order", 99)), str(item.get("segment") or "")))
    assembly_coverage = dict(assembly_coverage)
    assembly_coverage["segments"] = normalized_segments
    assembly_coverage["label"] = "Bandavirus L/M/S 覆盖度"
    assembly_coverage["segmented_view"] = "bandavirus_lms"
    return assembly_coverage

def _format_fraction_percent(value: object) -> str:
    numeric = _safe_float(value)
    if numeric is None:
        return "--"
    return f"{numeric * 100:.2f}%"

def _format_percent_value(value: object) -> str:
    numeric = _safe_float(value)
    if numeric is None:
        return "--"
    return f"{numeric:.2f}%"

def _read_norovirus_gene_coverage(report_dir: Path, sample_name: str) -> dict[str, str]:
    if not sample_name:
        return {}
    dual_typing = _read_tsv_rows(
        report_dir / f"{sample_name}_norovirus_reference_selection" / "typing" / "dual_typing.tsv"
    )
    columns = dual_typing.get("columns") or []
    rows = dual_typing.get("rows") or []
    if not columns or not rows:
        return {}
    gene_index = columns.index("gene") if "gene" in columns else -1
    coverage_index = columns.index("coverage") if "coverage" in columns else -1
    if gene_index < 0 or coverage_index < 0:
        return {}
    gene_coverage: dict[str, str] = {}
    for row in rows:
        if not isinstance(row, list):
            continue
        gene_name = str(row[gene_index] if gene_index < len(row) else "").strip().lower()
        coverage_value = row[coverage_index] if coverage_index < len(row) else ""
        if gene_name in {"rdrp", "vp1"}:
            gene_coverage[gene_name] = _format_percent_value(coverage_value)
    return gene_coverage

def _merge_assembly_coverage_summary(summary_table: dict, coverage_section: dict) -> dict:
    columns = list(summary_table.get("columns") or [])
    rows = [list(row) if isinstance(row, list) else row for row in (summary_table.get("rows") or [])]
    if not rows:
        return summary_table
    coverage_columns = [
        ("覆盖度", _format_fraction_percent(coverage_section.get("coverage_fraction"))),
        ("10x 百分比", _format_fraction_percent(coverage_section.get("coverage_10x_fraction"))),
        ("100x 百分比", _format_fraction_percent(coverage_section.get("coverage_100x_fraction"))),
    ]
    missing_columns = [label for label, _ in coverage_columns if label not in columns]
    if missing_columns:
        columns.extend(missing_columns)
        rows = [
            list(row) + [""] * len(missing_columns)
            for row in rows
        ]
    for row in rows:
        for label, value in coverage_columns:
            try:
                index = columns.index(label)
            except ValueError:
                continue
            if index >= len(row):
                row.extend([""] * (index + 1 - len(row)))
            row[index] = value
    return {"columns": columns, "rows": rows}

def _build_virus_fallback_assembly_summary(report_dir: Path, sample_name: str, coverage_section: dict) -> dict:
    fasta_path = _resolve_report_artifact_path(
        report_dir,
        [f"{sample_name}.final.fasta"] if sample_name else [],
        ["*.final.fasta", "*.consensus.fasta"],
    )
    fasta_summary = _read_fasta_assembly_summary(fasta_path) if fasta_path else {"columns": [], "rows": []}
    if not fasta_summary.get("rows"):
        return {"columns": [], "rows": []}
    first_row = fasta_summary["rows"][0]
    base_columns = fasta_summary.get("columns") or []
    value_by_column = {
        str(base_columns[index]): first_row[index] if index < len(first_row) else ""
        for index in range(len(base_columns))
    }
    norovirus_gene_coverage = _read_norovirus_gene_coverage(report_dir, sample_name)
    if norovirus_gene_coverage:
        return {
            "columns": ["样本名称", "Contig数量", "总长度(bp)", "最大片段长度(bp)", "N50长度(bp)", "RdRp覆盖度", "VP1覆盖度"],
            "rows": [[
                sample_name or str(value_by_column.get("结果文件") or fasta_path.name),
                str(value_by_column.get("Contig数量") or "--"),
                str(value_by_column.get("总长度(bp)") or "--"),
                str(value_by_column.get("最大Contig(bp)") or "--"),
                str(value_by_column.get("N50(bp)") or "--"),
                norovirus_gene_coverage.get("rdrp", "--"),
                norovirus_gene_coverage.get("vp1", "--"),
            ]],
        }
    return {
        "columns": ["样本名称", "Contig数量", "总长度(bp)", "最大片段长度(bp)", "N50长度(bp)"],
        "rows": [[
            sample_name or str(value_by_column.get("结果文件") or fasta_path.name),
            str(value_by_column.get("Contig数量") or "--"),
            str(value_by_column.get("总长度(bp)") or "--"),
            str(value_by_column.get("最大Contig(bp)") or "--"),
            str(value_by_column.get("N50(bp)") or "--"),
        ]],
    }
