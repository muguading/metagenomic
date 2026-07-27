from __future__ import annotations

import csv
import re
from pathlib import Path

from .admin_runtime import _resolve_workstation_module
from .parse_utils import _safe_float, _safe_int
from .phylogeny_reports import _attach_itol_payload, _read_pathosource_tree
from .table_io import _read_tsv_rows

def _is_pathosource_report_task(task: dict) -> bool:
    params = task.get("params") or {}
    workstation_key = _resolve_workstation_module(params.get("workstation_key")).get("key", "")
    pipeline_script = str(task.get("pipeline_script") or "").strip()
    return workstation_key == "pathosource" or pipeline_script == "PathoSource.py"

def _read_pathosource_distance_bins(path: Path) -> dict:
    raw = _read_tsv_rows(path)
    columns = raw.get("columns", [])
    rows = raw.get("rows", [])
    if not columns or not rows:
        return {"status": "empty", "columns": [], "rows": [], "bars": [], "total_pairs": 0}
    values = []
    for value in rows[0]:
        numeric = _safe_int(value)
        values.append(0 if numeric is None else numeric)
    total_pairs = sum(values)
    max_value = max(values) if values else 0
    bars = []
    for label, value in zip(columns, values):
        bars.append({
            "label": label,
            "value": value,
            "ratio": round((value / max_value) * 100, 2) if max_value else 0,
            "share": round((value / total_pairs) * 100, 2) if total_pairs else 0,
        })
    return {
        "status": "ready",
        "columns": columns,
        "rows": [values],
        "bars": bars,
        "total_pairs": total_pairs,
    }

def _read_pathosource_snp_matrix(path: Path, preview_size: int = 18) -> dict:
    if not path.is_file():
        return {
            "status": "empty",
            "summary": {},
            "preview": {"columns": [], "rows": []},
        }
    try:
        with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader, [])
            if len(header) < 2:
                return {"status": "empty", "summary": {}, "preview": {"columns": [], "rows": []}}
            sample_names = header[1:]
            preview_columns = ["样本"] + sample_names[:preview_size]
            preview_rows: list[list[object]] = []
            values: list[int] = []
            for row_index, row in enumerate(reader):
                if not row:
                    continue
                sample_name = row[0]
                numeric_cells = row[1:]
                if row_index < preview_size:
                    preview_rows.append([sample_name] + numeric_cells[:preview_size])
                for column_index, cell in enumerate(numeric_cells):
                    value = _safe_int(cell)
                    if value is None:
                        continue
                    if column_index < row_index:
                        values.append(value)
    except OSError:
        return {"status": "empty", "summary": {}, "preview": {"columns": [], "rows": []}}
    ordered_values = sorted(values)
    pair_count = len(ordered_values)
    positive_values = [value for value in ordered_values if value > 0]
    median_value = ordered_values[pair_count // 2] if pair_count else None
    return {
        "status": "ready",
        "summary": {
            "sample_count": len(sample_names),
            "pair_count": pair_count,
            "min_distance": positive_values[0] if positive_values else 0,
            "max_distance": ordered_values[-1] if ordered_values else 0,
            "median_distance": median_value if median_value is not None else 0,
            "zero_distance_pairs": pair_count - len(positive_values),
        },
        "preview": {
            "columns": preview_columns,
            "rows": preview_rows,
        },
    }

def _read_pathosource_ani(path: Path, preview_size: int = 16, top_pair_size: int = 20) -> dict:
    if not path.is_file():
        return {
            "status": "empty",
            "summary": {},
            "top_pairs": {"columns": [], "rows": []},
            "preview": {"columns": [], "rows": []},
        }
    try:
        with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            samples: list[str] = []
            sample_set: set[str] = set()
            pair_values: dict[tuple[str, str], float] = {}
            non_self_pairs: list[tuple[str, str, float, float | None, float | None]] = []
            for row in reader:
                ref_name = str(row.get("Ref_name") or "").strip()
                query_name = str(row.get("Query_name") or "").strip()
                ani_value = _safe_float(row.get("ANI"))
                align_ref = _safe_float(row.get("Align_fraction_ref"))
                align_query = _safe_float(row.get("Align_fraction_query"))
                if not ref_name or not query_name or ani_value is None:
                    continue
                if ref_name not in sample_set:
                    sample_set.add(ref_name)
                    samples.append(ref_name)
                if query_name not in sample_set:
                    sample_set.add(query_name)
                    samples.append(query_name)
                pair_values[(ref_name, query_name)] = ani_value
                if ref_name == query_name:
                    continue
                key = tuple(sorted((ref_name, query_name)))
                if (key[0], key[1]) in pair_values and key != (ref_name, query_name):
                    continue
                if ref_name < query_name:
                    non_self_pairs.append((ref_name, query_name, ani_value, align_ref, align_query))
    except OSError:
        return {
            "status": "empty",
            "summary": {},
            "top_pairs": {"columns": [], "rows": []},
            "preview": {"columns": [], "rows": []},
        }
    ani_values = sorted([item[2] for item in non_self_pairs])
    pair_count = len(ani_values)
    preview_samples = samples[:preview_size]
    preview_columns = ["样本"] + preview_samples
    preview_rows: list[list[object]] = []
    for sample in preview_samples:
        preview_row: list[object] = [sample]
        for target in preview_samples:
            value = pair_values.get((sample, target))
            if value is None:
                value = pair_values.get((target, sample))
            preview_row.append(f"{value:.2f}" if value is not None else "-")
        preview_rows.append(preview_row)
    top_pairs_rows = [
        [left, right, f"{ani:.2f}", f"{align_ref:.2f}" if align_ref is not None else "-", f"{align_query:.2f}" if align_query is not None else "-"]
        for left, right, ani, align_ref, align_query in sorted(non_self_pairs, key=lambda item: item[2], reverse=True)[:top_pair_size]
    ]
    return {
        "status": "ready",
        "summary": {
            "sample_count": len(samples),
            "pair_count": pair_count,
            "min_ani": round(ani_values[0], 2) if ani_values else None,
            "max_ani": round(ani_values[-1], 2) if ani_values else None,
            "median_ani": round(ani_values[pair_count // 2], 2) if ani_values else None,
        },
        "top_pairs": {
            "columns": ["样本A", "样本B", "ANI(%)", "Ref比对覆盖(%)", "Query比对覆盖(%)"],
            "rows": top_pairs_rows,
        },
        "preview": {
            "columns": preview_columns,
            "rows": preview_rows,
        },
    }

def _read_pathosource_mlst(path: Path) -> dict:
    if not path.is_file():
        return {
            "status": "empty",
            "columns": [],
            "rows": [],
            "summary": {"sample_count": 0, "scheme_count": 0, "st_count": 0, "dominant_st": "--"},
        }
    rows: list[list[str]] = []
    scheme_counter: dict[str, int] = {}
    st_counter: dict[str, int] = {}
    locus_count = 0
    locus_headers: list[str] = []
    try:
        with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            for row in reader:
                if len(row) < 3:
                    continue
                sample_name = re.sub(r"\.(raw\.)?f(ast)?a(sta)?$", "", str(row[0]).strip(), flags=re.IGNORECASE)
                scheme = str(row[1]).strip()
                st_value = str(row[2]).strip()
                loci = []
                raw_loci = row[3:]
                for index, locus in enumerate(raw_loci):
                    locus_text = str(locus or "").strip()
                    match = re.match(r"^\s*([^(]+?)\(([^()]*)\)\s*$", locus_text)
                    if len(locus_headers) <= index:
                        locus_headers.append(match.group(1).strip() if match else f"Locus{index + 1}")
                    loci.append(match.group(2).strip() if match else locus_text)
                locus_count = max(locus_count, len(loci))
                rows.append([sample_name, scheme, st_value, *loci])
                if scheme:
                    scheme_counter[scheme] = scheme_counter.get(scheme, 0) + 1
                if st_value:
                    st_counter[st_value] = st_counter.get(st_value, 0) + 1
    except OSError:
        return {
            "status": "empty",
            "columns": [],
            "rows": [],
            "summary": {"sample_count": 0, "scheme_count": 0, "st_count": 0, "dominant_st": "--"},
        }
    columns = ["样本名称", "MLST方案", "ST"] + [header or f"Locus{i}" for i, header in enumerate(locus_headers[:locus_count], start=1)]
    normalized_rows = [row + [""] * (len(columns) - len(row)) for row in rows]
    dominant_st = max(st_counter.items(), key=lambda item: item[1])[0] if st_counter else "--"
    return {
        "status": "ready" if normalized_rows else "empty",
        "columns": columns,
        "rows": normalized_rows,
        "summary": {
            "sample_count": len(normalized_rows),
            "scheme_count": len(scheme_counter),
            "st_count": len(st_counter),
            "dominant_st": dominant_st,
        },
    }

def _read_pathosource_mutations(path: Path, top_sample_size: int = 15, preview_size: int = 200) -> dict:
    if not path.is_file():
        return {
            "status": "empty",
            "summary": {
                "site_count": 0,
                "mutated_site_count": 0,
                "sample_count": 0,
                "mutated_sample_count": 0,
                "max_mutation_sample": "--",
                "max_mutation_count": 0,
            },
            "sample_bars": [],
            "table": {"columns": [], "rows": []},
        }
    try:
        with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader, [])
            if len(header) < 4:
                return {
                    "status": "empty",
                    "summary": {
                        "site_count": 0,
                        "mutated_site_count": 0,
                        "sample_count": 0,
                        "mutated_sample_count": 0,
                        "max_mutation_sample": "--",
                        "max_mutation_count": 0,
                    },
                    "sample_bars": [],
                    "table": {"columns": [], "rows": []},
                }
            sample_names = [str(item or "").strip() for item in header[3:]]
            sample_mutation_counter = {sample: 0 for sample in sample_names if sample}
            mutation_rows: list[list[object]] = []
            site_count = 0
            mutated_site_count = 0
            for row in reader:
                if len(row) < 3:
                    continue
                site_count += 1
                chrom = str(row[0] or "").strip()
                position = str(row[1] or "").strip()
                ref_base = str(row[2] or "").strip()
                mutated_samples: list[str] = []
                alt_counts: dict[str, int] = {}
                for sample_name, base in zip(sample_names, row[3:]):
                    sample = str(sample_name or "").strip()
                    observed = str(base or "").strip()
                    if not sample or not observed or observed == "-" or observed == ref_base:
                        continue
                    mutated_samples.append(sample)
                    sample_mutation_counter[sample] = sample_mutation_counter.get(sample, 0) + 1
                    alt_counts[observed] = alt_counts.get(observed, 0) + 1
                if not mutated_samples:
                    continue
                mutated_site_count += 1
                alt_summary = " / ".join(f"{allele}×{count}" for allele, count in sorted(alt_counts.items(), key=lambda item: (-item[1], item[0])))
                preview = "，".join(mutated_samples[:8])
                if len(mutated_samples) > 8:
                    preview = f"{preview} 等 {len(mutated_samples)} 个样本"
                mutation_rows.append([
                    chrom,
                    position,
                    ref_base,
                    len(mutated_samples),
                    alt_summary or "-",
                    preview or "-",
                ])
    except OSError:
        return {
            "status": "empty",
            "summary": {
                "site_count": 0,
                "mutated_site_count": 0,
                "sample_count": 0,
                "mutated_sample_count": 0,
                "max_mutation_sample": "--",
                "max_mutation_count": 0,
            },
            "sample_bars": [],
            "table": {"columns": [], "rows": []},
        }
    sample_items = sorted(sample_mutation_counter.items(), key=lambda item: (-item[1], item[0]))
    max_mutation_sample, max_mutation_count = sample_items[0] if sample_items else ("--", 0)
    mutated_sample_count = sum(1 for _, count in sample_items if count > 0)
    max_bar_value = max((count for _, count in sample_items), default=0)
    sample_bars = [
        {
            "label": sample,
            "value": count,
            "ratio": round((count / max_bar_value) * 100, 2) if max_bar_value else 0,
        }
        for sample, count in sample_items[:top_sample_size]
        if count > 0
    ]
    return {
        "status": "ready" if mutation_rows else "empty",
        "summary": {
            "site_count": site_count,
            "mutated_site_count": mutated_site_count,
            "sample_count": len(sample_names),
            "mutated_sample_count": mutated_sample_count,
            "max_mutation_sample": max_mutation_sample,
            "max_mutation_count": max_mutation_count,
        },
        "sample_bars": sample_bars,
        "table": {
            "columns": ["染色体", "变异位点位置", "参考位点碱基", "突变样本数", "替代碱基统计", "发生突变的样本"],
            "rows": mutation_rows[:preview_size],
        },
    }

def _build_pathosource_report_payload(
    *,
    task: dict,
    report_dir: Path,
    report_source: dict,
    sample_name: str,
    sample_display_name: str,
) -> dict:
    params = task.get("params") or {}
    cluster_table = _read_tsv_rows(report_dir / "Cluster.tsv")
    distance_bins = _read_pathosource_distance_bins(report_dir / "dis_bin.tsv")
    snp_matrix = _read_pathosource_snp_matrix(report_dir / "dis.mat.txt")
    ani_section = _read_pathosource_ani(report_dir / "Full_ANI.txt")
    mlst_section = _read_pathosource_mlst(report_dir / "mlst.txt")
    mutation_section = _read_pathosource_mutations(report_dir / "Mutate.tsv")
    grapetree_tree = _read_pathosource_tree(report_dir / "grapetree.nwk", "GrapeTree 最小生成树")
    mlst_tree = _read_pathosource_tree(report_dir / "mlst.nwk", "MLST 进化树")
    core_tree = _attach_itol_payload(
        _read_pathosource_tree(report_dir / "rmref.core.aln.contree", "核心 SNP 系统发育树"),
        report_dir / "rmref.core.aln.contree",
    )

    cluster_rows = cluster_table.get("rows", [])
    sample_total = max(
        snp_matrix.get("summary", {}).get("sample_count") or 0,
        mlst_section.get("summary", {}).get("sample_count") or 0,
    )
    cluster_sizes = []
    for row in cluster_rows:
        if len(row) > 1:
            size_value = _safe_int(row[1])
            if size_value is not None:
                cluster_sizes.append(size_value)
    largest_cluster = max(cluster_sizes) if cluster_sizes else 0
    singleton_clusters = sum(1 for size in cluster_sizes if size <= 1)
    tree_ready_count = sum(1 for item in [grapetree_tree, mlst_tree, core_tree] if item.get("status") == "ready")

    payload = {
        "task": {
            "id": task.get("id"),
            "name": task.get("name"),
            "status": task.get("status"),
            "owner": task.get("owner"),
            "group": task.get("owner_group", ""),
            "created_at": task.get("created_at"),
            "started_at": task.get("started_at"),
            "finished_at": task.get("finished_at"),
            "input_path": params.get("input_path", ""),
            "output_dir": params.get("output_dir", ""),
            "asm_type": params.get("msamethod", "") or "溯源分析",
            "method": params.get("treemethod", "") or "PathoSource",
            "analysis_target": "bacteria",
            "species": params.get("species", ""),
            "sample_name": sample_name,
            "sample_display_name": sample_display_name or task.get("name") or "溯源进化树",
            "samples": report_source.get("samples", []),
            "report_mode": report_source.get("mode", "single"),
            "report_kind": "pathosource_phylogeny",
            "workstation_key": "pathosource",
        },
        "overview_metrics": [
            {"key": "phylo_sample_count", "label": "纳入样本", "type": "single", "display": str(sample_total or "--")},
            {"key": "phylo_cluster_count", "label": "成簇数量", "type": "single", "display": str(len(cluster_rows) or "--")},
            {"key": "phylo_largest_cluster", "label": "最大簇规模", "type": "single", "display": str(largest_cluster or "--")},
            {"key": "phylo_pair_count", "label": "距离比较对", "type": "single", "display": str(snp_matrix.get('summary', {}).get('pair_count') or "--")},
            {"key": "phylo_tree_ready", "label": "可视化树数", "type": "single", "display": str(tree_ready_count or "--")},
        ],
        "sections": {
            "pathosource": {
                "status": "ready",
                "cluster": {
                    "table": cluster_table,
                    "summary": {
                        "cluster_count": len(cluster_rows),
                        "largest_cluster": largest_cluster,
                        "singleton_clusters": singleton_clusters,
                    },
                },
                "distance_bins": distance_bins,
                "snp_matrix": snp_matrix,
                "ani": ani_section,
                "mlst": mlst_section,
                "mutations": mutation_section,
                "trees": {
                    "grapetree": grapetree_tree,
                    "mlst": mlst_tree,
                    "core": core_tree,
                },
            },
        },
    }
    return payload
