from __future__ import annotations

import csv
import gzip
import json
import re
from pathlib import Path

from .community_reports import _cell_at, _find_column_index, _first_mapping_value, _first_row_mapping, _first_table_cell
from .knowledge_interpretation import SALMONELLA_SEROVAR_ALIAS_MAP, _normalize_serotype_lookup_text
from .parse_utils import _safe_float, _safe_int
from .report_formatters import _coerce_percent_value, _display_percent, _display_percent_points, _human_bp, _human_count
from .table_io import _read_tsv_rows

def _read_summary_metrics(path: Path) -> dict:
    if not path.is_file():
        return {}
    total_sum_len = 0
    total_reads = 0
    weighted_q20 = 0.0
    weighted_q30 = 0.0
    try:
        with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                value = _safe_int(row.get("sum_len"))
                reads_value = _safe_int(
                    row.get("sum_reads")
                    or row.get("reads")
                    or row.get("read_count")
                    or row.get("序列数量")
                    or row.get("sum_num")
                )
                if reads_value is not None:
                    total_reads += reads_value
                if value is not None:
                    total_sum_len += value
                    q20 = _safe_float(row.get("Q20(%)"))
                    q30 = _safe_float(row.get("Q30(%)"))
                    if q20 is not None:
                        weighted_q20 += q20 * value
                    if q30 is not None:
                        weighted_q30 += q30 * value
    except OSError:
        return {}
    return {
        "sum_len": total_sum_len or None,
        "sum_reads": total_reads or None,
        "q20_rate": round(weighted_q20 / total_sum_len, 2) if total_sum_len else None,
        "q30_rate": round(weighted_q30 / total_sum_len, 2) if total_sum_len else None,
    }

def _read_checkm_metrics(path: Path) -> dict:
    if not path.is_file():
        return {}
    try:
        with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="	")
            first = next(reader, None)
            if not first:
                return {}
            return {
                "contamination": first.get("污染率") or first.get("contamination") or None,
                "completeness": first.get("完整性") or first.get("completeness") or None,
                "species_name": first.get("物种名称") or first.get("species_name") or None,
                "mlst_species_name": first.get("mlst 物种名称") or first.get("mlst_species_name") or None,
            }
    except OSError:
        return {}

def _looks_like_report_directory(path: Path) -> bool:
    if not path.is_dir():
        return False
    marker_patterns = (
        "summary.tsv",
        "community_summary.json",
        "community_command_plan.tsv",
        "community_metadata_preview.tsv",
        "community_demux_summary.tsv",
        "community_taxonomy_preview.tsv",
        "Assem_info.tsv",
        "Assem_info1.tsv",
        "*.checkm.tsv",
        "*.fastp2.json",
        "*_2.list.txt",
        "*.mlst_Stat.txt",
        "Cluster.tsv",
        "dis_bin.tsv",
        "dis.mat.txt",
        "Full_ANI.txt",
        "grapetree.nwk",
        "mlst.txt",
        "rmref.core.aln.contree",
    )
    for pattern in marker_patterns:
        if "*" in pattern:
            if any(path.glob(pattern)):
                return True
        elif (path / pattern).exists():
            return True
    if (path / "flye_output").is_dir() or (path / "checkm2_out").is_dir():
        return True
    if (path / "nextclade_output" / "nextclade.tsv").is_file():
        return True
    return False

def _sample_sort_key(path: Path, preferred_name: str = "") -> tuple[int, int, str]:
    if preferred_name and path.name == preferred_name:
        preferred_rank = 0
    else:
        preferred_rank = 1
    ready_rank = 0 if _looks_like_report_directory(path) else 1
    return (preferred_rank, ready_rank, path.name.lower())

def _resolve_report_source(task: dict, selected_sample: str = "") -> dict:
    output_dir = str((task.get("params") or {}).get("output_dir", "")).strip()
    preferred_sample = str((task.get("params") or {}).get("sample_name", "")).strip()
    task_name = str(task.get("name") or task.get("id") or "当前任务").strip()

    if not output_dir:
        return {
            "available": False,
            "reason": f"{task_name} 缺少服务器输出目录 output_dir，无法定位结果。",
            "mode": "single",
            "root_dir": None,
            "report_dir": None,
            "samples": [],
            "selected_sample": "",
        }

    try:
        output_root = Path(output_dir).expanduser().resolve()
    except OSError:
        return {
            "available": False,
            "reason": f"{task_name} 的服务器输出目录无效: {output_dir}",
            "mode": "single",
            "root_dir": None,
            "report_dir": None,
            "samples": [],
            "selected_sample": "",
        }

    pipeline_script = str(task.get("pipeline_script") or "").strip()
    if str(task.get("demo_type") or "").strip() and pipeline_script.startswith("demo_data/"):
        demo_report_dir = (Path(__file__).resolve().parent.parent / pipeline_script).resolve()
        if demo_report_dir.is_dir() and _looks_like_report_directory(demo_report_dir):
            return {
                "available": True,
                "mode": "single",
                "root_dir": demo_report_dir,
                "report_dir": demo_report_dir,
                "samples": [],
                "selected_sample": "",
            }
        if pipeline_script == "demo_data/meta":
            meta_demo_report_dir = (Path(__file__).resolve().parent.parent / "demo_data" / "meta_1").resolve()
            if meta_demo_report_dir.is_dir() and _looks_like_report_directory(meta_demo_report_dir):
                return {
                    "available": True,
                    "mode": "single",
                    "root_dir": meta_demo_report_dir,
                    "report_dir": meta_demo_report_dir,
                    "samples": [],
                    "selected_sample": "",
                }

    if not output_root.is_dir():
        return {
            "available": False,
            "reason": f"{task_name} 的服务器结果目录不存在: {output_root}",
            "mode": "single",
            "root_dir": output_root,
            "report_dir": None,
            "samples": [],
            "selected_sample": "",
        }

    params = task.get("params") or {}
    workstation_key = str(params.get("workstation_key") or "").strip().lower()
    if not _looks_like_report_directory(output_root) and (
        task_name == "demo_tree"
        or workstation_key == "pathosource"
        or Path(pipeline_script).name == "PathoSource.py"
    ):
        try:
            input_root = Path(str(params.get("input_path") or "")).expanduser().resolve()
        except OSError:
            input_root = Path("")
        if input_root.is_dir() and _looks_like_report_directory(input_root):
            return {
                "available": True,
                "mode": "single",
                "root_dir": input_root,
                "report_dir": input_root,
                "samples": [],
                "selected_sample": "",
            }

    if not _looks_like_report_directory(output_root) and task_name == "demo_meta":
        meta_demo_report_dir = (Path(__file__).resolve().parent.parent / "demo_data" / "meta_1").resolve()
        if meta_demo_report_dir.is_dir() and _looks_like_report_directory(meta_demo_report_dir):
            return {
                "available": True,
                "mode": "single",
                "root_dir": meta_demo_report_dir,
                "report_dir": meta_demo_report_dir,
                "samples": [],
                "selected_sample": "",
            }

    fastq_analysis_roots = []
    if output_root.name == "fastq_analysis":
        fastq_analysis_roots.append(output_root)
    fastq_analysis_roots.append(output_root / "fastq_analysis")

    for fastq_analysis_root in fastq_analysis_roots:
        if not fastq_analysis_root.is_dir():
            continue
        sample_dirs = sorted(
            [child for child in fastq_analysis_root.iterdir() if child.is_dir() and not child.name.startswith(".")],
            key=lambda item: item.name.lower(),
        )
        if sample_dirs:
            explicit_selected_sample = str(selected_sample or "").strip()
            preferred_name = explicit_selected_sample or preferred_sample
            ready_dirs = [item for item in sample_dirs if _looks_like_report_directory(item)]
            selected_dir = None
            if explicit_selected_sample:
                selected_dir = next((item for item in sample_dirs if item.name == explicit_selected_sample), None)
            if selected_dir is None and preferred_sample:
                preferred_dir = next((item for item in ready_dirs if item.name == preferred_sample), None)
                if preferred_dir is not None:
                    selected_dir = preferred_dir
            if selected_dir is None and ready_dirs:
                selected_dir = sorted(ready_dirs, key=lambda item: _sample_sort_key(item, preferred_name))[0]
            if selected_dir is None:
                selected_dir = sorted(sample_dirs, key=lambda item: _sample_sort_key(item, preferred_name))[0]
            if not ready_dirs and str(task.get("status") or "").upper() in {"QUEUED", "RUNNING", "PAUSED"}:
                return {
                    "available": False,
                    "reason": f"{task_name} 已创建 {len(sample_dirs)} 个样本目录，但还没有样本产出可展示结果文件。",
                    "mode": "multi",
                    "root_dir": fastq_analysis_root,
                    "report_dir": None,
                    "samples": [item.name for item in sample_dirs],
                    "selected_sample": "",
                    "task_params": params,
                }
            return {
                "available": True,
                "mode": "multi",
                "root_dir": fastq_analysis_root,
                "report_dir": selected_dir,
                "samples": [item.name for item in sample_dirs],
                "selected_sample": selected_dir.name,
                "task_params": params,
            }

    if _looks_like_report_directory(output_root):
        return {
            "available": True,
            "mode": "single",
            "root_dir": output_root,
            "report_dir": output_root,
            "samples": [],
            "selected_sample": "",
        }

    return {
        "available": False,
        "reason": f"{task_name} 的服务器输出目录存在，但其中没有识别到结果文件: {output_root}",
        "mode": "single",
        "root_dir": output_root,
        "report_dir": None,
        "samples": [],
        "selected_sample": "",
    }

def _resolve_report_sample_name(task: dict, report_dir: Path) -> str:
    def _has_sample_artifacts(name: str) -> bool:
        sample = str(name or "").strip()
        if not sample:
            return False
        candidate_paths = [
            report_dir / f"{sample}.fastp2.json",
            report_dir / f"{sample}.final.fasta",
            report_dir / f"{sample}.assemble.result.tsv",
            report_dir / f"{sample}_serotype_result.tsv",
            report_dir / f"{sample}.mlst_Stat.txt",
            report_dir / f"{sample}_2.list.txt",
            report_dir / f"{sample}_2.list2.txt",
            report_dir / f"{sample}_ngs.per-base.bed.gz",
            report_dir / f"{sample}_ngs.per-base.bed",
            report_dir / f"{sample}_rsv_reference_selection" / "selection.tsv",
            report_dir / f"{sample}_denv_reference_selection" / "selection.tsv",
            report_dir / f"{sample}_hpiv_reference_selection" / "selection.tsv",
            report_dir / f"{sample}_hadv_reference_selection" / "selection.tsv",
            report_dir / f"{sample}_enterovirus_reference_selection" / "selection.tsv",
            report_dir / f"{sample}_hepatovirus_reference_selection" / "selection.tsv",
            report_dir / f"{sample}_bandavirus_reference_selection" / "selection.tsv",
            report_dir / f"{sample}_orthohantavirus_reference_selection" / "selection.tsv",
            report_dir / f"{sample}_astroviridae_reference_selection" / "selection.tsv",
            report_dir / f"{sample}_rhinovirus_reference_selection" / "selection.tsv",
        ]
        return any(path.exists() for path in candidate_paths)

    params = task.get("params") or {}
    explicit_sample_name = str(params.get("sample_name") or "").strip()
    if explicit_sample_name and _has_sample_artifacts(explicit_sample_name):
        return explicit_sample_name
    task_name = str(task.get("name") or params.get("task_name") or "").strip()
    input_path = str(params.get("input_path") or "").strip()
    input_name = Path(input_path).name if input_path else ""
    if input_name == "fastq" and (report_dir / "Men-IGT.fastp2.json").is_file():
        return "Men-IGT"
    if input_name == "meta" and (report_dir / "A01.fastp2.json").is_file():
        return "A01"
    for pattern in ("*.checkm.tsv", "*.fastp2.json", "*_bacgenome.html"):
        matches = sorted(report_dir.glob(pattern))
        if matches:
            name = matches[0].name
            if name.endswith(".checkm.tsv"):
                return name.removesuffix(".checkm.tsv")
            if name.endswith(".fastp2.json"):
                return name.removesuffix(".fastp2.json")
            if name.endswith("_bacgenome.html"):
                return name.removesuffix("_bacgenome.html")
    if explicit_sample_name:
        return explicit_sample_name
    if (report_dir / "nextclade_output" / "nextclade.tsv").is_file() and input_name:
        parts = Path(input_name).name.split(".")
        return parts[0] if parts else Path(input_name).stem
    if task_name and task_name not in {"demo_fastq", "demo_meta", "demo_tree"}:
        return task_name
    return Path(input_name).stem if input_name else ""

def _resolve_report_sample_display_name(task: dict, sample_name: str) -> str:
    params = task.get("params") or {}
    input_path = str(params.get("input_path") or "").strip()
    input_name = Path(input_path).name if input_path else ""
    task_name = str(task.get("name") or params.get("task_name") or "").strip()
    if task_name == "demo_meta" or (input_name == "meta" and sample_name == "A01"):
        return "A01"
    return sample_name

def _build_multi_sample_queue_summary(report_source: dict) -> dict:
    if str(report_source.get("mode") or "").strip() != "multi":
        return {}
    root_dir = report_source.get("root_dir")
    samples = [str(item or "").strip() for item in (report_source.get("samples") or []) if str(item or "").strip()]
    if not isinstance(root_dir, Path) or not samples:
        return {}

    sample_rows: list[dict] = []
    total_bases = 0
    total_reads = 0
    weighted_q20 = 0.0
    weighted_q30 = 0.0
    contig_count_total = 0
    plasmid_count_total = 0
    assembly_length_total = 0
    completeness_values: list[float] = []
    contamination_values: list[float] = []
    species_counts: dict[str, int] = {}
    ready_count = 0
    params = report_source.get("task_params") if isinstance(report_source.get("task_params"), dict) else {}
    analysis_target = str(params.get("analysis_target") or "").strip().lower()
    is_virus = analysis_target == "virus" or str(params.get("workstation_key") or "").strip().lower() == "virus"

    for sample in samples:
        report_dir = root_dir / sample
        if not report_dir.is_dir():
            continue
        summary_info = _read_summary_metrics(report_dir / "summary.tsv")
        checkm_info = _read_checkm_metrics(report_dir / f"{sample}.checkm.tsv")
        assembly_profile = _read_multi_sample_assembly_profile(report_dir, sample)
        has_result = any([
            _looks_like_report_directory(report_dir),
            (report_dir / "summary.tsv").is_file(),
            (report_dir / f"{sample}.fastp2.json").is_file(),
            any(report_dir.glob("*.fastp2.json")),
            (report_dir / f"{sample}.checkm.tsv").is_file(),
            (report_dir / "Assem_info.tsv").is_file(),
            (report_dir / "nextclade_output" / "nextclade.tsv").is_file(),
        ])
        if has_result:
            ready_count += 1

        bases = _safe_int(summary_info.get("sum_len")) or 0
        reads = _safe_int(summary_info.get("sum_reads")) or 0
        total_bases += bases
        total_reads += reads
        q20 = _safe_float(summary_info.get("q20_rate"))
        q30 = _safe_float(summary_info.get("q30_rate"))
        if bases and q20 is not None:
            weighted_q20 += q20 * bases
        if bases and q30 is not None:
            weighted_q30 += q30 * bases

        contigs = _safe_int(assembly_profile.get("contig_count")) or 0
        plasmids = _safe_int(assembly_profile.get("plasmid_count")) or 0
        assembly_length = _safe_int(assembly_profile.get("total_length")) or 0
        contig_count_total += contigs
        plasmid_count_total += plasmids
        assembly_length_total += assembly_length

        completeness = _safe_float(checkm_info.get("completeness"))
        contamination = _safe_float(checkm_info.get("contamination"))
        if completeness is not None:
            completeness_values.append(completeness)
        if contamination is not None:
            contamination_values.append(contamination)
        taxonomy_call = _read_multi_sample_taxonomy_call(report_dir, sample)
        typing_call = _read_multi_sample_typing_call(report_dir, sample, checkm_info, is_virus=is_virus)
        species_name = str(
            typing_call.get("species")
            or checkm_info.get("species_name")
            or taxonomy_call.get("species")
            or checkm_info.get("mlst_species_name")
            or (params.get("species") if is_virus and typing_call else "")
            or ""
        ).strip()
        if species_name:
            species_counts[species_name] = species_counts.get(species_name, 0) + 1
        resistance_count = len(_read_tsv_rows(report_dir / "Assem_abricate_CARD.tsv").get("rows") or [])
        virulence_count = len(_read_tsv_rows(report_dir / "Assem_abricate_VFDB.tsv").get("rows") or [])

        sample_rows.append({
            "sample": sample,
            "ready": has_result,
            "total_bases": bases or None,
            "total_reads": reads or None,
            "q20_rate": q20,
            "q30_rate": q30,
            "contig_count": contigs or None,
            "plasmid_count": plasmids or None,
            "assembly_length": assembly_length or None,
            "completeness": completeness,
            "contamination": contamination,
            "species_name": species_name,
            "taxonomy_ratio": taxonomy_call.get("ratio"),
            "typing": typing_call.get("typing") or "",
            "nextclade_typing": typing_call.get("nextclade_typing") or "-",
            "serotype": typing_call.get("serotype") or "",
            "coverage": _display_percent(typing_call.get("coverage")) if typing_call.get("coverage") not in {None, ""} else "",
            "coverage_1x": typing_call.get("coverage_1x") or (_display_percent(typing_call.get("coverage")) if typing_call.get("coverage") not in {None, ""} else ""),
            "coverage_10x": typing_call.get("coverage_10x") or "",
            "coverage_100x": typing_call.get("coverage_100x") or "",
            "mean_depth": typing_call.get("mean_depth") or "",
            "qc_status": typing_call.get("qc_status") or "",
            "resistance_count": resistance_count or None,
            "virulence_count": virulence_count or None,
            "note": typing_call.get("note") or "",
        })

    def _average(values: list[float]) -> float | None:
        return round(sum(values) / len(values), 2) if values else None

    species_rank = [
        {"name": name, "count": count}
        for name, count in sorted(species_counts.items(), key=lambda item: (-item[1], item[0].lower()))
    ]
    return {
        "sample_count": len(samples),
        "ready_count": ready_count,
        "total_bases": total_bases or None,
        "total_reads": total_reads or None,
        "q20_rate": round(weighted_q20 / total_bases, 2) if total_bases else None,
        "q30_rate": round(weighted_q30 / total_bases, 2) if total_bases else None,
        "contig_count_total": contig_count_total or None,
        "plasmid_count_total": plasmid_count_total or None,
        "assembly_length_total": assembly_length_total or None,
        "avg_completeness": _average(completeness_values),
        "avg_contamination": _average(contamination_values),
        "species_rank": species_rank[:5],
        "samples": sample_rows,
        "analysis_target": "virus" if is_virus else "bacteria",
        "table": _build_multi_sample_overview_table(sample_rows, is_virus=is_virus),
    }

def _read_multi_sample_taxonomy_call(report_dir: Path, sample_name: str) -> dict:
    table = _read_tsv_rows(report_dir / f"{sample_name}_2.list.txt")
    columns = table.get("columns") or []
    rows = table.get("rows") or []
    if not columns or not rows:
        return {}
    species_index = _find_column_index(columns, ["种", "species", "scientific_name", "name"])
    count_index = _find_column_index(columns, ["序列数量", "reads", "read_count", "count"])
    ratio_index = _find_column_index(columns, ["比例", "比例数值", "ratio", "abundance", "relative_abundance"])
    best: dict[str, object] = {}
    best_score = (-1.0, -1.0)
    for row in rows:
        if not isinstance(row, list):
            continue
        species = _cell_at(row, species_index)
        if not species or species == "-":
            continue
        reads = _safe_float(_cell_at(row, count_index)) or 0.0
        ratio = _coerce_percent_value(_cell_at(row, ratio_index))
        ratio_score = ratio if ratio is not None else 0.0
        if (reads, ratio_score) > best_score:
            best_score = (reads, ratio_score)
            best = {"species": species, "reads": int(reads) if reads else None, "ratio": ratio}
    return best

def _read_multi_sample_assembly_profile(report_dir: Path, sample_name: str) -> dict:
    profile = _read_assembly_profile(report_dir / "Assem_info.tsv")
    if _safe_int(profile.get("contig_count")) or _safe_int(profile.get("total_length")):
        return profile
    fasta_candidates = [
        report_dir / f"{sample_name}.final.fasta",
        report_dir / f"{sample_name}.consensus.fasta",
        report_dir / "tmp_combine.fa",
    ]
    fasta_path = next((path for path in fasta_candidates if path.is_file()), Path(""))
    fasta_summary = _read_fasta_assembly_summary(fasta_path)
    columns = fasta_summary.get("columns") or []
    rows = fasta_summary.get("rows") or []
    first = rows[0] if rows and isinstance(rows[0], list) else []
    return {
        "contig_count": _first_table_cell(first, columns, ["Contig数量"]),
        "plasmid_count": None,
        "total_count": _first_table_cell(first, columns, ["Contig数量"]),
        "total_length": _first_table_cell(first, columns, ["总长度(bp)"]),
    }


def _read_multi_sample_influenza_typing_call(report_dir: Path, sample_name: str) -> dict[str, str]:
    """Read the subtype selected by the workflow's influenza reference stage."""
    table = _read_tsv_rows(report_dir / "wf_flu" / "typing_summary.tsv")
    row = _first_row_mapping(table)
    if not row:
        return {}
    status = _first_mapping_value(row, ["status", "状态"])
    subtype_call = _first_mapping_value(row, ["subtype_call", "分型结果", "subtype"])
    ha_subtype = _first_mapping_value(row, ["ha_subtype", "HA亚型"])
    na_subtype = _first_mapping_value(row, ["na_subtype", "NA亚型"])
    if not _has_meaningful_serotype_value(subtype_call):
        subtype_call = "".join(
            value for value in (ha_subtype, na_subtype)
            if _has_meaningful_serotype_value(value)
        )
    if not _has_meaningful_serotype_value(subtype_call):
        return {}
    return {
        "species": _first_mapping_value(row, ["influenza_type", "流感类型", "virus_type"]),
        "typing": subtype_call,
        "status": status,
    }


def _read_multi_sample_influenza_nextclade_typing(report_dir: Path) -> str:
    """Return ready HA and NA Nextclade clades as ``HA|NA`` for batch display."""
    table = _read_tsv_rows(report_dir / "wf_flu" / "nextclade" / "segment_analysis.tsv")
    columns = table.get("columns") or []
    rows = table.get("rows") or []
    segment_index = _find_column_index(columns, ["segment", "片段"])
    clade_index = _find_column_index(columns, ["clade", "分支", "Nextclade分型"])
    status_index = _find_column_index(columns, ["status", "状态"])
    if segment_index is None or clade_index is None:
        return "-"

    values_by_segment: dict[str, list[str]] = {"HA": [], "NA": []}
    for row in rows:
        if not isinstance(row, list):
            continue
        segment = _cell_at(row, segment_index).upper()
        if segment not in values_by_segment:
            continue
        status = _cell_at(row, status_index).lower() if status_index is not None else "ready"
        clade = _cell_at(row, clade_index)
        if status != "ready" or not _has_meaningful_serotype_value(clade):
            continue
        if clade not in values_by_segment[segment]:
            values_by_segment[segment].append(clade)

    ha = "/".join(values_by_segment["HA"])
    na = "/".join(values_by_segment["NA"])
    return f"{ha or '-'}|{na or '-'}" if ha or na else "-"


def _read_multi_sample_typing_call(report_dir: Path, sample_name: str, checkm_info: dict, *, is_virus: bool) -> dict:
    result: dict[str, object] = {}
    serotype_table = _read_multi_sample_serotype_table(report_dir, sample_name)
    serotype_row = _first_row_mapping(serotype_table)
    if is_virus:
        nextclade_row = _first_row_mapping(_read_tsv_rows(report_dir / "nextclade_output" / "nextclade.tsv"))
        influenza_typing = _read_multi_sample_influenza_typing_call(report_dir, sample_name)
        influenza_nextclade_typing = _read_multi_sample_influenza_nextclade_typing(report_dir)
        species = _first_mapping_value(serotype_row, ["病毒类型", "物种", "species", "virus_type"])
        clade = _first_mapping_value(serotype_row, ["Nextclade分型", "大类分型", "大亚型", "分型结果", "基因型", "亚型", "clade", "type"])
        lineage = _first_mapping_value(serotype_row, ["Pango谱系", "S子亚型", "子亚型", "G分型", "P分型", "组合分型", "lineage", "subtype"])
        if influenza_typing:
            species = influenza_typing.get("species") or species
            clade = influenza_typing.get("typing") or clade
            lineage = ""
        hiv_typing = _read_multi_sample_hiv_typing_call(report_dir, sample_name)
        if not species:
            species = str(hiv_typing.get("species") or "").strip()
        if not clade:
            clade = str(hiv_typing.get("clade") or "").strip()
        if not lineage:
            lineage = str(hiv_typing.get("lineage") or "").strip()
        if not clade:
            clade = _first_mapping_value(nextclade_row, ["clade", "clade_display", "clade_who", "Nextclade_pango", "lineage", "genotype"])
        if not lineage:
            lineage = _first_mapping_value(nextclade_row, ["lineage", "genotype", "Nextclade_pango", "serotype"])
        coverage = _first_mapping_value(serotype_row, ["覆盖度", "全长覆盖度", "coverage"])
        if not coverage:
            coverage = _first_mapping_value(nextclade_row, ["coverage", "coverageScore"])
        depth_coverage = _read_multi_sample_depth_coverage(report_dir, sample_name)
        mean_depth = _first_mapping_value(serotype_row, ["平均深度", "mean_depth", "depth"])
        if not mean_depth:
            mean_depth = _read_multi_sample_mean_depth(report_dir, sample_name)
        typing_parts = []
        for item in (clade, lineage):
            item_text = str(item or "").strip()
            if item_text and item_text != "-" and item_text not in typing_parts:
                typing_parts.append(item_text)
        result.update({
            "species": species,
            "typing": " / ".join(typing_parts),
            "nextclade_typing": influenza_nextclade_typing,
            "coverage": coverage,
            "coverage_1x": depth_coverage.get("coverage_1x") or (_display_percent(coverage) if coverage else ""),
            "coverage_10x": depth_coverage.get("coverage_10x") or "",
            "coverage_100x": depth_coverage.get("coverage_100x") or "",
            "mean_depth": mean_depth,
            "qc_status": _first_mapping_value(serotype_row, ["QC状态", "qc_status"]) or _first_mapping_value(nextclade_row, ["qc.overallStatus"]),
            "note": _first_mapping_value(serotype_row, ["说明", "note"]),
        })
        return result

    mlst_row = _first_row_mapping(_read_tsv_rows(report_dir / f"{sample_name}.mlst_Stat.txt"))
    result.update({
        "species": str(checkm_info.get("species_name") or checkm_info.get("mlst_species_name") or "").strip(),
        "typing": _first_mapping_value(mlst_row, ["序列分型(ST)", "ST", "mlst_st"]),
        "serotype": _format_multi_sample_bacterial_serotype(serotype_row),
        "note": _first_mapping_value(serotype_row, ["说明", "note"]),
    })
    return result


def _read_multi_sample_hiv_typing_call(report_dir: Path, sample_name: str) -> dict[str, str]:
    candidates = [
        report_dir / f"{sample_name}_hiv_subtyping.json",
    ]
    candidates.extend(sorted(report_dir.glob("*_hiv_subtyping.json")))
    seen: set[Path] = set()
    for path in candidates:
        resolved = path.resolve()
        if resolved in seen or not path.is_file():
            continue
        seen.add(resolved)
        try:
            payload = json.loads(path.read_text(encoding="utf-8", errors="ignore"))
        except (OSError, json.JSONDecodeError):
            continue
        samples = payload.get("samples") if isinstance(payload, dict) else None
        if not isinstance(samples, list) or not samples:
            continue
        sample_payload = samples[0] if isinstance(samples[0], dict) else {}
        if not sample_payload:
            continue
        clade = str(sample_payload.get("predicted_group") or sample_payload.get("predicted_clade") or "").strip()
        lineage = str(sample_payload.get("predicted_clade") or "").strip()
        species = "HIV-1" if clade else ""
        return {"species": species, "clade": clade, "lineage": lineage}
    return {}


def _read_multi_sample_depth_coverage(report_dir: Path, sample_name: str) -> dict[str, str]:
    stats = _read_multi_sample_region_depth_stats(report_dir, sample_name) or _read_multi_sample_per_base_stats(report_dir, sample_name)
    total_bases = _safe_int(stats.get("total_bases"))
    if total_bases:
        return {
            "coverage_1x": _display_percent(stats.get("covered_1x") / total_bases),
            "coverage_10x": _display_percent(stats.get("covered_10x") / total_bases),
            "coverage_100x": _display_percent(stats.get("covered_100x") / total_bases),
        }
    return _read_multi_sample_depth_coverage_from_dist(report_dir, sample_name)


def _read_multi_sample_depth_coverage_from_dist(report_dir: Path, sample_name: str) -> dict[str, str]:
    candidates = [
        report_dir / f"{sample_name}.mosdepth.global.dist.txt",
        report_dir / f"{sample_name}.mosdepth.region.dist.txt",
        report_dir / "ref_map.mosdepth.global.dist.txt",
        report_dir / "ref_map.mosdepth.region.dist.txt",
    ]
    candidates.extend(sorted(report_dir.glob("*.mosdepth.global.dist.txt")))
    candidates.extend(sorted(report_dir.glob("*.mosdepth.region.dist.txt")))
    seen: set[Path] = set()
    for path in candidates:
        resolved = path.resolve()
        if resolved in seen or not path.is_file():
            continue
        seen.add(resolved)
        values: dict[int, float] = {}
        try:
            with path.open("r", encoding="utf-8", errors="ignore") as handle:
                for line in handle:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 3:
                        continue
                    if parts[0] != "total":
                        continue
                    depth = _safe_int(parts[1])
                    fraction = _safe_float(parts[2])
                    if depth in {1, 10, 100} and fraction is not None:
                        values[depth] = fraction
        except OSError:
            continue
        if values:
            return {
                "coverage_1x": _display_percent(values.get(1)) if 1 in values else "",
                "coverage_10x": _display_percent(values.get(10)) if 10 in values else "",
                "coverage_100x": _display_percent(values.get(100)) if 100 in values else "",
            }
    return {}


def _read_multi_sample_mean_depth(report_dir: Path, sample_name: str) -> str:
    candidates = [
        report_dir / f"{sample_name}.mosdepth.summary.txt",
        report_dir / "ref_map.mosdepth.summary.txt",
    ]
    candidates.extend(sorted(report_dir.glob("*.mosdepth.summary.txt")))
    seen: set[Path] = set()
    for path in candidates:
        resolved = path.resolve()
        if resolved in seen or not path.is_file():
            continue
        seen.add(resolved)
        table = _read_tsv_rows(path)
        columns = table.get("columns") or []
        rows = table.get("rows") or []
        if not columns or not rows:
            continue
        chrom_index = _find_column_index(columns, ["chrom", "chromosome", "region", "参考序列"])
        mean_index = _find_column_index(columns, ["mean", "平均深度", "depth", "mean_depth"])
        if mean_index < 0:
            continue
        selected_row = None
        for preferred_chrom in ("total_region", "total"):
            selected_row = next((row for row in rows if _cell_at(row, chrom_index) == preferred_chrom), None)
            if selected_row:
                break
        if selected_row is None:
            selected_row = rows[0]
        mean_value = _safe_float(_cell_at(selected_row, mean_index))
        if mean_value is not None:
            return f"{mean_value:.2f}"
    return _read_multi_sample_mean_depth_from_region_or_per_base(report_dir, sample_name)


def _read_multi_sample_mean_depth_from_region_or_per_base(report_dir: Path, sample_name: str) -> str:
    stats = _read_multi_sample_region_depth_stats(report_dir, sample_name) or _read_multi_sample_per_base_stats(report_dir, sample_name)
    mean_depth = _safe_float(stats.get("mean_depth"))
    return f"{mean_depth:.2f}" if mean_depth is not None else ""


def _read_multi_sample_region_depth_stats(report_dir: Path, sample_name: str) -> dict[str, float]:
    candidates = [
        report_dir / "ref.regions.bed",
        report_dir / "ref.regions.bed.gz",
        report_dir / f"{sample_name}.regions.bed.gz",
        report_dir / f"{sample_name}.regions.bed",
        report_dir / f"{sample_name}_ngs.regions.bed.gz",
        report_dir / f"{sample_name}_ngs.regions.bed",
        report_dir / "ref_map.regions.bed.gz",
        report_dir / "ref_map.regions.bed",
    ]
    candidates.extend(sorted(report_dir.glob("*.regions.bed.gz")))
    candidates.extend(sorted(report_dir.glob("*.regions.bed")))
    return _read_multi_sample_depth_stats_from_bed_candidates(candidates)


def _read_multi_sample_per_base_stats(report_dir: Path, sample_name: str) -> dict[str, float]:
    candidates = [
        report_dir / f"{sample_name}_ngs.per-base.bed.gz",
        report_dir / f"{sample_name}_ngs.per-base.bed",
        report_dir / f"{sample_name}.per-base.bed.gz",
        report_dir / f"{sample_name}.per-base.bed",
        report_dir / "ref_map.per-base.bed.gz",
        report_dir / "ref_map.per-base.bed",
    ]
    candidates.extend(sorted(report_dir.glob("*.per-base.bed.gz")))
    candidates.extend(sorted(report_dir.glob("*.per-base.bed")))
    return _read_multi_sample_depth_stats_from_bed_candidates(candidates)


def _read_multi_sample_depth_stats_from_bed_candidates(candidates: list[Path]) -> dict[str, float]:
    seen: set[Path] = set()
    for path in candidates:
        resolved = path.resolve()
        if resolved in seen or not path.is_file():
            continue
        seen.add(resolved)
        opener = gzip.open if path.suffix == ".gz" else open
        depth_total = 0.0
        span_total = 0
        covered_1x = 0
        covered_10x = 0
        covered_100x = 0
        try:
            with opener(path, "rt", encoding="utf-8", errors="ignore") as handle:
                for line in handle:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 4:
                        continue
                    depth = _safe_float(parts[3])
                    if depth is None:
                        continue
                    start = _safe_int(parts[1])
                    end = _safe_int(parts[2])
                    span = max(1, (end - start) if start is not None and end is not None else 1)
                    depth_total += depth * span
                    span_total += span
                    if depth >= 1:
                        covered_1x += span
                    if depth >= 10:
                        covered_10x += span
                    if depth >= 100:
                        covered_100x += span
        except OSError:
            continue
        if span_total:
            return {
                "total_bases": float(span_total),
                "mean_depth": depth_total / span_total,
                "covered_1x": float(covered_1x),
                "covered_10x": float(covered_10x),
                "covered_100x": float(covered_100x),
            }
    return {}


def _read_multi_sample_serotype_table(report_dir: Path, sample_name: str) -> dict:
    exact_path = report_dir / f"{sample_name}_serotype_result.tsv"
    exact_table = _read_tsv_rows(exact_path)
    exact_row = _first_row_mapping(exact_table)
    if _has_meaningful_serotype_value(_format_multi_sample_bacterial_serotype(exact_row)):
        return exact_table

    candidates = sorted(
        [
            path
            for path in report_dir.glob("*_serotype_result.tsv")
            if path.is_file() and path != exact_path
        ],
        key=lambda item: item.name.lower(),
    )
    fallback_table = exact_table if exact_row else {}
    for path in candidates:
        table = _read_tsv_rows(path)
        row = _first_row_mapping(table)
        if _has_meaningful_serotype_value(_format_multi_sample_bacterial_serotype(row)):
            return table
        if not _first_row_mapping(fallback_table) and row:
            fallback_table = table
    return fallback_table


def _has_meaningful_serotype_value(value: object) -> bool:
    text = str(value or "").strip()
    return bool(text) and text not in {"-", "--", "-:-", "-|-", "nan", "None"}


def _format_multi_sample_bacterial_serotype(serotype_row: dict[str, str]) -> str:
    ko_serotype = _first_mapping_value(serotype_row, [
        "KO血清型",
        "K/O血清型",
        "K-O血清型",
    ])
    if _has_meaningful_serotype_value(ko_serotype):
        return ko_serotype

    k_locus = _first_mapping_value(serotype_row, [
        "K_locus",
        "klebsiella_pneumo_complex__kaptive__K_locus",
    ])
    o_locus = _first_mapping_value(serotype_row, [
        "O_locus",
        "klebsiella_pneumo_complex__kaptive__O_locus",
    ])
    if _has_meaningful_serotype_value(k_locus) or _has_meaningful_serotype_value(o_locus):
        return f"{k_locus or '-'}|{o_locus or '-'}"

    annotation = _first_mapping_value(serotype_row, [
        "血清型注释信息(simple)",
        "血清型注释信息(details)",
        "亚型全称",
        "菌种",
    ])
    fallback = _first_mapping_value(serotype_row, [
        "血清型结果",
        "分型结果",
        "serotype",
        "血清型",
        "O抗原",
        "H抗原",
    ])
    subtype = _extract_salmonella_serovar_name(annotation or fallback)
    antigen_formula = _first_mapping_value(serotype_row, [
        "血清式",
        "抗原组成",
        "抗原式",
        "antigenic_formula",
        "antigen_formula",
    ]) or _extract_salmonella_antigen_formula(annotation)
    if subtype and antigen_formula:
        return f"{subtype} / {antigen_formula}"
    if subtype:
        return subtype
    if antigen_formula:
        return antigen_formula
    return fallback or annotation

def _extract_salmonella_serovar_name(value: object) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    prefix = re.split(r"\s*[\(（]", text, maxsplit=1)[0].strip(" ，,;；")
    for candidate in [prefix, text]:
        normalized_candidate = _normalize_serotype_lookup_text(candidate)
        for canonical, aliases in SALMONELLA_SEROVAR_ALIAS_MAP.items():
            normalized_aliases = {_normalize_serotype_lookup_text(item) for item in [canonical, *aliases]}
            if normalized_candidate in normalized_aliases:
                return canonical
    if "沙门" in prefix:
        return prefix.replace("沙门氏菌", "沙门菌")
    return prefix

def _extract_salmonella_antigen_formula(value: object) -> str:
    text = str(value or "").strip().replace("（", "(").replace("）", ")")
    if not text:
        return ""
    match = re.search(r",\s*\(([^()]*)\)\s*\)?$", text)
    if match:
        return match.group(1).strip()
    return ""

def _build_multi_sample_overview_table(sample_rows: list[dict], *, is_virus: bool) -> dict:
    if is_virus:
        columns = [
            {"key": "sample", "label": "样本"},
            {"key": "ready_label", "label": "状态"},
            {"key": "total_bases_label", "label": "测序量"},
            {"key": "q_label", "label": "Q20 / Q30"},
            {"key": "species_name", "label": "病毒/物种"},
            {"key": "typing", "label": "分型/谱系"},
            {"key": "nextclade_typing", "label": "Nextclade分型"},
            {"key": "coverage_1x", "label": "1x覆盖度"},
            {"key": "coverage_10x", "label": "10x覆盖度"},
            {"key": "coverage_100x", "label": "100x覆盖度"},
            {"key": "mean_depth", "label": "平均深度"},
            {"key": "assembly_label", "label": "组装"},
            {"key": "qc_status", "label": "分型QC"},
            {"key": "note", "label": "关注说明"},
        ]
    else:
        columns = [
            {"key": "sample", "label": "样本"},
            {"key": "ready_label", "label": "状态"},
            {"key": "total_bases_label", "label": "测序量"},
            {"key": "q_label", "label": "Q20 / Q30"},
            {"key": "species_name", "label": "物种鉴定"},
            {"key": "typing", "label": "MLST/ST"},
            {"key": "serotype", "label": "血清型注释"},
            {"key": "checkm_label", "label": "完整性/污染率"},
            {"key": "assembly_label", "label": "组装"},
            {"key": "rv_label", "label": "耐药/毒力"},
            {"key": "note", "label": "关注说明"},
        ]
    rows = []
    for item in sample_rows:
        q20 = _display_percent(item.get("q20_rate"))
        q30 = _display_percent(item.get("q30_rate"))
        completeness = _display_percent_points(item.get("completeness"))
        contamination = _display_percent_points(item.get("contamination"))
        rows.append({
            **item,
            "ready_label": "已生成" if item.get("ready") else "缺结果",
            "total_bases_label": _human_bp(item.get("total_bases")),
            "q_label": f"{q20} / {q30}",
            "checkm_label": f"{completeness} / {contamination}",
            "assembly_label": f"{_human_count(item.get('contig_count'))} contig / {_human_bp(item.get('assembly_length'))}",
            "rv_label": f"{_human_count(item.get('resistance_count'))} / {_human_count(item.get('virulence_count'))}",
        })
    return {"columns": columns, "rows": rows}

def _read_assembly_profile(path: Path) -> dict:
    if not path.is_file():
        return {"contig_count": None, "plasmid_count": None, "total_count": None, "total_length": None}
    contig_count = 0
    plasmid_count = 0
    total_count = 0
    total_length = 0
    try:
        with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
            reader = csv.reader(handle, delimiter="	")
            next(reader, None)
            for row in reader:
                if len(row) < 5:
                    continue
                total_count += 1
                length_value = _safe_int(row[1] if len(row) > 1 else None)
                if length_value is not None:
                    total_length += length_value
                genome_type = str(row[4]).lower()
                if "plasmid" in genome_type:
                    plasmid_count += 1
                else:
                    contig_count += 1
    except OSError:
        return {"contig_count": None, "plasmid_count": None, "total_count": None, "total_length": None}
    return {
        "contig_count": contig_count,
        "plasmid_count": plasmid_count,
        "total_count": total_count,
        "total_length": total_length or None,
    }

def _read_fasta_assembly_summary(path: Path) -> dict:
    if not path.is_file():
        return {"columns": [], "rows": []}
    contig_count = 0
    total_length = 0
    max_length = 0
    lengths: list[int] = []
    current_length = 0
    try:
        with path.open("r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if line.startswith(">"):
                    if current_length:
                        lengths.append(current_length)
                        total_length += current_length
                        max_length = max(max_length, current_length)
                        current_length = 0
                    contig_count += 1
                else:
                    current_length += len(line.strip())
            if current_length:
                lengths.append(current_length)
                total_length += current_length
                max_length = max(max_length, current_length)
    except OSError:
        return {"columns": [], "rows": []}
    if not contig_count:
        return {"columns": [], "rows": []}
    lengths.sort(reverse=True)
    half_total = total_length / 2 if total_length else 0
    cumulative = 0
    n50 = 0
    for length in lengths:
        cumulative += length
        if cumulative >= half_total:
            n50 = length
            break
    avg_length = round(total_length / contig_count, 2) if contig_count else 0
    return {
        "columns": ["结果文件", "Contig数量", "总长度(bp)", "最大Contig(bp)", "N50(bp)", "平均长度(bp)"],
        "rows": [[path.name, str(contig_count), str(total_length), str(max_length), str(n50), f"{avg_length:.2f}"]],
    }
