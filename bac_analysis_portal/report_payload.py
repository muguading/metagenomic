from __future__ import annotations

import csv
import json
from pathlib import Path

from .assembly_reports import (
    _build_virus_fallback_assembly_summary,
    _discover_cgview_assets,
    _fallback_assembly_profile,
    _merge_assembly_coverage_summary,
    _normalize_bandavirus_assembly_coverage,
    _read_contig_depth_relationship,
    _read_coverage_profile,
)
from .binning_reports import (
    _build_binning_quality_section,
    _build_binning_taxonomy_section,
    _build_meta_viral_assembly_section,
    _read_checkm2_quality,
)
from .community_reports import _build_community_report_payload, _is_community_report_task
from .export_checklist import build_export_checklist
from .knowledge_interpretation import (
    _build_kb_taxonomy_index,
    _build_knowledge_interpretation,
    _build_viral_serotype_knowledge_summary,
    _enrich_priority_serotype_rows,
)
from .meta_risk_reports import (
    _build_category_gene_relationship,
    _build_resistance_virulence_summary,
    _merge_annotation_with_summary,
    _read_gene_length_distribution,
    _read_meta_resistance_virulence,
    _read_mge_monitoring,
)
from .parse_utils import _safe_float, _safe_int
from .pathosource_reports import _build_pathosource_report_payload, _is_pathosource_report_task
from .public_health_support import _build_public_health_support
from .report_artifacts import (
    _read_gff_genome_features,
    _read_ncov_genome_features,
    _resolve_report_artifact_path,
    _simplify_virus_coverage_features,
)
from .report_cache import (
    _read_report_cache,
    _refresh_cached_report_payload,
    _report_cache_fingerprint,
    _write_report_cache,
)
from .report_formatters import (
    _coerce_percent_value,
    _display_percent,
    _display_percent_points,
    _human_bp,
    _human_count,
)
from .report_sources import (
    _build_multi_sample_queue_summary,
    _read_assembly_profile,
    _read_checkm_metrics,
    _read_fasta_assembly_summary,
    _read_summary_metrics,
    _resolve_report_sample_display_name,
    _resolve_report_sample_name,
    _resolve_report_source,
)
from .report_modeling import build_modeling_risk_section
from .runtime_paths import _resolve_runtime_database_root
from .serotype_reports import _build_serotype_section
from .table_io import _read_tsv_rows
from .task_manager import ValidationError
from .taxonomy_reports import (
    _build_taxonomy_abundance,
    _build_taxonomy_interpretation,
    _build_taxonomy_rarefaction,
    _build_taxonomy_risk_summary,
    _extract_dominant_fungus_taxonomy,
    _extract_dominant_species,
    _extract_dominant_virus_taxonomy,
    _read_assembly_taxonomy,
    _read_taxonomy_list,
)
from .typing_reports import (
    _build_neisseria_amr_section,
    _build_tb_amr_section,
    _build_tb_serotype_section,
    _read_mlst_result,
)
from .virus_report_templates import _build_virus_report_template_summary
from .workflow_closure import build_workflow_closure

def _build_report_payload(task: dict, selected_sample: str = "") -> dict:
    project_root = Path(__file__).resolve().parent.parent
    database_root = _resolve_runtime_database_root()
    params = task.get("params", {})
    report_source = _resolve_report_source(task, selected_sample)
    if not report_source.get("available"):
        raise ValidationError(str(report_source.get("reason") or "服务器结果目录尚未就绪。"))
    report_dir = report_source["report_dir"]
    sample_name = report_source.get("selected_sample") or _resolve_report_sample_name(task, report_dir)
    sample_display_name = _resolve_report_sample_display_name(task, sample_name)
    if _is_community_report_task(task):
        cache_fingerprint = _report_cache_fingerprint(task, report_dir, sample_name)
        cached_payload = _read_report_cache(report_dir, sample_name, cache_fingerprint)
        if cached_payload:
            return _refresh_cached_report_payload(
                cached_payload,
                task=task,
                report_dir=report_dir,
                report_source=report_source,
                sample_name=sample_name,
                sample_display_name=sample_display_name,
            )
        payload = _build_community_report_payload(
            task=task,
            report_dir=report_dir,
            report_source=report_source,
        )
        _write_report_cache(report_dir, sample_name, cache_fingerprint, payload)
        return payload
    if _is_pathosource_report_task(task):
        cache_fingerprint = _report_cache_fingerprint(task, report_dir, sample_name)
        cached_payload = _read_report_cache(report_dir, sample_name, cache_fingerprint)
        if cached_payload:
            return _refresh_cached_report_payload(
                cached_payload,
                task=task,
                report_dir=report_dir,
                report_source=report_source,
                sample_name=sample_name,
                sample_display_name=sample_display_name,
            )
        payload = _build_pathosource_report_payload(
            task=task,
            report_dir=report_dir,
            report_source=report_source,
            sample_name=sample_name,
            sample_display_name=sample_display_name,
        )
        _write_report_cache(report_dir, sample_name, cache_fingerprint, payload)
        return payload
    cache_fingerprint = _report_cache_fingerprint(task, report_dir, sample_name)
    cached_payload = _read_report_cache(report_dir, sample_name, cache_fingerprint)
    if cached_payload:
        return _refresh_cached_report_payload(
            cached_payload,
            task=task,
            report_dir=report_dir,
            report_source=report_source,
            sample_name=sample_name,
            sample_display_name=sample_display_name,
        )
    fastp_path = _resolve_report_artifact_path(
        report_dir,
        [f"{sample_name}.fastp2.json"] if sample_name else [],
        ["*.fastp2.json"],
    )
    coverage_path = _resolve_report_artifact_path(
        report_dir,
        ([ "ref.regions.bed", f"{sample_name}_ngs.per-base.bed.gz"] if sample_name else ["ref.regions.bed"]),
        ["ref.regions.bed", "*.regions.bed", "*_ngs.per-base.bed.gz", "*.regions.bed.gz", "*.per-base.bed.gz", "*.per-base.bed"],
    )
    summary_info = _read_summary_metrics(report_dir / "summary.tsv")
    checkm_info = _read_checkm_metrics(report_dir / f"{sample_name}.checkm.tsv") if sample_name else {}
    assembly_profile = _read_assembly_profile(report_dir / "Assem_info.tsv")
    contig_depth_relationship = _read_contig_depth_relationship(report_dir / "Assem_info.tsv")
    fastp_info = _read_fastp_metrics(fastp_path) if fastp_path else {}
    assembly_summary = _read_tsv_rows(report_dir / f"{sample_name}.assemble.result.tsv") if sample_name else {"columns": [], "rows": []}
    assembly_coverage = _read_coverage_profile(coverage_path) if coverage_path else {}
    is_virus_analysis = str(params.get("analysis_target", "bacteria")).strip() == "virus"
    species_hint = str(params.get("species") or "").strip().lower()
    hpiv_coverage_path = report_dir / "hpiv_coverage" / "hpiv.coverage.regions.bed"
    if ("parainfluenza" in species_hint or "副流感" in species_hint or species_hint == "hpiv") and hpiv_coverage_path.is_file():
        coverage_path = hpiv_coverage_path
        assembly_coverage = _read_coverage_profile(coverage_path)
    if is_virus_analysis and not assembly_summary.get("rows"):
        assembly_summary = _build_virus_fallback_assembly_summary(report_dir, sample_name, assembly_coverage)
    if is_virus_analysis and assembly_summary.get("rows") and assembly_coverage.get("status") == "ready":
        assembly_summary = _merge_assembly_coverage_summary(assembly_summary, assembly_coverage)
    nextclade_report_path = report_dir / "nextclade_output" / "nextclade.tsv"
    if assembly_coverage.get("status") == "ready":
        if nextclade_report_path.is_file() and ("sars-cov-2" in species_hint or "cov" in species_hint or "新冠" in species_hint):
            assembly_coverage["view_mode"] = "ncov_annotated"
            assembly_coverage["reference_name"] = "NC_045512.2"
            assembly_coverage["annotation_source"] = "database/virus/ncov/genomic.gff"
            assembly_coverage["annotation_label"] = "Wuhan-Hu-1 参考基因组 GFF"
            assembly_coverage["genome_features"] = _read_ncov_genome_features()
        elif (
            "respiratory syncytial virus" in species_hint
            or "rsv" in species_hint
            or "合胞病毒" in species_hint
            or "human metapneumovirus" in species_hint
            or "hmpv" in species_hint
            or "偏肺病毒" in species_hint
            or "dengue virus" in species_hint
            or species_hint == "denv"
            or "登革热" in species_hint
            or "zika virus" in species_hint
            or species_hint == "zikav"
            or "寨卡" in species_hint
            or "chikungunya" in species_hint
            or species_hint == "chikv"
            or "基孔肯雅" in species_hint
            or "parainfluenza" in species_hint
            or species_hint == "hpiv"
            or "副流感" in species_hint
            or "adenovirus" in species_hint
            or species_hint == "hadv"
            or "腺病毒" in species_hint
            or "norovirus" in species_hint
            or "norwalk" in species_hint
            or "诺如" in species_hint
            or "rhinovirus" in species_hint
            or "hrv" in species_hint
            or "鼻病毒" in species_hint
            or "human coronavirus" in species_hint
            or "seasonal coronavirus" in species_hint
            or "seasonal_hcov" in species_hint
            or "季节性冠状病毒" in species_hint
            or "229e" in species_hint
            or "nl63" in species_hint
            or "oc43" in species_hint
            or "hku1" in species_hint
            or "orthohantavirus" in species_hint
            or "hantavirus" in species_hint
            or "汉坦" in species_hint
            or "汉他" in species_hint
            or "ebola virus" in species_hint
            or "ebolavirus" in species_hint
            or "orthoebolavirus" in species_hint
            or "ebov" in species_hint
            or "埃博拉" in species_hint
        ):
            virus_gff_candidates = [report_dir / "ref" / "genes.gff"]
            for selection_dir in sorted(report_dir.glob("*_reference_selection")):
                virus_gff_candidates.extend(
                    [
                        selection_dir / "snpeff_reference.gff3",
                        selection_dir / "snpeff_vadr_ref" / "ref" / "genes.gff",
                    ]
                )
                selection_tsv = selection_dir / "selection.tsv"
                if not selection_tsv.is_file():
                    continue
                try:
                    with selection_tsv.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
                        selection_row = next(csv.DictReader(handle, delimiter="\t"), None) or {}
                except OSError:
                    selection_row = {}
                gff_text = str(selection_row.get("gff_path") or "").strip()
                if gff_text and gff_text != "nogtf":
                    virus_gff_candidates.append(Path(gff_text).expanduser())
            virus_gff_path = next((path for path in virus_gff_candidates if path.is_file() and path.stat().st_size > 0), virus_gff_candidates[0])
            virus_features = _read_gff_genome_features(virus_gff_path)
            virus_features = _simplify_virus_coverage_features(virus_features, species_hint)
            if virus_features:
                assembly_coverage["view_mode"] = "ncov_annotated"
                assembly_coverage["reference_name"] = str(params.get("species") or "Virus")
                try:
                    annotation_source = str(virus_gff_path.relative_to(report_dir))
                except ValueError:
                    annotation_source = virus_gff_path.name
                assembly_coverage["annotation_source"] = annotation_source
                assembly_coverage["annotation_label"] = "病毒参考基因组注释 GFF"
                assembly_coverage["genome_features"] = virus_features
    contig_annotation = _read_tsv_rows(report_dir / "flye_output" / "assembly_info.txt")
    assembly_profile = _fallback_assembly_profile(assembly_profile, contig_annotation, assembly_summary)
    checkm_quality = _read_checkm2_quality(report_dir / "checkm2_out" / "quality_report.tsv")
    gene_annotation_summary = _read_tsv_rows(report_dir / f"{sample_name}.genefun_summary.tsv") if sample_name else {"columns": [], "rows": []}
    gene_length_distribution = _read_gene_length_distribution(report_dir / f"{sample_name}_gene_raw_sum.tsv") if sample_name else {"status": "empty", "points": []}
    is_meta_method = str(params.get("method", "")).strip() == "meta"
    if is_meta_method:
        meta_rv = _read_meta_resistance_virulence(report_dir / "meta_plas_vf_card.tsv")
        rv_summary = meta_rv.get("summary", {"columns": [], "rows": []})
        virulence_elements = meta_rv.get("virulence_elements", {"columns": [], "rows": []})
        resistance_elements = meta_rv.get("resistance_elements", {"columns": [], "rows": []})
        virulence_relationship = _build_category_gene_relationship(
            virulence_elements,
            left_key="VF分类",
            right_key="基因名称",
            label="基因组/质粒与毒力基因关系图",
        )
        resistance_relationship = _build_category_gene_relationship(
            resistance_elements,
            left_key="耐药药物",
            right_key="基因名称",
            label="基因组/质粒与耐药基因关系图",
            split_delimiters=[";"],
        )
    else:
        rv_summary = _read_tsv_rows(report_dir / "Assem_info1.tsv")
        virulence_elements = _merge_annotation_with_summary(
            report_dir / "Assem_abricate_VFDB.tsv",
            report_dir / "VFDB_summary.tsv",
            mode="virulence",
        )
        virulence_relationship = _build_category_gene_relationship(
            virulence_elements,
            left_key="VF分类",
            right_key="基因名称",
            label="VF 分类与毒力基因关系图",
        )
        resistance_elements = _merge_annotation_with_summary(
            report_dir / "Assem_abricate_CARD.tsv",
            report_dir / "CARD_summary.tsv",
            mode="resistance",
        )
        resistance_relationship = _build_category_gene_relationship(
            resistance_elements,
            left_key="耐药药物",
            right_key="基因名称",
            label="耐药药物与耐药基因关系图",
            split_delimiters=[";"],
        )
    rv_overview = _build_resistance_virulence_summary(rv_summary, virulence_elements, resistance_elements)
    mge_monitoring = _read_mge_monitoring(report_dir, sample_name)
    taxonomy_index = _build_kb_taxonomy_index(str(database_root))
    species_taxonomy = _read_taxonomy_list(report_dir / f"{sample_name}_2.list.txt", terminal_column="种", taxonomy_index=taxonomy_index) if sample_name else {"rows": [], "rank_options": []}
    subspecies_taxonomy = _read_taxonomy_list(report_dir / f"{sample_name}_2.list2.txt", terminal_column="亚种", taxonomy_index=taxonomy_index) if sample_name else {"rows": [], "rank_options": []}
    taxonomy_abundance = _build_taxonomy_abundance(species_taxonomy, subspecies_taxonomy)
    taxonomy_risk_summary = _build_taxonomy_risk_summary(species_taxonomy, subspecies_taxonomy)
    taxonomy_interpretation = _build_taxonomy_interpretation(species_taxonomy, checkm_info)
    assembly_taxonomy = _read_assembly_taxonomy(
        report_dir / "Assem_info1.tsv",
        report_dir / f"{sample_name}_assem.kraken2.txt",
        taxonomy_index=taxonomy_index,
    ) if sample_name else {"columns": [], "rows": []}
    binning_quality = _build_binning_quality_section(report_dir / "bin_checkm2out" / "quality_report.tsv")
    binning_taxonomy = _build_binning_taxonomy_section(report_dir / "gtdbtk_out" / "gtdbtk.bac120.summary.tsv")
    viral_assembly = _build_meta_viral_assembly_section(report_dir / "viral_assembly") if is_meta_method else {"status": "empty", "summary": {}, "table": {"columns": [], "rows": []}}
    mlst_result = _read_mlst_result(report_dir / f"{sample_name}.mlst_Stat.txt", str(project_root)) if sample_name else {"columns": [], "rows": [], "gene_show_map": {}, "default_gene": ""}
    neisseria_amr_result = _build_neisseria_amr_section(report_dir, sample_name, mlst_result, checkm_info, project_root) if sample_name else {"status": "empty", "columns": [], "rows": []}
    tb_amr_result = _build_tb_amr_section(report_dir, sample_name, checkm_info) if sample_name else {"status": "empty", "columns": [], "rows": []}
    serotype_result = _build_serotype_section(report_dir, sample_name, checkm_info) if sample_name else {"status": "empty", "mode": "generic", "columns": [], "rows": []}
    tb_serotype_result = _build_tb_serotype_section(report_dir, sample_name, checkm_info) if sample_name else None
    if isinstance(tb_serotype_result, dict) and str(tb_serotype_result.get("status") or "").strip() != "":
        serotype_result = tb_serotype_result
    if str(serotype_result.get("mode") or "").strip() == "bandavirus_typing":
        assembly_coverage = _normalize_bandavirus_assembly_coverage(assembly_coverage, serotype_result)
    priority_serotype = _read_tsv_rows(report_dir / f"{sample_name}.pathonet_result.tsv") if sample_name else {"columns": [], "rows": []}
    if sample_name:
        priority_serotype = _enrich_priority_serotype_rows(str(project_root), priority_serotype, serotype_result, checkm_info)
    if is_meta_method:
        assembly_summary = _read_fasta_assembly_summary(report_dir / "tmp_combine.fa")
        assembly_coverage = {"status": "empty", "points": []}
        contig_annotation = _read_meta_annotation_subset(report_dir / "meta_plas_vf_card.tsv")
        contig_depth_relationship = {"status": "empty", "points": [], "length_depth_scatter": {"status": "empty", "points": []}}
        assembly_profile = _fallback_assembly_profile(assembly_profile, contig_annotation, assembly_summary)
    cgview_assets = _discover_cgview_assets(report_dir, sample_name) if not is_meta_method else {"status": "empty", "summary": {"map_count": 0}, "maps": []}
    dominant_species = _extract_dominant_species(species_taxonomy)
    dominant_virus_taxonomy = _extract_dominant_virus_taxonomy(species_taxonomy, subspecies_taxonomy)
    dominant_fungus_taxonomy = _extract_dominant_fungus_taxonomy(species_taxonomy, subspecies_taxonomy)
    public_health_support = _build_public_health_support(
        project_root,
        checkm_info.get("species_name") or dominant_species.get("species") or "",
        resistance_elements,
    )
    knowledge_interpretation = _build_knowledge_interpretation(
        project_root,
        species_taxonomy,
        subspecies_taxonomy,
        resistance_elements,
        virulence_elements,
        mge_monitoring,
        mlst_result,
        serotype_result,
        public_health_support,
        is_meta_method,
    )
    total_reads = _safe_int(summary_info.get("sum_reads"))
    if total_reads is None:
        total_reads = _safe_int(((fastp_info.get("summary") or {}).get("before_filtering") or {}).get("total_reads"))
    fastp_before_summary = (fastp_info.get("summary") or {}).get("before_filtering") or {}
    resolved_total_bases = summary_info.get("sum_len")
    if resolved_total_bases is None:
        resolved_total_bases = fastp_before_summary.get("total_bases")
    resolved_q20_rate = summary_info.get("q20_rate")
    if resolved_q20_rate is None:
        resolved_q20_rate = _coerce_percent_value(fastp_before_summary.get("q20_rate"))
    resolved_q30_rate = summary_info.get("q30_rate")
    if resolved_q30_rate is None:
        resolved_q30_rate = _coerce_percent_value(fastp_before_summary.get("q30_rate"))
    if isinstance(serotype_result, dict):
        mode = str(serotype_result.get("mode") or "").strip()
        if mode == "tb_profiler":
            serotype_result["knowledge_summary"] = serotype_result.get("knowledge_summary") or {"headline": "", "items": []}
        else:
            serotype_result["knowledge_summary"] = _build_viral_serotype_knowledge_summary(
                project_root,
                params.get("species", ""),
                serotype_result,
            )
            serotype_result["report_template"] = _build_virus_report_template_summary(project_root, serotype_result)
    influenza_overview_metrics: list[dict] = []
    monkeypox_overview_metrics: list[dict] = []
    rsv_overview_metrics: list[dict] = []
    hmpv_overview_metrics: list[dict] = []
    denv_overview_metrics: list[dict] = []
    zikav_overview_metrics: list[dict] = []
    chikv_overview_metrics: list[dict] = []
    ebola_overview_metrics: list[dict] = []
    hpiv_overview_metrics: list[dict] = []
    hiv_overview_metrics: list[dict] = []
    hadv_overview_metrics: list[dict] = []
    norovirus_overview_metrics: list[dict] = []
    enterovirus_overview_metrics: list[dict] = []
    hepatovirus_overview_metrics: list[dict] = []
    bandavirus_overview_metrics: list[dict] = []
    orthohantavirus_overview_metrics: list[dict] = []
    astroviridae_overview_metrics: list[dict] = []
    rhinovirus_overview_metrics: list[dict] = []
    seasonal_hcov_overview_metrics: list[dict] = []
    rotavirus_overview_metrics: list[dict] = []
    if str(serotype_result.get("mode") or "").strip() == "influenza_typing":
        segment_manifest = serotype_result.get("segment_manifest") if isinstance(serotype_result.get("segment_manifest"), dict) else {}
        segment_columns = segment_manifest.get("columns") if isinstance(segment_manifest.get("columns"), list) else []
        segment_rows = segment_manifest.get("rows") if isinstance(segment_manifest.get("rows"), list) else []
        coverage_segments = assembly_coverage.get("segments") if isinstance(assembly_coverage.get("segments"), list) else []
        coverage_by_name: dict[str, dict] = {}
        coverage_by_group: dict[str, dict] = {}
        for item in coverage_segments:
            if not isinstance(item, dict):
                continue
            coverage_name = str(item.get("name") or item.get("label") or "").strip()
            if not coverage_name:
                continue
            coverage_by_name[coverage_name] = item
            parts = coverage_name.split("_")
            if len(parts) >= 2 and parts[1]:
                coverage_by_group[parts[1]] = item
        segment_name_index = -1
        reference_id_index = -1
        segment_subtype_index = -1
        for index, value in enumerate(segment_columns):
            text = str(value or "").strip()
            if segment_name_index < 0 and text in {"segment_group", "片段"}:
                segment_name_index = index
            if reference_id_index < 0 and text in {"reference_id", "参考ID"}:
                reference_id_index = index
            if segment_subtype_index < 0 and (text in {"subtype", "亚型"} or text.lower() == "subtype"):
                segment_subtype_index = index
        segment_labels: list[str] = []
        for row in segment_rows:
            if not isinstance(row, list):
                continue
            segment_name = str(row[segment_name_index] if segment_name_index >= 0 and segment_name_index < len(row) else "").strip()
            reference_id = str(row[reference_id_index] if reference_id_index >= 0 and reference_id_index < len(row) else "").strip()
            if not segment_name:
                continue
            segment_subtype = str(row[segment_subtype_index] if segment_subtype_index >= 0 and segment_subtype_index < len(row) else "").strip()
            label = f"{segment_name}({segment_subtype})" if segment_subtype and segment_subtype != "-" else segment_name
            coverage_item = coverage_by_name.get(reference_id) or coverage_by_group.get(segment_name)
            mean_depth = _safe_float(coverage_item.get("mean_depth")) if isinstance(coverage_item, dict) else None
            if mean_depth is not None:
                label = f"{label}: {mean_depth:.2f}x"
            if label not in segment_labels:
                segment_labels.append(label)
        influenza_type = str(serotype_result.get("influenza_type") or "--").strip() or "--"
        ha_subtype = str(serotype_result.get("ha_subtype") or "--").strip() or "--"
        na_subtype = str(serotype_result.get("na_subtype") or "--").strip() or "--"
        influenza_overview_metrics = [
            {
                "key": "influenza_segments",
                "label": "组装情况",
                "type": "influenza_segments",
                "segment_count": len(segment_labels),
                "segments": segment_labels,
            },
            {
                "key": "influenza_species_estimation",
                "label": "物种预估",
                "type": "influenza_species_estimation",
                "influenza_type": influenza_type,
                "ha_subtype": ha_subtype,
                "na_subtype": na_subtype,
            },
        ]
    if str(serotype_result.get("mode") or "").strip() == "monkeypox_nextclade":
        monkeypox_overview_metrics = [
            {
                "key": "monkeypox_assembly_coverage",
                "label": "组装情况",
                "type": "paired",
                "items": [
                    {"label": "1x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_fraction"))},
                    {"label": "10x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_10x_fraction"))},
                    {"label": "100x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_100x_fraction"))},
                ],
            },
            {
                "key": "monkeypox_species_estimation",
                "label": "物种预估",
                "type": "paired",
                "items": [
                    {"label": "Nextclade Clade", "display": str(serotype_result.get("predicted_clade") or "--")},
                    {"label": "Lineage", "display": str(serotype_result.get("predicted_lineage") or "--")},
                ],
            },
        ]
    if str(serotype_result.get("mode") or "").strip() == "rsv_nextclade":
        rsv_overview_metrics = [
            {
                "key": "rsv_assembly_coverage",
                "label": "组装情况",
                "type": "paired",
                "items": [
                    {"label": "1x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_fraction"))},
                    {"label": "10x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_10x_fraction"))},
                    {"label": "100x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_100x_fraction"))},
                ],
            },
            {
                "key": "rsv_species_estimation",
                "label": "物种预估",
                "type": "paired",
                "items": [
                    {"label": "Nextclade Clade", "display": str(serotype_result.get("predicted_clade") or "--")},
                    {"label": "Lineage", "display": str(serotype_result.get("predicted_lineage") or "--")},
                ],
            },
        ]
    if str(serotype_result.get("mode") or "").strip() == "hmpv_nextclade":
        hmpv_overview_metrics = [
            {
                "key": "hmpv_assembly_coverage",
                "label": "组装情况",
                "type": "paired",
                "items": [
                    {"label": "1x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_fraction"))},
                    {"label": "10x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_10x_fraction"))},
                    {"label": "100x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_100x_fraction"))},
                ],
            },
            {
                "key": "hmpv_species_estimation",
                "label": "物种预估",
                "type": "paired",
                "items": [
                    {"label": "Nextclade Clade", "display": str(serotype_result.get("predicted_clade") or "--")},
                    {"label": "Lineage", "display": str(serotype_result.get("predicted_lineage") or "--")},
                ],
            },
        ]
    if str(serotype_result.get("mode") or "").strip() == "denv_nextclade":
        denv_overview_metrics = [
            {
                "key": "denv_assembly_coverage",
                "label": "组装情况",
                "type": "paired",
                "items": [
                    {"label": "1x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_fraction"))},
                    {"label": "10x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_10x_fraction"))},
                    {"label": "100x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_100x_fraction"))},
                ],
            },
            {
                "key": "denv_species_estimation",
                "label": "物种预估",
                "type": "paired",
                "items": [
                    {"label": "Nextclade Clade", "display": str(serotype_result.get("predicted_clade") or "--")},
                    {"label": "Lineage", "display": str(serotype_result.get("predicted_lineage") or "--")},
                ],
            },
        ]
    if str(serotype_result.get("mode") or "").strip() == "zikav_nextclade":
        zikav_overview_metrics = [
            {
                "key": "zikav_assembly_coverage",
                "label": "组装情况",
                "type": "paired",
                "items": [
                    {"label": "1x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_fraction"))},
                    {"label": "10x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_10x_fraction"))},
                    {"label": "100x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_100x_fraction"))},
                ],
            },
            {
                "key": "zikav_species_estimation",
                "label": "物种预估",
                "type": "paired",
                "items": [
                    {"label": "Nextclade Clade", "display": str(serotype_result.get("predicted_clade") or "--")},
                    {"label": "Lineage", "display": str(serotype_result.get("predicted_lineage") or "--")},
                ],
            },
        ]
    if str(serotype_result.get("mode") or "").strip() == "chikv_nextclade":
        chikv_overview_metrics = [
            {
                "key": "chikv_assembly_coverage",
                "label": "组装情况",
                "type": "paired",
                "items": [
                    {"label": "1x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_fraction"))},
                    {"label": "10x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_10x_fraction"))},
                    {"label": "100x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_100x_fraction"))},
                ],
            },
            {
                "key": "chikv_species_estimation",
                "label": "物种预估",
                "type": "paired",
                "items": [
                    {"label": "Nextclade Clade", "display": str(serotype_result.get("predicted_clade") or "--")},
                    {"label": "Lineage", "display": str(serotype_result.get("predicted_lineage") or "--")},
                ],
            },
        ]
    if str(serotype_result.get("mode") or "").strip() == "ebola_nextclade":
        ebola_overview_metrics = [
            {
                "key": "ebola_assembly_coverage",
                "label": "组装情况",
                "type": "paired",
                "items": [
                    {"label": "1x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_fraction"))},
                    {"label": "10x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_10x_fraction"))},
                    {"label": "100x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_100x_fraction"))},
                ],
            },
            {
                "key": "ebola_species_estimation",
                "label": "物种预估",
                "type": "paired",
                "items": [
                    {"label": "本地参考分型", "display": str(serotype_result.get("predicted_serotype") or "--")},
                    {"label": "Nextclade Clade", "display": str(serotype_result.get("predicted_clade") or "--")},
                    {"label": "Lineage / Genotype", "display": str(serotype_result.get("predicted_lineage") or "--")},
                    {"label": "参考序列", "display": str(serotype_result.get("reference_name") or "--")},
                ],
            },
        ]
    if str(serotype_result.get("mode") or "").strip() == "hpiv_typing":
        hpiv_overview_metrics = [
            {
                "key": "hpiv_assembly_coverage",
                "label": "组装情况",
                "type": "paired",
                "items": [
                    {"label": "1x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_fraction"))},
                    {"label": "10x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_10x_fraction"))},
                    {"label": "100x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_100x_fraction"))},
                ],
            },
            {
                "key": "hpiv_species_estimation",
                "label": "物种预估",
                "type": "paired",
                "items": [
                    {"label": "HPIV 亚型", "display": str(serotype_result.get("predicted_clade") or "--")},
                    {"label": "参考序列", "display": str(serotype_result.get("reference_name") or "--")},
                ],
            },
        ]
    if str(serotype_result.get("mode") or "").strip() == "hiv_resistance":
        hiv_overview_metrics = [
            {
                "key": "hiv_assembly_coverage",
                "label": "组装情况",
                "type": "paired",
                "items": [
                    {"label": "1x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_fraction"))},
                    {"label": "10x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_10x_fraction"))},
                    {"label": "100x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_100x_fraction"))},
                ],
            },
            {
                "key": "hiv_species_estimation",
                "label": "物种预估",
                "type": "paired",
                "items": [
                    {"label": "大亚型", "display": str(serotype_result.get("predicted_group") or "--")},
                    {"label": "子亚型", "display": str(serotype_result.get("predicted_clade") or "--")},
                    {"label": "重组判定", "display": str((serotype_result.get("summary_cards") or [{} ,{}, {} ,{}])[2].get("value") if isinstance(serotype_result.get("summary_cards"), list) and len(serotype_result.get("summary_cards")) >= 3 else "--")},
                    {"label": "代表株参考", "display": str((serotype_result.get("summary_cards") or [{} ,{}, {} ,{}])[3].get("value") if isinstance(serotype_result.get("summary_cards"), list) and len(serotype_result.get("summary_cards")) >= 4 else "--")},
                ],
            },
        ]
    if str(serotype_result.get("mode") or "").strip() == "hadv_typing":
        hadv_overview_metrics = [
            {
                "key": "hadv_assembly_coverage",
                "label": "组装情况",
                "type": "paired",
                "items": [
                    {"label": "1x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_fraction"))},
                    {"label": "10x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_10x_fraction"))},
                    {"label": "100x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_100x_fraction"))},
                ],
            },
            {
                "key": "hadv_species_estimation",
                "label": "物种预估",
                "type": "paired",
                "items": [
                    {"label": "HAdV 分型", "display": str(serotype_result.get("predicted_clade") or "--")},
                    {"label": "参考序列", "display": str(serotype_result.get("reference_name") or "--")},
                ],
            },
        ]
    if str(serotype_result.get("mode") or "").strip() == "norovirus_typing":
        norovirus_overview_metrics = [
            {
                "key": "norovirus_assembly_coverage",
                "label": "组装情况",
                "type": "paired",
                "items": [
                    {"label": "1x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_fraction"))},
                    {"label": "10x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_10x_fraction"))},
                    {"label": "100x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_100x_fraction"))},
                ],
            },
            {
                "key": "norovirus_species_estimation",
                "label": "物种预估",
                "type": "paired",
                "items": [
                    {"label": "双位点分型", "display": str(serotype_result.get("predicted_clade") or "--")},
                    {"label": "参考序列", "display": str(serotype_result.get("reference_name") or "--")},
                ],
            },
        ]
    if str(serotype_result.get("mode") or "").strip() == "enterovirus_typing":
        enterovirus_overview_metrics = [
            {
                "key": "rhinovirus_assembly_coverage",
                "label": "组装情况",
                "type": "paired",
                "items": [
                    {"label": "1x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_fraction"))},
                    {"label": "10x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_10x_fraction"))},
                    {"label": "100x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_100x_fraction"))},
                ],
            },
            {
                "key": "rhinovirus_species_estimation",
                "label": "物种预估",
                "type": "paired",
                "items": [
                    {"label": "VP1 分型", "display": str(serotype_result.get("predicted_clade") or "--")},
                    {"label": "大亚型", "display": str(serotype_result.get("predicted_group") or "--")},
                ],
            },
        ]
    if str(serotype_result.get("mode") or "").strip() == "hepatovirus_typing":
        hepatovirus_overview_metrics = [
            {
                "key": "hepatovirus_assembly_coverage",
                "label": "组装情况",
                "type": "paired",
                "items": [
                    {"label": "1x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_fraction"))},
                    {"label": "10x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_10x_fraction"))},
                    {"label": "100x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_100x_fraction"))},
                ],
            },
            {
                "key": "hepatovirus_species_estimation",
                "label": "物种预估",
                "type": "paired",
                "items": [
                    {"label": "大亚型", "display": str(serotype_result.get("predicted_group") or "--")},
                    {"label": "子亚型", "display": str(serotype_result.get("predicted_clade") or "--")},
                    {"label": "参考序列", "display": str(serotype_result.get("reference_name") or "--")},
                ],
            },
        ]
    if str(serotype_result.get("mode") or "").strip() == "bandavirus_typing":
        bandavirus_overview_metrics = [
            {
                "key": "bandavirus_assembly_coverage",
                "label": "组装情况",
                "type": "paired",
                "items": [
                    {"label": "1x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_fraction"))},
                    {"label": "10x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_10x_fraction"))},
                    {"label": "100x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_100x_fraction"))},
                ],
            },
            {
                "key": "bandavirus_species_estimation",
                "label": "物种预估",
                "type": "paired",
                "items": [
                    {"label": "大亚型", "display": str(serotype_result.get("predicted_group") or "--")},
                    {"label": "A_F 分型", "display": str(serotype_result.get("predicted_clade") or "--")},
                    {"label": "CJ 分型", "display": str(serotype_result.get("predicted_lineage") or "--")},
                ],
            },
        ]
    if str(serotype_result.get("mode") or "").strip() == "orthohantavirus_typing":
        orthohantavirus_overview_metrics = [
            {
                "key": "orthohantavirus_assembly_coverage",
                "label": "组装情况",
                "type": "paired",
                "items": [
                    {"label": "1x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_fraction"))},
                    {"label": "10x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_10x_fraction"))},
                    {"label": "100x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_100x_fraction"))},
                ],
            },
            {
                "key": "orthohantavirus_species_estimation",
                "label": "物种预估",
                "type": "paired",
                "items": [
                    {"label": "Orthohantavirus 分型", "display": str(serotype_result.get("predicted_clade") or "--")},
                    {"label": "S片段分型", "display": str(serotype_result.get("predicted_group") or "--")},
                ],
            },
        ]
    if str(serotype_result.get("mode") or "").strip() == "astroviridae_typing":
        astroviridae_overview_metrics = [
            {
                "key": "astroviridae_assembly_coverage",
                "label": "组装情况",
                "type": "paired",
                "items": [
                    {"label": "1x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_fraction"))},
                    {"label": "10x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_10x_fraction"))},
                    {"label": "100x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_100x_fraction"))},
                ],
            },
            {
                "key": "astroviridae_species_estimation",
                "label": "物种预估",
                "type": "paired",
                "items": [
                    {"label": "ORF2 分型", "display": str(serotype_result.get("predicted_clade") or "--")},
                    {"label": "病毒属", "display": str(serotype_result.get("predicted_group") or "--")},
                    {"label": "参考序列", "display": str(serotype_result.get("reference_name") or "--")},
                ],
            },
        ]
    if str(serotype_result.get("mode") or "").strip() == "rhinovirus_typing":
        rhinovirus_overview_metrics = [
            {
                "key": "rhinovirus_assembly_coverage",
                "label": "组装情况",
                "type": "paired",
                "items": [
                    {"label": "1x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_fraction"))},
                    {"label": "10x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_10x_fraction"))},
                    {"label": "100x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_100x_fraction"))},
                ],
            },
            {
                "key": "rhinovirus_species_estimation",
                "label": "物种预估",
                "type": "paired",
                "items": [
                    {"label": "VP1 分型", "display": str(serotype_result.get("predicted_clade") or "--")},
                    {"label": "参考序列", "display": str(serotype_result.get("reference_name") or "--")},
                ],
            },
        ]
    if str(serotype_result.get("mode") or "").strip() == "seasonal_hcov_typing":
        seasonal_hcov_overview_metrics = [
            {
                "key": "seasonal_hcov_assembly_coverage",
                "label": "组装情况",
                "type": "paired",
                "items": [
                    {"label": "1x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_fraction"))},
                    {"label": "10x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_10x_fraction"))},
                    {"label": "100x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_100x_fraction"))},
                ],
            },
            {
                "key": "seasonal_hcov_species_estimation",
                "label": "物种预估",
                "type": "paired",
                "items": [
                    {"label": "大类分型", "display": str(serotype_result.get("predicted_clade") or "--")},
                    {"label": "S 子亚型", "display": str(serotype_result.get("predicted_subtype") or "--")},
                    {"label": "参考序列", "display": str(serotype_result.get("reference_name") or "--")},
                ],
            },
        ]
    if str(serotype_result.get("mode") or "").strip() == "rotavirus_typing":
        rotavirus_overview_metrics = [
            {
                "key": "rotavirus_assembly_coverage",
                "label": "组装情况",
                "type": "paired",
                "items": [
                    {"label": "1x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_fraction"))},
                    {"label": "10x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_10x_fraction"))},
                    {"label": "100x 覆盖度", "display": _display_percent(assembly_coverage.get("coverage_100x_fraction"))},
                ],
            },
            {
                "key": "rotavirus_species_estimation",
                "label": "物种预估",
                "type": "paired",
                "items": [
                    {"label": "大组分型", "display": str(serotype_result.get("predicted_group") or "--")},
                    {"label": "组合分型", "display": str(serotype_result.get("predicted_subtype") or "--")},
                    {"label": "最优参考株", "display": str(serotype_result.get("reference_name") or "--")},
                ],
            },
        ]
    meta_overview_metrics = [
        {
            "key": "meta_sequencing",
            "label": "测序数据量 / 读取记录",
            "type": "paired",
            "items": [
                {"label": "总测序数据量", "display": _human_bp(summary_info.get("sum_len") or fastp_before_summary.get("total_bases"))},
                {"label": "读取记录", "display": _human_count(total_reads)},
            ],
        },
        {
            "key": "meta_assembly",
            "label": "组装情况",
            "type": "paired",
            "items": [
                {"label": "Contig数量", "display": _human_count(assembly_profile.get("contig_count"))},
                {"label": "总长度", "display": _human_bp(assembly_profile.get("total_length"))},
            ],
        },
        {
            "key": "meta_binning_quality",
            "label": "Binning完整性 / 污染率",
            "type": "paired",
            "items": [
                {"label": "平均完整性", "display": f"{binning_quality.get('summary', {}).get('avg_completeness', '--')}%" if binning_quality.get("summary", {}).get("avg_completeness") is not None else "--"},
                {"label": "平均污染率", "display": f"{binning_quality.get('summary', {}).get('avg_contamination', '--')}%" if binning_quality.get("summary", {}).get("avg_contamination") is not None else "--"},
            ],
        },
        {
            "key": "meta_species_mge",
            "label": "优势物种 / 移动元件",
            "type": "paired",
            "items": [
                {"label": "优势物种", "display": dominant_species.get("species") or "--"},
                {"label": "移动元件数量", "display": _human_count((mge_monitoring.get("overview") or {}).get("total_hits"))},
            ],
        },
    ]
    if dominant_virus_taxonomy:
        meta_overview_metrics.append(
            {
                "key": "virus_taxonomy",
                "label": "主导病毒 Taxonomy",
                "type": "paired",
                "items": [
                    {"label": "NCBI TaxID", "display": dominant_virus_taxonomy.get("taxid") or "--"},
                    {"label": "NCBI学名", "display": dominant_virus_taxonomy.get("scientific_name") or "--"},
                    {
                        "label": "科 / 属",
                        "display": " / ".join([
                            dominant_virus_taxonomy.get("family") or "-",
                            dominant_virus_taxonomy.get("genus") or "-",
                        ]),
                    },
                    {"label": "种", "display": dominant_virus_taxonomy.get("species_rank") or dominant_virus_taxonomy.get("species") or "--"},
                ],
            }
        )
    if dominant_fungus_taxonomy:
        meta_overview_metrics.append(
            {
                "key": "fungus_taxonomy",
                "label": "主导真菌 Taxonomy",
                "type": "paired",
                "items": [
                    {"label": "NCBI TaxID", "display": dominant_fungus_taxonomy.get("taxid") or "--"},
                    {"label": "NCBI学名", "display": dominant_fungus_taxonomy.get("scientific_name") or "--"},
                    {
                        "label": "科 / 属",
                        "display": " / ".join([
                            dominant_fungus_taxonomy.get("family") or "-",
                            dominant_fungus_taxonomy.get("genus") or "-",
                        ]),
                    },
                    {"label": "种", "display": dominant_fungus_taxonomy.get("species_rank") or dominant_fungus_taxonomy.get("species") or "--"},
                ],
            }
        )
    virus_overview_metrics = (
        influenza_overview_metrics
        or monkeypox_overview_metrics
        or rsv_overview_metrics
        or hmpv_overview_metrics
        or denv_overview_metrics
        or zikav_overview_metrics
        or chikv_overview_metrics
        or ebola_overview_metrics
        or hpiv_overview_metrics
        or hiv_overview_metrics
        or hadv_overview_metrics
        or norovirus_overview_metrics
        or enterovirus_overview_metrics
        or hepatovirus_overview_metrics
        or bandavirus_overview_metrics
        or orthohantavirus_overview_metrics
        or astroviridae_overview_metrics
        or rhinovirus_overview_metrics
        or seasonal_hcov_overview_metrics
        or rotavirus_overview_metrics
    )
    overview_metrics = meta_overview_metrics if is_meta_method else [
        {"key": "total_bases", "label": "总测序数据量", "type": "single", "value": resolved_total_bases, "display": _human_bp(resolved_total_bases), "unit": "bp"},
        *(virus_overview_metrics[:1] if virus_overview_metrics else []),
        *(
            []
            if virus_overview_metrics
            else [
                {
                    "key": "assembly_profile",
                    "label": "组装情况",
                    "type": "assembly_profile",
                    "contig_count": assembly_profile.get("contig_count"),
                    "plasmid_count": assembly_profile.get("plasmid_count"),
                    "total_count": assembly_profile.get("total_count"),
                    "total_length": assembly_profile.get("total_length"),
                }
            ]
        ),
        {
            "key": "q_metrics",
            "label": "Q20 / Q30",
            "type": "paired",
            "items": [
                {"label": "Q20", "display": _display_percent(resolved_q20_rate)},
                {"label": "Q30", "display": _display_percent(resolved_q30_rate)},
            ],
        },
        *(
            []
            if virus_overview_metrics
            else [
                {
                    "key": "checkm_metrics",
                    "label": "完整性 / 污染率",
                    "type": "paired",
                    "items": [
                        {"label": "完整性", "display": _display_percent_points(checkm_info.get("completeness"))},
                        {"label": "污染率", "display": _display_percent_points(checkm_info.get("contamination"))},
                    ],
                }
            ]
        ),
        *(virus_overview_metrics[1:] if virus_overview_metrics else []),
        *(
            []
            if virus_overview_metrics
            else [
                {
                    "key": "species_estimation",
                    "label": "物种预估",
                    "type": "paired",
                    "items": [
                        {"label": "物种名称", "display": checkm_info.get("species_name") or "--"},
                        {"label": "MLST 物种名称", "display": checkm_info.get("mlst_species_name") or "--"},
                    ],
                }
            ]
        ),
        *(
            [
                {
                    "key": "virus_taxonomy",
                    "label": "主导病毒 Taxonomy",
                    "type": "paired",
                    "items": [
                        {"label": "NCBI TaxID", "display": dominant_virus_taxonomy.get("taxid") or "--"},
                        {"label": "NCBI学名", "display": dominant_virus_taxonomy.get("scientific_name") or "--"},
                        {
                            "label": "科 / 属",
                            "display": " / ".join([
                                dominant_virus_taxonomy.get("family") or "-",
                                dominant_virus_taxonomy.get("genus") or "-",
                            ]),
                        },
                        {"label": "种", "display": dominant_virus_taxonomy.get("species_rank") or dominant_virus_taxonomy.get("species") or "--"},
                    ],
                }
            ] if dominant_virus_taxonomy else []
        ),
        *(
            [
                {
                    "key": "fungus_taxonomy",
                    "label": "主导真菌 Taxonomy",
                    "type": "paired",
                    "items": [
                        {"label": "NCBI TaxID", "display": dominant_fungus_taxonomy.get("taxid") or "--"},
                        {"label": "NCBI学名", "display": dominant_fungus_taxonomy.get("scientific_name") or "--"},
                        {
                            "label": "科 / 属",
                            "display": " / ".join([
                                dominant_fungus_taxonomy.get("family") or "-",
                                dominant_fungus_taxonomy.get("genus") or "-",
                            ]),
                        },
                        {"label": "种", "display": dominant_fungus_taxonomy.get("species_rank") or dominant_fungus_taxonomy.get("species") or "--"},
                    ],
                }
            ] if dominant_fungus_taxonomy else []
        ),
    ]
    overview_status = "ready" if overview_metrics else "empty"
    multi_sample_summary = _build_multi_sample_queue_summary(report_source)
    assembly_status = "ready" if (
        assembly_summary.get("rows")
        or contig_annotation.get("rows")
        or checkm_quality.get("rows")
        or gene_annotation_summary.get("rows")
        or assembly_coverage.get("status") == "ready"
        or contig_depth_relationship.get("status") == "ready"
        or gene_length_distribution.get("status") == "ready"
        or cgview_assets.get("status") == "ready"
    ) else "empty"
    resistance_virulence_status = "ready" if (
        rv_overview.get("status") == "ready"
        or rv_summary.get("rows")
        or virulence_elements.get("rows")
        or resistance_elements.get("rows")
    ) else "empty"
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
            "asm_type": params.get("asm_type", ""),
            "method": params.get("method", ""),
            "analysis_target": params.get("analysis_target", "bacteria"),
            "species": params.get("species", ""),
            "sample_name": sample_name,
            "sample_display_name": sample_display_name,
            "samples": report_source.get("samples", []),
            "report_mode": report_source.get("mode", "single"),
            "multi_sample_summary": multi_sample_summary,
        },
        "overview_metrics": overview_metrics,
        "sections": {
            "overview": {"status": overview_status, "multi_sample_summary": multi_sample_summary},
            "raw_qc": {
                "status": "ready" if fastp_info else "empty",
                "paired_end": {
                    "left": fastp_info.get("read1", {"status": "empty"}),
                    "right": fastp_info.get("read2", {"status": "empty"}),
                },
                "fastp": fastp_info.get("summary", {"status": "empty", "plots": []}),
            },
            "species_identification": {
                "status": "ready" if species_taxonomy.get("rows") or subspecies_taxonomy.get("rows") or assembly_taxonomy.get("rows") else "empty",
                "species": species_taxonomy,
                "subspecies": subspecies_taxonomy,
                "abundance": taxonomy_abundance,
                "risk_summary": taxonomy_risk_summary,
                "interpretation": taxonomy_interpretation,
                "rarefaction": _build_taxonomy_rarefaction(species_taxonomy, subspecies_taxonomy),
                "assembly_taxonomy": assembly_taxonomy,
            },
            "binning_results": {
                "status": "ready" if binning_quality.get("status") == "ready" or binning_taxonomy.get("status") == "ready" or viral_assembly.get("status") == "ready" else "empty",
                "quality": binning_quality,
                "taxonomy": binning_taxonomy,
                "viral_assembly": viral_assembly,
            },
            "assembly": {
                "status": assembly_status,
                "summary": assembly_summary,
                "coverage": assembly_coverage,
                "contig_annotation": contig_annotation,
                "contig_depth_relationship": contig_depth_relationship,
                "cgview": cgview_assets,
                "checkm": checkm_quality,
                "gene_annotation_summary": gene_annotation_summary,
                "gene_length_distribution": gene_length_distribution,
            },
            "resistance_virulence": {
                "status": resistance_virulence_status,
                "overview": rv_overview,
                "summary": rv_summary,
                "virulence_elements": virulence_elements,
                "virulence_relationship": virulence_relationship,
                "resistance_elements": resistance_elements,
                "resistance_relationship": resistance_relationship,
            },
            "mlst": {
                "status": "ready" if mlst_result.get("rows") else "empty",
                "columns": mlst_result.get("columns", []),
                "rows": mlst_result.get("rows", []),
                "gene_show_map": mlst_result.get("gene_show_map", {}),
                "default_gene": mlst_result.get("default_gene", ""),
                "knowledge_summary": mlst_result.get("knowledge_summary", {"headline": "", "items": []}),
                "title": mlst_result.get("title", ""),
                "tag_label": mlst_result.get("tag_label", ""),
                "empty_message": mlst_result.get("empty_message", ""),
                "detail_empty_message": mlst_result.get("detail_empty_message", ""),
                "generic_detail_note": mlst_result.get("generic_detail_note", ""),
                "neisseria_amr": neisseria_amr_result,
            },
            "tb_amr": tb_amr_result,
            "serotype": serotype_result,
            "priority_serotype": {
                "status": "ready" if priority_serotype.get("rows") else "empty",
                "columns": priority_serotype.get("columns", []),
                "rows": priority_serotype.get("rows", []),
            },
            "mge_monitoring": mge_monitoring,
            "public_health_support": public_health_support,
            "knowledge_interpretation": knowledge_interpretation,
            "modeling_risk": build_modeling_risk_section(project_root, sample_name=sample_name, report_dir=report_dir),
        },
    }
    payload["sections"]["workflow_closure"] = build_workflow_closure(payload)
    payload["sections"]["export_checklist"] = build_export_checklist(payload)
    _write_report_cache(report_dir, sample_name, cache_fingerprint, payload)
    return payload

def _read_meta_annotation_subset(path: Path) -> dict:
    raw = _read_tsv_rows(path)
    columns = raw.get("columns") or []
    rows = raw.get("rows") or []
    if not columns or not rows:
        return {"columns": [], "rows": []}
    selected_columns = [columns[0], columns[1], columns[-1]]
    selected_indexes = [0, 1, len(columns) - 1]
    return {
        "columns": ["Contig名称", "质粒基因组预测", "TPM"],
        "rows": [[row[index] if index < len(row) else "" for index in selected_indexes] for row in rows],
    }

def _read_fastp_metrics(path: Path) -> dict:
    if not path.is_file():
        return {}
    try:
        payload = json.loads(path.read_text(encoding="utf-8", errors="ignore"))
    except (OSError, json.JSONDecodeError):
        return {}
    summary = payload.get("summary", {})
    before = summary.get("before_filtering", {})
    after = summary.get("after_filtering", {})
    filtering = payload.get("filtering_result", {})
    duplication = payload.get("duplication", {})
    read1_before = payload.get("read1_before_filtering", {})
    read2_before = payload.get("read2_before_filtering", {})
    read1_after = payload.get("read1_after_filtering", {})
    read2_after = payload.get("read2_after_filtering", {})
    return {
        "summary": {
            "status": "ready",
            "sequencing": summary.get("sequencing", ""),
            "before_filtering": before,
            "after_filtering": after,
            "filtering_result": filtering,
            "duplication_rate": duplication.get("rate"),
            "insert_size": payload.get("insert_size", {}),
            "adapter_cutting": payload.get("adapter_cutting", {}),
            "base_distribution": {
                "read1": read1_after.get("content_curves", {}),
                "read2": read2_after.get("content_curves", {}),
            },
        },
        "read1": {
            "status": "ready",
            "label": "R1",
            "before": read1_before,
            "after": read1_after,
            "quality_curves": read1_before.get("quality_curves", {}),
            "content_curves": read1_before.get("content_curves", {}),
            "before_summary": _summarize_fastp_read_block(read1_before),
            "after_summary": _summarize_fastp_read_block(read1_before),
        },
        "read2": {
            "status": "ready",
            "label": "R2",
            "before": read2_before,
            "after": read2_after,
            "quality_curves": read2_before.get("quality_curves", {}),
            "content_curves": read2_before.get("content_curves", {}),
            "before_summary": _summarize_fastp_read_block(read2_before),
            "after_summary": _summarize_fastp_read_block(read2_before),
        },
    }

def _summarize_fastp_read_block(section: dict) -> dict:
    total_reads = _safe_int(section.get("total_reads"))
    total_bases = _safe_int(section.get("total_bases"))
    q20_bases = _safe_int(section.get("q20_bases"))
    q30_bases = _safe_int(section.get("q30_bases"))
    content_curves = section.get("content_curves", {}) or {}
    gc_curve = content_curves.get("GC", []) or []
    mean_length = round(total_bases / total_reads, 2) if total_reads and total_bases else None
    q20_rate = round(q20_bases / total_bases, 6) if q20_bases is not None and total_bases else None
    q30_rate = round(q30_bases / total_bases, 6) if q30_bases is not None and total_bases else None
    gc_content = round(sum((_safe_float(value) or 0) for value in gc_curve) / len(gc_curve), 6) if gc_curve else None
    return {
        "total_reads": total_reads,
        "total_bases": total_bases,
        "mean_length": mean_length,
        "q20_rate": q20_rate,
        "q30_rate": q30_rate,
        "gc_content": gc_content,
    }
