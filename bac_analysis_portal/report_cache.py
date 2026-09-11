from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path

from .report_artifacts import _resolve_report_artifact_path
from .report_modeling import build_modeling_risk_section
from .report_sources import _build_multi_sample_queue_summary
from .runtime_paths import _resolve_runtime_database_root

REPORT_CACHE_VERSION = 137


def _report_cache_dir(report_dir: Path) -> Path:
    return report_dir / ".portal_report_cache"

def _report_cache_path(report_dir: Path, sample_name: str) -> Path:
    safe_sample = str(sample_name or "default").strip() or "default"
    safe_sample = "".join(ch if ch.isalnum() or ch in {"-", "_", "."} else "_" for ch in safe_sample)
    return _report_cache_dir(report_dir) / f"report_payload_{safe_sample}.json"

def _report_cache_fingerprint(task: dict, report_dir: Path, sample_name: str) -> str:
    params = task.get("params") or {}
    database_root = _resolve_runtime_database_root()
    knowledge_base_root = database_root / "knowledge_base"
    sources = [
        Path(__file__).resolve(),
        report_dir / "community_summary.json",
        report_dir / "community_command_plan.tsv",
        report_dir / "community_metadata_preview.tsv",
        report_dir / "community_demux_summary.tsv",
        report_dir / "community_taxonomy_preview.tsv",
        report_dir / "taxonomy_export" / "taxonomy.tsv",
        report_dir / "ancombc2_export" / "q.jsonl",
        report_dir / "ancombc2_export" / "diff.jsonl",
        report_dir / "ancombc2_export" / "lfc.jsonl",
        report_dir / "microeco_beta" / "pcoa_plot.png",
        report_dir / "microeco_beta" / "nmds_plot.png",
        report_dir / "microeco_beta" / "within_group_distance_plot.png",
        report_dir / "microeco_beta" / "pcoa_scores.tsv",
        report_dir / "microeco_beta" / "nmds_scores.tsv",
        report_dir / "microeco_beta" / "permanova.tsv",
        report_dir / "microeco_beta" / "anosim.tsv",
        report_dir / "microeco_beta" / "betadisper.txt",
        report_dir / "microeco_beta" / "within_group_distance_stats.tsv",
        report_dir / "microeco_beta" / "group_sample_counts.tsv",
        report_dir / "microeco_beta" / "run_summary.txt",
        report_dir / "microeco_biomarker" / "lefse_diff.tsv",
        report_dir / "microeco_biomarker" / "lefse_barplot.png",
        report_dir / "microeco_biomarker" / "lefse_barplot.pdf",
        report_dir / "microeco_biomarker" / "rf_importance.tsv",
        report_dir / "microeco_biomarker" / "rf_importance.png",
        report_dir / "microeco_biomarker" / "rf_importance.pdf",
        report_dir / "microeco_biomarker" / "run_summary.txt",
        report_dir / "microeco_network" / "network_summary.tsv",
        report_dir / "microeco_network" / "node_table.tsv",
        report_dir / "microeco_network" / "edge_table.tsv",
        report_dir / "microeco_network" / "module_summary.tsv",
        report_dir / "microeco_network" / "role_summary.tsv",
        report_dir / "microeco_network" / "eigen_summary.tsv",
        report_dir / "microeco_network" / "run_summary.txt",
        report_dir / "summary.tsv",
        report_dir / "Assem_info.tsv",
        report_dir / "Assem_info1.tsv",
        report_dir / "binning_name.tsv",
        report_dir / "tmp_combine.fa",
        report_dir / "meta_plas_vf_card.tsv",
        report_dir / "viral_assembly" / "viral_summary.tsv",
        report_dir / "viral_assembly" / "viral_contig_summary.tsv",
        report_dir / "viral_assembly" / "viral_retained_contigs.fa",
        report_dir / "viral_assembly" / "megahit_output" / "final.contigs.fa",
        report_dir / "viral_assembly" / "virsorter2" / "final-viral-score.tsv",
        report_dir / "viral_assembly" / "checkv" / "contamination.tsv",
        report_dir / "viral_assembly" / "checkv" / "quality_summary.tsv",
        report_dir / "bin_vfdb.tsv",
        report_dir / "bin_card.tsv",
        report_dir / "binning_rgi_new.txt",
        report_dir / "staramr_result" / "resfinder.tsv",
        database_root / "who_bppl_2024_support.json",
        database_root / "china_cdc_bacteria_support.json",
        report_dir / "flye_output" / "assembly_info.txt",
        report_dir / "checkm2_out" / "quality_report.tsv",
        report_dir / "bin_checkm2out" / "quality_report.tsv",
        report_dir / "gtdbtk_out" / "gtdbtk.bac120.summary.tsv",
        report_dir / f"{sample_name}.fastp2.json",
        report_dir / f"{sample_name}.assemble.result.tsv",
        report_dir / f"{sample_name}_ngs.per-base.bed.gz",
        report_dir / f"{sample_name}_ngs.per-base.bed",
        report_dir / f"{sample_name}.regions.bed",
        report_dir / f"{sample_name}.genefun_summary.tsv",
        report_dir / f"{sample_name}_gene_raw_sum.tsv",
        report_dir / f"{sample_name}_2.list.txt",
        report_dir / f"{sample_name}_2.list2.txt",
        report_dir / f"{sample_name}_assem.kraken2.txt",
        report_dir / f"{sample_name}.mlst_Stat.txt",
        report_dir / f"{sample_name}.neisseria_amr_calls.csv",
        report_dir / "tb_analysis" / "tb_summary.json",
        report_dir / "tb_analysis" / "tb_catalogue_matches.tsv",
        report_dir / "tb_analysis" / "tbprofiler" / f"{sample_name}.results.json",
        report_dir / "tb_analysis" / "reference_call" / "snps.filt1.vcf",
        report_dir / "tb_analysis" / "reference_call" / "snps.anno.vcf",
        report_dir / "tb_analysis" / "reference_call" / f"{sample_name}.anno.tsv",
        report_dir / f"{sample_name}_serotype_result.tsv",
        report_dir / f"{sample_name}.pathonet_result.tsv",
        report_dir / f"{sample_name}.mge_risk_summary.tsv",
        report_dir / f"{sample_name}.integrated_mge_summary.tsv",
        report_dir / "nextclade_output" / "nextclade.tsv",
        report_dir / "nextclade_output" / "nextclade.json",
        report_dir / "nextclade_output" / "nextclade.csv",
        report_dir / "nextclade_output" / "nextclade.auspice.json",
        report_dir / "nextclade_output" / "nextclade.nwk",
        report_dir / "wf_flu" / "typing_summary.tsv",
        report_dir / "wf_flu" / "nextclade" / "segment_analysis.tsv",
        database_root / "virus" / "ncov" / "genomic.gff",
        report_dir / "genomes" / "ref.fa",
        report_dir / "genomes" / "ref.fa.fai",
        report_dir / "ref.mapping.bam",
        report_dir / "ref.mapping.bam.bai",
        report_dir / "ref" / "genes.gff",
        report_dir / "snps.raw.vcf",
        report_dir / "snps.filt1.vcf",
        report_dir / "snps.anno.vcf",
        report_dir / "snps.raw.mutation_table.json",
        report_dir / "snps.raw.mutation_table.tsv",
        report_dir / "snps.filt1.mutation_table.json",
        report_dir / "snps.filt1.mutation_table.tsv",
        report_dir / "snps.filt1.resistance_annotation.json",
        report_dir / "snps.filt1.resistance_annotation.tsv",
        report_dir / f"{sample_name}_rsv_reference_selection" / "selection.tsv",
        report_dir / f"{sample_name}_denv_reference_selection" / "selection.tsv",
        report_dir / f"{sample_name}_hpiv_reference_selection" / "selection.tsv",
        report_dir / f"{sample_name}_hadv_reference_selection" / "selection.tsv",
        report_dir / f"{sample_name}_hadv_reference_selection" / "phf_typing" / "phf_typing.tsv",
        report_dir / f"{sample_name}_enterovirus_reference_selection" / "selection.tsv",
        report_dir / f"{sample_name}_enterovirus_reference_selection" / "consensus_typing" / "consensus_typing.tsv",
        report_dir / f"{sample_name}_hepatovirus_reference_selection" / "selection.tsv",
        report_dir / f"{sample_name}_hepatovirus_reference_selection" / "consensus_typing" / "consensus_typing.tsv",
        report_dir / f"{sample_name}_bandavirus_reference_selection" / "selection.tsv",
        report_dir / f"{sample_name}_bandavirus_reference_selection" / "consensus_typing.tsv",
        report_dir / f"{sample_name}_bandavirus_reference_selection" / "selected_segments.tsv",
        report_dir / f"{sample_name}_bandavirus_reference_selection" / "cj_typing" / "selected_segments.tsv",
        report_dir / f"{sample_name}_orthohantavirus_reference_selection" / "selection.tsv",
        report_dir / f"{sample_name}_orthohantavirus_reference_selection" / "consensus_typing.tsv",
        report_dir / f"{sample_name}_orthohantavirus_reference_selection" / "selected_segments.tsv",
        report_dir / f"{sample_name}_orthohantavirus_reference_selection" / "snpeff_reference.gff3",
        report_dir / f"{sample_name}_orthohantavirus_reference_selection" / "snpeff_vadr_ref" / "ref" / "genes.gff",
        report_dir / f"{sample_name}_orthoebolavirus_reference_selection" / "selection.tsv",
        report_dir / f"{sample_name}_orthoebolavirus_reference_selection" / "consensus_typing" / "consensus_typing.tsv",
        report_dir / f"{sample_name}_astroviridae_reference_selection" / "selection.tsv",
        report_dir / f"{sample_name}_astroviridae_reference_selection" / "consensus_typing" / "consensus_typing.tsv",
        report_dir / f"{sample_name}_astroviridae_reference_selection" / "phylogeny" / "summary.tsv",
        report_dir / f"{sample_name}_rhinovirus_reference_selection" / "selection.tsv",
        report_dir / f"{sample_name}_rhinovirus_reference_selection" / "consensus_typing" / "consensus_typing.tsv",
        report_dir / f"{sample_name}_rhinovirus_reference_selection" / "phylogeny" / "summary.tsv",
        report_dir / "hpiv_coverage" / "hpiv.coverage.summary.tsv",
        report_dir / "hpiv_coverage" / "hpiv.coverage.regions.bed",
        database_root / "virus" / "rsv" / "nmdc_hrsv_variation_list.tsv",
        database_root / "virus" / "rsv" / "nmdc_hrsv_variation_list_zh.tsv",
        report_dir / "vadr" / f"{sample_name}.vadr.gff3",
        report_dir / f"{sample_name}_virus_typing" / "vadr" / f"{sample_name}.vadr.gff3",
        report_dir / "Cluster.tsv",
        report_dir / "dis_bin.tsv",
        report_dir / "dis.mat.txt",
        report_dir / "Full_ANI.txt",
        report_dir / "grapetree.nwk",
        report_dir / "mlst.nwk",
        report_dir / "mlst.txt",
        report_dir / "rmref.core.aln.contree",
    ]
    if knowledge_base_root.is_dir():
        for path in sorted(knowledge_base_root.rglob("*.json"), key=lambda item: str(item).lower()):
            sources.append(path)
    fallback_artifacts = [
        _resolve_report_artifact_path(report_dir, [f"{sample_name}.fastp2.json"] if sample_name else [], ["*.fastp2.json"]),
        _resolve_report_artifact_path(
            report_dir,
            [f"{sample_name}_ngs.per-base.bed.gz", f"{sample_name}_ngs.per-base.bed", f"{sample_name}.regions.bed"] if sample_name else [],
            ["*_ngs.per-base.bed.gz", "*_ngs.per-base.bed", "*.regions.bed.gz", "*.regions.bed", "*.per-base.bed.gz", "*.per-base.bed"],
        ),
    ]
    for artifact in fallback_artifacts:
        if artifact and artifact not in sources:
            sources.append(artifact)
    if sample_name:
        prokka_dir = report_dir / f"{sample_name}_prokka"
        if prokka_dir.is_dir():
            faa_path = prokka_dir / f"{sample_name}.faa"
            if faa_path.is_file():
                sources.append(faa_path)
            sources.extend(
                sorted(
                    [
                        child
                        for child in prokka_dir.iterdir()
                        if child.is_file() and child.suffix.lower() == ".gbk" and not child.name.startswith(".")
                    ],
                    key=lambda item: item.name.lower(),
                )
            )
    site_table = Path(__file__).resolve().parent.parent / "database" / "NM_mutate" / "neisseria_meningitidis_snp_amr_associations_literature_updated.csv"
    if site_table.is_file():
        sources.append(site_table)
    stamp = {
        "version": REPORT_CACHE_VERSION,
        "task_id": str(task.get("id") or ""),
        "task_status": str(task.get("status") or ""),
        "method": str(params.get("method") or ""),
        "sample_name": sample_name,
        "files": [],
    }
    for path in sources:
        try:
            if path.exists():
                stat = path.stat()
                stamp["files"].append({
                    "path": str(path.relative_to(report_dir)) if path.is_relative_to(report_dir) else str(path),
                    "size": stat.st_size,
                    "mtime_ns": stat.st_mtime_ns,
                })
        except OSError:
            continue
    return hashlib.sha256(json.dumps(stamp, ensure_ascii=False, sort_keys=True).encode("utf-8")).hexdigest()

def _read_report_cache(report_dir: Path, sample_name: str, fingerprint: str) -> dict | None:
    cache_path = _report_cache_path(report_dir, sample_name)
    if not cache_path.is_file():
        return None
    try:
        payload = json.loads(cache_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if int(payload.get("version") or 0) != REPORT_CACHE_VERSION:
        return None
    if str(payload.get("fingerprint") or "") != fingerprint:
        return None
    return payload.get("data") if isinstance(payload.get("data"), dict) else None

def _refresh_cached_report_payload(
    payload: dict,
    *,
    task: dict,
    report_dir: Path,
    report_source: dict,
    sample_name: str,
    sample_display_name: str,
) -> dict:
    refreshed = copy.deepcopy(payload)
    cached_task = refreshed.get("task") if isinstance(refreshed.get("task"), dict) else {}
    params = task.get("params") or {}
    root_dir = report_source.get("root_dir")
    root_dir_text = str(root_dir) if root_dir else str(params.get("output_dir") or report_dir)
    cached_task.update({
        "id": task.get("id"),
        "name": task.get("name"),
        "status": task.get("status"),
        "owner": task.get("owner"),
        "group": task.get("owner_group", ""),
        "created_at": task.get("created_at"),
        "started_at": task.get("started_at"),
        "finished_at": task.get("finished_at"),
        "input_path": params.get("input_path", ""),
        "output_dir": root_dir_text,
        "sample_name": sample_name,
        "sample_display_name": sample_display_name,
        "samples": report_source.get("samples", []),
        "report_mode": report_source.get("mode", "single"),
        "multi_sample_summary": _build_multi_sample_queue_summary(report_source),
    })
    refreshed["task"] = cached_task
    refreshed_sections = refreshed.get("sections") if isinstance(refreshed.get("sections"), dict) else {}
    overview_section = refreshed_sections.get("overview") if isinstance(refreshed_sections.get("overview"), dict) else {}
    overview_section["multi_sample_summary"] = cached_task.get("multi_sample_summary") or {}
    refreshed_sections["overview"] = overview_section
    refreshed_sections["modeling_risk"] = build_modeling_risk_section(
        Path(__file__).resolve().parent.parent,
        sample_name=sample_name,
        report_dir=report_dir,
    )
    refreshed["sections"] = refreshed_sections
    return refreshed

def _make_json_safe(value: object) -> object:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _make_json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_make_json_safe(item) for item in value]
    if isinstance(value, set):
        return [_make_json_safe(item) for item in sorted(value, key=lambda item: str(item))]
    return value

def _write_report_cache(report_dir: Path, sample_name: str, fingerprint: str, data: dict) -> None:
    cache_dir = _report_cache_dir(report_dir)
    cache_path = _report_cache_path(report_dir, sample_name)
    temp_path = cache_path.with_suffix(".tmp")
    try:
        cache_dir.mkdir(parents=True, exist_ok=True)
        payload = _make_json_safe({"version": REPORT_CACHE_VERSION, "fingerprint": fingerprint, "data": data})
        temp_path.write_text(
            json.dumps(payload, ensure_ascii=False),
            encoding="utf-8",
        )
        temp_path.replace(cache_path)
    except OSError:
        try:
            if temp_path.exists():
                temp_path.unlink()
        except OSError:
            pass
