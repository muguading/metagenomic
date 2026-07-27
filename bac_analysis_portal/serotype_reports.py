from __future__ import annotations

import csv
import json
import re
from pathlib import Path

from .igv_assets import (
    _discover_influenza_igv_assets,
    _discover_hiv_igv_assets,
    _discover_ncov_igv_assets,
    _discover_monkeypox_igv_assets,
    _discover_rsv_igv_assets,
    _discover_hmpv_igv_assets,
    _discover_denv_igv_assets,
    _discover_zikav_igv_assets,
    _discover_chikv_igv_assets,
    _discover_ebola_igv_assets,
    _discover_hpiv_igv_assets,
    _discover_norovirus_igv_assets,
    _discover_enterovirus_igv_assets,
    _discover_hepatovirus_igv_assets,
    _discover_bandavirus_igv_assets,
    _discover_orthohantavirus_igv_assets,
    _discover_astroviridae_igv_assets,
    _discover_rhinovirus_igv_assets,
    _discover_seasonal_hcov_igv_assets,
    _discover_hadv_igv_assets,
)
from .knowledge_interpretation import _build_viral_serotype_knowledge_summary
from .phylogeny_reports import (
    _build_astroviridae_gene_phylogeny,
    _build_enterovirus_gene_phylogeny,
    _build_nextclade_phylogeny_tree,
    _build_norovirus_gene_phylogeny,
    _build_rhinovirus_gene_phylogeny,
    _build_seasonal_hcov_spike_phylogeny,
)
from .report_artifacts import _collect_ncov_ngdc_matches
from .table_io import _read_tsv_rows
from .variant_reports import (
    _build_influenza_mutation_display_table,
    _build_ncov_mutation_display_table,
    _read_hadv_phf_snp_section,
    _read_hpiv_functional_annotation_table,
    _read_influenza_resistance_annotation_table,
    _read_influenza_variant_annotation_table,
    _read_monkeypox_variant_annotation_table,
    _read_ncov_variant_annotation_table,
    _read_rsv_nmdc_annotation_table,
    _read_rsv_variant_annotation_table,
)

def _sanitize_virus_demo_note(value: object) -> str:
    text = str(value or "").strip()
    if not text or text == "-":
        return ""
    blocked_fragments = [
        "当前仅内置流感病毒分型流程",
        "当前仅支持流感病毒分型流程",
        "当前仅内置细菌分型流程",
        "当前仅支持细菌分型流程",
    ]
    if any(fragment in text for fragment in blocked_fragments):
        return ""
    return text

def _is_bordetella_pertussis(checkm_info: dict) -> bool:
    candidates = [
        checkm_info.get("species_name"),
        checkm_info.get("mlst_species_name"),
    ]
    for value in candidates:
        text = str(value or "").strip().lower()
        if "bordetella" in text and "pertussis" in text:
            return True
    return False

def _read_influenza_typing_section(report_dir: Path, sample_name: str) -> dict | None:
    summary_path = report_dir / "wf_flu" / "typing_summary.tsv"
    if not summary_path.is_file():
        return None
    try:
        with summary_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            first = next(reader, None)
    except OSError:
        return None
    if not first:
        return None
    status = str(first.get("status") or "").strip() or "-"
    influenza_type = str(first.get("influenza_type") or "").strip() or "-"
    ha_subtype = str(first.get("ha_subtype") or "").strip() or "-"
    na_subtype = str(first.get("na_subtype") or "").strip() or "-"
    subtype_call = str(first.get("subtype_call") or "").strip() or "-"
    reference_path = str(first.get("reference_path") or "").strip()
    segment_manifest = _read_tsv_rows(report_dir / "wf_flu" / "reference_sets" / "final_segments.tsv")
    variant_annotation = _read_influenza_variant_annotation_table(report_dir)
    resistance_annotation = _read_influenza_resistance_annotation_table(report_dir)
    igv_view = _discover_influenza_igv_assets(report_dir, sample_name)
    raw_mutation_table = variant_annotation if str(variant_annotation.get("status") or "") == "ready" else (
        _read_tsv_rows(report_dir / f"{sample_name}.anno.tsv") if sample_name else {"columns": [], "rows": []}
    )
    mutation_table = _build_influenza_mutation_display_table(report_dir, raw_mutation_table)
    summary_columns = ["样本名称", "流感类型", "HA亚型", "NA亚型", "分型结果", "状态"]
    summary_rows = [[sample_name, influenza_type, ha_subtype, na_subtype, subtype_call, status]]
    segment_count = len(segment_manifest.get("rows") or [])
    mutation_rows = mutation_table.get("rows") or []
    mutation_columns = mutation_table.get("columns") or []
    mutation_segment_count = 0
    if mutation_rows and mutation_columns:
        segment_column = "染色体" if "染色体" in mutation_columns else (mutation_columns[0] if mutation_columns else "")
        if segment_column in mutation_columns:
            segment_index = mutation_columns.index(segment_column)
            mutation_segment_count = len(
                {
                    str(row[segment_index]).strip()
                    for row in mutation_rows
                    if row and segment_index < len(row) and str(row[segment_index]).strip()
                }
            )
    summary_cards = [
        {"label": "流感类型", "value": influenza_type},
        {"label": "HA 亚型", "value": ha_subtype},
        {"label": "NA 亚型", "value": na_subtype},
        {"label": "分型结果", "value": subtype_call},
        {"label": "突变位点数", "value": int(variant_annotation.get("total_variants") or len(mutation_rows)) if mutation_rows else "--"},
    ]
    note_parts = []
    if segment_count:
        note_parts.append(f"已拼装 {segment_count} 个 segment 参考片段")
    if reference_path:
        note_parts.append(f"参考集合: {Path(reference_path).name}")
    if str(variant_annotation.get("status") or "") == "ready" and mutation_rows:
        note_parts.append("已基于 VADR GFF3 与 consensus FASTA 生成 snpEff 变异注释表")
    if status == "screening_stop":
        note_parts.append("初筛未满足后续流感组装条件")
    return {
        "status": "ready",
        "mode": "influenza_typing",
        "influenza_type": influenza_type,
        "ha_subtype": ha_subtype,
        "na_subtype": na_subtype,
        "predicted_serotype": subtype_call,
        "reference_path": reference_path,
        "summary_cards": summary_cards,
        "segment_manifest": segment_manifest,
        "mutation_table": mutation_table,
        "mutation_summary": {
            "count": int(variant_annotation.get("total_variants") or len(mutation_rows)),
            "segment_count": mutation_segment_count,
            "columns": mutation_columns,
        },
        "variant_annotation": variant_annotation,
        "resistance_annotation": resistance_annotation,
        "igv": igv_view,
        "notes": "；".join(note_parts),
        "columns": summary_columns,
        "rows": summary_rows,
    }

def _read_hadv_phf_coverage_section(report_dir: Path, sample_name: str) -> dict:
    phf_dir = report_dir / f"{sample_name}_hadv_reference_selection" / "phf_typing"
    phf_table = _read_tsv_rows(phf_dir / "phf_typing.tsv")
    phf_columns = phf_table.get("columns", [])
    phf_rows = phf_table.get("rows", [])
    subject_index = phf_columns.index("subject") if "subject" in phf_columns else -1
    type_index = phf_columns.index("matched_type") if "matched_type" in phf_columns else -1
    subject_map: dict[str, str] = {}
    type_map: dict[str, str] = {}
    for row in phf_rows:
        if not isinstance(row, list) or not row:
            continue
        gene_name = str(row[1] if len(row) > 1 else "").strip().lower()
        if not gene_name:
            continue
        if subject_index >= 0 and subject_index < len(row):
            subject_map[gene_name] = str(row[subject_index] or "").strip()
        if type_index >= 0 and type_index < len(row):
            type_map[gene_name] = str(row[type_index] or "").strip()

    rows: list[list[object]] = []
    coverage_points: list[float] = []
    depth_points: list[float] = []
    x_values: list[str] = []
    for gene_name, display_name in [("penton", "Penton"), ("hexon", "Hexon"), ("fiber", "Fiber")]:
        coverage_path = phf_dir / f"{gene_name}.coverage.tsv"
        if not coverage_path.is_file():
            continue
        matched_row: dict[str, str] | None = None
        best_row: dict[str, str] | None = None
        try:
            with coverage_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                for row in reader:
                    current = dict(row)
                    current_name = str(current.get("#rname") or current.get("rname") or "").strip()
                    current_coverage = float(current.get("coverage") or 0.0)
                    current_depth = float(current.get("meandepth") or 0.0)
                    if best_row is None:
                        best_row = current
                    else:
                        best_score = (
                            float(best_row.get("coverage") or 0.0),
                            float(best_row.get("meandepth") or 0.0),
                        )
                        current_score = (current_coverage, current_depth)
                        if current_score > best_score:
                            best_row = current
                    if subject_map.get(gene_name) and current_name == subject_map.get(gene_name):
                        matched_row = current
        except OSError:
            continue
        selected = matched_row or best_row
        if not selected:
            continue
        coverage_value = float(selected.get("coverage") or 0.0)
        mean_depth_value = float(selected.get("meandepth") or 0.0)
        reference_name = str(selected.get("#rname") or selected.get("rname") or subject_map.get(gene_name) or "").strip()
        type_label = type_map.get(gene_name) or "-"
        rows.append([
            display_name,
            type_label,
            reference_name or "-",
            f"{coverage_value:.2f}%",
            f"{mean_depth_value:.2f}",
            str(selected.get("covbases") or "-"),
            str(selected.get("endpos") or "-"),
        ])
        coverage_points.append(round(coverage_value, 2))
        depth_points.append(round(mean_depth_value, 2))
        x_values.append(display_name)
    if not rows:
        return {"status": "empty", "columns": [], "rows": []}
    return {
        "status": "ready",
        "columns": ["基因", "命中分型", "命中参考", "覆盖度", "平均深度", "覆盖碱基数", "参考长度"],
        "rows": rows,
        "coverage_points": coverage_points,
        "depth_points": depth_points,
        "x_values": x_values,
    }

def _build_serotype_section(report_dir: Path, sample_name: str, checkm_info: dict) -> dict:
    nextclade_tsv = report_dir / "nextclade_output" / "nextclade.tsv"
    if nextclade_tsv.is_file():
        raw = _read_tsv_rows(nextclade_tsv)
        columns = raw.get("columns", [])
        rows = raw.get("rows", [])
        first_row = rows[0] if rows and isinstance(rows[0], list) else []

        def _value(column_name: str, fallback: str = "-") -> str:
            if column_name not in columns:
                return fallback
            index = columns.index(column_name)
            value = str(first_row[index] if index < len(first_row) else "").strip()
            return value or fallback

        seq_name = _value("seqName", sample_name or "-")
        seq_names: list[str] = []
        if "seqName" in columns:
            seq_name_index = columns.index("seqName")
            for row in rows:
                if isinstance(row, list) and seq_name_index < len(row):
                    current_seq_name = str(row[seq_name_index] or "").strip()
                    if current_seq_name:
                        seq_names.append(current_seq_name)
        clade = _value("clade")
        clade_display = _value("clade_display")
        clade_who = _value("clade_who")
        pango = _value("Nextclade_pango")
        qc_status = _value("qc.overallStatus")
        qc_score = _value("qc.overallScore")
        total_substitutions = _value("totalSubstitutions")
        total_deletions = _value("totalDeletions")
        total_insertions = _value("totalInsertions")
        total_aa_substitutions = _value("totalAminoacidSubstitutions")
        coverage = _value("coverage")
        missing = _value("missing")
        private_nuc = _value("privateNucMutations.unlabeledSubstitutions")
        aa_substitutions = _value("aaSubstitutions")
        serotype_summary = _read_tsv_rows(report_dir / f"{sample_name}_serotype_result.tsv")
        summary_columns = serotype_summary.get("columns", [])
        summary_rows = serotype_summary.get("rows", [])

        def _summary_value(column_name: str, fallback: str = "-") -> str:
            if column_name not in summary_columns:
                return fallback
            first_summary_row = summary_rows[0] if summary_rows and isinstance(summary_rows[0], list) else []
            index = summary_columns.index(column_name)
            value = str(first_summary_row[index] if index < len(first_summary_row) else "").strip()
            return value or fallback

        virus_type = _summary_value("病毒类型")
        nextclade_dataset = _summary_value("Nextclade数据集")
        monkeypox_lineage = _value("lineage")
        monkeypox_outbreak = _value("outbreak")
        rsv_lineage = next(
            (
                candidate
                for candidate in [
                    _value("lineage", ""),
                    _value("genotype", ""),
                    _value("Nextclade_pango", ""),
                ]
                if candidate and candidate != "-"
            ),
            "-",
        )
        hmpv_lineage = next(
            (
                candidate
                for candidate in [
                    _value("lineage", ""),
                    _value("genotype", ""),
                    _value("Nextclade_pango", ""),
                ]
                if candidate and candidate != "-"
            ),
            "-",
        )
        denv_lineage = next(
            (
                candidate
                for candidate in [
                    _value("lineage", ""),
                    _value("genotype", ""),
                    _value("serotype", ""),
                    _value("Nextclade_pango", ""),
                ]
                if candidate and candidate != "-"
            ),
            "-",
        )
        zikav_lineage = next(
            (
                candidate
                for candidate in [
                    _value("lineage", ""),
                    _value("genotype", ""),
                    _value("Nextclade_pango", ""),
                ]
                if candidate and candidate != "-"
            ),
            "-",
        )
        chikv_lineage = next(
            (
                candidate
                for candidate in [
                    _value("lineage", ""),
                    _value("genotype", ""),
                    _value("Nextclade_pango", ""),
                ]
                if candidate and candidate != "-"
            ),
            "-",
        )
        ebola_lineage = next(
            (
                candidate
                for candidate in [
                    _value("lineage", ""),
                    _value("genotype", ""),
                    _value("Nextclade_pango", ""),
                ]
                if candidate and candidate != "-"
            ),
            "-",
        )
        is_hmpv_nextclade = (
            "Human metapneumovirus" in virus_type
            or "人偏肺病毒" in virus_type
            or "metapneumovirus" in virus_type.lower()
            or "/hmpv" in nextclade_dataset.lower()
            or "nextclade_db/hmpv" in nextclade_dataset.lower()
        )
        if is_hmpv_nextclade:
            notes: list[str] = []
            summary_note = _sanitize_virus_demo_note(_summary_value("说明", ""))
            mutation_table = _read_rsv_variant_annotation_table(report_dir)
            if nextclade_dataset != "-":
                notes.append(f"数据集：{nextclade_dataset}")
            if summary_note:
                notes.append(summary_note)
            return {
                "status": "ready",
                "mode": "hmpv_nextclade",
                "predicted_clade": clade,
                "predicted_lineage": hmpv_lineage,
                "summary_cards": [
                    {"label": "Nextclade Clade", "value": clade},
                    {"label": "Lineage / Genotype", "value": hmpv_lineage},
                    {"label": "QC 状态", "value": qc_status},
                    {"label": "覆盖度", "value": coverage},
                ],
                "quality_metrics": [
                    {"label": "QC 分数", "value": qc_score},
                    {"label": "核苷酸替换", "value": total_substitutions},
                    {"label": "缺失数", "value": total_deletions},
                    {"label": "插入数", "value": total_insertions},
                    {"label": "氨基酸替换", "value": total_aa_substitutions},
                    {"label": "缺失区段", "value": missing},
                ],
                "sequence_name": seq_name,
                "notes": "；".join(notes),
                "mutation_table": mutation_table,
                "igv": _discover_hmpv_igv_assets(report_dir),
                "phylogeny_tree": _build_nextclade_phylogeny_tree(report_dir, seq_names or seq_name),
                "columns": columns,
                "rows": rows[:1],
            }
        is_denv_nextclade = (
            "dengue virus" in virus_type.lower()
            or "登革热" in virus_type
            or virus_type.lower().startswith("denv")
            or "/denv1" in nextclade_dataset.lower()
            or "/denv2" in nextclade_dataset.lower()
            or "/denv3" in nextclade_dataset.lower()
            or "/denv4" in nextclade_dataset.lower()
            or "nextclade_db/denv1" in nextclade_dataset.lower()
            or "nextclade_db/denv2" in nextclade_dataset.lower()
            or "nextclade_db/denv3" in nextclade_dataset.lower()
            or "nextclade_db/denv4" in nextclade_dataset.lower()
        )
        if is_denv_nextclade:
            notes: list[str] = []
            summary_note = _sanitize_virus_demo_note(_summary_value("说明", ""))
            selection_path = report_dir / f"{sample_name}_denv_reference_selection" / "selection.tsv"
            selection_table = _read_tsv_rows(selection_path)
            selection_columns = selection_table.get("columns", [])
            selection_rows = selection_table.get("rows", [])
            denv_type_hint = ""
            dataset_hint = nextclade_dataset.lower()
            virus_type_hint = virus_type.lower()
            clade_hint = clade.strip()
            virus_type_match = re.search(r"dengue virus\s*([1-4])", virus_type_hint)
            dataset_match = re.search(r"denv([1-4])", dataset_hint)
            clade_match = re.match(r"\s*([1-4])", clade_hint)
            for match in [virus_type_match, dataset_match, clade_match]:
                if match:
                    denv_type_hint = match.group(1)
                    break

            first_selection: list[object] = []
            if selection_rows and isinstance(selection_rows[0], list):
                first_selection = selection_rows[0]
                if "denv_type" in selection_columns:
                    denv_type_index = selection_columns.index("denv_type")
                    coverage_index = selection_columns.index("coverage") if "coverage" in selection_columns else -1
                    best_row = None
                    best_score = (-1, -1.0)
                    for row in selection_rows:
                        if not isinstance(row, list):
                            continue
                        row_type = str(row[denv_type_index] if denv_type_index < len(row) else "").strip()
                        type_score = 1 if denv_type_hint and row_type == denv_type_hint else 0
                        try:
                            coverage_score = float(row[coverage_index]) if coverage_index >= 0 and coverage_index < len(row) else 0.0
                        except (TypeError, ValueError):
                            coverage_score = 0.0
                        score = (type_score, coverage_score)
                        if score > best_score:
                            best_score = score
                            best_row = row
                    if isinstance(best_row, list):
                        first_selection = best_row

            def _denv_selection_value(column_name: str, fallback: str = "-") -> str:
                if column_name not in selection_columns:
                    return fallback
                index = selection_columns.index(column_name)
                value = str(first_selection[index] if index < len(first_selection) else "").strip()
                return value or fallback

            selected_denv_type = _denv_selection_value("denv_type")
            selected_reference = Path(_denv_selection_value("reference_path", "")).name or "-"
            mutation_table = _read_rsv_variant_annotation_table(report_dir)
            if nextclade_dataset != "-":
                notes.append(f"数据集：{nextclade_dataset}")
            if selected_denv_type != "-":
                notes.append(f"自动选择参考型别：DENV{selected_denv_type}")
            if selected_reference != "-":
                notes.append(f"参考序列：{selected_reference}")
            if summary_note:
                notes.append(summary_note)
            return {
                "status": "ready",
                "mode": "denv_nextclade",
                "predicted_serotype": f"DENV-{selected_denv_type}" if selected_denv_type not in {"", "-"} else "",
                "predicted_clade": clade,
                "predicted_lineage": denv_lineage,
                "summary_cards": [
                    {"label": "Nextclade Clade", "value": clade},
                    {"label": "Lineage / Genotype", "value": denv_lineage},
                    {"label": "QC 状态", "value": qc_status},
                    {"label": "覆盖度", "value": coverage},
                ],
                "quality_metrics": [
                    {"label": "QC 分数", "value": qc_score},
                    {"label": "核苷酸替换", "value": total_substitutions},
                    {"label": "缺失数", "value": total_deletions},
                    {"label": "插入数", "value": total_insertions},
                    {"label": "氨基酸替换", "value": total_aa_substitutions},
                    {"label": "缺失区段", "value": missing},
                ],
                "sequence_name": seq_name,
                "notes": "；".join(notes),
                "mutation_table": mutation_table,
                "igv": _discover_denv_igv_assets(report_dir),
                "phylogeny_tree": _build_nextclade_phylogeny_tree(report_dir, seq_names or seq_name),
                "columns": columns,
                "rows": rows[:1],
            }
        is_zikav_nextclade = (
            "zika virus" in virus_type.lower()
            or "寨卡" in virus_type
            or virus_type.lower().startswith("zik")
            or "/zikav" in nextclade_dataset.lower()
            or "nextclade_db/zikav" in nextclade_dataset.lower()
        )
        if is_zikav_nextclade:
            notes: list[str] = []
            summary_note = _sanitize_virus_demo_note(_summary_value("说明", ""))
            mutation_table = _read_rsv_variant_annotation_table(report_dir)
            if nextclade_dataset != "-":
                notes.append(f"数据集：{nextclade_dataset}")
            if summary_note:
                notes.append(summary_note)
            return {
                "status": "ready",
                "mode": "zikav_nextclade",
                "predicted_clade": clade,
                "predicted_lineage": zikav_lineage,
                "summary_cards": [
                    {"label": "Nextclade Clade", "value": clade},
                    {"label": "Lineage / Genotype", "value": zikav_lineage},
                    {"label": "QC 状态", "value": qc_status},
                    {"label": "覆盖度", "value": coverage},
                ],
                "quality_metrics": [
                    {"label": "QC 分数", "value": qc_score},
                    {"label": "核苷酸替换", "value": total_substitutions},
                    {"label": "缺失数", "value": total_deletions},
                    {"label": "插入数", "value": total_insertions},
                    {"label": "氨基酸替换", "value": total_aa_substitutions},
                    {"label": "缺失区段", "value": missing},
                ],
                "sequence_name": seq_name,
                "notes": "；".join(notes),
                "mutation_table": mutation_table,
                "igv": _discover_zikav_igv_assets(report_dir),
                "phylogeny_tree": _build_nextclade_phylogeny_tree(report_dir, seq_names or seq_name),
                "columns": columns,
                "rows": rows[:1],
            }
        is_chikv_nextclade = (
            "chikungunya virus" in virus_type.lower()
            or "chikungunya" in virus_type.lower()
            or "基孔肯雅" in virus_type
            or virus_type.lower().startswith("chikv")
            or "/chikv" in nextclade_dataset.lower()
            or "nextclade_db/chikv" in nextclade_dataset.lower()
        )
        if is_chikv_nextclade:
            notes: list[str] = []
            summary_note = _sanitize_virus_demo_note(_summary_value("说明", ""))
            mutation_table = _read_rsv_variant_annotation_table(report_dir)
            selected_reference = "reference.fasta" if (report_dir / "genomes" / "ref.fa").is_file() else "-"
            display_columns = ["样本名称", "病毒类型", "Clade", "Lineage / Genotype", "参考序列", "QC 状态", "覆盖度", "说明"]
            display_rows = [[
                sample_name or seq_name or "-",
                "Chikungunya virus",
                clade,
                chikv_lineage,
                selected_reference,
                qc_status,
                coverage,
                summary_note or "基于 CHIKV 固定参考与 Nextclade 数据集完成 clade / lineage 判读，并结合突变位点表和 IGV 比对结果组织展示。",
            ]]
            if nextclade_dataset != "-":
                notes.append(f"数据集：{nextclade_dataset}")
            if selected_reference != "-":
                notes.append(f"参考序列：{selected_reference}")
            if summary_note:
                notes.append(summary_note)
            return {
                "status": "ready",
                "mode": "chikv_nextclade",
                "predicted_clade": clade,
                "predicted_lineage": chikv_lineage,
                "reference_name": selected_reference,
                "summary_cards": [
                    {"label": "Nextclade Clade", "value": clade},
                    {"label": "Lineage / Genotype", "value": chikv_lineage},
                    {"label": "参考序列", "value": selected_reference},
                    {"label": "QC 状态", "value": qc_status},
                    {"label": "覆盖度", "value": coverage},
                ],
                "quality_metrics": [
                    {"label": "QC 分数", "value": qc_score},
                    {"label": "核苷酸替换", "value": total_substitutions},
                    {"label": "缺失数", "value": total_deletions},
                    {"label": "插入数", "value": total_insertions},
                    {"label": "氨基酸替换", "value": total_aa_substitutions},
                    {"label": "缺失区段", "value": missing},
                ],
                "sequence_name": seq_name,
                "notes": "；".join(notes),
                "mutation_table": mutation_table,
                "igv": _discover_chikv_igv_assets(report_dir),
                "phylogeny_tree": _build_nextclade_phylogeny_tree(report_dir, seq_names or seq_name),
                "columns": display_columns,
                "rows": display_rows,
            }
        is_ebola_nextclade = (
            "ebola virus" in virus_type.lower()
            or "ebolavirus" in virus_type.lower()
            or "orthoebolavirus" in virus_type.lower()
            or "埃博拉" in virus_type
            or "/ebola" in nextclade_dataset.lower()
            or "nextclade_db/ebola" in nextclade_dataset.lower()
        )
        if is_ebola_nextclade:
            notes: list[str] = []
            summary_note = _sanitize_virus_demo_note(_summary_value("说明", ""))
            selection_path = report_dir / f"{sample_name}_orthoebolavirus_reference_selection" / "selection.tsv"
            selection_table = _read_tsv_rows(selection_path)
            selection_columns = selection_table.get("columns", [])
            selection_rows = selection_table.get("rows", [])
            first_selection = selection_rows[0] if selection_rows and isinstance(selection_rows[0], list) else []

            def _ebola_selection_value(column_name: str, fallback: str = "-") -> str:
                if column_name not in selection_columns:
                    return fallback
                index = selection_columns.index(column_name)
                value = str(first_selection[index] if index < len(first_selection) else "").strip()
                return value or fallback

            selected_type = _ebola_selection_value("predicted_type")
            selected_reference = Path(_ebola_selection_value("reference_path", "")).name or "-"
            mutation_table = _read_rsv_variant_annotation_table(report_dir)
            display_columns = ["样本名称", "病毒类型", "本地参考分型", "Nextclade Clade", "Lineage / Genotype", "参考序列", "QC 状态", "覆盖度", "说明"]
            display_rows = [[
                sample_name or seq_name or "-",
                "Ebola virus",
                selected_type,
                clade,
                ebola_lineage,
                selected_reference,
                qc_status,
                coverage,
                summary_note or "基于 Orthoebolavirus 本地参考库完成最优参考筛选；EBOV 样本继续进入 Ebola Nextclade 数据集完成 clade / lineage 判读。",
            ]]
            if nextclade_dataset != "-":
                notes.append(f"数据集：{nextclade_dataset}")
            if selected_type != "-":
                notes.append(f"本地参考筛选：{selected_type}")
            if selected_reference != "-":
                notes.append(f"参考序列：{selected_reference}")
            if summary_note:
                notes.append(summary_note)
            return {
                "status": "ready",
                "mode": "ebola_nextclade",
                "predicted_serotype": selected_type if selected_type != "-" else "",
                "predicted_clade": clade,
                "predicted_lineage": ebola_lineage,
                "reference_name": selected_reference,
                "summary_cards": [
                    {"label": "Nextclade Clade", "value": clade},
                    {"label": "Lineage / Genotype", "value": ebola_lineage},
                    {"label": "本地参考分型", "value": selected_type},
                    {"label": "QC 状态", "value": qc_status},
                    {"label": "覆盖度", "value": coverage},
                ],
                "quality_metrics": [
                    {"label": "QC 分数", "value": qc_score},
                    {"label": "核苷酸替换", "value": total_substitutions},
                    {"label": "缺失数", "value": total_deletions},
                    {"label": "插入数", "value": total_insertions},
                    {"label": "氨基酸替换", "value": total_aa_substitutions},
                    {"label": "缺失区段", "value": missing},
                ],
                "sequence_name": seq_name,
                "notes": "；".join(notes),
                "mutation_table": mutation_table,
                "igv": _discover_ebola_igv_assets(report_dir),
                "phylogeny_tree": _build_nextclade_phylogeny_tree(report_dir, seq_names or seq_name),
                "orthoebolavirus_selection": selection_table,
                "nextclade_columns": columns,
                "nextclade_rows": rows[:1],
                "columns": display_columns,
                "rows": display_rows,
            }
        is_rsv_nextclade = (
            "Respiratory syncytial" in virus_type
            or "合胞病毒" in virus_type
            or "orthopneumovirus" in virus_type.lower()
            or "/rsv_" in nextclade_dataset.lower()
            or "nextclade_db/rsv_" in nextclade_dataset.lower()
        )
        if is_rsv_nextclade:
            notes: list[str] = []
            summary_note = _sanitize_virus_demo_note(_summary_value("说明", ""))
            selection_path = report_dir / f"{sample_name}_rsv_reference_selection" / "selection.tsv"
            selection_table = _read_tsv_rows(selection_path)
            selection_columns = selection_table.get("columns", [])
            selection_rows = selection_table.get("rows", [])
            rsv_type_hint = ""
            dataset_hint = nextclade_dataset.lower()
            virus_type_hint = virus_type.lower()
            clade_hint = clade.strip().upper()
            if "rsv_b" in dataset_hint or "respiratory syncytial virus b" in virus_type_hint or clade_hint.startswith("B"):
                rsv_type_hint = "B"
            elif "rsv_a" in dataset_hint or "respiratory syncytial virus a" in virus_type_hint or clade_hint.startswith("A"):
                rsv_type_hint = "A"

            first_selection: list[object] = []
            if selection_rows and isinstance(selection_rows[0], list):
                first_selection = selection_rows[0]
                if "rsv_type" in selection_columns:
                    rsv_type_index = selection_columns.index("rsv_type")
                    coverage_index = selection_columns.index("coverage") if "coverage" in selection_columns else -1
                    best_row = None
                    best_score = (-1, -1.0)
                    for row in selection_rows:
                        if not isinstance(row, list):
                            continue
                        row_type = str(row[rsv_type_index] if rsv_type_index < len(row) else "").strip()
                        type_score = 1 if rsv_type_hint and row_type == rsv_type_hint else 0
                        try:
                            coverage_score = float(row[coverage_index]) if coverage_index >= 0 and coverage_index < len(row) else 0.0
                        except (TypeError, ValueError):
                            coverage_score = 0.0
                        score = (type_score, coverage_score)
                        if score > best_score:
                            best_score = score
                            best_row = row
                    if isinstance(best_row, list):
                        first_selection = best_row

            def _selection_value(column_name: str, fallback: str = "-") -> str:
                if column_name not in selection_columns:
                    return fallback
                index = selection_columns.index(column_name)
                value = str(first_selection[index] if index < len(first_selection) else "").strip()
                return value or fallback

            selected_rsv_type = _selection_value("rsv_type")
            selected_reference = Path(_selection_value("reference_path", "")).name or "-"
            mutation_table = _read_rsv_variant_annotation_table(report_dir)
            nmdc_annotation = _read_rsv_nmdc_annotation_table(report_dir, mutation_table, selected_rsv_type or virus_type)
            if nextclade_dataset != "-":
                notes.append(f"数据集：{nextclade_dataset}")
            if selected_rsv_type != "-":
                notes.append(f"自动选择参考型别：RSV {selected_rsv_type}")
            if selected_reference != "-":
                notes.append(f"参考序列：{selected_reference}")
            if summary_note:
                notes.append(summary_note)
            return {
                "status": "ready",
                "mode": "rsv_nextclade",
                "predicted_serotype": f"RSV-{selected_rsv_type}" if selected_rsv_type not in {"", "-"} else "",
                "predicted_clade": clade,
                "predicted_lineage": rsv_lineage,
                "summary_cards": [
                    {"label": "Nextclade Clade", "value": clade},
                    {"label": "Lineage / Genotype", "value": rsv_lineage},
                    {"label": "QC 状态", "value": qc_status},
                    {"label": "覆盖度", "value": coverage},
                ],
                "quality_metrics": [
                    {"label": "QC 分数", "value": qc_score},
                    {"label": "核苷酸替换", "value": total_substitutions},
                    {"label": "缺失数", "value": total_deletions},
                    {"label": "插入数", "value": total_insertions},
                    {"label": "氨基酸替换", "value": total_aa_substitutions},
                    {"label": "缺失区段", "value": missing},
                ],
                "sequence_name": seq_name,
                "notes": "；".join(notes),
                "mutation_table": mutation_table,
                "nmdc_annotation": nmdc_annotation,
                "igv": _discover_rsv_igv_assets(report_dir),
                "phylogeny_tree": _build_nextclade_phylogeny_tree(report_dir, seq_names or seq_name),
                "columns": columns,
                "rows": rows[:1],
            }
        is_monkeypox_nextclade = (
            "Monkeypox" in virus_type
            or "猴痘" in virus_type
            or ("lineage" in columns and "outbreak" in columns and "Nextclade_pango" not in columns)
        )
        if is_monkeypox_nextclade:
            notes: list[str] = []
            dataset_path = _summary_value("Nextclade数据集")
            summary_note = _sanitize_virus_demo_note(_summary_value("说明", ""))
            mutation_table = _read_monkeypox_variant_annotation_table(report_dir)
            if dataset_path != "-":
                notes.append(f"数据集：{dataset_path}")
            if summary_note:
                notes.append(summary_note)
            return {
                "status": "ready",
                "mode": "monkeypox_nextclade",
                "predicted_clade": clade,
                "predicted_lineage": monkeypox_lineage,
                "predicted_outbreak": monkeypox_outbreak,
                "summary_cards": [
                    {"label": "Nextclade Clade", "value": clade},
                    {"label": "Lineage", "value": monkeypox_lineage},
                    {"label": "Outbreak", "value": monkeypox_outbreak},
                    {"label": "QC 状态", "value": qc_status},
                ],
                "quality_metrics": [
                    {"label": "QC 分数", "value": qc_score},
                    {"label": "覆盖度", "value": coverage},
                    {"label": "核苷酸替换", "value": total_substitutions},
                    {"label": "缺失数", "value": total_deletions},
                    {"label": "插入数", "value": total_insertions},
                    {"label": "氨基酸替换", "value": total_aa_substitutions},
                    {"label": "移码数", "value": _value("totalFrameShifts")},
                    {"label": "缺失区段", "value": missing},
                ],
                "sequence_name": seq_name,
                "notes": "；".join(notes),
                "mutation_table": mutation_table,
                "igv": _discover_monkeypox_igv_assets(report_dir),
                "phylogeny_tree": _build_nextclade_phylogeny_tree(report_dir, seq_names or seq_name),
                "columns": columns,
                "rows": rows[:1],
            }

        raw_ncov_variant_annotation = _read_ncov_variant_annotation_table(report_dir)
        ncov_variant_annotation = _build_ncov_mutation_display_table(report_dir, raw_ncov_variant_annotation)
        ncov_igv_view = _discover_ncov_igv_assets(report_dir)
        ngdc_mutation_knowledge = _collect_ncov_ngdc_matches(columns, first_row)
        notes: list[str] = []
        if clade_who != "-":
            notes.append(f"WHO 命名：{clade_who}")
        if clade_display != "-" and clade_display != clade:
            notes.append(f"显示分型：{clade_display}")
        note_text = "；".join(notes)
        mutation_preview = [item.strip() for item in private_nuc.split(",") if item.strip()][:12]
        aa_preview = [item.strip() for item in aa_substitutions.split(",") if item.strip()][:20]
        return {
            "status": "ready",
            "mode": "sars_cov_2_nextclade",
            "predicted_clade": clade,
            "pango_lineage": pango,
            "summary_cards": [
                {"label": "Nextclade Clade", "value": clade},
                {"label": "Pango 谱系", "value": pango},
                {"label": "QC 状态", "value": qc_status},
                {"label": "覆盖度", "value": coverage},
            ],
            "quality_metrics": [
                {"label": "QC 分数", "value": qc_score},
                {"label": "核苷酸替换", "value": total_substitutions},
                {"label": "缺失数", "value": total_deletions},
                {"label": "插入数", "value": total_insertions},
                {"label": "氨基酸替换", "value": total_aa_substitutions},
                {"label": "缺失区段", "value": missing},
            ],
            "sequence_name": seq_name,
            "notes": note_text,
            "mutation_preview": mutation_preview,
            "aa_mutation_preview": aa_preview,
            "mutation_knowledge": ngdc_mutation_knowledge,
            "variant_annotation": ncov_variant_annotation,
            "igv": ncov_igv_view,
            "phylogeny_tree": _build_nextclade_phylogeny_tree(report_dir, seq_names or seq_name),
            "columns": columns,
            "rows": rows[:1],
        }

    influenza_typing = _read_influenza_typing_section(report_dir, sample_name)
    if influenza_typing:
        return influenza_typing

    hadv_summary = _read_tsv_rows(report_dir / f"{sample_name}_serotype_result.tsv")
    hadv_columns = hadv_summary.get("columns", [])
    hadv_rows = hadv_summary.get("rows", [])
    if any(str(column).strip() == "HAdV分型" for column in hadv_columns):
        first_summary_row = hadv_rows[0] if hadv_rows and isinstance(hadv_rows[0], list) else []

        def _hadv_summary_value(column_name: str, fallback: str = "-") -> str:
            if column_name not in hadv_columns:
                return fallback
            index = hadv_columns.index(column_name)
            value = str(first_summary_row[index] if index < len(first_summary_row) else "").strip()
            return value or fallback

        selection_path = report_dir / f"{sample_name}_hadv_reference_selection" / "selection.tsv"
        phf_path = report_dir / f"{sample_name}_hadv_reference_selection" / "phf_typing" / "phf_typing.tsv"
        selection_table = _read_tsv_rows(selection_path)
        phf_table = _read_tsv_rows(phf_path)
        selection_columns = selection_table.get("columns", [])
        selection_rows = selection_table.get("rows", [])
        first_selection = selection_rows[0] if selection_rows and isinstance(selection_rows[0], list) else []

        def _hadv_selection_value(column_name: str, fallback: str = "-") -> str:
            if column_name not in selection_columns:
                return fallback
            index = selection_columns.index(column_name)
            value = str(first_selection[index] if index < len(first_selection) else "").strip()
            return value or fallback

        selected_hadv_type = _hadv_summary_value("HAdV分型", _hadv_selection_value("hadv_type"))
        penton_type = _hadv_summary_value("Penton分型", _hadv_selection_value("penton_type"))
        hexon_type = _hadv_summary_value("Hexon分型", _hadv_selection_value("hexon_type"))
        fiber_type = _hadv_summary_value("Fiber分型", _hadv_selection_value("fiber_type"))
        selected_reference = _hadv_summary_value("参考序列", Path(_hadv_selection_value("reference_path", "")).name or _hadv_selection_value("reference_name"))
        selected_gff = _hadv_summary_value("注释文件", Path(_hadv_selection_value("gff_path", "")).name or "-")
        coverage_value = _hadv_summary_value("覆盖度", _hadv_selection_value("coverage"))
        mean_depth_value = _hadv_summary_value("平均深度", _hadv_selection_value("mean_depth"))
        summary_note = _sanitize_virus_demo_note(_hadv_summary_value("说明", ""))
        notes: list[str] = []
        if selected_hadv_type != "-":
            notes.append(f"总分型：{selected_hadv_type}")
        if penton_type != "-" or hexon_type != "-" or fiber_type != "-":
            notes.append(f"PHF 分型：{penton_type}/{hexon_type}/{fiber_type}")
        if selected_reference != "-":
            notes.append(f"参考序列：{selected_reference}")
        if selected_gff != "-":
            notes.append(f"注释文件：{selected_gff}")
        if summary_note and summary_note != "-":
            notes.append(summary_note)
        mutation_table = _read_rsv_variant_annotation_table(report_dir)
        phf_coverage = _read_hadv_phf_coverage_section(report_dir, sample_name)
        phf_snp = _read_hadv_phf_snp_section(report_dir, sample_name, mutation_table)
        display_columns = hadv_columns
        display_rows = hadv_rows[:1]
        if not display_rows and selection_rows:
            display_columns = [
                column
                for column in selection_columns
                if str(column).strip() not in {"reference_path", "gff_path", "phf_summary_path"}
            ]
            selected_indexes = [selection_columns.index(column) for column in display_columns]
            display_rows = [
                [row[index] if index < len(row) else "" for index in selected_indexes]
                for row in selection_rows
                if isinstance(row, list)
            ]
        return {
            "status": "ready",
            "mode": "hadv_typing",
            "predicted_clade": selected_hadv_type,
            "reference_name": selected_reference,
            "summary_cards": [
                {"label": "HAdV 分型", "value": selected_hadv_type},
                {"label": "Penton", "value": penton_type},
                {"label": "Hexon", "value": hexon_type},
                {"label": "Fiber", "value": fiber_type},
                {"label": "参考序列", "value": selected_reference},
                {"label": "覆盖度", "value": coverage_value},
                *(
                    [
                        {"label": f"{str(item[0])} SNP", "value": int(item[3])}
                        for item in phf_snp.get("rows", [])
                        if isinstance(item, list) and len(item) > 3
                    ]
                ),
            ],
            "quality_metrics": [
                {"label": "平均深度", "value": mean_depth_value},
            ] if mean_depth_value != "-" else [],
            "sequence_name": sample_name or "-",
            "notes": "；".join(notes),
            "mutation_table": mutation_table,
            "phf_table": phf_table,
            "phf_coverage": phf_coverage,
            "phf_snp": phf_snp,
            "igv": _discover_hadv_igv_assets(report_dir),
            "columns": display_columns,
            "rows": display_rows,
        }

    norovirus_summary = _read_tsv_rows(report_dir / f"{sample_name}_serotype_result.tsv")
    norovirus_columns = norovirus_summary.get("columns", [])
    norovirus_rows = norovirus_summary.get("rows", [])
    norovirus_selection = _read_tsv_rows(report_dir / f"{sample_name}_norovirus_reference_selection" / "selection.tsv")
    norovirus_dual_typing = _read_tsv_rows(report_dir / f"{sample_name}_norovirus_reference_selection" / "typing" / "dual_typing.tsv")
    norovirus_consensus_typing = _read_tsv_rows(report_dir / f"{sample_name}_norovirus_reference_selection" / "consensus_typing" / "consensus_typing.tsv")
    has_norovirus_assets = (
        any(str(column).strip() in {"双位点分型", "RdRp分型"} for column in norovirus_columns)
        or bool(norovirus_selection.get("rows"))
        or bool(norovirus_dual_typing.get("rows"))
    )
    if has_norovirus_assets:
        selection_columns = norovirus_selection.get("columns", [])
        selection_rows = norovirus_selection.get("rows", [])
        first_selection = selection_rows[0] if selection_rows and isinstance(selection_rows[0], list) else []

        def _norovirus_selection_value(column_name: str, fallback: str = "-") -> str:
            if column_name not in selection_columns:
                return fallback
            index = selection_columns.index(column_name)
            value = str(first_selection[index] if index < len(first_selection) else "").strip()
            return value or fallback

        dual_columns = norovirus_dual_typing.get("columns", [])
        dual_rows = [
            row for row in norovirus_dual_typing.get("rows", [])
            if isinstance(row, list)
        ]

        def _norovirus_dual_value(gene_name: str, column_name: str, fallback: str = "-") -> str:
            if not dual_rows or column_name not in dual_columns or "gene" not in dual_columns:
                return fallback
            gene_index = dual_columns.index("gene")
            column_index = dual_columns.index(column_name)
            for row in dual_rows:
                gene_value = str(row[gene_index] if gene_index < len(row) else "").strip().lower()
                if gene_value != gene_name.lower():
                    continue
                value = str(row[column_index] if column_index < len(row) else "").strip()
                if value:
                    return value
            return fallback

        consensus_columns = norovirus_consensus_typing.get("columns", [])
        consensus_rows = norovirus_consensus_typing.get("rows", [])
        first_consensus = consensus_rows[0] if consensus_rows and isinstance(consensus_rows[0], list) else []

        def _norovirus_consensus_value(column_name: str, fallback: str = "-") -> str:
            if column_name not in consensus_columns:
                return fallback
            index = consensus_columns.index(column_name)
            value = str(first_consensus[index] if index < len(first_consensus) else "").strip()
            return value or fallback

        dual_type = _norovirus_selection_value("dual_type", _norovirus_consensus_value("dual_type", "-"))
        rdrp_type = _norovirus_selection_value("rdrp_type", _norovirus_dual_value("rdrp", "matched_type", _norovirus_consensus_value("rdrp_type", "-")))
        vp1_type = _norovirus_selection_value("vp1_type", _norovirus_dual_value("vp1", "matched_type", _norovirus_consensus_value("vp1_type", "-")))
        selected_reference = _norovirus_selection_value("reference_name", Path(_norovirus_selection_value("reference_path", "")).name or "-")
        selected_gff = Path(_norovirus_selection_value("gff_path", "")).name if _norovirus_selection_value("gff_path", "") not in {"", "-", "nogtf"} else "-"
        coverage_value = _norovirus_selection_value("coverage", "-")
        mean_depth_value = _norovirus_selection_value("mean_depth", "-")
        covered_bases = _norovirus_selection_value("covered_bases", "-")
        num_reads = _norovirus_selection_value("num_reads", "-")
        if coverage_value not in {"-", ""}:
            try:
                coverage_value = f"{float(coverage_value):.2f}%"
            except ValueError:
                pass
        if mean_depth_value not in {"-", ""}:
            try:
                mean_depth_value = f"{float(mean_depth_value):.2f}"
            except ValueError:
                pass
        summary_note = ""
        if norovirus_columns and norovirus_rows:
            first_summary_row = norovirus_rows[0] if isinstance(norovirus_rows[0], list) else []
            if "说明" in norovirus_columns:
                note_index = norovirus_columns.index("说明")
                summary_note = _sanitize_virus_demo_note(first_summary_row[note_index] if note_index < len(first_summary_row) else "")
        notes: list[str] = []
        if dual_type != "-":
            notes.append(f"双位点分型：{dual_type}")
        if rdrp_type != "-" or vp1_type != "-":
            notes.append(f"RdRp/VP1：{rdrp_type}/{vp1_type}")
        if selected_reference != "-":
            notes.append(f"参考序列：{selected_reference}")
        if selected_gff != "-":
            notes.append(f"注释文件：{selected_gff}")
        if summary_note and summary_note != "-":
            notes.append(summary_note)
        mutation_table = _read_rsv_variant_annotation_table(report_dir)
        gene_phylogeny = _build_norovirus_gene_phylogeny(report_dir)
        display_columns = ["样本名称", "病毒类型", "双位点分型", "RdRp分型", "VP1分型", "参考序列", "覆盖度", "平均深度", "说明"]
        display_rows = [[
            sample_name or "-",
            f"Norovirus {dual_type}" if dual_type != "-" else "Norovirus",
            dual_type,
            rdrp_type,
            vp1_type,
            selected_reference,
            coverage_value,
            mean_depth_value,
            summary_note or "基于 CDC RdRp/VP1 双位点分型后选择最优全基因组参考。",
        ]]
        return {
            "status": "ready",
            "mode": "norovirus_typing",
            "predicted_clade": dual_type,
            "reference_name": selected_reference,
            "summary_cards": [
                {"label": "双位点分型", "value": dual_type},
                {"label": "RdRp 分型", "value": rdrp_type},
                {"label": "VP1 分型", "value": vp1_type},
                {"label": "参考序列", "value": selected_reference},
                {"label": "覆盖度", "value": coverage_value},
            ],
            "quality_metrics": [
                {"label": "平均深度", "value": mean_depth_value},
                {"label": "覆盖碱基", "value": covered_bases},
                {"label": "支持 reads", "value": num_reads},
            ],
            "sequence_name": sample_name or "-",
            "notes": "；".join([item for item in notes if item]),
            "mutation_table": mutation_table,
            "typing_table": norovirus_dual_typing,
            "consensus_typing": norovirus_consensus_typing,
            "igv": _discover_norovirus_igv_assets(report_dir),
            "gene_phylogeny": gene_phylogeny,
            "columns": display_columns,
            "rows": display_rows,
        }

    enterovirus_summary = _read_tsv_rows(report_dir / f"{sample_name}_serotype_result.tsv")
    enterovirus_columns = enterovirus_summary.get("columns", [])
    enterovirus_rows = enterovirus_summary.get("rows", [])
    enterovirus_selection = _read_tsv_rows(report_dir / f"{sample_name}_enterovirus_reference_selection" / "selection.tsv")
    enterovirus_consensus_typing = _read_tsv_rows(report_dir / f"{sample_name}_enterovirus_reference_selection" / "consensus_typing" / "consensus_typing.tsv")
    enterovirus_phylogeny = _build_enterovirus_gene_phylogeny(report_dir)
    has_enterovirus_assets = (
        any(str(column).strip() in {"VP1分型", "VP1 分型"} for column in enterovirus_columns)
        or bool(enterovirus_selection.get("rows"))
        or bool(enterovirus_consensus_typing.get("rows"))
        or str(enterovirus_phylogeny.get("status") or "").strip().lower() == "ready"
    )
    if has_enterovirus_assets:
        selection_columns = enterovirus_selection.get("columns", [])
        selection_rows = enterovirus_selection.get("rows", [])
        first_selection = selection_rows[0] if selection_rows and isinstance(selection_rows[0], list) else []

        def _enterovirus_selection_value(column_name: str, fallback: str = "-") -> str:
            if column_name not in selection_columns:
                return fallback
            index = selection_columns.index(column_name)
            value = str(first_selection[index] if index < len(first_selection) else "").strip()
            return value or fallback

        consensus_columns = enterovirus_consensus_typing.get("columns", [])
        consensus_rows = enterovirus_consensus_typing.get("rows", [])
        first_consensus = consensus_rows[0] if consensus_rows and isinstance(consensus_rows[0], list) else []

        def _enterovirus_consensus_value(column_name: str, fallback: str = "-") -> str:
            if column_name not in consensus_columns:
                return fallback
            index = consensus_columns.index(column_name)
            value = str(first_consensus[index] if index < len(first_consensus) else "").strip()
            return value or fallback

        vp1_type = _enterovirus_selection_value("vp1_type", _enterovirus_consensus_value("vp1_type", "-"))
        big_group = _enterovirus_selection_value("big_group", _enterovirus_consensus_value("big_group", "-")).upper()
        selected_reference = _enterovirus_selection_value("reference_name", Path(_enterovirus_selection_value("reference_path", "")).name or "-")
        selected_gff = Path(_enterovirus_selection_value("gff_path", "")).name if _enterovirus_selection_value("gff_path", "") not in {"", "-", "nogtf"} else "-"
        coverage_value = _enterovirus_selection_value("coverage", "-")
        mean_depth_value = _enterovirus_selection_value("mean_depth", "-")
        covered_bases = _enterovirus_selection_value("covered_bases", "-")
        num_reads = _enterovirus_selection_value("num_reads", "-")
        anchor_accession = _enterovirus_selection_value("anchor_accession", "-")
        candidate_count = _enterovirus_selection_value("candidate_count", "-")
        dedup_candidate_count = _enterovirus_selection_value("dedup_candidate_count", "-")
        if coverage_value not in {"-", ""}:
            try:
                coverage_value = f"{float(coverage_value):.2f}%"
            except ValueError:
                pass
        if mean_depth_value not in {"-", ""}:
            try:
                mean_depth_value = f"{float(mean_depth_value):.2f}"
            except ValueError:
                pass
        summary_note = ""
        if enterovirus_columns and enterovirus_rows:
            first_summary_row = enterovirus_rows[0] if isinstance(enterovirus_rows[0], list) else []
            if "说明" in enterovirus_columns:
                note_index = enterovirus_columns.index("说明")
                summary_note = _sanitize_virus_demo_note(first_summary_row[note_index] if note_index < len(first_summary_row) else "")
        notes: list[str] = []
        if vp1_type != "-":
            notes.append(f"VP1 分型：{vp1_type}")
        if big_group != "-":
            notes.append(f"大亚型：{big_group}")
        if selected_reference != "-":
            notes.append(f"参考序列：{selected_reference}")
        if anchor_accession not in {"", "-"}:
            notes.append(f"批次锚点参考：{anchor_accession}")
        if candidate_count not in {"", "-"} and dedup_candidate_count not in {"", "-"}:
            notes.append(f"候选基因组：{candidate_count} 条，95% 去冗余后 {dedup_candidate_count} 条")
        if summary_note and summary_note != "-":
            notes.append(summary_note)
        mutation_table = _read_rsv_variant_annotation_table(report_dir)
        display_columns = ["样本名称", "病毒类型", "大亚型", "VP1分型", "参考序列", "覆盖度", "平均深度", "说明"]
        display_rows = [[
            sample_name or "-",
            f"Human enterovirus {big_group}" if big_group != "-" else "Human enterovirus",
            big_group,
            vp1_type,
            selected_reference,
            coverage_value,
            mean_depth_value,
            summary_note or "基于 EV-A/B/C/D 的 VP1 分型结果，在对应亚型全基因组候选集中选择覆盖度最优参考；组装后再按大亚型使用 VADR 生成 GFF。",
        ]]
        return {
            "status": "ready",
            "mode": "enterovirus_typing",
            "predicted_clade": vp1_type,
            "predicted_group": big_group,
            "reference_name": selected_reference,
            "summary_cards": [
                {"label": "VP1 分型", "value": vp1_type},
                {"label": "大亚型", "value": big_group},
                {"label": "参考序列", "value": selected_reference},
                {"label": "覆盖度", "value": coverage_value},
            ],
            "quality_metrics": [
                {"label": "平均深度", "value": mean_depth_value},
                {"label": "覆盖碱基", "value": covered_bases},
                {"label": "支持 reads", "value": num_reads},
                {"label": "去冗余候选数", "value": dedup_candidate_count},
            ],
            "sequence_name": sample_name or "-",
            "notes": "；".join([item for item in notes if item]),
            "mutation_table": mutation_table,
            "consensus_typing": enterovirus_consensus_typing,
            "igv": _discover_enterovirus_igv_assets(report_dir),
            "gene_phylogeny": enterovirus_phylogeny,
            "columns": display_columns,
            "rows": display_rows,
        }

    hepatovirus_summary = _read_tsv_rows(report_dir / f"{sample_name}_serotype_result.tsv")
    hepatovirus_columns = hepatovirus_summary.get("columns", [])
    hepatovirus_rows = hepatovirus_summary.get("rows", [])
    hepatovirus_selection = _read_tsv_rows(report_dir / f"{sample_name}_hepatovirus_reference_selection" / "selection.tsv")
    hepatovirus_consensus_typing = _read_tsv_rows(report_dir / f"{sample_name}_hepatovirus_reference_selection" / "consensus_typing" / "consensus_typing.tsv")
    hepatovirus_first_row = hepatovirus_rows[0] if hepatovirus_rows and isinstance(hepatovirus_rows[0], list) else []
    hepatovirus_virus_type = ""
    if "病毒类型" in hepatovirus_columns:
        virus_type_index = hepatovirus_columns.index("病毒类型")
        hepatovirus_virus_type = str(hepatovirus_first_row[virus_type_index] if virus_type_index < len(hepatovirus_first_row) else "").strip()
    has_hiv_resistance_markers = any(
        str(column).strip() in {"NRTI最高等级", "NNRTI最高等级", "PI最高等级", "INSTI最高等级"}
        for column in hepatovirus_columns
    )
    has_hepatovirus_assets = (
        any(str(column).strip() in {"大亚型", "子亚型", "HAV子亚型", "HAV 子亚型"} for column in hepatovirus_columns)
        or bool(hepatovirus_selection.get("rows"))
        or bool(hepatovirus_consensus_typing.get("rows"))
    ) and not (
        has_hiv_resistance_markers
        or "hiv" in hepatovirus_virus_type.lower()
        or "immunodeficiency" in hepatovirus_virus_type.lower()
    )
    if has_hepatovirus_assets:
        selection_columns = hepatovirus_selection.get("columns", [])
        selection_rows = hepatovirus_selection.get("rows", [])
        first_selection = selection_rows[0] if selection_rows and isinstance(selection_rows[0], list) else []

        def _hepatovirus_selection_value(column_name: str, fallback: str = "-") -> str:
            if column_name not in selection_columns:
                return fallback
            index = selection_columns.index(column_name)
            value = str(first_selection[index] if index < len(first_selection) else "").strip()
            return value or fallback

        consensus_columns = hepatovirus_consensus_typing.get("columns", [])
        consensus_rows = hepatovirus_consensus_typing.get("rows", [])
        first_consensus = consensus_rows[0] if consensus_rows and isinstance(consensus_rows[0], list) else []

        def _hepatovirus_consensus_value(column_name: str, fallback: str = "-") -> str:
            if column_name not in consensus_columns:
                return fallback
            index = consensus_columns.index(column_name)
            value = str(first_consensus[index] if index < len(first_consensus) else "").strip()
            return value or fallback

        broad_type = _hepatovirus_selection_value("broad_type", _hepatovirus_consensus_value("broad_type", "-")).upper()
        subtype = _hepatovirus_selection_value("subtype", _hepatovirus_consensus_value("subtype", "-")).upper()
        hav_subtype = _hepatovirus_selection_value("hav_subtype", _hepatovirus_consensus_value("hav_subtype", "-")).upper()
        if subtype == "-" and hav_subtype != "-":
            subtype = hav_subtype
        selected_reference = _hepatovirus_selection_value("reference_name", Path(_hepatovirus_selection_value("reference_path", "")).name or "-")
        selected_gff = Path(_hepatovirus_selection_value("gff_path", "")).name if _hepatovirus_selection_value("gff_path", "") not in {"", "-", "nogtf"} else "-"
        coverage_value = _hepatovirus_selection_value("coverage", "-")
        mean_depth_value = _hepatovirus_selection_value("mean_depth", "-")
        covered_bases = _hepatovirus_selection_value("covered_bases", "-")
        num_reads = _hepatovirus_selection_value("num_reads", "-")
        isolate_label = _hepatovirus_selection_value("isolate", "-")
        if coverage_value not in {"-", ""}:
            try:
                coverage_value = f"{float(coverage_value):.2f}%"
            except ValueError:
                pass
        if mean_depth_value not in {"-", ""}:
            try:
                mean_depth_value = f"{float(mean_depth_value):.2f}"
            except ValueError:
                pass
        summary_note = ""
        if hepatovirus_columns and hepatovirus_rows:
            first_summary_row = hepatovirus_rows[0] if isinstance(hepatovirus_rows[0], list) else []
            if "说明" in hepatovirus_columns:
                note_index = hepatovirus_columns.index("说明")
                summary_note = _sanitize_virus_demo_note(first_summary_row[note_index] if note_index < len(first_summary_row) else "")
        notes: list[str] = []
        if broad_type != "-":
            notes.append(f"大亚型：{broad_type}")
        if subtype != "-":
            notes.append(f"子亚型：{subtype}")
        if isolate_label not in {"", "-"}:
            notes.append(f"参考株：{isolate_label}")
        if selected_reference != "-":
            notes.append(f"参考序列：{selected_reference}")
        if selected_gff != "-":
            notes.append(f"注释文件：{selected_gff}")
        if summary_note and summary_note != "-":
            notes.append(summary_note)
        mutation_table = _read_rsv_variant_annotation_table(report_dir)
        species_label = _hepatovirus_selection_value("species_label", "")
        if not species_label or species_label == "-":
            species_label = {
                "HAV": "Hepatitis A virus",
                "HBV": "Hepatitis B virus",
                "HCV": "Hepatitis C virus",
                "HDV": "Hepatitis D virus",
                "HEV": "Hepatitis E virus",
            }.get(broad_type, "Hepatovirus")
        display_columns = ["样本名称", "病毒类型", "大亚型", "子亚型", "参考序列", "覆盖度", "平均深度", "说明"]
        display_rows = [[
            sample_name or "-",
            species_label,
            broad_type,
            subtype,
            selected_reference,
            coverage_value,
            mean_depth_value,
            summary_note or "先基于肝炎病毒 broad 参考库完成大亚型判定；再进入对应大亚型参考库选择覆盖度最优的子亚型参考。",
        ]]
        return {
            "status": "ready",
            "mode": "hepatovirus_typing",
            "predicted_clade": subtype if subtype != "-" else broad_type,
            "predicted_subtype": subtype,
            "predicted_group": broad_type,
            "reference_name": selected_reference,
            "summary_cards": [
                {"label": "大亚型", "value": broad_type},
                {"label": "子亚型", "value": subtype},
                {"label": "参考序列", "value": selected_reference},
                {"label": "覆盖度", "value": coverage_value},
            ],
            "quality_metrics": [
                {"label": "平均深度", "value": mean_depth_value},
                {"label": "覆盖碱基", "value": covered_bases},
                {"label": "支持 reads", "value": num_reads},
            ],
            "sequence_name": sample_name or "-",
            "notes": "；".join([item for item in notes if item]),
            "mutation_table": mutation_table,
            "consensus_typing": hepatovirus_consensus_typing,
            "igv": _discover_hepatovirus_igv_assets(report_dir),
            "columns": display_columns,
            "rows": display_rows,
        }

    bandavirus_summary = _read_tsv_rows(report_dir / f"{sample_name}_serotype_result.tsv")
    bandavirus_columns = bandavirus_summary.get("columns", [])
    bandavirus_rows = bandavirus_summary.get("rows", [])
    bandavirus_selection = _read_tsv_rows(report_dir / f"{sample_name}_bandavirus_reference_selection" / "selection.tsv")
    bandavirus_consensus_typing = _read_tsv_rows(report_dir / f"{sample_name}_bandavirus_reference_selection" / "consensus_typing.tsv")
    bandavirus_af_segments = _read_tsv_rows(report_dir / f"{sample_name}_bandavirus_reference_selection" / "selected_segments.tsv")
    bandavirus_cj_segments = _read_tsv_rows(report_dir / f"{sample_name}_bandavirus_reference_selection" / "cj_typing" / "selected_segments.tsv")
    has_bandavirus_assets = (
        any(str(column).strip() in {"A_F(LMS)", "CJ(LMS)"} for column in bandavirus_columns)
        or bool(bandavirus_selection.get("rows"))
        or bool(bandavirus_consensus_typing.get("rows"))
        or bool(bandavirus_af_segments.get("rows"))
        or bool(bandavirus_cj_segments.get("rows"))
    )
    if has_bandavirus_assets:
        summary_first_row = bandavirus_rows[0] if bandavirus_rows and isinstance(bandavirus_rows[0], list) else []
        selection_columns = bandavirus_selection.get("columns", [])
        selection_rows = bandavirus_selection.get("rows", [])
        first_selection = selection_rows[0] if selection_rows and isinstance(selection_rows[0], list) else []
        consensus_columns = bandavirus_consensus_typing.get("columns", [])
        consensus_rows = bandavirus_consensus_typing.get("rows", [])
        first_consensus = consensus_rows[0] if consensus_rows and isinstance(consensus_rows[0], list) else []

        def _bandavirus_summary_value(column_name: str, fallback: str = "-") -> str:
            if column_name not in bandavirus_columns:
                return fallback
            index = bandavirus_columns.index(column_name)
            value = str(summary_first_row[index] if index < len(summary_first_row) else "").strip()
            return value or fallback

        def _bandavirus_selection_value(column_name: str, fallback: str = "-") -> str:
            if column_name not in selection_columns:
                return fallback
            index = selection_columns.index(column_name)
            value = str(first_selection[index] if index < len(first_selection) else "").strip()
            return value or fallback

        def _bandavirus_consensus_value(column_name: str, fallback: str = "-") -> str:
            if column_name not in consensus_columns:
                return fallback
            index = consensus_columns.index(column_name)
            value = str(first_consensus[index] if index < len(first_consensus) else "").strip()
            return value or fallback

        broad_type = _bandavirus_selection_value("broad_type", _bandavirus_consensus_value("broad_type", "-")).upper()
        af_group = _bandavirus_selection_value("af_group", _bandavirus_consensus_value("af_group", "-")).upper()
        af_lms = _bandavirus_summary_value("A_F(LMS)", _bandavirus_selection_value("segment_groups", "-"))
        cj_lms = _bandavirus_summary_value("CJ(LMS)", "-")
        typing_result = _bandavirus_summary_value("分型结果", "-")
        reference_name = _bandavirus_summary_value("参考命中", "-")
        summary_note = _sanitize_virus_demo_note(_bandavirus_summary_value("说明", ""))
        reassortment_flag = _bandavirus_selection_value("reassortment_flag", _bandavirus_consensus_value("reassortment_flag", "no"))
        segment_groups = _bandavirus_selection_value("segment_groups", _bandavirus_consensus_value("segment_groups", "-"))
        notes: list[str] = []
        if broad_type not in {"", "-"}:
            notes.append(f"大亚型：{broad_type}")
        if af_lms not in {"", "-"}:
            notes.append(f"A_F(LMS)：{af_lms}")
        if cj_lms not in {"", "-"}:
            notes.append(f"CJ(LMS)：{cj_lms}")
        if segment_groups not in {"", "-"}:
            notes.append(f"三片段判定：{segment_groups}")
        if reassortment_flag == "yes":
            notes.append("L/M/S 最优分型不一致，提示疑似重组/重配")
        if summary_note and summary_note != "-":
            notes.append(summary_note)
        mutation_table = _read_rsv_variant_annotation_table(report_dir)
        display_columns = ["样本名称", "病毒类型", "大亚型", "A_F(LMS)", "CJ(LMS)", "分型结果", "参考命中", "说明"]
        display_rows = [[
            sample_name or "-",
            _bandavirus_summary_value("病毒类型", "Bandavirus"),
            broad_type,
            af_lms,
            cj_lms,
            typing_result,
            reference_name,
            summary_note or "先根据 Bandavirus reference_genomes 判定大亚型；若为 SFTSV，再结合 A_F 分型与 CJ 三片段分型结果完成参考选择和结果汇总。",
        ]]
        return {
            "status": "ready",
            "mode": "bandavirus_typing",
            "predicted_clade": af_group if af_group not in {"", "-"} else af_lms,
            "predicted_group": broad_type,
            "predicted_lineage": cj_lms,
            "reference_name": reference_name,
            "summary_cards": [
                {"label": "大亚型", "value": broad_type or "--"},
                {"label": "A_F(LMS)", "value": af_lms or "--"},
                {"label": "CJ(LMS)", "value": cj_lms or "--"},
                {"label": "重组提示", "value": "疑似重组/重配" if reassortment_flag == "yes" else "未见异常"},
            ],
            "quality_metrics": [
                {"label": "A_F 子亚型", "value": af_group or "--"},
                {"label": "三片段判定", "value": segment_groups or "--"},
                {"label": "A_F 片段数", "value": str(len(bandavirus_af_segments.get("rows") or []))},
                {"label": "CJ 片段数", "value": str(len(bandavirus_cj_segments.get("rows") or []))},
            ],
            "sequence_name": sample_name or "-",
            "notes": "；".join([item for item in notes if item]),
            "mutation_table": mutation_table,
            "bandavirus_selection": bandavirus_selection,
            "consensus_typing": bandavirus_consensus_typing,
            "af_segment_typing": bandavirus_af_segments,
            "cj_segment_typing": bandavirus_cj_segments,
            "igv": _discover_bandavirus_igv_assets(report_dir),
            "columns": display_columns,
            "rows": display_rows,
        }

    orthohantavirus_summary = _read_tsv_rows(report_dir / f"{sample_name}_serotype_result.tsv")
    orthohantavirus_columns = orthohantavirus_summary.get("columns", [])
    orthohantavirus_rows = orthohantavirus_summary.get("rows", [])
    orthohantavirus_selection = _read_tsv_rows(report_dir / f"{sample_name}_orthohantavirus_reference_selection" / "selection.tsv")
    orthohantavirus_consensus_typing = _read_tsv_rows(report_dir / f"{sample_name}_orthohantavirus_reference_selection" / "consensus_typing.tsv")
    orthohantavirus_segments = _read_tsv_rows(report_dir / f"{sample_name}_orthohantavirus_reference_selection" / "selected_segments.tsv")
    orthohantavirus_broad_typing = _read_tsv_rows(report_dir / f"{sample_name}_orthohantavirus_reference_selection" / "broad_typing" / "typing_summary.tsv")
    has_orthohantavirus_assets = (
        any(str(column).strip() in {"S分型", "LMS参考"} for column in orthohantavirus_columns)
        or bool(orthohantavirus_selection.get("rows"))
        or bool(orthohantavirus_consensus_typing.get("rows"))
        or bool(orthohantavirus_segments.get("rows"))
        or bool(orthohantavirus_broad_typing.get("rows"))
    )
    if has_orthohantavirus_assets:
        summary_first_row = orthohantavirus_rows[0] if orthohantavirus_rows and isinstance(orthohantavirus_rows[0], list) else []
        selection_columns = orthohantavirus_selection.get("columns", [])
        selection_rows = orthohantavirus_selection.get("rows", [])
        first_selection = selection_rows[0] if selection_rows and isinstance(selection_rows[0], list) else []
        consensus_columns = orthohantavirus_consensus_typing.get("columns", [])
        consensus_rows = orthohantavirus_consensus_typing.get("rows", [])
        first_consensus = consensus_rows[0] if consensus_rows and isinstance(consensus_rows[0], list) else []
        broad_columns = orthohantavirus_broad_typing.get("columns", [])
        broad_rows = orthohantavirus_broad_typing.get("rows", [])
        first_broad = broad_rows[0] if broad_rows and isinstance(broad_rows[0], list) else []

        def _orthohantavirus_summary_value(column_name: str, fallback: str = "-") -> str:
            if column_name not in orthohantavirus_columns:
                return fallback
            index = orthohantavirus_columns.index(column_name)
            value = str(summary_first_row[index] if index < len(summary_first_row) else "").strip()
            return value or fallback

        def _orthohantavirus_selection_value(column_name: str, fallback: str = "-") -> str:
            if column_name not in selection_columns:
                return fallback
            index = selection_columns.index(column_name)
            value = str(first_selection[index] if index < len(first_selection) else "").strip()
            return value or fallback

        def _orthohantavirus_consensus_value(column_name: str, fallback: str = "-") -> str:
            if column_name not in consensus_columns:
                return fallback
            index = consensus_columns.index(column_name)
            value = str(first_consensus[index] if index < len(first_consensus) else "").strip()
            return value or fallback

        def _orthohantavirus_broad_value(column_name: str, fallback: str = "-") -> str:
            if column_name not in broad_columns:
                return fallback
            index = broad_columns.index(column_name)
            value = str(first_broad[index] if index < len(first_broad) else "").strip()
            return value or fallback

        predicted_type = _orthohantavirus_selection_value("predicted_type", _orthohantavirus_consensus_value("predicted_type", "-")).upper()
        s_segment_type = _orthohantavirus_summary_value("S分型", _orthohantavirus_selection_value("s_segment_type", _orthohantavirus_consensus_value("s_segment_type", "-"))).upper()
        selected_segments = _orthohantavirus_summary_value("LMS参考", _orthohantavirus_selection_value("selected_segments", _orthohantavirus_consensus_value("selected_segments", "-")))
        typing_result = _orthohantavirus_summary_value("分型结果", predicted_type or "-")
        reference_name = _orthohantavirus_summary_value("参考命中", "-")
        summary_note = _sanitize_virus_demo_note(_orthohantavirus_summary_value("说明", ""))
        segment_count = _orthohantavirus_selection_value("segment_count", _orthohantavirus_consensus_value("segment_count", str(len(orthohantavirus_segments.get("rows") or []))))
        broad_segment_count = _orthohantavirus_broad_value("segment_count", segment_count or "--")
        broad_coverage_sum = _orthohantavirus_broad_value("coverage_sum", "--")
        broad_depth_sum = _orthohantavirus_broad_value("depth_sum", "--")
        broad_reads_sum = _orthohantavirus_broad_value("num_reads_sum", "--")
        notes: list[str] = []
        if predicted_type not in {"", "-"}:
            notes.append(f"分型：{predicted_type}")
        if s_segment_type not in {"", "-"}:
            notes.append(f"S片段：{s_segment_type}")
        if selected_segments not in {"", "-"}:
            notes.append(f"L/M/S参考：{selected_segments}")
        if broad_segment_count not in {"", "-"}:
            notes.append(f"broad 命中片段数：{broad_segment_count}")
        if broad_coverage_sum not in {"", "-"}:
            notes.append(f"broad coverage_sum：{broad_coverage_sum}")
        if summary_note and summary_note != "-":
            notes.append(summary_note)
        mutation_table = _read_influenza_variant_annotation_table(report_dir)
        display_columns = ["样本名称", "病毒类型", "S分型", "LMS参考", "分型结果", "参考命中", "说明"]
        display_rows = [[
            sample_name or "-",
            _orthohantavirus_summary_value("病毒类型", "Orthohantavirus"),
            s_segment_type,
            selected_segments,
            typing_result,
            reference_name,
            summary_note or "基于 Orthohantavirus 本地三片段参考库进行型别筛选，并优先以 S 片段分型结果辅助解释。",
        ]]
        return {
            "status": "ready",
            "mode": "orthohantavirus_typing",
            "predicted_clade": predicted_type,
            "predicted_group": s_segment_type,
            "reference_name": reference_name,
            "summary_cards": [
                {"label": "Orthohantavirus 分型", "value": predicted_type or "--"},
                {"label": "S片段分型", "value": s_segment_type or "--"},
                {"label": "L/M/S参考", "value": selected_segments or "--"},
                {"label": "参考片段数", "value": segment_count or "--"},
            ],
            "quality_metrics": [
                {"label": "broad 片段数", "value": broad_segment_count or "--"},
                {"label": "broad coverage_sum", "value": broad_coverage_sum or "--"},
                {"label": "broad depth_sum", "value": broad_depth_sum or "--"},
                {"label": "支持 reads", "value": broad_reads_sum or "--"},
            ],
            "sequence_name": sample_name or "-",
            "notes": "；".join([item for item in notes if item]),
            "orthohantavirus_selection": orthohantavirus_selection,
            "broad_typing": orthohantavirus_broad_typing,
            "consensus_typing": orthohantavirus_consensus_typing,
            "segment_typing": orthohantavirus_segments,
            "mutation_table": mutation_table,
            "igv": _discover_orthohantavirus_igv_assets(report_dir),
            "columns": display_columns,
            "rows": display_rows,
        }

    orthoebolavirus_summary = _read_tsv_rows(report_dir / f"{sample_name}_serotype_result.tsv")
    orthoebolavirus_columns = orthoebolavirus_summary.get("columns", [])
    orthoebolavirus_rows = orthoebolavirus_summary.get("rows", [])
    orthoebolavirus_selection = _read_tsv_rows(report_dir / f"{sample_name}_orthoebolavirus_reference_selection" / "selection.tsv")
    orthoebolavirus_consensus_typing = _read_tsv_rows(report_dir / f"{sample_name}_orthoebolavirus_reference_selection" / "consensus_typing" / "consensus_typing.tsv")
    has_orthoebolavirus_assets = (
        any(str(column).strip() == "Orthoebolavirus分型" for column in orthoebolavirus_columns)
        or bool(orthoebolavirus_selection.get("rows"))
        or bool(orthoebolavirus_consensus_typing.get("rows"))
    )
    if has_orthoebolavirus_assets:
        summary_first_row = orthoebolavirus_rows[0] if orthoebolavirus_rows and isinstance(orthoebolavirus_rows[0], list) else []
        selection_columns = orthoebolavirus_selection.get("columns", [])
        selection_rows = orthoebolavirus_selection.get("rows", [])
        first_selection = selection_rows[0] if selection_rows and isinstance(selection_rows[0], list) else []
        consensus_columns = orthoebolavirus_consensus_typing.get("columns", [])
        consensus_rows = orthoebolavirus_consensus_typing.get("rows", [])
        first_consensus = consensus_rows[0] if consensus_rows and isinstance(consensus_rows[0], list) else []

        def _orthoebolavirus_summary_value(column_name: str, fallback: str = "-") -> str:
            if column_name not in orthoebolavirus_columns:
                return fallback
            index = orthoebolavirus_columns.index(column_name)
            value = str(summary_first_row[index] if index < len(summary_first_row) else "").strip()
            return value or fallback

        def _orthoebolavirus_selection_value(column_name: str, fallback: str = "-") -> str:
            if column_name not in selection_columns:
                return fallback
            index = selection_columns.index(column_name)
            value = str(first_selection[index] if index < len(first_selection) else "").strip()
            return value or fallback

        def _orthoebolavirus_consensus_value(column_name: str, fallback: str = "-") -> str:
            if column_name not in consensus_columns:
                return fallback
            index = consensus_columns.index(column_name)
            value = str(first_consensus[index] if index < len(first_consensus) else "").strip()
            return value or fallback

        predicted_type = _orthoebolavirus_selection_value("predicted_type", _orthoebolavirus_consensus_value("predicted_type", "-")).upper()
        species_label = _orthoebolavirus_selection_value("species_label", _orthoebolavirus_summary_value("病毒类型", "Orthoebolavirus"))
        virus_name = _orthoebolavirus_selection_value("virus_name", _orthoebolavirus_summary_value("病毒名称", "-"))
        reference_name = _orthoebolavirus_selection_value("reference_name", _orthoebolavirus_summary_value("参考序列", "-"))
        coverage_value = _orthoebolavirus_selection_value("coverage", _orthoebolavirus_summary_value("覆盖度", "-"))
        mean_depth_value = _orthoebolavirus_selection_value("mean_depth", _orthoebolavirus_summary_value("平均深度", "-"))
        summary_note = _sanitize_virus_demo_note(_orthoebolavirus_selection_value("note", _orthoebolavirus_summary_value("说明", "")))
        mutation_table = _read_influenza_variant_annotation_table(report_dir)
        display_columns = ["样本名称", "病毒类型", "Orthoebolavirus分型", "病毒名称", "参考命中", "覆盖度", "平均深度", "说明"]
        display_rows = [[
            sample_name or "-",
            species_label,
            predicted_type or "-",
            virus_name or "-",
            reference_name or "-",
            coverage_value,
            mean_depth_value,
            summary_note or "基于 Orthoebolavirus 本地完整/编码完整参考基因组库进行最优参考筛选；EBOV 样本另接 Ebola Nextclade。",
        ]]
        return {
            "status": "ready",
            "mode": "orthoebolavirus_typing",
            "predicted_clade": predicted_type,
            "predicted_group": virus_name,
            "reference_name": reference_name,
            "summary_cards": [
                {"label": "Orthoebolavirus 分型", "value": predicted_type or "--"},
                {"label": "病毒名称", "value": virus_name or "--"},
                {"label": "参考序列", "value": reference_name or "--"},
                {"label": "覆盖度", "value": coverage_value or "--"},
            ],
            "quality_metrics": [
                {"label": "平均深度", "value": mean_depth_value or "--"},
                {"label": "覆盖碱基", "value": _orthoebolavirus_selection_value("covered_bases", "--")},
                {"label": "支持 reads", "value": _orthoebolavirus_selection_value("num_reads", "--")},
                {"label": "Accession", "value": _orthoebolavirus_selection_value("accession", "--")},
            ],
            "sequence_name": sample_name or "-",
            "notes": summary_note,
            "orthoebolavirus_selection": orthoebolavirus_selection,
            "consensus_typing": orthoebolavirus_consensus_typing,
            "mutation_table": mutation_table,
            "igv": _discover_ebola_igv_assets(report_dir),
            "columns": display_columns,
            "rows": display_rows,
        }

    astroviridae_summary = _read_tsv_rows(report_dir / f"{sample_name}_serotype_result.tsv")
    astroviridae_columns = astroviridae_summary.get("columns", [])
    astroviridae_rows = astroviridae_summary.get("rows", [])
    astroviridae_selection = _read_tsv_rows(report_dir / f"{sample_name}_astroviridae_reference_selection" / "selection.tsv")
    astroviridae_consensus_typing = _read_tsv_rows(report_dir / f"{sample_name}_astroviridae_reference_selection" / "consensus_typing" / "consensus_typing.tsv")
    astroviridae_phylogeny = _build_astroviridae_gene_phylogeny(report_dir)
    has_astroviridae_assets = (
        any(str(column).strip() in {"ORF2分型", "病毒属", "病毒种"} for column in astroviridae_columns)
        or bool(astroviridae_selection.get("rows"))
        or bool(astroviridae_consensus_typing.get("rows"))
        or str(astroviridae_phylogeny.get("status") or "").strip().lower() == "ready"
    )
    if has_astroviridae_assets:
        selection_columns = astroviridae_selection.get("columns", [])
        selection_rows = astroviridae_selection.get("rows", [])
        first_selection = selection_rows[0] if selection_rows and isinstance(selection_rows[0], list) else []

        def _astroviridae_selection_value(column_name: str, fallback: str = "-") -> str:
            if column_name not in selection_columns:
                return fallback
            index = selection_columns.index(column_name)
            value = str(first_selection[index] if index < len(first_selection) else "").strip()
            return value or fallback

        consensus_columns = astroviridae_consensus_typing.get("columns", [])
        consensus_rows = astroviridae_consensus_typing.get("rows", [])
        first_consensus = consensus_rows[0] if consensus_rows and isinstance(consensus_rows[0], list) else []

        def _astroviridae_consensus_value(column_name: str, fallback: str = "-") -> str:
            if column_name not in consensus_columns:
                return fallback
            index = consensus_columns.index(column_name)
            value = str(first_consensus[index] if index < len(first_consensus) else "").strip()
            return value or fallback

        orf2_type = _astroviridae_selection_value("orf2_type", _astroviridae_consensus_value("orf2_type", "-"))
        genus = _astroviridae_selection_value("genus", _astroviridae_consensus_value("genus", "-"))
        species = _astroviridae_selection_value("species", _astroviridae_consensus_value("species", "-"))
        selected_reference = _astroviridae_selection_value("reference_name", Path(_astroviridae_selection_value("reference_path", "")).name or "-")
        selected_gff = Path(_astroviridae_selection_value("gff_path", "")).name if _astroviridae_selection_value("gff_path", "") not in {"", "-", "nogtf"} else "-"
        coverage_value = _astroviridae_selection_value("coverage", "-")
        mean_depth_value = _astroviridae_selection_value("mean_depth", "-")
        covered_bases = _astroviridae_selection_value("covered_bases", "-")
        num_reads = _astroviridae_selection_value("num_reads", "-")
        anchor_accession = _astroviridae_selection_value("anchor_accession", "-")
        candidate_count = _astroviridae_selection_value("candidate_count", "-")
        dedup_candidate_count = _astroviridae_selection_value("dedup_candidate_count", "-")
        if coverage_value not in {"-", ""}:
            try:
                coverage_value = f"{float(coverage_value):.2f}%"
            except ValueError:
                pass
        if mean_depth_value not in {"-", ""}:
            try:
                mean_depth_value = f"{float(mean_depth_value):.2f}"
            except ValueError:
                pass
        summary_note = ""
        if astroviridae_columns and astroviridae_rows:
            first_summary_row = astroviridae_rows[0] if isinstance(astroviridae_rows[0], list) else []
            if "说明" in astroviridae_columns:
                note_index = astroviridae_columns.index("说明")
                summary_note = _sanitize_virus_demo_note(first_summary_row[note_index] if note_index < len(first_summary_row) else "")
        notes: list[str] = []
        if orf2_type != "-":
            notes.append(f"ORF2 分型：{orf2_type}")
        if genus != "-":
            notes.append(f"病毒属：{genus}")
        if species != "-":
            notes.append(f"病毒种：{species}")
        if selected_reference != "-":
            notes.append(f"参考序列：{selected_reference}")
        if selected_gff != "-":
            notes.append(f"注释文件：{selected_gff}")
        if anchor_accession not in {"", "-"}:
            notes.append(f"批次锚点参考：{anchor_accession}")
        if candidate_count not in {"", "-"} and dedup_candidate_count not in {"", "-"}:
            notes.append(f"候选基因组：{candidate_count} 条，95% 去冗余后 {dedup_candidate_count} 条")
        if summary_note and summary_note != "-":
            notes.append(summary_note)
        mutation_table = _read_rsv_variant_annotation_table(report_dir)
        display_columns = ["样本名称", "病毒种", "病毒属", "ORF2分型", "参考序列", "覆盖度", "平均深度", "说明"]
        display_rows = [[
            sample_name or "-",
            species if species != "-" else "Mamastrovirus/Avastrovirus",
            genus,
            orf2_type,
            selected_reference,
            coverage_value,
            mean_depth_value,
            summary_note or "基于 Astroviridae ORF2 分型结果，在对应亚型全基因组候选集中选择覆盖度最优参考；组装后结合 VADR 注释并提取 ORF2 构建系统发育树进行二次分类。",
        ]]
        return {
            "status": "ready",
            "mode": "astroviridae_typing",
            "predicted_clade": orf2_type,
            "predicted_group": genus,
            "predicted_lineage": species,
            "reference_name": selected_reference,
            "summary_cards": [
                {"label": "ORF2 分型", "value": orf2_type},
                {"label": "病毒属", "value": genus},
                {"label": "病毒种", "value": species},
                {"label": "参考序列", "value": selected_reference},
                {"label": "覆盖度", "value": coverage_value},
            ],
            "quality_metrics": [
                {"label": "平均深度", "value": mean_depth_value},
                {"label": "覆盖碱基", "value": covered_bases},
                {"label": "支持 reads", "value": num_reads},
                {"label": "去冗余候选数", "value": dedup_candidate_count},
            ],
            "sequence_name": sample_name or "-",
            "notes": "；".join([item for item in notes if item]),
            "mutation_table": mutation_table,
            "consensus_typing": astroviridae_consensus_typing,
            "igv": _discover_astroviridae_igv_assets(report_dir),
            "gene_phylogeny": astroviridae_phylogeny,
            "columns": display_columns,
            "rows": display_rows,
        }

    rhinovirus_summary = _read_tsv_rows(report_dir / f"{sample_name}_serotype_result.tsv")
    rhinovirus_columns = rhinovirus_summary.get("columns", [])
    rhinovirus_rows = rhinovirus_summary.get("rows", [])
    rhinovirus_selection = _read_tsv_rows(report_dir / f"{sample_name}_rhinovirus_reference_selection" / "selection.tsv")
    rhinovirus_consensus_typing = _read_tsv_rows(report_dir / f"{sample_name}_rhinovirus_reference_selection" / "consensus_typing" / "consensus_typing.tsv")
    has_rhinovirus_assets = (
        any(str(column).strip() in {"VP1分型", "物种组"} for column in rhinovirus_columns)
        or bool(rhinovirus_selection.get("rows"))
        or bool(rhinovirus_consensus_typing.get("rows"))
    )
    if has_rhinovirus_assets:
        selection_columns = rhinovirus_selection.get("columns", [])
        selection_rows = rhinovirus_selection.get("rows", [])
        first_selection = selection_rows[0] if selection_rows and isinstance(selection_rows[0], list) else []

        def _rhinovirus_selection_value(column_name: str, fallback: str = "-") -> str:
            if column_name not in selection_columns:
                return fallback
            index = selection_columns.index(column_name)
            value = str(first_selection[index] if index < len(first_selection) else "").strip()
            return value or fallback

        consensus_columns = rhinovirus_consensus_typing.get("columns", [])
        consensus_rows = rhinovirus_consensus_typing.get("rows", [])
        first_consensus = consensus_rows[0] if consensus_rows and isinstance(consensus_rows[0], list) else []

        def _rhinovirus_consensus_value(column_name: str, fallback: str = "-") -> str:
            if column_name not in consensus_columns:
                return fallback
            index = consensus_columns.index(column_name)
            value = str(first_consensus[index] if index < len(first_consensus) else "").strip()
            return value or fallback

        vp1_type = _rhinovirus_selection_value("vp1_type", _rhinovirus_consensus_value("vp1_type", "-"))
        species_group = _rhinovirus_selection_value("species_group", _rhinovirus_consensus_value("species_group", "-")).upper()
        selected_reference = _rhinovirus_selection_value("reference_name", Path(_rhinovirus_selection_value("reference_path", "")).name or "-")
        selected_gff = Path(_rhinovirus_selection_value("gff_path", "")).name if _rhinovirus_selection_value("gff_path", "") not in {"", "-", "nogtf"} else "-"
        coverage_value = _rhinovirus_selection_value("coverage", "-")
        mean_depth_value = _rhinovirus_selection_value("mean_depth", "-")
        covered_bases = _rhinovirus_selection_value("covered_bases", "-")
        num_reads = _rhinovirus_selection_value("num_reads", "-")
        if coverage_value not in {"-", ""}:
            try:
                coverage_value = f"{float(coverage_value):.2f}%"
            except ValueError:
                pass
        if mean_depth_value not in {"-", ""}:
            try:
                mean_depth_value = f"{float(mean_depth_value):.2f}"
            except ValueError:
                pass
        summary_note = ""
        if rhinovirus_columns and rhinovirus_rows:
            first_summary_row = rhinovirus_rows[0] if isinstance(rhinovirus_rows[0], list) else []
            if "说明" in rhinovirus_columns:
                note_index = rhinovirus_columns.index("说明")
                summary_note = _sanitize_virus_demo_note(first_summary_row[note_index] if note_index < len(first_summary_row) else "")
        notes: list[str] = []
        if vp1_type != "-":
            notes.append(f"VP1 分型：{vp1_type}")
        if species_group != "-":
            notes.append(f"物种组：{species_group}")
        if selected_reference != "-":
            notes.append(f"参考序列：{selected_reference}")
        if selected_gff != "-":
            notes.append(f"注释文件：{selected_gff}")
        if summary_note and summary_note != "-":
            notes.append(summary_note)
        mutation_table = _read_rsv_variant_annotation_table(report_dir)
        gene_phylogeny = _build_rhinovirus_gene_phylogeny(report_dir)
        display_columns = ["样本名称", "病毒类型", "VP1分型", "物种组", "参考序列", "覆盖度", "平均深度", "说明"]
        display_rows = [[
            sample_name or "-",
            f"Rhinovirus {species_group}" if species_group != "-" else "Rhinovirus",
            vp1_type,
            species_group,
            selected_reference,
            coverage_value,
            mean_depth_value,
            summary_note or "基于 VP1 分型后选择最优全基因组参考，并结合 VADR 注释生成鼻病毒报告。",
        ]]
        return {
            "status": "ready",
            "mode": "rhinovirus_typing",
            "predicted_clade": vp1_type,
            "predicted_group": species_group,
            "reference_name": selected_reference,
            "summary_cards": [
                {"label": "VP1 分型", "value": vp1_type},
                {"label": "物种组", "value": species_group},
                {"label": "参考序列", "value": selected_reference},
                {"label": "覆盖度", "value": coverage_value},
            ],
            "quality_metrics": [
                {"label": "平均深度", "value": mean_depth_value},
                {"label": "覆盖碱基", "value": covered_bases},
                {"label": "支持 reads", "value": num_reads},
            ],
            "sequence_name": sample_name or "-",
            "notes": "；".join([item for item in notes if item]),
            "mutation_table": mutation_table,
            "consensus_typing": rhinovirus_consensus_typing,
            "igv": _discover_rhinovirus_igv_assets(report_dir),
            "gene_phylogeny": gene_phylogeny,
            "columns": display_columns,
            "rows": display_rows,
        }

    seasonal_hcov_summary = _read_tsv_rows(report_dir / f"{sample_name}_serotype_result.tsv")
    seasonal_hcov_columns = seasonal_hcov_summary.get("columns", [])
    seasonal_hcov_rows = seasonal_hcov_summary.get("rows", [])
    seasonal_hcov_selection = _read_tsv_rows(report_dir / f"{sample_name}_seasonal_hcov_reference_selection" / "selection.tsv")
    seasonal_hcov_consensus_typing = _read_tsv_rows(report_dir / f"{sample_name}_seasonal_hcov_reference_selection" / "consensus_typing" / "consensus_typing.tsv")
    seasonal_hcov_phylogeny = _build_seasonal_hcov_spike_phylogeny(report_dir)
    has_seasonal_hcov_assets = (
        any(str(column).strip() in {"大类分型", "S子亚型"} for column in seasonal_hcov_columns)
        or bool(seasonal_hcov_selection.get("rows"))
        or bool(seasonal_hcov_consensus_typing.get("rows"))
        or str(seasonal_hcov_phylogeny.get("status") or "").strip().lower() == "ready"
    )
    if has_seasonal_hcov_assets:
        selection_columns = seasonal_hcov_selection.get("columns", [])
        selection_rows = seasonal_hcov_selection.get("rows", [])
        first_selection = selection_rows[0] if selection_rows and isinstance(selection_rows[0], list) else []

        def _seasonal_hcov_selection_value(column_name: str, fallback: str = "-") -> str:
            if column_name not in selection_columns:
                return fallback
            index = selection_columns.index(column_name)
            value = str(first_selection[index] if index < len(first_selection) else "").strip()
            return value or fallback

        summary_columns = seasonal_hcov_columns
        summary_rows = seasonal_hcov_rows
        first_summary = summary_rows[0] if summary_rows and isinstance(summary_rows[0], list) else []

        def _seasonal_hcov_summary_value(column_name: str, fallback: str = "-") -> str:
            if column_name not in summary_columns:
                return fallback
            index = summary_columns.index(column_name)
            value = str(first_summary[index] if index < len(first_summary) else "").strip()
            return value or fallback

        consensus_columns = seasonal_hcov_consensus_typing.get("columns", [])
        consensus_rows = seasonal_hcov_consensus_typing.get("rows", [])
        first_consensus = consensus_rows[0] if consensus_rows and isinstance(consensus_rows[0], list) else []

        def _seasonal_hcov_consensus_value(column_name: str, fallback: str = "-") -> str:
            if column_name not in consensus_columns:
                return fallback
            index = consensus_columns.index(column_name)
            value = str(first_consensus[index] if index < len(first_consensus) else "").strip()
            return value or fallback

        major_type = _seasonal_hcov_summary_value("大类分型", _seasonal_hcov_selection_value("hcov_type", "-"))
        spike_subtype = _seasonal_hcov_summary_value("S子亚型", _seasonal_hcov_consensus_value("subtype", "-"))
        selected_reference = _seasonal_hcov_summary_value("病毒类型", _seasonal_hcov_selection_value("reference_name", Path(_seasonal_hcov_selection_value("reference_path", "")).name or "-"))
        nearest_reference = _seasonal_hcov_summary_value("最近参考", _seasonal_hcov_consensus_value("nearest_tree_label", "-"))
        selected_gff = _seasonal_hcov_summary_value("注释文件", Path(_seasonal_hcov_selection_value("gff_path", "")).name if _seasonal_hcov_selection_value("gff_path", "") not in {"", "-", "nogtf"} else "-")
        coverage_value = _seasonal_hcov_summary_value("覆盖度", _seasonal_hcov_selection_value("coverage", "-"))
        mean_depth_value = _seasonal_hcov_summary_value("平均深度", _seasonal_hcov_selection_value("mean_depth", "-"))
        covered_bases = _seasonal_hcov_selection_value("covered_bases", "-")
        num_reads = _seasonal_hcov_selection_value("num_reads", "-")
        if coverage_value not in {"-", ""} and not str(coverage_value).endswith("%"):
            try:
                coverage_value = f"{float(coverage_value):.2f}%"
            except ValueError:
                pass
        if mean_depth_value not in {"-", ""}:
            try:
                mean_depth_value = f"{float(mean_depth_value):.2f}"
            except ValueError:
                pass
        summary_note = _sanitize_virus_demo_note(_seasonal_hcov_summary_value("说明", ""))
        notes: list[str] = []
        if major_type != "-":
            notes.append(f"大类分型：{major_type}")
        if spike_subtype != "-":
            notes.append(f"S 子亚型：{spike_subtype}")
        if nearest_reference != "-":
            notes.append(f"最近参考：{nearest_reference}")
        if selected_gff != "-":
            notes.append(f"注释文件：{selected_gff}")
        if summary_note and summary_note != "-":
            notes.append(summary_note)
        mutation_table = _read_rsv_variant_annotation_table(report_dir)
        display_columns = ["样本名称", "病毒类型", "大类分型", "S子亚型", "最近参考", "注释文件", "覆盖度", "平均深度", "说明"]
        display_rows = [[
            sample_name or "-",
            selected_reference,
            major_type,
            spike_subtype,
            nearest_reference,
            selected_gff,
            coverage_value,
            mean_depth_value,
            summary_note or "基于季节性冠状病毒参考选择、VADR 注释和 S 基因系统发育树完成型别与子亚型判定。",
        ]]
        return {
            "status": "ready",
            "mode": "seasonal_hcov_typing",
            "predicted_clade": major_type,
            "predicted_subtype": spike_subtype,
            "reference_name": selected_reference,
            "nearest_reference": nearest_reference,
            "summary_cards": [
                {"label": "大类分型", "value": major_type},
                {"label": "S 子亚型", "value": spike_subtype},
                {"label": "参考序列", "value": selected_reference},
                {"label": "最近参考", "value": nearest_reference},
                {"label": "覆盖度", "value": coverage_value},
            ],
            "quality_metrics": [
                {"label": "平均深度", "value": mean_depth_value},
                {"label": "覆盖碱基", "value": covered_bases},
                {"label": "支持 reads", "value": num_reads},
            ],
            "sequence_name": sample_name or "-",
            "notes": "；".join([item for item in notes if item]),
            "mutation_table": mutation_table,
            "consensus_typing": seasonal_hcov_consensus_typing,
            "igv": _discover_seasonal_hcov_igv_assets(report_dir),
            "gene_phylogeny": seasonal_hcov_phylogeny,
            "columns": display_columns,
            "rows": display_rows,
        }

    rotavirus_summary = _read_tsv_rows(report_dir / f"{sample_name}_serotype_result.tsv")
    rotavirus_columns = rotavirus_summary.get("columns", [])
    rotavirus_rows = rotavirus_summary.get("rows", [])
    rotavirus_selection = _read_tsv_rows(report_dir / f"{sample_name}_rotavirus_reference_selection" / "selection.tsv")
    rotavirus_consensus_typing = _read_tsv_rows(report_dir / f"{sample_name}_rotavirus_reference_selection" / "consensus_typing.tsv")
    rotavirus_group_typing = _read_tsv_rows(report_dir / f"{sample_name}_rotavirus_reference_selection" / "group_typing" / "group_typing.tsv")
    rotavirus_subtype_typing = _read_tsv_rows(report_dir / f"{sample_name}_rotavirus_reference_selection" / "subtype_typing" / "subtype_typing.tsv")
    has_rotavirus_assets = (
        any(str(column).strip() in {"大组分型", "G分型", "P分型", "组合分型"} for column in rotavirus_columns)
        or bool(rotavirus_selection.get("rows"))
        or bool(rotavirus_group_typing.get("rows"))
        or bool(rotavirus_subtype_typing.get("rows"))
        or bool(rotavirus_consensus_typing.get("rows"))
    )
    if has_rotavirus_assets:
        selection_columns = rotavirus_selection.get("columns", [])
        selection_rows = rotavirus_selection.get("rows", [])
        first_selection = selection_rows[0] if selection_rows and isinstance(selection_rows[0], list) else []

        def _rotavirus_selection_value(column_name: str, fallback: str = "-") -> str:
            if column_name not in selection_columns:
                return fallback
            index = selection_columns.index(column_name)
            value = str(first_selection[index] if index < len(first_selection) else "").strip()
            return value or fallback

        consensus_columns = rotavirus_consensus_typing.get("columns", [])
        consensus_rows = rotavirus_consensus_typing.get("rows", [])
        first_consensus = consensus_rows[0] if consensus_rows and isinstance(consensus_rows[0], list) else []

        def _rotavirus_consensus_value(column_name: str, fallback: str = "-") -> str:
            if column_name not in consensus_columns:
                return fallback
            index = consensus_columns.index(column_name)
            value = str(first_consensus[index] if index < len(first_consensus) else "").strip()
            return value or fallback

        group_columns = rotavirus_group_typing.get("columns", [])
        group_rows = rotavirus_group_typing.get("rows", [])
        selected_group_row: list | None = None
        is_selected_index = group_columns.index("is_selected") if "is_selected" in group_columns else -1
        for row in group_rows:
            if not isinstance(row, list):
                continue
            if is_selected_index >= 0 and is_selected_index < len(row) and str(row[is_selected_index] or "").strip().lower() == "yes":
                selected_group_row = row
                break
        if selected_group_row is None and group_rows and isinstance(group_rows[0], list):
            selected_group_row = group_rows[0]

        def _rotavirus_group_value(column_name: str, fallback: str = "-") -> str:
            if selected_group_row is None or column_name not in group_columns:
                return fallback
            index = group_columns.index(column_name)
            value = str(selected_group_row[index] if index < len(selected_group_row) else "").strip()
            return value or fallback

        summary_columns = rotavirus_columns
        summary_rows = rotavirus_rows
        first_summary = summary_rows[0] if summary_rows and isinstance(summary_rows[0], list) else []

        def _rotavirus_summary_value(column_name: str, fallback: str = "-") -> str:
            if column_name not in summary_columns:
                return fallback
            index = summary_columns.index(column_name)
            value = str(first_summary[index] if index < len(first_summary) else "").strip()
            return value or fallback

        group_type = _rotavirus_summary_value(
            "大组分型",
            _rotavirus_selection_value("group_type", _rotavirus_consensus_value("group_type", _rotavirus_group_value("selected_group_type", _rotavirus_group_value("group_type", "-")))),
        )
        g_genotype = _rotavirus_summary_value("G分型", _rotavirus_selection_value("g_genotype", _rotavirus_consensus_value("g_genotype", "-")))
        p_genotype = _rotavirus_summary_value("P分型", _rotavirus_selection_value("p_genotype", _rotavirus_consensus_value("p_genotype", "-")))
        subtype_combo = _rotavirus_summary_value("组合分型", _rotavirus_selection_value("subtype_combo", _rotavirus_consensus_value("subtype_combo", "-")))
        isolate_label = _rotavirus_selection_value("isolate_label", _rotavirus_group_value("isolate_label", "-"))
        reference_accession = _rotavirus_group_value("reference_name", Path(_rotavirus_selection_value("reference_path", "")).stem or "-")
        selected_reference = isolate_label if isolate_label != "-" else reference_accession
        segment_count = _rotavirus_group_value("segment_count", "-")
        full_length_coverage = _rotavirus_group_value("full_length_coverage_pct", "-")
        covered_bases = _rotavirus_group_value("covered_bases_sum", "-")
        reference_bases = _rotavirus_group_value("reference_bases_sum", "-")
        supported_reads = _rotavirus_group_value("num_reads_sum", "-")
        coverage_sum = _rotavirus_group_value("coverage_sum", "-")
        if full_length_coverage not in {"-", ""} and not str(full_length_coverage).endswith("%"):
            try:
                full_length_coverage = f"{float(full_length_coverage):.2f}%"
            except ValueError:
                pass
        if coverage_sum not in {"-", ""}:
            try:
                coverage_sum = f"{float(coverage_sum):.2f}"
            except ValueError:
                pass
        summary_note = _sanitize_virus_demo_note(_rotavirus_summary_value("说明", ""))
        if not summary_note:
            summary_note = "先按 A/B/C 大组完整参考覆盖度选择最优组别，再对 A 组样本结合 VP4/VP7 分型确定最优参考株。"
        notes: list[str] = []
        if group_type != "-":
            notes.append(f"大组分型：{group_type}")
        if subtype_combo != "-" and subtype_combo:
            notes.append(f"组合分型：{subtype_combo}")
        elif g_genotype != "-" or p_genotype != "-":
            notes.append(f"G/P 分型：{g_genotype} / {p_genotype}")
        if selected_reference != "-":
            notes.append(f"最优参考株：{selected_reference}")
        if reference_accession != "-" and reference_accession != selected_reference:
            notes.append(f"参考 accession：{reference_accession}")
        if summary_note and summary_note != "-":
            notes.append(summary_note)
        display_columns = ["样本名称", "病毒类型", "大组分型", "G分型", "P分型", "组合分型", "最优参考株", "参考片段数", "全长覆盖度", "说明"]
        display_rows = [[
            sample_name or "-",
            f"Human rotavirus {group_type}" if group_type != "-" else "Human rotavirus",
            group_type,
            g_genotype,
            p_genotype,
            subtype_combo,
            selected_reference,
            segment_count,
            full_length_coverage,
            summary_note,
        ]]
        return {
            "status": "ready",
            "mode": "rotavirus_typing",
            "predicted_clade": group_type,
            "predicted_lineage": subtype_combo if subtype_combo not in {"", "-"} else "-",
            "predicted_group": group_type,
            "predicted_subtype": subtype_combo if subtype_combo not in {"", "-"} else "-",
            "reference_name": selected_reference,
            "reference_accession": reference_accession,
            "summary_cards": [
                {"label": "大组分型", "value": group_type},
                {"label": "G 分型", "value": g_genotype},
                {"label": "P 分型", "value": p_genotype},
                {"label": "组合分型", "value": subtype_combo},
                {"label": "最优参考株", "value": selected_reference},
            ],
            "quality_metrics": [
                {"label": "全长覆盖度", "value": full_length_coverage},
                {"label": "命中片段数", "value": segment_count},
                {"label": "覆盖碱基", "value": covered_bases},
                {"label": "参考总长度", "value": reference_bases},
                {"label": "支持 reads", "value": supported_reads},
                {"label": "coverage_sum", "value": coverage_sum},
            ],
            "sequence_name": sample_name or "-",
            "notes": "；".join([item for item in notes if item]),
            "mutation_table": {"status": "empty", "columns": [], "rows": []},
            "group_typing": rotavirus_group_typing,
            "subtype_typing": rotavirus_subtype_typing,
            "consensus_typing": rotavirus_consensus_typing,
            "igv": {"status": "empty"},
            "columns": display_columns,
            "rows": display_rows,
        }

    hpiv_summary = _read_tsv_rows(report_dir / f"{sample_name}_serotype_result.tsv")
    hpiv_columns = hpiv_summary.get("columns", [])
    hpiv_rows = hpiv_summary.get("rows", [])
    if any(str(column).strip() == "HPIV亚型" for column in hpiv_columns):
        first_summary_row = hpiv_rows[0] if hpiv_rows and isinstance(hpiv_rows[0], list) else []

        def _hpiv_summary_value(column_name: str, fallback: str = "-") -> str:
            if column_name not in hpiv_columns:
                return fallback
            index = hpiv_columns.index(column_name)
            value = str(first_summary_row[index] if index < len(first_summary_row) else "").strip()
            return value or fallback

        selection_path = report_dir / f"{sample_name}_hpiv_reference_selection" / "selection.tsv"
        coverage_summary_path = report_dir / "hpiv_coverage" / "hpiv.coverage.summary.tsv"
        selection_table = _read_tsv_rows(selection_path)
        coverage_summary_table = _read_tsv_rows(coverage_summary_path)
        selection_columns = selection_table.get("columns", [])
        selection_rows = selection_table.get("rows", [])
        coverage_summary_columns = coverage_summary_table.get("columns", [])
        coverage_summary_rows = coverage_summary_table.get("rows", [])
        first_selection = selection_rows[0] if selection_rows and isinstance(selection_rows[0], list) else []
        first_coverage_summary = coverage_summary_rows[0] if coverage_summary_rows and isinstance(coverage_summary_rows[0], list) else []

        display_selection_columns = [
            column
            for column in selection_columns
            if str(column).strip() not in {"reference_path", "gff_path"}
        ]
        display_selection_rows = []
        if display_selection_columns and selection_rows:
            selected_indexes = [selection_columns.index(column) for column in display_selection_columns]
            for row in selection_rows:
                if not isinstance(row, list):
                    continue
                display_selection_rows.append([row[index] if index < len(row) else "" for index in selected_indexes])

        def _hpiv_selection_value(column_name: str, fallback: str = "-") -> str:
            if column_name not in selection_columns:
                return fallback
            index = selection_columns.index(column_name)
            value = str(first_selection[index] if index < len(first_selection) else "").strip()
            return value or fallback

        def _hpiv_coverage_summary_value(column_name: str, fallback: str = "-") -> str:
            if column_name not in coverage_summary_columns:
                return fallback
            index = coverage_summary_columns.index(column_name)
            value = str(first_coverage_summary[index] if index < len(first_coverage_summary) else "").strip()
            return value or fallback

        selected_hpiv_type = _hpiv_summary_value("HPIV亚型", _hpiv_selection_value("hpiv_type"))
        selected_reference = _hpiv_summary_value("参考序列", Path(_hpiv_selection_value("reference_path", "")).name or "-")
        selected_gff = _hpiv_summary_value("注释文件", Path(_hpiv_selection_value("gff_path", "")).name or "-")
        coverage_value = _hpiv_summary_value("覆盖度", _hpiv_coverage_summary_value("coverage", _hpiv_selection_value("coverage")))
        mean_depth_value = _hpiv_summary_value("平均深度", _hpiv_coverage_summary_value("mean_depth", _hpiv_selection_value("mean_depth")))
        summary_note = _sanitize_virus_demo_note(_hpiv_summary_value("说明", ""))
        coverage_reference_name = _hpiv_coverage_summary_value("coverage_reference_name", "")
        notes: list[str] = []
        if selected_hpiv_type != "-":
            notes.append(f"自动选择参考型别：HPIV{selected_hpiv_type}")
        if selected_reference != "-":
            notes.append(f"参考序列：{selected_reference}")
        if selected_gff != "-":
            notes.append(f"注释文件：{selected_gff}")
        if coverage_reference_name:
            notes.append(f"覆盖度参考：{coverage_reference_name}")
        if summary_note != "-":
            notes.append(summary_note)
        mutation_table = _read_rsv_variant_annotation_table(report_dir)
        functional_annotation = _read_hpiv_functional_annotation_table(report_dir, mutation_table, selected_hpiv_type)
        return {
            "status": "ready",
            "mode": "hpiv_typing",
            "predicted_clade": selected_hpiv_type,
            "reference_name": selected_reference,
            "summary_cards": [
                {"label": "HPIV 亚型", "value": selected_hpiv_type},
                {"label": "参考序列", "value": selected_reference},
                {"label": "注释文件", "value": selected_gff},
                {"label": "覆盖度", "value": coverage_value},
            ],
            "quality_metrics": [
                {"label": "平均深度", "value": mean_depth_value},
            ] if mean_depth_value != "-" else [],
            "sequence_name": sample_name or "-",
            "notes": "；".join(notes),
            "mutation_table": mutation_table,
            "functional_annotation": functional_annotation,
            "igv": _discover_hpiv_igv_assets(report_dir),
            "columns": display_selection_columns or hpiv_columns,
            "rows": display_selection_rows or hpiv_rows[:1],
        }

    hiv_summary = _read_tsv_rows(report_dir / f"{sample_name}_serotype_result.tsv")
    hiv_columns = hiv_summary.get("columns", [])
    hiv_rows = hiv_summary.get("rows", [])
    hiv_resistance = _read_tsv_rows(report_dir / f"{sample_name}_hiv_resistance.tsv")
    has_hiv_assets = (
        any(str(column).strip() in {"NRTI最高等级", "NNRTI最高等级", "PI最高等级", "INSTI最高等级"} for column in hiv_columns)
        or bool(hiv_resistance.get("rows"))
    )
    if has_hiv_assets:
        first_summary_row = hiv_rows[0] if hiv_rows and isinstance(hiv_rows[0], list) else []

        def _hiv_summary_value(column_name: str, fallback: str = "-") -> str:
            if column_name not in hiv_columns:
                return fallback
            index = hiv_columns.index(column_name)
            value = str(first_summary_row[index] if index < len(first_summary_row) else "").strip()
            return value or fallback

        notes: list[str] = []
        for label in ("候选父本", "输入序列", "耐药输入序列", "序列告警", "说明"):
            value = _hiv_summary_value(label, "-")
            if value != "-":
                prefix = "" if label == "说明" else f"{label}："
                notes.append(f"{prefix}{value}")

        display_columns = hiv_columns or hiv_resistance.get("columns", [])
        display_rows = hiv_rows[:1] if hiv_rows else []
        broad_type = _hiv_summary_value("大亚型", _hiv_summary_value("病毒类型", "HIV"))
        subtype = _hiv_summary_value("子亚型", "-")
        recombination = _hiv_summary_value("重组判定", "-")
        representative_reference = _hiv_summary_value("代表株参考", "-")
        mutation_table = _read_rsv_variant_annotation_table(report_dir)
        hiv_reference_root = report_dir / f"{sample_name}_hiv_reference_selection"
        broad_typing = _read_tsv_rows(hiv_reference_root / "broad_typing.tsv")
        subtype_reference_typing = _read_tsv_rows(hiv_reference_root / "representative_subtype_typing.tsv")
        reference_selection = _read_tsv_rows(hiv_reference_root / "selection.tsv")
        subtyping_json_path = report_dir / f"{sample_name}_hiv_subtyping.json"
        subtyping_payload: dict[str, object] = {}
        if subtyping_json_path.is_file():
            try:
                subtyping_payload = json.loads(subtyping_json_path.read_text(encoding="utf-8", errors="ignore"))
            except (OSError, json.JSONDecodeError):
                subtyping_payload = {}
        resistance_json_path = report_dir / f"{sample_name}_hiv_resistance.json"
        resistance_payload: dict[str, object] = {}
        if resistance_json_path.is_file():
            try:
                raw_resistance_payload = json.loads(resistance_json_path.read_text(encoding="utf-8", errors="ignore"))
                if isinstance(raw_resistance_payload, dict):
                    resistance_payload = raw_resistance_payload
            except (OSError, json.JSONDecodeError):
                resistance_payload = {}
        subtyping_sample = {}
        if isinstance(subtyping_payload.get("samples"), list) and subtyping_payload["samples"]:
            candidate = subtyping_payload["samples"][0]
            if isinstance(candidate, dict):
                subtyping_sample = candidate
        resistance_sample = {}
        if isinstance(resistance_payload.get("samples"), list) and resistance_payload["samples"]:
            candidate = resistance_payload["samples"][0]
            if isinstance(candidate, dict):
                resistance_sample = candidate
        resistance_algorithm = resistance_payload.get("algorithm") if isinstance(resistance_payload.get("algorithm"), dict) else {}
        bootscan_assets: dict[str, str] = {}
        raw_bootscan_assets = subtyping_sample.get("bootscan_assets")
        if isinstance(raw_bootscan_assets, dict):
            for key, value in raw_bootscan_assets.items():
                path_text = str(value or "").strip()
                if not path_text:
                    continue
                candidate_path = Path(path_text)
                if not candidate_path.is_absolute():
                    candidate_path = report_dir / candidate_path
                try:
                    if candidate_path.is_file() and candidate_path.is_relative_to(report_dir):
                        bootscan_assets[key] = str(candidate_path.relative_to(report_dir))
                except OSError:
                    continue
        candidate_parent_parts: list[str] = []
        for item in subtyping_sample.get("candidate_parents", []) if isinstance(subtyping_sample, dict) else []:
            if not isinstance(item, dict):
                continue
            group = str(item.get("group") or "").strip()
            if not group:
                continue
            fraction = item.get("fraction")
            try:
                fraction_text = f"{float(fraction):.3f}"
            except (TypeError, ValueError):
                fraction_text = ""
            candidate_parent_parts.append(f"{group}{f' ({fraction_text})' if fraction_text else ''}")
        broad_coverage = "-"
        broad_depth = "-"
        if broad_typing.get("columns") and broad_typing.get("rows"):
            broad_columns = broad_typing["columns"]
            broad_row = broad_typing["rows"][0] if broad_typing["rows"] and isinstance(broad_typing["rows"][0], list) else []
            if "coverage" in broad_columns:
                idx = broad_columns.index("coverage")
                broad_coverage = str(broad_row[idx] if idx < len(broad_row) else "").strip() or "-"
            if "mean_depth" in broad_columns:
                idx = broad_columns.index("mean_depth")
                broad_depth = str(broad_row[idx] if idx < len(broad_row) else "").strip() or "-"
        subtype_coverage = "-"
        subtype_depth = "-"
        if subtype_reference_typing.get("columns") and subtype_reference_typing.get("rows"):
            subtype_columns = subtype_reference_typing["columns"]
            subtype_row = subtype_reference_typing["rows"][0] if subtype_reference_typing["rows"] and isinstance(subtype_reference_typing["rows"][0], list) else []
            if "coverage" in subtype_columns:
                idx = subtype_columns.index("coverage")
                subtype_coverage = str(subtype_row[idx] if idx < len(subtype_row) else "").strip() or "-"
            if "mean_depth" in subtype_columns:
                idx = subtype_columns.index("mean_depth")
                subtype_depth = str(subtype_row[idx] if idx < len(subtype_row) else "").strip() or "-"
        knowledge_summary = _build_viral_serotype_knowledge_summary(
            str(Path(__file__).resolve().parents[1]),
            "Human immunodeficiency virus",
            {
                "mode": "hiv_resistance",
                "predicted_group": broad_type,
                "predicted_clade": subtype if subtype != "-" else broad_type,
            },
        )
        return {
            "status": "ready",
            "mode": "hiv_resistance",
            "predicted_group": broad_type,
            "predicted_clade": subtype if subtype != "-" else broad_type,
            "summary_cards": [
                {"label": "HIV 大亚型", "value": broad_type},
                {"label": "子亚型", "value": subtype},
                {"label": "重组判定", "value": recombination},
                {"label": "代表株参考", "value": representative_reference},
                {"label": "NRTI 最高等级", "value": _hiv_summary_value("NRTI最高等级", "-")},
                {"label": "NNRTI 最高等级", "value": _hiv_summary_value("NNRTI最高等级", "-")},
                {"label": "PI 最高等级", "value": _hiv_summary_value("PI最高等级", "-")},
                {"label": "INSTI 最高等级", "value": _hiv_summary_value("INSTI最高等级", "-")},
            ],
            "quality_metrics": [
                {"label": "大亚型覆盖度", "value": broad_coverage},
                {"label": "大亚型平均深度", "value": broad_depth},
                {"label": "子亚型覆盖度", "value": subtype_coverage},
                {"label": "子亚型平均深度", "value": subtype_depth},
            ],
            "sequence_name": sample_name or "-",
            "notes": "；".join(notes),
            "knowledge_summary": knowledge_summary,
            "mutation_panels": [
                {"label": "PR 突变", "value": _hiv_summary_value("PR突变", "-")},
                {"label": "RT 突变", "value": _hiv_summary_value("RT突变", "-")},
                {"label": "IN 突变", "value": _hiv_summary_value("IN突变", "-")},
                {"label": "候选父本", "value": "；".join(candidate_parent_parts) if candidate_parent_parts else _hiv_summary_value("候选父本", "-")},
            ],
            "reference_selection": reference_selection,
            "broad_typing": broad_typing,
            "subtype_reference_typing": subtype_reference_typing,
            "subtyping_summary": {
                "assignment_label": str(subtyping_sample.get("assignment_label") or recombination or "-"),
                "predicted_group": str(subtyping_sample.get("predicted_group") or broad_type or "-"),
                "predicted_clade": str(subtyping_sample.get("predicted_clade") or subtype or "-"),
                "pure_tree_support": str(((subtyping_sample.get("pure_tree") or {}) if isinstance(subtyping_sample, dict) else {}).get("support") or "-"),
                "overall_tree_support": str(((subtyping_sample.get("overall_tree") or {}) if isinstance(subtyping_sample, dict) else {}).get("support") or "-"),
                "candidate_parent_count": len(candidate_parent_parts),
            },
            "bootscan_assets": bootscan_assets,
            "resistance_table": hiv_resistance,
            "resistance_payload": {
                "algorithm": resistance_algorithm,
                "sample": resistance_sample,
            },
            "mutation_table": mutation_table,
            "igv": _discover_hiv_igv_assets(report_dir),
            "columns": display_columns,
            "rows": display_rows,
        }

    default_table = _read_tsv_rows(report_dir / f"{sample_name}_serotype_result.tsv")
    if not _is_bordetella_pertussis(checkm_info):
        return {
            "status": "ready" if default_table.get("rows") else "empty",
            "mode": "generic",
            "columns": default_table.get("columns", []),
            "rows": default_table.get("rows", []),
        }
    rrn_path = report_dir / f"{sample_name}.2037.tsv"
    scheme_table = _read_tsv_rows(report_dir / f"{sample_name}_scheme.tsv")
    scheme_summary = ""
    if scheme_table.get("columns") and "分型" in scheme_table["columns"]:
        type_index = scheme_table["columns"].index("分型")
        parts: list[str] = []
        for row in scheme_table.get("rows", []):
            value = str((row[type_index] if type_index < len(row) else "") or "").strip()
            if value and value != "-" and value not in parts:
                parts.append(value)
        scheme_summary = "/".join(parts)
    notes = [
        {
            "label": "rrn A2037G",
            "state": "success" if rrn_path.is_file() else "neutral",
            "text": "鉴定到了rrn A2037G" if rrn_path.is_file() else "未鉴定到rrn A2037G",
        },
        {
            "label": "抗原基因",
            "state": "success" if scheme_table.get("rows") else "neutral",
            "text": "检测到抗原基因" if scheme_table.get("rows") else "未检测到抗原基因",
        },
    ]
    status = "ready" if rrn_path.is_file() or scheme_table.get("rows") else "empty"
    return {
        "status": status,
        "mode": "bordetella_pertussis",
        "columns": scheme_table.get("columns", []),
        "rows": scheme_table.get("rows", []),
        "notes": notes,
        "scheme_summary": scheme_summary,
    }
