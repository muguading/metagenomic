from __future__ import annotations

import csv
import json
import re
from collections import Counter, defaultdict
from difflib import SequenceMatcher
from pathlib import Path

from .parse_utils import _safe_float, _safe_int
from .runtime_paths import _resolve_runtime_database_root
from .table_io import _read_tsv_rows

def _read_first_fasta_sequence(path: Path) -> str:
    if not path.is_file():
        return ""
    sequence_chunks: list[str] = []
    seen_header = False
    try:
        with path.open("r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                text = line.strip()
                if not text:
                    continue
                if text.startswith(">"):
                    if seen_header:
                        break
                    seen_header = True
                    continue
                if seen_header:
                    sequence_chunks.append(text)
    except OSError:
        return ""
    return "".join(sequence_chunks).upper()

def _read_named_fasta_sequence(path: Path, record_name: str) -> str:
    target = str(record_name or "").strip()
    if not path.is_file() or not target:
        return ""
    collecting = False
    sequence_chunks: list[str] = []
    try:
        with path.open("r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                text = line.strip()
                if not text:
                    continue
                if text.startswith(">"):
                    current_header = text[1:].strip()
                    current_id = current_header.split()[0]
                    if collecting:
                        break
                    collecting = target in {current_header, current_id}
                    continue
                if collecting:
                    sequence_chunks.append(text)
    except OSError:
        return ""
    return "".join(sequence_chunks).upper()

def _estimate_sequence_offset(reference_sequence: str, sample_sequence: str, kmer_size: int = 15) -> int | None:
    ref_seq = str(reference_sequence or "").upper()
    sample_seq = str(sample_sequence or "").upper()
    if len(ref_seq) < kmer_size or len(sample_seq) < kmer_size:
        return None
    sample_index: dict[str, list[int]] = defaultdict(list)
    for sample_pos in range(0, len(sample_seq) - kmer_size + 1):
        kmer = sample_seq[sample_pos:sample_pos + kmer_size]
        if "N" in kmer:
            continue
        if len(sample_index[kmer]) < 32:
            sample_index[kmer].append(sample_pos)
    offset_counter: Counter[int] = Counter()
    step = max(1, kmer_size // 2)
    for ref_pos in range(0, len(ref_seq) - kmer_size + 1, step):
        kmer = ref_seq[ref_pos:ref_pos + kmer_size]
        if "N" in kmer:
            continue
        for sample_pos in sample_index.get(kmer, []):
            offset_counter[sample_pos - ref_pos] += 1
    if not offset_counter:
        return None
    return offset_counter.most_common(1)[0][0]

def _count_substitutions_between_sequences(reference_sequence: str, sample_sequence: str) -> int:
    ref_seq = str(reference_sequence or "").upper()
    query_seq = str(sample_sequence or "").upper()
    if not ref_seq or not query_seq:
        return 0
    matcher = SequenceMatcher(None, ref_seq, query_seq, autojunk=False)
    snp_count = 0
    for tag, i1, i2, j1, j2 in matcher.get_opcodes():
        if tag != "replace":
            continue
        ref_slice = ref_seq[i1:i2]
        query_slice = query_seq[j1:j2]
        snp_count += sum(1 for left, right in zip(ref_slice, query_slice) if left != right)
    return snp_count

def _list_substitutions_between_sequences(reference_sequence: str, sample_sequence: str, sample_offset: int = 0) -> list[dict[str, object]]:
    ref_seq = str(reference_sequence or "").upper()
    query_seq = str(sample_sequence or "").upper()
    if not ref_seq or not query_seq:
        return []
    matcher = SequenceMatcher(None, ref_seq, query_seq, autojunk=False)
    results: list[dict[str, object]] = []
    for tag, i1, i2, j1, j2 in matcher.get_opcodes():
        if tag != "replace":
            continue
        ref_slice = ref_seq[i1:i2]
        query_slice = query_seq[j1:j2]
        for index, (left, right) in enumerate(zip(ref_slice, query_slice)):
            if left == right:
                continue
            ref_pos = i1 + index + 1
            sample_pos = sample_offset + j1 + index + 1
            results.append({
                "ref_pos": ref_pos,
                "sample_pos": sample_pos,
                "ref": left,
                "alt": right,
            })
    return results

def _build_hadv_mutation_lookup(mutation_table: dict) -> dict[int, list[dict[str, object]]]:
    columns = mutation_table.get("columns") if isinstance(mutation_table, dict) else []
    rows = mutation_table.get("rows") if isinstance(mutation_table, dict) else []
    if not isinstance(columns, list) or not isinstance(rows, list):
        return {}
    pos_idx = columns.index("位置") if "位置" in columns else -1
    ref_idx = columns.index("参考碱基") if "参考碱基" in columns else -1
    alt_idx = columns.index("突变碱基") if "突变碱基" in columns else -1
    depth_idx = columns.index("测序深度 / 突变频率") if "测序深度 / 突变频率" in columns else -1
    hgvs_c_idx = columns.index("HGVS.c") if "HGVS.c" in columns else -1
    hgvs_p_idx = columns.index("HGVS.p") if "HGVS.p" in columns else -1
    quality_idx = columns.index("质量分层") if "质量分层" in columns else -1
    lookup: dict[int, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        if not isinstance(row, list) or pos_idx < 0 or pos_idx >= len(row):
            continue
        try:
            pos_value = int(row[pos_idx])
        except (TypeError, ValueError):
            continue
        lookup[pos_value].append({
            "ref": str(row[ref_idx] if ref_idx >= 0 and ref_idx < len(row) else "").strip(),
            "alt": str(row[alt_idx] if alt_idx >= 0 and alt_idx < len(row) else "").strip(),
            "depth": str(row[depth_idx] if depth_idx >= 0 and depth_idx < len(row) else "").strip(),
            "hgvs_c": str(row[hgvs_c_idx] if hgvs_c_idx >= 0 and hgvs_c_idx < len(row) else "").strip(),
            "hgvs_p": str(row[hgvs_p_idx] if hgvs_p_idx >= 0 and hgvs_p_idx < len(row) else "").strip(),
            "quality": str(row[quality_idx] if quality_idx >= 0 and quality_idx < len(row) else "").strip(),
        })
    return lookup

def _read_hadv_phf_snp_section(report_dir: Path, sample_name: str, mutation_table: dict | None = None) -> dict:
    phf_dir = report_dir / f"{sample_name}_hadv_reference_selection" / "phf_typing"
    phf_table = _read_tsv_rows(phf_dir / "phf_typing.tsv")
    phf_columns = phf_table.get("columns", [])
    phf_rows = phf_table.get("rows", [])
    subject_index = phf_columns.index("subject") if "subject" in phf_columns else -1
    type_index = phf_columns.index("matched_type") if "matched_type" in phf_columns else -1
    sample_consensus = _read_first_fasta_sequence(report_dir / f"{sample_name}.consensus.fasta") or _read_first_fasta_sequence(report_dir / f"{sample_name}.final.fasta")
    if not sample_consensus:
        return {"status": "empty", "columns": [], "rows": []}
    db_root = Path(__file__).resolve().parent.parent / "database" / "virus" / "hadv"
    db_paths = {
        "penton": db_root / "blastn_db_penton" / "hadv_types_ref_penton.fa",
        "hexon": db_root / "blastn_db_hexon" / "hadv_types_ref_hexon.fa",
        "fiber": db_root / "blastn_db_fiber" / "hadv_types_ref_fiber.fa",
    }
    mutation_lookup = _build_hadv_mutation_lookup(mutation_table or {})
    rows: list[list[object]] = []
    values: list[int] = []
    x_values: list[str] = []
    detail_rows: list[list[object]] = []
    for row in phf_rows:
        if not isinstance(row, list) or len(row) < 2:
            continue
        gene_name = str(row[1] or "").strip().lower()
        if gene_name not in db_paths:
            continue
        subject = str(row[subject_index] if subject_index >= 0 and subject_index < len(row) else "").strip()
        matched_type = str(row[type_index] if type_index >= 0 and type_index < len(row) else "").strip()
        reference_sequence = _read_named_fasta_sequence(db_paths[gene_name], subject)
        if not reference_sequence:
            continue
        estimated_offset = _estimate_sequence_offset(reference_sequence, sample_consensus)
        if estimated_offset is None:
            continue
        start = max(0, estimated_offset - 80)
        end = min(len(sample_consensus), estimated_offset + len(reference_sequence) + 80)
        sample_window = sample_consensus[start:end]
        snp_count = _count_substitutions_between_sequences(reference_sequence, sample_window)
        substitutions = _list_substitutions_between_sequences(reference_sequence, sample_window, sample_offset=start)
        display_name = gene_name.capitalize()
        rows.append([display_name, matched_type or "-", subject or "-", int(snp_count)])
        values.append(int(snp_count))
        x_values.append(display_name)
        for item in substitutions:
            sample_pos = int(item.get("sample_pos") or 0)
            ref_base = str(item.get("ref") or "")
            alt_base = str(item.get("alt") or "")
            matching_mutations = mutation_lookup.get(sample_pos, [])
            matched_meta = next(
                (
                    meta for meta in matching_mutations
                    if str(meta.get("ref") or "").upper() == ref_base.upper()
                    and str(meta.get("alt") or "").upper() == alt_base.upper()
                ),
                matching_mutations[0] if matching_mutations else {},
            )
            detail_rows.append([
                display_name,
                matched_type or "-",
                subject or "-",
                int(item.get("ref_pos") or 0),
                sample_pos,
                ref_base,
                alt_base,
                str(matched_meta.get("depth") or "-"),
                str(matched_meta.get("quality") or "-"),
                str(matched_meta.get("hgvs_c") or "-"),
                str(matched_meta.get("hgvs_p") or "-"),
            ])
    if not rows:
        return {"status": "empty", "columns": [], "rows": []}
    return {
        "status": "ready",
        "columns": ["基因", "命中分型", "命中参考", "差异SNP数"],
        "rows": rows,
        "values": values,
        "x_values": x_values,
        "detail_columns": ["基因", "命中分型", "命中参考", "参考基因位点", "样本坐标", "Ref", "Alt", "Depth/MAF", "质量分层", "HGVS.c", "HGVS.p"],
        "detail_rows": detail_rows,
    }

def _read_variant_annotation_table(json_path: Path, tsv_path: Path, source_vcf_name: str) -> dict:
    if json_path.is_file():
        try:
            payload = json.loads(json_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            payload = {}
        rows = payload.get("rows") if isinstance(payload.get("rows"), list) else []
        if rows and isinstance(rows[0], dict):
            preferred = [
                "核苷酸突变",
                "变异类型",
                "质量分层",
                "位置",
                "参考碱基",
                "突变碱基",
                "质量值QUAL",
                "测序深度DP",
                "突变频率MAF",
                "注释效应",
                "影响等级",
                "基因",
                "HGVS.c",
                "HGVS.p",
                "氨基酸位点",
                "警告信息",
            ]
            columns = [column for column in preferred if any(column in row for row in rows)]
            seen = set(columns)
            for row in rows:
                for key in row.keys():
                    if key not in seen:
                        columns.append(key)
                        seen.add(key)
            list_rows = [[row.get(column, "") for column in columns] for row in rows]
            return {
                "status": str(payload.get("status") or "ready"),
                "columns": columns,
                "rows": list_rows,
                "total_variants": int(payload.get("total_variants") or len(rows)),
                "high_quality_variants": int(payload.get("high_quality_variants") or 0),
                "low_quality_variants": int(payload.get("low_quality_variants") or 0),
                "source_vcf": str(payload.get("source_vcf") or ""),
            }
    raw = _read_tsv_rows(tsv_path)
    columns = raw.get("columns", [])
    rows = raw.get("rows", [])
    quality_index = columns.index("质量分层") if "质量分层" in columns else -1
    high_quality = 0
    low_quality = 0
    if quality_index >= 0:
        for row in rows:
            value = str(row[quality_index] if quality_index < len(row) else "").strip()
            if value == "高质量突变":
                high_quality += 1
            elif value == "低质量突变":
                low_quality += 1
    return {
        "status": "ready" if rows else "empty",
        "columns": columns,
        "rows": rows,
        "total_variants": len(rows),
        "high_quality_variants": high_quality,
        "low_quality_variants": low_quality,
        "source_vcf": source_vcf_name,
    }

def _read_simple_annotated_variant_table(report_dir: Path, *, prefer_filtered: bool = False) -> dict:
    annotated_vcf_path = report_dir / "snps.anno.vcf"
    preferred_vcf_path = report_dir / ("snps.filt1.vcf" if prefer_filtered else "snps.raw.vcf")
    fallback_vcf_path = report_dir / ("snps.raw.vcf" if prefer_filtered else "snps.filt1.vcf")
    raw_vcf_path = preferred_vcf_path if preferred_vcf_path.is_file() else fallback_vcf_path
    if not raw_vcf_path.is_file():
        return {"status": "empty", "columns": [], "rows": [], "total_variants": 0, "source_vcf": ""}

    def _parse_kv_pairs(text: str) -> dict[str, str]:
        result: dict[str, str] = {}
        for item in str(text or "").split(";"):
            if not item:
                continue
            if "=" not in item:
                result[item] = ""
                continue
            key, value = item.split("=", 1)
            result[str(key).strip()] = str(value).strip()
        return result

    def _choose_ann_entry(entries: list[list[str]]) -> list[str]:
        if not entries:
            return []

        def _ann_score(parts: list[str]) -> tuple[int, int, int]:
            effect = str(parts[1] if len(parts) > 1 else "").strip().lower()
            impact = str(parts[2] if len(parts) > 2 else "").strip().upper()
            has_gene = 1 if str(parts[3] if len(parts) > 3 else "").strip() else 0
            impact_rank = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1}.get(impact, 0)
            non_modifier = 1 if effect and effect != "intergenic_region" and "upstream" not in effect and "downstream" not in effect else 0
            return (non_modifier, impact_rank, has_gene)

        return max(entries, key=_ann_score)

    annotated_lookup: dict[tuple[str, str, str, str], dict[str, str]] = {}
    if annotated_vcf_path.is_file():
        try:
            with annotated_vcf_path.open("r", encoding="utf-8", errors="ignore") as handle:
                for line in handle:
                    if not line or line.startswith("#"):
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 8:
                        continue
                    chrom, pos, _vid, ref, alt, _qual, _flt, info = parts[:8]
                    alt_value = str(alt or "").split(",")[0].strip()
                    info_map = _parse_kv_pairs(info)
                    ann_entries = []
                    ann_raw = str(info_map.get("ANN") or "").strip()
                    if ann_raw:
                        ann_entries = [entry.split("|") for entry in ann_raw.split(",") if entry]
                    best_ann = _choose_ann_entry(ann_entries)
                    annotated_lookup[(chrom, pos, ref, alt_value)] = {
                        "gene_name": str(best_ann[3] if len(best_ann) > 3 else "").strip() or chrom,
                        "effect": str(best_ann[1] if len(best_ann) > 1 else "").strip(),
                        "impact": str(best_ann[2] if len(best_ann) > 2 else "").strip(),
                        "hgvs_c": str(best_ann[9] if len(best_ann) > 9 else "").strip(),
                        "hgvs_p": str(best_ann[10] if len(best_ann) > 10 else "").strip(),
                    }
        except OSError:
            annotated_lookup = {}

    rows: list[list[object]] = []
    high_quality = 0
    low_quality = 0
    try:
        with raw_vcf_path.open("r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if not line or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 10:
                    continue
                chrom, pos, _vid, ref, alt, qual_raw, _filt, info, fmt, sample = parts[:10]
                alt_value = str(alt or "").split(",")[0].strip()
                info_map = _parse_kv_pairs(info)
                ann_info = annotated_lookup.get((chrom, pos, ref, alt_value), {})
                sample_fields = str(fmt or "").split(":")
                sample_values = str(sample or "").split(":")
                sample_map = {sample_fields[i]: sample_values[i] if i < len(sample_values) else "" for i in range(len(sample_fields))}

                qual_value = _safe_float(qual_raw) or 0.0
                depth_value = _safe_int(info_map.get("DP"))
                if depth_value is None:
                    depth_value = _safe_int(sample_map.get("DP")) or 0
                ao_value = _safe_float(str(info_map.get("AO") or "").split(",")[0].strip())
                if ao_value is None:
                    ao_value = _safe_float(str(sample_map.get("AO") or "").split(",")[0].strip()) or 0.0
                ro_value = _safe_float(str(info_map.get("RO") or "").split(",")[0].strip())
                if ro_value is None:
                    ro_value = _safe_float(str(sample_map.get("RO") or "").split(",")[0].strip()) or 0.0
                allele_depth = float(ao_value) + float(ro_value)
                if allele_depth <= 0 and depth_value > 0:
                    allele_depth = float(depth_value)
                maf_value = round((float(ao_value) / allele_depth), 6) if allele_depth > 0 else 0.0

                quality_label = "高质量突变" if qual_value > 10 and depth_value > 10 and maf_value > 0.1 else "低质量突变"
                if prefer_filtered and raw_vcf_path.name == "snps.filt1.vcf" and quality_label != "高质量突变":
                    quality_label = "高质量突变"
                if quality_label == "高质量突变":
                    high_quality += 1
                else:
                    low_quality += 1

                rows.append([
                    str(ann_info.get("gene_name") or chrom).strip() or chrom,
                    pos,
                    ref,
                    alt_value,
                    str(depth_value),
                    f"{maf_value:.4f}",
                    str(ann_info.get("hgvs_c") or "").strip(),
                    str(ann_info.get("hgvs_p") or "").strip(),
                    chrom,
                    quality_label,
                    str(ann_info.get("impact") or "").strip(),
                    str(ann_info.get("effect") or "").strip(),
                    "插入" if len(alt_value) > len(ref) else ("缺失" if len(ref) > len(alt_value) else "替换"),
                ])
    except OSError:
        return {"status": "empty", "columns": [], "rows": [], "total_variants": 0, "source_vcf": ""}

    columns = [
        "基因名",
        "位置",
        "参考碱基",
        "突变碱基",
        "测序深度DP",
        "突变频率MAF",
        "HGVS.c",
        "HGVS.p",
        "染色体",
        "质量分层",
        "影响等级",
        "注释效应",
        "变异类型",
    ]
    return {
        "status": "ready" if rows else "empty",
        "columns": columns,
        "rows": rows,
        "total_variants": len(rows),
        "high_quality_variants": high_quality,
        "low_quality_variants": low_quality,
        "source_vcf": str(report_dir / "snps.anno.vcf"),
    }

def _read_ncov_variant_annotation_table(report_dir: Path) -> dict:
    return _read_variant_annotation_table(
        report_dir / "snps.raw.mutation_table.json",
        report_dir / "snps.raw.mutation_table.tsv",
        str(report_dir / "snps.raw.ann.vcf"),
    )

def _read_influenza_variant_annotation_table(report_dir: Path) -> dict:
    return _read_variant_annotation_table(
        report_dir / "snps.filt1.mutation_table.json",
        report_dir / "snps.filt1.mutation_table.tsv",
        str(report_dir / "snps.anno.vcf"),
    )

def _build_influenza_mutation_display_table(report_dir: Path, fallback_table: dict) -> dict:
    json_path = report_dir / "snps.filt1.mutation_table.json"
    if json_path.is_file():
        try:
            payload = json.loads(json_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            payload = {}
        raw_rows = payload.get("rows") if isinstance(payload.get("rows"), list) else []
        if raw_rows and isinstance(raw_rows[0], dict):
            columns = [
                "基因名",
                "位置",
                "参考碱基",
                "突变碱基",
                "测序深度 / 突变频率",
                "HGVS.c",
                "HGVS.p",
                "染色体",
                "质量分层",
                "影响等级",
                "注释效应",
                "变异类型",
            ]
            rows: list[list[object]] = []
            for row in raw_rows:
                gene_name = str(row.get("染色体") or "-").strip() or "-"
                depth_value = str(row.get("测序深度DP") or "").strip()
                maf_raw = row.get("突变频率MAF")
                maf_value = ""
                try:
                    if maf_raw not in (None, ""):
                        maf_value = f"{float(maf_raw):.4f}"
                except (TypeError, ValueError):
                    maf_value = str(maf_raw or "").strip()
                rows.append([
                    gene_name,
                    row.get("位置", ""),
                    str(row.get("参考碱基") or "").strip(),
                    str(row.get("突变碱基") or "").strip(),
                    " / ".join([part for part in [depth_value, maf_value] if part]) or "-",
                    str(row.get("HGVS.c") or "").strip(),
                    str(row.get("HGVS.p") or "").strip(),
                    str(row.get("染色体") or "").strip(),
                    str(row.get("质量分层") or "").strip(),
                    str(row.get("影响等级") or "").strip(),
                    str(row.get("注释效应") or "").strip(),
                    str(row.get("变异类型") or "").strip(),
                ])
            return {"columns": columns, "rows": rows}
    fallback_columns = fallback_table.get("columns") or []
    fallback_rows = fallback_table.get("rows") or []
    if not fallback_columns or not fallback_rows:
        return fallback_table
    display_columns = [("基因名" if str(column) == "染色体" else str(column)) for column in fallback_columns]
    return {"columns": display_columns, "rows": fallback_rows}

def _build_ncov_mutation_display_table(report_dir: Path, fallback_table: dict) -> dict:
    json_path = report_dir / "snps.raw.mutation_table.json"
    if json_path.is_file():
        try:
            payload = json.loads(json_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            payload = {}
        raw_rows = payload.get("rows") if isinstance(payload.get("rows"), list) else []
        if raw_rows and isinstance(raw_rows[0], dict):
            columns = [
                "基因名",
                "位置",
                "参考碱基",
                "突变碱基",
                "测序深度 / 突变频率",
                "HGVS.c",
                "HGVS.p",
                "染色体",
                "质量分层",
                "影响等级",
                "注释效应",
                "变异类型",
            ]
            rows: list[list[object]] = []
            for row in raw_rows:
                gene_name = str(row.get("染色体") or "-").strip() or "-"
                depth_value = str(row.get("测序深度DP") or "").strip()
                maf_raw = row.get("突变频率MAF")
                maf_value = ""
                try:
                    if maf_raw not in (None, ""):
                        maf_value = f"{float(maf_raw):.4f}"
                except (TypeError, ValueError):
                    maf_value = str(maf_raw or "").strip()
                rows.append([
                    gene_name,
                    row.get("位置", ""),
                    str(row.get("参考碱基") or "").strip(),
                    str(row.get("突变碱基") or "").strip(),
                    " / ".join([part for part in [depth_value, maf_value] if part]) or "-",
                    str(row.get("HGVS.c") or "").strip(),
                    str(row.get("HGVS.p") or "").strip(),
                    str(row.get("染色体") or "").strip(),
                    str(row.get("质量分层") or "").strip(),
                    str(row.get("影响等级") or "").strip(),
                    str(row.get("注释效应") or "").strip(),
                    str(row.get("变异类型") or "").strip(),
                ])
            return {
                "status": str(payload.get("status") or "ready"),
                "columns": columns,
                "rows": rows,
                "total_variants": int(payload.get("total_variants") or len(rows)),
                "high_quality_variants": int(payload.get("high_quality_variants") or 0),
                "low_quality_variants": int(payload.get("low_quality_variants") or 0),
                "source_vcf": str(payload.get("source_vcf") or ""),
            }
    fallback_columns = fallback_table.get("columns") or []
    fallback_rows = fallback_table.get("rows") or []
    if not fallback_columns or not fallback_rows:
        return fallback_table
    display_columns = [("基因名" if str(column) == "染色体" else str(column)) for column in fallback_columns]
    return {
        "status": str(fallback_table.get("status") or ("ready" if fallback_rows else "empty")),
        "columns": display_columns,
        "rows": fallback_rows,
        "total_variants": int(fallback_table.get("total_variants") or len(fallback_rows)),
        "high_quality_variants": int(fallback_table.get("high_quality_variants") or 0),
        "low_quality_variants": int(fallback_table.get("low_quality_variants") or 0),
        "source_vcf": str(fallback_table.get("source_vcf") or ""),
    }

def _read_monkeypox_variant_annotation_table(report_dir: Path) -> dict:
    annotated_vcf_path = report_dir / "snps.anno.vcf"
    raw_vcf_path = report_dir / "snps.raw.vcf"
    if not raw_vcf_path.is_file():
        return {"status": "empty", "columns": [], "rows": [], "total_variants": 0, "source_vcf": ""}
    rows: list[dict[str, object]] = []

    def _parse_kv_pairs(text: str) -> dict[str, str]:
        result: dict[str, str] = {}
        for item in str(text or "").split(";"):
            if not item:
                continue
            if "=" not in item:
                result[item] = ""
                continue
            key, value = item.split("=", 1)
            result[str(key).strip()] = str(value).strip()
        return result

    def _choose_ann_entry(entries: list[list[str]]) -> list[str]:
        if not entries:
            return []
        def _ann_score(parts: list[str]) -> tuple[int, int, int]:
            effect = str(parts[1] if len(parts) > 1 else "").strip().lower()
            impact = str(parts[2] if len(parts) > 2 else "").strip().upper()
            has_gene = 1 if str(parts[3] if len(parts) > 3 else "").strip() else 0
            impact_rank = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1}.get(impact, 0)
            non_modifier = 1 if effect and effect != "intergenic_region" and "upstream" not in effect and "downstream" not in effect else 0
            return (non_modifier, impact_rank, has_gene)
        return max(entries, key=_ann_score)

    def _max_homopolymer_run(sequence: str | None) -> int:
        text = str(sequence or "").strip().upper()
        if not text:
            return 0
        best = 1
        current = 1
        for index in range(1, len(text)):
            if text[index] == text[index - 1]:
                current += 1
                best = max(best, current)
            else:
                current = 1
        return best

    def _is_long_poly_indel(ref: str | None, alt: str | None) -> bool:
        ref_text = str(ref or "").strip().upper()
        alt_text = str(alt or "").strip().upper()
        if not ref_text or not alt_text or len(ref_text) == len(alt_text):
            return False
        return max(_max_homopolymer_run(ref_text), _max_homopolymer_run(alt_text)) > 6

    annotated_lookup: dict[tuple[str, str, str, str], dict[str, str]] = {}
    if annotated_vcf_path.is_file():
        try:
            with annotated_vcf_path.open("r", encoding="utf-8", errors="ignore") as handle:
                for line in handle:
                    if not line or line.startswith("#"):
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 8:
                        continue
                    chrom, pos, _vid, ref, alt, _qual, _flt, info = parts[:8]
                    alt_value = str(alt or "").split(",")[0].strip()
                    info_map = _parse_kv_pairs(info)
                    ann_entries = []
                    ann_raw = str(info_map.get("ANN") or "").strip()
                    if ann_raw:
                        ann_entries = [entry.split("|") for entry in ann_raw.split(",") if entry]
                    best_ann = _choose_ann_entry(ann_entries)
                    annotated_lookup[(chrom, pos, ref, alt_value)] = {
                        "gene_name": str(best_ann[3] if len(best_ann) > 3 else "").strip() or chrom,
                        "effect": str(best_ann[1] if len(best_ann) > 1 else "").strip(),
                        "impact": str(best_ann[2] if len(best_ann) > 2 else "").strip(),
                        "hgvs_c": str(best_ann[9] if len(best_ann) > 9 else "").strip(),
                        "hgvs_p": str(best_ann[10] if len(best_ann) > 10 else "").strip(),
                    }
        except OSError:
            annotated_lookup = {}

    try:
        with raw_vcf_path.open("r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if not line or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 10:
                    continue
                chrom, pos, _vid, ref, alt, qual_raw, filt, info, fmt, sample = parts[:10]
                info_map = _parse_kv_pairs(info)
                alt_value = str(alt or "").split(",")[0].strip()
                ann_info = annotated_lookup.get((chrom, pos, ref, alt_value), {})
                gene_name = str(ann_info.get("gene_name") or chrom).strip() or chrom
                effect = str(ann_info.get("effect") or "").strip()
                impact = str(ann_info.get("impact") or "").strip()
                hgvs_c = str(ann_info.get("hgvs_c") or "").strip()
                hgvs_p = str(ann_info.get("hgvs_p") or "").strip()
                qual = _safe_float(qual_raw) or 0.0
                depth_value = _safe_int(info_map.get("DP")) or 0
                ao_value = _safe_float(str(info_map.get("AO") or "").split(",")[0].strip()) or 0.0
                ro_value = _safe_float(str(info_map.get("RO") or "").split(",")[0].strip()) or 0.0
                sample_fields = str(fmt or "").split(":")
                sample_values = str(sample or "").split(":")
                sample_map = {sample_fields[i]: sample_values[i] if i < len(sample_values) else "" for i in range(len(sample_fields))}
                if depth_value <= 0:
                    depth_value = _safe_int(sample_map.get("DP")) or 0
                if ao_value <= 0:
                    ao_value = _safe_float(str(sample_map.get("AO") or "").split(",")[0].strip()) or 0.0
                if ro_value <= 0:
                    ro_value = _safe_float(str(sample_map.get("RO") or "").split(",")[0].strip()) or 0.0
                maf_numeric = 0.0
                try:
                    allele_depth = float(ao_value) + float(ro_value)
                    if allele_depth <= 0 and depth_value > 0:
                        allele_depth = float(depth_value)
                    maf_numeric = round((float(ao_value) / allele_depth), 6) if allele_depth > 0 else 0.0
                except (TypeError, ValueError, ZeroDivisionError):
                    maf_numeric = 0.0
                is_poly_indel = _is_long_poly_indel(ref, alt_value)
                if is_poly_indel:
                    quality_label = "低质量突变"
                else:
                    quality_label = "高质量突变" if qual > 10 and depth_value > 10 and maf_numeric > 0.1 else "低质量突变"
                rows.append(
                    {
                        "基因名": gene_name or chrom,
                        "位置": int(pos) if str(pos).isdigit() else pos,
                        "参考碱基": ref,
                        "突变碱基": alt_value,
                        "测序深度 / 突变频率": " / ".join([part for part in [str(depth_value or ""), f"{maf_numeric:.4f}" if depth_value > 0 else ""] if part]) or "-",
                        "HGVS.c": hgvs_c,
                        "HGVS.p": hgvs_p,
                        "染色体": chrom,
                        "质量分层": quality_label,
                        "质量值QUAL": round(qual, 2),
                        "测序深度DP": depth_value,
                        "突变频率MAF": round(maf_numeric, 6),
                        "Poly位点": "是" if is_poly_indel else "否",
                        "影响等级": impact,
                        "注释效应": effect,
                        "变异类型": str(info_map.get("TYPE") or "").split(",")[0].strip(),
                        "FILTER": filt,
                    }
                )
    except OSError:
        return {"status": "empty", "columns": [], "rows": [], "total_variants": 0, "source_vcf": ""}
    columns = [
        "基因名",
        "位置",
        "参考碱基",
        "突变碱基",
        "测序深度 / 突变频率",
        "HGVS.c",
        "HGVS.p",
        "染色体",
        "质量分层",
        "Poly位点",
        "影响等级",
        "注释效应",
        "变异类型",
    ]
    list_rows = [[row.get(column, "") for column in columns] for row in rows]
    high_quality = sum(1 for row in rows if str(row.get("质量分层") or "").strip() == "高质量突变")
    return {
        "status": "ready" if rows else "empty",
        "columns": columns,
        "rows": list_rows,
        "total_variants": len(rows),
        "high_quality_variants": high_quality,
        "low_quality_variants": len(rows) - high_quality,
        "source_vcf": str(raw_vcf_path),
    }

def _read_rsv_variant_annotation_table(report_dir: Path) -> dict:
    table = _read_variant_annotation_table(
        report_dir / "snps.filt1.mutation_table.json",
        report_dir / "snps.filt1.mutation_table.tsv",
        str(report_dir / "snps.anno.vcf"),
    )
    table_rows = table.get("rows")
    if isinstance(table_rows, list) and table_rows:
        return table

    raw_vcf_path = report_dir / "snps.raw.vcf"
    filt_vcf_path = report_dir / "snps.filt1.vcf"
    annotated_vcf_path = report_dir / "snps.anno.vcf"
    if not raw_vcf_path.is_file():
        fallback = _read_simple_annotated_variant_table(report_dir, prefer_filtered=True)
        return fallback if fallback.get("status") == "ready" else table
    rows: list[dict[str, object]] = []

    def _parse_kv_pairs(text: str) -> dict[str, str]:
        result: dict[str, str] = {}
        for item in str(text or "").split(";"):
            if not item:
                continue
            if "=" not in item:
                result[item] = ""
                continue
            key, value = item.split("=", 1)
            result[str(key).strip()] = str(value).strip()
        return result

    def _choose_ann_entry(entries: list[list[str]]) -> list[str]:
        if not entries:
            return []

        def _ann_score(parts: list[str]) -> tuple[int, int, int]:
            effect = str(parts[1] if len(parts) > 1 else "").strip().lower()
            impact = str(parts[2] if len(parts) > 2 else "").strip().upper()
            has_gene = 1 if str(parts[3] if len(parts) > 3 else "").strip() else 0
            impact_rank = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1}.get(impact, 0)
            non_modifier = 1 if effect and effect != "intergenic_region" and "upstream" not in effect and "downstream" not in effect else 0
            return (non_modifier, impact_rank, has_gene)

        return max(entries, key=_ann_score)

    high_quality_keys: set[tuple[str, str, str, str]] = set()
    if filt_vcf_path.is_file():
        try:
            with filt_vcf_path.open("r", encoding="utf-8", errors="ignore") as handle:
                for line in handle:
                    if not line or line.startswith("#"):
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 5:
                        continue
                    chrom, pos, _vid, ref, alt = parts[:5]
                    high_quality_keys.add((chrom, pos, ref, str(alt or "").split(",")[0].strip()))
        except OSError:
            high_quality_keys = set()

    annotated_lookup: dict[tuple[str, str, str, str], dict[str, str]] = {}
    if annotated_vcf_path.is_file():
        try:
            with annotated_vcf_path.open("r", encoding="utf-8", errors="ignore") as handle:
                for line in handle:
                    if not line or line.startswith("#"):
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 8:
                        continue
                    chrom, pos, _vid, ref, alt, _qual, _flt, info = parts[:8]
                    alt_value = str(alt or "").split(",")[0].strip()
                    info_map = _parse_kv_pairs(info)
                    ann_entries = []
                    ann_raw = str(info_map.get("ANN") or "").strip()
                    if ann_raw:
                        ann_entries = [entry.split("|") for entry in ann_raw.split(",") if entry]
                    best_ann = _choose_ann_entry(ann_entries)
                    annotated_lookup[(chrom, pos, ref, alt_value)] = {
                        "gene_name": str(best_ann[3] if len(best_ann) > 3 else "").strip() or chrom,
                        "effect": str(best_ann[1] if len(best_ann) > 1 else "").strip(),
                        "impact": str(best_ann[2] if len(best_ann) > 2 else "").strip(),
                        "hgvs_c": str(best_ann[9] if len(best_ann) > 9 else "").strip(),
                        "hgvs_p": str(best_ann[10] if len(best_ann) > 10 else "").strip(),
                    }
        except OSError:
            annotated_lookup = {}

    try:
        with raw_vcf_path.open("r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if not line or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 10:
                    continue
                chrom, pos, _vid, ref, alt, qual_raw, filt, info, fmt, sample = parts[:10]
                alt_value = str(alt or "").split(",")[0].strip()
                info_map = _parse_kv_pairs(info)
                ann_info = annotated_lookup.get((chrom, pos, ref, alt_value), {})
                gene_name = str(ann_info.get("gene_name") or chrom).strip() or chrom
                effect = str(ann_info.get("effect") or "").strip()
                impact = str(ann_info.get("impact") or "").strip()
                hgvs_c = str(ann_info.get("hgvs_c") or "").strip()
                hgvs_p = str(ann_info.get("hgvs_p") or "").strip()
                qual = _safe_float(qual_raw) or 0.0
                depth_value = _safe_int(info_map.get("DP")) or 0
                ao_value = _safe_float(str(info_map.get("AO") or "").split(",")[0].strip()) or 0.0
                ro_value = _safe_float(str(info_map.get("RO") or "").split(",")[0].strip()) or 0.0
                sample_fields = str(fmt or "").split(":")
                sample_values = str(sample or "").split(":")
                sample_map = {sample_fields[i]: sample_values[i] if i < len(sample_values) else "" for i in range(len(sample_fields))}
                if depth_value <= 0:
                    depth_value = _safe_int(sample_map.get("DP")) or 0
                if ao_value <= 0:
                    ao_value = _safe_float(str(sample_map.get("AO") or "").split(",")[0].strip()) or 0.0
                if ro_value <= 0:
                    ro_value = _safe_float(str(sample_map.get("RO") or "").split(",")[0].strip()) or 0.0
                maf_numeric = 0.0
                try:
                    allele_depth = float(ao_value) + float(ro_value)
                    if allele_depth <= 0 and depth_value > 0:
                        allele_depth = float(depth_value)
                    maf_numeric = round((float(ao_value) / allele_depth), 6) if allele_depth > 0 else 0.0
                except (TypeError, ValueError, ZeroDivisionError):
                    maf_numeric = 0.0
                quality_label = "高质量突变" if (chrom, pos, ref, alt_value) in high_quality_keys else "低质量突变"
                variant_type = str(info_map.get("TYPE") or "").split(",")[0].strip()
                if not variant_type:
                    variant_type = "snp" if len(ref) == len(alt_value) else ("ins" if len(alt_value) > len(ref) else "del")
                rows.append(
                    {
                        "基因名": gene_name,
                        "位置": int(pos) if str(pos).isdigit() else pos,
                        "参考碱基": ref,
                        "突变碱基": alt_value,
                        "测序深度 / 突变频率": " / ".join([part for part in [str(depth_value or ""), f"{maf_numeric:.4f}" if depth_value > 0 else ""] if part]) or "-",
                        "HGVS.c": hgvs_c,
                        "HGVS.p": hgvs_p,
                        "染色体": chrom,
                        "质量分层": quality_label,
                        "质量值QUAL": round(qual, 2),
                        "测序深度DP": depth_value,
                        "突变频率MAF": round(maf_numeric, 6),
                        "影响等级": impact,
                        "注释效应": effect,
                        "变异类型": variant_type,
                        "FILTER": filt,
                    }
                )
    except OSError:
        return {"status": "empty", "columns": [], "rows": [], "total_variants": 0, "source_vcf": ""}
    columns = [
        "基因名",
        "位置",
        "参考碱基",
        "突变碱基",
        "测序深度 / 突变频率",
        "HGVS.c",
        "HGVS.p",
        "染色体",
        "质量分层",
        "影响等级",
        "注释效应",
        "变异类型",
    ]
    list_rows = [[row.get(column, "") for column in columns] for row in rows]
    high_quality = sum(1 for row in rows if str(row.get("质量分层") or "").strip() == "高质量突变")
    return {
        "status": "ready" if rows else "empty",
        "columns": columns,
        "rows": list_rows,
        "total_variants": len(rows),
        "high_quality_variants": high_quality,
        "low_quality_variants": len(rows) - high_quality,
        "source_vcf": str(annotated_vcf_path if annotated_vcf_path.is_file() else raw_vcf_path),
    }

def _normalize_rsv_subtype_label(value: object) -> str:
    text = str(value or "").strip().upper()
    if not text or text == "-":
        return ""
    if text in {"A", "B"}:
        return text
    if "RSV" in text:
        if "A" in text and "B" not in text:
            return "A"
        if "B" in text and "A" not in text:
            return "B"
    return text[-1] if text.endswith(("A", "B")) else ""

def _extract_aa_position_from_hgvs_p(value: object) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    matched = re.search(r"p\.[A-Za-z*]+(\d+)", text)
    if matched:
        return matched.group(1)
    matched = re.search(r"p\.(\d+)", text)
    if matched:
        return matched.group(1)
    return ""

def _normalize_hpiv_subtype_label(value: object) -> str:
    text = str(value or "").strip().upper().replace("HPIV-", "").replace("HPIV", "")
    if not text or text == "-":
        return ""
    if text in {"1", "2", "3", "4A", "4B"}:
        return text
    matched = re.search(r"\b([1-3]|4[A-B])\b", text)
    return matched.group(1) if matched else ""

def _read_hpiv_functional_annotation_table(report_dir: Path, mutation_table: dict, hpiv_subtype: str) -> dict:
    subtype = _normalize_hpiv_subtype_label(hpiv_subtype)
    if not subtype:
        return {"status": "empty", "columns": [], "rows": [], "total_hits": 0}
    db_path = Path(__file__).resolve().parent.parent / "database" / "virus" / "hpiv" / "HPIV_functional_sites_literature_curated.tsv"
    if not db_path.is_file():
        return {"status": "empty", "columns": [], "rows": [], "total_hits": 0}

    mutation_columns = mutation_table.get("columns") if isinstance(mutation_table, dict) else []
    mutation_rows = mutation_table.get("rows") if isinstance(mutation_table, dict) else []
    if not isinstance(mutation_columns, list) or not isinstance(mutation_rows, list) or not mutation_rows:
        return {"status": "empty", "columns": [], "rows": [], "total_hits": 0}

    gene_idx = mutation_columns.index("基因名") if "基因名" in mutation_columns else -1
    pos_idx = mutation_columns.index("位置") if "位置" in mutation_columns else -1
    hgvs_p_idx = mutation_columns.index("HGVS.p") if "HGVS.p" in mutation_columns else -1
    quality_idx = mutation_columns.index("质量分层") if "质量分层" in mutation_columns else -1
    if gene_idx < 0 or hgvs_p_idx < 0:
        return {"status": "empty", "columns": [], "rows": [], "total_hits": 0}

    mutation_lookup: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in mutation_rows:
        if not isinstance(row, list):
            continue
        gene_name = str(row[gene_idx] if gene_idx < len(row) else "").strip().upper()
        aa_position = _extract_aa_position_from_hgvs_p(row[hgvs_p_idx] if hgvs_p_idx < len(row) else "")
        if not gene_name or not aa_position:
            continue
        mutation_lookup.setdefault((gene_name, aa_position), []).append(
            {
                "gene": gene_name,
                "genome_position": str(row[pos_idx] if pos_idx >= 0 and pos_idx < len(row) else "").strip(),
                "hgvs_p": str(row[hgvs_p_idx] if hgvs_p_idx < len(row) else "").strip(),
                "quality": str(row[quality_idx] if quality_idx >= 0 and quality_idx < len(row) else "").strip() or "-",
            }
        )
    if not mutation_lookup:
        return {"status": "empty", "columns": [], "rows": [], "total_hits": 0}

    try:
        with db_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            db_rows = list(reader)
    except OSError:
        return {"status": "empty", "columns": [], "rows": [], "total_hits": 0}

    matched_rows: list[list[str]] = []
    target_virus_type = f"HPIV-{subtype}"
    for record in db_rows:
        virus_type = str(record.get("virus_type") or "").strip().upper()
        if virus_type != target_virus_type:
            continue
        gene_value = str(record.get("gene") or "").strip().upper()
        aa_position = str(record.get("aa_position_hint") or "").strip()
        if not gene_value or not aa_position:
            continue
        current_hits = mutation_lookup.get((gene_value, aa_position))
        if not current_hits:
            continue
        for hit in current_hits:
            matched_rows.append(
                [
                    subtype,
                    hit.get("gene", gene_value),
                    aa_position,
                    hit.get("genome_position", ""),
                    hit.get("hgvs_p", ""),
                    hit.get("quality", "-"),
                    str(record.get("site_or_mutation") or "").strip(),
                    str(record.get("functional_relevance_zh") or record.get("functional_relevance") or "").strip(),
                    str(record.get("observed_phenotype_zh") or record.get("observed_phenotype") or "").strip(),
                    str(record.get("evidence_type_zh") or record.get("evidence_type") or "").strip(),
                    str(record.get("supporting_reference") or "").strip(),
                    str(record.get("doi_or_pmid") or "").strip(),
                    str(record.get("confidence_note_zh") or record.get("confidence_note") or "").strip(),
                ]
            )

    matched_rows.sort(key=lambda row: (0 if row[5] == "高质量突变" else 1, row[1], _safe_int(row[2]) or 0, row[4]))
    columns = [
        "亚型",
        "基因",
        "氨基酸位点",
        "基因组位置",
        "我们的突变结果",
        "质量分层",
        "文献功能位点",
        "功能影响",
        "已知表型",
        "证据类型",
        "参考文献",
        "PMID/DOI",
        "证据说明",
    ]
    return {
        "status": "ready" if matched_rows else "empty",
        "columns": columns,
        "rows": matched_rows,
        "total_hits": len(matched_rows),
    }

def _read_rsv_nmdc_annotation_table(report_dir: Path, mutation_table: dict, rsv_subtype: str) -> dict:
    subtype = _normalize_rsv_subtype_label(rsv_subtype)
    if not subtype:
        return {"status": "empty", "columns": [], "rows": [], "total_hits": 0}
    database_root = _resolve_runtime_database_root()
    db_candidates = [
        database_root / "virus" / "rsv" / "nmdc_hrsv_variation_list_zh.tsv",
        database_root / "virus" / "rsv" / "nmdc_hrsv_variation_list.tsv",
    ]
    db_path = next((path for path in db_candidates if path.is_file()), None)
    if db_path is None:
        return {"status": "empty", "columns": [], "rows": [], "total_hits": 0}

    mutation_columns = mutation_table.get("columns") if isinstance(mutation_table, dict) else []
    mutation_rows = mutation_table.get("rows") if isinstance(mutation_table, dict) else []
    if not isinstance(mutation_columns, list) or not isinstance(mutation_rows, list) or not mutation_rows:
        return {"status": "empty", "columns": [], "rows": [], "total_hits": 0}

    gene_idx = mutation_columns.index("基因名") if "基因名" in mutation_columns else -1
    pos_idx = mutation_columns.index("位置") if "位置" in mutation_columns else -1
    hgvs_p_idx = mutation_columns.index("HGVS.p") if "HGVS.p" in mutation_columns else -1
    quality_idx = mutation_columns.index("质量分层") if "质量分层" in mutation_columns else -1
    if gene_idx < 0:
        return {"status": "empty", "columns": [], "rows": [], "total_hits": 0}

    mutation_lookup: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in mutation_rows:
        if not isinstance(row, list):
            continue
        gene_name = str(row[gene_idx] if gene_idx < len(row) else "").strip()
        aa_position = _extract_aa_position_from_hgvs_p(row[hgvs_p_idx] if hgvs_p_idx < len(row) else "")
        if not gene_name or not aa_position:
            continue
        mutation_lookup.setdefault((gene_name.upper(), aa_position), []).append(
            {
                "gene": gene_name,
                "genome_position": str(row[pos_idx] if pos_idx >= 0 and pos_idx < len(row) else "").strip(),
                "hgvs_p": str(row[hgvs_p_idx] if hgvs_p_idx < len(row) else "").strip(),
                "quality": str(row[quality_idx] if quality_idx >= 0 and quality_idx < len(row) else "").strip() or "-",
            }
        )
    if not mutation_lookup:
        return {"status": "empty", "columns": [], "rows": [], "total_hits": 0}

    try:
        with db_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            db_rows = list(reader)
    except OSError:
        return {"status": "empty", "columns": [], "rows": [], "total_hits": 0}

    def _clean_nmdc_value(value: object) -> str:
        text = str(value or "").strip()
        if not text or text in {"0", "0.0", "None", "none"}:
            return "-"
        if text.startswith("（") or text.startswith("("):
            return "-"
        return text

    matched_rows: list[list[str]] = []
    for record in db_rows:
        reference_value = str(record.get("参考株") or record.get("reference") or "").strip()
        reference_subtype = _normalize_rsv_subtype_label(reference_value)
        if reference_subtype and reference_subtype != subtype:
            continue
        gene_value = str(record.get("基因") or record.get("gene") or "").strip()
        aa_site_value = str(record.get("参考氨基酸位点") or record.get("ref_aminoacid_site") or "").strip()
        if not gene_value or not aa_site_value:
            continue
        current_hits = mutation_lookup.get((gene_value.upper(), aa_site_value))
        if not current_hits:
            continue
        for hit in current_hits:
            matched_rows.append(
                [
                    subtype,
                    hit.get("gene", gene_value),
                    aa_site_value,
                    hit.get("genome_position", ""),
                    hit.get("hgvs_p", ""),
                    hit.get("quality", "-"),
                    str(record.get("氨基酸变异") or record.get("amino_acid_variants") or "").strip(),
                    str(record.get("核酸变异") or record.get("basepair_variants") or "").strip(),
                    _clean_nmdc_value(record.get("抗体亲和") or record.get("anti_risk")),
                    _clean_nmdc_value(record.get("氨基酸突变") or record.get("matrix_risk")),
                    _clean_nmdc_value(record.get("氨基酸变异注释") or record.get("variant_amino_type")),
                    _clean_nmdc_value(record.get("结构相关关键词") or record.get("structure")),
                    _clean_nmdc_value(record.get("功能相关关键词") or record.get("function")),
                ]
            )

    matched_rows.sort(key=lambda row: (0 if row[5] == "高质量突变" else 1, row[1], _safe_int(row[2]) or 0, row[4]))
    columns = [
        "亚型",
        "基因",
        "突变位点",
        "基因组位置",
        "我们的突变结果",
        "质量分层",
        "数据库氨基酸变异",
        "数据库核酸变异",
        "抗体亲和",
        "氨基酸突变",
        "变异注释",
        "结构相关关键词",
        "功能相关关键词",
    ]
    return {
        "status": "ready" if matched_rows else "empty",
        "columns": columns,
        "rows": matched_rows,
        "total_hits": len(matched_rows),
        "source": str(db_path),
    }

def _read_influenza_resistance_annotation_table(report_dir: Path) -> dict:
    json_path = report_dir / "snps.filt1.resistance_annotation.json"
    tsv_path = report_dir / "snps.filt1.resistance_annotation.tsv"
    if json_path.is_file():
        try:
            payload = json.loads(json_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            payload = {}
        rows = payload.get("rows") if isinstance(payload.get("rows"), list) else []
        if rows and isinstance(rows[0], dict):
            columns = [
                "基因",
                "规则突变",
                "命中突变",
                "药物类别",
                "药物",
                "风险等级",
                "权威结论",
                "注释说明",
                "证据来源",
                "证据强度",
                "适用范围",
                "备注",
            ]
            return {
                "status": str(payload.get("status") or "ready"),
                "columns": columns,
                "rows": [[row.get(column, "") for column in columns] for row in rows],
                "total_hits": int(payload.get("total_hits") or len(rows)),
            }
    raw = _read_tsv_rows(tsv_path)
    return {
        "status": "ready" if raw.get("rows") is not None and raw.get("columns") else "empty",
        "columns": raw.get("columns", []),
        "rows": raw.get("rows", []),
        "total_hits": len(raw.get("rows", [])),
    }
