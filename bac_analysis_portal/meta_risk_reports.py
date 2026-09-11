from __future__ import annotations

from pathlib import Path

from .parse_utils import _safe_int
from .table_io import _read_tsv_rows

def _read_gene_length_distribution(path: Path) -> dict:
    raw = _read_tsv_rows(path)
    if not raw["columns"] or not raw["rows"]:
        return {"status": "empty", "points": []}
    x_values: list[str] = []
    points: list[int] = []
    for row in raw["rows"]:
        label = row[0] if row else ""
        value = _safe_int(row[1] if len(row) > 1 else None)
        if not label or value is None:
            continue
        x_values.append(str(label))
        points.append(value)
    if not points:
        return {"status": "empty", "points": []}
    return {
        "status": "ready",
        "label": "基因长度与数量分布",
        "x_label": raw["columns"][0] if raw["columns"] else "基因长度范围",
        "y_label": raw["columns"][1] if len(raw["columns"]) > 1 else "Gene数量",
        "x_values": x_values,
        "points": points,
        "max_count": max(points),
    }

def _build_resistance_virulence_summary(rv_summary: dict, virulence_elements: dict, resistance_elements: dict) -> dict:
    rv_rows = rv_summary.get("rows") or []
    rv_columns = rv_summary.get("columns") or []
    virulence_rows = virulence_elements.get("rows") or []
    resistance_rows = resistance_elements.get("rows") or []

    def _count_by_column(rows: list[list[str]], columns: list[str], target: str, split_delimiters: list[str] | None = None) -> list[dict]:
        if target not in columns:
            return []
        index = columns.index(target)
        counts: dict[str, int] = {}
        for row in rows:
            raw = str(row[index] if index < len(row) else "").strip()
            if not raw or raw == '-':
                continue
            parts = [raw]
            if split_delimiters:
                parts = [raw]
                for delimiter in split_delimiters:
                    expanded = []
                    for item in parts:
                        expanded.extend(item.split(delimiter))
                    parts = expanded
            for item in [part.strip() for part in parts if part.strip() and part.strip() != '-']:
                counts[item] = counts.get(item, 0) + 1
        return [
            {"label": label, "count": count}
            for label, count in sorted(counts.items(), key=lambda item: item[1], reverse=True)[:8]
        ]

    def _rv_metric(column_name: str) -> int:
        if column_name not in rv_columns:
            return 0
        idx = rv_columns.index(column_name)
        total = 0
        for row in rv_rows:
            try:
                total += int(float(str(row[idx]).strip() or '0'))
            except ValueError:
                continue
        return total

    virulence_top = _count_by_column(virulence_rows, virulence_elements.get("columns") or [], "VF分类")
    resistance_top = _count_by_column(resistance_rows, resistance_elements.get("columns") or [], "耐药药物", split_delimiters=[';'])
    virulence_genes = len({str(row[(virulence_elements.get("columns") or []).index("基因名称")]).strip() for row in virulence_rows if "基因名称" in (virulence_elements.get("columns") or []) and str(row[(virulence_elements.get("columns") or []).index("基因名称")]).strip() not in {'', '-'}}) if virulence_rows else 0
    resistance_genes = len({str(row[(resistance_elements.get("columns") or []).index("基因名称")]).strip() for row in resistance_rows if "基因名称" in (resistance_elements.get("columns") or []) and str(row[(resistance_elements.get("columns") or []).index("基因名称")]).strip() not in {'', '-'}}) if resistance_rows else 0

    resistance_note = f"当前共检出 {len(resistance_rows)} 条耐药元件记录、{resistance_genes} 个耐药基因，主要集中于“{resistance_top[0]['label']}”相关类别。" if resistance_top else f"当前共检出 {len(resistance_rows)} 条耐药元件记录、{resistance_genes} 个耐药基因。"
    virulence_note = f"当前共检出 {len(virulence_rows)} 条毒力元件记录、{virulence_genes} 个毒力基因，主要集中于“{virulence_top[0]['label']}”类别。" if virulence_top else f"当前共检出 {len(virulence_rows)} 条毒力元件记录、{virulence_genes} 个毒力基因。"

    return {
        "status": "ready" if rv_rows or virulence_rows or resistance_rows else "empty",
        "resistance": {
            "note": resistance_note,
            "hit_count": len(resistance_rows),
            "gene_count": resistance_genes,
            "top_categories": resistance_top,
            "summary_count": _rv_metric("耐药基因数量"),
        },
        "virulence": {
            "note": virulence_note,
            "hit_count": len(virulence_rows),
            "gene_count": virulence_genes,
            "top_categories": virulence_top,
            "summary_count": _rv_metric("毒力基因数量"),
        },
    }

def _meta_sequence_label(value: str) -> str:
    primary = str(value or "").strip().split(".", 1)[0].strip().lower()
    mapping = {
        "plasmid": "质粒",
        "chromosome": "染色体",
        "unclassified": "未分类",
    }
    return mapping.get(primary, primary or "未分类")

def _meta_taxonomy_label(value: str) -> str:
    parts = str(value or "").strip().split(".", 1)
    if len(parts) < 2:
        return "未分类"
    tax = parts[1].strip()
    return tax if tax and tax != "-" else "未分类"

def _normalize_meta_gene(value: str) -> str:
    return str(value or "").strip().lower()

def _extract_vf_category_from_product(product: str) -> str:
    text = str(product or "").strip()
    if not text:
        return "未分类"
    if "[" in text and "]" in text:
        inner = text[text.find("[") + 1:text.rfind("]")]
        if " - " in inner:
            tail = inner.split(" - ", 1)[1]
            if " (" in tail:
                tail = tail.split(" (", 1)[0]
            tail = tail.strip()
            if tail:
                return tail
    return "未分类"

def _load_meta_vf_categories(base_dir: Path) -> dict[tuple[str, str], str]:
    raw = _read_tsv_rows(base_dir / "bin_vfdb.tsv")
    columns = raw.get("columns") or []
    rows = raw.get("rows") or []
    if not columns or not rows:
        return {}
    index = {name: idx for idx, name in enumerate(columns)}
    result: dict[tuple[str, str], str] = {}
    for row in rows:
        try:
            file_name = str(row[index["#FILE"]]).strip().split("/")[-1].split(".")[0]
            sequence = str(row[index["SEQUENCE"]]).strip()
            gene = _normalize_meta_gene(row[index["GENE"]])
            product = str(row[index["PRODUCT"]]).strip()
        except Exception:
            continue
        if not file_name or not sequence or not gene:
            continue
        result[(f"{file_name}_{sequence}", gene)] = _extract_vf_category_from_product(product)
    return result

def _load_meta_card_categories(base_dir: Path) -> dict[tuple[str, str], str]:
    raw = _read_tsv_rows(base_dir / "bin_card.tsv")
    columns = raw.get("columns") or []
    rows = raw.get("rows") or []
    if not columns or not rows:
        return {}
    index = {name: idx for idx, name in enumerate(columns)}
    result: dict[tuple[str, str], str] = {}
    for row in rows:
        try:
            file_name = str(row[index["#FILE"]]).strip().split("/")[-1].split(".")[0]
            sequence = str(row[index["SEQUENCE"]]).strip()
            gene = _normalize_meta_gene(row[index["GENE"]])
            resistance = str(row[index["RESISTANCE"]]).strip() or str(row[index.get("PRODUCT", -1)]).strip()
        except Exception:
            continue
        if not file_name or not sequence or not gene:
            continue
        result[(f"{file_name}_{sequence}", gene)] = resistance or "未分类"
    return result

def _load_meta_rgi_categories(base_dir: Path) -> dict[tuple[str, str], str]:
    raw = _read_tsv_rows(base_dir / "binning_rgi_new.txt")
    columns = raw.get("columns") or []
    rows = raw.get("rows") or []
    if not columns or not rows:
        return {}
    index = {name: idx for idx, name in enumerate(columns)}
    result: dict[tuple[str, str], str] = {}
    for row in rows:
        try:
            cutoff = str(row[index["Cut_Off"]]).strip()
            contig = str(row[index["Contig"]]).strip()
            gene = _normalize_meta_gene(row[index["Best_Hit_ARO"]])
            drug_class = str(row[index.get("Drug Class", -1)]).strip()
            antibiotic = str(row[index.get("Antibiotic", -1)]).strip()
        except Exception:
            continue
        if cutoff not in {"Strict", "Perfect"} or not contig or not gene:
            continue
        result[(contig, gene)] = drug_class or antibiotic or "未分类"
    return result

def _load_meta_resfinder_categories(base_dir: Path) -> dict[tuple[str, str], str]:
    raw = _read_tsv_rows(base_dir / "staramr_result" / "resfinder.tsv")
    columns = raw.get("columns") or []
    rows = raw.get("rows") or []
    if not columns or not rows:
        return {}
    index = {name: idx for idx, name in enumerate(columns)}
    result: dict[tuple[str, str], str] = {}
    for row in rows:
        try:
            isolate_id = str(row[index["Isolate ID"]]).strip()
            contig = str(row[index["Contig"]]).strip()
            gene = _normalize_meta_gene(row[index["Gene"]])
            phenotype = str(row[index.get("CGE Predicted Phenotype", -1)]).strip() or str(row[index.get("Predicted Phenotype", -1)]).strip()
        except Exception:
            continue
        if not isolate_id or not contig or not gene:
            continue
        result[(f"{isolate_id}_{contig}", gene)] = phenotype or "未分类"
    return result

def _read_meta_resistance_virulence(path: Path) -> dict:
    raw = _read_tsv_rows(path)
    columns = raw.get("columns") or []
    rows = raw.get("rows") or []
    if not rows or not columns:
        return {
            "summary": {"columns": [], "rows": []},
            "virulence_elements": {"columns": [], "rows": []},
            "resistance_elements": {"columns": [], "rows": []},
        }

    base_dir = path.parent
    vf_categories = _load_meta_vf_categories(base_dir)
    card_categories = _load_meta_card_categories(base_dir)
    rgi_categories = _load_meta_rgi_categories(base_dir)
    resfinder_categories = _load_meta_resfinder_categories(base_dir)

    sample_depth_column = next((column for column in columns if column not in {
        "contig_name", "label", "Plasmid", "VF Gene", "ARG", "Name",
        "AR Gene(abricate)", "AR Gene(rgi)", "AR Gene(resfinder)",
        "D", "P", "C", "O", "F", "G", "S",
    }), "")
    records = [{columns[index]: row[index] if index < len(row) else "" for index in range(len(columns))} for row in rows]

    summary_columns = ["Bin名称", "物种名称", "毒力基因数量", "耐药基因数量", "质粒相关Contig数", "染色体相关Contig数", "平均深度"]
    summary_groups: dict[tuple[str, str], dict[str, object]] = {}
    virulence_columns = ["Contig名称", "Bin名称", "物种名称", "VF分类", "基因名称", "基因组/质粒", "门", "属", "种", "平均深度"]
    resistance_columns = ["Contig名称", "Bin名称", "物种名称", "耐药药物", "基因名称", "基因组/质粒", "软件支持", "门", "属", "种", "平均深度"]
    virulence_rows: list[list[str]] = []
    resistance_rows: list[list[str]] = []

    for record in records:
        contig_name = str(record.get("contig_name", "")).strip()
        bin_name = str(record.get("Name", "")).strip() or contig_name.split("_", 1)[0]
        species_name = str(record.get("S", "")).strip() or "-"
        seq_label = _meta_sequence_label(record.get("label", ""))
        depth_value = str(record.get(sample_depth_column, "")).strip() if sample_depth_column else ""
        genus_name = str(record.get("G", "")).strip() or "-"
        phylum_name = str(record.get("P", "")).strip() or "-"
        support_flags = [
            label for label, column_name in (
                ("abricate", "AR Gene(abricate)"),
                ("rgi", "AR Gene(rgi)"),
                ("resfinder", "AR Gene(resfinder)"),
            ) if str(record.get(column_name, "")).strip() == "+"
        ]
        group_key = (bin_name, species_name)
        group = summary_groups.setdefault(group_key, {
            "vf_genes": set(),
            "arg_genes": set(),
            "plasmid_contigs": set(),
            "chromosome_contigs": set(),
            "depth_values": [],
        })
        if depth_value not in {"", "-", "NA"}:
            try:
                group["depth_values"].append(float(depth_value))
            except ValueError:
                pass
        if seq_label == "质粒" and contig_name:
            group["plasmid_contigs"].add(contig_name)
        if seq_label == "染色体" and contig_name:
            group["chromosome_contigs"].add(contig_name)

        vf_gene = str(record.get("VF Gene", "")).strip()
        if vf_gene and vf_gene != "-":
            group["vf_genes"].add(vf_gene)
            vf_category = vf_categories.get((contig_name, _normalize_meta_gene(vf_gene)), "未分类")
            virulence_rows.append([
                contig_name,
                bin_name,
                species_name,
                vf_category,
                vf_gene,
                seq_label,
                phylum_name,
                genus_name,
                species_name,
                depth_value or "-",
            ])

        arg_gene = str(record.get("ARG", "")).strip()
        if arg_gene and arg_gene != "-":
            group["arg_genes"].add(arg_gene)
            normalized_arg = _normalize_meta_gene(arg_gene)
            resistance_labels: list[str] = []
            if str(record.get("AR Gene(abricate)", "")).strip() == "+":
                resistance_labels.append(card_categories.get((contig_name, normalized_arg), "未分类"))
            if str(record.get("AR Gene(rgi)", "")).strip() == "+":
                resistance_labels.append(rgi_categories.get((contig_name, normalized_arg), "未分类"))
            if str(record.get("AR Gene(resfinder)", "")).strip() == "+":
                resistance_labels.append(resfinder_categories.get((contig_name, normalized_arg), "未分类"))
            resistance_labels = [label for label in resistance_labels if label and label != "-"]
            resistance_rows.append([
                contig_name,
                bin_name,
                species_name,
                "; ".join(dict.fromkeys(resistance_labels)) if resistance_labels else "未分类",
                arg_gene,
                seq_label,
                " / ".join(support_flags) if support_flags else "-",
                phylum_name,
                genus_name,
                species_name,
                depth_value or "-",
            ])

    summary_rows: list[list[str]] = []
    for (bin_name, species_name), stats in sorted(summary_groups.items(), key=lambda item: item[0][0]):
        depth_values = stats["depth_values"]
        avg_depth = f"{(sum(depth_values) / len(depth_values)):.2f}" if depth_values else "-"
        summary_rows.append([
            bin_name,
            species_name,
            str(len(stats["vf_genes"])),
            str(len(stats["arg_genes"])),
            str(len(stats["plasmid_contigs"])),
            str(len(stats["chromosome_contigs"])),
            avg_depth,
        ])

    return {
        "summary": {"columns": summary_columns, "rows": summary_rows},
        "virulence_elements": {"columns": virulence_columns, "rows": virulence_rows},
        "resistance_elements": {"columns": resistance_columns, "rows": resistance_rows},
    }

def _risk_level_label(level: str) -> str:
    mapping = {
        "A": "Level A 高迁移风险",
        "B": "Level B 中迁移风险",
        "C": "Level C 低到中等风险",
        "D": "Level D 弱证据",
    }
    return mapping.get(str(level).strip().upper(), str(level).strip() or "未分级")

def _risk_statement_label(level: str) -> str:
    mapping = {
        "A": "高度提示可移动或潜在可转移",
        "B": "可能可被动员或与移动元件相关",
        "C": "属于移动元件相关邻域，不建议直接推断可转移",
        "D": "证据较弱，通常不应据此推断可移动",
    }
    return mapping.get(str(level).strip().upper(), "")

def _mge_type_label(value: str) -> str:
    text = str(value or "").strip()
    primary = text.split("|")[0].strip() if "|" in text else text
    normalized = primary.lower()
    mapping = {
        "plasmid": "质粒",
        "provirus": "前噬菌体",
        "virus": "病毒/噬菌体",
        "phage": "噬菌体",
        "ie": "整合/切除模块",
        "transfer": "转移模块",
        "boundary": "预测边界特征",
        "rrr_std": "RRR/STD（复制/重组/修复；稳定/转移/防御）",
        "stability/transfer/defense": "STD（稳定/转移/防御）",
        "replication/recombination/repair": "RRR（复制/重组/修复）",
        "mge": "移动元件相关",
    }
    if normalized in mapping:
        return mapping[normalized]
    if not normalized:
        return "未识别"
    return primary

def _read_mge_monitoring(report_dir: Path, sample_name: str) -> dict:
    if not sample_name:
        return {"status": "empty"}
    risk_raw = _read_tsv_rows(report_dir / f"{sample_name}.mge_risk_summary.tsv")
    integrated_raw = _read_tsv_rows(report_dir / f"{sample_name}.integrated_mge_summary.tsv")
    if not risk_raw.get("rows") and not integrated_raw.get("rows"):
        return {"status": "empty"}

    integrated_columns = integrated_raw.get("columns") or []
    integrated_records = [
        {integrated_columns[index]: row[index] if index < len(row) else "" for index in range(len(integrated_columns))}
        for row in integrated_raw.get("rows") or []
    ]
    integrated_table_columns = ["样本名称", "元件类型", "所在序列", "起始", "终止", "长度(bp)", "得分", "识别方法", "注释"]
    integrated_rows: list[list[str]] = []
    mge_type_by_sequence: dict[str, set[str]] = {}
    overall_mge_counts: dict[str, int] = {}
    for record in integrated_records:
        sequence_id = str(record.get("sequence_id", "")).strip()
        mge_type = _mge_type_label(str(record.get("mge_type", "")).strip())
        if sequence_id:
            integrated_rows.append([
                sample_name,
                mge_type,
                sequence_id,
                str(record.get("start", "")).strip(),
                str(record.get("end", "")).strip(),
                str(record.get("length", "")).strip(),
                str(record.get("score", "")).strip(),
                str(record.get("method", "")).strip(),
                str(record.get("annotation", "")).strip(),
            ])
        if not sequence_id:
            continue
        mge_type_by_sequence.setdefault(sequence_id, set()).add(mge_type)
        overall_mge_counts[mge_type] = overall_mge_counts.get(mge_type, 0) + 1

    risk_columns = risk_raw.get("columns") or []
    risk_records = [
        {risk_columns[index]: row[index] if index < len(row) else "" for index in range(len(risk_columns))}
        for row in risk_raw.get("rows") or []
    ]
    table_columns = [
        "样本名称",
        "基因类型",
        "基因名称",
        "产物/功能",
        "所在序列",
        "基因起始",
        "基因终止",
        "关联元件类型",
        "转移风险等级",
        "风险说明",
        "位于预测 MGE 边界",
        "最近核心模块",
        "最近核心模块类型",
        "最近距离(bp)",
        "最近距离(ORF)",
        "同序列核心类别",
        "存在第二核心模块或边界",
    ]
    resistance_rows: list[list[str]] = []
    virulence_rows: list[list[str]] = []
    resistance_identified_hits = 0
    virulence_identified_hits = 0
    risk_counts_all: dict[str, int] = {}
    risk_counts_by_type: dict[str, dict[str, int]] = {"ARG": {}, "VF": {}}
    mge_counts_by_type: dict[str, dict[str, int]] = {"ARG": {}, "VF": {}}

    for record in risk_records:
        gene_type = str(record.get("gene_type", "")).strip().upper()
        if gene_type not in {"ARG", "VF"}:
            continue
        level = str(record.get("risk_level", "")).strip().upper()
        risk_counts_all[level] = risk_counts_all.get(level, 0) + 1
        risk_counts_by_type[gene_type][level] = risk_counts_by_type[gene_type].get(level, 0) + 1
        sequence_id = str(record.get("sequence_id", "")).strip()
        related_mge_types = sorted(mge_type_by_sequence.get(sequence_id) or set())
        if not related_mge_types:
            related_mge_types = [
                _mge_type_label(str(record.get("nearest_core_type", "")).strip())
            ] if str(record.get("nearest_core_type", "")).strip() else ["未识别"]
        related_mge_types = [label for label in related_mge_types if label and label != "未识别"]
        if not related_mge_types:
            related_mge_types = ["未识别"]
        for label in related_mge_types:
            mge_counts_by_type[gene_type][label] = mge_counts_by_type[gene_type].get(label, 0) + 1
        row = [
            sample_name,
            "耐药基因" if gene_type == "ARG" else "毒力基因",
            str(record.get("gene_name", "")).strip(),
            str(record.get("product", "")).strip(),
            sequence_id,
            str(record.get("gene_start", "")).strip(),
            str(record.get("gene_end", "")).strip(),
            "、".join(related_mge_types),
            _risk_level_label(level),
            _risk_statement_label(level),
            "是" if str(record.get("within_mge_boundary", "")).strip().lower() == "yes" else "否",
            str(record.get("nearest_core_gene", "")).strip() or "-",
            _mge_type_label(str(record.get("nearest_core_type", "")).strip()),
            str(record.get("nearest_core_distance_bp", "")).strip() or "-",
            str(record.get("nearest_core_distance_orf", "")).strip() or "-",
            str(record.get("same_contig_core_categories", "")).strip() or "-",
            "是" if str(record.get("has_second_core_or_boundary", "")).strip().lower() == "yes" else "否",
        ]
        if gene_type == "ARG":
            resistance_rows.append(row)
            if any(label != "未识别" for label in related_mge_types):
                resistance_identified_hits += 1
        else:
            virulence_rows.append(row)
            if any(label != "未识别" for label in related_mge_types):
                virulence_identified_hits += 1

    def _sorted_counts(mapping: dict[str, int], order: list[str] | None = None) -> list[dict]:
        if order:
            items = [(key, mapping.get(key, 0)) for key in order if mapping.get(key, 0)]
        else:
            items = sorted(mapping.items(), key=lambda item: (-item[1], item[0]))
        return [{"label": key, "count": value} for key, value in items]

    def _build_block(title: str, rows: list[list[str]], key: str) -> dict:
        risk_counts = _sorted_counts(
            {_risk_level_label(level): count for level, count in risk_counts_by_type[key].items()},
            ["Level A 高迁移风险", "Level B 中迁移风险", "Level C 低到中等风险", "Level D 弱证据"],
        )
        mge_counts = _sorted_counts(mge_counts_by_type[key])
        high_risk = sum(risk_counts_by_type[key].get(level, 0) for level in ("A", "B"))
        relation_rows = []
        for row in rows:
            mge_types = str(row[7] if len(row) > 7 else "").strip()
            gene_name = str(row[2] if len(row) > 2 else "").strip()
            if not mge_types or not gene_name:
                continue
            filtered_types = "、".join([item for item in mge_types.split("、") if item.strip() and item.strip() != "未识别"])
            if not filtered_types:
                continue
            relation_rows.append([filtered_types, gene_name])
        note = (
            f"{title}共检出 {len(rows)} 条与移动元件相关记录，其中 {high_risk} 条处于 Level A/B，"
            f"提示存在较高或中等的潜在转移风险。"
            if rows else f"{title}未检出可展示的移动元件监测结果。"
        )
        return {
            "note": note,
            "hit_count": len(rows),
            "high_risk_count": high_risk,
            "risk_levels": risk_counts,
            "mge_types": mge_counts,
            "gene_mge_relationship": _build_category_gene_relationship(
                {"columns": ["移动元件类型", "基因名称"], "rows": relation_rows},
                left_key="移动元件类型",
                right_key="基因名称",
                label=f"{title}基因与移动元件类型关系图",
                split_delimiters=["、"],
            ),
            "columns": table_columns,
            "rows": rows,
        }

    overview_risk = _sorted_counts(
        {_risk_level_label(level): count for level, count in risk_counts_all.items()},
        ["Level A 高迁移风险", "Level B 中迁移风险", "Level C 低到中等风险", "Level D 弱证据"],
    )
    overview_mge = _sorted_counts(overall_mge_counts)
    return {
        "status": "ready" if integrated_rows or resistance_rows or virulence_rows else "empty",
        "overview": {
            "note": (
                f"当前共检出 {len(integrated_records)} 条移动元件记录，"
                f"其中与耐药/毒力基因相关的记录共 {len(risk_records)} 条，"
                f"Level A/B 共 {sum(risk_counts_all.get(level, 0) for level in ('A', 'B'))} 条。"
            ) if integrated_records or risk_records else "当前未检出可展示的移动元件监测记录。",
            "total_hits": len(integrated_records),
            "resistance_hits": resistance_identified_hits,
            "virulence_hits": virulence_identified_hits,
            "risk_levels": overview_risk,
            "mge_types": overview_mge,
        },
        "elements": {
            "columns": integrated_table_columns,
            "rows": integrated_rows,
        },
        "resistance": _build_block("耐药相关位点", resistance_rows, "ARG"),
        "virulence": _build_block("毒力相关位点", virulence_rows, "VF"),
    }

def _merge_annotation_with_summary(detail_path: Path, summary_path: Path, mode: str) -> dict:
    detail = _read_tsv_rows(detail_path)
    summary = _read_tsv_rows(summary_path)
    if not detail["columns"] and not summary["columns"]:
        return {"columns": [], "rows": []}

    detail_records = [
        {detail["columns"][index]: row[index] if index < len(row) else "" for index in range(len(detail["columns"]))}
        for row in detail["rows"]
    ]
    summary_records = [
        {summary["columns"][index]: row[index] if index < len(row) else "" for index in range(len(summary["columns"]))}
        for row in summary["rows"]
    ]
    summary_index: dict[tuple[str, str], dict] = {}
    for record in summary_records:
        key = (
            str(record.get("基因名称", "")).strip(),
            str(record.get("片段名称", "")).strip(),
        )
        summary_index[key] = record

    if mode == "virulence":
        columns = [
            "Contig名称", "物种名称", "taxid", "基因名称", "覆盖度%", "一致性%", "产物",
            "VF分类", "VF名称", "起始碱基", "终止碱基", "正负链",
            "覆盖度(>0)%", "覆盖度(>10)%", "覆盖度(>100)%", "平均深度", "最低深度", "最高深度",
        ]
    else:
        columns = [
            "Contig名称", "物种名称", "taxid", "基因名称", "覆盖度%", "一致性%", "产物",
            "耐药药物", "起始碱基", "终止碱基", "正负链",
            "覆盖度(>0)%", "覆盖度(>10)%", "覆盖度(>100)%", "平均深度", "最低深度", "最高深度",
        ]

    rows = []
    for record in detail_records:
        summary_record = summary_index.get(
            (
                str(record.get("基因名称", "")).strip(),
                str(record.get("Contig名称", "")).strip(),
            ),
            {},
        )
        merged = []
        for column in columns:
            if column in record:
                merged.append(record.get(column, ""))
            elif column in summary_record:
                merged.append(summary_record.get(column, ""))
            else:
                merged.append("")
        rows.append(merged)

    if rows:
        return {"columns": columns, "rows": rows}

    fallback_columns = summary["columns"] or detail["columns"]
    fallback_rows = summary["rows"] or detail["rows"]
    return {"columns": fallback_columns, "rows": fallback_rows}

def _build_category_gene_relationship(data: dict, left_key: str, right_key: str, label: str, split_delimiters: list[str] | None = None) -> dict:
    columns = data.get("columns") or []
    rows = data.get("rows") or []
    if not columns or not rows:
        return {"status": "empty", "nodes_left": [], "nodes_right": [], "links": []}
    left_index = columns.index(left_key) if left_key in columns else -1
    right_index = columns.index(right_key) if right_key in columns else -1
    if left_index < 0 or right_index < 0:
        return {"status": "empty", "nodes_left": [], "nodes_right": [], "links": []}
    split_delimiters = split_delimiters or []
    link_counts: dict[tuple[str, str], int] = {}
    left_totals: dict[str, int] = {}
    right_totals: dict[str, int] = {}
    for row in rows:
        left_value = str(row[left_index] if left_index < len(row) else "").strip()
        right_value = str(row[right_index] if right_index < len(row) else "").strip()
        if not left_value or not right_value or left_value == "-" or right_value == "-":
            continue
        left_values = [left_value]
        for delimiter in split_delimiters:
            if delimiter in left_value:
                left_values = [item.strip() for item in left_value.split(delimiter) if item.strip()]
                break
        for left_item in left_values:
            key = (left_item, right_value)
            link_counts[key] = link_counts.get(key, 0) + 1
            left_totals[left_item] = left_totals.get(left_item, 0) + 1
            right_totals[right_value] = right_totals.get(right_value, 0) + 1
    if not link_counts:
        return {"status": "empty", "nodes_left": [], "nodes_right": [], "links": []}
    top_left = sorted(left_totals.items(), key=lambda item: (-item[1], item[0]))[:8]
    top_right = sorted(right_totals.items(), key=lambda item: (-item[1], item[0]))[:12]
    keep_left = {name for name, _ in top_left}
    keep_right = {name for name, _ in top_right}
    links = [
        {"source": left, "target": right, "value": count}
        for (left, right), count in sorted(link_counts.items(), key=lambda item: (-item[1], item[0][0], item[0][1]))
        if left in keep_left and right in keep_right
    ]
    return {
        "status": "ready" if links else "empty",
        "label": label,
        "left_label": left_key,
        "right_label": right_key,
        "nodes_left": [{"name": name, "value": value} for name, value in top_left],
        "nodes_right": [{"name": name, "value": value} for name, value in top_right],
        "links": links,
    }
