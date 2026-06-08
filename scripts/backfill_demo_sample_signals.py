#!/usr/bin/env python3
from __future__ import annotations

import json
import sqlite3
from pathlib import Path


ROOT = Path("/Users/wuhhh/Desktop/徐老师/代码/metagenomic")
DB_PATH = ROOT / "bac_analysis_portal.sqlite3"


def extract_mlst_st(payload: object) -> str:
    sections = payload.get("sections", {}) if isinstance(payload, dict) else {}
    mlst_section = sections.get("mlst", {}) if isinstance(sections.get("mlst"), dict) else {}
    columns = mlst_section.get("columns") if isinstance(mlst_section.get("columns"), list) else []
    rows = mlst_section.get("rows") if isinstance(mlst_section.get("rows"), list) else []
    if not columns or not rows:
        return ""
    target_index = next(
        (idx for idx, label in enumerate(columns) if str(label or "").strip() in {"序列分型(ST)", "序列分型", "ST"}),
        -1,
    )
    if target_index < 0:
        return ""
    seen: list[str] = []
    for row in rows:
        if not isinstance(row, list) or target_index >= len(row):
            continue
        value = str(row[target_index] or "").strip()
        if not value or value == "-":
            continue
        normalized = value if value.upper().startswith("ST") else f"ST{value}"
        if normalized not in seen:
            seen.append(normalized)
    return " / ".join(seen[:3])


def extract_serotype(payload: object) -> str:
    sections = payload.get("sections", {}) if isinstance(payload, dict) else {}
    tables: list[dict] = []
    if isinstance(sections.get("priority_serotype"), dict):
        tables.append(sections["priority_serotype"])
    if isinstance(sections.get("serotype"), dict):
        tables.append(sections["serotype"])
    candidates = {"知识库命中血清型", "血清型", "血清群", "亚型预测", "predicted_serotype", "Serotype", "New_serotype"}
    for table in tables:
        columns = table.get("columns") if isinstance(table.get("columns"), list) else []
        rows = table.get("rows") if isinstance(table.get("rows"), list) else []
        if not columns or not rows:
            continue
        for idx, label in enumerate(columns):
            if str(label or "").strip() not in candidates:
                continue
            for row in rows:
                if not isinstance(row, list) or idx >= len(row):
                    continue
                value = str(row[idx] or "").strip()
                if value and value != "-":
                    return value
    return ""


def find_column_index(columns: list[object], candidates: tuple[str, ...]) -> int:
    labels = [str(label or "").strip() for label in (columns or [])]
    for candidate in candidates:
        if candidate in labels:
            return labels.index(candidate)
    return -1


def summarize_detected_genes(table: object, *, limit: int = 8) -> str:
    if not isinstance(table, dict):
        return ""
    columns = table.get("columns") if isinstance(table.get("columns"), list) else []
    rows = table.get("rows") if isinstance(table.get("rows"), list) else []
    gene_index = find_column_index(columns, ("基因名称", "gene_name", "Gene"))
    if gene_index < 0 or not rows:
        return ""
    seen: list[str] = []
    for row in rows:
        if not isinstance(row, list) or gene_index >= len(row):
            continue
        gene_name = str(row[gene_index] or "").strip()
        if not gene_name or gene_name == "-" or gene_name in seen:
            continue
        seen.append(gene_name)
    if not seen:
        return ""
    if len(seen) <= limit:
        return "、".join(seen)
    return f"{'、'.join(seen[:limit])} 等 {len(seen)} 项"


def summarize_mge_gene_hits(table: object, *, limit: int = 6) -> str:
    if not isinstance(table, dict):
        return ""
    columns = table.get("columns") if isinstance(table.get("columns"), list) else []
    rows = table.get("rows") if isinstance(table.get("rows"), list) else []
    gene_index = find_column_index(columns, ("基因名称", "gene_name", "Gene"))
    mge_index = find_column_index(columns, ("关联元件类型", "元件类型", "MGE类型"))
    if gene_index < 0 or mge_index < 0 or not rows:
        return ""
    seen: list[str] = []
    for row in rows:
        if not isinstance(row, list) or gene_index >= len(row) or mge_index >= len(row):
            continue
        gene_name = str(row[gene_index] or "").strip()
        mge_type = str(row[mge_index] or "").strip()
        if not gene_name or gene_name == "-":
            continue
        normalized_types = [item.strip() for item in mge_type.split("、") if item.strip() and item.strip() != "未识别"]
        if not normalized_types:
            continue
        summary = f"{gene_name}（{'/'.join(normalized_types)}）"
        if summary in seen:
            continue
        seen.append(summary)
    if not seen:
        return ""
    if len(seen) <= limit:
        return "、".join(seen)
    return f"{'、'.join(seen[:limit])} 等 {len(seen)} 项"


def build_signal_summary(report_dir: Path, sample_name: str) -> dict[str, str] | None:
    payload_path = report_dir / ".portal_report_cache" / f"report_payload_{sample_name}.json"
    if not payload_path.is_file():
        return None
    payload = json.loads(payload_path.read_text()).get("data", {})
    sections = payload.get("sections", {}) if isinstance(payload, dict) else {}
    rv_sections = sections.get("resistance_virulence", {}) if isinstance(sections.get("resistance_virulence"), dict) else {}
    mge_sections = sections.get("mge_monitoring", {}) if isinstance(sections.get("mge_monitoring"), dict) else {}
    return {
        "mlst_st": extract_mlst_st(payload),
        "serotype_result": extract_serotype(payload),
        "resistance_gene_hits": summarize_detected_genes(rv_sections.get("resistance_elements")),
        "virulence_gene_hits": summarize_detected_genes(rv_sections.get("virulence_elements")),
        "resistance_mge_hits": summarize_mge_gene_hits(mge_sections.get("resistance")),
        "virulence_mge_hits": summarize_mge_gene_hits(mge_sections.get("virulence")),
    }


def main() -> None:
    conn = sqlite3.connect(DB_PATH)
    cur = conn.cursor()
    rows = cur.execute(
        """
        select sample_key, sample_name, report_dir,
               mlst_st, serotype_result,
               resistance_gene_hits, virulence_gene_hits,
               resistance_mge_hits, virulence_mge_hits
          from sample_library
         where report_dir like ?
         order by sample_name, sample_key
        """,
        ("%demo_data%",),
    ).fetchall()
    updated = 0
    skipped = 0
    for row in rows:
        sample_key, sample_name, report_dir, *_ = row
        summary = build_signal_summary(Path(str(report_dir)).resolve(), str(sample_name))
        if not summary:
            skipped += 1
            print(f"SKIP\t{sample_name}\tmissing_payload")
            continue
        cur.execute(
            """
            update sample_library
               set mlst_st = ?,
                   serotype_result = ?,
                   resistance_gene_hits = ?,
                   virulence_gene_hits = ?,
                   resistance_mge_hits = ?,
                   virulence_mge_hits = ?,
                   updated_at = datetime('now')
             where sample_key = ?
            """,
            (
                summary["mlst_st"],
                summary["serotype_result"],
                summary["resistance_gene_hits"],
                summary["virulence_gene_hits"],
                summary["resistance_mge_hits"],
                summary["virulence_mge_hits"],
                sample_key,
            ),
        )
        updated += 1
        print(
            "UPDATED\t{sample}\t{mlst}\t{serotype}\t{arg}\t{vf}\t{arg_mge}\t{vf_mge}".format(
                sample=sample_name,
                mlst=summary["mlst_st"] or "-",
                serotype=summary["serotype_result"] or "-",
                arg=summary["resistance_gene_hits"] or "-",
                vf=summary["virulence_gene_hits"] or "-",
                arg_mge=summary["resistance_mge_hits"] or "-",
                vf_mge=summary["virulence_mge_hits"] or "-",
            )
        )
    conn.commit()
    conn.close()
    print(f"SUMMARY\tupdated={updated}\tskipped={skipped}")


if __name__ == "__main__":
    main()
