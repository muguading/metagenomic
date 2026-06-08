from __future__ import annotations

import json
import sqlite3
from datetime import date
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DB_PATH = PROJECT_ROOT / "bac_analysis_portal.sqlite3"
TARGET_TOTAL = 300
NOW_ISO = "2026-04-05T10:30:00+08:00"
HISTORY_PREFIX = "demo-history-plus"

SAMPLE_LIBRARY_COLUMNS = [
    "sample_key",
    "genome_id",
    "sample_name",
    "task_id",
    "task_name",
    "owner",
    "owner_group",
    "report_dir",
    "output_dir",
    "final_fasta_path",
    "species_name",
    "taxid",
    "mlst_species_name",
    "mlst_st",
    "serotype_result",
    "genome_length",
    "q20_rate",
    "q30_rate",
    "completeness",
    "contamination",
    "contig_count",
    "plasmid_count",
    "total_length",
    "resistance_count",
    "virulence_count",
    "resistance_gene_hits",
    "virulence_gene_hits",
    "resistance_mge_hits",
    "virulence_mge_hits",
    "description",
    "gender",
    "country",
    "location_json",
    "sample_type",
    "sequencing_method",
    "custom_metadata_json",
    "sample_alias",
    "sample_source",
    "collection_date",
    "host_info",
    "note",
    "library_scope",
    "visibility_scope",
    "source_submission_id",
    "imported_at",
    "updated_at",
]

SPECIES_SEEDS = {
    "Acinetobacter baumannii": "local::/Users/wuhhh/Downloads/ncbi_dataset (3)/ncbi_dataset/data/GCA_001272635.1/GCA_001272635.1_ASM127263v1_genomic.fna",
    "Escherichia coli": "local::/Users/wuhhh/Downloads/ncbi_dataset (8)/ncbi_dataset/data/GCF_000203795.2/GCF_000203795.2_v1.0_genomic.fna",
    "Klebsiella pneumoniae": "local::/Users/wuhhh/Downloads/ncbi_dataset (7)/ncbi_dataset/data/GCF_000864005.1/GCF_000864005.1_ViralProj15520_genomic.fna",
    "Staphylococcus aureus": "local::/Users/wuhhh/Downloads/ncbi_dataset (20)/ncbi_dataset/data/GCF_002290485.1/GCF_002290485.1_ASM229048v1_genomic.fna",
    "Pseudomonas aeruginosa": "local::/Users/wuhhh/Downloads/ncbi_dataset (3)/ncbi_dataset/data/GCF_001272635.1/GCF_001272635.1_ASM127263v1_genomic.fna",
    "Neisseria meningitidis": "local::/Users/wuhhh/Downloads/ncbi_dataset (20)/ncbi_dataset/data/GCF_022354085.1/GCF_022354085.1_ASM2235408v1_genomic.fna",
}

EXTRA_MONTHS = [
    "2025-04",
    "2025-05",
    "2025-06",
    "2025-07",
    "2025-08",
    "2025-09",
    "2025-10",
    "2025-11",
    "2025-12",
    "2026-01",
    "2026-02",
]

SPECIES_ROTATION = [
    "Acinetobacter baumannii",
    "Escherichia coli",
    "Klebsiella pneumoniae",
    "Staphylococcus aureus",
    "Pseudomonas aeruginosa",
    "Neisseria meningitidis",
]

LOCATION_POOL = [
    ("江苏", "南京", "鼓楼区", "南京示例医院"),
    ("上海", "上海", "浦东新区", "上海示例医院"),
    ("浙江", "杭州", "西湖区", "杭州示例医院"),
]


def choose_location(index: int) -> tuple[str, str, str, str]:
    return LOCATION_POOL[index % len(LOCATION_POOL)]


def month_day(year_month: str, offset: int) -> str:
    year, month = map(int, year_month.split("-"))
    day = min(26, 5 + offset * 3)
    return date(year, month, day).isoformat()


def build_sample_name(seed_name: str, year_month: str, seq: int) -> str:
    stem = seed_name if seed_name.startswith(("GCF_", "GCA_", "RNA")) else f"GCF_{seed_name}"
    return f"{stem}_{year_month.replace('-', '')}_{seq:02d}"


def upsert_row(conn: sqlite3.Connection, row: dict[str, str]) -> None:
    placeholders = ", ".join(["?"] * len(SAMPLE_LIBRARY_COLUMNS))
    assignments = ", ".join([f"{column}=excluded.{column}" for column in SAMPLE_LIBRARY_COLUMNS[1:]])
    sql = f"""
        INSERT INTO sample_library ({", ".join(SAMPLE_LIBRARY_COLUMNS)})
        VALUES ({placeholders})
        ON CONFLICT(sample_key) DO UPDATE SET
        {assignments}
    """
    conn.execute(sql, [row.get(column, "") for column in SAMPLE_LIBRARY_COLUMNS])


def get_next_seq(conn: sqlite3.Connection, scope: str, species_slug: str, year_month: str) -> int:
    prefix = (
        f"{HISTORY_PREFIX}::{species_slug}::{year_month.replace('-', '')}::"
        if scope == "main"
        else f"personal::{HISTORY_PREFIX}::{species_slug}::{year_month.replace('-', '')}::"
    )
    row = conn.execute(
        """
        SELECT sample_key
          FROM sample_library
         WHERE sample_key LIKE ?
         ORDER BY sample_key DESC
         LIMIT 1
        """,
        (f"{prefix}%",),
    ).fetchone()
    if not row:
        return 1
    return int(str(row[0]).rsplit("::", 1)[-1]) + 1


def build_sample_key(scope: str, species_slug: str, year_month: str, seq: int) -> str:
    base = f"{HISTORY_PREFIX}::{species_slug}::{year_month.replace('-', '')}::{seq:02d}"
    return base if scope == "main" else f"personal::{base}"


def main() -> None:
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    total = conn.execute("SELECT COUNT(*) FROM sample_library").fetchone()[0]
    if total >= TARGET_TOTAL:
        print(json.dumps({"created": 0, "updated": 0, "total": total}, ensure_ascii=False))
        return

    needed = TARGET_TOTAL - total
    created = 0
    updated = 0
    rotation_index = 0

    while needed >= 2:
        year_month = EXTRA_MONTHS[(created // 2) % len(EXTRA_MONTHS)]
        species_name = SPECIES_ROTATION[rotation_index % len(SPECIES_ROTATION)]
        rotation_index += 1
        species_slug = species_name.lower().replace(" ", "_")
        main_seed_key = SPECIES_SEEDS[species_name]
        main_seed = conn.execute("SELECT * FROM sample_library WHERE sample_key = ?", (main_seed_key,)).fetchone()
        personal_seed = conn.execute("SELECT * FROM sample_library WHERE sample_key = ?", (f"personal::{main_seed_key}",)).fetchone()
        if not main_seed or not personal_seed:
            continue

        main_seq = get_next_seq(conn, "main", species_slug, year_month)
        for scope, seed, seq in (
            ("main", main_seed, main_seq),
            ("personal", personal_seed, get_next_seq(conn, "personal", species_slug, year_month)),
        ):
            sample_key = build_sample_key(scope, species_slug, year_month, seq)
            existed = conn.execute("SELECT 1 FROM sample_library WHERE sample_key = ?", (sample_key,)).fetchone() is not None
            province, city, district, hospital = choose_location(seq + rotation_index + (1 if scope == "personal" else 0))
            collection_date = month_day(year_month, seq)
            location_payload = {
                "province": province,
                "city": city,
                "district": district,
                "detail": f"{hospital}{collection_date[-2:]}号采样点",
            }
            row = {column: seed[column] for column in SAMPLE_LIBRARY_COLUMNS}
            row.update(
                sample_key=sample_key,
                sample_name=build_sample_name(seed["sample_name"], year_month, seq),
                task_id=f"demo-history-extra-{year_month}-{species_slug}",
                task_name="扩充示例历史样本",
                location_json=json.dumps(location_payload, ensure_ascii=False),
                collection_date=collection_date,
                note="扩充示例历史样本，用于展示年度历史同期与专题预警。",
                description=f"{species_name} 扩充演示样本（{year_month}）",
                imported_at=NOW_ISO,
                updated_at=NOW_ISO,
                custom_metadata_json="[]",
            )
            upsert_row(conn, row)
            if existed:
                updated += 1
            else:
                created += 1
                needed -= 1
                if needed <= 0:
                    break

    if needed == 1:
        species_name = "Neisseria meningitidis"
        year_month = "2026-02"
        species_slug = species_name.lower().replace(" ", "_")
        main_seed = conn.execute("SELECT * FROM sample_library WHERE sample_key = ?", (SPECIES_SEEDS[species_name],)).fetchone()
        if main_seed:
            seq = get_next_seq(conn, "main", species_slug, year_month)
            sample_key = build_sample_key("main", species_slug, year_month, seq)
            province, city, district, hospital = choose_location(seq + 99)
            collection_date = month_day(year_month, seq)
            location_payload = {
                "province": province,
                "city": city,
                "district": district,
                "detail": f"{hospital}{collection_date[-2:]}号采样点",
            }
            row = {column: main_seed[column] for column in SAMPLE_LIBRARY_COLUMNS}
            row.update(
                sample_key=sample_key,
                sample_name=build_sample_name(main_seed["sample_name"], year_month, seq),
                task_id=f"demo-history-extra-{year_month}-{species_slug}",
                task_name="扩充示例历史样本",
                location_json=json.dumps(location_payload, ensure_ascii=False),
                collection_date=collection_date,
                note="扩充示例历史样本，用于展示年度历史同期与专题预警。",
                description=f"{species_name} 扩充演示样本（{year_month}）",
                imported_at=NOW_ISO,
                updated_at=NOW_ISO,
                custom_metadata_json="[]",
            )
            upsert_row(conn, row)
            created += 1

    conn.commit()
    final_total = conn.execute("SELECT COUNT(*) FROM sample_library").fetchone()[0]
    conn.close()
    print(json.dumps({"created": created, "updated": updated, "total": final_total}, ensure_ascii=False))


if __name__ == "__main__":
    main()
