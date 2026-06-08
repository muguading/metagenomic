from __future__ import annotations

import json
import sqlite3
from datetime import date
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DB_PATH = PROJECT_ROOT / "bac_analysis_portal.sqlite3"

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

MONTH_CONFIGS = [
    ("2025-04", {"Acinetobacter baumannii": 1, "Escherichia coli": 2, "Klebsiella pneumoniae": 1, "Staphylococcus aureus": 1, "Pseudomonas aeruginosa": 1, "Neisseria meningitidis": 1}),
    ("2025-05", {"Acinetobacter baumannii": 1, "Escherichia coli": 2, "Klebsiella pneumoniae": 1, "Staphylococcus aureus": 1, "Pseudomonas aeruginosa": 1, "Neisseria meningitidis": 1}),
    ("2025-06", {"Acinetobacter baumannii": 1, "Escherichia coli": 1, "Klebsiella pneumoniae": 2, "Staphylococcus aureus": 1, "Pseudomonas aeruginosa": 1, "Neisseria meningitidis": 1}),
    ("2025-07", {"Acinetobacter baumannii": 2, "Escherichia coli": 1, "Klebsiella pneumoniae": 2, "Staphylococcus aureus": 1, "Pseudomonas aeruginosa": 1, "Neisseria meningitidis": 1}),
    ("2025-08", {"Acinetobacter baumannii": 2, "Escherichia coli": 1, "Klebsiella pneumoniae": 2, "Staphylococcus aureus": 1, "Pseudomonas aeruginosa": 2, "Neisseria meningitidis": 1}),
    ("2025-09", {"Acinetobacter baumannii": 2, "Escherichia coli": 1, "Klebsiella pneumoniae": 2, "Staphylococcus aureus": 1, "Pseudomonas aeruginosa": 2, "Neisseria meningitidis": 1}),
    ("2025-10", {"Acinetobacter baumannii": 2, "Escherichia coli": 1, "Klebsiella pneumoniae": 2, "Staphylococcus aureus": 2, "Pseudomonas aeruginosa": 1, "Neisseria meningitidis": 1}),
    ("2025-11", {"Acinetobacter baumannii": 2, "Escherichia coli": 1, "Klebsiella pneumoniae": 2, "Staphylococcus aureus": 2, "Pseudomonas aeruginosa": 1, "Neisseria meningitidis": 1}),
    ("2025-12", {"Acinetobacter baumannii": 3, "Escherichia coli": 1, "Klebsiella pneumoniae": 2, "Staphylococcus aureus": 1, "Pseudomonas aeruginosa": 1, "Neisseria meningitidis": 1}),
    ("2026-01", {"Acinetobacter baumannii": 3, "Escherichia coli": 1, "Klebsiella pneumoniae": 2, "Staphylococcus aureus": 1, "Pseudomonas aeruginosa": 1, "Neisseria meningitidis": 2}),
    ("2026-02", {"Acinetobacter baumannii": 3, "Escherichia coli": 1, "Klebsiella pneumoniae": 2, "Staphylococcus aureus": 1, "Pseudomonas aeruginosa": 1, "Neisseria meningitidis": 2}),
]

LOCATION_POOL = [
    ("江苏", "南京", "鼓楼区", "南京示例医院"),
    ("上海", "上海", "浦东新区", "上海示例医院"),
    ("浙江", "杭州", "西湖区", "杭州示例医院"),
]


def month_day(year_month: str, offset: int) -> str:
    year, month = map(int, year_month.split("-"))
    day = min(24, 3 + offset * 4)
    return date(year, month, day).isoformat()


def choose_location(index: int) -> tuple[str, str, str, str]:
    return LOCATION_POOL[index % len(LOCATION_POOL)]


def build_sample_key(scope: str, species_slug: str, year_month: str, seq: int) -> str:
    base = f"demo-history::{scope}::{species_slug}::{year_month.replace('-', '')}::{seq:02d}"
    return base if scope == "main" else f"personal::{base}"


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


def main() -> None:
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    created = 0
    updated = 0
    now_iso = "2026-04-05T09:00:00+08:00"

    for species_name, main_seed_key in SPECIES_SEEDS.items():
        main_seed = conn.execute("SELECT * FROM sample_library WHERE sample_key = ?", (main_seed_key,)).fetchone()
        personal_seed = conn.execute("SELECT * FROM sample_library WHERE sample_key = ?", (f"personal::{main_seed_key}",)).fetchone()
        if not main_seed or not personal_seed:
            continue
        species_slug = species_name.lower().replace(" ", "_")
        for year_month, distribution in MONTH_CONFIGS:
            sample_count = int(distribution.get(species_name, 0))
            for seq in range(1, sample_count + 1):
                for scope, seed in (("main", main_seed), ("personal", personal_seed)):
                    sample_key = build_sample_key(scope, species_slug, year_month, seq)
                    existed = conn.execute("SELECT 1 FROM sample_library WHERE sample_key = ?", (sample_key,)).fetchone() is not None
                    province, city, district, hospital = choose_location(seq + len(year_month) + (1 if scope == "personal" else 0))
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
                        task_id=f"demo-history-{year_month}-{species_slug}",
                        task_name="自动生成示例历史样本",
                        location_json=json.dumps(location_payload, ensure_ascii=False),
                        collection_date=collection_date,
                        note="自动生成示例历史样本，用于展示历史同期对比体系。",
                        description=f"{species_name} 历史演示样本（{year_month}）",
                        imported_at=now_iso,
                        updated_at=now_iso,
                        custom_metadata_json="[]",
                    )
                    upsert_row(conn, row)
                    if existed:
                        updated += 1
                    else:
                        created += 1
    conn.commit()
    print(json.dumps({"created": created, "updated": updated}, ensure_ascii=False))


if __name__ == "__main__":
    main()
