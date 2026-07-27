from __future__ import annotations

import json
import random
from datetime import date, timedelta
from pathlib import Path

from bac_analysis_portal.store import PortalStore, utc_now_iso


PROJECT_ROOT = Path(__file__).resolve().parents[1]
TARGET_COUNT = 400
KEY_PREFIX = "modeling-test-salmonella-"
RANDOM_SEED = 20260615

PROFILES = [
    {
        "serotype": "Enteritidis",
        "sts": ["ST11", "ST183"],
        "weight": 22,
        "risk": "重点关注",
        "amr": ["aac(6')-Iaa", "blaTEM-1B", "sul2"],
        "vf": ["invA", "sipB", "sopE", "spvC"],
    },
    {
        "serotype": "Typhimurium",
        "sts": ["ST19", "ST34", "ST313"],
        "weight": 22,
        "risk": "异常信号",
        "amr": ["aac(6')-Iaa", "blaTEM-1B", "tet(B)", "sul2", "floR"],
        "vf": ["invA", "sipB", "sopE2", "spvB", "spvC"],
    },
    {
        "serotype": "Derby",
        "sts": ["ST40", "ST71"],
        "weight": 12,
        "risk": "常规",
        "amr": ["aac(6')-Iaa", "tet(A)"],
        "vf": ["invA", "sipB", "sopB"],
    },
    {
        "serotype": "London",
        "sts": ["ST155", "ST166"],
        "weight": 8,
        "risk": "重点关注",
        "amr": ["aac(6')-Iaa", "blaTEM-1B", "sul1"],
        "vf": ["invA", "sipB", "sopE2"],
    },
    {
        "serotype": "Rissen",
        "sts": ["ST469", "ST887"],
        "weight": 8,
        "risk": "重点关注",
        "amr": ["aac(6')-Iaa", "tet(A)", "sul2"],
        "vf": ["invA", "sipB", "sopB"],
    },
    {
        "serotype": "Agona",
        "sts": ["ST13", "ST103"],
        "weight": 7,
        "risk": "常规",
        "amr": ["aac(6')-Iaa"],
        "vf": ["invA", "sipB"],
    },
    {
        "serotype": "Thompson",
        "sts": ["ST26", "ST128"],
        "weight": 6,
        "risk": "常规",
        "amr": ["aac(6')-Iaa", "sul2"],
        "vf": ["invA", "sipB", "sopE"],
    },
    {
        "serotype": "Choleraesuis",
        "sts": ["ST145", "ST68"],
        "weight": 5,
        "risk": "异常信号",
        "amr": ["aac(6')-Iaa", "blaTEM-1B", "floR", "sul2"],
        "vf": ["invA", "sipB", "spvB", "spvC"],
    },
    {
        "serotype": "Typhi",
        "sts": ["ST1", "ST2"],
        "weight": 5,
        "risk": "异常信号",
        "amr": ["aac(6')-Iaa", "blaTEM-1B", "qnrS1"],
        "vf": ["invA", "sipB", "tviB", "viaB"],
    },
    {
        "serotype": "Paratyphi A",
        "sts": ["ST85", "ST129"],
        "weight": 5,
        "risk": "异常信号",
        "amr": ["aac(6')-Iaa", "qnrS1"],
        "vf": ["invA", "sipB", "tviA"],
    },
]

LOCATIONS = [
    ("上海", "上海", "黄浦区", "黄浦区疾控中心"),
    ("上海", "上海", "浦东新区", "浦东新区疾控中心"),
    ("江苏", "南京", "鼓楼区", "南京市示例实验室"),
    ("浙江", "杭州", "西湖区", "杭州市示例实验室"),
    ("安徽", "合肥", "蜀山区", "合肥市示例实验室"),
    ("山东", "青岛", "市南区", "青岛市示例实验室"),
]

SOURCE_OPTIONS = [
    ("临床", "粪便/肛拭子", "患者"),
    ("食品", "食品", "食品样本"),
    ("环境", "环境拭子", "环境样本"),
    ("动物源", "组织/活检", "动物样本"),
]


def metadata_item(key: str, label: str, value: str, field_type: str = "text") -> dict[str, str]:
    return {"key": key, "label": label, "type": field_type, "value": value}


def choose_subset(rng: random.Random, values: list[str], minimum: int = 1) -> list[str]:
    upper = max(minimum, len(values))
    size = rng.randint(minimum, upper)
    return sorted(rng.sample(values, size))


def build_record(index: int, rng: random.Random) -> dict[str, object]:
    profile = rng.choices(PROFILES, weights=[item["weight"] for item in PROFILES], k=1)[0]
    province, city, district, lab = LOCATIONS[index % len(LOCATIONS)]
    source_category, sample_source, patient_group = SOURCE_OPTIONS[index % len(SOURCE_OPTIONS)]
    collection_date = date(2023, 1, 1) + timedelta(days=(index * 7 + rng.randint(0, 20)) % 1240)
    st = rng.choice(profile["sts"])
    amr = choose_subset(rng, profile["amr"])
    vf = choose_subset(rng, profile["vf"])
    risk = str(profile["risk"])
    if index % 20 == 0:
        risk = "排除"
    elif risk == "常规" and len(amr) >= 3:
        risk = "重点关注"
    q30 = round(rng.uniform(91.5, 99.2), 1)
    completeness = round(rng.uniform(91.0, 99.8), 1)
    contamination = round(rng.uniform(0.1, 4.2), 1)
    batch = f"SAL-{collection_date.year}-{(index % 12) + 1:02d}"
    project = f"SALMONELLA-SURVEILLANCE-{collection_date.year}"
    metadata = [
        metadata_item("modeling_review_label", "建模复核标签", risk, "select"),
        metadata_item("case_id", "病例/事件编号", f"SAL-CASE-{index:04d}"),
        metadata_item("project_id", "项目 ID", project),
        metadata_item("batch_id", "批次 ID", batch),
        metadata_item("hospital_or_lab", "医院/实验室", lab),
        metadata_item("source_category", "来源类型", source_category, "select"),
        metadata_item("patient_group", "人群/来源组", patient_group),
        metadata_item("age_group", "年龄组", rng.choice(["儿童", "青年", "中年", "老年", "不适用"])),
        metadata_item("cgMLST cluster", "cgMLST cluster", f"SC-{profile['serotype'][:3].upper()}-{(index % 18) + 1:02d}"),
        metadata_item("SNP cluster", "SNP cluster", f"SNP-SAL-{(index % 25) + 1:02d}"),
        metadata_item("phylogenetic cluster", "系统发育簇", f"PHY-SAL-{(index % 14) + 1:02d}"),
        metadata_item("review_note", "复核说明", f"沙门菌模型训练测试数据：{risk}"),
    ]
    now = utc_now_iso()
    return {
        "sample_key": f"{KEY_PREFIX}{index:04d}",
        "genome_id": f"SAL-GENOME-{index:04d}",
        "sample_name": f"salmonella_model_train_{index:04d}",
        "task_id": "modeling-test-salmonella-seed",
        "task_name": "沙门菌建模训练测试样本",
        "owner": "admin",
        "owner_group": "training",
        "report_dir": "",
        "output_dir": "",
        "final_fasta_path": "",
        "species_name": "Salmonella enterica",
        "pathogen_type": "bacteria",
        "taxid": "28901",
        "mlst_species_name": "Salmonella enterica",
        "mlst_st": st,
        "serotype_result": profile["serotype"],
        "genome_length": str(rng.randint(4_650_000, 5_150_000)),
        "q20_rate": str(round(min(99.9, q30 + rng.uniform(0.3, 1.5)), 1)),
        "q30_rate": str(q30),
        "completeness": str(completeness),
        "contamination": str(contamination),
        "contig_count": str(rng.randint(28, 185)),
        "plasmid_count": str(rng.randint(0, 7)),
        "total_length": str(rng.randint(4_650_000, 5_150_000)),
        "resistance_count": str(len(amr)),
        "virulence_count": str(len(vf)),
        "resistance_gene_hits": ";".join(amr),
        "virulence_gene_hits": ";".join(vf),
        "resistance_mge_hits": rng.choice(["IncFIB(S)", "IncHI2", "IncI1-I", ""]),
        "virulence_mge_hits": rng.choice(["SPI-1", "SPI-2", "spv-region", ""]),
        "description": f"沙门菌建模测试样本，血清型 {profile['serotype']}，{st}",
        "gender": rng.choice(["男", "女", "不适用", ""]),
        "country": "中国",
        "location_json": json.dumps(
            {"province": province, "city": city, "district": district, "detail": lab},
            ensure_ascii=False,
        ),
        "sample_type": sample_source,
        "sequencing_method": rng.choice(["Illumina", "MGI", "Nanopore+Illumina"]),
        "custom_metadata_json": json.dumps(metadata, ensure_ascii=False),
        "sample_alias": f"SAL-{profile['serotype']}-{index:04d}",
        "sample_source": sample_source,
        "collection_date": collection_date.isoformat(),
        "host_info": patient_group,
        "note": "仅用于沙门菌模型训练与预测功能测试。",
        "library_scope": "main",
        "visibility_scope": "group",
        "source_submission_id": "",
        "imported_at": now,
        "updated_at": now,
    }


def main() -> None:
    rng = random.Random(RANDOM_SEED)
    store = PortalStore.from_project_root(PROJECT_ROOT)
    created = 0
    updated = 0
    for index in range(1, TARGET_COUNT + 1):
        key = f"{KEY_PREFIX}{index:04d}"
        try:
            store.get_sample_library_record(key)
            updated += 1
        except KeyError:
            created += 1
        store.upsert_sample_library_record(build_record(index, rng))
    total = sum(
        1
        for row in store.list_sample_library()
        if str(row.get("species_name") or "") == "Salmonella enterica"
        and str(row.get("sample_key") or "").startswith("modeling-test-")
    )
    print(json.dumps({"created": created, "updated": updated, "generated": TARGET_COUNT, "salmonella_modeling_total": total}, ensure_ascii=False))


if __name__ == "__main__":
    main()
