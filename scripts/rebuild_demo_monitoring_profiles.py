from __future__ import annotations

import json
import sqlite3
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DB_PATH = ROOT / "bac_analysis_portal.sqlite3"
NOW_ISO = "2026-04-05T12:10:00+08:00"

SAMPLE_COLUMNS = [
    "sample_key",
    "species_name",
    "taxid",
    "mlst_species_name",
    "mlst_st",
    "serotype_result",
    "resistance_count",
    "virulence_count",
    "resistance_gene_hits",
    "virulence_gene_hits",
    "resistance_mge_hits",
    "virulence_mge_hits",
    "description",
    "country",
    "location_json",
    "sample_type",
    "sequencing_method",
    "custom_metadata_json",
    "sample_source",
    "collection_date",
    "host_info",
    "note",
    "updated_at",
]

SPECIES_PROFILES = {
    "Acinetobacter baumannii": {
        "taxid": "470",
        "syndrome": "医院感染监测",
        "sources": ["ICU", "住院", "ICU", "住院", "ICU"],
        "specimen_category": "呼吸道",
        "sample_source": "呼吸道标本",
        "sample_type": "fasta",
        "sequencing_method": "spades",
        "host_info": "住院患者",
        "st_serotype": [("ST2", "KL2"), ("ST2", "KL2"), ("ST208", "KL2"), ("ST2", "KL9"), ("ST191", "KL2")],
        "arg": "OXA-23、blaADC-25、aph(3')-Ia",
        "vf": "ompA、bap、csuE",
        "arg_mge": "OXA-23（转座子）、aph(3')-Ia（整合子）",
        "vf_mge": "bap（质粒）",
        "resistance_count": "3",
        "virulence_count": "3",
    },
    "Klebsiella pneumoniae": {
        "taxid": "573",
        "syndrome": "医院感染监测",
        "sources": ["住院", "ICU", "住院", "ICU", "住院"],
        "specimen_category": "呼吸道",
        "sample_source": "呼吸道标本",
        "sample_type": "fasta",
        "sequencing_method": "spades",
        "host_info": "住院患者",
        "st_serotype": [("ST11", "K64"), ("ST11", "K64"), ("ST11", "K47"), ("ST15", "K64"), ("ST23", "K1")],
        "arg": "blaKPC-2、blaCTX-M-65、rmtB",
        "vf": "rmpA、iucA、ybtS",
        "arg_mge": "blaKPC-2（质粒）、rmtB（质粒）",
        "vf_mge": "rmpA（质粒）、iucA（质粒）",
        "resistance_count": "3",
        "virulence_count": "3",
    },
    "Pseudomonas aeruginosa": {
        "taxid": "287",
        "syndrome": "医院感染监测",
        "sources": ["ICU", "住院", "ICU", "住院"],
        "specimen_category": "呼吸道",
        "sample_source": "呼吸道标本",
        "sample_type": "fasta",
        "sequencing_method": "spades",
        "host_info": "住院患者",
        "st_serotype": [("ST235", "O11"), ("ST235", "O11"), ("ST244", "O11"), ("ST773", "O12")],
        "arg": "blaVIM-2、aacA4、fosA",
        "vf": "exoU、lasB、toxA",
        "arg_mge": "blaVIM-2（整合子）、aacA4（整合子）",
        "vf_mge": "exoU（基因岛）",
        "resistance_count": "3",
        "virulence_count": "3",
    },
    "Staphylococcus aureus": {
        "taxid": "1280",
        "syndrome": "败血症/血流感染",
        "sources": ["住院", "ICU", "住院", "门急诊"],
        "specimen_category": "血液",
        "sample_source": "血培养",
        "sample_type": "fasta",
        "sequencing_method": "spades",
        "host_info": "住院患者",
        "st_serotype": [("ST239", "Capsule type 5"), ("ST239", "Capsule type 5"), ("ST59", "Capsule type 8"), ("ST5", "Capsule type 5")],
        "arg": "mecA、ermC、tetK",
        "vf": "clfA、hla、eta",
        "arg_mge": "ermC（质粒）、tetK（质粒）",
        "vf_mge": "eta（噬菌体）",
        "resistance_count": "3",
        "virulence_count": "3",
    },
    "Escherichia coli": {
        "taxid": "562",
        "syndrome": "泌尿生殖道感染",
        "sources": ["门急诊", "住院", "社区监测", "门急诊"],
        "specimen_category": "尿液",
        "sample_source": "尿液",
        "sample_type": "fasta",
        "sequencing_method": "spades",
        "host_info": "门诊患者",
        "st_serotype": [("ST131", "O25b:H4"), ("ST131", "O25b:H4"), ("ST69", "O15:H1"), ("ST1193", "O75:H5"), ("ST95", "O1:K1")],
        "arg": "blaCTX-M-15、aac(6')-Ib-cr、qnrS1",
        "vf": "fimH、iutA、papC",
        "arg_mge": "blaCTX-M-15（质粒）、qnrS1（质粒）",
        "vf_mge": "iutA（质粒）",
        "resistance_count": "3",
        "virulence_count": "3",
    },
    "Neisseria meningitidis": {
        "taxid": "487",
        "syndrome": "脑膜炎/脑膜脑炎",
        "sources": ["发热门诊", "住院", "发热门诊", "住院"],
        "specimen_category": "脑脊液",
        "sample_source": "脑脊液",
        "sample_type": "fasta",
        "sequencing_method": "spades",
        "host_info": "侵袭性感染病例",
        "st_serotype": [("ST4821", "C"), ("ST4821", "B"), ("ST11", "C"), ("ST41/44", "B")],
        "arg": "penA、rpoB、gyrA",
        "vf": "fHbp、lbpA、lbpB",
        "arg_mge": "macB（复制/重组）、mtrR（修复）",
        "vf_mge": "fHbp（重组）、lbpA（重组）、lbpB（防御）",
        "resistance_count": "3",
        "virulence_count": "3",
    },
    "Streptococcus pneumoniae": {
        "taxid": "1313",
        "syndrome": "脑膜炎/脑膜脑炎",
        "sources": ["住院", "发热门诊", "住院"],
        "specimen_category": "脑脊液",
        "sample_source": "脑脊液",
        "sample_type": "fasta",
        "sequencing_method": "spades",
        "host_info": "侵袭性感染病例",
        "st_serotype": [("ST320", "19A"), ("ST271", "19F"), ("ST180", "3")],
        "arg": "pbp2x、ermB、tetM",
        "vf": "pspC、ply、cpsA",
        "arg_mge": "ermB（转座子）、tetM（转座子）",
        "vf_mge": "ply（重组）",
        "resistance_count": "3",
        "virulence_count": "3",
    },
    "Haemophilus influenzae": {
        "taxid": "727",
        "syndrome": "脑膜炎/脑膜脑炎",
        "sources": ["住院", "发热门诊"],
        "specimen_category": "脑脊液",
        "sample_source": "脑脊液",
        "sample_type": "fasta",
        "sequencing_method": "spades",
        "host_info": "侵袭性感染病例",
        "st_serotype": [("ST6", "b"), ("ST103", "非分型")],
        "arg": "blaTEM-1、ftsI、acrB",
        "vf": "iga1、hmw1A、ompP5",
        "arg_mge": "blaTEM-1（质粒）",
        "vf_mge": "hmw1A（重组）",
        "resistance_count": "3",
        "virulence_count": "3",
    },
    "Salmonella enterica": {
        "taxid": "28901",
        "syndrome": "腹泻/肠道感染",
        "sources": ["门急诊", "社区监测", "食品监测"],
        "specimen_category": "粪便/肛拭子",
        "sample_source": "粪便",
        "sample_type": "fasta",
        "sequencing_method": "spades",
        "host_info": "腹泻病例",
        "st_serotype": [("ST11", "Enteritidis"), ("ST34", "Typhimurium"), ("ST19", "Typhimurium")],
        "arg": "blaTEM-1、sul2、tetA",
        "vf": "invA、sipB、spvC",
        "arg_mge": "blaTEM-1（质粒）、tetA（质粒）",
        "vf_mge": "spvC（毒力质粒）",
        "resistance_count": "3",
        "virulence_count": "3",
    },
    "Shigella spp.": {
        "taxid": "620",
        "syndrome": "腹泻/肠道感染",
        "sources": ["门急诊", "社区监测", "门急诊"],
        "specimen_category": "粪便/肛拭子",
        "sample_source": "粪便",
        "sample_type": "fasta",
        "sequencing_method": "spades",
        "host_info": "腹泻病例",
        "st_serotype": [("ST100", "flexneri 2a"), ("ST152", "sonnei phase II")],
        "arg": "qnrS1、blaTEM-1、dfrA14",
        "vf": "ipaH、virA、icsA",
        "arg_mge": "qnrS1（质粒）、dfrA14（整合子）",
        "vf_mge": "virA（毒力质粒）",
        "resistance_count": "3",
        "virulence_count": "3",
    },
    "Vibrio parahaemolyticus": {
        "taxid": "691",
        "syndrome": "腹泻/肠道感染",
        "sources": ["门急诊", "食品监测", "社区监测"],
        "specimen_category": "粪便/肛拭子",
        "sample_source": "粪便",
        "sample_type": "fasta",
        "sequencing_method": "spades",
        "host_info": "腹泻病例",
        "st_serotype": [("ST3", "O3:K6"), ("ST36", "O4:K12")],
        "arg": "qnrVC6、tet(34)、blaCARB",
        "vf": "tdh、trh、toxR",
        "arg_mge": "qnrVC6（整合子）",
        "vf_mge": "tdh（重组）",
        "resistance_count": "3",
        "virulence_count": "3",
    },
    "Enterococcus faecium": {
        "taxid": "1352",
        "syndrome": "败血症/血流感染",
        "sources": ["住院", "ICU"],
        "specimen_category": "血液",
        "sample_source": "血培养",
        "sample_type": "fasta",
        "sequencing_method": "spades",
        "host_info": "住院患者",
        "st_serotype": [("ST78", "无血清型"), ("ST80", "无血清型"), ("ST17", "无血清型")],
        "arg": "vanA、vanB、aac(6')-Ii",
        "vf": "esp、hyl、acm",
        "arg_mge": "vanA（转座子）、vanB（整合元件）",
        "vf_mge": "esp（质粒）",
        "resistance_count": "3",
        "virulence_count": "3",
    },
    "Campylobacter jejuni": {
        "taxid": "197",
        "syndrome": "腹泻/肠道感染",
        "sources": ["门急诊", "社区监测", "食品监测"],
        "specimen_category": "粪便/肛拭子",
        "sample_source": "粪便",
        "sample_type": "fasta",
        "sequencing_method": "spades",
        "host_info": "腹泻病例",
        "st_serotype": [("ST21", "HS:2"), ("ST50", "HS:4 complex"), ("ST45", "HS:1/44")],
        "arg": "tetO、gyrA、blaOXA-61",
        "vf": "cadF、ciaB、cdtA",
        "arg_mge": "tetO（质粒）",
        "vf_mge": "",
        "resistance_count": "3",
        "virulence_count": "3",
    },
    "Enterobacter cloacae": {
        "taxid": "550",
        "syndrome": "医院感染监测",
        "sources": ["住院", "ICU", "住院"],
        "specimen_category": "呼吸道",
        "sample_source": "呼吸道标本",
        "sample_type": "fasta",
        "sequencing_method": "spades",
        "host_info": "住院患者",
        "st_serotype": [("ST114", "无血清型"), ("ST171", "无血清型"), ("ST78", "无血清型")],
        "arg": "blaACT-7、fosA2、qnrS1",
        "vf": "ompX、entB、fepA",
        "arg_mge": "qnrS1（质粒）",
        "vf_mge": "",
        "resistance_count": "3",
        "virulence_count": "3",
    },
}

HOSPITAL_LOCATIONS = [
    ("江苏", "南京", "鼓楼区", "南京示例医院ICU监测点"),
    ("江苏", "苏州", "工业园区", "苏州示例医院呼吸病区"),
    ("江苏", "无锡", "梁溪区", "无锡示例医院ICU监测点"),
    ("江苏", "常州", "天宁区", "常州示例医院感染病区"),
    ("上海", "上海", "黄浦区", "黄浦示例医院重症病区"),
    ("上海", "上海", "浦东新区", "浦东示例医院感染病区"),
    ("上海", "上海", "徐汇区", "徐汇示例医院ICU监测点"),
    ("上海", "上海", "静安区", "静安示例医院呼吸病区"),
    ("浙江", "杭州", "西湖区", "杭州示例医院ICU监测点"),
    ("浙江", "宁波", "鄞州区", "宁波示例医院重症病区"),
    ("浙江", "温州", "鹿城区", "温州示例医院感染病区"),
    ("浙江", "绍兴", "越城区", "绍兴示例医院ICU监测点"),
]

DIARRHEAL_LOCATIONS = [
    ("江苏", "南京", "玄武区", "南京社区腹泻监测点"),
    ("江苏", "苏州", "姑苏区", "苏州肠道门诊采样点"),
    ("江苏", "无锡", "滨湖区", "无锡肠道门诊采样点"),
    ("上海", "上海", "徐汇区", "上海肠道门诊采样点"),
    ("上海", "上海", "闵行区", "闵行食品相关采样点"),
    ("上海", "上海", "宝山区", "宝山腹泻监测点"),
    ("浙江", "杭州", "滨江区", "杭州肠道门诊采样点"),
    ("浙江", "温州", "鹿城区", "温州肠道门诊采样点"),
    ("浙江", "嘉兴", "南湖区", "嘉兴腹泻监测点"),
]

MENINGITIS_LOCATIONS = [
    ("上海", "上海", "黄浦区", "上海脑膜炎监测点"),
    ("江苏", "南京", "鼓楼区", "南京神经感染监测点"),
    ("浙江", "杭州", "西湖区", "杭州脑膜炎监测点"),
    ("江苏", "苏州", "工业园区", "苏州神经感染监测点"),
    ("浙江", "宁波", "海曙区", "宁波神经感染监测点"),
    ("上海", "上海", "杨浦区", "杨浦神经感染监测点"),
    ("江苏", "无锡", "梁溪区", "无锡脑膜炎监测点"),
]

ENV_LOCATIONS = [
    ("江苏", "苏州", "工业园区", "苏州污水处理点"),
    ("上海", "上海", "浦东新区", "浦东环境水体点"),
    ("浙江", "杭州", "滨江区", "杭州食品相关采样点"),
    ("浙江", "嘉兴", "南湖区", "嘉兴环境监测点"),
    ("江苏", "南京", "建邺区", "南京污水监测点"),
    ("上海", "上海", "宝山区", "宝山食品环境点"),
]

OTHER_LOCATIONS = [
    ("安徽", "合肥", "包河区", "合肥协作医院采样点"),
    ("北京", "北京市", "东城区", "国家示范协作点"),
]

SCENARIO_LIBRARY = {
    "hai_ab": ("Acinetobacter baumannii", HOSPITAL_LOCATIONS),
    "hai_kp": ("Klebsiella pneumoniae", HOSPITAL_LOCATIONS),
    "hai_pa": ("Pseudomonas aeruginosa", HOSPITAL_LOCATIONS),
    "hai_ec": ("Enterobacter cloacae", HOSPITAL_LOCATIONS),
    "bsi_sa": ("Staphylococcus aureus", HOSPITAL_LOCATIONS),
    "bsi_efm": ("Enterococcus faecium", HOSPITAL_LOCATIONS),
    "uti_ec": ("Escherichia coli", HOSPITAL_LOCATIONS),
    "uti_kp": ("Klebsiella pneumoniae", HOSPITAL_LOCATIONS),
    "men_nm": ("Neisseria meningitidis", MENINGITIS_LOCATIONS),
    "men_spn": ("Streptococcus pneumoniae", MENINGITIS_LOCATIONS),
    "men_hi": ("Haemophilus influenzae", MENINGITIS_LOCATIONS),
    "gut_sal": ("Salmonella enterica", DIARRHEAL_LOCATIONS),
    "gut_shi": ("Shigella spp.", DIARRHEAL_LOCATIONS),
    "gut_vp": ("Vibrio parahaemolyticus", DIARRHEAL_LOCATIONS),
    "gut_cj": ("Campylobacter jejuni", DIARRHEAL_LOCATIONS),
    "env_vp": ("Vibrio parahaemolyticus", ENV_LOCATIONS),
    "env_ec": ("Escherichia coli", ENV_LOCATIONS),
}

MONTH_SCENARIOS = {
    4: ["hai_ab", "hai_kp", "uti_ec", "bsi_sa", "men_nm", "gut_sal", "hai_pa", "gut_shi", "men_spn", "hai_ec", "bsi_efm", "uti_kp"],
    5: ["hai_ab", "hai_kp", "uti_ec", "bsi_sa", "men_nm", "gut_sal", "gut_shi", "hai_pa", "men_hi", "hai_ec", "bsi_efm", "uti_kp"],
    6: ["gut_sal", "gut_shi", "gut_vp", "gut_cj", "uti_ec", "hai_kp", "hai_ab", "env_ec", "gut_sal", "bsi_sa", "uti_kp", "hai_pa"],
    7: ["gut_sal", "gut_shi", "gut_vp", "env_vp", "gut_cj", "hai_kp", "hai_ab", "env_ec", "gut_sal", "bsi_sa", "uti_ec", "hai_ec"],
    8: ["gut_vp", "gut_sal", "gut_shi", "env_vp", "gut_cj", "hai_ab", "hai_kp", "env_ec", "men_nm", "bsi_sa", "uti_ec", "hai_pa"],
    9: ["gut_sal", "gut_shi", "gut_vp", "gut_cj", "uti_ec", "hai_kp", "hai_ab", "env_vp", "gut_sal", "men_spn", "bsi_sa", "uti_kp"],
    10: ["hai_ab", "hai_kp", "uti_ec", "bsi_sa", "men_nm", "men_spn", "gut_sal", "hai_pa", "bsi_efm", "hai_ec", "uti_kp", "men_hi"],
    11: ["hai_ab", "hai_kp", "men_nm", "men_spn", "hai_pa", "bsi_sa", "uti_ec", "bsi_efm", "men_hi", "hai_ec", "uti_kp", "hai_ab"],
    12: ["hai_ab", "hai_kp", "men_nm", "men_spn", "hai_pa", "bsi_sa", "uti_ec", "bsi_efm", "hai_ec", "men_hi", "uti_kp", "hai_ab"],
    1: ["hai_ab", "hai_kp", "men_nm", "men_spn", "hai_pa", "bsi_sa", "uti_ec", "bsi_efm", "hai_ec", "men_nm", "uti_kp", "men_hi"],
    2: ["hai_ab", "hai_kp", "men_nm", "men_spn", "hai_pa", "bsi_sa", "uti_ec", "bsi_efm", "hai_ec", "men_hi", "uti_kp", "hai_ab"],
    3: ["hai_ab", "hai_kp", "men_nm", "men_spn", "hai_pa", "bsi_sa", "uti_ec", "gut_sal", "hai_ec", "men_nm", "uti_kp", "bsi_efm"],
}

QUARTER_SCENARIOS = {
    2: [
        "gut_sal", "gut_shi", "gut_vp", "gut_cj",
        "uti_ec", "uti_kp",
        "hai_kp", "hai_ab", "bsi_sa",
        "env_ec", "env_vp",
        "men_nm",
    ],
    3: [
        "gut_vp", "gut_sal", "gut_shi", "gut_cj",
        "env_vp", "env_ec",
        "uti_ec",
        "hai_kp", "hai_ab",
        "bsi_sa",
        "men_nm",
        "gut_vp",
    ],
    4: [
        "hai_ab", "hai_kp", "hai_pa", "hai_ec",
        "bsi_sa", "bsi_efm",
        "men_nm", "men_spn", "men_hi",
        "uti_ec", "uti_kp",
        "hai_ab",
    ],
    1: [
        "hai_ab", "hai_kp", "hai_pa", "hai_ec",
        "men_nm", "men_spn", "men_hi",
        "bsi_sa", "bsi_efm",
        "uti_ec",
        "hai_kp",
        "men_nm",
    ],
}

QUARTER_LOCATION_BIAS = {
    2: ["江苏", "江苏", "上海", "浙江", "浙江"],
    3: ["上海", "浙江", "浙江", "上海", "江苏"],
    4: ["浙江", "江苏", "上海", "江苏", "浙江"],
    1: ["上海", "上海", "浙江", "江苏", "浙江"],
}


def build_metadata(
    *,
    case_id: str,
    patient_id: str,
    source: str,
    syndrome: str,
    collection_unit: str,
    submitting_unit: str,
    ward_department: str,
    specimen_category: str,
    location: tuple[str, str, str, str],
    cluster_status: str,
    epidemiology_link: str,
    traditional_result: str,
) -> str:
    province, city, district, detail = location
    payload = [
        {"key": "case_id", "label": "病例/事件编号", "type": "text", "value": case_id},
        {"key": "patient_id", "label": "患者/个案编号", "type": "text", "value": patient_id},
        {"key": "surveillance_source", "label": "监测来源", "type": "select", "value": source},
        {"key": "suspected_syndrome", "label": "疑似症候群", "type": "select", "value": syndrome},
        {"key": "submitting_unit", "label": "送检单位", "type": "text", "value": submitting_unit},
        {"key": "collection_unit", "label": "采样单位", "type": "text", "value": collection_unit},
        {"key": "ward_department", "label": "科室/病区", "type": "text", "value": ward_department},
        {"key": "specimen_category", "label": "标本类别", "type": "select", "value": specimen_category},
        {
            "key": "collection_site",
            "label": "采样地点",
            "type": "location",
            "value": {"province": province, "city": city, "district": district, "detail": detail},
        },
        {"key": "cluster_status", "label": "聚集性状态", "type": "select", "value": cluster_status},
        {"key": "epidemiology_link", "label": "流行病学关联", "type": "text", "value": epidemiology_link},
        {"key": "traditional_result", "label": "传统检测结果", "type": "text", "value": traditional_result},
    ]
    return json.dumps(payload, ensure_ascii=False)


def choose_profile(month: int, index: int) -> str:
    quarter = ((month - 1) // 3) + 1
    if index % 5 == 0:
        scenarios = QUARTER_SCENARIOS.get(quarter) or MONTH_SCENARIOS.get(month) or MONTH_SCENARIOS[3]
    else:
        scenarios = MONTH_SCENARIOS.get(month) or QUARTER_SCENARIOS.get(quarter) or MONTH_SCENARIOS[3]
    return scenarios[index % len(scenarios)]


def choose_location(pool: list[tuple[str, str, str, str]], idx: int, month: int) -> tuple[str, str, str, str]:
    quarter = ((month - 1) // 3) + 1
    province_bias = QUARTER_LOCATION_BIAS.get(quarter) or []
    if province_bias:
        target_province = province_bias[idx % len(province_bias)]
        filtered = [item for item in pool if item[0] == target_province]
        if filtered:
            location = filtered[(idx // max(1, len(province_bias))) % len(filtered)]
        else:
            location = pool[idx % len(pool)]
    else:
        location = pool[idx % len(pool)]
    if idx % 23 == 0:
        return OTHER_LOCATIONS[(idx // 17) % len(OTHER_LOCATIONS)]
    return location


def choose_source(profile: dict[str, object], idx: int) -> str:
    options = list(profile["sources"])
    return str(options[idx % len(options)])


def choose_department(syndrome: str, source: str) -> str:
    if syndrome == "医院感染监测":
        return "ICU" if source == "ICU" else "感染科"
    if syndrome == "脑膜炎/脑膜脑炎":
        return "神经内科"
    if syndrome == "败血症/血流感染":
        return "感染科 / 重症医学科"
    if syndrome == "腹泻/肠道感染":
        return "肠道门诊"
    if syndrome == "泌尿生殖道感染":
        return "泌尿科门诊"
    if syndrome == "环境异常事件":
        return "环境采样区"
    return "临床微生物室"


def choose_cluster_status(syndrome: str, source: str, idx: int) -> str:
    if syndrome == "环境异常事件":
        return "聚集性" if idx % 2 == 0 else "待判定"
    if syndrome == "腹泻/肠道感染":
        return "聚集性" if idx % 5 == 0 else "散发"
    if syndrome == "脑膜炎/脑膜脑炎":
        return "待判定" if idx % 4 == 0 else "散发"
    if source == "ICU":
        return "待判定" if idx % 3 == 0 else "散发"
    return "散发"


def tune_syndrome_for_species(species_name: str, source: str, idx: int) -> str:
    if species_name == "Klebsiella pneumoniae":
        return "泌尿生殖道感染" if idx % 7 == 0 else "医院感染监测"
    if species_name == "Staphylococcus aureus":
        return "医院感染监测" if idx % 5 == 0 else "败血症/血流感染"
    if species_name == "Escherichia coli":
        return "腹泻/肠道感染" if idx % 9 == 0 else "泌尿生殖道感染"
    if species_name == "Haemophilus influenzae":
        return "肺炎/呼吸道感染" if idx % 4 == 0 else "脑膜炎/脑膜脑炎"
    if species_name == "Neisseria meningitidis":
        return "不明原因感染" if idx % 11 == 0 else "脑膜炎/脑膜脑炎"
    return str(SPECIES_PROFILES[species_name]["syndrome"])


def tune_specimen_for_syndrome(syndrome: str, default_value: str) -> str:
    mapping = {
        "医院感染监测": "呼吸道",
        "脑膜炎/脑膜脑炎": "脑脊液",
        "肺炎/呼吸道感染": "呼吸道",
        "败血症/血流感染": "血液",
        "腹泻/肠道感染": "粪便/肛拭子",
        "泌尿生殖道感染": "尿液",
        "不明原因感染": default_value,
    }
    return mapping.get(syndrome, default_value)


def tune_sample_source_for_syndrome(syndrome: str, default_value: str) -> str:
    mapping = {
        "医院感染监测": "呼吸道标本",
        "脑膜炎/脑膜脑炎": "脑脊液",
        "肺炎/呼吸道感染": "呼吸道标本",
        "败血症/血流感染": "血培养",
        "腹泻/肠道感染": "粪便",
        "泌尿生殖道感染": "尿液",
    }
    return mapping.get(syndrome, default_value)


def summarize_gene_hits(raw: str, idx: int, *, mge: bool = False) -> str:
    items = [part.strip() for part in str(raw or "").replace("，", "、").split("、") if part.strip()]
    if not items:
        return ""
    if mge:
        if idx % 5 == 0:
            return ""
        keep = 1 if idx % 3 == 0 else min(2, len(items))
        return "、".join(items[:keep])
    keep = 2 if idx % 4 == 0 else min(3, len(items))
    return "、".join(items[:keep])


def rebuild() -> dict[str, int]:
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    rows = conn.execute(
        """
        select rowid as __rowid, sample_key, sample_name, collection_date
          from sample_library
         where library_scope='main'
         order by case
           when collection_date is null or collection_date='' then '9999-99-99'
           else collection_date
         end,
         sample_key
        """
    ).fetchall()
    by_month: dict[str, list[sqlite3.Row]] = defaultdict(list)
    for row in rows:
        collection_date = str(row["collection_date"] or "").strip() or "2026-03-23"
        by_month[collection_date[:7]].append(row)

    updated = 0
    for year_month in sorted(by_month):
        month_rows = by_month[year_month]
        month_num = int(year_month.split("-")[1])
        for idx, row in enumerate(month_rows):
            scenario_key = choose_profile(month_num, idx)
            species_name, location_pool = SCENARIO_LIBRARY[scenario_key]
            profile = SPECIES_PROFILES[species_name]
            location = choose_location(location_pool, idx, month_num)
            st, serotype = profile["st_serotype"][idx % len(profile["st_serotype"])]
            source = choose_source(profile, idx)
            syndrome = tune_syndrome_for_species(species_name, source, idx)
            department = choose_department(syndrome, source)
            cluster_status = choose_cluster_status(syndrome, source, idx)
            province, city, district, detail = location
            collection_date = str(row["collection_date"] or "").strip() or f"{year_month}-23"
            collection_unit = f"{city}疾控采样组"
            submitting_unit = "黄浦区疾控中心" if idx % 7 == 0 else f"{city}示例送检单位"
            case_id = f"DEMO-{year_month.replace('-', '')}-{idx + 1:03d}"
            patient_id = f"PAT-{year_month.replace('-', '')}-{idx + 1:03d}"
            metadata_json = build_metadata(
                case_id=case_id,
                patient_id=patient_id,
                source=source,
                syndrome=syndrome,
                collection_unit=collection_unit,
                submitting_unit=submitting_unit,
                ward_department=department,
                specimen_category=tune_specimen_for_syndrome(syndrome, str(profile["specimen_category"])),
                location=location,
                cluster_status=cluster_status,
                epidemiology_link=f"{year_month} {syndrome} 专题演示链路",
                traditional_result=f"{species_name} 传统检测提示阳性，建议与分型结果联合解读。",
            )
            location_json = json.dumps(
                {"province": province, "city": city, "district": district, "detail": detail},
                ensure_ascii=False,
            )
            description = f"{syndrome}演示样本（{species_name}，{year_month}）"
            note = f"示例样本：用于展示 {syndrome} 在 {year_month} 的监测波动与历史同期对比。"
            conn.execute(
                """
                update sample_library
                   set species_name=?,
                       taxid=?,
                       mlst_species_name=?,
                       mlst_st=?,
                       serotype_result=?,
                       resistance_count=?,
                       virulence_count=?,
                       resistance_gene_hits=?,
                       virulence_gene_hits=?,
                       resistance_mge_hits=?,
                       virulence_mge_hits=?,
                       description=?,
                       country='China',
                       location_json=?,
                       sample_type=?,
                       sequencing_method=?,
                       custom_metadata_json=?,
                       sample_source=?,
                       collection_date=?,
                       host_info=?,
                       note=?,
                       updated_at=?
                 where rowid=?
                """,
                (
                    species_name,
                    profile["taxid"],
                    species_name,
                    st,
                    serotype,
                    profile["resistance_count"],
                    profile["virulence_count"],
                    summarize_gene_hits(str(profile["arg"]), idx),
                    summarize_gene_hits(str(profile["vf"]), idx),
                    summarize_gene_hits(str(profile["arg_mge"]), idx, mge=True),
                    summarize_gene_hits(str(profile["vf_mge"]), idx, mge=True),
                    description,
                    location_json,
                    profile["sample_type"],
                    profile["sequencing_method"],
                    metadata_json,
                    tune_sample_source_for_syndrome(syndrome, str(profile["sample_source"])),
                    collection_date,
                    profile["host_info"],
                    note,
                    NOW_ISO,
                    row["__rowid"],
                ),
            )
            updated += 1
    conn.commit()
    conn.close()
    return {"updated": updated}


if __name__ == "__main__":
    print(json.dumps(rebuild(), ensure_ascii=False))
