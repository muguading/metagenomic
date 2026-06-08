from __future__ import annotations

import csv
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote

from Bio import AlignIO, Phylo, SeqIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from metagenomic_refactor.common import conda_env_path, run_command
from metagenomic_refactor.context import update_runtime_context

GENE_TYPING_MIN_READS = max(0, int(os.environ.get("META_GENE_TYPING_MIN_READS", "20") or 0))
ORTHOHANTAVIRUS_BROAD_MIN_READS = max(0, int(os.environ.get("META_ORTHOHANTAVIRUS_BROAD_MIN_READS", "20") or 0))
ORTHOHANTAVIRUS_BROAD_MIN_COVERAGE_SUM = max(0.0, float(os.environ.get("META_ORTHOHANTAVIRUS_BROAD_MIN_COVERAGE_SUM", "60") or 0.0))
_FASTA_INDEX_CACHE: dict[str, object] = {}
ENTEROVIRUS_REFERENCE_BATCH_SIZE = max(2, int(os.environ.get("META_ENTEROVIRUS_REFERENCE_BATCH_SIZE", "64") or 64))
ENTEROVIRUS_REFERENCE_PLAYOFF_SIZE = max(2, int(os.environ.get("META_ENTEROVIRUS_REFERENCE_PLAYOFF_SIZE", "24") or 24))
ENTEROVIRUS_REFERENCE_DEDUP_SEQ_ID = float(os.environ.get("META_ENTEROVIRUS_REFERENCE_DEDUP_SEQ_ID", "0.95") or 0.95)


def _parse_gff_attributes(field: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for chunk in str(field or "").strip().split(";"):
        item = chunk.strip()
        if not item or "=" not in item:
            continue
        key, value = item.split("=", 1)
        key = key.strip()
        value = unquote(value.strip())
        if key:
            attributes[key] = value
    return attributes


def prepare_snpeff_reference_gff(source_gff: Path, output_gff: Path) -> bool:
    if not source_gff.is_file() or source_gff.stat().st_size == 0:
        return False
    cds_rows: list[tuple[str, str, int, int, str, str, str, dict[str, str]]] = []
    gene_rows: list[tuple[str, str, int, int, str, str, dict[str, str]]] = []
    sequence_regions: list[str] = []
    seen_regions: set[str] = set()
    with source_gff.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith("##sequence-region"):
                if line not in seen_regions:
                    sequence_regions.append(line)
                    seen_regions.add(line)
                continue
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 9:
                continue
            seqid, source, feature_type, start, end, score, strand, phase, attrs = parts
            if feature_type not in {"CDS", "gene"}:
                continue
            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                continue
            attributes = _parse_gff_attributes(attrs)
            normalized_source = source or "RefSeq"
            normalized_score = score or "."
            normalized_strand = strand or "+"
            if feature_type == "CDS":
                cds_rows.append((seqid, normalized_source, start_i, end_i, normalized_score, normalized_strand, phase or "0", attributes))
            elif feature_type == "gene":
                gene_rows.append((seqid, normalized_source, start_i, end_i, normalized_score, normalized_strand, attributes))
    if not cds_rows and not gene_rows:
        shutil.copy2(source_gff, output_gff)
        return True
    output_gff.parent.mkdir(parents=True, exist_ok=True)
    with output_gff.open("w", encoding="utf-8") as handle:
        handle.write("##gff-version 3\n")
        for line in sequence_regions:
            handle.write(f"{line}\n")
        source_rows: list[tuple[str, str, int, int, str, str, str, dict[str, str]]] = cds_rows
        if not source_rows:
            source_rows = [
                (seqid, source, start_i, end_i, score, strand, "0", attributes)
                for seqid, source, start_i, end_i, score, strand, attributes in gene_rows
            ]
        for idx, (seqid, source, start_i, end_i, score, strand, phase, attributes) in enumerate(source_rows, start=1):
            cds_id_raw = (
                attributes.get("ID")
                or attributes.get("protein_id")
                or attributes.get("gene")
                or attributes.get("gene_name")
                or attributes.get("Name")
                or f"cds_{idx}"
            )
            cds_token = re.sub(r"[^A-Za-z0-9_.:-]+", "_", cds_id_raw)
            gene_id = f"gene_{cds_token}"
            transcript_id = f"rna_{cds_token}"
            exon_id = f"exon_{cds_token}"
            gene_name = (
                attributes.get("gene_name")
                or attributes.get("gene")
                or attributes.get("Name")
                or attributes.get("product")
                or cds_token
            )
            product = attributes.get("product") or gene_name
            gene_attrs = f"ID={gene_id};Name={gene_name}"
            transcript_attrs = f"ID={transcript_id};Parent={gene_id};Name={gene_name}"
            exon_attrs = f"ID={exon_id};Parent={transcript_id};Name={gene_name}"
            cds_attrs = f"ID={cds_token};Parent={transcript_id};Name={gene_name};product={product}"
            if attributes.get("protein_id"):
                cds_attrs += f";protein_id={attributes['protein_id']}"
            handle.write(f"{seqid}\t{source}\tgene\t{start_i}\t{end_i}\t{score}\t{strand}\t.\t{gene_attrs}\n")
            handle.write(f"{seqid}\t{source}\tmRNA\t{start_i}\t{end_i}\t{score}\t{strand}\t.\t{transcript_attrs}\n")
            handle.write(f"{seqid}\t{source}\texon\t{start_i}\t{end_i}\t{score}\t{strand}\t.\t{exon_attrs}\n")
            handle.write(f"{seqid}\t{source}\tCDS\t{start_i}\t{end_i}\t{score}\t{strand}\t{phase}\t{cds_attrs}\n")
    return True


INFLUENZA_A_LABELS = {
    "influenza a virus",
    "influenza a",
    "甲型流感病毒",
    "甲流",
}

INFLUENZA_B_LABELS = {
    "influenza b virus",
    "influenza b",
    "乙型流感病毒",
    "乙流",
}

INFLUENZA_UNIFIED_LABELS = {
    "influenza virus",
    "influenza",
    "流感病毒",
    "流感",
} | INFLUENZA_A_LABELS | INFLUENZA_B_LABELS

INFLUENZA_C_LABELS = {
    "influenza c virus",
    "influenza c",
    "丙型流感病毒",
    "丙流",
}

INFLUENZA_D_LABELS = {
    "influenza d virus",
    "influenza d",
    "丁型流感病毒",
    "丁流",
}

INFLUENZA_TYPE_TO_LABEL = {
    "A": "Influenza A virus",
    "B": "Influenza B virus",
    "C": "Influenza C virus",
    "D": "Influenza D virus",
}

SARS_COV_2_LABELS = {
    "sars-cov-2",
    "sars cov 2",
    "severe acute respiratory syndrome coronavirus 2",
    "2019-ncov",
    "covid-19 virus",
    "covid19 virus",
    "新型冠状病毒",
    "新冠病毒",
    "新冠",
}

MONKEYPOX_LABELS = {
    "monkeypox virus",
    "mpox virus",
    "human monkeypox virus",
    "猴痘病毒",
    "猴痘",
    "mpox",
    "hmpxv",
}

RSV_A_LABELS = {
    "respiratory syncytial virus a",
    "human respiratory syncytial virus a",
    "orthopneumovirus hominis a",
    "呼吸道合胞病毒a",
    "呼吸道合胞病毒 a",
    "rsv a",
    "rsv-a",
}

RSV_B_LABELS = {
    "respiratory syncytial virus b",
    "human respiratory syncytial virus b",
    "orthopneumovirus hominis b",
    "呼吸道合胞病毒b",
    "呼吸道合胞病毒 b",
    "rsv b",
    "rsv-b",
}

RSV_LABELS = {
    "respiratory syncytial virus",
    "human respiratory syncytial virus",
    "orthopneumovirus hominis",
    "呼吸道合胞病毒",
    "rsv",
    "h/rsv",
} | RSV_A_LABELS | RSV_B_LABELS

DENV_1_LABELS = {
    "dengue virus 1",
    "dengue virus type 1",
    "denv1",
    "denv-1",
    "登革热病毒1型",
    "登革热1型",
}

DENV_2_LABELS = {
    "dengue virus 2",
    "dengue virus type 2",
    "denv2",
    "denv-2",
    "登革热病毒2型",
    "登革热2型",
}

DENV_3_LABELS = {
    "dengue virus 3",
    "dengue virus type 3",
    "denv3",
    "denv-3",
    "登革热病毒3型",
    "登革热3型",
}

DENV_4_LABELS = {
    "dengue virus 4",
    "dengue virus type 4",
    "denv4",
    "denv-4",
    "登革热病毒4型",
    "登革热4型",
}

DENV_LABELS = {
    "dengue virus",
    "dengue",
    "登革热病毒",
    "登革热",
    "denv",
} | DENV_1_LABELS | DENV_2_LABELS | DENV_3_LABELS | DENV_4_LABELS

CHIKV_LABELS = {
    "chikungunya virus",
    "chikungunya",
    "chikv",
    "基孔肯雅病毒",
    "基孔肯雅",
}

HPIV_1_LABELS = {
    "human parainfluenza virus 1",
    "parainfluenza virus 1",
    "hpiv1",
    "hpiv-1",
    "副流感病毒1型",
    "副流感1型",
}

HPIV_2_LABELS = {
    "human parainfluenza virus 2",
    "parainfluenza virus 2",
    "hpiv2",
    "hpiv-2",
    "副流感病毒2型",
    "副流感2型",
}

HPIV_3_LABELS = {
    "human parainfluenza virus 3",
    "parainfluenza virus 3",
    "hpiv3",
    "hpiv-3",
    "副流感病毒3型",
    "副流感3型",
}

HPIV_4A_LABELS = {
    "human parainfluenza virus 4a",
    "parainfluenza virus 4a",
    "hpiv4a",
    "hpiv-4a",
    "副流感病毒4a型",
    "副流感4a型",
}

HPIV_4B_LABELS = {
    "human parainfluenza virus 4b",
    "parainfluenza virus 4b",
    "hpiv4b",
    "hpiv-4b",
    "副流感病毒4b型",
    "副流感4b型",
}

HPIV_LABELS = {
    "human parainfluenza virus",
    "parainfluenza virus",
    "副流感病毒",
    "hpiv",
} | HPIV_1_LABELS | HPIV_2_LABELS | HPIV_3_LABELS | HPIV_4A_LABELS | HPIV_4B_LABELS

HADV_LABELS = {
    "human adenovirus",
    "human adenovirus virus",
    "human mastadenovirus",
    "mastadenovirus hominis",
    "adenovirus",
    "hadv",
    "人腺病毒",
    "腺病毒",
}

NOROVIRUS_LABELS = {
    "norovirus",
    "human norovirus",
    "norwalk virus",
    "诺如病毒",
    "诺瓦克病毒",
    "noro",
}

HEPATOVIRUS_LABELS = {
    "hepatovirus",
    "hepatitis a virus",
    "hepatitis b virus",
    "hepatitis c virus",
    "hepatitis d virus",
    "hepatitis e virus",
    "hepatovirus ahepa",
    "hav",
    "hbv",
    "hcv",
    "hdv",
    "hev",
    "甲肝病毒",
    "乙肝病毒",
    "丙肝病毒",
    "丁肝病毒",
    "戊肝病毒",
    "甲型肝炎病毒",
    "乙型肝炎病毒",
    "丙型肝炎病毒",
    "丁型肝炎病毒",
    "戊型肝炎病毒",
    "甲型肝炎",
    "乙型肝炎",
    "丙型肝炎",
    "丁型肝炎",
    "戊型肝炎",
}

RHINOVIRUS_LABELS = {
    "human rhinovirus",
    "rhinovirus",
    "鼻病毒",
    "hrv",
    "rv",
    "human rhinovirus a",
    "human rhinovirus b",
    "human rhinovirus c",
    "rhinovirus a",
    "rhinovirus b",
    "rhinovirus c",
}

ENTEROVIRUS_LABELS = {
    "enterovirus",
    "human enterovirus",
    "肠道病毒",
    "人肠道病毒",
}

ENTEROVIRUS_A_LABELS = {
    "enterovirus a",
    "human enterovirus a",
    "ev-a",
    "eva",
    "肠道病毒a",
    "肠道病毒a组",
}

ENTEROVIRUS_B_LABELS = {
    "enterovirus b",
    "human enterovirus b",
    "ev-b",
    "evb",
    "肠道病毒b",
    "肠道病毒b组",
}

ENTEROVIRUS_C_LABELS = {
    "enterovirus c",
    "human enterovirus c",
    "ev-c",
    "evc",
    "肠道病毒c",
    "肠道病毒c组",
}

ENTEROVIRUS_D_LABELS = {
    "enterovirus d",
    "human enterovirus d",
    "ev-d",
    "evd",
    "肠道病毒d",
    "肠道病毒d组",
}

BANDAVIRUS_LABELS = {
    "bandavirus",
    "bandavirus dabieense",
    "bandavirus heartlandense",
    "bandavirus guertuense",
    "bandavirus bhanjanagarense",
    "bandavirus amblyommae",
    "bandavirus kismaayoense",
    "bandavirus razdanense",
    "bandavirus albatrossense",
    "bandavirus zwieselense",
    "severe fever with thrombocytopenia syndrome virus",
    "sftsv",
    "heartland virus",
    "guertu virus",
    "bhanja virus",
    "lone star virus",
    "kismaayo virus",
    "razdan virus",
    "hunter island virus",
    "zwiesel bat bandavirus",
    "班达病毒",
    "发热伴血小板减少综合征病毒",
}

ORTHOHANTAVIRUS_LABELS = {
    "orthohantavirus",
    "hantavirus",
    "orthohantavirus",
    "汉坦病毒",
    "汉他病毒",
}

SFTSV_LABELS = {
    "bandavirus dabieense",
    "severe fever with thrombocytopenia syndrome virus",
    "sftsv",
    "发热伴血小板减少综合征病毒",
}

ASTROVIRIDAE_LABELS = {
    "astroviridae",
    "astrovirus",
    "human astrovirus",
    "avian astrovirus",
    "mamastrovirus",
    "avastrovirus",
    "星状病毒",
}

MAMASTROVIRUS_LABELS = {
    "mamastrovirus",
    "human astrovirus",
    "mamastrovirus hominis",
    "mamastrovirus 1",
    "mamastrovirus 2",
    "mamastrovirus 3",
    "mamastrovirus 5",
    "mamastrovirus 6",
    "mamastrovirus 8",
    "mamastrovirus 9",
}

AVASTROVIRUS_LABELS = {
    "avastrovirus",
    "avian astrovirus",
    "avastrovirus 1",
    "avastrovirus 2",
    "avastrovirus 3",
}

ROTAVIRUS_A_LABELS = {
    "rotavirus a",
    "human rotavirus a",
    "a组轮状病毒",
    "轮状病毒a",
    "轮状病毒a组",
    "rva",
}

ROTAVIRUS_B_LABELS = {
    "rotavirus b",
    "human rotavirus b",
    "b组轮状病毒",
    "轮状病毒b",
    "轮状病毒b组",
    "rvb",
}

ROTAVIRUS_C_LABELS = {
    "rotavirus c",
    "human rotavirus c",
    "c组轮状病毒",
    "轮状病毒c",
    "轮状病毒c组",
    "rvc",
}

ROTAVIRUS_LABELS = {
    "rotavirus",
    "human rotavirus",
    "轮状病毒",
} | ROTAVIRUS_A_LABELS | ROTAVIRUS_B_LABELS | ROTAVIRUS_C_LABELS

SEASONAL_HCOV_229E_LABELS = {
    "human coronavirus 229e",
    "human coronavirus hcov-229e",
    "hcov-229e",
    "hcov 229e",
    "229e",
    "人冠状病毒229e",
}

SEASONAL_HCOV_NL63_LABELS = {
    "human coronavirus nl63",
    "human coronavirus hcov-nl63",
    "hcov-nl63",
    "hcov nl63",
    "nl63",
    "人冠状病毒nl63",
}

SEASONAL_HCOV_OC43_LABELS = {
    "human coronavirus oc43",
    "human coronavirus hcov-oc43",
    "hcov-oc43",
    "hcov oc43",
    "oc43",
    "人冠状病毒oc43",
}

SEASONAL_HCOV_HKU1_LABELS = {
    "human coronavirus hku1",
    "human coronavirus hcov-hku1",
    "hcov-hku1",
    "hcov hku1",
    "hku1",
    "人冠状病毒hku1",
}

SEASONAL_HCOV_LABELS = {
    "human coronavirus",
    "seasonal human coronavirus",
    "seasonal coronavirus",
    "hcov",
    "human cov",
    "人冠状病毒",
    "季节性冠状病毒",
} | SEASONAL_HCOV_229E_LABELS | SEASONAL_HCOV_NL63_LABELS | SEASONAL_HCOV_OC43_LABELS | SEASONAL_HCOV_HKU1_LABELS

SEASONAL_HCOV_TYPE_LABELS = {
    "HCoV-229E": SEASONAL_HCOV_229E_LABELS,
    "HCoV-NL63": SEASONAL_HCOV_NL63_LABELS,
    "HCoV-OC43": SEASONAL_HCOV_OC43_LABELS,
    "HCoV-HKU1": SEASONAL_HCOV_HKU1_LABELS,
}

HADV_MANUAL_PHF_COMBOS = {
    ("HAdV-D37", "HAdV-D22", "HAdV-D8"): "HAdV-D53",
    ("HAdV-B11", "HAdV-B11", "HAdV-B14"): "HAdV-B55",
    ("HAdV-D30", "HAdV-D30", "HAdV-D29"): "HAdV-D63",
    ("HAdV-B11", "HAdV-B14", "HAdV-B35"): "HAdV-B68",
    ("HAdV-D67", "HAdV-D45", "HAdV-D27"): "HAdV-D73",
    ("HAdV-B21", "HAdV-B21", "HAdV-B16"): "HAdV-B76",
    ("HAdV-B35", "HAdV-B34", "HAdV-B7"): "HAdV-B77",
    ("HAdV-B11", "HAdV-B11", "HAdV-B7"): "HAdV-B78",
    ("HAdV-B11", "HAdV-B34", "HAdV-B11"): "HAdV-B79",
    ("HAdV-D22", "HAdV-D28", "HAdV-D22"): "HAdV-D80",
    ("HAdV-D65", "HAdV-D48", "HAdV-D60"): "HAdV-D81",
    ("HAdV-D56", "HAdV-D15", "HAdV-D37"): "HAdV-D82",
    ("HAdV-D37", "HAdV-D19", "HAdV-D8"): "HAdV-D85",
    ("HAdV-D9", "HAdV-D25", "HAdV-D25"): "HAdV-D86",
    ("HAdV-D9", "HAdV-D15", "HAdV-D25"): "HAdV-D87",
    ("HAdV-D37", "HAdV-D37", "HAdV-D17"): "HAdV-D91",
    ("HAdV-D28", "HAdV-D37", "HAdV-D38"): "HAdV-D93",
    ("HAdV-D33", "HAdV-D15", "HAdV-D9"): "HAdV-D94",
    ("HAdV-D23", "HAdV-D32", "HAdV-D62"): "HAdV-D96",
    ("HAdV-D67", "HAdV-D28", "HAdV-D60"): "HAdV-D97",
    ("HAdV-D9", "HAdV-D46", "HAdV-D39"): "HAdV-D99",
    ("HAdV-C1", "HAdV-C1", "HAdV-C2"): "HAdV-C104",
    ("HAdV-B11", "HAdV-B11", "HAdV-B35"): "HAdV-B106",
    ("HAdV-C1", "HAdV-C2", "HAdV-C2"): "HAdV-C108",
    ("HAdV-D22", "HAdV-D19", "HAdV-D9"): "HAdV-D109",
    ("HAdV-D37", "HAdV-D9", "HAdV-D9"): "HAdV-D111",
    ("HAdV-D112", "HAdV-D112", "HAdV-D67"): "HAdV-D112",
    ("HAdV-D20", "HAdV-D42", "HAdV-D42"): "HAdV-D113",
}

AA3_TO_1 = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Gln": "Q",
    "Glu": "E",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
    "Ter": "*",
    "Stop": "*",
}


def _normalize_species_label(value: str) -> str:
    return str(value or "").strip().lower()


def _is_influenza_a(species: str) -> bool:
    normalized = _normalize_species_label(species)
    return normalized in INFLUENZA_A_LABELS


def _is_influenza_b(species: str) -> bool:
    normalized = _normalize_species_label(species)
    return normalized in INFLUENZA_B_LABELS


def _is_influenza(species: str) -> bool:
    normalized = _normalize_species_label(species)
    return normalized in INFLUENZA_UNIFIED_LABELS | INFLUENZA_C_LABELS | INFLUENZA_D_LABELS


def _is_sars_cov_2(species: str) -> bool:
    normalized = _normalize_species_label(species)
    return normalized in SARS_COV_2_LABELS or "sars-cov-2" in normalized or "新冠" in normalized


def _is_monkeypox(species: str) -> bool:
    normalized = _normalize_species_label(species)
    return normalized in MONKEYPOX_LABELS or "monkeypox" in normalized or "mpox" in normalized or "猴痘" in normalized


def _is_rsv(species: str) -> bool:
    normalized = _normalize_species_label(species)
    return normalized in RSV_LABELS or "respiratory syncytial" in normalized or normalized == "rsv" or "合胞病毒" in normalized


def _is_denv(species: str) -> bool:
    normalized = _normalize_species_label(species)
    return normalized in DENV_LABELS or "dengue" in normalized or "登革热" in normalized or normalized == "denv"


def _is_zika(species: str) -> bool:
    normalized = _normalize_species_label(species)
    return (
        normalized in {"zika virus", "zika", "zikv", "zikav", "寨卡病毒", "寨卡"}
        or "zika" in normalized
        or "寨卡" in normalized
    )


def _is_chikv(species: str) -> bool:
    normalized = _normalize_species_label(species)
    return normalized in CHIKV_LABELS or "chikungunya" in normalized or "基孔肯雅" in normalized


def _is_hpiv(species: str) -> bool:
    normalized = _normalize_species_label(species)
    return normalized in HPIV_LABELS or "parainfluenza" in normalized or "副流感" in normalized or normalized == "hpiv"


def _is_hadv(species: str) -> bool:
    normalized = _normalize_species_label(species)
    return normalized in HADV_LABELS or normalized.startswith("hadv") or "adenovirus" in normalized or "mastadenovirus" in normalized or "腺病毒" in normalized


def _is_norovirus(species: str) -> bool:
    normalized = _normalize_species_label(species)
    return normalized in NOROVIRUS_LABELS or "norovirus" in normalized or "norwalk" in normalized or "诺如" in normalized


def _is_hepatovirus(species: str) -> bool:
    normalized = _normalize_species_label(species)
    return (
        normalized in HEPATOVIRUS_LABELS
        or "hepatovirus" in normalized
        or re.search(r"\bhepatitis\s+[abcde]\s+virus\b", normalized) is not None
        or normalized in {"hav", "hbv", "hcv", "hdv", "hev"}
        or any(token in normalized for token in ("甲肝", "乙肝", "丙肝", "丁肝", "戊肝", "肝炎病毒"))
    )


def _is_hiv(species: str) -> bool:
    normalized = _normalize_species_label(species)
    return (
        normalized in {"hiv", "hiv-1", "hiv1", "hiv-2", "hiv2", "human immunodeficiency virus", "human immunodeficiency virus 1", "human immunodeficiency virus 2"}
        or "human immunodeficiency virus" in normalized
        or re.search(r"\bhiv(?:[-\s]?[12])?\b", normalized) is not None
        or "艾滋" in normalized
    )


def _is_rhinovirus(species: str) -> bool:
    normalized = _normalize_species_label(species)
    return normalized in RHINOVIRUS_LABELS or "rhinovirus" in normalized or "鼻病毒" in normalized


def _is_enterovirus(species: str) -> bool:
    normalized = _normalize_species_label(species)
    if not normalized or _is_rhinovirus(species):
        return False
    return (
        normalized in ENTEROVIRUS_LABELS
        or normalized in ENTEROVIRUS_A_LABELS
        or normalized in ENTEROVIRUS_B_LABELS
        or normalized in ENTEROVIRUS_C_LABELS
        or normalized in ENTEROVIRUS_D_LABELS
        or "enterovirus" in normalized
        or "肠道病毒" in normalized
    )


def _is_bandavirus(species: str) -> bool:
    normalized = _normalize_species_label(species)
    return (
        normalized in BANDAVIRUS_LABELS
        or "bandavirus" in normalized
        or "sftsv" in normalized
        or "thrombocytopenia syndrome virus" in normalized
        or "班达病毒" in normalized
        or "血小板减少综合征病毒" in normalized
    )


def _is_orthohantavirus(species: str) -> bool:
    normalized = _normalize_species_label(species)
    if not normalized or _is_bandavirus(species):
        return False
    return (
        normalized in ORTHOHANTAVIRUS_LABELS
        or "orthohantavirus" in normalized
        or "orthohantavirus" in normalized
        or "hantavirus" in normalized
        or "汉坦" in normalized
        or "汉他" in normalized
    )


def _is_orthoebolavirus(species: str) -> bool:
    normalized = _normalize_species_label(species)
    return (
        normalized in {
            "orthoebolavirus",
            "ebola virus",
            "ebolavirus",
            "ebov",
            "sudan virus",
            "sudv",
            "bundibugyo virus",
            "bdbv",
            "reston virus",
            "restv",
            "tai forest virus",
            "taï forest virus",
            "tafv",
            "bombali virus",
            "bomv",
            "埃博拉病毒",
            "正埃博拉病毒",
        }
        or "orthoebolavirus" in normalized
        or "ebolavirus" in normalized
        or "ebola" in normalized
        or "埃博拉" in normalized
    )


def _is_astroviridae(species: str) -> bool:
    normalized = _normalize_species_label(species)
    return (
        normalized in ASTROVIRIDAE_LABELS
        or normalized in MAMASTROVIRUS_LABELS
        or normalized in AVASTROVIRUS_LABELS
        or "astrovirus" in normalized
        or "mamastrovirus" in normalized
        or "avastrovirus" in normalized
        or "星状病毒" in normalized
    )


def _sample_skip_flag_path(pre: str) -> Path:
    return Path(f"{pre}.skip_remaining.txt")


def _clear_sample_skip_flag(pre: str) -> None:
    path = _sample_skip_flag_path(pre)
    if path.is_file():
        path.unlink()


def _write_sample_skip_flag(pre: str, note: str) -> None:
    _sample_skip_flag_path(pre).write_text(str(note or "").strip() + "\n", encoding="utf-8")


def _is_low_support_gene_typing_hit(hit: dict[str, object] | None) -> bool:
    if not isinstance(hit, dict):
        return False
    if str(hit.get("method") or "").strip() != "read_coverage":
        return False
    try:
        support_reads = float(hit.get("num_reads") or 0.0)
    except (TypeError, ValueError):
        support_reads = 0.0
    return support_reads > 0.0 and support_reads < float(GENE_TYPING_MIN_READS)


def _low_support_gene_typing_note(gene_name: str, hit: dict[str, object] | None) -> str:
    if not isinstance(hit, dict):
        return f"{gene_name} 分型命中支持 reads 过少，已跳过当前样本。"
    try:
        support_reads = int(float(hit.get("num_reads") or 0.0))
    except (TypeError, ValueError):
        support_reads = 0
    matched_type = str(hit.get("type") or "").strip()
    return (
        f"{gene_name} 分型命中支持 reads 过少"
        f"（{support_reads} < {GENE_TYPING_MIN_READS}）"
        f"{f'，当前命中为 {matched_type}' if matched_type else ''}，已跳过当前样本。"
    )


def _is_low_support_orthohantavirus_broad_result(result: dict[str, object] | None) -> bool:
    if not isinstance(result, dict):
        return False
    if str(result.get("method") or "").strip() != "read_coverage":
        return False
    try:
        segment_count = int(float(result.get("segment_count") or 0))
    except (TypeError, ValueError):
        segment_count = 0
    try:
        coverage_sum = float(result.get("coverage_sum") or 0.0)
    except (TypeError, ValueError):
        coverage_sum = 0.0
    try:
        support_reads = float(result.get("num_reads_sum") or 0.0)
    except (TypeError, ValueError):
        support_reads = 0.0
    if segment_count <= 0:
        return True
    if coverage_sum < ORTHOHANTAVIRUS_BROAD_MIN_COVERAGE_SUM:
        return True
    return support_reads > 0.0 and support_reads < float(ORTHOHANTAVIRUS_BROAD_MIN_READS)


def _low_support_orthohantavirus_broad_note(result: dict[str, object] | None, species: str = "") -> str:
    if not isinstance(result, dict):
        return f"汉坦病毒 broad 分型支持不足，species={species or '-'}，判定为非汉坦病毒或证据不足，不继续后续 L/M/S 细分。"
    label = str(result.get("label") or "").strip().upper()
    try:
        segment_count = int(float(result.get("segment_count") or 0))
    except (TypeError, ValueError):
        segment_count = 0
    try:
        coverage_sum = float(result.get("coverage_sum") or 0.0)
    except (TypeError, ValueError):
        coverage_sum = 0.0
    try:
        support_reads = int(float(result.get("num_reads_sum") or 0.0))
    except (TypeError, ValueError):
        support_reads = 0
    parts = [
        "汉坦病毒 broad 分型支持不足",
        f"segments={segment_count}",
        f"coverage_sum={coverage_sum:.2f}",
        f"reads={support_reads}",
    ]
    if label:
        parts.append(f"最高命中={label}")
    if species:
        parts.append(f"species={species}")
    parts.append("判定为非汉坦病毒或证据不足，不继续后续 L/M/S 细分")
    return "，".join(parts) + "。"


def _is_rotavirus(species: str) -> bool:
    normalized = _normalize_species_label(species)
    return normalized in ROTAVIRUS_LABELS or "rotavirus" in normalized or "轮状病毒" in normalized


def _is_seasonal_hcov(species: str) -> bool:
    normalized = _normalize_species_label(species)
    if not normalized or _is_sars_cov_2(species):
        return False
    return (
        normalized in SEASONAL_HCOV_LABELS
        or "seasonal coronavirus" in normalized
        or ("coronavirus" in normalized and any(token in normalized for token in ["229e", "nl63", "oc43", "hku1"]))
        or ("冠状病毒" in normalized and any(token in normalized for token in ["229e", "nl63", "oc43", "hku1"]))
    )


def _infer_seasonal_hcov_type_from_label(label: str) -> str:
    normalized = _normalize_species_label(label)
    for type_label, aliases in SEASONAL_HCOV_TYPE_LABELS.items():
        if normalized in aliases:
            return type_label
    if "229e" in normalized:
        return "HCoV-229E"
    if "nl63" in normalized:
        return "HCoV-NL63"
    if "oc43" in normalized:
        return "HCoV-OC43"
    if "hku1" in normalized:
        return "HCoV-HKU1"
    return ""


def _seasonal_hcov_species_label(type_label: str = "") -> str:
    normalized = str(type_label or "").strip()
    mapping = {
        "HCoV-229E": "Human coronavirus 229E",
        "HCoV-NL63": "Human coronavirus NL63",
        "HCoV-OC43": "Human coronavirus OC43",
        "HCoV-HKU1": "Human coronavirus HKU1",
    }
    return mapping.get(normalized, "Human coronavirus")


def _infer_rsv_type_from_label(label: str) -> str:
    normalized = _normalize_species_label(label)
    if normalized in RSV_A_LABELS or "respiratory syncytial virus a" in normalized or "rsv-a" in normalized or normalized.endswith(" rsv a"):
        return "A"
    if normalized in RSV_B_LABELS or "respiratory syncytial virus b" in normalized or "rsv-b" in normalized or normalized.endswith(" rsv b"):
        return "B"
    return "-"


def _infer_denv_type_from_label(label: str) -> str:
    normalized = _normalize_species_label(label)
    if normalized in DENV_1_LABELS or "dengue virus 1" in normalized or "denv-1" in normalized:
        return "1"
    if normalized in DENV_2_LABELS or "dengue virus 2" in normalized or "denv-2" in normalized:
        return "2"
    if normalized in DENV_3_LABELS or "dengue virus 3" in normalized or "denv-3" in normalized:
        return "3"
    if normalized in DENV_4_LABELS or "dengue virus 4" in normalized or "denv-4" in normalized:
        return "4"
    match = re.search(r"\bdenv[\s\-_]*([1-4])\b", normalized, re.IGNORECASE)
    if match:
        return match.group(1)
    match = re.search(r"登革热(?:病毒)?\s*([1-4])\s*型", normalized)
    if match:
        return match.group(1)
    return "-"


def _infer_hpiv_type_from_label(label: str) -> str:
    normalized = _normalize_species_label(label)
    if normalized in HPIV_1_LABELS or "parainfluenza virus 1" in normalized or "hpiv-1" in normalized:
        return "1"
    if normalized in HPIV_2_LABELS or "parainfluenza virus 2" in normalized or "hpiv-2" in normalized:
        return "2"
    if normalized in HPIV_3_LABELS or "parainfluenza virus 3" in normalized or "hpiv-3" in normalized:
        return "3"
    if normalized in HPIV_4A_LABELS or "parainfluenza virus 4a" in normalized or "hpiv-4a" in normalized:
        return "4A"
    if normalized in HPIV_4B_LABELS or "parainfluenza virus 4b" in normalized or "hpiv-4b" in normalized:
        return "4B"
    match = re.search(r"\bhpiv[\s\-_]*([1-3])\b", normalized, re.IGNORECASE)
    if match:
        return match.group(1)
    match = re.search(r"\bhpiv[\s\-_]*4[\s\-_]*([ab])\b", normalized, re.IGNORECASE)
    if match:
        return f"4{match.group(1).upper()}"
    match = re.search(r"副流感(?:病毒)?\s*([1-3])\s*型", normalized)
    if match:
        return match.group(1)
    match = re.search(r"副流感(?:病毒)?\s*4\s*([abAB])\s*型", normalized)
    if match:
        return f"4{match.group(1).upper()}"
    return "-"


def _extract_rhinovirus_species_group(value: str) -> str:
    text = str(value or "").strip().upper()
    if not text:
        return ""
    if text in {"A", "B", "C"}:
        return text
    match = re.search(r"RV[-_ ]?([ABC])", text)
    if match:
        return match.group(1)
    match = re.search(r"RHINOVIRUS[ _-]?([ABC])", text)
    if match:
        return match.group(1)
    if "HUMAN RHINOVIRUS A" in text or "RHINOVIRUS A" in text:
        return "A"
    if "HUMAN RHINOVIRUS B" in text or "RHINOVIRUS B" in text:
        return "B"
    if "HUMAN RHINOVIRUS C" in text or "RHINOVIRUS C" in text:
        return "C"
    return ""


def _normalize_rhinovirus_type(value: str) -> str:
    text = str(value or "").strip().upper()
    if not text:
        return ""
    direct = re.search(r"RV[-_ ]?([ABC])[-_ ]?([0-9]+)", text)
    if direct:
        return f"RV-{direct.group(1)}{direct.group(2)}"
    group = _extract_rhinovirus_species_group(text)
    digits = re.findall(r"\d+", text)
    if group and digits:
        return f"RV-{group}{digits[0]}"
    return text.replace(" ", "")


def _normalize_rhinovirus_type_from_row(row: dict[str, object]) -> str:
    normalized = _normalize_rhinovirus_type(str(row.get("normalized_type") or ""))
    if normalized:
        return normalized
    group = str(row.get("species_group") or "").strip().upper()
    type_label = str(row.get("type_label") or "").strip()
    if group and type_label:
        return _normalize_rhinovirus_type(f"RV-{group}{type_label}")
    return _normalize_rhinovirus_type(type_label)


def _rhinovirus_species_label(group: str) -> str:
    normalized = str(group or "").strip().upper()
    if normalized in {"A", "B", "C"}:
        return f"Rhinovirus {normalized}"
    return "Human rhinovirus"


def _infer_rotavirus_group_from_label(label: str) -> str:
    normalized = _normalize_species_label(label)
    if normalized in ROTAVIRUS_A_LABELS or "rotavirus a" in normalized or "轮状病毒a" in normalized or "a组轮状病毒" in normalized:
        return "A"
    if normalized in ROTAVIRUS_B_LABELS or "rotavirus b" in normalized or "轮状病毒b" in normalized or "b组轮状病毒" in normalized:
        return "B"
    if normalized in ROTAVIRUS_C_LABELS or "rotavirus c" in normalized or "轮状病毒c" in normalized or "c组轮状病毒" in normalized:
        return "C"
    return ""


def _rotavirus_species_label(group: str) -> str:
    normalized = str(group or "").strip().upper()
    if normalized in {"A", "B", "C"}:
        return f"Human rotavirus {normalized}"
    return "Rotavirus"


def _parse_rotavirus_combo(combo: str) -> tuple[str, str]:
    text = str(combo or "").strip()
    if not text:
        return "", ""
    match = re.search(r"(G\d+)\s*(P\[[^\]]+\])", text, flags=re.IGNORECASE)
    if match:
        return match.group(1).upper(), match.group(2).upper()
    g_match = re.search(r"(G\d+)", text, flags=re.IGNORECASE)
    p_match = re.search(r"(P\[[^\]]+\])", text, flags=re.IGNORECASE)
    return (
        g_match.group(1).upper() if g_match else "",
        p_match.group(1).upper() if p_match else "",
    )


def _infer_influenza_type_from_label(label: str) -> str:
    text = _normalize_species_label(label)
    if text in INFLUENZA_A_LABELS or "influenza a" in text or "甲流" in text:
        return "Influenza A virus"
    if text in INFLUENZA_B_LABELS or "influenza b" in text or "乙流" in text:
        return "Influenza B virus"
    if text in INFLUENZA_C_LABELS or "influenza c" in text or "丙流" in text:
        return "Influenza C virus"
    if text in INFLUENZA_D_LABELS or "influenza d" in text or "丁流" in text:
        return "Influenza D virus"
    return "-"


def _infer_norovirus_dual_type_from_label(label: str) -> str:
    text = str(label or "").strip()
    matched = re.search(r"(G[IVX]+\.[Pp][A-Za-z0-9 ._-]+_[GIVX]+\.[A-Za-z0-9 ._-]+)", text, flags=re.IGNORECASE)
    if matched:
        return matched.group(1).strip()
    return ""


def _infer_hepatovirus_broad_type_from_label(label: str) -> str:
    normalized = _normalize_species_label(label)
    if not normalized:
        return ""
    if "hepatitis a virus" in normalized or normalized == "hav" or "甲肝" in normalized:
        return "HAV"
    if "hepatitis b virus" in normalized or normalized == "hbv" or "乙肝" in normalized:
        return "HBV"
    if "hepatitis c virus" in normalized or normalized == "hcv" or "丙肝" in normalized:
        return "HCV"
    if "hepatitis d virus" in normalized or normalized == "hdv" or "丁肝" in normalized:
        return "HDV"
    if "hepatitis e virus" in normalized or normalized == "hev" or "戊肝" in normalized:
        return "HEV"
    matched = re.search(r"\bhepv[-\s]?([a-z0-9]+)\b", normalized, flags=re.IGNORECASE)
    if matched:
        suffix = str(matched.group(1) or "").strip().upper()
        return f"HepV-{suffix}" if suffix else ""
    if "phopivirus" in normalized or normalized == "phv":
        return "PhV"
    return ""


def _looks_like_typing_reference(path_value: str | Path) -> bool:
    text = _normalize_species_label(str(path_value or ""))
    if not text:
        return False
    return any(flag in text for flag in [
        "ha_subtypes",
        "na_subtypes",
        "type_refs",
        "blastdb",
        "insaflu",
    ])


def _project_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _build_vadr_env(project_root: Path, model_dir: Path) -> dict[str, str] | None:
    vadr_root = (project_root / "soft").resolve()
    vadr_scripts_dir = (vadr_root / "vadr").resolve()
    infernal_bin_dir = (vadr_root / "infernal" / "binaries").resolve()
    bio_easel_dir = (vadr_root / "Bio-Easel").resolve()
    alt_bio_easel_dir = (vadr_root / "Bio-Easel-ncov").resolve()
    if (alt_bio_easel_dir / "blib" / "lib").exists() and (alt_bio_easel_dir / "blib" / "arch").exists():
        bio_easel_dir = alt_bio_easel_dir
    sequip_dir = (vadr_root / "sequip").resolve()
    blast_bin_dir = (vadr_root / "ncbi-blast" / "bin").resolve()
    fasta_bin_dir = (vadr_root / "fasta" / "bin").resolve()
    minimap2_dir = (vadr_root / "minimap2").resolve()
    required_paths = [
        vadr_scripts_dir / "v-annotate.pl",
        vadr_scripts_dir / "miniscripts" / "annotate-tbl2gff.pl",
        infernal_bin_dir / "cmalign",
        bio_easel_dir / "blib" / "lib",
        bio_easel_dir / "blib" / "arch",
        sequip_dir,
        blast_bin_dir / "blastn",
        fasta_bin_dir / "fasta36",
        minimap2_dir / "minimap2",
        model_dir,
    ]
    if not all(path.exists() for path in required_paths):
        return None
    inherited_path = str(os.environ.get("PATH") or "").strip()
    inherited_perl5lib = str(os.environ.get("PERL5LIB") or "").strip()
    return {
        **os.environ,
        "VADRINSTALLDIR": str(vadr_root),
        "VADRSCRIPTSDIR": str(vadr_scripts_dir),
        "VADRCONFIGFILE": str((vadr_scripts_dir / "vadr.config").resolve()),
        "VADRMODELDIR": str(model_dir.resolve()),
        "VADRINFERNALDIR": str(infernal_bin_dir),
        "VADREASELDIR": str(infernal_bin_dir),
        "VADRHMMERDIR": str(infernal_bin_dir),
        "VADRBIOEASELDIR": str(bio_easel_dir),
        "VADRSEQUIPDIR": str(sequip_dir),
        "VADRBLASTDIR": str(blast_bin_dir),
        "VADRFASTADIR": str(fasta_bin_dir),
        "VADRMINIMAP2DIR": str(minimap2_dir),
        "PERL5LIB": os.pathsep.join(
            [
                str(vadr_scripts_dir),
                str(sequip_dir),
                str((bio_easel_dir / "blib" / "lib").resolve()),
                str((bio_easel_dir / "blib" / "arch").resolve()),
            ] + ([inherited_perl5lib] if inherited_perl5lib else [])
        ),
        "PATH": os.pathsep.join(
            [
                str(vadr_scripts_dir),
                str(blast_bin_dir),
                str(fasta_bin_dir),
                str(infernal_bin_dir),
                str(minimap2_dir),
            ] + ([inherited_path] if inherited_path else [])
        ),
    }


def _resolve_vadr_perl_bin() -> str:
    for candidate in ["/usr/bin/perl", "perl"]:
        try:
            completed = subprocess.run(
                [candidate, "-MInline", "-e", "print qq(ok\\n)"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )
            if completed.returncode == 0:
                return candidate
        except Exception:
            continue
    return "perl"


def _resolve_working_perl_with_module(env: dict[str, str] | None, module_name: str, candidates: list[str]) -> str:
    for candidate in candidates:
        if not candidate:
            continue
        try:
            completed = subprocess.run(
                [candidate, f"-M{module_name}", "-e", "print qq(ok\\n)"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
                env=env,
            )
            if completed.returncode == 0:
                return candidate
        except Exception:
            continue
    return candidates[0] if candidates else "perl"


def _resolve_flu_vadr_model_dir(project_root: Path) -> Path | None:
    primary = (project_root / "soft" / "vadr-models-flu").resolve()
    if primary.is_dir():
        return primary
    fallback = (project_root / "soft" / "vadr-models-flu-1.6.3-2").resolve()
    if fallback.is_dir():
        return fallback
    return None


def _resolve_monkeypox_vadr_model_dir(project_root: Path) -> Path | None:
    primary = (project_root / "soft" / "vadr-models-mpxv").resolve()
    if primary.is_dir():
        return primary
    return None


def _parse_info_field(info_text: str) -> dict[str, str]:
    info_map: dict[str, str] = {}
    for item in str(info_text or "").split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            info_map[key] = value
        else:
            info_map[item] = "true"
    return info_map


def _first_float(value: str | None, default: float = 0.0) -> float:
    if not value:
        return default
    token = str(value).split(",")[0]
    try:
        return float(token)
    except ValueError:
        return default


def _first_int(value: str | None, default: int = 0) -> int:
    if not value:
        return default
    token = str(value).split(",")[0]
    try:
        return int(float(token))
    except ValueError:
        return default


def _choose_best_snpeff_ann(ann_value: str | None) -> dict[str, str]:
    empty = {
        "annotation": "",
        "impact": "",
        "gene_name": "",
        "gene_id": "",
        "feature_id": "",
        "hgvs_c": "",
        "hgvs_p": "",
        "aa": "",
        "messages": "",
    }
    if not ann_value:
        return empty

    impact_rank = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3, "": 4}
    best_parts: list[str] | None = None
    best_key: tuple[int, int, int] | None = None
    for entry in str(ann_value).split(","):
        parts = entry.split("|")
        parts += [""] * (16 - len(parts))
        key = (
            impact_rank.get(parts[2], 9),
            0 if parts[5] == "transcript" else 1,
            0 if "intergenic" not in parts[1] else 1,
        )
        if best_key is None or key < best_key:
            best_key = key
            best_parts = parts[:16]
    if best_parts is None:
        return empty
    return {
        "annotation": best_parts[1],
        "impact": best_parts[2],
        "gene_name": best_parts[3],
        "gene_id": best_parts[4],
        "feature_id": best_parts[6],
        "hgvs_c": best_parts[9],
        "hgvs_p": best_parts[10],
        "aa": best_parts[13],
        "messages": best_parts[15],
    }


def _classify_variant_quality(qual: float, dp: int, maf: float) -> str:
    return "高质量突变" if qual > 10 and dp > 10 and maf > 0.1 else "低质量突变"


def _infer_influenza_resistance_gene(row: dict[str, object]) -> str:
    chrom = str(row.get("染色体") or "").upper()
    gene = str(row.get("基因") or "").upper()
    transcript = str(row.get("转录本ID") or "").upper()
    if "_NA" in chrom or "_NA_" in gene or "_NA_" in transcript:
        return "NA"
    if "_PA" in chrom or "_PA_" in gene or "_PA_" in transcript:
        return "PA"
    if "_MP" in chrom or "_MP_" in gene or "_MP_" in transcript:
        if "_715_982" in gene or "_715_982" in transcript:
            return "M2"
        return "M1"
    return ""


def _hgvs_p_to_short_aa_change(hgvs_p: object) -> str:
    text = str(hgvs_p or "").strip()
    matched = re.match(r"^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter|Stop)$", text)
    if not matched:
        return ""
    ref3, pos, alt3 = matched.groups()
    ref1 = AA3_TO_1.get(ref3, "")
    alt1 = AA3_TO_1.get(alt3, "")
    if not ref1 or not alt1:
        return ""
    return f"{ref1}{pos}{alt1}"


def _run_influenza_resistance_annotation(
    report_dir: Path,
    mutation_table_json: Path,
    logf=None,
) -> dict[str, str]:
    if not mutation_table_json.is_file() or mutation_table_json.stat().st_size == 0:
        return {"status": "missing", "note": "未找到流感 snpEff 变异注释 JSON", "table_json": "", "table_tsv": ""}

    rules_path = _project_root() / "database" / "virus" / "influenza_resistance_db_pack" / "influenza_resistance_rules.tsv"
    if not rules_path.is_file():
        return {"status": "missing", "note": "未找到流感耐药规则表", "table_json": "", "table_tsv": ""}

    try:
        payload = json.loads(mutation_table_json.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {"status": "failed", "note": "流感 snpEff 变异注释 JSON 解析失败", "table_json": "", "table_tsv": ""}
    raw_rows = payload.get("rows") if isinstance(payload.get("rows"), list) else []
    if not raw_rows or not isinstance(raw_rows[0], dict):
        return {"status": "empty", "note": "未找到可用于耐药注释的流感突变记录", "table_json": "", "table_tsv": ""}

    detected_by_gene: dict[str, dict[str, dict[str, object]]] = defaultdict(dict)
    for row in raw_rows:
        if str(row.get("质量分层") or "").strip() != "高质量突变":
            continue
        gene = _infer_influenza_resistance_gene(row)
        aa_change = _hgvs_p_to_short_aa_change(row.get("HGVS.p"))
        if not gene or not aa_change:
            continue
        detected_by_gene[gene][aa_change] = row

    output_json = report_dir / "snps.filt1.resistance_annotation.json"
    output_tsv = report_dir / "snps.filt1.resistance_annotation.tsv"
    if not detected_by_gene:
        empty_payload = {"status": "ready", "total_hits": 0, "rows": []}
        output_json.write_text(json.dumps(empty_payload, ensure_ascii=False, indent=2), encoding="utf-8")
        with output_tsv.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["基因", "规则突变", "命中突变", "药物类别", "药物", "风险等级", "权威结论", "注释说明", "证据来源", "证据强度", "适用范围", "备注"])
        return {"status": "ready", "note": "未命中流感耐药规则", "table_json": str(output_json.resolve()), "table_tsv": str(output_tsv.resolve())}

    result_rows: list[dict[str, object]] = []
    try:
        with rules_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for rule in reader:
                gene = str(rule.get("gene") or "").strip()
                mutation = str(rule.get("mutation") or "").strip()
                if not gene or not mutation or gene not in detected_by_gene:
                    continue
                required_mutations = [item.strip() for item in mutation.split("+") if item.strip()]
                if not required_mutations or not all(item in detected_by_gene[gene] for item in required_mutations):
                    continue
                matched_rows = [detected_by_gene[gene][item] for item in required_mutations]
                maf_values = []
                dp_values = []
                for matched_row in matched_rows:
                    try:
                        maf_values.append(float(matched_row.get("突变频率MAF") or 0))
                    except (TypeError, ValueError):
                        pass
                    try:
                        dp_values.append(int(float(matched_row.get("测序深度DP") or 0)))
                    except (TypeError, ValueError):
                        pass
                result_rows.append(
                    {
                        "基因": gene,
                        "规则突变": mutation,
                        "命中突变": ";".join(required_mutations),
                        "药物类别": str(rule.get("drug_class") or "").strip(),
                        "药物": str(rule.get("drugs") or "").strip(),
                        "风险等级": str(rule.get("effect_level") or "").strip(),
                        "权威结论": str(rule.get("cdc_or_who_call") or "").strip(),
                        "注释说明": str(rule.get("report_text") or "").strip(),
                        "证据来源": str(rule.get("evidence_source") or "").strip(),
                        "证据强度": str(rule.get("evidence_strength") or "").strip(),
                        "适用范围": str(rule.get("virus_scope") or "").strip(),
                        "备注": str(rule.get("notes") or "").strip(),
                        "平均突变频率MAF": round(sum(maf_values) / len(maf_values), 4) if maf_values else "",
                        "最小测序深度DP": min(dp_values) if dp_values else "",
                    }
                )
    except OSError:
        return {"status": "failed", "note": "读取流感耐药规则表失败", "table_json": "", "table_tsv": ""}

    fieldnames = [
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
        "平均突变频率MAF",
        "最小测序深度DP",
    ]
    with output_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(result_rows)
    output_json.write_text(
        json.dumps({"status": "ready", "total_hits": len(result_rows), "rows": result_rows}, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    return {
        "status": "ready",
        "note": f"已生成流感耐药突变注释表（命中 {len(result_rows)} 条）",
        "table_json": str(output_json.resolve()),
        "table_tsv": str(output_tsv.resolve()),
    }


def _run_influenza_snpeff_annotation(
    pre: str,
    consensus_fasta: Path,
    gff_path: Path,
    filtered_vcf: Path,
    blast_dir: Path,
    logf=None,
) -> dict[str, str]:
    if not gff_path.is_file() or gff_path.stat().st_size == 0:
        return {"status": "missing", "note": "未生成可用于 snpEff 的 VADR GFF3", "annotated_vcf": "", "table_json": ""}
    if not consensus_fasta.is_file() or consensus_fasta.stat().st_size == 0:
        return {"status": "missing", "note": "未找到流感 consensus FASTA，无法执行 snpEff 注释", "annotated_vcf": "", "table_json": ""}
    if not filtered_vcf.is_file() or filtered_vcf.stat().st_size == 0:
        return {"status": "missing", "note": "未找到 snps.filt1.vcf，无法执行 snpEff 注释", "annotated_vcf": "", "table_json": ""}

    project_root = _project_root()
    snpeff_jar = project_root / "snpEff" / "snpEff.jar"
    if not snpeff_jar.is_file():
        return {"status": "missing", "note": "未找到 snpEff.jar", "annotated_vcf": "", "table_json": ""}

    snpeff_dir = blast_dir / "snpeff_vadr_ref"
    ref_dir = snpeff_dir / "ref"
    ref_dir.mkdir(parents=True, exist_ok=True)
    sequences_fa = ref_dir / "sequences.fa"
    genes_gff = ref_dir / "genes.gff"
    shutil.copy2(consensus_fasta, sequences_fa)
    shutil.copy2(gff_path, genes_gff)

    config_path = snpeff_dir / "snpEff.config"
    config_path.write_text(f"ref.genome : {pre}_vadr\n", encoding="utf-8")
    predictor_bin = ref_dir / "snpEffectPredictor.bin"
    snpeff_jar_abs = snpeff_jar.resolve()
    config_path_abs = config_path.resolve()
    snpeff_dir_abs = snpeff_dir.resolve()
    filtered_vcf_abs = filtered_vcf.resolve()
    if not predictor_bin.is_file():
        build_cmd = " ".join(
            [
                "java",
                "-Xmx4g",
                "-jar",
                shlex.quote(str(snpeff_jar_abs)),
                "build",
                "-gff3",
                "-noCheckCds",
                "-noCheckProtein",
                "-c",
                shlex.quote(str(config_path_abs)),
                "-dataDir",
                shlex.quote(str(snpeff_dir_abs)),
                "ref",
            ]
        )
        run_command(build_cmd, logf=logf)
    if not predictor_bin.is_file():
        return {"status": "failed", "note": "snpEff 数据库构建失败", "annotated_vcf": "", "table_json": ""}

    annotated_vcf = filtered_vcf.with_name("snps.anno.vcf")
    sample_annotated_vcf = filtered_vcf.with_name(f"{pre}.anno.vcf")
    table_tsv = filtered_vcf.with_name(filtered_vcf.stem + ".mutation_table.tsv")
    table_json = filtered_vcf.with_name(filtered_vcf.stem + ".mutation_table.json")
    if annotated_vcf.exists() or annotated_vcf.is_symlink():
        annotated_vcf.unlink()
    ann_cmd = [
        "java",
        "-Xmx4g",
        "-jar",
        str(snpeff_jar_abs),
        "ann",
        "-noStats",
        "-c",
        str(config_path_abs),
        "-dataDir",
        str(snpeff_dir_abs),
        "ref",
        str(filtered_vcf_abs),
    ]
    if logf is not None:
        logf.write(f"\n[CMD] {' '.join(shlex.quote(part) for part in ann_cmd)} > {shlex.quote(str(annotated_vcf))}\n")
        logf.flush()
    with annotated_vcf.open("w", encoding="utf-8") as handle:
        completed = subprocess.run(
            ann_cmd,
            stdout=handle,
            stderr=logf,
            check=False,
        )
    if logf is not None:
        logf.write(f"[CMD_EXIT] code={completed.returncode}\n")
        logf.flush()
    if completed.returncode != 0:
        raise RuntimeError(f"snpEff ann 执行失败(returncode={completed.returncode})")
    if not annotated_vcf.is_file() or annotated_vcf.stat().st_size == 0:
        return {"status": "failed", "note": "snpEff 未生成注释 VCF", "annotated_vcf": "", "table_json": ""}
    if sample_annotated_vcf.name != annotated_vcf.name:
        try:
            if sample_annotated_vcf.exists() or sample_annotated_vcf.is_symlink():
                sample_annotated_vcf.unlink()
            os.symlink(annotated_vcf.name, sample_annotated_vcf)
        except OSError:
            shutil.copy2(annotated_vcf, sample_annotated_vcf)

    rows: list[dict[str, object]] = []
    total = 0
    high_quality = 0
    with annotated_vcf.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            total += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            chrom, pos, record_id, ref, alt, qual_raw, filt, info = fields[:8]
            info_map = _parse_info_field(info)
            ann = _choose_best_snpeff_ann(info_map.get("ANN"))
            qual = _first_float(qual_raw, 0.0)
            dp = _first_int(info_map.get("DP"), 0)
            ao = _first_float(info_map.get("AO"), 0.0)
            ro = _first_float(info_map.get("RO"), 0.0)
            allele_depth = ao + ro
            if allele_depth > 0:
                maf = round((ao / allele_depth), 6)
            else:
                maf = _first_float(info_map.get("AF"), -1.0)
                if maf < 0 and dp > 0:
                    maf = round((ao / float(dp)), 6)
                elif maf < 0:
                    maf = 0.0
            quality_label = _classify_variant_quality(qual, dp, maf)
            if quality_label == "高质量突变":
                high_quality += 1
            rows.append(
                {
                    "染色体": chrom,
                    "位置": int(pos),
                    "ID": record_id,
                    "参考碱基": ref,
                    "突变碱基": alt,
                    "核苷酸突变": f"{ref}{pos}{alt}",
                    "变异类型": info_map.get("TYPE", ""),
                    "质量值QUAL": qual,
                    "测序深度DP": dp,
                    "突变频率MAF": round(maf, 6),
                    "质量分层": quality_label,
                    "注释效应": ann["annotation"],
                    "影响等级": ann["impact"],
                    "基因": ann["gene_name"],
                    "基因ID": ann["gene_id"],
                    "转录本ID": ann["feature_id"],
                    "HGVS.c": ann["hgvs_c"],
                    "HGVS.p": ann["hgvs_p"],
                    "氨基酸位点": ann["aa"],
                    "警告信息": ann["messages"],
                    "FILTER": filt,
                }
            )

    fieldnames = [
        "染色体",
        "位置",
        "ID",
        "参考碱基",
        "突变碱基",
        "核苷酸突变",
        "变异类型",
        "质量值QUAL",
        "测序深度DP",
        "突变频率MAF",
        "质量分层",
        "注释效应",
        "影响等级",
        "基因",
        "基因ID",
        "转录本ID",
        "HGVS.c",
        "HGVS.p",
        "氨基酸位点",
        "警告信息",
        "FILTER",
    ]
    with table_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    with table_json.open("w", encoding="utf-8") as handle:
        json.dump(
            {
                "status": "ready",
                "source_vcf": str(annotated_vcf.absolute()),
                "total_variants": total,
                "high_quality_variants": high_quality,
                "low_quality_variants": total - high_quality,
                "rows": rows,
            },
            handle,
            ensure_ascii=False,
            indent=2,
        )
    return {
        "status": "ready",
        "note": "已生成流感 snpEff 变异注释表",
        "annotated_vcf": str(annotated_vcf.absolute()),
        "table_json": str(table_json.resolve()),
    }


def _run_vadr_flu_annotation(pre: str, final_fasta: Path, blast_dir: Path, threads: int = 2, logf=None) -> dict[str, str]:
    project_root = _project_root()
    model_dir = _resolve_flu_vadr_model_dir(project_root)
    if model_dir is None:
        return {"status": "missing", "note": "未找到流感 VADR 模型目录", "output_dir": "", "gff_path": ""}
    env = _build_vadr_env(project_root, model_dir)
    if env is None:
        return {"status": "missing", "note": "VADR 运行依赖不完整", "output_dir": "", "gff_path": ""}

    vadr_root = blast_dir / "vadr"
    vadr_root.mkdir(exist_ok=True)
    output_dir = vadr_root / f"{pre}_vadr"
    prefix = f"{output_dir.name}.vadr"
    pass_tbl = output_dir / f"{prefix}.pass.tbl"
    fail_tbl = output_dir / f"{prefix}.fail.tbl"
    gff_path = vadr_root / f"{pre}.vadr.gff3"
    gff_ready = gff_path.is_file() and gff_path.stat().st_size > 0
    if output_dir.is_dir() and (pass_tbl.is_file() or fail_tbl.is_file()):
        if pass_tbl.is_file() and not gff_ready:
            annotate_tbl2gff = project_root / "soft" / "vadr" / "miniscripts" / "annotate-tbl2gff.pl"
            try:
                run_command(
                    f"perl {shlex.quote(str(annotate_tbl2gff))} {shlex.quote(str(pass_tbl))} > {shlex.quote(str(gff_path))}",
                    logf=logf,
                    env=env,
                )
            except Exception:
                pass
            gff_ready = gff_path.is_file() and gff_path.stat().st_size > 0
        return {
            "status": "ready",
            "note": "已检测到现有 VADR 注释结果",
            "output_dir": str(output_dir.resolve()),
            "gff_path": str(gff_path.resolve()) if gff_ready else "",
        }

    cpu = max(1, min(int(threads or 1), 2))
    vadr_script = project_root / "soft" / "vadr" / "v-annotate.pl"
    cmd = " ".join(
        [
            "perl",
            shlex.quote(str(vadr_script)),
            "--split",
            "--cpu",
            str(cpu),
            "-r",
            "--atgonly",
            "--alt_fail",
            "extrant5,extrant3",
            "--xnocomp",
            "--nomisc",
            "--forcegene",
            "--mkey",
            "flu",
            "--mdir",
            shlex.quote(str(model_dir)),
            "-f",
            shlex.quote(str(final_fasta)),
            shlex.quote(str(output_dir)),
        ]
    )
    try:
        run_command(cmd, logf=logf, env=norovirus_env)
    except Exception as exc:
        return {
            "status": "failed",
            "note": f"VADR 注释失败: {exc}",
            "output_dir": str(output_dir.resolve()) if output_dir.exists() else "",
            "gff_path": "",
        }

    if pass_tbl.is_file():
        annotate_tbl2gff = project_root / "soft" / "vadr" / "miniscripts" / "annotate-tbl2gff.pl"
        try:
            run_command(
                f"perl {shlex.quote(str(annotate_tbl2gff))} {shlex.quote(str(pass_tbl))} > {shlex.quote(str(gff_path))}",
                logf=logf,
                env=env,
            )
        except Exception:
            pass
    gff_ready = gff_path.is_file() and gff_path.stat().st_size > 0
    return {
        "status": "ready",
        "note": "VADR 注释完成",
        "output_dir": str(output_dir.resolve()),
        "gff_path": str(gff_path.resolve()) if gff_ready else "",
    }


def _run_vadr_monkeypox_annotation(pre: str, input_fasta: Path, output_root: Path, logf=None) -> dict[str, str]:
    project_root = _project_root()
    model_dir = _resolve_monkeypox_vadr_model_dir(project_root)
    if model_dir is None:
        return {"status": "missing", "note": "未找到猴痘 VADR 模型目录", "output_dir": "", "gff_path": ""}
    env = _build_vadr_env(project_root, model_dir)
    if env is None:
        return {"status": "missing", "note": "VADR 运行依赖不完整", "output_dir": "", "gff_path": ""}

    vadr_root = output_root / "vadr"
    vadr_root.mkdir(parents=True, exist_ok=True)
    output_dir = vadr_root / f"{pre}_vadr"
    prefix = f"{output_dir.name}.vadr"
    pass_tbl = output_dir / f"{prefix}.pass.tbl"
    fail_tbl = output_dir / f"{prefix}.fail.tbl"
    gff_path = vadr_root / f"{pre}.vadr.gff3"
    gff_ready = gff_path.is_file() and gff_path.stat().st_size > 0
    if output_dir.is_dir() and (pass_tbl.is_file() or fail_tbl.is_file()):
        if pass_tbl.is_file() and not gff_ready:
            annotate_tbl2gff = project_root / "soft" / "vadr" / "miniscripts" / "annotate-tbl2gff.pl"
            try:
                run_command(
                    f"perl {shlex.quote(str(annotate_tbl2gff))} {shlex.quote(str(pass_tbl))} > {shlex.quote(str(gff_path))}",
                    logf=logf,
                    env=env,
                )
            except Exception:
                pass
            gff_ready = gff_path.is_file() and gff_path.stat().st_size > 0
        return {
            "status": "ready",
            "note": "已检测到现有猴痘 VADR 注释结果",
            "output_dir": str(output_dir.resolve()),
            "gff_path": str(gff_path.resolve()) if gff_ready else "",
        }

    vadr_script = project_root / "soft" / "vadr" / "v-annotate.pl"
    cmd = " ".join(
        [
            "perl",
            shlex.quote(str(vadr_script)),
            "--split",
            "--cpu",
            "1",
            "--glsearch",
            "--minimap2",
            "-s",
            "-r",
            "--nomisc",
            "--r_lowsimok",
            "--r_lowsimxd",
            "100",
            "--r_lowsimxl",
            "2000",
            "--alt_pass",
            "discontn,dupregin",
            "--s_overhang",
            "150",
            "--mkey",
            "mpxv",
            "--mdir",
            shlex.quote(str(model_dir)),
            "-f",
            shlex.quote(str(input_fasta)),
            shlex.quote(str(output_dir)),
        ]
    )
    try:
        run_command(cmd, logf=logf, env=env)
    except Exception as exc:
        return {
            "status": "failed",
            "note": f"猴痘 VADR 注释失败: {exc}",
            "output_dir": str(output_dir.resolve()) if output_dir.exists() else "",
            "gff_path": "",
        }

    if pass_tbl.is_file():
        annotate_tbl2gff = project_root / "soft" / "vadr" / "miniscripts" / "annotate-tbl2gff.pl"
        try:
            run_command(
                f"perl {shlex.quote(str(annotate_tbl2gff))} {shlex.quote(str(pass_tbl))} > {shlex.quote(str(gff_path))}",
                logf=logf,
                env=env,
            )
        except Exception:
            pass
    gff_ready = gff_path.is_file() and gff_path.stat().st_size > 0
    return {
        "status": "ready",
        "note": "猴痘 VADR 注释完成",
        "output_dir": str(output_dir.resolve()),
        "gff_path": str(gff_path.resolve()) if gff_ready else "",
    }


def prepare_monkeypox_reference_annotation(pre: str, reference_fasta: Path, output_root: Path, logf=None) -> Path | None:
    result = _run_vadr_monkeypox_annotation(pre, reference_fasta, output_root, logf=logf)
    gff_path = str(result.get("gff_path") or "").strip()
    if result.get("status") == "ready" and gff_path:
        gff_file = Path(gff_path)
        if gff_file.is_file() and gff_file.stat().st_size > 0:
            return gff_file
    return None


def _resolve_influenza_screening_sources() -> dict[str, Path]:
    project_root = _project_root()
    database_root = str(os.environ.get("META_DATABASE_ROOT") or "").strip()
    candidates = {
        "ha": [
            Path(str(os.environ.get("META_FLUA_HA_DB") or "")).expanduser() if str(os.environ.get("META_FLUA_HA_DB") or "").strip() else None,
            Path(database_root).expanduser() / "virus" / "influenza_a" / "ha_subtypes.dedup_98.fa" if database_root else None,
            project_root / "database" / "virus" / "influenza_a" / "ha_subtypes.dedup_98.fa",
            Path(database_root).expanduser() / "virus" / "influenza_a" / "ha_subtypes.fa" if database_root else None,
            project_root / "database" / "virus" / "influenza_a" / "ha_subtypes.fa",
        ],
        "na": [
            Path(str(os.environ.get("META_FLUA_NA_DB") or "")).expanduser() if str(os.environ.get("META_FLUA_NA_DB") or "").strip() else None,
            Path(database_root).expanduser() / "virus" / "influenza_a" / "na_subtypes.dedup_98.fa" if database_root else None,
            project_root / "database" / "virus" / "influenza_a" / "na_subtypes.dedup_98.fa",
            Path(database_root).expanduser() / "virus" / "influenza_a" / "na_subtypes.fa" if database_root else None,
            project_root / "database" / "virus" / "influenza_a" / "na_subtypes.fa",
        ],
        "irma": [
            Path(str(os.environ.get("META_FLU_IRMA_REFERENCE_SET") or "")).expanduser() if str(os.environ.get("META_FLU_IRMA_REFERENCE_SET") or "").strip() else None,
            Path(database_root).expanduser() / "virus" / "influenza" / "consensus_irma.fasta" if database_root else None,
            project_root / "database" / "virus" / "influenza" / "consensus_irma.fasta",
            project_root / "public" / "wf-flu-master" / "data" / "primer_schemes" / "V1" / "consensus_irma.fasta",
        ],
    }
    resolved: dict[str, Path] = {}
    for key, paths in candidates.items():
        found = next((path.resolve() for path in paths if path and path.is_file()), None)
        if found is None:
            raise FileNotFoundError(f"未找到流感参考数据库: {key}")
        resolved[key] = found
    return resolved


def _read_fasta_records(path: Path) -> list:
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        return list(SeqIO.parse(handle, "fasta"))


def _segment_group_from_reference_id(reference_id: str) -> str:
    text = str(reference_id or "").strip()
    parts = text.split("_")
    if len(parts) >= 2:
        return parts[1]
    return text


def _infer_subtype_from_header(text: str, prefix: str) -> str:
    header = str(text or "").strip()
    match = re.search(rf"(?<![A-Za-z0-9])({prefix.upper()}[0-9]+)(?![A-Za-z0-9])", header, re.IGNORECASE)
    if match:
        return match.group(1).upper()
    parts = [part.strip() for part in header.split("|")]
    for part in parts:
        if re.fullmatch(rf"{prefix.upper()}[0-9]+", part.upper()):
            return part.upper()
    return "-"


def _normalize_subtype_record(record, segment_group: str) -> tuple:
    raw_id = str(record.id or "").strip()
    subtype = _infer_subtype_from_header(raw_id, "H" if segment_group == "HA" else "N")
    if re.fullmatch(r"A_(HA|NA)_[A-Z0-9-]+__.+", raw_id):
        record.id = raw_id
        record.name = raw_id
        record.description = raw_id
        return record, subtype
    accession = re.sub(r"[^A-Za-z0-9._-]+", "_", raw_id.split("|")[0] or raw_id) or "unknown"
    normalized_id = f"A_{segment_group}_{subtype}__{accession}"
    record.id = normalized_id
    record.name = normalized_id
    record.description = normalized_id
    return record, subtype


def _build_influenza_type_reference_set(work_dir: Path) -> tuple[Path, dict]:
    sources = _resolve_influenza_screening_sources()
    work_dir.mkdir(parents=True, exist_ok=True)
    screening_fasta = work_dir / "irma_type_candidates.fa"
    metadata: dict[str, dict] = {}
    records = []

    irma_records = {str(record.id): record for record in _read_fasta_records(sources["irma"])}
    a_core = ["A_PB2", "A_PB1", "A_PA", "A_NP", "A_MP", "A_NS"]
    b_segments = ["B_PB2", "B_PB1", "B_PA", "B_NP", "B_MP", "B_NS", "B_HA", "B_NA"]
    for segment_id in a_core + b_segments:
        if segment_id not in irma_records:
            raise FileNotFoundError(f"IRMA 参考集中缺少片段: {segment_id}")
        record = irma_records[segment_id]
        records.append(record)
        metadata[segment_id] = {
            "reference_id": segment_id,
            "influenza_type": "Influenza A virus" if segment_id.startswith("A_") else "Influenza B virus",
            "segment_group": _segment_group_from_reference_id(segment_id),
            "subtype": "-" if "__" not in segment_id else segment_id.split("__", 1)[0].split("_")[-1],
            "source": "irma_consensus",
        }

    with screening_fasta.open("w", encoding="utf-8") as handle:
        SeqIO.write(records, handle, "fasta")
    return screening_fasta, metadata


def _build_influenza_subtype_reference_set(work_dir: Path) -> tuple[Path, dict]:
    sources = _resolve_influenza_screening_sources()
    work_dir.mkdir(parents=True, exist_ok=True)
    screening_fasta = work_dir / "subtype_candidates.fa"
    metadata: dict[str, dict] = {}
    records = []
    for record in _read_fasta_records(sources["ha"]):
        normalized, subtype = _normalize_subtype_record(record, "HA")
        records.append(normalized)
        metadata[normalized.id] = {
            "reference_id": normalized.id,
            "influenza_type": "Influenza A virus",
            "segment_group": "HA",
            "subtype": subtype,
            "source": "ha_subtype",
        }
    for record in _read_fasta_records(sources["na"]):
        normalized, subtype = _normalize_subtype_record(record, "NA")
        records.append(normalized)
        metadata[normalized.id] = {
            "reference_id": normalized.id,
            "influenza_type": "Influenza A virus",
            "segment_group": "NA",
            "subtype": subtype,
            "source": "na_subtype",
        }
    with screening_fasta.open("w", encoding="utf-8") as handle:
        SeqIO.write(records, handle, "fasta")
    return screening_fasta, metadata


def _run_influenza_screening_alignment(screening_fasta: Path, output_dir: Path, single_fastq: str = "", fq1: str = "", fq2: str = "", long_type: str = "", threads: int = 4, logf=None) -> tuple[Path, Path, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    bam_path = output_dir / "screening.bam"
    idxstats_path = output_dir / "screening.idxstats.tsv"
    depth_path = output_dir / "screening.depth.tsv"
    if Path('10239.1.fastq').is_file() and Path('10239.1.fastq').stat().st_size > 0:
        run_command(f"bwa index {shlex.quote(str(screening_fasta))}", logf=logf)
        if str(fq2 or "").strip():
            run_command(
                f"bwa mem -t {threads} {shlex.quote(str(screening_fasta))} 10239.1.fastq 10239.2.fastq | samtools sort -o {shlex.quote(str(bam_path))}",
                logf=logf,
            )
        else:
            run_command(
                f"bwa mem -t {threads}  {shlex.quote(str(screening_fasta))} 10239.1.fastq | samtools sort -o {shlex.quote(str(bam_path))}",
                logf=logf,
            )
    else:
        source_fastq = str(fq1 or single_fastq or "").strip()
        if not source_fastq:
            raise FileNotFoundError("未找到可用于流感初筛比对的 FASTQ 输入。")
        platform = "pb" if "pacbio" in _normalize_species_label(long_type) else "ont"
        run_command(
            f"minimap2 -ax map-{platform} {shlex.quote(str(screening_fasta))} {shlex.quote(source_fastq)} -t {threads} | samtools sort -o {shlex.quote(str(bam_path))}",
            logf=logf,
        )
    run_command(f"samtools index {shlex.quote(str(bam_path))}", logf=logf)
    run_command(f"samtools idxstats {shlex.quote(str(bam_path))} > {shlex.quote(str(idxstats_path))}", logf=logf)
    run_command(f"samtools depth -aa {shlex.quote(str(bam_path))} > {shlex.quote(str(depth_path))}", logf=logf)
    return bam_path, idxstats_path, depth_path


def _collect_depth_stats(depth_path: Path) -> dict[str, dict]:
    stats = defaultdict(lambda: {"covered_bases": 0, "depth_sum": 0.0, "length": 0})
    with depth_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            parts = line.rstrip().split("\t")
            if len(parts) < 3:
                continue
            contig = parts[0]
            depth = float(parts[2] or 0)
            stats[contig]["length"] += 1
            stats[contig]["depth_sum"] += depth
            if depth > 0:
                stats[contig]["covered_bases"] += 1
    return stats


def _collect_idxstats(idxstats_path: Path) -> dict[str, int]:
    mapped = {}
    with idxstats_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            parts = line.rstrip().split("\t")
            if len(parts) < 4:
                continue
            if parts[0] == "*":
                continue
            mapped[parts[0]] = int(float(parts[2] or 0))
    return mapped


def _write_influenza_segment_stats(output_dir: Path, metadata: dict, idxstats_path: Path, depth_path: Path) -> list[dict]:
    mapped = _collect_idxstats(idxstats_path)
    depth_stats = _collect_depth_stats(depth_path)
    rows = []
    for reference_id, meta in metadata.items():
        length = int(depth_stats.get(reference_id, {}).get("length", 0))
        covered = int(depth_stats.get(reference_id, {}).get("covered_bases", 0))
        depth_sum = float(depth_stats.get(reference_id, {}).get("depth_sum", 0.0))
        coverage_pct = round((covered / length) * 100, 2) if length else 0.0
        mean_depth = round(depth_sum / length, 2) if length else 0.0
        mapped_reads = int(mapped.get(reference_id, 0))
        rows.append({
            "reference_id": reference_id,
            "influenza_type": meta["influenza_type"],
            "segment_group": meta["segment_group"],
            "subtype": meta["subtype"],
            "source": meta["source"],
            "length": length,
            "covered_bases": covered,
            "coverage_pct": coverage_pct,
            "mean_depth": mean_depth,
            "mapped_reads": mapped_reads,
        })
    rows.sort(key=lambda item: (item["influenza_type"], item["segment_group"], item["coverage_pct"], item["mapped_reads"]), reverse=True)
    stats_path = output_dir / "segment_stats.tsv"
    with stats_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()) if rows else ["reference_id", "influenza_type", "segment_group", "subtype", "source", "length", "covered_bases", "coverage_pct", "mean_depth", "mapped_reads"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return rows


def _summarize_type_support(rows: list[dict], influenza_type: str) -> dict:
    best_by_group = {}
    for row in rows:
        if row["influenza_type"] != influenza_type:
            continue
        group = row["segment_group"]
        current = best_by_group.get(group)
        score = (float(row["coverage_pct"]), float(row["mean_depth"]), int(row["mapped_reads"]))
        if current is None or score > current[0]:
            best_by_group[group] = (score, row)
    selected_rows = [item[1] for item in best_by_group.values()]
    present_count = sum(1 for row in selected_rows if row["coverage_pct"] >= 20.0 or row["mapped_reads"] >= 20)
    coverage_sum = round(sum(float(row["coverage_pct"]) for row in selected_rows), 2)
    mapped_sum = sum(int(row["mapped_reads"]) for row in selected_rows)
    return {
        "rows": selected_rows,
        "present_count": present_count,
        "coverage_sum": coverage_sum,
        "mapped_sum": mapped_sum,
        "best_by_group": {row["segment_group"]: row for row in selected_rows},
    }


def _determine_influenza_type_from_segments(rows: list[dict]) -> str:
    a_support = _summarize_type_support(rows, "Influenza A virus")
    b_support = _summarize_type_support(rows, "Influenza B virus")
    if a_support["present_count"] < 4 and b_support["present_count"] < 4:
        return "Other"
    if a_support["present_count"] > b_support["present_count"]:
        return "Influenza A virus"
    if b_support["present_count"] > a_support["present_count"]:
        return "Influenza B virus"
    if a_support["coverage_sum"] > b_support["coverage_sum"] * 1.1:
        return "Influenza A virus"
    if b_support["coverage_sum"] > a_support["coverage_sum"] * 1.1:
        return "Influenza B virus"
    if a_support["mapped_sum"] > b_support["mapped_sum"]:
        return "Influenza A virus"
    if b_support["mapped_sum"] > a_support["mapped_sum"]:
        return "Influenza B virus"
    return "Other"


def _write_final_influenza_reference_set(output_dir: Path, rows: list[dict], metadata: dict, reference_fastas, influenza_type: str, sample_name: str) -> dict:
    fasta_list = reference_fastas if isinstance(reference_fastas, (list, tuple)) else [reference_fastas]
    records = {}
    for fasta_path in fasta_list:
        for record in _read_fasta_records(Path(fasta_path)):
            records[str(record.id)] = record
    output_dir.mkdir(parents=True, exist_ok=True)
    final_reference = output_dir / f"{sample_name}.final_segments.fa"
    manifest_path = output_dir / "final_segments.tsv"
    selected_records = []
    manifest_rows = []
    if influenza_type == "Influenza A virus":
        a_support = _summarize_type_support(rows, influenza_type)["best_by_group"]
        chosen_ha = a_support.get("HA")
        chosen_na = a_support.get("NA")
        if chosen_ha is None or chosen_na is None:
            raise RuntimeError("甲流样本未能稳定选出 HA/NA 最优参考。")
        ordered_ids = [
            "A_PB2",
            "A_PB1",
            "A_PA",
            "A_NP",
            chosen_ha["reference_id"],
            "A_MP",
            chosen_na["reference_id"],
            "A_NS",
        ]
        subtype_call = f"{chosen_ha['subtype']}{chosen_na['subtype']}".replace("-", "")
        ha_subtype = chosen_ha["subtype"]
        na_subtype = chosen_na["subtype"]
    else:
        ordered_ids = ["B_PB2", "B_PB1", "B_PA", "B_NP", "B_HA", "B_MP", "B_NA", "B_NS"]
        subtype_call = "-"
        ha_subtype = "-"
        na_subtype = "-"
    for reference_id in ordered_ids:
        record = records[reference_id]
        selected_records.append(record)
        manifest_rows.append({
            "segment_order": len(manifest_rows) + 1,
            "reference_id": reference_id,
            "influenza_type": influenza_type,
            "segment_group": metadata[reference_id]["segment_group"],
            "subtype": metadata[reference_id]["subtype"],
            "source": metadata[reference_id]["source"],
        })
    with final_reference.open("w", encoding="utf-8") as handle:
        SeqIO.write(selected_records, handle, "fasta")
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)
    return {
        "reference_path": str(final_reference.resolve()),
        "manifest_path": str(manifest_path.resolve()),
        "ha_subtype": ha_subtype,
        "na_subtype": na_subtype,
        "subtype_call": subtype_call,
    }


def _load_existing_influenza_reference_result(pre: str) -> dict | None:
    summary_path = Path("wf_flu") / "typing_summary.tsv"
    reference_dir = Path("wf_flu") / "reference_sets"
    manifest_path = reference_dir / "final_segments.tsv"
    reference_path = reference_dir / f"{pre}.final_segments.fa"
    type_stats_path = Path("wf_flu") / "screening" / "type" / "segment_stats.tsv"
    subtype_stats_path = Path("wf_flu") / "screening" / "subtype" / "segment_stats.tsv"
    if not summary_path.is_file() or summary_path.stat().st_size == 0:
        return None
    with summary_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        row = next(reader, None)
    if not row:
        return None
    status = str(row.get("status") or "").strip()
    influenza_type = str(row.get("influenza_type") or "").strip() or "Other"
    if status == "screening_stop":
        return {
            "status": "other",
            "influenza_type": "Other",
            "reference_path": "",
            "ha_subtype": "-",
            "na_subtype": "-",
            "subtype_call": "-",
            "type_stats_path": str(type_stats_path.resolve()) if type_stats_path.is_file() else "",
            "summary_path": str(summary_path.resolve()),
        }
    if status != "ready":
        return None
    resolved_reference = str(row.get("reference_path") or "").strip()
    if not resolved_reference and reference_path.is_file():
        resolved_reference = str(reference_path.resolve())
    if not resolved_reference:
        return None
    return {
        "status": "ready",
        "influenza_type": influenza_type,
        "reference_path": resolved_reference,
        "manifest_path": str(manifest_path.resolve()) if manifest_path.is_file() else "",
        "ha_subtype": str(row.get("ha_subtype") or "-").strip() or "-",
        "na_subtype": str(row.get("na_subtype") or "-").strip() or "-",
        "subtype_call": str(row.get("subtype_call") or "-").strip() or "-",
        "type_stats_path": str(type_stats_path.resolve()) if type_stats_path.is_file() else "",
        "subtype_stats_path": str(subtype_stats_path.resolve()) if subtype_stats_path.is_file() else "",
        "summary_path": str(summary_path.resolve()),
    }


def prepare_influenza_reference_set(
    pre: str,
    species: str = "",
    requested_ref: str = "",
    single_fastq: str = "",
    fq1: str = "",
    fq2: str = "",
    long_type: str = "",
    threads: int = 4,
    logf=None,
) -> dict:
    requested = str(requested_ref or "").strip()
    requested_path = Path(requested).expanduser() if requested and requested != "noref" else None
    if requested_path and requested_path.is_file() and not _looks_like_typing_reference(requested_path):
        inferred_type = _infer_influenza_type_from_label(species)
        return {
            "status": "custom_reference",
            "influenza_type": inferred_type if inferred_type != "-" else "Influenza A virus",
            "reference_path": str(requested_path.resolve()),
            "ha_subtype": "-",
            "na_subtype": "-",
            "subtype_call": "-",
        }

    existing = _load_existing_influenza_reference_result(pre)
    if existing is not None:
        return existing

    screening_dir = Path("wf_flu") / "screening"
    reference_dir = Path("wf_flu") / "reference_sets"
    type_reference_fasta, type_metadata = _build_influenza_type_reference_set(reference_dir)
    virus_log_path = Path("virus_ana.log")
    with virus_log_path.open("a", encoding="utf-8") as virus_logf:
        active_logf = virus_logf
        if logf is not None and hasattr(logf, "write"):
            try:
                logf.write("[INFO] 流感初筛日志同步写入 virus_ana.log\n")
                logf.flush()
            except Exception:
                pass
        _, type_idxstats_path, type_depth_path = _run_influenza_screening_alignment(
            type_reference_fasta,
            screening_dir / "type",
            single_fastq=single_fastq,
            fq1=fq1,
            fq2=fq2,
            long_type=long_type,
            threads=threads,
            logf=active_logf,
        )
    type_rows = _write_influenza_segment_stats(screening_dir / "type", type_metadata, type_idxstats_path, type_depth_path)
    influenza_type = _determine_influenza_type_from_segments(type_rows)
    summary_path = Path("wf_flu") / "typing_summary.tsv"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    if influenza_type == "Other":
        with summary_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sample", "influenza_type", "ha_subtype", "na_subtype", "subtype_call", "status"])
            writer.writerow([pre, "Other", "-", "-", "-", "screening_stop"])
        return {
            "status": "other",
            "influenza_type": "Other",
            "reference_path": "",
            "ha_subtype": "-",
            "na_subtype": "-",
            "subtype_call": "-",
            "type_stats_path": str((screening_dir / "type" / "segment_stats.tsv").resolve()),
            "summary_path": str(summary_path.resolve()),
        }

    subtype_reference_fasta = None
    subtype_rows = []
    final_info = {}
    if influenza_type == "Influenza A virus":
        subtype_reference_fasta, subtype_metadata = _build_influenza_subtype_reference_set(reference_dir)
        with virus_log_path.open("a", encoding="utf-8") as virus_logf:
            _, subtype_idxstats_path, subtype_depth_path = _run_influenza_screening_alignment(
                subtype_reference_fasta,
                screening_dir / "subtype",
                single_fastq=single_fastq,
                fq1=fq1,
                fq2=fq2,
                long_type=long_type,
                threads=threads,
                logf=virus_logf,
            )
        subtype_rows = _write_influenza_segment_stats(screening_dir / "subtype", subtype_metadata, subtype_idxstats_path, subtype_depth_path)
        combined_rows = type_rows + subtype_rows
        combined_metadata = {**type_metadata, **subtype_metadata}
        final_info = _write_final_influenza_reference_set(reference_dir, combined_rows, combined_metadata, [type_reference_fasta, subtype_reference_fasta], influenza_type, pre)
    else:
        final_info = _write_final_influenza_reference_set(reference_dir, type_rows, type_metadata, [type_reference_fasta], influenza_type, pre)
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "influenza_type", "ha_subtype", "na_subtype", "subtype_call", "status", "reference_path"])
        writer.writerow([pre, influenza_type, final_info["ha_subtype"], final_info["na_subtype"], final_info["subtype_call"], "ready", final_info["reference_path"]])
    return {
        "status": "ready",
        "influenza_type": influenza_type,
        "reference_path": final_info["reference_path"],
        "manifest_path": final_info["manifest_path"],
        "ha_subtype": final_info["ha_subtype"],
        "na_subtype": final_info["na_subtype"],
        "subtype_call": final_info["subtype_call"],
        "type_stats_path": str((screening_dir / "type" / "segment_stats.tsv").resolve()),
        "subtype_stats_path": str((screening_dir / "subtype" / "segment_stats.tsv").resolve()) if influenza_type == "Influenza A virus" else "",
        "summary_path": str(summary_path.resolve()),
    }


def resolve_influenza_reference(
    pre: str,
    species: str = "",
    requested_ref: str = "",
    single_fastq: str = "",
    fq1: str = "",
    fq2: str = "",
    long_type: str = "",
    threads: int = 4,
    logf=None,
) -> str:
    result = prepare_influenza_reference_set(
        pre=pre,
        species=species,
        requested_ref=requested_ref,
        single_fastq=single_fastq,
        fq1=fq1,
        fq2=fq2,
        long_type=long_type,
        threads=threads,
        logf=logf,
    )
    return str(result.get("reference_path") or "").strip()


def _load_taxonomy_candidates(path: Path) -> list[dict]:
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows: list[dict] = []
        for row in reader:
            if not row:
                continue
            label = str(row.get("种") or row.get("物种") or row.get("name") or "").strip()
            if not label:
                continue
            ratio = float(row.get("比例") or row.get("相对丰度") or 0.0)
            reads = int(float(row.get("序列数量") or row.get("new_est_reads") or 0))
            rows.append({"label": label, "ratio": ratio, "reads": reads})
    rows.sort(key=lambda item: (item["reads"], item["ratio"]), reverse=True)
    return rows


def detect_influenza_type(pre: str, species: str = "") -> str:
    explicit = _infer_influenza_type_from_label(species)
    if explicit in {"Influenza A virus", "Influenza B virus", "Influenza C virus", "Influenza D virus"}:
        return explicit

    candidates: list[dict] = []
    for filename in (f"{pre}_2.list.txt", f"{pre}.list.txt", f"{pre}.taxonomy_summary.tsv"):
        candidates.extend(_load_taxonomy_candidates(Path(filename)))
    candidates.sort(key=lambda item: (item["reads"], item["ratio"]), reverse=True)
    for candidate in candidates:
        inferred = _infer_influenza_type_from_label(candidate["label"])
        if inferred != "-":
            return inferred
    return "-"


def _resolve_db_path(env_name: str, default_path: str) -> Path:
    project_root = _project_root()
    database_root = str(os.environ.get("META_DATABASE_ROOT") or "").strip()
    raw = str(os.environ.get(env_name) or "").strip()
    if not raw:
        if database_root:
            if env_name == "META_FLU_TYPE_DB":
                raw = str(Path(database_root) / "virus" / "influenza" / "type_refs.fa")
            elif env_name == "META_FLUA_HA_DB":
                raw = str(Path(database_root) / "virus" / "influenza_a" / "ha_subtypes.fa")
            elif env_name == "META_FLUA_NA_DB":
                raw = str(Path(database_root) / "virus" / "influenza_a" / "na_subtypes.fa")
        if not raw:
            if env_name == "META_FLU_TYPE_DB":
                candidates = [
                    project_root / "database" / "virus" / "influenza" / "type_refs.fa",
                    project_root / "public" / "wf-flu-master" / "data" / "primer_schemes" / "V1" / "consensus_irma.fasta",
                ]
            elif env_name == "META_FLUA_HA_DB":
                candidates = [
                    project_root / "database" / "virus" / "influenza_a" / "ha_subtypes.fa",
                ]
            elif env_name == "META_FLUA_NA_DB":
                candidates = [
                    project_root / "database" / "virus" / "influenza_a" / "na_subtypes.fa",
                ]
            else:
                candidates = []
            found = next((path for path in candidates if path.is_file()), None)
            if found is not None:
                raw = str(found)
        if not raw:
            raw = default_path
    return Path(raw).expanduser().resolve()


def _resolve_nextclade_binary() -> str:
    configured = str(os.environ.get("META_NEXTCLADE_BIN") or "").strip()
    if configured:
        return configured
    found = shutil.which("nextclade")
    if found:
        return found
    ncov_nextclade = Path(conda_env_path("ncov", "bin", "nextclade"))
    if ncov_nextclade.is_file():
        return str(ncov_nextclade)
    return "nextclade"


def _resolve_nextclade_dataset(flu_type: str) -> Path:
    normalized = str(flu_type or "").strip().upper()
    env_map = {
        "A": "META_NEXTCLADE_FLUA_DATASET",
        "B": "META_NEXTCLADE_FLUB_DATASET",
        "C": "META_NEXTCLADE_FLUC_DATASET",
        "D": "META_NEXTCLADE_FLUD_DATASET",
    }
    raw = str(os.environ.get(env_map.get(normalized, ""), "") or "").strip()
    if raw:
        return Path(raw).expanduser().resolve()
    database_root = str(os.environ.get("META_DATABASE_ROOT") or "").strip()
    if database_root:
        return (Path(database_root) / "virus" / "nextclade" / f"influenza_{normalized.lower()}").expanduser().resolve()
    return Path(f"/data/deploy/meta_genome/database/virus/nextclade/influenza_{normalized.lower()}").expanduser().resolve()


def _resolve_sars_cov_2_nextclade_dataset() -> Path:
    raw = str(os.environ.get("META_NEXTCLADE_SC2_DATASET") or "").strip()
    if raw:
        return Path(raw).expanduser().resolve()
    database_root = str(os.environ.get("META_DATABASE_ROOT") or "").strip()
    if database_root:
        return (Path(database_root) / "virus" / "nextclade" / "sars-cov-2").expanduser().resolve()
    return Path("/data/deploy/meta_genome/database/virus/nextclade/sars-cov-2").expanduser().resolve()


def _resolve_monkeypox_nextclade_dataset() -> Path:
    raw = str(os.environ.get("META_NEXTCLADE_HMPXV_DATASET") or "").strip()
    if raw:
        return Path(raw).expanduser().resolve()
    database_root = str(os.environ.get("META_DATABASE_ROOT") or "").strip()
    project_root = Path(__file__).resolve().parents[1]
    candidates: list[Path] = []
    if database_root:
        candidates.extend([
            (Path(database_root) / "virus" / "nextclade" / "hMPXV").expanduser(),
            (Path(database_root) / "nextclade_db" / "hMPXV").expanduser(),
        ])
    candidates.extend([
        project_root / "database" / "virus" / "nextclade" / "hMPXV",
        project_root / "database" / "nextclade_db" / "hMPXV",
        Path("/data/deploy/meta_genome/database/virus/nextclade/hMPXV"),
    ])
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    return candidates[0].resolve()


def _resolve_hmpv_nextclade_dataset() -> Path:
    raw = str(os.environ.get("META_NEXTCLADE_HMPV_DATASET") or "").strip()
    if raw:
        return Path(raw).expanduser().resolve()
    database_root = str(os.environ.get("META_DATABASE_ROOT") or "").strip()
    project_root = Path(__file__).resolve().parents[1]
    candidates: list[Path] = []
    if database_root:
        candidates.extend([
            (Path(database_root) / "virus" / "nextclade" / "hmpv").expanduser(),
            (Path(database_root) / "nextclade_db" / "hmpv").expanduser(),
        ])
    candidates.extend([
        project_root / "database" / "virus" / "nextclade" / "hmpv",
        project_root / "database" / "nextclade_db" / "hmpv",
        Path("/data/deploy/meta_genome/database/nextclade_db/hmpv"),
        Path("/data/deploy/meta_genome/database/virus/nextclade/hmpv"),
    ])
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    return candidates[0].resolve()


def _resolve_denv_nextclade_dataset(denv_type: str) -> Path:
    normalized = str(denv_type or "").strip()
    suffix = f"denv{normalized}" if normalized in {"1", "2", "3", "4"} else "denv1"
    env_name = f"META_NEXTCLADE_DENV{normalized}_DATASET" if normalized in {"1", "2", "3", "4"} else "META_NEXTCLADE_DENV_DATASET"
    raw = str(os.environ.get(env_name) or "").strip()
    if raw:
        return Path(raw).expanduser().resolve()
    database_root = str(os.environ.get("META_DATABASE_ROOT") or "").strip()
    project_root = Path(__file__).resolve().parents[1]
    candidates: list[Path] = []
    if database_root:
        candidates.extend([
            (Path(database_root) / "virus" / "nextclade" / suffix).expanduser(),
            (Path(database_root) / "nextclade_db" / suffix).expanduser(),
        ])
    candidates.extend([
        project_root / "database" / "virus" / "nextclade" / suffix,
        project_root / "database" / "nextclade_db" / suffix,
        Path(f"/data/deploy/meta_genome/database/nextclade_db/{suffix}"),
        Path(f"/data/deploy/meta_genome/database/virus/nextclade/{suffix}"),
    ])
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    return candidates[0].resolve()


def _resolve_zika_nextclade_dataset() -> Path:
    raw = str(os.environ.get("META_NEXTCLADE_ZIKA_DATASET") or "").strip()
    if raw:
        return Path(raw).expanduser().resolve()
    database_root = str(os.environ.get("META_DATABASE_ROOT") or "").strip()
    project_root = Path(__file__).resolve().parents[1]
    candidates: list[Path] = []
    if database_root:
        candidates.extend([
            (Path(database_root) / "virus" / "nextclade" / "zikav").expanduser(),
            (Path(database_root) / "nextclade_db" / "zikav").expanduser(),
        ])
    candidates.extend([
        project_root / "database" / "virus" / "nextclade" / "zikav",
        project_root / "database" / "nextclade_db" / "zikav",
        Path("/data/deploy/meta_genome/database/nextclade_db/zikav"),
        Path("/data/deploy/meta_genome/database/virus/nextclade/zikav"),
    ])
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    return candidates[0].resolve()


def _resolve_chikv_nextclade_dataset() -> Path:
    raw = str(os.environ.get("META_NEXTCLADE_CHIKV_DATASET") or "").strip()
    if raw:
        return Path(raw).expanduser().resolve()
    database_root = str(os.environ.get("META_DATABASE_ROOT") or "").strip()
    project_root = Path(__file__).resolve().parents[1]
    candidates: list[Path] = []
    if database_root:
        candidates.extend([
            (Path(database_root) / "virus" / "nextclade" / "chikv").expanduser(),
            (Path(database_root) / "nextclade_db" / "chikv").expanduser(),
        ])
    candidates.extend([
        project_root / "database" / "virus" / "nextclade" / "chikv",
        project_root / "database" / "nextclade_db" / "chikv",
        Path("/data/deploy/meta_genome/database/nextclade_db/chikv"),
        Path("/data/deploy/meta_genome/database/virus/nextclade/chikv"),
    ])
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    return candidates[0].resolve()


def _resolve_ebola_nextclade_dataset() -> Path:
    raw = str(os.environ.get("META_NEXTCLADE_EBOLA_DATASET") or "").strip()
    if raw:
        return Path(raw).expanduser().resolve()
    database_root = str(os.environ.get("META_DATABASE_ROOT") or "").strip()
    project_root = Path(__file__).resolve().parents[1]
    candidates: list[Path] = []
    if database_root:
        candidates.extend([
            (Path(database_root) / "virus" / "nextclade" / "ebola").expanduser(),
            (Path(database_root) / "nextclade_db" / "ebola").expanduser(),
        ])
    candidates.extend([
        project_root / "database" / "virus" / "nextclade" / "ebola",
        project_root / "database" / "nextclade_db" / "ebola",
        Path("/data/deploy/meta_genome/database/nextclade_db/ebola"),
        Path("/data/deploy/meta_genome/database/virus/nextclade/ebola"),
    ])
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    return candidates[0].resolve()


def _resolve_rsv_nextclade_dataset(rsv_type: str) -> Path:
    normalized = str(rsv_type or "").strip().upper()
    suffix = "rsv_a" if normalized == "A" else "rsv_b"
    env_name = "META_NEXTCLADE_RSVA_DATASET" if normalized == "A" else "META_NEXTCLADE_RSVB_DATASET"
    raw = str(os.environ.get(env_name) or "").strip()
    if raw:
        return Path(raw).expanduser().resolve()
    database_root = str(os.environ.get("META_DATABASE_ROOT") or "").strip()
    project_root = Path(__file__).resolve().parents[1]
    candidates: list[Path] = []
    if database_root:
        candidates.extend([
            (Path(database_root) / "virus" / "nextclade" / suffix).expanduser(),
            (Path(database_root) / "nextclade_db" / suffix).expanduser(),
        ])
    candidates.extend([
        project_root / "database" / "virus" / "nextclade" / suffix,
        project_root / "database" / "nextclade_db" / suffix,
        Path(f"/data/deploy/meta_genome/database/nextclade_db/{suffix}"),
        Path(f"/data/deploy/meta_genome/database/virus/nextclade/{suffix}"),
    ])
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    return candidates[0].resolve()


def _resolve_rsv_reference_assets(rsv_type: str) -> dict[str, Path]:
    dataset_dir = _resolve_rsv_nextclade_dataset(rsv_type)
    return {
        "type": Path(""),
        "dataset_dir": dataset_dir,
        "reference_fasta": (dataset_dir / "reference.fasta").resolve(),
        "annotation_gff": (dataset_dir / "genome_annotation.gff3").resolve(),
    }


def _resolve_denv_reference_assets(denv_type: str) -> dict[str, Path]:
    normalized = str(denv_type or "").strip()
    suffix = f"denv{normalized}" if normalized in {"1", "2", "3", "4"} else "denv1"
    env_name = f"META_NEXTCLADE_DENV{normalized}_REF_DIR" if normalized in {"1", "2", "3", "4"} else ""
    raw = str(os.environ.get(env_name) or "").strip() if env_name else ""
    if raw:
        ref_dir = Path(raw).expanduser().resolve()
    else:
        database_root = str(os.environ.get("META_DATABASE_ROOT") or "").strip()
        project_root = Path(__file__).resolve().parents[1]
        candidates: list[Path] = []
        if database_root:
            candidates.extend([
                (Path(database_root) / "virus" / "nextclade" / suffix).expanduser(),
                (Path(database_root) / "nextclade_db" / suffix).expanduser(),
            ])
        candidates.extend([
            project_root / "database" / "virus" / "nextclade" / suffix,
            project_root / "database" / "nextclade_db" / suffix,
            Path(f"/data/deploy/meta_genome/database/nextclade_db/{suffix}"),
            Path(f"/data/deploy/meta_genome/database/virus/nextclade/{suffix}"),
        ])
        ref_dir = next((candidate.resolve() for candidate in candidates if candidate.exists()), candidates[0].resolve())
    return {
        "type": Path(""),
        "reference_dir": ref_dir,
        "reference_fasta": (ref_dir / "reference.fasta").resolve(),
        "annotation_gff": (ref_dir / "genome_annotation.gff3").resolve(),
        "dataset_dir": _resolve_denv_nextclade_dataset(normalized),
    }


def _resolve_hpiv_db_dir() -> Path:
    env_root = str(os.environ.get("META_HPIV_DB_DIR") or "").strip()
    if env_root:
        return Path(env_root).expanduser().resolve()
    database_root = str(os.environ.get("META_DATABASE_ROOT") or "").strip()
    project_root = Path(__file__).resolve().parents[1]
    candidates: list[Path] = []
    if database_root:
        candidates.append((Path(database_root) / "virus" / "hpiv").expanduser())
    candidates.extend([
        project_root / "database" / "virus" / "hpiv",
        Path("/data/deploy/meta_genome/database/virus/hpiv"),
    ])
    return next((candidate.resolve() for candidate in candidates if candidate.exists()), candidates[0].resolve())


def _resolve_hpiv_reference_assets(hpiv_type: str) -> dict[str, Path]:
    normalized = str(hpiv_type or "").strip().upper()
    suffix_map = {
        "1": "HPIV1",
        "2": "HPIV2",
        "3": "HPIV3",
        "4A": "HPIV4a",
        "4B": "HPIV4b",
    }
    suffix = suffix_map.get(normalized, "HPIV1")
    root_dir = _resolve_hpiv_db_dir()
    return {
        "type": Path(""),
        "reference_dir": root_dir,
        "reference_fasta": (root_dir / f"{suffix}.fna").resolve(),
        "annotation_gff": (root_dir / f"{suffix}.gff3").resolve(),
    }


def _resolve_hpiv_subtype_db_fasta(hpiv_type: str) -> Path:
    normalized = str(hpiv_type or "").strip().upper()
    suffix_map = {
        "1": "HPIV1",
        "2": "HPIV2",
        "3": "HPIV3",
        "4A": "HPIV4a",
        "4B": "HPIV4b",
    }
    suffix = suffix_map.get(normalized, "HPIV1")
    return (_resolve_hpiv_db_dir() / f"{suffix}_db.fasta").resolve()


def _resolve_hadv_db_dir() -> Path:
    env_root = str(os.environ.get("META_HADV_DB_DIR") or "").strip()
    if env_root:
        return Path(env_root).expanduser().resolve()
    database_root = str(os.environ.get("META_DATABASE_ROOT") or "").strip()
    project_root = Path(__file__).resolve().parents[1]
    candidates: list[Path] = []
    if database_root:
        candidates.append((Path(database_root) / "virus" / "hadv").expanduser())
    candidates.extend([
        project_root / "database" / "virus" / "hadv",
        Path("/data/deploy/meta_genome/database/virus/hadv"),
    ])
    return next((candidate.resolve() for candidate in candidates if candidate.exists()), candidates[0].resolve())


def _resolve_hadv_blastn_bin() -> str:
    configured = str(os.environ.get("META_BLASTN_BIN") or "").strip()
    if configured:
        return configured
    ncov_blastn = Path(conda_env_path("ncov", "bin", "blastn"))
    if ncov_blastn.is_file():
        return str(ncov_blastn)
    return "blastn"


def _resolve_hadv_reference_fasta() -> Path:
    return (_resolve_hadv_db_dir() / "reference_genomes" / "hadv_representative_genomes_expanded.fasta").resolve()


def _resolve_hadv_full_genome_fasta() -> Path:
    return (_resolve_hadv_db_dir() / "full_genomes" / "human_mastadenovirus_A_G_complete_genomes.fasta").resolve()


def _resolve_hadv_full_genome_manifest() -> Path:
    return (_resolve_hadv_db_dir() / "full_genomes" / "human_mastadenovirus_A_G_manifest.tsv").resolve()


def _resolve_rhinovirus_db_dir() -> Path:
    env_root = str(os.environ.get("META_RHINOVIRUS_DB_DIR") or "").strip()
    if env_root:
        return Path(env_root).expanduser().resolve()
    database_root = str(os.environ.get("META_DATABASE_ROOT") or "").strip()
    project_root = Path(__file__).resolve().parents[1]
    candidates: list[Path] = []
    if database_root:
        candidates.append((Path(database_root) / "virus" / "rhinovirus").expanduser())
    candidates.extend([
        project_root / "database" / "virus" / "rhinovirus",
        Path("/data/deploy/meta_genome/database/virus/rhinovirus"),
    ])
    return next((candidate.resolve() for candidate in candidates if candidate.exists()), candidates[0].resolve())


def _resolve_rhinovirus_reference_genome_manifest() -> Path:
    return (_resolve_rhinovirus_db_dir() / "reference_genomes" / "human_rhinovirus_representative_genomes.tsv").resolve()


def _resolve_rhinovirus_reference_genome_fasta() -> Path:
    return (_resolve_rhinovirus_db_dir() / "reference_genomes" / "human_rhinovirus_representative_genomes.fasta").resolve()


def _resolve_rhinovirus_reference_vp1_manifest() -> Path:
    return (_resolve_rhinovirus_db_dir() / "reference_genomes" / "human_rhinovirus_vp1_representative_genes.tsv").resolve()


def _resolve_rhinovirus_reference_vp1_fasta() -> Path:
    return (_resolve_rhinovirus_db_dir() / "reference_genomes" / "human_rhinovirus_vp1_representative_genes.fasta").resolve()


def _resolve_rhinovirus_full_genome_manifest() -> Path:
    return (_resolve_rhinovirus_db_dir() / "full_genomes" / "human_rhinovirus_A_B_C_with_vp1_manifest.tsv").resolve()


def _resolve_rhinovirus_full_genome_fasta() -> Path:
    return (_resolve_rhinovirus_db_dir() / "full_genomes" / "human_rhinovirus_A_B_C_with_vp1.fasta").resolve()


def _load_rhinovirus_manifest(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open("r", encoding="utf-8-sig", errors="ignore", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t") if row]


def _build_rhinovirus_sequence_map(fasta_path: Path) -> dict[str, SeqRecord]:
    sequence_map: dict[str, SeqRecord] = {}
    if not fasta_path.is_file():
        return sequence_map
    for record in SeqIO.parse(str(fasta_path), "fasta"):
        sequence_map[str(record.description).strip()] = record
        sequence_map[str(record.id).strip()] = record
    return sequence_map


def _build_rhinovirus_subject_meta() -> dict[str, dict[str, str]]:
    manifest_rows = _load_rhinovirus_manifest(_resolve_rhinovirus_reference_vp1_manifest())
    subject_meta: dict[str, dict[str, str]] = {}
    for row in manifest_rows:
        record_id = str(row.get("vp1_record_id") or "").strip()
        accession = str(row.get("accession") or "").strip()
        accession_root = str(row.get("accession_root") or accession.split(".", 1)[0]).strip()
        normalized_type = _normalize_rhinovirus_type_from_row(row)
        meta = {
            "species_group": str(row.get("species_group") or "").strip().upper(),
            "normalized_type": normalized_type,
            "type_label": str(row.get("type_label") or "").strip(),
            "accession": accession,
            "accession_root": accession_root,
            "gff_path": str(row.get("gff_path") or "").strip(),
        }
        for key in {record_id, accession, accession_root}:
            if key:
                subject_meta[key] = meta
    return subject_meta


def _score_rhinovirus_type_rows(rows: list[dict[str, object]]) -> dict[str, str]:
    by_type: dict[str, dict[str, object]] = {}
    for row in rows:
        meta = row.get("_meta") if isinstance(row.get("_meta"), dict) else {}
        normalized_type = _normalize_rhinovirus_type(str(meta.get("normalized_type") or ""))
        if not normalized_type:
            continue
        score = (
            float(row.get("coverage") or 0.0),
            float(row.get("mean_depth") or 0.0),
            float(row.get("covered_bases") or 0.0),
            float(row.get("num_reads") or 0.0),
        )
        bucket = by_type.setdefault(
            normalized_type,
            {
                "normalized_type": normalized_type,
                "species_group": str(meta.get("species_group") or ""),
                "subject": str(row.get("reference_name") or ""),
                "coverage": score[0],
                "mean_depth": score[1],
                "covered_bases": score[2],
                "num_reads": score[3],
                "reference_count": set(),
            },
        )
        if score > (
            float(bucket["coverage"]),
            float(bucket["mean_depth"]),
            float(bucket["covered_bases"]),
            float(bucket["num_reads"]),
        ):
            bucket["subject"] = str(row.get("reference_name") or "")
            bucket["coverage"] = score[0]
            bucket["mean_depth"] = score[1]
            bucket["covered_bases"] = score[2]
            bucket["num_reads"] = score[3]
        bucket["reference_count"].add(str(meta.get("accession_root") or str(row.get("reference_name") or "").split("_", 1)[0]))
    if not by_type:
        return {"type": "", "species_group": "", "subject": "", "coverage": "", "mean_depth": "", "covered_bases": "", "num_reads": "", "reference_count": "0"}
    ranked = sorted(
        by_type.values(),
        key=lambda item: (
            float(item["coverage"]),
            float(item["mean_depth"]),
            float(item["covered_bases"]),
            float(item["num_reads"]),
            len(item["reference_count"]),
        ),
        reverse=True,
    )
    best = ranked[0]
    return {
        "type": str(best["normalized_type"]),
        "species_group": str(best["species_group"]),
        "subject": str(best["subject"]),
        "coverage": f"{float(best['coverage']):.2f}",
        "mean_depth": f"{float(best['mean_depth']):.2f}",
        "covered_bases": f"{float(best['covered_bases']):.0f}",
        "num_reads": f"{float(best['num_reads']):.0f}",
        "reference_count": str(len(best["reference_count"])),
    }


def _run_rhinovirus_vp1_read_typing(
    db_fasta: Path,
    out_dir: Path,
    single_fastq: str = "",
    fq1: str = "",
    fq2: str = "",
    long_type: str = "",
    threads: int = 4,
    logf=None,
) -> dict[str, str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    bam_path = out_dir / "vp1.screening.bam"
    coverage_path = out_dir / "vp1.coverage.tsv"
    if not bam_path.is_file():
        if fq1:
            map_cmd = " ".join(
                [
                    "minimap2",
                    "-ax",
                    "sr",
                    shlex.quote(str(db_fasta)),
                    shlex.quote(str(fq1)),
                    *([shlex.quote(str(fq2))] if fq2 else []),
                    "|",
                    "samtools",
                    "sort",
                    "-o",
                    shlex.quote(str(bam_path)),
                ]
            )
        else:
            map_cmd = " ".join(
                [
                    "minimap2",
                    "-ax",
                    _choose_minimap2_preset(long_type),
                    shlex.quote(str(db_fasta)),
                    shlex.quote(str(single_fastq)),
                    "-t",
                    str(max(1, int(threads or 1))),
                    "|",
                    "samtools",
                    "sort",
                    "-o",
                    shlex.quote(str(bam_path)),
                ]
            )
        run_command(map_cmd, logf=logf)
    if not bam_path.is_file():
        return {"type": "", "species_group": "", "subject": "", "identity": "", "coverage": "", "mean_depth": "", "covered_bases": "", "num_reads": "", "reference_count": "0", "method": "read_coverage"}
    run_command(f"samtools index -@ {max(1, int(threads or 1))} {shlex.quote(str(bam_path))}", logf=logf)
    run_command(f"samtools coverage {shlex.quote(str(bam_path))} > {shlex.quote(str(coverage_path))}", logf=logf)
    subject_meta = _build_rhinovirus_subject_meta()
    rows = _parse_samtools_coverage_rows(coverage_path)
    enriched_rows = []
    for row in rows:
        subject = str(row.get("reference_name") or "").strip()
        meta = subject_meta.get(subject) or subject_meta.get(subject.split("_", 1)[0]) or {}
        enriched = dict(row)
        enriched["_meta"] = meta
        enriched_rows.append(enriched)
    best = _score_rhinovirus_type_rows(enriched_rows)
    return {
        "type": best["type"],
        "species_group": best["species_group"],
        "subject": best["subject"],
        "identity": "",
        "coverage": best["coverage"],
        "mean_depth": best["mean_depth"],
        "covered_bases": best["covered_bases"],
        "num_reads": best["num_reads"],
        "reference_count": best["reference_count"],
        "method": "read_coverage",
    }


def _run_rhinovirus_vp1_blast_typing(query_fasta: Path, db_fasta: Path, out_dir: Path, logf=None) -> dict[str, str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    out_tsv = out_dir / "vp1.blastn.tsv"
    blastn_bin = _resolve_hadv_blastn_bin()
    cmd = " ".join(
        [
            shlex.quote(blastn_bin),
            "-task blastn",
            "-query",
            shlex.quote(str(query_fasta)),
            "-subject",
            shlex.quote(str(db_fasta)),
            "-outfmt",
            shlex.quote("6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"),
            "-evalue 1e-20",
            "-max_target_seqs 50",
            "-dust no",
            "-num_threads 4",
            "-out",
            shlex.quote(str(out_tsv)),
        ]
    )
    run_command(cmd, logf=logf)
    if not out_tsv.is_file() or out_tsv.stat().st_size == 0:
        return {"type": "", "species_group": "", "subject": "", "identity": "", "coverage": "", "mean_depth": "", "covered_bases": "", "num_reads": "", "reference_count": "0", "method": "blastn"}

    subject_meta = _build_rhinovirus_subject_meta()
    ref_lengths: dict[str, int] = {}
    with db_fasta.open("r", encoding="utf-8", errors="ignore") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            ref_lengths[str(record.id)] = len(str(record.seq))

    by_type: dict[str, dict[str, object]] = {}
    with out_tsv.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            _qseqid, sseqid, pident, length, _mismatch, _gapopen, _qstart, _qend, _sstart, _send, _evalue, bitscore = parts[:12]
            meta = subject_meta.get(sseqid) or subject_meta.get(sseqid.split("_", 1)[0]) or {}
            normalized_type = _normalize_rhinovirus_type(str(meta.get("normalized_type") or ""))
            species_group = str(meta.get("species_group") or "").strip().upper()
            if not normalized_type:
                continue
            ref_len = ref_lengths.get(sseqid, 0)
            try:
                qcov_ref = min(100.0, (float(length) / ref_len) * 100) if ref_len else 0.0
                score = (float(bitscore), float(pident), qcov_ref, float(length))
            except ValueError:
                continue
            bucket = by_type.setdefault(
                normalized_type,
                {
                    "type": normalized_type,
                    "species_group": species_group,
                    "subject": sseqid,
                    "identity": str(pident),
                    "coverage": f"{qcov_ref:.2f}",
                    "covered_bases": str(length),
                    "reference_count": set(),
                    "_score": score,
                },
            )
            if score > bucket["_score"]:
                bucket["subject"] = sseqid
                bucket["identity"] = str(pident)
                bucket["coverage"] = f"{qcov_ref:.2f}"
                bucket["covered_bases"] = str(length)
                bucket["_score"] = score
                bucket["species_group"] = species_group
            bucket["reference_count"].add(str(meta.get("accession_root") or sseqid.split("_", 1)[0]))
    if not by_type:
        return {"type": "", "species_group": "", "subject": "", "identity": "", "coverage": "", "mean_depth": "", "covered_bases": "", "num_reads": "", "reference_count": "0", "method": "blastn"}
    ranked = sorted(
        by_type.values(),
        key=lambda item: (
            item["_score"][0],
            item["_score"][1],
            item["_score"][2],
            len(item["reference_count"]),
        ),
        reverse=True,
    )
    best = ranked[0]
    return {
        "type": str(best["type"]),
        "species_group": str(best["species_group"]),
        "subject": str(best["subject"]),
        "identity": str(best["identity"]),
        "coverage": str(best["coverage"]),
        "mean_depth": "",
        "covered_bases": str(best["covered_bases"]),
        "num_reads": "",
        "reference_count": str(len(best["reference_count"])),
        "method": "blastn",
    }


def _run_rhinovirus_vp1_typing(
    pre: str,
    query_fasta: Path | None = None,
    single_fastq: str = "",
    fq1: str = "",
    fq2: str = "",
    long_type: str = "",
    threads: int = 4,
    output_dir: Path | None = None,
    logf=None,
) -> dict[str, object]:
    typing_dir = (output_dir or (Path(f"{pre}_rhinovirus_reference_selection") / "typing")).resolve()
    typing_dir.mkdir(parents=True, exist_ok=True)
    db_fasta = _resolve_rhinovirus_reference_vp1_fasta()
    if not db_fasta.is_file():
        return {"vp1_type": "", "species_group": "", "summary_path": "", "hit": {"method": "missing_db"}}
    use_blast = query_fasta is not None and query_fasta.is_file() and query_fasta.stat().st_size > 0
    hit = (
        _run_rhinovirus_vp1_blast_typing(query_fasta, db_fasta, typing_dir, logf=logf)
        if use_blast
        else _run_rhinovirus_vp1_read_typing(
            db_fasta,
            typing_dir,
            single_fastq=single_fastq,
            fq1=fq1,
            fq2=fq2,
            long_type=long_type,
            threads=threads,
            logf=logf,
        )
    )
    hit_status = "ready"
    hit_note = ""
    if _is_low_support_gene_typing_hit(hit):
        hit_status = "low_support"
        hit_note = _low_support_gene_typing_note("VP1", hit)
        hit = dict(hit)
        hit["type"] = ""
        hit["species_group"] = ""
    summary_path = typing_dir / "vp1_typing.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "gene", "matched_type", "species_group", "subject", "identity", "coverage", "mean_depth", "covered_bases", "num_reads", "reference_count", "method", "status", "note"])
        writer.writerow([
            pre,
            "vp1",
            hit.get("type", ""),
            hit.get("species_group", ""),
            hit.get("subject", ""),
            hit.get("identity", ""),
            hit.get("coverage", ""),
            hit.get("mean_depth", ""),
            hit.get("covered_bases", ""),
            hit.get("num_reads", ""),
            hit.get("reference_count", "0"),
            hit.get("method", ""),
            hit_status,
            hit_note,
        ])
    return {
        "vp1_type": str(hit.get("type") or ""),
        "species_group": str(hit.get("species_group") or ""),
        "summary_path": str(summary_path.resolve()),
        "hit": hit,
        "status": hit_status,
        "note": hit_note,
    }


def _resolve_rhinovirus_vadr_model_dir(project_root: Path, species_group: str) -> Path | None:
    normalized = str(species_group or "").strip().upper()
    if normalized not in {"A", "B", "C"}:
        return None
    model_dir = (project_root / "soft" / "vadr-models-hrv" / f"hrv{normalized}").resolve()
    if model_dir.is_dir():
        return model_dir
    return None


def _run_vadr_rhinovirus_annotation(pre: str, input_fasta: Path, output_root: Path, species_group: str, logf=None) -> dict[str, str]:
    project_root = _project_root()
    normalized = str(species_group or "").strip().upper()
    model_dir = _resolve_rhinovirus_vadr_model_dir(project_root, normalized)
    if model_dir is None:
        return {"status": "missing", "note": f"未找到 Rhinovirus {normalized or '?'} VADR 模型目录", "output_dir": "", "gff_path": ""}
    env = _build_vadr_env(project_root, model_dir)
    if env is None:
        return {"status": "missing", "note": "VADR 运行依赖不完整", "output_dir": "", "gff_path": ""}
    rhinovirus_env = dict(env)
    alt_bio_easel_dir = project_root / "soft" / "Bio-Easel-ncov"
    alt_bio_lib = alt_bio_easel_dir / "blib" / "lib"
    alt_bio_arch = alt_bio_easel_dir / "blib" / "arch"
    if alt_bio_lib.is_dir() and alt_bio_arch.is_dir():
        rhinovirus_env["VADRBIOEASELDIR"] = str(alt_bio_easel_dir)
        rhinovirus_env["PERL5LIB"] = os.pathsep.join(
            [
                str(project_root / "soft" / "vadr"),
                str(project_root / "soft" / "sequip"),
                str(alt_bio_lib),
                str(alt_bio_arch),
            ]
        )
    ncov_bin_dir = conda_env_path("ncov", "bin")
    rhinovirus_env["PATH"] = os.pathsep.join(
        [
            ncov_bin_dir,
            "/usr/bin",
            "/bin",
            str(rhinovirus_env.get("PATH") or ""),
        ]
    )
    ncov_perl = f"{ncov_bin_dir}/perl"
    perl_candidates = [ncov_perl, "/usr/bin/perl", _resolve_vadr_perl_bin(), "perl"]
    perl_bin = _resolve_working_perl_with_module(rhinovirus_env, "Bio::Easel::MSA", perl_candidates)
    vadr_root = output_root / "vadr"
    vadr_root.mkdir(parents=True, exist_ok=True)
    output_dir = vadr_root / f"{pre}_vadr"
    prefix = f"{output_dir.name}.vadr"
    pass_tbl = output_dir / f"{prefix}.pass.tbl"
    fail_tbl = output_dir / f"{prefix}.fail.tbl"
    gff_path = vadr_root / f"{pre}.vadr.gff3"
    if output_dir.is_dir() and (pass_tbl.is_file() or fail_tbl.is_file()) and gff_path.is_file() and gff_path.stat().st_size > 0:
        return {"status": "ready", "note": "已检测到现有鼻病毒 VADR 注释结果", "output_dir": str(output_dir.resolve()), "gff_path": str(gff_path.resolve())}
    vadr_script = project_root / "soft" / "vadr" / "v-annotate.pl"
    cmd = " ".join(
        [
            shlex.quote(str(perl_bin)),
            shlex.quote(str(vadr_script)),
            "-f",
            "-r",
            "--ignore_exc",
            "--mkey",
            f"hrv{normalized}",
            "--mdir",
            shlex.quote(str(model_dir)),
            shlex.quote(str(input_fasta)),
            shlex.quote(str(output_dir)),
        ]
    )
    try:
            run_command(cmd, logf=logf, env=rhinovirus_env)
    except Exception as exc:
        return {"status": "failed", "note": f"鼻病毒 VADR 注释失败: {exc}", "output_dir": str(output_dir.resolve()) if output_dir.exists() else "", "gff_path": ""}
    annotate_tbl2gff = project_root / "soft" / "vadr" / "miniscripts" / "annotate-tbl2gff.pl"
    source_tbl = pass_tbl if pass_tbl.is_file() else fail_tbl
    if source_tbl.is_file():
        try:
            run_command(
                f"{shlex.quote(str(perl_bin))} {shlex.quote(str(annotate_tbl2gff))} {shlex.quote(str(source_tbl))} > {shlex.quote(str(gff_path))}",
                logf=logf,
                env=rhinovirus_env,
            )
        except Exception:
            pass
    return {
        "status": "ready" if gff_path.is_file() and gff_path.stat().st_size > 0 else "failed",
        "note": "鼻病毒 VADR 注释完成" if gff_path.is_file() and gff_path.stat().st_size > 0 else "鼻病毒 VADR 未生成可用 GFF",
        "output_dir": str(output_dir.resolve()),
        "gff_path": str(gff_path.resolve()) if gff_path.is_file() and gff_path.stat().st_size > 0 else "",
    }


def prepare_rhinovirus_sample_annotation(pre: str, sample_fasta: Path, output_root: Path, species_group: str, logf=None) -> Path | None:
    result = _run_vadr_rhinovirus_annotation(pre, sample_fasta, output_root, species_group, logf=logf)
    gff_path = str(result.get("gff_path") or "").strip()
    if result.get("status") == "ready" and gff_path:
        gff_file = Path(gff_path)
        if gff_file.is_file() and gff_file.stat().st_size > 0:
            return gff_file
    return None


def prepare_rhinovirus_reference_annotation(pre: str, reference_fasta: Path, output_root: Path, species_group: str, logf=None) -> Path | None:
    result = _run_vadr_rhinovirus_annotation(pre, reference_fasta, output_root, species_group, logf=logf)
    gff_path = str(result.get("gff_path") or "").strip()
    if result.get("status") == "ready" and gff_path:
        gff_file = Path(gff_path)
        if gff_file.is_file() and gff_file.stat().st_size > 0:
            return gff_file
    return None


def _find_rhinovirus_vp1_feature(gff_path: Path) -> dict[str, object] | None:
    if not gff_path.is_file():
        return None
    candidates: list[tuple[int, dict[str, object]]] = []
    with gff_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, _source, feature_type, start_text, end_text, _score, strand, _phase, attrs_text = parts[:9]
            attrs = _parse_norovirus_gff_attributes(attrs_text)
            product = str(attrs.get("product") or "").strip().upper()
            note = str(attrs.get("Note") or attrs.get("note") or "").strip().upper()
            label = " ".join([product, note, str(attrs.get("gene") or "").strip().upper()])
            if feature_type not in {"mature_protein_region_of_CDS", "mat_peptide"}:
                continue
            if "VP1" not in label and product not in {"1D", "1D (VP1)"}:
                continue
            try:
                start = int(start_text)
                end = int(end_text)
            except ValueError:
                continue
            priority = 0 if feature_type == "mat_peptide" else 1
            candidates.append(
                (
                    priority,
                    {
                        "seqid": seqid,
                        "feature_type": feature_type,
                        "start": start,
                        "end": end,
                        "strand": strand,
                        "label": product or "VP1",
                    },
                )
            )
    if not candidates:
        return None
    candidates.sort(key=lambda item: (item[0], int(item[1]["start"]), int(item[1]["end"])))
    return candidates[0][1]


def _extract_rhinovirus_vp1_region(sample_fasta: Path, gff_path: Path, sample_label: str) -> tuple[SeqRecord | None, dict[str, object]]:
    feature = _find_rhinovirus_vp1_feature(gff_path)
    if feature is None:
        return None, {"status": "missing", "note": "GFF 中未找到 VP1 区域"}
    records = list(SeqIO.parse(str(sample_fasta), "fasta"))
    if not records:
        return None, {"status": "missing", "note": f"{sample_fasta.name} 中未找到序列记录"}
    seqid = str(feature.get("seqid") or "").strip()
    record = next((item for item in records if str(item.id).strip() == seqid), None)
    if record is None:
        record = max(records, key=lambda item: len(str(item.seq)))
    start = int(feature["start"])
    end = int(feature["end"])
    sequence = str(record.seq)
    if start < 1 or end > len(sequence) or start > end:
        return None, {"status": "missing", "note": "VP1 坐标超出样本序列长度"}
    fragment = sequence[start - 1:end]
    if str(feature.get("strand") or "+") == "-":
        fragment = str(Seq(fragment).reverse_complement())
    fragment = fragment.strip().strip("Nn")
    if not fragment:
        return None, {"status": "missing", "note": "VP1 区域仅包含空白或 N"}
    record_id = f"{sample_label}_vp1_sample"
    return SeqRecord(Seq(fragment), id=record_id, description=f"{sample_label} VP1 sample"), {
        "status": "ready",
        "start": start,
        "end": end,
        "strand": str(feature.get("strand") or "+"),
        "seqid": str(record.id),
        "length": len(fragment),
    }


def _filter_rhinovirus_vp1_refs_by_group(species_group: str) -> tuple[list[dict[str, str]], list[SeqRecord]]:
    manifest_rows = _load_rhinovirus_manifest(_resolve_rhinovirus_reference_vp1_manifest())
    fasta_path = _resolve_rhinovirus_reference_vp1_fasta()
    if not fasta_path.is_file():
        return [], []
    target_group = str(species_group or "").strip().upper()
    selected_rows = [row for row in manifest_rows if str(row.get("species_group") or "").strip().upper() == target_group]
    if not selected_rows:
        return [], []
    sequence_map = _build_rhinovirus_sequence_map(fasta_path)
    selected_records: list[SeqRecord] = []
    ordered_rows: list[dict[str, str]] = []
    for row in selected_rows:
        record_id = str(row.get("vp1_record_id") or "").strip()
        template = sequence_map.get(record_id)
        if template is None:
            continue
        ordered_rows.append(
            {
                "gene": "VP1",
                "subtype": str(row.get("normalized_type") or row.get("type_label") or "").strip(),
                "accession": str(row.get("accession_root") or row.get("accession") or "").split(".", 1)[0],
                "backup_rank": "1",
            }
        )
        selected_records.append(
            SeqRecord(
                template.seq,
                id=record_id,
                name=record_id,
                description=str(template.description),
            )
        )
    return ordered_rows, selected_records


def build_rhinovirus_vp1_phylogeny_assets(pre: str, sample_fasta: Path, gff_path: Path, species_group: str = "", logf=None) -> dict[str, object]:
    output_root = Path(f"{pre}_rhinovirus_reference_selection") / "phylogeny"
    output_root.mkdir(parents=True, exist_ok=True)
    selection_path = Path(f"{pre}_rhinovirus_reference_selection") / "selection.tsv"
    vp1_type = ""
    selected_group = str(species_group or "").strip().upper()
    if selection_path.is_file():
        with selection_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
            row = next(csv.DictReader(handle, delimiter="\t"), None) or {}
        vp1_type = str(row.get("vp1_type") or "").strip()
        selected_group = selected_group or str(row.get("species_group") or "").strip().upper()
    result_row = {
        "gene": "VP1",
        "subtype": vp1_type,
        "genogroup": selected_group,
        "status": "missing",
        "note": "",
        "sample_region_path": "",
        "backup_fasta": "",
        "tree_path": "",
        "member_count": "0",
    }
    sample_record, sample_meta = _extract_rhinovirus_vp1_region(sample_fasta, gff_path, _sanitize_tree_label(pre))
    if sample_record is None:
        result_row["note"] = str(sample_meta.get("note") or "未提取到 VP1 区域")
    else:
        backup_rows, backup_records = _filter_rhinovirus_vp1_refs_by_group(selected_group)
        if not backup_records:
            result_row["note"] = f"VP1 参考库中未找到 Rhinovirus {selected_group or '?'} 记录"
        else:
            tree_result = _build_norovirus_gene_tree("vp1", sample_record, backup_rows, backup_records, output_root / "vp1", logf=logf)
            result_row["status"] = str(tree_result.get("status") or "missing")
            result_row["note"] = str(tree_result.get("note") or "")
            result_row["sample_region_path"] = str((output_root / "vp1" / "vp1.sample.fasta").resolve())
            result_row["backup_fasta"] = str(tree_result.get("backup_fasta") or "")
            result_row["tree_path"] = str(tree_result.get("tree_path") or "")
            result_row["member_count"] = str(tree_result.get("member_count") or "0")
    summary_path = output_root / "summary.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["gene", "subtype", "genogroup", "status", "note", "sample_region_path", "backup_fasta", "tree_path", "member_count"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(result_row)
    return {
        "status": "ready" if result_row.get("status") == "ready" else "missing",
        "summary_path": str(summary_path.resolve()),
        "rows": [result_row],
    }


def resolve_rhinovirus_reference(
    pre: str,
    species: str = "",
    requested_ref: str = "",
    single_fastq: str = "",
    fq1: str = "",
    fq2: str = "",
    long_type: str = "",
    threads: int = 4,
    query_fasta: str = "",
    logf=None,
) -> dict[str, str]:
    _clear_sample_skip_flag(pre)
    requested = str(requested_ref or "").strip()
    if requested and Path(requested).is_file():
        group = _extract_rhinovirus_species_group(species)
        return {
            "status": "ready",
            "vp1_type": _normalize_rhinovirus_type(species),
            "species_group": group,
            "species_label": _rhinovirus_species_label(group),
            "reference_path": str(Path(requested).resolve()),
            "gff_path": "nogtf",
            "summary_path": "",
            "typing_summary_path": "",
        }
    screening_dir = Path(f"{pre}_rhinovirus_reference_selection")
    screening_dir.mkdir(exist_ok=True)
    query_path = Path(query_fasta) if str(query_fasta or "").strip() else None
    typing_result = _run_rhinovirus_vp1_typing(
        pre,
        query_fasta=query_path,
        single_fastq=single_fastq,
        fq1=fq1,
        fq2=fq2,
        long_type=long_type,
        threads=threads,
        output_dir=screening_dir / "typing",
        logf=logf,
    )
    vp1_type = str(typing_result.get("vp1_type") or "").strip()
    species_group = str(typing_result.get("species_group") or "").strip().upper()
    summary_path = screening_dir / "selection.tsv"
    if not vp1_type or species_group not in {"A", "B", "C"}:
        status_text = "skipped" if str(typing_result.get("status") or "").strip() == "low_support" else "missing"
        note_text = str(typing_result.get("note") or "").strip() or f"未能根据 VP1 或 species={species or '-'} 确定鼻病毒分型"
        with summary_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sample", "vp1_type", "species_group", "status", "note", "typing_summary_path"])
            writer.writerow([pre, vp1_type, species_group, status_text, note_text, str(typing_result.get("summary_path") or "")])
        if status_text == "skipped":
            _write_sample_skip_flag(pre, note_text)
        return {
            "status": status_text,
            "vp1_type": vp1_type,
            "species_group": species_group,
            "species_label": species or "Human rhinovirus",
            "reference_path": "",
            "gff_path": "nogtf",
            "summary_path": str(summary_path.resolve()),
            "typing_summary_path": str(typing_result.get("summary_path") or ""),
        }
    manifest_rows = _load_rhinovirus_manifest(_resolve_rhinovirus_full_genome_manifest())
    fasta_path = _resolve_rhinovirus_full_genome_fasta()
    sequence_map = _build_rhinovirus_sequence_map(fasta_path)
    candidates = [
        row for row in manifest_rows
        if str(row.get("species_group") or "").strip().upper() == species_group
        and _normalize_rhinovirus_type_from_row(row) == vp1_type
    ]
    candidate_fasta = screening_dir / "candidate_references.fasta"
    with candidate_fasta.open("w", encoding="utf-8") as handle:
        for row in candidates:
            header = str(row.get("header") or row.get("source_header") or row.get("accession") or "").strip()
            record = sequence_map.get(header) or sequence_map.get(str(row.get("accession") or "").strip())
            if record is None:
                continue
            fasta_id = _sanitize_tree_label(str(row.get("accession_root") or row.get("accession") or record.id))
            handle.write(f">{fasta_id}\n")
            seq = str(record.seq)
            for index in range(0, len(seq), 80):
                handle.write(seq[index:index + 80] + "\n")
    if not candidate_fasta.is_file() or candidate_fasta.stat().st_size == 0:
        with summary_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sample", "vp1_type", "species_group", "status", "note", "typing_summary_path"])
            writer.writerow([pre, vp1_type, species_group, "missing", "未找到对应亚型的全基因组候选参考", str(typing_result.get("summary_path") or "")])
        return {
            "status": "missing",
            "vp1_type": vp1_type,
            "species_group": species_group,
            "species_label": _rhinovirus_species_label(species_group),
            "reference_path": "",
            "gff_path": "nogtf",
            "summary_path": str(summary_path.resolve()),
            "typing_summary_path": str(typing_result.get("summary_path") or ""),
        }
    bam_path = screening_dir / "candidate_references.bam"
    coverage_path = screening_dir / "candidate_references.coverage.tsv"
    if fq1:
        map_cmd = " ".join(
            [
                "minimap2",
                "-ax",
                "sr",
                shlex.quote(str(candidate_fasta)),
                shlex.quote(str(fq1)),
                *([shlex.quote(str(fq2))] if fq2 else []),
                "|",
                "samtools",
                "sort",
                "-o",
                shlex.quote(str(bam_path)),
            ]
        )
    else:
        map_cmd = " ".join(
            [
                "minimap2",
                "-ax",
                _choose_minimap2_preset(long_type),
                shlex.quote(str(candidate_fasta)),
                shlex.quote(str(single_fastq)),
                "-t",
                str(max(1, int(threads or 1))),
                "|",
                "samtools",
                "sort",
                "-o",
                shlex.quote(str(bam_path)),
            ]
        )
    run_command(map_cmd, logf=logf)
    run_command(f"samtools index -@ {max(1, int(threads or 1))} {shlex.quote(str(bam_path))}", logf=logf)
    run_command(f"samtools coverage {shlex.quote(str(bam_path))} > {shlex.quote(str(coverage_path))}", logf=logf)
    coverage_rows = _parse_samtools_coverage_rows(coverage_path)
    candidate_meta = { _sanitize_tree_label(str(row.get('accession_root') or row.get('accession') or '')): row for row in candidates }
    best = None
    for row in coverage_rows:
        ref_name = str(row.get("reference_name") or "").strip()
        meta = candidate_meta.get(ref_name)
        if meta is None:
            continue
        score = (
            float(row.get("coverage") or 0.0),
            float(row.get("mean_depth") or 0.0),
            float(row.get("covered_bases") or 0.0),
            float(row.get("num_reads") or 0.0),
        )
        if best is None or score > best["score"]:
            best = {"score": score, "meta": meta, "coverage_row": row}
    if best is None:
        return {
            "status": "missing",
            "vp1_type": vp1_type,
            "species_group": species_group,
            "species_label": _rhinovirus_species_label(species_group),
            "reference_path": "",
            "gff_path": "nogtf",
            "summary_path": "",
            "typing_summary_path": str(typing_result.get("summary_path") or ""),
        }
    best_meta = best["meta"]
    best_row = best["coverage_row"]
    best_record = sequence_map.get(str(best_meta.get("header") or best_meta.get("source_header") or "").strip()) or sequence_map.get(str(best_meta.get("accession") or "").strip())
    best_reference_fasta = screening_dir / f"{_sanitize_tree_label(str(best_meta.get('normalized_type') or vp1_type))}.reference.fasta"
    with best_reference_fasta.open("w", encoding="utf-8") as handle:
        handle.write(f">{str(best_meta.get('accession') or best_record.id if best_record else 'reference')}\n")
        sequence = str(best_record.seq if best_record is not None else "")
        for index in range(0, len(sequence), 80):
            handle.write(sequence[index:index + 80] + "\n")
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "vp1_type", "species_group", "coverage", "mean_depth", "covered_bases", "num_reads", "reference_name", "accession", "reference_path", "gff_path", "typing_summary_path"])
        writer.writerow([
            pre,
            vp1_type,
            species_group,
            f"{float(best_row.get('coverage') or 0.0):.6f}",
            f"{float(best_row.get('mean_depth') or 0.0):.6f}",
            f"{float(best_row.get('covered_bases') or 0.0):.0f}",
            f"{float(best_row.get('num_reads') or 0.0):.0f}",
            str(best_meta.get("source_title") or best_meta.get("header") or best_meta.get("accession") or ""),
            str(best_meta.get("accession") or ""),
            str(best_reference_fasta.resolve()),
            str(best_meta.get("gff_path") or "nogtf"),
            str(typing_result.get("summary_path") or ""),
        ])
    return {
        "status": "ready",
        "vp1_type": vp1_type,
        "species_group": species_group,
        "species_label": _rhinovirus_species_label(species_group),
        "reference_path": str(best_reference_fasta.resolve()),
        "gff_path": str(best_meta.get("gff_path") or "nogtf"),
        "summary_path": str(summary_path.resolve()),
        "typing_summary_path": str(typing_result.get("summary_path") or ""),
    }


def run_rhinovirus_consensus_typing(pre: str, consensus_fasta: Path, logf=None) -> dict[str, str]:
    if not consensus_fasta.is_file() or consensus_fasta.stat().st_size == 0:
        return {"status": "missing", "summary_path": "", "vp1_type": "", "species_group": ""}
    output_dir = Path(f"{pre}_rhinovirus_reference_selection") / "consensus_typing"
    result = _run_rhinovirus_vp1_typing(pre, query_fasta=consensus_fasta, output_dir=output_dir, logf=logf)
    summary_path = output_dir / "consensus_typing.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "vp1_type", "species_group", "typing_summary_path"])
        writer.writerow([pre, result.get("vp1_type", ""), result.get("species_group", ""), result.get("summary_path", "")])
    return {
        "status": "ready",
        "summary_path": str(summary_path.resolve()),
        "vp1_type": str(result.get("vp1_type") or ""),
        "species_group": str(result.get("species_group") or ""),
    }


def _resolve_enterovirus_db_dir() -> Path:
    env_root = str(os.environ.get("META_ENTEROVIRUS_DB_DIR") or "").strip()
    if env_root:
        return Path(env_root).expanduser().resolve()
    database_root = str(os.environ.get("META_DATABASE_ROOT") or "").strip()
    project_root = _project_root()
    candidates: list[Path] = []
    if database_root:
        candidates.append((Path(database_root) / "virus" / "enterovirus").expanduser())
    candidates.extend(
        [
            project_root / "database" / "virus" / "enterovirus",
            Path("/data/deploy/meta_genome/database/virus/enterovirus"),
        ]
    )
    for candidate in candidates:
        if candidate.is_dir():
            return candidate.resolve()
    return candidates[0].resolve()


def _resolve_enterovirus_reference_vp1_manifest() -> Path:
    return (_resolve_enterovirus_db_dir() / "reference_genomes" / "abcd_vp1" / "enterovirus_abcd_vp1.tsv").resolve()


def _resolve_enterovirus_reference_vp1_fasta() -> Path:
    return (_resolve_enterovirus_db_dir() / "reference_genomes" / "abcd_vp1" / "enterovirus_abcd_vp1.fasta").resolve()


def _resolve_enterovirus_expanded_complete_manifest() -> Path:
    return (_resolve_enterovirus_db_dir() / "reference_genomes" / "abcd_vp1" / "enterovirus_abcd_complete_genomes_expanded_manifest.tsv").resolve()


def _resolve_enterovirus_ev_typing_dir() -> Path:
    return (_resolve_enterovirus_db_dir() / "reference_genomes" / "ev_typing_out").resolve()


def _resolve_enterovirus_ev_typing_summary() -> Path:
    return (_resolve_enterovirus_ev_typing_dir() / "typing_summary.tsv").resolve()


def _resolve_enterovirus_ev_db_fasta() -> Path:
    return (_resolve_enterovirus_db_dir() / "reference_genomes" / "ev_db.fasta").resolve()


def _is_enterovirus_abcd_reference_row(row: dict[str, object]) -> bool:
    big_group = str(row.get("big_group") or "").strip().upper()
    abbrev = str(row.get("abbrev") or "").strip().upper()
    fasta_path = Path(str(row.get("fasta_path") or "").strip())
    vp1_fasta_path = Path(str(row.get("vp1_fasta_path") or "").strip())
    if big_group not in {"A", "B", "C", "D"}:
        return False
    if not abbrev or abbrev.startswith("RV-"):
        return False
    return fasta_path.is_file() and vp1_fasta_path.is_file()


def _load_enterovirus_abcd_manifest() -> list[dict[str, str]]:
    rows = _load_rhinovirus_manifest(_resolve_enterovirus_reference_vp1_manifest())
    return [row for row in rows if _is_enterovirus_abcd_reference_row(row)]


def _is_enterovirus_candidate_reference_row(row: dict[str, object]) -> bool:
    big_group = str(row.get("big_group") or "").strip().upper()
    abbrev = str(row.get("abbrev") or "").strip().upper()
    fasta_path = Path(str(row.get("fasta_path") or "").strip())
    if big_group not in {"A", "B", "C", "D"}:
        return False
    if not abbrev or abbrev.startswith("RV-"):
        return False
    return fasta_path.is_file()


def _load_enterovirus_candidate_reference_manifest() -> list[dict[str, str]]:
    combined_by_accession: dict[str, dict[str, str]] = {}
    for row in _load_enterovirus_abcd_manifest():
        accession = str(row.get("accession") or row.get("accession_full") or "").split(".", 1)[0]
        if accession:
            combined_by_accession[accession] = dict(row)

    expanded_manifest = _resolve_enterovirus_expanded_complete_manifest()
    if expanded_manifest.is_file():
        for row in _load_rhinovirus_manifest(expanded_manifest):
            if not _is_enterovirus_candidate_reference_row(row):
                continue
            accession = str(row.get("accession") or row.get("accession_full") or "").split(".", 1)[0]
            if not accession:
                continue
            merged = dict(combined_by_accession.get(accession, {}))
            merged.update(row)
            merged.setdefault("accession", accession)
            combined_by_accession[accession] = merged

    ev_typing_summary = _resolve_enterovirus_ev_typing_summary()
    ev_db_fasta = _resolve_enterovirus_ev_db_fasta()
    if ev_typing_summary.is_file() and ev_db_fasta.is_file():
        for row in _load_rhinovirus_manifest(ev_typing_summary):
            big_group = str(row.get("big_group") or "").strip().upper()
            subtype = _normalize_enterovirus_type(str(row.get("subtype") or ""))
            status = str(row.get("status") or "").strip().lower()
            accession_full = str(row.get("query_id") or "").strip()
            accession = accession_full.split(".", 1)[0] if accession_full else ""
            if big_group not in {"A", "B", "C", "D"}:
                continue
            if not subtype or subtype.startswith("RV-"):
                continue
            if status not in {"typed", "review"}:
                continue
            if not accession:
                continue
            local_gff = (_resolve_enterovirus_db_dir() / "reference_genomes" / "gff3" / f"{accession}.gff3").resolve()
            local_vadr_gff = (_resolve_enterovirus_db_dir() / "reference_genomes" / "gff3_vadr" / f"{accession}.gff3").resolve()
            gff_path = ""
            if local_gff.is_file():
                gff_path = str(local_gff)
            elif local_vadr_gff.is_file():
                gff_path = str(local_vadr_gff)
            merged = dict(combined_by_accession.get(accession, {}))
            merged.update(
                {
                    "accession": accession,
                    "accession_full": accession_full,
                    "abbrev": subtype,
                    "virus_name": str(row.get("virus_name") or "").strip(),
                    "big_group": big_group,
                    "header": accession_full,
                    "fasta_path": str(ev_db_fasta),
                    "gff_path": gff_path,
                    "available": "Complete genome",
                    "source_record_id": accession_full,
                    "evidence": "ev_typing_out_summary",
                    "typing_status": status,
                }
            )
            combined_by_accession[accession] = merged

    return sorted(
        combined_by_accession.values(),
        key=lambda item: (
            str(item.get("big_group") or "").strip(),
            _normalize_enterovirus_type_from_row(item),
            str(item.get("accession") or item.get("accession_full") or "").strip(),
        ),
    )


def _get_fasta_record_by_id(fasta_path: Path, record_id: str = "") -> SeqRecord | None:
    fasta_key = str(fasta_path.resolve())
    wanted = str(record_id or "").strip()
    if wanted:
        cache = _FASTA_INDEX_CACHE.get(fasta_key)
        if cache is None:
            cache = SeqIO.index(fasta_key, "fasta")
            _FASTA_INDEX_CACHE[fasta_key] = cache
        candidates = [wanted, wanted.split()[0], wanted.split(".", 1)[0]]
        for candidate in candidates:
            if candidate and candidate in cache:
                return cache[candidate]
        for key in cache.keys():
            text = str(key or "").strip()
            if text == wanted or text.split(".", 1)[0] == wanted.split(".", 1)[0]:
                return cache[key]
        return None
    return next(SeqIO.parse(str(fasta_path), "fasta"), None)


def _normalize_enterovirus_type(value: str) -> str:
    text = str(value or "").strip().upper()
    if not text:
        return ""
    direct = re.search(r"\b(?:EV|ENTEROVIRUS)[-_ ]?([ABCD])[-_ ]?(\d+)\b", text)
    if direct:
        return f"EV-{direct.group(1)}{direct.group(2)}"
    cva = re.search(r"\b(?:CVA|COXSACKIEVIRUS A)[-_ ]?(\d+)\b", text)
    if cva:
        return f"CVA{cva.group(1)}"
    cvb = re.search(r"\b(?:CVB|COXSACKIEVIRUS B)[-_ ]?(\d+)\b", text)
    if cvb:
        return f"CVB{cvb.group(1)}"
    echo = re.search(r"\b(?:ECHO(?:VIRUS)?|E)[-_ ]?(\d+)\b", text)
    if echo:
        return f"E{echo.group(1)}"
    polio = re.search(r"\b(?:PV|POLIOVIRUS)[-_ ]?([123])\b", text)
    if polio:
        return f"PV{polio.group(1)}"
    compact = re.sub(r"[\s_]+", "", text)
    return compact


def _normalize_enterovirus_type_from_row(row: dict[str, object]) -> str:
    for key in ["abbrev", "virus_name", "accession", "header"]:
        normalized = _normalize_enterovirus_type(str(row.get(key) or ""))
        if normalized:
            return normalized
    return ""


def _extract_enterovirus_big_group(value: str) -> str:
    text = str(value or "").strip().upper()
    if not text:
        return ""
    if text in {"A", "B", "C", "D"}:
        return text
    direct = re.search(r"\bEV[-_ ]?([ABCD])\b", text)
    if direct:
        return direct.group(1)
    if "ENTEROVIRUS A" in text:
        return "A"
    if "ENTEROVIRUS B" in text:
        return "B"
    if "ENTEROVIRUS C" in text:
        return "C"
    if "ENTEROVIRUS D" in text:
        return "D"
    normalized_type = _normalize_enterovirus_type(text)
    for row in _load_enterovirus_abcd_manifest():
        if _normalize_enterovirus_type_from_row(row) == normalized_type:
            return str(row.get("big_group") or "").strip().upper()
    return ""


def _enterovirus_species_label(group: str) -> str:
    normalized = str(group or "").strip().upper()
    if normalized in {"A", "B", "C", "D"}:
        return f"Human enterovirus {normalized}"
    return "Human enterovirus"


def _build_enterovirus_subject_meta() -> dict[str, dict[str, str]]:
    subject_meta: dict[str, dict[str, str]] = {}
    for row in _load_enterovirus_abcd_manifest():
        accession = str(row.get("accession") or row.get("accession_full") or "").strip()
        accession_root = accession.split(".", 1)[0]
        record_id = accession_root
        normalized_type = _normalize_enterovirus_type_from_row(row)
        meta = {
            "big_group": str(row.get("big_group") or "").strip().upper(),
            "normalized_type": normalized_type,
            "abbrev": str(row.get("abbrev") or "").strip(),
            "accession": accession,
            "accession_root": accession_root,
            "gff_path": str(row.get("gff_path") or "").strip(),
            "fasta_path": str(row.get("fasta_path") or "").strip(),
        }
        for key in {record_id, accession, accession_root, normalized_type}:
            if key:
                subject_meta[key] = meta
    return subject_meta


def _score_enterovirus_type_rows(rows: list[dict[str, object]]) -> dict[str, str]:
    by_type: dict[str, dict[str, object]] = {}
    for row in rows:
        meta = row.get("_meta") if isinstance(row.get("_meta"), dict) else {}
        normalized_type = _normalize_enterovirus_type(str(meta.get("normalized_type") or ""))
        if not normalized_type:
            continue
        score = (
            float(row.get("coverage") or 0.0),
            float(row.get("mean_depth") or 0.0),
            float(row.get("covered_bases") or 0.0),
            float(row.get("num_reads") or 0.0),
        )
        bucket = by_type.setdefault(
            normalized_type,
            {
                "normalized_type": normalized_type,
                "big_group": str(meta.get("big_group") or ""),
                "subject": str(row.get("reference_name") or ""),
                "coverage": score[0],
                "mean_depth": score[1],
                "covered_bases": score[2],
                "num_reads": score[3],
                "reference_count": set(),
            },
        )
        if score > (
            float(bucket["coverage"]),
            float(bucket["mean_depth"]),
            float(bucket["covered_bases"]),
            float(bucket["num_reads"]),
        ):
            bucket["subject"] = str(row.get("reference_name") or "")
            bucket["coverage"] = score[0]
            bucket["mean_depth"] = score[1]
            bucket["covered_bases"] = score[2]
            bucket["num_reads"] = score[3]
        bucket["reference_count"].add(str(meta.get("accession_root") or str(row.get("reference_name") or "").split("_", 1)[0]))
    if not by_type:
        return {"type": "", "big_group": "", "subject": "", "coverage": "", "mean_depth": "", "covered_bases": "", "num_reads": "", "reference_count": "0"}
    ranked = sorted(
        by_type.values(),
        key=lambda item: (
            float(item["coverage"]),
            float(item["mean_depth"]),
            float(item["covered_bases"]),
            float(item["num_reads"]),
            len(item["reference_count"]),
        ),
        reverse=True,
    )
    best = ranked[0]
    return {
        "type": str(best["normalized_type"]),
        "big_group": str(best["big_group"]),
        "subject": str(best["subject"]),
        "coverage": f"{float(best['coverage']):.2f}",
        "mean_depth": f"{float(best['mean_depth']):.2f}",
        "covered_bases": f"{float(best['covered_bases']):.0f}",
        "num_reads": f"{float(best['num_reads']):.0f}",
        "reference_count": str(len(best["reference_count"])),
    }


def _run_enterovirus_vp1_read_typing(
    db_fasta: Path,
    out_dir: Path,
    single_fastq: str = "",
    fq1: str = "",
    fq2: str = "",
    long_type: str = "",
    threads: int = 4,
    logf=None,
) -> dict[str, str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    bam_path = out_dir / "vp1.screening.bam"
    coverage_path = out_dir / "vp1.coverage.tsv"
    if not bam_path.is_file():
        if fq1:
            map_cmd = " ".join(
                [
                    "minimap2",
                    "-ax",
                    "sr",
                    shlex.quote(str(db_fasta)),
                    shlex.quote(str(fq1)),
                    *([shlex.quote(str(fq2))] if fq2 else []),
                    "|",
                    "samtools",
                    "sort",
                    "-o",
                    shlex.quote(str(bam_path)),
                ]
            )
        else:
            map_cmd = " ".join(
                [
                    "minimap2",
                    "-ax",
                    _choose_minimap2_preset(long_type),
                    shlex.quote(str(db_fasta)),
                    shlex.quote(str(single_fastq)),
                    "-t",
                    str(max(1, int(threads or 1))),
                    "|",
                    "samtools",
                    "sort",
                    "-o",
                    shlex.quote(str(bam_path)),
                ]
            )
        run_command(map_cmd, logf=logf)
    if not bam_path.is_file():
        return {"type": "", "big_group": "", "subject": "", "identity": "", "coverage": "", "mean_depth": "", "covered_bases": "", "num_reads": "", "reference_count": "0", "method": "read_coverage"}
    run_command(f"samtools index -@ {max(1, int(threads or 1))} {shlex.quote(str(bam_path))}", logf=logf)
    run_command(f"samtools coverage {shlex.quote(str(bam_path))} > {shlex.quote(str(coverage_path))}", logf=logf)
    subject_meta = _build_enterovirus_subject_meta()
    rows = _parse_samtools_coverage_rows(coverage_path)
    enriched_rows = []
    for row in rows:
        subject = str(row.get("reference_name") or "").strip()
        meta = subject_meta.get(subject) or subject_meta.get(subject.split("_", 1)[0]) or {}
        enriched = dict(row)
        enriched["_meta"] = meta
        enriched_rows.append(enriched)
    best = _score_enterovirus_type_rows(enriched_rows)
    return {
        "type": best["type"],
        "big_group": best["big_group"],
        "subject": best["subject"],
        "identity": "",
        "coverage": best["coverage"],
        "mean_depth": best["mean_depth"],
        "covered_bases": best["covered_bases"],
        "num_reads": best["num_reads"],
        "reference_count": best["reference_count"],
        "method": "read_coverage",
    }


def _run_enterovirus_vp1_blast_typing(query_fasta: Path, db_fasta: Path, out_dir: Path, logf=None) -> dict[str, str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    out_tsv = out_dir / "vp1.blastn.tsv"
    blastn_bin = _resolve_hadv_blastn_bin()
    cmd = " ".join(
        [
            shlex.quote(blastn_bin),
            "-task blastn",
            "-query",
            shlex.quote(str(query_fasta)),
            "-subject",
            shlex.quote(str(db_fasta)),
            "-outfmt",
            shlex.quote("6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"),
            "-evalue 1e-20",
            "-max_target_seqs 100",
            "-dust no",
            "-num_threads 4",
            "-out",
            shlex.quote(str(out_tsv)),
        ]
    )
    run_command(cmd, logf=logf)
    if not out_tsv.is_file() or out_tsv.stat().st_size == 0:
        return {"type": "", "big_group": "", "subject": "", "identity": "", "coverage": "", "mean_depth": "", "covered_bases": "", "num_reads": "", "reference_count": "0", "method": "blastn"}
    subject_meta = _build_enterovirus_subject_meta()
    ref_lengths: dict[str, int] = {}
    with db_fasta.open("r", encoding="utf-8", errors="ignore") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            ref_lengths[str(record.id)] = len(str(record.seq))
    by_type: dict[str, dict[str, object]] = {}
    with out_tsv.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            _qseqid, sseqid, pident, length, _mismatch, _gapopen, _qstart, _qend, _sstart, _send, _evalue, bitscore = parts[:12]
            meta = subject_meta.get(sseqid) or subject_meta.get(sseqid.split("_", 1)[0]) or {}
            normalized_type = _normalize_enterovirus_type(str(meta.get("normalized_type") or ""))
            big_group = str(meta.get("big_group") or "").strip().upper()
            if not normalized_type:
                continue
            ref_len = ref_lengths.get(sseqid, 0)
            try:
                qcov_ref = min(100.0, (float(length) / ref_len) * 100) if ref_len else 0.0
                score = (float(bitscore), float(pident), qcov_ref, float(length))
            except ValueError:
                continue
            bucket = by_type.setdefault(
                normalized_type,
                {
                    "type": normalized_type,
                    "big_group": big_group,
                    "subject": sseqid,
                    "identity": str(pident),
                    "coverage": f"{qcov_ref:.2f}",
                    "covered_bases": str(length),
                    "reference_count": set(),
                    "_score": score,
                },
            )
            if score > bucket["_score"]:
                bucket["subject"] = sseqid
                bucket["identity"] = str(pident)
                bucket["coverage"] = f"{qcov_ref:.2f}"
                bucket["covered_bases"] = str(length)
                bucket["_score"] = score
                bucket["big_group"] = big_group
            bucket["reference_count"].add(str(meta.get("accession_root") or sseqid.split("_", 1)[0]))
    if not by_type:
        return {"type": "", "big_group": "", "subject": "", "identity": "", "coverage": "", "mean_depth": "", "covered_bases": "", "num_reads": "", "reference_count": "0", "method": "blastn"}
    ranked = sorted(
        by_type.values(),
        key=lambda item: (
            item["_score"][0],
            item["_score"][1],
            item["_score"][2],
            len(item["reference_count"]),
        ),
        reverse=True,
    )
    best = ranked[0]
    return {
        "type": str(best["type"]),
        "big_group": str(best["big_group"]),
        "subject": str(best["subject"]),
        "identity": str(best["identity"]),
        "coverage": str(best["coverage"]),
        "mean_depth": "",
        "covered_bases": str(best["covered_bases"]),
        "num_reads": "",
        "reference_count": str(len(best["reference_count"])),
        "method": "blastn",
    }


def _run_enterovirus_vp1_typing(
    pre: str,
    query_fasta: Path | None = None,
    single_fastq: str = "",
    fq1: str = "",
    fq2: str = "",
    long_type: str = "",
    threads: int = 4,
    output_dir: Path | None = None,
    logf=None,
) -> dict[str, object]:
    typing_dir = (output_dir or (Path(f"{pre}_enterovirus_reference_selection") / "typing")).resolve()
    typing_dir.mkdir(parents=True, exist_ok=True)
    db_fasta = _resolve_enterovirus_reference_vp1_fasta()
    if not db_fasta.is_file():
        return {"vp1_type": "", "big_group": "", "summary_path": "", "hit": {"method": "missing_db"}}
    use_blast = query_fasta is not None and query_fasta.is_file() and query_fasta.stat().st_size > 0
    hit = (
        _run_enterovirus_vp1_blast_typing(query_fasta, db_fasta, typing_dir, logf=logf)
        if use_blast
        else _run_enterovirus_vp1_read_typing(
            db_fasta,
            typing_dir,
            single_fastq=single_fastq,
            fq1=fq1,
            fq2=fq2,
            long_type=long_type,
            threads=threads,
            logf=logf,
        )
    )
    hit_status = "ready"
    hit_note = ""
    if _is_low_support_gene_typing_hit(hit):
        hit_status = "low_support"
        hit_note = _low_support_gene_typing_note("VP1", hit)
        hit = dict(hit)
        hit["type"] = ""
        hit["big_group"] = ""
    summary_path = typing_dir / "vp1_typing.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "gene", "matched_type", "big_group", "subject", "identity", "coverage", "mean_depth", "covered_bases", "num_reads", "reference_count", "method", "status", "note"])
        writer.writerow([
            pre,
            "vp1",
            hit.get("type", ""),
            hit.get("big_group", ""),
            hit.get("subject", ""),
            hit.get("identity", ""),
            hit.get("coverage", ""),
            hit.get("mean_depth", ""),
            hit.get("covered_bases", ""),
            hit.get("num_reads", ""),
            hit.get("reference_count", "0"),
            hit.get("method", ""),
            hit_status,
            hit_note,
        ])
    return {
        "vp1_type": str(hit.get("type") or ""),
        "big_group": str(hit.get("big_group") or ""),
        "summary_path": str(summary_path.resolve()),
        "hit": hit,
        "status": hit_status,
        "note": hit_note,
    }


def _resolve_enterovirus_vadr_model_dir(project_root: Path, big_group: str) -> Path | None:
    normalized = str(big_group or "").strip().upper()
    if normalized not in {"A", "B", "C", "D"}:
        return None
    model_dir = (project_root / "soft" / "vadr-models-ev" / f"ev{normalized}").resolve()
    if model_dir.is_dir():
        return model_dir
    return None


def _run_vadr_enterovirus_annotation(pre: str, input_fasta: Path, output_root: Path, big_group: str, logf=None) -> dict[str, str]:
    project_root = _project_root()
    normalized = str(big_group or "").strip().upper()
    model_dir = _resolve_enterovirus_vadr_model_dir(project_root, normalized)
    if model_dir is None:
        return {"status": "missing", "note": f"未找到 Enterovirus {normalized or '?'} VADR 模型目录", "output_dir": "", "gff_path": ""}
    env = _build_vadr_env(project_root, model_dir)
    if env is None:
        return {"status": "missing", "note": "VADR 运行依赖不完整", "output_dir": "", "gff_path": ""}
    enterovirus_env = dict(env)
    alt_bio_easel_dir = project_root / "soft" / "Bio-Easel-ncov"
    alt_bio_lib = alt_bio_easel_dir / "blib" / "lib"
    alt_bio_arch = alt_bio_easel_dir / "blib" / "arch"
    if alt_bio_lib.is_dir() and alt_bio_arch.is_dir():
        enterovirus_env["VADRBIOEASELDIR"] = str(alt_bio_easel_dir)
        enterovirus_env["PERL5LIB"] = os.pathsep.join(
            [
                str(project_root / "soft" / "vadr"),
                str(project_root / "soft" / "sequip"),
                str(alt_bio_lib),
                str(alt_bio_arch),
            ]
        )
    ncov_bin_dir = conda_env_path("ncov", "bin")
    enterovirus_env["PATH"] = os.pathsep.join(
        [
            ncov_bin_dir,
            "/usr/bin",
            "/bin",
            str(enterovirus_env.get("PATH") or ""),
        ]
    )
    perl_candidates = [f"{ncov_bin_dir}/perl", "/usr/bin/perl", _resolve_vadr_perl_bin(), "perl"]
    perl_bin = _resolve_working_perl_with_module(enterovirus_env, "Bio::Easel::MSA", perl_candidates)
    vadr_root = output_root / "vadr"
    vadr_root.mkdir(parents=True, exist_ok=True)
    output_dir = vadr_root / f"{pre}_vadr"
    prefix = f"{output_dir.name}.vadr"
    pass_tbl = output_dir / f"{prefix}.pass.tbl"
    fail_tbl = output_dir / f"{prefix}.fail.tbl"
    gff_path = vadr_root / f"{pre}.vadr.gff3"
    if output_dir.is_dir() and (pass_tbl.is_file() or fail_tbl.is_file()) and gff_path.is_file() and gff_path.stat().st_size > 0:
        return {"status": "ready", "note": "已检测到现有肠道病毒 VADR 注释结果", "output_dir": str(output_dir.resolve()), "gff_path": str(gff_path.resolve())}
    vadr_script = project_root / "soft" / "vadr" / "v-annotate.pl"
    cmd = " ".join(
        [
            shlex.quote(str(perl_bin)),
            shlex.quote(str(vadr_script)),
            "-f",
            "-r",
            "--ignore_exc",
            "--mkey",
            f"ev{normalized}",
            "--mdir",
            shlex.quote(str(model_dir)),
            shlex.quote(str(input_fasta)),
            shlex.quote(str(output_dir)),
        ]
    )
    try:
        run_command(cmd, logf=logf, env=enterovirus_env)
    except Exception as exc:
        return {"status": "failed", "note": f"肠道病毒 VADR 注释失败: {exc}", "output_dir": str(output_dir.resolve()) if output_dir.exists() else "", "gff_path": ""}
    annotate_tbl2gff = project_root / "soft" / "vadr" / "miniscripts" / "annotate-tbl2gff.pl"
    source_tbl = pass_tbl if _tbl_has_feature_rows(pass_tbl) else fail_tbl
    if _tbl_has_feature_rows(source_tbl):
        try:
            run_command(
                f"{shlex.quote(str(perl_bin))} {shlex.quote(str(annotate_tbl2gff))} {shlex.quote(str(source_tbl))} > {shlex.quote(str(gff_path))}",
                logf=logf,
                env=enterovirus_env,
            )
        except Exception:
            pass
    return {
        "status": "ready" if gff_path.is_file() and gff_path.stat().st_size > 0 else "failed",
        "note": "肠道病毒 VADR 注释完成" if gff_path.is_file() and gff_path.stat().st_size > 0 else "肠道病毒 VADR 未生成可用 GFF",
        "output_dir": str(output_dir.resolve()),
        "gff_path": str(gff_path.resolve()) if gff_path.is_file() and gff_path.stat().st_size > 0 else "",
    }


def prepare_enterovirus_sample_annotation(pre: str, sample_fasta: Path, output_root: Path, big_group: str, logf=None) -> Path | None:
    result = _run_vadr_enterovirus_annotation(pre, sample_fasta, output_root, big_group, logf=logf)
    gff_path = str(result.get("gff_path") or "").strip()
    if result.get("status") == "ready" and gff_path:
        gff_file = Path(gff_path)
        if gff_file.is_file() and gff_file.stat().st_size > 0:
            return gff_file
    return None


def prepare_enterovirus_reference_annotation(pre: str, reference_fasta: Path, output_root: Path, big_group: str, logf=None) -> Path | None:
    result = _run_vadr_enterovirus_annotation(pre, reference_fasta, output_root, big_group, logf=logf)
    gff_path = str(result.get("gff_path") or "").strip()
    if result.get("status") == "ready" and gff_path:
        gff_file = Path(gff_path)
        if gff_file.is_file() and gff_file.stat().st_size > 0:
            return gff_file
    return None


def _tbl_has_feature_rows(tbl_path: Path) -> bool:
    if not tbl_path.is_file() or tbl_path.stat().st_size == 0:
        return False
    try:
        with tbl_path.open("r", encoding="utf-8", errors="ignore") as handle:
            for raw_line in handle:
                if raw_line.startswith(">Feature "):
                    return True
    except OSError:
        return False
    return False


def _find_enterovirus_vp1_feature(gff_path: Path) -> dict[str, object] | None:
    return _find_rhinovirus_vp1_feature(gff_path)


def _extract_enterovirus_vp1_region(sample_fasta: Path, gff_path: Path, sample_label: str) -> tuple[SeqRecord | None, dict[str, object]]:
    return _extract_rhinovirus_vp1_region(sample_fasta, gff_path, sample_label)


def _build_enterovirus_tree_member_label(subtype: str, member_label: str, suffix: str = "") -> str:
    subtype_text = str(subtype or "").strip()
    member_text = str(member_label or "").strip()
    parts = [item for item in [subtype_text, member_text, str(suffix or "").strip()] if item]
    return _sanitize_tree_label("_".join(parts))


def _resolve_enterovirus_vp1_fallback_gff(vp1_type: str, big_group: str) -> Path | None:
    normalized_type = _normalize_enterovirus_type(vp1_type)
    normalized_group = str(big_group or "").strip().upper()
    group_match: Path | None = None
    for row in _load_enterovirus_abcd_manifest():
        row_group = str(row.get("big_group") or "").strip().upper()
        if normalized_group and row_group != normalized_group:
            continue
        gff_path = Path(str(row.get("gff_path") or "").strip())
        if not gff_path.is_file() or gff_path.stat().st_size == 0:
            continue
        row_type = _normalize_enterovirus_type_from_row(row)
        if normalized_type and row_type == normalized_type:
            return gff_path
        if group_match is None:
            group_match = gff_path
    return group_match


def _filter_enterovirus_vp1_refs_by_group(big_group: str) -> tuple[list[dict[str, str]], list[SeqRecord]]:
    target_group = str(big_group or "").strip().upper()
    selected_records: list[SeqRecord] = []
    ordered_rows: list[dict[str, str]] = []
    for row in _load_enterovirus_abcd_manifest():
        if str(row.get("big_group") or "").strip().upper() != target_group:
            continue
        vp1_fasta_path = Path(str(row.get("vp1_fasta_path") or "").strip())
        if not vp1_fasta_path.is_file():
            continue
        record = next(SeqIO.parse(str(vp1_fasta_path), "fasta"), None)
        if record is None:
            continue
        subtype = _normalize_enterovirus_type_from_row(row)
        accession = str(row.get("accession") or row.get("accession_full") or "").split(".", 1)[0]
        tree_label = _build_enterovirus_tree_member_label(subtype, accession or record.id)
        ordered_rows.append(
            {
                "gene": "VP1",
                "subtype": subtype,
                "accession": accession,
                "backup_rank": "1",
                "tree_label": tree_label,
            }
        )
        selected_records.append(
            SeqRecord(
                record.seq,
                id=tree_label,
                name=tree_label,
                description=f"{subtype} {accession}".strip(),
            )
        )
    return ordered_rows, selected_records


def build_enterovirus_vp1_phylogeny_assets(pre: str, sample_fasta: Path, gff_path: Path, big_group: str = "", logf=None) -> dict[str, object]:
    output_root = Path(f"{pre}_enterovirus_reference_selection") / "phylogeny"
    output_root.mkdir(parents=True, exist_ok=True)
    selection_path = Path(f"{pre}_enterovirus_reference_selection") / "selection.tsv"
    vp1_type = ""
    selected_group = str(big_group or "").strip().upper()
    if selection_path.is_file():
        with selection_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
            row = next(csv.DictReader(handle, delimiter="\t"), None) or {}
        vp1_type = str(row.get("vp1_type") or "").strip()
        selected_group = selected_group or str(row.get("big_group") or "").strip().upper()
    result_row = {
        "gene": "VP1",
        "subtype": vp1_type,
        "genogroup": selected_group,
        "status": "missing",
        "note": "",
        "sample_region_path": "",
        "backup_fasta": "",
        "tree_path": "",
        "member_count": "0",
    }
    sample_record, sample_meta = _extract_enterovirus_vp1_region(sample_fasta, gff_path, _sanitize_tree_label(pre))
    if sample_record is None:
        fallback_gff = _resolve_enterovirus_vp1_fallback_gff(vp1_type, selected_group)
        if fallback_gff is not None and fallback_gff != gff_path:
            sample_record, sample_meta = _extract_enterovirus_vp1_region(sample_fasta, fallback_gff, _sanitize_tree_label(pre))
            if sample_record is not None:
                sample_meta = dict(sample_meta)
                sample_meta["note"] = f"样本 VADR 注释未定位 VP1，已回退到 {fallback_gff.name} 坐标"
    if sample_record is None:
        result_row["note"] = str(sample_meta.get("note") or "未提取到 VP1 区域")
    else:
        sample_label = _build_enterovirus_tree_member_label(vp1_type or selected_group or "Enterovirus", pre, "sample")
        sample_record.id = sample_label
        sample_record.name = sample_label
        sample_record.description = f"{vp1_type or selected_group or 'Enterovirus'} {pre} sample".strip()
        backup_rows, backup_records = _filter_enterovirus_vp1_refs_by_group(selected_group)
        if not backup_records:
            result_row["note"] = f"VP1 参考库中未找到 Enterovirus {selected_group or '?'} 记录"
        else:
            tree_result = _build_norovirus_gene_tree("vp1", sample_record, backup_rows, backup_records, output_root / "vp1", logf=logf)
            result_row["status"] = str(tree_result.get("status") or "missing")
            fallback_note = str(sample_meta.get("note") or "").strip()
            tree_note = str(tree_result.get("note") or "").strip()
            result_row["note"] = "；".join([item for item in [fallback_note, tree_note] if item])
            result_row["sample_region_path"] = str((output_root / "vp1" / "vp1.sample.fasta").resolve())
            result_row["backup_fasta"] = str(tree_result.get("backup_fasta") or "")
            result_row["tree_path"] = str(tree_result.get("tree_path") or "")
            result_row["member_count"] = str(tree_result.get("member_count") or "0")
    summary_path = output_root / "summary.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["gene", "subtype", "genogroup", "status", "note", "sample_region_path", "backup_fasta", "tree_path", "member_count"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(result_row)
    return {
        "status": "ready" if result_row.get("status") == "ready" else "missing",
        "summary_path": str(summary_path.resolve()),
        "rows": [result_row],
    }


def resolve_enterovirus_reference(
    pre: str,
    species: str = "",
    requested_ref: str = "",
    single_fastq: str = "",
    fq1: str = "",
    fq2: str = "",
    long_type: str = "",
    threads: int = 4,
    query_fasta: str = "",
    logf=None,
) -> dict[str, str]:
    _clear_sample_skip_flag(pre)
    requested = str(requested_ref or "").strip()
    if requested and Path(requested).is_file():
        big_group = _extract_enterovirus_big_group(species)
        return {
            "status": "ready",
            "vp1_type": _normalize_enterovirus_type(species),
            "big_group": big_group,
            "species_label": _enterovirus_species_label(big_group),
            "reference_path": str(Path(requested).resolve()),
            "gff_path": "nogtf",
            "summary_path": "",
            "typing_summary_path": "",
        }
    screening_dir = Path(f"{pre}_enterovirus_reference_selection")
    screening_dir.mkdir(exist_ok=True)
    query_path = Path(query_fasta) if str(query_fasta or "").strip() else None
    typing_result = _run_enterovirus_vp1_typing(
        pre,
        query_fasta=query_path,
        single_fastq=single_fastq,
        fq1=fq1,
        fq2=fq2,
        long_type=long_type,
        threads=threads,
        output_dir=screening_dir / "typing",
        logf=logf,
    )
    vp1_type = str(typing_result.get("vp1_type") or "").strip()
    big_group = str(typing_result.get("big_group") or "").strip().upper()
    summary_path = screening_dir / "selection.tsv"
    if not vp1_type or big_group not in {"A", "B", "C", "D"}:
        status_text = "skipped" if str(typing_result.get("status") or "").strip() == "low_support" else "missing"
        note_text = str(typing_result.get("note") or "").strip() or f"未能根据 VP1 或 species={species or '-'} 确定肠道病毒分型"
        with summary_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sample", "vp1_type", "big_group", "status", "note", "typing_summary_path"])
            writer.writerow([pre, vp1_type, big_group, status_text, note_text, str(typing_result.get("summary_path") or "")])
        if status_text == "skipped":
            _write_sample_skip_flag(pre, note_text)
        return {
            "status": status_text,
            "vp1_type": vp1_type,
            "big_group": big_group,
            "species_label": species or "Human enterovirus",
            "reference_path": "",
            "gff_path": "nogtf",
            "summary_path": str(summary_path.resolve()),
            "typing_summary_path": str(typing_result.get("summary_path") or ""),
        }
    manifest_rows = _load_enterovirus_candidate_reference_manifest()
    exact_candidates = [
        row for row in manifest_rows
        if str(row.get("big_group") or "").strip().upper() == big_group
        and _normalize_enterovirus_type_from_row(row) == vp1_type
    ]
    candidates = exact_candidates or [row for row in manifest_rows if str(row.get("big_group") or "").strip().upper() == big_group]
    candidate_records = _build_enterovirus_candidate_records(candidates)
    dedup_candidate_records = _deduplicate_enterovirus_candidate_records(
        candidate_records,
        screening_dir=screening_dir,
        threads=threads,
        logf=logf,
    )
    if not dedup_candidate_records:
        with summary_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sample", "vp1_type", "big_group", "status", "note", "typing_summary_path"])
            writer.writerow([pre, vp1_type, big_group, "missing", "未找到对应亚型的全基因组候选参考", str(typing_result.get("summary_path") or "")])
        return {
            "status": "missing",
            "vp1_type": vp1_type,
            "big_group": big_group,
            "species_label": _enterovirus_species_label(big_group),
            "reference_path": "",
            "gff_path": "nogtf",
            "summary_path": str(summary_path.resolve()),
            "typing_summary_path": str(typing_result.get("summary_path") or ""),
        }
    final_rows, audit_rows, anchor_record = _select_enterovirus_candidates_in_batches(
        dedup_candidate_records,
        screening_dir=screening_dir,
        single_fastq=single_fastq,
        fq1=fq1,
        fq2=fq2,
        long_type=long_type,
        threads=threads,
        logf=logf,
    )
    candidate_fasta = screening_dir / "candidate_references.fasta"
    candidate_meta = _write_reference_records(dedup_candidate_records, candidate_fasta)
    coverage_path = screening_dir / "candidate_references.coverage.tsv"
    coverage_rows = final_rows or audit_rows
    _write_aggregated_coverage_rows(coverage_path, coverage_rows)
    dedup_manifest_path = screening_dir / "candidate_references.dedup.tsv"
    with dedup_manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["representative_id", "accession", "member_count", "member_accessions"])
        for record in dedup_candidate_records:
            writer.writerow([record.get("fasta_id", ""), record.get("accession", ""), record.get("member_count", 1), record.get("member_accessions", "")])
    best = None
    for row in coverage_rows:
        ref_name = str(row.get("reference_name") or "").strip()
        meta = candidate_meta.get(ref_name)
        if meta is None:
            continue
        score = _coverage_row_score(row)
        if best is None or score > best["score"]:
            best = {"score": score, "meta": meta, "coverage_row": row}
    if best is None:
        return {
            "status": "missing",
            "vp1_type": vp1_type,
            "big_group": big_group,
            "species_label": _enterovirus_species_label(big_group),
            "reference_path": "",
            "gff_path": "nogtf",
            "summary_path": "",
            "typing_summary_path": str(typing_result.get("summary_path") or ""),
        }
    best_meta = best["meta"]
    best_row = best["coverage_row"]
    best_accession = str(best_meta.get("accession") or "").strip()
    anchor_accession = str(anchor_record.get("accession") or "").strip() if anchor_record else ""
    best_reference_fasta = screening_dir / f"{_sanitize_tree_label(vp1_type or best_accession)}.reference.fasta"
    with best_reference_fasta.open("w", encoding="utf-8") as handle:
        handle.write(f">{best_accession or 'reference'}\n")
        sequence = str(best_meta.get("sequence") or "")
        for index in range(0, len(sequence), 80):
            handle.write(sequence[index:index + 80] + "\n")
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "vp1_type", "big_group", "coverage", "mean_depth", "covered_bases", "num_reads", "reference_name", "accession", "reference_path", "gff_path", "typing_summary_path", "anchor_accession", "candidate_count", "dedup_candidate_count"])
        writer.writerow([
            pre,
            vp1_type,
            big_group,
            f"{float(best_row.get('coverage') or 0.0):.6f}",
            f"{float(best_row.get('mean_depth') or 0.0):.6f}",
            f"{float(best_row.get('covered_bases') or 0.0):.0f}",
            f"{float(best_row.get('num_reads') or 0.0):.0f}",
            str(best_meta.get("meta", {}).get("header") or best_accession),
            best_accession,
            str(best_reference_fasta.resolve()),
            str(best_meta.get("meta", {}).get("gff_path") or "nogtf"),
            str(typing_result.get("summary_path") or ""),
            anchor_accession,
            len(candidate_records),
            len(dedup_candidate_records),
        ])
    return {
        "status": "ready",
        "vp1_type": vp1_type,
        "big_group": big_group,
        "species_label": _enterovirus_species_label(big_group),
        "reference_path": str(best_reference_fasta.resolve()),
        "gff_path": str(best_meta.get("meta", {}).get("gff_path") or "nogtf"),
        "summary_path": str(summary_path.resolve()),
        "typing_summary_path": str(typing_result.get("summary_path") or ""),
    }


def run_enterovirus_consensus_typing(pre: str, consensus_fasta: Path, logf=None) -> dict[str, str]:
    if not consensus_fasta.is_file() or consensus_fasta.stat().st_size == 0:
        return {"status": "missing", "summary_path": "", "vp1_type": "", "big_group": ""}
    output_dir = Path(f"{pre}_enterovirus_reference_selection") / "consensus_typing"
    result = _run_enterovirus_vp1_typing(pre, query_fasta=consensus_fasta, output_dir=output_dir, logf=logf)
    summary_path = output_dir / "consensus_typing.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "vp1_type", "big_group", "typing_summary_path"])
        writer.writerow([pre, result.get("vp1_type", ""), result.get("big_group", ""), result.get("summary_path", "")])
    return {
        "status": "ready",
        "summary_path": str(summary_path.resolve()),
        "vp1_type": str(result.get("vp1_type") or ""),
        "big_group": str(result.get("big_group") or ""),
    }


def _resolve_bandavirus_db_dir() -> Path:
    env_root = str(os.environ.get("META_BANDAVIRUS_DB_DIR") or "").strip()
    if env_root:
        return Path(env_root).expanduser().resolve()
    database_root = str(os.environ.get("META_DATABASE_ROOT") or "").strip()
    project_root = _project_root()
    candidates: list[Path] = []
    if database_root:
        candidates.append((Path(database_root) / "virus" / "bandavirus").expanduser())
    candidates.extend(
        [
            project_root / "database" / "virus" / "bandavirus",
            Path("/data/deploy/meta_genome/database/virus/bandavirus"),
        ]
    )
    for candidate in candidates:
        if candidate.is_dir():
            return candidate.resolve()
    return candidates[0].resolve()


def _resolve_orthohantavirus_db_dir() -> Path:
    env_root = str(os.environ.get("META_ORTHOHANTAVIRUS_DB_DIR") or "").strip()
    if env_root:
        return Path(env_root).expanduser().resolve()
    database_root = str(os.environ.get("META_DATABASE_ROOT") or "").strip()
    project_root = _project_root()
    candidates: list[Path] = []
    if database_root:
        candidates.append((Path(database_root) / "virus" / "orthohantavirus").expanduser())
    candidates.extend(
        [
            project_root / "database" / "virus" / "orthohantavirus",
            Path("/data/deploy/meta_genome/database/virus/orthohantavirus"),
        ]
    )
    for candidate in candidates:
        if candidate.is_dir():
            return candidate.resolve()
    return candidates[0].resolve()


def _resolve_orthohantavirus_typing_manifest() -> Path:
    return (_resolve_orthohantavirus_db_dir() / "reference_genomes" / "orthohantavirus_typing_reference_segments_manifest.tsv").resolve()


def _resolve_orthoebolavirus_db_dir() -> Path:
    env_root = str(os.environ.get("META_ORTHOEBOLAVIRUS_DB_DIR") or "").strip()
    if env_root:
        return Path(env_root).expanduser().resolve()
    database_root = str(os.environ.get("META_DATABASE_ROOT") or "").strip()
    project_root = _project_root()
    candidates: list[Path] = []
    if database_root:
        candidates.append((Path(database_root) / "virus" / "Orthoebolavirus").expanduser())
    candidates.extend(
        [
            project_root / "database" / "virus" / "Orthoebolavirus",
            Path("/data/deploy/meta_genome/database/virus/Orthoebolavirus"),
        ]
    )
    for candidate in candidates:
        if candidate.is_dir():
            return candidate.resolve()
    return candidates[0].resolve()


def _resolve_orthoebolavirus_typing_manifest() -> Path:
    return (_resolve_orthoebolavirus_db_dir() / "reference_genomes" / "orthoebolavirus_typing_reference_genomes_manifest.tsv").resolve()


def _resolve_bandavirus_typing_manifest() -> Path:
    return (_resolve_bandavirus_db_dir() / "reference_genomes" / "bandavirus_typing_reference_segments_manifest.tsv").resolve()


def _resolve_bandavirus_grouped_csv() -> Path:
    return (_resolve_bandavirus_db_dir() / "bandavirus_db_grouped.csv").resolve()


def _resolve_bandavirus_grouped_fasta() -> Path:
    return (_resolve_bandavirus_db_dir() / "bandavirus_db_grouped.fasta").resolve()


def _resolve_bandavirus_af_manifest() -> Path:
    return (_resolve_bandavirus_db_dir() / "grouped_typing_refs" / "A_Fgroup" / "SFTSV_tree_group_segments_manifest.tsv").resolve()


def _resolve_bandavirus_af_summary() -> Path:
    return (_resolve_bandavirus_db_dir() / "grouped_typing_refs" / "A_Fgroup" / "SFTSV_tree_group_summary.tsv").resolve()


def _resolve_bandavirus_cj_manifest() -> Path:
    return (_resolve_bandavirus_db_dir() / "grouped_typing_refs" / "CJ_group" / "CJ_group_segment_genomes_manifest.tsv").resolve()


def _load_bandavirus_grouped_rows() -> list[dict[str, str]]:
    path = _resolve_bandavirus_grouped_csv()
    if not path.is_file():
        return []
    with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle) if row]


def _read_bandavirus_grouped_fasta_map(path: Path | None = None) -> dict[str, dict[str, str]]:
    fasta_path = path or _resolve_bandavirus_grouped_fasta()
    fasta_map: dict[str, dict[str, str]] = {}
    if not fasta_path.is_file():
        return fasta_map
    header = ""
    seq_lines: list[str] = []
    current_accession = ""
    current_segment = ""
    current_sample = ""
    for raw_line in fasta_path.open("r", encoding="utf-8", errors="ignore"):
        line = raw_line.rstrip("\n")
        if line.startswith(">"):
            if header and current_accession and current_segment:
                fasta_map[f"{current_accession}|{current_segment}"] = {
                    "header": header,
                    "sequence": "".join(seq_lines).upper(),
                    "accession": current_accession,
                    "segment": current_segment,
                    "sample": current_sample,
                }
            header = line[1:].strip()
            meta: dict[str, str] = {}
            for part in [item.strip() for item in header.split("|")]:
                if "=" not in part:
                    continue
                key, value = part.split("=", 1)
                meta[key.strip()] = value.strip()
            current_accession = str(meta.get("accession") or "").split(".", 1)[0]
            current_segment = str(meta.get("segment") or "").strip().upper()
            current_sample = str(meta.get("sample") or "").strip()
            seq_lines = []
        else:
            seq_lines.append(line.strip())
    if header and current_accession and current_segment:
        fasta_map[f"{current_accession}|{current_segment}"] = {
            "header": header,
            "sequence": "".join(seq_lines).upper(),
            "accession": current_accession,
            "segment": current_segment,
            "sample": current_sample,
        }
    return fasta_map


def _bandavirus_species_label(abbrev: str, fallback: str = "") -> str:
    normalized = str(abbrev or "").strip().upper()
    if normalized == "SFTSV":
        return "Bandavirus dabieense"
    for row in _load_rhinovirus_manifest(_resolve_bandavirus_typing_manifest()):
        if str(row.get("abbrev") or "").strip().upper() == normalized:
            species_label = str(row.get("species") or "").strip()
            if species_label:
                return species_label
    fallback_text = str(fallback or "").strip()
    return fallback_text or "Bandavirus"


def _bandavirus_is_sftsv(abbrev: str) -> bool:
    normalized = str(abbrev or "").strip().upper()
    return normalized == "SFTSV" or _normalize_species_label(abbrev) in SFTSV_LABELS


def _build_bandavirus_typing_records() -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for row in _load_rhinovirus_manifest(_resolve_bandavirus_typing_manifest()):
        fasta_path = Path(str(row.get("fasta_path") or "").strip())
        abbrev = str(row.get("abbrev") or "").strip().upper()
        segment = str(row.get("segment") or "").strip().upper()
        accession = str(row.get("accession") or "").strip().split(".", 1)[0]
        if not abbrev or segment not in {"L", "M", "S"} or not accession or not fasta_path.is_file():
            continue
        record = _get_fasta_record_by_id(fasta_path)
        if record is None:
            continue
        sequence = str(record.seq).strip().upper()
        if not sequence:
            continue
        records.append(
            {
                "fasta_id": f"{abbrev}__{segment}__{accession}",
                "accession": accession,
                "sequence": sequence,
                "reference_length": len(sequence),
                "meta": dict(row),
                "source_record_id": str(record.id),
            }
        )
    return records


def _orthohantavirus_species_label(abbrev: str, fallback: str = "") -> str:
    normalized = str(abbrev or "").strip().upper()
    for row in _load_rhinovirus_manifest(_resolve_orthohantavirus_typing_manifest()):
        candidate = str(row.get("abbrev_primary") or row.get("abbrev") or "").strip().upper()
        if candidate == normalized:
            species_label = str(row.get("species") or "").strip()
            if species_label:
                return species_label
    fallback_text = str(fallback or "").strip()
    return fallback_text or "Orthohantavirus"


def _build_orthohantavirus_typing_records(segment: str = "") -> list[dict[str, object]]:
    requested_segment = str(segment or "").strip().upper()
    records: list[dict[str, object]] = []
    for row in _load_rhinovirus_manifest(_resolve_orthohantavirus_typing_manifest()):
        fasta_path = Path(str(row.get("fasta_path") or "").strip())
        abbrev = str(row.get("abbrev_primary") or row.get("abbrev") or "").strip().upper()
        segment_label = str(row.get("segment") or "").strip().upper()
        accession = str(row.get("accession_root") or row.get("accession") or "").strip().split(".", 1)[0]
        if requested_segment and segment_label != requested_segment:
            continue
        if not abbrev or segment_label not in {"L", "M", "S"} or not accession or not fasta_path.is_file():
            continue
        record = _get_fasta_record_by_id(fasta_path)
        if record is None:
            continue
        sequence = str(record.seq).strip().upper()
        if not sequence:
            continue
        records.append(
            {
                "fasta_id": f"{abbrev}__{segment_label}__{accession}",
                "accession": accession,
                "sequence": sequence,
                "reference_length": len(sequence),
                "meta": dict(row),
                "source_record_id": str(record.id),
            }
        )
    return records


def _build_bandavirus_af_records(group: str = "", segment: str = "") -> list[dict[str, object]]:
    grouped_rows = _load_bandavirus_grouped_rows()
    grouped_index = {
        f"{str(row.get('Accession') or '').strip().split('.', 1)[0]}|{str(row.get('Assigned_Segment') or '').strip().upper()}": row
        for row in grouped_rows
    }
    fasta_map = _read_bandavirus_grouped_fasta_map()
    records: list[dict[str, object]] = []
    for row in _load_rhinovirus_manifest(_resolve_bandavirus_af_manifest()):
        tree_group = str(row.get("Tree_Group") or "").strip()
        tree_segment = str(row.get("Tree_Segment") or "").strip().upper()
        accession = str(row.get("Tree_Resolved_Accession") or row.get("Accession") or "").strip().split(".", 1)[0]
        if group and tree_group != group:
            continue
        if segment and tree_segment != str(segment).strip().upper():
            continue
        key = f"{accession}|{tree_segment}"
        fasta_meta = fasta_map.get(key)
        grouped_row = grouped_index.get(key, {})
        if not fasta_meta:
            continue
        sequence = str(fasta_meta.get("sequence") or "").strip().upper()
        if not sequence:
            continue
        meta = dict(row)
        meta.update(grouped_row)
        records.append(
            {
                "fasta_id": f"{tree_group}__{tree_segment}__{accession}",
                "accession": accession,
                "sequence": sequence,
                "reference_length": len(sequence),
                "meta": meta,
                "source_record_id": accession,
            }
        )
    return records


def _build_bandavirus_cj_records(genotype: str = "", segment: str = "") -> list[dict[str, object]]:
    grouped_rows = _load_bandavirus_grouped_rows()
    grouped_index = {
        f"{str(row.get('Accession') or '').strip().split('.', 1)[0]}|{str(row.get('Assigned_Segment') or '').strip().upper()}": row
        for row in grouped_rows
    }
    fasta_map = _read_bandavirus_grouped_fasta_map()
    records: list[dict[str, object]] = []
    for row in _load_rhinovirus_manifest(_resolve_bandavirus_cj_manifest()):
        cj_genotype = str(row.get("CJ_Genotype") or "").strip()
        cj_segment = str(row.get("CJ_Segment") or "").strip().upper()
        accession = str(row.get("Accession") or "").strip().split(".", 1)[0]
        if genotype and cj_genotype != genotype:
            continue
        if segment and cj_segment != str(segment).strip().upper():
            continue
        key = f"{accession}|{cj_segment}"
        fasta_meta = fasta_map.get(key)
        grouped_row = grouped_index.get(key, {})
        if not fasta_meta:
            continue
        sequence = str(fasta_meta.get("sequence") or "").strip().upper()
        if not sequence:
            continue
        meta = dict(row)
        meta.update(grouped_row)
        records.append(
            {
                "fasta_id": f"{cj_genotype}__{cj_segment}__{accession}",
                "accession": accession,
                "sequence": sequence,
                "reference_length": len(sequence),
                "meta": meta,
                "source_record_id": accession,
            }
        )
    return records


def _aggregate_bandavirus_rows(
    rows: list[dict[str, object]],
    meta_by_id: dict[str, dict[str, object]],
    label_key: str,
    segment_key: str,
) -> dict[str, str]:
    by_label: dict[str, dict[str, object]] = {}
    for row in rows:
        ref_name = str(row.get("reference_name") or "").strip()
        meta_entry = meta_by_id.get(ref_name) or {}
        meta = meta_entry.get("meta") if isinstance(meta_entry.get("meta"), dict) else {}
        label = str(meta.get(label_key) or "").strip()
        segment = str(meta.get(segment_key) or meta.get("segment") or "").strip().upper()
        if not label:
            continue
        score = _coverage_row_score(row)
        bucket = by_label.setdefault(
            label,
            {
                "label": label,
                "coverage_sum": 0.0,
                "depth_sum": 0.0,
                "covered_bases_sum": 0.0,
                "num_reads_sum": 0.0,
                "segment_hits": set(),
                "segment_best": {},
            },
        )
        previous = bucket["segment_best"].get(segment)
        if previous is None or score > previous["score"]:
            if previous is not None:
                bucket["coverage_sum"] -= float(previous["row"].get("coverage") or 0.0)
                bucket["depth_sum"] -= float(previous["row"].get("mean_depth") or 0.0)
                bucket["covered_bases_sum"] -= float(previous["row"].get("covered_bases") or 0.0)
                bucket["num_reads_sum"] -= float(previous["row"].get("num_reads") or 0.0)
            bucket["segment_best"][segment] = {"score": score, "row": row, "meta": meta_entry}
            bucket["coverage_sum"] += float(row.get("coverage") or 0.0)
            bucket["depth_sum"] += float(row.get("mean_depth") or 0.0)
            bucket["covered_bases_sum"] += float(row.get("covered_bases") or 0.0)
            bucket["num_reads_sum"] += float(row.get("num_reads") or 0.0)
            if segment:
                bucket["segment_hits"].add(segment)
    if not by_label:
        return {
            "label": "",
            "coverage_sum": "0",
            "depth_sum": "0",
            "covered_bases_sum": "0",
            "num_reads_sum": "0",
            "segment_count": "0",
        }
    ranked = sorted(
        by_label.values(),
        key=lambda item: (
            len(item["segment_hits"]),
            float(item["coverage_sum"]),
            float(item["depth_sum"]),
            float(item["covered_bases_sum"]),
            float(item["num_reads_sum"]),
        ),
        reverse=True,
    )
    best = ranked[0]
    return {
        "label": str(best["label"]),
        "coverage_sum": f"{float(best['coverage_sum']):.2f}",
        "depth_sum": f"{float(best['depth_sum']):.2f}",
        "covered_bases_sum": f"{float(best['covered_bases_sum']):.0f}",
        "num_reads_sum": f"{float(best['num_reads_sum']):.0f}",
        "segment_count": str(len(best["segment_hits"])),
    }


def _run_bandavirus_reference_typing(
    records: list[dict[str, object]],
    out_dir: Path,
    label_key: str,
    segment_key: str,
    query_fasta: Path | None = None,
    single_fastq: str = "",
    fq1: str = "",
    fq2: str = "",
    long_type: str = "",
    threads: int = 4,
    logf=None,
) -> dict[str, str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    candidate_fasta = out_dir / "candidates.fasta"
    meta_by_id = _write_reference_records(records, candidate_fasta)
    if not meta_by_id:
        return {"label": "", "coverage_sum": "0", "depth_sum": "0", "covered_bases_sum": "0", "num_reads_sum": "0", "segment_count": "0", "summary_path": ""}
    method = "read_coverage"
    if query_fasta is not None and query_fasta.is_file():
        method = "blastn"
        blast_path = out_dir / "typing.blastn.tsv"
        cmd = " ".join(
            [
                shlex.quote(_resolve_hadv_blastn_bin()),
                "-task blastn",
                "-query",
                shlex.quote(str(query_fasta)),
                "-subject",
                shlex.quote(str(candidate_fasta)),
                "-outfmt",
                shlex.quote("6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"),
                "-evalue 1e-20",
                "-max_target_seqs 100",
                "-dust no",
                "-num_threads",
                str(max(1, int(threads or 1))),
                "-out",
                shlex.quote(str(blast_path)),
            ]
        )
        run_command(cmd, logf=logf)
        query_lengths: dict[str, int] = {}
        for record in SeqIO.parse(str(query_fasta), "fasta"):
            query_lengths[str(record.id).strip()] = len(str(record.seq))
        rows: list[dict[str, object]] = []
        if blast_path.is_file():
            with blast_path.open("r", encoding="utf-8", errors="ignore") as handle:
                best_by_subject: dict[str, dict[str, object]] = {}
                for line in handle:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 12:
                        continue
                    qseqid, sseqid, _pident, length, *_rest, bitscore = parts
                    qlen = query_lengths.get(str(qseqid).strip(), 0)
                    try:
                        align_len = float(length)
                        coverage = min(100.0, (align_len / qlen) * 100.0) if qlen else 0.0
                        score = (coverage, float(bitscore), align_len)
                    except ValueError:
                        continue
                    previous = best_by_subject.get(sseqid)
                    if previous is None or score > previous["score"]:
                        best_by_subject[sseqid] = {
                            "score": score,
                            "row": {
                                "reference_name": sseqid,
                                "coverage": coverage,
                                "mean_depth": float(bitscore),
                                "covered_bases": align_len,
                                "num_reads": 1.0,
                            },
                        }
                rows = [item["row"] for item in best_by_subject.values()]
    else:
        rows = _run_multi_reference_coverage(
            candidate_fasta,
            out_dir / "typing",
            single_fastq=single_fastq,
            fq1=fq1,
            fq2=fq2,
            long_type=long_type,
            threads=threads,
            logf=logf,
        )
    summary = _aggregate_bandavirus_rows(rows, meta_by_id, label_key=label_key, segment_key=segment_key)
    summary_path = out_dir / "typing_summary.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["label", "segment_count", "coverage_sum", "depth_sum", "covered_bases_sum", "num_reads_sum", "method"])
        writer.writerow([
            summary.get("label", ""),
            summary.get("segment_count", "0"),
            summary.get("coverage_sum", "0"),
            summary.get("depth_sum", "0"),
            summary.get("covered_bases_sum", "0"),
            summary.get("num_reads_sum", "0"),
            method,
        ])
    summary["method"] = method
    summary["summary_path"] = str(summary_path.resolve())
    return summary


def _choose_bandavirus_segment_reference(
    records: list[dict[str, object]],
    output_dir: Path,
    query_fasta: Path | None = None,
    single_fastq: str = "",
    fq1: str = "",
    fq2: str = "",
    long_type: str = "",
    threads: int = 4,
    logf=None,
) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    candidate_fasta = output_dir / "segment_candidates.fasta"
    meta_by_id = _write_reference_records(records, candidate_fasta)
    if not meta_by_id:
        return {}
    if query_fasta is not None and query_fasta.is_file():
        blast_path = output_dir / "segment_candidates.blastn.tsv"
        cmd = " ".join(
            [
                shlex.quote(_resolve_hadv_blastn_bin()),
                "-task blastn",
                "-query",
                shlex.quote(str(query_fasta)),
                "-subject",
                shlex.quote(str(candidate_fasta)),
                "-outfmt",
                shlex.quote("6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"),
                "-evalue 1e-20",
                "-max_target_seqs 100",
                "-dust no",
                "-num_threads",
                str(max(1, int(threads or 1))),
                "-out",
                shlex.quote(str(blast_path)),
            ]
        )
        run_command(cmd, logf=logf)
        query_lengths = {str(record.id).strip(): len(str(record.seq)) for record in SeqIO.parse(str(query_fasta), "fasta")}
        rows: list[dict[str, object]] = []
        if blast_path.is_file():
            best_by_subject: dict[str, dict[str, object]] = {}
            with blast_path.open("r", encoding="utf-8", errors="ignore") as handle:
                for line in handle:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 12:
                        continue
                    qseqid, sseqid, _pident, length, *_rest, bitscore = parts
                    qlen = query_lengths.get(str(qseqid).strip(), 0)
                    try:
                        align_len = float(length)
                        coverage = min(100.0, (align_len / qlen) * 100.0) if qlen else 0.0
                        score = (coverage, float(bitscore), align_len)
                    except ValueError:
                        continue
                    previous = best_by_subject.get(sseqid)
                    if previous is None or score > previous["score"]:
                        best_by_subject[sseqid] = {
                            "score": score,
                            "row": {
                                "reference_name": sseqid,
                                "coverage": coverage,
                                "mean_depth": float(bitscore),
                                "covered_bases": align_len,
                                "num_reads": 1.0,
                            },
                        }
            rows = [item["row"] for item in best_by_subject.values()]
    else:
        rows = _run_multi_reference_coverage(
            candidate_fasta,
            output_dir / "segment_candidates",
            single_fastq=single_fastq,
            fq1=fq1,
            fq2=fq2,
            long_type=long_type,
            threads=threads,
            logf=logf,
        )
    best: dict[str, object] | None = None
    for row in rows:
        ref_name = str(row.get("reference_name") or "").strip()
        meta = meta_by_id.get(ref_name)
        if meta is None:
            continue
        score = _coverage_row_score(row)
        if best is None or score > best["score"]:
            best = {"score": score, "row": row, "meta": meta}
    return best or {}


def _write_bandavirus_reference_fasta(
    reference_rows: list[dict[str, object]],
    output_fasta: Path,
) -> str:
    with output_fasta.open("w", encoding="utf-8") as handle:
        for row in reference_rows:
            accession = str(row.get("accession") or "").strip()
            segment = str(row.get("segment") or "").strip().upper()
            sequence = str(row.get("sequence") or "").strip().upper()
            if not accession or segment not in {"L", "M", "S"} or not sequence:
                continue
            handle.write(f">{accession}_{segment}\n")
            for index in range(0, len(sequence), 80):
                handle.write(sequence[index:index + 80] + "\n")
    return str(output_fasta.resolve()) if output_fasta.is_file() and output_fasta.stat().st_size > 0 else ""


def _cleanup_bandavirus_reference_outputs(screening_dir: Path) -> None:
    if not screening_dir.is_dir():
        return
    for path in screening_dir.glob("*.reference.fasta"):
        if path.is_file():
            path.unlink()


def _format_bandavirus_segment_groups(segment_groups: dict[str, str]) -> str:
    ordered = []
    for segment in ["L", "M", "S"]:
        group = str(segment_groups.get(segment) or "").strip()
        if group:
            ordered.append(f"{segment}={group}")
    return ";".join(ordered)


def _format_bandavirus_segment_values(segment_values: dict[str, str]) -> str:
    return ";".join(str(segment_values.get(segment) or "-").strip() or "-" for segment in ["L", "M", "S"])


def _parse_bandavirus_segment_assignments(text: str) -> dict[str, str]:
    assignments: dict[str, str] = {}
    for chunk in str(text or "").split(";"):
        item = chunk.strip()
        if "=" not in item:
            continue
        segment, value = item.split("=", 1)
        normalized_segment = str(segment or "").strip().upper()
        normalized_value = str(value or "").strip()
        if normalized_segment in {"L", "M", "S"} and normalized_value:
            assignments[normalized_segment] = normalized_value
    return assignments


def _load_bandavirus_segment_manifest_labels(path: Path, value_key: str) -> tuple[dict[str, str], dict[str, str]]:
    segment_values: dict[str, str] = {}
    accessions: dict[str, str] = {}
    if not path.is_file():
        return segment_values, accessions
    with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if not row:
                continue
            segment = str(row.get("segment") or "").strip().upper()
            value = str(row.get(value_key) or "").strip()
            accession = str(row.get("accession") or "").strip()
            if segment in {"L", "M", "S"} and value:
                segment_values[segment] = value
            if segment in {"L", "M", "S"} and accession:
                accessions[segment] = accession
    return segment_values, accessions


def _run_bandavirus_serotype_typing(pre: str, species: str, final_fasta: Path, consensus_fasta: Path) -> None:
    columns = ["样本名称", "病毒类型", "A_F(LMS)", "CJ(LMS)", "分型结果", "参考命中", "说明"]
    query_fasta = consensus_fasta if consensus_fasta.is_file() and consensus_fasta.stat().st_size > 0 else final_fasta
    if not query_fasta.is_file() or query_fasta.stat().st_size == 0:
        _write_placeholder(pre, columns, [pre, species or "Bandavirus", "-", "-", "-", "-", "未检测到可用于班达病毒分型的组装/consensus 序列"])
        return

    screening_dir = Path(f"{pre}_bandavirus_reference_selection")
    selection_path = screening_dir / "selection.tsv"
    selection: dict[str, str] = {}
    if selection_path.is_file():
        with selection_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
            selection = next(csv.DictReader(handle, delimiter="\t"), None) or {}

    af_segments: dict[str, str] = {}
    af_accessions: dict[str, str] = {}
    segment_manifest_path = Path(str(selection.get("segment_manifest_path") or "").strip()) if selection else Path()
    if not segment_manifest_path.is_file():
        segment_manifest_path = screening_dir / "selected_segments.tsv"
    if segment_manifest_path.is_file():
        af_segments, af_accessions = _load_bandavirus_segment_manifest_labels(segment_manifest_path, "af_group")
    if not af_segments:
        af_segments = _parse_bandavirus_segment_assignments(selection.get("segment_groups", "")) if selection else {}
    if not af_segments:
        af_group = str(selection.get("af_group") or "").strip() if selection else ""
        if af_group:
            af_segments = {segment: af_group for segment in ["L", "M", "S"]}

    cj_segments: dict[str, str] = {}
    cj_accessions: dict[str, str] = {}
    cj_dir = screening_dir / "cj_typing"
    cj_dir.mkdir(parents=True, exist_ok=True)
    cj_manifest_path = cj_dir / "selected_segments.tsv"
    cj_rows: list[list[str]] = []
    for segment in ["L", "M", "S"]:
        best = _choose_bandavirus_segment_reference(
            _build_bandavirus_cj_records(segment=segment),
            cj_dir / f"segment_{segment}",
            query_fasta=query_fasta,
        )
        if not best:
            continue
        meta_entry = best.get("meta") if isinstance(best.get("meta"), dict) else {}
        meta = meta_entry.get("meta") if isinstance(meta_entry.get("meta"), dict) else {}
        genotype = str(meta.get("CJ_Genotype") or "").strip()
        accession = str(meta_entry.get("accession") or "").strip()
        score_row = best.get("row") if isinstance(best.get("row"), dict) else {}
        if genotype:
            cj_segments[segment] = genotype
        if accession:
            cj_accessions[segment] = accession
        cj_rows.append(
            [
                segment,
                genotype,
                accession,
                f"{float(score_row.get('coverage') or 0.0):.2f}",
                f"{float(score_row.get('mean_depth') or 0.0):.2f}",
                f"{float(score_row.get('covered_bases') or 0.0):.0f}",
                str(meta.get("CJ_Viral_Strain") or ""),
            ]
        )
    if cj_rows:
        with cj_manifest_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["segment", "cj_group", "accession", "coverage", "mean_depth", "covered_bases", "strain"])
            writer.writerows(cj_rows)

    broad_type = str(selection.get("broad_type") or "").strip() if selection else ""
    af_group = str(selection.get("af_group") or "").strip() if selection else ""
    reassortment_flag = str(selection.get("reassortment_flag") or "no").strip().lower() if selection else "no"
    note_parts: list[str] = []
    selection_note = str(selection.get("note") or "").strip() if selection else ""
    if selection_note:
        note_parts.append(selection_note)
    if broad_type:
        note_parts.append(f"大亚型: {broad_type}")
    if cj_rows:
        note_parts.append("CJ 分型基于 L/M/S 三段分别与 grouped bandavirus 参考比对得到")

    typing_result = broad_type or species or "Bandavirus"
    if af_group:
        typing_result = f"{typing_result} / A_F={af_group}"
    if reassortment_flag == "yes":
        typing_result = f"{typing_result} / 疑似重组或重配"

    af_reference_text = _format_bandavirus_segment_values(af_accessions)
    cj_reference_text = _format_bandavirus_segment_values(cj_accessions)
    reference_text = f"A_F参考={af_reference_text}; CJ参考={cj_reference_text}"

    _write_placeholder(
        pre,
        columns,
        [
            pre,
            _bandavirus_species_label(broad_type, species),
            _format_bandavirus_segment_values(af_segments),
            _format_bandavirus_segment_values(cj_segments),
            typing_result,
            reference_text,
            "；".join(part for part in note_parts if part),
        ],
    )


def resolve_orthohantavirus_reference(
    pre: str,
    species: str = "",
    requested_ref: str = "",
    single_fastq: str = "",
    fq1: str = "",
    fq2: str = "",
    long_type: str = "",
    threads: int = 4,
    query_fasta: str = "",
    logf=None,
) -> dict[str, str]:
    _clear_sample_skip_flag(pre)
    requested = str(requested_ref or "").strip()
    if requested and Path(requested).is_file():
        return {
            "status": "ready",
            "predicted_type": "",
            "s_segment_type": "",
            "selected_segments": "",
            "segment_count": "0",
            "species_label": species or "Orthohantavirus",
            "reference_path": str(Path(requested).resolve()),
            "gff_path": "nogtf",
            "summary_path": "",
            "typing_summary_path": "",
            "segment_manifest_path": "",
        }
    screening_dir = Path(f"{pre}_orthohantavirus_reference_selection")
    screening_dir.mkdir(exist_ok=True)
    _cleanup_bandavirus_reference_outputs(screening_dir)
    query_path = Path(query_fasta) if str(query_fasta or "").strip() else None
    typing_result = _run_bandavirus_reference_typing(
        _build_orthohantavirus_typing_records(),
        screening_dir / "broad_typing",
        label_key="abbrev_primary",
        segment_key="segment",
        query_fasta=query_path,
        single_fastq=single_fastq,
        fq1=fq1,
        fq2=fq2,
        long_type=long_type,
        threads=threads,
        logf=logf,
    )
    predicted_type = str(typing_result.get("label") or "").strip().upper()
    summary_path = screening_dir / "selection.tsv"
    low_support = _is_low_support_orthohantavirus_broad_result(typing_result)
    if not predicted_type or low_support:
        status_text = "skipped" if low_support else "missing"
        note_text = (
            _low_support_orthohantavirus_broad_note(typing_result, species)
            if low_support
            else f"未能确定汉坦病毒型别，species={species or '-'}"
        )
        with summary_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sample", "predicted_type", "s_segment_type", "selected_segments", "segment_count", "status", "note", "reference_path", "typing_summary_path", "segment_manifest_path"])
            writer.writerow([pre, "", "", "", "0", status_text, note_text, "", str(typing_result.get("summary_path") or ""), ""])
        return {
            "status": status_text,
            "predicted_type": "",
            "s_segment_type": "",
            "selected_segments": "",
            "segment_count": "0",
            "species_label": species or "Orthohantavirus",
            "reference_path": "",
            "gff_path": "nogtf",
            "summary_path": str(summary_path.resolve()),
            "typing_summary_path": str(typing_result.get("summary_path") or ""),
            "segment_manifest_path": "",
            "note": note_text,
        }

    selected_reference_rows: list[dict[str, object]] = []
    segment_manifest_rows: list[list[str]] = []
    selected_segments: list[str] = []
    s_segment_type = ""
    for segment_label in ["L", "M", "S"]:
        best = _choose_bandavirus_segment_reference(
            [record for record in _build_orthohantavirus_typing_records(segment=segment_label) if str(record.get("meta", {}).get("abbrev_primary") or "").strip().upper() == predicted_type],
            screening_dir / f"segment_{segment_label}",
            query_fasta=query_path,
            single_fastq=single_fastq,
            fq1=fq1,
            fq2=fq2,
            long_type=long_type,
            threads=threads,
            logf=logf,
        )
        if not best:
            continue
        meta_entry = best.get("meta") if isinstance(best.get("meta"), dict) else {}
        meta = meta_entry.get("meta") if isinstance(meta_entry.get("meta"), dict) else {}
        accession = str(meta_entry.get("accession") or "").strip()
        score_row = best.get("row") if isinstance(best.get("row"), dict) else {}
        selected_reference_rows.append(
            {
                "accession": accession,
                "segment": segment_label,
                "sequence": str(meta_entry.get("sequence") or ""),
            }
        )
        selected_segments.append(segment_label)
        if segment_label == "S":
            s_segment_type = str(meta.get("abbrev_primary") or predicted_type).strip().upper()
        segment_manifest_rows.append(
            [
                segment_label,
                str(meta.get("abbrev_primary") or predicted_type).strip().upper(),
                accession,
                f"{float(score_row.get('coverage') or 0.0):.2f}",
                f"{float(score_row.get('mean_depth') or 0.0):.2f}",
                f"{float(score_row.get('covered_bases') or 0.0):.0f}",
                f"{float(score_row.get('num_reads') or 0.0):.0f}",
                str(meta.get("header") or accession),
            ]
        )
    selected_segments_text = ",".join(selected_segments)
    reference_basename = f"{_sanitize_tree_label(predicted_type)}.reference.fasta"
    reference_path = _write_bandavirus_reference_fasta(selected_reference_rows, screening_dir / reference_basename)
    segment_manifest_path = screening_dir / "selected_segments.tsv"
    with segment_manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["segment", "typed_label", "accession", "coverage", "mean_depth", "covered_bases", "num_reads", "header"])
        writer.writerows(segment_manifest_rows)
    note = ""
    if len(selected_segments) < 3:
        note = f"仅为该型别找到 {len(selected_segments)} 个可用片段（{selected_segments_text or '-' }）"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "predicted_type", "s_segment_type", "selected_segments", "segment_count", "status", "note", "reference_path", "typing_summary_path", "segment_manifest_path"])
        writer.writerow([
            pre,
            predicted_type,
            s_segment_type or predicted_type,
            selected_segments_text,
            str(len(selected_segments)),
            "ready" if reference_path else "missing",
            note,
            reference_path,
            str(typing_result.get("summary_path") or ""),
            str(segment_manifest_path.resolve()),
        ])
    return {
        "status": "ready" if reference_path else "missing",
        "predicted_type": predicted_type,
        "s_segment_type": s_segment_type or predicted_type,
        "selected_segments": selected_segments_text,
        "segment_count": str(len(selected_segments)),
        "species_label": _orthohantavirus_species_label(predicted_type, species),
        "reference_path": reference_path,
        "gff_path": "nogtf",
        "summary_path": str(summary_path.resolve()),
        "typing_summary_path": str(typing_result.get("summary_path") or ""),
        "segment_manifest_path": str(segment_manifest_path.resolve()),
    }


def _load_orthohantavirus_segment_manifest() -> list[dict[str, str]]:
    manifest_path = _resolve_orthohantavirus_typing_manifest()
    if not manifest_path.is_file():
        return []
    rows: list[dict[str, str]] = []
    with manifest_path.open("r", encoding="utf-8-sig", errors="ignore", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row:
                rows.append({str(key): str(value or "").strip() for key, value in row.items()})
    return rows


def _orthohantavirus_manifest_lookup() -> dict[tuple[str, str], dict[str, str]]:
    lookup: dict[tuple[str, str], dict[str, str]] = {}
    for row in _load_orthohantavirus_segment_manifest():
        accession = str(row.get("accession") or row.get("accession_root") or "").strip().upper()
        segment = str(row.get("segment") or "").strip().upper()
        if accession and segment:
            lookup[(accession, segment)] = row
    return lookup


def _load_fasta_record_ids(fasta_path: Path) -> list[str]:
    if not fasta_path.is_file() or fasta_path.stat().st_size == 0:
        return []
    record_ids: list[str] = []
    with fasta_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            record_id = str(line[1:].strip().split()[0] if line[1:].strip() else "").strip()
            if record_id:
                record_ids.append(record_id)
    return record_ids


def _match_orthohantavirus_reference_seqid(record_ids: list[str], accession: str, segment: str, fallback: str) -> str:
    normalized_accession = str(accession or "").strip().upper()
    normalized_segment = str(segment or "").strip().upper()
    candidates: list[str] = []
    fallback_text = str(fallback or "").strip()
    if fallback_text:
        candidates.extend(
            [
                fallback_text,
                fallback_text.split()[0],
                fallback_text.replace("|", " ").split()[0],
            ]
        )
    if normalized_accession and normalized_segment:
        candidates.extend(
            [
                f"{normalized_accession}_{normalized_segment}",
                f"{normalized_segment}_{normalized_accession}",
            ]
        )
    normalized_records = [str(item or "").strip() for item in record_ids if str(item or "").strip()]
    for candidate in candidates:
        candidate_text = str(candidate or "").strip()
        if candidate_text and candidate_text in normalized_records:
            return candidate_text
    for record_id in normalized_records:
        upper_id = record_id.upper()
        if normalized_accession and normalized_accession in upper_id and normalized_segment and normalized_segment in upper_id:
            return record_id
    return fallback_text or (normalized_records[0] if normalized_records else "")


def _build_orthohantavirus_combined_gff(pre: str, reference_fasta: Path, logf=None) -> Path | None:
    screening_dir = Path(f"{pre}_orthohantavirus_reference_selection")
    segment_manifest_path = screening_dir / "selected_segments.tsv"
    if not segment_manifest_path.is_file() or segment_manifest_path.stat().st_size == 0:
        return None
    with segment_manifest_path.open("r", encoding="utf-8-sig", errors="ignore", newline="") as handle:
        segment_rows = [dict(row) for row in csv.DictReader(handle, delimiter="\t") if row]
    if not segment_rows:
        return None
    manifest_lookup = _orthohantavirus_manifest_lookup()
    record_ids = _load_fasta_record_ids(reference_fasta)
    if not record_ids:
        return None
    output_gff = screening_dir / "snpeff_reference.gff3"
    wrote_any = False
    with output_gff.open("w", encoding="utf-8") as out_handle:
        out_handle.write("##gff-version 3\n")
        for row in segment_rows:
            accession = str(row.get("accession") or "").strip().upper()
            segment = str(row.get("segment") or "").strip().upper()
            if not accession or segment not in {"L", "M", "S"}:
                continue
            manifest_row = manifest_lookup.get((accession, segment))
            if not manifest_row:
                continue
            gff_path = Path(str(manifest_row.get("gff_path") or "").strip())
            if not gff_path.is_file() or gff_path.stat().st_size == 0:
                continue
            record_id = _match_orthohantavirus_reference_seqid(
                record_ids,
                accession,
                segment,
                str(row.get("header") or manifest_row.get("record_id") or ""),
            )
            if not record_id:
                continue
            try:
                with gff_path.open("r", encoding="utf-8", errors="ignore") as in_handle:
                    for line in in_handle:
                        if not line:
                            continue
                        if line.startswith("##gff-version"):
                            continue
                        if line.startswith("##sequence-region"):
                            parts = line.rstrip("\n").split()
                            if len(parts) >= 4:
                                out_handle.write(f"##sequence-region {record_id} {parts[-2]} {parts[-1]}\n")
                            continue
                        if line.startswith("#"):
                            continue
                        fields = line.rstrip("\n").split("\t")
                        if len(fields) < 9:
                            continue
                        fields[0] = record_id
                        out_handle.write("\t".join(fields) + "\n")
                        wrote_any = True
            except OSError as exc:
                if logf is not None:
                    logf.write(f"[WARN] failed to read orthohantavirus gff {gff_path}: {exc}\n")
                    logf.flush()
                continue
    if not wrote_any:
        try:
            output_gff.unlink(missing_ok=True)
        except OSError:
            pass
        return None
    return output_gff


def run_orthohantavirus_consensus_typing(pre: str, consensus_fasta: Path, logf=None) -> dict[str, str]:
    if not consensus_fasta.is_file() or consensus_fasta.stat().st_size == 0:
        return {"status": "missing", "summary_path": "", "predicted_type": "", "s_segment_type": "", "selected_segments": "", "segment_count": "0"}
    selection_path = Path(f"{pre}_orthohantavirus_reference_selection") / "selection.tsv"
    selection: dict[str, str] = {}
    if selection_path.is_file():
        with selection_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
            selection = next(csv.DictReader(handle, delimiter="\t"), None) or {}
    if not selection:
        selection = resolve_orthohantavirus_reference(
            f"{pre}_consensus",
            species="Orthohantavirus",
            requested_ref="",
            query_fasta=str(consensus_fasta),
            logf=logf,
        )
    summary_path = Path(f"{pre}_orthohantavirus_reference_selection") / "consensus_typing.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "predicted_type", "s_segment_type", "selected_segments", "segment_count", "status", "selection_summary_path", "source"])
        writer.writerow([
            pre,
            selection.get("predicted_type", ""),
            selection.get("s_segment_type", ""),
            selection.get("selected_segments", ""),
            selection.get("segment_count", "0"),
            selection.get("status", "missing"),
            selection.get("summary_path", "") or str(selection_path.resolve() if selection_path.is_file() else ""),
            "read_selection" if selection_path.is_file() else "consensus_fallback",
        ])
    return {
        "status": str(selection.get("status") or ("ready" if str(selection.get("predicted_type") or "").strip() else "missing")),
        "summary_path": str(summary_path.resolve()),
        "predicted_type": str(selection.get("predicted_type") or ""),
        "s_segment_type": str(selection.get("s_segment_type") or ""),
        "selected_segments": str(selection.get("selected_segments") or ""),
        "segment_count": str(selection.get("segment_count") or "0"),
    }


def _orthoebolavirus_species_label(abbrev: str = "", fallback: str = "") -> str:
    normalized = str(abbrev or "").strip().upper()
    for row in _load_rhinovirus_manifest(_resolve_orthoebolavirus_typing_manifest()):
        candidate = str(row.get("abbrev") or "").strip().upper()
        if candidate == normalized:
            species_label = str(row.get("species") or "").strip()
            virus_name = str(row.get("virus_name") or "").strip()
            if species_label and virus_name:
                return f"{virus_name} ({species_label})"
            if species_label:
                return species_label
            if virus_name:
                return virus_name
    fallback_text = str(fallback or "").strip()
    return fallback_text or "Orthoebolavirus"


def _infer_orthoebolavirus_abbrev_from_label(label: str = "") -> str:
    normalized = str(label or "").strip().lower()
    if not normalized:
        return ""
    rules = [
        ("BOMV", ["bomv", "bombali"]),
        ("BDBV", ["bdbv", "bundibugyo"]),
        ("RESTV", ["restv", "reston"]),
        ("SUDV", ["sudv", "sudan"]),
        ("TAFV", ["tafv", "taï forest", "tai forest", "taiense"]),
        ("EBOV", ["ebov", "ebola virus", "zairense", "zaire ebola", "埃博拉"]),
    ]
    for abbrev, tokens in rules:
        if any(token in normalized for token in tokens):
            return abbrev
    return ""


def _build_orthoebolavirus_typing_records() -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for row in _load_rhinovirus_manifest(_resolve_orthoebolavirus_typing_manifest()):
        fasta_path = Path(str(row.get("fasta_path") or "").strip())
        abbrev = str(row.get("abbrev") or "").strip().upper()
        accession = str(row.get("accession") or "").strip()
        record_id = str(row.get("record_id") or "").strip()
        if not abbrev or not accession or not fasta_path.is_file():
            continue
        record = _get_fasta_record_by_id(fasta_path, record_id or None)
        if record is None:
            record = next(SeqIO.parse(str(fasta_path), "fasta"), None)
        if record is None:
            continue
        sequence = str(record.seq).strip().upper()
        if not sequence:
            continue
        records.append(
            {
                "fasta_id": record_id or f"{abbrev}_{accession.split('.', 1)[0]}",
                "accession": accession,
                "sequence": sequence,
                "reference_length": len(sequence),
                "meta": dict(row),
                "header": str(record.description or record.id).strip(),
                "abbrev": abbrev,
                "source_record_id": str(record.id),
            }
        )
    return records


def resolve_orthoebolavirus_reference(
    pre: str,
    species: str = "",
    requested_ref: str = "",
    single_fastq: str = "",
    fq1: str = "",
    fq2: str = "",
    long_type: str = "",
    threads: int = 4,
    query_fasta: str = "",
    output_dir: Path | None = None,
    logf=None,
) -> dict[str, str]:
    requested = str(requested_ref or "").strip()
    explicit_type = _infer_orthoebolavirus_abbrev_from_label(species)
    if requested and Path(requested).is_file():
        return {
            "status": "ready",
            "predicted_type": explicit_type,
            "species_label": _orthoebolavirus_species_label(explicit_type, species or "Orthoebolavirus"),
            "virus_name": "",
            "accession": "",
            "reference_path": str(Path(requested).resolve()),
            "gff_path": "nogtf",
            "summary_path": "",
            "nextclade_dataset": str(_resolve_ebola_nextclade_dataset()) if explicit_type == "EBOV" else "",
        }

    screening_dir = output_dir if output_dir is not None else Path(f"{pre}_orthoebolavirus_reference_selection")
    screening_dir.mkdir(exist_ok=True)
    final_summary_path = screening_dir / "selection.tsv"
    records = _build_orthoebolavirus_typing_records()
    if not records:
        with final_summary_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sample", "predicted_type", "species_label", "virus_name", "accession", "status", "note", "reference_path", "gff_path", "nextclade_dataset"])
            writer.writerow([pre, explicit_type, _orthoebolavirus_species_label(explicit_type, species), "", "", "missing", "未找到 Orthoebolavirus 本地参考库", "", "nogtf", ""])
        return {
            "status": "missing",
            "predicted_type": explicit_type,
            "species_label": _orthoebolavirus_species_label(explicit_type, species),
            "virus_name": "",
            "accession": "",
            "reference_path": "",
            "gff_path": "nogtf",
            "summary_path": str(final_summary_path.resolve()),
            "nextclade_dataset": "",
        }

    candidate_fasta = screening_dir / "candidate_references.fasta"
    meta_by_id = _write_reference_records(records, candidate_fasta)
    query_path = Path(query_fasta) if str(query_fasta or "").strip() else None
    if query_path is not None and query_path.is_file() and query_path.stat().st_size > 0:
        coverage_rows = _run_multi_reference_blast(query_path, candidate_fasta, screening_dir / "reference_selection", meta_by_id, threads=threads, logf=logf)
        coverage_rows.sort(key=_blast_row_score, reverse=True)
    else:
        coverage_rows = _run_multi_reference_coverage(
            candidate_fasta,
            screening_dir / "reference_selection",
            single_fastq=single_fastq,
            fq1=fq1,
            fq2=fq2,
            long_type=long_type,
            threads=threads,
            logf=logf,
        )
        coverage_rows.sort(key=_coverage_row_score, reverse=True)
    best_row = coverage_rows[0] if coverage_rows else {}
    best_meta = meta_by_id.get(str(best_row.get("reference_name") or "").strip()) if best_row else {}
    meta_row = dict(best_meta.get("meta") or {}) if isinstance(best_meta, dict) else {}
    predicted_type = str(best_meta.get("abbrev") or meta_row.get("abbrev") or explicit_type).strip().upper() if best_meta else explicit_type
    virus_name = str(meta_row.get("virus_name") or "").strip()
    species_label = _orthoebolavirus_species_label(predicted_type, species or "Orthoebolavirus")
    accession = str(best_meta.get("accession") or meta_row.get("accession") or "").strip() if best_meta else ""
    selected_reference_path = ""
    if best_meta:
        selected_reference_path = str((screening_dir / f"{_sanitize_tree_label(predicted_type or accession or 'orthoebolavirus')}.reference.fasta").resolve())
        if not _extract_named_fasta_record(candidate_fasta, str(best_meta.get("fasta_id") or ""), Path(selected_reference_path)):
            selected_reference_path = ""
    selected_gff_path = str(meta_row.get("gff_path") or "nogtf")
    nextclade_dataset = str(_resolve_ebola_nextclade_dataset()) if predicted_type == "EBOV" else ""
    note = "基于 Orthoebolavirus 本地完整/编码完整参考基因组库选择覆盖度最优参考。"
    if predicted_type == "EBOV":
        note = f"{note} EBOV 样本后续同时执行 Ebola Nextclade 分型。"
    with final_summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "predicted_type", "species_label", "virus_name", "isolate", "accession", "coverage", "mean_depth", "covered_bases", "num_reads", "reference_name", "reference_path", "gff_path", "status", "note", "nextclade_dataset"])
        writer.writerow([
            pre,
            predicted_type,
            species_label,
            virus_name,
            str(meta_row.get("isolate") or ""),
            accession,
            f"{float(best_row.get('coverage') or 0.0):.6f}",
            f"{float(best_row.get('mean_depth') or 0.0):.6f}" if str(best_row.get("mean_depth") or "").strip() else "",
            f"{float(best_row.get('covered_bases') or 0.0):.0f}",
            f"{float(best_row.get('num_reads') or 0.0):.0f}" if str(best_row.get("num_reads") or "").strip() else "",
            str(meta_row.get("header") or best_meta.get("header") or accession),
            selected_reference_path,
            selected_gff_path,
            "ready" if selected_reference_path else "missing",
            note if best_meta else "未获得 Orthoebolavirus 参考命中结果",
            nextclade_dataset,
        ])
    return {
        "status": "ready" if selected_reference_path else "missing",
        "predicted_type": predicted_type,
        "species_label": species_label,
        "virus_name": virus_name,
        "accession": accession,
        "reference_path": selected_reference_path,
        "gff_path": selected_gff_path,
        "summary_path": str(final_summary_path.resolve()),
        "nextclade_dataset": nextclade_dataset,
    }


def run_orthoebolavirus_consensus_typing(pre: str, consensus_fasta: Path, logf=None) -> dict[str, str]:
    if not consensus_fasta.is_file() or consensus_fasta.stat().st_size == 0:
        return {"status": "missing", "summary_path": "", "predicted_type": "", "species_label": ""}
    output_dir = Path(f"{pre}_orthoebolavirus_reference_selection") / "consensus_typing"
    selection = resolve_orthoebolavirus_reference(
        pre,
        species="Orthoebolavirus",
        requested_ref="",
        query_fasta=str(consensus_fasta),
        output_dir=output_dir,
        logf=logf,
    )
    summary_path = output_dir / "consensus_typing.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "predicted_type", "species_label", "virus_name", "accession", "selection_summary_path", "nextclade_dataset"])
        writer.writerow([
            pre,
            selection.get("predicted_type", ""),
            selection.get("species_label", ""),
            selection.get("virus_name", ""),
            selection.get("accession", ""),
            selection.get("summary_path", ""),
            selection.get("nextclade_dataset", ""),
        ])
    return {
        "status": "ready",
        "summary_path": str(summary_path.resolve()),
        "predicted_type": str(selection.get("predicted_type") or ""),
        "species_label": str(selection.get("species_label") or ""),
        "nextclade_dataset": str(selection.get("nextclade_dataset") or ""),
    }


def resolve_bandavirus_reference(
    pre: str,
    species: str = "",
    requested_ref: str = "",
    single_fastq: str = "",
    fq1: str = "",
    fq2: str = "",
    long_type: str = "",
    threads: int = 4,
    query_fasta: str = "",
    logf=None,
) -> dict[str, str]:
    requested = str(requested_ref or "").strip()
    if requested and Path(requested).is_file():
        return {
            "status": "ready",
            "broad_type": "SFTSV" if _is_bandavirus(species) and _bandavirus_is_sftsv(species) else "",
            "af_group": "",
            "segment_groups": "",
            "reassortment_flag": "no",
            "species_label": species or "Bandavirus",
            "reference_path": str(Path(requested).resolve()),
            "gff_path": "nogtf",
            "summary_path": "",
            "broad_summary_path": "",
            "subgroup_summary_path": "",
        }
    screening_dir = Path(f"{pre}_bandavirus_reference_selection")
    screening_dir.mkdir(exist_ok=True)
    _cleanup_bandavirus_reference_outputs(screening_dir)
    query_path = Path(query_fasta) if str(query_fasta or "").strip() else None
    broad_result = _run_bandavirus_reference_typing(
        _build_bandavirus_typing_records(),
        screening_dir / "broad_typing",
        label_key="abbrev",
        segment_key="segment",
        query_fasta=query_path,
        single_fastq=single_fastq,
        fq1=fq1,
        fq2=fq2,
        long_type=long_type,
        threads=threads,
        logf=logf,
    )
    broad_type = str(broad_result.get("label") or "").strip().upper()
    summary_path = screening_dir / "selection.tsv"
    if not broad_type:
        with summary_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sample", "broad_type", "af_group", "status", "note", "reference_path", "broad_summary_path", "subgroup_summary_path"])
            writer.writerow([pre, "", "", "missing", f"未能确定班达病毒大亚型，species={species or '-'}", "", str(broad_result.get("summary_path") or ""), ""])
        return {
            "status": "missing",
            "broad_type": "",
            "af_group": "",
            "segment_groups": "",
            "reassortment_flag": "no",
            "species_label": species or "Bandavirus",
            "reference_path": "",
            "gff_path": "nogtf",
            "summary_path": str(summary_path.resolve()),
            "broad_summary_path": str(broad_result.get("summary_path") or ""),
            "subgroup_summary_path": "",
        }
    if not _bandavirus_is_sftsv(broad_type):
        selected_rows: list[dict[str, object]] = []
        for record in _build_bandavirus_typing_records():
            meta = record.get("meta") if isinstance(record.get("meta"), dict) else {}
            if str(meta.get("abbrev") or "").strip().upper() != broad_type:
                continue
            selected_rows.append(
                {
                    "accession": str(record.get("accession") or ""),
                    "segment": str(meta.get("segment") or "").strip().upper(),
                    "sequence": str(record.get("sequence") or ""),
                }
            )
        selected_rows.sort(key=lambda item: {"L": 0, "M": 1, "S": 2}.get(str(item.get("segment") or ""), 9))
        reference_path = _write_bandavirus_reference_fasta(selected_rows, screening_dir / f"{_sanitize_tree_label(broad_type)}.reference.fasta")
        with summary_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sample", "broad_type", "af_group", "segment_groups", "reassortment_flag", "status", "note", "reference_path", "broad_summary_path", "subgroup_summary_path"])
            writer.writerow([pre, broad_type, "", "", "no", "ready" if reference_path else "missing", "", reference_path, str(broad_result.get("summary_path") or ""), ""])
        return {
            "status": "ready" if reference_path else "missing",
            "broad_type": broad_type,
            "af_group": "",
            "segment_groups": "",
            "reassortment_flag": "no",
            "species_label": _bandavirus_species_label(broad_type, species),
            "reference_path": reference_path,
            "gff_path": "nogtf",
            "summary_path": str(summary_path.resolve()),
            "broad_summary_path": str(broad_result.get("summary_path") or ""),
            "subgroup_summary_path": "",
        }
    subgroup_result = _run_bandavirus_reference_typing(
        _build_bandavirus_af_records(),
        screening_dir / "subgroup_typing",
        label_key="Tree_Group",
        segment_key="Tree_Segment",
        query_fasta=query_path,
        single_fastq=single_fastq,
        fq1=fq1,
        fq2=fq2,
        long_type=long_type,
        threads=threads,
        logf=logf,
    )
    af_group = str(subgroup_result.get("label") or "").strip()
    selected_reference_rows: list[dict[str, object]] = []
    segment_manifest_rows: list[list[str]] = []
    segment_best_groups: dict[str, str] = {}
    for segment in ["L", "M", "S"]:
        best = _choose_bandavirus_segment_reference(
            _build_bandavirus_af_records(segment=segment),
            screening_dir / f"segment_{segment}",
            query_fasta=query_path,
            single_fastq=single_fastq,
            fq1=fq1,
            fq2=fq2,
            long_type=long_type,
            threads=threads,
            logf=logf,
        )
        if not best:
            continue
        meta_entry = best.get("meta") if isinstance(best.get("meta"), dict) else {}
        meta = meta_entry.get("meta") if isinstance(meta_entry.get("meta"), dict) else {}
        segment_group = str(meta.get("Tree_Group") or "").strip()
        if segment_group:
            segment_best_groups[segment] = segment_group
        selected_reference_rows.append(
            {
                "accession": str(meta_entry.get("accession") or ""),
                "segment": segment,
                "sequence": str(meta_entry.get("sequence") or ""),
            }
        )
        score_row = best.get("row") if isinstance(best.get("row"), dict) else {}
        segment_manifest_rows.append(
            [
                segment,
                segment_group,
                str(meta_entry.get("accession") or ""),
                f"{float(score_row.get('coverage') or 0.0):.2f}",
                f"{float(score_row.get('mean_depth') or 0.0):.2f}",
                f"{float(score_row.get('covered_bases') or 0.0):.0f}",
                f"{float(score_row.get('num_reads') or 0.0):.0f}",
                str(meta.get("Tree_Sample_ID") or ""),
            ]
        )
    unique_groups = sorted({group for group in segment_best_groups.values() if group})
    reassortment_flag = "yes" if len(unique_groups) > 1 else "no"
    final_af_group = af_group
    if reassortment_flag == "yes":
        final_af_group = ""
    elif not final_af_group and len(unique_groups) == 1:
        final_af_group = unique_groups[0]
    segment_group_text = _format_bandavirus_segment_groups(segment_best_groups)
    note_text = ""
    if reassortment_flag == "yes":
        note_text = f"L/M/S 最优分型不一致（{segment_group_text}），提示疑似重组/重配。"
    reference_basename = f"SFTSV_{_sanitize_tree_label(final_af_group or 'reassortant')}.reference.fasta"
    reference_path = _write_bandavirus_reference_fasta(selected_reference_rows, screening_dir / reference_basename)
    segment_manifest_path = screening_dir / "selected_segments.tsv"
    with segment_manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["segment", "af_group", "accession", "coverage", "mean_depth", "covered_bases", "num_reads", "sample_id"])
        writer.writerows(segment_manifest_rows)
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "broad_type", "af_group", "segment_groups", "reassortment_flag", "status", "note", "reference_path", "broad_summary_path", "subgroup_summary_path", "segment_manifest_path"])
        writer.writerow([
            pre,
            broad_type,
            final_af_group,
            segment_group_text,
            reassortment_flag,
            "ready" if len(selected_reference_rows) == 3 and reference_path else "missing",
            note_text if note_text else ("" if final_af_group else "未能确定 SFTSV A_F 子亚型"),
            reference_path,
            str(broad_result.get("summary_path") or ""),
            str(subgroup_result.get("summary_path") or ""),
            str(segment_manifest_path.resolve()),
        ])
    return {
        "status": "ready" if len(selected_reference_rows) == 3 and reference_path else "missing",
        "broad_type": broad_type,
        "af_group": final_af_group,
        "segment_groups": segment_group_text,
        "reassortment_flag": reassortment_flag,
        "species_label": "Bandavirus dabieense",
        "reference_path": reference_path,
        "gff_path": "nogtf",
        "summary_path": str(summary_path.resolve()),
        "broad_summary_path": str(broad_result.get("summary_path") or ""),
        "subgroup_summary_path": str(subgroup_result.get("summary_path") or ""),
    }


def run_bandavirus_consensus_typing(pre: str, consensus_fasta: Path, logf=None) -> dict[str, str]:
    if not consensus_fasta.is_file() or consensus_fasta.stat().st_size == 0:
        return {"status": "missing", "summary_path": "", "broad_type": "", "af_group": "", "segment_groups": "", "reassortment_flag": "no"}
    selection_path = Path(f"{pre}_bandavirus_reference_selection") / "selection.tsv"
    selection: dict[str, str] = {}
    if selection_path.is_file():
        with selection_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
            selection = next(csv.DictReader(handle, delimiter="\t"), None) or {}
    if not selection:
        selection = resolve_bandavirus_reference(
            f"{pre}_consensus",
            species="Bandavirus",
            requested_ref="",
            query_fasta=str(consensus_fasta),
            logf=logf,
        )
    summary_path = Path(f"{pre}_bandavirus_reference_selection") / "consensus_typing.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "broad_type", "af_group", "segment_groups", "reassortment_flag", "selection_summary_path", "source"])
        writer.writerow([
            pre,
            selection.get("broad_type", ""),
            selection.get("af_group", ""),
            selection.get("segment_groups", ""),
            selection.get("reassortment_flag", "no"),
            selection.get("summary_path", "") or str(selection_path.resolve() if selection_path.is_file() else ""),
            "read_selection" if selection_path.is_file() else "consensus_fallback",
        ])
    return {
        "status": "ready" if str(selection.get("broad_type") or "").strip() else "missing",
        "summary_path": str(summary_path.resolve()),
        "broad_type": str(selection.get("broad_type") or ""),
        "af_group": str(selection.get("af_group") or ""),
        "segment_groups": str(selection.get("segment_groups") or ""),
        "reassortment_flag": str(selection.get("reassortment_flag") or "no"),
    }


def _run_orthohantavirus_serotype_typing(pre: str, species: str, final_fasta: Path, consensus_fasta: Path, logf=None) -> None:
    columns = ["样本名称", "病毒类型", "S分型", "LMS参考", "分型结果", "参考命中", "说明"]
    query_fasta = consensus_fasta if consensus_fasta.is_file() and consensus_fasta.stat().st_size > 0 else final_fasta
    if not query_fasta.is_file() or query_fasta.stat().st_size == 0:
        _write_placeholder(pre, columns, [pre, species or "Orthohantavirus", "-", "-", "-", "-", "未检测到可用于汉坦病毒分型的组装/consensus 序列"])
        return
    screening_dir = Path(f"{pre}_orthohantavirus_reference_selection")
    selection_path = screening_dir / "selection.tsv"
    selection: dict[str, str] = {}
    if selection_path.is_file():
        with selection_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
            selection = next(csv.DictReader(handle, delimiter="\t"), None) or {}
    if not selection:
        selection = resolve_orthohantavirus_reference(
            pre,
            species=species,
            requested_ref="",
            query_fasta=str(query_fasta),
            logf=None,
        )
    segment_columns: list[str] = []
    segment_rows: list[dict[str, str]] = []
    segment_manifest_path = screening_dir / "selected_segments.tsv"
    if segment_manifest_path.is_file():
        with segment_manifest_path.open("r", encoding="utf-8-sig", errors="ignore", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            segment_columns = list(reader.fieldnames or [])
            segment_rows = [dict(row) for row in reader if row]
    segment_accessions: dict[str, str] = {}
    for row in segment_rows:
        segment = str(row.get("segment") or "").strip().upper()
        accession = str(row.get("accession") or "").strip()
        if segment in {"L", "M", "S"} and accession:
            segment_accessions[segment] = accession
    selected_type = str(selection.get("predicted_type") or "").strip().upper() if selection else ""
    s_segment_type = str(selection.get("s_segment_type") or "").strip().upper() if selection else ""
    selected_segments = str(selection.get("selected_segments") or "").strip() if selection else ""
    segment_count = str(selection.get("segment_count") or "").strip() if selection else ""
    selection_status = str(selection.get("status") or "").strip().lower() if selection else ""
    note_parts: list[str] = []
    selection_note = str(selection.get("note") or "").strip() if selection else ""
    if selection_note:
        note_parts.append(selection_note)
    if selection_status in {"missing", "skipped"} or not (selected_type or s_segment_type):
        _write_placeholder(
            pre,
            columns,
            [
                pre,
                species or "Orthohantavirus",
                "-",
                "-",
                "非汉坦病毒/证据不足",
                "-",
                "；".join(part for part in note_parts if part) or "汉坦病毒 broad 分型支持不足，不继续后续分型。",
            ],
        )
        return
    if selected_segments:
        note_parts.append(f"可用参考片段：{selected_segments}")
    if segment_count:
        note_parts.append(f"片段数：{segment_count}")
    reference_fasta = Path("genomes") / "ref.fa"
    filtered_vcf = Path("snps.filt1.vcf")
    screening_dir = Path(f"{pre}_orthohantavirus_reference_selection")
    snpeff_note = ""
    try:
        combined_gff = _build_orthohantavirus_combined_gff(pre, reference_fasta)
        if combined_gff is not None and reference_fasta.is_file() and filtered_vcf.is_file():
            snpeff_result = _run_influenza_snpeff_annotation(
                pre,
                reference_fasta,
                combined_gff,
                filtered_vcf,
                screening_dir,
                logf=logf,
            )
            if snpeff_result.get("status") == "ready":
                snpeff_note = "已基于对应参考 GFF 生成汉坦病毒突变注释表。"
            elif snpeff_result.get("status") == "failed":
                snpeff_note = str(snpeff_result.get("note") or "").strip()
        elif not reference_fasta.is_file():
            snpeff_note = "未找到 genomes/ref.fa，无法生成汉坦病毒突变注释表。"
        elif not filtered_vcf.is_file():
            snpeff_note = "未找到 snps.filt1.vcf，无法生成汉坦病毒突变注释表。"
        else:
            snpeff_note = "未能为当前汉坦参考组合生成可用 GFF，无法继续突变注释。"
    except Exception as exc:
        snpeff_note = f"汉坦病毒突变注释失败: {exc}"
    if snpeff_note:
        note_parts.append(snpeff_note)
    typing_result = _orthohantavirus_species_label(selected_type or s_segment_type, species)
    reference_text = _format_bandavirus_segment_values(segment_accessions)
    _write_placeholder(
        pre,
        columns,
        [
            pre,
            typing_result or species or "Orthohantavirus",
            s_segment_type or selected_type or "-",
            selected_segments or "-",
            selected_type or s_segment_type or "-",
            reference_text or "-",
            "；".join(part for part in note_parts if part),
        ],
    )


def _resolve_astroviridae_db_dir() -> Path:
    env_root = str(os.environ.get("META_ASTROVIRIDAE_DB_DIR") or "").strip()
    if env_root:
        return Path(env_root).expanduser().resolve()
    database_root = str(os.environ.get("META_DATABASE_ROOT") or "").strip()
    project_root = _project_root()
    candidates: list[Path] = []
    if database_root:
        candidates.append((Path(database_root) / "virus" / "astroviridae").expanduser())
    candidates.extend(
        [
            project_root / "database" / "virus" / "astroviridae",
            Path("/data/deploy/meta_genome/database/virus/astroviridae"),
        ]
    )
    for candidate in candidates:
        if candidate.is_dir():
            return candidate.resolve()
    return candidates[0].resolve()


def _resolve_astroviridae_reference_orf2_manifest() -> Path:
    return (_resolve_astroviridae_db_dir() / "reference_genomes" / "astroviridae_orf2_capsid_references.tsv").resolve()


def _resolve_astroviridae_reference_orf2_fasta() -> Path:
    return (_resolve_astroviridae_db_dir() / "reference_genomes" / "astroviridae_orf2_capsid_references.fasta").resolve()


def _resolve_astroviridae_typing_reference_manifest() -> Path:
    return (_resolve_astroviridae_db_dir() / "reference_genomes" / "astroviridae_typing_reference_genomes_manifest.tsv").resolve()


def _resolve_astroviridae_typed_db_summary() -> Path:
    return (_resolve_astroviridae_db_dir() / "astroviridae_orf2_typing" / "typing_summary.tsv").resolve()


def _resolve_astroviridae_db_fasta() -> Path:
    return (_resolve_astroviridae_db_dir() / "astro_db.fasta").resolve()


def _resolve_astroviridae_db_metadata() -> Path:
    return (_resolve_astroviridae_db_dir() / "astro_db.csv").resolve()


def _normalize_astro_subtype(value: str) -> str:
    text = str(value or "").strip().upper()
    text = re.sub(r"\s*;\s*", ";", text)
    text = re.sub(r"\s+", " ", text)
    return text


def _extract_astro_genus(value: str) -> str:
    normalized = _normalize_species_label(value)
    if "mamastrovirus" in normalized:
        return "Mamastrovirus"
    if "avastrovirus" in normalized:
        return "Avastrovirus"
    return ""


def _astro_species_label(genus: str, species: str = "") -> str:
    species_text = str(species or "").strip()
    if species_text:
        return species_text
    genus_text = str(genus or "").strip()
    if genus_text:
        return genus_text
    return "Astroviridae"


def _load_astroviridae_reference_orf2_manifest() -> list[dict[str, str]]:
    manifest_path = _resolve_astroviridae_reference_orf2_manifest()
    if not manifest_path.is_file():
        return []
    return _load_rhinovirus_manifest(manifest_path)


def _load_astroviridae_typed_db_metadata() -> dict[str, dict[str, str]]:
    path = _resolve_astroviridae_db_metadata()
    if not path.is_file():
        return {}
    rows: dict[str, dict[str, str]] = {}
    with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            accession = str(row.get("Accession") or "").strip()
            if accession:
                rows[accession.split(".", 1)[0]] = {str(key): str(value or "") for key, value in row.items()}
    return rows


def _load_astroviridae_candidate_reference_manifest() -> list[dict[str, str]]:
    combined_by_accession: dict[str, dict[str, str]] = {}
    for row in _load_astroviridae_reference_orf2_manifest():
        accession = str(row.get("accession") or row.get("accession_full") or "").split(".", 1)[0]
        if accession:
            combined_by_accession[accession] = dict(row)

    typed_db_summary = _resolve_astroviridae_typed_db_summary()
    typed_db_fasta = _resolve_astroviridae_db_fasta()
    typed_db_meta = _load_astroviridae_typed_db_metadata()
    if typed_db_summary.is_file() and typed_db_fasta.is_file():
        for row in _load_rhinovirus_manifest(typed_db_summary):
            status = str(row.get("status") or "").strip().lower()
            accession_full = str(row.get("query_id") or "").strip()
            accession = accession_full.split(".", 1)[0] if accession_full else ""
            subtype = _normalize_astro_subtype(str(row.get("subtype") or ""))
            genus = str(row.get("genus") or "").strip()
            if status not in {"typed", "review"}:
                continue
            if not accession or not subtype or genus not in {"Mamastrovirus", "Avastrovirus"}:
                continue
            meta_row = typed_db_meta.get(accession, {})
            completeness = _normalize_species_label(str(meta_row.get("Nuc_Completeness") or ""))
            if completeness and "complete" not in completeness:
                continue
            merged = dict(combined_by_accession.get(accession, {}))
            merged.update(
                {
                    "accession": accession,
                    "accession_full": accession_full or accession,
                    "abbrev": subtype,
                    "genus": genus,
                    "species": str(row.get("species") or meta_row.get("Species") or "").strip(),
                    "virus_name": str(row.get("virus_name") or meta_row.get("Organism_Name") or "").strip(),
                    "isolate": str(meta_row.get("Isolate") or "").strip(),
                    "header": accession_full or accession,
                    "fasta_path": str(typed_db_fasta),
                    "gff_path": "",
                    "available_sequence": str(meta_row.get("Nuc_Completeness") or "complete").strip() or "complete",
                    "source_record_id": accession_full or accession,
                    "evidence": "astro_db_orf2_typing",
                    "typing_status": status,
                }
            )
            combined_by_accession[accession] = merged
    return sorted(
        combined_by_accession.values(),
        key=lambda item: (
            str(item.get("genus") or "").strip(),
            _normalize_astro_subtype(str(item.get("abbrev") or "")),
            str(item.get("accession") or item.get("accession_full") or "").strip(),
        ),
    )


def _build_astro_subject_meta() -> dict[str, dict[str, str]]:
    subject_meta: dict[str, dict[str, str]] = {}
    for row in _load_astroviridae_reference_orf2_manifest():
        accession = str(row.get("accession") or row.get("accession_full") or "").strip()
        accession_root = accession.split(".", 1)[0]
        header = str(row.get("header") or "").strip()
        meta = {
            "genus": str(row.get("genus") or "").strip(),
            "species": str(row.get("species") or "").strip(),
            "normalized_type": _normalize_astro_subtype(str(row.get("abbrev") or "")),
            "abbrev": str(row.get("abbrev") or "").strip(),
            "accession": accession,
            "accession_root": accession_root,
            "gff_path": str(row.get("gff_path") or "").strip(),
            "fasta_path": str(row.get("fasta_path") or "").strip(),
        }
        for key in {accession, accession_root, header, header.split()[0].split(".", 1)[0] if header else ""}:
            if key:
                subject_meta[key] = meta
    return subject_meta


def _score_astro_type_rows(rows: list[dict[str, object]]) -> dict[str, str]:
    by_type: dict[str, dict[str, object]] = {}
    for row in rows:
        meta = row.get("_meta") if isinstance(row.get("_meta"), dict) else {}
        normalized_type = _normalize_astro_subtype(str(meta.get("normalized_type") or ""))
        if not normalized_type:
            continue
        score = (
            float(row.get("coverage") or 0.0),
            float(row.get("mean_depth") or 0.0),
            float(row.get("covered_bases") or 0.0),
            float(row.get("num_reads") or 0.0),
        )
        bucket = by_type.setdefault(
            normalized_type,
            {
                "normalized_type": normalized_type,
                "genus": str(meta.get("genus") or ""),
                "species": str(meta.get("species") or ""),
                "subject": str(row.get("reference_name") or ""),
                "coverage": score[0],
                "mean_depth": score[1],
                "covered_bases": score[2],
                "num_reads": score[3],
                "reference_count": set(),
            },
        )
        if score > (
            float(bucket["coverage"]),
            float(bucket["mean_depth"]),
            float(bucket["covered_bases"]),
            float(bucket["num_reads"]),
        ):
            bucket["subject"] = str(row.get("reference_name") or "")
            bucket["coverage"] = score[0]
            bucket["mean_depth"] = score[1]
            bucket["covered_bases"] = score[2]
            bucket["num_reads"] = score[3]
            bucket["genus"] = str(meta.get("genus") or "")
            bucket["species"] = str(meta.get("species") or "")
        bucket["reference_count"].add(str(meta.get("accession_root") or str(row.get("reference_name") or "").split("_", 1)[0]))
    if not by_type:
        return {
            "type": "",
            "genus": "",
            "species": "",
            "subject": "",
            "coverage": "",
            "mean_depth": "",
            "covered_bases": "",
            "num_reads": "",
            "reference_count": "0",
        }
    ranked = sorted(
        by_type.values(),
        key=lambda item: (
            float(item["coverage"]),
            float(item["mean_depth"]),
            float(item["covered_bases"]),
            float(item["num_reads"]),
            len(item["reference_count"]),
        ),
        reverse=True,
    )
    best = ranked[0]
    return {
        "type": str(best["normalized_type"]),
        "genus": str(best["genus"]),
        "species": str(best["species"]),
        "subject": str(best["subject"]),
        "coverage": f"{float(best['coverage']):.2f}",
        "mean_depth": f"{float(best['mean_depth']):.2f}",
        "covered_bases": f"{float(best['covered_bases']):.0f}",
        "num_reads": f"{float(best['num_reads']):.0f}",
        "reference_count": str(len(best["reference_count"])),
    }


def _run_astro_orf2_read_typing(
    db_fasta: Path,
    out_dir: Path,
    single_fastq: str = "",
    fq1: str = "",
    fq2: str = "",
    long_type: str = "",
    threads: int = 4,
    logf=None,
) -> dict[str, str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    bam_path = out_dir / "orf2.screening.bam"
    coverage_path = out_dir / "orf2.coverage.tsv"
    if not bam_path.is_file():
        if fq1:
            map_cmd = " ".join(
                [
                    "minimap2",
                    "-ax",
                    "sr",
                    shlex.quote(str(db_fasta)),
                    shlex.quote(str(fq1)),
                    *([shlex.quote(str(fq2))] if fq2 else []),
                    "|",
                    "samtools",
                    "sort",
                    "-o",
                    shlex.quote(str(bam_path)),
                ]
            )
        else:
            map_cmd = " ".join(
                [
                    "minimap2",
                    "-ax",
                    _choose_minimap2_preset(long_type),
                    shlex.quote(str(db_fasta)),
                    shlex.quote(str(single_fastq)),
                    "-t",
                    str(max(1, int(threads or 1))),
                    "|",
                    "samtools",
                    "sort",
                    "-o",
                    shlex.quote(str(bam_path)),
                ]
            )
        run_command(map_cmd, logf=logf)
    if not bam_path.is_file():
        return {"type": "", "genus": "", "species": "", "subject": "", "identity": "", "coverage": "", "mean_depth": "", "covered_bases": "", "num_reads": "", "reference_count": "0", "method": "read_coverage"}
    run_command(f"samtools index -@ {max(1, int(threads or 1))} {shlex.quote(str(bam_path))}", logf=logf)
    run_command(f"samtools coverage {shlex.quote(str(bam_path))} > {shlex.quote(str(coverage_path))}", logf=logf)
    subject_meta = _build_astro_subject_meta()
    rows = _parse_samtools_coverage_rows(coverage_path)
    enriched_rows = []
    for row in rows:
        subject = str(row.get("reference_name") or "").strip()
        meta = subject_meta.get(subject) or subject_meta.get(subject.split("_", 1)[0]) or {}
        enriched = dict(row)
        enriched["_meta"] = meta
        enriched_rows.append(enriched)
    best = _score_astro_type_rows(enriched_rows)
    return {
        "type": best["type"],
        "genus": best["genus"],
        "species": best["species"],
        "subject": best["subject"],
        "identity": "",
        "coverage": best["coverage"],
        "mean_depth": best["mean_depth"],
        "covered_bases": best["covered_bases"],
        "num_reads": best["num_reads"],
        "reference_count": best["reference_count"],
        "method": "read_coverage",
    }


def _run_astro_orf2_blast_typing(query_fasta: Path, db_fasta: Path, out_dir: Path, logf=None) -> dict[str, str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    out_tsv = out_dir / "orf2.blastn.tsv"
    blastn_bin = _resolve_hadv_blastn_bin()
    cmd = " ".join(
        [
            shlex.quote(blastn_bin),
            "-task blastn",
            "-query",
            shlex.quote(str(query_fasta)),
            "-subject",
            shlex.quote(str(db_fasta)),
            "-outfmt",
            shlex.quote("6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"),
            "-evalue 1e-20",
            "-max_target_seqs 100",
            "-dust no",
            "-num_threads 4",
            "-out",
            shlex.quote(str(out_tsv)),
        ]
    )
    run_command(cmd, logf=logf)
    if not out_tsv.is_file() or out_tsv.stat().st_size == 0:
        return {"type": "", "genus": "", "species": "", "subject": "", "identity": "", "coverage": "", "mean_depth": "", "covered_bases": "", "num_reads": "", "reference_count": "0", "method": "blastn"}
    subject_meta = _build_astro_subject_meta()
    ref_lengths: dict[str, int] = {}
    with db_fasta.open("r", encoding="utf-8", errors="ignore") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            ref_lengths[str(record.id)] = len(str(record.seq))
    by_type: dict[str, dict[str, object]] = {}
    with out_tsv.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            _qseqid, sseqid, pident, length, _mismatch, _gapopen, _qstart, _qend, _sstart, _send, _evalue, bitscore = parts[:12]
            meta = subject_meta.get(sseqid) or subject_meta.get(sseqid.split("_", 1)[0]) or {}
            normalized_type = _normalize_astro_subtype(str(meta.get("normalized_type") or ""))
            genus = str(meta.get("genus") or "").strip()
            species = str(meta.get("species") or "").strip()
            if not normalized_type:
                continue
            ref_len = ref_lengths.get(sseqid, 0)
            try:
                qcov_ref = min(100.0, (float(length) / ref_len) * 100) if ref_len else 0.0
                score = (float(bitscore), float(pident), qcov_ref, float(length))
            except ValueError:
                continue
            bucket = by_type.setdefault(
                normalized_type,
                {
                    "type": normalized_type,
                    "genus": genus,
                    "species": species,
                    "subject": sseqid,
                    "identity": str(pident),
                    "coverage": f"{qcov_ref:.2f}",
                    "covered_bases": str(length),
                    "reference_count": set(),
                    "_score": score,
                },
            )
            if score > bucket["_score"]:
                bucket["subject"] = sseqid
                bucket["identity"] = str(pident)
                bucket["coverage"] = f"{qcov_ref:.2f}"
                bucket["covered_bases"] = str(length)
                bucket["_score"] = score
                bucket["genus"] = genus
                bucket["species"] = species
            bucket["reference_count"].add(str(meta.get("accession_root") or sseqid.split("_", 1)[0]))
    if not by_type:
        return {"type": "", "genus": "", "species": "", "subject": "", "identity": "", "coverage": "", "mean_depth": "", "covered_bases": "", "num_reads": "", "reference_count": "0", "method": "blastn"}
    ranked = sorted(
        by_type.values(),
        key=lambda item: (
            item["_score"][0],
            item["_score"][1],
            item["_score"][2],
            len(item["reference_count"]),
        ),
        reverse=True,
    )
    best = ranked[0]
    return {
        "type": str(best["type"]),
        "genus": str(best["genus"]),
        "species": str(best["species"]),
        "subject": str(best["subject"]),
        "identity": str(best["identity"]),
        "coverage": str(best["coverage"]),
        "mean_depth": "",
        "covered_bases": str(best["covered_bases"]),
        "num_reads": "",
        "reference_count": str(len(best["reference_count"])),
        "method": "blastn",
    }


def _run_astro_orf2_typing(
    pre: str,
    query_fasta: Path | None = None,
    single_fastq: str = "",
    fq1: str = "",
    fq2: str = "",
    long_type: str = "",
    threads: int = 4,
    output_dir: Path | None = None,
    logf=None,
) -> dict[str, object]:
    typing_dir = (output_dir or (Path(f"{pre}_astroviridae_reference_selection") / "typing")).resolve()
    typing_dir.mkdir(parents=True, exist_ok=True)
    db_fasta = _resolve_astroviridae_reference_orf2_fasta()
    if not db_fasta.is_file():
        return {"orf2_type": "", "genus": "", "species": "", "summary_path": "", "hit": {"method": "missing_db"}}
    use_blast = query_fasta is not None and query_fasta.is_file() and query_fasta.stat().st_size > 0
    hit = (
        _run_astro_orf2_blast_typing(query_fasta, db_fasta, typing_dir, logf=logf)
        if use_blast
        else _run_astro_orf2_read_typing(
            db_fasta,
            typing_dir,
            single_fastq=single_fastq,
            fq1=fq1,
            fq2=fq2,
            long_type=long_type,
            threads=threads,
            logf=logf,
        )
    )
    hit_status = "ready"
    hit_note = ""
    if _is_low_support_gene_typing_hit(hit):
        hit_status = "low_support"
        hit_note = _low_support_gene_typing_note("ORF2", hit)
        hit = dict(hit)
        hit["type"] = ""
        hit["genus"] = ""
        hit["species"] = ""
    summary_path = typing_dir / "orf2_typing.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "gene", "matched_type", "genus", "species", "subject", "identity", "coverage", "mean_depth", "covered_bases", "num_reads", "reference_count", "method", "status", "note"])
        writer.writerow([
            pre,
            "orf2",
            hit.get("type", ""),
            hit.get("genus", ""),
            hit.get("species", ""),
            hit.get("subject", ""),
            hit.get("identity", ""),
            hit.get("coverage", ""),
            hit.get("mean_depth", ""),
            hit.get("covered_bases", ""),
            hit.get("num_reads", ""),
            hit.get("reference_count", "0"),
            hit.get("method", ""),
            hit_status,
            hit_note,
        ])
    return {
        "orf2_type": str(hit.get("type") or ""),
        "genus": str(hit.get("genus") or ""),
        "species": str(hit.get("species") or ""),
        "summary_path": str(summary_path.resolve()),
        "hit": hit,
        "status": hit_status,
        "note": hit_note,
    }


def _resolve_astroviridae_vadr_model_dir(project_root: Path, genus: str) -> Path | None:
    normalized = str(genus or "").strip()
    model_name = ""
    if normalized == "Mamastrovirus":
        model_name = "mamastrovirus"
    elif normalized == "Avastrovirus":
        model_name = "avastrovirus"
    if not model_name:
        return None
    model_dir = (project_root / "soft" / "vadr-models-astro" / model_name).resolve()
    if model_dir.is_dir():
        return model_dir
    return None


def _run_vadr_astroviridae_annotation(pre: str, input_fasta: Path, output_root: Path, genus: str, logf=None) -> dict[str, str]:
    project_root = _project_root()
    normalized = str(genus or "").strip()
    model_dir = _resolve_astroviridae_vadr_model_dir(project_root, normalized)
    if model_dir is None:
        return {"status": "missing", "note": f"未找到 {normalized or 'Astroviridae'} VADR 模型目录", "output_dir": "", "gff_path": ""}
    env = _build_vadr_env(project_root, model_dir)
    if env is None:
        return {"status": "missing", "note": "VADR 运行依赖不完整", "output_dir": "", "gff_path": ""}
    astro_env = dict(env)
    ncov_bin_dir = conda_env_path("ncov", "bin")
    astro_env["PATH"] = os.pathsep.join([ncov_bin_dir, "/usr/bin", "/bin", str(astro_env.get("PATH") or "")])
    perl_candidates = [f"{ncov_bin_dir}/perl", "/usr/bin/perl", _resolve_vadr_perl_bin(), "perl"]
    perl_bin = _resolve_working_perl_with_module(astro_env, "Bio::Easel::MSA", perl_candidates)
    vadr_root = output_root / "vadr"
    vadr_root.mkdir(parents=True, exist_ok=True)
    output_dir = vadr_root / f"{pre}_vadr"
    prefix = f"{output_dir.name}.vadr"
    pass_tbl = output_dir / f"{prefix}.pass.tbl"
    fail_tbl = output_dir / f"{prefix}.fail.tbl"
    gff_path = vadr_root / f"{pre}.vadr.gff3"
    if output_dir.is_dir() and (pass_tbl.is_file() or fail_tbl.is_file()) and gff_path.is_file() and gff_path.stat().st_size > 0:
        return {"status": "ready", "note": "已检测到现有星状病毒 VADR 注释结果", "output_dir": str(output_dir.resolve()), "gff_path": str(gff_path.resolve())}
    vadr_script = project_root / "soft" / "vadr" / "v-annotate.pl"
    cmd = " ".join(
        [
            shlex.quote(str(perl_bin)),
            shlex.quote(str(vadr_script)),
            "-f",
            "-r",
            "--ignore_exc",
            "--mkey",
            shlex.quote(str(model_dir.name)),
            "--mdir",
            shlex.quote(str(model_dir)),
            shlex.quote(str(input_fasta)),
            shlex.quote(str(output_dir)),
        ]
    )
    try:
        run_command(cmd, logf=logf, env=astro_env)
    except Exception as exc:
        return {"status": "failed", "note": f"星状病毒 VADR 注释失败: {exc}", "output_dir": str(output_dir.resolve()) if output_dir.exists() else "", "gff_path": ""}
    annotate_tbl2gff = project_root / "soft" / "vadr" / "miniscripts" / "annotate-tbl2gff.pl"
    source_tbl = pass_tbl if _tbl_has_feature_rows(pass_tbl) else fail_tbl
    if _tbl_has_feature_rows(source_tbl):
        try:
            run_command(
                f"{shlex.quote(str(perl_bin))} {shlex.quote(str(annotate_tbl2gff))} {shlex.quote(str(source_tbl))} > {shlex.quote(str(gff_path))}",
                logf=logf,
                env=astro_env,
            )
        except Exception:
            pass
    return {
        "status": "ready" if gff_path.is_file() and gff_path.stat().st_size > 0 else "failed",
        "note": "星状病毒 VADR 注释完成" if gff_path.is_file() and gff_path.stat().st_size > 0 else "星状病毒 VADR 未生成可用 GFF",
        "output_dir": str(output_dir.resolve()),
        "gff_path": str(gff_path.resolve()) if gff_path.is_file() and gff_path.stat().st_size > 0 else "",
    }


def prepare_astroviridae_sample_annotation(pre: str, sample_fasta: Path, output_root: Path, genus: str, logf=None) -> Path | None:
    result = _run_vadr_astroviridae_annotation(pre, sample_fasta, output_root, genus, logf=logf)
    gff_path = str(result.get("gff_path") or "").strip()
    if result.get("status") == "ready" and gff_path:
        gff_file = Path(gff_path)
        if gff_file.is_file() and gff_file.stat().st_size > 0:
            return gff_file
    return None


def prepare_astroviridae_reference_annotation(pre: str, reference_fasta: Path, output_root: Path, genus: str, logf=None) -> Path | None:
    result = _run_vadr_astroviridae_annotation(pre, reference_fasta, output_root, genus, logf=logf)
    gff_path = str(result.get("gff_path") or "").strip()
    if result.get("status") == "ready" and gff_path:
        gff_file = Path(gff_path)
        if gff_file.is_file() and gff_file.stat().st_size > 0:
            return gff_file
    return None


def _parse_gff_attributes(text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in str(text or "").strip().split(";"):
        if not item or "=" not in item:
            continue
        key, value = item.split("=", 1)
        attrs[str(key)] = str(value)
    return attrs


def _normalize_astro_feature_label(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(text or "").strip().lower())


def _is_astro_orf2_annotation(product: str, gene: str, note: str) -> bool:
    raw_values = [str(product or "").strip().lower(), str(gene or "").strip().lower(), str(note or "").strip().lower()]
    normalized = [_normalize_astro_feature_label(value) for value in raw_values]
    product_l = raw_values[0]
    product_n = normalized[0]
    if any("orf2" in value or "orf 2" in value for value in raw_values):
        return True
    if any(value in {"orf2", "orftwo", "ofr2"} for value in normalized):
        return True
    if any("capsid" in value for value in raw_values):
        return True
    if ("structural polyprotein" in product_l or product_n == "structuralpolyprotein") and "non-structural" not in product_l:
        return True
    return False


def _score_astro_orf2_annotation(product: str, gene: str, note: str) -> int:
    product_l = str(product or "").strip().lower()
    gene_l = str(gene or "").strip().lower()
    note_l = str(note or "").strip().lower()
    product_n = _normalize_astro_feature_label(product_l)
    gene_n = _normalize_astro_feature_label(gene_l)
    note_n = _normalize_astro_feature_label(note_l)
    score = 0
    if gene_n in {"orf2", "ofr2"}:
        score += 10
    if "orf2" in product_l or "orf 2" in product_l or product_n == "ofr2":
        score += 8
    if "orf2" in note_l or "orf 2" in note_l:
        score += 4
    if "capsid" in product_l:
        score += 8
    elif "capsid" in note_l:
        score += 4
    if "structural polyprotein" in product_l or product_n == "structuralpolyprotein":
        score += 7
    if "non-structural" in product_l or "nonstructural" in product_n:
        score -= 8
    if gene_n.startswith("orf1") or "orf1" in product_n or "orf1" in note_n:
        score -= 6
    if "polyprotein1ab" in product_n or "rna-dependentrnapolymerase" in note_n:
        score -= 4
    return score


def _find_astro_orf2_feature(gff_path: Path) -> dict[str, object] | None:
    best: dict[str, object] | None = None
    with gff_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            seqid, _source, feature_type, start_text, end_text, _score, strand, _phase, attributes = parts
            attrs = _parse_gff_attributes(attributes)
            product = str(attrs.get("product") or "").strip()
            gene = str(attrs.get("gene") or "").strip()
            note = str(attrs.get("Note") or attrs.get("note") or "").strip()
            if not _is_astro_orf2_annotation(product, gene, note):
                continue
            feature_score = _score_astro_orf2_annotation(product, gene, note)
            if feature_score <= 0:
                continue
            try:
                start = int(start_text)
                end = int(end_text)
            except ValueError:
                continue
            candidate = {
                "seqid": seqid,
                "feature_type": feature_type,
                "start": start,
                "end": end,
                "strand": strand,
                "product": product,
                "gene": gene,
                "note": note,
                "score": feature_score,
            }
            if best is None:
                best = candidate
                continue
            best_score = int(best.get("score") or 0)
            if feature_score > best_score:
                best = candidate
                continue
            if feature_type == "CDS" and str(best.get("feature_type") or "") != "CDS":
                best = candidate
                continue
            candidate_len = end - start + 1
            best_len = int(best.get("end") or 0) - int(best.get("start") or 0) + 1
            if candidate_len > best_len:
                best = candidate
    return best


def _extract_astro_orf2_region(sample_fasta: Path, gff_path: Path, sample_label: str) -> tuple[SeqRecord | None, dict[str, object]]:
    feature = _find_astro_orf2_feature(gff_path)
    if feature is None:
        return None, {"status": "missing", "note": "GFF 中未找到 ORF2/capsid 区域"}
    records = list(SeqIO.parse(str(sample_fasta), "fasta"))
    if not records:
        return None, {"status": "missing", "note": f"{sample_fasta.name} 中未找到序列记录"}
    seqid = str(feature.get("seqid") or "").strip()
    record = next((item for item in records if str(item.id).strip() == seqid), None)
    if record is None:
        record = max(records, key=lambda item: len(str(item.seq)))
    start = int(feature["start"])
    end = int(feature["end"])
    sequence = str(record.seq)
    if start < 1 or end > len(sequence) or start > end:
        return None, {"status": "missing", "note": "ORF2 坐标超出样本序列长度"}
    fragment = sequence[start - 1:end]
    if str(feature.get("strand") or "+") == "-":
        fragment = str(Seq(fragment).reverse_complement())
    fragment = fragment.strip().strip("Nn")
    if not fragment:
        return None, {"status": "missing", "note": "ORF2 区域仅包含空白或 N"}
    record_id = f"{sample_label}_orf2_sample"
    return SeqRecord(Seq(fragment), id=record_id, description=f"{sample_label} ORF2 sample"), {
        "status": "ready",
        "start": start,
        "end": end,
        "strand": str(feature.get("strand") or "+"),
        "seqid": str(record.id),
        "length": len(fragment),
    }


def _build_astro_tree_member_label(subtype: str, member_label: str, suffix: str = "") -> str:
    parts = [item for item in [str(subtype or "").strip(), str(member_label or "").strip(), str(suffix or "").strip()] if item]
    return _sanitize_tree_label("_".join(parts))


def _resolve_astro_orf2_fallback_gff(orf2_type: str, genus: str) -> Path | None:
    normalized_type = _normalize_astro_subtype(orf2_type)
    normalized_genus = str(genus or "").strip()
    genus_match: Path | None = None
    for row in _load_astroviridae_reference_orf2_manifest():
        row_genus = str(row.get("genus") or "").strip()
        if normalized_genus and row_genus != normalized_genus:
            continue
        gff_path = Path(str(row.get("gff_path") or "").strip())
        if not gff_path.is_file() or gff_path.stat().st_size == 0:
            continue
        row_type = _normalize_astro_subtype(str(row.get("abbrev") or ""))
        if normalized_type and row_type == normalized_type:
            return gff_path
        if genus_match is None:
            genus_match = gff_path
    return genus_match


def _filter_astro_orf2_refs_by_group(orf2_type: str, genus: str) -> tuple[list[dict[str, str]], list[SeqRecord]]:
    target_type = _normalize_astro_subtype(orf2_type)
    target_genus = str(genus or "").strip()
    exact_rows: list[dict[str, str]] = []
    genus_rows: list[dict[str, str]] = []
    for row in _load_astroviridae_reference_orf2_manifest():
        row_genus = str(row.get("genus") or "").strip()
        if target_genus and row_genus != target_genus:
            continue
        orf2_fasta_path = Path(str(row.get("orf2_fasta_path") or "").strip())
        if not orf2_fasta_path.is_file():
            continue
        row_type = _normalize_astro_subtype(str(row.get("abbrev") or ""))
        if target_type and row_type == target_type:
            exact_rows.append(row)
        genus_rows.append(row)
    selected_rows = exact_rows or genus_rows
    selected_records: list[SeqRecord] = []
    ordered_rows: list[dict[str, str]] = []
    for row in selected_rows:
        orf2_fasta_path = Path(str(row.get("orf2_fasta_path") or "").strip())
        record = next(SeqIO.parse(str(orf2_fasta_path), "fasta"), None)
        if record is None:
            continue
        subtype = _normalize_astro_subtype(str(row.get("abbrev") or ""))
        accession = str(row.get("accession") or row.get("accession_full") or "").split(".", 1)[0]
        tree_label = _build_astro_tree_member_label(subtype or target_genus, accession or record.id)
        ordered_rows.append(
            {
                "gene": "ORF2",
                "subtype": subtype,
                "accession": accession,
                "backup_rank": "1",
                "tree_label": tree_label,
            }
        )
        selected_records.append(
            SeqRecord(
                record.seq,
                id=tree_label,
                name=tree_label,
                description=f"{subtype} {accession}".strip(),
            )
        )
    return ordered_rows, selected_records


def build_astroviridae_orf2_phylogeny_assets(pre: str, sample_fasta: Path, gff_path: Path, genus: str = "", logf=None) -> dict[str, object]:
    output_root = Path(f"{pre}_astroviridae_reference_selection") / "phylogeny"
    output_root.mkdir(parents=True, exist_ok=True)
    selection_path = Path(f"{pre}_astroviridae_reference_selection") / "selection.tsv"
    orf2_type = ""
    selected_genus = str(genus or "").strip()
    if selection_path.is_file():
        with selection_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
            row = next(csv.DictReader(handle, delimiter="\t"), None) or {}
        orf2_type = str(row.get("orf2_type") or "").strip()
        selected_genus = selected_genus or str(row.get("genus") or "").strip()
    result_row = {
        "gene": "ORF2",
        "subtype": orf2_type,
        "genogroup": selected_genus,
        "status": "missing",
        "note": "",
        "sample_region_path": "",
        "backup_fasta": "",
        "tree_path": "",
        "member_count": "0",
    }
    sample_record, sample_meta = _extract_astro_orf2_region(sample_fasta, gff_path, _sanitize_tree_label(pre))
    if sample_record is None:
        fallback_gff = _resolve_astro_orf2_fallback_gff(orf2_type, selected_genus)
        if fallback_gff is not None and fallback_gff != gff_path:
            sample_record, sample_meta = _extract_astro_orf2_region(sample_fasta, fallback_gff, _sanitize_tree_label(pre))
            if sample_record is not None:
                sample_meta = dict(sample_meta)
                sample_meta["note"] = f"样本 VADR 注释未定位 ORF2，已回退到 {fallback_gff.name} 坐标"
    if sample_record is None:
        result_row["note"] = str(sample_meta.get("note") or "未提取到 ORF2 区域")
    else:
        sample_label = _build_astro_tree_member_label(orf2_type or selected_genus or "Astroviridae", pre, "sample")
        sample_record.id = sample_label
        sample_record.name = sample_label
        sample_record.description = f"{orf2_type or selected_genus or 'Astroviridae'} {pre} sample".strip()
        backup_rows, backup_records = _filter_astro_orf2_refs_by_group(orf2_type, selected_genus)
        if not backup_records:
            result_row["note"] = f"ORF2 参考库中未找到 {selected_genus or 'Astroviridae'} 记录"
        else:
            tree_result = _build_norovirus_gene_tree("orf2", sample_record, backup_rows, backup_records, output_root / "orf2", logf=logf)
            result_row["status"] = str(tree_result.get("status") or "missing")
            fallback_note = str(sample_meta.get("note") or "").strip()
            tree_note = str(tree_result.get("note") or "").strip()
            result_row["note"] = "；".join([item for item in [fallback_note, tree_note] if item])
            result_row["sample_region_path"] = str((output_root / "orf2" / "orf2.sample.fasta").resolve())
            result_row["backup_fasta"] = str(tree_result.get("backup_fasta") or "")
            result_row["tree_path"] = str(tree_result.get("tree_path") or "")
            result_row["member_count"] = str(tree_result.get("member_count") or "0")
    summary_path = output_root / "summary.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["gene", "subtype", "genogroup", "status", "note", "sample_region_path", "backup_fasta", "tree_path", "member_count"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(result_row)
    return {
        "status": "ready" if result_row.get("status") == "ready" else "missing",
        "summary_path": str(summary_path.resolve()),
        "rows": [result_row],
    }


def resolve_astroviridae_reference(
    pre: str,
    species: str = "",
    requested_ref: str = "",
    single_fastq: str = "",
    fq1: str = "",
    fq2: str = "",
    long_type: str = "",
    threads: int = 4,
    query_fasta: str = "",
    logf=None,
) -> dict[str, str]:
    _clear_sample_skip_flag(pre)
    requested = str(requested_ref or "").strip()
    if requested and Path(requested).is_file():
        inferred_genus = _extract_astro_genus(species)
        return {
            "status": "ready",
            "orf2_type": _normalize_astro_subtype(species),
            "genus": inferred_genus,
            "species_label": _astro_species_label(inferred_genus, species),
            "reference_path": str(Path(requested).resolve()),
            "gff_path": "nogtf",
            "summary_path": "",
            "typing_summary_path": "",
        }
    screening_dir = Path(f"{pre}_astroviridae_reference_selection")
    screening_dir.mkdir(exist_ok=True)
    query_path = Path(query_fasta) if str(query_fasta or "").strip() else None
    typing_result = _run_astro_orf2_typing(
        pre,
        query_fasta=query_path,
        single_fastq=single_fastq,
        fq1=fq1,
        fq2=fq2,
        long_type=long_type,
        threads=threads,
        output_dir=screening_dir / "typing",
        logf=logf,
    )
    orf2_type = str(typing_result.get("orf2_type") or "").strip()
    genus = str(typing_result.get("genus") or "").strip()
    species_label = str(typing_result.get("species") or "").strip()
    summary_path = screening_dir / "selection.tsv"
    if not orf2_type or genus not in {"Mamastrovirus", "Avastrovirus"}:
        status_text = "skipped" if str(typing_result.get("status") or "").strip() == "low_support" else "missing"
        note_text = str(typing_result.get("note") or "").strip() or f"未能根据 ORF2 或 species={species or '-'} 确定星状病毒分型"
        with summary_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sample", "orf2_type", "genus", "species", "status", "note", "typing_summary_path"])
            writer.writerow([pre, orf2_type, genus, species_label, status_text, note_text, str(typing_result.get("summary_path") or "")])
        if status_text == "skipped":
            _write_sample_skip_flag(pre, note_text)
        return {
            "status": status_text,
            "orf2_type": orf2_type,
            "genus": genus,
            "species_label": species_label or species or "Astroviridae",
            "reference_path": "",
            "gff_path": "nogtf",
            "summary_path": str(summary_path.resolve()),
            "typing_summary_path": str(typing_result.get("summary_path") or ""),
        }
    manifest_rows = _load_astroviridae_candidate_reference_manifest()
    exact_candidates = [
        row for row in manifest_rows
        if str(row.get("genus") or "").strip() == genus
        and _normalize_astro_subtype(str(row.get("abbrev") or "")) == orf2_type
    ]
    species_candidates = [
        row for row in manifest_rows
        if str(row.get("genus") or "").strip() == genus
        and str(row.get("species") or "").strip() == species_label
    ]
    genus_candidates = [row for row in manifest_rows if str(row.get("genus") or "").strip() == genus]
    candidates = exact_candidates or species_candidates or genus_candidates
    candidate_records = _build_enterovirus_candidate_records(candidates)
    dedup_candidate_records = _deduplicate_enterovirus_candidate_records(
        candidate_records,
        screening_dir=screening_dir,
        threads=threads,
        logf=logf,
    )
    if not dedup_candidate_records:
        with summary_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sample", "orf2_type", "genus", "species", "status", "note", "typing_summary_path"])
            writer.writerow([pre, orf2_type, genus, species_label, "missing", "未找到对应亚型的全基因组候选参考", str(typing_result.get("summary_path") or "")])
        return {
            "status": "missing",
            "orf2_type": orf2_type,
            "genus": genus,
            "species_label": _astro_species_label(genus, species_label),
            "reference_path": "",
            "gff_path": "nogtf",
            "summary_path": str(summary_path.resolve()),
            "typing_summary_path": str(typing_result.get("summary_path") or ""),
        }
    final_rows, audit_rows, anchor_record = _select_enterovirus_candidates_in_batches(
        dedup_candidate_records,
        screening_dir=screening_dir,
        single_fastq=single_fastq,
        fq1=fq1,
        fq2=fq2,
        long_type=long_type,
        threads=threads,
        logf=logf,
    )
    candidate_fasta = screening_dir / "candidate_references.fasta"
    candidate_meta = _write_reference_records(dedup_candidate_records, candidate_fasta)
    coverage_path = screening_dir / "candidate_references.coverage.tsv"
    coverage_rows = final_rows or audit_rows
    _write_aggregated_coverage_rows(coverage_path, coverage_rows)
    dedup_manifest_path = screening_dir / "candidate_references.dedup.tsv"
    with dedup_manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["representative_id", "accession", "member_count", "member_accessions"])
        for record in dedup_candidate_records:
            writer.writerow([record.get("fasta_id", ""), record.get("accession", ""), record.get("member_count", 1), record.get("member_accessions", "")])
    best = None
    for row in coverage_rows:
        ref_name = str(row.get("reference_name") or "").strip()
        meta = candidate_meta.get(ref_name)
        if meta is None:
            continue
        score = _coverage_row_score(row)
        if best is None or score > best["score"]:
            best = {"score": score, "meta": meta, "coverage_row": row}
    if best is None:
        return {
            "status": "missing",
            "orf2_type": orf2_type,
            "genus": genus,
            "species_label": _astro_species_label(genus, species_label),
            "reference_path": "",
            "gff_path": "nogtf",
            "summary_path": "",
            "typing_summary_path": str(typing_result.get("summary_path") or ""),
        }
    best_meta = best["meta"]
    best_row = best["coverage_row"]
    best_accession = str(best_meta.get("accession") or "").strip()
    anchor_accession = str(anchor_record.get("accession") or "").strip() if anchor_record else ""
    best_reference_fasta = screening_dir / f"{_sanitize_tree_label(orf2_type or best_accession)}.reference.fasta"
    with best_reference_fasta.open("w", encoding="utf-8") as handle:
        handle.write(f">{best_accession or 'reference'}\n")
        sequence = str(best_meta.get("sequence") or "")
        for index in range(0, len(sequence), 80):
            handle.write(sequence[index:index + 80] + "\n")
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "orf2_type", "genus", "species", "coverage", "mean_depth", "covered_bases", "num_reads", "reference_name", "accession", "reference_path", "gff_path", "typing_summary_path", "anchor_accession", "candidate_count", "dedup_candidate_count"])
        writer.writerow([
            pre,
            orf2_type,
            genus,
            species_label,
            f"{float(best_row.get('coverage') or 0.0):.6f}",
            f"{float(best_row.get('mean_depth') or 0.0):.6f}",
            f"{float(best_row.get('covered_bases') or 0.0):.0f}",
            f"{float(best_row.get('num_reads') or 0.0):.0f}",
            str(best_meta.get("meta", {}).get("header") or best_accession),
            best_accession,
            str(best_reference_fasta.resolve()),
            str(best_meta.get("meta", {}).get("gff_path") or "nogtf"),
            str(typing_result.get("summary_path") or ""),
            anchor_accession,
            len(candidate_records),
            len(dedup_candidate_records),
        ])
    return {
        "status": "ready",
        "orf2_type": orf2_type,
        "genus": genus,
        "species_label": _astro_species_label(genus, species_label),
        "reference_path": str(best_reference_fasta.resolve()),
        "gff_path": str(best_meta.get("meta", {}).get("gff_path") or "nogtf"),
        "summary_path": str(summary_path.resolve()),
        "typing_summary_path": str(typing_result.get("summary_path") or ""),
    }


def run_astroviridae_consensus_typing(pre: str, consensus_fasta: Path, logf=None) -> dict[str, str]:
    if not consensus_fasta.is_file() or consensus_fasta.stat().st_size == 0:
        return {"status": "missing", "summary_path": "", "orf2_type": "", "genus": "", "species": ""}
    output_dir = Path(f"{pre}_astroviridae_reference_selection") / "consensus_typing"
    result = _run_astro_orf2_typing(pre, query_fasta=consensus_fasta, output_dir=output_dir, logf=logf)
    summary_path = output_dir / "consensus_typing.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "orf2_type", "genus", "species", "typing_summary_path"])
        writer.writerow([pre, result.get("orf2_type", ""), result.get("genus", ""), result.get("species", ""), result.get("summary_path", "")])
    return {
        "status": "ready",
        "summary_path": str(summary_path.resolve()),
        "orf2_type": str(result.get("orf2_type") or ""),
        "genus": str(result.get("genus") or ""),
        "species": str(result.get("species") or ""),
    }


def _resolve_rotavirus_db_dir() -> Path:
    env_root = str(os.environ.get("META_ROTAVIRUS_DB_DIR") or "").strip()
    if env_root:
        return Path(env_root).expanduser().resolve()
    database_root = str(os.environ.get("META_DATABASE_ROOT") or "").strip()
    project_root = _project_root()
    candidates: list[Path] = []
    if database_root:
        candidates.append((Path(database_root) / "virus" / "Rotavirus").expanduser())
    candidates.extend(
        [
            project_root / "database" / "virus" / "Rotavirus",
            Path("/data/deploy/meta_genome/database/virus/Rotavirus"),
        ]
    )
    return next((candidate.resolve() for candidate in candidates if candidate.exists()), candidates[0].resolve())


def _resolve_rotavirus_full_manifest() -> Path:
    return (_resolve_rotavirus_db_dir() / "full_genomes" / "human_rotavirus_A_B_C_manifest.tsv").resolve()


def _resolve_rotavirus_full_fasta() -> Path:
    return (_resolve_rotavirus_db_dir() / "full_genomes" / "human_rotavirus_A_B_C_complete_genomes.fasta").resolve()


def _resolve_rotavirus_a_subtype_manifest() -> Path:
    return (_resolve_rotavirus_db_dir() / "subtype_complete_genomes" / "rotavirus_a_subtype_complete_genomes_manifest.tsv").resolve()


def _resolve_rotavirus_a_subtype_summary() -> Path:
    return (_resolve_rotavirus_db_dir() / "subtype_complete_genomes" / "rotavirus_a_subtype_complete_genomes_summary.tsv").resolve()


def _resolve_rotavirus_a_subtype_fasta() -> Path:
    return (_resolve_rotavirus_db_dir() / "subtype_complete_genomes" / "rotavirus_a_subtype_complete_genomes.fasta").resolve()


def _resolve_rotavirus_group_bc_manifest() -> Path:
    return (_resolve_rotavirus_db_dir() / "group_bc_complete_genomes" / "rotavirus_bc_manifest.tsv").resolve()


def _resolve_rotavirus_group_bc_fasta() -> Path:
    return (_resolve_rotavirus_db_dir() / "group_bc_complete_genomes" / "rotavirus_b_bang373_complete_segments.fasta").resolve()


def _rotavirus_segment_sort_key(segment_label: str) -> tuple[int, str]:
    segment = str(segment_label or "").strip().upper()
    order = {
        "VP1": 1,
        "VP2": 2,
        "VP3": 3,
        "VP4": 4,
        "VP6": 5,
        "VP7": 6,
        "NSP1": 7,
        "NSP2": 8,
        "NSP3": 9,
        "NSP4": 10,
        "NSP5": 11,
    }
    return (order.get(segment, 999), segment)


def _load_rotavirus_manifest(path: Path) -> list[dict[str, str]]:
    return _load_rhinovirus_manifest(path)


def _choose_rotavirus_complete_isolate(rows: list[dict[str, str]], target_group: str = "") -> list[dict[str, str]]:
    filtered_rows = [
        row for row in rows
        if (not target_group or str(row.get("species_group") or "").strip().upper() == target_group)
        and str(row.get("isolate_label") or "").strip()
        and str(row.get("segment_label") or "").strip()
    ]
    by_isolate: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)
    for row in filtered_rows:
        isolate_label = str(row.get("isolate_label") or "").strip()
        segment_label = str(row.get("segment_label") or "").strip().upper()
        current = by_isolate[isolate_label].get(segment_label)
        if current is None:
            by_isolate[isolate_label][segment_label] = row
            continue
        current_nc = str(current.get("accession") or "").strip().startswith("NC_")
        next_nc = str(row.get("accession") or "").strip().startswith("NC_")
        current_len = int(current.get("sequence_length") or 0)
        next_len = int(row.get("sequence_length") or 0)
        if (next_nc and not current_nc) or (next_nc == current_nc and next_len > current_len):
            by_isolate[isolate_label][segment_label] = row
    ranked: list[tuple[int, int, str, list[dict[str, str]]]] = []
    for isolate_label, segments in by_isolate.items():
        selected = list(segments.values())
        nc_count = sum(1 for row in selected if str(row.get("accession") or "").strip().startswith("NC_"))
        ranked.append((len(selected), nc_count, isolate_label, sorted(selected, key=lambda row: _rotavirus_segment_sort_key(str(row.get("segment_label") or "")))))
    if not ranked:
        return []
    ranked.sort(key=lambda item: (item[0], item[1], item[2]), reverse=True)
    return ranked[0][3]


def _load_rotavirus_group_panel_rows() -> list[dict[str, str]]:
    panel_rows: list[dict[str, str]] = []
    subtype_rows = _load_rotavirus_manifest(_resolve_rotavirus_a_subtype_manifest())
    for row in subtype_rows:
        segment = str(row.get("segment_label") or "").strip().upper()
        if segment not in {"VP1", "VP2", "VP3", "VP4", "VP6", "VP7", "NSP1", "NSP2", "NSP3", "NSP4", "NSP5"}:
            continue
        enriched = dict(row)
        enriched["species_group"] = "A"
        enriched["group_label"] = "A"
        enriched["subtype_combo"] = str(row.get("combo") or "").strip()
        enriched["panel_source"] = "rotavirus_a_subtype_complete"
        enriched["fasta_source"] = str(_resolve_rotavirus_a_subtype_fasta())
        panel_rows.append(enriched)
    full_rows = _load_rotavirus_manifest(_resolve_rotavirus_full_manifest())
    for group in ["B", "C"]:
        selected = _choose_rotavirus_complete_isolate(full_rows, target_group=group)
        for row in selected:
            enriched = dict(row)
            enriched["species_group"] = group
            enriched["group_label"] = group
            enriched["subtype_combo"] = ""
            enriched["panel_source"] = f"rotavirus_{group.lower()}_complete"
            enriched["fasta_source"] = str(_resolve_rotavirus_full_fasta())
            panel_rows.append(enriched)
    return panel_rows


def _build_rotavirus_group_subject_meta() -> dict[str, dict[str, str]]:
    subject_meta: dict[str, dict[str, str]] = {}
    for row in _load_rotavirus_group_panel_rows():
        accession = str(row.get("accession") or "").strip()
        accession_root = accession.split(".", 1)[0]
        header = str(row.get("header") or "").strip()
        meta = {
            "species_group": str(row.get("species_group") or "").strip().upper(),
            "segment_label": str(row.get("segment_label") or "").strip().upper(),
            "subtype_combo": str(row.get("subtype_combo") or "").strip(),
            "isolate_label": str(row.get("isolate_label") or "").strip(),
            "accession": accession,
            "accession_root": accession_root,
            "header": header,
            "fasta_source": str(row.get("fasta_source") or "").strip(),
            "sequence_length": str(row.get("sequence_length") or "").strip(),
        }
        for key in {accession, accession_root, header, _sanitize_tree_label(accession_root), _sanitize_tree_label(accession)}:
            if key:
                subject_meta[key] = meta
    return subject_meta


def _write_rotavirus_group_panel_fasta(output_path: Path) -> list[dict[str, str]]:
    panel_rows = _load_rotavirus_group_panel_rows()
    by_fasta_source: dict[str, dict[str, SeqRecord]] = {}
    with output_path.open("w", encoding="utf-8") as handle:
        for row in panel_rows:
            fasta_source = str(row.get("fasta_source") or "").strip()
            sequence_map = by_fasta_source.setdefault(fasta_source, _build_rhinovirus_sequence_map(Path(fasta_source)))
            header = str(row.get("header") or "").strip()
            accession = str(row.get("accession") or "").strip()
            record = sequence_map.get(header) or sequence_map.get(accession)
            if record is None:
                continue
            fasta_id = _sanitize_tree_label(accession.split(".", 1)[0] or accession)
            handle.write(f">{fasta_id}\n")
            sequence = str(record.seq)
            for index in range(0, len(sequence), 80):
                handle.write(sequence[index:index + 80] + "\n")
    return panel_rows


def _rank_rotavirus_group_support(entries: list[dict[str, object]]) -> list[dict[str, str]]:
    by_group_isolate: dict[tuple[str, str], dict[str, object]] = {}
    for entry in entries:
        group = str(entry.get("species_group") or "").strip().upper()
        isolate = str(entry.get("isolate_label") or "").strip()
        segment = str(entry.get("segment_label") or "").strip().upper()
        if group not in {"A", "B", "C"} or not isolate or not segment:
            continue
        coverage = float(entry.get("coverage") or 0.0)
        depth_value = float(entry.get("depth_value") or 0.0)
        covered_bases = float(entry.get("covered_bases") or 0.0)
        support_value = float(entry.get("support_value") or 0.0)
        ref_length = float(entry.get("reference_bases") or 0.0)
        reference_name = str(entry.get("reference_name") or "").strip()
        score = (coverage, depth_value, covered_bases, support_value)
        bucket = by_group_isolate.setdefault(
            (group, isolate),
            {
                "species_group": group,
                "isolate_label": isolate,
                "segment_hits": set(),
                "coverage_sum": 0.0,
                "depth_sum": 0.0,
                "covered_bases_sum": 0.0,
                "num_reads_sum": 0.0,
                "reference_bases_sum": 0.0,
                "best_reference_name": reference_name,
                "best_row_score": score,
            },
        )
        if segment not in bucket["segment_hits"]:
            bucket["reference_bases_sum"] = float(bucket["reference_bases_sum"]) + ref_length
        bucket["segment_hits"].add(segment)
        bucket["coverage_sum"] = float(bucket["coverage_sum"]) + coverage
        bucket["depth_sum"] = float(bucket["depth_sum"]) + depth_value
        bucket["covered_bases_sum"] = float(bucket["covered_bases_sum"]) + (min(covered_bases, ref_length) if ref_length > 0 else covered_bases)
        bucket["num_reads_sum"] = float(bucket["num_reads_sum"]) + support_value
        if score > bucket["best_row_score"]:
            bucket["best_row_score"] = score
            bucket["best_reference_name"] = reference_name
    ranked = sorted(
        by_group_isolate.values(),
        key=lambda item: (
            len(item["segment_hits"]),
            ((float(item["covered_bases_sum"]) / float(item["reference_bases_sum"]) * 100.0) if float(item["reference_bases_sum"]) > 0 else 0.0),
            float(item["covered_bases_sum"]),
            float(item["depth_sum"]),
            float(item["num_reads_sum"]),
        ),
        reverse=True,
    )
    rows: list[dict[str, str]] = []
    for index, item in enumerate(ranked):
        ref_bases = float(item.get("reference_bases_sum") or 0.0)
        covered_bases = float(item.get("covered_bases_sum") or 0.0)
        rows.append(
            {
                "species_group": str(item["species_group"]),
                "isolate_label": str(item["isolate_label"]),
                "reference_name": str(item["best_reference_name"]),
                "segment_count": str(len(item["segment_hits"])),
                "coverage_sum": f"{float(item['coverage_sum']):.2f}",
                "full_length_coverage_pct": f"{((covered_bases / ref_bases) * 100.0) if ref_bases > 0 else 0.0:.2f}",
                "depth_sum": f"{float(item['depth_sum']):.2f}",
                "covered_bases_sum": f"{covered_bases:.0f}",
                "reference_bases_sum": f"{ref_bases:.0f}",
                "num_reads_sum": f"{float(item['num_reads_sum']):.0f}",
                "selected_group_type": str(item["species_group"]) if index == 0 else "",
                "is_selected": "yes" if index == 0 else "no",
            }
        )
    return rows


def _score_rotavirus_group_rows(rows: list[dict[str, object]]) -> dict[str, str]:
    entries: list[dict[str, object]] = []
    for row in rows:
        meta = row.get("_meta") if isinstance(row.get("_meta"), dict) else {}
        group = str(meta.get("species_group") or "").strip().upper()
        isolate = str(meta.get("isolate_label") or "").strip()
        segment = str(meta.get("segment_label") or "").strip().upper()
        if group not in {"A", "B", "C"} or not isolate or not segment:
            continue
        try:
            ref_length = max(0.0, float(meta.get("sequence_length") or 0.0))
        except ValueError:
            ref_length = 0.0
        entries.append(
            {
                "species_group": group,
                "isolate_label": isolate,
                "segment_label": segment,
                "coverage": float(row.get("coverage") or 0.0),
                "depth_value": float(row.get("mean_depth") or 0.0),
                "covered_bases": float(row.get("covered_bases") or 0.0),
                "support_value": float(row.get("num_reads") or 0.0),
                "reference_bases": ref_length,
                "reference_name": str(row.get("reference_name") or ""),
            }
        )
    ranked = _rank_rotavirus_group_support(entries)
    if not ranked:
        return {"species_group": "", "isolate_label": "", "reference_name": "", "segment_count": "0", "coverage_sum": "", "full_length_coverage_pct": "", "depth_sum": "", "covered_bases_sum": "", "reference_bases_sum": "", "num_reads_sum": "", "selected_group_type": "", "is_selected": "no"}
    return ranked[0]


def _run_rotavirus_group_read_typing(
    db_fasta: Path,
    out_dir: Path,
    single_fastq: str = "",
    fq1: str = "",
    fq2: str = "",
    long_type: str = "",
    threads: int = 4,
    logf=None,
) -> dict[str, object]:
    out_dir.mkdir(parents=True, exist_ok=True)
    bam_path = out_dir / "group.screening.bam"
    coverage_path = out_dir / "group.coverage.tsv"
    if fq1:
        map_cmd = " ".join(
            [
                "minimap2",
                "-ax",
                "sr",
                shlex.quote(str(db_fasta)),
                shlex.quote(str(fq1)),
                *([shlex.quote(str(fq2))] if fq2 else []),
                "|",
                "samtools",
                "sort",
                "-o",
                shlex.quote(str(bam_path)),
            ]
        )
    else:
        map_cmd = " ".join(
            [
                "minimap2",
                "-ax",
                _choose_minimap2_preset(long_type),
                shlex.quote(str(db_fasta)),
                shlex.quote(str(single_fastq)),
                "-t",
                str(max(1, int(threads or 1))),
                "|",
                "samtools",
                "sort",
                "-o",
                shlex.quote(str(bam_path)),
            ]
        )
    run_command(map_cmd, logf=logf)
    run_command(f"samtools index -@ {max(1, int(threads or 1))} {shlex.quote(str(bam_path))}", logf=logf)
    run_command(f"samtools coverage {shlex.quote(str(bam_path))} > {shlex.quote(str(coverage_path))}", logf=logf)
    subject_meta = _build_rotavirus_group_subject_meta()
    coverage_rows = _parse_samtools_coverage_rows(coverage_path)
    enriched_rows = []
    for row in coverage_rows:
        subject = str(row.get("reference_name") or "").strip()
        meta = subject_meta.get(subject) or subject_meta.get(subject.split("_", 1)[0]) or {}
        enriched = dict(row)
        enriched["_meta"] = meta
        enriched_rows.append(enriched)
    ranked_rows = _rank_rotavirus_group_support(
        [
            {
                "species_group": str((row.get("_meta") or {}).get("species_group") or "").strip().upper(),
                "isolate_label": str((row.get("_meta") or {}).get("isolate_label") or "").strip(),
                "segment_label": str((row.get("_meta") or {}).get("segment_label") or "").strip().upper(),
                "coverage": float(row.get("coverage") or 0.0),
                "depth_value": float(row.get("mean_depth") or 0.0),
                "covered_bases": float(row.get("covered_bases") or 0.0),
                "support_value": float(row.get("num_reads") or 0.0),
                "reference_bases": float((row.get("_meta") or {}).get("sequence_length") or 0.0),
                "reference_name": str(row.get("reference_name") or ""),
            }
            for row in enriched_rows
            if isinstance(row.get("_meta"), dict)
        ]
    )
    best = ranked_rows[0] if ranked_rows else {"species_group": "", "isolate_label": "", "reference_name": "", "segment_count": "0", "coverage_sum": "", "full_length_coverage_pct": "", "depth_sum": "", "covered_bases_sum": "", "reference_bases_sum": "", "num_reads_sum": "", "selected_group_type": "", "is_selected": "no"}
    best["method"] = "read_coverage"
    return {"hit": best, "rows": ranked_rows}


def _run_rotavirus_group_blast_typing(query_fasta: Path, db_fasta: Path, out_dir: Path, logf=None) -> dict[str, str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    out_tsv = out_dir / "group.blastn.tsv"
    blastn_bin = _resolve_hadv_blastn_bin()
    cmd = " ".join(
        [
            shlex.quote(blastn_bin),
            "-task blastn",
            "-query",
            shlex.quote(str(query_fasta)),
            "-subject",
            shlex.quote(str(db_fasta)),
            "-outfmt",
            shlex.quote("6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"),
            "-evalue 1e-20",
            "-max_target_seqs 100",
            "-dust no",
            "-num_threads 4",
            "-out",
            shlex.quote(str(out_tsv)),
        ]
    )
    run_command(cmd, logf=logf)
    if not out_tsv.is_file() or out_tsv.stat().st_size == 0:
        return {"hit": {"species_group": "", "isolate_label": "", "reference_name": "", "segment_count": "0", "coverage_sum": "", "full_length_coverage_pct": "", "depth_sum": "", "covered_bases_sum": "", "reference_bases_sum": "", "num_reads_sum": "", "selected_group_type": "", "is_selected": "no", "method": "blastn"}, "rows": []}
    subject_meta = _build_rotavirus_group_subject_meta()
    ref_lengths = {str(record.id): len(str(record.seq)) for record in SeqIO.parse(str(db_fasta), "fasta")}
    by_group_isolate: dict[tuple[str, str, str], dict[str, object]] = {}
    with out_tsv.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            _qseqid, sseqid, pident, length, _mismatch, _gapopen, _qstart, _qend, _sstart, _send, _evalue, bitscore = parts[:12]
            meta = subject_meta.get(sseqid) or subject_meta.get(sseqid.split("_", 1)[0]) or {}
            group = str(meta.get("species_group") or "").strip().upper()
            isolate = str(meta.get("isolate_label") or "").strip()
            segment = str(meta.get("segment_label") or "").strip().upper()
            if group not in {"A", "B", "C"} or not isolate or not segment:
                continue
            ref_len = ref_lengths.get(sseqid, 0)
            try:
                qcov_ref = min(100.0, (float(length) / ref_len) * 100) if ref_len else 0.0
                score = (float(bitscore), float(pident), qcov_ref, float(length))
            except ValueError:
                continue
            bucket = by_group_isolate.setdefault(
                (group, isolate, segment),
                {
                    "group": group,
                    "isolate": isolate,
                    "segment": segment,
                    "coverage": qcov_ref,
                    "bitscore": float(bitscore),
                    "length": float(length),
                    "reference_name": sseqid,
                },
            )
            if score > (float(bucket["bitscore"]), 0.0, float(bucket["coverage"]), float(bucket["length"])):
                bucket["coverage"] = qcov_ref
                bucket["bitscore"] = float(bitscore)
                bucket["length"] = float(length)
                bucket["reference_name"] = sseqid
    aggregate: dict[tuple[str, str], dict[str, object]] = {}
    for item in by_group_isolate.values():
        key = (str(item["group"]), str(item["isolate"]))
        bucket = aggregate.setdefault(
            key,
            {
                "species_group": str(item["group"]),
                "isolate_label": str(item["isolate"]),
                "segment_hits": set(),
                "coverage_sum": 0.0,
                "bitscore_sum": 0.0,
                "covered_bases_sum": 0.0,
                "reference_bases_sum": 0.0,
                "reference_name": str(item["reference_name"]),
            },
        )
        bucket["segment_hits"].add(str(item["segment"]))
        bucket["coverage_sum"] = float(bucket["coverage_sum"]) + float(item["coverage"])
        bucket["bitscore_sum"] = float(bucket["bitscore_sum"]) + float(item["bitscore"])
        ref_len = ref_lengths.get(str(item["reference_name"]), 0)
        bucket["reference_bases_sum"] = float(bucket["reference_bases_sum"]) + float(ref_len)
        bucket["covered_bases_sum"] = float(bucket["covered_bases_sum"]) + min(float(item["length"]), float(ref_len)) if ref_len else float(bucket["covered_bases_sum"]) + float(item["length"])
    entries = []
    for item in aggregate.values():
        entries.append(
            {
                "species_group": str(item["species_group"]),
                "isolate_label": str(item["isolate_label"]),
                "segment_label": ",".join(sorted(item["segment_hits"])),
                "coverage": float(item["coverage_sum"]),
                "depth_value": 0.0,
                "covered_bases": float(item["covered_bases_sum"]),
                "support_value": float(item["bitscore_sum"]),
                "reference_bases": float(item["reference_bases_sum"]),
                "reference_name": str(item["reference_name"]),
            }
        )
    ranked_rows = _rank_rotavirus_group_support(
        [
            {
                "species_group": str(item["group"]),
                "isolate_label": str(item["isolate"]),
                "segment_label": str(item["segment"]),
                "coverage": float(item["coverage"]),
                "depth_value": 0.0,
                "covered_bases": float(item["length"]),
                "support_value": float(item["bitscore"]),
                "reference_bases": float(ref_lengths.get(str(item["reference_name"]), 0)),
                "reference_name": str(item["reference_name"]),
            }
            for item in by_group_isolate.values()
        ]
    )
    best = ranked_rows[0] if ranked_rows else {"species_group": "", "isolate_label": "", "reference_name": "", "segment_count": "0", "coverage_sum": "", "full_length_coverage_pct": "", "depth_sum": "", "covered_bases_sum": "", "reference_bases_sum": "", "num_reads_sum": "", "selected_group_type": "", "is_selected": "no"}
    best["method"] = "blastn"
    return {"hit": best, "rows": ranked_rows}


def _run_rotavirus_group_typing(
    pre: str,
    query_fasta: Path | None = None,
    single_fastq: str = "",
    fq1: str = "",
    fq2: str = "",
    long_type: str = "",
    threads: int = 4,
    output_dir: Path | None = None,
    logf=None,
) -> dict[str, object]:
    typing_dir = (output_dir or (Path(f"{pre}_rotavirus_reference_selection") / "group_typing")).resolve()
    typing_dir.mkdir(parents=True, exist_ok=True)
    db_fasta = typing_dir / "rotavirus_group_panel.fasta"
    panel_rows = _write_rotavirus_group_panel_fasta(db_fasta)
    if not db_fasta.is_file() or db_fasta.stat().st_size == 0 or not panel_rows:
        return {"species_group": "", "isolate_label": "", "summary_path": "", "hit": {"method": "missing_db"}, "rows": []}
    use_blast = query_fasta is not None and query_fasta.is_file() and query_fasta.stat().st_size > 0
    result = (
        _run_rotavirus_group_blast_typing(query_fasta, db_fasta, typing_dir, logf=logf)
        if use_blast
        else _run_rotavirus_group_read_typing(
            db_fasta,
            typing_dir,
            single_fastq=single_fastq,
            fq1=fq1,
            fq2=fq2,
            long_type=long_type,
            threads=threads,
            logf=logf,
        )
    )
    hit = dict(result.get("hit") or {})
    all_rows = list(result.get("rows") or [])
    summary_path = typing_dir / "group_typing.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "group_type", "isolate_label", "reference_name", "segment_count", "coverage_sum", "full_length_coverage_pct", "depth_sum", "covered_bases_sum", "reference_bases_sum", "num_reads_sum", "selected_group_type", "is_selected", "method"])
        for row in all_rows:
            writer.writerow([
                pre,
                row.get("species_group", ""),
                row.get("isolate_label", ""),
                row.get("reference_name", ""),
                row.get("segment_count", ""),
                row.get("coverage_sum", ""),
                row.get("full_length_coverage_pct", ""),
                row.get("depth_sum", ""),
                row.get("covered_bases_sum", ""),
                row.get("reference_bases_sum", ""),
                row.get("num_reads_sum", ""),
                row.get("selected_group_type", ""),
                row.get("is_selected", "no"),
                hit.get("method", ""),
            ])
    return {
        "species_group": str(hit.get("species_group") or ""),
        "isolate_label": str(hit.get("isolate_label") or ""),
        "summary_path": str(summary_path.resolve()),
        "hit": hit,
        "rows": all_rows,
    }


def _load_rotavirus_a_subtype_subject_meta() -> dict[str, dict[str, str]]:
    subject_meta: dict[str, dict[str, str]] = {}
    for row in _load_rotavirus_manifest(_resolve_rotavirus_a_subtype_manifest()):
        accession = str(row.get("accession") or "").strip()
        accession_root = accession.split(".", 1)[0]
        header = str(row.get("header") or "").strip()
        combo = str(row.get("combo") or "").strip()
        g_genotype, p_genotype = _parse_rotavirus_combo(combo)
        meta = {
            "combo": combo,
            "g_genotype": str(row.get("g_genotype") or "").strip() or g_genotype,
            "p_genotype": str(row.get("p_genotype") or "").strip() or p_genotype,
            "segment_label": str(row.get("segment_label") or "").strip().upper(),
            "isolate_label": str(row.get("isolate_label") or "").strip(),
            "accession": accession,
            "accession_root": accession_root,
            "header": header,
        }
        for key in {accession, accession_root, header, _sanitize_tree_label(accession_root), _sanitize_tree_label(accession)}:
            if key:
                subject_meta[key] = meta
    return subject_meta


def _write_rotavirus_a_subtype_panel_fasta(output_path: Path) -> list[dict[str, str]]:
    panel_rows = [
        row for row in _load_rotavirus_manifest(_resolve_rotavirus_a_subtype_manifest())
        if str(row.get("segment_label") or "").strip().upper() in {"VP4", "VP7"}
    ]
    sequence_map = _build_rhinovirus_sequence_map(_resolve_rotavirus_a_subtype_fasta())
    with output_path.open("w", encoding="utf-8") as handle:
        for row in panel_rows:
            header = str(row.get("header") or "").strip()
            accession = str(row.get("accession") or "").strip()
            record = sequence_map.get(header) or sequence_map.get(accession)
            if record is None:
                continue
            fasta_id = _sanitize_tree_label(accession.split(".", 1)[0] or accession)
            handle.write(f">{fasta_id}\n")
            sequence = str(record.seq)
            for index in range(0, len(sequence), 80):
                handle.write(sequence[index:index + 80] + "\n")
    return panel_rows


def _score_rotavirus_a_subtype_rows(rows: list[dict[str, object]]) -> dict[str, str]:
    g_hits: dict[str, tuple[float, float, float, float]] = {}
    p_hits: dict[str, tuple[float, float, float, float]] = {}
    combo_hits: dict[str, dict[str, object]] = {}
    for row in rows:
        meta = row.get("_meta") if isinstance(row.get("_meta"), dict) else {}
        combo = str(meta.get("combo") or "").strip()
        segment = str(meta.get("segment_label") or "").strip().upper()
        g_genotype = str(meta.get("g_genotype") or "").strip()
        p_genotype = str(meta.get("p_genotype") or "").strip()
        if not combo or segment not in {"VP4", "VP7"}:
            continue
        score = (
            float(row.get("coverage") or 0.0),
            float(row.get("mean_depth") or 0.0),
            float(row.get("covered_bases") or 0.0),
            float(row.get("num_reads") or 0.0),
        )
        if segment == "VP7" and g_genotype:
            if g_genotype not in g_hits or score > g_hits[g_genotype]:
                g_hits[g_genotype] = score
        if segment == "VP4" and p_genotype:
            if p_genotype not in p_hits or score > p_hits[p_genotype]:
                p_hits[p_genotype] = score
        bucket = combo_hits.setdefault(
            combo,
            {
                "combo": combo,
                "g_genotype": g_genotype,
                "p_genotype": p_genotype,
                "vp4_score": (0.0, 0.0, 0.0, 0.0),
                "vp7_score": (0.0, 0.0, 0.0, 0.0),
                "isolate_label": str(meta.get("isolate_label") or ""),
            },
        )
        if segment == "VP4" and score > bucket["vp4_score"]:
            bucket["vp4_score"] = score
        if segment == "VP7" and score > bucket["vp7_score"]:
            bucket["vp7_score"] = score
    best_g = max(g_hits.items(), key=lambda item: item[1], default=("", (0, 0, 0, 0)))[0]
    best_p = max(p_hits.items(), key=lambda item: item[1], default=("", (0, 0, 0, 0)))[0]
    ranked_combos = sorted(
        combo_hits.values(),
        key=lambda item: (
            1 if item["g_genotype"] == best_g and item["p_genotype"] == best_p and best_g and best_p else 0,
            item["vp7_score"][0],
            item["vp4_score"][0],
            item["vp7_score"][1],
            item["vp4_score"][1],
            item["vp7_score"][2] + item["vp4_score"][2],
            item["vp7_score"][3] + item["vp4_score"][3],
        ),
        reverse=True,
    )
    best_combo = ranked_combos[0] if ranked_combos else {"combo": "", "g_genotype": best_g, "p_genotype": best_p, "isolate_label": ""}
    return {
        "combo": str(best_combo.get("combo") or ""),
        "g_genotype": str(best_combo.get("g_genotype") or best_g or ""),
        "p_genotype": str(best_combo.get("p_genotype") or best_p or ""),
        "isolate_label": str(best_combo.get("isolate_label") or ""),
    }


def _run_rotavirus_a_subtype_read_typing(
    db_fasta: Path,
    out_dir: Path,
    single_fastq: str = "",
    fq1: str = "",
    fq2: str = "",
    long_type: str = "",
    threads: int = 4,
    logf=None,
) -> dict[str, str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    bam_path = out_dir / "subtype.screening.bam"
    coverage_path = out_dir / "subtype.coverage.tsv"
    if fq1:
        map_cmd = " ".join(
            [
                "minimap2",
                "-ax",
                "sr",
                shlex.quote(str(db_fasta)),
                shlex.quote(str(fq1)),
                *([shlex.quote(str(fq2))] if fq2 else []),
                "|",
                "samtools",
                "sort",
                "-o",
                shlex.quote(str(bam_path)),
            ]
        )
    else:
        map_cmd = " ".join(
            [
                "minimap2",
                "-ax",
                _choose_minimap2_preset(long_type),
                shlex.quote(str(db_fasta)),
                shlex.quote(str(single_fastq)),
                "-t",
                str(max(1, int(threads or 1))),
                "|",
                "samtools",
                "sort",
                "-o",
                shlex.quote(str(bam_path)),
            ]
        )
    run_command(map_cmd, logf=logf)
    run_command(f"samtools index -@ {max(1, int(threads or 1))} {shlex.quote(str(bam_path))}", logf=logf)
    run_command(f"samtools coverage {shlex.quote(str(bam_path))} > {shlex.quote(str(coverage_path))}", logf=logf)
    subject_meta = _load_rotavirus_a_subtype_subject_meta()
    coverage_rows = _parse_samtools_coverage_rows(coverage_path)
    enriched_rows = []
    for row in coverage_rows:
        subject = str(row.get("reference_name") or "").strip()
        meta = subject_meta.get(subject) or subject_meta.get(subject.split("_", 1)[0]) or {}
        enriched = dict(row)
        enriched["_meta"] = meta
        enriched_rows.append(enriched)
    hit = _score_rotavirus_a_subtype_rows(enriched_rows)
    hit["method"] = "read_coverage"
    return hit


def _run_rotavirus_a_subtype_blast_typing(query_fasta: Path, db_fasta: Path, out_dir: Path, logf=None) -> dict[str, str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    out_tsv = out_dir / "subtype.blastn.tsv"
    blastn_bin = _resolve_hadv_blastn_bin()
    cmd = " ".join(
        [
            shlex.quote(blastn_bin),
            "-task blastn",
            "-query",
            shlex.quote(str(query_fasta)),
            "-subject",
            shlex.quote(str(db_fasta)),
            "-outfmt",
            shlex.quote("6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"),
            "-evalue 1e-20",
            "-max_target_seqs 100",
            "-dust no",
            "-num_threads 4",
            "-out",
            shlex.quote(str(out_tsv)),
        ]
    )
    run_command(cmd, logf=logf)
    if not out_tsv.is_file() or out_tsv.stat().st_size == 0:
        return {"combo": "", "g_genotype": "", "p_genotype": "", "isolate_label": "", "method": "blastn"}
    subject_meta = _load_rotavirus_a_subtype_subject_meta()
    ref_lengths = {str(record.id): len(str(record.seq)) for record in SeqIO.parse(str(db_fasta), "fasta")}
    by_reference: dict[str, dict[str, object]] = {}
    with out_tsv.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            _qseqid, sseqid, pident, length, _mismatch, _gapopen, _qstart, _qend, _sstart, _send, _evalue, bitscore = parts[:12]
            meta = subject_meta.get(sseqid) or subject_meta.get(sseqid.split("_", 1)[0]) or {}
            if str(meta.get("segment_label") or "").strip().upper() not in {"VP4", "VP7"}:
                continue
            ref_len = ref_lengths.get(sseqid, 0)
            try:
                qcov_ref = min(100.0, (float(length) / ref_len) * 100) if ref_len else 0.0
                score = (float(bitscore), float(pident), qcov_ref, float(length))
            except ValueError:
                continue
            bucket = by_reference.setdefault(
                sseqid,
                {
                    "reference_name": sseqid,
                    "coverage": qcov_ref,
                    "mean_depth": float(pident),
                    "covered_bases": float(length),
                    "num_reads": float(bitscore),
                    "_meta": meta,
                },
            )
            if score > (float(bucket["num_reads"]), float(bucket["mean_depth"]), float(bucket["coverage"]), float(bucket["covered_bases"])):
                bucket["coverage"] = qcov_ref
                bucket["mean_depth"] = float(pident)
                bucket["covered_bases"] = float(length)
                bucket["num_reads"] = float(bitscore)
    hit = _score_rotavirus_a_subtype_rows(list(by_reference.values()))
    hit["method"] = "blastn"
    return hit


def _run_rotavirus_a_subtype_typing(
    pre: str,
    query_fasta: Path | None = None,
    single_fastq: str = "",
    fq1: str = "",
    fq2: str = "",
    long_type: str = "",
    threads: int = 4,
    output_dir: Path | None = None,
    logf=None,
) -> dict[str, object]:
    typing_dir = (output_dir or (Path(f"{pre}_rotavirus_reference_selection") / "subtype_typing")).resolve()
    typing_dir.mkdir(parents=True, exist_ok=True)
    db_fasta = typing_dir / "rotavirus_a_vp4_vp7_panel.fasta"
    panel_rows = _write_rotavirus_a_subtype_panel_fasta(db_fasta)
    if not db_fasta.is_file() or db_fasta.stat().st_size == 0 or not panel_rows:
        return {"combo": "", "g_genotype": "", "p_genotype": "", "summary_path": "", "hit": {"method": "missing_db"}}
    use_blast = query_fasta is not None and query_fasta.is_file() and query_fasta.stat().st_size > 0
    hit = (
        _run_rotavirus_a_subtype_blast_typing(query_fasta, db_fasta, typing_dir, logf=logf)
        if use_blast
        else _run_rotavirus_a_subtype_read_typing(
            db_fasta,
            typing_dir,
            single_fastq=single_fastq,
            fq1=fq1,
            fq2=fq2,
            long_type=long_type,
            threads=threads,
            logf=logf,
        )
    )
    summary_path = typing_dir / "subtype_typing.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "group_type", "g_genotype", "p_genotype", "subtype_combo", "isolate_label", "method"])
        writer.writerow([pre, "A", hit.get("g_genotype", ""), hit.get("p_genotype", ""), hit.get("combo", ""), hit.get("isolate_label", ""), hit.get("method", "")])
    return {
        "combo": str(hit.get("combo") or ""),
        "g_genotype": str(hit.get("g_genotype") or ""),
        "p_genotype": str(hit.get("p_genotype") or ""),
        "isolate_label": str(hit.get("isolate_label") or ""),
        "summary_path": str(summary_path.resolve()),
        "hit": hit,
    }


def resolve_rotavirus_reference(
    pre: str,
    species: str = "",
    requested_ref: str = "",
    single_fastq: str = "",
    fq1: str = "",
    fq2: str = "",
    long_type: str = "",
    threads: int = 4,
    query_fasta: str = "",
    logf=None,
) -> dict[str, str]:
    requested = str(requested_ref or "").strip()
    explicit_group = _infer_rotavirus_group_from_label(species)
    if requested and Path(requested).is_file():
        return {
            "status": "ready",
            "group_type": explicit_group,
            "species_label": _rotavirus_species_label(explicit_group),
            "g_genotype": "",
            "p_genotype": "",
            "subtype_combo": "",
            "isolate_label": "",
            "reference_path": str(Path(requested).resolve()),
            "gff_path": "nogtf",
            "summary_path": "",
            "group_summary_path": "",
            "subtype_summary_path": "",
        }
    screening_dir = Path(f"{pre}_rotavirus_reference_selection")
    screening_dir.mkdir(exist_ok=True)
    query_path = Path(query_fasta) if str(query_fasta or "").strip() else None
    group_result = _run_rotavirus_group_typing(
        pre,
        query_fasta=query_path,
        single_fastq=single_fastq,
        fq1=fq1,
        fq2=fq2,
        long_type=long_type,
        threads=threads,
        output_dir=screening_dir / "group_typing",
        logf=logf,
    )
    group_type = str(group_result.get("species_group") or explicit_group or "").strip().upper()
    summary_path = screening_dir / "selection.tsv"
    if group_type not in {"A", "B", "C"}:
        with summary_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sample", "group_type", "status", "note", "group_summary_path"])
            writer.writerow([pre, group_type, "missing", f"未能确定轮状病毒大型别，species={species or '-'}", str(group_result.get("summary_path") or "")])
        return {
            "status": "missing",
            "group_type": group_type,
            "species_label": species or "Rotavirus",
            "g_genotype": "",
            "p_genotype": "",
            "subtype_combo": "",
            "isolate_label": "",
            "reference_path": "",
            "gff_path": "nogtf",
            "summary_path": str(summary_path.resolve()),
            "group_summary_path": str(group_result.get("summary_path") or ""),
            "subtype_summary_path": "",
        }
    if group_type == "A":
        subtype_result = _run_rotavirus_a_subtype_typing(
            pre,
            query_fasta=query_path,
            single_fastq=single_fastq,
            fq1=fq1,
            fq2=fq2,
            long_type=long_type,
            threads=threads,
            output_dir=screening_dir / "subtype_typing",
            logf=logf,
        )
        subtype_combo = str(subtype_result.get("combo") or "").strip()
        parsed_g, parsed_p = _parse_rotavirus_combo(subtype_combo)
        summary_rows = _load_rotavirus_manifest(_resolve_rotavirus_a_subtype_summary())
        selected_summary = next((row for row in summary_rows if str(row.get("combo") or "").strip() == subtype_combo), {})
        reference_path = ""
        if selected_summary:
            fasta_file = str(selected_summary.get("fasta_file") or "").strip()
            candidate_path = (_resolve_rotavirus_db_dir() / "subtype_complete_genomes" / fasta_file).resolve()
            if candidate_path.is_file():
                reference_path = str(candidate_path)
        with summary_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sample", "group_type", "g_genotype", "p_genotype", "subtype_combo", "isolate_label", "reference_path", "gff_path", "group_summary_path", "subtype_summary_path"])
            writer.writerow([
                pre,
                group_type,
                subtype_result.get("g_genotype", "") or parsed_g,
                subtype_result.get("p_genotype", "") or parsed_p,
                subtype_combo,
                selected_summary.get("isolate_label", subtype_result.get("isolate_label", "")),
                reference_path,
                "nogtf",
                str(group_result.get("summary_path") or ""),
                str(subtype_result.get("summary_path") or ""),
            ])
        return {
            "status": "ready" if reference_path else "missing",
            "group_type": group_type,
            "species_label": _rotavirus_species_label(group_type),
            "g_genotype": str(subtype_result.get("g_genotype") or parsed_g or ""),
            "p_genotype": str(subtype_result.get("p_genotype") or parsed_p or ""),
            "subtype_combo": subtype_combo,
            "isolate_label": str(selected_summary.get("isolate_label") or subtype_result.get("isolate_label") or ""),
            "reference_path": reference_path,
            "gff_path": "nogtf",
            "summary_path": str(summary_path.resolve()),
            "group_summary_path": str(group_result.get("summary_path") or ""),
            "subtype_summary_path": str(subtype_result.get("summary_path") or ""),
        }
    if group_type == "B":
        reference_path = str(_resolve_rotavirus_group_bc_fasta()) if _resolve_rotavirus_group_bc_fasta().is_file() else ""
        isolate_label = str(group_result.get("isolate_label") or "Bang373")
    else:
        full_rows = _load_rotavirus_manifest(_resolve_rotavirus_full_manifest())
        selected_rows = _choose_rotavirus_complete_isolate(full_rows, target_group="C")
        reference_path = ""
        isolate_label = str(group_result.get("isolate_label") or "")
        if selected_rows:
            full_fasta = _resolve_rotavirus_full_fasta()
            sequence_map = _build_rhinovirus_sequence_map(full_fasta)
            c_reference_fasta = screening_dir / "rotavirus_c.reference.fasta"
            with c_reference_fasta.open("w", encoding="utf-8") as handle:
                for row in selected_rows:
                    header = str(row.get("header") or "").strip()
                    accession = str(row.get("accession") or "").strip()
                    record = sequence_map.get(header) or sequence_map.get(accession)
                    if record is None:
                        continue
                    handle.write(f">{accession}\n")
                    sequence = str(record.seq)
                    for index in range(0, len(sequence), 80):
                        handle.write(sequence[index:index + 80] + "\n")
            if c_reference_fasta.is_file() and c_reference_fasta.stat().st_size > 0:
                reference_path = str(c_reference_fasta.resolve())
                isolate_label = str(selected_rows[0].get("isolate_label") or isolate_label)
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "group_type", "g_genotype", "p_genotype", "subtype_combo", "isolate_label", "reference_path", "gff_path", "group_summary_path", "subtype_summary_path"])
        writer.writerow([pre, group_type, "", "", "", isolate_label, reference_path, "nogtf", str(group_result.get("summary_path") or ""), ""])
    return {
        "status": "ready" if reference_path else "missing",
        "group_type": group_type,
        "species_label": _rotavirus_species_label(group_type),
        "g_genotype": "",
        "p_genotype": "",
        "subtype_combo": "",
        "isolate_label": isolate_label,
        "reference_path": reference_path,
        "gff_path": "nogtf",
        "summary_path": str(summary_path.resolve()),
        "group_summary_path": str(group_result.get("summary_path") or ""),
        "subtype_summary_path": "",
    }


def run_rotavirus_consensus_typing(pre: str, consensus_fasta: Path, logf=None) -> dict[str, str]:
    if not consensus_fasta.is_file() or consensus_fasta.stat().st_size == 0:
        return {"status": "missing", "summary_path": "", "group_type": "", "g_genotype": "", "p_genotype": "", "subtype_combo": ""}
    selection = resolve_rotavirus_reference(
        pre,
        species="Rotavirus",
        requested_ref="",
        query_fasta=str(consensus_fasta),
        logf=logf,
    )
    summary_path = Path(f"{pre}_rotavirus_reference_selection") / "consensus_typing.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "group_type", "g_genotype", "p_genotype", "subtype_combo", "selection_summary_path"])
        writer.writerow([
            pre,
            selection.get("group_type", ""),
            selection.get("g_genotype", ""),
            selection.get("p_genotype", ""),
            selection.get("subtype_combo", ""),
            selection.get("summary_path", ""),
        ])
    return {
        "status": "ready",
        "summary_path": str(summary_path.resolve()),
        "group_type": str(selection.get("group_type") or ""),
        "g_genotype": str(selection.get("g_genotype") or ""),
        "p_genotype": str(selection.get("p_genotype") or ""),
        "subtype_combo": str(selection.get("subtype_combo") or ""),
    }


def _resolve_seasonal_hcov_db_dir() -> Path:
    env_root = str(os.environ.get("META_SEASONAL_HCOV_DB_DIR") or "").strip()
    if env_root:
        return Path(env_root).expanduser().resolve()
    database_root = str(os.environ.get("META_DATABASE_ROOT") or "").strip()
    project_root = _project_root()
    candidates: list[Path] = []
    if database_root:
        candidates.append((Path(database_root) / "virus" / "seasonal_coronavirus").expanduser())
    candidates.extend(
        [
            project_root / "database" / "virus" / "seasonal_coronavirus",
            Path("/data/deploy/meta_genome/database/virus/seasonal_coronavirus"),
        ]
    )
    return next((candidate.resolve() for candidate in candidates if candidate.exists()), candidates[0].resolve())


def _resolve_seasonal_hcov_reference_manifest() -> Path:
    return (_resolve_seasonal_hcov_db_dir() / "reference_genomes" / "seasonal_hcov_reference_genomes_manifest.tsv").resolve()


def _resolve_seasonal_hcov_reference_fasta() -> Path:
    return (_resolve_seasonal_hcov_db_dir() / "reference_genomes" / "seasonal_hcov_reference_genomes.fasta").resolve()


def _resolve_seasonal_hcov_type_dir(hcov_type: str) -> Path:
    mapping = {
        "HCoV-229E": "HCoV_229E_genomes",
        "HCoV-NL63": "HCoV_NL63_genomes",
        "HCoV-OC43": "HCoV_OC43_genomes",
        "HCoV-HKU1": "HCoV_HKU1_genomes",
    }
    directory_name = mapping.get(str(hcov_type or "").strip(), "")
    return (_resolve_seasonal_hcov_db_dir() / directory_name).resolve() if directory_name else _resolve_seasonal_hcov_db_dir()


def _resolve_seasonal_hcov_reference_gff(hcov_type: str) -> Path:
    type_dir = _resolve_seasonal_hcov_type_dir(hcov_type)
    accession_map = {
        "HCoV-229E": "NC_002645.1",
        "HCoV-NL63": "NC_005831.2",
        "HCoV-OC43": "NC_006213.1",
        "HCoV-HKU1": "NC_006577.2",
    }
    accession = accession_map.get(str(hcov_type or "").strip(), "")
    return (type_dir / f"{accession}.gff3").resolve() if accession else type_dir / "missing.gff3"


def _resolve_seasonal_hcov_vadr_model_dir(project_root: Path, hcov_type: str) -> Path | None:
    mapping = {
        "HCoV-229E": "229E",
        "HCoV-NL63": "NL63",
        "HCoV-OC43": "OC43",
        "HCoV-HKU1": "HKU1",
    }
    model_key = mapping.get(str(hcov_type or "").strip(), "")
    if not model_key:
        return None
    model_dir = (project_root / "soft" / "vadr-models-hcov" / model_key).resolve()
    return model_dir if model_dir.is_dir() else None


def _load_seasonal_hcov_reference_manifest() -> list[dict[str, str]]:
    return _load_rhinovirus_manifest(_resolve_seasonal_hcov_reference_manifest())


def _infer_seasonal_hcov_type_from_ref_text(text: str) -> str:
    value = str(text or "").strip()
    if not value:
        return ""
    type_from_label = _infer_seasonal_hcov_type_from_label(value)
    if type_from_label:
        return type_from_label
    normalized = value.upper()
    accession_map = {
        "NC_002645": "HCoV-229E",
        "NC_005831": "HCoV-NL63",
        "NC_006213": "HCoV-OC43",
        "NC_006577": "HCoV-HKU1",
    }
    for accession_root, type_label in accession_map.items():
        if accession_root in normalized:
            return type_label
    return ""


def _run_vadr_seasonal_hcov_annotation(pre: str, input_fasta: Path, output_root: Path, hcov_type: str, logf=None) -> dict[str, str]:
    project_root = _project_root()
    model_dir = _resolve_seasonal_hcov_vadr_model_dir(project_root, hcov_type)
    if model_dir is None:
        return {"status": "missing", "note": f"未找到 {hcov_type or '?'} VADR 模型目录", "output_dir": "", "gff_path": ""}
    env = _build_vadr_env(project_root, model_dir)
    if env is None:
        return {"status": "missing", "note": "VADR 运行依赖不完整", "output_dir": "", "gff_path": ""}
    seasonal_env = dict(env)
    alt_bio_easel_dir = project_root / "soft" / "Bio-Easel-ncov"
    alt_bio_lib = alt_bio_easel_dir / "blib" / "lib"
    alt_bio_arch = alt_bio_easel_dir / "blib" / "arch"
    if alt_bio_lib.is_dir() and alt_bio_arch.is_dir():
        seasonal_env["VADRBIOEASELDIR"] = str(alt_bio_easel_dir)
        seasonal_env["PERL5LIB"] = os.pathsep.join(
            [
                str(project_root / "soft" / "vadr"),
                str(project_root / "soft" / "sequip"),
                str(alt_bio_lib),
                str(alt_bio_arch),
            ]
        )
    ncov_bin_dir = conda_env_path("ncov", "bin")
    seasonal_env["PATH"] = os.pathsep.join(
        [
            ncov_bin_dir,
            "/usr/bin",
            "/bin",
            str(seasonal_env.get("PATH") or ""),
        ]
    )
    perl_candidates = [f"{ncov_bin_dir}/perl", "/usr/bin/perl", _resolve_vadr_perl_bin(), "perl"]
    perl_bin = _resolve_working_perl_with_module(seasonal_env, "Bio::Easel::MSA", perl_candidates)
    vadr_root = output_root / "vadr"
    vadr_root.mkdir(parents=True, exist_ok=True)
    output_dir = vadr_root / f"{pre}_vadr"
    prefix = f"{output_dir.name}.vadr"
    pass_tbl = output_dir / f"{prefix}.pass.tbl"
    fail_tbl = output_dir / f"{prefix}.fail.tbl"
    gff_path = vadr_root / f"{pre}.vadr.gff3"
    if output_dir.is_dir() and (pass_tbl.is_file() or fail_tbl.is_file()) and gff_path.is_file() and gff_path.stat().st_size > 0:
        return {"status": "ready", "note": "已检测到现有季节性冠状病毒 VADR 注释结果", "output_dir": str(output_dir.resolve()), "gff_path": str(gff_path.resolve())}
    model_key = model_dir.name
    extra_args = ["-s", "--glsearch", "-r", "-f"]
    if model_key in {"HKU1", "OC43"}:
        extra_args.extend(["--alt_pass", "discontn"])
    vadr_script = project_root / "soft" / "vadr" / "v-annotate.pl"
    cmd_parts = [
        shlex.quote(str(perl_bin)),
        shlex.quote(str(vadr_script)),
        *extra_args,
        "--mkey",
        model_key,
        "--mdir",
        shlex.quote(str(model_dir)),
        shlex.quote(str(input_fasta)),
        shlex.quote(str(output_dir)),
    ]
    try:
        run_command(" ".join(cmd_parts), logf=logf, env=seasonal_env)
    except Exception as exc:
        return {"status": "failed", "note": f"季节性冠状病毒 VADR 注释失败: {exc}", "output_dir": str(output_dir.resolve()) if output_dir.exists() else "", "gff_path": ""}
    annotate_tbl2gff = project_root / "soft" / "vadr" / "miniscripts" / "annotate-tbl2gff.pl"
    source_tbl = pass_tbl if pass_tbl.is_file() else fail_tbl
    if source_tbl.is_file():
        try:
            run_command(
                f"{shlex.quote(str(perl_bin))} {shlex.quote(str(annotate_tbl2gff))} {shlex.quote(str(source_tbl))} > {shlex.quote(str(gff_path))}",
                logf=logf,
                env=seasonal_env,
            )
        except Exception:
            pass
    return {
        "status": "ready" if gff_path.is_file() and gff_path.stat().st_size > 0 else "failed",
        "note": "季节性冠状病毒 VADR 注释完成" if gff_path.is_file() and gff_path.stat().st_size > 0 else "季节性冠状病毒 VADR 未生成可用 GFF",
        "output_dir": str(output_dir.resolve()),
        "gff_path": str(gff_path.resolve()) if gff_path.is_file() and gff_path.stat().st_size > 0 else "",
    }


def prepare_seasonal_hcov_reference_annotation(pre: str, reference_fasta: Path, output_root: Path, hcov_type: str, logf=None) -> Path | None:
    result = _run_vadr_seasonal_hcov_annotation(pre, reference_fasta, output_root, hcov_type, logf=logf)
    gff_path = str(result.get("gff_path") or "").strip()
    if result.get("status") == "ready" and gff_path:
        gff_file = Path(gff_path)
        if gff_file.is_file() and gff_file.stat().st_size > 0:
            return gff_file
    return None


def prepare_seasonal_hcov_sample_annotation(pre: str, sample_fasta: Path, output_root: Path, hcov_type: str, logf=None) -> Path | None:
    result = _run_vadr_seasonal_hcov_annotation(pre, sample_fasta, output_root, hcov_type, logf=logf)
    gff_path = str(result.get("gff_path") or "").strip()
    if result.get("status") == "ready" and gff_path:
        gff_file = Path(gff_path)
        if gff_file.is_file() and gff_file.stat().st_size > 0:
            return gff_file
    return None


def _find_seasonal_hcov_spike_feature(gff_path: Path) -> dict[str, object] | None:
    if not gff_path.is_file():
        return None
    candidates: list[tuple[int, int, dict[str, object]]] = []
    with gff_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, _source, feature_type, start_text, end_text, _score, strand, _phase, attrs_text = parts[:9]
            attrs = _parse_norovirus_gff_attributes(attrs_text)
            gene_name = str(attrs.get("gene") or attrs.get("Name") or "").strip().upper()
            product = str(attrs.get("product") or attrs.get("Note") or attrs.get("note") or "").strip().upper()
            label = " ".join(
                [
                    gene_name,
                    product,
                    str(attrs.get("gbkey") or "").strip().upper(),
                ]
            )
            if feature_type not in {"CDS", "gene", "sequence_feature"}:
                continue
            if gene_name not in {"S", "SPIKE"} and "SPIKE" not in label:
                continue
            try:
                start = int(start_text)
                end = int(end_text)
            except ValueError:
                continue
            priority = {"CDS": 0, "gene": 1, "sequence_feature": 2}.get(feature_type, 9)
            candidates.append(
                (
                    priority,
                    -(end - start + 1),
                    {
                        "seqid": seqid,
                        "start": start,
                        "end": end,
                        "strand": strand,
                        "feature_type": feature_type,
                        "gene_name": gene_name or "S",
                        "product": product or "spike",
                    },
                )
            )
    if not candidates:
        return None
    candidates.sort(key=lambda item: (item[0], item[1], item[2]["start"], item[2]["end"]))
    return candidates[0][2]


def _extract_seasonal_hcov_spike_region(sample_fasta: Path, gff_path: Path, sample_label: str) -> tuple[SeqRecord | None, dict[str, object]]:
    feature = _find_seasonal_hcov_spike_feature(gff_path)
    if feature is None:
        return None, {"status": "missing", "note": "GFF 中未找到 S/spike 区域"}
    records = _read_fasta_records(sample_fasta)
    if not records:
        return None, {"status": "missing", "note": f"{sample_fasta.name} 中未找到序列记录"}
    seqid = str(feature.get("seqid") or "").strip()
    record = next((item for item in records if str(item.id).strip() == seqid), None)
    if record is None:
        record = max(records, key=lambda item: len(str(item.seq)))
    start = int(feature["start"])
    end = int(feature["end"])
    sequence = str(record.seq)
    if start < 1 or end > len(sequence) or start > end:
        return None, {"status": "missing", "note": "S/spike 坐标超出样本序列长度"}
    fragment = sequence[start - 1:end]
    if str(feature.get("strand") or "+") == "-":
        fragment = str(Seq(fragment).reverse_complement())
    fragment = fragment.strip().strip("Nn")
    if not fragment:
        return None, {"status": "missing", "note": "S/spike 区域仅包含空白或 N"}
    record_id = f"{sample_label}_spike_sample"
    return SeqRecord(Seq(fragment), id=record_id, description=f"{sample_label} spike sample"), {
        "status": "ready",
        "start": start,
        "end": end,
        "strand": str(feature.get("strand") or "+"),
        "seqid": str(record.id),
        "length": len(fragment),
    }


def _load_seasonal_hcov_spike_reference_records(hcov_type: str) -> tuple[list[dict[str, str]], list[SeqRecord]]:
    type_label = str(hcov_type or "").strip()
    if type_label == "HCoV-OC43":
        manifest_path = _resolve_seasonal_hcov_type_dir(type_label) / "HCoV_OC43_spike_genes_merged_with_references.tsv"
        if not manifest_path.is_file():
            return [], []
        rows = _load_rhinovirus_manifest(manifest_path)
        records: list[SeqRecord] = []
        ordered_rows: list[dict[str, str]] = []
        for row in rows:
            fasta_path = Path(str(row.get("fasta_path") or "")).expanduser()
            accession = str(row.get("accession") or "").strip()
            if not fasta_path.is_file() or not accession:
                continue
            start = str(row.get("spike_start") or "").strip()
            end = str(row.get("spike_end") or "").strip()
            strand = str(row.get("spike_strand") or "+").strip() or "+"
            subtype = str(row.get("tree_genotype") or "").strip() or "unresolved"
            tree_label = str(row.get("tree_label") or row.get("region_strain") or accession).strip()
            extracted, _meta = _extract_seasonal_hcov_spike_region(fasta_path, Path(str(row.get("gff_path") or "")), _sanitize_tree_label(tree_label))
            if extracted is None:
                continue
            record_id = _sanitize_tree_label(f"{accession}_{subtype}_{tree_label}")
            extracted.id = record_id
            extracted.name = record_id
            extracted.description = tree_label
            records.append(extracted)
            ordered_rows.append(
                {
                    "gene": "S",
                    "subtype": subtype,
                    "accession": accession,
                    "backup_rank": str(len(ordered_rows) + 1),
                    "tree_label": tree_label,
                }
            )
        return ordered_rows, records
    if type_label == "HCoV-HKU1":
        manifest_path = _resolve_seasonal_hcov_type_dir(type_label) / "HCoV_HKU1_spike_genes.tsv"
        if not manifest_path.is_file():
            return [], []
        rows = _load_rhinovirus_manifest(manifest_path)
        records: list[SeqRecord] = []
        ordered_rows: list[dict[str, str]] = []
        for row in rows:
            fasta_path = Path(str(row.get("fasta_path") or "")).expanduser()
            accession = str(row.get("accession") or "").strip()
            if not fasta_path.is_file() or not accession:
                continue
            extracted, _meta = _extract_seasonal_hcov_spike_region(fasta_path, Path(str(row.get("gff_path") or "")), _sanitize_tree_label(str(row.get("tree_label") or accession)))
            if extracted is None:
                continue
            subtype = str(row.get("subtype") or "").strip() or "unresolved"
            tree_label = str(row.get("tree_label") or accession).strip()
            record_id = _sanitize_tree_label(f"{accession}_{subtype}_{tree_label}")
            extracted.id = record_id
            extracted.name = record_id
            extracted.description = tree_label
            records.append(extracted)
            ordered_rows.append(
                {
                    "gene": "S",
                    "subtype": subtype,
                    "accession": accession,
                    "backup_rank": str(len(ordered_rows) + 1),
                    "tree_label": tree_label,
                }
            )
        return ordered_rows, records
    if type_label == "HCoV-NL63":
        manifest_path = _resolve_seasonal_hcov_type_dir(type_label) / "HCoV_NL63_spike_genes.tsv"
        if not manifest_path.is_file():
            return [], []
        rows = _load_rhinovirus_manifest(manifest_path)
        records: list[SeqRecord] = []
        ordered_rows: list[dict[str, str]] = []
        for row in rows:
            fasta_path = _resolve_seasonal_hcov_type_dir(type_label) / f"{str(row.get('accession') or '').strip()}.fasta"
            accession = str(row.get("accession") or "").strip()
            if not fasta_path.is_file() or not accession:
                continue
            extracted, _meta = _extract_seasonal_hcov_spike_region(fasta_path, Path(str(row.get("gff_path") or "")), _sanitize_tree_label(accession))
            if extracted is None:
                continue
            subtype = accession
            record_id = _sanitize_tree_label(f"{accession}_{subtype}")
            extracted.id = record_id
            extracted.name = record_id
            extracted.description = str(row.get("fasta_header") or accession)
            records.append(extracted)
            ordered_rows.append(
                {
                    "gene": "S",
                    "subtype": subtype,
                    "accession": accession,
                    "backup_rank": str(len(ordered_rows) + 1),
                    "tree_label": str(row.get("fasta_header") or accession),
                }
            )
        return ordered_rows, records
    if type_label == "HCoV-229E":
        manifest_path = _resolve_seasonal_hcov_type_dir(type_label) / "HCoV_229E_sample_spike_tree_annotations.tsv"
        if not manifest_path.is_file():
            return [], []
        rows = _load_rhinovirus_manifest(manifest_path)
        records: list[SeqRecord] = []
        ordered_rows: list[dict[str, str]] = []
        for row in rows:
            sample_name = str(row.get("sample_name") or "").strip()
            subtype = str(row.get("genogroup") or "").strip() or "unresolved"
            accessions = [item.strip() for item in str(row.get("accessions") or "").split(",") if item.strip()]
            spike_accession = next((item for item in accessions if item.startswith("MT797") and 676 <= int(item[5:]) <= 716), "")
            if not spike_accession:
                continue
            fasta_path = _resolve_seasonal_hcov_type_dir(type_label) / f"{spike_accession}.fasta"
            if not fasta_path.is_file():
                continue
            records_in_file = _read_fasta_records(fasta_path)
            if not records_in_file:
                continue
            template = records_in_file[0]
            record_id = _sanitize_tree_label(f"{spike_accession}_{subtype}_{sample_name}")
            cloned = SeqRecord(template.seq, id=record_id, name=record_id, description=str(template.description))
            records.append(cloned)
            ordered_rows.append(
                {
                    "gene": "S",
                    "subtype": subtype,
                    "accession": spike_accession,
                    "backup_rank": str(len(ordered_rows) + 1),
                    "tree_label": str(row.get("tree_label") or sample_name),
                }
            )
        return ordered_rows, records
    return [], []


def _build_seasonal_hcov_spike_tree(
    hcov_type: str,
    sample_record: SeqRecord,
    backup_rows: list[dict[str, str]],
    backup_records: list[SeqRecord],
    out_dir: Path,
    logf=None,
) -> dict[str, str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    sample_fasta = out_dir / "spike.sample.fasta"
    backup_fasta = out_dir / "spike.backup_refs.fasta"
    combined_fasta = out_dir / "spike.combined.fasta"
    aligned_fasta = out_dir / "spike.aligned.fasta"
    tree_path = out_dir / "spike.tree.nwk"
    members_tsv = out_dir / "spike.members.tsv"
    nearest_tsv = out_dir / "spike.nearest.tsv"

    SeqIO.write([sample_record], str(sample_fasta), "fasta")
    SeqIO.write(backup_records, str(backup_fasta), "fasta")
    SeqIO.write([sample_record, *backup_records], str(combined_fasta), "fasta")
    with members_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["tree_label", "member_type", "gene", "subtype", "accession", "backup_rank"])
        writer.writerow([sample_record.id, "sample", "S", "", "", ""])
        for row, record in zip(backup_rows, backup_records):
            writer.writerow([record.id, "backup_ref", "S", str(row.get("subtype") or ""), str(row.get("accession") or ""), str(row.get("backup_rank") or "")])

    if not aligned_fasta.is_file() or aligned_fasta.stat().st_size == 0:
        run_command(
            f"{shlex.quote('mafft')} --retree 1 --maxiterate 0 --quiet {shlex.quote(str(combined_fasta))} > {shlex.quote(str(aligned_fasta))}",
            logf=logf,
        )
    alignment = AlignIO.read(str(aligned_fasta), "fasta")
    if len(alignment) < 2:
        return {"status": "missing", "note": "S 基因对齐序列不足，无法建树", "tree_path": ""}
    calculator = DistanceCalculator("identity")
    distance_matrix = calculator.get_distance(alignment)
    nearest_row: dict[str, str] = {}
    sample_id = sample_record.id
    best_distance = None
    for row, record in zip(backup_rows, backup_records):
        try:
            distance_value = float(distance_matrix[sample_id, record.id])
        except Exception:
            continue
        if best_distance is None or distance_value < best_distance:
            best_distance = distance_value
            nearest_row = {
                "sample": sample_id,
                "hcov_type": hcov_type,
                "nearest_label": record.id,
                "nearest_tree_label": str(row.get("tree_label") or record.description or record.id),
                "nearest_accession": str(row.get("accession") or ""),
                "nearest_subtype": str(row.get("subtype") or ""),
                "distance": f"{distance_value:.6f}",
            }
    with nearest_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["sample", "hcov_type", "nearest_label", "nearest_tree_label", "nearest_accession", "nearest_subtype", "distance"],
            delimiter="\t",
        )
        writer.writeheader()
        if nearest_row:
            writer.writerow(nearest_row)

    if len(alignment) == 2:
        labels = [record.id for record in alignment]
        tree_path.write_text(f"({labels[0]}:0.1,{labels[1]}:0.1);", encoding="utf-8")
    else:
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(distance_matrix)
        Phylo.write(tree, str(tree_path), "newick")
    return {
        "status": "ready",
        "tree_path": str(tree_path.resolve()),
        "aligned_fasta": str(aligned_fasta.resolve()),
        "combined_fasta": str(combined_fasta.resolve()),
        "members_tsv": str(members_tsv.resolve()),
        "sample_fasta": str(sample_fasta.resolve()),
        "backup_fasta": str(backup_fasta.resolve()),
        "nearest_tsv": str(nearest_tsv.resolve()),
        "nearest_subtype": str(nearest_row.get("nearest_subtype") or ""),
        "nearest_accession": str(nearest_row.get("nearest_accession") or ""),
        "nearest_tree_label": str(nearest_row.get("nearest_tree_label") or ""),
        "nearest_distance": str(nearest_row.get("distance") or ""),
        "member_count": str(len(backup_records) + 1),
        "backup_count": str(len(backup_records)),
    }


def build_seasonal_hcov_spike_phylogeny_assets(pre: str, sample_fasta: Path, gff_path: Path, hcov_type: str, logf=None) -> dict[str, object]:
    output_root = Path(f"{pre}_seasonal_hcov_reference_selection") / "phylogeny"
    output_root.mkdir(parents=True, exist_ok=True)
    summary_path = output_root / "summary.tsv"
    consensus_typing_path = Path(f"{pre}_seasonal_hcov_reference_selection") / "consensus_typing" / "consensus_typing.tsv"
    result_row = {
        "gene": "S",
        "hcov_type": hcov_type,
        "subtype": "",
        "status": "missing",
        "note": "",
        "sample_region_path": "",
        "backup_fasta": "",
        "tree_path": "",
        "member_count": "0",
        "nearest_accession": "",
        "nearest_tree_label": "",
        "nearest_distance": "",
    }
    sample_record, sample_meta = _extract_seasonal_hcov_spike_region(sample_fasta, gff_path, _sanitize_tree_label(pre))
    if sample_record is None:
        result_row["note"] = str(sample_meta.get("note") or "未提取到 S 基因")
    else:
        backup_rows, backup_records = _load_seasonal_hcov_spike_reference_records(hcov_type)
        if not backup_records:
            result_row["note"] = f"未找到 {hcov_type or '?'} 的 S 基因参考库"
        else:
            tree_result = _build_seasonal_hcov_spike_tree(hcov_type, sample_record, backup_rows, backup_records, output_root / "spike", logf=logf)
            result_row["status"] = str(tree_result.get("status") or "missing")
            result_row["note"] = str(tree_result.get("note") or "")
            result_row["sample_region_path"] = str((output_root / "spike" / "spike.sample.fasta").resolve())
            result_row["backup_fasta"] = str(tree_result.get("backup_fasta") or "")
            result_row["tree_path"] = str(tree_result.get("tree_path") or "")
            result_row["member_count"] = str(tree_result.get("member_count") or "0")
            result_row["subtype"] = str(tree_result.get("nearest_subtype") or "")
            result_row["nearest_accession"] = str(tree_result.get("nearest_accession") or "")
            result_row["nearest_tree_label"] = str(tree_result.get("nearest_tree_label") or "")
            result_row["nearest_distance"] = str(tree_result.get("nearest_distance") or "")
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "gene",
                "hcov_type",
                "subtype",
                "status",
                "note",
                "sample_region_path",
                "backup_fasta",
                "tree_path",
                "member_count",
                "nearest_accession",
                "nearest_tree_label",
                "nearest_distance",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(result_row)
    consensus_typing_path.parent.mkdir(parents=True, exist_ok=True)
    with consensus_typing_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "hcov_type", "subtype", "nearest_accession", "nearest_tree_label", "nearest_distance", "summary_path"])
        writer.writerow([
            pre,
            hcov_type,
            result_row["subtype"],
            result_row["nearest_accession"],
            result_row["nearest_tree_label"],
            result_row["nearest_distance"],
            str(summary_path.resolve()),
        ])
    return {
        "status": "ready" if result_row.get("status") == "ready" else "missing",
        "summary_path": str(summary_path.resolve()),
        "consensus_typing_path": str(consensus_typing_path.resolve()),
        "rows": [result_row],
    }


def resolve_seasonal_hcov_reference(
    pre: str,
    species: str = "",
    requested_ref: str = "",
    single_fastq: str = "",
    fq1: str = "",
    fq2: str = "",
    long_type: str = "",
    threads: int = 4,
    logf=None,
) -> dict[str, str]:
    requested = str(requested_ref or "").strip()
    requested_path = Path(requested).expanduser() if requested else None
    explicit_type = _infer_seasonal_hcov_type_from_label(species)
    if requested_path and requested_path.is_file():
        inferred_type = explicit_type or _infer_seasonal_hcov_type_from_ref_text(requested_path.name)
        gff_path = ""
        if inferred_type:
            reference_gff = _resolve_seasonal_hcov_reference_gff(inferred_type)
            if reference_gff.is_file():
                gff_path = str(reference_gff)
        return {
            "status": "ready" if inferred_type else "missing",
            "hcov_type": inferred_type,
            "species_label": _seasonal_hcov_species_label(inferred_type) if inferred_type else (species or "Human coronavirus"),
            "reference_path": str(requested_path.resolve()),
            "gff_path": gff_path or "nogtf",
            "summary_path": "",
        }

    screening_dir = Path(f"{pre}_seasonal_hcov_reference_selection")
    screening_dir.mkdir(exist_ok=True)
    manifest_rows = _load_seasonal_hcov_reference_manifest()
    fasta_path = _resolve_seasonal_hcov_reference_fasta()
    sequence_map = _build_rhinovirus_sequence_map(fasta_path)
    candidate_fasta = screening_dir / "candidate_references.fasta"
    selected_rows = manifest_rows[:]
    with candidate_fasta.open("w", encoding="utf-8") as handle:
        for row in selected_rows:
            header = str(row.get("header") or row.get("accession") or "").strip()
            record = sequence_map.get(header) or sequence_map.get(str(row.get("accession") or "").strip())
            if record is None:
                continue
            fasta_id = _sanitize_tree_label(str(row.get("accession") or record.id))
            handle.write(f">{fasta_id}\n")
            seq = str(record.seq)
            for index in range(0, len(seq), 80):
                handle.write(seq[index:index + 80] + "\n")
    summary_path = screening_dir / "selection.tsv"
    if not candidate_fasta.is_file() or candidate_fasta.stat().st_size == 0:
        with summary_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sample", "hcov_type", "status", "note"])
            writer.writerow([pre, explicit_type, "missing", "未找到季节性冠状病毒参考候选序列"])
        return {
            "status": "missing",
            "hcov_type": explicit_type,
            "species_label": _seasonal_hcov_species_label(explicit_type) if explicit_type else (species or "Human coronavirus"),
            "reference_path": "",
            "gff_path": "nogtf",
            "summary_path": str(summary_path.resolve()),
        }
    bam_path = screening_dir / "candidate_references.bam"
    coverage_path = screening_dir / "candidate_references.coverage.tsv"
    if fq1:
        map_cmd = " ".join(
            [
                "minimap2",
                "-ax",
                "sr",
                shlex.quote(str(candidate_fasta)),
                shlex.quote(str(fq1)),
                *([shlex.quote(str(fq2))] if fq2 else []),
                "|",
                "samtools",
                "sort",
                "-o",
                shlex.quote(str(bam_path)),
            ]
        )
    else:
        map_cmd = " ".join(
            [
                "minimap2",
                "-ax",
                _choose_minimap2_preset(long_type),
                shlex.quote(str(candidate_fasta)),
                shlex.quote(str(single_fastq)),
                "-t",
                str(max(1, int(threads or 1))),
                "|",
                "samtools",
                "sort",
                "-o",
                shlex.quote(str(bam_path)),
            ]
        )
    run_command(map_cmd, logf=logf)
    run_command(f"samtools index -@ {max(1, int(threads or 1))} {shlex.quote(str(bam_path))}", logf=logf)
    run_command(f"samtools coverage {shlex.quote(str(bam_path))} > {shlex.quote(str(coverage_path))}", logf=logf)
    coverage_rows = _parse_samtools_coverage_rows(coverage_path)
    candidate_meta = {_sanitize_tree_label(str(row.get('accession') or '')): row for row in selected_rows}
    best = None
    for row in coverage_rows:
        ref_name = str(row.get("reference_name") or "").strip()
        meta = candidate_meta.get(ref_name)
        if meta is None:
            continue
        score = (
            float(row.get("coverage") or 0.0),
            float(row.get("mean_depth") or 0.0),
            float(row.get("covered_bases") or 0.0),
            float(row.get("num_reads") or 0.0),
        )
        if best is None or score > best["score"]:
            best = {"score": score, "meta": meta, "coverage_row": row}
    if best is None:
        with summary_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sample", "hcov_type", "status", "note"])
            writer.writerow([pre, explicit_type, "missing", "未能根据四种季节性冠状病毒参考完成分型"])
        return {
            "status": "missing",
            "hcov_type": explicit_type,
            "species_label": _seasonal_hcov_species_label(explicit_type) if explicit_type else (species or "Human coronavirus"),
            "reference_path": "",
            "gff_path": "nogtf",
            "summary_path": str(summary_path.resolve()),
        }
    best_meta = best["meta"]
    best_row = best["coverage_row"]
    best_record = sequence_map.get(str(best_meta.get("header") or "").strip()) or sequence_map.get(str(best_meta.get("accession") or "").strip())
    best_type = str(best_meta.get("virus") or "").strip()
    best_reference_fasta = screening_dir / f"{_sanitize_tree_label(best_type or pre)}.reference.fasta"
    with best_reference_fasta.open("w", encoding="utf-8") as handle:
        handle.write(f">{str(best_meta.get('accession') or best_record.id if best_record else 'reference')}\n")
        sequence = str(best_record.seq if best_record is not None else "")
        for index in range(0, len(sequence), 80):
            handle.write(sequence[index:index + 80] + "\n")
    reference_gff = _resolve_seasonal_hcov_reference_gff(best_type)
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "hcov_type", "coverage", "mean_depth", "covered_bases", "num_reads", "reference_name", "reference_accession", "reference_path", "gff_path"])
        writer.writerow([
            pre,
            best_type,
            f"{float(best_row.get('coverage') or 0.0):.6f}",
            f"{float(best_row.get('mean_depth') or 0.0):.6f}",
            f"{float(best_row.get('covered_bases') or 0.0):.0f}",
            f"{float(best_row.get('num_reads') or 0.0):.0f}",
            str(best_meta.get("species") or best_meta.get("header") or best_meta.get("accession") or ""),
            str(best_meta.get("accession") or ""),
            str(best_reference_fasta.resolve()),
            str(reference_gff.resolve()) if reference_gff.is_file() else "nogtf",
        ])
    return {
        "status": "ready",
        "hcov_type": best_type,
        "species_label": _seasonal_hcov_species_label(best_type),
        "reference_path": str(best_reference_fasta.resolve()),
        "gff_path": str(reference_gff.resolve()) if reference_gff.is_file() else "nogtf",
        "summary_path": str(summary_path.resolve()),
    }


def _resolve_hadv_typing_tsv() -> Path:
    return (_resolve_hadv_db_dir() / "hadv_typing.tsv").resolve()


def _resolve_hadv_gff_path(accession: str) -> Path:
    return (_resolve_hadv_db_dir() / "reference_genomes" / "gff3" / f"{str(accession or '').split('.')[0]}.gff3").resolve()


def _resolve_norovirus_db_dir() -> Path:
    env_root = str(os.environ.get("META_NOROVIRUS_DB_DIR") or "").strip()
    if env_root:
        return Path(env_root).expanduser().resolve()
    database_root = str(os.environ.get("META_DATABASE_ROOT") or "").strip()
    project_root = Path(__file__).resolve().parents[1]
    candidates: list[Path] = []
    if database_root:
        candidates.append((Path(database_root) / "virus" / "norovirus").expanduser())
    candidates.extend([
        project_root / "database" / "virus" / "norovirus",
        Path("/data/deploy/meta_genome/database/virus/norovirus"),
    ])
    return next((candidate.resolve() for candidate in candidates if candidate.exists()), candidates[0].resolve())


def _resolve_norovirus_cdc_ref_fasta(gene_name: str) -> Path:
    normalized = str(gene_name or "").strip().lower()
    return (_resolve_norovirus_db_dir() / "cdc_typing_refs" / f"cdc_norovirus_{normalized}_refs.fasta").resolve()


def _resolve_norovirus_cdc_ref_manifest(gene_name: str) -> Path:
    normalized = str(gene_name or "").strip().lower()
    return (_resolve_norovirus_db_dir() / "cdc_typing_refs" / f"cdc_norovirus_{normalized}_refs.tsv").resolve()


def _resolve_norovirus_full_genome_fasta() -> Path:
    return (_resolve_norovirus_db_dir() / "full_genomes" / "human_norovirus_complete_genomes.fasta").resolve()


def _resolve_norovirus_full_genome_manifest() -> Path:
    return (_resolve_norovirus_db_dir() / "full_genomes" / "human_norovirus_complete_genomes_manifest.tsv").resolve()


def _resolve_norovirus_backup_ref_fasta(gene_name: str) -> Path:
    normalized = str(gene_name or "").strip().lower()
    return (_resolve_norovirus_db_dir() / "cdc_typing_refs" / "backup_refs" / f"cdc_norovirus_{normalized}_backup_refs.fasta").resolve()


def _resolve_norovirus_backup_ref_manifest(gene_name: str) -> Path:
    normalized = str(gene_name or "").strip().lower()
    return (_resolve_norovirus_db_dir() / "cdc_typing_refs" / "backup_refs" / f"cdc_norovirus_{normalized}_backup_refs.tsv").resolve()


def _resolve_norovirus_fragment_typing_tsv() -> Path:
    return (_resolve_norovirus_db_dir() / "full_genomes" / "noro_fragment_typing.tsv").resolve()


def _hadv_species_label(type_label: str = "") -> str:
    normalized = str(type_label or "").strip()
    match = re.search(r"HAdV-([A-G])([0-9]{1,3})", normalized, flags=re.IGNORECASE)
    normalized = f"HAdV-{match.group(1).upper()}{match.group(2)}" if match else normalized
    return f"Human adenovirus {normalized}" if normalized else "Human adenovirus"


def _extract_hadv_type_label(text: str) -> str:
    value = str(text or "")
    matched = re.search(r"HAdV-([A-G])\s*([0-9]{1,3})", value, flags=re.IGNORECASE)
    if matched:
        return f"HAdV-{matched.group(1).upper()}{matched.group(2)}"
    matched = re.search(r"HAdV\s*([A-G])\s*([0-9]{1,3})", value, flags=re.IGNORECASE)
    if matched:
        return f"HAdV-{matched.group(1).upper()}{matched.group(2)}"
    matched = re.search(r"Human adenovirus(?:\s+type)?\s*([A-G])\s*([0-9]{1,3})", value, flags=re.IGNORECASE)
    if matched:
        return f"HAdV-{matched.group(1).upper()}{matched.group(2)}"
    matched = re.search(r"Human mastadenovirus\s+([A-G]).*?(?:type\s*)?([0-9]{1,3})", value, flags=re.IGNORECASE)
    if matched:
        return f"HAdV-{matched.group(1).upper()}{matched.group(2)}"
    matched = re.search(r"Human adenovirus(?:\s+type)?\s*([0-9]{1,3})", value, flags=re.IGNORECASE)
    if matched:
        number = matched.group(1)
        for row in _read_hadv_reference_manifest():
            type_label = _normalize_hadv_type_label(str(row.get("type_label") or "").strip())
            group, type_number = _hadv_type_parts(type_label)
            if group and type_number == number:
                return type_label
    return ""


def _normalize_hadv_type_label(type_label: str) -> str:
    group, number = _hadv_type_parts(type_label)
    return f"HAdV-{group}{number}" if group and number else str(type_label or "").strip()


def _hadv_type_parts(type_label: str) -> tuple[str, str]:
    matched = re.search(r"HAdV-([A-G])([0-9]{1,3})", str(type_label or ""), flags=re.IGNORECASE)
    if not matched:
        return "", ""
    return matched.group(1).upper(), matched.group(2)


def _hadv_header_matches_type_label(header: str, type_label: str) -> bool:
    group, number = _hadv_type_parts(type_label)
    if not group or not number:
        return False
    text = str(header or "")
    return bool(
        re.search(rf"HAdV-{group}\s*{number}(?!\d)", text, flags=re.IGNORECASE)
        or re.search(rf"Human adenovirus(?: type)?\s*{number}(?!\d)", text, flags=re.IGNORECASE)
        or re.search(rf"Human mastadenovirus\s+{group}\b.*?[/_\-\s]{number}(?:\[|/|_|-|\b)", text, flags=re.IGNORECASE)
    )


def _hadv_type_sort_key(type_label: str) -> tuple[str, int]:
    group, number = _hadv_type_parts(type_label)
    try:
        return group, int(number)
    except ValueError:
        return group, 999999


def _parse_hadv_gene_type(subject_id: str) -> str:
    return _extract_hadv_type_label(subject_id)


def _extract_hadv_type_number(value: str) -> str:
    group, number = _hadv_type_parts(value)
    if number:
        return number
    matched = re.search(r"([0-9]{1,3})", str(value or ""))
    return matched.group(1) if matched else ""


def _read_hadv_typing_rows() -> list[dict[str, str]]:
    typing_path = _resolve_hadv_typing_tsv()
    if not typing_path.is_file():
        return []
    with typing_path.open("r", encoding="utf-8-sig", errors="ignore", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def _hadv_typing_cell_matches(cell_value: str, observed_type: str) -> bool:
    observed_number = _extract_hadv_type_number(observed_type)
    if not observed_number:
        return False
    options = [item.strip() for item in str(cell_value or "").split(",") if item.strip()]
    return observed_number in options


def _read_hadv_reference_manifest() -> list[dict[str, str]]:
    manifest_path = (_resolve_hadv_db_dir() / "reference_genomes" / "hadv_representative_genomes_expanded.tsv").resolve()
    if not manifest_path.is_file():
        return []
    with manifest_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def _parse_hadv_phf_from_row(row: dict[str, str]) -> tuple[str, str, str]:
    def _from_header(header: str) -> str:
        return _extract_hadv_type_label(header)

    penton = _from_header(row.get("penton_header", ""))
    hexon = _from_header(row.get("hexon_header", ""))
    fiber = _from_header(row.get("fiber_header", ""))
    full_header = str(row.get("full_genome_header") or "")
    matched = re.search(r"\[P([0-9]+(?:/[0-9]+)?)H([0-9]+)F([0-9]+)\]", full_header, flags=re.IGNORECASE)
    group, _number = _hadv_type_parts(row.get("type_label", ""))
    if matched and group:
        penton = penton or f"HAdV-{group}{matched.group(1).split('/')[0]}"
        hexon = hexon or f"HAdV-{group}{matched.group(2)}"
        fiber = fiber or f"HAdV-{group}{matched.group(3)}"
    type_label = str(row.get("type_label") or "").strip()
    if type_label and not any([penton, hexon, fiber]):
        penton = hexon = fiber = type_label
    return penton, hexon, fiber


def _infer_hadv_type_from_phf(penton_type: str, hexon_type: str, fiber_type: str) -> str:
    query = {
        "Penton": str(penton_type or "").strip(),
        "Hexon": str(hexon_type or "").strip(),
        "Fiber": str(fiber_type or "").strip(),
    }
    if all(query.values()):
        for row in _read_hadv_typing_rows():
            if (
                _hadv_typing_cell_matches(row.get("Penton", ""), query["Penton"])
                and _hadv_typing_cell_matches(row.get("Hexon", ""), query["Hexon"])
                and _hadv_typing_cell_matches(row.get("Fiber", ""), query["Fiber"])
            ):
                return _normalize_hadv_type_label(str(row.get("Adenovirus Genotype") or "").strip())
    normalized_query = tuple(_normalize_hadv_type_label(value) for value in query.values())
    if normalized_query in HADV_MANUAL_PHF_COMBOS:
        return HADV_MANUAL_PHF_COMBOS[normalized_query]
    manifest_rows = _read_hadv_reference_manifest()
    for row in manifest_rows:
        if _parse_hadv_phf_from_row(row) == normalized_query:
            return _normalize_hadv_type_label(str(row.get("type_label") or "").strip())
    non_empty = [item for item in normalized_query if item]
    if non_empty and len(set(non_empty)) == 1:
        return _normalize_hadv_type_label(non_empty[0])
    return ""


def _read_hadv_full_genome_candidates(type_label: str) -> list[tuple[str, str, str]]:
    group, number = _hadv_type_parts(type_label)
    if not group or not number:
        return []
    fasta_path = _resolve_hadv_full_genome_fasta()
    manifest_path = _resolve_hadv_full_genome_manifest()
    manifest: dict[str, dict[str, str]] = {}
    if manifest_path.is_file():
        with manifest_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                accession = str(row.get("accession") or "").split(".")[0]
                if accession:
                    manifest[accession] = dict(row)
    candidates: list[tuple[str, str, str]] = []
    if fasta_path.is_file():
        with fasta_path.open("r", encoding="utf-8", errors="ignore") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                accession = str(record.id or "").split(".")[0]
                row = manifest.get(accession, {})
                row_group = str(row.get("species_group") or "").strip().upper()
                row_type = str(row.get("type_label") or "").strip()
                header = str(record.description or record.id)
                if (
                    (row_group == group and row_type == number)
                    or _hadv_header_matches_type_label(header, type_label)
                ):
                    candidates.append((str(record.id), header, str(record.seq)))
    if candidates:
        return candidates

    ref_fasta = _resolve_hadv_reference_fasta()
    if ref_fasta.is_file():
        with ref_fasta.open("r", encoding="utf-8", errors="ignore") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                header = str(record.description or record.id)
                if _hadv_header_matches_type_label(header, type_label) or _normalize_hadv_type_label(_extract_hadv_type_label(header)).lower() == _normalize_hadv_type_label(type_label).lower():
                    candidates.append((str(record.id), header, str(record.seq)))
    return candidates


def _write_hadv_candidate_fasta(candidates: list[tuple[str, str, str]], output_fasta: Path) -> None:
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with output_fasta.open("w", encoding="utf-8") as handle:
        for record_id, header, sequence in candidates:
            handle.write(f">{header or record_id}\n")
            for index in range(0, len(sequence), 80):
                handle.write(sequence[index:index + 80] + "\n")


def _extract_named_fasta_record(source_fasta: Path, record_name: str, output_fasta: Path) -> bool:
    if not source_fasta.is_file() or not record_name:
        return False
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with source_fasta.open("r", encoding="utf-8", errors="ignore") as in_handle:
        for record in SeqIO.parse(in_handle, "fasta"):
            identifiers = {str(record.id).strip(), str(record.name).strip(), str(record.description).strip()}
            if record_name in identifiers:
                with output_fasta.open("w", encoding="utf-8") as out_handle:
                    SeqIO.write(record, out_handle, "fasta")
                return output_fasta.is_file() and output_fasta.stat().st_size > 0
    return False


def _rsv_species_label(rsv_type: str) -> str:
    return f"Respiratory syncytial virus {str(rsv_type or '').strip().upper()}".strip()


def _denv_species_label(denv_type: str) -> str:
    normalized = str(denv_type or "").strip()
    return f"Dengue virus {normalized}".strip() if normalized in {"1", "2", "3", "4"} else "Dengue virus"


def _hpiv_species_label(hpiv_type: str) -> str:
    normalized = str(hpiv_type or "").strip().upper()
    return f"Human parainfluenza virus {normalized}".strip() if normalized in {"1", "2", "3", "4A", "4B"} else "Human parainfluenza virus"


def _norovirus_species_label(dual_type: str = "") -> str:
    normalized = str(dual_type or "").strip()
    return f"Norovirus {normalized}" if normalized else "Norovirus"


def _normalize_norovirus_subtype(value: str) -> str:
    return re.sub(r"\s+", " ", str(value or "").strip())


def _norovirus_subtype_match_key(value: str) -> str:
    normalized = _normalize_norovirus_subtype(value).upper()
    return re.sub(r"[\s_\-]+", "", normalized)


def _extract_norovirus_genogroup(type_label: str) -> str:
    matched = re.match(r"^(G[IVX]+)", _normalize_norovirus_subtype(type_label), flags=re.IGNORECASE)
    return matched.group(1).upper() if matched else ""


def _read_norovirus_manifest_rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open("r", encoding="utf-8-sig", errors="ignore", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t") if row]


def _read_norovirus_gene_manifest(gene_name: str) -> list[dict[str, str]]:
    return _read_norovirus_manifest_rows(_resolve_norovirus_cdc_ref_manifest(gene_name))


def _read_norovirus_full_genome_manifest_rows() -> list[dict[str, str]]:
    return _read_norovirus_manifest_rows(_resolve_norovirus_full_genome_manifest())


def _build_norovirus_subject_meta(gene_name: str) -> dict[str, dict[str, str]]:
    mapping: dict[str, dict[str, str]] = {}
    for row in _read_norovirus_gene_manifest(gene_name):
        accession = str(row.get("accession") or "").split(".", 1)[0]
        subtype = _normalize_norovirus_subtype(str(row.get("subtype") or ""))
        if not accession or not subtype:
            continue
        subject_id = f"{accession}_{subtype}_{str(gene_name or '').upper()}"
        enriched = dict(row)
        enriched["accession_core"] = accession
        enriched["subtype"] = subtype
        mapping[subject_id] = enriched
        mapping[accession] = enriched
    return mapping


def _parse_norovirus_subject_id(subject_id: str, gene_name: str = "") -> dict[str, str]:
    subject_meta = _build_norovirus_subject_meta(gene_name) if gene_name else {}
    matched = subject_meta.get(subject_id) or subject_meta.get(str(subject_id or "").split("_", 1)[0]) or {}
    if matched:
        return matched
    value = str(subject_id or "")
    suffix = f"_{str(gene_name or '').upper()}" if gene_name else ""
    core = value[:-len(suffix)] if suffix and value.endswith(suffix) else value
    accession = core.split("_", 1)[0]
    subtype = core[len(accession) + 1:] if accession and "_" in core else ""
    return {
        "accession": accession,
        "accession_core": accession.split(".", 1)[0],
        "subtype": _normalize_norovirus_subtype(subtype),
        "label": value,
    }


def _score_norovirus_type_rows(rows: list[dict[str, object]]) -> dict[str, str]:
    if not rows:
        return {"type": "", "subject": "", "coverage": "", "mean_depth": "", "covered_bases": "", "num_reads": "", "reference_count": "0"}
    by_type: dict[str, dict[str, object]] = {}
    for row in rows:
        subject = str(row.get("reference_name") or "").strip()
        meta = row.get("_meta", {})
        subtype = _normalize_norovirus_subtype(str(meta.get("subtype") or ""))
        if not subtype:
            continue
        bucket = by_type.setdefault(
            subtype,
            {
                "type": subtype,
                "subject": subject,
                "coverage": float(row.get("coverage") or 0.0),
                "mean_depth": float(row.get("mean_depth") or 0.0),
                "covered_bases": float(row.get("covered_bases") or 0.0),
                "num_reads": float(row.get("num_reads") or 0.0),
                "reference_count": set(),
            },
        )
        score = (
            float(row.get("coverage") or 0.0),
            float(row.get("mean_depth") or 0.0),
            float(row.get("covered_bases") or 0.0),
            float(row.get("num_reads") or 0.0),
        )
        best_score = (
            float(bucket["coverage"]),
            float(bucket["mean_depth"]),
            float(bucket["covered_bases"]),
            float(bucket["num_reads"]),
        )
        if score > best_score:
            bucket["subject"] = subject
            bucket["coverage"] = score[0]
            bucket["mean_depth"] = score[1]
            bucket["covered_bases"] = score[2]
            bucket["num_reads"] = score[3]
        bucket["reference_count"].add(str(meta.get("accession_core") or subject.split("_", 1)[0]))
    if not by_type:
        return {"type": "", "subject": "", "coverage": "", "mean_depth": "", "covered_bases": "", "num_reads": "", "reference_count": "0"}
    ranked = sorted(
        by_type.values(),
        key=lambda item: (
            float(item["coverage"]),
            float(item["mean_depth"]),
            float(item["covered_bases"]),
            float(item["num_reads"]),
            len(item["reference_count"]),
        ),
        reverse=True,
    )
    best = ranked[0]
    return {
        "type": str(best["type"]),
        "subject": str(best["subject"]),
        "coverage": f"{float(best['coverage']):.2f}",
        "mean_depth": f"{float(best['mean_depth']):.2f}",
        "covered_bases": f"{float(best['covered_bases']):.0f}",
        "num_reads": f"{float(best['num_reads']):.0f}",
        "reference_count": str(len(best["reference_count"])),
    }


def _run_norovirus_gene_read_typing(
    gene_name: str,
    db_fasta: Path,
    out_dir: Path,
    single_fastq: str = "",
    fq1: str = "",
    fq2: str = "",
    long_type: str = "",
    threads: int = 4,
    logf=None,
) -> dict[str, str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    bam_path = out_dir / f"{gene_name}.screening.bam"
    coverage_path = out_dir / f"{gene_name}.coverage.tsv"
    if not bam_path.is_file():
        if fq1:
            map_cmd = " ".join(
                [
                    "minimap2",
                    "-ax",
                    "sr",
                    shlex.quote(str(db_fasta)),
                    shlex.quote(str(fq1)),
                    *([shlex.quote(str(fq2))] if fq2 else []),
                    "|",
                    "samtools",
                    "sort",
                    "-o",
                    shlex.quote(str(bam_path)),
                ]
            )
        else:
            map_cmd = " ".join(
                [
                    "minimap2",
                    "-ax",
                    _choose_minimap2_preset(long_type),
                    shlex.quote(str(db_fasta)),
                    shlex.quote(str(single_fastq)),
                    "-t",
                    str(max(1, int(threads or 1))),
                    "|",
                    "samtools",
                    "sort",
                    "-o",
                    shlex.quote(str(bam_path)),
                ]
            )
        run_command(map_cmd, logf=logf)
    if not bam_path.is_file():
        return {"gene": gene_name, "type": "", "subject": "", "identity": "", "coverage": "", "mean_depth": "", "covered_bases": "", "num_reads": "", "method": "read_coverage"}
    run_command(f"samtools index -@ {max(1, int(threads or 1))} {shlex.quote(str(bam_path))}", logf=logf)
    run_command(f"samtools coverage {shlex.quote(str(bam_path))} > {shlex.quote(str(coverage_path))}", logf=logf)
    subject_meta = _build_norovirus_subject_meta(gene_name)
    rows = _parse_samtools_coverage_rows(coverage_path)
    enriched_rows = []
    for row in rows:
        subject = str(row.get("reference_name") or "").strip()
        meta = subject_meta.get(subject) or subject_meta.get(subject.split("_", 1)[0]) or {}
        if not meta:
            meta = _parse_norovirus_subject_id(subject, gene_name)
        enriched = dict(row)
        enriched["_meta"] = meta
        enriched_rows.append(enriched)
    best = _score_norovirus_type_rows(enriched_rows)
    return {
        "gene": gene_name,
        "type": best["type"],
        "subject": best["subject"],
        "identity": "",
        "coverage": best["coverage"],
        "mean_depth": best["mean_depth"],
        "covered_bases": best["covered_bases"],
        "num_reads": best["num_reads"],
        "reference_count": best["reference_count"],
        "method": "read_coverage",
    }


def _run_norovirus_gene_blast_typing(query_fasta: Path, gene_name: str, db_fasta: Path, out_dir: Path, logf=None) -> dict[str, str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    out_tsv = out_dir / f"{gene_name}.blastn.tsv"
    blastn_bin = _resolve_hadv_blastn_bin()
    cmd = " ".join(
        [
            shlex.quote(blastn_bin),
            "-task blastn",
            "-query",
            shlex.quote(str(query_fasta)),
            "-subject",
            shlex.quote(str(db_fasta)),
            "-outfmt",
            shlex.quote("6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"),
            "-evalue 1e-20",
            "-max_target_seqs 50",
            "-dust no",
            "-num_threads 4",
            "-out",
            shlex.quote(str(out_tsv)),
        ]
    )
    run_command(cmd, logf=logf)
    if not out_tsv.is_file() or out_tsv.stat().st_size == 0:
        return {"gene": gene_name, "type": "", "subject": "", "identity": "", "coverage": "", "mean_depth": "", "covered_bases": "", "num_reads": "", "reference_count": "0", "method": "blastn"}

    subject_meta = _build_norovirus_subject_meta(gene_name)
    ref_lengths: dict[str, int] = {}
    with db_fasta.open("r", encoding="utf-8", errors="ignore") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            ref_lengths[str(record.id)] = len(str(record.seq))

    by_type: dict[str, dict[str, object]] = {}
    with out_tsv.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            _qseqid, sseqid, pident, length, _mismatch, _gapopen, qstart, qend, sstart, send, evalue, bitscore = parts[:12]
            meta = subject_meta.get(sseqid) or subject_meta.get(sseqid.split("_", 1)[0]) or _parse_norovirus_subject_id(sseqid, gene_name)
            subtype = _normalize_norovirus_subtype(str(meta.get("subtype") or ""))
            if not subtype:
                continue
            ref_len = ref_lengths.get(sseqid, 0)
            try:
                qcov_ref = (float(length) / ref_len) * 100 if ref_len else 0.0
                score = (float(bitscore), float(pident), qcov_ref, float(length))
            except ValueError:
                continue
            bucket = by_type.setdefault(
                subtype,
                {
                    "type": subtype,
                    "subject": sseqid,
                    "identity": str(pident),
                    "coverage": f"{qcov_ref:.2f}",
                    "mean_depth": "",
                    "covered_bases": str(length),
                    "num_reads": "",
                    "reference_count": set(),
                    "_score": score,
                },
            )
            if score > bucket["_score"]:
                bucket["subject"] = sseqid
                bucket["identity"] = str(pident)
                bucket["coverage"] = f"{qcov_ref:.2f}"
                bucket["covered_bases"] = str(length)
                bucket["_score"] = score
            bucket["reference_count"].add(str(meta.get("accession_core") or sseqid.split("_", 1)[0]))
    if not by_type:
        return {"gene": gene_name, "type": "", "subject": "", "identity": "", "coverage": "", "mean_depth": "", "covered_bases": "", "num_reads": "", "reference_count": "0", "method": "blastn"}
    ranked = sorted(
        by_type.values(),
        key=lambda item: (
            item["_score"][0],
            item["_score"][1],
            item["_score"][2],
            len(item["reference_count"]),
        ),
        reverse=True,
    )
    best = ranked[0]
    return {
        "gene": gene_name,
        "type": str(best["type"]),
        "subject": str(best["subject"]),
        "identity": str(best["identity"]),
        "coverage": str(best["coverage"]),
        "mean_depth": "",
        "covered_bases": str(best["covered_bases"]),
        "num_reads": "",
        "reference_count": str(len(best["reference_count"])),
        "method": "blastn",
    }


def _build_norovirus_dual_type(rdrp_type: str, vp1_type: str) -> str:
    rdrp = _normalize_norovirus_subtype(rdrp_type)
    vp1 = _normalize_norovirus_subtype(vp1_type)
    if rdrp and vp1:
        return f"{rdrp}_{vp1}"
    return rdrp or vp1


def _run_norovirus_dual_typing(
    pre: str,
    query_fasta: Path | None = None,
    single_fastq: str = "",
    fq1: str = "",
    fq2: str = "",
    long_type: str = "",
    threads: int = 4,
    output_dir: Path | None = None,
    logf=None,
) -> dict[str, object]:
    typing_dir = (output_dir or (Path(f"{pre}_norovirus_reference_selection") / "typing")).resolve()
    typing_dir.mkdir(parents=True, exist_ok=True)
    gene_dbs = {
        "rdrp": _resolve_norovirus_cdc_ref_fasta("rdrp"),
        "vp1": _resolve_norovirus_cdc_ref_fasta("vp1"),
    }
    hits: dict[str, dict[str, str]] = {}
    use_blast = query_fasta is not None and query_fasta.is_file() and query_fasta.stat().st_size > 0
    for gene_name, db_fasta in gene_dbs.items():
        if not db_fasta.is_file():
            hits[gene_name] = {"gene": gene_name, "type": "", "subject": "", "identity": "", "coverage": "", "mean_depth": "", "covered_bases": "", "num_reads": "", "reference_count": "0", "method": "missing_db"}
            continue
        if use_blast:
            hits[gene_name] = _run_norovirus_gene_blast_typing(query_fasta, gene_name, db_fasta, typing_dir, logf=logf)
        else:
            hits[gene_name] = _run_norovirus_gene_read_typing(
                gene_name,
                db_fasta,
                typing_dir,
                single_fastq=single_fastq,
                fq1=fq1,
                fq2=fq2,
                long_type=long_type,
                threads=threads,
                logf=logf,
            )
    status_by_gene: dict[str, str] = {}
    note_by_gene: dict[str, str] = {}
    for gene_name in ["rdrp", "vp1"]:
        hit = hits.get(gene_name, {})
        if _is_low_support_gene_typing_hit(hit):
            status_by_gene[gene_name] = "low_support"
            note_by_gene[gene_name] = _low_support_gene_typing_note(gene_name.upper(), hit)
            adjusted = dict(hit)
            adjusted["type"] = ""
            hits[gene_name] = adjusted
        else:
            status_by_gene[gene_name] = "ready"
            note_by_gene[gene_name] = ""
    rdrp_type = hits.get("rdrp", {}).get("type", "")
    vp1_type = hits.get("vp1", {}).get("type", "")
    dual_type = _build_norovirus_dual_type(rdrp_type, vp1_type)
    summary_path = typing_dir / "dual_typing.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "gene", "matched_type", "subject", "identity", "coverage", "mean_depth", "covered_bases", "num_reads", "reference_count", "method", "status", "note"])
        for gene_name in ["rdrp", "vp1"]:
            hit = hits.get(gene_name, {})
            writer.writerow([
                pre,
                gene_name,
                hit.get("type", ""),
                hit.get("subject", ""),
                hit.get("identity", ""),
                hit.get("coverage", ""),
                hit.get("mean_depth", ""),
                hit.get("covered_bases", ""),
                hit.get("num_reads", ""),
                hit.get("reference_count", "0"),
                hit.get("method", ""),
                status_by_gene.get(gene_name, "ready"),
                note_by_gene.get(gene_name, ""),
            ])
    overall_status = "low_support" if any(value == "low_support" for value in status_by_gene.values()) else "ready"
    overall_note = "；".join([value for value in note_by_gene.values() if value])
    return {
        "rdrp_type": rdrp_type,
        "vp1_type": vp1_type,
        "dual_type": dual_type,
        "summary_path": str(summary_path.resolve()),
        "hits": hits,
        "status": overall_status,
        "note": overall_note,
    }


def _load_norovirus_fragment_typing_rows() -> list[dict[str, str]]:
    typing_path = _resolve_norovirus_fragment_typing_tsv()
    if not typing_path.is_file():
        return []
    with typing_path.open("r", encoding="utf-8-sig", errors="ignore", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t") if row]


def _load_noro_fasta_sequence_map() -> dict[str, tuple[str, str]]:
    root = _resolve_norovirus_db_dir() / "full_genomes"
    sequence_map: dict[str, tuple[str, str]] = {}
    for fasta_path in sorted(root.glob("noro*.fasta")):
        if fasta_path.name in {"noro_all_for_typing.fasta", "norovirus_human.fasta"}:
            continue
        with fasta_path.open("r", encoding="utf-8", errors="ignore") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                contig = str(record.id or "").strip()
                if contig and contig not in sequence_map:
                    sequence_map[contig] = (str(record.description or record.id), str(record.seq))
    return sequence_map


def _collect_norovirus_reference_candidates(rdrp_type: str, vp1_type: str) -> list[dict[str, str]]:
    rdrp_type = _normalize_norovirus_subtype(rdrp_type)
    vp1_type = _normalize_norovirus_subtype(vp1_type)
    rdrp_key = _norovirus_subtype_match_key(rdrp_type)
    vp1_key = _norovirus_subtype_match_key(vp1_type)
    typing_rows = _load_norovirus_fragment_typing_rows()
    sequence_map = _load_noro_fasta_sequence_map()
    exact_rows: list[dict[str, str]] = []
    partial_rows: list[dict[str, str]] = []
    for row in typing_rows:
        contig = str(row.get("contig名称") or "").strip()
        row_vp1 = _normalize_norovirus_subtype(str(row.get("vp1分型") or ""))
        row_rdrp = _normalize_norovirus_subtype(str(row.get("rdrp分型") or ""))
        row_vp1_key = _norovirus_subtype_match_key(row_vp1)
        row_rdrp_key = _norovirus_subtype_match_key(row_rdrp)
        if not contig or contig not in sequence_map:
            continue
        if row_vp1_key == vp1_key and row_rdrp_key == rdrp_key:
            exact_rows.append(row)
        elif (vp1_key and row_vp1_key == vp1_key) or (rdrp_key and row_rdrp_key == rdrp_key):
            partial_rows.append(row)
    selected_rows = exact_rows or partial_rows
    candidates_by_contig: dict[str, dict[str, str]] = {}
    for row in selected_rows:
        contig = str(row.get("contig名称") or "").strip()
        header, sequence = sequence_map.get(contig, ("", ""))
        if not sequence:
            continue
        if contig in candidates_by_contig:
            continue
        candidates_by_contig[contig] = {
            "accession": contig.split(".", 1)[0],
            "contig_name": contig,
            "header": header,
            "sequence": sequence,
            "rdrp_type": str(row.get("rdrp分型") or rdrp_type or ""),
            "vp1_type": str(row.get("vp1分型") or vp1_type or ""),
            "genogroup": _extract_norovirus_genogroup(str(row.get("vp1分型") or "")) or _extract_norovirus_genogroup(str(row.get("rdrp分型") or "")) or _extract_norovirus_genogroup(vp1_type) or _extract_norovirus_genogroup(rdrp_type),
            "genotype": str(row.get("vp1分型") or ""),
        }
    return list(candidates_by_contig.values())


def _write_norovirus_candidate_fasta(candidates: list[dict[str, str]], output_fasta: Path) -> None:
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with output_fasta.open("w", encoding="utf-8") as handle:
        for candidate in candidates:
            header = str(candidate.get("header") or candidate.get("accession") or "candidate")
            sequence = str(candidate.get("sequence") or "")
            handle.write(f">{header}\n")
            for index in range(0, len(sequence), 80):
                handle.write(sequence[index:index + 80] + "\n")


def _resolve_norovirus_vadr_model_dir(project_root: Path) -> Path | None:
    primary = (project_root / "soft" / "vadr-models-calici").resolve()
    if primary.is_dir():
        return primary
    return None


def _run_vadr_norovirus_annotation(pre: str, input_fasta: Path, output_root: Path, logf=None) -> dict[str, str]:
    project_root = _project_root()
    model_dir = _resolve_norovirus_vadr_model_dir(project_root)
    if model_dir is None:
        return {"status": "missing", "note": "未找到诺如病毒 VADR 模型目录", "output_dir": "", "gff_path": ""}
    env = _build_vadr_env(project_root, model_dir)
    if env is None:
        return {"status": "missing", "note": "VADR 运行依赖不完整", "output_dir": "", "gff_path": ""}
    norovirus_env = dict(env)
    alt_bio_easel_dir = project_root / "soft" / "Bio-Easel-ncov"
    alt_bio_lib = alt_bio_easel_dir / "blib" / "lib"
    alt_bio_arch = alt_bio_easel_dir / "blib" / "arch"
    if alt_bio_lib.is_dir() and alt_bio_arch.is_dir():
        norovirus_env["VADRBIOEASELDIR"] = str(alt_bio_easel_dir)
        norovirus_env["PERL5LIB"] = os.pathsep.join(
            [
                str(project_root / "soft" / "vadr"),
                str(project_root / "soft" / "sequip"),
                str(alt_bio_lib),
                str(alt_bio_arch),
            ]
        )
    norovirus_env["PATH"] = os.pathsep.join(
        [
            "/usr/bin",
            "/bin",
            str(norovirus_env.get("PATH") or ""),
        ]
    )
    ncov_perl = conda_env_path("ncov", "bin", "perl")
    perl_candidates = [ncov_perl, "/usr/bin/perl", _resolve_vadr_perl_bin(), "perl"]
    perl_bin = _resolve_working_perl_with_module(
        norovirus_env,
        "Bio::Easel::MSA",
        perl_candidates,
    )

    vadr_root = output_root / "vadr"
    vadr_root.mkdir(parents=True, exist_ok=True)
    output_dir = vadr_root / f"{pre}_vadr"
    prefix = f"{output_dir.name}.vadr"
    pass_tbl = output_dir / f"{prefix}.pass.tbl"
    fail_tbl = output_dir / f"{prefix}.fail.tbl"
    gff_path = vadr_root / f"{pre}.vadr.gff3"
    gff_ready = gff_path.is_file() and gff_path.stat().st_size > 0
    if output_dir.is_dir() and (pass_tbl.is_file() or fail_tbl.is_file()):
        if not gff_ready:
            annotate_tbl2gff = project_root / "soft" / "vadr" / "miniscripts" / "annotate-tbl2gff.pl"
            source_tbl = pass_tbl if pass_tbl.is_file() else fail_tbl
            try:
                run_command(
                    f"{shlex.quote(str(perl_bin))} {shlex.quote(str(annotate_tbl2gff))} {shlex.quote(str(source_tbl))} > {shlex.quote(str(gff_path))}",
                    logf=logf,
                    env=norovirus_env,
                )
            except Exception:
                pass
            gff_ready = gff_path.is_file() and gff_path.stat().st_size > 0
        return {
            "status": "ready",
            "note": "已检测到现有诺如病毒 VADR 注释结果",
            "output_dir": str(output_dir.resolve()),
            "gff_path": str(gff_path.resolve()) if gff_ready else "",
        }

    vadr_script = project_root / "soft" / "vadr" / "v-annotate.pl"
    cmd = " ".join(
        [
            shlex.quote(str(perl_bin)),
            shlex.quote(str(vadr_script)),
            "--split",
            "--cpu",
            "1",
            "--group",
            "Norovirus",
            "--nomisc",
            "--noprotid",
            "--mkey",
            "calici",
            "-r",
            "--mdir",
            shlex.quote(str(model_dir)),
            "-f",
            shlex.quote(str(input_fasta)),
            shlex.quote(str(output_dir)),
        ]
    )
    try:
        run_command(cmd, logf=logf, env=norovirus_env)
    except Exception as exc:
        return {
            "status": "failed",
            "note": f"诺如病毒 VADR 注释失败: {exc}",
            "output_dir": str(output_dir.resolve()) if output_dir.exists() else "",
            "gff_path": "",
        }

    annotate_tbl2gff = project_root / "soft" / "vadr" / "miniscripts" / "annotate-tbl2gff.pl"
    source_tbl = pass_tbl if pass_tbl.is_file() else fail_tbl
    if source_tbl.is_file():
        try:
            run_command(
                f"{shlex.quote(str(perl_bin))} {shlex.quote(str(annotate_tbl2gff))} {shlex.quote(str(source_tbl))} > {shlex.quote(str(gff_path))}",
                logf=logf,
                env=norovirus_env,
            )
        except Exception:
            pass
    gff_ready = gff_path.is_file() and gff_path.stat().st_size > 0
    return {
        "status": "ready",
        "note": "诺如病毒 VADR 注释完成",
        "output_dir": str(output_dir.resolve()),
        "gff_path": str(gff_path.resolve()) if gff_ready else "",
    }


def prepare_norovirus_reference_annotation(pre: str, reference_fasta: Path, output_root: Path, logf=None) -> Path | None:
    result = _run_vadr_norovirus_annotation(pre, reference_fasta, output_root, logf=logf)
    gff_path = str(result.get("gff_path") or "").strip()
    if result.get("status") == "ready" and gff_path:
        gff_file = Path(gff_path)
        if gff_file.is_file() and gff_file.stat().st_size > 0:
            return gff_file
    return None


def _parse_norovirus_gff_attributes(raw_text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for item in str(raw_text or "").strip().split(";"):
        token = item.strip()
        if not token:
            continue
        if "=" in token:
            key, value = token.split("=", 1)
        elif ":" in token:
            key, value = token.split(":", 1)
        else:
            continue
        attributes[key.strip()] = value.strip()
    return attributes


def _load_norovirus_typing_call(pre: str) -> dict[str, str]:
    selection_path = Path(f"{pre}_norovirus_reference_selection") / "selection.tsv"
    typing_path = Path(f"{pre}_norovirus_reference_selection") / "typing" / "dual_typing.tsv"
    result = {"dual_type": "", "rdrp_type": "", "vp1_type": ""}
    if selection_path.is_file():
        try:
            with selection_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"), None) or {}
            result["dual_type"] = str(row.get("dual_type") or "").strip()
            result["rdrp_type"] = str(row.get("rdrp_type") or "").strip()
            result["vp1_type"] = str(row.get("vp1_type") or "").strip()
        except OSError:
            pass
    if typing_path.is_file():
        try:
            with typing_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    gene_name = str(row.get("gene") or "").strip().lower()
                    matched_type = str(row.get("matched_type") or "").strip()
                    if gene_name in {"rdrp", "vp1"} and matched_type:
                        result[f"{gene_name}_type"] = matched_type
        except OSError:
            pass
    if not result["dual_type"]:
        result["dual_type"] = _build_norovirus_dual_type(result["rdrp_type"], result["vp1_type"])
    return result


def _find_norovirus_gene_feature(gff_path: Path, gene_name: str) -> dict[str, object] | None:
    target = str(gene_name or "").strip().lower()
    if not gff_path.is_file():
        return None
    candidates: list[tuple[int, dict[str, object]]] = []
    with gff_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, _source, feature_type, start_text, end_text, _score, strand, _phase, attrs_text = parts[:9]
            attrs = _parse_norovirus_gff_attributes(attrs_text)
            label = " ".join(
                [
                    str(attrs.get("gene") or ""),
                    str(attrs.get("product") or ""),
                    str(attrs.get("Name") or ""),
                ]
            ).upper()
            is_match = False
            priority = 9
            if target == "vp1":
                if "VP1" in label:
                    is_match = True
                    priority = 0 if feature_type == "CDS" else 1
                elif str(attrs.get("gene") or "").strip().upper() == "ORF2":
                    is_match = True
                    priority = 2
            elif target == "rdrp":
                if "RDRP" in label or "POLYMERASE" in label:
                    is_match = True
                    priority = 0 if feature_type == "mat_peptide" else 1
            if not is_match:
                continue
            try:
                start = int(start_text)
                end = int(end_text)
            except ValueError:
                continue
            candidates.append(
                (
                    priority,
                    {
                        "seqid": seqid,
                        "feature_type": feature_type,
                        "start": start,
                        "end": end,
                        "strand": strand,
                        "label": attrs.get("product") or attrs.get("gene") or target.upper(),
                    },
                )
            )
    if not candidates:
        return None
    candidates.sort(key=lambda item: (item[0], item[1]["start"], item[1]["end"]))
    return candidates[0][1]


def _extract_norovirus_gene_region(sample_fasta: Path, gff_path: Path, gene_name: str, sample_label: str) -> tuple[SeqRecord | None, dict[str, object]]:
    feature = _find_norovirus_gene_feature(gff_path, gene_name)
    if feature is None:
        return None, {"status": "missing", "note": f"GFF 中未找到 {gene_name.upper()} 区域"}
    records = list(SeqIO.parse(str(sample_fasta), "fasta"))
    if not records:
        return None, {"status": "missing", "note": f"{sample_fasta.name} 中未找到序列记录"}
    seqid = str(feature.get("seqid") or "").strip()
    record = next((item for item in records if str(item.id).strip() == seqid), None)
    if record is None:
        record = max(records, key=lambda item: len(str(item.seq)))
    start = int(feature["start"])
    end = int(feature["end"])
    sequence = str(record.seq)
    if start < 1 or end > len(sequence) or start > end:
        return None, {"status": "missing", "note": f"{gene_name.upper()} 坐标超出样本序列长度"}
    fragment = sequence[start - 1:end]
    if str(feature.get("strand") or "+") == "-":
        fragment = str(Seq(fragment).reverse_complement())
    fragment = fragment.strip().strip("Nn")
    if not fragment:
        return None, {"status": "missing", "note": f"{gene_name.upper()} 区域仅包含空白或 N"}
    record_id = f"{sample_label}_{gene_name.lower()}_sample"
    return SeqRecord(Seq(fragment), id=record_id, description=f"{sample_label} {gene_name.upper()} sample"), {
        "status": "ready",
        "start": start,
        "end": end,
        "strand": str(feature.get("strand") or "+"),
        "seqid": str(record.id),
        "length": len(fragment),
    }


def _filter_norovirus_backup_refs_by_genogroup(gene_name: str, genogroup: str) -> tuple[list[dict[str, str]], list[SeqRecord]]:
    manifest_path = _resolve_norovirus_backup_ref_manifest(gene_name)
    fasta_path = _resolve_norovirus_backup_ref_fasta(gene_name)
    if not manifest_path.is_file() or not fasta_path.is_file():
        return [], []
    genogroup_text = str(genogroup or "").strip().upper()
    with manifest_path.open("r", encoding="utf-8-sig", errors="ignore", newline="") as handle:
        manifest_rows = [dict(row) for row in csv.DictReader(handle, delimiter="\t") if row]
    selected_rows = [
        row for row in manifest_rows
        if _extract_norovirus_genogroup(str(row.get("subtype") or "")) == genogroup_text
    ]
    if not selected_rows:
        return [], []
    records_by_accession: dict[str, SeqRecord] = {}
    for record in SeqIO.parse(str(fasta_path), "fasta"):
        accession = str(record.id).split("_", 1)[0].strip()
        if accession and accession not in records_by_accession:
            records_by_accession[accession] = record
    ordered_rows: list[dict[str, str]] = []
    ordered_records: list[SeqRecord] = []
    for row in selected_rows:
        accession = str(row.get("accession") or "").strip()
        template = records_by_accession.get(accession)
        if template is None:
            continue
        subtype = str(row.get("subtype") or "").strip()
        label_gene = str(row.get("gene") or gene_name).strip()
        canonical_id = f"{accession}_{subtype}_{label_gene}"
        canonical_label = _sanitize_tree_label(canonical_id)
        cloned = SeqRecord(
            template.seq,
            id=canonical_label,
            name=canonical_label,
            description=f"{canonical_label} {template.description}".strip(),
        )
        ordered_rows.append(row)
        ordered_records.append(cloned)
    return ordered_rows, ordered_records


def _sanitize_tree_label(value: str) -> str:
    sanitized = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value or "").strip())
    return sanitized.strip("_") or "sample"


def _build_norovirus_gene_tree(
    gene_name: str,
    sample_record: SeqRecord,
    backup_rows: list[dict[str, str]],
    backup_records: list[SeqRecord],
    out_dir: Path,
    logf=None,
) -> dict[str, str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    sample_fasta = out_dir / f"{gene_name}.sample.fasta"
    backup_fasta = out_dir / f"{gene_name}.backup_refs.fasta"
    combined_fasta = out_dir / f"{gene_name}.combined.fasta"
    aligned_fasta = out_dir / f"{gene_name}.aligned.fasta"
    tree_path = out_dir / f"{gene_name}.tree.nwk"
    members_tsv = out_dir / f"{gene_name}.members.tsv"

    SeqIO.write([sample_record], str(sample_fasta), "fasta")
    SeqIO.write(backup_records, str(backup_fasta), "fasta")
    SeqIO.write([sample_record, *backup_records], str(combined_fasta), "fasta")
    with members_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["tree_label", "member_type", "gene", "subtype", "accession", "backup_rank"])
        writer.writerow([sample_record.id, "sample", gene_name.upper(), "", "", ""])
        for row, record in zip(backup_rows, backup_records):
            writer.writerow([
                record.id,
                "backup_ref",
                str(row.get("gene") or gene_name).strip(),
                str(row.get("subtype") or "").strip(),
                str(row.get("accession") or "").strip(),
                str(row.get("backup_rank") or "").strip(),
            ])

    if tree_path.is_file() and tree_path.stat().st_size > 0:
        return {
            "status": "ready",
            "tree_path": str(tree_path.resolve()),
            "aligned_fasta": str(aligned_fasta.resolve()),
            "combined_fasta": str(combined_fasta.resolve()),
            "members_tsv": str(members_tsv.resolve()),
            "sample_fasta": str(sample_fasta.resolve()),
            "backup_fasta": str(backup_fasta.resolve()),
            "member_count": str(len(backup_records) + 1),
            "backup_count": str(len(backup_records)),
        }
    if not aligned_fasta.is_file() or aligned_fasta.stat().st_size == 0:
        run_command(
            f"{shlex.quote('mafft')} --retree 1 --maxiterate 0 --quiet {shlex.quote(str(combined_fasta))} > {shlex.quote(str(aligned_fasta))}",
            logf=logf,
        )
    alignment = AlignIO.read(str(aligned_fasta), "fasta")
    if len(alignment) < 2:
        return {"status": "missing", "note": f"{gene_name.upper()} 对齐序列不足，无法建树", "tree_path": ""}
    if len(alignment) == 2:
        labels = [record.id for record in alignment]
        tree_text = f"({labels[0]}:0.1,{labels[1]}:0.1);"
        tree_path.write_text(tree_text, encoding="utf-8")
    else:
        calculator = DistanceCalculator("identity")
        distance_matrix = calculator.get_distance(alignment)
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(distance_matrix)
        Phylo.write(tree, str(tree_path), "newick")
    return {
        "status": "ready",
        "tree_path": str(tree_path.resolve()),
        "aligned_fasta": str(aligned_fasta.resolve()),
        "combined_fasta": str(combined_fasta.resolve()),
        "members_tsv": str(members_tsv.resolve()),
        "sample_fasta": str(sample_fasta.resolve()),
        "backup_fasta": str(backup_fasta.resolve()),
        "member_count": str(len(backup_records) + 1),
        "backup_count": str(len(backup_records)),
    }


def build_norovirus_gene_phylogeny_assets(pre: str, sample_fasta: Path, gff_path: Path, logf=None) -> dict[str, object]:
    output_root = Path(f"{pre}_norovirus_reference_selection") / "phylogeny"
    output_root.mkdir(parents=True, exist_ok=True)
    summary_path = output_root / "summary.tsv"
    typing_call = _load_norovirus_typing_call(pre)
    sample_label = _sanitize_tree_label(pre)
    results: list[dict[str, str]] = []
    for gene_name, subtype_key in [("vp1", "vp1_type"), ("rdrp", "rdrp_type")]:
        subtype = str(typing_call.get(subtype_key) or "").strip()
        genogroup = _extract_norovirus_genogroup(subtype)
        result_row = {
            "gene": gene_name.upper(),
            "subtype": subtype,
            "genogroup": genogroup,
            "status": "missing",
            "note": "",
            "sample_region_path": "",
            "backup_fasta": "",
            "tree_path": "",
            "member_count": "0",
        }
        if not sample_fasta.is_file():
            result_row["note"] = f"未找到样本序列文件 {sample_fasta.name}"
            results.append(result_row)
            continue
        if not gff_path.is_file():
            result_row["note"] = f"未找到 GFF 注释文件 {gff_path.name}"
            results.append(result_row)
            continue
        if not genogroup:
            result_row["note"] = f"未获得 {gene_name.upper()} 分型对应的大亚型"
            results.append(result_row)
            continue
        sample_record, sample_meta = _extract_norovirus_gene_region(sample_fasta, gff_path, gene_name, sample_label)
        if sample_record is None:
            result_row["note"] = str(sample_meta.get("note") or f"未提取到 {gene_name.upper()} 区域")
            results.append(result_row)
            continue
        backup_rows, backup_records = _filter_norovirus_backup_refs_by_genogroup(gene_name, genogroup)
        if not backup_records:
            result_row["sample_region_path"] = ""
            result_row["note"] = f"备用库中未找到 {gene_name.upper()} {genogroup} 参考序列"
            results.append(result_row)
            continue
        gene_dir = output_root / gene_name
        try:
            tree_result = _build_norovirus_gene_tree(gene_name, sample_record, backup_rows, backup_records, gene_dir, logf=logf)
        except Exception as exc:
            result_row["note"] = f"{gene_name.upper()} 建树失败: {exc}"
            results.append(result_row)
            continue
        result_row["status"] = str(tree_result.get("status") or "missing")
        result_row["note"] = str(tree_result.get("note") or "")
        result_row["sample_region_path"] = str((gene_dir / f"{gene_name}.sample.fasta").resolve())
        result_row["backup_fasta"] = str(tree_result.get("backup_fasta") or "")
        result_row["tree_path"] = str(tree_result.get("tree_path") or "")
        result_row["member_count"] = str(tree_result.get("member_count") or "0")
        results.append(result_row)
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["gene", "subtype", "genogroup", "status", "note", "sample_region_path", "backup_fasta", "tree_path", "member_count"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(results)
    return {
        "status": "ready" if any(row.get("status") == "ready" for row in results) else "missing",
        "summary_path": str(summary_path.resolve()),
        "rows": results,
    }


def resolve_norovirus_reference(
    pre: str,
    species: str = "",
    requested_ref: str = "",
    single_fastq: str = "",
    fq1: str = "",
    fq2: str = "",
    long_type: str = "",
    threads: int = 4,
    query_fasta: str = "",
    logf=None,
) -> dict[str, str]:
    _clear_sample_skip_flag(pre)
    requested = str(requested_ref or "").strip()
    if requested and Path(requested).is_file():
        return {
            "status": "ready",
            "dual_type": _infer_norovirus_dual_type_from_label(species),
            "rdrp_type": "",
            "vp1_type": "",
            "species_label": species or "Norovirus",
            "reference_path": str(Path(requested).resolve()),
            "gff_path": "nogtf",
            "summary_path": "",
            "typing_summary_path": "",
        }
    screening_dir = Path(f"{pre}_norovirus_reference_selection")
    screening_dir.mkdir(exist_ok=True)
    query_path = Path(query_fasta) if str(query_fasta or "").strip() else None
    typing_result = _run_norovirus_dual_typing(
        pre,
        query_fasta=query_path,
        single_fastq=single_fastq,
        fq1=fq1,
        fq2=fq2,
        long_type=long_type,
        threads=threads,
        output_dir=screening_dir / "typing",
        logf=logf,
    )
    rdrp_type = str(typing_result.get("rdrp_type") or "").strip()
    vp1_type = str(typing_result.get("vp1_type") or "").strip()
    dual_type = str(typing_result.get("dual_type") or "").strip()
    if not dual_type:
        dual_type = _infer_norovirus_dual_type_from_label(species)

    summary_path = screening_dir / "selection.tsv"
    if not dual_type:
        status_text = "skipped" if str(typing_result.get("status") or "").strip() == "low_support" else "missing"
        note_text = str(typing_result.get("note") or "").strip() or f"未能根据 RdRp/VP1 或 species={species or '-'} 确定诺如病毒分型"
        with summary_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sample", "dual_type", "rdrp_type", "vp1_type", "status", "note", "typing_summary_path"])
            writer.writerow([pre, "", rdrp_type, vp1_type, status_text, note_text, str(typing_result.get("summary_path") or "")])
        if status_text == "skipped":
            _write_sample_skip_flag(pre, note_text)
        return {
            "status": status_text,
            "dual_type": "",
            "rdrp_type": rdrp_type,
            "vp1_type": vp1_type,
            "species_label": species or "Norovirus",
            "reference_path": "",
            "gff_path": "nogtf",
            "summary_path": str(summary_path.resolve()),
            "typing_summary_path": str(typing_result.get("summary_path") or ""),
        }

    candidates = _collect_norovirus_reference_candidates(rdrp_type, vp1_type)
    candidate_fasta = screening_dir / "candidate_references.fasta"
    if candidates:
        _write_norovirus_candidate_fasta(candidates, candidate_fasta)
    else:
        with summary_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sample", "dual_type", "rdrp_type", "vp1_type", "status", "note", "typing_summary_path"])
            writer.writerow([pre, dual_type, rdrp_type, vp1_type, "missing", "未找到对应 Norovirus 分型候选参考库", str(typing_result.get("summary_path") or "")])
        return {
            "status": "missing",
            "dual_type": dual_type,
            "rdrp_type": rdrp_type,
            "vp1_type": vp1_type,
            "species_label": _norovirus_species_label(dual_type),
            "reference_path": "",
            "gff_path": "nogtf",
            "summary_path": str(summary_path.resolve()),
            "typing_summary_path": str(typing_result.get("summary_path") or ""),
        }

    bam_path = screening_dir / "candidate.screening.bam"
    coverage_path = screening_dir / "candidate.coverage.tsv"
    if not bam_path.is_file():
        if fq1:
            map_cmd = " ".join(
                [
                    "minimap2",
                    "-ax",
                    "sr",
                    shlex.quote(str(candidate_fasta)),
                    shlex.quote(str(fq1)),
                    *([shlex.quote(str(fq2))] if fq2 else []),
                    "|",
                    "samtools",
                    "sort",
                    "-o",
                    shlex.quote(str(bam_path)),
                ]
            )
        else:
            map_cmd = " ".join(
                [
                    "minimap2",
                    "-ax",
                    _choose_minimap2_preset(long_type),
                    shlex.quote(str(candidate_fasta)),
                    shlex.quote(str(single_fastq)),
                    "-t",
                    str(max(1, int(threads or 1))),
                    "|",
                    "samtools",
                    "sort",
                    "-o",
                    shlex.quote(str(bam_path)),
                ]
            )
        run_command(map_cmd, logf=logf)
    if bam_path.is_file():
        run_command(f"samtools index -@ {max(1, int(threads or 1))} {shlex.quote(str(bam_path))}", logf=logf)
        run_command(f"samtools coverage {shlex.quote(str(bam_path))} > {shlex.quote(str(coverage_path))}", logf=logf)
    rows = _parse_samtools_coverage_rows(coverage_path)
    rows.sort(
        key=lambda item: (
            float(item.get("coverage") or 0.0),
            float(item.get("mean_depth") or 0.0),
            float(item.get("covered_bases") or 0.0),
            float(item.get("num_reads") or 0.0),
        ),
        reverse=True,
    )
    if not rows:
        return {
            "status": "missing",
            "dual_type": dual_type,
            "rdrp_type": rdrp_type,
            "vp1_type": vp1_type,
            "species_label": _norovirus_species_label(dual_type),
            "reference_path": "",
            "gff_path": "nogtf",
            "summary_path": str(summary_path.resolve()),
            "typing_summary_path": str(typing_result.get("summary_path") or ""),
        }
    best = rows[0]
    best_reference_name = str(best.get("reference_name") or "").strip()
    best_reference_fasta = screening_dir / "norovirus.best_reference.fasta"
    if not _extract_named_fasta_record(candidate_fasta, best_reference_name, best_reference_fasta):
        shutil.copy2(candidate_fasta, best_reference_fasta)
    accession = best_reference_name.split()[0].split(".")[0]
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "dual_type", "rdrp_type", "vp1_type", "coverage", "mean_depth", "covered_bases", "num_reads", "reference_name", "reference_accession", "reference_path", "gff_path", "typing_summary_path"])
        writer.writerow([
            pre,
            dual_type,
            rdrp_type,
            vp1_type,
            f"{float(best.get('coverage') or 0.0):.6f}",
            f"{float(best.get('mean_depth') or 0.0):.6f}",
            f"{float(best.get('covered_bases') or 0.0):.0f}",
            f"{float(best.get('num_reads') or 0.0):.0f}",
            best_reference_name,
            accession,
            str(best_reference_fasta.resolve()),
            "nogtf",
            str(typing_result.get("summary_path") or ""),
        ])
    return {
        "status": "ready",
        "dual_type": dual_type,
        "rdrp_type": rdrp_type,
        "vp1_type": vp1_type,
        "species_label": _norovirus_species_label(dual_type),
        "reference_path": str(best_reference_fasta.resolve()),
        "gff_path": "nogtf",
        "summary_path": str(summary_path.resolve()),
        "typing_summary_path": str(typing_result.get("summary_path") or ""),
    }


def run_norovirus_consensus_typing(pre: str, consensus_fasta: Path, logf=None) -> dict[str, str]:
    if not consensus_fasta.is_file() or consensus_fasta.stat().st_size == 0:
        return {"status": "missing", "summary_path": "", "dual_type": "", "rdrp_type": "", "vp1_type": ""}
    output_dir = Path(f"{pre}_norovirus_reference_selection") / "consensus_typing"
    result = _run_norovirus_dual_typing(pre, query_fasta=consensus_fasta, output_dir=output_dir, logf=logf)
    summary_path = output_dir / "consensus_typing.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "dual_type", "rdrp_type", "vp1_type", "typing_summary_path"])
        writer.writerow([pre, result.get("dual_type", ""), result.get("rdrp_type", ""), result.get("vp1_type", ""), result.get("summary_path", "")])
    return {
        "status": "ready",
        "summary_path": str(summary_path.resolve()),
        "dual_type": str(result.get("dual_type") or ""),
        "rdrp_type": str(result.get("rdrp_type") or ""),
        "vp1_type": str(result.get("vp1_type") or ""),
    }


def _resolve_hepatovirus_db_dir() -> Path:
    env_root = str(os.environ.get("META_HEPATOVIRUS_DB_DIR") or "").strip()
    if env_root:
        return Path(env_root).expanduser().resolve()
    database_root = str(os.environ.get("META_DATABASE_ROOT") or "").strip()
    project_root = _project_root()
    candidates: list[Path] = []
    if database_root:
        candidates.append((Path(database_root) / "virus" / "Hepatovirus").expanduser())
    candidates.extend(
        [
            project_root / "database" / "virus" / "Hepatovirus",
            Path("/data/deploy/meta_genome/database/virus/Hepatovirus"),
        ]
    )
    return next((candidate.resolve() for candidate in candidates if candidate.exists()), candidates[0].resolve())


def _resolve_hepatovirus_reference_manifest() -> Path:
    db_dir = _resolve_hepatovirus_db_dir()
    combined_path = (db_dir / "broad_reference_genomes" / "hepatitis_broad_reference_genomes_manifest.tsv").resolve()
    if combined_path.exists():
        return combined_path
    return (db_dir / "reference_genomes" / "hepatovirus_typing_reference_genomes_manifest.tsv").resolve()


def _resolve_hav_subtype_manifest() -> Path:
    return (_resolve_hepatovirus_db_dir() / "HAV_subtypes" / "hav_subtype_complete_genomes_manifest.tsv").resolve()


def _resolve_hepatitis_type_reference_manifest(broad_type: str) -> Path:
    db_dir = _resolve_hepatovirus_db_dir()
    normalized = str(broad_type or "").strip().upper()
    if normalized == "HAV":
        return _resolve_hav_subtype_manifest()
    manifest_names = {
        "HBV": ("typingB_reference_genomes", "typingB_reference_genomes_manifest.tsv"),
        "HCV": ("typingC_reference_genomes", "typingC_reference_genomes_manifest.tsv"),
        "HDV": ("typingD_reference_genomes", "typingD_reference_genomes_manifest.tsv"),
        "HEV": ("typingE_reference_genomes", "typingE_reference_genomes_manifest.tsv"),
    }
    folder, filename = manifest_names.get(normalized, ("", ""))
    return (db_dir / folder / filename).resolve() if folder else Path("")


def _load_hepatovirus_reference_manifest() -> list[dict[str, str]]:
    return _load_rhinovirus_manifest(_resolve_hepatovirus_reference_manifest())


def _load_hav_subtype_manifest() -> list[dict[str, str]]:
    return _load_rhinovirus_manifest(_resolve_hav_subtype_manifest())


def _load_hepatitis_type_reference_manifest(broad_type: str) -> list[dict[str, str]]:
    manifest_path = _resolve_hepatitis_type_reference_manifest(broad_type)
    if not str(manifest_path) or not manifest_path.is_file():
        return []
    return _load_rhinovirus_manifest(manifest_path)


def _read_first_fasta_sequence(fasta_path: Path) -> tuple[str, str]:
    if not fasta_path.is_file() or fasta_path.stat().st_size == 0:
        return "", ""
    record = next(SeqIO.parse(str(fasta_path), "fasta"), None)
    if record is None:
        return "", ""
    return str(record.description or record.id).strip(), str(record.seq).upper()


def _hepatovirus_species_label(broad_type: str = "", fallback: str = "") -> str:
    normalized = str(broad_type or "").strip().upper()
    if normalized == "HAV":
        return "Hepatitis A virus"
    if normalized == "HBV":
        return "Hepatitis B virus"
    if normalized == "HCV":
        return "Hepatitis C virus"
    if normalized == "HDV":
        return "Hepatitis D virus"
    if normalized == "HEV":
        return "Hepatitis E virus"
    if normalized:
        return f"Hepatovirus ({normalized})"
    return str(fallback or "Hepatovirus").strip() or "Hepatovirus"


def _hepatovirus_row_species_label(row: dict[str, str], broad_type: str = "") -> str:
    explicit_label = str(row.get("species_label") or "").strip()
    if explicit_label:
        return explicit_label
    virus_name = str(row.get("virus_name") or "").strip()
    if ";" in virus_name:
        alias = str(virus_name.split(";")[-1] or "").strip()
        if alias:
            return alias
    species_name = str(row.get("species") or "").strip()
    return _hepatovirus_species_label(broad_type, species_name or virus_name or "Hepatovirus")


def _infer_hepatitis_subtype_from_row(row: dict[str, str], broad_type: str) -> str:
    normalized_broad = str(broad_type or "").strip().upper()
    genotype = str(row.get("genotype") or "").strip()
    if genotype:
        return genotype
    search_text = " ".join(
        str(row.get(key) or "")
        for key in ("virus_name", "header", "species", "abbrev")
    )
    genotype_match = re.search(r"\bgenotype\s*([0-9]+[A-Za-z]?|[A-Za-z])\b", search_text, flags=re.IGNORECASE)
    if genotype_match:
        return genotype_match.group(1)
    subtype_match = re.search(r"\bsubtype\s*([0-9]+[A-Za-z]?|[A-Za-z])\b", search_text, flags=re.IGNORECASE)
    if subtype_match:
        return subtype_match.group(1)
    abbrev = str(row.get("abbrev") or "").strip()
    if normalized_broad == "HBV" and abbrev.upper().startswith("HBV-"):
        return abbrev.split("-", 1)[1].strip()
    if normalized_broad == "HCV" and abbrev.upper().startswith("HCV"):
        return abbrev[3:].strip("-_ ")
    if normalized_broad == "HDV" and abbrev.upper().startswith("HDV"):
        return abbrev[3:].strip("-_ ")
    return abbrev


def _build_hepatovirus_records(rows: list[dict[str, str]], *, subtype: bool = False, subtype_broad_type: str = "") -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for row in rows:
        fasta_path = Path(str(row.get("fasta_path") or "").strip())
        if not fasta_path.is_file() or fasta_path.stat().st_size == 0:
            continue
        header, sequence = _read_first_fasta_sequence(fasta_path)
        accession = str(row.get("accession") or "").strip()
        if not accession or not sequence:
            continue
        header_token = str(header or "").split()[0].strip()
        if subtype:
            fasta_id = header_token or accession
            broad_type = str(subtype_broad_type or row.get("broad_type") or row.get("abbrev") or "").strip()
        else:
            broad_type = str(row.get("broad_type") or row.get("abbrev") or "").strip()
            fasta_id = header_token or accession
        records.append(
            {
                "fasta_id": fasta_id,
                "accession": accession,
                "sequence": sequence,
                "reference_length": len(sequence),
                "meta": dict(row),
                "header": header,
                "broad_type": broad_type,
                "subtype": _infer_hepatitis_subtype_from_row(row, broad_type) if subtype else "",
            }
        )
    return records


def _select_hepatovirus_reference_records(
    records: list[dict[str, object]],
    screening_dir: Path,
    file_prefix: str,
    single_fastq: str = "",
    fq1: str = "",
    fq2: str = "",
    long_type: str = "",
    threads: int = 4,
    query_fasta: Path | None = None,
    logf=None,
) -> dict[str, object]:
    candidate_fasta = screening_dir / f"{file_prefix}.candidate_references.fasta"
    coverage_prefix = screening_dir / file_prefix
    meta_by_id = _write_reference_records(records, candidate_fasta)
    if query_fasta is not None and query_fasta.is_file() and query_fasta.stat().st_size > 0:
        coverage_rows = _run_multi_reference_blast(
            query_fasta,
            candidate_fasta,
            coverage_prefix,
            meta_by_id,
            threads=threads,
            logf=logf,
        )
        coverage_rows.sort(key=_blast_row_score, reverse=True)
    else:
        coverage_rows = _run_multi_reference_coverage(
            candidate_fasta,
            coverage_prefix,
            single_fastq=single_fastq,
            fq1=fq1,
            fq2=fq2,
            long_type=long_type,
            threads=threads,
            logf=logf,
        )
        coverage_rows.sort(key=_coverage_row_score, reverse=True)
    coverage_path = coverage_prefix.with_suffix(".coverage.tsv")
    best_row = coverage_rows[0] if coverage_rows else {}
    best_meta = meta_by_id.get(str(best_row.get("reference_name") or "").strip()) if best_row else None
    return {
        "candidate_fasta": candidate_fasta,
        "coverage_path": coverage_path,
        "coverage_rows": coverage_rows,
        "best_row": best_row,
        "best_meta": best_meta,
    }


def resolve_hepatovirus_reference(
    pre: str,
    species: str = "",
    requested_ref: str = "",
    single_fastq: str = "",
    fq1: str = "",
    fq2: str = "",
    long_type: str = "",
    threads: int = 4,
    query_fasta: str = "",
    output_dir: Path | None = None,
    logf=None,
) -> dict[str, str]:
    requested = str(requested_ref or "").strip()
    explicit_broad_type = _infer_hepatovirus_broad_type_from_label(species)
    if requested and Path(requested).is_file():
        return {
            "status": "ready",
            "broad_type": explicit_broad_type,
            "subtype": "",
            "hav_subtype": "",
            "species_label": _hepatovirus_species_label(explicit_broad_type, species or "Hepatovirus"),
            "reference_path": str(Path(requested).resolve()),
            "gff_path": "nogtf",
            "summary_path": "",
            "broad_summary_path": "",
            "subtype_summary_path": "",
        }

    screening_dir = output_dir if output_dir is not None else Path(f"{pre}_hepatovirus_reference_selection")
    screening_dir.mkdir(exist_ok=True)
    broad_summary_path = screening_dir / "broad_typing.tsv"
    final_summary_path = screening_dir / "selection.tsv"
    query_path = Path(query_fasta) if str(query_fasta or "").strip() else None

    broad_rows = _load_hepatovirus_reference_manifest()
    broad_records = _build_hepatovirus_records(broad_rows, subtype=False)
    if not broad_records:
        with final_summary_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sample", "broad_type", "subtype", "hav_subtype", "status", "note", "broad_summary_path", "subtype_summary_path"])
            writer.writerow([pre, explicit_broad_type, "", "", "missing", "未找到 Hepatovirus 属级参考库", "", ""])
        return {
            "status": "missing",
            "broad_type": explicit_broad_type,
            "subtype": "",
            "hav_subtype": "",
            "species_label": _hepatovirus_species_label(explicit_broad_type, species or "Hepatovirus"),
            "reference_path": "",
            "gff_path": "nogtf",
            "summary_path": str(final_summary_path.resolve()),
            "broad_summary_path": "",
            "subtype_summary_path": "",
        }

    broad_pick = _select_hepatovirus_reference_records(
        broad_records,
        screening_dir,
        "broad_typing",
        single_fastq=single_fastq,
        fq1=fq1,
        fq2=fq2,
        long_type=long_type,
        threads=threads,
        query_fasta=query_path,
        logf=logf,
    )
    broad_best_row = broad_pick.get("best_row") or {}
    broad_best_meta = broad_pick.get("best_meta") or {}
    broad_meta_row = dict(broad_best_meta.get("meta") or {})
    broad_type = str(broad_best_meta.get("broad_type") or explicit_broad_type).strip()
    broad_species_label = _hepatovirus_row_species_label(broad_meta_row, broad_type) if broad_meta_row else _hepatovirus_species_label(broad_type, species or "Hepatovirus")
    broad_reference_path = ""
    if broad_best_meta:
        broad_reference_path = str((screening_dir / f"{_sanitize_tree_label(broad_type or str(broad_best_meta.get('accession') or 'hepatovirus'))}.reference.fasta").resolve())
        if not _extract_named_fasta_record(Path(str(broad_pick.get("candidate_fasta") or "")), str(broad_best_meta.get("fasta_id") or ""), Path(broad_reference_path)):
            broad_reference_path = ""
    broad_note = "基于肝炎病毒 broad 参考库选择命中最优的大亚型参考。"
    with broad_summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "broad_type", "species_label", "accession", "isolate", "coverage", "mean_depth", "covered_bases", "num_reads", "reference_name", "reference_path", "gff_path", "status", "note"])
        writer.writerow([
            pre,
            broad_type,
            broad_species_label,
            str(broad_best_meta.get("accession") or ""),
            str(broad_meta_row.get("isolate") or ""),
            f"{float(broad_best_row.get('coverage') or 0.0):.6f}",
            f"{float(broad_best_row.get('mean_depth') or 0.0):.6f}",
            f"{float(broad_best_row.get('covered_bases') or 0.0):.0f}",
            f"{float(broad_best_row.get('num_reads') or 0.0):.0f}",
            str(broad_meta_row.get("header") or broad_best_meta.get("accession") or ""),
            broad_reference_path,
            str(broad_meta_row.get("gff3_path") or "nogtf"),
            "ready" if broad_best_meta else "missing",
            broad_note if broad_best_meta else "未获得肝炎病毒 broad 大亚型命中结果",
        ])

    if not broad_best_meta:
        with final_summary_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sample", "broad_type", "subtype", "hav_subtype", "status", "note", "broad_summary_path", "subtype_summary_path"])
            writer.writerow([pre, broad_type, "", "", "missing", "未获得肝炎病毒 broad 大亚型命中结果", str(broad_summary_path.resolve()), ""])
        return {
            "status": "missing",
            "broad_type": broad_type,
            "subtype": "",
            "hav_subtype": "",
            "species_label": broad_species_label,
            "reference_path": "",
            "gff_path": "nogtf",
            "summary_path": str(final_summary_path.resolve()),
            "broad_summary_path": str(broad_summary_path.resolve()),
            "subtype_summary_path": "",
        }

    selected_reference_path = broad_reference_path
    selected_gff_path = str(broad_meta_row.get("gff3_path") or "nogtf")
    selected_coverage_row = broad_best_row
    selected_reference_name = str(broad_meta_row.get("header") or broad_best_meta.get("accession") or "")
    selected_accession = str(broad_best_meta.get("accession") or "")
    selected_isolate = str(broad_meta_row.get("isolate") or "")
    selected_species_label = broad_species_label
    subtype = ""
    hav_subtype = ""
    final_note = broad_note
    subtype_summary_resolved = ""

    normalized_broad_type = broad_type.upper()
    subtype_rows = _load_hepatitis_type_reference_manifest(normalized_broad_type)
    subtype_summary_path = screening_dir / ("hav_subtype_typing.tsv" if normalized_broad_type == "HAV" else f"{normalized_broad_type.lower()}_subtype_typing.tsv")
    if subtype_rows:
        subtype_records = _build_hepatovirus_records(subtype_rows, subtype=True, subtype_broad_type=normalized_broad_type)
        if subtype_records:
            subtype_pick = _select_hepatovirus_reference_records(
                subtype_records,
                screening_dir,
                "hav_subtype_typing" if normalized_broad_type == "HAV" else f"{normalized_broad_type.lower()}_subtype_typing",
                single_fastq=single_fastq,
                fq1=fq1,
                fq2=fq2,
                long_type=long_type,
                threads=threads,
                query_fasta=query_path,
                logf=logf,
            )
            subtype_best_row = subtype_pick.get("best_row") or {}
            subtype_best_meta = subtype_pick.get("best_meta") or {}
            subtype_meta_row = dict(subtype_best_meta.get("meta") or {})
            subtype = str(subtype_best_meta.get("subtype") or _infer_hepatitis_subtype_from_row(subtype_meta_row, normalized_broad_type)).strip()
            hav_subtype = subtype if normalized_broad_type == "HAV" else ""
            subtype_reference_path = ""
            if subtype_best_meta:
                subtype_reference_path = str((screening_dir / f"{_sanitize_tree_label(subtype or str(subtype_best_meta.get('accession') or normalized_broad_type.lower()))}.reference.fasta").resolve())
                if not _extract_named_fasta_record(Path(str(subtype_pick.get("candidate_fasta") or "")), str(subtype_best_meta.get("fasta_id") or ""), Path(subtype_reference_path)):
                    subtype_reference_path = ""
            with subtype_summary_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["sample", "broad_type", "subtype", "hav_subtype", "accession", "isolate", "coverage", "mean_depth", "covered_bases", "num_reads", "reference_name", "reference_path", "gff_path", "status", "note"])
                writer.writerow([
                    pre,
                    normalized_broad_type,
                    subtype,
                    hav_subtype,
                    str(subtype_best_meta.get("accession") or ""),
                    str(subtype_meta_row.get("isolate") or ""),
                    f"{float(subtype_best_row.get('coverage') or 0.0):.6f}",
                    f"{float(subtype_best_row.get('mean_depth') or 0.0):.6f}",
                    f"{float(subtype_best_row.get('covered_bases') or 0.0):.0f}",
                    f"{float(subtype_best_row.get('num_reads') or 0.0):.0f}",
                    str(subtype_meta_row.get("header") or subtype_best_meta.get("accession") or ""),
                    subtype_reference_path,
                    str(subtype_meta_row.get("gff3_path") or "nogtf"),
                    "ready" if subtype_best_meta else "missing",
                    f"基于 {normalized_broad_type} 亚型参考库选择覆盖度最优的子亚型参考。" if subtype_best_meta else f"{normalized_broad_type} 大亚型已命中，但未获得稳定的子亚型命中结果",
                ])
            subtype_summary_resolved = str(subtype_summary_path.resolve())
            if subtype_best_meta and subtype_reference_path:
                selected_reference_path = subtype_reference_path
                selected_gff_path = str(subtype_meta_row.get("gff3_path") or selected_gff_path)
                selected_coverage_row = subtype_best_row
                selected_reference_name = str(subtype_meta_row.get("header") or subtype_best_meta.get("accession") or "")
                selected_accession = str(subtype_best_meta.get("accession") or "")
                selected_isolate = str(subtype_meta_row.get("isolate") or "")
                selected_species_label = _hepatovirus_species_label(normalized_broad_type, selected_species_label)
                final_note = f"先基于肝炎病毒 broad 参考库完成大亚型判定；命中 {normalized_broad_type} 后，再在 {normalized_broad_type} 亚型参考库中选择覆盖度最优的子亚型参考。"

    with final_summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "species_label", "broad_type", "subtype", "hav_subtype", "isolate", "coverage", "mean_depth", "covered_bases", "num_reads", "reference_name", "accession", "reference_path", "gff_path", "status", "note", "broad_summary_path", "subtype_summary_path"])
        writer.writerow([
            pre,
            selected_species_label,
            broad_type,
            subtype,
            hav_subtype,
            selected_isolate,
            f"{float(selected_coverage_row.get('coverage') or 0.0):.6f}",
            f"{float(selected_coverage_row.get('mean_depth') or 0.0):.6f}",
            f"{float(selected_coverage_row.get('covered_bases') or 0.0):.0f}",
            f"{float(selected_coverage_row.get('num_reads') or 0.0):.0f}",
            selected_reference_name,
            selected_accession,
            selected_reference_path,
            selected_gff_path,
            "ready" if selected_reference_path else "missing",
            final_note,
            str(broad_summary_path.resolve()),
            subtype_summary_resolved,
        ])
    return {
        "status": "ready" if selected_reference_path else "missing",
        "broad_type": broad_type,
        "subtype": subtype,
        "hav_subtype": hav_subtype,
        "species_label": selected_species_label,
        "reference_path": selected_reference_path,
        "gff_path": selected_gff_path,
        "summary_path": str(final_summary_path.resolve()),
        "broad_summary_path": str(broad_summary_path.resolve()),
        "subtype_summary_path": subtype_summary_resolved,
    }


def run_hepatovirus_consensus_typing(pre: str, consensus_fasta: Path, logf=None) -> dict[str, str]:
    if not consensus_fasta.is_file() or consensus_fasta.stat().st_size == 0:
        return {"status": "missing", "summary_path": "", "broad_type": "", "subtype": "", "hav_subtype": ""}
    output_dir = Path(f"{pre}_hepatovirus_reference_selection") / "consensus_typing"
    selection = resolve_hepatovirus_reference(
        pre,
        species="Hepatovirus",
        requested_ref="",
        query_fasta=str(consensus_fasta),
        output_dir=output_dir,
        logf=logf,
    )
    summary_path = output_dir / "consensus_typing.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "broad_type", "subtype", "hav_subtype", "selection_summary_path"])
        writer.writerow([
            pre,
            selection.get("broad_type", ""),
            selection.get("subtype", ""),
            selection.get("hav_subtype", ""),
            selection.get("summary_path", ""),
        ])
    return {
        "status": "ready",
        "summary_path": str(summary_path.resolve()),
        "broad_type": str(selection.get("broad_type") or ""),
        "subtype": str(selection.get("subtype") or ""),
        "hav_subtype": str(selection.get("hav_subtype") or ""),
    }


def _choose_minimap2_preset(long_type: str) -> str:
    normalized = str(long_type or "").strip().lower()
    if normalized in {"pb", "pacbio", "hifi", "clr"}:
        return "map-pb"
    return "map-ont"


def _parse_samtools_coverage_tsv(path: Path) -> dict[str, float]:
    if not path.is_file() or path.stat().st_size == 0:
        return {"coverage": 0.0, "mean_depth": 0.0, "covered_bases": 0.0, "num_reads": 0.0}
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        row = handle.readline().rstrip("\n").split("\t")
    if not header or not row or len(row) != len(header):
        return {"coverage": 0.0, "mean_depth": 0.0, "covered_bases": 0.0, "num_reads": 0.0}
    values = dict(zip(header, row))
    return {
        "coverage": float(values.get("coverage") or 0.0),
        "mean_depth": float(values.get("meandepth") or 0.0),
        "covered_bases": float(values.get("covbases") or 0.0),
        "num_reads": float(values.get("numreads") or 0.0),
    }


def _parse_samtools_coverage_rows(path: Path) -> list[dict[str, object]]:
    if not path.is_file() or path.stat().st_size == 0:
        return []
    rows: list[dict[str, object]] = []
    with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if not row:
                continue
            reference_name = str(row.get("#rname") or row.get("rname") or row.get("name") or "").strip()
            if not reference_name or reference_name == "*":
                continue
            rows.append(
                {
                    "reference_name": reference_name,
                    "coverage": float(row.get("coverage") or 0.0),
                    "mean_depth": float(row.get("meandepth") or 0.0),
                    "covered_bases": float(row.get("covbases") or 0.0),
                    "num_reads": float(row.get("numreads") or 0.0),
                }
            )
    return rows


def _coverage_row_score(row: dict[str, object]) -> tuple[float, float, float, float]:
    return (
        float(row.get("coverage") or 0.0),
        float(row.get("mean_depth") or 0.0),
        float(row.get("covered_bases") or 0.0),
        float(row.get("num_reads") or 0.0),
    )


def _blast_row_score(row: dict[str, object]) -> tuple[float, float, float, float]:
    return (
        float(row.get("coverage") or 0.0),
        float(row.get("identity") or 0.0),
        float(row.get("bitscore") or 0.0),
        float(row.get("covered_bases") or 0.0),
    )


def _merged_interval_length(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    merged: list[tuple[int, int]] = []
    for start, end in sorted((min(a, b), max(a, b)) for a, b in intervals):
        if not merged or start > merged[-1][1] + 1:
            merged.append((start, end))
            continue
        previous_start, previous_end = merged[-1]
        merged[-1] = (previous_start, max(previous_end, end))
    return sum((end - start + 1) for start, end in merged)


def _write_reference_records(records: list[dict[str, object]], fasta_path: Path) -> dict[str, dict[str, object]]:
    meta_by_id: dict[str, dict[str, object]] = {}
    with fasta_path.open("w", encoding="utf-8") as handle:
        for record in records:
            fasta_id = str(record.get("fasta_id") or "").strip()
            sequence = str(record.get("sequence") or "").strip()
            if not fasta_id or not sequence:
                continue
            meta_by_id[fasta_id] = record
            handle.write(f">{fasta_id}\n")
            for index in range(0, len(sequence), 80):
                handle.write(sequence[index:index + 80] + "\n")
    return meta_by_id


def _run_multi_reference_blast(
    query_fasta: Path,
    reference_fasta: Path,
    output_prefix: Path,
    meta_by_id: dict[str, dict[str, object]],
    threads: int = 4,
    logf=None,
) -> list[dict[str, object]]:
    blast_path = output_prefix.with_suffix(".blastn.tsv")
    coverage_path = output_prefix.with_suffix(".coverage.tsv")
    blastn_bin = _resolve_hadv_blastn_bin()
    cmd = " ".join(
        [
            shlex.quote(blastn_bin),
            "-task blastn",
            "-query",
            shlex.quote(str(query_fasta)),
            "-subject",
            shlex.quote(str(reference_fasta)),
            "-outfmt",
            shlex.quote("6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"),
            "-evalue 1e-20",
            "-max_target_seqs 200",
            "-dust no",
            "-num_threads",
            str(max(1, int(threads or 1))),
            "-out",
            shlex.quote(str(blast_path)),
        ]
    )
    run_command(cmd, logf=logf)
    buckets: dict[str, dict[str, object]] = {}
    if blast_path.is_file() and blast_path.stat().st_size > 0:
        with blast_path.open("r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 12:
                    continue
                qseqid, sseqid, pident, length, _mismatch, _gapopen, _qstart, _qend, sstart, send, evalue, bitscore = parts[:12]
                subject = str(sseqid).strip()
                meta = meta_by_id.get(subject) or {}
                ref_len = int(meta.get("reference_length") or 0)
                try:
                    align_len = float(length)
                    identity = float(pident)
                    bit_score = float(bitscore)
                    start_i = int(float(sstart))
                    end_i = int(float(send))
                except ValueError:
                    continue
                bucket = buckets.setdefault(
                    subject,
                    {
                        "reference_name": subject,
                        "reference_length": ref_len,
                        "intervals": [],
                        "identity_weighted_sum": 0.0,
                        "aligned_bases_sum": 0.0,
                        "bitscore_sum": 0.0,
                        "best_evalue": "",
                        "query_ids": set(),
                    },
                )
                bucket["intervals"].append((start_i, end_i))
                bucket["identity_weighted_sum"] = float(bucket["identity_weighted_sum"]) + (identity * align_len)
                bucket["aligned_bases_sum"] = float(bucket["aligned_bases_sum"]) + align_len
                bucket["bitscore_sum"] = float(bucket["bitscore_sum"]) + bit_score
                bucket["query_ids"].add(str(qseqid).strip())
                previous_evalue = str(bucket.get("best_evalue") or "")
                try:
                    if not previous_evalue or float(evalue) < float(previous_evalue):
                        bucket["best_evalue"] = str(evalue)
                except ValueError:
                    if not previous_evalue:
                        bucket["best_evalue"] = str(evalue)
    rows: list[dict[str, object]] = []
    for subject, bucket in buckets.items():
        ref_len = int(bucket.get("reference_length") or 0)
        covered_bases = _merged_interval_length(bucket.get("intervals") or [])
        coverage = min(100.0, (covered_bases / ref_len) * 100.0) if ref_len else 0.0
        aligned_bases = float(bucket.get("aligned_bases_sum") or 0.0)
        identity = (float(bucket.get("identity_weighted_sum") or 0.0) / aligned_bases) if aligned_bases else 0.0
        rows.append(
            {
                "reference_name": subject,
                "identity": f"{identity:.6f}",
                "coverage": f"{coverage:.6f}",
                "mean_depth": "",
                "covered_bases": f"{covered_bases:.0f}",
                "num_reads": "",
                "bitscore": f"{float(bucket.get('bitscore_sum') or 0.0):.6f}",
                "evalue": str(bucket.get("best_evalue") or ""),
                "method": "blastn",
            }
        )
    rows = sorted(rows, key=_blast_row_score, reverse=True)
    with coverage_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["reference_name", "identity", "coverage", "mean_depth", "covered_bases", "num_reads", "bitscore", "evalue", "method"])
        for row in rows:
            writer.writerow([
                row.get("reference_name", ""),
                row.get("identity", ""),
                row.get("coverage", ""),
                row.get("mean_depth", ""),
                row.get("covered_bases", ""),
                row.get("num_reads", ""),
                row.get("bitscore", ""),
                row.get("evalue", ""),
                row.get("method", ""),
            ])
    return rows


def _run_multi_reference_coverage(
    reference_fasta: Path,
    output_prefix: Path,
    single_fastq: str = "",
    fq1: str = "",
    fq2: str = "",
    long_type: str = "",
    threads: int = 4,
    logf=None,
) -> list[dict[str, object]]:
    bam_path = output_prefix.with_suffix(".bam")
    coverage_path = output_prefix.with_suffix(".coverage.tsv")
    if fq1:
        map_cmd = " ".join(
            [
                "minimap2",
                "-ax",
                "sr",
                "-t",
                str(max(1, int(threads or 1))),
                shlex.quote(str(reference_fasta)),
                shlex.quote(str(fq1)),
                *([shlex.quote(str(fq2))] if fq2 else []),
                "|",
                "samtools",
                "sort",
                "-o",
                shlex.quote(str(bam_path)),
            ]
        )
    else:
        map_cmd = " ".join(
            [
                "minimap2",
                "-ax",
                _choose_minimap2_preset(long_type),
                shlex.quote(str(reference_fasta)),
                shlex.quote(str(single_fastq)),
                "-t",
                str(max(1, int(threads or 1))),
                "|",
                "samtools",
                "sort",
                "-o",
                shlex.quote(str(bam_path)),
            ]
        )
    run_command(map_cmd, logf=logf)
    run_command(f"samtools index -@ {max(1, int(threads or 1))} {shlex.quote(str(bam_path))}", logf=logf)
    run_command(f"samtools coverage {shlex.quote(str(bam_path))} > {shlex.quote(str(coverage_path))}", logf=logf)
    return _parse_samtools_coverage_rows(coverage_path)


def _write_aggregated_coverage_rows(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "#rname",
                "startpos",
                "endpos",
                "numreads",
                "covbases",
                "coverage",
                "meandepth",
                "stage",
                "batch_id",
                "is_anchor",
                "member_count",
                "member_accessions",
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    row.get("reference_name", ""),
                    1,
                    row.get("reference_length", ""),
                    row.get("num_reads", 0),
                    row.get("covered_bases", 0),
                    row.get("coverage", 0),
                    row.get("mean_depth", 0),
                    row.get("stage", ""),
                    row.get("batch_id", ""),
                    row.get("is_anchor", ""),
                    row.get("member_count", ""),
                    row.get("member_accessions", ""),
                ]
            )


def _chunk_records(records: list[dict[str, object]], size: int) -> list[list[dict[str, object]]]:
    if size <= 0:
        return [records]
    return [records[index:index + size] for index in range(0, len(records), size)]


def _resolve_mmseqs_bin() -> str:
    env_value = str(os.environ.get("META_MMSEQS_BIN") or "").strip()
    candidates = [env_value] if env_value else []
    candidates.extend(
        [
            conda_env_path("ncov", "bin", "mmseqs"),
            shutil.which("mmseqs") or "",
        ]
    )
    for candidate in candidates:
        if not candidate:
            continue
        path = Path(candidate).expanduser()
        if path.is_file():
            return str(path.resolve())
    return ""


def _build_enterovirus_candidate_records(candidates: list[dict[str, str]]) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for row in candidates:
        accession = str(row.get("accession") or row.get("accession_full") or "").split(".", 1)[0]
        fasta_path = Path(str(row.get("fasta_path") or "").strip())
        if not accession or not fasta_path.is_file():
            continue
        source_record_id = str(row.get("source_record_id") or row.get("accession_full") or accession)
        record = _get_fasta_record_by_id(fasta_path, source_record_id)
        if record is None:
            continue
        sequence = str(record.seq).strip().upper()
        if not sequence:
            continue
        records.append(
            {
                "fasta_id": _sanitize_tree_label(accession),
                "accession": accession,
                "sequence": sequence,
                "reference_length": len(sequence),
                "meta": row,
                "source_record_id": source_record_id,
            }
        )
    return records


def _deduplicate_enterovirus_candidate_records_exact(records: list[dict[str, object]]) -> list[dict[str, object]]:
    dedup_by_sequence: dict[str, dict[str, object]] = {}
    for record in records:
        sequence = str(record.get("sequence") or "")
        if not sequence:
            continue
        bucket = dedup_by_sequence.get(sequence)
        if bucket is None:
            bucket = dict(record)
            bucket["members"] = [record]
            dedup_by_sequence[sequence] = bucket
            continue
        bucket.setdefault("members", []).append(record)
    deduped: list[dict[str, object]] = []
    for item in dedup_by_sequence.values():
        members = list(item.get("members") or [])
        member_accessions = [str(member.get("accession") or "") for member in members if str(member.get("accession") or "").strip()]
        item["member_count"] = len(member_accessions)
        item["member_accessions"] = ",".join(sorted(member_accessions))
        deduped.append(item)
    deduped.sort(key=lambda item: (str(item.get("accession") or ""), str(item.get("fasta_id") or "")))
    return deduped


def _deduplicate_enterovirus_candidate_records(
    records: list[dict[str, object]],
    screening_dir: Path,
    threads: int = 4,
    logf=None,
) -> list[dict[str, object]]:
    if not records:
        return []
    mmseqs_bin = _resolve_mmseqs_bin()
    if not mmseqs_bin:
        return _deduplicate_enterovirus_candidate_records_exact(records)
    cluster_dir = screening_dir / "candidate_mmseqs_dedup"
    cluster_dir.mkdir(exist_ok=True)
    input_fasta = cluster_dir / "input.fasta"
    record_lookup = _write_reference_records(records, input_fasta)
    result_prefix = cluster_dir / "clusters"
    tmp_dir = cluster_dir / "tmp"
    cmd = " ".join(
        [
            shlex.quote(mmseqs_bin),
            "easy-cluster",
            shlex.quote(str(input_fasta)),
            shlex.quote(str(result_prefix)),
            shlex.quote(str(tmp_dir)),
            "--dbtype",
            "2",
            "--min-seq-id",
            f"{ENTEROVIRUS_REFERENCE_DEDUP_SEQ_ID:.3f}",
            "-c",
            f"{ENTEROVIRUS_REFERENCE_DEDUP_SEQ_ID:.3f}",
            "--cov-mode",
            "5",
            "--cluster-mode",
            "2",
            "--similarity-type",
            "2",
            "--threads",
            str(max(1, int(threads or 1))),
            "-v",
            "1",
        ]
    )
    try:
        run_command(cmd, logf=logf)
    except Exception:
        return _deduplicate_enterovirus_candidate_records_exact(records)
    cluster_tsv = cluster_dir / "clusters_cluster.tsv"
    if not cluster_tsv.is_file() or cluster_tsv.stat().st_size == 0:
        return _deduplicate_enterovirus_candidate_records_exact(records)
    members_by_rep: dict[str, list[dict[str, object]]] = defaultdict(list)
    assigned_ids: set[str] = set()
    with cluster_tsv.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            parts = raw_line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            rep_id = str(parts[0] or "").strip()
            member_id = str(parts[1] or "").strip()
            if not rep_id or not member_id:
                continue
            member = record_lookup.get(member_id)
            representative = record_lookup.get(rep_id)
            if member is None or representative is None:
                continue
            members_by_rep[rep_id].append(member)
            assigned_ids.add(member_id)
            assigned_ids.add(rep_id)
    deduped: list[dict[str, object]] = []
    for rep_id, members in members_by_rep.items():
        representative = dict(record_lookup[rep_id])
        member_accessions = sorted(
            {
                str(member.get("accession") or "").strip()
                for member in members
                if str(member.get("accession") or "").strip()
            }
        )
        representative["members"] = list(members)
        representative["member_count"] = len(member_accessions) or 1
        representative["member_accessions"] = ",".join(member_accessions) if member_accessions else str(representative.get("accession") or "")
        deduped.append(representative)
    for record in records:
        fasta_id = str(record.get("fasta_id") or "").strip()
        if not fasta_id or fasta_id in assigned_ids:
            continue
        singleton = dict(record)
        singleton["members"] = [record]
        singleton["member_count"] = 1
        singleton["member_accessions"] = str(record.get("accession") or "")
        deduped.append(singleton)
    deduped.sort(key=lambda item: (str(item.get("accession") or ""), str(item.get("fasta_id") or "")))
    return deduped


def _select_enterovirus_candidates_in_batches(
    records: list[dict[str, object]],
    screening_dir: Path,
    single_fastq: str = "",
    fq1: str = "",
    fq2: str = "",
    long_type: str = "",
    threads: int = 4,
    logf=None,
) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, object] | None]:
    if not records:
        return [], [], None
    batches_dir = screening_dir / "candidate_batches"
    batches_dir.mkdir(exist_ok=True)
    dedup_fasta = screening_dir / "candidate_references.dedup.fasta"
    dedup_meta = _write_reference_records(records, dedup_fasta)
    dedup_rows = _run_multi_reference_coverage(
        dedup_fasta,
        screening_dir / "candidate_references.dedup",
        single_fastq=single_fastq,
        fq1=fq1,
        fq2=fq2,
        long_type=long_type,
        threads=threads,
        logf=logf,
    )
    for row in dedup_rows:
        row["stage"] = "dedup_screen"
        row["batch_id"] = "dedup"
        meta = dedup_meta.get(str(row.get("reference_name") or "").strip(), {})
        row["member_count"] = meta.get("member_count", 1)
        row["member_accessions"] = meta.get("member_accessions", meta.get("accession", ""))
        row["reference_length"] = meta.get("reference_length", "")
    best_anchor_row = None
    for row in dedup_rows:
        if best_anchor_row is None or _coverage_row_score(row) > _coverage_row_score(best_anchor_row):
            best_anchor_row = row
    if best_anchor_row is None:
        return [], dedup_rows, None
    anchor_id = str(best_anchor_row.get("reference_name") or "").strip()
    anchor_record = dedup_meta.get(anchor_id)
    if anchor_record is None:
        return [], dedup_rows, None
    other_records = [record for record in records if str(record.get("fasta_id") or "").strip() != anchor_id]
    batch_size = max(2, ENTEROVIRUS_REFERENCE_BATCH_SIZE)
    batch_rows_by_ref: dict[str, dict[str, object]] = {}
    round_rows: list[dict[str, object]] = []
    for batch_index, batch_records in enumerate(_chunk_records(other_records, max(1, batch_size - 1)), start=1):
        round_records = [anchor_record] + batch_records
        batch_fasta = batches_dir / f"batch_{batch_index:03d}.fasta"
        batch_meta = _write_reference_records(round_records, batch_fasta)
        rows = _run_multi_reference_coverage(
            batch_fasta,
            batches_dir / f"batch_{batch_index:03d}",
            single_fastq=single_fastq,
            fq1=fq1,
            fq2=fq2,
            long_type=long_type,
            threads=threads,
            logf=logf,
        )
        for row in rows:
            ref_name = str(row.get("reference_name") or "").strip()
            meta = batch_meta.get(ref_name, {})
            row["stage"] = "batched_screen"
            row["batch_id"] = f"batch_{batch_index:03d}"
            row["is_anchor"] = "yes" if ref_name == anchor_id else "no"
            row["member_count"] = meta.get("member_count", 1)
            row["member_accessions"] = meta.get("member_accessions", meta.get("accession", ""))
            row["reference_length"] = meta.get("reference_length", "")
            round_rows.append(row)
            current = batch_rows_by_ref.get(ref_name)
            if current is None or _coverage_row_score(row) > _coverage_row_score(current):
                batch_rows_by_ref[ref_name] = dict(row)
    playoff_candidates = sorted(batch_rows_by_ref.values(), key=_coverage_row_score, reverse=True)[: max(2, ENTEROVIRUS_REFERENCE_PLAYOFF_SIZE)]
    playoff_records = [dedup_meta[str(row.get("reference_name") or "").strip()] for row in playoff_candidates if str(row.get("reference_name") or "").strip() in dedup_meta]
    final_rows: list[dict[str, object]] = []
    if playoff_records:
        playoff_fasta = screening_dir / "candidate_references.playoff.fasta"
        playoff_meta = _write_reference_records(playoff_records, playoff_fasta)
        playoff_rows = _run_multi_reference_coverage(
            playoff_fasta,
            screening_dir / "candidate_references.playoff",
            single_fastq=single_fastq,
            fq1=fq1,
            fq2=fq2,
            long_type=long_type,
            threads=threads,
            logf=logf,
        )
        for row in playoff_rows:
            ref_name = str(row.get("reference_name") or "").strip()
            meta = playoff_meta.get(ref_name, {})
            row["stage"] = "playoff"
            row["batch_id"] = "playoff"
            row["is_anchor"] = "yes" if ref_name == anchor_id else "no"
            row["member_count"] = meta.get("member_count", 1)
            row["member_accessions"] = meta.get("member_accessions", meta.get("accession", ""))
            row["reference_length"] = meta.get("reference_length", "")
            final_rows.append(row)
    else:
        final_rows = sorted(batch_rows_by_ref.values(), key=_coverage_row_score, reverse=True)
    return final_rows, dedup_rows + round_rows, anchor_record


def resolve_rsv_reference(
    pre: str,
    species: str = "",
    requested_ref: str = "",
    single_fastq: str = "",
    fq1: str = "",
    fq2: str = "",
    long_type: str = "",
    threads: int = 4,
    logf=None,
) -> dict[str, str]:
    requested = str(requested_ref or "").strip()
    explicit_type = _infer_rsv_type_from_label(species)
    if explicit_type in {"A", "B"}:
        assets = _resolve_rsv_reference_assets(explicit_type)
        return {
            "status": "ready",
            "rsv_type": explicit_type,
            "species_label": _rsv_species_label(explicit_type),
            "reference_path": str(assets["reference_fasta"]),
            "gff_path": str(assets["annotation_gff"]) if assets["annotation_gff"].is_file() else "nogtf",
            "dataset_dir": str(assets["dataset_dir"]),
            "summary_path": "",
        }
    if requested and Path(requested).is_file():
        lowered = requested.lower()
        if "rsv_a" in lowered:
            assets = _resolve_rsv_reference_assets("A")
            return {
                "status": "ready",
                "rsv_type": "A",
                "species_label": _rsv_species_label("A"),
                "reference_path": str(Path(requested).resolve()),
                "gff_path": str(assets["annotation_gff"]) if assets["annotation_gff"].is_file() else "nogtf",
                "dataset_dir": str(assets["dataset_dir"]),
                "summary_path": "",
            }
        if "rsv_b" in lowered:
            assets = _resolve_rsv_reference_assets("B")
            return {
                "status": "ready",
                "rsv_type": "B",
                "species_label": _rsv_species_label("B"),
                "reference_path": str(Path(requested).resolve()),
                "gff_path": str(assets["annotation_gff"]) if assets["annotation_gff"].is_file() else "nogtf",
                "dataset_dir": str(assets["dataset_dir"]),
                "summary_path": "",
            }

    screening_dir = Path(f"{pre}_rsv_reference_selection")
    screening_dir.mkdir(exist_ok=True)
    results: list[dict[str, object]] = []
    for rsv_type in ["A", "B"]:
        assets = _resolve_rsv_reference_assets(rsv_type)
        reference_fasta = assets["reference_fasta"]
        if not reference_fasta.is_file():
            continue
        type_dir = screening_dir / f"rsv_{rsv_type.lower()}"
        type_dir.mkdir(parents=True, exist_ok=True)
        bam_path = type_dir / "screening.bam"
        coverage_path = type_dir / "coverage.tsv"
        if not bam_path.is_file():
            if fq1:
                map_cmd = " ".join(
                    [
                        "minimap2",
                        "-ax",
                        "sr",
                        shlex.quote(str(reference_fasta)),
                        shlex.quote(str(fq1)),
                        *( [shlex.quote(str(fq2))] if fq2 else [] ),
                        "|",
                        "samtools",
                        "sort",
                        "-o",
                        shlex.quote(str(bam_path)),
                    ]
                )
            else:
                map_cmd = " ".join(
                    [
                        "minimap2",
                        "-ax",
                        _choose_minimap2_preset(long_type),
                        shlex.quote(str(reference_fasta)),
                        shlex.quote(str(single_fastq)),
                        "-t",
                        str(max(1, int(threads or 1))),
                        "|",
                        "samtools",
                        "sort",
                        "-o",
                        shlex.quote(str(bam_path)),
                    ]
                )
            run_command(map_cmd, logf=logf)
        if bam_path.is_file():
            run_command(f"samtools index -@ {max(1, int(threads or 1))} {shlex.quote(str(bam_path))}", logf=logf)
            run_command(f"samtools coverage {shlex.quote(str(bam_path))} > {shlex.quote(str(coverage_path))}", logf=logf)
        coverage_info = _parse_samtools_coverage_tsv(coverage_path)
        results.append(
            {
                "rsv_type": rsv_type,
                "species_label": _rsv_species_label(rsv_type),
                "reference_path": str(reference_fasta),
                "gff_path": str(assets["annotation_gff"]) if assets["annotation_gff"].is_file() else "nogtf",
                "dataset_dir": str(assets["dataset_dir"]),
                "coverage": coverage_info["coverage"],
                "mean_depth": coverage_info["mean_depth"],
                "covered_bases": coverage_info["covered_bases"],
                "num_reads": coverage_info["num_reads"],
            }
        )
    summary_path = screening_dir / "selection.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "rsv_type", "coverage", "mean_depth", "covered_bases", "num_reads", "reference_path", "gff_path", "dataset_dir"])
        for item in results:
            writer.writerow([
                pre,
                item["rsv_type"],
                f"{float(item['coverage']):.6f}",
                f"{float(item['mean_depth']):.6f}",
                f"{float(item['covered_bases']):.0f}",
                f"{float(item['num_reads']):.0f}",
                item["reference_path"],
                item["gff_path"],
                item["dataset_dir"],
            ])
    if not results:
        return {"status": "missing", "rsv_type": "-", "species_label": species or "Respiratory syncytial virus", "reference_path": "", "gff_path": "nogtf", "dataset_dir": "", "summary_path": str(summary_path)}
    results.sort(key=lambda item: (float(item["coverage"]), float(item["mean_depth"]), float(item["covered_bases"]), float(item["num_reads"])), reverse=True)
    best = results[0]
    return {
        "status": "ready",
        "rsv_type": str(best["rsv_type"]),
        "species_label": str(best["species_label"]),
        "reference_path": str(best["reference_path"]),
        "gff_path": str(best["gff_path"]),
        "dataset_dir": str(best["dataset_dir"]),
        "summary_path": str(summary_path.resolve()),
    }


def resolve_denv_reference(
    pre: str,
    species: str = "",
    requested_ref: str = "",
    single_fastq: str = "",
    fq1: str = "",
    fq2: str = "",
    long_type: str = "",
    threads: int = 4,
    logf=None,
) -> dict[str, str]:
    requested = str(requested_ref or "").strip()
    explicit_type = _infer_denv_type_from_label(species)
    if explicit_type in {"1", "2", "3", "4"}:
        assets = _resolve_denv_reference_assets(explicit_type)
        return {
            "status": "ready",
            "denv_type": explicit_type,
            "species_label": _denv_species_label(explicit_type),
            "reference_path": str(assets["reference_fasta"]),
            "gff_path": str(assets["annotation_gff"]) if assets["annotation_gff"].is_file() else "nogtf",
            "dataset_dir": str(assets["dataset_dir"]),
            "summary_path": "",
        }
    if requested and Path(requested).is_file():
        lowered = requested.lower()
        for denv_type in ["1", "2", "3", "4"]:
            if f"denv{denv_type}" in lowered:
                assets = _resolve_denv_reference_assets(denv_type)
                return {
                    "status": "ready",
                    "denv_type": denv_type,
                    "species_label": _denv_species_label(denv_type),
                    "reference_path": str(Path(requested).resolve()),
                    "gff_path": str(assets["annotation_gff"]) if assets["annotation_gff"].is_file() else "nogtf",
                    "dataset_dir": str(assets["dataset_dir"]),
                    "summary_path": "",
                }

    screening_dir = Path(f"{pre}_denv_reference_selection")
    screening_dir.mkdir(exist_ok=True)
    results: list[dict[str, object]] = []
    for denv_type in ["1", "2", "3", "4"]:
        assets = _resolve_denv_reference_assets(denv_type)
        reference_fasta = assets["reference_fasta"]
        if not reference_fasta.is_file():
            continue
        type_dir = screening_dir / f"denv{denv_type}"
        type_dir.mkdir(parents=True, exist_ok=True)
        bam_path = type_dir / "screening.bam"
        coverage_path = type_dir / "coverage.tsv"
        if not bam_path.is_file():
            if fq1:
                map_cmd = " ".join(
                    [
                        "minimap2",
                        "-ax",
                        "sr",
                        shlex.quote(str(reference_fasta)),
                        shlex.quote(str(fq1)),
                        *([shlex.quote(str(fq2))] if fq2 else []),
                        "|",
                        "samtools",
                        "sort",
                        "-o",
                        shlex.quote(str(bam_path)),
                    ]
                )
            else:
                map_cmd = " ".join(
                    [
                        "minimap2",
                        "-ax",
                        _choose_minimap2_preset(long_type),
                        shlex.quote(str(reference_fasta)),
                        shlex.quote(str(single_fastq)),
                        "-t",
                        str(max(1, int(threads or 1))),
                        "|",
                        "samtools",
                        "sort",
                        "-o",
                        shlex.quote(str(bam_path)),
                    ]
                )
            run_command(map_cmd, logf=logf)
        if bam_path.is_file():
            run_command(f"samtools index -@ {max(1, int(threads or 1))} {shlex.quote(str(bam_path))}", logf=logf)
            run_command(f"samtools coverage {shlex.quote(str(bam_path))} > {shlex.quote(str(coverage_path))}", logf=logf)
        coverage_info = _parse_samtools_coverage_tsv(coverage_path)
        results.append(
            {
                "denv_type": denv_type,
                "species_label": _denv_species_label(denv_type),
                "reference_path": str(reference_fasta),
                "gff_path": str(assets["annotation_gff"]) if assets["annotation_gff"].is_file() else "nogtf",
                "dataset_dir": str(assets["dataset_dir"]),
                "coverage": coverage_info["coverage"],
                "mean_depth": coverage_info["mean_depth"],
                "covered_bases": coverage_info["covered_bases"],
                "num_reads": coverage_info["num_reads"],
            }
        )
    summary_path = screening_dir / "selection.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "denv_type", "coverage", "mean_depth", "covered_bases", "num_reads", "reference_path", "gff_path", "dataset_dir"])
        for item in results:
            writer.writerow([
                pre,
                item["denv_type"],
                f"{float(item['coverage']):.6f}",
                f"{float(item['mean_depth']):.6f}",
                f"{float(item['covered_bases']):.0f}",
                f"{float(item['num_reads']):.0f}",
                item["reference_path"],
                item["gff_path"],
                item["dataset_dir"],
            ])
    if not results:
        return {"status": "missing", "denv_type": "-", "species_label": species or "Dengue virus", "reference_path": "", "gff_path": "nogtf", "dataset_dir": "", "summary_path": str(summary_path)}
    results.sort(key=lambda item: (float(item["coverage"]), float(item["mean_depth"]), float(item["covered_bases"]), float(item["num_reads"])), reverse=True)
    best = results[0]
    return {
        "status": "ready",
        "denv_type": str(best["denv_type"]),
        "species_label": str(best["species_label"]),
        "reference_path": str(best["reference_path"]),
        "gff_path": str(best["gff_path"]),
        "dataset_dir": str(best["dataset_dir"]),
        "summary_path": str(summary_path.resolve()),
    }


def resolve_hpiv_reference(
    pre: str,
    species: str = "",
    requested_ref: str = "",
    single_fastq: str = "",
    fq1: str = "",
    fq2: str = "",
    long_type: str = "",
    threads: int = 4,
    logf=None,
) -> dict[str, str]:
    requested = str(requested_ref or "").strip()
    explicit_type = _infer_hpiv_type_from_label(species)
    if explicit_type in {"1", "2", "3", "4A", "4B"}:
        assets = _resolve_hpiv_reference_assets(explicit_type)
        return {
            "status": "ready",
            "hpiv_type": explicit_type,
            "species_label": _hpiv_species_label(explicit_type),
            "reference_path": str(assets["reference_fasta"]),
            "gff_path": str(assets["annotation_gff"]) if assets["annotation_gff"].is_file() else "nogtf",
            "summary_path": "",
        }
    if requested and Path(requested).is_file():
        lowered = requested.lower()
        path_map = {
            "hpiv1": "1",
            "hpiv2": "2",
            "hpiv3": "3",
            "hpiv4a": "4A",
            "hpiv4b": "4B",
        }
        for token, hpiv_type in path_map.items():
            if token in lowered:
                assets = _resolve_hpiv_reference_assets(hpiv_type)
                return {
                    "status": "ready",
                    "hpiv_type": hpiv_type,
                    "species_label": _hpiv_species_label(hpiv_type),
                    "reference_path": str(Path(requested).resolve()),
                    "gff_path": str(assets["annotation_gff"]) if assets["annotation_gff"].is_file() else "nogtf",
                    "summary_path": "",
                }

    screening_dir = Path(f"{pre}_hpiv_reference_selection")
    screening_dir.mkdir(exist_ok=True)
    results: list[dict[str, object]] = []
    for hpiv_type in ["1", "2", "3", "4A", "4B"]:
        assets = _resolve_hpiv_reference_assets(hpiv_type)
        reference_fasta = assets["reference_fasta"]
        if not reference_fasta.is_file():
            continue
        type_dir = screening_dir / f"hpiv{str(hpiv_type).lower()}"
        type_dir.mkdir(parents=True, exist_ok=True)
        bam_path = type_dir / "screening.bam"
        coverage_path = type_dir / "coverage.tsv"
        if not bam_path.is_file():
            if fq1:
                map_cmd = " ".join(
                    [
                        "minimap2",
                        "-ax",
                        "sr",
                        shlex.quote(str(reference_fasta)),
                        shlex.quote(str(fq1)),
                        *([shlex.quote(str(fq2))] if fq2 else []),
                        "|",
                        "samtools",
                        "sort",
                        "-o",
                        shlex.quote(str(bam_path)),
                    ]
                )
            else:
                map_cmd = " ".join(
                    [
                        "minimap2",
                        "-ax",
                        _choose_minimap2_preset(long_type),
                        shlex.quote(str(reference_fasta)),
                        shlex.quote(str(single_fastq)),
                        "-t",
                        str(max(1, int(threads or 1))),
                        "|",
                        "samtools",
                        "sort",
                        "-o",
                        shlex.quote(str(bam_path)),
                    ]
                )
            run_command(map_cmd, logf=logf)
        if bam_path.is_file():
            run_command(f"samtools index -@ {max(1, int(threads or 1))} {shlex.quote(str(bam_path))}", logf=logf)
            run_command(f"samtools coverage {shlex.quote(str(bam_path))} > {shlex.quote(str(coverage_path))}", logf=logf)
        coverage_info = _parse_samtools_coverage_tsv(coverage_path)
        results.append(
            {
                "hpiv_type": hpiv_type,
                "species_label": _hpiv_species_label(hpiv_type),
                "reference_path": str(reference_fasta),
                "gff_path": str(assets["annotation_gff"]) if assets["annotation_gff"].is_file() else "nogtf",
                "coverage": coverage_info["coverage"],
                "mean_depth": coverage_info["mean_depth"],
                "covered_bases": coverage_info["covered_bases"],
                "num_reads": coverage_info["num_reads"],
            }
        )
    summary_path = screening_dir / "selection.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "hpiv_type", "coverage", "mean_depth", "covered_bases", "num_reads", "reference_path", "gff_path"])
        for item in results:
            writer.writerow([
                pre,
                item["hpiv_type"],
                f"{float(item['coverage']):.6f}",
                f"{float(item['mean_depth']):.6f}",
                f"{float(item['covered_bases']):.0f}",
                f"{float(item['num_reads']):.0f}",
                item["reference_path"],
                item["gff_path"],
            ])
    if not results:
        return {"status": "missing", "hpiv_type": "-", "species_label": species or "Human parainfluenza virus", "reference_path": "", "gff_path": "nogtf", "summary_path": str(summary_path)}
    results.sort(key=lambda item: (float(item["coverage"]), float(item["mean_depth"]), float(item["covered_bases"]), float(item["num_reads"])), reverse=True)
    best = results[0]
    return {
        "status": "ready",
        "hpiv_type": str(best["hpiv_type"]),
        "species_label": str(best["species_label"]),
        "reference_path": str(best["reference_path"]),
        "gff_path": str(best["gff_path"]),
        "summary_path": str(summary_path.resolve()),
    }


def _run_hadv_gene_blast_typing(query_fasta: Path, gene_name: str, db_fasta: Path, out_dir: Path, logf=None) -> dict[str, str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    out_tsv = out_dir / f"{gene_name}.blastn.tsv"
    blastn_bin = _resolve_hadv_blastn_bin()
    cmd = " ".join(
        [
            shlex.quote(blastn_bin),
            "-task blastn",
            "-query",
            shlex.quote(str(query_fasta)),
            "-db",
            shlex.quote(str(db_fasta)),
            "-outfmt",
            shlex.quote("6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore"),
            "-evalue 1e-20",
            "-max_target_seqs 10",
            "-dust no",
            "-num_threads 4",
            "-out",
            shlex.quote(str(out_tsv)),
        ]
    )
    run_command(cmd, logf=logf)
    best: dict[str, str] = {}
    if not out_tsv.is_file() or out_tsv.stat().st_size == 0:
        return {"gene": gene_name, "type": "", "subject": "", "identity": "", "coverage": "", "bitscore": "", "method": "blastn"}
    with out_tsv.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            qseqid, sseqid, pident, length, qlen, slen, qstart, qend, sstart, send, evalue, bitscore = parts[:12]
            try:
                coverage = 100.0 * int(length) / max(1, int(slen))
                score = (float(bitscore), float(pident), coverage, int(length))
            except (TypeError, ValueError):
                continue
            previous_score = best.get("_score")
            if previous_score is None or score > previous_score:
                best = {
                    "gene": gene_name,
                    "type": _parse_hadv_gene_type(sseqid),
                    "subject": sseqid,
                    "identity": pident,
                    "coverage": f"{coverage:.2f}",
                    "bitscore": bitscore,
                    "method": "blastn",
                    "_score": score,
                }
    best.pop("_score", None)
    return best or {"gene": gene_name, "type": "", "subject": "", "identity": "", "coverage": "", "bitscore": "", "method": "blastn"}


def _run_hadv_gene_read_typing(
    gene_name: str,
    db_fasta: Path,
    out_dir: Path,
    single_fastq: str = "",
    fq1: str = "",
    fq2: str = "",
    long_type: str = "",
    threads: int = 4,
    logf=None,
) -> dict[str, str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    bam_path = out_dir / f"{gene_name}.screening.bam"
    coverage_path = out_dir / f"{gene_name}.coverage.tsv"
    if not bam_path.is_file():
        if fq1:
            map_cmd = " ".join(
                [
                    "minimap2",
                    "-ax",
                    "sr",
                    shlex.quote(str(db_fasta)),
                    shlex.quote(str(fq1)),
                    *([shlex.quote(str(fq2))] if fq2 else []),
                    "|",
                    "samtools",
                    "sort",
                    "-o",
                    shlex.quote(str(bam_path)),
                ]
            )
        else:
            map_cmd = " ".join(
                [
                    "minimap2",
                    "-ax",
                    _choose_minimap2_preset(long_type),
                    shlex.quote(str(db_fasta)),
                    shlex.quote(str(single_fastq)),
                    "-t",
                    str(max(1, int(threads or 1))),
                    "|",
                    "samtools",
                    "sort",
                    "-o",
                    shlex.quote(str(bam_path)),
                ]
            )
        run_command(map_cmd, logf=logf)
    if not bam_path.is_file():
        return {"gene": gene_name, "type": "", "subject": "", "identity": "", "coverage": "", "bitscore": "", "method": "read_coverage"}
    run_command(f"samtools index -@ {max(1, int(threads or 1))} {shlex.quote(str(bam_path))}", logf=logf)
    run_command(f"samtools coverage {shlex.quote(str(bam_path))} > {shlex.quote(str(coverage_path))}", logf=logf)
    rows = _parse_samtools_coverage_rows(coverage_path)
    if not rows:
        return {"gene": gene_name, "type": "", "subject": "", "identity": "", "coverage": "", "bitscore": "", "method": "read_coverage"}
    rows.sort(
        key=lambda item: (
            float(item.get("coverage") or 0.0),
            float(item.get("mean_depth") or 0.0),
            float(item.get("covered_bases") or 0.0),
            float(item.get("num_reads") or 0.0),
        ),
        reverse=True,
    )
    best = rows[0]
    subject = str(best.get("reference_name") or "")
    return {
        "gene": gene_name,
        "type": _parse_hadv_gene_type(subject),
        "subject": subject,
        "identity": "",
        "coverage": f"{float(best.get('coverage') or 0.0):.2f}",
        "bitscore": "",
        "method": "read_coverage",
    }


def _run_hadv_phf_typing(
    pre: str,
    query_fasta: Path | None = None,
    single_fastq: str = "",
    fq1: str = "",
    fq2: str = "",
    long_type: str = "",
    threads: int = 4,
    logf=None,
) -> dict[str, object]:
    root_dir = _resolve_hadv_db_dir()
    typing_dir = Path(f"{pre}_hadv_reference_selection") / "phf_typing"
    gene_dbs = {
        "penton": root_dir / "blastn_db_penton" / "hadv_types_ref_penton.fa",
        "hexon": root_dir / "blastn_db_hexon" / "hadv_types_ref_hexon.fa",
        "fiber": root_dir / "blastn_db_fiber" / "hadv_types_ref_fiber.fa",
    }
    hits: dict[str, dict[str, str]] = {}
    use_blast = query_fasta is not None and query_fasta.is_file() and query_fasta.stat().st_size > 0
    for gene_name, db_fasta in gene_dbs.items():
        if not db_fasta.is_file():
            hits[gene_name] = {"gene": gene_name, "type": "", "subject": "", "identity": "", "coverage": "", "bitscore": "", "method": "missing_db"}
            continue
        if use_blast:
            hits[gene_name] = _run_hadv_gene_blast_typing(query_fasta, gene_name, db_fasta.resolve(), typing_dir, logf=logf)
        else:
            hits[gene_name] = _run_hadv_gene_read_typing(
                gene_name,
                db_fasta.resolve(),
                typing_dir,
                single_fastq=single_fastq,
                fq1=fq1,
                fq2=fq2,
                long_type=long_type,
                threads=threads,
                logf=logf,
            )
    penton_type = hits.get("penton", {}).get("type", "")
    hexon_type = hits.get("hexon", {}).get("type", "")
    fiber_type = hits.get("fiber", {}).get("type", "")
    hadv_type = _infer_hadv_type_from_phf(penton_type, hexon_type, fiber_type)
    summary_path = typing_dir / "phf_typing.tsv"
    typing_dir.mkdir(parents=True, exist_ok=True)
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "gene", "matched_type", "subject", "identity", "coverage", "bitscore", "method"])
        for gene_name in ["penton", "hexon", "fiber"]:
            hit = hits.get(gene_name, {})
            writer.writerow([pre, gene_name, hit.get("type", ""), hit.get("subject", ""), hit.get("identity", ""), hit.get("coverage", ""), hit.get("bitscore", ""), hit.get("method", "")])
    return {
        "penton_type": penton_type,
        "hexon_type": hexon_type,
        "fiber_type": fiber_type,
        "hadv_type": hadv_type,
        "summary_path": str(summary_path.resolve()),
        "hits": hits,
    }


def resolve_hadv_reference(
    pre: str,
    species: str = "",
    requested_ref: str = "",
    single_fastq: str = "",
    fq1: str = "",
    fq2: str = "",
    long_type: str = "",
    threads: int = 4,
    query_fasta: str = "",
    logf=None,
) -> dict[str, str]:
    requested = str(requested_ref or "").strip()
    if requested and Path(requested).is_file():
        return {
            "status": "ready",
            "hadv_type": _extract_hadv_type_label(requested),
            "species_label": species or "Human adenovirus",
            "reference_path": str(Path(requested).resolve()),
            "gff_path": "nogtf",
            "summary_path": "",
            "phf_summary_path": "",
        }
    query_path = Path(query_fasta) if str(query_fasta or "").strip() else None
    phf_result = _run_hadv_phf_typing(
        pre,
        query_fasta=query_path,
        single_fastq=single_fastq,
        fq1=fq1,
        fq2=fq2,
        long_type=long_type,
        threads=threads,
        logf=logf,
    )
    hadv_type = str(phf_result.get("hadv_type") or "").strip()
    if not hadv_type:
        hadv_type = _extract_hadv_type_label(species)
    screening_dir = Path(f"{pre}_hadv_reference_selection")
    screening_dir.mkdir(exist_ok=True)
    summary_path = screening_dir / "selection.tsv"
    if not hadv_type:
        with summary_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sample", "hadv_type", "penton_type", "hexon_type", "fiber_type", "status", "note"])
            writer.writerow([pre, "", phf_result.get("penton_type", ""), phf_result.get("hexon_type", ""), phf_result.get("fiber_type", ""), "missing", f"未能根据 PHF 组合或 species={species or '-'} 确定 HAdV 分型"])
        return {"status": "missing", "hadv_type": "", "species_label": species or "Human adenovirus", "reference_path": "", "gff_path": "nogtf", "summary_path": str(summary_path), "phf_summary_path": str(phf_result.get("summary_path") or "")}

    candidates = _read_hadv_full_genome_candidates(hadv_type)
    candidate_fasta = screening_dir / "candidate_references.fasta"
    if candidates:
        _write_hadv_candidate_fasta(candidates, candidate_fasta)
    else:
        with summary_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sample", "hadv_type", "penton_type", "hexon_type", "fiber_type", "status", "note"])
            writer.writerow([pre, hadv_type, phf_result.get("penton_type", ""), phf_result.get("hexon_type", ""), phf_result.get("fiber_type", ""), "missing", "未找到对应 HAdV 分型候选参考库"])
        return {"status": "missing", "hadv_type": hadv_type, "species_label": _hadv_species_label(hadv_type), "reference_path": "", "gff_path": "nogtf", "summary_path": str(summary_path), "phf_summary_path": str(phf_result.get("summary_path") or "")}

    bam_path = screening_dir / "candidate.screening.bam"
    coverage_path = screening_dir / "candidate.coverage.tsv"
    if not bam_path.is_file():
        if fq1:
            map_cmd = " ".join(
                [
                    "minimap2",
                    "-ax",
                    "sr",
                    shlex.quote(str(candidate_fasta)),
                    shlex.quote(str(fq1)),
                    *([shlex.quote(str(fq2))] if fq2 else []),
                    "|",
                    "samtools",
                    "sort",
                    "-o",
                    shlex.quote(str(bam_path)),
                ]
            )
        else:
            map_cmd = " ".join(
                [
                    "minimap2",
                    "-ax",
                    _choose_minimap2_preset(long_type),
                    shlex.quote(str(candidate_fasta)),
                    shlex.quote(str(single_fastq)),
                    "-t",
                    str(max(1, int(threads or 1))),
                    "|",
                    "samtools",
                    "sort",
                    "-o",
                    shlex.quote(str(bam_path)),
                ]
            )
        run_command(map_cmd, logf=logf)
    if bam_path.is_file():
        run_command(f"samtools index -@ {max(1, int(threads or 1))} {shlex.quote(str(bam_path))}", logf=logf)
        run_command(f"samtools coverage {shlex.quote(str(bam_path))} > {shlex.quote(str(coverage_path))}", logf=logf)
    rows = _parse_samtools_coverage_rows(coverage_path)
    rows.sort(
        key=lambda item: (
            float(item.get("coverage") or 0.0),
            float(item.get("mean_depth") or 0.0),
            float(item.get("covered_bases") or 0.0),
            float(item.get("num_reads") or 0.0),
        ),
        reverse=True,
    )
    if not rows:
        return {"status": "missing", "hadv_type": hadv_type, "species_label": _hadv_species_label(hadv_type), "reference_path": "", "gff_path": "nogtf", "summary_path": str(summary_path), "phf_summary_path": str(phf_result.get("summary_path") or "")}
    best = rows[0]
    best_reference_name = str(best.get("reference_name") or "").strip()
    best_reference_fasta = screening_dir / "hadv.best_reference.fasta"
    if not _extract_named_fasta_record(candidate_fasta, best_reference_name, best_reference_fasta):
        shutil.copy2(candidate_fasta, best_reference_fasta)
    accession = best_reference_name.split()[0].split(".")[0]
    gff_path = _resolve_hadv_gff_path(accession)
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "hadv_type", "penton_type", "hexon_type", "fiber_type", "coverage", "mean_depth", "covered_bases", "num_reads", "reference_name", "reference_path", "gff_path", "phf_summary_path"])
        writer.writerow([
            pre,
            hadv_type,
            phf_result.get("penton_type", ""),
            phf_result.get("hexon_type", ""),
            phf_result.get("fiber_type", ""),
            f"{float(best.get('coverage') or 0.0):.6f}",
            f"{float(best.get('mean_depth') or 0.0):.6f}",
            f"{float(best.get('covered_bases') or 0.0):.0f}",
            f"{float(best.get('num_reads') or 0.0):.0f}",
            best_reference_name,
            str(best_reference_fasta.resolve()),
            str(gff_path) if gff_path.is_file() else "nogtf",
            str(phf_result.get("summary_path") or ""),
        ])
    return {
        "status": "ready",
        "hadv_type": hadv_type,
        "species_label": _hadv_species_label(hadv_type),
        "reference_path": str(best_reference_fasta.resolve()),
        "gff_path": str(gff_path) if gff_path.is_file() else "nogtf",
        "summary_path": str(summary_path.resolve()),
        "phf_summary_path": str(phf_result.get("summary_path") or ""),
    }


def build_hpiv_coverage_assets(
    pre: str,
    species: str = "",
    single_fastq: str = "",
    fq1: str = "",
    fq2: str = "",
    long_type: str = "",
    threads: int = 4,
    logf=None,
) -> dict[str, object]:
    hpiv_type = _infer_hpiv_type_from_label(species)
    selection_path = Path(f"{pre}_hpiv_reference_selection") / "selection.tsv"
    if hpiv_type not in {"1", "2", "3", "4A", "4B"} and selection_path.is_file():
        with selection_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            rows = list(reader)
        rows.sort(key=lambda item: (float(item.get("coverage") or 0.0), float(item.get("mean_depth") or 0.0)), reverse=True)
        if rows:
            hpiv_type = str(rows[0].get("hpiv_type") or "").strip().upper()
    if hpiv_type not in {"1", "2", "3", "4A", "4B"}:
        return {"status": "missing", "note": "未能确定 HPIV 分型，无法构建覆盖度专用参考。"}

    db_fasta = _resolve_hpiv_subtype_db_fasta(hpiv_type)
    if not db_fasta.is_file():
        return {"status": "missing", "note": f"未找到 HPIV{hpiv_type} 的 _db.fasta 参考集合。"}

    source_fastq = str(fq1 or single_fastq or "").strip()
    if not source_fastq:
        return {"status": "missing", "note": "未找到可用于 HPIV 覆盖度重比对的 FASTQ 输入。"}

    coverage_root = Path("hpiv_coverage")
    coverage_root.mkdir(parents=True, exist_ok=True)
    multi_ref_bam = coverage_root / "screening.bam"
    multi_ref_coverage = coverage_root / "screening.coverage.tsv"
    if not multi_ref_bam.is_file():
        if fq1:
            map_cmd = " ".join(
                [
                    "minimap2",
                    "-ax",
                    "sr",
                    shlex.quote(str(db_fasta)),
                    shlex.quote(str(fq1)),
                    *([shlex.quote(str(fq2))] if fq2 else []),
                    "|",
                    "samtools",
                    "sort",
                    "-o",
                    shlex.quote(str(multi_ref_bam)),
                ]
            )
        else:
            map_cmd = " ".join(
                [
                    "minimap2",
                    "-ax",
                    _choose_minimap2_preset(long_type),
                    shlex.quote(str(db_fasta)),
                    shlex.quote(str(single_fastq)),
                    "-t",
                    str(max(1, int(threads or 1))),
                    "|",
                    "samtools",
                    "sort",
                    "-o",
                    shlex.quote(str(multi_ref_bam)),
                ]
            )
        run_command(map_cmd, logf=logf)
    if not multi_ref_bam.is_file():
        return {"status": "failed", "note": "HPIV 覆盖度参考初筛比对失败。"}

    run_command(f"samtools index -@ {max(1, int(threads or 1))} {shlex.quote(str(multi_ref_bam))}", logf=logf)
    run_command(f"samtools coverage {shlex.quote(str(multi_ref_bam))} > {shlex.quote(str(multi_ref_coverage))}", logf=logf)
    coverage_rows = _parse_samtools_coverage_rows(multi_ref_coverage)
    if not coverage_rows:
        return {"status": "failed", "note": "未获得 HPIV 覆盖度参考初筛结果。"}
    coverage_rows.sort(
        key=lambda item: (
            float(item.get("coverage") or 0.0),
            float(item.get("mean_depth") or 0.0),
            float(item.get("covered_bases") or 0.0),
            float(item.get("num_reads") or 0.0),
        ),
        reverse=True,
    )
    best = coverage_rows[0]
    best_reference_name = str(best.get("reference_name") or "").strip()
    best_reference_fasta = coverage_root / "hpiv.coverage.ref.fa"
    if not _extract_named_fasta_record(db_fasta, best_reference_name, best_reference_fasta):
        return {"status": "failed", "note": f"无法从 {db_fasta.name} 中提取覆盖度最优参考：{best_reference_name}。"}

    run_command(f"samtools faidx {shlex.quote(str(best_reference_fasta))}", logf=logf)
    coverage_bam = coverage_root / "hpiv.coverage.mapping.bam"
    coverage_bed = coverage_root / "hpiv.coverage.regions.bed"
    coverage_depth = coverage_root / "hpiv.coverage.depth.tsv"
    if fq1:
        remap_cmd = " ".join(
            [
                "minimap2",
                "-ax",
                "sr",
                shlex.quote(str(best_reference_fasta)),
                shlex.quote(str(fq1)),
                *([shlex.quote(str(fq2))] if fq2 else []),
                "|",
                "samtools",
                "sort",
                "-o",
                shlex.quote(str(coverage_bam)),
            ]
        )
    else:
        remap_cmd = " ".join(
            [
                "minimap2",
                "-ax",
                _choose_minimap2_preset(long_type),
                shlex.quote(str(best_reference_fasta)),
                shlex.quote(str(single_fastq)),
                "-t",
                str(max(1, int(threads or 1))),
                "|",
                "samtools",
                "sort",
                "-o",
                shlex.quote(str(coverage_bam)),
            ]
        )
    run_command(remap_cmd, logf=logf)
    if not coverage_bam.is_file():
        return {"status": "failed", "note": "HPIV 覆盖度专用参考重比对失败。"}
    run_command(f"samtools index -@ {max(1, int(threads or 1))} {shlex.quote(str(coverage_bam))}", logf=logf)
    run_command(f"samtools depth -aa {shlex.quote(str(coverage_bam))} > {shlex.quote(str(coverage_depth))}", logf=logf)
    if not coverage_depth.is_file() or coverage_depth.stat().st_size == 0:
        return {"status": "failed", "note": "未生成 HPIV 覆盖度深度文件。"}

    with coverage_depth.open("r", encoding="utf-8", errors="ignore") as in_handle, coverage_bed.open("w", encoding="utf-8") as out_handle:
        for line in in_handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            try:
                end_pos = int(parts[1])
            except (TypeError, ValueError):
                continue
            out_handle.write(f"{parts[0]}\t{max(0, end_pos - 1)}\t{end_pos}\t{parts[2]}\n")

    summary_path = coverage_root / "hpiv.coverage.summary.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "sample",
                "hpiv_type",
                "coverage_reference_name",
                "coverage_reference_path",
                "coverage",
                "mean_depth",
                "covered_bases",
                "num_reads",
            ]
        )
        writer.writerow(
            [
                pre,
                hpiv_type,
                best_reference_name,
                str(best_reference_fasta.resolve()),
                f"{float(best.get('coverage') or 0.0):.6f}",
                f"{float(best.get('mean_depth') or 0.0):.6f}",
                f"{float(best.get('covered_bases') or 0.0):.0f}",
                f"{float(best.get('num_reads') or 0.0):.0f}",
            ]
        )
    return {
        "status": "ready",
        "hpiv_type": hpiv_type,
        "coverage_reference_name": best_reference_name,
        "coverage_reference_path": str(best_reference_fasta.resolve()),
        "coverage_bam_path": str(coverage_bam.resolve()),
        "coverage_bed_path": str(coverage_bed.resolve()),
        "summary_path": str(summary_path.resolve()),
        "coverage": float(best.get("coverage") or 0.0),
        "mean_depth": float(best.get("mean_depth") or 0.0),
    }


def _nextclade_type_candidates(species: str) -> list[str]:
    inferred = _infer_influenza_type_from_label(species)
    if inferred == "Influenza A virus":
        return ["A"]
    if inferred == "Influenza B virus":
        return ["B"]
    if inferred == "Influenza C virus":
        return ["C"]
    if inferred == "Influenza D virus":
        return ["D"]
    return ["A", "B", "C", "D"]


def _parse_nextclade_result(tsv_path: Path) -> dict[str, str]:
    if not tsv_path.is_file() or tsv_path.stat().st_size == 0:
        return {}
    with tsv_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        first = next(reader, None) or {}
    return {str(key): str(value or "").strip() for key, value in first.items()}


def _nextclade_status_rank(value: str) -> int:
    mapping = {"good": 0, "mediocre": 1, "bad": 2}
    return mapping.get(str(value or "").strip().lower(), 3)


def _nextclade_score(value: str) -> float:
    try:
        return float(str(value or "").strip())
    except (TypeError, ValueError):
        return 10**9


def _choose_best_nextclade_result(results: list[dict[str, object]]) -> dict[str, object] | None:
    valid = [item for item in results if item.get("row")]
    if not valid:
        return None
    valid.sort(
        key=lambda item: (
            1 if str(item["row"].get("errors") or "").strip() else 0,
            _nextclade_status_rank(str(item["row"].get("qc.overallStatus") or "")),
            _nextclade_score(str(item["row"].get("qc.overallScore") or "")),
            str(item["row"].get("clade") or "") == "",
        )
    )
    return valid[0]


def _write_nextclade_taxonomy(pre: str, species_label: str, clade: str) -> None:
    species_path = Path(f"{pre}_2.list.txt")
    subspecies_path = Path(f"{pre}_2.list2.txt")
    with species_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["界", "属", "种", "比例", "序列数量"])
        writer.writerow(["病毒", "Influenza", species_label, "100", "1"])
    with subspecies_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["界", "属", "种", "亚种", "比例", "序列数量"])
        writer.writerow(["病毒", "Influenza", species_label, clade or "-", "100", "1"])


def nextclade_identify(pre: str, species: str) -> dict[str, str]:
    columns = ["样本名称", "病毒类型", "Nextclade数据集", "Clade", "QC状态", "QC分数", "说明"]
    final_fasta = Path(f"{pre}.final.fasta")
    out_dir = Path(f"{pre}_nextclade")
    out_dir.mkdir(exist_ok=True)
    summary_path = out_dir / "nextclade_identification.tsv"
    if not final_fasta.is_file() or final_fasta.stat().st_size == 0:
        with summary_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(columns)
            writer.writerow([pre, species or "-", "-", "-", "-", "-", "未检测到可用于 Nextclade 鉴定的 FASTA"])
        return {"virus_type": species or "-", "clade": "-", "status": "empty"}

    if not _is_influenza(species):
        with summary_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(columns)
            writer.writerow([pre, species or "-", "-", "-", "-", "-", "当前仅对流感病毒 FASTA 输入启用 Nextclade 鉴定"])
        return {"virus_type": species or "-", "clade": "-", "status": "skipped"}

    nextclade_bin = _resolve_nextclade_binary()
    log_path = out_dir / "nextclade.log"
    results: list[dict[str, object]] = []
    with log_path.open("w", encoding="utf-8") as logf:
        for flu_type in _nextclade_type_candidates(species):
            dataset_dir = _resolve_nextclade_dataset(flu_type)
            if not dataset_dir.exists():
                continue
            tsv_path = out_dir / f"influenza_{flu_type.lower()}.tsv"
            cmd = " ".join(
                [
                    shlex.quote(nextclade_bin),
                    "run",
                    "--input-dataset",
                    shlex.quote(str(dataset_dir)),
                    "--output-tsv",
                    shlex.quote(str(tsv_path)),
                    shlex.quote(str(final_fasta)),
                ]
            )
            try:
                run_command(cmd, logf=logf)
            except Exception as exc:
                logf.write(f"[NEXTCLADE_FAIL] dataset={dataset_dir} error={exc}\n")
                logf.flush()
                continue
            results.append(
                {
                    "type": flu_type,
                    "dataset": str(dataset_dir),
                    "tsv_path": str(tsv_path),
                    "row": _parse_nextclade_result(tsv_path),
                }
            )

    best = _choose_best_nextclade_result(results)
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(columns)
        if not best:
            writer.writerow([pre, species or "-", "-", "-", "-", "-", "未获得可用的 Nextclade 结果，请检查 nextclade 软件和数据集目录"])
            return {"virus_type": species or "-", "clade": "-", "status": "failed"}

        row = dict(best["row"])
        flu_type = str(best["type"])
        species_label = INFLUENZA_TYPE_TO_LABEL.get(flu_type, species or "-")
        clade = str(row.get("clade") or row.get("Nextclade_clade") or "").strip() or "-"
        qc_status = str(row.get("qc.overallStatus") or "").strip() or "-"
        qc_score = str(row.get("qc.overallScore") or "").strip() or "-"
        message = "Nextclade FASTA 鉴定完成"
        if str(row.get("errors") or "").strip():
            message = f"Nextclade 返回错误提示：{row['errors']}"
        writer.writerow([pre, species_label, str(best["dataset"]), clade, qc_status, qc_score, message])
        _write_nextclade_taxonomy(pre, species_label, clade)
        update_runtime_context(species=species_label)
        return {"virus_type": species_label, "clade": clade, "status": "ready"}


def _ensure_blast_db(source_path: Path, prefix_path: Path, logf=None) -> Path:
    nin_path = prefix_path.with_suffix(".nin")
    nhr_path = prefix_path.with_suffix(".nhr")
    nsq_path = prefix_path.with_suffix(".nsq")
    if nin_path.is_file() and nhr_path.is_file() and nsq_path.is_file():
        return prefix_path
    run_command(
        f"makeblastdb -in {shlex.quote(str(source_path))} -dbtype nucl -out {shlex.quote(str(prefix_path))}",
        logf=logf,
    )
    return prefix_path


def _extract_subtype(text: str, prefix: str) -> str:
    match = re.search(rf"{prefix}\s*([0-9]+)", str(text or ""), re.IGNORECASE)
    if match:
        return f"{prefix.upper()}{match.group(1)}"
    compact = re.search(rf"({prefix.upper()}[0-9]+)", str(text or ""), re.IGNORECASE)
    return compact.group(1).upper() if compact else "-"


def _blast_top_hit(query_fasta: Path, db_prefix: Path, out_path: Path, subtype_prefix: str, logf=None) -> dict:
    run_command(
        " ".join(
            [
                "blastn",
                "-query", shlex.quote(str(query_fasta)),
                "-db", shlex.quote(str(db_prefix)),
                "-outfmt", shlex.quote("6 qseqid sseqid pident length evalue bitscore qcovs"),
                "-max_target_seqs", "1",
                "-max_hsps", "1",
                "-out", shlex.quote(str(out_path)),
            ]
        ),
        logf=logf,
    )
    if not out_path.is_file() or out_path.stat().st_size == 0:
        return {
            "subtype": "-",
            "subject": "-",
            "identity": "-",
            "coverage": "-",
        }
    with out_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        first = next(reader, None)
    if not first or len(first) < 7:
        return {
            "subtype": "-",
            "subject": "-",
            "identity": "-",
            "coverage": "-",
        }
    subject = str(first[1] or "").strip()
    return {
        "subtype": _extract_subtype(subject, subtype_prefix),
        "subject": subject or "-",
        "identity": f"{float(first[2]):.2f}%" if first[2] else "-",
        "coverage": f"{float(first[6]):.2f}%" if first[6] else "-",
    }


def _blast_type_hit(query_fasta: Path, db_prefix: Path, out_path: Path, logf=None) -> dict:
    run_command(
        " ".join(
            [
                "blastn",
                "-query", shlex.quote(str(query_fasta)),
                "-db", shlex.quote(str(db_prefix)),
                "-outfmt", shlex.quote("6 qseqid sseqid pident length evalue bitscore qcovs"),
                "-max_target_seqs", "1",
                "-max_hsps", "1",
                "-out", shlex.quote(str(out_path)),
            ]
        ),
        logf=logf,
    )
    if not out_path.is_file() or out_path.stat().st_size == 0:
        return {"virus_type": "-", "subject": "-", "identity": "-", "coverage": "-"}
    with out_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        first = next(reader, None)
    if not first or len(first) < 7:
        return {"virus_type": "-", "subject": "-", "identity": "-", "coverage": "-"}
    subject = str(first[1] or "").strip()
    subject_lower = subject.lower()
    if "influenza_a" in subject_lower or "influenza a" in subject_lower or "甲流" in subject:
        virus_type = "Influenza A virus"
    elif "influenza_b" in subject_lower or "influenza b" in subject_lower or "乙流" in subject:
        virus_type = "Influenza B virus"
    else:
        virus_type = "-"
    return {
        "virus_type": virus_type,
        "subject": subject or "-",
        "identity": f"{float(first[2]):.2f}%" if first[2] else "-",
        "coverage": f"{float(first[6]):.2f}%" if first[6] else "-",
    }


def _write_placeholder(pre: str, columns: list[str], row: list[str]) -> None:
    out_path = Path(f"{pre}_serotype_result.tsv")
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(columns)
        writer.writerow(row)


def _resolve_sars_cov_2_nextclade_input_fasta(pre: str, final_fasta: Path) -> Path:
    preferred = Path("snps.raw.high_quality.consensus.fasta")
    if preferred.is_file() and preferred.stat().st_size > 0:
        return preferred
    legacy_preferred = Path(f"{pre}.high_quality.consensus.fasta")
    if legacy_preferred.is_file() and legacy_preferred.stat().st_size > 0:
        return legacy_preferred
    return final_fasta


def _resolve_monkeypox_nextclade_input_fasta(pre: str, final_fasta: Path) -> Path:
    preferred = Path(f"{pre}.consensus.fasta")
    if preferred.is_file() and preferred.stat().st_size > 0:
        return preferred
    return final_fasta


def _resolve_rsv_nextclade_input_fasta(pre: str, final_fasta: Path) -> Path:
    preferred = Path(f"{pre}.consensus.fasta")
    if preferred.is_file() and preferred.stat().st_size > 0:
        return preferred
    return final_fasta


def _is_hmpv(species: str) -> bool:
    normalized = str(species or "").strip().lower()
    return (
        normalized in {"human metapneumovirus", "metapneumovirus", "人偏肺病毒", "hmpv"}
        or "metapneumovirus" in normalized
        or "偏肺病毒" in normalized
    )


def _resolve_hmpv_nextclade_input_fasta(pre: str, final_fasta: Path) -> Path:
    preferred = Path(f"{pre}.consensus.fasta")
    if preferred.is_file() and preferred.stat().st_size > 0:
        return preferred
    return final_fasta


def _resolve_denv_nextclade_input_fasta(pre: str, final_fasta: Path) -> Path:
    preferred = Path(f"{pre}.consensus.fasta")
    if preferred.is_file() and preferred.stat().st_size > 0:
        return preferred
    return final_fasta


def _resolve_zika_nextclade_input_fasta(pre: str, final_fasta: Path) -> Path:
    preferred = Path(f"{pre}.consensus.fasta")
    if preferred.is_file() and preferred.stat().st_size > 0:
        return preferred
    return final_fasta


def _resolve_chikv_nextclade_input_fasta(pre: str, final_fasta: Path) -> Path:
    preferred = Path(f"{pre}.consensus.fasta")
    if preferred.is_file() and preferred.stat().st_size > 0:
        return preferred
    return final_fasta


def _resolve_ebola_nextclade_input_fasta(pre: str, final_fasta: Path) -> Path:
    preferred = Path(f"{pre}.consensus.fasta")
    if preferred.is_file() and preferred.stat().st_size > 0:
        return preferred
    return final_fasta


def _run_generic_nextclade_typing(
    pre: str,
    virus_label: str,
    dataset_dir: Path,
    input_fasta: Path,
    clade_keys: list[str],
    lineage_keys: list[str],
) -> bool:
    columns = ["样本名称", "病毒类型", "Nextclade分型", "Pango谱系", "Nextclade数据集", "QC状态", "QC分数", "说明"]
    if not dataset_dir.exists():
        _write_placeholder(pre, columns, [pre, virus_label, "-", "-", str(dataset_dir), "-", "-", f"未找到 {virus_label} Nextclade 数据集"])
        return True

    out_dir = Path("nextclade_output")
    out_dir.mkdir(exist_ok=True)
    compat_out_dir = Path(f"{pre}_nextclade")
    compat_out_dir.mkdir(exist_ok=True)
    tsv_path = out_dir / "nextclade.tsv"
    compat_tsv_path = compat_out_dir / "nextclade.tsv"
    log_path = out_dir / "nextclade.log"
    nextclade_bin = _resolve_nextclade_binary()
    with log_path.open("a", encoding="utf-8") as logf:
        cmd = " ".join(
            [
                shlex.quote(nextclade_bin),
                "run",
                "--input-dataset",
                shlex.quote(str(dataset_dir)),
                "--output-tsv",
                shlex.quote(str(tsv_path)),
                shlex.quote(str(input_fasta)),
            ]
        )
        try:
            run_command(cmd, logf=logf)
        except Exception as exc:
            _write_placeholder(pre, columns, [pre, virus_label, "-", "-", str(dataset_dir), "-", "-", f"Nextclade 运行失败: {exc}"])
            return True
    if tsv_path.is_file():
        try:
            shutil.copy2(tsv_path, compat_tsv_path)
        except OSError:
            pass

    row = _parse_nextclade_result(tsv_path)
    if not row:
        _write_placeholder(pre, columns, [pre, virus_label, "-", "-", str(dataset_dir), "-", "-", "未获得可用的 Nextclade 分型结果"])
        return True

    clade = "-"
    for key in clade_keys:
        value = str(row.get(key) or "").strip()
        if value:
            clade = value
            break
    lineage = "-"
    for key in lineage_keys:
        value = str(row.get(key) or "").strip()
        if value:
            lineage = value
            break
    qc_status = str(row.get("qc.overallStatus") or "").strip() or "-"
    qc_score = str(row.get("qc.overallScore") or "").strip() or "-"
    message = f"{virus_label} Nextclade 分型完成（输入序列：{input_fasta.name}）"
    if str(row.get("errors") or "").strip():
        message = f"Nextclade 返回错误提示：{row['errors']}"
    _write_placeholder(pre, columns, [pre, virus_label, clade, lineage, str(dataset_dir), qc_status, qc_score, message])
    return True


def _run_sars_cov_2_nextclade_typing(pre: str, species: str, final_fasta: Path) -> bool:
    dataset_dir = _resolve_sars_cov_2_nextclade_dataset()
    input_fasta = _resolve_sars_cov_2_nextclade_input_fasta(pre, final_fasta)
    return _run_generic_nextclade_typing(
        pre,
        "SARS-CoV-2",
        dataset_dir,
        input_fasta,
        ["clade", "Nextclade_clade"],
        ["pangoLineage", "Nextclade_pango"],
    )


def _run_monkeypox_nextclade_typing(pre: str, species: str, final_fasta: Path) -> bool:
    dataset_dir = _resolve_monkeypox_nextclade_dataset()
    input_fasta = _resolve_monkeypox_nextclade_input_fasta(pre, final_fasta)
    return _run_generic_nextclade_typing(
        pre,
        species or "Monkeypox virus",
        dataset_dir,
        input_fasta,
        ["clade", "Nextclade_clade"],
        ["lineage", "Nextclade_pango"],
    )


def _run_rsv_nextclade_typing(pre: str, species: str, final_fasta: Path) -> bool:
    rsv_type = _infer_rsv_type_from_label(species)
    if rsv_type not in {"A", "B"}:
        selection_path = Path(f"{pre}_rsv_reference_selection") / "selection.tsv"
        if selection_path.is_file():
            with selection_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                rows = list(reader)
            rows.sort(key=lambda item: (float(item.get("coverage") or 0.0), float(item.get("mean_depth") or 0.0)), reverse=True)
            if rows:
                rsv_type = str(rows[0].get("rsv_type") or "").strip().upper()
    if rsv_type not in {"A", "B"}:
        _write_placeholder(pre, ["样本名称", "病毒类型", "Nextclade分型", "Pango谱系", "Nextclade数据集", "QC状态", "QC分数", "说明"], [pre, species or "Respiratory syncytial virus", "-", "-", "-", "-", "-", "未能确定 RSV A/B 参考类型，无法执行 Nextclade"])
        return True
    dataset_dir = _resolve_rsv_nextclade_dataset(rsv_type)
    input_fasta = _resolve_rsv_nextclade_input_fasta(pre, final_fasta)
    return _run_generic_nextclade_typing(
        pre,
        _rsv_species_label(rsv_type),
        dataset_dir,
        input_fasta,
        ["clade", "Nextclade_clade"],
        ["lineage", "genotype"],
    )


def _run_hmpv_nextclade_typing(pre: str, species: str, final_fasta: Path) -> bool:
    dataset_dir = _resolve_hmpv_nextclade_dataset()
    input_fasta = _resolve_hmpv_nextclade_input_fasta(pre, final_fasta)
    return _run_generic_nextclade_typing(
        pre,
        species or "Human metapneumovirus",
        dataset_dir,
        input_fasta,
        ["clade", "Nextclade_clade"],
        ["lineage", "genotype"],
    )


def _run_denv_nextclade_typing(pre: str, species: str, final_fasta: Path) -> bool:
    denv_type = _infer_denv_type_from_label(species)
    if denv_type not in {"1", "2", "3", "4"}:
        selection_path = Path(f"{pre}_denv_reference_selection") / "selection.tsv"
        if selection_path.is_file():
            with selection_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                rows = list(reader)
            rows.sort(key=lambda item: (float(item.get("coverage") or 0.0), float(item.get("mean_depth") or 0.0)), reverse=True)
            if rows:
                denv_type = str(rows[0].get("denv_type") or "").strip()
    if denv_type not in {"1", "2", "3", "4"}:
        _write_placeholder(pre, ["样本名称", "病毒类型", "Nextclade分型", "Pango谱系", "Nextclade数据集", "QC状态", "QC分数", "说明"], [pre, species or "Dengue virus", "-", "-", "-", "-", "-", "未能确定 DENV 参考型别，无法执行 Nextclade"])
        return True
    dataset_dir = _resolve_denv_nextclade_dataset(denv_type)
    input_fasta = _resolve_denv_nextclade_input_fasta(pre, final_fasta)
    return _run_generic_nextclade_typing(
        pre,
        _denv_species_label(denv_type),
        dataset_dir,
        input_fasta,
        ["clade", "Nextclade_clade"],
        ["lineage", "genotype"],
    )


def _run_zika_nextclade_typing(pre: str, species: str, final_fasta: Path) -> bool:
    dataset_dir = _resolve_zika_nextclade_dataset()
    input_fasta = _resolve_zika_nextclade_input_fasta(pre, final_fasta)
    return _run_generic_nextclade_typing(
        pre,
        species or "Zika virus",
        dataset_dir,
        input_fasta,
        ["clade", "Nextclade_clade"],
        ["lineage", "genotype"],
    )


def _run_chikv_nextclade_typing(pre: str, species: str, final_fasta: Path) -> bool:
    dataset_dir = _resolve_chikv_nextclade_dataset()
    input_fasta = _resolve_chikv_nextclade_input_fasta(pre, final_fasta)
    return _run_generic_nextclade_typing(
        pre,
        species or "Chikungunya virus",
        dataset_dir,
        input_fasta,
        ["clade", "Nextclade_clade"],
        ["lineage", "genotype"],
    )


def _run_ebola_nextclade_typing(pre: str, species: str, final_fasta: Path) -> bool:
    predicted_type = _infer_orthoebolavirus_abbrev_from_label(species)
    selection_path = Path(f"{pre}_orthoebolavirus_reference_selection") / "selection.tsv"
    if selection_path.is_file():
        try:
            with selection_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"), None) or {}
            predicted_type = str(row.get("predicted_type") or predicted_type).strip().upper()
        except OSError:
            pass
    if predicted_type and predicted_type != "EBOV":
        return False
    dataset_dir = _resolve_ebola_nextclade_dataset()
    input_fasta = _resolve_ebola_nextclade_input_fasta(pre, final_fasta)
    return _run_generic_nextclade_typing(
        pre,
        "Ebola virus",
        dataset_dir,
        input_fasta,
        ["clade", "Nextclade_clade"],
        ["lineage", "genotype"],
    )


def _run_hpiv_reference_typing(pre: str, species: str) -> bool:
    selection_path = Path(f"{pre}_hpiv_reference_selection") / "selection.tsv"
    coverage_summary_path = Path("hpiv_coverage") / "hpiv.coverage.summary.tsv"
    hpiv_type = _infer_hpiv_type_from_label(species)
    reference_path = ""
    gff_path = "nogtf"
    coverage = "-"
    mean_depth = "-"
    coverage_reference_name = ""
    if selection_path.is_file():
        with selection_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            rows = list(reader)
        rows.sort(key=lambda item: (float(item.get("coverage") or 0.0), float(item.get("mean_depth") or 0.0)), reverse=True)
        if rows:
            best = rows[0]
            if hpiv_type not in {"1", "2", "3", "4A", "4B"}:
                hpiv_type = str(best.get("hpiv_type") or "").strip().upper()
            reference_path = str(best.get("reference_path") or "").strip()
            gff_path = str(best.get("gff_path") or "").strip() or "nogtf"
            coverage = f"{float(best.get('coverage') or 0.0):.2f}%"
            mean_depth = f"{float(best.get('mean_depth') or 0.0):.2f}"
    if coverage_summary_path.is_file():
        with coverage_summary_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            coverage_row = next(reader, None) or {}
        coverage_reference_name = str(coverage_row.get("coverage_reference_name") or "").strip()
        if coverage_reference_name:
            coverage = f"{float(coverage_row.get('coverage') or 0.0):.2f}%"
            mean_depth = f"{float(coverage_row.get('mean_depth') or 0.0):.2f}"
    if hpiv_type not in {"1", "2", "3", "4A", "4B"}:
        _write_placeholder(
            pre,
            ["样本名称", "病毒类型", "HPIV亚型", "参考序列", "注释文件", "覆盖度", "平均深度", "说明"],
            [pre, species or "Human parainfluenza virus", "-", "-", "-", "-", "-", "未能确定 HPIV 最优参考型别"],
        )
        return True
    _write_placeholder(
        pre,
        ["样本名称", "病毒类型", "HPIV亚型", "参考序列", "注释文件", "覆盖度", "平均深度", "说明"],
        [
            pre,
            _hpiv_species_label(hpiv_type),
            hpiv_type,
            Path(reference_path).name if reference_path else "-",
            Path(gff_path).name if gff_path and gff_path != "nogtf" else "-",
            coverage,
            mean_depth,
            (
                f"覆盖度基于 {coverage_reference_name} 计算，SNP 仍沿用原参考序列。"
                if coverage_reference_name
                else "基于最优参考基因组覆盖度完成 HPIV 分型，不执行 Nextclade。"
            ),
        ],
    )
    return True


def _run_hadv_reference_typing(pre: str, species: str) -> bool:
    selection_path = Path(f"{pre}_hadv_reference_selection") / "selection.tsv"
    if not selection_path.is_file():
        _write_placeholder(
            pre,
            ["样本名称", "病毒类型", "HAdV分型", "Penton分型", "Hexon分型", "Fiber分型", "参考序列", "注释文件", "覆盖度", "平均深度", "说明"],
            [pre, species or "Human adenovirus", "-", "-", "-", "-", "-", "-", "-", "-", "未找到 HAdV 参考选择结果"],
        )
        return True
    with selection_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        _write_placeholder(
            pre,
            ["样本名称", "病毒类型", "HAdV分型", "Penton分型", "Hexon分型", "Fiber分型", "参考序列", "注释文件", "覆盖度", "平均深度", "说明"],
            [pre, species or "Human adenovirus", "-", "-", "-", "-", "-", "-", "-", "-", "HAdV 参考选择结果为空"],
        )
        return True
    row = rows[0]
    hadv_type = str(row.get("hadv_type") or "").strip()
    _write_placeholder(
        pre,
        ["样本名称", "病毒类型", "HAdV分型", "Penton分型", "Hexon分型", "Fiber分型", "参考序列", "注释文件", "覆盖度", "平均深度", "说明"],
        [
            pre,
            _hadv_species_label(hadv_type),
            hadv_type or "-",
            str(row.get("penton_type") or "-"),
            str(row.get("hexon_type") or "-"),
            str(row.get("fiber_type") or "-"),
            Path(str(row.get("reference_path") or "")).name if row.get("reference_path") else "-",
            Path(str(row.get("gff_path") or "")).name if row.get("gff_path") and row.get("gff_path") != "nogtf" else "-",
            f"{float(row.get('coverage') or 0.0):.2f}%" if row.get("coverage") else "-",
            f"{float(row.get('mean_depth') or 0.0):.2f}" if row.get("mean_depth") else "-",
            "先基于 PHF 三基因完成分型，再从对应分型全基因组候选库中选择覆盖度最优参考。",
        ],
    )
    return True


def _run_norovirus_reference_typing(pre: str, species: str) -> bool:
    selection_path = Path(f"{pre}_norovirus_reference_selection") / "selection.tsv"
    typing_summary_path = Path(f"{pre}_norovirus_reference_selection") / "typing" / "dual_typing.tsv"
    consensus_typing_path = Path(f"{pre}_norovirus_reference_selection") / "consensus_typing" / "consensus_typing.tsv"
    if not selection_path.is_file():
        _write_placeholder(
            pre,
            ["样本名称", "病毒类型", "双位点分型", "RdRp分型", "VP1分型", "参考序列", "注释文件", "覆盖度", "平均深度", "说明"],
            [pre, species or "Norovirus", "-", "-", "-", "-", "-", "-", "-", "未找到 Norovirus 参考选择结果"],
        )
        return True

    with selection_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        _write_placeholder(
            pre,
            ["样本名称", "病毒类型", "双位点分型", "RdRp分型", "VP1分型", "参考序列", "注释文件", "覆盖度", "平均深度", "说明"],
            [pre, species or "Norovirus", "-", "-", "-", "-", "-", "-", "-", "Norovirus 参考选择结果为空"],
        )
        return True

    row = rows[0]
    status_value = str(row.get("status") or "").strip().lower()
    if status_value == "skipped":
        skip_note = str(row.get("note") or "").strip() or "RdRp/VP1 分型支持 reads 过少，已跳过当前样本。"
        _write_placeholder(
            pre,
            ["样本名称", "病毒类型", "双位点分型", "RdRp分型", "VP1分型", "参考序列", "注释文件", "覆盖度", "平均深度", "说明"],
            [pre, species or "Norovirus", "-", "-", "-", "-", "-", "-", "-", skip_note],
        )
        return True
    dual_type = str(row.get("dual_type") or "").strip()
    rdrp_type = str(row.get("rdrp_type") or "").strip()
    vp1_type = str(row.get("vp1_type") or "").strip()
    typing_note = "基于 CDC RdRp/VP1 双位点分型结果，从对应候选全基因组参考中选择覆盖度最优参考。"
    if typing_summary_path.is_file():
        try:
            with typing_summary_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
                typing_rows = list(csv.DictReader(handle, delimiter="\t"))
            matched_types = {
                str(item.get("gene") or "").strip().lower(): str(item.get("matched_type") or "").strip()
                for item in typing_rows
            }
            rdrp_type = matched_types.get("rdrp") or rdrp_type
            vp1_type = matched_types.get("vp1") or vp1_type
        except OSError:
            pass
    if consensus_typing_path.is_file():
        try:
            with consensus_typing_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
                consensus_row = next(csv.DictReader(handle, delimiter="\t"), None) or {}
            consensus_dual = str(consensus_row.get("dual_type") or "").strip()
            if consensus_dual and consensus_dual != dual_type:
                typing_note = f"{typing_note} consensus 双位点分型为 {consensus_dual}。"
        except OSError:
            pass

    species_label = _norovirus_species_label(dual_type)
    _write_placeholder(
        pre,
        ["样本名称", "病毒类型", "双位点分型", "RdRp分型", "VP1分型", "参考序列", "注释文件", "覆盖度", "平均深度", "说明"],
        [
            pre,
            species_label if species_label != "Norovirus" else (species or "Norovirus"),
            dual_type or "-",
            rdrp_type or "-",
            vp1_type or "-",
            Path(str(row.get("reference_path") or "")).name if row.get("reference_path") else "-",
            Path(str(row.get("gff_path") or "")).name if row.get("gff_path") and row.get("gff_path") != "nogtf" else "-",
            f"{float(row.get('coverage') or 0.0):.2f}%" if row.get("coverage") else "-",
            f"{float(row.get('mean_depth') or 0.0):.2f}" if row.get("mean_depth") else "-",
            typing_note,
        ],
    )
    return True


def _run_hepatovirus_reference_typing(pre: str, species: str) -> bool:
    selection_path = Path(f"{pre}_hepatovirus_reference_selection") / "selection.tsv"
    consensus_typing_path = Path(f"{pre}_hepatovirus_reference_selection") / "consensus_typing" / "consensus_typing.tsv"
    columns = ["样本名称", "病毒类型", "大亚型", "子亚型", "参考序列", "注释文件", "覆盖度", "平均深度", "说明"]
    if not selection_path.is_file():
        final_fasta = Path(f"{pre}.final.fasta")
        if final_fasta.is_file() and final_fasta.stat().st_size > 0:
            resolve_hepatovirus_reference(
                pre,
                species=species,
                requested_ref="",
                query_fasta=str(final_fasta),
            )
    if not selection_path.is_file():
        _write_placeholder(
            pre,
            columns,
            [pre, species or "Hepatovirus", "-", "-", "-", "-", "-", "-", "未找到 Hepatovirus 参考选择结果"],
        )
        return True
    with selection_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        _write_placeholder(
            pre,
            columns,
            [pre, species or "Hepatovirus", "-", "-", "-", "-", "-", "-", "Hepatovirus 参考选择结果为空"],
        )
        return True
    row = rows[0]
    broad_type = str(row.get("broad_type") or "").strip().upper()
    subtype = str(row.get("subtype") or row.get("hav_subtype") or "").strip().upper()
    species_label = str(row.get("species_label") or "").strip() or _hepatovirus_species_label(broad_type, species or "Hepatovirus")
    note = str(row.get("note") or "").strip() or "先基于肝炎病毒 broad 参考库完成大亚型判定；再进入对应大亚型参考库选择覆盖度最优的子亚型参考。"
    if consensus_typing_path.is_file():
        try:
            with consensus_typing_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
                consensus_row = next(csv.DictReader(handle, delimiter="\t"), None) or {}
            consensus_broad = str(consensus_row.get("broad_type") or "").strip().upper()
            consensus_subtype = str(consensus_row.get("subtype") or consensus_row.get("hav_subtype") or "").strip().upper()
            if consensus_broad and consensus_broad != broad_type:
                note = f"{note} consensus 大亚型为 {consensus_broad}。"
            elif consensus_subtype and consensus_subtype != subtype:
                note = f"{note} consensus 子亚型为 {consensus_subtype}。"
        except OSError:
            pass
    _write_placeholder(
        pre,
        columns,
        [
            pre,
            species_label,
            broad_type or "-",
            subtype or "-",
            Path(str(row.get("reference_path") or "")).name if row.get("reference_path") else "-",
            Path(str(row.get("gff_path") or "")).name if row.get("gff_path") and row.get("gff_path") != "nogtf" else "-",
            f"{float(row.get('coverage') or 0.0):.2f}%" if row.get("coverage") else "-",
            f"{float(row.get('mean_depth') or 0.0):.2f}" if row.get("mean_depth") else "-",
            note,
        ],
    )
    return True


def _run_orthoebolavirus_reference_typing(pre: str, species: str) -> bool:
    selection_path = Path(f"{pre}_orthoebolavirus_reference_selection") / "selection.tsv"
    consensus_typing_path = Path(f"{pre}_orthoebolavirus_reference_selection") / "consensus_typing" / "consensus_typing.tsv"
    columns = ["样本名称", "病毒类型", "Orthoebolavirus分型", "病毒名称", "参考序列", "注释文件", "覆盖度", "平均深度", "说明"]
    if not selection_path.is_file():
        final_fasta = Path(f"{pre}.final.fasta")
        if final_fasta.is_file() and final_fasta.stat().st_size > 0:
            resolve_orthoebolavirus_reference(
                pre,
                species=species,
                requested_ref="",
                query_fasta=str(final_fasta),
            )
    if not selection_path.is_file():
        _write_placeholder(
            pre,
            columns,
            [pre, species or "Orthoebolavirus", "-", "-", "-", "-", "-", "-", "未找到 Orthoebolavirus 参考选择结果"],
        )
        return True
    with selection_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        _write_placeholder(
            pre,
            columns,
            [pre, species or "Orthoebolavirus", "-", "-", "-", "-", "-", "-", "Orthoebolavirus 参考选择结果为空"],
        )
        return True
    row = rows[0]
    predicted_type = str(row.get("predicted_type") or "").strip().upper()
    species_label = str(row.get("species_label") or "").strip() or _orthoebolavirus_species_label(predicted_type, species or "Orthoebolavirus")
    virus_name = str(row.get("virus_name") or "").strip()
    note = str(row.get("note") or "").strip() or "基于 Orthoebolavirus 本地参考基因组库完成最优参考选择。"
    if consensus_typing_path.is_file():
        try:
            with consensus_typing_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
                consensus_row = next(csv.DictReader(handle, delimiter="\t"), None) or {}
            consensus_type = str(consensus_row.get("predicted_type") or "").strip().upper()
            if consensus_type and consensus_type != predicted_type:
                note = f"{note} consensus 复判为 {consensus_type}。"
        except OSError:
            pass
    _write_placeholder(
        pre,
        columns,
        [
            pre,
            species_label,
            predicted_type or "-",
            virus_name or "-",
            Path(str(row.get("reference_path") or "")).name if row.get("reference_path") else "-",
            Path(str(row.get("gff_path") or "")).name if row.get("gff_path") and row.get("gff_path") != "nogtf" else "-",
            f"{float(row.get('coverage') or 0.0):.2f}%" if row.get("coverage") else "-",
            f"{float(row.get('mean_depth') or 0.0):.2f}" if row.get("mean_depth") else "-",
            note,
        ],
    )
    return True


def _run_rhinovirus_reference_typing(pre: str, species: str) -> bool:
    selection_path = Path(f"{pre}_rhinovirus_reference_selection") / "selection.tsv"
    consensus_typing_path = Path(f"{pre}_rhinovirus_reference_selection") / "consensus_typing" / "consensus_typing.tsv"
    if not selection_path.is_file():
        _write_placeholder(
            pre,
            ["样本名称", "病毒类型", "VP1分型", "物种组", "参考序列", "注释文件", "覆盖度", "平均深度", "说明"],
            [pre, species or "Human rhinovirus", "-", "-", "-", "-", "-", "-", "未找到 Rhinovirus 参考选择结果"],
        )
        return True
    with selection_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        _write_placeholder(
            pre,
            ["样本名称", "病毒类型", "VP1分型", "物种组", "参考序列", "注释文件", "覆盖度", "平均深度", "说明"],
            [pre, species or "Human rhinovirus", "-", "-", "-", "-", "-", "-", "Rhinovirus 参考选择结果为空"],
        )
        return True
    row = rows[0]
    vp1_type = str(row.get("vp1_type") or "").strip()
    species_group = str(row.get("species_group") or "").strip().upper()
    typing_note = "基于 VP1 分型结果，在对应鼻病毒全基因组候选集中选择覆盖度最优参考，并使用 VADR hrv 模型生成样本注释。"
    if consensus_typing_path.is_file():
        try:
            with consensus_typing_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
                consensus_row = next(csv.DictReader(handle, delimiter="\t"), None) or {}
            consensus_type = str(consensus_row.get("vp1_type") or "").strip()
            if consensus_type and consensus_type != vp1_type:
                typing_note = f"{typing_note} consensus VP1 分型为 {consensus_type}。"
        except OSError:
            pass
    _write_placeholder(
        pre,
        ["样本名称", "病毒类型", "VP1分型", "物种组", "参考序列", "注释文件", "覆盖度", "平均深度", "说明"],
        [
            pre,
            _rhinovirus_species_label(species_group) if species_group else (species or "Human rhinovirus"),
            vp1_type or "-",
            species_group or "-",
            Path(str(row.get("reference_path") or "")).name if row.get("reference_path") else "-",
            Path(str(row.get("gff_path") or "")).name if row.get("gff_path") and row.get("gff_path") != "nogtf" else "-",
            f"{float(row.get('coverage') or 0.0):.2f}%" if row.get("coverage") else "-",
            f"{float(row.get('mean_depth') or 0.0):.2f}" if row.get("mean_depth") else "-",
            typing_note,
        ],
    )
    return True


def _run_enterovirus_reference_typing(pre: str, species: str) -> bool:
    selection_path = Path(f"{pre}_enterovirus_reference_selection") / "selection.tsv"
    consensus_typing_path = Path(f"{pre}_enterovirus_reference_selection") / "consensus_typing" / "consensus_typing.tsv"
    phylogeny_summary_path = Path(f"{pre}_enterovirus_reference_selection") / "phylogeny" / "summary.tsv"
    columns = ["样本名称", "病毒类型", "大亚型", "VP1分型", "参考序列", "注释文件", "覆盖度", "平均深度", "说明"]

    def _phylogeny_ready(path: Path) -> bool:
        if not path.is_file() or path.stat().st_size == 0:
            return False
        try:
            with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
                phylo_row = next(csv.DictReader(handle, delimiter="\t"), None) or {}
            return str(phylo_row.get("status") or "").strip().lower() == "ready"
        except OSError:
            return False

    if not selection_path.is_file():
        _write_placeholder(
            pre,
            columns,
            [pre, species or "Human enterovirus", "-", "-", "-", "-", "-", "-", "未找到 Enterovirus 参考选择结果"],
        )
        return True
    with selection_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        _write_placeholder(
            pre,
            columns,
            [pre, species or "Human enterovirus", "-", "-", "-", "-", "-", "-", "Enterovirus 参考选择结果为空"],
        )
        return True
    row = rows[0]
    status_value = str(row.get("status") or "").strip().lower()
    if status_value == "skipped":
        skip_note = str(row.get("note") or "").strip() or "VP1 分型支持 reads 过少，已跳过当前样本。"
        _write_placeholder(
            pre,
            columns,
            [pre, species or "Human enterovirus", "-", "-", "-", "-", "-", "-", skip_note],
        )
        return True
    vp1_type = str(row.get("vp1_type") or "").strip()
    big_group = str(row.get("big_group") or "").strip().upper()
    consensus_fasta = Path(f"{pre}.consensus.fasta")
    final_fasta = Path(f"{pre}.final.fasta")
    typing_input = consensus_fasta if consensus_fasta.is_file() and consensus_fasta.stat().st_size > 0 else final_fasta
    if big_group in {"A", "B", "C", "D"} and typing_input.is_file() and typing_input.stat().st_size > 0 and not _phylogeny_ready(phylogeny_summary_path):
        sample_gff = prepare_enterovirus_sample_annotation(pre, typing_input, Path.cwd(), big_group)
        phylogeny_gff = sample_gff if sample_gff is not None and sample_gff.is_file() else _resolve_enterovirus_vp1_fallback_gff(vp1_type, big_group)
        if phylogeny_gff is not None and phylogeny_gff.is_file():
            build_enterovirus_vp1_phylogeny_assets(pre, typing_input, phylogeny_gff, big_group)
    nearest_reference = ""
    note_suffix = ""
    if phylogeny_summary_path.is_file():
        try:
            with phylogeny_summary_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
                phylo_row = next(csv.DictReader(handle, delimiter="\t"), None) or {}
            nearest_reference = str(phylo_row.get("nearest_tree_label") or phylo_row.get("nearest_accession") or "").strip()
        except OSError:
            pass
    if consensus_typing_path.is_file():
        try:
            with consensus_typing_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
                consensus_row = next(csv.DictReader(handle, delimiter="\t"), None) or {}
            consensus_type = str(consensus_row.get("vp1_type") or "").strip()
            if consensus_type and consensus_type != vp1_type:
                note_suffix = f" consensus VP1 分型为 {consensus_type}。"
        except OSError:
            pass
    typing_note = "基于 EV-A/B/C/D 的 VP1 分型结果，在对应亚型全基因组候选集中选择覆盖度最优参考；组装后再按大亚型使用 VADR 生成 GFF，并提取样本 VP1 做子亚型树。"
    if nearest_reference:
        typing_note = f"{typing_note} VP1 树最近参考为 {nearest_reference}。{note_suffix}".strip()
    _write_placeholder(
        pre,
        columns,
        [
            pre,
            _enterovirus_species_label(big_group) if big_group else (species or "Human enterovirus"),
            big_group or "-",
            vp1_type or "-",
            Path(str(row.get("reference_path") or "")).name if row.get("reference_path") else "-",
            Path(str(row.get("gff_path") or "")).name if row.get("gff_path") and row.get("gff_path") != "nogtf" else "-",
            f"{float(row.get('coverage') or 0.0):.2f}%" if row.get("coverage") else "-",
            f"{float(row.get('mean_depth') or 0.0):.2f}" if row.get("mean_depth") else "-",
            typing_note,
        ],
    )
    return True


def _run_seasonal_hcov_reference_typing(pre: str, species: str) -> bool:
    selection_path = Path(f"{pre}_seasonal_hcov_reference_selection") / "selection.tsv"
    consensus_typing_path = Path(f"{pre}_seasonal_hcov_reference_selection") / "consensus_typing" / "consensus_typing.tsv"
    phylogeny_summary_path = Path(f"{pre}_seasonal_hcov_reference_selection") / "phylogeny" / "summary.tsv"
    columns = ["样本名称", "病毒类型", "大类分型", "S子亚型", "最近参考", "注释文件", "覆盖度", "平均深度", "说明"]

    def _read_local_hint(path: Path, max_chars: int = 4096) -> str:
        try:
            if path.is_file() and path.stat().st_size > 0:
                return path.read_text(encoding="utf-8", errors="ignore")[:max_chars]
        except OSError:
            return ""
        return ""

    def _resolve_postassembly_hcov_type() -> str:
        explicit_type = _infer_seasonal_hcov_type_from_label(species)
        if explicit_type:
            return explicit_type
        for hint_path in [
            Path("ref/genes.gff"),
            Path("ref.txt"),
            Path("genomes/ref.fa"),
            Path(f"{pre}.consensus.fasta"),
            Path(f"{pre}.final.fasta"),
        ]:
            inferred = _infer_seasonal_hcov_type_from_ref_text(f"{hint_path.name}\n{_read_local_hint(hint_path)}")
            if inferred:
                return inferred
        return ""

    row: dict[str, str] = {}
    if selection_path.is_file():
        with selection_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        if rows:
            row = rows[0]
    hcov_type = str(row.get("hcov_type") or "").strip() or _resolve_postassembly_hcov_type()
    consensus_fasta = Path(f"{pre}.consensus.fasta")
    final_fasta = Path(f"{pre}.final.fasta")
    typing_input = consensus_fasta if consensus_fasta.is_file() and consensus_fasta.stat().st_size > 0 else final_fasta
    if hcov_type and typing_input.is_file() and typing_input.stat().st_size > 0 and (not consensus_typing_path.is_file() or not phylogeny_summary_path.is_file()):
        sample_gff = prepare_seasonal_hcov_sample_annotation(pre, typing_input, Path.cwd(), hcov_type)
        if sample_gff is not None and sample_gff.is_file():
            build_seasonal_hcov_spike_phylogeny_assets(pre, typing_input, sample_gff, hcov_type)
    if not hcov_type:
        _write_placeholder(pre, columns, [pre, species or "Human coronavirus", "-", "-", "-", "-", "-", "-", "未能确定季节性冠状病毒型别，无法继续执行 VADR 和 S 子亚型鉴定"])
        return True
    spike_subtype = ""
    nearest_reference = ""
    nearest_distance = ""
    note_suffix = ""
    if phylogeny_summary_path.is_file():
        try:
            with phylogeny_summary_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
                phylo_row = next(csv.DictReader(handle, delimiter="\t"), None) or {}
            spike_subtype = str(phylo_row.get("subtype") or "").strip()
            nearest_reference = str(phylo_row.get("nearest_tree_label") or phylo_row.get("nearest_accession") or "").strip()
            nearest_distance = str(phylo_row.get("nearest_distance") or "").strip()
        except OSError:
            pass
    if consensus_typing_path.is_file():
        try:
            with consensus_typing_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
                consensus_row = next(csv.DictReader(handle, delimiter="\t"), None) or {}
            consensus_subtype = str(consensus_row.get("subtype") or "").strip()
            if consensus_subtype and consensus_subtype != spike_subtype:
                note_suffix = f" consensus S 子亚型为 {consensus_subtype}。"
        except OSError:
            pass
    gff_name = "-"
    gff_path_text = str(row.get("gff_path") or "").strip()
    if gff_path_text and gff_path_text != "nogtf":
        gff_name = Path(gff_path_text).name
    else:
        sample_vadr_gff = Path("vadr") / f"{pre}.vadr.gff3"
        if sample_vadr_gff.is_file() and sample_vadr_gff.stat().st_size > 0:
            gff_name = sample_vadr_gff.name
    typing_note = "先用四种季节性冠状病毒参考基因组判定病毒类型，再对 consensus 序列进行 VADR 注释并提取 S 基因，结合最近距离参考判断子亚型。"
    if not selection_path.is_file():
        typing_note = "当前为手动指定 ref/gtf 或无参考筛选结果场景，已直接基于样本 consensus/组装序列执行 VADR 注释和 S 基因子亚型鉴定。"
    if nearest_distance:
        typing_note = f"{typing_note} 最近参考距离 {nearest_distance}。{note_suffix}".strip()
    _write_placeholder(
        pre,
        columns,
        [
            pre,
            _seasonal_hcov_species_label(hcov_type) if hcov_type else (species or "Human coronavirus"),
            hcov_type or "-",
            spike_subtype or "-",
            nearest_reference or "-",
            gff_name,
            f"{float(row.get('coverage') or 0.0):.2f}%" if row.get("coverage") else "-",
            f"{float(row.get('mean_depth') or 0.0):.2f}" if row.get("mean_depth") else "-",
            typing_note,
        ],
    )
    return True


def _resolve_hiv_db_dir() -> Path:
    env_root = str(os.environ.get("META_HIV_DB_DIR") or "").strip()
    if env_root:
        return Path(env_root).expanduser().resolve()
    database_root = str(os.environ.get("META_DATABASE_ROOT") or "").strip()
    project_root = _project_root()
    candidates: list[Path] = []
    if database_root:
        candidates.append((Path(database_root) / "virus" / "HIV").expanduser())
    candidates.extend(
        [
            project_root / "database" / "virus" / "HIV",
            Path("/data/deploy/meta_genome/database/virus/HIV"),
        ]
    )
    return next((candidate.resolve() for candidate in candidates if candidate.exists()), candidates[0].resolve())


def _resolve_hiv_broad_reference_paths() -> dict[str, Path]:
    db_dir = _resolve_hiv_db_dir()
    return {
        "HIV-1": (db_dir / "HIV1.fasta").resolve(),
        "HIV-2": (db_dir / "HIV2.fasta").resolve(),
    }


def _resolve_hiv_rega_reference_fasta() -> Path:
    return (_resolve_hiv_db_dir() / "rega_reference_genomes" / "rega_hiv_reference_genomes.fasta").resolve()


def _resolve_hiv_rega_reference_manifest() -> Path:
    return (_resolve_hiv_db_dir() / "rega_reference_genomes" / "reference_manifest.tsv").resolve()


def _resolve_hiv_rega_supplement_manifest() -> Path:
    return (_resolve_hiv_db_dir() / "rega_reference_genomes" / "supplement_manifest.tsv").resolve()


def _normalize_hiv_subtype_label(group: str) -> str:
    raw = str(group or "").strip()
    if raw.startswith("Subtype "):
        return raw.split(" ", 1)[1].strip().replace("-", "_").upper()
    return raw.replace("-", "_").upper()


def _build_hiv_broad_reference_records() -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for broad_type, fasta_path in _resolve_hiv_broad_reference_paths().items():
        header, sequence = _read_first_fasta_sequence(fasta_path)
        if not sequence:
            continue
        fasta_id = broad_type.replace("-", "_")
        records.append(
            {
                "fasta_id": fasta_id,
                "accession": fasta_id,
                "sequence": sequence,
                "reference_length": len(sequence),
                "broad_type": broad_type,
                "subtype": "",
                "meta": {
                    "header": header or broad_type,
                    "species_label": broad_type,
                    "fasta_path": str(fasta_path.resolve()),
                },
            }
        )
    return records


def _load_hiv_subtype_representative_records() -> list[dict[str, object]]:
    combined_fasta = _resolve_hiv_rega_reference_fasta()
    manifest_path = _resolve_hiv_rega_reference_manifest()
    supplement_manifest_path = _resolve_hiv_rega_supplement_manifest()
    if not combined_fasta.is_file() or not manifest_path.is_file() or not supplement_manifest_path.is_file():
        return []
    fasta_records: dict[str, tuple[str, str]] = {}
    for record in SeqIO.parse(str(combined_fasta), "fasta"):
        record_id = str(record.id).strip()
        record_desc = str(record.description or record.id).strip()
        record_seq = str(record.seq).upper()
        fasta_records[record_id] = (record_desc, record_seq)
        fasta_records.setdefault(record_id.split(".", 1)[0], (record_desc, record_seq))

    representatives: dict[str, dict[str, object]] = {}

    def _maybe_add_record(accession: str, group: str, source_tag: str) -> None:
        accession = str(accession or "").strip()
        if not accession or accession not in fasta_records:
            return
        subtype = _normalize_hiv_subtype_label(group)
        if not subtype or subtype in representatives:
            return
        header, sequence = fasta_records[accession]
        representatives[subtype] = {
            "fasta_id": accession,
            "accession": accession,
            "sequence": sequence,
            "reference_length": len(sequence),
            "broad_type": "HIV-1",
            "subtype": subtype,
            "meta": {
                "header": header,
                "group": group,
                "source": source_tag,
                "fasta_path": str(combined_fasta.resolve()),
            },
        }

    with manifest_path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if str(row.get("download_status") or "").strip() != "downloaded":
                continue
            _maybe_add_record(str(row.get("note") or "").strip(), str(row.get("group") or "").strip(), "rega_manifest")
    with supplement_manifest_path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            _maybe_add_record(str(row.get("accession") or "").strip(), str(row.get("group") or "").strip(), "supplement_manifest")
    return [representatives[key] for key in sorted(representatives)]


def resolve_hiv_reference(
    pre: str,
    species: str = "",
    requested_ref: str = "",
    single_fastq: str = "",
    fq1: str = "",
    fq2: str = "",
    long_type: str = "",
    threads: int = 4,
    query_fasta: str = "",
    output_dir: Path | None = None,
    logf=None,
) -> dict[str, str]:
    requested = str(requested_ref or "").strip()
    if requested and Path(requested).is_file():
        broad_type = "HIV-2" if "hiv2" in Path(requested).name.lower() else "HIV-1"
        return {
            "status": "ready",
            "broad_type": broad_type,
            "subtype": "",
            "species_label": broad_type,
            "reference_path": str(Path(requested).resolve()),
            "gff_path": "nogtf",
            "summary_path": "",
            "broad_summary_path": "",
            "subtype_summary_path": "",
        }

    screening_dir = output_dir if output_dir is not None else Path(f"{pre}_hiv_reference_selection")
    screening_dir.mkdir(parents=True, exist_ok=True)
    broad_summary_path = screening_dir / "broad_typing.tsv"
    final_summary_path = screening_dir / "selection.tsv"
    query_path = Path(query_fasta) if str(query_fasta or "").strip() else None

    broad_records = _build_hiv_broad_reference_records()
    if not broad_records:
        with final_summary_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sample", "species_label", "broad_type", "subtype", "coverage", "mean_depth", "covered_bases", "num_reads", "reference_name", "reference_path", "gff_path", "status", "note", "broad_summary_path", "subtype_summary_path"])
            writer.writerow([pre, species or "HIV", "", "", "", "", "", "", "", "", "nogtf", "missing", "未找到 HIV-1/HIV-2 broad 参考库", "", ""])
        return {
            "status": "missing",
            "broad_type": "",
            "subtype": "",
            "species_label": species or "HIV",
            "reference_path": "",
            "gff_path": "nogtf",
            "summary_path": str(final_summary_path.resolve()),
            "broad_summary_path": "",
            "subtype_summary_path": "",
        }

    broad_pick = _select_hepatovirus_reference_records(
        broad_records,
        screening_dir,
        "broad_typing",
        single_fastq=single_fastq,
        fq1=fq1,
        fq2=fq2,
        long_type=long_type,
        threads=threads,
        query_fasta=query_path,
        logf=logf,
    )
    broad_best_row = broad_pick.get("best_row") or {}
    broad_best_meta = broad_pick.get("best_meta") or {}
    broad_meta_row = dict(broad_best_meta.get("meta") or {})
    broad_type = str(broad_best_meta.get("broad_type") or "").strip()
    broad_species_label = str(broad_meta_row.get("species_label") or broad_type or species or "HIV").strip()
    broad_reference_path = ""
    if broad_best_meta:
        broad_reference_path = str((screening_dir / f"{_sanitize_tree_label(broad_type or 'hiv')}.reference.fasta").resolve())
        if not _extract_named_fasta_record(Path(str(broad_pick.get("candidate_fasta") or "")), str(broad_best_meta.get("fasta_id") or ""), Path(broad_reference_path)):
            broad_reference_path = ""
    broad_note = "先基于 HIV-1/HIV-2 代表全基因组参考库按覆盖度/序列覆盖度判定大亚型。"
    with broad_summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "broad_type", "coverage", "mean_depth", "covered_bases", "num_reads", "reference_name", "reference_path", "status", "note"])
        writer.writerow([
            pre,
            broad_type,
            f"{float(broad_best_row.get('coverage') or 0.0):.6f}",
            f"{float(broad_best_row.get('mean_depth') or 0.0):.6f}" if str(broad_best_row.get("mean_depth") or "").strip() else "",
            f"{float(broad_best_row.get('covered_bases') or 0.0):.0f}" if str(broad_best_row.get("covered_bases") or "").strip() else "",
            f"{float(broad_best_row.get('num_reads') or 0.0):.0f}" if str(broad_best_row.get("num_reads") or "").strip() else "",
            str(broad_meta_row.get("header") or broad_type or ""),
            broad_reference_path,
            "ready" if broad_best_meta else "missing",
            broad_note if broad_best_meta else "未获得 HIV-1/HIV-2 broad 命中结果",
        ])

    if not broad_best_meta:
        with final_summary_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sample", "species_label", "broad_type", "subtype", "coverage", "mean_depth", "covered_bases", "num_reads", "reference_name", "reference_path", "gff_path", "status", "note", "broad_summary_path", "subtype_summary_path"])
            writer.writerow([pre, species or "HIV", "", "", "", "", "", "", "", "", "nogtf", "missing", "未获得 HIV-1/HIV-2 broad 命中结果", str(broad_summary_path.resolve()), ""])
        return {
            "status": "missing",
            "broad_type": "",
            "subtype": "",
            "species_label": species or "HIV",
            "reference_path": "",
            "gff_path": "nogtf",
            "summary_path": str(final_summary_path.resolve()),
            "broad_summary_path": str(broad_summary_path.resolve()),
            "subtype_summary_path": "",
        }

    selected_reference_path = broad_reference_path
    selected_row = broad_best_row
    selected_reference_name = str(broad_meta_row.get("header") or broad_type or "")
    subtype = ""
    subtype_summary_resolved = ""
    final_note = broad_note

    if broad_type == "HIV-1":
        subtype_records = _load_hiv_subtype_representative_records()
        subtype_summary_path = screening_dir / "representative_subtype_typing.tsv"
        if subtype_records:
            subtype_pick = _select_hepatovirus_reference_records(
                subtype_records,
                screening_dir,
                "representative_subtype_typing",
                single_fastq=single_fastq,
                fq1=fq1,
                fq2=fq2,
                long_type=long_type,
                threads=threads,
                query_fasta=query_path,
                logf=logf,
            )
            subtype_best_row = subtype_pick.get("best_row") or {}
            subtype_best_meta = subtype_pick.get("best_meta") or {}
            subtype_meta_row = dict(subtype_best_meta.get("meta") or {})
            subtype = str(subtype_best_meta.get("subtype") or "").strip()
            subtype_reference_path = ""
            if subtype_best_meta:
                subtype_reference_path = str((screening_dir / f"{_sanitize_tree_label(subtype or 'hiv1')}.reference.fasta").resolve())
                if not _extract_named_fasta_record(Path(str(subtype_pick.get("candidate_fasta") or "")), str(subtype_best_meta.get("fasta_id") or ""), Path(subtype_reference_path)):
                    subtype_reference_path = ""
            with subtype_summary_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["sample", "broad_type", "subtype", "coverage", "mean_depth", "covered_bases", "num_reads", "reference_name", "reference_path", "status", "note"])
                writer.writerow([
                    pre,
                    broad_type,
                    subtype,
                    f"{float(subtype_best_row.get('coverage') or 0.0):.6f}",
                    f"{float(subtype_best_row.get('mean_depth') or 0.0):.6f}" if str(subtype_best_row.get("mean_depth") or "").strip() else "",
                    f"{float(subtype_best_row.get('covered_bases') or 0.0):.0f}" if str(subtype_best_row.get("covered_bases") or "").strip() else "",
                    f"{float(subtype_best_row.get('num_reads') or 0.0):.0f}" if str(subtype_best_row.get("num_reads") or "").strip() else "",
                    str(subtype_meta_row.get("header") or subtype_best_meta.get("accession") or ""),
                    subtype_reference_path,
                    "ready" if subtype_best_meta else "missing",
                    "基于 HIV-1 各子亚型代表株覆盖度选择一致性生成参考。" if subtype_best_meta else "未获得 HIV-1 子亚型代表株稳定命中结果",
                ])
            subtype_summary_resolved = str(subtype_summary_path.resolve())
            if subtype_best_meta and subtype_reference_path:
                selected_reference_path = subtype_reference_path
                selected_row = subtype_best_row
                selected_reference_name = str(subtype_meta_row.get("header") or subtype_best_meta.get("accession") or "")
                final_note = "先基于 HIV-1/HIV-2 broad 参考库区分大亚型；命中 HIV-1 后，再按各子亚型代表株覆盖度选择参考生成一致性序列，随后执行 HIV-1 子亚型/重组与耐药分析。"

    with final_summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "species_label", "broad_type", "subtype", "coverage", "mean_depth", "covered_bases", "num_reads", "reference_name", "reference_path", "gff_path", "status", "note", "broad_summary_path", "subtype_summary_path"])
        writer.writerow([
            pre,
            broad_species_label,
            broad_type,
            subtype,
            f"{float(selected_row.get('coverage') or 0.0):.6f}",
            f"{float(selected_row.get('mean_depth') or 0.0):.6f}" if str(selected_row.get("mean_depth") or "").strip() else "",
            f"{float(selected_row.get('covered_bases') or 0.0):.0f}" if str(selected_row.get("covered_bases") or "").strip() else "",
            f"{float(selected_row.get('num_reads') or 0.0):.0f}" if str(selected_row.get("num_reads") or "").strip() else "",
            selected_reference_name,
            selected_reference_path,
            "nogtf",
            "ready" if selected_reference_path else "missing",
            final_note,
            str(broad_summary_path.resolve()),
            subtype_summary_resolved,
        ])
    return {
        "status": "ready" if selected_reference_path else "missing",
        "broad_type": broad_type,
        "subtype": subtype,
        "species_label": broad_species_label,
        "reference_path": selected_reference_path,
        "gff_path": "nogtf",
        "summary_path": str(final_summary_path.resolve()),
        "broad_summary_path": str(broad_summary_path.resolve()),
        "subtype_summary_path": subtype_summary_resolved,
    }


def _run_hiv_resistance_annotation(pre: str, species: str, input_fasta: Path) -> dict[str, object]:
    project_root = Path(__file__).resolve().parent.parent
    script_path = project_root / "scripts" / "annotate_hivdb_resistance.py"
    xml_path = project_root / "database" / "virus" / "HIV" / "HIVDB_10.2.xml"
    detail_tsv_path = Path(f"{pre}_hiv_resistance.tsv")
    detail_json_path = Path(f"{pre}_hiv_resistance.json")
    empty_result = {
        "status": "missing",
        "note": "",
        "class_best": {},
        "pr_mutations": "",
        "rt_mutations": "",
        "in_mutations": "",
        "alerts_text": "",
        "detail_tsv_path": str(detail_tsv_path.resolve()),
        "detail_json_path": str(detail_json_path.resolve()),
    }
    if not script_path.is_file() or not xml_path.is_file():
        empty_result["note"] = "未找到 HIVDB 注释脚本或 XML 规则库"
        return empty_result
    if not input_fasta.is_file() or input_fasta.stat().st_size == 0:
        empty_result["note"] = "未检测到可用于 HIV 耐药分析的 consensus/组装序列"
        return empty_result
    cmd = [
        sys.executable,
        str(script_path),
        "--json",
        "--xml",
        str(xml_path),
        "--fasta",
        str(input_fasta),
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as exc:
        empty_result["status"] = "failed"
        empty_result["note"] = (exc.stderr or exc.stdout or "").strip() or "HIVDB 注释脚本执行失败"
        return empty_result
    try:
        payload = json.loads(result.stdout)
    except json.JSONDecodeError:
        empty_result["status"] = "failed"
        empty_result["note"] = "HIVDB 注释输出不是有效 JSON"
        return empty_result

    samples = payload.get("samples") if isinstance(payload, dict) else None
    sample_entry = samples[0] if isinstance(samples, list) and samples else {}
    input_mutations = sample_entry.get("input_mutations") if isinstance(sample_entry, dict) else {}
    drug_results = sample_entry.get("drug_results") if isinstance(sample_entry, dict) else []
    sequence_alerts = sample_entry.get("sequence_alerts") if isinstance(sample_entry, dict) else []
    filtered_drug_results = [
        row for row in (drug_results if isinstance(drug_results, list) else [])
        if str(row.get("drug_class") or "").strip() in {"NRTI", "NNRTI", "PI", "INSTI"}
    ]
    class_best: dict[str, tuple[int, str]] = {}
    for row in filtered_drug_results:
        drug_class = str(row.get("drug_class") or "").strip()
        level = int(row.get("level") or 0)
        level_name = str(row.get("level_name") or "").strip() or "-"
        previous = class_best.get(drug_class)
        if previous is None or level > previous[0]:
            class_best[drug_class] = (level, level_name)
    pr_mutations = ",".join(input_mutations.get("PR") or []) if isinstance(input_mutations, dict) else ""
    rt_mutations = ",".join(input_mutations.get("RT") or []) if isinstance(input_mutations, dict) else ""
    in_mutations = ",".join(input_mutations.get("IN") or []) if isinstance(input_mutations, dict) else ""
    alerts_text = "；".join(sequence_alerts) if isinstance(sequence_alerts, list) else ""

    with detail_json_path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
    with detail_tsv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["drug_class", "drug", "fullname", "score", "level", "level_name", "sir", "triggered_rules", "result_comments"])
        for row in filtered_drug_results:
            writer.writerow([
                str(row.get("drug_class") or ""),
                str(row.get("drug") or ""),
                str(row.get("fullname") or ""),
                str(row.get("score") or ""),
                str(row.get("level") or ""),
                str(row.get("level_name") or ""),
                str(row.get("sir") or ""),
                "; ".join(row.get("triggered_rules") or []) if isinstance(row.get("triggered_rules"), list) else "",
                "；".join(row.get("result_comments") or []) if isinstance(row.get("result_comments"), list) else "",
            ])
    return {
        "status": "ready",
        "note": f"基于 HIVDB 10.2 XML 规则对 PR/RT/IN 突变进行药物耐药解释；输入序列为 {input_fasta.name}。",
        "class_best": class_best,
        "pr_mutations": pr_mutations,
        "rt_mutations": rt_mutations,
        "in_mutations": in_mutations,
        "alerts_text": alerts_text,
        "detail_tsv_path": str(detail_tsv_path.resolve()),
        "detail_json_path": str(detail_json_path.resolve()),
    }


def _run_hiv_typing_workflow(pre: str, species: str) -> bool:
    project_root = Path(__file__).resolve().parent.parent
    typing_script = project_root / "scripts" / "type_hiv_rega_like.py"
    consensus_fasta = Path(f"{pre}.consensus.fasta")
    hxb2_consensus_fasta = Path(f"{pre}.hxb2.consensus.fasta")
    final_fasta = Path(f"{pre}.final.fasta")
    input_fasta = consensus_fasta if consensus_fasta.is_file() and consensus_fasta.stat().st_size > 0 else final_fasta
    summary_path = Path(f"{pre}_serotype_result.tsv")
    selection_path = Path(f"{pre}_hiv_reference_selection") / "selection.tsv"
    subtype_json_path = Path(f"{pre}_hiv_subtyping.json")
    columns = [
        "样本名称",
        "病毒类型",
        "大亚型",
        "子亚型",
        "重组判定",
        "候选父本",
        "代表株参考",
        "输入序列",
        "耐药输入序列",
        "NRTI最高等级",
        "NNRTI最高等级",
        "PI最高等级",
        "INSTI最高等级",
        "PR突变",
        "RT突变",
        "IN突变",
        "序列告警",
        "说明",
    ]
    if not input_fasta.is_file() or input_fasta.stat().st_size == 0:
        _write_placeholder(pre, columns, [pre, species or "HIV", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "未检测到可用于 HIV 分型/耐药分析的 consensus/组装序列"])
        return False

    if not selection_path.is_file():
        resolve_hiv_reference(pre, species=species, query_fasta=str(input_fasta))

    selection_row: dict[str, str] = {}
    if selection_path.is_file():
        with selection_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
            selection_row = next(csv.DictReader(handle, delimiter="\t"), None) or {}
    broad_type = str(selection_row.get("broad_type") or "").strip() or ("HIV-2" if "hiv-2" in str(species or "").lower() else "HIV-1")
    representative_subtype = str(selection_row.get("subtype") or "").strip()
    representative_reference = Path(str(selection_row.get("reference_path") or "")).name if str(selection_row.get("reference_path") or "").strip() else "-"
    note_parts: list[str] = []
    if str(selection_row.get("note") or "").strip():
        note_parts.append(str(selection_row.get("note") or "").strip())

    subtype_name = representative_subtype or "-"
    recombination_label = "-"
    candidate_parents_text = "-"
    if broad_type == "HIV-1" and typing_script.is_file():
        subtype_output_dir = Path(f"{pre}_hiv_rega_like")
        cmd = [sys.executable, str(typing_script), "--json", "--output-dir", str(subtype_output_dir), str(input_fasta)]
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            payload = json.loads(result.stdout)
            with subtype_json_path.open("w", encoding="utf-8") as handle:
                json.dump(payload, handle, ensure_ascii=False, indent=2)
            samples = payload.get("samples") if isinstance(payload, dict) else None
            subtype_entry = samples[0] if isinstance(samples, list) and samples else {}
            subtype_name = str(subtype_entry.get("predicted_group") or representative_subtype or "").strip() or "-"
            recombination_label = str(subtype_entry.get("assignment_label") or "").strip() or "-"
            parents = subtype_entry.get("candidate_parents") if isinstance(subtype_entry, dict) else []
            parent_parts: list[str] = []
            for item in parents[:3] if isinstance(parents, list) else []:
                group = str(item.get("group") or "").strip()
                fraction = str(item.get("fraction") or "").strip()
                if group:
                    parent_parts.append(f"{group}({fraction})" if fraction else group)
            candidate_parents_text = "；".join(parent_parts) if parent_parts else "-"
            notes = subtype_entry.get("notes") if isinstance(subtype_entry, dict) else []
            if isinstance(notes, list):
                note_parts.extend(str(item).strip() for item in notes[:3] if str(item).strip())
        except (subprocess.CalledProcessError, json.JSONDecodeError) as exc:
            note_parts.append(f"HIV-1 子亚型/重组分析失败：{exc}")
    elif broad_type == "HIV-1":
        note_parts.append("未找到 HIV-1 子亚型分析脚本。")
    else:
        note_parts.append("当前判定为 HIV-2，不执行 HIV-1 子亚型/重组与 HIVDB 耐药分析。")

    resistance_result = {
        "status": "missing",
        "note": "",
        "class_best": {},
        "pr_mutations": "",
        "rt_mutations": "",
        "in_mutations": "",
        "alerts_text": "",
    }
    resistance_input_fasta = hxb2_consensus_fasta if hxb2_consensus_fasta.is_file() and hxb2_consensus_fasta.stat().st_size > 0 else input_fasta
    if broad_type == "HIV-1":
        if resistance_input_fasta != input_fasta:
            note_parts.append(f"耐药分析使用 HXB2 参考重建的一致性序列：{resistance_input_fasta.name}")
        else:
            note_parts.append(f"未检测到 HXB2 专用一致性序列，耐药分析回退使用：{resistance_input_fasta.name}")
        resistance_result = _run_hiv_resistance_annotation(pre, broad_type, resistance_input_fasta)
        if str(resistance_result.get("note") or "").strip():
            note_parts.append(str(resistance_result.get("note") or "").strip())
        if str(resistance_result.get("alerts_text") or "").strip():
            note_parts.append(f"序列告警：{resistance_result.get('alerts_text')}")

    class_best = resistance_result.get("class_best") if isinstance(resistance_result, dict) else {}
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(columns)
        writer.writerow([
            pre,
            broad_type or species or "HIV",
            broad_type or "-",
            subtype_name or "-",
            recombination_label or "-",
            candidate_parents_text or "-",
            representative_reference or "-",
            input_fasta.name,
            resistance_input_fasta.name if broad_type == "HIV-1" else "-",
            (class_best.get("NRTI") or (0, "-"))[1] if isinstance(class_best, dict) else "-",
            (class_best.get("NNRTI") or (0, "-"))[1] if isinstance(class_best, dict) else "-",
            (class_best.get("PI") or (0, "-"))[1] if isinstance(class_best, dict) else "-",
            (class_best.get("INSTI") or (0, "-"))[1] if isinstance(class_best, dict) else "-",
            str(resistance_result.get("pr_mutations") or "-"),
            str(resistance_result.get("rt_mutations") or "-"),
            str(resistance_result.get("in_mutations") or "-"),
            str(resistance_result.get("alerts_text") or "-"),
            " ".join(part for part in note_parts if part).strip() or "-",
        ])
    return True


def virus_typing(pre: str, species: str) -> None:
    columns = ["样本名称", "病毒类型", "HA亚型", "NA亚型", "分型结果", "HA参考命中", "NA参考命中", "说明"]
    final_fasta = Path(f"{pre}.final.fasta")
    consensus_fasta = Path(f"{pre}.consensus.fasta")
    if _is_bandavirus(species):
        _run_bandavirus_serotype_typing(pre, species, final_fasta, consensus_fasta)
        return
    if _is_orthohantavirus(species):
        _run_orthohantavirus_serotype_typing(pre, species, final_fasta, consensus_fasta)
        return
    if _is_orthoebolavirus(species):
        _run_orthoebolavirus_reference_typing(pre, species)
        if final_fasta.is_file() and final_fasta.stat().st_size > 0:
            _run_ebola_nextclade_typing(pre, species, final_fasta)
        return
    if _is_hepatovirus(species):
        _run_hepatovirus_reference_typing(pre, species)
        return
    if _is_hiv(species):
        _run_hiv_typing_workflow(pre, species)
        return
    if not final_fasta.is_file() or final_fasta.stat().st_size == 0:
        inferred = _infer_influenza_type_from_label(species)
        if inferred in {"Influenza C virus", "Influenza D virus"}:
            _write_placeholder(pre, columns, [pre, inferred, "-", "-", "-", "-", "-", "已完成物种鉴定；当前对丙流/丁流不继续进行组装与后续分型"])
            return
        _write_placeholder(pre, columns, [pre, species or "-", "-", "-", "-", "-", "-", "未检测到可用于分型的组装结果"])
        return

    if _is_sars_cov_2(species):
        _run_sars_cov_2_nextclade_typing(pre, species, final_fasta)
        return

    if _is_monkeypox(species):
        _run_monkeypox_nextclade_typing(pre, species, final_fasta)
        return

    if _is_rsv(species):
        _run_rsv_nextclade_typing(pre, species, final_fasta)
        return

    if _is_hmpv(species):
        _run_hmpv_nextclade_typing(pre, species, final_fasta)
        return

    if _is_denv(species):
        _run_denv_nextclade_typing(pre, species, final_fasta)
        return

    if _is_zika(species):
        _run_zika_nextclade_typing(pre, species, final_fasta)
        return

    if _is_chikv(species):
        _run_chikv_nextclade_typing(pre, species, final_fasta)
        return

    if _is_hpiv(species):
        _run_hpiv_reference_typing(pre, species)
        return

    if _is_hadv(species):
        _run_hadv_reference_typing(pre, species)
        return

    if _is_norovirus(species):
        _run_norovirus_reference_typing(pre, species)
        return

    if _is_rhinovirus(species):
        _run_rhinovirus_reference_typing(pre, species)
        return

    if _is_enterovirus(species):
        _run_enterovirus_reference_typing(pre, species)
        return

    if _is_seasonal_hcov(species):
        _run_seasonal_hcov_reference_typing(pre, species)
        return

    if not _is_influenza(species):
        _write_placeholder(pre, columns, [pre, species or "-", "-", "-", "-", "-", "-", "当前仅内置流感病毒分型流程"])
        return

    type_source = _resolve_db_path("META_FLU_TYPE_DB", "/data/deploy/meta_genome/database/virus/influenza/type_refs.fa")
    ha_source = _resolve_db_path("META_FLUA_HA_DB", "/data/deploy/meta_genome/database/virus/influenza_a/ha_subtypes.fa")
    na_source = _resolve_db_path("META_FLUA_NA_DB", "/data/deploy/meta_genome/database/virus/influenza_a/na_subtypes.fa")
    if not type_source.is_file():
        _write_placeholder(pre, columns, [pre, species or "-", "-", "-", "-", "-", "-", "未找到流感 A/B 判型参考数据库"])
        return

    blast_dir = Path(f"{pre}_virus_typing")
    blast_dir.mkdir(exist_ok=True)
    with open(blast_dir / "virus_typing.log", "w", encoding="utf-8") as logf:
        type_db = _ensure_blast_db(type_source, blast_dir / "influenza_type", logf=logf)
        type_hit = _blast_type_hit(final_fasta, type_db, blast_dir / "type.top.tsv", logf=logf)
        inferred_type = type_hit["virus_type"]
        if inferred_type == "-" and _is_influenza_a(species):
            inferred_type = "Influenza A virus"
        elif inferred_type == "-" and _is_influenza_b(species):
            inferred_type = "Influenza B virus"

        if inferred_type == "Influenza B virus":
            note = "流感 A/B 判型结果支持乙型流感病毒。"
            if type_hit["identity"] != "-" or type_hit["coverage"] != "-":
                note = f"流感 A/B 判型命中一致性 {type_hit['identity']}，覆盖度 {type_hit['coverage']}。"
            vadr_result = _run_vadr_flu_annotation(pre, final_fasta, blast_dir, logf=logf)
            snpeff_result = {"status": "missing", "note": "", "annotated_vcf": "", "table_json": ""}
            if vadr_result["status"] == "ready" and vadr_result.get("gff_path"):
                preferred_fasta = consensus_fasta if consensus_fasta.is_file() and consensus_fasta.stat().st_size > 0 else final_fasta
                snpeff_result = _run_influenza_snpeff_annotation(
                    pre,
                    preferred_fasta,
                    Path(vadr_result["gff_path"]),
                    Path("snps.filt1.vcf"),
                    blast_dir,
                    logf=logf,
                )
            if vadr_result["status"] == "ready":
                note += f" VADR 注释目录: {vadr_result['output_dir']}。"
            elif vadr_result["status"] == "failed":
                note += f" {vadr_result['note']}。"
            if snpeff_result["status"] == "ready":
                note += " 已生成 snpEff 变异注释表。"
            elif snpeff_result["status"] == "failed":
                note += f" {snpeff_result['note']}。"
            _write_placeholder(
                pre,
                columns,
                [pre, inferred_type, "-", "-", "B型流感", "-", "-", note],
            )
            return

        if inferred_type != "Influenza A virus":
            _write_placeholder(pre, columns, [pre, species or "-", "-", "-", "-", "-", "-", "未获得稳定的流感 A/B 判型结果"])
            return

        if not ha_source.is_file() or not na_source.is_file():
            _write_placeholder(pre, columns, [pre, inferred_type, "-", "-", "-", "-", "-", "已判定为甲型流感，但未找到 HA/NA 分型参考数据库"])
            return

        ha_db = _ensure_blast_db(ha_source, blast_dir / "ha_subtypes", logf=logf)
        na_db = _ensure_blast_db(na_source, blast_dir / "na_subtypes", logf=logf)
        ha_hit = _blast_top_hit(final_fasta, ha_db, blast_dir / "ha.top.tsv", "H", logf=logf)
        na_hit = _blast_top_hit(final_fasta, na_db, blast_dir / "na.top.tsv", "N", logf=logf)
        vadr_result = _run_vadr_flu_annotation(pre, final_fasta, blast_dir, logf=logf)
        snpeff_result = {"status": "missing", "note": "", "annotated_vcf": "", "table_json": ""}
        resistance_result = {"status": "missing", "note": "", "table_json": "", "table_tsv": ""}
        if vadr_result["status"] == "ready" and vadr_result.get("gff_path"):
            preferred_fasta = consensus_fasta if consensus_fasta.is_file() and consensus_fasta.stat().st_size > 0 else final_fasta
            snpeff_result = _run_influenza_snpeff_annotation(
                pre,
                preferred_fasta,
                Path(vadr_result["gff_path"]),
                Path("snps.filt1.vcf"),
                blast_dir,
                logf=logf,
            )
        if snpeff_result["status"] == "ready" and snpeff_result.get("table_json"):
            resistance_result = _run_influenza_resistance_annotation(
                Path.cwd(),
                Path(snpeff_result["table_json"]),
                logf=logf,
            )

    subtype = "-"
    if ha_hit["subtype"] != "-" and na_hit["subtype"] != "-":
        subtype = f"{ha_hit['subtype']}{na_hit['subtype']}"
    elif ha_hit["subtype"] != "-":
        subtype = ha_hit["subtype"]
    elif na_hit["subtype"] != "-":
        subtype = na_hit["subtype"]

    note_parts = []
    if type_hit["identity"] != "-" or type_hit["coverage"] != "-":
        note_parts.append(f"A/B 判型一致性 {type_hit['identity']}，覆盖度 {type_hit['coverage']}")
    if ha_hit["identity"] != "-" or ha_hit["coverage"] != "-":
        note_parts.append(f"HA 命中一致性 {ha_hit['identity']}，覆盖度 {ha_hit['coverage']}")
    if na_hit["identity"] != "-" or na_hit["coverage"] != "-":
        note_parts.append(f"NA 命中一致性 {na_hit['identity']}，覆盖度 {na_hit['coverage']}")
    if vadr_result["status"] == "ready":
        note_parts.append(f"VADR 注释目录 {vadr_result['output_dir']}")
    elif vadr_result["status"] == "failed":
        note_parts.append(vadr_result["note"])
    if snpeff_result["status"] == "ready":
        note_parts.append("已生成 snpEff 变异注释表")
    elif snpeff_result["status"] == "failed":
        note_parts.append(snpeff_result["note"])
    if resistance_result["status"] == "ready":
        note_parts.append("已生成耐药突变注释结果")
    elif resistance_result["status"] == "failed":
        note_parts.append(resistance_result["note"])
    note = "；".join(note_parts) if note_parts else "未获得稳定的 HA/NA 命中结果"

    _write_placeholder(
        pre,
        columns,
        [pre, "Influenza A virus", ha_hit["subtype"], na_hit["subtype"], subtype, ha_hit["subject"], na_hit["subject"], note],
    )
