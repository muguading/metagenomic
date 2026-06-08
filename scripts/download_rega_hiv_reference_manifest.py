#!/usr/bin/env python3
from __future__ import annotations

import csv
import re
from pathlib import Path


OUTPUT_DIR = Path("/Users/wuhhh/Desktop/徐老师/代码/metagenomic/database/virus/HIV/rega_reference_genomes")
MANIFEST_PATH = OUTPUT_DIR / "reference_manifest.tsv"


RAW_REFERENCES = {
    "Subtype A1": [
        "A1_UA_01_01UADN139_DQ823357",
        "A1_RU_00_RU00051_EF545108",
        "A1_UG_UG031_AB098330",
    ],
    "Subtype A2": [
        "A2_CD9797CDKTB48",
        "A2_CY9494CY01741",
    ],
    "Subtype B": [
        "B_US_1986_5019_86_AY835780",
        "B_US_90_US2_AY173953",
        "B_BR_03_BREPM1028_EF637053",
    ],
    "Subtype C": [
        "C_ZA_98_TV001_AY162223",
        "C_IN_x_VB39_EF694033",
        "C_IL_99_99ET1_AY255823",
    ],
    "Subtype D": [
        "D_UG_99_99UGG35093_AF484495",
        "D_UG_98_98UG57143_AF484514",
        "D_ZA_84_R2_AY773338",
    ],
    "Subtype F1": [
        "F1_ES_x_X1670_DQ979024",
        "F1_BR_89_BZ126_AY173957",
    ],
    "Subtype F2": [
        "F2_CM_95_MP255_AJ249236",
        "F2_CM_97_CM53657_AF377956",
    ],
    "Subtype G": [
        "G_NG_01_01NGPL0669_DQ168576",
        "G_ES_05_ES_EU786670",
        "G_GH_2003_GHNJ175_AB231893",
    ],
    "Subtype H": [
        "H_CF_90_056_AF005496",
        "H_BE_93_VI997_AF190128",
        "H_BE_93_VI991_AF190127",
    ],
    "Subtype J": [
        "J_SE_93_SE7887_AF082394",
        "J_SE_1994_SE7022_AF082395_1",
        "J_CD_97_J_97DC_KTB147_EF614151",
    ],
    "Subtype K": [
        "K_CD_97_EQTB11C_AJ249235",
        "K_CM_96_MP535_AJ249239",
    ],
    "CRF01_AE": ["01_AECF19090CF11697", "01_AETH90CM240"],
    "CRF02_AG": ["02_AGCM9797CMP807", "02_AGFR91DJ264"],
    "CRF03_AB": ["03_ABRU97KAL153", "03_ABRU98RU98001"],
    "CRF04_CPX": ["04_CPXG1R9197PVCH", "04_CPXGR9797PVMY"],
    "CRF05_DF": ["05_DFBEVI1310", "05_DFBE93VI961"],
    "CRF06_CPX": ["06_CPXAU96BFP90", "06_CPXML9595ML84"],
    "CRF07_BC": ["07_BCC1N9797CN001", "07_BCCN9898CN009"],
    "CRF08_BC": ["08_BCCN19797CNGX7F", "08_BCCN9898CN006"],
    "CRF09_CPX": ["09_cpxCI0000IC_10092AJ866553", "09_cpxGH9696GH2911AY093605"],
    "CRF10_CD": ["10_CDTZ19696TZBF071", "10_CDTZ9696TZBF110"],
    "CRF11_CPX": ["11_CPXF1R99MP1298", "11_CPXGRGR17"],
    "CRF12_BF": ["12_BFAR97A32989", "12_BFU1Y99URTR23"],
    "CRF13_CPX": ["13_CPXC1M961849", "13_CPXCM964164"],
    "CRF14_BG": ["14_BGE1S00X477", "14_BGES99X397"],
    "CRF18_CPX": ["18_cpxCM97CM53379AF377959", "18_cpxCU99CU14AY586541"],
    "CRF19_CPX": ["19_cpxCU99CU29AY588971", "19_cpxCU99CU7AY894994"],
    "CRF20_BG": ["20_BGCU99Cu103AY586545", "20_BGES99R77AY586544"],
    "CRF24_BG": ["24_BGCU03CB378AY900574", "24_BGCU03CB471AY900575"],
    "CRF25_CPX": ["25_cpxCM0606CM_BA_040EU693240", "25_cpxSA03J11451EU697908"],
    "CRF27_CPX": ["27_cpxCD9797CDKTB49AJ404325", "27_cpxFR0404CD_FR_KZSAM851091"],
    "CRF29_BF": ["29_BFBR01BREPM16704DQ085876", "29_BFBR99BREPM11948DQ085871"],
    "CRF31_BC": ["31_BCBR02110PAEF091932", "31_BCBR0404BR142AY727527"],
    "CRF35_AD": ["35_ADAF0505AF026EF158043", "35_ADAF0505AF094EF158040"],
    "CRF37_CPX": ["37_cpxCM0000CMNYU926EF116594", "37_cpxCM97CM53392AF377957"],
    "CRF39_BF": ["39_BFBR0303BRRJ103EU735534", "39_BFBR0303BRRJ327EU735536"],
    "CRF40_BF": ["40_BFBR0404BRRJ115EU735538", "40_BFBR0404BRSQ46EU735540"],
    "CRF42_BF": ["42_BF_LU_04_EU170142", "42_BF_LU_03_EU170151"],
    "CRF43_02G": ["43_02GSA03J11223EU697904", "43_02GSA03J11243EU697907"],
    "CRF47_BF": ["47BFES2008P1942GQ372987", "47BFES2008X24572FJ670529"],
}


ACCESSION_PATTERNS = [
    re.compile(r"[A-Z]{6}\d{9}(?:_\d+)?(?!\d)"),
    re.compile(r"[A-Z]{4}\d{8}(?:_\d+)?(?!\d)"),
    re.compile(r"[A-Z]{2}\d{6}(?:_\d+)?(?!\d)"),
    re.compile(r"[A-Z]\d{5}(?:_\d+)?(?!\d)"),
]


def extract_accession_candidates(label: str) -> list[str]:
    candidates: list[str] = []
    seen: set[str] = set()
    for pattern in ACCESSION_PATTERNS:
        for match in pattern.finditer(label):
            candidate = match.group(0)
            if "_" in candidate:
                candidate = candidate.split("_", 1)[0]
            if candidate not in seen:
                seen.add(candidate)
                candidates.append(candidate)
    # Prefer the rightmost candidate in the label.
    candidates.sort(key=label.rfind, reverse=True)
    return candidates


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    with MANIFEST_PATH.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["group", "label", "accession_candidates", "download_status", "note"])
        for group, labels in RAW_REFERENCES.items():
            for label in labels:
                candidates = extract_accession_candidates(label)
                writer.writerow(
                    [
                        group,
                        label,
                        ",".join(candidates),
                        "pending" if candidates else "unresolved",
                        "" if candidates else "No accession-like token found in REGA label",
                    ]
                )
    print(MANIFEST_PATH)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
