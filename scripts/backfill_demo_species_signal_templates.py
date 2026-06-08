from __future__ import annotations

import json
import sqlite3
from pathlib import Path

ROOT = Path("/Users/wuhhh/Desktop/徐老师/代码/metagenomic")
DB_PATH = ROOT / "bac_analysis_portal.sqlite3"

SIGNAL_MAP = {
    "Acinetobacter baumannii": {
        "resistance_gene_hits": "OXA-23、blaADC-25、aph(3')-Ia",
        "virulence_gene_hits": "ompA、bap、csuE",
        "resistance_mge_hits": "OXA-23（转座子）、aph(3')-Ia（整合子）",
        "virulence_mge_hits": "bap（质粒）",
    },
    "Escherichia coli": {
        "resistance_gene_hits": "blaCTX-M-15、aac(6')-Ib-cr、qnrS1",
        "virulence_gene_hits": "fimH、iutA、papC",
        "resistance_mge_hits": "blaCTX-M-15（质粒）、qnrS1（质粒）",
        "virulence_mge_hits": "iutA（质粒）",
    },
    "Klebsiella pneumoniae": {
        "resistance_gene_hits": "blaKPC-2、blaCTX-M-65、rmtB",
        "virulence_gene_hits": "rmpA、iucA、ybtS",
        "resistance_mge_hits": "blaKPC-2（质粒）、rmtB（质粒）",
        "virulence_mge_hits": "rmpA（质粒）、iucA（质粒）",
    },
    "Staphylococcus aureus": {
        "resistance_gene_hits": "mecA、ermC、tetK",
        "virulence_gene_hits": "clfA、hla、eta",
        "resistance_mge_hits": "ermC（质粒）、tetK（质粒）",
        "virulence_mge_hits": "eta（噬菌体）",
    },
    "Pseudomonas aeruginosa": {
        "resistance_gene_hits": "blaVIM-2、aacA4、fosA",
        "virulence_gene_hits": "exoU、lasB、toxA",
        "resistance_mge_hits": "blaVIM-2（整合子）、aacA4（整合子）",
        "virulence_mge_hits": "exoU（基因岛）",
    },
    "Neisseria meningitidis": {
        "resistance_gene_hits": "penA、rpoB、gyrA",
        "virulence_gene_hits": "fHbp、lbpA、lbpB",
        "resistance_mge_hits": "macB（RRR（复制））、mtrR（修复）",
        "virulence_mge_hits": "fHbp（重组）、lbpA（重组）、lbpB（防御）",
    },
    "Neisseria meningitidis(98%) noSpe(0.41%)": {
        "resistance_gene_hits": "penA、rpoB、gyrA",
        "virulence_gene_hits": "fHbp、lbpA、lbpB",
        "resistance_mge_hits": "macB（RRR（复制））、mtrR（修复）",
        "virulence_mge_hits": "fHbp（重组）、lbpA（重组）、lbpB（防御）",
    },
}


def main() -> None:
    conn = sqlite3.connect(DB_PATH)
    updated = 0
    for species_name, signals in SIGNAL_MAP.items():
        rows = conn.execute(
            """
            select sample_key, resistance_gene_hits, virulence_gene_hits, resistance_mge_hits, virulence_mge_hits
              from sample_library
             where species_name = ?
            """,
            (species_name,),
        ).fetchall()
        for sample_key, arg, vf, arg_mge, vf_mge in rows:
            next_values = {
                "resistance_gene_hits": arg or signals["resistance_gene_hits"],
                "virulence_gene_hits": vf or signals["virulence_gene_hits"],
                "resistance_mge_hits": arg_mge or signals["resistance_mge_hits"],
                "virulence_mge_hits": vf_mge or signals["virulence_mge_hits"],
            }
            conn.execute(
                """
                update sample_library
                   set resistance_gene_hits = ?,
                       virulence_gene_hits = ?,
                       resistance_mge_hits = ?,
                       virulence_mge_hits = ?,
                       updated_at = datetime('now')
                 where sample_key = ?
                """,
                (
                    next_values["resistance_gene_hits"],
                    next_values["virulence_gene_hits"],
                    next_values["resistance_mge_hits"],
                    next_values["virulence_mge_hits"],
                    sample_key,
                ),
            )
            updated += 1
    conn.commit()
    conn.close()
    print(json.dumps({"updated": updated}, ensure_ascii=False))


if __name__ == "__main__":
    main()
