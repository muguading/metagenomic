from __future__ import annotations

from pathlib import Path

from bac_analysis_portal import app as portal_app
from bac_analysis_portal.knowledge_base import load_knowledge_base_bundle


ROOT = Path(__file__).resolve().parents[1]


def test_hepatitis_typing_knowledge_base_links_hev_3a() -> None:
    bundle = load_knowledge_base_bundle(str(ROOT))

    assert bundle["summary"]["typing_rule_count"] >= 70
    rules = bundle["collections"]["typing_rules"]
    assert {rule.get("broad_type") for rule in rules if rule.get("level") == "broad_type"} == {"HAV", "HBV", "HCV", "HDV", "HEV"}
    hev_3a = next(
        rule
        for rule in rules
        if rule.get("broad_type") == "HEV"
        and rule.get("level") == "subtype"
        and str(rule.get("subtype") or rule.get("serotype") or "").lower() == "3a"
    )

    assert hev_3a["species"] == "Hepatitis E virus"
    assert "AF082843" in hev_3a["reference_accessions"]
    assert any(str(item.get("fasta_path") or "").endswith("AF082843.fasta") for item in hev_3a["references"])
    assert any(str(item.get("gff3_path") or "").endswith("AF082843.gff3") for item in hev_3a["references"])


def test_hepatovirus_report_summary_uses_broad_and_subtype_typing_rules() -> None:
    summary = portal_app._build_viral_serotype_knowledge_summary(
        str(ROOT),
        "Hepatitis E virus",
        {
            "mode": "hepatovirus_typing",
            "predicted_group": "HEV",
            "predicted_subtype": "3A",
            "predicted_clade": "3A",
        },
    )

    assert "HEV 大亚型" in summary["headline"]
    assert "HEV 3A 子亚型" in summary["headline"]
    assert len(summary["items"]) >= 2
    assert summary["items"][0]["serotype"] == "HEV"
    assert summary["items"][1]["serotype"].lower() == "3a"
    assert "AF082843" in summary["items"][1]["reference_accessions"]
