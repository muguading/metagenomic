from __future__ import annotations

import json
import sys
from pathlib import Path
from types import ModuleType

import pandas as pd
import pytest

from bac_analysis_portal.knowledge_base import load_knowledge_base_bundle

# PathoNet itself does not use pytaxonkit; stub the optional import so this unit
# test remains runnable in the lightweight web-test environment.
sys.modules.setdefault("pytaxonkit", ModuleType("pytaxonkit"))

from metagenomic_refactor.strain_typing import PathoNet, _load_pathonet_knowledge_base
from conftest import create_portal_user, login_as


ROOT = Path(__file__).resolve().parents[1]
KNOWLEDGE_BASE_BROWSER_SCRIPT = ROOT / "bac_analysis_portal" / "static" / "app.js"


def _write_rules(path: Path, *, serotypes: list[str], genes: list[str]) -> None:
    path.write_text(
        json.dumps(
            {
                "schema_version": "v1",
                "id": "test_pathonet_typing",
                "entries": [
                    {
                        "species": "vcholerae",
                        "serotype": serotypes,
                        "vfgene": genes,
                    }
                ],
            }
        ),
        encoding="utf-8",
    )


def test_default_pathonet_rules_are_registered_in_the_knowledge_base() -> None:
    rules = _load_pathonet_knowledge_base(
        ROOT / "database" / "knowledge_base" / "pathonet" / "pathonet_typing.json"
    )
    bundle = load_knowledge_base_bundle(str(ROOT))

    assert rules["vcholerae"] == {"serotype": ["O1", "O139"], "vfgene": ["ctxA", "ctxB"]}
    assert bundle["summary"]["pathonet_rule_count"] == len(rules)
    assert bundle["validation"]["pathonet_rules"] == []


def test_knowledge_base_page_exposes_pathonet_runtime_rules() -> None:
    script = KNOWLEDGE_BASE_BROWSER_SCRIPT.read_text(encoding="utf-8")

    assert 'key: "pathonet_rules", label: "PathoNet 运行规则"' in script
    assert 'summary.pathonet_rule_count' in script
    assert 'META_PATHONET_KNOWLEDGE_BASE' in script
    assert 'key: "vfgene", label: "重点毒力基因"' in script
    assert 'renderPathoNetRuleCrudToolbar(payload, canManagePathoNet)' in script
    assert 'openWorkbenchFormModal({' in script
    assert 'panelClass: "pathonet-rule-modal-panel"' in script
    assert 'formClass: "pathonet-rule-modal-form"' in script
    assert 'data-pathonet-modal-error' in script
    assert 'data-add-pathonet-rule' in script
    assert 'data-pathonet-edit' in script
    assert 'data-pathonet-delete' in script
    assert 'pathonet-rule-action-button is-edit' in script
    assert 'pathonet-rule-action-button is-delete' in script
    assert 'requestJson("/api/knowledge-base/pathonet-rules"' in script
    assert 'const canManagePathoNet = activeCollection === "pathonet_rules" && state.currentUser?.role === "admin";' in script
    assert 'confirmDangerAction({' in script
    assert 'renderPathoNetRuleCrudEditor' not in script


def test_pathonet_uses_user_runtime_rule_file(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    knowledge_base_path = tmp_path / "custom_pathonet_typing.json"
    _write_rules(knowledge_base_path, serotypes=["O777"], genes=["custom_vf"])
    monkeypatch.setenv("META_PATHONET_KNOWLEDGE_BASE", str(knowledge_base_path))
    monkeypatch.chdir(tmp_path)

    pd.DataFrame({"基因名称": ["custom_vf", "ctxA"]}).to_csv("sample.vfdb.tsv", sep="\t", index=False)
    pd.DataFrame({"血清型": ["O777"]}).to_csv("sample_serotype_result.tsv", sep="\t", index=False)

    PathoNet("sample", "vcholerae")
    result = pd.read_table("sample.pathonet_result.tsv")

    assert result.loc[0, "血清型"] == "O777(重点关注)"
    assert result.loc[0, "毒力基因"] == "custom_vf"


def test_pathonet_allows_user_rule_files_that_omit_unneeded_species(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    knowledge_base_path = tmp_path / "custom_pathonet_typing.json"
    _write_rules(knowledge_base_path, serotypes=["O777"], genes=["custom_vf"])
    monkeypatch.setenv("META_PATHONET_KNOWLEDGE_BASE", str(knowledge_base_path))
    monkeypatch.chdir(tmp_path)

    pd.DataFrame({"基因名称": ["hcp"]}).to_csv("sample.vfdb.tsv", sep="\t", index=False)
    pd.DataFrame({"血清型": ["HS:1"]}).to_csv("sample_serotype_result.tsv", sep="\t", index=False)

    PathoNet("sample", "campylobacter")
    result = pd.read_table("sample.pathonet_result.tsv")

    assert result.loc[0, "血清型"] == "HS:1"
    assert result.loc[0, "毒力基因"] == "-"


def test_pathonet_rejects_malformed_user_runtime_rule_file(tmp_path: Path) -> None:
    knowledge_base_path = tmp_path / "invalid_pathonet_typing.json"
    knowledge_base_path.write_text(
        json.dumps({"entries": [{"species": "vcholerae", "serotype": [], "vfgene": "ctxA"}]}),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="vfgene.*list of strings"):
        _load_pathonet_knowledge_base(knowledge_base_path)


def test_only_admin_can_save_pathonet_rules(release_app_client, tmp_path: Path) -> None:
    app, client = release_app_client
    database_root = tmp_path / "runtime_database"
    database_root.mkdir()
    app.config["PORTAL_STORE"].set_setting("database_root", str(database_root))
    payload = {
        "entries": [
            {
                "species": "custom_species",
                "serotype": ["S1", "S2"],
                "vfgene": ["vfA"],
            }
        ]
    }
    create_portal_user(app, username="analyst", password="analyst-pass")

    analyst_headers = login_as(client, "analyst", "analyst-pass")
    forbidden = client.put("/api/knowledge-base/pathonet-rules", json=payload, headers=analyst_headers)
    assert forbidden.status_code == 403

    admin_headers = login_as(client)
    response = client.put("/api/knowledge-base/pathonet-rules", json=payload, headers=admin_headers)

    assert response.status_code == 200
    bundle = response.get_json()
    assert bundle["collections"]["pathonet_rules"] == payload["entries"]
    assert bundle["summary"]["pathonet_rule_count"] == 1
    updated_payload = {
        "entries": [
            {
                "species": "custom_species",
                "serotype": ["S9"],
                "vfgene": ["vfB", "vfC"],
            }
        ]
    }
    updated = client.put("/api/knowledge-base/pathonet-rules", json=updated_payload, headers=admin_headers)
    assert updated.status_code == 200
    assert updated.get_json()["collections"]["pathonet_rules"] == updated_payload["entries"]

    deleted = client.put("/api/knowledge-base/pathonet-rules", json={"entries": []}, headers=admin_headers)
    assert deleted.status_code == 200
    assert deleted.get_json()["collections"]["pathonet_rules"] == []
    saved = json.loads(
        (database_root / "knowledge_base" / "pathonet" / "pathonet_typing.json").read_text(encoding="utf-8")
    )
    assert saved["entries"] == []
