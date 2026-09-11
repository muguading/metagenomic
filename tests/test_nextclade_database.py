from __future__ import annotations

import json
from pathlib import Path
from subprocess import CompletedProcess
from unittest.mock import Mock

import pytest

from conftest import create_portal_user, login_as

from bac_analysis_portal import nextclade_database
from bac_analysis_portal.nextclade_database import (
    NextcladeDatabaseManager,
    _download_and_replace,
    _match_remote_dataset,
    _resolve_dataset_root,
)
from bac_analysis_portal.task_manager import ValidationError


def _dataset_payload(
    *,
    name: str,
    accession: str,
    reference: str,
    tag: str,
    path: str = "",
    shortcuts: list[str] | None = None,
) -> dict[str, object]:
    return {
        "path": path,
        "shortcuts": shortcuts,
        "attributes": {
            "name": name,
            "reference accession": accession,
            "reference name": reference,
        },
        "version": {"tag": tag},
    }


def _write_dataset(root: Path, directory: str, payload: dict[str, object]) -> Path:
    target = root / directory
    target.mkdir(parents=True)
    (target / "pathogen.json").write_text(json.dumps(payload), encoding="utf-8")
    return target


def test_match_remote_dataset_prefers_exact_shortcut_when_reference_is_shared() -> None:
    local = _dataset_payload(
        name="Mpox virus (Clade IIb)",
        accession="NC_063383.1",
        reference="MPXV-M5312_HM12_Rivers",
        tag="old",
        shortcuts=["hMPXV"],
    )
    remote = [
        _dataset_payload(
            name="Mpox virus (All clades)",
            accession="NC_063383.1",
            reference="MPXV-M5312_HM12_Rivers",
            tag="new",
            path="nextstrain/mpox/all-clades",
            shortcuts=["MPXV"],
        ),
        _dataset_payload(
            name="Mpox virus (Clade IIb)",
            accession="NC_063383.1",
            reference="MPXV-M5312_HM12_Rivers",
            tag="new",
            path="nextstrain/mpox/clade-iib",
            shortcuts=["hMPXV"],
        ),
    ]

    matched = _match_remote_dataset(local, remote)

    assert matched is remote[1]


def test_resolve_dataset_root_supports_requested_and_legacy_layouts(tmp_path: Path) -> None:
    requested = tmp_path / "virus" / "nextclade_db"
    _write_dataset(requested, "sars-cov-2", {"version": {"tag": "v1"}})
    legacy = tmp_path / "nextclade_db"
    _write_dataset(legacy, "rsv_a", {"version": {"tag": "v1"}})

    assert _resolve_dataset_root(tmp_path) == requested.resolve()

    for child in requested.rglob("*"):
        if child.is_file():
            child.unlink()
    for child in sorted(requested.rglob("*"), reverse=True):
        if child.is_dir():
            child.rmdir()
    requested.rmdir()
    assert _resolve_dataset_root(tmp_path) == legacy.resolve()


def test_manager_check_reports_latest_and_outdated_datasets(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    dataset_root = tmp_path / "database" / "nextclade_db"
    _write_dataset(
        dataset_root,
        "sars-cov-2",
        _dataset_payload(
            name="SARS-CoV-2",
            accession="MN908947",
            reference="Wuhan-Hu-1/2019",
            tag="old",
            shortcuts=["sars-cov-2"],
        ),
    )
    _write_dataset(
        dataset_root,
        "rsv_a",
        _dataset_payload(
            name="RSV-A",
            accession="EPI_ISL_412866",
            reference="hRSV/A/England/397/2017",
            tag="same",
            shortcuts=["rsv_a"],
        ),
    )
    remote = [
        _dataset_payload(
            name="SARS-CoV-2",
            accession="MN908947",
            reference="Wuhan-Hu-1/2019",
            tag="new",
            path="nextstrain/sars-cov-2/wuhan-hu-1/orfs",
            shortcuts=["sars-cov-2"],
        ),
        _dataset_payload(
            name="RSV-A",
            accession="EPI_ISL_412866",
            reference="hRSV/A/England/397/2017",
            tag="same",
            path="nextstrain/rsv/a/EPI_ISL_412866",
            shortcuts=["rsv_a"],
        ),
    ]
    store = Mock()
    store.get_setting.side_effect = lambda key, default="": str(tmp_path) if key == "database_root" else default
    monkeypatch.setattr(nextclade_database, "_resolve_ncov_nextclade", lambda _store: Path("/fake/nextclade"))
    monkeypatch.setattr(nextclade_database, "_load_remote_datasets", lambda _executable: remote)
    monkeypatch.setattr(nextclade_database, "_nextclade_version", lambda _executable: "nextclade 3.test")

    result = NextcladeDatabaseManager(project_root=tmp_path, store=store).check()

    assert result["summary"] == {
        "total_count": 2,
        "latest_count": 1,
        "outdated_count": 1,
        "missing_count": 0,
        "unmatched_count": 0,
        "failed_count": 0,
    }
    assert [item["status"] for item in result["items"]] == ["latest", "outdated"]


def test_download_replaces_dataset_only_after_version_validation(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    dataset_root = tmp_path / "nextclade_db"
    target = _write_dataset(dataset_root, "virus-a", {"version": {"tag": "old"}})
    (target / "old.txt").write_text("old", encoding="utf-8")

    def fake_run(_executable: Path, arguments: list[str], *, timeout: int) -> CompletedProcess[str]:
        output_dir = Path(arguments[arguments.index("--output-dir") + 1])
        output_dir.mkdir()
        (output_dir / "pathogen.json").write_text(json.dumps({"version": {"tag": "new"}}), encoding="utf-8")
        (output_dir / "new.txt").write_text("new", encoding="utf-8")
        return CompletedProcess(["nextclade"], 0, "", "")

    monkeypatch.setattr(nextclade_database, "_run_nextclade", fake_run)

    _download_and_replace(
        Path("/fake/nextclade"),
        dataset_root,
        {"directory": "virus-a", "dataset_path": "remote/virus-a", "latest_tag": "new"},
    )

    assert (target / "new.txt").read_text(encoding="utf-8") == "new"
    assert not (target / "old.txt").exists()
    assert not list(dataset_root.glob(".*.backup-*"))


def test_download_keeps_original_when_downloaded_version_is_wrong(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    dataset_root = tmp_path / "nextclade_db"
    target = _write_dataset(dataset_root, "virus-a", {"version": {"tag": "old"}})
    (target / "old.txt").write_text("old", encoding="utf-8")

    def fake_run(_executable: Path, arguments: list[str], *, timeout: int) -> CompletedProcess[str]:
        output_dir = Path(arguments[arguments.index("--output-dir") + 1])
        output_dir.mkdir()
        (output_dir / "pathogen.json").write_text(json.dumps({"version": {"tag": "wrong"}}), encoding="utf-8")
        return CompletedProcess(["nextclade"], 0, "", "")

    monkeypatch.setattr(nextclade_database, "_run_nextclade", fake_run)

    with pytest.raises(ValidationError, match="下载版本校验失败"):
        _download_and_replace(
            Path("/fake/nextclade"),
            dataset_root,
            {"directory": "virus-a", "dataset_path": "remote/virus-a", "latest_tag": "new"},
        )

    assert (target / "old.txt").read_text(encoding="utf-8") == "old"


def test_nextclade_check_endpoint_is_admin_only(
    release_app_client, monkeypatch: pytest.MonkeyPatch
) -> None:
    app, client = release_app_client
    create_portal_user(app, username="analyst", password="analyst-pass", role="user", group_name="qa")
    user_headers = login_as(client, "analyst", "analyst-pass")

    forbidden = client.get("/api/admin/nextclade-datasets", headers=user_headers)
    client.post("/logout", headers=user_headers)

    monkeypatch.setattr(
        NextcladeDatabaseManager,
        "check",
        lambda self: {"items": [], "summary": {"total_count": 0}},
    )
    admin_headers = login_as(client)
    allowed = client.get("/api/admin/nextclade-datasets", headers=admin_headers)

    assert forbidden.status_code == 403
    assert allowed.status_code == 200
    assert allowed.get_json()["summary"]["total_count"] == 0


def test_scan_adds_missing_influenza_ha_and_na_datasets(tmp_path: Path) -> None:
    from bac_analysis_portal.nextclade_database import _scan_local_datasets

    dataset_root = tmp_path / "nextclade_db"
    dataset_root.mkdir()
    remote = [
        _dataset_payload(
            name="Influenza A H9 HA",
            accession="HA_REF",
            reference="H9 reference",
            tag="2026-08-01",
            path="nextstrain/flu/h9n2/ha/REF_H9",
            shortcuts=["flu_h9n2_ha"],
        ),
        _dataset_payload(
            name="Influenza A H9 NA",
            accession="NA_REF",
            reference="N2 reference",
            tag="2026-08-01",
            path="nextstrain/flu/h9n2/na/REF_N2",
            shortcuts=["flu_h9n2_na"],
        ),
    ]

    items = _scan_local_datasets(dataset_root, remote)

    assert [(item["display_name"], item["status"], item["updatable"]) for item in items] == [
        ("Influenza A H9 HA", "missing", True),
        ("Influenza A H9 NA", "missing", True),
    ]
    assert {item["dataset_path"] for item in items} == {
        "nextstrain/flu/h9n2/ha/REF_H9",
        "nextstrain/flu/h9n2/na/REF_N2",
    }


def test_download_creates_missing_dataset_directory(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    dataset_root = tmp_path / "nextclade_db"
    dataset_root.mkdir()

    def fake_run(_executable: Path, arguments: list[str], *, timeout: int) -> CompletedProcess[str]:
        output_dir = Path(arguments[arguments.index("--output-dir") + 1])
        output_dir.mkdir()
        (output_dir / "pathogen.json").write_text(json.dumps({"version": {"tag": "new"}}), encoding="utf-8")
        return CompletedProcess(["nextclade"], 0, "", "")

    monkeypatch.setattr(nextclade_database, "_run_nextclade", fake_run)
    _download_and_replace(
        Path("/fake/nextclade"),
        dataset_root,
        {"directory": "flu_h9n2_ha", "dataset_path": "nextstrain/flu/h9n2/ha", "latest_tag": "new"},
    )

    assert (dataset_root / "flu_h9n2_ha" / "pathogen.json").is_file()
