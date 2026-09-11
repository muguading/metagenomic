from __future__ import annotations

from pathlib import Path

import pytest

from conftest import login_as


pytestmark = pytest.mark.release


def _post_task(client, headers, payload):
    return client.post("/api/tasks", json=payload, headers=headers)


def test_release_workstation_pages_for_all_modules_are_accessible(release_app_client) -> None:
    _app, client = release_app_client
    headers = login_as(client)

    for module_key in ("bacteria", "virus", "metagenome", "community", "pathosource"):
        response = client.get(f"/workstation?module={module_key}", headers=headers)
        assert response.status_code == 200
        assert "快速提交分析任务" in response.get_data(as_text=True)


def test_release_task_submission_accepts_each_analysis_module(release_app_client, release_data) -> None:
    _app, client = release_app_client
    headers = login_as(client)
    data = release_data
    output_root = data["output_root"]

    cases = [
        {
            "workstation_key": "bacteria",
            "task_name": "release_bacteria",
            "input_path": str(data["input_dir"]),
            "output_dir": str(output_root),
            "asm_type": "shortasm",
            "method": "spades",
            "thread": 2,
        },
        {
            "workstation_key": "virus",
            "task_name": "release_virus",
            "input_path": str(data["fasta"]),
            "inputtype": "fasta",
            "output_dir": str(output_root),
            "asm_type": "shortref",
            "method": "bwa",
            "species": "SARS-CoV-2",
            "thread": 2,
        },
        {
            "workstation_key": "metagenome",
            "task_name": "release_metagenome",
            "input_path": str(data["input_dir"]),
            "output_dir": str(output_root),
            "asm_type": "shortasm",
            "method": "meta",
            "thread": 2,
        },
        {
            "workstation_key": "pathosource",
            "task_name": "release_pathosource",
            "input_path": str(data["input_dir"]),
            "output_dir": str(output_root),
            "species": "salmonella",
            "thread": 2,
        },
        {
            "workstation_key": "community",
            "task_name": "release_community",
            "input_path": str(data["input_dir"]),
            "output_dir": str(output_root),
            "community_metadata": str(data["metadata"]),
            "community_taxonomy": str(data["taxonomy"]),
            "community_group_column": "Group",
            "thread": 2,
        },
    ]

    for payload in cases:
        response = _post_task(client, headers, payload)
        assert response.status_code == 201, response.get_data(as_text=True)
        created = response.get_json()
        assert created["status"] == "QUEUED"
        assert created["params"]["task_name"] == payload["task_name"]
        assert Path(created["params"]["output_dir"]).name == payload["task_name"]


def test_release_task_submission_rejects_bad_enum_and_missing_reference(release_app_client, release_data) -> None:
    _app, client = release_app_client
    headers = login_as(client)
    payload = {
        "workstation_key": "bacteria",
        "task_name": "release_bad_enum",
        "input_path": str(release_data["input_dir"]),
        "output_dir": str(release_data["output_root"]),
        "asm_type": "shortasm",
        "method": "not-a-method",
        "ref": str(release_data["input_dir"] / "missing.fa"),
        "thread": 2,
    }

    response = _post_task(client, headers, payload)

    assert response.status_code == 400
    assert "不支持组装方法" in response.get_json()["error"]


def test_release_task_submission_rejects_nonexistent_input_path(release_app_client, release_data) -> None:
    _app, client = release_app_client
    headers = login_as(client)
    response = _post_task(
        client,
        headers,
        {
            "workstation_key": "bacteria",
            "task_name": "release_missing_input",
            "input_path": str(release_data["input_dir"] / "missing"),
            "output_dir": str(release_data["output_root"]),
            "asm_type": "shortasm",
            "method": "spades",
            "thread": 2,
        },
    )

    assert response.status_code == 400
    assert "input_path 不存在" in response.get_json()["error"]


def test_release_task_submission_rejects_excessive_thread_count(release_app_client, release_data) -> None:
    _app, client = release_app_client
    headers = login_as(client)
    response = _post_task(
        client,
        headers,
        {
            "workstation_key": "bacteria",
            "task_name": "release_too_many_threads",
            "input_path": str(release_data["input_dir"]),
            "output_dir": str(release_data["output_root"]),
            "asm_type": "shortasm",
            "method": "spades",
            "thread": 9999,
        },
    )

    assert response.status_code == 400
