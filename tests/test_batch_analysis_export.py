from __future__ import annotations

import json
import secrets
import zipfile

from bac_analysis_portal.application import create_app


def _login(client) -> dict[str, str]:
    client.get("/login")
    with client.session_transaction() as session:
        token = session["_csrf_token"]
    headers = {"X-CSRF-Token": token}
    response = client.post("/login", json={"username": "admin", "password": "admin123"}, headers=headers)
    assert response.status_code == 200
    token = secrets.token_urlsafe(32)
    with client.session_transaction() as session:
        session["_csrf_token"] = token
    return {"X-CSRF-Token": token}


def test_batch_analysis_export_zip_includes_manifest_and_missing_files(tmp_path, monkeypatch) -> None:
    monkeypatch.setenv("PORTAL_DB_PATH", str(tmp_path / "portal.sqlite3"))
    monkeypatch.setenv("BAC_ANALYSIS_TASK_ROOT", str(tmp_path / "tasks"))

    sample_dir = tmp_path / "outputs" / "fastq_analysis" / "S1"
    sample_dir.mkdir(parents=True)
    (sample_dir / "S1.final.fasta").write_text(">S1\nACGT\n", encoding="utf-8")
    (sample_dir / "S1_serotype_result.tsv").write_text("样本名称\t血清型\nS1\tK1\n", encoding="utf-8")

    task_dir = tmp_path / "tasks" / "task-001"
    task_dir.mkdir(parents=True)
    (task_dir / "task.json").write_text(
        json.dumps(
            {
                "id": "task-001",
                "name": "batch_test",
                "owner": "admin",
                "status": "SUCCEEDED",
                "created_at": "2026-07-03T00:00:00",
                "params": {"output_dir": str(tmp_path / "outputs")},
            }
        ),
        encoding="utf-8",
    )

    app = create_app({"TESTING": True, "SECRET_KEY": "test-secret"})
    client = app.test_client()
    headers = _login(client)
    response = client.post(
        "/api/tasks/batch-analysis-export",
        json={
            "samples": [{"task_id": "task-001", "sample_name": "S1"}],
            "artifact_types": ["fasta", "qc", "serotype"],
        },
        headers=headers,
    )

    assert response.status_code == 200
    assert response.mimetype == "application/zip"
    zip_path = tmp_path / "analysis_results.zip"
    zip_path.write_bytes(response.data)
    with zipfile.ZipFile(zip_path) as archive:
        names = set(archive.namelist())
        assert "fasta/S1.final.fasta" in names
        assert "serotype/S1_serotype_result.tsv" in names
        assert "manifest.tsv" in names
        assert "missing_files.tsv" in names
        missing = archive.read("missing_files.tsv").decode("utf-8-sig")
        assert "S1.fastp2.json" in missing
        assert "S1.fastp2.log" in missing


def test_batch_analysis_export_prefixes_shared_filenames(tmp_path, monkeypatch) -> None:
    monkeypatch.setenv("PORTAL_DB_PATH", str(tmp_path / "portal.sqlite3"))
    monkeypatch.setenv("BAC_ANALYSIS_TASK_ROOT", str(tmp_path / "tasks"))

    sample_dir = tmp_path / "outputs" / "fastq_analysis" / "S1"
    sample_dir.mkdir(parents=True)
    (sample_dir / "Assem_info.tsv").write_text("sample\tcontig\nS1\t1\n", encoding="utf-8")
    (sample_dir / "Assem_abricate_CARD.tsv").write_text("gene\nbla\n", encoding="utf-8")

    task_dir = tmp_path / "tasks" / "task-001"
    task_dir.mkdir(parents=True)
    (task_dir / "task.json").write_text(
        json.dumps(
            {
                "id": "task-001",
                "name": "batch_test",
                "owner": "admin",
                "status": "SUCCEEDED",
                "created_at": "2026-07-03T00:00:00",
                "params": {"output_dir": str(tmp_path / "outputs")},
            }
        ),
        encoding="utf-8",
    )

    app = create_app({"TESTING": True, "SECRET_KEY": "test-secret"})
    client = app.test_client()
    headers = _login(client)
    response = client.post(
        "/api/tasks/batch-analysis-export",
        json={
            "samples": [{"task_id": "task-001", "sample_name": "S1"}],
            "artifact_types": ["assembly", "resistance_virulence"],
        },
        headers=headers,
    )

    assert response.status_code == 200
    zip_path = tmp_path / "analysis_results.zip"
    zip_path.write_bytes(response.data)
    with zipfile.ZipFile(zip_path) as archive:
        names = set(archive.namelist())
        assert "assembly/S1_Assem_info.tsv" in names
        assert "resistance_virulence/S1_Assem_abricate_CARD.tsv" in names


def test_batch_analysis_export_prefixes_task_id_for_duplicate_sample_names(tmp_path, monkeypatch) -> None:
    monkeypatch.setenv("PORTAL_DB_PATH", str(tmp_path / "portal.sqlite3"))
    monkeypatch.setenv("BAC_ANALYSIS_TASK_ROOT", str(tmp_path / "tasks"))

    for task_id, output_name in (("task-a", "outputs_a"), ("task-b", "outputs_b")):
        sample_dir = tmp_path / output_name / "fastq_analysis" / "S1"
        sample_dir.mkdir(parents=True)
        (sample_dir / "Assem_info.tsv").write_text(f"sample\tcontig\n{task_id}\t1\n", encoding="utf-8")
        task_dir = tmp_path / "tasks" / task_id
        task_dir.mkdir(parents=True)
        (task_dir / "task.json").write_text(
            json.dumps(
                {
                    "id": task_id,
                    "name": task_id,
                    "owner": "admin",
                    "status": "SUCCEEDED",
                    "created_at": "2026-07-03T00:00:00",
                    "params": {"output_dir": str(tmp_path / output_name)},
                }
            ),
            encoding="utf-8",
        )

    app = create_app({"TESTING": True, "SECRET_KEY": "test-secret"})
    client = app.test_client()
    headers = _login(client)
    response = client.post(
        "/api/tasks/batch-analysis-export",
        json={
            "samples": [
                {"task_id": "task-a", "sample_name": "S1"},
                {"task_id": "task-b", "sample_name": "S1"},
            ],
            "artifact_types": ["assembly"],
        },
        headers=headers,
    )

    assert response.status_code == 200
    zip_path = tmp_path / "analysis_results.zip"
    zip_path.write_bytes(response.data)
    with zipfile.ZipFile(zip_path) as archive:
        names = set(archive.namelist())
        assert "assembly/task-a_S1_Assem_info.tsv" in names
        assert "assembly/task-b_S1_Assem_info.tsv" in names
        assert "assembly/S1_Assem_info.tsv" not in names
