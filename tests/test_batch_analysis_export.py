from __future__ import annotations

import json
import secrets
import zipfile
from io import BytesIO

import pytest

from bac_analysis_portal.application import create_app
from bac_analysis_portal.export_utils import _build_xlsx_bytes


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


def test_sample_meta_preview_and_mapping_support_xlsx(tmp_path, monkeypatch) -> None:
    monkeypatch.setenv("PORTAL_DB_PATH", str(tmp_path / "portal.sqlite3"))
    monkeypatch.setenv("BAC_ANALYSIS_TASK_ROOT", str(tmp_path / "tasks"))
    app = create_app({"TESTING": True, "SECRET_KEY": "test-secret"})
    client = app.test_client()
    headers = _login(client)
    xlsx = _build_xlsx_bytes("meta", ["测序名称", "样本名称"], [["SEQ-001", "Case-001"], ["SEQ-002", "Case-002"]])

    preview = client.post(
        "/api/export/sample-meta/preview",
        data={"file": (BytesIO(xlsx), "sample_meta.xlsx")},
        content_type="multipart/form-data",
        headers=headers,
    )

    assert preview.status_code == 200
    preview_body = preview.get_json()
    assert preview_body["headers"] == ["测序名称", "样本名称"]
    assert preview_body["row_count"] == 2
    mapping = client.post(
        "/api/export/sample-meta/mapping",
        json={
            "import_id": preview_body["import_id"],
            "sequencing_column": "测序名称",
            "sample_column": "样本名称",
        },
        headers=headers,
    )

    assert mapping.status_code == 200
    assert mapping.get_json()["mapping"] == {"SEQ-001": "Case-001", "SEQ-002": "Case-002"}


def test_sample_meta_mapping_normalizes_illumina_library_suffixes(tmp_path, monkeypatch) -> None:
    monkeypatch.setenv("PORTAL_DB_PATH", str(tmp_path / "portal.sqlite3"))
    monkeypatch.setenv("BAC_ANALYSIS_TASK_ROOT", str(tmp_path / "tasks"))
    app = create_app({"TESTING": True, "SECRET_KEY": "test-secret"})
    client = app.test_client()
    headers = _login(client)
    xlsx = _build_xlsx_bytes(
        "meta",
        ["测序名称", "样本名称"],
        [["26-2452", "Case-2452"], ["26-2456", "Case-2456"]],
    )

    preview = client.post(
        "/api/export/sample-meta/preview",
        data={"file": (BytesIO(xlsx), "sample_meta.xlsx")},
        content_type="multipart/form-data",
        headers=headers,
    )
    mapping = client.post(
        "/api/export/sample-meta/mapping",
        json={
            "import_id": preview.get_json()["import_id"],
            "sequencing_column": "测序名称",
            "sample_column": "样本名称",
        },
        headers=headers,
    )

    assert mapping.status_code == 200
    body = mapping.get_json()
    assert body["mapping"] == {"26-2452": "Case-2452", "26-2456": "Case-2456"}
    assert body["normalized_mapping"] == {"26-2452": "Case-2452", "26-2456": "Case-2456"}
    assert body["library_suffix_normalized_count"] == 2
    assert body["ambiguous_normalized_sequencing_names"] == []


@pytest.mark.parametrize(
    ("task_id", "task_name", "species"),
    [
        ("task-ncov", "ncov_export_test", "SARS-CoV-2"),
        ("task-mpox", "hmpxv_export_test", "Monkeypox virus"),
    ],
)
def test_batch_analysis_export_renames_ncov_or_monkeypox_fasta_file_and_first_contig(tmp_path, monkeypatch, task_id, task_name, species) -> None:
    monkeypatch.setenv("PORTAL_DB_PATH", str(tmp_path / "portal.sqlite3"))
    monkeypatch.setenv("BAC_ANALYSIS_TASK_ROOT", str(tmp_path / "tasks"))

    sample_dir = tmp_path / "outputs" / "fastq_analysis" / "SEQ-001"
    sample_dir.mkdir(parents=True)
    source_fasta = sample_dir / "SEQ-001.final.fasta"
    source_fasta.write_text(">original_contig original-description\nACGT\n>second_contig\nTGCA\n", encoding="utf-8")
    task_dir = tmp_path / "tasks" / task_id
    task_dir.mkdir(parents=True)
    (task_dir / "task.json").write_text(
        json.dumps(
            {
                "id": task_id,
                "name": task_name,
                "owner": "admin",
                "status": "SUCCEEDED",
                "created_at": "2026-07-03T00:00:00",
                "params": {"output_dir": str(tmp_path / "outputs"), "species": species},
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
            "samples": [{"task_id": task_id, "sample_name": "SEQ-001", "export_sample_name": "Case-001"}],
            "artifact_types": ["fasta"],
        },
        headers=headers,
    )

    assert response.status_code == 200
    with zipfile.ZipFile(BytesIO(response.data)) as archive:
        assert "fasta/Case-001.final.fasta" in archive.namelist()
        assert archive.read("fasta/Case-001.final.fasta").decode("utf-8") == ">Case-001 original-description\nACGT\n>second_contig\nTGCA\n"
    assert source_fasta.read_text(encoding="utf-8").startswith(">original_contig")
