from __future__ import annotations

import zipfile

import pytest

from conftest import login_as, write_task_record


pytestmark = pytest.mark.release


def test_release_succeeded_task_report_page_and_data_are_available(release_app_client, tmp_path) -> None:
    app, client = release_app_client
    headers = login_as(client)
    output_dir = tmp_path / "outputs"
    output_dir.mkdir()
    (output_dir / "summary.tsv").write_text("sample\treads\nS1\t10\n", encoding="utf-8")
    write_task_record(app, "release-report", status="SUCCEEDED", output_dir=output_dir)

    page = client.get("/tasks/release-report/result-page", headers=headers)
    data = client.get("/api/tasks/release-report/report-data", headers=headers)

    assert page.status_code == 200
    assert data.status_code == 200
    assert data.get_json()["task"]["id"] == "release-report"


def test_release_pending_report_data_returns_clear_reason_for_missing_results(release_app_client, tmp_path) -> None:
    app, client = release_app_client
    headers = login_as(client)
    output_dir = tmp_path / "outputs"
    (output_dir / "fastq_analysis" / "S1").mkdir(parents=True)
    write_task_record(app, "release-pending-report", status="RUNNING", output_dir=output_dir)

    task = client.get("/api/tasks/release-pending-report", headers=headers)

    assert task.status_code == 200
    assert task.get_json()["result_exists"] is False


def test_release_table_export_supports_csv_tsv_and_xlsx(release_app_client) -> None:
    _app, client = release_app_client
    headers = login_as(client)

    for export_format, mimetype in (
        ("csv", "text/csv"),
        ("tsv", "text/tab-separated-values"),
        ("xlsx", "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"),
    ):
        response = client.post(
            "/api/export/table",
            json={"title": "结果表", "format": export_format, "columns": ["样本", "Reads"], "rows": [["S1", "10"]]},
            headers=headers,
        )
        assert response.status_code == 200
        assert response.mimetype == mimetype
        assert response.data


def test_release_batch_zip_export_records_manifest_missing_files_and_duplicate_sample_names(
    release_app_client, tmp_path
) -> None:
    app, client = release_app_client
    headers = login_as(client)
    for task_id, root_name in (("release-export-a", "outputs-a"), ("release-export-b", "outputs-b")):
        sample_dir = tmp_path / root_name / "fastq_analysis" / "S1"
        sample_dir.mkdir(parents=True)
        (sample_dir / "S1.final.fasta").write_text(f">{task_id}\nACGT\n", encoding="utf-8")
        write_task_record(app, task_id, status="SUCCEEDED", output_dir=tmp_path / root_name)

    response = client.post(
        "/api/tasks/batch-analysis-export",
        json={
            "samples": [
                {"task_id": "release-export-a", "sample_name": "S1"},
                {"task_id": "release-export-b", "sample_name": "S1"},
            ],
            "artifact_types": ["fasta", "qc"],
        },
        headers=headers,
    )

    assert response.status_code == 200
    archive_path = tmp_path / "export.zip"
    archive_path.write_bytes(response.data)
    with zipfile.ZipFile(archive_path) as archive:
        names = set(archive.namelist())
        assert "fasta/release-export-a_S1.final.fasta" in names
        assert "fasta/release-export-b_S1.final.fasta" in names
        assert "manifest.tsv" in names
        assert "missing_files.tsv" in names
        missing = archive.read("missing_files.tsv").decode("utf-8-sig")
        assert "S1.fastp2.json" in missing
