from __future__ import annotations

import json
import secrets
from pathlib import Path

import pytest

from bac_analysis_portal.application import create_app
from bac_analysis_portal.export_routes import _find_analysis_artifacts
from bac_analysis_portal.portal_helpers import task_report_availability
from bac_analysis_portal.report_sources import _resolve_report_source
from bac_analysis_portal.task_manager import AnalysisTaskManager, ValidationError


@pytest.fixture
def app_client(tmp_path, monkeypatch):
    monkeypatch.setenv("PORTAL_DB_PATH", str(tmp_path / "portal.sqlite3"))
    monkeypatch.setenv("BAC_ANALYSIS_TASK_ROOT", str(tmp_path / "tasks"))
    app = create_app({"TESTING": True, "SECRET_KEY": "test-secret"})
    client = app.test_client()
    return app, client


def _login(client) -> dict[str, str]:
    client.get("/login")
    with client.session_transaction() as session:
        token = session["_csrf_token"]
    response = client.post(
        "/login",
        json={"username": "admin", "password": "admin123"},
        headers={"X-CSRF-Token": token},
    )
    assert response.status_code == 200
    token = secrets.token_urlsafe(32)
    with client.session_transaction() as session:
        session["_csrf_token"] = token
    return {"X-CSRF-Token": token}


def test_login_and_workstation_pages_are_accessible(app_client) -> None:
    _app, client = app_client

    login_page = client.get("/login")
    assert login_page.status_code == 200

    headers = _login(client)
    workstation = client.get("/workstation?module=community", headers=headers)

    assert workstation.status_code == 200
    assert "快速提交分析任务" in workstation.get_data(as_text=True)
    assert "群落分析" in workstation.get_data(as_text=True)


def test_task_submission_rejects_missing_input_path(app_client, tmp_path) -> None:
    _app, client = app_client
    headers = _login(client)

    response = client.post(
        "/api/tasks",
        json={
            "workstation_key": "pathosource",
            "pipeline_script": "PathoSource.py",
            "task_name": "missing_input",
            "input_path": str(tmp_path / "does-not-exist"),
            "output_dir": str(tmp_path / "outputs"),
            "thread": 2,
        },
        headers=headers,
    )

    assert response.status_code == 400
    assert "input_path 不存在" in response.get_json()["error"]


def test_fastq_batch_generation_reports_empty_directory(app_client, tmp_path) -> None:
    _app, client = app_client
    headers = _login(client)
    empty_fastq_dir = tmp_path / "empty-fastq"
    empty_fastq_dir.mkdir()

    response = client.post(
        "/api/batch-inputs/fastq-single",
        json={"task_name": "empty_fastq", "short_left": str(empty_fastq_dir), "short_right": ""},
        headers=headers,
    )

    assert response.status_code == 400
    assert "未识别到可组成样本的 fastq 文件" in response.get_json()["error"]


def test_table_export_rejects_unsupported_format(app_client) -> None:
    _app, client = app_client
    headers = _login(client)

    response = client.post(
        "/api/export/table",
        json={"title": "结果表", "format": "pdf", "columns": ["样本"], "rows": [["S1"]]},
        headers=headers,
    )

    assert response.status_code == 400
    assert "仅支持导出 csv、tsv 或 xlsx" in response.get_json()["error"]


def test_task_parameter_validation_rejects_invalid_thread_value(tmp_path) -> None:
    input_dir = tmp_path / "inputs"
    input_dir.mkdir()
    manager = AnalysisTaskManager(project_root=tmp_path, task_root=tmp_path / "tasks", python_executable="python")

    with pytest.raises(ValidationError, match="整数参数格式错误"):
        manager._normalize_payload(
            {
                "workstation_key": "bacteria",
                "task_name": "bad_thread",
                "input_path": str(input_dir),
                "output_dir": str(tmp_path / "outputs"),
                "asm_type": "shortasm",
                "method": "spades",
                "thread": "many",
            }
        )


def test_task_output_path_rejects_non_ascii_characters(tmp_path) -> None:
    input_dir = tmp_path / "inputs"
    input_dir.mkdir()
    manager = AnalysisTaskManager(project_root=tmp_path, task_root=tmp_path / "tasks", python_executable="python")

    with pytest.raises(ValidationError, match="输出路径不能包含非 ASCII 字符"):
        manager._normalize_payload(
            {
                "workstation_key": "bacteria",
                "task_name": "中文任务",
                "input_path": str(input_dir),
                "output_dir": str(tmp_path / "outputs"),
                "asm_type": "shortasm",
                "method": "spades",
                "thread": 2,
            }
        )


def test_task_output_path_appends_safe_task_directory(tmp_path) -> None:
    input_dir = tmp_path / "inputs"
    input_dir.mkdir()
    manager = AnalysisTaskManager(project_root=tmp_path, task_root=tmp_path / "tasks", python_executable="python")

    params = manager._normalize_payload(
        {
            "workstation_key": "bacteria",
            "task_name": "sample_001",
            "input_path": str(input_dir),
            "output_dir": str(tmp_path / "outputs"),
            "asm_type": "shortasm",
            "method": "spades",
            "thread": 2,
        }
    )

    assert Path(params["output_dir"]) == tmp_path / "outputs" / "sample_001"


def test_monitor_snapshot_only_counts_supported_sequence_files(tmp_path) -> None:
    input_dir = tmp_path / "watch"
    input_dir.mkdir()
    (input_dir / "S1_R1.fastq.gz").write_bytes(b"@r1\nACGT\n+\n!!!!\n")
    (input_dir / "S1_R2.fq").write_bytes(b"@r2\nACGT\n+\n!!!!\n")
    (input_dir / "notes.txt").write_text("ignore me", encoding="utf-8")
    manager = AnalysisTaskManager(project_root=tmp_path, task_root=tmp_path / "tasks", python_executable="python")

    snapshot = manager._build_monitor_snapshot(input_dir)

    assert snapshot["file_count"] == 2
    assert "S1_R1.fastq.gz" in snapshot["signature"]
    assert "notes.txt" not in snapshot["signature"]


def test_report_availability_detects_dynamic_result_directory(tmp_path) -> None:
    report_dir = tmp_path / "outputs"
    report_dir.mkdir()
    (report_dir / "summary.tsv").write_text("sample\treads\nS1\t10\n", encoding="utf-8")

    availability = task_report_availability(
        {"id": "task-1", "status": "SUCCEEDED", "params": {"output_dir": str(report_dir)}}
    )

    assert availability["available"] is True
    assert availability["mode"] == "dynamic"


def test_report_source_explains_missing_result_files(tmp_path) -> None:
    fastq_root = tmp_path / "outputs" / "fastq_analysis"
    (fastq_root / "S1").mkdir(parents=True)

    source = _resolve_report_source(
        {"id": "task-1", "name": "batch", "status": "RUNNING", "params": {"output_dir": str(tmp_path / "outputs")}}
    )

    assert source["available"] is False
    assert source["mode"] == "multi"
    assert "还没有样本产出可展示结果文件" in source["reason"]


def test_analysis_artifact_lookup_reports_found_and_missing_files(tmp_path) -> None:
    sample_dir = tmp_path / "S1"
    sample_dir.mkdir()
    fasta = sample_dir / "S1.final.fasta"
    fasta.write_text(">S1\nACGT\n", encoding="utf-8")

    found, missing = _find_analysis_artifacts(sample_dir, "S1", "fasta")
    qc_found, qc_missing = _find_analysis_artifacts(sample_dir, "S1", "qc")

    assert found == [fasta.resolve()]
    assert missing == []
    assert qc_found == []
    assert qc_missing == ["S1.fastp2.json", "S1.fastp2.log"]


def test_queued_task_serialization_exposes_status_and_progress(tmp_path) -> None:
    task_root = tmp_path / "tasks"
    task_dir = task_root / "task-001"
    task_dir.mkdir(parents=True)
    (task_dir / "task.json").write_text(
        json.dumps(
            {
                "id": "task-001",
                "name": "queued_task",
                "status": "QUEUED",
                "created_at": "2026-07-08T00:00:00+00:00",
                "log_path": str(task_dir / "pipeline.log"),
                "params": {"output_dir": str(tmp_path / "outputs")},
            }
        ),
        encoding="utf-8",
    )
    manager = AnalysisTaskManager(project_root=tmp_path, task_root=task_root, python_executable="python")

    task = manager.get_task("task-001")

    assert task["status"] == "QUEUED"
    assert task["progress"]["label"] == "等待任务启动"
    assert task["log_tail"] == ""
