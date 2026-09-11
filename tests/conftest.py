from __future__ import annotations

import json
import secrets
from pathlib import Path
from typing import Any

import pytest

from bac_analysis_portal.application import create_app
from bac_analysis_portal.app_services import APP_SERVICES_EXTENSION_KEY
from bac_analysis_portal.store import PortalStore
from bac_analysis_portal.task_manager import AnalysisTaskManager, utc_now_iso


@pytest.fixture
def release_app_client(tmp_path, monkeypatch):
    monkeypatch.setenv("PORTAL_DB_PATH", str(tmp_path / "portal.sqlite3"))
    monkeypatch.setenv("BAC_ANALYSIS_TASK_ROOT", str(tmp_path / "tasks"))

    def _no_start_runner(self: AnalysisTaskManager, task: dict[str, Any], task_file: Path) -> None:
        task["runner_pid"] = None
        task["status"] = "QUEUED"
        task_file.write_text(json.dumps(task, ensure_ascii=False, indent=2), encoding="utf-8")

    monkeypatch.setattr(AnalysisTaskManager, "_start_runner", _no_start_runner)
    app = create_app({"TESTING": True, "SECRET_KEY": "release-test-secret"})
    dummy_pipeline = tmp_path / "dummy_pipeline.py"
    dummy_pipeline.write_text("raise SystemExit(0)\n", encoding="utf-8")
    store = app.config["PORTAL_STORE"]
    assert isinstance(store, PortalStore)
    store.set_setting("workspace_root", str(tmp_path))
    store.set_setting("pipeline_script", str(dummy_pipeline))
    client = app.test_client()
    return app, client


@pytest.fixture
def release_data(tmp_path):
    input_dir = tmp_path / "inputs"
    input_dir.mkdir()
    fastq_r1 = input_dir / "S1_R1.fastq.gz"
    fastq_r2 = input_dir / "S1_R2.fastq.gz"
    fastq_r1.write_bytes(b"@S1/1\nACGT\n+\n!!!!\n")
    fastq_r2.write_bytes(b"@S1/2\nTGCA\n+\n!!!!\n")
    fasta = input_dir / "S1.final.fasta"
    fasta.write_text(">S1\nACGTACGT\n", encoding="utf-8")
    metadata = input_dir / "metadata.tsv"
    metadata.write_text("SampleID\tGroup\nS1\tCase\n", encoding="utf-8")
    taxonomy = input_dir / "taxonomy.tsv"
    taxonomy.write_text("Feature ID\tTaxon\nASV1\tBacteria;Firmicutes\n", encoding="utf-8")
    batch_tsv = input_dir / "batch.tsv"
    batch_tsv.write_text(
        f"样本名称\t三代数据\t二代数据左\t二代数据右\t物种信息\nS1\t\t{fastq_r1}\t{fastq_r2}\tEscherichia coli\n",
        encoding="utf-8",
    )
    return {
        "input_dir": input_dir,
        "fastq_r1": fastq_r1,
        "fastq_r2": fastq_r2,
        "fasta": fasta,
        "metadata": metadata,
        "taxonomy": taxonomy,
        "batch_tsv": batch_tsv,
        "output_root": tmp_path / "outputs",
    }


def login_as(client, username: str = "admin", password: str = "admin123") -> dict[str, str]:
    client.get("/login")
    with client.session_transaction() as session:
        token = session["_csrf_token"]
    response = client.post(
        "/login",
        json={"username": username, "password": password},
        headers={"X-CSRF-Token": token},
    )
    assert response.status_code == 200
    token = secrets.token_urlsafe(32)
    with client.session_transaction() as session:
        session["_csrf_token"] = token
    return {"X-CSRF-Token": token}


def create_portal_user(
    app,
    *,
    username: str,
    password: str = "password123",
    role: str = "user",
    group_name: str = "qa",
    allowed_modules: list[str] | None = None,
) -> dict[str, Any]:
    store = app.config["PORTAL_STORE"]
    assert isinstance(store, PortalStore)
    return store.create_user(
        username=username,
        password=password,
        role=role,
        group_name=group_name,
        allowed_modules=allowed_modules or ["bacteria"],
    )


def write_task_record(
    app,
    task_id: str,
    *,
    status: str = "SUCCEEDED",
    owner: str = "admin",
    owner_group: str = "",
    output_dir: Path | None = None,
    params: dict[str, Any] | None = None,
    log_text: str = "",
) -> Path:
    manager = app.extensions[APP_SERVICES_EXTENSION_KEY].task_lifecycle_service.task_manager
    task_dir = manager.task_root / task_id
    task_dir.mkdir(parents=True, exist_ok=True)
    log_path = task_dir / "pipeline.log"
    if log_text:
        log_path.write_text(log_text, encoding="utf-8")
    merged_params = {"output_dir": str(output_dir or task_dir / "outputs")}
    if params:
        merged_params.update(params)
    payload = {
        "id": task_id,
        "name": task_id,
        "owner": owner,
        "owner_group": owner_group,
        "status": status,
        "created_at": utc_now_iso(),
        "started_at": utc_now_iso() if status in {"RUNNING", "SUCCEEDED", "FAILED", "STOPPED"} else None,
        "finished_at": utc_now_iso() if status in {"SUCCEEDED", "FAILED", "STOPPED"} else None,
        "exit_code": 0 if status == "SUCCEEDED" else (1 if status == "FAILED" else None),
        "runner_pid": None,
        "pipeline_pid": 12345 if status in {"RUNNING", "PAUSED"} else None,
        "resume_pending": False,
        "requested_status": "",
        "params": merged_params,
        "command": [],
        "project_root": str(Path(app.config["PROJECT_ROOT"])),
        "log_path": str(log_path),
        "pipeline_script": str(Path(app.config["PROJECT_ROOT"]) / "Bac_assemble_260112_newformat.py"),
        "pipeline_python": "base",
        "conda_root": "",
    }
    (task_dir / "task.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    return task_dir
