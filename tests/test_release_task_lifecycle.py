from __future__ import annotations

import json

import pytest

from bac_analysis_portal.app_services import APP_SERVICES_EXTENSION_KEY
from bac_analysis_portal.task_manager import AnalysisTaskManager, read_json, write_json
from conftest import login_as, write_task_record


pytestmark = pytest.mark.release


def test_release_queued_task_can_be_stopped_through_lifecycle_api(release_app_client) -> None:
    app, client = release_app_client
    headers = login_as(client)
    write_task_record(app, "release-stop", status="QUEUED")

    response = client.post("/api/tasks/release-stop/stop", headers=headers)

    assert response.status_code == 200
    assert response.get_json()["status"] == "STOPPED"


def test_release_running_task_can_pause_and_resume_when_process_signals_succeed(
    release_app_client, monkeypatch
) -> None:
    app, client = release_app_client
    headers = login_as(client)
    services = app.extensions[APP_SERVICES_EXTENSION_KEY]
    manager = services.task_lifecycle_service.task_manager
    write_task_record(app, "release-pause-resume", status="RUNNING")

    monkeypatch.setattr(AnalysisTaskManager, "_signal_pipeline", lambda self, task, sig: True)
    monkeypatch.setattr(AnalysisTaskManager, "_resume_pipeline", lambda self, task: task.update({"status": "RUNNING", "resume_pending": False}) or True)

    paused = client.post("/api/tasks/release-pause-resume/pause", headers=headers)
    resumed = client.post("/api/tasks/release-pause-resume/resume", headers=headers)

    assert paused.status_code == 200
    assert paused.get_json()["status"] == "PAUSED"
    assert resumed.status_code == 200
    assert resumed.get_json()["status"] == "RUNNING"
    assert manager.get_task("release-pause-resume")["status"] == "RUNNING"


def test_release_queue_reconcile_respects_max_concurrency(release_app_client, tmp_path, monkeypatch) -> None:
    app, _client = release_app_client
    services = app.extensions[APP_SERVICES_EXTENSION_KEY]
    manager = services.task_lifecycle_service.task_manager
    write_task_record(app, "release-running", status="RUNNING")
    write_task_record(app, "release-queued-a", status="QUEUED")
    write_task_record(app, "release-queued-b", status="QUEUED")
    started: list[str] = []

    def record_start(self, task, task_file):
        started.append(str(task["id"]))
        task["status"] = "RUNNING"
        write_json(task_file, task)

    monkeypatch.setattr(AnalysisTaskManager, "_start_runner", record_start)
    manager.reconcile_queue(1)
    assert started == []

    running_file = manager.task_root / "release-running" / "task.json"
    running_task = read_json(running_file)
    running_task["status"] = "SUCCEEDED"
    write_json(running_file, running_task)
    manager.reconcile_queue(1)

    assert started == ["release-queued-a"]


def test_release_failed_task_exposes_user_readable_failure_diagnosis(release_app_client) -> None:
    app, client = release_app_client
    headers = login_as(client)
    write_task_record(
        app,
        "release-failed",
        status="FAILED",
        log_text="Traceback\nModuleNotFoundError: No module named 'Bio'\n",
    )

    response = client.get("/api/tasks/release-failed", headers=headers)

    assert response.status_code == 200
    diagnosis = response.get_json()["failure_diagnosis"]
    assert diagnosis["category"] == "dependency_missing"
    assert diagnosis["label"] == "运行环境依赖缺失"
    assert diagnosis["rerun_recommended"] is True
