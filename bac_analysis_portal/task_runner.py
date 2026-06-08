from __future__ import annotations

import csv
import os
import signal
import subprocess
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path

from .task_manager import AnalysisTaskManager, evaluate_species_review_gate, read_json, write_json
from .task_analytics import build_queue_analytics_snapshot, write_task_analytics_snapshot
from .store import PortalStore


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def main(task_file_arg: str) -> int:
    task_file = Path(task_file_arg).expanduser().resolve()
    task = read_json(task_file)
    if str(task.get("status") or "").upper() == "STOPPED" or str(task.get("requested_status") or "").upper() == "STOPPED":
        return 0
    log_path = Path(task["log_path"]).expanduser().resolve()
    log_path.parent.mkdir(parents=True, exist_ok=True)

    task["status"] = "RUNNING"
    task["started_at"] = task.get("started_at") or utc_now_iso()
    write_json(task_file, task)

    with log_path.open("a", encoding="utf-8") as log_handle:
        log_handle.write(f"[{utc_now_iso()}] 任务开始: {task['id']}\n")
        log_handle.flush()
        env = os.environ.copy()
        env["PYTHONUNBUFFERED"] = "1"
        for key, value in dict(task.get("env") or {}).items():
            if value is None:
                continue
            env[str(key)] = str(value)
        try:
            process = subprocess.Popen(
                task["command"],
                cwd=task["project_root"],
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                env=env,
                start_new_session=True,
            )
            task["pipeline_pid"] = process.pid
            write_json(task_file, task)
            while True:
                try:
                    exit_code = process.wait(timeout=5)
                    break
                except subprocess.TimeoutExpired:
                    _auto_pause_for_species_review(task_file, process, log_handle)
            latest = read_json(task_file)
            requested_status = str(latest.get("requested_status") or "").upper()
            if exit_code == 0:
                _postprocess_neisseria_amr(task, log_handle)
                _postprocess_tb_catalogue_and_lineage(task, log_handle)
            task["exit_code"] = exit_code
            task["status"] = "STOPPED" if requested_status == "STOPPED" else ("SUCCEEDED" if exit_code == 0 else "FAILED")
            task["finished_at"] = utc_now_iso()
            task["pipeline_pid"] = None
            task["runner_pid"] = None
            task["resume_pending"] = False
            task["requested_status"] = ""
            log_handle.write(f"\n[{utc_now_iso()}] 任务结束，退出码: {exit_code}\n")
            log_handle.flush()
            write_json(task_file, task)
            if task["status"] == "SUCCEEDED":
                _cache_task_analytics_snapshot(task_file, task, log_handle)
            _reconcile_after_task(task)
            return exit_code
        except Exception:
            task["status"] = "FAILED"
            task["finished_at"] = utc_now_iso()
            task["exit_code"] = -1
            task["pipeline_pid"] = None
            task["runner_pid"] = None
            task["resume_pending"] = False
            task["requested_status"] = ""
            log_handle.write(f"\n[{utc_now_iso()}] 调度器异常退出\n")
            log_handle.write(traceback.format_exc())
            log_handle.flush()
            write_json(task_file, task)
            _reconcile_after_task(task)
            return 1


def _cache_task_analytics_snapshot(task_file: Path, task: dict, log_handle) -> None:
    try:
        from .app import _build_report_payload

        report = _build_report_payload(task)
        snapshot = build_queue_analytics_snapshot(report)
        write_task_analytics_snapshot(task_file.parent, snapshot)
        log_handle.write(f"[{utc_now_iso()}] 已写入任务统计摘要缓存\n")
        log_handle.flush()
    except Exception:
        log_handle.write(f"[{utc_now_iso()}] 任务统计摘要缓存写入失败\n")
        log_handle.write(traceback.format_exc())
        log_handle.flush()


def _looks_like_neisseria_meningitidis(value: object) -> bool:
    text = str(value or "").strip().lower()
    return bool(text) and ("neisseria meningitidis" in text or "meningitidis" in text or "脑膜炎奈瑟" in text or "流脑" in text)


def _looks_like_mycobacterium_tuberculosis(value: object) -> bool:
    text = str(value or "").strip().lower()
    return bool(text) and (
        "mycobacterium tuberculosis" in text
        or "mycobacterium_tuberculosis" in text
        or "mycobacterium tuberculosis complex" in text
        or "m. tuberculosis" in text
        or "结核分枝杆菌" in text
        or "结核杆菌" in text
    )


def _read_first_tsv_row(path: Path) -> dict[str, str]:
    if not path.is_file():
        return {}
    try:
        with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            row = next(reader, None)
            return {str(k): str(v or "") for k, v in (row or {}).items()}
    except OSError:
        return {}


def _iter_report_dirs_for_samples(output_dir: Path):
    if not output_dir.is_dir():
        return
    seen: set[Path] = set()
    patterns = ("*.checkm.tsv", "*.mlst_Stat.txt", "*.fastp2.json")
    for pattern in patterns:
        for candidate in sorted(output_dir.rglob(pattern)):
            report_dir = candidate.parent.resolve()
            if report_dir in seen:
                continue
            seen.add(report_dir)
            sample_name = candidate.name
            for suffix in (".checkm.tsv", ".mlst_Stat.txt", ".fastp2.json"):
                if sample_name.endswith(suffix):
                    sample_name = sample_name.removesuffix(suffix)
                    break
            if sample_name:
                yield report_dir, sample_name


def _find_amr_gbk(report_dir: Path, sample_name: str) -> Path | None:
    for candidate in (
        report_dir / f"{sample_name}_prokka" / f"{sample_name}.gbk",
        report_dir / f"{sample_name}_prokka" / "main.gbk",
        report_dir / f"{sample_name}.gbk",
    ):
        if candidate.is_file():
            return candidate
    return None


def _postprocess_neisseria_amr(task: dict, log_handle) -> None:
    try:
        output_dir = Path(str((task.get("params") or {}).get("output_dir") or "")).expanduser().resolve()
    except OSError:
        return
    if not output_dir.is_dir():
        return

    project_root = Path(str(task.get("project_root") or "")).expanduser().resolve()
    script_path = project_root / "scripts" / "check_neisseria_meningitidis_amr_sites.py"
    site_table = project_root / "database" / "NM_mutate" / "neisseria_meningitidis_snp_amr_associations_literature_updated.csv"
    if not script_path.is_file() or not site_table.is_file():
        return

    for report_dir, sample_name in _iter_report_dirs_for_samples(output_dir):
        checkm_row = _read_first_tsv_row(report_dir / f"{sample_name}.checkm.tsv")
        if not any([
            _looks_like_neisseria_meningitidis(checkm_row.get("物种名称")),
            _looks_like_neisseria_meningitidis(checkm_row.get("species_name")),
            _looks_like_neisseria_meningitidis(checkm_row.get("mlst 物种名称")),
            _looks_like_neisseria_meningitidis(checkm_row.get("mlst_species_name")),
        ]):
            continue
        gbk_path = _find_amr_gbk(report_dir, sample_name)
        if gbk_path is None:
            continue
        output_path = report_dir / f"{sample_name}.neisseria_amr_calls.csv"
        command = [
            sys.executable,
            str(script_path),
            str(gbk_path),
            "--site-table",
            str(site_table),
            "--output",
            str(output_path),
        ]
        result = subprocess.run(
            command,
            cwd=str(project_root),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
        if result.returncode == 0:
            log_handle.write(f"[{utc_now_iso()}] 脑膜炎奈瑟耐药位点识别完成: {output_path}\n")
        else:
            log_handle.write(f"[{utc_now_iso()}] 脑膜炎奈瑟耐药位点识别失败: {sample_name}\n{result.stdout}\n")
        log_handle.flush()


def _postprocess_tb_catalogue_and_lineage(task: dict, log_handle) -> None:
    try:
        output_dir = Path(str((task.get("params") or {}).get("output_dir") or "")).expanduser().resolve()
    except OSError:
        return
    if not output_dir.is_dir():
        return
    project_root = Path(str(task.get("project_root") or "")).expanduser().resolve()
    script_path = project_root / "scripts" / "run_tb_sample_postprocess.py"
    if not script_path.is_file():
        return
    task_env = {str(key): str(value) for key, value in dict(task.get("env") or {}).items() if value is not None}
    env = os.environ.copy()
    env.update(task_env)
    threads = str((task.get("params") or {}).get("thread") or 8)
    long_type = str((task.get("params") or {}).get("long_type") or "").strip().lower()
    platform = "illumina"
    if "nano" in long_type:
        platform = "nanopore"
    elif "pacbio" in long_type or "hifi" in long_type:
        platform = "pacbio"
    for report_dir, sample_name in _iter_report_dirs_for_samples(output_dir):
        checkm_row = _read_first_tsv_row(report_dir / f"{sample_name}.checkm.tsv")
        if not any([
            _looks_like_mycobacterium_tuberculosis(checkm_row.get("物种名称")),
            _looks_like_mycobacterium_tuberculosis(checkm_row.get("species_name")),
            _looks_like_mycobacterium_tuberculosis(checkm_row.get("mlst 物种名称")),
            _looks_like_mycobacterium_tuberculosis(checkm_row.get("mlst_species_name")),
        ]):
            continue
        command = [
            sys.executable,
            str(script_path),
            "--report-dir",
            str(report_dir),
            "--sample",
            sample_name,
            "--project-root",
            str(project_root),
            "--threads",
            threads,
            "--platform",
            platform,
        ]
        result = subprocess.run(
            command,
            cwd=str(project_root),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            env=env,
            check=False,
        )
        if result.returncode == 0:
            log_handle.write(f"[{utc_now_iso()}] 结核分枝杆菌 H37Rv 有参 SNP 与 tb-profiler 分析完成: {sample_name}\n")
        else:
            log_handle.write(f"[{utc_now_iso()}] 结核分枝杆菌后处理失败: {sample_name}\n{result.stdout}\n")
        log_handle.flush()


def _reconcile_after_task(task: dict) -> None:
    try:
        project_root = Path(task["project_root"]).expanduser().resolve()
        store = PortalStore.from_project_root(project_root)
        max_concurrent = int(store.get_setting("max_concurrent_tasks", "2") or "2")
        from .task_manager import AnalysisTaskManager
        AnalysisTaskManager.from_project_root(project_root).reconcile_queue(max_concurrent)
    except Exception:
        return


def _signal_process_group(process: subprocess.Popen, sig: int) -> bool:
    try:
        os.killpg(os.getpgid(process.pid), sig)
        return True
    except ProcessLookupError:
        return False
    except OSError:
        try:
            os.kill(process.pid, sig)
            return True
        except OSError:
            return False


def _auto_pause_for_species_review(task_file: Path, process: subprocess.Popen, log_handle) -> None:
    latest = read_json(task_file)
    if str(latest.get("status") or "").upper() != "RUNNING":
        return
    result = evaluate_species_review_gate(latest, Path(str(latest.get("log_path") or "")).expanduser())
    if not result:
        return
    review_gate = latest.get("review_gate") or {}
    if result.get("decision") == "pass":
        evaluated = [str(item or "").strip() for item in (review_gate.get("evaluated_samples") or []) if str(item or "").strip()]
        sample_name = str(result.get("sample_name") or "").strip()
        if sample_name and sample_name not in evaluated:
            evaluated.append(sample_name)
            review_gate["evaluated_samples"] = evaluated
            latest["review_gate"] = review_gate
            write_json(task_file, latest)
        return
    if str(review_gate.get("state") or "").strip().lower() == "pending":
        return
    if not _signal_process_group(process, signal.SIGSTOP):
        return
    sample_name = str(result.get("sample_name") or "").strip()
    latest["status"] = "PAUSED"
    latest["review_gate"] = {
        **review_gate,
        "type": "species_identification",
        "state": "pending",
        "sample_name": sample_name,
        "summary": str(result.get("summary") or "").strip(),
        "reason_codes": list(result.get("reason_codes") or []),
        "evidence": list(result.get("evidence") or []),
        "classified_reads": result.get("classified_reads"),
        "bacteria_ratio": result.get("bacteria_ratio"),
        "dominant_species": result.get("dominant_species"),
        "secondary_species": result.get("secondary_species"),
        "triggered_at": utc_now_iso(),
        "approved_samples": list(review_gate.get("approved_samples") or []),
        "evaluated_samples": list(review_gate.get("evaluated_samples") or []),
    }
    write_json(task_file, latest)
    log_handle.write(f"\n[{utc_now_iso()}] 自动暂停：{latest['review_gate']['summary']}\n")
    for item in latest["review_gate"].get("evidence") or []:
        log_handle.write(f"[review] {item}\n")
    log_handle.flush()
    _reconcile_after_task(latest)


if __name__ == "__main__":
    import sys

    raise SystemExit(main(sys.argv[1]))
