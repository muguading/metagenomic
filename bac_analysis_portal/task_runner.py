from __future__ import annotations

import os
import subprocess
import traceback
from datetime import datetime, timezone
from pathlib import Path

from .task_manager import read_json, write_json


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def main(task_file_arg: str) -> int:
    task_file = Path(task_file_arg).expanduser().resolve()
    task = read_json(task_file)
    log_path = Path(task["log_path"]).expanduser().resolve()
    log_path.parent.mkdir(parents=True, exist_ok=True)

    task["status"] = "RUNNING"
    task["started_at"] = utc_now_iso()
    write_json(task_file, task)

    with log_path.open("a", encoding="utf-8") as log_handle:
        log_handle.write(f"[{utc_now_iso()}] 任务开始: {task['id']}\n")
        log_handle.flush()
        env = os.environ.copy()
        env["PYTHONUNBUFFERED"] = "1"
        try:
            process = subprocess.Popen(
                task["command"],
                cwd=task["project_root"],
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                env=env,
            )
            task["pipeline_pid"] = process.pid
            write_json(task_file, task)
            exit_code = process.wait()
            task["exit_code"] = exit_code
            task["status"] = "SUCCEEDED" if exit_code == 0 else "FAILED"
            task["finished_at"] = utc_now_iso()
            log_handle.write(f"\n[{utc_now_iso()}] 任务结束，退出码: {exit_code}\n")
            log_handle.flush()
            write_json(task_file, task)
            return exit_code
        except Exception:
            task["status"] = "FAILED"
            task["finished_at"] = utc_now_iso()
            task["exit_code"] = -1
            log_handle.write(f"\n[{utc_now_iso()}] 调度器异常退出\n")
            log_handle.write(traceback.format_exc())
            log_handle.flush()
            write_json(task_file, task)
            return 1


if __name__ == "__main__":
    import sys

    raise SystemExit(main(sys.argv[1]))

