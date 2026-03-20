from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from uuid import uuid4

ASM_METHOD_OPTIONS: dict[str, list[str]] = {
    "shortasm": ["spades", "masurca", "meta"],
    "longasm": ["flye", "miniasm", "wtdbg2", "canu", "unicycler", "raven"],
    "shortref": ["bwa"],
    "longref": ["minimap2"],
    "shortlongasm": ["unicycler"],
}


class ValidationError(ValueError):
    pass


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_json(path: Path, payload: dict[str, Any]) -> None:
    tmp_path = path.with_suffix(f"{path.suffix}.tmp")
    tmp_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    tmp_path.replace(path)


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


@dataclass
class AnalysisTaskManager:
    project_root: Path
    task_root: Path
    python_executable: str

    @classmethod
    def from_project_root(cls, project_root: Path) -> "AnalysisTaskManager":
        task_root = Path(os.environ.get("BAC_ANALYSIS_TASK_ROOT", project_root / "analysis_tasks")).expanduser()
        task_root.mkdir(parents=True, exist_ok=True)
        return cls(
            project_root=project_root,
            task_root=task_root,
            python_executable=os.environ.get("BAC_ANALYSIS_PYTHON", sys.executable),
        )

    def list_tasks(self, owner: str | None = None) -> list[dict[str, Any]]:
        tasks = []
        for path in self.task_root.glob("*/task.json"):
            task = read_json(path)
            if owner and task.get("owner") != owner:
                continue
            tasks.append(self._serialize_task(task, include_log=False))
        tasks.sort(key=lambda item: item.get("created_at", ""), reverse=True)
        return tasks

    def get_task(self, task_id: str, *, log_lines: int = 120, owner: str | None = None) -> dict[str, Any]:
        task_file = self.task_root / task_id / "task.json"
        if not task_file.is_file():
            raise KeyError(f"Task not found: {task_id}")
        task = read_json(task_file)
        if owner and task.get("owner") != owner:
            raise KeyError(f"Task not found: {task_id}")
        return self._serialize_task(task, include_log=True, log_lines=log_lines)

    def delete_task(self, task_id: str, *, owner: str | None = None) -> None:
        task_dir = self.task_root / task_id
        task_file = task_dir / "task.json"
        if not task_file.is_file():
            raise KeyError(f"Task not found: {task_id}")
        task = read_json(task_file)
        if owner and task.get("owner") != owner:
            raise KeyError(f"Task not found: {task_id}")
        shutil.rmtree(task_dir)

    def create_task(
        self,
        payload: dict[str, Any],
        *,
        owner: str,
        owner_group: str = "",
        pipeline_script: str,
        pipeline_python: str | None = None,
    ) -> dict[str, Any]:
        pipeline_script_path = Path(pipeline_script).expanduser()
        pipeline_script_path = (
            (self.project_root / pipeline_script_path).resolve()
            if not pipeline_script_path.is_absolute()
            else pipeline_script_path.resolve()
        )
        if not pipeline_script_path.is_file():
            raise ValidationError(f"分析脚本不存在: {pipeline_script_path}")

        task_id = datetime.now().strftime("%Y%m%d%H%M%S") + "_" + uuid4().hex[:8]
        task_dir = self.task_root / task_id
        task_dir.mkdir(parents=True, exist_ok=False)
        task_file = task_dir / "task.json"
        log_file = task_dir / "pipeline.log"

        params = self._normalize_payload(payload)
        pipeline_python_executable = self._resolve_runtime_python(pipeline_python)
        command = self._build_command(params, pipeline_script_path, pipeline_python_executable)
        task = {
            "id": task_id,
            "name": params["task_name"],
            "owner": owner,
            "owner_group": owner_group,
            "status": "QUEUED",
            "created_at": utc_now_iso(),
            "started_at": None,
            "finished_at": None,
            "exit_code": None,
            "runner_pid": None,
            "pipeline_pid": None,
            "params": params,
            "command": command,
            "project_root": str(self.project_root),
            "log_path": str(log_file),
            "pipeline_script": str(pipeline_script_path),
            "pipeline_python": pipeline_python_executable,
        }
        write_json(task_file, task)

        runner = subprocess.Popen(
            [self.python_executable, "-m", "bac_analysis_portal.task_runner", str(task_file)],
            cwd=str(self.project_root),
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            start_new_session=True,
        )
        task["runner_pid"] = runner.pid
        write_json(task_file, task)
        return self.get_task(task_id, owner=owner)

    def create_demo_task(self, *, owner: str, owner_group: str = "") -> dict[str, Any]:
        demo_input = (self.project_root / "demo_data" / "fastq").resolve()
        if not demo_input.exists():
            raise ValidationError(f"Demo 数据不存在: {demo_input}")

        task_id = datetime.now().strftime("%Y%m%d%H%M%S") + "_" + uuid4().hex[:8]
        task_dir = self.task_root / task_id
        task_dir.mkdir(parents=True, exist_ok=False)
        task_file = task_dir / "task.json"
        log_file = task_dir / "pipeline.log"
        output_dir = task_dir / "demo_output"
        output_dir.mkdir(parents=True, exist_ok=True)

        params = {
            "task_name": "demo_fastq",
            "input_path": str(demo_input),
            "inputtype": "fastq",
            "output_dir": str(output_dir),
            "thread": 10,
            "minlongfilt": "500",
            "Qfilt": "10",
            "barcodekit": "none",
            "method": "spades",
            "long_type": "Nanopore",
            "ref": "noref",
            "gtf": "nogtf",
            "genome_len": "4m",
            "asm_type": "shortasm",
            "polish_times": "1",
            "polish_soft": "medaka",
            "species": "False",
            "rmhost": "norm",
            "runflow": "All",
            "abun": "1",
            "rna": "0",
            "fake_pip": 1,
        }
        finished_at = utc_now_iso()
        task = {
            "id": task_id,
            "name": "demo_fastq",
            "owner": owner,
            "owner_group": owner_group,
            "status": "SUCCEEDED",
            "created_at": finished_at,
            "started_at": finished_at,
            "finished_at": finished_at,
            "exit_code": 0,
            "runner_pid": None,
            "pipeline_pid": None,
            "params": params,
            "command": ["demo_task", str(demo_input)],
            "project_root": str(self.project_root),
            "log_path": str(log_file),
            "pipeline_script": "demo_data/fastq",
        }
        log_file.write_text(
            f"[{finished_at}] Demo 任务已生成\n"
            f"[{finished_at}] 输入目录: {demo_input}\n"
            f"[{finished_at}] 状态: SUCCEEDED\n",
            encoding="utf-8",
        )
        write_json(task_file, task)
        return self.get_task(task_id, owner=owner)

    def _normalize_payload(self, payload: dict[str, Any]) -> dict[str, Any]:
        input_path = self._resolve_existing_path(payload.get("input_path"), "input_path")
        output_dir = self._resolve_output_path(payload.get("output_dir"))
        output_dir.mkdir(parents=True, exist_ok=True)
        asm_type = self._normalize_asm_type(payload.get("asm_type"))
        method = self._normalize_method(payload.get("method"), asm_type)

        return {
            "task_name": self._clean_str(payload.get("task_name")) or f"bac_assemble_{Path(input_path).stem}",
            "input_path": input_path,
            "inputtype": self._clean_str(payload.get("inputtype")) or "fastq",
            "output_dir": str(output_dir),
            "thread": self._clean_int(payload.get("thread"), default=10, minimum=1),
            "minlongfilt": self._clean_str(payload.get("minlongfilt")) or "500",
            "Qfilt": self._clean_str(payload.get("Qfilt")) or "10",
            "barcodekit": self._clean_str(payload.get("barcodekit")) or "none",
            "method": method,
            "long_type": self._clean_str(payload.get("long_type")) or "Nanopore",
            "ref": self._resolve_optional_path(payload.get("ref"), "noref"),
            "gtf": self._resolve_optional_path(payload.get("gtf"), "nogtf"),
            "genome_len": self._clean_str(payload.get("genome_len")) or "4m",
            "asm_type": asm_type,
            "polish_times": self._clean_str(payload.get("polish_times")) or "1",
            "polish_soft": self._clean_str(payload.get("polish_soft")) or "medaka",
            "species": self._clean_str(payload.get("species")) or "False",
            "rmhost": self._clean_str(payload.get("rmhost")) or "norm",
            "runflow": self._clean_str(payload.get("runflow")) or "All",
            "abun": self._clean_str(payload.get("abun")) or "1",
            "rna": self._clean_str(payload.get("rna")) or "0",
            "fake_pip": self._clean_int(payload.get("fake_pip"), default=0, minimum=0),
        }

    def _build_command(self, params: dict[str, Any], pipeline_script: Path, pipeline_python: str) -> list[str]:
        return [
            pipeline_python,
            "-u",
            str(pipeline_script),
            "--input",
            params["input_path"],
            "--inputtype",
            params["inputtype"],
            "--minlongfilt",
            params["minlongfilt"],
            "--Qfilt",
            params["Qfilt"],
            "--barcodekit",
            params["barcodekit"],
            "--thread",
            str(params["thread"]),
            "--output",
            params["output_dir"],
            "--fake_pip",
            str(params["fake_pip"]),
            "--method",
            params["method"],
            "--long_type",
            params["long_type"],
            "--ref",
            params["ref"],
            "--gtf",
            params["gtf"],
            "--genome_len",
            params["genome_len"],
            "--asm_type",
            params["asm_type"],
            "--polish_times",
            params["polish_times"],
            "--polish_soft",
            params["polish_soft"],
            "--species",
            params["species"],
            "--rmhost",
            params["rmhost"],
            "--runflow",
            params["runflow"],
            "--abun",
            params["abun"],
            "--rna",
            params["rna"],
        ]

    def _serialize_task(self, task: dict[str, Any], *, include_log: bool, log_lines: int = 120) -> dict[str, Any]:
        payload = {
            "id": task.get("id"),
            "name": task.get("name"),
            "owner": task.get("owner", ""),
            "owner_group": task.get("owner_group", ""),
            "status": task.get("status"),
            "created_at": task.get("created_at"),
            "started_at": task.get("started_at"),
            "finished_at": task.get("finished_at"),
            "exit_code": task.get("exit_code"),
            "runner_pid": task.get("runner_pid"),
            "pipeline_pid": task.get("pipeline_pid"),
            "params": task.get("params", {}),
            "command": task.get("command", []),
            "log_path": task.get("log_path"),
            "pipeline_script": task.get("pipeline_script", ""),
        }
        if include_log:
            payload["log_tail"] = self._read_log_tail(Path(task.get("log_path", "")), log_lines)
        return payload

    def _read_log_tail(self, path: Path, log_lines: int) -> str:
        if not path.is_file():
            return ""
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        return "\n".join(lines[-log_lines:])

    def _resolve_runtime_python(self, raw_value: str | None) -> str:
        candidate = self._clean_str(raw_value) or self.python_executable
        python_path = Path(candidate).expanduser()
        python_path = python_path.resolve() if python_path.is_absolute() else (self.project_root / python_path).resolve()
        if not python_path.is_file():
            raise ValidationError(f"运行环境 Python 不存在: {python_path}")
        return str(python_path)

    def _resolve_existing_path(self, raw_value: Any, field_name: str) -> str:
        text = self._clean_str(raw_value)
        if not text:
            raise ValidationError(f"{field_name} 不能为空")
        path = Path(text).expanduser()
        path = (self.project_root / path).resolve() if not path.is_absolute() else path.resolve()
        if not path.exists():
            raise ValidationError(f"{field_name} 不存在: {path}")
        return str(path)

    def _resolve_output_path(self, raw_value: Any) -> Path:
        text = self._clean_str(raw_value)
        if not text:
            raise ValidationError("output_dir 不能为空")
        path = Path(text).expanduser()
        return (self.project_root / path).resolve() if not path.is_absolute() else path.resolve()

    def _resolve_optional_path(self, raw_value: Any, placeholder: str) -> str:
        text = self._clean_str(raw_value)
        if not text:
            return placeholder
        path = Path(text).expanduser()
        path = (self.project_root / path).resolve() if not path.is_absolute() else path.resolve()
        if not path.exists():
            raise ValidationError(f"文件不存在: {path}")
        return str(path)

    def _clean_str(self, value: Any) -> str:
        return str(value).strip() if value is not None else ""

    def _normalize_asm_type(self, value: Any) -> str:
        text = self._clean_str(value) or "shortasm"
        if text == "longshortasm":
            text = "shortlongasm"
        if text not in ASM_METHOD_OPTIONS:
            raise ValidationError(f"不支持的组装类型: {text}")
        return text

    def _normalize_method(self, value: Any, asm_type: str) -> str:
        text = self._clean_str(value)
        options = ASM_METHOD_OPTIONS.get(asm_type, [])
        if not text:
            return options[0] if options else ""
        if text not in options:
            raise ValidationError(f"{asm_type} 不支持组装方法 {text}，可选: {', '.join(options)}")
        return text

    def _clean_int(self, value: Any, *, default: int, minimum: int) -> int:
        if value in (None, ""):
            return default
        try:
            parsed = int(value)
        except (TypeError, ValueError) as exc:
            raise ValidationError(f"整数参数格式错误: {value}") from exc
        if parsed < minimum:
            raise ValidationError(f"整数参数必须 >= {minimum}")
        return parsed
