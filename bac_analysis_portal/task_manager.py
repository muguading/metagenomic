from __future__ import annotations

import csv
import json
import os
import re
import shlex
import signal
import shutil
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from uuid import uuid4

from metagenomic_refactor.common import resolve_conda_env_name

from .task_analytics import read_task_analytics_snapshot

ASM_METHOD_OPTIONS: dict[str, list[str]] = {
    "shortasm": ["spades", "masurca", "meta"],
    "longasm": ["flye", "miniasm", "wtdbg2", "canu", "unicycler", "raven"],
    "shortref": ["bwa"],
    "longref": ["minimap2"],
    "shortlongasm": ["unicycler"],
}

ANALYSIS_TARGET_OPTIONS = {"bacteria", "virus"}

MONITOR_INPUT_EXTENSIONS = (
    ".fastq", ".fq", ".fastq.gz", ".fq.gz",
    ".fasta", ".fa", ".fna", ".fasta.gz", ".fa.gz", ".fna.gz",
)
SHORT_READ_R1_RE = re.compile(r"(?i)(?:^|[._-])(r?1)(?:$|[._-])")
SHORT_READ_R2_RE = re.compile(r"(?i)(?:^|[._-])(r?2)(?:$|[._-])")
AUTO_REVIEW_SUPPORTED_WORKSTATIONS = {"bacteria"}
AUTO_REVIEW_MIN_CLASSIFIED_READS = 100
AUTO_REVIEW_MIN_BACTERIA_RATIO = 40.0


class ValidationError(ValueError):
    pass


PROGRESS_LINE_RE = re.compile(
    r"task_step：(?P<step>\d+)/(?:\s*)?(?P<total_step>\d+)\s+样本进度：(?P<sample_index>\d+)/(?:\s*)?(?P<sample_total>\d+)\s+样本：(?P<sample>[^\t]+)\s+(?P<message>[^\t\r\n]+)"
)


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _is_demo_task_record(task: dict[str, Any]) -> bool:
    if bool(task.get("is_demo")):
        return True
    if str(task.get("demo_type") or "").strip():
        return True
    command = task.get("command") if isinstance(task.get("command"), list) else []
    if command and str(command[0] or "").strip() == "demo_task":
        return True
    task_name = str(task.get("name") or "").strip().lower()
    if task_name.startswith("demo_"):
        return True
    pipeline_script = str(task.get("pipeline_script") or "").strip().lower()
    return pipeline_script.startswith("demo_data/")


def default_runflow_for_selection(method: str, analysis_target: str = "bacteria") -> str:
    if str(analysis_target).strip() == "virus":
        return "基因组组装,物种鉴定,分型鉴定"
    return "基因组组装,病毒组装,物种鉴定,元件预测" if str(method).strip() == "meta" else "基因组组装,物种鉴定,耐药与毒力,mlst与血清型"


def _looks_like_influenza_task(params: dict[str, Any]) -> bool:
    if str(params.get("analysis_target") or "").strip().lower() != "virus":
        return False
    species = str(params.get("species") or "").strip().lower()
    return "influenza" in species or "流感" in species


def _resolve_conda_exe(conda_root: str | Path | None) -> str:
    raw_root = str(conda_root or "").strip()
    if raw_root:
        root = Path(raw_root).expanduser()
        candidates = (
            root / "bin" / "conda",
            root / "condabin" / "conda",
            root / "Scripts" / "conda.exe",
            root / "conda.exe",
        )
        for candidate in candidates:
            if candidate.is_file():
                return str(candidate)
        return str(candidates[0])
    return "conda"


def _conda_run_bin(conda_exe: str, env_name: str, executable: str) -> str:
    return " ".join(
        shlex.quote(part)
        for part in [conda_exe, "run", "-n", env_name, "--no-capture-output", executable]
    )


def _build_task_env(
    *,
    project_root: Path,
    params: dict[str, Any],
    database_root: str,
    conda_root: str = "",
    conda_env_map: dict[str, str],
) -> dict[str, str]:
    conda_exe = _resolve_conda_exe(conda_root)
    env = {
        "META_DATABASE_ROOT": str(database_root or ""),
        "META_CONDA_ROOT": str(conda_root or ""),
        "META_CONDA_EXE": conda_exe,
        "META_VIRSORTER2_BIN": _conda_run_bin(conda_exe, resolve_conda_env_name(conda_env_map.get("vfind", "genomad_aux")), "virsorter"),
        "META_CHECKV_BIN": _conda_run_bin(conda_exe, resolve_conda_env_name(conda_env_map.get("vfind", "genomad_aux")), "checkv"),
        "TB_PROFILER_BIN": _conda_run_bin(conda_exe, resolve_conda_env_name(conda_env_map.get("tb_profiler", "ncov")), "tb-profiler"),
        "TB_PROFILER_ENV": resolve_conda_env_name(conda_env_map.get("tb_profiler", "ncov")),
        "META_SISTR_HICAP_ENV": resolve_conda_env_name(conda_env_map.get("sistr_hicap", "sistr_hicap")),
        "COMMUNITY_QIIME_ENV": conda_env_map.get("qiime2", "qiime2-amplicon"),
        "COMMUNITY_MICROECO_ENV": conda_env_map.get("microeco", "microeco"),
    }
    if not _looks_like_influenza_task(params):
        return env

    vadr_root = (project_root / "soft").resolve()
    vadr_scripts_dir = vadr_root / "vadr"
    bio_easel_dir = vadr_root / "Bio-Easel"
    infernal_bin_dir = vadr_root / "infernal" / "binaries"
    sequip_dir = vadr_root / "sequip"
    blast_bin_dir = vadr_root / "ncbi-blast" / "bin"
    fasta_bin_dir = vadr_root / "fasta" / "bin"
    minimap2_dir = vadr_root / "minimap2"
    flu_model_dir = (vadr_root / "vadr-models-flu").resolve()
    if not flu_model_dir.is_dir():
        fallback_model_dir = (vadr_root / "vadr-models-flu-1.6.3-2").resolve()
        if fallback_model_dir.is_dir():
            flu_model_dir = fallback_model_dir
    required_paths = [
        vadr_scripts_dir / "v-annotate.pl",
        bio_easel_dir / "blib" / "lib",
        bio_easel_dir / "blib" / "arch",
        infernal_bin_dir / "cmalign",
        sequip_dir,
        blast_bin_dir / "blastn",
        fasta_bin_dir / "fasta36",
        minimap2_dir / "minimap2",
        flu_model_dir,
    ]
    if not all(path.exists() for path in required_paths):
        return env

    perl5lib_parts = [
        str(vadr_scripts_dir),
        str(sequip_dir),
        str((bio_easel_dir / "blib" / "lib").resolve()),
        str((bio_easel_dir / "blib" / "arch").resolve()),
    ]
    path_parts = [
        str(vadr_scripts_dir),
        str(blast_bin_dir),
        str(fasta_bin_dir),
        str(infernal_bin_dir),
        str(minimap2_dir),
    ]
    inherited_perl5lib = str(os.environ.get("PERL5LIB") or "").strip()
    inherited_path = str(os.environ.get("PATH") or "").strip()
    env.update(
        {
            "VADRINSTALLDIR": str(vadr_root),
            "VADRSCRIPTSDIR": str(vadr_scripts_dir),
            "VADRCONFIGFILE": str((vadr_scripts_dir / "vadr.config").resolve()),
            "VADRMODELDIR": str(flu_model_dir),
            "VADRINFERNALDIR": str(infernal_bin_dir),
            "VADREASELDIR": str(infernal_bin_dir),
            "VADRHMMERDIR": str(infernal_bin_dir),
            "VADRBIOEASELDIR": str(bio_easel_dir),
            "VADRSEQUIPDIR": str(sequip_dir),
            "VADRBLASTDIR": str(blast_bin_dir),
            "VADRFASTADIR": str(fasta_bin_dir),
            "VADRMINIMAP2DIR": str(minimap2_dir),
            "PERL5LIB": os.pathsep.join(perl5lib_parts + ([inherited_perl5lib] if inherited_perl5lib else [])),
            "PATH": os.pathsep.join(path_parts + ([inherited_path] if inherited_path else [])),
        }
    )
    return env


def write_json(path: Path, payload: dict[str, Any]) -> None:
    tmp_path = path.with_suffix(f"{path.suffix}.tmp")
    tmp_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    tmp_path.replace(path)


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _read_tsv_rows(path: Path) -> dict[str, list]:
    if not path.is_file():
        return {"columns": [], "rows": []}
    try:
        with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            rows = [row for row in reader if row]
    except OSError:
        return {"columns": [], "rows": []}
    if not rows:
        return {"columns": [], "rows": []}
    return {"columns": rows[0], "rows": rows[1:]}


def _safe_float(value: Any) -> float | None:
    try:
        if value in (None, ""):
            return None
        return float(str(value).strip().rstrip("%"))
    except (TypeError, ValueError):
        return None


def _safe_int(value: Any) -> int | None:
    try:
        if value in (None, ""):
            return None
        return int(float(str(value).strip()))
    except (TypeError, ValueError):
        return None


def _read_latest_progress_match(log_path: Path) -> re.Match[str] | None:
    if not log_path.is_file():
        return None
    latest_match: re.Match[str] | None = None
    for line in log_path.read_text(encoding="utf-8", errors="replace").splitlines():
        matched = PROGRESS_LINE_RE.search(line)
        if matched:
            latest_match = matched
    return latest_match


def _resolve_review_report_dir(task: dict[str, Any], sample_name: str) -> Path | None:
    params = task.get("params") or {}
    output_dir = str(params.get("output_dir") or "").strip()
    if not output_dir:
        return None
    output_root = Path(output_dir).expanduser()
    try:
        output_root = output_root.resolve()
    except OSError:
        return None
    fastq_analysis_root = output_root / "fastq_analysis"
    if sample_name and fastq_analysis_root.is_dir():
        sample_dir = fastq_analysis_root / sample_name
        if sample_dir.is_dir():
            return sample_dir
    if output_root.is_dir():
        return output_root
    return None


def _read_species_taxonomy_rows(path: Path) -> list[dict[str, Any]]:
    raw = _read_tsv_rows(path)
    if not raw["columns"] or not raw["rows"]:
        return []
    rows: list[dict[str, Any]] = []
    for row in raw["rows"]:
        record = {raw["columns"][index]: row[index] if index < len(row) else "" for index in range(len(raw["columns"]))}
        record["比例数值"] = _safe_float(record.get("比例")) or 0.0
        record["序列数量数值"] = _safe_int(record.get("序列数量")) or 0
        rows.append(record)
    return rows


def evaluate_species_review_gate(task: dict[str, Any], log_path: Path) -> dict[str, Any] | None:
    params = task.get("params") or {}
    workstation_key = str(params.get("workstation_key") or "").strip().lower()
    analysis_target = str(params.get("analysis_target") or "").strip().lower()
    method = str(params.get("method") or "").strip().lower()
    if workstation_key not in AUTO_REVIEW_SUPPORTED_WORKSTATIONS:
        return None
    if analysis_target != "bacteria" or method == "meta":
        return None

    latest_match = _read_latest_progress_match(log_path)
    if latest_match is None:
        return None
    step = int(latest_match.group("step"))
    if step < 2:
        return None
    sample_name = latest_match.group("sample").strip()
    if not sample_name:
        return None

    review_gate = task.get("review_gate") or {}
    approved_samples = {str(item or "").strip() for item in (review_gate.get("approved_samples") or []) if str(item or "").strip()}
    evaluated_samples = {str(item or "").strip() for item in (review_gate.get("evaluated_samples") or []) if str(item or "").strip()}
    if sample_name in approved_samples or sample_name in evaluated_samples:
        return None

    report_dir = _resolve_review_report_dir(task, sample_name)
    if report_dir is None:
        return None
    taxonomy_path = report_dir / f"{sample_name}_2.list.txt"
    rows = _read_species_taxonomy_rows(taxonomy_path)
    if not rows:
        return None

    valid_rows: list[dict[str, Any]] = []
    kingdom_reads: dict[str, int] = {"细菌": 0, "病毒": 0, "真菌": 0}
    total_reads = 0
    for row in rows:
        reads = int(row.get("序列数量数值") or 0)
        ratio = float(row.get("比例数值") or 0.0)
        species_name = str(row.get("种", "")).strip()
        genus_name = str(row.get("属", "")).strip() or (species_name.split()[0] if species_name else "")
        kingdom_name = str(row.get("界", "")).strip().lower()
        if "细菌" in kingdom_name or "bacteria" in kingdom_name or "eubacteria" in kingdom_name:
            kingdom_reads["细菌"] += reads
        elif "病毒" in kingdom_name or "virus" in kingdom_name or "viruses" in kingdom_name:
            kingdom_reads["病毒"] += reads
        elif "真菌" in kingdom_name or "fungi" in kingdom_name or "fungus" in kingdom_name or "mycota" in kingdom_name:
            kingdom_reads["真菌"] += reads
        total_reads += reads
        if not species_name or species_name == "-":
            continue
        valid_rows.append({
            "species": species_name,
            "genus": genus_name,
            "ratio": ratio,
            "reads": reads,
        })

    valid_rows.sort(key=lambda item: (item["reads"], item["ratio"]), reverse=True)
    if not valid_rows:
        return None

    top = valid_rows[0]
    second = valid_rows[1] if len(valid_rows) > 1 else None
    same_genus_rows = [row for row in valid_rows[1:] if row["genus"] and row["genus"] == top["genus"]]
    same_genus_count = len(same_genus_rows)
    same_genus_ratio = round(sum(row["ratio"] for row in same_genus_rows), 2)
    species_over_five = [row for row in valid_rows if row["ratio"] >= 5]
    top_two_ratio = round(top["ratio"] + (second["ratio"] if second else 0.0), 2)
    bacteria_ratio = round((kingdom_reads["细菌"] / total_reads) * 100, 2) if total_reads else 0.0
    non_bacteria_label = max(("病毒", "真菌"), key=lambda key: kingdom_reads[key])
    non_bacteria_ratio = round((kingdom_reads[non_bacteria_label] / total_reads) * 100, 2) if total_reads else 0.0

    same_genus_confusion = bool(second and second["genus"] == top["genus"] and same_genus_count > 0)
    low_bacteria = total_reads >= AUTO_REVIEW_MIN_CLASSIFIED_READS and bacteria_ratio < AUTO_REVIEW_MIN_BACTERIA_RATIO and non_bacteria_ratio > bacteria_ratio
    obvious_contamination = (
        total_reads >= AUTO_REVIEW_MIN_CLASSIFIED_READS
        and (top["ratio"] < 45 or len(species_over_five) >= 3 or (second and second["ratio"] >= 20 and top_two_ratio <= 85))
        and not same_genus_confusion
    )
    if not low_bacteria and not obvious_contamination:
        return {
            "sample_name": sample_name,
            "decision": "pass",
        }

    reason_codes: list[str] = []
    title_parts: list[str] = []
    evidence: list[str] = [
        f"当前样本：{sample_name}",
        f"主导物种：{top['species']}（{top['ratio']:.2f}%）",
        f"次高物种：{second['species']}（{second['ratio']:.2f}%）" if second else "次高物种：--",
        f"纳入统计的分类序列数：{total_reads}",
    ]
    if low_bacteria:
        reason_codes.append("low_bacteria")
        title_parts.append("细菌信号偏低")
        evidence.append(f"细菌分类占比约 {bacteria_ratio:.2f}% ，{non_bacteria_label}占比约 {non_bacteria_ratio:.2f}%")
    if obvious_contamination:
        reason_codes.append("contamination")
        title_parts.append("存在明显杂菌污染")
        evidence.append(f"占比 >=5% 的物种数为 {len(species_over_five)} 个，前两位物种合计占比约 {top_two_ratio:.2f}%")
    if same_genus_ratio:
        evidence.append(f"同属候选累计占比约 {same_genus_ratio:.2f}%")

    summary = "；".join(title_parts) + "。任务已在物种鉴定后暂停，等待人工确认是否继续后续组装与分型。"
    return {
        "sample_name": sample_name,
        "decision": "pause",
        "reason_codes": reason_codes,
        "summary": summary,
        "evidence": evidence,
        "bacteria_ratio": bacteria_ratio,
        "classified_reads": total_reads,
        "dominant_species": top["species"],
        "secondary_species": second["species"] if second else "",
    }


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
        self.refresh_monitored_tasks()
        tasks = []
        for path in self.task_root.glob("*/task.json"):
            task = read_json(path)
            if owner and task.get("owner") != owner:
                continue
            tasks.append(self._serialize_task(task, include_log=False))
        tasks.sort(key=lambda item: item.get("created_at", ""), reverse=True)
        return tasks

    def get_task(self, task_id: str, *, log_lines: int = 120, owner: str | None = None) -> dict[str, Any]:
        self.refresh_monitored_tasks()
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
        if not _is_demo_task_record(task):
            shutil.rmtree(task_dir)
            return

        output_dir_text = str((task.get("params") or {}).get("output_dir") or "").strip()
        preserved_output_dir: Path | None = None
        if output_dir_text:
            try:
                preserved_output_dir = Path(output_dir_text).expanduser().resolve()
            except OSError:
                preserved_output_dir = None
        for child in task_dir.iterdir():
            try:
                child_resolved = child.resolve()
            except OSError:
                child_resolved = child
            if preserved_output_dir is not None and child_resolved == preserved_output_dir:
                continue
            if child.is_dir():
                shutil.rmtree(child, ignore_errors=True)
            else:
                child.unlink(missing_ok=True)
        try:
            next(task_dir.iterdir())
        except StopIteration:
            task_dir.rmdir()

    def reconcile_queue(self, max_concurrent_tasks: int) -> None:
        limit = max(1, int(max_concurrent_tasks or 1))
        task_files = sorted(self.task_root.glob("*/task.json"), key=lambda item: item.parent.name)
        task_entries: list[tuple[Path, dict[str, Any]]] = []
        running_count = 0
        for task_file in task_files:
            try:
                task = read_json(task_file)
            except Exception:
                continue
            status = str(task.get("status") or "").upper()
            if status == "RUNNING":
                running_count += 1
            task_entries.append((task_file, task))

        available_slots = max(0, limit - running_count)
        if available_slots <= 0:
            return

        for task_file, task in sorted(task_entries, key=lambda item: item[1].get("created_at", "")):
            if available_slots <= 0:
                break
            status = str(task.get("status") or "").upper()
            if status != "QUEUED":
                continue
            if self._monitor_waiting(task):
                continue
            if task.get("resume_pending") and task.get("pipeline_pid"):
                if self._resume_pipeline(task):
                    write_json(task_file, task)
                    available_slots -= 1
                continue
            if task.get("runner_pid"):
                continue
            self._start_runner(task, task_file)
            available_slots -= 1

    def pause_task(self, task_id: str, *, owner: str | None = None, max_concurrent_tasks: int = 1) -> dict[str, Any]:
        task_file, task = self._load_task_for_change(task_id, owner=owner)
        status = str(task.get("status") or "").upper()
        if status != "RUNNING":
            raise ValidationError("只有运行中的任务可以暂停")
        if not self._signal_pipeline(task, signal.SIGSTOP):
            raise ValidationError("任务进程不存在，无法暂停")
        task["status"] = "PAUSED"
        write_json(task_file, task)
        self.reconcile_queue(max_concurrent_tasks)
        return self.get_task(task_id, owner=owner)

    def resume_task(self, task_id: str, *, owner: str | None = None, max_concurrent_tasks: int = 1) -> dict[str, Any]:
        task_file, task = self._load_task_for_change(task_id, owner=owner)
        status = str(task.get("status") or "").upper()
        if status != "PAUSED":
            raise ValidationError("只有已暂停的任务可以继续运行")
        review_gate = task.get("review_gate") or {}
        if str(review_gate.get("state") or "").strip().lower() == "pending":
            sample_name = str(review_gate.get("sample_name") or "").strip()
            approved_samples = [str(item or "").strip() for item in (review_gate.get("approved_samples") or []) if str(item or "").strip()]
            if sample_name and sample_name not in approved_samples:
                approved_samples.append(sample_name)
            review_gate["approved_samples"] = approved_samples
            review_gate["state"] = "approved"
            review_gate["approved_at"] = utc_now_iso()
            task["review_gate"] = review_gate
        task["status"] = "QUEUED"
        task["resume_pending"] = True
        write_json(task_file, task)
        self.reconcile_queue(max_concurrent_tasks)
        return self.get_task(task_id, owner=owner)

    def stop_task(self, task_id: str, *, owner: str | None = None, max_concurrent_tasks: int = 1) -> dict[str, Any]:
        task_file, task = self._load_task_for_change(task_id, owner=owner)
        status = str(task.get("status") or "").upper()
        if status == "QUEUED":
            task["status"] = "STOPPED"
            task["finished_at"] = utc_now_iso()
            task["exit_code"] = None
            task["runner_pid"] = None
            task["pipeline_pid"] = None
            task["resume_pending"] = False
            write_json(task_file, task)
            return self.get_task(task_id, owner=owner)
        if status not in {"RUNNING", "PAUSED"}:
            raise ValidationError("当前任务状态不支持停止")
        task["requested_status"] = "STOPPED"
        task["resume_pending"] = False
        write_json(task_file, task)
        if not self._signal_pipeline(task, signal.SIGTERM):
            task["status"] = "STOPPED"
            task["finished_at"] = utc_now_iso()
            task["runner_pid"] = None
            task["pipeline_pid"] = None
            write_json(task_file, task)
        self.reconcile_queue(max_concurrent_tasks)
        return self.get_task(task_id, owner=owner)

    def rerun_task(
        self,
        task_id: str,
        payload: dict[str, Any],
        *,
        owner: str | None = None,
        owner_group: str = "",
        pipeline_script: str,
        pipeline_python: str | None = None,
        max_concurrent_tasks: int = 1,
        database_root: str = "",
        conda_root: str = "",
        conda_envs: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        task_file, task = self._load_task_for_change(task_id, owner=owner)
        status = str(task.get("status") or "").upper()
        if status in {"RUNNING", "PAUSED", "QUEUED"}:
            raise ValidationError("当前任务仍在执行或等待中，不能重新运行")
        pipeline_script_path = Path(pipeline_script).expanduser()
        pipeline_script_path = (
            (self.project_root / pipeline_script_path).resolve()
            if not pipeline_script_path.is_absolute()
            else pipeline_script_path.resolve()
        )
        if not pipeline_script_path.is_file():
            raise ValidationError(f"分析脚本不存在: {pipeline_script_path}")
        params = self._normalize_payload(payload)
        pipeline_runtime_env = self._resolve_runtime_python(pipeline_python)
        conda_env_map = {str(key): resolve_conda_env_name(str(value).strip()) for key, value in dict(conda_envs or {}).items() if str(value).strip()}
        conda_exe = _resolve_conda_exe(conda_root)
        command = self._build_command(params, pipeline_script_path, pipeline_runtime_env, conda_exe)
        task["name"] = params["task_name"]
        task["owner_group"] = owner_group or str(task.get("owner_group") or "")
        task["status"] = "QUEUED"
        task["started_at"] = None
        task["finished_at"] = None
        task["exit_code"] = None
        task["runner_pid"] = None
        task["pipeline_pid"] = None
        task["resume_pending"] = False
        task["requested_status"] = ""
        task["params"] = params
        task["command"] = command
        task["env"] = _build_task_env(
            project_root=self.project_root,
            params=params,
            database_root=str(database_root or ""),
            conda_root=str(conda_root or ""),
            conda_env_map=conda_env_map,
        )
        task["pipeline_script"] = str(pipeline_script_path)
        task["pipeline_python"] = pipeline_runtime_env
        task["conda_root"] = str(conda_root or "")
        task["rerun_count"] = int(task.get("rerun_count") or 0) + 1
        task["rerun_at"] = utc_now_iso()
        task.pop("review_gate", None)
        write_json(task_file, task)
        self.reconcile_queue(max_concurrent_tasks)
        return self.get_task(task_id, owner=owner)

    def create_task(
        self,
        payload: dict[str, Any],
        *,
        owner: str,
        owner_group: str = "",
        pipeline_script: str,
        pipeline_python: str | None = None,
        max_concurrent_tasks: int = 1,
        database_root: str = "",
        conda_root: str = "",
        conda_envs: dict[str, Any] | None = None,
        extra_task_fields: dict[str, Any] | None = None,
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
        pipeline_runtime_env = self._resolve_runtime_python(pipeline_python)
        conda_env_map = {str(key): resolve_conda_env_name(str(value).strip()) for key, value in dict(conda_envs or {}).items() if str(value).strip()}
        conda_exe = _resolve_conda_exe(conda_root)
        command = self._build_command(params, pipeline_script_path, pipeline_runtime_env, conda_exe)
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
            "resume_pending": False,
            "requested_status": "",
            "params": params,
            "command": command,
            "env": _build_task_env(
                project_root=self.project_root,
                params=params,
                database_root=str(database_root or ""),
                conda_root=str(conda_root or ""),
                conda_env_map=conda_env_map,
            ),
            "project_root": str(self.project_root),
            "log_path": str(log_file),
            "pipeline_script": str(pipeline_script_path),
            "pipeline_python": pipeline_runtime_env,
            "conda_root": str(conda_root or ""),
        }
        if str(params.get("watch_mode") or "0") == "1":
            task["watch"] = self._build_initial_watch_state(params)
        if isinstance(extra_task_fields, dict):
            for key, value in extra_task_fields.items():
                if key in {"id", "params", "command", "log_path"}:
                    continue
                task[key] = value
        write_json(task_file, task)
        self.reconcile_queue(max_concurrent_tasks)
        return self.get_task(task_id, owner=owner)

    def _resolve_demo_workspace_root(
        self,
        workspace_root: Path,
        database_root: Path | None = None,
    ) -> Path:
        candidates: list[Path] = []

        def add_candidate(path: Path | None) -> None:
            if path is None:
                return
            resolved = path.expanduser().resolve()
            if resolved not in candidates:
                candidates.append(resolved)

        add_candidate(workspace_root)
        add_candidate(database_root)
        if database_root and database_root.name == "database":
            add_candidate(database_root.parent)
        add_candidate(self.project_root)

        for candidate in candidates:
            if (candidate / "demo_data").exists():
                return candidate
        return workspace_root

    def create_demo_task(
        self,
        *,
        owner: str,
        owner_group: str = "",
        demo_type: str = "fastq",
        workspace_root: str | Path | None = None,
        database_root: str | Path | None = None,
    ) -> dict[str, Any]:
        demo_type = str(demo_type or "fastq").strip().lower()
        workspace_root_path = Path(workspace_root).expanduser().resolve() if workspace_root else self.project_root
        database_root_path = Path(database_root).expanduser().resolve() if database_root else None
        workspace_root_path = self._resolve_demo_workspace_root(workspace_root_path, database_root_path)
        demo_definitions = {
            "fastq": {
                "task_name": "demo_fastq",
                "input_path": workspace_root_path / "demo_data" / "fastq",
                "output_path": workspace_root_path / "demo_data" / "fastq",
                "analysis_target": "bacteria",
                "inputtype": "fastq",
                "method": "spades",
                "asm_type": "shortasm",
                "runflow": default_runflow_for_selection("spades", "bacteria"),
                "pipeline_script": "demo_data/fastq",
            },
            "tb": {
                "task_name": "demo_tb",
                "input_path": workspace_root_path / "demo_data" / "tb_demo",
                "output_path": workspace_root_path / "demo_data" / "tb_demo",
                "analysis_target": "bacteria",
                "inputtype": "fastq",
                "method": "spades",
                "asm_type": "shortasm",
                "runflow": default_runflow_for_selection("spades", "bacteria"),
                "pipeline_script": "demo_data/tb_demo",
            },
            "meta": {
                "task_name": "demo_meta",
                "input_path": workspace_root_path / "demo_data" / "meta_1",
                "output_path": workspace_root_path / "demo_data" / "meta_1",
                "analysis_target": "bacteria",
                "inputtype": "fastq",
                "method": "meta",
                "asm_type": "shortasm",
                "runflow": default_runflow_for_selection("meta", "bacteria"),
                "pipeline_script": "demo_data/meta",
            },
            "ncov": {
                "task_name": "demo_ncov",
                "input_path": workspace_root_path / "demo_data" / "ncov" / "test.final.fasta",
                "output_path": workspace_root_path / "demo_data" / "ncov",
                "analysis_target": "virus",
                "inputtype": "fasta",
                "method": "spades",
                "asm_type": "shortref",
                "runflow": "分型鉴定",
                "pipeline_script": "demo_data/ncov",
                "workstation_key": "virus",
                "species": "SARS-CoV-2",
                "ref": workspace_root_path / "database" / "virus" / "ncov" / "ref.fna",
            },
            "flu": {
                "task_name": "demo_flu",
                "input_path": workspace_root_path / "demo_data" / "influenza" / "flu_samples.tsv",
                "output_path": workspace_root_path / "demo_data" / "flu_demo",
                "analysis_target": "virus",
                "inputtype": "fastq",
                "method": "freebayes",
                "asm_type": "shortref",
                "runflow": "物种鉴定,基因组组装,分型鉴定",
                "pipeline_script": "demo_data/flu_demo",
                "workstation_key": "virus",
                "species": "Influenza virus",
            },
            "rsv": {
                "task_name": "demo_rsv",
                "input_path": workspace_root_path / "demo_data" / "rsv_demo" / "ngs",
                "output_path": workspace_root_path / "demo_data" / "rsv_demo",
                "analysis_target": "virus",
                "inputtype": "fastq",
                "method": "freebayes",
                "asm_type": "shortref",
                "runflow": "物种鉴定,基因组组装,分型鉴定",
                "pipeline_script": "demo_data/rsv_demo",
                "workstation_key": "virus",
                "species": "Respiratory syncytial virus",
            },
            "hmpv": {
                "task_name": "demo_hmpv",
                "input_path": workspace_root_path / "demo_data" / "hmpv_demo",
                "output_path": workspace_root_path / "demo_data" / "hmpv_demo",
                "analysis_target": "virus",
                "inputtype": "fastq",
                "method": "freebayes",
                "asm_type": "shortref",
                "runflow": "物种鉴定,基因组组装,分型鉴定",
                "pipeline_script": "demo_data/hmpv_demo",
                "workstation_key": "virus",
                "species": "Human metapneumovirus",
            },
            "denv": {
                "task_name": "demo_denv",
                "input_path": workspace_root_path / "demo_data" / "denv_demo",
                "output_path": workspace_root_path / "demo_data" / "denv_demo",
                "analysis_target": "virus",
                "inputtype": "fastq",
                "method": "freebayes",
                "asm_type": "shortref",
                "runflow": "物种鉴定,基因组组装,分型鉴定",
                "pipeline_script": "demo_data/denv_demo",
                "workstation_key": "virus",
                "species": "Dengue virus",
            },
            "zikav": {
                "task_name": "demo_zikav",
                "input_path": workspace_root_path / "demo_data" / "zikav_demo_1" / "zikav_demo_1.final.fasta",
                "output_path": workspace_root_path / "demo_data" / "zikav_demo_1",
                "analysis_target": "virus",
                "inputtype": "fasta",
                "method": "spades",
                "asm_type": "shortref",
                "runflow": "分型鉴定",
                "pipeline_script": "demo_data/zikav_demo_1",
                "workstation_key": "virus",
                "species": "Zika virus",
            },
            "chikv": {
                "task_name": "demo_chikv",
                "input_path": workspace_root_path / "database" / "nextclade_db" / "chikv" / "reference.fasta",
                "output_path": workspace_root_path / "demo_data" / "chikv_demo_1",
                "analysis_target": "virus",
                "inputtype": "fasta",
                "method": "spades",
                "asm_type": "shortref",
                "runflow": "分型鉴定",
                "pipeline_script": "demo_data/chikv_demo_1",
                "workstation_key": "virus",
                "species": "Chikungunya virus",
            },
            "hpiv": {
                "task_name": "demo_hpiv",
                "input_path": workspace_root_path / "demo_data" / "hpiv_demo",
                "output_path": workspace_root_path / "demo_data" / "hpiv_demo",
                "analysis_target": "virus",
                "inputtype": "fastq",
                "method": "freebayes",
                "asm_type": "shortref",
                "runflow": "物种鉴定,基因组组装,分型鉴定",
                "pipeline_script": "demo_data/hpiv_demo",
                "workstation_key": "virus",
                "species": "Human parainfluenza virus",
            },
            "hadv": {
                "task_name": "demo_hadv",
                "input_path": workspace_root_path / "demo_data" / "hadv_demo",
                "output_path": workspace_root_path / "demo_data" / "hadv_demo",
                "analysis_target": "virus",
                "inputtype": "fastq",
                "method": "freebayes",
                "asm_type": "shortref",
                "runflow": "物种鉴定,基因组组装,分型鉴定",
                "pipeline_script": "demo_data/hadv_demo",
                "workstation_key": "virus",
                "species": "Human adenovirus",
            },
            "norovirus": {
                "task_name": "demo_norovirus",
                "input_path": workspace_root_path / "demo_data" / "norovirus_demo",
                "output_path": workspace_root_path / "demo_data" / "norovirus_demo",
                "analysis_target": "virus",
                "inputtype": "fastq",
                "method": "freebayes",
                "asm_type": "shortref",
                "runflow": "物种鉴定,基因组组装,分型鉴定",
                "pipeline_script": "demo_data/norovirus_demo",
                "workstation_key": "virus",
                "species": "Norovirus",
            },
            "enterovirus": {
                "task_name": "demo_enterovirus",
                "input_path": workspace_root_path / "demo_data" / "enterovirus",
                "output_path": workspace_root_path / "demo_data" / "enterovirus",
                "analysis_target": "virus",
                "inputtype": "fastq",
                "method": "freebayes",
                "asm_type": "shortref",
                "runflow": "物种鉴定,基因组组装,分型鉴定",
                "pipeline_script": "demo_data/enterovirus",
                "workstation_key": "virus",
                "species": "Human enterovirus",
            },
            "hepatovirus": {
                "task_name": "demo_hepatovirus",
                "input_path": workspace_root_path / "demo_data" / "HAV_demo",
                "output_path": workspace_root_path / "demo_data" / "HAV_demo",
                "analysis_target": "virus",
                "inputtype": "fastq",
                "method": "freebayes",
                "asm_type": "shortref",
                "runflow": "物种鉴定,基因组组装,分型鉴定",
                "pipeline_script": "demo_data/HAV_demo",
                "workstation_key": "virus",
                "species": "Hepatovirus",
            },
            "hiv": {
                "task_name": "demo_hiv",
                "input_path": workspace_root_path / "demo_data" / "hiv_demo",
                "output_path": workspace_root_path / "demo_data" / "hiv_demo",
                "analysis_target": "virus",
                "inputtype": "fastq",
                "method": "freebayes",
                "asm_type": "shortref",
                "runflow": "物种鉴定,基因组组装,分型鉴定",
                "pipeline_script": "demo_data/hiv_demo",
                "workstation_key": "virus",
                "species": "HIV-1",
            },
            "bandavirus": {
                "task_name": "demo_bandavirus",
                "input_path": workspace_root_path / "demo_data" / "sftsv_demo" / "sftsv_samples.tsv",
                "output_path": workspace_root_path / "demo_data" / "sftsv_demo",
                "analysis_target": "virus",
                "inputtype": "fastq",
                "method": "freebayes",
                "asm_type": "shortref",
                "runflow": "物种鉴定,基因组组装,分型鉴定",
                "pipeline_script": "demo_data/sftsv_demo",
                "workstation_key": "virus",
                "species": "Bandavirus dabieense",
            },
            "orthohantavirus": {
                "task_name": "demo_orthohantavirus",
                "input_path": workspace_root_path / "demo_data" / "orthohantavirus_demo" / "orthohantavirus_samples.tsv",
                "output_path": workspace_root_path / "demo_data" / "orthohantavirus_demo",
                "analysis_target": "virus",
                "inputtype": "fastq",
                "method": "freebayes",
                "asm_type": "shortref",
                "runflow": "物种鉴定,基因组组装,分型鉴定",
                "pipeline_script": "demo_data/orthohantavirus_demo",
                "workstation_key": "virus",
                "species": "Orthohantavirus",
            },
            "orthoebolavirus": {
                "task_name": "demo_orthoebolavirus",
                "input_path": workspace_root_path / "demo_data" / "orthoebolavirus_demo" / "orthoebolavirus_samples.tsv",
                "output_path": workspace_root_path / "demo_data" / "orthoebolavirus_demo",
                "analysis_target": "virus",
                "inputtype": "fastq",
                "method": "freebayes",
                "asm_type": "shortref",
                "runflow": "物种鉴定,基因组组装,分型鉴定",
                "pipeline_script": "demo_data/orthoebolavirus_demo",
                "workstation_key": "virus",
                "species": "Ebola virus",
            },
            "astroviridae": {
                "task_name": "demo_astroviridae",
                "input_path": workspace_root_path / "demo_data" / "astro_demo_paired",
                "output_path": workspace_root_path / "demo_data" / "astro_demo_paired",
                "analysis_target": "virus",
                "inputtype": "fastq",
                "method": "freebayes",
                "asm_type": "shortref",
                "runflow": "物种鉴定,基因组组装,分型鉴定",
                "pipeline_script": "demo_data/astro_demo_paired",
                "workstation_key": "virus",
                "species": "Human astrovirus",
            },
            "rhinovirus": {
                "task_name": "demo_rhinovirus",
                "input_path": workspace_root_path / "demo_data" / "rhinovirus",
                "output_path": workspace_root_path / "demo_data" / "rhinovirus_demo",
                "analysis_target": "virus",
                "inputtype": "fastq",
                "method": "freebayes",
                "asm_type": "shortref",
                "runflow": "物种鉴定,基因组组装,分型鉴定",
                "pipeline_script": "demo_data/rhinovirus",
                "workstation_key": "virus",
                "species": "Human rhinovirus",
            },
            "seasonal_hcov": {
                "task_name": "demo_seasonal_hcov",
                "input_path": workspace_root_path / "demo_data" / "seasonal_hcov",
                "output_path": workspace_root_path / "demo_data" / "seasonal_hcov",
                "analysis_target": "virus",
                "inputtype": "fastq",
                "method": "freebayes",
                "asm_type": "shortref",
                "runflow": "物种鉴定,基因组组装,分型鉴定",
                "pipeline_script": "demo_data/seasonal_hcov",
                "workstation_key": "virus",
                "species": "Human coronavirus",
            },
            "rotavirus": {
                "task_name": "demo_rotavirus",
                "input_path": workspace_root_path / "demo_data" / "rotavirus",
                "output_path": workspace_root_path / "demo_data" / "rotavirus",
                "analysis_target": "virus",
                "inputtype": "fastq",
                "method": "freebayes",
                "asm_type": "shortref",
                "runflow": "物种鉴定,基因组组装,分型鉴定",
                "pipeline_script": "demo_data/rotavirus",
                "workstation_key": "virus",
                "species": "Rotavirus",
            },
            "hmpxv": {
                "task_name": "demo_hmpxv",
                "input_path": workspace_root_path / "demo_data" / "hmpxv_demo" / "hmpxv_demo.final.fasta",
                "output_path": workspace_root_path / "demo_data" / "hmpxv_demo",
                "analysis_target": "virus",
                "inputtype": "fasta",
                "method": "spades",
                "asm_type": "shortref",
                "runflow": "分型鉴定",
                "pipeline_script": "demo_data/hmpxv_demo",
                "workstation_key": "virus",
                "species": "Monkeypox virus",
            },
            "tree": {
                "task_name": "demo_tree",
                "input_path": workspace_root_path / "demo_data" / "Tree",
                "output_path": workspace_root_path / "demo_data" / "Tree",
                "analysis_target": "bacteria",
                "inputtype": "fasta",
                "method": "snippy",
                "asm_type": "shortasm",
                "runflow": "物种鉴定",
                "pipeline_script": "PathoSource.py",
                "workstation_key": "pathosource",
                "pathosource_species": "Bordetella pertussis",
            },
            "community": {
                "task_name": "demo_community",
                "input_path": str(workspace_root_path / "demo_data" / "community_demo"),
                "metadata_path": str(workspace_root_path / "demo_data" / "community_demo" / "sample-metadata.tsv"),
                "taxonomy_path": str(workspace_root_path / "demo_data" / "community_demo" / "community_taxonomy.tsv"),
                "output_path": str(workspace_root_path / "demo_data" / "community_demo" / "output"),
                "analysis_target": "bacteria",
                "inputtype": "directory",
                "method": "community",
                "asm_type": "shortasm",
                "runflow": "群落分析",
                "pipeline_script": "CommunityAnalysis.py",
                "workstation_key": "community",
            },
        }
        demo_config = demo_definitions.get(demo_type)
        if demo_config is None:
            raise ValidationError(f"不支持的 Demo 类型: {demo_type}")

        demo_input = Path(demo_config["input_path"]).resolve()
        if not demo_input.exists():
            raise ValidationError(f"Demo 数据不存在: {demo_input}")

        task_id = datetime.now().strftime("%Y%m%d%H%M%S") + "_" + uuid4().hex[:8]
        task_dir = self.task_root / task_id
        task_dir.mkdir(parents=True, exist_ok=False)
        task_file = task_dir / "task.json"
        log_file = task_dir / "pipeline.log"
        if "output_path" in demo_config:
            output_dir = Path(demo_config["output_path"]).resolve()
        else:
            output_dir = task_dir / "demo_output"
            output_dir.mkdir(parents=True, exist_ok=True)

        params = {
            "task_name": str(demo_config["task_name"]),
            "input_path": str(demo_input),
            "analysis_target": str(demo_config["analysis_target"]),
            "inputtype": str(demo_config["inputtype"]),
            "output_dir": str(output_dir),
            "thread": 10,
            "minlongfilt": "500",
            "Qfilt": "10",
            "barcodekit": "none",
            "method": str(demo_config["method"]),
            "long_type": "Nanopore",
            "ref": "noref",
            "gtf": "nogtf",
            "genome_len": "4m",
            "asm_type": str(demo_config["asm_type"]),
            "polish_times": "1",
            "polish_soft": "medaka",
            "species": "False",
            "rmhost": "norm",
            "runflow": str(demo_config["runflow"]),
            "abun": "1",
            "rna": "0",
            "fake_pip": 1,
            "watch_mode": "0",
            "watch_stable_minutes": 30,
            "watch_poll_minutes": 5,
            "watch_max_samples": 0,
        }
        if demo_type == "ncov":
            params.update({
                "species": str(demo_config["species"]),
                "ref": str(Path(demo_config["ref"]).resolve()),
                "gtf": "nogtf",
                "fake_pip": 0,
                "workstation_key": "virus",
                "sample_name": "test",
            })
        elif demo_type == "tb":
            params.update({
                "species": "Mycobacterium tuberculosis",
                "fake_pip": 0,
                "sample_name": "tb_demo",
            })
        elif demo_type == "flu":
            params.update({
                "species": str(demo_config["species"]),
                "ref": "noref",
                "gtf": "nogtf",
                "fake_pip": 0,
                "workstation_key": "virus",
                "sample_name": "flu_demo",
            })
        elif demo_type == "rsv":
            params.update({
                "species": str(demo_config["species"]),
                "ref": "noref",
                "gtf": "nogtf",
                "fake_pip": 0,
                "workstation_key": "virus",
                "sample_name": "rsv_demo",
            })
        elif demo_type == "hmpv":
            hmpv_db_root = workspace_root_path / "database" / "nextclade_db" / "hmpv"
            if not hmpv_db_root.is_dir():
                hmpv_db_root = workspace_root_path / "database" / "virus" / "nextclade" / "hmpv"
            params.update({
                "species": str(demo_config["species"]),
                "ref": str((hmpv_db_root / "reference.fasta").resolve()) if (hmpv_db_root / "reference.fasta").is_file() else "noref",
                "gtf": str((hmpv_db_root / "genome_annotation.gff3").resolve()) if (hmpv_db_root / "genome_annotation.gff3").is_file() else "nogtf",
                "fake_pip": 0,
                "workstation_key": "virus",
                "sample_name": "hmpv_demo",
            })
        elif demo_type == "denv":
            params.update({
                "species": str(demo_config["species"]),
                "ref": "noref",
                "gtf": "nogtf",
                "fake_pip": 0,
                "workstation_key": "virus",
                "sample_name": "denv_demo",
            })
        elif demo_type == "zikav":
            zikav_db_root = workspace_root_path / "database" / "nextclade_db" / "zikav"
            params.update({
                "species": str(demo_config["species"]),
                "ref": str((zikav_db_root / "reference.fasta").resolve()) if (zikav_db_root / "reference.fasta").is_file() else "noref",
                "gtf": str((zikav_db_root / "genome_annotation.gff3").resolve()) if (zikav_db_root / "genome_annotation.gff3").is_file() else "nogtf",
                "fake_pip": 0,
                "workstation_key": "virus",
                "sample_name": "zikav_demo_1",
            })
        elif demo_type == "chikv":
            chikv_db_root = workspace_root_path / "database" / "nextclade_db" / "chikv"
            params.update({
                "species": str(demo_config["species"]),
                "ref": str((chikv_db_root / "reference.fasta").resolve()) if (chikv_db_root / "reference.fasta").is_file() else "noref",
                "gtf": str((chikv_db_root / "genome_annotation.gff3").resolve()) if (chikv_db_root / "genome_annotation.gff3").is_file() else "nogtf",
                "fake_pip": 0,
                "workstation_key": "virus",
                "sample_name": "chikv_demo_1",
            })
        elif demo_type == "hpiv":
            params.update({
                "species": str(demo_config["species"]),
                "ref": "noref",
                "gtf": "nogtf",
                "fake_pip": 0,
                "workstation_key": "virus",
                "sample_name": "hpiv_demo",
            })
        elif demo_type == "hadv":
            params.update({
                "species": str(demo_config["species"]),
                "ref": "noref",
                "gtf": "nogtf",
                "fake_pip": 0,
                "workstation_key": "virus",
                "sample_name": "hadv_demo",
            })
        elif demo_type == "norovirus":
            params.update({
                "species": str(demo_config["species"]),
                "ref": "noref",
                "gtf": "nogtf",
                "fake_pip": 0,
                "workstation_key": "virus",
                "sample_name": "norovirus_demo",
            })
        elif demo_type == "enterovirus":
            params.update({
                "species": str(demo_config["species"]),
                "ref": "noref",
                "gtf": "nogtf",
                "fake_pip": 0,
                "workstation_key": "virus",
                "sample_name": "ev_a_demo",
            })
        elif demo_type == "hepatovirus":
            params.update({
                "species": str(demo_config["species"]),
                "ref": "noref",
                "gtf": "nogtf",
                "fake_pip": 0,
                "workstation_key": "virus",
                "sample_name": "HAV_demo",
            })
        elif demo_type == "hiv":
            params.update({
                "species": str(demo_config["species"]),
                "ref": "noref",
                "gtf": "nogtf",
                "fake_pip": 0,
                "workstation_key": "virus",
                "sample_name": "hiv_demo",
            })
        elif demo_type == "bandavirus":
            params.update({
                "species": str(demo_config["species"]),
                "ref": "noref",
                "gtf": "nogtf",
                "fake_pip": 0,
                "workstation_key": "virus",
                "sample_name": "sftsv_demo",
            })
        elif demo_type == "orthohantavirus":
            params.update({
                "species": str(demo_config["species"]),
                "ref": "noref",
                "gtf": "nogtf",
                "fake_pip": 0,
                "workstation_key": "virus",
                "sample_name": "orthohantavirus_demo",
            })
        elif demo_type == "orthoebolavirus":
            params.update({
                "species": str(demo_config["species"]),
                "ref": "noref",
                "gtf": "nogtf",
                "fake_pip": 0,
                "workstation_key": "virus",
                "sample_name": "orthoebolavirus_demo",
                "genome_len": "19k",
                "rna": "1",
            })
        elif demo_type == "astroviridae":
            params.update({
                "species": str(demo_config["species"]),
                "ref": "noref",
                "gtf": "nogtf",
                "fake_pip": 0,
                "workstation_key": "virus",
                "sample_name": "astro_demo_paired",
            })
        elif demo_type == "rhinovirus":
            params.update({
                "species": str(demo_config["species"]),
                "ref": "noref",
                "gtf": "nogtf",
                "fake_pip": 0,
                "workstation_key": "virus",
                "sample_name": "rhinovirus_demo",
            })
        elif demo_type == "seasonal_hcov":
            params.update({
                "species": str(demo_config["species"]),
                "ref": "noref",
                "gtf": "nogtf",
                "fake_pip": 0,
                "workstation_key": "virus",
                "sample_name": "hcov_229e_demo",
            })
        elif demo_type == "rotavirus":
            params.update({
                "species": str(demo_config["species"]),
                "ref": "noref",
                "gtf": "nogtf",
                "fake_pip": 0,
                "workstation_key": "virus",
                "sample_name": "rotavirus_a_demo",
            })
        elif demo_type == "hmpxv":
            hmpxv_db_root = workspace_root_path / "database" / "virus" / "nextclade" / "hMPXV"
            if not hmpxv_db_root.is_dir():
                hmpxv_db_root = workspace_root_path / "database" / "nextclade_db" / "hMPXV"
            params.update({
                "species": str(demo_config["species"]),
                "ref": str((hmpxv_db_root / "reference.fasta").resolve()) if (hmpxv_db_root / "reference.fasta").is_file() else "noref",
                "gtf": str((hmpxv_db_root / "genome_annotation.gff3").resolve()) if (hmpxv_db_root / "genome_annotation.gff3").is_file() else "nogtf",
                "fake_pip": 0,
                "workstation_key": "virus",
                "sample_name": "hmpxv_demo",
            })
        elif demo_type == "tree":
            params.update({
                "species": "Bordetella pertussis",
                "thread": 8,
                "ref": "False",
                "gtf": "nogtf",
                "fake_pip": 0,
                "workstation_key": "pathosource",
                "cgmlstana": "yes",
                "gubbins": "yes",
                "msamethod": "snippy",
                "treemethod": "NJ",
                "Bootstrap": 500,
                "mltype": "MFP",
                "mode": "P",
                "cgmlst": "none",
            })
        elif demo_type == "community":
            params = {
                "task_name": str(demo_config["task_name"]),
                "input_path": str(demo_input),
                "analysis_target": "bacteria",
                "inputtype": "directory",
                "output_dir": str(output_dir),
                "thread": 8,
                "method": "community",
                "metadata": str(Path(demo_config["metadata_path"]).resolve()),
                "taxonomy": str(Path(demo_config["taxonomy_path"]).resolve()),
                "group_column": "SampleType",
                "taxonomy_level": "genus",
                "normalization": "relative",
                "analyses": "alpha,beta,lefse,ml,network",
                "alpha_metrics": "shannon,simpson,chao1",
                "beta_metric": "braycurtis",
                "ml_model": "random_forest",
                "watch_mode": "0",
                "watch_stable_minutes": 30,
                "watch_poll_minutes": 5,
                "watch_max_samples": 0,
                "community_mode": "amplicon",
                "workstation_key": "community",
            }
        finished_at = utc_now_iso()
        task = {
            "id": task_id,
            "name": str(demo_config["task_name"]),
            "is_demo": True,
            "demo_type": demo_type,
            "owner": owner,
            "owner_group": owner_group,
            "status": "SUCCEEDED",
            "created_at": finished_at,
            "started_at": finished_at,
            "finished_at": finished_at,
            "exit_code": 0,
            "runner_pid": None,
            "pipeline_pid": None,
            "resume_pending": False,
            "requested_status": "",
            "params": params,
            "command": ["demo_task", str(demo_input)],
            "project_root": str(workspace_root_path),
            "log_path": str(log_file),
            "pipeline_script": str(demo_config["pipeline_script"]),
        }
        log_file.write_text(
            f"[{finished_at}] Demo 任务已生成\n"
            f"[{finished_at}] 输入目录: {demo_input}\n"
            f"[{finished_at}] 状态: SUCCEEDED\n",
            encoding="utf-8",
        )
        write_json(task_file, task)
        return self.get_task(task_id, owner=owner)

    def _is_pathosource_payload(self, payload: dict[str, Any]) -> bool:
        workstation_key = self._clean_str(payload.get("workstation_key")).lower()
        pipeline_script = Path(str(payload.get("pipeline_script") or "")).name
        return workstation_key == "pathosource" or pipeline_script == "PathoSource.py"

    def _is_community_payload(self, payload: dict[str, Any]) -> bool:
        workstation_key = self._clean_str(payload.get("workstation_key")).lower()
        pipeline_script = Path(str(payload.get("pipeline_script") or "")).name
        return workstation_key == "community" or pipeline_script == "CommunityAnalysis.py"

    def _normalize_pathosource_payload(self, payload: dict[str, Any]) -> dict[str, Any]:
        input_path = self._resolve_existing_path(payload.get("input_path"), "input_path")
        task_name = self._clean_str(payload.get("task_name")) or f"pathosource_{Path(input_path).stem}"
        output_dir = self._resolve_task_output_path(payload.get("output_dir"), task_name)
        output_dir.mkdir(parents=True, exist_ok=True)
        meta = self._clean_str(payload.get("pathosource_meta"))
        ref = self._clean_str(payload.get("ref"))
        msamethod = self._clean_str(payload.get("pathosource_msamethod")) or "snippy"
        msa_alias = {"ska": "ska2", "ska2": "ska2", "roary": "roray", "roray": "roray"}
        msamethod = msa_alias.get(msamethod, msamethod)
        return {
            "task_name": task_name,
            "input_path": input_path,
            "output_dir": str(output_dir),
            "thread": self._clean_int(payload.get("thread"), default=10, minimum=1),
            "species": self._clean_str(payload.get("species")) or "salmonella",
            "ref": self._resolve_optional_path(ref, "False") if ref else "False",
            "meta": self._resolve_optional_path(meta, "") if meta else "",
            "cgmlstana": self._clean_str(payload.get("pathosource_cgmlstana")) or "no",
            "gubbins": self._clean_str(payload.get("pathosource_gubbins")) or "yes",
            "msamethod": msamethod,
            "treemethod": self._clean_str(payload.get("pathosource_treemethod")) or "ML",
            "Bootstrap": self._clean_int(payload.get("pathosource_bootstrap"), default=1000, minimum=0),
            "mltype": self._clean_str(payload.get("pathosource_mltype")) or "MFP",
            "mode": self._clean_str(payload.get("pathosource_mode")) or "P",
            "cgmlst": self._clean_str(payload.get("pathosource_cgmlstversion")) or "none",
            "watch_mode": "0",
            "watch_stable_minutes": self._clean_int(payload.get("watch_stable_minutes"), default=30, minimum=1),
            "watch_poll_minutes": self._clean_int(payload.get("watch_poll_minutes"), default=5, minimum=1),
            "watch_max_samples": self._clean_int(payload.get("watch_max_samples"), default=0, minimum=0),
            "workstation_key": "pathosource",
        }

    def _normalize_community_payload(self, payload: dict[str, Any]) -> dict[str, Any]:
        input_path = self._resolve_existing_path(payload.get("input_path"), "input_path")
        task_name = self._clean_str(payload.get("task_name")) or f"community_{Path(input_path).stem}"
        output_dir = self._resolve_task_output_path(payload.get("output_dir"), task_name)
        output_dir.mkdir(parents=True, exist_ok=True)
        metadata_path = self._resolve_existing_path(payload.get("community_metadata"), "community_metadata")
        input_path_obj = Path(input_path)
        amplicon_mode = (
            (input_path_obj / "demux.qzv").is_file()
            or (input_path_obj / "demux.qza").is_file()
            or input_path_obj.name.endswith(".qzv")
            or input_path_obj.name.endswith(".qza")
        )
        taxonomy_raw = self._clean_str(payload.get("community_taxonomy"))
        taxonomy_path = self._resolve_existing_path(taxonomy_raw, "community_taxonomy") if taxonomy_raw else ""
        group_column = self._clean_str(payload.get("community_group_column"))
        if not group_column:
            raise ValidationError("多样本群落分析必须提供分组列")
        if not amplicon_mode and not taxonomy_path:
            raise ValidationError("未检测到 demux.qzv/qza 时，群落分析必须提供 taxonomy 文件")
        analyses = [item.strip().lower() for item in str(payload.get("community_analyses") or "").split(",") if item.strip()]
        allowed_analyses = ["alpha", "beta", "lefse", "ml", "network"]
        normalized_analyses = [item for item in allowed_analyses if item in analyses]
        if not normalized_analyses:
            normalized_analyses = ["alpha", "beta"]
        merge_items = payload.get("community_merge_tasks") if isinstance(payload.get("community_merge_tasks"), list) else []
        return {
            "task_name": task_name,
            "input_path": input_path,
            "analysis_target": "bacteria",
            "inputtype": "directory",
            "output_dir": str(output_dir),
            "thread": self._clean_int(payload.get("thread"), default=8, minimum=1),
            "method": "community",
            "metadata": str(metadata_path),
            "taxonomy": str(taxonomy_path) if taxonomy_path else "",
            "group_column": group_column,
            "taxonomy_level": self._clean_str(payload.get("community_taxonomy_level")) or "genus",
            "normalization": self._clean_str(payload.get("community_normalization")) or "relative",
            "analyses": ",".join(normalized_analyses),
            "alpha_metrics": self._clean_str(payload.get("community_alpha_metrics")) or "shannon,simpson,chao1",
            "beta_metric": self._clean_str(payload.get("community_beta_metric")) or "braycurtis",
            "ml_model": self._clean_str(payload.get("community_ml_model")) or "random_forest",
            "watch_mode": "0",
            "watch_stable_minutes": self._clean_int(payload.get("watch_stable_minutes"), default=30, minimum=1),
            "watch_poll_minutes": self._clean_int(payload.get("watch_poll_minutes"), default=5, minimum=1),
            "watch_max_samples": self._clean_int(payload.get("watch_max_samples"), default=0, minimum=0),
            "community_mode": "amplicon" if amplicon_mode else "abundance",
            "community_input_source": self._clean_str(payload.get("community_input_source")) or "manual",
            "community_merge_tasks": merge_items,
            "community_generated_abundance": self._clean_str(payload.get("community_generated_abundance")),
            "community_generated_manifest": self._clean_str(payload.get("community_generated_manifest")),
            "workstation_key": "community",
        }

    def _normalize_payload(self, payload: dict[str, Any]) -> dict[str, Any]:
        if self._is_pathosource_payload(payload):
            return self._normalize_pathosource_payload(payload)
        if self._is_community_payload(payload):
            return self._normalize_community_payload(payload)
        workstation_key = self._clean_str(payload.get("workstation_key")).lower() or "bacteria"
        if workstation_key in {"pathogen", "single"}:
            workstation_key = "bacteria"
        input_path = self._resolve_existing_path(payload.get("input_path"), "input_path")
        task_name = self._clean_str(payload.get("task_name")) or f"bac_assemble_{Path(input_path).stem}"
        output_dir = self._resolve_task_output_path(payload.get("output_dir"), task_name)
        output_dir.mkdir(parents=True, exist_ok=True)
        if workstation_key == "virus":
            analysis_target = "virus"
        elif workstation_key == "metagenome":
            analysis_target = "bacteria"
        else:
            analysis_target = "bacteria"
        asm_type_raw = payload.get("asm_type")
        if analysis_target == "virus" and asm_type_raw in (None, ""):
            asm_type_raw = "shortref"
        if workstation_key == "metagenome":
            asm_type_raw = "shortasm"
        asm_type = self._normalize_asm_type(asm_type_raw)
        if workstation_key == "virus" and asm_type not in {"shortref", "longref"}:
            raise ValidationError("病毒分析仅支持短读长有参组装(shortref)或长读长有参组装(longref)")
        method = self._normalize_method(payload.get("method"), asm_type)
        if workstation_key == "metagenome":
            method = "meta"
        watch_mode = "1" if str(payload.get("watch_mode") or "0").strip() == "1" else "0"
        inputtype = self._clean_str(payload.get("inputtype")) or "fastq"
        if watch_mode == "1":
            input_candidate = Path(input_path)
            if not input_candidate.is_dir():
                raise ValidationError("监听任务要求输入路径必须是目录")
            if inputtype.lower() != "fastq":
                raise ValidationError("监听任务当前仅支持 fastq 输入目录")

        return {
            "task_name": task_name,
            "input_path": input_path,
            "analysis_target": analysis_target,
            "inputtype": inputtype,
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
            "species": "False" if workstation_key == "metagenome" else (self._clean_str(payload.get("species")) or "False"),
            "rmhost": self._normalize_rmhost(payload.get("rmhost")),
            "runflow": self._clean_str(payload.get("runflow")) or default_runflow_for_selection(method, analysis_target),
            "abun": self._clean_str(payload.get("abun")) or "1",
            "rna": self._clean_str(payload.get("rna")) or "0",
            "fake_pip": self._clean_int(payload.get("fake_pip"), default=0, minimum=0),
            "watch_mode": watch_mode,
            "watch_stable_minutes": self._clean_int(payload.get("watch_stable_minutes"), default=30, minimum=1),
            "watch_poll_minutes": self._clean_int(payload.get("watch_poll_minutes"), default=5, minimum=1),
            "watch_max_samples": self._clean_int(payload.get("watch_max_samples"), default=0, minimum=0),
            "workstation_key": workstation_key,
        }

    def _build_command(
        self,
        params: dict[str, Any],
        pipeline_script: Path,
        pipeline_python: str,
        conda_exe: str = "conda",
    ) -> list[str]:
        if pipeline_script.name == "PathoSource.py" or str(params.get("workstation_key") or "").lower() == "pathosource":
            command = [
                conda_exe,
                "run",
                "-n",
                pipeline_python,
                "--no-capture-output",
                "python",
                "-u",
                str(pipeline_script),
                "--input",
                params["input_path"],
                "--species",
                params["species"],
                "--threads",
                str(params["thread"]),
                "--output",
                params["output_dir"],
                "--ref",
                params["ref"],
                "--cgmlstana",
                params["cgmlstana"],
                "--gubbins",
                params.get("gubbins", "yes"),
                "--msamethod",
                params["msamethod"],
                "--treemethod",
                params["treemethod"],
                "--Bootstrap",
                str(params["Bootstrap"]),
                "--mltype",
                params["mltype"],
                "--mode",
                params["mode"],
                "--cgmlst",
                params["cgmlst"],
            ]
            if str(params.get("meta") or "").strip():
                command.extend(["--meta", params["meta"]])
            return command
        if pipeline_script.name == "CommunityAnalysis.py" or str(params.get("workstation_key") or "").lower() == "community":
            command = [
                conda_exe,
                "run",
                "-n",
                pipeline_python,
                "--no-capture-output",
                "python",
                "-u",
                str(pipeline_script),
                "--input",
                params["input_path"],
                "--metadata",
                params["metadata"],
                "--output",
                params["output_dir"],
                "--thread",
                str(params["thread"]),
                "--group-column",
                params["group_column"],
                "--taxonomy-level",
                params["taxonomy_level"],
                "--normalization",
                params["normalization"],
                "--analyses",
                params["analyses"],
                "--alpha-metrics",
                params["alpha_metrics"],
                "--beta-metric",
                params["beta_metric"],
                "--ml-model",
                params["ml_model"],
            ]
            if str(params.get("taxonomy") or "").strip():
                command[command.index("--output"):command.index("--output")] = ["--taxonomy", params["taxonomy"]]
            return command
        return [
            conda_exe,
            "run",
            "-n",
            pipeline_python,
            "--no-capture-output",
            "python",
            "-u",
            str(pipeline_script),
            "--input",
            params["input_path"],
            "--analysis_target",
            params["analysis_target"],
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
        log_path = Path(task.get("log_path", ""))
        task_dir = self.task_root / str(task.get("id") or "")
        progress = self._build_progress(task, log_path)
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
            "resume_pending": bool(task.get("resume_pending")),
            "requested_status": task.get("requested_status", ""),
            "params": task.get("params", {}),
            "command": task.get("command", []),
            "log_path": task.get("log_path"),
            "pipeline_script": task.get("pipeline_script", ""),
            "conda_root": task.get("conda_root", ""),
            "is_demo": bool(task.get("is_demo")),
            "demo_type": str(task.get("demo_type") or ""),
            "progress": progress,
            "last_log_line": self._read_last_log_line(log_path),
            "watch": task.get("watch", {}),
            "review_gate": task.get("review_gate", {}),
            "parent_task_id": task.get("parent_task_id", ""),
            "trigger_context": task.get("trigger_context", {}),
            "auto_pathosource": task.get("auto_pathosource", {}),
            "analytics_snapshot": read_task_analytics_snapshot(task_dir),
        }
        if include_log:
            payload["log_tail"] = self._read_log_tail(log_path, log_lines)
        return payload

    def update_task_fields(self, task_id: str, updates: dict[str, Any], *, owner: str | None = None) -> dict[str, Any]:
        task_file, task = self._load_task_for_change(task_id, owner=owner)
        for key, value in dict(updates or {}).items():
            if key in {"id", "log_path"}:
                continue
            task[key] = value
        write_json(task_file, task)
        return self._serialize_task(task, include_log=False)

    def _read_log_tail(self, path: Path, log_lines: int) -> str:
        if not path.is_file():
            return ""
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        return "\n".join(lines[-log_lines:])

    def _read_last_log_line(self, path: Path) -> str:
        if not path.is_file():
            return ""
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        for line in reversed(lines):
            text = line.strip()
            if text:
                return text
        return ""

    def _build_progress(self, task: dict[str, Any], log_path: Path) -> dict[str, Any]:
        status = str(task.get("status") or "").upper()
        review_gate = task.get("review_gate") or {}
        progress = {
            "overall_percent": 0.0,
            "sample_percent": 0.0,
            "sample_index": None,
            "sample_total": None,
            "sample_name": "",
            "step": None,
            "total_step": None,
            "message": "",
            "label": "等待调度",
        }
        watch = task.get("watch") or {}
        if status == "QUEUED" and self._monitor_waiting(task):
            progress["label"] = str(watch.get("label") or "监听中 · 等待测序结束")
            progress["message"] = str(watch.get("message") or "")
            progress["sample_total"] = int(watch.get("sample_count") or 0) or None
            progress["overall_percent"] = float(watch.get("progress_percent") or 0.0)
            progress["sample_percent"] = float(watch.get("progress_percent") or 0.0)
            return progress
        if status == "QUEUED" and not task.get("resume_pending"):
            progress["label"] = "等待任务启动"
            return progress
        if status in {"SUCCEEDED", "FAILED"} and not log_path.is_file():
            progress["overall_percent"] = 100.0 if status == "SUCCEEDED" else 0.0
            progress["sample_percent"] = 100.0 if status == "SUCCEEDED" else 0.0
            progress["label"] = "任务已完成" if status == "SUCCEEDED" else "任务已失败"
            return progress
        if not log_path.is_file():
            progress["label"] = "日志尚未生成"
            return progress

        latest_match: re.Match[str] | None = None
        for line in log_path.read_text(encoding="utf-8", errors="replace").splitlines():
            matched = PROGRESS_LINE_RE.search(line)
            if matched:
                latest_match = matched

        if latest_match is None:
            progress["overall_percent"] = 100.0 if status == "SUCCEEDED" else 0.0
            progress["sample_percent"] = 100.0 if status == "SUCCEEDED" else 0.0
            progress["label"] = "任务已完成" if status == "SUCCEEDED" else ("任务已暂停" if status == "PAUSED" else ("运行中，等待阶段日志" if status == "RUNNING" else ("任务已停止" if status == "STOPPED" else "任务已失败")))
            return progress

        step = int(latest_match.group("step"))
        total_step = max(int(latest_match.group("total_step")), 1)
        sample_index = max(int(latest_match.group("sample_index")), 1)
        sample_total = max(int(latest_match.group("sample_total")), 1)
        sample_name = latest_match.group("sample").strip()
        message = latest_match.group("message").strip()
        sample_percent = max(0.0, min(100.0, (step / total_step) * 100.0))
        overall_fraction = ((sample_index - 1) + (step / total_step)) / sample_total
        overall_percent = max(0.0, min(100.0, overall_fraction * 100.0))

        if status == "SUCCEEDED":
            sample_percent = 100.0
            overall_percent = 100.0

        label = f"样本 {sample_index}/{sample_total} · 步骤 {step}/{total_step}"
        if status == "PAUSED":
            if str(review_gate.get("state") or "").strip().lower() == "pending":
                label = f"待人工确认 · 样本 {sample_index}/{sample_total} · 步骤 {step}/{total_step}"
                message = str(review_gate.get("summary") or message).strip()
            else:
                label = f"已暂停 · 样本 {sample_index}/{sample_total} · 步骤 {step}/{total_step}"
        elif status == "QUEUED" and task.get("resume_pending"):
            label = f"等待恢复 · 样本 {sample_index}/{sample_total} · 步骤 {step}/{total_step}"
        elif status == "STOPPED":
            label = f"已停止 · 样本 {sample_index}/{sample_total} · 步骤 {step}/{total_step}"

        progress.update(
            {
                "overall_percent": round(overall_percent, 1),
                "sample_percent": round(sample_percent, 1),
                "sample_index": sample_index,
                "sample_total": sample_total,
                "sample_name": sample_name,
                "step": step,
                "total_step": total_step,
                "message": message,
                "label": label,
            }
        )
        return progress

    def refresh_monitored_tasks(self) -> None:
        for task_file in sorted(self.task_root.glob("*/task.json")):
            try:
                task = read_json(task_file)
            except Exception:
                continue
            updated = False
            if str(task.get("status") or "").upper() == "RUNNING":
                result = evaluate_species_review_gate(task, Path(str(task.get("log_path") or "")))
                review_gate = task.get("review_gate") or {}
                if result and result.get("decision") == "pass":
                    evaluated = [str(item or "").strip() for item in (review_gate.get("evaluated_samples") or []) if str(item or "").strip()]
                    sample_name = str(result.get("sample_name") or "").strip()
                    if sample_name and sample_name not in evaluated:
                        evaluated.append(sample_name)
                        review_gate["evaluated_samples"] = evaluated
                        task["review_gate"] = review_gate
                        updated = True
            if not self._monitor_waiting(task):
                if updated:
                    write_json(task_file, task)
                continue
            if self._refresh_single_monitored_task(task, task_file) or updated:
                write_json(task_file, task)

    def _refresh_single_monitored_task(self, task: dict[str, Any], task_file: Path) -> bool:
        watch = task.setdefault("watch", {})
        params = task.get("params") or {}
        poll_minutes = max(1, int(watch.get("poll_minutes") or params.get("watch_poll_minutes") or 5))
        now = datetime.now(timezone.utc)
        last_scan = self._parse_task_time(watch.get("last_scan_at"))
        if last_scan and (now - last_scan).total_seconds() < poll_minutes * 60:
            return False

        input_dir = Path(str(watch.get("input_dir") or params.get("watch_input_path") or params.get("input_path") or "")).expanduser()
        if not input_dir.exists() or not input_dir.is_dir():
            watch.update({
                "state": "invalid_input",
                "label": "监听失败 · 输入目录不存在",
                "message": f"输入目录不存在：{input_dir}",
                "last_scan_at": now.isoformat(),
                "progress_percent": 0,
            })
            return True

        snapshot = self._build_monitor_snapshot(input_dir)
        sample_rows = self._detect_monitored_samples(input_dir, params)
        sample_names = [row[0] for row in sample_rows]
        watch["last_scan_at"] = now.isoformat()
        watch["sample_count"] = len(sample_rows)
        watch["sample_preview"] = sample_names[:12]
        watch["file_count"] = snapshot["file_count"]
        watch["total_size"] = snapshot["total_size"]

        previous_signature = str(watch.get("snapshot_signature") or "")
        current_signature = str(snapshot.get("signature") or "")
        stable_minutes = max(1, int(watch.get("stable_minutes") or params.get("watch_stable_minutes") or 30))
        max_samples = max(0, int(watch.get("max_samples") or params.get("watch_max_samples") or 0))

        if snapshot["file_count"] <= 0:
            watch.update({
                "state": "waiting_data",
                "label": "监听中 · 等待测序数据写入",
                "message": "当前目录尚未发现可识别的测序文件。",
                "progress_percent": 4,
            })
            return True

        if current_signature != previous_signature:
            watch.update({
                "snapshot_signature": current_signature,
                "last_change_at": now.isoformat(),
                "state": "watching",
                "label": "监听中 · 检测到目录仍在变化",
                "message": f"已识别 {len(sample_rows)} 个样本，目录内文件仍有新增或大小变化。",
                "progress_percent": 26,
            })
            return True

        last_change = self._parse_task_time(watch.get("last_change_at")) or now
        quiet_minutes = max(0, int((now - last_change).total_seconds() // 60))
        if quiet_minutes < stable_minutes:
            remaining = max(0, stable_minutes - quiet_minutes)
            watch.update({
                "state": "stabilizing",
                "label": "监听中 · 等待目录稳定",
                "message": f"已识别 {len(sample_rows)} 个样本，最近 {quiet_minutes} 分钟无变化，还需稳定 {remaining} 分钟。",
                "progress_percent": min(92, max(30, int((quiet_minutes / stable_minutes) * 100))),
            })
            return True

        if not sample_rows:
            watch.update({
                "state": "no_samples",
                "label": "监听完成 · 未识别到可运行样本",
                "message": "目录已稳定，但未识别到可组成样本的 fastq 文件。",
                "progress_percent": 100,
            })
            return True

        if max_samples and len(sample_rows) > max_samples:
            watch.update({
                "state": "limit_exceeded",
                "label": "监听完成 · 样本数超出自动运行上限",
                "message": f"已识别 {len(sample_rows)} 个样本，超过自动运行上限 {max_samples}，任务保持排队等待人工处理。",
                "progress_percent": 100,
            })
            return True

        batch_path = self._write_monitored_batch_input(task_file.parent, sample_rows)
        params["watch_input_path"] = str(input_dir)
        params["input_path"] = str(batch_path)
        task["params"] = params
        task["command"] = self._build_command(
            params,
            Path(str(task.get("pipeline_script") or "")),
            str(task.get("pipeline_python") or "base"),
            _resolve_conda_exe(str(task.get("conda_root") or "")),
        )
        watch.update({
            "state": "ready",
            "label": "监听完成 · 已生成批量输入",
            "message": f"已识别 {len(sample_rows)} 个样本，目录稳定，准备启动分析。",
            "generated_batch_path": str(batch_path),
            "materialized_at": now.isoformat(),
            "progress_percent": 100,
        })
        return True

    def _monitor_waiting(self, task: dict[str, Any]) -> bool:
        if str(task.get("status") or "").upper() != "QUEUED":
            return False
        watch = task.get("watch") or {}
        if not watch.get("enabled"):
            return False
        return str(watch.get("state") or "") not in {"ready", "dispatched"}

    def _build_initial_watch_state(self, params: dict[str, Any]) -> dict[str, Any]:
        return {
            "enabled": True,
            "state": "waiting_data",
            "label": "监听中 · 等待测序数据写入",
            "message": "任务已创建，将持续检查输入目录内测序文件是否还在增长。",
            "stable_minutes": int(params.get("watch_stable_minutes") or 30),
            "poll_minutes": int(params.get("watch_poll_minutes") or 5),
            "max_samples": int(params.get("watch_max_samples") or 0),
            "input_dir": str(params.get("input_path") or ""),
            "progress_percent": 0,
            "sample_count": 0,
            "sample_preview": [],
            "snapshot_signature": "",
            "last_scan_at": "",
            "last_change_at": "",
            "generated_batch_path": "",
        }

    def _parse_task_time(self, value: Any) -> datetime | None:
        text = str(value or "").strip()
        if not text:
            return None
        try:
            return datetime.fromisoformat(text)
        except ValueError:
            return None

    def _build_monitor_snapshot(self, input_dir: Path) -> dict[str, Any]:
        files: list[tuple[str, int]] = []
        total_size = 0
        for path in sorted(input_dir.rglob("*")):
            if not path.is_file():
                continue
            lower_name = path.name.lower()
            if not lower_name.endswith(MONITOR_INPUT_EXTENSIONS):
                continue
            try:
                size = path.stat().st_size
            except OSError:
                continue
            rel_path = str(path.relative_to(input_dir))
            files.append((rel_path, size))
            total_size += size
        signature = json.dumps(files, ensure_ascii=False, separators=(",", ":"))
        return {
            "file_count": len(files),
            "total_size": total_size,
            "signature": signature,
        }

    def _detect_monitored_samples(self, input_dir: Path, params: dict[str, Any]) -> list[list[str]]:
        groups: dict[str, dict[str, str]] = {}
        for path in sorted(input_dir.rglob("*")):
            if not path.is_file():
                continue
            lower_name = path.name.lower()
            if not lower_name.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
                continue
            sample_name, field = self._classify_monitored_fastq(path)
            record = groups.setdefault(sample_name, {"sample_name": sample_name, "third_gen": "", "short_left": "", "short_right": ""})
            record[field] = str(path.resolve())
        rows: list[list[str]] = []
        species = str(params.get("species") or "").strip()
        for sample_name in sorted(groups):
            record = groups[sample_name]
            if not any([record["third_gen"], record["short_left"], record["short_right"]]):
                continue
            rows.append([sample_name, record["third_gen"], record["short_left"], record["short_right"], species])
        return rows

    def _classify_monitored_fastq(self, path: Path) -> tuple[str, str]:
        name = path.name
        lower_name = name.lower()
        for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
            if lower_name.endswith(suffix):
                name = name[: -len(suffix)]
                break
        field = "third_gen"
        sample_name = name
        if SHORT_READ_R1_RE.search(name):
            field = "short_left"
            sample_name = SHORT_READ_R1_RE.sub("_", name)
        elif SHORT_READ_R2_RE.search(name):
            field = "short_right"
            sample_name = SHORT_READ_R2_RE.sub("_", name)
        sample_name = re.sub(r"[._-]+$", "", sample_name).strip() or path.stem
        return sample_name, field

    def _write_monitored_batch_input(self, task_dir: Path, rows: list[list[str]]) -> Path:
        target = task_dir / "monitored_batch_input.tsv"
        with target.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["样本名称", "三代数据", "二代数据左", "二代数据右", "物种信息"])
            writer.writerows(rows)
        return target

    def _load_task_for_change(self, task_id: str, *, owner: str | None = None) -> tuple[Path, dict[str, Any]]:
        task_file = self.task_root / task_id / "task.json"
        if not task_file.is_file():
            raise KeyError(f"Task not found: {task_id}")
        task = read_json(task_file)
        if owner and task.get("owner") != owner:
            raise KeyError(f"Task not found: {task_id}")
        return task_file, task

    def _start_runner(self, task: dict[str, Any], task_file: Path) -> None:
        runner = subprocess.Popen(
            [self.python_executable, "-m", "bac_analysis_portal.task_runner", str(task_file)],
            cwd=str(self.project_root),
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            start_new_session=True,
        )
        task["runner_pid"] = runner.pid
        task["status"] = "QUEUED"
        task["resume_pending"] = False
        write_json(task_file, task)

    def _signal_pipeline(self, task: dict[str, Any], sig: int) -> bool:
        pid = int(task.get("pipeline_pid") or 0)
        if pid <= 0:
            return False
        try:
            os.killpg(os.getpgid(pid), sig)
            return True
        except ProcessLookupError:
            return False
        except OSError:
            try:
                os.kill(pid, sig)
                return True
            except OSError:
                return False

    def _resume_pipeline(self, task: dict[str, Any]) -> bool:
        if not self._signal_pipeline(task, signal.SIGCONT):
            return False
        task["status"] = "RUNNING"
        task["resume_pending"] = False
        return True

    def _resolve_runtime_python(self, raw_value: str | None) -> str:
        candidate = self._clean_str(raw_value) or "base"
        if "/" not in candidate and "\\" not in candidate:
            return candidate

        python_path = Path(candidate).expanduser()
        python_path = python_path.resolve() if python_path.is_absolute() else (self.project_root / python_path).resolve()
        name_lower = python_path.name.lower()
        if name_lower.endswith(".py") or "bac_assemble_260112_newformat.py" in name_lower:
            raise ValidationError("运行环境配置错误：这里应填写 Conda 环境名，或选择该环境中的 python 可执行文件，而不是分析脚本路径")
        if not python_path.is_file():
            raise ValidationError(f"运行环境 Python 不存在: {python_path}")
        parts = list(python_path.parts)
        if "envs" in parts:
            env_index = parts.index("envs")
            if env_index + 1 < len(parts):
                return parts[env_index + 1]
        if python_path.name.lower().startswith("python"):
            return "base"
        raise ValidationError("无法从所选路径识别 Conda 环境名，请直接填写环境名，或选择 envs/<环境名>/bin/python")

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

    def _resolve_task_output_path(self, raw_value: Any, task_name: str) -> Path:
        output_root = self._resolve_output_path(raw_value)
        if output_root.exists() and not output_root.is_dir():
            raise ValidationError(f"输出路径不是文件夹: {output_root}")
        safe_task_name = re.sub(r"[\\/]+", "_", self._clean_str(task_name))
        if not safe_task_name:
            raise ValidationError("task_name 不能为空")
        if output_root.name == safe_task_name:
            output_path = output_root
        else:
            output_path = output_root / safe_task_name
        self._validate_ascii_output_path(output_path)
        return output_path

    def _validate_ascii_output_path(self, path: Path) -> None:
        text = str(path)
        if text.isascii():
            return
        invalid_chars = "".join(dict.fromkeys(char for char in text if ord(char) > 127))
        preview = invalid_chars[:12]
        suffix = "..." if len(invalid_chars) > 12 else ""
        raise ValidationError(f"输出路径不能包含非 ASCII 字符，请去掉中文或全角字符: {preview}{suffix}")

    def _resolve_optional_path(self, raw_value: Any, placeholder: str) -> str:
        text = self._clean_str(raw_value)
        if not text:
            return placeholder
        path = Path(text).expanduser()
        path = (self.project_root / path).resolve() if not path.is_absolute() else path.resolve()
        if not path.exists():
            raise ValidationError(f"文件不存在: {path}")
        return str(path)

    def _normalize_rmhost(self, raw_value: Any) -> str:
        text = self._clean_str(raw_value) or "norm"
        if text == "norm":
            return "norm"
        path = Path(text).expanduser()
        path = (self.project_root / path).resolve() if not path.is_absolute() else path.resolve()
        if path.exists():
            return str(path)
        mmi_path = Path(f"{path}.mmi")
        bt2_matches = list(path.parent.glob(f"{path.name}*.bt2*")) if path.parent.exists() else []
        if mmi_path.exists() or bt2_matches:
            return str(path)
        raise ValidationError(f"宿主索引不存在: {path}")

    def _clean_str(self, value: Any) -> str:
        return str(value).strip() if value is not None else ""

    def _normalize_asm_type(self, value: Any) -> str:
        text = self._clean_str(value) or "shortasm"
        if text == "longshortasm":
            text = "shortlongasm"
        if text not in ASM_METHOD_OPTIONS:
            raise ValidationError(f"不支持的组装类型: {text}")
        return text

    def _normalize_analysis_target(self, value: Any) -> str:
        text = self._clean_str(value) or "bacteria"
        if text not in ANALYSIS_TARGET_OPTIONS:
            raise ValidationError(f"不支持的分析对象: {text}")
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
