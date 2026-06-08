from __future__ import annotations

import argparse
import csv
import functools
import json
import os
import subprocess
import zipfile
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from statistics import median
from typing import Iterable


def utc_now_text() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

def write_log(log_path: Path, message: str) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write(f"[{utc_now_text()}] {message}\n")


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, default=str), encoding="utf-8")


def write_tsv(path: Path, columns: list[str], rows: Iterable[Iterable[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(columns)
        for row in rows:
            writer.writerow(list(row))


def guess_delimiter(sample_text: str) -> str:
    if "\t" in sample_text:
        return "\t"
    if "," in sample_text:
        return ","
    return "\t"


def read_delimited_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    raw_text = path.read_text(encoding="utf-8", errors="ignore")
    lines = [line for line in raw_text.splitlines() if line.strip()]
    if not lines:
        raise ValueError(f"文件为空: {path}")
    delimiter = guess_delimiter(lines[0])
    reader = csv.DictReader(lines, delimiter=delimiter)
    fieldnames = [str(name or "").strip() for name in (reader.fieldnames or []) if str(name or "").strip()]
    if not fieldnames:
        raise ValueError(f"文件缺少表头: {path}")
    rows: list[dict[str, str]] = []
    for row in reader:
        normalized = {key: str(value or "").strip() for key, value in row.items() if key}
        if any(normalized.values()):
            rows.append(normalized)
    if not rows:
        raise ValueError(f"文件中没有可用记录: {path}")
    return fieldnames, rows


@dataclass
class MetadataSummary:
    sample_count: int
    sample_id_column: str
    group_column: str
    metadata_columns: list[str]
    preview_rows: list[dict[str, str]]
    group_counts: list[tuple[str, int]]


@dataclass
class TaxonomySummary:
    feature_count: int
    feature_column: str
    taxonomy_columns: list[str]
    preview_rows: list[dict[str, str]]


@dataclass
class DemuxSummary:
    visualization_path: str
    artifact_path: str
    artifact_available: bool
    provenance_type: str
    forward_reads_total: int
    reverse_reads_total: int
    min_forward_reads: int
    max_forward_reads: int
    median_forward_reads: int
    min_reverse_reads: int
    max_reverse_reads: int
    median_reverse_reads: int
    sample_count: int
    suggested_trim_left_f: int
    suggested_trim_left_r: int
    suggested_trunc_len_f: int
    suggested_trunc_len_r: int
    suggested_sampling_depth: int
    preview_rows: list[dict[str, str]]


def read_metadata(metadata_path: Path, group_column: str) -> MetadataSummary:
    fieldnames, rows = read_delimited_rows(metadata_path)
    sample_id_column = fieldnames[0]
    if group_column not in fieldnames:
        raise ValueError(f"元数据文件中未找到分组列：{group_column}")
    group_counts = Counter(str(row.get(group_column) or "未分组").strip() or "未分组" for row in rows)
    return MetadataSummary(
        sample_count=len(rows),
        sample_id_column=sample_id_column,
        group_column=group_column,
        metadata_columns=fieldnames,
        preview_rows=rows[:8],
        group_counts=group_counts.most_common(),
    )


def read_taxonomy(taxonomy_path: Path) -> TaxonomySummary:
    fieldnames, rows = read_delimited_rows(taxonomy_path)
    return TaxonomySummary(
        feature_count=len(rows),
        feature_column=fieldnames[0],
        taxonomy_columns=fieldnames,
        preview_rows=rows[:8],
    )


def collect_input_summary(input_path: Path) -> dict[str, object]:
    if input_path.is_file():
        return {
            "path_type": "file",
            "entry_count": 1,
            "entries": [input_path.name],
            "suffix": "".join(input_path.suffixes) or input_path.suffix,
        }
    entries = sorted(child.name for child in input_path.iterdir()) if input_path.is_dir() else []
    return {
        "path_type": "directory",
        "entry_count": len(entries),
        "entries": entries[:30],
        "suffix": "",
    }


def parse_analysis_list(text: str) -> list[str]:
    allowed = {"alpha", "beta", "lefse", "ml", "network"}
    ordered: list[str] = []
    for item in str(text or "").split(","):
        key = item.strip().lower()
        if key and key in allowed and key not in ordered:
            ordered.append(key)
    return ordered or ["alpha", "beta"]


def resolve_amplicon_files(input_path: Path) -> tuple[Path | None, Path | None]:
    if input_path.is_file():
        if input_path.name.endswith(".qzv"):
            sibling = input_path.with_suffix(".qza")
            return input_path, sibling if sibling.is_file() else None
        if input_path.name.endswith(".qza"):
            sibling = input_path.with_suffix(".qzv")
            return sibling if sibling.is_file() else None, input_path
        return None, None
    if not input_path.is_dir():
        return None, None
    demux_qzv = input_path / "demux.qzv"
    demux_qza = input_path / "demux.qza"
    return (demux_qzv if demux_qzv.is_file() else None, demux_qza if demux_qza.is_file() else None)


def is_amplicon_input(input_path: Path) -> bool:
    demux_qzv, demux_qza = resolve_amplicon_files(input_path)
    return demux_qzv is not None or demux_qza is not None


def _find_qzv_member(names: list[str], suffix: str) -> str:
    for name in names:
        if name.endswith(suffix):
            return name
    raise FileNotFoundError(f"未在 qzv 中找到 {suffix}")


def _read_qzv_text(path: Path, suffix: str) -> str:
    with zipfile.ZipFile(path) as archive:
        member = _find_qzv_member(archive.namelist(), suffix)
        return archive.read(member).decode("utf-8", errors="ignore")


def _read_qzv_tsv(path: Path, suffix: str) -> tuple[list[str], list[dict[str, str]]]:
    raw_text = _read_qzv_text(path, suffix)
    lines = [line for line in raw_text.splitlines() if line.strip()]
    reader = csv.DictReader(lines, delimiter="\t")
    fieldnames = [str(name or "").strip() for name in (reader.fieldnames or []) if str(name or "").strip()]
    rows = [{key: str(value or "").strip() for key, value in row.items() if key} for row in reader]
    return fieldnames, rows


def _read_qzv_quality_summary(path: Path, suffix: str) -> tuple[list[int], dict[str, list[float]]]:
    raw_text = _read_qzv_text(path, suffix)
    lines = [line for line in raw_text.splitlines() if line.strip()]
    reader = csv.reader(lines, delimiter="\t")
    rows = list(reader)
    positions = [int(str(item or "0").strip()) for item in rows[0][1:] if str(item or "").strip()]
    summary: dict[str, list[float]] = {}
    for row in rows[1:]:
        key = str(row[0] or "").strip()
        values = [float(str(item or "0").strip()) for item in row[1:1 + len(positions)]]
        summary[key] = values
    return positions, summary


def _extract_provenance_artifact_type(path: Path) -> str:
    try:
        text = _read_qzv_text(path, "provenance/artifacts/d27a741c-f7e9-48af-ad8a-a479bd89ec9e/metadata.yaml")
        for line in text.splitlines():
            if line.startswith("type:"):
                return line.split(":", 1)[1].strip()
    except Exception:
        pass
    try:
        with zipfile.ZipFile(path) as archive:
            for name in archive.namelist():
                if name.endswith("/metadata.yaml") and "/provenance/artifacts/" in name:
                    text = archive.read(name).decode("utf-8", errors="ignore")
                    for line in text.splitlines():
                        if line.startswith("type:"):
                            return line.split(":", 1)[1].strip()
    except Exception:
        pass
    return "SampleData[PairedEndSequencesWithQuality]"


def _suggest_trunc_len(positions: list[int], summary: dict[str, list[float]]) -> int:
    if not positions:
        return 240
    median_values = summary.get("50%") or summary.get("50")
    lower_quartile = summary.get("25%") or summary.get("25")
    candidate = positions[-1]
    for index, position in enumerate(positions):
        med = median_values[index] if median_values and index < len(median_values) else 0
        q1 = lower_quartile[index] if lower_quartile and index < len(lower_quartile) else med
        if med < 30 or q1 < 25:
            candidate = positions[max(index - 1, 0)]
            break
    return max(120, min(candidate, positions[-1]))


def _round_to_ten(value: int) -> int:
    return max(10, (value // 10) * 10)


def _project_root() -> Path:
    return Path(__file__).resolve().parent


def resolve_classifier_path() -> Path:
    configured = str(os.environ.get("COMMUNITY_QIIME_CLASSIFIER", "")).strip()
    if configured:
        candidate = Path(configured).expanduser()
        if not candidate.is_absolute():
            candidate = (_project_root() / candidate).resolve()
        if candidate.is_file():
            return candidate

    db_dir = _project_root() / "database" / "16s"
    preferred = [
        "silva-138-99-nb-classifier.qza",
        "2024.09.backbone.full-length.nb.sklearn-1.4.2.qza",
        "gtdb_classifier_r220.qza",
    ]
    for name in preferred:
        candidate = db_dir / name
        if candidate.is_file():
            return candidate
    for candidate in sorted(db_dir.glob("*.qza")):
        if candidate.name != "suboptimal-16S-rRNA-classifier.qza":
            return candidate
    raise FileNotFoundError(f"未找到可用的 16S classifier: {db_dir}")


def _qiime_shell_prefix() -> str:
    return (
        'export PATH="$CONDA_PREFIX/bin:$PATH"; '
        'export R_HOME="$CONDA_PREFIX/lib/R"; '
        'export R_LIBS_USER="$CONDA_PREFIX/lib/R/library"; '
        'export MPLCONFIGDIR="/tmp/mpl-qiime2-amplicon"; '
        'export NUMBA_CACHE_DIR="/tmp/numba-qiime2-amplicon"; '
        'mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"'
    )


def run_qiime_command(*, qiime_env: str, command: str, log_path: Path) -> None:
    wrapped = f"{_qiime_shell_prefix()}; {command}"
    write_log(log_path, f"执行命令: {command}")
    completed = subprocess.run(
        ["conda", "run", "-n", qiime_env, "/bin/zsh", "-lc", wrapped],
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.stdout.strip():
        write_log(log_path, completed.stdout.strip())
    if completed.stderr.strip():
        write_log(log_path, completed.stderr.strip())
    if completed.returncode != 0:
        raise subprocess.CalledProcessError(completed.returncode, completed.args, output=completed.stdout, stderr=completed.stderr)


def run_env_command(*, env_name: str, command: str, log_path: Path) -> None:
    wrapped = f"{_qiime_shell_prefix()}; {command}"
    write_log(log_path, f"执行命令: {command}")
    completed = subprocess.run(
        ["conda", "run", "-n", env_name, "/bin/zsh", "-lc", wrapped],
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.stdout.strip():
        write_log(log_path, completed.stdout.strip())
    if completed.stderr.strip():
        write_log(log_path, completed.stderr.strip())
    if completed.returncode != 0:
        raise subprocess.CalledProcessError(completed.returncode, completed.args, output=completed.stdout, stderr=completed.stderr)


def _path_has_content(path: Path) -> bool:
    if path.is_file():
        return path.stat().st_size > 0
    if path.is_dir():
        try:
            next(path.iterdir())
            return True
        except StopIteration:
            return False
    return False


def _targets_ready(targets: Iterable[object]) -> bool:
    normalized = [Path(str(item)).expanduser() for item in targets if str(item or "").strip()]
    if not normalized:
        return False
    return all(_path_has_content(path) for path in normalized)


@functools.lru_cache(maxsize=16)
def _conda_env_exists(env_name: str) -> bool:
    target = str(env_name or "").strip()
    if not target:
        return False
    completed = subprocess.run(
        ["conda", "env", "list", "--json"],
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode != 0:
        return False
    try:
        payload = json.loads(completed.stdout or "{}")
    except json.JSONDecodeError:
        return False
    env_paths = payload.get("envs") if isinstance(payload, dict) else []
    target_lower = target.lower()
    for raw_path in env_paths if isinstance(env_paths, list) else []:
        env_path = Path(str(raw_path or "")).expanduser()
        if env_path.name.lower() == target_lower:
            return True
    return False


@functools.lru_cache(maxsize=32)
def _r_packages_available(env_name: str, package_key: str) -> bool:
    target = str(env_name or "").strip()
    packages = [item.strip() for item in str(package_key or "").split(",") if item.strip()]
    if not target or not packages or not _conda_env_exists(target):
        return False
    package_expr = ", ".join([f"'{name}'" for name in packages])
    completed = subprocess.run(
        [
            "conda",
            "run",
            "-n",
            target,
            "--no-capture-output",
            "Rscript",
            "-e",
            f"pkgs <- c({package_expr}); ok <- all(vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)); quit(status=if (ok) 0 else 1)",
        ],
        text=True,
        capture_output=True,
        check=False,
    )
    return completed.returncode == 0


def resolve_runtime_env(preferred_env: str, *, fallback_envs: Iterable[str], required_r_packages: Iterable[str] | None = None) -> str:
    preferred = str(preferred_env or "").strip()
    candidates: list[str] = []
    for name in [preferred, *list(fallback_envs)]:
        normalized = str(name or "").strip()
        if normalized and normalized not in candidates:
            candidates.append(normalized)
    required = [item.strip() for item in (required_r_packages or []) if item and str(item).strip()]
    if not required:
        for env_name in candidates:
            if _conda_env_exists(env_name):
                return env_name
        return preferred or (candidates[0] if candidates else "")
    package_key = ",".join(required)
    for env_name in candidates:
        if _r_packages_available(env_name, package_key):
            return env_name
    for env_name in candidates:
        if _conda_env_exists(env_name):
            return env_name
    return preferred or (candidates[0] if candidates else "")


def summarize_demux_qzv(demux_qzv_path: Path, demux_qza_path: Path | None) -> DemuxSummary:
    _, count_rows = _read_qzv_tsv(demux_qzv_path, "data/per-sample-fastq-counts.tsv")
    forward_counts = [int(float(str(row.get("forward sequence count") or "0"))) for row in count_rows]
    reverse_counts = [int(float(str(row.get("reverse sequence count") or "0"))) for row in count_rows]
    forward_positions, forward_summary = _read_qzv_quality_summary(demux_qzv_path, "data/forward-seven-number-summaries.tsv")
    reverse_positions, reverse_summary = _read_qzv_quality_summary(demux_qzv_path, "data/reverse-seven-number-summaries.tsv")
    suggested_trunc_len_f = _round_to_ten(_suggest_trunc_len(forward_positions, forward_summary))
    suggested_trunc_len_r = _round_to_ten(_suggest_trunc_len(reverse_positions, reverse_summary))
    minimum_depth = min(forward_counts + reverse_counts) if (forward_counts and reverse_counts) else 300
    suggested_sampling_depth = _round_to_ten(max(100, minimum_depth - 20))
    preview_rows = [
        {
            "sample_id": str(row.get("sample ID") or "").strip(),
            "forward_count": str(row.get("forward sequence count") or "").strip(),
            "reverse_count": str(row.get("reverse sequence count") or "").strip(),
        }
        for row in count_rows[:8]
    ]
    return DemuxSummary(
        visualization_path=str(demux_qzv_path),
        artifact_path=str(demux_qza_path) if demux_qza_path else "",
        artifact_available=bool(demux_qza_path and demux_qza_path.is_file()),
        provenance_type=_extract_provenance_artifact_type(demux_qzv_path),
        forward_reads_total=sum(forward_counts),
        reverse_reads_total=sum(reverse_counts),
        min_forward_reads=min(forward_counts) if forward_counts else 0,
        max_forward_reads=max(forward_counts) if forward_counts else 0,
        median_forward_reads=int(median(forward_counts)) if forward_counts else 0,
        min_reverse_reads=min(reverse_counts) if reverse_counts else 0,
        max_reverse_reads=max(reverse_counts) if reverse_counts else 0,
        median_reverse_reads=int(median(reverse_counts)) if reverse_counts else 0,
        sample_count=len(count_rows),
        suggested_trim_left_f=0,
        suggested_trim_left_r=0,
        suggested_trunc_len_f=suggested_trunc_len_f,
        suggested_trunc_len_r=suggested_trunc_len_r,
        suggested_sampling_depth=suggested_sampling_depth,
        preview_rows=preview_rows,
    )


def build_amplicon_workflow(
    args: argparse.Namespace,
    metadata_path: Path,
    metadata_summary: MetadataSummary,
    demux_summary: DemuxSummary,
    output_dir: Path,
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]], list[dict[str, str]], list[str]]:
    qiime_env = str(os.environ.get("COMMUNITY_QIIME_ENV", "qiime2-amplicon")).strip() or "qiime2-amplicon"
    microeco_env_requested = str(os.environ.get("COMMUNITY_MICROECO_ENV", "microeco")).strip() or "microeco"
    microeco_env = resolve_runtime_env(
        microeco_env_requested,
        fallback_envs=["base", qiime_env],
        required_r_packages=["file2meco", "microeco", "ggplot2"],
    )
    biomarker_env = resolve_runtime_env(
        microeco_env_requested,
        fallback_envs=["base", qiime_env],
        required_r_packages=["file2meco", "microeco", "ggplot2", "randomForest"],
    )
    analyses = parse_analysis_list(args.analyses)
    project_root = _project_root()
    alpha_metrics = [item.strip() for item in str(args.alpha_metrics or "").split(",") if item.strip()]
    classifier = resolve_classifier_path()
    demux_qza = Path(demux_summary.artifact_path) if demux_summary.artifact_path else (output_dir.parent / "demux.qza")
    sample_metadata_qzv = output_dir / "sample-metadata.qzv"
    demux_output_qzv = output_dir / "demux.qzv"
    table_qza = output_dir / "table-dada2.qza"
    rep_seqs_qza = output_dir / "rep-seqs-dada2.qza"
    denoise_stats_qza = output_dir / "denoising-stats-dada2.qza"
    base_transition_qza = output_dir / "base-transition-stats-dada2.qza"
    denoise_stats_qzv = output_dir / "denoising-stats-dada2.qzv"
    table_qzv = output_dir / "table-dada2.qzv"
    rep_seqs_qzv = output_dir / "rep-seqs-dada2.qzv"
    feature_freq_qza = output_dir / "feature-frequencies.qza"
    sample_freq_qza = output_dir / "sample-frequencies.qza"
    taxonomy_qza = output_dir / "taxonomy.qza"
    taxonomy_qzv = output_dir / "taxonomy.qzv"
    taxonomy_export_dir = output_dir / "taxonomy_export"
    taxa_barplot_qzv = output_dir / "taxa-barplot.qzv"
    alpha_rarefaction_qzv = output_dir / "alpha-rarefaction.qzv"
    alpha_metric_outputs = [
        {
            "metric": metric,
            "qza": output_dir / f"alpha-{metric}.qza",
            "qzv": output_dir / f"alpha-{metric}.qzv",
            "export_dir": output_dir / f"alpha-{metric}_export",
        }
        for metric in alpha_metrics
    ]
    microeco_beta_dir = output_dir / "microeco_beta"
    microeco_biomarker_dir = output_dir / "microeco_biomarker"
    microeco_network_dir = output_dir / "microeco_network"
    biomarker_requested = "lefse" in analyses or "ml" in analyses
    network_requested = "network" in analyses

    modules = [
        {
            "key": "amplicon_overview",
            "label": "QIIME2 Amplicon 主流程",
            "status": "completed" if demux_summary.artifact_available else "needs-demux-artifact",
            "description": "按已验证的 QIIME2 流程真实执行 metadata、demux、DADA2、taxonomy、多样性分析，并为 microeco Biomarker 提供标准输入。",
            "runtime": qiime_env,
        },
        {
            "key": "amplicon_denoise",
            "label": "DADA2 去噪与特征表",
            "status": "completed" if demux_summary.artifact_available else "needs-demux-artifact",
            "description": f"基于 demux 质量概览，建议 trunc-len-f={demux_summary.suggested_trunc_len_f}、trunc-len-r={demux_summary.suggested_trunc_len_r}。",
            "runtime": qiime_env,
        },
        {
            "key": "amplicon_taxonomy",
            "label": "Taxonomy 注释与条形图",
            "status": "completed" if demux_summary.artifact_available else "needs-demux-artifact",
            "description": f"使用 {classifier.name} 进行 classify-sklearn 注释，并输出 taxa barplot。",
            "runtime": qiime_env,
        },
    ]
    if "alpha" in analyses:
        modules.append({
            "key": "amplicon_alpha",
            "label": "Alpha 多样性与 Rarefaction",
            "status": "completed" if demux_summary.artifact_available else "needs-demux-artifact",
            "description": f"基于建议采样深度 {demux_summary.suggested_sampling_depth} 计算 {', '.join(alpha_metrics) or 'alpha 指标'} 并输出 rarefaction。",
            "runtime": qiime_env,
        })
    if "beta" in analyses:
        modules.append({
            "key": "amplicon_beta",
            "label": "Beta 多样性",
            "status": "completed" if demux_summary.artifact_available else "needs-demux-artifact",
            "description": f"基于 microeco 计算 {args.beta_metric} 距离、PCoA / NMDS、PERMANOVA 与 betadisper。",
            "runtime": microeco_env,
        })
    if network_requested:
        modules.append({
            "key": "amplicon_network",
            "label": "网络分析",
            "status": "completed" if demux_summary.artifact_available else "needs-demux-artifact",
            "description": "基于 microeco trans_network 构建共现网络，输出节点/边属性、模块统计与网络拓扑摘要。",
            "runtime": microeco_env,
        })
    if biomarker_requested:
        biomarker_modules = []
        if "lefse" in analyses:
            biomarker_modules.append("LEfSe")
        if "ml" in analyses:
            biomarker_modules.append("RF")
        modules.append({
            "key": "amplicon_biomarker",
            "label": "Biomarker 分析",
            "status": "completed" if demux_summary.artifact_available else "needs-demux-artifact",
            "description": f"基于 microeco 运行 {' / '.join(biomarker_modules)}，直接输出差异特征与特征重要性排序。",
            "runtime": biomarker_env,
        })

    command_rows = [
        {"module": "Metadata 预览", "runtime": qiime_env, "env": qiime_env, "command": f"qiime metadata tabulate --m-input-file {metadata_path} --o-visualization {sample_metadata_qzv}", "targets": [sample_metadata_qzv]},
        {"module": "Demux 质控汇总", "runtime": qiime_env, "env": qiime_env, "command": f"qiime demux summarize --i-data {demux_qza} --o-visualization {demux_output_qzv}", "targets": [demux_output_qzv]},
        {
            "module": "DADA2 去噪",
            "runtime": qiime_env,
            "env": qiime_env,
            "command": (
                "qiime dada2 denoise-paired "
                f"--i-demultiplexed-seqs {demux_qza} "
                f"--p-trim-left-f {demux_summary.suggested_trim_left_f} --p-trim-left-r {demux_summary.suggested_trim_left_r} "
                f"--p-trunc-len-f {demux_summary.suggested_trunc_len_f} --p-trunc-len-r {demux_summary.suggested_trunc_len_r} "
                f"--o-table {table_qza} --o-representative-sequences {rep_seqs_qza} "
                f"--o-denoising-stats {denoise_stats_qza} --o-base-transition-stats {base_transition_qza} "
                f"--p-n-threads {max(1, int(args.thread))}"
            ),
            "targets": [table_qza, rep_seqs_qza, denoise_stats_qza, base_transition_qza],
        },
        {"module": "去噪统计", "runtime": qiime_env, "env": qiime_env, "command": f"qiime metadata tabulate --m-input-file {denoise_stats_qza} --o-visualization {denoise_stats_qzv}", "targets": [denoise_stats_qzv]},
        {
            "module": "特征表汇总",
            "runtime": qiime_env,
            "env": qiime_env,
            "command": (
                "qiime feature-table summarize "
                f"--i-table {table_qza} "
                f"--m-metadata-file {metadata_path} "
                f"--o-feature-frequencies {feature_freq_qza} --o-sample-frequencies {sample_freq_qza} "
                f"--o-summary {table_qzv}"
            ),
            "targets": [feature_freq_qza, sample_freq_qza, table_qzv],
        },
        {"module": "代表序列预览", "runtime": qiime_env, "env": qiime_env, "command": f"qiime feature-table tabulate-seqs --i-data {rep_seqs_qza} --o-visualization {rep_seqs_qzv}", "targets": [rep_seqs_qzv]},
        {
            "module": "分类注释",
            "runtime": qiime_env,
            "env": qiime_env,
            "command": (
                "qiime feature-classifier classify-sklearn "
                f"--i-reads {rep_seqs_qza} --i-classifier {classifier} --o-classification {taxonomy_qza}"
            ),
            "targets": [taxonomy_qza],
        },
        {"module": "Taxonomy 预览", "runtime": qiime_env, "env": qiime_env, "command": f"qiime metadata tabulate --m-input-file {taxonomy_qza} --o-visualization {taxonomy_qzv}", "targets": [taxonomy_qzv]},
        {"module": "Taxa 条形图", "runtime": qiime_env, "env": qiime_env, "command": f"qiime taxa barplot --i-table {table_qza} --i-taxonomy {taxonomy_qza} --m-metadata-file {metadata_path} --o-visualization {taxa_barplot_qzv}", "targets": [taxa_barplot_qzv]},
        {"module": "导出 taxonomy", "runtime": qiime_env, "env": qiime_env, "command": f"qiime tools export --input-path {taxonomy_qza} --output-path {taxonomy_export_dir}", "targets": [taxonomy_export_dir / "taxonomy.tsv"]},
    ]
    if "alpha" in analyses:
        command_rows.extend([
            {
                "module": "Alpha 多样性指标",
                "runtime": qiime_env,
                "env": qiime_env,
                "command": "\n".join(
                    [
                        " && ".join(
                            [
                                f"qiime diversity alpha --i-table {table_qza} --p-metric {item['metric']} --o-alpha-diversity {item['qza']}",
                                f"qiime metadata tabulate --m-input-file {item['qza']} --o-visualization {item['qzv']}",
                                f"qiime tools export --input-path {item['qza']} --output-path {item['export_dir']}",
                            ]
                        )
                        for item in alpha_metric_outputs
                    ]
                ),
                "targets": [
                    target
                    for item in alpha_metric_outputs
                    for target in (item["qza"], item["qzv"], item["export_dir"])
                ],
            },
            {
                "module": "Alpha 稀释曲线",
                "runtime": qiime_env,
                "env": qiime_env,
                "command": "qiime diversity alpha-rarefaction "
                f"--i-table {table_qza} "
                f"--p-max-depth {demux_summary.suggested_sampling_depth} "
                f"--m-metadata-file {metadata_path} --o-visualization {alpha_rarefaction_qzv}",
                "targets": [alpha_rarefaction_qzv],
            },
        ])
    if "beta" in analyses:
        command_rows.append({
            "module": "microeco Beta 多样性",
            "runtime": microeco_env,
            "env": microeco_env,
            "command": f"Rscript {project_root / 'scripts' / 'run_microeco_beta.R'} {table_qza} {taxonomy_qza} {metadata_path} {metadata_summary.group_column} {microeco_beta_dir} {args.beta_metric}",
            "targets": [microeco_beta_dir / "run_summary.txt", microeco_beta_dir / "pcoa_plot.png", microeco_beta_dir / "nmds_plot.png"],
        })
    if network_requested:
        command_rows.append({
            "module": "microeco 网络分析",
            "runtime": microeco_env,
            "env": microeco_env,
            "command": f"Rscript {project_root / 'scripts' / 'run_microeco_network.R'} {table_qza} {taxonomy_qza} {metadata_path} {metadata_summary.group_column} {microeco_network_dir} Genus spearman 0.6 0.05 0.0005",
            "targets": [
                microeco_network_dir / "run_summary.txt",
                microeco_network_dir / "network_summary.tsv",
                microeco_network_dir / "node_table.tsv",
                microeco_network_dir / "edge_table.tsv",
                microeco_network_dir / "module_summary.tsv",
            ],
        })
    if biomarker_requested:
        command_rows.append({
            "module": "microeco Biomarker",
            "runtime": biomarker_env,
            "env": biomarker_env,
            "command": f"Rscript {project_root / 'scripts' / 'run_microeco_biomarker.R'} {table_qza} {taxonomy_qza} {metadata_path} {metadata_summary.group_column} {microeco_biomarker_dir} Genus {args.ml_model}",
            "targets": [microeco_biomarker_dir / "run_summary.txt", microeco_biomarker_dir / "lefse_diff.tsv", microeco_biomarker_dir / "rf_importance.tsv"],
        })

    commands = [
        {
            **row,
            "status": (
                "ready"
                if demux_summary.artifact_available and _targets_ready(row.get("targets") or [])
                else ("planned" if demux_summary.artifact_available else "needs-demux-artifact")
            ),
        }
        for row in command_rows
    ]

    outputs = [
        {"label": "分析总览", "status": "ready", "path": str(output_dir / "community_summary.json")},
        {"label": "命令规划表", "status": "ready", "path": str(output_dir / "community_command_plan.tsv")},
        {"label": "元数据预览", "status": "ready", "path": str(output_dir / "community_metadata_preview.tsv")},
        {"label": "demux 汇总", "status": "ready", "path": str(output_dir / "community_demux_summary.tsv")},
        {"label": "教程流程清单", "status": "ready", "path": str(output_dir / "community_tutorial_workflow.tsv")},
        {"label": "sample-metadata.qzv", "status": "planned" if demux_summary.artifact_available else "needs-demux-artifact", "path": str(sample_metadata_qzv)},
        {"label": "demux.qzv", "status": "planned" if demux_summary.artifact_available else "needs-demux-artifact", "path": str(demux_output_qzv)},
        {"label": "table-dada2.qza", "status": "planned" if demux_summary.artifact_available else "needs-demux-artifact", "path": str(table_qza)},
        {"label": "rep-seqs-dada2.qza", "status": "planned" if demux_summary.artifact_available else "needs-demux-artifact", "path": str(rep_seqs_qza)},
        {"label": "table-dada2.qzv", "status": "planned" if demux_summary.artifact_available else "needs-demux-artifact", "path": str(table_qzv)},
        {"label": "taxonomy.qzv", "status": "planned" if demux_summary.artifact_available else "needs-demux-artifact", "path": str(taxonomy_qzv)},
        {"label": "taxa-barplot.qzv", "status": "planned" if demux_summary.artifact_available else "needs-demux-artifact", "path": str(taxa_barplot_qzv)},
    ]
    if "alpha" in analyses:
        outputs.extend([
            {"label": "alpha-rarefaction.qzv", "status": "planned" if demux_summary.artifact_available else "needs-demux-artifact", "path": str(alpha_rarefaction_qzv)},
            *[
                {"label": f"alpha-{item['metric']}.qzv", "status": "planned" if demux_summary.artifact_available else "needs-demux-artifact", "path": str(item["qzv"])}
                for item in alpha_metric_outputs
            ],
        ])
    if "beta" in analyses:
        outputs.append({"label": "microeco_beta", "status": "planned" if demux_summary.artifact_available else "needs-demux-artifact", "path": str(microeco_beta_dir)})
    if network_requested:
        outputs.append({"label": "microeco_network", "status": "planned" if demux_summary.artifact_available else "needs-demux-artifact", "path": str(microeco_network_dir)})
    if biomarker_requested:
        outputs.append({"label": "microeco_biomarker", "status": "planned" if demux_summary.artifact_available else "needs-demux-artifact", "path": str(microeco_biomarker_dir)})

    tutorial_steps = [
        {"step": "1", "stage": "Metadata", "tutorial_anchor": "metadata tabulate", "status": "planned"},
        {"step": "2", "stage": "Demux", "tutorial_anchor": "demux summarize", "status": "ready" if demux_summary.visualization_path else "missing"},
        {"step": "3", "stage": "Denoise", "tutorial_anchor": "dada2 denoise-paired", "status": "planned" if demux_summary.artifact_available else "needs-demux-artifact"},
        {"step": "4", "stage": "Feature table", "tutorial_anchor": "feature-table summarize / tabulate-seqs", "status": "completed" if demux_summary.artifact_available else "needs-demux-artifact"},
        {"step": "5", "stage": "Taxonomy", "tutorial_anchor": "classify-sklearn / taxa barplot", "status": "completed" if demux_summary.artifact_available else "needs-demux-artifact"},
    ]
    if "alpha" in analyses:
        tutorial_steps.append({"step": str(len(tutorial_steps) + 1), "stage": "Alpha", "tutorial_anchor": "alpha-rarefaction", "status": "completed" if demux_summary.artifact_available else "needs-demux-artifact"})
    if "beta" in analyses:
        tutorial_steps.append({"step": str(len(tutorial_steps) + 1), "stage": "Beta", "tutorial_anchor": "microeco beta", "status": "completed" if demux_summary.artifact_available else "needs-demux-artifact"})
    if network_requested:
        tutorial_steps.append({"step": str(len(tutorial_steps) + 1), "stage": "Network", "tutorial_anchor": "microeco trans_network", "status": "completed" if demux_summary.artifact_available else "needs-demux-artifact"})
    if biomarker_requested:
        tutorial_steps.append({"step": str(len(tutorial_steps) + 1), "stage": "Biomarker", "tutorial_anchor": "microeco lefse / rf", "status": "completed" if demux_summary.artifact_available else "needs-demux-artifact"})

    notes = [
        f"当前流程会按勾选模块规划模块；实际执行时会先检查各步骤最终输出是否已存在，存在则自动跳过，仅补跑缺失步骤。本次模块为 {', '.join(analyses)}。",
        f"当前 demux 可视化来自 {Path(demux_summary.visualization_path).name}，建议 trunc-len-f={demux_summary.suggested_trunc_len_f}、trunc-len-r={demux_summary.suggested_trunc_len_r}。",
        f"默认 classifier 为 {classifier.name}；Biomarker 模块默认使用 LEfSe 与 {args.ml_model}，Beta 模块使用 {args.beta_metric}，Network 模块使用 spearman 相关网络（COR_cut=0.6, P<=0.05），Alpha 模块计算 {', '.join(alpha_metrics) or 'shannon'}。",
    ]
    if microeco_env != microeco_env_requested:
        notes.append(f"未检测到可用的 microeco 专用环境 {microeco_env_requested}，Beta 已自动回退到 {microeco_env}。")
    if biomarker_env != microeco_env_requested:
        notes.append(f"Biomarker 所需 R 包未在 {microeco_env_requested} 中满足，已自动改用 {biomarker_env} 执行。")
    if not demux_summary.artifact_available:
        notes.append("当前输入目录只有 demux.qzv，可用于质量评估与参数推荐；若要实际运行 DADA2，请补充同目录 demux.qza。")
    return modules, commands, outputs, tutorial_steps, notes


def build_abundance_workflow(
    args: argparse.Namespace,
    metadata_summary: MetadataSummary,
    taxonomy_summary: TaxonomySummary,
    output_dir: Path,
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]], list[dict[str, str]], list[str]]:
    qiime_env = str(os.environ.get("COMMUNITY_QIIME_ENV", "qiime2")).strip() or "qiime2"
    microeco_env = str(os.environ.get("COMMUNITY_MICROECO_ENV", "microeco")).strip() or "microeco"
    analyses = parse_analysis_list(args.analyses)
    modules: list[dict[str, str]] = []
    commands: list[dict[str, str]] = []
    if "alpha" in analyses:
        modules.append({"key": "alpha", "label": "Alpha 多样性", "status": "planned", "description": f"基于 {args.alpha_metrics} 计算组内多样性。", "runtime": qiime_env})
        commands.append({"module": "Alpha 多样性", "runtime": qiime_env, "command": f"qiime diversity alpha --p-metric {args.alpha_metrics.split(',')[0].strip()} ...", "status": "planned"})
    if "beta" in analyses:
        modules.append({"key": "beta", "label": "Beta 多样性", "status": "planned", "description": f"基于 {args.beta_metric} 距离做 PCoA / NMDS。", "runtime": qiime_env})
        commands.append({"module": "Beta 多样性", "runtime": qiime_env, "command": f"qiime diversity beta-group-significance --p-method {args.beta_metric} ...", "status": "planned"})
    if "network" in analyses:
        modules.append({"key": "network", "label": "网络分析", "status": "planned", "description": "基于 microeco trans_network 构建共现网络并输出节点/边属性。", "runtime": microeco_env})
        commands.append({"module": "Network 共现网络", "runtime": microeco_env, "command": 'Rscript -e "library(microeco); # trans_network"' , "status": "planned"})
    if "lefse" in analyses:
        modules.append({"key": "lefse", "label": "Biomarker · LEfSe", "status": "planned", "description": "基于 microeco 筛选差异特征并输出 LDA 结果。", "runtime": microeco_env})
        commands.append({"module": "LEfSe 差异特征", "runtime": microeco_env, "command": 'Rscript -e "library(microeco); ..."', "status": "planned"})
    if "ml" in analyses:
        modules.append({"key": "ml", "label": "Biomarker · RF", "status": "planned", "description": f"基于 microeco 使用 {args.ml_model} 做分组判别与特征排序。", "runtime": microeco_env})
        commands.append({"module": "RF Biomarker", "runtime": microeco_env, "command": f'Rscript -e "library(microeco); # method = rf, model = {args.ml_model}"', "status": "planned"})
    outputs = [
        {"label": "分析总览", "status": "ready", "path": str(output_dir / "community_summary.json")},
        {"label": "命令规划表", "status": "ready", "path": str(output_dir / "community_command_plan.tsv")},
        {"label": "元数据预览", "status": "ready", "path": str(output_dir / "community_metadata_preview.tsv")},
        {"label": "taxonomy 预览", "status": "ready", "path": str(output_dir / "community_taxonomy_preview.tsv")},
    ] + [
        {"label": module["label"], "status": "planned", "path": str(output_dir / f"{module['key']}_results")}
        for module in modules
    ]
    tutorial_steps = [{"step": str(index + 1), "stage": module["label"], "tutorial_anchor": "abundance fallback", "status": "planned"} for index, module in enumerate(modules)]
    notes = [
        "当前输入不包含 demux.qzv/qza，因此回退到丰度表 + taxonomy 的群落分析骨架模式。",
        f"taxonomy 表共 {taxonomy_summary.feature_count} 条记录，统计层级为 {args.taxonomy_level}。",
    ]
    return modules, commands, outputs, tutorial_steps, notes


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Community analysis scaffold with QIIME2 amplicon workflow priority.")
    parser.add_argument("--input", required=True, help="输入目录、demux.qza/qzv 所在目录或丰度表目录")
    parser.add_argument("--metadata", required=True, help="样本元数据文件（csv/tsv）")
    parser.add_argument("--taxonomy", default="", help="丰度表模式下使用的 taxonomy 注释文件（csv/tsv）")
    parser.add_argument("--output", required=True, help="输出目录")
    parser.add_argument("--thread", type=int, default=8, help="线程数")
    parser.add_argument("--group-column", required=True, help="元数据中的分组列")
    parser.add_argument("--taxonomy-level", default="genus", help="统计层级，例如 phylum/genus/species")
    parser.add_argument("--normalization", default="relative", help="标准化方式，例如 relative/clr")
    parser.add_argument("--analyses", default="alpha,beta,lefse,ml,network", help="分析模块列表，逗号分隔")
    parser.add_argument("--alpha-metrics", default="shannon,simpson,chao1", help="alpha 多样性指标")
    parser.add_argument("--beta-metric", default="braycurtis", help="beta 多样性距离")
    parser.add_argument("--ml-model", default="random_forest", help="机器学习模型")
    return parser.parse_args()


def mark_output_statuses(outputs: list[dict[str, str]]) -> None:
    for item in outputs:
        path = Path(str(item.get("path") or "").strip())
        status = str(item.get("status") or "").strip().lower()
        if status in {"ready", "needs-demux-artifact"}:
            continue
        item["status"] = "ready" if path.exists() else "missing"


def execute_amplicon_workflow(
    *,
    qiime_env: str,
    commands: list[dict[str, str]],
    outputs: list[dict[str, str]],
    log_path: Path,
) -> None:
    for item in commands:
        command = str(item.get("command") or "").strip()
        if not command:
            item["status"] = "skipped"
            continue
        targets = item.get("targets") or []
        if _targets_ready(targets):
            write_log(log_path, f"跳过命令: {item.get('module') or '未命名步骤'}；检测到目标产物已存在")
            item["status"] = "skipped-existing"
            continue
        env_name = str(item.get("env") or qiime_env).strip() or qiime_env
        run_env_command(env_name=env_name, command=command, log_path=log_path)
        item["status"] = "ready" if _targets_ready(targets) else "missing"
    mark_output_statuses(outputs)


def main() -> int:
    args = parse_args()
    input_path = Path(args.input).expanduser().resolve()
    metadata_path = Path(args.metadata).expanduser().resolve()
    taxonomy_path = Path(args.taxonomy).expanduser().resolve() if str(args.taxonomy or "").strip() else None
    output_dir = Path(args.output).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    log_path = output_dir / "community_pipeline.log"

    if not input_path.exists():
        raise FileNotFoundError(f"输入路径不存在: {input_path}")
    if not metadata_path.is_file():
        raise FileNotFoundError(f"元数据文件不存在: {metadata_path}")

    write_log(log_path, "开始群落分析流程初始化")
    metadata_summary = read_metadata(metadata_path, args.group_column)
    input_summary = collect_input_summary(input_path)
    workflow_mode = "amplicon" if is_amplicon_input(input_path) else "abundance"

    metadata_preview_rows = [
        [str(index + 1), row.get(metadata_summary.sample_id_column, ""), row.get(metadata_summary.group_column, "")]
        for index, row in enumerate(metadata_summary.preview_rows)
    ]
    write_tsv(output_dir / "community_metadata_preview.tsv", ["序号", metadata_summary.sample_id_column, metadata_summary.group_column], metadata_preview_rows)

    summary_payload: dict[str, object] = {
        "report_kind": "community_meta_ecology",
        "generated_at": utc_now_text(),
        "workflow_mode": workflow_mode,
        "input_summary": input_summary,
        "metadata_summary": {
            "sample_count": metadata_summary.sample_count,
            "sample_id_column": metadata_summary.sample_id_column,
            "group_column": metadata_summary.group_column,
            "metadata_columns": metadata_summary.metadata_columns,
            "group_counts": [{"name": name, "count": count} for name, count in metadata_summary.group_counts[:12]],
        },
        "parameters": {
            "thread": args.thread,
            "taxonomy_level": args.taxonomy_level,
            "normalization": args.normalization,
            "alpha_metrics": [item.strip() for item in args.alpha_metrics.split(",") if item.strip()],
            "beta_metric": args.beta_metric,
            "ml_model": args.ml_model,
        },
        "source_inputs": {
            "metadata": str(metadata_path),
            "taxonomy": str(taxonomy_path) if taxonomy_path else "",
        },
    }

    if workflow_mode == "amplicon":
        demux_qzv_path, demux_qza_path = resolve_amplicon_files(input_path)
        if demux_qzv_path is None:
            raise FileNotFoundError(f"未在输入目录中找到 demux.qzv: {input_path}")
        demux_summary = summarize_demux_qzv(demux_qzv_path, demux_qza_path)
        modules, commands, outputs, tutorial_steps, notes = build_amplicon_workflow(args, metadata_path, metadata_summary, demux_summary, output_dir)
        qiime_env = str(os.environ.get("COMMUNITY_QIIME_ENV", "qiime2-amplicon")).strip() or "qiime2-amplicon"
        write_tsv(
            output_dir / "community_demux_summary.tsv",
            ["sample_id", "forward_count", "reverse_count"],
            ((row["sample_id"], row["forward_count"], row["reverse_count"]) for row in demux_summary.preview_rows),
        )
        summary_payload["demux_summary"] = {
            "visualization_path": demux_summary.visualization_path,
            "artifact_path": demux_summary.artifact_path,
            "artifact_available": demux_summary.artifact_available,
            "provenance_type": demux_summary.provenance_type,
            "sample_count": demux_summary.sample_count,
            "forward_reads_total": demux_summary.forward_reads_total,
            "reverse_reads_total": demux_summary.reverse_reads_total,
            "min_forward_reads": demux_summary.min_forward_reads,
            "median_forward_reads": demux_summary.median_forward_reads,
            "max_forward_reads": demux_summary.max_forward_reads,
            "min_reverse_reads": demux_summary.min_reverse_reads,
            "median_reverse_reads": demux_summary.median_reverse_reads,
            "max_reverse_reads": demux_summary.max_reverse_reads,
            "suggested_trim_left_f": demux_summary.suggested_trim_left_f,
            "suggested_trim_left_r": demux_summary.suggested_trim_left_r,
            "suggested_trunc_len_f": demux_summary.suggested_trunc_len_f,
            "suggested_trunc_len_r": demux_summary.suggested_trunc_len_r,
            "suggested_sampling_depth": demux_summary.suggested_sampling_depth,
        }
        summary_payload["tutorial_reference"] = {
            "name": "QIIME 2 gut-to-soil",
            "url": "https://amplicon-docs.qiime2.org/en/stable/tutorials/gut-to-soil.html",
            "workflow": [
                "metadata tabulate",
                "demux summarize",
                "dada2 denoise-paired",
                "feature-table summarize/tabulate-seqs",
                "classify-sklearn",
                "alpha-rarefaction",
                "taxa barplot",
                "taxa collapse",
                "microeco biomarker (lefse / rf)",
                "microeco trans_network",
            ],
        }
        write_log(log_path, f"已识别为 QIIME2 amplicon 模式，demux 可视化样本数 {demux_summary.sample_count}")
        write_log(log_path, f"建议 trunc-len-f={demux_summary.suggested_trunc_len_f}, trunc-len-r={demux_summary.suggested_trunc_len_r}, max-depth={demux_summary.suggested_sampling_depth}")
        if demux_summary.artifact_available:
            execute_amplicon_workflow(qiime_env=qiime_env, commands=commands, outputs=outputs, log_path=log_path)
            exported_taxonomy = output_dir / "taxonomy_export" / "taxonomy.tsv"
            if exported_taxonomy.is_file():
                taxonomy_summary = read_taxonomy(exported_taxonomy)
                summary_payload["taxonomy_summary"] = {
                    "feature_count": taxonomy_summary.feature_count,
                    "feature_column": taxonomy_summary.feature_column,
                    "taxonomy_columns": taxonomy_summary.taxonomy_columns,
                    "source_path": str(exported_taxonomy),
                }
            write_log(log_path, "QIIME2 amplicon 流程已真实执行完成")
        else:
            write_log(log_path, "缺少 demux.qza，当前仅生成流程规划与参数建议")
    else:
        if not taxonomy_path or not taxonomy_path.is_file():
            raise FileNotFoundError("当前输入未检测到 demux.qzv/qza，因此必须提供 taxonomy 文件")
        taxonomy_summary = read_taxonomy(taxonomy_path)
        taxonomy_preview_columns = [taxonomy_summary.feature_column] + taxonomy_summary.taxonomy_columns[1:min(len(taxonomy_summary.taxonomy_columns), 6)]
        taxonomy_preview_rows = [[row.get(column, "") for column in taxonomy_preview_columns] for row in taxonomy_summary.preview_rows]
        write_tsv(output_dir / "community_taxonomy_preview.tsv", taxonomy_preview_columns, taxonomy_preview_rows)
        modules, commands, outputs, tutorial_steps, notes = build_abundance_workflow(args, metadata_summary, taxonomy_summary, output_dir)
        summary_payload["taxonomy_summary"] = {
            "feature_count": taxonomy_summary.feature_count,
            "feature_column": taxonomy_summary.feature_column,
            "taxonomy_columns": taxonomy_summary.taxonomy_columns,
            "source_path": str(taxonomy_path),
        }
        write_log(log_path, f"已识别为丰度表模式，taxonomy 条目 {taxonomy_summary.feature_count}")

    write_tsv(
        output_dir / "community_command_plan.tsv",
        ["模块", "运行环境", "状态", "命令"],
        ((row["module"], row["runtime"], row.get("status") or "planned", row["command"]) for row in commands),
    )
    write_tsv(
        output_dir / "community_tutorial_workflow.tsv",
        ["步骤", "阶段", "教程锚点", "状态"],
        ((row["step"], row["stage"], row["tutorial_anchor"], row["status"]) for row in tutorial_steps),
    )

    summary_payload["modules"] = modules
    summary_payload["commands"] = commands
    summary_payload["outputs"] = outputs
    summary_payload["tutorial_steps"] = tutorial_steps
    summary_payload["notes"] = notes
    write_json(output_dir / "community_summary.json", summary_payload)

    write_log(log_path, f"已解析元数据，共纳入 {metadata_summary.sample_count} 个样本")
    write_log(log_path, f"已生成 {len(commands)} 条流程命令记录")
    write_log(log_path, "群落分析流程结果已写出")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
