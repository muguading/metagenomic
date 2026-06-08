from __future__ import annotations

import argparse
import csv
import gzip
import math
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

from metagenomic_refactor.common import conda_run_prefix


class MetaQuantError(RuntimeError):
    """Raised when the metagenome/metatranscriptome quantification workflow fails."""


@dataclass(frozen=True)
class MetaQuantSample:
    sample: str
    bins_dir: Path
    dna_fastq1: Path
    dna_fastq2: Path
    rna_fastq1: Path
    rna_fastq2: Path


@dataclass(frozen=True)
class MetaQuantConfig:
    outdir: Path
    threads: int = 16
    conda_env: str = "meta_quant"
    fastp_min_length: int = 50
    coverm_min_covered_fraction: float = 0.0
    coverm_min_read_percent_identity: float = 0.95
    coverm_min_read_aligned_percent: float = 0.75
    skip_ribodetector: bool = False
    force: bool = False


@dataclass(frozen=True)
class GeneRecord:
    gene_id: str
    bin_id: str
    source_gene_id: str
    nt_length: int


@dataclass(frozen=True)
class ScatterStats:
    point_count: int
    r2: float
    slope: float
    intercept: float


def _sanitize_sample_name(sample: str) -> str:
    return sample.strip().replace("/", "_").replace(" ", "_")


def _ensure_file(path: Path, label: str) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        raise MetaQuantError(f"{label}不存在: {resolved}")
    return resolved


def _ensure_dir(path: Path, label: str) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_dir():
        raise MetaQuantError(f"{label}不存在: {resolved}")
    return resolved


def _iter_bin_fastas(bins_dir: Path) -> list[Path]:
    fasta_files: list[Path] = []
    for pattern in ("*.fa", "*.fna", "*.fasta", "*.fa.gz", "*.fna.gz", "*.fasta.gz"):
        fasta_files.extend(sorted(bins_dir.glob(pattern)))
    if not fasta_files:
        raise MetaQuantError(f"目录下未找到bin FASTA文件: {bins_dir}")
    return fasta_files


def _opener(path: Path):
    return gzip.open if path.suffix == ".gz" else open


def _read_fasta_records(path: Path) -> Iterable[tuple[str, str]]:
    opener = _opener(path)
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        header = ""
        chunks: list[str] = []
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    yield header, "".join(chunks)
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
        if header:
            yield header, "".join(chunks)


def load_manifest(manifest_path: str | Path) -> list[MetaQuantSample]:
    manifest = _ensure_file(Path(manifest_path), "样本表")
    required = {"sample", "bins_dir", "dna_fastq1", "dna_fastq2", "rna_fastq1", "rna_fastq2"}
    samples: list[MetaQuantSample] = []
    with manifest.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            raise MetaQuantError(f"样本表缺少必要列，至少需要: {', '.join(sorted(required))}")
        for row_num, row in enumerate(reader, start=2):
            sample_raw = str(row.get("sample") or "").strip()
            if not sample_raw:
                continue
            bins_dir_raw = str(row.get("bins_dir") or "").strip()
            dna_fastq1_raw = str(row.get("dna_fastq1") or "").strip()
            dna_fastq2_raw = str(row.get("dna_fastq2") or "").strip()
            rna_fastq1_raw = str(row.get("rna_fastq1") or "").strip()
            rna_fastq2_raw = str(row.get("rna_fastq2") or "").strip()
            if not bins_dir_raw:
                raise MetaQuantError(f"第{row_num}行缺少bins_dir列")
            if not dna_fastq1_raw or not dna_fastq2_raw:
                raise MetaQuantError(f"第{row_num}行缺少宏基因组双端FASTQ")
            if not rna_fastq1_raw or not rna_fastq2_raw:
                raise MetaQuantError(f"第{row_num}行缺少宏转录组双端FASTQ")
            samples.append(
                MetaQuantSample(
                    sample=_sanitize_sample_name(sample_raw),
                    bins_dir=_ensure_dir(Path(bins_dir_raw), f"{sample_raw}的dRep bin目录"),
                    dna_fastq1=_ensure_file(Path(dna_fastq1_raw), f"{sample_raw}的宏基因组fastq1"),
                    dna_fastq2=_ensure_file(Path(dna_fastq2_raw), f"{sample_raw}的宏基因组fastq2"),
                    rna_fastq1=_ensure_file(Path(rna_fastq1_raw), f"{sample_raw}的宏转录组fastq1"),
                    rna_fastq2=_ensure_file(Path(rna_fastq2_raw), f"{sample_raw}的宏转录组fastq2"),
                )
            )
    if not samples:
        raise MetaQuantError(f"样本表中没有有效样本: {manifest}")
    return samples


def _tool_prefix(cfg: MetaQuantConfig) -> list[str]:
    if not cfg.conda_env:
        return []
    return conda_run_prefix(cfg.conda_env)


def _ribodetector_prefix() -> list[str]:
    return conda_run_prefix("Ribodetector")


def _coverm_prefix() -> list[str]:
    return conda_run_prefix("mag_aux")


def _tool_cmd(cfg: MetaQuantConfig, cmd: Sequence[str]) -> list[str]:
    return _tool_prefix(cfg) + list(cmd)


def _run_command(
    cfg: MetaQuantConfig,
    cmd: Sequence[str],
    stdout_log: Path,
    stderr_log: Path,
    cwd: Path | None = None,
) -> None:
    full_cmd = _tool_cmd(cfg, cmd)
    stdout_log.parent.mkdir(parents=True, exist_ok=True)
    stderr_log.parent.mkdir(parents=True, exist_ok=True)
    with stdout_log.open("w", encoding="utf-8") as stdout_handle, stderr_log.open(
        "w", encoding="utf-8"
    ) as stderr_handle:
        completed = subprocess.run(
            full_cmd,
            stdout=stdout_handle,
            stderr=stderr_handle,
            text=True,
            cwd=str(cwd) if cwd else None,
        )
    if completed.returncode != 0:
        raise MetaQuantError(
            f"命令执行失败: {' '.join(full_cmd)}\n请检查日志: {stdout_log} 和 {stderr_log}"
        )


def _run_prefixed_command(
    prefix: Sequence[str],
    cmd: Sequence[str],
    stdout_log: Path,
    stderr_log: Path,
    cwd: Path | None = None,
) -> None:
    full_cmd = list(prefix) + list(cmd)
    stdout_log.parent.mkdir(parents=True, exist_ok=True)
    stderr_log.parent.mkdir(parents=True, exist_ok=True)
    with stdout_log.open("w", encoding="utf-8") as stdout_handle, stderr_log.open(
        "w", encoding="utf-8"
    ) as stderr_handle:
        completed = subprocess.run(
            full_cmd,
            stdout=stdout_handle,
            stderr=stderr_handle,
            text=True,
            cwd=str(cwd) if cwd else None,
        )
    if completed.returncode != 0:
        raise MetaQuantError(
            f"命令执行失败: {' '.join(full_cmd)}\n请检查日志: {stdout_log} 和 {stderr_log}"
        )


def _run_shell_command(
    command: str,
    stdout_log: Path,
    stderr_log: Path,
    cwd: Path | None = None,
) -> None:
    stdout_log.parent.mkdir(parents=True, exist_ok=True)
    stderr_log.parent.mkdir(parents=True, exist_ok=True)
    with stdout_log.open("w", encoding="utf-8") as stdout_handle, stderr_log.open(
        "w", encoding="utf-8"
    ) as stderr_handle:
        completed = subprocess.run(
            ["bash", "-lc", f"set -euo pipefail; {command}"],
            stdout=stdout_handle,
            stderr=stderr_handle,
            text=True,
            cwd=str(cwd) if cwd else None,
        )
    if completed.returncode != 0:
        raise MetaQuantError(
            f"命令执行失败: {command}\n请检查日志: {stdout_log} 和 {stderr_log}"
        )


def _safe_clear_dir(path: Path) -> None:
    if not path.exists():
        return
    for item in path.iterdir():
        if item.is_dir() and not item.is_symlink():
            shutil.rmtree(item)
        else:
            item.unlink()


def _normalize_bin_name(path: Path) -> str:
    name = path.name
    for suffix in (".gz", ".fasta", ".fna", ".fa"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    return _sanitize_sample_name(name)


def _bin_gene_prediction_complete(bin_dir: Path, bin_id: str) -> bool:
    return all(
        path.is_file()
        for path in (
            bin_dir / f"{bin_id}.genes.fna",
            bin_dir / f"{bin_id}.proteins.faa",
            bin_dir / f"{bin_id}.genes.gff",
        )
    )


def _predict_genes_for_bins(sample: MetaQuantSample, cfg: MetaQuantConfig, sample_out: Path) -> tuple[Path, Path]:
    catalog_dir = sample_out / "gene_catalog"
    pred_dir = catalog_dir / "predictions"
    done_flag = catalog_dir / ".done"
    combined_cds = catalog_dir / "combined_cds.fna"
    gene_meta_tsv = catalog_dir / "gene_metadata.tsv"

    if cfg.force:
        _safe_clear_dir(catalog_dir)
    if done_flag.exists() and combined_cds.is_file() and gene_meta_tsv.is_file():
        return combined_cds, gene_meta_tsv

    pred_dir.mkdir(parents=True, exist_ok=True)
    gene_records: list[GeneRecord] = []
    combined_handle = combined_cds.open("w", encoding="utf-8")
    try:
        for fasta_path in _iter_bin_fastas(sample.bins_dir):
            bin_id = _normalize_bin_name(fasta_path)
            bin_dir = pred_dir / bin_id
            bin_dir.mkdir(parents=True, exist_ok=True)
            cds_fna = bin_dir / f"{bin_id}.genes.fna"
            proteins_faa = bin_dir / f"{bin_id}.proteins.faa"
            genes_gff = bin_dir / f"{bin_id}.genes.gff"
            if cfg.force or not _bin_gene_prediction_complete(bin_dir, bin_id):
                _run_command(
                    cfg,
                    [
                        "prodigal",
                        "-i",
                        str(fasta_path),
                        "-d",
                        str(cds_fna),
                        "-a",
                        str(proteins_faa),
                        "-o",
                        str(genes_gff),
                        "-f",
                        "gff",
                        "-p",
                        "meta",
                    ],
                    bin_dir / "prodigal.stdout.log",
                    bin_dir / "prodigal.stderr.log",
                )
            for source_gene_id, seq in _read_fasta_records(cds_fna):
                source_id = source_gene_id.split()[0]
                gene_id = f"{bin_id}|{source_id}"
                combined_handle.write(f">{gene_id}\n{seq}\n")
                gene_records.append(
                    GeneRecord(
                        gene_id=gene_id,
                        bin_id=bin_id,
                        source_gene_id=source_id,
                        nt_length=len(seq),
                    )
                )
    finally:
        combined_handle.close()

    if not gene_records:
        raise MetaQuantError(f"{sample.sample}未能从bin中预测到任何CDS: {sample.bins_dir}")

    with gene_meta_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["gene_id", "bin_id", "source_gene_id", "nt_length"])
        for record in gene_records:
            writer.writerow([record.gene_id, record.bin_id, record.source_gene_id, record.nt_length])

    done_flag.touch()
    return combined_cds, gene_meta_tsv


def _run_fastp_pair(
    cfg: MetaQuantConfig,
    fastq1: Path,
    fastq2: Path,
    out1: Path,
    out2: Path,
    json_path: Path,
    html_path: Path,
    log_prefix: Path,
) -> tuple[Path, Path]:
    if cfg.force:
        for stale in (out1, out2, json_path, html_path):
            if stale.exists():
                stale.unlink()
    if out1.is_file() and out2.is_file() and json_path.is_file():
        return out1, out2
    _run_command(
        cfg,
        [
            "fastp",
            "--in1",
            str(fastq1),
            "--in2",
            str(fastq2),
            "--out1",
            str(out1),
            "--out2",
            str(out2),
            "--thread",
            str(cfg.threads),
            "--length_required",
            str(cfg.fastp_min_length),
            "--detect_adapter_for_pe",
            "--json",
            str(json_path),
            "--html",
            str(html_path),
        ],
        log_prefix.with_suffix(".stdout.log"),
        log_prefix.with_suffix(".stderr.log"),
    )
    return out1, out2


def _run_ribodetector_pair(
    cfg: MetaQuantConfig,
    fastq1: Path,
    fastq2: Path,
    out1: Path,
    out2: Path,
    log_prefix: Path,
) -> tuple[Path, Path]:
    if cfg.force:
        for stale in (out1, out2):
            if stale.exists():
                stale.unlink()
    if out1.is_file() and out2.is_file():
        return out1, out2

    ribodetector_cmd = "ribodetector" if _server_has_gpu() else "ribodetector_cpu"
    _run_prefixed_command(
        _ribodetector_prefix(),
        [
            ribodetector_cmd,
            "-t",
            str(cfg.threads),
            "-e",
            "rrna",
            "-i",
            str(fastq1),
            str(fastq2),
            "-o",
            str(out1),
            str(out2),
        ],
        log_prefix.with_suffix(".stdout.log"),
        log_prefix.with_suffix(".stderr.log"),
    )
    return out1, out2


def _server_has_gpu() -> bool:
    try:
        completed = subprocess.run(
            ["nvidia-smi", "-L"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except OSError:
        return False
    return completed.returncode == 0 and bool(completed.stdout.strip())


def _run_coverm_pair(
    cfg: MetaQuantConfig,
    reference_fna: Path,
    fastq1: Path,
    fastq2: Path,
    output_tsv: Path,
    log_prefix: Path,
) -> Path:
    if cfg.force and output_tsv.exists():
        output_tsv.unlink()
    if output_tsv.is_file():
        return output_tsv
    _run_prefixed_command(
        _coverm_prefix(),
        [
            "coverm",
            "contig",
            "--reference",
            str(reference_fna),
            "--coupled",
            str(fastq1),
            str(fastq2),
            "--threads",
            str(cfg.threads),
            "--min-covered-fraction",
            str(cfg.coverm_min_covered_fraction),
            "--min-read-percent-identity",
            str(cfg.coverm_min_read_percent_identity),
            "--min-read-aligned-percent",
            str(cfg.coverm_min_read_aligned_percent),
            "--methods",
            "count",
            "tpm",
            "--output-file",
            str(output_tsv),
        ],
        log_prefix.with_suffix(".stdout.log"),
        log_prefix.with_suffix(".stderr.log"),
    )
    return output_tsv


def _extract_numeric(row: dict[str, str], keys: Sequence[str]) -> float:
    for key in keys:
        if key in row and str(row[key]).strip():
            return float(row[key])
    raise MetaQuantError(f"结果表缺少必要列: {', '.join(keys)}")


def load_coverm_table(path: str | Path) -> dict[str, dict[str, float]]:
    tsv_path = _ensure_file(Path(path), "CoverM结果表")
    with tsv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise MetaQuantError(f"CoverM结果表为空: {tsv_path}")
        lower_map = {field.lower(): field for field in reader.fieldnames}
        id_field = lower_map.get("contig") or lower_map.get("contig name") or reader.fieldnames[0]
        count_field = next(
            (field for field in reader.fieldnames if "count" in field.lower()),
            None,
        )
        tpm_field = next(
            (field for field in reader.fieldnames if field.lower() == "tpm" or field.lower().endswith(" tpm")),
            None,
        )
        if count_field is None:
            raise MetaQuantError(f"无法在CoverM结果中识别count列: {tsv_path}")
        if tpm_field is None:
            raise MetaQuantError(f"无法在CoverM结果中识别TPM列: {tsv_path}")
        result: dict[str, dict[str, float]] = {}
        for row in reader:
            gene_id = str(row.get(id_field) or "").strip()
            if not gene_id:
                continue
            result[gene_id] = {
                "count": _extract_numeric(row, [count_field]),
                "tpm": _extract_numeric(row, [tpm_field]),
            }
    if not result:
        raise MetaQuantError(f"CoverM结果表没有可用记录: {tsv_path}")
    return result


def merge_quant_tables(
    gene_metadata_tsv: str | Path,
    dna_coverm_tsv: str | Path,
    rna_coverm_tsv: str | Path,
    out_tsv: str | Path,
) -> Path:
    gene_meta_path = _ensure_file(Path(gene_metadata_tsv), "基因元数据表")
    dna_map = load_coverm_table(dna_coverm_tsv)
    rna_map = load_coverm_table(rna_coverm_tsv)
    out_path = Path(out_tsv).expanduser().resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, str | float | int]] = []
    total_dna_counts = sum(values["count"] for values in dna_map.values())
    if total_dna_counts <= 0:
        raise MetaQuantError("DNA count总和为0，无法计算FPM")

    with gene_meta_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise MetaQuantError(f"基因元数据表为空: {gene_meta_path}")
        for row in reader:
            gene_id = str(row.get("gene_id") or "").strip()
            if not gene_id:
                continue
            dna_values = dna_map.get(gene_id, {"count": 0.0, "tpm": 0.0})
            rna_values = rna_map.get(gene_id, {"count": 0.0, "tpm": 0.0})
            dna_count = float(dna_values["count"])
            rows.append(
                {
                    "gene_id": gene_id,
                    "bin_id": str(row.get("bin_id") or ""),
                    "source_gene_id": str(row.get("source_gene_id") or ""),
                    "nt_length": int(str(row.get("nt_length") or "0") or "0"),
                    "dna_count": round(dna_count, 6),
                    "dna_fpm": round(dna_count * 1_000_000.0 / total_dna_counts, 6),
                    "rna_count": round(float(rna_values["count"]), 6),
                    "rna_tpm": round(float(rna_values["tpm"]), 6),
                }
            )

    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=[
                "gene_id",
                "bin_id",
                "source_gene_id",
                "nt_length",
                "dna_count",
                "dna_fpm",
                "rna_count",
                "rna_tpm",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)
    return out_path


def compute_fpm_tpm_scatter_stats(gene_quant_tsv: str | Path) -> ScatterStats:
    quant_path = _ensure_file(Path(gene_quant_tsv), "gene_quant.tsv")
    xs: list[float] = []
    ys: list[float] = []
    with quant_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise MetaQuantError(f"gene_quant.tsv为空: {quant_path}")
        for row in reader:
            x = float(str(row.get("dna_fpm") or "0") or "0")
            y = float(str(row.get("rna_tpm") or "0") or "0")
            xs.append(math.log10(x + 1.0))
            ys.append(math.log10(y + 1.0))
    if not xs:
        raise MetaQuantError(f"gene_quant.tsv没有可用于绘图的数据: {quant_path}")

    n = len(xs)
    mean_x = sum(xs) / n
    mean_y = sum(ys) / n
    ss_xx = sum((x - mean_x) ** 2 for x in xs)
    ss_yy = sum((y - mean_y) ** 2 for y in ys)
    ss_xy = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))

    if ss_xx == 0:
        slope = 0.0
        intercept = mean_y
        r2 = 0.0
    else:
        slope = ss_xy / ss_xx
        intercept = mean_y - slope * mean_x
        if ss_yy == 0:
            r2 = 1.0
        else:
            corr = ss_xy / math.sqrt(ss_xx * ss_yy) if ss_xx > 0 and ss_yy > 0 else 0.0
            r2 = corr * corr
    return ScatterStats(
        point_count=n,
        r2=round(r2, 6),
        slope=round(slope, 6),
        intercept=round(intercept, 6),
    )


def render_fpm_tpm_scatter(gene_quant_tsv: str | Path, output_png: str | Path, output_stats_tsv: str | Path) -> tuple[Path, Path]:
    quant_path = _ensure_file(Path(gene_quant_tsv), "gene_quant.tsv")
    png_path = Path(output_png).expanduser().resolve()
    stats_path = Path(output_stats_tsv).expanduser().resolve()
    png_path.parent.mkdir(parents=True, exist_ok=True)
    stats_path.parent.mkdir(parents=True, exist_ok=True)

    xs: list[float] = []
    ys: list[float] = []
    with quant_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise MetaQuantError(f"gene_quant.tsv为空: {quant_path}")
        for row in reader:
            x = float(str(row.get("dna_fpm") or "0") or "0")
            y = float(str(row.get("rna_tpm") or "0") or "0")
            xs.append(math.log10(x + 1.0))
            ys.append(math.log10(y + 1.0))

    stats = compute_fpm_tpm_scatter_stats(quant_path)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise MetaQuantError("缺少matplotlib，无法绘制FPM/TPM散点图") from exc

    plt.figure(figsize=(6, 5))
    plt.scatter(xs, ys, s=10, alpha=0.55, edgecolors="none", color="#2f6c8f")
    if xs:
        min_x = min(xs)
        max_x = max(xs)
        line_x = [min_x, max_x]
        line_y = [stats.slope * x + stats.intercept for x in line_x]
        plt.plot(line_x, line_y, color="#c45b4d", linewidth=1.5)
    plt.xlabel("log10(DNA FPM + 1)")
    plt.ylabel("log10(RNA TPM + 1)")
    plt.title(f"FPM vs TPM (R²={stats.r2:.4f}, n={stats.point_count})")
    plt.tight_layout()
    plt.savefig(png_path, dpi=200)
    plt.close()

    with stats_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["point_count", stats.point_count])
        writer.writerow(["r2", f"{stats.r2:.6f}"])
        writer.writerow(["slope", f"{stats.slope:.6f}"])
        writer.writerow(["intercept", f"{stats.intercept:.6f}"])
        writer.writerow(["x_axis", "log10(DNA FPM + 1)"])
        writer.writerow(["y_axis", "log10(RNA TPM + 1)"])
    return png_path, stats_path


def _write_sample_summary(
    sample_out: Path,
    sample: MetaQuantSample,
    merged_tsv: Path,
    scatter_png: Path | None = None,
    scatter_stats_tsv: Path | None = None,
) -> Path:
    summary_path = sample_out / "summary.tsv"
    gene_count = 0
    bin_ids: set[str] = set()
    with merged_tsv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            gene_count += 1
            if row.get("bin_id"):
                bin_ids.add(str(row["bin_id"]))
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "bin_count", "gene_count", "quant_tsv", "scatter_png", "scatter_stats_tsv"])
        writer.writerow(
            [
                sample.sample,
                len(bin_ids),
                gene_count,
                str(merged_tsv),
                str(scatter_png) if scatter_png else "",
                str(scatter_stats_tsv) if scatter_stats_tsv else "",
            ]
        )
    return summary_path


def run_meta_quant_sample(sample: MetaQuantSample, cfg: MetaQuantConfig) -> Path:
    sample_out = cfg.outdir.expanduser().resolve() / sample.sample
    sample_out.mkdir(parents=True, exist_ok=True)

    combined_cds, gene_meta_tsv = _predict_genes_for_bins(sample, cfg, sample_out)

    dna_dir = sample_out / "dna_qc"
    rna_dir = sample_out / "rna_qc"
    dna_dir.mkdir(parents=True, exist_ok=True)
    rna_dir.mkdir(parents=True, exist_ok=True)

    dna_clean1, dna_clean2 = _run_fastp_pair(
        cfg,
        sample.dna_fastq1,
        sample.dna_fastq2,
        dna_dir / f"{sample.sample}.dna.clean.R1.fastq.gz",
        dna_dir / f"{sample.sample}.dna.clean.R2.fastq.gz",
        dna_dir / f"{sample.sample}.dna.fastp.json",
        dna_dir / f"{sample.sample}.dna.fastp.html",
        dna_dir / f"{sample.sample}.dna.fastp",
    )

    rna_fastp1, rna_fastp2 = _run_fastp_pair(
        cfg,
        sample.rna_fastq1,
        sample.rna_fastq2,
        rna_dir / f"{sample.sample}.rna.fastp.R1.fastq.gz",
        rna_dir / f"{sample.sample}.rna.fastp.R2.fastq.gz",
        rna_dir / f"{sample.sample}.rna.fastp.json",
        rna_dir / f"{sample.sample}.rna.fastp.html",
        rna_dir / f"{sample.sample}.rna.fastp",
    )
    if cfg.skip_ribodetector:
        rna_clean1, rna_clean2 = rna_fastp1, rna_fastp2
    else:
        rna_clean1, rna_clean2 = _run_ribodetector_pair(
            cfg,
            rna_fastp1,
            rna_fastp2,
            rna_dir / f"{sample.sample}.rna.non_rrna.R1.fastq",
            rna_dir / f"{sample.sample}.rna.non_rrna.R2.fastq",
            rna_dir / f"{sample.sample}.ribodetector",
        )

    quant_dir = sample_out / "quant"
    quant_dir.mkdir(parents=True, exist_ok=True)
    dna_coverm = _run_coverm_pair(
        cfg,
        combined_cds,
        dna_clean1,
        dna_clean2,
        quant_dir / f"{sample.sample}.dna.coverm.tsv",
        quant_dir / f"{sample.sample}.dna.coverm",
    )
    rna_coverm = _run_coverm_pair(
        cfg,
        combined_cds,
        rna_clean1,
        rna_clean2,
        quant_dir / f"{sample.sample}.rna.coverm.tsv",
        quant_dir / f"{sample.sample}.rna.coverm",
    )

    merged_tsv = merge_quant_tables(
        gene_meta_tsv,
        dna_coverm,
        rna_coverm,
        quant_dir / f"{sample.sample}.gene_quant.tsv",
    )
    scatter_png, scatter_stats_tsv = render_fpm_tpm_scatter(
        merged_tsv,
        quant_dir / f"{sample.sample}.fpm_tpm_scatter.png",
        quant_dir / f"{sample.sample}.fpm_tpm_r2.tsv",
    )
    _write_sample_summary(sample_out, sample, merged_tsv, scatter_png, scatter_stats_tsv)
    return merged_tsv


def run_meta_quant(samples: Sequence[MetaQuantSample], cfg: MetaQuantConfig) -> Path:
    outdir = cfg.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    run_summary = outdir / "run_summary.tsv"
    with run_summary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "status", "output_dir", "gene_quant_tsv"])
        for sample in samples:
            merged_tsv = run_meta_quant_sample(sample, cfg)
            writer.writerow([sample.sample, "done", str(outdir / sample.sample), str(merged_tsv)])
    return run_summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "基于dRep bins、宏基因组双端reads和宏转录组双端reads的基因丰度/表达定量模块"
        )
    )
    parser.add_argument(
        "--manifest",
        required=True,
        help="TSV样本表，包含sample/bins_dir/dna_fastq1/dna_fastq2/rna_fastq1/rna_fastq2",
    )
    parser.add_argument("--outdir", required=True, help="输出目录")
    parser.add_argument("--threads", type=int, default=16, help="线程数，默认16")
    parser.add_argument(
        "--conda-env",
        default="meta_quant",
        help="运行fastp和prodigal的Conda环境名，默认meta_quant；RiboDetector和CoverM使用各自固定环境",
    )
    parser.add_argument(
        "--fastp-min-length",
        type=int,
        default=50,
        help="fastp过滤的最短read长度，默认50",
    )
    parser.add_argument(
        "--coverm-min-covered-fraction",
        type=float,
        default=0.0,
        help="CoverM参数 --min-covered-fraction，默认0",
    )
    parser.add_argument(
        "--coverm-min-read-percent-identity",
        type=float,
        default=0.95,
        help="CoverM参数 --min-read-percent-identity，默认0.95",
    )
    parser.add_argument(
        "--coverm-min-read-aligned-percent",
        type=float,
        default=0.75,
        help="CoverM参数 --min-read-aligned-percent，默认0.75",
    )
    parser.add_argument(
        "--skip-ribodetector",
        action="store_true",
        help="已完成rRNA去除时跳过RiboDetector，直接使用fastp后的RNA reads做定量",
    )
    parser.add_argument("--force", action="store_true", help="覆盖已有结果并重跑")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    samples = load_manifest(args.manifest)
    cfg = MetaQuantConfig(
        outdir=Path(args.outdir),
        threads=args.threads,
        conda_env=args.conda_env,
        fastp_min_length=args.fastp_min_length,
        coverm_min_covered_fraction=args.coverm_min_covered_fraction,
        coverm_min_read_percent_identity=args.coverm_min_read_percent_identity,
        coverm_min_read_aligned_percent=args.coverm_min_read_aligned_percent,
        skip_ribodetector=args.skip_ribodetector,
        force=args.force,
    )
    summary = run_meta_quant(samples, cfg)
    print(f"宏基因组/宏转录组定量完成，结果汇总: {summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
