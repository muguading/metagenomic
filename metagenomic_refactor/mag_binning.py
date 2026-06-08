from __future__ import annotations

import argparse
import csv
import gzip
import shlex
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

from metagenomic_refactor.common import conda_run_prefix, get_conda_exe


class MagBinningError(RuntimeError):
    """Raised when the MAG binning workflow fails."""


@dataclass(frozen=True)
class MagSample:
    sample: str
    contigs: Path
    fastq1: Path
    fastq2: Path | None = None


@dataclass(frozen=True)
class MagBinningConfig:
    outdir: Path
    threads: int = 16
    min_contig_len: int = 1500
    semibin_env: str = "global"
    semibin_seq_type: str = "short_read"
    vamb_min_fasta: int = 200000
    score_threshold: float = 0.1
    force: bool = False


def _sanitize_sample_name(sample: str) -> str:
    return sample.strip().replace("/", "_").replace(" ", "_")


def _ensure_file(path: Path, label: str) -> Path:
    path = path.expanduser().resolve()
    if not path.is_file():
        raise MagBinningError(f"{label}不存在: {path}")
    return path


def load_manifest(manifest_path: str | Path) -> list[MagSample]:
    manifest = _ensure_file(Path(manifest_path), "样本表")
    samples: list[MagSample] = []
    with manifest.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        expected = {"sample", "contigs", "fastq1"}
        if reader.fieldnames is None or not expected.issubset(set(reader.fieldnames)):
            raise MagBinningError(
                f"样本表缺少必要列，至少需要: {', '.join(sorted(expected))}"
            )

        for row_num, row in enumerate(reader, start=2):
            sample_raw = (row.get("sample") or "").strip()
            contigs_raw = (row.get("contigs") or "").strip()
            fastq1_raw = (row.get("fastq1") or "").strip()
            fastq2_raw = (row.get("fastq2") or "").strip()
            if not sample_raw:
                continue
            if not contigs_raw:
                raise MagBinningError(f"第{row_num}行缺少contigs列")
            if not fastq1_raw:
                raise MagBinningError(f"第{row_num}行缺少fastq1列")

            samples.append(
                MagSample(
                    sample=_sanitize_sample_name(sample_raw),
                    contigs=_ensure_file(Path(contigs_raw), f"{sample_raw}的contigs文件"),
                    fastq1=_ensure_file(Path(fastq1_raw), f"{sample_raw}的fastq1文件"),
                    fastq2=_ensure_file(Path(fastq2_raw), f"{sample_raw}的fastq2文件")
                    if fastq2_raw
                    else None,
                )
            )
    if not samples:
        raise MagBinningError(f"样本表中没有有效样本: {manifest}")
    return samples


def _require_commands(commands: Sequence[str]) -> None:
    missing = [cmd for cmd in commands if shutil.which(cmd) is None]
    if missing:
        raise MagBinningError(f"缺少依赖命令: {', '.join(missing)}")


def _tool_cmd(cmd: Sequence[str]) -> list[str]:
    return conda_run_prefix("mag_aux") + list(cmd)


def _tool_cmd_text(cmd: Sequence[str]) -> str:
    return " ".join(shlex.quote(part) for part in _tool_cmd(cmd))


def _run_command(cmd: Sequence[str], stdout_log: Path, stderr_log: Path) -> None:
    full_cmd = _tool_cmd(cmd)
    stdout_log.parent.mkdir(parents=True, exist_ok=True)
    stderr_log.parent.mkdir(parents=True, exist_ok=True)
    with stdout_log.open("w", encoding="utf-8") as stdout_handle, stderr_log.open(
        "w", encoding="utf-8"
    ) as stderr_handle:
        completed = subprocess.run(full_cmd, stdout=stdout_handle, stderr=stderr_handle, text=True)
    if completed.returncode != 0:
        raise MagBinningError(
            f"命令执行失败: {' '.join(full_cmd)}\n请检查日志: {stdout_log} 和 {stderr_log}"
        )


def _run_shell_command(command: str, stdout_log: Path, stderr_log: Path) -> None:
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
        )
    if completed.returncode != 0:
        raise MagBinningError(
            f"命令执行失败: {command}\n请检查日志: {stdout_log} 和 {stderr_log}"
        )


def _read_fasta_headers(fasta_path: Path) -> Iterable[str]:
    opener = gzip.open if fasta_path.suffix == ".gz" else open
    with opener(fasta_path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                yield line[1:].strip().split()[0]


def _prepare_contigs(sample: MagSample, cfg: MagBinningConfig, sample_out: Path) -> Path:
    contig_dir = sample_out / "filtered_contigs"
    contig_dir.mkdir(parents=True, exist_ok=True)
    filtered_fasta = contig_dir / f"{sample.sample}.min{cfg.min_contig_len}.fa"
    done_flag = contig_dir / ".done"
    contig_path = Path(sample.contigs)

    if cfg.force:
        for stale in contig_dir.glob("*"):
            if stale.is_dir():
                shutil.rmtree(stale)
            else:
                stale.unlink()
        contig_dir.mkdir(parents=True, exist_ok=True)
    elif done_flag.exists() and filtered_fasta.is_file():
        return filtered_fasta

    opener = gzip.open if contig_path.suffix == ".gz" else open
    kept = 0
    with opener(contig_path, "rt", encoding="utf-8", errors="replace") as src, filtered_fasta.open(
        "w", encoding="utf-8"
    ) as dst:
        header = ""
        seq_chunks: list[str] = []
        for line in src:
            if line.startswith(">"):
                if header:
                    sequence = "".join(seq_chunks)
                    if len(sequence) >= cfg.min_contig_len:
                        dst.write(header)
                        dst.write(sequence)
                        dst.write("\n")
                        kept += 1
                header = line if line.endswith("\n") else f"{line}\n"
                seq_chunks = []
                continue
            seq_chunks.append(line.strip())
        if header:
            sequence = "".join(seq_chunks)
            if len(sequence) >= cfg.min_contig_len:
                dst.write(header)
                dst.write(sequence)
                dst.write("\n")
                kept += 1

    if kept == 0:
        raise MagBinningError(
            f"{sample.sample}过滤后没有保留任何长度 >= {cfg.min_contig_len} 的contig: {contig_path}"
        )
    done_flag.touch()
    return filtered_fasta


def bins_dir_to_contigs2bin(bins_dir: str | Path, out_tsv: str | Path) -> Path:
    bin_dir = Path(bins_dir).expanduser().resolve()
    out_path = Path(out_tsv).expanduser().resolve()
    fasta_files: list[Path] = []
    for pattern in ("*.fa", "*.fna", "*.fasta", "*.fa.gz", "*.fna.gz", "*.fasta.gz"):
        fasta_files.extend(sorted(bin_dir.glob(pattern)))
    if not fasta_files:
        raise MagBinningError(f"目录下未找到bin FASTA文件: {bin_dir}")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as handle:
        for fasta in fasta_files:
            bin_name = fasta.name
            for suffix in (".gz", ".fa", ".fna", ".fasta"):
                if bin_name.endswith(suffix):
                    bin_name = bin_name[: -len(suffix)]
            for contig_id in _read_fasta_headers(fasta):
                handle.write(f"{contig_id}\t{bin_name}\n")

    if out_path.stat().st_size == 0:
        raise MagBinningError(f"生成的contigs2bin文件为空: {out_path}")
    return out_path


def vamb_clusters_to_contigs2bin(clusters_tsv: str | Path, out_tsv: str | Path) -> Path:
    clusters_path = _ensure_file(Path(clusters_tsv), "VAMB聚类文件")
    out_path = Path(out_tsv).expanduser().resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    wrote = False
    with clusters_path.open("r", encoding="utf-8", errors="replace") as src, out_path.open(
        "w", encoding="utf-8"
    ) as dst:
        for line_number, line in enumerate(src, start=1):
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                continue
            if line_number == 1 and fields[0].lower() in {"clustername", "cluster", "bin"}:
                continue
            dst.write(f"{fields[1]}\t{fields[0]}\n")
            wrote = True

    if not wrote:
        raise MagBinningError(f"生成的VAMB contigs2bin文件为空: {out_path}")
    return out_path


def _prepare_bam(sample: MagSample, cfg: MagBinningConfig, sample_out: Path) -> tuple[Path, ...]:
    filtered_contigs = _prepare_contigs(sample, cfg, sample_out)
    mapping_dir = sample_out / "mapping"
    mapping_dir.mkdir(parents=True, exist_ok=True)
    done_flag = mapping_dir / ".done"
    assembly_link = mapping_dir / "assembly.fa"
    sorted_bam = mapping_dir / f"{sample.sample}.sorted.bam"

    if cfg.force and mapping_dir.exists():
        for stale in mapping_dir.glob("*"):
            if stale.is_dir():
                shutil.rmtree(stale)
            else:
                stale.unlink()
        mapping_dir.mkdir(parents=True, exist_ok=True)

    if done_flag.exists() and sorted_bam.is_file() and not cfg.force:
        return (sorted_bam,)

    if assembly_link.exists() or assembly_link.is_symlink():
        assembly_link.unlink()
    assembly_link.symlink_to(filtered_contigs)

    quoted_ref = shlex.quote(str(assembly_link))
    quoted_fastq1 = shlex.quote(str(sample.fastq1))
    quoted_bam = shlex.quote(str(sorted_bam))
    minimap2_cmd = _tool_cmd_text(["minimap2", "-ax", "sr", "-t", str(cfg.threads), str(assembly_link)])
    samtools_sort_cmd = _tool_cmd_text(["samtools", "sort", "-@", str(cfg.threads), "-o", str(sorted_bam)])
    if sample.fastq2 is not None:
        quoted_fastq2 = shlex.quote(str(sample.fastq2))
        mapping_cmd = (
            f"{minimap2_cmd} {quoted_fastq1} {quoted_fastq2} "
            f"| {samtools_sort_cmd}"
        )
    else:
        mapping_cmd = (
            f"{minimap2_cmd} {quoted_fastq1} "
            f"| {samtools_sort_cmd}"
        )
    _run_shell_command(
        mapping_cmd,
        mapping_dir / "mapping.stdout.log",
        mapping_dir / "mapping.stderr.log",
    )
    _run_command(
        ["samtools", "index", str(sorted_bam)],
        mapping_dir / "samtools_index.stdout.log",
        mapping_dir / "samtools_index.stderr.log",
    )
    done_flag.touch()
    return (sorted_bam,)


def _run_metabat2(sample: MagSample, cfg: MagBinningConfig, sample_out: Path) -> None:
    filtered_contigs = _prepare_contigs(sample, cfg, sample_out)
    bam_files = _prepare_bam(sample, cfg, sample_out)
    metabat_dir = sample_out / "metabat2"
    done_flag = metabat_dir / ".done"
    depth_tsv = metabat_dir / f"{sample.sample}.depth.tsv"
    bin_prefix = metabat_dir / f"{sample.sample}.bin"
    if done_flag.exists() and not cfg.force:
        return

    metabat_dir.mkdir(parents=True, exist_ok=True)
    cmd_depth = ["jgi_summarize_bam_contig_depths", "--outputDepth", str(depth_tsv)]
    cmd_depth.extend(str(bam) for bam in bam_files)
    _run_command(cmd_depth, metabat_dir / "jgi.stdout.log", metabat_dir / "jgi.stderr.log")

    cmd_metabat = [
        "metabat2",
        "-i",
        str(filtered_contigs),
        "-a",
        str(depth_tsv),
        "-o",
        str(bin_prefix),
        "-m",
        str(cfg.min_contig_len),
        "-t",
        str(cfg.threads),
    ]
    _run_command(
        cmd_metabat,
        metabat_dir / "metabat2.stdout.log",
        metabat_dir / "metabat2.stderr.log",
    )
    done_flag.touch()


def _run_semibin2(sample: MagSample, cfg: MagBinningConfig, sample_out: Path) -> None:
    filtered_contigs = _prepare_contigs(sample, cfg, sample_out)
    bam_files = _prepare_bam(sample, cfg, sample_out)
    semibin_dir = sample_out / "semibin2"
    done_flag = semibin_dir / ".done"
    if done_flag.exists() and not cfg.force:
        return

    semibin_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "SemiBin2",
        "single_easy_bin",
        "--input-fasta",
        str(filtered_contigs),
        "--output",
        str(semibin_dir),
        "--environment",
        cfg.semibin_env,
        "--threads",
        str(cfg.threads),
        "--sequencing-type",
        cfg.semibin_seq_type,
        "--input-bam",
    ]
    cmd.extend(str(bam) for bam in bam_files)
    _run_command(cmd, semibin_dir / "semibin2.stdout.log", semibin_dir / "semibin2.stderr.log")
    done_flag.touch()


def _vamb_cuda_available(sample_out: Path) -> bool:
    cuda_check_dir = sample_out / "vamb_cuda_check"
    cuda_check_dir.mkdir(parents=True, exist_ok=True)
    host_stdout = cuda_check_dir / "host.stdout.log"
    host_stderr = cuda_check_dir / "host.stderr.log"
    env_stdout = cuda_check_dir / "env.stdout.log"
    env_stderr = cuda_check_dir / "env.stderr.log"

    try:
        with host_stdout.open("w", encoding="utf-8") as stdout_handle, host_stderr.open(
            "w", encoding="utf-8"
        ) as stderr_handle:
            host_check = subprocess.run(
                ["nvidia-smi", "-L"],
                stdout=stdout_handle,
                stderr=stderr_handle,
                text=True,
            )
    except OSError as exc:
        host_stderr.write_text(f"{exc}\n", encoding="utf-8")
        return False
    if host_check.returncode != 0:
        return False

    cuda_probe = [
        *conda_run_prefix("mag_aux"),
        "python",
        "-c",
        (
            "import sys; "
            "import torch; "
            "sys.exit(0 if torch.cuda.is_available() else 1)"
        ),
    ]
    try:
        with env_stdout.open("w", encoding="utf-8") as stdout_handle, env_stderr.open(
            "w", encoding="utf-8"
        ) as stderr_handle:
            env_check = subprocess.run(
                cuda_probe,
                stdout=stdout_handle,
                stderr=stderr_handle,
                text=True,
            )
    except OSError as exc:
        env_stderr.write_text(f"{exc}\n", encoding="utf-8")
        return False
    return env_check.returncode == 0


def _run_vamb(sample: MagSample, cfg: MagBinningConfig, sample_out: Path) -> None:
    filtered_contigs = _prepare_contigs(sample, cfg, sample_out)
    vamb_dir = sample_out / "vamb"
    bam_dir = sample_out / "vamb_bams"
    done_flag = sample_out / "vamb.done"
    if done_flag.exists() and not cfg.force:
        return

    bam_files = _prepare_bam(sample, cfg, sample_out)
    if bam_dir.exists():
        shutil.rmtree(bam_dir)
    bam_dir.mkdir(parents=True, exist_ok=True)
    for idx, bam in enumerate(bam_files, start=1):
        target = bam_dir / f"{idx:03d}_{bam.name}"
        target.symlink_to(bam)
    if vamb_dir.exists() and cfg.force:
        shutil.rmtree(vamb_dir)
    elif vamb_dir.exists():
        raise MagBinningError(f"VAMB输出目录已存在，请使用force覆盖: {vamb_dir}")

    cmd = [
        "vamb",
        "bin",
        "default",
        "--outdir",
        str(vamb_dir),
        "--fasta",
        str(filtered_contigs),
        "--bamdir",
        str(bam_dir),
        "--minfasta",
        str(cfg.vamb_min_fasta),
        "-m",
        str(cfg.min_contig_len),
        "-p",
        str(cfg.threads),
        "-o",
    ]
    if _vamb_cuda_available(sample_out):
        cmd.append("--cuda")
    _run_command(cmd, sample_out / "vamb.stdout.log", sample_out / "vamb.stderr.log")
    done_flag.touch()


def _resolve_vamb_clusters(vamb_dir: Path) -> Path:
    split_file = vamb_dir / "vae_clusters_split.tsv"
    unsplit_file = vamb_dir / "vae_clusters_unsplit.tsv"
    if split_file.is_file():
        return split_file
    if unsplit_file.is_file():
        return unsplit_file
    raise MagBinningError(f"未找到VAMB聚类结果文件: {vamb_dir}")


def _run_dastool(sample: MagSample, cfg: MagBinningConfig, sample_out: Path) -> None:
    filtered_contigs = _prepare_contigs(sample, cfg, sample_out)
    dastool_dir = sample_out / "dastool"
    done_flag = dastool_dir / ".done"
    if done_flag.exists() and not cfg.force:
        return

    dastool_dir.mkdir(parents=True, exist_ok=True)
    metabat_tsv = bins_dir_to_contigs2bin(
        sample_out / "metabat2", dastool_dir / "metabat2_contigs2bin.tsv"
    )
    semibin_tsv = bins_dir_to_contigs2bin(
        sample_out / "semibin2" / "output_bins", dastool_dir / "semibin2_contigs2bin.tsv"
    )
    vamb_tsv = vamb_clusters_to_contigs2bin(
        _resolve_vamb_clusters(sample_out / "vamb"), dastool_dir / "vamb_contigs2bin.tsv"
    )

    cmd = [
        "DAS_Tool",
        "-i",
        ",".join((str(metabat_tsv), str(semibin_tsv), str(vamb_tsv))),
        "-l",
        "metabat2,semibin2,vamb",
        "-c",
        str(filtered_contigs),
        "-o",
        str(dastool_dir / sample.sample),
        "--write_bins",
        "--score_threshold",
        str(cfg.score_threshold),
        "--threads",
        str(cfg.threads),
    ]
    _run_command(cmd, dastool_dir / "dastool.stdout.log", dastool_dir / "dastool.stderr.log")
    done_flag.touch()


def _run_drep(sample: MagSample, cfg: MagBinningConfig, sample_out: Path) -> None:
    drep_dir = sample_out / "drep"
    derep_dir = drep_dir / "dereplicated_genomes"
    done_flag = drep_dir / ".done"
    if done_flag.exists() and not cfg.force:
        return

    drep_dir.mkdir(parents=True, exist_ok=True)
    derep_dir.mkdir(parents=True, exist_ok=True)
    dastool_bins_dir = find_dastool_bins_dir(sample_out)
    fasta_paths: list[Path] = []
    for pattern in ("*.fa", "*.fna", "*.fasta"):
        fasta_paths.extend(sorted(dastool_bins_dir.glob(pattern)))
    if not fasta_paths:
        raise MagBinningError(f"DASTool输出目录中没有可用于dRep的bin文件: {dastool_bins_dir}")
    if len(fasta_paths) == 1:
        for stale in derep_dir.glob("*"):
            if stale.is_dir():
                shutil.rmtree(stale)
            else:
                stale.unlink()
        shutil.copy2(fasta_paths[0], derep_dir / fasta_paths[0].name)
        done_flag.touch()
        return

    fasta_inputs = [str(path) for path in fasta_paths]

    cmd = [
        "dRep",
        "dereplicate",
        str(drep_dir),
        "-g",
        *fasta_inputs,
        "-comp",
        "10",
        "-con",
        "5",
        "--S_algorithm",
        "fastANI",
        "-sa",
        "0.95",
        "--ignoreGenomeQuality",
        "-p",
        str(cfg.threads),
    ]
    _run_command(cmd, drep_dir / "drep.stdout.log", drep_dir / "drep.stderr.log")
    done_flag.touch()


def find_dastool_bins_dir(sample_out: str | Path) -> Path:
    dastool_dir = Path(sample_out).expanduser().resolve() / "dastool"
    candidates = [
        dastool_dir / "DASTool_bins",
        dastool_dir / "bins",
    ]
    candidates.extend(sorted(dastool_dir.glob("*_DASTool_bins")))
    for candidate in candidates:
        if candidate.is_dir():
            fasta_files = list(candidate.glob("*.fa")) + list(candidate.glob("*.fasta")) + list(candidate.glob("*.fna"))
            if fasta_files:
                return candidate
    raise MagBinningError(f"未找到DASTool输出bin目录: {dastool_dir}")


def find_final_bins_dir(sample_out: str | Path) -> Path:
    sample_path = Path(sample_out).expanduser().resolve()
    drep_derep_dir = sample_path / "drep" / "dereplicated_genomes"
    if drep_derep_dir.is_dir():
        fasta_files = (
            list(drep_derep_dir.glob("*.fa"))
            + list(drep_derep_dir.glob("*.fna"))
            + list(drep_derep_dir.glob("*.fasta"))
        )
        if fasta_files:
            return drep_derep_dir
    return find_dastool_bins_dir(sample_path)


def export_legacy_binning_layout(sample_out: str | Path, legacy_root: str | Path) -> Path:
    sample_path = Path(sample_out).expanduser().resolve()
    legacy_root_path = Path(legacy_root).expanduser().resolve()
    derep_dir = legacy_root_path / "meta_drep_out" / "dereplicated_genomes"
    derep_dir.mkdir(parents=True, exist_ok=True)

    for stale in derep_dir.glob("*"):
        if stale.is_dir():
            shutil.rmtree(stale)
        else:
            stale.unlink()

    bins_dir = find_final_bins_dir(sample_path)
    copied = 0
    for fasta in sorted(bins_dir.glob("*")):
        if fasta.suffix.lower() not in {".fa", ".fna", ".fasta"}:
            continue
        shutil.copy2(fasta, derep_dir / fasta.name)
        copied += 1

    if copied == 0:
        raise MagBinningError(f"DASTool输出目录中没有可复制的bin文件: {bins_dir}")
    return derep_dir


def run_mag_binning(samples: Sequence[MagSample], config: MagBinningConfig) -> Path:
    _require_commands([get_conda_exe()])
    #-str 不能直接用expanduser
    outdir = Path(config.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    summary_path = outdir / "run_summary.tsv"

    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "status", "output_dir"])
        for sample in samples:
            sample_out = outdir / sample.sample
            sample_out.mkdir(parents=True, exist_ok=True)
            _run_metabat2(sample, config, sample_out)
            _run_semibin2(sample, config, sample_out)
            _run_vamb(sample, config, sample_out)
            _run_dastool(sample, config, sample_out)
            _run_drep(sample, config, sample_out)
            writer.writerow([sample.sample, "done", str(sample_out)])
    return summary_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="MAG分箱、DASTool精修与dRep去重复模块")
    parser.add_argument(
        "--manifest",
        required=True,
        help="TSV样本表，包含sample/contigs/fastq1三列，fastq2可选",
    )
    parser.add_argument("--outdir", required=True, help="输出目录")
    parser.add_argument("--threads", type=int, default=16, help="线程数，默认16")
    parser.add_argument("--min-contig-len", type=int, default=1500, help="最小contig长度，默认1500")
    parser.add_argument("--semibin-env", default="global", help="SemiBin2环境类型，默认global")
    parser.add_argument(
        "--semibin-seq-type",
        default="short_read",
        choices=["short_read", "long_read"],
        help="SemiBin2测序类型，默认short_read",
    )
    parser.add_argument(
        "--vamb-minfasta",
        type=int,
        default=200000,
        help="VAMB导出bin的最小总长度，默认200000",
    )
    parser.add_argument(
        "--score-threshold",
        type=float,
        default=0.1,
        help="DASTool score_threshold，默认0.1",
    )
    parser.add_argument("--force", action="store_true", help="覆盖已有结果并重跑")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    samples = load_manifest(args.manifest)
    config = MagBinningConfig(
        outdir=Path(args.outdir),
        threads=args.threads,
        min_contig_len=args.min_contig_len,
        semibin_env=args.semibin_env,
        semibin_seq_type=args.semibin_seq_type,
        vamb_min_fasta=args.vamb_minfasta,
        score_threshold=args.score_threshold,
        force=args.force,
    )
    summary_path = run_mag_binning(samples, config)
    print(f"MAG分箱流程完成，结果汇总: {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
