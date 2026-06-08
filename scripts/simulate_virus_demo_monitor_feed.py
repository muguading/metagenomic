#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import sys
import time
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FASTQ_SUFFIXES = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
DEFAULT_CHUNK_MB = 8.0
DEFAULT_INTERVAL_SECONDS = 8.0

DEMO_PRESETS: dict[str, list[Path]] = {
    "flu": [ROOT / "demo_data" / "flu_demo" / "10239.1.fastq", ROOT / "demo_data" / "flu_demo" / "10239.2.fastq"],
    "rsv": [ROOT / "demo_data" / "rsv_demo" / "10239.1.fastq", ROOT / "demo_data" / "rsv_demo" / "10239.2.fastq"],
    "hmpv": [ROOT / "demo_data" / "hmpv_demo" / "10239.1.fastq", ROOT / "demo_data" / "hmpv_demo" / "10239.2.fastq"],
    "hpiv": [ROOT / "demo_data" / "hpiv_demo" / "10239.1.fastq", ROOT / "demo_data" / "hpiv_demo" / "10239.2.fastq"],
    "hadv": [ROOT / "demo_data" / "hadv_demo" / "10239.1.fastq", ROOT / "demo_data" / "hadv_demo" / "10239.2.fastq"],
    "rhinovirus": [ROOT / "demo_data" / "rhinovirus_demo" / "10239.1.fastq", ROOT / "demo_data" / "rhinovirus_demo" / "10239.2.fastq"],
    "norovirus": [ROOT / "demo_data" / "norovirus_demo" / "10239.1.fastq", ROOT / "demo_data" / "norovirus_demo" / "10239.2.fastq"],
    "rotavirus": [ROOT / "demo_data" / "rotavirus_a_demo" / "10239.1.fastq", ROOT / "demo_data" / "rotavirus_a_demo" / "10239.2.fastq"],
    "denv": [ROOT / "demo_data" / "denv_demo" / "10239.1.fastq", ROOT / "demo_data" / "denv_demo" / "10239.2.fastq"],
    "zikav": [ROOT / "demo_data" / "zikav_demo_1" / "10239.1.fastq", ROOT / "demo_data" / "zikav_demo_1" / "10239.2.fastq"],
    "astro": [ROOT / "demo_data" / "astro_demo_paired" / "10239.1.fastq", ROOT / "demo_data" / "astro_demo_paired" / "10239.2.fastq"],
    "mpox": [ROOT / "demo_data" / "hmpxv_demo" / "hmpxv_demo.raw.fastq"],
    "ncov": [ROOT / "demo_data" / "ncov" / "All.fastq"],
}


@dataclass
class StreamJob:
    source: Path
    target: Path
    total_bytes: int
    reader: object
    writer: object
    written_bytes: int = 0
    finished: bool = False


def format_bytes(value: int) -> str:
    size = float(value)
    units = ["B", "KB", "MB", "GB", "TB"]
    for unit in units:
        if size < 1024 or unit == units[-1]:
            return f"{size:.1f}{unit}" if unit != "B" else f"{int(size)}B"
        size /= 1024
    return f"{value}B"


def log(message: str) -> None:
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}", flush=True)


def list_demos() -> int:
    log("可用病毒 demo 预设：")
    for key in sorted(DEMO_PRESETS):
        rels = ", ".join(str(path.relative_to(ROOT)) for path in DEMO_PRESETS[key])
        print(f"  - {key}: {rels}")
    return 0


def is_fastq_path(path: Path) -> bool:
    lower = path.name.lower()
    return any(lower.endswith(suffix) for suffix in FASTQ_SUFFIXES)


def collect_source_files(demos: list[str], sources: list[str], recursive: bool) -> list[Path]:
    files: list[Path] = []
    seen: set[str] = set()
    for demo in demos:
        if demo not in DEMO_PRESETS:
            raise SystemExit(f"未识别的 demo 预设：{demo}")
        for path in DEMO_PRESETS[demo]:
            if not path.is_file():
                raise SystemExit(f"demo 文件不存在：{path}")
            resolved = str(path.resolve())
            if resolved not in seen:
                files.append(path.resolve())
                seen.add(resolved)
    for raw in sources:
        path = Path(raw).expanduser().resolve()
        if path.is_file():
            if not is_fastq_path(path):
                raise SystemExit(f"不是 fastq/fq 文件：{path}")
            resolved = str(path)
            if resolved not in seen:
                files.append(path)
                seen.add(resolved)
            continue
        if not path.is_dir():
            raise SystemExit(f"路径不存在：{path}")
        candidates = path.rglob("*") if recursive else path.glob("*")
        matched = sorted(item.resolve() for item in candidates if item.is_file() and is_fastq_path(item))
        if not matched:
            raise SystemExit(f"目录内未找到 fastq/fq 文件：{path}")
        for item in matched:
            resolved = str(item)
            if resolved not in seen:
                files.append(item)
                seen.add(resolved)
    if not files:
        raise SystemExit("请至少通过 --demo 或 --source 提供一组病毒 fastq 文件。")
    return files


def ensure_targets(files: list[Path], dest_dir: Path, overwrite: bool) -> list[Path]:
    targets: list[Path] = []
    used_names: set[str] = set()
    for source in files:
        name = source.name
        if name in used_names:
            raise SystemExit(f"目标目录内会出现重名 fastq，当前脚本不自动改名：{name}")
        used_names.add(name)
        target = dest_dir / name
        if target.exists() and not overwrite:
            raise SystemExit(f"目标文件已存在，请先清理或加 --overwrite：{target}")
        targets.append(target)
    return targets


def open_jobs(files: list[Path], targets: list[Path], overwrite: bool) -> list[StreamJob]:
    jobs: list[StreamJob] = []
    for source, target in zip(files, targets):
        if overwrite and target.exists():
            target.unlink()
        target.parent.mkdir(parents=True, exist_ok=True)
        reader = source.open("rb")
        writer = target.open("wb")
        jobs.append(
            StreamJob(
                source=source,
                target=target,
                total_bytes=source.stat().st_size,
                reader=reader,
                writer=writer,
            )
        )
    return jobs


def close_jobs(jobs: list[StreamJob]) -> None:
    for job in jobs:
        try:
            job.reader.close()
        except Exception:
            pass
        try:
            job.writer.close()
        except Exception:
            pass


def stream_jobs(jobs: list[StreamJob], chunk_bytes: int, interval_seconds: float) -> None:
    active_jobs = [job for job in jobs if not job.finished]
    round_index = 0
    while active_jobs:
        round_index += 1
        changed = False
        for job in list(active_jobs):
            data = job.reader.read(chunk_bytes)
            if not data:
                job.writer.flush()
                os.fsync(job.writer.fileno())
                job.finished = True
                active_jobs.remove(job)
                log(f"已完成 {job.target.name} · {format_bytes(job.total_bytes)}")
                continue
            job.writer.write(data)
            job.writer.flush()
            os.fsync(job.writer.fileno())
            job.written_bytes += len(data)
            changed = True
            percent = (job.written_bytes / job.total_bytes) * 100 if job.total_bytes else 100.0
            log(
                f"第 {round_index} 轮写入 {job.target.name} · "
                f"{format_bytes(job.written_bytes)}/{format_bytes(job.total_bytes)} "
                f"({percent:.1f}%)"
            )
        if active_jobs and changed:
            time.sleep(interval_seconds)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="利用病毒 demo fastq 缓慢生成测序文件，模拟监听目录持续出数状态。",
    )
    parser.add_argument("--demo", action="append", default=[], help="使用内置病毒 demo 预设，可重复传入，如 flu、norovirus、zikav。")
    parser.add_argument("--source", action="append", default=[], help="额外指定 fastq 文件或目录，可重复传入。")
    parser.add_argument("--dest", required=False, help="模拟输出目录，也就是监听任务要盯住的目录。")
    parser.add_argument("--chunk-mb", type=float, default=DEFAULT_CHUNK_MB, help=f"每轮给每个 fastq 追加多少 MB，默认 {DEFAULT_CHUNK_MB}。")
    parser.add_argument("--interval-seconds", type=float, default=DEFAULT_INTERVAL_SECONDS, help=f"每轮写入后的等待秒数，默认 {DEFAULT_INTERVAL_SECONDS}。")
    parser.add_argument("--start-delay", type=float, default=0.0, help="正式开始写入前先等待多少秒，方便先启动监听。")
    parser.add_argument("--overwrite", action="store_true", help="如果目标目录已有同名文件，先覆盖。")
    parser.add_argument("--recursive", action="store_true", help="当 --source 给的是目录时，递归搜索 fastq 文件。")
    parser.add_argument("--list-demos", action="store_true", help="列出可用的 demo 预设并退出。")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.list_demos:
        return list_demos()

    if not args.dest:
        parser.error("未提供 --dest。")
    dest_dir = Path(args.dest).expanduser().resolve()
    dest_dir.mkdir(parents=True, exist_ok=True)

    if args.chunk_mb <= 0:
        raise SystemExit("--chunk-mb 必须大于 0。")
    if args.interval_seconds < 0:
        raise SystemExit("--interval-seconds 不能小于 0。")
    if args.start_delay < 0:
        raise SystemExit("--start-delay 不能小于 0。")

    files = collect_source_files(args.demo, args.source, args.recursive)
    targets = ensure_targets(files, dest_dir, args.overwrite)
    chunk_bytes = max(1, int(args.chunk_mb * 1024 * 1024))

    log(f"目标目录：{dest_dir}")
    log(f"本次共准备 {len(files)} 个 fastq 文件，按每轮每文件 {format_bytes(chunk_bytes)} 的速度写入。")
    for source, target in zip(files, targets):
        log(f"  {source.relative_to(ROOT) if source.is_relative_to(ROOT) else source} -> {target}")

    if args.start_delay > 0:
        log(f"等待 {args.start_delay:.1f} 秒后开始写入，方便你先启动监听任务。")
        time.sleep(args.start_delay)

    jobs = open_jobs(files, targets, args.overwrite)
    try:
        stream_jobs(jobs, chunk_bytes=chunk_bytes, interval_seconds=args.interval_seconds)
    finally:
        close_jobs(jobs)

    log("全部 fastq 已写入完成。现在只要目录继续保持静默，就会进入监听稳定窗口。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
