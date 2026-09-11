from __future__ import annotations

import re
from pathlib import Path

from .task_manager import SHORT_READ_R1_RE, SHORT_READ_R2_RE, ValidationError


def scan_fastq_directory_for_batch_rows(input_dir: Path, species: str) -> list[list[str]]:
    if not input_dir.is_dir():
        raise ValidationError(f"二代测序目录不存在: {input_dir}")
    groups: dict[str, dict[str, str]] = {}
    for path in sorted(input_dir.rglob("*")):
        if not path.is_file() or not path.name.lower().endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
            continue
        sample_name, field = classify_fastq_for_batch(path)
        groups.setdefault(sample_name, {"sample_name": sample_name, "third_gen": "", "short_left": "", "short_right": ""})[field] = str(path.resolve())
    return [[name, record["third_gen"], record["short_left"], record["short_right"], species] for name, record in sorted(groups.items()) if any([record["third_gen"], record["short_left"], record["short_right"]])]


def classify_fastq_for_batch(path: Path) -> tuple[str, str]:
    name = path.name
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
            break
    field, sample_name = "third_gen", name
    if SHORT_READ_R1_RE.search(name):
        field, sample_name = "short_left", SHORT_READ_R1_RE.sub("_", name)
    elif SHORT_READ_R2_RE.search(name):
        field, sample_name = "short_right", SHORT_READ_R2_RE.sub("_", name)
    return re.sub(r"[._-]+$", "", sample_name).strip() or path.stem, field
