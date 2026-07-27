from __future__ import annotations

import csv
from datetime import datetime
from pathlib import Path


BATCH_INPUT_HEADERS_ZH = ["样本名称", "三代数据", "二代数据左", "二代数据右", "物种信息"]


def write_batch_input(project_root: Path, rows: list[list[str]]) -> Path:
    batch_dir = project_root / "generated_batch_inputs"
    batch_dir.mkdir(parents=True, exist_ok=True)
    target = batch_dir / f"batch_input_{datetime.now().strftime('%Y%m%d%H%M%S%f')}.tsv"
    with target.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(BATCH_INPUT_HEADERS_ZH)
        writer.writerows(rows)
    return target.resolve()
