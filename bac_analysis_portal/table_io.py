from __future__ import annotations

import csv
from pathlib import Path

def _read_tsv_rows(path: Path) -> dict:
    if not path.is_file():
        return {"columns": [], "rows": []}
    try:
        with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
            reader = csv.reader(handle, delimiter="	")
            rows = list(reader)
    except OSError:
        return {"columns": [], "rows": []}
    if not rows:
        return {"columns": [], "rows": []}
    header = rows[0]
    data_rows = rows[1:]
    return {"columns": header, "rows": data_rows}


def write_tsv(path: Path, columns: list[str], rows: list[list[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(columns)
        for row in rows:
            writer.writerow(list(row))
