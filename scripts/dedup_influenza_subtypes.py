from __future__ import annotations

import argparse
import csv
import re
import shlex
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO


def infer_subtype(header: str, prefix: str) -> str:
    text = str(header or "").strip()
    match = re.search(rf"\b({prefix.upper()}[0-9]+)\b", text, re.IGNORECASE)
    if match:
        return match.group(1).upper()
    for part in text.split("|"):
        part = part.strip()
        if re.fullmatch(rf"{prefix.upper()}[0-9]+", part.upper()):
            return part.upper()
    return "-"


def normalize_record(record, segment_group: str):
    raw_id = str(record.id or "").strip()
    subtype = infer_subtype(raw_id, "H" if segment_group == "HA" else "N")
    accession = re.sub(r"[^A-Za-z0-9._-]+", "_", raw_id.split("|")[0] or raw_id) or "unknown"
    normalized_id = f"A_{segment_group}_{subtype}__{accession}"
    record.id = normalized_id
    record.name = normalized_id
    record.description = normalized_id
    return record, subtype


def run_mmseqs(mmseqs_cmd: str, input_fasta: Path, output_prefix: Path, tmp_dir: Path, min_seq_id: float) -> Path:
    cmd = (
        f"{mmseqs_cmd} easy-linclust "
        f"{shlex.quote(str(input_fasta))} "
        f"{shlex.quote(str(output_prefix))} "
        f"{shlex.quote(str(tmp_dir))} "
        f"--min-seq-id {min_seq_id} -c {min_seq_id} --cov-mode 0"
    )
    completed = subprocess.run(cmd, shell=True)
    if completed.returncode != 0:
        raise RuntimeError(f"mmseqs 去冗余失败: {cmd}")
    rep_fasta = Path(f"{output_prefix}_rep_seq.fasta")
    if not rep_fasta.is_file():
        raise FileNotFoundError(f"未找到 mmseqs 代表序列文件: {rep_fasta}")
    return rep_fasta


def dedup_file(input_fasta: Path, output_fasta: Path, summary_tsv: Path, segment_group: str, mmseqs_cmd: str, min_seq_id: float) -> None:
    grouped = defaultdict(list)
    for record in SeqIO.parse(str(input_fasta), "fasta"):
        normalized, subtype = normalize_record(record, segment_group)
        grouped[subtype].append(normalized)

    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    summary_rows = []
    final_records = []
    with tempfile.TemporaryDirectory(prefix=f"mmseqs_{segment_group.lower()}_") as tmp_root:
        tmp_root_path = Path(tmp_root)
        for subtype, records in sorted(grouped.items()):
            if len(records) == 1:
                final_records.extend(records)
                summary_rows.append({
                    "segment_group": segment_group,
                    "subtype": subtype,
                    "input_count": 1,
                    "output_count": 1,
                    "reduction_pct": 0.0,
                })
                continue
            subtype_dir = tmp_root_path / subtype
            subtype_dir.mkdir(parents=True, exist_ok=True)
            subtype_input = subtype_dir / "input.fa"
            with subtype_input.open("w", encoding="utf-8") as handle:
                SeqIO.write(records, handle, "fasta")
            rep_fasta = run_mmseqs(mmseqs_cmd, subtype_input, subtype_dir / "cluster", subtype_dir / "tmp", min_seq_id)
            reps = list(SeqIO.parse(str(rep_fasta), "fasta"))
            final_records.extend(reps)
            summary_rows.append({
                "segment_group": segment_group,
                "subtype": subtype,
                "input_count": len(records),
                "output_count": len(reps),
                "reduction_pct": round((1 - (len(reps) / len(records))) * 100, 2),
            })

    with output_fasta.open("w", encoding="utf-8") as handle:
        SeqIO.write(final_records, handle, "fasta")
    with summary_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["segment_group", "subtype", "input_count", "output_count", "reduction_pct"], delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)


def main() -> int:
    parser = argparse.ArgumentParser(description="Influenza HA/NA subtype reference deduplication with MMseqs2")
    parser.add_argument("--ha-input", type=Path, required=True)
    parser.add_argument("--na-input", type=Path, required=True)
    parser.add_argument("--ha-output", type=Path, required=True)
    parser.add_argument("--na-output", type=Path, required=True)
    parser.add_argument("--ha-summary", type=Path, required=True)
    parser.add_argument("--na-summary", type=Path, required=True)
    parser.add_argument("--mmseqs-cmd", default="conda run -n ncov --no-capture-output mmseqs")
    parser.add_argument("--min-seq-id", type=float, default=0.98)
    args = parser.parse_args()

    dedup_file(args.ha_input, args.ha_output, args.ha_summary, "HA", args.mmseqs_cmd, args.min_seq_id)
    dedup_file(args.na_input, args.na_output, args.na_summary, "NA", args.mmseqs_cmd, args.min_seq_id)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
