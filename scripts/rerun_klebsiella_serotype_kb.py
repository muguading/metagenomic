#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
from datetime import datetime
from pathlib import Path


ALIASES = {
    "ST": ("ST", "klebsiella_pneumo_complex__mlst__ST"),
    "virulence_score": ("virulence_score", "klebsiella_pneumo_complex__virulence_score__virulence_score"),
    "resistance_score": ("resistance_score", "klebsiella_pneumo_complex__resistance_score__resistance_score"),
    "Yersiniabactin": ("Yersiniabactin", "klebsiella__ybst__Yersiniabactin"),
    "Colibactin": ("Colibactin", "klebsiella__cbst__Colibactin"),
    "Bla_chr": ("Bla_chr", "klebsiella_pneumo_complex__amr__Bla_chr"),
    "SHV_mutations": ("SHV_mutations", "klebsiella_pneumo_complex__amr__SHV_mutations"),
    "wzi": ("wzi", "klebsiella_pneumo_complex__wzi__wzi"),
    "K_locus": ("K_locus", "klebsiella_pneumo_complex__kaptive__K_locus"),
    "O_locus": ("O_locus", "klebsiella_pneumo_complex__kaptive__O_locus"),
}

LEGACY_COLUMNS = [
    "样本名称",
    "ST",
    "毒力得分",
    "耐药得分",
    "耶尔森菌素",
    "大肠菌素",
    "氨苄类耐药SHV等位基因",
    "SHV耐药突变",
    "wzi荚膜预测",
    "KO血清型",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def value(row: dict[str, str], field: str) -> str:
    for column in ALIASES[field]:
        raw = str(row.get(column) or "").strip()
        if raw and raw.lower() != "nan":
            return raw
    for column, raw_value in row.items():
        if any(column.endswith(f"__{alias}") for alias in ALIASES[field]):
            raw = str(raw_value or "").strip()
            if raw and raw.lower() != "nan":
                return raw
    return "-"


def find_kleborate_table(output_dir: Path) -> Path:
    candidates = [
        output_dir / "klebsiella_pneumo_complex_output.txt",
        output_dir / "klebsiella_output.txt",
        output_dir / "results.txt",
    ]
    candidates.extend(sorted(output_dir.glob("*_output.txt")))
    candidates.extend(sorted(output_dir.glob("*.tsv")))
    candidates.extend(sorted(output_dir.glob("*.txt")))

    seen: set[Path] = set()
    for path in candidates:
        if path in seen or not path.is_file() or path.stat().st_size == 0:
            continue
        seen.add(path)
        try:
            rows = read_tsv(path)
        except Exception:
            continue
        if rows and (value(rows[0], "K_locus") != "-" or value(rows[0], "O_locus") != "-"):
            return path
    raise FileNotFoundError(f"未找到可解析的 kleborate 输出表: {output_dir}")


def legacy_row(raw: dict[str, str], sample: str) -> dict[str, str]:
    k_locus = value(raw, "K_locus")
    o_locus = value(raw, "O_locus")
    return {
        "样本名称": sample,
        "ST": value(raw, "ST"),
        "毒力得分": value(raw, "virulence_score"),
        "耐药得分": value(raw, "resistance_score"),
        "耶尔森菌素": value(raw, "Yersiniabactin"),
        "大肠菌素": value(raw, "Colibactin"),
        "氨苄类耐药SHV等位基因": value(raw, "Bla_chr"),
        "SHV耐药突变": value(raw, "SHV_mutations"),
        "wzi荚膜预测": value(raw, "wzi"),
        "KO血清型": f"{k_locus}|{o_locus}",
    }


def needs_rerun(fasta: Path, force: bool) -> bool:
    if force:
        return True
    sample = fasta.name.removesuffix(".final.fasta")
    result_path = fasta.with_name(f"{sample}_serotype_result.tsv")
    if not result_path.is_file() or result_path.stat().st_size == 0:
        return True
    try:
        rows = read_tsv(result_path)
    except Exception:
        return True
    if not rows or "KO血清型" not in rows[0]:
        return True
    ko = str(rows[0].get("KO血清型") or "").strip()
    return not ko or ko in {"-", "-|-"}


def rerun_one(fasta: Path, output_name: str, backup_suffix: str, force: bool) -> tuple[str, str]:
    sample = fasta.name.removesuffix(".final.fasta")
    if not needs_rerun(fasta, force):
        return sample, "skip:已有 KO血清型"

    sample_dir = fasta.parent
    output_dir = sample_dir / output_name
    if output_dir.exists():
        backup_dir = sample_dir / f"{output_name}.{backup_suffix}.bak"
        if backup_dir.exists():
            shutil.rmtree(backup_dir)
        output_dir.rename(backup_dir)

    cmd = ["kleborate", "-p", "kpsc", "-o", output_name, "-a", fasta.name]
    subprocess.run(cmd, cwd=sample_dir, check=True)

    table_path = find_kleborate_table(output_dir)
    raw_rows = read_tsv(table_path)
    row = legacy_row(raw_rows[0], sample)
    if row["KO血清型"] == "-|-":
        return sample, "fail:未解析到 K_locus/O_locus"

    result_path = sample_dir / f"{sample}_serotype_result.tsv"
    if result_path.exists():
        shutil.copy2(result_path, sample_dir / f"{result_path.name}.{backup_suffix}.bak")
    write_tsv(result_path, [row], LEGACY_COLUMNS)
    write_tsv(sample_dir / f"{sample}.keblo.tsv", [row], LEGACY_COLUMNS)
    return sample, f"ok:{row['KO血清型']}"


def main() -> int:
    parser = argparse.ArgumentParser(description="只补跑克雷伯菌 serotype_kb/klebsiella kleborate 结果。")
    parser.add_argument("--root", default=".", help="fastq_analysis 路径，默认当前目录")
    parser.add_argument("--force", action="store_true", help="即使已有 KO血清型 也重新生成")
    parser.add_argument("--output-dir", default="results", help="kleborate -o 输出目录名，默认 results")
    parser.add_argument("--sample", action="append", default=[], help="只跑指定样本名，可重复传入")
    args = parser.parse_args()

    root = Path(args.root).resolve()
    samples = set(args.sample)
    fastas = sorted(root.rglob("*.final.fasta"))
    if samples:
        fastas = [path for path in fastas if path.name.removesuffix(".final.fasta") in samples]

    backup_suffix = datetime.now().strftime("%Y%m%d_%H%M%S")
    if not fastas:
        print(f"未找到 *.final.fasta: {root}")
        return 1

    failed = 0
    for fasta in fastas:
        try:
            sample, status = rerun_one(fasta, args.output_dir, backup_suffix, args.force)
        except Exception as exc:
            failed += 1
            sample = fasta.name.removesuffix(".final.fasta")
            status = f"fail:{exc}"
        print(f"{sample}\t{status}\t{fasta.parent}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
