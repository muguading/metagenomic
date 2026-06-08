#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import re
import shutil
import subprocess
import tempfile
import zipfile
from pathlib import Path


SUBTYPES = ("HPIV1", "HPIV2", "HPIV3", "HPIV4a", "HPIV4b")
ACCESSION_RE = re.compile(r"^[A-Z]{1,4}_?\d+(?:\.\d+)?$", re.IGNORECASE)


def parse_args() -> argparse.Namespace:
    project_root = Path(__file__).resolve().parents[1]
    db_root = project_root / "database" / "virus" / "hpiv"
    parser = argparse.ArgumentParser(
        description="使用 datasets_macos 从 HPIV*_db.csv 的 Accession 列提取 accession，批量尝试下载对应 annotation。"
    )
    parser.add_argument("--db-root", type=Path, default=db_root, help="HPIV 数据库目录")
    parser.add_argument(
        "--datasets",
        type=Path,
        default=project_root / "scripts" / "datasets_macos",
        help="datasets_macos 可执行文件路径",
    )
    parser.add_argument(
        "--out-root",
        type=Path,
        default=db_root / "downloaded_annotations",
        help="下载输出目录",
    )
    parser.add_argument(
        "--summary",
        type=Path,
        default=db_root / "downloaded_annotations" / "hpiv_datasets_annotation_summary.tsv",
        help="汇总 TSV 路径",
    )
    parser.add_argument("--subtypes", nargs="*", default=list(SUBTYPES), help="要处理的 subtype")
    parser.add_argument("--chunk-size", type=int, default=100, help="每批 accession 数量")
    parser.add_argument("--limit-per-subtype", type=int, default=0, help="每个 subtype 只处理前 N 条 accession，0 为全部")
    parser.add_argument("--resume", action="store_true", help="已存在结果时跳过")
    return parser.parse_args()


def extract_accessions_from_csv(csv_path: Path, limit: int = 0) -> list[str]:
    accessions: list[str] = []
    seen: set[str] = set()
    with csv_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        reader = csv.DictReader(handle)
        accession_field = None
        fieldnames = [str(name or "").strip() for name in (reader.fieldnames or [])]
        for name in fieldnames:
            if name.lower() == "accession":
                accession_field = name
                break
        if accession_field is None:
            return accessions
        for row in reader:
            token = str((row or {}).get(accession_field) or "").strip()
            if not token:
                continue
            if not ACCESSION_RE.match(token):
                continue
            if token in seen:
                continue
            seen.add(token)
            accessions.append(token)
            if limit > 0 and len(accessions) >= limit:
                break
    return accessions


def build_fasta_accession_version_map(fasta_path: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    if not fasta_path.is_file():
        return mapping
    with fasta_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            token = line[1:].strip().split("|", 1)[0].strip().split()[0].strip()
            if not token or not ACCESSION_RE.match(token):
                continue
            base = token.split(".", 1)[0].strip()
            mapping.setdefault(base, token)
            mapping.setdefault(token, token)
    return mapping


def resolve_versioned_accessions(csv_accessions: list[str], fasta_path: Path) -> list[str]:
    version_map = build_fasta_accession_version_map(fasta_path)
    resolved: list[str] = []
    seen: set[str] = set()
    for accession in csv_accessions:
        token = version_map.get(accession, accession)
        if token in seen:
            continue
        seen.add(token)
        resolved.append(token)
    return resolved


def chunked(items: list[str], size: int) -> list[list[str]]:
    if size <= 0:
        return [items]
    return [items[index:index + size] for index in range(0, len(items), size)]


def run_datasets_download(datasets_bin: Path, accessions: list[str], zip_path: Path) -> tuple[bool, str]:
    with tempfile.NamedTemporaryFile("w", suffix=".txt", encoding="utf-8", delete=False) as handle:
        inputfile = Path(handle.name)
        handle.write("\n".join(accessions) + "\n")
    try:
        command = [
            str(datasets_bin),
            "download",
            "virus",
            "genome",
            "accession",
            "--inputfile",
            str(inputfile),
            "--include",
            "annotation",
            "--filename",
            str(zip_path),
            "--no-progressbar",
        ]
        completed = subprocess.run(command, capture_output=True, text=True, check=False)
        ok = completed.returncode == 0 and zip_path.is_file() and zip_path.stat().st_size > 0
        stderr = (completed.stderr or "").strip()
        stdout = (completed.stdout or "").strip()
        return ok, stderr or stdout
    finally:
        try:
            inputfile.unlink()
        except OSError:
            pass


def collect_annotation_candidates(extract_dir: Path) -> dict[str, list[Path]]:
    mapping: dict[str, list[Path]] = {}
    for path in extract_dir.rglob("*"):
        if not path.is_file():
            continue
        lower_name = path.name.lower()
        if not lower_name.endswith((".gff", ".gff3", ".gtf", ".gbff", ".jsonl", ".json")):
            continue
        for part in path.parts:
            token = part.strip()
            if ACCESSION_RE.match(token):
                mapping.setdefault(token, []).append(path)
        stem_token = path.stem.split("_", 1)[0].strip()
        if ACCESSION_RE.match(stem_token):
            mapping.setdefault(stem_token, []).append(path)
    return mapping


def best_annotation_file(paths: list[Path]) -> Path | None:
    if not paths:
        return None
    ranked = sorted(
        paths,
        key=lambda item: (
            0 if item.suffix.lower() in {".gff3", ".gff"} else 1 if item.suffix.lower() == ".gtf" else 2 if item.suffix.lower() == ".gbff" else 3,
            item.name.lower(),
        ),
    )
    return ranked[0]


def copy_annotation(src: Path, dest_dir: Path, accession: str) -> Path:
    dest_dir.mkdir(parents=True, exist_ok=True)
    suffix = src.suffix if src.suffix else ".txt"
    dest_path = dest_dir / f"{accession}{suffix}"
    shutil.copy2(src, dest_path)
    return dest_path


def main() -> None:
    args = parse_args()
    db_root = args.db_root.expanduser().resolve()
    datasets_bin = args.datasets.expanduser().resolve()
    out_root = args.out_root.expanduser().resolve()
    summary_path = args.summary.expanduser().resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)

    if not datasets_bin.is_file():
        raise FileNotFoundError(f"未找到 datasets 可执行文件: {datasets_bin}")

    summary_rows: list[dict[str, object]] = []

    for subtype in args.subtypes:
        subtype_name = str(subtype).strip()
        csv_path = db_root / f"{subtype_name}_db.csv"
        fasta_path = db_root / f"{subtype_name}_db.fasta"
        subtype_out_dir = out_root / subtype_name
        subtype_out_dir.mkdir(parents=True, exist_ok=True)

        if not csv_path.is_file():
            summary_rows.append(
                {
                    "subtype": subtype_name,
                    "accession": "-",
                    "status": "missing_csv",
                    "annotation_type": "",
                    "annotation_path": "",
                    "note": str(csv_path),
                }
            )
            continue

        csv_accessions = extract_accessions_from_csv(csv_path, args.limit_per_subtype)
        if not csv_accessions:
            summary_rows.append(
                {
                    "subtype": subtype_name,
                    "accession": "-",
                    "status": "no_accession",
                    "annotation_type": "",
                    "annotation_path": "",
                    "note": str(csv_path),
                }
            )
            continue
        accessions = resolve_versioned_accessions(csv_accessions, fasta_path)

        for batch_index, batch in enumerate(chunked(accessions, args.chunk_size), start=1):
            total_batches = math.ceil(len(accessions) / args.chunk_size) if args.chunk_size > 0 else 1
            with tempfile.TemporaryDirectory(prefix=f"{subtype_name}_datasets_") as tmp_dir:
                tmp_root = Path(tmp_dir)
                zip_path = tmp_root / f"{subtype_name}_batch{batch_index}.zip"
                ok, note = run_datasets_download(datasets_bin, batch, zip_path)
                if not ok:
                    for accession in batch:
                        summary_rows.append(
                            {
                                "subtype": subtype_name,
                                "accession": accession,
                                "status": "download_failed",
                                "annotation_type": "",
                                "annotation_path": "",
                                "note": f"batch {batch_index}/{total_batches}: {note}",
                            }
                        )
                    continue

                extract_dir = tmp_root / "unzipped"
                extract_dir.mkdir(parents=True, exist_ok=True)
                try:
                    with zipfile.ZipFile(zip_path, "r") as archive:
                        archive.extractall(extract_dir)
                except zipfile.BadZipFile:
                    for accession in batch:
                        summary_rows.append(
                            {
                                "subtype": subtype_name,
                                "accession": accession,
                                "status": "bad_zip",
                                "annotation_type": "",
                                "annotation_path": "",
                                "note": f"batch {batch_index}/{total_batches}",
                            }
                        )
                    continue

                candidates = collect_annotation_candidates(extract_dir)
                for accession in batch:
                    existing_matches = list(subtype_out_dir.glob(f"{accession}.*"))
                    if args.resume and existing_matches:
                        best_existing = sorted(existing_matches, key=lambda item: item.suffix.lower())[0]
                        summary_rows.append(
                            {
                                "subtype": subtype_name,
                                "accession": accession,
                                "status": "skipped_existing",
                                "annotation_type": best_existing.suffix.lower().lstrip("."),
                                "annotation_path": str(best_existing),
                                "note": "",
                            }
                        )
                        continue

                    matched_files = candidates.get(accession, [])
                    selected = best_annotation_file(matched_files)
                    if selected is None:
                        summary_rows.append(
                            {
                                "subtype": subtype_name,
                                "accession": accession,
                                "status": "no_annotation_file",
                                "annotation_type": "",
                                "annotation_path": "",
                                "note": f"batch {batch_index}/{total_batches}",
                            }
                        )
                        continue
                    copied = copy_annotation(selected, subtype_out_dir, accession)
                    suffix = copied.suffix.lower().lstrip(".")
                    status = "ok" if suffix in {"gff", "gff3"} else "non_gff_annotation"
                    summary_rows.append(
                        {
                            "subtype": subtype_name,
                            "accession": accession,
                            "status": status,
                            "annotation_type": suffix,
                            "annotation_path": str(copied),
                            "note": f"batch {batch_index}/{total_batches}",
                        }
                    )

    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["subtype", "accession", "status", "annotation_type", "annotation_path", "note"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"已完成 HPIV datasets annotation 下载尝试，输出目录: {out_root}")
    print(f"汇总文件: {summary_path}")


if __name__ == "__main__":
    main()
