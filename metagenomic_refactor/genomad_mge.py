from __future__ import annotations

import argparse
import csv
import os
import re
import shlex
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from metagenomic_refactor.common import conda_run_prefix


class GeNomadMGEError(RuntimeError):
    """Raised when the geNomad MGE workflow fails."""


@dataclass(frozen=True)
class GeNomadSample:
    sample: str
    fasta: Path


@dataclass(frozen=True)
class GeNomadConfig:
    outdir: Path
    database: Path
    threads: int = 16
    conda_env: str = "genomad_aux"
    force: bool = False


@dataclass(frozen=True)
class IntervalFeature:
    sequence_id: str
    start: int
    end: int
    label: str
    source: str
    feature_type: str
    metadata: str = ""


GENOMAD_DB_DEFAULT = "/data/deploy/meta_genome/database/genomad/genomad_db/"
MOBILEOG_DB_DEFAULT = "/data1/shanghai_pip/meta_genome/database/beatrix/mobileOG-db"
MOBILEOG_META_DEFAULT = "/data1/shanghai_pip/meta_genome/database/beatrix/mobileOG-db-beatrix-1.6-All.csv"

MGE_BOUNDARY_KEYWORDS = ("plasmid", "provirus", "phage", "virus")
CORE_IE_KEYWORDS = (
    "integrase",
    "excisionase",
    "recombinase",
    "transposase",
    "resolvase",
    "insertion sequence",
)
CORE_TRANSFER_KEYWORDS = (
    "conjug",
    "relaxase",
    "mobiliz",
    "type iv",
    "t4ss",
    "virb",
    "tra",
    "trb",
    "mob",
    "oriT".lower(),
)
CORE_PHAGE_KEYWORDS = (
    "capsid",
    "tail",
    "terminase",
    "portal",
    "baseplate",
    "head protein",
    "holin",
    "lysin",
    "tape measure",
    "phage",
)
RRR_STD_KEYWORDS = (
    "replication",
    "replicase",
    "partition",
    "toxin-antitoxin",
    "toxin",
    "antitoxin",
    "stability",
    "resistance",
)


def _sanitize_sample_name(sample: str) -> str:
    return sample.strip().replace("/", "_").replace(" ", "_")


def _ensure_file(path: Path, label: str) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        raise GeNomadMGEError(f"{label}不存在: {resolved}")
    return resolved


def _ensure_dir(path: Path, label: str) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_dir():
        raise GeNomadMGEError(f"{label}不存在: {resolved}")
    return resolved


def _runtime_method() -> str:
    try:
        from metagenomic_refactor.context import get_runtime_context
        return str(get_runtime_context().method or "").strip()
    except Exception:
        return ""


def _is_meta_method(method: str | None = None) -> bool:
    return str(_runtime_method() if method is None else method).strip() == "meta"


def _safe_int(value: str | int | None) -> int | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        return int(float(text))
    except ValueError:
        return None


def _interval_distance(start1: int, end1: int, start2: int, end2: int) -> int:
    if end1 < start2:
        return start2 - end1
    if end2 < start1:
        return start1 - end2
    return 0


def _parse_gff_attributes(raw: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in raw.split(";"):
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        attrs[key.strip()] = value.strip()
    return attrs


def load_manifest(manifest_path: str | Path) -> list[GeNomadSample]:
    manifest = _ensure_file(Path(manifest_path), "样本表")
    samples: list[GeNomadSample] = []
    with manifest.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample", "fasta"}
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            raise GeNomadMGEError(f"样本表缺少必要列，至少需要: {', '.join(sorted(required))}")
        for row_num, row in enumerate(reader, start=2):
            sample_raw = str(row.get("sample") or "").strip()
            fasta_raw = str(row.get("fasta") or "").strip()
            if not sample_raw:
                continue
            if not fasta_raw:
                raise GeNomadMGEError(f"第{row_num}行缺少fasta列")
            samples.append(
                GeNomadSample(
                    sample=_sanitize_sample_name(sample_raw),
                    fasta=_ensure_file(Path(fasta_raw), f"{sample_raw}的fasta文件"),
                )
            )
    if not samples:
        raise GeNomadMGEError(f"样本表中没有有效样本: {manifest}")
    return samples


def load_samples(
    *,
    manifest: str | None,
    sample: str | None,
    fasta: str | None,
) -> list[GeNomadSample]:
    if manifest:
        return load_manifest(manifest)
    if sample and fasta:
        return [
            GeNomadSample(
                sample=_sanitize_sample_name(sample),
                fasta=_ensure_file(Path(fasta), f"{sample}的fasta文件"),
            )
        ]
    raise GeNomadMGEError("请提供 --manifest，或同时提供 --sample 和 --fasta")


def _genomad_prefix() -> list[str]:
    return conda_run_prefix("genomad_aux")


def _tool_prefix(cfg: GeNomadConfig) -> list[str]:
    return conda_run_prefix(cfg.conda_env)


def _run_command(
    cmd: Sequence[str],
    stdout_log: Path,
    stderr_log: Path,
    cwd: Path | None = None,
) -> None:
    full_cmd = _genomad_prefix() + list(cmd)
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
        raise GeNomadMGEError(
            f"命令执行失败: {' '.join(full_cmd)}\n请检查日志: {stdout_log} 和 {stderr_log}"
        )


def _safe_clear_dir(path: Path) -> None:
    if not path.exists():
        return
    for item in path.iterdir():
        if item.is_dir() and not item.is_symlink():
            shutil.rmtree(item)
        else:
            item.unlink()


def _run_shell_command(command: str, stdout_log: Path, stderr_log: Path, cwd: Path | None = None) -> None:
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
        raise GeNomadMGEError(
            f"命令执行失败: {command}\n请检查日志: {stdout_log} 和 {stderr_log}"
        )


def _run_shell_command_with_cfg(
    cfg: GeNomadConfig,
    command: str,
    stdout_log: Path,
    stderr_log: Path,
    cwd: Path | None = None,
) -> None:
    prefix = shlex.join(_tool_prefix(cfg))
    wrapped = f"{prefix} bash -lc {shlex.quote(command)}"
    _run_shell_command(wrapped, stdout_log, stderr_log, cwd=cwd)


def _find_summary_table(sample_dir: Path, category: str) -> Path | None:
    candidates = [
        path
        for path in sample_dir.rglob("*.tsv")
        if path.is_file() and category in path.name.lower() and "summary" in path.name.lower()
    ]
    summary_only = [
        path for path in candidates
        if "genes" not in path.name.lower() and "gene" not in path.name.lower()
    ]
    filtered = summary_only or candidates
    if not filtered:
        return None
    filtered.sort(key=lambda path: (len(path.parts), len(path.name), str(path)))
    return filtered[0]


def _find_genes_table(sample_dir: Path, category: str) -> Path | None:
    candidates = [
        path
        for path in sample_dir.rglob("*.tsv")
        if path.is_file() and category in path.name.lower() and "genes" in path.name.lower()
    ]
    if not candidates:
        return None
    candidates.sort(key=lambda path: (len(path.parts), len(path.name), str(path)))
    return candidates[0]


def _extract_sequence_id(row: dict[str, str]) -> str:
    for key in ("seq_name", "sequence_name", "source_seq", "contig", "contig_name", "name", "virus_name", "plasmid_name"):
        value = str(row.get(key) or "").strip()
        if value:
            return value
    return ""


def _extract_length(row: dict[str, str]) -> str:
    for key in ("length", "seq_length", "size"):
        value = str(row.get(key) or "").strip()
        if value:
            return value
    return ""


def _extract_score(row: dict[str, str]) -> str:
    for key in ("score", "plasmid_score", "virus_score", "provirus_score", "fdr"):
        value = str(row.get(key) or "").strip()
        if value:
            return value
    return ""


def _extract_taxonomy(row: dict[str, str]) -> str:
    for key in ("taxonomy", "lineage", "taxname", "taxon"):
        value = str(row.get(key) or "").strip()
        if value:
            return value
    return ""


def _extract_coordinates(row: dict[str, str]) -> tuple[str, str]:
    for key in ("coordinates", "coord", "location"):
        value = str(row.get(key) or "").strip()
        if not value:
            continue
        match = re.search(r"(\d+)\D+(\d+)", value)
        if match:
            return match.group(1), match.group(2)
    start = str(row.get("start") or "").strip()
    end = str(row.get("end") or "").strip()
    if start and end:
        return start, end
    length = _extract_length(row)
    if length:
        return "1", length
    return start, end


def _load_gene_coordinate_map(path: Path | None) -> dict[str, tuple[str, str]]:
    if path is None or not path.is_file():
        return {}
    coords: dict[str, tuple[int, int]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return {}
        for row in reader:
            gene_id = str(row.get("gene") or "").strip()
            if not gene_id or "_" not in gene_id:
                continue
            seq_id = gene_id.rsplit("_", 1)[0]
            start_text = str(row.get("start") or "").strip()
            end_text = str(row.get("end") or "").strip()
            if not start_text or not end_text:
                continue
            try:
                start = int(start_text)
                end = int(end_text)
            except ValueError:
                continue
            current = coords.get(seq_id)
            if current is None:
                coords[seq_id] = (start, end)
            else:
                coords[seq_id] = (min(current[0], start), max(current[1], end))
    return {key: (str(value[0]), str(value[1])) for key, value in coords.items()}


def _load_summary_rows(
    path: Path,
    sample: str,
    mge_type: str,
    gene_coordinate_map: dict[str, tuple[str, str]] | None = None,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    gene_coordinate_map = gene_coordinate_map or {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return rows
        for row in reader:
            sequence_id = _extract_sequence_id(row)
            start, end = gene_coordinate_map.get(sequence_id, _extract_coordinates(row))
            rows.append(
                {
                    "sample": sample,
                    "mge_type": mge_type,
                    "sequence_id": sequence_id,
                    "length": _extract_length(row),
                    "start": start,
                    "end": end,
                    "score": _extract_score(row),
                    "topology": str(row.get("topology") or ""),
                    "taxonomy": _extract_taxonomy(row),
                    "source_tsv": str(path),
                }
            )
    return rows


def _load_mobileog_metadata(meta_csv: Path) -> dict[str, dict[str, str]]:
    metadata: dict[str, dict[str, str]] = {}
    with meta_csv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            key = str(row.get("mobileOG Entry Name") or "").strip()
            if key:
                metadata[key] = row
    return metadata


def _load_mobileog_rows(path: Path, sample: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    if not path.is_file():
        return rows
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows.append(
                {
                    "sample": sample,
                    "method": "mobileOG",
                    "mge_type": str(row.get("类型") or ""),
                    "sequence_id": str(row.get("序列名称") or ""),
                    "start": str(row.get("序列起始") or ""),
                    "end": str(row.get("序列终止") or ""),
                    "length": str(row.get("长度") or ""),
                    "score": str(row.get("比对得分") or ""),
                    "annotation": str(row.get("元件名称") or row.get("注释") or ""),
                    "source_tsv": str(path),
                }
            )
    return rows


def _load_genomad_combined_rows(path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    if not path.is_file():
        return rows
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows.append(
                {
                    "sample": str(row.get("sample") or ""),
                    "method": "geNomad",
                    "mge_type": str(row.get("mge_type") or ""),
                    "sequence_id": str(row.get("sequence_id") or ""),
                    "start": str(row.get("start") or ""),
                    "end": str(row.get("end") or ""),
                    "length": str(row.get("length") or ""),
                    "score": str(row.get("score") or ""),
                    "annotation": str(row.get("taxonomy") or ""),
                    "source_tsv": str(row.get("source_tsv") or path),
                }
            )
    return rows


def aggregate_genomad_results(sample_dir: str | Path, sample: str, outdir: str | Path) -> Path:
    sample_path = _ensure_dir(Path(sample_dir), "geNomad样本输出目录")
    out_path = Path(outdir).expanduser().resolve()
    out_path.mkdir(parents=True, exist_ok=True)

    plasmid_tsv = _find_summary_table(sample_path, "plasmid")
    provirus_tsv = _find_summary_table(sample_path, "provirus") or _find_summary_table(sample_path, "virus")
    plasmid_genes_tsv = _find_genes_table(sample_path, "plasmid")
    provirus_genes_tsv = _find_genes_table(sample_path, "provirus") or _find_genes_table(sample_path, "virus")

    plasmid_rows = (
        _load_summary_rows(plasmid_tsv, sample, "plasmid", _load_gene_coordinate_map(plasmid_genes_tsv))
        if plasmid_tsv else []
    )
    provirus_rows = (
        _load_summary_rows(provirus_tsv, sample, "provirus", _load_gene_coordinate_map(provirus_genes_tsv))
        if provirus_tsv else []
    )
    merged_rows = plasmid_rows + provirus_rows

    merged_tsv = out_path / f"{sample}.genomad_mge_summary.tsv"
    with merged_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=[
                "sample",
                "mge_type",
                "sequence_id",
                "length",
                "start",
                "end",
                "score",
                "topology",
                "taxonomy",
                "source_tsv",
            ],
        )
        writer.writeheader()
        writer.writerows(merged_rows)

    _write_category_table(out_path / f"{sample}.genomad_plasmid_summary.tsv", plasmid_rows)
    _write_category_table(out_path / f"{sample}.genomad_provirus_summary.tsv", provirus_rows)
    return merged_tsv


def _write_category_table(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=[
                "sample",
                "mge_type",
                "sequence_id",
                "length",
                "start",
                "end",
                "score",
                "topology",
                "taxonomy",
                "source_tsv",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


def _write_integrated_summary(path: Path, rows: list[dict[str, str]]) -> Path:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=[
                "sample",
                "method",
                "mge_type",
                "sequence_id",
                "start",
                "end",
                "length",
                "score",
                "annotation",
                "source_tsv",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)
    return path


def integrate_mge_tables(outdir: str | Path, sample: str) -> Path:
    out_path = Path(outdir).expanduser().resolve()
    genomad_summary = out_path / f"{sample}.genomad_mge_summary.tsv"
    genomad_dir = next(
        (
            path
            for path in sorted(out_path.rglob("genomad"))
            if path.is_dir() and path.parent.name == sample
        ),
        None,
    )
    if genomad_dir is not None:
        genomad_summary = aggregate_genomad_results(genomad_dir, sample, genomad_dir.parent)
    elif not genomad_summary.is_file():
        nested = next(
            (
                path
                for path in sorted(out_path.rglob(f"{sample}.genomad_mge_summary.tsv"))
                if path.is_file()
            ),
            None,
        )
        if nested is not None:
            genomad_summary = nested

    mobileog_summary = out_path / f"{sample}.mobileog.tsv"
    if not mobileog_summary.is_file():
        nested = next(
            (
                path
                for path in sorted(out_path.rglob(f"{sample}.mobileog.tsv"))
                if path.is_file()
            ),
            None,
        )
        if nested is not None:
            mobileog_summary = nested

    rows = _load_genomad_combined_rows(genomad_summary) + _load_mobileog_rows(mobileog_summary, sample)
    return _write_integrated_summary(out_path / f"{sample}.integrated_mge_summary.tsv", rows)


def _resolve_annotation_sources(pre: str, cwd: Path, cfg: GeNomadConfig) -> tuple[Path, Path, Path, Path]:
    if _is_meta_method():
        gff_path = _ensure_meta_orf_gff(pre, cwd, cfg)
        return (
            cwd / f"{pre}.integrated_mge_summary.tsv",
            cwd / "bin_card.tsv",
            cwd / "bin_vfdb.tsv",
            gff_path,
        )
    return (
        cwd / f"{pre}.integrated_mge_summary.tsv",
        cwd / f"{pre}.card.tsv",
        cwd / f"{pre}.vfdb.tsv",
        cwd / f"{pre}_prokka" / f"{pre}.gff",
    )


def _ensure_meta_orf_gff(pre: str, cwd: Path, cfg: GeNomadConfig) -> Path:
    candidates = [
        cwd / "tmp_combine.genes.gff",
        cwd / "tmp_combine.gff",
        cwd / f"{pre}.genes.gff",
        cwd / f"{pre}.gff",
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate

    fasta_path = cwd / "tmp_combine.fa"
    if not fasta_path.is_file():
        raise GeNomadMGEError(f"meta 方法未找到 ORF 预测所需的 fasta 文件: {fasta_path}")
    gff_path = cwd / "tmp_combine.genes.gff"
    _run_shell_command_with_cfg(
        cfg,
        f"prodigal -i {shlex.quote(str(fasta_path))} -o {shlex.quote(str(gff_path))} -f gff -p meta",
        cwd / "prodigal.stdout.log",
        cwd / "prodigal.stderr.log",
        cwd=cwd,
    )
    return gff_path


def _load_interval_features(path: Path, sample: str) -> list[IntervalFeature]:
    if not path.is_file():
        return []
    rows: list[IntervalFeature] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            sequence_id = str(row.get("sequence_id") or "").strip()
            start = _safe_int(row.get("start"))
            end = _safe_int(row.get("end"))
            if not sequence_id or start is None or end is None:
                continue
            annotation = str(row.get("annotation") or "").strip()
            rows.append(
                IntervalFeature(
                    sequence_id=sequence_id,
                    start=min(start, end),
                    end=max(start, end),
                    label=str(row.get("mge_type") or annotation or "MGE").strip(),
                    source=str(row.get("method") or "MGE").strip(),
                    feature_type=_classify_mge_feature_type(row),
                    metadata=annotation or sample,
                )
            )
    return rows


def _classify_mge_feature_type(row: dict[str, str]) -> str:
    text = " ".join(
        [
            str(row.get("mge_type") or ""),
            str(row.get("annotation") or ""),
            str(row.get("method") or ""),
        ]
    ).lower()
    if any(keyword in text for keyword in MGE_BOUNDARY_KEYWORDS):
        return "boundary"
    if any(keyword in text for keyword in CORE_IE_KEYWORDS):
        return "ie"
    if any(keyword in text for keyword in CORE_TRANSFER_KEYWORDS):
        return "transfer"
    if any(keyword in text for keyword in CORE_PHAGE_KEYWORDS):
        return "phage"
    if any(keyword in text for keyword in RRR_STD_KEYWORDS):
        return "rrr_std"
    return "mge"


def _load_gene_hits(path: Path, hit_type: str, is_meta: bool) -> list[IntervalFeature]:
    if not path.is_file():
        return []
    rows: list[IntervalFeature] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            feature = _parse_gene_hit_row(row, hit_type, is_meta)
            if feature is not None:
                rows.append(feature)
    return rows


def _parse_gene_hit_row(row: dict[str, str], hit_type: str, is_meta: bool) -> IntervalFeature | None:
    if is_meta:
        file_name = str(row.get("#FILE") or "").strip()
        sequence = str(row.get("SEQUENCE") or "").strip()
        if not file_name or not sequence:
            return None
        sequence_id = f"{Path(file_name).stem}_{sequence}"
        gene_name = str(row.get("GENE") or row.get("耐药基因") or row.get("毒力基因") or "").strip()
        product = str(row.get("PRODUCT") or row.get("毒力功能") or row.get("耐药机制") or "").strip()
        start = _safe_int(row.get("START"))
        end = _safe_int(row.get("END"))
    else:
        sequence_id = str(row.get("Contig名称") or "").strip()
        gene_name = str(row.get("基因名称") or "").strip()
        product = str(row.get("产物") or row.get("VF名称") or "").strip()
        start = _safe_int(row.get("起始碱基"))
        end = _safe_int(row.get("终止碱基"))
    if not sequence_id or not gene_name or start is None or end is None:
        return None
    return IntervalFeature(
        sequence_id=sequence_id,
        start=min(start, end),
        end=max(start, end),
        label=gene_name,
        source=hit_type,
        feature_type=hit_type.lower(),
        metadata=product,
    )


def _load_orf_features(path: Path) -> dict[str, list[IntervalFeature]]:
    features: dict[str, list[IntervalFeature]] = {}
    if not path.is_file():
        return features
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "CDS":
                continue
            start = _safe_int(fields[3])
            end = _safe_int(fields[4])
            if start is None or end is None:
                continue
            attrs = _parse_gff_attributes(fields[8])
            label = attrs.get("gene") or attrs.get("Name") or attrs.get("locus_tag") or attrs.get("ID") or "CDS"
            feature = IntervalFeature(
                sequence_id=fields[0].strip(),
                start=min(start, end),
                end=max(start, end),
                label=label,
                source="GFF",
                feature_type="orf",
                metadata=attrs.get("product", ""),
            )
            features.setdefault(feature.sequence_id, []).append(feature)
    for seq_id in features:
        features[seq_id].sort(key=lambda item: (item.start, item.end, item.label))
    return features


def _core_category(feature: IntervalFeature) -> str:
    text = f"{feature.label} {feature.metadata} {feature.source} {feature.feature_type}".lower()
    if feature.feature_type == "boundary":
        return "boundary"
    if any(keyword in text for keyword in CORE_IE_KEYWORDS):
        return "ie"
    if any(keyword in text for keyword in CORE_TRANSFER_KEYWORDS):
        return "transfer"
    if any(keyword in text for keyword in CORE_PHAGE_KEYWORDS):
        return "phage"
    if any(keyword in text for keyword in RRR_STD_KEYWORDS):
        return "rrr_std"
    if feature.source.lower() == "genomad" or feature.source.lower() == "mobileog":
        return "boundary"
    return ""


def _find_nearest_core(hit: IntervalFeature, mge_features: list[IntervalFeature]) -> tuple[IntervalFeature | None, int, int | None, list[str], bool]:
    same_contig = [feature for feature in mge_features if feature.sequence_id == hit.sequence_id]
    categories = {_core_category(feature) for feature in same_contig}
    categories.discard("")
    overlapping_boundaries = [
        feature for feature in same_contig
        if feature.feature_type == "boundary" and _interval_distance(hit.start, hit.end, feature.start, feature.end) == 0
    ]
    nearest_feature: IntervalFeature | None = None
    nearest_distance: int | None = None
    for feature in same_contig:
        if feature.start == hit.start and feature.end == hit.end and feature.label == hit.label:
            continue
        distance = _interval_distance(hit.start, hit.end, feature.start, feature.end)
        if nearest_distance is None or distance < nearest_distance:
            nearest_feature = feature
            nearest_distance = distance
    second_feature = len(categories) >= 2 or bool(overlapping_boundaries)
    return nearest_feature, nearest_distance if nearest_distance is not None else 10**9, None, sorted(categories), second_feature


def _orf_gap(hit: IntervalFeature, target: IntervalFeature | None, orf_map: dict[str, list[IntervalFeature]]) -> int | None:
    if target is None:
        return None
    orfs = orf_map.get(hit.sequence_id) or []
    if not orfs:
        return None
    hit_indexes = [idx for idx, orf in enumerate(orfs) if _interval_distance(hit.start, hit.end, orf.start, orf.end) == 0]
    target_indexes = [idx for idx, orf in enumerate(orfs) if _interval_distance(target.start, target.end, orf.start, orf.end) == 0]
    if not hit_indexes or not target_indexes:
        return None
    return min(abs(a - b) for a in hit_indexes for b in target_indexes)


def _risk_level(
    hit: IntervalFeature,
    nearest_feature: IntervalFeature | None,
    nearest_distance: int,
    orf_distance: int | None,
    categories: list[str],
    second_feature: bool,
    has_boundary_overlap: bool,
) -> tuple[str, str]:
    near_core = nearest_feature is not None and _core_category(nearest_feature) in {"ie", "transfer", "phage"}
    only_rrr_std = bool(categories) and set(categories).issubset({"rrr_std", "boundary"})
    if has_boundary_overlap or (near_core and nearest_distance <= 5000 and second_feature):
        return "A", "highly likely mobile / potentially transferable"
    if (near_core and nearest_distance <= 10000) or (orf_distance is not None and orf_distance <= 10):
        return "B", "possibly mobilizable / putatively associated with MGEs"
    if nearest_distance <= 25000 or only_rrr_std:
        return "C", "MGE-associated neighborhood"
    return "D", "weak evidence"


def summarize_mge_risk(pre: str, cwd: Path, cfg: GeNomadConfig) -> Path:
    integrated_path, card_path, vfdb_path, gff_path = _resolve_annotation_sources(pre, cwd, cfg)
    is_meta = _is_meta_method()
    mge_features = _load_interval_features(integrated_path, pre)
    arg_hits = _load_gene_hits(card_path, "ARG", is_meta)
    vf_hits = _load_gene_hits(vfdb_path, "VF", is_meta)
    orf_map = _load_orf_features(gff_path)

    rows: list[dict[str, str]] = []
    for hit in arg_hits + vf_hits:
        same_contig = [feature for feature in mge_features if feature.sequence_id == hit.sequence_id]
        has_boundary_overlap = any(
            feature.feature_type == "boundary" and _interval_distance(hit.start, hit.end, feature.start, feature.end) == 0
            for feature in same_contig
        )
        nearest_feature, nearest_distance, _, categories, second_feature = _find_nearest_core(hit, mge_features)
        orf_distance = _orf_gap(hit, nearest_feature, orf_map)
        level, statement = _risk_level(
            hit,
            nearest_feature,
            nearest_distance,
            orf_distance,
            categories,
            second_feature,
            has_boundary_overlap,
        )
        rows.append(
            {
                "sample": pre,
                "gene_type": hit.source,
                "gene_name": hit.label,
                "product": hit.metadata,
                "sequence_id": hit.sequence_id,
                "gene_start": str(hit.start),
                "gene_end": str(hit.end),
                "risk_level": level,
                "risk_statement": statement,
                "within_mge_boundary": "yes" if has_boundary_overlap else "no",
                "nearest_core_gene": nearest_feature.label if nearest_feature else "",
                "nearest_core_type": _core_category(nearest_feature) if nearest_feature else "",
                "nearest_core_distance_bp": "" if nearest_feature is None else str(nearest_distance),
                "nearest_core_distance_orf": "" if orf_distance is None else str(orf_distance),
                "same_contig_core_categories": ",".join(categories),
                "has_second_core_or_boundary": "yes" if second_feature else "no",
                "integrated_mge_tsv": str(integrated_path),
                "gene_source_tsv": str(card_path if hit.source == "ARG" else vfdb_path),
                "gff_path": str(gff_path),
            }
        )

    out_path = cwd / f"{pre}.mge_risk_summary.tsv"
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=[
                "sample",
                "gene_type",
                "gene_name",
                "product",
                "sequence_id",
                "gene_start",
                "gene_end",
                "risk_level",
                "risk_statement",
                "within_mge_boundary",
                "nearest_core_gene",
                "nearest_core_type",
                "nearest_core_distance_bp",
                "nearest_core_distance_orf",
                "same_contig_core_categories",
                "has_second_core_or_boundary",
                "integrated_mge_tsv",
                "gene_source_tsv",
                "gff_path",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)
    return out_path


def _resolve_mobileog_query_faa(pre: str, input_fasta: Path, cwd: Path) -> Path:
    del pre, cwd
    return input_fasta


def _run_mobileog(pre: str, input_fasta: Path, threads: int, cwd: Path, cfg: GeNomadConfig) -> Path | None:
    query_fasta = _resolve_mobileog_query_faa(pre, input_fasta, cwd)
    meta_csv = Path(os.environ.get("META_MOBILEOG_META", MOBILEOG_META_DEFAULT))
    db_path = os.environ.get("META_MOBILEOG_DB", MOBILEOG_DB_DEFAULT)
    if not query_fasta.is_file() or not meta_csv.is_file():
        return None
    blast_tsv = cwd / f"{pre}.mgeblast.tsv"
    out_tsv = cwd / f"{pre}.mobileog.tsv"
    _run_shell_command_with_cfg(
        cfg,
        (
            f"diamond blastx -q {query_fasta} --db {db_path} "
            f"--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore "
            f"--out {pre}.mgeblast.tsv --threads {threads}"
        ),
        cwd / "mobileog.stdout.log",
        cwd / "mobileog.stderr.log",
        cwd=cwd,
    )
    metadata = _load_mobileog_metadata(meta_csv)
    rows: list[dict[str, str]] = []
    if blast_tsv.is_file():
        with blast_tsv.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            for fields in reader:
                if len(fields) != 12:
                    continue
                pident = float(fields[2])
                evalue = float(fields[10])
                if pident <= 25 or evalue >= 1e-5:
                    continue
                subject = fields[1]
                entry_id = subject.split("|")[0]
                meta = metadata.get(entry_id, {})
                rows.append(
                    {
                        "序列名称": fields[0],
                        "参考基因组名称": subject,
                        "类型": "|".join(subject.split("|")[3:]),
                        "注释": str(meta.get("Manual Annotation") or ""),
                        "元件名称": str(meta.get("Name") or ""),
                        "相似性(%)": fields[2],
                        "长度": fields[3],
                        "差异数量": fields[4],
                        "空缺数量": fields[5],
                        "序列起始": fields[6],
                        "序列终止": fields[7],
                        "参考起始": fields[8],
                        "参考终止": fields[9],
                        "evalue": fields[10],
                        "比对得分": fields[11],
                    }
                )
    with out_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=[
                "序列名称",
                "参考基因组名称",
                "类型",
                "注释",
                "元件名称",
                "相似性(%)",
                "长度",
                "差异数量",
                "空缺数量",
                "序列起始",
                "序列终止",
                "参考起始",
                "参考终止",
                "evalue",
                "比对得分",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)
    return out_tsv


def _resolve_annoele_fasta_path(pre: str, cwd: Path) -> Path:
    try:
        from metagenomic_refactor.context import get_runtime_context
        runtime = get_runtime_context()
        method = str(runtime.method or "").strip()
    except Exception:
        method = ""
    if method == "meta":
        return cwd / "tmp_combine.fa"
    return cwd / f"{pre}.final.fasta"


def AnnoEle(pre: str, threads: int) -> Path:
    cwd = Path.cwd()
    fasta_path = _resolve_annoele_fasta_path(pre, cwd)
    if not fasta_path.is_file():
        raise GeNomadMGEError(f"未找到元件预测所需的fasta文件: {fasta_path}")

    db_path = Path(os.environ.get("META_GENOMAD_DB", GENOMAD_DB_DEFAULT))
    cfg = GeNomadConfig(
        outdir=cwd,
        database=db_path,
        threads=threads,
        conda_env=os.environ.get("META_MGE_ENV", "genomad_aux"),
    )
    cfg2 = GeNomadConfig(
        outdir=cwd,
        database=db_path,
        threads=threads,
        conda_env=os.environ.get("META_MGE_ENV", "genomad_aux"),
    )
    if db_path.is_dir():
        sample = GeNomadSample(sample=pre, fasta=fasta_path)
        run_genomad_sample(sample, cfg)

    _run_mobileog(pre, fasta_path, threads, cwd, cfg2)
    integrate_mge_tables(cwd, pre)
    return summarize_mge_risk(pre, cwd, cfg2)


def run_genomad_sample(sample: GeNomadSample, cfg: GeNomadConfig) -> Path:
    sample_out = cfg.outdir.expanduser().resolve() / sample.sample
    genomad_out = sample_out / "genomad"
    if cfg.force:
        _safe_clear_dir(sample_out)
    sample_out.mkdir(parents=True, exist_ok=True)

    summary_tsv = sample_out / f"{sample.sample}.genomad_mge_summary.tsv"
    if not cfg.force and summary_tsv.is_file():
        return summary_tsv

    _run_command(
        [
            "genomad",
            "end-to-end",
            str(sample.fasta),
            str(genomad_out),
            str(cfg.database),
            "--threads",
            str(cfg.threads),
        ],
        sample_out / "genomad.stdout.log",
        sample_out / "genomad.stderr.log",
    )
    return aggregate_genomad_results(genomad_out, sample.sample, sample_out)


def run_genomad(samples: Sequence[GeNomadSample], cfg: GeNomadConfig) -> Path:
    outdir = cfg.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    run_summary = outdir / "run_summary.tsv"
    with run_summary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "status", "output_dir", "mge_summary_tsv"])
        for sample in samples:
            summary_tsv = run_genomad_sample(sample, cfg)
            writer.writerow([sample.sample, "done", str(outdir / sample.sample), str(summary_tsv)])
    return run_summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="基于geNomad预测移动遗传元件(MGE)，包括质粒和前噬菌体(provirus)"
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--manifest", help="TSV样本表，包含sample和fasta两列")
    group.add_argument("--fasta", help="单个样本的装配fasta文件")
    parser.add_argument("--sample", help="单样本模式下的样本名")
    parser.add_argument("--outdir", required=True, help="输出目录")
    parser.add_argument("--database", required=True, help="geNomad数据库目录")
    parser.add_argument("--threads", type=int, default=16, help="线程数，默认16")
    parser.add_argument(
        "--conda-env",
        default="genomad_aux",
        help="除geNomad外其他软件(如mobileOG)使用的Conda环境名，默认genomad_aux",
    )
    parser.add_argument("--force", action="store_true", help="覆盖已有结果并重跑")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.fasta and not args.sample:
        raise SystemExit("单样本模式下请同时提供 --sample")
    samples = load_samples(manifest=args.manifest, sample=args.sample, fasta=args.fasta)
    cfg = GeNomadConfig(
        outdir=Path(args.outdir),
        database=_ensure_dir(Path(args.database), "geNomad数据库目录"),
        threads=args.threads,
        conda_env=args.conda_env,
        force=args.force,
    )
    summary = run_genomad(samples, cfg)
    print(f"geNomad MGE预测完成，结果汇总: {summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
