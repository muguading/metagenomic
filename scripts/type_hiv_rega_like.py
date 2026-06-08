#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import io
import json
import math
import re
import subprocess
import tempfile
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path

try:
    from Bio import Phylo
except ModuleNotFoundError as exc:
    raise SystemExit(
        "Biopython is required. Please run this script with the ncov environment Python: "
        "/opt/homebrew/Caskroom/mambaforge/base/envs/ncov/bin/python3"
    ) from exc


DEFAULT_REFERENCE_FASTA = Path(
    "/Users/wuhhh/Desktop/徐老师/代码/metagenomic/database/virus/HIV/rega_reference_genomes/rega_hiv_reference_genomes.fasta"
)
DEFAULT_REFERENCE_MANIFEST = Path(
    "/Users/wuhhh/Desktop/徐老师/代码/metagenomic/database/virus/HIV/rega_reference_genomes/reference_manifest.tsv"
)
DEFAULT_SUPPLEMENT_MANIFEST = Path(
    "/Users/wuhhh/Desktop/徐老师/代码/metagenomic/database/virus/HIV/rega_reference_genomes/supplement_manifest.tsv"
)
DEFAULT_RULES_JSON = Path(
    "/Users/wuhhh/Desktop/徐老师/代码/metagenomic/database/virus/HIV/rega_subtyping_rules.json"
)
DEFAULT_FASTTREE = Path("/opt/homebrew/Caskroom/mambaforge/base/envs/ncov/bin/FastTree")
DEFAULT_SEQKIT = Path("/opt/homebrew/Caskroom/mambaforge/base/envs/ncov/bin/seqkit")

PURE_GROUP_ORDER = ["A1", "A2", "B", "C", "D", "F1", "F2", "G", "H", "J", "K"]
BOOTSCAN_WINDOW = 400
BOOTSCAN_STEP = 40
MIN_WINDOW_FRACTION = 0.10
MIN_COMPARABLE_BASES = 200
MIN_IDENTITY_STRONG = 0.85
MIN_IDENTITY_WEAK = 0.80
MIN_MARGIN_STRONG = 0.015
MIN_MARGIN_WEAK = 0.008


@dataclass
class ReferenceRecord:
    header: str
    description: str
    accession: str
    group: str
    subtype: str
    is_crf: bool
    sequence: str


@dataclass
class GroupScore:
    group: str
    max_identity: float
    mean_identity: float
    support_count: int
    best_reference: str


@dataclass
class WindowAssignment:
    start: int
    end: int
    scope: str
    top_group: str
    top_identity: float
    second_group: str
    second_identity: float
    margin: float
    supported: bool


@dataclass
class TreePlacement:
    scope: str
    cluster_group: str
    cluster_mode: str
    support: float
    nearest_reference: str
    nearest_group: str
    member_count: int


@dataclass
class WindowProfile:
    start: int
    end: int
    midpoint: int
    scope: str
    supports: dict[str, float]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Approximate REGA/Stanford HIV-1 subtyping using local reference genomes and windowed similarity."
    )
    parser.add_argument("fasta", type=Path, help="Input HIV-1 consensus FASTA.")
    parser.add_argument("--reference-fasta", type=Path, default=DEFAULT_REFERENCE_FASTA)
    parser.add_argument("--reference-manifest", type=Path, default=DEFAULT_REFERENCE_MANIFEST)
    parser.add_argument("--supplement-manifest", type=Path, default=DEFAULT_SUPPLEMENT_MANIFEST)
    parser.add_argument("--rules-json", type=Path, default=DEFAULT_RULES_JSON)
    parser.add_argument("--fasttree-bin", type=Path, default=DEFAULT_FASTTREE)
    parser.add_argument("--seqkit-bin", type=Path, default=DEFAULT_SEQKIT)
    parser.add_argument("--output-dir", type=Path, default=None, help="Optional directory for CSV/SVG/PNG outputs.")
    parser.add_argument("--no-assets", action="store_true", help="Skip per-sample bootscan CSV/SVG generation for faster batch testing.")
    parser.add_argument("--json", action="store_true", help="Emit JSON.")
    return parser.parse_args()


def parse_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    seq_parts: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts).upper()))
                header = line[1:].strip() or "unnamed"
                seq_parts = []
                continue
            seq_parts.append(re.sub(r"[^A-Za-z-]", "", line))
    if header is not None:
        records.append((header, "".join(seq_parts).upper()))
    if not records:
        raise SystemExit(f"Empty or invalid FASTA: {path}")
    return records


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def run_mafft_alignment(reference_records: list[tuple[str, str]], query_records: list[tuple[str, str]]) -> dict[str, str]:
    mafft = shutil_which("mafft")
    if not mafft:
        raise SystemExit("mafft not found in PATH")
    with tempfile.TemporaryDirectory(prefix="hiv_rega_") as tmpdir:
        tmpdir_path = Path(tmpdir)
        combined_path = tmpdir_path / "combined.fasta"
        write_fasta(combined_path, [*reference_records, *query_records])
        result = subprocess.run(
            [mafft, "--retree", "1", "--maxiterate", "0", "--quiet", str(combined_path)],
            capture_output=True,
            text=True,
            check=True,
        )
    aligned: dict[str, str] = {}
    header: str | None = None
    seq_parts: list[str] = []
    for raw_line in result.stdout.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                aligned[header] = "".join(seq_parts).upper()
            header = line[1:].strip()
            seq_parts = []
            continue
        seq_parts.append(line)
    if header is not None:
        aligned[header] = "".join(seq_parts).upper()
    return aligned


def shutil_which(program: str) -> str:
    result = subprocess.run(["/usr/bin/which", program], capture_output=True, text=True)
    return result.stdout.strip()


def normalize_group_label(group: str) -> tuple[str, bool]:
    raw = str(group).strip()
    if raw.startswith("Subtype "):
        return raw.split(" ", 1)[1].strip(), False
    return raw, raw.upper().startswith("CRF")


def infer_subtype_from_group(group: str) -> str:
    normalized, _is_crf = normalize_group_label(group)
    return normalized.replace("-", "_").upper()


def resolve_fasta_header_by_accession(fasta_records: dict[str, str], accession: str) -> str:
    accession = accession.strip()
    if accession in fasta_records:
        return accession
    return next((header for header in fasta_records if header.startswith(accession + ".")), "")


def load_reference_metadata(
    reference_fasta: Path,
    reference_manifest: Path,
    supplement_manifest: Path,
) -> list[ReferenceRecord]:
    fasta_records = dict(parse_fasta(reference_fasta))
    accession_to_record: dict[str, ReferenceRecord] = {}

    with reference_manifest.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if str(row.get("download_status") or "").strip() != "downloaded":
                continue
            accession = str(row.get("note") or "").strip()
            header = resolve_fasta_header_by_accession(fasta_records, accession)
            if not header:
                continue
            group = str(row.get("group") or "").strip()
            subtype = infer_subtype_from_group(group)
            is_crf = subtype.upper().startswith("CRF")
            accession_to_record[accession] = ReferenceRecord(
                header=accession,
                description=header,
                accession=accession,
                group=group,
                subtype=subtype,
                is_crf=is_crf,
                sequence=fasta_records[header],
            )

    with supplement_manifest.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            accession = str(row.get("accession") or "").strip()
            header = resolve_fasta_header_by_accession(fasta_records, accession)
            if not header:
                continue
            group = str(row.get("group") or "").strip()
            subtype = str(row.get("subtype") or "").strip() or infer_subtype_from_group(group)
            accession_to_record[accession] = ReferenceRecord(
                header=accession,
                description=header,
                accession=accession,
                group=group,
                subtype=subtype.replace("-", "_").upper(),
                is_crf=group.upper().startswith("CRF"),
                sequence=fasta_records[header],
            )

    return list(accession_to_record.values())


def comparable_bases(seq_a: str, seq_b: str) -> tuple[int, int]:
    matches = 0
    comparable = 0
    for base_a, base_b in zip(seq_a.upper(), seq_b.upper()):
        if "-" in {base_a, base_b}:
            continue
        if "N" in {base_a, base_b}:
            continue
        comparable += 1
        if base_a == base_b:
            matches += 1
    return matches, comparable


def normalize_fasta_with_seqkit(seqkit_bin: Path, source_fasta: Path) -> list[tuple[str, str]]:
    result = subprocess.run(
        [str(seqkit_bin), "seq", "-w", "0", str(source_fasta)],
        capture_output=True,
        text=True,
        check=True,
    )
    with tempfile.TemporaryDirectory(prefix="hiv_seqkit_") as tmpdir:
        tmp_path = Path(tmpdir) / "normalized.fasta"
        tmp_path.write_text(result.stdout, encoding="utf-8")
        return parse_fasta(tmp_path)


def identity_for_pair(seq_a: str, seq_b: str) -> float:
    matches, comparable = comparable_bases(seq_a, seq_b)
    return (matches / comparable) if comparable else 0.0


def group_scores_for_query(
    aligned_query: str,
    aligned_references: dict[str, str],
    metadata: list[ReferenceRecord],
    scope: str,
) -> list[GroupScore]:
    bucket: dict[str, list[tuple[float, str]]] = defaultdict(list)
    for reference in metadata:
        if scope == "pure" and reference.is_crf:
            continue
        identity = identity_for_pair(aligned_query, aligned_references[reference.header])
        bucket[reference.subtype].append((identity, reference.header))

    scores: list[GroupScore] = []
    for group, values in bucket.items():
        ordered = sorted(values, key=lambda item: item[0], reverse=True)
        scores.append(
            GroupScore(
                group=group,
                max_identity=ordered[0][0],
                mean_identity=sum(item[0] for item in ordered) / len(ordered),
                support_count=len(values),
                best_reference=ordered[0][1],
            )
        )
    scores.sort(key=lambda item: (item.max_identity, item.mean_identity), reverse=True)
    return scores


def query_position_columns(aligned_query: str) -> list[int]:
    columns: list[int] = []
    for index, base in enumerate(aligned_query):
        if base == "-" or base.upper() == "N":
            continue
        columns.append(index)
    return columns


def slice_alignment_by_query_window(aligned_sequence: str, query_columns: list[int], start: int, end: int) -> str:
    selected = set(query_columns[start:end])
    return "".join(base for index, base in enumerate(aligned_sequence) if index in selected)


def assign_windows(
    aligned_query: str,
    aligned_references: dict[str, str],
    metadata: list[ReferenceRecord],
    scope: str,
    window_bp: int = BOOTSCAN_WINDOW,
    step_bp: int = BOOTSCAN_STEP,
) -> list[WindowAssignment]:
    columns = query_position_columns(aligned_query)
    if len(columns) < min(window_bp, MIN_COMPARABLE_BASES):
        return []
    windows: list[WindowAssignment] = []
    for start in range(0, max(len(columns) - window_bp + 1, 1), step_bp):
        end = min(start + window_bp, len(columns))
        if end - start < MIN_COMPARABLE_BASES:
            continue
        query_window = slice_alignment_by_query_window(aligned_query, columns, start, end)
        bucket: dict[str, list[float]] = defaultdict(list)
        for reference in metadata:
            if scope == "pure" and reference.is_crf:
                continue
            ref_window = slice_alignment_by_query_window(aligned_references[reference.header], columns, start, end)
            matches, comparable = comparable_bases(query_window, ref_window)
            if comparable < MIN_COMPARABLE_BASES:
                continue
            bucket[reference.subtype].append(matches / comparable)
        if not bucket:
            continue
        ranked = sorted(
            (
                (group, max(values), sum(values) / len(values))
                for group, values in bucket.items()
            ),
            key=lambda item: (item[1], item[2]),
            reverse=True,
        )
        top_group, top_identity, _top_mean = ranked[0]
        second_group, second_identity = ("", 0.0)
        if len(ranked) > 1:
            second_group, second_identity, _second_mean = ranked[1]
        margin = top_identity - second_identity
        supported = top_identity >= MIN_IDENTITY_WEAK and margin >= MIN_MARGIN_WEAK
        windows.append(
            WindowAssignment(
                start=start + 1,
                end=end,
                scope=scope,
                top_group=top_group,
                top_identity=top_identity,
                second_group=second_group,
                second_identity=second_identity,
                margin=margin,
                supported=supported,
            )
        )
    return windows


def build_window_profiles(
    aligned_query: str,
    aligned_references: dict[str, str],
    metadata: list[ReferenceRecord],
    scope: str,
    window_bp: int = BOOTSCAN_WINDOW,
    step_bp: int = BOOTSCAN_STEP,
    temperature: float = 250.0,
) -> list[WindowProfile]:
    columns = query_position_columns(aligned_query)
    if len(columns) < min(window_bp, MIN_COMPARABLE_BASES):
        return []
    profiles: list[WindowProfile] = []
    for start in range(0, max(len(columns) - window_bp + 1, 1), step_bp):
        end = min(start + window_bp, len(columns))
        if end - start < MIN_COMPARABLE_BASES:
            continue
        query_window = slice_alignment_by_query_window(aligned_query, columns, start, end)
        raw_scores: dict[str, float] = {}
        for reference in metadata:
            if scope == "pure" and reference.is_crf:
                continue
            ref_window = slice_alignment_by_query_window(aligned_references[reference.header], columns, start, end)
            matches, comparable = comparable_bases(query_window, ref_window)
            if comparable < MIN_COMPARABLE_BASES:
                continue
            identity = matches / comparable
            current = raw_scores.get(reference.subtype)
            if current is None or identity > current:
                raw_scores[reference.subtype] = identity
        if not raw_scores:
            continue
        max_identity = max(raw_scores.values())
        weights = {
            group: math.exp((identity - max_identity) * temperature)
            for group, identity in raw_scores.items()
        }
        weight_sum = sum(weights.values()) or 1.0
        supports = {group: (weight / weight_sum) * 100.0 for group, weight in weights.items()}
        profiles.append(
            WindowProfile(
                start=start + 1,
                end=end,
                midpoint=(start + end) // 2,
                scope=scope,
                supports=supports,
            )
        )
    return profiles


def summarize_window_parents(windows: list[WindowAssignment]) -> list[tuple[str, int, float]]:
    supported = [window for window in windows if window.supported and window.top_group]
    total = len(supported)
    if total == 0:
        return []
    counts = Counter(window.top_group for window in supported)
    ranked = sorted(counts.items(), key=lambda item: item[1], reverse=True)
    return [(group, count, count / total) for group, count in ranked]


def top_two(scores: list[GroupScore]) -> tuple[GroupScore | None, GroupScore | None]:
    top = scores[0] if scores else None
    second = scores[1] if len(scores) > 1 else None
    return top, second


def build_fasttree(
    fasttree_bin: Path,
    aligned: dict[str, str],
    references: list[ReferenceRecord],
    query_ids: list[str],
    scope: str,
) -> object:
    selected = [reference for reference in references if scope != "pure" or not reference.is_crf]
    records = [(reference.header, aligned[reference.header]) for reference in selected]
    records.extend((query_id, aligned[query_id]) for query_id in query_ids)
    with tempfile.TemporaryDirectory(prefix=f"hiv_fasttree_{scope}_") as tmpdir:
        aligned_fasta = Path(tmpdir) / f"{scope}.aligned.fasta"
        write_fasta(aligned_fasta, records)
        result = subprocess.run(
            [str(fasttree_bin), "-nt", str(aligned_fasta)],
            capture_output=True,
            text=True,
            check=True,
        )
    return Phylo.read(io.StringIO(result.stdout), "newick")


def tree_placement_from_tree(
    tree: object,
    aligned_query_id: str,
    references: list[ReferenceRecord],
    scope: str,
) -> TreePlacement:
    selected = [reference for reference in references if scope != "pure" or not reference.is_crf]
    group_map = {reference.header: reference.subtype for reference in selected}
    query_terminal = next(clade for clade in tree.get_terminals() if clade.name == aligned_query_id)
    nearest_reference = ""
    nearest_group = ""
    nearest_distance = None
    for terminal in tree.get_terminals():
        if terminal.name == aligned_query_id:
            continue
        distance = tree.distance(query_terminal, terminal)
        if nearest_distance is None or distance < nearest_distance:
            nearest_distance = distance
            nearest_reference = terminal.name or ""
            nearest_group = group_map.get(nearest_reference, "")

    ancestor_path = tree.get_path(query_terminal)
    for ancestor in reversed(ancestor_path):
        names = [terminal.name for terminal in ancestor.get_terminals() if terminal.name != aligned_query_id]
        if not names:
            continue
        groups = [group_map.get(name, "") for name in names if group_map.get(name, "")]
        if not groups:
            continue
        group_counts = Counter(groups)
        cluster_group, cluster_count = group_counts.most_common(1)[0]
        support = 0.0
        if getattr(ancestor, "confidence", None) is not None:
            try:
                support = float(ancestor.confidence)
            except (TypeError, ValueError):
                support = 0.0
        elif ancestor.name:
            try:
                support = float(str(ancestor.name))
            except ValueError:
                support = 0.0
        cluster_mode = "inside" if len(group_counts) == 1 else "near"
        return TreePlacement(
            scope=scope,
            cluster_group=cluster_group,
            cluster_mode=cluster_mode,
            support=support,
            nearest_reference=nearest_reference,
            nearest_group=nearest_group,
            member_count=len(names),
        )

    return TreePlacement(
        scope=scope,
        cluster_group=nearest_group,
        cluster_mode="nearest",
        support=0.0,
        nearest_reference=nearest_reference,
        nearest_group=nearest_group,
        member_count=1 if nearest_reference else 0,
    )


def classify_assignment(
    sequence_length: int,
    pure_tree: TreePlacement,
    all_tree: TreePlacement,
    pure_scores: list[GroupScore],
    all_scores: list[GroupScore],
    pure_windows: list[WindowAssignment],
    all_windows: list[WindowAssignment],
) -> tuple[str, str, list[str]]:
    pure_top, pure_second = top_two(pure_scores)
    all_top, all_second = top_two(all_scores)
    parent_summary = summarize_window_parents(all_windows)
    notes: list[str] = []
    if pure_tree.cluster_group:
        notes.append(f"pure tree {pure_tree.cluster_mode} {pure_tree.cluster_group} support={pure_tree.support:.2f}")
    if all_tree.cluster_group:
        notes.append(f"overall tree {all_tree.cluster_mode} {all_tree.cluster_group} support={all_tree.support:.2f}")
    if pure_top:
        notes.append(f"top pure {pure_top.group} identity={pure_top.max_identity:.4f}")
    if all_top:
        notes.append(f"top overall {all_top.group} identity={all_top.max_identity:.4f}")

    recomb_parents = [item for item in parent_summary if item[2] >= MIN_WINDOW_FRACTION]
    recombination_detected = len(recomb_parents) >= 2
    if recombination_detected:
        notes.append("window scan supports >=2 parental groups")

    pure_margin = (pure_top.max_identity - pure_second.max_identity) if pure_top and pure_second else (pure_top.max_identity if pure_top else 0.0)
    all_margin = (all_top.max_identity - all_second.max_identity) if all_top and all_second else (all_top.max_identity if all_top else 0.0)

    if sequence_length <= 800:
        if all_tree.cluster_group.upper().startswith("CRF") and (all_tree.cluster_mode == "inside" or (all_top and all_top.max_identity >= MIN_IDENTITY_WEAK and all_margin >= MIN_MARGIN_WEAK)):
            return "CRF", all_tree.cluster_group, notes
        if pure_tree.cluster_group and (pure_tree.cluster_mode == "inside" or (pure_top and pure_top.max_identity >= MIN_IDENTITY_WEAK and pure_margin >= MIN_MARGIN_WEAK)):
            return "PURE", pure_tree.cluster_group, notes
        return "NOT ASSIGNED", (all_top.group if all_top else (pure_top.group if pure_top else "")), notes

    if recombination_detected:
        parent_groups = [group for group, _count, _fraction in recomb_parents]
        has_crf = any(group.upper().startswith("CRF") for group in parent_groups)
        has_pure = any(not group.upper().startswith("CRF") for group in parent_groups)
        if has_crf and has_pure:
            return "CRF PURE RECOMBINANT", (all_top.group if all_top else parent_groups[0]), notes
        if not has_crf and len(parent_groups) >= 2:
            return "PURE RECOMBINANT", parent_groups[0], notes
        return "POTENTIAL RECOMBINANT", (all_top.group if all_top else parent_groups[0]), notes

    if all_tree.cluster_group.upper().startswith("CRF"):
        if all_tree.cluster_mode == "inside" and ((all_top and all_top.max_identity >= MIN_IDENTITY_STRONG) or all_tree.support >= 0.70):
            return "CRF", all_tree.cluster_group, notes
        if all_top and all_top.max_identity >= MIN_IDENTITY_WEAK:
            return "CRF-LIKE", all_tree.cluster_group, notes

    if pure_tree.cluster_group and pure_tree.cluster_mode == "inside" and ((pure_top and pure_top.max_identity >= MIN_IDENTITY_STRONG) or pure_tree.support >= 0.70):
        return "PURE", pure_tree.cluster_group, notes
    if pure_tree.cluster_group and pure_top and pure_top.max_identity >= MIN_IDENTITY_WEAK:
        return "PURE-LIKE", pure_tree.cluster_group, notes
    return "NOT ASSIGNED", (all_top.group if all_top else (pure_top.group if pure_top else "")), notes


def format_group_scores(scores: list[GroupScore], top_n: int = 5) -> list[dict[str, object]]:
    return [
        {
            "group": item.group,
            "max_identity": round(item.max_identity, 6),
            "mean_identity": round(item.mean_identity, 6),
            "support_count": item.support_count,
            "best_reference": item.best_reference,
        }
        for item in scores[:top_n]
    ]


def top_plot_groups(profiles: list[WindowProfile], preferred: list[str] | None = None, limit: int = 10) -> list[str]:
    totals: dict[str, float] = defaultdict(float)
    for profile in profiles:
        for group, support in profile.supports.items():
            totals[group] += support
    ordered = [group for group, _value in sorted(totals.items(), key=lambda item: item[1], reverse=True)]
    chosen: list[str] = []
    seen: set[str] = set()
    for group in preferred or []:
        if group in totals and group not in seen:
            chosen.append(group)
            seen.add(group)
    for group in ordered:
        if group not in seen:
            chosen.append(group)
            seen.add(group)
        if len(chosen) >= limit:
            break
    return chosen[:limit]


def svg_escape(text: str) -> str:
    return (
        str(text)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def color_for_group(group: str) -> str:
    palette = {
        "A1": "#ff2a2a",
        "A2": "#e84d8a",
        "B": "#2da1ff",
        "C": "#a16a2b",
        "D": "#f3a6d6",
        "F1": "#b8ff00",
        "F2": "#86d800",
        "G": "#68d37e",
        "H": "#ffbf00",
        "J": "#1fe0ea",
        "K": "#a759ff",
    }
    if group in palette:
        return palette[group]
    if group.startswith("CRF"):
        return "#435f88" if group == "CRF01_AE" else "#4e6f9f"
    return "#666666"


def write_bootscan_csv(path: Path, profiles: list[WindowProfile], groups: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["midpoint", "start", "end", *groups])
        for profile in profiles:
            writer.writerow([
                profile.midpoint,
                profile.start,
                profile.end,
                *[f"{profile.supports.get(group, 0.0):.6f}" for group in groups],
            ])


def build_polyline_points(
    profiles: list[WindowProfile],
    group: str,
    chart_left: float,
    chart_top: float,
    chart_width: float,
    chart_height: float,
    x_max: int,
) -> str:
    coords: list[str] = []
    for profile in profiles:
        x = chart_left + (profile.midpoint / x_max) * chart_width
        y = chart_top + chart_height - (profile.supports.get(group, 0.0) / 100.0) * chart_height
        coords.append(f"{x:.2f},{y:.2f}")
    return " ".join(coords)


def render_bootscan_svg(
    path: Path,
    profiles: list[WindowProfile],
    groups: list[str],
    sample_name: str,
    support_score: float,
    title_suffix: str,
    predicted_label: str,
    x_max: int,
) -> None:
    width = 1500
    height = 980
    chart_left = 120
    chart_top = 110
    chart_width = 1080
    chart_height = 680
    legend_x = 1230
    baseline_y = chart_top + chart_height - (70.0 / 100.0) * chart_height
    x_tick_max = max(10000, int(math.ceil(x_max / 1000.0) * 1000))

    svg: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<rect x="0" y="0" width="100%" height="74" fill="#11149d"/>',
        '<text x="750" y="47" font-size="28" font-family="Arial, Helvetica, sans-serif" fill="#ffffff" text-anchor="middle" font-weight="700">HIV-1 CRF/亚型重组分析</text>',
        f'<text x="750" y="110" font-size="34" font-family="Arial, Helvetica, sans-serif" fill="#111111" text-anchor="middle">Bootscan analysis{svg_escape(title_suffix)}</text>',
        f'<text x="34" y="{chart_top + 10}" font-size="28" font-family="Arial, Helvetica, sans-serif" fill="#111111">Support</text>',
    ]

    svg.append(f'<line x1="{chart_left}" y1="{chart_top}" x2="{chart_left}" y2="{chart_top + chart_height}" stroke="#222222" stroke-width="3"/>')
    svg.append(f'<line x1="{chart_left}" y1="{chart_top + chart_height}" x2="{chart_left + chart_width}" y2="{chart_top + chart_height}" stroke="#222222" stroke-width="3"/>')

    for support in range(0, 101, 20):
        y = chart_top + chart_height - (support / 100.0) * chart_height
        svg.append(f'<line x1="{chart_left - 10}" y1="{y}" x2="{chart_left + 10}" y2="{y}" stroke="#222222" stroke-width="2"/>')
        svg.append(f'<text x="{chart_left - 20}" y="{y + 8}" font-size="20" font-family="Arial, Helvetica, sans-serif" text-anchor="end">{support}</text>')

    for support in range(10, 100, 10):
        if support % 20 == 0:
            continue
        y = chart_top + chart_height - (support / 100.0) * chart_height
        svg.append(f'<line x1="{chart_left - 6}" y1="{y}" x2="{chart_left + 6}" y2="{y}" stroke="#222222" stroke-width="1.4"/>')

    for tick in range(0, x_tick_max + 1, 1000):
        x = chart_left + (tick / x_tick_max) * chart_width
        svg.append(f'<line x1="{x}" y1="{chart_top + chart_height - 10}" x2="{x}" y2="{chart_top + chart_height + 10}" stroke="#222222" stroke-width="2"/>')
        svg.append(f'<text x="{x}" y="{chart_top + chart_height + 40}" font-size="20" font-family="Arial, Helvetica, sans-serif" text-anchor="middle">{tick}</text>')

    for tick in range(500, x_tick_max, 500):
        if tick % 1000 == 0:
            continue
        x = chart_left + (tick / x_tick_max) * chart_width
        svg.append(f'<line x1="{x}" y1="{chart_top + chart_height - 6}" x2="{x}" y2="{chart_top + chart_height + 6}" stroke="#222222" stroke-width="1.4"/>')

    svg.append(f'<line x1="{chart_left}" y1="{baseline_y}" x2="{chart_left + chart_width}" y2="{baseline_y}" stroke="#ff5a5a" stroke-width="2" stroke-dasharray="8,6"/>')

    for group in groups:
        points = build_polyline_points(profiles, group, chart_left, chart_top, chart_width, chart_height, x_tick_max)
        if not points:
            continue
        svg.append(
            f'<polyline fill="none" stroke="{color_for_group(group)}" stroke-width="2.6" '
            f'stroke-linejoin="round" stroke-linecap="round" points="{points}"/>'
        )

    svg.append(f'<text x="{chart_left + chart_width / 2}" y="{chart_top + chart_height + 82}" font-size="26" font-family="Arial, Helvetica, sans-serif" text-anchor="middle">Nucleotide position</text>')
    svg.append(f'<text x="750" y="855" font-size="24" font-family="Arial, Helvetica, sans-serif" text-anchor="middle">样本: {svg_escape(sample_name)}   判定: {svg_escape(predicted_label)}</text>')
    svg.append(f'<text x="750" y="892" font-size="24" font-family="Arial, Helvetica, sans-serif" text-anchor="middle">Bootscan 集群支持: {support_score:.3f}</text>')
    svg.append(f'<text x="750" y="929" font-size="22" font-family="Arial, Helvetica, sans-serif" text-anchor="middle">Bootscan 分析采用窗口大小 {BOOTSCAN_WINDOW} 和步长 {BOOTSCAN_STEP} 执行。</text>')

    for index, group in enumerate(groups):
        y = 250 + index * 42
        svg.append(f'<line x1="{legend_x}" y1="{y}" x2="{legend_x + 32}" y2="{y}" stroke="{color_for_group(group)}" stroke-width="5"/>')
        svg.append(f'<text x="{legend_x + 40}" y="{y + 7}" font-size="22" font-family="Arial, Helvetica, sans-serif">{svg_escape(group.replace("CRF", "").replace("_", "_"))}</text>')

    svg.append("</svg>")
    path.write_text("\n".join(svg) + "\n", encoding="utf-8")


def maybe_render_png_from_svg(svg_path: Path) -> Path | None:
    try:
        subprocess.run(
            ["/usr/bin/qlmanage", "-t", "-s", "2400", "-o", str(svg_path.parent), str(svg_path)],
            capture_output=True,
            text=True,
            check=True,
        )
    except Exception:
        return None
    candidate = svg_path.parent / f"{svg_path.name}.png"
    return candidate if candidate.is_file() else None


def main() -> int:
    args = parse_args()
    _rules = json.loads(args.rules_json.read_text(encoding="utf-8"))
    if not args.fasttree_bin.is_file():
        raise SystemExit(f"FastTree binary not found: {args.fasttree_bin}")
    if not args.seqkit_bin.is_file():
        raise SystemExit(f"seqkit binary not found: {args.seqkit_bin}")
    references = load_reference_metadata(args.reference_fasta, args.reference_manifest, args.supplement_manifest)
    reference_records = [(record.header, record.sequence) for record in references]
    input_records = normalize_fasta_with_seqkit(args.seqkit_bin, args.fasta)
    query_records = [(f"query_{index}", sequence) for index, (_header, sequence) in enumerate(input_records, start=1)]
    query_name_lookup = {f"query_{index}": header for index, (header, _sequence) in enumerate(input_records, start=1)}

    aligned = run_mafft_alignment(reference_records, query_records)
    reference_aligned = {record.header: aligned[record.header] for record in references if record.header in aligned}
    query_ids = [query_id for query_id, _sequence in query_records]
    pure_tree_global = build_fasttree(args.fasttree_bin, aligned, references, query_ids, scope="pure")
    overall_tree_global = build_fasttree(args.fasttree_bin, aligned, references, query_ids, scope="all")

    samples: list[dict[str, object]] = []
    output_root = args.output_dir or (args.fasta.parent / f"{args.fasta.stem}_rega_like")
    output_root.mkdir(parents=True, exist_ok=True)
    for query_id, query_sequence in query_records:
        aligned_query = aligned[query_id]
        pure_scores = group_scores_for_query(aligned_query, reference_aligned, references, scope="pure")
        all_scores = group_scores_for_query(aligned_query, reference_aligned, references, scope="all")
        pure_windows = assign_windows(aligned_query, reference_aligned, references, scope="pure")
        all_windows = assign_windows(aligned_query, reference_aligned, references, scope="all")
        pure_profiles = build_window_profiles(aligned_query, reference_aligned, references, scope="pure")
        all_profiles = build_window_profiles(aligned_query, reference_aligned, references, scope="all")
        pure_tree = tree_placement_from_tree(pure_tree_global, query_id, references, scope="pure")
        all_tree = tree_placement_from_tree(overall_tree_global, query_id, references, scope="all")
        assignment_label, predicted_group, notes = classify_assignment(
            sequence_length=len(query_sequence.replace("N", "").replace("-", "")),
            pure_tree=pure_tree,
            all_tree=all_tree,
            pure_scores=pure_scores,
            all_scores=all_scores,
            pure_windows=pure_windows,
            all_windows=all_windows,
        )
        parent_summary = summarize_window_parents(all_windows)
        sample_name = query_name_lookup[query_id]
        preferred_plot_groups = []
        if predicted_group:
            preferred_plot_groups.append(predicted_group)
        if pure_tree.cluster_group and pure_tree.cluster_group not in preferred_plot_groups:
            preferred_plot_groups.append(pure_tree.cluster_group)
        pure_plot_groups = top_plot_groups(pure_profiles, preferred=preferred_plot_groups + PURE_GROUP_ORDER, limit=10)
        overall_plot_groups = top_plot_groups(all_profiles, preferred=preferred_plot_groups + PURE_GROUP_ORDER, limit=11)
        bootscan_assets = {
            "pure_csv": "",
            "overall_csv": "",
            "pure_svg": "",
            "overall_svg": "",
            "pure_png": "",
            "overall_png": "",
        }
        if not args.no_assets:
            sample_dir = output_root / sample_name
            sample_dir.mkdir(parents=True, exist_ok=True)
            pure_csv = sample_dir / "bootscan_pure.csv"
            overall_csv = sample_dir / "bootscan_overall.csv"
            pure_svg = sample_dir / "bootscan_pure.svg"
            overall_svg = sample_dir / "bootscan_overall.svg"
            write_bootscan_csv(pure_csv, pure_profiles, pure_plot_groups)
            write_bootscan_csv(overall_csv, all_profiles, overall_plot_groups)
            dominant_support = parent_summary[0][2] if parent_summary else 0.0
            render_bootscan_svg(
                pure_svg,
                pure_profiles,
                pure_plot_groups,
                sample_name=sample_name,
                support_score=dominant_support,
                title_suffix="",
                predicted_label=f"{assignment_label} / {predicted_group or '-'}",
                x_max=len(query_sequence.replace("N", "").replace("-", "")),
            )
            render_bootscan_svg(
                overall_svg,
                all_profiles,
                overall_plot_groups,
                sample_name=sample_name,
                support_score=dominant_support,
                title_suffix="",
                predicted_label=f"{assignment_label} / {predicted_group or '-'}",
                x_max=len(query_sequence.replace("N", "").replace("-", "")),
            )
            pure_png = maybe_render_png_from_svg(pure_svg)
            overall_png = maybe_render_png_from_svg(overall_svg)
            bootscan_assets = {
                "pure_csv": str(pure_csv),
                "overall_csv": str(overall_csv),
                "pure_svg": str(pure_svg),
                "overall_svg": str(overall_svg),
                "pure_png": str(pure_png) if pure_png else "",
                "overall_png": str(overall_png) if overall_png else "",
            }
        samples.append(
            {
                "sample": sample_name,
                "sequence_length": len(query_sequence.replace("N", "").replace("-", "")),
                "assignment_label": assignment_label,
                "predicted_group": predicted_group,
                "predicted_clade": predicted_group,
                "recombination_detected": len(parent_summary) >= 2,
                "candidate_parents": [
                    {"group": group, "supported_windows": count, "fraction": round(fraction, 6)}
                    for group, count, fraction in parent_summary
                ],
                "pure_tree": asdict(pure_tree),
                "overall_tree": asdict(all_tree),
                "pure_top_hits": format_group_scores(pure_scores),
                "overall_top_hits": format_group_scores(all_scores),
                "plot_groups": {
                    "pure": pure_plot_groups,
                    "overall": overall_plot_groups,
                },
                "bootscan_assets": bootscan_assets,
                "pure_window_supports": [asdict(item) for item in pure_windows[:50]],
                "overall_window_supports": [asdict(item) for item in all_windows[:50]],
                "notes": notes + [
                    "Approximation of the REGA/Stanford flow using MAFFT alignment, FastTree placements, and 400/40 window scans, not the official Stanford backend."
                ],
            }
        )

    payload = {
        "tool": "hiv_rega_like",
        "reference_count": len(references),
        "pure_reference_groups": PURE_GROUP_ORDER,
        "bootscan_window_bp": BOOTSCAN_WINDOW,
        "bootscan_step_bp": BOOTSCAN_STEP,
        "samples": samples,
    }

    if args.json:
        print(json.dumps(payload, ensure_ascii=False, indent=2))
        return 0

    for sample in samples:
        print(f"sample\t{sample['sample']}")
        print(f"sequence_length\t{sample['sequence_length']}")
        print(f"assignment_label\t{sample['assignment_label']}")
        print(f"predicted_group\t{sample['predicted_group']}")
        print(f"recombination_detected\t{sample['recombination_detected']}")
        parent_text = ",".join(
            f"{row['group']}:{row['supported_windows']}({row['fraction']:.3f})"
            for row in sample["candidate_parents"]
        ) or "-"
        print(f"candidate_parents\t{parent_text}")
        if sample["overall_top_hits"]:
            top = sample["overall_top_hits"][0]
            print(f"top_hit\t{top['group']} {top['max_identity']:.4f} {top['best_reference']}")
        print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
