from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "database/virus/orthohantavirus/reference_genomes/orthohantavirus_typing_reference_segments_manifest.tsv"
GFF_DIR = ROOT / "database/virus/orthohantavirus/reference_genomes/gff3_normalized"
OUTPUT = ROOT / "tmp/orthohantavirus_structure_comparison.png"

VIRUS_ORDER = [
    ("HTNV", "Hantaan virus", "Orthohantavirus hantanense", {"L": "X55901", "M": "M14627", "S": "M14626"}),
    ("SEOV", "Seoul virus", "Orthohantavirus seoulense", {"L": "X56492", "M": "S47716", "S": "AY273791"}),
    ("PUUV", "Puumala virus", "Orthohantavirus puumalaense", {"L": "MN832782", "M": "MN832783", "S": "MN832784"}),
    ("SNV", "Sin Nombre virus", "Orthohantavirus sinnombreense", {"L": "L37902", "M": "L37903", "S": "L37904"}),
    ("ANDV", "Andes virus", "Orthohantavirus andesense", {"L": "AF291704", "M": "AF291703", "S": "AF291702"}),
]

SEGMENT_ORDER = ("L", "M", "S")
SEGMENT_COLORS = {"L": "#20E3E6", "M": "#F4B2BF", "S": "#F012F6"}
FEATURE_COLORS = {"NSs": "#FFC857"}
SEGMENT_MAX = {"L": 7000, "M": 4000, "S": 2200}
SEGMENT_LABEL = {"L": "L", "M": "M", "S": "S"}


@dataclass
class Feature:
    start: int
    end: int
    label: str


@dataclass
class SegmentRecord:
    virus_name: str
    species: str
    segment: str
    accession: str
    record_id: str
    sequence_length: int
    features: list[Feature]


def load_manifest() -> list[dict[str, str]]:
    with MANIFEST.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def simplify_label(segment: str, product: str, cds_index: int) -> str:
    product_l = product.lower()
    if "non-structural" in product_l or "nss" in product_l:
        return "NSs"
    if segment == "L":
        return "L"
    if segment == "M":
        return "GPC"
    if segment == "S":
        return "N" if cds_index == 0 else "NSs"
    return product


def parse_gff(record_id: str, segment: str) -> list[Feature]:
    gff_path = GFF_DIR / f"{record_id}.gff3"
    features: list[Feature] = []
    cds_index = 0
    with gff_path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9 or parts[2] != "CDS":
                continue
            attrs = {}
            for item in parts[8].split(";"):
                if "=" in item:
                    key, value = item.split("=", 1)
                    attrs[key] = value
            label = simplify_label(segment, attrs.get("product", ""), cds_index)
            features.append(Feature(start=int(parts[3]), end=int(parts[4]), label=label))
            cds_index += 1
    return sorted(features, key=lambda feature: (feature.label == "NSs", feature.start))


def collect_records() -> list[SegmentRecord]:
    manifest = load_manifest()
    records: list[SegmentRecord] = []
    for abbrev, virus_name, species, accessions in VIRUS_ORDER:
        for segment in SEGMENT_ORDER:
            accession = accessions[segment]
            row = next(
                item
                for item in manifest
                if item["abbrev"] == abbrev and item["segment"] == segment and item["accession"] == accession
            )
            records.append(
                SegmentRecord(
                    virus_name=virus_name,
                    species=species,
                    segment=segment,
                    accession=row["accession"],
                    record_id=row["record_id"],
                    sequence_length=int(row["sequence_length"]),
                    features=parse_gff(row["record_id"], segment),
                )
            )
    return records


def add_segment_bar(ax, x0: float, y0: float, width: float, height: float, color: str) -> None:
    wedge = min(22, max(8, width * 0.022))
    points = [
        (x0, y0),
        (x0 + width, y0),
        (x0 + width - wedge, y0 - height),
        (x0, y0 - height),
    ]
    ax.add_patch(Polygon(points, closed=True, facecolor=color, edgecolor="#111111", linewidth=1.15, joinstyle="miter"))


def add_feature_bar(ax, x0: float, y0: float, width: float, height: float, color: str) -> None:
    wedge = min(16, max(6, width * 0.03))
    points = [
        (x0, y0),
        (x0 + width, y0),
        (x0 + width - wedge, y0 - height),
        (x0, y0 - height),
    ]
    ax.add_patch(Polygon(points, closed=True, facecolor=color, edgecolor="#111111", linewidth=0.9))


def draw(records: list[SegmentRecord]) -> None:
    groups: dict[str, dict[str, SegmentRecord]] = {}
    for record in records:
        groups.setdefault(record.virus_name, {})[record.segment] = record

    left_margin = -900
    title_y = 2.0
    row_gap = 1.85
    bar_height = 0.28
    panel_gap = 220

    panel_offsets: dict[str, float] = {}
    cursor = 0.0
    for segment in SEGMENT_ORDER:
        panel_offsets[segment] = cursor
        cursor += SEGMENT_MAX[segment] + panel_gap
    total_width = cursor - panel_gap

    total_height = len(VIRUS_ORDER) * row_gap + 0.6
    fig, ax = plt.subplots(figsize=(17.2, 8.8), dpi=220)
    ax.set_xlim(left_margin, total_width + 180)
    ax.set_ylim(-0.4, total_height + 1.25)
    ax.axis("off")

    ax.text(left_margin + 20, total_height + title_y, "Orthohantavirus Genome Structure Comparison", fontsize=24, fontweight="bold", ha="left", va="top")
    ax.text(
        left_margin + 20,
        total_height + title_y - 0.42,
        "Each virus is shown on a single row; L, M, and S segments are arranged side-by-side using segment-specific scales",
        fontsize=11.5,
        color="#555555",
        ha="left",
        va="top",
    )

    for segment in SEGMENT_ORDER:
        x0 = panel_offsets[segment]
        x1 = x0 + SEGMENT_MAX[segment]
        center = (x0 + x1) / 2
        ax.text(center, total_height + 0.72, f"{SEGMENT_LABEL[segment]} segment", fontsize=13.5, fontweight="bold", ha="center", va="bottom")
        ax.text(center, total_height + 0.48, f"0-{SEGMENT_MAX[segment]} nt", fontsize=10.2, color="#666666", ha="center", va="bottom")

    for index, (_abbrev, virus_name, species, _accessions) in enumerate(VIRUS_ORDER):
        y_mid = total_height - index * row_gap
        ax.text(left_margin + 20, y_mid + 0.33, virus_name, fontsize=17.5, fontweight="bold", ha="left", va="center")
        ax.text(left_margin + 20, y_mid + 0.08, species, fontsize=11.5, color="#555555", fontstyle="italic", ha="left", va="center")

        for segment in SEGMENT_ORDER:
            record = groups[virus_name][segment]
            panel_x = panel_offsets[segment]
            y_top = y_mid - 0.15
            add_segment_bar(ax, panel_x, y_top, record.sequence_length, bar_height, SEGMENT_COLORS[segment])
            ax.text(panel_x - 7, y_top - bar_height / 2, "3'", fontsize=9.7, ha="right", va="center")
            ax.text(panel_x + record.sequence_length + 20, y_top - bar_height / 2, "5'", fontsize=9.7, ha="left", va="center")
            ax.text(panel_x + record.sequence_length - 5, y_top + 0.11, str(record.sequence_length), fontsize=9.6, ha="center", va="bottom", color="#222222")
            ax.text(panel_x, y_top - bar_height - 0.13, record.accession, fontsize=10.6, ha="left", va="top", color="#333333")

            for feature in record.features:
                feature_width = feature.end - feature.start + 1
                feature_height = bar_height * (0.62 if feature.label == "NSs" else 1.0)
                feature_y = y_top - (bar_height - feature_height) / 2
                if feature.label == "NSs":
                    feature_y = y_top - feature_height - 0.01
                add_feature_bar(
                    ax,
                    panel_x + feature.start - 1,
                    feature_y,
                    feature_width,
                    feature_height,
                    FEATURE_COLORS.get(feature.label, SEGMENT_COLORS[segment]),
                )
                if feature_width >= 120:
                    ax.text(
                        panel_x + feature.start - 1 + feature_width / 2,
                        feature_y - feature_height / 2,
                        feature.label,
                        fontsize=8.8 if feature.label != "NSs" else 7.8,
                        fontstyle="italic" if feature.label in {"L", "GPC", "N"} else "normal",
                        ha="center",
                        va="center",
                    )
                elif feature.label == "NSs":
                    ax.text(panel_x + feature.end + 14, feature_y - feature_height / 2, "NSs", fontsize=7.8, ha="left", va="center", color="#7A5200")

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    draw(collect_records())
