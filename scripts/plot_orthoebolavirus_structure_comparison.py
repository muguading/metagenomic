from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon


ROOT = Path(__file__).resolve().parents[1]
REF_DIR = ROOT / "database/virus/Orthoebolavirus/reference_genomes"
MANIFEST = REF_DIR / "orthoebolavirus_typing_reference_genomes_manifest.tsv"
GFF_DIR = REF_DIR / "gff3_normalized"
OUTPUT = ROOT / "tmp/orthoebolavirus_structure_comparison.png"

PLOT_ORDER = [
    ("EBOV", "Zaire"),
    ("SUDV", "Sudan"),
    ("TAFV", "Taï Forest"),
    ("BDBV", "Bundibugyo"),
    ("RESTV", "Reston"),
]

GENE_ORDER = ["NP", "VP35", "VP40", "GP", "VP30", "VP24", "L"]
COLORS = {
    "NP": "#F012F6",
    "VP35": "#FF1717",
    "VP40": "#F3FF1B",
    "GP": "#F5B5C4",
    "VP30": "#F8E8C8",
    "VP24": "#FFA600",
    "L": "#22E3E6",
}


@dataclass
class Feature:
    start: int
    end: int
    label: str
    lane: int = 0


@dataclass
class VirusRecord:
    abbrev: str
    accession: str
    species: str
    virus_name: str
    record_id: str
    sequence_length: int
    features: list[Feature]


def load_manifest() -> dict[str, dict[str, str]]:
    with MANIFEST.open(newline="", encoding="utf-8") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        return {row["abbrev"]: row for row in rows}


def parse_attrs(attr_text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in attr_text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = unquote(value)
    return attrs


def parse_features(record_id: str) -> list[Feature]:
    features_by_gene: dict[str, list[tuple[int, int, str, str]]] = {}
    gp_edit_sites: list[int] = []
    with (GFF_DIR / f"{record_id}.gff3").open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) != 9:
                continue
            attrs = parse_attrs(columns[8])
            gene = attrs.get("gene", "")
            note = attrs.get("Note", "").lower()
            if columns[2] != "CDS":
                if gene in {"GP", "sGP"} and any(token in note for token in ("editing", "insertion", "polymerase slippage")):
                    gp_edit_sites.append(int(columns[3]))
                continue
            if gene not in GENE_ORDER and gene != "sGP":
                continue
            start = int(columns[3])
            end = int(columns[4])
            feature_id = attrs.get("ID", "")
            product = attrs.get("product", "")
            features_by_gene.setdefault(gene, []).append((start, end, feature_id, product))

    features: list[Feature] = []
    for gene in GENE_ORDER:
        parts = features_by_gene.get(gene, [])
        if not parts:
            continue
        if gene == "GP":
            structural = [
                part
                for part in parts
                if "nonstructural" not in part[3].lower()
                and any(token in part[3].lower() for token in ("spike", "virion", "structural glycoprotein"))
            ]
            structural = sorted(structural or parts, key=lambda part: (part[0], part[1]))
            if len(structural) == 1 and gp_edit_sites:
                start, end, _feature_id, _product = structural[0]
                split = min(site for site in gp_edit_sites if start < site < end)
                features.append(Feature(start, split, "GP", 0))
                features.append(Feature(split, end, "GP", 1))
            else:
                for index, part in enumerate(structural):
                    features.append(Feature(part[0], part[1], "GP", min(index, 1)))
            continue
        starts = [part[0] for part in parts]
        ends = [part[1] for part in parts]
        features.append(Feature(min(starts), max(ends), gene, 0))

    return features


def collect_records() -> list[VirusRecord]:
    manifest = load_manifest()
    records: list[VirusRecord] = []
    for abbrev, _label_en in PLOT_ORDER:
        row = manifest[abbrev]
        records.append(
            VirusRecord(
                abbrev=abbrev,
                accession=row["accession"],
                species=row["species"],
                virus_name=row["virus_name"],
                record_id=row["record_id"],
                sequence_length=int(row["sequence_length"]),
                features=parse_features(row["record_id"]),
            )
        )
    return records


def add_bar(ax, x0: float, y0: float, width: float, height: float, color: str, lw: float = 1.0) -> None:
    wedge = min(95, max(28, width * 0.08))
    points = [
        (x0, y0),
        (x0 + width, y0),
        (x0 + width + wedge, y0 - height / 2),
        (x0 + width, y0 - height),
        (x0, y0 - height),
    ]
    ax.add_patch(Polygon(points, closed=True, facecolor=color, edgecolor="#111111", linewidth=lw, joinstyle="miter"))


def draw(records: list[VirusRecord]) -> None:
    max_len = 20000
    left_margin = -1080
    row_gap = 2.0
    base_y = len(records) * row_gap + 0.2
    bar_height = 0.34

    fig, ax = plt.subplots(figsize=(16, 10), dpi=220)
    ax.set_xlim(left_margin, max_len + 450)
    ax.set_ylim(-0.8, base_y + 1.7)
    ax.axis("off")

    ax.text(left_margin + 25, base_y + 1.35, "genus Orthoebolavirus", fontsize=13.5, fontweight="bold", ha="left", va="center")
    ax.text(left_margin + 25, base_y + 1.0, "Genome structure comparison", fontsize=20, fontweight="bold", ha="left", va="center")

    for row_index, record in enumerate(records):
        _abbrev, label_en = PLOT_ORDER[row_index]
        y = base_y - row_index * row_gap
        ax.text(left_margin + 25, y + 0.56, label_en, fontsize=16, fontweight="bold", ha="left", va="center")
        ax.text(left_margin + 25, y + 0.26, f"species {record.species}", fontsize=12.5, fontweight="bold", fontstyle="italic", ha="left", va="center")
        ax.text(left_margin + 430, y - 0.05, f"{record.accession} {record.virus_name}", fontsize=12.2, ha="left", va="center")
        ax.text(record.sequence_length + 90, y - 0.42, f"{record.sequence_length:,}".replace(",", " "), fontsize=10.5, ha="center", va="bottom")

        axis_y = y - 0.52
        ax.plot([1, record.sequence_length], [axis_y, axis_y], color="#111111", lw=1.4)
        ax.text(-70, axis_y, "3'", fontsize=10.8, ha="right", va="center")
        ax.text(record.sequence_length + 120, axis_y, "5'", fontsize=10.8, ha="left", va="center")

        for feature in sorted(record.features, key=lambda item: (item.lane, item.start)):
            feature_y = axis_y - 0.02
            if feature.lane == 1:
                feature_y = axis_y - 0.52
                ax.plot([feature.start, feature.start], [axis_y - 0.05, feature_y + 0.02], color="#111111", lw=0.95)
            add_bar(ax, feature.start, feature_y, feature.end - feature.start + 1, bar_height, COLORS[feature.label])
            ax.text(
                feature.start + (feature.end - feature.start + 1) / 2,
                feature_y - bar_height / 2,
                feature.label,
                fontsize=9.8,
                fontstyle="italic",
                ha="center",
                va="center",
            )

    axis_y = -0.08
    ax.plot([1, max_len], [axis_y, axis_y], color="#999999", lw=0.9)
    for pos in range(0, max_len + 1, 200):
        height = 0.18
        if pos % 2000 == 0:
            height = 0.52
        elif pos % 1000 == 0:
            height = 0.34
        ax.plot([pos, pos], [axis_y, axis_y + height], color="#999999", lw=0.75)
    for pos in range(0, max_len + 1, 2000):
        label = "1" if pos == 0 else f"{pos:,}".replace(",", " ")
        ax.text(max(pos, 1), axis_y - 0.25, label, fontsize=9.5, ha="center", va="top", color="#333333")

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    draw(collect_records())
