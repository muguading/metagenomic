from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path


DEFAULT_CSV = Path("database/virus/bandavirus/bandavirus_db.csv")
DEFAULT_FASTA = Path("database/virus/bandavirus/bandavirus_db.fasta")
DEFAULT_OUT_CSV = Path("database/virus/bandavirus/bandavirus_db_grouped.csv")
DEFAULT_OUT_FASTA = Path("database/virus/bandavirus/bandavirus_db_grouped.fasta")
DEFAULT_OUT_SUMMARY = Path("database/virus/bandavirus/bandavirus_db_grouped_summary.tsv")
DEFAULT_OUT_INCOMPLETE = Path("database/virus/bandavirus/bandavirus_db_incomplete_samples.tsv")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Group bandavirus records by sample and assign L/M/S segments based on segment annotation and genome length."
    )
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV, help="Input metadata CSV")
    parser.add_argument("--fasta", type=Path, default=DEFAULT_FASTA, help="Input FASTA")
    parser.add_argument("--out-csv", type=Path, default=DEFAULT_OUT_CSV, help="Output grouped CSV")
    parser.add_argument("--out-fasta", type=Path, default=DEFAULT_OUT_FASTA, help="Output grouped FASTA")
    parser.add_argument("--summary-tsv", type=Path, default=DEFAULT_OUT_SUMMARY, help="Output summary TSV")
    parser.add_argument(
        "--incomplete-tsv",
        type=Path,
        default=DEFAULT_OUT_INCOMPLETE,
        help="Output incomplete-sample report TSV",
    )
    return parser.parse_args()


def normalize_text(value: str) -> str:
    text = re.sub(r"\s+", " ", str(value or "").strip())
    return text


def extract_title_sample_label(title: str) -> tuple[str, str]:
    text = normalize_text(title)
    if not text:
        return "", ""
    stop_terms = (
        r"(?:\s+segment\b|\s+glycoprotein gene\b|\s+glycoprotein precursor\b|\s+glycoprotein\b|"
        r"\s+polyprotein gene\b|\s+polymerase gene\b|\s+RNA-dependent RNA polymerase\b|"
        r"\s+nucleocapsid protein\b|\s+nonstructural protein\b|\s+complete genome\b|"
        r"\s+complete sequence\b|\s+complete cds\b)"
    )
    patterns = [
        rf"\bstrain\s+(.+?)(?:{stop_terms}|,|$)",
        rf"\bisolate\s+(.+?)(?:{stop_terms}|,|$)",
    ]
    for pattern in patterns:
        match = re.search(pattern, text, flags=re.IGNORECASE)
        if match:
            return normalize_text(match.group(1)), "specific"
    text = re.sub(
        rf"{stop_terms}.*$",
        "",
        text,
        flags=re.IGNORECASE,
    )
    text = re.sub(r",.*$", "", text).strip()
    return normalize_text(text), "stem"


def sample_key(row: dict[str, str]) -> tuple[str, str]:
    title_label, title_label_mode = extract_title_sample_label(row.get("GenBank_Title", ""))
    candidates = [
        ("organism_isolate_biosample_assembly", ["Organism_Name", "Isolate", "BioSample", "Assembly"]),
        ("organism_isolate_assembly", ["Organism_Name", "Isolate", "Assembly"]),
        ("organism_biosample_assembly", ["Organism_Name", "BioSample", "Assembly"]),
        ("organism_isolate", ["Organism_Name", "Isolate"]),
        ("organism_title_label_country_host_date", ["Organism_Name", "__TITLE_LABEL__", "Country", "Host", "Collection_Date"]),
        ("organism_title_label_host_date", ["Organism_Name", "__TITLE_LABEL__", "Host", "Collection_Date"]),
        ("species_isolate", ["Species", "Isolate"]),
        ("organism_biosample", ["Organism_Name", "BioSample"]),
        ("organism_assembly", ["Organism_Name", "Assembly"]),
    ]
    for label, fields in candidates:
        if "__TITLE_LABEL__" in fields and not title_label:
            continue
        if label == "organism_title_label_host_date" and title_label_mode != "specific":
            continue
        values = []
        for field in fields:
            if field == "__TITLE_LABEL__":
                values.append(title_label)
            else:
                values.append(normalize_text(row.get(field, "")))
        if all(values):
            return "|".join(values), label
    fallback = [normalize_text(row.get("Organism_Name", "")), normalize_text(row.get("Accession", ""))]
    return "|".join(fallback), "organism_accession"


def read_fasta_map(path: Path) -> dict[str, tuple[str, str]]:
    fasta_map: dict[str, tuple[str, str]] = {}
    header = ""
    seq_lines: list[str] = []
    accession = ""
    for raw_line in path.open(encoding="utf-8", errors="ignore"):
        line = raw_line.rstrip("\n")
        if line.startswith(">"):
            if header:
                fasta_map[accession] = (header, "".join(seq_lines).upper())
            header = line[1:].strip()
            accession = header.split()[0].split("|")[0].strip()
            accession = accession.split(".", 1)[0]
            seq_lines = []
        else:
            seq_lines.append(line.strip())
    if header:
        fasta_map[accession] = (header, "".join(seq_lines).upper())
    return fasta_map


def numeric_length(row: dict[str, str]) -> int:
    try:
        return int(str(row.get("Length", "")).strip())
    except ValueError:
        return 0


def segment_priority(segment: str) -> int:
    return {"L": 0, "M": 1, "S": 2}.get(segment, 9)


def segment_from_length(length: int) -> str:
    if length >= 5000:
        return "L"
    if length >= 2500:
        return "M"
    return "S"


def completeness_score(value: str) -> int:
    text = normalize_text(value).lower()
    if text == "complete":
        return 2
    if text:
        return 1
    return 0


def source_score(value: str) -> int:
    text = normalize_text(value)
    if text == "RefSeq":
        return 2
    if text == "GenBank":
        return 1
    return 0


def record_rank(row: dict[str, str]) -> tuple[int, int, int, str]:
    return (
        completeness_score(row.get("Nuc_Completeness", "")),
        source_score(row.get("GenBank_RefSeq", "")),
        numeric_length(row),
        normalize_text(row.get("Accession", "")),
    )


def assign_segments(group_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    rows = [dict(row) for row in group_rows]
    for row in rows:
        original_segment = normalize_text(row.get("Segment", "")).upper()
        row["Original_Segment"] = original_segment
        row["Length_Int"] = str(numeric_length(row))
        if original_segment in {"L", "M", "S"}:
            row["Assigned_Segment"] = original_segment
            row["Segment_Source"] = "original"
        else:
            row["Assigned_Segment"] = segment_from_length(numeric_length(row))
            row["Segment_Source"] = "length"

    best_by_segment: dict[str, dict[str, str]] = {}
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["Assigned_Segment"]].append(row)

    for segment in ("L", "M", "S"):
        candidates = grouped.get(segment, [])
        if not candidates:
            continue
        candidates.sort(key=record_rank, reverse=True)
        best = candidates[0]
        for candidate in candidates:
            candidate["Selected_For_Sample"] = "yes" if candidate is best else "no"
        best_by_segment[segment] = best

    selected_rows: list[dict[str, str]] = []
    for segment in ("L", "M", "S"):
        row = best_by_segment.get(segment)
        if row:
            selected_rows.append(row)
    return selected_rows


def build_fasta_header(row: dict[str, str]) -> str:
    parts = [
        f"sample={row['Sample_Key']}",
        f"segment={row['Assigned_Segment']}",
        f"accession={row['Accession']}",
        f"length={row['Length']}",
        f"organism={normalize_text(row.get('Organism_Name', ''))}",
    ]
    isolate = normalize_text(row.get("Isolate", ""))
    if isolate:
        parts.append(f"isolate={isolate}")
    return " | ".join(parts)


def wrap_sequence(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[i:i + width] for i in range(0, len(sequence), width))


def main() -> int:
    args = parse_args()
    csv_path = args.csv.expanduser().resolve()
    fasta_path = args.fasta.expanduser().resolve()
    out_csv = args.out_csv.expanduser().resolve()
    out_fasta = args.out_fasta.expanduser().resolve()
    summary_tsv = args.summary_tsv.expanduser().resolve()
    incomplete_tsv = args.incomplete_tsv.expanduser().resolve()

    out_csv.parent.mkdir(parents=True, exist_ok=True)

    fasta_map = read_fasta_map(fasta_path)
    with csv_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))

    usable_rows: list[dict[str, str]] = []
    for row in rows:
        accession = normalize_text(row.get("Accession", "")).split(".", 1)[0]
        if accession not in fasta_map:
            continue
        row = dict(row)
        row["Accession"] = accession
        key, key_source = sample_key(row)
        row["Sample_Key"] = key
        row["Sample_Key_Source"] = key_source
        usable_rows.append(row)

    sample_groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in usable_rows:
        sample_groups[row["Sample_Key"]].append(row)

    grouped_rows: list[dict[str, str]] = []
    summary_rows: list[dict[str, str]] = []
    incomplete_rows: list[dict[str, str]] = []
    fasta_chunks: list[str] = []

    for key, group in sorted(sample_groups.items()):
        selected = assign_segments(group)
        selected_segments = [row["Assigned_Segment"] for row in selected]
        status = "complete" if selected_segments == ["L", "M", "S"] else "incomplete"

        summary_rows.append(
            {
                "sample_key": key,
                "sample_key_source": group[0]["Sample_Key_Source"],
                "input_records": str(len(group)),
                "selected_records": str(len(selected)),
                "selected_segments": ",".join(selected_segments),
                "status": status,
            }
        )

        if status != "complete":
            incomplete_rows.append(
                {
                    "sample_key": key,
                    "sample_key_source": group[0]["Sample_Key_Source"],
                    "input_records": str(len(group)),
                    "selected_segments": ",".join(selected_segments),
                    "accessions": ",".join(normalize_text(row.get("Accession", "")) for row in group),
                }
            )
            continue

        selected.sort(key=lambda row: segment_priority(row["Assigned_Segment"]))
        for row in selected:
            header, sequence = fasta_map[row["Accession"]]
            row["Original_Header"] = header
            row["Sequence_Length_From_FASTA"] = str(len(sequence))
            row["Group_Status"] = status
            grouped_rows.append(row)
            fasta_chunks.append(f">{build_fasta_header(row)}\n{wrap_sequence(sequence)}")

    base_columns = list(rows[0].keys()) if rows else []
    extra_columns = [
        "Sample_Key",
        "Sample_Key_Source",
        "Original_Segment",
        "Assigned_Segment",
        "Segment_Source",
        "Original_Header",
        "Sequence_Length_From_FASTA",
        "Group_Status",
    ]
    out_columns = base_columns + [column for column in extra_columns if column not in base_columns]

    with out_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=out_columns)
        writer.writeheader()
        for row in grouped_rows:
            writer.writerow({column: row.get(column, "") for column in out_columns})

    out_fasta.write_text("\n".join(fasta_chunks).strip() + "\n", encoding="utf-8")

    with summary_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["sample_key", "sample_key_source", "input_records", "selected_records", "selected_segments", "status"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    with incomplete_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["sample_key", "sample_key_source", "input_records", "selected_segments", "accessions"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(incomplete_rows)

    counts = Counter(row["Assigned_Segment"] for row in grouped_rows)
    print(f"complete_samples\t{len(grouped_rows) // 3}")
    print(f"grouped_records\t{len(grouped_rows)}")
    print(f"segment_counts\tL={counts.get('L', 0)};M={counts.get('M', 0)};S={counts.get('S', 0)}")
    print(f"out_csv\t{out_csv}")
    print(f"out_fasta\t{out_fasta}")
    print(f"summary_tsv\t{summary_tsv}")
    print(f"incomplete_tsv\t{incomplete_tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
