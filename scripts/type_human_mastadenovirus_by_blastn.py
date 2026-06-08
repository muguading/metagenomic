from __future__ import annotations

import argparse
import csv
import re
import subprocess
from collections import defaultdict
from pathlib import Path


GENE_CONFIGS = {
    "fiber": Path("database/virus/hadv/blastn_db_fiber/hadv_types_ref_fiber.fa"),
    "hexon": Path("database/virus/hadv/blastn_db_hexon/hadv_types_ref_hexon.fa"),
    "penton": Path("database/virus/hadv/blastn_db_penton/hadv_types_ref_penton.fa"),
}

OUTFMT_FIELDS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
]


def _read_fasta(path: Path) -> dict[str, tuple[str, str]]:
    records: dict[str, tuple[str, str]] = {}
    header = ""
    seq_lines: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header:
                    seq_id = header.split()[0]
                    records[seq_id] = (header, "".join(seq_lines))
                header = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line.strip())
    if header:
        seq_id = header.split()[0]
        records[seq_id] = (header, "".join(seq_lines))
    return records


def _read_manifest(path: Path) -> dict[str, dict[str, str]]:
    mapping: dict[str, dict[str, str]] = {}
    if not path.exists():
        return mapping
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            accession = row.get("accession", "").strip()
            if accession:
                mapping[accession] = row
    return mapping


def _extract_hadv_type(subject_id: str) -> str:
    matched = re.search(r"(HAdV-[A-Z]\d+)", subject_id, flags=re.IGNORECASE)
    if matched:
        return matched.group(1).upper()
    return ""


def _run_blastn(blastn_bin: Path, query_fasta: Path, db_fasta: Path, out_path: Path) -> None:
    cmd = [
        str(blastn_bin),
        "-task",
        "blastn",
        "-query",
        str(query_fasta),
        "-db",
        str(db_fasta),
        "-outfmt",
        "6 " + " ".join(OUTFMT_FIELDS),
        "-evalue",
        "1e-20",
        "-max_target_seqs",
        "10",
        "-dust",
        "no",
        "-num_threads",
        "8",
        "-out",
        str(out_path),
    ]
    subprocess.run(cmd, check=True)


def _pick_best_hits(blast_tsv: Path) -> dict[str, dict[str, str]]:
    best_hits: dict[str, dict[str, str]] = {}
    with blast_tsv.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            values = line.split("\t")
            row = dict(zip(OUTFMT_FIELDS, values))
            qseqid = row["qseqid"]
            candidate = {
                "subject_id": row["sseqid"],
                "matched_type": _extract_hadv_type(row["sseqid"]),
                "pident": row["pident"],
                "align_length": row["length"],
                "mismatch": row["mismatch"],
                "gapopen": row["gapopen"],
                "qstart": row["qstart"],
                "qend": row["qend"],
                "sstart": row["sstart"],
                "send": row["send"],
                "evalue": row["evalue"],
                "bitscore": row["bitscore"],
            }
            previous = best_hits.get(qseqid)
            if previous is None:
                best_hits[qseqid] = candidate
                continue
            previous_score = (float(previous["bitscore"]), float(previous["pident"]), int(previous["align_length"]))
            candidate_score = (float(candidate["bitscore"]), float(candidate["pident"]), int(candidate["align_length"]))
            if candidate_score > previous_score:
                best_hits[qseqid] = candidate
    return best_hits


def _write_group_fastas(
    out_dir: Path,
    gene_name: str,
    records: dict[str, tuple[str, str]],
    best_hits: dict[str, dict[str, str]],
) -> None:
    gene_dir = out_dir / f"by_{gene_name}"
    gene_dir.mkdir(parents=True, exist_ok=True)
    grouped: dict[str, list[str]] = defaultdict(list)
    for query_id, (header, sequence) in records.items():
        type_label = best_hits.get(query_id, {}).get("matched_type") or "UNASSIGNED"
        grouped[type_label].append(f">{header}\n{_wrap_sequence(sequence)}\n")
    for type_label, entries in grouped.items():
        safe_label = re.sub(r"[^A-Za-z0-9._-]+", "_", type_label)
        out_path = gene_dir / f"{safe_label}.{gene_name}.fasta"
        out_path.write_text("".join(entries), encoding="utf-8")


def _wrap_sequence(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[index:index + width] for index in range(0, len(sequence), width))


def main() -> int:
    parser = argparse.ArgumentParser(description="Type Human mastadenovirus full genomes with fiber/hexon/penton blastn databases")
    parser.add_argument(
        "--query-fasta",
        type=Path,
        default=Path("database/virus/hadv/full_genomes/human_mastadenovirus_A_G_complete_genomes.fasta"),
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("database/virus/hadv/full_genomes/human_mastadenovirus_A_G_manifest.tsv"),
    )
    parser.add_argument(
        "--blastn-bin",
        type=Path,
        default=Path("/opt/homebrew/Caskroom/mambaforge/base/envs/ncov/bin/blastn"),
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("database/virus/hadv/full_genomes/typing_by_blastn"),
    )
    args = parser.parse_args()

    query_fasta = args.query_fasta.resolve()
    manifest_path = args.manifest.resolve()
    blastn_bin = args.blastn_bin.resolve()
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    records = _read_fasta(query_fasta)
    manifest = _read_manifest(manifest_path)
    gene_best_hits: dict[str, dict[str, dict[str, str]]] = {}

    for gene_name, db_path in GENE_CONFIGS.items():
        db_fasta = db_path.resolve()
        blast_out = out_dir / f"{gene_name}.blastn.tsv"
        _run_blastn(blastn_bin=blastn_bin, query_fasta=query_fasta, db_fasta=db_fasta, out_path=blast_out)
        best_hits = _pick_best_hits(blast_out)
        gene_best_hits[gene_name] = best_hits

        detail_tsv = out_dir / f"{gene_name}.typing.tsv"
        columns = [
            "query_id",
            "species_group",
            "species_name",
            "sequence_length",
            "matched_type",
            "subject_id",
            "pident",
            "align_length",
            "qstart",
            "qend",
            "sstart",
            "send",
            "evalue",
            "bitscore",
        ]
        with detail_tsv.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
            writer.writeheader()
            for query_id in sorted(records):
                row = manifest.get(query_id, {})
                hit = best_hits.get(query_id, {})
                writer.writerow(
                    {
                        "query_id": query_id,
                        "species_group": row.get("species_group", ""),
                        "species_name": row.get("species_name", ""),
                        "sequence_length": row.get("sequence_length", str(len(records[query_id][1]))),
                        "matched_type": hit.get("matched_type", ""),
                        "subject_id": hit.get("subject_id", ""),
                        "pident": hit.get("pident", ""),
                        "align_length": hit.get("align_length", ""),
                        "qstart": hit.get("qstart", ""),
                        "qend": hit.get("qend", ""),
                        "sstart": hit.get("sstart", ""),
                        "send": hit.get("send", ""),
                        "evalue": hit.get("evalue", ""),
                        "bitscore": hit.get("bitscore", ""),
                    }
                )
        _write_group_fastas(out_dir=out_dir, gene_name=gene_name, records=records, best_hits=best_hits)

    summary_columns = [
        "query_id",
        "species_group",
        "species_name",
        "sequence_length",
        "fiber_type",
        "hexon_type",
        "penton_type",
        "fiber_subject",
        "hexon_subject",
        "penton_subject",
    ]
    summary_path = out_dir / "fiber_hexon_penton.typing.summary.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_columns, delimiter="\t")
        writer.writeheader()
        for query_id in sorted(records):
            row = manifest.get(query_id, {})
            writer.writerow(
                {
                    "query_id": query_id,
                    "species_group": row.get("species_group", ""),
                    "species_name": row.get("species_name", ""),
                    "sequence_length": row.get("sequence_length", str(len(records[query_id][1]))),
                    "fiber_type": gene_best_hits["fiber"].get(query_id, {}).get("matched_type", ""),
                    "hexon_type": gene_best_hits["hexon"].get(query_id, {}).get("matched_type", ""),
                    "penton_type": gene_best_hits["penton"].get(query_id, {}).get("matched_type", ""),
                    "fiber_subject": gene_best_hits["fiber"].get(query_id, {}).get("subject_id", ""),
                    "hexon_subject": gene_best_hits["hexon"].get(query_id, {}).get("subject_id", ""),
                    "penton_subject": gene_best_hits["penton"].get(query_id, {}).get("subject_id", ""),
                }
            )

    print(summary_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
