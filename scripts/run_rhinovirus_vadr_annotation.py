#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import os
import shlex
import subprocess
from pathlib import Path


ROOT = Path("/Users/wuhhh/Desktop/徐老师/代码/metagenomic")
SOFT_DIR = ROOT / "soft"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Annotate rhinovirus genomes with VADR using hrvA/hrvB/hrvC models."
    )
    parser.add_argument(
        "--input-fasta",
        type=Path,
        help="Input FASTA file. Required unless --manifest is provided.",
    )
    parser.add_argument(
        "--group",
        choices=["A", "B", "C"],
        help="Rhinovirus species group for --input-fasta mode.",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        help="Optional TSV manifest with at least accession/species_group columns for batch mode.",
    )
    parser.add_argument(
        "--sequence-fasta",
        type=Path,
        help="Optional FASTA corresponding to --manifest. If omitted in batch mode, --input-fasta is used.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        required=True,
        help="Output directory.",
    )
    parser.add_argument(
        "--prefix",
        default="rhinovirus",
        help="Output prefix in single-group mode. Default: rhinovirus",
    )
    parser.add_argument(
        "--keep-existing",
        action="store_true",
        help="Reuse existing output directories if they already contain VADR tables.",
    )
    parser.add_argument(
        "--split",
        action="store_true",
        help="Pass --split to v-annotate.pl. Enabled automatically when --cpu > 1.",
    )
    parser.add_argument(
        "--cpu",
        type=int,
        default=1,
        help="CPU workers for v-annotate.pl. Default: 1",
    )
    return parser.parse_args()


def parse_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    header: str | None = None
    chunks: list[str] = []
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records[header] = "".join(chunks)
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        records[header] = "".join(chunks)
    return records


def build_vadr_env(model_dir: Path) -> dict[str, str]:
    vadr_root = SOFT_DIR.resolve()
    vadr_scripts_dir = (vadr_root / "vadr").resolve()
    infernal_bin_dir = (vadr_root / "infernal" / "binaries").resolve()
    bio_easel_dir = (vadr_root / "Bio-Easel").resolve()
    sequip_dir = (vadr_root / "sequip").resolve()
    blast_bin_dir = (vadr_root / "ncbi-blast" / "bin").resolve()
    fasta_bin_dir = (vadr_root / "fasta" / "bin").resolve()
    minimap2_dir = (vadr_root / "minimap2").resolve()

    required_paths = [
        vadr_scripts_dir / "v-annotate.pl",
        vadr_scripts_dir / "miniscripts" / "annotate-tbl2gff.pl",
        vadr_scripts_dir / "miniscripts" / "fasta-trim-terminal-ambigs.pl",
        infernal_bin_dir / "cmalign",
        bio_easel_dir / "blib" / "lib",
        bio_easel_dir / "blib" / "arch",
        sequip_dir,
        blast_bin_dir / "blastn",
        fasta_bin_dir / "fasta36",
        minimap2_dir / "minimap2",
        model_dir,
    ]
    missing = [str(path) for path in required_paths if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing VADR dependency paths: " + ", ".join(missing[:10]))

    inherited_path = str(os.environ.get("PATH") or "").strip()
    inherited_perl5lib = str(os.environ.get("PERL5LIB") or "").strip()
    return {
        **os.environ,
        "VADRINSTALLDIR": str(vadr_root),
        "VADRSCRIPTSDIR": str(vadr_scripts_dir),
        "VADRCONFIGFILE": str((vadr_scripts_dir / "vadr.config").resolve()),
        "VADRMODELDIR": str(model_dir.resolve()),
        "VADRINFERNALDIR": str(infernal_bin_dir),
        "VADREASELDIR": str(infernal_bin_dir),
        "VADRHMMERDIR": str(infernal_bin_dir),
        "VADRBIOEASELDIR": str(bio_easel_dir),
        "VADRSEQUIPDIR": str(sequip_dir),
        "VADRBLASTDIR": str(blast_bin_dir),
        "VADRFASTADIR": str(fasta_bin_dir),
        "VADRMINIMAP2DIR": str(minimap2_dir),
        "PERL5LIB": os.pathsep.join(
            [
                str(vadr_scripts_dir),
                str(sequip_dir),
                str((bio_easel_dir / "blib" / "lib").resolve()),
                str((bio_easel_dir / "blib" / "arch").resolve()),
            ] + ([inherited_perl5lib] if inherited_perl5lib else [])
        ),
        "PATH": os.pathsep.join(
            [
                str(vadr_scripts_dir),
                str(blast_bin_dir),
                str(fasta_bin_dir),
                str(infernal_bin_dir),
                str(minimap2_dir),
            ] + ([inherited_path] if inherited_path else [])
        ),
    }


def run_command(cmd: str, env: dict[str, str]) -> None:
    subprocess.run(cmd, shell=True, check=True, env=env)


def tbl_has_feature_rows(tbl_path: Path) -> bool:
    if not tbl_path.exists():
        return False
    with tbl_path.open() as handle:
        return any(line.startswith(">Feature ") for line in handle)


def split_gff_by_seqid(group_gff: Path, output_dir: Path) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    buckets: dict[str, list[str]] = {}
    headers: list[str] = []
    with group_gff.open() as handle:
        for raw_line in handle:
            if raw_line.startswith("#"):
                headers.append(raw_line)
                continue
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            seqid = fields[0]
            if seqid not in buckets:
                buckets[seqid] = list(headers)
            buckets[seqid].append(raw_line)

    output_paths: list[Path] = []
    for seqid, lines in buckets.items():
        out_path = output_dir / f"{seqid}.gff3"
        with out_path.open("w") as handle:
            handle.writelines(lines)
        output_paths.append(out_path)
    return sorted(output_paths)


def build_group_fasta_from_manifest(manifest_path: Path, sequence_fasta_path: Path, out_dir: Path) -> dict[str, Path]:
    records = parse_fasta(sequence_fasta_path)
    with manifest_path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    group_rows: dict[str, list[dict[str, str]]] = {"A": [], "B": [], "C": []}
    for row in rows:
        group = str(row.get("species_group") or "").strip()
        if group in group_rows:
            group_rows[group].append(row)

    outputs: dict[str, Path] = {}
    for group, items in group_rows.items():
        if not items:
            continue
        fasta_path = out_dir / f"hrv{group}.input.fasta"
        with fasta_path.open("w") as handle:
            for row in items:
                header = str(row.get("header") or "").strip()
                accession = str(row.get("accession") or "").strip()
                seq = records[header]
                handle.write(f">{accession.split('.', 1)[0]}\n")
                for i in range(0, len(seq), 80):
                    handle.write(seq[i : i + 80] + "\n")
        outputs[group] = fasta_path
    return outputs


def run_group_annotation(
    group: str,
    input_fasta: Path,
    out_dir: Path,
    keep_existing: bool,
    split: bool,
    cpu: int,
) -> dict[str, str]:
    model_key = f"hrv{group}"
    model_dir = ROOT / "soft" / "vadr-models-hrv" / model_key
    env = build_vadr_env(model_dir)

    out_dir.mkdir(parents=True, exist_ok=True)
    trimmed_fasta = out_dir / f"{model_key}.trimmed.fasta"
    vadr_dir = out_dir / f"{model_key}_vadr"
    prefix = f"{vadr_dir.name}.vadr"
    pass_tbl = vadr_dir / f"{prefix}.pass.tbl"
    fail_tbl = vadr_dir / f"{prefix}.fail.tbl"
    merged_gff = out_dir / f"{model_key}.vadr.gff3"
    split_dir = out_dir / "split_gff3"

    if not keep_existing or not (pass_tbl.exists() or fail_tbl.exists()):
        trim_script = ROOT / "soft" / "vadr" / "miniscripts" / "fasta-trim-terminal-ambigs.pl"
        run_command(
            f"perl {shlex.quote(str(trim_script))} --minlen 50 --maxlen 8000 "
            f"{shlex.quote(str(input_fasta))} > {shlex.quote(str(trimmed_fasta))}",
            env=env,
        )

        vadr_script = ROOT / "soft" / "vadr" / "v-annotate.pl"
        cmd_parts = [
            "perl",
            shlex.quote(str(vadr_script)),
            "-f",
            "-r",
            "--ignore_exc",
        ]
        if split or cpu > 1:
            cmd_parts.extend(["--split", "--cpu", str(max(1, cpu))])
        cmd_parts.extend(
            [
                "--mkey",
                model_key,
                "--mdir",
                shlex.quote(str(model_dir)),
                shlex.quote(str(trimmed_fasta)),
                shlex.quote(str(vadr_dir)),
            ]
        )
        cmd = " ".join(cmd_parts)
        run_command(cmd, env=env)

    source_tbl = pass_tbl if tbl_has_feature_rows(pass_tbl) else fail_tbl
    if not tbl_has_feature_rows(source_tbl):
        raise FileNotFoundError(f"No usable VADR feature table found: {source_tbl}")

    annotate_tbl2gff = ROOT / "soft" / "vadr" / "miniscripts" / "annotate-tbl2gff.pl"
    run_command(
        f"perl {shlex.quote(str(annotate_tbl2gff))} {shlex.quote(str(source_tbl))} > {shlex.quote(str(merged_gff))}",
        env=env,
    )
    split_paths = split_gff_by_seqid(merged_gff, split_dir)
    return {
        "group": group,
        "input_fasta": str(input_fasta),
        "trimmed_fasta": str(trimmed_fasta),
        "vadr_dir": str(vadr_dir),
        "feature_table": str(source_tbl),
        "merged_gff": str(merged_gff),
        "split_gff_dir": str(split_dir),
        "split_gff_count": str(len(split_paths)),
        "split": "yes" if (split or cpu > 1) else "no",
        "cpu": str(max(1, cpu)),
    }


def main() -> int:
    args = parse_args()
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    group_fastas: dict[str, Path] = {}
    if args.manifest:
        sequence_fasta = (args.sequence_fasta or args.input_fasta)
        if sequence_fasta is None:
            raise SystemExit("--manifest mode requires --sequence-fasta or --input-fasta")
        batch_dir = out_dir / "batch_inputs"
        batch_dir.mkdir(parents=True, exist_ok=True)
        group_fastas = build_group_fasta_from_manifest(args.manifest.resolve(), sequence_fasta.resolve(), batch_dir)
    else:
        if args.input_fasta is None or args.group is None:
            raise SystemExit("single-group mode requires both --input-fasta and --group")
        group_fastas = {args.group: args.input_fasta.resolve()}

    summary_rows: list[dict[str, str]] = []
    for group, fasta_path in sorted(group_fastas.items()):
        group_prefix = args.prefix if len(group_fastas) == 1 else f"{args.prefix}.hrv{group}"
        group_out_dir = out_dir / group_prefix
        summary_rows.append(
            run_group_annotation(
                group,
                fasta_path,
                group_out_dir,
                args.keep_existing,
                args.split,
                args.cpu,
            )
        )

    summary_path = out_dir / "vadr_annotation_summary.tsv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "group",
                "input_fasta",
                "trimmed_fasta",
                "vadr_dir",
                "feature_table",
                "merged_gff",
                "split_gff_dir",
                "split_gff_count",
                "split",
                "cpu",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    print(summary_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
