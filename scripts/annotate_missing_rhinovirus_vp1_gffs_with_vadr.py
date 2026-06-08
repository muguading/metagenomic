#!/usr/bin/env python3
from __future__ import annotations

import csv
import os
import shlex
import subprocess
from pathlib import Path


ROOT = Path("/Users/wuhhh/Desktop/徐老师/代码/metagenomic")
REFERENCE_DIR = ROOT / "database/virus/rhinovirus/reference_genomes"
REPRESENTATIVE_FASTA = REFERENCE_DIR / "human_rhinovirus_representative_genomes.fasta"
REPRESENTATIVE_MANIFEST = REFERENCE_DIR / "human_rhinovirus_representative_genomes.tsv"
ORIGINAL_GFF_DIR = REFERENCE_DIR / "gff3"
VADR_GFF_DIR = REFERENCE_DIR / "gff3_vadr"
VADR_WORK_DIR = REFERENCE_DIR / "vadr_annotation"
STATUS_TSV = VADR_GFF_DIR / "annotation_status.tsv"


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


def parse_gff_attributes(raw: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in raw.split(";"):
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        attrs[key] = value
    return attrs


def gff_has_vp1(gff_path: Path) -> bool:
    if not gff_path.exists():
        return False
    with gff_path.open() as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                continue
            feature_type = fields[2]
            attrs = parse_gff_attributes(fields[8])
            combined = f"{attrs.get('product', '')} {attrs.get('gene', '')} {attrs.get('Name', '')} {attrs.get('Note', '')}".upper()
            if feature_type in {"gene", "CDS", "mature_protein_region_of_CDS"} and ("VP1" in combined or "1D" in combined):
                return True
    return False


def build_vadr_env(model_dir: Path) -> dict[str, str]:
    vadr_root = (ROOT / "soft").resolve()
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
        raise FileNotFoundError("VADR dependencies missing: " + ", ".join(missing[:10]))
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
            ]
            + ([inherited_perl5lib] if inherited_perl5lib else [])
        ),
        "PATH": os.pathsep.join(
            [
                str(vadr_scripts_dir),
                str(blast_bin_dir),
                str(fasta_bin_dir),
                str(infernal_bin_dir),
                str(minimap2_dir),
            ]
            + ([inherited_path] if inherited_path else [])
        ),
    }


def run_command(cmd: str, env: dict[str, str], cwd: Path | None = None) -> None:
    subprocess.run(cmd, shell=True, check=True, env=env, cwd=str(cwd) if cwd else None)


def split_group_gff(group_gff: Path, output_map: dict[str, Path]) -> None:
    buckets: dict[str, list[str]] = {seqid: [] for seqid in output_map}
    current_headers: list[str] = []
    with group_gff.open() as handle:
        for raw_line in handle:
            if raw_line.startswith("#"):
                current_headers.append(raw_line)
                continue
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            seqid = fields[0]
            if seqid not in buckets:
                continue
            if not buckets[seqid]:
                buckets[seqid].extend(current_headers)
            buckets[seqid].append(raw_line)
    for seqid, out_path in output_map.items():
        content = buckets.get(seqid, [])
        if not content:
            continue
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with out_path.open("w") as handle:
            handle.writelines(content)


def tbl_has_feature_rows(tbl_path: Path) -> bool:
    if not tbl_path.exists():
        return False
    with tbl_path.open() as handle:
        for raw_line in handle:
            if raw_line.startswith(">Feature "):
                return True
    return False


def main() -> None:
    VADR_GFF_DIR.mkdir(parents=True, exist_ok=True)
    VADR_WORK_DIR.mkdir(parents=True, exist_ok=True)

    genomes = parse_fasta(REPRESENTATIVE_FASTA)
    with REPRESENTATIVE_MANIFEST.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    missing_by_group: dict[str, list[dict[str, str]]] = {"A": [], "B": [], "C": []}
    for row in rows:
        accession_root = row["accession_root"]
        original_gff = ORIGINAL_GFF_DIR / f"{accession_root}.gff3"
        supplement_gff = VADR_GFF_DIR / f"{accession_root}.gff3"
        if gff_has_vp1(original_gff) or gff_has_vp1(supplement_gff):
            continue
        missing_by_group[row["species_group"]].append(row)

    status_rows: list[dict[str, str]] = []

    for group, group_rows in missing_by_group.items():
        if not group_rows:
            continue
        model_key = f"hrv{group}"
        model_dir = ROOT / "soft" / "vadr-models-hrv" / model_key
        env = build_vadr_env(model_dir)

        group_dir = VADR_WORK_DIR / model_key
        group_dir.mkdir(parents=True, exist_ok=True)
        input_fasta = group_dir / f"{model_key}.input.fasta"
        trimmed_fasta = group_dir / f"{model_key}.trimmed.fasta"
        output_dir = group_dir / f"{model_key}_vadr"
        prefix = f"{output_dir.name}.vadr"
        pass_tbl = output_dir / f"{prefix}.pass.tbl"
        fail_tbl = output_dir / f"{prefix}.fail.tbl"
        merged_gff = group_dir / f"{model_key}.vadr.gff3"

        with input_fasta.open("w") as handle:
            for row in group_rows:
                header = f"{row['record_id']} {row['species_name']} type {row['normalized_type']} accession {row['accession']}"
                seq = genomes[header]
                handle.write(f">{row['accession_root']}\n")
                for i in range(0, len(seq), 80):
                    handle.write(seq[i : i + 80] + "\n")

        trim_script = ROOT / "soft" / "vadr" / "miniscripts" / "fasta-trim-terminal-ambigs.pl"
        run_command(
            f"perl {shlex.quote(str(trim_script))} --minlen 50 --maxlen 8000 {shlex.quote(str(input_fasta))} > {shlex.quote(str(trimmed_fasta))}",
            env=env,
        )

        vadr_script = ROOT / "soft" / "vadr" / "v-annotate.pl"
        cmd_parts = [
            "perl",
            shlex.quote(str(vadr_script)),
            "-f",
            "-r",
            "--ignore_exc",
            "--mkey",
            model_key,
            "--mdir",
            shlex.quote(str(model_dir)),
        ]
        cmd_parts.extend([shlex.quote(str(trimmed_fasta)), shlex.quote(str(output_dir))])
        run_command(" ".join(cmd_parts), env=env)

        annotate_tbl2gff = ROOT / "soft" / "vadr" / "miniscripts" / "annotate-tbl2gff.pl"
        source_tbl = pass_tbl if tbl_has_feature_rows(pass_tbl) else fail_tbl
        if not tbl_has_feature_rows(source_tbl):
            raise FileNotFoundError(f"VADR feature table not found for {model_key}: {source_tbl}")
        run_command(
            f"perl {shlex.quote(str(annotate_tbl2gff))} {shlex.quote(str(source_tbl))} > {shlex.quote(str(merged_gff))}",
            env=env,
        )

        output_map = {row["accession_root"]: VADR_GFF_DIR / f"{row['accession_root']}.gff3" for row in group_rows}
        split_group_gff(merged_gff, output_map)

        failed_accessions: set[str] = set()
        if fail_tbl.exists():
            with fail_tbl.open() as handle:
                for raw_line in handle:
                    line = raw_line.strip()
                    if not line or line.startswith(">"):
                        continue
                    failed_accessions.add(line.split("\t", 1)[0].split()[0])

        for row in group_rows:
            accession_root = row["accession_root"]
            out_gff = output_map[accession_root]
            status_rows.append(
                {
                    "accession_root": accession_root,
                    "accession": row["accession"],
                    "species_group": group,
                    "normalized_type": row["normalized_type"],
                    "status": "annotated" if out_gff.exists() else ("failed" if accession_root in failed_accessions else "missing"),
                    "vadr_gff_path": str(out_gff),
                }
            )

    with STATUS_TSV.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["accession_root", "accession", "species_group", "normalized_type", "status", "vadr_gff_path"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(status_rows)


if __name__ == "__main__":
    main()
