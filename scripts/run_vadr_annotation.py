#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations

import argparse
import os
import shlex
import subprocess
from pathlib import Path


ROOT = Path("/Users/wuhhh/Desktop/徐老师/代码/metagenomic")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run VADR annotation for an input FASTA using a given model directory."
    )
    parser.add_argument(
        "--input-fasta",
        "--input",
        dest="input_fasta",
        required=True,
        help="Input FASTA file to annotate",
    )
    parser.add_argument("--model-dir", required=True, help="VADR model directory, e.g. soft/vadr-models-ev/evB")
    parser.add_argument("--output-dir", required=True, help="Output directory for VADR results")
    parser.add_argument(
        "--mkey",
        default="",
        help="Optional model key. Defaults to model-dir basename, e.g. evB/hrvA/229E",
    )
    parser.add_argument("--minlen", type=int, default=50, help="Minimum sequence length for trim step")
    parser.add_argument("--maxlen", type=int, default=1000000, help="Maximum sequence length for trim step")
    parser.add_argument(
        "--extra-args",
        default="",
        help="Optional raw extra args appended to v-annotate.pl, e.g. '--glsearch --alt_pass discontn'",
    )
    return parser.parse_args()


def run_command(cmd: str, env: dict[str, str], cwd: Path | None = None) -> None:
    subprocess.run(cmd, shell=True, check=True, env=env, cwd=str(cwd) if cwd else None)


def run_command_to_file(cmd_parts: list[str], output_path: Path, env: dict[str, str], cwd: Path | None = None) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        subprocess.run(
            cmd_parts,
            check=True,
            env=env,
            cwd=str(cwd) if cwd else None,
            stdout=handle,
        )


def resolve_vadr_perl_bin() -> str:
    for candidate in ["/usr/bin/perl", "perl"]:
        try:
            completed = subprocess.run(
                [candidate, "-MInline", "-e", "print qq(ok\\n)"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )
            if completed.returncode == 0:
                return candidate
        except Exception:
            continue
    return "perl"


def resolve_working_perl_with_module(env: dict[str, str] | None, module_name: str, candidates: list[str]) -> str:
    for candidate in candidates:
        if not candidate:
            continue
        try:
            completed = subprocess.run(
                [candidate, f"-M{module_name}", "-e", "print qq(ok\\n)"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
                env=env,
            )
            if completed.returncode == 0:
                return candidate
        except Exception:
            continue
    return candidates[0] if candidates else "perl"


def build_vadr_env(model_dir: Path) -> tuple[dict[str, str], str]:
    vadr_root = (ROOT / "soft").resolve()
    vadr_scripts_dir = (vadr_root / "vadr").resolve()
    infernal_bin_dir = (vadr_root / "infernal" / "binaries").resolve()
    bio_easel_dir = (vadr_root / "Bio-Easel").resolve()
    alt_bio_easel_dir = (vadr_root / "Bio-Easel-ncov").resolve()
    if (alt_bio_easel_dir / "blib" / "lib").exists() and (alt_bio_easel_dir / "blib" / "arch").exists():
        bio_easel_dir = alt_bio_easel_dir
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
    env = {
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
                str(blast_bin_dir),
                str(fasta_bin_dir),
                str(infernal_bin_dir),
                str(minimap2_dir),
            ]
            + ([inherited_path] if inherited_path else [])
        ),
    }
    ncov_bin_dir = "/opt/homebrew/Caskroom/mambaforge/base/envs/ncov/bin"
    env["PATH"] = os.pathsep.join([ncov_bin_dir, "/usr/bin", "/bin", str(env.get("PATH") or "")])
    perl_candidates = [f"{ncov_bin_dir}/perl", "/usr/bin/perl", resolve_vadr_perl_bin(), "perl"]
    perl_bin = resolve_working_perl_with_module(env, "Bio::Easel::MSA", perl_candidates)
    return env, perl_bin


def tbl_has_feature_rows(tbl_path: Path) -> bool:
    if not tbl_path.exists():
        return False
    with tbl_path.open() as handle:
        for raw_line in handle:
            if raw_line.startswith(">Feature "):
                return True
    return False


def main() -> int:
    args = parse_args()
    input_fasta = Path(args.input_fasta).resolve()
    model_dir = Path(args.model_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    mkey = args.mkey.strip() or model_dir.name

    if not input_fasta.is_file():
        raise FileNotFoundError(f"Input FASTA not found: {input_fasta}")
    if not model_dir.is_dir():
        raise FileNotFoundError(f"Model directory not found: {model_dir}")

    env, perl_bin = build_vadr_env(model_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    trimmed_fasta = output_dir.parent / f"{output_dir.name}.{input_fasta.stem}.trimmed.fasta"
    prefix = f"{output_dir.name}.vadr"
    pass_tbl = output_dir / f"{prefix}.pass.tbl"
    fail_tbl = output_dir / f"{prefix}.fail.tbl"
    gff_path = output_dir / f"{input_fasta.stem}.vadr.gff3"

    trim_script = ROOT / "soft" / "vadr" / "miniscripts" / "fasta-trim-terminal-ambigs.pl"
    run_command_to_file(
        [
            perl_bin,
            str(trim_script),
            "--minlen",
            str(args.minlen),
            "--maxlen",
            str(args.maxlen),
            str(input_fasta),
        ],
        trimmed_fasta,
        env=env,
    )

    vadr_script = ROOT / "soft" / "vadr" / "v-annotate.pl"
    cmd_parts = [
        shlex.quote(perl_bin),
        shlex.quote(str(vadr_script)),
        "-f",
        "-r",
        "--ignore_exc",
        "--mkey",
        shlex.quote(mkey),
        "--mdir",
        shlex.quote(str(model_dir)),
    ]
    if args.extra_args.strip():
        cmd_parts.append(args.extra_args.strip())
    cmd_parts.extend([shlex.quote(str(trimmed_fasta)), shlex.quote(str(output_dir))])
    run_command(" ".join(cmd_parts), env=env)

    annotate_tbl2gff = ROOT / "soft" / "vadr" / "miniscripts" / "annotate-tbl2gff.pl"
    source_tbl = pass_tbl if tbl_has_feature_rows(pass_tbl) else fail_tbl
    if not tbl_has_feature_rows(source_tbl):
        raise FileNotFoundError(f"No usable VADR feature table found: {source_tbl}")
    run_command(
        f"{shlex.quote(perl_bin)} {shlex.quote(str(annotate_tbl2gff))} {shlex.quote(str(source_tbl))} > {shlex.quote(str(gff_path))}",
        env=env,
    )

    print(f"trimmed_fasta\t{trimmed_fasta}")
    print(f"vadr_output_dir\t{output_dir}")
    print(f"vadr_gff3\t{gff_path}")
    print(f"vadr_pass_tbl\t{pass_tbl}")
    print(f"vadr_fail_tbl\t{fail_tbl}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
