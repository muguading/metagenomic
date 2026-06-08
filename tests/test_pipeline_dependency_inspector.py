from __future__ import annotations

from pathlib import Path

from metagenomic_refactor.pipeline_dependency_inspector import (
    build_conda_env_yaml,
    inspect_python_files,
    render_unknown_markdown,
    summarize_commands,
    write_outputs,
)


def test_inspects_shell_pipeline_and_conda_wrapped_command(tmp_path: Path) -> None:
    script = tmp_path / "pipeline.py"
    script.write_text(
        """
import os
import subprocess

subprocess.run(f"fastp -i {fq1} -o clean.fq && seqkit stat clean.fq | awk '{{print $1}}'", shell=True)
subprocess.run("conda run --no-capture-output -n VFind checkv end_to_end contigs.fa checkv_out", shell=True)
subprocess.run(["bwa", "index", "ref.fa"], check=True)
os.system("bash -lc 'minimap2 -ax sr ref.fa reads.fq | samtools sort -o out.bam'")
""",
        encoding="utf-8",
    )

    result = inspect_python_files([script])
    rows = summarize_commands(result.commands)
    commands = {row["command"]: row for row in rows}

    assert {"fastp", "seqkit", "awk", "checkv", "bwa", "minimap2", "samtools"} <= set(commands)
    assert commands["checkv"]["conda_envs"] == ["VFind"]
    assert commands["fastp"]["group"] == "Read QC"
    assert commands["samtools"]["conda_package"] == "samtools"


def test_dynamic_commands_are_reported_as_skipped(tmp_path: Path) -> None:
    script = tmp_path / "dynamic.py"
    script.write_text(
        """
import subprocess

cmd = build_command()
subprocess.run(cmd, shell=True)
""",
        encoding="utf-8",
    )

    result = inspect_python_files([script])

    assert result.commands == []
    assert result.skipped_dynamic_calls == [
        {"source": str(script), "line": 5, "reason": "dynamic command expression"}
    ]


def test_outputs_env_and_unknown_report(tmp_path: Path) -> None:
    script = tmp_path / "pipeline.py"
    out_dir = tmp_path / "out"
    script.write_text(
        """
import subprocess

subprocess.run("megahit -1 R1.fq -2 R2.fq -o out && mysterytool run x", shell=True)
""",
        encoding="utf-8",
    )

    result = inspect_python_files([script])
    rows = summarize_commands(result.commands)

    env_yml = build_conda_env_yaml(rows, "demo")
    unknown = render_unknown_markdown(rows, [])
    write_outputs(result, out_dir, "demo")

    assert "name: demo" in env_yml
    assert "  - megahit" in env_yml
    assert "`mysterytool`" in unknown
    assert (out_dir / "commands.tsv").exists()
    assert (out_dir / "software_groups.md").exists()
    assert (out_dir / "env.yml").exists()
    assert (out_dir / "unknown_software.md").exists()
    assert (out_dir / "report.json").exists()
