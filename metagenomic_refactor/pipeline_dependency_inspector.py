from __future__ import annotations

import argparse
import ast
import json
import shlex
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence


SHELL_CONTROL_TOKENS = {"|", "||", "&&", ";", "\n"}
SHELL_BUILTINS = {
    "alias",
    "break",
    "case",
    "cd",
    "continue",
    "do",
    "done",
    "echo",
    "elif",
    "else",
    "esac",
    "eval",
    "exec",
    "exit",
    "export",
    "fi",
    "for",
    "function",
    "if",
    "printf",
    "pwd",
    "read",
    "return",
    "set",
    "shift",
    "then",
    "trap",
    "type",
    "ulimit",
    "unset",
    "while",
}
REDIRECT_TOKENS = {">", ">>", "<", "2>", "2>>", "&>", ">&"}
COMMAND_WRAPPERS = {"time", "env", "nohup", "command"}


@dataclass(frozen=True)
class SoftwareInfo:
    package: str | None
    group: str
    channel: str = "bioconda"
    note: str = ""


SOFTWARE_CATALOG: dict[str, SoftwareInfo] = {
    "abricate": SoftwareInfo("abricate", "AMR/VF annotation"),
    "bwa": SoftwareInfo("bwa", "Read mapping"),
    "samtools": SoftwareInfo("samtools", "Read mapping"),
    "minimap2": SoftwareInfo("minimap2", "Read mapping"),
    "mosdepth": SoftwareInfo("mosdepth", "Depth/coverage"),
    "bedtools": SoftwareInfo("bedtools", "Genome intervals"),
    "bcftools": SoftwareInfo("bcftools", "Variant calling"),
    "freebayes": SoftwareInfo("freebayes", "Variant calling"),
    "snpeff": SoftwareInfo("snpeff", "Variant annotation"),
    "java": SoftwareInfo("openjdk", "Runtime", channel="conda-forge"),
    "fastp": SoftwareInfo("fastp", "Read QC"),
    "seqkit": SoftwareInfo("seqkit", "Sequence utilities"),
    "rasusa": SoftwareInfo("rasusa", "Read QC"),
    "kneaddata": SoftwareInfo("kneaddata", "Host/background removal"),
    "kraken2": SoftwareInfo("kraken2", "Taxonomy profiling"),
    "bracken": SoftwareInfo("bracken", "Taxonomy profiling"),
    "blastn": SoftwareInfo("blast", "Similarity search"),
    "blastx": SoftwareInfo("blast", "Similarity search"),
    "makeblastdb": SoftwareInfo("blast", "Similarity search"),
    "diamond": SoftwareInfo("diamond", "Similarity search"),
    "taxonkit": SoftwareInfo("taxonkit", "Taxonomy utilities"),
    "megahit": SoftwareInfo("megahit", "Assembly"),
    "spades.py": SoftwareInfo("spades", "Assembly"),
    "spades": SoftwareInfo("spades", "Assembly"),
    "flye": SoftwareInfo("flye", "Assembly"),
    "miniasm": SoftwareInfo("miniasm", "Assembly"),
    "gfatools": SoftwareInfo("gfatools", "Assembly"),
    "canu": SoftwareInfo("canu", "Assembly"),
    "unicycler": SoftwareInfo("unicycler", "Assembly"),
    "raven": SoftwareInfo("raven-assembler", "Assembly"),
    "masurca": SoftwareInfo("masurca", "Assembly"),
    "prokka": SoftwareInfo("prokka", "Genome annotation"),
    "checkm2": SoftwareInfo("checkm2", "MAG quality"),
    "gtdbtk": SoftwareInfo("gtdbtk", "MAG taxonomy"),
    "coverm": SoftwareInfo("coverm", "MAG abundance"),
    "rgi": SoftwareInfo("rgi", "AMR/VF annotation"),
    "staramr": SoftwareInfo("staramr", "AMR/VF annotation"),
    "virsorter": SoftwareInfo("virsorter", "Virus/MGE analysis"),
    "checkv": SoftwareInfo("checkv", "Virus/MGE analysis"),
    "genomad": SoftwareInfo("genomad", "Virus/MGE analysis"),
    "v-annotate.pl": SoftwareInfo("vadr", "Virus annotation"),
    "fasta_generate_regions.py": SoftwareInfo(None, "Local/custom scripts", note="Looks like a project-local helper script."),
    "perbase": SoftwareInfo("perbase", "Depth/coverage"),
    "plasflow.py": SoftwareInfo("plasflow", "Plasmid analysis"),
    "PlasFlow.py": SoftwareInfo("plasflow", "Plasmid analysis"),
    "genovi": SoftwareInfo("genovi", "Visualization/reporting"),
    "ruby": SoftwareInfo("ruby", "Runtime", channel="conda-forge"),
    "perl": SoftwareInfo("perl", "Runtime", channel="conda-forge"),
    "python": SoftwareInfo("python", "Runtime", channel="conda-forge"),
    "python3": SoftwareInfo("python", "Runtime", channel="conda-forge"),
}


SYSTEM_COMMANDS = {
    "awk",
    "bash",
    "cat",
    "cp",
    "cut",
    "find",
    "grep",
    "gunzip",
    "gzip",
    "head",
    "ln",
    "mkdir",
    "mv",
    "rename",
    "rm",
    "sed",
    "sh",
    "sort",
    "tar",
    "touch",
    "uniq",
    "wc",
    "xargs",
}


@dataclass
class CommandUse:
    command: str
    source: Path
    line: int
    raw: str
    conda_env: str | None = None


@dataclass
class InspectionResult:
    commands: list[CommandUse] = field(default_factory=list)
    skipped_dynamic_calls: list[dict[str, str | int]] = field(default_factory=list)


def _basename(command: str) -> str:
    return Path(command).name


def _literal_from_node(node: ast.AST) -> str | list[str] | None:
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value
    if isinstance(node, ast.JoinedStr):
        chunks: list[str] = []
        for part in node.values:
            if isinstance(part, ast.Constant) and isinstance(part.value, str):
                chunks.append(part.value)
            elif isinstance(part, ast.FormattedValue):
                chunks.append("{expr}")
        return "".join(chunks)
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
        left = _literal_from_node(node.left)
        right = _literal_from_node(node.right)
        if isinstance(left, str) and isinstance(right, str):
            return left + right
    if isinstance(node, (ast.List, ast.Tuple)):
        values: list[str] = []
        for element in node.elts:
            literal = _literal_from_node(element)
            if not isinstance(literal, str):
                return None
            values.append(literal)
        return values
    return None


def _called_name(node: ast.AST) -> str:
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        parent = _called_name(node.value)
        return f"{parent}.{node.attr}" if parent else node.attr
    return ""


def _is_command_call(call: ast.Call) -> bool:
    name = _called_name(call.func)
    return name in {
        "os.system",
        "subprocess.call",
        "subprocess.check_call",
        "subprocess.check_output",
        "subprocess.Popen",
        "subprocess.run",
        "run_cmd",
        "run_command",
    } or name.endswith(".run_cmd") or name.endswith(".run_command")


def inspect_python_files(paths: Iterable[Path]) -> InspectionResult:
    result = InspectionResult()
    for path in paths:
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call) or not _is_command_call(node) or not node.args:
                continue
            raw = _literal_from_node(node.args[0])
            if raw is None:
                result.skipped_dynamic_calls.append(
                    {"source": str(path), "line": getattr(node, "lineno", 0), "reason": "dynamic command expression"}
                )
                continue
            line = getattr(node, "lineno", 0)
            if isinstance(raw, list):
                result.commands.extend(_commands_from_argv(raw, path, line))
            else:
                result.commands.extend(_commands_from_shell(raw, path, line))
    result.commands = sorted(result.commands, key=lambda item: (item.command, str(item.source), item.line, item.raw))
    return result


def _commands_from_argv(argv: Sequence[str], source: Path, line: int) -> list[CommandUse]:
    if not argv:
        return []
    command = _basename(argv[0])
    if command in {"bash", "sh"} and len(argv) >= 3 and argv[1] in {"-c", "-lc"}:
        return _commands_from_shell(argv[2], source, line)
    conda = _unwrap_conda_run(list(argv))
    if conda:
        inner, env_name = conda
        return [CommandUse(command=_basename(inner[0]), source=source, line=line, raw=" ".join(argv), conda_env=env_name)]
    return [CommandUse(command=command, source=source, line=line, raw=" ".join(argv))]


def _shell_tokens(command: str) -> list[str]:
    lexer = shlex.shlex(command, posix=True, punctuation_chars="|&;<>")
    lexer.whitespace_split = True
    lexer.commenters = ""
    return list(lexer)


def _commands_from_shell(command: str, source: Path, line: int) -> list[CommandUse]:
    try:
        tokens = _shell_tokens(command)
    except ValueError:
        return [CommandUse(command="<unparseable-shell>", source=source, line=line, raw=command)]

    commands: list[CommandUse] = []
    segment: list[str] = []
    for token in tokens + [";"]:
        if token in SHELL_CONTROL_TOKENS:
            commands.extend(_commands_from_shell_segment(segment, source, line, command))
            segment = []
        else:
            segment.append(token)
    return commands


def _commands_from_shell_segment(tokens: list[str], source: Path, line: int, raw: str) -> list[CommandUse]:
    tokens = _strip_leading_wrappers(tokens)
    if not tokens:
        return []

    command = _basename(tokens[0])
    if command in SHELL_BUILTINS or command in REDIRECT_TOKENS:
        return []
    if command in {"bash", "sh"} and len(tokens) >= 3 and tokens[1] in {"-c", "-lc"}:
        return _commands_from_shell(tokens[2], source, line)

    conda = _unwrap_conda_run(tokens)
    if conda:
        inner, env_name = conda
        return [CommandUse(command=_basename(inner[0]), source=source, line=line, raw=raw, conda_env=env_name)]

    return [CommandUse(command=command, source=source, line=line, raw=raw)]


def _strip_leading_wrappers(tokens: list[str]) -> list[str]:
    stripped = list(tokens)
    while stripped:
        command = _basename(stripped[0])
        if "=" in stripped[0] and not stripped[0].startswith("-"):
            stripped.pop(0)
            continue
        if command in COMMAND_WRAPPERS:
            stripped.pop(0)
            while stripped and stripped[0].startswith("-"):
                stripped.pop(0)
            continue
        break
    return stripped


def _unwrap_conda_run(tokens: list[str]) -> tuple[list[str], str | None] | None:
    if not tokens or _basename(tokens[0]) != "conda" or len(tokens) < 3 or tokens[1] != "run":
        return None
    env_name: str | None = None
    index = 2
    options_with_value = {"-n", "--name", "-p", "--prefix"}
    while index < len(tokens):
        token = tokens[index]
        if token in options_with_value:
            env_name = tokens[index + 1] if index + 1 < len(tokens) else None
            index += 2
            continue
        if token.startswith("-"):
            index += 1
            continue
        return tokens[index:], env_name
    return None


def known_software(command: str) -> SoftwareInfo | None:
    normalized = command.lower()
    return SOFTWARE_CATALOG.get(command) or SOFTWARE_CATALOG.get(normalized)


def summarize_commands(commands: Sequence[CommandUse]) -> list[dict[str, object]]:
    by_command: dict[str, list[CommandUse]] = defaultdict(list)
    for item in commands:
        if item.command in SHELL_BUILTINS or not item.command:
            continue
        by_command[item.command].append(item)

    rows: list[dict[str, object]] = []
    for command, uses in sorted(by_command.items()):
        info = known_software(command)
        rows.append(
            {
                "command": command,
                "count": len(uses),
                "group": info.group if info else ("System utilities" if command in SYSTEM_COMMANDS else "Unknown"),
                "conda_package": info.package if info else (None if command not in SYSTEM_COMMANDS else command),
                "conda_envs": sorted({use.conda_env for use in uses if use.conda_env}),
                "locations": [f"{use.source}:{use.line}" for use in uses],
            }
        )
    return rows


def build_conda_env_yaml(rows: Sequence[dict[str, object]], env_name: str) -> str:
    dependencies: dict[str, str] = {"python": "conda-forge"}
    for row in rows:
        command = str(row["command"])
        info = known_software(command)
        if info and info.package:
            dependencies[info.package] = info.channel

    channels = ["conda-forge", "bioconda", "defaults"]
    lines = [f"name: {env_name}", "channels:"]
    lines.extend(f"  - {channel}" for channel in channels)
    lines.append("dependencies:")
    for package in sorted(dependencies):
        lines.append(f"  - {package}")
    return "\n".join(lines) + "\n"


def render_groups_markdown(rows: Sequence[dict[str, object]]) -> str:
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[str(row["group"])].append(row)

    lines = ["# Software Group Suggestions", ""]
    for group in sorted(grouped):
        lines.append(f"## {group}")
        for row in sorted(grouped[group], key=lambda item: str(item["command"])):
            package = row["conda_package"] or "unknown"
            lines.append(f"- `{row['command']}` -> `{package}` ({row['count']} use(s))")
        lines.append("")
    return "\n".join(lines)


def render_unknown_markdown(rows: Sequence[dict[str, object]], skipped: Sequence[dict[str, str | int]]) -> str:
    unknown = [row for row in rows if row["group"] == "Unknown"]
    lines = ["# Missing / Unknown Software Report", ""]
    if not unknown and not skipped:
        lines.append("No unknown software or dynamic command calls detected.")
        return "\n".join(lines) + "\n"
    if unknown:
        lines.append("## Unknown executables")
        for row in unknown:
            lines.append(f"- `{row['command']}` at {', '.join(row['locations'])}")
        lines.append("")
    if skipped:
        lines.append("## Dynamic command calls skipped")
        for item in skipped:
            lines.append(f"- {item['source']}:{item['line']} ({item['reason']})")
        lines.append("")
    return "\n".join(lines)


def write_outputs(result: InspectionResult, output_dir: Path, env_name: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    rows = summarize_commands(result.commands)
    (output_dir / "commands.tsv").write_text(_commands_tsv(rows), encoding="utf-8")
    (output_dir / "software_groups.md").write_text(render_groups_markdown(rows), encoding="utf-8")
    (output_dir / "env.yml").write_text(build_conda_env_yaml(rows, env_name), encoding="utf-8")
    (output_dir / "unknown_software.md").write_text(render_unknown_markdown(rows, result.skipped_dynamic_calls), encoding="utf-8")
    (output_dir / "report.json").write_text(
        json.dumps({"commands": rows, "skipped_dynamic_calls": result.skipped_dynamic_calls}, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )


def _commands_tsv(rows: Sequence[dict[str, object]]) -> str:
    lines = ["command\tcount\tgroup\tconda_package\tconda_envs\tlocations"]
    for row in rows:
        lines.append(
            "\t".join(
                [
                    str(row["command"]),
                    str(row["count"]),
                    str(row["group"]),
                    str(row["conda_package"] or ""),
                    ",".join(row["conda_envs"]),
                    ";".join(row["locations"]),
                ]
            )
        )
    return "\n".join(lines) + "\n"


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Inspect Python pipeline scripts and propose conda dependencies.")
    parser.add_argument("scripts", nargs="+", type=Path, help="Python pipeline script(s) to inspect.")
    parser.add_argument("-o", "--output-dir", type=Path, default=Path("dependency_inspector_out"))
    parser.add_argument("--env-name", default="pipeline-tools")
    args = parser.parse_args(argv)

    result = inspect_python_files(args.scripts)
    write_outputs(result, args.output_dir, args.env_name)
    rows = summarize_commands(result.commands)
    unknown_count = sum(1 for row in rows if row["group"] == "Unknown")
    print(f"Detected {len(rows)} unique commands; unknown={unknown_count}; outputs={args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
