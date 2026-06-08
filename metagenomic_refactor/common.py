from __future__ import annotations

import os
import shlex
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import IO, Sequence

import glob


@dataclass
class CommandResult:
    command: str
    returncode: int
    elapsed_seconds: float


class CommandExecutionError(RuntimeError):
    def __init__(self, command: str, returncode: int):
        super().__init__(f"命令执行失败(returncode={returncode}): {command}")
        self.command = command
        self.returncode = returncode


CONDA_ENV_ALIASES = {
    "VFind": "genomad_aux",
    "geNomad": "genomad_aux",
    "genomad": "genomad_aux",
    "TB_ONT": "meta_main",
    "hamronization": "amr_aux",
    "RGI": "amr_aux",
    "RGI_new": "amr_aux",
    "hostile": "host_filter",
    "kneaddata": "host_filter",
    "BASALT": "mag_aux",
    "mag_binning": "mag_aux",
    "coverm": "mag_aux",
    "gtdbtk": "mag_aux",
    "GTDBtk": "mag_aux",
    "TB": "ncov",
    "PathoSource": "ncov",
    "medaka": "longread_aux",
    "clair3": "longread_aux",
}


def resolve_conda_env_name(env_name: str) -> str:
    return CONDA_ENV_ALIASES.get(str(env_name), str(env_name))


def get_conda_root() -> str:
    raw_root = str(os.environ.get("META_CONDA_ROOT") or "").strip()
    if raw_root:
        return str(Path(raw_root).expanduser())

    raw_exe = str(os.environ.get("META_CONDA_EXE") or os.environ.get("CONDA_EXE") or "").strip()
    if raw_exe and raw_exe != "conda":
        exe_path = Path(raw_exe).expanduser()
        if exe_path.parent.name in {"bin", "condabin", "Scripts"}:
            return str(exe_path.parent.parent)

    which_conda = shutil.which("conda")
    if which_conda:
        exe_path = Path(which_conda)
        if exe_path.parent.name in {"bin", "condabin", "Scripts"}:
            return str(exe_path.parent.parent)

    raw_prefix = str(os.environ.get("CONDA_PREFIX") or "").strip()
    if raw_prefix:
        prefix_path = Path(raw_prefix).expanduser()
        if prefix_path.parent.name == "envs":
            return str(prefix_path.parent.parent)
        return str(prefix_path)

    return ""


def get_conda_exe() -> str:
    raw_exe = str(os.environ.get("META_CONDA_EXE") or "").strip()
    if raw_exe:
        return str(Path(raw_exe).expanduser()) if raw_exe != "conda" else "conda"

    root = get_conda_root()
    if root:
        root_path = Path(root).expanduser()
        candidates = (
            root_path / "bin" / "conda",
            root_path / "condabin" / "conda",
            root_path / "Scripts" / "conda.exe",
            root_path / "conda.exe",
        )
        for candidate in candidates:
            if candidate.is_file():
                return str(candidate)
        return str(candidates[0])

    raw_conda_exe = str(os.environ.get("CONDA_EXE") or "").strip()
    return str(Path(raw_conda_exe).expanduser()) if raw_conda_exe else "conda"


def conda_run_prefix(env_name: str, *, no_capture: bool = True) -> list[str]:
    env_name = resolve_conda_env_name(env_name)
    prefix = [get_conda_exe(), "run"]
    if no_capture:
        prefix.append("--no-capture-output")
    prefix.extend(["-n", str(env_name)])
    return prefix


def conda_run_command(env_name: str, command: str | Sequence[str], *, no_capture: bool = True) -> str:
    prefix = " ".join(shlex.quote(part) for part in conda_run_prefix(env_name, no_capture=no_capture))
    if isinstance(command, str):
        return f"{prefix} {command}".strip()
    return f"{prefix} {' '.join(shlex.quote(str(part)) for part in command)}".strip()


def conda_env_path(env_name: str, *parts: str) -> str:
    env_name = resolve_conda_env_name(env_name)
    root = get_conda_root()
    if root:
        return str(Path(root).expanduser() / "envs" / str(env_name) / Path(*parts))
    return str(Path("envs") / str(env_name) / Path(*parts))


def conda_base_bin(command_name: str) -> str:
    root = get_conda_root()
    if root:
        return str(Path(root).expanduser() / "bin" / command_name)
    return command_name


def is_fasta(filename):
    from Bio import SeqIO

    try:
        for fmt in ("fasta", "fasta-blast", "fasta-pearson"):
            try:
                with open(filename, "r", encoding="utf-8", errors="ignore") as handle:
                    fasta = SeqIO.parse(handle, fmt)
                    if any(fasta):
                        return True
            except ValueError:
                continue
        return False
    except OSError:
        return False


def is_fastq(file_path):
    if file_path != 0:
        suflist = ["fastq", "fq", "fastq.gz", "fq.gz"]
        tlist = [i for i in suflist if str(file_path).endswith(i)]
        return len(tlist) > 0
    return False


def run_command(
    cmd: str,
    logf: IO[str] | None = None,
    check: bool = True,
    cwd: str | None = None,
    env: dict[str, str] | None = None,
) -> CommandResult:
    start = time.time()
    if logf is not None:
        logf.write(f"\n[CMD] {cmd}\n")
        logf.flush()
    merged_env = None
    if env:
        merged_env = os.environ.copy()
        merged_env.update({str(key): str(value) for key, value in env.items()})
    completed = subprocess.run(cmd, shell=True, stdout=logf, stderr=logf, cwd=cwd, env=merged_env)
    elapsed_seconds = time.time() - start
    if logf is not None:
        logf.write(f"[CMD_EXIT] code={completed.returncode} elapsed={elapsed_seconds:.2f}s\n")
        logf.flush()
    if check and completed.returncode != 0:
        raise CommandExecutionError(cmd, completed.returncode)
    return CommandResult(command=cmd, returncode=completed.returncode, elapsed_seconds=elapsed_seconds)


def run_cmd(cmd, logf=None):
    return run_command(cmd, logf=logf, check=True)


def copy_pattern(patterns, dest):
    dest = Path(dest)
    dest.mkdir(parents=True, exist_ok=True)

    for p in patterns:
        for f in glob.glob(p):
            try:
                shutil.copy(f, dest)
            except Exception:
                pass


def format_seconds(seconds):
    total_seconds = int(seconds)
    hours = total_seconds // 3600
    minutes = (total_seconds % 3600) // 60
    seconds = total_seconds % 60
    return f"{hours}小时{minutes}分钟{seconds}秒"


def file_exists(path):
    return bool(path) and str(path) != "0" and os.path.isfile(path)


def is_fastq_input(infile, fastq1, fastq2):
    return is_fastq(infile) or is_fastq(fastq1) or is_fastq(fastq2)


def check_input(infile, intp):
    if intp == "barcode_fastq":
        if os.path.isdir(infile):
            for i in os.listdir(infile):
                if i.startswith("barcode"):
                    return "bardir"
            print("文件夹下未检测到barcode文件夹")
            raise SystemExit()
        print("请输入barcode_fastq文件夹")
        raise SystemExit()

    if intp == "fastq":
        if os.path.isdir(infile):
            return "fqdir"
        if os.path.isfile(infile):
            return "fqfile"

    elif ".cfg" in intp:
        if os.path.isdir(infile):
            return "f5dir"
        print("请输入f5dir文件夹")
        raise SystemExit()

    elif "@" in intp:
        if os.path.isdir(infile):
            return "pod5"
        print("请输入pod5文件夹")
        raise SystemExit()

    elif intp == "fasta":
        if os.path.isdir(infile):
            return "fadir"
        return "fafile"

    return None


def get_free_gpu_memory():
    try:
        result = subprocess.run(
            ["nvidia-smi", "--query-gpu=memory.free", "--format=csv,noheader,nounits"],
            stdout=subprocess.PIPE,
        )
        memory_free = result.stdout.decode("utf-8").strip()
        return memory_free or "noGPU"
    except Exception:
        return "noGPU"


def basecaller(inputtype, params, inf, ofn="basecaller_outputs"):
    freem = get_free_gpu_memory()
    if freem == "noGPU":
        return

    freem = int(freem)
    if not os.path.isdir(ofn):
        os.makedirs(ofn)

    if os.path.isfile(inf):
        os.makedirs(f"{ofn}/tmp")
        subprocess.run(f"cp {inf} {ofn}/tmp", shell=True)
        inf = f"{ofn}/tmp"

    if inputtype == "fast5":
        print("下机数据为fast5模式")
        print(f"guppy_basecaller -r -i {inf} -s {ofn} -x auto -c {params}")
        if freem / 1.2 > 5500:
            subprocess.run(f"guppy_basecaller -r -i {inf} -s {ofn} -x auto -c {params}", shell=True)
        else:
            subprocess.run(
                f"guppy_basecaller -r -i {inf} -s {ofn} -x auto -c {params} --chunk_size 1000",
                shell=True,
            )
        subprocess.run(f"seqkit seq {ofn}/pass/*.fastq > {ofn}/basecaller.fastq", shell=True)
        return

    print("下机数据为pod5模式")
    qdict = {"fast@": 8, "hac@": 9, "sup@": 10}
    configdir = "/home/dell/biosoft/dorado"
    paramlist = params.split("_")
    if len(paramlist) == 4:
        _, _, _, modeln = paramlist
        qthod = qdict.get(modeln)
        subprocess.run(
            f"{configdir}/bin/dorado basecaller -r {configdir}/{params}v3.3 {inf} --min-qscore {qthod} --emit-fastq > {ofn}/basecaller.fastq",
            shell=True,
        )
    else:
        _, _, _, cspeed, modeln = paramlist
        qthod = qdict.get(modeln)
        if cspeed == "400bps":
            subprocess.run(
                f"{configdir}/bin/dorado basecaller -r {configdir}/{params}v4.1.0 {inf} --min-qscore {qthod} --emit-fastq  > {ofn}/basecaller.fastq",
                shell=True,
            )
            if os.path.getsize(f"{ofn}/basecaller.fastq") == 0:
                subprocess.run(
                    f"{configdir}/bin/dorado basecaller -r {configdir}/{params}v4.2.0 {inf} --min-qscore {qthod} --emit-fastq  > {ofn}/basecaller.fastq",
                    shell=True,
                )
        else:
            subprocess.run(
                f"{configdir}/bin/dorado basecaller -r {configdir}/{params}v4.1.0 {inf} --min-qscore {qthod} --emit-fastq  > {ofn}/basecaller.fastq",
                shell=True,
            )
