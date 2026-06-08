from __future__ import annotations

import gzip


def checkfile(infile: str) -> str:
    if not infile.endswith("gz"):
        with open(infile, "r") as handle:
            line = handle.readline()
    else:
        with gzip.open(infile, "rb") as handle:
            line = handle.readline().decode(errors="ignore")
    if line and line[0] == ">":
        return "fasta"
    if line and line[0] == "@":
        return "fastq"
    raise ValueError(f"无法识别输入文件格式: {infile}")


def select_cols(row):
    if row[1] == "none":
        return row[:1]
    return row


def format_seconds(seconds: float) -> str:
    total_seconds = int(seconds)
    hours = total_seconds // 3600
    minutes = (total_seconds % 3600) // 60
    remain_seconds = total_seconds % 60
    return f"{hours}小时{minutes}分钟{remain_seconds}秒"
