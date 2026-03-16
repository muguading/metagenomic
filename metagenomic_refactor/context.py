from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass
class RuntimeContext:
    ofn: str
    runflow: str
    method: str
    rmhost: str
    tspeabun: str
    nt: int = 0
    krdb: str = ""
    wkdir: str = ""
    long_type: str = ""
    species: str = ""
    genome_len: str = ""
    ref: str = ""
    gtf: str = ""
    vfmeta: Any = None


runtime: RuntimeContext | None = None


def set_runtime_context(ctx: RuntimeContext) -> None:
    global runtime
    runtime = ctx


def get_runtime_context() -> RuntimeContext:
    if runtime is None:
        raise RuntimeError("Runtime context has not been initialized.")
    return runtime


def update_runtime_context(**kwargs) -> None:
    global runtime
    if runtime is None:
        raise RuntimeError("Runtime context has not been initialized.")
    for key, value in kwargs.items():
        if hasattr(runtime, key):
            setattr(runtime, key, value)
