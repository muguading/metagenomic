from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from pathosource_refactor.config import ResourcePaths


@dataclass
class RuntimeContext:
    input_path: str
    run_path: str
    species: str
    threads: int
    meta: str | None
    ref: str
    mode: str
    cgmlstana: str
    gubbins: str
    msamethod: str
    treemethod: str
    bootstrap: int
    mltype: str
    cgmlstversion: str
    speciesdb: pd.DataFrame
    fametadb: pd.DataFrame
    resources: ResourcePaths


runtime: RuntimeContext | None = None


def set_runtime_context(ctx: RuntimeContext) -> None:
    global runtime
    runtime = ctx


def get_runtime_context() -> RuntimeContext:
    if runtime is None:
        raise RuntimeError("Runtime context has not been initialized.")
    return runtime


def update_runtime_context(**kwargs) -> None:
    current = get_runtime_context()
    for key, value in kwargs.items():
        if hasattr(current, key):
            setattr(current, key, value)
