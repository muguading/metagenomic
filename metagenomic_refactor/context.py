from __future__ import annotations

from dataclasses import dataclass


@dataclass
class RuntimeContext:
    ofn: str
    runflow: str
    method: str
    rmhost: str
    tspeabun: str
    krdb: str = ""
    wkdir: str = ""


runtime: RuntimeContext | None = None


def set_runtime_context(ctx: RuntimeContext) -> None:
    global runtime
    runtime = ctx


def get_runtime_context() -> RuntimeContext:
    if runtime is None:
        raise RuntimeError("Runtime context has not been initialized.")
    return runtime
