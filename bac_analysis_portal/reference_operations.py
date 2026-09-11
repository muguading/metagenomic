from __future__ import annotations

from dataclasses import dataclass
from typing import Callable


@dataclass(frozen=True)
class ReferenceOperations:
    list_download_jobs: Callable[[str], list[dict[str, object]]]
    delete_download_job: Callable[[str, str], None]
    start_download_job: Callable[..., dict[str, object]]
    import_from_source: Callable[..., dict[str, object]]
    register_record: Callable[..., dict[str, object]]
    batch_import: Callable[..., dict[str, object]]
    precheck_batch: Callable[..., dict[str, object]]
    store_precheck: Callable[..., dict[str, object]]
    load_precheck: Callable[..., tuple[str, bytes, str] | None]
