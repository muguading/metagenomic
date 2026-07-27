from __future__ import annotations

import time
from pathlib import Path

from bac_analysis_portal.store import PortalStore


def make_store(tmp_path: Path) -> PortalStore:
    store = PortalStore(db_path=tmp_path / "portal.sqlite3", project_root=tmp_path)
    store.initialize()
    return store


def test_only_one_owner_holds_active_service_lease(tmp_path: Path) -> None:
    first = make_store(tmp_path)
    second = PortalStore(db_path=first.db_path, project_root=tmp_path)

    assert first.acquire_service_lease("queue-maintenance", "first", 60) is True
    assert second.acquire_service_lease("queue-maintenance", "second", 60) is False
    assert first.renew_service_lease("queue-maintenance", "first", 60) is True
    assert second.renew_service_lease("queue-maintenance", "second", 60) is False


def test_expired_or_released_lease_can_be_acquired(tmp_path: Path) -> None:
    first = make_store(tmp_path)
    second = PortalStore(db_path=first.db_path, project_root=tmp_path)

    assert first.acquire_service_lease("queue-maintenance", "first", 1) is True
    with first.connect() as conn:
        conn.execute("UPDATE service_leases SET expires_at = ?", (time.time() - 1,))
    assert second.acquire_service_lease("queue-maintenance", "second", 60) is True
    second.release_service_lease("queue-maintenance", "second")
    assert first.acquire_service_lease("queue-maintenance", "first", 60) is True
