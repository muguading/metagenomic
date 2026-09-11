from __future__ import annotations

from dataclasses import dataclass

from .store import PortalStore
from .task_manager import AnalysisTaskManager


@dataclass(frozen=True)
class QueueMaintenanceService:
    store: PortalStore
    task_manager: AnalysisTaskManager

    def acquire_lease(self, owner_id: str, *, ttl_seconds: int = 180) -> bool:
        return self.store.acquire_service_lease("queue-maintenance", owner_id, ttl_seconds)

    def renew_lease(self, owner_id: str, *, ttl_seconds: int = 180) -> bool:
        return self.store.renew_service_lease("queue-maintenance", owner_id, ttl_seconds)

    def release_lease(self, owner_id: str) -> None:
        self.store.release_service_lease("queue-maintenance", owner_id)

    def maintain_once(self) -> None:
        self.task_manager.refresh_monitored_tasks()
        self.task_manager.reconcile_queue(
            max(1, int(self.store.get_setting("max_concurrent_tasks", "2") or "2"))
        )
