from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime

from .store import PortalStore


@dataclass(frozen=True)
class AdminAuditService:
    store: PortalStore

    def query(self, *, username: str = "", module: str = "", action: str = "", outcome: str = "", search: str = "", limit: int = 500) -> dict:
        items = self.store.list_audit_logs(
            username=username, module=module, action=action, outcome=outcome, search=search, limit=limit
        )
        today_prefix = datetime.utcnow().date().isoformat()
        return {
            "items": items,
            "summary": {
                "total": len(items),
                "failed": sum(1 for item in items if str(item.get("outcome") or "") == "failed"),
                "today": sum(1 for item in items if str(item.get("created_at") or "").startswith(today_prefix)),
                "users": len({str(item.get("username") or "") for item in items if str(item.get("username") or "").strip()}),
            },
            "facets": {
                "modules": sorted({str(item.get("module") or "") for item in items if str(item.get("module") or "").strip()}),
                "actions": sorted({str(item.get("action") or "") for item in items if str(item.get("action") or "").strip()}),
                "users": sorted({str(item.get("username") or "") for item in items if str(item.get("username") or "").strip()}),
            },
        }
