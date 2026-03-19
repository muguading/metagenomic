from __future__ import annotations

import logging

from sqlalchemy.orm import Session

from .models import AuditLog


def configure_audit_logger(name: str = "genome_db.audit") -> logging.Logger:
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler()
        formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s")
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
        logger.propagate = False
    return logger


class AuditTrail:
    def __init__(self, logger: logging.Logger | None = None) -> None:
        self.logger = logger or configure_audit_logger()

    def log(
        self,
        session: Session,
        *,
        operation: str,
        genome_id: str | None,
        operator: str | None,
        status: str,
        details: str | None = None,
    ) -> AuditLog:
        audit_entry = AuditLog(
            operation=operation,
            genome_id=genome_id,
            operator=operator,
            status=status,
            details=details,
        )
        session.add(audit_entry)
        self.logger.info(
            "operation=%s genome_id=%s operator=%s status=%s details=%s",
            operation,
            genome_id,
            operator,
            status,
            details,
        )
        return audit_entry
