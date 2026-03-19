from __future__ import annotations

import json
from datetime import datetime

from sqlalchemy import DateTime, Index, Integer, String, Text, UniqueConstraint, func
from sqlalchemy.orm import Mapped, mapped_column

from .database import Base


class Genome(Base):
    __tablename__ = "genomes"
    __table_args__ = (
        UniqueConstraint("genome_id", "submitter", name="uq_genomes_genome_id_submitter"),
        UniqueConstraint("genome_file_path", name="uq_genomes_genome_file_path"),
        Index("ix_genomes_species_name", "species_name"),
        Index("ix_genomes_taxid", "taxid"),
        Index("ix_genomes_submitter", "submitter"),
    )

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    genome_id: Mapped[str] = mapped_column(String(128), nullable=False)
    sample_name: Mapped[str] = mapped_column(String(255), nullable=False)
    species_name: Mapped[str] = mapped_column(String(255), nullable=False)
    taxid: Mapped[int] = mapped_column(Integer, nullable=False)
    genome_length: Mapped[int] = mapped_column(Integer, nullable=False)
    genome_file_path: Mapped[str] = mapped_column(String(1024), nullable=False)
    submitter: Mapped[str] = mapped_column(String(255), nullable=False)
    submit_time: Mapped[datetime] = mapped_column(DateTime(timezone=True), nullable=False, server_default=func.now())
    last_modified_time: Mapped[datetime] = mapped_column(
        DateTime(timezone=True),
        nullable=False,
        server_default=func.now(),
        onupdate=func.now(),
    )
    description: Mapped[str | None] = mapped_column(Text, nullable=True)
    gender: Mapped[str | None] = mapped_column(String(64), nullable=True)
    country: Mapped[str | None] = mapped_column(String(255), nullable=True)
    location_json: Mapped[str | None] = mapped_column(Text, nullable=True)
    collection_time: Mapped[datetime | None] = mapped_column(DateTime(timezone=True), nullable=True)
    sample_type: Mapped[str | None] = mapped_column(String(255), nullable=True)
    sequencing_method: Mapped[str | None] = mapped_column(String(255), nullable=True)
    custom_metadata: Mapped[str | None] = mapped_column(Text, nullable=True)

    def to_dict(self) -> dict[str, object]:
        return {
            "id": self.id,
            "genome_id": self.genome_id,
            "sample_name": self.sample_name,
            "species_name": self.species_name,
            "taxid": self.taxid,
            "genome_length": self.genome_length,
            "genome_file_path": self.genome_file_path,
            "submitter": self.submitter,
            "submit_time": self.submit_time.isoformat() if self.submit_time else None,
            "last_modified_time": self.last_modified_time.isoformat() if self.last_modified_time else None,
            "description": self.description,
            "gender": self.gender,
            "country": self.country,
            "location": json.loads(self.location_json) if self.location_json else None,
            "collection_time": self.collection_time.isoformat() if self.collection_time else None,
            "sample_type": self.sample_type,
            "sequencing_method": self.sequencing_method,
            "custom_metadata": json.loads(self.custom_metadata) if self.custom_metadata else [],
        }


class AuditLog(Base):
    __tablename__ = "audit_log"
    __table_args__ = (
        Index("ix_audit_log_genome_id", "genome_id"),
        Index("ix_audit_log_operation", "operation"),
        Index("ix_audit_log_action_time", "action_time"),
    )

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    operation: Mapped[str] = mapped_column(String(64), nullable=False)
    genome_id: Mapped[str | None] = mapped_column(String(128), nullable=True)
    operator: Mapped[str | None] = mapped_column(String(255), nullable=True)
    status: Mapped[str] = mapped_column(String(32), nullable=False)
    details: Mapped[str | None] = mapped_column(Text, nullable=True)
    action_time: Mapped[datetime] = mapped_column(DateTime(timezone=True), nullable=False, server_default=func.now())

    def to_dict(self) -> dict[str, object]:
        return {
            "operation": self.operation,
            "genome_id": self.genome_id,
            "operator": self.operator,
            "status": self.status,
            "details": self.details,
            "action_time": self.action_time.isoformat() if self.action_time else None,
        }


class User(Base):
    __tablename__ = "users"
    __table_args__ = (
        UniqueConstraint("username", name="uq_users_username"),
        Index("ix_users_role", "role"),
    )

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    username: Mapped[str] = mapped_column(String(255), nullable=False)
    password_hash: Mapped[str] = mapped_column(String(512), nullable=False)
    role: Mapped[str] = mapped_column(String(32), nullable=False, default="user")
    display_name: Mapped[str | None] = mapped_column(String(255), nullable=True)
    email: Mapped[str | None] = mapped_column(String(255), nullable=True)
    created_time: Mapped[datetime] = mapped_column(DateTime(timezone=True), nullable=False, server_default=func.now())
    last_modified_time: Mapped[datetime] = mapped_column(
        DateTime(timezone=True),
        nullable=False,
        server_default=func.now(),
        onupdate=func.now(),
    )
    last_login_time: Mapped[datetime | None] = mapped_column(DateTime(timezone=True), nullable=True)

    def to_dict(self) -> dict[str, object]:
        return {
            "username": self.username,
            "role": self.role,
            "display_name": self.display_name,
            "email": self.email,
            "created_time": self.created_time.isoformat() if self.created_time else None,
            "last_modified_time": self.last_modified_time.isoformat() if self.last_modified_time else None,
            "last_login_time": self.last_login_time.isoformat() if self.last_login_time else None,
        }


class MetadataTemplate(Base):
    __tablename__ = "metadata_templates"
    __table_args__ = (
        UniqueConstraint("submitter", "field_key", name="uq_metadata_templates_submitter_field_key"),
        Index("ix_metadata_templates_submitter_position", "submitter", "position"),
    )

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    submitter: Mapped[str] = mapped_column(String(255), nullable=False)
    field_key: Mapped[str] = mapped_column(String(255), nullable=False)
    label: Mapped[str] = mapped_column(String(255), nullable=False)
    field_type: Mapped[str] = mapped_column(String(32), nullable=False)
    options_json: Mapped[str | None] = mapped_column(Text, nullable=True)
    position: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    created_time: Mapped[datetime] = mapped_column(DateTime(timezone=True), nullable=False, server_default=func.now())
    last_modified_time: Mapped[datetime] = mapped_column(
        DateTime(timezone=True),
        nullable=False,
        server_default=func.now(),
        onupdate=func.now(),
    )

    def to_dict(self) -> dict[str, object]:
        return {
            "key": self.field_key,
            "label": self.label,
            "type": self.field_type,
            "options": json.loads(self.options_json) if self.options_json else [],
            "position": self.position,
        }
