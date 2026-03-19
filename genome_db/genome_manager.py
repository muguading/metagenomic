from __future__ import annotations

import json
from dataclasses import dataclass
from math import ceil
from pathlib import Path

from sqlalchemy import Select, func, select
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm import Session, sessionmaker

from .audit_logger import AuditTrail
from .database import DEFAULT_DB_URL, init_db, session_scope
from .models import AuditLog, Genome, MetadataTemplate
from .validators import (
    ValidationError,
    validate_custom_metadata,
    validate_fasta_file,
    validate_optional_datetime,
    validate_optional_location,
    validate_optional_text,
    validate_required_text,
    validate_taxid,
)


@dataclass(slots=True)
class SearchResult:
    items: list[dict[str, object]]
    total: int
    page: int
    page_size: int
    pages: int


@dataclass(slots=True)
class AuditSearchResult:
    items: list[dict[str, object]]
    total: int
    page: int
    page_size: int
    pages: int


@dataclass(slots=True)
class DashboardDataResult:
    items: list[dict[str, object]]
    templates: list[dict[str, object]]


class DuplicateGenomeError(ValueError):
    pass


class GenomeManager:
    def __init__(
        self,
        database_url: str = DEFAULT_DB_URL,
        *,
        session_factory: sessionmaker[Session] | None = None,
        audit_trail: AuditTrail | None = None,
    ) -> None:
        self.session_factory = session_factory or init_db(database_url)
        self.audit_trail = audit_trail or AuditTrail()

    def add_genome(
        self,
        *,
        genome_id: str,
        sample_name: str,
        species_name: str,
        taxid: int,
        genome_file_path: str,
        submitter: str,
        description: str | None = None,
        gender: object = None,
        country: object = None,
        location: object = None,
        collection_time: object = None,
        sample_type: object = None,
        sequencing_method: object = None,
        custom_metadata: object = None,
    ) -> dict[str, object]:
        try:
            normalized = self._normalize_input(
                genome_id=genome_id,
                sample_name=sample_name,
                species_name=species_name,
                taxid=taxid,
                genome_file_path=genome_file_path,
                submitter=submitter,
                description=description,
                gender=gender,
                country=country,
                location=location,
                collection_time=collection_time,
                sample_type=sample_type,
                sequencing_method=sequencing_method,
                custom_metadata=custom_metadata,
            )
        except Exception as exc:
            self._log_failure(
                operation="add_genome",
                genome_id=genome_id,
                operator=submitter,
                details=str(exc),
            )
            raise

        with session_scope(self.session_factory) as session:
            try:
                self._ensure_no_duplicate(
                    session,
                    genome_id=normalized["genome_id"],
                    genome_file_path=normalized["genome_file_path"],
                    submitter=normalized["submitter"],
                )
                genome = Genome(**normalized)
                session.add(genome)
                session.flush()
                metadata_items = genome.to_dict().get("custom_metadata", [])
                self._sync_metadata_templates(session, submitter=genome.submitter, metadata_items=metadata_items, prune_missing=False)
                self.audit_trail.log(
                    session,
                    operation="add_genome",
                    genome_id=genome.genome_id,
                    operator=genome.submitter,
                    status="SUCCESS",
                    details=f"Added genome for species={genome.species_name}",
                )
                return genome.to_dict()
            except IntegrityError as exc:
                self._log_failure(
                    operation="add_genome",
                    genome_id=normalized.get("genome_id"),
                    operator=normalized.get("submitter"),
                    details=str(exc.orig),
                )
                raise DuplicateGenomeError("Duplicate genome record detected") from exc
            except Exception as exc:
                self._log_failure(
                    operation="add_genome",
                    genome_id=normalized.get("genome_id"),
                    operator=normalized.get("submitter"),
                    details=str(exc),
                )
                raise

    def update_genome(self, record_id: int, *, operator: str | None = None, **updates: object) -> dict[str, object]:
        with session_scope(self.session_factory) as session:
            before: dict[str, object] = {}
            try:
                genome = self._get_genome_model(session, record_id)
                before = genome.to_dict()
                if "genome_id" in updates and updates["genome_id"] != genome.genome_id:
                    raise ValidationError("genome_id cannot be changed")
                self._apply_updates(genome, updates)
                session.flush()
                self._sync_metadata_templates(
                    session,
                    submitter=genome.submitter,
                    metadata_items=genome.to_dict().get("custom_metadata", []),
                    prune_missing=True,
                )
                self.audit_trail.log(
                    session,
                    operation="update_genome",
                    genome_id=genome.genome_id,
                    operator=operator or genome.submitter,
                    status="SUCCESS",
                    details=f"Updated fields: {', '.join(sorted(updates)) or 'none'}",
                )
                return genome.to_dict()
            except IntegrityError as exc:
                self._log_failure(
                    operation="update_genome",
                    genome_id=before.get("genome_id") if before else None,
                    operator=operator or before.get("submitter"),
                    details=str(exc.orig),
                )
                raise DuplicateGenomeError("Update would create a duplicate genome record") from exc
            except Exception as exc:
                self._log_failure(
                    operation="update_genome",
                    genome_id=before.get("genome_id") if before else None,
                    operator=operator or before.get("submitter"),
                    details=str(exc),
                )
                raise

    def delete_genome(self, record_id: int, *, operator: str | None = None) -> None:
        with session_scope(self.session_factory) as session:
            genome_id: str | None = None
            submitter: str | None = None
            try:
                genome = self._get_genome_model(session, record_id)
                genome_id = genome.genome_id
                submitter = genome.submitter
                session.delete(genome)
                session.flush()
                self.audit_trail.log(
                    session,
                    operation="delete_genome",
                    genome_id=genome_id,
                    operator=operator or submitter,
                    status="SUCCESS",
                    details="Genome record deleted",
                )
            except Exception as exc:
                self._log_failure(
                    operation="delete_genome",
                    genome_id=genome_id,
                    operator=operator or submitter,
                    details=str(exc),
                )
                raise

    def get_genome(self, record_id: int, *, operator: str | None = None) -> dict[str, object]:
        with session_scope(self.session_factory) as session:
            genome_id: str | None = None
            try:
                genome = self._get_genome_model(session, record_id)
                genome_id = genome.genome_id
                self.audit_trail.log(
                    session,
                    operation="get_genome",
                    genome_id=genome.genome_id,
                    operator=operator,
                    status="SUCCESS",
                    details="Genome record retrieved",
                )
                return genome.to_dict()
            except Exception as exc:
                self._log_failure(
                    operation="get_genome",
                    genome_id=genome_id,
                    operator=operator,
                    details=str(exc),
                )
                raise

    def search_genomes(
        self,
        *,
        species_name: str | None = None,
        taxid: int | None = None,
        submitter: str | None = None,
        custom_logic: str = "and",
        custom_filters: object = None,
        page: int = 1,
        page_size: int = 20,
        operator: str | None = None,
    ) -> SearchResult:
        if page < 1 or page_size < 1:
            self._log_failure(
                operation="search_genomes",
                genome_id=None,
                operator=operator,
                details="page and page_size must be positive integers",
            )
            raise ValidationError("page and page_size must be positive integers")

        with session_scope(self.session_factory) as session:
            try:
                stmt = self._build_search_query(species_name=species_name, taxid=taxid, submitter=submitter)
                rows = session.scalars(stmt.order_by(Genome.submit_time.desc(), Genome.id.desc())).all()
                normalized_filters = self._validate_search_filters(custom_filters)
                if normalized_filters:
                    rows = [
                        row for row in rows if self._matches_custom_filters(
                            row.to_dict(),
                            normalized_filters,
                            custom_logic=custom_logic,
                        )
                    ]
                total = len(rows)
                page_rows = rows[(page - 1) * page_size : page * page_size]
                self.audit_trail.log(
                    session,
                    operation="search_genomes",
                    genome_id=None,
                    operator=operator,
                    status="SUCCESS",
                    details=(
                        f"filters species_name={species_name}, taxid={taxid}, submitter={submitter}, "
                        f"custom_logic={custom_logic}, custom_filter_count={len(normalized_filters)}"
                    ),
                )
                return SearchResult(
                    items=[row.to_dict() for row in page_rows],
                    total=total,
                    page=page,
                    page_size=page_size,
                    pages=ceil(total / page_size) if total else 0,
                )
            except Exception as exc:
                self._log_failure(
                    operation="search_genomes",
                    genome_id=None,
                    operator=operator,
                    details=str(exc),
                )
                raise

    def list_audit_logs(
        self,
        *,
        genome_id: str | None = None,
        operation: str | None = None,
        page: int = 1,
        page_size: int = 20,
        operator: str | None = None,
    ) -> AuditSearchResult:
        if page < 1 or page_size < 1:
            raise ValidationError("page and page_size must be positive integers")

        with session_scope(self.session_factory) as session:
            stmt = select(AuditLog)
            if genome_id:
                stmt = stmt.where(AuditLog.genome_id == genome_id.strip())
            if operation:
                stmt = stmt.where(AuditLog.operation.ilike(f"%{operation.strip()}%"))

            total = session.scalar(select(func.count()).select_from(stmt.subquery())) or 0
            rows = session.scalars(
                stmt.order_by(AuditLog.action_time.desc(), AuditLog.id.desc())
                .offset((page - 1) * page_size)
                .limit(page_size)
            ).all()
            return AuditSearchResult(
                items=[row.to_dict() for row in rows],
                total=total,
                page=page,
                page_size=page_size,
                pages=ceil(total / page_size) if total else 0,
            )

    def list_metadata_templates(self, *, submitter: str | None = None) -> list[dict[str, object]]:
        with session_scope(self.session_factory) as session:
            owner = submitter.strip() if submitter else None
            templates = self._ensure_metadata_templates(session, submitter=owner)
            return [template.to_dict() for template in templates]

    def get_dashboard_data(self, *, submitter: str | None = None, operator: str | None = None) -> DashboardDataResult:
        with session_scope(self.session_factory) as session:
            try:
                stmt = self._build_search_query(species_name=None, taxid=None, submitter=submitter)
                rows = session.scalars(stmt.order_by(Genome.submit_time.desc(), Genome.id.desc())).all()
                templates = self._ensure_metadata_templates(session, submitter=submitter.strip() if submitter else None)
                self.audit_trail.log(
                    session,
                    operation="get_dashboard_data",
                    genome_id=None,
                    operator=operator,
                    status="SUCCESS",
                    details=f"dashboard submitter={submitter}",
                )
                return DashboardDataResult(
                    items=[row.to_dict() for row in rows],
                    templates=[template.to_dict() for template in templates],
                )
            except Exception as exc:
                self._log_failure(
                    operation="get_dashboard_data",
                    genome_id=None,
                    operator=operator,
                    details=str(exc),
                )
                raise

    def _normalize_input(
        self,
        *,
        genome_id: str,
        sample_name: str,
        species_name: str,
        taxid: int,
        genome_file_path: str,
        submitter: str,
        description: str | None,
        gender: object,
        country: object,
        location: object,
        collection_time: object,
        sample_type: object,
        sequencing_method: object,
        custom_metadata: object,
    ) -> dict[str, object]:
        resolved_path, genome_length = validate_fasta_file(genome_file_path)
        metadata_items = validate_custom_metadata(custom_metadata)
        normalized_location = validate_optional_location(location, "location")
        return {
            "genome_id": validate_required_text(genome_id, "genome_id"),
            "sample_name": validate_required_text(sample_name, "sample_name"),
            "species_name": validate_required_text(species_name, "species_name"),
            "taxid": validate_taxid(int(taxid)),
            "genome_length": genome_length,
            "genome_file_path": resolved_path,
            "submitter": validate_required_text(submitter, "submitter"),
            "description": description.strip() if isinstance(description, str) and description.strip() else None,
            "gender": validate_optional_text(gender, "gender"),
            "country": validate_optional_text(country, "country"),
            "location_json": json.dumps(normalized_location, ensure_ascii=False) if normalized_location else None,
            "collection_time": validate_optional_datetime(collection_time, "collection_time"),
            "sample_type": validate_optional_text(sample_type, "sample_type"),
            "sequencing_method": validate_optional_text(sequencing_method, "sequencing_method"),
            "custom_metadata": json.dumps(metadata_items, ensure_ascii=False) if metadata_items else None,
        }

    def _sync_metadata_templates(
        self,
        session: Session,
        *,
        submitter: str,
        metadata_items: list[dict[str, object]],
        prune_missing: bool = False,
    ) -> None:
        existing = {
            template.field_key: template
            for template in session.scalars(
                select(MetadataTemplate).where(MetadataTemplate.submitter == submitter).order_by(MetadataTemplate.position.asc())
            ).all()
        }
        next_position = 0
        seen_keys: set[str] = set()
        for item in metadata_items:
            field_key = str(item.get("key", "")).strip()
            if not field_key or field_key in seen_keys:
                continue
            seen_keys.add(field_key)
            options = list(item.get("options", [])) if item.get("type") == "select" else []
            template = existing.get(field_key)
            if template is None:
                session.add(
                    MetadataTemplate(
                        submitter=submitter,
                        field_key=field_key,
                        label=str(item.get("label", "")).strip(),
                        field_type=str(item.get("type", "")).strip(),
                        options_json=json.dumps(options, ensure_ascii=False) if options else None,
                        position=next_position,
                    )
                )
            else:
                template.label = str(item.get("label", "")).strip()
                template.field_type = str(item.get("type", "")).strip()
                template.options_json = json.dumps(options, ensure_ascii=False) if options else None
                template.position = next_position
            next_position += 1

        if prune_missing:
            removed_keys = set(existing) - seen_keys
            if removed_keys:
                self._remove_metadata_fields(session, submitter=submitter, removed_keys=removed_keys)
            return

        for template in session.scalars(
            select(MetadataTemplate).where(MetadataTemplate.submitter == submitter).order_by(MetadataTemplate.position.asc())
        ).all():
            if template.field_key in seen_keys:
                continue
            template.position = next_position
            next_position += 1

    def _remove_metadata_fields(self, session: Session, *, submitter: str, removed_keys: set[str]) -> None:
        if not removed_keys:
            return

        genomes = session.scalars(select(Genome).where(Genome.submitter == submitter)).all()
        for genome in genomes:
            if not genome.custom_metadata:
                continue
            try:
                items = validate_custom_metadata(json.loads(genome.custom_metadata))
            except Exception:
                continue
            filtered_items = [item for item in items if str(item.get("key", "")).strip() not in removed_keys]
            if len(filtered_items) != len(items):
                genome.custom_metadata = json.dumps(filtered_items, ensure_ascii=False) if filtered_items else None

        for template in session.scalars(
            select(MetadataTemplate).where(
                (MetadataTemplate.submitter == submitter) & (MetadataTemplate.field_key.in_(sorted(removed_keys)))
            )
        ).all():
            session.delete(template)

    def _ensure_metadata_templates(self, session: Session, *, submitter: str | None) -> list[MetadataTemplate]:
        stmt = select(MetadataTemplate)
        if submitter:
            stmt = stmt.where(MetadataTemplate.submitter == submitter)
        templates = session.scalars(stmt.order_by(MetadataTemplate.position.asc(), MetadataTemplate.id.asc())).all()
        if templates or submitter is None:
            return templates

        rows = session.scalars(
            select(Genome.custom_metadata)
            .where((Genome.submitter == submitter) & (Genome.custom_metadata.is_not(None)))
            .order_by(Genome.submit_time.asc(), Genome.id.asc())
        ).all()
        synthesized: list[dict[str, object]] = []
        seen_keys: set[str] = set()
        for raw in rows:
            if not raw:
                continue
            try:
                items = validate_custom_metadata(json.loads(raw))
            except Exception:
                continue
            for item in items:
                field_key = str(item.get("key", "")).strip()
                if not field_key or field_key in seen_keys:
                    continue
                seen_keys.add(field_key)
                synthesized.append(item)
        if synthesized:
            self._sync_metadata_templates(session, submitter=submitter, metadata_items=synthesized)
            session.flush()
            templates = session.scalars(
                select(MetadataTemplate)
                .where(MetadataTemplate.submitter == submitter)
                .order_by(MetadataTemplate.position.asc(), MetadataTemplate.id.asc())
            ).all()
        return templates

    def _apply_updates(self, genome: Genome, updates: dict[str, object]) -> None:
        if not updates:
            return

        if "sample_name" in updates:
            genome.sample_name = validate_required_text(str(updates["sample_name"]), "sample_name")
        if "species_name" in updates:
            genome.species_name = validate_required_text(str(updates["species_name"]), "species_name")
        if "taxid" in updates:
            genome.taxid = validate_taxid(int(updates["taxid"]))
        if "submitter" in updates:
            genome.submitter = validate_required_text(str(updates["submitter"]), "submitter")
        if "description" in updates:
            value = updates["description"]
            genome.description = value.strip() if isinstance(value, str) and value.strip() else None
        if "gender" in updates:
            genome.gender = validate_optional_text(updates["gender"], "gender")
        if "country" in updates:
            genome.country = validate_optional_text(updates["country"], "country")
        if "location" in updates:
            normalized_location = validate_optional_location(updates["location"], "location")
            genome.location_json = json.dumps(normalized_location, ensure_ascii=False) if normalized_location else None
        if "collection_time" in updates:
            genome.collection_time = validate_optional_datetime(updates["collection_time"], "collection_time")
        if "sample_type" in updates:
            genome.sample_type = validate_optional_text(updates["sample_type"], "sample_type")
        if "sequencing_method" in updates:
            genome.sequencing_method = validate_optional_text(updates["sequencing_method"], "sequencing_method")
        if "custom_metadata" in updates:
            metadata_items = validate_custom_metadata(updates["custom_metadata"])
            genome.custom_metadata = json.dumps(metadata_items, ensure_ascii=False) if metadata_items else None
        if "genome_file_path" in updates:
            resolved_path, genome_length = validate_fasta_file(str(updates["genome_file_path"]))
            if resolved_path != genome.genome_file_path:
                genome.genome_file_path = resolved_path
                genome.genome_length = genome_length

    def _get_genome_model(self, session: Session, record_id: int) -> Genome:
        stmt = select(Genome).where(Genome.id == record_id)
        genome = session.scalar(stmt)
        if genome is None:
            raise KeyError(f"Genome not found: {record_id}")
        return genome

    def _ensure_no_duplicate(self, session: Session, *, genome_id: str, genome_file_path: str, submitter: str) -> None:
        duplicate_id = session.scalar(
            select(Genome).where((Genome.genome_id == genome_id) & (Genome.submitter == submitter))
        )
        duplicate_path = session.scalar(select(Genome).where(Genome.genome_file_path == genome_file_path))
        if duplicate_id is not None and duplicate_path is not None:
            raise DuplicateGenomeError(
                f"Duplicate genome_id={genome_id!r} for submitter={submitter!r} and genome_file_path={genome_file_path!r}"
            )
        if duplicate_id is not None:
            raise DuplicateGenomeError(f"Duplicate genome_id for this account: {genome_id!r}")
        if duplicate_path is not None:
            raise DuplicateGenomeError(f"Duplicate genome_file_path detected: {genome_file_path!r}")

    def _validate_search_filters(self, value: object) -> list[dict[str, str]]:
        if value in (None, "", []):
            return []
        if not isinstance(value, list):
            raise ValidationError("custom_filters must be a list")

        allowed_operators = {"equals", "not_equals", "contains", "not_contains", "empty", "not_empty"}
        normalized: list[dict[str, str]] = []
        for index, item in enumerate(value, start=1):
            if not isinstance(item, dict):
                raise ValidationError(f"custom_filters[{index}] must be an object")
            key = validate_required_text(str(item.get("key", "")), f"custom_filters[{index}].key")
            operator = validate_required_text(str(item.get("operator", "")), f"custom_filters[{index}].operator").lower()
            if operator not in allowed_operators:
                raise ValidationError(f"custom_filters[{index}].operator is invalid")
            normalized.append(
                {
                    "key": key,
                    "label": str(item.get("label", "")).strip(),
                    "type": str(item.get("type", "")).strip().lower(),
                    "operator": operator,
                    "value": str(item.get("value", "")).strip(),
                }
            )
        return normalized

    def _matches_custom_filters(
        self,
        genome_data: dict[str, object],
        filters: list[dict[str, str]],
        *,
        custom_logic: str,
    ) -> bool:
        normalized_logic = str(custom_logic or "and").strip().lower()
        if normalized_logic not in {"and", "or"}:
            raise ValidationError("custom_logic must be and or or")
        results = [self._matches_single_custom_filter(genome_data, item) for item in filters]
        return all(results) if normalized_logic == "and" else any(results)

    def _matches_single_custom_filter(self, genome_data: dict[str, object], search_filter: dict[str, str]) -> bool:
        search_terms = self._search_terms_for_filter(genome_data, search_filter)
        operator = search_filter["operator"]
        needle = search_filter["value"].strip().casefold()

        if operator == "empty":
            return not search_terms
        if operator == "not_empty":
            return bool(search_terms)
        if not needle:
            return False
        if operator == "equals":
            return any(term == needle for term in search_terms)
        if operator == "not_equals":
            return all(term != needle for term in search_terms)
        if operator == "contains":
            return any(needle in term for term in search_terms)
        if operator == "not_contains":
            return all(needle not in term for term in search_terms)
        return False

    def _search_terms_for_filter(self, genome_data: dict[str, object], search_filter: dict[str, str]) -> list[str]:
        key = str(search_filter.get("key", "")).strip()
        if key.startswith("standard:"):
            return self._standard_search_terms(key.replace("standard:", "", 1), genome_data)
        metadata_items = genome_data.get("custom_metadata", []) if isinstance(genome_data, dict) else []
        metadata_item = next(
            (item for item in metadata_items if str(item.get("key", "")).strip() == key.replace("meta:", "")),
            None,
        )
        return self._metadata_search_terms(metadata_item)

    def _standard_search_terms(self, field_key: str, genome_data: dict[str, object]) -> list[str]:
        def normalize(text: object) -> str:
            return str(text or "").strip().casefold()

        if field_key == "location":
            raw_value = genome_data.get("location")
            if not isinstance(raw_value, dict):
                return []
            parts = [
                normalize(raw_value.get("province")),
                normalize(raw_value.get("city")),
                normalize(raw_value.get("district")),
                normalize(raw_value.get("detail")),
            ]
            joined = " / ".join(part for part in parts if part)
            return [part for part in parts if part] + ([joined] if joined else [])
        value = normalize(genome_data.get(field_key))
        return [value] if value else []

    def _metadata_search_terms(self, item: dict[str, object] | None) -> list[str]:
        if not item:
            return []
        field_type = str(item.get("type", "")).strip().lower()
        raw_value = item.get("value")

        def normalize(text: object) -> str:
            return str(text or "").strip().casefold()

        if field_type == "location" and isinstance(raw_value, dict):
            parts = [
                normalize(raw_value.get("province")),
                normalize(raw_value.get("city")),
                normalize(raw_value.get("district")),
                normalize(raw_value.get("detail")),
            ]
            joined = " / ".join(part for part in parts if part)
            return [part for part in parts if part] + ([joined] if joined else [])
        if field_type == "file":
            terms = [normalize(item.get("filename")), normalize(raw_value)]
            return [term for term in terms if term]
        value = normalize(raw_value)
        return [value] if value else []

    def is_genome_file_referenced(self, genome_file_path: str) -> bool:
        path = str(Path(genome_file_path).expanduser().resolve())
        with session_scope(self.session_factory) as session:
            existing = session.scalar(select(Genome).where(Genome.genome_file_path == path))
            return existing is not None

    def is_metadata_file_referenced(self, file_path: str) -> bool:
        path = str(Path(file_path).expanduser().resolve())
        with session_scope(self.session_factory) as session:
            rows = session.scalars(select(Genome.custom_metadata).where(Genome.custom_metadata.is_not(None))).all()
        for raw in rows:
            if not raw:
                continue
            try:
                items = json.loads(raw)
            except json.JSONDecodeError:
                continue
            for item in items:
                if isinstance(item, dict) and item.get("type") == "file" and item.get("value") == path:
                    return True
        return False

    def _build_search_query(
        self,
        *,
        species_name: str | None,
        taxid: int | None,
        submitter: str | None,
    ) -> Select[tuple[Genome]]:
        stmt = select(Genome)
        if species_name:
            stmt = stmt.where(Genome.species_name.ilike(f"%{species_name.strip()}%"))
        if taxid is not None:
            stmt = stmt.where(Genome.taxid == validate_taxid(int(taxid)))
        if submitter:
            stmt = stmt.where(Genome.submitter.ilike(f"%{submitter.strip()}%"))
        return stmt

    def _log_failure(self, *, operation: str, genome_id: str | None, operator: str | None, details: str) -> None:
        with session_scope(self.session_factory) as audit_session:
            self.audit_trail.log(
                audit_session,
                operation=operation,
                genome_id=genome_id,
                operator=operator,
                status="FAILED",
                details=details,
            )
