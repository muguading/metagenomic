from __future__ import annotations

import os
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path

from sqlalchemy import create_engine, inspect, text
from sqlalchemy.orm import DeclarativeBase, Session, sessionmaker


PACKAGE_DIR = Path(__file__).resolve().parent
DEFAULT_DB_DIR = PACKAGE_DIR / "data"
DEFAULT_DB_PATH = DEFAULT_DB_DIR / "genome_db.sqlite3"
DEFAULT_DB_URL = os.environ.get("GENOME_DB_URL", f"sqlite:///{DEFAULT_DB_PATH}")


class Base(DeclarativeBase):
    pass


def create_session_factory(database_url: str = DEFAULT_DB_URL) -> sessionmaker[Session]:
    if database_url.startswith("sqlite:///"):
        db_path = Path(database_url.replace("sqlite:///", "", 1)).expanduser()
        db_path.parent.mkdir(parents=True, exist_ok=True)
    connect_args = {"check_same_thread": False} if database_url.startswith("sqlite") else {}
    engine = create_engine(database_url, future=True, echo=False, connect_args=connect_args)
    return sessionmaker(bind=engine, autoflush=False, expire_on_commit=False, future=True)


def init_db(database_url: str = DEFAULT_DB_URL) -> sessionmaker[Session]:
    from .models import AuditLog, Genome, MetadataTemplate, User

    session_factory = create_session_factory(database_url)
    engine = session_factory.kw["bind"]
    Base.metadata.create_all(engine, tables=[Genome.__table__, AuditLog.__table__, User.__table__, MetadataTemplate.__table__])
    _upgrade_sqlite_schema(engine)
    return session_factory


def _upgrade_sqlite_schema(engine) -> None:
    if engine.dialect.name != "sqlite":
        return
    inspector = inspect(engine)
    _upgrade_genomes_schema(engine, inspector)
    _upgrade_metadata_templates_schema(engine, inspector)
    if "users" not in inspector.get_table_names():
        return
    columns = {column["name"] for column in inspector.get_columns("users")}
    with engine.begin() as conn:
        if "display_name" not in columns:
            conn.execute(text("ALTER TABLE users ADD COLUMN display_name VARCHAR(255)"))
        if "email" not in columns:
            conn.execute(text("ALTER TABLE users ADD COLUMN email VARCHAR(255)"))
        if "last_modified_time" not in columns:
            conn.execute(text("ALTER TABLE users ADD COLUMN last_modified_time DATETIME"))
            conn.execute(text("UPDATE users SET last_modified_time = CURRENT_TIMESTAMP WHERE last_modified_time IS NULL"))


def _upgrade_genomes_schema(engine, inspector) -> None:
    if "genomes" not in inspector.get_table_names():
        return
    columns = {column["name"] for column in inspector.get_columns("genomes")}
    with engine.begin() as conn:
        if "custom_metadata" not in columns:
            conn.execute(text("ALTER TABLE genomes ADD COLUMN custom_metadata TEXT"))
        if "gender" not in columns:
            conn.execute(text("ALTER TABLE genomes ADD COLUMN gender VARCHAR(64)"))
        if "country" not in columns:
            conn.execute(text("ALTER TABLE genomes ADD COLUMN country VARCHAR(255)"))
        if "location_json" not in columns:
            conn.execute(text("ALTER TABLE genomes ADD COLUMN location_json TEXT"))
        if "collection_time" not in columns:
            conn.execute(text("ALTER TABLE genomes ADD COLUMN collection_time DATETIME"))
        if "sample_type" not in columns:
            conn.execute(text("ALTER TABLE genomes ADD COLUMN sample_type VARCHAR(255)"))
        if "sequencing_method" not in columns:
            conn.execute(text("ALTER TABLE genomes ADD COLUMN sequencing_method VARCHAR(255)"))
    unique_constraints = {item["name"]: tuple(item["column_names"]) for item in inspector.get_unique_constraints("genomes")}
    desired = ("genome_id", "submitter")
    if unique_constraints.get("uq_genomes_genome_id_submitter") == desired:
        return

    with engine.begin() as conn:
        conn.execute(
            text(
                """
                CREATE TABLE genomes_new (
                    id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
                    genome_id VARCHAR(128) NOT NULL,
                    sample_name VARCHAR(255) NOT NULL,
                    species_name VARCHAR(255) NOT NULL,
                    taxid INTEGER NOT NULL,
                    genome_length INTEGER NOT NULL,
                    genome_file_path VARCHAR(1024) NOT NULL,
                    submitter VARCHAR(255) NOT NULL,
                    submit_time DATETIME NOT NULL DEFAULT CURRENT_TIMESTAMP,
                    last_modified_time DATETIME NOT NULL DEFAULT CURRENT_TIMESTAMP,
                    description TEXT,
                    gender VARCHAR(64),
                    country VARCHAR(255),
                    location_json TEXT,
                    collection_time DATETIME,
                    sample_type VARCHAR(255),
                    sequencing_method VARCHAR(255),
                    custom_metadata TEXT,
                    CONSTRAINT uq_genomes_genome_id_submitter UNIQUE (genome_id, submitter),
                    CONSTRAINT uq_genomes_genome_file_path UNIQUE (genome_file_path)
                )
                """
            )
        )
        conn.execute(
            text(
                """
                INSERT INTO genomes_new (
                    id, genome_id, sample_name, species_name, taxid, genome_length,
                    genome_file_path, submitter, submit_time, last_modified_time, description,
                    gender, country, location_json, collection_time, sample_type, sequencing_method, custom_metadata
                )
                SELECT
                    id, genome_id, sample_name, species_name, taxid, genome_length,
                    genome_file_path, submitter, submit_time, last_modified_time, description,
                    gender, country, location_json, collection_time, sample_type, sequencing_method, custom_metadata
                FROM genomes
                """
            )
        )
        conn.execute(text("DROP TABLE genomes"))
        conn.execute(text("ALTER TABLE genomes_new RENAME TO genomes"))
        conn.execute(text("CREATE INDEX ix_genomes_species_name ON genomes (species_name)"))
        conn.execute(text("CREATE INDEX ix_genomes_taxid ON genomes (taxid)"))
        conn.execute(text("CREATE INDEX ix_genomes_submitter ON genomes (submitter)"))


def _upgrade_metadata_templates_schema(engine, inspector) -> None:
    if "metadata_templates" in inspector.get_table_names():
        return
    with engine.begin() as conn:
        conn.execute(
            text(
                """
                CREATE TABLE metadata_templates (
                    id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
                    submitter VARCHAR(255) NOT NULL,
                    field_key VARCHAR(255) NOT NULL,
                    label VARCHAR(255) NOT NULL,
                    field_type VARCHAR(32) NOT NULL,
                    options_json TEXT,
                    position INTEGER NOT NULL DEFAULT 0,
                    created_time DATETIME NOT NULL DEFAULT CURRENT_TIMESTAMP,
                    last_modified_time DATETIME NOT NULL DEFAULT CURRENT_TIMESTAMP,
                    CONSTRAINT uq_metadata_templates_submitter_field_key UNIQUE (submitter, field_key)
                )
                """
            )
        )
        conn.execute(text("CREATE INDEX ix_metadata_templates_submitter_position ON metadata_templates (submitter, position)"))


@contextmanager
def session_scope(session_factory: sessionmaker[Session]) -> Iterator[Session]:
    session = session_factory()
    try:
        yield session
        session.commit()
    except Exception:
        session.rollback()
        raise
    finally:
        session.close()
