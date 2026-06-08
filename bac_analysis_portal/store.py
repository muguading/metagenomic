from __future__ import annotations

import json
import sqlite3
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from werkzeug.security import check_password_hash, generate_password_hash

SAMPLE_LIBRARY_COLUMNS = """
sample_key, genome_id, sample_name, task_id, task_name, owner, owner_group, report_dir, output_dir,
final_fasta_path, species_name, pathogen_type, taxid, mlst_species_name, mlst_st, serotype_result, genome_length, q20_rate, q30_rate, completeness,
contamination, contig_count, plasmid_count, total_length, resistance_count,
virulence_count, resistance_gene_hits, virulence_gene_hits, resistance_mge_hits, virulence_mge_hits, description, gender, country, location_json, sample_type, sequencing_method,
custom_metadata_json, sample_alias, sample_source, collection_date, host_info, note,
library_scope, visibility_scope, source_submission_id, imported_at, updated_at
"""

SUBMISSION_COLUMNS = """
request_id, personal_sample_key, sample_name, owner, owner_group, payload_json,
status, review_note, reviewed_by, reviewed_at, created_at, updated_at
"""

VERSION_LOG_COLUMNS = """
event_id, version_label, request_id, sample_key, sample_name, owner, owner_group,
action, summary, operator, payload_json, created_at
"""

RELEASE_VERSION_COLUMNS = """
release_id, version_label, scope, summary, note, operator, change_count,
added_count, updated_count, deleted_count, changes_json, snapshot_json, created_at
"""

HOST_DATABASE_COLUMNS = """
host_key, db_category, host_name, genome_name, taxid, source_type, source_label, source_accession, source_url, description,
fasta_path, index_prefix, index_status, index_message, owner, created_at, updated_at
"""

REFERENCE_PANEL_COLUMNS = """
panel_key, db_category, panel_type, species_name, species_slug, selected_host_keys_json, genome_count,
schema_dir, cgmlst_dir, loci_list_path, results_alleles_path, threshold, status, message, owner, created_at, updated_at
"""

AUDIT_LOG_COLUMNS = """
event_id, username, role, group_name, method, path, endpoint, module, action, target_type, target_id,
status_code, outcome, ip_address, user_agent, request_summary, response_summary, created_at
"""

BATCH_IMPORT_RUN_COLUMNS = """
batch_id, import_type, category, filename, operator, role, group_name, status,
precheck_json, result_json, imported_count, skipped_count, issue_count, content_hash,
created_at, updated_at, completed_at
"""

DEFAULT_ALLOWED_MODULES = ["bacteria", "virus", "metagenome", "community", "pathosource"]
DEFAULT_ALLOWED_VIRUSES = [
    "ncov",
    "flu",
    "rsv",
    "hmpv",
    "hpiv",
    "hadv",
    "rhinovirus",
    "seasonal_hcov",
    "mpox",
    "denv",
    "zikav",
    "chikv",
    "bandavirus",
    "orthohantavirus",
    "orthoebolavirus",
    "norovirus",
    "rotavirus",
    "astroviridae",
    "enterovirus",
    "hepatovirus",
    "hiv",
    "sapovirus",
]


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


@dataclass
class PortalStore:
    db_path: Path
    project_root: Path

    @classmethod
    def from_project_root(cls, project_root: Path) -> "PortalStore":
        db_path = project_root / "bac_analysis_portal.sqlite3"
        store = cls(db_path=db_path, project_root=project_root)
        store.initialize()
        return store

    def initialize(self) -> None:
        with self.connect() as conn:
            conn.execute(
                """
                CREATE TABLE IF NOT EXISTS users (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    username TEXT NOT NULL UNIQUE,
                    password_hash TEXT NOT NULL,
                    role TEXT NOT NULL,
                    group_name TEXT NOT NULL DEFAULT '',
                    allowed_modules TEXT NOT NULL DEFAULT '',
                    module_expirations TEXT NOT NULL DEFAULT '',
                    account_expires_at TEXT NOT NULL DEFAULT '',
                    display_name TEXT,
                    created_at TEXT NOT NULL,
                    updated_at TEXT NOT NULL
                )
                """
            )
            columns = {row["name"] for row in conn.execute("PRAGMA table_info(users)").fetchall()}
            if "group_name" not in columns:
                conn.execute("ALTER TABLE users ADD COLUMN group_name TEXT NOT NULL DEFAULT ''")
            if "allowed_modules" not in columns:
                conn.execute("ALTER TABLE users ADD COLUMN allowed_modules TEXT NOT NULL DEFAULT ''")
            if "allowed_viruses" not in columns:
                conn.execute("ALTER TABLE users ADD COLUMN allowed_viruses TEXT NOT NULL DEFAULT ''")
            if "module_expirations" not in columns:
                conn.execute("ALTER TABLE users ADD COLUMN module_expirations TEXT NOT NULL DEFAULT ''")
            if "account_expires_at" not in columns:
                conn.execute("ALTER TABLE users ADD COLUMN account_expires_at TEXT NOT NULL DEFAULT ''")
            conn.execute("UPDATE users SET role = 'group_admin' WHERE role = 'group'")
            conn.execute(
                """
                CREATE TABLE IF NOT EXISTS settings (
                    key TEXT PRIMARY KEY,
                    value TEXT NOT NULL,
                    updated_at TEXT NOT NULL
                )
                """
            )
            conn.execute(
                """
                CREATE TABLE IF NOT EXISTS sample_library (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    sample_key TEXT NOT NULL UNIQUE,
                    genome_id TEXT NOT NULL DEFAULT '',
                    sample_name TEXT NOT NULL,
                    task_id TEXT NOT NULL,
                    task_name TEXT NOT NULL,
                    owner TEXT NOT NULL,
                    owner_group TEXT NOT NULL DEFAULT '',
                    report_dir TEXT NOT NULL,
                    output_dir TEXT NOT NULL DEFAULT '',
                    final_fasta_path TEXT NOT NULL,
                    species_name TEXT NOT NULL DEFAULT '',
                    pathogen_type TEXT NOT NULL DEFAULT '',
                    taxid TEXT NOT NULL DEFAULT '',
                    mlst_species_name TEXT NOT NULL DEFAULT '',
                    mlst_st TEXT NOT NULL DEFAULT '',
                    serotype_result TEXT NOT NULL DEFAULT '',
                    genome_length TEXT NOT NULL DEFAULT '',
                    q20_rate TEXT NOT NULL DEFAULT '',
                    q30_rate TEXT NOT NULL DEFAULT '',
                    completeness TEXT NOT NULL DEFAULT '',
                    contamination TEXT NOT NULL DEFAULT '',
                    contig_count TEXT NOT NULL DEFAULT '',
                    plasmid_count TEXT NOT NULL DEFAULT '',
                    total_length TEXT NOT NULL DEFAULT '',
                    resistance_count TEXT NOT NULL DEFAULT '',
                    virulence_count TEXT NOT NULL DEFAULT '',
                    resistance_gene_hits TEXT NOT NULL DEFAULT '',
                    virulence_gene_hits TEXT NOT NULL DEFAULT '',
                    resistance_mge_hits TEXT NOT NULL DEFAULT '',
                    virulence_mge_hits TEXT NOT NULL DEFAULT '',
                    description TEXT NOT NULL DEFAULT '',
                    gender TEXT NOT NULL DEFAULT '',
                    country TEXT NOT NULL DEFAULT '',
                    location_json TEXT NOT NULL DEFAULT '',
                    sample_type TEXT NOT NULL DEFAULT '',
                    sequencing_method TEXT NOT NULL DEFAULT '',
                    custom_metadata_json TEXT NOT NULL DEFAULT '[]',
                    sample_alias TEXT NOT NULL DEFAULT '',
                    sample_source TEXT NOT NULL DEFAULT '',
                    collection_date TEXT NOT NULL DEFAULT '',
                    host_info TEXT NOT NULL DEFAULT '',
                    note TEXT NOT NULL DEFAULT '',
                    library_scope TEXT NOT NULL DEFAULT 'main',
                    visibility_scope TEXT NOT NULL DEFAULT 'group',
                    source_submission_id TEXT NOT NULL DEFAULT '',
                    imported_at TEXT NOT NULL,
                    updated_at TEXT NOT NULL
                )
                """
            )
            sample_columns = {
                row["name"]
                for row in conn.execute("PRAGMA table_info(sample_library)").fetchall()
            }
            if "pathogen_type" not in sample_columns:
                conn.execute("ALTER TABLE sample_library ADD COLUMN pathogen_type TEXT NOT NULL DEFAULT ''")
            conn.execute(
                """
                CREATE TABLE IF NOT EXISTS sample_library_submissions (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    request_id TEXT NOT NULL UNIQUE,
                    personal_sample_key TEXT NOT NULL,
                    sample_name TEXT NOT NULL DEFAULT '',
                    owner TEXT NOT NULL DEFAULT '',
                    owner_group TEXT NOT NULL DEFAULT '',
                    payload_json TEXT NOT NULL DEFAULT '{}',
                    status TEXT NOT NULL DEFAULT 'pending',
                    review_note TEXT NOT NULL DEFAULT '',
                    reviewed_by TEXT NOT NULL DEFAULT '',
                    reviewed_at TEXT NOT NULL DEFAULT '',
                    created_at TEXT NOT NULL,
                    updated_at TEXT NOT NULL
                )
                """
            )
            conn.execute(
                """
                CREATE TABLE IF NOT EXISTS sample_library_metadata_templates (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    field_key TEXT NOT NULL UNIQUE,
                    label TEXT NOT NULL DEFAULT '',
                    field_type TEXT NOT NULL DEFAULT 'text',
                    options_json TEXT,
                    config_json TEXT NOT NULL DEFAULT '{}',
                    position INTEGER NOT NULL DEFAULT 0,
                    updated_at TEXT NOT NULL
                )
                """
            )
            conn.execute(
                """
                CREATE TABLE IF NOT EXISTS sample_library_version_logs (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    event_id TEXT NOT NULL UNIQUE,
                    version_label TEXT NOT NULL DEFAULT '',
                    request_id TEXT NOT NULL DEFAULT '',
                    sample_key TEXT NOT NULL DEFAULT '',
                    sample_name TEXT NOT NULL DEFAULT '',
                    owner TEXT NOT NULL DEFAULT '',
                    owner_group TEXT NOT NULL DEFAULT '',
                    action TEXT NOT NULL DEFAULT 'publish',
                    summary TEXT NOT NULL DEFAULT '',
                    operator TEXT NOT NULL DEFAULT '',
                    payload_json TEXT NOT NULL DEFAULT '{}',
                    created_at TEXT NOT NULL
                )
                """
            )
            conn.execute(
                """
                CREATE TABLE IF NOT EXISTS sample_library_release_versions (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    release_id TEXT NOT NULL UNIQUE,
                    version_label TEXT NOT NULL DEFAULT '',
                    scope TEXT NOT NULL DEFAULT 'main',
                    summary TEXT NOT NULL DEFAULT '',
                    note TEXT NOT NULL DEFAULT '',
                    operator TEXT NOT NULL DEFAULT '',
                    change_count INTEGER NOT NULL DEFAULT 0,
                    added_count INTEGER NOT NULL DEFAULT 0,
                    updated_count INTEGER NOT NULL DEFAULT 0,
                    deleted_count INTEGER NOT NULL DEFAULT 0,
                    changes_json TEXT NOT NULL DEFAULT '{}',
                    snapshot_json TEXT NOT NULL DEFAULT '[]',
                    created_at TEXT NOT NULL
                )
                """
            )
            metadata_template_columns = {row["name"] for row in conn.execute("PRAGMA table_info(sample_library_metadata_templates)").fetchall()}
            if "config_json" not in metadata_template_columns:
                conn.execute("ALTER TABLE sample_library_metadata_templates ADD COLUMN config_json TEXT NOT NULL DEFAULT '{}'")
            conn.execute(
                """
                CREATE TABLE IF NOT EXISTS host_database (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    host_key TEXT NOT NULL UNIQUE,
                    db_category TEXT NOT NULL DEFAULT 'host',
                    host_name TEXT NOT NULL DEFAULT '',
                    genome_name TEXT NOT NULL DEFAULT '',
                    taxid TEXT NOT NULL DEFAULT '',
                    source_type TEXT NOT NULL DEFAULT 'local',
                    source_label TEXT NOT NULL DEFAULT '',
                    source_accession TEXT NOT NULL DEFAULT '',
                    source_url TEXT NOT NULL DEFAULT '',
                    description TEXT NOT NULL DEFAULT '',
                    fasta_path TEXT NOT NULL DEFAULT '',
                    index_prefix TEXT NOT NULL DEFAULT '',
                    index_status TEXT NOT NULL DEFAULT 'pending',
                    index_message TEXT NOT NULL DEFAULT '',
                    owner TEXT NOT NULL DEFAULT '',
                    created_at TEXT NOT NULL,
                    updated_at TEXT NOT NULL
                )
                """
            )
            conn.execute(
                """
                CREATE TABLE IF NOT EXISTS reference_panels (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    panel_key TEXT NOT NULL UNIQUE,
                    db_category TEXT NOT NULL DEFAULT 'pathogen',
                    panel_type TEXT NOT NULL DEFAULT 'cgmlst',
                    species_name TEXT NOT NULL DEFAULT '',
                    species_slug TEXT NOT NULL DEFAULT '',
                    selected_host_keys_json TEXT NOT NULL DEFAULT '[]',
                    genome_count INTEGER NOT NULL DEFAULT 0,
                    schema_dir TEXT NOT NULL DEFAULT '',
                    cgmlst_dir TEXT NOT NULL DEFAULT '',
                    loci_list_path TEXT NOT NULL DEFAULT '',
                    results_alleles_path TEXT NOT NULL DEFAULT '',
                    threshold TEXT NOT NULL DEFAULT '0.95',
                    status TEXT NOT NULL DEFAULT 'pending',
                    message TEXT NOT NULL DEFAULT '',
                    owner TEXT NOT NULL DEFAULT '',
                    created_at TEXT NOT NULL,
                    updated_at TEXT NOT NULL
                )
                """
            )
            conn.execute(
                """
                CREATE TABLE IF NOT EXISTS audit_logs (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    event_id TEXT NOT NULL UNIQUE,
                    username TEXT NOT NULL DEFAULT '',
                    role TEXT NOT NULL DEFAULT '',
                    group_name TEXT NOT NULL DEFAULT '',
                    method TEXT NOT NULL DEFAULT '',
                    path TEXT NOT NULL DEFAULT '',
                    endpoint TEXT NOT NULL DEFAULT '',
                    module TEXT NOT NULL DEFAULT '',
                    action TEXT NOT NULL DEFAULT '',
                    target_type TEXT NOT NULL DEFAULT '',
                    target_id TEXT NOT NULL DEFAULT '',
                    status_code INTEGER NOT NULL DEFAULT 0,
                    outcome TEXT NOT NULL DEFAULT '',
                    ip_address TEXT NOT NULL DEFAULT '',
                    user_agent TEXT NOT NULL DEFAULT '',
                    request_summary TEXT NOT NULL DEFAULT '',
                    response_summary TEXT NOT NULL DEFAULT '',
                    created_at TEXT NOT NULL
                )
                """
            )
            conn.execute(
                """
                CREATE TABLE IF NOT EXISTS batch_import_runs (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    batch_id TEXT NOT NULL UNIQUE,
                    import_type TEXT NOT NULL DEFAULT '',
                    category TEXT NOT NULL DEFAULT '',
                    filename TEXT NOT NULL DEFAULT '',
                    operator TEXT NOT NULL DEFAULT '',
                    role TEXT NOT NULL DEFAULT '',
                    group_name TEXT NOT NULL DEFAULT '',
                    status TEXT NOT NULL DEFAULT 'precheck_failed',
                    precheck_json TEXT NOT NULL DEFAULT '{}',
                    result_json TEXT NOT NULL DEFAULT '{}',
                    imported_count INTEGER NOT NULL DEFAULT 0,
                    skipped_count INTEGER NOT NULL DEFAULT 0,
                    issue_count INTEGER NOT NULL DEFAULT 0,
                    content_hash TEXT NOT NULL DEFAULT '',
                    created_at TEXT NOT NULL,
                    updated_at TEXT NOT NULL,
                    completed_at TEXT NOT NULL DEFAULT ''
                )
                """
            )
            host_columns = {row["name"] for row in conn.execute("PRAGMA table_info(host_database)").fetchall()}
            host_extra_columns = {
                "db_category": "TEXT NOT NULL DEFAULT 'host'",
                "genome_name": "TEXT NOT NULL DEFAULT ''",
                "taxid": "TEXT NOT NULL DEFAULT ''",
                "source_label": "TEXT NOT NULL DEFAULT ''",
                "source_accession": "TEXT NOT NULL DEFAULT ''",
                "source_url": "TEXT NOT NULL DEFAULT ''",
                "description": "TEXT NOT NULL DEFAULT ''",
                "index_prefix": "TEXT NOT NULL DEFAULT ''",
                "index_status": "TEXT NOT NULL DEFAULT 'pending'",
                "index_message": "TEXT NOT NULL DEFAULT ''",
                "owner": "TEXT NOT NULL DEFAULT ''",
            }
            for column_name, column_type in host_extra_columns.items():
                if column_name not in host_columns:
                    conn.execute(f"ALTER TABLE host_database ADD COLUMN {column_name} {column_type}")
            panel_columns = {row["name"] for row in conn.execute("PRAGMA table_info(reference_panels)").fetchall()}
            panel_extra_columns = {
                "db_category": "TEXT NOT NULL DEFAULT 'pathogen'",
                "panel_type": "TEXT NOT NULL DEFAULT 'cgmlst'",
                "species_name": "TEXT NOT NULL DEFAULT ''",
                "species_slug": "TEXT NOT NULL DEFAULT ''",
                "selected_host_keys_json": "TEXT NOT NULL DEFAULT '[]'",
                "genome_count": "INTEGER NOT NULL DEFAULT 0",
                "schema_dir": "TEXT NOT NULL DEFAULT ''",
                "cgmlst_dir": "TEXT NOT NULL DEFAULT ''",
                "loci_list_path": "TEXT NOT NULL DEFAULT ''",
                "results_alleles_path": "TEXT NOT NULL DEFAULT ''",
                "threshold": "TEXT NOT NULL DEFAULT '0.95'",
                "status": "TEXT NOT NULL DEFAULT 'pending'",
                "message": "TEXT NOT NULL DEFAULT ''",
                "owner": "TEXT NOT NULL DEFAULT ''",
            }
            for column_name, column_type in panel_extra_columns.items():
                if column_name not in panel_columns:
                    conn.execute(f"ALTER TABLE reference_panels ADD COLUMN {column_name} {column_type}")
            sample_columns = {row["name"] for row in conn.execute("PRAGMA table_info(sample_library)").fetchall()}
            sample_extra_columns = {
                "genome_id": "TEXT NOT NULL DEFAULT ''",
                "taxid": "TEXT NOT NULL DEFAULT ''",
                "mlst_st": "TEXT NOT NULL DEFAULT ''",
                "serotype_result": "TEXT NOT NULL DEFAULT ''",
                "resistance_gene_hits": "TEXT NOT NULL DEFAULT ''",
                "virulence_gene_hits": "TEXT NOT NULL DEFAULT ''",
                "resistance_mge_hits": "TEXT NOT NULL DEFAULT ''",
                "virulence_mge_hits": "TEXT NOT NULL DEFAULT ''",
                "genome_length": "TEXT NOT NULL DEFAULT ''",
                "description": "TEXT NOT NULL DEFAULT ''",
                "gender": "TEXT NOT NULL DEFAULT ''",
                "country": "TEXT NOT NULL DEFAULT ''",
                "location_json": "TEXT NOT NULL DEFAULT ''",
                "sample_type": "TEXT NOT NULL DEFAULT ''",
                "sequencing_method": "TEXT NOT NULL DEFAULT ''",
                "custom_metadata_json": "TEXT NOT NULL DEFAULT '[]'",
                "sample_alias": "TEXT NOT NULL DEFAULT ''",
                "sample_source": "TEXT NOT NULL DEFAULT ''",
                "collection_date": "TEXT NOT NULL DEFAULT ''",
                "host_info": "TEXT NOT NULL DEFAULT ''",
                "note": "TEXT NOT NULL DEFAULT ''",
                "library_scope": "TEXT NOT NULL DEFAULT 'main'",
                "visibility_scope": "TEXT NOT NULL DEFAULT 'group'",
                "source_submission_id": "TEXT NOT NULL DEFAULT ''",
            }
            for column_name, column_type in sample_extra_columns.items():
                if column_name not in sample_columns:
                    conn.execute(f"ALTER TABLE sample_library ADD COLUMN {column_name} {column_type}")
            admin = conn.execute("SELECT username FROM users WHERE role = 'admin' LIMIT 1").fetchone()
            if admin is None:
                now = utc_now_iso()
                conn.execute(
                    """
                    INSERT INTO users (username, password_hash, role, display_name, created_at, updated_at)
                    VALUES (?, ?, ?, ?, ?, ?)
                    """,
                    ("admin", generate_password_hash("admin123"), "admin", "System Admin", now, now),
                )
            default_workspace_root = str(self.project_root)
            default_script = "Bac_assemble_260112_newformat.py"
            conn.execute(
                """
                INSERT INTO settings (key, value, updated_at)
                VALUES (?, ?, ?)
                ON CONFLICT(key) DO NOTHING
                """,
                ("workspace_root", default_workspace_root, utc_now_iso()),
            )
            conn.execute(
                """
                INSERT INTO settings (key, value, updated_at)
                VALUES (?, ?, ?)
                ON CONFLICT(key) DO NOTHING
                """,
                ("pipeline_script", default_script, utc_now_iso()),
            )
            conn.execute(
                """
                INSERT INTO settings (key, value, updated_at)
                VALUES (?, ?, ?)
                ON CONFLICT(key) DO NOTHING
                """,
                ("pipeline_python", sys.executable, utc_now_iso()),
            )
            conn.execute(
                """
                INSERT INTO settings (key, value, updated_at)
                VALUES (?, ?, ?)
                ON CONFLICT(key) DO NOTHING
                """,
                ("database_root", str(self.project_root), utc_now_iso()),
            )
            conn.execute(
                """
                INSERT INTO settings (key, value, updated_at)
                VALUES (?, ?, ?)
                ON CONFLICT(key) DO NOTHING
                """,
                ("max_concurrent_tasks", "2", utc_now_iso()),
            )

    def connect(self) -> sqlite3.Connection:
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row
        return conn

    def authenticate(self, username: str, password: str) -> dict[str, Any] | None:
        with self.connect() as conn:
            row = conn.execute(
                "SELECT username, password_hash, role, group_name, allowed_modules, allowed_viruses, module_expirations, account_expires_at, display_name, created_at, updated_at FROM users WHERE username = ?",
                (username,),
            ).fetchone()
        if row is None or not check_password_hash(row["password_hash"], password):
            return None
        serialized = self._serialize_user_row(row)
        if serialized.get("is_expired"):
            return None
        return serialized

    def get_user(self, username: str) -> dict[str, Any]:
        with self.connect() as conn:
            row = conn.execute(
                "SELECT username, role, group_name, allowed_modules, allowed_viruses, module_expirations, account_expires_at, display_name, created_at, updated_at FROM users WHERE username = ?",
                (username,),
            ).fetchone()
        if row is None:
            raise KeyError(f"用户不存在: {username}")
        return self._serialize_user_row(row)

    def list_users(self) -> list[dict[str, Any]]:
        with self.connect() as conn:
            rows = conn.execute(
                "SELECT username, role, group_name, allowed_modules, allowed_viruses, module_expirations, account_expires_at, display_name, created_at, updated_at FROM users ORDER BY role DESC, group_name ASC, username ASC"
            ).fetchall()
        return [self._serialize_user_row(row) for row in rows]

    def create_user(
        self,
        username: str,
        password: str,
        role: str,
        display_name: str | None = None,
        group_name: str | None = None,
        allowed_modules: list[str] | None = None,
        allowed_viruses: list[str] | None = None,
        module_expirations: dict[str, str] | None = None,
        account_expires_at: str | None = None,
    ) -> dict[str, Any]:
        role = self._normalize_role(role)
        if role not in {"admin", "group_admin", "user"}:
            raise ValueError("role must be admin, group_admin or user")
        now = utc_now_iso()
        normalized_group = "" if role == "admin" else (group_name or "").strip()
        normalized_modules = self._normalize_allowed_modules(allowed_modules, role=role, default_if_empty=["bacteria"])
        normalized_viruses = self._normalize_allowed_viruses(allowed_viruses, default_if_empty=DEFAULT_ALLOWED_VIRUSES)
        normalized_module_expirations = self._normalize_module_expirations(module_expirations, granted_modules=normalized_modules, role=role)
        normalized_expiration = self._normalize_account_expires_at(account_expires_at)
        if role != "admin" and not normalized_group:
            raise ValueError("group_name is required for non-admin users")
        try:
            with self.connect() as conn:
                conn.execute(
                    """
                    INSERT INTO users (username, password_hash, role, group_name, allowed_modules, allowed_viruses, module_expirations, account_expires_at, display_name, created_at, updated_at)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    (username, generate_password_hash(password), role, normalized_group, json.dumps(normalized_modules, ensure_ascii=False), json.dumps(normalized_viruses, ensure_ascii=False), json.dumps(normalized_module_expirations, ensure_ascii=False), normalized_expiration, display_name or "", now, now),
                )
        except sqlite3.IntegrityError as exc:
            raise ValueError(f"用户名已存在: {username}") from exc
        return self.get_user(username)

    def update_user(
        self,
        username: str,
        *,
        new_username: str | None = None,
        role: str | None = None,
        group_name: str | None = None,
        display_name: str | None = None,
        new_password: str | None = None,
        allowed_modules: list[str] | None = None,
        allowed_viruses: list[str] | None = None,
        module_expirations: dict[str, str] | None = None,
        account_expires_at: str | None = None,
    ) -> dict[str, Any]:
        current = self.get_user_with_hash(username)
        next_username = (new_username or current["username"]).strip()
        if not next_username:
            raise ValueError("用户名不能为空")
        next_role = self._normalize_role(role or current["role"])
        if next_role not in {"admin", "group_admin", "user"}:
            raise ValueError("role must be admin, group_admin or user")
        next_group_name = group_name.strip() if group_name is not None else (current.get("group_name") or "")
        if next_role == "admin":
            next_group_name = ""
        elif not next_group_name:
            raise ValueError("group_name is required for non-admin users")
        next_allowed_modules = self._normalize_allowed_modules(
            allowed_modules if allowed_modules is not None else current.get("allowed_modules"),
            role=next_role,
            default_if_empty=current.get("granted_modules") or current.get("allowed_modules"),
        )
        next_allowed_viruses = self._normalize_allowed_viruses(
            allowed_viruses if allowed_viruses is not None else current.get("allowed_viruses"),
            default_if_empty=current.get("granted_viruses") or current.get("allowed_viruses") or DEFAULT_ALLOWED_VIRUSES,
        )
        next_module_expirations = self._normalize_module_expirations(
            module_expirations if module_expirations is not None else current.get("module_expirations"),
            granted_modules=next_allowed_modules,
            role=next_role,
        )
        next_account_expires_at = (
            self._normalize_account_expires_at(account_expires_at)
            if account_expires_at is not None
            else self._normalize_account_expires_at(current.get("account_expires_at"))
        )
        next_display_name = display_name if display_name is not None else current["display_name"]
        next_hash = generate_password_hash(new_password) if new_password else current["password_hash"]
        try:
            with self.connect() as conn:
                conn.execute(
                    """
                    UPDATE users
                    SET username = ?, role = ?, group_name = ?, allowed_modules = ?, allowed_viruses = ?, module_expirations = ?, account_expires_at = ?, display_name = ?, password_hash = ?, updated_at = ?
                    WHERE username = ?
                    """,
                    (next_username, next_role, next_group_name, json.dumps(next_allowed_modules, ensure_ascii=False), json.dumps(next_allowed_viruses, ensure_ascii=False), json.dumps(next_module_expirations, ensure_ascii=False), next_account_expires_at, next_display_name, next_hash, utc_now_iso(), username),
                )
        except sqlite3.IntegrityError as exc:
            raise ValueError(f"用户名已存在: {next_username}") from exc
        return self.get_user(next_username)

    def delete_user(self, username: str) -> None:
        with self.connect() as conn:
            row = conn.execute("SELECT username, role FROM users WHERE username = ?", (username,)).fetchone()
            if row is None:
                raise KeyError(f"用户不存在: {username}")
            if row["role"] == "admin":
                admin_count = conn.execute("SELECT COUNT(*) AS count FROM users WHERE role = 'admin'").fetchone()["count"]
                if admin_count <= 1:
                    raise ValueError("不能删除最后一个管理员")
            conn.execute("DELETE FROM users WHERE username = ?", (username,))

    def get_user_with_hash(self, username: str) -> dict[str, Any]:
        with self.connect() as conn:
            row = conn.execute("SELECT * FROM users WHERE username = ?", (username,)).fetchone()
        if row is None:
            raise KeyError(f"用户不存在: {username}")
        return dict(row)

    def get_setting(self, key: str, default: str = "") -> str:
        with self.connect() as conn:
            row = conn.execute("SELECT value FROM settings WHERE key = ?", (key,)).fetchone()
        return row["value"] if row else default

    def set_setting(self, key: str, value: str) -> None:
        with self.connect() as conn:
            conn.execute(
                """
                INSERT INTO settings (key, value, updated_at)
                VALUES (?, ?, ?)
                ON CONFLICT(key) DO UPDATE SET value = excluded.value, updated_at = excluded.updated_at
                """,
                (key, value, utc_now_iso()),
            )

    def record_audit_log(self, payload: dict[str, Any]) -> dict[str, Any]:
        now = payload.get("created_at") or utc_now_iso()
        event_id = str(payload.get("event_id") or f"audit-{now}-{hash(json.dumps(payload, ensure_ascii=False, sort_keys=True))}").strip()
        record = {
            "event_id": event_id,
            "username": str(payload.get("username") or "").strip(),
            "role": str(payload.get("role") or "").strip(),
            "group_name": str(payload.get("group_name") or "").strip(),
            "method": str(payload.get("method") or "").strip(),
            "path": str(payload.get("path") or "").strip(),
            "endpoint": str(payload.get("endpoint") or "").strip(),
            "module": str(payload.get("module") or "").strip(),
            "action": str(payload.get("action") or "").strip(),
            "target_type": str(payload.get("target_type") or "").strip(),
            "target_id": str(payload.get("target_id") or "").strip(),
            "status_code": int(payload.get("status_code") or 0),
            "outcome": str(payload.get("outcome") or "").strip(),
            "ip_address": str(payload.get("ip_address") or "").strip(),
            "user_agent": str(payload.get("user_agent") or "").strip(),
            "request_summary": str(payload.get("request_summary") or "").strip(),
            "response_summary": str(payload.get("response_summary") or "").strip(),
            "created_at": str(now),
        }
        with self.connect() as conn:
            conn.execute(
                """
                INSERT OR REPLACE INTO audit_logs (
                    event_id, username, role, group_name, method, path, endpoint, module, action, target_type, target_id,
                    status_code, outcome, ip_address, user_agent, request_summary, response_summary, created_at
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    record["event_id"],
                    record["username"],
                    record["role"],
                    record["group_name"],
                    record["method"],
                    record["path"],
                    record["endpoint"],
                    record["module"],
                    record["action"],
                    record["target_type"],
                    record["target_id"],
                    record["status_code"],
                    record["outcome"],
                    record["ip_address"],
                    record["user_agent"],
                    record["request_summary"],
                    record["response_summary"],
                    record["created_at"],
                ),
            )
        return record

    def list_audit_logs(
        self,
        *,
        username: str = "",
        module: str = "",
        action: str = "",
        outcome: str = "",
        search: str = "",
        limit: int = 500,
    ) -> list[dict[str, Any]]:
        conditions = []
        params: list[Any] = []
        if username.strip():
            conditions.append("username = ?")
            params.append(username.strip())
        if module.strip():
            conditions.append("module = ?")
            params.append(module.strip())
        if action.strip():
            conditions.append("action = ?")
            params.append(action.strip())
        if outcome.strip():
            conditions.append("outcome = ?")
            params.append(outcome.strip())
        if search.strip():
            keyword = f"%{search.strip()}%"
            conditions.append("(path LIKE ? OR request_summary LIKE ? OR response_summary LIKE ? OR target_id LIKE ?)")
            params.extend([keyword, keyword, keyword, keyword])
        where_clause = f"WHERE {' AND '.join(conditions)}" if conditions else ""
        capped_limit = max(1, min(int(limit or 500), 2000))
        with self.connect() as conn:
            rows = conn.execute(
                f"""
                SELECT {AUDIT_LOG_COLUMNS}
                FROM audit_logs
                {where_clause}
                ORDER BY created_at DESC, id DESC
                LIMIT ?
                """,
                (*params, capped_limit),
            ).fetchall()
        return [dict(row) for row in rows]

    def create_batch_import_run(self, payload: dict[str, Any]) -> dict[str, Any]:
        now = str(payload.get("created_at") or utc_now_iso())
        batch_id = str(payload.get("batch_id") or "").strip()
        if not batch_id:
            raise ValueError("batch_id 不能为空")
        precheck = payload.get("precheck")
        result = payload.get("result")
        status = str(payload.get("status") or "precheck_failed").strip() or "precheck_failed"
        record = {
            "batch_id": batch_id,
            "import_type": str(payload.get("import_type") or "").strip(),
            "category": str(payload.get("category") or "").strip(),
            "filename": str(payload.get("filename") or "").strip(),
            "operator": str(payload.get("operator") or "").strip(),
            "role": str(payload.get("role") or "").strip(),
            "group_name": str(payload.get("group_name") or "").strip(),
            "status": status,
            "precheck_json": json.dumps(precheck if isinstance(precheck, dict) else {}, ensure_ascii=False),
            "result_json": json.dumps(result if isinstance(result, dict) else {}, ensure_ascii=False),
            "imported_count": int(payload.get("imported_count") or 0),
            "skipped_count": int(payload.get("skipped_count") or 0),
            "issue_count": int(payload.get("issue_count") or 0),
            "content_hash": str(payload.get("content_hash") or "").strip(),
            "created_at": now,
            "updated_at": now,
            "completed_at": str(payload.get("completed_at") or "").strip(),
        }
        with self.connect() as conn:
            conn.execute(
                """
                INSERT INTO batch_import_runs (
                    batch_id, import_type, category, filename, operator, role, group_name, status,
                    precheck_json, result_json, imported_count, skipped_count, issue_count, content_hash,
                    created_at, updated_at, completed_at
                )
                VALUES (
                    :batch_id, :import_type, :category, :filename, :operator, :role, :group_name, :status,
                    :precheck_json, :result_json, :imported_count, :skipped_count, :issue_count, :content_hash,
                    :created_at, :updated_at, :completed_at
                )
                ON CONFLICT(batch_id) DO UPDATE SET
                    import_type=excluded.import_type,
                    category=excluded.category,
                    filename=excluded.filename,
                    operator=excluded.operator,
                    role=excluded.role,
                    group_name=excluded.group_name,
                    status=excluded.status,
                    precheck_json=excluded.precheck_json,
                    result_json=excluded.result_json,
                    imported_count=excluded.imported_count,
                    skipped_count=excluded.skipped_count,
                    issue_count=excluded.issue_count,
                    content_hash=excluded.content_hash,
                    updated_at=excluded.updated_at,
                    completed_at=excluded.completed_at
                """,
                record,
            )
        return self.get_batch_import_run(batch_id)

    def update_batch_import_run(self, batch_id: str, payload: dict[str, Any]) -> dict[str, Any]:
        normalized_id = str(batch_id or "").strip()
        if not normalized_id:
            raise ValueError("batch_id 不能为空")
        current = self.get_batch_import_run(normalized_id)
        precheck = payload.get("precheck")
        result = payload.get("result")
        record = {
            **current,
            "status": str(payload.get("status") or current.get("status") or "").strip(),
            "precheck_json": json.dumps(precheck, ensure_ascii=False) if isinstance(precheck, dict) else current.get("precheck_json", "{}"),
            "result_json": json.dumps(result, ensure_ascii=False) if isinstance(result, dict) else current.get("result_json", "{}"),
            "imported_count": int(payload.get("imported_count") if payload.get("imported_count") is not None else current.get("imported_count") or 0),
            "skipped_count": int(payload.get("skipped_count") if payload.get("skipped_count") is not None else current.get("skipped_count") or 0),
            "issue_count": int(payload.get("issue_count") if payload.get("issue_count") is not None else current.get("issue_count") or 0),
            "updated_at": utc_now_iso(),
            "completed_at": str(payload.get("completed_at") if payload.get("completed_at") is not None else current.get("completed_at") or "").strip(),
        }
        with self.connect() as conn:
            conn.execute(
                """
                UPDATE batch_import_runs
                SET status=:status,
                    precheck_json=:precheck_json,
                    result_json=:result_json,
                    imported_count=:imported_count,
                    skipped_count=:skipped_count,
                    issue_count=:issue_count,
                    updated_at=:updated_at,
                    completed_at=:completed_at
                WHERE batch_id=:batch_id
                """,
                record,
            )
        return self.get_batch_import_run(normalized_id)

    def get_batch_import_run(self, batch_id: str) -> dict[str, Any]:
        with self.connect() as conn:
            row = conn.execute(
                f"""
                SELECT {BATCH_IMPORT_RUN_COLUMNS}
                FROM batch_import_runs
                WHERE batch_id = ?
                """,
                (str(batch_id or "").strip(),),
            ).fetchone()
        if row is None:
            raise KeyError(f"批量导入批次不存在: {batch_id}")
        return dict(row)

    def list_batch_import_runs(self, *, import_type: str = "", category: str = "", limit: int = 100) -> list[dict[str, Any]]:
        conditions = []
        params: list[Any] = []
        if import_type.strip():
            conditions.append("import_type = ?")
            params.append(import_type.strip())
        if category.strip():
            conditions.append("category = ?")
            params.append(category.strip())
        where_clause = f"WHERE {' AND '.join(conditions)}" if conditions else ""
        safe_limit = max(1, min(int(limit or 100), 500))
        with self.connect() as conn:
            rows = conn.execute(
                f"""
                SELECT {BATCH_IMPORT_RUN_COLUMNS}
                FROM batch_import_runs
                {where_clause}
                ORDER BY created_at DESC, id DESC
                LIMIT ?
                """,
                (*params, safe_limit),
            ).fetchall()
        return [dict(row) for row in rows]

    def list_sample_library(self) -> list[dict[str, Any]]:
        with self.connect() as conn:
            rows = conn.execute(
                """
                SELECT
                    """ + SAMPLE_LIBRARY_COLUMNS + """
                FROM sample_library
                ORDER BY updated_at DESC, sample_name ASC
                """
            ).fetchall()
        return [dict(row) for row in rows]

    def list_host_database(self, category: str = "host") -> list[dict[str, Any]]:
        with self.connect() as conn:
            rows = conn.execute(
                """
                SELECT
                    """ + HOST_DATABASE_COLUMNS + """
                FROM host_database
                WHERE db_category = ?
                ORDER BY updated_at DESC, genome_name ASC, host_name ASC
                """,
                (category,),
            ).fetchall()
        return [dict(row) for row in rows]

    def get_host_database_record(self, host_key: str) -> dict[str, Any]:
        with self.connect() as conn:
            row = conn.execute(
                """
                SELECT
                    """ + HOST_DATABASE_COLUMNS + """
                FROM host_database
                WHERE host_key = ?
                """,
                (host_key,),
            ).fetchone()
        if row is None:
            raise KeyError(f"宿主记录不存在: {host_key}")
        return dict(row)

    def upsert_host_database_record(self, record: dict[str, Any]) -> dict[str, Any]:
        host_key = str(record.get("host_key") or "").strip()
        if not host_key:
            raise ValueError("host_key 不能为空")
        now = utc_now_iso()
        payload = {
            "host_key": host_key,
            "db_category": str(record.get("db_category") or "host").strip() or "host",
            "host_name": str(record.get("host_name") or "").strip(),
            "genome_name": str(record.get("genome_name") or record.get("host_name") or "").strip(),
            "taxid": str(record.get("taxid") or "").strip(),
            "source_type": str(record.get("source_type") or "local").strip() or "local",
            "source_label": str(record.get("source_label") or "").strip(),
            "source_accession": str(record.get("source_accession") or "").strip(),
            "source_url": str(record.get("source_url") or "").strip(),
            "description": str(record.get("description") or "").strip(),
            "fasta_path": str(record.get("fasta_path") or "").strip(),
            "index_prefix": str(record.get("index_prefix") or "").strip(),
            "index_status": str(record.get("index_status") or "pending").strip() or "pending",
            "index_message": str(record.get("index_message") or "").strip(),
            "owner": str(record.get("owner") or "").strip(),
        }
        with self.connect() as conn:
            conn.execute(
                """
                INSERT INTO host_database (
                    host_key, db_category, host_name, genome_name, taxid, source_type, source_label, source_accession, source_url, description,
                    fasta_path, index_prefix, index_status, index_message, owner, created_at, updated_at
                )
                VALUES (
                    :host_key, :db_category, :host_name, :genome_name, :taxid, :source_type, :source_label, :source_accession, :source_url, :description,
                    :fasta_path, :index_prefix, :index_status, :index_message, :owner, :created_at, :updated_at
                )
                ON CONFLICT(host_key) DO UPDATE SET
                    db_category=excluded.db_category,
                    host_name=excluded.host_name,
                    genome_name=excluded.genome_name,
                    taxid=excluded.taxid,
                    source_type=excluded.source_type,
                    source_label=excluded.source_label,
                    source_accession=excluded.source_accession,
                    source_url=excluded.source_url,
                    description=excluded.description,
                    fasta_path=excluded.fasta_path,
                    index_prefix=excluded.index_prefix,
                    index_status=excluded.index_status,
                    index_message=excluded.index_message,
                    owner=excluded.owner,
                    updated_at=excluded.updated_at
                """,
                {
                    **payload,
                    "created_at": now,
                    "updated_at": now,
                },
            )
        return self.get_host_database_record(host_key)

    def update_host_database_index(self, host_key: str, *, index_prefix: str | None = None, index_status: str | None = None, index_message: str | None = None) -> dict[str, Any]:
        current = self.get_host_database_record(host_key)
        next_values = {
            "index_prefix": current.get("index_prefix", "") if index_prefix is None else str(index_prefix).strip(),
            "index_status": current.get("index_status", "pending") if index_status is None else (str(index_status).strip() or "pending"),
            "index_message": current.get("index_message", "") if index_message is None else str(index_message).strip(),
        }
        with self.connect() as conn:
            conn.execute(
                """
                UPDATE host_database
                SET index_prefix = ?, index_status = ?, index_message = ?, updated_at = ?
                WHERE host_key = ?
                """,
                (
                    next_values["index_prefix"],
                    next_values["index_status"],
                    next_values["index_message"],
                    utc_now_iso(),
                    host_key,
                ),
            )
        return self.get_host_database_record(host_key)

    def update_host_database_record(
        self,
        host_key: str,
        *,
        host_name: str | None = None,
        genome_name: str | None = None,
        taxid: str | None = None,
        source_label: str | None = None,
        source_accession: str | None = None,
        source_url: str | None = None,
        description: str | None = None,
        fasta_path: str | None = None,
        index_prefix: str | None = None,
        index_status: str | None = None,
        index_message: str | None = None,
    ) -> dict[str, Any]:
        current = self.get_host_database_record(host_key)
        next_values = {
            "host_name": current.get("host_name", "") if host_name is None else str(host_name).strip(),
            "genome_name": current.get("genome_name", "") if genome_name is None else str(genome_name).strip(),
            "taxid": current.get("taxid", "") if taxid is None else str(taxid).strip(),
            "source_label": current.get("source_label", "") if source_label is None else str(source_label).strip(),
            "source_accession": current.get("source_accession", "") if source_accession is None else str(source_accession).strip(),
            "source_url": current.get("source_url", "") if source_url is None else str(source_url).strip(),
            "description": current.get("description", "") if description is None else str(description).strip(),
            "fasta_path": current.get("fasta_path", "") if fasta_path is None else str(fasta_path).strip(),
            "index_prefix": current.get("index_prefix", "") if index_prefix is None else str(index_prefix).strip(),
            "index_status": current.get("index_status", "pending") if index_status is None else (str(index_status).strip() or "pending"),
            "index_message": current.get("index_message", "") if index_message is None else str(index_message).strip(),
        }
        with self.connect() as conn:
            conn.execute(
                """
                UPDATE host_database
                SET host_name = ?, genome_name = ?, taxid = ?, source_label = ?, source_accession = ?, source_url = ?, description = ?,
                    fasta_path = ?, index_prefix = ?, index_status = ?, index_message = ?, updated_at = ?
                WHERE host_key = ?
                """,
                (
                    next_values["host_name"],
                    next_values["genome_name"],
                    next_values["taxid"],
                    next_values["source_label"],
                    next_values["source_accession"],
                    next_values["source_url"],
                    next_values["description"],
                    next_values["fasta_path"],
                    next_values["index_prefix"],
                    next_values["index_status"],
                    next_values["index_message"],
                    utc_now_iso(),
                    host_key,
                ),
            )
        return self.get_host_database_record(host_key)

    def delete_host_database_record(self, host_key: str) -> None:
        with self.connect() as conn:
            conn.execute("DELETE FROM host_database WHERE host_key = ?", (host_key,))

    def list_reference_panels(self, category: str = "pathogen", panel_type: str = "cgmlst") -> list[dict[str, Any]]:
        with self.connect() as conn:
            rows = conn.execute(
                """
                SELECT
                    """ + REFERENCE_PANEL_COLUMNS + """
                FROM reference_panels
                WHERE db_category = ? AND panel_type = ?
                ORDER BY updated_at DESC, species_name ASC
                """,
                (category, panel_type),
            ).fetchall()
        return [dict(row) for row in rows]

    def get_reference_panel(self, panel_key: str) -> dict[str, Any]:
        with self.connect() as conn:
            row = conn.execute(
                """
                SELECT
                    """ + REFERENCE_PANEL_COLUMNS + """
                FROM reference_panels
                WHERE panel_key = ?
                """,
                (panel_key,),
            ).fetchone()
        if row is None:
            raise KeyError(f"panel 不存在: {panel_key}")
        return dict(row)

    def upsert_reference_panel(self, record: dict[str, Any]) -> dict[str, Any]:
        panel_key = str(record.get("panel_key") or "").strip()
        if not panel_key:
            raise ValueError("panel_key 不能为空")
        now = utc_now_iso()
        selected_host_keys = record.get("selected_host_keys_json")
        if isinstance(selected_host_keys, (list, tuple)):
            selected_host_keys = json.dumps([str(item or "").strip() for item in selected_host_keys if str(item or "").strip()], ensure_ascii=False)
        payload = {
            "panel_key": panel_key,
            "db_category": str(record.get("db_category") or "pathogen").strip() or "pathogen",
            "panel_type": str(record.get("panel_type") or "cgmlst").strip() or "cgmlst",
            "species_name": str(record.get("species_name") or "").strip(),
            "species_slug": str(record.get("species_slug") or "").strip(),
            "selected_host_keys_json": str(selected_host_keys or "[]").strip() or "[]",
            "genome_count": int(record.get("genome_count") or 0),
            "schema_dir": str(record.get("schema_dir") or "").strip(),
            "cgmlst_dir": str(record.get("cgmlst_dir") or "").strip(),
            "loci_list_path": str(record.get("loci_list_path") or "").strip(),
            "results_alleles_path": str(record.get("results_alleles_path") or "").strip(),
            "threshold": str(record.get("threshold") or "0.95").strip() or "0.95",
            "status": str(record.get("status") or "pending").strip() or "pending",
            "message": str(record.get("message") or "").strip(),
            "owner": str(record.get("owner") or "").strip(),
        }
        with self.connect() as conn:
            conn.execute(
                """
                INSERT INTO reference_panels (
                    panel_key, db_category, panel_type, species_name, species_slug, selected_host_keys_json, genome_count,
                    schema_dir, cgmlst_dir, loci_list_path, results_alleles_path, threshold, status, message, owner, created_at, updated_at
                )
                VALUES (
                    :panel_key, :db_category, :panel_type, :species_name, :species_slug, :selected_host_keys_json, :genome_count,
                    :schema_dir, :cgmlst_dir, :loci_list_path, :results_alleles_path, :threshold, :status, :message, :owner, :created_at, :updated_at
                )
                ON CONFLICT(panel_key) DO UPDATE SET
                    db_category=excluded.db_category,
                    panel_type=excluded.panel_type,
                    species_name=excluded.species_name,
                    species_slug=excluded.species_slug,
                    selected_host_keys_json=excluded.selected_host_keys_json,
                    genome_count=excluded.genome_count,
                    schema_dir=excluded.schema_dir,
                    cgmlst_dir=excluded.cgmlst_dir,
                    loci_list_path=excluded.loci_list_path,
                    results_alleles_path=excluded.results_alleles_path,
                    threshold=excluded.threshold,
                    status=excluded.status,
                    message=excluded.message,
                    owner=excluded.owner,
                    updated_at=excluded.updated_at
                """,
                {
                    **payload,
                    "created_at": now,
                    "updated_at": now,
                },
            )
        return self.get_reference_panel(panel_key)

    def list_sample_library_by_scope(self, scope: str) -> list[dict[str, Any]]:
        with self.connect() as conn:
            rows = conn.execute(
                """
                SELECT
                    """ + SAMPLE_LIBRARY_COLUMNS + """
                FROM sample_library
                WHERE library_scope = ?
                ORDER BY updated_at DESC, sample_name ASC
                """,
                (scope,),
            ).fetchall()
        return [dict(row) for row in rows]

    def list_sample_library_metadata_templates(self) -> list[dict[str, Any]]:
        with self.connect() as conn:
            rows = conn.execute(
                """
                SELECT field_key, label, field_type, options_json, config_json, position, updated_at
                FROM sample_library_metadata_templates
                ORDER BY position ASC, id ASC
                """
            ).fetchall()
        items: list[dict[str, Any]] = []
        for row in rows:
            raw_options = str(row["options_json"] or "").strip()
            raw_config = str(row["config_json"] or "").strip()
            try:
                options = json.loads(raw_options) if raw_options else []
            except Exception:
                options = []
            try:
                config = json.loads(raw_config) if raw_config else {}
            except Exception:
                config = {}
            if not isinstance(config, dict):
                config = {}
            items.append(
                {
                    "key": row["field_key"],
                    "label": row["label"] or row["field_key"],
                    "type": row["field_type"] or "text",
                    "options": options if isinstance(options, list) else [],
                    "group": str(config.get("group") or "").strip(),
                    "requirement": str(config.get("requirement") or "").strip(),
                    "placeholder": str(config.get("placeholder") or "").strip(),
                    "help_text": str(config.get("help_text") or "").strip(),
                    "dictionary_name": str(config.get("dictionary_name") or "").strip(),
                    "position": row["position"],
                    "updated_at": row["updated_at"],
                }
            )
        return items

    def upsert_sample_library_metadata_templates(self, items: list[dict[str, Any]]) -> list[dict[str, Any]]:
        normalized: list[dict[str, Any]] = []
        for position, item in enumerate(items):
            field_key = str(item.get("key") or "").strip()
            if not field_key:
                continue
            field_type = str(item.get("type") or "text").strip() or "text"
            options = [str(option).strip() for option in (item.get("options") or []) if str(option).strip()] if isinstance(item.get("options"), list) else []
            normalized.append(
                {
                    "field_key": field_key,
                    "label": str(item.get("label") or field_key).strip(),
                    "field_type": field_type,
                    "options_json": json.dumps(options, ensure_ascii=False) if options else None,
                    "config_json": json.dumps(
                        {
                            "group": str(item.get("group") or "").strip(),
                            "requirement": str(item.get("requirement") or "").strip(),
                            "placeholder": str(item.get("placeholder") or "").strip(),
                            "help_text": str(item.get("help_text") or "").strip(),
                            "dictionary_name": str(item.get("dictionary_name") or "").strip(),
                        },
                        ensure_ascii=False,
                    ),
                    "position": position,
                    "updated_at": utc_now_iso(),
                }
            )
        with self.connect() as conn:
            keep_keys = [payload["field_key"] for payload in normalized]
            if keep_keys:
                placeholders = ",".join(["?"] * len(keep_keys))
                conn.execute(
                    f"DELETE FROM sample_library_metadata_templates WHERE field_key NOT IN ({placeholders})",
                    keep_keys,
                )
            else:
                conn.execute("DELETE FROM sample_library_metadata_templates")
            for payload in normalized:
                conn.execute(
                    """
                    INSERT INTO sample_library_metadata_templates (field_key, label, field_type, options_json, config_json, position, updated_at)
                    VALUES (:field_key, :label, :field_type, :options_json, :config_json, :position, :updated_at)
                    ON CONFLICT(field_key) DO UPDATE SET
                        label=excluded.label,
                        field_type=excluded.field_type,
                        options_json=excluded.options_json,
                        config_json=excluded.config_json,
                        position=excluded.position,
                        updated_at=excluded.updated_at
                    """,
                    payload,
                )
        return self.list_sample_library_metadata_templates()

    def get_sample_library_record(self, sample_key: str) -> dict[str, Any]:
        with self.connect() as conn:
            row = conn.execute(
                """
                SELECT
                    """ + SAMPLE_LIBRARY_COLUMNS + """
                FROM sample_library
                WHERE sample_key = ?
                """,
                (sample_key,),
            ).fetchone()
        if row is None:
            raise KeyError(f"样本不存在: {sample_key}")
        return dict(row)

    def upsert_sample_library_record(self, record: dict[str, Any]) -> dict[str, Any]:
        sample_key = str(record.get("sample_key") or "").strip()
        if not sample_key:
            raise ValueError("sample_key 不能为空")
        now = utc_now_iso()
        payload = {
            "sample_key": sample_key,
            "genome_id": str(record.get("genome_id") or "").strip(),
            "sample_name": str(record.get("sample_name") or "").strip(),
            "task_id": str(record.get("task_id") or "").strip(),
            "task_name": str(record.get("task_name") or "").strip(),
            "owner": str(record.get("owner") or "").strip(),
            "owner_group": str(record.get("owner_group") or "").strip(),
            "report_dir": str(record.get("report_dir") or "").strip(),
            "output_dir": str(record.get("output_dir") or "").strip(),
            "final_fasta_path": str(record.get("final_fasta_path") or "").strip(),
            "species_name": str(record.get("species_name") or "").strip(),
            "pathogen_type": str(record.get("pathogen_type") or "").strip(),
            "taxid": str(record.get("taxid") or "").strip(),
            "mlst_species_name": str(record.get("mlst_species_name") or "").strip(),
            "mlst_st": str(record.get("mlst_st") or "").strip(),
            "serotype_result": str(record.get("serotype_result") or "").strip(),
            "genome_length": str(record.get("genome_length") or "").strip(),
            "q20_rate": str(record.get("q20_rate") or "").strip(),
            "q30_rate": str(record.get("q30_rate") or "").strip(),
            "completeness": str(record.get("completeness") or "").strip(),
            "contamination": str(record.get("contamination") or "").strip(),
            "contig_count": str(record.get("contig_count") or "").strip(),
            "plasmid_count": str(record.get("plasmid_count") or "").strip(),
            "total_length": str(record.get("total_length") or "").strip(),
            "resistance_count": str(record.get("resistance_count") or "").strip(),
            "virulence_count": str(record.get("virulence_count") or "").strip(),
            "resistance_gene_hits": str(record.get("resistance_gene_hits") or "").strip(),
            "virulence_gene_hits": str(record.get("virulence_gene_hits") or "").strip(),
            "resistance_mge_hits": str(record.get("resistance_mge_hits") or "").strip(),
            "virulence_mge_hits": str(record.get("virulence_mge_hits") or "").strip(),
            "description": str(record.get("description") or "").strip(),
            "gender": str(record.get("gender") or "").strip(),
            "country": str(record.get("country") or "").strip(),
            "location_json": str(record.get("location_json") or "").strip(),
            "sample_type": str(record.get("sample_type") or "").strip(),
            "sequencing_method": str(record.get("sequencing_method") or "").strip(),
            "custom_metadata_json": str(record.get("custom_metadata_json") or "[]").strip() or "[]",
            "sample_alias": str(record.get("sample_alias") or "").strip(),
            "sample_source": str(record.get("sample_source") or "").strip(),
            "collection_date": str(record.get("collection_date") or "").strip(),
            "host_info": str(record.get("host_info") or "").strip(),
            "note": str(record.get("note") or "").strip(),
            "library_scope": str(record.get("library_scope") or "main").strip() or "main",
            "visibility_scope": str(record.get("visibility_scope") or "group").strip() or "group",
            "source_submission_id": str(record.get("source_submission_id") or "").strip(),
        }
        with self.connect() as conn:
            conn.execute(
                """
                INSERT INTO sample_library (
                    sample_key, genome_id, sample_name, task_id, task_name, owner, owner_group, report_dir, output_dir,
                    final_fasta_path, species_name, pathogen_type, taxid, mlst_species_name, mlst_st, serotype_result, genome_length, q20_rate, q30_rate, completeness,
                    contamination, contig_count, plasmid_count, total_length, resistance_count,
                    virulence_count, resistance_gene_hits, virulence_gene_hits, resistance_mge_hits, virulence_mge_hits, description, gender, country, location_json, sample_type, sequencing_method,
                    custom_metadata_json, sample_alias, sample_source, collection_date, host_info, note, library_scope, visibility_scope, source_submission_id,
                    imported_at, updated_at
                )
                VALUES (
                    :sample_key, :genome_id, :sample_name, :task_id, :task_name, :owner, :owner_group, :report_dir, :output_dir,
                    :final_fasta_path, :species_name, :pathogen_type, :taxid, :mlst_species_name, :mlst_st, :serotype_result, :genome_length, :q20_rate, :q30_rate, :completeness,
                    :contamination, :contig_count, :plasmid_count, :total_length, :resistance_count,
                    :virulence_count, :resistance_gene_hits, :virulence_gene_hits, :resistance_mge_hits, :virulence_mge_hits, :description, :gender, :country, :location_json, :sample_type, :sequencing_method,
                    :custom_metadata_json, :sample_alias, :sample_source, :collection_date, :host_info, :note, :library_scope, :visibility_scope, :source_submission_id,
                    :imported_at, :updated_at
                )
                ON CONFLICT(sample_key) DO UPDATE SET
                    genome_id=excluded.genome_id,
                    sample_name=excluded.sample_name,
                    task_id=excluded.task_id,
                    task_name=excluded.task_name,
                    owner=excluded.owner,
                    owner_group=excluded.owner_group,
                    report_dir=excluded.report_dir,
                    output_dir=excluded.output_dir,
                    final_fasta_path=excluded.final_fasta_path,
                    species_name=excluded.species_name,
                    pathogen_type=excluded.pathogen_type,
                    taxid=excluded.taxid,
                    mlst_species_name=excluded.mlst_species_name,
                    mlst_st=excluded.mlst_st,
                    serotype_result=excluded.serotype_result,
                    genome_length=excluded.genome_length,
                    q20_rate=excluded.q20_rate,
                    q30_rate=excluded.q30_rate,
                    completeness=excluded.completeness,
                    contamination=excluded.contamination,
                    contig_count=excluded.contig_count,
                    plasmid_count=excluded.plasmid_count,
                    total_length=excluded.total_length,
                    resistance_count=excluded.resistance_count,
                    virulence_count=excluded.virulence_count,
                    resistance_gene_hits=excluded.resistance_gene_hits,
                    virulence_gene_hits=excluded.virulence_gene_hits,
                    resistance_mge_hits=excluded.resistance_mge_hits,
                    virulence_mge_hits=excluded.virulence_mge_hits,
                    description=excluded.description,
                    gender=excluded.gender,
                    country=excluded.country,
                    location_json=excluded.location_json,
                    sample_type=excluded.sample_type,
                    sequencing_method=excluded.sequencing_method,
                    custom_metadata_json=excluded.custom_metadata_json,
                    sample_alias=excluded.sample_alias,
                    sample_source=excluded.sample_source,
                    collection_date=excluded.collection_date,
                    host_info=excluded.host_info,
                    note=excluded.note,
                    library_scope=excluded.library_scope,
                    visibility_scope=excluded.visibility_scope,
                    source_submission_id=excluded.source_submission_id,
                    updated_at=excluded.updated_at
                """,
                {
                    **payload,
                    "imported_at": now,
                    "updated_at": now,
                },
            )
            row = conn.execute(
                """
                SELECT
                    """ + SAMPLE_LIBRARY_COLUMNS + """
                FROM sample_library
                WHERE sample_key = ?
                """,
                (sample_key,),
            ).fetchone()
        return dict(row) if row else payload

    def update_sample_library_record(
        self,
        sample_key: str,
        *,
        genome_id: str | None = None,
        pathogen_type: str | None = None,
        sample_alias: str | None = None,
        taxid: str | None = None,
        mlst_st: str | None = None,
        serotype_result: str | None = None,
        resistance_gene_hits: str | None = None,
        virulence_gene_hits: str | None = None,
        resistance_mge_hits: str | None = None,
        virulence_mge_hits: str | None = None,
        description: str | None = None,
        sample_source: str | None = None,
        collection_date: str | None = None,
        gender: str | None = None,
        country: str | None = None,
        host_info: str | None = None,
        location_json: str | None = None,
        sample_type: str | None = None,
        sequencing_method: str | None = None,
        genome_length: str | None = None,
        note: str | None = None,
        visibility_scope: str | None = None,
        custom_metadata_json: str | None = None,
    ) -> dict[str, Any]:
        current = self.get_sample_library_record(sample_key)
        next_values = {
            "genome_id": current.get("genome_id", "") if genome_id is None else str(genome_id).strip(),
            "pathogen_type": current.get("pathogen_type", "") if pathogen_type is None else str(pathogen_type).strip(),
            "sample_alias": current.get("sample_alias", "") if sample_alias is None else str(sample_alias).strip(),
            "taxid": current.get("taxid", "") if taxid is None else str(taxid).strip(),
            "mlst_st": current.get("mlst_st", "") if mlst_st is None else str(mlst_st).strip(),
            "serotype_result": current.get("serotype_result", "") if serotype_result is None else str(serotype_result).strip(),
            "resistance_gene_hits": current.get("resistance_gene_hits", "") if resistance_gene_hits is None else str(resistance_gene_hits).strip(),
            "virulence_gene_hits": current.get("virulence_gene_hits", "") if virulence_gene_hits is None else str(virulence_gene_hits).strip(),
            "resistance_mge_hits": current.get("resistance_mge_hits", "") if resistance_mge_hits is None else str(resistance_mge_hits).strip(),
            "virulence_mge_hits": current.get("virulence_mge_hits", "") if virulence_mge_hits is None else str(virulence_mge_hits).strip(),
            "description": current.get("description", "") if description is None else str(description).strip(),
            "sample_source": current.get("sample_source", "") if sample_source is None else str(sample_source).strip(),
            "collection_date": current.get("collection_date", "") if collection_date is None else str(collection_date).strip(),
            "gender": current.get("gender", "") if gender is None else str(gender).strip(),
            "country": current.get("country", "") if country is None else str(country).strip(),
            "host_info": current.get("host_info", "") if host_info is None else str(host_info).strip(),
            "location_json": current.get("location_json", "") if location_json is None else str(location_json).strip(),
            "sample_type": current.get("sample_type", "") if sample_type is None else str(sample_type).strip(),
            "sequencing_method": current.get("sequencing_method", "") if sequencing_method is None else str(sequencing_method).strip(),
            "genome_length": current.get("genome_length", "") if genome_length is None else str(genome_length).strip(),
            "note": current.get("note", "") if note is None else str(note).strip(),
            "visibility_scope": current.get("visibility_scope", "group") if visibility_scope is None else (str(visibility_scope).strip() or "group"),
            "custom_metadata_json": current.get("custom_metadata_json", "[]") if custom_metadata_json is None else (str(custom_metadata_json).strip() or "[]"),
        }
        with self.connect() as conn:
            conn.execute(
                """
                UPDATE sample_library
                SET genome_id = ?, pathogen_type = ?, sample_alias = ?, taxid = ?, mlst_st = ?, serotype_result = ?, resistance_gene_hits = ?, virulence_gene_hits = ?, resistance_mge_hits = ?, virulence_mge_hits = ?, description = ?, sample_source = ?, collection_date = ?,
                    gender = ?, country = ?, host_info = ?, location_json = ?, sample_type = ?, sequencing_method = ?,
                    genome_length = ?, note = ?, visibility_scope = ?, custom_metadata_json = ?, updated_at = ?
                WHERE sample_key = ?
                """,
                (
                    next_values["genome_id"],
                    next_values["pathogen_type"],
                    next_values["sample_alias"],
                    next_values["taxid"],
                    next_values["mlst_st"],
                    next_values["serotype_result"],
                    next_values["resistance_gene_hits"],
                    next_values["virulence_gene_hits"],
                    next_values["resistance_mge_hits"],
                    next_values["virulence_mge_hits"],
                    next_values["description"],
                    next_values["sample_source"],
                    next_values["collection_date"],
                    next_values["gender"],
                    next_values["country"],
                    next_values["host_info"],
                    next_values["location_json"],
                    next_values["sample_type"],
                    next_values["sequencing_method"],
                    next_values["genome_length"],
                    next_values["note"],
                    next_values["visibility_scope"],
                    next_values["custom_metadata_json"],
                    utc_now_iso(),
                    sample_key,
                ),
            )
        return self.get_sample_library_record(sample_key)

    def delete_sample_library_record(self, sample_key: str) -> None:
        with self.connect() as conn:
            row = conn.execute("SELECT sample_key FROM sample_library WHERE sample_key = ?", (sample_key,)).fetchone()
            if row is None:
                raise KeyError(f"样本不存在: {sample_key}")
            conn.execute("DELETE FROM sample_library WHERE sample_key = ?", (sample_key,))

    def list_sample_library_submissions(self, status: str | None = None) -> list[dict[str, Any]]:
        query = """
            SELECT
                """ + SUBMISSION_COLUMNS + """
            FROM sample_library_submissions
        """
        params: list[Any] = []
        if status:
            query += " WHERE status = ?"
            params.append(status)
        query += " ORDER BY created_at DESC"
        with self.connect() as conn:
            rows = conn.execute(query, params).fetchall()
        return [dict(row) for row in rows]

    def get_sample_library_submission(self, request_id: str) -> dict[str, Any]:
        with self.connect() as conn:
            row = conn.execute(
                """
                SELECT
                    """ + SUBMISSION_COLUMNS + """
                FROM sample_library_submissions
                WHERE request_id = ?
                """,
                (request_id,),
            ).fetchone()
        if row is None:
            raise KeyError(f"审核申请不存在: {request_id}")
        return dict(row)

    def find_pending_submission_for_sample(self, personal_sample_key: str) -> dict[str, Any] | None:
        with self.connect() as conn:
            row = conn.execute(
                """
                SELECT
                    """ + SUBMISSION_COLUMNS + """
                FROM sample_library_submissions
                WHERE personal_sample_key = ? AND status = 'pending'
                ORDER BY created_at DESC
                LIMIT 1
                """,
                (personal_sample_key,),
            ).fetchone()
        return dict(row) if row else None

    def upsert_sample_library_submission(self, submission: dict[str, Any]) -> dict[str, Any]:
        request_id = str(submission.get("request_id") or "").strip()
        if not request_id:
            raise ValueError("request_id 不能为空")
        payload = {
            "request_id": request_id,
            "personal_sample_key": str(submission.get("personal_sample_key") or "").strip(),
            "sample_name": str(submission.get("sample_name") or "").strip(),
            "owner": str(submission.get("owner") or "").strip(),
            "owner_group": str(submission.get("owner_group") or "").strip(),
            "payload_json": str(submission.get("payload_json") or "{}").strip() or "{}",
            "status": str(submission.get("status") or "pending").strip() or "pending",
            "review_note": str(submission.get("review_note") or "").strip(),
            "reviewed_by": str(submission.get("reviewed_by") or "").strip(),
            "reviewed_at": str(submission.get("reviewed_at") or "").strip(),
            "created_at": str(submission.get("created_at") or utc_now_iso()),
            "updated_at": utc_now_iso(),
        }
        with self.connect() as conn:
            conn.execute(
                """
                INSERT INTO sample_library_submissions (
                    request_id, personal_sample_key, sample_name, owner, owner_group, payload_json,
                    status, review_note, reviewed_by, reviewed_at, created_at, updated_at
                )
                VALUES (
                    :request_id, :personal_sample_key, :sample_name, :owner, :owner_group, :payload_json,
                    :status, :review_note, :reviewed_by, :reviewed_at, :created_at, :updated_at
                )
                ON CONFLICT(request_id) DO UPDATE SET
                    personal_sample_key=excluded.personal_sample_key,
                    sample_name=excluded.sample_name,
                    owner=excluded.owner,
                    owner_group=excluded.owner_group,
                    payload_json=excluded.payload_json,
                    status=excluded.status,
                    review_note=excluded.review_note,
                    reviewed_by=excluded.reviewed_by,
                    reviewed_at=excluded.reviewed_at,
                    updated_at=excluded.updated_at
                """,
                payload,
            )
        return self.get_sample_library_submission(request_id)

    def list_sample_library_version_logs(self, limit: int = 50) -> list[dict[str, Any]]:
        safe_limit = max(1, min(int(limit or 50), 200))
        with self.connect() as conn:
            rows = conn.execute(
                f"""
                SELECT
                    {VERSION_LOG_COLUMNS}
                FROM sample_library_version_logs
                ORDER BY created_at DESC
                LIMIT ?
                """,
                (safe_limit,),
            ).fetchall()
        return [dict(row) for row in rows]

    def create_sample_library_version_log(self, payload: dict[str, Any]) -> dict[str, Any]:
        event_id = str(payload.get("event_id") or "").strip()
        if not event_id:
            raise ValueError("event_id 不能为空")
        record = {
            "event_id": event_id,
            "version_label": str(payload.get("version_label") or "").strip(),
            "request_id": str(payload.get("request_id") or "").strip(),
            "sample_key": str(payload.get("sample_key") or "").strip(),
            "sample_name": str(payload.get("sample_name") or "").strip(),
            "owner": str(payload.get("owner") or "").strip(),
            "owner_group": str(payload.get("owner_group") or "").strip(),
            "action": str(payload.get("action") or "publish").strip() or "publish",
            "summary": str(payload.get("summary") or "").strip(),
            "operator": str(payload.get("operator") or "").strip(),
            "payload_json": str(payload.get("payload_json") or "{}").strip() or "{}",
            "created_at": str(payload.get("created_at") or utc_now_iso()),
        }
        with self.connect() as conn:
            conn.execute(
                """
                INSERT INTO sample_library_version_logs (
                    event_id, version_label, request_id, sample_key, sample_name, owner, owner_group,
                    action, summary, operator, payload_json, created_at
                )
                VALUES (
                    :event_id, :version_label, :request_id, :sample_key, :sample_name, :owner, :owner_group,
                    :action, :summary, :operator, :payload_json, :created_at
                )
                """,
                record,
            )
        return record

    def list_sample_library_release_versions(self, scope: str = "main", limit: int = 30) -> list[dict[str, Any]]:
        safe_limit = max(1, min(int(limit or 30), 200))
        with self.connect() as conn:
            rows = conn.execute(
                f"""
                SELECT
                    {RELEASE_VERSION_COLUMNS}
                FROM sample_library_release_versions
                WHERE scope = ?
                ORDER BY created_at DESC
                LIMIT ?
                """,
                (scope, safe_limit),
            ).fetchall()
        return [dict(row) for row in rows]

    def create_sample_library_release_version(self, payload: dict[str, Any]) -> dict[str, Any]:
        release_id = str(payload.get("release_id") or "").strip()
        if not release_id:
            raise ValueError("release_id 不能为空")
        record = {
            "release_id": release_id,
            "version_label": str(payload.get("version_label") or "").strip(),
            "scope": str(payload.get("scope") or "main").strip() or "main",
            "summary": str(payload.get("summary") or "").strip(),
            "note": str(payload.get("note") or "").strip(),
            "operator": str(payload.get("operator") or "").strip(),
            "change_count": int(payload.get("change_count") or 0),
            "added_count": int(payload.get("added_count") or 0),
            "updated_count": int(payload.get("updated_count") or 0),
            "deleted_count": int(payload.get("deleted_count") or 0),
            "changes_json": str(payload.get("changes_json") or "{}").strip() or "{}",
            "snapshot_json": str(payload.get("snapshot_json") or "[]").strip() or "[]",
            "created_at": str(payload.get("created_at") or utc_now_iso()),
        }
        with self.connect() as conn:
            conn.execute(
                """
                INSERT INTO sample_library_release_versions (
                    release_id, version_label, scope, summary, note, operator, change_count,
                    added_count, updated_count, deleted_count, changes_json, snapshot_json, created_at
                )
                VALUES (
                    :release_id, :version_label, :scope, :summary, :note, :operator, :change_count,
                    :added_count, :updated_count, :deleted_count, :changes_json, :snapshot_json, :created_at
                )
                """,
                record,
            )
        return record

    def _serialize_user_row(self, row: sqlite3.Row) -> dict[str, Any]:
        role = self._normalize_role(row["role"])
        granted_modules = self._normalize_allowed_modules(row["allowed_modules"], role=role, default_if_empty=DEFAULT_ALLOWED_MODULES)
        module_expirations = self._normalize_module_expirations(
            row["module_expirations"] if "module_expirations" in row.keys() else "",
            granted_modules=granted_modules,
            role=role,
        )
        account_expires_at = self._normalize_account_expires_at(row["account_expires_at"] if "account_expires_at" in row.keys() else "")
        is_expired = self._is_account_expired(account_expires_at)
        allowed_modules = [] if is_expired else self._filter_active_modules(granted_modules, module_expirations, role=role)
        granted_viruses = self._normalize_allowed_viruses(
            row["allowed_viruses"] if "allowed_viruses" in row.keys() else "",
            default_if_empty=DEFAULT_ALLOWED_VIRUSES,
        )
        allowed_viruses = [] if is_expired else list(granted_viruses)
        return {
            "username": row["username"],
            "role": role,
            "group_name": row["group_name"] or "",
            "allowed_modules": allowed_modules,
            "granted_modules": granted_modules,
            "allowed_viruses": allowed_viruses,
            "granted_viruses": granted_viruses,
            "module_expirations": module_expirations,
            "module_access": [
                {
                    "module": module,
                    "expires_at": module_expirations.get(module, ""),
                    "is_expired": self._is_module_expired(module_expirations.get(module, "")),
                    "enabled": module in allowed_modules,
                }
                for module in granted_modules
            ],
            "account_expires_at": account_expires_at,
            "is_expired": is_expired,
            "display_name": row["display_name"] or "",
            "created_at": row["created_at"],
            "updated_at": row["updated_at"],
        }

    def _normalize_role(self, role: str) -> str:
        text = str(role or "").strip()
        return "group_admin" if text == "group" else text

    def _normalize_allowed_modules(self, value: Any, *, role: str = "", default_if_empty: Any = None) -> list[str]:
        normalized_role = self._normalize_role(role)
        if normalized_role == "admin":
            return list(DEFAULT_ALLOWED_MODULES)
        modules: list[str] = []

        def _append_module(raw_module: Any) -> None:
            key = str(raw_module or "").strip().lower()
            if not key:
                return
            if key == "pathogen":
                for legacy_key in ("bacteria", "virus", "metagenome"):
                    if legacy_key not in modules:
                        modules.append(legacy_key)
                return
            if key == "single":
                for legacy_key in ("bacteria", "virus"):
                    if legacy_key not in modules:
                        modules.append(legacy_key)
                return
            if key in DEFAULT_ALLOWED_MODULES and key not in modules:
                modules.append(key)

        if isinstance(value, str):
            text = value.strip()
            if text:
                try:
                    parsed = json.loads(text)
                except Exception:
                    parsed = [item.strip() for item in text.split(",") if item.strip()]
            else:
                parsed = default_if_empty
        else:
            parsed = value if value is not None else default_if_empty
        if isinstance(parsed, (list, tuple, set)):
            for item in parsed:
                _append_module(item)
        if not modules:
            fallback = default_if_empty if default_if_empty is not None else ["bacteria"]
            if isinstance(fallback, str):
                try:
                    fallback = json.loads(fallback)
                except Exception:
                    fallback = [item.strip() for item in fallback.split(",") if item.strip()]
            if isinstance(fallback, (list, tuple, set)):
                for item in fallback:
                    _append_module(item)
        if not modules:
            modules = ["bacteria"]
        return modules

    def _normalize_allowed_viruses(self, value: Any, *, default_if_empty: Any = None) -> list[str]:
        viruses: list[str] = []

        def _append_virus(raw_virus: Any) -> None:
            key = str(raw_virus or "").strip().lower()
            if key in DEFAULT_ALLOWED_VIRUSES and key not in viruses:
                viruses.append(key)

        if isinstance(value, str):
            text = value.strip()
            if text:
                try:
                    parsed = json.loads(text)
                except Exception:
                    parsed = [item.strip() for item in text.split(",") if item.strip()]
            else:
                parsed = default_if_empty
        else:
            parsed = value if value is not None else default_if_empty
        if isinstance(parsed, (list, tuple, set)):
            for item in parsed:
                _append_virus(item)
        elif parsed is not None:
            _append_virus(parsed)
        if not viruses:
            fallback = default_if_empty if default_if_empty is not None else DEFAULT_ALLOWED_VIRUSES
            if isinstance(fallback, str):
                try:
                    fallback = json.loads(fallback)
                except Exception:
                    fallback = [item.strip() for item in fallback.split(",") if item.strip()]
            if isinstance(fallback, (list, tuple, set)):
                for item in fallback:
                    _append_virus(item)
        if not viruses:
            viruses = list(DEFAULT_ALLOWED_VIRUSES)
        return viruses

    def _normalize_account_expires_at(self, value: Any) -> str:
        text = str(value or "").strip()
        if not text:
            return ""
        normalized = text.replace(" ", "T")
        if len(normalized) == 10:
            normalized = f"{normalized}T23:59:59"
        try:
            parsed = datetime.fromisoformat(normalized.replace("Z", "+00:00"))
        except ValueError as exc:
            raise ValueError("账号有效期格式不正确，请使用 YYYY-MM-DD 或 YYYY-MM-DDTHH:MM") from exc
        if parsed.tzinfo is None:
            parsed = parsed.replace(tzinfo=timezone.utc)
        else:
            parsed = parsed.astimezone(timezone.utc)
        return parsed.isoformat()

    def _is_account_expired(self, value: Any) -> bool:
        expires_at = self._normalize_account_expires_at(value)
        if not expires_at:
            return False
        parsed = datetime.fromisoformat(expires_at.replace("Z", "+00:00"))
        return parsed <= datetime.now(timezone.utc)

    def _normalize_module_expirations(self, value: Any, *, granted_modules: list[str], role: str = "") -> dict[str, str]:
        normalized_role = self._normalize_role(role)
        if normalized_role == "admin":
            return {}
        parsed: Any = value
        if isinstance(value, str):
            text = value.strip()
            if not text:
                parsed = {}
            else:
                try:
                    parsed = json.loads(text)
                except Exception:
                    parsed = {}
        expirations: dict[str, str] = {}
        if isinstance(parsed, dict):
            for module in granted_modules:
                if module not in DEFAULT_ALLOWED_MODULES:
                    continue
                if module in parsed:
                    expirations[module] = self._normalize_account_expires_at(parsed.get(module))
        return expirations

    def _is_module_expired(self, value: Any) -> bool:
        expires_at = self._normalize_account_expires_at(value)
        if not expires_at:
            return False
        parsed = datetime.fromisoformat(expires_at.replace("Z", "+00:00"))
        return parsed <= datetime.now(timezone.utc)

    def _filter_active_modules(self, granted_modules: list[str], module_expirations: dict[str, str], *, role: str = "") -> list[str]:
        if self._normalize_role(role) == "admin":
            return list(DEFAULT_ALLOWED_MODULES)
        active: list[str] = []
        for module in granted_modules:
            if module not in DEFAULT_ALLOWED_MODULES:
                continue
            if self._is_module_expired(module_expirations.get(module, "")):
                continue
            if module not in active:
                active.append(module)
        return active
