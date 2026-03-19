from __future__ import annotations

import sqlite3
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from werkzeug.security import check_password_hash, generate_password_hash


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
                    display_name TEXT,
                    created_at TEXT NOT NULL,
                    updated_at TEXT NOT NULL
                )
                """
            )
            conn.execute(
                """
                CREATE TABLE IF NOT EXISTS settings (
                    key TEXT PRIMARY KEY,
                    value TEXT NOT NULL,
                    updated_at TEXT NOT NULL
                )
                """
            )
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

    def connect(self) -> sqlite3.Connection:
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row
        return conn

    def authenticate(self, username: str, password: str) -> dict[str, Any] | None:
        with self.connect() as conn:
            row = conn.execute(
                "SELECT username, password_hash, role, display_name, created_at, updated_at FROM users WHERE username = ?",
                (username,),
            ).fetchone()
        if row is None or not check_password_hash(row["password_hash"], password):
            return None
        return self._serialize_user_row(row)

    def get_user(self, username: str) -> dict[str, Any]:
        with self.connect() as conn:
            row = conn.execute(
                "SELECT username, role, display_name, created_at, updated_at FROM users WHERE username = ?",
                (username,),
            ).fetchone()
        if row is None:
            raise KeyError(f"用户不存在: {username}")
        return self._serialize_user_row(row)

    def list_users(self) -> list[dict[str, Any]]:
        with self.connect() as conn:
            rows = conn.execute(
                "SELECT username, role, display_name, created_at, updated_at FROM users ORDER BY role DESC, username ASC"
            ).fetchall()
        return [self._serialize_user_row(row) for row in rows]

    def create_user(self, username: str, password: str, role: str, display_name: str | None = None) -> dict[str, Any]:
        if role not in {"admin", "user"}:
            raise ValueError("role must be admin or user")
        now = utc_now_iso()
        with self.connect() as conn:
            conn.execute(
                """
                INSERT INTO users (username, password_hash, role, display_name, created_at, updated_at)
                VALUES (?, ?, ?, ?, ?, ?)
                """,
                (username, generate_password_hash(password), role, display_name or "", now, now),
            )
        return self.get_user(username)

    def update_user(
        self,
        username: str,
        *,
        role: str | None = None,
        display_name: str | None = None,
        new_password: str | None = None,
    ) -> dict[str, Any]:
        current = self.get_user_with_hash(username)
        next_role = role or current["role"]
        next_display_name = display_name if display_name is not None else current["display_name"]
        next_hash = generate_password_hash(new_password) if new_password else current["password_hash"]
        with self.connect() as conn:
            conn.execute(
                """
                UPDATE users
                SET role = ?, display_name = ?, password_hash = ?, updated_at = ?
                WHERE username = ?
                """,
                (next_role, next_display_name, next_hash, utc_now_iso(), username),
            )
        return self.get_user(username)

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

    def _serialize_user_row(self, row: sqlite3.Row) -> dict[str, Any]:
        return {
            "username": row["username"],
            "role": row["role"],
            "display_name": row["display_name"] or "",
            "created_at": row["created_at"],
            "updated_at": row["updated_at"],
        }
