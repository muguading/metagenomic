from __future__ import annotations

import sqlite3
import sys
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
                    group_name TEXT NOT NULL DEFAULT '',
                    display_name TEXT,
                    created_at TEXT NOT NULL,
                    updated_at TEXT NOT NULL
                )
                """
            )
            columns = {row["name"] for row in conn.execute("PRAGMA table_info(users)").fetchall()}
            if "group_name" not in columns:
                conn.execute("ALTER TABLE users ADD COLUMN group_name TEXT NOT NULL DEFAULT ''")
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
                ("max_concurrent_tasks", "2", utc_now_iso()),
            )

    def connect(self) -> sqlite3.Connection:
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row
        return conn

    def authenticate(self, username: str, password: str) -> dict[str, Any] | None:
        with self.connect() as conn:
            row = conn.execute(
                "SELECT username, password_hash, role, group_name, display_name, created_at, updated_at FROM users WHERE username = ?",
                (username,),
            ).fetchone()
        if row is None or not check_password_hash(row["password_hash"], password):
            return None
        return self._serialize_user_row(row)

    def get_user(self, username: str) -> dict[str, Any]:
        with self.connect() as conn:
            row = conn.execute(
                "SELECT username, role, group_name, display_name, created_at, updated_at FROM users WHERE username = ?",
                (username,),
            ).fetchone()
        if row is None:
            raise KeyError(f"用户不存在: {username}")
        return self._serialize_user_row(row)

    def list_users(self) -> list[dict[str, Any]]:
        with self.connect() as conn:
            rows = conn.execute(
                "SELECT username, role, group_name, display_name, created_at, updated_at FROM users ORDER BY role DESC, group_name ASC, username ASC"
            ).fetchall()
        return [self._serialize_user_row(row) for row in rows]

    def create_user(
        self,
        username: str,
        password: str,
        role: str,
        display_name: str | None = None,
        group_name: str | None = None,
    ) -> dict[str, Any]:
        role = self._normalize_role(role)
        if role not in {"admin", "group_admin", "user"}:
            raise ValueError("role must be admin, group_admin or user")
        now = utc_now_iso()
        normalized_group = "" if role == "admin" else (group_name or "").strip()
        if role != "admin" and not normalized_group:
            raise ValueError("group_name is required for non-admin users")
        try:
            with self.connect() as conn:
                conn.execute(
                    """
                    INSERT INTO users (username, password_hash, role, group_name, display_name, created_at, updated_at)
                    VALUES (?, ?, ?, ?, ?, ?, ?)
                    """,
                    (username, generate_password_hash(password), role, normalized_group, display_name or "", now, now),
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
        next_display_name = display_name if display_name is not None else current["display_name"]
        next_hash = generate_password_hash(new_password) if new_password else current["password_hash"]
        try:
            with self.connect() as conn:
                conn.execute(
                    """
                    UPDATE users
                    SET username = ?, role = ?, group_name = ?, display_name = ?, password_hash = ?, updated_at = ?
                    WHERE username = ?
                    """,
                    (next_username, next_role, next_group_name, next_display_name, next_hash, utc_now_iso(), username),
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

    def _serialize_user_row(self, row: sqlite3.Row) -> dict[str, Any]:
        return {
            "username": row["username"],
            "role": self._normalize_role(row["role"]),
            "group_name": row["group_name"] or "",
            "display_name": row["display_name"] or "",
            "created_at": row["created_at"],
            "updated_at": row["updated_at"],
        }

    def _normalize_role(self, role: str) -> str:
        text = str(role or "").strip()
        return "group_admin" if text == "group" else text
