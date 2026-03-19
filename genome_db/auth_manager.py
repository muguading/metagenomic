from __future__ import annotations

import os
from dataclasses import dataclass
from datetime import datetime

from sqlalchemy import select, update
from sqlalchemy.orm import Session, sessionmaker
from werkzeug.security import check_password_hash, generate_password_hash

from .database import session_scope
from .models import AuditLog, Genome, User
from .validators import ValidationError, validate_required_text


DEFAULT_ADMIN_USERNAME = os.environ.get("GENOME_DB_ADMIN_USERNAME", "admin")
DEFAULT_ADMIN_PASSWORD = os.environ.get("GENOME_DB_ADMIN_PASSWORD", "admin123")


@dataclass(slots=True)
class AuthenticatedUser:
    username: str
    role: str


class AuthenticationError(ValueError):
    pass


class AuthorizationError(PermissionError):
    pass


class UserManager:
    def __init__(self, session_factory: sessionmaker[Session]) -> None:
        self.session_factory = session_factory
        self.ensure_default_admin()

    def ensure_default_admin(self) -> None:
        with session_scope(self.session_factory) as session:
            existing = session.scalar(select(User).where(User.username == DEFAULT_ADMIN_USERNAME))
            if existing is None:
                session.add(
                    User(
                        username=DEFAULT_ADMIN_USERNAME,
                        password_hash=generate_password_hash(DEFAULT_ADMIN_PASSWORD),
                        role="admin",
                    )
                )

    def authenticate(self, username: str, password: str) -> AuthenticatedUser:
        username = validate_required_text(username, "username")
        password = validate_required_text(password, "password")
        with session_scope(self.session_factory) as session:
            user = session.scalar(select(User).where(User.username == username))
            if user is None or not check_password_hash(user.password_hash, password):
                raise AuthenticationError("Invalid username or password")
            user.last_login_time = datetime.utcnow()
            session.flush()
            return AuthenticatedUser(username=user.username, role=user.role)

    def create_user(self, *, username: str, password: str, role: str) -> dict[str, object]:
        return self._create_user(username=username, password=password, role=role, display_name=None, email=None)

    def register_user(
        self,
        *,
        username: str,
        password: str,
        display_name: str | None = None,
        email: str | None = None,
    ) -> dict[str, object]:
        return self._create_user(
            username=username,
            password=password,
            role="user",
            display_name=display_name,
            email=email,
        )

    def _create_user(
        self,
        *,
        username: str,
        password: str,
        role: str,
        display_name: str | None,
        email: str | None,
    ) -> dict[str, object]:
        username = validate_required_text(username, "username")
        password = validate_required_text(password, "password")
        role = role.strip().lower()
        if role not in {"admin", "user"}:
            raise ValidationError("role must be 'admin' or 'user'")

        with session_scope(self.session_factory) as session:
            existing = session.scalar(select(User).where(User.username == username))
            if existing is not None:
                raise ValidationError(f"User already exists: {username}")
            user = User(
                username=username,
                password_hash=generate_password_hash(password),
                role=role,
                display_name=display_name.strip() if isinstance(display_name, str) and display_name.strip() else None,
                email=email.strip() if isinstance(email, str) and email.strip() else None,
            )
            session.add(user)
            session.flush()
            return user.to_dict()

    def list_users(self) -> list[dict[str, object]]:
        with session_scope(self.session_factory) as session:
            users = session.scalars(select(User).order_by(User.role.desc(), User.username.asc())).all()
            return [user.to_dict() for user in users]

    def get_user_profile(self, username: str) -> dict[str, object]:
        with session_scope(self.session_factory) as session:
            user = self._get_user(session, username)
            return user.to_dict()

    def update_user_profile(
        self,
        *,
        username: str,
        new_username: str,
        display_name: str | None,
        email: str | None,
    ) -> dict[str, object]:
        new_username = validate_required_text(new_username, "username")
        with session_scope(self.session_factory) as session:
            user = self._get_user(session, username)
            if new_username != username:
                existing = session.scalar(select(User).where(User.username == new_username))
                if existing is not None:
                    raise ValidationError(f"User already exists: {new_username}")
                user.username = new_username
                session.execute(
                    update(Genome)
                    .where(Genome.submitter == username)
                    .values(submitter=new_username)
                )
                session.execute(
                    update(AuditLog)
                    .where(AuditLog.operator == username)
                    .values(operator=new_username)
                )
            user.display_name = display_name.strip() if isinstance(display_name, str) and display_name.strip() else None
            user.email = email.strip() if isinstance(email, str) and email.strip() else None
            user.last_modified_time = datetime.utcnow()
            session.flush()
            return user.to_dict()

    def change_password(self, *, username: str, current_password: str, new_password: str) -> None:
        current_password = validate_required_text(current_password, "current_password")
        new_password = validate_required_text(new_password, "new_password")
        with session_scope(self.session_factory) as session:
            user = self._get_user(session, username)
            if not check_password_hash(user.password_hash, current_password):
                raise AuthenticationError("Current password is incorrect")
            user.password_hash = generate_password_hash(new_password)
            user.last_modified_time = datetime.utcnow()
            session.flush()

    def admin_reset_password(self, *, target_username: str, new_password: str) -> dict[str, object]:
        new_password = validate_required_text(new_password, "new_password")
        with session_scope(self.session_factory) as session:
            user = self._get_user(session, target_username)
            if user.role != "user":
                raise ValidationError("Administrator can only reset passwords for ordinary users")
            user.password_hash = generate_password_hash(new_password)
            user.last_modified_time = datetime.utcnow()
            session.flush()
            return user.to_dict()

    def _get_user(self, session: Session, username: str) -> User:
        user = session.scalar(select(User).where(User.username == username))
        if user is None:
            raise ValidationError(f"User not found: {username}")
        return user
