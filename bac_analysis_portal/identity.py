from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class UserIdentity:
    username: str
    role: str
    group_name: str

    @classmethod
    def from_mapping(cls, value: dict[str, str]) -> "UserIdentity":
        return cls(
            username=str(value.get("username") or "").strip(),
            role=str(value.get("role") or "").strip(),
            group_name=str(value.get("group_name") or "").strip(),
        )
