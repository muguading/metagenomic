from __future__ import annotations

from datetime import datetime
from pathlib import Path
from uuid import uuid4


class ValidationError(ValueError):
    pass


def validate_required_text(value: str, field_name: str) -> str:
    normalized = value.strip()
    if not normalized:
        raise ValidationError(f"{field_name} cannot be empty")
    return normalized


def validate_taxid(taxid: int) -> int:
    if taxid <= 0:
        raise ValidationError("taxid must be a positive integer")
    return taxid


def validate_optional_text(value: object, field_name: str) -> str | None:
    if value in (None, ""):
        return None
    text = str(value).strip()
    if not text:
        return None
    return text


def validate_optional_datetime(value: object, field_name: str) -> datetime | None:
    if value in (None, ""):
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        return datetime.fromisoformat(text)
    except ValueError as exc:
        raise ValidationError(f"{field_name} must be a valid datetime") from exc


def validate_optional_location(value: object, field_name: str = "location") -> dict[str, str] | None:
    if value in (None, "", {}):
        return None
    if not isinstance(value, dict):
        raise ValidationError(f"{field_name} must be an object")
    province = str(value.get("province", "")).strip()
    city = str(value.get("city", "")).strip()
    district = str(value.get("district", "")).strip()
    detail = str(value.get("detail", "")).strip()
    if not any([province, city, district, detail]):
        return None
    if any([province, city, district, detail]) and not all([province, city, district]):
        raise ValidationError(f"{field_name} must include province, city, and district")
    return {
        "province": province,
        "city": city,
        "district": district,
        "detail": detail,
    }


def validate_fasta_file(file_path: str) -> tuple[str, int]:
    path = Path(file_path).expanduser().resolve()
    if not path.is_file():
        raise ValidationError(f"FASTA file does not exist: {path}")

    total_length = 0
    header_seen = False
    sequence_seen = False
    allowed = set("ACGTRYSWKMBDHVNacgtryswkmbdhvn-.*")

    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if len(line) == 1:
                    raise ValidationError(f"Invalid FASTA header at line {line_number}: missing sequence identifier")
                header_seen = True
                continue
            if not header_seen:
                raise ValidationError(f"Invalid FASTA format at line {line_number}: sequence data before header")
            invalid = set(line) - allowed
            if invalid:
                invalid_chars = "".join(sorted(invalid))
                raise ValidationError(f"Invalid FASTA sequence at line {line_number}: unsupported characters {invalid_chars!r}")
            total_length += len(line.replace("-", "").replace(".", "").replace("*", ""))
            sequence_seen = True

    if not header_seen or not sequence_seen:
        raise ValidationError("FASTA file is empty or missing sequence records")

    return str(path), total_length


def validate_custom_metadata(value: object) -> list[dict[str, object]]:
    if value in (None, "", []):
        return []
    if not isinstance(value, list):
        raise ValidationError("custom_metadata must be a list")

    normalized: list[dict[str, object]] = []
    seen_keys: set[str] = set()
    seen_labels: set[str] = set()
    for index, item in enumerate(value, start=1):
        if not isinstance(item, dict):
            raise ValidationError(f"custom_metadata item {index} must be an object")

        label = validate_required_text(str(item.get("label", "")), f"custom_metadata[{index}].label")
        field_type = validate_required_text(str(item.get("type", "")), f"custom_metadata[{index}].type").lower()
        if field_type not in {"text", "select", "file", "datetime", "location", "country"}:
            raise ValidationError(f"custom_metadata[{index}].type must be text, select, file, datetime, location, or country")
        normalized_label = " ".join(label.split()).casefold()
        if normalized_label in seen_labels:
            raise ValidationError(f"custom_metadata[{index}].label is duplicated: {label!r}")
        seen_labels.add(normalized_label)

        field_key = str(item.get("key") or f"meta_{uuid4().hex[:8]}")
        if field_key in seen_keys:
            raise ValidationError(f"custom_metadata[{index}].key is duplicated")
        seen_keys.add(field_key)
        normalized_item: dict[str, object] = {
            "key": field_key,
            "label": label,
            "type": field_type,
        }

        if field_type == "text":
            normalized_item["value"] = str(item.get("value", "")).strip()
        elif field_type == "country":
            normalized_item["value"] = str(item.get("value", "")).strip()
        elif field_type == "datetime":
            raw_value = str(item.get("value", "")).strip()
            if raw_value:
                try:
                    datetime.fromisoformat(raw_value)
                except ValueError as exc:
                    raise ValidationError(f"custom_metadata[{index}].value must be a valid datetime") from exc
            normalized_item["value"] = raw_value
        elif field_type == "select":
            raw_options = item.get("options", [])
            if isinstance(raw_options, str):
                raw_options = [part.strip() for part in raw_options.split(",")]
            if not isinstance(raw_options, list):
                raise ValidationError(f"custom_metadata[{index}].options must be a list")
            options = [str(option).strip() for option in raw_options if str(option).strip()]
            if not options:
                raise ValidationError(f"custom_metadata[{index}].options cannot be empty")
            selected = str(item.get("value", "")).strip()
            if selected and selected not in options:
                raise ValidationError(f"custom_metadata[{index}].value must be one of the configured options")
            normalized_item["options"] = options
            normalized_item["value"] = selected
        elif field_type == "file":
            path_value = str(item.get("value", "")).strip()
            if path_value:
                path = Path(path_value).expanduser().resolve()
                if not path.is_file():
                    raise ValidationError(f"custom_metadata[{index}].file does not exist: {path}")
                normalized_item["value"] = str(path)
                normalized_item["filename"] = str(item.get("filename", "")).strip() or path.name
            else:
                normalized_item["value"] = ""
                if item.get("filename"):
                    normalized_item["filename"] = str(item.get("filename")).strip()
        else:
            raw_location = item.get("value", {})
            if raw_location in (None, ""):
                raw_location = {}
            if not isinstance(raw_location, dict):
                raise ValidationError(f"custom_metadata[{index}].value must be an object for location type")
            province = str(raw_location.get("province", "")).strip()
            city = str(raw_location.get("city", "")).strip()
            district = str(raw_location.get("district", "")).strip()
            detail = str(raw_location.get("detail", "")).strip()
            if any([province, city, district, detail]) and not all([province, city, district]):
                raise ValidationError(f"custom_metadata[{index}].location must include province, city, and district")
            normalized_item["value"] = {
                "province": province,
                "city": city,
                "district": district,
                "detail": detail,
            }

        normalized.append(normalized_item)

    return normalized
