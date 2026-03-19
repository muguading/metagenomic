from __future__ import annotations

import csv
import io
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET


EXPECTED_COLUMNS = [
    "genome_id",
    "sample_name",
    "species_name",
    "taxid",
    "genome_file_path",
    "gender",
    "country",
    "location_province",
    "location_city",
    "location_district",
    "location_detail",
    "collection_time",
    "sample_type",
    "sequencing_method",
    "description",
]


class BulkImportFormatError(ValueError):
    pass


def parse_bulk_import_rows(filename: str, file_bytes: bytes) -> list[dict[str, str]]:
    suffix = Path(filename or "").suffix.lower()
    if suffix == ".csv":
        rows = _parse_delimited_rows(file_bytes, delimiter=",")
    elif suffix == ".tsv":
        rows = _parse_delimited_rows(file_bytes, delimiter="\t")
    elif suffix == ".xlsx":
        rows = _parse_xlsx_rows(file_bytes)
    else:
        raise BulkImportFormatError("Only .xlsx, .csv, and .tsv files are supported")
    _validate_headers(rows)
    return rows


def _parse_delimited_rows(file_bytes: bytes, *, delimiter: str) -> list[dict[str, str]]:
    text = file_bytes.decode("utf-8-sig")
    reader = csv.DictReader(io.StringIO(text), delimiter=delimiter)
    return [_normalize_row(row) for row in reader]


def _parse_xlsx_rows(file_bytes: bytes) -> list[dict[str, str]]:
    with zipfile.ZipFile(io.BytesIO(file_bytes)) as archive:
        shared_strings = _read_shared_strings(archive)
        sheet_xml = archive.read("xl/worksheets/sheet1.xml")
        root = ET.fromstring(sheet_xml)
        namespace = {"main": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
        rows: list[list[str]] = []
        for row in root.findall(".//main:sheetData/main:row", namespace):
            cells: list[str] = []
            current_column = 0
            for cell in row.findall("main:c", namespace):
                reference = cell.attrib.get("r", "")
                column_index = _column_index_from_reference(reference)
                while current_column < column_index:
                    cells.append("")
                    current_column += 1
                value = _read_xlsx_cell(cell, shared_strings, namespace)
                cells.append(value)
                current_column += 1
            rows.append(cells)
    if not rows:
        return []
    header = rows[0]
    data_rows = rows[1:]
    normalized: list[dict[str, str]] = []
    for values in data_rows:
        row = {header[index]: values[index] if index < len(values) else "" for index in range(len(header))}
        normalized.append(_normalize_row(row))
    return normalized


def _read_shared_strings(archive: zipfile.ZipFile) -> list[str]:
    try:
        raw = archive.read("xl/sharedStrings.xml")
    except KeyError:
        return []
    root = ET.fromstring(raw)
    namespace = {"main": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    values: list[str] = []
    for item in root.findall("main:si", namespace):
        parts = [node.text or "" for node in item.findall(".//main:t", namespace)]
        values.append("".join(parts))
    return values


def _read_xlsx_cell(cell: ET.Element, shared_strings: list[str], namespace: dict[str, str]) -> str:
    cell_type = cell.attrib.get("t", "")
    value_node = cell.find("main:v", namespace)
    if value_node is None or value_node.text is None:
        inline = cell.find("main:is/main:t", namespace)
        return (inline.text or "").strip() if inline is not None else ""
    raw = value_node.text
    if cell_type == "s":
        index = int(raw)
        return shared_strings[index].strip() if 0 <= index < len(shared_strings) else ""
    return raw.strip()


def _column_index_from_reference(reference: str) -> int:
    letters = "".join(char for char in reference if char.isalpha()).upper()
    index = 0
    for char in letters:
        index = index * 26 + (ord(char) - 64)
    return max(index - 1, 0)


def _normalize_row(row: dict[str, object]) -> dict[str, str]:
    return {str(key or "").strip(): str(value or "").strip() for key, value in row.items() if key is not None}


def _validate_headers(rows: list[dict[str, str]]) -> None:
    if not rows:
        raise BulkImportFormatError("Import file is empty")
    headers = set(rows[0].keys())
    missing = [column for column in EXPECTED_COLUMNS if column not in headers]
    if missing:
        raise BulkImportFormatError(f"Import file is missing required columns: {', '.join(missing)}")
