from __future__ import annotations

import csv
import io
import zipfile
from xml.sax.saxutils import escape as xml_escape

from .import_templates import _extract_batch_upload_headers, _parse_database_batch_upload
from .task_manager import ValidationError

def read_sample_meta_upload(filename: str, content: bytes) -> tuple[list[str], list[dict[str, str]]]:
    """Read the headers and rows used to map exported analysis sample names."""
    return _extract_batch_upload_headers(filename, content), _parse_database_batch_upload(filename, content)


def _normalize_export_columns(columns: object) -> list[str]:
    if not isinstance(columns, list):
        raise ValidationError("导出表头格式不正确")
    return [str(value or "") for value in columns]


def _normalize_export_rows(rows: object) -> list[list[str]]:
    if not isinstance(rows, list):
        raise ValidationError("导出表格内容格式不正确")
    normalized: list[list[str]] = []
    for row in rows:
        if not isinstance(row, list):
            raise ValidationError("导出表格内容格式不正确")
        normalized.append([str(value or "") for value in row])
    return normalized


def _normalize_export_sheets(sheets: object) -> list[dict[str, object]]:
    if not sheets:
        return []
    if not isinstance(sheets, list):
        raise ValidationError("导出 sheet 格式不正确")
    normalized: list[dict[str, object]] = []
    for index, sheet in enumerate(sheets, start=1):
        if not isinstance(sheet, dict):
            raise ValidationError("导出 sheet 格式不正确")
        title = str(sheet.get("title") or f"Sheet{index}").strip() or f"Sheet{index}"
        normalized.append({
            "title": title,
            "columns": _normalize_export_columns(sheet.get("columns", [])),
            "rows": _normalize_export_rows(sheet.get("rows", [])),
        })
    return normalized


def _sanitize_export_filename(name: str) -> str:
    safe = "".join(char if char.isalnum() or char in {"-", "_"} else "_" for char in name.strip())
    return safe.strip("_") or "export"


def _build_delimited_bytes(columns: list[str], rows: list[list[str]], delimiter: str) -> bytes:
    buffer = io.StringIO()
    writer = csv.writer(buffer, delimiter=delimiter, lineterminator="\n")
    if columns:
        writer.writerow(columns)
    writer.writerows(rows)
    return buffer.getvalue().encode("utf-8-sig")


def _xlsx_column_name(index: int) -> str:
    result = ""
    current = index + 1
    while current > 0:
        current, remainder = divmod(current - 1, 26)
        result = chr(65 + remainder) + result
    return result


def _sanitize_xlsx_sheet_name(name: str, used: set[str]) -> str:
    invalid = set("[]:*?/\\")
    base = "".join("_" if char in invalid else char for char in str(name or "").strip()).strip("'") or "Sheet"
    base = base[:31] or "Sheet"
    candidate = base
    counter = 2
    while candidate.lower() in used:
        suffix = f"_{counter}"
        candidate = f"{base[:31 - len(suffix)]}{suffix}"
        counter += 1
    used.add(candidate.lower())
    return candidate


def _build_xlsx_sheet_xml(columns: list[str], rows: list[list[str]]) -> str:
    def build_row_xml(row_number: int, values: list[str]) -> str:
        cells = []
        for index, value in enumerate(values):
            reference = f"{_xlsx_column_name(index)}{row_number}"
            text = xml_escape(str(value or ""))
            cells.append(
                f'<c r="{reference}" t="inlineStr"><is><t xml:space="preserve">{text}</t></is></c>'
            )
        return f'<row r="{row_number}">{"".join(cells)}</row>'

    sheet_rows = []
    current_row = 1
    if columns:
        sheet_rows.append(build_row_xml(current_row, columns))
        current_row += 1
    for row in rows:
        sheet_rows.append(build_row_xml(current_row, row))
        current_row += 1
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        f'<sheetData>{"".join(sheet_rows)}</sheetData>'
        '</worksheet>'
    )


def _build_xlsx_workbook_bytes(sheets: list[dict[str, object]]) -> bytes:
    normalized_sheets = sheets or [{"title": "Sheet1", "columns": [], "rows": []}]
    used_names: set[str] = set()
    sheet_meta = [
        {
            "name": _sanitize_xlsx_sheet_name(str(sheet.get("title") or f"Sheet{index}"), used_names),
            "columns": sheet.get("columns") if isinstance(sheet.get("columns"), list) else [],
            "rows": sheet.get("rows") if isinstance(sheet.get("rows"), list) else [],
            "index": index,
        }
        for index, sheet in enumerate(normalized_sheets, start=1)
    ]
    workbook_sheets_xml = "".join(
        f'<sheet name="{xml_escape(item["name"])}" sheetId="{item["index"]}" r:id="rId{item["index"]}"/>'
        for item in sheet_meta
    )
    workbook_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        f'<sheets>{workbook_sheets_xml}</sheets>'
        '</workbook>'
    )
    worksheet_overrides = "".join(
        f'<Override PartName="/xl/worksheets/sheet{item["index"]}.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
        for item in sheet_meta
    )
    content_types_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
        '<Default Extension="xml" ContentType="application/xml"/>'
        '<Override PartName="/xl/workbook.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
        f'{worksheet_overrides}'
        '<Override PartName="/xl/styles.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>'
        '</Types>'
    )
    rels_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" '
        'Target="xl/workbook.xml"/>'
        '</Relationships>'
    )
    workbook_relationships = "".join(
        f'<Relationship Id="rId{item["index"]}" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" '
        f'Target="worksheets/sheet{item["index"]}.xml"/>'
        for item in sheet_meta
    )
    workbook_rels_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        f'{workbook_relationships}'
        f'<Relationship Id="rId{len(sheet_meta) + 1}" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" '
        'Target="styles.xml"/>'
        '</Relationships>'
    )
    styles_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        '<fonts count="1"><font><sz val="11"/><name val="Calibri"/></font></fonts>'
        '<fills count="1"><fill><patternFill patternType="none"/></fill></fills>'
        '<borders count="1"><border/></borders>'
        '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>'
        '<cellXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/></cellXfs>'
        '<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>'
        '</styleSheet>'
    )

    output = io.BytesIO()
    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("[Content_Types].xml", content_types_xml)
        archive.writestr("_rels/.rels", rels_xml)
        archive.writestr("xl/workbook.xml", workbook_xml)
        archive.writestr("xl/_rels/workbook.xml.rels", workbook_rels_xml)
        archive.writestr("xl/styles.xml", styles_xml)
        for item in sheet_meta:
            archive.writestr(
                f'xl/worksheets/sheet{item["index"]}.xml',
                _build_xlsx_sheet_xml(item["columns"], item["rows"]),
            )
    return output.getvalue()


def _build_xlsx_bytes(title: str, columns: list[str], rows: list[list[str]]) -> bytes:
    sheet_name = _sanitize_export_filename(title)[:31] or "Sheet1"

    def build_row_xml(row_number: int, values: list[str]) -> str:
        cells = []
        for index, value in enumerate(values):
            reference = f"{_xlsx_column_name(index)}{row_number}"
            text = xml_escape(str(value or ""))
            cells.append(
                f'<c r="{reference}" t="inlineStr"><is><t xml:space="preserve">{text}</t></is></c>'
            )
        return f'<row r="{row_number}">{"".join(cells)}</row>'

    sheet_rows = []
    current_row = 1
    if columns:
        sheet_rows.append(build_row_xml(current_row, columns))
        current_row += 1
    for row in rows:
        sheet_rows.append(build_row_xml(current_row, row))
        current_row += 1

    workbook_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        f'<sheets><sheet name="{xml_escape(sheet_name)}" sheetId="1" r:id="rId1"/></sheets>'
        '</workbook>'
    )
    sheet_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        f'<sheetData>{"".join(sheet_rows)}</sheetData>'
        '</worksheet>'
    )
    content_types_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
        '<Default Extension="xml" ContentType="application/xml"/>'
        '<Override PartName="/xl/workbook.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
        '<Override PartName="/xl/worksheets/sheet1.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
        '<Override PartName="/xl/styles.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>'
        '</Types>'
    )
    rels_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" '
        'Target="xl/workbook.xml"/>'
        '</Relationships>'
    )
    workbook_rels_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" '
        'Target="worksheets/sheet1.xml"/>'
        '<Relationship Id="rId2" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" '
        'Target="styles.xml"/>'
        '</Relationships>'
    )
    styles_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        '<fonts count="1"><font><sz val="11"/><name val="Calibri"/></font></fonts>'
        '<fills count="1"><fill><patternFill patternType="none"/></fill></fills>'
        '<borders count="1"><border/></borders>'
        '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>'
        '<cellXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/></cellXfs>'
        '<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>'
        '</styleSheet>'
    )

    output = io.BytesIO()
    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("[Content_Types].xml", content_types_xml)
        archive.writestr("_rels/.rels", rels_xml)
        archive.writestr("xl/workbook.xml", workbook_xml)
        archive.writestr("xl/_rels/workbook.xml.rels", workbook_rels_xml)
        archive.writestr("xl/worksheets/sheet1.xml", sheet_xml)
        archive.writestr("xl/styles.xml", styles_xml)
    return output.getvalue()
