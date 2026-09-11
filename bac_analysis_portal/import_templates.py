from __future__ import annotations

import csv
import hashlib
import io
import time
import uuid
import zipfile
from datetime import datetime
from pathlib import Path
from typing import Any
from xml.etree import ElementTree as ET
from xml.sax.saxutils import escape as xml_escape

from .filesystem_helpers import _resolve_optional_existing_path
from .sample_library_manager import SampleLibraryManager
from .task_manager import ValidationError

BATCH_IMPORT_PRECHECK_TTL_SECONDS = 30 * 60

BATCH_INPUT_HEADER_MAP = {
    "样本名称": "sample_name",
    "三代数据": "third_gen",
    "三代测序": "third_gen",
    "二代数据左": "short_left",
    "二代测序左": "short_left",
    "二代数据右": "short_right",
    "二代测序右": "short_right",
    "物种信息": "species",
    "sample_name": "sample_name",
    "third_gen": "third_gen",
    "short_left": "short_left",
    "short_right": "short_right",
    "species": "species",
}

DATABASE_IMPORT_COLUMNS: list[tuple[str, str]] = [
    ("sample_name", "样本名称"),
    ("genome_id", "Genome ID"),
    ("taxid", "TaxID"),
    ("final_fasta_path", "Final FASTA 路径"),
    ("species_name", "物种预估"),
    ("mlst_species_name", "MLST 物种名称"),
    ("mlst_st", "MLST ST"),
    ("serotype_result", "血清型"),
    ("q20_rate", "Q20"),
    ("q30_rate", "Q30"),
    ("completeness", "完整性"),
    ("contamination", "污染率"),
    ("contig_count", "Contig 数"),
    ("plasmid_count", "质粒数"),
    ("resistance_count", "耐药条目数"),
    ("virulence_count", "毒力条目数"),
    ("resistance_gene_hits", "耐药基因摘要"),
    ("virulence_gene_hits", "毒力基因摘要"),
    ("resistance_mge_hits", "耐药基因-MGE"),
    ("virulence_mge_hits", "毒力基因-MGE"),
    ("sample_source", "样本来源"),
    ("collection_date", "采样日期"),
    ("gender", "性别"),
    ("country", "国家 / 地区"),
    ("host_info", "宿主 / 样本背景"),
    ("province", "省 / Province"),
    ("city", "市 / City"),
    ("district", "区县 / District"),
    ("location_detail", "位置详情"),
    ("sample_type", "输入类型"),
    ("sequencing_method", "测序 / 组装方法"),
    ("description", "描述"),
    ("note", "备注"),
    ("case_id", "病例/事件编号"),
    ("patient_id", "患者/个案编号"),
    ("surveillance_source", "监测来源"),
    ("suspected_syndrome", "疑似症候群"),
    ("submitting_unit", "送检单位"),
    ("collection_unit", "采样单位"),
    ("ward_department", "科室/病区"),
    ("specimen_category", "标本类别"),
    ("cluster_status", "聚集性状态"),
    ("epidemiology_link", "流行病学关联"),
    ("traditional_result", "传统检测结果"),
]

REFERENCE_IMPORT_COLUMNS = [
    ("host_name", "物种名称"),
    ("genome_name", "基因组名称"),
    ("taxid", "TaxID"),
    ("source_accession", "NCBI编号"),
    ("source_label", "来源标签"),
    ("fasta_path", "FASTA路径"),
    ("description", "备注"),
]


def _database_import_example_row() -> dict[str, str]:
    return {
        "sample_name": "Men-IGT",
        "genome_id": "Men-IGT",
        "taxid": "",
        "final_fasta_path": "/path/to/Men-IGT.final.fasta",
        "species_name": "Neisseria meningitidis",
        "mlst_species_name": "neisseria",
        "mlst_st": "ST4821",
        "serotype_result": "B",
        "q20_rate": "96.49%",
        "q30_rate": "94.07%",
        "completeness": "82.47%",
        "contamination": "10.42%",
        "contig_count": "613",
        "plasmid_count": "71",
        "resistance_count": "8",
        "virulence_count": "82",
        "resistance_gene_hits": "OXA-23、blaADC-25、aph(3')-Ia",
        "virulence_gene_hits": "ompA、bap、csuE",
        "resistance_mge_hits": "OXA-23（转座子）",
        "virulence_mge_hits": "bap（质粒）",
        "sample_source": "血液",
        "collection_date": "2026-03-23",
        "gender": "",
        "country": "China",
        "host_info": "人源",
        "province": "浙江",
        "city": "杭州",
        "district": "西湖区",
        "location_detail": "某医院",
        "sample_type": "fastq",
        "sequencing_method": "spades",
        "description": "示例样本",
        "note": "可按需补充",
        "case_id": "HP-CDC-2026-001",
        "patient_id": "PID-20260323-01",
        "surveillance_source": "ICU",
        "suspected_syndrome": "脑膜炎/脑膜脑炎",
        "submitting_unit": "黄浦区疾控中心",
        "collection_unit": "某三甲医院神经内科",
        "ward_department": "ICU",
        "specimen_category": "脑脊液",
        "cluster_status": "待判定",
        "epidemiology_link": "与 2026 春季 ICU 监测事件关联",
        "traditional_result": "培养/PCR结果待补录",
    }

def _database_import_headers() -> list[str]:
    return [item[0] for item in DATABASE_IMPORT_COLUMNS]

def _reference_import_headers() -> list[str]:
    return [item[0] for item in REFERENCE_IMPORT_COLUMNS]

def _reference_import_example_row(category: str) -> dict[str, str]:
    if category == "host":
        return {
            "host_name": "Homo sapiens",
            "genome_name": "GRCh38",
            "taxid": "9606",
            "source_accession": "GCF_000001405.40",
            "source_label": "本地库 / 人源参考",
            "fasta_path": "/data/reference/human/GRCh38.fa",
            "description": "用于去宿主分析",
        }
    return {
        "host_name": "Mycobacterium tuberculosis",
        "genome_name": "H37Rv",
        "taxid": "1773",
        "source_accession": "GCF_000195955.2",
        "source_label": "本地库 / 参考菌株",
        "fasta_path": "/data/reference/pathogen/H37Rv.fa",
        "description": "病原参考基因组示例",
    }

def _reference_import_labels(category: str) -> list[str]:
    species_label = "宿主名称" if category == "host" else "病原名称"
    return [species_label if key == "host_name" else label for key, label in REFERENCE_IMPORT_COLUMNS]

def _build_reference_import_template_text(category: str, delimiter: str) -> str:
    output = io.StringIO()
    writer = csv.writer(output, delimiter=delimiter, lineterminator="\n")
    headers = _reference_import_headers()
    labels = _reference_import_labels(category)
    example = _reference_import_example_row(category)
    writer.writerow(headers)
    writer.writerow(labels)
    writer.writerow([example.get(header, "") for header in headers])
    return output.getvalue()

def _build_reference_import_template_xlsx(category: str) -> bytes:
    headers = _reference_import_headers()
    labels = _reference_import_labels(category)
    example = _reference_import_example_row(category)
    rows = [headers, labels, [example.get(header, "") for header in headers]]

    def build_sheet_xml() -> str:
        row_xml: list[str] = []
        for row_index, row in enumerate(rows, start=1):
            cell_xml: list[str] = []
            for col_index, value in enumerate(row, start=1):
                ref = f"{_xlsx_column_name(col_index)}{row_index}"
                text = xml_escape(str(value or ""))
                cell_xml.append(f'<c r="{ref}" t="inlineStr"><is><t>{text}</t></is></c>')
            row_xml.append(f'<row r="{row_index}">{"".join(cell_xml)}</row>')
        return (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
            "<sheetData>"
            + "".join(row_xml)
            + "</sheetData></worksheet>"
        )

    workbook_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        '<sheets><sheet name="references" sheetId="1" r:id="rId1"/></sheets></workbook>'
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
    root_rels_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" '
        'Target="xl/workbook.xml"/>'
        '<Relationship Id="rId2" '
        'Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" '
        'Target="docProps/core.xml"/>'
        '<Relationship Id="rId3" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties" '
        'Target="docProps/app.xml"/>'
        '</Relationships>'
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
        '<Override PartName="/docProps/core.xml" '
        'ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>'
        '<Override PartName="/docProps/app.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>'
        '</Types>'
    )
    styles_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        '<fonts count="1"><font><sz val="11"/><name val="Calibri"/></font></fonts>'
        '<fills count="1"><fill><patternFill patternType="none"/></fill></fills>'
        '<borders count="1"><border/></borders>'
        '<cellStyleXfs count="1"><xf/></cellStyleXfs>'
        '<cellXfs count="1"><xf xfId="0"/></cellXfs>'
        '<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>'
        '</styleSheet>'
    )
    now = datetime.utcnow().replace(microsecond=0).isoformat() + "Z"
    core_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" '
        'xmlns:dc="http://purl.org/dc/elements/1.1/" '
        'xmlns:dcterms="http://purl.org/dc/terms/" '
        'xmlns:dcmitype="http://purl.org/dc/dcmitype/" '
        'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">'
        '<dc:title>reference_database_import_template</dc:title>'
        '<dc:creator>bac_analysis_portal</dc:creator>'
        f'<dcterms:created xsi:type="dcterms:W3CDTF">{now}</dcterms:created>'
        f'<dcterms:modified xsi:type="dcterms:W3CDTF">{now}</dcterms:modified>'
        '</cp:coreProperties>'
    )
    app_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties" '
        'xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes">'
        '<Application>bac_analysis_portal</Application>'
        '</Properties>'
    )

    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("[Content_Types].xml", content_types_xml)
        archive.writestr("_rels/.rels", root_rels_xml)
        archive.writestr("xl/workbook.xml", workbook_xml)
        archive.writestr("xl/_rels/workbook.xml.rels", workbook_rels_xml)
        archive.writestr("xl/worksheets/sheet1.xml", build_sheet_xml())
        archive.writestr("xl/styles.xml", styles_xml)
        archive.writestr("docProps/core.xml", core_xml)
        archive.writestr("docProps/app.xml", app_xml)
    return buffer.getvalue()

def _build_database_import_template_text(delimiter: str) -> str:
    output = io.StringIO()
    writer = csv.writer(output, delimiter=delimiter, lineterminator="\n")
    headers = _database_import_headers()
    labels = [item[1] for item in DATABASE_IMPORT_COLUMNS]
    example = _database_import_example_row()
    writer.writerow(headers)
    writer.writerow(labels)
    writer.writerow([example.get(header, "") for header in headers])
    return output.getvalue()

def _xlsx_column_name(index: int) -> str:
    result = ""
    current = index
    while current > 0:
        current, remainder = divmod(current - 1, 26)
        result = chr(65 + remainder) + result
    return result

def _build_database_import_template_guide_rows(metadata_templates: list[dict[str, Any]] | None = None) -> list[list[str]]:
    template_map = {
        str(item.get("key") or "").strip(): item
        for item in (metadata_templates or [])
        if str(item.get("key") or "").strip()
    }
    rows: list[list[str]] = [
        ["模块", "字段 key", "显示名称", "填写要求", "标准词表/推荐项", "填写说明", "示例"],
        ["模板说明", "-", "-", "-", "-", "XLSX 模板包含两页：samples 为可直接导入的数据页；说明与词表示例 为填写说明页。", "-"],
        ["模板说明", "-", "-", "-", "-", "CSV/TSV 模板保持纯数据结构，适合程序批量生成或脚本处理。", "-"],
        ["模板说明", "-", "-", "-", "-", "导入时至少保证 sample_name、final_fasta_path 以及主档必填字段有效。", "-"],
        ["模板说明", "-", "-", "-", "-", "采样地点请分别填写 province/city/district/location_detail；location_detail 优先选标准点位。", "-"],
    ]
    example_row = _database_import_example_row()
    for key, label in DATABASE_IMPORT_COLUMNS:
        template = template_map.get(key, {})
        requirement = str(template.get("requirement") or "").strip()
        requirement_text = "必填" if requirement == "required" else "推荐" if requirement == "recommended" else "选填"
        options = template.get("options") if isinstance(template.get("options"), list) else []
        option_text = "；".join([str(item).strip() for item in options if str(item).strip()][:8])
        dictionary_name = str(template.get("dictionary_name") or "").strip()
        help_text = str(template.get("help_text") or "").strip()
        if dictionary_name:
            option_text = f"{dictionary_name}{'：' if option_text else ''}{option_text}"
        rows.append([
            "字段说明",
            key,
            label,
            requirement_text,
            option_text or "-",
            help_text or "-",
            str(example_row.get(key, "") or "-"),
        ])
    return rows

def _build_database_import_template_xlsx(metadata_templates: list[dict[str, Any]] | None = None) -> bytes:
    headers = _database_import_headers()
    labels = [item[1] for item in DATABASE_IMPORT_COLUMNS]
    example = _database_import_example_row()
    rows = [headers, labels, [example.get(header, "") for header in headers]]
    guide_rows = _build_database_import_template_guide_rows(metadata_templates)

    def build_sheet_xml(sheet_rows: list[list[str]]) -> str:
        row_xml: list[str] = []
        for row_index, row in enumerate(sheet_rows, start=1):
            cell_xml: list[str] = []
            for col_index, value in enumerate(row, start=1):
                ref = f"{_xlsx_column_name(col_index)}{row_index}"
                text = xml_escape(str(value or ""))
                cell_xml.append(
                    f'<c r="{ref}" t="inlineStr"><is><t>{text}</t></is></c>'
                )
            row_xml.append(f'<row r="{row_index}">{"".join(cell_xml)}</row>')
        return (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
            "<sheetData>"
            + "".join(row_xml)
            + "</sheetData></worksheet>"
        )

    workbook_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        '<sheets>'
        '<sheet name="samples" sheetId="1" r:id="rId1"/>'
        '<sheet name="说明与词表示例" sheetId="2" r:id="rId2"/>'
        '</sheets></workbook>'
    )
    workbook_rels_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" '
        'Target="worksheets/sheet1.xml"/>'
        '<Relationship Id="rId2" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" '
        'Target="worksheets/sheet2.xml"/>'
        '<Relationship Id="rId3" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" '
        'Target="styles.xml"/>'
        '</Relationships>'
    )
    root_rels_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" '
        'Target="xl/workbook.xml"/>'
        '<Relationship Id="rId2" '
        'Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" '
        'Target="docProps/core.xml"/>'
        '<Relationship Id="rId3" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties" '
        'Target="docProps/app.xml"/>'
        '</Relationships>'
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
        '<Override PartName="/xl/worksheets/sheet2.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
        '<Override PartName="/xl/styles.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>'
        '<Override PartName="/docProps/core.xml" '
        'ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>'
        '<Override PartName="/docProps/app.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>'
        '</Types>'
    )
    styles_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        '<fonts count="1"><font><sz val="11"/><name val="Calibri"/></font></fonts>'
        '<fills count="1"><fill><patternFill patternType="none"/></fill></fills>'
        '<borders count="1"><border/></borders>'
        '<cellStyleXfs count="1"><xf/></cellStyleXfs>'
        '<cellXfs count="1"><xf xfId="0"/></cellXfs>'
        '<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>'
        '</styleSheet>'
    )
    now = datetime.utcnow().replace(microsecond=0).isoformat() + "Z"
    core_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" '
        'xmlns:dc="http://purl.org/dc/elements/1.1/" '
        'xmlns:dcterms="http://purl.org/dc/terms/" '
        'xmlns:dcmitype="http://purl.org/dc/dcmitype/" '
        'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">'
        '<dc:title>sample_database_import_template</dc:title>'
        '<dc:creator>bac_analysis_portal</dc:creator>'
        f'<dcterms:created xsi:type="dcterms:W3CDTF">{now}</dcterms:created>'
        f'<dcterms:modified xsi:type="dcterms:W3CDTF">{now}</dcterms:modified>'
        '</cp:coreProperties>'
    )
    app_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties" '
        'xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes">'
        '<Application>bac_analysis_portal</Application>'
        '</Properties>'
    )

    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("[Content_Types].xml", content_types_xml)
        archive.writestr("_rels/.rels", root_rels_xml)
        archive.writestr("xl/workbook.xml", workbook_xml)
        archive.writestr("xl/_rels/workbook.xml.rels", workbook_rels_xml)
        archive.writestr("xl/worksheets/sheet1.xml", build_sheet_xml(rows))
        archive.writestr("xl/worksheets/sheet2.xml", build_sheet_xml(guide_rows))
        archive.writestr("xl/styles.xml", styles_xml)
        archive.writestr("docProps/core.xml", core_xml)
        archive.writestr("docProps/app.xml", app_xml)
    return buffer.getvalue()

def _build_database_import_metadata_items(
    row: dict[str, str],
    metadata_templates: list[dict[str, object]],
) -> list[dict[str, object]]:
    template_map = {
        str(item.get("key") or "").strip(): item
        for item in metadata_templates
        if isinstance(item, dict) and str(item.get("key") or "").strip()
    }
    items: list[dict[str, object]] = []
    for key in (
        "case_id",
        "patient_id",
        "surveillance_source",
        "suspected_syndrome",
        "submitting_unit",
        "collection_unit",
        "ward_department",
        "specimen_category",
        "cluster_status",
        "epidemiology_link",
        "traditional_result",
    ):
        value = str(row.get(key) or "").strip()
        if not value:
            continue
        template = template_map.get(key) or {}
        items.append(
            {
                "key": key,
                "label": str(template.get("label") or key).strip(),
                "type": str(template.get("type") or "text").strip() or "text",
                "options": template.get("options") or [],
                "value": value,
            }
        )

    province = str(row.get("province") or "").strip()
    city = str(row.get("city") or "").strip()
    district = str(row.get("district") or "").strip()
    detail = str(row.get("location_detail") or "").strip()
    if any([province, city, district, detail]):
        template = template_map.get("collection_site") or {}
        items.append(
            {
                "key": "collection_site",
                "label": str(template.get("label") or "采样地点").strip(),
                "type": "location",
                "options": [],
                "value": {
                    "province": province,
                    "city": city,
                    "district": district,
                    "detail": detail,
                },
            }
        )
    return items

def _parse_database_batch_upload(filename: str, content: bytes) -> list[dict[str, str]]:
    suffix = Path(str(filename or "")).suffix.lower()
    if suffix == ".csv":
        return _parse_database_text_table(content.decode("utf-8-sig", errors="ignore"), ",")
    if suffix == ".tsv":
        return _parse_database_text_table(content.decode("utf-8-sig", errors="ignore"), "\t")
    if suffix == ".xlsx":
        return _parse_database_xlsx_table(content)
    raise ValidationError("批量导入只支持 xlsx、csv 或 tsv")

def _extract_batch_upload_headers(filename: str, content: bytes) -> list[str]:
    suffix = Path(str(filename or "")).suffix.lower()
    if suffix in {".csv", ".tsv"}:
        delimiter = "," if suffix == ".csv" else "\t"
        reader = csv.reader(io.StringIO(content.decode("utf-8-sig", errors="ignore")), delimiter=delimiter)
        for row in reader:
            headers = [str(item or "").strip().replace("\ufeff", "") for item in row]
            if any(headers):
                return [item for item in headers if item]
        return []
    if suffix == ".xlsx":
        namespace = {"main": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
        with zipfile.ZipFile(io.BytesIO(content)) as archive:
            try:
                sheet_xml = archive.read("xl/worksheets/sheet1.xml")
            except KeyError as exc:
                raise ValidationError("xlsx 模板缺少第一个工作表") from exc
            shared_strings: list[str] = []
            try:
                shared_root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
                for item in shared_root.findall("main:si", namespace):
                    shared_strings.append("".join(node.text or "" for node in item.iterfind(".//main:t", namespace)))
            except KeyError:
                shared_strings = []
            root = ET.fromstring(sheet_xml)
            first_row = root.find(".//main:sheetData/main:row", namespace)
            if first_row is None:
                return []
            headers: list[str] = []
            for cell in first_row.findall("main:c", namespace):
                cell_type = cell.attrib.get("t", "")
                if cell_type == "inlineStr":
                    value = "".join(node.text or "" for node in cell.iterfind(".//main:t", namespace))
                else:
                    raw = cell.findtext("main:v", default="", namespaces=namespace)
                    if cell_type == "s" and raw:
                        try:
                            value = shared_strings[int(raw)]
                        except (ValueError, IndexError):
                            value = raw
                    else:
                        value = raw or ""
                headers.append(str(value or "").strip())
            return [item for item in headers if item]
    raise ValidationError("批量导入只支持 xlsx、csv 或 tsv")

def _parse_database_text_table(text: str, delimiter: str) -> list[dict[str, str]]:
    reader = csv.DictReader(io.StringIO(text), delimiter=delimiter)
    if not reader.fieldnames:
        raise ValidationError("导入文件缺少表头")
    return [
        {str(key or "").strip(): str(value or "").strip() for key, value in row.items()}
        for row in reader
    ]

def _parse_database_xlsx_table(content: bytes) -> list[dict[str, str]]:
    namespace = {"main": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    with zipfile.ZipFile(io.BytesIO(content)) as archive:
        try:
            sheet_xml = archive.read("xl/worksheets/sheet1.xml")
        except KeyError as exc:
            raise ValidationError("xlsx 模板缺少第一个工作表") from exc
        shared_strings: list[str] = []
        try:
            shared_root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            for item in shared_root.findall("main:si", namespace):
                text = "".join(node.text or "" for node in item.iterfind(".//main:t", namespace))
                shared_strings.append(text)
        except KeyError:
            shared_strings = []
        root = ET.fromstring(sheet_xml)
        rows: list[list[str]] = []
        for row in root.findall(".//main:sheetData/main:row", namespace):
            current: list[str] = []
            for cell in row.findall("main:c", namespace):
                cell_type = cell.attrib.get("t", "")
                value = ""
                if cell_type == "inlineStr":
                    value = "".join(node.text or "" for node in cell.iterfind(".//main:t", namespace))
                else:
                    raw = cell.findtext("main:v", default="", namespaces=namespace)
                    if cell_type == "s" and raw:
                        try:
                            value = shared_strings[int(raw)]
                        except (ValueError, IndexError):
                            value = raw
                    else:
                        value = raw or ""
                current.append(str(value).strip())
            rows.append(current)
    if not rows:
        raise ValidationError("xlsx 文件为空")
    headers = [str(item or "").strip() for item in rows[0]]
    if not any(headers):
        raise ValidationError("xlsx 文件缺少表头")
    parsed_rows: list[dict[str, str]] = []
    for raw_row in rows[1:]:
        padded = raw_row + [""] * max(0, len(headers) - len(raw_row))
        parsed_rows.append({headers[i]: str(padded[i] or "").strip() for i in range(len(headers)) if headers[i]})
    if parsed_rows and parsed_rows[0].get("sample_name") == "样本名称":
        parsed_rows = parsed_rows[1:]
    return parsed_rows

def _build_precheck_summary(
    *,
    filename: str,
    headers: list[str],
    rows: list[dict[str, str]],
    required_columns: list[str],
    valid_rows: int,
    empty_rows: int,
    issues: list[dict[str, str]],
) -> dict[str, object]:
    missing_columns = [column for column in required_columns if column not in headers]
    issue_limit = 20
    can_import = not missing_columns and valid_rows > 0 and not issues
    status = "ok" if can_import else ("error" if missing_columns or valid_rows == 0 else "warning")
    return {
        "status": status,
        "can_import": can_import,
        "filename": filename,
        "headers": headers,
        "required_columns": required_columns,
        "missing_columns": missing_columns,
        "total_rows": len(rows),
        "valid_rows": valid_rows,
        "empty_rows": empty_rows,
        "issue_count": len(issues),
        "issues": issues[:issue_limit],
        "truncated": len(issues) > issue_limit,
    }

def _prune_batch_import_precheck_cache(cache: dict[str, dict[str, object]]) -> None:
    now = time.time()
    expired_ids = [
        precheck_id
        for precheck_id, entry in cache.items()
        if now - float(entry.get("created_at") or 0) > BATCH_IMPORT_PRECHECK_TTL_SECONDS
    ]
    for precheck_id in expired_ids:
        cache.pop(precheck_id, None)

def _store_batch_import_precheck(
    cache: dict[str, dict[str, object]],
    *,
    kind: str,
    category: str,
    owner: str,
    filename: str,
    content: bytes,
    result: dict[str, object],
    batch_id: str = "",
) -> dict[str, object]:
    _prune_batch_import_precheck_cache(cache)
    if not result.get("can_import"):
        return result
    precheck_id = uuid.uuid4().hex
    cache[precheck_id] = {
        "kind": kind,
        "category": category,
        "owner": owner,
        "filename": filename,
        "content": content,
        "created_at": time.time(),
        "content_hash": hashlib.sha256(content).hexdigest(),
        "result": result,
        "batch_id": batch_id,
    }
    enriched = dict(result)
    enriched["precheck_id"] = precheck_id
    if batch_id:
        enriched["batch_id"] = batch_id
    enriched["precheck_ttl_seconds"] = BATCH_IMPORT_PRECHECK_TTL_SECONDS
    return enriched

def _load_batch_import_precheck(
    cache: dict[str, dict[str, object]],
    *,
    precheck_id: str,
    kind: str,
    category: str,
    owner: str,
) -> tuple[str, bytes, str] | None:
    normalized_id = str(precheck_id or "").strip()
    if not normalized_id:
        return None
    _prune_batch_import_precheck_cache(cache)
    entry = cache.get(normalized_id)
    if not entry:
        raise ValidationError("预检结果已过期，请重新选择文件完成预检")
    if str(entry.get("kind") or "") != kind:
        raise ValidationError("预检结果与当前导入类型不匹配，请重新预检")
    if str(entry.get("category") or "") != category:
        raise ValidationError("预检结果与当前数据库类型不匹配，请重新预检")
    if str(entry.get("owner") or "") != owner:
        raise PermissionError("不能使用其他用户生成的预检结果")
    result = entry.get("result")
    if not isinstance(result, dict) or not result.get("can_import"):
        raise ValidationError("预检未通过，不能进入正式导入")
    return str(entry.get("filename") or ""), bytes(entry.get("content") or b""), str(entry.get("batch_id") or "")

def _precheck_database_batch_upload(
    *,
    project_root: Path,
    sample_library_manager: SampleLibraryManager,
    filename: str,
    content: bytes,
) -> dict[str, object]:
    headers = _extract_batch_upload_headers(filename, content)
    rows = _parse_database_batch_upload(filename, content)
    required_columns = ["sample_name", "final_fasta_path"]
    metadata_templates = sample_library_manager.list_metadata_templates()
    valid_rows = 0
    empty_rows = 0
    issues: list[dict[str, str]] = []
    for index, row in enumerate(rows, start=2):
        sample_name = str(row.get("sample_name") or "").strip()
        final_fasta_path = str(row.get("final_fasta_path") or "").strip()
        if sample_name == "样本名称":
            continue
        if not sample_name and not final_fasta_path:
            empty_rows += 1
            continue
        row_errors: list[str] = []
        if not sample_name:
            row_errors.append("缺少 sample_name")
        if not final_fasta_path:
            row_errors.append("缺少 final_fasta_path")
        if final_fasta_path:
            try:
                _resolve_optional_existing_path(project_root, final_fasta_path)
            except Exception as exc:
                row_errors.append(str(exc))
        try:
            metadata_items = _build_database_import_metadata_items(row, metadata_templates)
            sample_library_manager.validate_metadata_items(metadata_items)
        except Exception as exc:
            row_errors.append(str(exc))
        if row_errors:
            issues.append({"row": str(index), "name": sample_name or "-", "reason": "；".join(row_errors)})
        else:
            valid_rows += 1
    return _build_precheck_summary(
        filename=filename,
        headers=headers,
        rows=rows,
        required_columns=required_columns,
        valid_rows=valid_rows,
        empty_rows=empty_rows,
        issues=issues,
    )

def _precheck_reference_batch_upload(
    *,
    project_root: Path,
    category: str,
    filename: str,
    content: bytes,
) -> dict[str, object]:
    headers = _extract_batch_upload_headers(filename, content)
    rows = _parse_database_batch_upload(filename, content)
    required_columns = ["host_name", "genome_name", "fasta_path"]
    valid_rows = 0
    empty_rows = 0
    issues: list[dict[str, str]] = []
    display_name = "宿主" if category == "host" else "病原"
    for index, row in enumerate(rows, start=2):
        host_name = str(row.get("host_name") or "").strip()
        genome_name = str(row.get("genome_name") or "").strip()
        fasta_path = str(row.get("fasta_path") or "").strip()
        if host_name in {"宿主名称", "病原名称", "物种名称"} and genome_name == "基因组名称":
            continue
        if not host_name and not genome_name and not fasta_path:
            empty_rows += 1
            continue
        row_errors: list[str] = []
        if not host_name:
            row_errors.append(f"缺少{display_name}名称")
        if not genome_name:
            row_errors.append("缺少 genome_name")
        if not fasta_path:
            row_errors.append("缺少 fasta_path")
        else:
            try:
                _resolve_optional_existing_path(project_root, fasta_path)
            except Exception as exc:
                row_errors.append(str(exc))
        if row_errors:
            issues.append({"row": str(index), "name": host_name or "-", "reason": "；".join(row_errors)})
        else:
            valid_rows += 1
    return _build_precheck_summary(
        filename=filename,
        headers=headers,
        rows=rows,
        required_columns=required_columns,
        valid_rows=valid_rows,
        empty_rows=empty_rows,
        issues=issues,
    )

def _normalize_batch_input_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    normalized_rows: list[dict[str, str]] = []
    for row in rows:
        normalized: dict[str, str] = {}
        for raw_key, raw_value in row.items():
            key = BATCH_INPUT_HEADER_MAP.get(str(raw_key or "").strip())
            if not key:
                continue
            normalized[key] = str(raw_value or "").strip()
        if normalized:
            normalized_rows.append(normalized)
    return normalized_rows
