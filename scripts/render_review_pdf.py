from __future__ import annotations

import re
import sys
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import Paragraph, SimpleDocTemplate, Spacer


FONT_PATHS = [
    "/System/Library/Fonts/STHeiti Light.ttc",
    "/System/Library/Fonts/STHeiti Medium.ttc",
    "/System/Library/Fonts/Supplemental/Arial Unicode.ttf",
]


def register_font() -> str:
    for path in FONT_PATHS:
        if Path(path).exists():
            try:
                pdfmetrics.registerFont(TTFont("ReviewCN", path))
                return "ReviewCN"
            except Exception:
                continue
    return "Helvetica"


def build_styles(font_name: str):
    styles = getSampleStyleSheet()
    styles.add(
        ParagraphStyle(
            name="ReviewTitle",
            parent=styles["Title"],
            fontName=font_name,
            fontSize=19,
            leading=26,
            alignment=TA_CENTER,
            textColor=colors.HexColor("#20354a"),
            spaceAfter=10,
        )
    )
    styles.add(
        ParagraphStyle(
            name="ReviewHeading1",
            parent=styles["Heading1"],
            fontName=font_name,
            fontSize=14,
            leading=20,
            textColor=colors.HexColor("#223c57"),
            spaceBefore=10,
            spaceAfter=8,
        )
    )
    styles.add(
        ParagraphStyle(
            name="ReviewBody",
            parent=styles["BodyText"],
            fontName=font_name,
            fontSize=10.5,
            leading=18,
            alignment=TA_JUSTIFY,
            textColor=colors.HexColor("#222222"),
            firstLineIndent=21,
            spaceAfter=6,
        )
    )
    styles.add(
        ParagraphStyle(
            name="ReviewRef",
            parent=styles["BodyText"],
            fontName=font_name,
            fontSize=9.5,
            leading=15,
            textColor=colors.HexColor("#333333"),
            spaceAfter=4,
        )
    )
    return styles


def escape_text(text: str) -> str:
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
    )


def inline_format(text: str) -> str:
    text = escape_text(text)
    text = re.sub(r"`([^`]+)`", r"<font face='Courier'>\1</font>", text)
    return text


def markdown_to_story(markdown_text: str, styles):
    story = []
    lines = markdown_text.splitlines()
    in_refs = False

    for raw_line in lines:
        line = raw_line.strip()
        if not line:
            story.append(Spacer(1, 3))
            continue
        if line.startswith("# "):
            story.append(Paragraph(inline_format(line[2:].strip()), styles["ReviewTitle"]))
            story.append(Spacer(1, 6))
            continue
        if line.startswith("## "):
            if line == "## 参考文献":
                in_refs = True
            story.append(Paragraph(inline_format(line[3:].strip()), styles["ReviewHeading1"]))
            continue
        if re.match(r"^\d+\.\s", line) and in_refs:
            story.append(Paragraph(inline_format(line), styles["ReviewRef"]))
            continue
        story.append(Paragraph(inline_format(line), styles["ReviewBody"]))

    return story


def main() -> int:
    if len(sys.argv) != 3:
        print("Usage: python3 render_review_pdf.py <input.md> <output.pdf>")
        return 1

    input_path = Path(sys.argv[1])
    output_path = Path(sys.argv[2])
    font_name = register_font()
    styles = build_styles(font_name)
    markdown_text = input_path.read_text(encoding="utf-8")
    story = markdown_to_story(markdown_text, styles)

    doc = SimpleDocTemplate(
        str(output_path),
        pagesize=A4,
        leftMargin=20 * mm,
        rightMargin=20 * mm,
        topMargin=18 * mm,
        bottomMargin=18 * mm,
        title="脑膜炎奈瑟菌基因组研究进展综述",
        author="OpenAI Codex",
    )
    doc.build(story)
    print(output_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
