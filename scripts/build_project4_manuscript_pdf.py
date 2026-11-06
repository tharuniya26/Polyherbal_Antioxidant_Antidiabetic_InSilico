"""Build the Project 4 supervisor-review manuscript PDF from its Markdown source.

The renderer intentionally supports only the Markdown features used by the
manuscript: headings, paragraphs, emphasis, links, pipe tables, numbered and
bulleted lists, and local PNG figures. It does not alter scientific content.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from html import escape
from pathlib import Path
import re

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    Image,
    KeepTogether,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = ROOT / "manuscript" / "project4_insilico_manuscript_draft.md"
DEFAULT_OUTPUT = ROOT / "manuscript" / "Tharuniya_Anpalakan_Polyherbal_Computational_Manuscript.pdf"
SHORT_TITLE = "Exploratory phytochemical pose geometry at human ntMGAM"


@dataclass(frozen=True)
class FontFamily:
    regular: str
    bold: str
    italic: str
    bold_italic: str


def _register_serif_fonts() -> FontFamily:
    """Register Times New Roman when available, with a Unicode serif fallback."""
    candidates = [
        (
            Path("C:/Windows/Fonts/times.ttf"),
            Path("C:/Windows/Fonts/timesbd.ttf"),
            Path("C:/Windows/Fonts/timesi.ttf"),
            Path("C:/Windows/Fonts/timesbi.ttf"),
        ),
        (
            Path("/usr/share/fonts/truetype/dejavu/DejaVuSerif.ttf"),
            Path("/usr/share/fonts/truetype/dejavu/DejaVuSerif-Bold.ttf"),
            Path("/usr/share/fonts/truetype/dejavu/DejaVuSerif-Italic.ttf"),
            Path("/usr/share/fonts/truetype/dejavu/DejaVuSerif-BoldItalic.ttf"),
        ),
    ]
    for regular, bold, italic, bold_italic in candidates:
        if all(path.is_file() for path in (regular, bold, italic, bold_italic)):
            names = FontFamily("ManuscriptSerif", "ManuscriptSerif-Bold", "ManuscriptSerif-Italic", "ManuscriptSerif-BoldItalic")
            for name, path in zip(names.__dict__.values(), (regular, bold, italic, bold_italic), strict=True):
                pdfmetrics.registerFont(TTFont(name, str(path)))
            pdfmetrics.registerFontFamily(
                "ManuscriptSerif",
                normal=names.regular,
                bold=names.bold,
                italic=names.italic,
                boldItalic=names.bold_italic,
            )
            return names
    raise FileNotFoundError("No supported Unicode serif font family was found")


def _inline(text: str) -> str:
    """Translate the manuscript's inline Markdown subset to ReportLab markup."""
    text = escape(text, quote=True)
    text = re.sub(
        r"\[([^\]]+)\]\((https?://[^)]+)\)",
        lambda match: f'<a href="{match.group(2)}" color="#1f4e79">{match.group(1)}</a>',
        text,
    )
    text = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", text)
    text = re.sub(r"(?<!\*)\*([^*]+?)\*(?!\*)", r"<i>\1</i>", text)
    text = re.sub(r"`([^`]+)`", r'<font name="Courier">\1</font>', text)
    return text


def _styles(fonts: FontFamily) -> dict[str, ParagraphStyle]:
    base = getSampleStyleSheet()
    return {
        "title": ParagraphStyle(
            "ManuscriptTitle",
            parent=base["Title"],
            fontName=fonts.bold,
            fontSize=17,
            leading=21,
            alignment=TA_CENTER,
            textColor=colors.HexColor("#16324F"),
            spaceAfter=5 * mm,
        ),
        "author": ParagraphStyle(
            "Author",
            parent=base["Normal"],
            fontName=fonts.regular,
            fontSize=10.5,
            leading=14,
            alignment=TA_CENTER,
            spaceAfter=4 * mm,
        ),
        "h1": ParagraphStyle(
            "H1",
            parent=base["Heading1"],
            fontName=fonts.bold,
            fontSize=14,
            leading=17,
            textColor=colors.HexColor("#16324F"),
            spaceBefore=5 * mm,
            spaceAfter=2.5 * mm,
            keepWithNext=True,
        ),
        "h2": ParagraphStyle(
            "H2",
            parent=base["Heading2"],
            fontName=fonts.bold,
            fontSize=11.5,
            leading=14,
            textColor=colors.HexColor("#28536B"),
            spaceBefore=3.5 * mm,
            spaceAfter=1.5 * mm,
            keepWithNext=True,
        ),
        "body": ParagraphStyle(
            "Body",
            parent=base["BodyText"],
            fontName=fonts.regular,
            fontSize=9.5,
            leading=13,
            alignment=TA_JUSTIFY,
            firstLineIndent=5 * mm,
            spaceAfter=2.5 * mm,
            allowWidows=0,
            allowOrphans=0,
        ),
        "body_code": ParagraphStyle(
            "BodyWithCode",
            parent=base["BodyText"],
            fontName=fonts.regular,
            fontSize=9.2,
            leading=12.8,
            alignment=TA_LEFT,
            firstLineIndent=5 * mm,
            spaceAfter=2.5 * mm,
            allowWidows=0,
            allowOrphans=0,
        ),
        "abstract": ParagraphStyle(
            "AbstractBody",
            parent=base["BodyText"],
            fontName=fonts.regular,
            fontSize=9,
            leading=12.2,
            alignment=TA_JUSTIFY,
            leftIndent=7 * mm,
            rightIndent=7 * mm,
            spaceAfter=2 * mm,
        ),
        "list": ParagraphStyle(
            "List",
            parent=base["BodyText"],
            fontName=fonts.regular,
            fontSize=9.2,
            leading=12.4,
            leftIndent=7 * mm,
            firstLineIndent=-4 * mm,
            spaceAfter=1.5 * mm,
        ),
        "reference": ParagraphStyle(
            "Reference",
            parent=base["BodyText"],
            fontName=fonts.regular,
            fontSize=8.2,
            leading=10.6,
            leftIndent=6 * mm,
            firstLineIndent=-6 * mm,
            alignment=TA_LEFT,
            spaceAfter=1.2 * mm,
        ),
        "caption": ParagraphStyle(
            "Caption",
            parent=base["BodyText"],
            fontName=fonts.regular,
            fontSize=8.2,
            leading=10.5,
            alignment=TA_JUSTIFY,
            leftIndent=4 * mm,
            rightIndent=4 * mm,
            spaceBefore=1.2 * mm,
            spaceAfter=3.5 * mm,
        ),
        "table": ParagraphStyle(
            "TableCell",
            parent=base["BodyText"],
            fontName=fonts.regular,
            fontSize=7.4,
            leading=9.3,
            alignment=TA_LEFT,
        ),
        "table_header": ParagraphStyle(
            "TableHeader",
            parent=base["BodyText"],
            fontName=fonts.bold,
            fontSize=7.4,
            leading=9.3,
            textColor=colors.white,
            alignment=TA_LEFT,
        ),
    }


def _pipe_cells(line: str) -> list[str]:
    return [cell.strip() for cell in line.strip().strip("|").split("|")]


def _table(lines: list[str], styles: dict[str, ParagraphStyle], available_width: float) -> Table:
    raw = [_pipe_cells(line) for line in lines]
    raw = [raw[0], *raw[2:]]  # remove Markdown alignment separator
    column_count = max(len(row) for row in raw)
    rows = []
    for row_index, row in enumerate(raw):
        row += [""] * (column_count - len(row))
        style = styles["table_header"] if row_index == 0 else styles["table"]
        rows.append([Paragraph(_inline(cell), style) for cell in row])
    weights = [1.0] * column_count
    if column_count == 3:
        weights = [0.23, 0.29, 0.48]
    elif column_count == 4:
        weights = [0.20, 0.20, 0.24, 0.36]
    total = sum(weights)
    widths = [available_width * weight / total for weight in weights]
    table = Table(rows, colWidths=widths, repeatRows=1, hAlign="CENTER")
    table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#28536B")),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("GRID", (0, 0), (-1, -1), 0.35, colors.HexColor("#9AA6AE")),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F2F5F7")]),
                ("LEFTPADDING", (0, 0), (-1, -1), 4),
                ("RIGHTPADDING", (0, 0), (-1, -1), 4),
                ("TOPPADDING", (0, 0), (-1, -1), 3),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
            ]
        )
    )
    return table


def _figure(path: Path, max_width: float, max_height: float) -> Image:
    image = Image(str(path))
    scale = min(max_width / image.imageWidth, max_height / image.imageHeight)
    image.drawWidth = image.imageWidth * scale
    image.drawHeight = image.imageHeight * scale
    image.hAlign = "CENTER"
    return image


def _footer(canvas, document) -> None:  # ReportLab callback signature
    canvas.saveState()
    width, _ = A4
    canvas.setStrokeColor(colors.HexColor("#B7C2C9"))
    canvas.setLineWidth(0.35)
    canvas.line(22 * mm, 15 * mm, width - 22 * mm, 15 * mm)
    canvas.setFont("Helvetica", 7.5)
    canvas.setFillColor(colors.HexColor("#58646C"))
    canvas.drawString(22 * mm, 10.5 * mm, "Supervisor-review draft • Not peer-reviewed")
    canvas.drawRightString(width - 22 * mm, 10.5 * mm, f"Page {canvas.getPageNumber()}")
    canvas.restoreState()


def _story(source: Path, styles: dict[str, ParagraphStyle], available_width: float) -> list:
    lines = source.read_text(encoding="utf-8").splitlines()
    story: list = []
    in_abstract = False
    in_references = False
    index = 0
    while index < len(lines):
        line = lines[index].rstrip()
        stripped = line.strip()
        if not stripped:
            index += 1
            continue

        if stripped.startswith("# "):
            story.append(Paragraph(_inline(stripped[2:]), styles["title"]))
            index += 1
            author_lines: list[str] = []
            while index < len(lines) and not lines[index].startswith("## "):
                if lines[index].strip():
                    author_lines.append(lines[index].strip().rstrip("  "))
                index += 1
            if author_lines:
                story.append(Paragraph("<br/>".join(_inline(item) for item in author_lines), styles["author"]))
            continue

        if stripped.startswith("## "):
            heading = stripped[3:]
            in_abstract = heading == "Abstract"
            in_references = heading == "References"
            story.append(Paragraph(_inline(heading), styles["h1"]))
            index += 1
            continue

        if stripped.startswith("### "):
            story.append(Paragraph(_inline(stripped[4:]), styles["h2"]))
            index += 1
            continue

        image_match = re.fullmatch(r"!\[([^\]]*)\]\(([^)]+)\)", stripped)
        if image_match:
            image_path = (source.parent / image_match.group(2)).resolve()
            if not image_path.is_file():
                raise FileNotFoundError(f"Figure not found: {image_path}")
            figure = _figure(image_path, available_width, 102 * mm)
            next_index = index + 1
            while next_index < len(lines) and not lines[next_index].strip():
                next_index += 1
            grouped = [figure]
            if next_index < len(lines) and lines[next_index].lstrip().startswith("**Figure "):
                grouped.append(Paragraph(_inline(lines[next_index].strip()), styles["caption"]))
                index = next_index + 1
            else:
                index += 1
            story.append(KeepTogether(grouped))
            continue

        if stripped.startswith("|") and index + 1 < len(lines) and re.match(r"^\s*\|?\s*:?-+", lines[index + 1]):
            table_lines = [stripped, lines[index + 1].strip()]
            index += 2
            while index < len(lines) and lines[index].strip().startswith("|"):
                table_lines.append(lines[index].strip())
                index += 1
            story.extend([_table(table_lines, styles, available_width), Spacer(1, 3 * mm)])
            continue

        numbered = re.match(r"^(\d+)\.\s+(.+)$", stripped)
        bullet = re.match(r"^[-*]\s+(.+)$", stripped)
        if numbered:
            style = styles["reference"] if in_references else styles["list"]
            story.append(Paragraph(f"{numbered.group(1)}. {_inline(numbered.group(2))}", style))
            index += 1
            continue
        if bullet:
            story.append(Paragraph(f"• {_inline(bullet.group(1))}", styles["list"]))
            index += 1
            continue

        paragraph_lines = [stripped]
        index += 1
        while index < len(lines):
            candidate = lines[index].strip()
            if not candidate:
                break
            if candidate.startswith(("#", "|", "![", "- ", "* ")) or re.match(r"^\d+\.\s+", candidate):
                break
            paragraph_lines.append(candidate)
            index += 1
        text = " ".join(paragraph_lines)
        if in_abstract:
            style = styles["abstract"]
        elif "`" in text:
            style = styles["body_code"]
        else:
            style = styles["body"]
        story.append(Paragraph(_inline(text), style))

    return story


def build(source: Path, output: Path) -> None:
    fonts = _register_serif_fonts()
    styles = _styles(fonts)
    output.parent.mkdir(parents=True, exist_ok=True)
    document = SimpleDocTemplate(
        str(output),
        pagesize=A4,
        rightMargin=22 * mm,
        leftMargin=22 * mm,
        topMargin=20 * mm,
        bottomMargin=20 * mm,
        title="Reproducible pose-geometry evaluation of six literature-associated phytochemicals",
        author="Tharuniya Anpalakan",
        subject="Supervisor-review computational research manuscript; not peer-reviewed",
        creator="Project 4 reproducible ReportLab renderer",
        pageCompression=1,
    )
    story = _story(source, styles, A4[0] - 44 * mm)
    document.build(story, onFirstPage=_footer, onLaterPages=_footer)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    build(args.source.resolve(), args.output.resolve())
    print(args.output.resolve())


if __name__ == "__main__":
    main()
