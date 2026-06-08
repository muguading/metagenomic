from __future__ import annotations

import argparse
import csv
import html
import http.client
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from html.parser import HTMLParser
from pathlib import Path


CDC_URL = "https://calicivirustypingtool.cdc.gov/becerance.cgi"
EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def _request(url: str, retries: int = 5, sleep_seconds: float = 1.5) -> str:
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        request = urllib.request.Request(
            url,
            headers={"User-Agent": "metagenomic-cdc-norovirus-ref-downloader/1.0"},
        )
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                return response.read().decode("utf-8", errors="ignore")
        except (urllib.error.URLError, TimeoutError, http.client.IncompleteRead, ConnectionError) as exc:
            last_error = exc
            if attempt == retries:
                break
            print(f"[retry] request failed on attempt {attempt}/{retries}: {exc}", file=sys.stderr)
            time.sleep(sleep_seconds * attempt)
    assert last_error is not None
    raise last_error


def _batched(items: list[str], size: int) -> list[list[str]]:
    return [items[index:index + size] for index in range(0, len(items), size)]


def efetch_fasta(id_batch: list[str], email: str = "", api_key: str = "") -> str:
    params = {
        "db": "nuccore",
        "id": ",".join(id_batch),
        "rettype": "fasta",
        "retmode": "text",
    }
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key
    url = f"{EUTILS_BASE}/efetch.fcgi?{urllib.parse.urlencode(params)}"
    return _request(url)


class TableParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.tables: list[list[list[dict[str, object]]]] = []
        self._current_table: list[list[dict[str, object]]] | None = None
        self._current_row: list[dict[str, object]] | None = None
        self._current_cell: dict[str, object] | None = None
        self._capture_text = False
        self._capture_link = False

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        attr_map = dict(attrs)
        if tag == "table":
            self._current_table = []
        elif tag == "tr" and self._current_table is not None:
            self._current_row = []
        elif tag in {"td", "th"} and self._current_row is not None:
            self._current_cell = {
                "tag": tag,
                "text": "",
                "links": [],
                "rowspan": int(attr_map.get("rowspan") or "1"),
                "colspan": int(attr_map.get("colspan") or "1"),
            }
            self._capture_text = True
        elif tag == "a" and self._current_cell is not None:
            href = attr_map.get("href") or ""
            if href:
                cast_links = self._current_cell["links"]
                assert isinstance(cast_links, list)
                cast_links.append(href)
            self._capture_link = True

    def handle_endtag(self, tag: str) -> None:
        if tag in {"td", "th"} and self._current_cell is not None and self._current_row is not None:
            self._current_cell["text"] = html.unescape(str(self._current_cell["text"])).strip()
            self._current_row.append(self._current_cell)
            self._current_cell = None
            self._capture_text = False
            self._capture_link = False
        elif tag == "tr" and self._current_row is not None and self._current_table is not None:
            self._current_table.append(self._current_row)
            self._current_row = None
        elif tag == "table" and self._current_table is not None:
            self.tables.append(self._current_table)
            self._current_table = None

    def handle_data(self, data: str) -> None:
        if self._current_cell is not None and (self._capture_text or self._capture_link):
            current = str(self._current_cell["text"])
            self._current_cell["text"] = current + data


def _expand_rows(table_rows: list[list[dict[str, object]]]) -> list[list[dict[str, object]]]:
    expanded: list[list[dict[str, object]]] = []
    active_rowspans: list[tuple[int, dict[str, object]]] = []
    for row in table_rows:
        output_row: list[dict[str, object]] = []
        col = 0
        cell_iter = iter(row)
        next_cell = next(cell_iter, None)
        while next_cell is not None or any(span > 0 for span, _ in active_rowspans[col:]):
            while col < len(active_rowspans) and active_rowspans[col][0] > 0:
                span_left, span_cell = active_rowspans[col]
                output_row.append(span_cell)
                active_rowspans[col] = (span_left - 1, span_cell)
                col += 1
            if next_cell is None:
                continue
            cell = next_cell
            rowspan = int(cell.get("rowspan") or 1)
            colspan = int(cell.get("colspan") or 1)
            for _ in range(colspan):
                output_row.append(cell)
                while len(active_rowspans) <= col:
                    active_rowspans.append((0, {}))
                active_rowspans[col] = (rowspan - 1, cell)
                col += 1
            next_cell = next(cell_iter, None)
        expanded.append(output_row)
    return expanded


def _extract_accession(cell: dict[str, object]) -> str:
    links = cell.get("links") or []
    for link in links:
        matched = re.search(r"/nuccore/([A-Z0-9_.]+)", str(link), flags=re.IGNORECASE)
        if matched:
            return matched.group(1)
    text = str(cell.get("text") or "")
    matched = re.search(r"\b([A-Z]{1,4}_?\d+(?:\.\d+)?)\b", text)
    return matched.group(1) if matched else ""


def _normalize_subtype(text: str) -> str:
    return re.sub(r"\s*\(.*?\)\s*$", "", text.replace("\xa0", " ").strip())


def _parse_cdc_table(html_text: str) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    parser = TableParser()
    parser.feed(html_text)
    if not parser.tables:
        raise ValueError("No tables found in CDC page")

    reference_table: list[list[dict[str, object]]] | None = None
    for table in parser.tables:
        flattened = " ".join(str(cell.get("text") or "") for row in table for cell in row)
        if "Norovirus (262)" in flattened and "polymerase (RdRp)" in flattened and "capsid (VP1)" in flattened:
            reference_table = table
            break
    if reference_table is None:
        raise ValueError("Could not locate norovirus reference table")

    rows = _expand_rows(reference_table)
    rdrp_rows: list[dict[str, str]] = []
    vp1_rows: list[dict[str, str]] = []
    seen_rdrp: set[tuple[str, str]] = set()
    seen_vp1: set[tuple[str, str]] = set()

    for row in rows[2:]:
        if len(row) < 8:
            continue
        rdrp_subtype = _normalize_subtype(str(row[1].get("text") or ""))
        rdrp_acc = _extract_accession(row[2])
        rdrp_label = str(row[2].get("text") or "").strip()
        if rdrp_subtype and rdrp_acc:
            key = (rdrp_subtype, rdrp_acc)
            if key not in seen_rdrp:
                seen_rdrp.add(key)
                rdrp_rows.append(
                    {
                        "source": "CDC calicivirus typing tool",
                        "gene": "RdRp",
                        "subtype": rdrp_subtype,
                        "accession": rdrp_acc,
                        "label": rdrp_label,
                        "source_url": CDC_URL,
                    }
                )

        vp1_subtype = _normalize_subtype(str(row[5].get("text") or ""))
        vp1_acc = _extract_accession(row[6])
        vp1_label = str(row[6].get("text") or "").strip()
        if vp1_subtype and vp1_acc:
            key = (vp1_subtype, vp1_acc)
            if key not in seen_vp1:
                seen_vp1.add(key)
                vp1_rows.append(
                    {
                        "source": "CDC calicivirus typing tool",
                        "gene": "VP1",
                        "subtype": vp1_subtype,
                        "accession": vp1_acc,
                        "label": vp1_label,
                        "source_url": CDC_URL,
                    }
                )
    return rdrp_rows, vp1_rows


def _parse_fasta_records(fasta_text: str) -> dict[str, tuple[str, str]]:
    records: dict[str, tuple[str, str]] = {}
    header = ""
    chunks: list[str] = []

    def _flush() -> None:
        nonlocal header, chunks
        if not header:
            return
        accession = header.split()[0].lstrip(">")
        records[accession.split(".", 1)[0]] = (header.lstrip(">"), "".join(chunks))
        header = ""
        chunks = []

    for line in fasta_text.splitlines():
        if line.startswith(">"):
            _flush()
            header = line.strip()
            chunks = []
        elif header:
            chunks.append(line.strip())
    _flush()
    return records


def _write_fasta(path: Path, ordered_rows: list[dict[str, str]], sequences: dict[str, tuple[str, str]]) -> None:
    entries: list[str] = []
    for row in ordered_rows:
        accession = row["accession"]
        key = accession.split(".", 1)[0]
        if key not in sequences:
            continue
        original_header, sequence = sequences[key]
        header = f"{accession}_{row['subtype']}_{row['gene']} {original_header}"
        wrapped = "\n".join(sequence[index:index + 80] for index in range(0, len(sequence), 80))
        entries.append(f">{header}\n{wrapped}")
    path.write_text("\n".join(entries).strip() + "\n", encoding="utf-8")


def _write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["source", "gene", "subtype", "accession", "label", "source_url"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description="Download CDC norovirus typing reference sequences for RdRp and VP1")
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("database/virus/norovirus/cdc_typing_refs"),
        help="Output directory",
    )
    parser.add_argument("--email", default="", help="Optional email for NCBI E-utilities")
    parser.add_argument("--api-key", default="", help="Optional NCBI API key")
    args = parser.parse_args()

    out_dir = args.out_dir.expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    html_path = out_dir / "cdc_calicivirus_reference_sequences.html"
    rdrp_manifest = out_dir / "cdc_norovirus_rdrp_refs.tsv"
    vp1_manifest = out_dir / "cdc_norovirus_vp1_refs.tsv"
    rdrp_fasta = out_dir / "cdc_norovirus_rdrp_refs.fasta"
    vp1_fasta = out_dir / "cdc_norovirus_vp1_refs.fasta"

    html_text = _request(CDC_URL)
    html_path.write_text(html_text, encoding="utf-8")

    rdrp_rows, vp1_rows = _parse_cdc_table(html_text)
    _write_manifest(rdrp_manifest, rdrp_rows)
    _write_manifest(vp1_manifest, vp1_rows)

    all_accessions = sorted(
        {
            row["accession"].split(".", 1)[0]
            for row in rdrp_rows + vp1_rows
            if row.get("accession")
        }
    )
    sequences: dict[str, tuple[str, str]] = {}
    for batch in _batched(all_accessions, 200):
        fasta_text = efetch_fasta(batch, email=args.email, api_key=args.api_key)
        time.sleep(0.34)
        sequences.update(_parse_fasta_records(fasta_text))

    _write_fasta(rdrp_fasta, rdrp_rows, sequences)
    _write_fasta(vp1_fasta, vp1_rows, sequences)

    print(f"[ok] wrote {rdrp_manifest}", file=sys.stderr)
    print(f"[ok] wrote {vp1_manifest}", file=sys.stderr)
    print(f"[ok] wrote {rdrp_fasta}", file=sys.stderr)
    print(f"[ok] wrote {vp1_fasta}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
