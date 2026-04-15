#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FermFooDb ACE downloader v1（纯标准库版）
======================================

目标
----
从 FermFooDb 的 Browse -> Activity -> ACE-inhibitory 入口抓取：
1. 列表页 HTML / TXT
2. ACE 列表原始表 fermfoodb_ace_list.csv
3. 详情页 HTML / TXT
4. 详情页原始表 fermfoodb_ace_detail.csv
5. crawl_summary.json / failed_detail_pages.csv

默认入口
--------
https://webs.iiitd.edu.in/raghava/fermfoodb/browse_sub4.php?col=12&token=Ace-inhibitory

默认输出
--------
DB/raw/ace/databases/FermFooDb/experimental/
├── meta/
│   ├── crawl_summary.json
│   └── failed_detail_pages.csv
├── raw_html/
│   ├── list_pages/
│   └── detail_pages/
├── raw_txt/
│   ├── list_pages/
│   └── detail_cards/
└── raw_tables_csv/
    ├── fermfoodb_ace_list.csv
    └── fermfoodb_ace_detail.csv

说明
----
- 这是 raw 层抓取脚本，不做标准化
- 默认按列表页顺序抓全量
- 详情页字段尽量保守提取，提不到就留空
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import re
import socket
import time
from dataclasses import dataclass
from html.parser import HTMLParser
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from urllib.error import HTTPError, URLError
from urllib.parse import parse_qs, urljoin, urlparse, urlunparse
from urllib.request import Request, urlopen


BASE_URL = "https://webs.iiitd.edu.in/raghava/fermfoodb/"
START_URL = "https://webs.iiitd.edu.in/raghava/fermfoodb/browse_sub4.php?col=12&token=Ace-inhibitory"

USER_AGENT = (
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
    "AppleWebKit/537.36 (KHTML, like Gecko) "
    "Chrome/124.0 Safari/537.36"
)

LIST_COLUMNS = [
    "list_page_no",
    "list_page_url",
    "row_no_on_page",
    "fmdb_id",
    "detail_url",
    "pubmed_id",
    "title",
    "peptide_sequence",
    "length_of_peptide",
    "food_matrix",
    "protein",
    "ph",
    "temperature",
    "incubation_time",
    "activity",
    "experiment",
    "model",
    "assay_for_activity_measurement",
    "culture",
    "hydrolysis",
    "method_of_analysis",
    "mz_ratio",
    "mass",
    "ic50_raw",
]

DETAIL_COLUMNS = [
    "fmdb_id",
    "detail_url",
    "detail_title",
    "peptide_sequence",
    "length_of_peptide",
    "pubmed_id",
    "food_matrix",
    "protein_name",
    "ph",
    "temperature",
    "incubation_time",
    "activity",
    "experiment",
    "model",
    "assay_for_activity_measurement",
    "starter_culture",
    "hydrolysis",
    "method_of_analysis",
    "mz_ratio",
    "mass",
    "ic50_raw",
    "detail_text_path",
    "detail_html_path",
]

FAILED_DETAIL_COLUMNS = [
    "fmdb_id",
    "detail_url",
    "error_type",
    "error_message",
]

MISSING_TOKENS = {"", "-", "--", "na", "n/a", "none", "null"}


def infer_project_root(script_file: Path) -> Path:
    try:
        return script_file.resolve().parents[2]
    except IndexError:
        return script_file.resolve().parent


def normalize_spaces(text: str) -> str:
    text = text.replace("\ufeff", " ")
    text = text.replace("\r", " ").replace("\n", " ")
    return re.sub(r"\s+", " ", text).strip()


def clean_text(value: object) -> Optional[str]:
    if value is None:
        return None
    text = normalize_spaces(str(value))
    if text.lower() in MISSING_TOKENS or text == "":
        return None
    return text


def canonicalize_url(url: str) -> str:
    parsed = urlparse(url)
    scheme = parsed.scheme.lower()
    netloc = parsed.netloc.lower()
    path = parsed.path or "/"
    path = re.sub(r"/+", "/", path)
    query = "&".join(
        f"{k}={v}"
        for k, values in sorted(parse_qs(parsed.query, keep_blank_values=True).items())
        for v in sorted(values)
    )
    return urlunparse((scheme, netloc, path, "", query, ""))


def safe_filename(prefix: str, value: str) -> str:
    digest = hashlib.md5(value.encode("utf-8")).hexdigest()[:10]
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", value)[:60]
    value = value.strip("_") or "item"
    return f"{prefix}_{value}_{digest}"


def save_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def write_csv(path: Path, rows: List[Dict[str, object]], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def fetch_url(url: str, timeout: int = 20, sleep_sec: float = 0.0) -> Dict[str, object]:
    if sleep_sec > 0:
        time.sleep(sleep_sec)

    req = Request(
        url,
        headers={
            "User-Agent": USER_AGENT,
            "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
            "Accept-Language": "en-US,en;q=0.9",
            "Connection": "close",
        },
        method="GET",
    )

    started_at = time.time()
    try:
        with urlopen(req, timeout=timeout) as resp:
            content_type = resp.headers.get("Content-Type", "")
            final_url = resp.geturl()
            status = getattr(resp, "status", 200)
            raw = resp.read()
            elapsed = time.time() - started_at

            text = None
            encodings = []
            m = re.search(r"charset=([A-Za-z0-9._-]+)", content_type or "", flags=re.I)
            if m:
                encodings.append(m.group(1))
            encodings.extend(["utf-8", "utf-8-sig", "gb18030", "gbk", "latin-1"])

            for enc in encodings:
                try:
                    text = raw.decode(enc, errors="replace")
                    break
                except Exception:
                    continue
            if text is None:
                text = raw.decode("utf-8", errors="replace")

            return {
                "ok": True,
                "status": status,
                "final_url": final_url,
                "content_type": content_type,
                "elapsed_sec": round(elapsed, 4),
                "text": text,
                "error_type": "",
                "error_message": "",
            }

    except HTTPError as e:
        return {
            "ok": False,
            "status": getattr(e, "code", None),
            "final_url": url,
            "content_type": "",
            "elapsed_sec": round(time.time() - started_at, 4),
            "text": "",
            "error_type": "HTTPError",
            "error_message": str(e),
        }
    except URLError as e:
        return {
            "ok": False,
            "status": None,
            "final_url": url,
            "content_type": "",
            "elapsed_sec": round(time.time() - started_at, 4),
            "text": "",
            "error_type": "URLError",
            "error_message": str(e),
        }
    except socket.timeout as e:
        return {
            "ok": False,
            "status": None,
            "final_url": url,
            "content_type": "",
            "elapsed_sec": round(time.time() - started_at, 4),
            "text": "",
            "error_type": "Timeout",
            "error_message": str(e),
        }
    except Exception as e:
        return {
            "ok": False,
            "status": None,
            "final_url": url,
            "content_type": "",
            "elapsed_sec": round(time.time() - started_at, 4),
            "text": "",
            "error_type": e.__class__.__name__,
            "error_message": str(e),
        }


class TableAndLinkParser(HTMLParser):
    """
    轻量 HTML 解析器：
    - 提取页面纯文本
    - 提取链接
    - 提取 table -> rows -> cells
    """
    def __init__(self):
        super().__init__()
        self.links: List[Tuple[str, str]] = []
        self.text_parts: List[str] = []

        self.tables: List[List[List[str]]] = []
        self._in_table = False
        self._current_table: List[List[str]] = []
        self._in_row = False
        self._current_row: List[str] = []
        self._in_cell = False
        self._current_cell_parts: List[str] = []

        self._tag_stack: List[str] = []
        self._current_link_href: Optional[str] = None
        self._current_link_text_parts: List[str] = []

    def handle_starttag(self, tag, attrs):
        tag = tag.lower()
        self._tag_stack.append(tag)
        attrs_dict = dict(attrs)

        if tag == "a":
            self._current_link_href = attrs_dict.get("href")
            self._current_link_text_parts = []

        if tag == "table":
            self._in_table = True
            self._current_table = []

        if self._in_table and tag == "tr":
            self._in_row = True
            self._current_row = []

        if self._in_table and tag in ("td", "th"):
            self._in_cell = True
            self._current_cell_parts = []

    def handle_endtag(self, tag):
        tag = tag.lower()

        if tag == "a":
            if self._current_link_href:
                link_text = normalize_spaces(" ".join(self._current_link_text_parts))
                self.links.append((self._current_link_href, link_text))
            self._current_link_href = None
            self._current_link_text_parts = []

        if self._in_table and tag in ("td", "th") and self._in_cell:
            cell_text = normalize_spaces(" ".join(self._current_cell_parts))
            self._current_row.append(cell_text)
            self._in_cell = False
            self._current_cell_parts = []

        if self._in_table and tag == "tr" and self._in_row:
            if any(clean_text(x) for x in self._current_row):
                self._current_table.append(self._current_row)
            self._current_row = []
            self._in_row = False

        if tag == "table" and self._in_table:
            if self._current_table:
                self.tables.append(self._current_table)
            self._current_table = []
            self._in_table = False

        if self._tag_stack and self._tag_stack[-1] == tag:
            self._tag_stack.pop()
        else:
            for i in range(len(self._tag_stack) - 1, -1, -1):
                if self._tag_stack[i] == tag:
                    del self._tag_stack[i]
                    break

    def handle_data(self, data):
        text = normalize_spaces(html.unescape(data))
        if not text:
            return
        self.text_parts.append(text)
        if self._current_link_href is not None:
            self._current_link_text_parts.append(text)
        if self._in_cell:
            self._current_cell_parts.append(text)

    @property
    def plain_text(self) -> str:
        return normalize_spaces(" ".join(self.text_parts))


def parse_html(html_text: str) -> TableAndLinkParser:
    parser = TableAndLinkParser()
    parser.feed(html_text)
    return parser


def extract_total_entries(text: str) -> Optional[int]:
    patterns = [
        r"total number entries retrieved from this search are\s+(\d+)",
        r"total number of entries.*?(\d+)",
        r"entries retrieved.*?(\d+)",
    ]
    low = text.lower()
    for pat in patterns:
        m = re.search(pat, low, flags=re.I)
        if m:
            return int(m.group(1))
    return None


def choose_list_table(tables: List[List[List[str]]]) -> Optional[List[List[str]]]:
    """
    找最像列表结果的表。
    规则：表头里包含 FMDB_ID / Peptide / IC50 / Activity 任意若干关键字。
    """
    best = None
    best_score = -1
    for table in tables:
        if not table:
            continue
        header = " | ".join(table[0]).lower()
        score = 0
        for key in ["fmdb", "peptide", "sequence", "activity", "ic50", "food", "protein"]:
            if key in header:
                score += 1
        if score > best_score and len(table) >= 2:
            best = table
            best_score = score
    return best


def normalize_header_cell(text: str) -> str:
    text = clean_text(text) or ""
    text = text.replace(" ", "_")
    text = text.replace("-", "_")
    text = re.sub(r"[^A-Za-z0-9_]+", "", text)
    return text.lower()


def build_header_index(header_row: List[str]) -> Dict[str, int]:
    return {normalize_header_cell(x): i for i, x in enumerate(header_row)}


def get_cell(row: List[str], idx_map: Dict[str, int], *candidates: str) -> str:
    for cand in candidates:
        key = normalize_header_cell(cand)
        if key in idx_map:
            idx = idx_map[key]
            if idx < len(row):
                return row[idx]
    return ""


def extract_detail_url_from_links(links: List[Tuple[str, str]], fmdb_id: str) -> str:
    fmdb_id = clean_text(fmdb_id) or ""
    for href, link_text in links:
        abs_url = canonicalize_url(urljoin(BASE_URL, href))
        if "display_sub.php" in abs_url and fmdb_id and (fmdb_id in abs_url or fmdb_id == clean_text(link_text)):
            return abs_url
    # 退路：按 details 参数名找
    for href, _ in links:
        abs_url = canonicalize_url(urljoin(BASE_URL, href))
        if "display_sub.php" in abs_url and "details=" in abs_url:
            return abs_url
    return ""


def parse_list_page(
    page_no: int,
    page_url: str,
    html_text: str,
) -> Tuple[List[Dict[str, object]], Optional[int], List[str]]:
    parser = parse_html(html_text)
    total_entries = extract_total_entries(parser.plain_text)
    table = choose_list_table(parser.tables)
    warnings = []

    if table is None or len(table) < 2:
        warnings.append("list_table_not_found")
        return [], total_entries, warnings

    header = table[0]
    idx_map = build_header_index(header)
    rows: List[Dict[str, object]] = []

    for row_no, row in enumerate(table[1:], start=1):
        fmdb_id = get_cell(row, idx_map, "FMDB_ID", "FMDB ID", "ID")
        if not clean_text(fmdb_id):
            continue

        detail_url = extract_detail_url_from_links(parser.links, fmdb_id)

        rows.append({
            "list_page_no": page_no,
            "list_page_url": page_url,
            "row_no_on_page": row_no,
            "fmdb_id": clean_text(fmdb_id) or "",
            "detail_url": detail_url,
            "pubmed_id": clean_text(get_cell(row, idx_map, "PubMed ID", "PubMed_ID")) or "",
            "title": clean_text(get_cell(row, idx_map, "Title")) or "",
            "peptide_sequence": clean_text(get_cell(row, idx_map, "Peptide_Sequence", "Peptide Sequence")) or "",
            "length_of_peptide": clean_text(get_cell(row, idx_map, "Length of peptide", "Length_of_peptide")) or "",
            "food_matrix": clean_text(get_cell(row, idx_map, "Food_Matrix", "Food Matrix")) or "",
            "protein": clean_text(get_cell(row, idx_map, "Protein")) or "",
            "ph": clean_text(get_cell(row, idx_map, "pH", "PH")) or "",
            "temperature": clean_text(get_cell(row, idx_map, "Temperature")) or "",
            "incubation_time": clean_text(get_cell(row, idx_map, "Incubation Time", "Incubation_Time")) or "",
            "activity": clean_text(get_cell(row, idx_map, "Activity")) or "",
            "experiment": clean_text(get_cell(row, idx_map, "Experiment")) or "",
            "model": clean_text(get_cell(row, idx_map, "Model")) or "",
            "assay_for_activity_measurement": clean_text(get_cell(row, idx_map, "Assay for Activity Measurement", "Assay_for_Activity_Measurement")) or "",
            "culture": clean_text(get_cell(row, idx_map, "Culture")) or "",
            "hydrolysis": clean_text(get_cell(row, idx_map, "Hydrolysis")) or "",
            "method_of_analysis": clean_text(get_cell(row, idx_map, "Method of analysis", "Method_of_analysis")) or "",
            "mz_ratio": clean_text(get_cell(row, idx_map, "M_Z ratio", "M_Z_Ratio", "M/Z ratio")) or "",
            "mass": clean_text(get_cell(row, idx_map, "Mass")) or "",
            "ic50_raw": clean_text(get_cell(row, idx_map, "IC50")) or "",
        })

    return rows, total_entries, warnings


def next_page_candidates(current_url: str, page_no: int) -> List[str]:
    """
    尝试构造分页链接。
    FermFooDb 常见分页可能带 page / pageno / start 等参数。
    这里保守尝试几类。
    """
    parsed = urlparse(current_url)
    q = parse_qs(parsed.query, keep_blank_values=True)

    candidates = []

    for key in ["page", "pageno", "Page", "PageNo"]:
        q2 = {k: list(v) for k, v in q.items()}
        q2[key] = [str(page_no + 1)]
        query = "&".join(f"{k}={v}" for k, vals in sorted(q2.items()) for v in vals)
        candidates.append(urlunparse((parsed.scheme, parsed.netloc, parsed.path, "", query, "")))

    for key in ["start", "offset"]:
        q2 = {k: list(v) for k, v in q.items()}
        q2[key] = [str((page_no) * 50)]
        query = "&".join(f"{k}={v}" for k, vals in sorted(q2.items()) for v in vals)
        candidates.append(urlunparse((parsed.scheme, parsed.netloc, parsed.path, "", query, "")))

    # 也尝试直接附加 ?page=
    if "page" not in q:
        q2 = {k: list(v) for k, v in q.items()}
        q2["page"] = [str(page_no + 1)]
        query = "&".join(f"{k}={v}" for k, vals in sorted(q2.items()) for v in vals)
        candidates.append(urlunparse((parsed.scheme, parsed.netloc, parsed.path, "", query, "")))

    dedup = []
    seen = set()
    for c in candidates:
        cc = canonicalize_url(c)
        if cc not in seen:
            seen.add(cc)
            dedup.append(cc)
    return dedup


def discover_explicit_next_link(links: List[Tuple[str, str]], current_url: str) -> Optional[str]:
    keywords = {"next", ">", ">>", "next page"}
    for href, text in links:
        link_text = (clean_text(text) or "").lower()
        if link_text in keywords or "next" in link_text:
            return canonicalize_url(urljoin(current_url, href))
    return None


def parse_detail_pairs(plain_text: str) -> Dict[str, str]:
    """
    详情页通常是标签: 值 的成对结构。
    这里做保守抽取。
    """
    fields = {
        "detail_title": "",
        "peptide_sequence": "",
        "length_of_peptide": "",
        "pubmed_id": "",
        "food_matrix": "",
        "protein_name": "",
        "ph": "",
        "temperature": "",
        "incubation_time": "",
        "activity": "",
        "experiment": "",
        "model": "",
        "assay_for_activity_measurement": "",
        "starter_culture": "",
        "hydrolysis": "",
        "method_of_analysis": "",
        "mz_ratio": "",
        "mass": "",
        "ic50_raw": "",
    }

    patterns = {
        "detail_title": [r"Title\s*:\s*(.+?)\s*(?:PubMed ID|Peptide Sequence|$)"],
        "peptide_sequence": [r"Peptide Sequence\s*:\s*(.+?)\s*(?:Length of peptide|PubMed ID|Food Matrix|$)"],
        "length_of_peptide": [r"Length of peptide\s*:\s*(.+?)\s*(?:Food Matrix|Protein Name|PubMed ID|$)"],
        "pubmed_id": [r"PubMed ID\s*:\s*(.+?)\s*(?:Title|Peptide Sequence|Food Matrix|$)"],
        "food_matrix": [r"Food Matrix\s*:\s*(.+?)\s*(?:Protein Name|pH|Temperature|$)"],
        "protein_name": [r"Protein Name\s*:\s*(.+?)\s*(?:pH|Temperature|Incubation Time|$)"],
        "ph": [r"pH\s*:\s*(.+?)\s*(?:Temperature|Incubation Time|Activity|$)"],
        "temperature": [r"Temperature\s*:\s*(.+?)\s*(?:Incubation Time|Activity|Experiment|$)"],
        "incubation_time": [r"Incubation Time\s*:\s*(.+?)\s*(?:Activity|Experiment|Model|$)"],
        "activity": [r"Activity\s*:\s*(.+?)\s*(?:Experiment|Model|Assay for Activity Measurement|$)"],
        "experiment": [r"Experiment\s*:\s*(.+?)\s*(?:Model|Assay for Activity Measurement|Starter Culture|$)"],
        "model": [r"Model\s*:\s*(.+?)\s*(?:Assay for Activity Measurement|Starter Culture|Hydrolysis|$)"],
        "assay_for_activity_measurement": [r"Assay for Activity Measurement\s*:\s*(.+?)\s*(?:Starter Culture|Hydrolysis|Method of Analysis|$)"],
        "starter_culture": [r"(?:Starter Culture|Culture)\s*:\s*(.+?)\s*(?:Hydrolysis|Method of Analysis|M/Z ratio|$)"],
        "hydrolysis": [r"Hydrolysis\s*:\s*(.+?)\s*(?:Method of Analysis|M/Z ratio|Mass|$)"],
        "method_of_analysis": [r"Method of Analysis\s*:\s*(.+?)\s*(?:M/Z ratio|M_Z ratio|Mass|IC50|$)"],
        "mz_ratio": [r"(?:M/Z ratio|M_Z ratio)\s*:\s*(.+?)\s*(?:Mass|IC50|$)"],
        "mass": [r"Mass\s*:\s*(.+?)\s*(?:IC50|$)"],
        "ic50_raw": [r"IC50\s*:\s*(.+?)\s*$"],
    }

    text = plain_text
    for key, pats in patterns.items():
        for pat in pats:
            m = re.search(pat, text, flags=re.I)
            if m:
                fields[key] = normalize_spaces(m.group(1))
                break

    return fields


def parse_detail_page(fmdb_id: str, detail_url: str, html_text: str) -> Dict[str, object]:
    parser = parse_html(html_text)
    plain_text = parser.plain_text
    pairs = parse_detail_pairs(plain_text)

    return {
        "fmdb_id": fmdb_id,
        "detail_url": detail_url,
        "detail_title": clean_text(pairs["detail_title"]) or "",
        "peptide_sequence": clean_text(pairs["peptide_sequence"]) or "",
        "length_of_peptide": clean_text(pairs["length_of_peptide"]) or "",
        "pubmed_id": clean_text(pairs["pubmed_id"]) or "",
        "food_matrix": clean_text(pairs["food_matrix"]) or "",
        "protein_name": clean_text(pairs["protein_name"]) or "",
        "ph": clean_text(pairs["ph"]) or "",
        "temperature": clean_text(pairs["temperature"]) or "",
        "incubation_time": clean_text(pairs["incubation_time"]) or "",
        "activity": clean_text(pairs["activity"]) or "",
        "experiment": clean_text(pairs["experiment"]) or "",
        "model": clean_text(pairs["model"]) or "",
        "assay_for_activity_measurement": clean_text(pairs["assay_for_activity_measurement"]) or "",
        "starter_culture": clean_text(pairs["starter_culture"]) or "",
        "hydrolysis": clean_text(pairs["hydrolysis"]) or "",
        "method_of_analysis": clean_text(pairs["method_of_analysis"]) or "",
        "mz_ratio": clean_text(pairs["mz_ratio"]) or "",
        "mass": clean_text(pairs["mass"]) or "",
        "ic50_raw": clean_text(pairs["ic50_raw"]) or "",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="FermFooDb ACE downloader v1（纯标准库版）")
    parser.add_argument("--project-root", type=str, default=None, help="PepDB 项目根目录")
    parser.add_argument("--start-url", type=str, default=START_URL, help="ACE 列表入口 URL")
    parser.add_argument("--timeout", type=int, default=20, help="请求超时秒数")
    parser.add_argument("--sleep-sec", type=float, default=0.3, help="请求间隔秒数")
    parser.add_argument("--max-pages", type=int, default=100, help="最多抓多少个列表页")
    parser.add_argument("--max-detail", type=int, default=None, help="最多抓多少个详情页，默认抓全量")
    args = parser.parse_args()

    script_file = Path(__file__)
    project_root = Path(args.project_root) if args.project_root else infer_project_root(script_file)

    base_dir = project_root / "DB" / "raw" / "ace" / "databases" / "FermFooDb" / "experimental"
    meta_dir = base_dir / "meta"
    list_html_dir = base_dir / "raw_html" / "list_pages"
    detail_html_dir = base_dir / "raw_html" / "detail_pages"
    list_txt_dir = base_dir / "raw_txt" / "list_pages"
    detail_txt_dir = base_dir / "raw_txt" / "detail_cards"
    raw_tables_dir = base_dir / "raw_tables_csv"

    crawl_summary_path = meta_dir / "crawl_summary.json"
    failed_detail_path = meta_dir / "failed_detail_pages.csv"
    list_csv_path = raw_tables_dir / "fermfoodb_ace_list.csv"
    detail_csv_path = raw_tables_dir / "fermfoodb_ace_detail.csv"

    meta_dir.mkdir(parents=True, exist_ok=True)
    raw_tables_dir.mkdir(parents=True, exist_ok=True)

    # 1) 抓列表页
    list_rows: List[Dict[str, object]] = []
    visited_list_urls = set()
    failed_list_urls: List[Dict[str, str]] = []
    total_entries_reported: Optional[int] = None

    current_url = canonicalize_url(args.start_url)
    page_no = 1

    while current_url and page_no <= args.max_pages and current_url not in visited_list_urls:
        result = fetch_url(current_url, timeout=args.timeout, sleep_sec=args.sleep_sec)
        visited_list_urls.add(current_url)

        if not result["ok"]:
            failed_list_urls.append({
                "page_no": str(page_no),
                "url": current_url,
                "error_type": result["error_type"],
                "error_message": result["error_message"],
            })
            break

        page_html = result["text"]
        page_url = canonicalize_url(result["final_url"])
        parser_obj = parse_html(page_html)
        plain_text = parser_obj.plain_text

        html_path = list_html_dir / f"list_page_{page_no:03d}.html"
        txt_path = list_txt_dir / f"list_page_{page_no:03d}.txt"
        save_text(html_path, page_html)
        save_text(txt_path, plain_text)

        parsed_rows, total_entries_found, warnings = parse_list_page(page_no, page_url, page_html)
        if total_entries_reported is None and total_entries_found is not None:
            total_entries_reported = total_entries_found

        # 如果这一页没有解析到任何记录，停止
        if not parsed_rows:
            break

        list_rows.extend(parsed_rows)

        explicit_next = discover_explicit_next_link(parser_obj.links, page_url)
        next_url = None

        if explicit_next and explicit_next not in visited_list_urls:
            next_url = explicit_next
        else:
            for cand in next_page_candidates(page_url, page_no):
                if cand not in visited_list_urls:
                    # 试探下一页是否与当前不同
                    next_url = cand
                    break

        current_url = next_url
        page_no += 1

    write_csv(list_csv_path, list_rows, LIST_COLUMNS)

    # 2) 抓详情页
    detail_rows: List[Dict[str, object]] = []
    failed_detail_rows: List[Dict[str, object]] = []

    detail_seen = set()
    detail_targets = [
        (clean_text(r.get("fmdb_id")) or "", clean_text(r.get("detail_url")) or "")
        for r in list_rows
        if clean_text(r.get("fmdb_id")) and clean_text(r.get("detail_url"))
    ]

    if args.max_detail is not None:
        detail_targets = detail_targets[: args.max_detail]

    for idx, (fmdb_id, detail_url) in enumerate(detail_targets, start=1):
        detail_url = canonicalize_url(detail_url)
        if detail_url in detail_seen:
            continue
        detail_seen.add(detail_url)

        result = fetch_url(detail_url, timeout=args.timeout, sleep_sec=args.sleep_sec)
        if not result["ok"]:
            failed_detail_rows.append({
                "fmdb_id": fmdb_id,
                "detail_url": detail_url,
                "error_type": result["error_type"],
                "error_message": result["error_message"],
            })
            continue

        detail_html = result["text"]
        parser_obj = parse_html(detail_html)
        plain_text = parser_obj.plain_text

        stem = safe_filename("detail", fmdb_id or detail_url)
        html_path = detail_html_dir / f"{stem}.html"
        txt_path = detail_txt_dir / f"{stem}.txt"
        save_text(html_path, detail_html)
        save_text(txt_path, plain_text)

        detail_row = parse_detail_page(fmdb_id, detail_url, detail_html)
        detail_row["detail_text_path"] = str(txt_path)
        detail_row["detail_html_path"] = str(html_path)
        detail_rows.append(detail_row)

    write_csv(detail_csv_path, detail_rows, DETAIL_COLUMNS)
    write_csv(failed_detail_path, failed_detail_rows, FAILED_DETAIL_COLUMNS)

    summary = {
        "source_name": "FermFooDb",
        "task": "ace",
        "start_url": args.start_url,
        "list_pages_fetched": len(visited_list_urls),
        "list_rows_parsed": len(list_rows),
        "total_entries_reported_on_site": total_entries_reported,
        "detail_targets": len(detail_targets),
        "detail_rows_parsed": len(detail_rows),
        "failed_detail_count": len(failed_detail_rows),
        "failed_list_count": len(failed_list_urls),
        "list_csv_path": str(list_csv_path),
        "detail_csv_path": str(detail_csv_path),
    }
    crawl_summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")

    print("=" * 80)
    print("FermFooDb ACE downloader v1 完成")
    print(f"项目根目录：{project_root}")
    print(f"列表页抓取数：{len(visited_list_urls)}")
    print(f"列表记录数：{len(list_rows)}")
    print(f"站点报告总条数：{total_entries_reported}")
    print(f"详情页目标数：{len(detail_targets)}")
    print(f"详情页成功数：{len(detail_rows)}")
    print(f"详情页失败数：{len(failed_detail_rows)}")
    print(f"列表 CSV：{list_csv_path}")
    print(f"详情 CSV：{detail_csv_path}")
    print(f"摘要 JSON：{crawl_summary_path}")
    print("=" * 80)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())