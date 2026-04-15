#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FermFooDb ACE downloader v2（纯标准库版）
======================================

目标
----
从 FermFooDb 的 Browse -> Activity -> ACE-inhibitory 入口抓取：
1. 列表页 HTML / TXT
2. ACE 列表原始表 fermfoodb_ace_list.csv
3. 详情页 HTML / TXT
4. 详情页原始表 fermfoodb_ace_detail.csv
5. crawl_summary.json / failed_detail_pages.csv

和 v1 的区别
------------
- 默认只抓 ACE 结果入口这一页，不自动翻页
- 按 FMDB_ID 去重
- 详情页链接直接按 FMDB_ID 构造，避免 detail_targets=0
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import re
import socket
import time
from collections import OrderedDict
from html.parser import HTMLParser
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from urllib.error import HTTPError, URLError
from urllib.parse import urljoin
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
    "pubmed_id",
    "peptide_sequence",
    "length_of_peptide",
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


class TableParser(HTMLParser):
    """
    轻量解析：
    - 页面纯文本
    - 所有表格
    """

    def __init__(self):
        super().__init__()
        self.text_parts: List[str] = []

        self.tables: List[List[List[str]]] = []
        self._in_table = False
        self._current_table: List[List[str]] = []

        self._in_row = False
        self._current_row: List[str] = []

        self._in_cell = False
        self._current_cell_parts: List[str] = []

    def handle_starttag(self, tag, attrs):
        tag = tag.lower()

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

    def handle_data(self, data):
        text = normalize_spaces(html.unescape(data))
        if not text:
            return
        self.text_parts.append(text)
        if self._in_cell:
            self._current_cell_parts.append(text)

    @property
    def plain_text(self) -> str:
        return normalize_spaces(" ".join(self.text_parts))


def parse_html_tables(html_text: str) -> TableParser:
    parser = TableParser()
    parser.feed(html_text)
    return parser


def extract_total_entries(text: str) -> Optional[int]:
    low = text.lower()
    patterns = [
        r"total number entries retrieved from this search are\s+(\d+)",
        r"total number of entries.*?(\d+)",
        r"entries retrieved.*?(\d+)",
    ]
    for pat in patterns:
        m = re.search(pat, low, flags=re.I)
        if m:
            return int(m.group(1))
    return None


def choose_list_table(tables: List[List[List[str]]]) -> Optional[List[List[str]]]:
    best = None
    best_score = -1
    for table in tables:
        if not table or len(table) < 2:
            continue
        header = " | ".join(table[0]).lower()
        score = 0
        for key in ["fmdb", "peptide", "sequence", "activity", "ic50", "food", "protein"]:
            if key in header:
                score += 1
        if score > best_score:
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


def fmdb_id_to_detail_url(fmdb_id: str) -> str:
    return urljoin(BASE_URL, f"display_sub.php?details={fmdb_id}")


def parse_list_page(page_no: int, page_url: str, html_text: str) -> Tuple[List[Dict[str, object]], Optional[int], List[str]]:
    parser = parse_html_tables(html_text)
    total_entries = extract_total_entries(parser.plain_text)
    table = choose_list_table(parser.tables)
    warnings: List[str] = []

    if table is None:
        warnings.append("list_table_not_found")
        return [], total_entries, warnings

    header = table[0]
    idx_map = build_header_index(header)

    rows: List[Dict[str, object]] = []
    for row_no, row in enumerate(table[1:], start=1):
        fmdb_id = clean_text(get_cell(row, idx_map, "FMDB_ID", "FMDB ID", "ID")) or ""
        if not fmdb_id:
            continue

        rows.append({
            "list_page_no": page_no,
            "list_page_url": page_url,
            "row_no_on_page": row_no,
            "fmdb_id": fmdb_id,
            "detail_url": fmdb_id_to_detail_url(fmdb_id),
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
            "assay_for_activity_measurement": clean_text(
                get_cell(row, idx_map, "Assay for Activity Measurement", "Assay_for_Activity_Measurement")
            ) or "",
            "culture": clean_text(get_cell(row, idx_map, "Culture")) or "",
            "hydrolysis": clean_text(get_cell(row, idx_map, "Hydrolysis")) or "",
            "method_of_analysis": clean_text(get_cell(row, idx_map, "Method of analysis", "Method_of_analysis")) or "",
            "mz_ratio": clean_text(get_cell(row, idx_map, "M_Z ratio", "M_Z_Ratio", "M/Z ratio")) or "",
            "mass": clean_text(get_cell(row, idx_map, "Mass")) or "",
            "ic50_raw": clean_text(get_cell(row, idx_map, "IC50")) or "",
        })

    return rows, total_entries, warnings


def deduplicate_list_rows(rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
    dedup = OrderedDict()
    for row in rows:
        fmdb_id = clean_text(row.get("fmdb_id")) or ""
        if not fmdb_id:
            continue
        if fmdb_id not in dedup:
            dedup[fmdb_id] = row
    return list(dedup.values())


def parse_detail_pairs(plain_text: str) -> Dict[str, str]:
    fields = {
        "detail_title": "",
        "pubmed_id": "",
        "peptide_sequence": "",
        "length_of_peptide": "",
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
        "pubmed_id": [r"PubMed ID\s*:\s*(.+?)\s*(?:Title|Peptide Sequence|$)"],
        "peptide_sequence": [r"Peptide Sequence\s*:\s*(.+?)\s*(?:Length of peptide|PubMed ID|Food Matrix|$)"],
        "length_of_peptide": [r"Length of peptide\s*:\s*(.+?)\s*(?:Food Matrix|Protein Name|PubMed ID|$)"],
        "food_matrix": [r"Food Matrix\s*:\s*(.+?)\s*(?:Protein Name|pH|Temperature|$)"],
        "protein_name": [r"Protein Name\s*:\s*(.+?)\s*(?:pH|Temperature|Incubation Time|$)"],
        "ph": [r"pH\s*:\s*(.+?)\s*(?:Temperature|Incubation Time|Activity|$)"],
        "temperature": [r"Temperature\s*:\s*(.+?)\s*(?:Incubation Time|Activity|Experiment|$)"],
        "incubation_time": [r"Incubation Time\s*:\s*(.+?)\s*(?:Activity|Experiment|Model|$)"],
        "activity": [r"Activity\s*:\s*(.+?)\s*(?:Experiment|Model|Assay for Activity Measurement|$)"],
        "experiment": [r"Experiment\s*:\s*(.+?)\s*(?:Model|Assay for Activity Measurement|Starter Culture|$)"],
        "model": [r"Model\s*:\s*(.+?)\s*(?:Assay for Activity Measurement|Starter Culture|Hydrolysis|$)"],
        "assay_for_activity_measurement": [r"Assay for Activity Measurement\s*:\s*(.+?)\s*(?:Starter Culture|Hydrolysis|Method of Analysis|$)"],
        "starter_culture": [r"(?:Starter Culture|Culture)\s*:\s*(.+?)\s*(?:Hydrolysis|Method of Analysis|M/Z ratio|M_Z ratio|$)"],
        "hydrolysis": [r"Hydrolysis\s*:\s*(.+?)\s*(?:Method of Analysis|M/Z ratio|M_Z ratio|Mass|$)"],
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
    parser = parse_html_tables(html_text)
    plain_text = parser.plain_text
    pairs = parse_detail_pairs(plain_text)

    return {
        "fmdb_id": fmdb_id,
        "detail_url": detail_url,
        "detail_title": clean_text(pairs["detail_title"]) or "",
        "pubmed_id": clean_text(pairs["pubmed_id"]) or "",
        "peptide_sequence": clean_text(pairs["peptide_sequence"]) or "",
        "length_of_peptide": clean_text(pairs["length_of_peptide"]) or "",
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
    parser = argparse.ArgumentParser(description="FermFooDb ACE downloader v2（纯标准库版）")
    parser.add_argument("--project-root", type=str, default=None, help="PepDB 项目根目录")
    parser.add_argument("--start-url", type=str, default=START_URL, help="ACE 列表入口 URL")
    parser.add_argument("--timeout", type=int, default=20, help="请求超时秒数")
    parser.add_argument("--sleep-sec", type=float, default=0.3, help="请求间隔秒数")
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

    # 1) 抓 ACE 列表页（只抓这一页）
    result = fetch_url(args.start_url, timeout=args.timeout, sleep_sec=0.0)
    if not result["ok"]:
        raise RuntimeError(f"列表页抓取失败: {result['error_type']} | {result['error_message']}")

    page_html = result["text"]
    page_url = result["final_url"]

    list_html_path = list_html_dir / "list_page_001.html"
    list_txt_path = list_txt_dir / "list_page_001.txt"
    parser_obj = parse_html_tables(page_html)
    plain_text = parser_obj.plain_text
    save_text(list_html_path, page_html)
    save_text(list_txt_path, plain_text)

    list_rows_raw, total_entries_reported, warnings = parse_list_page(
        page_no=1,
        page_url=page_url,
        html_text=page_html,
    )
    list_rows = deduplicate_list_rows(list_rows_raw)
    write_csv(list_csv_path, list_rows, LIST_COLUMNS)

    # 2) 抓详情页
    detail_targets = [
        (clean_text(r.get("fmdb_id")) or "", clean_text(r.get("detail_url")) or "")
        for r in list_rows
        if clean_text(r.get("fmdb_id")) and clean_text(r.get("detail_url"))
    ]

    if args.max_detail is not None:
        detail_targets = detail_targets[: args.max_detail]

    detail_rows: List[Dict[str, object]] = []
    failed_detail_rows: List[Dict[str, object]] = []

    for fmdb_id, detail_url in detail_targets:
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
        detail_parser = parse_html_tables(detail_html)
        detail_text = detail_parser.plain_text

        html_path = detail_html_dir / f"{fmdb_id}.html"
        txt_path = detail_txt_dir / f"{fmdb_id}.txt"
        save_text(html_path, detail_html)
        save_text(txt_path, detail_text)

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
        "list_pages_fetched": 1,
        "list_rows_parsed_before_dedup": len(list_rows_raw),
        "list_rows_parsed_after_dedup": len(list_rows),
        "total_entries_reported_on_site": total_entries_reported,
        "warnings": warnings,
        "detail_targets": len(detail_targets),
        "detail_rows_parsed": len(detail_rows),
        "failed_detail_count": len(failed_detail_rows),
        "list_csv_path": str(list_csv_path),
        "detail_csv_path": str(detail_csv_path),
    }
    crawl_summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")

    print("=" * 80)
    print("FermFooDb ACE downloader v2 完成")
    print(f"项目根目录：{project_root}")
    print(f"列表页抓取数：1")
    print(f"列表记录数（去重前）：{len(list_rows_raw)}")
    print(f"列表记录数（去重后）：{len(list_rows)}")
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