# -*- coding: utf-8 -*-
"""
BIOPEP-UWM experimental（ACE inhibitor）原始下载脚本 v3
=========================================================

【脚本目标】
1. 抓取 BIOPEP-UWM experimental 列表页原始 HTML / TXT；
2. 从列表页中筛出 Name = "ACE inhibitor" 的条目；
3. 针对每个 ACE 条目，抓取更适合保存与查阅的报告页：peptidedatacard.php?zm_ID=...
4. 为列表页与详情页分别输出 CSV，方便本地 Excel 查阅；
5. 当前阶段只做“原始下载与归档”，不做清洗、标准化、单位换算、去重。

【这版相对前版的关键修正】
- 不再依赖脆弱的 HTML 表格结构解析；
- 列表页按“Peptide Data + 行文本模式”解析；
- 详情页改抓 peptidedatacard.php 报告页，因为该页把序列、活性、文献、交叉库引用集中在一个页面；
- 若首页解析不到任何记录，脚本会立即报错停止，避免继续空跑几百页；
- 翻页 URL 改为稳定候选方案，优先尝试不带后缀的 pageNum_result1，再尝试兼容旧站点的 bbb / aaa 形式；
- 控制台实时打印进度。

【技术约束】
- 仅使用 Python 标准库；
- 默认适配 Windows 路径；
- CSV 统一输出为 UTF-8-SIG，方便 Excel 打开。
"""

from __future__ import annotations

import csv
import json
import os
import re
import ssl
import sys
import time
import traceback
from dataclasses import dataclass, asdict
from html import unescape
from typing import Dict, List, Optional, Tuple
from urllib.parse import urljoin
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError


# ============================================================
# 一、前部配置区（路径、抓取行为、调试开关）
# ============================================================

# 是否优先使用下面手动配置的路径
PRIORITY_USE_MANUAL_PATHS = True

# PepDB 项目根目录
MANUAL_PROJECT_ROOT = r"E:\MYS\PepDB"

# 输出目录：原始下载结果将放在这里
MANUAL_OUTPUT_ROOT = os.path.join(
    MANUAL_PROJECT_ROOT,
    "DB",
    "raw",
    "ace",
    "databases",
    "BIOPEP_UWM",
    "experimental",
)

# 站点配置
BASE_URL = "https://biochemia.uwm.edu.pl/biopep/"
LIST_URL = urljoin(BASE_URL, "peptide_data.php")
DETAIL_REPORT_URL_TEMPLATE = urljoin(BASE_URL, "peptidedatacard.php?zm_ID={zm_id}")
TARGET_ACTIVITY_NAME = "ACE inhibitor"
TARGET_ACTIVITY_CODE = "ah"

# 网络配置
REQUEST_TIMEOUT = 15
REQUEST_RETRY_TIMES = 2
REQUEST_RETRY_SLEEP_SECONDS = 2.0
REQUEST_SLEEP_BETWEEN_PAGES = 0.6
REQUEST_SLEEP_BETWEEN_DETAILS = 0.3
USER_AGENT = (
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
    "AppleWebKit/537.36 (KHTML, like Gecko) "
    "Chrome/124.0 Safari/537.36"
)

# SSL 证书失败时，是否自动回退为“不校验证书继续下载”
ALLOW_INSECURE_SSL_FALLBACK = True

# 覆盖 / 跳过策略
OVERWRITE_EXISTING_FILES = False
SKIP_ALREADY_DOWNLOADED_DETAIL_HTML = True

# 调试开关：正式跑保持 None
DEBUG_MAX_LIST_PAGES: Optional[int] = None
DEBUG_MAX_DETAIL_RECORDS: Optional[int] = None

# 是否保存 HTML / TXT
SAVE_RAW_HTML = True
SAVE_LIST_PAGE_TXT = True
SAVE_DETAIL_TXT = True


# ============================================================
# 二、数据结构
# ============================================================

@dataclass
class ListRecord:
    list_page_number: int
    list_page_url: str
    list_row_index: int
    peptide_id: int
    name: str
    sequence: str
    chemical_mass: str
    monoisotopic_mass: str
    activity_value_raw: str
    measure_type_raw: str
    detail_report_url: str


@dataclass
class DetailRecord:
    peptide_id: int
    detail_report_url: str
    name: str
    sequence: str
    function_text: str
    number_of_residues: str
    activity_code: str
    activity_name: str
    chemical_mass: str
    monoisotopic_mass: str
    activity_measure_label: str
    activity_value_text: str
    bibliographic_authors: str
    bibliographic_title: str
    bibliographic_year: str
    bibliographic_source_type: str
    bibliographic_raw: str
    additional_information_raw: str
    smiles: str
    inchi: str
    inchikey: str
    database_reference_raw: str
    raw_text_path: str
    raw_html_path: str


# ============================================================
# 三、路径工具
# ============================================================


def get_output_root() -> str:
    if PRIORITY_USE_MANUAL_PATHS:
        return MANUAL_OUTPUT_ROOT
    return MANUAL_OUTPUT_ROOT



def build_output_paths(output_root: str) -> Dict[str, str]:
    return {
        "root": output_root,
        "raw_html_root": os.path.join(output_root, "raw_html"),
        "raw_html_list": os.path.join(output_root, "raw_html", "list_pages"),
        "raw_html_detail": os.path.join(output_root, "raw_html", "detail_pages"),
        "raw_txt_root": os.path.join(output_root, "raw_txt"),
        "raw_txt_list": os.path.join(output_root, "raw_txt", "list_pages"),
        "raw_txt_detail": os.path.join(output_root, "raw_txt", "detail_cards"),
        "csv_root": os.path.join(output_root, "raw_tables_csv"),
        "meta_root": os.path.join(output_root, "meta"),
        "log_file": os.path.join(output_root, "meta", "run_log.txt"),
        "summary_json": os.path.join(output_root, "meta", "crawl_summary.json"),
        "readme_txt": os.path.join(output_root, "meta", "README_抓取说明.txt"),
        "debug_page1_txt": os.path.join(output_root, "meta", "debug_page1_text.txt"),
        "debug_page1_html": os.path.join(output_root, "meta", "debug_page1_raw.html"),
        "debug_page2_txt": os.path.join(output_root, "meta", "debug_page2_text.txt"),
        "debug_page2_html": os.path.join(output_root, "meta", "debug_page2_raw.html"),
        "ace_list_csv": os.path.join(output_root, "raw_tables_csv", "biopep_uwm_experimental_ace_list.csv"),
        "ace_detail_csv": os.path.join(output_root, "raw_tables_csv", "biopep_uwm_experimental_ace_detail.csv"),
        "failed_detail_csv": os.path.join(output_root, "meta", "failed_detail_urls.csv"),
    }



def ensure_directories(paths: Dict[str, str]) -> None:
    for key in [
        "root",
        "raw_html_root",
        "raw_html_list",
        "raw_html_detail",
        "raw_txt_root",
        "raw_txt_list",
        "raw_txt_detail",
        "csv_root",
        "meta_root",
    ]:
        os.makedirs(paths[key], exist_ok=True)


# ============================================================
# 四、通用文本 / 文件工具
# ============================================================


def write_text(path: str, text: str, encoding: str = "utf-8-sig") -> None:
    with open(path, "w", encoding=encoding, newline="") as f:
        f.write(text)



def write_json(path: str, data: dict) -> None:
    with open(path, "w", encoding="utf-8-sig") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)



def write_csv(path: str, rows: List[Dict[str, object]], fieldnames: List[str]) -> None:
    with open(path, "w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})



def append_log(log_file: str, text: str) -> None:
    with open(log_file, "a", encoding="utf-8-sig") as f:
        f.write(text.rstrip() + "\n")



def log_and_print(log_file: str, text: str) -> None:
    print(text, flush=True)
    append_log(log_file, text)



def normalize_whitespace(text: str) -> str:
    text = text.replace("\xa0", " ")
    text = text.replace("\r", "\n")
    text = re.sub(r"[ \t]+", " ", text)
    text = re.sub(r"\n{3,}", "\n\n", text)
    return text.strip()



def slugify_filename(text: str) -> str:
    text = re.sub(r"[^0-9A-Za-z_\-\.]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("_")
    return text or "untitled"



def safe_int(text: str) -> Optional[int]:
    try:
        return int(text)
    except Exception:
        return None



def html_to_text(html: str) -> str:
    """
    把原始 HTML 尽量稳定地转成“可读文本”。

    说明：
    - 这里不是做严格 DOM 解析，而是为了后续用正则从文本里抽字段；
    - 对 BIOPEP 这种老站，很多页面结构很朴素，用文本模式反而比死扣表格标签更稳。
    """
    text = html
    text = re.sub(r"(?is)<script[^>]*>.*?</script>", " ", text)
    text = re.sub(r"(?is)<style[^>]*>.*?</style>", " ", text)
    text = re.sub(r"(?is)<!--.*?-->", " ", text)

    # 先给常见块级标签补换行，尽量保留视觉结构
    text = re.sub(r"(?is)<br\s*/?>", "\n", text)
    text = re.sub(r"(?is)</(p|div|tr|table|ul|ol|li|h\d)>", "\n", text)
    text = re.sub(r"(?is)</(td|th)>", " ", text)

    # 再去掉其他标签
    text = re.sub(r"(?is)<[^>]+>", " ", text)
    text = unescape(text)

    # BIOPEP 页面里常见的 IC_{50} / EC_{50} / IC 50 统一处理
    text = text.replace("IC_{50}", "IC50")
    text = text.replace("EC_{50}", "EC50")
    text = text.replace("IC 50", "IC50")
    text = text.replace("EC 50", "EC50")

    # 一些视觉上有意义的模式，主动补换行，方便后续 split
    text = re.sub(r"\bPeptide Data\b", "\nPeptide Data ", text)
    text = re.sub(r"\bFunction:\b", "\nFunction:\n", text)
    text = re.sub(r"\bBibliographic data:\b", "\nBibliographic data:\n", text)
    text = re.sub(r"\bAdditional information:\b", "\nAdditional information:\n", text)
    text = re.sub(r"\bDatabase reference:\b", "\nDatabase reference:\n", text)
    text = re.sub(r"\bID\s+(\d+)\b", r"\nID \1", text)

    return normalize_whitespace(text)


# ============================================================
# 五、网络请求工具
# ============================================================


def _decode_html_bytes(raw: bytes, resp) -> str:
    encodings = []
    try:
        charset = resp.headers.get_content_charset() or ""
        if charset:
            encodings.append(charset)
    except Exception:
        pass
    encodings.extend(["utf-8", "utf-8-sig", "latin1"])

    for enc in encodings:
        try:
            return raw.decode(enc)
        except Exception:
            continue
    return raw.decode("utf-8", errors="replace")



def _is_ssl_cert_verify_error(exc: BaseException) -> bool:
    current = exc
    seen = set()
    while current is not None and id(current) not in seen:
        seen.add(id(current))
        if isinstance(current, ssl.SSLCertVerificationError):
            return True
        msg = str(current).lower()
        if "certificate verify failed" in msg or "ssl: cert" in msg:
            return True
        current = getattr(current, "reason", None) or getattr(current, "__cause__", None)
    return False



def fetch_html(url: str, timeout: int = REQUEST_TIMEOUT) -> str:
    headers = {"User-Agent": USER_AGENT}
    req = Request(url, headers=headers)

    last_exc: Optional[BaseException] = None
    for attempt in range(1, REQUEST_RETRY_TIMES + 1):
        try:
            try:
                with urlopen(req, timeout=timeout) as resp:
                    raw = resp.read()
                    return _decode_html_bytes(raw, resp)
            except Exception as exc:
                if ALLOW_INSECURE_SSL_FALLBACK and _is_ssl_cert_verify_error(exc):
                    ctx = ssl.create_default_context()
                    ctx.check_hostname = False
                    ctx.verify_mode = ssl.CERT_NONE
                    with urlopen(req, timeout=timeout, context=ctx) as resp:
                        raw = resp.read()
                        return _decode_html_bytes(raw, resp)
                raise
        except Exception as exc:
            last_exc = exc
            if attempt < REQUEST_RETRY_TIMES:
                time.sleep(REQUEST_RETRY_SLEEP_SECONDS)
            else:
                break

    raise RuntimeError(f"下载失败：{url}\n最后一次错误：{last_exc}") from last_exc


# ============================================================
# 六、列表页 URL 与解析逻辑
# ============================================================


def build_list_page_candidate_urls(page_number: int) -> List[str]:
    """
    构造某一页的候选 URL。

    说明：
    - BIOPEP 列表第一页直接是 peptide_data.php
    - 其余分页常见写法是 pageNum_result1={page_index_from_0}
    - 旧站点上偶尔能看到 pageNum_result1=3bbb / 372aaa 之类形式；
      这些后缀本质上会被 PHP 当作前面的整数部分来处理。
    - 为了兼容，脚本把这几种形式都当候选。
    """
    if page_number == 1:
        return [LIST_URL]

    offset = page_number - 1
    return [
        f"{LIST_URL}?pageNum_result1={offset}",
        f"{LIST_URL}?pageNum_result1={offset}bbb",
        f"{LIST_URL}?pageNum_result1={offset}aaa",
    ]



def parse_total_pages_from_text(text: str) -> Optional[int]:
    matches = re.findall(r"Page\s+(\d+)\s*/\s*(\d+)", text)
    if not matches:
        return None
    # 一般最后一个匹配最接近页面底部的真实页码提示
    current, total = matches[-1]
    return safe_int(total)



def try_fetch_list_page(page_number: int) -> Tuple[str, str]:
    """
    尝试下载指定页的列表页，并返回：
    - 实际使用的 URL
    - HTML 内容

    选择逻辑：
    - 若页面文本里含有 Page X / Y，且 X 与目标页码一致，则优先认为成功；
    - 否则退而求其次，只要拿到一个非空 HTML，也返回给上层，由上层决定是否可解析。
    """
    candidates = build_list_page_candidate_urls(page_number)
    fallback: Optional[Tuple[str, str]] = None

    for url in candidates:
        html = fetch_html(url)
        if html and fallback is None:
            fallback = (url, html)

        text = html_to_text(html)
        m = re.search(r"Page\s+(\d+)\s*/\s*(\d+)", text)
        if m:
            current_page = safe_int(m.group(1))
            if current_page == page_number:
                return url, html

    if fallback is not None:
        return fallback

    raise RuntimeError(f"列表页下载失败：page={page_number}")



def strip_html_tags(html_fragment: str) -> str:
    """把一个 HTML 片段转成简洁文本。"""
    frag = html_fragment
    frag = re.sub(r"(?is)<br\s*/?>", " ", frag)
    frag = re.sub(r"(?is)</?(sub|sup)[^>]*>", "", frag)
    frag = re.sub(r"(?is)<[^>]+>", " ", frag)
    frag = unescape(frag)
    frag = frag.replace("IC 50", "IC50").replace("EC 50", "EC50")
    return normalize_whitespace(frag)



def parse_activity_cell_text(cell_text: str) -> Tuple[str, str]:
    """把最后一列类似 `29.00 IC50` / `4.40 EC50` 拆成数值与类型。"""
    text = normalize_whitespace(cell_text)
    text = text.replace("IC 50", "IC50").replace("EC 50", "EC50")
    m = re.search(r"(-?\d+(?:\.\d+)?)\s*((?:IC|EC)50)\b", text, re.IGNORECASE)
    if not m:
        return text, ""
    return m.group(1), m.group(2).upper()



def parse_list_rows_from_html(html: str) -> List[Dict[str, str]]:
    """
    直接从列表页 HTML 解析数据行。

    这一版不再依赖 html_to_text 之后的松散文本结构，而是直接按：
    - `<tr class="info"> ... </tr>` 抓数据行
    - 每行再按 `<td>...</td>` 抓 7 个单元格

    这样能稳定处理 BIOPEP 首页这类 `Peptide / Data` 被拆行的问题。
    """
    records: List[Dict[str, str]] = []

    row_blocks = re.findall(r"(?is)<tr\s+class=['\"]info['\"][^>]*>(.*?)</tr>", html)
    for row_html in row_blocks:
        cells = re.findall(r"(?is)<td\b[^>]*>(.*?)</td>", row_html)
        if len(cells) < 7:
            continue

        link_cell = cells[0]
        peptide_id_text = strip_html_tags(cells[1])
        name = strip_html_tags(cells[2])
        sequence = strip_html_tags(cells[3])
        chemical_mass = strip_html_tags(cells[4])
        monoisotopic_mass = strip_html_tags(cells[5])
        activity_cell_text = strip_html_tags(cells[6])
        activity_value_raw, measure_type_raw = parse_activity_cell_text(activity_cell_text)

        peptide_id = safe_int(peptide_id_text)
        if peptide_id is None:
            continue

        href_match = re.search(r"href=['\"']([^'\"']*zm_ID=(\d+)[^'\"']*)['\"']", link_cell, re.IGNORECASE)
        detail_list_url = ""
        if href_match:
            detail_list_url = urljoin(BASE_URL, href_match.group(1))

        records.append({
            "peptide_id": str(peptide_id),
            "name": name,
            "sequence": sequence,
            "chemical_mass": chemical_mass,
            "monoisotopic_mass": monoisotopic_mass,
            "activity_value_raw": activity_value_raw,
            "measure_type_raw": measure_type_raw,
            "detail_list_url": detail_list_url,
        })

    return records



def parse_list_rows_from_text(text: str) -> List[Dict[str, str]]:
    """
    文本模式兜底解析。

    注意：BIOPEP 首页会把 `Peptide` 和 `Data` 在纯文本里拆开，
    因此这条路径只作为 HTML 解析失败时的备用，不再作为主解析逻辑。
    """
    records: List[Dict[str, str]] = []

    text = re.sub(r"Peptide\s+Data", "Peptide Data", text, flags=re.IGNORECASE)
    parts = re.split(r"\bPeptide Data\b", text)
    if len(parts) <= 1:
        return records

    row_pattern = re.compile(
        r"^\s*(?P<peptide_id>\d+)\s+"
        r"(?P<name>.+?)\s+"
        r"(?P<sequence>[A-Za-z0-9~\-\*\.]+)\s+"
        r"(?P<chemical_mass>\d+(?:\.\d+)?)\s+"
        r"(?P<monoisotopic_mass>\d+(?:\.\d+)?)\s+"
        r"(?P<activity_value_raw>\d+(?:\.\d+)?)\s+"
        r"(?P<measure_type_raw>(?:IC|EC)50)\b",
        re.IGNORECASE,
    )

    for seg in parts[1:]:
        seg = normalize_whitespace(seg)
        if not seg:
            continue
        seg = re.split(r"\b(?:First\b|Search:|Page\s+\d+\s*/\s*\d+)\b", seg, maxsplit=1)[0].strip()
        if not seg:
            continue
        m = row_pattern.search(seg)
        if not m:
            continue
        rec = {k: v.strip() for k, v in m.groupdict().items()}
        rec["detail_list_url"] = ""
        records.append(rec)

    return records


# ============================================================
# 七、详情页解析逻辑（基于 peptidedatacard 报告页）
# ============================================================


def extract_first(pattern: str, text: str, flags: int = re.IGNORECASE | re.DOTALL) -> str:
    m = re.search(pattern, text, flags)
    if not m:
        return ""
    return normalize_whitespace(m.group(1))



def parse_detail_report_text(text: str, peptide_id: int, detail_report_url: str, raw_text_path: str, raw_html_path: str) -> DetailRecord:
    # 大块 section
    function_text = extract_first(r"Function:\s*(.*?)\s*Number of residues", text)
    number_of_residues = extract_first(r"Number of residues\s*(.*?)\s*Activity code", text)
    activity_code = extract_first(r"Activity code\s*(.*?)\s*Activity\s*:", text)
    activity_name = extract_first(r"Activity\s*:\s*(.*?)\s*Chemical mass", text)
    chemical_mass = extract_first(r"Chemical mass\s*(\d+(?:\.\d+)?)", text)
    monoisotopic_mass = extract_first(r"Monoisotopic mass\s*(\d+(?:\.\d+)?)", text)

    # 名称与序列
    name = extract_first(r"ID\s+\d+\s+Name\s+(.*?)\s+sequence", text)
    sequence = extract_first(r"sequence\s+(.*?)\s+\*\s*\*\s*\*", text)

    # 活性数值区块：可能是 IC50 / EC50 等
    measure_match = re.search(
        r"\b((?:IC|EC)50)\s*:\s*(.*?)\s*\*\s*\*\s*\*",
        text,
        re.IGNORECASE | re.DOTALL,
    )
    if measure_match:
        activity_measure_label = normalize_whitespace(measure_match.group(1)).upper()
        activity_value_text = normalize_whitespace(measure_match.group(2))
    else:
        activity_measure_label = ""
        activity_value_text = ""

    bibliographic_raw = extract_first(r"Bibliographic data:\s*(.*?)\s*\*\s*\*\s*\*", text)
    bibliographic_authors = extract_first(r"Authors\s*(.*?)\s*Title", bibliographic_raw)
    bibliographic_title = extract_first(r"Title\s*(.*?)\s*Year\s+Source", bibliographic_raw)
    bibliographic_year = extract_first(r"Year\s+Source\s*(\d{4})", bibliographic_raw)
    bibliographic_source_type = extract_first(r"Year\s+Source\s*\d{4}\s+(.+)$", bibliographic_raw)

    additional_information_raw = extract_first(r"Additional information:\s*(.*?)\s*\*\s*\*\s*\*", text)
    smiles = extract_first(r"SMILES:\s*(.*?)\s*InChI=", additional_information_raw)
    inchi = extract_first(r"(InChI=[^\n]+)", additional_information_raw)
    inchikey = extract_first(r"InChIKey:\s*(.*?)\s*(?:Information concerning|$)", additional_information_raw)

    database_reference_raw = extract_first(r"Database reference:\s*(.*)$", text)

    return DetailRecord(
        peptide_id=peptide_id,
        detail_report_url=detail_report_url,
        name=name,
        sequence=sequence,
        function_text=function_text,
        number_of_residues=number_of_residues,
        activity_code=activity_code,
        activity_name=activity_name,
        chemical_mass=chemical_mass,
        monoisotopic_mass=monoisotopic_mass,
        activity_measure_label=activity_measure_label,
        activity_value_text=activity_value_text,
        bibliographic_authors=bibliographic_authors,
        bibliographic_title=bibliographic_title,
        bibliographic_year=bibliographic_year,
        bibliographic_source_type=bibliographic_source_type,
        bibliographic_raw=bibliographic_raw,
        additional_information_raw=additional_information_raw,
        smiles=smiles,
        inchi=inchi,
        inchikey=inchikey,
        database_reference_raw=database_reference_raw,
        raw_text_path=raw_text_path,
        raw_html_path=raw_html_path,
    )


# ============================================================
# 八、抓取主流程：列表页
# ============================================================


def save_debug_snapshot_if_needed(paths: Dict[str, str], page_number: int, html: str, text: str) -> None:
    if page_number == 1:
        write_text(paths["debug_page1_html"], html, encoding="utf-8")
        write_text(paths["debug_page1_txt"], text)
    elif page_number == 2:
        write_text(paths["debug_page2_html"], html, encoding="utf-8")
        write_text(paths["debug_page2_txt"], text)



def crawl_list_pages(paths: Dict[str, str]) -> Tuple[List[ListRecord], Dict[str, int]]:
    log_file = paths["log_file"]
    all_records: List[ListRecord] = []
    total_pages_seen = 0
    total_rows_parsed = 0
    target_rows_found = 0

    # 先抓第 1 页，确定总页数，并验证解析是否正常
    page_number = 1
    page_url, html = try_fetch_list_page(page_number)
    text = html_to_text(html)
    rows = parse_list_rows_from_html(html) or parse_list_rows_from_text(text)
    save_debug_snapshot_if_needed(paths, page_number, html, text)

    if SAVE_RAW_HTML:
        write_text(os.path.join(paths["raw_html_list"], f"list_page_{page_number:04d}.html"), html, encoding="utf-8")
    if SAVE_LIST_PAGE_TXT:
        write_text(os.path.join(paths["raw_txt_list"], f"list_page_{page_number:04d}.txt"), text)

    if not rows:
        raise RuntimeError(
            "首页列表页解析结果为 0 行，脚本已停止。\n"
            "请查看 meta/debug_page1_text.txt 与 meta/debug_page1_raw.html 排查页面结构。"
        )

    total_pages = parse_total_pages_from_text(text) or 1
    total_pages_seen = max(total_pages_seen, total_pages)

    ace_rows_this_page = 0
    for row_index, row in enumerate(rows, start=1):
        total_rows_parsed += 1
        if normalize_whitespace(row["name"]) != TARGET_ACTIVITY_NAME:
            continue
        ace_rows_this_page += 1
        target_rows_found += 1
        peptide_id = safe_int(row["peptide_id"]) or 0
        detail_report_url = DETAIL_REPORT_URL_TEMPLATE.format(zm_id=peptide_id)
        all_records.append(
            ListRecord(
                list_page_number=page_number,
                list_page_url=page_url,
                list_row_index=row_index,
                peptide_id=peptide_id,
                name=row["name"],
                sequence=row["sequence"],
                chemical_mass=row["chemical_mass"],
                monoisotopic_mass=row["monoisotopic_mass"],
                activity_value_raw=row["activity_value_raw"],
                measure_type_raw=row["measure_type_raw"],
                detail_report_url=detail_report_url,
            )
        )

    log_and_print(log_file, f"[列表页] page={page_number} rows={len(rows)} target={ace_rows_this_page} url={page_url}")

    max_pages_to_run = total_pages
    if DEBUG_MAX_LIST_PAGES is not None:
        max_pages_to_run = min(max_pages_to_run, DEBUG_MAX_LIST_PAGES)

    # 再抓后续页
    for page_number in range(2, max_pages_to_run + 1):
        time.sleep(REQUEST_SLEEP_BETWEEN_PAGES)
        page_url, html = try_fetch_list_page(page_number)
        text = html_to_text(html)
        rows = parse_list_rows_from_html(html) or parse_list_rows_from_text(text)
        save_debug_snapshot_if_needed(paths, page_number, html, text)

        if SAVE_RAW_HTML:
            write_text(os.path.join(paths["raw_html_list"], f"list_page_{page_number:04d}.html"), html, encoding="utf-8")
        if SAVE_LIST_PAGE_TXT:
            write_text(os.path.join(paths["raw_txt_list"], f"list_page_{page_number:04d}.txt"), text)

        # 这里不再立即报错；若个别页失败，保留调试材料，但继续跑
        ace_rows_this_page = 0
        for row_index, row in enumerate(rows, start=1):
            total_rows_parsed += 1
            if normalize_whitespace(row["name"]) != TARGET_ACTIVITY_NAME:
                continue
            ace_rows_this_page += 1
            target_rows_found += 1
            peptide_id = safe_int(row["peptide_id"]) or 0
            detail_report_url = DETAIL_REPORT_URL_TEMPLATE.format(zm_id=peptide_id)
            all_records.append(
                ListRecord(
                    list_page_number=page_number,
                    list_page_url=page_url,
                    list_row_index=row_index,
                    peptide_id=peptide_id,
                    name=row["name"],
                    sequence=row["sequence"],
                    chemical_mass=row["chemical_mass"],
                    monoisotopic_mass=row["monoisotopic_mass"],
                    activity_value_raw=row["activity_value_raw"],
                    measure_type_raw=row["measure_type_raw"],
                    detail_report_url=detail_report_url,
                )
            )

        log_and_print(log_file, f"[列表页] page={page_number} rows={len(rows)} target={ace_rows_this_page} url={page_url}")

    return all_records, {
        "list_pages_crawled": max_pages_to_run,
        "list_total_pages_seen": total_pages_seen,
        "list_total_rows_parsed": total_rows_parsed,
        "target_rows_found": target_rows_found,
    }


# ============================================================
# 九、抓取主流程：详情页（报告页）
# ============================================================


def crawl_detail_pages(paths: Dict[str, str], list_records: List[ListRecord]) -> Tuple[List[DetailRecord], List[Dict[str, str]]]:
    log_file = paths["log_file"]
    detail_records: List[DetailRecord] = []
    failed_rows: List[Dict[str, str]] = []

    max_details_to_run = len(list_records)
    if DEBUG_MAX_DETAIL_RECORDS is not None:
        max_details_to_run = min(max_details_to_run, DEBUG_MAX_DETAIL_RECORDS)

    for idx, rec in enumerate(list_records[:max_details_to_run], start=1):
        peptide_id = rec.peptide_id
        detail_report_url = rec.detail_report_url

        html_filename = f"detail_report_{peptide_id}.html"
        txt_filename = f"detail_report_{peptide_id}.txt"
        html_path = os.path.join(paths["raw_html_detail"], html_filename)
        txt_path = os.path.join(paths["raw_txt_detail"], txt_filename)

        try:
            if SKIP_ALREADY_DOWNLOADED_DETAIL_HTML and os.path.exists(html_path) and not OVERWRITE_EXISTING_FILES:
                with open(html_path, "r", encoding="utf-8", errors="replace") as f:
                    html = f.read()
            else:
                html = fetch_html(detail_report_url)
                if SAVE_RAW_HTML and (OVERWRITE_EXISTING_FILES or not os.path.exists(html_path)):
                    write_text(html_path, html, encoding="utf-8")
                time.sleep(REQUEST_SLEEP_BETWEEN_DETAILS)

            text = html_to_text(html)
            if SAVE_DETAIL_TXT and (OVERWRITE_EXISTING_FILES or not os.path.exists(txt_path)):
                write_text(txt_path, text)

            detail_record = parse_detail_report_text(
                text=text,
                peptide_id=peptide_id,
                detail_report_url=detail_report_url,
                raw_text_path=txt_path,
                raw_html_path=html_path,
            )
            detail_records.append(detail_record)
            log_and_print(log_file, f"[详情页] {idx}/{max_details_to_run} success id={peptide_id} url={detail_report_url}")
        except Exception as exc:
            failed_rows.append({
                "peptide_id": str(peptide_id),
                "detail_report_url": detail_report_url,
                "error": str(exc),
            })
            log_and_print(log_file, f"[详情页] {idx}/{max_details_to_run} FAILED id={peptide_id} url={detail_report_url} error={exc}")

    return detail_records, failed_rows


# ============================================================
# 十、README 与汇总输出
# ============================================================


def build_readme_text() -> str:
    return f"""BIOPEP-UWM experimental ACE inhibitor 原始抓取说明
===========================================

一、抓取对象
- 站点：BIOPEP-UWM experimental
- 列表页入口：{LIST_URL}
- 目标活动名称：{TARGET_ACTIVITY_NAME}
- 目标活动代码：{TARGET_ACTIVITY_CODE}
- 详情页抓取对象：peptidedatacard.php?zm_ID=...

二、目录说明
1. raw_html/list_pages
   - 每个列表页的原始 HTML
2. raw_html/detail_pages
   - 每个详情报告页的原始 HTML
3. raw_txt/list_pages
   - 每个列表页转出的便于人工阅读的 TXT
4. raw_txt/detail_cards
   - 每个详情报告页转出的便于人工阅读的 TXT
5. raw_tables_csv
   - biopep_uwm_experimental_ace_list.csv：ACE 条目列表索引表
   - biopep_uwm_experimental_ace_detail.csv：详情报告页字段汇总表
6. meta
   - run_log.txt：运行日志
   - crawl_summary.json：抓取摘要
   - failed_detail_urls.csv：失败详情页记录
   - debug_page1_raw.html / debug_page1_text.txt：首页调试快照

三、当前不做的事情
- 不做数据清洗
- 不做标准化
- 不做单位换算
- 不做跨库去重
- 不做与 AHTPDB 合并

四、推荐查看顺序
1. 先看 raw_tables_csv/biopep_uwm_experimental_ace_detail.csv
2. 再看 raw_tables_csv/biopep_uwm_experimental_ace_list.csv
3. 若发现字段缺失，再看 raw_txt/detail_cards
4. 若仍有疑问，再看 raw_html 原文
"""


# ============================================================
# 十一、主函数
# ============================================================


def main() -> int:
    output_root = get_output_root()
    paths = build_output_paths(output_root)
    ensure_directories(paths)

    # 每次运行前，重置日志文件，避免混淆旧记录
    write_text(paths["log_file"], "", encoding="utf-8-sig")
    write_text(paths["readme_txt"], build_readme_text())

    started_at = time.strftime("%Y-%m-%d %H:%M:%S")

    header = [
        "=" * 80,
        f"开始抓取 BIOPEP-UWM experimental ACE inhibitor 原始数据（v4）",
        f"输出目录：{output_root}",
        f"列表页入口：{LIST_URL}",
        f"目标活动名称：{TARGET_ACTIVITY_NAME}",
        f"目标活动代码：{TARGET_ACTIVITY_CODE}",
        "=" * 80,
    ]
    for line in header:
        log_and_print(paths["log_file"], line)

    try:
        list_records, list_summary = crawl_list_pages(paths)
        list_rows_for_csv = [asdict(x) for x in list_records]
        write_csv(
            paths["ace_list_csv"],
            list_rows_for_csv,
            fieldnames=[
                "list_page_number",
                "list_page_url",
                "list_row_index",
                "peptide_id",
                "name",
                "sequence",
                "chemical_mass",
                "monoisotopic_mass",
                "activity_value_raw",
                "measure_type_raw",
                "detail_report_url",
            ],
        )

        detail_records, failed_rows = crawl_detail_pages(paths, list_records)
        detail_rows_for_csv = [asdict(x) for x in detail_records]
        write_csv(
            paths["ace_detail_csv"],
            detail_rows_for_csv,
            fieldnames=[
                "peptide_id",
                "detail_report_url",
                "name",
                "sequence",
                "function_text",
                "number_of_residues",
                "activity_code",
                "activity_name",
                "chemical_mass",
                "monoisotopic_mass",
                "activity_measure_label",
                "activity_value_text",
                "bibliographic_authors",
                "bibliographic_title",
                "bibliographic_year",
                "bibliographic_source_type",
                "bibliographic_raw",
                "additional_information_raw",
                "smiles",
                "inchi",
                "inchikey",
                "database_reference_raw",
                "raw_text_path",
                "raw_html_path",
            ],
        )

        if failed_rows:
            write_csv(
                paths["failed_detail_csv"],
                failed_rows,
                fieldnames=["peptide_id", "detail_report_url", "error"],
            )

        ended_at = time.strftime("%Y-%m-%d %H:%M:%S")
        summary = {
            "started_at": started_at,
            "ended_at": ended_at,
            "output_root": output_root,
            "target_activity_name": TARGET_ACTIVITY_NAME,
            "target_activity_code": TARGET_ACTIVITY_CODE,
            **list_summary,
            "detail_success_count": len(detail_records),
            "detail_failed_count": len(failed_rows),
            "list_csv": paths["ace_list_csv"],
            "detail_csv": paths["ace_detail_csv"],
        }
        write_json(paths["summary_json"], summary)

        log_and_print(paths["log_file"], "抓取完成。")
        log_and_print(paths["log_file"], json.dumps(summary, ensure_ascii=False, indent=2))
        print("\n可直接查看的两个核心 CSV：", flush=True)
        print(paths["ace_list_csv"], flush=True)
        print(paths["ace_detail_csv"], flush=True)
        return 0

    except KeyboardInterrupt:
        log_and_print(paths["log_file"], "用户中断运行。")
        return 1
    except Exception as exc:
        ended_at = time.strftime("%Y-%m-%d %H:%M:%S")
        append_log(paths["log_file"], "程序运行失败： " + str(exc))
        append_log(paths["log_file"], traceback.format_exc())
        print(f"程序运行失败： {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
