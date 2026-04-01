# -*- coding: utf-8 -*-
"""
BIOPEP-UWM experimental 库 ACE inhibitor 原始下载脚本
---------------------------------------------------
用途：
1. 从 BIOPEP-UWM experimental 主列表逐页抓取原始页面；
2. 只筛出 Name 为 "ACE inhibitor" 的条目；
3. 进入每条详情页，下载原始 HTML，并提取为便于查阅的 TXT；
4. 额外生成两张 CSV：
   - 一张是列表页筛出的 ACE 条目索引表
   - 一张是详情页解析出的字段汇总表

重要说明：
- 这个脚本当前目标只是“把数据下载下来并整理成便于查阅的原始表”，
  不做后续清洗、标准化、单位换算、去重合并等工作。
- 为了兼容你当前本机环境，脚本只使用 Python 标准库，不依赖 pandas / bs4 / requests。
- 脚本默认把数据放到：
  E:/MYS/PepDB/DB/raw/ace/databases/BIOPEP_UWM/experimental
- 所有 CSV 统一用 UTF-8-SIG 编码，方便直接用 Excel 打开。

建议运行环境：
- Python 3.9+
- Windows / macOS / Linux 均可，但默认路径按你的 Windows 项目目录写死在前部配置里。
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
from dataclasses import dataclass
from html.parser import HTMLParser
from typing import Dict, List, Optional, Tuple
from urllib.parse import urljoin
from urllib.request import Request, urlopen
from urllib.error import URLError, HTTPError


# ============================================================
# 一、前部路径与行为配置（便于后续改路径、改活动类型、改调试行为）
# ============================================================

# 是否优先使用手动路径
PRIORITY_USE_MANUAL_PATHS = True

# 你的 PepDB 项目根目录（默认就是你当前长期项目）
MANUAL_PROJECT_ROOT = r"E:\MYS\PepDB"

# 默认输出根目录：BIOPEP-UWM experimental 原始数据目录
MANUAL_OUTPUT_ROOT = os.path.join(
    MANUAL_PROJECT_ROOT,
    "DB",
    "raw",
    "ace",
    "databases",
    "BIOPEP_UWM",
    "experimental",
)

# 目标站点配置
BASE_URL = "https://biochemia.uwm.edu.pl/biopep/"
LIST_URL = urljoin(BASE_URL, "peptide_data.php")
TARGET_ACTIVITY_NAME = "ACE inhibitor"
TARGET_ACTIVITY_CODE = "ah"  # BIOPEP-UWM 活性字典中 ACE inhibitor 的代码

# 网络行为配置
REQUEST_TIMEOUT = 30
REQUEST_RETRY_TIMES = 3
REQUEST_RETRY_SLEEP_SECONDS = 3
REQUEST_SLEEP_BETWEEN_PAGES = 0.8
REQUEST_SLEEP_BETWEEN_DETAILS = 0.6

# SSL 相关配置：
# 1. 默认先按正常 HTTPS 证书校验访问；
# 2. 若目标站点证书链在本机环境下无法通过校验，且这里允许回退，
#    则脚本会自动切换到“不校验证书”的模式继续抓取。
#
# 说明：
# - 这是为了解决某些站点证书链较旧、或本机 Python/Anaconda 证书根链不完整时的下载失败问题；
# - 当前目标只是原始数据采集，不涉及登录态、表单提交或敏感信息，因此这种回退在这里可接受；
# - 若你后续修复了本机证书环境，也可以把下面开关改回 False。
ALLOW_INSECURE_SSL_FALLBACK = True
USER_AGENT = (
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
    "AppleWebKit/537.36 (KHTML, like Gecko) "
    "Chrome/124.0 Safari/537.36"
)

# 文件覆盖 / 断点续跑配置
OVERWRITE_EXISTING_FILES = False
SKIP_ALREADY_DOWNLOADED_DETAIL_HTML = True

# 调试配置：正式跑请保持 None
DEBUG_MAX_LIST_PAGES: Optional[int] = None
DEBUG_MAX_DETAIL_RECORDS: Optional[int] = None

# 是否保存原始 HTML
SAVE_RAW_HTML = True

# 是否为每条详情额外保存 TXT 文本版
SAVE_DETAIL_TXT = True

# 是否为每个列表页额外保存 TXT 文本版
SAVE_LIST_PAGE_TXT = True


# ============================================================
# 二、路径相关工具
# ============================================================

def get_output_root() -> str:
    """返回输出根目录。当前脚本默认使用手动路径。"""
    if PRIORITY_USE_MANUAL_PATHS:
        return MANUAL_OUTPUT_ROOT
    return MANUAL_OUTPUT_ROOT


def build_output_paths(output_root: str) -> Dict[str, str]:
    """统一管理输出目录，避免路径散落在代码里。"""
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
        "ace_list_csv": os.path.join(output_root, "raw_tables_csv", "biopep_uwm_experimental_ace_list.csv"),
        "ace_detail_csv": os.path.join(output_root, "raw_tables_csv", "biopep_uwm_experimental_ace_detail.csv"),
    }


def ensure_directories(paths: Dict[str, str]) -> None:
    """创建所有需要的目录。"""
    dir_keys = [
        "root",
        "raw_html_root",
        "raw_html_list",
        "raw_html_detail",
        "raw_txt_root",
        "raw_txt_list",
        "raw_txt_detail",
        "csv_root",
        "meta_root",
    ]
    for key in dir_keys:
        os.makedirs(paths[key], exist_ok=True)


# ============================================================
# 三、通用工具函数
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
            clean_row = {k: row.get(k, "") for k in fieldnames}
            writer.writerow(clean_row)


def append_log(log_file: str, text: str) -> None:
    with open(log_file, "a", encoding="utf-8-sig") as f:
        f.write(text.rstrip() + "\n")


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


def safe_float(text: str) -> Optional[float]:
    try:
        return float(text)
    except Exception:
        return None


def _decode_html_bytes(raw: bytes, resp) -> str:
    """统一处理网页字节流解码，尽量兼容站点编码声明不规范的情况。"""
    content_type = ""
    try:
        content_type = resp.headers.get_content_charset() or ""
    except Exception:
        content_type = ""

    candidate_encodings = []
    if content_type:
        candidate_encodings.append(content_type)
    candidate_encodings.extend(["utf-8", "utf-8-sig", "latin1"])

    for enc in candidate_encodings:
        try:
            return raw.decode(enc)
        except Exception:
            continue
    return raw.decode("utf-8", errors="replace")


def _is_ssl_cert_verify_error(exc: BaseException) -> bool:
    """尽量稳妥地判断是否属于证书校验失败。"""
    current = exc
    seen = set()
    while current is not None and id(current) not in seen:
        seen.add(id(current))
        if isinstance(current, ssl.SSLCertVerificationError):
            return True
        message = str(current).lower()
        if "certificate verify failed" in message or "ssl: cert" in message:
            return True
        current = getattr(current, "reason", None) or getattr(current, "__cause__", None)
    return False


def _urlopen_with_optional_ssl_fallback(req: Request, timeout: int):
    """
    先用默认 SSL 校验访问；若命中证书校验失败且允许回退，则切换到不校验证书模式。
    返回已打开的响应对象。
    """
    try:
        return urlopen(req, timeout=timeout)
    except Exception as exc:
        if not ALLOW_INSECURE_SSL_FALLBACK:
            raise
        if not _is_ssl_cert_verify_error(exc):
            raise

        insecure_context = ssl._create_unverified_context()
        return urlopen(req, timeout=timeout, context=insecure_context)


def fetch_html(url: str, timeout: int = REQUEST_TIMEOUT) -> str:
    """
    下载网页 HTML。
    只用标准库 urllib，避免额外依赖。

    额外兼容：
    - 若本机 SSL 证书环境与目标站点证书链不兼容，
      可按前部配置自动回退到“不校验证书”模式继续抓取。
    """
    last_error = None
    for attempt in range(1, REQUEST_RETRY_TIMES + 1):
        try:
            req = Request(url, headers={"User-Agent": USER_AGENT})
            with _urlopen_with_optional_ssl_fallback(req, timeout=timeout) as resp:
                raw = resp.read()
                return _decode_html_bytes(raw, resp)
        except (URLError, HTTPError, TimeoutError, OSError, ssl.SSLError) as exc:
            last_error = exc
            if attempt < REQUEST_RETRY_TIMES:
                time.sleep(REQUEST_RETRY_SLEEP_SECONDS)
            else:
                raise RuntimeError(f"下载失败：{url}\n最后一次错误：{exc}") from exc
    raise RuntimeError(f"下载失败：{url}\n最后一次错误：{last_error}")


# ============================================================
# 四、HTML 解析器：列表页表格 + 链接 + 通用文本提取
# ============================================================

class TableAndLinkParser(HTMLParser):
    """
    一个尽量轻量、但足够实用的 HTML 解析器。

    作用：
    1. 解析列表页中的表格行 / 单元格；
    2. 记录所有链接，便于找到“Next”“Last”“Peptide Data”等链接；
    3. 顺便生成一个带换行的纯文本版本，方便保存 TXT。
    """

    BLOCK_TAGS = {"p", "div", "br", "tr", "table", "hr", "li", "ul", "ol", "h1", "h2", "h3", "h4", "h5", "h6"}

    def __init__(self) -> None:
        super().__init__()
        self.rows: List[List[Dict[str, object]]] = []
        self.links: List[Dict[str, str]] = []

        self._in_tr = False
        self._in_td = False
        self._current_row: List[Dict[str, object]] = []
        self._current_cell_text_parts: List[str] = []
        self._current_cell_links: List[str] = []

        self._current_anchor_href: Optional[str] = None
        self._current_anchor_text_parts: List[str] = []

        self._text_parts: List[str] = []

    def handle_starttag(self, tag: str, attrs: List[Tuple[str, Optional[str]]]) -> None:
        attrs_dict = dict(attrs)

        if tag in self.BLOCK_TAGS:
            self._text_parts.append("\n")

        if tag == "tr":
            self._in_tr = True
            self._current_row = []

        elif tag in ("td", "th"):
            self._in_td = True
            self._current_cell_text_parts = []
            self._current_cell_links = []

        elif tag == "a":
            href = attrs_dict.get("href")
            self._current_anchor_href = href
            self._current_anchor_text_parts = []
            if self._in_td and href:
                self._current_cell_links.append(href)

        elif tag == "sub":
            # 让 IC50 / EC50 在纯文本里更接近可读形式
            self._text_parts.append("_")

    def handle_endtag(self, tag: str) -> None:
        if tag == "a":
            anchor_text = normalize_whitespace("".join(self._current_anchor_text_parts))
            href = self._current_anchor_href
            if href and anchor_text:
                self.links.append({"text": anchor_text, "href": href})
            self._current_anchor_href = None
            self._current_anchor_text_parts = []

        elif tag in ("td", "th"):
            text = normalize_whitespace("".join(self._current_cell_text_parts))
            self._current_row.append({
                "text": text,
                "links": list(self._current_cell_links),
            })
            self._in_td = False
            self._current_cell_text_parts = []
            self._current_cell_links = []

        elif tag == "tr":
            self._in_tr = False
            if self._current_row:
                self.rows.append(self._current_row)
            self._current_row = []
            self._text_parts.append("\n")

        elif tag in self.BLOCK_TAGS:
            self._text_parts.append("\n")

    def handle_data(self, data: str) -> None:
        if not data:
            return
        self._text_parts.append(data)
        if self._in_td:
            self._current_cell_text_parts.append(data)
        if self._current_anchor_href is not None:
            self._current_anchor_text_parts.append(data)

    def get_text(self) -> str:
        text = "".join(self._text_parts)
        text = text.replace("\xa0", " ")
        lines = [normalize_whitespace(line) for line in text.splitlines()]
        lines = [line for line in lines if line]
        return "\n".join(lines).strip()


# ============================================================
# 五、列表页 / 详情页字段结构
# ============================================================

@dataclass
class ListRecord:
    list_page_no: int
    list_page_url: str
    peptide_row_id: str
    name: str
    sequence: str
    chem_mass: str
    monois_mass: str
    activity_value_raw: str
    assay_type_raw: str
    detail_url: str


@dataclass
class DetailRecord:
    detail_id: str
    name: str
    sequence: str
    function_text: str
    residue_count_raw: str
    activity_code: str
    activity_name: str
    chemical_mass: str
    monoisotopic_mass: str
    activity_measure_label: str
    activity_measure_value_raw: str
    authors: str
    title: str
    year: str
    source: str
    smiles: str
    inchi: str
    inchikey: str
    database_reference_text: str
    raw_text_file: str
    raw_html_file: str
    detail_url: str


# ============================================================
# 六、列表页解析
# ============================================================

def parse_total_pages_from_text(text: str) -> Optional[int]:
    m = re.search(r"Page:\s*\d+\s*/\s*(\d+)", text)
    if m:
        return safe_int(m.group(1))
    m = re.search(r"Page\s+\d+\s*/\s*(\d+)", text)
    if m:
        return safe_int(m.group(1))
    return None


def find_next_page_url(parser: TableAndLinkParser, current_url: str) -> Optional[str]:
    """从当前页所有链接里找到文字为 Next 的链接。"""
    for item in parser.links:
        if normalize_whitespace(item.get("text", "")).lower() == "next":
            href = item.get("href", "")
            if href:
                return urljoin(current_url, href)
    return None


def parse_list_records_from_parser(parser: TableAndLinkParser, current_url: str, page_no: int) -> List[ListRecord]:
    """
    从列表页中提取记录。

    优先策略：
    1. 直接按表格 <tr>/<td> 解析；
    2. 若站点页面结构略变，表格解析失败，则退回到“纯文本行 + 链接顺序”解析。

    这样做是为了提高脚本在站点小改版后的可用性。
    """
    records: List[ListRecord] = []

    # ------------------------------
    # A. 优先使用表格行解析
    # ------------------------------
    for row in parser.rows:
        cell_texts = [normalize_whitespace(str(cell.get("text", ""))) for cell in row]
        if len(cell_texts) < 8:
            continue
        if not cell_texts[0]:
            continue
        if "Peptide Data" not in cell_texts[0]:
            continue

        detail_links = row[0].get("links", []) if row else []
        if not detail_links:
            continue

        detail_url = urljoin(current_url, str(detail_links[0]))

        peptide_row_id = cell_texts[1]
        name = cell_texts[2]
        sequence = cell_texts[3]
        chem_mass = cell_texts[4]
        monois_mass = cell_texts[5]
        activity_value_raw = cell_texts[6]
        assay_type_raw = cell_texts[7]

        record = ListRecord(
            list_page_no=page_no,
            list_page_url=current_url,
            peptide_row_id=peptide_row_id,
            name=name,
            sequence=sequence,
            chem_mass=chem_mass,
            monois_mass=monois_mass,
            activity_value_raw=activity_value_raw,
            assay_type_raw=assay_type_raw,
            detail_url=detail_url,
        )
        records.append(record)

    if records:
        return records

    # ------------------------------
    # B. 回退：按纯文本行解析
    # ------------------------------
    text_lines = [x.strip() for x in parser.get_text().splitlines() if x.strip()]
    peptide_links = [
        urljoin(current_url, item["href"])
        for item in parser.links
        if normalize_whitespace(item.get("text", "")) == "Peptide Data"
    ]

    line_candidates = [line for line in text_lines if line.startswith("Peptide Data ")]
    pair_count = min(len(line_candidates), len(peptide_links))

    for idx in range(pair_count):
        line = line_candidates[idx]
        detail_url = peptide_links[idx]

        # 典型格式：
        # Peptide Data 10989 ACE inhibitor WP 301.3396 301.1422 217.00 IC_50
        parts = line.split()
        if len(parts) < 8:
            continue
        if parts[0] != "Peptide" or parts[1] != "Data":
            continue

        peptide_row_id = parts[2]
        assay_type_raw = parts[-1].replace("{", "").replace("}", "")
        activity_value_raw = parts[-2]
        monois_mass = parts[-3]
        chem_mass = parts[-4]
        sequence = parts[-5]
        name = " ".join(parts[3:-5]).strip()

        records.append(ListRecord(
            list_page_no=page_no,
            list_page_url=current_url,
            peptide_row_id=peptide_row_id,
            name=name,
            sequence=sequence,
            chem_mass=chem_mass,
            monois_mass=monois_mass,
            activity_value_raw=activity_value_raw,
            assay_type_raw=assay_type_raw,
            detail_url=detail_url,
        ))

    return records


# ============================================================
# 七、详情页解析
# ============================================================

def extract_block_after_label(lines: List[str], label: str, stop_labels: set) -> str:
    """
    给定一个标签行，如 “Authors” 或 “Title”，取其后面的连续内容，
    直到遇到下一个停止标签。
    """
    for idx, line in enumerate(lines):
        if line == label:
            collected: List[str] = []
            for j in range(idx + 1, len(lines)):
                next_line = lines[j].strip()
                if not next_line:
                    continue
                if next_line in stop_labels:
                    break
                # 页面里常出现 * * * 分隔线，这里遇到也停止
                if re.fullmatch(r"\*+( \*+)*", next_line):
                    break
                collected.append(next_line)
            return " ".join(collected).strip()
    return ""


def extract_single_line_value(lines: List[str], prefix: str) -> str:
    for line in lines:
        if line.startswith(prefix):
            return normalize_whitespace(line[len(prefix):].strip())
    return ""


def extract_next_nonempty_line(lines: List[str], label: str) -> str:
    for idx, line in enumerate(lines):
        if line == label:
            for j in range(idx + 1, len(lines)):
                val = lines[j].strip()
                if val:
                    return val
    return ""


def parse_year_and_source(lines: List[str]) -> Tuple[str, str]:
    """
    页面里 Year / Source 常见样式是：
    Year Source
    2006 Journal

    有时也可能拆成多行，因此这里尽量宽松解析。
    """
    for idx, line in enumerate(lines):
        if line == "Year Source":
            for j in range(idx + 1, len(lines)):
                nxt = lines[j].strip()
                if not nxt:
                    continue
                m = re.match(r"^(\d{4})\s+(.+)$", nxt)
                if m:
                    return m.group(1), m.group(2).strip()
                if re.fullmatch(r"\d{4}", nxt):
                    year = nxt
                    source = ""
                    for k in range(j + 1, len(lines)):
                        after = lines[k].strip()
                        if after:
                            source = after
                            break
                    return year, source
    return "", ""


def extract_database_reference_text(lines: List[str]) -> str:
    for idx, line in enumerate(lines):
        if line == "Database reference:":
            collected: List[str] = []
            for j in range(idx + 1, len(lines)):
                nxt = lines[j].strip()
                if not nxt:
                    continue
                collected.append(nxt)
            return " | ".join(collected)
    return ""


def parse_detail_record_from_text(detail_url: str, raw_text: str, raw_text_file: str, raw_html_file: str) -> DetailRecord:
    """从详情页纯文本中提取主要字段。"""
    lines = [normalize_whitespace(x) for x in raw_text.splitlines()]
    lines = [x for x in lines if x]

    stop_labels = {
        "Authors", "Title", "Year Source", "Additional information:", "Database reference:",
        "Name", "sequence", "Function:", "Number of residues", "Activity code", "Activity :",
        "IC50 :", "EC50 :", "IC 50 :", "EC 50 :",
    }

    # 1）ID / Name / sequence / Function / Number of residues / Activity code / Activity
    detail_id = extract_single_line_value(lines, "ID ")
    name = extract_single_line_value(lines, "Name ")
    sequence = extract_next_nonempty_line(lines, "sequence")
    function_text = extract_next_nonempty_line(lines, "Function:")
    residue_count_raw = extract_next_nonempty_line(lines, "Number of residues")
    activity_code = extract_next_nonempty_line(lines, "Activity code")
    activity_name = extract_next_nonempty_line(lines, "Activity :")

    # 2）Chemical mass / Monoisotopic mass
    chemical_mass = ""
    monoisotopic_mass = ""
    for line in lines:
        m = re.search(r"Chemical mass\s+([0-9.]+)\s+Monoisotopic mass\s+([0-9.]+)", line)
        if m:
            chemical_mass = m.group(1)
            monoisotopic_mass = m.group(2)
            break

    # 3）活性值标签与数值（优先识别 IC50，再识别 EC50）
    activity_measure_label = ""
    activity_measure_value_raw = ""
    for label in ["IC50 :", "IC 50 :", "EC50 :", "EC 50 :"]:
        value = extract_next_nonempty_line(lines, label)
        if value:
            activity_measure_label = label.replace(" ", "")
            activity_measure_value_raw = value
            break

    # 4）文献信息
    authors = extract_block_after_label(lines, "Authors", stop_labels)
    title = extract_block_after_label(lines, "Title", stop_labels)
    year, source = parse_year_and_source(lines)

    # 5）结构字符串
    smiles = extract_single_line_value(lines, "SMILES:")

    inchi = ""
    for line in lines:
        if line.startswith("InChI="):
            inchi = line
            break

    inchikey = extract_single_line_value(lines, "InChIKey:")

    # 6）交叉库参考
    database_reference_text = extract_database_reference_text(lines)

    return DetailRecord(
        detail_id=detail_id,
        name=name,
        sequence=sequence,
        function_text=function_text,
        residue_count_raw=residue_count_raw,
        activity_code=activity_code,
        activity_name=activity_name,
        chemical_mass=chemical_mass,
        monoisotopic_mass=monoisotopic_mass,
        activity_measure_label=activity_measure_label,
        activity_measure_value_raw=activity_measure_value_raw,
        authors=authors,
        title=title,
        year=year,
        source=source,
        smiles=smiles,
        inchi=inchi,
        inchikey=inchikey,
        database_reference_text=database_reference_text,
        raw_text_file=raw_text_file,
        raw_html_file=raw_html_file,
        detail_url=detail_url,
    )


# ============================================================
# 八、核心抓取流程
# ============================================================

def save_list_page_files(paths: Dict[str, str], page_no: int, url: str, html: str, text: str) -> Tuple[str, str]:
    html_file = os.path.join(paths["raw_html_list"], f"list_page_{page_no:04d}.html")
    txt_file = os.path.join(paths["raw_txt_list"], f"list_page_{page_no:04d}.txt")

    if SAVE_RAW_HTML:
        if OVERWRITE_EXISTING_FILES or (not os.path.exists(html_file)):
            write_text(html_file, html, encoding="utf-8")

    if SAVE_LIST_PAGE_TXT:
        header = f"BIOPEP-UWM experimental 列表页\n页码: {page_no}\nURL: {url}\n\n"
        if OVERWRITE_EXISTING_FILES or (not os.path.exists(txt_file)):
            write_text(txt_file, header + text, encoding="utf-8-sig")

    return html_file, txt_file


def save_detail_files(paths: Dict[str, str], detail_id: str, html: str, text: str) -> Tuple[str, str]:
    base_name = f"detail_{slugify_filename(detail_id)}"
    html_file = os.path.join(paths["raw_html_detail"], base_name + ".html")
    txt_file = os.path.join(paths["raw_txt_detail"], base_name + ".txt")

    if SAVE_RAW_HTML:
        if OVERWRITE_EXISTING_FILES or (not os.path.exists(html_file)):
            write_text(html_file, html, encoding="utf-8")

    if SAVE_DETAIL_TXT:
        if OVERWRITE_EXISTING_FILES or (not os.path.exists(txt_file)):
            write_text(txt_file, text, encoding="utf-8-sig")

    return html_file, txt_file


def crawl_list_pages(paths: Dict[str, str]) -> Tuple[List[ListRecord], Dict[str, object]]:
    """
    从 peptide_data.php 首页开始，沿着 Next 链接一路抓到底。

    这样做的好处是：
    - 不依赖我们手工猜分页参数格式；
    - 以后站点如果只改了分页参数写法，只要 Next 链接还在，脚本还能继续用。
    """
    all_target_records: List[ListRecord] = []
    visited_urls = set()

    current_url = LIST_URL
    current_page_no = 1
    total_pages_seen = None
    page_count = 0
    all_rows_count = 0

    while current_url and current_url not in visited_urls:
        if DEBUG_MAX_LIST_PAGES is not None and page_count >= DEBUG_MAX_LIST_PAGES:
            break

        visited_urls.add(current_url)
        page_count += 1

        html = fetch_html(current_url)
        parser = TableAndLinkParser()
        parser.feed(html)
        page_text = parser.get_text()

        if total_pages_seen is None:
            total_pages_seen = parse_total_pages_from_text(page_text)

        save_list_page_files(paths, current_page_no, current_url, html, page_text)

        page_records = parse_list_records_from_parser(parser, current_url, current_page_no)
        all_rows_count += len(page_records)

        target_records = [r for r in page_records if normalize_whitespace(r.name) == TARGET_ACTIVITY_NAME]
        all_target_records.extend(target_records)

        append_log(
            paths["log_file"],
            f"[列表页] page={current_page_no} rows={len(page_records)} target={len(target_records)} url={current_url}"
        )

        next_url = find_next_page_url(parser, current_url)
        if not next_url:
            break

        current_url = next_url
        current_page_no += 1
        time.sleep(REQUEST_SLEEP_BETWEEN_PAGES)

    summary = {
        "list_pages_crawled": page_count,
        "list_total_pages_seen": total_pages_seen,
        "list_total_rows_parsed": all_rows_count,
        "target_rows_found": len(all_target_records),
    }
    return all_target_records, summary


def crawl_detail_pages(paths: Dict[str, str], list_records: List[ListRecord]) -> Tuple[List[DetailRecord], Dict[str, object]]:
    detail_records: List[DetailRecord] = []
    success_count = 0
    failed_count = 0

    for idx, list_record in enumerate(list_records, start=1):
        if DEBUG_MAX_DETAIL_RECORDS is not None and idx > DEBUG_MAX_DETAIL_RECORDS:
            break

        try:
            # 先下载 HTML
            html = fetch_html(list_record.detail_url)
            parser = TableAndLinkParser()
            parser.feed(html)
            raw_text = parser.get_text()

            # 先从文本里解析 detail_id；如果暂时拿不到，就用列表页 row_id 当备用
            detail_id_match = re.search(r"\bID\s+(\d+)\b", raw_text)
            detail_id = detail_id_match.group(1) if detail_id_match else str(list_record.peptide_row_id)

            html_file_candidate = os.path.join(paths["raw_html_detail"], f"detail_{slugify_filename(detail_id)}.html")
            if SKIP_ALREADY_DOWNLOADED_DETAIL_HTML and os.path.exists(html_file_candidate) and not OVERWRITE_EXISTING_FILES:
                # 已存在时，直接复用文件内容，减少重复写入；但为了简洁，这里仍按当前下载结果继续解析
                pass

            raw_html_file, raw_text_file = save_detail_files(paths, detail_id, html, raw_text)
            detail_record = parse_detail_record_from_text(
                detail_url=list_record.detail_url,
                raw_text=raw_text,
                raw_text_file=raw_text_file,
                raw_html_file=raw_html_file,
            )

            # 补一层交叉核对：把列表页已有字段也一起带进详情 CSV 的最终导出时再合并
            detail_records.append(detail_record)
            success_count += 1

            append_log(
                paths["log_file"],
                f"[详情页] ok index={idx}/{len(list_records)} detail_id={detail_record.detail_id} url={list_record.detail_url}"
            )
        except Exception as exc:
            failed_count += 1
            append_log(
                paths["log_file"],
                f"[详情页] FAILED index={idx}/{len(list_records)} url={list_record.detail_url} error={exc}"
            )
            append_log(paths["log_file"], traceback.format_exc())

        time.sleep(REQUEST_SLEEP_BETWEEN_DETAILS)

    summary = {
        "detail_success_count": success_count,
        "detail_failed_count": failed_count,
    }
    return detail_records, summary


# ============================================================
# 九、导出 CSV
# ============================================================

def export_list_csv(paths: Dict[str, str], list_records: List[ListRecord]) -> None:
    rows = []
    for r in list_records:
        rows.append({
            "列表页页码": r.list_page_no,
            "列表页URL": r.list_page_url,
            "列表页行ID": r.peptide_row_id,
            "名称": r.name,
            "序列": r.sequence,
            "化学质量": r.chem_mass,
            "单同位素质量": r.monois_mass,
            "列表页活性值原文": r.activity_value_raw,
            "列表页测定类型原文": r.assay_type_raw,
            "详情页URL": r.detail_url,
        })

    fieldnames = [
        "列表页页码",
        "列表页URL",
        "列表页行ID",
        "名称",
        "序列",
        "化学质量",
        "单同位素质量",
        "列表页活性值原文",
        "列表页测定类型原文",
        "详情页URL",
    ]
    write_csv(paths["ace_list_csv"], rows, fieldnames)


def export_detail_csv(paths: Dict[str, str], list_records: List[ListRecord], detail_records: List[DetailRecord]) -> None:
    """
    最终详情 CSV 同时保留：
    - 列表页已有字段（方便回溯分页来源）
    - 详情页解析字段（方便你直接查阅）
    """
    list_map: Dict[str, ListRecord] = {}
    for item in list_records:
        # 用详情 URL 作为最稳定的合并键
        list_map[item.detail_url] = item

    rows: List[Dict[str, object]] = []
    for d in detail_records:
        src = list_map.get(d.detail_url)
        rows.append({
            "列表页页码": src.list_page_no if src else "",
            "列表页URL": src.list_page_url if src else "",
            "列表页行ID": src.peptide_row_id if src else "",
            "列表页名称": src.name if src else "",
            "列表页序列": src.sequence if src else "",
            "列表页化学质量": src.chem_mass if src else "",
            "列表页单同位素质量": src.monois_mass if src else "",
            "列表页活性值原文": src.activity_value_raw if src else "",
            "列表页测定类型原文": src.assay_type_raw if src else "",
            "详情ID": d.detail_id,
            "详情名称": d.name,
            "详情序列": d.sequence,
            "功能说明": d.function_text,
            "残基数原文": d.residue_count_raw,
            "活性代码": d.activity_code,
            "活性名称": d.activity_name,
            "详情化学质量": d.chemical_mass,
            "详情单同位素质量": d.monoisotopic_mass,
            "活性指标标签": d.activity_measure_label,
            "活性指标原文": d.activity_measure_value_raw,
            "作者": d.authors,
            "题目": d.title,
            "年份": d.year,
            "来源类型": d.source,
            "SMILES": d.smiles,
            "InChI": d.inchi,
            "InChIKey": d.inchikey,
            "交叉库参考原文": d.database_reference_text,
            "详情页URL": d.detail_url,
            "原始TXT文件": d.raw_text_file,
            "原始HTML文件": d.raw_html_file,
        })

    fieldnames = [
        "列表页页码",
        "列表页URL",
        "列表页行ID",
        "列表页名称",
        "列表页序列",
        "列表页化学质量",
        "列表页单同位素质量",
        "列表页活性值原文",
        "列表页测定类型原文",
        "详情ID",
        "详情名称",
        "详情序列",
        "功能说明",
        "残基数原文",
        "活性代码",
        "活性名称",
        "详情化学质量",
        "详情单同位素质量",
        "活性指标标签",
        "活性指标原文",
        "作者",
        "题目",
        "年份",
        "来源类型",
        "SMILES",
        "InChI",
        "InChIKey",
        "交叉库参考原文",
        "详情页URL",
        "原始TXT文件",
        "原始HTML文件",
    ]
    write_csv(paths["ace_detail_csv"], rows, fieldnames)


# ============================================================
# 十、README 与汇总信息输出
# ============================================================

def build_readme_text(paths: Dict[str, str], summary: Dict[str, object]) -> str:
    return f"""BIOPEP-UWM experimental ACE inhibitor 原始抓取说明
===============================

一、默认输出位置
{paths['root']}

二、目录说明
1. raw_html/list_pages
   - 列表页原始 HTML
2. raw_html/detail_pages
   - 详情页原始 HTML
3. raw_txt/list_pages
   - 列表页纯文本版本，便于快速查看
4. raw_txt/detail_cards
   - 每条详情页纯文本版本，便于快速核对
5. raw_tables_csv
   - biopep_uwm_experimental_ace_list.csv
     列表页中筛到的 ACE inhibitor 条目索引表
   - biopep_uwm_experimental_ace_detail.csv
     详情页主要字段汇总表
6. meta
   - run_log.txt
   - crawl_summary.json
   - 本说明文件

三、本次抓取摘要
{json.dumps(summary, ensure_ascii=False, indent=2)}

四、当前脚本目标
- 只做 BIOPEP-UWM experimental 中 ACE inhibitor 原始数据下载
- 不做清洗、标准化、单位换算、去重或多库融合

五、后续扩展建议
- 若以后要抓其他活性，只需修改 TARGET_ACTIVITY_NAME / TARGET_ACTIVITY_CODE
- 若以后要做标准化，可在此脚本产物基础上新建 standardized 脚本，而不是直接改 raw 数据
"""


# ============================================================
# 十一、主程序
# ============================================================

def main() -> int:
    output_root = get_output_root()
    paths = build_output_paths(output_root)
    ensure_directories(paths)

    append_log(paths["log_file"], "=" * 80)
    append_log(paths["log_file"], "开始抓取 BIOPEP-UWM experimental ACE inhibitor 原始数据")
    append_log(paths["log_file"], f"输出目录：{paths['root']}")
    append_log(paths["log_file"], f"列表页入口：{LIST_URL}")
    append_log(paths["log_file"], f"目标活动名称：{TARGET_ACTIVITY_NAME}")
    append_log(paths["log_file"], f"目标活动代码：{TARGET_ACTIVITY_CODE}")

    started_at = time.strftime("%Y-%m-%d %H:%M:%S")

    # 第一步：抓列表页，筛出 ACE inhibitor
    list_records, list_summary = crawl_list_pages(paths)
    export_list_csv(paths, list_records)

    # 第二步：抓详情页
    detail_records, detail_summary = crawl_detail_pages(paths, list_records)
    export_detail_csv(paths, list_records, detail_records)

    ended_at = time.strftime("%Y-%m-%d %H:%M:%S")

    final_summary = {
        "started_at": started_at,
        "ended_at": ended_at,
        "output_root": paths["root"],
        "target_activity_name": TARGET_ACTIVITY_NAME,
        "target_activity_code": TARGET_ACTIVITY_CODE,
        **list_summary,
        **detail_summary,
        "list_csv": paths["ace_list_csv"],
        "detail_csv": paths["ace_detail_csv"],
    }

    write_json(paths["summary_json"], final_summary)
    write_text(paths["readme_txt"], build_readme_text(paths, final_summary), encoding="utf-8-sig")

    append_log(paths["log_file"], "抓取完成。")
    append_log(paths["log_file"], json.dumps(final_summary, ensure_ascii=False, indent=2))

    print("=" * 80)
    print("BIOPEP-UWM experimental ACE inhibitor 原始抓取完成")
    print(f"输出根目录: {paths['root']}")
    print(f"列表索引 CSV: {paths['ace_list_csv']}")
    print(f"详情汇总 CSV: {paths['ace_detail_csv']}")
    print(f"运行摘要 JSON: {paths['summary_json']}")
    print("=" * 80)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("用户中断。")
        sys.exit(130)
    except Exception as exc:
        print("程序运行失败：", exc)
        traceback.print_exc()
        sys.exit(1)
