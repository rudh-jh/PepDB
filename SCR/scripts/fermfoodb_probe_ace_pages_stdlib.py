#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FermFooDb ACE 页面探测脚本（纯标准库版）
=====================================

用途
----
1. 探测 FermFooDb 站点是否可访问
2. 小规模抓取站内页面，观察链接结构与 ACE 相关关键词命中
3. 为后续 ACE 定向 downloader 提供页面样本和结构线索

默认输出
--------
DB/raw/ace/databases/FermFooDb/experimental/
├── meta/
│   ├── fermfoodb_probe_summary.json
│   ├── discovered_pages.csv
│   ├── keyword_hits.csv
│   └── errors.csv
├── raw_html/
│   ├── list_pages/
│   └── detail_pages/
└── raw_txt/
    ├── list_pages/
    └── detail_cards/

说明
----
- 这不是正式 downloader，只是 probe
- 不做标准化，不做去重，不做单位换算
- 默认只探测少量页面
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
from collections import Counter, deque
from html.parser import HTMLParser
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
from urllib.error import HTTPError, URLError
from urllib.parse import parse_qs, urljoin, urlparse, urlunparse
from urllib.request import Request, urlopen


DEFAULT_BASE_URLS = [
    "https://webs.iiitd.edu.in/raghava/fermfoodb/",
    "http://webs.iiitd.edu.in/raghava/fermfoodb/",
]

DEFAULT_KEYWORDS = [
    "ace",
    "ace inhibitory",
    "ace-inhibitory",
    "angiotensin",
    "angiotensin converting enzyme",
    "angiotensin-converting enzyme",
    "antihypertensive",
    "hypertension",
    "blood pressure",
]

DEFAULT_PAGE_KEYWORDS = [
    "peptide",
    "activity",
    "ic50",
    "function",
    "matrix",
    "culture",
    "protein",
    "fermentation",
    "hydrolysis",
]

USER_AGENT = (
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
    "AppleWebKit/537.36 (KHTML, like Gecko) "
    "Chrome/124.0 Safari/537.36"
)


def infer_project_root(script_file: Path) -> Path:
    try:
        return script_file.resolve().parents[2]
    except IndexError:
        return script_file.resolve().parent


def normalize_spaces(text: str) -> str:
    return re.sub(r"\s+", " ", text.replace("\r", " ").replace("\n", " ")).strip()


def safe_filename_from_url(url: str, prefix: str = "page") -> str:
    h = hashlib.md5(url.encode("utf-8")).hexdigest()[:12]
    parsed = urlparse(url)
    path = parsed.path.strip("/").replace("/", "_")
    if not path:
        path = "root"
    path = re.sub(r"[^A-Za-z0-9._-]+", "_", path)[:60]
    return f"{prefix}_{path}_{h}"


def is_html_content_type(content_type: Optional[str]) -> bool:
    if not content_type:
        return True
    low = content_type.lower()
    return ("html" in low) or ("text/" in low)


def same_host(url: str, hosts: Set[str]) -> bool:
    parsed = urlparse(url)
    return parsed.netloc.lower() in hosts


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
    fragment = ""
    return urlunparse((scheme, netloc, path, "", query, fragment))


def contains_any_keyword(text: str, keywords: List[str]) -> List[str]:
    low = text.lower()
    hits = []
    for kw in keywords:
        if kw.lower() in low:
            hits.append(kw)
    return hits


class LinkParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links: List[Tuple[str, str]] = []
        self.title_parts: List[str] = []
        self._in_title = False
        self.text_parts: List[str] = []

    def handle_starttag(self, tag, attrs):
        attrs_dict = dict(attrs)
        tag_low = tag.lower()
        if tag_low == "a":
            href = attrs_dict.get("href")
            if href:
                self.links.append((href, attrs_dict.get("title", "")))
        elif tag_low == "title":
            self._in_title = True

    def handle_endtag(self, tag):
        if tag.lower() == "title":
            self._in_title = False

    def handle_data(self, data):
        text = normalize_spaces(html.unescape(data))
        if not text:
            return
        if self._in_title:
            self.title_parts.append(text)
        self.text_parts.append(text)

    @property
    def title(self) -> str:
        return normalize_spaces(" ".join(self.title_parts))

    @property
    def text(self) -> str:
        return normalize_spaces(" ".join(self.text_parts))


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
                "url": url,
                "final_url": final_url,
                "status": status,
                "content_type": content_type,
                "elapsed_sec": round(elapsed, 4),
                "text": text,
                "error_type": "",
                "error_message": "",
            }

    except HTTPError as e:
        return {
            "ok": False,
            "url": url,
            "final_url": url,
            "status": getattr(e, "code", None),
            "content_type": "",
            "elapsed_sec": round(time.time() - started_at, 4),
            "text": "",
            "error_type": "HTTPError",
            "error_message": str(e),
        }
    except URLError as e:
        return {
            "ok": False,
            "url": url,
            "final_url": url,
            "status": None,
            "content_type": "",
            "elapsed_sec": round(time.time() - started_at, 4),
            "text": "",
            "error_type": "URLError",
            "error_message": str(e),
        }
    except socket.timeout as e:
        return {
            "ok": False,
            "url": url,
            "final_url": url,
            "status": None,
            "content_type": "",
            "elapsed_sec": round(time.time() - started_at, 4),
            "text": "",
            "error_type": "Timeout",
            "error_message": str(e),
        }
    except Exception as e:
        return {
            "ok": False,
            "url": url,
            "final_url": url,
            "status": None,
            "content_type": "",
            "elapsed_sec": round(time.time() - started_at, 4),
            "text": "",
            "error_type": e.__class__.__name__,
            "error_message": str(e),
        }


def parse_html_page(url: str, html_text: str) -> Dict[str, object]:
    parser = LinkParser()
    parser.feed(html_text)

    abs_links = []
    for href, title_attr in parser.links:
        abs_url = urljoin(url, href)
        abs_url = canonicalize_url(abs_url)
        abs_links.append({
            "url": abs_url,
            "title_attr": normalize_spaces(title_attr),
        })

    return {
        "title": parser.title,
        "text": parser.text,
        "links": abs_links,
    }


def classify_page(url: str, title: str, text: str) -> str:
    blob = " ".join([url, title, text[:4000]]).lower()

    detail_patterns = [
        "detail", "peptide", "peptide sequence", "ic50", "assay", "matrix",
        "culture", "protein", "doi", "entry", "record",
    ]
    list_patterns = [
        "search", "advance", "browse", "activity", "result", "results",
        "page=", "list", "table", "hydrolysis", "experiment",
    ]

    detail_score = sum(1 for p in detail_patterns if p in blob)
    list_score = sum(1 for p in list_patterns if p in blob)

    parsed = urlparse(url)
    query_keys = set(parse_qs(parsed.query).keys())

    if {"id", "seqid", "pid", "detail"} & query_keys:
        detail_score += 2
    if {"page", "search", "act", "activity", "keyword", "type"} & query_keys:
        list_score += 2

    if detail_score > list_score:
        return "detail_like"
    if list_score > detail_score:
        return "list_like"
    return "unknown"


def save_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def save_probe_page(base_dir: Path, page_type: str, url: str, html_text: str, plain_text: str) -> Tuple[Path, Path]:
    stem = safe_filename_from_url(url, prefix=page_type)
    if page_type == "detail_like":
        html_path = base_dir / "raw_html" / "detail_pages" / f"{stem}.html"
        txt_path = base_dir / "raw_txt" / "detail_cards" / f"{stem}.txt"
    else:
        html_path = base_dir / "raw_html" / "list_pages" / f"{stem}.html"
        txt_path = base_dir / "raw_txt" / "list_pages" / f"{stem}.txt"

    save_text(html_path, html_text)
    save_text(txt_path, plain_text)
    return html_path, txt_path


def write_csv(path: Path, rows: List[Dict[str, object]], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def build_seed_urls(start_url: str) -> List[str]:
    """
    给 probe 一些常见种子页面，提升更快命中 ACE / Browse / Search 的概率
    """
    seeds = [
        start_url,
        urljoin(start_url, "index.html"),
        urljoin(start_url, "search.php"),
        urljoin(start_url, "browse.php"),
        urljoin(start_url, "advance.php"),
        urljoin(start_url, "activity.php"),
    ]
    normalized = []
    seen = set()
    for x in seeds:
        cx = canonicalize_url(x)
        if cx not in seen:
            seen.add(cx)
            normalized.append(cx)
    return normalized


def main() -> int:
    parser = argparse.ArgumentParser(description="FermFooDb ACE 页面探测脚本（纯标准库版）")
    parser.add_argument(
        "--project-root",
        type=str,
        default=None,
        help="PepDB 项目根目录，例如 E:\\MYS\\PepDB",
    )
    parser.add_argument(
        "--base-url",
        type=str,
        default=None,
        help="单个入口 URL。不给时会按内置候选入口依次尝试。",
    )
    parser.add_argument(
        "--max-pages",
        type=int,
        default=30,
        help="最多探测多少个页面，默认 30。",
    )
    parser.add_argument(
        "--max-depth",
        type=int,
        default=2,
        help="BFS 最大深度，默认 2。",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=20,
        help="单次请求超时秒数，默认 20。",
    )
    parser.add_argument(
        "--sleep-sec",
        type=float,
        default=0.3,
        help="相邻请求间隔秒数，默认 0.3。",
    )
    parser.add_argument(
        "--save-limit",
        type=int,
        default=12,
        help="最多保存多少个样本页面 HTML/TXT，默认 12。",
    )
    args = parser.parse_args()

    script_file = Path(__file__)
    project_root = Path(args.project_root) if args.project_root else infer_project_root(script_file)

    base_dir = project_root / "DB" / "raw" / "ace" / "databases" / "FermFooDb" / "experimental"
    meta_dir = base_dir / "meta"
    meta_dir.mkdir(parents=True, exist_ok=True)

    summary_path = meta_dir / "fermfoodb_probe_summary.json"
    discovered_csv = meta_dir / "discovered_pages.csv"
    hits_csv = meta_dir / "keyword_hits.csv"
    errors_csv = meta_dir / "errors.csv"

    candidate_base_urls = [args.base_url] if args.base_url else DEFAULT_BASE_URLS[:]

    start_success = None
    start_trials = []

    for candidate in candidate_base_urls:
        trial = fetch_url(candidate, timeout=args.timeout, sleep_sec=0.0)
        start_trials.append({
            "candidate_url": candidate,
            "ok": trial["ok"],
            "status": trial["status"],
            "final_url": trial["final_url"],
            "content_type": trial["content_type"],
            "error_type": trial["error_type"],
            "error_message": trial["error_message"],
            "elapsed_sec": trial["elapsed_sec"],
        })
        if trial["ok"] and is_html_content_type(trial["content_type"]):
            start_success = trial
            break

    discovered_rows: List[Dict[str, object]] = []
    hit_rows: List[Dict[str, object]] = []
    error_rows: List[Dict[str, object]] = []

    if not start_success:
        summary = {
            "probe_ok": False,
            "message": "all candidate base urls failed",
            "candidate_base_urls": candidate_base_urls,
            "start_trials": start_trials,
            "max_pages": args.max_pages,
            "max_depth": args.max_depth,
            "timeout": args.timeout,
            "sleep_sec": args.sleep_sec,
        }
        summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
        write_csv(
            errors_csv,
            [
                {
                    "stage": "start_probe",
                    "url": x["candidate_url"],
                    "status": x["status"],
                    "error_type": x["error_type"],
                    "error_message": x["error_message"],
                    "elapsed_sec": x["elapsed_sec"],
                }
                for x in start_trials
            ],
            ["stage", "url", "status", "error_type", "error_message", "elapsed_sec"],
        )
        print("=" * 80)
        print("FermFooDb probe 未找到可访问入口")
        print(f"已尝试：{candidate_base_urls}")
        print(f"summary: {summary_path}")
        print(f"errors : {errors_csv}")
        print("=" * 80)
        return 1

    start_url = canonicalize_url(start_success["final_url"])
    start_host = urlparse(start_url).netloc.lower()
    allowed_hosts = {start_host}

    visited: Set[str] = set()
    q = deque()

    for seed in build_seed_urls(start_url):
        q.append((seed, 0))

    saved_count = 0
    total_requests = 0

    while q and len(visited) < args.max_pages:
        current_url, depth = q.popleft()
        current_url = canonicalize_url(current_url)

        if current_url in visited:
            continue
        if depth > args.max_depth:
            continue
        if not same_host(current_url, allowed_hosts):
            continue

        result = fetch_url(current_url, timeout=args.timeout, sleep_sec=args.sleep_sec)
        total_requests += 1
        visited.add(current_url)

        if not result["ok"]:
            error_rows.append({
                "stage": "crawl",
                "url": current_url,
                "status": result["status"],
                "error_type": result["error_type"],
                "error_message": result["error_message"],
                "elapsed_sec": result["elapsed_sec"],
            })
            continue

        if not is_html_content_type(result["content_type"]):
            discovered_rows.append({
                "url": current_url,
                "final_url": result["final_url"],
                "depth": depth,
                "status": result["status"],
                "content_type": result["content_type"],
                "page_type": "non_html",
                "title": "",
                "text_len": 0,
                "link_count": 0,
                "keyword_hit_count": 0,
                "keyword_hits": "",
                "html_path": "",
                "txt_path": "",
                "elapsed_sec": result["elapsed_sec"],
            })
            continue

        parsed = parse_html_page(result["final_url"], result["text"])
        text_blob = parsed["text"]
        title = parsed["title"]
        page_type = classify_page(result["final_url"], title, text_blob)
        keyword_hits = contains_any_keyword(" ".join([result["final_url"], title, text_blob]), DEFAULT_KEYWORDS)
        page_keyword_hits = contains_any_keyword(" ".join([title, text_blob]), DEFAULT_PAGE_KEYWORDS)

        html_path = ""
        txt_path = ""
        if saved_count < args.save_limit:
            saved_html_path, saved_txt_path = save_probe_page(
                base_dir=base_dir,
                page_type=page_type,
                url=result["final_url"],
                html_text=result["text"],
                plain_text=text_blob,
            )
            html_path = str(saved_html_path)
            txt_path = str(saved_txt_path)
            saved_count += 1

        discovered_rows.append({
            "url": current_url,
            "final_url": result["final_url"],
            "depth": depth,
            "status": result["status"],
            "content_type": result["content_type"],
            "page_type": page_type,
            "title": title,
            "text_len": len(text_blob),
            "link_count": len(parsed["links"]),
            "keyword_hit_count": len(keyword_hits),
            "keyword_hits": " | ".join(keyword_hits),
            "html_path": html_path,
            "txt_path": txt_path,
            "elapsed_sec": result["elapsed_sec"],
        })

        if keyword_hits or page_keyword_hits:
            hit_rows.append({
                "url": result["final_url"],
                "depth": depth,
                "page_type": page_type,
                "title": title,
                "keyword_hits": " | ".join(keyword_hits),
                "page_keyword_hits": " | ".join(page_keyword_hits),
                "text_preview": text_blob[:500],
            })

        for link in parsed["links"]:
            nxt = link["url"]
            if not nxt.startswith(("http://", "https://")):
                continue
            if not same_host(nxt, allowed_hosts):
                continue
            if canonicalize_url(nxt) in visited:
                continue
            if depth + 1 <= args.max_depth:
                q.append((nxt, depth + 1))

    page_type_counter = Counter(str(x["page_type"]) for x in discovered_rows)
    keyword_hit_pages = sum(1 for x in discovered_rows if int(x["keyword_hit_count"]) > 0)

    summary = {
        "probe_ok": True,
        "base_url_candidates": candidate_base_urls,
        "start_trials": start_trials,
        "selected_start_url": start_url,
        "allowed_hosts": sorted(allowed_hosts),
        "max_pages": args.max_pages,
        "max_depth": args.max_depth,
        "timeout": args.timeout,
        "sleep_sec": args.sleep_sec,
        "save_limit": args.save_limit,
        "total_requests": total_requests,
        "visited_pages": len(discovered_rows),
        "saved_sample_pages": saved_count,
        "keyword_hit_pages": keyword_hit_pages,
        "page_type_counter": dict(page_type_counter),
        "error_count": len(error_rows),
        "ace_keywords": DEFAULT_KEYWORDS,
        "page_keywords": DEFAULT_PAGE_KEYWORDS,
    }

    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")

    write_csv(
        discovered_csv,
        discovered_rows,
        [
            "url", "final_url", "depth", "status", "content_type", "page_type",
            "title", "text_len", "link_count", "keyword_hit_count",
            "keyword_hits", "html_path", "txt_path", "elapsed_sec",
        ],
    )

    write_csv(
        hits_csv,
        hit_rows if hit_rows else [{
            "url": "",
            "depth": "",
            "page_type": "",
            "title": "",
            "keyword_hits": "",
            "page_keyword_hits": "",
            "text_preview": "",
        }],
        ["url", "depth", "page_type", "title", "keyword_hits", "page_keyword_hits", "text_preview"],
    )

    write_csv(
        errors_csv,
        error_rows if error_rows else [{
            "stage": "",
            "url": "",
            "status": "",
            "error_type": "",
            "error_message": "",
            "elapsed_sec": "",
        }],
        ["stage", "url", "status", "error_type", "error_message", "elapsed_sec"],
    )

    print("=" * 80)
    print("FermFooDb probe 完成")
    print(f"项目根目录: {project_root}")
    print(f"起始 URL   : {start_url}")
    print(f"访问页面数 : {len(discovered_rows)}")
    print(f"关键词命中 : {keyword_hit_pages}")
    print(f"保存样本数 : {saved_count}")
    print(f"错误数     : {len(error_rows)}")
    print(f"summary    : {summary_path}")
    print(f"pages csv  : {discovered_csv}")
    print(f"hits csv   : {hits_csv}")
    print(f"errors csv : {errors_csv}")
    print("=" * 80)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())