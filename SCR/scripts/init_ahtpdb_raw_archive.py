# -*- coding: utf-8 -*-
"""
脚本名称：
    init_ahtpdb_raw_archive_retry.py

脚本定位：
    PepDB 项目中用于“原始数据归档层（raw layer）”的 AHTPDB 试点归档脚本（重试增强版）。

与上一版相比，本版新增：
    1. 下载失败自动重试
    2. 支持指数退避（每次失败后延长等待时间）
    3. 支持只补下载某一个数据文件，例如只补 small peptides
    4. 支持跳过页面或跳过数据集，避免每次整包重跑
    5. 检索日志中记录每次重试尝试的细节，便于后续排查
    6. 对官方下载的 txt 原始表，自动额外生成 csv 版本，便于后续查看和处理

当前设计目标：
    先把 AHTPDB 的 raw 层归档尽量补完整，不在 raw 层做字段标准化。

推荐保存位置：
    E:\\MYS\\PepDB\\SCR\\scripts\\init_ahtpdb_raw_archive.py

推荐运行方式：
    1. 全量补跑（默认跳过已存在文件）
       python E:\\MYS\\PepDB\\SCR\\scripts\\init_ahtpdb_raw_archive.py

    2. 只补下载 small peptides
       python E:\\MYS\\PepDB\\SCR\\scripts\\init_ahtpdb_raw_archive.py --skip-pages --dataset-key small

    3. 只重试数据集下载，并把超时放长
       python E:\\MYS\\PepDB\\SCR\\scripts\\init_ahtpdb_raw_archive.py --skip-pages --timeout 60 --retry-count 4

重要原则：
    - 原始数据永不覆盖（默认不覆盖已存在文件）
    - raw 层只做原始保存，不做标准化
    - 官方原始 txt 需要保留
    - csv 只是为了后续查看和处理方便而额外生成，不替代原始 txt
    - 失败会记录到检索日志中，便于人工检查
    - Windows 兼容优先，同时兼容命令行和 PyCharm
"""

from __future__ import annotations

import argparse
import csv
import json
import ssl
import sys
import time
import urllib.request
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


# =====================================================================
# 1. AHTPDB 官方来源配置
# =====================================================================
AHTPDB_SOURCE_ID = "ace_db_ahtpdb_001"
AHTPDB_SOURCE_NAME = "AHTPDB"
AHTPDB_TASK = "ace"
AHTPDB_SOURCE_TYPE = "database"
AHTPDB_EVIDENCE_LEVEL = "curated_db"

AHTPDB_PAGE_URLS: Dict[str, str] = {
    "home": "https://webs.iiitd.edu.in/raghava/ahtpdb/cp_index.php",
    "information": "https://webs.iiitd.edu.in/raghava/ahtpdb/info2.php",
    "download": "https://webs.iiitd.edu.in/raghava/ahtpdb/download.php",
    "guide": "https://webs.iiitd.edu.in/raghava/ahtpdb/main_help.php",
    "publications": "https://webs.iiitd.edu.in/raghava/ahtpdb/pub1.php",
}

AHTPDB_DATASET_URLS: Dict[str, str] = {
    "small": "https://webs.iiitd.edu.in/raghava/ahtpdb/downloads/small.txt",
    "long": "https://webs.iiitd.edu.in/raghava/ahtpdb/downloads/long.txt",
    "ic50": "https://webs.iiitd.edu.in/raghava/ahtpdb/downloads/pepic50.txt",
}

AHTPDB_DATASET_FILENAMES: Dict[str, str] = {
    "small": "ahtpdb_small_peptides_2_5_residues.txt",
    "long": "ahtpdb_long_peptides_6_plus_residues.txt",
    "ic50": "ahtpdb_peptides_with_ic50.txt",
}

AHTPDB_PAGE_EXPORT_FILENAMES: Dict[str, str] = {
    "home": "ahtpdb_home.html",
    "information": "ahtpdb_information.html",
    "download": "ahtpdb_download.html",
    "guide": "ahtpdb_guide.html",
    "publications": "ahtpdb_publications.html",
}


# =====================================================================
# 2. 文本模板
# =====================================================================
def build_ahtpdb_readme(now_str: str) -> str:
    """构造 AHTPDB 原始归档目录 README 内容。"""
    return f"""# AHTPDB 原始归档说明

## 目录用途
本目录用于保存 AHTPDB 与 ACE 抑制肽任务相关的原始获取材料。

## 当前阶段
当前处于原始数据归档阶段。
本目录中的文件只代表“原始来源材料”或“原始保存版本”，不代表标准化结果。

## 当前脚本做了什么
1. 保存 AHTPDB 官方页面原始 HTML
2. 保存 AHTPDB 官方下载文本文件（txt）
3. 将官方下载的 txt 额外转换为 csv，便于后续查看与处理
4. 记录获取时间与来源 URL
5. 更新 PepDB 的 meta 清单文件

## 子目录说明
- page_exports/: 页面原始导出内容（HTML 等）
- raw_tables/: 官方下载的原始文本表格（txt）
- raw_tables_csv/: 由官方 txt 自动转换生成的 csv
- screenshots/: 页面截图（需后续手动补充）
- retrieval_logs/: 获取过程日志、URL 清单、运行摘要

## 当前归档原则
1. 原始 txt 文件不覆盖（除非脚本显式开启覆盖）
2. 原始字段不改写
3. csv 只是对原始 txt 的格式转换，不替代原始 txt
4. 后续清洗、去重、标准化结果不得回写本目录
5. 本目录内容应可追溯到来源页面或来源说明

## 当前不应直接做什么
- 不应把本目录文件直接当作最终分析输入
- 不应在本目录直接修改值、统一单位、去重
- 不应把手工整理后的结果覆盖原文件

## 来源信息
- task: ace
- source_name: AHTPDB
- source_type: database
- archived_at: {now_str}
- archived_by: 当前脚本运行者

## 备注
AHTPDB 当前 raw_tables 目录中的文件，是官方原始保存版本。
raw_tables_csv 目录中的文件，是由 txt 自动转换生成的 csv。
后续 standardized / worksets / analysis 层，应另行生成新文件。
"""


def build_source_note(now_str: str) -> str:
    """构造来源说明文件内容。"""
    lines = [
        "AHTPDB 原始来源说明",
        "=" * 60,
        f"记录时间: {now_str}",
        f"source_id: {AHTPDB_SOURCE_ID}",
        f"task: {AHTPDB_TASK}",
        f"source_type: {AHTPDB_SOURCE_TYPE}",
        f"source_name: {AHTPDB_SOURCE_NAME}",
        f"evidence_level: {AHTPDB_EVIDENCE_LEVEL}",
        "",
        "官方页面 URL:",
    ]
    for name, url in AHTPDB_PAGE_URLS.items():
        lines.append(f"- {name}: {url}")

    lines.append("")
    lines.append("官方下载数据 URL:")
    for name, url in AHTPDB_DATASET_URLS.items():
        lines.append(f"- {name}: {url}")

    lines.extend(
        [
            "",
            "当前说明:",
            "1. 本目录保存的是 AHTPDB 原始页面与官方下载文本。",
            "2. 当前脚本不会在 raw 层做字段标准化。",
            "3. 当前脚本会保留官方 txt，并额外生成 csv。",
            "4. 后续 standardized / worksets / analysis 层应基于本目录另行生成。",
            "5. screenshots 目录默认留空，建议后续手工补充页面结构截图。",
            "6. 本版脚本支持对失败文件进行重试与单文件补下载。",
        ]
    )
    return "\n".join(lines) + "\n"


# =====================================================================
# 3. 通用工具函数
# =====================================================================
def ensure_dir(path: Path) -> None:
    """若目录不存在则创建。"""
    path.mkdir(parents=True, exist_ok=True)


def write_text(path: Path, content: str, overwrite: bool = False) -> Tuple[bool, str]:
    """写入文本文件；默认不覆盖。"""
    if path.exists() and not overwrite:
        return False, "skip_exists"

    ensure_dir(path.parent)
    path.write_text(content, encoding="utf-8")
    return True, "written"


def build_ssl_context(insecure: bool = False) -> ssl.SSLContext:
    """
    构造 SSL 上下文。

    参数：
        insecure:
            False -> 使用默认安全上下文
            True  -> 不校验证书（仅在确有网络兼容问题时临时使用）
    """
    if insecure:
        return ssl._create_unverified_context()  # noqa: SLF001
    return ssl.create_default_context()


def fetch_url_bytes(url: str, timeout: int = 30, insecure: bool = False) -> bytes:
    """下载 URL 内容并返回原始字节。"""
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": (
                "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
                "AppleWebKit/537.36 (KHTML, like Gecko) "
                "Chrome/123.0.0.0 Safari/537.36"
            )
        },
    )
    context = build_ssl_context(insecure=insecure)
    with urllib.request.urlopen(request, timeout=timeout, context=context) as response:
        return response.read()


def download_file(
    url: str,
    save_path: Path,
    overwrite: bool = False,
    timeout: int = 30,
    insecure: bool = False,
) -> Tuple[bool, str, int]:
    """下载单个文件。"""
    if save_path.exists() and not overwrite:
        return False, "skip_exists", save_path.stat().st_size

    content = fetch_url_bytes(url, timeout=timeout, insecure=insecure)
    ensure_dir(save_path.parent)
    save_path.write_bytes(content)
    return True, "downloaded", len(content)


def download_file_with_retry(
    *,
    url: str,
    save_path: Path,
    overwrite: bool,
    timeout: int,
    retry_count: int,
    retry_delay: float,
    backoff: float,
    insecure: bool,
) -> Tuple[bool, str, int, List[dict]]:
    """
    带重试的下载。

    返回：
        (是否成功, 最终状态, 字节数, 尝试日志列表)

    说明：
        - retry_count=4 表示最多尝试 4 次
        - 若文件已存在且 overwrite=False，则直接 skip，不再发起网络请求
    """
    if save_path.exists() and not overwrite:
        return True, "skip_exists", save_path.stat().st_size, [
            {
                "attempt": 0,
                "status": "skip_exists",
                "size_bytes": save_path.stat().st_size,
            }
        ]

    attempt_logs: List[dict] = []
    total_attempts = max(1, retry_count)

    for attempt in range(1, total_attempts + 1):
        try:
            written, status, size = download_file(
                url=url,
                save_path=save_path,
                overwrite=True if attempt > 1 else overwrite,
                timeout=timeout,
                insecure=insecure,
            )
            attempt_logs.append(
                {
                    "attempt": attempt,
                    "status": status,
                    "size_bytes": size,
                }
            )
            return True, status, size, attempt_logs
        except Exception as exc:  # noqa: BLE001
            attempt_logs.append(
                {
                    "attempt": attempt,
                    "status": "error",
                    "error": str(exc),
                }
            )
            if attempt < total_attempts:
                sleep_seconds = retry_delay * (backoff ** (attempt - 1))
                time.sleep(sleep_seconds)
            else:
                return False, "error", 0, attempt_logs

    return False, "error", 0, attempt_logs


def convert_tab_txt_to_csv(
    txt_path: Path,
    csv_path: Path,
    overwrite: bool = False,
) -> Tuple[bool, str, int]:
    """
    将以 Tab 为主分隔的 txt 表格转换为 csv。

    参数：
        txt_path: 原始 txt 文件路径
        csv_path: 输出 csv 文件路径
        overwrite: 是否覆盖已存在 csv

    返回：
        (是否成功, 状态, 写入总行数)

    状态含义：
        - written: 成功写入
        - skip_exists: csv 已存在且未要求覆盖
        - missing_txt: txt 文件不存在
    """
    if not txt_path.exists():
        return False, "missing_txt", 0

    if csv_path.exists() and not overwrite:
        return True, "skip_exists", 0

    ensure_dir(csv_path.parent)

    row_count = 0
    with txt_path.open("r", encoding="utf-8", errors="replace", newline="") as fin, \
         csv_path.open("w", encoding="utf-8-sig", newline="") as fout:
        reader = csv.reader(fin, delimiter="\t")
        writer = csv.writer(fout)

        for raw_row in reader:
            if not raw_row:
                continue

            cleaned = [cell.strip() for cell in raw_row]

            if all(cell == "" for cell in cleaned):
                continue

            writer.writerow(cleaned)
            row_count += 1

    return True, "written", row_count


def read_csv_rows(path: Path) -> List[dict]:
    """读取 CSV 为字典列表；若文件不存在则返回空列表。"""
    if not path.exists():
        return []

    with path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        return list(reader)


def write_csv_rows(path: Path, fieldnames: List[str], rows: List[dict]) -> None:
    """覆盖写入 CSV。"""
    ensure_dir(path.parent)
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def append_unique_source_manifest_row(csv_path: Path, row: dict) -> str:
    """向 source_manifest.csv 追加一条来源记录；以 source_id 去重。"""
    fieldnames = [
        "source_id",
        "task",
        "source_type",
        "source_name",
        "source_url_or_reference",
        "evidence_level",
        "access_date",
        "raw_storage_path",
        "notes",
    ]

    rows = read_csv_rows(csv_path)
    existing_ids = {r.get("source_id", "").strip() for r in rows}

    if row["source_id"] in existing_ids:
        return "skip_exists"

    rows.append(row)
    write_csv_rows(csv_path, fieldnames, rows)
    return "appended"


def _extract_existing_file_index(rows: Iterable[dict]) -> int:
    """从 file_catalog 现有 file_id 中提取最大编号。"""
    max_idx = 0
    prefix = "ace_db_ahtpdb_file_"

    for row in rows:
        file_id = row.get("file_id", "").strip()
        if not file_id.startswith(prefix):
            continue
        suffix = file_id.replace(prefix, "")
        if suffix.isdigit():
            max_idx = max(max_idx, int(suffix))
    return max_idx


def append_unique_file_catalog_rows(csv_path: Path, new_rows: List[dict]) -> Tuple[int, int]:
    """向 file_catalog.csv 追加若干文件记录；以 relative_path 去重。"""
    fieldnames = [
        "file_id",
        "source_id",
        "task",
        "relative_path",
        "file_type",
        "content_role",
        "created_or_downloaded_at",
        "generated_by",
        "version_tag",
        "notes",
    ]

    rows = read_csv_rows(csv_path)
    existing_paths = {r.get("relative_path", "").strip().replace("\\", "/") for r in rows}
    next_idx = _extract_existing_file_index(rows)

    appended = 0
    skipped = 0

    for row in new_rows:
        rel_path = row["relative_path"].strip().replace("\\", "/")
        if rel_path in existing_paths:
            skipped += 1
            continue

        next_idx += 1
        row["file_id"] = f"ace_db_ahtpdb_file_{next_idx:03d}"
        rows.append(row)
        existing_paths.add(rel_path)
        appended += 1

    write_csv_rows(csv_path, fieldnames, rows)
    return appended, skipped


def relative_to_root(root_dir: Path, target_path: Path) -> str:
    """将绝对路径转换为相对项目根目录路径，统一使用 / 分隔。"""
    return str(target_path.relative_to(root_dir)).replace("\\", "/")


def make_file_row(
    *,
    root_dir: Path,
    file_path: Path,
    file_type: str,
    content_role: str,
    created_or_downloaded_at: str,
    generated_by: str,
    notes: str,
) -> dict:
    """构造 file_catalog 单条记录。"""
    return {
        "file_id": "",
        "source_id": AHTPDB_SOURCE_ID,
        "task": AHTPDB_TASK,
        "relative_path": relative_to_root(root_dir, file_path),
        "file_type": file_type,
        "content_role": content_role,
        "created_or_downloaded_at": created_or_downloaded_at,
        "generated_by": generated_by,
        "version_tag": "v1",
        "notes": notes,
    }


# =====================================================================
# 4. AHTPDB 归档主逻辑
# =====================================================================
def run_ahtpdb_raw_archive(
    root_dir: Path,
    overwrite: bool = False,
    timeout: int = 30,
    retry_count: int = 4,
    retry_delay: float = 3.0,
    backoff: float = 1.8,
    insecure: bool = False,
    skip_pages: bool = False,
    skip_datasets: bool = False,
    page_key: str = "all",
    dataset_key: str = "all",
) -> None:
    """执行 AHTPDB 原始归档流程（重试增强版）。"""
    now = datetime.now()
    now_str = now.strftime("%Y-%m-%d %H:%M:%S")
    access_date = now.strftime("%Y-%m-%d")
    timestamp_tag = now.strftime("%Y%m%d_%H%M%S")

    db_dir = root_dir / "DB"
    meta_dir = db_dir / "meta"
    ahtpdb_dir = db_dir / "raw" / "ace" / "databases" / "AHTPDB"
    page_exports_dir = ahtpdb_dir / "page_exports"
    raw_tables_dir = ahtpdb_dir / "raw_tables"
    raw_tables_csv_dir = ahtpdb_dir / "raw_tables_csv"
    screenshots_dir = ahtpdb_dir / "screenshots"
    retrieval_logs_dir = ahtpdb_dir / "retrieval_logs"

    for path in [
        db_dir,
        meta_dir,
        ahtpdb_dir,
        page_exports_dir,
        raw_tables_dir,
        raw_tables_csv_dir,
        screenshots_dir,
        retrieval_logs_dir,
    ]:
        ensure_dir(path)

    print("=" * 80)
    print("AHTPDB 原始归档开始（重试增强版）")
    print(f"项目根目录: {root_dir}")
    print(f"AHTPDB 原始目录: {ahtpdb_dir}")
    print(f"timeout={timeout}s, retry_count={retry_count}, retry_delay={retry_delay}s, backoff={backoff}")
    print(f"skip_pages={skip_pages}, skip_datasets={skip_datasets}, page_key={page_key}, dataset_key={dataset_key}")
    print("=" * 80)

    log_records: List[dict] = []
    file_rows_to_add: List[dict] = []
    summary = {
        "page_success": 0,
        "page_skip": 0,
        "page_error": 0,
        "dataset_success": 0,
        "dataset_skip": 0,
        "dataset_error": 0,
        "csv_success": 0,
        "csv_skip": 0,
        "csv_error": 0,
    }

    # 1) README 与来源说明
    readme_path = ahtpdb_dir / "README.md"
    source_note_path = ahtpdb_dir / "source_note.txt"

    for file_path, content, content_role, generated_by, notes in [
        (
            readme_path,
            build_ahtpdb_readme(now_str),
            "source_readme",
            "script:init_ahtpdb_raw_archive.py",
            "AHTPDB 原始归档目录说明",
        ),
        (
            source_note_path,
            build_source_note(now_str),
            "retrieval_note",
            "script:init_ahtpdb_raw_archive.py",
            "AHTPDB 来源说明、URL 与获取口径说明",
        ),
    ]:
        written, status = write_text(file_path, content, overwrite=overwrite)
        print(f"[{status}] {file_path}")
        log_records.append(
            {
                "step": "write_text",
                "path": str(file_path),
                "status": status,
                "timestamp": now_str,
            }
        )
        file_rows_to_add.append(
            make_file_row(
                root_dir=root_dir,
                file_path=file_path,
                file_type=file_path.suffix.replace(".", "") or "txt",
                content_role=content_role,
                created_or_downloaded_at=now_str,
                generated_by=generated_by,
                notes=notes,
            )
        )

    # 2) 页面下载
    if not skip_pages:
        selected_page_items = AHTPDB_PAGE_URLS.items() if page_key == "all" else [(page_key, AHTPDB_PAGE_URLS[page_key])]
        for cur_page_key, url in selected_page_items:
            filename = AHTPDB_PAGE_EXPORT_FILENAMES[cur_page_key]
            save_path = page_exports_dir / filename
            ok, status, size, attempt_logs = download_file_with_retry(
                url=url,
                save_path=save_path,
                overwrite=overwrite,
                timeout=timeout,
                retry_count=retry_count,
                retry_delay=retry_delay,
                backoff=backoff,
                insecure=insecure,
            )
            if ok and status != "skip_exists":
                summary["page_success"] += 1
                print(f"[{status}] {save_path} ({size} bytes)")
            elif status == "skip_exists":
                summary["page_skip"] += 1
                print(f"[{status}] {save_path} ({size} bytes)")
            else:
                summary["page_error"] += 1
                last_error = attempt_logs[-1].get("error", "unknown error") if attempt_logs else "unknown error"
                print(f"[error] 页面下载失败: {url}\n        原因: {last_error}")

            log_records.append(
                {
                    "step": "download_page",
                    "page_key": cur_page_key,
                    "url": url,
                    "path": str(save_path),
                    "status": status,
                    "size_bytes": size,
                    "attempt_logs": attempt_logs,
                    "timestamp": now_str,
                }
            )

            if ok:
                file_rows_to_add.append(
                    make_file_row(
                        root_dir=root_dir,
                        file_path=save_path,
                        file_type="html",
                        content_role=f"page_export_{cur_page_key}",
                        created_or_downloaded_at=now_str,
                        generated_by="script:urllib_request_retry",
                        notes=f"AHTPDB 官方页面导出，来源 URL: {url}",
                    )
                )

    # 3) 数据集下载 + 自动生成 csv
    if not skip_datasets:
        selected_dataset_items = AHTPDB_DATASET_URLS.items() if dataset_key == "all" else [(dataset_key, AHTPDB_DATASET_URLS[dataset_key])]
        for cur_dataset_key, url in selected_dataset_items:
            filename = AHTPDB_DATASET_FILENAMES[cur_dataset_key]
            save_path = raw_tables_dir / filename

            ok, status, size, attempt_logs = download_file_with_retry(
                url=url,
                save_path=save_path,
                overwrite=overwrite,
                timeout=timeout,
                retry_count=retry_count,
                retry_delay=retry_delay,
                backoff=backoff,
                insecure=insecure,
            )

            if ok and status != "skip_exists":
                summary["dataset_success"] += 1
                print(f"[{status}] {save_path} ({size} bytes)")
            elif status == "skip_exists":
                summary["dataset_skip"] += 1
                print(f"[{status}] {save_path} ({size} bytes)")
            else:
                summary["dataset_error"] += 1
                last_error = attempt_logs[-1].get("error", "unknown error") if attempt_logs else "unknown error"
                print(f"[error] 数据集下载失败: {url}\n        原因: {last_error}")

            log_records.append(
                {
                    "step": "download_dataset",
                    "dataset_key": cur_dataset_key,
                    "url": url,
                    "path": str(save_path),
                    "status": status,
                    "size_bytes": size,
                    "attempt_logs": attempt_logs,
                    "timestamp": now_str,
                }
            )

            if ok:
                # 记录 txt 文件
                file_rows_to_add.append(
                    make_file_row(
                        root_dir=root_dir,
                        file_path=save_path,
                        file_type="txt",
                        content_role="raw_table",
                        created_or_downloaded_at=now_str,
                        generated_by="script:urllib_request_retry",
                        notes=f"AHTPDB 官方下载文本（dataset_key={cur_dataset_key}），来源 URL: {url}",
                    )
                )

                # ===== 新增：txt 下载成功/已存在后，自动生成 csv =====
                csv_save_path = raw_tables_csv_dir / (Path(filename).stem + ".csv")

                try:
                    csv_ok, csv_status, csv_row_count = convert_tab_txt_to_csv(
                        txt_path=save_path,
                        csv_path=csv_save_path,
                        overwrite=overwrite,
                    )

                    if csv_ok and csv_status == "written":
                        summary["csv_success"] += 1
                        print(f"[csv_written] {csv_save_path} | rows={csv_row_count}")
                    elif csv_ok and csv_status == "skip_exists":
                        summary["csv_skip"] += 1
                        print(f"[csv_skip_exists] {csv_save_path}")
                    else:
                        summary["csv_error"] += 1
                        print(f"[csv_error] {csv_save_path} | status={csv_status}")

                    log_records.append(
                        {
                            "step": "convert_dataset_txt_to_csv",
                            "dataset_key": cur_dataset_key,
                            "source_txt_path": str(save_path),
                            "csv_path": str(csv_save_path),
                            "status": csv_status,
                            "row_count": csv_row_count,
                            "timestamp": now_str,
                        }
                    )

                    if csv_ok:
                        file_rows_to_add.append(
                            make_file_row(
                                root_dir=root_dir,
                                file_path=csv_save_path,
                                file_type="csv",
                                content_role="raw_table_csv",
                                created_or_downloaded_at=now_str,
                                generated_by="script:txt_to_csv_after_download",
                                notes=(
                                    f"AHTPDB 原始 txt 自动转换生成的 csv（dataset_key={cur_dataset_key}），"
                                    f"对应原文件: {save_path.name}"
                                ),
                            )
                        )

                except Exception as exc:  # noqa: BLE001
                    summary["csv_error"] += 1
                    print(f"[csv_error] {csv_save_path} | reason={exc}")
                    log_records.append(
                        {
                            "step": "convert_dataset_txt_to_csv",
                            "dataset_key": cur_dataset_key,
                            "source_txt_path": str(save_path),
                            "csv_path": str(csv_save_path),
                            "status": "error",
                            "error": str(exc),
                            "timestamp": now_str,
                        }
                    )

    # 4) 检索日志
    retrieval_log_path = retrieval_logs_dir / f"ahtpdb_retrieval_log_{timestamp_tag}.json"
    retrieval_summary = {
        "script": "init_ahtpdb_raw_archive.py",
        "mode": "retry_enhanced",
        "source_id": AHTPDB_SOURCE_ID,
        "task": AHTPDB_TASK,
        "source_name": AHTPDB_SOURCE_NAME,
        "run_at": now_str,
        "project_root": str(root_dir),
        "overwrite": overwrite,
        "timeout": timeout,
        "retry_count": retry_count,
        "retry_delay": retry_delay,
        "backoff": backoff,
        "insecure": insecure,
        "skip_pages": skip_pages,
        "skip_datasets": skip_datasets,
        "page_key": page_key,
        "dataset_key": dataset_key,
        "summary": summary,
        "log_records": log_records,
    }
    retrieval_log_path.write_text(json.dumps(retrieval_summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"[written] {retrieval_log_path}")

    file_rows_to_add.append(
        make_file_row(
            root_dir=root_dir,
            file_path=retrieval_log_path,
            file_type="json",
            content_role="retrieval_log",
            created_or_downloaded_at=now_str,
            generated_by="script:init_ahtpdb_raw_archive.py",
            notes="AHTPDB 试点归档运行日志（重试增强版）",
        )
    )

    # 5) 更新 source_manifest.csv
    source_manifest_path = meta_dir / "source_manifest.csv"
    source_manifest_status = append_unique_source_manifest_row(
        source_manifest_path,
        {
            "source_id": AHTPDB_SOURCE_ID,
            "task": AHTPDB_TASK,
            "source_type": AHTPDB_SOURCE_TYPE,
            "source_name": AHTPDB_SOURCE_NAME,
            "source_url_or_reference": AHTPDB_PAGE_URLS["home"],
            "evidence_level": AHTPDB_EVIDENCE_LEVEL,
            "access_date": access_date,
            "raw_storage_path": relative_to_root(root_dir, ahtpdb_dir),
            "notes": (
                "AHTPDB 试点归档来源；当前阶段完成原始页面、官方下载 txt 与自动转换 csv 的归档，"
                "尚未进入标准化与工作集构建。"
            ),
        },
    )
    print(f"[{source_manifest_status}] {source_manifest_path}")

    # 6) 更新 file_catalog.csv
    file_catalog_path = meta_dir / "file_catalog.csv"
    appended_count, skipped_count = append_unique_file_catalog_rows(file_catalog_path, file_rows_to_add)
    print(f"[file_catalog] appended={appended_count}, skipped={skipped_count}: {file_catalog_path}")

    print("=" * 80)
    print("AHTPDB 原始归档完成（重试增强版）")
    print("运行摘要：")
    print(f"- 页面：success={summary['page_success']}, skip={summary['page_skip']}, error={summary['page_error']}")
    print(f"- 数据集：success={summary['dataset_success']}, skip={summary['dataset_skip']}, error={summary['dataset_error']}")
    print(f"- CSV：success={summary['csv_success']}, skip={summary['csv_skip']}, error={summary['csv_error']}")
    print("下一步建议：")
    print("1. 检查 raw_tables 下是否保留原始 txt")
    print("2. 检查 raw_tables_csv 下是否已自动生成对应 csv")
    print("3. 再检查 file_catalog.csv 中是否出现 raw_table_csv 类型记录")
    print("=" * 80)


# =====================================================================
# 5. 命令行入口
# =====================================================================
def parse_args(argv: List[str]) -> argparse.Namespace:
    """解析命令行参数。"""
    parser = argparse.ArgumentParser(description="PepDB - AHTPDB 原始归档脚本（重试增强版）")
    parser.add_argument(
        "--root-dir",
        default=r"E:\MYS\PepDB",
        help="PepDB 项目根目录，默认 E:\\MYS\\PepDB",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="若提供该参数，则覆盖已存在文件；默认不覆盖。",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=30,
        help="单次网络请求超时时间（秒），默认 30。",
    )
    parser.add_argument(
        "--retry-count",
        type=int,
        default=4,
        help="单个文件最大尝试次数，默认 4。",
    )
    parser.add_argument(
        "--retry-delay",
        type=float,
        default=3.0,
        help="首次重试等待秒数，默认 3.0。",
    )
    parser.add_argument(
        "--backoff",
        type=float,
        default=1.8,
        help="重试等待倍数，默认 1.8。",
    )
    parser.add_argument(
        "--insecure",
        action="store_true",
        help="若提供该参数，则使用不校验证书的 SSL 上下文。仅在网络兼容问题下临时使用。",
    )
    parser.add_argument(
        "--skip-pages",
        action="store_true",
        help="跳过页面下载，只处理说明文件与/或数据集。",
    )
    parser.add_argument(
        "--skip-datasets",
        action="store_true",
        help="跳过数据集下载，只处理说明文件与/或页面。",
    )
    parser.add_argument(
        "--page-key",
        choices=["all", "home", "information", "download", "guide", "publications"],
        default="all",
        help="只处理某个页面；默认 all。",
    )
    parser.add_argument(
        "--dataset-key",
        choices=["all", "small", "long", "ic50"],
        default="all",
        help="只处理某个数据集；默认 all。",
    )
    return parser.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> None:
    """主函数。"""
    args = parse_args(argv or sys.argv[1:])
    root_dir = Path(args.root_dir)
    run_ahtpdb_raw_archive(
        root_dir=root_dir,
        overwrite=args.overwrite,
        timeout=args.timeout,
        retry_count=args.retry_count,
        retry_delay=args.retry_delay,
        backoff=args.backoff,
        insecure=args.insecure,
        skip_pages=args.skip_pages,
        skip_datasets=args.skip_datasets,
        page_key=args.page_key,
        dataset_key=args.dataset_key,
    )


if __name__ == "__main__":
    main()