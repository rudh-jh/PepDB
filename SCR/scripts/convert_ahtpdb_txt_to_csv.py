# -*- coding: utf-8 -*-
"""
脚本名称：
    convert_ahtpdb_txt_to_csv.py

脚本作用：
    将 AHTPDB 当前已下载的 3 个原始 txt 文件转换为 csv 文件。
    这些 txt 文件本质上是以 Tab 为主分隔的表格文本，因此可直接解析后另存为 CSV。

适用场景：
    1. 你已经完成 AHTPDB 原始下载
    2. 不想重跑下载脚本
    3. 只想把 raw_tables 下现有 txt 转成 csv

默认输入目录：
    E:\MYS\PepDB\DB\raw\ace\databases\AHTPDB\raw_tables

默认输出目录：
    E:\MYS\PepDB\DB\raw\ace\databases\AHTPDB\raw_tables_csv

注意：
    1. 本脚本不会删除原始 txt
    2. 本脚本不会修改原始 txt
    3. 本脚本只做格式转换，不做单位统一、不做去重、不做字段标准化
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import List


DEFAULT_ROOT = Path(r"E:\MYS\PepDB")

TXT_FILES = [
    "ahtpdb_small_peptides_2_5_residues.txt",
    "ahtpdb_long_peptides_6_plus_residues.txt",
    "ahtpdb_peptides_with_ic50.txt",
]


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def parse_tab_text_file(txt_path: Path) -> List[List[str]]:
    """
    读取 AHTPDB txt 文件并解析为二维列表。

    当前 AHTPDB 的 txt 文件可以按 Tab 读取。
    同时这里会对每个单元格做 strip()，尽量去掉首尾空白。
    """
    rows: List[List[str]] = []
    with txt_path.open("r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        for raw_row in reader:
            if not raw_row:
                continue
            cleaned = [cell.strip() for cell in raw_row]
            if all(cell == "" for cell in cleaned):
                continue
            rows.append(cleaned)
    return rows


def write_csv(csv_path: Path, rows: List[List[str]]) -> None:
    ensure_dir(csv_path.parent)
    with csv_path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="将 AHTPDB txt 转为 csv")
    parser.add_argument(
        "--root-dir",
        default=str(DEFAULT_ROOT),
        help="PepDB 项目根目录，默认 E:\\MYS\\PepDB",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="若提供该参数，则允许覆盖已存在 csv 文件",
    )
    args = parser.parse_args()

    root_dir = Path(args.root_dir)
    txt_dir = root_dir / "DB" / "raw" / "ace" / "databases" / "AHTPDB" / "raw_tables"
    csv_dir = root_dir / "DB" / "raw" / "ace" / "databases" / "AHTPDB" / "raw_tables_csv"

    print("=" * 70)
    print("AHTPDB txt -> csv 转换开始")
    print(f"输入目录: {txt_dir}")
    print(f"输出目录: {csv_dir}")
    print("=" * 70)

    ensure_dir(csv_dir)

    for txt_name in TXT_FILES:
        txt_path = txt_dir / txt_name
        csv_name = txt_path.stem + ".csv"
        csv_path = csv_dir / csv_name

        if not txt_path.exists():
            print(f"[missing] {txt_path}")
            continue

        if csv_path.exists() and not args.overwrite:
            print(f"[skip_exists] {csv_path}")
            continue

        rows = parse_tab_text_file(txt_path)
        if not rows:
            print(f"[warning] 未解析到有效行: {txt_path}")
            continue

        write_csv(csv_path, rows)
        print(f"[written] {csv_path} | rows={len(rows)}")

    print("=" * 70)
    print("AHTPDB txt -> csv 转换完成")
    print("=" * 70)


if __name__ == "__main__":
    main()
