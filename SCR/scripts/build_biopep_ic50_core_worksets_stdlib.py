#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BIOPEP-UWM ACE IC50 核心工作集构建脚本（纯标准库版）
=================================================

输入
----
DB/standardized/ace/biopep_uwm/biopep_uwm_master_enriched.csv

输出
----
DB/worksets/ace/biopep_uwm/
├── biopep_uwm_ic50_core_record_level.csv
├── biopep_uwm_ic50_core_sequence_level.csv
├── biopep_uwm_ic50_short_2_3_record_level.csv
├── biopep_uwm_ic50_short_2_3_sequence_level.csv
├── biopep_uwm_excluded_ec50_records.csv
├── biopep_uwm_excluded_noncanonical_records.csv
├── biopep_uwm_excluded_length_conflict_records.csv
└── biopep_uwm_workset_summary.csv

筛选逻辑
--------
核心 IC50 record-level 保留条件：
1. ic50_parse_status ∈ {exact_molar, molar_threshold}
2. ic50_uM 非空
3. 不是 EC50 记录
4. notes 不含 non_canonical_or_contains_unknown_residue
5. length_consistency_flag != conflict

说明
----
- 这是 workset 层，不修改 standardized 原表
- 排除集单独导出，便于后续写数据说明
- sequence-level 按 sequence 聚合
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple


INPUT_FILE = "biopep_uwm_master_enriched.csv"

CORE_RECORD_FILE = "biopep_uwm_ic50_core_record_level.csv"
CORE_SEQUENCE_FILE = "biopep_uwm_ic50_core_sequence_level.csv"
SHORT_RECORD_FILE = "biopep_uwm_ic50_short_2_3_record_level.csv"
SHORT_SEQUENCE_FILE = "biopep_uwm_ic50_short_2_3_sequence_level.csv"

EXCLUDED_EC50_FILE = "biopep_uwm_excluded_ec50_records.csv"
EXCLUDED_NONCANONICAL_FILE = "biopep_uwm_excluded_noncanonical_records.csv"
EXCLUDED_LENGTH_CONFLICT_FILE = "biopep_uwm_excluded_length_conflict_records.csv"

SUMMARY_FILE = "biopep_uwm_workset_summary.csv"

MISSING_TOKENS = {
    "",
    "na",
    "n/a",
    "nan",
    "none",
    "null",
    "-",
    "--",
    "nd",
    "n.d.",
}

RESOLVED_STATUSES = {"exact_molar", "molar_threshold"}


# =========================
# 通用工具
# =========================
def infer_project_root(script_file: Path) -> Path:
    try:
        return script_file.resolve().parents[2]
    except IndexError:
        return script_file.resolve().parent


def strip_bom(text: str) -> str:
    return text.lstrip("\ufeff")


def normalize_spaces(text: str) -> str:
    return " ".join(text.replace("\r", " ").replace("\n", " ").split()).strip()


def is_missing(value: object) -> bool:
    if value is None:
        return True
    text = normalize_spaces(str(value)).lower()
    return text in MISSING_TOKENS


def clean_text(value: object) -> Optional[str]:
    if is_missing(value):
        return None
    return normalize_spaces(str(value))


def to_float_safe(value: object) -> Optional[float]:
    text = clean_text(value)
    if text is None:
        return None
    text = text.replace(",", "")
    try:
        x = float(text)
    except ValueError:
        return None
    if math.isnan(x) or math.isinf(x):
        return None
    return x


def to_int_safe(value: object) -> Optional[int]:
    x = to_float_safe(value)
    if x is None:
        return None
    if abs(x - round(x)) < 1e-9:
        return int(round(x))
    return None


def write_csv(path: Path, rows: List[Dict[str, object]], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def read_csv_rows(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"找不到输入文件：{path}")
    with open(path, "r", encoding="utf-8-sig", errors="replace", newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError(f"CSV 无表头：{path}")
        rows = []
        for row in reader:
            clean_row = {}
            for k, v in row.items():
                clean_row[strip_bom(k)] = v
            rows.append(clean_row)
    return rows


def split_notes(notes: object) -> List[str]:
    text = clean_text(notes)
    if text is None:
        return []
    return [x.strip() for x in text.split(";") if x.strip()]


def has_note(row: Dict[str, str], token: str) -> bool:
    return token in split_notes(row.get("notes"))


def contains_case_insensitive(text: object, keyword: str) -> bool:
    t = clean_text(text)
    if t is None:
        return False
    return keyword.lower() in t.lower()


def safe_min(values: List[float]) -> Optional[float]:
    return min(values) if values else None


def safe_max(values: List[float]) -> Optional[float]:
    return max(values) if values else None


def safe_mean(values: List[float]) -> Optional[float]:
    return statistics.mean(values) if values else None


def safe_median(values: List[float]) -> Optional[float]:
    return statistics.median(values) if values else None


def sorted_unique_join(values: List[object], sep: str = " | ") -> str:
    cleaned = []
    seen = set()
    for v in values:
        t = clean_text(v)
        if t is None:
            continue
        if t not in seen:
            seen.add(t)
            cleaned.append(t)
    cleaned.sort()
    return sep.join(cleaned)


def year_int(value: object) -> Optional[int]:
    text = clean_text(value)
    if text is None:
        return None
    digits = "".join(ch for ch in text if ch.isdigit())
    if len(digits) >= 4:
        y = int(digits[:4])
        if 1800 <= y <= 2100:
            return y
    return None


# =========================
# 记录判定
# =========================
def is_resolved_ic50_row(row: Dict[str, str]) -> bool:
    status = clean_text(row.get("ic50_parse_status"))
    ic50_uM = to_float_safe(row.get("ic50_uM"))
    return (status in RESOLVED_STATUSES) and (ic50_uM is not None)


def is_ec50_row(row: Dict[str, str]) -> bool:
    candidates = [
        row.get("measure_type_raw"),
        row.get("activity_measure_label_raw"),
        row.get("ic50_raw"),
        row.get("activity_name_raw"),
        row.get("activity_label_raw"),
    ]
    return any(contains_case_insensitive(x, "EC50") for x in candidates)


def is_noncanonical_row(row: Dict[str, str]) -> bool:
    seq = clean_text(row.get("sequence")) or ""
    if "{" in seq or "}" in seq:
        return True
    return has_note(row, "non_canonical_or_contains_unknown_residue")


def is_length_conflict_row(row: Dict[str, str]) -> bool:
    return clean_text(row.get("length_consistency_flag")) == "conflict"


def build_exclusion_reasons(row: Dict[str, str]) -> List[str]:
    reasons = []

    if not is_resolved_ic50_row(row):
        reasons.append("unresolved_ic50")

    if is_ec50_row(row):
        reasons.append("ec50_record")

    if is_noncanonical_row(row):
        reasons.append("noncanonical_or_modified_residue")

    if is_length_conflict_row(row):
        reasons.append("length_conflict")

    return reasons


def is_core_record(row: Dict[str, str]) -> bool:
    if not is_resolved_ic50_row(row):
        return False
    if is_ec50_row(row):
        return False
    if is_noncanonical_row(row):
        return False
    if is_length_conflict_row(row):
        return False
    return True


# =========================
# 导出行加工
# =========================
def add_workset_fields(row: Dict[str, str], keep_flag: int, exclusion_reasons: List[str]) -> Dict[str, object]:
    out = dict(row)
    out["workset_keep_flag"] = keep_flag
    out["workset_exclusion_reasons"] = "; ".join(exclusion_reasons) if exclusion_reasons else ""
    return out


def core_record_level_rows(rows: List[Dict[str, str]]) -> List[Dict[str, object]]:
    out = []
    for row in rows:
        reasons = build_exclusion_reasons(row)
        keep = 1 if is_core_record(row) else 0
        if keep == 1:
            out.append(add_workset_fields(row, keep, reasons))
    return out


def excluded_rows(rows: List[Dict[str, str]], mode: str) -> List[Dict[str, object]]:
    out = []
    for row in rows:
        reasons = build_exclusion_reasons(row)

        match = False
        if mode == "ec50":
            match = "ec50_record" in reasons
        elif mode == "noncanonical":
            match = "noncanonical_or_modified_residue" in reasons
        elif mode == "length_conflict":
            match = "length_conflict" in reasons
        else:
            raise ValueError(f"未知 mode: {mode}")

        if match:
            out.append(add_workset_fields(row, 0, reasons))
    return out


# =========================
# sequence-level 聚合
# =========================
SEQUENCE_LEVEL_COLUMNS = [
    "sequence",
    "peptide_length",
    "record_count",
    "record_id_list",
    "source_record_id_list",
    "source_name_set",
    "ic50_parse_status_set",
    "ic50_relation_set",
    "ic50_uM_min",
    "ic50_uM_max",
    "ic50_uM_mean",
    "ic50_uM_median",
    "ic50_spread_ratio_max_min",
    "bibliographic_year_min",
    "bibliographic_year_max",
    "detail_report_url_list",
    "notes_set",
]


def build_sequence_level(rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
    groups: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    for row in rows:
        seq = clean_text(row.get("sequence"))
        if seq is None:
            continue
        groups[seq].append(row)

    out = []
    for seq, members in sorted(groups.items(), key=lambda kv: (len(kv[0]), kv[0])):
        ic50_values = []
        years = []

        for m in members:
            x = to_float_safe(m.get("ic50_uM"))
            if x is not None:
                ic50_values.append(x)

            y = year_int(m.get("bibliographic_year_raw"))
            if y is not None:
                years.append(y)

        ic50_min = safe_min(ic50_values)
        ic50_max = safe_max(ic50_values)
        if ic50_min is not None and ic50_max is not None and ic50_min > 0:
            spread_ratio = ic50_max / ic50_min
        else:
            spread_ratio = None

        row = {
            "sequence": seq,
            "peptide_length": to_int_safe(members[0].get("peptide_length")),
            "record_count": len(members),
            "record_id_list": sorted_unique_join([m.get("record_id") for m in members]),
            "source_record_id_list": sorted_unique_join([m.get("source_record_id") for m in members]),
            "source_name_set": sorted_unique_join([m.get("source_name") for m in members]),
            "ic50_parse_status_set": sorted_unique_join([m.get("ic50_parse_status") for m in members]),
            "ic50_relation_set": sorted_unique_join([m.get("ic50_relation") for m in members]),
            "ic50_uM_min": ic50_min,
            "ic50_uM_max": ic50_max,
            "ic50_uM_mean": safe_mean(ic50_values),
            "ic50_uM_median": safe_median(ic50_values),
            "ic50_spread_ratio_max_min": spread_ratio,
            "bibliographic_year_min": min(years) if years else None,
            "bibliographic_year_max": max(years) if years else None,
            "detail_report_url_list": sorted_unique_join([m.get("detail_report_url") for m in members]),
            "notes_set": sorted_unique_join([m.get("notes") for m in members]),
        }
        out.append(row)

    out.sort(key=lambda r: (
        to_int_safe(r.get("peptide_length")) if r.get("peptide_length") is not None else 10**9,
        str(r.get("sequence") or ""),
    ))
    return out


# =========================
# summary
# =========================
def build_summary_rows(
    all_rows: List[Dict[str, str]],
    core_rows: List[Dict[str, object]],
    short_rows: List[Dict[str, object]],
    ec50_rows: List[Dict[str, object]],
    noncanonical_rows: List[Dict[str, object]],
    conflict_rows: List[Dict[str, object]],
    core_seq_rows: List[Dict[str, object]],
    short_seq_rows: List[Dict[str, object]],
) -> List[Dict[str, object]]:
    summary = []

    total_rows = len(all_rows)
    resolved_rows = sum(1 for r in all_rows if is_resolved_ic50_row(r))
    ec50_count = len(ec50_rows)
    noncanonical_count = len(noncanonical_rows)
    conflict_count = len(conflict_rows)
    core_count = len(core_rows)
    short_count = len(short_rows)

    summary.extend([
        {"metric": "input_rows_total", "value": total_rows},
        {"metric": "resolved_ic50_rows", "value": resolved_rows},
        {"metric": "excluded_ec50_rows", "value": ec50_count},
        {"metric": "excluded_noncanonical_rows", "value": noncanonical_count},
        {"metric": "excluded_length_conflict_rows", "value": conflict_count},
        {"metric": "core_record_level_rows", "value": core_count},
        {"metric": "core_sequence_level_rows", "value": len(core_seq_rows)},
        {"metric": "short_2_3_record_level_rows", "value": short_count},
        {"metric": "short_2_3_sequence_level_rows", "value": len(short_seq_rows)},
    ])

    length_counter = Counter(
        to_int_safe(r.get("peptide_length"))
        for r in core_rows
        if to_int_safe(r.get("peptide_length")) is not None
    )
    for k in sorted(length_counter):
        summary.append({"metric": f"core_length::{k}", "value": length_counter[k]})

    relation_counter = Counter(
        clean_text(r.get("ic50_relation")) or ""
        for r in core_rows
    )
    for k, v in relation_counter.items():
        summary.append({"metric": f"core_ic50_relation::{k}", "value": v})

    return summary


# =========================
# 主程序
# =========================
def main() -> int:
    parser = argparse.ArgumentParser(description="构建 BIOPEP-UWM ACE IC50 核心工作集（纯标准库版）")
    parser.add_argument(
        "--project-root",
        type=str,
        default=None,
        help="PepDB 项目根目录，例如 E:\\MYS\\PepDB",
    )
    parser.add_argument(
        "--input-path",
        type=str,
        default=None,
        help="输入 enriched 主表路径；默认 DB/standardized/ace/biopep_uwm/biopep_uwm_master_enriched.csv",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="输出 worksets 目录；默认 DB/worksets/ace/biopep_uwm",
    )
    args = parser.parse_args()

    script_file = Path(__file__)
    project_root = Path(args.project_root) if args.project_root else infer_project_root(script_file)

    input_path = Path(args.input_path) if args.input_path else (
        project_root / "DB" / "standardized" / "ace" / "biopep_uwm" / INPUT_FILE
    )
    output_dir = Path(args.output_dir) if args.output_dir else (
        project_root / "DB" / "worksets" / "ace" / "biopep_uwm"
    )

    all_rows = read_csv_rows(input_path)
    if not all_rows:
        raise ValueError(f"输入表为空：{input_path}")

    record_columns = list(all_rows[0].keys()) + ["workset_keep_flag", "workset_exclusion_reasons"]

    core_rows = core_record_level_rows(all_rows)
    core_seq_rows = build_sequence_level(core_rows)

    short_rows = [
        r for r in core_rows
        if to_int_safe(r.get("peptide_length")) in {2, 3}
    ]
    short_seq_rows = build_sequence_level(short_rows)

    ec50_rows = excluded_rows(all_rows, "ec50")
    noncanonical_rows = excluded_rows(all_rows, "noncanonical")
    conflict_rows = excluded_rows(all_rows, "length_conflict")

    write_csv(output_dir / CORE_RECORD_FILE, core_rows, record_columns)
    write_csv(output_dir / CORE_SEQUENCE_FILE, core_seq_rows, SEQUENCE_LEVEL_COLUMNS)

    write_csv(output_dir / SHORT_RECORD_FILE, short_rows, record_columns)
    write_csv(output_dir / SHORT_SEQUENCE_FILE, short_seq_rows, SEQUENCE_LEVEL_COLUMNS)

    write_csv(output_dir / EXCLUDED_EC50_FILE, ec50_rows, record_columns)
    write_csv(output_dir / EXCLUDED_NONCANONICAL_FILE, noncanonical_rows, record_columns)
    write_csv(output_dir / EXCLUDED_LENGTH_CONFLICT_FILE, conflict_rows, record_columns)

    summary_rows = build_summary_rows(
        all_rows=all_rows,
        core_rows=core_rows,
        short_rows=short_rows,
        ec50_rows=ec50_rows,
        noncanonical_rows=noncanonical_rows,
        conflict_rows=conflict_rows,
        core_seq_rows=core_seq_rows,
        short_seq_rows=short_seq_rows,
    )
    write_csv(output_dir / SUMMARY_FILE, summary_rows, ["metric", "value"])

    print("=" * 80)
    print("BIOPEP-UWM IC50 核心工作集构建完成")
    print(f"项目根目录：{project_root}")
    print(f"输入表：{input_path}")
    print(f"输出目录：{output_dir}")
    print(f"输入记录数：{len(all_rows)}")
    print(f"核心 record-level：{len(core_rows)}")
    print(f"核心 sequence-level：{len(core_seq_rows)}")
    print(f"二肽/三肽 record-level：{len(short_rows)}")
    print(f"二肽/三肽 sequence-level：{len(short_seq_rows)}")
    print(f"排除 EC50：{len(ec50_rows)}")
    print(f"排除 非标准残基：{len(noncanonical_rows)}")
    print(f"排除 长度冲突：{len(conflict_rows)}")
    print("=" * 80)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())