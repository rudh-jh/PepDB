#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MBPDB ACE core workset 构建脚本（纯标准库版）
==========================================

输入
----
DB/standardized/ace/mbpdb/mbpdb_ace_master_clean.csv

输出
----
DB/worksets/ace/mbpdb/
├── mbpdb_ace_core_record_level.csv
├── mbpdb_ace_core_sequence_level.csv
├── mbpdb_ace_short_2_3_record_level.csv
├── mbpdb_ace_short_2_3_sequence_level.csv
├── mbpdb_ace_excluded_missing_ic50_records.csv
├── mbpdb_ace_excluded_noncanonical_records.csv
└── mbpdb_ace_workset_summary.csv

说明
----
- 这是 workset 层，不修改 standardized 原表
- MBPDB 当前已经是 ACE-inhibitory 子集，因此不再做 function 主过滤
- 核心保留规则相对简单：
  1) ic50_parse_status ∈ {exact_molar, molar_threshold}
  2) ic50_uM 非空
  3) sequence 非空且为 canonical
- 再额外导出二肽/三肽子集
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional


INPUT_FILE = "mbpdb_ace_master_clean.csv"

CORE_RECORD_FILE = "mbpdb_ace_core_record_level.csv"
CORE_SEQUENCE_FILE = "mbpdb_ace_core_sequence_level.csv"
SHORT_RECORD_FILE = "mbpdb_ace_short_2_3_record_level.csv"
SHORT_SEQUENCE_FILE = "mbpdb_ace_short_2_3_sequence_level.csv"
EXCLUDED_MISSING_IC50_FILE = "mbpdb_ace_excluded_missing_ic50_records.csv"
EXCLUDED_NONCANONICAL_FILE = "mbpdb_ace_excluded_noncanonical_records.csv"
SUMMARY_FILE = "mbpdb_ace_workset_summary.csv"

RESOLVED_STATUSES = {"exact_molar", "molar_threshold"}
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
    return normalize_spaces(str(value)).lower() in MISSING_TOKENS


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


def write_csv(path: Path, rows: List[Dict[str, object]], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def split_notes(notes: object) -> List[str]:
    text = clean_text(notes)
    if text is None:
        return []
    return [x.strip() for x in text.split(";") if x.strip()]


def has_note(row: Dict[str, str], token: str) -> bool:
    return token in split_notes(row.get("notes"))


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


def year_int_from_text(value: object) -> Optional[int]:
    text = clean_text(value)
    if text is None:
        return None
    m = __import__("re").findall(r"\b(19\d{2}|20\d{2})\b", text)
    if not m:
        return None
    years = [int(x) for x in m]
    return min(years) if years else None


# =========================
# 核心判定
# =========================
def is_resolved_ic50_row(row: Dict[str, str]) -> bool:
    status = clean_text(row.get("ic50_parse_status"))
    ic50_uM = to_float_safe(row.get("ic50_uM"))
    return (status in RESOLVED_STATUSES) and (ic50_uM is not None)


def is_noncanonical_row(row: Dict[str, str]) -> bool:
    seq = clean_text(row.get("sequence")) or ""
    if not seq:
        return True
    if has_note(row, "non_canonical_or_contains_unknown_residue"):
        return True
    return False


def build_exclusion_reasons(row: Dict[str, str]) -> List[str]:
    reasons = []
    if not is_resolved_ic50_row(row):
        reasons.append("missing_or_unresolved_ic50")
    if is_noncanonical_row(row):
        reasons.append("noncanonical_or_missing_sequence")
    return reasons


def is_core_record(row: Dict[str, str]) -> bool:
    if not is_resolved_ic50_row(row):
        return False
    if is_noncanonical_row(row):
        return False
    return True


def add_workset_fields(row: Dict[str, str], keep_flag: int, exclusion_reasons: List[str]) -> Dict[str, object]:
    out = dict(row)
    out["workset_keep_flag"] = keep_flag
    out["workset_exclusion_reasons"] = "; ".join(exclusion_reasons) if exclusion_reasons else ""
    return out


# =========================
# record-level
# =========================
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
        if mode == "missing_ic50":
            match = "missing_or_unresolved_ic50" in reasons
        elif mode == "noncanonical":
            match = "noncanonical_or_missing_sequence" in reasons
        else:
            raise ValueError(f"未知 mode: {mode}")

        if match:
            out.append(add_workset_fields(row, 0, reasons))
    return out


# =========================
# sequence-level
# =========================
SEQUENCE_LEVEL_COLUMNS = [
    "sequence",
    "peptide_length",
    "record_count",
    "record_id_list",
    "source_record_id_list",
    "species_set",
    "protein_id_set",
    "protein_description_set",
    "doi_set",
    "ic50_relation_set",
    "ic50_parse_status_set",
    "ic50_uM_min",
    "ic50_uM_max",
    "ic50_uM_mean",
    "ic50_uM_median",
    "ic50_spread_ratio_max_min",
    "bibliographic_year_min",
    "bibliographic_year_max",
    "title_set",
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

            y = year_int_from_text(m.get("title_raw"))
            if y is None:
                y = year_int_from_text(m.get("doi_raw"))
            if y is not None:
                years.append(y)

        ic50_min = safe_min(ic50_values)
        ic50_max = safe_max(ic50_values)
        spread_ratio = None
        if ic50_min is not None and ic50_max is not None and ic50_min > 0:
            spread_ratio = ic50_max / ic50_min

        row = {
            "sequence": seq,
            "peptide_length": to_int_safe(members[0].get("peptide_length")),
            "record_count": len(members),
            "record_id_list": sorted_unique_join([m.get("record_id") for m in members]),
            "source_record_id_list": sorted_unique_join([m.get("source_record_id") for m in members]),
            "species_set": sorted_unique_join([m.get("species_raw") for m in members]),
            "protein_id_set": sorted_unique_join([m.get("protein_id_raw") for m in members]),
            "protein_description_set": sorted_unique_join([m.get("protein_description_raw") for m in members]),
            "doi_set": sorted_unique_join([m.get("doi_raw") for m in members]),
            "ic50_relation_set": sorted_unique_join([m.get("ic50_relation") for m in members]),
            "ic50_parse_status_set": sorted_unique_join([m.get("ic50_parse_status") for m in members]),
            "ic50_uM_min": ic50_min,
            "ic50_uM_max": ic50_max,
            "ic50_uM_mean": safe_mean(ic50_values),
            "ic50_uM_median": safe_median(ic50_values),
            "ic50_spread_ratio_max_min": spread_ratio,
            "bibliographic_year_min": min(years) if years else None,
            "bibliographic_year_max": max(years) if years else None,
            "title_set": sorted_unique_join([m.get("title_raw") for m in members]),
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
    missing_ic50_rows: List[Dict[str, object]],
    noncanonical_rows: List[Dict[str, object]],
    core_seq_rows: List[Dict[str, object]],
    short_seq_rows: List[Dict[str, object]],
) -> List[Dict[str, object]]:
    summary = []

    total_rows = len(all_rows)
    resolved_rows = sum(1 for r in all_rows if is_resolved_ic50_row(r))
    core_count = len(core_rows)
    short_count = len(short_rows)

    summary.extend([
        {"metric": "input_rows_total", "value": total_rows},
        {"metric": "resolved_ic50_rows", "value": resolved_rows},
        {"metric": "excluded_missing_or_unresolved_ic50_rows", "value": len(missing_ic50_rows)},
        {"metric": "excluded_noncanonical_or_missing_sequence_rows", "value": len(noncanonical_rows)},
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

    species_counter = Counter(
        clean_text(r.get("species_raw")) or ""
        for r in core_rows
        if clean_text(r.get("species_raw")) is not None
    )
    for k, v in species_counter.items():
        summary.append({"metric": f"core_species::{k}", "value": v})

    relation_counter = Counter(
        clean_text(r.get("ic50_relation")) or ""
        for r in core_rows
    )
    for k, v in relation_counter.items():
        summary.append({"metric": f"core_ic50_relation::{k}", "value": v})

    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description="构建 MBPDB ACE core workset（纯标准库版）")
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
        help="输入 standardized 主表路径；默认 DB/standardized/ace/mbpdb/mbpdb_ace_master_clean.csv",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="输出 worksets 目录；默认 DB/worksets/ace/mbpdb",
    )
    args = parser.parse_args()

    script_file = Path(__file__)
    project_root = Path(args.project_root) if args.project_root else infer_project_root(script_file)

    input_path = Path(args.input_path) if args.input_path else (
        project_root / "DB" / "standardized" / "ace" / "mbpdb" / INPUT_FILE
    )
    output_dir = Path(args.output_dir) if args.output_dir else (
        project_root / "DB" / "worksets" / "ace" / "mbpdb"
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

    missing_ic50_rows = excluded_rows(all_rows, "missing_ic50")
    noncanonical_rows = excluded_rows(all_rows, "noncanonical")

    write_csv(output_dir / CORE_RECORD_FILE, core_rows, record_columns)
    write_csv(output_dir / CORE_SEQUENCE_FILE, core_seq_rows, SEQUENCE_LEVEL_COLUMNS)
    write_csv(output_dir / SHORT_RECORD_FILE, short_rows, record_columns)
    write_csv(output_dir / SHORT_SEQUENCE_FILE, short_seq_rows, SEQUENCE_LEVEL_COLUMNS)
    write_csv(output_dir / EXCLUDED_MISSING_IC50_FILE, missing_ic50_rows, record_columns)
    write_csv(output_dir / EXCLUDED_NONCANONICAL_FILE, noncanonical_rows, record_columns)

    summary_rows = build_summary_rows(
        all_rows=all_rows,
        core_rows=core_rows,
        short_rows=short_rows,
        missing_ic50_rows=missing_ic50_rows,
        noncanonical_rows=noncanonical_rows,
        core_seq_rows=core_seq_rows,
        short_seq_rows=short_seq_rows,
    )
    write_csv(output_dir / SUMMARY_FILE, summary_rows, ["metric", "value"])

    print("=" * 80)
    print("MBPDB ACE core workset 构建完成")
    print(f"项目根目录：{project_root}")
    print(f"输入表：{input_path}")
    print(f"输出目录：{output_dir}")
    print(f"输入记录数：{len(all_rows)}")
    print(f"核心 record-level：{len(core_rows)}")
    print(f"核心 sequence-level：{len(core_seq_rows)}")
    print(f"二肽/三肽 record-level：{len(short_rows)}")
    print(f"二肽/三肽 sequence-level：{len(short_seq_rows)}")
    print(f"排除 missing/unresolved IC50：{len(missing_ic50_rows)}")
    print(f"排除 noncanonical/missing sequence：{len(noncanonical_rows)}")
    print("=" * 80)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())