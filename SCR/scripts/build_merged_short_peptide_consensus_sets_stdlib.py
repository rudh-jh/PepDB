#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
构建 merged 短肽共识分层表（纯标准库版）
====================================

输入
----
DB/worksets/ace/merged/ace_short_2_3_merged_sequence_level_v0_1.csv

输出
----
DB/worksets/ace/merged/
├── ace_short_2_3_consensus_annotated_v0_1.csv
├── ace_short_2_3_consensus_tier_A.csv
├── ace_short_2_3_consensus_tier_B.csv
├── ace_short_2_3_singletons.csv
├── ace_short_2_3_conflicts.csv
├── ace_short_2_3_highconf_union.csv
├── ace_short_2_3_conflict_details.csv
└── ace_short_2_3_consensus_summary_v0_1.csv

分层逻辑
--------
Tier_A = cross_source_consistent
Tier_B = within_range_consistent
Tier_C = singleton
Tier_X = conflicting / unknown / 其它未覆盖状态

附加说明
--------
- 本脚本只做 sequence-level 的证据分层
- 不修改原 merged 输入表
- 会额外给出 highconf_union = Tier_A + Tier_B
- 会额外给出 conflict_details，便于后续人工排查冲突短肽
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional


INPUT_FILE = "ace_short_2_3_merged_sequence_level_v0_1.csv"

ANNOTATED_FILE = "ace_short_2_3_consensus_annotated_v0_1.csv"
TIER_A_FILE = "ace_short_2_3_consensus_tier_A.csv"
TIER_B_FILE = "ace_short_2_3_consensus_tier_B.csv"
TIER_C_FILE = "ace_short_2_3_singletons.csv"
TIER_X_FILE = "ace_short_2_3_conflicts.csv"
HIGHCONF_FILE = "ace_short_2_3_highconf_union.csv"
CONFLICT_DETAILS_FILE = "ace_short_2_3_conflict_details.csv"
SUMMARY_FILE = "ace_short_2_3_consensus_summary_v0_1.csv"

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

BASE_COLUMNS = [
    "sequence",
    "peptide_length",
    "source_count",
    "source_name_set",
    "record_count_total",
    "record_count_ahtpdb",
    "record_count_biopep",
    "merged_record_id_list",
    "source_record_id_list",
    "ic50_relation_set",
    "ic50_status_set",
    "ic50_uM_min",
    "ic50_uM_max",
    "ic50_uM_mean",
    "ic50_uM_median",
    "ic50_spread_ratio_max_min",
    "cross_source_overlap_flag",
    "stability_flag",
    "stability_note",
    "bibliographic_year_min",
    "bibliographic_year_max",
    "database_reference_raw_set",
]

ANNOTATED_EXTRA_COLUMNS = [
    "consensus_tier",
    "consensus_keep_flag",
    "consensus_rank",
    "consensus_note",
    "priority_score",
]

CONFLICT_DETAIL_COLUMNS = [
    "sequence",
    "peptide_length",
    "source_count",
    "source_name_set",
    "record_count_total",
    "record_count_ahtpdb",
    "record_count_biopep",
    "ic50_uM_min",
    "ic50_uM_max",
    "ic50_uM_median",
    "ic50_spread_ratio_max_min",
    "cross_source_overlap_flag",
    "stability_flag",
    "stability_note",
    "conflict_reason",
    "conflict_severity",
    "merged_record_id_list",
    "source_record_id_list",
    "database_reference_raw_set",
]


def infer_project_root(script_file: Path) -> Path:
    try:
        return script_file.resolve().parents[2]
    except IndexError:
        return script_file.resolve().parent


def strip_bom(text: str) -> str:
    return text.lstrip("\ufeff")


def normalize_spaces(text: str) -> str:
    text = text.replace("\ufeff", " ")
    text = text.replace("\r", " ").replace("\n", " ")
    text = " ".join(text.split())
    return text.strip()


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


def truthy_flag(value: object) -> bool:
    text = clean_text(value)
    if text is None:
        return False
    return text.lower() in {"1", "true", "yes", "y"}


def read_csv_rows(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"找不到输入文件：{path}")

    with open(path, "r", encoding="utf-8-sig", errors="replace", newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError(f"CSV 无表头：{path}")
        rows: List[Dict[str, str]] = []
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


def classify_tier(stability_flag: Optional[str]) -> str:
    flag = clean_text(stability_flag) or ""

    if flag == "cross_source_consistent":
        return "Tier_A"
    if flag == "within_range_consistent":
        return "Tier_B"
    if flag == "singleton":
        return "Tier_C"
    return "Tier_X"


def tier_rank(tier: str) -> int:
    mapping = {
        "Tier_A": 1,
        "Tier_B": 2,
        "Tier_C": 3,
        "Tier_X": 9,
    }
    return mapping.get(tier, 99)


def build_consensus_note(row: Dict[str, str], tier: str) -> str:
    source_count = to_int_safe(row.get("source_count"))
    record_count = to_int_safe(row.get("record_count_total"))
    spread_ratio = to_float_safe(row.get("ic50_spread_ratio_max_min"))
    stability_flag = clean_text(row.get("stability_flag")) or ""

    if tier == "Tier_A":
        if spread_ratio is not None:
            return f"跨来源一致；source_count={source_count}; record_count={record_count}; spread_ratio={spread_ratio:.4f}"
        return f"跨来源一致；source_count={source_count}; record_count={record_count}"

    if tier == "Tier_B":
        if spread_ratio is not None:
            return f"单源或局部多源范围内一致；source_count={source_count}; record_count={record_count}; spread_ratio={spread_ratio:.4f}"
        return f"单源或局部多源范围内一致；source_count={source_count}; record_count={record_count}"

    if tier == "Tier_C":
        return f"单条证据；source_count={source_count}; record_count={record_count}"

    if spread_ratio is not None:
        return f"存在冲突或状态未明；stability_flag={stability_flag}; source_count={source_count}; record_count={record_count}; spread_ratio={spread_ratio:.4f}"
    return f"存在冲突或状态未明；stability_flag={stability_flag}; source_count={source_count}; record_count={record_count}"


def build_priority_score(row: Dict[str, str], tier: str) -> float:
    """
    一个很轻量的优先级分数：
    - Tier_A 最高
    - Tier_B 次之
    - Tier_C 再次
    - Tier_X 最低
    同 tier 内，source_count 越高、record_count 越高、spread_ratio 越低越靠前
    """
    base = {
        "Tier_A": 3000.0,
        "Tier_B": 2000.0,
        "Tier_C": 1000.0,
        "Tier_X": 0.0,
    }.get(tier, 0.0)

    source_count = float(to_int_safe(row.get("source_count")) or 0)
    record_count = float(to_int_safe(row.get("record_count_total")) or 0)
    spread_ratio = to_float_safe(row.get("ic50_spread_ratio_max_min"))

    if spread_ratio is None or spread_ratio <= 0:
        spread_bonus = 0.0
    else:
        spread_bonus = max(0.0, 100.0 - min(spread_ratio, 100.0))

    return base + source_count * 20.0 + record_count * 5.0 + spread_bonus


def annotate_rows(rows: List[Dict[str, str]]) -> List[Dict[str, object]]:
    out: List[Dict[str, object]] = []

    for row in rows:
        annotated = dict(row)
        tier = classify_tier(row.get("stability_flag"))
        keep_flag = 1 if tier in {"Tier_A", "Tier_B"} else 0
        priority = build_priority_score(row, tier)

        annotated["consensus_tier"] = tier
        annotated["consensus_keep_flag"] = keep_flag
        annotated["consensus_rank"] = tier_rank(tier)
        annotated["consensus_note"] = build_consensus_note(row, tier)
        annotated["priority_score"] = round(priority, 4)

        out.append(annotated)

    out.sort(
        key=lambda r: (
            int(r.get("consensus_rank", 99)),
            -(float(r.get("priority_score", 0.0))),
            to_int_safe(r.get("peptide_length")) or 999,
            str(r.get("sequence") or ""),
        )
    )
    return out


def build_conflict_reason(row: Dict[str, str]) -> str:
    stability_flag = clean_text(row.get("stability_flag")) or ""
    source_count = to_int_safe(row.get("source_count")) or 0
    spread_ratio = to_float_safe(row.get("ic50_spread_ratio_max_min"))
    record_count = to_int_safe(row.get("record_count_total")) or 0

    if stability_flag == "conflicting":
        if source_count >= 2 and spread_ratio is not None and spread_ratio > 10:
            return "cross_source_large_spread"
        if spread_ratio is not None and spread_ratio > 10:
            return "large_within_sequence_spread"
        return "conflicting_flag_from_merged_input"

    if stability_flag == "unknown":
        if spread_ratio is None:
            return "spread_ratio_missing"
        return "unknown_flag_from_merged_input"

    if record_count > 1 and spread_ratio is not None and spread_ratio > 10:
        return "manual_conflict_override_large_spread"

    return "other_unresolved_case"


def build_conflict_severity(row: Dict[str, str]) -> str:
    spread_ratio = to_float_safe(row.get("ic50_spread_ratio_max_min"))
    source_count = to_int_safe(row.get("source_count")) or 0

    if spread_ratio is None:
        return "unknown"

    if source_count >= 2 and spread_ratio >= 100:
        return "very_high"
    if spread_ratio >= 100:
        return "high"
    if spread_ratio >= 10:
        return "medium"
    return "low"


def build_conflict_details(rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
    out = []
    for row in rows:
        detail = {
            "sequence": row.get("sequence"),
            "peptide_length": row.get("peptide_length"),
            "source_count": row.get("source_count"),
            "source_name_set": row.get("source_name_set"),
            "record_count_total": row.get("record_count_total"),
            "record_count_ahtpdb": row.get("record_count_ahtpdb"),
            "record_count_biopep": row.get("record_count_biopep"),
            "ic50_uM_min": row.get("ic50_uM_min"),
            "ic50_uM_max": row.get("ic50_uM_max"),
            "ic50_uM_median": row.get("ic50_uM_median"),
            "ic50_spread_ratio_max_min": row.get("ic50_spread_ratio_max_min"),
            "cross_source_overlap_flag": row.get("cross_source_overlap_flag"),
            "stability_flag": row.get("stability_flag"),
            "stability_note": row.get("stability_note"),
            "conflict_reason": build_conflict_reason(row),
            "conflict_severity": build_conflict_severity(row),
            "merged_record_id_list": row.get("merged_record_id_list"),
            "source_record_id_list": row.get("source_record_id_list"),
            "database_reference_raw_set": row.get("database_reference_raw_set"),
        }
        out.append(detail)

    out.sort(
        key=lambda r: (
            {"very_high": 1, "high": 2, "medium": 3, "low": 4, "unknown": 5}.get(str(r["conflict_severity"]), 9),
            -(to_float_safe(r.get("ic50_spread_ratio_max_min")) or -1.0),
            str(r.get("sequence") or ""),
        )
    )
    return out


def build_summary_rows(annotated_rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
    summary: List[Dict[str, object]] = []

    total = len(annotated_rows)
    tier_counter = Counter(clean_text(r.get("consensus_tier")) or "" for r in annotated_rows)
    stability_counter = Counter(clean_text(r.get("stability_flag")) or "" for r in annotated_rows)
    overlap_counter = Counter(str(int(bool(truthy_flag(r.get("cross_source_overlap_flag"))))) for r in annotated_rows)

    summary.append({"metric": "input_sequence_rows_total", "value": total})

    for tier in ["Tier_A", "Tier_B", "Tier_C", "Tier_X"]:
        summary.append({"metric": f"consensus_tier::{tier}", "value": tier_counter.get(tier, 0)})

    summary.append({
        "metric": "highconf_union_rows",
        "value": tier_counter.get("Tier_A", 0) + tier_counter.get("Tier_B", 0),
    })

    for k, v in stability_counter.items():
        summary.append({"metric": f"stability_flag::{k}", "value": v})

    for k, v in overlap_counter.items():
        label = "cross_source_overlap_yes" if k == "1" else "cross_source_overlap_no"
        summary.append({"metric": label, "value": v})

    # 按长度统计
    for tier in ["Tier_A", "Tier_B", "Tier_C", "Tier_X"]:
        subset = [r for r in annotated_rows if clean_text(r.get("consensus_tier")) == tier]
        length_counter = Counter(to_int_safe(r.get("peptide_length")) for r in subset if to_int_safe(r.get("peptide_length")) is not None)
        for length in sorted(length_counter):
            summary.append({
                "metric": f"{tier}::length::{length}",
                "value": length_counter[length],
            })

    # spread ratio 粗略分桶
    buckets = {
        "le_2": lambda x: x is not None and x <= 2,
        "gt_2_le_5": lambda x: x is not None and 2 < x <= 5,
        "gt_5_le_10": lambda x: x is not None and 5 < x <= 10,
        "gt_10_le_100": lambda x: x is not None and 10 < x <= 100,
        "gt_100": lambda x: x is not None and x > 100,
        "missing": lambda x: x is None,
    }
    bucket_counter = Counter()
    for row in annotated_rows:
        x = to_float_safe(row.get("ic50_spread_ratio_max_min"))
        for name, func in buckets.items():
            if func(x):
                bucket_counter[name] += 1
                break

    for name, value in bucket_counter.items():
        summary.append({"metric": f"spread_bucket::{name}", "value": value})

    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description="构建 merged 短肽共识分层表（纯标准库版）")
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
        help="输入路径；默认 DB/worksets/ace/merged/ace_short_2_3_merged_sequence_level_v0_1.csv",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="输出目录；默认 DB/worksets/ace/merged",
    )
    args = parser.parse_args()

    script_file = Path(__file__)
    project_root = Path(args.project_root) if args.project_root else infer_project_root(script_file)

    input_path = Path(args.input_path) if args.input_path else (
        project_root / "DB" / "worksets" / "ace" / "merged" / INPUT_FILE
    )
    output_dir = Path(args.output_dir) if args.output_dir else (
        project_root / "DB" / "worksets" / "ace" / "merged"
    )

    rows = read_csv_rows(input_path)
    if not rows:
        raise ValueError(f"输入表为空：{input_path}")

    annotated_rows = annotate_rows(rows)

    tier_a_rows = [r for r in annotated_rows if clean_text(r.get("consensus_tier")) == "Tier_A"]
    tier_b_rows = [r for r in annotated_rows if clean_text(r.get("consensus_tier")) == "Tier_B"]
    tier_c_rows = [r for r in annotated_rows if clean_text(r.get("consensus_tier")) == "Tier_C"]
    tier_x_rows = [r for r in annotated_rows if clean_text(r.get("consensus_tier")) == "Tier_X"]
    highconf_rows = [r for r in annotated_rows if clean_text(r.get("consensus_tier")) in {"Tier_A", "Tier_B"}]

    conflict_detail_rows = build_conflict_details(tier_x_rows)
    summary_rows = build_summary_rows(annotated_rows)

    annotated_columns = BASE_COLUMNS + ANNOTATED_EXTRA_COLUMNS

    write_csv(output_dir / ANNOTATED_FILE, annotated_rows, annotated_columns)
    write_csv(output_dir / TIER_A_FILE, tier_a_rows, annotated_columns)
    write_csv(output_dir / TIER_B_FILE, tier_b_rows, annotated_columns)
    write_csv(output_dir / TIER_C_FILE, tier_c_rows, annotated_columns)
    write_csv(output_dir / TIER_X_FILE, tier_x_rows, annotated_columns)
    write_csv(output_dir / HIGHCONF_FILE, highconf_rows, annotated_columns)
    write_csv(output_dir / CONFLICT_DETAILS_FILE, conflict_detail_rows, CONFLICT_DETAIL_COLUMNS)
    write_csv(output_dir / SUMMARY_FILE, summary_rows, ["metric", "value"])

    print("=" * 80)
    print("merged 短肽共识分层构建完成")
    print(f"项目根目录：{project_root}")
    print(f"输入文件：{input_path}")
    print(f"输出目录：{output_dir}")
    print(f"总 sequence-level 条数：{len(annotated_rows)}")
    print(f"Tier_A：{len(tier_a_rows)}")
    print(f"Tier_B：{len(tier_b_rows)}")
    print(f"Tier_C：{len(tier_c_rows)}")
    print(f"Tier_X：{len(tier_x_rows)}")
    print(f"HighConf Union：{len(highconf_rows)}")
    print("=" * 80)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())