#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ACE short 2/3 merged consensus tiers v0.1（纯标准库版）
====================================================

输入
----
DB/worksets/ace/merged/ace_short_2_3_merged_sequence_level_v0_3.csv

输出
----
DB/worksets/ace/merged/
├── ace_short_2_3_consensus_tiers_v0_1.csv
├── ace_short_2_3_consensus_high_confidence_v0_1.csv
├── ace_short_2_3_consensus_conflicting_v0_1.csv
└── ace_short_2_3_consensus_summary_v0_1.csv

分层规则
--------
Tier A:
    - source_count >= 3
    - stability_flag in {cross_source_consistent, within_range_consistent}
    - spread_ratio <= 10
Tier B:
    - source_count == 2
    - stability_flag in {cross_source_consistent, within_range_consistent}
    - spread_ratio <= 10
Tier C:
    - source_count == 1
    - 有可用 ic50_uM_median
Tier D:
    - source_count >= 2
    - stability_flag == conflicting
Tier E:
    - 其余情况（unknown / 信息不足 / spread_ratio 缺失等）

说明
----
- 这是 sequence-level 共识层，不回写 merged 原表
- 重点服务于二肽/三肽高置信工作集与冲突分析
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional

csv.field_size_limit(10**8)

INPUT_FILE = "ace_short_2_3_merged_sequence_level_v0_3.csv"

CONSENSUS_OUTPUT = "ace_short_2_3_consensus_tiers_v0_1.csv"
HIGH_CONF_OUTPUT = "ace_short_2_3_consensus_high_confidence_v0_1.csv"
CONFLICT_OUTPUT = "ace_short_2_3_consensus_conflicting_v0_1.csv"
SUMMARY_OUTPUT = "ace_short_2_3_consensus_summary_v0_1.csv"

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

CONSENSUS_COLUMNS = [
    "sequence",
    "peptide_length",
    "source_count",
    "source_name_set",
    "record_count_total",
    "record_count_ahtpdb",
    "record_count_biopep",
    "record_count_mbpdb",
    "record_count_fermfooddb",
    "ic50_uM_min",
    "ic50_uM_max",
    "ic50_uM_mean",
    "ic50_uM_median",
    "ic50_spread_ratio_max_min",
    "cross_source_overlap_flag",
    "stability_flag",
    "stability_note",
    "consensus_tier",
    "consensus_score",
    "consensus_note",
    "high_confidence_flag",
    "conflict_flag",
    "source_richness_flag",
    "bibliographic_year_min",
    "bibliographic_year_max",
    "merged_record_id_list",
    "source_record_id_list",
    "species_set",
    "protein_id_set",
    "doi_set",
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


def compute_consensus_score(
    source_count: int,
    stability_flag: str,
    spread_ratio: Optional[float],
    record_count_total: int,
) -> int:
    """
    粗粒度打分，便于后续排序
    """
    score = 0

    # 来源数
    score += min(source_count, 4) * 20

    # 记录条数
    score += min(record_count_total, 10) * 2

    # 稳定性
    if stability_flag == "cross_source_consistent":
        score += 25
    elif stability_flag == "within_range_consistent":
        score += 18
    elif stability_flag == "singleton":
        score += 8
    elif stability_flag == "unknown":
        score += 3
    elif stability_flag == "conflicting":
        score -= 20

    # 数值范围
    if spread_ratio is not None:
        if spread_ratio <= 3:
            score += 15
        elif spread_ratio <= 10:
            score += 8
        elif spread_ratio <= 30:
            score -= 5
        else:
            score -= 12

    return score


def assign_consensus_tier(row: Dict[str, str]) -> Dict[str, object]:
    source_count = to_int_safe(row.get("source_count")) or 0
    record_count_total = to_int_safe(row.get("record_count_total")) or 0
    spread_ratio = to_float_safe(row.get("ic50_spread_ratio_max_min"))
    stability_flag = clean_text(row.get("stability_flag")) or ""
    ic50_median = to_float_safe(row.get("ic50_uM_median"))

    tier = "Tier E"
    note = ""
    high_conf = 0
    conflict_flag = 0
    source_richness_flag = 1 if source_count >= 3 else 0

    consistent_flags = {"cross_source_consistent", "within_range_consistent"}

    if (
        source_count >= 3
        and stability_flag in consistent_flags
        and spread_ratio is not None
        and spread_ratio <= 10
        and ic50_median is not None
    ):
        tier = "Tier A"
        note = ">=3 sources and consistent within acceptable spread"
        high_conf = 1

    elif (
        source_count == 2
        and stability_flag in consistent_flags
        and spread_ratio is not None
        and spread_ratio <= 10
        and ic50_median is not None
    ):
        tier = "Tier B"
        note = "2 sources and consistent within acceptable spread"
        high_conf = 1

    elif (
        source_count == 1
        and ic50_median is not None
    ):
        tier = "Tier C"
        note = "single-source evidence with usable IC50"

    elif (
        source_count >= 2
        and stability_flag == "conflicting"
    ):
        tier = "Tier D"
        note = "multi-source evidence but conflicting IC50 range"
        conflict_flag = 1

    else:
        tier = "Tier E"
        note = "unknown or insufficient evidence pattern"

    score = compute_consensus_score(
        source_count=source_count,
        stability_flag=stability_flag,
        spread_ratio=spread_ratio,
        record_count_total=record_count_total,
    )

    out = dict(row)
    out["consensus_tier"] = tier
    out["consensus_score"] = score
    out["consensus_note"] = note
    out["high_confidence_flag"] = high_conf
    out["conflict_flag"] = conflict_flag
    out["source_richness_flag"] = source_richness_flag
    return out


def build_summary_rows(consensus_rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
    summary: List[Dict[str, object]] = []
    summary.append({"metric": "input_sequence_rows_total", "value": len(consensus_rows)})

    tier_counter = Counter(clean_text(r.get("consensus_tier")) or "" for r in consensus_rows)
    for k, v in tier_counter.items():
        summary.append({"metric": f"consensus_tier::{k}", "value": v})

    stability_counter = Counter(clean_text(r.get("stability_flag")) or "" for r in consensus_rows)
    for k, v in stability_counter.items():
        summary.append({"metric": f"stability_flag::{k}", "value": v})

    high_conf_count = sum(1 for r in consensus_rows if to_int_safe(r.get("high_confidence_flag")) == 1)
    summary.append({"metric": "high_confidence_rows", "value": high_conf_count})

    conflict_count = sum(1 for r in consensus_rows if to_int_safe(r.get("conflict_flag")) == 1)
    summary.append({"metric": "conflicting_rows", "value": conflict_count})

    overlap_yes = sum(1 for r in consensus_rows if to_int_safe(r.get("cross_source_overlap_flag")) == 1)
    overlap_no = sum(1 for r in consensus_rows if to_int_safe(r.get("cross_source_overlap_flag")) == 0)
    summary.append({"metric": "cross_source_overlap_yes", "value": overlap_yes})
    summary.append({"metric": "cross_source_overlap_no", "value": overlap_no})

    source_count_counter = Counter(to_int_safe(r.get("source_count")) or 0 for r in consensus_rows)
    for k in sorted(source_count_counter):
        summary.append({"metric": f"source_count::{k}", "value": source_count_counter[k]})

    length_counter = Counter(to_int_safe(r.get("peptide_length")) or 0 for r in consensus_rows)
    for k in sorted(length_counter):
        summary.append({"metric": f"peptide_length::{k}", "value": length_counter[k]})

    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description="构建 ACE short 2/3 merged consensus tiers v0.1（纯标准库版）")
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
        help="输入 merged short sequence-level 文件路径；默认 DB/worksets/ace/merged/ace_short_2_3_merged_sequence_level_v0_3.csv",
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

    input_rows = read_csv_rows(input_path)
    if not input_rows:
        raise ValueError(f"输入表为空：{input_path}")

    consensus_rows = [assign_consensus_tier(r) for r in input_rows]

    # 排序：tier 优先，其次 score 降序，再按 sequence
    tier_rank = {"Tier A": 1, "Tier B": 2, "Tier C": 3, "Tier D": 4, "Tier E": 5}
    consensus_rows.sort(
        key=lambda r: (
            tier_rank.get(clean_text(r.get("consensus_tier")) or "Tier E", 99),
            -(to_int_safe(r.get("consensus_score")) or -999999),
            str(r.get("sequence") or ""),
        )
    )

    high_conf_rows = [
        r for r in consensus_rows
        if clean_text(r.get("consensus_tier")) in {"Tier A", "Tier B"}
    ]

    conflict_rows = [
        r for r in consensus_rows
        if clean_text(r.get("consensus_tier")) == "Tier D"
    ]

    summary_rows = build_summary_rows(consensus_rows)

    write_csv(output_dir / CONSENSUS_OUTPUT, consensus_rows, CONSENSUS_COLUMNS)
    write_csv(output_dir / HIGH_CONF_OUTPUT, high_conf_rows, CONSENSUS_COLUMNS)
    write_csv(output_dir / CONFLICT_OUTPUT, conflict_rows, CONSENSUS_COLUMNS)
    write_csv(output_dir / SUMMARY_OUTPUT, summary_rows, ["metric", "value"])

    print("=" * 80)
    print("ACE short 2/3 consensus tiers v0.1 构建完成")
    print(f"项目根目录：{project_root}")
    print(f"输入文件：{input_path}")
    print(f"总序列数：{len(consensus_rows)}")
    print(f"高置信（Tier A/B）：{len(high_conf_rows)}")
    print(f"冲突层（Tier D）：{len(conflict_rows)}")
    print(f"输出：{output_dir / CONSENSUS_OUTPUT}")
    print(f"输出：{output_dir / HIGH_CONF_OUTPUT}")
    print(f"输出：{output_dir / CONFLICT_OUTPUT}")
    print(f"输出：{output_dir / SUMMARY_OUTPUT}")
    print("=" * 80)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())