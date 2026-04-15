#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ACE short 2/3 benchmark candidates v0.1（纯标准库版，修正版）
==========================================================

输入
----
DB/worksets/ace/merged/ace_short_2_3_consensus_tiers_v0_1.csv

输出
----
DB/worksets/ace/merged/
├── ace_short_2_3_benchmark_candidates_all_v0_1.csv
├── ace_short_2_3_benchmark_high_confidence_core_v0_1.csv
├── ace_short_2_3_benchmark_single_source_candidates_v0_1.csv
├── ace_short_2_3_benchmark_conflict_review_pool_v0_1.csv
└── ace_short_2_3_benchmark_candidates_summary_v0_1.csv

用途
----
- 把共识分层结果整理成更适合 benchmark / 论文 / 工作集使用的三层结构
- 对原来的 Tier D / Tier E 分别映射到：
  * Tier_D_conflicting
  * Tier_E_insufficient

规则
----
1) high_confidence_core
   - 原 Tier A / Tier B
2) single_source_candidates
   - 原 Tier C
3) conflict_review_pool
   - 原 Tier D / Tier E
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional

csv.field_size_limit(10**8)

INPUT_FILE = "ace_short_2_3_consensus_tiers_v0_1.csv"

ALL_OUTPUT = "ace_short_2_3_benchmark_candidates_all_v0_1.csv"
HIGH_CONF_OUTPUT = "ace_short_2_3_benchmark_high_confidence_core_v0_1.csv"
SINGLE_SOURCE_OUTPUT = "ace_short_2_3_benchmark_single_source_candidates_v0_1.csv"
CONFLICT_OUTPUT = "ace_short_2_3_benchmark_conflict_review_pool_v0_1.csv"
SUMMARY_OUTPUT = "ace_short_2_3_benchmark_candidates_summary_v0_1.csv"

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

OUTPUT_COLUMNS = [
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
    "refined_tier",
    "benchmark_bucket",
    "benchmark_priority",
    "benchmark_note",
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


def assign_benchmark_bucket(row: Dict[str, str]) -> Dict[str, object]:
    out = dict(row)

    tier = clean_text(row.get("consensus_tier")) or ""
    stability_flag = clean_text(row.get("stability_flag")) or ""
    source_count = to_int_safe(row.get("source_count")) or 0
    spread_ratio = to_float_safe(row.get("ic50_spread_ratio_max_min"))
    score = to_int_safe(row.get("consensus_score")) or 0

    refined_tier = tier
    bucket = ""
    priority = 0
    note = ""

    if tier in {"Tier A", "Tier B"}:
        refined_tier = tier.replace(" ", "_")
        bucket = "high_confidence_core"
        priority = 1 if tier == "Tier A" else 2
        note = "high-consensus short peptide candidate"

    elif tier == "Tier C":
        refined_tier = "Tier_C_single_source"
        bucket = "single_source_candidates"
        priority = 3
        note = "single-source usable IC50 candidate"

    elif tier == "Tier D":
        refined_tier = "Tier_D_conflicting"
        bucket = "conflict_review_pool"
        priority = 4
        note = "multi-source conflicting short peptide; requires review"

    elif tier == "Tier E":
        refined_tier = "Tier_E_insufficient"
        bucket = "conflict_review_pool"
        priority = 5
        note = "insufficient or unclear evidence pattern"

    else:
        refined_tier = "Tier_E_insufficient"
        bucket = "conflict_review_pool"
        priority = 5
        note = "fallback bucket for unknown tier"

    if spread_ratio is not None:
        note = f"{note}; spread_ratio={spread_ratio:.4f}; source_count={source_count}; score={score}"
    else:
        note = f"{note}; spread_ratio=NA; source_count={source_count}; score={score}"

    out["refined_tier"] = refined_tier
    out["benchmark_bucket"] = bucket
    out["benchmark_priority"] = priority
    out["benchmark_note"] = note
    return out


def build_summary_rows(rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
    summary: List[Dict[str, object]] = []

    summary.append({"metric": "input_rows_total", "value": len(rows)})

    bucket_counter = Counter(clean_text(r.get("benchmark_bucket")) or "" for r in rows)
    for k, v in bucket_counter.items():
        summary.append({"metric": f"benchmark_bucket::{k}", "value": v})

    refined_counter = Counter(clean_text(r.get("refined_tier")) or "" for r in rows)
    for k, v in refined_counter.items():
        summary.append({"metric": f"refined_tier::{k}", "value": v})

    consensus_counter = Counter(clean_text(r.get("consensus_tier")) or "" for r in rows)
    for k, v in consensus_counter.items():
        summary.append({"metric": f"consensus_tier::{k}", "value": v})

    stability_counter = Counter(clean_text(r.get("stability_flag")) or "" for r in rows)
    for k, v in stability_counter.items():
        summary.append({"metric": f"stability_flag::{k}", "value": v})

    source_count_counter = Counter(to_int_safe(r.get("source_count")) or 0 for r in rows)
    for k in sorted(source_count_counter):
        summary.append({"metric": f"source_count::{k}", "value": source_count_counter[k]})

    high_conf_rows = [
        r for r in rows
        if clean_text(r.get("benchmark_bucket")) == "high_confidence_core"
    ]
    single_rows = [
        r for r in rows
        if clean_text(r.get("benchmark_bucket")) == "single_source_candidates"
    ]
    conflict_rows = [
        r for r in rows
        if clean_text(r.get("benchmark_bucket")) == "conflict_review_pool"
    ]

    summary.append({"metric": "high_confidence_core_rows", "value": len(high_conf_rows)})
    summary.append({"metric": "single_source_candidates_rows", "value": len(single_rows)})
    summary.append({"metric": "conflict_review_pool_rows", "value": len(conflict_rows)})

    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description="构建 ACE short 2/3 benchmark candidates v0.1（纯标准库版，修正版）")
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
        help="输入 consensus tiers 文件路径；默认 DB/worksets/ace/merged/ace_short_2_3_consensus_tiers_v0_1.csv",
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

    benchmark_rows = [assign_benchmark_bucket(r) for r in input_rows]

    benchmark_rows.sort(
        key=lambda r: (
            to_int_safe(r.get("benchmark_priority")) if r.get("benchmark_priority") is not None else 999,
            -(to_int_safe(r.get("consensus_score")) or -999999),
            str(r.get("sequence") or ""),
        )
    )

    high_conf_rows = [
        r for r in benchmark_rows
        if clean_text(r.get("benchmark_bucket")) == "high_confidence_core"
    ]
    single_source_rows = [
        r for r in benchmark_rows
        if clean_text(r.get("benchmark_bucket")) == "single_source_candidates"
    ]
    conflict_rows = [
        r for r in benchmark_rows
        if clean_text(r.get("benchmark_bucket")) == "conflict_review_pool"
    ]

    summary_rows = build_summary_rows(benchmark_rows)

    write_csv(output_dir / ALL_OUTPUT, benchmark_rows, OUTPUT_COLUMNS)
    write_csv(output_dir / HIGH_CONF_OUTPUT, high_conf_rows, OUTPUT_COLUMNS)
    write_csv(output_dir / SINGLE_SOURCE_OUTPUT, single_source_rows, OUTPUT_COLUMNS)
    write_csv(output_dir / CONFLICT_OUTPUT, conflict_rows, OUTPUT_COLUMNS)
    write_csv(output_dir / SUMMARY_OUTPUT, summary_rows, ["metric", "value"])

    print("=" * 80)
    print("ACE short 2/3 benchmark candidates v0.1 构建完成（修正版）")
    print(f"项目根目录：{project_root}")
    print(f"输入文件：{input_path}")
    print(f"总记录数：{len(benchmark_rows)}")
    print(f"高置信核心层：{len(high_conf_rows)}")
    print(f"单源候选层：{len(single_source_rows)}")
    print(f"冲突复核层：{len(conflict_rows)}")
    print(f"输出：{output_dir / ALL_OUTPUT}")
    print(f"输出：{output_dir / HIGH_CONF_OUTPUT}")
    print(f"输出：{output_dir / SINGLE_SOURCE_OUTPUT}")
    print(f"输出：{output_dir / CONFLICT_OUTPUT}")
    print(f"输出：{output_dir / SUMMARY_OUTPUT}")
    print("=" * 80)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())