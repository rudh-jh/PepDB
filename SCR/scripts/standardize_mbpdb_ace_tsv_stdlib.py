#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MBPDB ACE TSV 标准化脚本（纯标准库版，修正版）
=============================================

输入
----
默认从以下目录自动寻找一个 TSV：
DB/raw/ace/databases/MBPDB/manual_exports/

输出
----
DB/standardized/ace/mbpdb/
├── mbpdb_ace_master_clean.csv
├── mbpdb_ace_parse_log.csv
└── mbpdb_ace_standardize_summary.csv

DB/analysis/ace/mbpdb/标准化检查/
├── 字段缺失统计.csv
├── 解析状态统计.csv
├── 重复sequence检查.csv
└── 物种分布统计.csv

说明
----
- 当前输入应为 MBPDB 的 ACE-inhibitory 手工导出 TSV
- 不做网页抓取，不做详情页 enrich
- IC50 列标题是 "IC50 (μM)"，默认统一解析为 uM
- 对一个单元格内出现多个 IC50 数值的情况，默认取中位数作为 ic50_value/ic50_uM，
  同时在 parse_log 和 notes 中保留痕迹
- 修复点：
  1) 不再把 1) / 2) / 3) 这种编号当成 IC50 数值
  2) 自动跳过空白异常行
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import statistics
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Tuple


INPUT_GLOB_PATTERNS = [
    "mbpdb_ace_inhibitory_export*.tsv",
    "*.tsv",
]

REQUIRED_COLUMNS = [
    "Protein ID",
    "Peptide",
    "Protein description",
    "Species",
    "Intervals",
    "Function",
    "Additional details",
    "IC50 (μM)",
    "Inhibition type",
    "Inhibited microorganisms",
    "PTM",
    "Title",
    "Authors",
    "Abstract",
    "DOI",
]

MASTER_COLUMNS = [
    "record_id",
    "task",
    "source_type",
    "source_name",
    "source_record_id",
    "target",
    "sequence",
    "peptide_length",
    "activity_label_raw",
    "ic50_raw",
    "ic50_value",
    "ic50_unit",
    "ic50_relation",
    "ic50_uM",
    "ic50_parse_status",
    "ic50_parse_note",
    "protein_id_raw",
    "protein_description_raw",
    "species_raw",
    "intervals_raw",
    "function_raw",
    "additional_details_raw",
    "inhibition_type_raw",
    "inhibited_microorganisms_raw",
    "ptm_raw",
    "title_raw",
    "authors_raw",
    "abstract_raw",
    "doi_raw",
    "notes",
]

CANONICAL_AA = set("ACDEFGHIKLMNPQRSTVWY")
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


def normalize_spaces(text: str) -> str:
    text = text.replace("\ufeff", " ")
    text = text.replace("\r", " ").replace("\n", " ")
    return re.sub(r"\s+", " ", text).strip()


def strip_bom(text: str) -> str:
    return text.lstrip("\ufeff")


def is_missing(value: object) -> bool:
    if value is None:
        return True
    text = normalize_spaces(str(value)).lower()
    return text in MISSING_TOKENS


def clean_text(value: object) -> Optional[str]:
    if is_missing(value):
        return None
    return normalize_spaces(str(value))


def clean_sequence(value: object) -> Optional[str]:
    text = clean_text(value)
    if text is None:
        return None
    text = text.upper().replace(" ", "")
    return text if text else None


def sequence_is_canonical(seq: Optional[str]) -> bool:
    if not seq:
        return False
    return set(seq).issubset(CANONICAL_AA)


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


def write_csv(path: Path, rows: List[Dict[str, object]], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def read_tsv_rows(path: Path) -> List[Dict[str, Optional[str]]]:
    if not path.exists():
        raise FileNotFoundError(f"找不到输入文件：{path}")

    with open(path, "r", encoding="utf-8-sig", errors="replace", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV 无表头：{path}")

        fieldnames = [strip_bom(x) for x in reader.fieldnames]
        missing = [c for c in REQUIRED_COLUMNS if c not in fieldnames]
        if missing:
            raise ValueError(
                f"文件 {path.name} 缺少必要列：{missing}\n"
                f"当前列：{fieldnames}"
            )

        rows: List[Dict[str, Optional[str]]] = []
        for row in reader:
            clean_row: Dict[str, Optional[str]] = {}
            for key, value in row.items():
                clean_key = strip_bom(key)
                clean_row[clean_key] = clean_text(value)
            rows.append(clean_row)
    return rows


def resolve_input_tsv(input_path: Optional[str], project_root: Path) -> Path:
    if input_path:
        return Path(input_path)

    manual_dir = project_root / "DB" / "raw" / "ace" / "databases" / "MBPDB" / "manual_exports"
    if not manual_dir.exists():
        raise FileNotFoundError(f"找不到目录：{manual_dir}")

    candidates: List[Path] = []
    for pattern in INPUT_GLOB_PATTERNS:
        candidates.extend(sorted(manual_dir.glob(pattern)))

    seen = set()
    unique_candidates = []
    for p in candidates:
        rp = p.resolve()
        if rp not in seen:
            seen.add(rp)
            unique_candidates.append(p)

    if not unique_candidates:
        raise FileNotFoundError(
            f"在 {manual_dir} 下没有找到 TSV 文件。"
            f"建议放置如 mbpdb_ace_inhibitory_export_2026-04-14.tsv"
        )

    return unique_candidates[0]


def build_source_record_id(row: Dict[str, Optional[str]]) -> str:
    protein_id = clean_text(row.get("Protein ID")) or "NA"
    peptide = clean_sequence(row.get("Peptide")) or "NA"
    intervals = clean_text(row.get("Intervals")) or "NA"
    intervals = re.sub(r"[^A-Za-z0-9._-]+", "_", intervals)
    return f"MBPDB_{protein_id}_{peptide}_{intervals}"


def parse_ic50_um(raw_value: object) -> Dict[str, object]:
    """
    解析 TSV 中的 "IC50 (μM)" 列。

    常见格式：
    - 52.8
    - < 20
    - 1) 315.0, 2) 205.0, 3) 315.0
    - 315.0; 205.0; 315.0
    - 空值

    规则：
    - 单值：exact_molar
    - 阈值：molar_threshold
    - 多值：取中位数作为 ic50_value / ic50_uM
    - 关键修复：先移除 1) / 2) / 3) 这种编号，避免把编号当成数值
    """
    result = {
        "ic50_value": None,
        "ic50_unit": None,
        "ic50_relation": None,
        "ic50_uM": None,
        "ic50_parse_status": "missing",
        "ic50_parse_note": None,
        "ic50_numeric_count": 0,
        "ic50_numeric_min": None,
        "ic50_numeric_max": None,
        "ic50_numeric_median": None,
        "ic50_numeric_values_serialized": None,
    }

    text = clean_text(raw_value)
    if text is None:
        return result

    work_text = text
    relation = "="
    if work_text.startswith("<"):
        relation = "<"
    elif work_text.startswith(">"):
        relation = ">"
    elif work_text.startswith("≤"):
        relation = "≤"
    elif work_text.startswith("≥"):
        relation = "≥"

    # 去掉列表编号，例如：1) 315.0, 2) 205.0, 3) 315.0
    work_text = re.sub(r'(?<!\d)\b\d+\s*[\)\.]\s*', ' ', work_text)

    # 统一常见分隔符
    work_text = work_text.replace("；", ";").replace("，", ",")
    work_text = work_text.replace("/", " / ")

    # 只抓真正数值；不再把编号抓进去
    numeric_strings = re.findall(r'(?<![A-Za-z])[-+]?\d+(?:\.\d+)?(?!\s*[\)\.])', work_text)

    values = []
    for x in numeric_strings:
        try:
            values.append(float(x))
        except ValueError:
            continue

    if not values:
        result["ic50_parse_status"] = "unparsed_text"
        result["ic50_parse_note"] = "text_present_but_no_numeric_value_found"
        return result

    values_sorted = sorted(values)
    median_value = statistics.median(values_sorted)

    result["ic50_value"] = median_value
    result["ic50_unit"] = "uM"
    result["ic50_relation"] = relation
    result["ic50_uM"] = median_value
    result["ic50_parse_status"] = "exact_molar" if relation == "=" else "molar_threshold"
    result["ic50_numeric_count"] = len(values_sorted)
    result["ic50_numeric_min"] = min(values_sorted)
    result["ic50_numeric_max"] = max(values_sorted)
    result["ic50_numeric_median"] = median_value
    result["ic50_numeric_values_serialized"] = " | ".join(str(v) for v in values_sorted)

    if len(values_sorted) == 1:
        result["ic50_parse_note"] = "single_numeric_uM_value_from_tsv"
    else:
        result["ic50_parse_note"] = f"multi_value_median_used_n={len(values_sorted)}"

    return result


def is_effectively_empty_row(row: Dict[str, Optional[str]]) -> bool:
    """
    过滤掉导出 TSV 中偶发的空白/异常行。
    只要核心字段几乎全空，就不进入主表。
    """
    core_fields = [
        row.get("Protein ID"),
        row.get("Peptide"),
        row.get("Function"),
        row.get("IC50 (μM)"),
        row.get("Title"),
        row.get("Authors"),
        row.get("DOI"),
    ]
    return all(clean_text(x) is None for x in core_fields)


def build_notes(row: Dict[str, Optional[str]], sequence: Optional[str], ic50_info: Dict[str, object]) -> Optional[str]:
    notes = []

    function_raw = clean_text(row.get("Function"))
    if function_raw and function_raw != "ACE-inhibitory":
        notes.append("function_not_exactly_ACE-inhibitory")

    if sequence is None:
        notes.append("sequence_missing")
    elif not sequence_is_canonical(sequence):
        notes.append("non_canonical_or_contains_unknown_residue")

    if ic50_info["ic50_parse_status"] == "missing":
        notes.append("ic50_missing")
    elif ic50_info["ic50_parse_status"] == "unparsed_text":
        notes.append("ic50_needs_manual_review")

    if isinstance(ic50_info.get("ic50_numeric_count"), int) and ic50_info["ic50_numeric_count"] > 1:
        notes.append("multiple_ic50_values_in_single_cell")

    if clean_text(row.get("PTM")):
        notes.append("ptm_present")

    if clean_text(row.get("Inhibition type")):
        notes.append("inhibition_type_present")

    if not clean_text(row.get("DOI")):
        notes.append("doi_missing")

    if not clean_text(row.get("Title")):
        notes.append("title_missing")

    return "; ".join(notes) if notes else None


def build_master_row(row: Dict[str, Optional[str]], idx: int) -> Tuple[Dict[str, object], Dict[str, object]]:
    sequence = clean_sequence(row.get("Peptide"))
    peptide_length = len(sequence) if sequence else None
    ic50_raw = clean_text(row.get("IC50 (μM)"))
    ic50_info = parse_ic50_um(ic50_raw)

    source_record_id = build_source_record_id(row)

    master_row = {
        "record_id": f"mbpdb_ace_{idx:06d}",
        "task": "ace",
        "source_type": "database",
        "source_name": "MBPDB",
        "source_record_id": source_record_id,
        "target": "ACE",
        "sequence": sequence,
        "peptide_length": peptide_length,
        "activity_label_raw": clean_text(row.get("Function")),
        "ic50_raw": ic50_raw,
        "ic50_value": ic50_info["ic50_value"],
        "ic50_unit": ic50_info["ic50_unit"],
        "ic50_relation": ic50_info["ic50_relation"],
        "ic50_uM": ic50_info["ic50_uM"],
        "ic50_parse_status": ic50_info["ic50_parse_status"],
        "ic50_parse_note": ic50_info["ic50_parse_note"],
        "protein_id_raw": clean_text(row.get("Protein ID")),
        "protein_description_raw": clean_text(row.get("Protein description")),
        "species_raw": clean_text(row.get("Species")),
        "intervals_raw": clean_text(row.get("Intervals")),
        "function_raw": clean_text(row.get("Function")),
        "additional_details_raw": clean_text(row.get("Additional details")),
        "inhibition_type_raw": clean_text(row.get("Inhibition type")),
        "inhibited_microorganisms_raw": clean_text(row.get("Inhibited microorganisms")),
        "ptm_raw": clean_text(row.get("PTM")),
        "title_raw": clean_text(row.get("Title")),
        "authors_raw": clean_text(row.get("Authors")),
        "abstract_raw": clean_text(row.get("Abstract")),
        "doi_raw": clean_text(row.get("DOI")),
        "notes": build_notes(row, sequence, ic50_info),
    }

    parse_log = {
        "record_id": master_row["record_id"],
        "source_record_id": source_record_id,
        "sequence": sequence,
        "peptide_length": peptide_length,
        "ic50_raw": ic50_raw,
        "ic50_value": ic50_info["ic50_value"],
        "ic50_unit": ic50_info["ic50_unit"],
        "ic50_relation": ic50_info["ic50_relation"],
        "ic50_uM": ic50_info["ic50_uM"],
        "ic50_parse_status": ic50_info["ic50_parse_status"],
        "ic50_parse_note": ic50_info["ic50_parse_note"],
        "ic50_numeric_count": ic50_info["ic50_numeric_count"],
        "ic50_numeric_min": ic50_info["ic50_numeric_min"],
        "ic50_numeric_max": ic50_info["ic50_numeric_max"],
        "ic50_numeric_median": ic50_info["ic50_numeric_median"],
        "ic50_numeric_values_serialized": ic50_info["ic50_numeric_values_serialized"],
        "function_raw": clean_text(row.get("Function")),
        "species_raw": clean_text(row.get("Species")),
        "notes": master_row["notes"],
    }

    return master_row, parse_log


def count_missing(rows: List[Dict[str, object]], columns: List[str]) -> List[Dict[str, object]]:
    out = []
    total = len(rows)
    for col in columns:
        missing_count = 0
        for row in rows:
            val = row.get(col)
            if val is None or str(val).strip() == "":
                missing_count += 1
        out.append({
            "field_name": col,
            "total_rows": total,
            "missing_count": missing_count,
            "non_missing_count": total - missing_count,
            "missing_ratio": round(missing_count / total, 6) if total else 0.0,
        })
    return out


def build_counter_rows(counter: Counter, key_name: str, value_name: str) -> List[Dict[str, object]]:
    return [{key_name: k, value_name: v} for k, v in counter.most_common()]


def build_duplicate_sequence_rows(master_rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
    counter = Counter(str(x.get("sequence") or "") for x in master_rows if x.get("sequence"))
    out = []
    for seq, cnt in counter.most_common():
        if seq and cnt > 1:
            out.append({"sequence": seq, "duplicate_count": cnt})
    return out


def build_species_rows(master_rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
    counter = Counter(str(x.get("species_raw") or "") for x in master_rows if x.get("species_raw"))
    return build_counter_rows(counter, "species_raw", "count")


def build_summary_rows(
    input_path: Path,
    input_rows: List[Dict[str, Optional[str]]],
    master_rows: List[Dict[str, object]],
    skipped_empty_rows: int,
) -> List[Dict[str, object]]:
    summary = [
        {"metric": "input_filename", "value": input_path.name},
        {"metric": "input_rows", "value": len(input_rows)},
        {"metric": "skipped_empty_rows", "value": skipped_empty_rows},
        {"metric": "master_rows", "value": len(master_rows)},
    ]

    parse_counter = Counter(str(x.get("ic50_parse_status") or "") for x in master_rows)
    for k, v in parse_counter.items():
        summary.append({"metric": f"ic50_parse_status::{k}", "value": v})

    multi_ic50_rows = sum(
        1 for x in master_rows
        if clean_text(x.get("notes")) and "multiple_ic50_values_in_single_cell" in str(x.get("notes"))
    )
    summary.append({"metric": "rows_with_multiple_ic50_values_in_single_cell", "value": multi_ic50_rows})

    noncanonical_rows = sum(
        1 for x in master_rows
        if clean_text(x.get("notes")) and "non_canonical_or_contains_unknown_residue" in str(x.get("notes"))
    )
    summary.append({"metric": "rows_with_noncanonical_sequence", "value": noncanonical_rows})

    species_counter = Counter(str(x.get("species_raw") or "") for x in master_rows if x.get("species_raw"))
    for k, v in species_counter.most_common():
        summary.append({"metric": f"species::{k}", "value": v})

    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description="MBPDB ACE TSV 标准化脚本（纯标准库版，修正版）")
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
        help="MBPDB ACE TSV 文件路径；默认自动到 DB/raw/ace/databases/MBPDB/manual_exports 下查找。",
    )
    parser.add_argument(
        "--output-path",
        type=str,
        default=None,
        help="主表输出路径；默认 DB/standardized/ace/mbpdb/mbpdb_ace_master_clean.csv",
    )
    args = parser.parse_args()

    script_file = Path(__file__)
    project_root = Path(args.project_root) if args.project_root else infer_project_root(script_file)

    input_path = resolve_input_tsv(args.input_path, project_root)
    output_path = Path(args.output_path) if args.output_path else (
        project_root / "DB" / "standardized" / "ace" / "mbpdb" / "mbpdb_ace_master_clean.csv"
    )

    standardized_dir = output_path.parent
    analysis_dir = project_root / "DB" / "analysis" / "ace" / "mbpdb" / "标准化检查"

    parse_log_path = standardized_dir / "mbpdb_ace_parse_log.csv"
    summary_path = standardized_dir / "mbpdb_ace_standardize_summary.csv"
    missing_stats_path = analysis_dir / "字段缺失统计.csv"
    parse_status_stats_path = analysis_dir / "解析状态统计.csv"
    duplicate_sequence_path = analysis_dir / "重复sequence检查.csv"
    species_stats_path = analysis_dir / "物种分布统计.csv"

    input_rows = read_tsv_rows(input_path)

    master_rows: List[Dict[str, object]] = []
    parse_logs: List[Dict[str, object]] = []
    skipped_empty_rows = 0

    for row in input_rows:
        if is_effectively_empty_row(row):
            skipped_empty_rows += 1
            continue

        idx = len(master_rows) + 1
        master_row, parse_log = build_master_row(row, idx)
        master_rows.append(master_row)
        parse_logs.append(parse_log)

    write_csv(output_path, master_rows, MASTER_COLUMNS)
    write_csv(parse_log_path, parse_logs, list(parse_logs[0].keys()) if parse_logs else [])

    missing_stats = count_missing(master_rows, MASTER_COLUMNS)
    write_csv(
        missing_stats_path,
        missing_stats,
        ["field_name", "total_rows", "missing_count", "non_missing_count", "missing_ratio"],
    )

    parse_counter = Counter(str(x.get("ic50_parse_status") or "") for x in master_rows)
    parse_status_rows = build_counter_rows(parse_counter, "ic50_parse_status", "count")
    write_csv(parse_status_stats_path, parse_status_rows, ["ic50_parse_status", "count"])

    duplicate_sequence_rows = build_duplicate_sequence_rows(master_rows)
    if not duplicate_sequence_rows:
        duplicate_sequence_rows = [{"sequence": "", "duplicate_count": 0}]
    write_csv(duplicate_sequence_path, duplicate_sequence_rows, ["sequence", "duplicate_count"])

    species_rows = build_species_rows(master_rows)
    if not species_rows:
        species_rows = [{"species_raw": "", "count": 0}]
    write_csv(species_stats_path, species_rows, ["species_raw", "count"])

    summary_rows = build_summary_rows(input_path, input_rows, master_rows, skipped_empty_rows)
    write_csv(summary_path, summary_rows, ["metric", "value"])

    print("=" * 80)
    print("MBPDB ACE TSV 标准化完成（修正版）")
    print(f"项目根目录：{project_root}")
    print(f"输入文件：{input_path}")
    print(f"主表输出：{output_path}")
    print(f"解析日志：{parse_log_path}")
    print(f"检查目录：{analysis_dir}")
    print(f"input_rows={len(input_rows)} skipped_empty_rows={skipped_empty_rows} master_rows={len(master_rows)}")
    print("=" * 80)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())