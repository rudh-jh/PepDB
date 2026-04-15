#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FermFooDb ACE core workset 构建脚本（纯标准库版）
==============================================

输入
----
DB/standardized/ace/fermfoodb/fermfoodb_ace_master_rescued.csv

输出
----
DB/worksets/ace/fermfoodb/
├── fermfoodb_ace_core_strict_record_level.csv
├── fermfoodb_ace_core_strict_sequence_level.csv
├── fermfoodb_ace_core_expanded_record_level.csv
├── fermfoodb_ace_core_expanded_sequence_level.csv
├── fermfoodb_ace_short_2_3_strict_record_level.csv
├── fermfoodb_ace_short_2_3_strict_sequence_level.csv
├── fermfoodb_ace_short_2_3_expanded_record_level.csv
├── fermfoodb_ace_short_2_3_expanded_sequence_level.csv
├── fermfoodb_ace_excluded_records.csv
└── fermfoodb_ace_workset_summary.csv

说明
----
- strict: 只保留 exact_molar / molar_threshold
- expanded: 在 strict 基础上再纳入 rescued_numeric_no_unit
- 要求 sequence 为 canonical，且 ic50_uM 非空
- 不回写 standardized 原表
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional

csv.field_size_limit(10**8)

INPUT_FILE = "fermfoodb_ace_master_rescued.csv"

STRICT_RECORD_FILE = "fermfoodb_ace_core_strict_record_level.csv"
STRICT_SEQUENCE_FILE = "fermfoodb_ace_core_strict_sequence_level.csv"
EXPANDED_RECORD_FILE = "fermfoodb_ace_core_expanded_record_level.csv"
EXPANDED_SEQUENCE_FILE = "fermfoodb_ace_core_expanded_sequence_level.csv"

SHORT23_STRICT_RECORD_FILE = "fermfoodb_ace_short_2_3_strict_record_level.csv"
SHORT23_STRICT_SEQUENCE_FILE = "fermfoodb_ace_short_2_3_strict_sequence_level.csv"
SHORT23_EXPANDED_RECORD_FILE = "fermfoodb_ace_short_2_3_expanded_record_level.csv"
SHORT23_EXPANDED_SEQUENCE_FILE = "fermfoodb_ace_short_2_3_expanded_sequence_level.csv"

EXCLUDED_FILE = "fermfoodb_ace_excluded_records.csv"
SUMMARY_FILE = "fermfoodb_ace_workset_summary.csv"

STRICT_STATUSES = {"exact_molar", "molar_threshold"}
EXPANDED_STATUSES = {"exact_molar", "molar_threshold", "rescued_numeric_no_unit"}

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
CANONICAL_AA = set("ACDEFGHIKLMNPQRSTVWY")


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


def clean_sequence(value: object) -> Optional[str]:
    text = clean_text(value)
    if text is None:
        return None
    text = text.upper().replace(" ", "").replace("-", "")
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


def sorted_unique_join(values: List[object], sep: str = " | ") -> str:
    seen = set()
    out = []
    for v in values:
        t = clean_text(v)
        if t is None:
            continue
        if t not in seen:
            seen.add(t)
            out.append(t)
    out.sort()
    return sep.join(out)


def year_from_any_text(value: object) -> Optional[int]:
    text = clean_text(value)
    if text is None:
        return None
    import re
    years = re.findall(r"\b(19\d{2}|20\d{2})\b", text)
    if not years:
        return None
    years_int = [int(y) for y in years]
    return max(years_int)


def safe_min(values: List[float]) -> Optional[float]:
    return min(values) if values else None


def safe_max(values: List[float]) -> Optional[float]:
    return max(values) if values else None


def safe_mean(values: List[float]) -> Optional[float]:
    return statistics.mean(values) if values else None


def safe_median(values: List[float]) -> Optional[float]:
    return statistics.median(values) if values else None


def is_noncanonical_row(row: Dict[str, str]) -> bool:
    seq = clean_sequence(row.get("sequence"))
    if seq is None:
        return True
    if not sequence_is_canonical(seq):
        return True
    if has_note(row, "non_canonical_or_contains_unknown_residue"):
        return True
    return False


def is_short_2_3_row(row: Dict[str, object]) -> bool:
    return to_int_safe(row.get("peptide_length")) in {2, 3}


def build_exclusion_reasons(row: Dict[str, str]) -> List[str]:
    reasons = []

    seq = clean_sequence(row.get("sequence"))
    if seq is None:
        reasons.append("sequence_missing")
    elif not sequence_is_canonical(seq):
        reasons.append("noncanonical_sequence")

    peptide_length = to_int_safe(row.get("peptide_length"))
    if peptide_length is None:
        reasons.append("peptide_length_missing")

    status = clean_text(row.get("ic50_parse_status"))
    ic50_uM = to_float_safe(row.get("ic50_uM"))
    if status not in (STRICT_STATUSES | {"rescued_numeric_no_unit"}):
        reasons.append("ic50_status_not_usable")
    if ic50_uM is None:
        reasons.append("ic50_uM_missing")

    return reasons


def keep_row(row: Dict[str, str], mode: str) -> bool:
    seq = clean_sequence(row.get("sequence"))
    if seq is None:
        return False
    if not sequence_is_canonical(seq):
        return False

    ic50_uM = to_float_safe(row.get("ic50_uM"))
    if ic50_uM is None:
        return False

    status = clean_text(row.get("ic50_parse_status"))
    if mode == "strict":
        return status in STRICT_STATUSES
    if mode == "expanded":
        return status in EXPANDED_STATUSES
    raise ValueError(f"未知 mode: {mode}")


def add_workset_fields(row: Dict[str, str], tier: str, keep_flag: int, exclusion_reasons: List[str]) -> Dict[str, object]:
    out = dict(row)
    out["workset_tier"] = tier
    out["workset_keep_flag"] = keep_flag
    out["workset_exclusion_reasons"] = "; ".join(exclusion_reasons) if exclusion_reasons else ""
    return out


def build_record_level(rows: List[Dict[str, str]], mode: str) -> List[Dict[str, object]]:
    out = []
    for row in rows:
        reasons = build_exclusion_reasons(row)
        if keep_row(row, mode):
            out.append(add_workset_fields(row, mode, 1, reasons))
    return out


def build_excluded_rows(rows: List[Dict[str, str]]) -> List[Dict[str, object]]:
    out = []
    for row in rows:
        reasons = build_exclusion_reasons(row)
        if reasons:
            out.append(add_workset_fields(row, "excluded", 0, reasons))
    return out


SEQUENCE_LEVEL_COLUMNS = [
    "sequence",
    "peptide_length",
    "record_count",
    "source_record_id_list",
    "record_id_list",
    "activity_label_set",
    "food_matrix_set",
    "protein_set",
    "culture_set",
    "hydrolysis_set",
    "experiment_set",
    "model_set",
    "assay_set",
    "method_of_analysis_set",
    "pubmed_id_set",
    "title_set",
    "ic50_status_set",
    "ic50_relation_set",
    "ic50_uM_min",
    "ic50_uM_max",
    "ic50_uM_mean",
    "ic50_uM_median",
    "ic50_spread_ratio_max_min",
    "stability_flag",
    "stability_note",
    "bibliographic_year_min",
    "bibliographic_year_max",
]


def build_stability_flag(record_count: int, spread_ratio: Optional[float]) -> Tuple[str, str]:
    if record_count <= 1:
        return "singleton", "only_one_record_in_fermfoodb"
    if spread_ratio is None:
        return "unknown", "spread_ratio_unavailable"
    if spread_ratio <= 10:
        return "within_range_consistent", f"spread_ratio={spread_ratio:.4f}"
    return "conflicting", f"spread_ratio={spread_ratio:.4f}"


def build_sequence_level(rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
    groups: Dict[str, List[Dict[str, object]]] = defaultdict(list)

    for row in rows:
        seq = clean_sequence(row.get("sequence"))
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

            for field in ["title_raw", "pubmed_id_raw"]:
                y = year_from_any_text(m.get(field))
                if y is not None:
                    years.append(y)

        ic50_min = safe_min(ic50_values)
        ic50_max = safe_max(ic50_values)
        spread_ratio = None
        if ic50_min is not None and ic50_max is not None and ic50_min > 0:
            spread_ratio = ic50_max / ic50_min

        stability_flag, stability_note = build_stability_flag(
            record_count=len(members),
            spread_ratio=spread_ratio,
        )

        row = {
            "sequence": seq,
            "peptide_length": to_int_safe(members[0].get("peptide_length")),
            "record_count": len(members),
            "source_record_id_list": sorted_unique_join([m.get("source_record_id") for m in members]),
            "record_id_list": sorted_unique_join([m.get("record_id") for m in members]),
            "activity_label_set": sorted_unique_join([m.get("activity_label_raw") for m in members]),
            "food_matrix_set": sorted_unique_join([m.get("food_matrix_raw") for m in members]),
            "protein_set": sorted_unique_join([m.get("protein_raw") for m in members]),
            "culture_set": sorted_unique_join([m.get("culture_raw") for m in members]),
            "hydrolysis_set": sorted_unique_join([m.get("hydrolysis_raw") for m in members]),
            "experiment_set": sorted_unique_join([m.get("experiment_raw") for m in members]),
            "model_set": sorted_unique_join([m.get("model_raw") for m in members]),
            "assay_set": sorted_unique_join([m.get("assay_raw") for m in members]),
            "method_of_analysis_set": sorted_unique_join([m.get("method_of_analysis_raw") for m in members]),
            "pubmed_id_set": sorted_unique_join([m.get("pubmed_id_raw") for m in members]),
            "title_set": sorted_unique_join([m.get("title_raw") for m in members]),
            "ic50_status_set": sorted_unique_join([m.get("ic50_parse_status") for m in members]),
            "ic50_relation_set": sorted_unique_join([m.get("ic50_relation") for m in members]),
            "ic50_uM_min": ic50_min,
            "ic50_uM_max": ic50_max,
            "ic50_uM_mean": safe_mean(ic50_values),
            "ic50_uM_median": safe_median(ic50_values),
            "ic50_spread_ratio_max_min": spread_ratio,
            "stability_flag": stability_flag,
            "stability_note": stability_note,
            "bibliographic_year_min": min(years) if years else None,
            "bibliographic_year_max": max(years) if years else None,
        }
        out.append(row)

    out.sort(key=lambda r: (
        to_int_safe(r.get("peptide_length")) if r.get("peptide_length") is not None else 10**9,
        str(r.get("sequence") or ""),
    ))
    return out


def build_summary_rows(
    input_rows: List[Dict[str, str]],
    strict_record_rows: List[Dict[str, object]],
    strict_seq_rows: List[Dict[str, object]],
    expanded_record_rows: List[Dict[str, object]],
    expanded_seq_rows: List[Dict[str, object]],
    short23_strict_record_rows: List[Dict[str, object]],
    short23_strict_seq_rows: List[Dict[str, object]],
    short23_expanded_record_rows: List[Dict[str, object]],
    short23_expanded_seq_rows: List[Dict[str, object]],
    excluded_rows: List[Dict[str, object]],
) -> List[Dict[str, object]]:
    summary = [
        {"metric": "input_rows_total", "value": len(input_rows)},
        {"metric": "excluded_rows_total", "value": len(excluded_rows)},
        {"metric": "strict_record_level_rows", "value": len(strict_record_rows)},
        {"metric": "strict_sequence_level_rows", "value": len(strict_seq_rows)},
        {"metric": "expanded_record_level_rows", "value": len(expanded_record_rows)},
        {"metric": "expanded_sequence_level_rows", "value": len(expanded_seq_rows)},
        {"metric": "short_2_3_strict_record_level_rows", "value": len(short23_strict_record_rows)},
        {"metric": "short_2_3_strict_sequence_level_rows", "value": len(short23_strict_seq_rows)},
        {"metric": "short_2_3_expanded_record_level_rows", "value": len(short23_expanded_record_rows)},
        {"metric": "short_2_3_expanded_sequence_level_rows", "value": len(short23_expanded_seq_rows)},
    ]

    input_status_counter = Counter(
        clean_text(r.get("ic50_parse_status")) or ""
        for r in input_rows
    )
    for k, v in input_status_counter.items():
        summary.append({"metric": f"input_parse_status::{k}", "value": v})

    strict_status_counter = Counter(
        clean_text(r.get("ic50_parse_status")) or ""
        for r in strict_record_rows
    )
    for k, v in strict_status_counter.items():
        summary.append({"metric": f"strict_parse_status::{k}", "value": v})

    expanded_status_counter = Counter(
        clean_text(r.get("ic50_parse_status")) or ""
        for r in expanded_record_rows
    )
    for k, v in expanded_status_counter.items():
        summary.append({"metric": f"expanded_parse_status::{k}", "value": v})

    strict_length_counter = Counter(
        to_int_safe(r.get("peptide_length"))
        for r in strict_record_rows
        if to_int_safe(r.get("peptide_length")) is not None
    )
    for k in sorted(strict_length_counter):
        summary.append({"metric": f"strict_length::{k}", "value": strict_length_counter[k]})

    expanded_length_counter = Counter(
        to_int_safe(r.get("peptide_length"))
        for r in expanded_record_rows
        if to_int_safe(r.get("peptide_length")) is not None
    )
    for k in sorted(expanded_length_counter):
        summary.append({"metric": f"expanded_length::{k}", "value": expanded_length_counter[k]})

    strict_stability_counter = Counter(clean_text(r.get("stability_flag")) or "" for r in strict_seq_rows)
    for k, v in strict_stability_counter.items():
        summary.append({"metric": f"strict_seq_stability::{k}", "value": v})

    expanded_stability_counter = Counter(clean_text(r.get("stability_flag")) or "" for r in expanded_seq_rows)
    for k, v in expanded_stability_counter.items():
        summary.append({"metric": f"expanded_seq_stability::{k}", "value": v})

    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description="构建 FermFooDb ACE core worksets（纯标准库版）")
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
        help="输入 rescued 主表路径；默认 DB/standardized/ace/fermfoodb/fermfoodb_ace_master_rescued.csv",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="输出 worksets 目录；默认 DB/worksets/ace/fermfoodb",
    )
    args = parser.parse_args()

    script_file = Path(__file__)
    project_root = Path(args.project_root) if args.project_root else infer_project_root(script_file)

    input_path = Path(args.input_path) if args.input_path else (
        project_root / "DB" / "standardized" / "ace" / "fermfoodb" / INPUT_FILE
    )
    output_dir = Path(args.output_dir) if args.output_dir else (
        project_root / "DB" / "worksets" / "ace" / "fermfoodb"
    )

    input_rows = read_csv_rows(input_path)
    if not input_rows:
        raise ValueError(f"输入表为空：{input_path}")

    record_columns = list(input_rows[0].keys()) + [
        "workset_tier",
        "workset_keep_flag",
        "workset_exclusion_reasons",
    ]

    strict_record_rows = build_record_level(input_rows, "strict")
    strict_seq_rows = build_sequence_level(strict_record_rows)

    expanded_record_rows = build_record_level(input_rows, "expanded")
    expanded_seq_rows = build_sequence_level(expanded_record_rows)

    short23_strict_record_rows = [r for r in strict_record_rows if is_short_2_3_row(r)]
    short23_strict_seq_rows = build_sequence_level(short23_strict_record_rows)

    short23_expanded_record_rows = [r for r in expanded_record_rows if is_short_2_3_row(r)]
    short23_expanded_seq_rows = build_sequence_level(short23_expanded_record_rows)

    excluded_rows = build_excluded_rows(input_rows)

    write_csv(output_dir / STRICT_RECORD_FILE, strict_record_rows, record_columns)
    write_csv(output_dir / STRICT_SEQUENCE_FILE, strict_seq_rows, SEQUENCE_LEVEL_COLUMNS)
    write_csv(output_dir / EXPANDED_RECORD_FILE, expanded_record_rows, record_columns)
    write_csv(output_dir / EXPANDED_SEQUENCE_FILE, expanded_seq_rows, SEQUENCE_LEVEL_COLUMNS)

    write_csv(output_dir / SHORT23_STRICT_RECORD_FILE, short23_strict_record_rows, record_columns)
    write_csv(output_dir / SHORT23_STRICT_SEQUENCE_FILE, short23_strict_seq_rows, SEQUENCE_LEVEL_COLUMNS)
    write_csv(output_dir / SHORT23_EXPANDED_RECORD_FILE, short23_expanded_record_rows, record_columns)
    write_csv(output_dir / SHORT23_EXPANDED_SEQUENCE_FILE, short23_expanded_seq_rows, SEQUENCE_LEVEL_COLUMNS)

    write_csv(output_dir / EXCLUDED_FILE, excluded_rows, record_columns)

    summary_rows = build_summary_rows(
        input_rows=input_rows,
        strict_record_rows=strict_record_rows,
        strict_seq_rows=strict_seq_rows,
        expanded_record_rows=expanded_record_rows,
        expanded_seq_rows=expanded_seq_rows,
        short23_strict_record_rows=short23_strict_record_rows,
        short23_strict_seq_rows=short23_strict_seq_rows,
        short23_expanded_record_rows=short23_expanded_record_rows,
        short23_expanded_seq_rows=short23_expanded_seq_rows,
        excluded_rows=excluded_rows,
    )
    write_csv(output_dir / SUMMARY_FILE, summary_rows, ["metric", "value"])

    print("=" * 80)
    print("FermFooDb ACE core worksets 构建完成")
    print(f"项目根目录：{project_root}")
    print(f"输入表：{input_path}")
    print(f"输出目录：{output_dir}")
    print(f"输入记录数：{len(input_rows)}")
    print(f"strict record-level：{len(strict_record_rows)}")
    print(f"strict sequence-level：{len(strict_seq_rows)}")
    print(f"expanded record-level：{len(expanded_record_rows)}")
    print(f"expanded sequence-level：{len(expanded_seq_rows)}")
    print(f"short 2/3 strict record-level：{len(short23_strict_record_rows)}")
    print(f"short 2/3 strict sequence-level：{len(short23_strict_seq_rows)}")
    print(f"short 2/3 expanded record-level：{len(short23_expanded_record_rows)}")
    print(f"short 2/3 expanded sequence-level：{len(short23_expanded_seq_rows)}")
    print(f"excluded rows：{len(excluded_rows)}")
    print("=" * 80)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())