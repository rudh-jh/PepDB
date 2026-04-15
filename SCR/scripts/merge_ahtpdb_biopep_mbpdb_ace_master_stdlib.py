#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AHTPDB + BIOPEP-UWM + MBPDB ACE merged master 构建脚本（纯标准库版）
==================================================================

输入
----
1) DB/standardized/ace/ahtpdb/ahtpdb_master_clean_um_expanded.csv
2) DB/worksets/ace/biopep_uwm/biopep_uwm_ic50_core_record_level.csv
3) DB/worksets/ace/mbpdb/mbpdb_ace_core_record_level.csv

输出
----
1) DB/standardized/ace/merged/ace_master_merged_record_level_v0_2.csv
2) DB/worksets/ace/merged/ace_short_2_3_merged_record_level_v0_2.csv
3) DB/worksets/ace/merged/ace_short_2_3_merged_sequence_level_v0_2.csv
4) DB/worksets/ace/merged/ace_merge_summary_v0_2.csv

说明
----
- 当前阶段先做 merged record-level / short 2-3 record-level / short 2-3 sequence-level
- 不回写源表
- 自动把多行文本字段压成单行，避免 GitHub 预览和后续 CSV 复用混乱
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional


AHTPDB_INPUT = "ahtpdb_master_clean_um_expanded.csv"
BIOPEP_INPUT = "biopep_uwm_ic50_core_record_level.csv"
MBPDB_INPUT = "mbpdb_ace_core_record_level.csv"

MERGED_RECORD_OUTPUT = "ace_master_merged_record_level_v0_2.csv"
SHORT_RECORD_OUTPUT = "ace_short_2_3_merged_record_level_v0_2.csv"
SHORT_SEQUENCE_OUTPUT = "ace_short_2_3_merged_sequence_level_v0_2.csv"
SUMMARY_OUTPUT = "ace_merge_summary_v0_2.csv"

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

MERGED_RECORD_COLUMNS = [
    "merged_record_id",
    "task",
    "source_name",
    "source_level",
    "source_record_id",
    "sequence",
    "peptide_length",
    "sequence_valid_flag",
    "noncanonical_flag",
    "ic50_raw",
    "ic50_relation",
    "ic50_value",
    "ic50_unit",
    "ic50_uM",
    "ic50_status",
    "assay_raw",
    "source_context_raw",
    "species_raw",
    "protein_id_raw",
    "protein_description_raw",
    "bibliographic_raw",
    "doi_raw",
    "database_reference_raw",
    "is_dipeptide",
    "is_tripeptide",
    "merge_note",
]

SHORT_SEQUENCE_COLUMNS = [
    "sequence",
    "peptide_length",
    "source_count",
    "source_name_set",
    "record_count_total",
    "record_count_ahtpdb",
    "record_count_biopep",
    "record_count_mbpdb",
    "merged_record_id_list",
    "source_record_id_list",
    "species_set",
    "protein_id_set",
    "doi_set",
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


def join_nonmissing(values: List[object], sep: str = " | ") -> str:
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


# =========================================================
# AHTPDB 映射
# =========================================================
def keep_ahtpdb_row(row: Dict[str, str]) -> bool:
    sequence_valid = truthy_flag(row.get("sequence_valid"))
    noncanonical = truthy_flag(row.get("sequence_has_noncanonical_char"))
    ic50_uM = to_float_safe(row.get("ic50_uM"))
    ic50_type = clean_text(row.get("ic50_type"))

    if not sequence_valid:
        return False
    if noncanonical:
        return False
    if ic50_uM is None:
        return False
    if ic50_type not in {"exact", "threshold"}:
        return False
    return True


def map_ahtpdb_row(row: Dict[str, str], idx: int) -> Dict[str, object]:
    sequence = clean_text(row.get("sequence_clean"))
    peptide_length = to_int_safe(row.get("len_clean"))
    ic50_type = clean_text(row.get("ic50_type"))

    if ic50_type == "exact":
        ic50_status = "exact_molar"
    elif ic50_type == "threshold":
        ic50_status = "molar_threshold"
    else:
        ic50_status = ic50_type

    assay_parts = [
        row.get("method_raw"),
        row.get("assay_raw"),
        row.get("mice_raw"),
    ]

    source_context_parts = [
        row.get("source_raw"),
        row.get("table_from"),
    ]

    return {
        "merged_record_id": f"ace_merged_{idx:06d}",
        "task": "ace",
        "source_name": "AHTPDB",
        "source_level": "ahtpdb_standardized",
        "source_record_id": clean_text(row.get("id")),
        "sequence": sequence,
        "peptide_length": peptide_length,
        "sequence_valid_flag": 1 if truthy_flag(row.get("sequence_valid")) else 0,
        "noncanonical_flag": 1 if truthy_flag(row.get("sequence_has_noncanonical_char")) else 0,
        "ic50_raw": clean_text(row.get("ic50_raw")),
        "ic50_relation": clean_text(row.get("ic50_relation")),
        "ic50_value": to_float_safe(row.get("ic50_value")),
        "ic50_unit": clean_text(row.get("ic50_unit_raw")),
        "ic50_uM": to_float_safe(row.get("ic50_uM")),
        "ic50_status": ic50_status,
        "assay_raw": join_nonmissing(assay_parts),
        "source_context_raw": join_nonmissing(source_context_parts),
        "species_raw": "",
        "protein_id_raw": "",
        "protein_description_raw": "",
        "bibliographic_raw": "",
        "doi_raw": "",
        "database_reference_raw": "",
        "is_dipeptide": 1 if truthy_flag(row.get("is_dipeptide")) else 0,
        "is_tripeptide": 1 if truthy_flag(row.get("is_tripeptide")) else 0,
        "merge_note": "mapped_from_ahtpdb_master_clean_um_expanded",
    }


# =========================================================
# BIOPEP 映射
# =========================================================
def keep_biopep_row(row: Dict[str, str]) -> bool:
    sequence = clean_text(row.get("sequence"))
    peptide_length = to_int_safe(row.get("peptide_length"))
    ic50_uM = to_float_safe(row.get("ic50_uM"))
    ic50_status = clean_text(row.get("ic50_parse_status"))

    if sequence is None:
        return False
    if peptide_length is None:
        return False
    if ic50_uM is None:
        return False
    if ic50_status not in {"exact_molar", "molar_threshold"}:
        return False
    return True


def map_biopep_row(row: Dict[str, str], idx: int) -> Dict[str, object]:
    peptide_length = to_int_safe(row.get("peptide_length"))

    return {
        "merged_record_id": f"ace_merged_{idx:06d}",
        "task": "ace",
        "source_name": "BIOPEP-UWM",
        "source_level": "biopep_core_workset",
        "source_record_id": clean_text(row.get("source_record_id")),
        "sequence": clean_text(row.get("sequence")),
        "peptide_length": peptide_length,
        "sequence_valid_flag": 1,
        "noncanonical_flag": 0,
        "ic50_raw": clean_text(row.get("ic50_raw")),
        "ic50_relation": clean_text(row.get("ic50_relation")),
        "ic50_value": to_float_safe(row.get("ic50_value")),
        "ic50_unit": clean_text(row.get("ic50_unit")),
        "ic50_uM": to_float_safe(row.get("ic50_uM")),
        "ic50_status": clean_text(row.get("ic50_parse_status")),
        "assay_raw": join_nonmissing([
            row.get("activity_measure_label_raw"),
            row.get("measure_type_raw"),
            row.get("function_text_raw"),
        ]),
        "source_context_raw": join_nonmissing([
            row.get("source_name"),
            row.get("source_type"),
            row.get("target"),
            row.get("evidence_type"),
        ]),
        "species_raw": "",
        "protein_id_raw": "",
        "protein_description_raw": "",
        "bibliographic_raw": clean_text(row.get("bibliographic_raw")) or "",
        "doi_raw": "",
        "database_reference_raw": clean_text(row.get("database_reference_raw")) or "",
        "is_dipeptide": 1 if peptide_length == 2 else 0,
        "is_tripeptide": 1 if peptide_length == 3 else 0,
        "merge_note": "mapped_from_biopep_ic50_core_record_level",
    }


# =========================================================
# MBPDB 映射
# =========================================================
def keep_mbpdb_row(row: Dict[str, str]) -> bool:
    sequence = clean_text(row.get("sequence"))
    peptide_length = to_int_safe(row.get("peptide_length"))
    ic50_uM = to_float_safe(row.get("ic50_uM"))
    ic50_status = clean_text(row.get("ic50_parse_status"))

    if sequence is None:
        return False
    if peptide_length is None:
        return False
    if ic50_uM is None:
        return False
    if ic50_status not in {"exact_molar", "molar_threshold"}:
        return False
    return True


def map_mbpdb_row(row: Dict[str, str], idx: int) -> Dict[str, object]:
    peptide_length = to_int_safe(row.get("peptide_length"))
    source_record_id = clean_text(row.get("source_record_id"))
    species_raw = clean_text(row.get("species_raw"))
    protein_id_raw = clean_text(row.get("protein_id_raw"))
    protein_description_raw = clean_text(row.get("protein_description_raw"))

    return {
        "merged_record_id": f"ace_merged_{idx:06d}",
        "task": "ace",
        "source_name": "MBPDB",
        "source_level": "mbpdb_core_workset",
        "source_record_id": source_record_id,
        "sequence": clean_text(row.get("sequence")),
        "peptide_length": peptide_length,
        "sequence_valid_flag": 1,
        "noncanonical_flag": 0,
        "ic50_raw": clean_text(row.get("ic50_raw")),
        "ic50_relation": clean_text(row.get("ic50_relation")),
        "ic50_value": to_float_safe(row.get("ic50_value")),
        "ic50_unit": clean_text(row.get("ic50_unit")),
        "ic50_uM": to_float_safe(row.get("ic50_uM")),
        "ic50_status": clean_text(row.get("ic50_parse_status")),
        "assay_raw": join_nonmissing([
            row.get("function_raw"),
            row.get("inhibition_type_raw"),
            row.get("additional_details_raw"),
        ]),
        "source_context_raw": join_nonmissing([
            row.get("source_name"),
            row.get("source_type"),
            row.get("target"),
            row.get("intervals_raw"),
        ]),
        "species_raw": species_raw or "",
        "protein_id_raw": protein_id_raw or "",
        "protein_description_raw": protein_description_raw or "",
        "bibliographic_raw": join_nonmissing([
            row.get("title_raw"),
            row.get("authors_raw"),
            row.get("abstract_raw"),
        ]),
        "doi_raw": clean_text(row.get("doi_raw")) or "",
        "database_reference_raw": "",
        "is_dipeptide": 1 if peptide_length == 2 else 0,
        "is_tripeptide": 1 if peptide_length == 3 else 0,
        "merge_note": "mapped_from_mbpdb_ace_core_record_level",
    }


# =========================================================
# sequence-level 聚合
# =========================================================
def build_stability_flag(source_count: int, record_count: int, spread_ratio: Optional[float]) -> (str, str):
    if record_count <= 1:
        return "singleton", "only_one_record_after_merge"
    if spread_ratio is None:
        return "unknown", "spread_ratio_unavailable"
    if source_count >= 2 and spread_ratio <= 10:
        return "cross_source_consistent", f"source_count={source_count}; spread_ratio={spread_ratio:.4f}"
    if spread_ratio <= 10:
        return "within_range_consistent", f"source_count={source_count}; spread_ratio={spread_ratio:.4f}"
    return "conflicting", f"source_count={source_count}; spread_ratio={spread_ratio:.4f}"


def build_short_sequence_level(rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
    groups: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    for row in rows:
        seq = clean_text(row.get("sequence"))
        if seq is None:
            continue
        groups[seq].append(row)

    out = []
    for seq, members in groups.items():
        source_names = sorted({
            clean_text(m.get("source_name")) for m in members
            if clean_text(m.get("source_name")) is not None
        })
        source_count = len(source_names)

        ic50_values = [
            to_float_safe(m.get("ic50_uM"))
            for m in members
            if to_float_safe(m.get("ic50_uM")) is not None
        ]
        ic50_min = safe_min(ic50_values)
        ic50_max = safe_max(ic50_values)
        spread_ratio = None
        if ic50_min is not None and ic50_max is not None and ic50_min > 0:
            spread_ratio = ic50_max / ic50_min

        years = []
        for m in members:
            for field in ["bibliographic_raw", "doi_raw", "source_context_raw"]:
                y = year_from_any_text(m.get(field))
                if y is not None:
                    years.append(y)

        stability_flag, stability_note = build_stability_flag(
            source_count=source_count,
            record_count=len(members),
            spread_ratio=spread_ratio,
        )

        out.append({
            "sequence": seq,
            "peptide_length": to_int_safe(members[0].get("peptide_length")),
            "source_count": source_count,
            "source_name_set": join_nonmissing([m.get("source_name") for m in members]),
            "record_count_total": len(members),
            "record_count_ahtpdb": sum(1 for m in members if clean_text(m.get("source_name")) == "AHTPDB"),
            "record_count_biopep": sum(1 for m in members if clean_text(m.get("source_name")) == "BIOPEP-UWM"),
            "record_count_mbpdb": sum(1 for m in members if clean_text(m.get("source_name")) == "MBPDB"),
            "merged_record_id_list": join_nonmissing([m.get("merged_record_id") for m in members]),
            "source_record_id_list": join_nonmissing([m.get("source_record_id") for m in members]),
            "species_set": join_nonmissing([m.get("species_raw") for m in members]),
            "protein_id_set": join_nonmissing([m.get("protein_id_raw") for m in members]),
            "doi_set": join_nonmissing([m.get("doi_raw") for m in members]),
            "ic50_relation_set": join_nonmissing([m.get("ic50_relation") for m in members]),
            "ic50_status_set": join_nonmissing([m.get("ic50_status") for m in members]),
            "ic50_uM_min": ic50_min,
            "ic50_uM_max": ic50_max,
            "ic50_uM_mean": safe_mean(ic50_values),
            "ic50_uM_median": safe_median(ic50_values),
            "ic50_spread_ratio_max_min": spread_ratio,
            "cross_source_overlap_flag": 1 if source_count >= 2 else 0,
            "stability_flag": stability_flag,
            "stability_note": stability_note,
            "bibliographic_year_min": min(years) if years else None,
            "bibliographic_year_max": max(years) if years else None,
            "database_reference_raw_set": join_nonmissing([m.get("database_reference_raw") for m in members]),
        })

    out.sort(key=lambda r: (
        to_int_safe(r.get("peptide_length")) if r.get("peptide_length") is not None else 10**9,
        str(r.get("sequence") or "")
    ))
    return out


def build_summary_rows(
    ahtpdb_input_count: int,
    biopep_input_count: int,
    mbpdb_input_count: int,
    ahtpdb_kept_count: int,
    biopep_kept_count: int,
    mbpdb_kept_count: int,
    merged_rows: List[Dict[str, object]],
    short_rows: List[Dict[str, object]],
    short_seq_rows: List[Dict[str, object]],
) -> List[Dict[str, object]]:
    summary = [
        {"metric": "ahtpdb_input_rows", "value": ahtpdb_input_count},
        {"metric": "biopep_input_rows", "value": biopep_input_count},
        {"metric": "mbpdb_input_rows", "value": mbpdb_input_count},
        {"metric": "ahtpdb_kept_rows", "value": ahtpdb_kept_count},
        {"metric": "biopep_kept_rows", "value": biopep_kept_count},
        {"metric": "mbpdb_kept_rows", "value": mbpdb_kept_count},
        {"metric": "merged_record_level_rows", "value": len(merged_rows)},
        {"metric": "merged_short_2_3_record_level_rows", "value": len(short_rows)},
        {"metric": "merged_short_2_3_sequence_level_rows", "value": len(short_seq_rows)},
    ]

    source_counter = Counter(clean_text(r.get("source_name")) or "" for r in merged_rows)
    for k, v in source_counter.items():
        summary.append({"metric": f"merged_source::{k}", "value": v})

    short_source_counter = Counter(clean_text(r.get("source_name")) or "" for r in short_rows)
    for k, v in short_source_counter.items():
        summary.append({"metric": f"short_source::{k}", "value": v})

    stability_counter = Counter(clean_text(r.get("stability_flag")) or "" for r in short_seq_rows)
    for k, v in stability_counter.items():
        summary.append({"metric": f"short_seq_stability::{k}", "value": v})

    overlap_counter = Counter(int(bool(to_int_safe(r.get("cross_source_overlap_flag")))) for r in short_seq_rows)
    for k, v in overlap_counter.items():
        metric = "short_seq_cross_source_overlap_yes" if k == 1 else "short_seq_cross_source_overlap_no"
        summary.append({"metric": metric, "value": v})

    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description="Merge AHTPDB + BIOPEP-UWM + MBPDB ACE master（纯标准库版）")
    parser.add_argument(
        "--project-root",
        type=str,
        default=None,
        help="PepDB 项目根目录，例如 E:\\MYS\\PepDB",
    )
    parser.add_argument(
        "--ahtpdb-input",
        type=str,
        default=None,
        help="AHTPDB 输入路径；默认 DB/standardized/ace/ahtpdb/ahtpdb_master_clean_um_expanded.csv",
    )
    parser.add_argument(
        "--biopep-input",
        type=str,
        default=None,
        help="BIOPEP 输入路径；默认 DB/worksets/ace/biopep_uwm/biopep_uwm_ic50_core_record_level.csv",
    )
    parser.add_argument(
        "--mbpdb-input",
        type=str,
        default=None,
        help="MBPDB 输入路径；默认 DB/worksets/ace/mbpdb/mbpdb_ace_core_record_level.csv",
    )
    args = parser.parse_args()

    script_file = Path(__file__)
    project_root = Path(args.project_root) if args.project_root else infer_project_root(script_file)

    ahtpdb_input = Path(args.ahtpdb_input) if args.ahtpdb_input else (
        project_root / "DB" / "standardized" / "ace" / "ahtpdb" / AHTPDB_INPUT
    )
    biopep_input = Path(args.biopep_input) if args.biopep_input else (
        project_root / "DB" / "worksets" / "ace" / "biopep_uwm" / BIOPEP_INPUT
    )
    mbpdb_input = Path(args.mbpdb_input) if args.mbpdb_input else (
        project_root / "DB" / "worksets" / "ace" / "mbpdb" / MBPDB_INPUT
    )

    merged_standardized_dir = project_root / "DB" / "standardized" / "ace" / "merged"
    merged_workset_dir = project_root / "DB" / "worksets" / "ace" / "merged"

    ahtpdb_rows = read_csv_rows(ahtpdb_input)
    biopep_rows = read_csv_rows(biopep_input)
    mbpdb_rows = read_csv_rows(mbpdb_input)

    merged_rows: List[Dict[str, object]] = []
    merged_idx = 1

    ahtpdb_kept_count = 0
    for row in ahtpdb_rows:
        if keep_ahtpdb_row(row):
            merged_rows.append(map_ahtpdb_row(row, merged_idx))
            merged_idx += 1
            ahtpdb_kept_count += 1

    biopep_kept_count = 0
    for row in biopep_rows:
        if keep_biopep_row(row):
            merged_rows.append(map_biopep_row(row, merged_idx))
            merged_idx += 1
            biopep_kept_count += 1

    mbpdb_kept_count = 0
    for row in mbpdb_rows:
        if keep_mbpdb_row(row):
            merged_rows.append(map_mbpdb_row(row, merged_idx))
            merged_idx += 1
            mbpdb_kept_count += 1

    merged_rows.sort(key=lambda r: (
        to_int_safe(r.get("peptide_length")) if r.get("peptide_length") is not None else 10**9,
        str(r.get("sequence") or ""),
        str(r.get("source_name") or ""),
        str(r.get("source_record_id") or ""),
    ))

    for i, row in enumerate(merged_rows, start=1):
        row["merged_record_id"] = f"ace_merged_{i:06d}"

    short_rows = [
        row for row in merged_rows
        if to_int_safe(row.get("peptide_length")) in {2, 3}
    ]

    short_seq_rows = build_short_sequence_level(short_rows)

    summary_rows = build_summary_rows(
        ahtpdb_input_count=len(ahtpdb_rows),
        biopep_input_count=len(biopep_rows),
        mbpdb_input_count=len(mbpdb_rows),
        ahtpdb_kept_count=ahtpdb_kept_count,
        biopep_kept_count=biopep_kept_count,
        mbpdb_kept_count=mbpdb_kept_count,
        merged_rows=merged_rows,
        short_rows=short_rows,
        short_seq_rows=short_seq_rows,
    )

    write_csv(merged_standardized_dir / MERGED_RECORD_OUTPUT, merged_rows, MERGED_RECORD_COLUMNS)
    write_csv(merged_workset_dir / SHORT_RECORD_OUTPUT, short_rows, MERGED_RECORD_COLUMNS)
    write_csv(merged_workset_dir / SHORT_SEQUENCE_OUTPUT, short_seq_rows, SHORT_SEQUENCE_COLUMNS)
    write_csv(merged_workset_dir / SUMMARY_OUTPUT, summary_rows, ["metric", "value"])

    print("=" * 80)
    print("AHTPDB + BIOPEP-UWM + MBPDB ACE merged master 构建完成")
    print(f"项目根目录：{project_root}")
    print(f"AHTPDB 输入：{ahtpdb_input}")
    print(f"BIOPEP 输入：{biopep_input}")
    print(f"MBPDB 输入：{mbpdb_input}")
    print(f"AHTPDB 输入行数：{len(ahtpdb_rows)}；保留：{ahtpdb_kept_count}")
    print(f"BIOPEP 输入行数：{len(biopep_rows)}；保留：{biopep_kept_count}")
    print(f"MBPDB 输入行数：{len(mbpdb_rows)}；保留：{mbpdb_kept_count}")
    print(f"merged record-level：{len(merged_rows)}")
    print(f"merged short 2/3 record-level：{len(short_rows)}")
    print(f"merged short 2/3 sequence-level：{len(short_seq_rows)}")
    print(f"输出：{merged_standardized_dir / MERGED_RECORD_OUTPUT}")
    print(f"输出：{merged_workset_dir / SHORT_RECORD_OUTPUT}")
    print(f"输出：{merged_workset_dir / SHORT_SEQUENCE_OUTPUT}")
    print(f"输出：{merged_workset_dir / SUMMARY_OUTPUT}")
    print("=" * 80)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())