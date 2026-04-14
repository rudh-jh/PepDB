#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BIOPEP-UWM experimental ACE 原始表标准化脚本（纯标准库版）
====================================================

本脚本面向 PepDB 当前阶段的 BIOPEP-UWM ACE raw -> standardized 过渡：
1. 读取 raw_tables_csv 下的两张原始表：
   - biopep_uwm_experimental_ace_list.csv
   - biopep_uwm_experimental_ace_detail.csv
2. 以 peptide_id 为键做合并
3. 输出一张 BIOPEP 标准化中间主表
4. 输出解析日志与若干基础检查表

设计原则
--------
- 只做“工程化标准化”，不冒进做复杂生物学判定
- 原始字段必须保留
- 自动解析字段单独输出
- 当前默认不强行做单位换算
- 当前默认不做跨库去重，也不与 AHTPDB 合并
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Tuple

LIST_FILE = "biopep_uwm_experimental_ace_list.csv"
DETAIL_FILE = "biopep_uwm_experimental_ace_detail.csv"

LIST_REQUIRED_COLUMNS = [
    "list_page_number",
    "list_page_url",
    "list_row_index",
    "peptide_id",
    "name",
    "sequence",
    "chemical_mass",
    "monoisotopic_mass",
    "activity_value_raw",
    "measure_type_raw",
    "detail_report_url",
]

DETAIL_REQUIRED_COLUMNS = [
    "peptide_id",
    "detail_report_url",
    "name",
    "sequence",
    "function_text",
    "number_of_residues",
    "activity_code",
    "activity_name",
    "chemical_mass",
    "monoisotopic_mass",
    "activity_measure_label",
    "activity_value_text",
    "bibliographic_authors",
    "bibliographic_title",
    "bibliographic_year",
    "bibliographic_source_type",
    "bibliographic_raw",
    "additional_information_raw",
    "smiles",
    "inchi",
    "inchikey",
    "database_reference_raw",
    "raw_text_path",
    "raw_html_path",
]

MASTER_COLUMNS = [
    "record_id",
    "task",
    "sequence",
    "peptide_length",
    "source_type",
    "source_name",
    "source_record_id",
    "evidence_type",
    "target",
    "activity_label_raw",
    "ic50_raw",
    "ic50_value",
    "ic50_unit",
    "ic50_relation",
    "ic50_uM",
    "ic50_parse_status",
    "ic50_parse_note",
    "name_raw",
    "function_text_raw",
    "activity_code_raw",
    "activity_name_raw",
    "activity_measure_label_raw",
    "measure_type_raw",
    "number_of_residues_raw",
    "chemical_mass",
    "monoisotopic_mass",
    "bibliographic_authors_raw",
    "bibliographic_title_raw",
    "bibliographic_year_raw",
    "bibliographic_source_type_raw",
    "bibliographic_raw",
    "database_reference_raw",
    "additional_information_raw",
    "smiles",
    "inchi",
    "inchikey",
    "detail_report_url",
    "list_page_number",
    "list_page_url",
    "list_row_index",
    "raw_text_path",
    "raw_html_path",
    "sequence_from",
    "length_consistency_flag",
    "merge_source",
    "notes",
]

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
MOLAR_UNIT_TO_UM_FACTOR = {
    "nm": 1e-3,
    "um": 1.0,
    "mm": 1e3,
    "m": 1e6,
}


# =========================
# 通用辅助函数
# =========================
def infer_project_root(script_file: Path) -> Path:
    try:
        return script_file.resolve().parents[2]
    except IndexError:
        return script_file.resolve().parent


def normalize_spaces(text: str) -> str:
    return re.sub(r"\s+", " ", text.strip())


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


def to_float_safe(value: object) -> Optional[float]:
    text = clean_text(value)
    if text is None:
        return None
    text = text.replace(",", "")
    try:
        return float(text)
    except ValueError:
        return None


def read_csv_rows(path: Path, required_columns: List[str]) -> List[Dict[str, Optional[str]]]:
    if not path.exists():
        raise FileNotFoundError(f"找不到输入文件：{path}")

    with open(path, "r", encoding="utf-8-sig", errors="replace", newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError(f"CSV 无表头：{path}")

        fieldnames = [strip_bom(x) for x in reader.fieldnames]
        missing = [c for c in required_columns if c not in fieldnames]
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


def write_csv(path: Path, rows: List[Dict[str, object]], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def clean_sequence(seq: object) -> Optional[str]:
    text = clean_text(seq)
    if text is None:
        return None
    text = text.upper().replace(" ", "")
    return text if text else None


def is_valid_canonical_sequence(seq: Optional[str]) -> bool:
    if not seq:
        return False
    return set(seq).issubset(CANONICAL_AA)


def compute_length(seq: Optional[str]) -> Optional[int]:
    if not seq:
        return None
    return len(seq)


def choose_first(*values: object) -> Optional[str]:
    for value in values:
        text = clean_text(value)
        if text is not None:
            return text
    return None


def normalize_unit_text(text: str) -> str:
    text = text.strip()
    text = text.replace("μ", "u").replace("µ", "u")
    text = text.replace("UM", "uM").replace("NM", "nM").replace("MM", "mM")
    text = text.replace("ug/ml", "ug/mL").replace("UG/ML", "ug/mL")
    text = text.replace("mg/ml", "mg/mL").replace("MG/ML", "mg/mL")
    text = text.replace("ng/ml", "ng/mL").replace("NG/ML", "ng/mL")
    return text


def compose_activity_label_raw(list_row: Dict[str, Optional[str]], detail_row: Dict[str, Optional[str]]) -> Optional[str]:
    parts = [
        detail_row.get("activity_name"),
        list_row.get("name"),
        detail_row.get("function_text"),
    ]
    seen = []
    for x in parts:
        x = clean_text(x)
        if x and x not in seen:
            seen.append(x)
    if not seen:
        return None
    return " | ".join(seen)


def compose_ic50_raw(list_row: Dict[str, Optional[str]], detail_row: Dict[str, Optional[str]]) -> Optional[str]:
    detail_value = clean_text(detail_row.get("activity_value_text"))
    detail_measure = clean_text(detail_row.get("activity_measure_label"))
    list_value = clean_text(list_row.get("activity_value_raw"))
    list_measure = clean_text(list_row.get("measure_type_raw"))

    if detail_value and detail_measure:
        return f"{detail_value} [{detail_measure}]"
    if detail_value:
        return detail_value
    if list_value and list_measure:
        return f"{list_value} [{list_measure}]"
    if list_value:
        return list_value
    return None


def parse_ic50(raw_value: object) -> Dict[str, object]:
    result = {
        "ic50_value": None,
        "ic50_unit": None,
        "ic50_relation": None,
        "ic50_uM": None,
        "ic50_parse_status": "missing",
        "ic50_parse_note": None,
    }

    text = clean_text(raw_value)
    if text is None:
        return result

    work_text = re.sub(r"\s*\[[^\]]+\]\s*$", "", text)
    work_text = normalize_unit_text(work_text)

    m_plain = re.fullmatch(r"([0-9]+(?:\.[0-9]+)?)", work_text)
    if m_plain:
        result["ic50_value"] = float(m_plain.group(1))
        result["ic50_relation"] = "="
        result["ic50_parse_status"] = "numeric_no_unit"
        result["ic50_parse_note"] = "numeric_value_present_but_unit_missing"
        return result

    if "%" in work_text:
        m_percent = re.search(r"([0-9]+(?:\.[0-9]+)?)", work_text)
        if m_percent:
            result["ic50_value"] = float(m_percent.group(1))
        result["ic50_parse_status"] = "percent_like"
        result["ic50_parse_note"] = "percent_or_percent_like_text"
        return result

    m_range = re.fullmatch(
        r"([0-9]+(?:\.[0-9]+)?)\s*[-~]\s*([0-9]+(?:\.[0-9]+)?)\s*(nM|uM|mM|M)",
        work_text,
    )
    if m_range:
        result["ic50_parse_status"] = "range"
        result["ic50_parse_note"] = "range_value_not_collapsed_to_single_number"
        result["ic50_unit"] = m_range.group(3)
        return result

    m_mass = re.fullmatch(r"([<>≤≥])?\s*([0-9]+(?:\.[0-9]+)?)\s*(ng/mL|ug/mL|mg/mL)", work_text)
    if m_mass:
        result["ic50_value"] = float(m_mass.group(2))
        result["ic50_relation"] = m_mass.group(1) if m_mass.group(1) else "="
        result["ic50_unit"] = m_mass.group(3)
        result["ic50_parse_status"] = "mass_unit"
        result["ic50_parse_note"] = "mass_concentration_kept_without_uM_conversion"
        return result

    m_threshold = re.fullmatch(r"([<>≤≥])\s*([0-9]+(?:\.[0-9]+)?)\s*(nM|uM|mM|M)", work_text)
    if m_threshold:
        relation = m_threshold.group(1)
        value = float(m_threshold.group(2))
        unit = m_threshold.group(3)
        factor = MOLAR_UNIT_TO_UM_FACTOR[unit.lower()]
        result["ic50_value"] = value
        result["ic50_unit"] = unit
        result["ic50_relation"] = relation
        result["ic50_uM"] = value * factor
        result["ic50_parse_status"] = "molar_threshold"
        result["ic50_parse_note"] = "molar_threshold_converted_to_uM"
        return result

    m_exact = re.fullmatch(r"([0-9]+(?:\.[0-9]+)?)\s*(nM|uM|mM|M)", work_text)
    if m_exact:
        value = float(m_exact.group(1))
        unit = m_exact.group(2)
        factor = MOLAR_UNIT_TO_UM_FACTOR[unit.lower()]
        result["ic50_value"] = value
        result["ic50_unit"] = unit
        result["ic50_relation"] = "="
        result["ic50_uM"] = value * factor
        result["ic50_parse_status"] = "exact_molar"
        result["ic50_parse_note"] = "exact_molar_value_converted_to_uM"
        return result

    m_any = re.search(r"([0-9]+(?:\.[0-9]+)?)", work_text)
    if m_any:
        result["ic50_value"] = float(m_any.group(1))
        result["ic50_parse_status"] = "ambiguous_text"
        result["ic50_parse_note"] = "number_detected_but_not_safely_normalized"
        return result

    result["ic50_parse_status"] = "unparsed_text"
    result["ic50_parse_note"] = "text_present_but_no_safe_numeric_parse"
    return result


def build_merge_source(list_row: Optional[Dict[str, Optional[str]]], detail_row: Optional[Dict[str, Optional[str]]]) -> str:
    if list_row and detail_row:
        return "list+detail"
    if list_row:
        return "list_only"
    if detail_row:
        return "detail_only"
    return "none"


def build_master_row(
    peptide_id: str,
    idx: int,
    list_row: Optional[Dict[str, Optional[str]]],
    detail_row: Optional[Dict[str, Optional[str]]],
) -> Tuple[Dict[str, object], Dict[str, object]]:
    list_row = list_row or {}
    detail_row = detail_row or {}

    seq_list = clean_sequence(list_row.get("sequence"))
    seq_detail = clean_sequence(detail_row.get("sequence"))
    sequence = seq_list or seq_detail
    sequence_from = "list.sequence" if seq_list else ("detail.sequence" if seq_detail else "missing")

    seq_valid = is_valid_canonical_sequence(sequence)
    seq_note = None
    if sequence and not seq_valid:
        seq_note = "non_canonical_or_contains_unknown_residue"

    peptide_length = compute_length(sequence)
    number_of_residues_raw = clean_text(detail_row.get("number_of_residues"))
    residues_from_detail = to_float_safe(number_of_residues_raw)
    if residues_from_detail is not None:
        residues_from_detail = int(residues_from_detail)

    if peptide_length is None and residues_from_detail is not None:
        peptide_length = residues_from_detail

    if peptide_length is None and residues_from_detail is None:
        length_consistency_flag = "missing_both"
    elif peptide_length is not None and residues_from_detail is None:
        length_consistency_flag = "sequence_only"
    elif peptide_length is None and residues_from_detail is not None:
        length_consistency_flag = "detail_only"
    elif peptide_length == residues_from_detail:
        length_consistency_flag = "consistent"
    else:
        length_consistency_flag = "conflict"

    ic50_raw = compose_ic50_raw(list_row, detail_row)
    ic50_parsed = parse_ic50(ic50_raw)

    chemical_mass = choose_first(list_row.get("chemical_mass"), detail_row.get("chemical_mass"))
    monoisotopic_mass = choose_first(list_row.get("monoisotopic_mass"), detail_row.get("monoisotopic_mass"))
    detail_report_url = choose_first(detail_row.get("detail_report_url"), list_row.get("detail_report_url"))

    notes_list = []
    if seq_note:
        notes_list.append(seq_note)
    if sequence is None:
        notes_list.append("sequence_missing")
    if length_consistency_flag == "conflict":
        notes_list.append("sequence_length_conflicts_with_detail_number_of_residues")
    if ic50_parsed["ic50_parse_status"] in {"numeric_no_unit", "ambiguous_text", "unparsed_text"}:
        notes_list.append("ic50_needs_manual_review")
    if not clean_text(detail_row.get("bibliographic_raw")):
        notes_list.append("bibliography_missing")
    notes = "; ".join(notes_list) if notes_list else None

    master_row: Dict[str, object] = {
        "record_id": f"biopep_ace_{idx:06d}",
        "task": "ace",
        "sequence": sequence,
        "peptide_length": peptide_length,
        "source_type": "database",
        "source_name": "BIOPEP-UWM",
        "source_record_id": peptide_id,
        "evidence_type": "experimental",
        "target": "ACE",
        "activity_label_raw": compose_activity_label_raw(list_row, detail_row),
        "ic50_raw": ic50_raw,
        "ic50_value": ic50_parsed["ic50_value"],
        "ic50_unit": ic50_parsed["ic50_unit"],
        "ic50_relation": ic50_parsed["ic50_relation"],
        "ic50_uM": ic50_parsed["ic50_uM"],
        "ic50_parse_status": ic50_parsed["ic50_parse_status"],
        "ic50_parse_note": ic50_parsed["ic50_parse_note"],
        "name_raw": choose_first(detail_row.get("name"), list_row.get("name")),
        "function_text_raw": detail_row.get("function_text"),
        "activity_code_raw": detail_row.get("activity_code"),
        "activity_name_raw": detail_row.get("activity_name"),
        "activity_measure_label_raw": detail_row.get("activity_measure_label"),
        "measure_type_raw": list_row.get("measure_type_raw"),
        "number_of_residues_raw": number_of_residues_raw,
        "chemical_mass": chemical_mass,
        "monoisotopic_mass": monoisotopic_mass,
        "bibliographic_authors_raw": detail_row.get("bibliographic_authors"),
        "bibliographic_title_raw": detail_row.get("bibliographic_title"),
        "bibliographic_year_raw": detail_row.get("bibliographic_year"),
        "bibliographic_source_type_raw": detail_row.get("bibliographic_source_type"),
        "bibliographic_raw": detail_row.get("bibliographic_raw"),
        "database_reference_raw": detail_row.get("database_reference_raw"),
        "additional_information_raw": detail_row.get("additional_information_raw"),
        "smiles": detail_row.get("smiles"),
        "inchi": detail_row.get("inchi"),
        "inchikey": detail_row.get("inchikey"),
        "detail_report_url": detail_report_url,
        "list_page_number": list_row.get("list_page_number"),
        "list_page_url": list_row.get("list_page_url"),
        "list_row_index": list_row.get("list_row_index"),
        "raw_text_path": detail_row.get("raw_text_path"),
        "raw_html_path": detail_row.get("raw_html_path"),
        "sequence_from": sequence_from,
        "length_consistency_flag": length_consistency_flag,
        "merge_source": build_merge_source(
            list_row if list_row else None,
            detail_row if detail_row else None,
        ),
        "notes": notes,
    }

    parse_log: Dict[str, object] = {
        "record_id": master_row["record_id"],
        "source_record_id": peptide_id,
        "ic50_raw": ic50_raw,
        "ic50_value": ic50_parsed["ic50_value"],
        "ic50_unit": ic50_parsed["ic50_unit"],
        "ic50_relation": ic50_parsed["ic50_relation"],
        "ic50_uM": ic50_parsed["ic50_uM"],
        "ic50_parse_status": ic50_parsed["ic50_parse_status"],
        "ic50_parse_note": ic50_parsed["ic50_parse_note"],
        "sequence": sequence,
        "sequence_from": sequence_from,
        "peptide_length": peptide_length,
        "number_of_residues_raw": number_of_residues_raw,
        "length_consistency_flag": length_consistency_flag,
        "notes": notes,
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


def build_counter_table(counter: Counter, key_name: str, value_name: str) -> List[Dict[str, object]]:
    return [{key_name: k, value_name: v} for k, v in counter.most_common()]


def build_duplicate_id_table(master_rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
    counter = Counter(str(x.get("source_record_id") or "") for x in master_rows)
    out = []
    for pid, cnt in counter.most_common():
        if pid and cnt > 1:
            out.append({"source_record_id": pid, "duplicate_count": cnt})
    return out


def build_sequence_duplicate_table(master_rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
    counter = Counter(str(x.get("sequence") or "") for x in master_rows if x.get("sequence"))
    out = []
    for seq, cnt in counter.most_common():
        if seq and cnt > 1:
            out.append({"sequence": seq, "duplicate_count": cnt})
    return out


def build_summary_rows(
    list_rows: List[Dict[str, Optional[str]]],
    detail_rows: List[Dict[str, Optional[str]]],
    master_rows: List[Dict[str, object]],
) -> List[Dict[str, object]]:
    merge_counter = Counter(str(x.get("merge_source") or "") for x in master_rows)
    parse_counter = Counter(str(x.get("ic50_parse_status") or "") for x in master_rows)
    length_counter = Counter(str(x.get("length_consistency_flag") or "") for x in master_rows)

    summary = [
        {"metric": "list_rows", "value": len(list_rows)},
        {"metric": "detail_rows", "value": len(detail_rows)},
        {"metric": "master_rows", "value": len(master_rows)},
    ]
    for k, v in merge_counter.items():
        summary.append({"metric": f"merge_source::{k}", "value": v})
    for k, v in parse_counter.items():
        summary.append({"metric": f"ic50_parse_status::{k}", "value": v})
    for k, v in length_counter.items():
        summary.append({"metric": f"length_consistency::{k}", "value": v})
    return summary


def standardize(list_rows: List[Dict[str, Optional[str]]], detail_rows: List[Dict[str, Optional[str]]]):
    list_map = {}
    for row in list_rows:
        pid = clean_text(row.get("peptide_id"))
        if pid:
            list_map[pid] = row

    detail_map = {}
    for row in detail_rows:
        pid = clean_text(row.get("peptide_id"))
        if pid:
            detail_map[pid] = row

    all_ids = sorted(set(list_map.keys()) | set(detail_map.keys()), key=lambda x: int(x) if x.isdigit() else x)

    master_rows: List[Dict[str, object]] = []
    parse_logs: List[Dict[str, object]] = []
    for idx, peptide_id in enumerate(all_ids, start=1):
        master_row, parse_log = build_master_row(
            peptide_id=peptide_id,
            idx=idx,
            list_row=list_map.get(peptide_id),
            detail_row=detail_map.get(peptide_id),
        )
        master_rows.append(master_row)
        parse_logs.append(parse_log)

    return master_rows, parse_logs


def main() -> int:
    parser = argparse.ArgumentParser(description="BIOPEP-UWM experimental ACE raw -> standardized 标准化脚本（纯标准库版）")
    parser.add_argument(
        "--project-root",
        type=str,
        default=None,
        help="PepDB 项目根目录，例如 E:\\MYS\\PepDB",
    )
    parser.add_argument(
        "--input-dir",
        type=str,
        default=None,
        help="BIOPEP raw_tables_csv 目录；默认自动推断到 DB/raw/ace/databases/BIOPEP_UWM/experimental/raw_tables_csv",
    )
    parser.add_argument(
        "--output-path",
        type=str,
        default=None,
        help="标准化主表输出路径；默认输出到 DB/standardized/ace/biopep_uwm/biopep_uwm_master_clean.csv",
    )
    args = parser.parse_args()

    script_file = Path(__file__)
    project_root = Path(args.project_root) if args.project_root else infer_project_root(script_file)

    input_dir = Path(args.input_dir) if args.input_dir else (
        project_root / "DB" / "raw" / "ace" / "databases" / "BIOPEP_UWM" / "experimental" / "raw_tables_csv"
    )
    output_path = Path(args.output_path) if args.output_path else (
        project_root / "DB" / "standardized" / "ace" / "biopep_uwm" / "biopep_uwm_master_clean.csv"
    )

    standardized_dir = output_path.parent
    analysis_dir = project_root / "DB" / "analysis" / "ace" / "biopep_uwm" / "标准化检查"

    parse_log_path = standardized_dir / "biopep_uwm_parse_log.csv"
    summary_path = standardized_dir / "biopep_uwm_standardize_summary.csv"
    missing_stats_path = analysis_dir / "字段缺失统计.csv"
    parse_status_stats_path = analysis_dir / "解析状态统计.csv"
    duplicate_id_path = analysis_dir / "重复peptide_id检查.csv"
    duplicate_seq_path = analysis_dir / "重复sequence检查.csv"
    length_consistency_path = analysis_dir / "长度一致性统计.csv"

    list_path = input_dir / LIST_FILE
    detail_path = input_dir / DETAIL_FILE

    list_rows = read_csv_rows(list_path, LIST_REQUIRED_COLUMNS)
    detail_rows = read_csv_rows(detail_path, DETAIL_REQUIRED_COLUMNS)

    master_rows, parse_logs = standardize(list_rows, detail_rows)

    write_csv(output_path, master_rows, MASTER_COLUMNS)
    write_csv(parse_log_path, parse_logs, list(parse_logs[0].keys()) if parse_logs else [])

    missing_stats = count_missing(master_rows, MASTER_COLUMNS)
    write_csv(missing_stats_path, missing_stats, ["field_name", "total_rows", "missing_count", "non_missing_count", "missing_ratio"])

    parse_counter = Counter(str(x.get("ic50_parse_status") or "") for x in master_rows)
    parse_status_rows = build_counter_table(parse_counter, "ic50_parse_status", "count")
    write_csv(parse_status_stats_path, parse_status_rows, ["ic50_parse_status", "count"])

    duplicate_id_rows = build_duplicate_id_table(master_rows)
    if not duplicate_id_rows:
        duplicate_id_rows = [{"source_record_id": "", "duplicate_count": 0}]
    write_csv(duplicate_id_path, duplicate_id_rows, ["source_record_id", "duplicate_count"])

    duplicate_seq_rows = build_sequence_duplicate_table(master_rows)
    if not duplicate_seq_rows:
        duplicate_seq_rows = [{"sequence": "", "duplicate_count": 0}]
    write_csv(duplicate_seq_path, duplicate_seq_rows, ["sequence", "duplicate_count"])

    length_counter = Counter(str(x.get("length_consistency_flag") or "") for x in master_rows)
    length_rows = build_counter_table(length_counter, "length_consistency_flag", "count")
    write_csv(length_consistency_path, length_rows, ["length_consistency_flag", "count"])

    summary_rows = build_summary_rows(list_rows, detail_rows, master_rows)
    write_csv(summary_path, summary_rows, ["metric", "value"])

    print("=" * 80)
    print("BIOPEP-UWM ACE 标准化完成")
    print(f"项目根目录：{project_root}")
    print(f"输入目录：{input_dir}")
    print(f"主表输出：{output_path}")
    print(f"解析日志：{parse_log_path}")
    print(f"检查目录：{analysis_dir}")
    print(f"list_rows={len(list_rows)} detail_rows={len(detail_rows)} master_rows={len(master_rows)}")
    print("=" * 80)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())