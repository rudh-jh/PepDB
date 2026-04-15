#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FermFooDb ACE 标准化脚本（纯标准库版）
====================================

输入
----
1) DB/raw/ace/databases/FermFooDb/experimental/raw_tables_csv/fermfoodb_ace_list.csv
2) DB/raw/ace/databases/FermFooDb/experimental/raw_tables_csv/fermfoodb_ace_detail.csv

输出
----
DB/standardized/ace/fermfoodb/
├── fermfoodb_ace_master_clean.csv
├── fermfoodb_ace_parse_log.csv
└── fermfoodb_ace_standardize_summary.csv

DB/analysis/ace/fermfoodb/标准化检查/
├── 字段缺失统计.csv
├── 解析状态统计.csv
├── 重复sequence检查.csv
├── activity分布统计.csv
└── food_matrix分布统计.csv

说明
----
- 主表以 list.csv 为主，detail.csv 作为按 fmdb_id 的补充信息源
- 重点对 ic50_raw 做解析和单位统一
- 对于摩尔浓度单位，统一转换为 uM
- 对于质量浓度单位（mg/mL, ug/mL, ng/mL），若 mass 可解析为 Da，则换算为 uM
- 不可解析或不可换算的 IC50 保留 raw，并标记 parse_status
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

# 允许长文本字段
csv.field_size_limit(10**8)

INPUT_LIST = "fermfoodb_ace_list.csv"
INPUT_DETAIL = "fermfoodb_ace_detail.csv"

MASTER_OUTPUT = "fermfoodb_ace_master_clean.csv"
PARSE_LOG_OUTPUT = "fermfoodb_ace_parse_log.csv"
SUMMARY_OUTPUT = "fermfoodb_ace_standardize_summary.csv"

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
    "food_matrix_raw",
    "protein_raw",
    "culture_raw",
    "hydrolysis_raw",
    "experiment_raw",
    "model_raw",
    "assay_raw",
    "method_of_analysis_raw",
    "mz_ratio_raw",
    "mass_raw",
    "mass_da",
    "pubmed_id_raw",
    "title_raw",
    "ph_raw",
    "temperature_raw",
    "incubation_time_raw",
    "detail_url",
    "notes",
]

PARSE_LOG_COLUMNS = [
    "record_id",
    "source_record_id",
    "sequence",
    "peptide_length",
    "ic50_raw",
    "mass_raw",
    "mass_da",
    "ic50_value",
    "ic50_unit",
    "ic50_relation",
    "ic50_uM",
    "ic50_parse_status",
    "ic50_parse_note",
    "ic50_numeric_count",
    "ic50_numeric_values_serialized",
    "activity_label_raw",
    "notes",
]


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
    # 去掉常见分隔符，但不强行删除字母外字符，后续 notes 会标记
    text = text.replace("-", "")
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


def prefer(*values: object) -> Optional[str]:
    for v in values:
        t = clean_text(v)
        if t is not None:
            return t
    return None


def write_csv(path: Path, rows: List[Dict[str, object]], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def read_csv_rows(path: Path) -> List[Dict[str, Optional[str]]]:
    if not path.exists():
        raise FileNotFoundError(f"找不到输入文件：{path}")

    with open(path, "r", encoding="utf-8-sig", errors="replace", newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError(f"CSV 无表头：{path}")
        rows: List[Dict[str, Optional[str]]] = []
        for row in reader:
            clean_row: Dict[str, Optional[str]] = {}
            for k, v in row.items():
                clean_row[strip_bom(k)] = clean_text(v)
            rows.append(clean_row)
    return rows


def resolve_input_path(input_path: Optional[str], default_path: Path) -> Path:
    return Path(input_path) if input_path else default_path


def build_detail_index(detail_rows: List[Dict[str, Optional[str]]]) -> Dict[str, Dict[str, Optional[str]]]:
    """
    detail 可能存在重复 FMDB_ID，优先保留“字段更丰富”的那条
    """
    detail_index: Dict[str, Dict[str, Optional[str]]] = {}
    detail_score: Dict[str, int] = {}

    for row in detail_rows:
        fmdb_id = clean_text(row.get("fmdb_id"))
        if not fmdb_id:
            continue

        score = sum(
            1
            for k, v in row.items()
            if k not in {"fmdb_id", "detail_url", "detail_text_path", "detail_html_path"} and clean_text(v) is not None
        )

        if fmdb_id not in detail_index or score > detail_score[fmdb_id]:
            detail_index[fmdb_id] = row
            detail_score[fmdb_id] = score

    return detail_index


def parse_mass_da(raw_mass: object, raw_mz_ratio: object) -> Tuple[Optional[float], Optional[str]]:
    """
    先用 mass，取首个数值作为 Da；
    如果 mass 缺失，再尝试 m/z ratio
    """
    mass_text = clean_text(raw_mass)
    if mass_text:
        nums = re.findall(r"[-+]?\d+(?:\.\d+)?", mass_text)
        if nums:
            try:
                return float(nums[0]), "mass_numeric_extracted"
            except ValueError:
                pass

    mz_text = clean_text(raw_mz_ratio)
    if mz_text:
        nums = re.findall(r"[-+]?\d+(?:\.\d+)?", mz_text)
        if nums:
            try:
                return float(nums[0]), "mz_numeric_extracted_as_fallback"
            except ValueError:
                pass

    return None, None


def normalize_unit_text(text: str) -> str:
    x = text.lower()
    x = x.replace("μ", "u").replace("µ", "u")
    x = x.replace("μmol", "umol").replace("µmol", "umol")
    x = x.replace("μm", "um").replace("µm", "um")
    x = x.replace("μg", "ug").replace("µg", "ug")
    x = x.replace("l-1", "/l")
    x = x.replace(" ", "")
    return x


def detect_relation(text: str) -> str:
    text = text.strip()
    if text.startswith("≤"):
        return "≤"
    if text.startswith("≥"):
        return "≥"
    if text.startswith("<"):
        return "<"
    if text.startswith(">"):
        return ">"
    return "="


def parse_ic50(raw_value: object, mass_da: Optional[float]) -> Dict[str, object]:
    """
    解析并归一化 ic50_raw

    支持示例：
    - 2.15±0.02 uM
    - 720uM
    - 300.1umol/l
    - 173.3 umol/l
    - >1000umol/l
    - 1.2 mg/ml  （若 mass_da 可用）
    - 250 ug/ml  （若 mass_da 可用）
    """
    result = {
        "ic50_value": None,
        "ic50_unit": None,
        "ic50_relation": None,
        "ic50_uM": None,
        "ic50_parse_status": "missing",
        "ic50_parse_note": None,
        "ic50_numeric_count": 0,
        "ic50_numeric_values_serialized": None,
    }

    text = clean_text(raw_value)
    if text is None:
        return result

    relation = detect_relation(text)
    work = normalize_unit_text(text)

    # 去掉前导关系符
    work = re.sub(r"^[<>≤≥]+\s*", "", work)

    # 把 ± 误差项抹掉：2.15±0.02um -> 2.15um
    work = re.sub(r"(\d+(?:\.\d+)?)\s*±\s*\d+(?:\.\d+)?", r"\1", work)

    # 去掉括号内说明，但保留简单数值
    work = re.sub(r"\(([^)]*)\)", " ", work)
    work = normalize_spaces(work)

    numeric_strings = re.findall(r"[-+]?\d+(?:\.\d+)?", work)
    values: List[float] = []
    for x in numeric_strings:
        try:
            values.append(float(x))
        except ValueError:
            continue

    if not values:
        result["ic50_parse_status"] = "unparsed_text"
        result["ic50_parse_note"] = "text_present_but_no_numeric_value_found"
        result["ic50_relation"] = relation
        return result

    # 单元格中若有多个数，取中位数
    median_value = statistics.median(values)
    result["ic50_numeric_count"] = len(values)
    result["ic50_numeric_values_serialized"] = " | ".join(str(v) for v in values)
    result["ic50_relation"] = relation

    # 单位检测：从高优先级到低优先级
    unit_key = None

    # 质量浓度
    if "mg/ml" in work or "mgperml" in work:
        unit_key = "mg/ml"
    elif "ug/ml" in work or "ugperml" in work:
        unit_key = "ug/ml"
    elif "ng/ml" in work or "ngperml" in work:
        unit_key = "ng/ml"

    # 摩尔浓度
    elif "mmol/l" in work or re.search(r"(?<![a-z])mm\b", work):
        unit_key = "mM"
    elif "umol/l" in work or re.search(r"(?<![a-z])um\b", work):
        unit_key = "uM"
    elif "nmol/l" in work or re.search(r"(?<![a-z])nm\b", work):
        unit_key = "nM"
    elif "pmol/l" in work or re.search(r"(?<![a-z])pm\b", work):
        unit_key = "pM"
    elif "mol/l" in work:
        unit_key = "M"

    # 没写单位但有数值
    if unit_key is None:
        result["ic50_value"] = median_value
        result["ic50_unit"] = None
        result["ic50_uM"] = None
        result["ic50_parse_status"] = "numeric_no_unit"
        result["ic50_parse_note"] = f"numeric_value_present_but_unit_missing_n={len(values)}"
        return result

    # 摩尔浓度统一换算到 uM
    if unit_key == "uM":
        ic50_uM = median_value
        parse_note = "parsed_direct_uM_like_unit"
    elif unit_key == "mM":
        ic50_uM = median_value * 1000.0
        parse_note = "converted_from_mM_to_uM"
    elif unit_key == "nM":
        ic50_uM = median_value / 1000.0
        parse_note = "converted_from_nM_to_uM"
    elif unit_key == "pM":
        ic50_uM = median_value / 1_000_000.0
        parse_note = "converted_from_pM_to_uM"
    elif unit_key == "M":
        ic50_uM = median_value * 1_000_000.0
        parse_note = "converted_from_M_to_uM"

    # 质量浓度：若有分子量则换算
    elif unit_key == "mg/ml":
        if mass_da and mass_da > 0:
            # 1 mg/mL = 1 g/L；uM = (g/L / Da) * 1e6
            ic50_uM = median_value * (1_000_000.0 / mass_da)
            parse_note = "converted_from_mg_per_ml_to_uM_using_mass_da"
        else:
            result["ic50_value"] = median_value
            result["ic50_unit"] = unit_key
            result["ic50_parse_status"] = "mass_concentration_needs_mass"
            result["ic50_parse_note"] = "mg_per_ml_detected_but_mass_da_missing"
            return result
    elif unit_key == "ug/ml":
        if mass_da and mass_da > 0:
            # 1 ug/mL = 0.001 g/L；uM = (0.001 / Da) * 1e6 = 1000 / Da
            ic50_uM = median_value * (1000.0 / mass_da)
            parse_note = "converted_from_ug_per_ml_to_uM_using_mass_da"
        else:
            result["ic50_value"] = median_value
            result["ic50_unit"] = unit_key
            result["ic50_parse_status"] = "mass_concentration_needs_mass"
            result["ic50_parse_note"] = "ug_per_ml_detected_but_mass_da_missing"
            return result
    elif unit_key == "ng/ml":
        if mass_da and mass_da > 0:
            # 1 ng/mL = 1e-6 g/L；uM = (1e-6 / Da) * 1e6 = 1 / Da
            ic50_uM = median_value * (1.0 / mass_da)
            parse_note = "converted_from_ng_per_ml_to_uM_using_mass_da"
        else:
            result["ic50_value"] = median_value
            result["ic50_unit"] = unit_key
            result["ic50_parse_status"] = "mass_concentration_needs_mass"
            result["ic50_parse_note"] = "ng_per_ml_detected_but_mass_da_missing"
            return result
    else:
        result["ic50_value"] = median_value
        result["ic50_unit"] = unit_key
        result["ic50_parse_status"] = "unparsed_text"
        result["ic50_parse_note"] = "unit_detected_but_conversion_rule_missing"
        return result

    result["ic50_value"] = median_value
    result["ic50_unit"] = unit_key
    result["ic50_uM"] = ic50_uM
    result["ic50_parse_status"] = "exact_molar" if relation == "=" else "molar_threshold"
    if len(values) == 1:
        result["ic50_parse_note"] = parse_note
    else:
        result["ic50_parse_note"] = f"{parse_note}; multi_value_median_used_n={len(values)}"
    return result


def merge_list_and_detail(
    list_row: Dict[str, Optional[str]],
    detail_row: Optional[Dict[str, Optional[str]]],
) -> Dict[str, Optional[str]]:
    detail_row = detail_row or {}

    return {
        "fmdb_id": prefer(list_row.get("fmdb_id"), detail_row.get("fmdb_id")),
        "detail_url": prefer(list_row.get("detail_url"), detail_row.get("detail_url")),
        "pubmed_id_raw": prefer(detail_row.get("pubmed_id"), list_row.get("pubmed_id")),
        "title_raw": prefer(detail_row.get("detail_title"), list_row.get("title")),
        "sequence_raw": prefer(detail_row.get("peptide_sequence"), list_row.get("peptide_sequence")),
        "length_raw": prefer(detail_row.get("length_of_peptide"), list_row.get("length_of_peptide")),
        "food_matrix_raw": prefer(detail_row.get("food_matrix"), list_row.get("food_matrix")),
        "protein_raw": prefer(detail_row.get("protein_name"), list_row.get("protein")),
        "ph_raw": prefer(detail_row.get("ph"), list_row.get("ph")),
        "temperature_raw": prefer(detail_row.get("temperature"), list_row.get("temperature")),
        "incubation_time_raw": prefer(detail_row.get("incubation_time"), list_row.get("incubation_time")),
        "activity_label_raw": prefer(detail_row.get("activity"), list_row.get("activity")),
        "experiment_raw": prefer(detail_row.get("experiment"), list_row.get("experiment")),
        "model_raw": prefer(detail_row.get("model"), list_row.get("model")),
        "assay_raw": prefer(detail_row.get("assay_for_activity_measurement"), list_row.get("assay_for_activity_measurement")),
        "culture_raw": prefer(detail_row.get("starter_culture"), list_row.get("culture")),
        "hydrolysis_raw": prefer(detail_row.get("hydrolysis"), list_row.get("hydrolysis")),
        "method_of_analysis_raw": prefer(detail_row.get("method_of_analysis"), list_row.get("method_of_analysis")),
        "mz_ratio_raw": prefer(detail_row.get("mz_ratio"), list_row.get("mz_ratio")),
        "mass_raw": prefer(detail_row.get("mass"), list_row.get("mass")),
        "ic50_raw": prefer(detail_row.get("ic50_raw"), list_row.get("ic50_raw")),
    }


def build_notes(
    merged_raw: Dict[str, Optional[str]],
    sequence: Optional[str],
    peptide_length: Optional[int],
    ic50_info: Dict[str, object],
    detail_hit: bool,
    mass_note: Optional[str],
) -> Optional[str]:
    notes: List[str] = []

    if not detail_hit:
        notes.append("detail_missing_or_not_joined")

    activity = clean_text(merged_raw.get("activity_label_raw"))
    if activity and activity.lower() not in {"ace-inhibitory", "ace inhibitory", "aceinhibitory"}:
        notes.append("activity_label_not_exact_ace_inhibitory")

    if sequence is None:
        notes.append("sequence_missing")
    elif not sequence_is_canonical(sequence):
        notes.append("non_canonical_or_contains_unknown_residue")

    if peptide_length is None:
        notes.append("peptide_length_missing")
    elif sequence is not None and peptide_length != len(sequence):
        notes.append("peptide_length_conflicts_with_sequence_length")

    parse_status = str(ic50_info.get("ic50_parse_status") or "")
    if parse_status == "missing":
        notes.append("ic50_missing")
    elif parse_status in {"unparsed_text", "numeric_no_unit", "mass_concentration_needs_mass"}:
        notes.append("ic50_needs_manual_review")

    if mass_note:
        notes.append(mass_note)

    if clean_text(merged_raw.get("pubmed_id_raw")) is None:
        notes.append("pubmed_id_missing")

    return "; ".join(notes) if notes else None


def build_master_row(
    list_row: Dict[str, Optional[str]],
    detail_row: Optional[Dict[str, Optional[str]]],
    idx: int,
) -> Tuple[Dict[str, object], Dict[str, object]]:
    merged_raw = merge_list_and_detail(list_row, detail_row)
    detail_hit = detail_row is not None

    sequence = clean_sequence(merged_raw.get("sequence_raw"))

    peptide_length = to_int_safe(merged_raw.get("length_raw"))
    if peptide_length is None and sequence is not None:
        peptide_length = len(sequence)

    mass_da, mass_note = parse_mass_da(merged_raw.get("mass_raw"), merged_raw.get("mz_ratio_raw"))
    ic50_info = parse_ic50(merged_raw.get("ic50_raw"), mass_da=mass_da)

    notes = build_notes(
        merged_raw=merged_raw,
        sequence=sequence,
        peptide_length=peptide_length,
        ic50_info=ic50_info,
        detail_hit=detail_hit,
        mass_note=mass_note,
    )

    master_row = {
        "record_id": f"fermfoodb_ace_{idx:06d}",
        "task": "ace",
        "source_type": "database",
        "source_name": "FermFooDb",
        "source_record_id": clean_text(merged_raw.get("fmdb_id")),
        "target": "ACE",
        "sequence": sequence,
        "peptide_length": peptide_length,
        "activity_label_raw": clean_text(merged_raw.get("activity_label_raw")),
        "ic50_raw": clean_text(merged_raw.get("ic50_raw")),
        "ic50_value": ic50_info["ic50_value"],
        "ic50_unit": ic50_info["ic50_unit"],
        "ic50_relation": ic50_info["ic50_relation"],
        "ic50_uM": ic50_info["ic50_uM"],
        "ic50_parse_status": ic50_info["ic50_parse_status"],
        "ic50_parse_note": ic50_info["ic50_parse_note"],
        "food_matrix_raw": clean_text(merged_raw.get("food_matrix_raw")),
        "protein_raw": clean_text(merged_raw.get("protein_raw")),
        "culture_raw": clean_text(merged_raw.get("culture_raw")),
        "hydrolysis_raw": clean_text(merged_raw.get("hydrolysis_raw")),
        "experiment_raw": clean_text(merged_raw.get("experiment_raw")),
        "model_raw": clean_text(merged_raw.get("model_raw")),
        "assay_raw": clean_text(merged_raw.get("assay_raw")),
        "method_of_analysis_raw": clean_text(merged_raw.get("method_of_analysis_raw")),
        "mz_ratio_raw": clean_text(merged_raw.get("mz_ratio_raw")),
        "mass_raw": clean_text(merged_raw.get("mass_raw")),
        "mass_da": mass_da,
        "pubmed_id_raw": clean_text(merged_raw.get("pubmed_id_raw")),
        "title_raw": clean_text(merged_raw.get("title_raw")),
        "ph_raw": clean_text(merged_raw.get("ph_raw")),
        "temperature_raw": clean_text(merged_raw.get("temperature_raw")),
        "incubation_time_raw": clean_text(merged_raw.get("incubation_time_raw")),
        "detail_url": clean_text(merged_raw.get("detail_url")),
        "notes": notes,
    }

    parse_log = {
        "record_id": master_row["record_id"],
        "source_record_id": master_row["source_record_id"],
        "sequence": master_row["sequence"],
        "peptide_length": master_row["peptide_length"],
        "ic50_raw": master_row["ic50_raw"],
        "mass_raw": master_row["mass_raw"],
        "mass_da": master_row["mass_da"],
        "ic50_value": master_row["ic50_value"],
        "ic50_unit": master_row["ic50_unit"],
        "ic50_relation": master_row["ic50_relation"],
        "ic50_uM": master_row["ic50_uM"],
        "ic50_parse_status": master_row["ic50_parse_status"],
        "ic50_parse_note": master_row["ic50_parse_note"],
        "ic50_numeric_count": ic50_info["ic50_numeric_count"],
        "ic50_numeric_values_serialized": ic50_info["ic50_numeric_values_serialized"],
        "activity_label_raw": master_row["activity_label_raw"],
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


def build_food_matrix_rows(master_rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
    counter = Counter(str(x.get("food_matrix_raw") or "") for x in master_rows if x.get("food_matrix_raw"))
    return build_counter_rows(counter, "food_matrix_raw", "count")


def build_activity_rows(master_rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
    counter = Counter(str(x.get("activity_label_raw") or "") for x in master_rows if x.get("activity_label_raw"))
    return build_counter_rows(counter, "activity_label_raw", "count")


def build_summary_rows(
    list_rows: List[Dict[str, Optional[str]]],
    detail_rows: List[Dict[str, Optional[str]]],
    detail_index: Dict[str, Dict[str, Optional[str]]],
    master_rows: List[Dict[str, object]],
    input_list_path: Path,
    input_detail_path: Path,
) -> List[Dict[str, object]]:
    summary = [
        {"metric": "input_list_filename", "value": input_list_path.name},
        {"metric": "input_detail_filename", "value": input_detail_path.name},
        {"metric": "input_list_rows", "value": len(list_rows)},
        {"metric": "input_detail_rows", "value": len(detail_rows)},
        {"metric": "unique_detail_fmdb_id", "value": len(detail_index)},
        {"metric": "master_rows", "value": len(master_rows)},
    ]

    detail_joined_rows = sum(
        1
        for row in master_rows
        if row.get("detail_url")
    )
    summary.append({"metric": "rows_with_detail_url", "value": detail_joined_rows})

    parse_counter = Counter(str(x.get("ic50_parse_status") or "") for x in master_rows)
    for k, v in parse_counter.items():
        summary.append({"metric": f"ic50_parse_status::{k}", "value": v})

    unit_counter = Counter(str(x.get("ic50_unit") or "") for x in master_rows if x.get("ic50_unit"))
    for k, v in unit_counter.items():
        summary.append({"metric": f"ic50_unit::{k}", "value": v})

    detail_hit_counter = Counter(
        "detail_present" if ("detail_missing_or_not_joined" not in str(x.get("notes") or "")) else "detail_missing"
        for x in master_rows
    )
    for k, v in detail_hit_counter.items():
        summary.append({"metric": f"detail_join::{k}", "value": v})

    noncanonical_rows = sum(
        1 for x in master_rows
        if "non_canonical_or_contains_unknown_residue" in str(x.get("notes") or "")
    )
    summary.append({"metric": "rows_with_noncanonical_sequence", "value": noncanonical_rows})

    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description="FermFooDb ACE 标准化脚本（纯标准库版）")
    parser.add_argument(
        "--project-root",
        type=str,
        default=None,
        help="PepDB 项目根目录，例如 E:\\MYS\\PepDB",
    )
    parser.add_argument(
        "--list-input",
        type=str,
        default=None,
        help="list CSV 路径；默认 raw_tables_csv/fermfoodb_ace_list.csv",
    )
    parser.add_argument(
        "--detail-input",
        type=str,
        default=None,
        help="detail CSV 路径；默认 raw_tables_csv/fermfoodb_ace_detail.csv",
    )
    args = parser.parse_args()

    script_file = Path(__file__)
    project_root = Path(args.project_root) if args.project_root else infer_project_root(script_file)

    list_input = resolve_input_path(
        args.list_input,
        project_root / "DB" / "raw" / "ace" / "databases" / "FermFooDb" / "experimental" / "raw_tables_csv" / INPUT_LIST,
    )
    detail_input = resolve_input_path(
        args.detail_input,
        project_root / "DB" / "raw" / "ace" / "databases" / "FermFooDb" / "experimental" / "raw_tables_csv" / INPUT_DETAIL,
    )

    standardized_dir = project_root / "DB" / "standardized" / "ace" / "fermfoodb"
    analysis_dir = project_root / "DB" / "analysis" / "ace" / "fermfoodb" / "标准化检查"

    master_output = standardized_dir / MASTER_OUTPUT
    parse_log_output = standardized_dir / PARSE_LOG_OUTPUT
    summary_output = standardized_dir / SUMMARY_OUTPUT

    missing_stats_output = analysis_dir / "字段缺失统计.csv"
    parse_stats_output = analysis_dir / "解析状态统计.csv"
    duplicate_output = analysis_dir / "重复sequence检查.csv"
    activity_output = analysis_dir / "activity分布统计.csv"
    food_matrix_output = analysis_dir / "food_matrix分布统计.csv"

    list_rows = read_csv_rows(list_input)
    detail_rows = read_csv_rows(detail_input)
    detail_index = build_detail_index(detail_rows)

    master_rows: List[Dict[str, object]] = []
    parse_logs: List[Dict[str, object]] = []

    for idx, list_row in enumerate(list_rows, start=1):
        fmdb_id = clean_text(list_row.get("fmdb_id"))
        detail_row = detail_index.get(fmdb_id) if fmdb_id else None
        master_row, parse_log = build_master_row(list_row, detail_row, idx)
        master_rows.append(master_row)
        parse_logs.append(parse_log)

    write_csv(master_output, master_rows, MASTER_COLUMNS)
    write_csv(parse_log_output, parse_logs, PARSE_LOG_COLUMNS)

    missing_stats = count_missing(master_rows, MASTER_COLUMNS)
    write_csv(
        missing_stats_output,
        missing_stats,
        ["field_name", "total_rows", "missing_count", "non_missing_count", "missing_ratio"],
    )

    parse_counter = Counter(str(x.get("ic50_parse_status") or "") for x in master_rows)
    parse_rows = build_counter_rows(parse_counter, "ic50_parse_status", "count")
    write_csv(parse_stats_output, parse_rows, ["ic50_parse_status", "count"])

    duplicate_rows = build_duplicate_sequence_rows(master_rows)
    if not duplicate_rows:
        duplicate_rows = [{"sequence": "", "duplicate_count": 0}]
    write_csv(duplicate_output, duplicate_rows, ["sequence", "duplicate_count"])

    activity_rows = build_activity_rows(master_rows)
    if not activity_rows:
        activity_rows = [{"activity_label_raw": "", "count": 0}]
    write_csv(activity_output, activity_rows, ["activity_label_raw", "count"])

    food_matrix_rows = build_food_matrix_rows(master_rows)
    if not food_matrix_rows:
        food_matrix_rows = [{"food_matrix_raw": "", "count": 0}]
    write_csv(food_matrix_output, food_matrix_rows, ["food_matrix_raw", "count"])

    summary_rows = build_summary_rows(
        list_rows=list_rows,
        detail_rows=detail_rows,
        detail_index=detail_index,
        master_rows=master_rows,
        input_list_path=list_input,
        input_detail_path=detail_input,
    )
    write_csv(summary_output, summary_rows, ["metric", "value"])

    print("=" * 80)
    print("FermFooDb ACE 标准化完成")
    print(f"项目根目录：{project_root}")
    print(f"list 输入：{list_input}")
    print(f"detail 输入：{detail_input}")
    print(f"master 输出：{master_output}")
    print(f"parse log：{parse_log_output}")
    print(f"summary：{summary_output}")
    print(f"input_list_rows={len(list_rows)} input_detail_rows={len(detail_rows)} master_rows={len(master_rows)}")
    print("=" * 80)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())