#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FermFooDb ACE IC50 rescue / enrich 脚本（纯标准库版）
===================================================

输入
----
1) DB/standardized/ace/fermfoodb/fermfoodb_ace_master_clean.csv
2) DB/raw/ace/databases/FermFooDb/experimental/raw_tables_csv/fermfoodb_ace_list.csv
3) DB/raw/ace/databases/FermFooDb/experimental/raw_tables_csv/fermfoodb_ace_detail.csv

输出
----
DB/standardized/ace/fermfoodb/
├── fermfoodb_ace_master_rescued.csv
├── fermfoodb_ace_ic50_rescue_log.csv
└── fermfoodb_ace_ic50_rescue_summary.csv

DB/analysis/ace/fermfoodb/IC50补救检查/
├── 补救前后解析状态统计.csv
├── 补救动作统计.csv
└── 字段回填统计.csv

说明
----
- 不覆盖原始 standardized 主表
- 以 source_record_id / fmdb_id 为主键，把 raw list/detail 重新补回 standardized
- 优先挽救：
  1) numeric_no_unit
  2) mass_concentration_needs_mass
  3) 部分写法怪异的 unparsed_text
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

csv.field_size_limit(10**8)

MASTER_INPUT = "fermfoodb_ace_master_clean.csv"
RAW_LIST_INPUT = "fermfoodb_ace_list.csv"
RAW_DETAIL_INPUT = "fermfoodb_ace_detail.csv"

MASTER_OUTPUT = "fermfoodb_ace_master_rescued.csv"
LOG_OUTPUT = "fermfoodb_ace_ic50_rescue_log.csv"
SUMMARY_OUTPUT = "fermfoodb_ace_ic50_rescue_summary.csv"

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

# 单位 Da = g/mol；这里用单体残基质量 + H2O 估算肽分子量
AA_MONOISOTOPIC_RESIDUE_MASS = {
    "A": 71.03711,
    "C": 103.00919,
    "D": 115.02694,
    "E": 129.04259,
    "F": 147.06841,
    "G": 57.02146,
    "H": 137.05891,
    "I": 113.08406,
    "K": 128.09496,
    "L": 113.08406,
    "M": 131.04049,
    "N": 114.04293,
    "P": 97.05276,
    "Q": 128.05858,
    "R": 156.10111,
    "S": 87.03203,
    "T": 101.04768,
    "V": 99.06841,
    "W": 186.07931,
    "Y": 163.06333,
}
WATER_MASS = 18.01056

BASE_COLUMNS = [
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

EXTRA_COLUMNS = [
    "mass_da_source",
    "rescue_changed_flag",
    "rescue_action",
]

MASTER_OUTPUT_COLUMNS = BASE_COLUMNS + EXTRA_COLUMNS

LOG_COLUMNS = [
    "record_id",
    "source_record_id",
    "sequence",
    "ic50_raw_before",
    "ic50_parse_status_before",
    "ic50_uM_before",
    "mass_da_before",
    "ic50_raw_after",
    "ic50_parse_status_after",
    "ic50_uM_after",
    "mass_da_after",
    "mass_da_source_after",
    "rescue_action",
    "pubmed_backfilled",
    "title_backfilled",
    "protein_backfilled",
    "food_matrix_backfilled",
    "culture_backfilled",
    "hydrolysis_backfilled",
    "experiment_backfilled",
    "model_backfilled",
    "assay_backfilled",
    "method_backfilled",
    "mz_ratio_backfilled",
    "mass_backfilled",
    "notes_after",
]


def infer_project_root(script_file: Path) -> Path:
    try:
        return script_file.resolve().parents[2]
    except IndexError:
        return script_file.resolve().parent


def normalize_spaces(text: str) -> str:
    text = text.replace("\ufeff", " ").replace("\r", " ").replace("\n", " ")
    return re.sub(r"\s+", " ", text).strip()


def strip_bom(text: str) -> str:
    return text.lstrip("\ufeff")


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


def read_csv_rows(path: Path) -> List[Dict[str, Optional[str]]]:
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


def resolve_path(user_path: Optional[str], default_path: Path) -> Path:
    return Path(user_path) if user_path else default_path


def normalize_unit_text(text: str) -> str:
    x = text.lower()
    x = x.replace("μ", "u").replace("µ", "u")
    x = x.replace("μmol", "umol").replace("µmol", "umol")
    x = x.replace("μm", "um").replace("µm", "um")
    x = x.replace("μg", "ug").replace("µg", "ug")
    x = x.replace("mcg", "ug")
    x = x.replace("l-1", "/l")
    x = x.replace("perml", "/ml")
    x = x.replace("perl", "/l")
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


def estimate_mass_da_from_sequence(seq: Optional[str]) -> Optional[float]:
    if not sequence_is_canonical(seq):
        return None
    total = WATER_MASS
    for aa in seq:
        total += AA_MONOISOTOPIC_RESIDUE_MASS[aa]
    return total


def parse_mass_da(raw_mass: object, raw_mz_ratio: object) -> Tuple[Optional[float], Optional[str]]:
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


def build_raw_index(rows: List[Dict[str, Optional[str]]], key_field: str) -> Dict[str, Dict[str, Optional[str]]]:
    """
    同键可能多条，保留“信息最丰富”的那条。
    """
    out: Dict[str, Dict[str, Optional[str]]] = {}
    score_map: Dict[str, int] = {}

    for row in rows:
        key = clean_text(row.get(key_field))
        if key is None:
            continue
        score = sum(1 for _, v in row.items() if clean_text(v) is not None)
        if key not in out or score > score_map[key]:
            out[key] = row
            score_map[key] = score

    return out


def prefer_nonempty(*values: object) -> Optional[str]:
    for v in values:
        t = clean_text(v)
        if t is not None:
            return t
    return None


def append_note(notes: Optional[str], token: str) -> str:
    existing = [x.strip() for x in str(notes or "").split(";") if x.strip()]
    if token not in existing:
        existing.append(token)
    return "; ".join(existing)


def parse_ic50_advanced(raw_value: object, mass_da: Optional[float]) -> Dict[str, object]:
    """
    更保守但更宽容的 IC50 解析器
    """
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

    relation = detect_relation(text)
    work = normalize_unit_text(text)

    # 去前导关系符
    work = re.sub(r"^[<>≤≥]+\s*", "", work)

    # 把 2.15±0.02um -> 2.15um
    work = re.sub(r"(\d+(?:\.\d+)?)\s*±\s*\d+(?:\.\d+)?", r"\1", work)

    # 兼容 like 2.15+/-0.02uM
    work = re.sub(r"(\d+(?:\.\d+)?)\s*\+/-\s*\d+(?:\.\d+)?", r"\1", work)

    # 去括号
    work = re.sub(r"\(([^)]*)\)", " ", work)
    work = normalize_spaces(work)

    nums = re.findall(r"[-+]?\d+(?:\.\d+)?", work)
    if not nums:
        result["ic50_relation"] = relation
        result["ic50_parse_status"] = "unparsed_text"
        result["ic50_parse_note"] = "no_numeric_token_found"
        return result

    values = []
    for x in nums:
        try:
            values.append(float(x))
        except ValueError:
            continue

    if not values:
        result["ic50_relation"] = relation
        result["ic50_parse_status"] = "unparsed_text"
        result["ic50_parse_note"] = "numeric_conversion_failed"
        return result

    numeric_value = statistics.median(values)
    result["ic50_value"] = numeric_value
    result["ic50_relation"] = relation

    unit_key = None

    # 质量浓度
    if "mg/ml" in work:
        unit_key = "mg/ml"
    elif "ug/ml" in work:
        unit_key = "ug/ml"
    elif "ng/ml" in work:
        unit_key = "ng/ml"

    # 摩尔浓度
    elif "mmol/l" in work or re.search(r"(?<![a-z])mm(?![a-z])", work):
        unit_key = "mM"
    elif "umol/l" in work or "micromolar" in work or re.search(r"(?<![a-z])um(?![a-z])", work):
        unit_key = "uM"
    elif "nmol/l" in work or re.search(r"(?<![a-z])nm(?![a-z])", work):
        unit_key = "nM"
    elif "pmol/l" in work or re.search(r"(?<![a-z])pm(?![a-z])", work):
        unit_key = "pM"
    elif "mol/l" in work:
        unit_key = "M"

    # 没单位
    if unit_key is None:
        result["ic50_parse_status"] = "numeric_no_unit"
        result["ic50_parse_note"] = "numeric_value_present_but_unit_missing"
        return result

    # 换算到 uM
    if unit_key == "uM":
        ic50_uM = numeric_value
        note = "parsed_uM_like_unit"
    elif unit_key == "mM":
        ic50_uM = numeric_value * 1000.0
        note = "converted_from_mM_to_uM"
    elif unit_key == "nM":
        ic50_uM = numeric_value / 1000.0
        note = "converted_from_nM_to_uM"
    elif unit_key == "pM":
        ic50_uM = numeric_value / 1_000_000.0
        note = "converted_from_pM_to_uM"
    elif unit_key == "M":
        ic50_uM = numeric_value * 1_000_000.0
        note = "converted_from_M_to_uM"
    elif unit_key == "mg/ml":
        if mass_da and mass_da > 0:
            ic50_uM = numeric_value * (1_000_000.0 / mass_da)
            note = "converted_from_mg_per_ml_to_uM_using_mass_da"
        else:
            result["ic50_unit"] = unit_key
            result["ic50_parse_status"] = "mass_concentration_needs_mass"
            result["ic50_parse_note"] = "mg_per_ml_but_mass_da_missing"
            return result
    elif unit_key == "ug/ml":
        if mass_da and mass_da > 0:
            ic50_uM = numeric_value * (1000.0 / mass_da)
            note = "converted_from_ug_per_ml_to_uM_using_mass_da"
        else:
            result["ic50_unit"] = unit_key
            result["ic50_parse_status"] = "mass_concentration_needs_mass"
            result["ic50_parse_note"] = "ug_per_ml_but_mass_da_missing"
            return result
    elif unit_key == "ng/ml":
        if mass_da and mass_da > 0:
            ic50_uM = numeric_value * (1.0 / mass_da)
            note = "converted_from_ng_per_ml_to_uM_using_mass_da"
        else:
            result["ic50_unit"] = unit_key
            result["ic50_parse_status"] = "mass_concentration_needs_mass"
            result["ic50_parse_note"] = "ng_per_ml_but_mass_da_missing"
            return result
    else:
        result["ic50_unit"] = unit_key
        result["ic50_parse_status"] = "unparsed_text"
        result["ic50_parse_note"] = "unit_detected_but_conversion_rule_missing"
        return result

    result["ic50_unit"] = unit_key
    result["ic50_uM"] = ic50_uM
    result["ic50_parse_status"] = "exact_molar" if relation == "=" else "molar_threshold"
    result["ic50_parse_note"] = note
    return result


def build_sequence_context(rows: List[Dict[str, str]]) -> Dict[str, Dict[str, object]]:
    """
    为 numeric_no_unit 准备 sequence 级上下文：
    只收集已经解析成功的 exact/molar_threshold 记录
    """
    bucket: Dict[str, List[float]] = defaultdict(list)
    for row in rows:
        seq = clean_sequence(row.get("sequence"))
        status = clean_text(row.get("ic50_parse_status"))
        ic50_uM = to_float_safe(row.get("ic50_uM"))
        if seq and status in {"exact_molar", "molar_threshold"} and ic50_uM is not None:
            bucket[seq].append(ic50_uM)

    out = {}
    for seq, vals in bucket.items():
        out[seq] = {
            "count": len(vals),
            "median_uM": statistics.median(vals),
            "min_uM": min(vals),
            "max_uM": max(vals),
        }
    return out


def maybe_rescue_numeric_no_unit(
    row: Dict[str, object],
    seq_context: Dict[str, Dict[str, object]],
) -> Tuple[Optional[float], Optional[str], Optional[str]]:
    """
    仅对 numeric_no_unit 做保守补救：
    如果同 sequence 已有至少 2 条 resolved 记录，就把裸数字当作 uM
    """
    seq = clean_sequence(row.get("sequence"))
    raw = clean_text(row.get("ic50_raw"))
    if seq is None or raw is None:
        return None, None, None

    nums = re.findall(r"[-+]?\d+(?:\.\d+)?", raw)
    if not nums:
        return None, None, None

    try:
        numeric_value = statistics.median([float(x) for x in nums])
    except Exception:
        return None, None, None

    ctx = seq_context.get(seq)
    if not ctx:
        return None, None, None

    if int(ctx["count"]) < 2:
        return None, None, None

    note = (
        f"rescued_numeric_no_unit_assumed_uM_from_sequence_context"
        f"; context_count={ctx['count']}"
        f"; context_median_uM={ctx['median_uM']}"
    )
    return numeric_value, "uM", note


def count_backfill(before: object, after: object) -> int:
    return int(clean_text(before) is None and clean_text(after) is not None)


def main() -> int:
    parser = argparse.ArgumentParser(description="FermFooDb ACE IC50 rescue / enrich 脚本（纯标准库版）")
    parser.add_argument("--project-root", type=str, default=None, help="PepDB 项目根目录")
    parser.add_argument("--master-input", type=str, default=None, help="标准化主表路径")
    parser.add_argument("--raw-list-input", type=str, default=None, help="raw list 表路径")
    parser.add_argument("--raw-detail-input", type=str, default=None, help="raw detail 表路径")
    args = parser.parse_args()

    script_file = Path(__file__)
    project_root = Path(args.project_root) if args.project_root else infer_project_root(script_file)

    master_input = resolve_path(
        args.master_input,
        project_root / "DB" / "standardized" / "ace" / "fermfoodb" / MASTER_INPUT,
    )
    raw_list_input = resolve_path(
        args.raw_list_input,
        project_root / "DB" / "raw" / "ace" / "databases" / "FermFooDb" / "experimental" / "raw_tables_csv" / RAW_LIST_INPUT,
    )
    raw_detail_input = resolve_path(
        args.raw_detail_input,
        project_root / "DB" / "raw" / "ace" / "databases" / "FermFooDb" / "experimental" / "raw_tables_csv" / RAW_DETAIL_INPUT,
    )

    standardized_dir = project_root / "DB" / "standardized" / "ace" / "fermfoodb"
    analysis_dir = project_root / "DB" / "analysis" / "ace" / "fermfoodb" / "IC50补救检查"

    master_output = standardized_dir / MASTER_OUTPUT
    log_output = standardized_dir / LOG_OUTPUT
    summary_output = standardized_dir / SUMMARY_OUTPUT

    parse_compare_output = analysis_dir / "补救前后解析状态统计.csv"
    action_stats_output = analysis_dir / "补救动作统计.csv"
    backfill_stats_output = analysis_dir / "字段回填统计.csv"

    master_rows = read_csv_rows(master_input)
    raw_list_rows = read_csv_rows(raw_list_input)
    raw_detail_rows = read_csv_rows(raw_detail_input)

    raw_list_index = build_raw_index(raw_list_rows, "fmdb_id")
    raw_detail_index = build_raw_index(raw_detail_rows, "fmdb_id")

    seq_context = build_sequence_context(master_rows)

    rescued_rows: List[Dict[str, object]] = []
    rescue_logs: List[Dict[str, object]] = []

    before_counter = Counter()
    after_counter = Counter()
    action_counter = Counter()
    backfill_counter = Counter()

    for row in master_rows:
        out = dict(row)
        source_record_id = clean_text(out.get("source_record_id"))
        raw_list = raw_list_index.get(source_record_id or "")
        raw_detail = raw_detail_index.get(source_record_id or "")

        # 记录补救前
        status_before = clean_text(out.get("ic50_parse_status")) or ""
        before_counter[status_before] += 1

        mass_da_before = to_float_safe(out.get("mass_da"))
        ic50_uM_before = to_float_safe(out.get("ic50_uM"))

        pubmed_before = out.get("pubmed_id_raw")
        title_before = out.get("title_raw")
        protein_before = out.get("protein_raw")
        food_before = out.get("food_matrix_raw")
        culture_before = out.get("culture_raw")
        hydrolysis_before = out.get("hydrolysis_raw")
        experiment_before = out.get("experiment_raw")
        model_before = out.get("model_raw")
        assay_before = out.get("assay_raw")
        method_before = out.get("method_of_analysis_raw")
        mz_before = out.get("mz_ratio_raw")
        mass_before = out.get("mass_raw")

        # ===== 1) 字段回填 =====
        out["pubmed_id_raw"] = prefer_nonempty(
            out.get("pubmed_id_raw"),
            (raw_detail or {}).get("pubmed_id"),
            (raw_list or {}).get("pubmed_id"),
        )
        out["title_raw"] = prefer_nonempty(
            out.get("title_raw"),
            (raw_detail or {}).get("detail_title"),
            (raw_list or {}).get("title"),
        )
        out["protein_raw"] = prefer_nonempty(
            out.get("protein_raw"),
            (raw_detail or {}).get("protein_name"),
            (raw_list or {}).get("protein"),
        )
        out["food_matrix_raw"] = prefer_nonempty(
            out.get("food_matrix_raw"),
            (raw_detail or {}).get("food_matrix"),
            (raw_list or {}).get("food_matrix"),
        )
        out["culture_raw"] = prefer_nonempty(
            out.get("culture_raw"),
            (raw_detail or {}).get("starter_culture"),
            (raw_list or {}).get("culture"),
        )
        out["hydrolysis_raw"] = prefer_nonempty(
            out.get("hydrolysis_raw"),
            (raw_detail or {}).get("hydrolysis"),
            (raw_list or {}).get("hydrolysis"),
        )
        out["experiment_raw"] = prefer_nonempty(
            out.get("experiment_raw"),
            (raw_detail or {}).get("experiment"),
            (raw_list or {}).get("experiment"),
        )
        out["model_raw"] = prefer_nonempty(
            out.get("model_raw"),
            (raw_detail or {}).get("model"),
            (raw_list or {}).get("model"),
        )
        out["assay_raw"] = prefer_nonempty(
            out.get("assay_raw"),
            (raw_detail or {}).get("assay_for_activity_measurement"),
            (raw_list or {}).get("assay_for_activity_measurement"),
        )
        out["method_of_analysis_raw"] = prefer_nonempty(
            out.get("method_of_analysis_raw"),
            (raw_detail or {}).get("method_of_analysis"),
            (raw_list or {}).get("method_of_analysis"),
        )
        out["mz_ratio_raw"] = prefer_nonempty(
            out.get("mz_ratio_raw"),
            (raw_detail or {}).get("mz_ratio"),
            (raw_list or {}).get("mz_ratio"),
        )
        out["mass_raw"] = prefer_nonempty(
            out.get("mass_raw"),
            (raw_detail or {}).get("mass"),
            (raw_list or {}).get("mass"),
        )
        out["detail_url"] = prefer_nonempty(
            out.get("detail_url"),
            (raw_detail or {}).get("detail_url"),
            (raw_list or {}).get("detail_url"),
        )
        out["activity_label_raw"] = prefer_nonempty(
            out.get("activity_label_raw"),
            (raw_detail or {}).get("activity"),
            (raw_list or {}).get("activity"),
        )
        out["ph_raw"] = prefer_nonempty(
            out.get("ph_raw"),
            (raw_detail or {}).get("ph"),
            (raw_list or {}).get("ph"),
        )
        out["temperature_raw"] = prefer_nonempty(
            out.get("temperature_raw"),
            (raw_detail or {}).get("temperature"),
            (raw_list or {}).get("temperature"),
        )
        out["incubation_time_raw"] = prefer_nonempty(
            out.get("incubation_time_raw"),
            (raw_detail or {}).get("incubation_time"),
            (raw_list or {}).get("incubation_time"),
        )
        out["ic50_raw"] = prefer_nonempty(
            out.get("ic50_raw"),
            (raw_detail or {}).get("ic50_raw"),
            (raw_list or {}).get("ic50_raw"),
        )

        # ===== 2) 质量补救 =====
        mass_da_source = clean_text(out.get("mass_da_source"))
        mass_da = to_float_safe(out.get("mass_da"))

        if mass_da is None:
            mass_da, mass_note = parse_mass_da(out.get("mass_raw"), out.get("mz_ratio_raw"))
            if mass_da is not None:
                mass_da_source = mass_note

        if mass_da is None:
            seq = clean_sequence(out.get("sequence"))
            est_mass = estimate_mass_da_from_sequence(seq)
            if est_mass is not None:
                mass_da = est_mass
                mass_da_source = "estimated_from_sequence"

        out["mass_da"] = mass_da
        out["mass_da_source"] = mass_da_source

        # ===== 3) IC50 重新解析 =====
        parsed = parse_ic50_advanced(out.get("ic50_raw"), mass_da=mass_da)

        rescue_action = []
        changed = 0

        status_old = clean_text(out.get("ic50_parse_status"))
        ic50_old = to_float_safe(out.get("ic50_uM"))

        # 直接采用 improved parser 的结果，如果更好
        old_resolved = status_old in {"exact_molar", "molar_threshold"} and ic50_old is not None
        new_resolved = parsed["ic50_parse_status"] in {"exact_molar", "molar_threshold"} and parsed["ic50_uM"] is not None

        if (not old_resolved and new_resolved) or (old_resolved and new_resolved and parsed["ic50_uM"] != ic50_old):
            out["ic50_value"] = parsed["ic50_value"]
            out["ic50_unit"] = parsed["ic50_unit"]
            out["ic50_relation"] = parsed["ic50_relation"]
            out["ic50_uM"] = parsed["ic50_uM"]
            out["ic50_parse_status"] = parsed["ic50_parse_status"]
            out["ic50_parse_note"] = parsed["ic50_parse_note"]
            changed = 1
            rescue_action.append("reparsed_ic50_with_improved_rules")

        # numeric_no_unit 的 sequence-context 补救
        status_now = clean_text(out.get("ic50_parse_status"))
        if status_now == "numeric_no_unit":
            rescued_value, rescued_unit, rescued_note = maybe_rescue_numeric_no_unit(out, seq_context)
            if rescued_value is not None:
                out["ic50_value"] = rescued_value
                out["ic50_unit"] = rescued_unit
                out["ic50_relation"] = out.get("ic50_relation") or "="
                out["ic50_uM"] = rescued_value
                out["ic50_parse_status"] = "rescued_numeric_no_unit"
                out["ic50_parse_note"] = rescued_note
                changed = 1
                rescue_action.append("rescued_numeric_no_unit_from_sequence_context")

        # 质量估算补救
        if mass_da_source == "estimated_from_sequence":
            rescue_action.append("mass_da_estimated_from_sequence")

        # notes 更新
        notes_after = clean_text(out.get("notes")) or ""
        if changed:
            notes_after = append_note(notes_after, "ic50_rescued_or_updated")
        if mass_da_source == "estimated_from_sequence":
            notes_after = append_note(notes_after, "mass_da_estimated_from_sequence")
        out["notes"] = notes_after

        out["rescue_changed_flag"] = 1 if changed else 0
        out["rescue_action"] = "; ".join(dict.fromkeys(rescue_action)) if rescue_action else ""

        status_after = clean_text(out.get("ic50_parse_status")) or ""
        after_counter[status_after] += 1

        # 统计
        action_counter[out["rescue_action"] or "no_action"] += 1
        backfill_counter["pubmed_id_raw"] += count_backfill(pubmed_before, out.get("pubmed_id_raw"))
        backfill_counter["title_raw"] += count_backfill(title_before, out.get("title_raw"))
        backfill_counter["protein_raw"] += count_backfill(protein_before, out.get("protein_raw"))
        backfill_counter["food_matrix_raw"] += count_backfill(food_before, out.get("food_matrix_raw"))
        backfill_counter["culture_raw"] += count_backfill(culture_before, out.get("culture_raw"))
        backfill_counter["hydrolysis_raw"] += count_backfill(hydrolysis_before, out.get("hydrolysis_raw"))
        backfill_counter["experiment_raw"] += count_backfill(experiment_before, out.get("experiment_raw"))
        backfill_counter["model_raw"] += count_backfill(model_before, out.get("model_raw"))
        backfill_counter["assay_raw"] += count_backfill(assay_before, out.get("assay_raw"))
        backfill_counter["method_of_analysis_raw"] += count_backfill(method_before, out.get("method_of_analysis_raw"))
        backfill_counter["mz_ratio_raw"] += count_backfill(mz_before, out.get("mz_ratio_raw"))
        backfill_counter["mass_raw"] += count_backfill(mass_before, out.get("mass_raw"))

        rescue_logs.append({
            "record_id": out.get("record_id"),
            "source_record_id": out.get("source_record_id"),
            "sequence": out.get("sequence"),
            "ic50_raw_before": row.get("ic50_raw"),
            "ic50_parse_status_before": status_before,
            "ic50_uM_before": ic50_uM_before,
            "mass_da_before": mass_da_before,
            "ic50_raw_after": out.get("ic50_raw"),
            "ic50_parse_status_after": out.get("ic50_parse_status"),
            "ic50_uM_after": out.get("ic50_uM"),
            "mass_da_after": out.get("mass_da"),
            "mass_da_source_after": out.get("mass_da_source"),
            "rescue_action": out.get("rescue_action"),
            "pubmed_backfilled": count_backfill(pubmed_before, out.get("pubmed_id_raw")),
            "title_backfilled": count_backfill(title_before, out.get("title_raw")),
            "protein_backfilled": count_backfill(protein_before, out.get("protein_raw")),
            "food_matrix_backfilled": count_backfill(food_before, out.get("food_matrix_raw")),
            "culture_backfilled": count_backfill(culture_before, out.get("culture_raw")),
            "hydrolysis_backfilled": count_backfill(hydrolysis_before, out.get("hydrolysis_raw")),
            "experiment_backfilled": count_backfill(experiment_before, out.get("experiment_raw")),
            "model_backfilled": count_backfill(model_before, out.get("model_raw")),
            "assay_backfilled": count_backfill(assay_before, out.get("assay_raw")),
            "method_backfilled": count_backfill(method_before, out.get("method_of_analysis_raw")),
            "mz_ratio_backfilled": count_backfill(mz_before, out.get("mz_ratio_raw")),
            "mass_backfilled": count_backfill(mass_before, out.get("mass_raw")),
            "notes_after": out.get("notes"),
        })

        rescued_rows.append(out)

    # 输出 rescued master / log
    write_csv(master_output, rescued_rows, MASTER_OUTPUT_COLUMNS)
    write_csv(log_output, rescue_logs, LOG_COLUMNS)

    # summary
    summary_rows: List[Dict[str, object]] = [
        {"metric": "input_master_rows", "value": len(master_rows)},
        {"metric": "raw_list_rows", "value": len(raw_list_rows)},
        {"metric": "raw_detail_rows", "value": len(raw_detail_rows)},
        {"metric": "rescued_master_rows", "value": len(rescued_rows)},
    ]

    for k, v in before_counter.items():
        summary_rows.append({"metric": f"before_parse_status::{k}", "value": v})
    for k, v in after_counter.items():
        summary_rows.append({"metric": f"after_parse_status::{k}", "value": v})
    for k, v in action_counter.items():
        summary_rows.append({"metric": f"rescue_action::{k}", "value": v})
    for k, v in backfill_counter.items():
        summary_rows.append({"metric": f"backfill::{k}", "value": v})

    write_csv(summary_output, summary_rows, ["metric", "value"])

    # analysis
    parse_compare_rows = []
    all_statuses = sorted(set(before_counter) | set(after_counter))
    for status in all_statuses:
        parse_compare_rows.append({
            "ic50_parse_status": status,
            "before_count": before_counter.get(status, 0),
            "after_count": after_counter.get(status, 0),
            "delta": after_counter.get(status, 0) - before_counter.get(status, 0),
        })
    write_csv(parse_compare_output, parse_compare_rows, ["ic50_parse_status", "before_count", "after_count", "delta"])

    action_rows = [{"rescue_action": k, "count": v} for k, v in action_counter.most_common()]
    write_csv(action_stats_output, action_rows, ["rescue_action", "count"])

    backfill_rows = [{"field_name": k, "backfilled_count": v} for k, v in backfill_counter.items()]
    write_csv(backfill_stats_output, backfill_rows, ["field_name", "backfilled_count"])

    print("=" * 80)
    print("FermFooDb ACE IC50 rescue / enrich 完成")
    print(f"项目根目录：{project_root}")
    print(f"输入主表：{master_input}")
    print(f"输入 raw list：{raw_list_input}")
    print(f"输入 raw detail：{raw_detail_input}")
    print(f"输出主表：{master_output}")
    print(f"输出日志：{log_output}")
    print(f"输出摘要：{summary_output}")
    print("=" * 80)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())