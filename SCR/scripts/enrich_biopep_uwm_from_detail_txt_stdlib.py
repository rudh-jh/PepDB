#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BIOPEP-UWM detail txt 二次富化脚本（纯标准库版）
================================================

用途
----
在 biopep_uwm_master_clean.csv 的基础上：
1. 读取每条记录对应的 raw_txt/detail_cards/detail_report_*.txt
2. 从富文本中回填：
   - ic50_value / ic50_unit / ic50_relation / ic50_uM / ic50_parse_status / ic50_parse_note
   - bibliographic_authors_raw / bibliographic_title_raw / bibliographic_year_raw
   - bibliographic_source_type_raw / bibliographic_raw
   - additional_information_raw
   - database_reference_raw
3. 输出 enriched 主表、log、以及若干检查表

设计原则
--------
- 不覆盖第一阶段 clean 主表
- 只做保守解析
- 原始字段尽量保留
- 回填优先级：detail txt > 第一阶段弱解析结果
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Tuple

INPUT_MASTER = "biopep_uwm_master_clean.csv"
OUTPUT_MASTER = "biopep_uwm_master_enriched.csv"
OUTPUT_LOG = "biopep_uwm_text_enrich_log.csv"

MOLAR_UNIT_TO_UM_FACTOR = {
    "nm": 1e-3,
    "um": 1.0,
    "mm": 1e3,
    "m": 1e6,
}

WEAK_PARSE_STATUSES = {
    "",
    "missing",
    "numeric_no_unit",
    "ambiguous_text",
    "unparsed_text",
    "percent_like",
    "range",
    "mass_unit",
}

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


# =========================
# 通用工具
# =========================
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
    text = re.sub(r"\s+", " ", text)
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
        return float(text)
    except ValueError:
        return None


def write_csv(path: Path, rows: List[Dict[str, object]], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


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


def normalize_unit_token(unit: str) -> str:
    unit = unit.strip()
    unit = unit.replace("μ", "u").replace("µ", "u")
    unit = unit.replace("micromole", "uM")
    unit = unit.replace("micromolar", "uM")
    unit = unit.replace("Micromole", "uM")
    unit = unit.replace("Micromolar", "uM")
    unit = unit.replace("UM", "uM").replace("NM", "nM").replace("MM", "mM")
    return unit


def unit_to_uM_factor(unit: str) -> Optional[float]:
    key = unit.strip().lower()
    return MOLAR_UNIT_TO_UM_FACTOR.get(key)


def canonical_unit_display(unit: str) -> str:
    unit = normalize_unit_token(unit)
    if unit.lower() == "um":
        return "uM"
    if unit.lower() == "nm":
        return "nM"
    if unit.lower() == "mm":
        return "mM"
    if unit == "M":
        return "M"
    return unit


def nonempty_count(*values: object) -> int:
    return sum(0 if is_missing(v) else 1 for v in values)


def first_nonmissing(*values: object) -> Optional[str]:
    for v in values:
        t = clean_text(v)
        if t is not None:
            return t
    return None


def safe_read_text(path: Path) -> Optional[str]:
    if not path.exists():
        return None
    for enc in ("utf-8-sig", "utf-8", "cp1252", "latin-1"):
        try:
            return path.read_text(encoding=enc, errors="replace")
        except Exception:
            continue
    return None


# =========================
# 路径解析
# =========================
def resolve_detail_txt_path(row: Dict[str, str], project_root: Path) -> Tuple[Optional[Path], str]:
    raw_text_path = clean_text(row.get("raw_text_path"))
    source_record_id = clean_text(row.get("source_record_id"))

    if raw_text_path:
        p = Path(raw_text_path)
        if p.exists():
            return p, "raw_text_path_exists"

        # 尝试把绝对路径退化成相对 DB/... 路径
        raw_lower = raw_text_path.lower().replace("/", "\\")
        marker = "\\db\\raw\\ace\\databases\\biopep_uwm\\experimental\\raw_txt\\detail_cards\\"
        idx = raw_lower.find(marker)
        if idx != -1:
            relative_tail = raw_text_path[idx + 1:].replace("\\", "/")
            candidate = project_root / relative_tail
            if candidate.exists():
                return candidate, "raw_text_path_rebased_from_project_root"

    if source_record_id:
        candidate = (
            project_root
            / "DB"
            / "raw"
            / "ace"
            / "databases"
            / "BIOPEP_UWM"
            / "experimental"
            / "raw_txt"
            / "detail_cards"
            / f"detail_report_{source_record_id}.txt"
        )
        if candidate.exists():
            return candidate, "fallback_from_source_record_id"

    return None, "detail_txt_not_found"


# =========================
# 文本区块提取
# =========================
def extract_section(text: str, start_label: str, end_labels: List[str]) -> Optional[str]:
    idx = text.find(start_label)
    if idx == -1:
        return None
    start = idx + len(start_label)
    end = len(text)
    for lab in end_labels:
        pos = text.find(lab, start)
        if pos != -1 and pos < end:
            end = pos
    value = normalize_spaces(text[start:end])
    return value or None


def extract_bibliography_raw(text: str) -> Optional[str]:
    return extract_section(
        text,
        "Bibliographic data:",
        ["Additional information:", "SMILES:", "Database reference:"]
    )


def extract_additional_information(text: str) -> Optional[str]:
    return extract_section(
        text,
        "Additional information:",
        ["SMILES:", "InChI=", "InChIKey:", "Information concerning", "Database reference:"]
    )


def extract_database_reference(text: str) -> Optional[str]:
    return extract_section(
        text,
        "Database reference:",
        []
    )


def extract_authors_from_biblio(biblio: str) -> Optional[str]:
    m = re.search(r"Authors\s+(.+?)\s+Title\s+", biblio, flags=re.IGNORECASE)
    if m:
        return normalize_spaces(m.group(1))
    return None


def extract_year_from_biblio(biblio: str) -> Optional[str]:
    # 优先 Year Source 之后的 4 位数
    m = re.search(r"Year\s+Source\s+(\d{4})", biblio, flags=re.IGNORECASE)
    if m:
        return m.group(1)

    # 退路：取最后一个 1900-2099 年份
    years = re.findall(r"\b(19\d{2}|20\d{2})\b", biblio)
    if years:
        return years[-1]
    return None


def extract_source_type_from_biblio(biblio: str) -> Optional[str]:
    m = re.search(r"Year\s+Source\s+\d{4}\s+([A-Za-z][A-Za-z \-/]+)$", biblio, flags=re.IGNORECASE)
    if m:
        return normalize_spaces(m.group(1))

    low = biblio.lower()
    if "journal" in low:
        return "Journal"
    if "book" in low:
        return "Book"
    if "patent" in low:
        return "Patent"
    if "database" in low:
        return "Database"
    return None


def extract_title_from_biblio(biblio: str) -> Optional[str]:
    m = re.search(r"Title\s+(.+?)\s+Year\s+Source\s+\d{4}", biblio, flags=re.IGNORECASE)
    if not m:
        return None

    title_part = normalize_spaces(m.group(1))

    # 尽量去掉末尾夹带的期刊卷期页码
    title_part = re.sub(r"\s+[A-Z][A-Za-z\.\-]+.*?\(\d{4}\)\s*$", "", title_part)
    title_part = normalize_spaces(title_part)

    return title_part or None


# =========================
# IC50 解析
# =========================
def parse_ic50_from_text(detail_text: str) -> Dict[str, object]:
    result = {
        "ic50_raw_from_text": None,
        "ic50_value": None,
        "ic50_unit": None,
        "ic50_relation": None,
        "ic50_uM": None,
        "ic50_parse_status": "missing",
        "ic50_parse_note": None,
    }

    patterns = [
        # IC50 : 29.00 µM
        r"\bIC50\s*:?\s*([<>≤≥]?)\s*([0-9]+(?:\.[0-9]+)?)\s*(µM|μM|uM|nM|mM|M|micromole|micromolar)\b",
        # The IC50 ... was 29.00 µM
        r"\bIC50\b.*?\b(?:was|is|of)\s*([<>≤≥]?)\s*([0-9]+(?:\.[0-9]+)?)\s*(µM|μM|uM|nM|mM|M|micromole|micromolar)\b",
    ]

    for pat in patterns:
        m = re.search(pat, detail_text, flags=re.IGNORECASE)
        if m:
            relation = m.group(1) if m.group(1) else "="
            value = float(m.group(2))
            unit_raw = m.group(3)
            unit = canonical_unit_display(unit_raw)
            factor = unit_to_uM_factor(unit)
            if factor is None:
                continue

            result["ic50_raw_from_text"] = normalize_spaces(m.group(0))
            result["ic50_value"] = value
            result["ic50_unit"] = unit
            result["ic50_relation"] = relation
            result["ic50_uM"] = value * factor
            result["ic50_parse_status"] = "exact_molar" if relation == "=" else "molar_threshold"
            result["ic50_parse_note"] = "parsed_from_detail_txt"
            return result

    # 如果只有数字，没有单位，不强行升级
    m_num = re.search(r"\bIC50\s*:?\s*([<>≤≥]?)\s*([0-9]+(?:\.[0-9]+)?)\b", detail_text, flags=re.IGNORECASE)
    if m_num:
        result["ic50_raw_from_text"] = normalize_spaces(m_num.group(0))
        result["ic50_value"] = float(m_num.group(2))
        result["ic50_relation"] = m_num.group(1) if m_num.group(1) else "="
        result["ic50_parse_status"] = "numeric_no_unit"
        result["ic50_parse_note"] = "detail_txt_numeric_but_unit_missing"
        return result

    result["ic50_parse_status"] = "missing"
    result["ic50_parse_note"] = "detail_txt_ic50_not_found"
    return result


def should_upgrade_ic50(row: Dict[str, str], parsed: Dict[str, object]) -> bool:
    current_status = clean_text(row.get("ic50_parse_status")) or ""
    current_uM = to_float_safe(row.get("ic50_uM"))
    current_unit = clean_text(row.get("ic50_unit"))

    new_status = clean_text(parsed.get("ic50_parse_status")) or ""
    new_uM = parsed.get("ic50_uM")
    new_unit = clean_text(parsed.get("ic50_unit"))

    if new_status in {"exact_molar", "molar_threshold"}:
        if current_status in WEAK_PARSE_STATUSES:
            return True
        if current_uM is None:
            return True
        if current_unit is None and new_unit is not None:
            return True

    return False


# =========================
# 备注更新
# =========================
def clean_notes_tokens(notes: Optional[str]) -> List[str]:
    text = clean_text(notes)
    if text is None:
        return []
    return [x.strip() for x in text.split(";") if x.strip()]


def update_notes(old_notes: Optional[str], row_after: Dict[str, object], detail_found: bool) -> Optional[str]:
    tokens = clean_notes_tokens(old_notes)

    # 若已成功回填 bibliography，去掉旧的 bibliography_missing
    if clean_text(row_after.get("bibliographic_raw")):
        tokens = [x for x in tokens if x != "bibliography_missing"]

    # 若 IC50 已升为 exact_molar / molar_threshold，去掉旧的 ic50_needs_manual_review
    if clean_text(row_after.get("ic50_parse_status")) in {"exact_molar", "molar_threshold"}:
        tokens = [x for x in tokens if x != "ic50_needs_manual_review"]

    if not detail_found and "detail_txt_missing" not in tokens:
        tokens.append("detail_txt_missing")

    return "; ".join(tokens) if tokens else None


# =========================
# 行处理
# =========================
def enrich_one_row(row: Dict[str, str], project_root: Path) -> Tuple[Dict[str, object], Dict[str, object]]:
    out = dict(row)

    path, path_status = resolve_detail_txt_path(row, project_root)
    detail_found = path is not None
    detail_text_raw = safe_read_text(path) if path else None
    detail_text = normalize_spaces(detail_text_raw) if detail_text_raw else None

    before_ic50_status = clean_text(row.get("ic50_parse_status")) or ""
    before_biblio_count = nonempty_count(
        row.get("bibliographic_authors_raw"),
        row.get("bibliographic_title_raw"),
        row.get("bibliographic_year_raw"),
        row.get("bibliographic_source_type_raw"),
        row.get("bibliographic_raw"),
    )

    ic50_changed = False
    bibliography_changed = False
    additional_info_changed = False
    db_ref_changed = False

    ic50_raw_from_text = None
    detail_txt_note = None

    if detail_text:
        parsed_ic50 = parse_ic50_from_text(detail_text)
        ic50_raw_from_text = parsed_ic50["ic50_raw_from_text"]

        if should_upgrade_ic50(row, parsed_ic50):
            out["ic50_value"] = parsed_ic50["ic50_value"]
            out["ic50_unit"] = parsed_ic50["ic50_unit"]
            out["ic50_relation"] = parsed_ic50["ic50_relation"]
            out["ic50_uM"] = parsed_ic50["ic50_uM"]
            out["ic50_parse_status"] = parsed_ic50["ic50_parse_status"]
            out["ic50_parse_note"] = parsed_ic50["ic50_parse_note"]
            # 只在原始 ic50_raw 缺失或很弱时补一份
            if is_missing(out.get("ic50_raw")) and parsed_ic50["ic50_raw_from_text"]:
                out["ic50_raw"] = parsed_ic50["ic50_raw_from_text"]
            ic50_changed = True

        bibliographic_raw = extract_bibliography_raw(detail_text)
        additional_information_raw = extract_additional_information(detail_text)
        database_reference_raw = extract_database_reference(detail_text)

        if bibliographic_raw:
            parsed_authors = extract_authors_from_biblio(bibliographic_raw)
            parsed_title = extract_title_from_biblio(bibliographic_raw)
            parsed_year = extract_year_from_biblio(bibliographic_raw)
            parsed_source_type = extract_source_type_from_biblio(bibliographic_raw)

            before = (
                clean_text(out.get("bibliographic_authors_raw")),
                clean_text(out.get("bibliographic_title_raw")),
                clean_text(out.get("bibliographic_year_raw")),
                clean_text(out.get("bibliographic_source_type_raw")),
                clean_text(out.get("bibliographic_raw")),
            )

            out["bibliographic_authors_raw"] = first_nonmissing(parsed_authors, out.get("bibliographic_authors_raw"))
            out["bibliographic_title_raw"] = first_nonmissing(parsed_title, out.get("bibliographic_title_raw"))
            out["bibliographic_year_raw"] = first_nonmissing(parsed_year, out.get("bibliographic_year_raw"))
            out["bibliographic_source_type_raw"] = first_nonmissing(parsed_source_type, out.get("bibliographic_source_type_raw"))
            out["bibliographic_raw"] = first_nonmissing(bibliographic_raw, out.get("bibliographic_raw"))

            after = (
                clean_text(out.get("bibliographic_authors_raw")),
                clean_text(out.get("bibliographic_title_raw")),
                clean_text(out.get("bibliographic_year_raw")),
                clean_text(out.get("bibliographic_source_type_raw")),
                clean_text(out.get("bibliographic_raw")),
            )
            bibliography_changed = before != after

        if additional_information_raw and is_missing(out.get("additional_information_raw")):
            out["additional_information_raw"] = additional_information_raw
            additional_info_changed = True

        if database_reference_raw and is_missing(out.get("database_reference_raw")):
            out["database_reference_raw"] = database_reference_raw
            db_ref_changed = True

        detail_txt_note = "detail_txt_read_ok"
    else:
        detail_txt_note = "detail_txt_unavailable"

    out["notes"] = update_notes(row.get("notes"), out, detail_found)

    after_biblio_count = nonempty_count(
        out.get("bibliographic_authors_raw"),
        out.get("bibliographic_title_raw"),
        out.get("bibliographic_year_raw"),
        out.get("bibliographic_source_type_raw"),
        out.get("bibliographic_raw"),
    )

    log_row = {
        "record_id": out.get("record_id"),
        "source_record_id": out.get("source_record_id"),
        "raw_text_path_in_master": row.get("raw_text_path"),
        "resolved_detail_txt_path": str(path) if path else "",
        "detail_txt_path_status": path_status,
        "detail_txt_found": 1 if detail_found else 0,
        "detail_txt_note": detail_txt_note,
        "ic50_before_status": before_ic50_status,
        "ic50_after_status": out.get("ic50_parse_status"),
        "ic50_raw_from_text": ic50_raw_from_text,
        "ic50_changed": 1 if ic50_changed else 0,
        "bibliography_before_nonmissing_count": before_biblio_count,
        "bibliography_after_nonmissing_count": after_biblio_count,
        "bibliography_changed": 1 if bibliography_changed else 0,
        "additional_information_changed": 1 if additional_info_changed else 0,
        "database_reference_changed": 1 if db_ref_changed else 0,
        "notes_after": out.get("notes"),
    }

    return out, log_row


# =========================
# 统计输出
# =========================
def build_metric_rows(pairs: List[Tuple[str, object]]) -> List[Dict[str, object]]:
    return [{"metric": k, "value": v} for k, v in pairs]


def build_counter_rows(counter: Counter, key_name: str, value_name: str) -> List[Dict[str, object]]:
    return [{key_name: k, value_name: v} for k, v in counter.most_common()]


def build_unresolved_ic50_rows(rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
    out = []
    for row in rows:
        status = clean_text(row.get("ic50_parse_status")) or ""
        if status not in {"exact_molar", "molar_threshold"}:
            out.append({
                "record_id": row.get("record_id"),
                "source_record_id": row.get("source_record_id"),
                "sequence": row.get("sequence"),
                "ic50_raw": row.get("ic50_raw"),
                "ic50_value": row.get("ic50_value"),
                "ic50_unit": row.get("ic50_unit"),
                "ic50_relation": row.get("ic50_relation"),
                "ic50_uM": row.get("ic50_uM"),
                "ic50_parse_status": row.get("ic50_parse_status"),
                "ic50_parse_note": row.get("ic50_parse_note"),
                "raw_text_path": row.get("raw_text_path"),
                "detail_report_url": row.get("detail_report_url"),
                "notes": row.get("notes"),
            })
    return out


def build_length_conflict_rows(rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
    out = []
    for row in rows:
        if clean_text(row.get("length_consistency_flag")) == "conflict":
            out.append({
                "record_id": row.get("record_id"),
                "source_record_id": row.get("source_record_id"),
                "sequence": row.get("sequence"),
                "peptide_length": row.get("peptide_length"),
                "number_of_residues_raw": row.get("number_of_residues_raw"),
                "length_consistency_flag": row.get("length_consistency_flag"),
                "raw_text_path": row.get("raw_text_path"),
                "detail_report_url": row.get("detail_report_url"),
                "notes": row.get("notes"),
            })
    return out


def summarize_ic50_backfill(before_rows: List[Dict[str, str]], after_rows: List[Dict[str, object]], logs: List[Dict[str, object]]) -> List[Dict[str, object]]:
    before_status = Counter(clean_text(r.get("ic50_parse_status")) or "" for r in before_rows)
    after_status = Counter(clean_text(r.get("ic50_parse_status")) or "" for r in after_rows)

    before_resolved = sum(1 for r in before_rows if (clean_text(r.get("ic50_parse_status")) or "") in {"exact_molar", "molar_threshold"})
    after_resolved = sum(1 for r in after_rows if (clean_text(r.get("ic50_parse_status")) or "") in {"exact_molar", "molar_threshold"})
    changed = sum(int(log.get("ic50_changed") or 0) for log in logs)
    detail_found = sum(int(log.get("detail_txt_found") or 0) for log in logs)

    rows = build_metric_rows([
        ("total_rows", len(after_rows)),
        ("detail_txt_found_rows", detail_found),
        ("ic50_changed_rows", changed),
        ("resolved_ic50_rows_before", before_resolved),
        ("resolved_ic50_rows_after", after_resolved),
        ("resolved_ic50_rows_gain", after_resolved - before_resolved),
    ])

    for k, v in before_status.items():
        rows.append({"metric": f"before_status::{k}", "value": v})
    for k, v in after_status.items():
        rows.append({"metric": f"after_status::{k}", "value": v})

    return rows


def summarize_bibliography_backfill(before_rows: List[Dict[str, str]], after_rows: List[Dict[str, object]], logs: List[Dict[str, object]]) -> List[Dict[str, object]]:
    fields = [
        "bibliographic_authors_raw",
        "bibliographic_title_raw",
        "bibliographic_year_raw",
        "bibliographic_source_type_raw",
        "bibliographic_raw",
        "additional_information_raw",
        "database_reference_raw",
    ]

    out = []
    for field in fields:
        before_nonmissing = sum(0 if is_missing(r.get(field)) else 1 for r in before_rows)
        after_nonmissing = sum(0 if is_missing(r.get(field)) else 1 for r in after_rows)
        out.append({
            "field_name": field,
            "before_nonmissing_count": before_nonmissing,
            "after_nonmissing_count": after_nonmissing,
            "gain": after_nonmissing - before_nonmissing,
        })

    out.append({
        "field_name": "__bibliography_changed_rows__",
        "before_nonmissing_count": "",
        "after_nonmissing_count": sum(int(log.get("bibliography_changed") or 0) for log in logs),
        "gain": "",
    })
    return out


# =========================
# 主程序
# =========================
def main() -> int:
    parser = argparse.ArgumentParser(description="BIOPEP-UWM 从 detail txt 二次富化主表（纯标准库版）")
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
        help="输入主表路径；默认 DB/standardized/ace/biopep_uwm/biopep_uwm_master_clean.csv",
    )
    parser.add_argument(
        "--output-path",
        type=str,
        default=None,
        help="输出 enriched 主表路径；默认 DB/standardized/ace/biopep_uwm/biopep_uwm_master_enriched.csv",
    )
    args = parser.parse_args()

    script_file = Path(__file__)
    project_root = Path(args.project_root) if args.project_root else infer_project_root(script_file)

    input_path = Path(args.input_path) if args.input_path else (
        project_root / "DB" / "standardized" / "ace" / "biopep_uwm" / INPUT_MASTER
    )
    output_path = Path(args.output_path) if args.output_path else (
        project_root / "DB" / "standardized" / "ace" / "biopep_uwm" / OUTPUT_MASTER
    )

    standardized_dir = output_path.parent
    analysis_dir = project_root / "DB" / "analysis" / "ace" / "biopep_uwm" / "标准化检查"

    log_path = standardized_dir / OUTPUT_LOG
    ic50_stats_path = analysis_dir / "IC50文本回填统计.csv"
    biblio_stats_path = analysis_dir / "文献信息回填统计.csv"
    unresolved_ic50_path = analysis_dir / "仍未解析IC50记录.csv"
    length_conflict_path = analysis_dir / "长度冲突明细.csv"

    before_rows = read_csv_rows(input_path)
    if not before_rows:
        raise ValueError(f"输入主表为空：{input_path}")

    master_columns = list(before_rows[0].keys())

    after_rows: List[Dict[str, object]] = []
    logs: List[Dict[str, object]] = []

    for row in before_rows:
        enriched_row, log_row = enrich_one_row(row, project_root)
        after_rows.append(enriched_row)
        logs.append(log_row)

    write_csv(output_path, after_rows, master_columns)
    write_csv(log_path, logs, list(logs[0].keys()) if logs else [])

    ic50_stat_rows = summarize_ic50_backfill(before_rows, after_rows, logs)
    write_csv(ic50_stats_path, ic50_stat_rows, ["metric", "value"])

    biblio_stat_rows = summarize_bibliography_backfill(before_rows, after_rows, logs)
    write_csv(
        biblio_stats_path,
        biblio_stat_rows,
        ["field_name", "before_nonmissing_count", "after_nonmissing_count", "gain"],
    )

    unresolved_rows = build_unresolved_ic50_rows(after_rows)
    if not unresolved_rows:
        unresolved_rows = [{
            "record_id": "",
            "source_record_id": "",
            "sequence": "",
            "ic50_raw": "",
            "ic50_value": "",
            "ic50_unit": "",
            "ic50_relation": "",
            "ic50_uM": "",
            "ic50_parse_status": "",
            "ic50_parse_note": "",
            "raw_text_path": "",
            "detail_report_url": "",
            "notes": "",
        }]
    write_csv(
        unresolved_ic50_path,
        unresolved_rows,
        [
            "record_id", "source_record_id", "sequence", "ic50_raw", "ic50_value",
            "ic50_unit", "ic50_relation", "ic50_uM", "ic50_parse_status",
            "ic50_parse_note", "raw_text_path", "detail_report_url", "notes",
        ],
    )

    conflict_rows = build_length_conflict_rows(after_rows)
    if not conflict_rows:
        conflict_rows = [{
            "record_id": "",
            "source_record_id": "",
            "sequence": "",
            "peptide_length": "",
            "number_of_residues_raw": "",
            "length_consistency_flag": "",
            "raw_text_path": "",
            "detail_report_url": "",
            "notes": "",
        }]
    write_csv(
        length_conflict_path,
        conflict_rows,
        [
            "record_id", "source_record_id", "sequence", "peptide_length",
            "number_of_residues_raw", "length_consistency_flag",
            "raw_text_path", "detail_report_url", "notes",
        ],
    )

    detail_found_count = sum(int(x.get("detail_txt_found") or 0) for x in logs)
    ic50_changed_count = sum(int(x.get("ic50_changed") or 0) for x in logs)
    biblio_changed_count = sum(int(x.get("bibliography_changed") or 0) for x in logs)

    print("=" * 80)
    print("BIOPEP-UWM detail txt 二次富化完成")
    print(f"项目根目录：{project_root}")
    print(f"输入主表：{input_path}")
    print(f"输出主表：{output_path}")
    print(f"日志文件：{log_path}")
    print(f"分析目录：{analysis_dir}")
    print(f"总记录数：{len(after_rows)}")
    print(f"找到 detail txt 的记录数：{detail_found_count}")
    print(f"IC50 被升级/回填的记录数：{ic50_changed_count}")
    print(f"bibliography 被回填的记录数：{biblio_changed_count}")
    print("=" * 80)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())