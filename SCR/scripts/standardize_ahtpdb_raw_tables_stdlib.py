#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AHTPDB raw_tables 三表清洗与标准化脚本（纯标准库版）
====================================================

为什么要写这个版本
------------------
你本机当前报错不是业务逻辑有问题，而是 NumPy / pandas / pyarrow 的二进制兼容问题。
为了让脚本在“只有 Python 标准库”的环境中也能直接运行，这个版本完全不依赖 pandas。

本脚本服务的目标（按当前 PepDB 阶段收束）
---------------------------------------
1. 读取 AHTPDB 的三张原始表：
   - ahtpdb_small_peptides_2_5_residues.txt
   - ahtpdb_long_peptides_6_plus_residues.txt
   - ahtpdb_peptides_with_ic50.txt
2. 以 small + long 作为主体记录表
3. 用 peptides_with_ic50 表补充 / 优先活性信息
4. 只重点清洗最关键字段：
   - seq（序列）
   - len（长度）
   - molwt（分子量）
   - ic50（活性）
5. 输出一张可直接用于后续分析的总表 CSV

刻意不做的事情
--------------
- 不做复杂数据库建模
- 不删除异常记录，只加标记列
- 不合并重复序列，只统计重复和冲突
- 不把可疑单位（如 uM/L）直接改写成 uM
- 不在第一版中强行把 mg/mL、ug/mL 全部换算成 uM
- 不对 source / mice / method / assay 做过度结构化

总表可直接支持的分析
--------------------
- 长度分布
- 二肽/三肽数量
- 精确 IC50 样本量
- 单位分布
- 重复序列与冲突值
- N/C 端位点富集
- 强活性 vs 全体富集
- 顺序翻转对比

默认输入 / 输出路径
--------------------
假定脚本放在：PepDB/SCR/scripts/
则默认：
- 输入目录：PepDB/DB/raw/ace/databases/AHTPDB/raw_tables/
- 输出文件：PepDB/DB/standardized/ace/ahtpdb/ahtpdb_master_clean.csv

运行示例
--------
1) 若脚本已放入 PepDB/SCR/scripts/：
   python standardize_ahtpdb_raw_tables_stdlib.py

2) 指定项目根目录：
   python standardize_ahtpdb_raw_tables_stdlib.py --project-root E:\\MYS\\PepDB

3) 指定输出文件：
   python standardize_ahtpdb_raw_tables_stdlib.py --output-path E:\\MYS\\PepDB\\DB\\standardized\\ace\\ahtpdb\\ahtpdb_master_clean.csv
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# =========================
# 1. 基础常量
# =========================
SMALL_FILE = "ahtpdb_small_peptides_2_5_residues.txt"
LONG_FILE = "ahtpdb_long_peptides_6_plus_residues.txt"
IC50_FILE = "ahtpdb_peptides_with_ic50.txt"

EXPECTED_COLUMNS = [
    "id",
    "seq",
    "len",
    "molwt",
    "ic50",
    "source",
    "mice",
    "method",
    "assay",
    "bitter",
    "pi",
    "bp",
]

MISSING_TOKENS = {
    "",
    "nd",
    "n.d.",
    "na",
    "n/a",
    "nan",
    "none",
    "null",
    "-",
    "--",
}

CANONICAL_AA = set("ACDEFGHIKLMNPQRSTVWY")

# 仅用于把“明确摩尔浓度点值 / 阈值”统一到 uM。
MOLAR_UNIT_TO_UM_FACTOR = {
    "nm": 1e-3,
    "um": 1.0,
    "mm": 1e3,
    "m": 1e6,
}


# =========================
# 2. 通用辅助函数
# =========================
def infer_project_root(script_file: Path) -> Path:
    """
    根据脚本位置推断 PepDB 项目根目录。

    约定：脚本放在 PepDB/SCR/scripts/
    那么 script_file.parents[2] 就是 PepDB 根目录。
    """
    try:
        return script_file.resolve().parents[2]
    except IndexError:
        return script_file.resolve().parent


def normalize_spaces(text: str) -> str:
    """去掉首尾空格，并把中间连续空白压缩成单个空格。"""
    return re.sub(r"\s+", " ", text.strip())


def is_missing(value: object) -> bool:
    """判断一个值是否应视为缺失。"""
    if value is None:
        return True
    text = str(value).strip().lower()
    return text in MISSING_TOKENS


def clean_text(value: object) -> Optional[str]:
    """
    文本轻清洗：
    - 缺失返回 None
    - 非缺失返回“去首尾空格 + 压缩多余空格”的字符串
    """
    if is_missing(value):
        return None
    return normalize_spaces(str(value))


def to_float_safe(value: object) -> Optional[float]:
    """尝试把一个值转成 float；失败返回 None。"""
    text = clean_text(value)
    if text is None:
        return None
    text = text.replace(",", "")
    try:
        return float(text)
    except ValueError:
        return None


# =========================
# 3. 原始表读取
# =========================
def _split_line_to_fields(line: str) -> List[str]:
    """
    将原始文本的一行切成字段。

    读取策略：
    1. 若包含 tab，则优先按 tab 切分。
    2. 否则按“两个及以上空白”切分。

    之所以不按单个空格切，是因为 method / assay 等字段内部可能含有单个空格。
    """
    raw = line.rstrip("\n\r")
    if "\t" in raw:
        return [cell.strip() for cell in raw.split("\t")]
    return [cell.strip() for cell in re.split(r"\s{2,}", raw.strip())]


def read_raw_table(path: Path, table_from: str) -> List[Dict[str, Optional[str]]]:
    """
    读取 AHTPDB 原始表，返回“行字典列表”。

    设计原则：
    - 不依赖 pandas
    - 尽量兼容 tab 分隔和多空格分隔
    - 所有值先按字符串保留，避免自动类型推断带来的误伤
    """
    if not path.exists():
        raise FileNotFoundError(f"找不到输入文件：{path}")

    with open(path, "r", encoding="utf-8", errors="replace") as f:
        lines = [ln for ln in f if ln.strip()]

    if not lines:
        raise ValueError(f"文件为空：{path}")

    header = _split_line_to_fields(lines[0])
    if header != EXPECTED_COLUMNS:
        raise ValueError(
            f"文件 {path.name} 的列名不是预期格式。\n"
            f"当前列名：{header}\n"
            f"预期列名：{EXPECTED_COLUMNS}"
        )

    rows: List[Dict[str, Optional[str]]] = []
    for line_no, line in enumerate(lines[1:], start=2):
        parts = _split_line_to_fields(line)

        # 如果因为某些异常格式导致字段数不对，优先尝试保底修正；修不了就报错。
        if len(parts) != len(EXPECTED_COLUMNS):
            raise ValueError(
                f"文件 {path.name} 第 {line_no} 行字段数异常：实际 {len(parts)}，预期 {len(EXPECTED_COLUMNS)}。\n"
                f"原始行内容：{line.strip()}"
            )

        row: Dict[str, Optional[str]] = {col: clean_text(val) for col, val in zip(EXPECTED_COLUMNS, parts)}
        row["table_from"] = table_from
        rows.append(row)

    return rows


# =========================
# 4. Sequence / Length / MolWt 处理
# =========================
def clean_sequence(seq: object) -> Optional[str]:
    """
    清洗氨基酸序列：
    - 缺失返回 None
    - 转大写
    - 去空格
    - 不改顺序
    - 不自动替换可疑字符（例如 J）
    """
    text = clean_text(seq)
    if text is None:
        return None
    text = text.upper().replace(" ", "")
    return text if text else None


def is_valid_canonical_sequence(seq: Optional[str]) -> bool:
    """判断是否为仅包含标准 20 种氨基酸的序列。"""
    if not seq:
        return False
    return set(seq).issubset(CANONICAL_AA)


def compute_length(seq: Optional[str]) -> Optional[int]:
    """根据清洗后的序列重算长度。"""
    if not seq:
        return None
    return len(seq)


def parse_molwt(value: object) -> Optional[float]:
    """第一版只做最基础的数值提取。"""
    return to_float_safe(value)


# =========================
# 5. IC50 解析与标准化
# =========================
def normalize_unit_text(text: str) -> str:
    """
    把单位字符串规范化成便于解析的形式。

    例如：
    - μM / µM -> uM
    - ug/ml / μg/ml / µg/ml -> ug/mL
    - mg/ml -> mg/mL
    """
    text = text.strip()
    text = text.replace("μ", "u").replace("µ", "u")
    text = text.replace("UM", "uM").replace("uM", "uM")
    text = text.replace("NM", "nM").replace("MM", "mM")
    text = text.replace("ug/ml", "ug/mL").replace("UG/ML", "ug/mL")
    text = text.replace("mg/ml", "mg/mL").replace("MG/ML", "mg/mL")
    text = text.replace("ng/ml", "ng/mL").replace("NG/ML", "ng/mL")
    return text


def parse_ic50(raw_value: object) -> Dict[str, object]:
    """
    解析 IC50 原始文本。

    输出字段说明：
    - ic50_type:
        exact           单一数值 + 明确摩尔浓度单位 + 无关系符
        threshold       带 > / < 等关系符的摩尔浓度
        range           数值范围，例如 0.32-14 uM
        mass_unit       质量浓度，例如 ug/mL, mg/mL
        percent         百分比抑制率，不是精确 IC50
        suspicious_unit 可疑单位，例如 uM/L
        dirty           明显脏值或难以可靠解释的文本
        missing         缺失
    - ic50_exact_flag:
        只有 exact 才为 1，其余都为 0
    - ic50_uM:
        仅对 exact / threshold 且单位为明确摩尔浓度时做换算
        其余情况留空
    """
    result = {
        "ic50_type": "missing",
        "ic50_relation": None,
        "ic50_value": None,
        "ic50_low": None,
        "ic50_high": None,
        "ic50_unit_raw": None,
        "ic50_uM": None,
        "ic50_exact_flag": 0,
        "ic50_parse_note": None,
    }

    text = clean_text(raw_value)
    if text is None:
        return result

    result["ic50_unit_raw"] = text
    norm = normalize_unit_text(text)
    compact = norm.replace(" ", "")

    # 百分比：不是 IC50 点值
    if "%" in compact:
        result["ic50_type"] = "percent"
        m = re.search(r"([0-9]+(?:\.[0-9]+)?)", compact)
        if m:
            result["ic50_value"] = float(m.group(1))
        result["ic50_parse_note"] = "percent_inhibition_or_percent_like_value"
        return result

    # 可疑单位：如 uM/L、mM/L
    if re.search(r"(?:nM|uM|mM|M)/L", compact):
        result["ic50_type"] = "suspicious_unit"
        num = re.search(r"([<>≤≥])?\s*([0-9]+(?:\.[0-9]+)?)", norm)
        if num:
            result["ic50_relation"] = num.group(1) if num.group(1) else "="
            result["ic50_value"] = float(num.group(2))
        result["ic50_parse_note"] = "unit_like_molar_per_liter_kept_as_suspicious"
        return result

    # 明显重复数字脏值：例如 26.2 26.2
    m_dup = re.fullmatch(r"([0-9]+(?:\.[0-9]+)?)\s+\1", norm)
    if m_dup:
        result["ic50_type"] = "dirty"
        result["ic50_value"] = float(m_dup.group(1))
        result["ic50_parse_note"] = "duplicated_number_text"
        return result

    # 范围值：例如 0.32-14 uM
    m_range = re.fullmatch(
        r"([0-9]+(?:\.[0-9]+)?)\s*[-~]\s*([0-9]+(?:\.[0-9]+)?)\s*(nM|uM|mM|M)",
        norm,
    )
    if m_range:
        result["ic50_type"] = "range"
        result["ic50_low"] = float(m_range.group(1))
        result["ic50_high"] = float(m_range.group(2))
        result["ic50_unit_raw"] = m_range.group(3)
        result["ic50_parse_note"] = "range_value_not_treated_as_exact"
        return result

    # 质量浓度：ug/mL、mg/mL、ng/mL
    m_mass = re.fullmatch(r"([<>≤≥])?\s*([0-9]+(?:\.[0-9]+)?)\s*(ng/mL|ug/mL|mg/mL)", norm)
    if m_mass:
        result["ic50_type"] = "mass_unit"
        result["ic50_relation"] = m_mass.group(1) if m_mass.group(1) else "="
        result["ic50_value"] = float(m_mass.group(2))
        result["ic50_unit_raw"] = m_mass.group(3)
        result["ic50_parse_note"] = "mass_concentration_kept_without_uM_conversion"
        return result

    # 明确摩尔浓度阈值：>1500 uM
    m_threshold = re.fullmatch(r"([<>≤≥])\s*([0-9]+(?:\.[0-9]+)?)\s*(nM|uM|mM|M)", norm)
    if m_threshold:
        relation = m_threshold.group(1)
        value = float(m_threshold.group(2))
        unit = m_threshold.group(3)
        factor = MOLAR_UNIT_TO_UM_FACTOR[unit.lower()]

        result["ic50_type"] = "threshold"
        result["ic50_relation"] = relation
        result["ic50_value"] = value
        result["ic50_unit_raw"] = unit
        result["ic50_uM"] = value * factor
        result["ic50_parse_note"] = "molar_threshold_value"
        return result

    # 明确摩尔浓度点值：15.7 uM / 0.43 mM
    m_exact = re.fullmatch(r"([0-9]+(?:\.[0-9]+)?)\s*(nM|uM|mM|M)", norm)
    if m_exact:
        value = float(m_exact.group(1))
        unit = m_exact.group(2)
        factor = MOLAR_UNIT_TO_UM_FACTOR[unit.lower()]

        result["ic50_type"] = "exact"
        result["ic50_relation"] = "="
        result["ic50_value"] = value
        result["ic50_unit_raw"] = unit
        result["ic50_uM"] = value * factor
        result["ic50_exact_flag"] = 1
        result["ic50_parse_note"] = "exact_molar_value"
        return result

    # 其余未覆盖情况：一律保留但标为 dirty
    result["ic50_type"] = "dirty"
    result["ic50_parse_note"] = "unparsed_or_ambiguous_text"
    return result


# =========================
# 6. 合并与构建总表
# =========================
def coalesce_text(preferred: object, fallback: object) -> Optional[str]:
    """优先取 preferred；为空则退回 fallback。"""
    a = clean_text(preferred)
    b = clean_text(fallback)
    return a if a is not None else b


def deduplicate_by_id_keep_first(rows: List[Dict[str, Optional[str]]], tag: str) -> Tuple[List[Dict[str, Optional[str]]], int]:
    """
    按 id 去重，只保留每个 id 的第一条记录。

    说明：
    - 这里不是按 sequence 去重
    - 只是防止原表里同一个 id 重复出现导致后面 one-to-many 合并
    """
    seen = set()
    out: List[Dict[str, Optional[str]]] = []
    dup_count = 0
    for row in rows:
        row_id = row.get("id")
        if row_id in seen:
            dup_count += 1
            continue
        seen.add(row_id)
        out.append(row)
    if dup_count > 0:
        print(f"[警告] {tag} 中检测到 {dup_count} 个重复 id。脚本将保留每个 id 的第一条记录。")
    return out, dup_count


def build_master_table(base_rows: List[Dict[str, Optional[str]]], ic50_rows: List[Dict[str, Optional[str]]]) -> List[Dict[str, object]]:
    """
    把 small + long 作为主体表，再用 ic50 表补充信息，最终生成一张总表。

    处理原则：
    1) 记录主体来自 small + long
    2) ic50 表只做补充 / 校正用途，不新增主体记录
    3) 若 base 与 ic50 表的 ic50 文本都存在且不同，则优先使用 ic50 表的值
    4) 其他字段优先保留 base 表值，缺失时再回填 ic50 表值
    """
    ic50_map = {row["id"]: row for row in ic50_rows if row.get("id") is not None}

    merged_rows: List[Dict[str, object]] = []

    for base in base_rows:
        extra = ic50_map.get(base.get("id"))
        row: Dict[str, object] = {}

        # 主体识别信息
        row["id"] = base.get("id")
        row["table_from"] = base.get("table_from")
        row["appears_in_ic50_table"] = 1 if extra else 0

        # sequence / len / molwt：优先用 base，缺失时回填 ic50 表
        row["sequence_raw"] = coalesce_text(base.get("seq"), extra.get("seq") if extra else None)
        row["len_raw"] = coalesce_text(base.get("len"), extra.get("len") if extra else None)
        row["molwt_raw"] = coalesce_text(base.get("molwt"), extra.get("molwt") if extra else None)

        # IC50：优先取 ic50 专门表的值；若没有，再用 base 表中的值
        row["ic50_raw_base"] = clean_text(base.get("ic50"))
        row["ic50_raw_ic50tbl"] = clean_text(extra.get("ic50")) if extra else None
        row["ic50_raw"] = coalesce_text(
            extra.get("ic50") if extra else None,
            base.get("ic50"),
        )

        base_ic50 = clean_text(base.get("ic50"))
        extra_ic50 = clean_text(extra.get("ic50")) if extra else None
        row["ic50_merge_conflict_flag"] = 1 if (base_ic50 is not None and extra_ic50 is not None and base_ic50 != extra_ic50) else 0

        # 其他原始上下文字段：保留原样为主，缺失时从 ic50 表回填
        for col in ["source", "mice", "method", "assay", "bitter", "pi", "bp"]:
            row[f"{col}_raw"] = coalesce_text(base.get(col), extra.get(col) if extra else None)

        # Sequence 处理
        row["sequence_clean"] = clean_sequence(row["sequence_raw"])
        seq_valid = is_valid_canonical_sequence(row["sequence_clean"])
        row["sequence_valid"] = 1 if seq_valid else 0
        row["sequence_has_noncanonical_char"] = 0 if row["sequence_clean"] is None else (0 if seq_valid else 1)

        # Length 处理
        row["len_value_raw"] = to_float_safe(row["len_raw"])
        row["len_clean"] = compute_length(row["sequence_clean"])
        if row["len_value_raw"] is None or row["len_clean"] is None:
            row["len_match"] = None
        else:
            row["len_match"] = 1 if int(row["len_value_raw"]) == int(row["len_clean"]) else 0

        # MolWt 处理
        row["molwt_value"] = parse_molwt(row["molwt_raw"])

        # IC50 处理
        row.update(parse_ic50(row["ic50_raw"]))

        merged_rows.append(row)

    # ===== 重复序列、冲突值、顺序翻转 =====
    seq_counter = Counter(
        row["sequence_clean"]
        for row in merged_rows
        if isinstance(row.get("sequence_clean"), str) and row.get("sequence_clean")
    )

    exact_values_by_seq: Dict[str, set] = defaultdict(set)
    for row in merged_rows:
        seq = row.get("sequence_clean")
        if isinstance(seq, str) and row.get("ic50_exact_flag") == 1 and row.get("ic50_uM") is not None:
            exact_values_by_seq[seq].add(row["ic50_uM"])

    existing_sequences = set(seq_counter.keys())

    for row in merged_rows:
        seq = row.get("sequence_clean")
        if isinstance(seq, str) and seq:
            dup_count = seq_counter.get(seq, 0)
            row["sequence_dup_count"] = dup_count
            row["is_duplicate_sequence"] = 1 if dup_count > 1 else 0
            row["reverse_sequence"] = seq[::-1]
            row["reverse_exists_flag"] = 1 if row["reverse_sequence"] in existing_sequences else 0
            row["n_terminal_aa"] = seq[0] if len(seq) >= 1 else None
            row["c_terminal_aa"] = seq[-1] if len(seq) >= 1 else None
            exact_unique_count = len(exact_values_by_seq.get(seq, set()))
            row["exact_ic50_unique_count_by_sequence"] = exact_unique_count
            row["exact_ic50_conflict_flag"] = 1 if exact_unique_count > 1 else 0
        else:
            row["sequence_dup_count"] = 0
            row["is_duplicate_sequence"] = 0
            row["reverse_sequence"] = None
            row["reverse_exists_flag"] = 0
            row["n_terminal_aa"] = None
            row["c_terminal_aa"] = None
            row["exact_ic50_unique_count_by_sequence"] = 0
            row["exact_ic50_conflict_flag"] = 0

        row["is_dipeptide"] = 1 if row.get("len_clean") == 2 else 0
        row["is_tripeptide"] = 1 if row.get("len_clean") == 3 else 0

    return merged_rows


# =========================
# 7. 输出与摘要
# =========================
OUTPUT_COLUMNS = [
    # 主体识别信息
    "id",
    "table_from",
    "appears_in_ic50_table",
    # sequence / length / molwt
    "sequence_raw",
    "sequence_clean",
    "sequence_valid",
    "sequence_has_noncanonical_char",
    "len_raw",
    "len_value_raw",
    "len_clean",
    "len_match",
    "molwt_raw",
    "molwt_value",
    # activity
    "ic50_raw_base",
    "ic50_raw_ic50tbl",
    "ic50_raw",
    "ic50_merge_conflict_flag",
    "ic50_type",
    "ic50_relation",
    "ic50_value",
    "ic50_low",
    "ic50_high",
    "ic50_unit_raw",
    "ic50_uM",
    "ic50_exact_flag",
    "ic50_parse_note",
    # 原始上下文字段
    "source_raw",
    "mice_raw",
    "method_raw",
    "assay_raw",
    "bitter_raw",
    "pi_raw",
    "bp_raw",
    # 后续分析辅助列
    "sequence_dup_count",
    "is_duplicate_sequence",
    "exact_ic50_unique_count_by_sequence",
    "exact_ic50_conflict_flag",
    "reverse_sequence",
    "reverse_exists_flag",
    "n_terminal_aa",
    "c_terminal_aa",
    "is_dipeptide",
    "is_tripeptide",
]


def write_csv(rows: List[Dict[str, object]], output_path: Path) -> None:
    """把总表写出为 UTF-8 BOM CSV，方便 Windows/Excel 打开。"""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=OUTPUT_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            out_row = {}
            for col in OUTPUT_COLUMNS:
                val = row.get(col)
                out_row[col] = "" if val is None else val
            writer.writerow(out_row)


def print_summary(rows: List[Dict[str, object]]) -> None:
    """
    在控制台打印一份简洁摘要，方便你确认脚本处理结果是否大体合理。
    这就是你当前阶段够用的“轻量 QC”。
    """
    total_rows = len(rows)
    valid_seq = sum(1 for r in rows if r.get("sequence_valid") == 1)
    dipep = sum(1 for r in rows if r.get("is_dipeptide") == 1)
    tripep = sum(1 for r in rows if r.get("is_tripeptide") == 1)
    exact_ic50 = sum(1 for r in rows if r.get("ic50_exact_flag") == 1)
    dup_seq = sum(1 for r in rows if r.get("is_duplicate_sequence") == 1)
    exact_conflict = sum(1 for r in rows if r.get("exact_ic50_conflict_flag") == 1)

    print("\n========== AHTPDB 总表摘要 ==========")
    print(f"总记录数: {total_rows}")
    print(f"有效标准氨基酸序列数: {valid_seq}")
    print(f"二肽记录数: {dipep}")
    print(f"三肽记录数: {tripep}")
    print(f"精确 IC50 记录数: {exact_ic50}")
    print(f"重复序列记录数: {dup_seq}")
    print(f"存在 exact IC50 冲突的记录数: {exact_conflict}")

    ic50_type_counts = Counter(r.get("ic50_type") or "missing" for r in rows)
    print("\nIC50 类型分布:")
    for k, v in ic50_type_counts.most_common():
        print(f"  {k}: {v}")

    unit_counts = Counter(r.get("ic50_unit_raw") or "missing" for r in rows)
    print("\nIC50 原始单位 / 文本前 15 类分布:")
    for k, v in unit_counts.most_common(15):
        print(f"  {k}: {v}")
    print("====================================\n")


# =========================
# 8. 命令行入口
# =========================
def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="清洗并标准化 AHTPDB raw_tables 三张表，输出一张总表 CSV（纯标准库版）")
    parser.add_argument(
        "--project-root",
        type=str,
        default=None,
        help="PepDB 项目根目录，例如 E:/MYS/PepDB。若不传，则尝试根据脚本位置自动推断。",
    )
    parser.add_argument(
        "--raw-dir",
        type=str,
        default=None,
        help="AHTPDB raw_tables 目录。若不传，则默认使用 <project_root>/DB/raw/ace/databases/AHTPDB/raw_tables",
    )
    parser.add_argument(
        "--output-path",
        type=str,
        default=None,
        help="输出 CSV 文件路径。若不传，则默认输出到 <project_root>/DB/standardized/ace/ahtpdb/ahtpdb_master_clean.csv",
    )
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()

    script_file = Path(__file__)
    project_root = Path(args.project_root) if args.project_root else infer_project_root(script_file)

    raw_dir = (
        Path(args.raw_dir)
        if args.raw_dir
        else project_root / "DB" / "raw" / "ace" / "databases" / "AHTPDB" / "raw_tables"
    )

    output_path = (
        Path(args.output_path)
        if args.output_path
        else project_root / "DB" / "standardized" / "ace" / "ahtpdb" / "ahtpdb_master_clean.csv"
    )

    # 读取三张原表
    small_rows = read_raw_table(raw_dir / SMALL_FILE, table_from="small")
    long_rows = read_raw_table(raw_dir / LONG_FILE, table_from="long")
    ic50_rows = read_raw_table(raw_dir / IC50_FILE, table_from="ic50")

    # small + long 作为主体
    base_rows = small_rows + long_rows

    # 仅按 id 去掉重复，避免后续 one-to-many 合并问题
    base_rows, _ = deduplicate_by_id_keep_first(base_rows, tag="small+long 主体表")
    ic50_rows, _ = deduplicate_by_id_keep_first(ic50_rows, tag="ic50 表")

    # 构建总表
    master_rows = build_master_table(base_rows=base_rows, ic50_rows=ic50_rows)

    # 输出
    write_csv(master_rows, output_path)
    print(f"[完成] 已输出总表 CSV：{output_path}")

    # 打印摘要
    print_summary(master_rows)


if __name__ == "__main__":
    main()
