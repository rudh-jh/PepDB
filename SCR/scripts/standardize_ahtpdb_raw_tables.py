#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AHTPDB raw_tables 三表清洗与标准化脚本（简化版，面向“总表”）

适用场景
--------
本脚本专门服务于 PepDB 当前阶段的目标：
- 读取 AHTPDB 的三张原始表：small / long / ic50
- 以 small + long 作为主体记录表
- 用 peptides_with_ic50 表补充 / 校正活性相关信息
- 对最关键字段进行清洗：sequence、length、molwt、ic50
- 输出一张可直接用于后续分析的 CSV 总表

本脚本刻意避免“过度工程化”：
- 不拆成复杂多表数据库
- 不删除异常值，而是尽量保留并加标记列
- 不在第一版中强行把所有质量浓度换算为 uM
- 不对 source / mice / method / assay 做过度结构化

后续你可以基于这张总表继续做：
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

如果你的脚本不放在上述位置，也可以使用命令行参数覆盖路径。

运行示例
--------
1) 若脚本已放入 PepDB/SCR/scripts/：
   python standardize_ahtpdb_raw_tables.py

2) 指定项目根目录：
   python standardize_ahtpdb_raw_tables.py --project-root E:\\MYS\\PepDB

3) 指定输出文件：
   python standardize_ahtpdb_raw_tables.py --output-path E:\\MYS\\PepDB\\DB\\standardized\\ace\\ahtpdb\\ahtpdb_master_clean.csv
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, Optional

import pandas as pd


# =========================
# 1. 基础常量设置
# =========================
# AHTPDB 原始三张表的默认文件名。
# 如果后续官方文件名变化，只需要修改这里。
SMALL_FILE = "ahtpdb_small_peptides_2_5_residues.txt"
LONG_FILE = "ahtpdb_long_peptides_6_plus_residues.txt"
IC50_FILE = "ahtpdb_peptides_with_ic50.txt"

# 原始表的标准列顺序。当前 AHTPDB 的三张表都是这套表头。
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

# 常见缺失值表达。注意：这里仅用于“工作中的判定”。
# 输出时仍会保留原始文本列，例如 sequence_raw、ic50_raw。
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

# 标准 20 种氨基酸单字母代码。
CANONICAL_AA = set("ACDEFGHIKLMNPQRSTVWY")

# 用于把不同写法统一成标准单位名称。
# 第一版只对“明确摩尔浓度点值 / 阈值”做 uM 统一，其他类型先保留。
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
    尝试根据脚本位置推断 PepDB 项目根目录。

    约定：脚本放在 PepDB/SCR/scripts/
    那么 script_file.parents[2] 就是 PepDB 根目录。

    如果未来你把脚本挪走，建议直接使用 --project-root 显式传入。
    """
    try:
        return script_file.resolve().parents[2]
    except IndexError:
        return script_file.resolve().parent


def normalize_spaces(text: str) -> str:
    """去掉首尾空格，并把中间连续空白压缩成单个空格。"""
    return re.sub(r"\s+", " ", text.strip())


def is_missing(value: object) -> bool:
    """
    判断一个值是否应视为缺失。

    说明：
    - pandas 的 NA / NaN 视为缺失
    - 常见的 ND / N.A. / 空字符串等也视为缺失
    """
    if pd.isna(value):
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
    """
    尝试把一个值转成 float。
    转换失败时返回 None，而不是报错。
    """
    text = clean_text(value)
    if text is None:
        return None
    # 兼容少量千位分隔或奇怪空格
    text = text.replace(",", "")
    try:
        return float(text)
    except ValueError:
        return None


def read_raw_table(path: Path, table_from: str) -> pd.DataFrame:
    """
    读取 AHTPDB 原始表。

    设计原则：
    - 优先按 TAB 读取，因为 AHTPDB 官方导出通常是制表符分隔
    - 若列数与预期不一致，再退回到“两个及以上空白分隔”读取
    - 所有列先按字符串读入，避免 pandas 自动把 ND、科学计数法等误处理
    """
    if not path.exists():
        raise FileNotFoundError(f"找不到输入文件：{path}")

    # 第一轮：按 tab 读取
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, engine="python")

    # 若列数不对，再尝试更宽松的空白分隔
    if list(df.columns) != EXPECTED_COLUMNS:
        df = pd.read_csv(path, sep=r"\s{2,}|\t", dtype=str, keep_default_na=False, engine="python")

    # 如果仍然不是预期表头，直接报错，防止悄悄读错列
    if list(df.columns) != EXPECTED_COLUMNS:
        raise ValueError(
            f"文件 {path.name} 的列名不是预期格式。\n"
            f"当前列名：{list(df.columns)}\n"
            f"预期列名：{EXPECTED_COLUMNS}"
        )

    # 加一个来源表标记，后续合并时很有用
    df["table_from"] = table_from

    # 统一做轻清洗，保证后续字符串操作更稳
    for col in EXPECTED_COLUMNS:
        df[col] = df[col].map(clean_text)

    return df


# =========================
# 3. Sequence / Length / MolWt 处理
# =========================
def clean_sequence(seq: object) -> Optional[str]:
    """
    清洗氨基酸序列：
    - 缺失返回 None
    - 转大写
    - 去空格
    - 不改序列顺序
    - 不自动替换可疑字符（例如 J），因为这样会引入主观假设
    """
    text = clean_text(seq)
    if text is None:
        return None
    text = text.upper().replace(" ", "")
    return text if text else None


def is_valid_canonical_sequence(seq: Optional[str]) -> bool:
    """
    判断是否为“仅包含标准 20 种氨基酸”的序列。
    说明：
    - 只要出现 J / B / Z / X / 数字 / 符号，就判为 False
    - 这样做比较保守，但适合后续二肽/三肽规律分析
    """
    if not seq:
        return False
    return set(seq).issubset(CANONICAL_AA)


def compute_length(seq: Optional[str]) -> Optional[int]:
    """根据清洗后的序列重算长度。"""
    if not seq:
        return None
    return len(seq)


def parse_molwt(value: object) -> Optional[float]:
    """
    解析分子量列。
    第一版只做最基础的数值提取，不做“按序列回算分子量”的复杂校验。
    这样既能保留总表的简洁性，也能为后续需要时留下接口。
    """
    return to_float_safe(value)


# =========================
# 4. IC50 解析与标准化
# =========================
def normalize_unit_text(text: str) -> str:
    """
    把单位字符串规范化成便于解析的形式。

    例如：
    - μM / µM -> uM
    - ug/ml / μg/ml / µg/ml -> ug/mL
    - mg/ml -> mg/mL

    注意：
    这里只做“字符串层面的统一写法”，并不意味着一定会换算。
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

    目标不是“吃掉所有边角情况”，而是：
    1) 尽量把最常见、最重要的情况识别出来
    2) 让“精确 IC50”判定清楚可靠
    3) 其余情况一律保留原始值，并分到对应类型，便于后续筛选

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

    # ---------- 4.1 百分比：不是 IC50 点值 ----------
    if "%" in compact:
        result["ic50_type"] = "percent"
        m = re.search(r"([0-9]+(?:\.[0-9]+)?)", compact)
        if m:
            result["ic50_value"] = float(m.group(1))
        result["ic50_parse_note"] = "percent_inhibition_or_percent_like_value"
        return result

    # ---------- 4.2 可疑单位：如 uM/L、mM/L ----------
    if re.search(r"(?:nM|uM|mM|M)/L", compact):
        result["ic50_type"] = "suspicious_unit"
        num = re.search(r"([<>≤≥])?\s*([0-9]+(?:\.[0-9]+)?)", norm)
        if num:
            result["ic50_relation"] = num.group(1) if num.group(1) else "="
            result["ic50_value"] = float(num.group(2))
        result["ic50_parse_note"] = "unit_like_molar_per_liter_kept_as_suspicious"
        return result

    # ---------- 4.3 明显重复数字脏值：例如 26.2 26.2 ----------
    m_dup = re.fullmatch(r"([0-9]+(?:\.[0-9]+)?)\s+\1", norm)
    if m_dup:
        result["ic50_type"] = "dirty"
        result["ic50_value"] = float(m_dup.group(1))
        result["ic50_parse_note"] = "duplicated_number_text"
        return result

    # ---------- 4.4 范围值：例如 0.32-14 uM ----------
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

    # ---------- 4.5 质量浓度：ug/mL、mg/mL、ng/mL ----------
    m_mass = re.fullmatch(r"([<>≤≥])?\s*([0-9]+(?:\.[0-9]+)?)\s*(ng/mL|ug/mL|mg/mL)", norm)
    if m_mass:
        result["ic50_type"] = "mass_unit"
        result["ic50_relation"] = m_mass.group(1) if m_mass.group(1) else "="
        result["ic50_value"] = float(m_mass.group(2))
        result["ic50_unit_raw"] = m_mass.group(3)
        result["ic50_parse_note"] = "mass_concentration_kept_without_uM_conversion"
        return result

    # ---------- 4.6 明确摩尔浓度阈值：>1500 uM ----------
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

    # ---------- 4.7 明确摩尔浓度点值：15.7 uM / 0.43 mM ----------
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

    # ---------- 4.8 其余未覆盖情况：一律保留但标为 dirty ----------
    result["ic50_type"] = "dirty"
    result["ic50_parse_note"] = "unparsed_or_ambiguous_text"
    return result


# =========================
# 5. 构建总表的核心流程
# =========================
def coalesce_text(preferred: object, fallback: object) -> Optional[str]:
    """
    优先取 preferred；若其为空，则退回 fallback。
    适合合并 base 表与 ic50 表的同名文本字段。
    """
    a = clean_text(preferred)
    b = clean_text(fallback)
    return a if a is not None else b


def build_master_table(base_df: pd.DataFrame, ic50_df: pd.DataFrame) -> pd.DataFrame:
    """
    把 small + long 作为主体表，再用 ic50 表补充信息，最终生成一张总表。

    处理原则：
    1) 记录主体来自 small + long
    2) ic50 表只做“补充 / 校正”用途，不新增主体记录
    3) 若 base 与 ic50 表的 ic50 文本都存在且不同，则优先使用 ic50 表的值
       因为它本身是“带活性值的专门表”，但同时保留冲突标记
    4) 其他字段一般优先保留 base 表值，缺失时再回填 ic50 表值
    """
    # 为了避免同名列冲突，先给 ic50 表全部加后缀
    ic50_renamed = ic50_df.rename(columns={col: f"{col}_ic50tbl" for col in ic50_df.columns})

    merged = base_df.merge(
        ic50_renamed,
        left_on="id",
        right_on="id_ic50tbl",
        how="left",
        validate="one_to_one",
    )

    # 基础合并标记
    merged["appears_in_ic50_table"] = merged["id_ic50tbl"].notna().astype(int)

    # 合并后，以 base 表为主体生成最终字段。
    # sequence / len / molwt 等主体信息一般先取 base，再补 ic50 表。
    merged["sequence_raw"] = merged.apply(lambda r: coalesce_text(r.get("seq"), r.get("seq_ic50tbl")), axis=1)
    merged["len_raw"] = merged.apply(lambda r: coalesce_text(r.get("len"), r.get("len_ic50tbl")), axis=1)
    merged["molwt_raw"] = merged.apply(lambda r: coalesce_text(r.get("molwt"), r.get("molwt_ic50tbl")), axis=1)

    # IC50 的合并规则：优先取 ic50 专门表中的值；若没有，再用 base 表中的值。
    merged["ic50_raw_base"] = merged["ic50"]
    merged["ic50_raw_ic50tbl"] = merged["ic50_ic50tbl"]
    merged["ic50_raw"] = merged.apply(lambda r: coalesce_text(r.get("ic50_ic50tbl"), r.get("ic50")), axis=1)

    # 记录 base 与 ic50 表的 IC50 文本是否冲突。
    def ic50_merge_conflict(row: pd.Series) -> int:
        base = clean_text(row.get("ic50"))
        extra = clean_text(row.get("ic50_ic50tbl"))
        if base is None or extra is None:
            return 0
        return int(base != extra)

    merged["ic50_merge_conflict_flag"] = merged.apply(ic50_merge_conflict, axis=1)

    # 其他原始上下文字段：保留原样为主，缺失时从 ic50 表回填
    for col in ["source", "mice", "method", "assay", "bitter", "pi", "bp"]:
        merged[f"{col}_raw"] = merged.apply(
            lambda r, c=col: coalesce_text(r.get(c), r.get(f"{c}_ic50tbl")), axis=1
        )

    # ========== Sequence 处理 ==========
    merged["sequence_clean"] = merged["sequence_raw"].map(clean_sequence)
    merged["sequence_valid"] = merged["sequence_clean"].map(is_valid_canonical_sequence).astype(int)
    merged["sequence_has_noncanonical_char"] = (1 - merged["sequence_valid"]).where(merged["sequence_clean"].notna(), 0)

    # ========== Length 处理 ==========
    merged["len_value_raw"] = merged["len_raw"].map(to_float_safe)
    merged["len_clean"] = merged["sequence_clean"].map(compute_length)

    def length_match(row: pd.Series) -> Optional[int]:
        raw_len = row.get("len_value_raw")
        clean_len = row.get("len_clean")
        if raw_len is None or pd.isna(raw_len) or clean_len is None or pd.isna(clean_len):
            return None
        return int(int(raw_len) == int(clean_len))

    merged["len_match"] = merged.apply(length_match, axis=1)

    # ========== MolWt 处理 ==========
    merged["molwt_value"] = merged["molwt_raw"].map(parse_molwt)

    # ========== IC50 处理 ==========
    ic50_parsed = merged["ic50_raw"].map(parse_ic50)
    ic50_parsed_df = pd.DataFrame(list(ic50_parsed))
    merged = pd.concat([merged, ic50_parsed_df], axis=1)

    # ========== 重复序列、冲突值、顺序翻转 ==========
    # 注意：这些列不是删除依据，只是后续分析的辅助信息。
    merged["sequence_dup_count"] = merged.groupby("sequence_clean")["sequence_clean"].transform("count")
    merged["sequence_dup_count"] = merged["sequence_dup_count"].fillna(0).astype(int)
    merged["is_duplicate_sequence"] = (merged["sequence_dup_count"] > 1).astype(int)

    # reverse sequence：为后续“顺序翻转对比”准备
    merged["reverse_sequence"] = merged["sequence_clean"].map(lambda x: x[::-1] if isinstance(x, str) else None)
    existing_sequences = set(merged["sequence_clean"].dropna().tolist())
    merged["reverse_exists_flag"] = merged["reverse_sequence"].map(lambda x: int(x in existing_sequences) if isinstance(x, str) else 0)

    # exact IC50 冲突：同一序列若对应多个不同 exact uM，则记为冲突
    exact_only = merged[(merged["sequence_clean"].notna()) & (merged["ic50_exact_flag"] == 1)].copy()
    exact_unique_counts = (
        exact_only.groupby("sequence_clean")["ic50_uM"]
        .nunique(dropna=True)
        .rename("exact_ic50_unique_count_by_sequence")
        .reset_index()
    )
    merged = merged.merge(exact_unique_counts, on="sequence_clean", how="left")
    merged["exact_ic50_unique_count_by_sequence"] = merged["exact_ic50_unique_count_by_sequence"].fillna(0).astype(int)
    merged["exact_ic50_conflict_flag"] = (merged["exact_ic50_unique_count_by_sequence"] > 1).astype(int)

    # ========== 一些简单分析友好列 ==========
    merged["n_terminal_aa"] = merged["sequence_clean"].map(lambda x: x[0] if isinstance(x, str) and len(x) >= 1 else None)
    merged["c_terminal_aa"] = merged["sequence_clean"].map(lambda x: x[-1] if isinstance(x, str) and len(x) >= 1 else None)
    merged["is_dipeptide"] = merged["len_clean"].map(lambda x: int(x == 2) if pd.notna(x) else 0)
    merged["is_tripeptide"] = merged["len_clean"].map(lambda x: int(x == 3) if pd.notna(x) else 0)

    # ========== 最终列顺序 ==========
    # 优先把真正分析常用的列放前面，其他辅助列放后面。
    ordered_cols = [
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

    # 只保留存在的列，避免未来上游文件小变动时因为列顺序强绑定而报错
    ordered_cols = [c for c in ordered_cols if c in merged.columns]
    merged = merged[ordered_cols].copy()

    # 为了更适合保存 CSV，把所有 NA 统一留空，而不是写成 <NA>
    merged = merged.where(pd.notna(merged), None)

    return merged


# =========================
# 6. 输出与简单统计
# =========================
def print_summary(df: pd.DataFrame) -> None:
    """
    在控制台打印一份简洁摘要，方便你确认脚本处理结果是否大体合理。
    这部分就相当于“轻量 QC”。
    """
    total_rows = len(df)
    valid_seq = int((df["sequence_valid"] == 1).sum()) if "sequence_valid" in df else 0
    dipep = int((df["is_dipeptide"] == 1).sum()) if "is_dipeptide" in df else 0
    tripep = int((df["is_tripeptide"] == 1).sum()) if "is_tripeptide" in df else 0
    exact_ic50 = int((df["ic50_exact_flag"] == 1).sum()) if "ic50_exact_flag" in df else 0
    dup_seq = int((df["is_duplicate_sequence"] == 1).sum()) if "is_duplicate_sequence" in df else 0
    exact_conflict = int((df["exact_ic50_conflict_flag"] == 1).sum()) if "exact_ic50_conflict_flag" in df else 0

    print("\n========== AHTPDB 总表摘要 ==========")
    print(f"总记录数: {total_rows}")
    print(f"有效标准氨基酸序列数: {valid_seq}")
    print(f"二肽记录数: {dipep}")
    print(f"三肽记录数: {tripep}")
    print(f"精确 IC50 记录数: {exact_ic50}")
    print(f"重复序列记录数: {dup_seq}")
    print(f"存在 exact IC50 冲突的记录数: {exact_conflict}")

    if "ic50_type" in df.columns:
        print("\nIC50 类型分布:")
        type_counts = df["ic50_type"].fillna("missing").value_counts(dropna=False)
        for k, v in type_counts.items():
            print(f"  {k}: {v}")

    if "ic50_unit_raw" in df.columns:
        print("\nIC50 原始单位 / 文本前 15 类分布:")
        unit_counts = df["ic50_unit_raw"].fillna("missing").value_counts(dropna=False).head(15)
        for k, v in unit_counts.items():
            print(f"  {k}: {v}")
    print("====================================\n")


# =========================
# 7. 命令行入口
# =========================
def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="清洗并标准化 AHTPDB raw_tables 三张表，输出一张总表 CSV")
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

    output_path.parent.mkdir(parents=True, exist_ok=True)

    # ---------- 读取三张原表 ----------
    small_df = read_raw_table(raw_dir / SMALL_FILE, table_from="small")
    long_df = read_raw_table(raw_dir / LONG_FILE, table_from="long")
    ic50_df = read_raw_table(raw_dir / IC50_FILE, table_from="ic50")

    # ---------- small + long 作为主体 ----------
    base_df = pd.concat([small_df, long_df], ignore_index=True)

    # 简单去掉 base 中“完全重复行”。
    # 说明：
    # - 这不是序列层面的去重
    # - 只是防止极端情况下原表出现完全一样的重复记录
    base_df = base_df.drop_duplicates().copy()

    # ic50 表也做同样处理，避免 merge 时意外 one-to-many
    ic50_df = ic50_df.drop_duplicates().copy()

    # 如果 ic50 表里 id 有重复，保留第一条，同时提醒用户。
    # 这里不直接报错，是为了让脚本尽量跑完；但控制台会提示你复核。
    if ic50_df["id"].duplicated().any():
        dup_count = int(ic50_df["id"].duplicated().sum())
        print(f"[警告] ic50 表中检测到 {dup_count} 个重复 id。脚本将保留每个 id 的第一条记录。")
        ic50_df = ic50_df.drop_duplicates(subset=["id"], keep="first").copy()

    if base_df["id"].duplicated().any():
        dup_count = int(base_df["id"].duplicated().sum())
        print(f"[警告] small+long 主体表中检测到 {dup_count} 个重复 id。脚本将保留每个 id 的第一条记录。")
        base_df = base_df.drop_duplicates(subset=["id"], keep="first").copy()

    # ---------- 构建总表 ----------
    master_df = build_master_table(base_df=base_df, ic50_df=ic50_df)

    # ---------- 输出 ----------
    master_df.to_csv(output_path, index=False, encoding="utf-8-sig")
    print(f"[完成] 已输出总表 CSV：{output_path}")

    # ---------- 打印摘要 ----------
    print_summary(master_df)


if __name__ == "__main__":
    main()
