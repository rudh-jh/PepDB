#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
脚本名称：
    expand_ahtpdb_unit_conversion_stdlib.py

脚本用途：
    对当前已经标准化完成的 AHTPDB 总表：
        DB/standardized/ace/ahtpdb/ahtpdb_master_clean.csv
    进行“可安全单位换算记录”的再处理。

核心目标：
    1）识别当前总表中仍未转成 uM、但其实可以安全换算的记录
    2）对这些旧记录进行单位转换
    3）将“旧记录剔除 + 转换后的新记录加入”
    4）输出一个新的、已合并更新的总表
    5）输出一张转换日志表，方便复核

当前脚本主要处理哪类数据：
    只自动处理这一类：
        - ic50_type == mass_unit
        - ic50_relation == "="
        - ic50_value 可解析为数值
        - molwt_value 可解析为正数
        - ic50_unit_raw 为 mg/mL 或 ug/mL（大小写、μ/µ/u 写法会自动兼容）

为什么只处理这类：
    因为这类记录的换算逻辑最明确、最稳妥：
        - mg/mL -> uM：uM = (mg/mL * 1,000,000) / 分子量(Da)
        - ug/mL -> uM：uM = (ug/mL * 1,000) / 分子量(Da)

为什么暂时不自动处理 μM/L 这类：
    因为当前标准化表中这类记录已被刻意标成 suspicious_unit，
    例如 11.6 μM/L 这类写法，很像可疑录入，不能在未人工复核前直接当成 μM。
    所以本脚本默认不动这类记录。

处理后的更新策略：
    - 不修改原始 raw 含义字段，例如 ic50_raw、ic50_unit_raw、ic50_value 保留原样
    - 只更新标准化结果字段，使其能进入后续“高置信标签”口径：
        * ic50_type -> exact
        * ic50_relation -> "="
        * ic50_uM -> 换算后的数值
        * ic50_exact_flag -> 1
        * ic50_parse_note -> exact_mass_value_converted_to_uM_by_molwt
    - 同时新增若干追踪字段，便于后续回查：
        * unit_conversion_applied
        * unit_conversion_from_unit
        * unit_conversion_to_unit
        * unit_conversion_formula
        * unit_conversion_note

输出文件：
    默认输出到：
        DB/standardized/ace/ahtpdb/

    生成：
        1. ahtpdb_master_clean_um_expanded.csv
           -> 更新后的新版总表（建议你后续分析用这张）

        2. ahtpdb_unit_conversion_log.csv
           -> 所有被转换记录的明细表

    如果加了参数 --overwrite-master：
        - 会先备份原表为：
            ahtpdb_master_clean.backup_before_unit_conversion.csv
        - 再用新版结果覆盖 ahtpdb_master_clean.csv

设计原则：
    - 纯 Python 标准库实现
    - 兼容 Windows 编码问题
    - 注释尽量详细
    - 不依赖 pandas / numpy / pyarrow
"""

from __future__ import annotations

import argparse
import csv
import shutil
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


# =========================================================
# 一、基础工具函数
# =========================================================

def 清洗字符串(value: Any) -> str:
    """
    将任意值安全转成字符串并去掉首尾空格。
    None -> ""
    """
    if value is None:
        return ""
    return str(value).strip()


def 判定为真(value: Any) -> bool:
    """
    将 CSV 中常见真假写法统一转为布尔值。
    支持：
        1 / 0
        true / false
        yes / no
        y / n
        t / f
    """
    text = 清洗字符串(value).lower()
    return text in {"1", "true", "t", "yes", "y"}


def 安全转浮点(value: Any) -> Optional[float]:
    """
    安全转 float。
    若为空或非法，返回 None。
    """
    text = 清洗字符串(value)
    if text == "":
        return None
    try:
        return float(text)
    except Exception:
        return None


def 格式化浮点(value: Optional[float], digits: int = 10) -> str:
    """
    将浮点数格式化成适合写回 CSV 的字符串。
    None -> ""
    """
    if value is None:
        return ""
    return f"{value:.{digits}g}"


def 规范化单位(unit_text: str) -> str:
    """
    对单位文本做统一化，便于判断是否可换算。

    处理逻辑：
    - 去掉首尾空格
    - 全部转小写
    - 将 μ / µ 统一替换成 u
    - 去掉内部空格

    例如：
        'mg/mL'   -> 'mg/ml'
        'μg/ml'   -> 'ug/ml'
        'µg/mL'   -> 'ug/ml'
        ' u g / m l ' 不会出现这种极端格式，但也尽量容错
    """
    text = 清洗字符串(unit_text).lower()
    text = text.replace("μ", "u").replace("µ", "u")
    text = text.replace(" ", "")
    return text


def 计算uM(value: float, unit_norm: str, molwt: float) -> Optional[float]:
    """
    根据原始质量浓度和分子量，计算换算后的 uM。

    参数：
        value:
            原始数值，例如 0.06（对应 0.06 mg/mL）
        unit_norm:
            已经规范化后的单位，例如 'mg/ml' / 'ug/ml'
        molwt:
            分子量，单位按 Da / g/mol 处理

    公式：
        mg/mL -> uM:
            uM = (mg/mL × 1,000,000) / MW

        ug/mL -> uM:
            uM = (ug/mL × 1,000) / MW

    返回：
        若可计算，返回换算后的 uM
        否则返回 None
    """
    if molwt <= 0:
        return None

    if unit_norm == "mg/ml":
        return (value * 1_000_000.0) / molwt

    if unit_norm == "ug/ml":
        return (value * 1_000.0) / molwt

    return None


# =========================================================
# 二、路径与 CSV 读写函数
# =========================================================

def 自动推断项目根目录(脚本路径: Path) -> Path:
    """
    若脚本位于 PepDB/SCR/scripts/ 下，则自动推断项目根目录为 PepDB/
    """
    return 脚本路径.resolve().parent.parent.parent


def 获取路径(project_root: Path) -> Tuple[Path, Path, Path]:
    """
    返回：
        1）输入主表路径
        2）更新版主表输出路径
        3）转换日志表输出路径
    """
    base_dir = project_root / "DB" / "standardized" / "ace" / "ahtpdb"
    input_csv = base_dir / "ahtpdb_master_clean.csv"
    output_csv = base_dir / "ahtpdb_master_clean_um_expanded.csv"
    log_csv = base_dir / "ahtpdb_unit_conversion_log.csv"
    return input_csv, output_csv, log_csv


def 读取CSV(filepath: Path) -> Tuple[List[Dict[str, str]], List[str], str]:
    """
    读取 CSV，并自动尝试多种编码与分隔符。

    返回：
        rows: 数据行（字典列表）
        fieldnames: 原始列名顺序
        used_encoding: 成功读取所使用的编码

    为什么要返回 fieldnames：
    因为我们后面要尽量保持原始主表列顺序不变，再把新加列追加到最后。
    """
    编码候选 = ["utf-8-sig", "utf-8", "gb18030", "gbk", "cp936", "latin1"]
    异常列表 = []

    for encoding in 编码候选:
        try:
            with filepath.open("r", encoding=encoding, newline="") as f:
                sample = f.read(4096)
                f.seek(0)

                try:
                    dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
                except Exception:
                    dialect = csv.excel

                reader = csv.DictReader(f, dialect=dialect)
                fieldnames = reader.fieldnames or []

                rows: List[Dict[str, str]] = []
                for row in reader:
                    cleaned = {}
                    for k, v in row.items():
                        cleaned[清洗字符串(k)] = 清洗字符串(v)
                    rows.append(cleaned)

                return rows, [清洗字符串(x) for x in fieldnames], encoding

        except Exception as e:
            异常列表.append(f"{encoding}: {type(e).__name__} -> {e}")

    raise RuntimeError(
        "无法成功读取输入 CSV 文件，请检查编码或文件完整性。\n"
        + "\n".join(异常列表)
    )


def 写出CSV(filepath: Path, fieldnames: List[str], rows: List[Dict[str, Any]]) -> None:
    """
    将字典列表写出为 CSV。
    使用 utf-8-sig，方便 Excel 打开中文。
    """
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with filepath.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def 检查必要列(rows: List[Dict[str, str]], required_cols: List[str]) -> None:
    """
    检查输入主表是否包含必要列。
    """
    if not rows:
        raise ValueError("输入表为空，无法继续处理。")

    actual = set(rows[0].keys())
    missing = [c for c in required_cols if c not in actual]
    if missing:
        raise ValueError(
            "输入表缺少以下必要列：\n"
            + ", ".join(missing)
            + "\n\n实际读取到的列名为：\n"
            + ", ".join(rows[0].keys())
        )


# =========================================================
# 三、核心判定与转换函数
# =========================================================

def 是否可自动转换(row: Dict[str, str]) -> Tuple[bool, str]:
    """
    判断一条记录是否满足“可安全自动转换”的条件。

    当前严格规则：
        1）ic50_type == mass_unit
        2）ic50_relation == "="
        3）ic50_value 可解析为正数
        4）molwt_value 可解析为正数
        5）ic50_unit_raw 规范化后为 mg/ml 或 ug/ml

    返回：
        (是否可转, 原因说明)
    """
    ic50_type = 清洗字符串(row.get("ic50_type", "")).lower()
    ic50_relation = 清洗字符串(row.get("ic50_relation", ""))
    ic50_value = 安全转浮点(row.get("ic50_value", ""))
    molwt_value = 安全转浮点(row.get("molwt_value", ""))
    unit_norm = 规范化单位(row.get("ic50_unit_raw", ""))

    if ic50_type != "mass_unit":
        return False, "ic50_type不是mass_unit"

    if ic50_relation != "=":
        return False, "ic50_relation不是等号"

    if ic50_value is None or ic50_value <= 0:
        return False, "ic50_value不是正数"

    if molwt_value is None or molwt_value <= 0:
        return False, "molwt_value不是正数"

    if unit_norm not in {"mg/ml", "ug/ml"}:
        return False, "单位不是mg/mL或ug/mL"

    return True, "满足自动转换条件"


def 转换并更新记录(row: Dict[str, str]) -> Tuple[Dict[str, str], Dict[str, Any]]:
    """
    对一条满足条件的旧记录进行单位转换，并返回：
        1）更新后的新记录
        2）一条转换日志

    注意：
    - 原始字段尽量保留不改，例如 ic50_raw / ic50_unit_raw / ic50_value
    - 只更新“标准化结果”相关字段
    """
    new_row = dict(row)  # 复制一份，避免直接改原数据

    original_type = 清洗字符串(row.get("ic50_type", ""))
    original_unit_raw = 清洗字符串(row.get("ic50_unit_raw", ""))
    original_ic50_uM = 清洗字符串(row.get("ic50_uM", ""))

    ic50_value = 安全转浮点(row.get("ic50_value", ""))
    molwt_value = 安全转浮点(row.get("molwt_value", ""))
    unit_norm = 规范化单位(original_unit_raw)

    converted_uM = 计算uM(ic50_value, unit_norm, molwt_value)

    if converted_uM is None:
        raise ValueError("当前记录不应进入转换函数：无法完成uM计算。")

    # 写回标准化结果字段
    # 说明：
    # - ic50_type 改为 exact，这样你后续跑标签质量分析与高置信训练集时能纳入
    # - ic50_relation 保持 "="
    # - ic50_uM 写入换算后的数值
    # - ic50_exact_flag 置为 1
    # - ic50_parse_note 更新为可追溯说明
    new_row["ic50_type"] = "exact"
    new_row["ic50_relation"] = "="
    new_row["ic50_uM"] = 格式化浮点(converted_uM, digits=12)
    new_row["ic50_exact_flag"] = "1"
    new_row["ic50_parse_note"] = "exact_mass_value_converted_to_uM_by_molwt"

    # 追加追踪字段（如果原表里没有这些列，后面会自动加到表尾）
    new_row["unit_conversion_applied"] = "1"
    new_row["unit_conversion_from_unit"] = original_unit_raw
    new_row["unit_conversion_to_unit"] = "uM"

    if unit_norm == "mg/ml":
        formula = "uM = (mg/mL * 1000000) / molwt_value"
        note = "mass_unit_mg_per_mL_converted_to_uM"
    else:
        formula = "uM = (ug/mL * 1000) / molwt_value"
        note = "mass_unit_ug_per_mL_converted_to_uM"

    new_row["unit_conversion_formula"] = formula
    new_row["unit_conversion_note"] = note

    # 构造转换日志
    log_row = {
        "原始记录ID": 清洗字符串(row.get("id", "")),
        "序列": 清洗字符串(row.get("sequence_clean", "")),
        "来源原表": 清洗字符串(row.get("table_from", "")),
        "原始IC50文本": 清洗字符串(row.get("ic50_raw", "")),
        "原始IC50类型": original_type,
        "原始IC50单位": original_unit_raw,
        "原始IC50数值": 格式化浮点(ic50_value, digits=12),
        "分子量": 格式化浮点(molwt_value, digits=12),
        "转换后IC50_uM": 格式化浮点(converted_uM, digits=12),
        "转换公式": formula,
        "转换说明": note,
        "sequence_valid": 清洗字符串(row.get("sequence_valid", "")),
        "len_match": 清洗字符串(row.get("len_match", "")),
        "原始ic50_uM字段": original_ic50_uM,
        "转换后是否可作为精确值": "是",
    }

    return new_row, log_row


# =========================================================
# 四、主流程
# =========================================================

def main():
    parser = argparse.ArgumentParser(
        description="扩展 AHTPDB 标准化总表中的可自动单位转换记录（纯标准库版）。"
    )
    parser.add_argument(
        "--project-root",
        type=str,
        default="",
        help="PepDB 项目根目录，例如 E:/MYS/PepDB。若不传，则按脚本位置自动推断。"
    )
    parser.add_argument(
        "--overwrite-master",
        action="store_true",
        help="若提供该参数，则会先备份原 ahtpdb_master_clean.csv，再用新结果覆盖它。默认不覆盖。"
    )
    args = parser.parse_args()

    script_path = Path(__file__)

    # 1）确定项目根目录
    if 清洗字符串(args.project_root) != "":
        project_root = Path(args.project_root).resolve()
    else:
        project_root = 自动推断项目根目录(script_path)

    # 2）路径
    input_csv, output_csv, log_csv = 获取路径(project_root)

    print("=" * 100)
    print("AHTPDB 单位转换扩展脚本（纯标准库版）")
    print("=" * 100)
    print(f"项目根目录：{project_root}")
    print(f"输入主表：{input_csv}")
    print(f"更新版主表输出：{output_csv}")
    print(f"转换日志输出：{log_csv}")
    print(f"是否覆盖原主表：{'是' if args.overwrite_master else '否'}")
    print("-" * 100)

    if not input_csv.exists():
        raise FileNotFoundError(f"未找到输入文件：{input_csv}")

    # 3）读取主表
    rows, original_fieldnames, used_encoding = 读取CSV(input_csv)
    print(f"成功读取记录数：{len(rows)}")
    print(f"读取编码：{used_encoding}")

    检查必要列(
        rows,
        [
            "id",
            "table_from",
            "sequence_clean",
            "sequence_valid",
            "len_match",
            "molwt_value",
            "ic50_type",
            "ic50_relation",
            "ic50_value",
            "ic50_unit_raw",
            "ic50_uM",
            "ic50_exact_flag",
            "ic50_parse_note",
        ]
    )

    print("检测到的列名如下：")
    print(original_fieldnames)
    print("-" * 100)

    # 4）遍历记录，执行转换
    updated_rows: List[Dict[str, str]] = []
    conversion_logs: List[Dict[str, Any]] = []

    total_rows = 0
    convertible_rows = 0
    converted_rows = 0
    skipped_rows = 0

    for row in rows:
        total_rows += 1

        can_convert, reason = 是否可自动转换(row)

        if can_convert:
            convertible_rows += 1
            try:
                new_row, log_row = 转换并更新记录(row)
                updated_rows.append(new_row)
                conversion_logs.append(log_row)
                converted_rows += 1
            except Exception as e:
                # 理论上不应发生；若发生，则保留原记录不动，并记为跳过
                old_row = dict(row)
                old_row["unit_conversion_applied"] = "0"
                old_row["unit_conversion_from_unit"] = ""
                old_row["unit_conversion_to_unit"] = ""
                old_row["unit_conversion_formula"] = ""
                old_row["unit_conversion_note"] = f"conversion_failed: {e}"
                updated_rows.append(old_row)
                skipped_rows += 1
        else:
            # 不能转的记录原样保留
            old_row = dict(row)
            old_row["unit_conversion_applied"] = "0"
            old_row["unit_conversion_from_unit"] = ""
            old_row["unit_conversion_to_unit"] = ""
            old_row["unit_conversion_formula"] = ""
            old_row["unit_conversion_note"] = f"not_converted: {reason}"
            updated_rows.append(old_row)

    # 5）准备输出列顺序
    #    尽量保持原始列顺序不变，再把新增追踪列加到末尾
    extra_fields = [
        "unit_conversion_applied",
        "unit_conversion_from_unit",
        "unit_conversion_to_unit",
        "unit_conversion_formula",
        "unit_conversion_note",
    ]

    output_fieldnames = list(original_fieldnames)
    for col in extra_fields:
        if col not in output_fieldnames:
            output_fieldnames.append(col)

    # 6）写出更新版主表
    写出CSV(output_csv, output_fieldnames, updated_rows)

    # 7）写出转换日志
    log_fieldnames = [
        "原始记录ID",
        "序列",
        "来源原表",
        "原始IC50文本",
        "原始IC50类型",
        "原始IC50单位",
        "原始IC50数值",
        "分子量",
        "转换后IC50_uM",
        "转换公式",
        "转换说明",
        "sequence_valid",
        "len_match",
        "原始ic50_uM字段",
        "转换后是否可作为精确值",
    ]
    写出CSV(log_csv, log_fieldnames, conversion_logs)

    # 8）如需覆盖原主表，则先备份再覆盖
    if args.overwrite_master:
        backup_path = input_csv.with_name("ahtpdb_master_clean.backup_before_unit_conversion.csv")
        shutil.copy2(input_csv, backup_path)
        shutil.copy2(output_csv, input_csv)
        print(f"已备份原主表：{backup_path}")
        print(f"已覆盖原主表：{input_csv}")

    # 9）控制台摘要
    print("-" * 100)
    print("处理完成。")
    print(f"总记录数：{total_rows}")
    print(f"满足自动转换条件的记录数：{convertible_rows}")
    print(f"成功完成转换的记录数：{converted_rows}")
    print(f"转换过程中失败并保留原样的记录数：{skipped_rows}")
    print("-" * 100)
    print(f"更新版主表：{output_csv}")
    print(f"转换日志表：{log_csv}")
    print("=" * 100)


if __name__ == "__main__":
    main()