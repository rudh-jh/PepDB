#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
脚本名称：
    build_highconf_short_peptide_worksets_stdlib.py

脚本用途：
    从已经标准化完成的 AHTPDB 总表 ahtpdb_master_clean.csv 中，
    自动构建“高置信短肽工作集”，用于后续：
    1）二肽 / 三肽构效分析
    2）N/C 端位点富集分析
    3）顺序翻转对比
    4）后续可解释预测模型的训练集准备

这次改动说明：
    你希望把脚本所依赖的“输入文件 / 输出目录”封装到脚本前部，
    以后只改几行路径变量即可，不用翻主程序逻辑。

    因此本脚本新增了“路径配置区”：
        - MANUAL_PROJECT_ROOT
        - MANUAL_INPUT_CSV
        - MANUAL_OUTPUT_DIR
        - PRIORITY_USE_MANUAL_PATHS

    当 PRIORITY_USE_MANUAL_PATHS = True 时：
        脚本优先使用文件前部你手动配置好的路径。

    这样后续你如果换表、换目录，只需要改前面这几项，
    其余业务逻辑保持不变。

脚本会生成三张表：
    1）高置信二三肽记录表.csv
       - 一行一条记录
       - 作用：保留所有高置信二肽/三肽记录，方便做记录级分析

    2）高置信二三肽序列聚合表.csv
       - 一行一个 sequence_clean
       - 作用：把同一序列的多条记录聚合起来，看 IC50 冲突、来源差异、assay 差异

    3）高置信二三肽稳定序列表.csv
       - 在“序列聚合表”基础上，进一步筛出“稳定序列”
       - 作用：作为后续更干净的二/三肽构效分析与第一版预测建模输入

脚本特点：
    - 纯 Python 标准库实现
    - 不依赖 pandas / numpy / pyarrow
    - 输出 CSV 编码为 utf-8-sig，方便 Excel 直接打开
    - 输出表名、列名全部为中文
    - 注释尽量详细，便于你后续自己改
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path
from statistics import median
from typing import Any, Dict, Iterable, List, Optional, Tuple


# =========================================================
# 一、路径配置区（你后续最常改这里）
# =========================================================
# 使用说明：
# 1）如果 PRIORITY_USE_MANUAL_PATHS = True，则脚本优先使用下面的手动路径。
# 2）如果你把它改成 False，则脚本会退回到“按项目根目录自动推导输入/输出路径”的逻辑。
# 3）你以后最常改的通常只有：
#       - MANUAL_INPUT_CSV
#       - MANUAL_OUTPUT_DIR
#    这样就不用去改主逻辑。

PRIORITY_USE_MANUAL_PATHS = True

# PepDB 项目根目录（可选）
MANUAL_PROJECT_ROOT = r"E:\MYS\PepDB"

# 输入文件：标准化后的 AHTPDB 总表
MANUAL_INPUT_CSV = r"E:\MYS\PepDB\DB\standardized\ace\ahtpdb\ahtpdb_master_clean_um_expanded.csv"

# 输出目录：高置信短肽工作集目录
MANUAL_OUTPUT_DIR = r"E:\MYS\PepDB\DB\worksets\ace\ahtpdb\高置信短肽工作集"


# =========================================================
# 二、基础工具函数
# =========================================================

def 清洗字符串(value: Any) -> str:
    """
    将任意对象安全转成字符串，并去掉首尾空白字符。
    None -> ""
    """
    if value is None:
        return ""
    return str(value).strip()


def 判定为真(value: Any) -> bool:
    """
    将 CSV 中常见的真假表达统一为布尔值。

    兼容：
    - 1 / 0
    - True / False
    - yes / no
    - y / n
    - t / f
    """
    text = 清洗字符串(value).lower()
    return text in {"1", "true", "t", "yes", "y"}


def 安全转浮点(value: Any) -> Optional[float]:
    """
    安全地把某个值转成浮点数。
    若为空、非法、不可解析，则返回 None。
    """
    text = 清洗字符串(value)
    if text == "":
        return None
    try:
        return float(text)
    except Exception:
        return None


def 格式化浮点(value: Optional[float], digits: int = 6) -> str:
    """
    将浮点数格式化为适合写入 CSV 的字符串。
    - None -> ""
    - 否则保留一定有效位数，尽量避免过长小数和夸张科学计数法
    """
    if value is None:
        return ""
    return f"{value:.{digits}g}"


def 百分比字符串(分子: int, 分母: int) -> str:
    """
    将分子/分母转成百分比字符串，保留两位小数。
    """
    if 分母 == 0:
        return "0.00%"
    return f"{(分子 / 分母) * 100:.2f}%"


def 去重并排序(values: Iterable[str], sep: str = " | ") -> str:
    """
    对一组字符串：
    1）去空
    2）去重
    3）排序
    4）用 sep 拼接

    常用于：
    - 来源列表
    - assay 列表
    - 记录 ID 列表
    """
    cleaned = [清洗字符串(v) for v in values if 清洗字符串(v) != ""]
    return sep.join(sorted(set(cleaned)))


def 对数10(value: Optional[float]) -> Optional[float]:
    """
    安全计算 log10。
    只有 value > 0 时才计算，否则返回 None。
    """
    if value is None or value <= 0:
        return None
    try:
        return math.log10(value)
    except Exception:
        return None


# =========================================================
# 三、路径与 CSV 读写函数
# =========================================================

def 自动推断项目根目录(脚本路径: Path) -> Path:
    """
    如果脚本放在 PepDB/SCR/scripts/ 目录下，
    则项目根目录可以自动推断为 PepDB。
    """
    return 脚本路径.resolve().parent.parent.parent


def 获取输入输出路径(project_root: Path) -> Tuple[Path, Path]:
    """
    根据项目根目录，计算默认：
    1）输入 CSV 路径
    2）输出目录路径

    注意：
    这个函数主要用于“你关闭手动路径优先”时的后备逻辑。
    """
    input_csv = project_root / "DB" / "standardized" / "ace" / "ahtpdb" / "ahtpdb_master_clean.csv"
    output_dir = project_root / "DB" / "worksets" / "ace" / "ahtpdb" / "高置信短肽工作集"
    return input_csv, output_dir


def 解析最终输入输出路径(args_project_root: str, script_path: Path) -> Tuple[Path, Path, Path]:
    """
    统一解析最终实际使用的：
        - 项目根目录
        - 输入文件路径
        - 输出目录路径

    优先级规则：
    1）如果 PRIORITY_USE_MANUAL_PATHS = True，则优先使用脚本前部配置的手动路径
    2）否则，如果命令行传了 --project-root，则按该根目录推导默认输入输出路径
    3）否则，自动按脚本位置推断项目根目录，再推导默认输入输出路径
    """
    if PRIORITY_USE_MANUAL_PATHS:
        if 清洗字符串(MANUAL_PROJECT_ROOT) != "":
            project_root = Path(MANUAL_PROJECT_ROOT).resolve()
        elif 清洗字符串(args_project_root) != "":
            project_root = Path(args_project_root).resolve()
        else:
            project_root = 自动推断项目根目录(script_path)

        input_csv = Path(MANUAL_INPUT_CSV).resolve()
        output_dir = Path(MANUAL_OUTPUT_DIR).resolve()
        return project_root, input_csv, output_dir

    if 清洗字符串(args_project_root) != "":
        project_root = Path(args_project_root).resolve()
    else:
        project_root = 自动推断项目根目录(script_path)

    input_csv, output_dir = 获取输入输出路径(project_root)
    return project_root, input_csv, output_dir


def 读取CSV(filepath: Path) -> List[Dict[str, str]]:
    """
    读取 CSV 文件，并自动尝试多种常见编码。

    为什么要这样做：
    你的 CSV 可能来自 Windows 环境、Excel、脚本导出等不同来源，
    编码不一定稳定是 utf-8-sig，所以这里采用“多编码尝试”的方式。

    尝试顺序：
    1. utf-8-sig
    2. utf-8
    3. gb18030
    4. gbk
    5. cp936
    6. latin1（最后兜底）

    同时还会尝试自动识别分隔符（逗号、制表符、分号）。
    """
    编码候选列表 = ["utf-8-sig", "utf-8", "gb18030", "gbk", "cp936", "latin1"]
    异常列表 = []

    for 编码 in 编码候选列表:
        try:
            with filepath.open("r", encoding=编码, newline="") as f:
                sample = f.read(4096)
                f.seek(0)

                try:
                    dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
                except Exception:
                    dialect = csv.excel

                reader = csv.DictReader(f, dialect=dialect)
                rows = []
                for row in reader:
                    cleaned_row = {}
                    for k, v in row.items():
                        cleaned_row[清洗字符串(k)] = 清洗字符串(v)
                    rows.append(cleaned_row)

                print(f"成功使用编码读取文件：{编码}")
                return rows

        except Exception as e:
            异常列表.append(f"{编码}: {type(e).__name__} -> {e}")

    raise RuntimeError(
        "无法成功读取输入 CSV 文件，请检查编码或文件是否损坏。\n"
        + "\n".join(异常列表)
    )


def 写出CSV(filepath: Path, fieldnames: List[str], rows: List[Dict[str, Any]]) -> None:
    """
    将结果写出为 CSV。

    特点：
    - 自动创建父目录
    - 使用 utf-8-sig，方便 Excel 正常显示中文
    - 如果某行多了字段，自动忽略，不报错
    """
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with filepath.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def 检查必要列(rows: List[Dict[str, str]], required_cols: List[str]) -> None:
    """
    检查输入表中是否包含后续分析必须使用的列。
    如果缺列，则尽早报错，避免后面在中间流程才崩。
    """
    if not rows:
        raise ValueError("输入 CSV 没有任何数据行。")

    actual_cols = set(rows[0].keys())
    missing = [c for c in required_cols if c not in actual_cols]
    if missing:
        raise ValueError(
            "输入表缺少以下必要列：\n"
            + ", ".join(missing)
            + "\n\n实际读取到的列名为：\n"
            + ", ".join(rows[0].keys())
        )


# =========================================================
# 四、核心筛选与判定函数
# =========================================================

def 是否高置信短肽记录(row: Dict[str, str]) -> bool:
    """
    判定一条记录是否属于“高置信短肽记录”。

    当前筛选条件：
    1）sequence_valid = 1
    2）len_match = 1
    3）ic50_exact_flag = 1
    4）ic50_type = exact
    5）ic50_uM 可解析为数值
    6）is_dipeptide = 1 或 is_tripeptide = 1
    """
    sequence_valid = 判定为真(row.get("sequence_valid", ""))
    len_match = 判定为真(row.get("len_match", ""))
    ic50_exact_flag = 判定为真(row.get("ic50_exact_flag", ""))
    ic50_type = 清洗字符串(row.get("ic50_type", "")).lower()
    ic50_uM = 安全转浮点(row.get("ic50_uM", ""))

    is_dipeptide = 判定为真(row.get("is_dipeptide", ""))
    is_tripeptide = 判定为真(row.get("is_tripeptide", ""))

    return (
        sequence_valid
        and len_match
        and ic50_exact_flag
        and ic50_type == "exact"
        and ic50_uM is not None
        and (is_dipeptide or is_tripeptide)
    )


def 短肽类型(row: Dict[str, str]) -> str:
    """
    返回该记录属于二肽、三肽还是其他。
    按当前规则，理论上只会进入二肽或三肽。
    """
    if 判定为真(row.get("is_dipeptide", "")):
        return "二肽"
    if 判定为真(row.get("is_tripeptide", "")):
        return "三肽"
    return "其他"


def 稳定性等级(精确值列表: List[float]) -> Tuple[str, str, Optional[float]]:
    """
    根据同一序列的多个精确 IC50 值，判断其稳定性等级。

    返回：
        (稳定性等级, 判定说明, 最大最小倍数)

    规则（第一版，偏保守）：
    1）如果只有 1 个精确值：
        -> 稳定序列
    2）如果有多个精确值，但 max/min <= 3：
        -> 稳定序列
    3）如果 3 < max/min <= 10：
        -> 条件敏感序列
    4）如果 max/min > 10：
        -> 高冲突序列
    """
    if not 精确值列表:
        return "无法判定", "无精确IC50值", None

    unique_values = sorted(set(精确值列表))
    if len(unique_values) == 1:
        return "稳定序列", "仅有1个精确IC50值", 1.0

    min_v = min(unique_values)
    max_v = max(unique_values)

    if min_v <= 0:
        return "无法判定", "存在非正IC50值，无法计算倍数", None

    ratio = max_v / min_v

    if ratio <= 3:
        return "稳定序列", "多个精确IC50值，但最大最小倍数<=3", ratio
    if ratio <= 10:
        return "条件敏感序列", "多个精确IC50值，且3<最大最小倍数<=10", ratio
    return "高冲突序列", "多个精确IC50值，且最大最小倍数>10", ratio


# =========================================================
# 五、构建表1：高置信二三肽记录表
# =========================================================

def 构建高置信二三肽记录表(rows: List[Dict[str, str]]) -> List[Dict[str, Any]]:
    """
    从总表中筛选出高置信短肽记录，并重命名为中文列。

    这张表是“记录级”表：
    - 一行 = 一条高置信二/三肽记录
    - 不做序列级聚合
    - 方便后续做记录级分析与回溯原记录
    """
    result = []

    for row in rows:
        if not 是否高置信短肽记录(row):
            continue

        seq = 清洗字符串(row.get("sequence_clean", ""))
        ic50_uM = 安全转浮点(row.get("ic50_uM", ""))
        log_ic50 = 对数10(ic50_uM)

        result.append({
            "原始记录ID": 清洗字符串(row.get("id", "")),
            "来源原表": 清洗字符串(row.get("table_from", "")),
            "序列": seq,
            "短肽类型": 短肽类型(row),
            "肽长": 清洗字符串(row.get("len_clean", "")),
            "N端氨基酸": 清洗字符串(row.get("n_terminal_aa", "")),
            "C端氨基酸": 清洗字符串(row.get("c_terminal_aa", "")),
            "精确IC50_uM": 格式化浮点(ic50_uM),
            "log10_IC50_uM": 格式化浮点(log_ic50),
            "来源原文": 清洗字符串(row.get("source_raw", "")),
            "Assay原文": 清洗字符串(row.get("assay_raw", "")),
            "方法原文": 清洗字符串(row.get("method_raw", "")),
            "模型原文": 清洗字符串(row.get("mice_raw", "")),
            "原始IC50文本": 清洗字符串(row.get("ic50_raw", "")),
            "是否存在反向序列": "是" if 判定为真(row.get("reverse_exists_flag", "")) else "否",
            "反向序列": 清洗字符串(row.get("reverse_sequence", "")),
            "该序列总重复记录数": 清洗字符串(row.get("sequence_dup_count", "")),
            "该序列精确IC50唯一值个数": 清洗字符串(row.get("exact_ic50_unique_count_by_sequence", "")),
            "该序列是否存在精确IC50冲突": "是" if 判定为真(row.get("exact_ic50_conflict_flag", "")) else "否",
        })

    def 排序键(x: Dict[str, Any]):
        peptide_rank = 0 if x["短肽类型"] == "二肽" else 1
        ic50 = 安全转浮点(x["精确IC50_uM"])
        ic50 = ic50 if ic50 is not None else float("inf")
        return (peptide_rank, x["序列"], ic50, x["原始记录ID"])

    result.sort(key=排序键)
    return result


# =========================================================
# 六、构建表2：高置信二三肽序列聚合表
# =========================================================

def 几何均值(values: List[float]) -> Optional[float]:
    """
    计算几何均值。
    只有在所有值都 > 0 时才计算，否则返回 None。
    """
    if not values:
        return None
    if any(v <= 0 for v in values):
        return None

    try:
        log_sum = sum(math.log(v) for v in values)
        return math.exp(log_sum / len(values))
    except Exception:
        return None


def 构建高置信二三肽序列聚合表(rows: List[Dict[str, str]]) -> List[Dict[str, Any]]:
    """
    将高置信短肽记录按 sequence_clean 聚合成“序列级”表。

    这张表可以回答：
    - 同一序列有多少条高置信记录？
    - 同一序列有多少个不同精确 IC50？
    - 这个序列整体上稳不稳？
    - 来源和 assay 是否很多样？
    """
    grouped = defaultdict(list)

    for row in rows:
        if 是否高置信短肽记录(row):
            seq = 清洗字符串(row.get("sequence_clean", ""))
            if seq != "":
                grouped[seq].append(row)

    result = []

    for seq, subrows in grouped.items():
        values = []
        ids = []
        sources = []
        assays = []
        methods = []
        models = []

        for r in subrows:
            v = 安全转浮点(r.get("ic50_uM", ""))
            if v is not None:
                values.append(v)

            ids.append(清洗字符串(r.get("id", "")))
            sources.append(清洗字符串(r.get("source_raw", "")))
            assays.append(清洗字符串(r.get("assay_raw", "")))
            methods.append(清洗字符串(r.get("method_raw", "")))
            models.append(清洗字符串(r.get("mice_raw", "")))

        unique_values = sorted(set(values))
        min_v = min(unique_values) if unique_values else None
        max_v = max(unique_values) if unique_values else None
        med_v = median(unique_values) if unique_values else None
        gmean_v = 几何均值(unique_values)
        log10_med_v = 对数10(med_v)

        等级, 说明, ratio = 稳定性等级(unique_values)

        first = subrows[0]
        短肽类别 = 短肽类型(first)

        result.append({
            "序列": seq,
            "短肽类型": 短肽类别,
            "肽长": 清洗字符串(first.get("len_clean", "")),
            "N端氨基酸": 清洗字符串(first.get("n_terminal_aa", "")),
            "C端氨基酸": 清洗字符串(first.get("c_terminal_aa", "")),
            "记录数": len(subrows),
            "精确IC50唯一值个数": len(unique_values),
            "精确IC50最小值_uM": 格式化浮点(min_v),
            "精确IC50最大值_uM": 格式化浮点(max_v),
            "精确IC50中位数_uM": 格式化浮点(med_v),
            "精确IC50几何均值_uM": 格式化浮点(gmean_v),
            "log10_中位数IC50_uM": 格式化浮点(log10_med_v),
            "最大最小倍数": 格式化浮点(ratio),
            "稳定性等级": 等级,
            "稳定性判定说明": 说明,
            "来源种类数": len(set([s for s in sources if s != ""])),
            "Assay种类数": len(set([a for a in assays if a != ""])),
            "方法种类数": len(set([m for m in methods if m != ""])),
            "模型种类数": len(set([m for m in models if m != ""])),
            "来源列表": 去重并排序(sources),
            "Assay列表": 去重并排序(assays),
            "方法列表": 去重并排序(methods),
            "模型列表": 去重并排序(models),
            "记录ID列表": 去重并排序(ids),
            "是否存在反向序列": "是" if 判定为真(first.get("reverse_exists_flag", "")) else "否",
            "反向序列": 清洗字符串(first.get("reverse_sequence", "")),
        })

    稳定性排序映射 = {
        "稳定序列": 0,
        "条件敏感序列": 1,
        "高冲突序列": 2,
        "无法判定": 3,
    }

    def 排序键(x: Dict[str, Any]):
        peptide_rank = 0 if x["短肽类型"] == "二肽" else 1
        stable_rank = 稳定性排序映射.get(x["稳定性等级"], 99)
        med = 安全转浮点(x["精确IC50中位数_uM"])
        med = med if med is not None else float("inf")
        return (peptide_rank, stable_rank, med, x["序列"])

    result.sort(key=排序键)
    return result


# =========================================================
# 七、构建表3：高置信二三肽稳定序列表
# =========================================================

def 构建高置信二三肽稳定序列表(聚合表: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """
    在“高置信二三肽序列聚合表”的基础上，进一步筛出稳定序列。

    当前纳入规则：
        只保留 稳定性等级 = 稳定序列 的行。

    这张表通常是你后续：
    - N/C 端富集
    - 强活性 vs 全体比较
    - 顺序翻转分析
    - 第一版预测模型
    最适合优先使用的核心样本表。
    """
    result = []

    for row in 聚合表:
        if 清洗字符串(row.get("稳定性等级", "")) != "稳定序列":
            continue

        med = 安全转浮点(row.get("精确IC50中位数_uM", ""))
        log_med = 对数10(med)

        if med is None:
            活性层级 = "无法判定"
        elif med <= 50:
            活性层级 = "强活性"
        elif med <= 200:
            活性层级 = "中等活性"
        else:
            活性层级 = "较弱活性"

        result.append({
            "序列": row.get("序列", ""),
            "短肽类型": row.get("短肽类型", ""),
            "肽长": row.get("肽长", ""),
            "N端氨基酸": row.get("N端氨基酸", ""),
            "C端氨基酸": row.get("C端氨基酸", ""),
            "精确IC50中位数_uM": row.get("精确IC50中位数_uM", ""),
            "log10_中位数IC50_uM": 格式化浮点(log_med),
            "精确IC50几何均值_uM": row.get("精确IC50几何均值_uM", ""),
            "精确IC50唯一值个数": row.get("精确IC50唯一值个数", ""),
            "记录数": row.get("记录数", ""),
            "最大最小倍数": row.get("最大最小倍数", ""),
            "活性分层建议": 活性层级,
            "来源种类数": row.get("来源种类数", ""),
            "Assay种类数": row.get("Assay种类数", ""),
            "来源列表": row.get("来源列表", ""),
            "Assay列表": row.get("Assay列表", ""),
            "是否存在反向序列": row.get("是否存在反向序列", ""),
            "反向序列": row.get("反向序列", ""),
            "稳定性判定说明": row.get("稳定性判定说明", ""),
        })

    活性排序映射 = {
        "强活性": 0,
        "中等活性": 1,
        "较弱活性": 2,
        "无法判定": 3,
    }

    def 排序键(x: Dict[str, Any]):
        peptide_rank = 0 if x["短肽类型"] == "二肽" else 1
        activity_rank = 活性排序映射.get(x["活性分层建议"], 99)
        med = 安全转浮点(x["精确IC50中位数_uM"])
        med = med if med is not None else float("inf")
        return (peptide_rank, activity_rank, med, x["序列"])

    result.sort(key=排序键)
    return result


# =========================================================
# 八、主程序
# =========================================================

def main():
    parser = argparse.ArgumentParser(
        description="从 AHTPDB 总表构建高置信短肽工作集（纯标准库版）。"
    )
    parser.add_argument(
        "--project-root",
        type=str,
        default="",
        help="PepDB 项目根目录，例如 E:/MYS/PepDB。仅在你关闭手动路径优先时有用。"
    )
    args = parser.parse_args()

    script_path = Path(__file__)

    project_root, input_csv, output_dir = 解析最终输入输出路径(args.project_root, script_path)

    print("=" * 90)
    print("构建高置信短肽工作集脚本（纯标准库版）")
    print("=" * 90)
    print(f"项目根目录：{project_root}")
    print(f"输入文件：{input_csv}")
    print(f"输出目录：{output_dir}")
    print(f"是否优先使用文件前部手动路径：{'是' if PRIORITY_USE_MANUAL_PATHS else '否'}")
    print("-" * 90)

    if not input_csv.exists():
        raise FileNotFoundError(
            f"未找到输入文件：{input_csv}\n"
            f"请先确认 ahtpdb_master_clean.csv 已存在，或修改脚本前部的 MANUAL_INPUT_CSV。"
        )

    rows = 读取CSV(input_csv)
    print(f"成功读取记录数：{len(rows)}")

    检查必要列(
        rows,
        [
            "id",
            "table_from",
            "sequence_clean",
            "sequence_valid",
            "len_clean",
            "len_match",
            "ic50_type",
            "ic50_exact_flag",
            "ic50_uM",
            "source_raw",
            "assay_raw",
            "method_raw",
            "mice_raw",
            "ic50_raw",
            "reverse_sequence",
            "reverse_exists_flag",
            "n_terminal_aa",
            "c_terminal_aa",
            "sequence_dup_count",
            "exact_ic50_unique_count_by_sequence",
            "exact_ic50_conflict_flag",
            "is_dipeptide",
            "is_tripeptide",
        ],
    )

    print("检测到的输入列名如下：")
    print(list(rows[0].keys()))
    print("-" * 90)

    表1_记录表 = 构建高置信二三肽记录表(rows)
    表2_聚合表 = 构建高置信二三肽序列聚合表(rows)
    表3_稳定表 = 构建高置信二三肽稳定序列表(表2_聚合表)

    写出CSV(
        output_dir / "高置信二三肽记录表.csv",
        [
            "原始记录ID", "来源原表", "序列", "短肽类型", "肽长",
            "N端氨基酸", "C端氨基酸",
            "精确IC50_uM", "log10_IC50_uM",
            "来源原文", "Assay原文", "方法原文", "模型原文",
            "原始IC50文本",
            "是否存在反向序列", "反向序列",
            "该序列总重复记录数", "该序列精确IC50唯一值个数", "该序列是否存在精确IC50冲突",
        ],
        表1_记录表,
    )

    写出CSV(
        output_dir / "高置信二三肽序列聚合表.csv",
        [
            "序列", "短肽类型", "肽长", "N端氨基酸", "C端氨基酸",
            "记录数", "精确IC50唯一值个数",
            "精确IC50最小值_uM", "精确IC50最大值_uM", "精确IC50中位数_uM",
            "精确IC50几何均值_uM", "log10_中位数IC50_uM", "最大最小倍数",
            "稳定性等级", "稳定性判定说明",
            "来源种类数", "Assay种类数", "方法种类数", "模型种类数",
            "来源列表", "Assay列表", "方法列表", "模型列表", "记录ID列表",
            "是否存在反向序列", "反向序列",
        ],
        表2_聚合表,
    )

    写出CSV(
        output_dir / "高置信二三肽稳定序列表.csv",
        [
            "序列", "短肽类型", "肽长", "N端氨基酸", "C端氨基酸",
            "精确IC50中位数_uM", "log10_中位数IC50_uM", "精确IC50几何均值_uM",
            "精确IC50唯一值个数", "记录数", "最大最小倍数",
            "活性分层建议", "来源种类数", "Assay种类数",
            "来源列表", "Assay列表",
            "是否存在反向序列", "反向序列",
            "稳定性判定说明",
        ],
        表3_稳定表,
    )

    二肽记录数 = sum(1 for r in 表1_记录表 if 清洗字符串(r.get("短肽类型", "")) == "二肽")
    三肽记录数 = sum(1 for r in 表1_记录表 if 清洗字符串(r.get("短肽类型", "")) == "三肽")

    二肽稳定数 = sum(1 for r in 表3_稳定表 if 清洗字符串(r.get("短肽类型", "")) == "二肽")
    三肽稳定数 = sum(1 for r in 表3_稳定表 if 清洗字符串(r.get("短肽类型", "")) == "三肽")

    print("结果表已输出：")
    print(f"1. {output_dir / '高置信二三肽记录表.csv'}")
    print(f"2. {output_dir / '高置信二三肽序列聚合表.csv'}")
    print(f"3. {output_dir / '高置信二三肽稳定序列表.csv'}")
    print("-" * 90)
    print("摘要：")
    print(f"高置信二三肽记录数：{len(表1_记录表)}")
    print(f"  其中二肽记录数：{二肽记录数}")
    print(f"  其中三肽记录数：{三肽记录数}")
    print(f"高置信二三肽唯一序列数：{len(表2_聚合表)}")
    print(f"高置信二三肽稳定序列数：{len(表3_稳定表)}")
    print(f"  其中稳定二肽序列数：{二肽稳定数}")
    print(f"  其中稳定三肽序列数：{三肽稳定数}")
    print("=" * 90)
    print("工作集构建完成。")


if __name__ == "__main__":
    main()
