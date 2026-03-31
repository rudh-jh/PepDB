#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
脚本名称：
    analyze_ahtpdb_label_quality_stdlib.py

脚本用途：
    读取已经标准化完成的 AHTPDB 总表（ahtpdb_master_clean.csv），
    对其中的“标签质量”做一轮基础分析，并输出若干张中文结果表，
    方便后续定义高置信训练集、判断重复序列冲突、开展构效分析和预测任务。

这次改动说明：
    你希望把“输入文件/输出目录”放到脚本前部单独封装，
    方便以后直接改路径，而不需要在主逻辑里翻找。

    因此本脚本新增了“路径配置区”：
        - MANUAL_PROJECT_ROOT
        - MANUAL_INPUT_CSV
        - MANUAL_OUTPUT_DIR
        - PRIORITY_USE_MANUAL_PATHS

    你以后最常改的通常只有：
        1）MANUAL_PROJECT_ROOT
        2）MANUAL_INPUT_CSV
        3）MANUAL_OUTPUT_DIR

    其余分析逻辑不变。

脚本特点：
    - 纯 Python 标准库实现，不依赖 pandas / numpy / pyarrow
    - 输出 CSV 使用 utf-8-sig，方便 Windows 下直接用 Excel 打开
    - 输出表名、列名均使用中文
    - 只做“标签质量分析”，不做建模、不画图

默认输出结果：
    1. IC50类型分布表.csv
    2. IC50单位分布表.csv
    3. 标签分级统计表.csv
    4. 重复序列次数分布表.csv
    5. 精确IC50冲突明细表.csv
    6. 冲突与实验来源关系表.csv
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Any, Iterable, Tuple


# =========================================================
# 一、路径配置区（你后续最常改这里）
# =========================================================
# 使用说明：
# 1）如果 PRIORITY_USE_MANUAL_PATHS = True，脚本会优先使用下面这几项手动路径。
# 2）如果你把它改成 False，脚本就会回到“自动推断项目根目录 + 由项目根目录计算输入/输出路径”的逻辑。
# 3）命令行里如果显式传入 --project-root，也仍然可以覆盖自动推断逻辑。
# 4）推荐你长期维护时，直接改 MANUAL_INPUT_CSV 和 MANUAL_OUTPUT_DIR，就最省事。

PRIORITY_USE_MANUAL_PATHS = True

# PepDB 项目根目录（可选）
MANUAL_PROJECT_ROOT = r"E:\MYS\PepDB"

# 输入文件：标准化后的 AHTPDB 总表
MANUAL_INPUT_CSV = r"E:\MYS\PepDB\DB\standardized\ace\ahtpdb\ahtpdb_master_clean_um_expanded.csv"

# 输出目录：标签质量分析结果目录
MANUAL_OUTPUT_DIR = r"E:\MYS\PepDB\DB\analysis\ace\ahtpdb\标签质量分析"


# =========================================================
# 二、通用小工具函数
# =========================================================

def 清洗字符串(value: Any) -> str:
    """
    将任意值安全地转成字符串，并做基本清洗。

    处理逻辑：
    - None -> 空字符串
    - 去掉首尾空白字符
    - 保持原始文本主体不变
    """
    if value is None:
        return ""
    return str(value).strip()


def 判定为真(value: Any) -> bool:
    """
    将表中的各种真假写法统一判定为布尔值。

    支持：
    - 1 / 0
    - True / False
    - yes / no
    - y / n
    - t / f
    """
    text = 清洗字符串(value).lower()
    return text in {"1", "true", "t", "yes", "y"}


def 安全转浮点(value: Any):
    """
    尝试将字段转成浮点数。
    如果无法转换，返回 None。
    """
    text = 清洗字符串(value)
    if text == "":
        return None
    try:
        return float(text)
    except Exception:
        return None


def 去重并排序字符串(values: Iterable[str], 分隔符: str = " | ") -> str:
    """
    对一组字符串去重、排序后拼接。
    适用于 assay/source/id 等列表字段。
    """
    cleaned = []
    for v in values:
        s = 清洗字符串(v)
        if s != "":
            cleaned.append(s)
    unique_sorted = sorted(set(cleaned))
    return 分隔符.join(unique_sorted)


def 百分比字符串(分子: int, 分母: int) -> str:
    """
    将分子/分母转成百分比字符串，保留两位小数。
    """
    if 分母 == 0:
        return "0.00%"
    return f"{(分子 / 分母) * 100:.2f}%"


# =========================================================
# 三、CSV 读写函数
# =========================================================

def 读取CSV(filepath: Path) -> List[Dict[str, str]]:
    """
    读取 CSV 文件，返回“字典列表”。

    这个增强版函数专门用于解决 Windows 环境下常见的 CSV 编码问题。
    它会依次尝试：
        - utf-8-sig
        - utf-8
        - gb18030
        - gbk
        - cp936
        - latin1（兜底）

    同时还会尝试自动识别分隔符（逗号 / 制表符 / 分号）。
    """
    编码候选列表 = [
        "utf-8-sig",
        "utf-8",
        "gb18030",
        "gbk",
        "cp936",
        "latin1",
    ]

    最后异常信息 = []

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

                临时结果 = []
                for row in reader:
                    cleaned_row = {}
                    for k, v in row.items():
                        key = 清洗字符串(k)
                        value = 清洗字符串(v)
                        cleaned_row[key] = value
                    临时结果.append(cleaned_row)

                print(f"成功使用编码读取文件：{编码}")
                return 临时结果

        except UnicodeDecodeError as e:
            最后异常信息.append(f"{编码}: UnicodeDecodeError -> {e}")
            continue
        except Exception as e:
            最后异常信息.append(f"{编码}: {type(e).__name__} -> {e}")
            continue

    raise RuntimeError(
        "无法成功读取输入 CSV 文件，已尝试多种编码但均失败。\n"
        "请检查文件是否损坏，或是否是 Excel 另存为时使用了特殊编码。\n"
        "已尝试的编码及报错如下：\n"
        + "\n".join(最后异常信息)
    )


def 写出CSV(filepath: Path, fieldnames: List[str], rows: List[Dict[str, Any]]) -> None:
    """
    将字典列表写出为 CSV。

    说明：
    - 自动创建父目录
    - 使用 utf-8-sig，便于 Excel 正确显示中文
    - extrasaction='ignore'：如果 rows 里有多余字段，自动忽略
    """
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with filepath.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def 检查必要列是否存在(rows: List[Dict[str, str]], 必要列列表: List[str]) -> None:
    """
    检查输入表是否包含后续分析所需的关键列。
    如果缺列，则尽早报错，避免后面在别处才崩。
    """
    if not rows:
        raise ValueError("输入表为空，无法检查列名。")

    实际列名 = set(rows[0].keys())
    缺失列 = [col for col in 必要列列表 if col not in 实际列名]

    if 缺失列:
        raise ValueError(
            "输入表缺少以下必要列：\n"
            + ", ".join(缺失列)
            + "\n\n实际读取到的列名为：\n"
            + ", ".join(rows[0].keys())
        )


# =========================================================
# 四、路径解析函数
# =========================================================

def 自动推断项目根目录(脚本路径: Path) -> Path:
    """
    根据脚本所在位置，自动推断 PepDB 项目根目录。
    约定：如果脚本放在 PepDB/SCR/scripts/ 目录下，则项目根目录为 PepDB/。
    """
    return 脚本路径.resolve().parent.parent.parent


def 获取输入输出路径(project_root: Path) -> Tuple[Path, Path]:
    """
    根据项目根目录，计算默认输入文件和默认输出目录。

    注意：
    这个函数只有在“不开启手动路径优先”时才会真正发挥作用。
    因为你这次希望路径能在文件前部直接改，所以脚本优先使用路径配置区。
    """
    input_csv = project_root / "DB" / "standardized" / "ace" / "ahtpdb" / "ahtpdb_master_clean.csv"
    output_dir = project_root / "DB" / "analysis" / "ace" / "ahtpdb" / "标签质量分析"
    return input_csv, output_dir


def 解析最终输入输出路径(args_project_root: str, script_path: Path) -> Tuple[Path, Path, Path]:
    """
    统一解析最终实际使用的：
        - 项目根目录
        - 输入文件路径
        - 输出目录路径

    优先级规则：
    1）如果 PRIORITY_USE_MANUAL_PATHS = True，则优先采用文件前部手动配置的路径
    2）否则，如果命令行传了 --project-root，则按该根目录推导默认输入输出路径
    3）否则，自动按脚本位置推断项目根目录，再推导默认输入输出路径

    这么封装的目的：
    让你以后修改路径时，只需要改前面的三个变量；
    主程序和统计逻辑都不用动。
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


# =========================================================
# 五、标签分级判定
# =========================================================

def 判定标签等级(row: Dict[str, str]) -> str:
    """
    将每条记录划分为 A/B/C 三个标签等级。

    A级高置信：
        - sequence_valid = 1
        - len_match = 1
        - ic50_exact_flag = 1
        - ic50_type = exact
        - ic50_uM 非空且能转成数值

    B级中置信：
        - 不是 A 级
        - 但 ic50_type 属于以下之一：
            threshold / range / mass_unit / mass_conc / mass_concentration

    C级低置信：
        - 其余情况
        - 包括 suspicious_unit / dirty / missing / percent 等
    """
    sequence_valid = 判定为真(row.get("sequence_valid", ""))
    len_match = 判定为真(row.get("len_match", ""))
    ic50_exact_flag = 判定为真(row.get("ic50_exact_flag", ""))
    ic50_type = 清洗字符串(row.get("ic50_type", "")).lower()
    ic50_uM = 安全转浮点(row.get("ic50_uM", ""))

    if sequence_valid and len_match and ic50_exact_flag and ic50_type == "exact" and ic50_uM is not None:
        return "A级高置信"

    if ic50_type in {"threshold", "range", "mass_unit", "mass_conc", "mass_concentration"}:
        return "B级中置信"

    return "C级低置信"


def 标签等级判定说明(等级: str) -> str:
    """
    为每个标签等级返回中文说明，便于结果表直接阅读。
    """
    if 等级 == "A级高置信":
        return "有效序列 + 长度一致 + 精确IC50点值 + 已统一到uM"
    if 等级 == "B级中置信":
        return "非精确点值，但仍具参考意义，如阈值/范围值/质量浓度"
    return "可疑单位、脏值、缺失值或其他低置信标签"


# =========================================================
# 六、各结果表构建函数
# =========================================================

def 生成IC50类型分布表(rows: List[Dict[str, str]]) -> List[Dict[str, Any]]:
    counter = Counter()
    for row in rows:
        ic50_type = 清洗字符串(row.get("ic50_type", ""))
        if ic50_type == "":
            ic50_type = "空白/未标注"
        counter[ic50_type] += 1

    total = sum(counter.values())
    result = []
    for ic50_type, count in counter.most_common():
        result.append({
            "IC50类型": ic50_type,
            "样本数": count,
            "占比": 百分比字符串(count, total),
        })
    return result


def 生成IC50单位分布表(rows: List[Dict[str, str]]) -> List[Dict[str, Any]]:
    counter = Counter()
    for row in rows:
        unit = 清洗字符串(row.get("ic50_unit_raw", ""))
        if unit == "":
            unit = "空白/未标注"
        counter[unit] += 1

    total = sum(counter.values())
    result = []
    for unit, count in counter.most_common():
        result.append({
            "IC50原始单位": unit,
            "样本数": count,
            "占比": 百分比字符串(count, total),
        })
    return result


def 生成标签分级统计表(rows: List[Dict[str, str]]) -> List[Dict[str, Any]]:
    bucket = {
        "A级高置信": [],
        "B级中置信": [],
        "C级低置信": [],
    }

    for row in rows:
        level = 判定标签等级(row)
        bucket[level].append(row)

    result = []
    for level in ["A级高置信", "B级中置信", "C级低置信"]:
        sub = bucket[level]
        total_count = len(sub)
        valid_seq_count = sum(1 for r in sub if 判定为真(r.get("sequence_valid", "")))
        di_count = sum(1 for r in sub if 判定为真(r.get("is_dipeptide", "")))
        tri_count = sum(1 for r in sub if 判定为真(r.get("is_tripeptide", "")))

        result.append({
            "标签等级": level,
            "判定规则": 标签等级判定说明(level),
            "样本数": total_count,
            "占全部样本比例": 百分比字符串(total_count, len(rows)),
            "有效序列样本数": valid_seq_count,
            "二肽样本数": di_count,
            "三肽样本数": tri_count,
        })

    return result


def 生成重复序列次数分布表(rows: List[Dict[str, str]]) -> List[Dict[str, Any]]:
    seq_to_dup = {}

    for row in rows:
        seq = 清洗字符串(row.get("sequence_clean", ""))
        if seq == "":
            continue

        dup_count_raw = 清洗字符串(row.get("sequence_dup_count", ""))
        try:
            dup_count = int(float(dup_count_raw)) if dup_count_raw != "" else 1
        except Exception:
            dup_count = 1

        if seq not in seq_to_dup:
            seq_to_dup[seq] = dup_count
        else:
            seq_to_dup[seq] = max(seq_to_dup[seq], dup_count)

    dup_counter = Counter(seq_to_dup.values())

    result = []
    for dup_count in sorted(dup_counter.keys()):
        unique_seq_num = dup_counter[dup_count]
        total_record_num = dup_count * unique_seq_num
        result.append({
            "重复次数": dup_count,
            "对应唯一序列数": unique_seq_num,
            "对应记录总数": total_record_num,
        })
    return result


def 生成精确IC50冲突明细表(rows: List[Dict[str, str]]) -> List[Dict[str, Any]]:
    seq_rows = defaultdict(list)
    for row in rows:
        seq = 清洗字符串(row.get("sequence_clean", ""))
        if seq != "":
            seq_rows[seq].append(row)

    result = []

    for seq, subrows in seq_rows.items():
        exact_values = []
        ids = []
        sources = []
        assays = []

        for r in subrows:
            ids.append(清洗字符串(r.get("id", "")))
            sources.append(清洗字符串(r.get("source_raw", "")))
            assays.append(清洗字符串(r.get("assay_raw", "")))

            exact_flag = 判定为真(r.get("ic50_exact_flag", ""))
            ic50_uM = 安全转浮点(r.get("ic50_uM", ""))

            if exact_flag and ic50_uM is not None:
                exact_values.append(f"{ic50_uM:.6g}")

        unique_exact_values = sorted(set(exact_values))
        explicit_conflict = any(判定为真(r.get("exact_ic50_conflict_flag", "")) for r in subrows)
        implicit_conflict = len(unique_exact_values) >= 2

        if explicit_conflict or implicit_conflict:
            dup_count_candidates = []
            for r in subrows:
                raw_dup = 清洗字符串(r.get("sequence_dup_count", ""))
                if raw_dup != "":
                    try:
                        dup_count_candidates.append(int(float(raw_dup)))
                    except Exception:
                        pass

            if dup_count_candidates:
                dup_count = max(dup_count_candidates)
            else:
                dup_count = len(subrows)

            result.append({
                "序列": seq,
                "重复记录数": dup_count,
                "精确IC50唯一值个数": len(unique_exact_values),
                "精确IC50值列表_uM": " | ".join(unique_exact_values),
                "来源列表": 去重并排序字符串(sources),
                "Assay列表": 去重并排序字符串(assays),
                "记录ID列表": 去重并排序字符串(ids),
            })

    result.sort(key=lambda x: (-int(x["精确IC50唯一值个数"]), -int(x["重复记录数"]), x["序列"]))
    return result


def 生成冲突与实验来源关系表(rows: List[Dict[str, str]]) -> List[Dict[str, Any]]:
    seq_rows = defaultdict(list)
    for row in rows:
        seq = 清洗字符串(row.get("sequence_clean", ""))
        if seq != "":
            seq_rows[seq].append(row)

    result = []

    for seq, subrows in seq_rows.items():
        exact_values = []
        sources = []
        assays = []

        for r in subrows:
            exact_flag = 判定为真(r.get("ic50_exact_flag", ""))
            ic50_uM = 安全转浮点(r.get("ic50_uM", ""))

            if exact_flag and ic50_uM is not None:
                exact_values.append(f"{ic50_uM:.6g}")

            source = 清洗字符串(r.get("source_raw", ""))
            assay = 清洗字符串(r.get("assay_raw", ""))

            if source != "":
                sources.append(source)
            if assay != "":
                assays.append(assay)

        unique_exact_values = sorted(set(exact_values))
        unique_sources = sorted(set(sources))
        unique_assays = sorted(set(assays))

        if len(unique_exact_values) >= 2:
            result.append({
                "序列": seq,
                "冲突精确IC50值列表_uM": " | ".join(unique_exact_values),
                "来源种类数": len(unique_sources),
                "来源列表": " | ".join(unique_sources),
                "Assay种类数": len(unique_assays),
                "Assay列表": " | ".join(unique_assays),
                "是否跨来源": "是" if len(unique_sources) >= 2 else "否",
                "是否跨Assay": "是" if len(unique_assays) >= 2 else "否",
            })

    result.sort(
        key=lambda x: (
            0 if x["是否跨Assay"] == "是" else 1,
            0 if x["是否跨来源"] == "是" else 1,
            -int(x["Assay种类数"]),
            -int(x["来源种类数"]),
            x["序列"],
        )
    )
    return result


# =========================================================
# 七、主程序
# =========================================================

def main():
    parser = argparse.ArgumentParser(
        description="读取 AHTPDB 总表，输出标签质量分析结果表（纯标准库版）。"
    )
    parser.add_argument(
        "--project-root",
        type=str,
        default="",
        help="PepDB 项目根目录，例如 E:/MYS/PepDB。仅在你关闭手动路径优先时有用。"
    )
    args = parser.parse_args()

    script_path = Path(__file__)

    # 统一解析最终实际使用的项目根目录、输入文件和输出目录
    project_root, input_csv, output_dir = 解析最终输入输出路径(args.project_root, script_path)

    print("=" * 80)
    print("AHTPDB 标签质量分析脚本（纯标准库版）")
    print("=" * 80)
    print(f"项目根目录：{project_root}")
    print(f"输入文件：{input_csv}")
    print(f"输出目录：{output_dir}")
    print(f"是否优先使用文件前部手动路径：{'是' if PRIORITY_USE_MANUAL_PATHS else '否'}")
    print("-" * 80)

    if not input_csv.exists():
        raise FileNotFoundError(
            f"未找到输入文件：{input_csv}\n"
            f"请确认你已经先生成 ahtpdb_master_clean.csv，"
            f"或者修改脚本前部的 MANUAL_INPUT_CSV。"
        )

    rows = 读取CSV(input_csv)
    print(f"成功读取记录数：{len(rows)}")

    if rows:
        print("检测到的列名如下：")
        print(list(rows[0].keys()))

    if len(rows) == 0:
        raise ValueError("输入 CSV 没有数据行，无法继续分析。")

    检查必要列是否存在(
        rows,
        [
            "sequence_clean",
            "sequence_valid",
            "len_match",
            "ic50_type",
            "ic50_exact_flag",
            "ic50_uM",
            "sequence_dup_count",
            "exact_ic50_conflict_flag",
            "source_raw",
            "assay_raw",
            "is_dipeptide",
            "is_tripeptide",
        ],
    )

    表1 = 生成IC50类型分布表(rows)
    表2 = 生成IC50单位分布表(rows)
    表3 = 生成标签分级统计表(rows)
    表4 = 生成重复序列次数分布表(rows)
    表5 = 生成精确IC50冲突明细表(rows)
    表6 = 生成冲突与实验来源关系表(rows)

    写出CSV(
        output_dir / "IC50类型分布表.csv",
        ["IC50类型", "样本数", "占比"],
        表1,
    )

    写出CSV(
        output_dir / "IC50单位分布表.csv",
        ["IC50原始单位", "样本数", "占比"],
        表2,
    )

    写出CSV(
        output_dir / "标签分级统计表.csv",
        ["标签等级", "判定规则", "样本数", "占全部样本比例", "有效序列样本数", "二肽样本数", "三肽样本数"],
        表3,
    )

    写出CSV(
        output_dir / "重复序列次数分布表.csv",
        ["重复次数", "对应唯一序列数", "对应记录总数"],
        表4,
    )

    写出CSV(
        output_dir / "精确IC50冲突明细表.csv",
        ["序列", "重复记录数", "精确IC50唯一值个数", "精确IC50值列表_uM", "来源列表", "Assay列表", "记录ID列表"],
        表5,
    )

    写出CSV(
        output_dir / "冲突与实验来源关系表.csv",
        ["序列", "冲突精确IC50值列表_uM", "来源种类数", "来源列表", "Assay种类数", "Assay列表", "是否跨来源", "是否跨Assay"],
        表6,
    )

    print("-" * 80)
    print("结果表已输出：")
    print(f"1. {output_dir / 'IC50类型分布表.csv'}")
    print(f"2. {output_dir / 'IC50单位分布表.csv'}")
    print(f"3. {output_dir / '标签分级统计表.csv'}")
    print(f"4. {output_dir / '重复序列次数分布表.csv'}")
    print(f"5. {output_dir / '精确IC50冲突明细表.csv'}")
    print(f"6. {output_dir / '冲突与实验来源关系表.csv'}")
    print("-" * 80)

    print("简要摘要：")
    print(f"总记录数：{len(rows)}")
    print(f"IC50类型种类数：{len(表1)}")
    print(f"IC50单位种类数：{len(表2)}")
    print(f"存在精确IC50冲突的序列数：{len(表5)}")
    print("=" * 80)
    print("分析完成。")


if __name__ == "__main__":
    main()
