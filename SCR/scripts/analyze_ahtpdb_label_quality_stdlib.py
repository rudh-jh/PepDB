#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
脚本名称：
    analyze_ahtpdb_label_quality_stdlib.py

脚本用途：
    读取已经标准化完成的 AHTPDB 总表（ahtpdb_master_clean.csv），
    对其中的“标签质量”做一轮基础分析，并输出若干张中文结果表，
    方便后续你去定义高置信训练集、判断重复序列冲突、开展构效分析和预测任务。

为什么写这个脚本：
    你当前最需要的，不是立刻上模型，而是先搞清楚：
    1）手上到底有多少“可用的精确 IC50 标签”
    2）同一序列是否有重复、是否存在多个冲突的精确 IC50
    3）冲突和 assay/source 是否有关
    4）哪些样本可作为高置信标签，哪些只能作为中低置信证据

脚本特点：
    - 纯 Python 标准库实现，不依赖 pandas / numpy / pyarrow
    - 输出 CSV 使用 utf-8-sig，方便 Windows 下直接用 Excel 打开
    - 输出表名、列名均使用中文
    - 只做“标签质量分析”，不做建模、不画图

默认输入文件：
    DB/standardized/ace/ahtpdb/ahtpdb_master_clean.csv

默认输出目录：
    DB/analysis/ace/ahtpdb/标签质量分析/

输出结果表：
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


# =========================
# 一、通用小工具函数
# =========================

def 清洗字符串(value: Any) -> str:
    """
    将任意值安全地转成字符串，并做基本清洗。

    处理逻辑：
    - None -> 空字符串
    - 去掉首尾空白字符
    - 保持原始文本主体不变

    注意：
    这里不做大小写强转，不做内容修改，只做最轻量清洗。
    """
    if value is None:
        return ""
    return str(value).strip()


def 判定为真(value: Any) -> bool:
    """
    将表中的各种“真假/0/1/yes/no/true/false”风格文本，统一判定为布尔值。

    这是因为 CSV 中很多列可能是：
    - 1 / 0
    - True / False
    - yes / no
    - y / n
    - 空字符串

    这里统一做兼容，避免后面判断出错。
    """
    text = 清洗字符串(value).lower()
    return text in {"1", "true", "t", "yes", "y"}


def 安全转浮点(value: Any):
    """
    尝试将字段转成浮点数。
    如果无法转换，返回 None。

    说明：
    后续在处理 ic50_uM、pi 等字段时会用到。
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

    适用于：
    - assay 列表
    - source 列表
    - id 列表
    - 精确 IC50 值列表

    注意：
    这里先去空、再去重、再排序。
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
    将“分子/分母”转成百分比字符串，保留两位小数。

    例如：
        25 / 100 -> '25.00%'
    """
    if 分母 == 0:
        return "0.00%"
    return f"{(分子 / 分母) * 100:.2f}%"


def 读取CSV(filepath: Path) -> List[Dict[str, str]]:
    """
    读取 CSV 文件，返回“字典列表”。

    这个增强版函数专门用来解决 Windows 环境下常见的 CSV 编码问题。
    因为你的 ahtpdb_master_clean.csv 可能不是 utf-8-sig，而可能是：
        - utf-8-sig
        - utf-8
        - gb18030
        - gbk
        - cp936

    所以这里会按顺序逐个尝试这些编码，哪个能成功读，就用哪个。

    额外增强：
    1）自动尝试多种常见编码
    2）尝试自动识别分隔符（逗号 / 制表符 / 分号）
    3）如果前几种编码都失败，最后用 latin1 兜底读取，避免脚本直接崩掉
       但如果走到 latin1，说明文件编码很可能不规范，后续要特别留意中文乱码
    """

    rows: List[Dict[str, str]] = []

    # 按常见程度排列编码尝试顺序
    # 你的场景里，最有可能是 utf-8-sig / utf-8 / gb18030 / gbk / cp936
    编码候选列表 = [
        "utf-8-sig",
        "utf-8",
        "gb18030",
        "gbk",
        "cp936",
        "latin1",   # 兜底，不推荐长期依赖，但可以避免直接报错
    ]

    最后异常信息 = []

    for 编码 in 编码候选列表:
        try:
            with filepath.open("r", encoding=编码, newline="") as f:
                # 先读取一小段样本，用于自动识别分隔符
                sample = f.read(4096)
                f.seek(0)

                # 尝试自动识别 CSV 分隔符
                # 你的文件大概率是逗号分隔，但也可能被保存成分号或制表符
                try:
                    dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
                except Exception:
                    # 如果识别失败，默认按逗号处理
                    dialect = csv.excel

                reader = csv.DictReader(f, dialect=dialect)

                临时结果 = []
                for row in reader:
                    # 将每行的键和值做基础清洗
                    cleaned_row = {}
                    for k, v in row.items():
                        key = 清洗字符串(k)
                        value = 清洗字符串(v)
                        cleaned_row[key] = value
                    临时结果.append(cleaned_row)

                # 如果成功走到这里，说明该编码可用
                print(f"成功使用编码读取文件：{编码}")
                return 临时结果

        except UnicodeDecodeError as e:
            最后异常信息.append(f"{编码}: UnicodeDecodeError -> {e}")
            continue
        except Exception as e:
            最后异常信息.append(f"{编码}: {type(e).__name__} -> {e}")
            continue

    # 如果所有编码都失败，抛出更详细的错误信息，方便排查
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
    - extrasaction='ignore'，即如果 rows 里有多余字段，不报错，自动忽略
    """
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with filepath.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


# =========================
# 二、路径相关函数
# =========================

def 自动推断项目根目录(脚本路径: Path) -> Path:
    """
    根据脚本所在位置，自动推断 PepDB 项目根目录。

    约定：
        如果脚本放在 PepDB/SCR/scripts/ 目录下，
        那么项目根目录就是 scripts 的上上级目录，也就是 PepDB/

    例如：
        E:/MYS/PepDB/SCR/scripts/analyze_xxx.py
        -> 项目根目录推断为 E:/MYS/PepDB

    如果用户显式传了 --project-root，则优先使用用户给的路径。
    """
    return 脚本路径.resolve().parent.parent.parent


def 获取输入输出路径(project_root: Path) -> Tuple[Path, Path]:
    """
    根据项目根目录，计算输入文件和输出目录。

    输入：
        DB/standardized/ace/ahtpdb/ahtpdb_master_clean.csv

    输出目录：
        DB/analysis/ace/ahtpdb/标签质量分析/
    """
    input_csv = project_root / "DB" / "standardized" / "ace" / "ahtpdb" / "ahtpdb_master_clean.csv"
    output_dir = project_root / "DB" / "analysis" / "ace" / "ahtpdb" / "标签质量分析"
    return input_csv, output_dir


# =========================
# 三、标签等级判定函数
# =========================

def 判定标签等级(row: Dict[str, str]) -> str:
    """
    将每条记录划分为 A/B/C 三个标签等级。

    这是整个脚本里最核心的规则之一。

    规则说明：
    -----------------------------------
    A级高置信：
        - sequence_valid = 1
        - len_match = 1
        - ic50_exact_flag = 1
        - ic50_type = exact
        - ic50_uM 非空、可转为数值

    B级中置信：
        - 不是 A 级
        - 但 ic50_type 属于以下之一：
            threshold
            range
            mass_unit
            mass_conc
            mass_concentration

    C级低置信：
        - 其余情况全部归为 C 级
        - 包括 suspicious_unit / dirty / missing / percent 等
    -----------------------------------

    注意：
    这里的分级是为了“标签质量分析”和“后续定义训练集”，
    不是生物学上的最终真值分级。
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
    为每个标签等级返回一段中文说明，便于结果表直接阅读。
    """
    if 等级 == "A级高置信":
        return "有效序列 + 长度一致 + 精确IC50点值 + 已统一到uM"
    if 等级 == "B级中置信":
        return "非精确点值，但仍具参考意义，如阈值/范围值/质量浓度"
    return "可疑单位、脏值、缺失值或其他低置信标签"


# =========================
# 四、结果表一：IC50类型分布表
# =========================

def 生成IC50类型分布表(rows: List[Dict[str, str]]) -> List[Dict[str, Any]]:
    """
    统计 ic50_type 各类型的样本数与占比。

    这个表能帮助你快速判断：
    - 有多少 exact
    - 有多少 threshold/range
    - 有多少 suspicious_unit / dirty

    这对于后续定义“高置信监督集”非常重要。
    """
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
            "占比": 百分比字符串(count, total)
        })
    return result


# =========================
# 五、结果表二：IC50单位分布表
# =========================

def 生成IC50单位分布表(rows: List[Dict[str, str]]) -> List[Dict[str, Any]]:
    """
    统计 ic50_unit_raw 各单位的样本数与占比。

    这个表的作用：
    - 帮你了解 raw 数据的单位异质性
    - 为后续决定“哪些单位可以纳入统一分析”提供依据
    """
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
            "占比": 百分比字符串(count, total)
        })
    return result


# =========================
# 六、结果表三：标签分级统计表
# =========================

def 生成标签分级统计表(rows: List[Dict[str, str]]) -> List[Dict[str, Any]]:
    """
    统计 A/B/C 三个标签等级的数量。

    除了总样本数之外，还额外统计：
    - 有效序列样本数
    - 二肽样本数
    - 三肽样本数

    这样你后面能直接判断：
    - 高置信标签里有多少二肽/三肽
    - 哪一层最适合作为训练集
    """
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


# =========================
# 七、结果表四：重复序列次数分布表
# =========================

def 生成重复序列次数分布表(rows: List[Dict[str, str]]) -> List[Dict[str, Any]]:
    """
    按 sequence_dup_count 统计重复序列分布。

    说明：
    你的总表中已经有 sequence_dup_count 这一列，
    它表示“同一 sequence_clean 在整张表中出现了多少次”。

    这里我们关心的是：
    - 重复 1 次的序列有多少条
    - 重复 2 次的序列有多少条
    - ...
    同时也统计：
    - 对应的总记录数是多少

    注意：
    这个表按“唯一序列”视角统计，而不是按记录视角。
    """
    # 用 sequence_clean 去重，只保留每个序列的一个代表
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

        # 同一序列理论上应该对应同一个 sequence_dup_count
        # 若表中出现不一致，这里以最大值为准，尽量避免低估重复程度
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
            "对应记录总数": total_record_num
        })
    return result


# =========================
# 八、结果表五：精确IC50冲突明细表
# =========================

def 生成精确IC50冲突明细表(rows: List[Dict[str, str]]) -> List[Dict[str, Any]]:
    """
    生成“同一序列存在多个不同精确 IC50” 的冲突明细表。

    这是你后续定义训练集时最关键的结果表之一。

    判定逻辑：
    - sequence_clean 非空
    - exact_ic50_conflict_flag = 1
      或者
    - 同一序列拥有两个及以上不同的精确 ic50_uM

    输出内容包括：
    - 序列
    - 重复记录数
    - 精确IC50唯一值个数
    - 精确IC50值列表
    - 来源列表
    - assay 列表
    - 记录ID列表
    """
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
                # 为了避免 15.700000 这类显示问题，统一格式化
                exact_values.append(f"{ic50_uM:.6g}")

        unique_exact_values = sorted(set(exact_values))

        # 以两种方式判断冲突：
        # 1）表中显式标记 exact_ic50_conflict_flag=1
        # 2）同一序列确实有两个及以上不同精确值
        explicit_conflict = any(判定为真(r.get("exact_ic50_conflict_flag", "")) for r in subrows)
        implicit_conflict = len(unique_exact_values) >= 2

        if explicit_conflict or implicit_conflict:
            # 重复记录数尽量优先使用表中已有的 sequence_dup_count
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

    # 冲突更严重的放前面：
    # 先按“精确IC50唯一值个数”降序，再按“重复记录数”降序，再按序列排序
    result.sort(key=lambda x: (-int(x["精确IC50唯一值个数"]), -int(x["重复记录数"]), x["序列"]))
    return result


# =========================
# 九、结果表六：冲突与实验来源关系表
# =========================

def 生成冲突与实验来源关系表(rows: List[Dict[str, str]]) -> List[Dict[str, Any]]:
    """
    汇总“发生精确 IC50 冲突的序列”与 assay/source 的关系。

    你后续最想回答的问题之一是：
    “同一序列的多个 IC50 值，到底是数据库噪音，还是来源/实验条件不同导致的？”

    所以这个表按“冲突序列”汇总，输出：
    - 序列
    - 冲突精确IC50值
    - assay 数量 / assay 列表
    - source 数量 / source 列表
    - 是否跨 assay
    - 是否跨来源
    """
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

    # 优先把“跨 assay + 跨来源”的复杂情况放前面
    result.sort(
        key=lambda x: (
            0 if x["是否跨Assay"] == "是" else 1,
            0 if x["是否跨来源"] == "是" else 1,
            -int(x["Assay种类数"]),
            -int(x["来源种类数"]),
            x["序列"]
        )
    )
    return result


# =========================
# 十、主流程
# =========================

def main():
    parser = argparse.ArgumentParser(
        description="读取 AHTPDB 总表，输出标签质量分析结果表（纯标准库版）。"
    )
    parser.add_argument(
        "--project-root",
        type=str,
        default="",
        help="PepDB 项目根目录，例如 E:/MYS/PepDB。若不传，则尝试按脚本位置自动推断。"
    )
    args = parser.parse_args()

    script_path = Path(__file__)

    # 1）确定项目根目录
    if 清洗字符串(args.project_root) != "":
        project_root = Path(args.project_root).resolve()
    else:
        project_root = 自动推断项目根目录(script_path)

    # 2）计算输入文件和输出目录
    input_csv, output_dir = 获取输入输出路径(project_root)

    print("=" * 80)
    print("AHTPDB 标签质量分析脚本（纯标准库版）")
    print("=" * 80)
    print(f"项目根目录：{project_root}")
    print(f"输入文件：{input_csv}")
    print(f"输出目录：{output_dir}")
    print("-" * 80)

    # 3）检查输入文件是否存在
    if not input_csv.exists():
        raise FileNotFoundError(
            f"未找到输入文件：{input_csv}\n"
            f"请确认你已经先生成 ahtpdb_master_clean.csv，"
            f"或者使用 --project-root 显式指定 PepDB 根目录。"
        )

    # 4）读取总表
    rows = 读取CSV(input_csv)
    print(f"成功读取记录数：{len(rows)}")
    if rows:
        print("检测到的列名如下：")
        print(list(rows[0].keys()))

    if len(rows) == 0:
        raise ValueError("输入 CSV 没有数据行，无法继续分析。")

    # 5）生成结果表
    表1 = 生成IC50类型分布表(rows)
    表2 = 生成IC50单位分布表(rows)
    表3 = 生成标签分级统计表(rows)
    表4 = 生成重复序列次数分布表(rows)
    表5 = 生成精确IC50冲突明细表(rows)
    表6 = 生成冲突与实验来源关系表(rows)

    # 6）写出结果表
    写出CSV(
        output_dir / "IC50类型分布表.csv",
        ["IC50类型", "样本数", "占比"],
        表1
    )

    写出CSV(
        output_dir / "IC50单位分布表.csv",
        ["IC50原始单位", "样本数", "占比"],
        表2
    )

    写出CSV(
        output_dir / "标签分级统计表.csv",
        ["标签等级", "判定规则", "样本数", "占全部样本比例", "有效序列样本数", "二肽样本数", "三肽样本数"],
        表3
    )

    写出CSV(
        output_dir / "重复序列次数分布表.csv",
        ["重复次数", "对应唯一序列数", "对应记录总数"],
        表4
    )

    写出CSV(
        output_dir / "精确IC50冲突明细表.csv",
        ["序列", "重复记录数", "精确IC50唯一值个数", "精确IC50值列表_uM", "来源列表", "Assay列表", "记录ID列表"],
        表5
    )

    写出CSV(
        output_dir / "冲突与实验来源关系表.csv",
        ["序列", "冲突精确IC50值列表_uM", "来源种类数", "来源列表", "Assay种类数", "Assay列表", "是否跨来源", "是否跨Assay"],
        表6
    )

    # 7）控制台打印摘要，方便你快速确认是否跑通
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