
#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
脚本名称：
    analyze_highconf_short_peptide_visuals_stdlib.py

脚本用途：
    读取 PepDB 中已经构建好的三张“高置信短肽工作集”表：
        1）高置信二三肽记录表.csv
        2）高置信二三肽序列聚合表.csv
        3）高置信二三肽稳定序列表.csv

    并自动输出你当前最需要的一组“图 + 表”分析结果，帮助你从直观层面看清：
        - 长度分布
        - 二肽 / 三肽数量
        - 精确 IC50 样本量
        - 重复序列与冲突值
        - N / C 端位点频率
        - 强活性 vs 全体富集
        - 顺序翻转对比
        - 位置化构效分析

为什么这份脚本坚持用“图 + 表”同时输出：
    1）图更适合快速发现规律、直观看趋势
    2）表更适合核对样本、追踪具体序列、写方法和结果部分
    3）你后面做建模时，需要图发现问题、表精确定位问题

脚本特点：
    - 纯 Python 标准库实现
    - 不依赖 pandas / numpy / matplotlib / seaborn
    - 图统一输出为 SVG，方便直接双击打开，也方便后续插入 PPT / 文档
    - 表统一输出为 CSV，编码为 utf-8-sig，便于 Excel 打开
    - 代码里保留了详细中文注释，后续你自己改也容易
    - 输入文件路径和输出目录放在文件前面集中配置，方便你后期修改

默认输出目录：
    DB/analysis/ace/ahtpdb/高置信短肽图表分析/

输出结果分为两类：
    一、表格（CSV）
    二、图形（SVG）

建议你后续怎么用：
    - 探索规律：先看图
    - 核对细节：再看表
    - 写汇报 / 论文：主文放图，附表留详细结果
"""

from __future__ import annotations

import argparse
import csv
import html
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple


# ======================================================================================
# 一、路径配置区（你后续最常改的是这里）
# ======================================================================================

# 是否优先使用下面手动指定的路径。
# True  = 完全按下面这些路径读写，不再依赖脚本自动推断
# False = 若你不想手工写死路径，则由脚本按 PepDB/SCR/scripts/ 的位置自动推断
PRIORITY_USE_MANUAL_PATHS = True

# PepDB 项目根目录（如果你想手动指定项目根目录，就改这里）
MANUAL_PROJECT_ROOT = r"E:\MYS\PepDB"

# 三张输入表
MANUAL_INPUT_RECORDS_CSV = r"E:\MYS\PepDB\DB\worksets\ace\ahtpdb\高置信短肽工作集\高置信二三肽记录表.csv"
MANUAL_INPUT_AGGREGATED_CSV = r"E:\MYS\PepDB\DB\worksets\ace\ahtpdb\高置信短肽工作集\高置信二三肽序列聚合表.csv"
MANUAL_INPUT_STABLE_CSV = r"E:\MYS\PepDB\DB\worksets\ace\ahtpdb\高置信短肽工作集\高置信二三肽稳定序列表.csv"

# 输出目录
MANUAL_OUTPUT_DIR = r"E:\MYS\PepDB\DB\analysis\ace\ahtpdb\高置信短肽图表分析"


# ======================================================================================
# 二、基础工具函数
# ======================================================================================

AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")


def 清洗字符串(value: Any) -> str:
    """
    将任意值安全转成字符串，并去掉首尾空白。

    说明：
    - None -> ""
    - 其他对象 -> str(value).strip()
    """
    if value is None:
        return ""
    return str(value).strip()


def 安全转浮点(value: Any) -> Optional[float]:
    """
    尝试把某个值转成 float。
    如果为空或无法解析，返回 None。
    """
    text = 清洗字符串(value)
    if text == "":
        return None
    try:
        return float(text)
    except Exception:
        return None


def 安全转整数(value: Any) -> Optional[int]:
    """
    尝试把某个值转成 int。
    对于 '2.0' 这类字符串也能兼容。
    """
    text = 清洗字符串(value)
    if text == "":
        return None
    try:
        return int(float(text))
    except Exception:
        return None


def 判定为真(value: Any) -> bool:
    """
    统一处理 CSV 中常见真假写法：
        1 / 0
        true / false
        yes / no
        y / n
        是 / 否

    返回布尔值。
    """
    text = 清洗字符串(value).lower()
    return text in {"1", "true", "t", "yes", "y", "是"}


def 格式化浮点(value: Optional[float], digits: int = 6) -> str:
    """
    将 float 格式化为适合写入 CSV 的字符串。
    None -> ""
    """
    if value is None:
        return ""
    return f"{value:.{digits}g}"


def 百分比字符串(num: int, den: int) -> str:
    """
    将分子/分母转换成百分比字符串。
    """
    if den == 0:
        return "0.00%"
    return f"{(num / den) * 100:.2f}%"


def 对数10(value: Optional[float]) -> Optional[float]:
    """
    安全计算 log10。
    只有 value > 0 时才返回结果，否则返回 None。
    """
    if value is None or value <= 0:
        return None
    try:
        return math.log10(value)
    except Exception:
        return None


def 去重并排序(values: Iterable[str], sep: str = " | ") -> str:
    """
    对一组字符串：
    1）去空
    2）去重
    3）排序
    4）拼接

    常用于来源列表、Assay 列表等。
    """
    cleaned = [清洗字符串(v) for v in values if 清洗字符串(v) != ""]
    return sep.join(sorted(set(cleaned)))


def 短肽类型排序值(peptide_type: str) -> int:
    """
    排序辅助：
        二肽 -> 0
        三肽 -> 1
        其他 -> 9
    """
    t = 清洗字符串(peptide_type)
    if t == "二肽":
        return 0
    if t == "三肽":
        return 1
    return 9


def 输出文件安全名(name: str) -> str:
    """
    清理文件名中不适合在 Windows 路径里直接使用的字符。
    """
    bad = ['\\', '/', ':', '*', '?', '"', '<', '>', '|']
    result = name
    for ch in bad:
        result = result.replace(ch, "_")
    return result


# ======================================================================================
# 三、路径解析与 CSV 读写函数
# ======================================================================================

def 自动推断项目根目录(脚本路径: Path) -> Path:
    """
    如果脚本放在 PepDB/SCR/scripts/ 目录下，
    那么项目根目录就是 scripts 的上上级目录，即 PepDB。
    """
    return 脚本路径.resolve().parent.parent.parent


def 获取默认路径(project_root: Path) -> Tuple[Path, Path, Path, Path]:
    """
    根据项目根目录，构造默认输入/输出路径。
    """
    workset_dir = project_root / "DB" / "worksets" / "ace" / "ahtpdb" / "高置信短肽工作集"
    records_csv = workset_dir / "高置信二三肽记录表.csv"
    aggregated_csv = workset_dir / "高置信二三肽序列聚合表.csv"
    stable_csv = workset_dir / "高置信二三肽稳定序列表.csv"
    output_dir = project_root / "DB" / "analysis" / "ace" / "ahtpdb" / "高置信短肽图表分析"
    return records_csv, aggregated_csv, stable_csv, output_dir


def 解析路径(args_project_root: str, script_path: Path) -> Tuple[Path, Path, Path, Path, Path]:
    """
    统一决定项目根目录、三张输入表、输出目录。

    优先级：
        1）如果 PRIORITY_USE_MANUAL_PATHS = True，则优先使用脚本顶部手写路径
        2）否则如果命令行给了 --project-root，则按它推断默认路径
        3）否则按脚本位置自动推断 PepDB 根目录
    """
    if PRIORITY_USE_MANUAL_PATHS:
        project_root = Path(MANUAL_PROJECT_ROOT).resolve() if 清洗字符串(MANUAL_PROJECT_ROOT) else 自动推断项目根目录(script_path)
        records_csv = Path(MANUAL_INPUT_RECORDS_CSV).resolve()
        aggregated_csv = Path(MANUAL_INPUT_AGGREGATED_CSV).resolve()
        stable_csv = Path(MANUAL_INPUT_STABLE_CSV).resolve()
        output_dir = Path(MANUAL_OUTPUT_DIR).resolve()
        return project_root, records_csv, aggregated_csv, stable_csv, output_dir

    if 清洗字符串(args_project_root):
        project_root = Path(args_project_root).resolve()
    else:
        project_root = 自动推断项目根目录(script_path)

    records_csv, aggregated_csv, stable_csv, output_dir = 获取默认路径(project_root)
    return project_root, records_csv, aggregated_csv, stable_csv, output_dir


def 读取CSV(filepath: Path) -> Tuple[List[Dict[str, str]], List[str], str]:
    """
    读取 CSV 文件，并自动尝试多种常见编码。

    返回：
        rows        -> 数据行（字典列表）
        fieldnames  -> 读取到的列名
        encoding    -> 最终成功读取所用编码

    为什么要写成这样：
    你前面已经遇到过 Windows / Excel 编码问题。
    为了让脚本更稳，这里会按顺序尝试：
        utf-8-sig
        utf-8
        gb18030
        gbk
        cp936
        latin1
    """
    encodings = ["utf-8-sig", "utf-8", "gb18030", "gbk", "cp936", "latin1"]
    errors = []

    for enc in encodings:
        try:
            with filepath.open("r", encoding=enc, newline="") as f:
                sample = f.read(4096)
                f.seek(0)

                try:
                    dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
                except Exception:
                    dialect = csv.excel

                reader = csv.DictReader(f, dialect=dialect)
                fieldnames = [清洗字符串(x) for x in (reader.fieldnames or [])]

                rows: List[Dict[str, str]] = []
                for row in reader:
                    cleaned = {}
                    for k, v in row.items():
                        cleaned[清洗字符串(k)] = 清洗字符串(v)
                    rows.append(cleaned)

                return rows, fieldnames, enc
        except Exception as e:
            errors.append(f"{enc}: {type(e).__name__} -> {e}")

    raise RuntimeError(
        f"无法成功读取文件：{filepath}\n"
        + "\n".join(errors)
    )


def 写出CSV(filepath: Path, fieldnames: List[str], rows: List[Dict[str, Any]]) -> None:
    """
    写出 CSV，统一使用 utf-8-sig，便于 Excel 打开。
    """
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with filepath.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def 写出文本(filepath: Path, text: str) -> None:
    """
    写出普通文本文件。
    """
    filepath.parent.mkdir(parents=True, exist_ok=True)
    filepath.write_text(text, encoding="utf-8-sig")


def 检查必要列(rows: List[Dict[str, str]], required_cols: List[str], file_label: str) -> None:
    """
    检查输入表中是否含有必要列。
    如果缺列，尽早报错，避免后续算到一半才崩。
    """
    if not rows:
        raise ValueError(f"{file_label} 没有任何数据行。")

    actual = set(rows[0].keys())
    missing = [c for c in required_cols if c not in actual]

    if missing:
        raise ValueError(
            f"{file_label} 缺少以下必要列：\n"
            + ", ".join(missing)
            + "\n\n实际读取到的列名为：\n"
            + ", ".join(rows[0].keys())
        )


# ======================================================================================
# 四、SVG 绘图基础函数（纯标准库版）
# ======================================================================================

def svg_escape(text: Any) -> str:
    """
    对文字做 XML / SVG 转义，防止特殊字符破坏 SVG。
    """
    return html.escape(清洗字符串(text))


def 颜色插值(c1: Tuple[int, int, int], c2: Tuple[int, int, int], t: float) -> str:
    """
    在线性空间内进行颜色插值。
    t 取值范围建议 0~1。
    返回形如 '#RRGGBB' 的颜色字符串。
    """
    t = max(0.0, min(1.0, t))
    r = round(c1[0] + (c2[0] - c1[0]) * t)
    g = round(c1[1] + (c2[1] - c1[1]) * t)
    b = round(c1[2] + (c2[2] - c1[2]) * t)
    return f"#{r:02x}{g:02x}{b:02x}"


def 单向配色(value: float, vmax: float) -> str:
    """
    用于非负数值（例如频率、计数）的颜色映射。
    0 -> 近白色
    vmax -> 深蓝色
    """
    if vmax <= 0:
        return "#f7fbff"
    ratio = max(0.0, min(1.0, value / vmax))
    return 颜色插值((247, 251, 255), (8, 81, 156), ratio)


def 双向配色(value: float, vabs: float) -> str:
    """
    用于带正负号的数值（例如 log2FC）的颜色映射。
    负值 -> 蓝
    0   -> 白
    正值 -> 红
    """
    if vabs <= 0:
        return "#ffffff"

    x = max(-vabs, min(vabs, value))
    if x >= 0:
        return 颜色插值((255, 255, 255), (203, 24, 29), x / vabs)
    else:
        return 颜色插值((33, 113, 181), (255, 255, 255), (x + vabs) / vabs)


def 保存SVG(filepath: Path, width: int, height: int, elements: List[str]) -> None:
    """
    将 SVG 元素列表保存为完整 SVG 文件。
    """
    filepath.parent.mkdir(parents=True, exist_ok=True)
    content = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect x="0" y="0" width="100%" height="100%" fill="white"/>'
    ]
    content.extend(elements)
    content.append("</svg>")
    filepath.write_text("\n".join(content), encoding="utf-8")


def 画标题(elements: List[str], title: str, width: int) -> None:
    """
    在 SVG 顶部中央写标题。
    """
    elements.append(
        f'<text x="{width/2:.1f}" y="34" text-anchor="middle" '
        f'font-size="22" font-weight="bold" font-family="Arial, Helvetica, sans-serif">'
        f'{svg_escape(title)}</text>'
    )


def 画脚注(elements: List[str], note: str, width: int, height: int) -> None:
    """
    在底部增加简短说明。
    """
    if 清洗字符串(note) == "":
        return
    elements.append(
        f'<text x="{width/2:.1f}" y="{height-12}" text-anchor="middle" '
        f'font-size="11" fill="#666" font-family="Arial, Helvetica, sans-serif">'
        f'{svg_escape(note)}</text>'
    )


def 绘制柱状图(
    filepath: Path,
    labels: List[str],
    values: List[float],
    title: str,
    x_label: str = "",
    y_label: str = "",
    color: str = "#2E86AB",
    value_digits: int = 4,
    rotate_x: bool = False,
    note: str = ""
) -> None:
    """
    绘制简单竖向柱状图。

    适用场景：
        - 长度分布
        - 精确 IC50 样本量
        - 重复次数分布
        - 冲突倍数分箱统计
    """
    width = 1200
    height = 760
    left = 100
    right = 40
    top = 70
    bottom = 150 if rotate_x else 120
    plot_w = width - left - right
    plot_h = height - top - bottom

    elements: List[str] = []
    画标题(elements, title, width)

    max_value = max(values) if values else 1.0
    if max_value <= 0:
        max_value = 1.0

    # 网格线与 y 轴刻度
    for i in range(6):
        y = top + plot_h - i * plot_h / 5
        tick_val = max_value * i / 5
        elements.append(f'<line x1="{left}" y1="{y:.1f}" x2="{left+plot_w}" y2="{y:.1f}" stroke="#dddddd" stroke-width="1"/>')
        elements.append(
            f'<text x="{left-12}" y="{y+4:.1f}" text-anchor="end" font-size="12" fill="#444" '
            f'font-family="Arial, Helvetica, sans-serif">{svg_escape(格式化浮点(tick_val, 6))}</text>'
        )

    # 坐标轴
    elements.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top+plot_h}" stroke="#333" stroke-width="1.5"/>')
    elements.append(f'<line x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" y2="{top+plot_h}" stroke="#333" stroke-width="1.5"/>')

    n = len(labels) if labels else 1
    slot = plot_w / max(n, 1)
    bar_w = min(54, slot * 0.65)

    for i, (lab, val) in enumerate(zip(labels, values)):
        cx = left + i * slot + slot / 2
        bar_h = 0 if max_value == 0 else (val / max_value) * plot_h
        x = cx - bar_w / 2
        y = top + plot_h - bar_h

        elements.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w:.1f}" height="{bar_h:.1f}" fill="{color}" opacity="0.88"/>')
        elements.append(
            f'<text x="{cx:.1f}" y="{y-6:.1f}" text-anchor="middle" font-size="11" fill="#222" '
            f'font-family="Arial, Helvetica, sans-serif">{svg_escape(格式化浮点(val, value_digits))}</text>'
        )

        if rotate_x:
            elements.append(
                f'<text x="{cx:.1f}" y="{top+plot_h+18:.1f}" transform="rotate(-45 {cx:.1f},{top+plot_h+18:.1f})" '
                f'text-anchor="end" font-size="12" fill="#222" font-family="Arial, Helvetica, sans-serif">{svg_escape(lab)}</text>'
            )
        else:
            elements.append(
                f'<text x="{cx:.1f}" y="{top+plot_h+22:.1f}" text-anchor="middle" font-size="12" fill="#222" '
                f'font-family="Arial, Helvetica, sans-serif">{svg_escape(lab)}</text>'
            )

    # 坐标轴标题
    if x_label:
        elements.append(
            f'<text x="{left + plot_w/2:.1f}" y="{height-58:.1f}" text-anchor="middle" font-size="14" '
            f'fill="#222" font-family="Arial, Helvetica, sans-serif">{svg_escape(x_label)}</text>'
        )
    if y_label:
        elements.append(
            f'<text x="24" y="{top + plot_h/2:.1f}" transform="rotate(-90 24,{top + plot_h/2:.1f})" '
            f'text-anchor="middle" font-size="14" fill="#222" font-family="Arial, Helvetica, sans-serif">{svg_escape(y_label)}</text>'
        )

    画脚注(elements, note, width, height)
    保存SVG(filepath, width, height, elements)


def 绘制分组柱状图(
    filepath: Path,
    categories: List[str],
    series: Dict[str, List[float]],
    title: str,
    x_label: str = "",
    y_label: str = "",
    note: str = ""
) -> None:
    """
    绘制分组柱状图。

    适用场景：
        - 不同数据层级的长度分布
        - 二肽 / 三肽数量对比
        - 精确 IC50 样本量对比
    """
    width = 1280
    height = 780
    left = 100
    right = 60
    top = 90
    bottom = 140
    plot_w = width - left - right
    plot_h = height - top - bottom

    colors = ["#2E86AB", "#F18F01", "#A23B72", "#4CAF50", "#C73E1D", "#6A4C93"]

    elements: List[str] = []
    画标题(elements, title, width)

    all_values = []
    for vals in series.values():
        all_values.extend(vals)
    max_value = max(all_values) if all_values else 1.0
    if max_value <= 0:
        max_value = 1.0

    # 网格线
    for i in range(6):
        y = top + plot_h - i * plot_h / 5
        tick_val = max_value * i / 5
        elements.append(f'<line x1="{left}" y1="{y:.1f}" x2="{left+plot_w}" y2="{y:.1f}" stroke="#dddddd" stroke-width="1"/>')
        elements.append(
            f'<text x="{left-12}" y="{y+4:.1f}" text-anchor="end" font-size="12" fill="#444" '
            f'font-family="Arial, Helvetica, sans-serif">{svg_escape(格式化浮点(tick_val, 6))}</text>'
        )

    elements.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top+plot_h}" stroke="#333" stroke-width="1.5"/>')
    elements.append(f'<line x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" y2="{top+plot_h}" stroke="#333" stroke-width="1.5"/>')

    cat_n = max(1, len(categories))
    group_slot = plot_w / cat_n
    series_names = list(series.keys())
    group_bar_area = min(120, group_slot * 0.8)
    bar_w = group_bar_area / max(1, len(series_names))

    for i, cat in enumerate(categories):
        group_center = left + i * group_slot + group_slot / 2
        group_left = group_center - group_bar_area / 2

        for j, sname in enumerate(series_names):
            val = series[sname][i]
            bar_h = (val / max_value) * plot_h if max_value > 0 else 0
            x = group_left + j * bar_w + 4
            y = top + plot_h - bar_h
            color = colors[j % len(colors)]

            elements.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w-8:.1f}" height="{bar_h:.1f}" fill="{color}" opacity="0.88"/>')
            elements.append(
                f'<text x="{x + (bar_w-8)/2:.1f}" y="{y-6:.1f}" text-anchor="middle" font-size="10" fill="#222" '
                f'font-family="Arial, Helvetica, sans-serif">{svg_escape(格式化浮点(val, 5))}</text>'
            )

        elements.append(
            f'<text x="{group_center:.1f}" y="{top+plot_h+24:.1f}" text-anchor="middle" font-size="12" fill="#222" '
            f'font-family="Arial, Helvetica, sans-serif">{svg_escape(cat)}</text>'
        )

    # 图例
    legend_x = left + 20
    legend_y = 54
    for idx, sname in enumerate(series_names):
        color = colors[idx % len(colors)]
        x = legend_x + idx * 180
        elements.append(f'<rect x="{x}" y="{legend_y}" width="18" height="12" fill="{color}" opacity="0.88"/>')
        elements.append(
            f'<text x="{x+24}" y="{legend_y+11}" font-size="12" fill="#222" '
            f'font-family="Arial, Helvetica, sans-serif">{svg_escape(sname)}</text>'
        )

    if x_label:
        elements.append(
            f'<text x="{left + plot_w/2:.1f}" y="{height-62:.1f}" text-anchor="middle" font-size="14" fill="#222" '
            f'font-family="Arial, Helvetica, sans-serif">{svg_escape(x_label)}</text>'
        )
    if y_label:
        elements.append(
            f'<text x="24" y="{top + plot_h/2:.1f}" transform="rotate(-90 24,{top + plot_h/2:.1f})" '
            f'text-anchor="middle" font-size="14" fill="#222" font-family="Arial, Helvetica, sans-serif">{svg_escape(y_label)}</text>'
        )

    画脚注(elements, note, width, height)
    保存SVG(filepath, width, height, elements)


def 绘制热图(
    filepath: Path,
    row_labels: List[str],
    col_labels: List[str],
    matrix: List[List[float]],
    title: str,
    diverging: bool = False,
    note: str = "",
    value_digits: int = 3
) -> None:
    """
    绘制热图。

    适用场景：
        - N / C 端位点频率热图
        - 二肽 / 三肽位置频率热图
        - 强活性 vs 全体的位点富集热图

    参数说明：
        diverging = False
            用于非负数据，例如频率
        diverging = True
            用于有正负号的数据，例如 log2FC
    """
    width = max(1100, 160 + len(col_labels) * 46)
    height = max(420, 160 + len(row_labels) * 56)
    left = 150
    top = 110
    bottom = 120
    right = 100

    plot_w = width - left - right
    plot_h = height - top - bottom

    cell_w = plot_w / max(1, len(col_labels))
    cell_h = plot_h / max(1, len(row_labels))

    elements: List[str] = []
    画标题(elements, title, width)

    if diverging:
        all_values = [abs(v) for row in matrix for v in row]
        vmax = max(all_values) if all_values else 1.0
        vmax = max(vmax, 1e-12)
    else:
        all_values = [v for row in matrix for v in row]
        vmax = max(all_values) if all_values else 1.0
        vmax = max(vmax, 1e-12)

    # 单元格
    for i, rlab in enumerate(row_labels):
        y = top + i * cell_h

        # 行标签
        elements.append(
            f'<text x="{left-12}" y="{y + cell_h/2 + 5:.1f}" text-anchor="end" font-size="13" fill="#222" '
            f'font-family="Arial, Helvetica, sans-serif">{svg_escape(rlab)}</text>'
        )

        for j, clab in enumerate(col_labels):
            x = left + j * cell_w
            value = matrix[i][j]
            color = 双向配色(value, vmax) if diverging else 单向配色(value, vmax)

            elements.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{cell_w:.1f}" height="{cell_h:.1f}" fill="{color}" stroke="#ffffff" stroke-width="1"/>')

            # 在单元格中央写数值
            text_color = "#111" if not diverging else ("#111" if abs(value) < vmax * 0.55 else "#fff")
            elements.append(
                f'<text x="{x + cell_w/2:.1f}" y="{y + cell_h/2 + 4:.1f}" text-anchor="middle" font-size="10" fill="{text_color}" '
                f'font-family="Arial, Helvetica, sans-serif">{svg_escape(格式化浮点(value, value_digits))}</text>'
            )

    # 列标签
    for j, clab in enumerate(col_labels):
        x = left + j * cell_w + cell_w / 2
        elements.append(
            f'<text x="{x:.1f}" y="{top-12:.1f}" text-anchor="middle" font-size="12" fill="#222" '
            f'font-family="Arial, Helvetica, sans-serif">{svg_escape(clab)}</text>'
        )

    # 图例
    legend_x = width - 72
    legend_y = top
    legend_h = plot_h
    steps = 100
    for k in range(steps):
        t = k / (steps - 1)
        if diverging:
            val = (1 - 2 * t) * vmax  # 从 +vmax 到 -vmax
            color = 双向配色(val, vmax)
        else:
            val = (1 - t) * vmax
            color = 单向配色(val, vmax)
        y = legend_y + k * legend_h / steps
        elements.append(f'<rect x="{legend_x}" y="{y:.1f}" width="18" height="{legend_h/steps + 1:.2f}" fill="{color}" stroke="none"/>')

    if diverging:
        legend_labels = [(legend_y, vmax), (legend_y + legend_h/2, 0.0), (legend_y + legend_h, -vmax)]
    else:
        legend_labels = [(legend_y, vmax), (legend_y + legend_h, 0.0)]

    for y, val in legend_labels:
        elements.append(
            f'<text x="{legend_x + 26}" y="{y+4:.1f}" font-size="11" fill="#222" '
            f'font-family="Arial, Helvetica, sans-serif">{svg_escape(格式化浮点(val, 4))}</text>'
        )

    画脚注(elements, note, width, height)
    保存SVG(filepath, width, height, elements)


def 绘制发散条形图(
    filepath: Path,
    labels: List[str],
    values: List[float],
    title: str,
    x_label: str = "log2FC",
    note: str = ""
) -> None:
    """
    绘制横向发散条形图。

    适用场景：
        - 强活性 vs 全体的富集差异
        - 某些位置化富集的 Top residue 比较

    说明：
        正值用红色（在强活性中更富集）
        负值用蓝色（在全体中更常见或在强活性中相对欠富集）
    """
    n = len(labels)
    width = 1200
    height = max(520, 120 + n * 30)
    left = 220
    right = 80
    top = 90
    bottom = 80

    plot_w = width - left - right
    plot_h = height - top - bottom

    elements: List[str] = []
    画标题(elements, title, width)

    vmax = max([abs(v) for v in values], default=1.0)
    vmax = max(vmax, 1e-6)

    zero_x = left + plot_w / 2

    # 网格与刻度
    ticks = [-vmax, -vmax/2, 0.0, vmax/2, vmax]
    for tv in ticks:
        x = zero_x + (tv / vmax) * (plot_w / 2)
        elements.append(f'<line x1="{x:.1f}" y1="{top}" x2="{x:.1f}" y2="{top+plot_h}" stroke="#e0e0e0" stroke-width="1"/>')
        elements.append(
            f'<text x="{x:.1f}" y="{top+plot_h+22:.1f}" text-anchor="middle" font-size="12" fill="#444" '
            f'font-family="Arial, Helvetica, sans-serif">{svg_escape(格式化浮点(tv, 4))}</text>'
        )

    elements.append(f'<line x1="{zero_x:.1f}" y1="{top}" x2="{zero_x:.1f}" y2="{top+plot_h}" stroke="#444" stroke-width="1.5"/>')

    row_h = plot_h / max(1, n)

    for i, (lab, val) in enumerate(zip(labels, values)):
        cy = top + i * row_h + row_h / 2
        bar_len = abs(val) / vmax * (plot_w / 2)
        bar_h = min(18, row_h * 0.65)

        if val >= 0:
            x = zero_x
            color = "#C81D25"
            text_anchor = "start"
            text_x = x + bar_len + 8
        else:
            x = zero_x - bar_len
            color = "#1F78B4"
            text_anchor = "end"
            text_x = x - 8

        elements.append(f'<rect x="{x:.1f}" y="{cy-bar_h/2:.1f}" width="{bar_len:.1f}" height="{bar_h:.1f}" fill="{color}" opacity="0.88"/>')
        elements.append(
            f'<text x="{left-12}" y="{cy+4:.1f}" text-anchor="end" font-size="12" fill="#222" '
            f'font-family="Arial, Helvetica, sans-serif">{svg_escape(lab)}</text>'
        )
        elements.append(
            f'<text x="{text_x:.1f}" y="{cy+4:.1f}" text-anchor="{text_anchor}" font-size="11" fill="#222" '
            f'font-family="Arial, Helvetica, sans-serif">{svg_escape(格式化浮点(val, 4))}</text>'
        )

    if x_label:
        elements.append(
            f'<text x="{left + plot_w/2:.1f}" y="{height-24:.1f}" text-anchor="middle" font-size="14" fill="#222" '
            f'font-family="Arial, Helvetica, sans-serif">{svg_escape(x_label)}</text>'
        )

    画脚注(elements, note, width, height)
    保存SVG(filepath, width, height, elements)


def 绘制顺序翻转散点图(
    filepath: Path,
    pair_rows: List[Dict[str, Any]],
    title: str,
    note: str = ""
) -> None:
    """
    绘制“顺序翻转对比散点图”。

    横轴：
        序列A 的中位数 IC50（log10）
    纵轴：
        序列B 的中位数 IC50（log10）

    为什么用 log10：
        IC50 往往跨数量级，用原值画图会严重挤在一侧；
        用 log10 更适合比较相对差异。
    """
    # 提取有效点
    points = []
    for row in pair_rows:
        xa = 安全转浮点(row.get("序列A中位数IC50_uM", ""))
        yb = 安全转浮点(row.get("序列B中位数IC50_uM", ""))
        if xa is None or yb is None or xa <= 0 or yb <= 0:
            continue
        lx = math.log10(xa)
        ly = math.log10(yb)
        label = f'{row.get("序列A", "")}↔{row.get("序列B", "")}'
        diff = abs(lx - ly)
        points.append((label, lx, ly, diff))

    width = 900
    height = 760
    left = 100
    right = 60
    top = 90
    bottom = 100
    plot_w = width - left - right
    plot_h = height - top - bottom

    elements: List[str] = []
    画标题(elements, title, width)

    if not points:
        elements.append('<text x="450" y="380" text-anchor="middle" font-size="20" fill="#444">没有可用于绘图的顺序翻转配对样本</text>')
        画脚注(elements, note, width, height)
        保存SVG(filepath, width, height, elements)
        return

    xs = [p[1] for p in points]
    ys = [p[2] for p in points]
    vmin = min(xs + ys)
    vmax = max(xs + ys)
    if abs(vmax - vmin) < 1e-9:
        vmax = vmin + 1

    pad = (vmax - vmin) * 0.08
    vmin -= pad
    vmax += pad

    def sx(v: float) -> float:
        return left + (v - vmin) / (vmax - vmin) * plot_w

    def sy(v: float) -> float:
        return top + plot_h - (v - vmin) / (vmax - vmin) * plot_h

    # 网格与刻度
    for i in range(6):
        frac = i / 5
        val = vmin + frac * (vmax - vmin)
        x = sx(val)
        y = sy(val)
        elements.append(f'<line x1="{x:.1f}" y1="{top}" x2="{x:.1f}" y2="{top+plot_h}" stroke="#efefef" stroke-width="1"/>')
        elements.append(f'<line x1="{left}" y1="{y:.1f}" x2="{left+plot_w}" y2="{y:.1f}" stroke="#efefef" stroke-width="1"/>')
        elements.append(
            f'<text x="{x:.1f}" y="{top+plot_h+22:.1f}" text-anchor="middle" font-size="12" fill="#444" '
            f'font-family="Arial, Helvetica, sans-serif">{svg_escape(格式化浮点(val, 4))}</text>'
        )
        elements.append(
            f'<text x="{left-12}" y="{y+4:.1f}" text-anchor="end" font-size="12" fill="#444" '
            f'font-family="Arial, Helvetica, sans-serif">{svg_escape(格式化浮点(val, 4))}</text>'
        )

    # 坐标轴框线
    elements.append(f'<rect x="{left}" y="{top}" width="{plot_w}" height="{plot_h}" fill="none" stroke="#333" stroke-width="1.4"/>')

    # 45度对角线：表示两者活性完全相等
    elements.append(
        f'<line x1="{sx(vmin):.1f}" y1="{sy(vmin):.1f}" x2="{sx(vmax):.1f}" y2="{sy(vmax):.1f}" '
        f'stroke="#C81D25" stroke-width="1.2" stroke-dasharray="6,4"/>'
    )

    # 点
    top_labels = sorted(points, key=lambda x: x[3], reverse=True)[:10]
    top_label_names = {p[0] for p in top_labels}

    for label, lx, ly, diff in points:
        x = sx(lx)
        y = sy(ly)
        elements.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4.4" fill="#2E86AB" opacity="0.8"/>')
        if label in top_label_names:
            elements.append(
                f'<text x="{x+8:.1f}" y="{y-8:.1f}" font-size="10" fill="#222" '
                f'font-family="Arial, Helvetica, sans-serif">{svg_escape(label)}</text>'
            )

    # 轴标题
    elements.append(
        f'<text x="{left + plot_w/2:.1f}" y="{height-28:.1f}" text-anchor="middle" font-size="14" fill="#222" '
        f'font-family="Arial, Helvetica, sans-serif">序列A 的 log10(IC50_uM)</text>'
    )
    elements.append(
        f'<text x="24" y="{top + plot_h/2:.1f}" transform="rotate(-90 24,{top + plot_h/2:.1f})" '
        f'text-anchor="middle" font-size="14" fill="#222" font-family="Arial, Helvetica, sans-serif">序列B 的 log10(IC50_uM)</text>'
    )

    画脚注(elements, note, width, height)
    保存SVG(filepath, width, height, elements)


# ======================================================================================
# 五、数据理解与辅助判定函数
# ======================================================================================

def 是否强活性_稳定表(row: Dict[str, str]) -> bool:
    """
    判定稳定序列表中的一条序列是否属于“强活性”。

    优先使用：
        活性分层建议 == 强活性

    如果该字段缺失，则回退为：
        精确IC50中位数_uM <= 50
    """
    label = 清洗字符串(row.get("活性分层建议", ""))
    if label != "":
        return label == "强活性"

    value = 安全转浮点(row.get("精确IC50中位数_uM", ""))
    return value is not None and value <= 50


def 是否强活性_记录表(row: Dict[str, str]) -> bool:
    """
    判定记录表中的一条记录是否属于“强活性”。

    这里用单条记录的 精确IC50_uM <= 50 作为第一版简单判定。
    """
    value = 安全转浮点(row.get("精确IC50_uM", ""))
    return value is not None and value <= 50


def 是否强活性_聚合表(row: Dict[str, str]) -> bool:
    """
    判定序列聚合表中的一条序列是否属于“强活性”。

    这里使用 精确IC50中位数_uM <= 50 作为第一版简单判定。
    """
    value = 安全转浮点(row.get("精确IC50中位数_uM", ""))
    return value is not None and value <= 50


def 获取肽长(row: Dict[str, str]) -> Optional[int]:
    """
    统一从不同表中读取肽长。
    """
    return 安全转整数(row.get("肽长", ""))


def 获取短肽类型(row: Dict[str, str]) -> str:
    """
    统一读取短肽类型。
    """
    return 清洗字符串(row.get("短肽类型", ""))


def 获取稳定性等级(row: Dict[str, str]) -> str:
    """
    统一读取稳定性等级。
    """
    return 清洗字符串(row.get("稳定性等级", ""))


# ======================================================================================
# 六、分析 1：长度分布
# ======================================================================================

def 生成长度分布统计表(records_rows, aggregated_rows, stable_rows) -> List[Dict[str, Any]]:
    """
    生成长度分布统计表。

    表的目的：
        从三个层级同时看“肽长 = 2 / 3”的分布情况：
        - 高置信记录层
        - 唯一序列层
        - 稳定序列层

    虽然你现在工作集本质上就是二肽和三肽，
    但把三个层级放在一起看，能直观看出：
        - 去重带来的压缩程度
        - 稳定筛选后各长度的保留程度
    """
    layer_map = {
        "高置信记录": records_rows,
        "唯一序列": aggregated_rows,
        "稳定序列": stable_rows,
    }

    all_lengths = sorted({
        x for rows in layer_map.values()
        for x in [获取肽长(r) for r in rows]
        if x is not None
    })

    result = []
    for layer_name, rows in layer_map.items():
        total = len(rows)
        counter = Counter()
        for row in rows:
            length = 获取肽长(row)
            if length is not None:
                counter[length] += 1

        for length in all_lengths:
            count = counter.get(length, 0)
            result.append({
                "数据层级": layer_name,
                "肽长": length,
                "样本数": count,
                "占该层级比例": 百分比字符串(count, total)
            })

    return result


def 生成长度分布图数据(长度分布表: List[Dict[str, Any]]) -> Tuple[List[str], Dict[str, List[float]]]:
    """
    将长度分布统计表转成“分组柱状图”所需的数据格式。
    """
    lengths = sorted(set(安全转整数(r["肽长"]) for r in 长度分布表 if 安全转整数(r["肽长"]) is not None))
    categories = [f"{x}肽" for x in lengths]

    layers = ["高置信记录", "唯一序列", "稳定序列"]
    series = {layer: [] for layer in layers}

    for layer in layers:
        mapping = {}
        for row in 长度分布表:
            if row["数据层级"] == layer:
                length = 安全转整数(row["肽长"])
                count = 安全转浮点(row["样本数"]) or 0
                mapping[length] = count
        for length in lengths:
            series[layer].append(mapping.get(length, 0))

    return categories, series


# ======================================================================================
# 七、分析 2：二肽 / 三肽数量
# ======================================================================================

def 生成二肽三肽数量统计表(stable_rows) -> List[Dict[str, Any]]:
    """
    生成“二肽 / 三肽数量统计表”。

    这里不再按数据层级比较，而是按“活性层级”比较稳定集中的二肽与三肽数量。
    这样做的意义：
        - 更贴近后续建模的核心数据集
        - 更容易看出三肽是否更偏向强活性
    """
    groups = {
        "全部稳定序列": stable_rows,
        "强活性稳定序列": [r for r in stable_rows if 是否强活性_稳定表(r)],
        "中等活性稳定序列": [r for r in stable_rows if 清洗字符串(r.get("活性分层建议", "")) == "中等活性"],
        "较弱活性稳定序列": [r for r in stable_rows if 清洗字符串(r.get("活性分层建议", "")) == "较弱活性"],
    }

    result = []
    for group_name, rows in groups.items():
        total = len(rows)
        di = sum(1 for r in rows if 获取短肽类型(r) == "二肽")
        tri = sum(1 for r in rows if 获取短肽类型(r) == "三肽")

        result.append({
            "分析分组": group_name,
            "二肽数量": di,
            "二肽占比": 百分比字符串(di, total),
            "三肽数量": tri,
            "三肽占比": 百分比字符串(tri, total),
            "总数": total,
        })

    return result


def 生成二肽三肽数量图数据(数量表: List[Dict[str, Any]]) -> Tuple[List[str], Dict[str, List[float]]]:
    """
    将二肽 / 三肽数量统计表转成分组柱状图所需的数据格式。
    """
    categories = [row["分析分组"] for row in 数量表]
    series = {
        "二肽": [安全转浮点(row["二肽数量"]) or 0 for row in 数量表],
        "三肽": [安全转浮点(row["三肽数量"]) or 0 for row in 数量表],
    }
    return categories, series


# ======================================================================================
# 八、分析 3：精确 IC50 样本量
# ======================================================================================

def 生成精确IC50样本量统计表(records_rows, aggregated_rows, stable_rows) -> List[Dict[str, Any]]:
    """
    生成“精确 IC50 样本量统计表”。

    这里虽然三张工作集都已经是高置信 / 精确 IC50 导向的表，
    但它们分别对应：
        - 记录层
        - 唯一序列层
        - 稳定序列层
    所以把三层的样本量写清楚，后续做方法部分和图注很方便。
    """
    layers = [
        ("高置信记录", records_rows, 是否强活性_记录表),
        ("唯一序列", aggregated_rows, 是否强活性_聚合表),
        ("稳定序列", stable_rows, 是否强活性_稳定表),
    ]

    result = []
    for layer_name, rows, strong_fn in layers:
        total = len(rows)
        di = sum(1 for r in rows if 获取短肽类型(r) == "二肽")
        tri = sum(1 for r in rows if 获取短肽类型(r) == "三肽")
        strong = sum(1 for r in rows if strong_fn(r))

        result.append({
            "数据层级": layer_name,
            "全部样本量": total,
            "二肽样本量": di,
            "三肽样本量": tri,
            "强活性样本量": strong,
            "强活性占比": 百分比字符串(strong, total),
        })

    return result


def 生成精确IC50样本量图数据(样本量表: List[Dict[str, Any]]) -> Tuple[List[str], Dict[str, List[float]]]:
    """
    将“精确 IC50 样本量统计表”转成分组柱状图数据。
    """
    categories = [row["数据层级"] for row in 样本量表]
    series = {
        "全部样本量": [安全转浮点(row["全部样本量"]) or 0 for row in 样本量表],
        "强活性样本量": [安全转浮点(row["强活性样本量"]) or 0 for row in 样本量表],
    }
    return categories, series


# ======================================================================================
# 九、分析 4：重复序列与冲突值
# ======================================================================================

def 生成重复序列次数分布表(aggregated_rows) -> List[Dict[str, Any]]:
    """
    统计“一个唯一序列对应多少条高置信记录”。

    这反映的是：
        - 某些序列在数据库中是否被重复报道很多次
        - 后续做聚合、稳定性分析时，数据支持是否足够丰富
    """
    counts = Counter()
    total = len(aggregated_rows)

    for row in aggregated_rows:
        n = 安全转整数(row.get("记录数", ""))
        if n is not None:
            counts[n] += 1

    result = []
    for record_count in sorted(counts.keys()):
        seq_num = counts[record_count]
        result.append({
            "重复记录数": record_count,
            "唯一序列数": seq_num,
            "占全部唯一序列比例": 百分比字符串(seq_num, total),
        })

    return result


def 冲突倍数分箱(value: Optional[float]) -> str:
    """
    将“最大最小倍数”映射到更容易展示的分箱。
    """
    if value is None:
        return "缺失"
    if value <= 1:
        return "1（无冲突）"
    if value <= 1.5:
        return "(1,1.5]"
    if value <= 2:
        return "(1.5,2]"
    if value <= 3:
        return "(2,3]"
    if value <= 5:
        return "(3,5]"
    if value <= 10:
        return "(5,10]"
    return ">10"


def 生成冲突倍数分箱统计表(aggregated_rows) -> List[Dict[str, Any]]:
    """
    对“最大最小倍数”做分箱统计。

    这张表比直接列原始倍数更适合做整体图，
    因为它能快速回答：
        - 多少序列几乎无冲突
        - 多少序列只是轻微波动
        - 多少序列属于高冲突
    """
    counter = Counter()
    total = len(aggregated_rows)

    for row in aggregated_rows:
        ratio = 安全转浮点(row.get("最大最小倍数", ""))
        bucket = 冲突倍数分箱(ratio)
        counter[bucket] += 1

    order = ["1（无冲突）", "(1,1.5]", "(1.5,2]", "(2,3]", "(3,5]", "(5,10]", ">10", "缺失"]
    result = []
    for bucket in order:
        if bucket in counter:
            count = counter[bucket]
            result.append({
                "冲突倍数分箱": bucket,
                "序列数": count,
                "占全部唯一序列比例": 百分比字符串(count, total),
            })

    return result


def 生成高冲突序列明细表(aggregated_rows) -> List[Dict[str, Any]]:
    """
    输出“高冲突序列明细表”。

    这张表的目标不是看整体，而是让你快速定位：
        - 哪些序列冲突最严重
        - 它们是不是跨来源 / 跨 Assay
        - 后续哪些序列该排除、哪些值得单独研究

    这里收录规则：
        - 稳定性等级不是“稳定序列”
          或
        - 最大最小倍数 > 3
    """
    result = []
    for row in aggregated_rows:
        ratio = 安全转浮点(row.get("最大最小倍数", ""))
        stable_level = 获取稳定性等级(row)

        if stable_level != "稳定序列" or (ratio is not None and ratio > 3):
            result.append({
                "序列": row.get("序列", ""),
                "短肽类型": row.get("短肽类型", ""),
                "记录数": row.get("记录数", ""),
                "精确IC50唯一值个数": row.get("精确IC50唯一值个数", ""),
                "精确IC50最小值_uM": row.get("精确IC50最小值_uM", ""),
                "精确IC50最大值_uM": row.get("精确IC50最大值_uM", ""),
                "精确IC50中位数_uM": row.get("精确IC50中位数_uM", ""),
                "最大最小倍数": row.get("最大最小倍数", ""),
                "稳定性等级": row.get("稳定性等级", ""),
                "稳定性判定说明": row.get("稳定性判定说明", ""),
                "来源种类数": row.get("来源种类数", ""),
                "Assay种类数": row.get("Assay种类数", ""),
                "来源列表": row.get("来源列表", ""),
                "Assay列表": row.get("Assay列表", ""),
            })

    result.sort(
        key=lambda x: (
            - (安全转浮点(x["最大最小倍数"]) or -1),
            - (安全转整数(x["记录数"]) or 0),
            短肽类型排序值(x["短肽类型"]),
            x["序列"]
        )
    )
    return result


# ======================================================================================
# 十、分析 5：N / C 端位点频率
# ======================================================================================

def 位点频率统计(rows: List[Dict[str, str]], residue_col: str, group_name: str) -> List[Dict[str, Any]]:
    """
    对某一组样本，在指定端点位点上统计氨基酸频率。

    residue_col 只能是：
        - N端氨基酸
        - C端氨基酸
    """
    counter = Counter()
    total = 0

    for row in rows:
        aa = 清洗字符串(row.get(residue_col, ""))
        if aa in AA_ORDER:
            counter[aa] += 1
            total += 1

    result = []
    for aa in AA_ORDER:
        count = counter.get(aa, 0)
        freq = (count / total) if total > 0 else 0.0
        result.append({
            "分析分组": group_name,
            "氨基酸": aa,
            "频数": count,
            "频率": 格式化浮点(freq, 6),
        })
    return result


def 生成端点位点频率表(stable_rows, residue_col: str) -> List[Dict[str, Any]]:
    """
    生成 N 端或 C 端位点频率表。

    输出分组固定为四类：
        - 全部-二肽
        - 全部-三肽
        - 强活性-二肽
        - 强活性-三肽
    """
    all_di = [r for r in stable_rows if 获取短肽类型(r) == "二肽"]
    all_tri = [r for r in stable_rows if 获取短肽类型(r) == "三肽"]
    strong_di = [r for r in all_di if 是否强活性_稳定表(r)]
    strong_tri = [r for r in all_tri if 是否强活性_稳定表(r)]

    result = []
    result.extend(位点频率统计(all_di, residue_col, "全部-二肽"))
    result.extend(位点频率统计(all_tri, residue_col, "全部-三肽"))
    result.extend(位点频率统计(strong_di, residue_col, "强活性-二肽"))
    result.extend(位点频率统计(strong_tri, residue_col, "强活性-三肽"))
    return result


def 端点位点热图数据(频率表: List[Dict[str, Any]]) -> Tuple[List[str], List[str], List[List[float]]]:
    """
    将 N / C 端频率表转换成热图所需矩阵。
    """
    groups = ["全部-二肽", "全部-三肽", "强活性-二肽", "强活性-三肽"]
    matrix = []

    for group in groups:
        mapping = {row["氨基酸"]: 安全转浮点(row["频率"]) or 0.0 for row in 频率表 if row["分析分组"] == group}
        matrix.append([mapping.get(aa, 0.0) for aa in AA_ORDER])

    return groups, AA_ORDER, matrix


# ======================================================================================
# 十一、分析 6：强活性 vs 全体富集（N / C 端）
# ======================================================================================

def 富集统计_单端点(stable_rows: List[Dict[str, str]], residue_col: str) -> List[Dict[str, Any]]:
    """
    对某个端点（N 端或 C 端），计算“强活性 vs 全体”的富集。

    这里先在稳定集上做“整体比较”，不再区分二肽和三肽。
    为什么这样做：
        - 图会更简洁
        - 你能更快看到整体信号
    二肽 / 三肽的更细粒度比较，放到后面的“位置化构效分析”里。

    log2FC 使用轻微平滑：
        freq_s = (count + 0.5) / (total + 0.5 * 20)
        freq_a = (count + 0.5) / (total + 0.5 * 20)
        log2FC = log2(freq_s / freq_a)
    """
    all_rows = list(stable_rows)
    strong_rows = [r for r in stable_rows if 是否强活性_稳定表(r)]

    all_total = 0
    strong_total = 0
    all_counter = Counter()
    strong_counter = Counter()

    for row in all_rows:
        aa = 清洗字符串(row.get(residue_col, ""))
        if aa in AA_ORDER:
            all_counter[aa] += 1
            all_total += 1

    for row in strong_rows:
        aa = 清洗字符串(row.get(residue_col, ""))
        if aa in AA_ORDER:
            strong_counter[aa] += 1
            strong_total += 1

    result = []
    pseudo = 0.5
    denom_all = all_total + pseudo * len(AA_ORDER)
    denom_strong = strong_total + pseudo * len(AA_ORDER)

    for aa in AA_ORDER:
        ca = all_counter.get(aa, 0)
        cs = strong_counter.get(aa, 0)

        fa = ca / all_total if all_total > 0 else 0.0
        fs = cs / strong_total if strong_total > 0 else 0.0

        fa_s = (ca + pseudo) / denom_all if denom_all > 0 else 0.0
        fs_s = (cs + pseudo) / denom_strong if denom_strong > 0 else 0.0

        log2fc = math.log2(fs_s / fa_s) if fa_s > 0 and fs_s > 0 else 0.0

        result.append({
            "氨基酸": aa,
            "全体频数": ca,
            "全体频率": 格式化浮点(fa, 6),
            "强活性频数": cs,
            "强活性频率": 格式化浮点(fs, 6),
            "频率差值_强减全体": 格式化浮点(fs - fa, 6),
            "log2FC_强_vs_全体": 格式化浮点(log2fc, 6),
        })

    result.sort(key=lambda x: abs(安全转浮点(x["log2FC_强_vs_全体"]) or 0), reverse=True)
    return result


# ======================================================================================
# 十二、分析 7：顺序翻转对比
# ======================================================================================

def 生成顺序翻转对比表(stable_rows: List[Dict[str, str]]) -> List[Dict[str, Any]]:
    """
    生成“顺序翻转对比表”。

    思路：
        在稳定序列表中，如果某个序列声明存在反向序列，
        并且它的反向序列本身也在稳定表中，
        那么就把它们组成一对。

    为避免重复：
        使用排序后的 (seq, reverse_seq) 作为唯一键去重。

    这张表非常适合后续做：
        - 顺序敏感性分析
        - 成对排序学习
        - 解释性比较
    """
    seq_map = {}
    for row in stable_rows:
        seq = 清洗字符串(row.get("序列", ""))
        if seq != "":
            seq_map[seq] = row

    seen = set()
    result = []

    for row in stable_rows:
        seq = 清洗字符串(row.get("序列", ""))
        rev = 清洗字符串(row.get("反向序列", ""))
        has_rev = 判定为真(row.get("是否存在反向序列", ""))

        if not seq or not rev or not has_rev:
            continue
        if rev not in seq_map:
            continue
        if seq == rev:
            continue

        key = tuple(sorted([seq, rev]))
        if key in seen:
            continue
        seen.add(key)

        row_a = seq_map[key[0]]
        row_b = seq_map[key[1]]

        a_val = 安全转浮点(row_a.get("精确IC50中位数_uM", ""))
        b_val = 安全转浮点(row_b.get("精确IC50中位数_uM", ""))

        if a_val is None or b_val is None or a_val <= 0 or b_val <= 0:
            fold = None
            log_gap = None
            stronger = "无法判定"
        else:
            fold = max(a_val, b_val) / min(a_val, b_val) if min(a_val, b_val) > 0 else None
            log_gap = abs(math.log10(a_val) - math.log10(b_val))
            if abs(a_val - b_val) < 1e-12:
                stronger = "近似一致"
            elif a_val < b_val:
                stronger = key[0]
            else:
                stronger = key[1]

        both_strong = "是" if (是否强活性_稳定表(row_a) and 是否强活性_稳定表(row_b)) else "否"

        result.append({
            "序列A": key[0],
            "序列B": key[1],
            "短肽类型": row_a.get("短肽类型", ""),
            "序列A中位数IC50_uM": 格式化浮点(a_val, 6),
            "序列B中位数IC50_uM": 格式化浮点(b_val, 6),
            "较强序列": stronger,
            "倍数差": 格式化浮点(fold, 6),
            "绝对log10差值": 格式化浮点(log_gap, 6),
            "序列A活性分层": row_a.get("活性分层建议", ""),
            "序列B活性分层": row_b.get("活性分层建议", ""),
            "是否同为强活性": both_strong,
        })

    result.sort(
        key=lambda x: (
            - (安全转浮点(x["倍数差"]) or -1),
            短肽类型排序值(x["短肽类型"]),
            x["序列A"],
            x["序列B"]
        )
    )
    return result


# ======================================================================================
# 十三、分析 8：位置化构效分析
# ======================================================================================

def 位置频率统计(rows: List[Dict[str, str]], peptide_type: str) -> List[Dict[str, Any]]:
    """
    对指定短肽类型（二肽或三肽）做“位置频率统计”。

    二肽：
        位置1, 位置2
    三肽：
        N端, 中间位, C端

    注意：
        这里直接用“序列”字段拆位点，不依赖别的列。
        这样更稳，也更适合后续你替换输入数据。
    """
    selected = [r for r in rows if 获取短肽类型(r) == peptide_type and 清洗字符串(r.get("序列", "")) != ""]

    if peptide_type == "二肽":
        positions = [("位置1", 0), ("位置2", 1)]
    else:
        positions = [("N端", 0), ("中间位", 1), ("C端", 2)]

    result = []
    for pos_name, idx in positions:
        counter = Counter()
        total = 0
        for row in selected:
            seq = 清洗字符串(row.get("序列", ""))
            if len(seq) > idx:
                aa = seq[idx]
                if aa in AA_ORDER:
                    counter[aa] += 1
                    total += 1
        for aa in AA_ORDER:
            count = counter.get(aa, 0)
            freq = count / total if total > 0 else 0.0
            result.append({
                "短肽类型": peptide_type,
                "位置": pos_name,
                "氨基酸": aa,
                "频数": count,
                "频率": 格式化浮点(freq, 6),
            })

    return result


def 位置频率热图数据(位置频率表: List[Dict[str, Any]], peptide_type: str) -> Tuple[List[str], List[str], List[List[float]]]:
    """
    将位置频率表转成热图矩阵。
    """
    if peptide_type == "二肽":
        positions = ["位置1", "位置2"]
    else:
        positions = ["N端", "中间位", "C端"]

    matrix = []
    for pos in positions:
        mapping = {
            row["氨基酸"]: 安全转浮点(row["频率"]) or 0.0
            for row in 位置频率表
            if row["短肽类型"] == peptide_type and row["位置"] == pos
        }
        matrix.append([mapping.get(aa, 0.0) for aa in AA_ORDER])

    return positions, AA_ORDER, matrix


def 位置富集统计(stable_rows: List[Dict[str, str]], peptide_type: str) -> List[Dict[str, Any]]:
    """
    计算指定短肽类型（只看稳定集）中，
    “强活性 vs 全体”的位置化富集。

    输出列：
        - 位置
        - 氨基酸
        - 全体频数 / 频率
        - 强活性频数 / 频率
        - 差值
        - log2FC

    这是后续做“位置化构效分析”最核心的一张表。
    """
    all_rows = [r for r in stable_rows if 获取短肽类型(r) == peptide_type]
    strong_rows = [r for r in all_rows if 是否强活性_稳定表(r)]

    if peptide_type == "二肽":
        positions = [("位置1", 0), ("位置2", 1)]
    else:
        positions = [("N端", 0), ("中间位", 1), ("C端", 2)]

    result = []
    pseudo = 0.5

    for pos_name, idx in positions:
        all_counter = Counter()
        strong_counter = Counter()
        all_total = 0
        strong_total = 0

        for row in all_rows:
            seq = 清洗字符串(row.get("序列", ""))
            if len(seq) > idx and seq[idx] in AA_ORDER:
                all_counter[seq[idx]] += 1
                all_total += 1

        for row in strong_rows:
            seq = 清洗字符串(row.get("序列", ""))
            if len(seq) > idx and seq[idx] in AA_ORDER:
                strong_counter[seq[idx]] += 1
                strong_total += 1

        denom_all = all_total + pseudo * len(AA_ORDER)
        denom_strong = strong_total + pseudo * len(AA_ORDER)

        for aa in AA_ORDER:
            ca = all_counter.get(aa, 0)
            cs = strong_counter.get(aa, 0)
            fa = ca / all_total if all_total > 0 else 0.0
            fs = cs / strong_total if strong_total > 0 else 0.0
            fa_s = (ca + pseudo) / denom_all if denom_all > 0 else 0.0
            fs_s = (cs + pseudo) / denom_strong if denom_strong > 0 else 0.0
            log2fc = math.log2(fs_s / fa_s) if fa_s > 0 and fs_s > 0 else 0.0

            result.append({
                "短肽类型": peptide_type,
                "位置": pos_name,
                "氨基酸": aa,
                "全体频数": ca,
                "全体频率": 格式化浮点(fa, 6),
                "强活性频数": cs,
                "强活性频率": 格式化浮点(fs, 6),
                "频率差值_强减全体": 格式化浮点(fs - fa, 6),
                "log2FC_强_vs_全体": 格式化浮点(log2fc, 6),
            })

    return result


def 位置富集热图数据(富集表: List[Dict[str, Any]], peptide_type: str) -> Tuple[List[str], List[str], List[List[float]]]:
    """
    将位置富集表转成热图矩阵。
    """
    if peptide_type == "二肽":
        positions = ["位置1", "位置2"]
    else:
        positions = ["N端", "中间位", "C端"]

    matrix = []
    for pos in positions:
        mapping = {
            row["氨基酸"]: 安全转浮点(row["log2FC_强_vs_全体"]) or 0.0
            for row in 富集表
            if row["短肽类型"] == peptide_type and row["位置"] == pos
        }
        matrix.append([mapping.get(aa, 0.0) for aa in AA_ORDER])

    return positions, AA_ORDER, matrix


# ======================================================================================
# 十四、结果索引说明
# ======================================================================================

def 生成结果索引文本() -> str:
    """
    生成一个简单的 Markdown 索引文件，帮助你快速知道每个输出文件是干什么的。
    """
    return """# 高置信短肽图表分析结果索引

本目录下的结果分为两类：

## 一、表格（CSV）
- 长度分布统计表.csv  
  看三个层级（高置信记录 / 唯一序列 / 稳定序列）的肽长分布
- 二肽三肽数量统计表.csv  
  看稳定集内部，不同活性层级下二肽 / 三肽数量
- 精确IC50样本量统计表.csv  
  看三个层级的样本总量、强活性样本量
- 重复序列次数分布表.csv  
  看一个唯一序列对应多少条高置信记录
- 冲突倍数分箱统计表.csv  
  看最大最小倍数的整体分布
- 高冲突序列明细表.csv  
  找出最不稳定、最值得单独复核的序列
- N端位点频率表.csv / C端位点频率表.csv  
  看端点残基在不同分组中的频率
- 强活性_vs_全体_N端富集表.csv / 强活性_vs_全体_C端富集表.csv  
  看强活性集合相对全体的端点富集差异
- 顺序翻转对比表.csv  
  看序列翻转后，活性变化大不大
- 二肽位置频率表.csv / 三肽位置频率表.csv  
  看位置化频率
- 二肽_强活性_vs_全体_位置富集表.csv / 三肽_强活性_vs_全体_位置富集表.csv  
  看位置化富集差异

## 二、图形（SVG）
- 长度分布分组柱状图.svg
- 二肽三肽数量分组柱状图.svg
- 精确IC50样本量分组柱状图.svg
- 重复序列次数分布柱状图.svg
- 冲突倍数分箱柱状图.svg
- N端位点频率热图.svg
- C端位点频率热图.svg
- 强活性_vs_全体_N端富集条形图.svg
- 强活性_vs_全体_C端富集条形图.svg
- 顺序翻转对比散点图.svg
- 二肽位置频率热图.svg
- 三肽位置频率热图.svg
- 二肽_强活性_vs_全体_位置富集热图.svg
- 三肽_强活性_vs_全体_位置富集热图.svg

## 建议的查看顺序
1. 先看：
   - 长度分布
   - 二肽三肽数量
   - 精确IC50样本量
2. 再看：
   - 重复序列次数
   - 冲突倍数
   - 高冲突序列明细
3. 然后重点看：
   - N / C 端位点热图
   - 强活性 vs 全体富集图
   - 顺序翻转对比图
   - 二肽 / 三肽位置化构效热图

## 最适合后续提炼计算机创新点的部分
- 顺序翻转对比表与散点图：适合做 pairwise ranking / 顺序敏感学习
- 高冲突序列明细表：适合做标签不确定性感知学习
- 位置化富集热图：适合做位置感知、可解释模型
"""


# ======================================================================================
# 十五、主程序
# ======================================================================================

def main():
    parser = argparse.ArgumentParser(
        description="基于高置信短肽工作集，自动生成图表分析结果（纯标准库版）。"
    )
    parser.add_argument(
        "--project-root",
        type=str,
        default="",
        help="PepDB 项目根目录，例如 E:/MYS/PepDB。若 PRIORITY_USE_MANUAL_PATHS=True，则此参数会被忽略。"
    )
    args = parser.parse_args()

    script_path = Path(__file__)
    project_root, records_csv, aggregated_csv, stable_csv, output_dir = 解析路径(args.project_root, script_path)

    tables_dir = output_dir / "表格"
    figs_dir = output_dir / "图形"

    print("=" * 100)
    print("高置信短肽工作集图表分析脚本（纯标准库版）")
    print("=" * 100)
    print(f"项目根目录：{project_root}")
    print(f"记录表：{records_csv}")
    print(f"聚合表：{aggregated_csv}")
    print(f"稳定表：{stable_csv}")
    print(f"输出目录：{output_dir}")
    print("-" * 100)

    for p in [records_csv, aggregated_csv, stable_csv]:
        if not p.exists():
            raise FileNotFoundError(f"未找到输入文件：{p}")

    # 读取三张表
    records_rows, records_fields, enc1 = 读取CSV(records_csv)
    aggregated_rows, aggregated_fields, enc2 = 读取CSV(aggregated_csv)
    stable_rows, stable_fields, enc3 = 读取CSV(stable_csv)

    print(f"记录表读取成功：{len(records_rows)} 行，编码 = {enc1}")
    print(f"聚合表读取成功：{len(aggregated_rows)} 行，编码 = {enc2}")
    print(f"稳定表读取成功：{len(stable_rows)} 行，编码 = {enc3}")
    print("-" * 100)

    # 检查必要列
    检查必要列(
        records_rows,
        ["序列", "短肽类型", "肽长", "N端氨基酸", "C端氨基酸", "精确IC50_uM", "是否存在反向序列", "反向序列"],
        "高置信二三肽记录表"
    )
    检查必要列(
        aggregated_rows,
        ["序列", "短肽类型", "肽长", "记录数", "最大最小倍数", "稳定性等级", "精确IC50中位数_uM"],
        "高置信二三肽序列聚合表"
    )
    检查必要列(
        stable_rows,
        ["序列", "短肽类型", "肽长", "N端氨基酸", "C端氨基酸", "精确IC50中位数_uM", "是否存在反向序列", "反向序列"],
        "高置信二三肽稳定序列表"
    )

    tables_dir.mkdir(parents=True, exist_ok=True)
    figs_dir.mkdir(parents=True, exist_ok=True)

    # =========================
    # 1）长度分布
    # =========================
    length_table = 生成长度分布统计表(records_rows, aggregated_rows, stable_rows)
    写出CSV(
        tables_dir / "长度分布统计表.csv",
        ["数据层级", "肽长", "样本数", "占该层级比例"],
        length_table
    )
    length_categories, length_series = 生成长度分布图数据(length_table)
    绘制分组柱状图(
        figs_dir / "长度分布分组柱状图.svg",
        length_categories,
        length_series,
        title="长度分布分组柱状图",
        x_label="肽长",
        y_label="样本数",
        note="同一长度在不同数据层级中的保留情况，可直观看到去重与稳定筛选带来的变化。"
    )

    # =========================
    # 2）二肽 / 三肽数量
    # =========================
    di_tri_table = 生成二肽三肽数量统计表(stable_rows)
    写出CSV(
        tables_dir / "二肽三肽数量统计表.csv",
        ["分析分组", "二肽数量", "二肽占比", "三肽数量", "三肽占比", "总数"],
        di_tri_table
    )
    di_tri_categories, di_tri_series = 生成二肽三肽数量图数据(di_tri_table)
    绘制分组柱状图(
        figs_dir / "二肽三肽数量分组柱状图.svg",
        di_tri_categories,
        di_tri_series,
        title="二肽与三肽数量分组柱状图",
        x_label="分析分组",
        y_label="序列数",
        note="这里聚焦稳定序列层，更适合后续构效分析与第一版建模。"
    )

    # =========================
    # 3）精确 IC50 样本量
    # =========================
    exact_count_table = 生成精确IC50样本量统计表(records_rows, aggregated_rows, stable_rows)
    写出CSV(
        tables_dir / "精确IC50样本量统计表.csv",
        ["数据层级", "全部样本量", "二肽样本量", "三肽样本量", "强活性样本量", "强活性占比"],
        exact_count_table
    )
    exact_categories, exact_series = 生成精确IC50样本量图数据(exact_count_table)
    绘制分组柱状图(
        figs_dir / "精确IC50样本量分组柱状图.svg",
        exact_categories,
        exact_series,
        title="精确 IC50 样本量分组柱状图",
        x_label="数据层级",
        y_label="样本量",
        note="“全部样本量”反映可用监督规模；“强活性样本量”反映高价值样本规模。"
    )

    # =========================
    # 4）重复序列与冲突值
    # =========================
    dup_table = 生成重复序列次数分布表(aggregated_rows)
    写出CSV(
        tables_dir / "重复序列次数分布表.csv",
        ["重复记录数", "唯一序列数", "占全部唯一序列比例"],
        dup_table
    )
    绘制柱状图(
        figs_dir / "重复序列次数分布柱状图.svg",
        labels=[str(r["重复记录数"]) for r in dup_table],
        values=[安全转浮点(r["唯一序列数"]) or 0 for r in dup_table],
        title="重复序列次数分布柱状图",
        x_label="同一序列对应的高置信记录数",
        y_label="唯一序列数",
        rotate_x=False,
        note="横轴越靠右，说明该序列被独立报告的次数越多。"
    )

    conflict_bin_table = 生成冲突倍数分箱统计表(aggregated_rows)
    写出CSV(
        tables_dir / "冲突倍数分箱统计表.csv",
        ["冲突倍数分箱", "序列数", "占全部唯一序列比例"],
        conflict_bin_table
    )
    绘制柱状图(
        figs_dir / "冲突倍数分箱柱状图.svg",
        labels=[r["冲突倍数分箱"] for r in conflict_bin_table],
        values=[安全转浮点(r["序列数"]) or 0 for r in conflict_bin_table],
        title="冲突倍数分箱柱状图",
        x_label="最大最小倍数分箱",
        y_label="序列数",
        rotate_x=True,
        note="最大最小倍数越大，说明同一序列在不同记录中的标签波动越大。"
    )

    high_conflict_table = 生成高冲突序列明细表(aggregated_rows)
    写出CSV(
        tables_dir / "高冲突序列明细表.csv",
        [
            "序列", "短肽类型", "记录数", "精确IC50唯一值个数",
            "精确IC50最小值_uM", "精确IC50最大值_uM", "精确IC50中位数_uM",
            "最大最小倍数", "稳定性等级", "稳定性判定说明",
            "来源种类数", "Assay种类数", "来源列表", "Assay列表"
        ],
        high_conflict_table
    )

    # =========================
    # 5）N / C 端位点频率
    # =========================
    n_table = 生成端点位点频率表(stable_rows, "N端氨基酸")
    c_table = 生成端点位点频率表(stable_rows, "C端氨基酸")

    写出CSV(
        tables_dir / "N端位点频率表.csv",
        ["分析分组", "氨基酸", "频数", "频率"],
        n_table
    )
    写出CSV(
        tables_dir / "C端位点频率表.csv",
        ["分析分组", "氨基酸", "频数", "频率"],
        c_table
    )

    n_rows, n_cols, n_matrix = 端点位点热图数据(n_table)
    c_rows, c_cols, c_matrix = 端点位点热图数据(c_table)

    绘制热图(
        figs_dir / "N端位点频率热图.svg",
        row_labels=n_rows,
        col_labels=n_cols,
        matrix=n_matrix,
        title="N 端位点频率热图",
        diverging=False,
        note="颜色越深，表示该分组中该氨基酸在 N 端出现越频繁。"
    )
    绘制热图(
        figs_dir / "C端位点频率热图.svg",
        row_labels=c_rows,
        col_labels=c_cols,
        matrix=c_matrix,
        title="C 端位点频率热图",
        diverging=False,
        note="颜色越深，表示该分组中该氨基酸在 C 端出现越频繁。"
    )

    # =========================
    # 6）强活性 vs 全体富集（N / C 端）
    # =========================
    n_enrich_table = 富集统计_单端点(stable_rows, "N端氨基酸")
    c_enrich_table = 富集统计_单端点(stable_rows, "C端氨基酸")

    写出CSV(
        tables_dir / "强活性_vs_全体_N端富集表.csv",
        ["氨基酸", "全体频数", "全体频率", "强活性频数", "强活性频率", "频率差值_强减全体", "log2FC_强_vs_全体"],
        n_enrich_table
    )
    写出CSV(
        tables_dir / "强活性_vs_全体_C端富集表.csv",
        ["氨基酸", "全体频数", "全体频率", "强活性频数", "强活性频率", "频率差值_强减全体", "log2FC_强_vs_全体"],
        c_enrich_table
    )

    # 只画 Top 15 绝对变化最大的氨基酸，防止图过长
    top_n = 15
    n_plot_rows = n_enrich_table[:top_n]
    c_plot_rows = c_enrich_table[:top_n]

    绘制发散条形图(
        figs_dir / "强活性_vs_全体_N端富集条形图.svg",
        labels=[r["氨基酸"] for r in n_plot_rows],
        values=[安全转浮点(r["log2FC_强_vs_全体"]) or 0 for r in n_plot_rows],
        title="强活性 vs 全体：N 端富集差异条形图",
        x_label="log2FC (强活性 / 全体)",
        note="红色表示在强活性集合中更富集，蓝色表示相对欠富集。"
    )
    绘制发散条形图(
        figs_dir / "强活性_vs_全体_C端富集条形图.svg",
        labels=[r["氨基酸"] for r in c_plot_rows],
        values=[安全转浮点(r["log2FC_强_vs_全体"]) or 0 for r in c_plot_rows],
        title="强活性 vs 全体：C 端富集差异条形图",
        x_label="log2FC (强活性 / 全体)",
        note="如果某些 C 端残基持续为正值，通常意味着它们更值得进入后续解释性分析。"
    )

    # =========================
    # 7）顺序翻转对比
    # =========================
    reverse_table = 生成顺序翻转对比表(stable_rows)
    写出CSV(
        tables_dir / "顺序翻转对比表.csv",
        ["序列A", "序列B", "短肽类型", "序列A中位数IC50_uM", "序列B中位数IC50_uM", "较强序列", "倍数差", "绝对log10差值", "序列A活性分层", "序列B活性分层", "是否同为强活性"],
        reverse_table
    )
    绘制顺序翻转散点图(
        figs_dir / "顺序翻转对比散点图.svg",
        pair_rows=reverse_table,
        title="顺序翻转对比散点图",
        note="虚线表示两者活性完全相等；偏离越大，说明序列顺序越关键。"
    )

    # =========================
    # 8）位置化构效分析
    # =========================
    di_pos_table = 位置频率统计(stable_rows, "二肽")
    tri_pos_table = 位置频率统计(stable_rows, "三肽")

    写出CSV(
        tables_dir / "二肽位置频率表.csv",
        ["短肽类型", "位置", "氨基酸", "频数", "频率"],
        di_pos_table
    )
    写出CSV(
        tables_dir / "三肽位置频率表.csv",
        ["短肽类型", "位置", "氨基酸", "频数", "频率"],
        tri_pos_table
    )

    di_rows, di_cols, di_matrix = 位置频率热图数据(di_pos_table, "二肽")
    tri_rows, tri_cols, tri_matrix = 位置频率热图数据(tri_pos_table, "三肽")

    绘制热图(
        figs_dir / "二肽位置频率热图.svg",
        row_labels=di_rows,
        col_labels=di_cols,
        matrix=di_matrix,
        title="二肽位置频率热图",
        diverging=False,
        note="用于看二肽两个位置上最常出现的氨基酸类型。"
    )
    绘制热图(
        figs_dir / "三肽位置频率热图.svg",
        row_labels=tri_rows,
        col_labels=tri_cols,
        matrix=tri_matrix,
        title="三肽位置频率热图",
        diverging=False,
        note="用于看三肽 N端 / 中间位 / C端 的整体组成规律。"
    )

    di_enrich_table = 位置富集统计(stable_rows, "二肽")
    tri_enrich_table = 位置富集统计(stable_rows, "三肽")

    写出CSV(
        tables_dir / "二肽_强活性_vs_全体_位置富集表.csv",
        ["短肽类型", "位置", "氨基酸", "全体频数", "全体频率", "强活性频数", "强活性频率", "频率差值_强减全体", "log2FC_强_vs_全体"],
        di_enrich_table
    )
    写出CSV(
        tables_dir / "三肽_强活性_vs_全体_位置富集表.csv",
        ["短肽类型", "位置", "氨基酸", "全体频数", "全体频率", "强活性频数", "强活性频率", "频率差值_强减全体", "log2FC_强_vs_全体"],
        tri_enrich_table
    )

    di_e_rows, di_e_cols, di_e_matrix = 位置富集热图数据(di_enrich_table, "二肽")
    tri_e_rows, tri_e_cols, tri_e_matrix = 位置富集热图数据(tri_enrich_table, "三肽")

    绘制热图(
        figs_dir / "二肽_强活性_vs_全体_位置富集热图.svg",
        row_labels=di_e_rows,
        col_labels=di_e_cols,
        matrix=di_e_matrix,
        title="二肽：强活性 vs 全体 位置富集热图",
        diverging=True,
        note="红色表示在强活性二肽中更富集，蓝色表示相对欠富集。"
    )
    绘制热图(
        figs_dir / "三肽_强活性_vs_全体_位置富集热图.svg",
        row_labels=tri_e_rows,
        col_labels=tri_e_cols,
        matrix=tri_e_matrix,
        title="三肽：强活性 vs 全体 位置富集热图",
        diverging=True,
        note="如果某些 C 端位点持续为红，通常是非常值得深入解释的残基模式。"
    )

    # =========================
    # 9）结果索引文件
    # =========================
    写出文本(output_dir / "结果索引.md", 生成结果索引文本())

    # 控制台摘要
    print("-" * 100)
    print("所有分析结果已输出。")
    print(f"表格目录：{tables_dir}")
    print(f"图形目录：{figs_dir}")
    print(f"索引文件：{output_dir / '结果索引.md'}")
    print("-" * 100)
    print("建议你优先查看：")
    print("1. 图形/强活性_vs_全体_C端富集条形图.svg")
    print("2. 图形/顺序翻转对比散点图.svg")
    print("3. 图形/三肽_强活性_vs_全体_位置富集热图.svg")
    print("4. 表格/高冲突序列明细表.csv")
    print("=" * 100)
    print("分析完成。")


if __name__ == "__main__":
    main()
