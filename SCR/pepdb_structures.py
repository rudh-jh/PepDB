# -*- coding: utf-8 -*-
"""
脚本名称：
    init_pepdb_structure.py

脚本作用：
    初始化 PepDB 项目的第一阶段骨架。
    该脚本主要完成以下事情：

    1. 创建项目所需的目录结构
    2. 创建基础 README / Markdown 说明文件
    3. 创建 meta 层的空白 CSV 清单文件
    4. 创建 SCR/src 下的基础 Python 包结构
    5. 所有操作默认“只创建，不覆盖”

为什么现在先做这个：
    你的项目当前重点不是立刻建模，而是先把：
    - 原始数据归档体系
    - 标准化层
    - 工作集层
    - 分析层
    - 元数据说明层
    的骨架搭起来。

    所以这个脚本属于：
    阶段 1：项目骨架与说明文档

输入：
    无需外部输入文件。
    默认使用固定根目录：E:\\MYS\\PepDB

输出：
    在 E:\\MYS\\PepDB 下创建或补齐目录与基础模板文件。

重要原则：
    1. 不覆盖已存在文件
    2. 不删除任何现有目录
    3. Windows 路径优先
    4. 便于以后扩展到 DPPIV 等任务
"""

from pathlib import Path
from datetime import datetime
import csv


# =========================
# 1. 基础路径设置
# =========================
# 这里明确整个项目的根目录。
# 如果以后项目根目录变了，只需要修改这里。
ROOT_DIR = Path(r"E:\MYS\PepDB")

# 数据目录
DB_DIR = ROOT_DIR / "DB"

# Python 工程目录
SCR_DIR = ROOT_DIR / "SCR"

# 说明文档目录（项目级）
DOCS_DIR = ROOT_DIR / "docs"

# 当前时间字符串，用于写入模板说明
NOW_STR = datetime.now().strftime("%Y-%m-%d %H:%M:%S")


# =========================
# 2. 目录清单
# =========================
# 说明：
# 这里列出需要创建的目录。
# 如果目录已经存在，脚本会直接跳过，不会报错。
DIRS_TO_CREATE = [
    # ---------- 项目级文档目录 ----------
    DOCS_DIR,

    # ---------- DB 主结构 ----------
    DB_DIR,
    DB_DIR / "raw",
    DB_DIR / "standardized",
    DB_DIR / "worksets",
    DB_DIR / "analysis",
    DB_DIR / "meta",

    # ---------- raw / ace ----------
    DB_DIR / "raw" / "ace",
    DB_DIR / "raw" / "ace" / "databases",
    DB_DIR / "raw" / "ace" / "databases" / "AHTPDB",
    DB_DIR / "raw" / "ace" / "databases" / "BIOPEP_UWM",
    DB_DIR / "raw" / "ace" / "databases" / "FermFooDb",
    DB_DIR / "raw" / "ace" / "databases" / "DFBP",
    DB_DIR / "raw" / "ace" / "literature",
    DB_DIR / "raw" / "ace" / "literature" / "search_records",
    DB_DIR / "raw" / "ace" / "literature" / "pdfs",
    DB_DIR / "raw" / "ace" / "literature" / "supplementary",
    DB_DIR / "raw" / "ace" / "literature" / "manual_extraction_raw",

    # ---------- standardized / ace ----------
    DB_DIR / "standardized" / "ace",

    # ---------- worksets / ace ----------
    DB_DIR / "worksets" / "ace",

    # ---------- analysis / ace ----------
    DB_DIR / "analysis" / "ace",

    # ---------- SCR 主结构 ----------
    SCR_DIR,
    SCR_DIR / "src",
    SCR_DIR / "scripts",
    SCR_DIR / "configs",
    SCR_DIR / "docs",
    SCR_DIR / "logs",
    SCR_DIR / "tests",

    # ---------- Python 包结构 ----------
    SCR_DIR / "src" / "pepdb",
    SCR_DIR / "src" / "pepdb" / "io",
    SCR_DIR / "src" / "pepdb" / "parsers",
    SCR_DIR / "src" / "pepdb" / "standardizers",
    SCR_DIR / "src" / "pepdb" / "worksets",
    SCR_DIR / "src" / "pepdb" / "analysis",
    SCR_DIR / "src" / "pepdb" / "reporting",
    SCR_DIR / "src" / "pepdb" / "utils",
]


# =========================
# 3. 文件模板内容
# =========================
# 说明：
# 这里为重要文件准备模板内容。
# 所有文件都只在“尚不存在时”创建。
ROOT_README = f"""# PepDB

## 项目定位
PepDB 是一个围绕生物活性肽数据整理与分析构建的科研数据工程项目。

当前主线任务：
- ACE inhibitory peptides（ACE 抑制肽）

当前研究重点：
- 二肽
- 三肽
- 尤其关注带 IC50 的子集

## 当前优先级
1. 原始数据收集与归档
2. 基础统计分析
3. 标准化与清洗
4. 构建高质量工作集
5. 建模

## 目录说明
- `DB/`：数据目录
- `SCR/`：Python 项目目录
- `docs/`：项目级架构与说明文档

## 重要原则
1. 原始数据永不覆盖
2. 后续整理结果另存
3. 重要文件必须可追溯
4. 当前不急于建模
5. 当前主线聚焦 ACE，后续保留扩展到 DPPIV 的接口

## 初始化信息
- 初始化时间：{NOW_STR}
- 初始化脚本：SCR/init_pepdb_structure.py
"""

DB_README = f"""# DB 目录说明

本目录用于存放 PepDB 项目的全部数据资产。

## 子目录说明

### raw/
保存原始数据，永不覆盖，永不直接修改。
包括但不限于：
- 数据库原始下载文件
- 网页快照
- 文献检索记录
- PDF
- supplementary
- 手工提取表格的原始版本

### standardized/
保存字段统一后的中间表。
目标是先统一字段结构，而不是直接追求最终清洗结果。

### worksets/
保存研究直接使用的工作集，例如：
- 二肽子集
- 三肽子集
- 二三肽合并子集
- 二三肽且带 IC50 的子集
- 高置信子集

### analysis/
保存分析结果，包括：
- 统计表
- 图片
- Markdown 报告
- 中间分析产物

### meta/
保存元数据与说明文件，包括：
- source_manifest.csv
- file_catalog.csv
- field_dictionary.md
- project_decisions.md

## 重要原则
1. 原始层 raw 永不覆盖
2. 任何处理结果都应另存
3. 每个重要文件都应能追溯来源、时间、脚本和用途

## 初始化信息
- 初始化时间：{NOW_STR}
- 初始化脚本：SCR/init_pepdb_structure.py
"""

SCR_README = f"""# SCR 目录说明

本目录用于存放 PepDB 项目的 Python 工程代码。

## 子目录说明
- `src/`：核心包代码
- `scripts/`：任务入口脚本
- `configs/`：配置文件
- `docs/`：代码相关文档
- `logs/`：运行日志
- `tests/`：测试代码
- `envs/`：虚拟环境目录（当前放在 SCR 下）

## 代码组织原则
1. 尽量模块化
2. Windows 兼容优先
3. 同时兼容命令行和 PyCharm
4. 不把所有逻辑堆到单一巨型脚本里
5. 便于扩展更多数据库和更多任务类型

## 当前包结构
`src/pepdb/` 下预留以下模块：
- io
- parsers
- standardizers
- worksets
- analysis
- reporting
- utils

## 初始化信息
- 初始化时间：{NOW_STR}
- 初始化脚本：SCR/init_pepdb_structure.py
"""

ARCHITECTURE_MD = f"""# PepDB 架构说明

## 当前主线
当前只聚焦 ACE inhibitory peptides。

## 总体数据流
原始数据层 -> 标准化层 -> 工作集层 -> 分析层 -> 建模层

## 分层说明

### 1. 原始数据层（raw）
目标：
- 完整保存数据库和文献原始证据
- 不覆盖，不直接修改

### 2. 标准化层（standardized）
目标：
- 统一字段
- 保留原始字段
- 记录解析状态

### 3. 工作集层（worksets）
目标：
- 根据研究问题构建稳定数据子集

### 4. 分析层（analysis）
目标：
- 自动输出图、表、Markdown 报告
- 优先保证表格完整性

### 5. 建模层（未来）
目标：
- 回归优先
- 排序次之

## 扩展原则
虽然当前主线是 ACE，但目录与脚本设计应保留扩展到 DPPIV 等任务的能力。

## 初始化信息
- 初始化时间：{NOW_STR}
- 初始化脚本：SCR/init_pepdb_structure.py
"""

PROJECT_DECISIONS_MD = f"""# 项目决策记录

> 本文件用于记录 PepDB 项目中的关键结构决策，避免后续反复改口径。

## 当前已确认决策

### 决策 1：当前主线只做 ACE 抑制肽
- 状态：已确认
- 说明：未来可扩展到 DPPIV，但当前不混做多个任务。

### 决策 2：当前阶段不以建模为优先
- 状态：已确认
- 当前优先级：
  1. 原始数据收集与归档
  2. 基础统计分析
  3. 标准化与清洗
  4. 工作集构建
  5. 建模

### 决策 3：主库不限长度，但重点关注二肽/三肽
- 状态：已确认
- 说明：后续重点工作集将突出二肽、三肽和带 IC50 的子集。

### 决策 4：原始数据永不覆盖
- 状态：已确认
- 说明：raw 层只保存原始证据，不做覆盖式修改。

### 决策 5：envs 放在 SCR 下
- 状态：已确认
- 路径：E:\\MYS\\PepDB\\SCR\\envs

## 后续可继续补充
建议以后每新增一个重要规则，就在此文件追加记录：
- 决策背景
- 采取方案
- 为什么这样做
- 对哪些脚本/目录产生影响

## 初始化信息
- 初始化时间：{NOW_STR}
- 初始化脚本：SCR/init_pepdb_structure.py
"""

FIELD_DICTIONARY_MD = f"""# 字段字典（初稿）

> 本文件用于记录标准化层和工作集层计划采用的关键字段。
> 当前是初稿，后续会逐步补充。

## 基础字段建议

| 字段名 | 含义 | 备注 |
|---|---|---|
| record_id | 记录唯一标识 | 标准化后生成 |
| task | 任务名称 | 当前主线为 ace |
| sequence | 肽序列 | 建议保留原始大小写信息来源 |
| peptide_length | 肽长度 | 可自动计算 |
| source_type | 来源类型 | 数据库 / 文献 / supplementary / 手工提取 |
| source_name | 来源名称 | 如 AHTPDB、BIOPEP-UWM、某论文 |
| source_record_id | 来源内部编号 | 若原始来源提供 |
| evidence_type | 证据类型 | 实验验证 / 预测 / 对接 / 其他 |
| target | 作用靶点 | 当前主线一般为 ACE |
| activity_label_raw | 原始活性描述 | 保留原始文本 |
| ic50_raw | 原始 IC50 文本 | 不覆盖 |
| ic50_value | 解析出的数值部分 | 自动解析字段 |
| ic50_unit | 原始或标准单位 | 如 uM / mg/mL |
| ic50_relation | 比较关系 | =, <, <=, > 等 |
| ic50_uM | 统一换算后的 uM 数值 | 若可换算 |
| ic50_parse_status | IC50 解析状态 | success / failed / partial |
| source_file | 来源文件路径 | 相对路径优先 |
| notes | 备注 | 用于补充解释 |

## 说明
1. 原始字段必须保留
2. 自动解析字段要明确标注
3. 统一后的值不代表最终生物学结论，只代表当前工程处理结果

## 初始化信息
- 初始化时间：{NOW_STR}
- 初始化脚本：SCR/init_pepdb_structure.py
"""

RAW_ACE_README = f"""# raw/ace 目录说明

本目录用于保存 ACE 抑制肽任务相关的原始数据。

## 子目录建议
- `databases/`：公共数据库原始下载内容
- `literature/search_records/`：检索式、检索结果、筛选记录
- `literature/pdfs/`：文献 PDF
- `literature/supplementary/`：补充材料
- `literature/manual_extraction_raw/`：手工提取表格的原始版本

## 原则
1. 不直接修改原始文件
2. 尽量保留下载时间、来源链接、来源说明
3. 后续清洗与整理结果不要回写到本层

## 初始化信息
- 初始化时间：{NOW_STR}
- 初始化脚本：SCR/init_pepdb_structure.py
"""

STANDARDIZED_ACE_README = f"""# standardized/ace 目录说明

本目录用于存放 ACE 抑制肽数据的标准化中间表。

## 本层目标
- 统一字段
- 保留原始值
- 记录解析状态
- 为后续工作集构建做准备

## 注意
本层不是最终清洁真值层，而是工程化统一字段层。

## 初始化信息
- 初始化时间：{NOW_STR}
- 初始化脚本：SCR/init_pepdb_structure.py
"""

WORKSETS_ACE_README = f"""# worksets/ace 目录说明

本目录用于存放研究直接使用的 ACE 工作集。

## 未来典型工作集
- 二肽子集
- 三肽子集
- 二三肽合并子集
- 二三肽且带 IC50 的子集
- 高置信子集

## 初始化信息
- 初始化时间：{NOW_STR}
- 初始化脚本：SCR/init_pepdb_structure.py
"""

ANALYSIS_ACE_README = f"""# analysis/ace 目录说明

本目录用于保存 ACE 抑制肽相关分析结果。

## 未来输出类型
- 统计表
- 图片
- Markdown 报告
- 中间分析结果

## 输出原则
1. 表格完整性优先
2. 图尽量生成
3. 同步生成 Markdown 报告
4. 每个重要输出目录尽量带 README

## 初始化信息
- 初始化时间：{NOW_STR}
- 初始化脚本：SCR/init_pepdb_structure.py
"""

CONFIG_README = f"""# configs 目录说明

本目录用于存放配置文件，例如：
- 路径配置
- 字段映射
- 阈值参数
- 分析参数

建议后续将“强活性阈值”等设计为可配置项，而不是硬编码在脚本中。

## 初始化信息
- 初始化时间：{NOW_STR}
"""

SCR_DOCS_README = f"""# SCR/docs 目录说明

本目录用于存放代码工程相关说明文档，例如：
- 模块说明
- 处理流程说明
- 脚本使用说明
- 开发记录

## 初始化信息
- 初始化时间：{NOW_STR}
"""

LOGS_README = f"""# logs 目录说明

本目录用于保存脚本运行日志。

建议后续：
- 每个批处理脚本单独输出日志
- 日志文件名带日期时间
- 避免把重要结果只放日志里而不落盘成表

## 初始化信息
- 初始化时间：{NOW_STR}
"""

TESTS_README = f"""# tests 目录说明

本目录预留给后续测试代码使用。
当前阶段不是最优先，但建议逐步为关键模块补基本测试。

## 初始化信息
- 初始化时间：{NOW_STR}
"""

SCRIPTS_README = f"""# scripts 目录说明

本目录用于存放任务入口脚本。

未来建议在这里放置：
- 初始化项目骨架
- 下载数据库原始数据
- 构建标准化表
- 构建工作集
- 运行分析
- 生成报告

## 初始化信息
- 初始化时间：{NOW_STR}
"""


# =========================
# 4. 文件路径清单
# =========================
FILES_TO_CREATE = {
    ROOT_DIR / "README.md": ROOT_README,
    DB_DIR / "README.md": DB_README,
    SCR_DIR / "README.md": SCR_README,
    DOCS_DIR / "architecture.md": ARCHITECTURE_MD,
    DB_DIR / "meta" / "project_decisions.md": PROJECT_DECISIONS_MD,
    DB_DIR / "meta" / "field_dictionary.md": FIELD_DICTIONARY_MD,
    DB_DIR / "raw" / "ace" / "README.md": RAW_ACE_README,
    DB_DIR / "standardized" / "ace" / "README.md": STANDARDIZED_ACE_README,
    DB_DIR / "worksets" / "ace" / "README.md": WORKSETS_ACE_README,
    DB_DIR / "analysis" / "ace" / "README.md": ANALYSIS_ACE_README,
    SCR_DIR / "configs" / "README.md": CONFIG_README,
    SCR_DIR / "docs" / "README.md": SCR_DOCS_README,
    SCR_DIR / "logs" / "README.md": LOGS_README,
    SCR_DIR / "tests" / "README.md": TESTS_README,
    SCR_DIR / "scripts" / "README.md": SCRIPTS_README,
}


# =========================
# 5. Python 包初始化文件
# =========================
INIT_FILES = {
    SCR_DIR / "src" / "pepdb" / "__init__.py": '"""PepDB 主包。"""\n',
    SCR_DIR / "src" / "pepdb" / "io" / "__init__.py": '"""I/O 相关模块。"""\n',
    SCR_DIR / "src" / "pepdb" / "parsers" / "__init__.py": '"""各数据源解析器模块。"""\n',
    SCR_DIR / "src" / "pepdb" / "standardizers" / "__init__.py": '"""标准化处理模块。"""\n',
    SCR_DIR / "src" / "pepdb" / "worksets" / "__init__.py": '"""工作集构建模块。"""\n',
    SCR_DIR / "src" / "pepdb" / "analysis" / "__init__.py": '"""统计分析模块。"""\n',
    SCR_DIR / "src" / "pepdb" / "reporting" / "__init__.py": '"""报告生成模块。"""\n',
    SCR_DIR / "src" / "pepdb" / "utils" / "__init__.py": '"""通用工具函数模块。"""\n',
}


# =========================
# 6. CSV 清单模板
# =========================
CSV_FILES = {
    DB_DIR / "meta" / "source_manifest.csv": [
        "source_id",
        "task",
        "source_type",
        "source_name",
        "source_url_or_reference",
        "evidence_level",
        "access_date",
        "raw_storage_path",
        "notes",
    ],
    DB_DIR / "meta" / "file_catalog.csv": [
        "file_id",
        "source_id",
        "task",
        "relative_path",
        "file_type",
        "content_role",
        "created_or_downloaded_at",
        "generated_by",
        "version_tag",
        "notes",
    ],
}


# =========================
# 7. 工具函数
# =========================
def ensure_dir(path: Path) -> None:
    """
    创建目录（若不存在）。

    参数：
        path: 需要创建的目录路径

    说明：
        parents=True 表示如果上级目录不存在，也一并创建。
        exist_ok=True 表示目录已存在时不报错。
    """
    path.mkdir(parents=True, exist_ok=True)


def write_text_if_missing(path: Path, content: str) -> bool:
    """
    若目标文件不存在，则写入文本内容。

    参数：
        path: 文件路径
        content: 要写入的文本内容

    返回：
        True  -> 本次新建了文件
        False -> 文件已存在，因此未覆盖
    """
    if path.exists():
        return False

    ensure_dir(path.parent)
    path.write_text(content, encoding="utf-8")
    return True


def write_csv_if_missing(path: Path, headers: list) -> bool:
    """
    若目标 CSV 文件不存在，则写入表头。

    参数：
        path: CSV 文件路径
        headers: 表头列表

    返回：
        True  -> 本次新建了文件
        False -> 文件已存在，因此未覆盖
    """
    if path.exists():
        return False

    ensure_dir(path.parent)
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerow(headers)
    return True


# =========================
# 8. 主执行逻辑
# =========================
def main():
    """
    主函数。

    执行顺序：
    1. 创建目录
    2. 创建 Markdown / README 文件
    3. 创建 __init__.py 文件
    4. 创建 meta CSV 文件
    5. 打印执行结果摘要
    """
    print("=" * 70)
    print("PepDB 项目骨架初始化开始")
    print(f"根目录：{ROOT_DIR}")
    print("=" * 70)

    created_dirs = []
    created_files = []
    skipped_files = []

    # 1) 创建目录
    for d in DIRS_TO_CREATE:
        existed_before = d.exists()
        ensure_dir(d)
        if not existed_before:
            created_dirs.append(d)

    # 2) 创建说明文件
    for file_path, content in FILES_TO_CREATE.items():
        if write_text_if_missing(file_path, content):
            created_files.append(file_path)
        else:
            skipped_files.append(file_path)

    # 3) 创建 __init__.py
    for file_path, content in INIT_FILES.items():
        if write_text_if_missing(file_path, content):
            created_files.append(file_path)
        else:
            skipped_files.append(file_path)

    # 4) 创建 CSV 模板
    for file_path, headers in CSV_FILES.items():
        if write_csv_if_missing(file_path, headers):
            created_files.append(file_path)
        else:
            skipped_files.append(file_path)

    # 5) 输出结果摘要
    print("\n[一] 新创建的目录：")
    if created_dirs:
        for item in created_dirs:
            print(f"  [DIR]  {item}")
    else:
        print("  没有新增目录（可能之前都已存在）。")

    print("\n[二] 新创建的文件：")
    if created_files:
        for item in created_files:
            print(f"  [FILE] {item}")
    else:
        print("  没有新增文件（可能之前都已存在）。")

    print("\n[三] 已存在、因此未覆盖的文件：")
    if skipped_files:
        for item in skipped_files:
            print(f"  [SKIP] {item}")
    else:
        print("  没有需要跳过的文件。")

    print("\n" + "=" * 70)
    print("PepDB 项目骨架初始化完成")
    print("=" * 70)
    print("建议下一步：")
    print("1. 先检查 DB/meta 下生成的说明文件和 CSV 表头是否符合你的想法")
    print("2. 然后开始进入阶段 2：原始数据归档系统")
    print("3. 优先从 AHTPDB / BIOPEP-UWM / 文献检索记录的 raw 层归档开始")
    print("=" * 70)


if __name__ == "__main__":
    main()