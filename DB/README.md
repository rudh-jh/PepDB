# DB 目录说明

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
- 初始化时间：2026-03-24 10:58:49
- 初始化脚本：SCR/init_pepdb_structure.py
