# AHTPDB 原始归档说明

## 目录用途
本目录用于保存 AHTPDB 与 ACE 抑制肽任务相关的原始获取材料。

## 当前阶段
当前处于原始数据归档阶段。
本目录中的文件只代表“原始来源材料”或“原始保存版本”，不代表标准化结果。

## 当前脚本做了什么
1. 保存 AHTPDB 官方页面原始 HTML
2. 保存 AHTPDB 官方下载文本文件
3. 记录获取时间与来源 URL
4. 更新 PepDB 的 meta 清单文件

## 子目录说明
- page_exports/: 页面原始导出内容（HTML 等）
- raw_tables/: 官方下载的原始文本表格
- screenshots/: 页面截图（需后续手动补充）
- retrieval_logs/: 获取过程日志、URL 清单、运行摘要

## 当前归档原则
1. 原始文件不覆盖（除非脚本显式开启覆盖）
2. 原始字段不改写
3. 后续清洗、去重、标准化结果不得回写本目录
4. 本目录内容应可追溯到来源页面或来源说明

## 当前不应直接做什么
- 不应把本目录文件直接当作最终分析输入
- 不应在本目录直接修改值、统一单位、去重
- 不应把手工整理后的结果覆盖原文件

## 来源信息
- task: ace
- source_name: AHTPDB
- source_type: database
- archived_at: 2026-03-24 13:59:20
- archived_by: 当前脚本运行者

## 备注
AHTPDB 当前公开提供的下载入口包括：
- Small Peptides (2-5 residues)
- Long Peptides (6-16 or more residues)
- Peptides with IC50

本目录中的 raw_tables 仅保存其官方下载文本，不代表后续标准化字段结构。
