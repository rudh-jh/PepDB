# 字段字典（初稿）

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
- 初始化时间：2026-03-24 10:58:49
- 初始化脚本：SCR/init_pepdb_structure.py
