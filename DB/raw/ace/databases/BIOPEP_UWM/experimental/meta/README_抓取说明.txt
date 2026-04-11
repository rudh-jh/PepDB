BIOPEP-UWM experimental ACE inhibitor 原始抓取说明
===========================================

一、抓取对象
- 站点：BIOPEP-UWM experimental
- 列表页入口：https://biochemia.uwm.edu.pl/biopep/peptide_data.php
- 目标活动名称：ACE inhibitor
- 目标活动代码：ah
- 详情页抓取对象：peptidedatacard.php?zm_ID=...

二、目录说明
1. raw_html/list_pages
   - 每个列表页的原始 HTML
2. raw_html/detail_pages
   - 每个详情报告页的原始 HTML
3. raw_txt/list_pages
   - 每个列表页转出的便于人工阅读的 TXT
4. raw_txt/detail_cards
   - 每个详情报告页转出的便于人工阅读的 TXT
5. raw_tables_csv
   - biopep_uwm_experimental_ace_list.csv：ACE 条目列表索引表
   - biopep_uwm_experimental_ace_detail.csv：详情报告页字段汇总表
6. meta
   - run_log.txt：运行日志
   - crawl_summary.json：抓取摘要
   - failed_detail_urls.csv：失败详情页记录
   - debug_page1_raw.html / debug_page1_text.txt：首页调试快照

三、当前不做的事情
- 不做数据清洗
- 不做标准化
- 不做单位换算
- 不做跨库去重
- 不做与 AHTPDB 合并

四、推荐查看顺序
1. 先看 raw_tables_csv/biopep_uwm_experimental_ace_detail.csv
2. 再看 raw_tables_csv/biopep_uwm_experimental_ace_list.csv
3. 若发现字段缺失，再看 raw_txt/detail_cards
4. 若仍有疑问，再看 raw_html 原文
