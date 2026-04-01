BIOPEP-UWM experimental ACE inhibitor 原始抓取说明
===============================

一、默认输出位置
E:\MYS\PepDB\DB\raw\ace\databases\BIOPEP_UWM\experimental

二、目录说明
1. raw_html/list_pages
   - 列表页原始 HTML
2. raw_html/detail_pages
   - 详情页原始 HTML
3. raw_txt/list_pages
   - 列表页纯文本版本，便于快速查看
4. raw_txt/detail_cards
   - 每条详情页纯文本版本，便于快速核对
5. raw_tables_csv
   - biopep_uwm_experimental_ace_list.csv
     列表页中筛到的 ACE inhibitor 条目索引表
   - biopep_uwm_experimental_ace_detail.csv
     详情页主要字段汇总表
6. meta
   - run_log.txt
   - crawl_summary.json
   - 本说明文件

三、本次抓取摘要
{
  "started_at": "2026-04-01 17:28:35",
  "ended_at": "2026-04-01 17:43:22",
  "output_root": "E:\\MYS\\PepDB\\DB\\raw\\ace\\databases\\BIOPEP_UWM\\experimental",
  "target_activity_name": "ACE inhibitor",
  "target_activity_code": "ah",
  "list_pages_crawled": 374,
  "list_total_pages_seen": 374,
  "list_total_rows_parsed": 0,
  "target_rows_found": 0,
  "detail_success_count": 0,
  "detail_failed_count": 0,
  "list_csv": "E:\\MYS\\PepDB\\DB\\raw\\ace\\databases\\BIOPEP_UWM\\experimental\\raw_tables_csv\\biopep_uwm_experimental_ace_list.csv",
  "detail_csv": "E:\\MYS\\PepDB\\DB\\raw\\ace\\databases\\BIOPEP_UWM\\experimental\\raw_tables_csv\\biopep_uwm_experimental_ace_detail.csv"
}

四、当前脚本目标
- 只做 BIOPEP-UWM experimental 中 ACE inhibitor 原始数据下载
- 不做清洗、标准化、单位换算、去重或多库融合

五、后续扩展建议
- 若以后要抓其他活性，只需修改 TARGET_ACTIVITY_NAME / TARGET_ACTIVITY_CODE
- 若以后要做标准化，可在此脚本产物基础上新建 standardized 脚本，而不是直接改 raw 数据
