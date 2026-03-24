# SCR 目录说明

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
- 初始化时间：2026-03-24 10:58:49
- 初始化脚本：SCR/init_pepdb_structure.py
