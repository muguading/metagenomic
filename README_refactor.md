# Metagenomic Pipeline Refactor Blueprint

这个目录提供了一个适用于宏基因组分析脚本重构的模块化骨架，目标是把单一超长脚本拆解为：

- `main.py` 负责入口与 workflow 调度
- `config/config.yaml` 负责参数配置
- `metagenomic_pipeline/modules/` 负责具体分析步骤
- `metagenomic_pipeline/utils/` 负责公共工具
- `logs/` 负责统一日志输出

## 推荐目录结构

```text
project/
├── main.py
├── config/
│   └── config.yaml
├── logs/
├── metagenomic_pipeline/
│   ├── __init__.py
│   ├── cli.py
│   ├── pipeline.py
│   ├── modules/
│   │   ├── __init__.py
│   │   ├── base.py
│   │   ├── qc.py
│   │   ├── host_filter.py
│   │   ├── taxonomy.py
│   │   ├── assembly.py
│   │   └── report.py
│   └── utils/
│       ├── __init__.py
│       ├── command.py
│       ├── config_loader.py
│       ├── logger.py
│       └── models.py
└── Bac_assemble_260112_newformat.py
```

## 每个模块建议职责

- `qc.py`: 质控、fastp、fastqc、质控结果表格整理
- `host_filter.py`: 去宿主与宿主污染率统计
- `taxonomy.py`: Kraken2 / Bracken / Krona 以及物种丰度汇总
- `assembly.py`: 组装、抛光、组装质量统计
- `report.py`: 汇总各模块输出，生成最终报告
- `utils/command.py`: 统一运行外部命令
- `utils/logger.py`: 统一日志格式与文件输出
- `utils/models.py`: workflow 共享上下文

## 模块调用关系

```text
main.py
  -> cli.py 读取配置
  -> pipeline.py 创建 context
  -> QCModule.run()
  -> HostFilterModule.run()
  -> TaxonomyModule.run()
  -> AssemblyModule.run()
  -> ReportModule.run()
```

上一步的输出通过 `PipelineContext` 传递给下一步，避免全局变量污染。

## 从旧脚本迁移时的建议

1. 先把旧脚本中的公共函数迁移到 `utils/`
2. 再把 `QC_func()`、Kraken2/Bracken 相关函数、组装函数分别迁入对应模块
3. 最后把 `main_process()` 中的分支逻辑收敛到 `pipeline.py`
4. 保持命令行参数名和原有 shell 命令不变，先保证行为一致，再做清理

## 与当前旧脚本的对应关系

从现有 `Bac_assemble_260112_newformat.py` 来看，比较自然的拆分点包括：

- `QC_func()` -> `modules/qc.py`
- `proc_kra()`、`exreadsID()`、Kraken2/Bracken 汇总逻辑 -> `modules/taxonomy.py`
- `assem_vfdr()`、组装主流程 -> `modules/assembly.py`
- `combine_func()`、汇总输出 -> `modules/report.py`
- `run_cmd()`、`copy_pattern()`、`get_logger()`、格式判断函数 -> `utils/`

## main.py 工作流示例

```python
from metagenomic_pipeline.cli import main

if __name__ == "__main__":
    main()
```

## 日志系统示例

```python
logger = setup_logger(
    name="pipeline",
    log_dir=Path("logs"),
    sample_id="sampleA",
)
logger.info("Start taxonomy analysis")
```

## 配置驱动的好处

- 参数集中管理，不再分散在 1000 多行脚本里
- 更容易为不同项目准备不同配置文件
- 更容易添加新模块，例如 `annotation.py`、`amr.py`、`virulence.py`
