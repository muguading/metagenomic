# Pathogen Workbench

面向公共卫生场景的病原微生物基因组与宏基因组分析工作台。项目把本地/服务器端分析流程、样本库、任务队列、报告浏览、群落分析和病毒分型能力收在一个 Flask + pywebview 应用里，同时保留可直接在命令行运行的批处理入口。

这不是一个“装完 pip 包就能跑完整分析”的轻量 demo。Web 工作台可以快速启动；真实生信分析还需要 conda 环境、外部命令行工具和本地数据库资源。这个边界必须提前说清楚，否则 README 就是在浪费后来者的时间。

## 功能概览

- 病原微生物分析工作台：任务创建、队列管理、日志查看、结果报告、权限与用户管理。
- 细菌/宏基因组流程：质控、宿主处理、Kraken2/Bracken 物种鉴定、组装、注释、耐药/毒力/移动元件分析、MLST/血清型等节点。
- 病毒分析流程：常见呼吸道、肠道、虫媒、自然疫源性、血源性病毒的分型、Nextclade/参考库比对和报告解释。
- 群落分析：支持 abundance 表与 metadata 输入，输出 alpha/beta/LEfSe/机器学习等分析结果。
- 样本与参考库管理：样本库、批量导入预检查、宿主/病原参考库与 cgMLST panel 管理。
- 桌面封装：通过 pywebview 将 Flask 服务封装为 macOS/Windows 桌面应用。

## 项目结构

```text
.
├── bac_analysis_portal/          # Flask Web 工作台、桌面壳、任务调度、前端静态资源
├── metagenomic_refactor/         # 重构后的分析模块：质控、分类、组装、注释、病毒、MAG 等
├── scripts/                      # 数据库构建、demo 数据下载、后处理、打包脚本
├── tests/                        # pytest 回归测试
├── 示例数据/群落分析示例/          # 群落分析演示输入与期望输出
├── Bac_assemble_260112_newformat.py
│                                # 兼容旧参数的主分析入口
├── CommunityAnalysis.py          # 群落分析入口
├── run_bac_analysis_desktop.py   # 桌面应用启动入口
├── run_metagenome_analysis.sh    # 命令行批处理入口
├── requirements-web.txt          # Web/桌面端 Python 依赖
└── requirements-viral-assembly.txt
                                 # 病毒组装节点的软件与数据库说明
```

## 快速启动 Web 工作台

建议使用 Python 3.10+。

```bash
python -m venv .venv_web
./.venv_web/bin/python -m pip install -U pip
./.venv_web/bin/python -m pip install -r requirements-web.txt
./.venv_web/bin/python -m bac_analysis_portal.app
```

启动后访问：

```text
http://127.0.0.1:5055
```

默认演示账号：

```text
admin / admin123
```

如果只是查看界面、任务创建逻辑、报告页和群落分析演示，这一步已经够了。提交真实测序数据前，请先完成下面的生信运行环境配置。

## 启动桌面应用

```bash
./.venv_web/bin/python run_bac_analysis_desktop.py
```

桌面入口会先显示连接页，可选择连接远端服务器，也可以在本机启动服务后自动进入系统。

## 运行宏基因组批处理

命令行入口：

```bash
./run_metagenome_analysis.sh <batch_input.tsv> [output_dir]
```

批处理 TSV 至少包含：

```text
样本名称	三代数据	二代数据左	二代数据右	物种信息
sample_001	0	/path/to/sample_001_R1.fastq.gz	/path/to/sample_001_R2.fastq.gz	nolevel
```

常用环境变量：

```bash
export CONDA_ROOT=/home/hpcdc/miniconda3
export CONDA_ENV=meta_main
export THREADS=10
export RUNFLOW="基因组组装,病毒组装,物种鉴定,元件预测"
export RMHOST=norm
export ABUN=1
```

真实分析依赖的数据库与软件较重。项目代码支持通过环境变量覆盖默认路径：

```bash
export META_DATABASE_ROOT=/path/to/database
export META_KRAKEN_DB=/path/to/kraken2_db
export META_VIRUS_KRAKEN_DB=/path/to/virus_kraken2_db
export META_VIRSORTER2_DB=/path/to/virsorter2_db
export META_CHECKV_DB=/path/to/checkv-db
export META_GENOMAD_DB=/path/to/genomad_db
export META_MOBILEOG_DB=/path/to/mobileOG-db
export META_MOBILEOG_META=/path/to/mobileOG-db-beatrix.csv
```

病毒宏基因组组装节点的详细依赖见 `requirements-viral-assembly.txt`。

## 群落分析示例

仓库提供了一套轻量演示数据：

```text
示例数据/群落分析示例/
├── input/
├── community_metadata.tsv
└── expected_output/
```

在工作台中推荐填写：

```text
输入路径：示例数据/群落分析示例/input
元数据文件：示例数据/群落分析示例/community_metadata.tsv
分组列：Group
统计层级：genus
标准化方式：relative
分析模块：alpha,beta,lefse,ml
```

## 测试

安装测试依赖：

```bash
./.venv_web/bin/python -m pip install -r requirements-dev.txt
```

运行测试：

```bash
./.venv_web/bin/python -m pytest tests
```

当前测试主要覆盖重构模块、流程分支、数据库导入预检查、病毒分型知识库和若干回归问题。外部生信工具链的端到端运行需要在配置完整数据库与 conda 环境后单独验证。

## 打包

macOS：

```bash
bash scripts/build_mac_desktop_app.sh
```

打包并打开：

```bash
bash scripts/build_mac_desktop_app.sh --open
```

Windows 建议在 Windows 10/11 上构建：

```powershell
python -m pip install -r requirements-web.txt
powershell -ExecutionPolicy Bypass -File .\build_windows_desktop_app.ps1
```

如果只生成可执行目录，不生成安装包：

```powershell
powershell -ExecutionPolicy Bypass -File .\build_windows_desktop_app.ps1 -SkipInstaller
```

## 开发约定

- Web 工作台入口在 `bac_analysis_portal/app.py:create_app()`。
- 任务运行由 `bac_analysis_portal/task_runner.py` 负责，它会读取任务 JSON、启动分析命令、写入日志并更新状态。
- 旧主入口 `Bac_assemble_260112_newformat.py` 仍保留参数兼容，具体节点逐步下沉到 `metagenomic_refactor/`。
- 新增分析能力时，优先在 `metagenomic_refactor/` 中补模块和测试，不要继续把逻辑堆回旧主脚本。
- 数据库、大体积中间结果、虚拟环境、打包产物不应进入 Git；它们属于部署资产，不属于源码。
- 仓库边界与资产归属见 `docs/repository-boundaries.md`。新增二进制、数据库、压缩包、运行输出或打包应用前，先确认它是否属于源码。

## 许可证

当前仓库未声明开源许可证。公开发布前请先补充 `LICENSE`，并确认数据库、第三方工具、示例数据和报告模板的再分发权限。
