---
name: virus-pipeline-extension
description: Extend this metagenomic project with a new virus database, typing workflow, knowledge base, demo task, and report page. Use when adding a new virus genus/species pipeline or when refactoring existing virus report integration.
---

# Virus Pipeline Extension Skill

用于在当前项目里新增一个病毒分析能力，并把它完整接到:

- 参考数据库
- 分型流程
- Demo 数据
- 知识库
- Portal 结果页面

目标不是“能跑一点”，而是把一个新病毒像 `Hepatovirus`、`Orthohantavirus`、`Astroviridae`、`Rotavirus` 一样，接成可演示、可报告、可维护的一整条链。

## 0. 默认工作方式

以后如果要新增一个新病毒，读完这份 `skill.md` 之后，不要直接改代码。

必须先产出一份“接入步骤大纲”，等用户确认后，再开始更新软件里的病毒流程。

这份大纲至少要包含：

1. 病毒基础信息
   - ICTV 对应章节
   - 主参考来源
   - 是否有公认 genotype / subtype / lineage
2. 分型方法
   - broad typing 怎么做
   - subtype typing 怎么做
   - 用全基因组、特定基因还是特定片段
   - 是否需要系统发育、VADR、Nextclade、双位点或多片段
3. 数据库方案
   - typing.xlsx 怎么组织
   - 参考 FASTA / GFF 来源
   - 是否需要子亚型完整基因组库
4. 知识库方案
   - pathogen profile
   - serotype/subtype associations
   - aliases
   - 科研判读摘要里预期展示什么
5. 软件改动范围
   - 主流程文件
   - demo 数据
   - portal 结果页
   - 测试与验收方式

只有在用户确认这份大纲之后，才进入真正的软件更新阶段。

## 1. 先确定病毒类型

先判断你要接入的是哪一类：

1. `Nextclade` 类  
   例如 `RSV / hMPV / DENV / ZIKV / mpox`
2. 参考筛选分型类  
   例如 `Hepatovirus / HPIV / HAdV / Enterovirus / Rhinovirus / Seasonal HCoV`
3. 多层分型类  
   例如 `Hepatovirus (broad -> subtype/genotype)`、`Bandavirus (broad + segments)`、`Orthohantavirus (broad + S segment)`
4. 特殊注释/专位点类  
   例如 `Astroviridae (ORF2 + VADR)`、`Rotavirus (group + G/P)`

先选“最像的现有病毒”，后续全部照它的模式接，不要发明新流派。

### 1.1 如何判断是不是 Nextclade 类

不要靠主观印象判断，优先直接查 `nextclade` 数据集。

推荐做法：

```bash
conda run -n ncov nextclade dataset list
```

或者在已有 `ncov` 环境里直接执行：

```bash
nextclade dataset list
```

判断原则：

1. 如果 `nextclade dataset list` 里已经有目标病毒或其公认数据集  
   则优先归到 `Nextclade` 类
2. 如果没有现成数据集  
   再判断它是不是参考筛选分型类、多层分型类或特殊注释类
3. 不要在 `nextclade` 已有现成库时，重新手搓一套平行分型体系

## 2. 数据库层

### 2.0 来源原则

`database` 里的 `typing.xlsx` 默认优先来自 `ICTV`。

推荐顺序：

1. 先查该病毒在 `ICTV Report` 对应章节
2. 主参考优先取 `Member Species` 表格
3. 如果 `ICTV` 直接给了 genotype / subtype / lineage 划分，就直接按 `ICTV` 组织 `typing.xlsx`
4. 如果 `ICTV` 没直接给子分型参考，但给了分型方法，就按它给的方法去查原始文献
5. 如果 `ICTV` 连子分型方法都没给，再基于文献补充分型策略

例如 `Hepatovirus` 的主参考应来自 ICTV 章节：

- [Genus: Hepatovirus | ICTV](https://ictv.global/report/chapter/picornaviridae/picornaviridae/hepatovirus)

这个页面可作为肝炎病毒分型设计的上游依据，但落地时不要把 `Hepatovirus` 错理解成只有 `HAV`。

当前项目里的肝炎病毒流程要求是：

- 先做 `HAV / HBV / HCV / HDV / HEV` 的 broad typing
- 再进入各自独立参考库做 subtype / genotype 细分
- 报告字段统一写成 `大亚型 / 子亚型`
- `HAV` 可以继续细到 `IA / IB / IIA / IIB / IIIA / IIIB`
- `HBV / HCV / HDV / HEV` 按各自 genotype / subtype 体系输出

不要再把肝炎病毒流程写成 “只支持 HAV 子亚型”。

### 2.1 参考基因组

把病毒参考放到 `database/virus/<VirusName>/` 下，至少准备：

- typing 参考 FASTA
- 对应 GFF/GFF3
- manifest
- 合并 FASTA

现成脚本参考：

- [download_hepatovirus_typing_refs.py](/Users/wuhhh/Desktop/徐老师/代码/metagenomic/scripts/download_hepatovirus_typing_refs.py)
- [download_orthohantavirus_typing_refs.py](/Users/wuhhh/Desktop/徐老师/代码/metagenomic/scripts/download_orthohantavirus_typing_refs.py)
- [download_astroviridae_typing_refs.py](/Users/wuhhh/Desktop/徐老师/代码/metagenomic/scripts/download_astroviridae_typing_refs.py)
- [download_human_rotavirus_complete_genomes.py](/Users/wuhhh/Desktop/徐老师/代码/metagenomic/scripts/download_human_rotavirus_complete_genomes.py)
- [download_rotavirus_a_subtype_complete_genomes.py](/Users/wuhhh/Desktop/徐老师/代码/metagenomic/scripts/download_rotavirus_a_subtype_complete_genomes.py)

如果需要补全长病毒基因组，优先用 `NCBI datasets` 下载，不要先用手工网页下载零散 FASTA。

要求：

- FASTA header 尽量保留 accession 版本号，例如 `.1`
- FASTA contig 名必须和 GFF contig 一致
- 后续会被拷到 `genomes/ref.fa` 供 `snpEff` 或可视化使用，不能只图“能比对”

### 2.2 子分型库

如果病毒有二级或三级分型：

- 像肝炎病毒这样按 broad type 分开建子亚型/基因型库
- 不要把 broad typing 和 subtype typing 混在一套参考里
- 如果 `ICTV` 已直接给出子亚型列表，优先照 `ICTV`
- 如果 `ICTV` 只给出“如何分型”的原则，就按该原则查文献，再决定是用全基因组、特定基因还是特定片段

对 `Hepatovirus` 的具体要求是：

- broad 库单独维护，用来区分 `HAV/HBV/HCV/HDV/HEV`
- subtype/genotype 库按大亚型拆开维护，不共用一套混合参考
- `HAV` 使用 `HAV_subtypes`
- `HBV/HCV/HDV/HEV` 使用各自 `typingB/C/D/E_reference_genomes`

### 2.3 文献补充规则

当 `ICTV` 不能直接给出完整子分型参考时，补文献时按这个顺序：

1. 先确认 `ICTV` 有没有写明分型依据  
   例如：全基因组、VP1、ORF2、某 segment、系统发育阈值等
2. 再查原始或综述文献，确认：
   - 分型位点/区域
   - 命名体系
   - 是否已有公认参考株
3. 最后才下载参考序列构建本地数据库

不要先下载一批序列再倒推“怎么分型”，顺序反了通常会把数据库做乱。

## 3. 主流程接入

核心改动通常在：

- [metagenomic_refactor/virus_analysis.py](/Users/wuhhh/Desktop/徐老师/代码/metagenomic/metagenomic_refactor/virus_analysis.py)
- [metagenomic_refactor/assembly.py](/Users/wuhhh/Desktop/徐老师/代码/metagenomic/metagenomic_refactor/assembly.py)

### 3.1 新增分型函数

最少要补：

- `resolve_*_reference()`
- `run_*_consensus_typing()`
- `_run_*_reference_typing()` 或等价入口
- `_build_serotype_section()` 对应结果汇总分支

### 3.2 输出结构要求

必须区分两类结果：

1. `read-based` 参考筛选
2. `consensus-based` 复判

不要让 consensus 结果覆盖 read-based 结果。  
正确做法参考其他病毒：

- 主目录保留 `selection.tsv` / `coverage.tsv`
- `consensus_typing/` 单独放复判结果

### 3.2.1 BLAST 覆盖度规则

如果 subtype/genotype 识别依赖 `blast`，覆盖度计算不能偷懒按“最佳单条 HSP”处理。

必须遵守：

- 同一参考上的多段 HSP 要先按参考坐标合并区间
- 覆盖碱基数按合并后的非重叠长度计算
- 序列里有 `N` 时，`blast` 被打成分段命中仍要能正确累计覆盖度
- 不能因为分段命中就把真实最佳参考错判成次优参考

这是肝炎病毒子亚型识别里的硬性规则，不是可选优化。

### 3.3 双端输入

如果是二代双端：

- 优先使用真实存在的 `fq1/fq2`
- 不要误把不存在的 `*.final.fastq` 继续往下传

### 3.4 参考拷贝与注释

如果后面要做 `snpEff`：

- 选中参考后复制到 `genomes/`
- FASTA header 和 GFF contig 名必须一致
- 特别检查 accession 版本号是否丢失

## 4. Demo 数据与 Demo 任务

至少准备一套公开 demo 数据，放到 `demo_data/<virus>_demo` 或等价目录。

脚本参考：

- [download_hepatovirus_demo_data.sh](/Users/wuhhh/Desktop/徐老师/代码/metagenomic/scripts/download_hepatovirus_demo_data.sh)
- [download_orthohantavirus_demo_data.sh](/Users/wuhhh/Desktop/徐老师/代码/metagenomic/scripts/download_orthohantavirus_demo_data.sh)
- [download_astroviridae_demo_data.sh](/Users/wuhhh/Desktop/徐老师/代码/metagenomic/scripts/download_astroviridae_demo_data.sh)
- [download_rotavirus_demo_data.sh](/Users/wuhhh/Desktop/徐老师/代码/metagenomic/scripts/download_rotavirus_demo_data.sh)

### 4.1 下载脚本怎么写

下载脚本统一放到 `scripts/download_<virus>_demo_data.sh`。

推荐结构直接照：

- [download_zikav_demo_data.sh](/Users/wuhhh/Desktop/徐老师/代码/metagenomic/scripts/download_zikav_demo_data.sh)
- [download_chikv_demo_data.sh](/Users/wuhhh/Desktop/徐老师/代码/metagenomic/scripts/download_chikv_demo_data.sh)

要求：

1. 脚本开头固定：
   - `#!/usr/bin/env bash`
   - `set -euo pipefail`
2. 统一用：
   - `PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"`
   - `OUT_DIR="${PROJECT_ROOT}/demo_data/<virus>_demo"`
3. `download()` 函数至少要有：
   - 已存在文件时 `skip`
   - `curl -sS -L --fail --retry 3 --retry-delay 2`
4. 如果下载文件体积不小，优先写成 `.part` 临时文件后再 `mv` 到正式文件，避免半截文件冒充下载完成
5. 下载完成后打印：
   - `[ok] <virus> demo data ready in ${OUT_DIR}`

### 4.2 Demo 数据选择规则

不要只要是公开数据就拿来当 demo，必须先确认它适合这个项目的演示流程。

最低要求：

1. 如果是二代数据，优先选真实 `paired-end WGS`
2. 不要选明显过短的 reads 充当标准 demo  
   例如只有 `35 bp / 36 bp` 这类数据，通常不适合当前项目里默认的组装和分型演示
3. 优先选：
   - Illumina MiSeq / NextSeq / NovaSeq 等常见平台
   - 有明确 `WGS` 或等价全基因组策略
   - 文件体积适中、下载成本可接受的数据
4. 不要一味选最大 run  
   demo 数据的目标是“能稳定演示”，不是“把网络和磁盘吃满”
5. 下载前最好先查：
   - `instrument_model`
   - `library_strategy`
   - `base_count`
   - `fastq_ftp`

如果需要从 `ENA` 挑 run，推荐先查元数据再决定，例如：

```bash
curl -sS 'https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=study_accession=%22<STUDY>%22&fields=run_accession,instrument_model,library_strategy,base_count,fastq_ftp&format=tsv&limit=50'
```

如果对读长有怀疑，不要只信项目页面。  
下载后至少抽查第一条 read 的实际长度，再决定这套数据要不要保留为 demo。

### 4.3 `cmd1.sh` 怎么写

在 `demo_data/<virus>_demo/` 下必须放一份 `cmd1.sh`，用于复现实例运行。

Nextclade 类病毒可直接参考：

- [demo_data/zikav_demo/cmd1.sh](/Users/wuhhh/Desktop/徐老师/代码/metagenomic/demo_data/zikav_demo/cmd1.sh)
- [demo_data/chikv_demo/cmd1.sh](/Users/wuhhh/Desktop/徐老师/代码/metagenomic/demo_data/chikv_demo/cmd1.sh)

要求：

1. 入口统一调用：
   - `/Users/wuhhh/Desktop/徐老师/代码/metagenomic/Bac_assemble_260112_newformat.py`
2. `--input` 优先指向样本表 TSV，而不是手写散乱参数
3. 二代 Nextclade 类病毒默认参数风格：
   - `--analysis_target virus`
   - `--inputtype fastq`
   - `--thread 8`
   - `--method freebayes`
   - `--asm_type shortref`
4. `--ref` 和 `--gtf` 要直接指向该病毒的数据集参考：
   - 例如 `database/nextclade_db/<virus>/reference.fasta`
   - 例如 `database/nextclade_db/<virus>/genome_annotation.gff3`
5. `--species` 必须是真实物种名  
   例如 `Chikungunya virus`、`Zika virus`
6. `--runflow` 必须是真正的流程串  
   通常写：
   - `"物种鉴定,基因组组装,分型鉴定"`
7. 不要把 `--species` 和 `--runflow` 写反
8. RNA 病毒通常要确认：
   - `--rna 1`
9. `--genome_len` 要和病毒量级匹配，不要复制别的病毒后忘改

生成完以后至少做一次：

```bash
bash -n demo_data/<virus>_demo/cmd1.sh
```

### 4.4 样本表规则

`demo_data/<virus>_demo/` 里同时要有 `<virus>_samples.tsv`。

要求：

1. 表头沿用现有格式：
   - `样本名称`
   - `三代数据`
   - `二代数据左`
   - `二代数据右`
   - `物种信息`
2. 二代双端必须填真实存在的 `fq1/fq2`
3. 不要让样本表继续引用旧文件名或已经废弃的短 reads
4. 如果更换了 demo run，`download` 脚本、`samples.tsv`、`cmd1.sh` 三处要一起更新

Portal demo 接入位置：

- [bac_analysis_portal/task_manager.py](/Users/wuhhh/Desktop/徐老师/代码/metagenomic/bac_analysis_portal/task_manager.py)
- [bac_analysis_portal/static/app.js](/Users/wuhhh/Desktop/徐老师/代码/metagenomic/bac_analysis_portal/static/app.js)
- [bac_analysis_portal/templates/index.html](/Users/wuhhh/Desktop/徐老师/代码/metagenomic/bac_analysis_portal/templates/index.html)

要补：

- `demo_type`
- 对应权限映射
- 首页 demo 按钮
- demo 输出目录

另外，凡是新增病毒流程，不能只补后端分型，必须同步补“提交入口”。

至少检查：

- 提交任务页 `物种信息` 下拉是否出现该病毒
- 如果该病毒有 broad / subtype 关系，是否需要同时提供“总入口”和“具体物种”
- 监控任务或筛选器里的病毒选择器是否同步出现
- demo 按钮、权限勾选、allowed virus 映射是否同步更新

对肝炎病毒的具体要求：

- 提交页要能选 `Hepatovirus`
- 也要能直接选 `Hepatitis A/B/C/D/E virus`
- 文案不能继续写成只有“甲肝病毒”，否则会误导成 HAV-only 流程

## 5. 知识库接入

知识库文件放在：

- `database/knowledge_base/pathogens/<pathogen>.json`
- `database/knowledge_base/typing/*.json`

HAV 参考：

- [hepatitis_a_virus.json](/Users/wuhhh/Desktop/徐老师/代码/metagenomic/database/knowledge_base/pathogens/hepatitis_a_virus.json)

至少包含：

- `species`
- `common_name`
- `aliases`
- `pathogen_type`
- `clinical_significance`
- `public_health_significance`
- `syndrome_associations`
- `serotype_associations`
- `serotype_panels`
- `interpretation_notes`

关键原则：

1. `aliases` 要覆盖真实现场名称  
   例如 `Hepatitis A virus / Hepatovirus A / Hepatovirus ahepa / 中文名 / 缩写`
2. `serotype_associations` 要写 broad type 和 subtype
3. 不要只写知识库条目，还要验证后端真的能命中

如果是像肝炎病毒这样有“统一大类 + 各自 subtype/genotype”的体系，不能只靠 `pathogens/*.json`。

必须额外做：

- 单独的 `typing` 知识库文件
- broad type 规则
- subtype / genotype 规则
- reference accession、FASTA、GFF 的关联
- 报告里 broad 与 subtype 的双层命中摘要

对 `Hepatovirus`，知识库必须覆盖：

- `HAV/HBV/HCV/HDV/HEV` 五个 broad type
- 各自 subtype / genotype
- broad 与 subtype 的关联命中
- 报告中 `knowledge_summary` 能同时解释大亚型和子亚型

后端命中逻辑主要在：

- [bac_analysis_portal/app.py](/Users/wuhhh/Desktop/徐老师/代码/metagenomic/bac_analysis_portal/app.py)
  - `_build_kb_serotype_index()`
  - `_lookup_kb_serotype_association()`
  - `_build_viral_serotype_knowledge_summary()`
  - `_build_knowledge_interpretation()`

如果病毒现场名称和知识库主名不一致，优先：

- 补 `aliases`
- 在病毒专用分支里做 species 归一化

## 6. 结果页面接入

主要文件：

- [bac_analysis_portal/app.py](/Users/wuhhh/Desktop/徐老师/代码/metagenomic/bac_analysis_portal/app.py)
- [bac_analysis_portal/static/report_runtime.js](/Users/wuhhh/Desktop/徐老师/代码/metagenomic/bac_analysis_portal/static/report_runtime.js)
- [bac_analysis_portal/templates/result_report.html](/Users/wuhhh/Desktop/徐老师/代码/metagenomic/bac_analysis_portal/templates/result_report.html)

### 6.1 后端

在 `_build_serotype_section()` 里新增专用 `mode`，例如：

- `hepatovirus_typing`
- `orthohantavirus_typing`
- `astroviridae_typing`
- `rotavirus_typing`

返回内容至少包括：

- `mode`
- `predicted_group / clade / lineage / subtype`
- `reference_name`
- `summary_cards`
- `mutation_table`
- `knowledge_summary`
- `igv`

### 6.2 前端

在 `report_runtime.js` 里补：

- `is*TypingReport()`
- `extract*ResearchInterpretation()` 或复用通用病毒解析
- `buildResearchScene()` 的专用分支
- 专用 serotype 页面内容

要求：

- 病毒报告不要落回细菌模板
- 如果专用分支失效，也要有“通用病毒兜底模板”
- 页面只展示命中的知识库结果，不要展示内部流程说明

### 6.3 缓存

如果改了结果页结构或 payload：

- 同步提升 `REPORT_CACHE_VERSION`

位置：

- [bac_analysis_portal/app.py](/Users/wuhhh/Desktop/徐老师/代码/metagenomic/bac_analysis_portal/app.py)

## 7. 验证清单

每次新增病毒后，至少做这几步：

1. 参考库检查
   - FASTA / GFF 数量一致
   - header 与 contig 一致
2. 流程检查
   - read-based 与 consensus-based 结果分目录
   - 双端输入不会误用不存在的单端文件
3. 结果页检查
   - `mode` 正确
   - 专用摘要生效
   - 不出现细菌模板词汇
4. 知识库检查
   - `knowledge_summary` 非空
   - 命中 broad type 和 subtype
   - 如果有独立 `typing` 知识库，确认 broad 规则和 subtype 规则都可命中
5. Demo 检查
   - demo 任务能跑
   - demo 结果页能打开
6. 提交入口检查
   - 提交页物种下拉可见
   - 监控/筛选器可见
   - demo 与权限标签不误导

命令建议：

```bash
python3 -m py_compile bac_analysis_portal/app.py
node --check bac_analysis_portal/static/report_runtime.js
node --check bac_analysis_portal/static/app.js
```

必要时补定向测试，例如：

- [tests/test_hepatovirus_reference_selection.py](/Users/wuhhh/Desktop/徐老师/代码/metagenomic/tests/test_hepatovirus_reference_selection.py)
- `blast` 分段 HSP 覆盖度累计测试
- `typing` 知识库 broad/subtype 命中测试

## 8. 最小接入顺序

如果要快速新增一个病毒，推荐严格按这个顺序：

1. 先建 typing 参考库
2. 再接主流程分型
3. 再补 demo 数据
4. 再接 portal `mode`
5. 最后补知识库和科研摘要

不要反过来先做结果页，否则你会一直和空结果、错模板、旧缓存打架。
