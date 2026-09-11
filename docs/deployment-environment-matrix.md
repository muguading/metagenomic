# 对外部署环境矩阵

## 结论

- [COMPUTED] 全量离线部署脚本固定要求 **19 个 Linux Conda 环境**；这是当前正式打包脚本的 `REQUIRED_ENVS` 清单，而不是建议数量。
- [COMPUTED] 另需 Ubuntu 系统层运行环境、Portal 状态目录和参考数据库目录；它们不属于 Conda 环境。
- [KNOWN] 所有 Conda 环境必须在 Linux x86_64 上创建或打包；不能将 macOS 环境复制到 Ubuntu。

## 系统层（不计入 19 个环境）

| 层级 | 必需内容 | 用途 |
| --- | --- | --- |
| Ubuntu | Ubuntu 22.04、`systemd`、`bash`、`python3`、`nginx`、`zstd`、`tar`、`sha256sum` | 全量离线安装脚本直接检查或调用。 |
| Conda 运行时 | Miniconda/Miniforge、`conda-pack`（仅构建离线包时） | 创建、打包与解包环境。 |
| Portal 数据目录 | 状态库、任务目录、输出目录 | SQLite 状态、任务 JSON、分析结果。 |
| 参考数据库目录 | Kraken2、病毒 Kraken2、VirSorter2、CheckV、geNomad、MobileOG 和 Nextclade 等资产 | 完整分析时由 `META_*` 变量定位。 |

## 19 个 Conda 环境

| # | 环境名 | 主要软件 / 责任范围 | 证据状态 |
| ---: | --- | --- | --- |
| 1 | `web_runtime` | Python 3.10、Gunicorn，以及 `requirements-release.lock` 中的 Flask、SQLAlchemy、pywebview、pandas、scikit-learn、xgboost、lightgbm、PyInstaller 等。 | 版本锁文件完整；已在测试服务器完成安装与导入验证。 |
| 2 | `meta_main` | Python 3.10、pandas、Biopython、matplotlib、openpyxl、perl、ruby、OpenJDK；fastp/fastqc/seqkit/pigz；Kraken2/Bracken/Krona/taxonkit；BLAST/DIAMOND/bwa/samtools/minimap2/mosdepth/bedtools/bcftools/freebayes；SPAdes/Flye/Unicycler/miniasm/gfatools/canu/raven；Prokka/ABRicate/mlst/minced/MAFFT/Nextclade。 | `envs/main_pipeline_env.yml` 已提供可创建清单。 |
| 3 | `genomad_aux` | VirSorter2、CheckV、geNomad、chewBBACA。 | 运行时映射定义了职责；需在发布前冻结 YAML 与版本。 |
| 4 | `amr_aux` | hamronize、ResFinder、AMRFinder、RGI，以及 pmga/salty/lissero 等耐药与分型补充工具。 | 运行时映射定义了职责；需在发布前冻结 YAML 与版本。 |
| 5 | `host_filter` | hostile、kneaddata。 | 运行时映射定义了职责；需在发布前冻结 YAML 与版本。 |
| 6 | `mag_aux` | megahit、BASALT、CoverM、GTDB-Tk；分箱基础软件还包括 bwa、samtools、MetaBAT2、SemiBin2、VAMB、DAS Tool。 | `envs/mag_binning_env.yml` 覆盖后半部分；其余软件需冻结 YAML。 |
| 7 | `cm210` | CheckM2、ectyper。 | 运行时映射定义了职责；需在发布前冻结 YAML 与版本。 |
| 8 | `sistr_hicap` | SISTR、HICAP。 | 运行时映射定义了职责；需在发布前冻结 YAML 与版本。 |
| 9 | `qiime2` | QIIME 2 的 alpha/beta 多样性、样本分类器与群落统计组件。 | 运行时映射定义了职责；需在发布前冻结 YAML 与版本。 |
| 10 | `microeco` | R、microeco、LEfSe 与生态统计/可视化依赖。 | 运行时映射定义了职责；需在发布前冻结 YAML 与版本。 |
| 11 | `genovi` | genovi、COG 图谱与统计依赖。 | 运行时映射定义了职责；需在发布前冻结 YAML 与版本。 |
| 12 | `plasflow` | PlasFlow。 | 运行时映射定义了职责；需在发布前冻结 YAML 与版本。 |
| 13 | `longread_aux` | medaka、Clair3。 | 运行时映射定义了职责；需在发布前冻结 YAML 与版本。 |
| 14 | `chewie` | chewBBACA 相关工具。 | 运行时映射定义了职责；需在发布前冻结 YAML 与版本。 |
| 15 | `ncov` | TB-Profiler、Nextclade 和病毒分型所需运行组件。 | 运行时映射与 Nextclade 查找逻辑共同表明此职责；需在发布前冻结 YAML 与版本。 |
| 16 | `choleraefinder` | CholeraeFinder。 | 运行时映射定义了职责；需在发布前冻结 YAML 与版本。 |
| 17 | `report_env` | R 与报告生成/结果整合依赖。 | 运行时映射定义了职责；需在发布前冻结 YAML 与版本。 |
| 18 | `Vlib` | 旧版病毒流程兼容依赖。 | 运行时映射定义了职责；需在发布前冻结 YAML 与版本。 |
| 19 | `tNGS` | 旧版病毒流程补充依赖。 | 运行时映射定义了职责；需在发布前冻结 YAML 与版本。 |

## 发布前必须补齐的资产

- [COMPUTED] 生产环境变量必须配置 `PORTAL_MODE`、`PORTAL_SECRET_KEY`、`PORTAL_INITIAL_ADMIN_PASSWORD`、`PORTAL_DB_PATH`、`BAC_ANALYSIS_TASK_ROOT`。
- [COMPUTED] 完整分析还必须配置并验证 `META_DATABASE_ROOT`、`META_KRAKEN_DB`、`META_VIRUS_KRAKEN_DB`、`META_VIRSORTER2_DB`、`META_CHECKV_DB`、`META_GENOMAD_DB`、`META_MOBILEOG_DB`、`META_MOBILEOG_META`。
- [INFERRED] 特殊 Conda 环境目前缺少逐环境的版本化 YAML，因此在生成正式离线包前应先补齐并锁定每个环境的 YAML/显式包清单；否则无法可靠复现相同软件版本。置信度：HIGH。

## 测试服务器验证（2026-07-28）

- [COMPUTED] 服务器为 Ubuntu 22.04.5 x86_64、4 vCPU、3.6 GiB 内存、无 swap。
- [COMPUTED] 初始可用磁盘为 15 GiB；安装 Miniconda 与已验证的 Web 环境后，可用磁盘为 12 GiB。
- [COMPUTED] 项目的全量打包脚本声明 80 GiB 最低可用空间，因此该服务器不能用于构建全量数据库或完整离线分析包。
- [COMPUTED] 已安装 `/home/ubuntu/miniconda3`，并创建 `pathogen-web-validate`；锁定依赖安装成功、`pip check` 成功，Flask、SQLAlchemy、pywebview、scikit-learn、xgboost、lightgbm 导入成功。
- [COMPUTED] TUNA 的 Anaconda 主仓库元数据下载约 1.16 秒，官方仓库约 10.73 秒；已写入用户级 `/home/ubuntu/.condarc` 并启用严格渠道优先级、libmamba 和 TUNA 镜像。

## 发布验证命令

```bash
/opt/pathogen-workbench/conda/envs/web_runtime/bin/python \
  /opt/pathogen-workbench/app/scripts/check_deployment.py \
  --project-root /opt/pathogen-workbench/app \
  --env-file /opt/pathogen-workbench/app/deployment/portal.env \
  --profile full-analysis
```

该命令只有在 19 个环境、关键二进制、生产变量和数据库路径均存在时才应通过。
