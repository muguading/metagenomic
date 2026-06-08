群落分析示例数据说明

1. 当前群落分析脚本 `CommunityAnalysis.py` 对输入数据的硬性要求是：
- `--input` 可以是目录或文件；
- `--metadata` 必须是 csv/tsv 元数据文件；
- 元数据第一列会被识别为样本 ID；
- 元数据中必须包含提交任务时填写的分组列，例如 `Group`。

2. 本目录提供了一套可直接用于软件演示和联调的示例：
- `input/`：多样本输入目录，占位为每个样本一个丰度表文件；
- `community_metadata.tsv`：群落分析元数据文件；
- `expected_output/`：使用当前脚本运行后生成的示例输出。

3. 推荐在软件中这样填写：
- 输入路径：`示例数据/群落分析示例/input`
- 元数据文件：`示例数据/群落分析示例/community_metadata.tsv`
- 分组列：`Group`
- 统计层级：`genus`
- 标准化方式：`relative`
- 分析模块：`alpha,beta,lefse,ml`
