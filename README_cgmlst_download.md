# cgMLST.org Panel 下载工具

这个目录包含多个用于从 cgMLST.org 下载全基因组 MLST 方案的脚本。

## 文件说明

| 文件 | 说明 | 推荐度 |
|------|------|--------|
| `test_cgmlst_api.py` | API 连接测试脚本 | ⭐⭐⭐ 先运行 |
| `download_cgmlst_panels.py` | 基础 Python 下载脚本 | ⭐⭐⭐ |
| `download_cgmlst_enhanced.py` | 增强版下载脚本（并发、断点续传） | ⭐⭐⭐⭐⭐ |
| `download_cgmlst.sh` | Shell 版本（使用 curl） | ⭐⭐⭐⭐ |

## 快速开始

### 第一步：测试 API 连接

```bash
python test_cgmlst_api.py
```

这个脚本会测试所有可能的 API 端点，并显示：
- 哪些端点可用
- 返回的数据结构
- 下载链接是否有效

### 第二步：选择下载方式

#### 方式 A：增强版 Python 脚本（推荐）

```bash
# 基础使用
python download_cgmlst_enhanced.py

# 指定参数
python download_cgmlst_enhanced.py \
    --output /path/to/save \
    --threads 5 \
    --delay 0.5

# 强制重新下载（不跳过已存在的）
python download_cgmlst_enhanced.py --no-skip
```

**特点：**
- 自动并发下载
- 断点续传支持
- 智能重试机制
- 详细的下载报告

#### 方式 B：基础 Python 脚本

```bash
python download_cgmlst_panels.py
```

**特点：**
- 简单直接
- 单线程顺序下载
- 适合测试或小规模下载

#### 方式 C：Shell 脚本（使用 curl）

```bash
# 添加执行权限
chmod +x download_cgmlst.sh

# 测试连接
./download_cgmlst.sh --test

# 开始下载
./download_cgmlst.sh
```

**特点：**
- 不依赖 Python
- 使用系统 curl
- 自动检测 jq 工具

## 下载的数据结构

下载完成后，每个 scheme 会有一个独立目录：

```
cgmlst_panels/
├── 0001_Escherichia_coli_ecoli_cgMLST/
│   ├── scheme_metadata.json    # 基础元数据
│   ├── scheme_details.json     # 详细描述
│   ├── alleles.fasta          # 等位基因序列 (FASTA)
│   ├── profiles.tsv           # ST 分型表 (TSV)
│   └── loci_list.txt          # 基因位点列表
├── 0002_Staphylococcus_aureus_saureus_cgMLST/
│   └── ...
├── download_report.json       # 下载统计报告
└── ...
```

## 常见问题

### Q: 脚本无法连接网站？

A: 可能的解决方案：
1. 检查网络连接
2. 尝试使用代理：
   ```bash
   export HTTP_PROXY=http://proxy.example.com:8080
   export HTTPS_PROXY=http://proxy.example.com:8080
   ```
3. 使用浏览器开发者工具查看实际 API 端点

### Q: 返回 429 Too Many Requests？

A: 增加请求延迟：
```bash
python download_cgmlst_enhanced.py --delay 3 --threads 1
```

### Q: 下载速度太慢？

A: 可以调整并发数：
```bash
python download_cgmlst_enhanced.py --threads 10 --delay 0.5
```

**注意：** 并发过高可能导致被暂时封禁 IP，建议从低并发开始测试。

### Q: 如何只下载特定物种？

A: 下载完成后使用 shell 命令过滤：
```bash
# 只保留包含 Salmonella 的目录
find cgmlst_panels -maxdepth 1 -type d ! -name '*Salmonella*' ! -name 'cgmlst_panels' -exec rm -rf {} +
```

### Q: 下载中断如何恢复？

A: 增强版脚本支持自动断点续传：
```bash
# 重新运行即可，会自动跳过已完成的
python download_cgmlst_enhanced.py
```

## cgMLST 数据使用说明

下载的 panel 可用于：

1. **细菌分型**
   ```bash
   # 使用 chewBBACA
   chewBBACA.py AlleleCall -i input/ -g alleles.fasta -o output/
   ```

2. **构建自定义数据库**
   ```bash
   # 合并所有 alleles
   cat */alleles.fasta > all_cgmlst_alleles.fasta
   makeblastdb -in all_cgmlst_alleles.fasta -dbtype nucl
   ```

3. **MLST 分型**
   ```bash
   # 使用 mlst 工具
   mlst --scheme your_scheme genome.fasta
   ```

## 相关工具推荐

| 工具 | 用途 | 链接 |
|------|------|------|
| chewBBACA | cgMLST 等位基因调用 | https://github.com/B-UMMI/chewBBACA |
| MLST | 传统 7 基因 MLST | https://github.com/tseemann/mlst |
| Blast | 序列比对 | NCBI BLAST+ |
| seqkit | 序列处理 | https://bioinf.shenwei.me/seqkit/ |

## 技术支持

cgMLST.org 官方资源：
- 网站：https://www.cgmlst.org
- NCS：https://www.cgmlst.org/ncs

如果遇到 API 变更，请检查：
1. 官方网站的 API 文档
2. 浏览器开发者工具中的网络请求
3. 更新脚本中的 `BASE_URL` 配置

## 许可说明

下载的数据遵循 cgMLST.org 的使用条款。请确保：
- 遵守数据共享协议
- 在发表论文时引用原始数据来源
- 不用于商业用途（如有限制）
