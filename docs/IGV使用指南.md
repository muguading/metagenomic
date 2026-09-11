# IGV 使用与复用指南（本软件）

> 目标：说明本软件如何把 IGV.js 嵌入病毒分析报告，并给出可直接迁移到其他项目的资产约定、后端接口和前端集成步骤。
>
> 证据标记：**[KNOWN]** 表示已由本仓库当前源码直接核对；**[COMMON]** 表示通用的生物信息学/前端实践；**[INFERRED]** 表示基于源码行为作出的解释。未标记内容不构成事实陈述。

## 1. 本软件中的 IGV 是什么

- **[KNOWN]** 本软件在报告页嵌入的是 **IGV.js**，本地脚本路径为 `/public/igv/igv.min.js`；前端调用 `igv.createBrowser()` 创建浏览器实例，而不是要求用户另行安装桌面版 IGV。
- **[KNOWN]** IGV 用于展示样本 reads 对当前参考序列的 BAM 比对，并可选加载 GFF/GFF3 基因注释轨道。
- **[KNOWN]** 报告只有在 `section.igv.status === "ready"` 时才显示“IGV 比对结果”卡片；资产不齐时，后端返回 `{"status": "empty"}`，该卡片不显示。
- **[INFERRED]** 因此，IGV 是对分型、覆盖度和突变表的可视化复核层，而不是独立的分型流程。**置信度：HIGH**。

## 2. 报告使用者：如何查看和核查位点

1. **[KNOWN]** 在已完成的病毒分析任务中打开结果报告，向下找到“IGV 比对结果”卡片。
2. **[KNOWN]** 点击“加载 IGV”；报告初次打开时不会立即创建 IGV iframe，以减少初始页面负载。
3. **[KNOWN]** 在突变表中点击某个位点所在的单元格或行；前端按照“染色体”和“位置”两列组成 `contig:position`，并让 IGV 跳转到该坐标。
4. **[KNOWN]** IGV 卡片包含导航栏和标尺，能够移动、缩放和输入坐标进行查看。
5. **[COMMON]** 对关键变异可结合深度轨道、reads 支持、正反向 reads、一致性以及邻近低复杂度区域进行人工复核；IGV 图像应与变异表和流程质控结果一并解读。
6. **[COMMON]** 记录复核结果时，应保存样本名、参考序列/contig、坐标、参考/替代碱基、查看倍率和结论，以便他人重现查看位置。

### 当前界面的联动规则

- **[KNOWN]** 首次点选突变表位点会触发 IGV 的按需加载，并将该位点作为初始坐标。
- **[KNOWN]** 已加载时，前端优先调用 iframe 内部的 `jump(locus)`，其后调用 `browser.search(locus)`；备用路径使用 `postMessage` 和搜索框事件。
- **[KNOWN]** 若位点 contig 不存在于当前参考序列，iframe 会拒绝跳转并显示“当前 IGV 参考序列不包含该 chromosome”的调试状态。
- **[KNOWN]** 前端对未就绪的 iframe 会延迟重试跳转，最大重试次数为 6 次。

## 3. 生成 IGV 所需的结果资产

### 最小资产约定

| 用途 | 推荐相对路径 | 是否必需 | 说明 |
|---|---|---:|---|
| 参考序列 | `genomes/ref.fa` | 是 | **[KNOWN]** 当前大多数发现函数使用该路径；其他项目可改名，但要同步传入配置。 |
| 参考索引 | `genomes/ref.fa.fai` | 是 | **[KNOWN]** 作为 `reference.indexURL`。 |
| 比对文件 | `ref.mapping.bam` | 是 | **[KNOWN]** 作为 alignment track 的 `url`。 |
| BAM 索引 | `ref.mapping.bam.bai` | 是 | **[KNOWN]** 作为 alignment track 的 `indexURL`。 |
| 基因注释 | `ref/genes.gff` 或 `vadr/*.vadr.gff3` | 视模式而定 | **[KNOWN]** 有值时以 annotation track 加载；部分模式允许缺失。 |

- **[KNOWN]** 流感、HIV、肠道病毒、肝炎病毒、Bandavirus、Orthohantavirus、Astrovirus 和 HAdV 的实现允许在基础四个资产存在时省略 GFF；其余模式是否需要 GFF 由各自的 `_discover_*_igv_assets()` 函数决定。
- **[KNOWN]** 新冠使用数据库中的 `ref.fna`、`ref.fna.fai` 和 `genomic.gff`，通过 `__ncov_*` 虚拟资产名暴露；猴痘优先使用本次 VADR GFF3，缺失时尝试数据库注释。
- **[COMMON]** BAM、BAI、FASTA 和 FAI 必须来自同一参考序列和同一坐标命名体系；否则坐标跳转或 reads 展示可能失败。

### 后端返回的数据结构

**[KNOWN]** 当前报告数据在 `sections.<模块>.igv` 中使用如下契约；其他项目只要提供同名字段即可复用同一类前端逻辑。

```json
{
  "status": "ready",
  "reference_asset": "genomes/ref.fa",
  "reference_index_asset": "genomes/ref.fa.fai",
  "bam_asset": "ref.mapping.bam",
  "bam_index_asset": "ref.mapping.bam.bai",
  "gff_asset": "ref/genes.gff",
  "viewer_label": "参考比对视图",
  "note": "点击上方变异位点后，IGV 会自动跳转到对应位置。",
  "alignment_visibility_window": 5000,
  "alignment_sampling_window_size": 50,
  "alignment_sampling_depth": 50,
  "alignment_downsample_reads": true,
  "alignment_height": 220
}
```

- **[KNOWN]** `gff_asset` 为可选字段；其他字段中的采样和高度也是可选字段，前端对缺失值使用上例中的默认值。
- **[KNOWN]** ZIKV、CHIKV 和 Ebola 的发现函数覆盖默认值为：`visibility_window=3000`、`sampling_window_size=25`、`sampling_depth=25`、`height=180`。
- **[COMMON]** 大深度数据优先设置采样窗口、采样深度和较小轨道高度，再根据实际复核需要扩大显示；这能降低浏览器渲染负担。

## 4. 本软件的实现结构

```text
分析流程产物
  ├── genomes/ref.fa + genomes/ref.fa.fai
  ├── ref.mapping.bam + ref.mapping.bam.bai
  └── ref/genes.gff 或 vadr/*.gff3
          ↓
igv_assets.py: _discover_*_igv_assets()
          ↓
serotype_reports.py: 在报告 section 中填入 igv 配置
          ↓
report-data API: 返回 section.igv
          ↓
report_runtime.js: 按需构建 iframe srcdoc
          ↓
igv.createBrowser(): 加载 reference、BAM 和可选 GFF/GFF3
```

- **[KNOWN]** 资产发现逻辑位于 `bac_analysis_portal/igv_assets.py`，报告数据装配位于 `bac_analysis_portal/serotype_reports.py`。
- **[KNOWN]** iframe 内容、IGV options 和表格联动逻辑位于 `bac_analysis_portal/static/report_runtime.js` 的 `buildInfluenzaIgvSrcdoc()`、`initializeDeferredIgvEmbed()` 和 `initializeInfluenzaIgvLink()`。
- **[KNOWN]** 文件由 `GET /api/tasks/<task_id>/report-asset/<path:asset_name>` 提供，当前接口要求登录、选择可见任务、解析所选样本对应报告目录，并限制普通资产必须位于报告目录内。

## 5. 迁移到其他项目：最小可用实现

### 步骤 A：在分析流程产出索引文件

**[COMMON]** 在报告生成前检查下列文件均存在且非空：`ref.fa`、`ref.fa.fai`、`ref.mapping.bam`、`ref.mapping.bam.bai`；注释轨道启用时也检查 GFF/GFF3。当前工程的发现函数正是以这些文件检查决定 `ready` 或 `empty`。

```python
from pathlib import Path


def discover_igv_assets(report_dir: Path) -> dict:
    required = {
        "reference_asset": report_dir / "genomes" / "ref.fa",
        "reference_index_asset": report_dir / "genomes" / "ref.fa.fai",
        "bam_asset": report_dir / "ref.mapping.bam",
        "bam_index_asset": report_dir / "ref.mapping.bam.bai",
    }
    if not report_dir.is_dir() or not all(path.is_file() for path in required.values()):
        return {"status": "empty"}

    result = {
        "status": "ready",
        **{key: str(path.relative_to(report_dir)) for key, path in required.items()},
        "viewer_label": "参考比对视图",
        "note": "点击突变位点后，IGV 会跳转到对应位置。",
    }
    gff = report_dir / "ref" / "genes.gff"
    if gff.is_file() and gff.stat().st_size > 0:
        result["gff_asset"] = str(gff.relative_to(report_dir))
    return result
```

- **[KNOWN]** 以上模式与本项目 `_discover_*_igv_assets()` 的通用行为一致；它是可复制的最小示例，特定项目可替换文件名与注释来源。

### 步骤 B：提供受控的报告资产 URL

**[COMMON]** 将每个相对资产名转换成 URL，并在服务端验证路径仍位于该任务的报告目录；不要直接接受任意本地文件路径。

```javascript
function buildReportAssetUrl(taskId, assetName) {
  return `/api/tasks/${encodeURIComponent(taskId)}/report-asset/${
    String(assetName).split("/").filter(Boolean).map(encodeURIComponent).join("/")
  }`;
}
```

- **[KNOWN]** 本软件的同名函数使用该编码思路，并保留当前报告所选样本参数。
- **[KNOWN]** 当前资产接口为 BAM/BAI、FASTA/FAI、GFF/GFF3 设置了相应响应类型，并关闭缓存。

### 步骤 C：创建 IGV.js 浏览器

```javascript
const options = {
  showNavigation: true,
  showRuler: true,
  loadDefaultGenomes: false,
  reference: {
    fastaURL: buildReportAssetUrl(task.id, igvView.reference_asset),
    indexURL: buildReportAssetUrl(task.id, igvView.reference_index_asset),
  },
  tracks: [{
    name: igvView.viewer_label || "参考比对视图",
    type: "alignment",
    format: "bam",
    url: buildReportAssetUrl(task.id, igvView.bam_asset),
    indexURL: buildReportAssetUrl(task.id, igvView.bam_index_asset),
    displayMode: "SQUISHED",
    height: igvView.alignment_height || 220,
    visibilityWindow: igvView.alignment_visibility_window || 5000,
    samplingWindowSize: igvView.alignment_sampling_window_size || 50,
    samplingDepth: igvView.alignment_sampling_depth || 50,
    downsampleReads: igvView.alignment_downsample_reads !== false,
  }],
};

if (igvView.gff_asset) {
  options.tracks.push({
    name: "注释",
    type: "annotation",
    format: "gff3",
    url: buildReportAssetUrl(task.id, igvView.gff_asset),
    displayMode: "EXPANDED",
    visibilityWindow: 300000000,
  });
}

const browser = await igv.createBrowser(document.getElementById("igv-div"), options);
```

- **[KNOWN]** 该示例保留了本软件的核心 reference、alignment、annotation、采样和显示配置；当前工程还将 `showSoftClips`、`showInsertionText` 和 `showMismatches` 均设为 `false`。
- **[COMMON]** 将 IGV.js 放在项目自己的静态资源目录中，可以避免运行时依赖外部 CDN；部署时仍应验证该静态文件能被浏览器成功请求。

### 步骤 D：让突变表跳转至 IGV

```javascript
function locusFromVariant(row) {
  return `${String(row.contig).trim()}:${String(row.position).trim()}`;
}

function jumpToVariant(browser, locus) {
  if (!locus || !browser || typeof browser.search !== "function") return;
  browser.search(locus);
}
```

- **[KNOWN]** 本软件根据中文列名“染色体”和“位置”构造同样的 `contig:position` 坐标，并优先通过 `browser.search()` 完成跳转。
- **[COMMON]** 变异表的 contig 值必须与 FASTA 头部中的参考序列名一致；这是联动成功的首要数据一致性检查。

## 6. 运行与故障排查

| 现象 | 核查顺序 | 处理方式 |
|---|---|---|
| 报告没有 IGV 卡片 | **[KNOWN]** 检查 `section.igv.status`；检查发现函数要求的文件。 | **[KNOWN]** 让后端返回完整的 `ready` 配置；缺少基础资产时卡片不会渲染。 |
| 卡片显示但空白 | **[COMMON]** 在浏览器网络面板检查 `igv.min.js`、FASTA、FAI、BAM、BAI 与 GFF 请求。 | **[COMMON]** 修正静态资源路径、任务资产 URL 或文件可读性后重新加载。 |
| 点击位点不跳转 | **[KNOWN]** 查看前端 IGV 调试状态，确认已构造 `contig:position`。 | **[COMMON]** 核对突变表 contig、FASTA 标题和 BAM 参考字典是否完全一致。 |
| IGV 很慢或卡顿 | **[KNOWN]** 当前工程可通过 `alignment_*` 字段配置下采样和轨道高度。 | **[COMMON]** 先降低 sampling depth/window，再用窄区间复核目标位点。 |
| 没有注释轨道 | **[KNOWN]** 当前实现仅在 `gff_asset` 有值时添加 annotation track。 | **[COMMON]** 生成有效的 GFF/GFF3，或在报告配置中省略该字段以只显示比对轨道。 |

## 7. 交付前验证清单

- **[COMMON]** 用 `samtools quickcheck` 或等效工具验证 BAM/BAM 索引可读；工具和版本应与项目环境一致。
- **[COMMON]** 用 FASTA 索引工具或流程质控确认 `.fai` 与 FASTA 对应。
- **[KNOWN]** 请求报告数据后确认 `section.igv.status` 为 `ready`，并确认四个基础资产名非空。
- **[KNOWN]** 在浏览器中加载一次 IGV，点击至少一个突变表位点，并确认浏览器导航坐标变化。
- **[COMMON]** 在高深度样本和含多个 contig 的样本上各做一次人工验收，以覆盖性能与坐标命名两类常见问题。

## 8. 本工程的源码索引

- **[KNOWN]** 资产发现：`/Users/wuhhh/Desktop/徐老师/代码/metagenomic/bac_analysis_portal/igv_assets.py`
- **[KNOWN]** 报告 section 装配：`/Users/wuhhh/Desktop/徐老师/代码/metagenomic/bac_analysis_portal/serotype_reports.py`
- **[KNOWN]** IGV iframe、延迟加载、位点联动：`/Users/wuhhh/Desktop/徐老师/代码/metagenomic/bac_analysis_portal/static/report_runtime.js`
- **[KNOWN]** 报告资产接口：`/Users/wuhhh/Desktop/徐老师/代码/metagenomic/bac_analysis_portal/task_report_routes.py`

**[INFERRED] 结论：** 其他项目只要稳定地产出“参考 FASTA + FAI、BAM + BAI、可选 GFF/GFF3”，并把这些相对路径填入 `section.igv` 契约，就能复用本软件的 IGV 报告模式。**置信度：HIGH**。
