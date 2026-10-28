// Descriptive guidance only. Actual analysis steps remain controlled by the task form.
const MODULE_GUIDES = {
  bacteria: {
    label: "细菌分析",
    description: "适用于纯培养或目标单菌样本。",
    steps: [
      ["质控与预处理", "检查测序质量，按需过滤和去宿主。"],
      ["基因组组装", "根据输入类型选择组装或参考比对。"],
      ["物种鉴定", "结合参考数据库确认物种。"],
      ["分型与功能分析", "按运行方案执行分型、耐药及毒力分析。"],
    ],
    outputs: ["组装序列与质量统计", "物种、MLST 与血清型", "耐药及毒力分析结果"],
  },
  virus: {
    label: "病毒分析",
    description: "适用于目标病毒样本的参考驱动分析。",
    steps: [
      ["输入与参考选择", "指定病毒类型、测序数据及参考序列。"],
      ["质控与参考比对", "检查测序数据并比对目标参考。"],
      ["一致性序列生成", "生成目标病毒的一致性序列。"],
      ["病毒分型", "按病毒类型执行对应分型流程。"],
    ],
    outputs: ["一致性序列与覆盖统计", "病毒鉴定与分型结果", "样本分析报告"],
  },
  metagenome: {
    label: "宏基因组分析",
    description: "适用于混合样本及复杂微生物背景。",
    steps: [
      ["质控与背景处理", "按需过滤测序数据并处理宿主背景。"],
      ["物种识别", "识别混合样本中的微生物组成。"],
      ["宏基因组组装", "根据数据类型组织混合样本组装。"],
      ["结果汇总", "汇总物种及组装分析结果。"],
    ],
    outputs: ["物种列表与丰度统计", "宏基因组组装结果", "样本质量与分析报告"],
  },
  community: {
    label: "群落分析",
    description: "适用于多样本群落项目及分组比较。",
    steps: [
      ["数据与分组", "提供丰度数据、分类信息及样本元数据。"],
      ["多样性分析", "按方案计算 alpha / beta 多样性。"],
      ["差异与标志物", "比较分组差异并筛选候选特征。"],
      ["共现网络", "按需开展共现网络分析。"],
    ],
    outputs: ["多样性统计与图表", "差异特征与标志物结果", "共现网络分析结果"],
  },
  pathosource: {
    label: "溯源分析",
    description: "适用于菌株关联分析及分子溯源。",
    steps: [
      ["菌株与参考选择", "组织待比较菌株及参考序列。"],
      ["比对与变异处理", "按方案执行序列比对及变异处理。"],
      ["系统发育与距离", "构建系统发育树并计算菌株距离。"],
      ["关联分析", "结合分型与聚类结果开展关联判读。"],
    ],
    outputs: ["系统发育树与距离矩阵", "MLST / cgMLST 分型结果", "聚类与菌株关联结果"],
  },
};

export function getSubmissionModuleGuide(key) {
  return Object.hasOwn(MODULE_GUIDES, key) ? MODULE_GUIDES[key] : null;
}

function escapeGuideText(value) {
  return String(value).replace(/[&<>"']/g, (character) => ({
    "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;",
  })[character]);
}

export function renderSubmissionModuleGuide(root, key) {
  const guide = getSubmissionModuleGuide(key);
  if (!root) return guide;
  root.hidden = !guide;
  if (!guide) return null;
  if (root.dataset?.moduleGuide === key) return guide;
  if (root.dataset) root.dataset.moduleGuide = key;
  root.querySelector("[data-guide-title]").textContent = `${guide.label} · 流程说明`;
  root.querySelector("[data-guide-description]").textContent = guide.description;
  root.querySelector("[data-guide-steps]").innerHTML = guide.steps
    .map(([title, description]) => `<li><strong>${escapeGuideText(title)}</strong><p>${escapeGuideText(description)}</p></li>`)
    .join("");
  root.querySelector("[data-guide-outputs]").innerHTML = guide.outputs
    .map((item) => `<li>${escapeGuideText(item)}</li>`).join("");
  return guide;
}
