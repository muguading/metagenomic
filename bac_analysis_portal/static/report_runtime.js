function escapeHtml(value) {
  return String(value ?? "")
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;")
    .replaceAll("'", "&#39;");
}

function formatDateTime(value) {
  if (!value) return "-";
  const date = new Date(value);
  return Number.isNaN(date.getTime()) ? String(value) : date.toLocaleString();
}

function formatBases(value) {
  const numeric = Number(value);
  if (!Number.isFinite(numeric)) return "-";
  if (numeric >= 1e9) return `${(numeric / 1e9).toFixed(2)} Gb`;
  if (numeric >= 1e6) return `${(numeric / 1e6).toFixed(2)} Mb`;
  if (numeric >= 1e3) return `${(numeric / 1e3).toFixed(2)} Kb`;
  return `${numeric} bp`;
}

function formatRate(value) {
  const numeric = Number(value);
  if (!Number.isFinite(numeric)) return "-";
  if (numeric <= 1) return `${(numeric * 100).toFixed(2)}%`;
  return `${numeric.toFixed(2)}%`;
}

function formatChartValue(value, formatter) {
  const numeric = Number(value);
  if (!Number.isFinite(numeric)) return "-";
  if (formatter === "percent") return formatRate(numeric);
  if (formatter === "bases") return formatBases(numeric);
  return numeric.toFixed(numeric >= 100 ? 0 : 2);
}

function parsePercentDisplay(value) {
  const numeric = Number(String(value ?? "").replace("%", "").trim());
  return Number.isFinite(numeric) ? numeric : null;
}

function getQualityCardState(q20, q30) {
  if (q20 !== null && q30 !== null && q20 < 90 && q30 < 80) return "danger";
  return "success";
}

function getQMetricState(label, value) {
  if (value === null) return "neutral";
  if (label === "Q20") return value >= 90 ? "success" : "danger";
  if (label === "Q30") return value >= 80 ? "success" : "danger";
  return "neutral";
}

function getCompletenessState(value) {
  if (value === null) return "neutral";
  if (value > 90) return "success";
  if (value >= 50) return "warning";
  return "danger";
}

function getContaminationState(value) {
  if (value === null) return "neutral";
  if (value < 5) return "success";
  if (value <= 20) return "warning";
  return "danger";
}

function clampPercent(value) {
  const numeric = Number(value);
  if (!Number.isFinite(numeric)) return 0;
  return Math.max(0, Math.min(100, numeric));
}

function buildChartInsight(text) {
  if (!text) return "";
  return `
    <div class="chart-insight" role="note" aria-label="图表判读">
      <span class="chart-insight-label">图表判读</span>
      <p>${escapeHtml(text)}</p>
    </div>
  `;
}

function summarizeRawQc(section, readLabel) {
  const q20 = Number(section?.before_summary?.q20_rate);
  const q30 = Number(section?.before_summary?.q30_rate);
  const gc = Number(section?.before_summary?.gc_content);
  const meanLength = Number(section?.before_summary?.mean_length);
  const qualityPieces = [];
  if (Number.isFinite(q20) && Number.isFinite(q30)) {
    if (q20 >= 0.9 && q30 >= 0.8) {
      qualityPieces.push(`${readLabel} 原始测序质量整体达标`);
    } else if (q20 < 0.85 || q30 < 0.75) {
      qualityPieces.push(`${readLabel} 原始测序质量偏低，建议重点关注过滤后结果`);
    } else {
      qualityPieces.push(`${readLabel} 原始测序质量接近阈值，建议结合后续质控结果综合判读`);
    }
  }
  if (Number.isFinite(gc)) {
    if (gc >= 0.35 && gc <= 0.65) {
      qualityPieces.push(`GC 比例处于常见范围`);
    } else {
      qualityPieces.push(`GC 比例偏离常见范围`);
    }
  }
  if (Number.isFinite(meanLength)) {
    qualityPieces.push(`平均读长约为 ${meanLength.toFixed(0)} bp`);
  }
  return qualityPieces.join("，") || "";
}

function summarizeInsertSize(fastp) {
  const peak = Number(fastp?.insert_size?.peak);
  const unknown = Number(fastp?.insert_size?.unknown);
  if (!Number.isFinite(peak) && !Number.isFinite(unknown)) return "";
  const summary = [];
  if (Number.isFinite(peak)) {
    summary.push(`插入片段主峰位于 ${peak.toFixed(0)} bp`);
  }
  if (Number.isFinite(unknown)) {
    if (unknown > 500000) {
      summary.push(`未识别配对 reads 较多`);
    } else {
      summary.push(`未识别配对 reads 处于可接受范围`);
    }
  }
  return summary.join("，");
}

function summarizeBaseDistribution(baseCurves, readLabel) {
  const gcValues = Array.isArray(baseCurves?.GC) ? baseCurves.GC : [];
  if (!gcValues.length) return "";
  const gcPercent = gcValues.map((value) => Number(value) * 100).filter(Number.isFinite);
  if (!gcPercent.length) return "";
  const min = Math.min(...gcPercent);
  const max = Math.max(...gcPercent);
  const mean = gcPercent.reduce((sum, value) => sum + value, 0) / gcPercent.length;
  const spread = max - min;
  if (spread <= 8) {
    return `${readLabel} 碱基组成整体平稳，GC 比例均值约为 ${mean.toFixed(1)}%，沿读长波动较小。`;
  }
  return `${readLabel} 碱基组成沿读长存在一定波动，GC 比例均值约为 ${mean.toFixed(1)}%，建议结合接头与过滤结果综合判断。`;
}

function summarizeCoverage(section) {
  const mean = Number(section?.mean_depth);
  const max = Number(section?.max_depth);
  const contigCount = Number(section?.contig_count);
  const notes = [];
  if (Number.isFinite(mean)) {
    if (mean >= 100) {
      notes.push(`整体覆盖深度充足`);
    } else if (mean >= 30) {
      notes.push(`整体覆盖深度基本可用`);
    } else {
      notes.push(`整体覆盖深度偏低`);
    }
  }
  if (Number.isFinite(mean) && Number.isFinite(max) && mean > 0) {
    const ratio = max / mean;
    if (ratio >= 8) {
      notes.push(`局部存在明显高深度峰值`);
    } else {
      notes.push(`覆盖深度分布相对平稳`);
    }
  }
  if (Number.isFinite(contigCount)) {
    notes.push(`共覆盖 ${contigCount.toFixed(0)} 条组装片段`);
  }
  return notes.join("，");
}

function summarizeGeneLengthDistribution(section) {
  const labels = Array.isArray(section?.x_values) ? section.x_values : [];
  const points = Array.isArray(section?.points) ? section.points.map((value) => Number(value) || 0) : [];
  if (!labels.length || !points.length) return "";
  const maxValue = Math.max(...points);
  const maxIndex = points.indexOf(maxValue);
  const total = points.reduce((sum, value) => sum + value, 0);
  const topThree = [...points].sort((a, b) => b - a).slice(0, 3).reduce((sum, value) => sum + value, 0);
  const concentration = total > 0 ? (topThree / total) * 100 : 0;
  const dominantLabel = labels[maxIndex] || "主要长度区间";
  if (concentration >= 55) {
    return `基因长度主要集中于 ${dominantLabel} 等少数区间，长度分布呈明显集中趋势。`;
  }
  return `基因长度以 ${dominantLabel} 区间最为常见，但整体仍呈较分散分布。`;
}

function summarizeContigDepth(section) {
  const points = Array.isArray(section?.points) ? section.points : [];
  if (!points.length) return "";
  const grouped = { "基因组": [], "质粒": [] };
  points.forEach((point) => {
    const type = point?.type === "质粒" ? "质粒" : "基因组";
    const depth = Number(point?.depth);
    if (Number.isFinite(depth)) grouped[type].push(depth);
  });
  const avgGenome = grouped["基因组"].length ? grouped["基因组"].reduce((a, b) => a + b, 0) / grouped["基因组"].length : null;
  const avgPlasmid = grouped["质粒"].length ? grouped["质粒"].reduce((a, b) => a + b, 0) / grouped["质粒"].length : null;
  if (avgGenome === null && avgPlasmid === null) return "";
  if (avgGenome !== null && avgPlasmid !== null) {
    if (avgPlasmid > avgGenome * 1.3) {
      return `质粒序列平均深度高于基因组背景，提示部分质粒可能存在较高拷贝特征。`;
    }
    if (avgGenome > avgPlasmid * 1.3) {
      return `基因组序列整体深度高于质粒，提示质粒丰度相对有限。`;
    }
    return `基因组与质粒序列平均深度接近，整体丰度差异不明显。`;
  }
  return avgPlasmid !== null ? `当前结果以质粒序列深度信息为主。` : `当前结果以基因组序列深度信息为主。`;
}

function summarizeRelationship(section) {
  const leftNodes = Array.isArray(section?.nodes_left) ? section.nodes_left : [];
  const rightNodes = Array.isArray(section?.nodes_right) ? section.nodes_right : [];
  const links = Array.isArray(section?.links) ? section.links : [];
  if (!leftNodes.length || !rightNodes.length || !links.length) return "";
  const dominantLeft = leftNodes[0];
  const dominantRight = rightNodes[0];
  const dominantLink = [...links].sort((a, b) => (Number(b.value) || 0) - (Number(a.value) || 0))[0];
  const leftName = dominantLeft?.name || section?.left_label || "主要分类";
  const rightName = dominantRight?.name || section?.right_label || "主要基因";
  if (dominantLink?.source && dominantLink?.target) {
    return `${leftName} 为当前最主要的分类来源，${dominantLink.source} 与 ${dominantLink.target} 的关联最为集中。`;
  }
  return `${leftName} 与 ${rightName} 为当前结果中的主要关系主体。`;
}

function getMetricByKey(metrics, key) {
  return Array.isArray(metrics) ? metrics.find((metric) => metric?.key === key) : null;
}

function buildSummaryCard({ icon, label, title, body, meta, state = "neutral", target = "" }) {
  return `
    <article class="summary-brief-card${target ? " is-jumpable" : ""}" data-state="${escapeHtml(state)}"${target ? ` data-summary-target="${escapeHtml(target)}" tabindex="0" role="button"` : ""}>
      <div class="summary-brief-head">
        <span class="summary-brief-icon" aria-hidden="true">${escapeHtml(icon || "•")}</span>
        <span class="summary-brief-label">${escapeHtml(label)}</span>
      </div>
      <h3>${escapeHtml(title || "-")}</h3>
      <p>${escapeHtml(body || "-")}</p>
      ${meta ? `<span class="summary-brief-meta">${escapeHtml(meta)}</span>` : ""}
    </article>
  `;
}

function countPrioritySerotypeHits(section) {
  const columns = Array.isArray(section?.columns) ? section.columns : [];
  const rows = Array.isArray(section?.rows) ? section.rows : [];
  if (!columns.length || !rows.length) return 0;
  const serotypeIndex = columns.indexOf("血清型");
  const virulenceIndex = columns.indexOf("毒力基因");
  if (serotypeIndex < 0 || virulenceIndex < 0) return 0;
  return rows.filter((row) => {
    const serotype = String(row?.[serotypeIndex] ?? "").trim();
    const virulence = String(row?.[virulenceIndex] ?? "").trim();
    const validSerotype = serotype && serotype !== "-";
    const validVirulence = virulence && virulence !== "-";
    return validSerotype && validVirulence;
  }).length;
}

function buildExecutiveSummary(data) {
  const metrics = Array.isArray(data?.overview_metrics) ? data.overview_metrics : [];
  const sections = data?.sections || {};

  const qMetric = getMetricByKey(metrics, "q_metrics");
  const q20 = parsePercentDisplay(qMetric?.items?.[0]?.display);
  const q30 = parsePercentDisplay(qMetric?.items?.[1]?.display);
  const qualityState = getQualityCardState(q20, q30);
  const qualityTitle = qualityState === "danger" ? "测序质量需重点关注" : "测序质量总体合格";
  const qualityBody = qualityState === "danger"
    ? `Q20 ${qMetric?.items?.[0]?.display || "--"}、Q30 ${qMetric?.items?.[1]?.display || "--"}，建议优先复核原始质控与 fastp 结果。`
    : `Q20 ${qMetric?.items?.[0]?.display || "--"}、Q30 ${qMetric?.items?.[1]?.display || "--"}，原始测序质量达到当前报告判读阈值。`;

  const assemblyMetric = getMetricByKey(metrics, "assembly_profile");
  const contigCount = Number(assemblyMetric?.contig_count);
  const plasmidCount = Number(assemblyMetric?.plasmid_count);
  const checkmMetric = getMetricByKey(metrics, "checkm_metrics");
  const completeness = parsePercentDisplay(checkmMetric?.items?.[0]?.display);
  const contamination = parsePercentDisplay(checkmMetric?.items?.[1]?.display);
  let assemblyState = "neutral";
  if (Number.isFinite(contigCount) && contigCount > 200) {
    assemblyState = "danger";
  } else if ((Number.isFinite(completeness) && completeness < 90) || (Number.isFinite(contamination) && contamination >= 5)) {
    assemblyState = "warning";
  } else if (Number.isFinite(contigCount)) {
    assemblyState = "success";
  }
  const assemblyTitle = assemblyState === "danger"
    ? "组装碎片化程度偏高"
    : assemblyState === "warning"
      ? "组装结果基本可用，建议结合质量指标判读"
      : "组装结果整体稳健";
  const assemblyBody = `当前共检出 ${Number.isFinite(contigCount) ? contigCount : "--"} 条 Contig、${Number.isFinite(plasmidCount) ? plasmidCount : "--"} 条质粒；完整性 ${checkmMetric?.items?.[0]?.display || "--"}，污染率 ${checkmMetric?.items?.[1]?.display || "--"}。`;

  const speciesMetric = getMetricByKey(metrics, "species_estimation");
  const speciesName = speciesMetric?.items?.[0]?.display || "--";
  const mlstSpecies = speciesMetric?.items?.[1]?.display || "--";
  const speciesTitle = speciesName && speciesName !== "--" ? speciesName : "未形成稳定物种预估";
  const speciesBody = mlstSpecies && mlstSpecies !== "--"
    ? `综合 CheckM 与 MLST 结果，当前样本更接近 ${mlstSpecies} 对应的分型背景。`
    : "当前报告未提供稳定的 MLST 物种辅助信息，建议以分类与分型结果综合判断。";

  const resistanceRows = sections?.resistance_virulence?.resistance_elements?.rows?.length || 0;
  const virulenceRows = sections?.resistance_virulence?.virulence_elements?.rows?.length || 0;
  const priorityRows = countPrioritySerotypeHits(sections?.priority_serotype || {});
  let focusState = "neutral";
  let focusTitle = "建议按模块顺序判读";
  let focusBody = "当前未检出需优先跳读的高风险模块，可依次查看组装、分型与耐药毒力结果。";
  if (priorityRows > 0) {
    focusState = "danger";
    focusTitle = "建议优先查看关注毒力血清型";
    focusBody = `当前已检出 ${priorityRows} 条重点关注结果，建议先结合血清型与毒力信息进行判读。`;
  } else if (resistanceRows > 0 || virulenceRows > 0) {
    focusState = resistanceRows > 0 && virulenceRows > 0 ? "warning" : "neutral";
    focusTitle = "建议重点查看耐药与毒力分析";
    focusBody = `当前耐药元件 ${resistanceRows} 条、毒力元件 ${virulenceRows} 条，建议结合关系图与汇总表共同判读。`;
  }

  return [
    { icon: "Q", label: "Sequencing", title: qualityTitle, body: qualityBody, meta: "原始数据质控 + fastp", state: qualityState, target: "#section-raw-qc" },
    { icon: "A", label: "Assembly", title: assemblyTitle, body: assemblyBody, meta: "组装情况 + CheckM", state: assemblyState, target: "#section-assembly" },
    { icon: "T", label: "Taxonomy", title: speciesTitle, body: speciesBody, meta: "物种预估", state: "neutral", target: "#section-species" },
    { icon: "!", label: "Priority", title: focusTitle, body: focusBody, meta: "优先阅读建议", state: focusState, target: priorityRows > 0 ? "#section-priority-serotype" : "#section-rv" },
  ];
}

function renderExecutiveSummary(data) {
  const container = document.getElementById("report-summary-grid");
  if (!container) return;
  container.innerHTML = buildExecutiveSummary(data).map((card) => buildSummaryCard(card)).join("");
  container.querySelectorAll("[data-summary-target]").forEach((card) => {
    const jump = () => {
      const selector = card.dataset.summaryTarget || "";
      const target = selector ? document.querySelector(selector) : null;
      if (!target) return;
      target.scrollIntoView({ behavior: "smooth", block: "start" });
      if (selector.startsWith("#")) {
        history.replaceState(null, "", selector);
      }
    };
    card.addEventListener("click", jump);
    card.addEventListener("keydown", (event) => {
      if (event.key === "Enter" || event.key === " ") {
        event.preventDefault();
        jump();
      }
    });
  });
}

function buildMetricCards(metrics) {
  const totalMetric = metrics.find((metric) => metric.key === "total_bases");
  const qMetric = metrics.find((metric) => metric.key === "q_metrics");
  const mergedMetrics = metrics.filter((metric) => metric.key !== "q_metrics").map((metric) => {
    if (metric.key === "total_bases" && totalMetric && qMetric) {
      return { ...metric, type: "quality_overview", items: qMetric.items || [] };
    }
    return metric;
  });

  return mergedMetrics.map((metric) => {
    if (metric.type === "quality_overview") {
      const q20 = parsePercentDisplay(metric.items?.[0]?.display);
      const q30 = parsePercentDisplay(metric.items?.[1]?.display);
      const cardState = getQualityCardState(q20, q30);
      return `
        <article class="metric-card quality-overview-card metric-state-${cardState}">
          <span class="metric-label">${escapeHtml(metric.label)}</span>
          <strong>${escapeHtml(metric.display ?? "--")}</strong>
          <div class="paired-metric-grid compact-paired-grid">
            ${(metric.items || []).map((item) => {
              const value = parsePercentDisplay(item.display);
              const itemState = getQMetricState(item.label, value);
              return `
              <div class="paired-metric-item inline-metric-item metric-state-${itemState}">
                <span>${escapeHtml(item.label)}</span>
                <strong>${escapeHtml(item.display ?? "--")}</strong>
                <div class="metric-progress-track" aria-hidden="true">
                  <div class="metric-progress-fill metric-progress-${itemState}" style="width:${clampPercent(value)}%"></div>
                </div>
              </div>
            `;}).join("")}
          </div>
          <div class="metric-status-note">${cardState === "danger" ? "测序质量差" : "测序质量合格"}</div>
        </article>
      `;
    }
    if (metric.type === "paired") {
      const items = (metric.items || []).map((item) => {
        let state = "neutral";
        const value = parsePercentDisplay(item.display);
        if (metric.key === "checkm_metrics") {
          state = item.label.includes("完整") ? getCompletenessState(value) : getContaminationState(value);
        }
        return { ...item, state, numericValue: value };
      });
      return `
        <article class="metric-card paired-metric-card">
          <span class="metric-label">${escapeHtml(metric.label)}</span>
          <div class="paired-metric-grid">
            ${items.map((item) => `
              <div class="paired-metric-item metric-state-${item.state}">
                <span>${escapeHtml(item.label)}</span>
                <strong>${escapeHtml(item.display ?? "--")}</strong>
                ${metric.key === "checkm_metrics" ? `
                  <div class="metric-progress-track" aria-hidden="true">
                    <div class="metric-progress-fill metric-progress-${item.state}" style="width:${clampPercent(item.numericValue)}%"></div>
                  </div>
                ` : ""}
              </div>
            `).join("")}
          </div>
        </article>
      `;
    }
    if (metric.type === "assembly_profile") {
      const contig = Number(metric.contig_count ?? 0);
      const plasmid = Number(metric.plasmid_count ?? 0);
      const total = Math.max(Number(metric.total_count ?? 0), contig + plasmid, 1);
      const contigState = contig > 200 ? "danger" : "success";
      const contigNote = contig > 200 ? "组装片段过多" : "组装质量良好";
      const contigRatio = Math.max(0, Math.min(1, contig / total));
      const plasmidRatio = Math.max(0, Math.min(1, plasmid / total));
      const contigDash = `${(contigRatio * 251.2).toFixed(2)} 251.2`;
      const plasmidDash = `${(plasmidRatio * 251.2).toFixed(2)} 251.2`;
      return `
        <article class="metric-card assembly-metric-card metric-state-${contigState}">
          <span class="metric-label">${escapeHtml(metric.label)}</span>
          <div class="assembly-metric-layout">
            <div class="donut-chart" aria-label="组装情况环形图">
              <svg viewBox="0 0 120 120" role="img">
                <circle class="donut-track" cx="60" cy="60" r="40"></circle>
                <circle class="donut-segment donut-segment-contig" cx="60" cy="60" r="40" stroke-dasharray="${contigDash}"></circle>
                <circle class="donut-segment donut-segment-plasmid" cx="60" cy="60" r="40" stroke-dasharray="${plasmidDash}" stroke-dashoffset="-${(contigRatio * 251.2).toFixed(2)}"></circle>
              </svg>
              <div class="donut-total">${escapeHtml(String(metric.total_count ?? "--"))}</div>
            </div>
            <div class="assembly-metric-values">
              <div class="paired-metric-item">
                <span>总长度</span>
                <strong>${escapeHtml(formatBases(metric.total_length))}</strong>
              </div>
              <div class="paired-metric-item">
                <span>Contig</span>
                <strong>${escapeHtml(String(metric.contig_count ?? "--"))}</strong>
              </div>
              <div class="paired-metric-item">
                <span>质粒</span>
                <strong>${escapeHtml(String(metric.plasmid_count ?? "--"))}</strong>
              </div>
            </div>
          </div>
          <div class="metric-status-note">${contigNote}</div>
        </article>
      `;
    }
    return `
      <article class="metric-card">
        <span class="metric-label">${escapeHtml(metric.label)}</span>
        <strong>${escapeHtml(metric.display ?? "--")}</strong>
      </article>
    `;
  }).join("");
}


function compareTableValues(left, right, direction) {
  const leftNumber = Number(left);
  const rightNumber = Number(right);
  const bothNumeric = Number.isFinite(leftNumber) && Number.isFinite(rightNumber);
  if (bothNumeric) {
    return direction === "asc" ? leftNumber - rightNumber : rightNumber - leftNumber;
  }
  const leftText = String(left ?? "");
  const rightText = String(right ?? "");
  return direction === "asc"
    ? leftText.localeCompare(rightText, "zh-CN", { numeric: true })
    : rightText.localeCompare(leftText, "zh-CN", { numeric: true });
}

function extractNumericValue(value) {
  if (typeof value === "number" && Number.isFinite(value)) return value;
  const text = String(value ?? "").trim();
  if (!text) return null;
  const normalized = text.replaceAll(",", "");
  const match = normalized.match(/-?\d+(?:\.\d+)?/);
  if (!match) return null;
  const numeric = Number(match[0]);
  return Number.isFinite(numeric) ? numeric : null;
}

function inferColumnKinds(columns, rows) {
  return columns.map((column, index) => {
    let numericCount = 0;
    let valueCount = 0;
    rows.forEach((row) => {
      const value = Array.isArray(row) ? row[index] : row[column];
      const text = String(value ?? "").trim();
      if (!text || text === "-") return;
      valueCount += 1;
      if (extractNumericValue(text) !== null) numericCount += 1;
    });
    return valueCount > 0 && numericCount / valueCount >= 0.8 ? "number" : "text";
  });
}

function matchNumericFilter(value, rawFilter) {
  const numeric = extractNumericValue(value);
  if (numeric === null) return false;
  const filter = String(rawFilter ?? "").trim().replace(/\s+/g, "");
  if (!filter) return true;

  const rangeMatch = filter.match(/^(-?\d+(?:\.\d+)?)\s*(?:-|~)\s*(-?\d+(?:\.\d+)?)$/);
  if (rangeMatch) {
    const min = Number(rangeMatch[1]);
    const max = Number(rangeMatch[2]);
    if (!Number.isFinite(min) || !Number.isFinite(max)) return false;
    return numeric >= Math.min(min, max) && numeric <= Math.max(min, max);
  }

  const compareMatch = filter.match(/^(<=|>=|=|<|>)(-?\d+(?:\.\d+)?)$/);
  if (compareMatch) {
    const operator = compareMatch[1];
    const target = Number(compareMatch[2]);
    if (!Number.isFinite(target)) return false;
    if (operator === "<") return numeric < target;
    if (operator === "<=") return numeric <= target;
    if (operator === ">") return numeric > target;
    if (operator === ">=") return numeric >= target;
    return numeric === target;
  }

  const exact = Number(filter);
  if (Number.isFinite(exact)) {
    return numeric === exact;
  }
  return String(value ?? "").toLowerCase().includes(String(rawFilter ?? "").trim().toLowerCase());
}

function renderTableCellContent(value) {
  const text = String(value ?? "-");
  if (text.length <= 30) {
    return `<span>${escapeHtml(text)}</span>`;
  }
  return `<span title="${escapeHtml(text)}">${escapeHtml(`${text.slice(0, 30)}...`)}</span>`;
}

const tableExportStore = new WeakMap();
let currentReportData = null;

function slugifyFilename(value) {
  const normalized = String(value ?? "").trim().replace(/[^\w\u4e00-\u9fff-]+/g, "_");
  return normalized.replace(/^_+|_+$/g, "") || "export";
}

function downloadBlob(blob, filename) {
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = filename;
  document.body.appendChild(link);
  link.click();
  link.remove();
  window.setTimeout(() => URL.revokeObjectURL(url), 1000);
}

function renderTableExportToolbar() {
  return `
    <div class="table-export-toolbar">
      <div class="table-export-menu">
        <button class="table-export-button table-export-trigger" type="button">导出表格</button>
        <div class="table-export-dropdown">
          <button class="table-export-option" type="button" data-table-export-format="xlsx">XLSX</button>
          <button class="table-export-option" type="button" data-table-export-format="tsv">TSV</button>
          <button class="table-export-option" type="button" data-table-export-format="csv">CSV</button>
        </div>
      </div>
    </div>
  `;
}

function renderChartExportToolbar() {
  return `
    <div class="chart-export-toolbar">
      <div class="chart-export-menu">
        <button class="chart-export-button chart-export-trigger" type="button">下载图表</button>
        <div class="chart-export-dropdown">
          <button class="chart-export-option" type="button" data-chart-export-format="svg">SVG</button>
          <button class="chart-export-option" type="button" data-chart-export-format="png">PNG</button>
          <button class="chart-export-option" type="button" data-chart-export-format="pdf">PDF</button>
        </div>
      </div>
    </div>
  `;
}

function getChartExportFilename(container) {
  const title = container.querySelector('.mini-chart-title')?.textContent?.trim()
    || container.closest('.result-card')?.querySelector('.card-title-stack h3, .card-head h3')?.textContent?.trim()
    || 'chart';
  return slugifyFilename(title);
}

function serializeChartSvg(svg) {
  const clone = svg.cloneNode(true);
  clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
  clone.setAttribute('xmlns:xlink', 'http://www.w3.org/1999/xlink');
  if (!clone.getAttribute('viewBox')) {
    const width = svg.viewBox?.baseVal?.width || svg.clientWidth || 1200;
    const height = svg.viewBox?.baseVal?.height || svg.clientHeight || 480;
    clone.setAttribute('viewBox', `0 0 ${width} ${height}`);
  }
  return new XMLSerializer().serializeToString(clone);
}

async function requestChartExport(container, format) {
  const svg = container.querySelector('svg');
  if (!svg) return;
  const filename = getChartExportFilename(container);
  const svgText = serializeChartSvg(svg);
  const svgBlob = new Blob([svgText], { type: 'image/svg+xml;charset=utf-8' });
  if (format === 'svg') {
    downloadBlob(svgBlob, `${filename}.svg`);
    return;
  }
  if (format === 'png') {
    const url = URL.createObjectURL(svgBlob);
    try {
      const img = new Image();
      img.decoding = 'async';
      await new Promise((resolve, reject) => {
        img.onload = resolve;
        img.onerror = reject;
        img.src = url;
      });
      const viewBox = svg.viewBox?.baseVal;
      const width = Math.max(1, Math.round(viewBox?.width || svg.clientWidth || img.width || 1200));
      const height = Math.max(1, Math.round(viewBox?.height || svg.clientHeight || img.height || 480));
      const canvas = document.createElement('canvas');
      canvas.width = width * 2;
      canvas.height = height * 2;
      const ctx = canvas.getContext('2d');
      if (!ctx) throw new Error('图表导出失败');
      ctx.scale(2, 2);
      ctx.fillStyle = '#fbf8f1';
      ctx.fillRect(0, 0, width, height);
      ctx.drawImage(img, 0, 0, width, height);
      const pngBlob = await new Promise((resolve) => canvas.toBlob(resolve, 'image/png'));
      if (!pngBlob) throw new Error('PNG 导出失败');
      downloadBlob(pngBlob, `${filename}.png`);
    } finally {
      URL.revokeObjectURL(url);
    }
    return;
  }
  const title = container.querySelector('.mini-chart-title')?.textContent?.trim() || '图表';
  const printWindow = window.open('', '_blank', 'noopener,noreferrer');
  if (!printWindow) return;
  printWindow.document.write(`<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><title>${escapeHtml(title)}</title><style>
    @page { size: A4 landscape; margin: 12mm; }
    body { margin: 0; font-family: "Source Sans 3", "PingFang SC", sans-serif; background: #fffaf2; color: #243447; }
    .sheet { padding: 16px 18px; }
    h1 { margin: 0 0 10px; font-size: 20px; }
    p { margin: 0 0 16px; color: #526170; }
    .frame { border: 1px solid rgba(36,52,71,0.12); border-radius: 16px; padding: 12px; background: #fffdf8; }
    svg { width: 100%; height: auto; display: block; }
  </style></head><body><div class="sheet"><h1>${escapeHtml(title)}</h1><p>图表导出副本</p><div class="frame">${svgText}</div></div><script>window.addEventListener('load',function(){setTimeout(function(){window.print();},300);});</script></body></html>`);
  printWindow.document.close();
}

function bindChartExportButtons(root = document) {
  root.querySelectorAll('.mini-chart-card').forEach((card) => {
    if (card.dataset.chartEnhanced === 'true') return;
    if (!card.querySelector('svg')) return;
    const title = card.querySelector('.mini-chart-title');
    if (!title) return;
    const bar = document.createElement('div');
    bar.className = 'chart-titlebar';
    title.replaceWith(bar);
    bar.appendChild(title);
    bar.insertAdjacentHTML('beforeend', renderChartExportToolbar());
    card.dataset.chartEnhanced = 'true';
  });
  root.querySelectorAll('[data-chart-export-format]').forEach((button) => {
    if (button.dataset.bound === 'true') return;
    button.dataset.bound = 'true';
    button.addEventListener('click', async (event) => {
      event.stopPropagation();
      try {
        await requestChartExport(button.closest('.mini-chart-card'), button.dataset.chartExportFormat || 'svg');
      } catch (error) {
        console.error(error);
      }
    });
  });
}

function renderInteractiveTableSummary(columns, columnKinds, state, totalRows, filteredCount) {
  const activeFilters = columns
    .map((column, index) => ({ column, value: String(state.filters[index] || "").trim(), kind: columnKinds[index], index }))
    .filter((item) => item.value);
  const sortLabel = state.sortIndex >= 0
    ? `${columns[state.sortIndex]} · ${state.sortDirection === "asc" ? "升序" : "降序"}`
    : "未排序";
  return `
    <div class="table-statebar">
      <div class="table-state-toprow">
        <div class="table-state-meta">
          <span class="table-state-pill">显示 ${filteredCount} / ${totalRows} 行</span>
          <span class="table-state-pill">${escapeHtml(sortLabel)}</span>
          ${activeFilters.length ? `<span class="table-state-pill">已筛选 ${activeFilters.length} 列</span>` : ""}
        </div>
        <div class="table-view-switch" role="tablist" aria-label="表格视图切换">
          <button class="table-view-button ${state.viewMode === "key" ? "active" : ""}" type="button" data-table-view="key">关键列</button>
          <button class="table-view-button ${state.viewMode === "all" ? "active" : ""}" type="button" data-table-view="all">全部列</button>
        </div>
      </div>
      <div class="table-filter-chip-row">
        ${activeFilters.length ? activeFilters.map((item) => `
          <button class="table-filter-chip" type="button" data-clear-filter-index="${item.index}">
            <span>${escapeHtml(item.column)}</span>
            <strong>${escapeHtml(item.value)}</strong>
          </button>
        `).join("") : `<span class="table-filter-placeholder">当前未设置筛选条件</span>`}
      </div>
    </div>
  `;
}

function getInteractiveTableSpec(tableId, columns) {
  const aliasSets = {
    "assembly-species-table": ["序列名称", "基因组/质粒", "平均深度", "物种名称", "界", "门", "属", "种"],
    "contig-annotation-table": ["seq_name", "length", "cov.", "circ.", "repeat", "mult.", "alt_group"],
    "rv-summary-table": ["序列名称", "物种名称", "平均深度", "基因组/质粒", "毒力基因", "耐药基因", "血清型"],
    "virulence-table": ["Contig名称", "物种名称", "基因名称", "VF分类", "覆盖度%", "一致性%", "平均深度"],
    "resistance-table": ["Contig名称", "物种名称", "基因名称", "耐药药物", "覆盖度%", "一致性%", "平均深度"],
  };
  const preferred = aliasSets[tableId] || [];
  const normalized = columns.map((column) => String(column || "").toLowerCase());
  const keyIndexes = preferred
    .map((alias) => normalized.findIndex((column) => column.includes(String(alias).toLowerCase())))
    .filter((index, position, list) => index >= 0 && list.indexOf(index) === position);
  return { keyIndexes: keyIndexes.length ? keyIndexes : columns.slice(0, Math.min(columns.length, 6)).map((_, index) => index) };
}

function getTableCellTone(column, value) {
  const text = String(value ?? "-").trim();
  const numeric = extractNumericValue(text);
  if (!text || text === "-") return { tone: "", render: `<span>${escapeHtml(text || "-")}</span>` };
  if (["基因组/质粒", "是否成环", "VF分类", "耐药药物", "物种名称", "属", "种"].includes(column)) {
    return { tone: "tag", render: `<span class="table-inline-tag">${escapeHtml(text)}</span>` };
  }
  if (column.includes("完整性")) {
    const tone = getCompletenessState(numeric);
    return { tone, render: `<span class="table-inline-chip table-inline-chip-${tone}">${escapeHtml(text)}</span>` };
  }
  if (column.includes("污染率")) {
    const tone = getContaminationState(numeric);
    return { tone, render: `<span class="table-inline-chip table-inline-chip-${tone}">${escapeHtml(text)}</span>` };
  }
  if (column.includes("覆盖度") && column.includes("%")) {
    const tone = numeric >= 95 ? "success" : numeric >= 80 ? "warning" : "danger";
    return { tone, render: `<span class="table-inline-chip table-inline-chip-${tone}">${escapeHtml(text)}</span>` };
  }
  if (column.includes("一致性")) {
    const tone = numeric >= 95 ? "success" : numeric >= 85 ? "warning" : "danger";
    return { tone, render: `<span class="table-inline-chip table-inline-chip-${tone}">${escapeHtml(text)}</span>` };
  }
  if (column.includes("平均深度")) {
    const tone = numeric >= 30 ? "success" : numeric >= 10 ? "warning" : "danger";
    return { tone, render: `<span class="table-inline-chip table-inline-chip-${tone}">${escapeHtml(text)}</span>` };
  }
  return { tone: "", render: renderTableCellContent(value) };
}

async function requestTableExport(container, format) {
  const shell = document.querySelector(".report-shell");
  const payload = tableExportStore.get(container);
  if (!shell || !payload) return;
  const response = await fetch(shell.dataset.exportEndpoint, {
    method: "POST",
    credentials: "same-origin",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({
      title: payload.title,
      filename: slugifyFilename(payload.title),
      columns: payload.columns,
      rows: payload.rows,
      format,
    }),
  });
  if (!response.ok) {
    const data = await response.json().catch(() => ({}));
    throw new Error(data.error || "表格导出失败");
  }
  const blob = await response.blob();
  downloadBlob(blob, `${slugifyFilename(payload.title)}.${format}`);
}

function bindTableExportButtons(container, title, columns, rows) {
  tableExportStore.set(container, { title, columns, rows });
  container.querySelectorAll("[data-table-export-format]").forEach((button) => {
    button.addEventListener("click", async () => {
      try {
        await requestTableExport(container, button.dataset.tableExportFormat || "csv");
      } catch (error) {
        console.error(error);
      }
    });
  });
}

async function fetchReportCssText() {
  const stylesheet = document.querySelector('link[href*="report.css"]');
  if (!stylesheet) return "";
  const response = await fetch(stylesheet.href, { credentials: "same-origin" });
  return response.ok ? response.text() : "";
}

async function fetchReportJsText() {
  const runtimeScript = document.querySelector('script[src*="report_runtime.js"]');
  if (!runtimeScript?.src) return "";
  const response = await fetch(runtimeScript.src, { credentials: "same-origin" });
  return response.ok ? response.text() : "";
}

async function buildReportExportHtml(printMode = false, interactiveMode = false) {
  const clone = document.documentElement.cloneNode(true);
  const shell = clone.querySelector(".report-shell");
  const taskName = document.querySelector(".report-shell")?.dataset.taskName || document.getElementById("report-task-name")?.textContent || "analysis_report";
  const taskStatus = document.querySelector(".report-status")?.textContent?.trim() || "未知";
  const exportedAt = new Intl.DateTimeFormat("zh-CN", {
    year: "numeric",
    month: "2-digit",
    day: "2-digit",
    hour: "2-digit",
    minute: "2-digit",
    second: "2-digit",
    hour12: false,
  }).format(new Date());
  clone.body?.classList.add("report-export-document");
  if (printMode) clone.body?.classList.add("report-print-export");
  if (interactiveMode) clone.body?.classList.add("report-interactive-export");
  clone.querySelectorAll("script").forEach((node) => node.remove());
  clone.querySelectorAll(".report-export-toolbar, .table-export-toolbar, .report-back").forEach((node) => node.remove());
  if (!interactiveMode) {
    clone.querySelectorAll(".table-filter-input, .table-statebar, .table-view-switch").forEach((node) => node.remove());
    clone.querySelectorAll(".table-sort-button").forEach((button) => {
      const replacement = clone.ownerDocument.createElement("span");
      replacement.textContent = button.querySelector("span")?.textContent || button.textContent || "";
      button.replaceWith(replacement);
    });
    clone.querySelectorAll(".mlst-hostgene-button").forEach((button) => {
      const replacement = clone.ownerDocument.createElement("span");
      replacement.textContent = button.textContent || "";
      button.replaceWith(replacement);
    });
    clone.querySelectorAll("td span[title]").forEach((span) => {
      span.textContent = span.getAttribute("title") || span.textContent || "";
      span.removeAttribute("title");
    });
  }
  clone.querySelectorAll('link[href*="report.css"]').forEach((node) => node.remove());
  const cssText = await fetchReportCssText();
  if (cssText) {
    const style = clone.ownerDocument.createElement("style");
    style.textContent = cssText;
    clone.querySelector("head")?.appendChild(style);
  }
  if (shell) {
    const masthead = clone.ownerDocument.createElement("section");
    masthead.className = "report-document-masthead";
    masthead.innerHTML = `
      <div class="report-document-headline">
        <span class="report-document-label">${printMode ? "打印版报告" : interactiveMode ? "交互导出版" : "归档导出版"}</span>
        <h2>细菌组装分析结果归档副本</h2>
        <p>本副本用于结果存档、离线查看与纸质输出。表格与图形已整理为正式报告阅读顺序。</p>
      </div>
      <dl class="report-document-meta">
        <div><dt>任务名称</dt><dd>${escapeHtml(taskName)}</dd></div>
        <div><dt>任务状态</dt><dd>${escapeHtml(taskStatus)}</dd></div>
        <div><dt>导出时间</dt><dd>${escapeHtml(exportedAt)}</dd></div>
        <div><dt>导出形式</dt><dd>${printMode ? "PDF打印" : interactiveMode ? "HTML交互文档" : "Word归档文档"}</dd></div>
      </dl>
    `;
    shell.insertBefore(masthead, shell.querySelector(".report-layout"));
  }
  if (interactiveMode && currentReportData) {
    if (shell) {
      shell.dataset.reportEndpoint = "";
      shell.dataset.exportEndpoint = "";
    }
    const dataScript = clone.ownerDocument.createElement("script");
    dataScript.textContent = `window.__EMBEDDED_REPORT_DATA__ = ${JSON.stringify(currentReportData)};`;
    clone.querySelector("body")?.appendChild(dataScript);
    const runtimeText = await fetchReportJsText();
    if (runtimeText) {
      const script = clone.ownerDocument.createElement("script");
      script.textContent = runtimeText;
      clone.querySelector("body")?.appendChild(script);
    }
  }
  if (printMode) {
    const script = clone.ownerDocument.createElement("script");
    script.textContent = "window.addEventListener('load', function(){ setTimeout(function(){ window.print(); }, 300); });";
    clone.querySelector("body")?.appendChild(script);
  }
  return `<!doctype html>
${clone.outerHTML}`;
}

async function exportReportPage(format) {
  const shell = document.querySelector(".report-shell");
  const baseName = slugifyFilename(shell?.dataset.taskName || "analysis_report");
  if (format === "pdf") {
    const html = await buildReportExportHtml(true, false);
    const blob = new Blob([html], { type: "text/html;charset=utf-8" });
    const url = URL.createObjectURL(blob);
    window.open(url, "_blank");
    window.setTimeout(() => URL.revokeObjectURL(url), 60000);
    return;
  }
  const html = await buildReportExportHtml(false, format === "html");
  if (format === "word") {
    downloadBlob(new Blob([html], { type: "application/msword;charset=utf-8" }), `${baseName}.doc`);
    return;
  }
  downloadBlob(new Blob([html], { type: "text/html;charset=utf-8" }), `${baseName}.html`);
}

function renderInteractiveContigTable(container, columns, rows, tableId = "") {
  const columnKinds = inferColumnKinds(columns, rows);
  const spec = getInteractiveTableSpec(tableId, columns);
  const state = {
    sortIndex: -1,
    sortDirection: "asc",
    filters: columns.map(() => ""),
    scrollLeft: 0,
    scrollTop: 0,
    viewMode: spec.keyIndexes?.length ? "key" : "all",
  };

  const applyState = (focusIndex = null, caretPosition = null) => {
    const previousFrame = container.querySelector(".table-frame");
    if (previousFrame) {
      state.scrollLeft = previousFrame.scrollLeft;
      state.scrollTop = previousFrame.scrollTop;
    }

    const filteredRows = rows.filter((row) => columns.every((column, index) => {
      const keyword = state.filters[index]?.trim();
      if (!keyword) return true;
      const value = Array.isArray(row) ? row[index] : row[column];
      if (columnKinds[index] === "number") {
        return matchNumericFilter(value, keyword);
      }
      return String(value ?? "").toLowerCase().includes(keyword.toLowerCase());
    }));

    if (state.sortIndex >= 0) {
      filteredRows.sort((leftRow, rightRow) => {
        const left = Array.isArray(leftRow) ? leftRow[state.sortIndex] : leftRow[columns[state.sortIndex]];
        const right = Array.isArray(rightRow) ? rightRow[state.sortIndex] : rightRow[columns[state.sortIndex]];
        return compareTableValues(left, right, state.sortDirection);
      });
    }

    const visibleIndexes = state.viewMode === "key" && spec.keyIndexes?.length
      ? spec.keyIndexes
      : columns.map((_, index) => index);
    const exportTitle = container.dataset.exportTitle || "结果表";
    container.innerHTML = `
      ${renderTableExportToolbar()}
      ${renderInteractiveTableSummary(columns, columnKinds, state, rows.length, filteredRows.length)}
      <div class="table-frame interactive-table-frame tall-table-frame">
        <table class="report-table report-table-interactive">
          <thead>
            <tr>
              ${visibleIndexes.map((index) => {
                const column = columns[index];
                return `
                <th>
                  <div class="table-head-stack">
                    <span class="table-head-kicker">${index === 0 ? "主索引" : (columnKinds[index] === "number" ? "数值列" : "文本列")}</span>
                    <button class="table-sort-button" type="button" data-sort-index="${index}">
                      <span>${escapeHtml(column)}</span>
                      <span class="table-sort-indicator">${state.sortIndex === index ? (state.sortDirection === "asc" ? "▲" : "▼") : "↕"}</span>
                    </button>
                  </div>
                  <div class="table-filter-wrap">
                    <input class="table-filter-input" type="text" placeholder="${columnKinds[index] === "number" ? "如 >100 或 50-200" : "筛选"}" value="${escapeHtml(state.filters[index] || "")}" data-filter-index="${index}">
                  </div>
                </th>`;
              }).join("")}
            </tr>
          </thead>
          <tbody>
            ${filteredRows.map((row) => `<tr>${visibleIndexes.map((index) => {
              const column = columns[index];
              const value = Array.isArray(row) ? row[index] : row[column];
              const cell = getTableCellTone(column, value);
              return `<td class="${cell.tone ? `table-cell-${cell.tone}` : ""}">${cell.render}</td>`;
            }).join("")}</tr>`).join("")}
          </tbody>
        </table>
      </div>
    `;
    bindTableExportButtons(container, exportTitle, columns, filteredRows.map((row) => (
      Array.isArray(row) ? row : columns.map((column) => row[column] ?? "")
    )));

    container.querySelectorAll('.table-sort-button').forEach((button) => {
      button.addEventListener('click', () => {
        const index = Number(button.dataset.sortIndex);
        if (state.sortIndex === index) {
          state.sortDirection = state.sortDirection === 'asc' ? 'desc' : 'asc';
        } else {
          state.sortIndex = index;
          state.sortDirection = 'asc';
        }
        applyState();
      });
    });

    container.querySelectorAll('.table-filter-input').forEach((input) => {
      input.addEventListener('input', () => {
        const index = Number(input.dataset.filterIndex);
        state.filters[index] = input.value;
        applyState(index, input.selectionStart ?? input.value.length);
      });
    });

    container.querySelectorAll('[data-clear-filter-index]').forEach((button) => {
      button.addEventListener('click', () => {
        const index = Number(button.dataset.clearFilterIndex);
        state.filters[index] = '';
        state.viewMode = state.viewMode || 'key';
        applyState(index, 0);
      });
    });

    container.querySelectorAll('[data-table-view]').forEach((button) => {
      button.addEventListener('click', () => {
        state.viewMode = button.dataset.tableView === 'all' ? 'all' : 'key';
        applyState();
      });
    });

    const nextFrame = container.querySelector(".table-frame");
    if (nextFrame) {
      nextFrame.scrollLeft = state.scrollLeft;
      nextFrame.scrollTop = state.scrollTop;
      nextFrame.addEventListener("scroll", () => {
        state.scrollLeft = nextFrame.scrollLeft;
        state.scrollTop = nextFrame.scrollTop;
      }, { passive: true });
    }

    if (focusIndex !== null) {
      const target = container.querySelector(`.table-filter-input[data-filter-index="${focusIndex}"]`);
      if (target) {
        target.focus();
        const nextPos = Math.min(caretPosition ?? target.value.length, target.value.length);
        target.setSelectionRange(nextPos, nextPos);
      }
    }
  };

  applyState();
}

function applyTableTone(container, containerId) {
  const tones = {
    "assembly-species-table": "taxonomy",
    "contig-annotation-table": "assembly",
    "assembly-summary-table": "assembly",
    "checkm-table": "quality",
    "gene-annotation-summary-table": "gene",
    "rv-summary-table": "risk",
    "virulence-table": "virulence",
    "resistance-table": "resistance",
    "mlst-table": "typing",
    "serotype-table": "typing",
    "priority-serotype-table": "risk",
    "taxonomy-source-table": "taxonomy",
    "taxonomy-summary-table": "taxonomy",
  };
  const tone = tones[containerId] || "neutral";
  [
    "table-tone-neutral",
    "table-tone-taxonomy",
    "table-tone-assembly",
    "table-tone-quality",
    "table-tone-gene",
    "table-tone-risk",
    "table-tone-virulence",
    "table-tone-resistance",
    "table-tone-typing",
  ].forEach((className) => container.classList.remove(className));
  container.classList.add(`table-tone-${tone}`);
}

function buildTableCard(containerId, title, columns, rows) {
  const container = document.getElementById(containerId);
  if (!container) return;
  container.dataset.exportTitle = title;
  applyTableTone(container, containerId);
  if (!Array.isArray(rows) || rows.length === 0) {
    container.innerHTML = `
      <div class="empty-table-state">
        <strong>${escapeHtml(title)} 暂无数据</strong>
        <p class="empty-copy">本节未检出可展示结果。</p>
      </div>
    `;
    return;
  }
  if (['assembly-species-table', 'contig-annotation-table', 'rv-summary-table', 'virulence-table', 'resistance-table'].includes(containerId)) {
    renderInteractiveContigTable(container, columns, rows, containerId);
    return;
  }
  container.innerHTML = `
    ${renderTableExportToolbar()}
    <div class="table-frame">
      <table class="report-table">
        <thead><tr>${columns.map((column) => `<th>${escapeHtml(column)}</th>`).join("")}</tr></thead>
        <tbody>
          ${rows.map((row) => `<tr>${columns.map((column, index) => {
            const value = Array.isArray(row) ? row[index] : row[column];
            return `<td>${renderTableCellContent(value)}</td>`;
          }).join("")}</tr>`).join("")}
        </tbody>
      </table>
    </div>
  `;
  bindTableExportButtons(container, title, columns, rows.map((row) => (
    Array.isArray(row) ? row : columns.map((column) => row[column] ?? "")
  )));
}

function fillTaskMeta(task) {
  document.getElementById("report-task-name").textContent = task.name || task.id || "-";
  document.getElementById("report-task-meta").textContent = `任务编号：${task.id || "-"}`;
  document.getElementById("report-sample-title").textContent = task.sample_name || task.name || task.id || "分析结果";
  document.getElementById("report-sample-copy").textContent = `创建时间：${formatDateTime(task.created_at)}；开始时间：${formatDateTime(task.started_at)}；结束时间：${formatDateTime(task.finished_at)}。`;
  document.getElementById("meta-owner").textContent = task.owner || "-";
  document.getElementById("meta-group").textContent = task.group || "-";
  document.getElementById("meta-asm-type").textContent = task.asm_type || "-";
  document.getElementById("meta-method").textContent = task.method || "-";
  document.getElementById("meta-input").textContent = task.input_path || "-";
  document.getElementById("meta-output").textContent = task.output_dir || "-";
}

function renderChartLegend(seriesList) {
  if (!Array.isArray(seriesList) || seriesList.length === 0) return "";
  return `
    <div class="chart-legend" aria-hidden="true">
      ${seriesList.map((series) => `
        <span class="chart-legend-item">
          <i class="chart-legend-swatch" style="--legend-color:${escapeHtml(series.color)}"></i>
          <span>${escapeHtml(series.label)}</span>
        </span>
      `).join("")}
    </div>
  `;
}

function createSeriesSvg(seriesList, options = {}) {
  const width = options.width || 1100;
  const height = options.height || 320;
  const padX = options.padX ?? 72;
  const padTop = options.padTop ?? 20;
  const padBottom = options.padBottom ?? 58;
  const innerWidth = width - padX * 2;
  const innerHeight = height - padTop - padBottom;
  const maxLength = Math.max(...seriesList.map((series) => (series.values || []).length), 0);
  const allValues = seriesList.flatMap((series) => (series.values || []).map((value) => Number(value) || 0));
  const computedMax = Math.max(...allValues, 1);
  const minValue = options.minValue ?? 0;
  const roundedMax = options.maxValue ?? (computedMax <= 1 ? 1 : Math.ceil(computedMax * 1.1));
  const yTicks = 4;
  const grid = Array.from({ length: yTicks + 1 }, (_, index) => {
    const ratio = index / yTicks;
    const y = padTop + innerHeight - innerHeight * ratio;
    const value = minValue + (roundedMax - minValue) * ratio;
    return `
      <line class="chart-grid-line" x1="${padX}" y1="${y}" x2="${width - padX}" y2="${y}"></line>
      <text class="chart-axis-label y-axis" x="${padX - 10}" y="${y + 4}">${escapeHtml(formatChartValue(value, options.yFormatter))}</text>
    `;
  }).join("");
  const xValues = Array.isArray(options.xValues) ? options.xValues : [];
  const xTickValues = options.xTicks || [1, Math.max(1, Math.round((maxLength + 1) / 2)), Math.max(1, maxLength)];
  const xTicks = [0, 0.5, 1].map((ratio, index) => {
    const x = padX + innerWidth * ratio;
    return `
      <line class="chart-axis-tick" x1="${x}" y1="${height - padBottom}" x2="${x}" y2="${height - padBottom + 6}"></line>
      <text class="chart-axis-label x-axis" x="${x}" y="${height - 26}">${escapeHtml(String(xTickValues[index] ?? ""))}</text>
    `;
  }).join("");
  const paths = seriesList.map((series) => {
    const values = series.values || [];
    const coords = values.map((value, pointIndex) => {
      const x = padX + (maxLength <= 1 ? 0 : (pointIndex / (maxLength - 1)) * innerWidth);
      const rawNumber = Number(value) || 0;
      const scaledValue = options.yFormatter === "percent" ? rawNumber * 100 : rawNumber;
      const normalized = (scaledValue - minValue) / Math.max(roundedMax - minValue, 1);
      const y = padTop + innerHeight - normalized * innerHeight;
      return { x, y };
    });
    const points = coords.map(({ x, y }) => `${x},${y}`).join(" ");
    const areaPoints = [
      `${padX},${height - padBottom}`,
      ...coords.map(({ x, y }) => `${x},${y}`),
      `${padX + innerWidth},${height - padBottom}`,
    ].join(" ");
    return `
      <polygon class="chart-area-fill" fill="${series.color}" fill-opacity="${seriesList.length > 2 ? "0.04" : "0.08"}" points="${areaPoints}"></polygon>
      <polyline class="chart-series-line" fill="none" stroke="${series.color}" stroke-width="${seriesList.length > 2 ? "2.1" : "2.4"}" points="${points}"></polyline>
    `;
  }).join("");
  const focusDots = seriesList.map((series) => (
    `<circle class="chart-focus-dot" fill="${series.color}" r="4" cx="${padX}" cy="${height - padBottom}"></circle>`
  )).join("");
  return `
    <div class="chart-shell">
      ${renderChartLegend(seriesList)}
      <div
        class="interactive-chart"
        data-chart-padx="${padX}"
        data-chart-padtop="${padTop}"
        data-chart-padbottom="${padBottom}"
        data-chart-innerwidth="${innerWidth}"
        data-chart-innerheight="${innerHeight}"
        data-chart-maxlength="${maxLength}"
        data-chart-minvalue="${minValue}"
        data-chart-maxvalue="${roundedMax}"
        data-chart-xlabel="${escapeHtml(options.xLabel || "位置") }"
        data-chart-ylabel="${escapeHtml(options.yLabel || "数值") }"
        data-chart-yformatter="${escapeHtml(options.yFormatter || "number") }"
        data-chart-xvalues="${escapeHtml(JSON.stringify(xValues))}"
        data-chart-series="${escapeHtml(JSON.stringify(seriesList.map((series) => ({ label: series.label, color: series.color, values: series.values || [] }))))}"
      >
        <div class="chart-canvas" style="--chart-height:${height}px">
          <svg class="sparkline-svg" viewBox="0 0 ${width} ${height}" preserveAspectRatio="xMidYMid meet" role="img" aria-label="${escapeHtml(options.label || "结果曲线")}">
            ${grid}
            <line class="chart-axis-line" x1="${padX}" y1="${padTop}" x2="${padX}" y2="${height - padBottom}"></line>
            <line class="chart-axis-line" x1="${padX}" y1="${height - padBottom}" x2="${width - padX}" y2="${height - padBottom}"></line>
            ${xTicks}
            ${paths}
            <line class="chart-focus-line" x1="${padX}" y1="${padTop}" x2="${padX}" y2="${height - padBottom}"></line>
            ${focusDots}
            <text class="chart-axis-title" x="${width / 2}" y="${height - 6}">${escapeHtml(options.xLabel || "位置")}</text>
            <text class="chart-axis-title chart-axis-title-y" x="24" y="${height / 2}">${escapeHtml(options.yLabel || "数值")}</text>
          </svg>
        </div>
        <div class="chart-tooltip" hidden></div>
      </div>
    </div>
  `;
}

function renderStatsWithCharts(containerId, section, readLabel) {
  const container = document.getElementById(containerId);
  if (!container) return;
  container.classList.remove("empty-box");
  const qualitySvg = createSeriesSvg([
    { label: `${readLabel} mean`, color: "#4e6177", values: section.quality_curves?.mean || [] },
  ], { label: `${readLabel} 数据质量`, width: 820, height: 280, xLabel: "碱基位置", yLabel: "质量分数", maxValue: 45 });
  const cards = [
    { label: "总 reads", value: String(section.before_summary?.total_reads ?? "-") },
    { label: "总数据量", value: formatBases(section.before_summary?.total_bases) },
    { label: "平均长度", value: String(section.before_summary?.mean_length ?? "-") },
    { label: "Q20", value: formatRate(section.before_summary?.q20_rate) },
    { label: "Q30", value: formatRate(section.before_summary?.q30_rate) },
    { label: "GC", value: formatRate(section.before_summary?.gc_content) },
  ];
  container.innerHTML = `
    <div class="mini-stat-grid">
      ${cards.map((item) => `
        <div class="mini-stat-card">
          <span>${escapeHtml(item.label)}</span>
          <strong>${escapeHtml(item.value)}</strong>
        </div>
      `).join("")}
    </div>
    <div class="mini-chart-stack">
      <div class="mini-chart-card">
        <span class="mini-chart-title">数据质量</span>
        ${buildChartInsight(summarizeRawQc(section, readLabel))}
        ${qualitySvg}
      </div>
    </div>
  `;
}

function renderAdapterTab(panel, adapter) {
  panel.innerHTML = `
    <div class="adapter-grid">
      <div class="adapter-card"><span>接头去除 reads</span><strong>${escapeHtml(String(adapter.adapter_trimmed_reads ?? "-"))}</strong></div>
      <div class="adapter-card"><span>接头去除碱基</span><strong>${escapeHtml(formatBases(adapter.adapter_trimmed_bases))}</strong></div>
      <div class="adapter-card"><span>Read1 接头序列</span><code>${escapeHtml(adapter.read1_adapter_sequence || "-")}</code></div>
      <div class="adapter-card"><span>Read2 接头序列</span><code>${escapeHtml(adapter.read2_adapter_sequence || "-")}</code></div>
    </div>
  `;
}

function renderFastpTabs(fastp) {
  const insertPanel = document.getElementById("fastp-tab-insert-size");
  const basePanel = document.getElementById("fastp-tab-base-content");
  const adapterPanel = document.getElementById("fastp-tab-adapter");
  if (insertPanel) {
    insertPanel.innerHTML = `
      <div class="mini-chart-card">
        <span class="mini-chart-title">插入片段长度</span>
        ${buildChartInsight(summarizeInsertSize(fastp))}
        ${createSeriesSvg([{ label: "insert size", color: "#4e6177", values: fastp.insert_size?.histogram || [] }], { label: "插入片段长度", width: 820, height: 240, xLabel: "插入片段长度区间", yLabel: "read 数" })}
        <p class="empty-copy">峰值：${escapeHtml(String(fastp.insert_size?.peak ?? "-"))}；未知配对：${escapeHtml(String(fastp.insert_size?.unknown ?? "-"))}</p>
      </div>
    `;
  }
  if (basePanel) {
    basePanel.innerHTML = `
      <div class="mini-chart-card">
        <div class="subreport-tabs" role="tablist" aria-label="碱基分布切换">
          <button class="subreport-tab-button active" type="button" data-base-tab="read1">R1</button>
          <button class="subreport-tab-button" type="button" data-base-tab="read2">R2</button>
        </div>
        <div class="subreport-tab-panel active" data-base-panel="read1">
          <span class="mini-chart-title">R1 碱基分布</span>
          ${buildChartInsight(summarizeBaseDistribution(fastp.base_distribution?.read1 || {}, "R1"))}
          ${createSeriesSvg([
            { label: "A", color: "#8a6654", values: fastp.base_distribution?.read1?.A || [] },
            { label: "T", color: "#7a7158", values: fastp.base_distribution?.read1?.T || [] },
            { label: "C", color: "#5d7c83", values: fastp.base_distribution?.read1?.C || [] },
            { label: "G", color: "#6d6481", values: fastp.base_distribution?.read1?.G || [] },
            { label: "GC", color: "#3e546f", values: fastp.base_distribution?.read1?.GC || [] },
          ], { label: "R1 碱基分布", width: 820, height: 240, xLabel: "碱基位置", yLabel: "比例", yFormatter: "percent", maxValue: 100 })}
        </div>
        <div class="subreport-tab-panel" data-base-panel="read2">
          <span class="mini-chart-title">R2 碱基分布</span>
          ${buildChartInsight(summarizeBaseDistribution(fastp.base_distribution?.read2 || {}, "R2"))}
          ${createSeriesSvg([
            { label: "A", color: "#8a6654", values: fastp.base_distribution?.read2?.A || [] },
            { label: "T", color: "#7a7158", values: fastp.base_distribution?.read2?.T || [] },
            { label: "C", color: "#5d7c83", values: fastp.base_distribution?.read2?.C || [] },
            { label: "G", color: "#6d6481", values: fastp.base_distribution?.read2?.G || [] },
            { label: "GC", color: "#3e546f", values: fastp.base_distribution?.read2?.GC || [] },
          ], { label: "R2 碱基分布", width: 820, height: 240, xLabel: "碱基位置", yLabel: "比例", yFormatter: "percent", maxValue: 100 })}
        </div>
      </div>
    `;
  }
  if (adapterPanel) {
    renderAdapterTab(adapterPanel, fastp.adapter_cutting || {});
  }
}

function initializeBaseTabs(root = document) {
  root.querySelectorAll(".subreport-tabs").forEach((tabGroup) => {
    if (tabGroup.dataset.initialized === "true") return;
    tabGroup.dataset.initialized = "true";
    const container = tabGroup.closest(".mini-chart-card");
    const buttons = Array.from(tabGroup.querySelectorAll(".subreport-tab-button"));
    const panels = Array.from(container?.querySelectorAll(".subreport-tab-panel") || []);
    buttons.forEach((button) => {
      button.addEventListener("click", () => {
        const target = button.dataset.baseTab;
        buttons.forEach((item) => item.classList.toggle("active", item === button));
        panels.forEach((panel) => {
          panel.classList.toggle("active", panel.dataset.basePanel === target);
        });
        initializeInteractiveCharts(container || document);
      });
    });
  });
}

function initializeInteractiveCharts(root = document) {
  root.querySelectorAll(".interactive-chart").forEach((chart) => {
    if (chart.dataset.initialized === "true") return;
    chart.dataset.initialized = "true";
    const svg = chart.querySelector(".sparkline-svg");
    const tooltip = chart.querySelector(".chart-tooltip");
    const focusLine = chart.querySelector(".chart-focus-line");
    const dots = Array.from(chart.querySelectorAll(".chart-focus-dot"));
    const seriesList = JSON.parse(chart.dataset.chartSeries || "[]");
    const maxLength = Number(chart.dataset.chartMaxlength || 0);
    const padX = Number(chart.dataset.chartPadx || 0);
    const padTop = Number(chart.dataset.chartPadtop || 0);
    const padBottom = Number(chart.dataset.chartPadbottom || 0);
    const innerWidth = Number(chart.dataset.chartInnerwidth || 0);
    const innerHeight = Number(chart.dataset.chartInnerheight || 0);
    const minValue = Number(chart.dataset.chartMinvalue || 0);
    const maxValue = Number(chart.dataset.chartMaxvalue || 1);
    const yFormatter = chart.dataset.chartYformatter || "number";
    const xLabel = chart.dataset.chartXlabel || "位置";
    const xValues = JSON.parse(chart.dataset.chartXvalues || "[]");
    if (!svg || !tooltip || maxLength < 1) return;

    const update = (clientX) => {
      const rect = svg.getBoundingClientRect();
      const relativeX = Math.min(Math.max(clientX - rect.left, 0), rect.width);
      const ratio = rect.width <= 0 ? 0 : relativeX / rect.width;
      const index = Math.min(Math.max(Math.round(ratio * Math.max(maxLength - 1, 0)), 0), Math.max(maxLength - 1, 0));
      const chartX = padX + (maxLength <= 1 ? 0 : (index / (maxLength - 1)) * innerWidth);
      focusLine.setAttribute("x1", chartX);
      focusLine.setAttribute("x2", chartX);
      const rows = [];
      dots.forEach((dot, seriesIndex) => {
        const series = seriesList[seriesIndex] || {};
        const rawValue = Number(series.values?.[index] ?? 0);
        const scaledValue = yFormatter === "percent" ? rawValue * 100 : rawValue;
        const normalized = (scaledValue - minValue) / Math.max(maxValue - minValue, 1);
        const y = padTop + innerHeight - normalized * innerHeight;
        dot.setAttribute("cx", chartX);
        dot.setAttribute("cy", y);
        rows.push(`<span><i style="background:${series.color}"></i>${escapeHtml(series.label)}: ${escapeHtml(formatChartValue(rawValue, yFormatter))}</span>`);
      });
      const xValue = xValues[index] ?? (index + 1);
      tooltip.innerHTML = `<strong>${escapeHtml(`${xLabel}: ${xValue}`)}</strong>${rows.join("")}`;
      tooltip.hidden = false;
    };

    svg.addEventListener("mousemove", (event) => update(event.clientX));
    svg.addEventListener("mouseenter", (event) => update(event.clientX));
    svg.addEventListener("mouseleave", () => {
      tooltip.hidden = true;
    });
  });
}

function renderRiskBars(items, type) {
  return `
    <div class="risk-bar-list risk-bar-list-${type}">
      ${items.map((item) => {
        const width = Math.max(Math.min(Number(item.ratio) || 0, 100), 2);
        return `
          <div class="risk-bar-item">
            <div class="risk-bar-meta">
              <span>${escapeHtml(item.label)}</span>
              <strong>${escapeHtml(String(item.reads))} reads</strong>
            </div>
            <div class="risk-bar-track" aria-hidden="true">
              <div class="risk-bar-fill" style="width:${width}%"></div>
            </div>
            <div class="risk-bar-foot">
              <span>占比 ${escapeHtml(formatRate((Number(item.ratio) || 0) / 100))}</span>
              <span>分类项 ${escapeHtml(String(item.records || 0))}</span>
            </div>
          </div>
        `;
      }).join("")}
    </div>
  `;
}

function renderKingdomCompareBars(items) {
  const totalReads = items.reduce((sum, item) => sum + (Number(item.reads) || 0), 0);
  return `
    <div class="kingdom-compare-list">
      ${items.map((item) => {
        const reads = Number(item.reads) || 0;
        const ratio = totalReads ? (reads / totalReads) * 100 : Number(item.ratio) || 0;
        const width = Math.max(Math.min(ratio, 100), reads > 0 ? 8 : 2);
        return `
          <div class="kingdom-compare-item">
            <div class="kingdom-compare-meta">
              <span class="kingdom-compare-label">${escapeHtml(item.label)}</span>
              <strong>${escapeHtml(String(reads))} reads</strong>
            </div>
            <div class="kingdom-compare-track" aria-hidden="true">
              <div class="kingdom-compare-fill" style="width:${width}%"></div>
            </div>
            <div class="kingdom-compare-foot">
              <span>占比 ${escapeHtml(formatRate(ratio / 100))}</span>
              <span>${escapeHtml(String(item.records || 0))} 项分类</span>
            </div>
          </div>
        `;
      }).join("")}
    </div>
  `;
}

function renderRvCategoryBars(items, type) {
  return `
    <div class="risk-bar-list risk-bar-list-${type}">
      ${items.map((item) => {
        const width = Math.max(Math.min((Number(item.count) || 0) * 10, 100), 4);
        return `
          <div class="risk-bar-item">
            <div class="risk-bar-meta">
              <span>${escapeHtml(item.label)}</span>
              <strong>${escapeHtml(String(item.count))} 项</strong>
            </div>
            <div class="risk-bar-track" aria-hidden="true">
              <div class="risk-bar-fill" style="width:${width}%"></div>
            </div>
          </div>
        `;
      }).join("")}
    </div>
  `;
}

function renderResistanceVirulenceOverview(section) {
  const resistanceContainer = document.getElementById("resistance-summary-panel");
  const virulenceContainer = document.getElementById("virulence-summary-panel");
  if (!resistanceContainer || !virulenceContainer) return;
  const resistance = section?.resistance || {};
  const virulence = section?.virulence || {};

  resistanceContainer.classList.remove("empty-box");
  virulenceContainer.classList.remove("empty-box");

  resistanceContainer.innerHTML = `
    ${buildChartInsight(resistance.note || "当前未检出可汇总的耐药元件结果。")}
    <div class="risk-summary-panel">
      <div class="risk-summary-head">
        <span class="mini-chart-title">耐药药物类别分布</span>
        <span class="risk-summary-count">${escapeHtml(String(resistance.gene_count || 0))} 个基因</span>
      </div>
      <div class="mini-stat-grid">
        <div class="mini-stat-card"><span>元件记录</span><strong>${escapeHtml(String(resistance.hit_count || 0))}</strong></div>
        <div class="mini-stat-card"><span>汇总基因数</span><strong>${escapeHtml(String(resistance.summary_count || 0))}</strong></div>
      </div>
      ${(resistance.top_categories || []).length ? renderRvCategoryBars(resistance.top_categories || [], "hazard") : `<p class="empty-copy">未检出可展示的耐药药物类别分布。</p>`}
    </div>
  `;

  virulenceContainer.innerHTML = `
    ${buildChartInsight(virulence.note || "当前未检出可汇总的毒力元件结果。")}
    <div class="risk-summary-panel">
      <div class="risk-summary-head">
        <span class="mini-chart-title">VF 分类分布</span>
        <span class="risk-summary-count">${escapeHtml(String(virulence.gene_count || 0))} 个基因</span>
      </div>
      <div class="mini-stat-grid">
        <div class="mini-stat-card"><span>元件记录</span><strong>${escapeHtml(String(virulence.hit_count || 0))}</strong></div>
        <div class="mini-stat-card"><span>汇总基因数</span><strong>${escapeHtml(String(virulence.summary_count || 0))}</strong></div>
      </div>
      ${(virulence.top_categories || []).length ? renderRvCategoryBars(virulence.top_categories || [], "pathogenicity") : `<p class="empty-copy">未检出可展示的 VF 分类分布。</p>`}
    </div>
  `;
}

function renderTaxonomyRiskSummary(section) {
  const summaryContainer = document.getElementById("taxonomy-risk-summary");
  const hazardContainer = document.getElementById("taxonomy-hazard-chart");
  if (!summaryContainer || !hazardContainer) return;
  const kingdomSummary = Array.isArray(section?.kingdom_summary) ? section.kingdom_summary : [];
  const pathogenicity = Array.isArray(section?.pathogenicity) ? section.pathogenicity : [];
  const hazard = Array.isArray(section?.hazard) ? section.hazard : [];
  const narrative = section?.narrative || "";

  if (!pathogenicity.length && !hazard.length) {
    summaryContainer.innerHTML = `
      <div class="empty-box coverage-empty">
        <p>当前序列物种鉴定结果未提供可汇总的致病性与危害等级信息。</p>
      </div>
    `;
    hazardContainer.innerHTML = `
      <div class="empty-box coverage-empty">
        <p>未检出可展示的危害程度等级分布。</p>
      </div>
    `;
    return;
  }

  summaryContainer.classList.remove("empty-box");
  hazardContainer.classList.remove("empty-box");

  summaryContainer.innerHTML = `
    <div class="taxonomy-summary-stack">
      ${buildChartInsight(narrative)}
      <div class="taxonomy-summary-feature">
        <div class="risk-summary-head taxonomy-summary-feature-head">
          <div>
            <span class="mini-chart-title">细菌 / 病毒 / 真菌对比</span>
            <p class="taxonomy-summary-caption">基于 3.1 序列物种鉴定结果，对不同大类微生物的序列数量与分类记录进行汇总。</p>
          </div>
          <span class="risk-summary-count">${escapeHtml(String(section?.total_records || 0))} 条分类记录</span>
        </div>
        ${renderKingdomCompareBars(kingdomSummary)}
      </div>
      <div class="taxonomy-summary-subpanel">
        <div class="risk-summary-head">
          <span class="mini-chart-title">致病性分布</span>
          <span class="risk-summary-count">${escapeHtml(String(pathogenicity.length))} 类</span>
        </div>
        ${pathogenicity.length ? renderRiskBars(pathogenicity, "pathogenicity") : `<p class="empty-copy">未检出可展示的致病性结果。</p>`}
      </div>
    </div>
  `;

  hazardContainer.innerHTML = `
    <div class="risk-summary-panel">
      <div class="risk-summary-head">
        <span class="mini-chart-title">危害程度等级</span>
        <span class="risk-summary-count">${escapeHtml(String(hazard.length))} 类</span>
      </div>
      ${hazard.length ? renderRiskBars(hazard, "hazard") : `<p class="empty-copy">未检出可展示的危害程度等级结果。</p>`}
    </div>
  `;
}

function renderSpeciesIdentification(section) {
  const sourceContainer = document.getElementById("taxonomy-source-table");
  const tagsContainer = document.getElementById("taxonomy-rank-tags");
  const tableContainer = document.getElementById("taxonomy-summary-table");
  if (!sourceContainer || !tagsContainer || !tableContainer) return;
  applyTableTone(sourceContainer, "taxonomy-source-table");
  applyTableTone(tableContainer, "taxonomy-summary-table");

  const datasets = {
    species: section?.species || { rows: [], rank_options: [], terminal_column: "种" },
    subspecies: section?.subspecies || { rows: [], rank_options: [], terminal_column: "亚种" },
  };
  const state = {
    tab: datasets.species.rows?.length ? "species" : "subspecies",
    rank: null,
  };

  const aggregateRows = (rows, rank) => {
    const groups = new Map();
    rows.forEach((row) => {
      const key = String(row?.[rank] ?? "").trim() || "未注释";
      const current = groups.get(key) || { ratio: 0, reads: 0 };
      current.ratio += Number(row?.["比例数值"] || 0);
      current.reads += Number(row?.["序列数量数值"] || 0);
      groups.set(key, current);
    });
    return Array.from(groups.entries())
      .map(([name, values]) => [rank, name, `${values.ratio.toFixed(2)}%`, String(values.reads)])
      .sort((left, right) => Number(right[3]) - Number(left[3]));
  };

  const renderSourceTable = (dataset) => {
    const hiddenColumns = new Set(["界", "门", "纲", "目", "科", "比例数值", "序列数量数值"]);
    const rows = dataset.rows || [];
    if (!rows.length) {
      sourceContainer.innerHTML = `
        <div class="empty-table-state">
          <strong>未检出分类明细结果</strong>
          <p class="empty-copy">未检出对应的种或亚种分类结果文件。</p>
        </div>
      `;
      return;
    }
    const sample = rows[0] || {};
    const columns = Object.keys(sample).filter((column) => !hiddenColumns.has(column));
    sourceContainer.dataset.exportTitle = state.tab === "species" ? "序列物种鉴定_种" : "序列物种鉴定_亚种";
    renderInteractiveContigTable(
      sourceContainer,
      columns,
      rows.map((row) => columns.map((column) => row?.[column] ?? "")),
    );
  };

  const render = () => {
    const dataset = datasets[state.tab] || { rows: [], rank_options: [], terminal_column: "" };
    const rankOptions = dataset.rank_options || [];
    if (!state.rank || !rankOptions.includes(state.rank)) {
      state.rank = rankOptions[rankOptions.length - 1] || null;
    }

    document.querySelectorAll('[data-taxonomy-tab]').forEach((button) => {
      button.classList.toggle('active', button.dataset.taxonomyTab === state.tab);
    });

    renderSourceTable(dataset);

    tagsContainer.innerHTML = rankOptions.map((rank) => `
      <button class="taxonomy-rank-tag ${state.rank === rank ? 'active' : ''}" type="button" data-taxonomy-rank="${escapeHtml(rank)}">${escapeHtml(rank)}</button>
    `).join('');

    if (!dataset.rows?.length || !state.rank) {
      tableContainer.innerHTML = `
        <div class="empty-table-state">
          <strong>未检出序列物种鉴定结果</strong>
          <p class="empty-copy">未检出对应的种或亚种分类结果文件。</p>
        </div>
      `;
    } else {
      tableContainer.dataset.exportTitle = `${state.tab === "species" ? "种" : "亚种"}_${state.rank}_分类聚合`;
      renderInteractiveContigTable(
        tableContainer,
        ["分类水平", "分类名称", "比例", "序列数量"],
        aggregateRows(dataset.rows, state.rank),
      );
    }

    tagsContainer.querySelectorAll('[data-taxonomy-rank]').forEach((button) => {
      button.addEventListener('click', () => {
        state.rank = button.dataset.taxonomyRank;
        render();
      });
    });
  };

  document.querySelectorAll('[data-taxonomy-tab]').forEach((button) => {
    button.onclick = () => {
      state.tab = button.dataset.taxonomyTab;
      state.rank = null;
      render();
    };
  });

  render();
}

function renderTaxonomyAbundance(section) {
  const container = document.getElementById("taxonomy-abundance-chart");
  const topnSelect = document.getElementById("taxonomy-abundance-topn");
  const tag = document.getElementById("taxonomy-abundance-tag");
  if (!container) return;
  const ranks = Array.isArray(section?.ranks) ? section.ranks : [];
  if (!ranks.length) {
    container.innerHTML = `
      <div class="empty-box coverage-empty">
        <p>未检出可展示的分类丰度数据。</p>
      </div>
    `;
    return;
  }
  const state = {
    topn: Number(topnSelect?.value || 10),
  };

  const buildSegments = (segments) => {
    const topn = Math.max(1, state.topn || 10);
    const topItems = segments.slice(0, topn);
    const otherItems = segments.slice(topn);
    const merged = topItems.map((segment, index) => ({
      ...segment,
      color: segment.color || ["#526a86", "#76834f", "#8a6654", "#6d6481", "#4e7b75", "#9b7a3f", "#8a4d47", "#5d7c83", "#7b6d5a", "#697789"][index % 10] || '#526a86',
    }));
    if (otherItems.length) {
      merged.push({
        name: 'Other',
        ratio: otherItems.reduce((sum, item) => sum + (Number(item.ratio) || 0), 0),
        reads: otherItems.reduce((sum, item) => sum + (Number(item.reads) || 0), 0),
        color: '#8d8d8d',
      });
    }
    return merged;
  };

  const render = () => {
    container.classList.remove("empty-box");
    if (tag) tag.textContent = `Top ${state.topn} + Other`;
    container.innerHTML = `
      <div class="abundance-stack-list">
        ${ranks.map((rankBlock) => {
          const segments = buildSegments(Array.isArray(rankBlock.segments) ? rankBlock.segments : []);
          return `
            <article class="abundance-rank-card">
              <div class="abundance-rank-head">
                <div>
                  <span class="abundance-rank-label">${escapeHtml(rankBlock.rank)}</span>
                  <strong>分类丰度堆积</strong>
                </div>
                <span class="abundance-rank-total">总占比 ${escapeHtml(formatRate(rankBlock.total_ratio))}</span>
              </div>
              <div class="abundance-bar" role="img" aria-label="${escapeHtml(rankBlock.rank)}分类丰度堆积图">
                ${segments.map((segment) => `
                  <span
                    class="abundance-segment"
                    style="width:${Math.max(Number(segment.ratio) || 0, 0)}%; --segment-color:${escapeHtml(segment.color || '#526a86')}"
                    title="${escapeHtml(`${segment.name} | ${formatRate(segment.ratio)} | ${segment.reads} reads`)}"
                  ></span>
                `).join("")}
              </div>
              <div class="abundance-legend">
                ${segments.map((segment) => `
                  <div class="abundance-legend-item">
                    <i style="--legend-color:${escapeHtml(segment.color || '#526a86')}"></i>
                    <span class="abundance-legend-name">${escapeHtml(segment.name)}</span>
                    <strong>${escapeHtml(formatRate(segment.ratio))}</strong>
                  </div>
                `).join("")}
              </div>
            </article>
          `;
        }).join("")}
      </div>
    `;
  };

  if (topnSelect && topnSelect.dataset.bound !== 'true') {
    topnSelect.dataset.bound = 'true';
    topnSelect.addEventListener('change', () => {
      state.topn = Number(topnSelect.value || 10);
      render();
    });
  }

  render();
}

function renderMlstSection(section) {
  const tableContainer = document.getElementById("mlst-table");
  const detailContainer = document.getElementById("mlst-gene-show");
  if (!tableContainer || !detailContainer) return;

  const columns = Array.isArray(section?.columns) ? section.columns : [];
  const rows = Array.isArray(section?.rows) ? section.rows : [];
  const geneShowMap = section?.gene_show_map || {};
  const hostGeneIdIndex = columns.indexOf("Host Gene ID");
  const hostGeneDisplayIndex = columns.indexOf("Host Gene 展示");
  const defaultGene = section?.default_gene || (hostGeneIdIndex >= 0 ? (rows[0]?.[hostGeneIdIndex] ?? "") : "");

  if (!rows.length || !columns.length) {
    buildTableCard("mlst-table", "MLST 分析结果", [], []);
    detailContainer.innerHTML = `
      <div class="empty-box">
        <p>未检出 MLST 结果文件或 host gene 比对文件。</p>
      </div>
    `;
    return;
  }

  const displayColumns = columns.filter((column) => column !== "Host Gene ID" && column !== "Host Gene 展示");
  tableContainer.dataset.exportTitle = "MLST分析结果";
  tableContainer.innerHTML = `
    ${renderTableExportToolbar()}
    <div class="table-frame">
      <table class="report-table mlst-report-table">
        <thead><tr>${displayColumns.map((column) => `<th>${escapeHtml(column)}</th>`).join("")}</tr></thead>
        <tbody>
          ${rows.map((row) => `
            <tr class="mlst-clickable-row" data-host-gene-id="${escapeHtml(hostGeneIdIndex >= 0 ? (row[hostGeneIdIndex] ?? "") : "")}">
              ${displayColumns.map((column) => {
                const index = columns.indexOf(column);
                const value = row[index] ?? "";
                if (column === "Host Gene") {
                  const displayValue = hostGeneDisplayIndex >= 0 ? (row[hostGeneDisplayIndex] ?? value) : value;
                  return `<td><button class="mlst-hostgene-button" type="button" data-host-gene-id="${escapeHtml(hostGeneIdIndex >= 0 ? (row[hostGeneIdIndex] ?? "") : "")}">${escapeHtml(String(displayValue))}</button></td>`;
                }
                return `<td>${renderTableCellContent(value)}</td>`;
              }).join("")}
            </tr>
          `).join("")}
        </tbody>
      </table>
    </div>
  `;
  bindTableExportButtons(
    tableContainer,
    "MLST分析结果",
    displayColumns,
    rows.map((row) => displayColumns.map((column) => row[columns.indexOf(column)] ?? "")),
  );

  const renderDetail = (hostGeneId) => {
    const detail = geneShowMap[hostGeneId] || "未检出对应的 host gene 比对详情。";
    detailContainer.classList.remove("empty-box");
    detailContainer.innerHTML = `
      <div class="gene-show-meta">当前 Host Gene：${escapeHtml(hostGeneId || "-")}</div>
      <pre class="gene-show-pre">${escapeHtml(detail)}</pre>
    `;
    tableContainer.querySelectorAll('tbody tr').forEach((rowNode) => {
      rowNode.classList.toggle('mlst-row-active', (rowNode.dataset.hostGeneId || '') === hostGeneId);
    });
  };

  tableContainer.querySelectorAll('[data-host-gene-id]').forEach((node) => {
    node.addEventListener('click', () => renderDetail(node.dataset.hostGeneId || ''));
  });
  renderDetail(defaultGene);
}

function renderAssemblyCoverage(section) {
  const container = document.getElementById("assembly-coverage-chart");
  if (!container) return;
  const points = Array.isArray(section?.points) ? section.points : [];
  if (!points.length) {
    container.innerHTML = `
      <div class="empty-box coverage-empty">
        <p>未检出覆盖度结果文件。</p>
      </div>
    `;
    return;
  }
  const chartSvg = createSeriesSvg(
    [{ label: "测序深度", color: "#4e6177", values: points }],
    {
      label: section.label || "基因组覆盖度",
      width: 1120,
      height: 420,
      padX: 44,
      padTop: 10,
      padBottom: 42,
      xLabel: section.x_label || "基因组位置",
      yLabel: section.y_label || "测序深度",
      xTicks: section.x_ticks || [1, Math.max(1, Math.round(points.length / 2)), points.length],
      xValues: section.x_values || [],
    },
  );
  container.classList.remove("empty-box");
  container.innerHTML = `
    <div class="coverage-summary-bar">
      <div class="coverage-summary-item">
        <span>总位点</span>
        <strong>${escapeHtml(String(section.total_bases ?? "-"))}</strong>
      </div>
      <div class="coverage-summary-item">
        <span>平均深度</span>
        <strong>${escapeHtml(String(section.mean_depth ?? "-"))}</strong>
      </div>
      <div class="coverage-summary-item">
        <span>最大深度</span>
        <strong>${escapeHtml(String(section.max_depth ?? "-"))}</strong>
      </div>
      <div class="coverage-summary-item">
        <span>Contig 数</span>
        <strong>${escapeHtml(String(section.contig_count ?? "-"))}</strong>
      </div>
    </div>
    <div class="mini-chart-card coverage-chart-card">
      <span class="mini-chart-title">${escapeHtml(section.label || "基因组覆盖度")}</span>
      ${buildChartInsight(summarizeCoverage(section))}
      ${chartSvg}
    </div>
  `;
}


function renderBarSvg(values, options = {}) {
  const width = options.width || 1100;
  const height = options.height || 360;
  const padX = options.padX ?? 60;
  const padTop = options.padTop ?? 32;
  const padBottom = options.padBottom ?? 88;
  const innerWidth = width - padX * 2;
  const innerHeight = height - padTop - padBottom;
  const maxValue = Math.max(...values.map((value) => Number(value) || 0), 1);
  const yTicks = 4;
  const valueLabelStep = values.length <= 12 ? 1 : Math.ceil(values.length / 8);
  const grid = Array.from({ length: yTicks + 1 }, (_, index) => {
    const ratio = index / yTicks;
    const y = padTop + innerHeight - innerHeight * ratio;
    const value = maxValue * ratio;
    return `
      <line class="chart-grid-line" x1="${padX}" y1="${y}" x2="${width - padX}" y2="${y}"></line>
      <text class="chart-axis-label y-axis" x="${padX - 10}" y="${y + 4}">${escapeHtml(formatChartValue(value))}</text>
    `;
  }).join("");
  const barGap = 8;
  const barWidth = Math.max((innerWidth - Math.max(values.length - 1, 0) * barGap) / Math.max(values.length, 1), 10);
  const bars = values.map((value, index) => {
    const numeric = Number(value) || 0;
    const barHeight = (numeric / Math.max(maxValue, 1)) * innerHeight;
    const x = padX + index * (barWidth + barGap);
    const y = padTop + innerHeight - barHeight;
    const showValueLabel = index % valueLabelStep === 0 || index === values.length - 1;
    return `
      <rect class="chart-bar" x="${x}" y="${y}" width="${barWidth}" height="${barHeight}" rx="6" ry="6">
        <title>${escapeHtml(`${options.xValues?.[index] ?? index + 1}: ${numeric}`)}</title>
      </rect>
      ${showValueLabel ? `<text class="chart-bar-value" x="${x + barWidth / 2}" y="${Math.max(y - 10, padTop - 8)}">${escapeHtml(String(numeric))}</text>` : ""}
      <text class="chart-axis-label x-axis chart-rotated-label" x="${x + barWidth / 2}" y="${height - 30}">${escapeHtml(String(options.xValues?.[index] ?? index + 1))}</text>
    `;
  }).join("");
  return `
    <div class="mini-chart-card coverage-chart-card">
      <span class="mini-chart-title">${escapeHtml(options.label || "柱状分布图")}</span>
      <div class="chart-canvas" style="--chart-height:${height}px">
        <svg class="sparkline-svg bar-chart-svg" viewBox="0 0 ${width} ${height}" preserveAspectRatio="xMidYMid meet" role="img" aria-label="${escapeHtml(options.label || "柱状分布图")}">
          ${grid}
          <line class="chart-axis-line" x1="${padX}" y1="${padTop}" x2="${padX}" y2="${height - padBottom}"></line>
          <line class="chart-axis-line" x1="${padX}" y1="${height - padBottom}" x2="${width - padX}" y2="${height - padBottom}"></line>
          ${bars}
          <text class="chart-axis-title" x="${width / 2}" y="${height - 8}">${escapeHtml(options.xLabel || "分类")}</text>
          <text class="chart-axis-title chart-axis-title-y" x="24" y="${height / 2}">${escapeHtml(options.yLabel || "数量")}</text>
        </svg>
      </div>
    </div>
  `;
}


function renderGeneLengthDistribution(section) {
  const container = document.getElementById("gene-length-distribution-chart");
  if (!container) return;
  const points = Array.isArray(section?.points) ? section.points : [];
  if (!points.length) {
    container.innerHTML = `
      <div class="empty-box coverage-empty">
        <p>未检出基因长度统计文件。</p>
      </div>
    `;
    return;
  }
  container.classList.remove("empty-box");
  container.innerHTML = `
    ${buildChartInsight(summarizeGeneLengthDistribution(section))}
    ${renderBarSvg(points, {
    label: section.label || "基因长度与数量分布",
    width: 1120,
    height: 420,
    padX: 48,
    padTop: 24,
    padBottom: 96,
    xLabel: section.x_label || "基因长度范围",
    yLabel: section.y_label || "Gene数量",
    xValues: section.x_values || [],
  })}
  `;
}


function renderContigDepthRelationship(section) {
  const container = document.getElementById("contig-depth-relationship-chart");
  if (!container) return;
  const points = Array.isArray(section?.points) ? section.points : [];
  if (!points.length) {
    container.innerHTML = `
      <div class="empty-box coverage-empty">
        <p>未检出可展示的序列类型与平均深度关系数据。</p>
      </div>
    `;
    return;
  }
  container.classList.remove("empty-box");
  const width = 980;
  const height = 320;
  const padLeft = 84;
  const padRight = 38;
  const padTop = 26;
  const padBottom = 58;
  const innerWidth = width - padLeft - padRight;
  const innerHeight = height - padTop - padBottom;
  const groups = ["基因组", "质粒"];
  const maxDepth = Math.max(...points.map((point) => Number(point.depth) || 0), 1);
  const roundedMax = Math.ceil(maxDepth * 1.08);
  const xPositions = { "基因组": padLeft + innerWidth * 0.28, "质粒": padLeft + innerWidth * 0.72 };
  const yTicks = 4;
  const grid = Array.from({ length: yTicks + 1 }, (_, index) => {
    const ratio = index / yTicks;
    const y = padTop + innerHeight - innerHeight * ratio;
    const value = roundedMax * ratio;
    return `
      <line class="chart-grid-line" x1="${padLeft}" y1="${y}" x2="${width - padRight}" y2="${y}"></line>
      <text class="chart-axis-label y-axis" x="${padLeft - 10}" y="${y + 4}">${escapeHtml(formatChartValue(value))}</text>
    `;
  }).join("");
  const dots = points.map((point, index) => {
    const baseX = xPositions[point.type] || xPositions["基因组"];
    const jitter = ((index % 11) - 5) * 7;
    const depth = Number(point.depth) || 0;
    const y = padTop + innerHeight - (depth / Math.max(roundedMax, 1)) * innerHeight;
    const fill = point.type === "质粒" ? "#a88262" : "#4e6177";
    return `<circle cx="${baseX + jitter}" cy="${y}" r="5.5" fill="${fill}" fill-opacity="0.82"><title>${escapeHtml(`${point.name} | ${point.type} | 平均深度 ${depth}`)}</title></circle>`;
  }).join("");
  const labels = groups.map((group) => `<text class="chart-axis-label x-axis" x="${xPositions[group]}" y="${height - 28}">${group}</text>`).join("");
  container.innerHTML = `
    <div class="mini-chart-card relation-chart-card">
      <span class="mini-chart-title">${escapeHtml(section.label || "基因组/质粒与平均深度关系图")}</span>
      ${buildChartInsight(summarizeContigDepth(section))}
      <div class="chart-canvas" style="--chart-height:${height}px">
        <svg class="sparkline-svg" viewBox="0 0 ${width} ${height}" preserveAspectRatio="xMidYMid meet" role="img" aria-label="${escapeHtml(section.label || "基因组/质粒与平均深度关系图")}">
          ${grid}
          <line class="chart-axis-line" x1="${padLeft}" y1="${padTop}" x2="${padLeft}" y2="${height - padBottom}"></line>
          <line class="chart-axis-line" x1="${padLeft}" y1="${height - padBottom}" x2="${width - padRight}" y2="${height - padBottom}"></line>
          ${labels}
          ${dots}
          <text class="chart-axis-title" x="${width / 2}" y="${height - 6}">${escapeHtml(section.x_label || "序列类型")}</text>
          <text class="chart-axis-title chart-axis-title-y" x="24" y="${height / 2}">${escapeHtml(section.y_label || "平均深度")}</text>
        </svg>
      </div>
    </div>
  `;
}

function summarizeContigLengthDepthScatter(section) {
  const points = Array.isArray(section?.points) ? section.points : [];
  if (!points.length) {
    return "未检出可用于长度与平均深度联合展示的 Contig 数据。";
  }
  const genome = points.filter((point) => point.type === "基因组");
  const plasmid = points.filter((point) => point.type === "质粒");
  const maxLengthPoint = points.reduce((best, current) => ((Number(current.length) || 0) > (Number(best.length) || 0) ? current : best), points[0]);
  const maxDepthPoint = points.reduce((best, current) => ((Number(current.depth) || 0) > (Number(best.depth) || 0) ? current : best), points[0]);
  const plasmidMeanDepth = plasmid.length
    ? plasmid.reduce((sum, item) => sum + (Number(item.depth) || 0), 0) / plasmid.length
    : 0;
  const genomeMeanDepth = genome.length
    ? genome.reduce((sum, item) => sum + (Number(item.depth) || 0), 0) / genome.length
    : 0;
  if (plasmid.length && genome.length && plasmidMeanDepth > genomeMeanDepth * 1.35) {
    return `质粒序列整体平均深度高于基因组背景，提示部分质粒可能具有更高拷贝特征；最长Contig为 ${maxLengthPoint?.name || "-"}。`;
  }
  return `Contig 长度与平均深度分布整体可见，最长Contig为 ${maxLengthPoint?.name || "-"}，最高深度序列为 ${maxDepthPoint?.name || "-"}。`;
}

function renderContigLengthDepthScatter(section) {
  const container = document.getElementById("contig-length-depth-scatter-chart");
  if (!container) return;
  const points = Array.isArray(section?.points) ? section.points : [];
  if (!points.length) {
    container.innerHTML = `
      <div class="empty-box coverage-empty">
        <p>未检出可展示的 Contig 长度与平均测序深度数据。</p>
      </div>
    `;
    return;
  }
  container.classList.remove("empty-box");
  const width = 920;
  const height = 320;
  const padLeft = 86;
  const padRight = 38;
  const padTop = 24;
  const padBottom = 66;
  const innerWidth = width - padLeft - padRight;
  const innerHeight = height - padTop - padBottom;
  const maxLength = Math.max(...points.map((point) => Number(point.length) || 0), 1);
  const maxDepth = Math.max(...points.map((point) => Number(point.depth) || 0), 1);
  const roundedLength = Math.ceil(maxLength * 1.04);
  const roundedDepth = Math.ceil(maxDepth * 1.08);
  const xTicks = 4;
  const yTicks = 4;
  const gridX = Array.from({ length: xTicks + 1 }, (_, index) => {
    const ratio = index / xTicks;
    const x = padLeft + innerWidth * ratio;
    const value = roundedLength * ratio;
    return `
      <line class="chart-grid-line" x1="${x}" y1="${padTop}" x2="${x}" y2="${height - padBottom}"></line>
      <text class="chart-axis-label x-axis" x="${x}" y="${height - 28}">${escapeHtml(formatChartValue(value))}</text>
    `;
  }).join("");
  const gridY = Array.from({ length: yTicks + 1 }, (_, index) => {
    const ratio = index / yTicks;
    const y = padTop + innerHeight - innerHeight * ratio;
    const value = roundedDepth * ratio;
    return `
      <line class="chart-grid-line" x1="${padLeft}" y1="${y}" x2="${width - padRight}" y2="${y}"></line>
      <text class="chart-axis-label y-axis" x="${padLeft - 10}" y="${y + 4}">${escapeHtml(formatChartValue(value))}</text>
    `;
  }).join("");
  const dots = points.map((point) => {
    const length = Number(point.length) || 0;
    const depth = Number(point.depth) || 0;
    const x = padLeft + (length / Math.max(roundedLength, 1)) * innerWidth;
    const y = padTop + innerHeight - (depth / Math.max(roundedDepth, 1)) * innerHeight;
    const fill = point.type === "质粒" ? "#a88262" : "#4e6177";
    return `<circle cx="${x}" cy="${y}" r="5" fill="${fill}" fill-opacity="0.82"><title>${escapeHtml(`${point.name} | ${point.type} | 长度 ${formatChartValue(length)} bp | 平均深度 ${depth}`)}</title></circle>`;
  }).join("");
  container.innerHTML = `
    <div class="mini-chart-card relation-chart-card">
      <span class="mini-chart-title">${escapeHtml(section.label || "Contig长度与平均测序深度散点图")}</span>
      ${buildChartInsight(summarizeContigLengthDepthScatter(section))}
      <div class="chart-canvas" style="--chart-height:${height}px">
        <svg class="sparkline-svg" viewBox="0 0 ${width} ${height}" preserveAspectRatio="xMidYMid meet" role="img" aria-label="${escapeHtml(section.label || "Contig长度与平均测序深度散点图")}">
          ${gridX}
          ${gridY}
          <line class="chart-axis-line" x1="${padLeft}" y1="${padTop}" x2="${padLeft}" y2="${height - padBottom}"></line>
          <line class="chart-axis-line" x1="${padLeft}" y1="${height - padBottom}" x2="${width - padRight}" y2="${height - padBottom}"></line>
          ${dots}
          <text class="chart-axis-title" x="${width / 2}" y="${height - 6}">${escapeHtml(section.x_label || "Contig长度(bp)")}</text>
          <text class="chart-axis-title chart-axis-title-y" x="24" y="${height / 2}">${escapeHtml(section.y_label || "平均测序深度")}</text>
        </svg>
      </div>
    </div>
  `;
}

function renderCategoryGeneRelationship(containerId, section) {
  const container = document.getElementById(containerId);
  if (!container) return;
  const leftNodes = Array.isArray(section?.nodes_left) ? section.nodes_left : [];
  const rightNodes = Array.isArray(section?.nodes_right) ? section.nodes_right : [];
  const links = Array.isArray(section?.links) ? section.links : [];
  if (!leftNodes.length || !rightNodes.length || !links.length) {
    container.innerHTML = `
      <div class="empty-box coverage-empty">
        <p>未检出可展示的关系数据。</p>
      </div>
    `;
    return;
  }
  container.classList.remove("empty-box");
  const width = 1040;
  const height = Math.max(300, Math.max(leftNodes.length, rightNodes.length) * 34 + 70);
  const leftX = 176;
  const rightX = width - 176;
  const leftYStep = leftNodes.length > 1 ? (height - 80) / (leftNodes.length - 1) : 0;
  const rightYStep = rightNodes.length > 1 ? (height - 80) / (rightNodes.length - 1) : 0;
  const leftPos = Object.fromEntries(leftNodes.map((node, index) => [node.name, 40 + index * leftYStep]));
  const rightPos = Object.fromEntries(rightNodes.map((node, index) => [node.name, 40 + index * rightYStep]));
  const maxLink = Math.max(...links.map((link) => Number(link.value) || 0), 1);
  const linkSvg = links.map((link) => {
    const y1 = leftPos[link.source];
    const y2 = rightPos[link.target];
    if (y1 === undefined || y2 === undefined) return '';
    const strokeWidth = 1.6 + ((Number(link.value) || 0) / maxLink) * 5.4;
    return `<path d="M ${leftX} ${y1} C ${leftX + 160} ${y1}, ${rightX - 160} ${y2}, ${rightX} ${y2}" class="relation-link" style="stroke-width:${strokeWidth}px"><title>${escapeHtml(`${link.source} -> ${link.target}: ${link.value}`)}</title></path>`;
  }).join('');
  const leftSvg = leftNodes.map((node) => `
    <g transform="translate(0 ${leftPos[node.name]})">
      <text class="relation-label relation-label-left" x="154" y="4">${escapeHtml(node.name)}</text>
      <rect class="relation-node relation-node-left" x="160" y="-10" width="16" height="20" rx="6"></rect>
      <text class="relation-value relation-value-left" x="150" y="-14">${escapeHtml(String(node.value))}</text>
    </g>
  `).join('');
  const rightSvg = rightNodes.map((node) => `
    <g transform="translate(0 ${rightPos[node.name]})">
      <rect class="relation-node relation-node-right" x="${rightX}" y="-10" width="16" height="20" rx="6"></rect>
      <text class="relation-label relation-label-right" x="${rightX + 22}" y="4">${escapeHtml(node.name)}</text>
      <text class="relation-value relation-value-right" x="${rightX + 22}" y="-14">${escapeHtml(String(node.value))}</text>
    </g>
  `).join('');
  container.innerHTML = `
    <div class="mini-chart-card relation-chart-card">
      <span class="mini-chart-title">${escapeHtml(section.label || '关系图')}</span>
      ${buildChartInsight(summarizeRelationship(section))}
      <div class="relationship-chart-head">
        <span>${escapeHtml(section.left_label || '左侧分类')}</span>
        <span>${escapeHtml(section.right_label || '右侧基因')}</span>
      </div>
      <div class="chart-canvas" style="--chart-height:${height}px">
        <svg class="sparkline-svg relation-svg" viewBox="0 0 ${width} ${height}" preserveAspectRatio="xMidYMid meet" role="img" aria-label="${escapeHtml(section.label || '关系图')}">
          ${linkSvg}
          ${leftSvg}
          ${rightSvg}
        </svg>
      </div>
    </div>
  `;
}

function renderRawQc(sections) {
  const left = sections?.raw_qc?.paired_end?.left || {};
  const right = sections?.raw_qc?.paired_end?.right || {};
  const fastp = sections?.raw_qc?.fastp || {};
  renderStatsWithCharts("raw-qc-left", left, "R1");
  renderStatsWithCharts("raw-qc-right", right, "R2");
  const summaryItems = [
    { label: "测序模式", value: String(fastp.sequencing || "-") },
    { label: "过滤后 Reads", value: String(fastp.filtering_result?.passed_filter_reads ?? "-") },
    { label: "过短 Reads", value: String(fastp.filtering_result?.too_short_reads ?? "-") },
    { label: "N 过多 Reads", value: String(fastp.filtering_result?.too_many_N_reads ?? "-") },
    { label: "重复率", value: formatRate(fastp.duplication_rate) },
  ];
  const summaryBox = document.getElementById("fastp-summary");
  if (summaryBox) {
    summaryBox.classList.remove("empty-box");
    const [primaryItem, ...secondaryItems] = summaryItems;
    const passed = Number(fastp.filtering_result?.passed_filter_reads || 0);
    const tooShort = Number(fastp.filtering_result?.too_short_reads || 0);
    const tooManyN = Number(fastp.filtering_result?.too_many_N_reads || 0);
    const duplication = Number(fastp.duplication_rate || 0);
    const totalObserved = passed + tooShort + tooManyN;
    const lossRate = totalObserved ? (tooShort + tooManyN) / totalObserved : 0;
    const riskState = duplication >= 0.5 || lossRate >= 0.2 ? 'danger' : duplication >= 0.3 || lossRate >= 0.1 ? 'warning' : 'success';
    const riskLabel = riskState === 'danger' ? '需重点关注' : riskState === 'warning' ? '轻度风险' : '质控平稳';
    summaryBox.innerHTML = `
      <div class="qc-summary-layout">
        <article class="qc-summary-primary qc-summary-primary-${riskState}">
          <div class="qc-summary-primary-head">
            <span class="qc-summary-kicker">质控总览</span>
            <span class="qc-summary-badge qc-summary-badge-${riskState}">${escapeHtml(riskLabel)}</span>
          </div>
          <strong>${escapeHtml(primaryItem.value)}</strong>
          <p>${escapeHtml(primaryItem.label)}</p>
        </article>
        <div class="qc-summary-grid">
          ${secondaryItems.map((item) => `
            <article class="qc-summary-metric">
              <span>${escapeHtml(item.label)}</span>
              <strong>${escapeHtml(item.value)}</strong>
            </article>
          `).join("")}
        </div>
      </div>
    `;
  }
  renderFastpTabs(fastp);
  initializeBaseTabs();
  bindChartExportButtons();
}

async function loadReport() {
  const shell = document.querySelector(".report-shell");
  if (!shell) return;
  let data = window.__EMBEDDED_REPORT_DATA__ || currentReportData;
  if (!data) {
    const response = await fetch(shell.dataset.reportEndpoint, { credentials: "same-origin" });
    data = await response.json();
    if (!response.ok) throw new Error(data.error || "结果数据加载失败");
  }
  currentReportData = data;
  fillTaskMeta(data.task || {});
  renderExecutiveSummary(data);
  document.getElementById("overview-metrics").innerHTML = buildMetricCards(data.overview_metrics || []);
  renderRawQc(data.sections || {});
  renderTaxonomyRiskSummary(data.sections?.species_identification?.risk_summary || {});
  renderSpeciesIdentification(data.sections?.species_identification || {});
  renderTaxonomyAbundance(data.sections?.species_identification?.abundance || {});
  buildTableCard("assembly-species-table", "组装物种鉴定", data.sections?.species_identification?.assembly_taxonomy?.columns || [], data.sections?.species_identification?.assembly_taxonomy?.rows || []);
  buildTableCard("assembly-summary-table", "组装后信息统计", data.sections?.assembly?.summary?.columns || [], data.sections?.assembly?.summary?.rows || []);
  renderAssemblyCoverage(data.sections?.assembly?.coverage || {});
  buildTableCard("contig-annotation-table", "各个 Contig 注释结果", data.sections?.assembly?.contig_annotation?.columns || [], data.sections?.assembly?.contig_annotation?.rows || []);
  renderContigDepthRelationship(data.sections?.assembly?.contig_depth_relationship || {});
  renderContigLengthDepthScatter(data.sections?.assembly?.contig_depth_relationship?.length_depth_scatter || {});
  buildTableCard("checkm-table", "CheckM 统计结果", data.sections?.assembly?.checkm?.columns || [], data.sections?.assembly?.checkm?.rows || []);
  buildTableCard("gene-annotation-summary-table", "基因注释统计", data.sections?.assembly?.gene_annotation_summary?.columns || [], data.sections?.assembly?.gene_annotation_summary?.rows || []);
  renderGeneLengthDistribution(data.sections?.assembly?.gene_length_distribution || {});
  renderResistanceVirulenceOverview(data.sections?.resistance_virulence?.overview || {});
  buildTableCard("rv-summary-table", "耐药毒力结果汇总", data.sections?.resistance_virulence?.summary?.columns || [], data.sections?.resistance_virulence?.summary?.rows || []);
  buildTableCard("virulence-table", "毒力元件", data.sections?.resistance_virulence?.virulence_elements?.columns || [], data.sections?.resistance_virulence?.virulence_elements?.rows || []);
  renderCategoryGeneRelationship("virulence-relationship-chart", data.sections?.resistance_virulence?.virulence_relationship || {});
  buildTableCard("resistance-table", "耐药元件", data.sections?.resistance_virulence?.resistance_elements?.columns || [], data.sections?.resistance_virulence?.resistance_elements?.rows || []);
  renderCategoryGeneRelationship("resistance-relationship-chart", data.sections?.resistance_virulence?.resistance_relationship || {});
  renderMlstSection(data.sections?.mlst || {});
  buildTableCard("serotype-table", "血清型鉴定", data.sections?.serotype?.columns || [], data.sections?.serotype?.rows || []);
  buildTableCard("priority-serotype-table", "关注毒力血清型", data.sections?.priority_serotype?.columns || [], data.sections?.priority_serotype?.rows || []);
  initializeInteractiveCharts();
  bindChartExportButtons();
}

function initializeReportNav() {
  const groups = Array.from(document.querySelectorAll('[data-nav-group]'));
  const topLinks = Array.from(document.querySelectorAll('.report-nav > .report-nav-group > a.report-nav-link[href^="#"]'));
  const subLinks = Array.from(document.querySelectorAll('.report-subnav a[href^="#"]'));
  const toggles = Array.from(document.querySelectorAll('[data-nav-toggle]'));
  const sectionNodes = Array.from(document.querySelectorAll('.report-content .report-section[id]'));

  const openGroup = (group) => {
    const toggle = group.querySelector('[data-nav-toggle]');
    const subnav = group.querySelector('.report-subnav');
    if (!toggle || !subnav) return;
    group.classList.add('is-open');
    toggle.setAttribute('aria-expanded', 'true');
    subnav.hidden = false;
  };

  const closeGroup = (group) => {
    const toggle = group.querySelector('[data-nav-toggle]');
    const subnav = group.querySelector('.report-subnav');
    if (!toggle || !subnav) return;
    group.classList.remove('is-open');
    toggle.setAttribute('aria-expanded', 'false');
    subnav.hidden = true;
  };

  const clearActive = () => {
    [...topLinks, ...subLinks, ...toggles].forEach((node) => node.classList.remove('is-active'));
  };

  const syncActiveState = () => {
    const threshold = 168;
    let current = sectionNodes[0]?.id || '';
    sectionNodes.forEach((section) => {
      if (section.getBoundingClientRect().top <= threshold) {
        current = section.id;
      }
    });
    clearActive();
    let matchedLink = document.querySelector(`.report-subnav a[href="#${current}"]`);
    if (matchedLink) {
      const group = matchedLink.closest('[data-nav-group]');
      groups.forEach((item) => {
        if (item === group) openGroup(item);
        else closeGroup(item);
      });
      matchedLink.classList.add('is-active');
      group?.querySelector('[data-nav-toggle]')?.classList.add('is-active');
      return;
    }
    const topLink = document.querySelector(`.report-nav > .report-nav-group > a.report-nav-link[href="#${current}"]`);
    if (topLink) {
      groups.forEach((item) => closeGroup(item));
      topLink.classList.add('is-active');
      return;
    }
    const toggle = toggles.find((node) => node.dataset.navSection === current);
    if (toggle) {
      const group = toggle.closest('[data-nav-group]');
      groups.forEach((item) => {
        if (item === group) openGroup(item);
        else closeGroup(item);
      });
      toggle.classList.add('is-active');
    }
  };

  groups.forEach((group) => {
    const toggle = group.querySelector('[data-nav-toggle]');
    const subnav = group.querySelector('.report-subnav');
    if (!toggle || !subnav) return;
    closeGroup(group);
    toggle.addEventListener('click', () => {
      groups.forEach((item) => {
        if (item === group) openGroup(item);
        else closeGroup(item);
      });
      clearActive();
      toggle.classList.add('is-active');
      const targetId = toggle.dataset.navSection;
      if (targetId) {
        const section = document.getElementById(targetId);
        if (section) {
          section.scrollIntoView({ behavior: 'smooth', block: 'start' });
          history.replaceState(null, '', `#${targetId}`);
          window.setTimeout(syncActiveState, 80);
        }
      }
    });
  });

  [...topLinks, ...subLinks].forEach((link) => {
    link.addEventListener('click', () => {
      window.setTimeout(syncActiveState, 40);
    });
  });

  window.addEventListener('scroll', syncActiveState, { passive: true });
  window.addEventListener('hashchange', syncActiveState);
  syncActiveState();
}


document.addEventListener("DOMContentLoaded", () => {
  document.querySelectorAll(".report-tab-button").forEach((button) => {
    button.addEventListener("click", () => {
      const target = button.dataset.reportTab;
      document.querySelectorAll(".report-tab-button").forEach((item) => item.classList.toggle("active", item === button));
      document.querySelectorAll(".report-tab-panel").forEach((panel) => {
        panel.classList.toggle("active", panel.dataset.reportPanel === target);
      });
      initializeInteractiveCharts();
    });
  });
  document.querySelectorAll("[data-report-export-format]").forEach((button) => {
    button.addEventListener("click", async () => {
      try {
        await exportReportPage(button.dataset.reportExportFormat || "html");
      } catch (error) {
        console.error(error);
      }
    });
  });
  initializeReportNav();
  loadReport().catch((error) => {
    console.error(error);
  });
});
