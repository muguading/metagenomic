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

function buildMetricCards(metrics) {
  return metrics.map((metric) => `
    <article class="metric-card">
      <span class="metric-label">${escapeHtml(metric.label)}</span>
      <strong>${escapeHtml(metric.display ?? (metric.value == null || metric.value === "" ? "--" : `${metric.value}${metric.unit ? ` ${metric.unit}` : ""}`))}</strong>
      <div class="metric-note">接口键名：${escapeHtml(metric.key)}</div>
    </article>
  `).join("");
}

function buildTableCard(containerId, title, columns, rows) {
  const container = document.getElementById(containerId);
  if (!container) return;
  if (!Array.isArray(rows) || rows.length === 0) {
    container.innerHTML = `
      <div class="empty-table-state">
        <strong>${escapeHtml(title)} 暂无数据</strong>
        <p class="empty-copy">当前模块已预留接口与表格结构，后续接入真实结果后会在这里自动展示。</p>
      </div>
    `;
    return;
  }

  container.innerHTML = `
    <div class="table-frame">
      <table class="report-table">
        <thead>
          <tr>${columns.map((column) => `<th>${escapeHtml(column)}</th>`).join("")}</tr>
        </thead>
        <tbody>
          ${rows.map((row) => `<tr>${columns.map((column, index) => {
            const value = Array.isArray(row) ? row[index] : row[column];
            return `<td>${escapeHtml(value ?? "-")}</td>`;
          }).join("")}</tr>`).join("")}
        </tbody>
      </table>
    </div>
  `;
}

function fillTaskMeta(task) {
  document.getElementById("report-task-name").textContent = task.name || task.id || "-";
  document.getElementById("report-task-meta").textContent = `任务编号：${task.id || "-"}`;
  document.getElementById("report-sample-title").textContent = task.name || task.id || "分析结果";
  document.getElementById("report-sample-copy").textContent = `创建时间：${formatDateTime(task.created_at)}；开始时间：${formatDateTime(task.started_at)}；结束时间：${formatDateTime(task.finished_at)}。`;
  document.getElementById("meta-owner").textContent = task.owner || "-";
  document.getElementById("meta-group").textContent = task.group || "-";
  document.getElementById("meta-asm-type").textContent = task.asm_type || "-";
  document.getElementById("meta-method").textContent = task.method || "-";
  document.getElementById("meta-input").textContent = task.input_path || "-";
  document.getElementById("meta-output").textContent = task.output_dir || "-";
}

function renderKeyValueGrid(containerId, title, items) {
  const container = document.getElementById(containerId);
  if (!container) return;
  if (!items.length) {
    container.innerHTML = `<p class="empty-copy">${escapeHtml(title)} 暂无数据。</p>`;
    return;
  }
  container.innerHTML = `
    <div class="mini-stat-grid">
      ${items.map((item) => `
        <div class="mini-stat-card">
          <span>${escapeHtml(item.label)}</span>
          <strong>${escapeHtml(item.value)}</strong>
        </div>
      `).join("")}
    </div>
  `;
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

function renderRawQc(sections) {
  const left = sections?.raw_qc?.paired_end?.left || {};
  const right = sections?.raw_qc?.paired_end?.right || {};
  const fastp = sections?.raw_qc?.fastp || {};
  renderKeyValueGrid("raw-qc-left", "左端测序数据", [
    { label: "过滤前 reads", value: String(left.before_summary?.total_reads ?? "-") },
    { label: "过滤前数据量", value: formatBases(left.before_summary?.total_bases) },
    { label: "过滤后平均长度", value: String(left.after_summary?.mean_length ?? "-") },
    { label: "过滤后 Q30", value: formatRate(left.after_summary?.q30_rate) },
    { label: "过滤后 GC", value: formatRate(left.after_summary?.gc_content) },
  ]);
  renderKeyValueGrid("raw-qc-right", "右端测序数据", [
    { label: "过滤前 reads", value: String(right.before_summary?.total_reads ?? "-") },
    { label: "过滤前数据量", value: formatBases(right.before_summary?.total_bases) },
    { label: "过滤后平均长度", value: String(right.after_summary?.mean_length ?? "-") },
    { label: "过滤后 Q30", value: formatRate(right.after_summary?.q30_rate) },
    { label: "过滤后 GC", value: formatRate(right.after_summary?.gc_content) },
  ]);

  renderKeyValueGrid("fastp-summary", "fastp 质控摘要", [
    { label: "测序模式", value: String(fastp.sequencing || "-") },
    { label: "过滤后 reads", value: String(fastp.filtering_result?.passed_filter_reads ?? "-") },
    { label: "过短 reads", value: String(fastp.filtering_result?.too_short_reads ?? "-") },
    { label: "N 过多 reads", value: String(fastp.filtering_result?.too_many_N_reads ?? "-") },
    { label: "重复率", value: formatRate(fastp.duplication_rate) },
  ]);

  const visual = document.getElementById("fastp-visual");
  if (visual) {
    visual.innerHTML = `
      <div class="chart-placeholder-text">
        <strong>FASTP JSON 已接入</strong>
        <span>当前已读取 ${escapeHtml(String(fastp.filtering_result?.passed_filter_reads ?? "-"))} 条过滤后 reads</span>
        <span>后续可在这里接质量曲线、碱基组成曲线和插入片段分布图。</span>
      </div>
    `;
  }
}

async function loadReport() {
  const shell = document.querySelector(".report-shell");
  if (!shell) return;
  const endpoint = shell.dataset.reportEndpoint;
  if (!endpoint) return;
  const response = await fetch(endpoint, { credentials: "same-origin" });
  const data = await response.json();
  if (!response.ok) {
    throw new Error(data.error || "结果数据加载失败");
  }

  fillTaskMeta(data.task || {});
  document.getElementById("overview-metrics").innerHTML = buildMetricCards(data.overview_metrics || []);
  renderRawQc(data.sections || {});
  buildTableCard("species-table", "物种鉴定", data.sections?.species_identification?.columns || [], data.sections?.species_identification?.rows || []);
  buildTableCard("assembly-summary-table", "组装后信息统计", data.sections?.assembly?.summary?.columns || [], data.sections?.assembly?.summary?.rows || []);
  buildTableCard("contig-annotation-table", "各个 Contig 注释结果", data.sections?.assembly?.contig_annotation?.columns || [], data.sections?.assembly?.contig_annotation?.rows || []);
  buildTableCard("checkm-table", "CheckM 统计结果", data.sections?.assembly?.checkm?.columns || [], data.sections?.assembly?.checkm?.rows || []);
  buildTableCard("rv-summary-table", "耐药毒力结果汇总", data.sections?.resistance_virulence?.summary?.columns || [], data.sections?.resistance_virulence?.summary?.rows || []);
  buildTableCard("virulence-table", "毒力元件", data.sections?.resistance_virulence?.virulence_elements?.columns || [], data.sections?.resistance_virulence?.virulence_elements?.rows || []);
  buildTableCard("resistance-table", "耐药元件", data.sections?.resistance_virulence?.resistance_elements?.columns || [], data.sections?.resistance_virulence?.resistance_elements?.rows || []);
  buildTableCard("mlst-table", "MLST 分析结果", data.sections?.mlst?.columns || [], data.sections?.mlst?.rows || []);
  buildTableCard("serotype-table", "血清型鉴定", data.sections?.serotype?.columns || [], data.sections?.serotype?.rows || []);
  buildTableCard("priority-serotype-table", "关注毒力血清型", data.sections?.priority_serotype?.columns || [], data.sections?.priority_serotype?.rows || []);
}

document.addEventListener("DOMContentLoaded", () => {
  loadReport().catch((error) => {
    console.error(error);
  });
});
