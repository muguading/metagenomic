let runtime = null;
let initialized = false;

const elements = () => runtime.selectors.elements();
const state = () => runtime.selectors.state();
const escapeHtml = (value) => runtime.ui.escapeHtml(value);

export function init(nextRuntime) {
  runtime = nextRuntime;
  if (initialized) return;
  initialized = true;
  const dom = elements();
  dom.refreshAuditButton?.addEventListener("click", () => refresh());
  runtime.actions.bindImeSafeInput(dom.auditSearch, commitSearch);
  dom.auditUserFilter?.addEventListener("change", onFacetChange);
  dom.auditModuleFilter?.addEventListener("change", onFacetChange);
  dom.auditActionFilter?.addEventListener("change", onFacetChange);
  dom.auditOutcomeFilter?.addEventListener("change", onFacetChange);
  dom.auditList?.addEventListener("click", handleTableClick);
  dom.auditList?.addEventListener("input", handleTableInput);
}

export function activate() {
  return load();
}

export function refresh() {
  return load(true);
}

async function load(force = false) {
  const appState = state();
  if (appState.currentUser?.role !== "admin") return;
  const data = await runtime.api.requestJson("/api/admin/audit-logs");
  appState.auditTrail.items = Array.isArray(data.items) ? data.items : [];
  appState.auditTrail.summary = data.summary || null;
  appState.auditTrail.facets = data.facets || { users: [], modules: [], actions: [] };
  renderFilters(force);
  renderSummary();
  renderLogs();
}

function renderFilters(force = false) {
  const appState = state();
  const dom = elements();
  const facets = appState.auditTrail.facets || {};
  const renderOptions = (items, value, emptyLabel) => [
    `<option value="">${escapeHtml(emptyLabel)}</option>`,
    ...(Array.isArray(items) ? items : []).map((item) => `<option value="${escapeHtml(item)}" ${item === value ? "selected" : ""}>${escapeHtml(item)}</option>`),
  ].join("");
  if (dom.auditUserFilter && (force || dom.auditUserFilter.options.length <= 1)) {
    dom.auditUserFilter.innerHTML = renderOptions(facets.users, appState.auditTrail.filters.username, "全部用户");
  }
  if (dom.auditModuleFilter && (force || dom.auditModuleFilter.options.length <= 1)) {
    dom.auditModuleFilter.innerHTML = renderOptions(facets.modules, appState.auditTrail.filters.module, "全部模块");
  }
  if (dom.auditActionFilter && (force || dom.auditActionFilter.options.length <= 1)) {
    dom.auditActionFilter.innerHTML = renderOptions(facets.actions, appState.auditTrail.filters.action, "全部动作");
  }
}

function renderSummary() {
  const dom = elements();
  if (!dom.auditSummary) return;
  const summary = state().auditTrail.summary || { total: 0, failed: 0, today: 0, users: 0 };
  const cards = [
    { label: "最近事件", value: summary.total, note: "当前加载到的审计记录数" },
    { label: "今日事件", value: summary.today, note: "当天产生的非查看行为" },
    { label: "失败事件", value: summary.failed, note: "返回 4xx / 5xx 的操作" },
    { label: "活跃用户", value: summary.users, note: "产生审计事件的用户数" },
  ];
  dom.auditSummary.innerHTML = cards.map((card) => `
    <article class="server-metric-card audit-summary-card">
      <span>${escapeHtml(card.label)}</span>
      <strong>${escapeHtml(String(card.value))}</strong>
      <small>${escapeHtml(card.note)}</small>
    </article>
  `).join("");
}

function filteredRows() {
  const auditTrail = state().auditTrail;
  const filters = auditTrail.filters || {};
  return (auditTrail.items || []).filter((item) => {
    const searchTarget = [item.path, item.target_id, item.request_summary, item.response_summary, item.username, item.module, item.action].join(" ").toLowerCase();
    if (filters.username && String(item.username || "") !== filters.username) return false;
    if (filters.module && String(item.module || "") !== filters.module) return false;
    if (filters.action && String(item.action || "") !== filters.action) return false;
    if (filters.outcome && String(item.outcome || "") !== filters.outcome) return false;
    if (filters.path && !String(item.path || "").toLowerCase().includes(String(filters.path).toLowerCase())) return false;
    if (filters.target_id && !String(item.target_id || "").toLowerCase().includes(String(filters.target_id).toLowerCase())) return false;
    return !filters.search || searchTarget.includes(String(filters.search).toLowerCase());
  });
}

function sortedRows(rows) {
  const auditTrail = state().auditTrail;
  const key = auditTrail.sortKey || "created_at";
  const direction = auditTrail.sortDirection === "asc" ? 1 : -1;
  return [...rows].sort((left, right) => String(left?.[key] || "").localeCompare(String(right?.[key] || ""), "zh-CN") * direction);
}

function renderLogs(focusKey = "", caretPosition = null) {
  const appState = state();
  const dom = elements();
  if (!dom.auditList) return;
  const columns = [
    { key: "created_at", label: "时间" },
    { key: "username", label: "用户" },
    { key: "module", label: "模块" },
    { key: "action", label: "动作" },
    { key: "target_id", label: "对象编号" },
    { key: "path", label: "路径" },
    { key: "outcome", label: "结果" },
  ];
  const rows = sortedRows(filteredRows());
  dom.auditList.innerHTML = `
    <div class="database-table-frame">
      <table class="database-table report-table table-tone-assembly">
        <thead><tr>
          ${columns.map((column) => `
            <th><div class="table-head-stack queue-table-head-stack">
              <button class="database-sort-button table-sort-button ${appState.auditTrail.sortKey === column.key ? "active" : ""}" type="button" data-audit-sort="${escapeHtml(column.key)}">
                <span class="queue-table-head-label">${escapeHtml(column.label)}</span>
                <span>${appState.auditTrail.sortKey === column.key ? (appState.auditTrail.sortDirection === "asc" ? "↑" : "↓") : "↕"}</span>
              </button>
              ${runtime.ui.renderTableFilterInput({
                value: appState.auditTrail.filters[column.key] || "",
                dataName: "data-audit-filter",
                dataValue: column.key,
                clearDataName: "data-clear-audit-filter",
              })}
            </div></th>
          `).join("")}
          <th><span class="queue-table-head-label">请求摘要</span></th>
        </tr></thead>
        <tbody>
          ${rows.length ? rows.map((row) => `
            <tr>
              <td${runtime.ui.renderMobileCellAttributes("时间")}>${escapeHtml(runtime.ui.formatDate(row.created_at) || row.created_at || "-")}</td>
              <td${runtime.ui.renderMobileCellAttributes("用户")}><strong>${escapeHtml(row.username || "-")}</strong></td>
              <td${runtime.ui.renderMobileCellAttributes("模块")}>${escapeHtml(row.module || "-")}</td>
              <td${runtime.ui.renderMobileCellAttributes("动作")}>${escapeHtml(row.action || "-")}</td>
              <td title="${escapeHtml(row.target_id || "-")}"${runtime.ui.renderMobileCellAttributes("对象编号")}>${escapeHtml(runtime.ui.truncateText(row.target_id || "-", 20))}</td>
              <td title="${escapeHtml(row.path || "-")}"${runtime.ui.renderMobileCellAttributes("路径")}>${escapeHtml(runtime.ui.truncateText(row.path || "-", 36))}</td>
              <td${runtime.ui.renderMobileCellAttributes("结果")}><span class="status-chip ${row.outcome === "failed" ? "failed" : "running"}">${row.outcome === "failed" ? "失败" : "成功"}</span></td>
              <td title="${escapeHtml(row.request_summary || row.response_summary || "-")}"${runtime.ui.renderMobileCellAttributes("请求摘要")}>${escapeHtml(runtime.ui.truncateText(row.request_summary || row.response_summary || "-", 54))}</td>
            </tr>
          `).join("") : `<tr><td colspan="8" class="database-empty-cell">${runtime.ui.renderEmptyState({
            title: "当前没有符合条件的审计事件",
            reason: appState.auditTrail.items.length ? "审计记录存在，但当前搜索、用户、模块、动作或结果筛选没有命中。" : "系统暂时还没有记录到可展示的操作事件。",
            action: appState.auditTrail.items.length ? "放宽上方筛选条件，或清空表头列筛选后再查看。" : "完成一次登录、提交、导入或管理操作后，这里会自动出现审计轨迹。",
            className: "queue-empty database-empty-state",
          })}</td></tr>`}
        </tbody>
      </table>
    </div>
  `;
  if (focusKey) {
    const target = dom.auditList.querySelector(`[data-audit-filter="${CSS.escape(focusKey)}"]`);
    if (target instanceof HTMLInputElement) {
      target.focus();
      const caret = caretPosition ?? target.value.length;
      target.setSelectionRange(caret, caret);
    }
  }
}

function handleTableClick(event) {
  const auditTrail = state().auditTrail;
  const clearButton = event.target.closest("[data-clear-audit-filter]");
  if (clearButton) {
    const key = String(clearButton.dataset.clearAuditFilter || "").trim();
    if (!key) return;
    auditTrail.filters[key] = "";
    renderLogs(key, 0);
    return;
  }
  const button = event.target.closest("[data-audit-sort]");
  if (!button) return;
  const key = String(button.dataset.auditSort || "").trim();
  if (!key) return;
  if (auditTrail.sortKey === key) {
    auditTrail.sortDirection = auditTrail.sortDirection === "asc" ? "desc" : "asc";
  } else {
    auditTrail.sortKey = key;
    auditTrail.sortDirection = key === "created_at" ? "desc" : "asc";
  }
  renderLogs();
}

function handleTableInput(event) {
  const input = event.target.closest("[data-audit-filter]");
  if (!(input instanceof HTMLInputElement) || event.isComposing || input.dataset.imeComposing === "1") return;
  const key = String(input.dataset.auditFilter || "").trim();
  if (!key) return;
  state().auditTrail.filters[key] = input.value;
  renderLogs(key, input.selectionStart ?? input.value.length);
}

function commitSearch(input = elements().auditSearch) {
  state().auditTrail.filters.search = String(input?.value || "");
  renderLogs();
}

function onFacetChange() {
  const dom = elements();
  const filters = state().auditTrail.filters;
  filters.username = dom.auditUserFilter?.value || "";
  filters.module = dom.auditModuleFilter?.value || "";
  filters.action = dom.auditActionFilter?.value || "";
  filters.outcome = dom.auditOutcomeFilter?.value || "";
  renderLogs();
}
