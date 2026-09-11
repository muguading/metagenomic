import {
  clearProjectLocalState,
  hideAutoProject,
  readArchivedProjectKeys,
  readHiddenProjectKeys,
  readProjectDefinitions,
  readProjectMilestoneOverrides,
  readProjectOwnerOverrides,
  removeProjectDefinition,
  writeArchivedProjectKeys,
  writeHiddenProjectKeys,
  writeProjectDefinitions,
  writeProjectMilestoneOverrides,
  writeProjectOwnerOverrides,
} from "./storage.js";

let state = null;
let elements = null;
let clearScopedValidation;
let appendFieldValidationMessage;
let markFieldInvalid;
let withSubmittingState;
let getDatabaseMetadataDisplayValue;
let showModalElement;
let hideModalElement;
let setActiveTab;
let scheduleSegmentedControlsSync;
let openDatabaseSampleModal;
let parseLocationJson;
let escapeHtml;
let formatDate;
let formatDateTime;
let renderMobileCellAttributes;
let showToast;

export function createProjectManagement(runtime) {
  state = runtime.selectors.state();
  elements = runtime.selectors.elements();
  ({ clearScopedValidation, appendFieldValidationMessage, markFieldInvalid, withSubmittingState, getDatabaseMetadataDisplayValue, showModalElement, hideModalElement, setActiveTab, scheduleSegmentedControlsSync, openDatabaseSampleModal, parseLocationJson } = runtime.actions);
  ({ escapeHtml, formatDate, renderMobileCellAttributes, showToast } = runtime.ui);
  formatDateTime = formatDate;
  return Object.freeze({
    appendProjectMilestoneChildRow,
    closeProjectCreateModal,
    closeProjectDeleteModal,
    closeProjectMilestoneModal,
    closeProjectOwnerModal,
    handleOpenProjectCreateModal,
    onDeleteProjectMilestone,
    onSubmitProjectCreate,
    onSubmitProjectDelete,
    onSubmitProjectMilestone,
    onSubmitProjectOwner,
    renderProjectManagement,
  });
}

const PROJECT_STAGE_CONFIG = [
  { key: "intake", label: "样本入组", tone: "steel" },
  { key: "analysis", label: "分析处理中", tone: "cobalt" },
  { key: "review", label: "结果复核", tone: "amber" },
  { key: "delivery", label: "报告交付", tone: "teal" },
];

function getSelectedProjectCandidateRows() {
  const selectedKeys = new Set((state.databaseSelectedSamples || []).map((item) => String(item || "").trim()).filter(Boolean));
  if (!selectedKeys.size) return [];
  return (Array.isArray(state.databaseRecords) ? state.databaseRecords : []).filter((row) => {
    const sampleKey = String(row?.sample_key || "").trim();
    return sampleKey && selectedKeys.has(sampleKey) && String(row?.library_scope || "main") === "main";
  });
}

function buildDefaultProjectName(rows = []) {
  const first = rows[0] || {};
  const submittingUnit = getDatabaseMetadataDisplayValue(first, "submitting_unit") || first.sample_name || "未命名项目";
  const anchorDate = startOfDay(first.collection_date || first.imported_at || first.created_at || "") || new Date();
  return `${submittingUnit} ${anchorDate.getFullYear()}年${anchorDate.getMonth() + 1}月项目`;
}

function bindProjectCreatePanelScroll() {
  const panel = elements.projectCreateForm;
  if (!(panel instanceof HTMLElement) || panel.dataset.scrollBound === "1") return;
  panel.dataset.scrollBound = "1";
  panel.addEventListener("wheel", (event) => {
    const deltaY = Number(event.deltaY || 0);
    if (!deltaY) return;
    const maxScrollTop = panel.scrollHeight - panel.clientHeight;
    if (maxScrollTop <= 0) return;
    const nextScrollTop = Math.max(0, Math.min(maxScrollTop, panel.scrollTop + deltaY));
    if (nextScrollTop === panel.scrollTop) return;
    panel.scrollTop = nextScrollTop;
    event.preventDefault();
    event.stopPropagation();
  }, { passive: false });
}

function handleOpenProjectCreateModal(event) {
  event.preventDefault();
  event.stopPropagation();
  const dropdown = elements.openProjectCreateModalButton?.closest("details");
  openProjectCreateModal();
  window.setTimeout(() => {
    dropdown?.removeAttribute("open");
  }, 0);
}

function openProjectCreateModal() {
  clearScopedValidation(elements.projectCreateForm);
  const rows = getSelectedProjectCandidateRows();
  const names = rows.slice(0, 4).map((row) => row.sample_name || row.sample_key).filter(Boolean);
  const remainder = Math.max(0, rows.length - names.length);
  const submittingUnits = Array.from(new Set(rows.map((row) => String(getDatabaseMetadataDisplayValue(row, "submitting_unit") || "").trim()).filter(Boolean)));
  const leadUnit = submittingUnits[0] || "当前送检单位";
  const syndromes = Array.from(new Set(rows.map((row) => String(getDatabaseMetadataDisplayValue(row, "suspected_syndrome") || "").trim()).filter(Boolean)));
  const sampleDates = rows
    .map((row) => startOfDay(row.collection_date || row.imported_at || row.created_at || ""))
    .filter(Boolean)
    .sort((left, right) => left.getTime() - right.getTime());
  const dateRangeLabel = sampleDates.length
    ? `${formatDate(sampleDates[0])} - ${formatDate(sampleDates[sampleDates.length - 1])}`
    : "未记录采样日期";
  if (elements.projectCreateSampleCount) {
    elements.projectCreateSampleCount.textContent = `${rows.length} 份样本`;
  }
  if (elements.projectCreateSampleMetrics) {
    elements.projectCreateSampleMetrics.innerHTML = rows.length
      ? `
        <span class="project-create-metric-chip"><small>送检单位</small><strong>${escapeHtml(leadUnit)}</strong></span>
        <span class="project-create-metric-chip"><small>症候群</small><strong>${escapeHtml(syndromes[0] || "未标注")}</strong></span>
        <span class="project-create-metric-chip"><small>时间范围</small><strong>${escapeHtml(dateRangeLabel)}</strong></span>
      `
      : "";
  }
  if (elements.projectCreateSampleSummary) {
    if (!rows.length) {
      elements.projectCreateSampleSummary.textContent = "当前还没有勾选样本。请先在样本列表里勾选需要纳入项目的样本，再提交创建。";
    } else {
      const tail = remainder > 0 ? `其余 ${remainder} 份样本将一并并入该项目。` : "当前选中的样本将全部并入该项目。";
      elements.projectCreateSampleSummary.textContent = `建议优先确认项目名称与总负责人，后续可在项目管理中继续维护交付节点、负责人分工与项目甘特图。${tail}`;
    }
  }
  if (elements.projectCreateSamplePreview) {
    if (!rows.length) {
      elements.projectCreateSamplePreview.innerHTML = "";
    } else {
      elements.projectCreateSamplePreview.innerHTML = `
        <span class="project-create-preview-label">代表样本</span>
        <div class="project-create-preview-list">
          ${names.map((name) => `<span class="project-create-preview-chip">${escapeHtml(name)}</span>`).join("")}
          ${remainder > 0 ? `<span class="project-create-preview-chip is-muted">+${remainder} 份</span>` : ""}
        </div>
      `;
    }
  }
  if (elements.projectCreateNameInput instanceof HTMLInputElement) {
    elements.projectCreateNameInput.value = rows.length ? buildDefaultProjectName(rows) : "";
  }
  if (elements.projectCreateOwnerInput instanceof HTMLInputElement) {
    elements.projectCreateOwnerInput.value = "";
  }
  if (elements.submitProjectCreateButton instanceof HTMLButtonElement) {
    elements.submitProjectCreateButton.disabled = rows.length === 0;
  }
  bindProjectCreatePanelScroll();
  if (elements.projectCreateForm instanceof HTMLElement) {
    elements.projectCreateForm.scrollTop = 0;
  }
  showModalElement(elements.projectCreateModal);
}

function closeProjectCreateModal() {
  clearScopedValidation(elements.projectCreateForm);
  hideModalElement(elements.projectCreateModal);
  elements.projectCreateForm?.reset();
  if (elements.projectCreateSampleCount) {
    elements.projectCreateSampleCount.textContent = "0 份样本";
  }
  if (elements.submitProjectCreateButton instanceof HTMLButtonElement) {
    elements.submitProjectCreateButton.disabled = false;
  }
  if (elements.projectCreateSampleSummary) {
    elements.projectCreateSampleSummary.textContent = "当前尚未选择样本。";
  }
  if (elements.projectCreateSampleMetrics) {
    elements.projectCreateSampleMetrics.innerHTML = "";
  }
  if (elements.projectCreateSamplePreview) {
    elements.projectCreateSamplePreview.innerHTML = "";
  }
}

function performDeleteProject(projectKey = "") {
  const project = getProjectByKey(projectKey);
  if (!project) return false;
  if (project.explicitProject) {
    removeProjectDefinition(project.key);
  } else {
    hideAutoProject(project.key);
  }
  clearProjectLocalState(project.key);
  if (state.projectManagementSelectedKey === project.key) {
    state.projectManagementSelectedKey = "";
    state.projectManagementStageKey = "";
  }
  showToast(project.explicitProject ? "项目已删除" : "项目已从项目管理中移除");
  renderProjectManagement();
  return true;
}

function saveProjectArchivedState(projectKey = "", archived = true) {
  const normalizedKey = String(projectKey || "").trim();
  if (!normalizedKey) return;
  const archivedKeys = new Set(readArchivedProjectKeys());
  if (archived) archivedKeys.add(normalizedKey);
  else archivedKeys.delete(normalizedKey);
  writeArchivedProjectKeys(Array.from(archivedKeys));
}

function isProjectArchived(projectKey = "") {
  const normalizedKey = String(projectKey || "").trim();
  if (!normalizedKey) return false;
  return new Set(readArchivedProjectKeys()).has(normalizedKey);
}

function ensureProjectEditable(projectKey = "", actionLabel = "修改项目") {
  const project = getProjectByKey(projectKey);
  if (!project) {
    showToast("未找到当前项目。", true);
    return null;
  }
  if (project.archived) {
    showToast(`项目已归档，暂不支持${actionLabel}。`, true);
    return null;
  }
  return project;
}

function performArchiveProject(projectKey = "", archived = true) {
  const project = getProjectByKey(projectKey);
  if (!project) return false;
  saveProjectArchivedState(project.key, archived);
  if (archived && state.projectManagementSection === "gantt" && state.projectManagementSelectedKey === project.key) {
    state.projectManagementSection = "overview";
    state.projectManagementStageKey = "";
  }
  showToast(archived ? "项目已归档，已从甘特图中隐藏" : "项目已取消归档");
  renderProjectManagement();
  return true;
}

function openProjectDeleteModal(projectKey = "", mode = "delete") {
  const project = getProjectByKey(projectKey);
  if (!project) return;
  state.currentProjectDeleteKey = project.key;
  state.currentProjectDeleteMode = mode === "archive" ? "archive" : mode === "unarchive" ? "unarchive" : "delete";
  if (elements.projectDeleteProjectName) {
    elements.projectDeleteProjectName.textContent = project.name || "当前项目";
  }
  if (elements.projectDeleteTitle) {
    elements.projectDeleteTitle.textContent = state.currentProjectDeleteMode === "archive"
      ? "归档项目"
      : state.currentProjectDeleteMode === "unarchive"
        ? "取消归档"
        : "删除项目";
  }
  if (elements.projectDeleteSummary) {
    elements.projectDeleteSummary.textContent = state.currentProjectDeleteMode === "archive"
      ? "确认后项目将被归档，并从甘特图中隐藏。"
      : state.currentProjectDeleteMode === "unarchive"
        ? "确认后项目将恢复到可编辑状态并重新进入甘特图。"
        : "确认后将移除当前项目定义或将该自动聚合项目从项目管理中隐藏。";
  }
  if (elements.projectDeleteKicker) {
    elements.projectDeleteKicker.textContent = state.currentProjectDeleteMode === "archive"
      ? "归档操作"
      : state.currentProjectDeleteMode === "unarchive"
        ? "恢复操作"
        : "删除操作";
  }
  if (elements.projectDeleteActionLabel) {
    if (state.currentProjectDeleteMode === "archive") {
      elements.projectDeleteActionLabel.textContent = "归档项目";
    } else if (state.currentProjectDeleteMode === "unarchive") {
      elements.projectDeleteActionLabel.textContent = "取消归档";
      } else {
        elements.projectDeleteActionLabel.textContent = project.explicitProject ? "删除项目定义" : "移除自动聚合项目";
      }
  }
  if (elements.submitProjectDeleteButton) {
    elements.submitProjectDeleteButton.textContent = state.currentProjectDeleteMode === "archive"
      ? "确认归档"
      : state.currentProjectDeleteMode === "unarchive"
        ? "确认恢复"
        : "确认删除";
    elements.submitProjectDeleteButton.classList.toggle("danger", state.currentProjectDeleteMode !== "unarchive");
  }
  if (elements.projectDeleteNote) {
    if (state.currentProjectDeleteMode === "archive") {
      elements.projectDeleteNote.textContent = "确认后项目将被归档：不会再出现在甘特图中，且项目负责人、交付节点和阶段备注都将变为只读。";
    } else if (state.currentProjectDeleteMode === "unarchive") {
      elements.projectDeleteNote.textContent = "确认后项目将恢复到可编辑状态，并重新出现在项目甘特图中。";
    } else {
      elements.projectDeleteNote.textContent = project.explicitProject
        ? "确认后将从项目管理中删除该客户项目定义，项目相关的本地负责人、节点配置和交付状态会一起移除。"
        : "确认后将把当前自动聚合项目从项目管理中隐藏，后续仍可通过样本变化重新生成新的自动聚合项目。";
    }
  }
  showModalElement(elements.projectDeleteModal);
}

function closeProjectDeleteModal() {
  hideModalElement(elements.projectDeleteModal);
  state.currentProjectDeleteKey = "";
  state.currentProjectDeleteMode = "delete";
}

function onSubmitProjectDelete(event) {
  event.preventDefault();
  const projectKey = String(state.currentProjectDeleteKey || "").trim();
  const mode = String(state.currentProjectDeleteMode || "delete");
  if (!projectKey) {
    closeProjectDeleteModal();
    return;
  }
  let handled = false;
  if (mode === "archive") {
    handled = performArchiveProject(projectKey, true);
  } else if (mode === "unarchive") {
    handled = performArchiveProject(projectKey, false);
  } else {
    handled = performDeleteProject(projectKey);
  }
  closeProjectDeleteModal();
  if (!handled) {
    showToast("未找到当前项目，当前操作未执行。", true);
  }
}

async function onSubmitProjectCreate(event) {
  event.preventDefault();
  clearScopedValidation(elements.projectCreateForm);
  const rows = getSelectedProjectCandidateRows();
  if (!rows.length) {
    appendFieldValidationMessage(elements.projectCreateSampleSummary?.parentElement || elements.projectCreateForm, "project-create-samples", "请先在样本列表里勾选需要纳入项目的样本。");
    showToast("当前没有可用于创建项目的样本，请重新勾选。", "warning");
    return;
  }
  const projectName = String(elements.projectCreateNameInput?.value || "").trim();
  const ownerName = String(elements.projectCreateOwnerInput?.value || "").trim();
  if (!projectName) {
    markFieldInvalid(elements.projectCreateNameInput, "project-create-name", "请填写项目名称，便于后续在项目管理中识别。");
    elements.projectCreateNameInput?.focus();
    return;
  }
  await withSubmittingState(elements.submitProjectCreateButton, "创建中...", async () => {
    const sampleKeys = rows.map((row) => String(row.sample_key || "").trim()).filter(Boolean);
    const firstSyndrome = getDatabaseMetadataDisplayValue(rows[0], "suspected_syndrome") || "未标注症候群";
    const definitions = readProjectDefinitions();
    const key = `manual__${Date.now()}`;
    definitions.unshift({
      key,
      project_name: projectName,
      owner_name: ownerName,
      syndrome: firstSyndrome,
      sample_keys: sampleKeys,
      created_at: new Date().toISOString(),
      created_by: String(state.currentUser?.username || "admin"),
    });
    writeProjectDefinitions(definitions);
    const hiddenKeys = new Set(readHiddenProjectKeys());
    hiddenKeys.delete(key);
    writeHiddenProjectKeys(Array.from(hiddenKeys));
    saveProjectArchivedState(key, false);
    state.projectManagementSelectedKey = key;
    state.projectManagementStageKey = "";
    closeProjectCreateModal();
    setActiveTab("project-tab");
    renderProjectManagement();
    showToast(`已创建项目：${projectName}`);
  });
}

function getProjectOwnerState(projectKey = "", fallbackOwner = "") {
  const overrides = readProjectOwnerOverrides();
  const stored = overrides[String(projectKey || "")];
  const ownerName = String(stored?.ownerName || "").trim();
  return ownerName || String(fallbackOwner || "").trim() || "待分配";
}

function saveProjectOwnerState(projectKey = "", ownerName = "") {
  const normalizedKey = String(projectKey || "").trim();
  if (!normalizedKey) return;
  const overrides = readProjectOwnerOverrides();
  overrides[normalizedKey] = {
    ownerName: String(ownerName || "").trim(),
    modifiedAt: new Date().toISOString(),
    modifiedBy: String(state.currentUser?.username || "admin"),
  };
  writeProjectOwnerOverrides(overrides);
}

function buildProjectDefaultMilestones({ sampleCount, taskLinkedCount, reportReadyCount, latestLog, progress, ownerName }) {
  const milestoneOwner = String(ownerName || "").trim() || "待分配";
  return [
    {
      id: "intake",
      tone: "steel",
      label: "样本入组",
      owner: milestoneOwner,
      weight: 20,
      completion: sampleCount > 0 ? 100 : 0,
      note: "",
      detail: `${sampleCount} 份样本纳入项目`,
    },
    {
      id: "analysis",
      tone: "cobalt",
      label: "任务关联",
      owner: milestoneOwner,
      weight: 28,
      completion: sampleCount ? Math.round((taskLinkedCount / sampleCount) * 100) : 0,
      note: "",
      detail: `${taskLinkedCount}/${sampleCount} 份样本已关联分析任务`,
    },
    {
      id: "review",
      tone: "amber",
      label: "报告整理",
      owner: milestoneOwner,
      weight: 24,
      completion: sampleCount ? Math.round((reportReadyCount / sampleCount) * 100) : 0,
      note: "",
      detail: `${reportReadyCount}/${sampleCount} 份样本已生成报告目录`,
    },
    {
      id: "release",
      tone: "sage",
      label: "主库发布",
      owner: milestoneOwner,
      weight: 13,
      completion: latestLog ? 100 : 0,
      note: "",
      detail: latestLog ? `${latestLog.version_label || "已发布"} · ${formatDate(latestLog.created_at) || "-"}` : "尚未形成样本发布记录",
    },
    {
      id: "delivery",
      tone: "teal",
      label: "客户交付",
      owner: milestoneOwner,
      weight: 15,
      completion: Math.max(0, Math.min(100, progress)),
      note: "",
      detail: progress >= 90 ? "可进入客户交付复核" : "待样本分析与报告整理完成后交付",
    },
  ];
}

function normalizeProjectMilestone(item, fallbackIndex = 0) {
  const weightValue = Number(item?.weight);
  const completionValue = Number(item?.completion);
  const label = String(item?.label || "").trim() || `节点${fallbackIndex + 1}`;
  const tone = String(item?.tone || "").trim() || PROJECT_STAGE_CONFIG[fallbackIndex % PROJECT_STAGE_CONFIG.length]?.tone || "steel";
  const id = String(item?.id || `${label}-${fallbackIndex + 1}`).trim();
  const weight = Number.isFinite(weightValue) ? Math.max(1, Math.min(100, Math.round(weightValue))) : 20;
  const completion = Number.isFinite(completionValue) ? Math.max(0, Math.min(100, Math.round(completionValue))) : 0;
  const status = completion >= 100 ? "done" : completion > 0 ? "active" : "pending";
  const startDate = startOfDay(item?.startDate || item?.start_date || "");
  const endDate = startOfDay(item?.endDate || item?.end_date || "");
  const children = Array.isArray(item?.children)
    ? item.children
        .map((child, childIndex) => {
          const childCompletion = Number(child?.completion);
          const childLabel = String(child?.label || "").trim();
          if (!childLabel) return null;
          return {
            id: String(child?.id || `${id}-child-${childIndex + 1}`),
            label: childLabel,
            note: String(child?.note || "").trim(),
            completion: Number.isFinite(childCompletion) ? Math.max(0, Math.min(100, Math.round(childCompletion))) : 0,
          };
        })
        .filter(Boolean)
    : [];
  return {
    id,
    tone,
    label,
    owner: String(item?.owner || "").trim(),
    weight,
    completion,
    status,
    startDate: startDate ? startDate.toISOString() : "",
    endDate: endDate ? endDate.toISOString() : "",
    note: String(item?.note || "").trim(),
    detail: String(item?.detail || "").trim(),
    children,
  };
}

function renderProjectMilestoneChildrenRows(children = []) {
  if (!(elements.projectMilestoneChildren instanceof HTMLElement)) return;
  const normalizedChildren = Array.isArray(children) ? children : [];
  elements.projectMilestoneChildren.innerHTML = normalizedChildren.length
    ? normalizedChildren.map((child, index) => `
        <div class="project-milestone-child-row">
          <label>
            <span>子节点名称</span>
            <input data-project-milestone-child-label value="${escapeHtml(String(child?.label || ""))}" placeholder="例如：中试样本准备">
          </label>
          <label>
            <span>完成（%）</span>
            <input data-project-milestone-child-completion type="number" min="0" max="100" step="1" value="${escapeHtml(String(child?.completion ?? ""))}" placeholder="0">
          </label>
          <label class="wide">
            <span>备注</span>
            <input data-project-milestone-child-note value="${escapeHtml(String(child?.note || ""))}" placeholder="补充该步骤的说明、阻塞项或输出结果">
          </label>
          <button class="ghost-button compact danger" type="button" data-project-milestone-child-remove="${index}">删除</button>
        </div>
      `).join("")
    : `<div class="field-note">当前节点尚未拆分子节点，可点击“新增子节点”补充详细步骤。</div>`;
}

function appendProjectMilestoneChildRow(child = {}) {
  if (!(elements.projectMilestoneChildren instanceof HTMLElement)) return;
  const emptyNote = elements.projectMilestoneChildren.querySelector(".field-note");
  emptyNote?.remove();
  const row = document.createElement("div");
  row.className = "project-milestone-child-row";
  row.innerHTML = `
    <label>
      <span>子节点名称</span>
      <input data-project-milestone-child-label value="${escapeHtml(String(child?.label || ""))}" placeholder="例如：中试样本准备">
    </label>
    <label>
      <span>完成（%）</span>
      <input data-project-milestone-child-completion type="number" min="0" max="100" step="1" value="${escapeHtml(String(child?.completion ?? ""))}" placeholder="0">
    </label>
    <label class="wide">
      <span>备注</span>
      <input data-project-milestone-child-note value="${escapeHtml(String(child?.note || ""))}" placeholder="补充该步骤的说明、阻塞项或输出结果">
    </label>
    <button class="ghost-button compact danger" type="button" data-project-milestone-child-remove="1">删除</button>
  `;
  elements.projectMilestoneChildren.appendChild(row);
}

function getProjectMilestoneState(projectKey, defaults) {
  const overrides = readProjectMilestoneOverrides();
  const stored = overrides[String(projectKey || "")];
  if (!stored || !Array.isArray(stored.milestones)) {
    return defaults.map((item, index) => normalizeProjectMilestone(item, index));
  }
  const normalized = stored.milestones
    .map((item, index) => normalizeProjectMilestone(item, index))
    .filter((item) => String(item.label || "").trim());
  return normalized.length ? normalized : defaults.map((item, index) => normalizeProjectMilestone(item, index));
}

function saveProjectMilestoneState(projectKey, milestones) {
  const normalizedKey = String(projectKey || "").trim();
  if (!normalizedKey) return;
  const overrides = readProjectMilestoneOverrides();
  overrides[normalizedKey] = {
    milestones: (Array.isArray(milestones) ? milestones : []).map((item, index) => normalizeProjectMilestone(item, index)),
    modifiedAt: new Date().toISOString(),
    modifiedBy: String(state.currentUser?.username || "admin"),
  };
  writeProjectMilestoneOverrides(overrides);
}

function buildProjectStageSegments(startDate, totalDays, milestones) {
  const normalizedMilestones = (Array.isArray(milestones) ? milestones : []).map((item, index) => normalizeProjectMilestone(item, index));
  const totalWeight = Math.max(1, normalizedMilestones.reduce((sum, item) => sum + Math.max(1, Number(item.weight) || 0), 0));
  const safeTotalDays = Math.max(totalDays, normalizedMilestones.length || 1);
  let cursor = startOfDay(startDate) || startOfDay(new Date());
  const segments = normalizedMilestones.map((item, index) => {
    const explicitStart = startOfDay(item.startDate || "");
    const explicitEnd = startOfDay(item.endDate || "");
    if (explicitStart && explicitEnd && explicitEnd.getTime() >= explicitStart.getTime()) {
      const segmentStart = explicitStart.getTime() >= cursor.getTime() ? explicitStart : cursor;
      const segmentEnd = explicitEnd.getTime() >= segmentStart.getTime() ? explicitEnd : segmentStart;
      cursor = addDays(segmentEnd, 1);
      return {
        key: item.id,
        label: item.label,
        tone: item.tone,
        start: segmentStart,
        end: segmentEnd,
        progressEnd: addDays(segmentStart, Math.max(0, Math.ceil(((daysBetween(segmentStart, segmentEnd) + 1) * item.completion) / 100) - 1)),
        completion: item.completion,
      };
    }
    const remainingWeight = normalizedMilestones.slice(index).reduce((sum, row) => sum + Math.max(1, Number(row.weight) || 0), 0);
    const remainingDays = Math.max(1, safeTotalDays - index - normalizedMilestones.slice(0, index).reduce((sum, row) => sum + (row.durationDays || 0), 0));
    const suggestedDays = index === normalizedMilestones.length - 1
      ? remainingDays
      : Math.max(1, Math.round((safeTotalDays * item.weight) / totalWeight));
    const maxDays = index === normalizedMilestones.length - 1 ? remainingDays : Math.max(1, remainingDays - (normalizedMilestones.length - index - 1));
    const durationDays = Math.min(Math.max(1, suggestedDays), maxDays);
    item.durationDays = durationDays;
    const segmentStart = startOfDay(cursor) || startOfDay(new Date());
    const segmentEnd = addDays(segmentStart, durationDays - 1);
    cursor = addDays(segmentEnd, 1);
    return {
      key: item.id,
      label: item.label,
      tone: item.tone,
      start: segmentStart,
      end: segmentEnd,
      progressEnd: addDays(segmentStart, Math.max(0, Math.ceil((durationDays * item.completion) / 100) - 1)),
      completion: item.completion,
    };
  });
  return segments;
}

function formatProjectMonthLabel(value) {
  const date = value instanceof Date ? value : new Date(value);
  if (Number.isNaN(date.getTime())) return "";
  return `${date.getFullYear()}-${String(date.getMonth() + 1).padStart(2, "0")}`;
}

function formatProjectDateRange(start, end) {
  const startLabel = formatDate(start);
  const endLabel = formatDate(end);
  if (!startLabel && !endLabel) return "-";
  if (!startLabel) return endLabel;
  if (!endLabel || startLabel === endLabel) return startLabel;
  return `${startLabel} - ${endLabel}`;
}

function startOfDay(value) {
  const date = value instanceof Date ? new Date(value) : new Date(value);
  if (Number.isNaN(date.getTime())) return null;
  date.setHours(0, 0, 0, 0);
  return date;
}

function endOfDay(value) {
  const date = startOfDay(value);
  if (!date) return null;
  date.setHours(23, 59, 59, 999);
  return date;
}

function addDays(date, days) {
  const next = new Date(date);
  next.setDate(next.getDate() + days);
  return next;
}

function daysBetween(start, end) {
  const startDate = startOfDay(start);
  const endDate = startOfDay(end);
  if (!startDate || !endDate) return 0;
  return Math.max(0, Math.round((endDate.getTime() - startDate.getTime()) / 86400000));
}

function buildProjectTimelineMonths(projects) {
  if (!projects.length) return [];
  const start = startOfDay(new Date(Math.min(...projects.map((item) => item.timelineStart.getTime()))));
  const end = startOfDay(new Date(Math.max(...projects.map((item) => item.timelineEnd.getTime()))));
  if (!start || !end) return [];
  const cursor = new Date(start.getFullYear(), start.getMonth(), 1);
  const limit = new Date(end.getFullYear(), end.getMonth(), 1);
  const monthStarts = [];
  while (cursor.getTime() <= limit.getTime()) {
    monthStarts.push(new Date(cursor));
    cursor.setMonth(cursor.getMonth() + 1);
  }
  return monthStarts;
}

function buildProjectManagementProjectRecord(group, versionLogs, releaseVersions, options = {}) {
    const explicitName = String(options.projectName || "").trim();
    const explicitOwnerName = String(options.ownerName || "").trim();
    const archived = Boolean(options.archived);
    const sampleCount = group.rows.length;
    const dates = group.rows
      .map((item) => startOfDay(item.collection_date || item.imported_at || item.created_at || ""))
      .filter(Boolean)
      .sort((a, b) => a.getTime() - b.getTime());
    const startDate = dates[0] || startOfDay(new Date());
    const endDate = dates[dates.length - 1] || startDate;
    const typedCount = group.rows.filter((item) => String(item.mlst_st || "").trim() || String(item.serotype_result || "").trim()).length;
    const molecularCount = group.rows.filter((item) => (
      String(item.resistance_gene_hits || "").trim()
      || String(item.virulence_gene_hits || "").trim()
      || String(item.resistance_mge_hits || "").trim()
      || String(item.virulence_mge_hits || "").trim()
    )).length;
    const completedCount = group.rows.filter((item) => String(item.metadata_completion_status || "").trim() === "complete").length;
    const metadataRatio = sampleCount ? completedCount / sampleCount : 0;
    const typedRatio = sampleCount ? typedCount / sampleCount : 0;
    const molecularRatio = sampleCount ? molecularCount / sampleCount : 0;
    const inferredProgress = Math.max(
      12,
      Math.min(
        100,
        Math.round((metadataRatio * 0.35 + typedRatio * 0.3 + molecularRatio * 0.2 + 0.15) * 100),
      ),
    );
    const area = group.rows
      .map((item) => {
        const location = parseLocationJson(item.location_json);
        return [location.province, location.city].filter(Boolean).join(" / ");
      })
      .filter(Boolean)
      .sort((a, b) => a.localeCompare(b))[0] || "区域待补充";
    const taskNames = Array.from(new Set(group.rows.map((item) => String(item.task_name || "").trim()).filter(Boolean)));
    const taskLinkedCount = group.rows.filter((item) => String(item.task_id || "").trim() || String(item.task_name || "").trim()).length;
    const reportReadyCount = group.rows.filter((item) => String(item.report_dir || "").trim()).length;
    const sampleKeySet = new Set(group.rows.map((item) => String(item.sample_key || "").trim()).filter(Boolean));
    const projectLogs = versionLogs.filter((item) => sampleKeySet.has(String(item.sample_key || "").trim()));
    const latestLog = projectLogs.length ? projectLogs[0] : null;
    const latestRelease = releaseVersions[0] || null;
    const inferredOwnerName = group.syndrome.includes("脑膜")
      ? "侵袭感染项目组"
      : group.syndrome.includes("腹泻") || group.syndrome.includes("肠道")
        ? "肠道病项目组"
        : group.syndrome.includes("环境")
          ? "环境监测项目组"
          : group.syndrome.includes("医院感染")
            ? "院感项目组"
            : "综合项目组";
    const ownerName = getProjectOwnerState(group.key, explicitOwnerName || inferredOwnerName);
    const deliveryEnd = addDays(endDate, 4);
    const totalDurationDays = Math.max(7, daysBetween(startDate, deliveryEnd) + 1);
    const blockerParts = [];
    if (completedCount < sampleCount) blockerParts.push("主档补录未完成");
    if (taskLinkedCount < sampleCount) blockerParts.push("部分样本未关联任务");
    if (reportReadyCount < sampleCount) blockerParts.push("部分样本尚未形成报告");
    if (!latestLog) blockerParts.push("未形成主库发布记录");
    const blockerText = blockerParts.length ? blockerParts.slice(0, 2).join("；") : "当前无明显阻塞项";
    const defaultMilestones = buildProjectDefaultMilestones({
      sampleCount,
      taskLinkedCount,
      reportReadyCount,
      latestLog,
      progress: inferredProgress,
      ownerName,
    });
    const milestones = getProjectMilestoneState(group.key, defaultMilestones);
    const milestoneWeightSum = Math.max(1, milestones.reduce((sum, item) => sum + Math.max(1, Number(item.weight) || 0), 0));
    const progress = Math.max(
      0,
      Math.min(
        100,
        Math.round(milestones.reduce((sum, item) => sum + ((Math.max(1, Number(item.weight) || 0) / milestoneWeightSum) * Math.max(0, Math.min(100, Number(item.completion) || 0))), 0)),
      ),
    );
    const currentMilestone = milestones.find((item) => item.completion < 100) || milestones[milestones.length - 1];
    const stageLabel = progress >= 90
      ? "交付复核"
      : currentMilestone?.label || "样本入组";
    const customerStatus = progress >= 90
      ? "待客户确认"
      : progress >= 68
        ? "交付准备中"
        : progress >= 42
          ? "分析推进中"
          : "样本收集中";
    const stageSegments = buildProjectStageSegments(startDate, totalDurationDays, milestones);
    const timelineEnd = stageSegments.length ? stageSegments[stageSegments.length - 1].end : deliveryEnd;
    const samples = group.rows
      .slice()
      .sort((a, b) => String(a.collection_date || a.imported_at || "").localeCompare(String(b.collection_date || b.imported_at || "")))
      .map((item) => {
        const itemLogs = projectLogs.filter((entry) => String(entry.sample_key || "").trim() === String(item.sample_key || "").trim());
        const itemLatestLog = itemLogs[0] || null;
        return {
          sampleKey: String(item.sample_key || ""),
          sampleName: String(item.sample_name || item.sample_key || "-"),
          date: item.collection_date || item.imported_at || "",
          species: item.species_name || item.mlst_species_name || "-",
          taskName: item.task_name || "-",
          mlst: item.mlst_st || "-",
          serotype: item.serotype_result || "-",
          releaseLabel: itemLatestLog?.version_label || latestRelease?.version_label || "未发布",
          taskStatus: String(item.report_dir || "").trim()
            ? "报告已生成"
            : (String(item.task_id || "").trim() || String(item.task_name || "").trim())
              ? "分析处理中"
              : "待关联任务",
        };
      });
    return {
      key: group.key,
      name: explicitName || `${group.customer} · ${group.monthLabel}`,
      shortName: explicitName || group.customer,
      monthLabel: group.monthLabel,
      syndrome: group.syndrome,
      area,
      sampleCount,
      metadataRatio,
      typedRatio,
      molecularRatio,
      progress,
      stageLabel,
      ownerName,
      plannedDeliveryDate: timelineEnd,
      blockerText,
      customerStatus,
      taskNames,
      taskLinkedCount,
      reportReadyCount,
      projectLogs,
      latestLog,
      latestRelease,
      milestones,
      samples,
      startDate,
      endDate,
      timelineStart: startDate,
      timelineEnd,
      stageSegments,
      explicitProject: Boolean(options.explicitProject),
      archived,
    };
}

function buildProjectManagementProjects() {
  const rows = Array.isArray(state.databaseRecords) ? state.databaseRecords : [];
  const scopedRows = rows.filter((row) => String(row?.library_scope || "main") === "main");
  const versionLogs = Array.isArray(state.databaseVersionLogs) ? state.databaseVersionLogs : [];
  const releaseVersions = Array.isArray(state.databaseReleaseVersions) ? state.databaseReleaseVersions : [];
  const rowMap = new Map(scopedRows.map((row) => [String(row.sample_key || "").trim(), row]));
  const projectDefinitions = readProjectDefinitions();
  const hiddenProjectKeys = new Set(readHiddenProjectKeys());
  const archivedProjectKeys = new Set(readArchivedProjectKeys());
  const explicitSampleKeys = new Set();
  const explicitGroups = projectDefinitions
    .map((definition) => {
      const sampleKeys = Array.isArray(definition.sample_keys) ? definition.sample_keys.map((item) => String(item || "").trim()).filter(Boolean) : [];
      const groupRows = sampleKeys.map((key) => rowMap.get(key)).filter(Boolean);
      if (!groupRows.length) return null;
      sampleKeys.forEach((key) => explicitSampleKeys.add(key));
      const firstDate = startOfDay(groupRows[0]?.collection_date || groupRows[0]?.imported_at || "");
      const monthLabel = firstDate ? `${firstDate.getFullYear()}-${String(firstDate.getMonth() + 1).padStart(2, "0")}` : "未定时间";
      return {
        key: String(definition.key || "").trim(),
        customer: String(definition.project_name || "").trim() || "未命名项目",
        monthLabel,
        syndrome: String(definition.syndrome || getDatabaseMetadataDisplayValue(groupRows[0], "suspected_syndrome") || "未标注症候群").trim(),
        rows: groupRows,
        projectName: String(definition.project_name || "").trim(),
        ownerName: String(definition.owner_name || "").trim(),
        explicitProject: true,
        archived: archivedProjectKeys.has(String(definition.key || "").trim()),
      };
    })
    .filter(Boolean);
  const grouped = new Map();
  scopedRows.forEach((row) => {
    const sampleKey = String(row.sample_key || "").trim();
    if (explicitSampleKeys.has(sampleKey)) return;
    const submittingUnit = getDatabaseMetadataDisplayValue(row, "submitting_unit") || "未标注送检单位";
    const syndrome = getDatabaseMetadataDisplayValue(row, "suspected_syndrome") || "未标注症候群";
    const collectionDate = startOfDay(row.collection_date || row.imported_at || row.created_at || "");
    const monthLabel = collectionDate ? `${collectionDate.getFullYear()}-${String(collectionDate.getMonth() + 1).padStart(2, "0")}` : "未定时间";
    const projectKey = `${submittingUnit}__${monthLabel}`;
    if (hiddenProjectKeys.has(projectKey)) return;
    if (!grouped.has(projectKey)) {
      grouped.set(projectKey, {
        key: projectKey,
        customer: submittingUnit,
        monthLabel,
        syndrome,
        rows: [],
        archived: archivedProjectKeys.has(projectKey),
      });
    }
    grouped.get(projectKey).rows.push(row);
  });
  const projects = explicitGroups
    .concat(Array.from(grouped.values()))
    .map((group) => buildProjectManagementProjectRecord(group, versionLogs, releaseVersions, {
      projectName: group.projectName,
      ownerName: group.ownerName,
      explicitProject: group.explicitProject,
      archived: group.archived,
    }))
    .sort((a, b) => {
      const endDelta = b.timelineEnd.getTime() - a.timelineEnd.getTime();
      if (endDelta !== 0) return endDelta;
      return b.sampleCount - a.sampleCount;
    });
  return projects;
}

function formatProjectMilestoneStepScore(project, segment) {
  const milestones = Array.isArray(project?.milestones) ? project.milestones : [];
  const milestone = milestones.find((item) => String(item?.id || "") === String(segment?.key || "")) || null;
  const children = Array.isArray(milestone?.children) ? milestone.children : [];
  if (children.length) {
    const completedCount = children.filter((child) => Number(child?.completion || 0) >= 100).length;
    return `${completedCount}/${children.length}`;
  }
  return Number(segment?.completion || 0) >= 100 ? "1/1" : "0/1";
}

function formatStandaloneMilestoneStepScore(milestone) {
  const children = Array.isArray(milestone?.children) ? milestone.children : [];
  if (children.length) {
    const completedCount = children.filter((child) => Number(child?.completion || 0) >= 100).length;
    return `${completedCount}/${children.length}`;
  }
  return Number(milestone?.completion || 0) >= 100 ? "1/1" : "0/1";
}

function getProjectOverallScoreMeta(project) {
  const milestones = Array.isArray(project?.milestones) ? project.milestones : [];
  if (!milestones.length) {
    return { completedCount: 0, totalCount: 0, label: "0/0" };
  }
  const completedCount = milestones.filter((item) => Number(item?.completion || 0) >= 100).length;
  return {
    completedCount,
    totalCount: milestones.length,
    label: `${completedCount}/${milestones.length}`,
  };
}

const PROJECT_GANTT_ZOOM_MIN = 0.35;
const PROJECT_GANTT_ZOOM_MAX = 1.75;
const PROJECT_GANTT_ZOOM_STEP = 0.15;
const PROJECT_GANTT_BASE_DAY_WIDTH = 92;
const PROJECT_GANTT_META_WIDTH = 220;
const PROJECT_GANTT_TRACK_GAP = 14;
const PROJECT_GANTT_SHELL_PADDING_X = 24;
const PROJECT_GANTT_ROW_PADDING_X = 12;

function clampProjectGanttZoom(value) {
  const numeric = Number(value);
  if (!Number.isFinite(numeric)) return 1;
  return Math.max(PROJECT_GANTT_ZOOM_MIN, Math.min(PROJECT_GANTT_ZOOM_MAX, Math.round(numeric * 100) / 100));
}

function calculateProjectGanttFitZoom(totalDays, shellWidth) {
  const safeDays = Math.max(1, Number(totalDays) || 1);
  const availableTrackWidth = calculateProjectGanttFitTrackWidth(shellWidth);
  return clampProjectGanttZoom(availableTrackWidth / (safeDays * PROJECT_GANTT_BASE_DAY_WIDTH));
}

function calculateProjectGanttFitTrackWidth(shellWidth) {
  const safeShellWidth = Math.max(0, Number(shellWidth) || 0);
  return Math.max(
    240,
    safeShellWidth
      - PROJECT_GANTT_SHELL_PADDING_X
      - PROJECT_GANTT_ROW_PADDING_X
      - PROJECT_GANTT_META_WIDTH
      - PROJECT_GANTT_TRACK_GAP
  );
}

function getProjectGanttTrackWidth(totalDaysForTrack, zoom, fitMode, shellWidth = 0) {
  if (fitMode) {
    return calculateProjectGanttFitTrackWidth(shellWidth);
  }
  return Math.max(720, Math.min(24000, totalDaysForTrack * PROJECT_GANTT_BASE_DAY_WIDTH * zoom));
}

function updateProjectGanttViewport(totalDaysForTrack, options = {}) {
  if (!elements.projectManagementContent) return;
  const shell = elements.projectManagementContent.querySelector(".project-gantt-shell-interactive");
  if (!(shell instanceof HTMLElement)) return;
  const zoom = clampProjectGanttZoom(options.zoom ?? state.projectManagementGanttZoom);
  const fitMode = Boolean(options.fitMode);
  const trackWidth = getProjectGanttTrackWidth(totalDaysForTrack, zoom, fitMode, shell.clientWidth || 0);
  const summaryValue = elements.projectManagementContent.querySelector("[data-project-gantt-zoom-value]");
  const summaryNote = elements.projectManagementContent.querySelector("[data-project-gantt-zoom-note]");
  const range = elements.projectManagementContent.querySelector("[data-project-gantt-zoom-range]");
  const fitButton = elements.projectManagementContent.querySelector("[data-project-gantt-zoom-fit]");
  const zoomOutButton = elements.projectManagementContent.querySelector("[data-project-gantt-zoom-out]");
  const zoomInButton = elements.projectManagementContent.querySelector("[data-project-gantt-zoom-in]");

  shell.style.setProperty("--project-gantt-track-width", `${trackWidth}px`);
  shell.dataset.fitMode = fitMode ? "true" : "false";
  if (fitMode && options.resetScroll !== false) {
    shell.scrollLeft = 0;
  }
  if (summaryValue) {
    summaryValue.textContent = fitMode ? "适配全宽" : `${Math.round(zoom * 100)}%`;
  }
  if (summaryNote) {
    summaryNote.textContent = fitMode ? "当前宽度已完整显示整个时间轴。" : "缩小时优先看全跨度，放大后便于查看单个阶段细节。";
  }
  if (range instanceof HTMLInputElement && document.activeElement !== range) {
    range.value = String(Math.max(PROJECT_GANTT_ZOOM_MIN, Math.min(PROJECT_GANTT_ZOOM_MAX, zoom)));
  }
  if (fitButton instanceof HTMLElement) {
    fitButton.classList.toggle("is-active", fitMode);
  }
  if (zoomOutButton instanceof HTMLButtonElement) {
    zoomOutButton.disabled = zoom <= PROJECT_GANTT_ZOOM_MIN;
  }
  if (zoomInButton instanceof HTMLButtonElement) {
    zoomInButton.disabled = zoom >= PROJECT_GANTT_ZOOM_MAX;
  }
}

function renderProjectManagement() {
  if (!elements.projectManagementSummary || !elements.projectManagementContent) return;
  const projects = buildProjectManagementProjects();
  if (!projects.length) {
    elements.projectManagementSummary.innerHTML = "";
    elements.projectManagementContent.innerHTML = `
      <section class="database-alert-rule-block project-management-block">
        <div class="panel-head compact">
          <div>
            <h3>项目进度概览</h3>
            <p class="field-note">当前主数据库尚未形成可用于项目聚合的样本批次。</p>
          </div>
        </div>
      </section>
    `;
    return;
  }
  const archivedProjects = projects.filter((item) => item.archived);
  const editableProjects = projects.filter((item) => !item.archived);
  const activeProjects = editableProjects.filter((item) => item.progress < 90).length;
  const completedProjects = editableProjects.filter((item) => item.progress >= 90).length;
  const overallProgress = Math.round(projects.reduce((sum, item) => sum + item.progress, 0) / projects.length);
  const overallCompletedMilestones = projects.reduce((sum, item) => sum + getProjectOverallScoreMeta(item).completedCount, 0);
  const overallMilestones = projects.reduce((sum, item) => sum + getProjectOverallScoreMeta(item).totalCount, 0);
  const nextDelivery = editableProjects.slice().sort((a, b) => a.timelineEnd.getTime() - b.timelineEnd.getTime())[0] || null;
  const totalSamples = projects.reduce((sum, item) => sum + item.sampleCount, 0);
  const activeSection = state.projectManagementSection === "gantt" ? "gantt" : "overview";
  if (!state.projectManagementSelectedKey || !projects.some((item) => item.key === state.projectManagementSelectedKey)) {
    state.projectManagementSelectedKey = projects[0].key;
  }
  let selectedProject = projects.find((item) => item.key === state.projectManagementSelectedKey) || projects[0];
  const timelineProjects = projects.filter((item) => !item.archived);
  if (activeSection === "gantt" && timelineProjects.length && selectedProject.archived) {
    selectedProject = timelineProjects[0];
    state.projectManagementSelectedKey = selectedProject.key;
    state.projectManagementStageKey = "";
  }
  const resolvedStageKey = String(state.projectManagementStageKey || "").trim();
  const selectedStageSegment = selectedProject.stageSegments.find((segment) => String(segment.key || "") === resolvedStageKey) || null;
  const selectedStageLabel = selectedStageSegment?.label || "";
  const selectedStageMilestone = resolvedStageKey
    ? selectedProject.milestones.find((item) => String(item.id || "") === resolvedStageKey) || null
    : null;
  const monthStarts = buildProjectTimelineMonths(timelineProjects);
  const chartStart = timelineProjects.length
    ? startOfDay(new Date(Math.min(...timelineProjects.map((item) => item.timelineStart.getTime()))))
    : startOfDay(new Date());
  const chartEnd = timelineProjects.length
    ? endOfDay(new Date(Math.max(...timelineProjects.map((item) => item.timelineEnd.getTime()))))
    : endOfDay(new Date());
  const chartWidth = 1120;
  const leftGutter = 248;
  const topGutter = 62;
  const rowHeight = 60;
  const chartHeight = topGutter + timelineProjects.length * rowHeight + 26;
  const plotWidth = chartWidth - leftGutter - 24;
  const totalMs = Math.max(86400000, chartEnd.getTime() - chartStart.getTime());
  const toX = (date) => {
    const value = startOfDay(date) || chartStart;
    const delta = Math.max(0, Math.min(totalMs, value.getTime() - chartStart.getTime()));
    return leftGutter + (delta / totalMs) * plotWidth;
  };
  const toPercent = (date) => {
    const value = startOfDay(date) || chartStart;
    const delta = Math.max(0, Math.min(totalMs, value.getTime() - chartStart.getTime()));
    return (delta / totalMs) * 100;
  };
  const monthBands = monthStarts.map((month) => {
    const next = new Date(month);
    next.setMonth(month.getMonth() + 1);
    const startPercent = toPercent(month);
    const endPercent = Math.min(100, toPercent(next));
    return {
      label: formatProjectMonthLabel(month),
      startPercent,
      widthPercent: Math.max(4, endPercent - startPercent),
    };
  });
  const totalDaysForTrack = Math.max(1, daysBetween(chartStart, chartEnd) + 1);
  const ganttZoom = clampProjectGanttZoom(state.projectManagementGanttZoom);
  if (ganttZoom !== state.projectManagementGanttZoom) {
    state.projectManagementGanttZoom = ganttZoom;
  }
  const sliderZoom = Math.max(PROJECT_GANTT_ZOOM_MIN, Math.min(PROJECT_GANTT_ZOOM_MAX, ganttZoom));
  const ganttTrackWidth = state.projectManagementGanttFitMode
    ? getProjectGanttTrackWidth(totalDaysForTrack, ganttZoom, true, elements.projectManagementContent.querySelector(".project-gantt-shell-interactive")?.clientWidth || 0)
    : getProjectGanttTrackWidth(totalDaysForTrack, ganttZoom, false, 0);
  const ganttZoomPercent = state.projectManagementGanttFitMode ? null : Math.round(ganttZoom * 100);
  const selectedProjectCurrentSegment = selectedProject.stageSegments.find((segment) => Number(segment.completion || 0) < 100) || selectedProject.stageSegments[selectedProject.stageSegments.length - 1];
  const ganttLegendItems = Array.from(
    new Map(
      timelineProjects.flatMap((project) =>
        project.stageSegments.map((segment) => {
          const stage = PROJECT_STAGE_CONFIG.find((item) => item.key === segment.key) || PROJECT_STAGE_CONFIG.find((item) => item.tone === segment.tone) || PROJECT_STAGE_CONFIG[0];
          const segmentLabel = String(segment.label || "").trim() || stage.label;
          return [
            `${segment.key || segmentLabel || stage.tone}`,
            {
              key: segment.key || segmentLabel || stage.tone,
              label: segmentLabel,
              tone: String(segment.tone || "").trim() || stage.tone,
            },
          ];
        })
      )
    ).values()
  );
  const stageMatchesSample = (sample, stageKey) => {
    if (!stageKey) return true;
    const taskStatus = String(sample?.taskStatus || "").trim();
    const releaseLabel = String(sample?.releaseLabel || "").trim();
    if (stageKey === "intake") return true;
    if (stageKey === "analysis") return taskStatus === "分析处理中" || taskStatus === "报告已生成";
    if (stageKey === "review") return taskStatus === "报告已生成";
    if (stageKey === "release") return releaseLabel && releaseLabel !== "未发布";
    if (stageKey === "delivery") return releaseLabel && releaseLabel !== "未发布";
    return true;
  };
  const filteredMilestones = resolvedStageKey
    ? selectedProject.milestones.filter((item) => String(item.id || "") === resolvedStageKey)
    : selectedProject.milestones;
  const filteredSamples = selectedProject.samples.filter((item) => stageMatchesSample(item, resolvedStageKey));
  const filteredProjectLogs = resolvedStageKey
    ? ((resolvedStageKey === "release" || resolvedStageKey === "delivery") ? selectedProject.projectLogs : [])
    : selectedProject.projectLogs;
  const releaseLogCards = Array.from(
    new Map(
      filteredProjectLogs.map((item) => [
        `${item.version_label || item.version_id || item.created_at}`,
        item,
      ])
    ).values()
  ).slice(0, 6);

  elements.projectManagementSummary.innerHTML = `
    <article class="knowledge-base-card project-summary-card project-summary-card-primary">
      <span>在管项目</span>
      <strong>${projects.length} 个</strong>
      <p>覆盖 ${totalSamples} 份样本，按送检单位与月份聚合项目批次。</p>
    </article>
    <article class="knowledge-base-card project-summary-card">
      <span>项目状态</span>
      <strong>${activeProjects} 个进行中</strong>
      <p>已完成 ${completedProjects} 个，已归档 ${archivedProjects.length} 个，当前项目整体节点完成 ${overallCompletedMilestones}/${overallMilestones}。</p>
    </article>
    <article class="knowledge-base-card project-summary-card">
      <span>当前最急交付</span>
      <strong>${escapeHtml(nextDelivery?.shortName || (editableProjects.length ? "-" : "暂无在制项目"))}</strong>
      <p>${escapeHtml(nextDelivery?.monthLabel || "-")} · ${escapeHtml(nextDelivery?.stageLabel || "-")} · ${escapeHtml(nextDelivery ? `预计 ${formatDate(nextDelivery?.timelineEnd) || "-"} 前交付。` : "当前仅剩归档项目。")}</p>
    </article>
    <article class="knowledge-base-card project-summary-card">
      <span>主要专题</span>
      <strong>${escapeHtml(projects[0]?.syndrome || "未标注症候群")}</strong>
      <p>当前样本量最高的项目批次主要围绕 ${escapeHtml(projects[0]?.syndrome || "未标注症候群")} 展开。</p>
    </article>
  `;

  elements.projectManagementContent.innerHTML = `
    <div class="project-management-tabs" role="tablist" aria-label="项目管理视图切换">
      <button class="tab-button database-section-button ${activeSection === "overview" ? "active" : ""}" type="button" data-project-section="overview" aria-selected="${activeSection === "overview" ? "true" : "false"}">项目批次概览</button>
      <button class="tab-button database-section-button ${activeSection === "gantt" ? "active" : ""}" type="button" data-project-section="gantt" aria-selected="${activeSection === "gantt" ? "true" : "false"}">项目甘特图</button>
    </div>
    <div class="project-management-panel ${activeSection === "overview" ? "" : "hidden"}" data-project-panel="overview">
      <div class="project-management-overview-grid">
        <section class="database-alert-rule-block project-management-block project-management-block-rail">
          <div class="panel-head compact">
            <div>
              <h3>项目批次概览</h3>
              <p class="field-note">将同一送检单位在同一月份的样本视作同一项目批次，用于跟踪样本入组、分析、复核和交付进度。</p>
            </div>
          </div>
          <div class="project-management-lanes">
            ${projects.slice(0, 8).map((project) => `
              <article class="project-lane-card tone-${escapeHtml(project.progress >= 90 ? "teal" : project.progress >= 68 ? "amber" : project.progress >= 42 ? "cobalt" : "steel")} ${project.key === selectedProject.key ? "is-active" : ""}" data-project-card="${escapeHtml(project.key)}" tabindex="0" role="button">
                <span>${escapeHtml(project.monthLabel)}</span>
                <strong>${escapeHtml(project.shortName)}</strong>
                <p>${escapeHtml(project.syndrome)} · ${escapeHtml(project.area)}</p>
                <p>${escapeHtml(`${project.sampleCount} 份样本 · ${project.archived ? "已归档" : `当前阶段 ${project.stageLabel}`} · 节点完成 ${getProjectOverallScoreMeta(project).label}`)}</p>
              </article>
            `).join("")}
          </div>
        </section>
        <section class="database-alert-rule-block project-management-block project-management-block-detail">
          <div class="panel-head compact">
            <div>
              <h3>项目协同详情</h3>
              <p class="field-note">围绕当前选中项目查看样本清单、任务关联、主库发布和客户交付节点，便于按项目而不是按单样本推进工作。</p>
            </div>
            <div class="project-detail-head-actions">
              <button class="ghost-button compact" type="button" data-project-archive="${escapeHtml(selectedProject.key)}" data-project-archive-mode="${selectedProject.archived ? "unarchive" : "archive"}">
                ${selectedProject.archived ? "取消归档" : "归档项目"}
              </button>
              <button class="project-delete-button" type="button" data-project-delete="${escapeHtml(selectedProject.key)}">
                <span aria-hidden="true">⌫</span>
                <span>${selectedProject.explicitProject ? "删除项目" : "移除项目"}</span>
              </button>
              ${selectedStageSegment ? `
              <div class="project-stage-filter-banner">
                <span>当前阶段视角</span>
                <strong>${escapeHtml(selectedStageLabel)}</strong>
                <button class="ghost-button compact" type="button" data-project-stage-clear="1">查看全部阶段</button>
              </div>
              ` : ""}
            </div>
          </div>
          <div class="project-detail-grid">
        <article class="project-detail-card">
          <span>当前项目</span>
          <strong>${escapeHtml(selectedProject.name)}</strong>
          <p>${escapeHtml(selectedProject.sampleCount + " 份样本")} · ${escapeHtml(selectedProject.syndrome)} · ${escapeHtml(selectedProject.area)}${selectedProject.archived ? " · 已归档" : ""}</p>
        </article>
        <article class="project-detail-card">
          <span>总负责人</span>
          <strong>${escapeHtml(selectedProject.ownerName)}</strong>
          <p>${escapeHtml(selectedProject.taskNames.slice(0, 2).join(" / ") || "待补充任务")} · ${escapeHtml(`${selectedProject.taskLinkedCount}/${selectedProject.sampleCount} 份样本已关联任务`)}</p>
          <div class="project-detail-card-actions">
            ${selectedProject.archived ? `<span class="field-note">归档项目不可修改负责人</span>` : ""}
            ${selectedProject.archived ? "" : `
            <button class="ghost-button compact" type="button" data-project-owner-edit="${escapeHtml(selectedProject.key)}">修改负责人</button>
            `}
          </div>
        </article>
        <article class="project-detail-card">
          <span>计划交付日</span>
          <strong>${escapeHtml(formatDate(selectedProject.plannedDeliveryDate) || "-")}</strong>
          <p>${escapeHtml(`当前阶段：${selectedProject.stageLabel}`)} · ${escapeHtml(selectedProject.customerStatus)}</p>
        </article>
        <article class="project-detail-card">
          <span>当前阻塞项</span>
          <strong>${escapeHtml(selectedProject.blockerText)}</strong>
          <p>${escapeHtml(`客户状态：${selectedProject.customerStatus}`)}</p>
        </article>
          </div>
          <div class="project-detail-layout">
            <section class="project-detail-section">
              <div class="panel-head compact">
                <div>
                  <h3>交付节点</h3>
                  <p class="field-note">${selectedProject.archived ? "当前项目已归档，仅保留节点只读查看，不再支持编辑，也不会出现在甘特图中。" : "按项目推进链路展示当前阶段，便于客户项目协同管理。"}</p>
                </div>
                <div class="database-alert-topic-actions">
                  ${selectedProject.archived ? "" : `
                  <button class="database-alert-rule-button is-primary" type="button" data-project-milestone-create="${escapeHtml(selectedProject.key)}">新增节点</button>
                  `}
                </div>
              </div>
              <div class="project-milestone-list">
                ${filteredMilestones.map((item) => `
                  <article class="project-milestone-card tone-${item.status}">
                    <div class="project-milestone-head">
                      <div>
                        <strong>${escapeHtml(item.label)}</strong>
                        <p>${escapeHtml(item.detail || "未补充节点说明")}</p>
                      </div>
                      <div class="project-milestone-actions">
                        ${selectedProject.archived ? "" : `
                        <button class="ghost-button compact" type="button" data-project-milestone-edit="${escapeHtml(item.id)}">修改</button>
                        <button class="ghost-button compact danger" type="button" data-project-milestone-delete="${escapeHtml(item.id)}">删除</button>
                        `}
                      </div>
                    </div>
                    <div class="project-milestone-meta">
                      <span>节点负责人 ${escapeHtml(item.owner || selectedProject.ownerName || "待分配")}</span>
                      ${item.startDate && item.endDate ? `<span>时间区段 ${escapeHtml(formatProjectDateRange(item.startDate, item.endDate))}</span>` : ""}
                      <span>节点权重 ${escapeHtml(String(item.weight))}%</span>
                      <span>节点进度 ${escapeHtml(formatStandaloneMilestoneStepScore(item))}</span>
                    </div>
                    ${item.note ? `<p class="project-milestone-note">${escapeHtml(item.note)}</p>` : ""}
                    ${item.children?.length ? `
                      <div class="project-milestone-children-preview">
                        ${item.children.map((child) => `
                          <div class="project-milestone-child-chip">
                            <span>${escapeHtml(child.label)}</span>
                            <strong>${escapeHtml(String(child.completion))}%</strong>
                          </div>
                        `).join("")}
                      </div>
                    ` : ""}
                  </article>
                `).join("")}
              </div>
            </section>
            <section class="project-detail-section">
              <div class="panel-head compact">
                <div>
                  <h3>报告版本</h3>
                  <p class="field-note">${selectedStageSegment ? "当前仅展示与所选阶段相关的版本记录。" : "展示当前项目关联的主库发布与版本记录。"}</p>
                </div>
              </div>
              <div class="project-release-log-list">
                ${releaseLogCards.length ? releaseLogCards.map((item) => `
                  <article class="project-release-log-card">
                    <div>
                      <strong>${escapeHtml(item.version_label || item.version_id || "未命名版本")}</strong>
                      <p>${escapeHtml(item.action_label || item.action || "版本发布")} · ${escapeHtml(formatDateTime(item.created_at) || "-")}</p>
                    </div>
                    <span>${escapeHtml(item.sample_name || item.sample_key || "项目关联样本")}</span>
                  </article>
                `).join("") : `<div class="field-note">当前阶段暂无关联版本记录。</div>`}
              </div>
            </section>
          </div>
        </section>
      </div>
      <section class="database-alert-rule-block project-management-block project-management-block-wide">
        <div class="panel-head compact">
          <div>
            <h3>样本清单</h3>
            <p class="field-note">${selectedStageSegment ? `当前仅展示与“${escapeHtml(selectedStageLabel)}”相关的样本、任务和版本状态。` : "当前项目下的样本、任务和版本状态一览，便于从项目视角统一核查分析与交付准备情况。"}</p>
          </div>
        </div>
        <div class="database-table-shell project-sample-table-shell">
          <div class="database-table-frame">
            <table class="database-table report-table table-tone-assembly project-sample-table">
              <thead>
                <tr>
                  <th>样本</th>
                  <th>采样时间</th>
                  <th>物种</th>
                  <th>任务状态</th>
                  <th>MLST</th>
                  <th>血清型</th>
                  <th>发布版本</th>
                  <th>操作</th>
                </tr>
              </thead>
              <tbody>
                ${filteredSamples.map((item) => `
                  <tr>
                    <td${renderMobileCellAttributes("样本")}>${escapeHtml(item.sampleName)}</td>
                    <td${renderMobileCellAttributes("采样时间")}>${escapeHtml(formatDate(item.date) || "-")}</td>
                    <td${renderMobileCellAttributes("物种")}>${escapeHtml(item.species)}</td>
                    <td${renderMobileCellAttributes("任务状态")}>${escapeHtml(item.taskStatus)}</td>
                    <td${renderMobileCellAttributes("MLST")}>${escapeHtml(item.mlst)}</td>
                    <td${renderMobileCellAttributes("血清型")}>${escapeHtml(item.serotype)}</td>
                    <td${renderMobileCellAttributes("发布版本")}>${escapeHtml(item.releaseLabel)}</td>
                    <td${renderMobileCellAttributes("操作", "actions")}><button class="ghost-button compact" type="button" data-project-sample-key="${escapeHtml(item.sampleKey)}">查看样本</button></td>
                  </tr>
                `).join("") || `<tr><td colspan="8" class="project-sample-table-empty">当前阶段暂无关联样本。</td></tr>`}
              </tbody>
            </table>
          </div>
        </div>
      </section>
    </div>
    <div class="project-management-panel ${activeSection === "gantt" ? "" : "hidden"}" data-project-panel="gantt">
      <section class="database-alert-rule-block project-management-block">
        <div class="panel-head compact">
          <div>
            <h3>项目甘特图</h3>
            <p class="field-note">按项目批次展示样本入组、分析处理中、结果复核和报告交付各阶段，归档项目不会出现在甘特图中。</p>
          </div>
        </div>
        ${timelineProjects.length ? `
        <div class="project-gantt-legend">
          ${ganttLegendItems.map((stage) => `
            <span class="project-gantt-legend-chip tone-${escapeHtml(stage.tone)}">${escapeHtml(stage.label)}</span>
          `).join("")}
        </div>
        <div class="project-gantt-toolbar">
          <div class="project-gantt-toolbar-copy">
            <span>时间缩放</span>
            <strong data-project-gantt-zoom-value>${state.projectManagementGanttFitMode ? "适配全宽" : `${ganttZoomPercent}%`}</strong>
            <p data-project-gantt-zoom-note>${state.projectManagementGanttFitMode ? "当前宽度已完整显示整个时间轴。" : "缩小时优先看全跨度，放大后便于查看单个阶段细节。"}</p>
          </div>
          <div class="project-gantt-toolbar-actions" role="group" aria-label="甘特图缩放控制">
            <button class="project-gantt-toolbar-button" type="button" data-project-gantt-zoom-out="1" ${ganttZoom <= PROJECT_GANTT_ZOOM_MIN ? "disabled" : ""}>缩小</button>
            <div class="project-gantt-zoom-control">
              <input class="project-gantt-zoom-range" type="range" min="${PROJECT_GANTT_ZOOM_MIN}" max="${PROJECT_GANTT_ZOOM_MAX}" step="${PROJECT_GANTT_ZOOM_STEP}" value="${sliderZoom}" data-project-gantt-zoom-range="1" aria-label="甘特图时间缩放">
            </div>
            <button class="project-gantt-toolbar-button" type="button" data-project-gantt-zoom-in="1" ${ganttZoom >= PROJECT_GANTT_ZOOM_MAX ? "disabled" : ""}>放大</button>
            <button class="project-gantt-toolbar-button is-fit ${state.projectManagementGanttFitMode ? "is-active" : ""}" type="button" data-project-gantt-zoom-fit="1">适配全宽</button>
          </div>
        </div>
        <div class="project-gantt-shell project-gantt-shell-interactive" style="--project-gantt-track-width:${ganttTrackWidth}px;">
          <div class="project-gantt-header">
            <div class="project-gantt-header-meta">
              <span>项目</span>
              <strong>当前阶段 / 进度</strong>
            </div>
            <div class="project-gantt-header-scale">
              ${monthBands.map((band, index) => `
                <div class="project-gantt-header-band ${index % 2 === 0 ? "is-even" : "is-odd"}" style="left:${band.startPercent.toFixed(3)}%; width:${band.widthPercent.toFixed(3)}%;">
                  <span>${escapeHtml(band.label)}</span>
                </div>
              `).join("")}
            </div>
          </div>
          <div class="project-gantt-rows">
            ${timelineProjects.map((project) => `
              <article class="project-gantt-row ${project.key === selectedProject.key ? "is-active" : ""}" data-project-gantt-select="${escapeHtml(project.key)}" tabindex="0" role="button">
                <div class="project-gantt-row-meta">
                  <strong>${escapeHtml(project.shortName)}</strong>
                  <p>${escapeHtml(`${project.sampleCount} 份样本 · ${project.syndrome}`)}</p>
                  <span>${escapeHtml(`${project.stageLabel} · ${getProjectOverallScoreMeta(project).label}`)}</span>
                </div>
                <div class="project-gantt-row-track">
                  <div class="project-gantt-row-bands">
                    ${monthBands.map((band, index) => `
                      <div class="project-gantt-row-band ${index % 2 === 0 ? "is-even" : "is-odd"}" style="left:${band.startPercent.toFixed(3)}%; width:${band.widthPercent.toFixed(3)}%;"></div>
                    `).join("")}
                  </div>
                  ${project.stageSegments.map((segment) => {
                    const stage = PROJECT_STAGE_CONFIG.find((item) => item.key === segment.key) || PROJECT_STAGE_CONFIG.find((item) => item.tone === segment.tone) || PROJECT_STAGE_CONFIG[0];
                    const segmentLabel = String(segment.label || "").trim() || stage.label;
                    const startPercent = toPercent(segment.start);
                    const endPercent = Math.min(100, toPercent(addDays(segment.end, 1)));
                    const widthPercent = Math.max(0.12, endPercent - startPercent);
                    const widthPx = ganttTrackWidth * (widthPercent / 100);
                    const progressWidth = Math.max(0, Math.min(100, Number(segment.completion) || 0));
                    const showLabel = widthPx >= 156;
                    const showValue = widthPx >= 252;
                    const compactLabel = !showLabel && widthPx >= 102;
                    const compactValue = !showValue && widthPx >= 96;
                    const compactSegmentLabel = compactLabel
                      ? (segmentLabel.length > 4 ? `${segmentLabel.slice(0, 4)}…` : segmentLabel)
                      : "";
                    const segmentScoreLabel = formatProjectMilestoneStepScore(project, segment);
                    return `
                      <button class="project-gantt-segment tone-${escapeHtml(stage.tone)} ${resolvedStageKey === String(segment.key || "") ? "is-active" : ""}" type="button" data-project-gantt-segment="${escapeHtml(project.key)}" data-project-gantt-stage-key="${escapeHtml(String(segment.key || ""))}" style="left:${startPercent.toFixed(3)}%; width:${widthPercent.toFixed(3)}%;" title="${escapeHtml(`${project.shortName} · ${segmentLabel} · ${formatProjectDateRange(segment.start, segment.end)} · 节点进度 ${segmentScoreLabel}`)}">
                        <span class="project-gantt-segment-progress" style="width:${progressWidth}%;"></span>
                        ${showLabel ? `<span class="project-gantt-segment-label">${escapeHtml(segmentLabel)}</span>` : compactLabel ? `<span class="project-gantt-segment-label is-compact">${escapeHtml(compactSegmentLabel)}</span>` : ""}
                        ${showValue ? `<span class="project-gantt-segment-value">${escapeHtml(segmentScoreLabel)}</span>` : compactValue ? `<span class="project-gantt-segment-value is-compact">${escapeHtml(segmentScoreLabel)}</span>` : ""}
                      </button>
                    `;
                  }).join("")}
                </div>
              </article>
            `).join("")}
          </div>
        </div>
        <div class="project-gantt-focus-card">
          <div class="project-gantt-focus-head">
            <div>
              <p class="section-kicker">Current Project Focus</p>
              <h4>${escapeHtml(selectedProject.name)}</h4>
              <p class="field-note">当前项目的阶段推进、计划交付窗口和节点完成情况会随着甘特图选择联动更新。</p>
            </div>
            <span class="project-gantt-focus-badge">${escapeHtml(getProjectOverallScoreMeta(selectedProject).label)}</span>
          </div>
          <div class="project-gantt-focus-grid">
            <article class="project-gantt-focus-item">
              <span>当前阶段</span>
              <strong>${escapeHtml(selectedProject.stageLabel)}</strong>
              <p>${escapeHtml(selectedProjectCurrentSegment?.label || "待确定")}</p>
            </article>
            <article class="project-gantt-focus-item">
              <span>计划窗口</span>
              <strong>${escapeHtml(formatProjectDateRange(selectedProject.timelineStart, selectedProject.timelineEnd))}</strong>
              <p>${escapeHtml(`计划交付日：${formatDate(selectedProject.plannedDeliveryDate) || "-"}`)}</p>
            </article>
            <article class="project-gantt-focus-item">
              <span>总负责人</span>
              <strong>${escapeHtml(selectedProject.ownerName)}</strong>
              <p>${escapeHtml(selectedProject.customerStatus)}</p>
            </article>
            <article class="project-gantt-focus-item">
              <span>节点负责人</span>
              <strong>${escapeHtml(selectedStageSegment?.owner || selectedProject.ownerName || "待分配")}</strong>
              <p>${escapeHtml(selectedStageSegment?.label || "当前未选中特定节点")}</p>
            </article>
          </div>
          ${selectedStageMilestone ? `
            <div class="project-gantt-stage-detail">
              <div class="project-gantt-stage-detail-head">
                <div>
                  <span>节点备注</span>
                  <strong>${escapeHtml(selectedStageMilestone.label || selectedStageLabel || "当前阶段")}</strong>
                </div>
              </div>
              <label class="project-gantt-stage-editor">
                <span>节点备注</span>
                <textarea data-project-stage-note rows="3" placeholder="补充当前节点的执行要点、阻塞项或交付说明">${escapeHtml(selectedStageMilestone.note || "")}</textarea>
              </label>
              <div class="project-gantt-stage-children">
                <strong>子节点步骤</strong>
                ${selectedStageMilestone.children?.length ? `
                  <div class="project-gantt-stage-child-list">
                    ${selectedStageMilestone.children.map((child, index) => `
                      <article class="project-gantt-stage-child-card" data-project-stage-child-index="${index}">
                        <div>
                          <span>${escapeHtml(child.label)}</span>
                          <textarea data-project-stage-child-note rows="2" placeholder="补充该子节点的执行备注">${escapeHtml(child.note || "")}</textarea>
                        </div>
                        <label class="project-gantt-stage-child-progress">
                          <small>完成度</small>
                          <input data-project-stage-child-completion type="number" min="0" max="100" step="1" value="${escapeHtml(String(child.completion ?? 0))}">
                        </label>
                      </article>
                    `).join("")}
                  </div>
                ` : `<div class="field-note">当前节点尚未拆分子节点步骤。</div>`}
              </div>
              <div class="project-gantt-stage-actions">
                <button class="primary-button" type="button" data-project-stage-save="${escapeHtml(selectedProject.key)}" data-project-stage-save-id="${escapeHtml(selectedStageMilestone.id)}">保存当前节点</button>
              </div>
            </div>
          ` : ""}
          <div class="project-gantt-stage-list">
            ${selectedProject.stageSegments.map((segment) => {
              const stage = PROJECT_STAGE_CONFIG.find((item) => item.key === segment.key) || PROJECT_STAGE_CONFIG.find((item) => item.tone === segment.tone) || PROJECT_STAGE_CONFIG[0];
              const segmentLabel = String(segment.label || "").trim() || stage.label;
              const segmentScoreLabel = formatProjectMilestoneStepScore(selectedProject, segment);
              return `
                <article class="project-gantt-stage-card tone-${escapeHtml(stage.tone)} ${resolvedStageKey === String(segment.key || "") ? "is-active" : ""}">
                  <div>
                    <strong>${escapeHtml(segmentLabel)}</strong>
                    <p>${escapeHtml(formatProjectDateRange(segment.start, segment.end))}</p>
                    <p>${escapeHtml(`负责人：${segment.owner || selectedProject.ownerName || "待分配"}`)}</p>
                  </div>
                  <span>${escapeHtml(segmentScoreLabel)}</span>
                </article>
              `;
            }).join("")}
          </div>
        </div>
        ` : `
        <div class="field-note" style="padding: 12px 0 4px;">
          当前没有可显示的在制项目。所有项目都已归档，或尚未形成可进入甘特图的项目批次。
        </div>
        `}
      </section>
    </div>
  `;

  if (activeSection === "gantt" && timelineProjects.length && state.projectManagementGanttAutoFitPending) {
    const shell = elements.projectManagementContent.querySelector(".project-gantt-shell-interactive");
    if (shell instanceof HTMLElement && shell.clientWidth > 0) {
      state.projectManagementGanttFitMode = true;
      state.projectManagementGanttZoom = calculateProjectGanttFitZoom(totalDaysForTrack, shell.clientWidth);
      state.projectManagementGanttAutoFitPending = false;
      renderProjectManagement();
      return;
    }
  }

  elements.projectManagementContent.querySelectorAll("[data-project-section]").forEach((node) => {
    node.addEventListener("click", () => {
      state.projectManagementSection = String(node.getAttribute("data-project-section") || "overview");
      renderProjectManagement();
    });
  });
  scheduleSegmentedControlsSync(elements.projectManagementContent);
  elements.projectManagementContent.querySelectorAll("[data-project-card]").forEach((node) => {
    const selectProject = () => {
      state.projectManagementSelectedKey = String(node.getAttribute("data-project-card") || "");
      state.projectManagementStageKey = "";
      renderProjectManagement();
    };
    node.addEventListener("click", selectProject);
    node.addEventListener("keydown", (event) => {
      if (event.key === "Enter" || event.key === " ") {
        event.preventDefault();
        selectProject();
      }
    });
  });
  elements.projectManagementContent.querySelectorAll("[data-project-gantt-select]").forEach((node) => {
    const selectProject = () => {
      state.projectManagementSelectedKey = String(node.getAttribute("data-project-gantt-select") || "");
      state.projectManagementStageKey = "";
      renderProjectManagement();
    };
    node.addEventListener("click", selectProject);
    node.addEventListener("keydown", (event) => {
      if (event.key === "Enter" || event.key === " ") {
        event.preventDefault();
        selectProject();
      }
    });
  });
  elements.projectManagementContent.querySelectorAll("[data-project-gantt-segment]").forEach((node) => {
    node.addEventListener("click", (event) => {
      event.stopPropagation();
      state.projectManagementSelectedKey = String(node.getAttribute("data-project-gantt-segment") || "");
      const stageKey = String(node.getAttribute("data-project-gantt-stage-key") || "");
      state.projectManagementStageKey = state.projectManagementStageKey === stageKey ? "" : stageKey;
      renderProjectManagement();
    });
  });
  elements.projectManagementContent.querySelectorAll("[data-project-stage-clear]").forEach((node) => {
    node.addEventListener("click", () => {
      state.projectManagementStageKey = "";
      renderProjectManagement();
    });
  });
  elements.projectManagementContent.querySelectorAll("[data-project-gantt-zoom-out]").forEach((node) => {
    node.addEventListener("click", () => {
      state.projectManagementGanttFitMode = false;
      state.projectManagementGanttZoom = clampProjectGanttZoom(state.projectManagementGanttZoom - PROJECT_GANTT_ZOOM_STEP);
      state.projectManagementGanttAutoFitPending = false;
      renderProjectManagement();
    });
  });
  elements.projectManagementContent.querySelectorAll("[data-project-gantt-zoom-in]").forEach((node) => {
    node.addEventListener("click", () => {
      state.projectManagementGanttFitMode = false;
      state.projectManagementGanttZoom = clampProjectGanttZoom(state.projectManagementGanttZoom + PROJECT_GANTT_ZOOM_STEP);
      state.projectManagementGanttAutoFitPending = false;
      renderProjectManagement();
    });
  });
  elements.projectManagementContent.querySelectorAll("[data-project-gantt-zoom-range]").forEach((node) => {
    let pendingFrame = 0;
    const commitZoom = () => {
      if (pendingFrame) {
        window.cancelAnimationFrame(pendingFrame);
        pendingFrame = 0;
      }
      state.projectManagementGanttFitMode = false;
      state.projectManagementGanttZoom = clampProjectGanttZoom(node.value);
      state.projectManagementGanttAutoFitPending = false;
      renderProjectManagement();
    };
    node.addEventListener("input", () => {
      state.projectManagementGanttFitMode = false;
      state.projectManagementGanttZoom = clampProjectGanttZoom(node.value);
      state.projectManagementGanttAutoFitPending = false;
      if (pendingFrame) {
        window.cancelAnimationFrame(pendingFrame);
      }
      pendingFrame = window.requestAnimationFrame(() => {
        pendingFrame = 0;
        updateProjectGanttViewport(totalDaysForTrack, {
          zoom: state.projectManagementGanttZoom,
          fitMode: false,
          resetScroll: false,
        });
      });
    });
    node.addEventListener("change", commitZoom);
    node.addEventListener("pointerup", commitZoom);
  });
  elements.projectManagementContent.querySelectorAll("[data-project-gantt-zoom-fit]").forEach((node) => {
    node.addEventListener("click", () => {
      const shell = elements.projectManagementContent.querySelector(".project-gantt-shell-interactive");
      state.projectManagementGanttFitMode = true;
      state.projectManagementGanttZoom = calculateProjectGanttFitZoom(totalDaysForTrack, shell?.clientWidth || 0);
      state.projectManagementGanttAutoFitPending = false;
      renderProjectManagement();
    });
  });
  elements.projectManagementContent.querySelectorAll(".project-gantt-shell-interactive").forEach((shell) => {
    if (state.projectManagementGanttFitMode) {
      shell.scrollLeft = 0;
    }
    let dragging = false;
    let startX = 0;
    let startScroll = 0;
    const onPointerMove = (event) => {
      if (!dragging) return;
      const delta = event.clientX - startX;
      shell.scrollLeft = startScroll - delta;
    };
    const endDrag = () => {
      if (!dragging) return;
      dragging = false;
      shell.classList.remove("is-dragging");
      window.removeEventListener("pointermove", onPointerMove);
      window.removeEventListener("pointerup", endDrag);
    };
    shell.addEventListener("pointerdown", (event) => {
      const target = event.target;
      if (!(target instanceof HTMLElement)) return;
      if (target.closest("button, input, textarea, select, a")) return;
      dragging = true;
      startX = event.clientX;
      startScroll = shell.scrollLeft;
      shell.classList.add("is-dragging");
      window.addEventListener("pointermove", onPointerMove);
      window.addEventListener("pointerup", endDrag);
    });
  });
  if (activeSection === "gantt") {
    updateProjectGanttViewport(totalDaysForTrack, {
      zoom: state.projectManagementGanttZoom,
      fitMode: state.projectManagementGanttFitMode,
      resetScroll: state.projectManagementGanttFitMode,
    });
  }
  elements.projectManagementContent.querySelectorAll("[data-project-sample-key]").forEach((node) => {
    node.addEventListener("click", async () => {
      const sampleKey = String(node.getAttribute("data-project-sample-key") || "");
      if (!sampleKey) return;
      await openDatabaseSampleModal(sampleKey);
    });
  });
  elements.projectManagementContent.querySelectorAll("[data-project-milestone-create]").forEach((node) => {
    node.addEventListener("click", () => {
      openProjectMilestoneModal(String(node.getAttribute("data-project-milestone-create") || ""), "");
    });
  });
  elements.projectManagementContent.querySelectorAll("[data-project-owner-edit]").forEach((node) => {
    node.addEventListener("click", () => {
      openProjectOwnerModal(String(node.getAttribute("data-project-owner-edit") || ""));
    });
  });
  elements.projectManagementContent.querySelectorAll("[data-project-archive]").forEach((node) => {
    node.addEventListener("click", () => {
      openProjectDeleteModal(
        String(node.getAttribute("data-project-archive") || ""),
        String(node.getAttribute("data-project-archive-mode") || "archive"),
      );
    });
  });
  elements.projectManagementContent.querySelectorAll("[data-project-delete]").forEach((node) => {
    node.addEventListener("click", () => {
      openProjectDeleteModal(String(node.getAttribute("data-project-delete") || ""));
    });
  });
  elements.projectManagementContent.querySelectorAll("[data-project-milestone-edit]").forEach((node) => {
    node.addEventListener("click", () => {
      openProjectMilestoneModal(selectedProject.key, String(node.getAttribute("data-project-milestone-edit") || ""));
    });
  });
  elements.projectManagementContent.querySelectorAll("[data-project-milestone-delete]").forEach((node) => {
    node.addEventListener("click", () => {
      deleteProjectMilestone(selectedProject.key, String(node.getAttribute("data-project-milestone-delete") || ""));
    });
  });
  elements.projectManagementContent.querySelectorAll("[data-project-stage-save]").forEach((node) => {
    node.addEventListener("click", () => {
      const projectKey = String(node.getAttribute("data-project-stage-save") || "");
      const milestoneId = String(node.getAttribute("data-project-stage-save-id") || "");
      const stageDetail = node.closest(".project-gantt-stage-detail");
      if (!(stageDetail instanceof HTMLElement)) return;
      const note = String(stageDetail.querySelector("[data-project-stage-note]")?.value || "").trim();
      const children = Array.from(stageDetail.querySelectorAll("[data-project-stage-child-index]")).map((row) => {
        const index = Number(row.getAttribute("data-project-stage-child-index") || 0);
        return {
          index,
          note: String(row.querySelector("[data-project-stage-child-note]")?.value || "").trim(),
          completion: Math.max(0, Math.min(100, Number(row.querySelector("[data-project-stage-child-completion]")?.value || 0))),
        };
      });
      const project = getProjectByKey(projectKey);
      const milestone = project?.milestones.find((item) => item.id === milestoneId);
      const nextChildren = Array.isArray(milestone?.children)
        ? milestone.children.map((child, index) => {
            const edited = children.find((item) => item.index === index);
            return edited ? { ...child, note: edited.note, completion: edited.completion } : child;
          })
        : [];
      saveProjectStageExecutionPanel(projectKey, milestoneId, { note, children: nextChildren });
    });
  });
}

function getProjectByKey(projectKey = "") {
  const normalizedKey = String(projectKey || "").trim();
  if (!normalizedKey) return null;
  return buildProjectManagementProjects().find((item) => item.key === normalizedKey) || null;
}

function openProjectMilestoneModal(projectKey = "", milestoneId = "") {
  const project = ensureProjectEditable(projectKey, milestoneId ? "修改交付节点" : "新增交付节点");
  if (!project || !elements.projectMilestoneModal) return;
  clearScopedValidation(elements.projectMilestoneForm);
  const milestone = project.milestones.find((item) => item.id === milestoneId) || null;
  state.currentProjectMilestoneEditor = {
    projectKey: project.key,
    milestoneId: milestone?.id || "",
  };
  if (elements.projectMilestoneTitle) {
    elements.projectMilestoneTitle.textContent = milestone ? "修改交付节点" : "新增交付节点";
  }
  if (elements.projectMilestoneProjectName) {
    elements.projectMilestoneProjectName.textContent = `${project.name} · 当前项目共 ${project.milestones.length} 个节点`;
  }
  if (elements.projectMilestoneNameInput instanceof HTMLInputElement) {
    elements.projectMilestoneNameInput.value = milestone?.label || "";
  }
  if (elements.projectMilestoneOwnerInput instanceof HTMLInputElement) {
    elements.projectMilestoneOwnerInput.value = milestone?.owner || project.ownerName || "";
  }
  if (elements.projectMilestoneStartDateInput instanceof HTMLInputElement) {
    elements.projectMilestoneStartDateInput.value = milestone?.startDate ? String(milestone.startDate).slice(0, 10) : "";
  }
  if (elements.projectMilestoneEndDateInput instanceof HTMLInputElement) {
    elements.projectMilestoneEndDateInput.value = milestone?.endDate ? String(milestone.endDate).slice(0, 10) : "";
  }
  if (elements.projectMilestoneNoteInput instanceof HTMLTextAreaElement) {
    elements.projectMilestoneNoteInput.value = milestone?.note || "";
  }
  if (elements.projectMilestoneWeightInput instanceof HTMLInputElement) {
    elements.projectMilestoneWeightInput.value = milestone ? String(milestone.weight || "") : "";
  }
  if (elements.projectMilestoneCompletionInput instanceof HTMLInputElement) {
    elements.projectMilestoneCompletionInput.value = milestone ? String(milestone.completion || "") : "";
  }
  renderProjectMilestoneChildrenRows(milestone?.children || []);
  if (elements.projectMilestoneDeleteButton instanceof HTMLButtonElement) {
    elements.projectMilestoneDeleteButton.classList.toggle("hidden", !milestone);
  }
  showModalElement(elements.projectMilestoneModal);
}

function closeProjectMilestoneModal() {
  clearScopedValidation(elements.projectMilestoneForm);
  state.currentProjectMilestoneEditor = null;
  hideModalElement(elements.projectMilestoneModal);
}

function openProjectOwnerModal(projectKey = "") {
  const project = ensureProjectEditable(projectKey, "修改项目负责人");
  if (!project || !elements.projectOwnerModal) return;
  clearScopedValidation(elements.projectOwnerForm);
  state.currentProjectOwnerEditor = { projectKey: project.key };
  if (elements.projectOwnerProjectName) {
    elements.projectOwnerProjectName.textContent = `${project.name} · 当前总负责人 ${project.ownerName || "待分配"}`;
  }
  if (elements.projectOwnerNameInput instanceof HTMLInputElement) {
    elements.projectOwnerNameInput.value = project.ownerName || "";
  }
  showModalElement(elements.projectOwnerModal);
}

function closeProjectOwnerModal() {
  clearScopedValidation(elements.projectOwnerForm);
  state.currentProjectOwnerEditor = null;
  hideModalElement(elements.projectOwnerModal);
}

function onSubmitProjectOwner(event) {
  event.preventDefault();
  clearScopedValidation(elements.projectOwnerForm);
  const editor = state.currentProjectOwnerEditor;
  if (!editor?.projectKey) return;
  if (!ensureProjectEditable(editor.projectKey, "修改项目负责人")) return;
  const ownerName = String(elements.projectOwnerNameInput?.value || "").trim();
  if (!ownerName) {
    markFieldInvalid(elements.projectOwnerNameInput, "project-owner-name", "请填写项目总负责人。");
    elements.projectOwnerNameInput?.focus();
    return;
  }
  return withSubmittingState(elements.submitProjectOwnerButton, "保存中...", async () => {
    saveProjectOwnerState(editor.projectKey, ownerName);
    closeProjectOwnerModal();
    renderProjectManagement();
  });
}

function upsertProjectMilestone(projectKey, nextMilestone, milestoneId = "") {
  const project = ensureProjectEditable(projectKey, milestoneId ? "修改交付节点" : "新增交付节点");
  if (!project) return;
  const milestones = project.milestones.slice();
  const normalized = normalizeProjectMilestone({
    ...nextMilestone,
    id: milestoneId || `custom-${Date.now()}`,
  }, milestones.length);
  const index = milestones.findIndex((item) => item.id === milestoneId);
  if (index >= 0) milestones[index] = normalized;
  else milestones.push(normalized);
  saveProjectMilestoneState(project.key, milestones);
  renderProjectManagement();
}

function deleteProjectMilestone(projectKey, milestoneId = "") {
  const project = ensureProjectEditable(projectKey, "删除交付节点");
  if (!project || !milestoneId) return;
  const milestones = project.milestones.filter((item) => item.id !== milestoneId);
  if (!milestones.length) return;
  saveProjectMilestoneState(project.key, milestones);
  if (state.currentProjectMilestoneEditor?.projectKey === project.key && state.currentProjectMilestoneEditor?.milestoneId === milestoneId) {
    closeProjectMilestoneModal();
  }
  renderProjectManagement();
}

function saveProjectStageExecutionPanel(projectKey = "", milestoneId = "", payload = {}) {
  const project = ensureProjectEditable(projectKey, "保存阶段备注");
  if (!project || !milestoneId) return;
  const milestone = project.milestones.find((item) => item.id === milestoneId);
  if (!milestone) return;
  upsertProjectMilestone(project.key, {
    ...milestone,
    note: String(payload.note ?? milestone.note ?? "").trim(),
    children: Array.isArray(payload.children) ? payload.children : milestone.children,
  }, milestoneId);
}

function onDeleteProjectMilestone() {
  const editor = state.currentProjectMilestoneEditor;
  if (!editor?.projectKey || !editor?.milestoneId) return;
  deleteProjectMilestone(editor.projectKey, editor.milestoneId);
}

function onSubmitProjectMilestone(event) {
  event.preventDefault();
  clearScopedValidation(elements.projectMilestoneForm);
  const editor = state.currentProjectMilestoneEditor;
  if (!editor?.projectKey) return;
  if (!ensureProjectEditable(editor.projectKey, editor.milestoneId ? "修改交付节点" : "新增交付节点")) return;
  const label = String(elements.projectMilestoneNameInput?.value || "").trim();
  const owner = String(elements.projectMilestoneOwnerInput?.value || "").trim();
  const startDate = String(elements.projectMilestoneStartDateInput?.value || "").trim();
  const endDate = String(elements.projectMilestoneEndDateInput?.value || "").trim();
  const note = String(elements.projectMilestoneNoteInput?.value || "").trim();
  const weight = Number(elements.projectMilestoneWeightInput?.value || 0);
  const completion = Number(elements.projectMilestoneCompletionInput?.value || 0);
  const children = elements.projectMilestoneChildren instanceof HTMLElement
    ? Array.from(elements.projectMilestoneChildren.querySelectorAll(".project-milestone-child-row"))
        .map((row, index) => {
          const label = String(row.querySelector("[data-project-milestone-child-label]")?.value || "").trim();
          const note = String(row.querySelector("[data-project-milestone-child-note]")?.value || "").trim();
          const completionValue = Number(row.querySelector("[data-project-milestone-child-completion]")?.value || 0);
          if (!label) return null;
          return {
            id: `${editor.milestoneId || "custom"}-child-${index + 1}`,
            label,
            note,
            completion: Number.isFinite(completionValue) ? Math.max(0, Math.min(100, Math.round(completionValue))) : 0,
          };
        })
        .filter(Boolean)
    : [];
  if (!label) {
    markFieldInvalid(elements.projectMilestoneNameInput, "project-milestone-name", "请填写节点名称。");
    elements.projectMilestoneNameInput?.focus();
    return;
  }
  if (!Number.isFinite(weight) || weight <= 0) {
    markFieldInvalid(elements.projectMilestoneWeightInput, "project-milestone-weight", "请填写大于 0 的节点权重。");
    elements.projectMilestoneWeightInput?.focus();
    return;
  }
  if (!Number.isFinite(completion) || completion < 0 || completion > 100) {
    markFieldInvalid(elements.projectMilestoneCompletionInput, "project-milestone-completion", "节点完成百分比需在 0 到 100 之间。");
    elements.projectMilestoneCompletionInput?.focus();
    return;
  }
  if ((startDate && !endDate) || (!startDate && endDate)) {
    markFieldInvalid(elements.projectMilestoneStartDateInput, "project-milestone-start-date", "请同时填写节点开始日期和结束日期，或都留空。");
    markFieldInvalid(elements.projectMilestoneEndDateInput, "project-milestone-end-date", "请同时填写节点开始日期和结束日期，或都留空。");
    return;
  }
  if (startDate && endDate && new Date(endDate).getTime() < new Date(startDate).getTime()) {
    markFieldInvalid(elements.projectMilestoneEndDateInput, "project-milestone-end-date", "节点结束日期不能早于开始日期。");
    elements.projectMilestoneEndDateInput?.focus();
    return;
  }
  return withSubmittingState(elements.submitProjectMilestoneButton, "保存中...", async () => {
    upsertProjectMilestone(editor.projectKey, {
      label,
      owner,
      startDate,
      endDate,
      note,
      weight,
      completion,
      children,
      detail: note || `${label}${startDate && endDate ? ` · ${formatProjectDateRange(startDate, endDate)}` : ""} · 当前完成 ${Math.round(completion)}%，节点权重 ${Math.round(weight)}%。`,
    }, editor.milestoneId || "");
    closeProjectMilestoneModal();
  });
}
