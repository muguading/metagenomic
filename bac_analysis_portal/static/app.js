const state = {
  tasks: [],
  selectedTaskId: "",
  activeTab: "submission-tab",
  pollTimer: null,
  asmOptions: {},
  currentUser: null,
  serverStatus: null,
  queueControls: {
    status: "ALL",
    search: "",
    sort: "created_desc",
  },
  taskDetailOpen: false,
  resultViewerOpen: false,
  rebuildModalOpen: false,
  taskDetailView: "result",
  currentTaskDetail: null,
  currentResultTask: null,
  taskDetailScrollTop: {
    result: 0,
    log: 0,
    meta: 0,
  },
  submissionPanelHome: null,
  batchInput: {
    mode: "single",
    open: false,
    rows: [],
    target: null,
  },
  pathBrowser: {
    open: false,
    selector: "input",
    mode: "project",
    root: "",
    relativePath: "",
    items: [],
    selectedItem: null,
    drag: null,
    currentPath: "",
    withinRoot: true,
  },
};

const elements = {
  form: document.getElementById("task-form"),
  submitButton: document.getElementById("submit-button"),
  demoButton: document.getElementById("demo-button"),
  heroStartButton: document.getElementById("hero-start-button"),
  heroQueueButton: document.getElementById("hero-queue-button"),
  refreshButton: document.getElementById("refresh-button"),
  refreshServerButton: document.getElementById("refresh-server-button"),
  logoutButton: document.getElementById("logout-button"),
  currentUserChip: document.getElementById("current-user-chip"),
  tabButtons: Array.from(document.querySelectorAll("[data-tab-target]")),
  tabPanels: Array.from(document.querySelectorAll(".tab-panel")),
  adminTabButton: document.getElementById("admin-tab-button"),
  chooseInputPathButton: document.getElementById("choose-input-path"),
  openBatchInputButton: document.getElementById("open-batch-input"),
  singleInputModeButton: document.getElementById("single-input-mode"),
  batchInputModeButton: document.getElementById("batch-input-mode"),
  inputPathNote: document.getElementById("input-path-note"),
  chooseOutputDirButton: document.getElementById("choose-output-dir"),
  asmType: document.getElementById("asm_type"),
  asmTypeNote: document.getElementById("asm-type-note"),
  method: document.getElementById("method"),
  longType: document.getElementById("long_type"),
  longTypeWrap: document.getElementById("long-type-wrap"),
  taskList: document.getElementById("task-list"),
  queueSummary: document.getElementById("queue-summary"),
  queueSearch: document.getElementById("queue-search"),
  queueSort: document.getElementById("queue-sort"),
  queueFilterButtons: Array.from(document.querySelectorAll("[data-filter-status]")),
  serverSummary: document.getElementById("server-summary"),
  serverMetrics: document.getElementById("server-metrics"),
  taskDetail: document.getElementById("task-detail"),
  taskStatusChip: document.getElementById("task-status-chip"),
  taskViewResultButton: document.getElementById("task-view-result"),
  taskViewLogButton: document.getElementById("task-view-log"),
  taskViewMetaButton: document.getElementById("task-view-meta"),
  taskDetailModal: document.getElementById("task-detail-modal"),
  resultViewerModal: document.getElementById("result-viewer-modal"),
  resultViewerBackdrop: document.getElementById("result-viewer-backdrop"),
  resultViewerBody: document.getElementById("result-viewer-body"),
  resultViewerStatusChip: document.getElementById("result-viewer-status-chip"),
  resultViewerSubtitle: document.getElementById("result-viewer-subtitle"),
  resultViewerBackButton: document.getElementById("result-viewer-back"),
  closeResultViewerButton: document.getElementById("close-result-viewer"),
  taskDetailBackdrop: document.getElementById("task-detail-backdrop"),
  closeTaskDetailButton: document.getElementById("close-task-detail"),
  rebuildTaskButton: document.getElementById("rebuild-task-button"),
  rebuildModal: document.getElementById("rebuild-modal"),
  rebuildBackdrop: document.getElementById("rebuild-backdrop"),
  rebuildModalPanel: document.getElementById("rebuild-modal-panel"),
  closeRebuildModalButton: document.getElementById("close-rebuild-modal"),
  rebuildPanelSlot: document.getElementById("rebuild-panel-slot"),
  submissionPanel: document.getElementById("submission-panel"),
  submissionTab: document.getElementById("submission-tab"),
  adminPanel: document.getElementById("admin-panel"),
  settingsForm: document.getElementById("settings-form"),
  adminWorkspaceRoot: document.getElementById("admin_workspace_root"),
  adminPipelineScript: document.getElementById("admin_pipeline_script"),
  adminPipelinePython: document.getElementById("admin_pipeline_python"),
  adminMaxConcurrentTasks: document.getElementById("admin_max_concurrent_tasks"),
  chooseWorkspaceRootButton: document.getElementById("choose-workspace-root"),
  choosePipelineScriptButton: document.getElementById("choose-pipeline-script"),
  choosePipelinePythonButton: document.getElementById("choose-pipeline-python"),
  userForm: document.getElementById("user-form"),
  userList: document.getElementById("user-list"),
  newGroupName: document.getElementById("new_group_name"),
  batchInputModal: document.getElementById("batch-input-modal"),
  batchInputBackdrop: document.getElementById("batch-input-backdrop"),
  closeBatchInputButton: document.getElementById("close-batch-input"),
  addBatchRowButton: document.getElementById("add-batch-row"),
  confirmBatchInputButton: document.getElementById("confirm-batch-input"),
  batchInputRows: document.getElementById("batch-input-rows"),
  batchInputPreview: document.getElementById("batch-input-preview"),
  pathBrowserModal: document.getElementById("path-browser-modal"),
  pathBrowserPanel: document.getElementById("path-browser-panel"),
  pathBrowserBackdrop: document.getElementById("path-browser-backdrop"),
  pathBrowserDragHandle: document.getElementById("path-browser-drag-handle"),
  closePathBrowserButton: document.getElementById("close-path-browser"),
  browserGoUpButton: document.getElementById("browser-go-up"),
  browserNewFolderButton: document.getElementById("browser-new-folder"),
  browserRenameButton: document.getElementById("browser-rename"),
  browserSelectCurrentButton: document.getElementById("browser-select-current"),
  browserCurrentPath: document.getElementById("browser-current-path"),
  browserBreadcrumbs: document.getElementById("browser-breadcrumbs"),
  browserSelectedName: document.getElementById("browser-selected-name"),
  browserSelectedPath: document.getElementById("browser-selected-path"),
  browserShortcutProject: document.getElementById("browser-shortcut-project"),
  browserShortcutHome: document.getElementById("browser-shortcut-home"),
  browserShortcutDesktop: document.getElementById("browser-shortcut-desktop"),
  browserShortcutWorkspace: document.getElementById("browser-shortcut-workspace"),
  browserList: document.getElementById("browser-list"),
  browserSubtitle: document.getElementById("path-browser-subtitle"),
  toast: document.getElementById("toast"),
};

document.addEventListener("DOMContentLoaded", async () => {
  initializeBatchInput();
  bindEvents();
  await loadSession();
  applyInitialTabFromUrl();
  await loadAsmOptions();
  await loadTasks();
  await loadServerStatus();
  if (state.currentUser?.role === "admin") {
    await Promise.all([loadAdminSettings(), loadUsers()]);
  }
  startPolling();
});

function applyInitialTabFromUrl() {
  const params = new URLSearchParams(window.location.search);
  const tab = params.get("tab");
  if (tab === "queue") {
    setActiveTab("queue-tab");
    return;
  }
  if (tab === "server") {
    setActiveTab("server-tab");
    return;
  }
  if (tab === "admin" && state.currentUser?.role === "admin") {
    setActiveTab("admin-tab");
  }
}

function bindEvents() {
  elements.form?.addEventListener("submit", onSubmitTask);
  elements.demoButton?.addEventListener("click", onCreateDemoTask);
  elements.tabButtons.forEach((button) => {
    button.addEventListener("click", () => setActiveTab(button.dataset.tabTarget));
  });
  elements.heroStartButton?.addEventListener("click", () => {
    setActiveTab("submission-tab");
    elements.chooseInputPathButton?.focus();
  });
  elements.heroQueueButton?.addEventListener("click", () => setActiveTab("queue-tab"));
  elements.chooseInputPathButton?.addEventListener("click", () => openPathBrowser("input"));
  elements.openBatchInputButton?.addEventListener("click", openBatchInputModal);
  elements.singleInputModeButton?.addEventListener("click", () => setInputMode("single"));
  elements.batchInputModeButton?.addEventListener("click", () => setInputMode("batch"));
  elements.chooseOutputDirButton?.addEventListener("click", () => openPathBrowser("output"));
  elements.asmType?.addEventListener("change", syncAsmTypeState);
  elements.refreshButton?.addEventListener("click", refreshQueueAndModal);
  elements.queueSearch?.addEventListener("input", onQueueSearchChange);
  elements.queueSort?.addEventListener("change", onQueueSortChange);
  elements.queueFilterButtons.forEach((button) => {
    button.addEventListener("click", () => onQueueFilterChange(button.dataset.filterStatus || "ALL"));
  });
  elements.refreshServerButton?.addEventListener("click", loadServerStatus);
  elements.logoutButton?.addEventListener("click", logout);
  elements.settingsForm?.addEventListener("submit", onSaveSettings);
  elements.userForm?.addEventListener("submit", onCreateUser);
  elements.closeBatchInputButton?.addEventListener("click", closeBatchInputModal);
  elements.batchInputBackdrop?.addEventListener("click", closeBatchInputModal);
  elements.addBatchRowButton?.addEventListener("click", addBatchInputRow);
  elements.confirmBatchInputButton?.addEventListener("click", onConfirmBatchInput);
  elements.chooseWorkspaceRootButton?.addEventListener("click", () => openPathBrowser("workspace_root"));
  elements.choosePipelineScriptButton?.addEventListener("click", () => openPathBrowser("script_file"));
  elements.choosePipelinePythonButton?.addEventListener("click", () => openPathBrowser("pipeline_python"));
  elements.closeTaskDetailButton?.addEventListener("click", closeTaskDetail);
  elements.closeResultViewerButton?.addEventListener("click", closeResultViewer);
  elements.resultViewerBackButton?.addEventListener("click", closeResultViewerToQueue);
  elements.resultViewerBackdrop?.addEventListener("click", closeResultViewerToQueue);
  elements.rebuildTaskButton?.addEventListener("click", onRebuildCurrentTask);
  elements.closeRebuildModalButton?.addEventListener("click", closeRebuildModal);
  elements.rebuildBackdrop?.addEventListener("click", closeRebuildModal);
  elements.taskViewResultButton?.addEventListener("click", () => setTaskDetailView("result"));
  elements.taskViewLogButton?.addEventListener("click", () => setTaskDetailView("log"));
  elements.taskViewMetaButton?.addEventListener("click", () => setTaskDetailView("meta"));
  elements.taskDetailBackdrop?.addEventListener("click", closeTaskDetail);
  elements.closePathBrowserButton?.addEventListener("click", closePathBrowser);
  elements.pathBrowserBackdrop?.addEventListener("click", closePathBrowser);
  elements.pathBrowserDragHandle?.addEventListener("pointerdown", startPathBrowserDrag);
  elements.browserGoUpButton?.addEventListener("click", onBrowserGoUp);
  elements.browserNewFolderButton?.addEventListener("click", onBrowserNewFolder);
  elements.browserRenameButton?.addEventListener("click", onBrowserRename);
  elements.browserSelectCurrentButton?.addEventListener("click", onSelectCurrentDirectory);
  elements.browserShortcutProject?.addEventListener("click", () => openBrowserShortcut("project"));
  elements.browserShortcutHome?.addEventListener("click", () => openBrowserShortcut("home"));
  elements.browserShortcutDesktop?.addEventListener("click", () => openBrowserShortcut("desktop"));
  elements.browserShortcutWorkspace?.addEventListener("click", () => openBrowserShortcut("workspace"));
  document.addEventListener("keydown", (event) => {
    if (event.key === "Escape" && state.batchInput.open) {
      closeBatchInputModal();
      return;
    }
    if (event.key === "Escape" && state.rebuildModalOpen) {
      closeRebuildModal();
      return;
    }
    if (event.key === "Escape" && state.resultViewerOpen) {
      closeResultViewerToQueue();
      return;
    }
    if (event.key === "Escape" && state.taskDetailOpen) {
      closeTaskDetail();
      return;
    }
    if (event.key === "Escape" && state.pathBrowser.open) {
      closePathBrowser();
    }
  });
  document.addEventListener("pointermove", onPathBrowserPointerMove);
  document.addEventListener("pointerup", stopPathBrowserDrag);
  window.addEventListener("resize", syncPathBrowserWithinViewport);
}

async function loadSession() {
  const user = await requestJson("/api/session");
  state.currentUser = user;
  elements.currentUserChip.textContent = `${user.username} (${buildRoleLabel(user.role)})`;
  const isAdmin = user.role === "admin";
  elements.adminPanel.classList.toggle("hidden", !isAdmin);
  elements.adminTabButton.classList.toggle("hidden", !isAdmin);
  if (!isAdmin && state.activeTab === "admin-tab") {
    setActiveTab("submission-tab");
  }
}

async function onSubmitTask(event) {
  event.preventDefault();
  const payload = collectFormPayload();
  elements.submitButton.disabled = true;
  elements.submitButton.textContent = "启动中...";
  try {
    const created = await requestJson("/api/tasks", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload),
    });
    showToast(`任务已提交：${created.name}`);
    state.currentTaskDetail = null;
    state.selectedTaskId = created.id;
    state.taskDetailView = "result";
    await loadTasks();
    if (state.rebuildModalOpen) {
      closeRebuildModal();
    }
    setActiveTab("queue-tab");
    await openResultViewerFromTaskId(created.id);
  } finally {
    elements.submitButton.disabled = false;
    elements.submitButton.textContent = "启动任务";
  }
}

async function onCreateDemoTask() {
  elements.demoButton.disabled = true;
  elements.demoButton.textContent = "生成中...";
  try {
    const created = await requestJson("/api/tasks/demo", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
    });
    showToast(`Demo任务已生成：${created.name}`);
    state.currentTaskDetail = null;
    state.selectedTaskId = created.id;
    state.taskDetailView = "result";
    await loadTasks();
    setActiveTab("queue-tab");
    await openResultViewerFromTaskId(created.id);
  } finally {
    elements.demoButton.disabled = false;
    elements.demoButton.textContent = "Demo任务";
  }
}

function collectFormPayload() {
  return {
    task_name: document.getElementById("task_name").value.trim(),
    input_path: document.getElementById("input_path").value.trim(),
    output_dir: document.getElementById("output_dir").value.trim(),
    inputtype: document.getElementById("inputtype").value.trim(),
    thread: document.getElementById("thread").value.trim(),
    minlongfilt: document.getElementById("minlongfilt").value.trim(),
    Qfilt: document.getElementById("Qfilt").value.trim(),
    barcodekit: document.getElementById("barcodekit").value.trim(),
    method: document.getElementById("method").value.trim(),
    long_type: document.getElementById("long_type").value.trim(),
    genome_len: document.getElementById("genome_len").value.trim(),
    asm_type: document.getElementById("asm_type").value.trim(),
    polish_times: document.getElementById("polish_times").value.trim(),
    polish_soft: document.getElementById("polish_soft").value.trim(),
    species: document.getElementById("species").value.trim(),
    runflow: document.getElementById("runflow").value.trim(),
    rmhost: document.getElementById("rmhost").value.trim(),
    abun: document.getElementById("abun").value.trim(),
    rna: document.getElementById("rna").value.trim(),
    fake_pip: document.getElementById("fake_pip").value.trim(),
    ref: document.getElementById("ref").value.trim(),
    gtf: document.getElementById("gtf").value.trim(),
  };
}

function initializeBatchInput() {
  if (!state.batchInput.rows.length) {
    state.batchInput.rows = [createEmptyBatchRow()];
  }
  setInputMode("single");
}

function createEmptyBatchRow() {
  return {
    sample_name: "",
    third_gen: "",
    short_left: "",
    short_right: "",
    species: "",
  };
}

function setInputMode(mode) {
  state.batchInput.mode = mode;
  const isBatch = mode === "batch";
  elements.singleInputModeButton?.classList.toggle("active", !isBatch);
  elements.batchInputModeButton?.classList.toggle("active", isBatch);
  elements.chooseInputPathButton?.classList.toggle("hidden", isBatch);
  elements.openBatchInputButton?.classList.toggle("hidden", !isBatch);
  if (elements.inputPathNote) {
    elements.inputPathNote.textContent = isBatch
      ? "点击“编辑批量表”填写多个样本，确认后自动生成 TSV 并回填到输入路径。"
      : "从当前项目工作区中选择输入文件或输入目录。";
  }
}

function openBatchInputModal() {
  if (!state.batchInput.rows.length) {
    state.batchInput.rows = [createEmptyBatchRow()];
  }
  state.batchInput.open = true;
  elements.batchInputModal?.classList.remove("hidden");
  elements.batchInputModal?.setAttribute("aria-hidden", "false");
  renderBatchInputRows();
}

function closeBatchInputModal() {
  state.batchInput.open = false;
  state.batchInput.target = null;
  elements.batchInputModal?.classList.add("hidden");
  elements.batchInputModal?.setAttribute("aria-hidden", "true");
}

function addBatchInputRow() {
  state.batchInput.rows.push(createEmptyBatchRow());
  renderBatchInputRows();
}

function renderBatchInputRows() {
  if (!elements.batchInputRows) return;
  elements.batchInputRows.replaceChildren();
  state.batchInput.rows.forEach((row, index) => {
    const rowNode = document.createElement("div");
    rowNode.className = "batch-input-row";
    rowNode.innerHTML = `
      <input class="batch-cell-input" data-field="sample_name" data-row-index="${index}" value="${escapeHtml(row.sample_name || "")}" placeholder="例如：sample_01">
      <label class="batch-path-cell">
        <input class="batch-cell-input batch-path-input" value="${escapeHtml(row.third_gen || "")}" readonly placeholder="选择三代数据文件">
        <button class="ghost-button batch-picker-button" type="button" data-field="third_gen" data-row-index="${index}">选择</button>
      </label>
      <label class="batch-path-cell">
        <input class="batch-cell-input batch-path-input" value="${escapeHtml(row.short_left || "")}" readonly placeholder="选择二代左文件">
        <button class="ghost-button batch-picker-button" type="button" data-field="short_left" data-row-index="${index}">选择</button>
      </label>
      <label class="batch-path-cell">
        <input class="batch-cell-input batch-path-input" value="${escapeHtml(row.short_right || "")}" readonly placeholder="选择二代右文件">
        <button class="ghost-button batch-picker-button" type="button" data-field="short_right" data-row-index="${index}">选择</button>
      </label>
      <input class="batch-cell-input" data-field="species" data-row-index="${index}" value="${escapeHtml(row.species || "")}" placeholder="例如：Escherichia coli">
      <button class="ghost-button batch-remove-button" type="button" data-row-index="${index}">删除</button>
    `;
    rowNode.querySelectorAll('.batch-cell-input[data-field]').forEach((input) => {
      input.addEventListener('input', (event) => {
        const target = event.currentTarget;
        updateBatchRowValue(Number(target.dataset.rowIndex), target.dataset.field, target.value);
      });
    });
    rowNode.querySelectorAll('.batch-picker-button').forEach((button) => {
      button.addEventListener('click', () => {
        openBatchFieldPicker(Number(button.dataset.rowIndex), button.dataset.field);
      });
    });
    rowNode.querySelector('.batch-remove-button')?.addEventListener('click', () => {
      removeBatchInputRow(index);
    });
    elements.batchInputRows.appendChild(rowNode);
  });
}

function updateBatchRowValue(index, field, value) {
  if (!state.batchInput.rows[index] || !field) return;
  state.batchInput.rows[index][field] = value;
}

function removeBatchInputRow(index) {
  state.batchInput.rows.splice(index, 1);
  if (!state.batchInput.rows.length) {
    state.batchInput.rows = [createEmptyBatchRow()];
  }
  renderBatchInputRows();
}

function openBatchFieldPicker(index, field) {
  state.batchInput.target = { index, field };
  openPathBrowser("batch_input_file");
}

function setBatchFieldValue(value) {
  const target = state.batchInput.target;
  if (!target || !state.batchInput.rows[target.index]) return;
  state.batchInput.rows[target.index][target.field] = value || "";
  state.batchInput.target = null;
  renderBatchInputRows();
}

async function onConfirmBatchInput() {
  const rows = state.batchInput.rows.map((row) => ({
    sample_name: String(row.sample_name || "").trim(),
    third_gen: String(row.third_gen || "").trim(),
    short_left: String(row.short_left || "").trim(),
    short_right: String(row.short_right || "").trim(),
    species: String(row.species || "").trim(),
  })).filter((row) => Object.values(row).some(Boolean));
  const data = await requestJson('/api/batch-inputs', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ rows }),
  });
  document.getElementById('input_path').value = data.path || data.absolute_path || '';
  if (elements.batchInputPreview) {
    elements.batchInputPreview.value = data.path || data.absolute_path || '';
  }
  setInputMode('batch');
  closeBatchInputModal();
  showToast(`已生成批量 TSV：${data.path}`);
}

async function loadTasks() {
  const data = await requestJson("/api/tasks");
  state.tasks = data.items || [];
  renderTaskList();
}

async function loadServerStatus() {
  const data = await requestJson("/api/server-status");
  state.serverStatus = data;
  renderServerStatus();
}

async function refreshQueueAndModal() {
  await loadTasks();
  if (state.selectedTaskId && state.taskDetailOpen) {
    await loadTaskDetail(state.selectedTaskId);
  }
}

async function loadAdminSettings() {
  const data = await requestJson("/api/admin/settings");
  elements.adminWorkspaceRoot.value = data.workspace_root || "";
  elements.adminPipelineScript.value = data.pipeline_script || "";
  elements.adminPipelinePython.value = data.pipeline_python || "";
  if (elements.adminMaxConcurrentTasks) {
    elements.adminMaxConcurrentTasks.value = String(data.max_concurrent_tasks || 2);
  }
}

async function loadUsers() {
  const data = await requestJson("/api/admin/users");
  renderUserList(data.items || []);
}

async function loadAsmOptions() {
  const data = await requestJson("/api/asm-options");
  state.asmOptions = data.items || {};
  syncAsmTypeState();
}

function syncAsmTypeState() {
  const asmType = elements.asmType.value || "shortasm";
  const options = state.asmOptions[asmType] || [];
  elements.method.innerHTML = options
    .map((item) => `<option value="${escapeHtml(item)}">${escapeHtml(item)}</option>`)
    .join("");
  elements.method.value = options[0] || "";
  elements.asmTypeNote.textContent = buildAsmTypeNote(asmType, options);
  const showLongType = ["longasm", "longref", "shortlongasm"].includes(asmType);
  elements.longTypeWrap.classList.toggle("hidden", !showLongType);
  elements.longType.disabled = !showLongType;
  if (!showLongType) {
    elements.longType.value = "Nanopore";
  }
}

function buildAsmTypeNote(asmType, options) {
  const textMap = {
    shortasm: "short 代表短读长，asm 代表无参组装。",
    longasm: "long 代表长读长，asm 代表无参组装。",
    shortref: "short 代表短读长，ref 代表有参组装。",
    longref: "long 代表长读长，ref 代表有参组装。",
    shortlongasm: "shortlongasm 代表短长读长混合组装，asm 代表无参组装。",
  };
  const suffix = options.length ? `当前可选方法：${options.join("、")}` : "当前没有可用方法。";
  return `${textMap[asmType] || ""} ${suffix}`;
}

function renderTaskList() {
  elements.taskList.replaceChildren();
  const visibleTasks = getVisibleTasks();
  renderQueueSummary(visibleTasks);
  syncQueueControls();
  if (!visibleTasks.length) {
    const empty = document.createElement("div");
    empty.className = "empty-state queue-empty";
    empty.innerHTML = `
      <strong>当前没有符合条件的任务</strong>
      <p>你可以调整状态筛选、搜索关键字或排序方式后再查看。</p>
    `;
    elements.taskList.appendChild(empty);
    return;
  }

  const shell = document.createElement("div");
  shell.className = "task-table-shell";
  const table = document.createElement("table");
  table.className = "task-table";
  table.innerHTML = `
    <thead>
      <tr>
        <th>任务名称</th>
        <th>归属用户</th>
        <th>组装方式</th>
        <th>组装软件</th>
        <th>运行时间</th>
        <th>操作</th>
      </tr>
    </thead>
  `;
  const tbody = document.createElement("tbody");

  visibleTasks.forEach((task) => {
    const row = document.createElement("tr");
    row.className = "task-table-row";
    if (task.id === state.selectedTaskId) {
      row.classList.add("active");
    }
    row.innerHTML = `
      <td>
        <div class="task-name-cell">
          <button class="task-name-button" type="button">${escapeHtml(task.name || task.id)}</button>
          <div class="task-name-meta">
            <span class="mini-chip ${statusClassName(task.status)}">${escapeHtml(task.status || "未知")}</span>
            <span>${escapeHtml(task.id || "-")}</span>
          </div>
        </div>
      </td>
      <td>${escapeHtml(task.owner || "-")}</td>
      <td>${escapeHtml(getAsmTypeLabel(task.params?.asm_type))}</td>
      <td>${escapeHtml(task.params?.method || "-")}</td>
      <td>${escapeHtml(getTaskRunTimeLabel(task))}</td>
      <td>
        <div class="task-table-actions">
          <button class="ghost-button task-action-button" type="button" data-view="log">日志</button>
          <button class="ghost-button task-action-button" type="button" data-view="meta">任务详情</button>
          ${canRebuildTask(task) ? '<button class="ghost-button task-action-button" type="button" data-view="rebuild">重建任务</button>' : ""}
          ${canDeleteTask(task) ? '<button class="ghost-button task-action-button task-danger-button" type="button" data-view="delete">删除任务</button>' : ""}
        </div>
      </td>
    `;

    row.addEventListener("click", (event) => {
      if (event.target.closest(".task-action-button")) {
        return;
      }
      window.location.href = `/tasks/${encodeURIComponent(task.id)}/result-page`;
    });

    row.querySelector(".task-name-button")?.addEventListener("click", (event) => {
      event.stopPropagation();
      window.location.href = `/tasks/${encodeURIComponent(task.id)}/result-page`;
    });

    row.querySelectorAll(".task-action-button").forEach((button) => {
      button.addEventListener("click", async (event) => {
        event.stopPropagation();
        state.currentTaskDetail = null;
        state.selectedTaskId = task.id;
        if ((button.dataset.view || "") === "delete") {
          await onDeleteTask(task);
          return;
        }
        if ((button.dataset.view || "") === "rebuild") {
          await openRebuildFromTaskId(task.id);
          return;
        }
        setTaskDetailView(button.dataset.view || "log");
        renderTaskList();
        await loadTaskDetail(task.id);
        openTaskDetail();
      });
    });

    tbody.appendChild(row);
  });

  table.appendChild(tbody);
  shell.appendChild(table);
  elements.taskList.appendChild(shell);
}

function getVisibleTasks() {
  const search = state.queueControls.search.trim().toLowerCase();
  const filtered = state.tasks.filter((task) => {
    if (state.queueControls.status !== "ALL" && String(task.status || "") !== state.queueControls.status) {
      return false;
    }
    if (!search) {
      return true;
    }
    const haystack = [
      task.name,
      task.id,
      task.owner,
      task.params?.method,
      task.params?.asm_type,
    ].join(" ").toLowerCase();
    return haystack.includes(search);
  });

  return filtered.sort((left, right) => compareTasks(left, right, state.queueControls.sort));
}

function compareTasks(left, right, sortMode) {
  if (sortMode === "created_asc") {
    return compareText(left.created_at, right.created_at);
  }
  if (sortMode === "name_asc") {
    return compareText(left.name || left.id, right.name || right.id);
  }
  if (sortMode === "name_desc") {
    return compareText(right.name || right.id, left.name || left.id);
  }
  if (sortMode === "status") {
    const order = { RUNNING: 0, QUEUED: 1, FAILED: 2, SUCCEEDED: 3 };
    const statusDelta = (order[left.status] ?? 9) - (order[right.status] ?? 9);
    if (statusDelta !== 0) {
      return statusDelta;
    }
    return compareText(right.created_at, left.created_at);
  }
  return compareText(right.created_at, left.created_at);
}

function compareText(left, right) {
  return String(left || "").localeCompare(String(right || ""), "zh-CN");
}

function onQueueFilterChange(status) {
  state.queueControls.status = status || "ALL";
  renderTaskList();
}

function onQueueSearchChange(event) {
  state.queueControls.search = String(event.target.value || "");
  renderTaskList();
}

function onQueueSortChange(event) {
  state.queueControls.sort = String(event.target.value || "created_desc");
  renderTaskList();
}

function syncQueueControls() {
  elements.queueFilterButtons.forEach((button) => {
    button.classList.toggle("active", (button.dataset.filterStatus || "ALL") === state.queueControls.status);
  });
  if (elements.queueSearch && elements.queueSearch.value !== state.queueControls.search) {
    elements.queueSearch.value = state.queueControls.search;
  }
  if (elements.queueSort && elements.queueSort.value !== state.queueControls.sort) {
    elements.queueSort.value = state.queueControls.sort;
  }
}

function getAsmTypeLabel(asmType) {
  const labelMap = {
    shortasm: "短读长无参组装",
    longasm: "长读长无参组装",
    shortref: "短读长有参组装",
    longref: "长读长有参组装",
    shortlongasm: "混合组装",
  };
  return labelMap[String(asmType || "").trim()] || String(asmType || "-");
}

function getTaskRunTimeLabel(task) {
  return formatDate(task.started_at || task.created_at || "");
}

function canDeleteTask(_task) {
  return state.currentUser?.role === "admin" || state.currentUser?.role === "group_admin";
}

function canRebuildTask(task) {
  return state.currentUser?.role === "admin" || state.currentUser?.username === task.owner;
}

async function onDeleteTask(task) {
  const confirmed = window.confirm(`确认删除任务“${task.name || task.id}”吗？该任务目录和日志会一起删除。`);
  if (!confirmed) return;
  await requestJson(`/api/tasks/${encodeURIComponent(task.id)}`, {
    method: "DELETE",
    headers: { "Content-Type": "application/json" },
  });
  if (state.selectedTaskId === task.id) {
    state.selectedTaskId = "";
    state.currentTaskDetail = null;
    if (state.taskDetailOpen) {
      closeTaskDetail();
    }
  }
  showToast(`已删除任务：${task.name || task.id}`);
  await loadTasks();
}

async function loadTaskDetail(taskId) {
  const task = await requestJson(`/api/tasks/${encodeURIComponent(taskId)}`);
  renderTaskDetail(task);
}

async function openResultViewerFromTaskId(taskId) {
  const task = await requestJson(`/api/tasks/${encodeURIComponent(taskId)}`);
  openResultViewer(task);
}

function openResultViewer(task) {
  state.currentResultTask = task;
  state.selectedTaskId = task.id;
  state.resultViewerOpen = true;
  if (elements.resultViewerStatusChip) {
    elements.resultViewerStatusChip.textContent = task.status || "未知";
    elements.resultViewerStatusChip.className = `status-chip ${statusClassName(task.status)}`;
  }
  if (elements.resultViewerSubtitle) {
    const label = task.result_name || task.name || task.id || "任务结果";
    elements.resultViewerSubtitle.textContent = `当前任务：${label}`;
  }
  renderResultViewer(task);
  elements.resultViewerModal?.classList.remove("hidden");
  elements.resultViewerModal?.setAttribute("aria-hidden", "false");
}

function closeResultViewer() {
  state.resultViewerOpen = false;
  elements.resultViewerModal?.classList.add("hidden");
  elements.resultViewerModal?.setAttribute("aria-hidden", "true");
}

function closeResultViewerToQueue() {
  closeResultViewer();
  setActiveTab("queue-tab");
}

function renderResultViewer(task) {
  if (!elements.resultViewerBody) return;
  if (task.result_exists && task.result_url) {
    elements.resultViewerBody.classList.remove("empty-state");
    elements.resultViewerBody.innerHTML = `
      <div class="result-viewer-frame-wrap">
        <iframe class="task-result-frame result-viewer-frame" src="${escapeHtml(task.result_url)}" title="结果预览"></iframe>
      </div>
    `;
    return;
  }
  elements.resultViewerBody.classList.add("empty-state");
  elements.resultViewerBody.innerHTML = `
    <div class="task-result-placeholder result-viewer-placeholder">
      <strong>当前没有可展示的结果页面</strong>
      <p>结果窗口已经打开，后续可以在这里显示真实分析结果；示例任务会优先读取输入目录中的 *_bacgenome.html 文件。</p>
    </div>
  `;
}

function renderTaskDetail(task) {
  const currentScrollNode = elements.taskDetail.querySelector('.task-detail-scroll');
  if (currentScrollNode) {
    state.taskDetailScrollTop[state.taskDetailView] = currentScrollNode.scrollTop;
  }
  state.currentTaskDetail = task;
  state.selectedTaskId = task.id;
  elements.taskStatusChip.textContent = task.status || "未知";
  elements.taskStatusChip.className = `status-chip ${statusClassName(task.status)}`;
  syncTaskDetailViewButtons();
  const command = Array.isArray(task.command) ? task.command.join(" ") : "";
  const params = task.params || {};
  elements.taskDetail.classList.remove("empty-state");

  const logHtml = `
    <div class="task-detail-scroll">
      <div class="detail-block">
        <h3>日志尾部</h3>
        <p class="field-note">显示最近日志输出，适合快速确认当前进展和报错位置。</p>
        <pre>${escapeHtml(task.log_tail || "暂无日志输出")}</pre>
      </div>
    </div>
  `;

  const metaHtml = `
    <div class="task-detail-scroll">
      <div class="detail-grid">
        <div><span>任务 ID</span><strong>${escapeHtml(task.id || "-")}</strong></div>
        <div><span>归属用户</span><strong>${escapeHtml(task.owner || "-")}</strong></div>
        <div><span>创建时间</span><strong>${escapeHtml(formatDate(task.created_at))}</strong></div>
        <div><span>开始时间</span><strong>${escapeHtml(formatDate(task.started_at))}</strong></div>
        <div><span>结束时间</span><strong>${escapeHtml(formatDate(task.finished_at))}</strong></div>
        <div><span>退出码</span><strong>${escapeHtml(task.exit_code ?? "-")}</strong></div>
        <div><span>输出目录</span><strong>${escapeHtml(params.output_dir || "-")}</strong></div>
        <div><span>脚本入口</span><strong>${escapeHtml(task.pipeline_script || "-")}</strong></div>
      </div>
      <div class="detail-block">
        <h3>核心参数</h3>
        <table class="detail-table">
          <tbody>
            ${Object.entries(params).map(([key, value]) => `
              <tr>
                <td>${escapeHtml(key)}</td>
                <td>${escapeHtml(String(value ?? "-"))}</td>
              </tr>
            `).join("")}
          </tbody>
        </table>
      </div>
      <div class="detail-block">
        <h3>执行命令</h3>
        <p class="field-note">用于复现实验或排查参数传递问题。</p>
        <pre>${escapeHtml(command)}</pre>
      </div>
    </div>
  `;

  if (state.taskDetailView === "meta") {
    elements.taskDetail.innerHTML = metaHtml;
  } else {
    elements.taskDetail.innerHTML = logHtml;
  }
  const nextScrollNode = elements.taskDetail.querySelector('.task-detail-scroll');
  if (nextScrollNode) {
    nextScrollNode.scrollTop = state.taskDetailScrollTop[state.taskDetailView] || 0;
    nextScrollNode.addEventListener('scroll', () => {
      state.taskDetailScrollTop[state.taskDetailView] = nextScrollNode.scrollTop;
    }, { passive: true });
  }
  renderTaskList();
}

function openTaskDetail() {
  state.taskDetailOpen = true;
  elements.taskDetailModal.classList.remove("hidden");
  elements.taskDetailModal.setAttribute("aria-hidden", "false");
}

function closeTaskDetail() {
  const currentScrollNode = elements.taskDetail.querySelector('.task-detail-scroll');
  if (currentScrollNode) {
    state.taskDetailScrollTop[state.taskDetailView] = currentScrollNode.scrollTop;
  }
  state.taskDetailOpen = false;
  elements.taskDetailModal.classList.add("hidden");
  elements.taskDetailModal.setAttribute("aria-hidden", "true");
}

function setTaskDetailView(view) {
  const currentScrollNode = elements.taskDetail.querySelector('.task-detail-scroll');
  if (currentScrollNode) {
    state.taskDetailScrollTop[state.taskDetailView] = currentScrollNode.scrollTop;
  }
  state.taskDetailView = view;
  syncTaskDetailViewButtons();
  if (view === "result") {
    if (state.currentTaskDetail?.id) {
      window.location.href = `/tasks/${encodeURIComponent(state.currentTaskDetail.id)}/result-page`;
    }
    return;
  }
  if (state.currentTaskDetail) {
    renderTaskDetail(state.currentTaskDetail);
  }
}

function syncTaskDetailViewButtons() {
  const mapping = [
    [elements.taskViewResultButton, "result"],
    [elements.taskViewLogButton, "log"],
    [elements.taskViewMetaButton, "meta"],
  ];
  mapping.forEach(([button, view]) => {
    if (!button) return;
    button.classList.toggle("active", state.taskDetailView === view);
  });
}

async function onRebuildCurrentTask() {
  if (!state.currentTaskDetail) {
    showToast("请先选择一个任务。", true);
    return;
  }
  await openRebuildModal(state.currentTaskDetail);
}

async function openRebuildFromTaskId(taskId) {
  const task = await requestJson(`/api/tasks/${encodeURIComponent(taskId)}`);
  await openRebuildModal(task);
}

async function openRebuildModal(task) {
  if (!task) {
    showToast("任务信息不存在。", true);
    return;
  }
  populateFormFromTask(task);
  if (state.taskDetailOpen) {
    closeTaskDetail();
  }
  if (state.resultViewerOpen) {
    closeResultViewer();
  }
  ensureSubmissionPanelHome();
  if (elements.submissionPanel && elements.rebuildPanelSlot && elements.submissionPanel.parentElement !== elements.rebuildPanelSlot) {
    elements.rebuildPanelSlot.appendChild(elements.submissionPanel);
  }
  state.rebuildModalOpen = true;
  elements.rebuildModal?.classList.remove("hidden");
  elements.rebuildModal?.setAttribute("aria-hidden", "false");
  elements.submitButton.textContent = "启动任务";
  showToast("已回填任务参数，可继续修改后再提交。");
}

function closeRebuildModal() {
  state.rebuildModalOpen = false;
  restoreSubmissionPanelHome();
  elements.rebuildModal?.classList.add("hidden");
  elements.rebuildModal?.setAttribute("aria-hidden", "true");
}

function ensureSubmissionPanelHome() {
  if (state.submissionPanelHome || !elements.submissionPanel || !elements.submissionTab) {
    return;
  }
  const marker = document.createElement("div");
  marker.id = "submission-panel-home";
  elements.submissionTab.insertBefore(marker, elements.submissionPanel);
  state.submissionPanelHome = marker;
}

function restoreSubmissionPanelHome() {
  if (!state.submissionPanelHome || !elements.submissionPanel) {
    return;
  }
  state.submissionPanelHome.replaceWith(elements.submissionPanel);
  state.submissionPanelHome = null;
}

function populateFormFromTask(task) {
  const params = task.params || {};
  const inputPath = String(params.input_path || "");
  document.getElementById("task_name").value = String(params.task_name || task.name || "");
  document.getElementById("input_path").value = inputPath;
  setInputMode(inputPath.endsWith(".tsv") ? "batch" : "single");
  document.getElementById("output_dir").value = String(params.output_dir || "");
  document.getElementById("inputtype").value = String(params.inputtype || "fastq");
  document.getElementById("thread").value = String(params.thread ?? 10);
  document.getElementById("minlongfilt").value = String(params.minlongfilt ?? "500");
  document.getElementById("Qfilt").value = String(params.Qfilt ?? "10");
  document.getElementById("barcodekit").value = String(params.barcodekit ?? "none");
  document.getElementById("asm_type").value = String(params.asm_type || "shortasm");
  syncAsmTypeState();
  document.getElementById("method").value = String(params.method || document.getElementById("method").value || "");
  document.getElementById("long_type").value = String(params.long_type || "Nanopore");
  document.getElementById("genome_len").value = String(params.genome_len ?? "4m");
  document.getElementById("polish_times").value = String(params.polish_times ?? "1");
  document.getElementById("polish_soft").value = String(params.polish_soft ?? "medaka");
  document.getElementById("species").value = String(params.species ?? "False");
  document.getElementById("runflow").value = String(params.runflow ?? "All");
  document.getElementById("rmhost").value = String(params.rmhost ?? "norm");
  document.getElementById("abun").value = String(params.abun ?? "1");
  document.getElementById("rna").value = String(params.rna ?? "0");
  document.getElementById("fake_pip").value = String(params.fake_pip ?? 0);
  document.getElementById("ref").value = params.ref && params.ref !== "noref" ? String(params.ref) : "";
  document.getElementById("gtf").value = params.gtf && params.gtf !== "nogtf" ? String(params.gtf) : "";
}

function renderUserList(items) {
  elements.userList.replaceChildren();
  if (!items.length) {
    elements.userList.innerHTML = '<div class="empty-state queue-empty"><strong>暂无用户</strong></div>';
    return;
  }
  items.forEach((user) => {
    const row = document.createElement("article");
    row.className = "user-card";
    row.innerHTML = `
      <div class="user-card-head">
        <strong>${escapeHtml(user.username)}</strong>
        <span class="mini-chip ${user.role === "admin" ? "running" : "queued"}">${escapeHtml(user.role)}</span>
      </div>
      <p>${escapeHtml(user.display_name || "未设置显示名")}</p>
      <p>${escapeHtml(user.group_name || "未分组")}</p>
      <p>${escapeHtml(formatDate(user.created_at))}</p>
      <div class="user-card-actions">
        <button class="ghost-button user-edit-button" type="button">修改</button>
        <button class="ghost-button task-danger-button user-delete-button" type="button">删除</button>
      </div>
    `;
    row.querySelector(".user-edit-button")?.addEventListener("click", async () => {
      const nextUsername = window.prompt("输入用户名", user.username || "");
      if (nextUsername === null) return;
      const nextGroupName = window.prompt("输入用户组，管理员可留空", user.group_name || "");
      if (nextGroupName === null) return;
      const displayName = window.prompt("输入显示名", user.display_name || "");
      if (displayName === null) return;
      const role = window.prompt("输入角色：admin / group_admin / user", user.role);
      if (role === null) return;
      const password = window.prompt("如需重置密码请输入新密码，留空则不修改", "");
      await requestJson(`/api/admin/users/${encodeURIComponent(user.username)}`, {
        method: "PUT",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          username: nextUsername,
          group_name: nextGroupName,
          display_name: displayName,
          role: role,
          password: password || "",
        }),
      });
      showToast(`已更新用户：${nextUsername || user.username}`);
      await loadUsers();
    });
    row.querySelector(".user-delete-button")?.addEventListener("click", async () => {
      const confirmed = window.confirm(`确认删除用户“${user.username}”吗？`);
      if (!confirmed) return;
      await requestJson(`/api/admin/users/${encodeURIComponent(user.username)}`, {
        method: "DELETE",
        headers: { "Content-Type": "application/json" },
      });
      showToast(`已删除用户：${user.username}`);
      await loadUsers();
    });
    elements.userList.appendChild(row);
  });
}

function setActiveTab(tabId) {
  state.activeTab = tabId;
  elements.tabButtons.forEach((button) => {
    button.classList.toggle("active", button.dataset.tabTarget === tabId);
  });
  elements.tabPanels.forEach((panel) => {
    panel.classList.toggle("active", panel.id === tabId);
  });
  if (tabId === "server-tab") {
    loadServerStatus().catch((error) => console.error(error));
  }
}

async function openPathBrowser(selector) {
  if (selector === "script_file" && !elements.adminWorkspaceRoot.value.trim()) {
    showToast("请先选择部署基准目录。", true);
    return;
  }
  state.pathBrowser.open = true;
  state.pathBrowser.selector = selector;
  state.pathBrowser.mode = ["workspace_root", "script_file", "pipeline_python"].includes(selector) ? "admin" : "project";
  state.pathBrowser.root = selector === "script_file"
    ? (elements.adminWorkspaceRoot.value.trim() || "")
    : "";
  state.pathBrowser.relativePath = "";
  state.pathBrowser.selectedItem = null;
  resetPathBrowserPosition();
  elements.pathBrowserModal.classList.remove("hidden");
  elements.pathBrowserModal.setAttribute("aria-hidden", "false");
  elements.browserSubtitle.textContent = describeBrowserSelector(selector);
  const adminBrowse = state.pathBrowser.mode === "admin";
  elements.browserNewFolderButton.disabled = adminBrowse;
  elements.browserRenameButton.disabled = adminBrowse;
  elements.browserSelectCurrentButton.disabled = false;
  elements.browserShortcutWorkspace.classList.toggle("hidden", !elements.adminWorkspaceRoot.value.trim());
  syncBrowserShortcuts();
  await loadBrowserDirectory(state.pathBrowser.mode === "admin" ? state.pathBrowser.root : "");
}

function closePathBrowser() {
  state.pathBrowser.open = false;
  state.pathBrowser.drag = null;
  elements.pathBrowserPanel?.classList.remove("is-dragging");
  elements.pathBrowserModal.classList.add("hidden");
  elements.pathBrowserModal.setAttribute("aria-hidden", "true");
}

async function loadBrowserDirectory(relativePath) {
  const params = new URLSearchParams({ selector: state.pathBrowser.selector });
  if (state.pathBrowser.mode === "admin") {
    if (state.pathBrowser.root) {
      params.set("root", state.pathBrowser.root);
    }
    if (relativePath) {
      params.set("path", relativePath);
    }
  } else {
    params.set("path", relativePath || "");
  }
  const endpoint = state.pathBrowser.mode === "admin" ? "/api/admin/filesystem" : "/api/filesystem";
  const data = await requestJson(`${endpoint}?${params.toString()}`);
  state.pathBrowser.root = data.root || state.pathBrowser.root;
  state.pathBrowser.currentPath = data.current_path || "";
  state.pathBrowser.withinRoot = data.within_root !== false;
  state.pathBrowser.relativePath = state.pathBrowser.mode === "admin"
    ? (data.current_path || "")
    : (data.relative_path || "");
  state.pathBrowser.items = data.items || [];
  state.pathBrowser.selectedItem = null;
  elements.browserCurrentPath.textContent = data.current_path || "/";
  const parentPath = state.pathBrowser.mode === "admin" ? (data.parent_path || "") : (data.parent_relative_path || "");
  elements.browserGoUpButton.disabled = !parentPath || data.current_path === '/';
  elements.browserGoUpButton.dataset.targetPath = parentPath || "";
  elements.browserNewFolderButton.disabled = state.pathBrowser.mode === "admin" || !state.pathBrowser.withinRoot;
  elements.browserRenameButton.disabled = !state.pathBrowser.selectedItem || state.pathBrowser.mode === "admin" || !state.pathBrowser.withinRoot;
  syncBrowserSelection();
  renderBrowserBreadcrumbs();
  syncBrowserShortcuts();
  renderBrowserList(data.items || []);
}

function renderBrowserList(items) {
  elements.browserList.replaceChildren();
  if (!items.length) {
    const empty = document.createElement("div");
    empty.className = "empty-state queue-empty";
    empty.innerHTML = "<strong>当前目录为空</strong><p>请返回上级目录或切换其他位置。</p>";
    elements.browserList.appendChild(empty);
    return;
  }

  items.forEach((item) => {
    const button = document.createElement("button");
    button.type = "button";
    button.className = `browser-item ${item.type}`;
    button.innerHTML = `
      <strong>${escapeHtml(item.name)}</strong>
      <span>${item.type === "directory" ? "目录" : "文件"}</span>
      <code>${escapeHtml(item.path)}</code>
    `;
    button.addEventListener("click", () => {
      state.pathBrowser.selectedItem = item;
      syncBrowserSelection();
      renderBrowserList(state.pathBrowser.items);
    });
    button.addEventListener("dblclick", async () => {
      await onOpenBrowserItem(item);
    });
    button.classList.toggle("active", state.pathBrowser.selectedItem?.path === item.path);
    elements.browserList.appendChild(button);
  });
}

async function onOpenBrowserItem(item) {
  if (item.type === "directory") {
    await loadBrowserDirectory(item.path);
    return;
  }
  if (state.pathBrowser.selector === "input") {
    setSelectedPath("input_path", item.path);
    closePathBrowser();
    return;
  }
  if (state.pathBrowser.selector === "batch_input_file") {
    setBatchFieldValue(item.path);
    closePathBrowser();
    return;
  }
  if (state.pathBrowser.selector === "script_file") {
    setAdminScriptPath(item.path);
    closePathBrowser();
    return;
  }
  if (state.pathBrowser.selector === "pipeline_python") {
    setAdminPipelinePython(item.path);
    closePathBrowser();
  }
}

async function onBrowserGoUp() {
  await loadBrowserDirectory(elements.browserGoUpButton.dataset.targetPath || "");
}

function startPathBrowserDrag(event) {
  if (!state.pathBrowser.open) return;
  const interactive = event.target.closest("button, input, textarea, select, a");
  if (interactive) return;
  const panel = elements.pathBrowserPanel;
  if (!panel) return;
  const rect = panel.getBoundingClientRect();
  state.pathBrowser.drag = {
    offsetX: event.clientX - rect.left,
    offsetY: event.clientY - rect.top,
  };
  panel.classList.add("is-positioned");
  panel.classList.add("is-dragging");
}

function onPathBrowserPointerMove(event) {
  if (!state.pathBrowser.drag || !state.pathBrowser.open) return;
  const panel = elements.pathBrowserPanel;
  if (!panel) return;
  const margin = 12;
  const width = rectWidth(panel);
  const height = rectHeight(panel);
  const maxLeft = Math.max(margin, window.innerWidth - width - margin);
  const maxTop = Math.max(margin, window.innerHeight - height - margin);
  const left = clamp(event.clientX - state.pathBrowser.drag.offsetX, margin, maxLeft);
  const top = clamp(event.clientY - state.pathBrowser.drag.offsetY, margin, maxTop);
  panel.style.left = `${left}px`;
  panel.style.top = `${top}px`;
  panel.style.transform = "none";
}

function stopPathBrowserDrag() {
  if (!state.pathBrowser.drag) return;
  state.pathBrowser.drag = null;
  elements.pathBrowserPanel?.classList.remove("is-dragging");
}

function resetPathBrowserPosition() {
  const panel = elements.pathBrowserPanel;
  if (!panel) return;
  panel.classList.remove("is-positioned");
  panel.style.left = "";
  panel.style.top = "";
  panel.style.transform = "";
}

function syncPathBrowserWithinViewport() {
  if (!state.pathBrowser.open) return;
  const panel = elements.pathBrowserPanel;
  if (!panel) return;
  if (!panel.style.left && !panel.style.top) return;
  const margin = 12;
  const width = rectWidth(panel);
  const height = rectHeight(panel);
  const left = clamp(parseFloat(panel.style.left || "0"), margin, Math.max(margin, window.innerWidth - width - margin));
  const top = clamp(parseFloat(panel.style.top || "0"), margin, Math.max(margin, window.innerHeight - height - margin));
  panel.style.left = `${left}px`;
  panel.style.top = `${top}px`;
}

function onSelectCurrentDirectory() {
  if (state.pathBrowser.selector === "workspace_root") {
    const path = state.pathBrowser.selectedItem?.type === "directory"
      ? state.pathBrowser.selectedItem.path
      : state.pathBrowser.relativePath;
    elements.adminWorkspaceRoot.value = path || "";
    elements.adminPipelineScript.value = "";
    closePathBrowser();
    return;
  }
  if (state.pathBrowser.selector === "script_file" && state.pathBrowser.selectedItem?.type === "file") {
    setAdminScriptPath(state.pathBrowser.selectedItem.path);
    closePathBrowser();
    return;
  }
  if (state.pathBrowser.selector === "pipeline_python" && state.pathBrowser.selectedItem?.type === "file") {
    setAdminPipelinePython(state.pathBrowser.selectedItem.path);
    closePathBrowser();
    return;
  }
  if (state.pathBrowser.selector === "input" && state.pathBrowser.selectedItem?.type === "file") {
    setSelectedPath("input_path", state.pathBrowser.selectedItem.path);
    closePathBrowser();
    return;
  }
  if (state.pathBrowser.selector === "batch_input_file") {
    if (state.pathBrowser.selectedItem?.type === "file") {
      setBatchFieldValue(state.pathBrowser.selectedItem.path);
      closePathBrowser();
      return;
    }
    showToast("批量输入列需要选择文件。", true);
    return;
  }
  if (state.pathBrowser.selectedItem?.type === "directory") {
    const fieldId = state.pathBrowser.selector === "output" ? "output_dir" : "input_path";
    setSelectedPath(fieldId, state.pathBrowser.selectedItem.path);
    closePathBrowser();
    return;
  }
  const fieldId = state.pathBrowser.selector === "output" ? "output_dir" : "input_path";
  setSelectedPath(fieldId, state.pathBrowser.relativePath || ".");
  closePathBrowser();
}

function setSelectedPath(fieldId, value) {
  document.getElementById(fieldId).value = value || ".";
}

function syncBrowserSelection() {
  const selected = state.pathBrowser.selectedItem;
  elements.browserSelectedName.textContent = selected
    ? `${selected.name} (${selected.type === "directory" ? "目录" : "文件"})`
    : "未选择";
  elements.browserSelectedPath.value = selected
    ? selected.path
    : (state.pathBrowser.relativePath || "");
  elements.browserRenameButton.disabled = !selected || state.pathBrowser.mode === "admin" || !state.pathBrowser.withinRoot;
}

async function onBrowserNewFolder() {
  if (state.pathBrowser.mode === "admin") {
    showToast("部署目录浏览器暂不支持在系统路径中新建文件夹。", true);
    return;
  }
  const name = window.prompt("输入新文件夹名称");
  if (!name) return;
  await requestJson("/api/filesystem/mkdir", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({
      path: state.pathBrowser.relativePath || "",
      name: name.trim(),
    }),
  });
  showToast(`已创建文件夹：${name.trim()}`);
  await loadBrowserDirectory(state.pathBrowser.relativePath || "");
}

async function onBrowserRename() {
  if (state.pathBrowser.mode === "admin") {
    showToast("部署目录浏览器暂不支持在系统路径中重命名。", true);
    return;
  }
  const selected = state.pathBrowser.selectedItem;
  if (!selected) {
    showToast("请先选择一个文件或目录。", true);
    return;
  }
  const name = window.prompt("输入新名称", selected.name);
  if (!name || name.trim() === selected.name) return;
  await requestJson("/api/filesystem/rename", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({
      path: selected.path,
      name: name.trim(),
    }),
  });
  showToast(`已重命名为：${name.trim()}`);
  await loadBrowserDirectory(state.pathBrowser.relativePath || "");
}

async function logout() {
  await requestJson("/logout", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
  });
  window.location.href = "/login";
}

async function onSaveSettings(event) {
  event.preventDefault();
  const data = await requestJson("/api/admin/settings", {
    method: "PUT",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({
      workspace_root: elements.adminWorkspaceRoot.value.trim(),
      pipeline_script: elements.adminPipelineScript.value.trim(),
      pipeline_python: elements.adminPipelinePython.value.trim(),
      max_concurrent_tasks: elements.adminMaxConcurrentTasks?.value?.trim() || "2",
    }),
  });
  elements.adminWorkspaceRoot.value = data.workspace_root || "";
  elements.adminPipelineScript.value = data.pipeline_script || "";
  elements.adminPipelinePython.value = data.pipeline_python || "";
  if (elements.adminMaxConcurrentTasks) {
    elements.adminMaxConcurrentTasks.value = String(data.max_concurrent_tasks || 2);
  }
  showToast("脚本路径已更新");
}

async function onCreateUser(event) {
  event.preventDefault();
  await requestJson("/api/admin/users", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({
      username: document.getElementById("new_username").value.trim(),
      display_name: document.getElementById("new_display_name").value.trim(),
      group_name: document.getElementById("new_group_name").value.trim(),
      password: document.getElementById("new_password").value,
      role: document.getElementById("new_role").value,
    }),
  });
  elements.userForm.reset();
  showToast("用户已创建");
  await loadUsers();
}

function startPolling() {
  window.clearInterval(state.pollTimer);
  state.pollTimer = window.setInterval(async () => {
    try {
      await loadTasks();
      if (state.activeTab === "server-tab") {
        await loadServerStatus();
      }
      if (state.selectedTaskId && state.taskDetailOpen) {
        await loadTaskDetail(state.selectedTaskId);
      }
      if (state.currentUser?.role === "admin" && state.activeTab === "admin-tab") {
        await loadAdminSettings();
      }
    } catch (error) {
      console.error(error);
    }
  }, 5000);
}

async function requestJson(url, options = {}) {
  const response = await fetch(url, options);
  const data = await response.json();
  if (!response.ok) {
    showToast(data.error || "请求失败", true);
    throw new Error(data.error || "请求失败");
  }
  return data;
}

function formatDate(value) {
  if (!value) return "-";
  const date = new Date(value);
  return Number.isNaN(date.getTime()) ? String(value) : date.toLocaleString();
}

function renderQueueSummary(items = state.tasks) {
  const total = items.length;
  const running = items.filter((task) => task.status === "RUNNING").length;
  const queued = items.filter((task) => task.status === "QUEUED").length;
  const failed = items.filter((task) => task.status === "FAILED").length;
  const succeeded = items.filter((task) => task.status === "SUCCEEDED").length;
  elements.queueSummary.innerHTML = `
    <article><span>总任务</span><strong>${total}</strong></article>
    <article><span>运行中</span><strong>${running}</strong></article>
    <article><span>已完成</span><strong>${succeeded}</strong></article>
    <article><span>失败</span><strong>${failed}</strong></article>
    <article><span>排队中</span><strong>${queued}</strong></article>
  `;
}

function renderServerStatus() {
  if (!elements.serverSummary || !elements.serverMetrics) return;
  const status = state.serverStatus;
  if (!status) {
    elements.serverSummary.innerHTML = "";
    elements.serverMetrics.innerHTML = `
      <div class="empty-state queue-empty">
        <strong>暂无服务器信息</strong>
        <p>请点击刷新状态重新获取当前服务器资源情况。</p>
      </div>
    `;
    return;
  }

  elements.serverSummary.innerHTML = `
    <article><span>主机名</span><strong>${escapeHtml(status.hostname || "-")}</strong></article>
    <article><span>操作系统</span><strong>${escapeHtml(status.platform || "-")}</strong></article>
    <article><span>CPU 核心</span><strong>${escapeHtml(String(status.cpu_count || "-"))}</strong></article>
    <article><span>总内存</span><strong>${escapeHtml(status.memory?.total_human || "-")}</strong></article>
    <article><span>磁盘剩余</span><strong>${escapeHtml(status.disk?.free_human || "-")}</strong></article>
  `;

  elements.serverMetrics.innerHTML = `
    <article class="server-card">
      <div class="server-card-head">
        <div>
          <p class="section-kicker">Hardware</p>
          <h3>硬件概况</h3>
        </div>
        <span class="server-update-time">更新于 ${escapeHtml(formatDate(status.sampled_at))}</span>
      </div>
      <div class="server-detail-grid">
        <div><span>主机名</span><strong>${escapeHtml(status.hostname || "-")}</strong></div>
        <div><span>系统平台</span><strong>${escapeHtml(status.platform || "-")}</strong></div>
        <div><span>架构</span><strong>${escapeHtml(status.machine || "-")}</strong></div>
        <div><span>Python</span><strong>${escapeHtml(status.python || "-")}</strong></div>
        <div><span>CPU 逻辑核心</span><strong>${escapeHtml(String(status.cpu_count || "-"))}</strong></div>
        <div><span>部署磁盘路径</span><strong>${escapeHtml(status.disk?.path || "-")}</strong></div>
      </div>
    </article>
    <article class="server-card">
      <div class="server-card-head">
        <div>
          <p class="section-kicker">Usage</p>
          <h3>资源占用</h3>
        </div>
      </div>
      <div class="server-usage-grid">
        ${renderUsageMetric("CPU 使用率", status.cpu?.percent, status.cpu?.detail)}
        ${renderUsageMetric("内存使用率", status.memory?.percent, `${status.memory?.used_human || "-"} / ${status.memory?.total_human || "-"}`)}
        ${renderUsageMetric("存储使用率", status.disk?.used_percent, `${status.disk?.used_human || "-"} / ${status.disk?.total_human || "-"}`)}
      </div>
      <div class="server-detail-grid">
        <div><span>内存剩余</span><strong>${escapeHtml(status.memory?.free_human || "-")}</strong></div>
        <div><span>磁盘剩余</span><strong>${escapeHtml(status.disk?.free_human || "-")}</strong></div>
        <div><span>1 分钟负载</span><strong>${escapeHtml(status.cpu?.load_average || "-")}</strong></div>
        <div><span>采样方式</span><strong>${escapeHtml(status.cpu?.source || "-")}</strong></div>
      </div>
    </article>
  `;
}

function renderUsageMetric(label, percent, detail) {
  const normalized = Number.isFinite(Number(percent)) ? Math.max(0, Math.min(100, Number(percent))) : 0;
  return `
    <div class="usage-metric-card">
      <span>${escapeHtml(label)}</span>
      <strong>${escapeHtml(formatPercent(percent))}</strong>
      <div class="usage-meter" aria-hidden="true">
        <div class="usage-meter-fill" style="width:${normalized}%"></div>
      </div>
      <small>${escapeHtml(detail || "-")}</small>
    </div>
  `;
}

function formatPercent(value) {
  const numeric = Number(value);
  if (!Number.isFinite(numeric)) return "-";
  return `${numeric.toFixed(1)}%`;
}

function buildRoleLabel(role) {
  if (role === "admin") return "管理员";
  if (role === "group_admin") return "group 管理";
  return "普通用户";
}

function escapeHtml(value) {
  return String(value)
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;")
    .replaceAll("'", "&#39;");
}

function statusClassName(status) {
  return String(status || "").toLowerCase();
}

function showToast(message, isError = false) {
  elements.toast.textContent = message;
  elements.toast.style.background = isError ? "rgba(143, 51, 47, 0.96)" : "rgba(39, 68, 109, 0.94)";
  elements.toast.classList.remove("hidden");
  window.clearTimeout(showToast.timer);
  showToast.timer = window.setTimeout(() => elements.toast.classList.add("hidden"), 2600);
}

function clamp(value, min, max) {
  return Math.min(Math.max(value, min), max);
}

function rectWidth(node) {
  return node.getBoundingClientRect().width || node.offsetWidth || 0;
}

function rectHeight(node) {
  return node.getBoundingClientRect().height || node.offsetHeight || 0;
}

function describeBrowserSelector(selector) {
  if (selector === "input") return "选择输入文件或输入目录。";
  if (selector === "batch_input_file") return "为批量输入行选择文件。";
  if (selector === "output") return "选择输出目录。";
  if (selector === "workspace_root") return "选择部署基准目录。";
  if (selector === "script_file") return "从部署基准目录中选择脚本文件。";
  if (selector === "pipeline_python") return "选择用于运行 Bac_assemble_260112_newformat.py 的 Python 可执行文件。";
  return "浏览路径。";
}

function setAdminPipelinePython(absolutePath) {
  if (elements.adminPipelinePython) {
    elements.adminPipelinePython.value = absolutePath || "";
  }
}

function setAdminScriptPath(absolutePath) {
  const workspaceRoot = elements.adminWorkspaceRoot.value.trim();
  if (!workspaceRoot) {
    showToast("请先选择部署基准目录。", true);
    return;
  }
  if (!absolutePath.startsWith(workspaceRoot)) {
    showToast("脚本文件必须位于部署基准目录内。", true);
    return;
  }
  const relative = absolutePath.slice(workspaceRoot.length).replace(/^[/\\\\]+/, "");
  elements.adminPipelineScript.value = relative || "";
}

function renderBrowserBreadcrumbs() {
  elements.browserBreadcrumbs.replaceChildren();
  const currentPath = state.pathBrowser.currentPath || state.pathBrowser.root || "/";
  const normalized = String(currentPath || "/");
  const parts = normalized.split("/").filter(Boolean);
  const crumbs = [];

  if (state.pathBrowser.mode === "project") {
    let current = normalized.startsWith("/") ? "/" : "";
    crumbs.push({ label: "/", value: "/" });
    parts.forEach((part) => {
      current = current === "/" ? `/${part}` : current ? `${current}/${part}` : part;
      const isProjectRoot = current === state.pathBrowser.root;
      crumbs.push({ label: isProjectRoot ? "项目根目录" : part, value: current });
    });
  } else {
    let current = normalized.startsWith("/") ? "/" : "";
    crumbs.push({ label: normalized.startsWith("/") ? "/" : "根目录", value: normalized.startsWith("/") ? "/" : state.pathBrowser.root || "" });
    parts.forEach((part) => {
      current = current === "/" ? `/${part}` : current ? `${current}/${part}` : part;
      crumbs.push({ label: part, value: current });
    });
  }

  crumbs.forEach((crumb, index) => {
    const button = document.createElement("button");
    button.type = "button";
    button.className = "browser-crumb";
    button.textContent = crumb.label;
    button.addEventListener("click", async () => {
      await loadBrowserDirectory(crumb.value);
    });
    elements.browserBreadcrumbs.appendChild(button);
    if (index < crumbs.length - 1) {
      const divider = document.createElement("span");
      divider.className = "browser-crumb-divider";
      divider.textContent = "/";
      elements.browserBreadcrumbs.appendChild(divider);
    }
  });
}

async function openBrowserShortcut(type) {
  if (type === "project") {
    state.pathBrowser.mode = "project";
    state.pathBrowser.root = "";
    await loadBrowserDirectory("");
    return;
  }
  if (type === "workspace") {
    const workspaceRoot = elements.adminWorkspaceRoot.value.trim();
    if (!workspaceRoot) {
      showToast("当前还没有部署基准目录。", true);
      return;
    }
    state.pathBrowser.mode = "admin";
    state.pathBrowser.root = workspaceRoot;
    await loadBrowserDirectory(workspaceRoot);
    return;
  }
  if (type === "home") {
    state.pathBrowser.mode = "admin";
    state.pathBrowser.root = "/";
    await loadBrowserDirectory("/Users");
    return;
  }
  if (type === "desktop") {
    state.pathBrowser.mode = "admin";
    state.pathBrowser.root = "/";
    await loadBrowserDirectory("/Users/wuhhh/Desktop");
  }
}

function syncBrowserShortcuts() {
  const current = state.pathBrowser.currentPath || state.pathBrowser.relativePath || "";
  const shortcuts = [
    [elements.browserShortcutProject, state.pathBrowser.mode === "project" && current === state.pathBrowser.root],
    [elements.browserShortcutWorkspace, elements.adminWorkspaceRoot.value.trim() && state.pathBrowser.mode === "admin" && current === elements.adminWorkspaceRoot.value.trim()],
    [elements.browserShortcutHome, state.pathBrowser.mode === "admin" && current === "/Users"],
    [elements.browserShortcutDesktop, state.pathBrowser.mode === "admin" && current === "/Users/wuhhh/Desktop"],
  ];
  shortcuts.forEach(([node, active]) => {
    if (!node) return;
    node.classList.toggle("active", Boolean(active));
  });
}
