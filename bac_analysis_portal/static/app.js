const state = {
  tasks: [],
  selectedTaskId: "",
  activeTab: "submission-tab",
  pollTimer: null,
  asmOptions: {},
  currentUser: null,
  taskDetailOpen: false,
  pathBrowser: {
    open: false,
    selector: "input",
    mode: "project",
    root: "",
    relativePath: "",
    items: [],
    selectedItem: null,
  },
};

const elements = {
  form: document.getElementById("task-form"),
  submitButton: document.getElementById("submit-button"),
  heroStartButton: document.getElementById("hero-start-button"),
  heroQueueButton: document.getElementById("hero-queue-button"),
  refreshButton: document.getElementById("refresh-button"),
  logoutButton: document.getElementById("logout-button"),
  currentUserChip: document.getElementById("current-user-chip"),
  tabButtons: Array.from(document.querySelectorAll("[data-tab-target]")),
  tabPanels: Array.from(document.querySelectorAll(".tab-panel")),
  adminTabButton: document.getElementById("admin-tab-button"),
  chooseInputPathButton: document.getElementById("choose-input-path"),
  chooseOutputDirButton: document.getElementById("choose-output-dir"),
  asmType: document.getElementById("asm_type"),
  asmTypeNote: document.getElementById("asm-type-note"),
  method: document.getElementById("method"),
  longType: document.getElementById("long_type"),
  longTypeWrap: document.getElementById("long-type-wrap"),
  taskList: document.getElementById("task-list"),
  queueSummary: document.getElementById("queue-summary"),
  taskDetail: document.getElementById("task-detail"),
  taskStatusChip: document.getElementById("task-status-chip"),
  taskDetailModal: document.getElementById("task-detail-modal"),
  taskDetailBackdrop: document.getElementById("task-detail-backdrop"),
  closeTaskDetailButton: document.getElementById("close-task-detail"),
  adminPanel: document.getElementById("admin-panel"),
  settingsForm: document.getElementById("settings-form"),
  adminWorkspaceRoot: document.getElementById("admin_workspace_root"),
  adminPipelineScript: document.getElementById("admin_pipeline_script"),
  chooseWorkspaceRootButton: document.getElementById("choose-workspace-root"),
  choosePipelineScriptButton: document.getElementById("choose-pipeline-script"),
  userForm: document.getElementById("user-form"),
  userList: document.getElementById("user-list"),
  pathBrowserModal: document.getElementById("path-browser-modal"),
  pathBrowserBackdrop: document.getElementById("path-browser-backdrop"),
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
  bindEvents();
  await loadSession();
  await loadAsmOptions();
  await loadTasks();
  if (state.currentUser?.role === "admin") {
    await Promise.all([loadAdminSettings(), loadUsers()]);
  }
  startPolling();
});

function bindEvents() {
  elements.form?.addEventListener("submit", onSubmitTask);
  elements.tabButtons.forEach((button) => {
    button.addEventListener("click", () => setActiveTab(button.dataset.tabTarget));
  });
  elements.heroStartButton?.addEventListener("click", () => {
    setActiveTab("submission-tab");
    elements.chooseInputPathButton?.focus();
  });
  elements.heroQueueButton?.addEventListener("click", () => setActiveTab("queue-tab"));
  elements.chooseInputPathButton?.addEventListener("click", () => openPathBrowser("input"));
  elements.chooseOutputDirButton?.addEventListener("click", () => openPathBrowser("output"));
  elements.asmType?.addEventListener("change", syncAsmTypeState);
  elements.refreshButton?.addEventListener("click", refreshQueueAndModal);
  elements.logoutButton?.addEventListener("click", logout);
  elements.settingsForm?.addEventListener("submit", onSaveSettings);
  elements.userForm?.addEventListener("submit", onCreateUser);
  elements.chooseWorkspaceRootButton?.addEventListener("click", () => openPathBrowser("workspace_root"));
  elements.choosePipelineScriptButton?.addEventListener("click", () => openPathBrowser("script_file"));
  elements.closeTaskDetailButton?.addEventListener("click", closeTaskDetail);
  elements.taskDetailBackdrop?.addEventListener("click", closeTaskDetail);
  elements.closePathBrowserButton?.addEventListener("click", closePathBrowser);
  elements.pathBrowserBackdrop?.addEventListener("click", closePathBrowser);
  elements.browserGoUpButton?.addEventListener("click", onBrowserGoUp);
  elements.browserNewFolderButton?.addEventListener("click", onBrowserNewFolder);
  elements.browserRenameButton?.addEventListener("click", onBrowserRename);
  elements.browserSelectCurrentButton?.addEventListener("click", onSelectCurrentDirectory);
  elements.browserShortcutProject?.addEventListener("click", () => openBrowserShortcut("project"));
  elements.browserShortcutHome?.addEventListener("click", () => openBrowserShortcut("home"));
  elements.browserShortcutDesktop?.addEventListener("click", () => openBrowserShortcut("desktop"));
  elements.browserShortcutWorkspace?.addEventListener("click", () => openBrowserShortcut("workspace"));
  document.addEventListener("keydown", (event) => {
    if (event.key === "Escape" && state.taskDetailOpen) {
      closeTaskDetail();
      return;
    }
    if (event.key === "Escape" && state.pathBrowser.open) {
      closePathBrowser();
    }
  });
}

async function loadSession() {
  const user = await requestJson("/api/session");
  state.currentUser = user;
  elements.currentUserChip.textContent = `${user.username} (${user.role === "admin" ? "管理员" : "普通用户"})`;
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
    state.selectedTaskId = created.id;
    await loadTasks();
    setActiveTab("queue-tab");
    await loadTaskDetail(created.id);
    openTaskDetail();
  } finally {
    elements.submitButton.disabled = false;
    elements.submitButton.textContent = "启动任务";
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

async function loadTasks() {
  const data = await requestJson("/api/tasks");
  state.tasks = data.items || [];
  renderTaskList();
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
  renderQueueSummary();
  if (!state.tasks.length) {
    const empty = document.createElement("div");
    empty.className = "empty-state queue-empty";
    empty.innerHTML = `
      <strong>当前没有任务</strong>
      <p>先在“提交任务”标签中填写输入路径和输出目录，再启动第一个分析任务。</p>
    `;
    elements.taskList.appendChild(empty);
    return;
  }

  state.tasks.forEach((task) => {
    const card = document.createElement("button");
    card.type = "button";
    card.className = "task-card";
    if (task.id === state.selectedTaskId) {
      card.classList.add("active");
    }
    card.innerHTML = `
      <div class="task-card-head">
        <strong>${escapeHtml(task.name || task.id)}</strong>
        <span class="mini-chip ${statusClassName(task.status)}">${escapeHtml(task.status || "-")}</span>
      </div>
      <dl class="task-card-meta">
        <div><dt>归属用户</dt><dd>${escapeHtml(task.owner || "-")}</dd></div>
        <div><dt>输入</dt><dd>${escapeHtml(task.params?.input_path || "-")}</dd></div>
        <div><dt>输出</dt><dd>${escapeHtml(task.params?.output_dir || "-")}</dd></div>
        <div><dt>创建</dt><dd>${escapeHtml(formatDate(task.created_at))}</dd></div>
      </dl>
    `;
    card.addEventListener("click", async () => {
      state.selectedTaskId = task.id;
      renderTaskList();
      await loadTaskDetail(task.id);
      openTaskDetail();
    });
    elements.taskList.appendChild(card);
  });
}

async function loadTaskDetail(taskId) {
  const task = await requestJson(`/api/tasks/${encodeURIComponent(taskId)}`);
  renderTaskDetail(task);
}

function renderTaskDetail(task) {
  state.selectedTaskId = task.id;
  elements.taskStatusChip.textContent = task.status || "未知";
  elements.taskStatusChip.className = `status-chip ${statusClassName(task.status)}`;
  const command = Array.isArray(task.command) ? task.command.join(" ") : "";
  const params = task.params || {};
  elements.taskDetail.classList.remove("empty-state");
  elements.taskDetail.innerHTML = `
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
    <div class="detail-block">
      <h3>日志尾部</h3>
      <p class="field-note">显示最近日志输出，适合快速确认当前进展和报错位置。</p>
      <pre>${escapeHtml(task.log_tail || "暂无日志输出")}</pre>
    </div>
  `;
  renderTaskList();
}

function openTaskDetail() {
  state.taskDetailOpen = true;
  elements.taskDetailModal.classList.remove("hidden");
  elements.taskDetailModal.setAttribute("aria-hidden", "false");
}

function closeTaskDetail() {
  state.taskDetailOpen = false;
  elements.taskDetailModal.classList.add("hidden");
  elements.taskDetailModal.setAttribute("aria-hidden", "true");
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
      <p>${escapeHtml(formatDate(user.created_at))}</p>
      <div class="user-card-actions">
        <button class="ghost-button" type="button">修改</button>
      </div>
    `;
    row.querySelector("button").addEventListener("click", async () => {
      const displayName = window.prompt("输入显示名", user.display_name || "");
      if (displayName === null) return;
      const role = window.prompt("输入角色：admin 或 user", user.role);
      if (role === null) return;
      const password = window.prompt("如需重置密码请输入新密码，留空则不修改", "");
      await requestJson(`/api/admin/users/${encodeURIComponent(user.username)}`, {
        method: "PUT",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          display_name: displayName,
          role: role,
          password: password || "",
        }),
      });
      showToast(`已更新用户：${user.username}`);
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
}

async function openPathBrowser(selector) {
  if (selector === "script_file" && !elements.adminWorkspaceRoot.value.trim()) {
    showToast("请先选择部署基准目录。", true);
    return;
  }
  state.pathBrowser.open = true;
  state.pathBrowser.selector = selector;
  state.pathBrowser.mode = ["workspace_root", "script_file"].includes(selector) ? "admin" : "project";
  state.pathBrowser.root = selector === "script_file"
    ? (elements.adminWorkspaceRoot.value.trim() || "")
    : "";
  state.pathBrowser.relativePath = "";
  state.pathBrowser.selectedItem = null;
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
  state.pathBrowser.relativePath = state.pathBrowser.mode === "admin"
    ? (data.current_path || "")
    : (data.relative_path || "");
  state.pathBrowser.items = data.items || [];
  state.pathBrowser.selectedItem = null;
  elements.browserCurrentPath.textContent = data.current_path || "/";
  const parentPath = state.pathBrowser.mode === "admin" ? (data.parent_path || "") : (data.parent_relative_path || "");
  elements.browserGoUpButton.disabled = !parentPath && !data.current_path;
  elements.browserGoUpButton.dataset.targetPath = parentPath || "";
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
  if (state.pathBrowser.selector === "script_file") {
    setAdminScriptPath(item.path);
    closePathBrowser();
  }
}

async function onBrowserGoUp() {
  await loadBrowserDirectory(elements.browserGoUpButton.dataset.targetPath || "");
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
  if (state.pathBrowser.selector === "input" && state.pathBrowser.selectedItem?.type === "file") {
    setSelectedPath("input_path", state.pathBrowser.selectedItem.path);
    closePathBrowser();
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
  elements.browserRenameButton.disabled = !selected;
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
    }),
  });
  elements.adminWorkspaceRoot.value = data.workspace_root || "";
  elements.adminPipelineScript.value = data.pipeline_script || "";
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

function renderQueueSummary() {
  const total = state.tasks.length;
  const running = state.tasks.filter((task) => task.status === "RUNNING").length;
  const queued = state.tasks.filter((task) => task.status === "QUEUED").length;
  const failed = state.tasks.filter((task) => task.status === "FAILED").length;
  const succeeded = state.tasks.filter((task) => task.status === "SUCCEEDED").length;
  elements.queueSummary.innerHTML = `
    <article><span>总任务</span><strong>${total}</strong></article>
    <article><span>运行中</span><strong>${running}</strong></article>
    <article><span>已完成</span><strong>${succeeded}</strong></article>
    <article><span>失败</span><strong>${failed}</strong></article>
    <article><span>排队中</span><strong>${queued}</strong></article>
  `;
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

function describeBrowserSelector(selector) {
  if (selector === "input") return "选择输入文件或输入目录。";
  if (selector === "output") return "选择输出目录。";
  if (selector === "workspace_root") return "选择部署基准目录。";
  if (selector === "script_file") return "从部署基准目录中选择脚本文件。";
  return "浏览路径。";
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
  const fullPath = state.pathBrowser.relativePath || state.pathBrowser.root || "/";
  const normalized = String(fullPath || "/");
  const isAbsolute = normalized.startsWith("/");
  const parts = normalized.split("/").filter(Boolean);
  const crumbs = [];

  if (state.pathBrowser.mode === "project") {
    crumbs.push({ label: "项目根目录", value: "" });
    let current = "";
    parts.forEach((part) => {
      current = current ? `${current}/${part}` : part;
      crumbs.push({ label: part, value: current });
    });
  } else {
    let current = isAbsolute ? "/" : "";
    crumbs.push({ label: isAbsolute ? "/" : "根目录", value: isAbsolute ? "/" : state.pathBrowser.root || "" });
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
  const current = state.pathBrowser.relativePath || "";
  const shortcuts = [
    [elements.browserShortcutProject, state.pathBrowser.mode === "project" && current === ""],
    [elements.browserShortcutWorkspace, elements.adminWorkspaceRoot.value.trim() && state.pathBrowser.mode === "admin" && current === elements.adminWorkspaceRoot.value.trim()],
    [elements.browserShortcutHome, state.pathBrowser.mode === "admin" && current === "/Users"],
    [elements.browserShortcutDesktop, state.pathBrowser.mode === "admin" && current === "/Users/wuhhh/Desktop"],
  ];
  shortcuts.forEach(([node, active]) => {
    if (!node) return;
    node.classList.toggle("active", Boolean(active));
  });
}
