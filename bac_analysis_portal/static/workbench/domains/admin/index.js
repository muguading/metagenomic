let runtime = null;
let initialized = false;
let nextcladeDatasetStatus = null;

export function init(nextRuntime) {
  runtime = nextRuntime;
  if (initialized) return;
  initialized = true;
  const elements = runtime.selectors.elements();
  elements.adminSectionTabs.forEach((button) => {
    button.addEventListener("click", () => activate(button.dataset.adminSection || "admin-settings-section"));
  });
  elements.settingsForm?.addEventListener("submit", runtime.actions.onSaveSettings);
  elements.userForm?.addEventListener("submit", runtime.actions.onCreateUser);
  elements.editUserForm?.addEventListener("submit", runtime.actions.onSaveUserEdit);
  elements.adminMonitorForm?.addEventListener("submit", runtime.actions.onSubmitAdminMonitor);
  elements.adminPathosourceTriggerForm?.addEventListener("submit", runtime.actions.onSaveAdminPathosourceTriggerRules);
  elements.checkNextcladeDatasetsButton?.addEventListener("click", onCheckNextcladeDatasets);
  elements.updateNextcladeDatasetsButton?.addEventListener("click", onUpdateNextcladeDatasets);
  elements.runOnlineUpdateButton?.addEventListener("click", runtime.actions.onRunOnlineUpdate);
  elements.runOfflineUpdateButton?.addEventListener("click", runtime.actions.onRunOfflineUpdate);
  bindPathBrowserTargets(elements);
}

export function activate(sectionId = "") {
  const remembered = sectionId || window.localStorage.getItem("bac-admin-section") || "admin-settings-section";
  const elements = runtime.selectors.elements();
  elements.adminSectionTabs.forEach((button) => {
    const active = button.dataset.adminSection === remembered;
    button.classList.toggle("active", active);
    button.setAttribute("aria-selected", active ? "true" : "false");
  });
  elements.adminSectionPanels.forEach((panel) => {
    panel.classList.toggle("active", panel.id === remembered);
  });
  window.localStorage.setItem("bac-admin-section", remembered);
  runtime.actions.scheduleSegmentedControlsSync(document.querySelector(".admin-section-tabs") || document);
}

export async function refresh() {
  await Promise.all([refreshSettings(true), refreshUsers()]);
}

export async function refreshSettings(force = false) {
  if (!force && runtime.actions.isAdminSettingsDirty()) return;
  const data = await runtime.api.requestJson("/api/admin/settings");
  const elements = runtime.selectors.elements();
  const state = runtime.selectors.state();
  elements.adminWorkspaceRoot.value = data.workspace_root || "";
  elements.adminPipelineScript.value = data.pipeline_script || "";
  if (elements.adminCondaRoot) elements.adminCondaRoot.value = data.conda_root || "";
  if (elements.adminDatabaseRoot) elements.adminDatabaseRoot.value = data.database_root || "";
  if (elements.adminMaxConcurrentTasks) elements.adminMaxConcurrentTasks.value = String(data.max_concurrent_tasks || 2);
  state.detectedCondaEnvs = Array.isArray(data.detected_conda_envs) ? data.detected_conda_envs : [];
  runtime.actions.rebuildAdminPipelineEnvSelectOptions(data.pipeline_python || "");
  runtime.actions.applyAdminCondaEnvValues(data.conda_envs || {});
  runtime.actions.renderAdminDetectedCondaEnvs();
  state.adminSettingsSnapshot = runtime.actions.getAdminSettingsDraft();
  if (elements.adminMonitorOutputRoot && !elements.adminMonitorOutputRoot.value.trim()) {
    elements.adminMonitorOutputRoot.value = data.workspace_root || "";
  }
  runtime.actions.applyAdminPathosourceTriggerRules(data.pathosource_trigger_rules || runtime.constants.adminPathosourceTriggerDefaults);
  runtime.actions.syncAdminMonitorModuleState();
  runtime.actions.renderAdminMonitorTaskList();
}

export async function refreshUsers() {
  const data = await runtime.api.requestJson("/api/admin/users");
  const state = runtime.selectors.state();
  state.adminUsers = data.items || [];
  runtime.actions.renderReviewerReadinessPanel();
  runtime.actions.renderUserList(state.adminUsers);
}

function setNextcladeDatasetPending(message) {
  const elements = runtime.selectors.elements();
  if (!elements.nextcladeDatasetResult) return;
  elements.nextcladeDatasetResult.className = "admin-update-result nextclade-dataset-result";
  elements.nextcladeDatasetResult.innerHTML = `
    <strong>${runtime.ui.escapeHtml(message)}</strong>
    <p>正在通过 ncov 环境读取本地版本和在线兼容版本。</p>
  `;
}

function nextcladeDatasetStatusLabel(status) {
  return {
    latest: "最新",
    outdated: "可更新",
    missing: "待下载",
    unmatched: "未匹配",
    invalid: "元数据异常",
    update_failed: "更新失败",
  }[status] || "未知";
}

function renderNextcladeDatasetStatus(data, isError = false) {
  const elements = runtime.selectors.elements();
  const escapeHtml = runtime.ui.escapeHtml;
  const target = elements.nextcladeDatasetResult;
  const updateButton = elements.updateNextcladeDatasetsButton;
  if (!target) return;
  if (isError) {
    nextcladeDatasetStatus = null;
    target.className = "admin-update-result nextclade-dataset-result error";
    target.innerHTML = `
      <strong>检测失败</strong>
      <p>${escapeHtml(data?.error || "Nextclade 数据库版本检测失败。")}</p>
    `;
    updateButton?.classList.add("hidden");
    if (updateButton) updateButton.disabled = true;
    return;
  }

  const items = Array.isArray(data?.items) ? data.items : [];
  const summary = data?.summary || {};
  const updatableItems = items.filter((item) => item?.updatable);
  const hasFailures = Number(summary.failed_count || 0) > 0;
  nextcladeDatasetStatus = data;
  target.className = `admin-update-result nextclade-dataset-result${hasFailures ? " error" : updatableItems.length ? "" : " success"}`;
  target.innerHTML = `
    <div class="nextclade-dataset-summary">
      <div>
        <strong>${updatableItems.length ? `发现 ${updatableItems.length} 个待获取或更新的数据集` : "本地数据集已完成检测"}</strong>
        <p>${escapeHtml(data?.database_root || "-")} · ${escapeHtml(data?.nextclade_version || "Nextclade")}</p>
      </div>
      <span>共 ${Number(summary.total_count || items.length)} · 最新 ${Number(summary.latest_count || 0)} · 待下载 ${Number(summary.missing_count || 0)} · 未匹配 ${Number(summary.unmatched_count || 0)} · 异常 ${Number(summary.failed_count || 0)}</span>
    </div>
    <div class="nextclade-dataset-table-shell">
      <table class="nextclade-dataset-table">
        <thead>
          <tr><th>数据集</th><th>本地版本</th><th>最新版本</th><th>状态</th></tr>
        </thead>
        <tbody>
          ${items.map((item) => `
            <tr>
              <td><strong>${escapeHtml(item.display_name || item.directory || "-")}</strong><small>${escapeHtml(item.directory || "-")}</small></td>
              <td><code>${escapeHtml(item.local_tag || "-")}</code></td>
              <td><code>${escapeHtml(item.latest_tag || "-")}</code></td>
              <td><span class="nextclade-dataset-status is-${escapeHtml(item.status || "unknown")}">${escapeHtml(nextcladeDatasetStatusLabel(item.status))}</span>${item.status === "latest" ? "" : `<small>${escapeHtml(item.message || "")}</small>`}</td>
            </tr>
          `).join("")}
        </tbody>
      </table>
    </div>
  `;
  updateButton?.classList.toggle("hidden", !updatableItems.length);
  if (updateButton) updateButton.disabled = !updatableItems.length;
}

async function onCheckNextcladeDatasets(event) {
  const button = event?.currentTarget;
  await runtime.actions.withSubmittingState(button, "检测中...", async () => {
    setNextcladeDatasetPending("正在检测 Nextclade 数据库版本");
    try {
      const data = await runtime.api.requestJson("/api/admin/nextclade-datasets");
      renderNextcladeDatasetStatus(data);
    } catch (error) {
      renderNextcladeDatasetStatus({ error: error.message || "检测失败" }, true);
    }
  });
}

async function onUpdateNextcladeDatasets(event) {
  const elements = runtime.selectors.elements();
  const updatableItems = (nextcladeDatasetStatus?.items || []).filter((item) => item?.updatable);
  if (!updatableItems.length) {
    runtime.ui.showToast("请先检测数据库版本。", "warning");
    return;
  }
  const confirmed = await runtime.actions.confirmDangerAction({
    title: "更新 Nextclade 数据库",
    message: `确认获取或更新 ${updatableItems.length} 个 Nextclade 数据集吗？`,
    impact: "缺失数据集将下载；已有数据集将替换为最新兼容版本",
    detail: updatableItems.map((item) => `${item.display_name || item.directory}：${item.local_tag || "-"} → ${item.latest_tag || "-"}`).join("；"),
    confirmLabel: "确认更新",
    tone: "warning",
  });
  if (!confirmed) return;

  const button = event?.currentTarget;
  const checkButton = elements.checkNextcladeDatasetsButton;
  await runtime.actions.withSubmittingState(button, "更新中...", async () => {
    if (checkButton) checkButton.disabled = true;
    setNextcladeDatasetPending("正在下载并替换 Nextclade 数据集");
    try {
      const data = await runtime.api.requestJson("/api/admin/nextclade-datasets/update", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ directories: updatableItems.map((item) => item.directory) }),
      });
      renderNextcladeDatasetStatus(data);
      if (Number(data.update_failed_count || 0)) {
        runtime.ui.showToast(`已更新 ${Number(data.updated_count || 0)} 个，失败 ${Number(data.update_failed_count || 0)} 个`, "warning");
      } else {
        runtime.ui.showToast(`已更新 ${Number(data.updated_count || 0)} 个 Nextclade 数据集`);
      }
    } catch (error) {
      renderNextcladeDatasetStatus({ error: error.message || "更新失败" }, true);
    } finally {
      if (checkButton) checkButton.disabled = false;
    }
  });
}

function bindPathBrowserTargets(elements) {
  const open = (selector) => () => runtime.pathBrowser.open(selector, { mode: "admin" });
  const buttonTargets = [
    [elements.chooseWorkspaceRootButton, "workspace_root"],
    [elements.choosePipelineScriptButton, "script_file"],
    [elements.chooseCondaRootButton, "conda_root"],
    [elements.chooseDatabaseRootButton, "database_root"],
    [elements.chooseAdminMonitorInputDirButton, "admin_monitor_input_dir"],
    [elements.chooseAdminMonitorOutputRootButton, "admin_monitor_output_root"],
    [elements.chooseAdminMonitorInputDirBacteriaButton, "admin_monitor_input_dir_bacteria"],
    [elements.chooseAdminMonitorInputDirVirusButton, "admin_monitor_input_dir_virus"],
    [elements.chooseAdminMonitorInputDirMetagenomeButton, "admin_monitor_input_dir_metagenome"],
    [elements.chooseAdminMonitorOutputRootBacteriaButton, "admin_monitor_output_root_bacteria"],
    [elements.chooseAdminMonitorOutputRootVirusButton, "admin_monitor_output_root_virus"],
    [elements.chooseAdminMonitorOutputRootMetagenomeButton, "admin_monitor_output_root_metagenome"],
    [elements.chooseOfflineUpdateSourceButton, "offline_update_source"],
  ];
  buttonTargets.forEach(([button, selector]) => button?.addEventListener("click", open(selector)));

  const directoryTargets = new Map([
    ["database_root", elements.adminDatabaseRoot],
    ["admin_monitor_input_dir", elements.adminMonitorInputDir],
    ["admin_monitor_output_root", elements.adminMonitorOutputRoot],
    ["admin_monitor_input_dir_bacteria", elements.adminMonitorModuleInputBacteria],
    ["admin_monitor_input_dir_virus", elements.adminMonitorModuleInputVirus],
    ["admin_monitor_input_dir_metagenome", elements.adminMonitorModuleInputMetagenome],
    ["admin_monitor_output_root_bacteria", elements.adminMonitorModuleOutputBacteria],
    ["admin_monitor_output_root_virus", elements.adminMonitorModuleOutputVirus],
    ["admin_monitor_output_root_metagenome", elements.adminMonitorModuleOutputMetagenome],
    ["offline_update_source", elements.offlineUpdateSource],
  ]);
  directoryTargets.forEach((input, selector) => {
    runtime.pathBrowser.registerTarget(selector, ({ currentPath }) => {
      if (input) input.value = currentPath || "";
      runtime.actions.rememberPathBrowserLocation(selector, currentPath || "");
      runtime.actions.closePathBrowser();
    });
  });
  runtime.pathBrowser.registerTarget("workspace_root", ({ currentPath }) => {
    elements.adminWorkspaceRoot.value = currentPath || "";
    elements.adminPipelineScript.value = "";
    runtime.actions.rememberPathBrowserLocation("workspace_root", currentPath || "");
    runtime.actions.closePathBrowser();
  });
  runtime.pathBrowser.registerTarget("conda_root", async ({ currentPath }) => {
    if (elements.adminCondaRoot) elements.adminCondaRoot.value = currentPath || "";
    runtime.actions.rememberPathBrowserLocation("conda_root", currentPath || "");
    runtime.actions.closePathBrowser();
    await runtime.actions.refreshAdminCondaEnvsForRoot({ announce: true });
  });
  runtime.pathBrowser.registerTarget("script_file", ({ selectedItem }) => {
    if (selectedItem?.type !== "file") return;
    runtime.actions.setAdminScriptPath(selectedItem.path);
    runtime.actions.closePathBrowser();
  });
  runtime.pathBrowser.registerTarget("pipeline_python", ({ selectedItem }) => {
    if (selectedItem?.type !== "file") return;
    runtime.actions.setAdminPipelinePython(selectedItem.path);
    runtime.actions.closePathBrowser();
  });
}
