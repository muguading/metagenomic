const state = {
  page: 1,
  pageSize: 8,
  totalPages: 1,
  databasePage: 1,
  databasePageSize: 12,
  databaseTotalPages: 1,
  auditPage: 1,
  auditPageSize: 10,
  auditTotalPages: 1,
  selectedGenomeId: null,
  mode: "create",
  currentUser: null,
  activeTab: "genome-tab",
  activeFormTab: "basic-tab",
  metadataTemplates: [],
  searchMetadataTemplates: [],
  monitoringData: [],
  monitoringTemplates: [],
  monitoringWidgets: [],
  monitoringDrag: null,
  surveillanceMapScope: "china",
  surveillanceMapDrilldownCode: "",
  surveillanceMapDrilldownName: "",
  surveillanceMapParentCode: "",
  surveillanceMapParentName: "",
};

let COUNTRY_OPTIONS = [];
let CHINA_REGION_TREE = {};
const MONITORING_MAX_WIDGETS = 10;
const MONITORING_COLUMNS = 12;
const MONITORING_ROW_HEIGHT = 42;
const MONITORING_DEFAULT_WIDGET = { w: 4, h: 6 };
let WORLD_GEOJSON = null;
let CHINA_GEOJSON = null;
const monitoringLeafletMaps = new Map();
const CHINA_DRILLDOWN_CACHE = new Map();

const elements = {
  genomeForm: document.getElementById("genome-form"),
  resetFormButton: document.getElementById("reset-form"),
  submitButton: document.getElementById("submit-button"),
  deleteButton: document.getElementById("delete-button"),
  fastaFileInput: document.getElementById("fasta_file"),
  selectFileButton: document.getElementById("select-file-button"),
  uploadDropzone: document.getElementById("upload-dropzone"),
  uploadStatus: document.getElementById("upload-status"),
  formTitle: document.getElementById("form-title"),
  formTabButtons: Array.from(document.querySelectorAll("[data-form-tab-target]")),
  formTabPanels: Array.from(document.querySelectorAll(".form-tab-panel")),
  genomeModal: document.getElementById("genome-modal"),
  bulkImportModal: document.getElementById("bulk-import-modal"),
  openGenomeModalButton: document.getElementById("open-genome-modal"),
  openBulkImportModalButton: document.getElementById("open-bulk-import-modal"),
  closeGenomeModalButton: document.getElementById("close-genome-modal"),
  closeBulkImportModalButton: document.getElementById("close-bulk-import-modal"),
  bulkImportFile: document.getElementById("bulk-import-file"),
  submitBulkImportButton: document.getElementById("submit-bulk-import"),
  bulkImportStatus: document.getElementById("bulk-import-status"),
  bulkImportResult: document.getElementById("bulk-import-result"),
  detailCard: document.getElementById("detail-card"),
  detailStatus: document.getElementById("detail-status"),
  genomeList: document.getElementById("genome-list"),
  searchForm: document.getElementById("search-form"),
  refreshButton: document.getElementById("refresh-button"),
  toggleSearchPanelButton: document.getElementById("toggle-search-panel"),
  addSearchFilterButton: document.getElementById("add-search-filter"),
  searchCustomLogic: document.getElementById("search_custom_logic"),
  searchFilterList: document.getElementById("search-filter-list"),
  searchFilterEmptyState: document.getElementById("search-filter-empty-state"),
  refreshMonitoring: document.getElementById("refresh-monitoring"),
  resetMonitoringLayout: document.getElementById("reset-monitoring-layout"),
  addMonitoringWidget: document.getElementById("add-monitoring-widget"),
  createMonitoringWidget: document.getElementById("create-monitoring-widget"),
  monitoringTitle: document.getElementById("monitoring_title"),
  monitoringChartType: document.getElementById("monitoring_chart_type"),
  monitoringField: document.getElementById("monitoring_field"),
  monitoringMapScopeWrap: document.getElementById("monitoring_map_scope_wrap"),
  monitoringMapScope: document.getElementById("monitoring_map_scope"),
  surveillanceChinaButton: document.getElementById("surveillance-china-button"),
  surveillanceWorldButton: document.getElementById("surveillance-world-button"),
  surveillanceKpis: document.getElementById("surveillance-kpis"),
  surveillanceMap: document.getElementById("surveillance-map"),
  surveillanceMapBack: document.getElementById("surveillance-map-back"),
  surveillanceMapCaption: document.getElementById("surveillance-map-caption"),
  surveillanceTimeline: document.getElementById("surveillance-timeline"),
  surveillanceCountry: document.getElementById("surveillance-country"),
  surveillanceSampleType: document.getElementById("surveillance-sample-type"),
  surveillanceSequencing: document.getElementById("surveillance-sequencing"),
  surveillanceRecent: document.getElementById("surveillance-recent"),
  monitoringCount: document.getElementById("monitoring-count"),
  monitoringBoard: document.getElementById("monitoring-board"),
  monitoringEmptyState: document.getElementById("monitoring-empty-state"),
  refreshDatabase: document.getElementById("refresh-database"),
  prevDatabasePage: document.getElementById("prev-database-page"),
  nextDatabasePage: document.getElementById("next-database-page"),
  databasePageIndicator: document.getElementById("database-page-indicator"),
  databaseCount: document.getElementById("database-count"),
  databaseTableBody: document.getElementById("database-table-body"),
  databaseEmptyState: document.getElementById("database-empty-state"),
  prevPage: document.getElementById("prev-page"),
  nextPage: document.getElementById("next-page"),
  pageIndicator: document.getElementById("page-indicator"),
  genomeCount: document.getElementById("genome-count"),
  auditForm: document.getElementById("audit-form"),
  auditCount: document.getElementById("audit-count"),
  auditList: document.getElementById("audit-list"),
  refreshAudit: document.getElementById("refresh-audit"),
  prevAuditPage: document.getElementById("prev-audit-page"),
  nextAuditPage: document.getElementById("next-audit-page"),
  auditPageIndicator: document.getElementById("audit-page-indicator"),
  logoutButton: document.getElementById("logout-button"),
  currentUserChip: document.getElementById("current-user-chip"),
  tabButtons: Array.from(document.querySelectorAll("[data-tab-target]")),
  tabPanels: Array.from(document.querySelectorAll(".tab-panel")),
  adminTabButton: document.getElementById("admin-tab-button"),
  profileForm: document.getElementById("profile-form"),
  refreshProfile: document.getElementById("refresh-profile"),
  passwordForm: document.getElementById("password-form"),
  adminPanel: document.getElementById("admin-panel"),
  userForm: document.getElementById("user-form"),
  refreshUsers: document.getElementById("refresh-users"),
  userList: document.getElementById("user-list"),
  addMetaFieldButton: document.getElementById("add-meta-field"),
  metaFieldList: document.getElementById("meta-field-list"),
  metaEmptyState: document.getElementById("meta-empty-state"),
  toast: document.getElementById("toast"),
  genomeCardTemplate: document.getElementById("genome-card-template"),
  auditItemTemplate: document.getElementById("audit-item-template"),
  userItemTemplate: document.getElementById("user-item-template"),
  metaFieldTemplate: document.getElementById("meta-field-template"),
  searchFilterTemplate: document.getElementById("search-filter-template"),
  monitoringWidgetTemplate: document.getElementById("monitoring-widget-template"),
};

document.addEventListener("DOMContentLoaded", async () => {
  try {
    await loadReferenceData();
    await loadSession();
    loadMonitoringWidgetsFromStorage();
    bindEvents();
    await Promise.all([loadMetadataTemplates(), loadSearchMetadataTemplates(), loadGenomes(), loadAuditLogs()]);
    if (state.currentUser.role === "admin") {
      await loadUsers();
    }
  } catch (error) {
    showToast(error.message || "Failed to initialize page", true);
  }
});

async function loadReferenceData() {
  try {
    const [countries, chinaRegions] = await Promise.all([
      requestJson("/static/data/countries.json"),
      requestJson("/static/data/china_regions.json"),
    ]);
    COUNTRY_OPTIONS = Array.isArray(countries) ? countries : [];
    CHINA_REGION_TREE = chinaRegions && typeof chinaRegions === "object" ? chinaRegions : {};
  } catch (error) {
    COUNTRY_OPTIONS = [];
    CHINA_REGION_TREE = {};
    console.error("Failed to load reference data", error);
    showToast("Country and China region reference data failed to load.", true);
  }
  try {
    WORLD_GEOJSON = await requestJson("/static/data/world_countries.geojson");
  } catch (error) {
    WORLD_GEOJSON = null;
    console.error("Failed to load world map boundaries", error);
  }
  try {
    CHINA_GEOJSON = await requestJson("/static/data/china_full.geojson");
  } catch (error) {
    CHINA_GEOJSON = null;
    console.error("Failed to load China map boundaries", error);
  }
}

function bindEvents() {
  elements.tabButtons.forEach((button) => {
    button.addEventListener("click", () => setActiveTab(button.dataset.tabTarget));
  });
  elements.formTabButtons.forEach((button) => {
    button.addEventListener("click", () => setActiveFormTab(button.dataset.formTabTarget));
  });
  elements.openGenomeModalButton.addEventListener("click", async () => {
    await loadMetadataTemplates(state.currentUser?.username);
    resetForm({ preserveDetail: true });
    openGenomeModal();
  });
  elements.openBulkImportModalButton?.addEventListener("click", openBulkImportModal);
  elements.closeGenomeModalButton.addEventListener("click", closeGenomeModal);
  elements.closeBulkImportModalButton?.addEventListener("click", closeBulkImportModal);
  elements.genomeModal.addEventListener("click", (event) => {
    if (event.target.dataset.closeModal === "genome") {
      closeGenomeModal();
    }
  });
  elements.bulkImportModal?.addEventListener("click", (event) => {
    if (event.target.dataset.closeModal === "bulk-import") {
      closeBulkImportModal();
    }
  });
  elements.genomeForm.addEventListener("submit", onSubmitGenomeForm);
  elements.bulkImportFile?.addEventListener("change", onBulkImportFileSelected);
  elements.submitBulkImportButton?.addEventListener("click", onSubmitBulkImport);
  elements.resetFormButton.addEventListener("click", resetForm);
  elements.deleteButton.addEventListener("click", onDeleteGenome);
  elements.selectFileButton.addEventListener("click", () => elements.fastaFileInput.click());
  elements.fastaFileInput.addEventListener("change", onFileSelected);
  document.getElementById("country").addEventListener("change", () => syncBuiltinLocationField({ preserveValue: false }));
  document.getElementById("builtin_location_province").addEventListener("change", () => {
    populateBuiltinCityOptions(document.getElementById("builtin_location_province").value, "");
    populateBuiltinDistrictOptions(
      document.getElementById("builtin_location_province").value,
      document.getElementById("builtin_location_city").value,
      ""
    );
  });
  document.getElementById("builtin_location_city").addEventListener("change", () => {
    populateBuiltinDistrictOptions(
      document.getElementById("builtin_location_province").value,
      document.getElementById("builtin_location_city").value,
      ""
    );
  });
  elements.uploadDropzone.addEventListener("click", (event) => {
    if (event.target.tagName !== "BUTTON") {
      elements.fastaFileInput.click();
    }
  });
  elements.uploadDropzone.addEventListener("dragover", onDragOver);
  elements.uploadDropzone.addEventListener("dragleave", onDragLeave);
  elements.uploadDropzone.addEventListener("drop", onDropFile);
  elements.uploadDropzone.addEventListener("keydown", (event) => {
    if (event.key === "Enter" || event.key === " ") {
      event.preventDefault();
      elements.fastaFileInput.click();
    }
  });
  elements.searchForm.addEventListener("submit", onSearch);
  elements.toggleSearchPanelButton.addEventListener("click", toggleSearchPanel);
  elements.addSearchFilterButton.addEventListener("click", () => {
    renderSearchFilter();
    updateSearchFilterEmptyState();
  });
  elements.refreshMonitoring.addEventListener("click", async () => {
    await loadMonitoringData();
    renderSurveillanceScreen();
    renderMonitoringBoard();
  });
  elements.resetMonitoringLayout.addEventListener("click", resetMonitoringLayout);
  elements.addMonitoringWidget.addEventListener("click", addMonitoringWidgetFromBuilder);
  elements.createMonitoringWidget.addEventListener("click", addMonitoringWidgetFromBuilder);
  elements.monitoringChartType.addEventListener("change", syncMonitoringBuilderState);
  elements.surveillanceChinaButton?.addEventListener("click", () => {
    state.surveillanceMapScope = "china";
    resetSurveillanceDrilldown();
    renderSurveillanceScreen();
  });
  elements.surveillanceWorldButton?.addEventListener("click", () => {
    state.surveillanceMapScope = "world";
    resetSurveillanceDrilldown();
    renderSurveillanceScreen();
  });
  elements.surveillanceMapBack?.addEventListener("click", () => {
    stepBackSurveillanceDrilldown();
    renderSurveillanceScreen();
  });
  elements.refreshButton.addEventListener("click", () => loadGenomes());
  elements.refreshDatabase.addEventListener("click", () => loadDatabaseTable());
  elements.prevDatabasePage.addEventListener("click", () => changeDatabasePage(-1));
  elements.nextDatabasePage.addEventListener("click", () => changeDatabasePage(1));
  elements.prevPage.addEventListener("click", () => changePage(-1));
  elements.nextPage.addEventListener("click", () => changePage(1));
  elements.auditForm.addEventListener("submit", onAuditSearch);
  elements.refreshAudit.addEventListener("click", () => loadAuditLogs());
  elements.prevAuditPage.addEventListener("click", () => changeAuditPage(-1));
  elements.nextAuditPage.addEventListener("click", () => changeAuditPage(1));
  elements.logoutButton.addEventListener("click", logout);
  elements.profileForm.addEventListener("submit", onUpdateProfile);
  elements.refreshProfile.addEventListener("click", () => loadProfile());
  elements.passwordForm.addEventListener("submit", onChangePassword);
  elements.addMetaFieldButton.addEventListener("click", () => {
    renderMetaField(createEmptyMetaField());
    updateMetaEmptyState();
  });
  document.addEventListener("keydown", onGlobalKeydown);
  document.addEventListener("pointermove", onMonitoringPointerMove);
  document.addEventListener("pointerup", onMonitoringPointerUp);
  window.addEventListener("resize", () => {
    if (state.activeTab === "monitoring-tab") {
      renderSurveillanceScreen();
      renderMonitoringBoard();
    }
  });
  if (elements.userForm) {
    elements.userForm.addEventListener("submit", onCreateUser);
    elements.refreshUsers.addEventListener("click", () => loadUsers());
  }
  document.getElementById("search_submitter").addEventListener("change", onSearchSubmitterChanged);
  syncMonitoringBuilderState();
}

function onGlobalKeydown(event) {
  if (event.key === "Escape" && !elements.genomeModal.classList.contains("hidden")) {
    closeGenomeModal();
    return;
  }
  if (event.key === "Escape" && elements.bulkImportModal && !elements.bulkImportModal.classList.contains("hidden")) {
    closeBulkImportModal();
    return;
  }
  if (event.key === "Escape") {
    const expandedWidget = state.monitoringWidgets.find((item) => item.expanded);
    if (expandedWidget) {
      expandedWidget.expanded = false;
      saveMonitoringWidgetsToStorage();
      renderMonitoringBoard();
    }
  }
}

async function loadSession() {
  const user = await requestJson("/api/session");
  applyCurrentUser(user);
  fillProfileForm(user);
}

function applyCurrentUser(user) {
  state.currentUser = user;
  elements.currentUserChip.textContent = `${user.username} (${user.role})`;
  document.getElementById("submitter").value = user.username;
  populateCountryOptions(document.getElementById("country"), document.getElementById("country")?.value || "");
  syncBuiltinLocationField({ preserveValue: true });
  if (user.role === "admin") {
    elements.adminPanel.classList.remove("hidden");
    elements.adminTabButton.classList.remove("hidden");
    document.getElementById("search_submitter").placeholder = "alice";
  } else {
    elements.adminPanel.classList.add("hidden");
    elements.adminTabButton.classList.add("hidden");
    if (state.activeTab === "users-tab") {
      setActiveTab("genome-tab");
    }
    document.getElementById("search_submitter").value = user.username;
    document.getElementById("search_submitter").readOnly = true;
  }
}

function setActiveTab(tabId) {
  state.activeTab = tabId;
  elements.tabButtons.forEach((button) => {
    const isActive = button.dataset.tabTarget === tabId;
    button.classList.toggle("active", isActive);
    button.setAttribute("aria-selected", isActive ? "true" : "false");
  });
  elements.tabPanels.forEach((panel) => {
    panel.classList.toggle("active", panel.id === tabId);
  });
  if (tabId === "database-tab") {
    loadDatabaseTable();
  }
  if (tabId === "monitoring-tab") {
    loadMonitoringData().then(() => {
      populateMonitoringFieldOptions();
      renderSurveillanceScreen();
      renderMonitoringBoard();
    }).catch((error) => {
      showToast(error.message || "Failed to load monitoring dashboard", true);
    });
  }
}

function setActiveFormTab(tabId) {
  state.activeFormTab = tabId;
  elements.formTabButtons.forEach((button) => {
    const isActive = button.dataset.formTabTarget === tabId;
    button.classList.toggle("active", isActive);
  });
  elements.formTabPanels.forEach((panel) => {
    panel.classList.toggle("active", panel.id === tabId);
  });
}

async function loadProfile() {
  const user = await requestJson("/api/profile");
  applyCurrentUser(user);
  fillProfileForm(user);
}

async function loadMetadataTemplates(submitter = state.currentUser?.username) {
  const params = new URLSearchParams();
  if (submitter) params.set("submitter", submitter);
  const suffix = params.toString() ? `?${params.toString()}` : "";
  const data = await requestJson(`/api/metadata-templates${suffix}`);
  state.metadataTemplates = data.items || [];
}

async function loadSearchMetadataTemplates(submitter = getSearchTemplateSubmitter()) {
  const params = new URLSearchParams();
  if (submitter) params.set("submitter", submitter);
  const suffix = params.toString() ? `?${params.toString()}` : "";
  const data = await requestJson(`/api/metadata-templates${suffix}`);
  state.searchMetadataTemplates = getBuiltinStandardFields().concat(data.items || []);
  rerenderSearchFilters();
}

async function loadMonitoringData(submitter = getMonitoringSubmitter()) {
  const params = new URLSearchParams();
  if (submitter) params.set("submitter", submitter);
  const suffix = params.toString() ? `?${params.toString()}` : "";
  const data = await requestJson(`/api/dashboard-data${suffix}`);
  state.monitoringData = data.items || [];
  state.monitoringTemplates = data.templates || [];
  populateMonitoringFieldOptions();
}

function fillProfileForm(user) {
  document.getElementById("profile_username").value = user.username || "";
  document.getElementById("profile_display_name").value = user.display_name || "";
  document.getElementById("profile_email").value = user.email || "";
}

async function logout() {
  await requestJson("/logout", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
  });
  window.location.href = "/login";
}

async function loadGenomes() {
  const params = new URLSearchParams({
    page: String(state.page),
    page_size: String(state.pageSize),
  });
  const speciesName = document.getElementById("search_species_name").value.trim();
  const taxid = document.getElementById("search_taxid").value.trim();
  const submitter = document.getElementById("search_submitter").value.trim();
  const customLogic = elements.searchCustomLogic?.value || "and";
  const customFilters = collectSearchFilters();
  if (speciesName) params.set("species_name", speciesName);
  if (taxid) params.set("taxid", taxid);
  if (submitter) params.set("submitter", submitter);
  if (customFilters.length) {
    params.set("custom_logic", customLogic);
    params.set("custom_filters", JSON.stringify(customFilters));
  }

  const data = await requestJson(`/api/genomes?${params.toString()}`);
  state.totalPages = Math.max(1, data.pages || 1);
  renderGenomeList(data.items, data.total);
}

async function loadDatabaseTable() {
  const params = new URLSearchParams({
    page: String(state.databasePage),
    page_size: String(state.databasePageSize),
  });
  const submitter = state.currentUser?.role === "admin" ? "" : state.currentUser?.username || "";
  if (submitter) {
    params.set("submitter", submitter);
  }
  const data = await requestJson(`/api/genomes?${params.toString()}`);
  state.databaseTotalPages = Math.max(1, data.pages || 1);
  renderDatabaseTable(data.items || [], data.total || 0);
}

async function loadAuditLogs() {
  const params = new URLSearchParams({
    page: String(state.auditPage),
    page_size: String(state.auditPageSize),
  });
  const genomeId = document.getElementById("audit_genome_id").value.trim();
  const operation = document.getElementById("audit_operation").value.trim();
  if (genomeId) params.set("genome_id", genomeId);
  if (operation) params.set("operation", operation);
  const data = await requestJson(`/api/audit-logs?${params.toString()}`);
  state.auditTotalPages = Math.max(1, data.pages || 1);
  renderAuditList(data.items, data.total);
}

async function loadUsers() {
  if (state.currentUser.role !== "admin") return;
  const data = await requestJson("/api/users");
  renderUserList(data.items);
}

function renderGenomeList(items, total) {
  elements.genomeList.replaceChildren();
  elements.genomeCount.textContent = `${total} genomes`;
  elements.pageIndicator.textContent = `Page ${state.page} / ${state.totalPages}`;
  elements.prevPage.disabled = state.page <= 1;
  elements.nextPage.disabled = state.page >= state.totalPages;
  if (!items.length) {
    const empty = document.createElement("div");
    empty.className = "detail-card empty-state";
    empty.textContent = "No genomes matched the current filters.";
    elements.genomeList.appendChild(empty);
    return;
  }
  items.forEach((item) => {
    const node = elements.genomeCardTemplate.content.firstElementChild.cloneNode(true);
    node.dataset.genomeRecordId = String(item.id);
    node.querySelector("h3").textContent = item.genome_id;
    node.querySelector(".species-tag").textContent = item.species_name;
    node.querySelector(".genome-card-meta").textContent = `${item.sample_name} | TaxID ${item.taxid} | ${item.submitter}`;
    node.querySelector(".mini-stats").innerHTML = `
      <div><dt>Length</dt><dd>${formatNumber(item.genome_length)}</dd></div>
      <div><dt>Submitted</dt><dd>${formatDate(item.submit_time)}</dd></div>
    `;
    node.addEventListener("click", () => selectGenome(item.id));
    node.querySelector(".genome-edit-button").addEventListener("click", (event) => {
      event.stopPropagation();
      editGenome(item.id);
    });
    if (state.selectedGenomeId === item.id) node.classList.add("active");
    elements.genomeList.appendChild(node);
  });
}

function renderDatabaseTable(items, total) {
  elements.databaseTableBody.replaceChildren();
  elements.databaseCount.textContent = `${total} samples`;
  elements.databasePageIndicator.textContent = `Page ${state.databasePage} / ${state.databaseTotalPages}`;
  elements.prevDatabasePage.disabled = state.databasePage <= 1;
  elements.nextDatabasePage.disabled = state.databasePage >= state.databaseTotalPages;
  elements.databaseEmptyState.classList.toggle("hidden", items.length > 0);
  if (!items.length) {
    return;
  }
  items.forEach((item) => {
    const row = document.createElement("tr");
    row.innerHTML = `
      <td>${escapeHtml(item.genome_id)}</td>
      <td>${escapeHtml(item.sample_name)}</td>
      <td>${escapeHtml(item.species_name)}</td>
      <td>${escapeHtml(item.taxid)}</td>
      <td>${escapeHtml(item.submitter)}</td>
      <td>${escapeHtml(formatNumber(item.genome_length))}</td>
      <td>${escapeHtml(formatDate(item.submit_time))}</td>
    `;
    row.addEventListener("click", async () => {
      setActiveTab("genome-tab");
      await selectGenome(item.id);
    });
    elements.databaseTableBody.appendChild(row);
  });
}

function getMonitoringSubmitter() {
  if (state.currentUser?.role !== "admin") {
    return state.currentUser?.username || "";
  }
  return "";
}

function getMonitoringStorageKey() {
  return `genome_db_monitoring_layout::${state.currentUser?.username || "anonymous"}`;
}

function loadMonitoringWidgetsFromStorage() {
  try {
    const raw = window.localStorage.getItem(getMonitoringStorageKey());
    const parsed = raw ? JSON.parse(raw) : [];
    state.monitoringWidgets = Array.isArray(parsed) ? parsed.slice(0, MONITORING_MAX_WIDGETS) : [];
  } catch {
    state.monitoringWidgets = [];
  }
}

function getBuiltinStandardFields() {
  return [
    { key: "standard:species_name", label: "Species Name", type: "text" },
    { key: "standard:sample_name", label: "Sample Name", type: "text" },
    { key: "standard:submitter", label: "Submitter", type: "text" },
    { key: "standard:taxid", label: "TaxID", type: "text" },
    { key: "standard:genome_length", label: "Genome Length", type: "text" },
    { key: "standard:gender", label: "Gender", type: "select" },
    { key: "standard:country", label: "Country", type: "country" },
    { key: "standard:location", label: "Location", type: "location" },
    { key: "standard:collection_time", label: "Collection Time", type: "datetime" },
    { key: "standard:sample_type", label: "Sample Type", type: "text" },
    { key: "standard:sequencing_method", label: "Sequencing Method", type: "text" },
    { key: "standard:submit_time", label: "Submit Time", type: "datetime" },
    { key: "standard:last_modified_time", label: "Last Modified", type: "datetime" },
  ];
}

function saveMonitoringWidgetsToStorage() {
  window.localStorage.setItem(getMonitoringStorageKey(), JSON.stringify(state.monitoringWidgets.slice(0, MONITORING_MAX_WIDGETS)));
}

function populateMonitoringFieldOptions() {
  const options = getMonitoringFieldOptions();
  elements.monitoringField.innerHTML = options
    .map((item) => `<option value="${escapeHtml(item.key)}">${escapeHtml(item.label)}</option>`)
    .join("");
}

function getMonitoringFieldOptions() {
  const standardFields = getBuiltinStandardFields();
  const templateFields = (state.monitoringTemplates || []).map((item) => ({
    key: `meta:${item.key}`,
    label: item.label || item.key || "Metadata",
    type: item.type || "text",
    metaKey: item.key,
  }));
  return standardFields.concat(templateFields);
}

function renderSurveillanceScreen() {
  if (!elements.surveillanceKpis) return;
  destroyLeafletMapByKey("surveillance-screen");
  const total = state.monitoringData.length;
  const countries = new Set(state.monitoringData.map((item) => String(item.country || "").trim()).filter(Boolean));
  const chinaRegions = new Set(
    state.monitoringData
      .map((item) => item.location?.province || "")
      .map((value) => String(value || "").trim())
      .filter(Boolean)
  );
  const latestDate = getLatestMonitoringDate();
  const sampledTypes = new Set(state.monitoringData.map((item) => String(item.sample_type || "").trim()).filter(Boolean));
  const kpis = [
    { label: "Total Samples", value: formatNumber(total), note: "Records in current monitoring scope" },
    { label: "Countries", value: formatNumber(countries.size), note: "Distinct country coverage" },
    { label: "China Regions", value: formatNumber(chinaRegions.size), note: "Distinct province-level coverage" },
    { label: "Latest Activity", value: latestDate ? latestDate.slice(0, 10) : "-", note: `${sampledTypes.size} sample types observed` },
  ];
  elements.surveillanceKpis.innerHTML = kpis.map((item) => `
    <article class="surveillance-kpi">
      <span class="surveillance-kpi-label">${escapeHtml(item.label)}</span>
      <strong class="surveillance-kpi-value">${escapeHtml(item.value)}</strong>
      <div class="surveillance-kpi-note">${escapeHtml(item.note)}</div>
    </article>
  `).join("");

  const timelineSeries = buildSurveillanceTimelineSeries();
  const countrySeries = buildMonitoringSeries(getMonitoringFieldDefinition("standard:country"), "bar");
  const sampleTypeSeries = buildMonitoringSeries(getMonitoringFieldDefinition("standard:sample_type"), "pie");
  const sequencingSeries = buildMonitoringSeries(getMonitoringFieldDefinition("standard:sequencing_method"), "bar");
  elements.surveillanceTimeline.innerHTML = timelineSeries.labels.length
    ? renderLineChart(timelineSeries, "Collection Timeline")
    : '<div class="empty-state">No temporal data available.</div>';
  elements.surveillanceCountry.innerHTML = countrySeries.labels.length
    ? renderBarChart(countrySeries, "Country Distribution")
    : '<div class="empty-state">No country data available.</div>';
  elements.surveillanceSampleType.innerHTML = sampleTypeSeries.labels.length
    ? renderPieChart(sampleTypeSeries, "Sample Type")
    : '<div class="empty-state">No sample type data available.</div>';
  elements.surveillanceSequencing.innerHTML = sequencingSeries.labels.length
    ? renderBarChart(sequencingSeries, "Sequencing Method")
    : '<div class="empty-state">No sequencing method data available.</div>';
  renderSurveillanceRecentTable();
  renderSurveillanceMap();
  elements.surveillanceChinaButton?.classList.toggle("active", state.surveillanceMapScope === "china");
  elements.surveillanceWorldButton?.classList.toggle("active", state.surveillanceMapScope === "world");
  elements.surveillanceMapBack?.classList.toggle("hidden", !state.surveillanceMapDrilldownCode);
  if (elements.surveillanceMapCaption) {
    elements.surveillanceMapCaption.textContent = state.surveillanceMapScope === "china"
      ? describeChinaDrilldownCaption(
          state.surveillanceMapDrilldownName,
          state.surveillanceMapParentName,
        )
      : "Country-level heat distribution. Use trackpad or mouse wheel to zoom.";
  }
}

function resetSurveillanceDrilldown() {
  state.surveillanceMapDrilldownCode = "";
  state.surveillanceMapDrilldownName = "";
  state.surveillanceMapParentCode = "";
  state.surveillanceMapParentName = "";
}

function stepBackSurveillanceDrilldown() {
  if (state.surveillanceMapParentCode) {
    state.surveillanceMapDrilldownCode = state.surveillanceMapParentCode;
    state.surveillanceMapDrilldownName = state.surveillanceMapParentName;
    state.surveillanceMapParentCode = "";
    state.surveillanceMapParentName = "";
    return;
  }
  resetSurveillanceDrilldown();
}

function stepBackDrilldown(target) {
  if (target.mapDrilldownParentCode) {
    target.mapDrilldownCode = target.mapDrilldownParentCode;
    target.mapDrilldownName = target.mapDrilldownParentName;
    target.mapDrilldownParentCode = "";
    target.mapDrilldownParentName = "";
    return;
  }
  target.mapDrilldownCode = "";
  target.mapDrilldownName = "";
}

function describeChinaDrilldownCaption(currentName, parentName) {
  if (parentName && currentName) {
    return `${parentName} / ${currentName} district-level heat distribution. Use trackpad or mouse wheel to zoom.`;
  }
  if (currentName) {
    return `${currentName} city-level or district-level heat distribution. Use trackpad or mouse wheel to zoom.`;
  }
  return "Province-level heat distribution. Use trackpad or mouse wheel to zoom.";
}

function buildSurveillanceTimelineSeries() {
  const preferredField = state.monitoringData.some((item) => item.collection_time)
    ? getMonitoringFieldDefinition("standard:collection_time")
    : getMonitoringFieldDefinition("standard:submit_time");
  return preferredField ? aggregateMonitoringTimeSeries(preferredField) : { labels: [], values: [] };
}

function getLatestMonitoringDate() {
  const values = state.monitoringData
    .flatMap((item) => [item.collection_time, item.submit_time])
    .filter(Boolean)
    .map((value) => new Date(value))
    .filter((date) => !Number.isNaN(date.getTime()))
    .sort((a, b) => b.getTime() - a.getTime());
  return values[0] ? values[0].toISOString() : "";
}

function renderSurveillanceRecentTable() {
  const rows = [...state.monitoringData]
    .sort((a, b) => {
      const aTime = new Date(a.collection_time || a.submit_time || 0).getTime();
      const bTime = new Date(b.collection_time || b.submit_time || 0).getTime();
      return bTime - aTime;
    })
    .slice(0, 8);
  if (!rows.length) {
    elements.surveillanceRecent.innerHTML = '<div class="empty-state">No recent records available.</div>';
    return;
  }
  elements.surveillanceRecent.innerHTML = `
    <table class="surveillance-table">
      <thead>
        <tr>
          <th>Genome</th>
          <th>Country</th>
          <th>Collection</th>
          <th>Type</th>
        </tr>
      </thead>
      <tbody>
        ${rows.map((item) => `
          <tr>
            <td>${escapeHtml(item.genome_id || "-")}</td>
            <td>${escapeHtml(item.country || "-")}</td>
            <td>${escapeHtml(formatDate(item.collection_time || item.submit_time))}</td>
            <td>${escapeHtml(item.sample_type || "-")}</td>
          </tr>
        `).join("")}
      </tbody>
    </table>
  `;
}

function renderSurveillanceMap() {
  if (!elements.surveillanceMap) return;
  elements.surveillanceMap.innerHTML = "";
  const field = getMonitoringFieldDefinition(state.surveillanceMapScope === "china" ? "standard:location" : "standard:country");
  if (!field) {
    elements.surveillanceMap.innerHTML = '<div class="empty-state">Map field is unavailable.</div>';
    return;
  }
  const mapContainer = document.createElement("div");
  mapContainer.className = "viz-geo-map";
  elements.surveillanceMap.appendChild(mapContainer);
  const widget = {
    id: "surveillance-screen",
    mapScope: state.surveillanceMapScope,
    mapDrilldownCode: state.surveillanceMapDrilldownCode,
    mapDrilldownName: state.surveillanceMapDrilldownName,
    mapDrilldownParentCode: state.surveillanceMapParentCode,
    mapDrilldownParentName: state.surveillanceMapParentName,
  };
  window.requestAnimationFrame(() => initializeMonitoringMap(mapContainer, field, widget, { allowDrilldown: true, mapKey: "surveillance-screen", onDrilldown: syncSurveillanceDrilldown }));
}

function syncSurveillanceDrilldown(widget) {
  state.surveillanceMapDrilldownCode = widget.mapDrilldownCode || "";
  state.surveillanceMapDrilldownName = widget.mapDrilldownName || "";
  state.surveillanceMapParentCode = widget.mapDrilldownParentCode || "";
  state.surveillanceMapParentName = widget.mapDrilldownParentName || "";
}

function addMonitoringWidgetFromBuilder() {
  if (state.monitoringWidgets.length >= MONITORING_MAX_WIDGETS) {
    showToast("You can keep at most 10 charts on the monitoring board.", true);
    return;
  }
  const field = getMonitoringFieldDefinition(elements.monitoringField.value);
  if (!field) {
    showToast("Please choose a field for the chart.", true);
    return;
  }
  const position = getNextMonitoringPosition();
  const widget = {
    id: createWidgetId(),
    title: elements.monitoringTitle.value.trim() || `${field.label} ${capitalize(elements.monitoringChartType.value)}`,
    chartType: elements.monitoringChartType.value,
    fieldKey: field.key,
    mapScope: elements.monitoringMapScope.value,
    collapsed: false,
    expanded: false,
    mapDrilldownCode: "",
    mapDrilldownName: "",
    mapDrilldownParentCode: "",
    mapDrilldownParentName: "",
    x: position.x,
    y: position.y,
    w: position.w,
    h: position.h,
  };
  state.monitoringWidgets.push(widget);
  saveMonitoringWidgetsToStorage();
  renderMonitoringBoard();
}

function createWidgetId() {
  return window.crypto?.randomUUID ? crypto.randomUUID() : `widget_${Date.now()}_${Math.random().toString(16).slice(2, 8)}`;
}

function syncMonitoringBuilderState() {
  const isMap = elements.monitoringChartType.value === "map";
  elements.monitoringMapScopeWrap.classList.toggle("hidden", !isMap);
}

function getNextMonitoringPosition() {
  return getNextMonitoringPositionForSize(MONITORING_DEFAULT_WIDGET.w, MONITORING_DEFAULT_WIDGET.h);
}

function rectsOverlap(a, b) {
  return !(a.x + a.w <= b.x || b.x + b.w <= a.x || a.y + a.h <= b.y || b.y + b.h <= a.y);
}

function renderMonitoringBoard() {
  destroyMonitoringLeafletMaps();
  elements.monitoringBoard.replaceChildren();
  elements.monitoringCount.textContent = `${state.monitoringWidgets.length} / ${MONITORING_MAX_WIDGETS} widgets`;
  elements.monitoringEmptyState.classList.toggle("hidden", state.monitoringWidgets.length > 0);
  if (!state.monitoringWidgets.length) {
    elements.monitoringBoard.style.height = "0px";
    return;
  }

  const boardWidth = elements.monitoringBoard.clientWidth || elements.monitoringBoard.offsetWidth || 1200;
  const colWidth = boardWidth / MONITORING_COLUMNS;
  let boardRows = 8;

  state.monitoringWidgets.forEach((widget) => {
    boardRows = Math.max(boardRows, widget.y + widget.h + 1);
    const node = elements.monitoringWidgetTemplate.content.firstElementChild.cloneNode(true);
    const titleText = node.querySelector(".monitoring-widget-title-text");
    const typeChip = node.querySelector(".monitoring-widget-type");
    const titleInput = node.querySelector(".monitoring-widget-title-input");
    const chartSelect = node.querySelector(".monitoring-widget-chart-select");
    const fieldSelect = node.querySelector(".monitoring-widget-field-select");
    const mapScopeWrap = node.querySelector(".monitoring-widget-map-scope-wrap");
    const mapScopeSelect = node.querySelector(".monitoring-widget-map-scope-select");
    const controlsWrap = node.querySelector(".monitoring-widget-controls");
    const mapToolbar = node.querySelector(".monitoring-widget-map-toolbar");
    const mapBackButton = node.querySelector(".monitoring-widget-map-back");
    const mapLevelLabel = node.querySelector(".monitoring-widget-map-level");
    const expandButton = node.querySelector(".monitoring-widget-expand");
    const toggleButton = node.querySelector(".monitoring-widget-toggle");
    const chartWrap = node.querySelector(".monitoring-widget-chart");
    node.dataset.widgetId = widget.id;
    if (!widget.expanded) {
      node.style.left = `${widget.x * colWidth}px`;
      node.style.top = `${widget.y * MONITORING_ROW_HEIGHT}px`;
      node.style.width = `${widget.w * colWidth}px`;
      node.style.height = `${widget.h * MONITORING_ROW_HEIGHT}px`;
    }

    titleText.textContent = widget.title || "Untitled Chart";
    typeChip.textContent = widget.chartType;
    titleInput.value = widget.title || "";
    chartSelect.value = widget.chartType || "bar";
    mapScopeSelect.value = widget.mapScope || "world";
    populateMonitoringWidgetFieldOptions(fieldSelect, widget.fieldKey);
    mapScopeWrap.classList.toggle("hidden", widget.chartType !== "map");
    controlsWrap.classList.toggle("hidden", Boolean(widget.collapsed));
    node.classList.toggle("is-collapsed", Boolean(widget.collapsed));
    node.classList.toggle("is-expanded", Boolean(widget.expanded));
    mapToolbar.classList.toggle("hidden", widget.chartType !== "map");
    mapLevelLabel.textContent = widget.mapDrilldownName || (widget.mapScope === "china" ? "China" : "World");
    mapBackButton.classList.toggle("hidden", !widget.mapDrilldownCode);
    toggleButton.textContent = widget.collapsed ? "Expand" : "Collapse";
    expandButton.textContent = widget.expanded ? "Restore" : "Enlarge";

    titleInput.addEventListener("input", () => {
      widget.title = titleInput.value.trim();
      titleText.textContent = widget.title || "Untitled Chart";
      saveMonitoringWidgetsToStorage();
    });
    chartSelect.addEventListener("change", () => {
      widget.chartType = chartSelect.value;
      typeChip.textContent = widget.chartType;
      mapScopeWrap.classList.toggle("hidden", widget.chartType !== "map");
      if (widget.chartType !== "map") {
        widget.mapDrilldownCode = "";
        widget.mapDrilldownName = "";
        widget.mapDrilldownParentCode = "";
        widget.mapDrilldownParentName = "";
      }
      saveMonitoringWidgetsToStorage();
      renderMonitoringBoard();
    });
    fieldSelect.addEventListener("change", () => {
      widget.fieldKey = fieldSelect.value;
      if (!widget.title) {
        const field = getMonitoringFieldDefinition(widget.fieldKey);
        titleText.textContent = field?.label || "Untitled Chart";
      }
      saveMonitoringWidgetsToStorage();
      renderMonitoringChart(chartWrap, widget);
    });
    mapScopeSelect.addEventListener("change", () => {
      widget.mapScope = mapScopeSelect.value;
      widget.mapDrilldownCode = "";
      widget.mapDrilldownName = "";
      widget.mapDrilldownParentCode = "";
      widget.mapDrilldownParentName = "";
      saveMonitoringWidgetsToStorage();
      renderMonitoringBoard();
    });
    expandButton.addEventListener("click", () => {
      const nextExpanded = !widget.expanded;
      state.monitoringWidgets.forEach((item) => {
        item.expanded = item.id === widget.id ? nextExpanded : false;
      });
      saveMonitoringWidgetsToStorage();
      renderMonitoringBoard();
    });
    toggleButton.addEventListener("click", () => {
      widget.collapsed = !widget.collapsed;
      saveMonitoringWidgetsToStorage();
      renderMonitoringBoard();
    });
    mapBackButton.addEventListener("click", () => {
      stepBackDrilldown(widget);
      saveMonitoringWidgetsToStorage();
      renderMonitoringBoard();
    });
    node.querySelector(".monitoring-widget-remove").addEventListener("click", () => {
      state.monitoringWidgets = state.monitoringWidgets.filter((item) => item.id !== widget.id);
      compactMonitoringWidgets();
      saveMonitoringWidgetsToStorage();
      renderMonitoringBoard();
    });
    node.querySelector(".monitoring-widget-drag-handle").addEventListener("pointerdown", (event) => startMonitoringDrag(event, widget.id, "move"));
    node.querySelector(".monitoring-widget-resize-handle").addEventListener("pointerdown", (event) => startMonitoringDrag(event, widget.id, "resize"));

    elements.monitoringBoard.appendChild(node);
    renderMonitoringChart(chartWrap, widget);
  });

  elements.monitoringBoard.style.height = `${boardRows * MONITORING_ROW_HEIGHT}px`;
}

function populateMonitoringWidgetFieldOptions(selectNode, selectedKey) {
  selectNode.innerHTML = getMonitoringFieldOptions()
    .map((item) => `<option value="${escapeHtml(item.key)}">${escapeHtml(item.label)}</option>`)
    .join("");
  selectNode.value = selectedKey || "";
}

function getMonitoringFieldDefinition(fieldKey) {
  return getMonitoringFieldOptions().find((item) => item.key === fieldKey) || null;
}

function renderMonitoringChart(container, widget) {
  const field = getMonitoringFieldDefinition(widget.fieldKey);
  if (!field) {
    container.innerHTML = '<div class="empty-state">Choose a field to render this widget.</div>';
    return;
  }
  const series = buildMonitoringSeries(field, widget.chartType);
  if (!series.labels.length) {
    container.innerHTML = '<div class="empty-state">No data available for this field.</div>';
    return;
  }
  if (widget.chartType === "pie") {
    container.innerHTML = renderPieChart(series, widget.title || field.label);
    return;
  }
  if (widget.chartType === "line") {
    container.innerHTML = renderLineChart(series, widget.title || field.label);
    return;
  }
  if (widget.chartType === "radar") {
    container.innerHTML = renderRadarChart(series, widget.title || field.label);
    return;
  }
  if (widget.chartType === "map") {
    container.innerHTML = renderMapChart(field, widget, widget.title || field.label);
    const mapContainer = container.querySelector(".viz-geo-map");
    window.requestAnimationFrame(() => initializeMonitoringMap(mapContainer, field, widget));
    return;
  }
  container.innerHTML = renderBarChart(series, widget.title || field.label);
}

function buildMonitoringSeries(field, chartType) {
  if (chartType === "line") {
    return aggregateMonitoringTimeSeries(field);
  }
  const counts = new Map();
  state.monitoringData.forEach((genome) => {
    getMonitoringFieldValues(genome, field).forEach((value) => {
      const key = String(value || "").trim();
      if (!key) return;
      counts.set(key, (counts.get(key) || 0) + 1);
    });
  });
  const entries = Array.from(counts.entries())
    .sort((a, b) => b[1] - a[1])
    .slice(0, chartType === "radar" ? 6 : 8);
  return {
    labels: entries.map(([label]) => label),
    values: entries.map(([, value]) => value),
  };
}

function aggregateMonitoringTimeSeries(field) {
  const counts = new Map();
  state.monitoringData.forEach((genome) => {
    getMonitoringFieldValues(genome, field).forEach((value) => {
      const bucket = toTimeBucket(value);
      if (!bucket) return;
      counts.set(bucket, (counts.get(bucket) || 0) + 1);
    });
  });
  const entries = Array.from(counts.entries()).sort((a, b) => a[0].localeCompare(b[0]));
  return { labels: entries.map(([label]) => label), values: entries.map(([, value]) => value) };
}

function toTimeBucket(value) {
  const date = new Date(value);
  if (Number.isNaN(date.getTime())) return "";
  const year = date.getFullYear();
  const month = `${date.getMonth() + 1}`.padStart(2, "0");
  return `${year}-${month}`;
}

function getMonitoringFieldValues(genome, field) {
  if (!genome || !field) return [];
  if (field.key.startsWith("standard:")) {
    const value = genome[field.key.replace("standard:", "")];
    return normalizeMonitoringValue(value, field.type);
  }
  const metaKey = field.key.replace("meta:", "");
  const item = (genome.custom_metadata || []).find((entry) => entry.key === metaKey);
  return normalizeMonitoringValue(item?.value, item?.type || field.type);
}

function normalizeMonitoringValue(value, type) {
  if (value == null) return [];
  if (type === "location" && typeof value === "object") {
    const province = value.province || "";
    const city = value.city || "";
    const district = value.district || "";
    const detail = value.detail || "";
    const joined = [province, city, district].filter(Boolean).join(" / ");
    return [joined || province || city || district || detail].filter(Boolean);
  }
  if (Array.isArray(value)) {
    return value.map((item) => String(item).trim()).filter(Boolean);
  }
  return [String(value).trim()].filter(Boolean);
}

function renderBarChart(series, title) {
  const max = Math.max(...series.values, 1);
  const bars = series.labels.map((label, index) => {
    const width = Math.max(10, (series.values[index] / max) * 100);
    return `
      <div class="viz-bar-row">
        <span class="viz-label">${escapeHtml(label)}</span>
        <div class="viz-bar-track"><span class="viz-bar-fill" style="width:${width}%"></span></div>
        <strong>${series.values[index]}</strong>
      </div>
    `;
  }).join("");
  return `<div class="viz-shell"><h4>${escapeHtml(title)}</h4><div class="viz-bar-chart">${bars}</div></div>`;
}

function renderPieChart(series, title) {
  const total = series.values.reduce((sum, value) => sum + value, 0) || 1;
  let current = 0;
  const segments = series.values.map((value, index) => {
    const start = current / total * Math.PI * 2;
    current += value;
    const end = current / total * Math.PI * 2;
    const path = describeArcSlice(70, 70, 58, start, end);
    return `<path d="${path}" fill="${chartColor(index)}"></path>`;
  }).join("");
  const legend = series.labels.map((label, index) => `
    <div class="viz-legend-item"><span class="viz-dot" style="background:${chartColor(index)}"></span>${escapeHtml(label)} <strong>${series.values[index]}</strong></div>
  `).join("");
  return `<div class="viz-shell"><h4>${escapeHtml(title)}</h4><div class="viz-pie-layout"><svg viewBox="0 0 140 140" class="viz-pie-chart">${segments}</svg><div class="viz-legend">${legend}</div></div></div>`;
}

function renderLineChart(series, title) {
  const max = Math.max(...series.values, 1);
  const points = series.values.map((value, index) => {
    const x = series.values.length === 1 ? 10 : (index / (series.values.length - 1)) * 100;
    const y = 90 - (value / max) * 70;
    return `${x},${y}`;
  }).join(" ");
  const markers = series.values.map((value, index) => {
    const x = series.values.length === 1 ? 10 : (index / (series.values.length - 1)) * 100;
    const y = 90 - (value / max) * 70;
    return `<circle cx="${x}" cy="${y}" r="2.8" fill="${chartColor(index)}"></circle>`;
  }).join("");
  const axisLabels = series.labels.map((label, index) => {
    const x = series.values.length === 1 ? 10 : (index / (series.values.length - 1)) * 100;
    return `<text x="${x}" y="98" text-anchor="middle">${escapeHtml(shortenLabel(label, 8))}</text>`;
  }).join("");
  return `<div class="viz-shell"><h4>${escapeHtml(title)}</h4><svg viewBox="0 0 100 100" class="viz-line-chart"><polyline fill="none" stroke="${chartColor(0)}" stroke-width="2.2" points="${points}"></polyline>${markers}<g class="viz-axis-labels">${axisLabels}</g></svg></div>`;
}

function renderRadarChart(series, title) {
  const count = Math.max(series.labels.length, 3);
  const max = Math.max(...series.values, 1);
  const points = series.values.map((value, index) => {
    const angle = (Math.PI * 2 * index) / count - Math.PI / 2;
    const radius = 18 + (value / max) * 26;
    const x = 50 + Math.cos(angle) * radius;
    const y = 50 + Math.sin(angle) * radius;
    return `${x},${y}`;
  }).join(" ");
  const spokes = Array.from({ length: count }, (_, index) => {
    const angle = (Math.PI * 2 * index) / count - Math.PI / 2;
    const x = 50 + Math.cos(angle) * 40;
    const y = 50 + Math.sin(angle) * 40;
    const label = shortenLabel(series.labels[index] || `Axis ${index + 1}`, 10);
    const lx = 50 + Math.cos(angle) * 46;
    const ly = 50 + Math.sin(angle) * 46;
    return `<line x1="50" y1="50" x2="${x}" y2="${y}"></line><text x="${lx}" y="${ly}" text-anchor="middle">${escapeHtml(label)}</text>`;
  }).join("");
  return `<div class="viz-shell"><h4>${escapeHtml(title)}</h4><svg viewBox="0 0 100 100" class="viz-radar-chart"><polygon points="50,12 86,36 72,84 28,84 14,36" class="viz-radar-grid"></polygon><g class="viz-radar-spokes">${spokes}</g><polygon points="${points}" class="viz-radar-area"></polygon></svg></div>`;
}

function renderMapChart(field, widget, title) {
  const mapId = `monitoring-map-${widget.id}`;
  const legend = renderMapLegend(field, widget);
  return `<div class="viz-shell"><h4>${escapeHtml(title)}</h4><div class="viz-map-layout"><div id="${escapeHtml(mapId)}" class="viz-geo-map"></div><div class="viz-legend">${legend}</div></div></div>`;
}

function extractProvinceFromValue(value) {
  const text = String(value || "").trim();
  if (!text) return "";
  return text.split("/")[0].trim().replace(/\s+/g, "");
}

function normalizeCountryForMap(value) {
  const text = String(value || "").trim();
  if (!text) return "";
  if (text === "中国" || text === "中华人民共和国") return "China";
  if (text === "英国") return "United Kingdom";
  if (text === "美国") return "United States";
  if (text === "日本") return "Japan";
  if (text === "韩国") return "South Korea";
  return text;
}

function renderMapLegend(field, widget) {
  const scope = widget.mapScope || "world";
  const counts = buildMapCounts(field, widget);
  const max = Math.max(...counts.values(), 0);
  if (!max) {
    return '<div class="field-note">No regional data available for this map.</div>';
  }
  const bins = buildHeatBins(max);
  return bins
    .map((bin) => `<div class="viz-legend-item"><span class="viz-dot" style="background:${heatColor(bin.sample, max, scope)}"></span>${escapeHtml(bin.label)}</div>`)
    .join("");
}

async function initializeMonitoringMap(container, field, widget, options = {}) {
  if (!window.L || !container) return;
  const scope = widget.mapScope || "world";
  const mapKey = options.mapKey || widget.id;
  const allowDrilldown = options.allowDrilldown ?? true;
  const geojson = scope === "china" ? await getChinaGeoJsonForWidget(widget) : WORLD_GEOJSON;
  if (!geojson) {
    container.innerHTML = '<div class="empty-state">Map boundary data is unavailable.</div>';
    return;
  }

  const counts = buildMapCounts(field, widget);
  const max = Math.max(...counts.values(), 1);
  destroyLeafletMapByKey(mapKey);
  const map = window.L.map(container, {
    zoomControl: true,
    attributionControl: false,
    dragging: true,
    scrollWheelZoom: true,
    touchZoom: true,
    doubleClickZoom: true,
    boxZoom: false,
    keyboard: false,
    tap: false,
    zoomSnap: 0.25,
  });
  monitoringLeafletMaps.set(mapKey, map);

  const layer = window.L.geoJSON(geojson, {
    style: (feature) => {
      const areaName = getFeatureAreaName(feature, scope);
      const value = counts.get(areaName) || 0;
      return {
        color: value ? "#4d6570" : "#b7c2c8",
        weight: value ? 1.2 : 0.8,
        fillColor: heatColor(value, max, scope),
        fillOpacity: value ? 0.92 : 0.12,
      };
    },
    onEachFeature: (feature, featureLayer) => {
      const areaName = getFeatureAreaName(feature, scope);
      const value = counts.get(areaName) || 0;
      featureLayer.bindTooltip(`${areaName || "Unknown"}: ${value} samples`, { sticky: true });
      if (allowDrilldown && scope === "china" && feature?.properties?.adcode) {
        featureLayer.on("mouseover", () => {
          featureLayer.setStyle({ weight: 1.4, color: "#284b63" });
        });
        featureLayer.on("mouseout", () => {
          featureLayer.setStyle({
            color: value ? "#6e7f89" : "#b7c2c8",
            weight: 1,
          });
        });
        featureLayer.on("click", () => {
          const currentLevel = getGeoJsonFeatureLevel(geojson);
          if (currentLevel === "district") return;
          if (!widget.mapDrilldownCode) {
            widget.mapDrilldownCode = String(feature.properties.adcode || "");
            widget.mapDrilldownName = areaName;
          } else if (!widget.mapDrilldownParentCode && currentLevel === "city") {
            widget.mapDrilldownParentCode = widget.mapDrilldownCode || "";
            widget.mapDrilldownParentName = widget.mapDrilldownName || "";
            widget.mapDrilldownCode = String(feature.properties.adcode || "");
            widget.mapDrilldownName = areaName;
          } else {
            return;
          }
          options.onDrilldown?.(widget);
          if (mapKey === "surveillance-screen") {
            renderSurveillanceScreen();
          } else {
            saveMonitoringWidgetsToStorage();
            renderMonitoringBoard();
          }
        });
      }
    },
  }).addTo(map);

  const bounds = layer.getBounds();
  if (bounds.isValid()) {
    map.fitBounds(bounds.pad(0.04));
  }
  window.setTimeout(() => map.invalidateSize(), 0);
}

function destroyMonitoringLeafletMaps() {
  Array.from(monitoringLeafletMaps.entries()).forEach(([key, map]) => {
    if (String(key) === "surveillance-screen") return;
    map.remove();
    monitoringLeafletMaps.delete(key);
  });
}

function destroyLeafletMapByKey(mapKey) {
  const existing = monitoringLeafletMaps.get(mapKey);
  if (!existing) return;
  existing.remove();
  monitoringLeafletMaps.delete(mapKey);
}

function buildMapCounts(field, widget) {
  const scope = widget.mapScope || "world";
  const counts = new Map();
  state.monitoringData.forEach((genome) => {
    getMonitoringFieldValues(genome, field).forEach((value) => {
      const areaName = scope === "china"
        ? extractChinaAreaForWidget(value, widget)
        : normalizeCountryForMap(value);
      if (!areaName) return;
      counts.set(areaName, (counts.get(areaName) || 0) + 1);
    });
  });
  return counts;
}

async function getChinaGeoJsonForWidget(widget) {
  if (!widget.mapDrilldownCode) {
    return CHINA_GEOJSON;
  }
  if (CHINA_DRILLDOWN_CACHE.has(widget.mapDrilldownCode)) {
    return CHINA_DRILLDOWN_CACHE.get(widget.mapDrilldownCode);
  }
  try {
    const localPaths = [
      `/static/data/china_drilldown/${encodeURIComponent(widget.mapDrilldownCode)}.geojson`,
      `/static/data/china_drilldown/cities/${encodeURIComponent(widget.mapDrilldownCode)}.geojson`,
    ];
    let geojson = null;
    for (const path of localPaths) {
      try {
        geojson = await requestJson(path);
        break;
      } catch {
        continue;
      }
    }
    if (!geojson) {
      geojson = await requestJson(`https://geo.datav.aliyun.com/areas_v3/bound/${encodeURIComponent(widget.mapDrilldownCode)}_full.json`);
    }
    CHINA_DRILLDOWN_CACHE.set(widget.mapDrilldownCode, geojson);
    return geojson;
  } catch (error) {
    console.error("Failed to load China drilldown boundary", widget.mapDrilldownCode, error);
    return null;
  }
}

function extractChinaAreaForWidget(value, widget) {
  const parts = extractChinaLocationParts(value);
  if (!parts.province) return "";
  if (!widget.mapDrilldownCode) {
    return parts.province;
  }
  if (widget.mapDrilldownParentCode) {
    if (parts.province !== widget.mapDrilldownParentName || parts.city !== widget.mapDrilldownName) {
      return "";
    }
    return parts.district || "";
  }
  if (parts.province !== widget.mapDrilldownName) {
    return "";
  }
  if (getChinaDrilldownMode(widget.mapDrilldownCode) === "district") {
    return parts.district || "";
  }
  return parts.city || "";
}

function getGeoJsonFeatureLevel(geojson) {
  const firstFeature = geojson?.features?.[0];
  return String(firstFeature?.properties?.level || "").trim().toLowerCase();
}

function extractChinaLocationParts(value) {
  const text = String(value || "").trim();
  if (!text) {
    return { province: "", city: "", district: "" };
  }
  const segments = text.split("/").map((item) => item.trim()).filter(Boolean);
  return {
    province: (segments[0] || "").replace(/\s+/g, ""),
    city: (segments[1] || "").replace(/\s+/g, ""),
    district: (segments[2] || "").replace(/\s+/g, ""),
  };
}

function getChinaDrilldownMode(code) {
  return ["110000", "120000", "310000", "500000"].includes(String(code || "")) ? "district" : "city";
}

function getFeatureAreaName(feature, scope) {
  const rawName = String(
    feature?.properties?.name ??
    feature?.properties?.NAME ??
    feature?.properties?.fullname ??
    feature?.properties?.行政区 ??
    ""
  ).trim();
  if (!rawName) return "";
  return scope === "china" ? rawName.replace(/\s+/g, "") : normalizeCountryForMap(rawName);
}

function heatColor(value, max, scope) {
  if (!value) return scope === "china" ? "#e2f0ef" : "#f7e6db";
  const ratio = Math.min(1, value / Math.max(max, 1));
  if (scope === "china") {
    if (ratio >= 0.8) return "#0b4f46";
    if (ratio >= 0.6) return "#0f7b6c";
    if (ratio >= 0.4) return "#29a08f";
    if (ratio >= 0.2) return "#6bcab9";
    return "#bdebe3";
  }
  if (ratio >= 0.8) return "#8f2f14";
  if (ratio >= 0.6) return "#b5481f";
  if (ratio >= 0.4) return "#d96f3d";
  if (ratio >= 0.2) return "#ef9b73";
  return "#f7d3bf";
}

function buildHeatBins(max) {
  if (max <= 1) {
    return [
      { sample: 1, label: "1 sample" },
    ];
  }
  const steps = [0.2, 0.4, 0.6, 0.8, 1];
  let previous = 1;
  return steps.map((ratio, index) => {
    const upper = Math.max(previous, Math.ceil(max * ratio));
    const label = index === 0 ? `1-${upper} samples` : `${previous + 1}-${upper} samples`;
    const sample = upper;
    previous = upper;
    return { sample, label };
  });
}

function describeArcSlice(cx, cy, r, startAngle, endAngle) {
  const start = polarToCartesian(cx, cy, r, endAngle);
  const end = polarToCartesian(cx, cy, r, startAngle);
  const largeArcFlag = endAngle - startAngle <= Math.PI ? "0" : "1";
  return `M ${cx} ${cy} L ${start.x} ${start.y} A ${r} ${r} 0 ${largeArcFlag} 0 ${end.x} ${end.y} Z`;
}

function polarToCartesian(cx, cy, r, angle) {
  return { x: cx + Math.cos(angle) * r, y: cy + Math.sin(angle) * r };
}

function chartColor(index) {
  const colors = ["#0f7b6c", "#d96f3d", "#2f5d8a", "#8f4c7d", "#4d7c0f", "#c05621", "#2563eb", "#be185d"];
  return colors[index % colors.length];
}

function shortenLabel(value, maxLength) {
  const text = String(value || "");
  return text.length > maxLength ? `${text.slice(0, Math.max(1, maxLength - 1))}…` : text;
}

function capitalize(value) {
  const text = String(value || "");
  return text ? `${text[0].toUpperCase()}${text.slice(1)}` : text;
}

function compactMonitoringWidgets() {
  const widgets = [...state.monitoringWidgets];
  state.monitoringWidgets = [];
  widgets.forEach((widget) => {
    const next = getNextMonitoringPositionForSize(widget.w, widget.h);
    state.monitoringWidgets.push({ ...widget, x: next.x, y: next.y });
  });
}

function getNextMonitoringPositionForSize(width, height) {
  const widgets = state.monitoringWidgets || [];
  let y = 0;
  while (y < 200) {
    for (let x = 0; x <= MONITORING_COLUMNS - width; x += 1) {
      if (!widgets.some((item) => rectsOverlap({ x, y, w: width, h: height }, item))) {
        return { x, y, w: width, h: height };
      }
    }
    y += 1;
  }
  return { x: 0, y: widgets.length * height, w: width, h: height };
}

function startMonitoringDrag(event, widgetId, mode) {
  event.preventDefault();
  const widget = state.monitoringWidgets.find((item) => item.id === widgetId);
  if (!widget) return;
  if (widget.expanded) return;
  const boardRect = elements.monitoringBoard.getBoundingClientRect();
  state.monitoringDrag = {
    mode,
    widgetId,
    startX: event.clientX,
    startY: event.clientY,
    originX: widget.x,
    originY: widget.y,
    originW: widget.w,
    originH: widget.h,
    boardWidth: boardRect.width || 1200,
  };
}

function onMonitoringPointerMove(event) {
  if (!state.monitoringDrag) return;
  const widget = state.monitoringWidgets.find((item) => item.id === state.monitoringDrag.widgetId);
  if (!widget) return;
  const colWidth = (state.monitoringDrag.boardWidth || 1200) / MONITORING_COLUMNS;
  const deltaCols = Math.round((event.clientX - state.monitoringDrag.startX) / colWidth);
  const deltaRows = Math.round((event.clientY - state.monitoringDrag.startY) / MONITORING_ROW_HEIGHT);
  if (state.monitoringDrag.mode === "move") {
    widget.x = clamp(state.monitoringDrag.originX + deltaCols, 0, MONITORING_COLUMNS - widget.w);
    widget.y = Math.max(0, state.monitoringDrag.originY + deltaRows);
  } else {
    widget.w = clamp(state.monitoringDrag.originW + deltaCols, 3, MONITORING_COLUMNS - widget.x);
    widget.h = Math.max(4, state.monitoringDrag.originH + deltaRows);
  }
  renderMonitoringBoard();
}

function onMonitoringPointerUp() {
  if (!state.monitoringDrag) return;
  state.monitoringDrag = null;
  saveMonitoringWidgetsToStorage();
}

function resetMonitoringLayout() {
  compactMonitoringWidgets();
  saveMonitoringWidgetsToStorage();
  renderMonitoringBoard();
}

function clamp(value, min, max) {
  return Math.min(Math.max(value, min), max);
}

function renderAuditList(items, total) {
  elements.auditList.replaceChildren();
  elements.auditCount.textContent = `${total} records`;
  elements.auditPageIndicator.textContent = `Page ${state.auditPage} / ${state.auditTotalPages}`;
  elements.prevAuditPage.disabled = state.auditPage <= 1;
  elements.nextAuditPage.disabled = state.auditPage >= state.auditTotalPages;
  if (!items.length) {
    const empty = document.createElement("div");
    empty.className = "detail-card empty-state";
    empty.textContent = "No audit records matched the current filters.";
    elements.auditList.appendChild(empty);
    return;
  }
  items.forEach((item) => {
    const node = elements.auditItemTemplate.content.firstElementChild.cloneNode(true);
    node.querySelector(".audit-operation").textContent = item.operation;
    const status = node.querySelector(".audit-status");
    status.textContent = item.status;
    status.classList.add(item.status.toLowerCase());
    node.querySelector(".audit-details").textContent = item.details || "No details provided.";
    node.querySelector(".audit-meta").textContent = `${item.genome_id || "global"} | ${item.operator || "system"} | ${formatDate(item.action_time)}`;
    elements.auditList.appendChild(node);
  });
}

function renderUserList(items) {
  elements.userList.replaceChildren();
  if (!items.length) {
    const empty = document.createElement("div");
    empty.className = "detail-card empty-state";
    empty.textContent = "No users found.";
    elements.userList.appendChild(empty);
    return;
  }
  items.forEach((item) => {
    const node = elements.userItemTemplate.content.firstElementChild.cloneNode(true);
    node.querySelector(".user-name").textContent = item.username;
    node.querySelector(".user-role").textContent = item.role;
    node.querySelector(".user-meta").textContent = `Created ${formatDate(item.created_time)} | Last login ${formatDate(item.last_login_time)}`;
    elements.userList.appendChild(node);
  });
}

async function selectGenome(recordId) {
  const genome = await requestJson(`/api/genomes/${encodeURIComponent(recordId)}`);
  state.selectedGenomeId = genome.id;
  await loadMetadataTemplates(genome.submitter);
  renderDetail(genome);
  await loadGenomes();
}

async function editGenome(recordId) {
  const genome = await requestJson(`/api/genomes/${encodeURIComponent(recordId)}`);
  state.selectedGenomeId = genome.id;
  state.mode = "edit";
  await loadMetadataTemplates(genome.submitter);
  fillForm(genome);
  renderDetail(genome);
  setActiveTab("genome-tab");
  openGenomeModal();
  await loadGenomes();
}

function openGenomeModal() {
  elements.genomeModal.classList.remove("hidden");
  elements.genomeModal.setAttribute("aria-hidden", "false");
}

function closeGenomeModal() {
  elements.genomeModal.classList.add("hidden");
  elements.genomeModal.setAttribute("aria-hidden", "true");
}

function openBulkImportModal() {
  elements.bulkImportModal?.classList.remove("hidden");
  elements.bulkImportModal?.setAttribute("aria-hidden", "false");
  if (elements.bulkImportFile) {
    elements.bulkImportFile.value = "";
  }
  if (elements.bulkImportStatus) {
    elements.bulkImportStatus.textContent = "No file selected.";
  }
}

function closeBulkImportModal() {
  elements.bulkImportModal?.classList.add("hidden");
  elements.bulkImportModal?.setAttribute("aria-hidden", "true");
}

function renderDetail(genome) {
  elements.detailStatus.textContent = "Loaded";
  elements.detailCard.classList.remove("empty-state");
  const orderedMetadata = mergeMetadataTemplates(state.metadataTemplates, genome.custom_metadata || []);
  const metadataMarkup = renderMetadataDetail(orderedMetadata);
  const locationParts = [genome.location?.province || "", genome.location?.city || "", genome.location?.district || ""].filter(Boolean);
  const locationText = locationParts.join(" / ") || "-";
  const locationDetail = genome.location?.detail ? `<br><span class="detail-path">${escapeHtml(genome.location.detail)}</span>` : "";
  elements.detailCard.innerHTML = `
    <h3>${genome.genome_id}</h3>
    <p>${genome.description || "No description provided."}</p>
    <dl class="detail-grid">
      <div><dt>Sample</dt><dd>${genome.sample_name}</dd></div>
      <div><dt>Species</dt><dd>${genome.species_name}</dd></div>
      <div><dt>TaxID</dt><dd>${genome.taxid}</dd></div>
      <div><dt>Length</dt><dd>${formatNumber(genome.genome_length)}</dd></div>
      <div><dt>Gender</dt><dd>${escapeHtml(genome.gender || "-")}</dd></div>
      <div><dt>Country</dt><dd>${escapeHtml(genome.country || "-")}</dd></div>
      <div><dt>Location</dt><dd>${escapeHtml(locationText)}${locationDetail}</dd></div>
      <div><dt>Collection Time</dt><dd>${escapeHtml(formatDate(genome.collection_time) || "-")}</dd></div>
      <div><dt>Sample Type</dt><dd>${escapeHtml(genome.sample_type || "-")}</dd></div>
      <div><dt>Sequencing Method</dt><dd>${escapeHtml(genome.sequencing_method || "-")}</dd></div>
      <div><dt>Submitter</dt><dd>${genome.submitter}</dd></div>
      <div><dt>Submitted</dt><dd>${formatDate(genome.submit_time)}</dd></div>
      <div><dt>Modified</dt><dd>${formatDate(genome.last_modified_time)}</dd></div>
      <div><dt>Path</dt><dd><code>${genome.genome_file_path}</code></dd></div>
    </dl>
    ${metadataMarkup}
  `;
}

function fillForm(genome) {
  document.getElementById("edit-genome-id").value = genome.id;
  document.getElementById("genome_id").value = genome.genome_id;
  document.getElementById("genome_id").disabled = true;
  document.getElementById("sample_name").value = genome.sample_name;
  document.getElementById("species_name").value = genome.species_name;
  document.getElementById("taxid").value = genome.taxid;
  document.getElementById("genome_file_path").value = genome.genome_file_path;
  document.getElementById("submitter").value = genome.submitter;
  document.getElementById("gender").value = genome.gender || "";
  populateCountryOptions(document.getElementById("country"), genome.country || "");
  document.getElementById("collection_time").value = toDatetimeLocalValue(genome.collection_time);
  document.getElementById("sample_type").value = genome.sample_type || "";
  document.getElementById("sequencing_method").value = genome.sequencing_method || "";
  syncBuiltinLocationField({ preserveValue: false });
  document.getElementById("builtin_location_province").value = genome.location?.province || "";
  populateBuiltinCityOptions(genome.location?.province || "", genome.location?.city || "");
  populateBuiltinDistrictOptions(genome.location?.province || "", genome.location?.city || "", genome.location?.district || "");
  document.getElementById("builtin_location_detail").value = genome.location?.detail || "";
  elements.uploadStatus.textContent = `Current FASTA: ${genome.genome_file_path}`;
  document.getElementById("description").value = genome.description || "";
  elements.formTitle.textContent = "Edit Genome";
  elements.submitButton.textContent = "Update Genome";
  elements.deleteButton.classList.remove("hidden");
  populateMetaFields(genome.custom_metadata || [], state.metadataTemplates);
  setActiveFormTab("basic-tab");
}

function resetForm(options = {}) {
  const preserveDetail = options.preserveDetail ?? false;
  state.mode = "create";
  elements.genomeForm.reset();
  document.getElementById("edit-genome-id").value = "";
  document.getElementById("genome_id").disabled = false;
  document.getElementById("submitter").value = state.currentUser.username;
  document.getElementById("gender").value = "";
  populateCountryOptions(document.getElementById("country"), "");
  document.getElementById("collection_time").value = "";
  document.getElementById("sample_type").value = "";
  document.getElementById("sequencing_method").value = "";
  syncBuiltinLocationField({ preserveValue: false });
  elements.formTitle.textContent = "Add Genome";
  elements.submitButton.textContent = "Add Genome";
  elements.deleteButton.classList.add("hidden");
  elements.fastaFileInput.value = "";
  elements.uploadStatus.textContent = "No file uploaded yet.";
  elements.uploadDropzone.classList.remove("is-active");
  populateMetaFields([], state.metadataTemplates);
  setActiveFormTab("basic-tab");
  if (!preserveDetail) {
    state.selectedGenomeId = null;
    elements.detailStatus.textContent = "Idle";
    elements.detailCard.className = "detail-card empty-state";
    elements.detailCard.textContent = "Select a genome card to inspect metadata and edit it.";
  }
  loadGenomes();
}

async function onSubmitGenomeForm(event) {
  event.preventDefault();
  const metadataFields = collectMetaFields();
  const metadataValidationError = validateMetaFieldsBeforeSubmit(metadataFields);
  if (metadataValidationError) {
    setActiveFormTab("meta-tab");
    showToast(metadataValidationError, true);
    return;
  }
  const isEdit = state.mode === "edit";
  const payload = {
    genome_id: document.getElementById("genome_id").value.trim(),
    sample_name: document.getElementById("sample_name").value.trim(),
    species_name: document.getElementById("species_name").value.trim(),
    taxid: Number(document.getElementById("taxid").value),
    genome_file_path: document.getElementById("genome_file_path").value.trim(),
    gender: document.getElementById("gender").value.trim(),
    country: document.getElementById("country").value.trim(),
    location: collectBuiltinLocation(),
    collection_time: document.getElementById("collection_time").value.trim(),
    sample_type: document.getElementById("sample_type").value.trim(),
    sequencing_method: document.getElementById("sequencing_method").value.trim(),
    description: document.getElementById("description").value.trim(),
    custom_metadata: metadataFields,
  };
  let response;
  if (isEdit) {
    response = await requestJson(`/api/genomes/${encodeURIComponent(document.getElementById("edit-genome-id").value)}`, {
      method: "PUT",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload),
    });
    showToast("Genome updated successfully.");
    closeGenomeModal();
  } else {
    response = await requestJson("/api/genomes", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload),
    });
    showToast("Genome added successfully.");
    closeGenomeModal();
    resetForm({ preserveDetail: true });
  }
  await Promise.all([loadGenomes(), loadAuditLogs()]);
  await loadMetadataTemplates(document.getElementById("submitter").value.trim() || state.currentUser?.username);
  if (response?.id) {
    await loadGenomeDetail(response.id);
  }
}

async function onDeleteGenome() {
  const recordId = document.getElementById("edit-genome-id").value;
  const genomeId = document.getElementById("genome_id").value;
  if (!recordId) return;
  const confirmed = window.confirm(`Delete genome ${genomeId}? This only removes the metadata record.`);
  if (!confirmed) return;
  await requestJson(`/api/genomes/${encodeURIComponent(recordId)}`, {
    method: "DELETE",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({}),
  });
  showToast(`Deleted ${genomeId}.`);
  resetForm();
  closeGenomeModal();
  await loadAuditLogs();
}

async function loadGenomeDetail(recordId) {
  const genome = await requestJson(`/api/genomes/${encodeURIComponent(recordId)}`);
  state.selectedGenomeId = genome.id;
  await loadMetadataTemplates(genome.submitter);
  renderDetail(genome);
  await loadGenomes();
}

async function onCreateUser(event) {
  event.preventDefault();
  const payload = {
    username: document.getElementById("new_username").value.trim(),
    password: document.getElementById("new_password").value,
    role: document.getElementById("new_role").value,
  };
  await requestJson("/api/users", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(payload),
  });
  elements.userForm.reset();
  showToast(`Created user ${payload.username}`);
  await loadUsers();
}

async function onUpdateProfile(event) {
  event.preventDefault();
  const payload = {
    username: document.getElementById("profile_username").value.trim(),
    display_name: document.getElementById("profile_display_name").value.trim(),
    email: document.getElementById("profile_email").value.trim(),
  };
  const updated = await requestJson("/api/profile", {
    method: "PUT",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(payload),
  });
  applyCurrentUser(updated);
  fillProfileForm(updated);
  if (state.mode === "create") {
    document.getElementById("submitter").value = updated.username;
  }
  if (state.currentUser.role !== "admin") {
    document.getElementById("search_submitter").value = updated.username;
  }
  showToast("Profile updated successfully.");
  await Promise.all([loadGenomes(), loadAuditLogs()]);
}

async function onChangePassword(event) {
  event.preventDefault();
  const payload = {
    current_password: document.getElementById("current_password").value,
    new_password: document.getElementById("new_profile_password").value,
    confirm_password: document.getElementById("confirm_profile_password").value,
  };
  await requestJson("/api/profile/password", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(payload),
  });
  elements.passwordForm.reset();
  showToast("Password changed successfully.");
}

async function onFileSelected(event) {
  const [file] = event.target.files;
  if (!file) return;
  await uploadFasta(file);
}

function onDragOver(event) {
  event.preventDefault();
  elements.uploadDropzone.classList.add("is-active");
}

function onDragLeave(event) {
  event.preventDefault();
  elements.uploadDropzone.classList.remove("is-active");
}

async function onDropFile(event) {
  event.preventDefault();
  elements.uploadDropzone.classList.remove("is-active");
  const [file] = event.dataTransfer.files;
  if (!file) return;
  await uploadFasta(file);
}

async function uploadFasta(file) {
  const formData = new FormData();
  formData.append("file", file);
  elements.uploadStatus.textContent = `Uploading ${file.name}...`;
  const response = await fetch("/api/upload-fasta", { method: "POST", body: formData });
  const data = await response.json();
  if (!response.ok) {
    elements.uploadStatus.textContent = "Upload failed.";
    showToast(data.error || "Upload failed", true);
    throw new Error(data.error || "Upload failed");
  }
  document.getElementById("genome_file_path").value = data.stored_path;
  elements.uploadStatus.textContent = `Uploaded ${data.filename} | validated length ${formatNumber(data.genome_length)} bp`;
  showToast(`Uploaded ${data.filename}`);
}

function createEmptyMetaField() {
  const suffix = window.crypto?.randomUUID ? crypto.randomUUID() : `${Date.now()}_${Math.random().toString(16).slice(2, 8)}`;
  return { key: `meta_${suffix}`, label: "", type: "text", value: "", options: [] };
}

function populateMetaFields(items, templates = []) {
  elements.metaFieldList.replaceChildren();
  const merged = mergeMetadataTemplates(templates, items);
  merged.forEach((item) => renderMetaField(item));
  updateMetaEmptyState();
}

function mergeMetadataTemplates(templates, values) {
  const valueMap = new Map(
    (Array.isArray(values) ? values : []).map((item) => [metadataTemplateKey(item), item])
  );
  const merged = [];
  (Array.isArray(templates) ? templates : []).forEach((template) => {
    const existing = valueMap.get(metadataTemplateKey(template));
    merged.push({
      key: existing?.key || template.key || createEmptyMetaField().key,
      label: existing?.label || template.label || "",
      type: existing?.type || template.type || "text",
      options: existing?.options || template.options || [],
      value: existing?.value || "",
      filename: existing?.filename || "",
    });
    valueMap.delete(metadataTemplateKey(template));
  });
  valueMap.forEach((item) => {
    merged.push({
      key: item.key || createEmptyMetaField().key,
      label: item.label || "",
      type: item.type || "text",
      options: item.options || [],
      value: item.value || "",
      filename: item.filename || "",
    });
  });
  return merged;
}

function renderMetaField(item) {
  const node = elements.metaFieldTemplate.content.firstElementChild.cloneNode(true);
  node.dataset.key = item.key || `meta_${Date.now()}`;
  const labelInput = node.querySelector(".meta-label-input");
  const typeSelect = node.querySelector(".meta-type-select");
  const optionsWrap = node.querySelector(".meta-options-wrap");
  const optionsInput = node.querySelector(".meta-options-input");
  const textWrap = node.querySelector(".meta-text-wrap");
  const valueInput = node.querySelector(".meta-value-input");
  const countryWrap = node.querySelector(".meta-country-wrap");
  const countryInput = node.querySelector(".meta-country-input");
  const datetimeWrap = node.querySelector(".meta-datetime-wrap");
  const datetimeInput = node.querySelector(".meta-datetime-input");
  const selectWrap = node.querySelector(".meta-select-wrap");
  const selectInput = node.querySelector(".meta-select-input");
  const locationWrap = node.querySelector(".meta-location-wrap");
  const locationProvinceInput = node.querySelector(".meta-location-province");
  const locationCityInput = node.querySelector(".meta-location-city");
  const locationDistrictInput = node.querySelector(".meta-location-district");
  const locationDetailInput = node.querySelector(".meta-location-detail");
  const locationNote = node.querySelector(".meta-location-note");
  const fileWrap = node.querySelector(".meta-file-wrap");
  const fileName = node.querySelector(".meta-file-name");
  const uploadButton = node.querySelector(".meta-file-upload-button");

  labelInput.value = item.label || "";
  typeSelect.value = item.type || "text";
  optionsInput.value = Array.isArray(item.options) ? item.options.join(", ") : "";
  valueInput.value = item.type === "text" ? item.value || "" : "";
  populateCountryOptions(countryInput, item.type === "country" ? item.value || "" : "");
  datetimeInput.value = item.type === "datetime" ? item.value || "" : "";
  locationProvinceInput.value = item.type === "location" ? item.value?.province || "" : "";
  locationCityInput.value = item.type === "location" ? item.value?.city || "" : "";
  locationDistrictInput.value = item.type === "location" ? item.value?.district || "" : "";
  locationDetailInput.value = item.type === "location" ? item.value?.detail || "" : "";
  node.dataset.fileValue = item.type === "file" ? item.value || "" : "";
  node.dataset.fileName = item.type === "file" ? item.filename || item.value || "" : "";
  fileName.textContent = node.dataset.fileName || "No file selected.";

  function refreshTypeState() {
    const currentType = typeSelect.value;
    optionsWrap.classList.toggle("hidden", currentType !== "select");
    textWrap.classList.toggle("hidden", currentType !== "text");
    countryWrap.classList.toggle("hidden", currentType !== "country");
    datetimeWrap.classList.toggle("hidden", currentType !== "datetime");
    selectWrap.classList.toggle("hidden", currentType !== "select");
    locationWrap.classList.toggle("hidden", currentType !== "location");
    fileWrap.classList.toggle("hidden", currentType !== "file");
    if (currentType === "select") {
      const options = optionsInput.value
        .split(",")
        .map((part) => part.trim())
        .filter(Boolean);
      selectInput.innerHTML = options.map((option) => `<option value="${escapeHtml(option)}">${escapeHtml(option)}</option>`).join("");
      if (!options.length) {
        selectInput.innerHTML = '<option value="">No options yet</option>';
      }
      if (item.type === "select" && item.value) {
        selectInput.value = item.value;
      }
    }
    if (currentType === "location") {
      syncLocationField(node, { preserveValue: true });
    }
  }

  optionsInput.addEventListener("input", refreshTypeState);
  typeSelect.addEventListener("change", refreshTypeState);
  countryInput.addEventListener("change", syncAllLocationFields);
  uploadButton.addEventListener("click", async () => {
    const picker = document.createElement("input");
    picker.type = "file";
    picker.addEventListener("change", async () => {
      const [file] = picker.files || [];
      if (!file) return;
      const uploaded = await uploadMetadataFile(file, fileName);
      node.dataset.fileValue = uploaded.stored_path;
      node.dataset.fileName = uploaded.filename;
      fileName.textContent = uploaded.filename;
    });
    picker.click();
  });
  node.querySelector(".meta-remove-button").addEventListener("click", () => {
    node.remove();
    updateMetaEmptyState();
  });
  node.querySelector(".meta-move-up-button").addEventListener("click", () => moveMetaField(node, -1));
  node.querySelector(".meta-move-down-button").addEventListener("click", () => moveMetaField(node, 1));

  refreshTypeState();
  locationNote.textContent = "Select country as China to enable province/city/district linkage.";
  locationProvinceInput.addEventListener("change", () => {
    populateCityOptions(node, locationProvinceInput.value, "");
    populateDistrictOptions(node, locationProvinceInput.value, locationCityInput.value, "");
  });
  locationCityInput.addEventListener("change", () => {
    populateDistrictOptions(node, locationProvinceInput.value, locationCityInput.value, "");
  });
  elements.metaFieldList.appendChild(node);
  syncAllLocationFields();
}

function populateCountryOptions(selectNode, selectedValue) {
  if (!selectNode) return;
  const values = COUNTRY_OPTIONS.map((item) => item.value);
  const optionMarkup = COUNTRY_OPTIONS.map(
    (item) => `<option value="${escapeHtml(item.value)}">${escapeHtml(item.label || item.value)}</option>`
  );
  if (selectedValue && !values.includes(selectedValue)) {
    optionMarkup.push(`<option value="${escapeHtml(selectedValue)}">${escapeHtml(selectedValue)}</option>`);
  }
  selectNode.innerHTML = ['<option value="">Select Country</option>'].concat(optionMarkup).join("");
  selectNode.value = selectedValue || "";
}

function syncBuiltinLocationField(options = {}) {
  const preserveValue = options.preserveValue ?? false;
  const countryValue = document.getElementById("country")?.value.trim() || "";
  const isChina = isChinaCountry(countryValue);
  const provinceSelect = document.getElementById("builtin_location_province");
  const citySelect = document.getElementById("builtin_location_city");
  const districtSelect = document.getElementById("builtin_location_district");
  const detailInput = document.getElementById("builtin_location_detail");
  const note = document.getElementById("builtin_location_note");
  if (!provinceSelect || !citySelect || !districtSelect || !detailInput || !note) return;

  provinceSelect.disabled = !isChina;
  citySelect.disabled = !isChina;
  districtSelect.disabled = !isChina;
  detailInput.disabled = !isChina;
  note.textContent = isChina
    ? "China selected: province/city/district linkage is enabled."
    : "Location linkage is only enabled when country is China. Leave location blank for other countries.";

  if (!isChina) {
    provinceSelect.innerHTML = '<option value="">Province / 省</option>';
    citySelect.innerHTML = '<option value="">City / 市</option>';
    districtSelect.innerHTML = '<option value="">District / 区县</option>';
    provinceSelect.value = "";
    citySelect.value = "";
    districtSelect.value = "";
    detailInput.value = "";
    return;
  }

  const provinceValue = preserveValue ? provinceSelect.value : "";
  const cityValue = preserveValue ? citySelect.value : "";
  const districtValue = preserveValue ? districtSelect.value : "";
  populateBuiltinProvinceOptions(provinceValue);
  populateBuiltinCityOptions(document.getElementById("builtin_location_province").value, cityValue);
  populateBuiltinDistrictOptions(
    document.getElementById("builtin_location_province").value,
    document.getElementById("builtin_location_city").value,
    districtValue
  );
}

function populateBuiltinProvinceOptions(selectedValue) {
  const provinceSelect = document.getElementById("builtin_location_province");
  const provinces = Object.keys(CHINA_REGION_TREE);
  provinceSelect.innerHTML = ['<option value="">Province / 省</option>']
    .concat(provinces.map((province) => `<option value="${escapeHtml(province)}">${escapeHtml(province)}</option>`))
    .join("");
  provinceSelect.value = provinces.includes(selectedValue) ? selectedValue : "";
}

function populateBuiltinCityOptions(province, selectedValue) {
  const citySelect = document.getElementById("builtin_location_city");
  const cities = province && CHINA_REGION_TREE[province] ? Object.keys(CHINA_REGION_TREE[province]) : [];
  citySelect.innerHTML = ['<option value="">City / 市</option>']
    .concat(cities.map((city) => `<option value="${escapeHtml(city)}">${escapeHtml(city)}</option>`))
    .join("");
  citySelect.value = cities.includes(selectedValue) ? selectedValue : "";
}

function populateBuiltinDistrictOptions(province, city, selectedValue) {
  const districtSelect = document.getElementById("builtin_location_district");
  const districts = province && city && CHINA_REGION_TREE[province]?.[city] ? CHINA_REGION_TREE[province][city] : [];
  districtSelect.innerHTML = ['<option value="">District / 区县</option>']
    .concat(districts.map((district) => `<option value="${escapeHtml(district)}">${escapeHtml(district)}</option>`))
    .join("");
  districtSelect.value = districts.includes(selectedValue) ? selectedValue : "";
}

function collectBuiltinLocation() {
  return {
    province: document.getElementById("builtin_location_province").value.trim(),
    city: document.getElementById("builtin_location_city").value.trim(),
    district: document.getElementById("builtin_location_district").value.trim(),
    detail: document.getElementById("builtin_location_detail").value.trim(),
  };
}

function collectMetaFields() {
  return Array.from(elements.metaFieldList.children).map((node) => {
    const type = node.querySelector(".meta-type-select").value;
    const base = {
      key: node.dataset.key,
      label: node.querySelector(".meta-label-input").value.trim(),
      type,
    };
    if (type === "text") {
      return { ...base, value: node.querySelector(".meta-value-input").value.trim() };
    }
    if (type === "country") {
      return { ...base, value: node.querySelector(".meta-country-input").value.trim() };
    }
    if (type === "datetime") {
      return { ...base, value: node.querySelector(".meta-datetime-input").value.trim() };
    }
    if (type === "select") {
      const options = node.querySelector(".meta-options-input").value.split(",").map((part) => part.trim()).filter(Boolean);
      return { ...base, options, value: node.querySelector(".meta-select-input").value.trim() };
    }
    if (type === "location") {
      return {
        ...base,
        value: {
          province: node.querySelector(".meta-location-province").value.trim(),
          city: node.querySelector(".meta-location-city").value.trim(),
          district: node.querySelector(".meta-location-district").value.trim(),
          detail: node.querySelector(".meta-location-detail").value.trim(),
        },
      };
    }
    return { ...base, value: node.dataset.fileValue || "", filename: node.dataset.fileName || "" };
  });
}

function validateMetaFieldsBeforeSubmit(items) {
  const seenLabels = new Set();
  for (const item of items) {
    const label = String(item?.label || "").trim();
    const normalizedLabel = label.replace(/\s+/g, " ").toLocaleLowerCase();
    if (!label) {
      return "Field Label cannot be empty.";
    }
    if (seenLabels.has(normalizedLabel)) {
      return `Field Label duplicated: ${label}`;
    }
    seenLabels.add(normalizedLabel);
  }
  return "";
}

function updateMetaEmptyState() {
  elements.metaEmptyState.classList.toggle("hidden", elements.metaFieldList.children.length > 0);
}

function moveMetaField(node, direction) {
  if (direction < 0) {
    const previous = node.previousElementSibling;
    if (previous) {
      elements.metaFieldList.insertBefore(node, previous);
    }
    return;
  }
  const next = node.nextElementSibling;
  if (next) {
    elements.metaFieldList.insertBefore(next, node);
  }
}

function syncAllLocationFields() {
  Array.from(elements.metaFieldList.children).forEach((node) => {
    if (node.querySelector(".meta-type-select")?.value === "location") {
      syncLocationField(node, { preserveValue: true });
    }
  });
}

function syncLocationField(node, options = {}) {
  const preserveValue = options.preserveValue ?? false;
  const countryValue = getSelectedCountryValue();
  const isChina = isChinaCountry(countryValue);
  const provinceSelect = node.querySelector(".meta-location-province");
  const citySelect = node.querySelector(".meta-location-city");
  const districtSelect = node.querySelector(".meta-location-district");
  const detailInput = node.querySelector(".meta-location-detail");
  const note = node.querySelector(".meta-location-note");

  provinceSelect.disabled = !isChina;
  citySelect.disabled = !isChina;
  districtSelect.disabled = !isChina;
  detailInput.disabled = !isChina;
  note.textContent = isChina
    ? "China selected: province/city/district linkage is enabled."
    : "Location linkage is only enabled when country is China. Leave this field blank for other countries.";

  if (!isChina) {
    provinceSelect.innerHTML = '<option value="">Province / 省</option>';
    citySelect.innerHTML = '<option value="">City / 市</option>';
    districtSelect.innerHTML = '<option value="">District / 区县</option>';
    detailInput.value = "";
    provinceSelect.value = "";
    citySelect.value = "";
    districtSelect.value = "";
    return;
  }

  const provinceValue = preserveValue ? provinceSelect.value : "";
  const cityValue = preserveValue ? citySelect.value : "";
  const districtValue = preserveValue ? districtSelect.value : "";
  populateProvinceOptions(node, provinceValue);
  populateCityOptions(node, provinceSelect.value, cityValue);
  populateDistrictOptions(node, provinceSelect.value, citySelect.value, districtValue);
}

function populateProvinceOptions(node, selectedValue) {
  const provinceSelect = node.querySelector(".meta-location-province");
  const provinces = Object.keys(CHINA_REGION_TREE);
  provinceSelect.innerHTML = ['<option value="">Province / 省</option>']
    .concat(provinces.map((province) => `<option value="${escapeHtml(province)}">${escapeHtml(province)}</option>`))
    .join("");
  provinceSelect.value = provinces.includes(selectedValue) ? selectedValue : "";
}

function populateCityOptions(node, province, selectedValue) {
  const citySelect = node.querySelector(".meta-location-city");
  const cities = province && CHINA_REGION_TREE[province] ? Object.keys(CHINA_REGION_TREE[province]) : [];
  citySelect.innerHTML = ['<option value="">City / 市</option>']
    .concat(cities.map((city) => `<option value="${escapeHtml(city)}">${escapeHtml(city)}</option>`))
    .join("");
  citySelect.value = cities.includes(selectedValue) ? selectedValue : "";
}

function populateDistrictOptions(node, province, city, selectedValue) {
  const districtSelect = node.querySelector(".meta-location-district");
  const districts = province && city && CHINA_REGION_TREE[province]?.[city] ? CHINA_REGION_TREE[province][city] : [];
  districtSelect.innerHTML = ['<option value="">District / 区县</option>']
    .concat(districts.map((district) => `<option value="${escapeHtml(district)}">${escapeHtml(district)}</option>`))
    .join("");
  districtSelect.value = districts.includes(selectedValue) ? selectedValue : "";
}

function getSelectedCountryValue() {
  const countryField = Array.from(elements.metaFieldList.children).find(
    (node) => node.querySelector(".meta-type-select")?.value === "country"
  );
  return countryField ? countryField.querySelector(".meta-country-input")?.value.trim() || "" : "";
}

function isChinaCountry(value) {
  const normalized = String(value || "").trim().toLowerCase();
  return ["中国", "china", "pr china", "people's republic of china", "中华人民共和国"].includes(normalized);
}

async function uploadMetadataFile(file, statusNode) {
  const formData = new FormData();
  formData.append("file", file);
  statusNode.textContent = `Uploading ${file.name}...`;
  const response = await fetch("/api/upload-metadata-file", { method: "POST", body: formData });
  const data = await response.json();
  if (!response.ok) {
    statusNode.textContent = "Upload failed.";
    showToast(data.error || "Metadata file upload failed", true);
    throw new Error(data.error || "Metadata file upload failed");
  }
  showToast(`Uploaded ${data.filename}`);
  return data;
}

function onBulkImportFileSelected() {
  const file = elements.bulkImportFile?.files?.[0];
  elements.bulkImportStatus.textContent = file
    ? `Selected ${file.name}`
    : "No file selected.";
}

async function onSubmitBulkImport() {
  const file = elements.bulkImportFile?.files?.[0];
  if (!file) {
    showToast("Please choose an import file first.", true);
    return;
  }
  const formData = new FormData();
  formData.append("file", file);
  elements.bulkImportStatus.textContent = `Importing ${file.name}...`;
  const response = await fetch("/api/bulk-import-genomes", {
    method: "POST",
    body: formData,
  });
  const data = await response.json();
  if (!response.ok) {
    elements.bulkImportStatus.textContent = "Import failed.";
    showToast(data.error || "Bulk import failed", true);
    throw new Error(data.error || "Bulk import failed");
  }
  elements.bulkImportStatus.textContent = `Finished: ${data.success_count} succeeded, ${data.failure_count} failed.`;
  renderBulkImportResult(data);
  await Promise.all([loadGenomes(), loadDatabaseTable(), loadAuditLogs()]);
  if (state.activeTab === "monitoring-tab") {
    await loadMonitoringData();
    renderSurveillanceScreen();
    renderMonitoringBoard();
  }
  showToast(`Bulk import completed: ${data.success_count} succeeded, ${data.failure_count} failed.`);
}

function renderBulkImportResult(data) {
  const items = Array.isArray(data?.results) ? data.results : [];
  if (!items.length) {
    elements.bulkImportResult.className = "detail-card empty-state";
    elements.bulkImportResult.textContent = "No row-level result returned.";
    return;
  }
  elements.bulkImportResult.className = "detail-card";
  elements.bulkImportResult.innerHTML = `
    <div class="results-meta">
      <span>Total ${escapeHtml(data.total || 0)} rows</span>
      <span>${escapeHtml(data.success_count || 0)} success / ${escapeHtml(data.failure_count || 0)} failed</span>
    </div>
    <div class="table-shell">
      <table class="data-table">
        <thead>
          <tr>
            <th>Row</th>
            <th>Genome ID</th>
            <th>Status</th>
            <th>Message</th>
          </tr>
        </thead>
        <tbody>
          ${items.map((item) => `
            <tr>
              <td>${escapeHtml(item.row_number || "-")}</td>
              <td>${escapeHtml(item.genome_id || "-")}</td>
              <td>${escapeHtml(item.status || "-")}</td>
              <td>${escapeHtml(item.error || "Imported successfully")}</td>
            </tr>
          `).join("")}
        </tbody>
      </table>
    </div>
  `;
}

function renderMetadataDetail(items) {
  if (!Array.isArray(items) || !items.length) return "";
  const entries = items
    .map((item) => {
      const label = escapeHtml(item.label || "Metadata");
      if (item.type === "location") {
        const parts = [
          item.value?.province || "",
          item.value?.city || "",
          item.value?.district || "",
        ].filter(Boolean);
        const locationLine = escapeHtml(parts.join(" / ") || "-");
        const detailLine = item.value?.detail ? `<br><span class="detail-path">${escapeHtml(item.value.detail)}</span>` : "";
        return `<div><dt>${label}</dt><dd>${locationLine}${detailLine}</dd></div>`;
      }
      if (item.type === "datetime") {
        const formatted = item.value ? formatDate(item.value) : "-";
        return `<div><dt>${label}</dt><dd>${escapeHtml(formatted)}</dd></div>`;
      }
      if (item.type === "file" && item.value) {
        const fileName = escapeHtml(item.filename || item.value.split("/").pop() || "Attached file");
        const filePath = escapeHtml(item.value);
        return `<div><dt>${label}</dt><dd><code>${fileName}</code><br><span class="detail-path">${filePath}</span></dd></div>`;
      }
      const value = escapeHtml(item.value || "-");
      return `<div><dt>${label}</dt><dd>${value}</dd></div>`;
    })
    .join("");
  return `<div class="detail-section"><h4>Custom Metadata</h4><dl class="detail-grid">${entries}</dl></div>`;
}

function onSearch(event) {
  event.preventDefault();
  state.page = 1;
  loadGenomes();
}

function toggleSearchPanel() {
  const isHidden = elements.searchForm.classList.contains("hidden");
  elements.searchForm.classList.toggle("hidden", !isHidden);
  elements.searchForm.setAttribute("aria-hidden", isHidden ? "false" : "true");
  elements.toggleSearchPanelButton.textContent = isHidden ? "Hide Search" : "Open Search";
  elements.toggleSearchPanelButton.setAttribute("aria-expanded", isHidden ? "true" : "false");
}

async function onSearchSubmitterChanged() {
  try {
    await loadSearchMetadataTemplates();
  } catch (error) {
    showToast(error.message || "Failed to load search metadata filters", true);
  }
}

function changePage(delta) {
  const next = state.page + delta;
  if (next < 1 || next > state.totalPages) return;
  state.page = next;
  loadGenomes();
}

function changeDatabasePage(delta) {
  const next = state.databasePage + delta;
  if (next < 1 || next > state.databaseTotalPages) return;
  state.databasePage = next;
  loadDatabaseTable();
}

function getSearchTemplateSubmitter() {
  const submitter = document.getElementById("search_submitter")?.value.trim();
  return submitter || state.currentUser?.username || "";
}

function rerenderSearchFilters() {
  const currentFilters = collectSearchFilters();
  elements.searchFilterList.replaceChildren();
  currentFilters.forEach((item) => renderSearchFilter(item));
  updateSearchFilterEmptyState();
}

function renderSearchFilter(item = null) {
  const node = elements.searchFilterTemplate.content.firstElementChild.cloneNode(true);
  const fieldSelect = node.querySelector(".search-filter-field");
  const operatorSelect = node.querySelector(".search-filter-operator");
  const valueWrap = node.querySelector(".search-filter-value-wrap");
  const valueInput = node.querySelector(".search-filter-value");
  populateSearchFieldOptions(fieldSelect, item?.key || "");

  function refreshFieldState() {
    const template = getSearchTemplateByKey(fieldSelect.value);
    populateSearchOperatorOptions(operatorSelect, template?.type, item?.operator || "");
    if (template?.type === "country") {
      const selectedValue = item?.value && typeof item.value === "string" ? item.value : "";
      valueWrap.innerHTML = `
        <span>Value</span>
        <select class="search-filter-value"></select>
      `;
      populateCountryOptions(valueWrap.querySelector(".search-filter-value"), selectedValue);
    } else {
      valueWrap.innerHTML = `
        <span>Value</span>
        <input class="search-filter-value" placeholder="Enter search value">
      `;
      valueWrap.querySelector(".search-filter-value").value = item?.value || "";
    }
    syncSearchValueVisibility(node);
  }

  fieldSelect.addEventListener("change", () => {
    item = { key: fieldSelect.value, operator: "", value: "" };
    refreshFieldState();
  });
  operatorSelect.addEventListener("change", () => syncSearchValueVisibility(node));
  node.querySelector(".search-filter-remove-button").addEventListener("click", () => {
    node.remove();
    updateSearchFilterEmptyState();
  });

  refreshFieldState();
  elements.searchFilterList.appendChild(node);
}

function populateSearchFieldOptions(selectNode, selectedKey) {
  const options = ['<option value="">Select Field</option>']
    .concat(
      state.searchMetadataTemplates.map(
        (item) => {
          const key = item.key && String(item.key).startsWith("standard:") ? item.key : `meta:${item.key || ""}`;
          return `<option value="${escapeHtml(key)}">${escapeHtml(item.label || item.key || "Metadata")}</option>`;
        }
      )
    )
    .join("");
  selectNode.innerHTML = options;
  selectNode.value = selectedKey || "";
}

function populateSearchOperatorOptions(selectNode, fieldType, selectedOperator) {
  const operators = getSearchOperators(fieldType);
  selectNode.innerHTML = operators
    .map((item) => `<option value="${escapeHtml(item.value)}">${escapeHtml(item.label)}</option>`)
    .join("");
  selectNode.value = operators.some((item) => item.value === selectedOperator) ? selectedOperator : operators[0]?.value || "contains";
}

function getSearchOperators(fieldType) {
  if (fieldType === "select" || fieldType === "country") {
    return [
      { value: "equals", label: "is" },
      { value: "not_equals", label: "is not" },
      { value: "empty", label: "is empty" },
      { value: "not_empty", label: "is not empty" },
    ];
  }
  if (fieldType === "location") {
    return [
      { value: "contains", label: "contains" },
      { value: "not_contains", label: "does not contain" },
      { value: "equals", label: "is" },
      { value: "not_equals", label: "is not" },
      { value: "empty", label: "is empty" },
      { value: "not_empty", label: "is not empty" },
    ];
  }
  return [
    { value: "contains", label: "contains" },
    { value: "not_contains", label: "does not contain" },
    { value: "equals", label: "is" },
    { value: "not_equals", label: "is not" },
    { value: "empty", label: "is empty" },
    { value: "not_empty", label: "is not empty" },
  ];
}

function syncSearchValueVisibility(node) {
  const operator = node.querySelector(".search-filter-operator")?.value || "contains";
  const valueWrap = node.querySelector(".search-filter-value-wrap");
  const valueInput = valueWrap.querySelector(".search-filter-value");
  const hideValue = operator === "empty" || operator === "not_empty";
  valueInput.disabled = hideValue;
  valueWrap.classList.toggle("hidden", hideValue);
  if (hideValue) {
    valueInput.value = "";
  }
}

function collectSearchFilters() {
  return Array.from(elements.searchFilterList.children)
    .map((node) => {
      const key = node.querySelector(".search-filter-field")?.value || "";
      const operator = node.querySelector(".search-filter-operator")?.value || "";
      const value = node.querySelector(".search-filter-value")?.value?.trim?.() || "";
      const template = getSearchTemplateByKey(key);
      if (!key || !operator || !template) return null;
      return {
        key,
        label: template.label || "",
        type: template.type || "text",
        operator,
        value,
      };
    })
    .filter(Boolean);
}

function getSearchTemplateByKey(key) {
  return state.searchMetadataTemplates.find((item) => {
    const templateKey = item.key && String(item.key).startsWith("standard:") ? item.key : `meta:${item.key || ""}`;
    return templateKey === key;
  }) || null;
}

function updateSearchFilterEmptyState() {
  elements.searchFilterEmptyState.classList.toggle("hidden", elements.searchFilterList.children.length > 0);
}

function onAuditSearch(event) {
  event.preventDefault();
  state.auditPage = 1;
  loadAuditLogs();
}

function changeAuditPage(delta) {
  const next = state.auditPage + delta;
  if (next < 1 || next > state.auditTotalPages) return;
  state.auditPage = next;
  loadAuditLogs();
}

async function requestJson(url, options = {}) {
  const response = await fetch(url, options);
  const data = await response.json();
  if (!response.ok) {
    showToast(data.error || "Request failed", true);
    throw new Error(data.error || "Request failed");
  }
  return data;
}

function formatDate(value) {
  if (!value) return "-";
  return new Date(value).toLocaleString();
}

function toDatetimeLocalValue(value) {
  if (!value) return "";
  const date = new Date(value);
  if (Number.isNaN(date.getTime())) return "";
  const pad = (part) => String(part).padStart(2, "0");
  return `${date.getFullYear()}-${pad(date.getMonth() + 1)}-${pad(date.getDate())}T${pad(date.getHours())}:${pad(date.getMinutes())}`;
}

function formatNumber(value) {
  return new Intl.NumberFormat().format(value);
}

function showToast(message, isError = false) {
  elements.toast.textContent = message;
  elements.toast.style.background = isError ? "rgba(166, 54, 42, 0.95)" : "rgba(16, 32, 39, 0.92)";
  elements.toast.classList.remove("hidden");
  window.clearTimeout(showToast.timer);
  showToast.timer = window.setTimeout(() => elements.toast.classList.add("hidden"), 2800);
}

function escapeHtml(value) {
  return String(value)
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;")
    .replaceAll("'", "&#39;");
}

function metadataTemplateKey(item) {
  return `${String(item?.label || "").trim().toLowerCase()}::${String(item?.type || "").trim().toLowerCase()}`;
}
