const KEYS = {
  milestones: "bac-project-milestone-overrides",
  owners: "bac-project-owner-overrides",
  definitions: "bac-project-definitions",
  hidden: "bac-project-hidden",
  archived: "bac-project-archived",
};

function readObject(key) {
  try {
    const parsed = JSON.parse(window.localStorage.getItem(key) || "{}");
    return parsed && typeof parsed === "object" && !Array.isArray(parsed) ? parsed : {};
  } catch (_error) {
    return {};
  }
}

function readList(key, objects = false) {
  try {
    const parsed = JSON.parse(window.localStorage.getItem(key) || "[]");
    if (!Array.isArray(parsed)) return [];
    return objects
      ? parsed.filter((item) => item && typeof item === "object")
      : parsed.map((item) => String(item || "").trim()).filter(Boolean);
  } catch (_error) {
    return [];
  }
}

function write(key, value) {
  window.localStorage.setItem(key, JSON.stringify(value));
}

export const readProjectMilestoneOverrides = () => readObject(KEYS.milestones);
export const writeProjectMilestoneOverrides = (value) => write(KEYS.milestones, value || {});
export const readProjectOwnerOverrides = () => readObject(KEYS.owners);
export const writeProjectOwnerOverrides = (value) => write(KEYS.owners, value || {});
export const readProjectDefinitions = () => readList(KEYS.definitions, true);
export const writeProjectDefinitions = (value) => write(KEYS.definitions, Array.isArray(value) ? value : []);
export const readHiddenProjectKeys = () => readList(KEYS.hidden);
export const writeHiddenProjectKeys = (value) => write(KEYS.hidden, Array.isArray(value) ? value : []);
export const readArchivedProjectKeys = () => readList(KEYS.archived);
export const writeArchivedProjectKeys = (value) => write(KEYS.archived, Array.isArray(value) ? value : []);

export function clearProjectLocalState(projectKey = "") {
  const key = String(projectKey || "").trim();
  if (!key) return;
  const milestones = readProjectMilestoneOverrides();
  const owners = readProjectOwnerOverrides();
  if (Object.prototype.hasOwnProperty.call(milestones, key)) {
    delete milestones[key];
    writeProjectMilestoneOverrides(milestones);
  }
  if (Object.prototype.hasOwnProperty.call(owners, key)) {
    delete owners[key];
    writeProjectOwnerOverrides(owners);
  }
}

export function removeProjectDefinition(projectKey = "") {
  const key = String(projectKey || "").trim();
  if (!key) return;
  writeProjectDefinitions(readProjectDefinitions().filter((item) => String(item?.key || "").trim() !== key));
  const archived = new Set(readArchivedProjectKeys());
  if (archived.delete(key)) writeArchivedProjectKeys(Array.from(archived));
}

export function hideAutoProject(projectKey = "") {
  const key = String(projectKey || "").trim();
  if (!key) return;
  const hidden = new Set(readHiddenProjectKeys());
  hidden.add(key);
  writeHiddenProjectKeys(Array.from(hidden));
}
