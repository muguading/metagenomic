import { createProjectManagement } from "./management.js";

let manager = null;
let initialized = false;

export function init(runtime) {
  manager = createProjectManagement(runtime);
  if (initialized) return;
  initialized = true;
  const elements = runtime.selectors.elements();
  elements.closeProjectMilestoneModalButton?.addEventListener("click", manager.closeProjectMilestoneModal);
  elements.cancelProjectMilestoneButton?.addEventListener("click", manager.closeProjectMilestoneModal);
  elements.projectMilestoneBackdrop?.addEventListener("click", manager.closeProjectMilestoneModal);
  elements.projectMilestoneDeleteButton?.addEventListener("click", manager.onDeleteProjectMilestone);
  elements.projectMilestoneForm?.addEventListener("submit", manager.onSubmitProjectMilestone);
  elements.closeProjectOwnerModalButton?.addEventListener("click", manager.closeProjectOwnerModal);
  elements.cancelProjectOwnerButton?.addEventListener("click", manager.closeProjectOwnerModal);
  elements.projectOwnerBackdrop?.addEventListener("click", manager.closeProjectOwnerModal);
  elements.projectOwnerForm?.addEventListener("submit", manager.onSubmitProjectOwner);
  elements.closeProjectDeleteModalButton?.addEventListener("click", manager.closeProjectDeleteModal);
  elements.cancelProjectDeleteButton?.addEventListener("click", manager.closeProjectDeleteModal);
  elements.projectDeleteBackdrop?.addEventListener("click", manager.closeProjectDeleteModal);
  elements.projectDeleteForm?.addEventListener("submit", manager.onSubmitProjectDelete);
  elements.openProjectCreateModalButton?.addEventListener("click", manager.handleOpenProjectCreateModal);
  elements.closeProjectCreateModalButton?.addEventListener("click", manager.closeProjectCreateModal);
  elements.cancelProjectCreateButton?.addEventListener("click", manager.closeProjectCreateModal);
  elements.projectCreateBackdrop?.addEventListener("click", manager.closeProjectCreateModal);
  elements.projectCreateForm?.addEventListener("submit", manager.onSubmitProjectCreate);
  elements.addProjectMilestoneChildButton?.addEventListener("click", manager.appendProjectMilestoneChildRow);
  elements.projectMilestoneChildren?.addEventListener("click", (event) => {
    const target = event.target;
    if (!(target instanceof HTMLElement)) return;
    target.closest("[data-project-milestone-child-remove]")?.closest(".project-milestone-child-row")?.remove();
  });
}

export function activate() {
  manager?.renderProjectManagement();
}

export function refresh() {
  manager?.renderProjectManagement();
}
