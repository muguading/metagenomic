(function () {
  function onReady(fn) {
    if (document.readyState === "loading") {
      document.addEventListener("DOMContentLoaded", fn, { once: true });
    } else {
      fn();
    }
  }

  function hideComplexActions(root) {
    root.querySelectorAll('[data-task-action="rerun"], [data-task-action="rebuild"], [data-task-action="delete"]').forEach((node) => {
      node.style.display = "none";
    });
    root.querySelectorAll('[data-database-delete], [data-database-submit-main], .database-action-dropdown-danger').forEach((node) => {
      node.style.display = "none";
    });
  }

  function forceDemoView() {
    document.body.classList.add("demo-lite");

    const allowedTabs = new Set(["queue-tab", "database-tab", "admin-tab"]);
    document.querySelectorAll(".app-nav-button[data-tab-target]").forEach((button) => {
      const target = String(button.getAttribute("data-tab-target") || "");
      if (!allowedTabs.has(target)) {
        button.style.display = "none";
      }
    });

    const visiblePanel = document.querySelector(".tab-panel.active");
    if (!visiblePanel || !allowedTabs.has(visiblePanel.id)) {
      const databaseButton = document.querySelector('.app-nav-button[data-tab-target="database-tab"]');
      if (databaseButton instanceof HTMLElement) {
        databaseButton.click();
      }
    }

    const reportPreset = document.getElementById("database-report-preset");
    if (reportPreset instanceof HTMLSelectElement) {
      reportPreset.value = "general";
      reportPreset.dispatchEvent(new Event("change", { bubbles: true }));
    }

    const reportSectionButton = document.querySelector('[data-database-section="database-report-panel"]');
    if (reportSectionButton instanceof HTMLElement) {
      reportSectionButton.addEventListener("dblclick", (event) => event.preventDefault());
    }

    hideComplexActions(document);
  }

  onReady(() => {
    forceDemoView();
    const observer = new MutationObserver(() => hideComplexActions(document));
    observer.observe(document.body, { childList: true, subtree: true });
  });
})();
