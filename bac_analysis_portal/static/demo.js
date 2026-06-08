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

    const allowedTabs = new Set(["submission-tab", "queue-tab", "server-tab", "database-tab", "admin-tab"]);
    document.querySelectorAll(".app-nav-button[data-tab-target]").forEach((button) => {
      const target = String(button.getAttribute("data-tab-target") || "");
      if (!allowedTabs.has(target)) {
        button.style.display = "none";
      }
    });

    const visiblePanel = document.querySelector(".tab-panel.active");
    if (!visiblePanel || !allowedTabs.has(visiblePanel.id)) {
      const queueButton = document.querySelector('.app-nav-button[data-tab-target="queue-tab"]');
      if (queueButton instanceof HTMLElement) {
        queueButton.click();
      }
    }

    hideComplexActions(document);
  }

  onReady(() => {
    forceDemoView();
    const observer = new MutationObserver(() => hideComplexActions(document));
    observer.observe(document.body, { childList: true, subtree: true });
  });
})();
