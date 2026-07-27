from __future__ import annotations

import json
from pathlib import Path
import re


ROOT = Path(__file__).parents[1]
STATIC = ROOT / "bac_analysis_portal" / "static"
TEMPLATES = ROOT / "bac_analysis_portal" / "templates"


def read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def test_frontend_growth_budgets() -> None:
    budgets = json.loads(read(ROOT / "frontend-size-budgets.json"))
    for relative_path, maximum_lines in budgets.items():
        source = read(ROOT / relative_path)
        assert len(source.splitlines()) <= maximum_lines, relative_path


def test_workbench_uses_module_entrypoint_and_domain_partials() -> None:
    index = read(TEMPLATES / "index.html")
    assert '<script type="module" src="{{ url_for(\'static\', filename=\'app.js\') }}' in index
    assert "20260612-workbench-modules-03" in index
    for partial in (
        "workbench/tabs/project.html",
        "workbench/tabs/admin.html",
        "workbench/tabs/audit.html",
        "workbench/modals/project.html",
        "workbench/modals/admin-users.html",
    ):
        assert f'{{% include "{partial}" %}}' in index


def test_domain_modules_do_not_import_each_other() -> None:
    domain_root = STATIC / "workbench" / "domains"
    offenders = []
    for module_path in domain_root.rglob("*.js"):
        source = read(module_path)
        if re.search(r"""from\s+["'][^"']*/domains/""", source):
            offenders.append(str(module_path.relative_to(ROOT)))
    assert offenders == []


def test_migrated_domain_event_bindings_stay_out_of_entrypoint() -> None:
    entrypoint = read(STATIC / "app.js")
    for binding in (
        "elements.settingsForm?.addEventListener",
        "elements.userForm?.addEventListener",
        "elements.projectCreateForm?.addEventListener",
        "elements.refreshAuditButton?.addEventListener",
        "function loadAuditLogs",
        "function renderAuditLogs",
        "function handleAuditTableClick",
        "function loadAdminSettings",
        "function loadUsers",
        "function setActiveAdminSection",
        'state.pathBrowser.selector === "workspace_root"',
        "PROJECT_MILESTONE_OVERRIDE_STORAGE_KEY",
        "function renderProjectManagement",
        "function buildProjectManagementProjects",
        "function openProjectMilestoneModal",
    ):
        assert binding not in entrypoint
    assert 'navigation.register("project-tab", projectDomain)' in entrypoint
    assert 'navigation.register("admin-tab", adminDomain)' in entrypoint
    assert 'navigation.register("audit-tab", auditDomain)' in entrypoint
    assert "await navigation.activate(tabId)" in entrypoint


def test_path_browser_dispatches_registered_domain_targets() -> None:
    entrypoint = read(STATIC / "app.js")
    adapter = read(STATIC / "workbench" / "core" / "path-browser.js")
    admin = read(STATIC / "workbench" / "domains" / "admin" / "index.js")
    assert "registerTarget(selector, handler)" in adapter
    assert "await pathBrowser.handleSelection" in entrypoint
    assert 'runtime.pathBrowser.registerTarget("workspace_root"' in admin
    assert 'runtime.pathBrowser.registerTarget("conda_root"' in admin


def test_project_domain_owns_management_implementation() -> None:
    management = read(STATIC / "workbench" / "domains" / "project" / "management.js")
    project_entry = read(STATIC / "workbench" / "domains" / "project" / "index.js")
    assert "function buildProjectManagementProjects" in management
    assert "function renderProjectManagement" in management
    assert "function openProjectMilestoneModal" in management
    assert "createProjectManagement(runtime)" in project_entry
