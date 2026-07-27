from pathlib import Path


ROOT = Path(__file__).parents[1]
APP_JS = (ROOT / "bac_analysis_portal" / "static" / "app.js").read_text(encoding="utf-8")


def test_database_navigation_uses_cache_and_defers_hidden_views() -> None:
    assert "databaseLoadedAt: 0" in APP_JS
    assert "databaseLoadPromise: null" in APP_JS
    assert "databaseCacheAge > 60000" in APP_JS
    assert 'state.databaseSection === "database-monitor-panel"' in APP_JS
    assert 'state.databaseSection === "database-report-panel"' in APP_JS
    assert 'state.activeTab === "host-tab"' in APP_JS
    assert "const pageRows = sortedRows.slice(pageStart, pageStart + pageSize);" in APP_JS
    assert "${pageRows.length ? pageRows.map((row) => `" in APP_JS


def test_monitoring_does_not_unconditionally_render_report() -> None:
    monitoring = APP_JS.split("async function renderDatabaseMonitoring()", 1)[1].split("function renderDatabaseReport()", 1)[0]
    assert 'if (state.databaseSection === "database-report-panel")' in monitoring
    assert "\n  renderDatabaseReport();\n" not in monitoring
