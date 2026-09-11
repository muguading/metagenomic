from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).parents[1]
REPORT_CSS = (ROOT / "bac_analysis_portal" / "static" / "report.css").read_text(encoding="utf-8")
REPORT_RUNTIME = (ROOT / "bac_analysis_portal" / "static" / "report_runtime.js").read_text(encoding="utf-8")


def test_report_metric_cells_use_compact_neutral_surface() -> None:
    assert "--report-radius-lg: 18px;" in REPORT_CSS
    assert "--report-radius-md: 14px;" in REPORT_CSS
    assert ".mini-stat-grid" in REPORT_CSS
    assert "gap: 8px;" in REPORT_CSS
    assert "border-radius: 10px;" in REPORT_CSS
    assert "background: color-mix(in srgb, var(--report-paper) 94%, var(--report-soft) 6%);" in REPORT_CSS


def test_report_tab_controls_expose_interaction_and_reduced_motion_states() -> None:
    assert ".report-tab-button,\n.subreport-tab-button" in REPORT_CSS
    assert "min-height: 44px;" in REPORT_CSS
    assert ".report-tab-button:hover:not(.active)" in REPORT_CSS
    assert ".report-tab-button:active" in REPORT_CSS
    assert ".report-tab-button:focus-visible" in REPORT_CSS
    assert "@media (prefers-reduced-motion: reduce)" in REPORT_CSS


def test_variant_quality_switches_have_spacing_and_accessible_state() -> None:
    assert REPORT_RUNTIME.count('class="mini-stat-grid variant-summary-grid"') == 2
    assert REPORT_RUNTIME.count('role="group" aria-label=') >= 2
    assert 'data-nextclade-variant-tab="high" aria-pressed="true"' in REPORT_RUNTIME
    assert 'data-rsv-variant-tab="high" aria-pressed="true"' in REPORT_RUNTIME
    assert 'data-nextclade-variant-tab="low" aria-pressed="false"' in REPORT_RUNTIME
    assert 'data-rsv-variant-tab="low" aria-pressed="false"' in REPORT_RUNTIME
    assert 'button.setAttribute("aria-pressed", String(button.getAttribute("data-nextclade-variant-tab") === tabKey));' in REPORT_RUNTIME
    assert 'button.setAttribute("aria-pressed", String(button.getAttribute("data-rsv-variant-tab") === tabKey));' in REPORT_RUNTIME
    assert "margin: 16px 0 12px;" in REPORT_CSS
    assert "padding-top: 10px;" in REPORT_CSS

def test_coverage_depth_controls_share_a_compact_toolbar_system() -> None:
    assert 'class="subreport-tabs ncov-depth-tabs ncov-depth-mode" role="group"' in REPORT_RUNTIME
    assert REPORT_RUNTIME.count("ncov-depth-button") == 3
    assert 'data-ncov-depth-mode="raw" aria-pressed="${depthMode === "raw"}"' in REPORT_RUNTIME
    assert 'aria-label="返回全基因组范围"' in REPORT_RUNTIME
    assert ".ncov-plot-toolbar .ncov-depth-tabs" in REPORT_CSS
    assert "min-height: 40px;" in REPORT_CSS
    assert "min-width: 96px;" in REPORT_CSS
    assert ".ncov-coverage-reset.table-export-button:focus-visible" in REPORT_CSS
    assert "background: #3d5b79;" in REPORT_CSS
    assert "color: #fff;" in REPORT_CSS
    assert "background: #f6f7f8;" in REPORT_CSS
    assert ".report-topbar" in REPORT_CSS
    assert "width: min(100%, 1280px);" in REPORT_CSS

