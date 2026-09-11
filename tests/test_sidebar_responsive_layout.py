from __future__ import annotations

import unittest
from pathlib import Path


class SidebarResponsiveLayoutTests(unittest.TestCase):
    def test_sidebar_uses_adaptive_menu_height_and_local_scrolling(self) -> None:
        styles = (Path(__file__).parents[1] / "bac_analysis_portal" / "static" / "styles.css").read_text(encoding="utf-8")

        self.assertIn("grid-template-rows: auto auto minmax(0, 1fr)", styles)
        self.assertIn("overflow-y: hidden", styles)
        self.assertIn(".app-nav::-webkit-scrollbar", styles)
        self.assertIn("justify-content: flex-start", styles)
        self.assertIn(".tab-bar.app-nav", styles)
        self.assertIn("flex-wrap: nowrap", styles)
        self.assertIn("overflow-y: auto", styles)
        self.assertIn("background: transparent", styles)
        self.assertIn("flex: 0 0 52px", styles)
        self.assertIn("min-height: 52px", styles)
        self.assertIn("max-height: 52px", styles)


if __name__ == "__main__":
    unittest.main()
