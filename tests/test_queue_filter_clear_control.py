from __future__ import annotations

import unittest
from pathlib import Path


class QueueFilterClearControlTests(unittest.TestCase):
    def test_clear_button_uses_a_compact_svg_icon_and_accessible_states(self) -> None:
        root = Path(__file__).parents[1]
        runtime = (root / "bac_analysis_portal" / "static" / "app.js").read_text(encoding="utf-8")
        styles = (root / "bac_analysis_portal" / "static" / "styles.css").read_text(encoding="utf-8")

        self.assertIn('<svg class="table-filter-clear-icon"', runtime)
        self.assertIn('viewBox="0 0 16 16"', runtime)
        self.assertIn('aria-label="清空筛选"', runtime)
        self.assertNotIn('>×</button>', runtime)
        self.assertIn('.table-filter-clear-icon', styles)
        self.assertIn('grid-template-areas: "control";', styles)
        self.assertIn('.table-filter-wrap > .queue-table-filter-input', styles)
        self.assertIn('grid-area: control;', styles)
        self.assertIn('align-self: center;', styles)
        self.assertIn('justify-self: end;', styles)
        self.assertIn('margin-right: 5px;', styles)
        clear_button_block = styles.split('.table-filter-clear-button {', 1)[1].split('}', 1)[0]
        self.assertIn('border: none;', clear_button_block)
        self.assertIn('width: 28px;', clear_button_block)
        self.assertNotIn('top:', clear_button_block)
        self.assertNotIn('transform:', clear_button_block)
        self.assertNotIn('.table-filter-clear-button:hover', styles)
        self.assertNotIn('.table-filter-clear-button:active', styles)
        self.assertIn('.table-filter-clear-button:focus-visible', styles)
        self.assertIn('visibility: hidden;', styles)


if __name__ == "__main__":
    unittest.main()
