from __future__ import annotations

import unittest
from pathlib import Path


class QueueControlAffordanceTests(unittest.TestCase):
    def test_sort_and_filter_controls_expose_interaction_states(self) -> None:
        root = Path(__file__).parents[1]
        runtime = (root / "bac_analysis_portal" / "static" / "app.js").read_text(encoding="utf-8")
        styles = (root / "bac_analysis_portal" / "static" / "styles.css").read_text(encoding="utf-8")

        self.assertIn('data-queue-sort-menu', runtime)
        self.assertIn('data-queue-filter-menu', runtime)
        self.assertIn('.queue-control-trigger::after', styles)
        self.assertIn('.queue-control-menu:not([open]) .queue-control-trigger:hover', styles)
        self.assertIn('.queue-control-trigger:active', styles)
        self.assertIn('.queue-control-trigger:focus-visible', styles)
        self.assertGreaterEqual(styles.count('min-height: 44px'), 2)
        self.assertIn('padding: 5px 28px 5px 9px;', styles)
        self.assertIn('@media (hover: hover) and (pointer: fine)', styles)
        self.assertIn('@media (prefers-reduced-motion: reduce)', styles)
        self.assertIn('transform: translateY(-25%) rotate(225deg);', styles)


if __name__ == "__main__":
    unittest.main()
