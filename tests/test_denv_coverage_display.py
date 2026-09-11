from pathlib import Path
import unittest

from bac_analysis_portal.serotype_reports import _build_serotype_section, _format_nextclade_coverage_percent


class DenvCoverageDisplayTests(unittest.TestCase):
    def test_fraction_is_shown_as_percent_with_two_decimal_places(self):
        self.assertEqual(_format_nextclade_coverage_percent("0.9875105643722415"), "98.75%")
        self.assertEqual(_format_nextclade_coverage_percent("1"), "100.00%")

    def test_denv_summary_card_and_table_share_percent_display(self):
        section = _build_serotype_section(Path("demo_data/denv_demo"), "denv_demo", {})
        card = next(item for item in section["summary_cards"] if item["label"] == "覆盖度")
        coverage_index = section["columns"].index("coverage")

        self.assertEqual(section["mode"], "denv_nextclade")
        self.assertEqual(card["value"], "98.75%")
        self.assertEqual(section["rows"][0][coverage_index], "98.75%")


if __name__ == "__main__":
    unittest.main()
