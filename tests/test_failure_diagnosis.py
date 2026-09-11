from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from bac_analysis_portal.failure_diagnosis import build_failure_diagnosis


class FailureDiagnosisTests(unittest.TestCase):
    def _diagnose(self, log_text: str, status: str = "FAILED") -> dict:
        with TemporaryDirectory() as directory:
            log_path = Path(directory) / "pipeline.log"
            log_path.write_text(log_text, encoding="utf-8")
            return build_failure_diagnosis({"status": status}, log_path)

    def test_missing_python_module_is_classified_as_environment_failure(self) -> None:
        diagnosis = self._diagnose("ModuleNotFoundError: No module named 'pandas'")
        self.assertEqual(diagnosis["category"], "dependency_missing")
        self.assertEqual(diagnosis["label"], "运行环境依赖缺失")
        self.assertTrue(diagnosis["rerun_recommended"])
        self.assertIn("pandas", diagnosis["evidence"][0])

    def test_input_failure_warns_against_negative_interpretation(self) -> None:
        diagnosis = self._diagnose("Input FASTQ file does not exist")
        self.assertEqual(diagnosis["category"], "input_invalid")
        self.assertIn("不能据此作出阴性或排除判断", diagnosis["impact"])

    def test_unknown_failure_requires_manual_confirmation(self) -> None:
        diagnosis = self._diagnose("process exited unexpectedly")
        self.assertEqual(diagnosis["category"], "unknown")
        self.assertFalse(diagnosis["rerun_recommended"])

    def test_stopped_task_has_explicit_disposition_guidance(self) -> None:
        diagnosis = self._diagnose("", status="STOPPED")
        self.assertEqual(diagnosis["category"], "manually_stopped")
        self.assertEqual(diagnosis["confidence"], "confirmed")


if __name__ == "__main__":
    unittest.main()
