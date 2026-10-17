from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from markdup_lookup_recheck import Recheck
from resume_markdup_lookup import Recovery


class RecoveryTest(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.recovery = object.__new__(Recovery)
        self.recovery.root = Path(temp.name)
        (self.recovery.root / "logs").mkdir()
        self.calls = []

    def run_record(self, instance, label, capacity, argv, expected, check, ranks):
        self.calls.append((label, capacity, argv, expected, check, ranks))
        (instance.root / "logs" / (label + ".log")).write_text("verified retry\n")

    def test_only_reread_gets_a_new_session_label(self):
        argv = ["frozen/stats", "--basic", "-i", "saved.bam"]
        with patch.object(Recheck, "run", autospec=True, side_effect=self.run_record):
            self.recovery.run("D3_validate", "8G", argv, {"records": 10}, "sorted")
        self.assertEqual(self.calls, [("D3_validate_retry1", "8G", argv, {"records": 10}, "sorted", 6)])
        self.assertEqual((self.recovery.root / "logs/D3_validate.log").read_text(), "verified retry\n")

    def test_controls_keep_original_configuration(self):
        argv = ["frozen/markdup", "--no-dump"]
        with patch.object(Recheck, "run", autospec=True, side_effect=self.run_record):
            self.recovery.run("D1_new", "8G", argv, {"records": 10}, "markdup")
        self.assertEqual(self.calls[0][0], "D1_new")
        self.assertFalse((self.recovery.root / "logs/D3_validate.log").exists())

    def test_failed_retry_is_not_published_as_validated(self):
        with patch.object(Recheck, "run", side_effect=RuntimeError("incomplete output")):
            with self.assertRaisesRegex(RuntimeError, "incomplete"):
                self.recovery.run("D3_validate", "8G", [], {"records": 10}, "sorted")
        self.assertFalse((self.recovery.root / "logs/D3_validate.log").exists())


if __name__ == "__main__":
    unittest.main()
