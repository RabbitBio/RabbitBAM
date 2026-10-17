from pathlib import Path
import csv
import json
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from confirm_markdup_markers import Confirmation, accepted_first, summarize


class ConfirmationTest(unittest.TestCase):
    def test_median_not_best_sample(self):
        result = summarize([dict(processing_max_s=value, e2e_max_s=value + 7)
                            for value in (103.5, 106.0, 101.0)])
        self.assertEqual(result["processing_median_s"], 103.5)
        self.assertEqual(result["processing_range_s"], [101.0, 106.0])
        self.assertEqual(result["mpi_interval_median_s"], 110.5)
        self.assertTrue(result["target_met"])

    def test_slow_samples_are_not_dropped(self):
        result = summarize([dict(processing_max_s=value, e2e_max_s=value + 7)
                            for value in (100, 140, 160)])
        self.assertEqual(result["processing_median_s"], 140)
        self.assertFalse(result["target_met"])

    def test_reject_incomplete_or_extra_samples(self):
        for count in (0, 1, 2, 4):
            with self.assertRaisesRegex(RuntimeError, "exactly three"):
                summarize([dict(processing_max_s=100, e2e_max_s=107)] * count)

    def test_wait_for_first_sample(self):
        with tempfile.TemporaryDirectory() as directory:
            self.assertIsNone(accepted_first(Path(directory), "frozen"))

    def test_reject_different_first_binary(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "binaries.sha256").write_text("old  /frozen/markdup-new\n")
            (root / "results.tsv").write_text("label\n")
            with self.assertRaisesRegex(RuntimeError, "different binary"):
                accepted_first(root, "new")

    def test_only_two_complete_output_confirmations(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "logs").mkdir()
            (root / "program").write_bytes(b"frozen")
            check = object.__new__(Confirmation)
            check.root, check.sw = root, "/frozen/confirm"
            check.args = SimpleNamespace(first_root=root, new_app=root / "program")
            calls, states, rows = [], [], []
            check.state = lambda phase, **values: states.append(phase)

            def run(label, capacity, argv, expected, validation):
                calls.append((label, capacity, argv, expected, validation))
                rows.append(dict(processing_max_s=104, e2e_max_s=111))
                with (root / "results.tsv").open("w") as stream:
                    writer = csv.DictWriter(stream, list(rows[0]), delimiter="\t")
                    writer.writeheader()
                    writer.writerows(rows)
                (root / "logs" / (label + ".log")).write_text(
                    "boundary halo=fallback\n" + "lookup=indexed index_peak=1\n" * 6 +
                    "".join("[md-window-reuse rank=%d]\n" % rank for rank in range(6)))

            check.run = run
            with patch("confirm_markdup_markers.accepted_first",
                       return_value=dict(processing_max_s=103, e2e_max_s=110)):
                check.checks()
            self.assertEqual([call[0] for call in calls], ["D3_confirm_r2", "D3_confirm_r3"])
            for _, capacity, argv, _, validation in calls:
                self.assertEqual((capacity, validation), ("24G", "markdup"))
                self.assertIn("--no-dump", argv)
                self.assertEqual(argv[argv.index("-m") + 1], "8G")
            self.assertEqual(states[-1], "COMPLETE")
            summary = json.loads((root / "confirmation_summary.json").read_text())
            self.assertEqual(summary["processing_median_s"], 104)


if __name__ == "__main__":
    unittest.main()
