from pathlib import Path
import tempfile
import unittest

from markdup_marker_recheck import MarkerCheck


class MarkerConfiguration(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.check = object.__new__(MarkerCheck)
        self.check.root = Path(temp.name)
        (self.check.root / "logs").mkdir()
        self.check.sw = "/frozen/markers"
        self.calls = []
        self.log = "boundary halo=fallback\n" + "lookup=indexed index_peak=100\n" * 6
        self.log += "".join("[md-boundary-state rank=%d] mode=incremental\n"
                            "[md-window-reuse rank=%d] marker_visits=100 marker_retained=10 marker_peak=5\n"
                            % (i, i) for i in range(6))

        def run(label, capacity, argv, expected, check):
            self.calls.append((label, capacity, argv))
            (self.check.root / "logs" / (label + ".log")).write_text(self.log)

        self.check.run = run

    def test_full_output_and_capacity(self):
        self.check.markdup("D3", dump=True)
        self.assertEqual(self.calls[-1][1], "24G")
        self.assertNotIn("--no-dump", self.calls[-1][2])

    def test_missing_counter_rejected(self):
        self.log = self.log.replace("[md-window-reuse rank=5]", "[missing rank=5]")
        with self.assertRaisesRegex(RuntimeError, "counters"):
            self.check.markdup("D3")

    def test_no_compaction_rejected(self):
        self.log = self.log.replace("marker_retained=10 ", "marker_retained=100 ")
        with self.assertRaisesRegex(RuntimeError, "did not exercise"):
            self.check.markdup("D3")

    def test_normal_path(self):
        self.log = "boundary halo=enabled\n" + "lookup=linear index_peak=0\n" * 6
        self.check.markdup("D1")
        self.assertEqual(self.calls[-1][1], "8G")


if __name__ == "__main__":
    unittest.main()
