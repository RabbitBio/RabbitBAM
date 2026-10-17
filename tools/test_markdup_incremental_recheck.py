from pathlib import Path
import tempfile
import unittest

from markdup_incremental_recheck import IncrementalCheck


class IncrementalConfiguration(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.check = object.__new__(IncrementalCheck)
        self.check.root = Path(temp.name)
        (self.check.root / "logs").mkdir()
        self.check.sw = "/frozen/incremental"
        self.calls = []
        self.log = "boundary halo=fallback\n" + "lookup=indexed index_peak=200\n" * 6
        self.log += "".join("[md-boundary-state rank=%d] mode=incremental\n" % i for i in range(6))

        def run(label, capacity, argv, expected, check):
            self.calls.append((label, capacity, argv))
            (self.check.root / "logs" / (label + ".log")).write_text(self.log)

        self.check.run = run

    def test_fallback_complete_storage(self):
        self.check.markdup("D3", dump=True)
        self.assertEqual(self.calls[-1][1], "24G")
        self.assertNotIn("--no-dump", self.calls[-1][2])

    def test_reject_missing_incremental_rank(self):
        self.log = self.log.replace("[md-boundary-state rank=5] mode=incremental\n", "")
        with self.assertRaisesRegex(RuntimeError, "all ranks"):
            self.check.markdup("D3")

    def test_normal_keeps_linear_lookup(self):
        self.log = "boundary halo=enabled\n" + "lookup=linear index_peak=0\n" * 6
        self.check.markdup("D1")
        self.assertEqual(self.calls[-1][1], "8G")

    def test_reject_incremental_normal_path(self):
        self.log = "boundary halo=enabled\n" + "lookup=linear index_peak=0\n" * 6
        self.log += "[md-boundary-state rank=0] mode=incremental\n"
        with self.assertRaisesRegex(RuntimeError, "normal halo path"):
            self.check.markdup("D2")


if __name__ == "__main__":
    unittest.main()
