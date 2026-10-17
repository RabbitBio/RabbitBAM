import tempfile
from pathlib import Path
from types import SimpleNamespace
import unittest

from markdup_lookup_recheck import EXPECTED, INPUTS, Recheck


class ConfigurationTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.check = object.__new__(Recheck)
        self.check.root = Path(self.temp.name)
        (self.check.root / "logs").mkdir()
        self.check.sw = "/frozen/check"
        self.calls = []
        self.log = "boundary halo=enabled\n" + "lookup=linear index_peak=0\n" * 6

        def run(label, capacity, argv, expected, check):
            self.calls.append((label, capacity, argv, expected, check))
            (self.check.root / "logs" / (label + ".log")).write_text(self.log)

        self.check.run = run

    def test_normal_configuration(self):
        self.check.markdup("D1")
        label, capacity, argv, expected, check = self.calls[-1]
        self.assertEqual((label, capacity, check), ("D1_new", "8G", "markdup"))
        self.assertEqual(expected, EXPECTED["D1"])
        self.assertIn(INPUTS["D1"], argv)
        self.assertIn("--no-dump", argv)
        self.assertEqual(argv[argv.index("-m") + 1], "8G")

    def test_normal_path_rejects_index(self):
        self.log = "boundary halo=enabled\n" + "lookup=indexed index_peak=2\n" * 6
        with self.assertRaisesRegex(RuntimeError, "normal path"):
            self.check.markdup("D2")

    def test_fallback_saved_output(self):
        self.log = "boundary halo=fallback\n" + "lookup=indexed index_peak=300\n" * 6
        output = self.check.markdup("D3", dump=True)
        self.assertEqual(output, "/frozen/check/out/D3_new.bam")
        self.assertEqual(self.calls[-1][1], "24G")
        self.assertNotIn("--no-dump", self.calls[-1][2])

    def test_fallback_requires_index(self):
        with self.assertRaisesRegex(RuntimeError, "fallback index"):
            self.check.markdup("D3")


if __name__ == "__main__":
    unittest.main()
