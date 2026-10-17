from pathlib import Path
import tempfile
import unittest
from unittest.mock import Mock

import probe_nonhuman_d3 as probe
from resume_markdup_markers import StartupWait, wait_for_startup


class SchedulerStartup(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.audit = Path(temp.name) / "scheduler.jsonl"
        self.pause = Mock()

    def call(self, responses, attempts=30, script="bjobs 42 2>/dev/null || true\n"):
        remote = Mock(side_effect=responses)
        result = StartupWait(remote, self.audit, self.pause, attempts)(
            "swls", ["bash", "-s"], script=script)
        return result, remote

    def test_starting_then_actual_running_status(self):
        result, remote = self.call(["42 STARTING user\n", "42 RUN user\n"])
        self.assertEqual(result, "42 RUN user\n")
        self.assertEqual(remote.call_count, 2)
        self.pause.assert_called_once_with(10)
        self.assertIn("STARTING", self.audit.read_text())

    def test_failure_not_hidden(self):
        result, _ = self.call(["42 STARTING user\n", "42 EXIT user\n"])
        self.assertEqual(result, "42 EXIT user\n")

    def test_bounded_startup_wait(self):
        with self.assertRaisesRegex(RuntimeError, "STARTING persisted"):
            self.call(["42 STARTING user\n"] * 3, attempts=3)
        self.assertEqual(self.pause.call_count, 2)

    def test_non_starting_states_unchanged(self):
        for status in ("PEND", "RUN", "EXITING", "DONE", "EXIT"):
            value = "42 " + status + " user\n"
            result, remote = self.call([value])
            self.assertEqual(result, value)
            self.assertEqual(remote.call_count, 1)

    def test_other_queries_unchanged(self):
        value = "42 STARTING user\n"
        result, remote = self.call([value], script="cat log\n")
        self.assertEqual(result, value)
        self.assertEqual(remote.call_count, 1)
        self.pause.assert_not_called()

    def test_frozen_remote_restored_on_error(self):
        original = probe.remote
        with self.assertRaisesRegex(ValueError, "test"):
            with wait_for_startup(self.audit):
                self.assertIsInstance(probe.remote, StartupWait)
                raise ValueError("test")
        self.assertIs(probe.remote, original)


if __name__ == "__main__":
    unittest.main()
