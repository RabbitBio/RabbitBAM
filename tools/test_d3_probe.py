#!/usr/bin/env python3
"""Focused checks of the one-off probe's log acceptance, with no scheduler calls."""

import unittest

from probe_nonhuman_d3 import completed_log, metrics


def log():
    return "\n".join([
        "sdk_count ranks=2 records=100 core=1.2",
        "[swbam-e2e rank=1] seconds=10.1 status=0",
        "[swbam-e2e rank=0] seconds=10.0 status=0",
        "[swbam-processing rank=0] seconds=2.0 status=0",
        "[swbam-processing rank=1] seconds=2.1 status=0",
        "[ram-storage rank=0] peak=123 output_bytes=456 temp_read_bytes=0 temp_write_bytes=0",
        "[ram-storage rank=1] peak=124 output_bytes=456 temp_read_bytes=0 temp_write_bytes=0",
        "[ram-prepare rank=0] seconds=3.0",
    ])


class ProbeTests(unittest.TestCase):
    def test_complete_sample(self):
        result = metrics(log(), 2, 100)
        self.assertEqual(result["processing_max_s"], 2.1)
        self.assertEqual(result["pool_peak_max_bytes"], 124)

    def test_no_missing_or_duplicate_rank(self):
        with self.assertRaises(RuntimeError):
            metrics(log().replace("processing rank=1", "processing rank=0"), 2, 100)

    def test_record_count(self):
        with self.assertRaises(RuntimeError):
            metrics(log(), 2, 101)

    def test_in_memory_and_temporary_io(self):
        with self.assertRaises(RuntimeError):
            metrics(log(), 2, 100, "memory")
        with self.assertRaises(RuntimeError):
            metrics(log().replace("temp_write_bytes=0", "temp_write_bytes=1"), 2, 100)

    def test_markdup_ram_spill_and_reference(self):
        data = log().replace("temp_write_bytes=0", "temp_write_bytes=1") + "\nduplicates marked=7\n"
        result = metrics(data, 2, 100, "markdup", marked=7)
        self.assertEqual(result["temp_write_bytes"], 2)
        with self.assertRaises(RuntimeError):
            metrics(data, 2, 100, "markdup", marked=8)

    def test_timing_checked_by_rank_not_log_order(self):
        data = log().replace("processing rank=0] seconds=2.0", "processing rank=0] seconds=10.05")
        with self.assertRaises(RuntimeError):
            metrics(data, 2, 100)

    def test_log_visibility_waits_for_all_rank_fields(self):
        snapshots = iter(["", log().split("[ram-storage rank=1]")[0], log()])
        sleeps = []
        text, result = completed_log(lambda: next(snapshots), 2, {"records": 100},
                                     attempts=3, pause=sleeps.append)
        self.assertEqual(text, log())
        self.assertEqual(result["processing_max_s"], 2.1)
        self.assertEqual(sleeps, [3, 3])

    def test_persistently_incomplete_log_is_not_accepted(self):
        with self.assertRaises(RuntimeError):
            completed_log(lambda: "", 2, {"records": 100}, attempts=2, pause=lambda _: None)

    def test_application_error_is_not_retried_as_log_delay(self):
        sleeps = []
        with self.assertRaises(RuntimeError):
            completed_log(lambda: "ERROR: failed", 2, {"records": 100}, pause=sleeps.append)
        self.assertFalse(sleeps)


if __name__ == "__main__":
    unittest.main()
