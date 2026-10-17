#!/usr/bin/env python3
"""Configuration and timing checks for the isolated new-D3 driver."""

import argparse
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import bench_samtools_d3_mouse as driver


class MouseTests(unittest.TestCase):
    def setUp(self):
        driver.base.DATASETS[driver.DATASET] = dict(label=driver.DATASET, records=driver.RECORDS,
                                                   bytes=driver.INPUT_BYTES, sort_mib=98304)
        self.bench = driver.MouseBenchmark.__new__(driver.MouseBenchmark)
        self.bench.args = argparse.Namespace(input="/data/mouse.bam", samtools="samtools")
        self.bench.root = Path("/bench")

    def test_scope_and_threads(self):
        self.assertEqual(len(driver.base.APPS), 8)
        self.assertNotIn("stats", driver.base.APPS)
        self.assertEqual(len(driver.LIBRARY_THREADS), 12)
        self.assertEqual(driver.LIBRARY_THREADS[0], 0)
        self.assertEqual(driver.LIBRARY_THREADS[-1], 127)

    def test_new_species_not_old_human_d3(self):
        self.assertNotEqual(driver.DATASET, "D3")
        self.assertEqual(str(self.bench.input(driver.DATASET)), "/data/mouse.bam")
        with self.assertRaises(RuntimeError):
            self.bench.input("D3")

    def test_collate_and_sort(self):
        command = self.bench.command(driver.DATASET, "collate", 96, "out.bam", "temp")
        self.assertEqual(command[command.index("-n") + 1], "66")
        command = self.bench.command(driver.DATASET, "sort", 24, "out.bam", "temp")
        self.assertEqual(command[command.index("-m") + 1], "4096M")
        self.assertEqual(command[command.index("-l") + 1], "1")

    def test_sort_must_prove_no_spill(self):
        with tempfile.TemporaryDirectory() as directory:
            stderr = Path(directory) / "sort.stderr"
            stderr.write_text("[bam_sort_core] merging from 0 files and 24 in-memory blocks...\n")
            self.bench.check_sort(driver.DATASET, stderr)
            for text in ("", "merging from 2 files and 24 in-memory blocks"):
                stderr.write_text(text)
                with self.assertRaises(RuntimeError):
                    self.bench.check_sort(driver.DATASET, stderr)

    def test_workflow_uses_stage_sum_not_runner_time(self):
        with tempfile.TemporaryDirectory() as directory:
            self.bench.root = Path(directory)
            (self.bench.root / "logs").mkdir()
            (self.bench.root / "logs" / ("workflow." + driver.DATASET + ".r1.json")).write_text(
                json.dumps(dict(stage_times_s=[1, 2, 3, 4], wall_s="12.0")))
            with patch.object(driver.base.Benchmark, "workflow"):
                self.bench.workflow(driver.DATASET, 1)
            row = json.loads((self.bench.root / "workflows.jsonl").read_text())
            self.assertEqual(row["command_wall_sum_s"], 10)
            self.assertEqual(row["runner_wall_s"], 12)


if __name__ == "__main__":
    unittest.main()
