#!/usr/bin/env python3
"""Small real-Samtools regression for the A4 retuning driver."""

import argparse
import csv
import json
import os
from pathlib import Path
import subprocess
import tempfile
from unittest.mock import Mock, patch

import bench_samtools_a4_retune as driver


class SmokeBenchmark(driver.Benchmark):
    def wait_resources(self, **kwargs):
        # The tiny fixture does not need the production 100-GiB disk reserve.
        pass


def test_resource_guard():
    bench = driver.Benchmark.__new__(driver.Benchmark)
    bench.args = argparse.Namespace(resource_wait_seconds=1800)
    bench.root = Path("/bench")
    bench.log = Mock()
    with patch.object(driver.shutil, "disk_usage", return_value=argparse.Namespace(free=1024 * driver.GIB)), \
            patch.object(driver.os, "getloadavg", return_value=(0, 0, 0)), \
            patch.object(driver, "memory_available", return_value=182 * driver.GIB), \
            patch.object(driver.time, "sleep") as sleep:
        for dataset in ("D1", "D2", "D3"):
            bench.wait_resources(memory_bytes=driver.required_memory(dataset, "workflow"))
        sleep.assert_not_called()
        with patch.object(driver, "memory_available", return_value=80 * driver.GIB), \
                patch.object(driver.time, "monotonic", side_effect=[0, 0, 1800]):
            try:
                bench.wait_resources(memory_bytes=driver.required_memory("D1", "sort"))
            except RuntimeError as error:
                assert "resource wait timed out" in str(error)
            else:
                raise AssertionError("resource waits must have a deadline")
        sleep.assert_called_once_with(60)
        assert "required=" in bench.log.call_args[0][0]
        with patch.object(driver.shutil, "disk_usage", return_value=argparse.Namespace(free=0)):
            try:
                bench.wait_resources()
            except RuntimeError as error:
                assert "insufficient free disk" in str(error)
            else:
                raise AssertionError("disk protection must remain enabled")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samtools", required=True)
    parser.add_argument("--sam", required=True)
    parser.add_argument("--libdeflate", required=True)
    parser.add_argument("--directory", required=True)
    args = parser.parse_args()
    driver.self_test()
    test_resource_guard()
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = args.libdeflate + ":" + env.get("LD_LIBRARY_PATH", "")
    with tempfile.TemporaryDirectory(prefix="retune-smoke-", dir=args.directory) as directory:
        base = Path(directory)
        src = base / "data" / "D2-smoke" / "initial" / "D2-smoke.coordinate.bam"
        src.parent.mkdir(parents=True)
        subprocess.run([args.samtools, "sort", "--no-PG", "-m", "16M", "-l", "1",
                        "-o", str(src), args.sam], env=env, check=True)
        records = int(subprocess.check_output([args.samtools, "view", "-c", str(src)],
                                             env=env, universal_newlines=True))
        driver.DATASETS["D2"] = dict(label="D2-smoke", bytes=src.stat().st_size,
                                     records=records, sort_mib=192)
        driver.CANDIDATES = [2, 4]
        driver.OLD_THREADS = dict((app, 4) for app in driver.APPS)
        bench_args = argparse.Namespace(root=str(base / "run"), data=str(base / "data"),
                                        samtools=args.samtools, libdeflate=args.libdeflate,
                                        workflows=True)
        bench = SmokeBenchmark(bench_args)
        try:
            for index, app in enumerate(driver.APPS):
                if app == "bam2sam":
                    for repeat in range(1, 4):
                        bench.workflow("D2", repeat)
                    for name in ("name.bam", "fixmate.bam", "fixmate.sorted.bam"):
                        (bench.root / "prepared" / "D2" / name).unlink()
                bench.tune("D2", app, index)
            with open(str(bench.root / "results.tsv")) as source:
                rows = list(csv.DictReader(source, delimiter="\t"))
            assert set(row["app"] for row in rows) == set(driver.APPS + ["workflow"])
            assert "stats" not in set(row["app"] for row in rows)
            commands = [json.loads(line) for line in
                        (bench.root / "commands.jsonl").read_text().splitlines()]
            collates = [entry["argv"] for entry in commands if entry["argv"][1] == "collate"]
            assert collates
            assert all(command[command.index("-n") + 1] == "66" for command in collates)
            for app in driver.APPS:
                selected = bench.selected[("D2", app)]
                selected_rows = [row for row in rows if row["phase"] == "confirm"
                                 and row["app"] == app and int(row["threads"]) == selected]
                assert len(selected_rows) == 3
            workflows = [row for row in rows if row["app"] == "workflow"]
            assert len(workflows) == 3
            assert all(int(row["records_checked"]) == records for row in workflows)
            assert not list((bench.root / "work").iterdir())
            assert (bench.root / "prepared" / "D2" / "input.sam").exists()
            sort_log = base / "sort.stderr"
            sort_log.write_text("[bam_sort_core] merging from 3 files and 4 in-memory blocks...\n")
            bench.check_sort("D3", sort_log)
            try:
                bench.check_sort("D2", sort_log)
            except RuntimeError:
                pass
            else:
                raise AssertionError("D2 spill should be rejected")
            sort_log.write_text("")
            try:
                bench.check_sort("D3", sort_log)
            except RuntimeError:
                pass
            else:
                raise AssertionError("D3 in-memory sort should be rejected")
            resumed_args = argparse.Namespace(**vars(bench_args))
            resumed_args.root = str(base / "resumed")
            resumed_args.resume_from = str(bench.root)
            bench.writer_file.close()
            bench.lock.close()
            resumed = SmokeBenchmark(resumed_args)
            try:
                assert resumed.cached_rows
                assert any(key[1] == "collate" for key in resumed.cached_rows)
                assert not any(key[1] in ("workflow", "stats")
                               for key in resumed.cached_rows)
                imported = len(resumed.cached_rows)
                resumed.tune("D2", "flagstat", 0)
                assert len(resumed.cached_rows) == imported
                assert not (resumed.root / "commands.jsonl").exists()
                assert resumed.selected[("D2", "flagstat")] == bench.selected[("D2", "flagstat")]
                resumed.tune("D2", "collate", 2)
                assert not (resumed.root / "commands.jsonl").exists()
                for app in ("fixmate", "sort", "markdup"):
                    resumed.selected[("D2", app)] = bench.selected[("D2", app)]
                resumed.workflow("D2", 1)
                resumed_commands = [json.loads(line) for line in
                                    (resumed.root / "commands.jsonl").read_text().splitlines()]
                assert all(entry["argv"][entry["argv"].index("-n") + 1] == "66"
                           for entry in resumed_commands if entry["argv"][1] == "collate")
                assert not (resumed.root / "prepared" / "D2" / "name.bam").exists()
                print("RESUME PASSED: matching-bin collate reused; workflows regenerated", flush=True)
            finally:
                resumed.writer_file.close()
                resumed.lock.close()
            plan_path = bench.root / "plan.json"
            old_plan = json.loads(plan_path.read_text())
            old_plan["collate_bins"] = 64
            plan_path.write_text(json.dumps(old_plan) + "\n")
            resumed_args.root = str(base / "resumed_old_bins")
            old_bins = SmokeBenchmark(resumed_args)
            try:
                assert old_bins.cached_rows
                assert not any(key[1] in ("collate", "workflow", "stats")
                               for key in old_bins.cached_rows)
                old_bins.tune("D2", "collate", 2)
                assert (old_bins.root / "commands.jsonl").exists()
                print("OLD-BIN RESUME PASSED: changed-bin collate not reused", flush=True)
            finally:
                old_bins.writer_file.close()
                old_bins.lock.close()
            print("SMOKE PASSED: 8 apps, 3 workflows, {} records, {} result rows".format(
                records, len(rows)), flush=True)
        finally:
            bench.writer_file.close()
            bench.lock.close()


if __name__ == "__main__":
    main()
