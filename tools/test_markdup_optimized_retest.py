#!/usr/bin/env python3
"""Offline plan, timing and external-run acceptance checks for the selective retest."""

import copy
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import markdup_optimized_retest as retest


ROOT = "/test/frozen"


def external_log():
    lines = ["OptimizedSortMPI finished. ranks=6 records=100",
             "Complete the 4.3 OptimizedSortMPI cost 1.0",
             "MPI BAM sort selected mode=external",
             "  external_actual_sum runs=18",
             "  external_memory_sum tracked_peak_max=8589934592"]
    for rank in range(6):
        lines.extend([
            "[swbam-processing rank=%d] seconds=2.0 status=0" % rank,
            "[swbam-e2e rank=%d] seconds=9.0 status=0" % rank,
            "[ram-storage rank=%d] capacity=68719476736 peak=200 output_bytes=456 "
            "temp_read_bytes=100 temp_write_bytes=100" % rank])
    return "\n".join(lines)


class RetestTests(unittest.TestCase):
    def test_only_affected_series(self):
        counts = []
        for case in retest.SPECS:
            plan = retest.build_plan(ROOT, case)
            counts.append(len(plan))
            self.assertEqual(len({row["label"] for row in plan}), len(plan))
            apps = [row for row in plan if row["kind"] == "app"]
            self.assertEqual({row["app"] for row in apps}, set() if case == "D2_external" else {"markdup"})
            self.assertEqual(len(apps), 0 if case == "D2_external" else 3)
            self.assertEqual(len([row for row in plan if row["kind"] == "workflow"]), 12)
            for row in plan:
                if row["app"] == "markdup":
                    self.assertTrue(row["argv"][0].endswith("/RabbitBAM-MPI-markdup"))
                    self.assertIn("--stream", row["argv"])
                    self.assertEqual(row["argv"][row["argv"].index("-l") + 1], "300")
                for path in row["cleanup"]:
                    self.assertTrue(path.startswith(ROOT + "/out/workflow_r"))
            self.assertEqual(len([row for row in plan if row["check"] == "external"]),
                             6 if case == "D2_external" else 0)
        self.assertEqual(sum(counts), 84)

    def test_real_round_dependencies_and_readback(self):
        for case in retest.SPECS:
            plan = retest.build_plan(ROOT, case)
            for repeat in (1, 2, 3):
                stages = [row for row in plan if row["repeat"] == repeat and row["kind"] == "workflow"]
                self.assertEqual([row["app"] for row in stages], ["collate", "fixmate", "sort", "markdup"])
                for before, after in zip(stages, stages[1:]):
                    output = before["argv"][before["argv"].index("-o") + 1]
                    self.assertEqual(after["argv"][after["argv"].index("-i") + 1], output)
                self.assertTrue(all("--no-dump" not in row["argv"] for row in stages))
                checks = [row for row in plan if row["repeat"] == repeat and row["app"] == "workflow_markdup"]
                self.assertEqual(len(checks), 1)
                self.assertEqual(checks[0]["check"], "duplicates")

    def test_binaries_and_workspaces_are_frozen(self):
        for case in retest.SPECS:
            frozen = retest.binaries(case)
            self.assertEqual(frozen["RabbitBAM-MPI-markdup"][1], retest.NEW_SHA)
            self.assertEqual(frozen["RabbitBAM-MPI"][1], retest.BASE_SHA)
            self.assertEqual(frozen["RabbitBAM-MPI-fixmate"][1],
                retest.MOUSE_FM_SHA if case == "D3" else retest.BASE_SHA)
        self.assertEqual(retest.SPECS["D2_external"]["memories"], ("8G", "8G"))
        self.assertEqual(retest.SPECS["D2_external"]["pools"][0], "64G")

    def test_actual_bounded_external_runs_required(self):
        sample = next(row for row in retest.build_plan(ROOT, "D2_external") if row["check"] == "external")
        spec = dict(records=100, marked=7)
        self.assertEqual(retest.sample_metrics(external_log(), sample, spec)["runs"], 18)
        for old, new in [("mode=external", "mode=memory"), ("runs=18", "runs=6"),
                         ("temp_read_bytes=100", "temp_read_bytes=0"),
                         ("tracked_peak_max=8589934592", "tracked_peak_max=8589934593")]:
            with self.assertRaises(RuntimeError):
                retest.sample_metrics(external_log().replace(old, new), sample, spec)

    def test_sum_actual_stages_then_median(self):
        rows = []
        for repeat, times in enumerate(((1, 2, 3, 40), (10, 20, 30, 4), (3, 5, 7, 9)), 1):
            for app, seconds in zip(("collate", "fixmate", "sort", "markdup"), times):
                rows.append(dict(kind="workflow", app=app, repeat=repeat,
                    label="%s%d" % (app, repeat), processing_max_s=seconds, e2e_max_s=seconds + 7))
            for app in ("workflow_sort", "workflow_markdup"):
                rows.append(dict(kind="validation", app=app, repeat=repeat, label=app + str(repeat)))
        result = retest.summaries(rows)
        self.assertEqual([row["processing_total_s"] for row in result["workflows"]], [46, 64, 24])
        self.assertEqual(result["workflow_processing_median_s"], 46)
        incomplete = [row for row in rows if row["app"] != "workflow_markdup"]
        self.assertNotIn("workflow_processing_median_s", retest.summaries(incomplete))
        with self.assertRaises(RuntimeError):
            retest.summaries(rows + [copy.deepcopy(rows[0])])

    def test_submission_manifest_nodes_prefix_and_starting(self):
        for nodes, prefix in ((None, None), ("43334-43339", "mdopt_test_")):
            with self.subTest(nodes=nodes), tempfile.TemporaryDirectory() as temporary:
                driver = retest.probe.Probe.__new__(retest.probe.Probe)
                driver.root, driver.sw = Path(temporary), ROOT
                driver.args = SimpleNamespace(swls="unused", resume=False, retry_no_start=False,
                    nodes=nodes, job_prefix=prefix, poll_seconds=0)
                driver.state = lambda *args, **kwargs: None
                states = iter(("123 STARTING", "123 DONE"))

                def remote(host, argv, timeout=120, script=None):
                    if script.startswith("set -eu\ncd "):
                        return "Job <123> has been submitted"
                    if script.startswith("bjobs "):
                        return next(states)
                    return external_log()

                with patch.object(retest.probe, "remote", side_effect=remote):
                    driver.run("sample", "64G", ["/frozen/app", "sort"],
                        dict(records=100), "external")
                command = json.loads((driver.root / "jobs.jsonl").read_text())["argv"]
                if nodes:
                    self.assertEqual(command[command.index("-node") + 1], nodes)
                else:
                    self.assertNotIn("-node", command)
                self.assertEqual(command[command.index("-J") + 1],
                    "mdopt_test_sample" if prefix else "d3mouse_sample")
                self.assertIsNone(driver.current_job)


if __name__ == "__main__":
    unittest.main()
