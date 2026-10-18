import csv
from pathlib import Path
import tempfile
import unittest

import probe_nonhuman_d3 as probe
from ram_d3_mouse_matrix import EXPECTED, INPUTS, Matrix, RECORDS, SDK_POOLS
from ram_d3_mouse_matrix import build_plan, sample_metrics, summaries


def log(sample):
    ranks = sample["ranks"]
    lines = ["sdk_count ranks=%d records=%d core=1.0" % (ranks, RECORDS),
             "selected mode=memory", "kept_records=134417603", "duplicates marked=33171753",
             "SN\tis sorted:\t1", "SN\treads duplicated:\t33171753",
             "boundary halo=fallback"]
    for rank in range(ranks):
        lines += ["[swbam-processing rank=%d] storage=ram seconds=2.0 status=0" % rank,
                  "[swbam-e2e rank=%d] storage=ram seconds=9.0 status=0" % rank,
                  "[ram-storage rank=%d] capacity=%d peak=1024 output_bytes=123 "
                  "temp_read_bytes=0 temp_write_bytes=0" % (rank, int(sample["capacity"][:-1]) * 1024 ** 3),
                  "[md-window-reuse rank=%d] marker_visits=10 marker_retained=1 marker_peak=1" % rank]
    return "\n".join(lines) + "\n"


class Plan(unittest.TestCase):
    def setUp(self):
        self.plan = build_plan("/frozen/D3")

    def test_counts_and_unique_samples(self):
        self.assertEqual(len(self.plan), 74)
        self.assertEqual(len({row["label"] for row in self.plan}), 74)
        for kind, count in (("sdk", 30), ("app", 24), ("workflow", 12), ("validation", 8)):
            self.assertEqual(sum(row["kind"] == kind for row in self.plan), count)

    def test_resource_points_and_order(self):
        for repeat, order in ((1, [1, 2, 4, 6, 12]), (2, [12, 6, 4, 2, 1]), (3, [1, 2, 4, 6, 12])):
            rows = [row for row in self.plan if row["kind"] == "sdk" and row["repeat"] == repeat]
            self.assertEqual([row["ranks"] for row in rows if row["app"] == "count"], order)
            for row in rows:
                self.assertEqual(row["capacity"], SDK_POOLS[row["ranks"]][row["app"] == "bam2bam"])
        self.assertTrue(all(row["ranks"] == 6 for row in self.plan if row["kind"] != "sdk"))

    def test_fixed_versions_and_memory_policy(self):
        for row in self.plan:
            if row["kind"] not in ("app", "workflow"):
                continue
            argv = row["argv"]
            if row["app"] == "fixmate":
                self.assertTrue(argv[0].endswith("/RabbitBAM-MPI-fixmate"))
                self.assertIn("-m", argv)
            elif row["app"] == "markdup":
                self.assertTrue(argv[0].endswith("/RabbitBAM-MPI-markdup"))
                self.assertEqual(argv[argv.index("-m") + 1], "8G")
            else:
                self.assertTrue(argv[0].endswith("/RabbitBAM-MPI"))
            if row["app"] in ("sort", "collate"):
                self.assertEqual(argv[argv.index("-m") + 1], "48G")
                self.assertEqual(row["check"], "memory")
            if row["app"] == "collate":
                self.assertEqual(argv[argv.index("-n") + 1], "66")

    def test_workflow_really_reloads_each_previous_output(self):
        for repeat in (1, 2, 3):
            stages = [row for row in self.plan if row["kind"] == "workflow" and row["repeat"] == repeat]
            self.assertEqual([row["app"] for row in stages], ["collate", "fixmate", "sort", "markdup"])
            source = INPUTS["bam"]
            for row in stages:
                argv = row["argv"]
                self.assertNotIn("--no-dump", argv)
                self.assertEqual(argv[argv.index("-i") + 1], source)
                source = argv[argv.index("-o") + 1]
            checks = [row["check"] for row in self.plan if row["kind"] == "validation" and
                      row["repeat"] == repeat and row["app"].startswith("workflow_")]
            self.assertEqual(checks, ["sorted", "duplicates"])

    def test_first_writer_outputs_saved(self):
        writers = [row for row in self.plan if row["kind"] == "sdk" and row["app"] == "bam2bam"]
        self.assertEqual([(row["ranks"], row["repeat"]) for row in writers if "--no-dump" not in row["argv"]], [(1, 1)])
        conversions = [row for row in self.plan if row["kind"] == "app" and row["app"] == "sam2bam"]
        self.assertNotIn("--no-dump", conversions[0]["argv"])
        self.assertTrue(all("--no-dump" in row["argv"] for row in conversions[1:]))


class Acceptance(unittest.TestCase):
    def setUp(self):
        self.plan = build_plan("/frozen/D3")

    def sample(self, name):
        return next(row for row in self.plan if row["app"] == name)

    def test_log_acceptance_and_core(self):
        sample = self.sample("bam2bam")
        result = sample_metrics(log(sample), sample)
        self.assertEqual(result["processing_max_s"], 2)
        self.assertEqual(result["core_s"], 1)

    def test_full_output_not_optional(self):
        sample = self.sample("bam2bam")
        with self.assertRaisesRegex(RuntimeError, "complete RAM output"):
            sample_metrics(log(sample).replace("output_bytes=123", "output_bytes=0"), sample)

    def test_capacity_and_unique_rank(self):
        sample = self.sample("count")
        with self.assertRaisesRegex(RuntimeError, "capacity"):
            sample_metrics(log(sample).replace("capacity=17179869184", "capacity=1"), sample)
        with self.assertRaises(RuntimeError):
            sample_metrics(log(sample).replace("processing rank=0", "processing rank=1"), sample)

    def test_external_sort_and_temp_traffic_rejected(self):
        sample = self.sample("sort")
        for value in (log(sample).replace("selected mode=memory", "selected mode=external"),
                      log(sample).replace("temp_write_bytes=0", "temp_write_bytes=1")):
            with self.assertRaises(RuntimeError):
                sample_metrics(value, sample)

    def test_markdup_spill_allowed_but_path_and_counts_required(self):
        sample = self.sample("markdup")
        value = log(sample).replace("temp_write_bytes=0", "temp_write_bytes=1")
        self.assertEqual(sample_metrics(value, sample)["temp_write_bytes"], 6)
        for bad in (value.replace("marked=33171753", "marked=1"),
                    value.replace("halo=fallback", "halo=enabled"),
                    value.replace("[md-window-reuse rank=5]", "[missing rank=5]")):
            with self.assertRaises(RuntimeError):
                sample_metrics(bad, sample)

    def test_reread_duplicates_required(self):
        sample = self.sample("workflow_markdup")
        with self.assertRaisesRegex(RuntimeError, "duplicate total"):
            sample_metrics(log(sample).replace("reads duplicated:\t33171753", "reads duplicated:\t1"), sample)


class Summary(unittest.TestCase):
    def test_round_sum_then_median_not_sum_of_stage_medians(self):
        rows = []
        for repeat, times in enumerate(((1, 100, 1, 1), (100, 1, 1, 1), (100, 100, 1, 1)), 1):
            for app, value in zip(("collate", "fixmate", "sort", "markdup"), times):
                rows.append(dict(kind="workflow", app=app, repeat=repeat, ranks=6, job="42",
                                 processing_max_s=value, e2e_max_s=value + 7))
            rows.append(dict(kind="validation", app="workflow_markdup", repeat=repeat, ranks=6,
                             processing_max_s=1000, e2e_max_s=1007))
            rows.append(dict(kind="validation", app="workflow_sort", repeat=repeat, ranks=6,
                             processing_max_s=1000, e2e_max_s=1007))
        result = summaries(rows)
        self.assertEqual(result["workflow_processing_median_s"], 103)
        self.assertEqual(result["workflow_mpi_interval_median_s"], 131)

    def test_no_partial_series_median(self):
        row = dict(kind="sdk", app="count", repeat=1, ranks=6, core_s=1,
                   processing_max_s=2, e2e_max_s=9)
        self.assertEqual(summaries([row])["series"], [])
        with self.assertRaisesRegex(RuntimeError, "duplicate repetitions"):
            summaries([row, row, row])

    def test_full_mock_matrix_does_not_submit_cluster_jobs(self):
        with tempfile.TemporaryDirectory() as directory:
            matrix = object.__new__(Matrix)
            matrix.root = Path(directory)
            matrix.sw = "/frozen/D3"
            (matrix.root / "logs").mkdir()
            plan = build_plan(matrix.sw)
            by_label = {row["label"]: row for row in plan}
            matrix.initialize = lambda: plan
            states = []
            matrix.state = lambda phase, **values: states.append(phase)

            def run(label, capacity, argv, expected, check, ranks):
                self.assertEqual(expected, EXPECTED)
                sample = by_label[label]
                self.assertEqual(argv, sample["argv"])
                value = log(sample)
                (matrix.root / "logs" / (label + ".log")).write_text(value)
                result = probe.metrics(value, ranks, RECORDS, check, EXPECTED["kept"], EXPECTED["marked"])
                result.update(label=label, job="42", ranks=ranks, capacity=capacity)
                with (matrix.root / "results.tsv").open("a", newline="") as stream:
                    writer = csv.DictWriter(stream, list(result), delimiter="\t")
                    if stream.tell() == 0:
                        writer.writeheader()
                    writer.writerow(result)

            matrix.run = run
            matrix.execute()
            self.assertEqual(states[-1], "COMPLETE")
            self.assertEqual(states.count("SAMPLE_ACCEPTED"), 74)
            with (matrix.root / "accepted_samples.tsv").open() as stream:
                self.assertEqual(len(list(csv.DictReader(stream, delimiter="\t"))), 74)


if __name__ == "__main__":
    unittest.main()
