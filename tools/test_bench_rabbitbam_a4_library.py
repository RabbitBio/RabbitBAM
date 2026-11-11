#!/usr/bin/env python3
"""Offline checks for the isolated RabbitBAM tuning experiment."""
import unittest
from unittest.mock import Mock

import bench_rabbitbam_a4_library as bench


class Tests(unittest.TestCase):
    def test_search_limits_and_distinct_configs(self):
        for mode in ("count", "copy"):
            configs = bench.screen_configs(mode)
            self.assertEqual(len(configs), len(set(configs)))
            for nr, nw in configs:
                self.assertGreater(nr, 0)
                self.assertLessEqual(nr + nw, 120)
                self.assertEqual(nw == 0, mode == "count")
                for r, w in bench.refine_configs(mode, (nr, nw)):
                    self.assertGreater(r, 0)
                    self.assertLessEqual(r + w, 120)
                    self.assertEqual(w == 0, mode == "count")
                    self.assertNotEqual((r, w), (nr, nw))

    def test_ranking_uses_matching_phase_medians(self):
        rows = []
        for config, times in [((2, 8), [1, 100, 100]), ((4, 16), [10, 11, 12])]:
            for seconds in times:
                rows.append(dict(dataset="D1", mode="copy", phase="confirm", nr=config[0],
                                 nw=config[1], wall_s=seconds))
        rows.append(dict(dataset="D1", mode="copy", phase="screen", nr=2, nw=8, wall_s=0.1))
        self.assertEqual(bench.rank_configs(rows, "D1", "copy", {"confirm"})[0], (4, 16))
        self.assertEqual(bench.rank_configs(rows, "D2", "copy", {"confirm"}), [])

    def test_full_record_count_required(self):
        text = "BENCH records=250003 batches=4 nr=4 nw=16 parser_threads=4\n"
        self.assertEqual(bench.parse_metrics(text, 250003), (250003, 4))
        for invalid in ("", text + text, text.replace("250003", "250002")):
            with self.assertRaises(RuntimeError):
                bench.parse_metrics(invalid, 250003)

    def test_confirmation_and_formal_runs_are_separate(self):
        runner = bench.Benchmark.__new__(bench.Benchmark)
        runner.rows, calls = [], []
        runner.log, runner.summary = Mock(), Mock()

        def sample(dataset, mode, config, phase, repeat):
            nr, nw = config
            seconds = abs(nr - 8) + 10
            if phase == "confirm":
                seconds = {6: [11, 12, 13], 8: [5, 25, 25], 10: [15, 16, 17]}[nr][repeat - 1]
            calls.append((config, phase, repeat))
            runner.rows.append(dict(dataset=dataset, mode=mode, phase=phase,
                                    nr=nr, nw=nw, wall_s=seconds))

        runner.run_sample = sample
        runner.tune("D1", "count")
        confirm = [entry for entry in calls if entry[1] == "confirm"]
        formal = [entry for entry in calls if entry[1] == "formal"]
        self.assertEqual(len(confirm), 9)
        self.assertEqual([entry[0] for entry in confirm[:3]],
                         list(reversed([entry[0] for entry in confirm[3:6]])))
        self.assertEqual(formal, [((6, 0), "formal", r) for r in (1, 2, 3)])
        runner.summary.assert_called_once_with()


if __name__ == "__main__":
    unittest.main()
