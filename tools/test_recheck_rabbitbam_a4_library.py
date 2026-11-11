#!/usr/bin/env python3
import unittest

import recheck_rabbitbam_a4_library as check


class Tests(unittest.TestCase):
    def test_single_node_does_not_oversubscribe(self):
        self.assertEqual(check.policy_command("node0", 8, 0),
                         ["numactl", "--cpunodebind=0", "--membind=0"])
        self.assertIsNone(check.policy_command("node1", 64, 0))
        self.assertIsNone(check.policy_command("node0", 24, 96))
        self.assertEqual(check.policy_command("active", 120, 0), [])

    def test_time_includes_io_faults_and_switches(self):
        value = check.metrics("10.2\t99\t4\t220000\t0\t123\t456\t7\t0\t80000\t0\n")
        self.assertEqual(value["input_blocks"], "0")
        self.assertEqual(value["voluntary_switches"], "456")
        self.assertEqual(value["major_faults"], "0")
        with self.assertRaises(RuntimeError):
            check.metrics("10.2\t99\t4\t220000\t0\t123\t456\t7\t0\t80000\t1\n")

    def test_groups_do_not_mix_policy_or_stage(self):
        rows = []
        for phase, policy, times in (("confirm", "baseline", [10, 12, 11]),
                                     ("confirm", "active", [20, 24, 22]),
                                     ("formal", "baseline", [13, 14, 15])):
            for seconds in times:
                rows.append(dict(dataset="D1", mode="count", phase=phase, policy=policy,
                                 nr=8, nw=0, wall_s=seconds, user_s=90, sys_s=4,
                                 input_blocks=0, voluntary_switches=10, involuntary_switches=5))
        result = check.summarize(rows)
        self.assertEqual(len(result), 3)
        medians = {(r["phase"], r["policy"]): r["median_s"] for r in result}
        self.assertEqual(medians[("confirm", "baseline")], 11)
        self.assertEqual(medians[("confirm", "active")], 22)
        self.assertEqual(medians[("formal", "baseline")], 14)

    def test_low_thread_candidate_added_only_when_needed(self):
        rows = [dict(dataset="D3", mode="count", phase="confirm", nr=n, nw=0, wall_s=n / 10)
                for n in (64, 80, 120)]
        rows.append(dict(dataset="D3", mode="count", phase="screen", nr=24, nw=0, wall_s=9.63))
        self.assertEqual(check.candidate_configs(rows, "D3", "count"),
                         [(64, 0), (80, 0), (120, 0), (24, 0)])
        rows.append(dict(dataset="D3", mode="count", phase="confirm", nr=8, nw=0, wall_s=5))
        self.assertEqual(check.candidate_configs(rows, "D3", "count"), [(8, 0), (64, 0), (80, 0)])


if __name__ == "__main__":
    unittest.main()
