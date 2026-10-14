#!/usr/bin/env python3
"""Capacity gates and formal D2 samples, without real cluster submissions."""

import unittest

import test_ram_d1_matrix as fixture


class D2MatrixTest(unittest.TestCase):
    def run_matrix(self, overrides):
        return fixture.MatrixTest().run_matrix(overrides, dataset='D2')

    def test_full_matrix_after_all_probes(self):
        status, files = self.run_matrix({'PRIOR_MATRIX_STATE': 'COMPLETE'})
        self.assertEqual(status, 0, files['diagnostic'])
        rows = [line.split('\t') for line in files['results.tsv'].splitlines()[1:]]
        self.assertEqual(len(rows), 84)
        for kind, count in (('probe', 6), ('probe_validation', 4), ('sdk', 30),
                            ('app', 24), ('workflow', 12), ('validation', 8)):
            self.assertEqual(sum(row[0] == kind for row in rows), count)
        self.assertEqual([row[0] for row in rows[:10]],
                         ['probe', 'probe_validation', 'probe', 'probe_validation',
                          'probe', 'probe', 'probe_validation', 'probe', 'probe', 'probe_validation'])
        self.assertEqual(files['state'], 'COMPLETE\n')
        self.assertIn('GATES PASSED', files['capacity_gate.txt'])
        self.assertEqual(len(files['prior_probes.tsv'].splitlines()), 8)
        self.assertEqual(files['workflow_summary.txt'],
                         'workflow_e2e_median_s=60.000000\nworkflow_processing_median_s=52.000000\n')
        self.assertIn('probe\tsort\t1\t6\t100\t16G', files['results.tsv'])
        self.assertIn('probe\tbam2bam\t1\t1\t100\t76G', files['results.tsv'])
        self.assertNotIn('--no-dump', files['commands.sh'].splitlines()[0])
        self.assertIn('SWBAM_SORT_MEMORY_DIAG=1', files['commands.sh'].splitlines()[0])
        self.assertIn('-x SWBAM_SORT_MEMORY_DIAG=1', files['commands.sh'].splitlines()[0])
        self.assertNotIn(' env ', files['commands.sh'].splitlines()[0])
        for ranks, read, write in (('1', '36G', '76G'), ('2', '20G', '40G'),
                                   ('4', '12G', '24G'), ('6', '8G', '16G'), ('12', '4G', '8G')):
            self.assertTrue(all(row[5] == read for row in rows if row[:2] == ['sdk', 'count'] and row[3] == ranks))
            self.assertTrue(all(row[5] == write for row in rows if row[:2] == ['sdk', 'bam2bam'] and row[3] == ranks))
        for row in rows:
            if row[0] in ('app', 'workflow') and row[1] == 'sort':
                self.assertEqual(row[5], '16G')

    def test_sort_failure_stops_before_formal_matrix(self):
        status, files = self.run_matrix({'FAIL_LABEL': 'probe_sort_n6_r1'})
        self.assertNotEqual(status, 0)
        self.assertEqual(len(files['commands.sh'].splitlines()), 1)
        self.assertFalse(files['capacity_gate.txt'])

    def test_single_rank_writer_failure_stops_matrix(self):
        status, files = self.run_matrix({'FAIL_LABEL': 'probe_bam2bam_n1_r1'})
        self.assertNotEqual(status, 0)
        self.assertEqual(len(files['commands.sh'].splitlines()), 9)
        self.assertNotIn('sdk\tcount', files['results.tsv'])

    def test_capacity_only_does_not_submit_formal_samples(self):
        status, files = self.run_matrix({'BENCH_STAGE': 'capacity'})
        self.assertEqual(status, 0, files['diagnostic'])
        self.assertEqual(files['state'], 'CAPACITY_COMPLETE\n')
        self.assertEqual(len(files['commands.sh'].splitlines()), 10)

    def test_mismatched_baseline_rejected_without_submission(self):
        status, files = self.run_matrix({'BAD_BASELINE_HASH': '1'})
        self.assertNotEqual(status, 0)
        self.assertFalse(files['commands.sh'])
        self.assertIn('baseline_gate', files['state'])

    def test_bad_baseline_count_rejected_without_submission(self):
        status, files = self.run_matrix({'BAD_RECORDS': 'count'})
        self.assertNotEqual(status, 0)
        self.assertFalse(files['commands.sh'])

    def test_external_sort_rejected(self):
        status, files = self.run_matrix({'BAD_MODE': 'probe_sort_n6_r1'})
        self.assertNotEqual(status, 0)
        self.assertEqual(len(files['commands.sh'].splitlines()), 1)

    def test_ram_temp_data_rejected_in_this_in_memory_matrix(self):
        status, files = self.run_matrix({'BAD_TEMP': 'probe_sort_n6_r1'})
        self.assertNotEqual(status, 0)
        self.assertEqual(len(files['commands.sh'].splitlines()), 1)

    def test_markdup_ram_deferred_temp_data_allowed(self):
        for label in ('probe_markdup_n6_r1', 'app_markdup_n6_r1',
                      'workflow_markdup_n6_r1'):
            with self.subTest(label=label):
                status, files = self.run_matrix({'BAD_TEMP': label})
                self.assertEqual(status, 0, files['diagnostic'])
                self.assertEqual(files['state'], 'COMPLETE\n')
                self.assertEqual(len(files['results.tsv'].splitlines()), 85)

    def test_markdup_temp_data_does_not_bypass_output_validation(self):
        status, files = self.run_matrix({
            'BAD_TEMP': 'probe_markdup_n6_r1',
            'BAD_DUPLICATES': 'probe_validation_markdup_n6_r1',
        })
        self.assertNotEqual(status, 0)
        self.assertEqual(len(files['commands.sh'].splitlines()), 4)

    def test_wrong_duplicate_count_rejected(self):
        status, files = self.run_matrix({'BAD_DUPLICATES': 'probe_validation_markdup_n6_r1'})
        self.assertNotEqual(status, 0)
        self.assertEqual(len(files['commands.sh'].splitlines()), 4)

    def test_d1_failure_prevents_formal_matrix(self):
        status, files = self.run_matrix({'PRIOR_MATRIX_STATE': 'FAILED app_sort job=123'})
        self.assertNotEqual(status, 0)
        self.assertEqual(len(files['commands.sh'].splitlines()), 10)
        self.assertIn('waiting_for_prior_matrix', files['state'])

    def test_resume_reuses_verified_probes(self):
        status, files = self.run_matrix({'FAIL_LABEL': 'probe_sort_n6_r1', 'RECOVER': True})
        self.assertEqual(status, 0, files['diagnostic'])
        self.assertEqual(len(files['commands.sh'].splitlines()), 84)
        self.assertEqual(len(files['results.tsv'].splitlines()), 85)

    def test_separate_markdup_keeps_baseline_and_other_commands(self):
        status, files = self.run_matrix({'BENCH_MARKDUP_APP': 'separate'})
        self.assertEqual(status, 0, files['diagnostic'])
        self.assertEqual(len(files['results.tsv'].splitlines()), 85)
        commands = files['commands.sh'].splitlines()
        self.assertEqual(sum('RabbitBAM-MPI-markdup markdup' in line for line in commands), 6)
        self.assertIn('RabbitBAM-MPI markdup', commands[2])
        self.assertEqual(sum('RabbitBAM-MPI sort' in line for line in commands), 7)
        self.assertIn('RabbitBAM-MPI-markdup', files['markdup.sha256'])


if __name__ == '__main__':
    unittest.main()
