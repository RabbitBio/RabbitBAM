#!/usr/bin/env python3
"""Validate the D1 matrix using a fake scheduler, without submitting jobs."""

import os
from pathlib import Path
import subprocess
import tempfile
import unittest


MOCK_SUBMIT = r'''#!/usr/bin/env python3
import os
from pathlib import Path
import sys
args = sys.argv[1:]
label = args[args.index('-J') + 1][len('d1ram_'):]
log = Path(args[args.index('-o') + 1])
ranks = int(args[args.index('--ram-ranks') + 1])
output_indices = [i for i, value in enumerate(args) if value == '-o']
if len(output_indices) > 1 and '--no-dump' not in args:
    Path(args[output_indices[1] + 1]).write_bytes(b'complete mock output')
prefix = 'sdk_count' if label.startswith('sdk_count_') else 'OptimizedMockMPI finished.'
dataset = args[args.index('-J') + 1][:2]
expected = 365079492 if dataset == 'd2' else 177798571
kept = 328059484 if dataset == 'd2' else 170359591
marked = 8894488 if dataset == 'd2' else 50084605
duplicated = marked - 1 if os.environ.get('BAD_DUPLICATES') == label else marked
records = expected - 1 if os.environ.get('BAD_RECORDS') == label else expected
mode = 'external' if os.environ.get('BAD_MODE') == label else 'memory'
lines = [f'{prefix} ranks={ranks} total_records={records} kept_records={kept} core=4.0',
         'Complete the 4.3 operator cost 4.000000',
         f'selected mode={mode}', f'  duplicates marked={marked}',
         'SN\tis sorted:\t1', f'SN\treads duplicated:\t{duplicated}\t# mock']
for rank in range(ranks):
    index = 0 if os.environ.get('DUPLICATE_RANKS') == label else rank
    processing = 8 + rank
    if os.environ.get('PROCESSING_MAX_OTHER_RANK') == label:
        processing = 9.5 if rank == 0 else 6.0
    if os.environ.get('BAD_PROCESSING') == label:
        processing = 999.0
    if os.environ.get('OMIT_PROCESSING') != label:
        processing_rank = 0 if os.environ.get('DUPLICATE_PROCESSING') == label else rank
        lines += [f'[swbam-processing rank={processing_rank}] storage=ram seconds={processing} status=0']
    lines += [f'[swbam-e2e rank={index}] storage=ram seconds={10 + rank}.0 status=0',
              f'[ram-prepare rank={rank}] seconds={100 + rank}.0',
              f'[ram-storage rank={rank}] peak={1000 + rank} output_bytes=2000 temp_read_bytes=0 '
              f'temp_write_bytes={1 if os.environ.get("BAD_TEMP") == label else 0}']
if os.environ.get('FAIL_LABEL') == label:
    lines += ['ERROR: mock capacity failure']
log.write_text('\n'.join(lines) + '\n')
print('Job <100> has been submitted to queue <q_share>')
'''


class MatrixTest(unittest.TestCase):
    def run_matrix(self, overrides, dataset='D1'):
        overrides = dict(overrides)
        recover = overrides.pop('RECOVER', False)
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            root = base / 'matrix'
            (root / 'bin').mkdir(parents=True)
            for name in ('RabbitBAM-MPI', 'swbam-sdk-optimized-count',
                         'swbam-sdk-optimized-bam2bam'):
                path = root / 'bin' / name
                path.write_text('#!/bin/sh\nexit 0\n')
                path.chmod(0o700)
            scheduler = base / 'scheduler'
            scheduler.mkdir()
            programs = {
                'bsub': MOCK_SUBMIT,
                'bjobs': '#!/bin/sh\nprintf "100 %s\\n" "${MOCK_STATUS:-DONE}"\n',
                'sleep': '#!/bin/sh\nexit 0\n',
            }
            for name, contents in programs.items():
                path = scheduler / name
                path.write_text(contents)
                path.chmod(0o700)
            inputs = []
            for name in ('input.bam', 'input.sam', 'name.bam', 'fixmate.bam', 'sorted.bam'):
                path = base / name
                path.write_bytes(b'mock input')
                inputs.append(str(path))
            env = dict(os.environ, PATH=str(scheduler) + ':' + os.environ['PATH'])
            env.update(overrides)
            if env.get('BENCH_MARKDUP_APP') == 'separate':
                markdup = root / 'bin' / 'RabbitBAM-MPI-markdup'
                markdup.write_text('#!/bin/sh\n# New markdup scheduler\nexit 0\n')
                markdup.chmod(0o700)
                env['BENCH_MARKDUP_APP'] = str(markdup)
            if dataset == 'D2':
                baseline = base / 'prior_probe'
                (baseline / 'logs').mkdir(parents=True)
                hashes = subprocess.run(['sha256sum', 'bin/RabbitBAM-MPI',
                                         'bin/swbam-sdk-optimized-count', 'bin/swbam-sdk-optimized-bam2bam'],
                                        cwd=root, check=True, stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE, universal_newlines=True).stdout
                if env.get('BAD_BASELINE_HASH'):
                    hashes = hashes.replace(hashes[:64], '0' * 64)
                (baseline / 'binaries.sha256').write_text(hashes)
                stats = subprocess.run(['stat', '-c', '%n\\t%s bytes'] + inputs,
                                       check=True, stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE, universal_newlines=True).stdout
                (baseline / 'configuration.txt').write_text('ranks=6\n' + stats)
                names = ('count', 'sdk_copy', 'sdk_validate', 'flagstat', 'filter', 'collate', 'fixmate')
                for name in names:
                    subprocess.run([str(scheduler / 'bsub'), '-J', 'd2ram_' + name,
                                    '-o', str(baseline / 'logs' / (name + '.log')), '--ram-ranks', '6'],
                                   env=env, check=True, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, universal_newlines=True)
                (baseline / 'results.tsv').write_text('name\tjob_id\n' + ''.join(name + '\t100\n' for name in names))
                env['BENCH_D2_BASELINE_ROOT'] = str(baseline)
                if env.get('PRIOR_MATRIX_STATE'):
                    prior = base / 'prior_matrix'
                    prior.mkdir()
                    (prior / 'state').write_text(env['PRIOR_MATRIX_STATE'] + '\n')
                    env['BENCH_WAIT_FOR'] = str(prior)
            script = Path(__file__).resolve().with_name('ram_' + dataset.lower() + '_matrix.sh')
            command = ['bash', str(script), str(root)] + inputs
            result = subprocess.run(command, env=env, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE, universal_newlines=True, timeout=30)
            if recover:
                self.assertNotEqual(result.returncode, 0)
                label = 'probe_sort_n6_r1' if dataset == 'D2' else 'sdk_count_n1_r1'
                log = root / 'logs' / (label + '.log')
                log.write_text(log.read_text().replace('ERROR: mock capacity failure\n', ''))
                env.pop('FAIL_LABEL')
                if env.pop('CHANGE_MARKDUP_BINARY', False):
                    markdup.write_text('#!/bin/sh\n# Changed after freezing\nexit 0\n')
                env['BENCH_RESUME'] = '1'
                result = subprocess.run(command, env=env, stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE, universal_newlines=True, timeout=30)
            files = {name: (root / name).read_text() if (root / name).exists() else ''
                     for name in ('state', 'commands.sh', 'results.tsv', 'workflows.tsv', 'workflow_summary.txt',
                                  'capacity_gate.txt', 'prior_probes.tsv', 'markdup.sha256')}
            files['diagnostic'] = result.stdout + result.stderr
            return result.returncode, files

    def test_complete_scope_and_workflow_totals(self):
        status, files = self.run_matrix({})
        self.assertEqual(status, 0, files['diagnostic'])
        self.assertEqual(files['state'], 'COMPLETE\n')
        rows = [line.split('\t') for line in files['results.tsv'].splitlines()[1:]]
        self.assertEqual(len(rows), 74)
        for kind, count in (('sdk', 30), ('app', 24), ('workflow', 12), ('validation', 8)):
            self.assertEqual(sum(row[0] == kind for row in rows), count)
        self.assertEqual(files['commands.sh'].count('--no-dump'), 34)
        self.assertEqual(files['workflow_summary.txt'],
                         'workflow_e2e_median_s=60.000000\nworkflow_processing_median_s=52.000000\n')
        workflow_rows = [line.split('\t') for line in files['workflows.tsv'].splitlines()[1:]]
        self.assertEqual(len(workflow_rows), 3)
        self.assertTrue(all(float(row[5]) == 60.0 for row in workflow_rows))
        self.assertTrue(all(float(row[-1]) == 52.0 for row in workflow_rows))
        self.assertTrue(all(len(row) == 12 for row in rows))
        self.assertIn('processing_max_s', files['results.tsv'].splitlines()[0])
        for repeat, points in (('1', ['1', '2', '4', '6', '12']),
                               ('2', ['12', '6', '4', '2', '1']),
                               ('3', ['1', '2', '4', '6', '12'])):
            self.assertEqual([row[3] for row in rows if row[0:3] == ['sdk', 'count', repeat]], points)
        self.assertIn('sdk\tcount\t1\t1\t100\t12G\t10.000000\t4.000000', files['results.tsv'])

    def test_bad_record_count_stops(self):
        status, files = self.run_matrix({'BAD_RECORDS': 'sdk_count_n1_r1'})
        self.assertNotEqual(status, 0)
        self.assertIn('FAILED sdk_count_n1_r1 job=100', files['state'])
        self.assertEqual(len(files['commands.sh'].splitlines()), 1)
        self.assertEqual(len(files['results.tsv'].splitlines()), 1)

    def test_capacity_error_stops(self):
        status, files = self.run_matrix({'FAIL_LABEL': 'sdk_count_n1_r1'})
        self.assertNotEqual(status, 0)
        self.assertEqual(len(files['commands.sh'].splitlines()), 1)

    def test_scheduler_exit_stops(self):
        status, files = self.run_matrix({'MOCK_STATUS': 'EXIT'})
        self.assertNotEqual(status, 0)
        self.assertEqual(len(files['commands.sh'].splitlines()), 1)

    def test_duplicate_rank_ids_rejected(self):
        status, files = self.run_matrix({'DUPLICATE_RANKS': 'sdk_count_n2_r1'})
        self.assertNotEqual(status, 0)
        self.assertIn('FAILED sdk_count_n2_r1', files['state'])

    def test_processing_max_is_independent_of_e2e_max(self):
        status, files = self.run_matrix({'PROCESSING_MAX_OTHER_RANK': 'sdk_count_n2_r1'})
        self.assertEqual(status, 0, files['diagnostic'])
        rows = [line.split('\t') for line in files['results.tsv'].splitlines()[1:]]
        sample = next(row for row in rows if row[:4] == ['sdk', 'count', '1', '2'])
        self.assertEqual(float(sample[6]), 11.0)
        self.assertEqual(float(sample[11]), 9.5)

    def test_missing_processing_timer_rejected(self):
        status, files = self.run_matrix({'OMIT_PROCESSING': 'sdk_count_n1_r1'})
        self.assertNotEqual(status, 0)
        self.assertIn('FAILED sdk_count_n1_r1', files['state'])

    def test_duplicate_processing_rank_rejected(self):
        status, files = self.run_matrix({'DUPLICATE_PROCESSING': 'sdk_count_n2_r1'})
        self.assertNotEqual(status, 0)
        self.assertIn('FAILED sdk_count_n2_r1', files['state'])

    def test_processing_longer_than_e2e_rejected(self):
        status, files = self.run_matrix({'BAD_PROCESSING': 'sdk_count_n1_r1'})
        self.assertNotEqual(status, 0)
        self.assertIn('FAILED sdk_count_n1_r1', files['state'])

    def test_external_sort_not_imported(self):
        status, files = self.run_matrix({'BAD_MODE': 'app_sort_n6_r1'})
        self.assertNotEqual(status, 0)
        self.assertIn('FAILED app_sort_n6_r1', files['state'])

    def test_resume_without_duplicate_submission(self):
        status, files = self.run_matrix({'FAIL_LABEL': 'sdk_count_n1_r1', 'RECOVER': True})
        self.assertEqual(status, 0)
        self.assertEqual(len(files['commands.sh'].splitlines()), 74)
        self.assertEqual(len(files['results.tsv'].splitlines()), 75)

    def test_markdup_and_workflow_only_scope(self):
        status, files = self.run_matrix({'BENCH_STAGE': 'markdup-workflow',
                                         'BENCH_MARKDUP_APP': 'separate'})
        self.assertEqual(status, 0, files['diagnostic'])
        rows = [line.split('\t') for line in files['results.tsv'].splitlines()[1:]]
        self.assertEqual(len(rows), 21)
        self.assertEqual(sum(row[0] == 'app' for row in rows), 3)
        self.assertTrue(all(row[1] == 'markdup' for row in rows if row[0] == 'app'))
        self.assertEqual(sum(row[0] == 'workflow' for row in rows), 12)
        self.assertEqual(sum(row[0] == 'validation' for row in rows), 6)
        commands = files['commands.sh'].splitlines()
        self.assertEqual(sum('RabbitBAM-MPI-markdup markdup' in line for line in commands), 6)
        self.assertEqual(sum('RabbitBAM-MPI collate' in line for line in commands), 3)
        self.assertEqual(files['state'], 'COMPLETE\n')
        self.assertIn('RabbitBAM-MPI-markdup', files['markdup.sha256'])
        self.assertIn('workflow_processing_median_s=52.000000', files['workflow_summary.txt'])

    def test_markdup_version_change_rejected_on_resume(self):
        status, files = self.run_matrix({'BENCH_MARKDUP_APP': 'separate',
                                         'FAIL_LABEL': 'sdk_count_n1_r1',
                                         'RECOVER': True, 'CHANGE_MARKDUP_BINARY': '1'})
        self.assertNotEqual(status, 0)
        self.assertEqual(len(files['commands.sh'].splitlines()), 1)
        self.assertIn('FAILED', files['diagnostic'])


if __name__ == '__main__':
    unittest.main()
