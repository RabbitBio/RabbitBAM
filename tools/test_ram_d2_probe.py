#!/usr/bin/env python3
"""Test probe orchestration with a mock scheduler, never submit real jobs."""

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
log = Path(args[args.index('-o') + 1])
name = args[args.index('-J') + 1][len('ram_D2_'):]
if name == 'sdk_copy':
    second_output = args.index('-o', args.index('-o') + 1)
    Path(args[second_output + 1]).write_bytes(b'mock output')
ranks = int(args[args.index('-N') + 1])
if os.environ.get('BAD_RANKS') == name:
    ranks -= 1
lines = ['total_records=365079492 kept_records=328059484 marked=8894488',
         'selected mode=memory', 'SN\tis sorted:\t1']
for rank in range(ranks):
    if os.environ.get('OMIT_PROCESSING') != name:
        processing = 999.0 if os.environ.get('BAD_PROCESSING') == name else 8 + rank
        processing_rank = 0 if os.environ.get('DUPLICATE_PROCESSING') == name else rank
        lines += [f'[swbam-processing rank={processing_rank}] seconds={processing} status=0']
    lines += [f'[swbam-e2e rank={rank}] seconds={10 + rank}.0 status=0',
              f'[ram-prepare rank={rank}] seconds={100 + rank}.0',
              f'[ram-storage rank={rank}] peak={1000 + rank} output_bytes=2000']
if os.environ.get('FAIL_NAME') == name:
    lines += ['ERROR: mock capacity failure']
if not os.environ.get('OMIT_FOOTER'):
    lines += ['Job 100 has been finished.']
log.write_text('\n'.join(lines) + '\n')
print('Job <100> has been submitted to queue <q_share>')
'''


class ProbeTest(unittest.TestCase):
    def run_probe(self, overrides):
        overrides = dict(overrides)
        recover = overrides.pop('RECOVER', False)
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            root = base / 'probe'
            (root / 'bin').mkdir(parents=True)
            for name in ('RabbitBAM-MPI', 'swbam-sdk-optimized-count',
                         'swbam-sdk-optimized-bam2bam'):
                binary = root / 'bin' / name
                binary.write_text('#!/bin/sh\nexit 0\n')
                binary.chmod(0o700)
            scheduler = base / 'scheduler'
            scheduler.mkdir()
            bsub = scheduler / 'bsub'
            bsub.write_text(MOCK_SUBMIT)
            bsub.chmod(0o700)
            bjobs = scheduler / 'bjobs'
            bjobs.write_text('#!/bin/sh\nprintf "100 %s\\n" "${MOCK_STATUS:-DONE}"\n')
            bjobs.chmod(0o700)
            sleep = scheduler / 'sleep'
            sleep.write_text('#!/bin/sh\nexit 0\n')
            sleep.chmod(0o700)
            inputs = []
            for name in ('input.bam', 'input.sam', 'name.bam',
                         'fixmate.bam', 'sorted.bam'):
                path = base / name
                path.write_bytes(b'mock input')
                inputs.append(str(path))
            env = dict(os.environ, PATH=str(scheduler) + ':' + os.environ['PATH'])
            env.update(overrides)
            script = Path(__file__).resolve().with_name('ram_d2_probe.sh')
            result = subprocess.run(['bash', str(script), str(root)] + inputs,
                                    env=env, capture_output=True, text=True,
                                    timeout=15)
            if recover:
                self.assertNotEqual(result.returncode, 0)
                log = root / 'logs' / 'count.log'
                log.write_text(log.read_text().replace('ERROR: mock capacity failure\n', ''))
                env.pop('FAIL_NAME', None)
                env['PROBE_RESUME'] = '1'
                result = subprocess.run(['bash', str(script), str(root)] + inputs,
                                        env=env, capture_output=True, text=True,
                                        timeout=15)
            return (result.returncode, (root / 'state').read_text(),
                    (root / 'results.tsv').read_text(),
                    (root / 'commands.sh').read_text())

    def test_complete_and_metrics(self):
        status, state, rows, commands = self.run_probe({})
        self.assertEqual(status, 0)
        self.assertEqual(state, 'COMPLETE\n')
        self.assertEqual(len(rows.splitlines()), 12)
        self.assertIn('count\t100\t8G\t15.000000\t105.000000\t1005\t2000', rows)
        self.assertIn('processing_max_s', rows.splitlines()[0])
        self.assertTrue(all(len(row.split('\t')) == 8 for row in rows.splitlines()))
        self.assertEqual(rows.splitlines()[1].split('\t')[-1], '13.000000')
        self.assertEqual(commands.count('--ram-ranks 6'), 11)
        self.assertEqual(commands.count('--no-dump'), 7)

    def test_capacity_error_stops(self):
        status, state, rows, commands = self.run_probe({'FAIL_NAME': 'count'})
        self.assertNotEqual(status, 0)
        self.assertIn('FAILED count job=100', state)
        self.assertEqual(len(rows.splitlines()), 1)
        self.assertEqual(commands.count('--ram-ranks 6'), 1)

    def test_single_rank_capacity_gates(self):
        status, state, rows, commands = self.run_probe({'PROBE_SDK_SINGLE_RANK': '1'})
        self.assertEqual(status, 0)
        self.assertEqual(state, 'COMPLETE\n')
        self.assertEqual(len(rows.splitlines()), 14)
        self.assertIn('count_n1\t100\t36G\t10.000000', rows)
        self.assertIn('sdk_copy_n1\t100\t76G\t10.000000', rows)
        self.assertEqual(commands.count('--ram-ranks 1 '), 2)
        self.assertEqual(commands.count('--ram-ranks 6'), 11)
        self.assertIn('--ram-capacity 16G', commands)

    def test_done_without_footer(self):
        status, state, rows, commands = self.run_probe({'OMIT_FOOTER': '1'})
        self.assertEqual(status, 0)
        self.assertEqual(state, 'COMPLETE\n')
        self.assertEqual(len(rows.splitlines()), 12)

    def test_scheduler_exit_rejects_success_log(self):
        status, state, rows, commands = self.run_probe({'MOCK_STATUS': 'EXIT'})
        self.assertNotEqual(status, 0)
        self.assertIn('FAILED count job=100', state)
        self.assertEqual(len(rows.splitlines()), 1)

    def test_resume_reuses_submitted_job(self):
        status, state, rows, commands = self.run_probe({'FAIL_NAME': 'count', 'RECOVER': True})
        self.assertEqual(status, 0)
        self.assertEqual(state, 'COMPLETE\n')
        self.assertEqual(len(rows.splitlines()), 12)
        self.assertEqual(commands.count('ram_D2_count '), 1)
        self.assertEqual(commands.count('--ram-ranks 6'), 11)

    def test_missing_rank_stops(self):
        status, state, rows, commands = self.run_probe({'BAD_RANKS': 'count'})
        self.assertNotEqual(status, 0)
        self.assertIn('FAILED count job=100', state)
        self.assertEqual(len(rows.splitlines()), 1)

    def test_invalid_processing_timer_stops(self):
        for flag in ('OMIT_PROCESSING', 'BAD_PROCESSING', 'DUPLICATE_PROCESSING'):
            with self.subTest(flag=flag):
                status, state, rows, commands = self.run_probe({flag: 'count'})
                self.assertNotEqual(status, 0)
                self.assertIn('FAILED count job=100', state)
                self.assertEqual(len(rows.splitlines()), 1)


if __name__ == '__main__':
    unittest.main()
