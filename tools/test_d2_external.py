#!/usr/bin/env python3
"""Validate external-mode gates, budgets, and real workflow aggregation."""

import argparse
import csv
import json
import os
from pathlib import Path
import shlex
import subprocess
import tempfile
import unittest

import bench_samtools_d2_external as a4


MOCK_SUBMIT = r'''#!/usr/bin/env python3
import os
from pathlib import Path
import sys
args = sys.argv[1:]
label = args[args.index('-J')+1][len('d2ram_external_'):]
if os.environ.get('SUBMISSION_FAILURE') == label:
    sys.exit('mock scheduler submission failure')
log = Path(args[args.index('-o')+1])
outputs = [i for i, arg in enumerate(args) if arg == '-o']
if len(outputs)>1 and '--no-dump' not in args:
    Path(args[outputs[1]+1]).write_bytes(b'complete output')
external = any(label.startswith(kind+'_'+app) for kind in ('probe','app','workflow') for app in ('collate','sort'))
mode = 'external' if external and os.environ.get('BAD_MODE') != label else 'memory'
records = 365079491 if os.environ.get('BAD_RECORDS') == label else 365079492
lines = [f'OptimizedMockMPI finished. ranks=6 records={records}',
         'Complete the 4.3 operator cost 4.0', f'selected mode={mode}',
         'external_actual_sum runs=18', 'external_memory_sum tracked_peak_max=8589934592',
         '  duplicates marked=8894488', 'SN\tis sorted:\t1',
         'SN\treads duplicated:\t8894488']
for rank in range(6):
    spill = 100 if external and os.environ.get('NO_TEMP') != label else 0
    lines += [f'[swbam-processing rank={rank}] storage=ram seconds=13 status=0',
              f'[swbam-e2e rank={rank}] storage=ram seconds=20 status=0',
              f'[ram-prepare rank={rank}] seconds=100',
              f'[ram-storage rank={rank}] peak=1000 output_bytes=2000 temp_read_bytes={spill} temp_write_bytes={spill}']
if os.environ.get('CAPACITY_FAILURE') == label:
    lines.append('ERROR: RAM pool capacity exhausted')
log.write_text('\n'.join(lines)+'\n')
print('Job <100> has been submitted to queue <q_share>')
'''


class ExternalTest(unittest.TestCase):
    def run_swls(self, overrides, resume=False):
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            root = base / 'run'
            (root / 'bin').mkdir(parents=True)
            for name in ('RabbitBAM-MPI', 'RabbitBAM-MPI-markdup'):
                path = root / 'bin' / name
                path.write_text('#!/bin/sh\nexit 0\n')
                path.chmod(0o700)
            scheduler = base / 'scheduler'
            scheduler.mkdir()
            for name, contents in (
                    ('bsub', MOCK_SUBMIT),
                    ('bjobs', '#!/bin/sh\nprintf "100 DONE\\n"\n'),
                    ('sleep', '#!/bin/sh\nexit 0\n')):
                path = scheduler / name
                path.write_text(contents)
                path.chmod(0o700)
            inputs = []
            for name in ('input.bam', 'fixmate.bam'):
                path = base / name
                path.write_bytes(b'mock input')
                inputs.append(str(path))
            env = dict(os.environ, PATH=str(scheduler)+':'+os.environ['PATH'])
            env.update(overrides)
            script = Path(__file__).with_name('ram_d2_external.sh')
            command = ['bash', str(script), str(root)] + inputs
            result = subprocess.run(command,
                                    env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    universal_newlines=True, timeout=30)
            before_resume = ''
            if resume:
                self.assertNotEqual(result.returncode, 0)
                before_resume = (root / 'results.tsv').read_text()
                env.pop('SUBMISSION_FAILURE', None)
                env['BENCH_RESUME'] = '1'
                result = subprocess.run(command, env=env, stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE, universal_newlines=True, timeout=30)
            files = {name: (root/name).read_text() if (root/name).exists() else ''
                     for name in ('state','commands.sh','results.tsv','workflows.tsv',
                                  'workflow_summary.txt','capacity_gate.txt')}
            files['before_resume_results.tsv'] = before_resume
            return result, files

    def test_complete_matrix(self):
        result, files = self.run_swls({})
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        rows = list(csv.DictReader(files['results.tsv'].splitlines(), delimiter='\t'))
        self.assertEqual(len(rows), 32)
        self.assertEqual(sum(row['kind']=='app' for row in rows), 6)
        self.assertEqual(sum(row['kind']=='workflow' for row in rows), 12)
        self.assertEqual(files['state'], 'COMPLETE\n')
        self.assertIn('GATES PASSED', files['capacity_gate.txt'])
        self.assertEqual(files['workflow_summary.txt'],
                         'workflow_e2e_median_s=80.000000\nworkflow_processing_median_s=52.000000\n')
        self.assertTrue(all(row['capacity']=='64G' and row['temp_read_bytes']=='600'
                            for row in rows if row['app'] in ('sort','collate')
                            and row['kind'] in ('app','probe','workflow')))
        self.assertEqual(next(row['capacity'] for row in rows
                              if row['kind']=='validation' and row['app']=='collate'), '12G')
        self.assertTrue(all(row['capacity']=='8G' for row in rows
                            if row['kind'] in ('validation', 'probe_validation')
                            and row['app'] != 'collate'))
        commands = files['commands.sh'].splitlines()
        self.assertEqual(len(commands), 32)
        self.assertTrue(all('--no-dump' not in line for line in commands if 'workflow_' in line))
        for app in ('collate', 'sort'):
            for repeat in (1, 2, 3):
                label = 'd2ram_external_app_{}_r{}'.format(app, repeat)
                command = next(shlex.split(line) for line in commands if label in line)
                self.assertEqual('-n' in command, app == 'collate')
                if app == 'collate':
                    self.assertEqual(command[command.index('-n')+1], '66')
                self.assertEqual('--no-dump' in command, repeat != 1)
                self.assertEqual(command[command.index('-m')+1], '8G')

    def test_capacity_failure_stops_before_formal_samples(self):
        result, files = self.run_swls({'CAPACITY_FAILURE':'probe_collate'})
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(len(files['commands.sh'].splitlines()), 1)
        self.assertFalse(files['capacity_gate.txt'])

    def test_resume_reuses_all_six_passed_probes(self):
        result, files = self.run_swls({'SUBMISSION_FAILURE': 'app_collate_r1'}, resume=True)
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        prior = files['before_resume_results.tsv']
        self.assertEqual(len(list(csv.DictReader(prior.splitlines(), delimiter='\t'))), 6)
        self.assertTrue(files['results.tsv'].startswith(prior))
        self.assertEqual(len(files['results.tsv'].splitlines()), 33)
        self.assertEqual(result.stdout.count('REUSE verified probe'), 6)
        commands = files['commands.sh'].splitlines()
        self.assertEqual(len(commands), 33)  # Includes the failed submission attempt.
        self.assertEqual(sum('-J d2ram_external_probe_' in line and
                             '-J d2ram_external_probe_validation_' not in line
                             for line in commands), 4)
        self.assertEqual(sum('-J d2ram_external_probe_validation_' in line for line in commands), 2)
        self.assertEqual(files['state'], 'COMPLETE\n')

    def test_memory_path_rejected(self):
        result, files = self.run_swls({'BAD_MODE':'probe_sort'})
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(len(files['commands.sh'].splitlines()), 3)
        self.assertFalse(files['capacity_gate.txt'])

    def test_resume_preserves_finished_collate_sample(self):
        result, files = self.run_swls({'SUBMISSION_FAILURE': 'validation_collate_r1'}, resume=True)
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        prior = files['before_resume_results.tsv']
        self.assertEqual(len(list(csv.DictReader(prior.splitlines(), delimiter='\t'))), 7)
        self.assertTrue(files['results.tsv'].startswith(prior))
        self.assertEqual(len(files['results.tsv'].splitlines()), 33)
        self.assertIn('REUSE verified app_collate_r1', result.stdout)
        self.assertEqual(sum('-J d2ram_external_app_collate_r1 ' in line
                             for line in files['commands.sh'].splitlines()), 1)
        self.assertEqual(files['state'], 'COMPLETE\n')

    def test_external_without_actual_temp_io_rejected(self):
        result, files = self.run_swls({'NO_TEMP':'probe_collate'})
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(files['capacity_gate.txt'])

    def test_wrong_count_rejected(self):
        result, files = self.run_swls({'BAD_RECORDS':'probe_validation_markdup'})
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(files['capacity_gate.txt'])

    def test_a4_sort_budget_and_mode(self):
        for threads in (24,48,96,127):
            memory = a4.sort_memory(49152, threads)
            self.assertLessEqual(memory*threads, 49152)
            self.assertLess(49152-memory*threads, threads)
        self.assertEqual(a4.sort_runs('[bam_sort_core] merging from 5 files and 24 in-memory blocks...'), 5)
        for message in ('', 'merging from 0 files and 24 in-memory blocks...'):
            with self.assertRaises(RuntimeError):
                a4.sort_runs(message)
        bench = a4.ExternalBenchmark.__new__(a4.ExternalBenchmark)
        bench.args = argparse.Namespace(samtools='/bin/samtools', input='/data/input.bam', sort_mib=49152)
        command = bench.command('D2','sort',24,Path('/out.bam'),Path('/tmp/run'))
        self.assertEqual(command[command.index('-m')+1], '2048M')
        command = bench.command('D2','collate',127,Path('/out.bam'),Path('/tmp/run'))
        self.assertEqual(command[command.index('-n')+1], '66')
        self.assertNotIn('-m', command)


if __name__ == '__main__':
    unittest.main()
