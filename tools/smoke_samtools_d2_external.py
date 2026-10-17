#!/usr/bin/env python3
"""Run the external driver on a generated paired-read fixture using Samtools."""

import argparse
import csv
import json
import os
from pathlib import Path
import subprocess

import bench_samtools_d2_external as driver


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', required=True)
    parser.add_argument('--samtools', required=True)
    parser.add_argument('--libdeflate', required=True)
    args = parser.parse_args()
    root = Path(args.directory)
    root.mkdir(parents=True, exist_ok=True)
    sam, bam = root/'input.sam', root/'input.bam'
    pairs = 100000
    sequence, quality = 'ACGT'*37+'AC', 'I'*150
    with open(str(sam), 'w') as output:
        output.write('@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:2000000\n')
        for i in range(pairs):
            position = i*10+1
            output.write('p{0}\t99\tchr1\t{1}\t60\t150M\t=\t{2}\t152\t{3}\t{4}\n'.format(
                i, position, position+2, sequence, quality))
            output.write('p{0}\t147\tchr1\t{1}\t60\t150M\t=\t{2}\t-152\t{3}\t{4}\n'.format(
                i, position+2, position, sequence, quality))
    env = dict(os.environ, LD_LIBRARY_PATH=args.libdeflate+':'+os.environ.get('LD_LIBRARY_PATH',''))
    subprocess.run([args.samtools,'view','--no-PG','-b','-o',str(bam),str(sam)], env=env, check=True)
    config = argparse.Namespace(root=str(root/'run'), input=str(bam), samtools=args.samtools,
                                libdeflate=args.libdeflate, records=2*pairs, duplicates=0,
                                input_bytes=bam.stat().st_size, sort_mib=24, candidates=[24],
                                repeats=1, resource_wait_seconds=30, temp_sample_interval=0.02)
    bench = driver.ExternalBenchmark(config)
    try:
        bench.run()
    finally:
        bench.writer_file.close()
        bench.lock.close()
    with open(str(root/'run/results.tsv')) as source:
        rows = list(csv.DictReader(source, delimiter='\t'))
    assert len(rows)==9, len(rows)
    assert all(int(row['sort_runs'])>0 for row in rows if row['app']=='sort')
    assert all(int(row['temp_files_peak'])>0 for row in rows if row['app']=='collate')
    assert (root/'run/ALL_DONE').exists()
    assert len(list(csv.DictReader((root/'run/workflows.tsv').read_text().splitlines(), delimiter='\t')))==1
    workflows = [row for row in rows if row['phase']=='workflow']
    assert len(workflows)==4
    summary = json.loads((root/'run/summary.json').read_text())
    assert abs(summary['workflow_median_s']-sum(float(row['wall_s']) for row in workflows))<0.000001
    assert int(workflows[-1]['records_checked'])==2*pairs
    assert int(workflows[-1]['duplicates_checked'])==0
    print('PASS: external sort, real collate temps, full workflow, {} records'.format(2*pairs))


if __name__=='__main__':
    main()
