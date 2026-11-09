#!/usr/bin/env python3
"""Create fresh SAMtools references and frozen per-stage common inputs on A4."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import time


DATA = {
    'D1': ('D1-WES-9G', 177798571, 50084605, 'b7bc9220bc5e7051b51b758e7018907644b8ba323d67dec82e63acd2c80bb046'),
    'D2': ('D2-WGS-25G', 365079492, 8894488, '1cbe4e17224d89ef331aeaa655835a0ec467ac94cf89511ddb07e6e9c5d041a6'),
    'D3': ('D3-mouse-SRR7768905', 165373770, 33171753, '472593083f7e98acb40a84ec917341adb9578a64a89bf5b0ecc3c1669a83dc8c'),
}


def sha256(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 << 20), b''):
            result.update(block)
    return result.hexdigest()


def write(path, value):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def shell(argv):
    return ' '.join(shlex.quote(str(value)) for value in argv)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--samtools', default='/home/user_home/wzs/software/samtools-1.20/bin/samtools')
    parser.add_argument('--data', type=Path, default=Path('/home/bigssd/wzs_data'))
    parser.add_argument('--digest', type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    root.mkdir(parents=True, exist_ok=True)
    if (root/'reference_manifest.json').exists() or (root/'commands.jsonl').exists():
        raise RuntimeError('refusing to overwrite an existing reference run')
    (root/'logs').mkdir(exist_ok=True)
    manifest = dict(samtools=args.samtools, samtools_sha256=sha256(args.samtools),
                    version=subprocess.check_output([args.samtools,'--version'], universal_newlines=True),
                    digest_sha256=sha256(args.digest), cases={})

    def state(phase, **extra):
        value = dict(phase=phase, time=time.strftime('%F %T %z'), **extra)
        write(root/'state.json', value)
        print(json.dumps(value), flush=True)

    def run(label, argv, output=None):
        with (root/'commands.jsonl').open('a') as stream:
            stream.write(json.dumps(dict(label=label, argv=argv))+'\n')
        with (root/'logs'/(label+'.stderr')).open('w') as err:
            if output is None:
                subprocess.run(argv, check=True, stderr=err)
            else:
                with output.open('w') as out:
                    subprocess.run(argv, check=True, stdout=out, stderr=err)

    def digest(label, path, grouped=False):
        group_path = root/(label+'.groups')
        pipeline = [args.samtools,'view','--no-PG','-h','-@','8',str(path)]
        checker = [str(args.digest)] + (['--groups',str(group_path)] if grouped else [])
        command = shell(pipeline)+' | '+shell(checker)
        output = root/(label+'.digest.json')
        run(label+'.digest', ['bash','-o','pipefail','-c',command], output)
        result = json.loads(output.read_text())
        if grouped:
            temporary = root/(label+'.sort-tmp')
            temporary.mkdir()
            unique = root/(label+'.unique.txt')
            command = shell(['sort','-u','-S','4G','--parallel=8','-T',str(temporary),str(group_path)])+' | wc -l'
            run(label+'.groups', ['env','LC_ALL=C','bash','-o','pipefail','-c',command], unique)
            result['unique_qnames'] = int(unique.read_text().strip())
            result['continuity_errors'] = result['qname_groups']-result['unique_qnames']
            if result['continuity_errors'] or result['read_order_errors']:
                raise RuntimeError(label+': invalid QNAME grouping/read order')
            group_path.unlink()
            temporary.rmdir()
            write(output,result)
        return result

    try:
        state('INITIALIZING')
        for case,(dataset,records,marked,expected_sha) in DATA.items():
            folder = root/case
            folder.mkdir()
            source = args.data/dataset/'initial'/(dataset+'.coordinate.bam')
            state('HASHING_INPUT',case=case,path=str(source))
            if sha256(source) != expected_sha:
                raise RuntimeError(case+': initial BAM differs from the frozen experiment input')
            original_digest = digest(case+'.input',source)
            if original_digest['records'] != records:
                raise RuntimeError('initial record count mismatch')
            item = dict(input=str(source),input_sha256=expected_sha,input_digest=original_digest, stages={})
            manifest['cases'][case] = item
            previous, previous_sha = source, expected_sha
            for app in ('collate','fixmate','sort','markdup'):
                state('REFERENCE',case=case,app=app)
                available = {line.split(':')[0]:int(line.split()[1])*1024
                             for line in Path('/proc/meminfo').read_text().splitlines() if line.startswith('MemAvailable:')}
                space = os.statvfs(root)
                if space.f_bavail*space.f_frsize < 160*1024**3:
                    raise RuntimeError('less than 160 GiB free; no automatic directory relocation')
                if app == 'sort' and case != 'D2' and available['MemAvailable'] < 128*1024**3:
                    raise RuntimeError('in-memory reference sort requires at least 128 GiB available RAM')
                target = folder/(app+'.bam')
                s = args.samtools
                if app == 'collate':
                    argv=[s,app,'--no-PG','-@','24','-n','66','-l','1','-T',str(folder/'tmp.collate'),'-o',str(target),str(previous)]
                elif app == 'fixmate':
                    argv=[s,app,'--no-PG','-@','24','-m','--output-fmt','BAM,level=1',str(previous),str(target)]
                elif app == 'sort':
                    threads,memory=('96','512M') if case=='D2' else ('24','4G')
                    argv=[s,app,'--no-PG','-@',threads,'-m',memory,'-l','1','-T',str(folder/'tmp.sort'),'-o',str(target),str(previous)]
                else:
                    argv=[s,app,'--no-PG','-@','24','-l','300','--output-fmt','BAM,level=1',str(previous),str(target)]
                label=case+'.'+app
                run(label,argv)
                run(label+'.quickcheck',[s,'quickcheck','-v',str(target)])
                result=digest(label,target,app=='collate')
                if result['records'] != records or (app in ('sort','markdup') and result['coordinate_errors']):
                    raise RuntimeError(label+': record count or coordinate order mismatch')
                if app=='markdup' and result['duplicates'] != marked:
                    raise RuntimeError(label+': unexpected DUP total')
                if app=='sort':
                    stderr=(root/'logs'/(label+'.stderr')).read_text()
                    match=re.search(r'merging from (\d+) files',stderr)
                    if not match or (int(match[1])==0)!=(case!='D2'):
                        raise RuntimeError(label+': unexpected memory/external sort mode')
                output_sha=sha256(target)
                item['stages'][app]=dict(input=str(previous),input_sha256=previous_sha,
                                        output=str(target),output_sha256=output_sha,bytes=target.stat().st_size,
                                        argv=argv,digest=result)
                write(root/'reference_manifest.json',manifest)
                previous,previous_sha=target,output_sha
                state('REFERENCE_PASSED',case=case,app=app)
            write(folder/'READY.json',item)
        state('COMPLETE',cases=list(DATA),scope='fresh SAMtools four-command references; not SWBAM validation')
    except Exception as error:
        state('FAILED',error=str(error))
        raise


if __name__=='__main__':
    main()
