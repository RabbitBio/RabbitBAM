#!/usr/bin/env python3
"""Validate standalone operations and a real four-command workflow with common BAMs."""
import argparse
import hashlib
import json
from pathlib import Path
import shlex
import subprocess
import time

import probe_nonhuman_d3 as probe
from markdup_optimized_retest import sample_metrics


REFERENCE = '/home/bigssd/wzs_data/samtools_bench/same_input_correctness_20261109'
BUILD = Path(__file__).resolve().parents[1]/'build_bounded_io/validation_same_input_20261109/app/RabbitBAM-MPI'
INITIAL = {'D1':'D1-WES-9G', 'D2':'D2-WGS-25G', 'D3':'D3-mouse-SRR7768905'}


def plan(case, app, original, common, out):
    external = case == 'D2'
    memory = ('8G','8G') if external else (('24G','32G') if case=='D1' else ('48G','48G'))
    pool = ('64G','24G','64G','24G') if external else (('8G','12G','8G','8G') if case=='D1' else ('16G','24G','16G','24G'))
    ordering = 'external' if external else 'memory'
    result=[]
    def add(label, operation, source, target, i, extra):
        result.append(dict(label=label, app=operation, source=source, target=target, capacity=pool[i],
                           kind='validation_operation', repeat=1,
                           check=ordering if operation in ('sort','collate') else ('markdup' if operation=='markdup' else 'write'),
                           argv=[app,operation,'-i',source,'-o',target,'--compress-level','1']+extra))
    collated=out+'/collate.bam'
    add('shared_collate','collate',original,collated,0,['-n','66','-m',memory[0]])
    add('standalone_fixmate','fixmate',common+'/collate.bam',out+'/standalone_fixmate.bam',1,['-m'])
    add('standalone_sort','sort',common+'/fixmate.bam',out+'/standalone_sort.bam',2,['-m',memory[1]])
    add('standalone_markdup','markdup',common+'/sort.bam',out+'/standalone_markdup.bam',3,['--stream','-l','300','-m','8G'])
    fixed,ordered=out+'/workflow_fixmate.bam',out+'/workflow_sort.bam'
    add('workflow_fixmate','fixmate',collated,fixed,1,['-m'])
    add('workflow_sort','sort',fixed,ordered,2,['-m',memory[1]])
    add('workflow_markdup','markdup',ordered,out+'/workflow_markdup.bam',3,['--stream','-l','300','-m','8G'])
    if case=='D3':
        source=probe.SW_DATA+'/swbam_bench/ram_D3_mouse_probe_20261016_023500/out/full/records.sam'
        target=out+'/integrated_sam2bam.bam'
        result.append(dict(label='integrated_sam2bam',app='sam2bam',source=source,target=target,
                           capacity='40G',kind='validation_operation',repeat=1,check='write',
                           argv=[app,'run_all','-i',source,'-o',target,'--compress-level','1']))
    return result


def checksum(host, path):
    return probe.remote(host,['sha256sum',path],timeout=7200).split()[0]


def transfer(source_host, source, target_host, target, expected, root):
    probe.remote(target_host,['mkdir','-p',str(Path(target).parent)])
    present=probe.remote(target_host,['bash','-s'],script='test ! -f '+shlex.quote(target)+' || printf present\n')
    if present:
        if checksum(target_host,target)!=expected:
            raise RuntimeError('existing destination has wrong SHA-256: '+target)
        return
    cache=root/'transfer'/Path(target).name
    cache.parent.mkdir(exist_ok=True)
    def copy(src,dst):
        subprocess.run(['rsync','-a','--partial','-e',' '.join(probe.SSH),src,dst],check=True,timeout=4*3600)
    copy(source_host+':'+source,str(cache))
    copy(str(cache),target_host+':'+target)
    if checksum(target_host,target)!=expected:
        raise RuntimeError('transferred BAM SHA-256 differs: '+target)
    cache.unlink()


class Validation(probe.Probe):
    def state(self, phase, **extra):
        value=dict(phase=phase,case=self.args.case,time=time.strftime('%F %T %z'),**extra)
        temporary=self.root/'state.tmp'
        temporary.write_text(json.dumps(value,indent=2)+'\n')
        temporary.replace(self.root/'state.json')
        print(json.dumps(value),flush=True)

    def digest(self, sample, reference):
        label=sample['label']
        target=self.a4+'/out/'+label+'.bam'
        digest_prefix=self.a4+'/out/'+label
        source_hash=checksum(self.args.swls,sample['target'])
        self.state('TRANSFERRING_OUTPUT',label=label)
        transfer(self.args.swls,sample['target'],self.args.a4,target,source_hash,self.root)
        grouped=sample['app']=='collate'
        command=' '.join(map(shlex.quote,[self.args.samtools,'view','--no-PG','-h','-@','8',target]))
        checker=[self.args.reference_root+'/bin/sam_correctness_digest']
        if grouped: checker+=['--groups',digest_prefix+'.groups']
        script='set -euo pipefail\nexport LC_ALL=C\n'+command+' | '+' '.join(map(shlex.quote,checker))+' > '+shlex.quote(digest_prefix+'.json')+'\n'
        if grouped:
            script+='mkdir -p '+shlex.quote(digest_prefix+'.sort-tmp')+'\n'
            script+='sort -u -S 4G --parallel=8 -T '+shlex.quote(digest_prefix+'.sort-tmp')+' '+shlex.quote(digest_prefix+'.groups')+' | wc -l > '+shlex.quote(digest_prefix+'.unique')+'\n'
        (self.root/(label+'.digest.command.sh')).write_text(script)
        self.state('DIGESTING_OUTPUT',label=label)
        probe.remote(self.args.a4,['bash','-s'],script=script,timeout=4*3600)
        actual=json.loads(probe.remote(self.args.a4,['cat',digest_prefix+'.json']))
        if grouped:
            actual['unique_qnames']=int(probe.remote(self.args.a4,['cat',digest_prefix+'.unique']).strip())
            actual['continuity_errors']=actual['qname_groups']-actual['unique_qnames']
        (self.root/(label+'.digest.json')).write_text(json.dumps(actual,indent=2)+'\n')
        expected=reference['input_digest'] if sample['app']=='sam2bam' else reference['stages'][sample['app']]['digest']
        failures=[]
        for key in ('records','duplicates','header_without_pg_sha256','record_multiset_sha256_sum',
                    'without_dup_multiset_sha256_sum','duplicate_members_sha256_sum'):
            if actual[key]!=expected[key]: failures.append(key)
        if sample['app'] in ('sort','markdup') and actual['coordinate_errors']:
            failures.append('coordinate_errors')
        if grouped and (actual['continuity_errors'] or actual['read_order_errors']):
            failures.append('QNAME grouping/read order')
        if sample['app']=='sam2bam' and actual['ordered_sha256']!=expected['ordered_sha256']:
            failures.append('ordered SAM-to-BAM records')
        comparison=dict(label=label,passed=not failures,failures=failures,sha256=source_hash,
                        actual=actual,expected=expected,
                        ordered_equal=actual['ordered_sha256']==expected['ordered_sha256'],
                        comparison='full canonical SAM record multiset and DUP membership fingerprints')
        (self.root/(label+'.comparison.json')).write_text(json.dumps(comparison,indent=2)+'\n')
        if failures: raise RuntimeError(label+': cross-tool differences: '+', '.join(failures))
        # The SWLS output and all fingerprints remain; discard only this verified A4 copy.
        if not target.startswith(self.a4+'/out/'):
            raise RuntimeError('invalid validation-copy cleanup path')
        probe.remote(self.args.a4,['rm','--',target])
        self.state('SEMANTIC_COMPARISON_PASSED',label=label)

    def execute(self):
        case=self.args.case
        self.state('CHECKING_REFERENCE')
        reference=json.loads(probe.remote(self.args.a4,['cat',self.args.reference_root+'/'+case+'/READY.json']))
        frozen=self.root/'reference.json'
        if frozen.exists() and json.loads(frozen.read_text())!=reference:
            raise RuntimeError('reference manifest changed')
        frozen.write_text(json.dumps(reference,indent=2)+'\n')
        app=self.sw+'/bin/RabbitBAM-MPI'
        original=probe.SW_DATA+'/'+INITIAL[case]+'/initial/'+INITIAL[case]+'.coordinate.bam'
        common=self.sw+'/common'
        samples=plan(case,app,original,common,self.sw+'/out')
        local_sha=hashlib.sha256(self.args.app.read_bytes()).hexdigest()
        scripts=('validate_same_input.py','probe_nonhuman_d3.py','markdup_optimized_retest.py','markdup_halo_metrics.py',
                 'sam_correctness_digest.cpp','samtools_same_input_reference.py')
        config=dict(case=case,binary_sha256=local_sha,reference_root=self.args.reference_root,samples=samples,
                    script_sha256={name:hashlib.sha256((Path(__file__).parent/name).read_bytes()).hexdigest() for name in scripts})
        path=self.root/'plan.json'
        if self.args.resume:
            if json.loads(path.read_text())!=config: raise RuntimeError('frozen plan or binary changed')
        else:
            probe.remote(self.args.swls,['bash','-s'],script='set -eu\ntest ! -e '+shlex.quote(self.sw)+'\nmkdir -p '+' '.join(shlex.quote(self.sw+'/'+name) for name in ('bin','logs','out','common','sessions'))+'\n')
            probe.copy(str(self.args.app),self.args.swls+':'+app)
            path.write_text(json.dumps(config,indent=2)+'\n')
        if checksum(self.args.swls,app)!=local_sha: raise RuntimeError('SWBAM binary changed')
        if checksum(self.args.swls,original)!=reference['input_sha256']:
            raise RuntimeError('initial BAM differs between platforms')
        for operation in ('collate','fixmate','sort'):
            item=reference['stages'][operation]
            self.state('MIRRORING_COMMON_INPUT',app=operation)
            transfer(self.args.a4,item['output'],self.args.swls,common+'/'+operation+'.bam',item['output_sha256'],self.root)
        probe.remote(self.args.a4,['mkdir','-p',self.a4+'/out'])
        spec=dict(records=reference['input_digest']['records'],marked=reference['stages']['markdup']['digest']['duplicates'])
        accepted=[]
        for sample in samples:
            self.run(sample['label'],sample['capacity'],sample['argv'],spec,sample['check'])
            log=(self.root/'logs'/(sample['label']+'.log')).read_text()
            performance=sample_metrics(log,sample,spec)
            # This is a correctness run; elapsed times never replace formal samples.
            actual_source_hash=checksum(self.args.swls,sample['source'])
            independent=sample['label'].startswith('standalone_') or sample['app']=='collate'
            if independent and actual_source_hash!=reference['stages'][sample['app']]['input_sha256']:
                raise RuntimeError('standalone input is not byte-identical: '+sample['label'])
            self.digest(sample,reference)
            accepted.append(dict(label=sample['label'],input_sha256=actual_source_hash,performance=performance))
            (self.root/'accepted.json').write_text(json.dumps(accepted,indent=2)+'\n')
        self.state('COMPLETE',accepted=len(accepted),scope='4 same-input standalone operations + actual 4-command workflow; collate shared; D3 SAM-to-BAM integration gate')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--case',choices=INITIAL,required=True)
    parser.add_argument('--app',type=Path,default=BUILD)
    parser.add_argument('--reference-root',default=REFERENCE)
    parser.add_argument('--samtools',default='/home/user_home/wzs/software/samtools-1.20/bin/samtools')
    parser.add_argument('--swls',default='swls')
    parser.add_argument('--a4',default='wzs@a4.rabbitbio.org')
    parser.add_argument('--resume',action='store_true')
    args=parser.parse_args()
    args.nodes=None; args.retry_no_start=False; args.poll_seconds=60; args.job_prefix='sameinput_'+args.case+'_'
    run=Validation(args)
    try:
        run.execute()
    except Exception as error:
        run.state('FAILED',error=str(error),note='No automatic resubmission; old reference/performance results preserved.')
        raise


if __name__=='__main__':
    main()
