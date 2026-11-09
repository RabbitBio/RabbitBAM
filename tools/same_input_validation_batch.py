#!/usr/bin/env python3
"""Wait for fresh references and SWLS connectivity, then run each correctness case once."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

import probe_nonhuman_d3 as probe


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--app',type=Path,required=True)
    parser.add_argument('--reference-root',required=True)
    parser.add_argument('--swls',default='swls')
    parser.add_argument('--a4',default='wzs@a4.rabbitbio.org')
    parser.add_argument('--wait-seconds',type=int,default=14400)
    parser.add_argument('--resume',action='store_true')
    args=parser.parse_args()
    root=args.root.resolve()
    root.mkdir(parents=True,exist_ok=True)
    scripts=Path(__file__).resolve().parent
    config=dict(app=str(args.app.resolve()),sha256=hashlib.sha256(args.app.read_bytes()).hexdigest(),
                reference_root=args.reference_root,cases=['D1','D2','D3'],
                validation_script_sha256=hashlib.sha256((scripts/'validate_same_input.py').read_bytes()).hexdigest())
    config_path=root/'batch.json'
    if config_path.exists():
        if not args.resume or json.loads(config_path.read_text())!=config:
            raise RuntimeError('existing batch requires explicit resume with the same program and scripts')
    else: config_path.write_text(json.dumps(config,indent=2)+'\n')

    def state(phase,**extra):
        value=dict(phase=phase,time=time.strftime('%F %T %z'),**extra)
        temporary=root/'state.tmp'
        temporary.write_text(json.dumps(value,indent=2)+'\n')
        temporary.replace(root/'state.json')
        print(json.dumps(value),flush=True)

    try:
        for case in config['cases']:
            child=root.parent/(root.name+'_'+case)
            if (child/'state.json').exists() and json.loads((child/'state.json').read_text())['phase']=='COMPLETE':
                if not args.resume: raise RuntimeError('completed case exists; explicit resume required')
                continue
            deadline=time.monotonic()+args.wait_seconds
            while True:
                try:
                    probe.remote(args.swls,['true'],timeout=20)
                except Exception as error:
                    state('WAITING_SWLS',case=case,error=str(error))
                else:
                    status=json.loads(probe.remote(args.a4,['cat',args.reference_root+'/state.json']))
                    if status['phase']=='FAILED': raise RuntimeError('A4 reference failed: '+str(status))
                    ready=probe.remote(args.a4,['bash','-s'],script='test ! -f '+args.reference_root+'/'+case+'/READY.json || printf ready\n')
                    if ready: break
                    state('WAITING_REFERENCE',case=case,reference_status=status)
                if time.monotonic()>=deadline:
                    raise RuntimeError('prerequisites not ready within the fixed wait limit; no automatic retries/submissions')
                time.sleep(120)
            state('VALIDATING',case=case)
            command=['python3','-u','-B',str(scripts/'validate_same_input.py'),'--root',str(child),
                     '--case',case,'--app',config['app'],'--reference-root',args.reference_root,
                     '--swls',args.swls,'--a4',args.a4]
            if args.resume and (child/'plan.json').exists(): command+=['--resume']
            with (root/(case+'.driver.log')).open('a') as log:
                subprocess.run(command,check=True,stdout=log,stderr=subprocess.STDOUT)
            state('CASE_COMPLETE',case=case)
        state('COMPLETE',cases=config['cases'])
    except Exception as error:
        state('FAILED',error=str(error),note='No automatic re-submission; inspect case state and existing job IDs.')
        raise


if __name__=='__main__': main()
