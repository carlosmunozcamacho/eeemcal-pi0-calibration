# SPDX-License-Identifier: LGPL-3.0-or-later
"""Slurm submission only: extraction array followed by one global solve."""
import math
import os
from pathlib import Path
import shlex
import subprocess
import sys
import yaml
from .io import read_json,write_json,digest

SITE=dict(account=None,partition=None,qos=None,launcher=[],python=None,files_per_shard=25,files_per_task=50,concurrency=16,
          extraction=dict(cpus=2,memory='8G',time='04:00:00'),analysis=dict(cpus=2,memory='16G',time='01:00:00'))


def site_config(path):
    import copy
    from .config import merge
    c=copy.deepcopy(SITE)
    if path: merge(c,yaml.safe_load(Path(path).read_text()) or {})
    c['python']=c['python'] or sys.executable
    if not isinstance(c['launcher'],list): raise ValueError('launcher must be an argv list')
    if c['files_per_task']%c['files_per_shard'] or min(c['files_per_task'],c['files_per_shard'],c['concurrency'])<1:
        raise ValueError('files_per_task must be a positive multiple of files_per_shard')
    for name in ['analysis','extraction']:
        if c[name]['cpus']<1: raise ValueError('Stage CPUs must be positive')
    return c


def plan(out,inputs,compact,cache,work,cfg,site):
    out=Path(out);shards=[]
    for i,start in enumerate(range(0,len(inputs),site['files_per_shard'])):
        files=inputs[start:start+site['files_per_shard']]; target=Path(cache)/'shards'/f'shard_{i:06d}'/'candidates.npz'
        shards.append(dict(index=i,inputs=files,output=str(target),signature=digest(dict(inputs=files,config=cfg))))
    per=site['files_per_task']//site['files_per_shard']
    tasks=[list(range(i,min(i+per,len(shards)))) for i in range(0,len(shards),per)]
    p=dict(schema_version=1,output=str(out),cache=str(cache) if cache else None,work=str(work),compact=compact,
           inputs=inputs,config=cfg,site=site,shards=shards,tasks=tasks)
    write_json(out/'batch/work_plan.json',p)
    manifests=out/'batch/task_manifests';manifests.mkdir(exist_ok=True)
    for i,task in enumerate(tasks): (manifests/f'{i:06d}.txt').write_text('\n'.join(f for j in task for f in shards[j]['inputs'])+'\n')
    return p


def complete(shard):
    path=Path(shard['output']);done=path.parent/'complete.json'
    return path.exists() and path.with_suffix('.json').exists() and done.exists() and read_json(done).get('signature')==shard['signature']


def query(submissions,runner=subprocess.run):
    ids=sorted(set(str(x['job_id']) for x in submissions))
    if not ids: return {}
    r=runner(['squeue','-r','-h','-j',','.join(ids),'-o','%i|%T'],capture_output=True,text=True,check=True)
    return dict(line.strip().split('|',1) for line in r.stdout.splitlines() if '|' in line)


def script(p,kind):
    site=p['site'];res=site['extraction' if kind=='extract' else 'analysis']
    command=[site['python'],'-m','eeemcal_pi0','--stage','worker','--plan',str(Path(p['output'])/'batch/work_plan.json'),'--worker-kind',kind]
    args=shlex.join(site['launcher']+command)
    # Array index is deliberately passed through the environment, never interpolated in shell source.
    text='#!/bin/bash\nset -euo pipefail\n'
    text+=f'export OMP_NUM_THREADS={res["cpus"]} OPENBLAS_NUM_THREADS={res["cpus"]} MKL_NUM_THREADS={res["cpus"]} PYTHONDONTWRITEBYTECODE=1\n'
    text+='exec '+args+'\n'
    path=Path(p['output'])/'batch'/f'{kind}.sh';path.write_text(text);return path


def submit(p,dry=False,runner=subprocess.run):
    out=Path(p['output']);record=out/'batch/submissions.json'
    history=read_json(record) if record.exists() else []
    if (out/'complete.json').exists(): return dict(status='already_complete',submissions=history)
    states=query(history,runner) if history and not dry else {}
    owned_analysis=[x for x in history if x['kind']=='analyze']
    if any(states.get(str(x['job_id'])) in ['RUNNING','COMPLETING','CONFIGURING'] for x in owned_analysis):
        return dict(status='analysis_running',submissions=history)
    active_tasks={}
    for x in history:
        if x['kind']!='extract': continue
        for task in x['tasks']:
            key=f'{x["job_id"]}_{task}'
            if key in states: active_tasks[task]=key
    missing=[i for i,t in enumerate(p['tasks']) if not all(complete(p['shards'][j]) for j in t)]
    retry=[i for i in missing if i not in active_tasks];deps=[active_tasks[i] for i in missing if i in active_tasks]
    proposed=[]
    def invoke(kind,extra,tasks=None):
        res=p['site']['extraction' if kind=='extract' else 'analysis'];path=script(p,kind)
        cmd=['sbatch','--parsable',f'--job-name=pi0-{kind}',f'--cpus-per-task={res["cpus"]}',f'--mem={res["memory"]}',f'--time={res["time"]}',
             '--nodes=1','--ntasks=1',f'--output={out}/logs/{kind}-%A_%a.log']
        for key in ['account','partition','qos']:
            if p['site'][key]: cmd.append(f'--{key}={p["site"][key]}')
        cmd+=extra+[str(path)];proposed.append(cmd)
        if dry: return 'DRY_RUN_ARRAY'
        value=runner(cmd,capture_output=True,text=True,check=True).stdout.strip().split(';')[0]
        if not value.isdigit(): raise ValueError('Unexpected sbatch job ID: '+value)
        history.append(dict(kind=kind,job_id=value,tasks=tasks or [],command=cmd));write_json(record,history)
        return value
    # Replace only recorded workflow-owned pending analysis dependencies.
    for x in owned_analysis:
        jid=str(x['job_id'])
        if jid in states and not dry: runner(['scancel',jid],capture_output=True,text=True,check=True)
    if retry:
        array=','.join(map(str,retry))+f'%{p["site"]["concurrency"]}'
        deps.append(invoke('extract',['--array='+array],retry))
    extra=['--dependency=afterok:'+':'.join(deps)] if deps else []
    invoke('analyze',extra)
    write_json(out/'batch/proposed_submissions.json',proposed)
    return dict(status='dry_run' if dry else 'submitted',commands=proposed,submissions=history)


def worker(plan_path,kind,task=None):
    from .config import cache_signature
    from .candidates import CandidateSink
    from .extract import extract
    from .analysis import run
    from .selection import discover
    p=read_json(plan_path);cfg=p['config'];out=Path(p['output'])
    if kind=='extract':
        task=int(os.environ['SLURM_ARRAY_TASK_ID']) if task is None else int(task)
        for j in p['tasks'][task]:
            s=p['shards'][j]
            if complete(s): continue
            target=Path(s['output']);target.parent.mkdir(parents=True,exist_ok=True)
            # Only known incomplete products in this run-owned shard may be replaced.
            for file in [target,target.with_suffix('.json'),target.with_name(target.stem+'.partial.npz')]:
                if file.exists(): file.unlink()
            extract(s['inputs'],target,cfg,CandidateSink(cfg['input']['candidate_floor_GeV'],cache_signature(cfg)))
            write_json(target.parent/'complete.json',dict(signature=s['signature'],schema_version=1,complete=True))
        return 0
    if p['compact']: sources=discover(p['compact'])
    else:
        if not all(complete(s) for s in p['shards']): raise ValueError('Extraction incomplete; global analysis refused')
        cache=Path(p['cache']);sources=[Path(s['output']) for s in p['shards']]
        write_json(cache/'index.json',dict(schema_version=1,complete=True,shards=[str(x.relative_to(cache)) for x in sources]))
    result=run(sources,out,p['work'],cfg)
    if result==0: write_json(out/'complete.json',dict(schema_version=1,complete=True))
    return result
