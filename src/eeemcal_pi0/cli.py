# SPDX-License-Identifier: LGPL-3.0-or-later
"""Portable producer CLI; submission returns immediately."""
import argparse
import fcntl
import importlib.metadata
import json
from pathlib import Path
import platform
import sys
import traceback
import yaml
from . import __version__
from .config import resolve
from .io import prepare,own_subdir,manifest,write_json,read_json,overlaps
from . import batch


def parser():
    p=argparse.ArgumentParser(description='Candidate relative pi0 channel calibration; no automatic deployment')
    source=p.add_mutually_exclusive_group();source.add_argument('--input',type=Path);source.add_argument('--compact',type=Path)
    p.add_argument('--output',type=Path);p.add_argument('--config',type=Path);p.add_argument('--set',action='append',default=[],metavar='SECTION.KEY=YAML')
    p.add_argument('--cache-dir',type=Path);p.add_argument('--work-dir',type=Path)
    p.add_argument('--executor',choices=['local','slurm'],default='local');p.add_argument('--batch-config',type=Path)
    p.add_argument('--stage',choices=['all','extract','analyze','plots','status','worker'],default='all')
    policy=p.add_mutually_exclusive_group();policy.add_argument('--resume',action='store_true');policy.add_argument('--overwrite',action='store_true')
    p.add_argument('--dry-run',action='store_true');p.add_argument('--plan',type=Path);p.add_argument('--worker-kind',choices=['extract','analyze']);p.add_argument('--task',type=int)
    return p


def main(argv=None):
    p=parser();args=p.parse_args(argv);out=None;owned=False
    try:
        if args.stage=='worker':
            out=Path(read_json(args.plan)['output']);owned=True
            return batch.worker(args.plan,args.worker_kind,args.task)
        if args.output is None: p.error('--output is required')
        out=args.output.resolve()
        if args.stage=='status':
            submissions=out/'batch/submissions.json'
            print(json.dumps(batch.query(read_json(submissions) if submissions.exists() else []),indent=2));return 0
        if args.stage=='plots':
            from .reporting import render
            render(out,yaml.safe_load((out/'config_resolved.yaml').read_text()));return 0
        if not (args.input or args.compact): p.error('Provide --input or --compact')
        if args.stage=='analyze' and not args.compact: p.error('--stage analyze requires --compact')
        if args.stage=='extract' and not args.input: p.error('--stage extract requires --input')
        cfg=resolve(args.config,args.set);site=batch.site_config(args.batch_config)
        inputs=manifest(args.input) if args.input else [];compact=str(args.compact.resolve()) if args.compact else None
        protected=inputs+[x for x in [compact,args.input,args.config,args.batch_config,cfg['input']['inventory'],cfg['input']['expected_manifest'],cfg['input']['identity_map'],cfg['normalization']['metadata']] if x]
        cache=args.cache_dir.resolve() if args.cache_dir else out/'cache';work=args.work_dir.resolve() if args.work_dir else out/'work'
        if overlaps(cache,work): raise ValueError('Cache and work directories must not overlap')
        for path in [cache,work]:
            if path==out or path in out.parents: raise ValueError('Cache/work may be children of output or separate, never its ancestors')
        identity=dict(inputs=inputs,compact=compact,config=cfg,cache=str(cache),work=str(work),site=site)
        owner=prepare(out,identity,protected,args.resume or args.overwrite,False)
        owned=True
        own_subdir(work,owner,protected)
        if not compact: own_subdir(cache,owner,protected)
        if args.resume and (out/'complete.json').exists():
            print('Already complete:',out);return 0
        (out/'config_resolved.yaml').write_text(yaml.safe_dump(cfg,sort_keys=False))
        if not (out/'inputs.txt').exists(): (out/'inputs.txt').write_text('\n'.join(inputs or [compact])+'\n')
        versions={x:importlib.metadata.version(x) for x in ['numpy','scipy','matplotlib','uproot','awkward','PyYAML']}
        write_json(out/'provenance.json',dict(schema_version=1,producer_version=__version__,python=platform.python_version(),dependencies=versions,
            run_id=owner['run_id'],input_manifest=str(args.input.resolve()) if args.input else None,compact=compact,cache=str(cache),work=str(work),
            baseline_id=cfg['input']['baseline_id'],geometry_id=cfg['input']['geometry_id'],cellid_encoding=cfg['input']['cellid_encoding']))
        lock=(out/'batch/submit.lock').open('w')
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        if args.overwrite:
            previous=out/'batch/submissions.json'
            if previous.exists() and batch.query(read_json(previous)):
                raise ValueError('Cannot overwrite while workflow jobs are queued/running')
            # Replace only standard result products; retain broad caches and all unknown files.
            for name in ['gains.csv','summary.json','report.md','complete.json']:
                file=out/name
                if file.exists(): file.unlink()
            diagnostic_names=['channel_status.csv','crystal_number_mapping.csv','convergence.csv','eigenvalues.csv','mass_spectra.csv','cutflow.csv','mass_fits.json','closure.csv','selected_pairs.csv']
            plot_names=['pi0_mass_training','pi0_mass_before_after_selection','pi0_mass_before_after_calibration','gain_histogram','gain_map','gain_vs_crystal_number','support_map','support_histogram','gain_stability','convergence','eigenvalues']
            standard=[out/'diagnostics'/name for name in diagnostic_names]+[out/'plots'/(name+'.'+fmt) for name in plot_names for fmt in ['png','pdf','svg']]
            for file in standard:
                if file.is_file(): file.unlink()
        planfile=out/'batch/work_plan.json'
        plan=read_json(planfile) if (args.resume or args.overwrite) and planfile.exists() else batch.plan(out,inputs,compact,cache,work,cfg,site)
        if args.executor=='slurm' or args.dry_run:
            if args.stage!='all': raise ValueError('Use local stages inside external jobs, or all for Slurm orchestration')
            print(json.dumps(batch.submit(plan,args.dry_run),indent=2));return 0
        if args.resume and not args.overwrite and (out/'complete.json').exists(): print('Already complete:',out);return 0
        if args.stage in ['all','extract'] and inputs:
            for task in range(len(plan['tasks'])): batch.worker(planfile,'extract',task)
            if args.stage=='extract':
                write_json(cache/'index.json',dict(schema_version=1,complete=True,shards=[str(Path(s['output']).relative_to(cache)) for s in plan['shards']]))
                print('Extraction complete:',cache);return 0
        return batch.worker(planfile,'analyze')
    except Exception as e:
        print(f'pi0-calibrate: {e}',file=sys.stderr)
        if owned and out and (out/'.pi0-owner.json').exists() and not isinstance(e,(FileExistsError,BlockingIOError)):
            # Do not replace a completed result with a failed invocation.
            if not (out/'summary.json').exists():
                write_json(out/'summary.json',dict(schema_version=1,execution_status='FAIL',numerical_status='NOT_EVALUATED',physics_qa_status='NOT_EVALUATED',normalization=dict(status='unavailable'),reasons=[str(e)]))
                (out/'report.md').write_text('# Producer did not complete\n\n'+str(e)+'\nNo constants are approved.\n')
            (out/'logs/error.txt').write_text(traceback.format_exc())
        return 2
