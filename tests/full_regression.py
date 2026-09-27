# SPDX-License-Identifier: LGPL-3.0-or-later
"""One compact-only regression; all production paths supplied explicitly.

Called within one allocated job. On success, archive includes this validation.
No central ROOT reading, no scheduler submission, no result retuning.
"""
import argparse
import csv
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import traceback
import numpy as np
from eeemcal_pi0.cli import main as producer
from eeemcal_pi0.model import Sample
from eeemcal_pi0.io import read_json,write_json


def compare(output,work,reference):
    s=read_json(output/'summary.json');ref=read_json(reference/'summary.json');r=ref['fits']['10']
    old=Sample.load(reference/'sample.npz');new=Sample.load(work/'sample.npz')
    checks={}
    def exact(name,a,b):
        np.testing.assert_array_equal(a,b,err_msg=name);checks[name]='exact'
    def close(name,a,b,rtol,atol):
        np.testing.assert_allclose(a,b,rtol=rtol,atol=atol,err_msg=name)
        checks[name]=dict(rtol=rtol,atol=atol,max_abs=float(np.max(abs(np.asarray(a)-np.asarray(b)))))
    for key in ['inputs','file_index','entry','cluster_index','train','cell_id']:
        exact('sample_'+key,new.arrays[key],old.arrays[key])
    for key in ['energy','position','vertex','angle_factor','original_m2']:
        close('sample_'+key,new.arrays[key],old.arrays[key],3e-14,1e-18)
    for q in range(2):
        exact(f'c{q}_indices',new.showers[q].indices,old.showers[q].indices)
        exact(f'c{q}_indptr',new.showers[q].indptr,old.showers[q].indptr)
        close(f'c{q}_contributions',new.showers[q].data,old.showers[q].data,3e-14,1e-18)
    with (output/'gains.csv').open() as f:table=list(csv.DictReader(f))
    with (reference/'gains_min10.csv').open() as f:oldtable=list(csv.DictReader(f))
    exact('gain_cellIDs',[int(x['cellID']) for x in table],[int(x['cellID']) for x in oldtable])
    exact('support_counts',[int(x['meaningful_support']) for x in table],[int(x['usage']) for x in oldtable])
    close('individual_gains',[float(x['relative_gain']) for x in table],[float(x['cumulative_gain']) for x in oldtable],2e-5,2e-7)
    exact('training_window_pairs',s['training_selected'],ref['training_selected'])
    exact('files',s['counts']['files'],ref['files']);exact('events',s['counts']['processed_events'],ref['events'])
    exact('active',s['active_parameters'],r['active']);exact('retained',s['final_matrix']['retained_modes'],r['retained_modes'])
    exact('half_stable_below10',s['stability']['below10'],r['stability']['worst_of_halves']['below10'])
    close('window',s['window']['window_GeV'],ref['window_GeV'],2e-5,2e-7)
    for label in ['before','after']:
        for key in ['peak_GeV','sigma_GeV']:
            close(label+'_'+key,s[label][key],r[label][key],2e-5,2e-7)
    with (reference/'stability_min10.csv').open() as f:st=list(csv.DictReader(f))
    close('individual_half_changes',[float(x['half_change']) for x in table],[float(x['worst_half_relative_change']) for x in st],2e-4,2e-5)
    if s['termination']!='converged': raise AssertionError('New fit did not converge')
    return dict(status='PASS',checks=checks,counts=s['counts'],training_selected=s['training_selected'],
                active=s['active_parameters'],retained=s['final_matrix']['retained_modes'],gain_quantiles=s['gain_quantiles'],
                half_stable10=s['stability']['below10'],before=s['before'],after=s['after'],normalization=s['normalization'])


def archive(source,path):
    path.parent.mkdir(parents=True,exist_ok=True)
    def filtered(info):
        parts=Path(info.name).parts
        roots={'CITATION.cff','CHANGELOG.md','RELEASE_NOTES_v0.1.0.md','.github','README.md','CONTRIBUTING.md','.gitignore','LICENSE','COPYING','NOTICE.md',
               'THIRD_PARTY_NOTICES.md','pyproject.toml','MANIFEST.in','docs','examples','schemas','src','tests'}
        if len(parts)>1 and parts[1] not in roots: return None
        if any(x in ['__pycache__','build','dist','.git','.pytest_cache','.venv','venv'] or x.endswith('.egg-info') for x in parts): return None
        if info.name.endswith(('.pyc','.root','.npz','.log','.out')): return None
        info.uid=info.gid=0;info.uname=info.gname=''
        return info
    with tarfile.open(path,'w:gz') as tar:tar.add(source,arcname='eeemcal-pi0-0.1.0',filter=filtered)


def main():
    p=argparse.ArgumentParser();p.add_argument('--compact',type=Path,required=True);p.add_argument('--reference',type=Path,required=True)
    p.add_argument('--config',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True)
    p.add_argument('--source',type=Path,required=True);p.add_argument('--archive',type=Path,required=True);args=p.parse_args()
    report=args.source/'docs/full_regression.md';comparison=args.source/'docs/full_regression.json'
    try:
        code=producer(['--compact',str(args.compact),'--config',str(args.config),'--output',str(args.output),'--work-dir',str(args.work),'--executor','local'])
        if code:raise RuntimeError(f'Producer exit {code}; see output/report.md')
        result=compare(args.output,args.work,args.reference);write_json(comparison,result)
        b,a=result['before'],result['after']
        text=f'''# Full compact regression: PASS

Compared against the saved detailed relaxed-B/min10 result. All event/pair IDs,
train flags, cellIDs, sparse structure, per-cell support and training-window
pair counts match exactly. Numerical tolerances and differences are recorded
in full_regression.json; no selection or solver retuning was performed.

Files: {result['counts']['files']}; reconstructed events: {result['counts']['processed_events']}.
Training-window pairs: {result['training_selected']}; active/retained: {result['active']}/{result['retained']}.
Median gain: {result['gain_quantiles'][2]:.12g}; half-stable below10%: {result['half_stable10']}.
Held-out peak: {b['peak_GeV']*1000:.6f} → {a['peak_GeV']*1000:.6f} MeV;
relative width: {b['sigma_GeV']/b['peak_GeV']:.8%} → {a['sigma_GeV']/a['peak_GeV']:.8%}.

Only existing broad compact inputs and saved reference products were read.
This reproduces the fixed-shower calculation; it does not validate real data,
a full reconstruction feedback loop or detector systematic errors.
Normalization was supplied as existing approximate
regression metadata, independently of the gain calculation.
'''
        report.write_text(text)
        (args.output/'regression_validation.md').write_text(text);write_json(args.output/'regression_comparison.json',result)
        archive(args.source,args.archive)
        print(json.dumps(dict(status='PASS',report=str(report),archive=str(args.archive))),flush=True)
        return 0
    except Exception as e:
        result=dict(status='FAIL',reason=str(e),traceback=traceback.format_exc());write_json(comparison,result)
        report.write_text('# Full compact regression: FAIL\n\n'+str(e)+'\n\nNo validated archive was produced. See full_regression.json and analysis output.\n')
        print(json.dumps(result),flush=True);return 1

if __name__=='__main__':raise SystemExit(main())
