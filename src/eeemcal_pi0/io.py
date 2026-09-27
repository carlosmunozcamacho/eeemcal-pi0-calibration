# SPDX-License-Identifier: LGPL-3.0-or-later
# clean/JSON helpers adapt calibration/pi0/src/dis_analysis.py; ownership/path logic is new.
"""Small atomic records and explicit ownership. Never remove input products."""
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import uuid
import numpy as np


def clean(x):
    if isinstance(x,dict): return {str(k):clean(v) for k,v in x.items()}
    if isinstance(x,(list,tuple,np.ndarray)): return [clean(v) for v in x]
    if isinstance(x,np.generic): x=x.item()
    if isinstance(x,float) and not np.isfinite(x): return None
    if isinstance(x,Path): return str(x)
    return x


def write_json(path,x):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+'.tmp'); tmp.write_text(json.dumps(clean(x),indent=2,allow_nan=False)+'\n');tmp.replace(path)


def read_json(path): return json.loads(Path(path).read_text())


def digest(x): return hashlib.sha256(json.dumps(clean(x),sort_keys=True).encode()).hexdigest()


def table(path,rows,fields=None):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    rows=list(rows)
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields or list(rows[0]) if rows else fields or [])
        w.writeheader();w.writerows(rows)


def manifest(path):
    path=Path(path).resolve(); result=[]
    for line in path.read_text().splitlines():
        s=line.strip()
        if not s or s.startswith('#'): continue
        if not s.startswith('root://'):
            p=Path(s).expanduser();s=os.path.abspath(path.parent/p if not p.is_absolute() else p)
        result.append(s)
    if not result: raise ValueError('Empty manifest')
    if len(result)!=len(set(result)): raise ValueError('Duplicate manifest entries; remove aliases/copies explicitly')
    return result


def overlaps(a,b):
    a,b=Path(a).resolve(),Path(b).resolve()
    return a==b or a in b.parents or b in a.parents


def prepare(output,identity,protected,resume=False,overwrite=False):
    output=Path(output).resolve()
    for p in protected:
        if p and not str(p).startswith('root://') and overlaps(output,p): raise ValueError(f'Input/output collision: {p}')
    marker=output/'.pi0-owner.json'; token=digest(identity)
    if output.exists():
        if not marker.exists(): raise FileExistsError('Output exists without producer ownership: '+str(output))
        old=read_json(marker)
        if old['identity']!=token: raise ValueError('Run identity/configuration differs; choose a new output directory')
        if not (resume or overwrite): raise FileExistsError('Output exists; use --resume or --overwrite')
        if overwrite:
            # The entire directory was created exclusively for this run. Reject symlink root.
            if Path(output).is_symlink(): raise ValueError('Refuse symlink output')
            shutil.rmtree(output)
    if not output.exists():
        output.mkdir(parents=True)
        write_json(marker,dict(schema_version=1,run_id=str(uuid.uuid4()),identity=token))
    for sub in ['logs','batch','diagnostics','plots']: (output/sub).mkdir(exist_ok=True)
    return read_json(marker)


def own_subdir(path,owner,protected):
    path=Path(path).resolve()
    for p in protected:
        if p and not str(p).startswith('root://') and overlaps(path,p): raise ValueError('Cache/work overlaps input: '+str(p))
    marker=path/'.pi0-owner.json'
    if path.exists() and not marker.exists(): raise FileExistsError('Cache/work directory is not owned by this run: '+str(path))
    if marker.exists() and read_json(marker)['run_id']!=owner['run_id']: raise ValueError('Cache/work owned by a different run')
    path.mkdir(parents=True,exist_ok=True);write_json(marker,owner)
    return path
