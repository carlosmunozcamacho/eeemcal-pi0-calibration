# SPDX-License-Identifier: LGPL-3.0-or-later
# Adapted from calibration/pi0/src/broad_candidates.py (select), sidis_compact.py
# (pair_shard/merge), and selection_B_calibration.py (build_sample).
"""Streaming mass-independent pair choice and bounded-shard global assembly."""
from itertools import combinations
from pathlib import Path
import csv
import json
import numpy as np
from scipy import sparse
from .candidates import load,hit,truth_proxy
from .diagnostics import pair_features
from .model import Sample,angular_factor,evaluate
from .io import read_json,write_json,manifest


def discover(path):
    p=Path(path)
    if p.is_file(): return [p.resolve()]
    index=p/'index.json'
    if index.exists():
        d=read_json(index)
        files=[(p/x).resolve() for x in d['shards']]
        if not d.get('complete') or not all(x.exists() and x.with_suffix('.json').exists() for x in files):
            raise ValueError('Incomplete compact production; analysis refused')
        return files
    root=p/'shards' if (p/'shards').is_dir() else p
    files=sorted(root.glob('shard_*/candidates.npz'))
    if not files and (p/'candidates.npz').exists(): files=[p/'candidates.npz']
    if not files: raise ValueError('No broad compact shards found')
    return [x.resolve() for x in files]


def compatible(a,cfg):
    if int(a.get('candidate_schema',-1))!=1: raise ValueError('Unsupported broad compact schema')
    if float(a['candidate_floor_GeV'])>cfg['selection']['energy_GeV']: raise ValueError('Requested photon energy is below stored candidate floor; cache cannot support relaxation')
    if cfg['selection']['neutral'] and ('track_matched' not in a or not bool(a.get('track_matching_available',True))):
        raise ValueError('Compact sample lacks required track matching')
    if a['cell_id'].dtype!=np.uint64 or a['hit_cell_id'].dtype!=np.uint64: raise ValueError('cellIDs must be stored as uint64')
    metadata=json.loads(str(a['cache_metadata'])) if 'cache_metadata' in a else None
    if metadata:
        for k in ['geometry_id','cellid_encoding','baseline_id','clusters','hits']:
            if metadata['input'][k]!=cfg['input'][k]: raise ValueError('Incompatible cache '+k)
        if metadata['split']!=cfg['split'] or metadata['vertex']!=cfg['vertex']: raise ValueError('Cache split/vertex differs; existing assignments cannot be changed')
    elif not (cfg['input']['legacy_compatibility_acknowledged'] and cfg['input']['geometry_id'] and cfg['input']['baseline_id'] and cfg['input']['cellid_encoding']):
        raise ValueError('Legacy cache lacks baseline/geometry declaration; explicit legacy_compatibility_acknowledged required')
    elif not cfg['input']['expected_manifest']:
        raise ValueError('Legacy compact reuse requires input.expected_manifest to reject missing production shards')
    elif cfg['split']!=dict(training_fraction=.8,seed='pi0-v1',convention='blake2b-mod5'):
        raise ValueError('Legacy cache preserves its original pi0-v1 mod5 split')


def choose(a,cfg,edges):
    hist=np.zeros(len(edges)-1,dtype=np.int64); pair_count=0; event_count=0
    selected=[]; c=cfg['selection']
    for ev,(lo,hi) in enumerate(zip(a['event_offsets'][:-1],a['event_offsets'][1:])):
        if hi-lo<2: continue
        # Broad preselection: all cached neutral (when enabled) pairs, all events.
        broad=[k for k in range(lo,hi) if not c['neutral'] or not a['track_matched'][k]]
        best=None; bestkey=None
        if len(broad)>=2: event_count+=1
        for i,j in combinations(broad,2):
            e=a['energy'][[i,j]];pos=a['position'][[i,j]];v=a['vertex'][i]
            mass=np.sqrt(float(np.prod(e)*angular_factor(pos[None],v[None])[0]))
            pair_count+=1
            b=np.searchsorted(edges,mass,side='right')-1
            if mass==edges[-1]: b=len(hist)-1
            if 0<=b<len(hist): hist[b]+=1
            if min(e)<c['energy_GeV']: continue
            feature=pair_features(hit(a,i),hit(a,j),e,pos)
            if c['asymmetry'] is not None and feature['asymmetry']>=c['asymmetry']: continue
            if feature['overlap']>=c['overlap'] or (c['different_seeds'] and feature['same_seed']): continue
            key=(-float(e.sum()),int(a['cluster_index'][i]),int(a['cluster_index'][j]))
            if bestkey is None or key<bestkey: bestkey=key;best=(i,j)
        if best is not None: selected.append((ev,*best))
    return selected,dict(preselection_pairs=pair_count,preselection_events=event_count,preselection_hist=hist.tolist())


def pair_shard(source,out,cfg,edges):
    a=load(source);compatible(a,cfg)
    info=read_json(source.with_suffix('.json'))
    if info.get('development_limit_per_file') is not None or info.get('relation_errors',0) or info.get('nonconserving_hit_assignments',0):
        raise ValueError('Incomplete/invalid extraction shard: '+str(source))
    if 'files' in info and isinstance(info['files'],list) and info['files'] and isinstance(info['files'][0],dict):
        if any(x['events_read']!=x['events_in_file'] for x in info['files']): raise ValueError('Partial RECO extraction')
    for fi in range(len(a['inputs'])):
        entries=a['event_entry'][a['event_file_index']==fi]
        if len(entries)!=len(np.unique(entries)): raise ValueError('Duplicate compact event identity')
    if info.get('files') and isinstance(info['files'][0],dict):
        if [x['path'] for x in info['files']]!=a['inputs'].tolist(): raise ValueError('Compact input bookkeeping mismatch')
        if sum(x['events_read'] for x in info['files'])!=len(a['event_entry']): raise ValueError('Compact event bookkeeping mismatch')
    chosen,counts=choose(a,cfg,edges); records=[];parts=[[[],[],[0]],[[],[],[0]]]
    for ev,i,j in chosen:
        r=dict(file_index=int(a['event_file_index'][ev]),entry=int(a['event_entry'][ev]),
               cluster_index=a['cluster_index'][[i,j]],train=bool(a['event_train'][ev]),
               energy=a['energy'][[i,j]],position=a['position'][[i,j]],vertex=a['vertex'][i],vertex_source=a['vertex_source'][i])
        lo,hi=a['event_offsets'][ev:ev+2]
        r['truth_label']=truth_proxy(a,lo,hi).pair(i-lo,j-lo)['truth_label'] if 'truth_sim' in a else 'not_available'
        records.append(r)
        for q,k in enumerate([i,j]):
            ids,c=hit(a,k); col=np.searchsorted(a['cell_id'],ids)
            if np.any(col>=len(a['cell_id'])) or not np.array_equal(a['cell_id'][col],ids): raise ValueError('Unknown contribution cellID')
            d,ix,ptr=parts[q];d.extend(c);ix.extend(col);ptr.append(len(d))
    if records: b={k:np.array([r[k] for r in records]) for k in records[0]}
    else:
        b=dict(file_index=np.empty(0,dtype=int),entry=np.empty(0,dtype=int),cluster_index=np.empty((0,2),dtype=int),
               train=np.empty(0,dtype=bool),energy=np.empty((0,2)),position=np.empty((0,2,3)),vertex=np.empty((0,3)),
               vertex_source=np.empty(0,dtype='U64'),truth_label=np.empty(0,dtype='U64'))
    b.update(cell_id=a['cell_id'],cell_position=a['cell_position'],inputs=a['inputs'])
    b['angle_factor']=angular_factor(b['position'],b['vertex']);b['original_m2']=np.prod(b['energy'],axis=1)*b['angle_factor']
    matrices=[]
    for d,ix,ptr in parts:
        m=sparse.csr_matrix((d,np.array(ix,dtype=int),ptr),shape=(len(records),len(a['cell_id'])))
        m.sum_duplicates();m.eliminate_zeros();matrices.append(m)
    Sample(b,tuple(matrices)).save(out)
    return dict(events=len(a['event_entry']),training_events=int(a['event_train'].sum()),pairs=len(records),files=a['inputs'].tolist(),
                shapes={k:list(v.shape[1:]) for k,v in b.items() if k not in ['cell_id','cell_position','inputs']},
                dtypes={k:v.dtype.str for k,v in b.items() if k not in ['cell_id','cell_position','inputs']},
                nnz=[m.nnz for m in matrices],**counts)


def build_sample(sources,work,cfg):
    work=Path(work);pairs=work/'pairs';pairs.mkdir(exist_ok=True)
    edges=np.linspace(*cfg['peak']['fit_range'],cfg['peak']['bins']+1)
    metadata=[];geometry={};inputs=[];dtypes={}
    for i,source in enumerate(sources):
        out=pairs/f'{i:06d}.npz';d=pair_shard(source,out,cfg,edges);metadata.append(d);inputs.extend(d['files'])
        for k,v in d['dtypes'].items(): dtypes[k]=np.promote_types(dtypes.get(k,np.dtype(v)),np.dtype(v))
        with np.load(out) as a:
            for cell,pos in zip(a['cell_id'],a['cell_position']):
                if int(cell) in geometry: np.testing.assert_allclose(geometry[int(cell)],pos,atol=.01,rtol=0)
                geometry[int(cell)]=pos
        print(f'SELECT shard={i+1}/{len(sources)} pairs={d["pairs"]}',flush=True)
    if len(inputs)!=len(set(inputs)): raise ValueError('Duplicate canonical input identities across shards')
    if cfg['input']['expected_manifest'] and inputs!=manifest(cfg['input']['expected_manifest']):
        raise ValueError('Compact inputs differ from expected complete manifest')
    inventory=cfg['input']['inventory']
    if inventory:
        with open(inventory) as f: inv=list(csv.DictReader(f))
        if len({int(r['cellID']) for r in inv})!=len(inv): raise ValueError('Duplicate inventory cellID')
        expected={int(r['cellID']):np.array([float(r.get(k) or 'nan') for k in ['x_mm','y_mm','z_mm']]) for r in inv}
        if not set(geometry)<=set(expected): raise ValueError('Observed cell outside configured detector inventory')
        for cell,pos in geometry.items():
            good=np.isfinite(expected[cell]);np.testing.assert_allclose(expected[cell][good],pos[good],atol=.01,rtol=0)
            expected[cell]=pos
        geometry=expected
    cells=np.array(sorted(geometry),dtype=np.uint64);n=sum(d['pairs'] for d in metadata)
    buffers=work/'buffers';buffers.mkdir(exist_ok=True)
    def alloc(name,dtype,shape):
        if np.prod(shape)==0: return np.empty(shape,dtype=dtype)
        return np.lib.format.open_memmap(buffers/(name+'.npy'),mode='w+',dtype=dtype,shape=shape)
    shapes=metadata[0]['shapes'];a={k:alloc(k,dtypes[k],(n,*v)) for k,v in shapes.items()}
    nnz=np.sum([d['nnz'] for d in metadata],axis=0)
    data=[alloc(f'c{q}_data',float,(int(nnz[q]),)) for q in range(2)]
    indices=[alloc(f'c{q}_indices',np.int32,(int(nnz[q]),)) for q in range(2)]
    ptr=[alloc(f'c{q}_indptr',np.int64,(n+1,)) for q in range(2)]
    row=0;offset=0;used=[0,0]
    for i,d in enumerate(metadata):
        s=Sample.load(pairs/f'{i:06d}.npz')
        for k in shapes: a[k][row:row+s.n]=s.arrays[k]+offset if k=='file_index' else s.arrays[k]
        remap=np.searchsorted(cells,s.arrays['cell_id'])
        for q,c in enumerate(s.showers):
            end=used[q]+c.nnz;data[q][used[q]:end]=c.data;indices[q][used[q]:end]=remap[c.indices]
            ptr[q][row:row+s.n]=c.indptr[:-1]+used[q];used[q]=end
        row+=s.n;offset+=len(d['files'])
    for q in range(2): ptr[q][-1]=used[q]
    a.update(cell_id=cells,cell_position=np.array([geometry[int(c)] for c in cells]),inputs=np.array(inputs))
    sample=Sample(a,tuple(sparse.csr_matrix((data[q],indices[q],ptr[q]),shape=(n,len(cells)),copy=False) for q in range(2)))
    energy,m2=evaluate(sample,np.ones(sample.k))
    np.testing.assert_allclose(energy,a['energy'],rtol=2e-15);np.testing.assert_allclose(m2,a['original_m2'],rtol=3e-15)
    counts=dict(files=len(inputs),processed_events=sum(d['events'] for d in metadata),available_reconstructed_events=sum(d['events'] for d in metadata),
                training_events=sum(d['training_events'] for d in metadata),selected_pairs=n,
                preselection_pairs=sum(d['preselection_pairs'] for d in metadata),preselection_events=sum(d['preselection_events'] for d in metadata),
                preselection_hist=np.sum([d['preselection_hist'] for d in metadata],axis=0).tolist(),hist_edges=edges.tolist(),
                coverage_denominator='configured_inventory' if inventory else 'observed_channels',channels=len(cells))
    counts['validation_events']=counts['processed_events']-counts['training_events']
    return sample,counts
