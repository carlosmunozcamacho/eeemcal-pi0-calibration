# SPDX-License-Identifier: LGPL-3.0-or-later
# New bounded reader informed by calibration/pi0/src/sample_luminosity.py metadata semantics.
# HepMC3 is a format reference, not copied implementation; see THIRD_PARTY_NOTICES.md.
"""Normalization only; never supplies weights to the unweighted calibration."""
from pathlib import Path
import yaml
import numpy as np

XS_TO_PB={'pb':1.,'nb':1e3,'ub':1e6,'mb':1e9,'b':1e12}
LUMI_TO_PB={'pb^-1':1.,'nb^-1':1e-3,'fb^-1':1e3}


def hepmc_final(path):
    """One final representative event, only attributes and weights, no particles."""
    import uproot
    import awkward as ak
    with uproot.open(path) as f:
        t=f['hepmc3_tree'];i=t.num_entries-1
        if i<0: raise ValueError('Empty EVGEN tree')
        def read(member):
            candidates=[k for k in t.keys() if (k.endswith('.'+member) or k.endswith('/'+member)) and k.startswith('hepmc3_event')]
            if len(candidates)!=1: raise ValueError('Ambiguous/missing HepMC '+member)
            return ak.to_list(t[candidates[0]].array(entry_start=i,entry_stop=i+1))[0]
        ids,names,values=read('attribute_id'),read('attribute_name'),read('attribute_string')
        records=[v for cid,n,v in zip(ids,names,values) if cid==0 and n=='GenCrossSection']
        if len(records)!=1: raise ValueError('Missing/ambiguous event GenCrossSection')
        fields=records[0].split();weights=read('weights')
        return dict(cross_section=dict(value=float(fields[0]),unit='pb'),error_pb=float(fields[1]),
                    serialized=records[0],representative_entry=i,representative_weights=weights,
                    weights='unit' if weights==[1.] else 'unknown',source=str(path))


def exposure(pb,instantaneous):
    return dict(pb_inverse=pb,nb_inverse=pb*1000,fb_inverse=pb/1000,
                equivalent_times=[dict(instantaneous_cm2_s=x,seconds=pb*1e36/x,
                    readable=f'{pb*1e36/x/3600:.3g} h' if pb*1e36/x>=3600 else f'{pb*1e36/x/60:.3g} min' if pb*1e36/x>=60 else f'{pb*1e36/x:.3g} s') for x in instantaneous])


def calculate(metadata,instantaneous):
    if not metadata: return dict(status='unavailable',reasons=['No normalization metadata supplied'])
    d=metadata;result=dict(source=d.get('source'),convention=d.get('convention'),assumptions=d.get('assumptions',[]),status='unavailable',reasons=[])
    if not d.get('source') or not d.get('convention'):
        result['reasons']=['Explicit source and convention are required'];return result
    if d.get('kind')=='data' or d.get('kind')=='reported_exposure':
        values={}
        for key in ['luminosity','delivered_luminosity','recorded_luminosity']:
            if key in d:
                v=d[key];values[key]=float(v['value'])*LUMI_TO_PB[v['unit']]
        result.update(metadata=d,luminosities_pb_inverse=values)
        key=d.get('use','recorded_luminosity' if d.get('kind')=='data' else 'luminosity')
        if key not in values: result['reasons']=['Requested luminosity convention missing'];return result
        result.update(status=d.get('precision','approximate'),**exposure(values[key],instantaneous));return result
    components=[]
    for c in d.get('components',[]):
        r=dict(c,status='unavailable')
        if c.get('evgen_file'):
            try:
                g=hepmc_final(c['evgen_file']);r['generator_metadata']=g
                c=dict(c,cross_section=g['cross_section'])
                if g['weights']!='unit': c['weights']='unknown'
            except Exception as e:
                r['reason']='Generator metadata unavailable: '+str(e);components.append(r);continue
        if c.get('weights')!='unit': r['reason']='Weighted/unknown bookkeeping: no supported exposure formula'
        elif not c.get('unit_weight_convention_confirmed'): r['reason']='Representative unit weight alone does not prove run convention'
        elif not c.get('generated_events_represented') or not c.get('count_convention'): r['reason']='Explicit represented generated count and loss/filter convention required; RECO counts are not substituted'
        elif not c.get('cross_section'): r['reason']='No cross section'
        else:
            sigma=float(c['cross_section']['value'])*XS_TO_PB[c['cross_section']['unit']]
            if not np.isfinite(sigma) or sigma<=0: raise ValueError('Invalid cross section')
            r.update(status=c.get('precision','approximate'),cross_section_pb=sigma,
                     pb_inverse=float(c['generated_events_represented'])/sigma)
        components.append(r)
    result['components']=components
    ids=[c.get('piece_id') for c in components]
    ensembles={c.get('ensemble') for c in components}
    if not components or any(c['status']=='unavailable' for c in components): result['reasons']=['Some component normalizations are unavailable']
    elif len(components)>1 and (not d.get('disjoint_generated_pieces') or None in ids or len(set(ids))!=len(ids)):
        result['reasons']=['Overlapping or undocumented generated-event pieces; no combined exposure']
    elif len(ensembles)!=1 or None in ensembles:
        result['reasons']=['Different/unspecified generation ensembles; explicit external stitching interpretation required']
    else:
        result.update(status='approximate' if any(c['status']!='exact' for c in components) else 'exact',
                      **exposure(sum(c['pb_inverse'] for c in components),instantaneous))
    return result


def resolve(cfg):
    try:
        path=cfg['metadata'];d=yaml.safe_load(Path(path).read_text()) if path else None
        if d:
            for c in d.get('components',[]):
                p=c.get('evgen_file')
                if p and not p.startswith('root://') and not Path(p).is_absolute(): c['evgen_file']=str((Path(path).parent/p).resolve())
        return calculate(d,cfg['instantaneous_cm2_s'])
    except Exception as e: return dict(status='unavailable',reasons=[str(e)])
