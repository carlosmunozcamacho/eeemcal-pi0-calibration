# SPDX-License-Identifier: LGPL-3.0-or-later
# Adapted from calibration/pi0/src/selection_B_calibration.py (main/final_matrix),
# dis_analysis.py (run_fit/matrix_report), and support_threshold_study.py (stability).
"""One global training decision and unweighted constrained fit."""
from pathlib import Path
import numpy as np
from .model import Sample,fit,evaluate,equations
from .diagnostics import meaningful_usage
from .peak_diagnostics import fit_peak,window_yields
from .selection import build_sample
from .normalization import resolve as normalization
from .io import write_json,table


def peak(mass,c):
    return fit_peak(mass,bins=c['bins'],background_form=c['background_form'],fit_range=c['fit_range'],criteria=c['criteria'])


def matrix(s,g,rows,active,c):
    if not len(active): return dict(retained_modes=0,discarded_modes=0,condition=None,eigenvalues=[],projection=[])
    C,*_=equations(s,g,rows,active,chunk=c['chunk']);eig,vec=np.linalg.eigh(C)
    keep=(eig>0)&(eig>c['rcond']*eig[-1])
    return dict(retained_modes=int(keep.sum()),discarded_modes=int((~keep).sum()),
                condition=float(eig[-1]/eig[keep][0]) if keep.any() else None,
                eigenvalues=eig.tolist(),projection=np.sum(vec[:,keep]**2,axis=1).tolist())


def run_fit(s,rows,active,c):
    settings={k:c[k] for k in ['rcond','iterations','max_step','step_tolerance','bias_tolerance','chunk']}
    g,_,_,history,reason=fit(s,rows,active_indices=active,**settings)
    return g,dict(termination=reason,history=history,retained_modes=history[-1].get('retained_modes',0))


def run(sources,out,work,cfg):
    out=Path(out);diag=out/'diagnostics';s,counts=build_sample(sources,work,cfg);a=s.arrays
    s.save(Path(work)/'sample.npz')
    (out/'inputs.txt').write_text('\n'.join(map(str,a['inputs']))+'\n')
    c=cfg['solver'];qa=cfg['qa'];train=a['train'].astype(bool);mass=np.sqrt(a['original_m2']);held=np.flatnonzero(~train)
    result=dict(schema_version=1,counts=counts,execution_status='PASS',numerical_status='FAIL',physics_qa_status='NOT_EVALUATED',
                reasons=[],normalization=normalization(cfg['normalization']),candidate_constants_only=True)
    if any(not cfg['input'][k] for k in ['geometry_id','baseline_id','cellid_encoding']):
        result['reasons'].append('Some compatibility identifiers are unspecified; equality of persisted positions alone does not prove baseline/encoding compatibility')
    components=result['normalization'].get('components',[])
    result['counts']['generated_events_represented']=sum(x['generated_events_represented'] for x in components) if components and 'pb_inverse' in result['normalization'] else None
    result['counts'].update(training_candidates=int(train.sum()),validation_candidates=len(held))
    gains=np.ones(s.k);use=np.zeros(s.k,dtype=int);active=np.array([],dtype=int);projection=np.full(s.k,np.nan);halfchange=np.full(s.k,np.nan)
    fits={};history=[];halves=[];before=after=None;half_valid=False
    try:
        p=peak(mass[train],cfg['peak']);fits['training']=p;result['training_peak']=p
        if not p['identified']: raise ValueError('No usable training pi0 peak under configured identification criteria')
        width=cfg['peak']['window_sigma'];window=window_yields(p,p['peak_GeV']-width*p['sigma_GeV'],p['peak_GeV']+width*p['sigma_GeV'])
        lo,hi=window['window_GeV'];rows=np.flatnonzero(train&(mass>=lo)&(mass<=hi))
        result.update(window=window,training_selected=len(rows))
        use=meaningful_usage(s,rows,c['meaningful_fraction']);active=np.flatnonzero(use>=c['min_usage'])
        result['initial_matrix']=matrix(s,gains,rows,active,c)
        if not len(active): raise ValueError('No support-qualified channels')
        gains,fitinfo=run_fit(s,rows,active,c);history=fitinfo['history'];result.update(fitinfo)
        final=matrix(s,gains,rows,active,c);result['final_matrix']=final;projection[active]=final['projection']
        result['numerical_status']=('PASS' if final['discarded_modes']==0 else 'WARNING') if fitinfo['termination']=='converged' else 'FAIL'
        if fitinfo['termination']!='converged': result['reasons'].append('Full fit did not converge: '+fitinfo['termination'])
        if final['discarded_modes']: result['reasons'].append('Active parameters include unresolved linear combinations')
        if qa['halves']:
            perm=np.random.default_rng(qa['half_seed']).permutation(rows);maps=[]
            for name,rr in [('A',perm[::2]),('B',perm[1::2])]:
                try:
                    g,h=run_fit(s,rr,active,c);maps.append(g);halves.append(dict(half=name,termination=h['termination'],retained_modes=h['retained_modes']))
                except (ValueError,np.linalg.LinAlgError) as e:
                    halves.append(dict(half=name,termination='failed',reason=str(e),retained_modes=0))
            if len(maps)==2:
                halfchange[active]=np.maximum(abs(maps[0][active]/gains[active]-1),abs(maps[1][active]/gains[active]-1))
            half_valid=len(maps)==2 and all(x['termination']=='converged' and x['retained_modes']==len(active) for x in halves)
            if not half_valid: result['reasons'].append('Half-fit convergence or retained rank failed; stability acceptance disabled')
        else: result['reasons'].append('Training-half QA not evaluated')
        result['half_states']=halves
        result['stability']=dict(status='PASS' if half_valid else 'NOT_EVALUATED' if not qa['halves'] else 'FAIL',
            median=float(np.nanmedian(halfchange[active])) if np.isfinite(halfchange).any() else None,
            p90=float(np.nanquantile(halfchange[active],.9)) if np.isfinite(halfchange).any() else None,
            below5=int(np.sum(halfchange[active]<.05)) if half_valid else None,
            below10=int(np.sum(halfchange[active]<.1)) if half_valid else None,
            fraction_below5=float(np.mean(halfchange[active]<.05)) if half_valid else None,
            fraction_below10=float(np.mean(halfchange[active]<.1)) if half_valid else None)
        before=peak(mass[held],cfg['peak']);after_mass=np.sqrt(evaluate(s,gains,held)[1]);after=peak(after_mass,cfg['peak'])
        fits.update(heldout_before=before,heldout_after=after);result.update(before=before,after=after)
        change=(after['sigma_GeV']/after['peak_GeV'])/(before['sigma_GeV']/before['peak_GeV'])-1
        result['relative_width_change']=change
        if not before['identified'] or not after['identified']: result['reasons'].append('Held-out peak QA inconclusive')
        if change>qa['width_warning']: result['reasons'].append('Descriptive relative mass width increase exceeds configured warning threshold; no significance test')
        if np.any(halfchange[active]>=qa['stability_fail']): result['reasons'].append('Some candidate channels exceed half-stability limit')
        result['physics_qa_status']='WARNING' if result['reasons'] else 'PASS'
        if qa['closure']:
            inject=np.ones(s.k);inject[active]=np.clip(1+np.random.default_rng(qa['closure_seed']).normal(0,.05,len(active)),.8,1.2)
            energy,m2=evaluate(s,inject);distorted=Sample(dict(a,energy=energy,original_m2=m2),tuple(x.multiply(inject).tocsr() for x in s.showers))
            recovered,ch=run_fit(distorted,rows,active,c);error=recovered[active]*inject[active]/gains[active]-1
            result['closure']=dict(termination=ch['termination'],rms=float(np.sqrt(np.mean(error**2))),max_abs=float(np.max(abs(error))),interpretation='Algebraic fixed-cluster closure against nominal map; not precision')
            table(diag/'closure.csv',[dict(cellID=int(a['cell_id'][k]),injected=inject[k],recovered=recovered[k],nominal=gains[k],relative_error=error[i]) for i,k in enumerate(active)])
        if qa['detailed_events']:
            table(diag/'selected_pairs.csv',[dict(input=a['inputs'][a['file_index'][i]],entry=int(a['entry'][i]),cluster1=int(a['cluster_index'][i,0]),cluster2=int(a['cluster_index'][i,1]),train=bool(train[i]),mass_GeV=float(mass[i])) for i in range(s.n)])
    except (ValueError,np.linalg.LinAlgError) as e:
        result['reasons'].append(str(e));result['numerical_status']='FAIL';result['physics_qa_status']='NOT_EVALUATED'
    result.update(active_parameters=len(active),expected_channels=s.k,unsupported=s.k-len(active),
                  gain_quantiles=np.quantile(gains[active],[0,.25,.5,.75,1]).tolist() if len(active) else [],
                  coverage={str(n):int(np.sum(use>=n)) for n in [1,5,10,20,50,100]})
    active_mask=np.zeros(s.k,dtype=bool);active_mask[active]=True
    channels=[]
    for k,cell in enumerate(a['cell_id']):
        flags=['candidate'] if active_mask[k] else ['unsupported','unchanged']
        if use[k]<qa['low_support']: flags.append('low_support')
        if active_mask[k]:
            if np.isfinite(halfchange[k]) and halfchange[k]>=qa['stability_fail']: flags.append('unstable')
            elif np.isfinite(halfchange[k]) and halfchange[k]>=qa['stability_warn']: flags.append('stability_warning')
            if not half_valid: flags.append('stability_not_accepted' if qa['halves'] else 'stability_not_evaluated')
            if projection[k]<qa['projection_min']: flags.append('rank_sensitive')
            if result.get('final_matrix',{}).get('discarded_modes',0): flags.append('global_rank_truncated')
        if result['numerical_status']=='FAIL' or result.get('termination','converged')!='converged': flags.append('numerical_failure')
        pos=a['cell_position'][k]
        channels.append(dict(cellID=int(cell),crystal_number=k+1,x_mm=pos[0],y_mm=pos[1],z_mm=pos[2],relative_gain=format(gains[k],'.17g'),
            meaningful_support=int(use[k]),active=bool(active_mask[k]),half_change=halfchange[k] if np.isfinite(halfchange[k]) else '',
            retained_projection=projection[k] if np.isfinite(projection[k]) else '',status=';'.join(flags)))
    table(out/'gains.csv',channels);table(diag/'channel_status.csv',channels)
    table(diag/'crystal_number_mapping.csv',[{k:r[k] for k in ['crystal_number','cellID','x_mm','y_mm','z_mm']} for r in channels])
    keys=list(dict.fromkeys(k for h in history for k in h));table(diag/'convergence.csv',history,keys)
    eig=result.get('final_matrix',{}).get('eigenvalues',[]);cut=c['rcond']*max(eig,default=0)
    table(diag/'eigenvalues.csv',[dict(mode=i+1,eigenvalue=v,cutoff=cut,retained=v>max(0,cut)) for i,v in enumerate(eig)],['mode','eigenvalue','cutoff','retained'])
    edges=np.array(counts['hist_edges']);selected=np.histogram(mass,edges)[0];hb=np.histogram(mass[held],edges)[0]
    ha=np.histogram(np.sqrt(evaluate(s,gains,held)[1]),edges)[0]
    table(diag/'mass_spectra.csv',[dict(low_GeV=l,high_GeV=h,preselection_pairs=pre,selected_pairs=sel,heldout_before=b,heldout_after=af,training=int(np.histogram(mass[train],edges)[0][i])) for i,(l,h,pre,sel,b,af) in enumerate(zip(edges[:-1],edges[1:],counts['preselection_hist'],selected,hb,ha))])
    table(diag/'cutflow.csv',[dict(stage=k,count=counts[k],unit='pairs' if k.endswith('pairs') else 'events') for k in ['processed_events','preselection_events','preselection_pairs','selected_pairs']])
    write_json(diag/'mass_fits.json',fits);write_json(out/'summary.json',result)
    from .reporting import render
    render(out,cfg)
    return 3 if result['numerical_status']=='FAIL' else 0
