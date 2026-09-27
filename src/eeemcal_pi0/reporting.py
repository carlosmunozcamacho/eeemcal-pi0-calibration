# SPDX-License-Identifier: LGPL-3.0-or-later
"""Regenerate every figure from persisted numerical products, without refitting."""
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from .io import read_json
from .model import M0


def records(p):
    with open(p) as f: return list(csv.DictReader(f))


def render(out,cfg):
    out=Path(out);d=out/'diagnostics';s=read_json(out/'summary.json');fits=read_json(d/'mass_fits.json')
    rows=records(out/'gains.csv');m=records(d/'mass_spectra.csv')
    def column(rows,key): return np.array([float(r[key]) if r[key] else np.nan for r in rows])
    lo=column(m,'low_GeV');hi=column(m,'high_GeV');edges=np.r_[lo,hi[-1]];mid=(lo+hi)/2
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':240,'pdf.fonttype':42})
    def figure(): return plt.subplots(figsize=(7.2,4.8))
    def save(fig,name):
        fig.tight_layout()
        for fmt in cfg['qa']['formats']: fig.savefig(out/'plots'/f'{name}.{fmt}',bbox_inches='tight')
        plt.close(fig)
    fig,ax=figure();ax.stairs(column(m,'training'),edges,label='Selected training pairs',color='black')
    if 'training' in fits:
        p=fits['training'];ax.plot(mid,p['fit_expected'],label='Signal + background',color='tab:blue');ax.plot(mid,p['background_expected'],'--',label='Background',color='tab:orange')
        if 'window' in s:
            ax.axvspan(*s['window']['window_GeV'],color='tab:blue',alpha=.10,label='Frozen training window')
            ax.set_title(f"Peak {p['peak_GeV']*1000:.2f} MeV; σ {p['sigma_GeV']*1000:.2f} MeV\nSignal {s['window']['signal']:.0f}, purity {s['window']['purity']:.1%} in window")
    ax.set(xlabel=r'$m_{\gamma\gamma}$ [GeV]',ylabel='Training pairs / bin');ax.legend(fontsize=9);save(fig,'pi0_mass_training')
    fig,ax=figure();ax.stairs(column(m,'preselection_pairs'),edges,label='All cached neutral pairs/event' if cfg['selection']['neutral'] else 'All cached pairs/event')
    ax.stairs(column(m,'selected_pairs'),edges,label='Final selection: ≤1 pair/event')
    ax.set(xlabel=r'$m_{\gamma\gamma}$ [GeV]',ylabel='Pairs / bin',title='Same input events; preselection permits multiple pairs')
    ax.legend(fontsize=9);save(fig,'pi0_mass_before_after_selection')
    fig,ax=figure()
    for key,label in [('heldout_before','Before'),('heldout_after','After')]:
        p=fits.get(key);desc=label
        if p: desc+=f": μ={p['peak_GeV']*1000:.2f}, σ={p['sigma_GeV']*1000:.2f} MeV; σ/μ={p['sigma_GeV']/p['peak_GeV']:.3%}"
        ax.stairs(column(m,key),edges,label=desc)
    ax.axvline(M0,color='gray',linestyle=':',label='Nominal π⁰ mass')
    ax.set(xlabel=r'$m_{\gamma\gamma}$ [GeV]',ylabel='Held-out pairs / bin',title='Identical held-out pairs; no training-window cut')
    ax.legend(fontsize=8);save(fig,'pi0_mass_before_after_calibration')
    gain=column(rows,'relative_gain');use=column(rows,'meaningful_support');change=column(rows,'half_change')
    active=np.array([r['active']=='True' for r in rows]);flag=np.array([any(x in r['status'] for x in ['unsupported','unstable','rank_sensitive','failure','low_support','stability_warning','not_accepted','not_evaluated','global_rank_truncated']) for r in rows])
    fig,ax=figure();ax.hist(gain[active],bins=cfg['qa']['histogram_bins'],histtype='step',range=None)
    ax.set(xlabel='Candidate relative energy multiplier G',ylabel='Active crystals',title='Full coefficient range, including tails');save(fig,'gain_histogram')
    fig,ax=figure();ax.scatter(column(rows,'crystal_number'),gain,s=7,c='tab:blue');ax.scatter(column(rows,'crystal_number')[flag],gain[flag],s=17,c='tab:red',marker='x',label='QA flagged / unsupported')
    ax.set(xlabel='Crystal number (cellID ascending)',ylabel='Relative gain G');ax.legend();save(fig,'gain_vs_crystal_number')
    xy=np.column_stack([column(rows,'x_mm'),column(rows,'y_mm')]);known=np.isfinite(xy).all(axis=1)
    for values,label,name in [(gain,'Candidate relative gain G','gain_map'),(use,'Meaningful independent training events','support_map')]:
        fig,ax=figure();vmin=float(values.min(initial=1e99));vmax=float(values.max(initial=-1e99))
        if vmax<=vmin: vmin-=.01;vmax+=.01
        im=ax.scatter(xy[known,0],xy[known,1],c=values[known],s=13,marker='s',cmap='viridis',norm=Normalize(vmin,vmax))
        ax.scatter(xy[known&flag,0],xy[known&flag,1],s=12,c='crimson',marker='x',linewidths=.5,label='QA flagged / unsupported')
        fig.colorbar(im,ax=ax,label=label,ticks=np.linspace(vmin,vmax,6))
        ax.set(aspect='equal',xlabel='x [mm]',ylabel='y [mm]');ax.legend(fontsize=8);save(fig,name)
    fig,ax=figure();ax.hist(use,bins=cfg['qa']['histogram_bins'],histtype='step');ax.set(xlabel='Meaningful independent training uses',ylabel='Crystals');save(fig,'support_histogram')
    fig,ax=figure();ax.scatter(use[active],100*change[active],s=7);ax.axhline(cfg['qa']['stability_fail']*100,color='tab:red',linestyle='--')
    ax.set(xlabel='Meaningful training uses',ylabel='Worst half/full gain change [%]',title='Reproducibility, not a gain uncertainty' if cfg['qa']['halves'] else 'Training-half QA not evaluated');save(fig,'gain_stability')
    history=records(d/'convergence.csv');fig,ax=figure()
    if history:
        for key,label in [('max_gain_step','Maximum gain step'),('m2_bias_GeV2','Squared-mass bias [GeV²]'),('m2_residual_std_GeV2','Residual spread [GeV²]')]: ax.plot(column(history,'iteration'),column(history,key),marker='.',label=label)
        ax.set_yscale('symlog',linthresh=1e-7);ax.legend(fontsize=9)
    ax.set(xlabel='Iteration',title=s.get('termination','No numerical fit'));save(fig,'convergence')
    eigen=records(d/'eigenvalues.csv');fig,ax=figure()
    if eigen:
        ax.plot(column(eigen,'mode'),column(eigen,'eigenvalue'));ax.axhline(float(eigen[0]['cutoff']),color='tab:red',linestyle='--',label='Configured cutoff');ax.set_yscale('symlog',linthresh=1e-12);ax.legend()
    ax.set(xlabel='Mode (ascending eigenvalue)',ylabel='Calibration matrix eigenvalue');save(fig,'eigenvalues')
    n=s['normalization'];counts=s['counts'];final=s.get('final_matrix',{})
    lines=['# Pi0 candidate channel calibration','',f"Execution: **{s['execution_status']}**. Numerical: **{s['numerical_status']}**. Physics QA: **{s['physics_qa_status']}**. Normalization: **{n['status']}**.",'',
           f"{counts['files']:,} input files; {counts['processed_events']:,} processed reconstructed events; {counts['selected_pairs']:,} selected pairs. Training/validation event counts: {counts['training_events']:,}/{counts['validation_events']:,}. Training-window pairs: {s.get('training_selected',0):,}.",'',
           f"Broad preselection contains {counts['preselection_pairs']:,} pairs from {counts['preselection_events']:,} distinct events at the stored floor; final selection chooses at most one pair per event by largest energy sum and deterministic cluster-index tie. Selection never uses truth or mass. Nominal vertex fallback is explicitly configured; cache vertex/source flags are preserved.",'',
           f"Active crystals: {s['active_parameters']}/{s['expected_channels']} ({counts['coverage_denominator']}). Retained/discarded modes: {final.get('retained_modes',0)}/{final.get('discarded_modes',0)}. Condition: {final.get('condition')}. Gain min/Q25/median/Q75/max: {s['gain_quantiles']}.",'',
           f"Half checks: {s.get('stability',{'status':'NOT_EVALUATED'})}. Failed or disabled checks never establish accepted stability. Discarded modes are combinations of gains; retained projection is not a statistical uncertainty.",'']
    if before:=s.get('before'):
        after=s['after'];lines += [f"Held-out peak: {before['peak_GeV']*1000:.3f} → {after['peak_GeV']*1000:.3f} MeV; σ: {before['sigma_GeV']*1000:.3f} → {after['sigma_GeV']*1000:.3f} MeV; σ/peak: {before['sigma_GeV']/before['peak_GeV']:.4%} → {after['sigma_GeV']/after['peak_GeV']:.4%}. Descriptive relative-width change: {s['relative_width_change']:+.3%}. Absolute width scales with peak position. No paired significance test or standalone calorimeter energy-resolution claim is made.",'']
    if 'pb_inverse' in n: lines += [f"Equivalent luminosity: {n['nb_inverse']:.7g} nb⁻¹ = {n['pb_inverse']:.7g} pb⁻¹ = {n['fb_inverse']:.7g} fb⁻¹.",'']+[f"- At {t['instantaneous_cm2_s']:.3g} cm⁻² s⁻¹: {t['seconds']:.7g} s ({t['readable']})." for t in n['equivalent_times']]+['']
    else: lines += ['Luminosity unavailable: '+'; '.join(n.get('reasons',[]))+'. Calibration does not require normalization.','']
    lines+=['Normalization assumptions and component bookkeeping are in summary.json. Beam times are ideal luminosity equivalents; elapsed time, recorded/delivered luminosity and live time are distinct. No generated exposure is inferred from selected pairs.','',
            'Reasons: '+('; '.join(s['reasons']) or 'No configured numerical/QA warnings; constants remain candidates.')+'.','',
            'All values are relative to the declared baseline, with frozen shower membership and directions. Unity unsupported values mean unchanged, never calibrated. Candidate constants require detector review before application; this run neither approves nor deploys them.','',
            'Products: gains.csv, summary.json, config_resolved.yaml, provenance.json, diagnostics/ numerical plot inputs, plots/. Regenerate figures with `pi0-calibrate --stage plots --output THIS_DIRECTORY`.']
    (out/'report.md').write_text('\n'.join(lines)+'\n')
