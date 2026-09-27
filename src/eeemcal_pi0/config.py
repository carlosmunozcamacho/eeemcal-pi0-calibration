# SPDX-License-Identifier: LGPL-3.0-or-later
# Partition rule adapted from calibration/pi0/src/extract.py; other configuration code is new.
"""Validated physics settings. Site/scheduler settings are separate."""
import copy
import hashlib
from pathlib import Path
import yaml

DEFAULT = dict(
    input=dict(clusters='EcalEndcapNClusters', hits='EcalEndcapNRecHits',
               event_header='EventHeader', primary_vertices='PrimaryVertices',
               track_matches='EcalEndcapNTrackClusterMatches', tracks='CentralCKFTracks',
               mc_particles='MCParticles', associations='EcalEndcapNClusterAssociations',
               chunk=128, candidate_floor_GeV=.15, mode='auto',
               geometry_id=None, cellid_encoding=None, baseline_id=None, inventory=None, expected_manifest=None, identity_map=None,
               legacy_compatibility_acknowledged=False),
    vertex=dict(allow_nominal=True, nominal_mm=[0.,0.,0.], max_xy_mm=20., max_z_mm=300.),
    split=dict(training_fraction=.8, seed='pi0-v1', convention='blake2b-mod5'),
    selection=dict(energy_GeV=.2, neutral=True, asymmetry=.8, overlap=.2,
                   different_seeds=True, ranking='largest_energy_sum'),
    peak=dict(fit_range=[.04,.24], bins=50, background_form='exponential',
              criteria=dict(entries=60, signal=15, significance=4.,
                            mean_bounds=[.076,.194], sigma_bounds=[.0041,.0279]), window_sigma=2.),
    solver=dict(min_usage=10, meaningful_fraction=.05, rcond=.001, iterations=24,
                max_step=.1, step_tolerance=1e-4, bias_tolerance=1e-6, chunk=2048),
    qa=dict(halves=True, half_seed=20260924, stability_warn=.05, stability_fail=.1,
            low_support=20, projection_min=.99, width_warning=.05, closure=False,
            closure_seed=41256, detailed_events=False, formats=['png','pdf'], histogram_bins=60),
    normalization=dict(metadata=None, instantaneous_cm2_s=[1e33,1e34]))


def merge(base, override):
    for k,v in override.items():
        if k not in base: raise ValueError(f'Unknown configuration key: {k}')
        if isinstance(base[k], dict):
            if not isinstance(v,dict): raise ValueError(f'{k} must be a mapping')
            merge(base[k],v)
        else: base[k]=v
    return base


def resolve(path=None, overrides=()):
    cfg=copy.deepcopy(DEFAULT)
    if path:
        cfg=merge(cfg,yaml.safe_load(Path(path).read_text()) or {})
    for override in overrides:
        key, value=override.split('=',1); node=cfg; parts=key.split('.')
        for part in parts[:-1]: node=node[part]
        merge(node,{parts[-1]:yaml.safe_load(value)})
    for section,key in [('input','inventory'),('input','expected_manifest'),('input','identity_map'),('normalization','metadata')]:
        if cfg[section][key]:
            p=Path(cfg[section][key]).expanduser()
            if not p.is_absolute(): p=(Path(path).resolve().parent if path else Path.cwd())/p
            cfg[section][key]=str(p.resolve())
    validate(cfg)
    return cfg


def validate(c):
    if c['input']['mode'] not in ['auto','simulation','data']: raise ValueError('input.mode must be auto, simulation or data')
    if c['selection']['ranking']!='largest_energy_sum': raise ValueError('Only mass-independent largest_energy_sum ranking is supported')
    if c['peak']['background_form'] not in ['exponential','linear']: raise ValueError('Unsupported background form')
    if c['split']['convention'] not in ['blake2b-mod5','blake2b-fraction']: raise ValueError('Unknown split convention')
    if c['split']['convention']=='blake2b-mod5' and c['split']['training_fraction']!=.8: raise ValueError('mod5 preserves the legacy 80/20 partition; choose blake2b-fraction for other fractions')
    for section,keys in [('solver',['min_usage','iterations','chunk','max_step','step_tolerance','bias_tolerance']),('input',['chunk']),('peak',['bins','window_sigma']),('qa',['histogram_bins'])]:
        for k in keys:
            if not isinstance(c[section][k],(float,int)) or c[section][k]<=0: raise ValueError(f'{section}.{k} must be positive')
    for section,k in [('solver','min_usage'),('solver','iterations'),('solver','chunk'),('input','chunk'),('peak','bins'),('qa','histogram_bins')]:
        if not isinstance(c[section][k],int): raise ValueError(f'{section}.{k} must be an integer')
    for section,k in [('solver','rcond'),('solver','meaningful_fraction'),('split','training_fraction')]:
        if not 0<c[section][k]<1: raise ValueError(f'Invalid {section}.{k}')
    for k in ['asymmetry','overlap']:
        if c['selection'][k] is not None and not 0<c['selection'][k]<=1: raise ValueError(f'Invalid selection.{k}')
    if c['selection']['overlap'] is None: raise ValueError('overlap requires an explicit upper limit')
    if c['input']['candidate_floor_GeV']<0 or c['selection']['energy_GeV']<0: raise ValueError('Negative energy floor')
    if len(c['vertex']['nominal_mm'])!=3: raise ValueError('nominal_mm requires xyz')
    if c['peak']['fit_range'][0]>=c['peak']['fit_range'][1]: raise ValueError('Invalid mass range')
    if not set(c['qa']['formats'])<= {'png','pdf','svg'}: raise ValueError('Unsupported plot format')
    if not c['qa']['stability_warn']<=c['qa']['stability_fail']: raise ValueError('Stability thresholds reversed')
    if any(x<=0 for x in c['normalization']['instantaneous_cm2_s']): raise ValueError('Reference luminosities must be positive')


def partition(identity, entry, cfg):
    n=int.from_bytes(hashlib.blake2b(f"{identity}:{entry}:{cfg['seed']}".encode(),digest_size=8).digest(),'big')
    return n%5!=0 if cfg['convention']=='blake2b-mod5' else n/(2**64)<cfg['training_fraction']


def cache_signature(c):
    return dict(producer_version='0.1.0',input={k:v for k,v in c['input'].items() if k not in ['inventory','expected_manifest','identity_map','legacy_compatibility_acknowledged']},
                vertex=c['vertex'],split=c['split'])
