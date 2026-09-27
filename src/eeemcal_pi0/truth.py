# SPDX-License-Identifier: LGPL-3.0-or-later
# Track checks adapt calibration/pi0/src/mc_diagnostics.py (fields/EventTruth); optional wiring is new.
"""Truth diagnostics optional; a missing track relation never disables the veto."""
import numpy as np
from .mc_diagnostics import fields


def truth_setup(tree,mapping,cfg):
    names={'EcalEndcapNClusters':cfg['input']['clusters'], 'MCParticles':cfg['input']['mc_particles'],
           'EcalEndcapNClusterAssociations':cfg['input']['associations'],
           'EcalEndcapNTrackClusterMatches':cfg['input']['track_matches'], 'CentralCKFTracks':cfg['input']['tracks']}
    ff=fields()
    for k,v in ff.items():
        for old,new in names.items(): v=v.replace(old,new)
        ff[k]=v
    tracks={k:v for k,v in ff.items() if k.startswith('track_')}
    track_ok=all(v in tree for v in tracks.values())
    if cfg['selection']['neutral'] and not track_ok:
        raise ValueError('Missing required track-match relations for neutral veto: '+str([v for v in tracks.values() if v not in tree]))
    truth={k:v for k,v in ff.items() if not k.startswith('track_')}
    truth_ok=cfg['input']['mode']!='data' and all(v in tree for v in truth.values())
    extra=dict(tracks if track_ok else {})
    if truth_ok: extra.update(truth)
    canonical={i:next((old for old,new in names.items() if n==new),n) for i,n in mapping.items()}
    return extra,truth_ok,track_ok,canonical


class RecoOnly:
    def __init__(self,a,mapping,energy):
        self.n_pi0=0; self.momentum=np.zeros((0,3)); self.track_matched=np.zeros(len(energy),dtype=bool)
        keys=['track_cluster_index','track_cluster_collectionID','track_track_index','track_track_collectionID']
        values=[a.get(k,[]) for k in keys]
        if len({len(x) for x in values})!=1: raise ValueError('Invalid track-match relation lengths')
        for i,c,t,tc in zip(*values):
            if not (0<=i<len(energy) and mapping.get(int(c))=='EcalEndcapNClusters' and t>=0 and mapping.get(int(tc))=='CentralCKFTracks'):
                raise ValueError('Invalid event-local track-match ObjectID')
            self.track_matched[i]=True
    def pair(self,i,j):
        return dict(truth_sim=np.array([-1,-1]),truth_pdg=np.zeros(2,dtype=int),truth_weight=np.zeros(2),
                    truth_weight_sum=np.zeros(2),truth_photon=np.array([-1,-1]),truth_parent=np.array([-1,-1]),
                    track_matched=self.track_matched[[i,j]])
