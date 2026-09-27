# SPDX-License-Identifier: LGPL-3.0-or-later
# Adapted from calibration/pi0/src/broad_candidates.py (CandidateSink/load/hit/truth_proxy).
"""Candidate-level DIS compact storage and mass-independent offline pair choice.

Uses the existing extractor, pair_features and EventTruth.pair without changing
relation, overlap, seed, vertex or truth-label semantics. No calibration fits.
"""
import argparse
import json
from itertools import combinations
from pathlib import Path
import numpy as np
from .diagnostics import pair_features
from .mc_diagnostics import EventTruth
from .model import angular_factor


class CandidateSink:
    def __init__(self, floor=.15, metadata=None):
        self.metadata = metadata or {}
        self.floor = floor
        self.rows, self.ids, self.contributions = [], [], []
        self.hit_offsets, self.event_offsets = [0], [0]

    def add(self, event_index, energy, position, vertex, source, valid, hits, truth):
        for j in np.flatnonzero(valid & (energy >= self.floor)):
            # Use the established truth diagnostic to retain all pair-independent fields.
            t = truth.pair(j, j)
            row = dict(candidate_event=event_index, cluster_index=int(j), energy=energy[j],
                       position=position[j], vertex=vertex, vertex_source=source)
            for key in ['truth_sim', 'truth_pdg', 'truth_weight', 'truth_weight_sum',
                        'truth_photon', 'truth_parent', 'track_matched']:
                row[key] = t[key][0]
            origin = int(t['truth_photon'][0])
            row['truth_photon_momentum'] = truth.momentum[origin] if origin >= 0 else np.zeros(3)
            ids, c = hits[j]
            # Preserve original relation order, including duplicate IDs, for exact seed ties.
            self.ids.extend(ids); self.contributions.extend(c)
            self.hit_offsets.append(len(self.ids))
            row['seed_cell_id'] = pair_features(hits[j], hits[j], [energy[j]]*2,
                                                position[[j, j]])['seed_cell_id'][0]
            self.rows.append(row)
        self.event_offsets.append(len(self.rows))

    def finish(self, paths, output, coords, events, decays, report):
        out = Path(output)
        out.parent.mkdir(parents=True, exist_ok=True)
        if out.exists(): raise FileExistsError(out)
        # A zero-candidate shard is allowed and still retains all event identities.
        keys = self.rows[0].keys() if self.rows else []
        a = {k: np.array([r[k] for r in self.rows]) for k in keys}
        for records in [events, decays]:
            if records:
                a.update({k: np.array([r[k] for r in records]) for k in records[0]})
        cells = np.array(sorted(coords), dtype=np.uint64)
        a['cache_metadata'] = np.array(json.dumps(self.metadata))
        a.update(candidate_schema=np.array(1), candidate_floor_GeV=np.array(self.floor),
                 inputs=np.array(paths), cell_id=cells,
                 cell_position=np.array([coords[int(k)] for k in cells]),
                 event_offsets=np.array(self.event_offsets, dtype=np.int64),
                 hit_offsets=np.array(self.hit_offsets, dtype=np.int64),
                 hit_cell_id=np.array(self.ids, dtype=np.uint64),
                 hit_contribution=np.array(self.contributions, dtype=np.float64))
        if self.rows: a['seed_cell_id'] = a['seed_cell_id'].astype(np.uint64)
        else:
            a.update(energy=np.empty(0), position=np.empty((0,3)), vertex=np.empty((0,3)),
                     cluster_index=np.empty(0,dtype=np.int64), vertex_source=np.empty(0,dtype='U64'),
                     track_matched=np.empty(0,dtype=bool))
        if not self.metadata.get('truth_available', True):
            a = {k:v for k,v in a.items() if not k.startswith(('truth_', 'decay_'))}
        a['track_matching_available'] = np.array(self.metadata.get('track_matching_available', True))
        temp = out.with_name(out.stem+'.partial.npz')
        np.savez_compressed(temp, **a)
        temp.rename(out)
        report.update(candidate_schema=1, retained_candidates=len(self.rows),
                      events=len(events), selection=f'valid candidates E >= {self.floor} GeV; no pair choice or mass cut')
        out.with_suffix('.json').write_text(json.dumps(report, indent=2)+'\n')
        print(json.dumps(dict(output=str(out), candidates=len(self.rows), events=len(events))), flush=True)


def load(path):
    with np.load(path, allow_pickle=False) as f:
        a = {k: f[k] for k in f.files}
    assert int(a['candidate_schema']) == 1
    return a


def hit(a, k):
    lo, hi = a['hit_offsets'][k:k+2]
    return a['hit_cell_id'][lo:hi], a['hit_contribution'][lo:hi]


def truth_proxy(a, lo, hi):
    """Reconstruct just the EventTruth.pair inputs, preserving original MC IDs."""
    t = EventTruth.__new__(EventTruth)
    for attr, key in [('sim','truth_sim'), ('weight','truth_weight'), ('weight_sum','truth_weight_sum'),
                      ('origin','truth_photon'), ('parent','truth_parent'), ('track_matched','track_matched')]:
        setattr(t, attr, a[key][lo:hi])
    n = max(int(t.sim.max(initial=-1)), int(t.origin.max(initial=-1)))+1
    t.pdg = np.zeros(n, dtype=int); t.momentum = np.zeros((n, 3))
    for k in range(lo, hi):
        sim, origin = int(a['truth_sim'][k]), int(a['truth_photon'][k])
        if sim >= 0: t.pdg[sim] = a['truth_pdg'][k]
        if origin >= 0: t.momentum[origin] = a['truth_photon_momentum'][k]
    return t

