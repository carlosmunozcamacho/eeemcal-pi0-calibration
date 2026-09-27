# SPDX-License-Identifier: LGPL-3.0-or-later
# Adapted from calibration/pi0/src/diagnostics.py (pair_features/event_ids/meaningful_usage).
"""Reusable pair diagnostics and exact-cellID merging for statistics shards.

No truth or mass is used by the data-like selection definitions below.
Energy thresholds are trial settings, not detector-qualified photon cuts.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy import sparse
from .model import Sample


def pair_features(hit1, hit2, energy, position):
    """Sharing uses assigned energies, with repeated cellIDs summed exactly."""
    def by_cell(hit):
        result = {}
        for cell, value in zip(*hit):
            result[int(cell)] = result.get(int(cell), 0.) + float(value)
        return result
    c1, c2 = by_cell(hit1), by_cell(hit2)
    common = c1.keys() & c2.keys()
    overlap = sum(min(c1[k], c2[k]) for k in common) / min(sum(c1.values()), sum(c2.values()))
    seed1, seed2 = max(c1, key=c1.get), max(c2, key=c2.get)
    return dict(overlap=overlap, shared_cells=sum(c1[k] > 0 and c2[k] > 0 for k in common),
                same_seed=seed1 == seed2, seed_cell_id=np.array([seed1, seed2], dtype=np.uint64),
                separation_mm=float(np.linalg.norm(position[0]-position[1])),
                asymmetry=float(abs(energy[0]-energy[1])/sum(energy)))


def event_ids(sample, rows):
    a = sample.arrays
    return np.rec.fromarrays([a['file_index'][rows], a['entry'][rows]], names='file,entry')


def meaningful_usage(sample, rows, min_fraction=.05):
    """Independent-event support from >=5% of either original shower energy.

    Tiny shared tails cannot qualify a channel. This support definition changes
    only which channels may be fitted; all cells stay in model denominators.
    """
    rows = np.asarray(rows, dtype=int)
    fractions = []
    for c in sample.showers:
        c = c[rows]
        total = np.asarray(c.sum(axis=1)).ravel()
        f = c.multiply((1 / total)[:, None]).tocsr()
        f.data = (f.data >= min_fraction).astype(float)
        f.eliminate_zeros()
        fractions.append(f)
    present = (fractions[0] + fractions[1]).tocoo()
    ids = event_ids(sample, rows)
    _, event = np.unique(ids, return_inverse=True)
    unique = sparse.coo_matrix((np.ones(len(present.data)), (event[present.row], present.col)),
                               shape=(int(event.max())+1 if len(event) else 0, sample.k)).tocsr()
    unique.data[:] = 1
    return np.asarray(unique.sum(axis=0)).ravel().astype(int)

