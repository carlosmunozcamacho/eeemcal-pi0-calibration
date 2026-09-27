# SPDX-License-Identifier: LGPL-3.0-or-later
# Adapted from calibration/pi0/src/extract.py. Axel's C++ prototype informed the earlier
# topology/mass workflow; see THIRD_PARTY_NOTICES.md for the distinct relation-reader provenance.
"""Read checked PODIO relations into broad candidate shards.

Assigned energies are used directly. Optional truth is separate from required
track matching; no pair ranking or mass selection occurs during extraction.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
from itertools import combinations
import awkward as ak
import numpy as np
from scipy import sparse
import uproot
from .model import Sample, angular_factor, evaluate

CL = 'EcalEndcapNClusters'
RH = 'EcalEndcapNRecHits'


def leaf(collection, member):
    return f'{collection}/{collection}.{member}'


def collection_map(f):
    m = f['podio_metadata']
    names = ak.to_list(m['events___CollectionTypeInfo.name'].array())[0]
    ids = ak.to_list(m['events___CollectionTypeInfo.collectionID'].array())[0]
    return dict(zip(map(int, ids), names))


def open_input(path):
    return uproot.open(path)


def select_vertex(a, mapping, vertices, cfg):
    """The original pilot vertex rule, shared by baseline and all-pair mode."""
    vertex, source = np.asarray(cfg['nominal_mm'], dtype=float), 'nominal_IP'
    choices = []
    for idx, cid in zip(a.get('vidx', []), a.get('vcid', [])):
        name = mapping.get(int(cid))
        if name not in vertices or idx < 0 or idx >= len(a[name+'.ndf']):
            raise ValueError('Invalid primary vertex relation')
        v = np.array([a[name+'.position.'+d][idx] for d in 'xyz'])
        ndf, chi2 = a[name+'.ndf'][idx], a[name+'.chi2'][idx]
        if np.all(np.isfinite(v)) and ndf > 0 and np.isfinite(chi2) and chi2 >= 0 and np.all(np.abs(v[:2]) < cfg['max_xy_mm']) and abs(v[2]) < cfg['max_z_mm']:
            choices.append((-ndf, chi2/ndf, idx, name, v))
    if choices:
        best = sorted(choices, key=lambda x: x[:4])[0]
        vertex, source = best[-1], best[-2]
    if source == 'nominal_IP' and not cfg['allow_nominal']:
        raise ValueError('No usable reconstructed vertex and nominal fallback disabled')
    return vertex, source


def extract(paths, output, config, candidate_sink):
    chunk, limit, diagnostics = config['input']['chunk'], None, True
    CL, RH = config['input']['clusters'], config['input']['hits']
    from .truth import RecoOnly, truth_setup
    from .config import partition
    if candidate_sink is not None and not diagnostics:
        raise ValueError('Candidate output requires diagnostics')
    identities = json.loads(Path(config['input']['identity_map']).read_text()) if config['input']['identity_map'] else {}
    canonical = [identities.get(path,path) for path in paths]
    if len(set(canonical)) != len(canonical): raise ValueError('Duplicate canonical source identities')
    start = time.monotonic()
    report = {'files': [], 'exactly_two': 0, 'accepted': 0, 'invalid_pair': 0,
              'development_limit_per_file': limit, 'diagnostics_enabled': diagnostics,
              'clusters_checked': 0, 'shared_hit_occurrences': 0,
              'relation_errors': 0, 'nonconserving_hit_assignments': 0,
              'max_assigned_to_hit_ratio': 0., 'max_contribution_to_hit_ratio': 0.,
              'vertex_sources': {}, 'vertex_rule': 'PrimaryVertices subset, ndf>0, finite chi2>=0, |x,y|<20 mm, |z|<300 mm; highest ndf then lowest chi2/ndf; otherwise nominal IP',
              'selection': 'broad valid candidate cache, no pair choice or mass cut'}
    rows, coords, ratios = [], {}, []
    pair_extras, event_rows, decay_rows = [], [], []
    if diagnostics:
        from .diagnostics import pair_features
        from .mc_diagnostics import EventTruth, fields as truth_fields
        report.update(statistics_schema=1, selection='broad candidates; no pair choice or mass cut')
    payload = [[], []]
    for fi, path in enumerate(paths):
        with open_input(path) as f:
            t = f['events']
            mapping = collection_map(f)
            vertices = [n for n in mapping.values() if leaf(n, 'ndf') in t and leaf(n, 'position.z') in t]
            fields = {f'cl_{m}': leaf(CL, m) for m in ['energy', 'position.x', 'position.y', 'position.z', 'nhits', 'hits_begin', 'hits_end', 'hitContributions_begin', 'hitContributions_end']}
            fields.update({f'rh_{m}': leaf(RH, m) for m in ['energy', 'cellID', 'position.x', 'position.y', 'position.z']})
            fields.update(idx=leaf('_'+CL+'_hits', 'index'), cid=leaf('_'+CL+'_hits', 'collectionID'),
                          contrib='_'+CL+'_hitContributions', event=leaf(config['input']['event_header'], 'eventNumber'), run=leaf(config['input']['event_header'], 'runNumber'))
            primary = config['input']['primary_vertices']+'_objIdx'
            if primary in t:
                fields.update(vidx=leaf(primary, 'index'), vcid=leaf(primary, 'collectionID'))
            for name in vertices:
                for member in ['position.x', 'position.y', 'position.z', 'ndf', 'chi2']:
                    fields[name+'.'+member] = leaf(name, member)
            extra, truth_enabled, track_available, canonical_mapping = truth_setup(t, mapping, config)
            fields.update(extra)
            candidate_sink.metadata['truth_available'] = candidate_sink.metadata.get('truth_available', True) and truth_enabled
            candidate_sink.metadata['track_matching_available'] = candidate_sink.metadata.get('track_matching_available', True) and track_available
            missing = [p for p in fields.values() if p not in t]
            if missing:
                raise ValueError(f'Missing required persisted information: {missing}')
            total = min(t.num_entries, limit) if limit else t.num_entries
            info = {'path': canonical[fi], 'reader_path': path, 'events_in_file': t.num_entries, 'events_read': total,
                    'size_bytes': Path(path).stat().st_size if Path(path).is_file() else None,
                    'collections': {n: n in mapping.values() for n in [CL, RH, 'EcalEndcapNRawHits', 'EcalEndcapNHits', 'EcalEndcapNProtoClusters', 'MCParticles', 'PrimaryVertices']},
                    'collection_ids': {str(k): v for k, v in mapping.items() if 'EcalEndcapN' in v or 'Vertic' in v}}
            report['files'].append(info)
            before = len(rows)
            for begin in range(0, total, chunk):
                # Separate explicit branches avoid ambiguous nested branch aliases.
                arrays = {alias: ak.to_list(t[p].array(entry_start=begin, entry_stop=min(begin+chunk, total))) for alias, p in fields.items()}
                for ev in range(len(arrays['cl_energy'])):
                    a = {k: v[ev] for k, v in arrays.items()}
                    if not track_available:
                        for rel in ['cluster','track']:
                            for member in ['index','collectionID']: a['track_'+rel+'_'+member]=[]
                    hit_ids = np.asarray(a['rh_cellID'], dtype=np.uint64)
                    hit_e = np.asarray(a['rh_energy'], dtype=float)
                    xyz = np.column_stack([a['rh_position.'+d] for d in 'xyz'])
                    for cell, pos in zip(hit_ids, xyz):
                        cell = int(cell)
                        if cell in coords and not np.allclose(coords[cell], pos, rtol=0, atol=.01):
                            raise ValueError(f'Inconsistent persisted position for cell {cell}')
                        coords[cell] = pos
                    cluster_hits, assigned = [], np.zeros(len(hit_e))
                    hit_uses = np.zeros(len(hit_e), dtype=int)
                    for j, energy in enumerate(a['cl_energy']):
                        hb, he = a['cl_hits_begin'][j], a['cl_hits_end'][j]
                        cb, ce = a['cl_hitContributions_begin'][j], a['cl_hitContributions_end'][j]
                        if not (0 <= hb <= he <= len(a['idx']) == len(a['cid']) and 0 <= cb <= ce <= len(a['contrib'])):
                            raise ValueError('Cluster relation offsets out of bounds')
                        idx = np.asarray(a['idx'][hb:he], dtype=int)
                        ids = a['cid'][hb:he]
                        c = np.asarray(a['contrib'][cb:ce], dtype=float)
                        if not (len(idx) == len(c) == a['cl_nhits'][j]) or any(mapping.get(int(cid)) != RH for cid in ids):
                            raise ValueError('Cluster hit/contribution length or target collection mismatch')
                        if np.any(idx < 0) or np.any(idx >= len(hit_e)) or np.any(c < 0) or not np.all(np.isfinite(c)):
                            raise ValueError('Invalid event-local hit relation or assigned energy')
                        if np.any(c > hit_e[idx] * (1+2e-5) + 1e-8):
                            raise ValueError('Assigned contribution exceeds linked hit energy')
                        np.add.at(assigned, idx, c)
                        np.add.at(hit_uses, idx, 1)
                        if len(c):
                            report['max_contribution_to_hit_ratio'] = max(report['max_contribution_to_hit_ratio'], float(np.max(c/np.maximum(hit_e[idx], 1e-30))))
                        if c.sum() > 0:
                            ratios.append(float(energy/c.sum()))
                        cluster_hits.append((hit_ids[idx], c))
                        report['clusters_checked'] += 1
                    report['shared_hit_occurrences'] += int(np.sum(hit_uses > 1))
                    report['nonconserving_hit_assignments'] += int(np.sum(assigned > hit_e*(1+2e-5)+1e-8))
                    if len(hit_e):
                        report['max_assigned_to_hit_ratio'] = max(report['max_assigned_to_hit_ratio'], float(np.max(assigned/np.maximum(hit_e, 1e-30))))
                    multiplicity = len(cluster_hits)
                    if multiplicity == 2:
                        report['exactly_two'] += 1
                    if not diagnostics and multiplicity != 2:
                        continue
                    pos = np.column_stack([a['cl_position.'+d] for d in 'xyz'])
                    energy = np.asarray(a['cl_energy'], dtype=float)
                    vertex, source = select_vertex(a, mapping, vertices, config['vertex'])
                    valid = (np.all(np.isfinite(pos), axis=1) & np.isfinite(energy) & (energy > 0)
                             & np.array([c.sum() > 0 for _, c in cluster_hits], dtype=bool)
                             & (np.linalg.norm(pos-vertex, axis=1) > 0))
                    entry = begin+ev
                    # Fixed deterministic split, assigned before any peak/window tuning.
                    train = partition(canonical[fi], entry, config['split'])
                    counts = {f'candidate_count_{tag}': int(np.sum(valid & (energy >= threshold)))
                              for tag, threshold in [('020', .2), ('030', .3), ('050', .5)]}
                    if diagnostics:
                        truth = EventTruth(a, canonical_mapping, energy) if truth_enabled else RecoOnly(a, canonical_mapping, energy)
                        counts['neutral_count_030'] = int(np.sum(valid & (energy >= .3) & ~truth.track_matched))
                        decays = []  # Optional per-candidate truth is sufficient; no generator census.
                        decay_rows.extend(decays)
                        event_rows.append(dict(event_file_index=fi, event_entry=entry, event_run=a['run'][0],
                                               event_number=a['event'][0], event_train=train,
                                               event_ncluster=multiplicity, event_nvalid=int(valid.sum()),
                                               event_vertex_source=source,
                                               **{'event_'+k: v for k, v in counts.items()}))
                    if multiplicity == 2 and not valid.all():
                        report['invalid_pair'] += 1
                    if candidate_sink is not None:
                        candidate_sink.add(len(event_rows)-1, energy, pos, vertex, source,
                                           valid, cluster_hits, truth)
                        continue
            info['accepted_pairs'] = len(rows)-before
            print(json.dumps(info), flush=True)
    if report['nonconserving_hit_assignments']:
        raise ValueError('Nonconserving assigned hit energies')
    if candidate_sink is not None:
        candidate_sink.finish(canonical, output, coords, event_rows, decay_rows, report)
        return
