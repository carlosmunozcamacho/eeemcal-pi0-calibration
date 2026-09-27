# SPDX-License-Identifier: LGPL-3.0-or-later
# Adapted from calibration/pi0/src/mc_diagnostics.py (fields/EventTruth).
"""Truth controls from persisted ObjectIDs, never inputs to data-like gains.

The 25.10.0 file links EcalEndcapNClusterAssociations to MCParticles and
EcalEndcapNClusters. Parent/daughter ObjectIDs target MCParticles. Stored
track-match weights can be zero: existence of the relation defines a match.
"""
import numpy as np


def fields():
    def leaf(c, m): return f'{c}/{c}.{m}'
    result = {'mc_'+m: leaf('MCParticles', m) for m in
              ['PDG', 'generatorStatus', 'parents_begin', 'parents_end', 'daughters_begin', 'daughters_end',
               'momentum.x', 'momentum.y', 'momentum.z', 'vertex.x', 'vertex.y', 'vertex.z']}
    for relation in ['parents', 'daughters']:
        for member in ['index', 'collectionID']:
            result['mc_'+relation+'_'+member] = leaf('_MCParticles_'+relation, member)
    for member in ['recID', 'simID', 'weight']:
        result['assoc_'+member] = leaf('EcalEndcapNClusterAssociations', member)
    for relation in ['rec', 'sim']:
        for member in ['index', 'collectionID']:
            result['assoc_'+relation+'_'+member] = leaf('_EcalEndcapNClusterAssociations_'+relation, member)
    for relation in ['cluster', 'track']:
        for member in ['index', 'collectionID']:
            result['track_'+relation+'_'+member] = leaf('_EcalEndcapNTrackClusterMatches_'+relation, member)
    return result


class EventTruth:
    def __init__(self, a, mapping, energy):
        self.a, self.mapping = a, mapping
        self.pdg = np.asarray(a['mc_PDG'], dtype=int)
        self.status = np.asarray(a['mc_generatorStatus'], dtype=int)
        n, nc = len(self.pdg), len(energy)
        self.momentum = np.column_stack([a['mc_momentum.'+d] for d in 'xyz'])
        self.vertex = np.column_stack([a['mc_vertex.'+d] for d in 'xyz'])
        self.relations = {}
        for rel in ['parents', 'daughters']:
            idx = np.asarray(a['mc_'+rel+'_index'], dtype=int)
            cid = a['mc_'+rel+'_collectionID']
            begins, ends = np.asarray(a['mc_'+rel+'_begin']), np.asarray(a['mc_'+rel+'_end'])
            if (len(idx) != len(cid) or len(begins) != n or len(ends) != n or
                    np.any(idx < 0) or np.any(idx >= n) or np.any(begins > ends) or np.any(ends > len(idx)) or
                    any(mapping.get(int(c)) != 'MCParticles' for c in set(cid))):
                raise ValueError('Invalid MC ancestry ObjectID/offset')
            self.relations[rel] = (idx, begins, ends)
        generated = np.flatnonzero((self.pdg == 111) & (self.status > 0))
        self.n_pi0 = len(generated)
        self.decays, self.photon_parent = {}, {}
        for parent in generated:
            daughters = self.related('daughters', parent)
            # Direct reciprocal two-photon decay, using emitted photons.
            if len(daughters) != 2 or len(set(daughters)) != 2:
                continue
            if not all(self.pdg[d] == 22 and self.status[d] > 0 and parent in self.related('parents', d) for d in daughters):
                continue
            self.decays[int(parent)] = tuple(map(int, daughters))
            for photon in daughters: self.photon_parent[int(photon)] = int(parent)
        weights = [{} for _ in energy]
        count = len(a['assoc_weight'])
        keys = ['assoc_rec_index', 'assoc_sim_index', 'assoc_rec_collectionID', 'assoc_sim_collectionID', 'assoc_recID', 'assoc_simID']
        if any(len(a[k]) != count for k in keys):
            raise ValueError('Inconsistent cluster-MC association lengths')
        for z in range(count):
            rec, sim, weight = int(a['assoc_rec_index'][z]), int(a['assoc_sim_index'][z]), float(a['assoc_weight'][z])
            if not (0 <= rec < nc and 0 <= sim < n and np.isfinite(weight) and weight >= 0 and
                    mapping.get(int(a['assoc_rec_collectionID'][z])) == 'EcalEndcapNClusters' and
                    mapping.get(int(a['assoc_sim_collectionID'][z])) == 'MCParticles' and
                    rec == a['assoc_recID'][z] and sim == a['assoc_simID'][z]):
                raise ValueError('Invalid cluster-MC association or redundant ID mismatch')
            weights[rec][sim] = weights[rec].get(sim, 0.) + weight
        self.sim = np.full(nc, -1, dtype=int)
        self.weight = np.zeros(nc)
        self.weight_sum = np.zeros(nc)
        self.origin = np.full(nc, -1, dtype=int)
        self.parent = np.full(nc, -1, dtype=int)
        for i, w in enumerate(weights):
            total = sum(w.values()); self.weight_sum[i] = total
            if total <= 0: continue
            sim = max(w, key=lambda k: (w[k], -k))
            self.sim[i], self.weight[i] = sim, w[sim]/total
            origins = self.photon_ancestors(sim)
            if len(origins) == 1:
                self.origin[i] = next(iter(origins))
                self.parent[i] = self.photon_parent[self.origin[i]]
            elif len(origins) > 1:
                self.origin[i] = -2  # Explicit ambiguous ancestry.
        self.track_matched = np.zeros(nc, dtype=bool)
        ti, tc = a['track_cluster_index'], a['track_cluster_collectionID']
        if len(ti) != len(tc) or len(ti) != len(a['track_track_index']) or len(ti) != len(a['track_track_collectionID']):
            raise ValueError('Invalid track-match relation lengths')
        for i, c, track, track_c in zip(ti, tc, a['track_track_index'], a['track_track_collectionID']):
            if not (0 <= i < nc and mapping.get(int(c)) == 'EcalEndcapNClusters' and track >= 0 and
                    mapping.get(int(track_c)) == 'CentralCKFTracks'):
                raise ValueError('Unexpected persisted track-match target')
            self.track_matched[i] = True

    def related(self, relation, i):
        idx, begin, end = self.relations[relation]
        return idx[begin[i]:end[i]]

    def photon_ancestors(self, i):
        found, seen, pending = set(), set(), [int(i)]
        while pending:
            k = pending.pop()
            if k in seen: continue
            seen.add(k)
            if k in self.photon_parent:
                found.add(k)
            else:
                pending.extend(map(int, self.related('parents', k)))
        return found

    def pair(self, i, j):
        indices = [i, j]
        sim = self.sim[indices]
        origin, parent = self.origin[indices], self.parent[indices]
        pdg = np.array([self.pdg[k] if k >= 0 else 0 for k in sim])
        weight = self.weight[indices]
        truth_m2 = -1.
        if np.any(sim < 0): label = 'unmatched'
        elif np.any(weight < .8) or np.any(origin == -2): label = 'ambiguous'
        elif origin[0] >= 0 and origin[1] >= 0 and origin[0] != origin[1] and parent[0] == parent[1]:
            label = 'distinct_photons_same_pi0'
            p = self.momentum[origin]
            e = np.linalg.norm(p, axis=1)
            truth_m2 = float(np.prod(e)*np.sum((p[0]/e[0]-p[1]/e[1])**2))
        elif origin[0] >= 0 and origin[0] == origin[1]: label = 'same_photon_shower'
        elif sim[0] == sim[1]:
            label = 'same_electron_shower' if abs(pdg[0]) == 11 else ('same_photon_shower' if pdg[0] == 22 else 'same_hadron_shower')
        elif np.any(np.abs(pdg) > 100): label = 'hadron_background'
        elif np.any(np.abs(pdg) == 11): label = 'electron_background'
        else: label = 'photon_or_other_combinatorial'
        return dict(truth_label=label, truth_sim=sim, truth_pdg=pdg, truth_weight=weight,
                    truth_weight_sum=self.weight_sum[indices], truth_photon=origin,
                    truth_parent=parent, truth_m2=truth_m2, track_matched=self.track_matched[indices])

    def decay_records(self, fi, entry, train, energy, position, vertex, valid):
        from .model import angular_factor
        records = []
        for parent, photons in self.decays.items():
            photons = np.asarray(photons, dtype=int)
            p, v = self.momentum[photons], self.vertex[photons]
            e = np.linalg.norm(p, axis=1)
            if np.any(e <= 0) or not np.all(np.isfinite(p)):
                raise ValueError('Invalid generated decay-photon momentum')
            m2 = float(np.prod(e)*np.sum((p[0]/e[0]-p[1]/e[1])**2))
            clusters = [np.flatnonzero(valid & (self.origin == photon) & (self.weight >= .8)) for photon in photons]
            best = [max(c, key=lambda k: (energy[k], -k)) if len(c) else -1 for c in clusters]
            reconstructed = all(i >= 0 for i in best) and best[0] != best[1]
            reco_m2 = (float(np.prod(energy[best])*angular_factor(position[best][None], vertex[None])[0])
                       if reconstructed else -1.)
            records.append(dict(decay_file_index=fi, decay_entry=entry, decay_train=train,
                                decay_parent=parent, decay_photon=photons, decay_momentum=p, decay_vertex=v,
                                decay_truth_m2=m2, decay_nclusters=[len(c) for c in clusters],
                                decay_best_cluster=best, decay_reco_m2=reco_m2))
        return records
