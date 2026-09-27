# SPDX-License-Identifier: LGPL-3.0-or-later
# Adapted from calibration/pi0/tests/test_model.py (toy/PhysicsTests).
"""Focused physics/numerics tests; run with Python's unittest, no pytest needed."""
import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np
from scipy import sparse
from eeemcal_pi0.model import Sample, M0, evaluate, equations, fit, angular_factor


def toy(n=800, seed=12):
    rng = np.random.default_rng(seed)
    k = 7  # last channel deliberately never illuminated
    energies = rng.uniform(.7, 8., (n, 2))
    contributions = []
    for j in range(2):
        c = np.zeros((n, k))
        for i in range(n):
            channels = rng.choice(k-1, 3, replace=False)
            c[i, channels] = rng.dirichlet(np.ones(3))*energies[i, j]/1.12
        contributions.append(c)
    angle = M0**2/np.prod(energies, axis=1)
    theta = np.arccos(1-angle/2)
    position = np.zeros((n, 2, 3))
    position[:, :, 2] = -1850
    position[:, 1, 0] = 1850*np.tan(theta)
    arrays = dict(energy=energies, angle_factor=angle, position=position,
                  vertex=np.zeros((n, 3)), cell_id=np.arange(k, dtype=np.uint64)+2**54)
    return Sample(arrays, tuple(sparse.csr_matrix(c) for c in contributions))


class PhysicsTests(unittest.TestCase):
    def test_derivative_symmetric_asymmetric_multicell_shared(self):
        a = dict(energy=np.array([[2., 2.], [.6, 9.], [3., 1.]]),
                 angle_factor=np.array([.003, .008, .01]), cell_id=np.arange(4, dtype=np.uint64))
        c1 = sparse.csr_matrix([[1., 1., 0., 0.], [.1, .4, .1, 0.], [.2, 0., 2.8, 0.]])
        c2 = sparse.csr_matrix([[0., 0., 1., 1.], [0., 4., 0., 5.], [0., .3, .5, .2]])
        sample = Sample(a, (c1, c2))
        gains = np.array([.9, 1.1, 1.2, .95])
        _, m2, J = evaluate(sample, gains, jacobian=True)
        np.testing.assert_allclose(np.asarray(J.sum(axis=1)).ravel(), 2*m2, rtol=2e-15)
        for k in range(4):
            plus, minus = gains.copy(), gains.copy()
            plus[k] *= 1+1e-6; minus[k] *= 1-1e-6
            fd = (evaluate(sample, plus)[1]-evaluate(sample, minus)[1])/2e-6
            np.testing.assert_allclose(fd, J[:, k].toarray().ravel(), rtol=1e-8, atol=1e-11)
        # Dropping fitted columns must not redistribute the frozen energy.
        C, D, L, B = equations(sample, gains, np.arange(3), np.array([0, 2]), chunk=2)
        np.testing.assert_allclose(L, np.asarray(J[:, [0, 2]].sum(axis=0)).ravel())

    def test_gain_recovery_without_injected_values_in_fit(self):
        nominal = toy()
        injected = np.array([.8, 1.2, .92, 1.12, .88, 1.05, 1.])
        energy, _ = evaluate(nominal, injected)
        distorted = Sample(dict(nominal.arrays, energy=energy),
                           tuple(c.multiply(injected).tocsr() for c in nominal.showers))
        result = fit(distorted, np.arange(distorted.n), min_usage=20, rcond=1e-8,
                     iterations=40, step_tolerance=1e-10, bias_tolerance=1e-13)
        gains, usage, active, history, reason = result
        np.testing.assert_allclose(gains[:6], 1/injected[:6], rtol=2e-7, atol=1e-8)
        self.assertEqual(gains[6], 1.)
        self.assertEqual(usage[6], 0)
        self.assertEqual(len(active), 6)
        self.assertEqual(reason, 'converged')
        self.assertLess(history[-1]['m2_residual_std_GeV2'], 1e-10)

    def test_unity_and_serialization_preserve_energy_geometry_cellids(self):
        sample = toy(n=12)
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)/'sample.npz'; sample.save(p)
            loaded = Sample.load(p)
        np.testing.assert_array_equal(loaded.arrays['cell_id'], sample.arrays['cell_id'])
        self.assertEqual(loaded.arrays['cell_id'].dtype, np.uint64)
        energy, m2 = evaluate(loaded, np.ones(loaded.k))
        np.testing.assert_allclose(energy, sample.arrays['energy'], rtol=2e-15)
        np.testing.assert_allclose(m2, M0**2, rtol=3e-15)
        np.testing.assert_allclose(angular_factor(sample.arrays['position'], sample.arrays['vertex']),
                                   sample.arrays['angle_factor'], rtol=1e-10)

    def test_chunk_and_shard_equations_match(self):
        sample = toy(n=37)
        gains = np.array([.8, 1.2, .92, 1.12, .88, 1.05, 1.])
        rows, active = np.arange(sample.n), np.arange(6)
        full = equations(sample, gains, rows, active, chunk=100)
        chunked = equations(sample, gains, rows, active, chunk=4)
        shards = [equations(sample, gains, r, active, chunk=5) for r in [rows[::2], rows[1::2]]]
        for a, b, c, d in zip(full, chunked, *shards):
            np.testing.assert_allclose(a, b, rtol=2e-14, atol=1e-16)
            np.testing.assert_allclose(a, c+d, rtol=2e-14, atol=1e-16)


if __name__ == '__main__':
    unittest.main(verbosity=2)
