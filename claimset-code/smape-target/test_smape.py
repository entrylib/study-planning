"""Checks of the code against exact results.

    python3 -m unittest test_smape
"""
import math
import unittest

import numpy as np

import dists
import minimiser
import theory
from measures import MEASURES, TestSet


def brute(name, y, f):
    e = y[:, None] - f[None, :]
    if name == "RMSE":
        return np.mean(e ** 2, axis=0) - y.var()
    if name == "MAE":
        return np.mean(np.abs(e), axis=0)
    if name == "MAPE":
        return 100 * np.mean(np.abs(e) / np.abs(y)[:, None], axis=0)
    d = np.abs(y)[:, None] + np.abs(f)[None, :]
    return 200 * np.mean(np.where(e == 0, 0.0, np.abs(e) / np.where(d == 0, 1.0, d)), axis=0)


class TestMeasures(unittest.TestCase):
    def test_fast_measures_equal_their_definitions(self):
        g = np.random.default_rng(3)
        for y in (g.lognormal(4.6, 0.5, 37), g.normal(0.5, 1, 40), np.where(g.random(30) < 0.4, 0.0, g.lognormal(4.6, 1, 30))):
            f = np.concatenate([[0.0], np.quantile(y, np.linspace(0, 1, 23)), [y.mean(), 1e-3, -0.7]])
            t = TestSet(y)
            for m, fn in MEASURES.items():
                if m == "MAPE" and np.any(y == 0):
                    self.assertTrue(np.all(np.isinf(fn(t, f))))
                    continue
                np.testing.assert_allclose(fn(t, f), brute(m, y, f), rtol=1e-9, atol=1e-9, err_msg=m)


class TestMinimiser(unittest.TestCase):
    def test_the_minimiser_is_an_outcome_and_the_search_finds_it(self):
        g = np.random.default_rng(7)
        for n, s in ((7, 1.0), (40, 3.0), (600, 0.5), (1500, 4.0)):
            x = np.sort(g.normal(0, s, (6, n)), axis=1)
            found = minimiser.log_minimisers(x)
            for i in range(6):
                fine = np.linspace(x[i, 0], x[i, -1], 20001)                 # no forecast between outcomes does better
                s_fine = np.tanh(np.abs(x[i][None, :] - fine[:, None]) / 2).sum(axis=1).min()
                s_found = np.tanh(np.abs(x[i] - found[i]) / 2).sum()
                self.assertLessEqual(s_found, s_fine + 1e-9)
                self.assertIn(found[i], x[i])
            if n > minimiser.EXACT_UP_TO:
                np.testing.assert_array_equal(found, minimiser.exact(x))


class TestTheory(unittest.TestCase):
    def test_moments_against_quadrature(self):
        for d in [dists.lognormal(1.0), dists.log_laplace(0.25), dists.log_skewed(4.0, 0.5, 1), dists.log_skewed(1.0, 0.5, -1),
                  dists.bimodal(0.35), dists.inverse_gaussian(20.0), dists.log_gumbel(0.2, 1), dists.log_gumbel(0.2, -1)]:
            gr = theory.Grid(d)
            self.assertAlmostEqual(gr.mass, 1.0, places=6, msg=d.name)
            for m, p in (("mean", 1), ("second", 2), ("inv", -1)):
                if d.finite[m]:
                    self.assertAlmostEqual(float(np.exp(p * gr.x) @ gr.gw) / d.moments[m], 1.0, places=5, msg=(d.name, m))

    def test_samplers_match_densities(self):
        g = np.random.default_rng(11)
        levels = np.array([0.05, 0.25, 0.5, 0.75, 0.95])
        for d in dists.pl1()[:-11] + [dists.log_gumbel(2.0, 1), dists.log_laplace(1.5)]:
            T = theory.theory(d)
            y = d.sample(g, 400000)
            np.testing.assert_allclose(T.level(np.quantile(y, levels)), levels, atol=4e-3, err_msg=d.name)

    def test_fft_equals_direct(self):
        T = theory.theory(dists.log_skewed(1.0, 1.0, 1))
        gr = T.grid
        for j in (100000, 250000, 400000):
            self.assertAlmostEqual(T.S_plus[j], gr.smape_at(gr.x[j])[0], places=7)
            self.assertAlmostEqual(T.S_plus[j], T.smape_exact(gr.x[j]), places=5)

    def test_log_symmetric_targets_are_the_median(self):           # H5
        for d in (dists.lognormal(0.3), dists.lognormal(2.0), dists.log_laplace(0.8)):
            T = theory.theory(d)
            self.assertEqual(len(T.minima_plus), 1)
            self.assertAlmostEqual(T.minima_plus[0][0], math.log(100), places=4)

    def test_mape_minimiser_of_a_lognormal(self):                   # the 1/y-weighted lognormal has median e^(mu - s^2)
        T = theory.theory(dists.lognormal(0.5))
        self.assertAlmostEqual(T.mape_minimiser() / (100 * math.exp(-0.25)), 1.0, places=4)
        self.assertEqual(theory.theory(dists.log_skewed(1.0, 1.0, -1)).mape_minimiser(), 0.0)    # E[1/Y] infinite

    def test_zero_threshold_balances(self):                         # H8
        base = dists.lognormal(0.5)
        T = theory.theory(base)
        p = T.zero_threshold()
        Tz = theory.theory(dists.ZeroInflated(base, p))
        at = Tz.smape_minima()
        self.assertAlmostEqual(at[0][1], at[1][1], places=6)
        self.assertLess(p, 0.5)

    def test_bimodal_equal_modes_has_two_equal_minima(self):        # H7
        T = theory.theory(dists.bimodal(0.5))
        self.assertEqual(len(T.minima_plus), 2)
        self.assertAlmostEqual(T.minima_plus[0][1], T.minima_plus[1][1], places=5)
        self.assertAlmostEqual(T.minima_plus[0][0] + T.minima_plus[1][0], 2 * math.log(100), places=4)

    def test_flat_approximation_improves_with_spread(self):          # H14
        errs = []
        for s in (2.0, 6.0, 20.0):
            T = theory.theory(dists.log_gumbel(s, 1))
            c, v = T.minima_plus[0]
            errs.append(abs(v - T.flat_approximation(c)))
        self.assertTrue(errs[0] > errs[1] > errs[2])
        self.assertLess(errs[2], 0.05)

    def test_asymptotic_relative_efficiency(self):                  # H11, H14 (iii)
        v_s, v_m, _ = theory.theory(dists.lognormal(1.0)).asymptotic_variances()
        self.assertAlmostEqual(v_m / v_s, 0.660, places=2)
        v_s, v_m, _ = theory.theory(dists.lognormal(20.0)).asymptotic_variances()
        self.assertAlmostEqual(v_m / v_s * 20 ** 3, 7.1, delta=0.15)

    def test_infinite_moments(self):                                # H12
        T = theory.theory(dists.log_skewed(1.0, 1.0, 1))            # E Y infinite
        e = T.expected(np.array([50.0, 100.0]))
        self.assertIsNone(e["RMSE"][1])
        self.assertTrue(np.all(np.isinf(e["MAE"][0])) and np.all(np.isfinite(e["MAE"][1])))
        T = theory.theory(dists.log_skewed(4.0, 1.0, 1))            # E Y finite, E Y^2 infinite
        e = T.expected(np.array([50.0, T.mean_plus, 300.0]))
        self.assertTrue(np.all(np.isinf(e["RMSE"][0])))
        self.assertEqual(int(np.argmin(e["RMSE"][1])), 1)


if __name__ == "__main__":
    unittest.main()
