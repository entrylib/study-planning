"""Checks of the code against what is known exactly.

    python3 -m unittest test_rankings
"""
import math
import unittest

import numpy as np

import reference
import settings
from measures import MEASURES, TestSet


def brute(name, y, f):
    e = y[:, None] - f[None, :]
    if name == "RMSE":
        return np.mean(e ** 2, axis=0) - y.var()
    if name == "MAE":
        return np.mean(np.abs(e), axis=0)
    if name == "MAPE":
        return 100 * np.mean(np.abs(e) / y[:, None], axis=0)
    d = np.abs(y)[:, None] + np.abs(f)[None, :]
    return 200 * np.mean(np.where(e == 0, 0.0, np.abs(e) / np.where(d == 0, 1.0, d)), axis=0)


class TestMeasures(unittest.TestCase):
    def test_the_fast_measures_equal_their_definitions(self):
        g = np.random.default_rng(3)
        for y in (g.lognormal(4.6, 0.5, 37), g.poisson(1.5, 40).astype(float), g.poisson(0.3, 25).astype(float)):
            f = np.concatenate([[0.0], np.quantile(y, np.linspace(0, 1, 23)), [y.mean(), 1e-3]])
            t = TestSet(y)
            for m, fn in MEASURES.items():
                if m == "MAPE" and np.any(y == 0):
                    continue
                np.testing.assert_allclose(fn(t, f), brute(m, y, f), rtol=1e-9, atol=1e-9, err_msg=m)


class TestReference(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ref, reference.REF = reference.REF, 400_000

    @classmethod
    def tearDownClass(cls):
        reference.REF = cls.ref

    def test_lognormal_preferred_candidates_match_theory(self):
        s = 0.5
        c, _ = reference.continuous(lambda g, n: g.lognormal(math.log(100.0), s, n), 5)
        val = lambda m: c.values[c.preferred[m]]
        lev = lambda m: c.level[c.preferred[m]]
        self.assertAlmostEqual(val("RMSE"), 100 * math.exp(s * s / 2), delta=0.5)        # the mean
        self.assertAlmostEqual(val("MAE"), 100.0, delta=0.5)                            # the median
        self.assertAlmostEqual(val("MAPE"), 100 * math.exp(-s * s), delta=0.5)          # exp(mu - s^2)
        self.assertAlmostEqual(lev("MAPE"), 0.5 * (1 + math.erf(-s / math.sqrt(2))), delta=0.005)   # level Phi(-s)
        self.assertAlmostEqual(val("sMAPE"), 100.0, delta=0.5)   # sMAPE depends on log(y / f) only: the median

    def test_counts_are_exact(self):
        name, _, pmf = settings.counts()[0]                      # Poisson 1.5
        c = reference.counts(pmf)
        self.assertAlmostEqual(c.mean, 1.5, places=12)
        self.assertAlmostEqual(c.values[c.preferred["RMSE"]], 1.5, places=12)
        self.assertEqual(c.values[c.preferred["MAE"]], 1.0)
        self.assertAlmostEqual(c.zero_share, math.exp(-1.5), places=12)

    def test_a_zero_forecast_can_be_the_smape_minimiser(self):
        c = reference.counts(settings.counts()[8][2])            # intermittent Poisson 0.65: 52% zeros
        self.assertEqual(c.values[c.preferred["sMAPE"]], 0.0)


class TestSettings(unittest.TestCase):
    def test_continuous_outcomes_are_positive(self):
        g = np.random.default_rng(9)
        for name, s in settings.continuous():
            self.assertTrue(np.all(s(g, 200_000) > 0), name)

    def test_intermittent_settings_have_more_than_half_zeros(self):
        for name, _, pmf in settings.counts():
            if name.startswith("intermittent"):
                self.assertGreater(float(pmf(np.array([0]))[0]), 0.5, name)


if __name__ == "__main__":
    unittest.main()
