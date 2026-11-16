"""Checks of the oracle's code against what is known exactly.

In one test set, the candidate with the lowest RMSE is the sample mean and the candidate with the lowest MAE is the
sample median (squared error is minimised at the mean, absolute error at any median; a tie is resolved at the
midpoint, which for MAE is the usual sample median). Both scores are convex in the candidate, so the exact minimum
lies within one grid step of the best candidate on the grid: the grid search must return the sample mean and the
sample median to within one grid step. And the vectorised measures must equal their definitions computed one candidate at
a time.

    python3 -m unittest test_oracle
"""
import unittest

import numpy as np

from measures import MEASURES
from oracle import grid, select, settings


def brute(name, y, f):
    e = y[:, None] - f[None, :]
    if name == "RMSE":
        return np.sqrt(np.mean(e ** 2, axis=0))
    if name == "MAE":
        return np.mean(np.abs(e), axis=0)
    if name == "MAPE":
        keep = y != 0
        return 100 * np.mean(np.abs(e[keep]) / np.abs(y[keep])[:, None], axis=0)
    d = np.abs(y)[:, None] + np.abs(f)[None, :]
    return 200 * np.mean(np.where(e == 0, 0.0, np.abs(e) / np.where(d == 0, 1.0, d)), axis=0)


class TestSelection(unittest.TestCase):
    def test_rmse_selects_the_sample_mean_and_mae_the_sample_median(self):
        g = np.random.default_rng(2026)
        checked = 0
        for name, kind, t, draw in settings():
            for n in (10, 11, 30, 50, 100, 1000):            # odd and even sizes
                for _ in range(20):
                    y = draw(g, n).astype(float)
                    f, step = grid(kind, t)
                    for measure, target in (("RMSE", np.mean(y)), ("MAE", np.median(y))):
                        choice, edge = select(MEASURES[measure], y, f, step, kind)
                        self.assertFalse(edge)
                        self.assertLessEqual(abs(choice - target), step + 1e-9,
                                             f"{measure} on {name}, n = {n}: {choice} vs {target}")
                        checked += 1
        self.assertEqual(checked, 18 * 6 * 20 * 2)


class TestMeasures(unittest.TestCase):
    def test_the_vectorised_measures_equal_their_definitions(self):
        g = np.random.default_rng(7)
        for name, kind, t, draw in settings():
            y = draw(g, 37).astype(float)
            f, _ = grid(kind, t)
            for m, fn in MEASURES.items():
                if m == "MAPE" and kind != "continuous":
                    continue
                np.testing.assert_allclose(fn(y, f), brute(m, y, f), rtol=1e-9, atol=1e-12, err_msg=f"{m} on {name}")


if __name__ == "__main__":
    unittest.main()
