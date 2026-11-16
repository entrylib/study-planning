"""The candidate forecasts of one outcome distribution, fixed before any test set is drawn, and their expected scores.

Continuous distributions: everything comes from one reference sample of REF outcomes, drawn once with its own seed
(the supervisor's decision DC2: "whatever works to get quantiles"). Counts: everything is exact, from the probability
mass function.

Candidates: the quantiles at levels 0.01, 0.02, ..., 0.99 (the 0.50 quantile is the median), the mean, and the values
that minimise the expected MAPE (continuous only) and the expected sMAPE. Candidates with equal values are merged into
one, keeping all their labels."""
import math

import numpy as np

from measures import MEASURES, TestSet

REF = 10_000_000
LEVELS = np.round(np.arange(1, 100) / 100, 2)


def _golden(fn, lo, hi, tol):
    phi = (math.sqrt(5) - 1) / 2
    a, b = lo, hi
    c, d = b - phi * (b - a), a + phi * (b - a)
    fc, fd = fn(c), fn(d)
    while b - a > tol:
        if fc <= fd:
            b, d, fd = d, c, fc
            c = b - phi * (b - a)
            fc = fn(c)
        else:
            a, c, fc = c, d, fd
            d = a + phi * (b - a)
            fd = fn(d)
    return (a + b) / 2


def _minimise(score, grid, tol):
    """Minimise a one-dimensional expected score: the best point of a grid, then golden-section search between its
    neighbours."""
    s = score(grid)
    k = int(np.argmin(s))
    lo, hi = grid[max(k - 1, 0)], grid[min(k + 1, grid.size - 1)]
    x = _golden(lambda x: float(score(np.array([x]))[0]), lo, hi, tol)
    # the grid's best point may be lower still: sMAPE jumps at a forecast of exactly 0 (a zero outcome scores 0 there
    # and 200 just above it), which a search of the smooth part cannot find
    return x if float(score(np.array([x]))[0]) < s[k] else float(grid[k])


class Candidates:
    """labels: list of label lists; values: array; level: the share of the distribution at or below each value."""

    def __init__(self, labels, values, level, expected, preferred, mean, median, zero_share=None):
        self.labels, self.values, self.level = labels, values, level
        self.expected, self.preferred, self.mean, self.median = expected, preferred, mean, median
        self.zero_share = zero_share


def _merge(pairs):
    pairs = sorted(pairs, key=lambda p: p[1])
    labels, values = [], []
    for lab, v in pairs:
        if values and abs(v - values[-1]) <= 1e-6 * max(1.0, abs(v)):   # a minimiser found to within the search
            labels[-1].append(lab)
        else:
            labels.append([lab])
            values.append(v)
    return labels, np.array(values)


def continuous(sampler, seed, mape_infinite=False):
    """mape_infinite: the distribution has positive density at zero, so the expected MAPE of every candidate is
    infinite; it is recorded as inf and MAPE has no minimiser and no preferred candidate."""
    g = np.random.default_rng(seed)
    ref = TestSet(sampler(g, REF))
    q = np.quantile(ref.y, LEVELS)
    pairs = [(f"q{p:.2f}", float(v)) for p, v in zip(LEVELS, q)]
    pairs[49] = ("median (q0.50)", pairs[49][1])
    pairs.append(("mean", ref.mean))
    tol = (q[-1] - q[0]) * 1e-7                             # the search grid is the 99 quantiles themselves
    if not mape_infinite:
        pairs.append(("MAPE minimiser", _minimise(lambda f: MEASURES["MAPE"](ref, f), q, tol)))
    pairs.append(("sMAPE minimiser", _minimise(lambda f: MEASURES["sMAPE"](ref, f, chunk=1), q, tol)))
    labels, values = _merge(pairs)
    level = np.searchsorted(ref.y, values, side="right") / ref.n
    expected = {m: (fn(ref, values, chunk=1) if m == "sMAPE" else fn(ref, values)) for m, fn in MEASURES.items()}
    expected["RMSE"] = np.sqrt(expected["RMSE"] + float(ref.y.var()))
    if mape_infinite:
        expected["MAPE"] = np.full(values.size, np.inf)
    preferred = {m: int(np.argmin(s)) for m, s in expected.items() if np.isfinite(s).any()}
    return Candidates(labels, values, level, expected, preferred, ref.mean, float(np.median(ref.y))), ref


def counts(pmf):
    k = np.arange(0, 2000)
    p = pmf(k)
    keep = np.cumsum(p) < 1 - 1e-16
    keep[: int(np.argmin(keep)) + 1] = True                 # up to the first value where the tail is negligible
    k, p = k[keep].astype(float), p[keep] / p[keep].sum()
    cdf = np.cumsum(p)
    q = k[np.searchsorted(cdf, LEVELS - 1e-12, side="left")]
    mean = float(k @ p)
    pairs = [(f"q{lv:.2f}", float(v)) for lv, v in zip(LEVELS, q)]
    pairs[49] = ("median (q0.50)", pairs[49][1])
    pairs.append(("mean", mean))

    class Exact:                                            # the distribution itself as a weighted "test set"
        vals, counts, n = k, p, 1.0
    smape = lambda f: MEASURES["sMAPE"](Exact, f)
    hi = float(k[np.searchsorted(cdf, 0.999)])
    grid = np.linspace(0.0, max(hi, 1.0), 4001)
    pairs.append(("sMAPE minimiser", _minimise(smape, grid, 1e-6)))
    labels, values = _merge(pairs)
    level = np.array([float(p[k <= v].sum()) for v in values])
    expected = {"RMSE": np.sqrt(((k[:, None] - values[None, :]) ** 2).T @ p),
                "MAE": np.abs(k[:, None] - values[None, :]).T @ p,
                "sMAPE": smape(values)}
    preferred = {m: int(np.argmin(s)) for m, s in expected.items()}
    median = float(k[np.searchsorted(cdf, 0.5 - 1e-12)])
    return Candidates(labels, values, level, expected, preferred, mean, median, zero_share=float(p[0]))
