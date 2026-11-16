"""The forecast that minimises the sMAPE of a test set of positive outcomes.

In c = ln f, the sMAPE of a test set is 200/n * sum tanh(|x_i - c| / 2) with x_i = ln y_i. Each term is concave in c
on each side of x_i, so the sum is concave between consecutive outcomes and its minimum over c is at an outcome: the
minimiser is one of the test set's outcomes (the lowest-scoring one; ties go to the lower). Small test sets are searched
at every outcome; large ones by a grid, a finer grid around the best two grid points, and then every outcome within the
finer grid's best cells."""
import numpy as np

EXACT_UP_TO = 500
GRID = 41


def _score(x, c):
    """x: (B, n) sorted log outcomes; c: (B, m) log forecasts -> (B, m) sums of tanh(|x - c| / 2)."""
    return np.tanh(np.abs(x[:, None, :] - c[:, :, None]) / 2).sum(axis=2)


def exact(x):
    """Search every outcome. x: (B, n) sorted; returns (B,) log minimisers."""
    s = _score(x, x)
    return x[np.arange(x.shape[0]), np.argmin(s, axis=1)]


def _cells(x, c, s, keep):
    """The `keep` best grid points of each row, with the cell between their neighbours: (B, keep, 2) bounds."""
    order = np.argsort(s, axis=1)[:, :keep]
    lo = np.take_along_axis(c, np.maximum(order - 1, 0), axis=1)
    hi = np.take_along_axis(c, np.minimum(order + 1, c.shape[1] - 1), axis=1)
    return np.stack([lo, hi], axis=2)


def search(x, keep=2):
    """Grid, finer grid, then every outcome in the best cells. x: (B, n) sorted; returns (B,) log minimisers."""
    B, n = x.shape
    lo, hi = x[:, n // 100], x[:, n - 1 - n // 100]
    c = lo[:, None] + (hi - lo)[:, None] * np.linspace(0, 1, GRID)[None, :]
    cells = _cells(x, c, _score(x, c), keep)
    best_c, best_s = np.empty(B), np.full(B, np.inf)
    for k in range(keep):
        a, b = cells[:, k, 0], cells[:, k, 1]
        c2 = a[:, None] + (b - a)[:, None] * np.linspace(0, 1, GRID)[None, :]
        cells2 = _cells(x, c2, _score(x, c2), 2)
        for k2 in range(2):
            a2, b2 = cells2[:, k2, 0], cells2[:, k2, 1]
            for i in range(B):                       # every outcome within the cell, scored exactly
                j0, j1 = np.searchsorted(x[i], [a2[i], b2[i]])
                cand = x[i, max(j0 - 1, 0): min(j1 + 1, n)]
                s = np.tanh(np.abs(x[i][None, :] - cand[:, None]) / 2).sum(axis=1)
                j = int(np.argmin(s))
                if s[j] < best_s[i]:
                    best_s[i], best_c[i] = s[j], cand[j]
    return best_c


def log_minimisers(x, batch=None):
    """x: (R, n) log outcomes (any order) -> (R,) log sMAPE minimisers."""
    x = np.sort(x, axis=1)
    R, n = x.shape
    fn = exact if n <= EXACT_UP_TO else search
    batch = batch or max(1, int(2e7 // (n * (n if n <= EXACT_UP_TO else GRID))))
    return np.concatenate([fn(x[i:i + batch]) for i in range(0, R, batch)])
