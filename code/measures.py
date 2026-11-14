"""The four measures, evaluated for a whole grid of candidate point forecasts at once, on the same outcomes.

Each function takes the outcomes `y` (1-d array) and the candidates `f` (1-d array, sorted ascending) and returns one
score per candidate. RMSE, MAE and MAPE use sorted outcomes and cumulative sums (exact, O(n log n + G log n)); sMAPE
has no such form and is evaluated directly, in chunks, or by counts when the outcomes are integers.
"""
import numpy as np


def rmse(y, f):
    m1, m2 = y.mean(), (y * y).mean()
    return np.sqrt(np.maximum(m2 - 2 * f * m1 + f * f, 0.0))


def mae(y, f):
    ys = np.sort(y)
    cs = np.concatenate([[0.0], np.cumsum(ys)])
    n, total = ys.size, cs[-1]
    k = np.searchsorted(ys, f, side="right")            # outcomes <= f
    below = cs[k]
    return (f * k - below + (total - below) - f * (n - k)) / n


def mape(y, f):
    """100 * mean |y - f| / |y| over outcomes other than zero (the zeros are left out: count them with
    `mape_left_out`)."""
    ys = np.sort(y[y != 0])
    z = 1.0 / np.abs(ys)
    cz = np.concatenate([[0.0], np.cumsum(z)])
    n, ztot = ys.size, cz[-1]
    k = np.searchsorted(ys, f, side="right")            # y <= f: |y - f| / y = f / y - 1
    zle = cz[k]
    return 100 * (f * zle - k + (n - k) - f * (ztot - zle)) / n


def mape_left_out(y):
    return int(np.sum(y == 0))


def smape(y, f, chunk=40):
    """200 * mean |y - f| / (|y| + |f|), a pair with y = f (including y = f = 0) contributing 0."""
    vals, counts = np.unique(y, return_counts=True)
    if vals.size < 0.2 * y.size:                         # integer outcomes: evaluate by counts
        w, ys = counts / counts.sum(), vals
        out = np.empty(f.size)
        for i in range(0, f.size, 512):
            ff = f[i:i + 512, None]
            num, den = np.abs(ys - ff), np.abs(ys) + np.abs(ff)
            term = np.where(num == 0, 0.0, num / np.where(den == 0, 1.0, den))
            out[i:i + 512] = 200 * term @ w
        return out
    out = np.empty(f.size)
    for i in range(0, f.size, chunk):
        ff = f[i:i + chunk, None]
        num, den = np.abs(y - ff), np.abs(y) + np.abs(ff)
        term = np.where(num == 0, 0.0, num / np.where(den == 0, 1.0, den))
        out[i:i + chunk] = 200 * term.mean(axis=1)
    return out


MEASURES = {"RMSE": rmse, "MAE": mae, "MAPE": mape, "sMAPE": smape}
