"""RMSE, MAE, MAPE and sMAPE of many fixed candidate forecasts f over one test set y, at once.

Each function takes y sorted ascending (with its cumulative sums where needed) and returns one score per candidate.
Scores are the measures as defined in the study's terms; RMSE is returned as its square, which ranks candidates the
same way (the square root is increasing)."""
import numpy as np


class TestSet:
    def __init__(self, y):
        self.y = np.sort(np.asarray(y, dtype=float))
        self.n = self.y.size
        self.cs = np.concatenate([[0.0], np.cumsum(self.y)])
        self.mean = float(self.y.mean())
        positive = self.y > 0
        self.inv = np.concatenate([[0.0], np.cumsum(np.where(positive, 1.0 / np.where(positive, self.y, 1.0), 0.0))])
        self.vals, self.counts = np.unique(self.y, return_counts=True)


def mse(t, f):
    """Mean squared error minus the test set's variance: (mean(y) - f)^2. It orders candidates as RMSE does."""
    return (t.mean - f) ** 2


def mae(t, f):
    k = np.searchsorted(t.y, f, side="left")
    return (f * k - t.cs[k] + (t.cs[-1] - t.cs[k]) - f * (t.n - k)) / t.n


def mape(t, f):
    """100 mean |y - f| / y; for positive outcomes only (the caller does not compute it on counts)."""
    k = np.searchsorted(t.y, f, side="left")
    return 100 * (f * t.inv[k] - k + (t.n - k) - f * (t.inv[-1] - t.inv[k])) / t.n


def smape(t, f, chunk=16):
    """200 mean |y - f| / (|y| + |f|), a pair with y = f contributing 0; computed over the distinct values of y."""
    v, w = t.vals, t.counts / t.n
    out = np.empty(f.size)
    for i in range(0, f.size, chunk):
        ff = f[i:i + chunk, None]
        num, den = np.abs(v - ff), np.abs(v) + np.abs(ff)
        out[i:i + chunk] = 200 * (np.where(num == 0, 0.0, num / np.where(den == 0, 1.0, den)) @ w)
    return out


MEASURES = {"RMSE": mse, "MAE": mae, "MAPE": mape, "sMAPE": smape}
