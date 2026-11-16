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
        nonzero = self.y != 0
        self.has_zero = not bool(np.all(nonzero))
        absy = np.abs(np.where(nonzero, self.y, 1.0))
        self.inv = np.concatenate([[0.0], np.cumsum(np.where(nonzero, 1.0 / absy, 0.0))])          # sums of 1 / |y|
        self.sgn = np.concatenate([[0.0], np.cumsum(np.where(nonzero, np.sign(self.y), 0.0))])     # sums of y / |y|
        self.vals, self.counts = np.unique(self.y, return_counts=True)


def mse(t, f):
    """Mean squared error minus the test set's variance: (mean(y) - f)^2. It orders candidates as RMSE does."""
    return (t.mean - f) ** 2


def mae(t, f):
    k = np.searchsorted(t.y, f, side="left")
    return (f * k - t.cs[k] + (t.cs[-1] - t.cs[k]) - f * (t.n - k)) / t.n


def mape(t, f):
    """100 mean |y - f| / |y|. An outcome equal to zero makes it infinite for every forecast (returned as inf: MAPE
    fails, and the failure is shown, not hidden). Outcomes of either sign are allowed."""
    if t.has_zero:
        return np.full(np.shape(f), np.inf)
    k = np.searchsorted(t.y, f, side="left")
    # below f: (f - y) / |y| = f / |y| - sign(y); at or above f: (y - f) / |y| = sign(y) - f / |y|
    return 100 * (f * t.inv[k] - t.sgn[k] + (t.sgn[-1] - t.sgn[k]) - f * (t.inv[-1] - t.inv[k])) / t.n


def smape(t, f, chunk=16):
    """200 mean |y - f| / (|y| + |f|), a pair with y = f contributing 0; computed over the distinct values of y."""
    v, w = t.vals, t.counts / t.n
    out = np.empty(f.size)
    for i in range(0, f.size, chunk):
        ff = f[i:i + chunk, None]
        num, den = np.abs(v - ff), np.abs(v) + np.abs(ff)
        out[i:i + chunk] = 200 * (np.where(num == 0, 0.0, num / np.where(den == 0, 1.0, den)) @ w)
    return out


def log_mae(t, f):
    """100 mean |ln y - ln f| (positive outcomes and forecasts): the small-error limit of sMAPE (H9), unbounded."""
    lf = np.log(np.asarray(f, dtype=float))
    return 100 * np.abs(np.log(t.y)[None, :] - lf[:, None]).mean(axis=1)


MEASURES = {"RMSE": mse, "MAE": mae, "MAPE": mape, "sMAPE": smape}
