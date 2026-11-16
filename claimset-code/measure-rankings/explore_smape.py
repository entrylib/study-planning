"""Exploration (not a planned test): where sMAPE's target lies, relative to the median, for normal outcomes of
different levels and standard deviations, and for inverse Gaussian outcomes, whose long right tail makes ln y skewed.

For each distribution, one reference sample of 2,000,000 outcomes; the sMAPE minimiser is found by a grid of the
sample's quantiles at levels 0.005 to 0.995, then golden-section search between the best point's neighbours. The
MAPE minimiser (the weighted median with weights 1 / y) is given beside it. Normal outcomes are drawn as
level x (1 + CV x Z) with the same Z for every level, so that a change of level alone can be seen. Outcomes at or
below zero are counted and left out.

    python3 explore_smape.py --seed 7
"""
import argparse
import json
import math
import os

import numpy as np

from measures import TestSet, smape
from reference import _golden

N = 2_000_000


def skew(x):
    d = x - x.mean()
    return float(np.mean(d ** 3) / np.mean(d ** 2) ** 1.5)


def describe(name, y, extra):
    nonpos = int(np.sum(y <= 0))
    y = np.sort(y[y > 0])
    t = TestSet(y)
    levels = np.linspace(0.005, 0.995, 199)
    q = np.quantile(y, levels)
    s = smape(t, q, chunk=2)
    k = int(np.argmin(s))
    best = _golden(lambda f: float(smape(t, np.array([f]), chunk=1)[0]), q[max(k - 1, 0)], q[min(k + 1, q.size - 1)],
                   (q[-1] - q[0]) * 1e-7)
    w = np.cumsum(1 / y)
    mape_min = float(y[np.searchsorted(w, w[-1] / 2)])
    med = float(np.median(y))
    level = lambda v: float(np.searchsorted(y, v, side="right") / y.size)
    cond = {"distribution": name, **extra}
    out = [("mean of the outcomes", float(y.mean())), ("median of the outcomes", med),
           ("sMAPE minimiser", best), ("quantile level of the sMAPE minimiser", level(best)),
           ("sMAPE minimiser divided by the median", best / med),
           ("MAPE minimiser", mape_min), ("quantile level of the MAPE minimiser", level(mape_min)),
           ("skewness of ln y", skew(np.log(y))), ("outcomes at or below zero", nonpos)]
    return [{"variable": v, "value": x, "conditions": cond} for v, x in out]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    g = np.random.default_rng(a.seed)
    z = g.standard_normal(N)
    results = []
    for level in (100.0, 1000.0):
        for cv in (0.01, 0.05, 0.1, 0.15, 0.2, 0.25):
            results += describe("normal", level * (1 + cv * z), {"level": level, "coefficient of variation": cv})
    for lam in (1000.0, 200.0, 100.0, 50.0, 20.0, 10.0):          # inverse Gaussian, mean 100, shape lam
        y = g.wald(100.0, lam, N)
        results += describe("inverse Gaussian", y, {"mean": 100.0, "shape": lam,
                                                    "coefficient of variation": math.sqrt(100.0 / lam)})
    path = os.environ.get("CLAIMSET_RESULTS")
    if path:
        with open(path, "w") as fh:
            json.dump(results, fh)
    else:
        print(json.dumps(results, indent=1))


if __name__ == "__main__":
    main()
