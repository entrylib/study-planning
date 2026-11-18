"""Exploratory computation made on 18 November 2026, before the hypotheses of version 3 were drafted: do three formulas
derived on paper agree with quadrature? (1) the signal-to-noise ratio of the score difference between two forecasts at
fixed quantile levels, for sMAPE and for the absolute log error, as the spread grows; (2) the large-sample relative
efficiency of the sMAPE minimiser and its rate for wide lognormal and log-Laplace actuals; (3) the bound on how far the
sMAPE minimiser of a test set is from minimising the mean absolute log error. No hypothesis is decided here.
Writes its results to $CLAIMSET_RESULTS when the lab runs it."""
import json
import os

import numpy as np
from numpy.polynomial.legendre import leggauss

LN2 = np.log(2)
X, W = leggauss(200)


def integ(fn, a, b, pieces):
    e = np.linspace(a, b, pieces + 1)
    tot = 0.0
    for lo, hi in zip(e[:-1], e[1:]):
        u = 0.5 * (hi - lo) * X + 0.5 * (hi + lo)
        tot += 0.5 * (hi - lo) * np.sum(W * fn(u))
    return tot


def norm(u):
    return np.exp(-u * u / 2) / np.sqrt(2 * np.pi)


def gumbel(u):
    return np.exp(-(u + np.exp(-u)))


def laplace(u):
    return 0.5 * np.exp(-np.abs(u))


def sech(u):
    return 1 / np.cosh(np.minimum(np.abs(u), 300.0))


def main():
    out = []
    # (1) signal-to-noise ratios
    for name, h, lo, hi in (("normal", norm, -12.0, 12.0), ("Gumbel", gumbel, -6.0, 30.0)):
        for c1, c2 in ((0.0, 0.5), (0.0, 1.0), (-0.5, 1.5)):
            for s in (2, 10, 50):
                pts = sorted({lo, hi, c1, c2})

                def parts(fn):
                    return sum(integ(fn, a, b, max(40, 40 * s)) for a, b in zip(pts[:-1], pts[1:]))
                def d(e):
                    return np.tanh(s * np.abs(e - c2) / 2) - np.tanh(s * np.abs(e - c1) / 2)
                def dl(e):
                    return np.abs(e - c2) - np.abs(e - c1)
                m = parts(lambda e: d(e) * h(e))
                v = parts(lambda e: d(e) ** 2 * h(e)) - m * m
                ml = parts(lambda e: dl(e) * h(e))
                vl = parts(lambda e: dl(e) ** 2 * h(e)) - ml * ml
                h1, h2 = float(h(np.array(c1))), float(h(np.array(c2)))
                formula = 2 * LN2 / np.sqrt(2 * LN2 - 1) * abs(h1 - h2) / np.sqrt(h1 + h2) / np.sqrt(s)
                cond = {"error distribution": name, "c1": c1, "c2": c2, "sigma": s}
                out += [("signal-to-noise ratio, sMAPE", abs(m) / np.sqrt(v), cond),
                        ("signal-to-noise ratio, sMAPE, over its large-spread formula", abs(m) / np.sqrt(v) / formula, cond),
                        ("signal-to-noise ratio, absolute log error", abs(ml) / np.sqrt(vl), cond)]
    # (2) large-sample relative efficiency
    rates = (("lognormal", norm, lambda s: 6 * LN2 ** 2 * np.sqrt(2 * np.pi) / s ** 3), ("log-Laplace", laplace, lambda s: 3 / s))
    for name, h, rate in rates:
        for s in (0.25, 1, 6, 50, 200):
            def g(u):
                return h(u / s) / s
            top = 40 * s + 60
            a = 2 * integ(lambda u: g(u) * sech(u / 2) ** 2 * np.tanh(u / 2), 0, top, 4000)
            b = 2 * integ(lambda u: g(u) * sech(u / 2) ** 4, 0, top, 4000)
            g0 = float(g(np.array(0.0)))
            are = (2 * g0 - a) ** 2 / (4 * g0 ** 2 * b)
            cond = {"family": name, "spread": s}
            out += [("large-sample relative efficiency, formula", are, cond),
                    ("large-sample relative efficiency over its large-spread rate", are / rate(s), cond)]
    # (3) the bound on random test sets
    rng = np.random.default_rng(1)
    worst, violations = -np.inf, 0
    for rep in range(400):
        n = int(rng.integers(5, 40))
        s = float(rng.choice([0.05, 0.2, 0.5, 1.0, 2.0]))
        x = s * (rng.gamma(2.0, 1.0, n) - 2.0) if rep % 2 else s * rng.standard_normal(n)
        grid = np.linspace(x.min(), x.max(), 20001)
        t = grid[np.tanh(np.abs(x[None, :] - grid[:, None]) / 2).mean(1).argmin()]
        lhs = np.abs(x - t).mean() - np.abs(x - np.median(x)).mean()
        rhs = (np.abs(x - t) ** 3).mean() / 12
        worst = max(worst, lhs - rhs)
        violations += lhs > rhs + 1e-9
    out += [("violations of the bound on the log-MAE of the sMAPE minimiser", int(violations), {"test sets": 400}),
            ("largest excess of the log-MAE over the median's, minus the bound", worst, {"test sets": 400})]
    results = [{"variable": v, "value": float(x) if not isinstance(x, int) else x, "conditions": c} for v, x, c in out]
    for r in results:
        print(r["variable"], r["conditions"], f"{r['value']:.6g}")
    path = os.environ.get("CLAIMSET_RESULTS")
    if path:
        with open(path, "w") as fh:
            json.dump(results, fh)


if __name__ == "__main__":
    main()
