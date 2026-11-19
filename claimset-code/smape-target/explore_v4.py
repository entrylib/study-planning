"""Exploratory computation made on 19 November 2026, before the hypotheses of version 4 were drafted, after a reader's
notes on version 3: (1) the large-sample relative efficiency (ARE) of the sMAPE minimiser to the sample median on
log-symmetric shapes, at small and large spreads, and the index h(0)E[Z^2]/E|Z| beside the kurtosis; (2) narrow-core
mixtures, where the ARE exceeds 1; (3) how the asymptotic variance of the sMAPE minimiser grows with the spread for
normal, Laplace and uniform shapes of ln y. No hypothesis is decided here.
Writes its results to $CLAIMSET_RESULTS when the lab runs it."""
import json
import math
import os

import numpy as np
from numpy.polynomial.legendre import leggauss

X, W = leggauss(100)


def integ(fn, a, b, pieces):
    e = np.linspace(a, b, pieces + 1)
    tot = 0.0
    for lo, hi in zip(e[:-1], e[1:]):
        u = 0.5 * (hi - lo) * X + 0.5 * (hi + lo)
        tot += 0.5 * (hi - lo) * np.sum(W * fn(u))
    return tot


def sech(u):
    return 1 / np.cosh(np.minimum(np.abs(u), 300.0))


def normal(z):
    return np.exp(-z * z / 2) / math.sqrt(2 * math.pi)


def t_density(nu):
    c = math.exp(math.lgamma((nu + 1) / 2) - math.lgamma(nu / 2)) / math.sqrt(nu * math.pi)
    return lambda z: c * (1 + z * z / nu) ** (-(nu + 1) / 2)


def gennorm(b):
    return lambda z: b / (2 * math.gamma(1 / b)) * np.exp(-np.abs(z) ** b)


SHAPES = {"uniform": (lambda z: np.where(np.abs(z) < 1, 0.5, 0.0), 1.0),
          "normal": (normal, 40.0), "logistic": (lambda z: 0.25 * sech(z / 2) ** 2, 80.0),
          "t6": (t_density(6), 4000.0), "Laplace": (lambda z: 0.5 * np.exp(-np.abs(z)), 80.0),
          "gennorm 0.75": (gennorm(0.75), 400.0), "gennorm 0.5": (gennorm(0.5), 3000.0)}


def are(h, top, s):
    """(2g - A)^2 / (4 g^2 B) for ln y = ln m + s Z, Z with symmetric density h, by quadrature on z > 0 (split near 0)."""
    edges = [0.0, min(top, 1e-3), min(top, 1.0), top]
    def part(fn):
        return sum(integ(fn, a, b, 400) for a, b in zip(edges[:-1], edges[1:]) if b > a)
    a_ = 2 * part(lambda z: h(z) * sech(s * z / 2) ** 2 * np.tanh(s * z / 2))
    b_ = 2 * part(lambda z: h(z) * sech(s * z / 2) ** 4)
    h0 = float(h(np.array(0.0)))
    return (2 * h0 / s - a_) ** 2 / (4 * (h0 / s) ** 2 * b_), a_, b_


def main():
    out = []
    for name, (h, top) in SHAPES.items():
        h0 = float(h(np.array(0.0)))
        m1 = 2 * integ(lambda z: z * h(z), 0, top, 4000)
        m2 = 2 * integ(lambda z: z * z * h(z), 0, top, 4000)
        out.append(("index h(0)E[Z^2]/E|Z|", h0 * m2 / m1, {"shape": name}))
        for s in (0.05, 0.25, 1.0, 4.0, 16.0):
            out.append(("ARE", are(h, top, s)[0], {"shape": name, "spread": s}))
    for q, s1, s2 in ((0.5, 0.1, 10.0), (0.3, 0.05, 20.0), (0.5, 0.5, 5.0)):
        def h(z, q=q, s1=s1, s2=s2):
            return q * normal(z / s1) / s1 + (1 - q) * normal(z / s2) / s2
        out.append(("ARE", are(h, 12 * s2, 1.0)[0], {"shape": f"normal mixture {q} sd {s1} and sd {s2}", "spread": 1.0}))
    for name in ("normal", "Laplace", "uniform"):
        h, top = SHAPES[name]
        h0 = float(h(np.array(0.0)))
        for s in (2.0, 4.0, 8.0, 16.0):
            r, a_, b_ = are(h, top, s)
            v = b_ / (2 * h0 / s - a_) ** 2
            out.append(("asymptotic variance of the log sMAPE minimiser", v, {"shape": name, "spread": s}))
    res = [{"variable": v, "value": float(x), "conditions": c} for v, x, c in out]
    for r in res:
        print(r["variable"], r["conditions"], f"{r['value']:.6g}")
    path = os.environ.get("CLAIMSET_RESULTS")
    if path:
        with open(path, "w") as fh:
            json.dump(res, fh)


if __name__ == "__main__":
    main()
