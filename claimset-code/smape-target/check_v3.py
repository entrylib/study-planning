"""Numerical checks of the formulas of the six hypotheses of version 3 (the plan's deciding run).

A proof decides a formal hypothesis; this script checks each formula it proves on concrete numbers. Each result is a
count of violations or a discrepancy that the proof says is zero, or a ratio that it says tends to 1.

    python3 check_v3.py --seed 8 [--rn7 claimset/smape-target/lab/runs/RN7/results.json]

Writes its results to $CLAIMSET_RESULTS when the lab runs it."""
import argparse
import json
import math
import os

import numpy as np
from numpy.polynomial.legendre import leggauss

import dists
import minimiser
import theory
from msmape import EPS, m_term, s_term

LN2 = math.log(2)
X, W = leggauss(200)
HERE = os.path.dirname(os.path.abspath(__file__))
RN7 = os.path.join(HERE, "..", "..", "claimset", "smape-target", "lab", "runs", "RN7", "results.json")


def integ(fn, a, b, pieces):
    e = np.linspace(a, b, pieces + 1)
    tot = 0.0
    for lo, hi in zip(e[:-1], e[1:]):
        u = 0.5 * (hi - lo) * X + 0.5 * (hi + lo)
        tot += 0.5 * (hi - lo) * np.sum(W * fn(u))
    return tot


def norm(u):
    return np.exp(-u * u / 2) / np.sqrt(2 * np.pi)


def norm_cdf(c):
    return 0.5 * (1 + math.erf(c / math.sqrt(2)))


def gumbel(u):
    return np.exp(-(u + np.exp(-u)))


def gumbel_cdf(c):
    return math.exp(-math.exp(-c))


def laplace(u):
    return 0.5 * np.exp(-np.abs(u))


def sech(u):
    return 1 / np.cosh(np.minimum(np.abs(u), 300.0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=8)
    ap.add_argument("--rn7", default=RN7)
    a = ap.parse_args()
    g = np.random.default_rng(a.seed)
    res = []

    def add(h, what, value, **cond):
        res.append({"variable": what, "value": value, "conditions": {"hypothesis": h, **cond}})

    # ---- H3 (c): the log-MAE of the sMAPE minimiser of a test set exceeds the median's by at most mean|z|^3 / 12
    violations, worst = 0, -math.inf
    for rep in range(2000):
        n = int(g.integers(5, 41))
        s = float(g.choice([0.05, 0.2, 0.5, 1.0, 2.0]))
        x = np.sort(s * (g.gamma(2.0, 1.0, n) - 2.0) if rep % 2 else s * g.standard_normal(n))
        grid = np.linspace(x[0], x[-1], 20001)
        t = float(grid[np.tanh(np.abs(x[None, :] - grid[:, None]) / 2).mean(1).argmin()])
        t_exact = float(minimiser.exact(x[None, :])[0])                 # the minimiser is one of the outcomes
        for c in (t, t_exact):
            lhs = np.abs(x - c).mean() - np.abs(x - np.median(x)).mean()
            rhs = (np.abs(x - c) ** 3).mean() / 12
            worst = max(worst, lhs - rhs)
            violations += bool(lhs > rhs + 1e-12)
    add("H3", "violations of the bound on the log-MAE of the sMAPE minimiser", int(violations), test_sets=2000)
    add("H3", "largest excess of the minimiser's log-MAE over the median's, minus the bound", float(worst), test_sets=2000)
    x5 = np.array([0.0, 0.0, 1.0, 4.0, 4.0])
    add("H3", "sum of sMAPE terms, forecast at the median, actuals 1, 1, e, e^4, e^4", float(np.tanh(np.abs(x5 - 1) / 2).sum()))
    add("H3", "sum of sMAPE terms, forecast at 1, actuals 1, 1, e, e^4, e^4", float(np.tanh(np.abs(x5) / 2).sum()))

    # ---- H4: at the target the weights below and above balance; the ordinary median: on which side the target lies
    six = [dists.lognormal(1.0), dists.log_skewed(4.0, 1.0, 1), dists.log_skewed(4.0, 1.0, -1),
           dists.inverse_gaussian(20.0), dists.log_gumbel(1.0, 1), dists.log_gumbel(1.0, -1)]
    bal, wrong_side, many = 0.0, 0, 0
    for d in six:
        T = theory.theory(d)
        many += len(T.minima_plus) != 1
        c = T.minima_plus[0][0]
        gr = T.grid
        total = float((0.5 / np.cosh((gr.x - c) / 2) ** 2) @ gr.gw)
        bal = max(bal, abs(T.balance_exact(c)) / total)
        D, _ = T.one_step()                                             # weight above minus weight below, at the median
        if abs(D) > 1e-9:
            wrong_side += (D > 0) != (c > T.median_x)
    add("H4", "largest difference of the weights below and above the target, over their sum", float(bal), distributions=6)
    add("H34", "distributions with more than one local minimum", int(many), distributions=6)
    add("H34", "distributions where the target is not on the side of the median that carries more weight", int(wrong_side),
        distributions=6)
    sym = dists.pl2()
    dev = max(abs(theory.theory(d).minima_plus[0][0] - theory.theory(d).median_x) for d in sym)
    add("H34", "largest distance in logs between the target and the median, lognormal and log-Laplace", float(dev),
        distributions=len(sym))

    # ---- H35: the Suilin term
    n = 200000
    y = g.choice([-1, 1], n) * 10.0 ** g.uniform(-4, 4, n) * g.random(n)
    f = g.choice([-1, 1], n) * 10.0 ** g.uniform(-4, 4, n) * g.random(n)
    y[: n // 10] = 0.0
    f[n // 10: n // 5] = 0.0
    near = slice(n // 5, n // 5 + n // 10)
    y[near] = g.uniform(0, 0.5, n // 10)
    f[near] = 0.5 - y[near] + g.normal(0, 1e-3, n // 10)
    s, m = s_term(y, f), m_term(y, f)
    tot = np.abs(y) + np.abs(f)
    far, close = tot >= 0.5, tot < 0.5
    add("H35", "largest error of the identity of (a)", float(np.max(np.abs(m[far] - s[far] * tot[far] / (tot[far] + EPS)))), pairs=n)
    add("H35", "largest error of the identity of (b)", float(np.max(np.abs(m[close] - np.abs(y - f)[close] / (0.5 + EPS)))), pairs=n)
    add("H35", "values of the Suilin term outside [0, 1)", int(np.sum((m < 0) | (m >= 1))), pairs=n)
    add("H35", "Suilin term of a zero forecast of a zero actual", float(m_term(0.0, 0.0)))

    # ---- H39: signal-to-noise ratios of a comparison of two forecasts
    cases = (("normal", norm, norm_cdf, -12.0, 12.0), ("Gumbel", gumbel, gumbel_cdf, -6.0, 30.0))
    for name, h, cdf, lo, hi in cases:
        for c1, c2 in ((0.0, 0.5), (0.0, 1.0), (-0.5, 1.5)):
            pts = sorted({lo, hi, c1, c2})
            log_ratios = []
            for sg in (2, 10, 50, 200):
                def parts(fn):
                    return sum(integ(fn, p, q, max(40, 40 * sg)) for p, q in zip(pts[:-1], pts[1:]))

                def snr(dfn):
                    mean = parts(lambda e: dfn(e) * h(e))
                    var = parts(lambda e: dfn(e) ** 2 * h(e)) - mean * mean
                    return abs(mean) / math.sqrt(var)
                cond = {"error_distribution": name, "c1": c1, "c2": c2, "sigma": sg}
                if sg <= 50:
                    log_ratios.append(snr(lambda e: sg * (np.abs(e - c2) - np.abs(e - c1))))
                if sg == 200:
                    h1, h2 = float(h(np.array(c1))), float(h(np.array(c2)))
                    formula = 2 * LN2 / math.sqrt(2 * LN2 - 1) * abs(h1 - h2) / math.sqrt(h1 + h2) / math.sqrt(sg)
                    r = snr(lambda e: np.tanh(sg * np.abs(e - c2) / 2) - np.tanh(sg * np.abs(e - c1) / 2))
                    add("H39", "signal-to-noise ratio of sMAPE over its large-spread formula", float(r / formula), **cond)
                    rho = math.exp(-sg * (c2 - c1))

                    def d_abs(e):                                       # (|y - f2| - |y - f1|) / f2
                        mid = 1 + rho - 2 * np.exp(np.minimum(sg * (e - c2), 0.0))
                        return np.where(e <= c1, 1 - rho, np.where(e >= c2, -(1 - rho), mid))
                    H2 = cdf(c2)
                    limit = abs(2 * H2 - 1) / (2 * math.sqrt(H2 * (1 - H2)))
                    add("H39", "signal-to-noise ratio of the absolute error over its limit", float(snr(d_abs) / limit), **cond)
            spread = (max(log_ratios) - min(log_ratios)) / max(log_ratios)
            add("H39", "relative spread of the absolute log error's signal-to-noise ratio over sigma = 2, 10, 50",
                float(spread), error_distribution=name, c1=c1, c2=c2)
    mu = math.log(100)
    for sg in (0.5, 1.0, 2.0):
        for c1, c2 in ((0.0, 0.5), (0.0, 1.0), (-0.5, 1.5)):
            f1, f2 = math.exp(mu + sg * c1), math.exp(mu + sg * c2)

            def d_sq(e):
                yy = np.exp(mu + sg * e)
                return (yy - f2) ** 2 - (yy - f1) ** 2
            pts = sorted({-14.0, 24.0, c1, c2})
            mean = sum(integ(lambda e: d_sq(e) * norm(e), p, q, 400) for p, q in zip(pts[:-1], pts[1:]))
            second = sum(integ(lambda e: d_sq(e) ** 2 * norm(e), p, q, 400) for p, q in zip(pts[:-1], pts[1:]))
            ey = math.exp(mu + sg * sg / 2)
            sdy = ey * math.sqrt(math.exp(sg * sg) - 1)
            formula = abs(2 * ey - f1 - f2) / (2 * sdy)
            add("H39", "signal-to-noise ratio of the squared error over its formula, lognormal",
                float(abs(mean) / math.sqrt(second - mean * mean) / formula), c1=c1, c2=c2, sigma=sg)
    for sg in (4.0, 6.0):
        ey = math.exp(mu + sg * sg / 2)
        f1, f2 = math.exp(mu), math.exp(mu + sg)
        ratio = abs(2 * ey - f1 - f2) / (2 * ey * math.sqrt(math.exp(sg * sg) - 1))
        add("H39", "exp(sigma^2 / 2) times the squared error's ratio, lognormal, formula", float(math.exp(sg * sg / 2) * ratio),
            c1=0.0, c2=1.0, sigma=sg)

    # ---- H40: the large-sample relative efficiency
    def are(h, spread):
        def gx(u):
            return h(u / spread) / spread
        top = 40 * spread + 60
        aa = 2 * integ(lambda u: gx(u) * sech(u / 2) ** 2 * np.tanh(u / 2), 0, top, 4000)
        bb = 2 * integ(lambda u: gx(u) * sech(u / 2) ** 4, 0, top, 4000)
        g0 = float(gx(np.array(0.0)))
        return (2 * g0 - aa) ** 2 / (4 * g0 ** 2 * bb), aa, g0
    recorded = {}
    if os.path.exists(a.rn7):
        with open(a.rn7) as fh:
            recorded = {r["conditions"]["setting"]: r["value"] for r in json.load(fh)
                        if r["variable"] == "asymptotic relative efficiency"}
    worst_rec, not_below, compared = 0.0, 0, 0
    for name, h, label in (("lognormal", norm, "lognormal sigma={:g}"), ("log-Laplace", laplace, "log-Laplace b={:g}")):
        for sp in dists.SPREADS:
            value, aa, g0 = are(h, sp)
            not_below += not (aa < 2 * g0)
            add("H40", "large-sample relative efficiency, formula", float(value), family=name, spread=sp)
            key = label.format(sp)
            if key in recorded:
                compared += 1
                worst_rec = max(worst_rec, abs(value - recorded[key]))
        add("H40", "distance of the large-sample relative efficiency from 1 at spread 0.01", float(abs(are(h, 0.01)[0] - 1)), family=name)
    add("H40", "distributions where A is not below 2 g at the median", int(not_below), distributions=16)
    add("H40", "distributions compared with the asymptotic values recorded by RN7", int(compared))
    add("H40", "largest difference from the asymptotic values recorded by RN7", float(worst_rec), distributions=compared)
    add("H40", "large-sample relative efficiency times sigma^3 over 6 (ln 2)^2 sqrt(2 pi), lognormal",
        float(are(norm, 200.0)[0] * 200.0 ** 3 / (6 * LN2 ** 2 * math.sqrt(2 * math.pi))), spread=200)
    add("H40", "large-sample relative efficiency times b over 3, log-Laplace", float(are(laplace, 200.0)[0] * 200.0 / 3), spread=200)

    for r in res:
        print(r["conditions"].get("hypothesis"), "|", r["variable"], {k: v for k, v in r["conditions"].items() if k != "hypothesis"},
              f"{r['value']:.6g}")
    path = os.environ.get("CLAIMSET_RESULTS")
    if path:
        with open(path, "w") as fh:
            json.dump(res, fh)


if __name__ == "__main__":
    main()
