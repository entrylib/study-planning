"""Numerical experiment of version 4 (plans of draft v4): checks the formulas of H47 to H50 by quadrature and on
simulated test sets.

A. Asymptotic variances (H49, H50 (a)-(c)). For each distribution of the actual and each of RMSE, MAE, MAPE and sMAPE:
   V from the formulas of H49, computed by quadrature of the density of ln y; and V estimated by simulation, n times the
   variance of the log minimiser over R test sets of n actuals. Where sMAPE and MAE both target the median, their ratio
   is the asymptotic relative efficiency (ARE), also estimated by simulation.
B. Growth with the spread (H50 (a)-(e)): the general formulas against the closed forms for lognormal and log-uniform
   actuals; V of sMAPE times sigma^-5 (normal shape), sigma^-3 (Laplace shape) and divided by sigma e^(2 sigma) (uniform
   shape) at growing spreads, against their limits.
C. Efficiency (H48): (a) (ARE - 1)/(sigma^2/2) at small spreads against E[Z^2] - E|Z|/h(0); (b) sigma ARE at large
   spreads against the bound 3L^2/(2h(0)^3) for Lipschitz shapes; (c) ARE of narrow-core mixtures against the limit
   1/(q + (1-q) E_k[sech^4(x/2)]).
D. Two separated groups (H47): expected sMAPE on a fine grid of ln f; its concavity across the gap, where its minimum
   is, the weight balance at its maximum in the gap; and, on simulated test sets, how often the minimiser lies in the
   target's group. A mixture that is not separated checks the local-maximum condition of H47 (e).

Usage: python3 experiment_v4.py [--seed 4] [--reps 2000] [--n 1000] [--quick]
Writes its results to $CLAIMSET_RESULTS when the lab runs it; prints each check and exits 1 if any fails."""
import argparse
import json
import math
import os
import sys

import numpy as np
from numpy.polynomial.legendre import leggauss

import minimiser

GX, GW = leggauss(20)


def sech(u):
    return 1 / np.cosh(np.minimum(np.abs(u), 300.0))


class Dist:
    """A positive actual through the density g of x = ln y on [lo, hi], with breakpoints where g is not smooth.
    Gauss-Legendre panels; an integrand with a kink at c (|x - c|, a sign) is integrated with c as a breakpoint."""

    def __init__(self, name, g, lo, hi, sample, breaks=(), panels=3000):
        self.name, self.g, self.sample = name, g, sample
        self.lo, self.hi, self.breaks, self.panels = lo, hi, tuple(breaks), panels
        self.x, self.w = self.nodes()

    def nodes(self, *extra):
        edges = sorted({self.lo, self.hi, *[b for b in (*self.breaks, *extra) if self.lo < b < self.hi]})
        xs, ws = [], []
        for a, b in zip(edges[:-1], edges[1:]):
            e = np.linspace(a, b, max(20, int(self.panels * (b - a) / (self.hi - self.lo))) + 1)
            mid, half = (e[:-1] + e[1:]) / 2, (e[1:] - e[:-1]) / 2
            xs.append((mid[:, None] + half[:, None] * GX[None, :]).ravel())
            ws.append((half[:, None] * GW[None, :]).ravel())
        x = np.concatenate(xs)
        return x, np.concatenate(ws) * self.g(x)

    def E(self, values):
        return float(np.sum(self.w * values))

    def E_at(self, c, fn):
        """E[fn(x - c)], with c a breakpoint."""
        x, w = self.nodes(c)
        return float(np.sum(w * fn(x - c)))

    def quantile_where(self, weight, level=0.5):
        """The point c with a share `level` of the weighted mass of weight(x) g(x) below it (by bisection)."""
        tot = self.E(weight(self.x))
        a, b = float(self.lo), float(self.hi)
        for _ in range(100):
            c = (a + b) / 2
            x, w = self.nodes(c)
            if float(np.sum(w * weight(x) * (x < c))) < level * tot:
                a = c
            else:
                b = c
        return (a + b) / 2

    def smape_target(self):
        """ln of the sMAPE target: on a grid, the minimiser of E tanh(|x - c|/2); then the root of its derivative
        E[sech^2((x - c)/2) sign(c - x)]/2 next to it, by bisection."""
        cs = np.linspace(self.lo, self.hi, 801)
        s = [self.E(np.tanh(np.abs(self.x - c) / 2)) for c in cs]
        k = int(np.argmin(s))
        a, b = cs[max(k - 1, 0)], cs[min(k + 1, len(cs) - 1)]
        def slope(c):
            return self.E_at(c, lambda z: sech(z / 2) ** 2 * -np.sign(z))
        if slope(a) >= 0 or slope(b) <= 0:
            return float(cs[k])
        for _ in range(100):
            c = (a + b) / 2
            if slope(c) < 0:
                a = c
            else:
                b = c
        return (a + b) / 2

    def variances(self, measures=("RMSE", "MAE", "MAPE", "sMAPE")):
        """V of each measure by the formulas of H49 (inf where the condition fails); ln targets."""
        x, out = self.x, {}
        if "RMSE" not in measures:
            pass
        elif self.finite_square:
            ey, ey2 = self.E(np.exp(x)), self.E(np.exp(2 * x))
            out["RMSE"] = (ey2 / ey ** 2 - 1, math.log(ey))
        else:
            out["RMSE"] = (math.inf, math.nan)
        med = 0.0 if self.symmetric else self.quantile_where(lambda t: np.ones_like(t))
        out["MAE"] = (1 / (4 * float(self.g(np.array(med))) ** 2), med)
        if "MAPE" not in measures:
            pass
        elif self.finite_inverse_square:
            c = self.quantile_where(lambda t: np.exp(-t))
            out["MAPE"] = (self.E_at(c, lambda z: np.exp(-2 * z)) / (4 * float(self.g(np.array(c))) ** 2), c)
        else:
            out["MAPE"] = (math.inf, math.nan)
        c = 0.0 if self.symmetric else self.smape_target()
        d = 2 * float(self.g(np.array(c))) - self.E_at(c, lambda z: sech(z / 2) ** 2 * np.tanh(np.abs(z) / 2))
        out["sMAPE"] = (self.E_at(c, lambda z: sech(z / 2) ** 4) / d ** 2 if d > 0 else math.inf, c)
        return out

    symmetric = False                                  # g symmetric about 0 and unimodal: both targets are 0
    finite_square = True                               # E[y^2] < inf
    finite_inverse_square = True                       # E[y^-2] < inf


def lognormal(s):
    d = Dist(f"lognormal sigma {s}", lambda x: np.exp(-x * x / (2 * s * s)) / (s * math.sqrt(2 * math.pi)),
             -12 * s, 12 * s, lambda rng, size: rng.normal(0, s, size), breaks=(0.0,))
    d.symmetric = True
    return d


def log_laplace(b):
    d = Dist(f"log-Laplace scale {b}", lambda x: np.exp(-np.abs(x) / b) / (2 * b), -40 * b, 40 * b,
             lambda rng, size: rng.laplace(0, b, size), breaks=(0.0,))
    d.finite_square = d.finite_inverse_square = b < 0.5      # E[y^2] = E[y^-2] = 1/(1 - 4b^2) when b < 1/2
    d.symmetric = True
    return d


def log_uniform(a):
    d = Dist(f"log-uniform half-width {a}", lambda x: np.where(np.abs(x) < a, 0.5 / a, 0.0), -a, a,
             lambda rng, size: rng.uniform(-a, a, size), breaks=(0.0,))
    d.symmetric = True
    return d


def gamma(k):
    lg = math.lgamma(k)
    d = Dist(f"gamma shape {k}", lambda x: np.exp(k * x - np.exp(np.minimum(x, 700)) - lg),
             math.log(1e-300) / k if k < 1 else -60 / k - 5, math.log(k + 40 * math.sqrt(k) + 60),
             lambda rng, size: np.log(rng.gamma(k, 1.0, size)))
    d.finite_inverse_square = k > 2
    return d


def inverse_gaussian(mu, lam):
    def g(x):
        y = np.exp(x)
        return np.sqrt(lam / (2 * math.pi * y)) * np.exp(-lam * (y - mu) ** 2 / (2 * mu * mu * y))
    return Dist(f"inverse Gaussian mean {mu} shape {lam}", g, math.log(mu) - 40, math.log(mu) + 8,
                lambda rng, size: np.log(rng.wald(mu, lam, size)))


def weibull(k):
    d = Dist(f"Weibull shape {k}", lambda x: k * np.exp(k * x - np.exp(np.minimum(k * x, 700))), -60 / k, 5 / k,
             lambda rng, size: np.log(rng.weibull(k, size)))
    d.finite_inverse_square = False                    # E[y^-2] = Gamma(1 - 2/k): infinite for k <= 2
    return d


def gennorm(beta, s=1.0):
    c = beta / (2 * math.gamma(1 / beta) * s)
    top = s * (60.0 ** (1 / beta))
    def sample(rng, size):
        return s * rng.choice([-1.0, 1.0], size) * rng.gamma(1 / beta, 1.0, size) ** (1 / beta)
    d = Dist(f"log generalised normal power {beta} spread {s}", lambda x: c * np.exp(-np.abs(x / s) ** beta),
             -top, top, sample, breaks=(0.0, -s * 1e-3, s * 1e-3), panels=20000)
    d.finite_square = d.finite_inverse_square = beta >= 1   # tails heavier than exponential: no E[y^2], E[y^-2]
    d.symmetric = True
    return d


def mixture(q, eps, k_sd):
    """ln y - ln m: normal sd eps with probability q, else normal sd k_sd."""
    def g(x):
        return (q * np.exp(-x * x / (2 * eps * eps)) / (eps * math.sqrt(2 * math.pi))
                + (1 - q) * np.exp(-x * x / (2 * k_sd * k_sd)) / (k_sd * math.sqrt(2 * math.pi)))
    def sample(rng, size):
        return np.where(rng.random(size) < q, rng.normal(0, eps, size), rng.normal(0, k_sd, size))
    lo = 12 * k_sd
    d = Dist(f"normal mixture q {q} sd {eps} and {k_sd}", g, -lo, lo, sample,
             breaks=(0.0, -12 * eps, 12 * eps, -eps, eps), panels=40000)
    d.symmetric = True
    return d


# --- sample minimisers, one per row of x = ln y (R, n) -------------------------------------------------------------

def sample_logs(x, measures=("RMSE", "MAE", "MAPE", "sMAPE")):
    out = {}
    if "RMSE" in measures:
        out["RMSE"] = np.log(np.exp(x).mean(axis=1))
    if "MAE" in measures:
        out["MAE"] = np.median(x, axis=1)
    if "MAPE" in measures:
        xs = np.sort(x, axis=1)
        cw = np.cumsum(np.exp(-(xs - xs[:, :1])), axis=1)   # weights 1/y, scaled by the row's smallest y
        j = np.argmax(cw >= cw[:, -1:] / 2, axis=1)        # weighted median
        out["MAPE"] = xs[np.arange(x.shape[0]), j]
    if "sMAPE" in measures:
        out["sMAPE"] = minimiser.log_minimisers(x)
    return out


def var_and_se(v):
    """Sample variance of v and its standard error, from the fourth central moment."""
    d = v - v.mean()
    s2 = float(np.mean(d * d)) * len(v) / (len(v) - 1)
    m4 = float(np.mean(d ** 4))
    return s2, math.sqrt(max(m4 - s2 * s2, 0.0) / len(v))


SIZES = (6, 12, 24, 48, 100, 250, 1000)


def trigamma(x):
    """psi_1(x) for x > 0: recurrence up to 6, then the asymptotic series."""
    t = 0.0
    while x < 6:
        t += 1 / (x * x)
        x += 1
    return t + 1 / x + 1 / (2 * x ** 2) + 1 / (6 * x ** 3) - 1 / (30 * x ** 5) + 1 / (42 * x ** 7) - 1 / (30 * x ** 9)


def exact_median_variance(d, n, points=400001):
    """Exact variance of the sample median of x = ln y for a test set of n actuals (the mean of the two middle order
    statistics when n is even, the middle one when odd), from the joint density of order statistics, on a fine grid:
    a uniform grid over the support, joined with a second one around the median, 40 times as wide as the median's
    standard deviation for n actuals (so that a sharp peak at the median is resolved)."""
    x0 = np.linspace(d.lo, d.hi, points)
    w0 = d.g(x0)
    F0 = np.cumsum(w0) / w0.sum()
    med0 = float(np.interp(0.5, F0, x0))
    half = 20 / (2 * max(float(d.g(np.array(med0))), 1e-12) * math.sqrt(n))
    x = np.unique(np.concatenate([x0, np.linspace(max(d.lo, med0 - half), min(d.hi, med0 + half), points)]))
    dx = np.gradient(x)
    g = d.g(x) * dx
    g = g / g.sum()
    F = np.cumsum(g) - g / 2
    med = float(np.interp(0.5, F, x))
    u = x - med
    if n % 2:
        k = (n + 1) // 2
        lw = (k - 1) * np.log(np.maximum(F, 1e-300)) + (n - k) * np.log(np.maximum(1 - F, 1e-300))
        w = np.exp(lw - lw.max()) * g                      # one factor: its largest value can be divided out
        w = w / w.sum()
        m1 = float(np.sum(w * u))
        return float(np.sum(w * u * u)) - m1 * m1
    k = n // 2                                   # X_(k) = a < b = X_(k+1)
    # F(a)^(k-1) (1 - F(b))^(n-k-1), scaled by 2^(n-2) so that it is about 1 at the median; the constant factor cancels
    # when divided by the mass. Where an exponent is clipped, the pair's product is negligible.
    la = np.clip((k - 1) * np.log(np.maximum(2 * F, 1e-300)), -745, 700)
    lb = np.clip((n - k - 1) * np.log(np.maximum(2 * (1 - F), 1e-300)), -745, 700)
    A = np.exp(la) * g
    B = np.exp(lb) * g
    def tail(v):                                 # sum over j > i of v_j
        c = np.cumsum(v[::-1])[::-1]
        return np.concatenate([c[1:], [0.0]])
    s0, s1, s2 = tail(B), tail(B * u), tail(B * u * u)
    mass = float(np.sum(A * s0))
    ea, eb = float(np.sum(A * u * s0)) / mass, float(np.sum(A * s1)) / mass
    ea2, eb2, eab = float(np.sum(A * u * u * s0)) / mass, float(np.sum(A * s2)) / mass, float(np.sum(A * u * s1)) / mass
    mean = (ea + eb) / 2
    return (ea2 + 2 * eab + eb2) / 4 - mean * mean


def exact_rmse_variance(name, n):
    """Exact variance of ln of the sample mean: gamma (shape k): the sum is gamma with shape nk, so trigamma(nk);
    inverse Gaussian (mean mu, shape lam): the sample mean is inverse Gaussian with shape n lam, by quadrature."""
    if name.startswith("gamma shape "):
        return trigamma(n * float(name.split()[-1]))
    if name.startswith("inverse Gaussian "):
        w = name.split()
        dm = inverse_gaussian(float(w[3]), n * float(w[5]))
        mass = dm.E(np.ones_like(dm.x))
        m1 = dm.E(dm.x) / mass
        return dm.E((dm.x - m1) ** 2) / mass
    return None


def extremity(d):
    """E|x - E x|^(1/2) / sd(x)^(1/2) for x = ln y (the supervisor's measure, TM46)."""
    m = d.E(d.x) / d.E(np.ones_like(d.x))
    sd = math.sqrt(d.E((d.x - m) ** 2) / d.E(np.ones_like(d.x)))
    return d.E_at(m, lambda z: np.sqrt(np.abs(z))) / d.E(np.ones_like(d.x)) / math.sqrt(sd)


def part_a(rng, reps, n, quick):
    settings = [lognormal(s) for s in (0.5, 1.0, 2.0)] + [log_laplace(b) for b in (0.25, 0.5, 1.0)] \
        + [log_uniform(a) for a in (0.5, 1.0, 2.0)] + [gamma(k) for k in (20, 5, 1, 0.5)] \
        + [inverse_gaussian(1.0, lam) for lam in (1.0, 0.2, 0.05)] + [weibull(k) for k in (2.0, 1.0, 0.5)] \
        + [gennorm(0.5, 1.0), mixture(0.5, 0.05, 1.0)]
    sizes = tuple(m for m in SIZES if m < n) + (n,)
    if quick:
        settings, sizes = settings[::4], (6, 24, n)
    results, checks = [], []
    for d in settings:
        vf = d.variances()
        if d.symmetric:
            c = d.smape_target()
            checks.append((f"A {d.name}: sMAPE target found numerically at ln y = {c:.2e}, the median", abs(c) < 1e-6))
            results.append(("extremity", extremity(d), {"distribution": d.name}))
        both_median = abs(vf["sMAPE"][1] - vf["MAE"][1]) < 1e-6
        if both_median:
            are_f = vf["MAE"][0] / vf["sMAPE"][0]
            results.append(("asymptotic relative efficiency by formula", are_f, {"distribution": d.name}))
        for size in sizes:
            logs = sample_logs(d.sample(rng, (reps, size)))
            for m in ("RMSE", "MAE", "MAPE", "sMAPE"):
                v, _ = vf[m]
                s2, se = var_and_se(logs[m])
                cond = {"distribution": d.name, "measure": m, "n": size, "test sets": reps}
                results.append(("asymptotic variance by simulation", size * s2, cond))
                results.append(("standard error of the simulated variance", size * se, cond))
                if size == n:
                    results.append(("asymptotic variance by formula", v, {k: cond[k] for k in ("distribution", "measure")}))
                if math.isfinite(v):
                    results.append(("n times variance over formula", size * s2 / v, cond))
                exact = None
                if m == "MAE":
                    exact = exact_median_variance(d, size)
                elif m == "RMSE":
                    exact = exact_rmse_variance(d.name, size)
                if exact is not None:
                    results.append(("exact variance of the log minimiser times n", size * exact, cond))
                    checks.append((f"E {d.name} {m} n {size}: simulated n Var {size * s2:.4g} vs exact {size * exact:.4g}",
                                   abs(s2 - exact) <= 4 * se + 0.01 * exact))
                if size == n and math.isfinite(v):
                    ok = abs(n * s2 - v) <= 4 * n * se + 0.05 * v
                    text = f"A {d.name} {m}: simulated V {n * s2:.4g} vs formula {v:.4g}"
                    if not ok:                   # the plan's second chance: 4n actuals per test set, same tolerance
                        n4 = 4 * n
                        s2b, seb = var_and_se(sample_logs(d.sample(rng, (reps, n4)), (m,))[m])
                        results.append(("asymptotic variance by simulation", n4 * s2b, {**cond, "n": n4}))
                        results.append(("standard error of the simulated variance", n4 * seb, {**cond, "n": n4}))
                        ok = abs(n4 * s2b - v) <= 4 * n4 * seb + 0.05 * v
                        text += f"; with {n4} actuals {n4 * s2b:.4g}"
                    checks.append((text, ok))
            if both_median:
                are_s = var_and_se(logs["MAE"])[0] / var_and_se(logs["sMAPE"])[0]
                results.append(("relative efficiency by simulation", are_s,
                                {"distribution": d.name, "n": size, "test sets": reps}))
                if size == n:
                    checks.append((f"A {d.name}: simulated ARE {are_s:.4g} on the same side of 1 as {are_f:.4g}",
                                   abs(are_f - 1) < 0.05 or (are_s - 1) * (are_f - 1) > 0))
        print(d.name, {m: f"{vf[m][0]:.4g}" for m in vf}, flush=True)
    return results, checks


def part_b():
    results, checks = [], []
    for s in (0.5, 1.0, 2.0):
        d, vf = lognormal(s), None
        vf = d.variances()
        closed = {"RMSE": math.exp(s * s) - 1, "MAE": math.pi * s * s / 2, "MAPE": math.pi * s * s * math.exp(s * s) / 2}
        for m, v in closed.items():
            checks.append((f"B lognormal {s} {m}: general formula {vf[m][0]:.6g} vs closed form {v:.6g}",
                           abs(vf[m][0] / v - 1) < 1e-4))
    for a in (0.5, 1.0, 2.0, 4.0):
        vf = log_uniform(a).variances()
        t = math.tanh(a / 2)
        closed = {"RMSE": a / math.tanh(a) - 1, "MAE": a * a, "MAPE": a * math.tanh(a),
                  "sMAPE": 2 * a * t * (3 - t * t) / (3 * (1 - t * t) ** 2)}
        for m, v in closed.items():
            checks.append((f"B log-uniform {a} {m}: general formula {vf[m][0]:.6g} vs closed form {v:.6g}",
                           abs(vf[m][0] / v - 1) < 1e-4))
    h0 = 1 / math.sqrt(2 * math.pi)
    lim_n = h0 / (24 * math.log(2) ** 2 * h0 ** 2)
    lim_l = 0.5 / (6 * 0.25)
    for s in (10.0, 20.0, 50.0, 100.0, 200.0):
        vn = lognormal(s).variances(("sMAPE",))["sMAPE"][0] / s ** 5
        vl = log_laplace(s).variances(("sMAPE",))["sMAPE"][0] / s ** 3
        results.append(("sMAPE asymptotic variance over sigma^5", vn, {"shape": "normal", "spread": s}))
        results.append(("sMAPE asymptotic variance over sigma^3", vl, {"shape": "Laplace", "spread": s}))
        if s == 200.0:
            checks.append((f"B normal shape: V/sigma^5 {vn:.5g} at 200 vs limit {lim_n:.5g}", abs(vn / lim_n - 1) < 0.03))
            checks.append((f"B Laplace shape: V/sigma^3 {vl:.5g} at 200 vs limit {lim_l:.5g}", abs(vl / lim_l - 1) < 0.03))
    for a in (2.0, 4.0, 8.0, 16.0):
        v = log_uniform(a).variances(("sMAPE",))["sMAPE"][0] / (a * math.exp(2 * a))
        results.append(("sMAPE asymptotic variance over sigma e^(2 sigma)", v, {"shape": "uniform", "spread": a}))
        if a == 16.0:
            checks.append((f"B uniform shape: V/(sigma e^(2 sigma)) {v:.5g} at 16 vs limit 1/12", abs(v * 12 - 1) < 0.03))
    return results, checks


SHAPES = {  # density h of Z, its Lipschitz constant (None: not Lipschitz at 0), sampler-free
    "normal": (lambda z: np.exp(-z * z / 2) / math.sqrt(2 * math.pi), math.exp(-0.5) / math.sqrt(2 * math.pi)),
    "logistic": (lambda z: 0.25 * sech(z / 2) ** 2, 1 / (6 * math.sqrt(3))),
    "Student t 6": (lambda z: math.exp(math.lgamma(3.5) - math.lgamma(3)) / math.sqrt(6 * math.pi)
                    * (1 + z * z / 6) ** -3.5, None),
    "Laplace": (lambda z: 0.5 * np.exp(-np.abs(z)), 0.5),
    "uniform": (lambda z: np.where(np.abs(z) < 1, 0.5, 0.0), None),
    "generalised normal 0.75": (lambda z: 0.75 / (2 * math.gamma(1 / 0.75)) * np.exp(-np.abs(z) ** 0.75), None),
}


def shape_dist(name, s):
    h, _ = SHAPES[name]
    top = {"normal": 40, "logistic": 80, "Student t 6": 4000, "Laplace": 80, "uniform": 1,
           "generalised normal 0.75": 400}[name] * s
    rings = [r for k in range(-3, 8) for r in (10.0 ** k, 3 * 10.0 ** k) if r < top]  # fine panels near 0
    d = Dist(f"{name} spread {s}", lambda x: h(x / s) / s, -top, top, None,
             breaks=(0.0, -s, s, *rings, *[-r for r in rings]), panels=20000)
    d.symmetric = True
    return d


def are_of(d):
    vf = d.variances(("MAE", "sMAPE"))
    return vf["MAE"][0] / vf["sMAPE"][0]


def t6_lipschitz():
    z = np.linspace(0, 10, 200001)
    h = SHAPES["Student t 6"][0](z)
    return float(np.max(np.abs(np.diff(h)) / np.diff(z)))


def part_c():
    results, checks = [], []
    for name, (h, _) in SHAPES.items():
        d1 = shape_dist(name, 1.0)
        h0 = float(h(np.array(0.0)))
        m1, m2 = d1.E(np.abs(d1.x)), d1.E(d1.x ** 2)
        lim = m2 - m1 / h0
        results.append(("peakedness index", h0 * m2 / m1, {"shape": name}))
        results.append(("extremity", extremity(d1), {"shape": name}))
        r = (are_of(shape_dist(name, 0.01)) - 1) / (0.01 ** 2 / 2)
        results.append(("(ARE - 1)/(sigma^2/2)", r, {"shape": name, "spread": 0.01}))
        checks.append((f"C(a) {name}: (ARE-1)/(sigma^2/2) {r:.4g} at 0.01 vs {lim:.4g}",
                       abs(r - lim) <= 0.02 * max(1.0, abs(lim))))
    for name in ("normal", "logistic", "Student t 6", "Laplace"):
        h, L = SHAPES[name]
        L = L if L is not None else t6_lipschitz()
        h0 = float(h(np.array(0.0)))
        bound = 3 * L * L / (2 * h0 ** 3)
        for s in (10.0, 100.0, 1000.0):
            v = s * are_of(shape_dist(name, s))
            results.append(("sigma times ARE", v, {"shape": name, "spread": s}))
            checks.append((f"C(b) {name}: sigma ARE {v:.4g} at {s} within bound {bound:.4g}", v <= bound))
    for q, k in ((0.5, 1.0), (0.3, 2.0), (0.8, 0.5)):
        kd = Dist("k", lambda x: np.exp(-x * x / (2 * k * k)) / (k * math.sqrt(2 * math.pi)), -12 * k, 12 * k, None)
        lim = 1 / (q + (1 - q) * kd.E(sech(kd.x / 2) ** 4))
        for eps in (1e-2, 1e-3):
            dm = mixture(q, eps, k)
            a = are_of(dm)
            results.append(("ARE of a narrow-core mixture", a, {"q": q, "core sd": eps, "wide sd": k, "limit": lim}))
            results.append(("extremity", extremity(dm), {"q": q, "core sd": eps, "wide sd": k}))
            if eps == 1e-3:
                checks.append((f"C(c) q {q} wide sd {k}: ARE {a:.5g} at core sd 0.001 vs limit {lim:.5g}",
                               abs(a / lim - 1) < 0.01))
    return results, checks


def expected_smape_grid(x, w, cs):
    """E tanh(|X - c|/2) for each c: X takes values x with probabilities w."""
    return np.array([np.sum(w * np.tanh(np.abs(x - c) / 2)) for c in cs])


def groups(l1, u1, l2, u2, p, points=4000):
    """Two uniform groups on the log scale (or points when l = u), as atoms on a fine grid."""
    def part(l, u, mass):
        if l == u:
            return np.array([l]), np.array([mass])
        e = np.linspace(l, u, points + 1)
        return (e[:-1] + e[1:]) / 2, np.full(points, mass / points)
    x1, w1 = part(l1, u1, p)
    x2, w2 = part(l2, u2, 1 - p)
    return np.concatenate([x1, x2]), np.concatenate([w1, w2])


def part_d(rng, quick):
    results, checks = [], []
    cases = [(0.0, 0.5, 2.0, 4.0, 0.5), (0.0, 0.5, 2.0, 4.0, 0.6), (0.0, 0.5, 2.0, 4.0, 0.4),
             (-1.0, 0.0, 1.0, 1.2, 0.3), (-1.0, 0.0, 1.0, 1.2, 0.5), (0.0, 0.0, 2.0, 2.0, 0.5),
             (0.0, 0.0, 2.0, 2.0, 0.55), (0.0, 0.0, 2.0, 2.0, 0.45), (-1.0, 1.0, 1.5, 1.7, 0.5)]
    for l1, u1, l2, u2, p in cases:
        x, w = groups(l1, u1, l2, u2, p)
        cond = {"groups": f"[{l1}, {u1}] and [{l2}, {u2}]", "p": p}
        gap = np.linspace(u1, l2, 2001)[1:-1]
        sg = expected_smape_grid(x, w, gap)
        sec = np.diff(sg, 2)
        checks.append((f"D(a) {cond}: concave across the gap (max second difference {sec.max():.3g})", sec.max() < 0))
        checks.append((f"D(a) {cond}: every gap forecast worse than an end",
                       bool(np.all(sg > min(expected_smape_grid(x, w, [u1, l2])) - 1e-15))))
        below = np.linspace(l1 - 2, l1, 401)[:-1]
        above = np.linspace(u2, u2 + 2, 401)[1:]
        se = expected_smape_grid(x, w, [l1, u2])
        checks.append((f"D(b) {cond}: worse below the lowest and above the highest group",
                       bool(np.all(expected_smape_grid(x, w, below) > se[0]) and
                            np.all(expected_smape_grid(x, w, above) > se[1]))))
        allc = np.unique(np.concatenate([np.linspace(l1, u1, 801), np.linspace(l2, u2, 801), gap]))
        sa = expected_smape_grid(x, w, allc)
        tgt = allc[np.argmin(sa)]
        in_group = 1 if l1 <= tgt <= u1 else 2 if l2 <= tgt <= u2 else 0
        results.append(("ln sMAPE target", float(tgt), cond))
        checks.append((f"D(b) {cond}: target {tgt:.4f} lies in group {in_group}", in_group > 0))
        if l1 == u1 and l2 == u2:
            vals = expected_smape_grid(x, w, [l1, l2])
            want = "both" if p == 0.5 else ("low" if p > 0.5 else "high")
            got = "both" if abs(vals[0] - vals[1]) < 1e-12 else ("low" if vals[0] < vals[1] else "high")
            checks.append((f"D(d) {cond}: target {got}, expected {want}", got == want))
        k = int(np.argmax(sg))
        f = math.exp(gap[k])
        y = np.exp(x)
        wt = w * y / (y + f) ** 2
        bal = (wt[x < gap[k]].sum() - wt[x > gap[k]].sum()) / wt.sum()
        results.append(("weight balance at the maximum in the gap", float(bal), cond))
        if 0 < k < len(gap) - 1:
            checks.append((f"D(e) {cond}: weights balance at the gap maximum ({bal:.2e})", abs(bal) < 1e-3))
        if not quick:
            for n in (50, 500):
                xs = np.where(rng.random((400, n)) < p,
                              rng.uniform(l1, u1, (400, n)) if u1 > l1 else l1,
                              rng.uniform(l2, u2, (400, n)) if u2 > l2 else l2)
                mins = minimiser.log_minimisers(xs)
                g_of = np.where(mins <= u1 + 1e-12, 1, 2)
                share = float(np.mean(g_of == in_group))
                results.append(("share of test-set minimisers in the target's group", share, {**cond, "n": n}))
                # an estimate, not a check: H47 says nothing about finite test sets (supervisor, 2026-11-20)
    # H47 (e), not separated: equal normal groups at +-1.5 with sd 0.4 on the log scale, mirror images about 0
    d = Dist("mirror mixture", lambda t: 0.5 * (np.exp(-(t - 1.5) ** 2 / 0.32) + np.exp(-(t + 1.5) ** 2 / 0.32))
             / (0.4 * math.sqrt(2 * math.pi)), -6, 6, None, breaks=(0.0,))
    g0 = float(d.g(np.array(0.0)))
    a = d.E(sech(d.x / 2) ** 2 * np.tanh(np.abs(d.x) / 2))
    s0, s1 = (d.E(np.tanh(np.abs(d.x - c) / 2)) for c in (0.0, 0.05))
    results.append(("2g(ln m) - E[sech^2 tanh] at the median", 2 * g0 - a, {"distribution": "mirror mixture"}))
    checks.append((f"D(e) mirror mixture: 2g {2 * g0:.4g} < {a:.4g} and the median is a local maximum",
                   2 * g0 < a and s0 > s1))
    return results, checks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=4)
    ap.add_argument("--reps", type=int, default=2000)
    ap.add_argument("--n", type=int, default=1000)
    ap.add_argument("--quick", action="store_true", help="fewer settings and test sets, to test the code")
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    out, checks = [], []
    for part in (lambda: part_a(rng, args.reps, args.n, args.quick), part_b, part_c, lambda: part_d(rng, args.quick)):
        r, c = part()
        out += r
        checks += c
    for text, ok in checks:
        print("PASS" if ok else "FAIL", text)
    res = [{"variable": v, "value": float(x), "conditions": c} for v, x, c in out if math.isfinite(x)]  # inf: no V
    res.append({"variable": "checks failed", "value": sum(not ok for _, ok in checks),
                "conditions": {"checks": len(checks)}})
    path = os.environ.get("CLAIMSET_RESULTS")
    if path:
        with open(path, "w") as fh:
            json.dump(res, fh)
    print(f"{sum(ok for _, ok in checks)} of {len(checks)} checks pass")
    sys.exit(0 if all(ok for _, ok in checks) else 1)


if __name__ == "__main__":
    main()
