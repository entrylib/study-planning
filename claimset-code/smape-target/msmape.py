"""msMAPE (Suilin 2017; Hewamalage, Bergmeir and Bandara 2021, Eq. 15) beside sMAPE: the terms, their expected values
under a distribution of the actual, and the forecasts that minimise them, computed from the distribution itself.

    s(y, f) = |y - f| / (|y| + |f|)                      the sMAPE term (0 when y = f = 0)
    m(y, f) = |y - f| / max(|y| + |f| + eps, 0.5 + eps)  the msMAPE term, eps = 0.1

Both measures are 200 times the mean of their term. Distributions: positive continuous ones by the density of
x = ln y (as dists.py describes them), zero-inflated ones (a share p0 of zeros), and Poisson counts.
"""
import math

import numpy as np

from dists import Positive, _normal_pdf
from theory import _gauss

EPS = 0.1
PANELS = 400


def s_term(y, f):
    y, f = np.asarray(y, float), np.asarray(f, float)
    den = np.abs(y) + np.abs(f)
    return np.where(den > 0, np.abs(y - f) / np.where(den > 0, den, 1.0), 0.0)


def m_term(y, f, eps=EPS):
    y, f = np.asarray(y, float), np.asarray(f, float)
    return np.abs(y - f) / np.maximum(np.abs(y) + np.abs(f) + eps, 0.5 + eps)


TERMS = {"sMAPE": s_term, "msMAPE": m_term}


def lognormal_at(median, sigma):
    """A lognormal distribution with the given median (dists.lognormal has median 100)."""
    mu = math.log(median)
    return Positive(f"lognormal median={median:g} sigma={sigma:g}", lambda g, n: g.lognormal(mu, sigma, n),
                    lambda x: _normal_pdf(x, mu, sigma), (mu - 9 * sigma, mu + 9 * sigma),
                    {"mean": True, "second": True, "inv": True}, {}, mu, "lognormal",
                    {"median": median, "sigma": sigma})


class Continuous:
    """A positive continuous distribution, with a share p0 of zeros added."""

    def __init__(self, d, p0=0.0, label=None):
        self.d, self.p0 = d, p0
        self.name = label or (d.name + (f" with zeros p0={p0:g}" if p0 else ""))
        lo, hi = d.support
        self.mass = sum(_gauss(d.pdf_x, a, b, PANELS) for a, b in zip(*self._pieces([lo, hi])))

    def _pieces(self, cuts):
        lo, hi = self.d.support
        c = sorted({min(max(x, lo), hi) for x in cuts} | {lo, hi} |
                   {k for k in self.d.kinks if lo < k < hi})
        return c[:-1], c[1:]

    def expected(self, measure, f):
        """Expected term of `measure` for the forecast f >= 0."""
        term, d = TERMS[measure], self.d
        cuts = [math.log(f)] if f > 0 else []
        if f < 0.5:
            cuts.append(math.log(0.5 - f))                     # where |y| + |f| = 0.5 (msMAPE switches)
        fn = lambda x: term(np.exp(x), f) * d.pdf_x(x)
        pos = sum(_gauss(fn, a, b, PANELS) for a, b in zip(*self._pieces(cuts))) / self.mass
        return float(self.p0 * term(0.0, f) + (1 - self.p0) * pos)

    def level(self, f):
        if f <= 0:
            return self.p0
        lo, _ = self.d.support
        x = math.log(f)
        pos = sum(_gauss(self.d.pdf_x, a, b, PANELS) for a, b in zip(*self._pieces([x]))
                  if b <= x) / self.mass if x > lo else 0.0
        return self.p0 + (1 - self.p0) * pos

    def median(self):
        if self.p0 >= 0.5:
            return 0.0
        lo, hi = self.d.support
        a, b = math.exp(lo), math.exp(hi)
        for _ in range(200):
            mid = math.sqrt(a * b)
            a, b = (mid, b) if self.level(mid) < 0.5 else (a, mid)
        return math.sqrt(a * b)


class Poisson:
    def __init__(self, mean):
        self.mean, self.name = mean, f"Poisson mean={mean:g}"
        k = np.arange(0, int(mean + 40 * math.sqrt(mean) + 60))
        logp = -mean + k * math.log(mean) - np.array([math.lgamma(i + 1) for i in k])
        self.k, self.p = k.astype(float), np.exp(logp)
        self.p0 = float(self.p[0])

    def expected(self, measure, f):
        return float(TERMS[measure](self.k, f) @ self.p)

    def level(self, f):
        return float(self.p[self.k <= f + 1e-12].sum())

    def median(self):
        return float(self.k[np.searchsorted(np.cumsum(self.p), 0.5)])


def _golden(fn, a, b, tol=1e-10):
    g = (math.sqrt(5) - 1) / 2
    c, d = b - g * (b - a), a + g * (b - a)
    fc, fd = fn(c), fn(d)
    while b - a > tol * (1 + abs(a) + abs(b)):
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - g * (b - a)
            fc = fn(c)
        else:
            a, c, fc = c, d, fd
            d = a + g * (b - a)
            fd = fn(d)
    return (a + b) / 2


def target(dist, measure, hi):
    """(forecast, expected score) minimising the expected term over forecasts in [0, hi]: a search on a grid of 0 and
    4000 log-spaced points, refined by golden section around each of the three best grid points (the expected term can
    have several local minima, and a minimum at a kink of the term, at an outcome of a count, is kept as found)."""
    grid = np.concatenate([[0.0], np.geomspace(1e-4, hi, 4000)])
    if isinstance(dist, Poisson):
        grid = np.union1d(grid, dist.k[dist.k <= hi])
    vals = np.array([dist.expected(measure, f) for f in grid])
    best = (float(grid[np.argmin(vals)]), float(vals.min()))
    for i in np.argsort(vals)[:3]:
        a, b = grid[max(i - 1, 0)], grid[min(i + 1, grid.size - 1)]
        if b > a:
            f = _golden(lambda t: dist.expected(measure, t), a, b)
            v = dist.expected(measure, f)
            if v < best[1]:
                best = (f, v)
    return best


def zero_threshold(base, measure, lo=0.0, hi=0.99, tol=1e-6):
    """The share of zeros at and above which the forecast 0 minimises the expected term (bisection on p0)."""
    up = math.exp(base.support[1])
    for _ in range(60):
        mid = (lo + hi) / 2
        f, _ = target(Continuous(base, mid), measure, up)
        lo, hi = (lo, mid) if f == 0.0 else (mid, hi)
        if hi - lo < tol:
            break
    return (lo + hi) / 2
