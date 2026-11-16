"""The outcome distributions of the study, each with a sampler and its exact description.

A positive distribution is described by the density of x = ln y on a grid (`logpdf_x` returns the density of x), its
support in x, and which of its moments are finite, decided from closed forms (hypothesis H12): E Y, E Y^2 and E[1/Y].
A zero-inflated distribution is a positive one with a share p0 of outcomes equal to 0. The straddling settings
(outcomes of both signs) are normal and are described by the density of y itself.

Families (y = 100 exp(x) unless said otherwise):
  inverse Gaussian(mean 100, shape lam)            all moments and E[1/Y] finite
  log-skewed: x = s sd (G - k)/sqrt(k), G ~ Gamma(k, 1), s = +1 (right) or -1 (left): ln y has standard deviation sd
              and skewness 2s/sqrt(k). E Y^p finite iff p s sd/sqrt(k) < 1; E[1/Y] finite iff -s sd/sqrt(k) < 1.
  lognormal(sigma)                                  all finite
  bimodal: x ~ N(-1, 0.25^2) with probability p, N(+1, 0.25^2) otherwise          all finite
  log-Laplace(b): x = b L, L standard Laplace       E Y^p finite iff p b < 1; E[1/Y] finite iff b < 1
  log-Gumbel(sigma, s): x = s sigma E, E standard Gumbel (mode 0, median -ln ln 2): s = +1 right-skewed, -1 left.
              E Y^p = 100^p Gamma(1 - p s sigma), finite iff p s sigma < 1; E[1/Y] finite iff s sigma > -1
  normal straddling zero: y ~ N(mu, 1)              E[1/|Y|] infinite (density at 0)
"""
import math

import numpy as np

LN100 = math.log(100.0)
SQ2PI = math.sqrt(2 * math.pi)


class Positive:
    """kind "positive". logpdf_x(x): density of x = ln y; support: (lo, hi) in x where the density is not negligible
    (beyond it the mass is below about 1e-15); finite: {"mean", "second", "inv"}; moments: closed forms where known
    ({"mean": E Y, "second": E Y^2, "inv": E[1/Y]}, a value only when finite); mode_x: the mode of ln y."""
    kind = "positive"
    p0 = 0.0

    def __init__(self, name, sample, pdf_x, support, finite, moments, mode_x, family, params, kinks=()):
        self.name, self._sample, self.pdf_x, self.support = name, sample, pdf_x, support
        self.kinks = tuple(kinks)                     # points of ln y where the density is not smooth
        self.finite, self.moments, self.mode_x = finite, moments, mode_x
        self.family, self.params = family, params

    def sample(self, g, n):
        return self._sample(g, n)


class ZeroInflated(Positive):
    """A positive distribution with a share p0 of zero outcomes."""
    kind = "zero-inflated"

    def __init__(self, base, p0):
        self.__dict__.update(base.__dict__)
        self.base, self.p0 = base, p0
        self.name = f"{base.name} with zeros p0={p0:g}"

    def sample(self, g, n):
        y = self.base.sample(g, n)
        return np.where(g.random(n) < self.p0, 0.0, y)


class Straddling:
    """Normal outcomes N(mu, sd^2) of both signs."""
    kind = "straddling"
    p0 = 0.0

    def __init__(self, mu, sd=1.0):
        self.mu, self.sd = mu, sd
        self.name = f"normal straddling zero mean={mu:g}"
        self.family, self.params = "normal straddling zero", {"mean": mu, "sd": sd}
        self.finite = {"mean": True, "second": True, "inv": False}
        self.moments = {"mean": mu, "second": mu * mu + sd * sd}
        self.support = (mu - 12 * sd, mu + 12 * sd)

    def sample(self, g, n):
        return g.normal(self.mu, self.sd, n)

    def pdf_y(self, y):
        return np.exp(-0.5 * ((y - self.mu) / self.sd) ** 2) / (self.sd * SQ2PI)


# ------------------------------------------------------------------ families

def _normal_pdf(x, m, s):
    return np.exp(-0.5 * ((x - m) / s) ** 2) / (s * SQ2PI)


def _exp(t):
    return math.exp(t) if t < 700 else math.inf       # a finite moment too large for a float is stored as inf


def lognormal(sigma):
    m = {"mean": 100 * _exp(sigma ** 2 / 2), "second": 1e4 * _exp(2 * sigma ** 2), "inv": _exp(sigma ** 2 / 2) / 100}
    return Positive(f"lognormal sigma={sigma:g}", lambda g, n: g.lognormal(LN100, sigma, n),
                    lambda x: _normal_pdf(x, LN100, sigma), (LN100 - 9 * sigma, LN100 + 9 * sigma),
                    {"mean": True, "second": True, "inv": True}, m, LN100, "lognormal", {"sigma": sigma})


def log_laplace(b):
    fin = {"mean": b < 1, "second": 2 * b < 1, "inv": b < 1}
    m = {}
    if fin["mean"]:
        m["mean"] = 100 / (1 - b * b)
        m["inv"] = 1 / (100 * (1 - b * b))
    if fin["second"]:
        m["second"] = 1e4 / (1 - 4 * b * b)
    return Positive(f"log-Laplace b={b:g}", lambda g, n: 100.0 * np.exp(g.laplace(0.0, b, n)),
                    lambda x: np.exp(-np.abs(x - LN100) / b) / (2 * b), (LN100 - 36 * b, LN100 + 36 * b),
                    fin, m, LN100, "log-Laplace", {"b": b}, kinks=(LN100,))


def log_skewed(k, sd, s):
    """x - ln 100 = s sd (G - k)/sqrt(k), G ~ Gamma(k, 1)."""
    a = sd / math.sqrt(k)                       # x - LN100 = s a (G - k)

    def mgf(t):                                 # E exp(t (x - LN100)) = exp(-t s a k) (1 - t s a)^(-k), if t s a < 1
        return math.exp(-t * s * a * k) * (1 - t * s * a) ** (-k) if t * s * a < 1 else math.inf
    fin = {"mean": s * a < 1, "second": 2 * s * a < 1, "inv": -s * a < 1}
    m = {}
    if fin["mean"]:
        m["mean"] = 100 * mgf(1)
    if fin["second"]:
        m["second"] = 1e4 * mgf(2)
    if fin["inv"]:
        m["inv"] = mgf(-1) / 100
    lgk = math.lgamma(k)

    def pdf_x(x):
        gg = k + s * (x - LN100) / a
        out = np.zeros_like(x)
        pos = gg > 0
        out[pos] = np.exp((k - 1) * np.log(gg[pos]) - gg[pos] - lgk) / a
        if k == 1:                              # the exponential density is 1 at G = 0, the end of the support
            out[np.abs(gg) <= 1e-9] = 1.0 / a
        return out
    lo_g, hi_g = 0.0, k + 40 + 12 * math.sqrt(k)           # the Gamma(k) mass beyond hi_g is below 1e-15
    ends = sorted([LN100 + s * a * (lo_g - k), LN100 + s * a * (hi_g - k)])
    mode_g = max(k - 1, 0.0)
    side = "right" if s > 0 else "left"
    return Positive(f"log-skewed {side} sd={sd:g} skew={2 * s / math.sqrt(k):+g}",
                    lambda g, n: 100.0 * np.exp(s * a * (g.gamma(k, 1.0, n) - k)), pdf_x, tuple(ends), fin, m,
                    LN100 + s * a * (mode_g - k), "log-skewed", {"sd": sd, "skewness": 2 * s / math.sqrt(k)})


def bimodal(p):
    def pdf_x(x):
        return p * _normal_pdf(x, LN100 - 1, 0.25) + (1 - p) * _normal_pdf(x, LN100 + 1, 0.25)
    e = math.exp(0.25 ** 2 / 2)
    m = {"mean": 100 * e * (p * math.exp(-1) + (1 - p) * math.e),
         "second": 1e4 * math.exp(2 * 0.25 ** 2) * (p * math.exp(-2) + (1 - p) * math.exp(2)),
         "inv": e * (p * math.e + (1 - p) * math.exp(-1)) / 100}
    what = {0.5: "equal modes", 0.65: "left mode higher", 0.35: "right mode higher"}.get(p, f"p={p:g}")
    mode = LN100 - 1 if p >= 0.5 else LN100 + 1
    return Positive(f"bimodal {what}", lambda g, n: 100.0 * np.exp(np.where(g.random(n) < p, g.normal(-1.0, 0.25, n),
                                                                             g.normal(1.0, 0.25, n))),
                    pdf_x, (LN100 - 4, LN100 + 4), {"mean": True, "second": True, "inv": True}, m, mode,
                    "bimodal", {"p": p})


def inverse_gaussian(lam, mu=100.0):
    def pdf_x(x):
        y = np.exp(x)
        return np.sqrt(lam / (2 * math.pi * y ** 3)) * np.exp(-lam * (y - mu) ** 2 / (2 * mu * mu * y)) * y
    # mode of ln y: d/dx [ln f_Y(e^x) + x] = 0 -> -1/2 - lam (y^2 - mu^2) / (2 mu^2 y) = 0 -> lam y^2 + mu^2 y - lam mu^2 = 0
    y_mode = (-mu * mu + math.sqrt(mu ** 4 + 4 * lam * lam * mu * mu)) / (2 * lam)
    cv = math.sqrt(mu / lam)
    lo = math.log(mu) - 6 - 3 * cv                 # the lower tail is exp(-lam/(2y)): negligible well before this
    hi = math.log(mu) + math.log(1 + 80 * cv * cv + 40 * cv)
    return Positive(f"inverse Gaussian shape={lam:g}", lambda g, n: g.wald(mu, lam, n), pdf_x, (lo, hi),
                    {"mean": True, "second": True, "inv": True},
                    {"mean": mu, "second": mu * mu + mu ** 3 / lam, "inv": 1 / mu + 1 / lam}, math.log(y_mode),
                    "inverse Gaussian", {"shape": lam})


EULER = 0.5772156649015329


def log_gumbel(sigma, s):
    def pdf_x(x):
        e = s * (x - LN100) / sigma
        return np.exp(-(e + np.exp(-e))) / sigma
    fin = {"mean": s * sigma < 1, "second": 2 * s * sigma < 1, "inv": -s * sigma < 1}
    m = {}
    if fin["mean"]:
        m["mean"] = 100 * math.gamma(1 - s * sigma)
    if fin["second"]:
        m["second"] = 1e4 * math.gamma(1 - 2 * s * sigma)
    if fin["inv"]:
        m["inv"] = math.gamma(1 + s * sigma) / 100
    lo_e, hi_e = -math.log(40 * math.log(10)), 40.0   # standard Gumbel: P(E < lo_e) = exp(-e^{-lo_e}) ~ 1e-40; P(E > 40) ~ 4e-18
    ends = sorted([LN100 + s * sigma * lo_e, LN100 + s * sigma * hi_e])
    side = "right" if s > 0 else "left"
    return Positive(f"log-Gumbel {side} sigma={sigma:g}",
                    lambda g, n: 100.0 * np.exp(s * sigma * g.gumbel(0.0, 1.0, n)), pdf_x, tuple(ends), fin, m,
                    LN100, "log-Gumbel", {"sigma": sigma, "side": side})


# ------------------------------------------------------------------ the settings of each plan

def pl1():
    """The selection settings: those of measure-rankings PL5, and zero-inflated lognormal outcomes (H8)."""
    out = [inverse_gaussian(lam) for lam in (1000.0, 200.0, 100.0, 50.0, 20.0, 10.0)]
    for sd in (0.5, 1.0):
        for s in (-1, 1):
            out += [log_skewed(k, sd, s) for k in (1.0, 4.0, 16.0)]
        out.append(lognormal(sd))
    out += [bimodal(p) for p in (0.5, 0.65, 0.35)]
    out += [Straddling(mu) for mu in (0.0, 0.5, 1.0)]
    out += [ZeroInflated(lognormal(sigma), p0) for sigma, p0s in ZERO_SHARES.items() for p0 in p0s]
    return out


# p0 on both sides of the threshold (200 - m+)/(400 - m+) of H8, fixed in advance from the expected sMAPE of the positive
# part (theory.py: lognormal sigma=0.5, m+ = 38.377, threshold 0.4469; sigma=1, m+ = 69.943, threshold 0.3940)
ZERO_SHARES = {0.5: (0.25, 0.42, 0.47, 0.6), 1.0: (0.25, 0.37, 0.42, 0.6)}

SPREADS = (0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0)


def pl2():
    """Log-symmetric, strictly unimodal settings, where sMAPE and MAE both target the median (H5)."""
    return [lognormal(s) for s in SPREADS] + [log_laplace(b) for b in SPREADS]


FLAT_SPREADS = (0.1, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0)


def pl3():
    """Families whose log-spread grows: symmetric (lognormal) and skewed both ways (log-Gumbel)."""
    out = []
    for s in FLAT_SPREADS:
        out += [lognormal(s), log_gumbel(s, 1), log_gumbel(s, -1)]
    return out
