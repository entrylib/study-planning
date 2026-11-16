"""The outcome distributions. Continuous ones are positive (the supervisor's decision DC3): normal and Laplace sit at a
level of 100 with standard deviations 1 to 15 (normal) and 1 to 7.5 (Laplace); the outlier settings keep their outliers
at least six standard deviations above zero. Counts may be zero; MAPE is not computed on them."""
import math

import numpy as np


def _contaminated(g, n):
    """95% N(100, 5^2), 5% N(100, 15^2): symmetric outliers three times as spread."""
    return 100.0 + 5.0 * g.standard_normal(n) * np.where(g.random(n) < 0.05, 3.0, 1.0)


def _spikes(g, n):
    """N(100, 5^2), with 5% of outcomes three times as large: upward outliers, as from promotions."""
    return g.normal(100.0, 5.0, n) * np.where(g.random(n) < 0.05, 3.0, 1.0)


def continuous():
    """(name, sampler(g, n)) for the continuous settings."""
    # level 100; the standard deviation sets the coefficient of variation. Laplace stops at 7.5: its exponential tail
    # puts a share 0.5 exp(-sqrt(2) / CV) of outcomes at or below zero (3e-9 at CV 0.075, 4e-7 at 0.1)
    out = [(f"normal sd={sd:g}", lambda g, n, sd=sd: g.normal(100.0, sd, n)) for sd in (1, 5, 10, 15)]
    out += [(f"Laplace sd={sd:g}", lambda g, n, sd=sd: g.laplace(100.0, sd / math.sqrt(2), n)) for sd in (1, 2.5, 5, 7.5)]
    for s in (0.1, 0.25, 0.5, 1.0):
        out.append((f"lognormal s={s}", lambda g, n, s=s: g.lognormal(math.log(100.0), s, n)))
    for b in (0.25, 0.5):
        out.append((f"log-Laplace b={b}", lambda g, n, b=b: 100.0 * np.exp(g.laplace(0.0, b, n))))
    for k in (2, 5, 20):
        out.append((f"gamma shape={k}", lambda g, n, k=k: g.gamma(k, 100.0 / k, n)))
    out += [("normal with symmetric outliers", _contaminated), ("normal with upward spikes", _spikes)]
    return out


def _pois(lam):
    def pmf(k):
        return np.exp(k * math.log(lam) - lam - np.array([math.lgamma(x + 1) for x in k]))
    return pmf


def _negbin(mean, size=1.0):
    p = size / (size + mean)

    def pmf(k):
        return np.exp(np.array([math.lgamma(x + size) - math.lgamma(size) - math.lgamma(x + 1) for x in k])
                      + size * math.log(p) + k * math.log(1 - p))
    return pmf


def counts():
    """(name, sampler(g, n), pmf(k)) for the count settings: Poisson, and intermittent ones (more than half zeros)."""
    out = []
    for lam in (1.5, 3.5, 10.0):
        out.append((f"Poisson lambda={lam}", lambda g, n, lam=lam: g.poisson(lam, n).astype(float), _pois(lam)))
    for m in (0.3, 0.5, 0.6):
        out.append((f"intermittent negative binomial mean={m}",
                    lambda g, n, m=m: g.negative_binomial(1.0, 1.0 / (1.0 + m), n).astype(float), _negbin(m)))
    for lam in (0.3, 0.5, 0.65):
        out.append((f"intermittent Poisson lambda={lam}",
                    lambda g, n, lam=lam: g.poisson(lam, n).astype(float), _pois(lam)))
    return out


# ------------------------------------------------------------------ the second set (plan PL2)

def _log_skewed(k, sd, sign):
    """y = 100 exp(X), X = sign * sd * (G - k) / sqrt(k), G ~ Gamma(k, 1): ln y has standard deviation sd and
    skewness sign * 2 / sqrt(k) (right-skewed for sign +1, left-skewed for -1)."""
    def draw(g, n):
        return 100.0 * np.exp(sign * sd * (g.gamma(k, 1.0, n) - k) / math.sqrt(k))
    return draw


def _bimodal(p):
    """y = 100 exp(X), X ~ N(-1, 0.25^2) with probability p, N(+1, 0.25^2) otherwise: two modes, at about 37 and 272."""
    def draw(g, n):
        return 100.0 * np.exp(np.where(g.random(n) < p, g.normal(-1.0, 0.25, n), g.normal(1.0, 0.25, n)))
    return draw


def second():
    """(name, sampler(g, n), flags) for the settings of the second plan. flags["mape_infinite"]: the distribution
    has positive density at zero, so the expected MAPE of every forecast is infinite."""
    out = []
    for lam in (1000.0, 200.0, 100.0, 50.0, 20.0, 10.0):          # mean 100, CV sqrt(100 / lam): 0.32 to 3.16
        out.append((f"inverse Gaussian shape={lam:g}", lambda g, n, lam=lam: g.wald(100.0, lam, n), {}))
    for sd in (0.5, 1.0):                                          # ln y: standard deviation sd, skewness -2 .. +2
        for sign, side in ((-1, "left"), (1, "right")):
            for k in (1.0, 4.0, 16.0):
                out.append((f"log-skewed {side} sd={sd:g} skew={sign * 2 / math.sqrt(k):+g}", _log_skewed(k, sd, sign), {}))
        out.append((f"log-skewed symmetric sd={sd:g} skew=0", lambda g, n, sd=sd: g.lognormal(math.log(100.0), sd, n), {}))
    for p, what in ((0.5, "equal modes"), (0.65, "left mode higher"), (0.35, "right mode higher")):
        out.append((f"bimodal {what}", _bimodal(p), {}))
    for mu in (0.0, 0.5, 1.0):                                     # standard deviation 1: 50%, 31%, 16% negative
        out.append((f"normal straddling zero mean={mu:g}", lambda g, n, mu=mu: g.normal(mu, 1.0, n), {"mape_infinite": True}))
    return out
