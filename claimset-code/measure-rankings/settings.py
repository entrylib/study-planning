"""The outcome distributions. Continuous ones are positive (the supervisor's decision DC3): normal and Laplace sit at a
level of 100 with standard deviation 5, twenty standard deviations above zero; the outlier settings keep their outliers
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
    out = [("normal", lambda g, n: g.normal(100.0, 5.0, n)),
           ("Laplace", lambda g, n: g.laplace(100.0, 5.0 / math.sqrt(2), n))]
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
