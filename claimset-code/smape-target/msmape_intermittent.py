"""Which forecast does msMAPE reward on intermittent demand, beside sMAPE and the Suilin loss?

Part 1, exact: the forecast that minimises the expected term of each measure on each setting, and for the lognormal
demand sizes the share of zeros at which each measure switches to the forecast 0. Part 2, simulated: on test sets of
n actuals, the forecast with the lowest value of each measure.

    python3 msmape_intermittent.py --seed 9

Measures: sMAPE; msMAPE of Smyl (2025) as Kolassa (2026, Definition 4) states it, with nu = 0.1, 0.5 and 1; the Suilin
loss. Settings: Poisson counts (means 0.1, 0.3, 0.5, 0.7, 1, 1.5, 3), negative binomial counts (size 0.5; means 0.5, 1,
2) and lognormal demand sizes (median 1, sigma 1) with a share 0.2, 0.4, 0.6 or 0.8 of zeros.
Writes its results to $CLAIMSET_RESULTS and a figure to $CLAIMSET_RUN_DIR/figures when the lab runs it."""
import argparse
import json
import math
import os

import numpy as np

import msmape
from msmape import Continuous, Poisson, lognormal_at, m_term, s_term, target, zero_threshold
from theory import _gauss

NUS = (0.1, 0.5, 1.0)
SIZES = (12, 52, 365)
REPS = 2000


def ms_term(nu):
    def term(y, f):
        y, f = np.asarray(y, float), np.asarray(f, float)
        return np.abs(f - y) / (f + np.maximum(y, nu))
    return term


MEASURES = {"sMAPE": s_term, **{f"msMAPE nu={nu:g}": ms_term(nu) for nu in NUS}, "Suilin": m_term}
msmape.TERMS.update(MEASURES)                                  # target() and zero_threshold() look the terms up there


class NegBin(Poisson):
    """Negative binomial counts with the given mean and size (variance mean + mean^2 / size)."""

    def __init__(self, mean, size):
        self.mean, self.size, self.name = mean, size, f"negative binomial mean={mean:g} size={size:g}"
        k = np.arange(0, int(mean + 60 * math.sqrt(mean + mean * mean / size) + 200))
        logp = (np.array([math.lgamma(i + size) - math.lgamma(i + 1) for i in k]) - math.lgamma(size)
                + size * math.log(size / (size + mean)) + k * math.log(mean / (size + mean)))
        self.k, self.p = k.astype(float), np.exp(logp)
        self.p0 = float(self.p[0])

    def sample(self, g, shape):
        return g.negative_binomial(self.size, self.size / (self.size + self.mean), shape).astype(float)


class Pois(Poisson):
    def sample(self, g, shape):
        return g.poisson(self.mean, shape).astype(float)


class ZeroLognormal(Continuous):
    """Lognormal demand sizes with a share p0 of zeros; the quadrature is also cut where a term has a kink."""

    def expected(self, measure, f):
        term, d = msmape.TERMS[measure], self.d
        cuts = [math.log(v) for v in (*NUS, f, 0.5 - f) if v > 0]
        fn = lambda x: term(np.exp(x), f) * d.pdf_x(x)
        pos = sum(_gauss(fn, a, b, msmape.PANELS) for a, b in zip(*self._pieces(cuts))) / self.mass
        return float(self.p0 * term(0.0, f) + (1 - self.p0) * pos)

    def sample(self, g, shape):
        y = g.lognormal(0.0, 1.0, shape)
        y[g.random(shape) < self.p0] = 0.0
        return y


def settings():
    out = [("Poisson", Pois(m)) for m in (0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 3.0)]
    out += [("negative binomial", NegBin(m, 0.5)) for m in (0.5, 1.0, 2.0)]
    out += [("lognormal sizes with zeros", ZeroLognormal(lognormal_at(1.0, 1.0), p0, f"lognormal sizes with zeros p0={p0:g}"))
            for p0 in (0.2, 0.4, 0.6, 0.8)]
    return out


def minimisers(y, term, counts):
    """The forecast with the lowest mean term in each test set (rows of y); ties go to the lowest forecast."""
    reps, n = y.shape
    out = np.empty(reps)
    if counts:
        top = int(y.max()) + 1
        cand = np.round(np.arange(0, top * 20 + 1) / 20.0, 2)                 # step 0.05, the integers among them
        values = np.arange(0, top + 1, dtype=float)
        hist = np.stack([(y == v).sum(axis=1) for v in values], axis=1)       # reps x values
        table = term(values[:, None], cand[None, :])                           # values x candidates
        return cand[np.argmin(hist @ table, axis=1)]
    grid = np.concatenate([[0.0], np.geomspace(1e-3, 1e3, 400)])
    for a in range(0, reps, 50):
        yy = y[a:a + 50]
        cand = np.sort(np.concatenate([np.broadcast_to(grid, (yy.shape[0], grid.size)), yy], axis=1), axis=1)
        score = term(yy[:, :, None], cand[:, None, :]).mean(axis=1)
        out[a:a + 50] = cand[np.arange(yy.shape[0]), np.argmin(score, axis=1)]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=9)
    a = ap.parse_args()
    res, exact = [], {}

    def add(variable, value, **cond):
        res.append({"variable": variable, "value": value, "conditions": cond})

    sets = settings()
    for group, d in sets:
        counts = isinstance(d, Poisson)
        hi = float(d.k[np.searchsorted(np.cumsum(d.p), 1 - 1e-9)]) + 5 if counts else math.exp(d.d.support[1])
        add("share of zero outcomes", float(d.p0), setting=d.name, group=group)
        add("median of the outcomes", float(d.median()), setting=d.name, group=group)
        if counts:
            add("mode of the outcomes", float(d.k[np.argmax(d.p)]), setting=d.name, group=group)
            add("mean of the outcomes", float(d.k @ d.p), setting=d.name, group=group)
        else:
            add("mean of the outcomes", float((1 - d.p0) * math.exp(0.5)), setting=d.name, group=group)
        for meas in MEASURES:
            f, v = target(d, meas, hi)
            if f < 1e-9:
                f = 0.0
            exact[(d.name, meas)] = f
            add("target", f, setting=d.name, group=group, measure=meas)
            add("quantile level of the target", d.level(f), setting=d.name, group=group, measure=meas)
            add("expected score at the target", 200 * v, setting=d.name, group=group, measure=meas)
            print(f"{d.name:<38} {meas:<15} target {f:.4g}  level {d.level(f):.3f}  score {200 * v:.2f}", flush=True)
    for meas in MEASURES:
        z = zero_threshold(lognormal_at(1.0, 1.0), meas)
        add("zero threshold", z, setting="lognormal sizes with zeros", group="lognormal sizes with zeros", measure=meas)
        print(f"zero threshold, {meas}: {z:.4f}", flush=True)

    for si, (group, d) in enumerate(sets):
        counts = isinstance(d, Poisson)
        for n in SIZES:
            g = np.random.default_rng([a.seed, si, n])
            y = d.sample(g, (REPS, n))
            for meas, term in MEASURES.items():
                m = minimisers(y, term, counts)
                cond = dict(setting=d.name, group=group, measure=meas, n=n)
                add("share of test sets whose minimiser is 0", float(np.mean(m == 0.0)), **cond)
                add("median of the test-set minimisers", float(np.median(m)), **cond)
                if counts:
                    add("share of test sets whose minimiser is the target",
                        float(np.mean(np.abs(m - exact[(d.name, meas)]) < 1e-9)), **cond)
        print("simulated", d.name, flush=True)

    out = os.environ.get("CLAIMSET_RUN_DIR") or "."
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        os.makedirs(os.path.join(out, "figures"), exist_ok=True)
        fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharex=True)
        for ax, grp in zip(axes, ("Poisson", "negative binomial", "lognormal sizes with zeros")):
            for meas in MEASURES:
                pts = sorted((d.p0, exact[(d.name, meas)]) for gname, d in sets if gname == grp)
                ax.plot([p for p, _ in pts], [t for _, t in pts], marker="o", label=meas, alpha=0.8)
            ax.set_title(grp)
            ax.set_xlabel("share of zero actuals")
        axes[0].set_ylabel("forecast with the lowest expected score")
        axes[0].legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(out, "figures", "msmape_intermittent.png"), dpi=150)
    except ImportError:
        pass
    path = os.environ.get("CLAIMSET_RESULTS")
    if path:
        with open(path, "w") as fh:
            json.dump([r for r in res if isinstance(r["value"], (int, float)) and math.isfinite(r["value"])], fh)


if __name__ == "__main__":
    main()
