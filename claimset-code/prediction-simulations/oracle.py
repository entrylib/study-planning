"""Part A, the oracle design: which point forecast does each measure select when outcomes come from a known
distribution?

For each setting, test-set size n and replication: draw n outcomes, score every candidate on a fine grid by RMSE,
MAE, MAPE and sMAPE (measures.py), and record each measure's selected forecast (the candidate with the lowest score).
Over the replications, report the mean selected forecast with its 99% confidence interval. Per hypothesis, count
how many of its intervals (settings x sizes) fall on each side of its target, at the Bonferroni level for that many
intervals (1% for the hypothesis as a whole), allowing half a grid step for the grid (a selection cannot be closer to
a target than the grid allows).

    python3 oracle.py --seed 1 --sizes 10 30 50 100 500 1000 100000 --reps-small 1000 --reps-large 200
"""
import argparse
import json
import math
import os

import numpy as np

from measures import MEASURES, mape_left_out

Z99 = 2.5758293035489004                                 # two-sided 99%
STEP_REL = 0.005                                         # continuous grid: multiples of the median
STEP_ABS = 0.01                                          # count grid: absolute


def settings():
    out = [("normal", "continuous", dict(mean=100.0, median=100.0, symmetric=True),
            lambda g, n: g.normal(100.0, 10.0, n)),
           ("Laplace", "continuous", dict(mean=100.0, median=100.0, symmetric=True),      # sd 10: b = 10 / sqrt(2)
            lambda g, n: g.laplace(100.0, 10.0 / math.sqrt(2), n)),
           ("Student t df=3", "continuous", dict(mean=100.0, median=100.0, symmetric=True),  # sd 10: scale 10 / sqrt(3)
            lambda g, n: 100.0 + 10.0 / math.sqrt(3) * g.standard_t(3, n))]
    for s in (0.1, 0.25, 0.5, 1.0):
        out.append((f"lognormal s={s}", "continuous", dict(mean=100 * math.exp(s * s / 2), median=100.0, s=s),
                    lambda g, n, s=s: g.lognormal(math.log(100.0), s, n)))
    for b in (0.25, 0.5):                                 # y = 100 exp(Laplace(0, b)): mean 100 / (1 - b^2)
        out.append((f"log-Laplace b={b}", "continuous", dict(mean=100 / (1 - b * b), median=100.0, b=b),
                    lambda g, n, b=b: 100.0 * np.exp(g.laplace(0.0, b, n))))
    for lam, med in ((1.5, 1.0), (3.5, 3.0), (10.0, 10.0)):
        out.append((f"Poisson lambda={lam}", "count", dict(mean=lam, median=med),
                    lambda g, n, lam=lam: g.poisson(lam, n).astype(float)))
    for m in (0.3, 0.5, 0.6):
        out.append((f"negative binomial mean={m}", "intermittent", dict(mean=m, median=0.0, p0=1 / (1 + m)),
                    lambda g, n, m=m: g.negative_binomial(1, 1 / (1 + m), n).astype(float)))
    for lam in (0.3, 0.5, 0.65):
        out.append((f"Poisson lambda={lam}", "intermittent", dict(mean=lam, median=0.0, p0=math.exp(-lam)),
                    lambda g, n, lam=lam: g.poisson(lam, n).astype(float)))
    return out


def grid(kind, t):
    if kind == "continuous":
        return t["median"] * np.round(np.arange(0.2, 3.0 + 1e-9, STEP_REL), 6), t["median"] * STEP_REL
    top = max(3 * t["mean"], t["mean"] + 6 * math.sqrt(t["mean"]) + 2)
    return np.round(np.arange(0.0, top + 1e-9, STEP_ABS), 6), STEP_ABS


TIE = 1e-9                                               # scores within this relative distance of the lowest tie


def select(fn, y, f, step, kind, tries=20):
    """The candidate with the lowest score, extending the grid (by its own step) whenever the choice lands on an
    edge: the upper edge always, the lower edge down to zero (forecasts are not negative; for counts zero is the
    floor). When several candidates tie for the lowest score (MAE is flat between the two middle observations of an
    even sample), the midpoint of the tied ones is taken: for MAE, the usual sample median. Returns (choice, whether it
    is still on an edge after the extensions)."""
    for _ in range(tries):
        s = fn(y, f)
        low = s.min()
        tied = np.flatnonzero(s <= low + TIE * max(abs(low), 1e-300))
        lo, hi = int(tied[0]), int(tied[-1])
        lower_open = f[0] > 0 and kind == "continuous"
        if hi == f.size - 1:
            f = np.concatenate([f, f[-1] + step * np.arange(1, f.size + 1)])
        elif lo == 0 and lower_open:
            new = f[0] - step * np.arange(f.size, 0, -1)
            f = np.concatenate([new[new > 0], f])
        else:
            return float((f[lo] + f[hi]) / 2), False
    return float((f[lo] + f[hi]) / 2), True


def z_for(k, alpha=0.01):
    """Two-sided z of a Bonferroni interval, one of k at family level alpha."""
    lo, hi = 0.0, 10.0
    target = 1 - alpha / (2 * k)
    for _ in range(100):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if 0.5 * (1 + math.erf(mid / math.sqrt(2))) < target else (lo, mid)
    return (lo + hi) / 2


def run(seed, sizes, reps_small, reps_large, large=100_000):
    results, per = [], {}
    for name, _kind, t, _draw in settings():
        results += [{"variable": "mean of the outcome distribution", "value": t["mean"], "conditions": {"setting": name}},
                    {"variable": "median of the outcome distribution", "value": t["median"],
                     "conditions": {"setting": name}}]
    for n in sizes:
        reps = reps_large if n >= large else reps_small
        results += one_size(seed, n, reps, per)
    results += summaries(per)
    return results


def one_size(seed, n, reps, per):
    results = []
    for i, (name, kind, t, draw) in enumerate(settings()):
        f, step = grid(kind, t)
        sel = {m: [] for m in MEASURES}
        edge = {m: 0 for m in MEASURES}
        zeros, left_out = [], []
        for r in range(reps):
            g = np.random.default_rng([seed + r, i, n])
            y = draw(g, n)
            zeros.append(float(np.mean(y == 0)))
            left_out.append(mape_left_out(y))
            for m, fn in MEASURES.items():
                if m == "MAPE" and kind != "continuous":
                    continue                              # MAPE leaves out zeros: not studied on counts
                choice, at_edge = select(fn, y, f, step, kind)
                sel[m].append(choice)
                edge[m] += at_edge
        for m, v in sel.items():
            if not v:
                continue
            v = np.asarray(v)
            mean, se = float(v.mean()), float(v.std(ddof=1)) / math.sqrt(v.size)
            per[(name, m, n)] = dict(mean=mean, se=se, step=step, t=t, kind=kind, name=name)
            cond = {"setting": name, "measure": m, "test-set size": n}
            results.append({"variable": "selected forecast", "value": mean,
                            "interval": {"level": 0.99, "lower": mean - Z99 * se, "upper": mean + Z99 * se},
                            "conditions": cond})
            results.append({"variable": "selections at the edge of the grid", "value": edge[m], "conditions": cond})
        sc = {"setting": name, "test-set size": n}
        results += [{"variable": "share of zero outcomes", "value": float(np.mean(zeros)), "conditions": sc},
                    {"variable": "outcomes left out of MAPE (zeros)", "value": int(np.sum(left_out)), "conditions": sc}]
        if t.get("symmetric"):                           # RMSE and MAE alike: MAE minus RMSE, per replication
            d = np.asarray(sel["MAE"]) - np.asarray(sel["RMSE"])
            mean, se = float(d.mean()), float(d.std(ddof=1)) / math.sqrt(d.size)
            per[("H3", name, n)] = dict(mean=mean, se=se, step=step)
            results.append({"variable": "MAE minus RMSE selected forecast", "value": mean,
                            "interval": {"level": 0.99, "lower": mean - Z99 * se, "upper": mean + Z99 * se},
                            "conditions": sc})
    return results



def side(p, target, z):
    """-1 if the interval (mean +/- z se) lies wholly below the target, +1 wholly above, 0 if it reaches it (half a
    grid step allowed on each side)."""
    h = p["step"] / 2
    if p["mean"] + z * p["se"] + h < target:
        return -1
    if p["mean"] - z * p["se"] - h > target:
        return 1
    return 0


def summaries(per):
    """Per hypothesis: how many settings fall on each side of its target."""
    out = []

    def count(hyp, measure, kinds, target, names=None):
        rows = [p for k, p in per.items() if len(k) == 3 and k[0] != "H3" and k[1] == measure and p["kind"] in kinds
                and (names is None or p["name"].startswith(names))]
        z = z_for(len(rows))
        sides = [side(p, p["t"][target], z) for p in rows]
        for label, v in (("intervals", len(sides)), ("intervals that reach the target", sides.count(0)),
                         ("intervals that lie below the target", sides.count(-1)),
                         ("intervals that lie above the target", sides.count(1)), ("Bonferroni z", z)):
            out.append({"variable": label, "value": v, "conditions": {"hypothesis": hyp}})

    allkinds = ("continuous", "count", "intermittent")
    count("RMSE selects the mean", "RMSE", allkinds, "mean")
    count("MAE selects the median", "MAE", allkinds, "median")
    skewed = ("lognormal", "log-Laplace")
    count("MAPE selects below the median", "MAPE", ("continuous",), "median", skewed)
    count("sMAPE selects above the mean", "sMAPE", ("continuous",), "mean", skewed)
    count("sMAPE selects above the median", "sMAPE", ("continuous",), "median", skewed)
    count("MAE selects zero", "MAE", ("intermittent",), "median")
    count("sMAPE selects zero", "sMAPE", ("intermittent",), "median")
    rows = [p for k, p in per.items() if k[0] == "H3"]
    z = z_for(len(rows))
    reach = sum(p["mean"] - z * p["se"] - p["step"] <= 0 <= p["mean"] + z * p["se"] + p["step"] for p in rows)
    hyp = "RMSE and MAE select the same on symmetric outcomes"
    out += [{"variable": "intervals", "value": len(rows), "conditions": {"hypothesis": hyp}},
            {"variable": "intervals that reach the target", "value": int(reach), "conditions": {"hypothesis": hyp}},
            {"variable": "Bonferroni z", "value": z, "conditions": {"hypothesis": hyp}}]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--sizes", type=int, nargs="+", default=[10, 30, 50, 100, 500, 1000, 100_000])
    ap.add_argument("--reps-small", type=int, default=1000, help="replications for test sets under 100,000")
    ap.add_argument("--reps-large", type=int, default=200, help="replications for test sets of 100,000 or more")
    a = ap.parse_args()
    results = run(a.seed, a.sizes, a.reps_small, a.reps_large)
    path = os.environ.get("CLAIMSET_RESULTS")
    if path:
        with open(path, "w") as fh:
            json.dump(results, fh)
    else:
        print(json.dumps(results, indent=1))


if __name__ == "__main__":
    main()
