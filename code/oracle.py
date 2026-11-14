"""Part A, the oracle design: which point forecast does each measure select when outcomes come from a known
distribution?

For each setting and replication: draw `n` outcomes, score every candidate on a fine grid by RMSE, MAE, MAPE and
sMAPE (measures.py), and record each measure's selected forecast (the candidate with the lowest score). Over the
replications, report the mean selected forecast with its 99% confidence interval (mean +/- z_0.995 sd / sqrt(R)),
and, per hypothesis, how many settings fall on each side of the target, allowing half a grid step for the grid
(a selection cannot be closer to a target than the grid allows).

    python3 oracle.py --seed 1 --reps 200 --n 100000
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
    out = [("normal", "continuous", dict(mean=100.0, median=100.0), lambda g, n: g.normal(100.0, 10.0, n))]
    for s in (0.1, 0.25, 0.5, 1.0):
        out.append((f"lognormal s={s}", "continuous", dict(mean=100 * math.exp(s * s / 2), median=100.0, s=s),
                    lambda g, n, s=s: g.lognormal(math.log(100.0), s, n)))
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


def run(seed, reps, n):
    results, per = [], {}
    for i, (name, kind, t, draw) in enumerate(settings()):
        f, step = grid(kind, t)
        sel = {m: [] for m in MEASURES}
        edge = {m: 0 for m in MEASURES}
        zeros, left_out = [], []
        for r in range(reps):
            g = np.random.default_rng([seed + r, i])
            y = draw(g, n)
            zeros.append(float(np.mean(y == 0)))
            left_out.append(mape_left_out(y))
            for m, fn in MEASURES.items():
                if m == "MAPE" and kind != "continuous":
                    continue                              # MAPE leaves out zeros: not studied on counts
                k = int(np.argmin(fn(y, f)))
                sel[m].append(f[k])
                edge[m] += k in (0, f.size - 1)
        for m, v in sel.items():
            if not v:
                continue
            v = np.asarray(v)
            mean, half = float(v.mean()), Z99 * float(v.std(ddof=1)) / math.sqrt(v.size)
            per[(name, m)] = dict(mean=mean, lower=mean - half, upper=mean + half, step=step, t=t, kind=kind)
            results.append({"variable": "selected forecast", "value": mean,
                            "interval": {"level": 0.99, "lower": mean - half, "upper": mean + half},
                            "conditions": {"setting": name, "measure": m}})
            results.append({"variable": "selections at the edge of the grid", "value": edge[m],
                            "conditions": {"setting": name, "measure": m}})
        results += [{"variable": "mean of the outcome distribution", "value": t["mean"], "conditions": {"setting": name}},
                    {"variable": "median of the outcome distribution", "value": t["median"], "conditions": {"setting": name}},
                    {"variable": "share of zero outcomes", "value": float(np.mean(zeros)), "conditions": {"setting": name}},
                    {"variable": "outcomes left out of MAPE (zeros)", "value": int(np.sum(left_out)),
                     "conditions": {"setting": name}}]
        if name == "normal":                              # H3: MAE minus RMSE selection, per replication
            d = np.asarray(sel["MAE"]) - np.asarray(sel["RMSE"])
            mean, half = float(d.mean()), Z99 * float(d.std(ddof=1)) / math.sqrt(d.size)
            per["H3"] = dict(mean=mean, lower=mean - half, upper=mean + half, step=step)
            results.append({"variable": "MAE minus RMSE selected forecast", "value": mean,
                            "interval": {"level": 0.99, "lower": mean - half, "upper": mean + half},
                            "conditions": {"setting": name}})
    results += summaries(per)
    return results


def side(p, target):
    """-1 if the interval lies wholly below the target, +1 wholly above, 0 if it reaches it (half a grid step
    allowed on each side)."""
    h = p["step"] / 2
    if p["upper"] + h < target:
        return -1
    if p["lower"] - h > target:
        return 1
    return 0


def summaries(per):
    """Per hypothesis: how many settings fall on each side of its target."""
    out = []

    def count(hyp, measure, kinds, target, names=None):
        rows = [p for (name, m), p in per.items() if isinstance(name, str) and m == measure and p["kind"] in kinds
                and (names is None or name.startswith(names))]
        sides = [side(p, p["t"][target]) for p in rows]
        for label, v in (("settings", len(sides)), ("settings whose interval reaches the target", sides.count(0)),
                         ("settings whose interval lies below the target", sides.count(-1)),
                         ("settings whose interval lies above the target", sides.count(1))):
            out.append({"variable": label, "value": v, "conditions": {"hypothesis": hyp}})

    allkinds = ("continuous", "count", "intermittent")
    count("RMSE selects the mean", "RMSE", allkinds, "mean")
    count("MAE selects the median", "MAE", allkinds, "median")
    count("MAPE selects below the median", "MAPE", ("continuous",), "median", "lognormal")
    count("sMAPE selects above the mean", "sMAPE", ("continuous",), "mean", "lognormal")
    count("sMAPE selects above the median", "sMAPE", ("continuous",), "median", "lognormal")
    count("MAE selects zero", "MAE", ("intermittent",), "median")
    count("sMAPE selects zero", "sMAPE", ("intermittent",), "median")
    p = per["H3"]
    out.append({"variable": "settings whose interval reaches the target",
                "value": int(p["lower"] - p["step"] <= 0 <= p["upper"] + p["step"]),
                "conditions": {"hypothesis": "RMSE and MAE select the same on normal outcomes"}})
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--reps", type=int, default=200)
    ap.add_argument("--n", type=int, default=100_000)
    a = ap.parse_args()
    results = run(a.seed, a.reps, a.n)
    path = os.environ.get("CLAIMSET_RESULTS")
    if path:
        with open(path, "w") as fh:
            json.dump(results, fh)
    else:
        print(json.dumps(results, indent=1))


if __name__ == "__main__":
    main()
