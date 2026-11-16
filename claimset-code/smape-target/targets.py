"""Plan PL1: the forecast each measure targets, derived exactly and checked on test sets.

For each outcome distribution (dists.pl1()), theory.py derives, from the distribution itself: every local minimum of
expected sMAPE (H4-H8, H13), the one-step quantile level of H6, the zero threshold of H8, which expected scores are
infinite (H12), and the target of each measure. The candidates are then fixed: the quantiles at levels 0.01 to 0.99,
the mean (when finite), the forecast 0, the log-mode, the minimiser of expected MAPE (when it exists) and every local
minimiser of expected sMAPE; candidates with equal values are merged. For each test-set size n and each of R test sets,
every candidate is scored by every measure and the one with the lowest score is counted as that measure's choice (a tie
of k candidates counts 1/k each; a test set in which every score is infinite chooses none, and is counted).

Outputs (in the run's directory when run under the lab, else in --out): results.json, candidates.csv,
selection_shares.csv, figures/*.png.

    python3 targets.py --seed 3 --reps 1000
"""
import argparse
import csv
import json
import math
import os

import numpy as np

import dists
import theory
from measures import MEASURES, TestSet

Z99 = 2.5758293035489004
TIE = 1e-10
SIZES = [6, 8, 12, 13, 14, 18, 24, 36, 48, 100, 250, 500, 1000, 10000]
PERCENTILES = (5, 25, 50, 75, 95)
LEVELS = np.round(np.arange(1, 100) / 100, 2)
MERGE = 1e-4                 # candidates closer than this, relative, are one (minimisers are found to about 1e-5)


def wilson(k, n, z=Z99):
    p = k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return max(0.0, c - h), min(1.0, c + h)


def choose(scores):
    low = scores.min()
    if not np.isfinite(low):
        return None
    tied = scores <= low + TIE * max(abs(low), 1e-300)
    return tied / tied.sum()


class Candidates:
    def __init__(self, labels, values, level):
        self.labels, self.values, self.level = labels, values, level


def candidates(d, T, smape_minima, mape_min):
    q = T.quantiles(LEVELS)
    pairs = [(f"q{p:.2f}", float(v)) for p, v in zip(LEVELS, q)]
    pairs[49] = ("median (q0.50)", pairs[49][1])
    if d.finite["mean"]:
        pairs.append(("mean", float(d.moments["mean"]) * (1 - d.p0)))
    pairs.append(("zero", 0.0))
    if d.kind != "straddling":
        pairs.append(("log-mode", math.exp(d.mode_x)))
    if mape_min is not None and mape_min != 0.0:
        pairs.append(("MAPE minimiser", mape_min))
    for r, (f, _) in enumerate(smape_minima, 1):
        pairs.append(("sMAPE minimiser" if r == 1 else f"sMAPE local minimiser {r}", float(f)))
    pairs.sort(key=lambda p: p[1])
    labels, values = [], []
    for lab, v in pairs:
        if values and abs(v - values[-1]) <= MERGE * max(1.0, abs(v)):
            labels[-1].append(lab)
        else:
            labels.append([lab])
            values.append(v)
    values = np.array(values)
    for i, lab in enumerate(labels):              # a group holding the forecast 0 is exactly 0 (MAPE scores it 100)
        if "zero" in lab:
            values[i] = 0.0
    return Candidates(labels, values, T.level(values))


def derived(d, T, cand, smape_minima, mape_min):
    """Results computed from the distribution itself."""
    cond = {"setting": d.name}
    res = [{"variable": "median of the outcome distribution", "value": float(T.quantiles([0.5])[0]), "conditions": cond}]
    for m in ("mean", "second", "inv"):
        res.append({"variable": f"moment finite: {dict(mean='E Y', second='E Y^2', inv='E[1/Y]')[m]}",
                    "value": bool(d.finite[m]), "conditions": cond})
    if d.finite["mean"]:
        res.append({"variable": "mean of the outcome distribution", "value": float(d.moments["mean"]) * (1 - d.p0),
                    "conditions": cond})
    if d.p0:
        res.append({"variable": "share of zero outcomes", "value": d.p0, "conditions": cond})
    for r, (f, v) in enumerate(smape_minima, 1):
        lv = float(T.level(np.array([f]))[0])
        c = {**cond, "rank": r}
        res += [{"variable": "derived local minimum of expected sMAPE", "value": float(f), "conditions": c},
                {"variable": "quantile level of a derived local minimum of expected sMAPE", "value": lv, "conditions": c},
                {"variable": "expected sMAPE at a derived local minimum", "value": float(v), "conditions": c}]
    if d.kind != "straddling":
        D, tau1 = T.one_step()
        res += [{"variable": "pull from above minus pull from below at the median", "value": D, "conditions": cond},
                {"variable": "one-step quantile level of the sMAPE target", "value": tau1, "conditions": cond},
                {"variable": "zero threshold of sMAPE", "value": T.zero_threshold(), "conditions": cond},
                {"variable": "quantile level of the log-mode", "value": float(T.level(np.array([math.exp(d.mode_x)]))[0]),
                 "conditions": cond}]
    res.append({"variable": "derived MAPE target", "value": mape_min if mape_min is not None else "none",
                "conditions": cond})
    exp = T.expected(cand.values)
    for m in MEASURES:
        sc, key = exp[m]
        mc = {**cond, "measure": m}
        res.append({"variable": "candidates whose expected score is infinite", "value": int(np.sum(~np.isfinite(sc))),
                    "conditions": mc})
        if key is None or not np.isfinite(key).any():
            res.append({"variable": "value of the preferred candidate", "value": "none", "conditions": mc})
            continue
        i = int(np.argmin(np.where(np.isfinite(key), key, np.inf)))
        lab = " = ".join(cand.labels[i])
        res += [{"variable": "value of the preferred candidate", "value": float(cand.values[i]),
                 "conditions": {**mc, "candidate": lab}},
                {"variable": "quantile level of the preferred candidate", "value": float(cand.level[i]),
                 "conditions": {**mc, "candidate": lab}}]
    return res, exp


def simulate(idx, d, cand, preferred, seed, sizes, reps):
    results, shares = [], []
    for n in sizes:
        g = np.random.default_rng([seed, idx, n])
        tally = {m: np.zeros(cand.values.size) for m in MEASURES}
        lev = {m: np.full(reps, np.nan) for m in MEASURES}
        infinite = {m: 0 for m in MEASURES}
        for r in range(reps):
            t = TestSet(d.sample(g, n))
            for m in MEASURES:
                w = choose(MEASURES[m](t, cand.values))
                if w is None:
                    infinite[m] += 1
                    continue
                tally[m] += w
                lev[m][r] = w @ cand.level
        for m in MEASURES:
            cond = {"setting": d.name, "measure": m, "test-set size": n}
            if infinite[m]:
                results.append({"variable": "test sets in which every candidate's score is infinite",
                                "value": infinite[m], "conditions": cond})
            x = lev[m][np.isfinite(lev[m])]
            if x.size >= 2:
                mu, se = float(x.mean()), float(x.std(ddof=1)) / math.sqrt(x.size)
                results.append({"variable": "mean quantile level of the chosen candidates", "value": mu,
                                "interval": {"level": 0.99, "lower": mu - Z99 * se, "upper": mu + Z99 * se},
                                "conditions": cond})
                for p, v in zip(PERCENTILES, np.percentile(x, PERCENTILES)):
                    results.append({"variable": "percentile of the chosen candidates' quantile levels", "value": float(v),
                                    "conditions": {**cond, "percentile": p}})
            if m in preferred:
                k = float(tally[m][preferred[m]])
                lo, hi = wilson(k, reps)
                results.append({"variable": "share of test sets choosing the preferred candidate", "value": k / reps,
                                "interval": {"level": 0.99, "lower": lo, "upper": hi}, "conditions": cond})
            for i in np.flatnonzero(tally[m]):
                shares.append([d.name, m, n, " = ".join(cand.labels[i]), float(cand.values[i]), float(cand.level[i]),
                               float(tally[m][i]) / reps, reps])
    return results, shares


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=3)
    ap.add_argument("--sizes", type=int, nargs="+", default=SIZES)
    ap.add_argument("--reps", type=int, default=1000)
    ap.add_argument("--only", nargs="*", help="run only the settings whose names start with these")
    ap.add_argument("--out", default=".")
    a = ap.parse_args()
    out = os.environ.get("CLAIMSET_RUN_DIR", a.out)
    results, shares, rows = [], [], []
    for idx, d in enumerate(dists.pl1()):
        if a.only and not d.name.startswith(tuple(a.only)):
            continue
        T = theory.theory(d)
        sm = T.smape_minima()
        mape_min = T.mape_minimiser()
        cand = candidates(d, T, sm, mape_min)
        res, exp = derived(d, T, cand, sm, mape_min)
        results += res
        preferred = {}
        for m in MEASURES:
            key = exp[m][1]
            if key is not None and np.isfinite(key).any():
                preferred[m] = int(np.argmin(np.where(np.isfinite(key), key, np.inf)))
        for i in range(cand.values.size):
            rows.append([d.name, " = ".join(cand.labels[i]), float(cand.values[i]), float(cand.level[i])]
                        + [float(exp[m][0][i]) for m in MEASURES])
        r, sh = simulate(idx, d, cand, preferred, a.seed, a.sizes, a.reps)
        results += r
        shares += sh
        print(f"{d.name}: done", flush=True)
    with open(os.path.join(out, "candidates.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["setting", "candidate", "value", "quantile level"] + [f"expected {m}" for m in MEASURES])
        w.writerows(rows)
    with open(os.path.join(out, "selection_shares.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["setting", "measure", "test-set size", "candidate", "value", "quantile level", "share", "test sets"])
        w.writerows(shares)
    path = os.environ.get("CLAIMSET_RESULTS") or os.path.join(out, "results.json")
    with open(path, "w") as fh:
        json.dump(results, fh)
    import plot_targets
    plot_targets.main(out)


if __name__ == "__main__":
    main()
