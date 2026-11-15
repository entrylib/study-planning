"""Which fixed candidate forecast does each error measure rank first, in the limit and in test sets of n outcomes?

For each outcome distribution (settings.py) the candidates are fixed first (reference.py): the quantiles at levels
0.01 to 0.99, the mean, and the minimisers of the expected MAPE and sMAPE. Their expected scores give each measure's
preferred candidate. Then, for each test-set size n and each of R replications, a test set of n outcomes is drawn,
every candidate is scored by every measure, and the candidate with the lowest score is counted as that measure's
choice (a tie between k candidates counts 1/k for each). No forecast is estimated from a test set.

Outputs (in the run's directory when run under the lab, else in --out):
  results.json         per setting and measure: the preferred candidate (value, quantile level); per setting, measure
                       and n: the mean value and mean quantile level of the chosen candidates with 99% intervals, and
                       the share of test sets choosing the preferred candidate with its 99% Wilson interval
  candidates.csv       every candidate: labels, value, quantile level, expected score under each measure
  selection_shares.csv per setting, measure and n: the share of test sets choosing each candidate (non-zero shares)
  figures/*.png        one figure per setting (plot.py)

    python3 simulate.py --seed 1 --ref-seed 1001 --sizes 10 30 50 100 500 1000 100000 --reps-small 1000 --reps-large 200
"""
import argparse
import csv
import json
import math
import os

import numpy as np

import reference
import settings
from measures import MEASURES, TestSet

Z99 = 2.5758293035489004
TIE = 1e-10


def wilson(k, n, z=Z99):
    p = k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return max(0.0, c - h), min(1.0, c + h)


def choose(scores):
    """Weights of the chosen candidates: 1 for the lowest score, 1/k each for a tie of k."""
    low = scores.min()
    tied = scores <= low + TIE * max(abs(low), 1e-300)
    return tied / tied.sum()


def one_setting(idx, name, kind, cand, sampler, seed, sizes, reps_small, reps_large):
    measures = [m for m in MEASURES if not (kind == "count" and m == "MAPE")]
    results, shares = [], []
    for m in measures:
        i = cand.preferred[m]
        cond = {"setting": name, "measure": m}
        results += [{"variable": "value of the preferred candidate", "value": float(cand.values[i]),
                     "conditions": {**cond, "candidate": " = ".join(cand.labels[i])}},
                    {"variable": "quantile level of the preferred candidate", "value": float(cand.level[i]),
                     "conditions": {**cond, "candidate": " = ".join(cand.labels[i])}}]
    nonpositive = 0
    for n in sizes:
        reps = reps_large if n >= 100_000 else reps_small
        g = np.random.default_rng([seed, idx, n])
        tally = {m: np.zeros(cand.values.size) for m in measures}
        chosen_value = {m: np.empty(reps) for m in measures}
        chosen_level = {m: np.empty(reps) for m in measures}
        for r in range(reps):
            y = sampler(g, n)
            if kind == "continuous":
                nonpositive += int(np.sum(y <= 0))
            t = TestSet(y)
            for m in measures:
                w = choose(MEASURES[m](t, cand.values))
                tally[m] += w
                chosen_value[m][r] = w @ cand.values
                chosen_level[m][r] = w @ cand.level
        for m in measures:
            cond = {"setting": name, "measure": m, "test-set size": n}
            stats = [("mean value of the chosen candidates", chosen_value[m])]
            if kind == "continuous":                    # on counts many levels share one value: the value only
                stats.append(("mean quantile level of the chosen candidates", chosen_level[m]))
            for var, x in stats:
                mu, se = float(x.mean()), float(x.std(ddof=1)) / math.sqrt(reps)
                results.append({"variable": var, "value": mu,
                                "interval": {"level": 0.99, "lower": mu - Z99 * se, "upper": mu + Z99 * se},
                                "conditions": cond})
            k = float(tally[m][cand.preferred[m]])
            lo, hi = wilson(k, reps)
            results.append({"variable": "share of test sets choosing the preferred candidate", "value": k / reps,
                            "interval": {"level": 0.99, "lower": lo, "upper": hi}, "conditions": cond})
            for i in np.flatnonzero(tally[m]):
                shares.append([name, m, n, " = ".join(cand.labels[i]), float(cand.values[i]), float(cand.level[i]),
                               float(tally[m][i]) / reps])
    if kind == "continuous":
        results.append({"variable": "outcomes at or below zero", "value": nonpositive, "conditions": {"setting": name}})
    return results, shares


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--ref-seed", type=int, default=1001, help="seed of the reference samples of continuous settings")
    ap.add_argument("--sizes", type=int, nargs="+", default=[10, 30, 50, 100, 500, 1000, 100_000])
    ap.add_argument("--reps-small", type=int, default=1000)
    ap.add_argument("--reps-large", type=int, default=200)
    ap.add_argument("--only", nargs="*", help="run only the settings whose names start with these")
    ap.add_argument("--out", default=".")
    a = ap.parse_args()
    out = os.environ.get("CLAIMSET_RUN_DIR", a.out)
    all_settings = [(name, "continuous", s, None) for name, s in settings.continuous()] + \
                   [(name, "count", s, pmf) for name, s, pmf in settings.counts()]
    results, shares, cands = [], [], []
    for idx, (name, kind, sampler, pmf) in enumerate(all_settings):
        if a.only and not name.startswith(tuple(a.only)):
            continue
        if kind == "continuous":
            cand, _ = reference.continuous(sampler, [a.ref_seed, idx])
        else:
            cand = reference.counts(pmf)
        results += [{"variable": "mean of the outcome distribution", "value": cand.mean, "conditions": {"setting": name}},
                    {"variable": "median of the outcome distribution", "value": cand.median, "conditions": {"setting": name}}]
        if cand.zero_share is not None:
            results.append({"variable": "share of zero outcomes", "value": cand.zero_share, "conditions": {"setting": name}})
        for i in range(cand.values.size):
            cands.append([name, " = ".join(cand.labels[i]), float(cand.values[i]), float(cand.level[i])]
                         + [float(cand.expected[m][i]) if m in cand.expected else "" for m in MEASURES])
        r, sh = one_setting(idx, name, kind, cand, sampler, a.seed, a.sizes, a.reps_small, a.reps_large)
        results += r
        shares += sh
        print(f"{name}: done", flush=True)
    with open(os.path.join(out, "candidates.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["setting", "candidate", "value", "quantile level"] + [f"expected {m}" for m in MEASURES])
        w.writerows(cands)
    with open(os.path.join(out, "selection_shares.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["setting", "measure", "test-set size", "candidate", "value", "quantile level", "share"])
        w.writerows(shares)
    path = os.environ.get("CLAIMSET_RESULTS") or os.path.join(out, "results.json")
    with open(path, "w") as fh:
        json.dump(results, fh)
    import plot
    plot.main(out)


if __name__ == "__main__":
    main()
