"""Plan PL3: the flatness and lack of sensitivity of sMAPE as the spread of ln y grows (H14), on lognormal outcomes
(ln y symmetric) and log-Gumbel outcomes (ln y skewed right or left), sigma = 0.1 to 10 (dists.pl3()).

From each distribution (theory.py):
  - the quantile level of the sMAPE target, beside the median (0.5) and the log-mode;
  - expected sMAPE at the target, and the large-spread approximation 200 - 800 ln 2 g(c) of H14 there;
  - for forecasts at quantile levels 0.1, 0.25, 0.75 and 0.9: their excess expected score over each measure's own
    target, relative to the score at the target, for sMAPE, MAE (where E Y is finite) and log-MAE (100 E|ln y - ln f|,
    unbounded, the small-error limit of sMAPE, H9); for sMAPE also the absolute excess;
  - expected sMAPE over quantile levels 0.01 to 0.99 (expected_smape_curves.csv).
From test sets:
  - ranking agreement (lognormal outcomes, where all three measures target the median): the share of test sets of
    n outcomes in which the measure scores the median lower than the forecast at quantile level tau (ties count as
    disagreement), with its 99% Wilson interval;
  - the quantile level of the sMAPE minimiser of a test set (minimiser.py): its mean (with 99% interval) and 5th,
    25th, 50th, 75th and 95th percentiles, on every setting.

    python3 flatness.py --seed 5 --reps 10000
"""
import argparse
import csv
import json
import math
import os

import numpy as np

import dists
import minimiser
import theory
from measures import TestSet, log_mae, mae, smape

Z99 = 2.5758293035489004
TAUS = (0.1, 0.25, 0.75, 0.9)
AGREE_SIZES = (12, 48, 250, 1000)
MIN_SIZES = (48, 1000)
CURVE_LEVELS = np.round(np.arange(1, 100) / 100, 2)
CHUNK = 500


def wilson(k, n, z=Z99):
    p = k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return max(0.0, c - h), min(1.0, c + h)


def log_mae_expected(T, f):
    gr = T.grid
    return np.array([100 * float(np.abs(gr.x - math.log(ff)) @ gr.gw) for ff in f])


def distribution_side(d, T, cond):
    res, curve = [], []
    c_t, s_t = T.minima_plus[0]
    target = math.exp(c_t)
    res += [{"variable": "quantile level of the sMAPE target", "value": float(T.level(np.array([target]))[0]),
             "conditions": cond},
            {"variable": "quantile level of the log-mode", "value": float(T.level(np.array([math.exp(d.mode_x)]))[0]),
             "conditions": cond},
            {"variable": "expected sMAPE at the target", "value": s_t, "conditions": cond},
            {"variable": "large-spread approximation of expected sMAPE at the target",
             "value": float(T.flat_approximation(c_t)), "conditions": cond},
            {"variable": "local minima of expected sMAPE", "value": len(T.minima_plus), "conditions": cond}]
    median = float(T.quantiles([0.5])[0])
    f = T.quantiles(np.array(TAUS))
    exp = T.expected(np.concatenate([[target, median], f]))
    sm = exp["sMAPE"][0]
    lm = log_mae_expected(T, np.concatenate([[median], f]))
    for i, tau in enumerate(TAUS):
        c = {**cond, "quantile level of the forecast": tau}
        res += [{"variable": "excess expected score over the target, relative", "value": float(sm[2 + i] / sm[0] - 1),
                 "conditions": {**c, "measure": "sMAPE"}},
                {"variable": "excess expected sMAPE over the target", "value": float(sm[2 + i] - sm[0]),
                 "conditions": c},
                {"variable": "excess expected score over the target, relative", "value": float(lm[1 + i] / lm[0] - 1),
                 "conditions": {**c, "measure": "log-MAE"}}]
        if d.finite["mean"]:
            ma = exp["MAE"][0]
            res.append({"variable": "excess expected score over the target, relative",
                        "value": float(ma[2 + i] / ma[1] - 1), "conditions": {**c, "measure": "MAE"}})
        else:
            res.append({"variable": "excess expected score over the target, relative", "value": "none",
                        "conditions": {**c, "measure": "MAE", "reason": "E Y infinite: every expected MAE is infinite"}})
    fq = T.quantiles(CURVE_LEVELS)
    for lv, v in zip(CURVE_LEVELS, T.expected(fq)["sMAPE"][0]):
        curve.append([d.name, float(lv), float(v)])
    return res, curve, median, f


def agreement(d, median, f, g, reps, cond):
    res = []
    fs = np.concatenate([[median], f])
    for n in AGREE_SIZES:
        wins = {m: np.zeros(len(TAUS)) for m in ("sMAPE", "MAE", "log-MAE")}
        for i in range(reps):
            t = TestSet(d.sample(g, n))
            for m, fn in (("sMAPE", smape), ("MAE", mae), ("log-MAE", log_mae)):
                s = fn(t, fs)
                wins[m] += s[1:] > s[0]                 # the median scores strictly lower: agreement
        for m, w in wins.items():
            for j, tau in enumerate(TAUS):
                lo, hi = wilson(w[j], reps)
                res.append({"variable": "ranking agreement", "value": float(w[j] / reps),
                            "interval": {"level": 0.99, "lower": lo, "upper": hi},
                            "conditions": {**cond, "measure": m, "test-set size": n, "quantile level of the forecast": tau}})
    return res


def minimiser_levels(d, T, g, reps, cond):
    res = []
    for n in MIN_SIZES:
        lv = np.empty(reps)
        for i in range(0, reps, CHUNK):
            r = min(CHUNK, reps - i)
            x = np.log(d.sample(g, r * n).reshape(r, n))
            lv[i:i + r] = T.grid.level_x(minimiser.log_minimisers(x))
        mu, se = float(lv.mean()), float(lv.std(ddof=1)) / math.sqrt(reps)
        c = {**cond, "test-set size": n}
        res.append({"variable": "mean quantile level of the sMAPE minimiser of a test set", "value": mu,
                    "interval": {"level": 0.99, "lower": mu - Z99 * se, "upper": mu + Z99 * se}, "conditions": c})
        for p, v in zip((5, 25, 50, 75, 95), np.percentile(lv, (5, 25, 50, 75, 95))):
            res.append({"variable": "percentile of the quantile level of the sMAPE minimiser of a test set",
                        "value": float(v), "conditions": {**c, "percentile": p}})
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=5)
    ap.add_argument("--reps", type=int, default=10000)
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--out", default=".")
    a = ap.parse_args()
    out = os.environ.get("CLAIMSET_RUN_DIR", a.out)
    results, curves = [], []
    for idx, d in enumerate(dists.pl3()):
        if a.only and not d.name.startswith(tuple(a.only)):
            continue
        T = theory.theory(d)
        cond = {"setting": d.name}
        r, cv, median, f = distribution_side(d, T, cond)
        results += r
        curves += cv
        if d.family == "lognormal":
            results += agreement(d, median, f, np.random.default_rng([a.seed, idx, 1]), a.reps, cond)
        results += minimiser_levels(d, T, np.random.default_rng([a.seed, idx, 2]), a.reps, cond)
        print(f"{d.name}: done", flush=True)
    with open(os.path.join(out, "expected_smape_curves.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["setting", "quantile level of the forecast", "expected sMAPE"])
        w.writerows(curves)
    path = os.environ.get("CLAIMSET_RESULTS") or os.path.join(out, "results.json")
    with open(path, "w") as fh:
        json.dump(results, fh)
    import plot_flatness
    plot_flatness.main(out)


if __name__ == "__main__":
    main()
