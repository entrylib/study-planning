"""Plan PL2: where sMAPE and MAE both target the median (ln y symmetric and strictly unimodal, H5), which measure's
minimiser over a test set varies less?

For each setting (dists.pl2()), each test-set size n and each of R test sets: the sample median (MAE's minimiser) and the
sMAPE minimiser of the test set (minimiser.py, exact). Over test sets: the variance of the log of each, their ratio, the
relative efficiency RE(n) = Var(ln median) / Var(ln sMAPE minimiser) with a 99% bootstrap percentile interval, and
the ratio of their mean squared errors about the true median. From the distribution (theory.py, H11): the large-sample
values n Var(ln median), n Var(ln sMAPE minimiser) and their ratio, the asymptotic relative efficiency.

    python3 efficiency.py --seed 4 --reps 10000
"""
import argparse
import json
import math
import multiprocessing
import os

import numpy as np

import dists
import minimiser
import theory

SIZES = [6, 8, 12, 13, 14, 18, 24, 36, 48, 100, 250, 500, 1000, 10000]
BOOT = 2000
CHUNK = 500


def boot_ratio(a, b, g, boots=BOOT):
    """99% percentile interval of Var(a) / Var(b) over resampled test sets (paired)."""
    R = a.size
    out = np.empty(boots)
    for k in range(boots):
        i = g.integers(0, R, R)
        out[k] = a[i].var(ddof=1) / b[i].var(ddof=1)
    return float(np.percentile(out, 0.5)), float(np.percentile(out, 99.5))


def one_setting(job):
    """Every result of one setting; run in a worker process. Each (setting, n) draws from its own seed, so the
    results do not depend on how settings are shared among workers."""
    idx, seed, sizes, reps = job
    d = dists.pl2()[idx]
    results = []
    T = theory.theory(d)
    v_s, v_m, s2 = T.asymptotic_variances()
    cond = {"setting": d.name}
    results += [
        {"variable": "large-sample n Var of ln(sample median)", "value": v_m, "conditions": cond},
        {"variable": "large-sample n Var of ln(sMAPE minimiser)", "value": v_s, "conditions": cond},
        {"variable": "asymptotic relative efficiency", "value": v_m / v_s, "conditions": cond},
        {"variable": "quadrature mass of the density", "value": T.mass, "conditions": cond},
        {"variable": "quantile level of the sMAPE target", "value": float(T.level(np.array([math.exp(T.minima_plus[0][0])]))[0]),
         "conditions": cond},
    ]
    true = T.median_x
    for n in sizes:
        g = np.random.default_rng([seed, idx, n])
        med, sm = np.empty(reps), np.empty(reps)
        for i in range(0, reps, CHUNK):            # test sets drawn in chunks, one after another
            r = min(CHUNK, reps - i)
            x = np.log(d.sample(g, r * n).reshape(r, n))
            med[i:i + r] = np.median(x, axis=1)
            sm[i:i + r] = minimiser.log_minimisers(x)
        c = {**cond, "test-set size": n}
        vm, vs = float(med.var(ddof=1)), float(sm.var(ddof=1))
        lo, hi = boot_ratio(med, sm, np.random.default_rng([seed, idx, n, 1]))
        results += [
            {"variable": "variance of ln(sample median)", "value": vm, "conditions": c},
            {"variable": "variance of ln(sMAPE minimiser)", "value": vs, "conditions": c},
            {"variable": "relative efficiency", "value": vm / vs,
             "interval": {"level": 0.99, "lower": lo, "upper": hi, "method": "bootstrap percentile, 2000 resamples"},
             "conditions": c},
            {"variable": "ratio of mean squared errors about the true median",
             "value": float(np.mean((med - true) ** 2) / np.mean((sm - true) ** 2)), "conditions": c},
            {"variable": "mean of ln(sMAPE minimiser) minus ln(median)", "value": float(sm.mean() - true),
             "conditions": c},
        ]
        print(f"{d.name} n={n}: RE={vm / vs:.3f} [{lo:.3f}, {hi:.3f}]", flush=True)
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=4)
    ap.add_argument("--sizes", type=int, nargs="+", default=SIZES)
    ap.add_argument("--reps", type=int, default=10000)
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--workers", type=int, default=os.cpu_count() or 1,
                    help="settings run in parallel processes; the results do not depend on it")
    ap.add_argument("--out", default=".")
    a = ap.parse_args()
    out = os.environ.get("CLAIMSET_RUN_DIR", a.out)
    jobs = [(idx, a.seed, a.sizes, a.reps) for idx, d in enumerate(dists.pl2())
            if not a.only or d.name.startswith(tuple(a.only))]
    with multiprocessing.Pool(min(a.workers, len(jobs))) as pool:
        results = [r for part in pool.map(one_setting, jobs, chunksize=1) for r in part]   # in setting order
    path = os.environ.get("CLAIMSET_RESULTS") or os.path.join(out, "results.json")
    with open(path, "w") as fh:
        json.dump(results, fh)
    import plot_efficiency
    plot_efficiency.main(out)


if __name__ == "__main__":
    main()
