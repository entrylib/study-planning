"""Part C: summing RMSE over series on different scales, against averaging RMSSE.

Each collection has 20 series of 120 observations, the last 12 held out:
- 19 small series, noisy local level in logs: log y_t = l_t + e_t, l_t = l_{t-1} + u_t, e_t ~ N(0, 0.25^2),
  u_t ~ N(0, 0.075^2), l_0 = log(100);
- 1 large series, a random walk in logs: log y_t = log y_{t-1} + w_t, w_t ~ N(0, 0.05^2), log y_0 = log(10000).
Two methods forecast every series: naive (the last observation) and SES with alpha fixed at 0.1 on the original
scale (level started at the first observation). Per series: RMSE over the 12 held-out values, and RMSSE (RMSE divided
by the in-sample RMSE of the one-step naive forecast).

Per collection, which method each criterion ranks first: summed RMSE; the large series alone; the majority of the 20
series (the method with the lower RMSE on at least 11 of them; a 10-10 tie has no majority and counts as a
disagreement); mean RMSSE. Over the collections: the share in which summed
RMSE agrees with the large series (H9), and in which mean RMSSE agrees with the majority (H10), each with its 99%
Wilson interval.

    python3 aggregation.py --seed 1 --collections 500
"""
import argparse
import json
import math
import os

import numpy as np

Z99 = 2.5758293035489004
T, H, N_SMALL = 120, 12, 19


def ses_fixed(y, alpha=0.1):
    lev = y[:, 0].copy()
    for t in range(1, y.shape[1]):
        lev = lev + alpha * (y[:, t] - lev)
    return lev


def collection(g):
    u = g.normal(0, 0.075, (N_SMALL, T))
    small = np.exp(math.log(100.0) + np.cumsum(u, axis=1) + g.normal(0, 0.25, (N_SMALL, T)))
    large = np.exp(math.log(10000.0) + np.cumsum(g.normal(0, 0.05, T)))[None, :]
    return np.vstack([small, large])


def scores(y):
    train, test = y[:, :T - H], y[:, T - H:]
    naive = train[:, -1:]
    ses = ses_fixed(train)[:, None]
    rmse_naive = np.sqrt(((test - naive) ** 2).mean(1))
    rmse_ses = np.sqrt(((test - ses) ** 2).mean(1))
    insample = np.sqrt((np.diff(train, axis=1) ** 2).mean(1))
    return rmse_naive, rmse_ses, rmse_naive / insample, rmse_ses / insample


def wilson(k, n, z=Z99):
    p = k / n
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return centre - half, centre + half


def run(seed, collections):
    agree_sum_large = agree_rmsse_majority = 0
    sum_picks_naive = majority_picks_ses = 0
    for c in range(collections):
        g = np.random.default_rng([seed, c])
        rn, rs, qn, qs = scores(collection(g))
        sum_first = "naive" if rn.sum() < rs.sum() else "SES"
        large_first = "naive" if rn[-1] < rs[-1] else "SES"
        wins = int(np.sum(rs < rn))                       # series on which SES has the lower RMSE, of 20
        majority_first = "SES" if wins >= 11 else "naive" if wins <= 9 else "tie"   # a tie agrees with nothing
        rmsse_first = "naive" if qn.mean() < qs.mean() else "SES"
        agree_sum_large += sum_first == large_first
        agree_rmsse_majority += rmsse_first == majority_first
        sum_picks_naive += sum_first == "naive"
        majority_picks_ses += majority_first == "SES"
    out = []
    for label, k, hyp in (("share of collections where summed RMSE ranks as the large series", agree_sum_large, "H9"),
                          ("share of collections where mean RMSSE ranks as the majority of series",
                           agree_rmsse_majority, "H10")):
        lo, hi = wilson(k, collections)
        out.append({"variable": label, "value": k / collections, "interval": {"level": 0.99, "lower": lo, "upper": hi},
                    "conditions": {"hypothesis": hyp}})
    out += [{"variable": "collections", "value": collections},
            {"variable": "collections where summed RMSE ranks naive first", "value": sum_picks_naive},
            {"variable": "collections where most series rank SES first", "value": majority_picks_ses}]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--collections", type=int, default=500)
    a = ap.parse_args()
    results = run(a.seed, a.collections)
    path = os.environ.get("CLAIMSET_RESULTS")
    if path:
        with open(path, "w") as fh:
            json.dump(results, fh)
    else:
        print(json.dumps(results, indent=1))


if __name__ == "__main__":
    main()
