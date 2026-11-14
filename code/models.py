"""Part B: SES on logs, its median and its mean forecast, scored by RMSE and MAE.

Series: log y_t = l_t + e_t, l_t = l_{t-1} + u_t, e_t ~ N(0, s^2), u_t ~ N(0, (0.3 s)^2), l_0 = log(100), 120
observations, the last 12 held out. SES is fitted to log y on the first 108: alpha in {0.05, ..., 0.95} with the
lowest one-step sum of squared errors (level started at the first observation); sigma^2 is the mean squared one-step
error. Forecasts for h = 1..12: median exp(l_T), mean exp(l_T + v_h / 2), v_h = sigma^2 (1 + (h - 1) alpha^2).

Per series: the difference in MSE (mean forecast minus median forecast) divided by exp(2 l_T), and the difference in
MAE (median forecast minus mean forecast) divided by exp(l_T), so that series are comparable. Over the series of each
noise level: the mean difference with its 99.5% confidence interval (z; 99.5% per level keeps the two levels at 1%).

    python3 models.py --seed 1 --series 2000
"""
import argparse
import json
import math
import os

import numpy as np

Z995 = 2.807033768343811                                 # two-sided 99.5%
ALPHAS = np.round(np.arange(0.05, 0.951, 0.05), 2)
T, H = 120, 12


def simulate(g, n, s):
    u = g.normal(0, 0.3 * s, (n, T))
    level = math.log(100.0) + np.cumsum(u, axis=1)
    return level + g.normal(0, s, (n, T))                # log y


def ses(x):
    """Fit SES to each row of x (log series): (alpha, final level, sigma^2), alpha chosen on ALPHAS by SSE."""
    n, m = x.shape
    lev = np.repeat(x[:, :1], ALPHAS.size, axis=1)       # n x A
    sse = np.zeros((n, ALPHAS.size))
    for t in range(1, m):
        err = x[:, t:t + 1] - lev
        sse += err ** 2
        lev = lev + ALPHAS * err
    k = np.argmin(sse, axis=1)
    rows = np.arange(n)
    return ALPHAS[k], lev[rows, k], sse[rows, k] / (m - 1), k


def run(seed, series):
    results = []
    for j, s in enumerate((0.25, 0.5)):
        g = np.random.default_rng([seed, j])
        x = simulate(g, series, s)
        train, test = x[:, :T - H], np.exp(x[:, T - H:])
        alpha, lt, sigma2, k = ses(train)
        h = np.arange(1, H + 1)
        v = sigma2[:, None] * (1 + (h - 1) * alpha[:, None] ** 2)
        med = np.exp(lt)[:, None] * np.ones(H)
        mean = np.exp(lt[:, None] + v / 2)
        scale = np.exp(lt)
        d_mse = (((test - mean) ** 2).mean(1) - ((test - med) ** 2).mean(1)) / scale ** 2
        d_mae = (np.abs(test - med).mean(1) - np.abs(test - mean).mean(1)) / scale
        for label, d in (("scaled MSE of the mean forecast minus that of the median forecast", d_mse),
                         ("scaled MAE of the median forecast minus that of the mean forecast", d_mae)):
            m, half = float(d.mean()), Z995 * float(d.std(ddof=1)) / math.sqrt(d.size)
            results.append({"variable": label, "value": m,
                            "interval": {"level": 0.995, "lower": m - half, "upper": m + half},
                            "conditions": {"noise s": s}})
        results += [{"variable": "series", "value": series, "conditions": {"noise s": s}},
                    {"variable": "mean fitted alpha", "value": float(alpha.mean()), "conditions": {"noise s": s}},
                    {"variable": "series whose alpha is at the edge of its grid",
                     "value": int(np.sum((k == 0) | (k == ALPHAS.size - 1))), "conditions": {"noise s": s}}]
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--series", type=int, default=2000)
    a = ap.parse_args()
    results = run(a.seed, a.series)
    path = os.environ.get("CLAIMSET_RESULTS")
    if path:
        with open(path, "w") as fh:
            json.dump(results, fh)
    else:
        print(json.dumps(results, indent=1))


if __name__ == "__main__":
    main()
