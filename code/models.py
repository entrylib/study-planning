"""Part B: ETS(A,N,N) on logs (simple exponential smoothing, statsforecast's AutoETS), its median and its mean
forecast, scored by RMSE and MAE, over in-sample sizes.

Series: log y_t = l_t + e_t, l_t = l_{t-1} + u_t, u_t ~ N(0, (0.3 s)^2), l_0 = log(100), and e_t of standard
deviation s from one of three noise families (--noises):
  normal    e_t ~ N(0, s^2);
  Laplace   e_t ~ Laplace(0, s / sqrt(2)), heavier tails than normal;
  outliers  a contaminated normal: e_t ~ N(0, c^2) with probability 0.95 and N(0, (4c)^2) with probability 0.05,
            c = s / sqrt(1.75), so that 5% of observations are outliers four times as spread as the rest.
Tails heavier than Laplace on the log scale (Student t, Cauchy) would give y no finite mean, and so no mean forecast
to compare; the contaminated normal gives rare large outliers with a finite mean. All three are symmetric, so exp(l_t)
stays the median of y_t. The model's mean forecast below assumes normal noise, as ETS(A,N,N) does; under the other two
families it is the forecast an analyst would use, not the exact mean. The first
`n` observations are the sample the model is fitted to (n in --sizes), the next 12 are held out. ETS(A,N,N) is fitted
to log y by maximum likelihood (statsforecast AutoETS(model="ANN"), its defaults otherwise). Forecasts for h = 1..12:
the median forecast exp(m_h) and the mean forecast exp(m_h + v_h / 2), where m_h is the model's point forecast of
log y and v_h its forecast variance, read from its 68.27% prediction interval (half-width = one standard deviation).

Per series: the difference in MSE (mean forecast minus median forecast) divided by exp(2 m_1), and the difference in
MAE (median forecast minus mean forecast) divided by exp(m_1), so that series are comparable. Per noise family, noise
level and size: the mean difference with its 99% confidence interval. Per hypothesis (MSE: the mean forecast is
better; MAE: the median forecast is better), how many of its intervals (families x levels x sizes) lie wholly below
zero, wholly above zero
or reach it, at the Bonferroni level for that many intervals (1% for the hypothesis as a whole).

    python3 models.py --seed 1 --series 2000 --sizes 10 30 50 100 500 1000 --noises normal Laplace outliers
"""
import argparse
import json
import math
import os
import warnings

import numpy as np

warnings.filterwarnings("ignore")
from statsforecast.models import AutoETS  # noqa: E402

H = 12
ONE_SD = 68.26894921370858                               # the interval whose half-width is one standard deviation


NOISES = ("normal", "Laplace", "outliers")


def noise(g, family, s, shape):
    """Noise of standard deviation s from one family."""
    if family == "normal":
        return g.normal(0, s, shape)
    if family == "Laplace":
        return g.laplace(0, s / math.sqrt(2), shape)
    if family == "outliers":
        c = s / math.sqrt(0.95 + 0.05 * 16)
        return g.normal(0, c, shape) * np.where(g.random(shape) < 0.05, 4.0, 1.0)
    raise ValueError(family)


def simulate(g, n, s, length, family="normal"):
    u = g.normal(0, 0.3 * s, (n, length))
    level = math.log(100.0) + np.cumsum(u, axis=1)
    return level + noise(g, family, s, (n, length))      # log y


def fit_forecast(x):
    """ETS(A,N,N) on one log series: (point forecasts m_h, variances v_h, alpha), or None if the fit fails."""
    try:
        m = AutoETS(model="ANN").fit(x)
        fc = m.predict(h=H, level=[ONE_SD])
    except Exception:  # noqa: BLE001 - counted as a failure, reported
        return None
    mean = np.asarray(fc["mean"])
    sd = (np.asarray(fc[f"hi-{ONE_SD}"]) - np.asarray(fc[f"lo-{ONE_SD}"])) / 2
    return mean, sd ** 2, float(m.model_["par"][0])


def run(seed, series, sizes, noises=NOISES):
    results, diffs = [], {"MSE": [], "MAE": []}
    for family in noises:
        for j, s in enumerate((0.25, 0.5)):
            for size in sizes:
                g = np.random.default_rng([seed, NOISES.index(family), j, size])
                x = simulate(g, series, s, size + H, family)
                d_mse, d_mae, alphas, failed = [], [], [], 0
                for row in x:
                    out = fit_forecast(row[:size])
                    if out is None:
                        failed += 1
                        continue
                    m, v, alpha = out
                    test = np.exp(row[size:])
                    med, mean = np.exp(m), np.exp(m + v / 2)
                    scale = math.exp(m[0])
                    d_mse.append((np.mean((test - mean) ** 2) - np.mean((test - med) ** 2)) / scale ** 2)
                    d_mae.append((np.mean(np.abs(test - med)) - np.mean(np.abs(test - mean))) / scale)
                    alphas.append(alpha)
                cond = {"noise": family, "noise s": s, "in-sample size": size}
                for key, label, d in (("MSE", "scaled MSE of the mean forecast minus that of the median forecast", d_mse),
                                      ("MAE", "scaled MAE of the median forecast minus that of the mean forecast", d_mae)):
                    d = np.asarray(d)
                    mu, se = float(d.mean()), float(d.std(ddof=1)) / math.sqrt(d.size)
                    diffs[key].append((mu, se))
                    results.append({"variable": label, "value": mu,
                                    "interval": {"level": 0.99, "lower": mu - Z99 * se, "upper": mu + Z99 * se},
                                    "conditions": cond})
                a = np.asarray(alphas)
                results += [{"variable": "series", "value": series, "conditions": cond},
                            {"variable": "series whose fit failed", "value": failed, "conditions": cond},
                            {"variable": "mean fitted alpha", "value": float(a.mean()), "conditions": cond},
                            {"variable": "series whose alpha is at a bound (below 0.001 or above 0.999)",
                             "value": int(np.sum((a < 0.001) | (a > 0.999))), "conditions": cond}]
    for key, hyp in (("MSE", "RMSE ranks the mean forecast first"), ("MAE", "MAE ranks the median forecast first")):
        rows = diffs[key]
        z = z_for(len(rows))
        below = sum(mu + z * se < 0 for mu, se in rows)
        above = sum(mu - z * se > 0 for mu, se in rows)
        for label, v in (("intervals", len(rows)), ("intervals that lie below zero", below),
                         ("intervals that lie above zero", above), ("intervals that reach zero", len(rows) - below - above),
                         ("Bonferroni z", z)):
            results.append({"variable": label, "value": v, "conditions": {"hypothesis": hyp}})
    return results


Z99 = 2.5758293035489004


def z_for(k, alpha=0.01):
    """Two-sided z of a Bonferroni interval, one of k at family level alpha."""
    lo, hi = 0.0, 10.0
    for _ in range(100):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if 0.5 * (1 + math.erf(mid / math.sqrt(2))) < 1 - alpha / (2 * k) else (lo, mid)
    return (lo + hi) / 2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--series", type=int, default=2000)
    ap.add_argument("--sizes", type=int, nargs="+", default=[10, 30, 50, 100, 500, 1000])
    ap.add_argument("--noises", nargs="+", choices=NOISES, default=list(NOISES))
    a = ap.parse_args()
    results = run(a.seed, a.series, a.sizes, a.noises)
    path = os.environ.get("CLAIMSET_RESULTS")
    if path:
        with open(path, "w") as fh:
            json.dump(results, fh)
    else:
        print(json.dumps(results, indent=1))


if __name__ == "__main__":
    main()
