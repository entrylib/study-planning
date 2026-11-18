"""Numerical checks of the signal-to-noise hypothesis (H39) at two spreads, so that convergence to each limit shows:
the deciding run of its second plan. The first plan's run (check_v3.py at sigma = 200) left the absolute error's ratio
3.5% from its limit in one of six cases, against a tolerance of 3%; a limit is checked here by its distance shrinking.

    python3 check_h39.py

Writes its results to $CLAIMSET_RESULTS when the lab runs it."""
import json
import math
import os

import numpy as np

from check_v3 import LN2, gumbel, gumbel_cdf, integ, norm, norm_cdf


def main():
    res = []

    def add(what, value, **cond):
        res.append({"variable": what, "value": value, "conditions": {"hypothesis": "H39", **cond}})

    for name, h, cdf, lo, hi in (("normal", norm, norm_cdf, -12.0, 12.0), ("Gumbel", gumbel, gumbel_cdf, -6.0, 30.0)):
        for c1, c2 in ((0.0, 0.5), (0.0, 1.0), (-0.5, 1.5)):
            # the integrands change over a width 1/sigma at c1 and c2 only: fine pieces there, coarse ones elsewhere
            def parts(fn, sg):
                w = 60.0 / sg
                pts = sorted({lo, hi, c1 - w, c1, c1 + w, c2 - w, c2, c2 + w})
                tot = 0.0
                for p, q in zip(pts[:-1], pts[1:]):
                    near = (q - p) <= w + 1e-12
                    tot += integ(fn, p, q, 400 if near else 200)
                return tot

            def snr(dfn, sg):
                mean = parts(lambda e: dfn(e) * h(e), sg)
                var = parts(lambda e: dfn(e) ** 2 * h(e), sg) - mean * mean
                return abs(mean) / math.sqrt(var)
            h1, h2, H2 = float(h(np.array(c1))), float(h(np.array(c2))), cdf(c2)
            limit = abs(2 * H2 - 1) / (2 * math.sqrt(H2 * (1 - H2)))
            log_ratios = []
            for sg in (2, 10, 50, 200, 2000):
                cond = {"error_distribution": name, "c1": c1, "c2": c2, "sigma": sg}
                if sg <= 50:
                    log_ratios.append(snr(lambda e: sg * (np.abs(e - c2) - np.abs(e - c1)), sg))
                    continue
                formula = 2 * LN2 / math.sqrt(2 * LN2 - 1) * abs(h1 - h2) / math.sqrt(h1 + h2) / math.sqrt(sg)
                r = snr(lambda e: np.tanh(sg * np.abs(e - c2) / 2) - np.tanh(sg * np.abs(e - c1) / 2), sg)
                add("distance of sMAPE's signal-to-noise ratio over its large-spread formula from 1", float(abs(r / formula - 1)), **cond)
                rho = math.exp(-sg * (c2 - c1))

                def d_abs(e):                                           # (|y - f2| - |y - f1|) / f2
                    mid = 1 + rho - 2 * np.exp(np.minimum(sg * (e - c2), 0.0))
                    return np.where(e <= c1, 1 - rho, np.where(e >= c2, -(1 - rho), mid))
                add("distance of the absolute error's signal-to-noise ratio over its limit from 1",
                    float(abs(snr(d_abs, sg) / limit - 1)), **cond)
            add("relative spread of the absolute log error's signal-to-noise ratio over sigma = 2, 10, 50",
                float((max(log_ratios) - min(log_ratios)) / max(log_ratios)), error_distribution=name, c1=c1, c2=c2)
    mu = math.log(100)
    for sg in (0.5, 1.0, 2.0):
        for c1, c2 in ((0.0, 0.5), (0.0, 1.0), (-0.5, 1.5)):
            f1, f2 = math.exp(mu + sg * c1), math.exp(mu + sg * c2)

            def d_sq(e):
                yy = np.exp(mu + sg * e)
                return (yy - f2) ** 2 - (yy - f1) ** 2
            pts = sorted({-14.0, 24.0, c1, c2})
            mean = sum(integ(lambda e: d_sq(e) * norm(e), p, q, 400) for p, q in zip(pts[:-1], pts[1:]))
            second = sum(integ(lambda e: d_sq(e) ** 2 * norm(e), p, q, 400) for p, q in zip(pts[:-1], pts[1:]))
            ey = math.exp(mu + sg * sg / 2)
            formula = abs(2 * ey - f1 - f2) / (2 * ey * math.sqrt(math.exp(sg * sg) - 1))
            add("distance of the squared error's signal-to-noise ratio over its formula from 1, lognormal",
                float(abs(abs(mean) / math.sqrt(second - mean * mean) / formula - 1)), c1=c1, c2=c2, sigma=sg)
    for sg in (4.0, 6.0, 8.0):
        ey = math.exp(mu + sg * sg / 2)
        f1, f2 = math.exp(mu), math.exp(mu + sg)
        ratio = abs(2 * ey - f1 - f2) / (2 * ey * math.sqrt(math.exp(sg * sg) - 1))
        add("distance of exp(sigma^2 / 2) times the squared error's ratio from 1, lognormal, formula",
            float(abs(math.exp(sg * sg / 2) * ratio - 1)), c1=0.0, c2=1.0, sigma=sg)
    for r in res:
        print(r["variable"], {k: v for k, v in r["conditions"].items() if k != "hypothesis"}, f"{r['value']:.4g}")
    path = os.environ.get("CLAIMSET_RESULTS")
    if path:
        with open(path, "w") as fh:
            json.dump(res, fh)


if __name__ == "__main__":
    main()
