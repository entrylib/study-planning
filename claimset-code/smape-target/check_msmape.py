"""Numerical checks of the msMAPE hypothesis (parts a, b, c) on random pairs of actuals and forecasts: the deciding run
of its formal plan. Writes its results to $CLAIMSET_RESULTS when the lab runs it."""
import argparse
import json
import os

import numpy as np

from msmape import EPS, m_term, s_term


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--pairs", type=int, default=200000)
    a = ap.parse_args()
    g = np.random.default_rng(a.seed)
    n = a.pairs
    # actuals and forecasts of any sign, over many scales, with exact zeros and values near the switch at 0.5
    scale = 10.0 ** g.uniform(-4, 4, n)
    y = g.choice([-1, 1], n) * scale * g.random(n)
    f = g.choice([-1, 1], n) * 10.0 ** g.uniform(-4, 4, n) * g.random(n)
    y[: n // 10] = 0.0
    f[n // 10: n // 5] = 0.0
    near = slice(n // 5, n // 5 + n // 10)
    y[near] = g.uniform(0, 0.5, n // 10)
    f[near] = 0.5 - y[near] + g.normal(0, 1e-3, n // 10)
    s, m = s_term(y, f), m_term(y, f)
    tot = np.abs(y) + np.abs(f)
    far, close = tot >= 0.5, tot < 0.5
    factor = tot[far] / (tot[far] + EPS)
    out = [
        ("maximum absolute difference between m and s (|y|+|f|)/(|y|+|f|+eps) where |y|+|f| >= 0.5",
         float(np.max(np.abs(m[far] - s[far] * factor)))),
        ("violations of 5/6 <= (|y|+|f|)/(|y|+|f|+eps) <= 1 where |y|+|f| >= 0.5",
         int(np.sum((factor < 5 / 6 - 1e-12) | (factor > 1)))),
        ("maximum absolute difference between m and |y-f|/(0.5+eps) where |y|+|f| < 0.5",
         float(np.max(np.abs(m[close] - np.abs(y[close] - f[close]) / (0.5 + EPS))))),
        ("violations of 0 <= m < 1", int(np.sum((m < 0) | (m >= 1)))),
        ("msMAPE term of a zero forecast of a zero actual", float(m_term(0.0, 0.0))),
        ("maximum absolute difference between m(0, f) and |f|/(0.5+eps) for 0 < |f| < 0.5",
         float(np.max(np.abs(m_term(0.0, np.linspace(1e-6, 0.4999, 1001)) - np.linspace(1e-6, 0.4999, 1001) / 0.6)))),
        ("pairs with |y|+|f| >= 0.5", int(far.sum())),
        ("pairs with |y|+|f| < 0.5", int(close.sum())),
    ]
    results = [{"variable": v, "value": x, "conditions": {"pairs": n, "seed": a.seed}} for v, x in out]
    for r in results:
        print(f"{r['variable']} = {r['value']}")
    path = os.environ.get("CLAIMSET_RESULTS")
    if path:
        with open(path, "w") as fh:
            json.dump(results, fh)


if __name__ == "__main__":
    main()
