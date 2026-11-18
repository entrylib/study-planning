"""Exploratory: does Smyl's msMAPE (constant nu) reward a different forecast than sMAPE on Poisson counts?
For each mean and nu: the forecast with the lowest expected term (searched on 0 and a grid of step 0.001 up to 12),
its expected score (times 200), and the same for sMAPE. Writes its results to $CLAIMSET_RESULTS when the lab runs it."""
import json
import math
import os

import numpy as np

K = np.arange(0, 80)
GRID = np.concatenate([[0.0], np.arange(1, 12001) / 1000.0])


def pmf(lam):
    return np.array([math.exp(-lam + k * math.log(lam) - math.lgamma(k + 1)) for k in K])


def expected(lam, nu):
    """Expected term of every forecast of GRID; nu None: the sMAPE term."""
    p = pmf(lam)
    f = GRID[:, None]
    k = K[None, :].astype(float)
    den = f + (k if nu is None else np.maximum(k, nu))
    term = np.where(den > 0, np.abs(f - k) / np.where(den > 0, den, 1.0), 0.0)
    return term @ p


def main():
    res = []
    for lam in (0.1, 0.3, 0.5, 0.7, 0.9, 1.0, 1.5, 3.0):
        row = {}
        for label, nu in (("sMAPE", None), ("msMAPE nu=0.1", 0.1), ("msMAPE nu=0.5", 0.5), ("msMAPE nu=1", 1.0)):
            e = expected(lam, nu)
            i = int(np.argmin(e))
            row[label] = (float(GRID[i]), 200 * float(e[i]), 200 * float(e[0]), 200 * float(e[1000]))
            cond = {"mean": lam, "measure": label, "share_of_zeros": math.exp(-lam)}
            res += [{"variable": "target", "value": float(GRID[i]), "conditions": cond},
                    {"variable": "expected score at the target", "value": 200 * float(e[i]), "conditions": cond},
                    {"variable": "expected score of the forecast 0", "value": 200 * float(e[0]), "conditions": cond},
                    {"variable": "expected score of the forecast 1", "value": 200 * float(e[1000]), "conditions": cond}]
        print(f"mean {lam:<4} zeros {math.exp(-lam):.2f} | " + " | ".join(
            f"{k}: target {v[0]:g} ({v[1]:.1f}; at 0: {v[2]:.1f}, at 1: {v[3]:.1f})" for k, v in row.items()))
    path = os.environ.get("CLAIMSET_RESULTS")
    if path:
        with open(path, "w") as fh:
            json.dump(res, fh)


if __name__ == "__main__":
    main()
