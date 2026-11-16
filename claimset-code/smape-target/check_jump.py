"""Plan PL6: numerical check of H31 (actuals of both signs: the sMAPE target jumps across zero with the mean).

For normal actuals with standard deviation 1 and mean mu just below and just above 0, computed from the distribution
by quadrature (theory.StraddlingTheory): every local minimum of the expected sMAPE, the sign of the global minimiser
(+1 or -1), its distance from 0, and the expected sMAPE term at forecasts very close to 0 (H31 says they score close to
1, the worst value).

    python3 check_jump.py
"""
import json
import os

import numpy as np

import dists
import theory

MEANS = (-0.5, -0.2, -0.05, -0.01, 0.01, 0.05, 0.2, 0.5)


def main():
    out = os.environ.get("CLAIMSET_RUN_DIR", ".")
    res = []
    for mu in MEANS:
        T = theory.theory(dists.Straddling(mu))
        mins = T.smape_minima()
        f0, v0 = mins[0]
        c = {"hypothesis": "H31", "mean": mu}
        res += [{"variable": "sign of the global minimiser of expected sMAPE", "value": int(np.sign(f0)), "conditions": c},
                {"variable": "distance of the global minimiser from 0", "value": float(abs(f0)), "conditions": c},
                {"variable": "lowest expected sMAPE term", "value": float(v0) / 200, "conditions": c},
                {"variable": "local minima of expected sMAPE", "value": len(mins), "conditions": c},
                {"variable": "smallest distance of any local minimiser from 0", "value": float(min(abs(f) for f, _ in mins)),
                 "conditions": c}]
        for f in (-1e-3, 1e-3):
            res.append({"variable": "expected sMAPE term at a forecast near 0", "value": float(T.smape_at(f)[0]) / 200,
                        "conditions": {**c, "forecast": f}})
        print(f"mu={mu:+}: global minimiser {f0:+.4f} (expected term {v0 / 200:.4f}); minima {[round(f, 4) for f, _ in mins]}")
    path = os.environ.get("CLAIMSET_RESULTS") or os.path.join(out, "results.json")
    with open(path, "w") as fh:
        json.dump(res, fh, indent=1)


if __name__ == "__main__":
    main()
