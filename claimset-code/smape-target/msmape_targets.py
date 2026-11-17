"""The forecasts that minimise expected msMAPE and expected sMAPE, computed from the distributions themselves (no
sampling): away from zero, near zero (zero-inflated lognormal outcomes), and on Poisson counts; and the share of zeros
at which each measure switches to the forecast 0. Writes its results to $CLAIMSET_RESULTS when the lab runs it.

Settings (agreed with the supervisor):
  away from zero  lognormal sigma 0.5 and 1, inverse Gaussian shape 10 and 100, log-skewed sd 1 with skewness -1 and +1,
                  all at median about 100 (dists.py)
  near zero       lognormal sigma 1 at median 1 and at median 100, with shares of zeros 0, 0.2, 0.4, 0.6 and 0.8
  counts          Poisson with means 0.3, 0.5, 1, 1.5, 3 and 10
"""
import json
import math
import os

from dists import inverse_gaussian, log_skewed, lognormal
from msmape import Continuous, Poisson, lognormal_at, target, zero_threshold

MEASURES = ("sMAPE", "msMAPE")


def settings():
    away = [Continuous(d) for d in (lognormal(0.5), lognormal(1.0), inverse_gaussian(10.0), inverse_gaussian(100.0),
                                    log_skewed(4.0, 1.0, -1), log_skewed(4.0, 1.0, 1))]
    near = [Continuous(lognormal_at(med, 1.0), p0) for med in (1.0, 100.0) for p0 in (0.0, 0.2, 0.4, 0.6, 0.8)]
    counts = [Poisson(m) for m in (0.3, 0.5, 1.0, 1.5, 3.0, 10.0)]
    return [("away from zero", d) for d in away] + [("near zero", d) for d in near] + [("counts", d) for d in counts]


def main():
    results = []

    def add(variable, value, **cond):
        results.append({"variable": variable, "value": value, "conditions": cond})

    for group, d in settings():
        hi = math.exp(d.d.support[1]) if isinstance(d, Continuous) else float(d.k[-1])
        found = {}
        for meas in MEASURES:
            f, v = target(d, meas, hi)
            found[meas] = f
            add("target", f, setting=d.name, group=group, measure=meas)
            add("quantile level of the target", d.level(f), setting=d.name, group=group, measure=meas)
            add("expected score at the target", 200 * v, setting=d.name, group=group, measure=meas)
            print(f"{group:15s} {d.name:45s} {meas:7s} target {f:.6g} level {d.level(f):.4f} score {200 * v:.4f}",
                  flush=True)
        a, b = found["msMAPE"], found["sMAPE"]
        add("relative difference of the msMAPE target from the sMAPE target",
            (a - b) / b if b > 0 else (0.0 if a == 0 else math.inf), setting=d.name, group=group)
        add("median of the outcomes", d.median(), setting=d.name, group=group)
        add("share of zero outcomes", float(d.p0), setting=d.name, group=group)
    for med in (1.0, 100.0):
        base = lognormal_at(med, 1.0)
        for meas in MEASURES:
            z = zero_threshold(base, meas)
            add("zero threshold", z, setting=f"lognormal median={med:g} sigma=1", group="near zero", measure=meas)
            print(f"zero threshold, lognormal median {med:g} sigma 1, {meas}: {z:.5f}", flush=True)
    results = [r for r in results if isinstance(r["value"], (int, float)) and math.isfinite(r["value"])]
    path = os.environ.get("CLAIMSET_RESULTS")
    if path:
        with open(path, "w") as fh:
            json.dump(results, fh)


if __name__ == "__main__":
    main()
