"""Plan PL4: numerical checks of the formulas of the formal hypotheses H1 to H14.

A proof decides a formal hypothesis; this script checks each formula it proves on concrete numbers, so that an error in a
proof that changes a formula shows. Each result is a discrepancy (or a count of violations) that the proof says is zero,
or a quantity it says converges to a stated value.

    python3 check_formulas.py --seed 6
"""
import argparse
import json
import math
import os

import numpy as np

import dists
import minimiser
import theory

LN2 = math.log(2)


def s_term(y, f):
    num, den = np.abs(y - f), np.abs(y) + np.abs(f)
    return 200 * np.where(num == 0, 0.0, num / np.where(den == 0, 1.0, den))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=6)
    ap.add_argument("--out", default=".")
    a = ap.parse_args()
    out = os.environ.get("CLAIMSET_RUN_DIR", a.out)
    g = np.random.default_rng(a.seed)
    res = []

    def add(h, what, value, **cond):
        res.append({"variable": what, "value": value, "conditions": {"hypothesis": h, **cond}})

    # H1-H3, H9, H10: identities and inequalities on random positive pairs
    y = np.exp(g.uniform(-8, 8, 200000))
    f = np.exp(g.uniform(-8, 8, 200000))
    r, z = y / f, np.log(y) - np.log(f)
    s = s_term(y, f)
    add("H1", "maximum absolute difference from the ratio form", float(np.max(np.abs(s - 200 * np.abs(r - 1) / (r + 1)))))
    add("H2", "maximum absolute difference between the ratio and its reciprocal",
        float(np.max(np.abs(200 * np.abs(r - 1) / (r + 1) - 200 * np.abs(1 / r - 1) / (1 / r + 1)))))
    add("H3", "maximum absolute difference from 200 tanh(|z|/2)", float(np.max(np.abs(s - 200 * np.tanh(np.abs(z) / 2)))))
    # on |z| <= 12 the steps of 200 tanh(|z|/2) stay far above floating-point resolution (beyond about 30 it rounds to 200)
    t = np.linspace(0, 12, 120001)
    h = 200 * np.tanh(t / 2)
    add("H3", "violations of monotonicity, concavity or the bound 200 on 0 <= |z| <= 12",
        int(np.sum(np.diff(h) <= 0) + np.sum(np.diff(h, 2) > 1e-12) + np.sum(h > 200)))
    small = np.abs(z) < 0.5
    gap = 100 * np.abs(z) - s
    add("H9", "violations of 0 <= 100|z| - s <= 25|z|^3/3", int(np.sum((gap < -1e-9) | (gap > 25 / 3 * np.abs(z) ** 3 + 1e-9))))
    zz = np.abs(z[small & (np.abs(z) > 0.05)])
    add("H9", "maximum of |s - 100|z| + 25|z|^3/3| / |z|^5 for 0.05 < |z| < 0.5",
        float(np.max(np.abs(200 * np.tanh(zz / 2) - 100 * zz + 25 / 3 * zz ** 3) / zz ** 5)))
    e = y * g.uniform(0.001, 0.999, y.size)
    add("H10", "violations of s(y, y+e) < s(y, y-e)", int(np.sum(s_term(y, y + e) >= s_term(y, y - e))))
    k = np.exp(g.uniform(0.01, 5, y.size))
    add("H10", "maximum absolute difference between s(y, ky) and s(y, y/k)",
        float(np.max(np.abs(s_term(y, k * y) - s_term(y, y / k)))))

    # H4: the derivative formula against a numerical derivative of the exact expected sMAPE
    worst = 0.0
    for d in (dists.inverse_gaussian(20.0), dists.log_skewed(4.0, 1.0, 1), dists.bimodal(0.35), dists.lognormal(1.5)):
        T = theory.theory(d)
        gr = T.grid
        for lv in (0.2, 0.5, 0.8):
            c = float(gr.quantile_x(lv))
            hstep = 1e-4
            num = (T.smape_exact(c + hstep) - T.smape_exact(c - hstep)) / (2 * hstep)   # dS/dc = f dS/df
            formula = 200 * T.balance_exact(c)
            worst = max(worst, abs(num - formula) / max(1.0, abs(formula)))
    add("H4", "maximum relative difference between the formula for S' and a numerical derivative", worst)

    # H5: log-symmetric unimodal outcomes are minimised at the median
    devs = [abs(theory.theory(d).minima_plus[0][0] - math.log(100))
            for d in [dists.lognormal(s_) for s_ in (0.1, 0.5, 1, 3)] + [dists.log_laplace(b) for b in (0.25, 1, 3)]]
    counts = [len(theory.theory(d).minima_plus) for d in [dists.lognormal(1.0), dists.log_laplace(1.0)]]
    add("H5", "maximum distance in ln f between the minimiser and the median", float(max(devs)))
    add("H5", "number of local minima, lognormal sigma=1 and log-Laplace b=1", int(max(counts)))

    # H6: balance at the minimiser, the direction of the shift from D, the one-step level
    mism, bal, step_err = 0, 0.0, 0.0
    for d in dists.pl1():
        if d.kind != "positive":
            continue
        T = theory.theory(d)
        if len(T.minima_plus) != 1:
            continue
        gr = T.grid
        c = T.minima_plus[0][0]
        bal = max(bal, abs(T.balance_exact(c)))
        D, tau1 = T.one_step()
        lvl = float(gr.level_x(c))
        if abs(D) > 1e-9 and (D > 0) != (lvl > 0.5):
            mism += 1
        step_err = max(step_err, abs(tau1 - lvl))
    add("H6", "maximum balance residual at the minimiser (pull below minus pull above)", bal)
    add("H6", "settings with one local minimum where the sign of D misplaces the target", mism)
    add("H6", "maximum distance between the one-step and the exact quantile level", step_err)

    # H7: the exact criterion and the minima of the equal-mode mixture
    T = theory.theory(dists.bimodal(0.5))
    gr = T.grid
    zx = gr.x - math.log(100)
    crit = 2 * float(gr.density_x(math.log(100))) - float((np.tanh(np.abs(zx) / 2) / np.cosh(zx / 2) ** 2) @ gr.gw)
    add("H7", "2 g(m) - E[sech^2(Z/2) tanh(|Z|/2)] for the equal-mode bimodal setting (negative: median a local maximum)", crit)
    add("H7", "number of local minima of the equal-mode bimodal setting", len(T.minima_plus))
    add("H7", "difference between the two minimal expected sMAPEs, equal-mode bimodal setting",
        abs(T.minima_plus[0][1] - T.minima_plus[1][1]))
    left = theory.theory(dists.bimodal(0.65)).minima_plus[0][0] < math.log(100)
    right = theory.theory(dists.bimodal(0.35)).minima_plus[0][0] > math.log(100)
    add("H7", "unequal weights: global minimiser on the heavier side (both settings)", bool(left and right))

    # H8: the threshold against the share of zeros at which 0 and the positive minimiser tie
    worst = 0.0
    for sg in (0.3, 0.5, 1.0, 2.0):
        base = dists.lognormal(sg)
        Tb = theory.theory(base)
        lo, hi = 0.0, 1.0
        for _ in range(60):                     # bisection on S(0) - min_{f>0} S(f), computed separately
            p = (lo + hi) / 2
            if 200 * (1 - p) < 200 * p + (1 - p) * Tb.m_plus:
                hi = p
            else:
                lo = p
        worst = max(worst, abs(lo - Tb.zero_threshold()))
        add("H8", "zero threshold", Tb.zero_threshold(), setting=base.name)
    add("H8", "maximum distance between the formula and the bisected threshold", worst)
    ys = np.concatenate([g.normal(0, 1, 1000), [0.0] * 10])
    fs = np.concatenate([-ys[:500], np.zeros(500), [1.0] * 10])
    opp = ((np.sign(ys) * np.sign(fs) < 0) | ((ys == 0) != (fs == 0)))
    add("H8", "opposite-sign or one-zero pairs not scoring 200", int(np.sum(np.abs(s_term(ys[opp], fs[opp]) - 200) > 1e-9)))

    # H11: large-sample variance against simulation (lognormal sigma=1), with a 99% bootstrap interval
    d = dists.lognormal(1.0)
    T = theory.theory(d)
    v_s, v_m, _ = T.asymptotic_variances()
    for n, reps in ((1000, 4000), (10000, 2000)):
        sm = np.concatenate([minimiser.log_minimisers(np.log(d.sample(g, n * 500).reshape(500, n)))
                             for _ in range(reps // 500)])
        ratio = n * sm.var(ddof=1) / v_s
        boots = [n * sm[g.integers(0, sm.size, sm.size)].var(ddof=1) / v_s for _ in range(1000)]
        res.append({"variable": "simulated n Var(ln sMAPE minimiser) / large-sample value, lognormal sigma=1",
                    "value": float(ratio), "interval": {"level": 0.99, "lower": float(np.percentile(boots, 0.5)),
                                                         "upper": float(np.percentile(boots, 99.5))},
                    "conditions": {"hypothesis": "H11", "test-set size": n, "test sets": reps}})

    # H12: the MAPE target (1/y-weighted median) against direct minimisation; RMSE difference formula
    d = dists.lognormal(0.5)
    T = theory.theory(d)
    fgrid = np.exp(np.linspace(math.log(50), math.log(120), 4001))
    mp = T.expected(fgrid)["MAPE"][0]
    add("H12", "relative distance between the 1/y-weighted median and the minimiser of expected MAPE",
        float(abs(fgrid[int(np.argmin(mp))] / T.mape_minimiser() - 1)))
    d = dists.log_skewed(4.0, 0.5, 1)                            # E Y^2 finite: check the difference formula
    T = theory.theory(d)
    f1, f2 = 80.0, 130.0
    rm = T.expected(np.array([f1, f2]))["RMSE"][0] ** 2
    add("H12", "difference of expected squared errors minus (f-g)(f+g-2EY)",
        float((rm[0] - rm[1]) - (f1 - f2) * (f1 + f2 - 2 * d.moments["mean"])))

    # H13: bounds and limits
    ys = g.normal(0, 3, 100000)
    fs = g.normal(0, 3, 100000)
    st = s_term(ys, fs)
    add("H13", "sMAPE terms outside [0, 200]", int(np.sum((st < 0) | (st > 200))))
    T = theory.theory(dists.lognormal(1.0))
    add("H13", "expected sMAPE at f = 1e-6 and f = 1e10 (both near 200)",
        float(min(T.smape_exact(math.log(1e-6)), T.smape_exact(math.log(1e10)))))

    # H14: the expansion's error times sigma^3 stays bounded; sigma^3 RE tends to 7.23 (normal)
    gg = np.linspace(-6, 40, 200001)
    g2 = np.gradient(np.gradient(np.exp(-(gg + np.exp(-gg))), gg), gg)
    add("H14", "the proof's bound on the error times sigma^3, 600 zeta(3) max|g''|, log-Gumbel",
        float(600 * 1.2020569031595942 * np.max(np.abs(g2))))
    for sg in (5.0, 10.0, 20.0, 40.0):
        T = theory.theory(dists.log_gumbel(sg, 1))
        cs = T.grid.quantile_x(np.array([0.1, 0.3, 0.5, 0.7, 0.9]))
        err = max(abs(T.smape_exact(c) - T.flat_approximation(c)) for c in cs)
        add("H14", "maximum |S - (200 - 800 ln2 g)| times sigma^3, log-Gumbel right", float(err * sg ** 3), sigma=sg)
    for sg in (5.0, 10.0, 20.0, 40.0):
        v_s, v_m, _ = theory.theory(dists.lognormal(sg)).asymptotic_variances()
        add("H14", "sigma^3 times the asymptotic relative efficiency, lognormal", float(v_m / v_s * sg ** 3), sigma=sg)
    g0 = 1 / math.sqrt(2 * math.pi)
    add("H14", "limit of sigma^3 times the relative efficiency for normal log-outcomes, from the proof",
        3 * (800 * LN2) ** 2 / (32e4 * g0))

    path = os.environ.get("CLAIMSET_RESULTS") or os.path.join(out, "results.json")
    with open(path, "w") as fh:
        json.dump(res, fh, indent=1)
    for r_ in res:
        print(r_["conditions"]["hypothesis"], r_["variable"], "=", r_["value"])


if __name__ == "__main__":
    main()
