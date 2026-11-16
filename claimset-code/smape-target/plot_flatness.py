"""Figures of plan PL3 (flatness of sMAPE, H14)."""
import csv
import json
import os

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

INK, MUTED, GRID, SURFACE = "#0b0b0b", "#6b6a66", "#e6e5e1", "#fcfcfb"
MEASURE = {"sMAPE": "#1c5cab", "MAE": "#c2571a", "RMSE": "#2f7d4f", "MAPE": "#8a5cc2"}
FAMILY = {"lognormal": "#1c5cab", "log-Gumbel right": "#c2571a", "log-Gumbel left": "#2f7d4f"}


def _style(ax):
    ax.set_facecolor(SURFACE)
    ax.grid(color=GRID, lw=0.6)
    ax.tick_params(colors=MUTED, labelsize=7)
    for sp in ax.spines.values():
        sp.set_visible(False)


def _family(setting):
    return setting.rsplit(" sigma=", 1)[0]


def _sigma(setting):
    return float(setting.rsplit("sigma=", 1)[1])


def _series(res, variable, family, **cond):
    pts = []
    for r in res:
        c = r["conditions"]
        if r["variable"] != variable or _family(c["setting"]) != family or not isinstance(r["value"], (int, float)):
            continue
        if all(c.get(k) == v for k, v in cond.items()):
            pts.append((_sigma(c["setting"]), r["value"], r.get("interval", {}).get("lower"),
                        r.get("interval", {}).get("upper")))
    return sorted(pts)


def main(out):
    res = json.load(open(os.environ.get("CLAIMSET_RESULTS") or os.path.join(out, "results.json")))
    fig_dir = os.path.join(out, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    # 1. where the target lies: from the median (narrow) to the log-mode (wide)
    fig, ax = plt.subplots(figsize=(7, 4), facecolor=SURFACE)
    _style(ax)
    for fam, col in FAMILY.items():
        t = _series(res, "quantile level of the sMAPE target", fam)
        m = _series(res, "quantile level of the log-mode", fam)
        if t:
            ax.plot([p[0] for p in t], [p[1] for p in t], color=col, lw=1.6, marker="o", ms=3, label=f"{fam}: sMAPE target")
            ax.plot([p[0] for p in m], [p[1] for p in m], color=col, lw=1, ls="--", label=f"{fam}: log-mode")
    ax.axhline(0.5, color=MUTED, ls=":", lw=0.9)
    ax.text(1.0, 0.5, " median", transform=ax.get_yaxis_transform(), fontsize=7, color=MUTED, va="center")
    ax.set_xscale("log")
    ax.set_xlabel("spread of ln y (sigma)", fontsize=8, color=INK)
    ax.set_ylabel("quantile level", fontsize=8, color=INK)
    ax.legend(fontsize=7, frameon=False)
    ax.set_title("The sMAPE target moves from the median to the mode of ln y as the spread grows", fontsize=9,
                 color=INK, loc="left")
    fig.tight_layout()
    fig.savefig(os.path.join(fig_dir, "flatness_target.png"), dpi=150)
    plt.close(fig)

    # 2. expected sMAPE over the quantile level of the forecast
    rows = list(csv.DictReader(open(os.path.join(out, "expected_smape_curves.csv"))))
    fams = [f for f in FAMILY if any(_family(r["setting"]) == f for r in rows)]
    fig, axes = plt.subplots(1, len(fams), figsize=(4.2 * len(fams), 3.6), sharey=True, facecolor=SURFACE)
    for ax, fam in zip(np.atleast_1d(axes), fams):
        _style(ax)
        sets = sorted({r["setting"] for r in rows if _family(r["setting"]) == fam}, key=_sigma)
        cmap = plt.get_cmap("Blues")
        for k, s in enumerate(sets):
            pts = [(float(r["quantile level of the forecast"]), float(r["expected sMAPE"])) for r in rows if r["setting"] == s]
            ax.plot(*zip(*pts), color=cmap(0.35 + 0.65 * k / max(1, len(sets) - 1)), lw=1.2, label=f"sigma={_sigma(s):g}")
        ax.set_title(fam, fontsize=8, color=INK, loc="left")
        ax.set_xlabel("quantile level of the forecast", fontsize=8, color=INK)
    np.atleast_1d(axes)[0].set_ylabel("expected sMAPE", fontsize=8, color=INK)
    np.atleast_1d(axes)[-1].legend(fontsize=6, frameon=False, loc="lower right")
    fig.suptitle("Expected sMAPE flattens towards 200 as the spread of ln y grows", fontsize=9, color=INK)
    fig.tight_layout()
    fig.savefig(os.path.join(fig_dir, "flatness_expected_sMAPE.png"), dpi=150)
    plt.close(fig)

    # 3. relative excess of forecasts away from each measure's target
    offsets = (-1.0, -0.5, 0.5, 1.0)
    fig, axes = plt.subplots(1, 4, figsize=(13, 3.4), sharey=True, facecolor=SURFACE)
    for ax, k in zip(axes, offsets):
        _style(ax)
        for m, col in MEASURE.items():
            pts = _series(res, "excess expected score over the target, relative", "lognormal", measure=m,
                          **{"offset in sigma": k})
            pts = [p for p in pts if p[1] > 0]
            if pts:
                ax.plot([p[0] for p in pts], [p[1] for p in pts], color=col, marker="o", ms=3, lw=1.4, label=m)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(f"forecast = target x exp({k:+g} sigma)", fontsize=8, color=INK, loc="left")
        ax.set_xlabel("sigma of ln y", fontsize=8, color=INK)
    axes[0].set_ylabel("excess expected score / score at the target", fontsize=8, color=INK)
    axes[-1].legend(fontsize=7, frameon=False)
    fig.suptitle("Lognormal outcomes: how much worse a forecast away from each measure's target scores", fontsize=9,
                 color=INK)
    fig.tight_layout()
    fig.savefig(os.path.join(fig_dir, "flatness_sensitivity.png"), dpi=150)
    plt.close(fig)

    # 4. ranking agreement on test sets
    sizes = sorted({r["conditions"]["test-set size"] for r in res if r["variable"] == "ranking agreement"})
    if sizes:
        fig, axes = plt.subplots(2, len(sizes), figsize=(3.3 * len(sizes), 5.6), sharex=True, sharey=True,
                                 facecolor=SURFACE)
        for row, k in enumerate((-0.5, 0.5)):
            for col_i, n in enumerate(sizes):
                ax = axes[row, col_i]
                _style(ax)
                for m, col in MEASURE.items():
                    pts = _series(res, "ranking agreement", "lognormal", measure=m,
                                  **{"test-set size": n, "offset in sigma": k})
                    if not pts:
                        continue
                    s_, v, lo, hi = map(np.array, zip(*pts))
                    ax.fill_between(s_, lo, hi, color=col, alpha=0.18, lw=0)
                    ax.plot(s_, v, color=col, lw=1.4, marker="o", ms=2.5, label=m)
                ax.axhline(0.5, color=MUTED, ls=":", lw=0.9)
                ax.set_xscale("log")
                ax.set_ylim(0, 1.0)
                ax.set_title(f"n={n}, target vs target x exp({k:+g} sigma)", fontsize=8, color=INK, loc="left")
        axes[0, -1].legend(fontsize=7, frameon=False, loc="lower left")
        fig.supxlabel("sigma of ln y (lognormal outcomes)", fontsize=8, color=INK)
        fig.supylabel("share of test sets scoring the target lower", fontsize=8, color=INK)
        fig.suptitle("Ranking agreement: whether a test set prefers a measure's target to a forecast away from it "
                     "(99% intervals; 0.5 is chance)", fontsize=9, color=INK)
        fig.tight_layout()
        fig.savefig(os.path.join(fig_dir, "flatness_agreement.png"), dpi=150)
        plt.close(fig)

    # 5. the sMAPE minimiser of a test set
    msz = sorted({r["conditions"]["test-set size"] for r in res
                  if r["variable"] == "mean quantile level of the sMAPE minimiser of a test set"})
    if msz:
        fig, axes = plt.subplots(1, len(fams), figsize=(4.2 * len(fams), 3.6), sharey=True, facecolor=SURFACE)
        for ax, fam in zip(np.atleast_1d(axes), fams):
            _style(ax)
            for k, n in enumerate(msz):
                p5 = _series(res, "percentile of the quantile level of the sMAPE minimiser of a test set", fam,
                             **{"test-set size": n, "percentile": 5})
                p95 = _series(res, "percentile of the quantile level of the sMAPE minimiser of a test set", fam,
                              **{"test-set size": n, "percentile": 95})
                p50 = _series(res, "percentile of the quantile level of the sMAPE minimiser of a test set", fam,
                              **{"test-set size": n, "percentile": 50})
                col = ["#3987e5", "#1c5cab"][k % 2]
                ax.fill_between([p[0] for p in p5], [p[1] for p in p5], [p[1] for p in p95], color=col, alpha=0.18, lw=0)
                ax.plot([p[0] for p in p50], [p[1] for p in p50], color=col, lw=1.3, label=f"n={n}: median, 5-95%")
            t = _series(res, "quantile level of the sMAPE target", fam)
            ax.plot([p[0] for p in t], [p[1] for p in t], color=INK, ls="--", lw=1, label="sMAPE target")
            ax.set_xscale("log")
            ax.set_ylim(0, 1)
            ax.set_title(fam, fontsize=8, color=INK, loc="left")
            ax.set_xlabel("sigma of ln y", fontsize=8, color=INK)
        np.atleast_1d(axes)[0].set_ylabel("quantile level of the test set's sMAPE minimiser", fontsize=8, color=INK)
        np.atleast_1d(axes)[-1].legend(fontsize=7, frameon=False)
        fig.suptitle("Where the sMAPE minimiser of a test set falls, and how widely it spreads", fontsize=9, color=INK)
        fig.tight_layout()
        fig.savefig(os.path.join(fig_dir, "flatness_minimiser.png"), dpi=150)
        plt.close(fig)


if __name__ == "__main__":
    import sys
    main(sys.argv[1] if len(sys.argv) > 1 else os.environ.get("CLAIMSET_RUN_DIR", "."))
